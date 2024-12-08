/*
 * $Id$
 *
 * Copyright and User License
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright 2006-2019 CERN and INFN
 * 
 *
 * Please consult the LICENSE file for the license 
 *
 * DISCLAIMER
 * ~~~~~~~~~~
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
 * NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY, OF
 * SATISFACTORY QUALITY, AND FITNESS FOR A PARTICULAR PURPOSE
 * OR USE ARE DISCLAIMED. THE COPYRIGHT HOLDERS AND THE
 * AUTHORS MAKE NO REPRESENTATION THAT THE SOFTWARE AND
 * MODIFICATIONS THEREOF, WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADE SECRET OR OTHER PROPRIETARY RIGHT.
 *
 * LIMITATION OF LIABILITY
 * ~~~~~~~~~~~~~~~~~~~~~~~
 * THE COPYRIGHT HOLDERS AND THE AUTHORS SHALL HAVE NO
 * LIABILITY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
 * CONSEQUENTIAL, EXEMPLARY, OR PUNITIVE DAMAGES OF ANY
 * CHARACTER INCLUDING, WITHOUT LIMITATION, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES, LOSS OF USE, DATA OR PROFITS,
 * OR BUSINESS INTERRUPTION, HOWEVER CAUSED AND ON ANY THEORY
 * OF CONTRACT, WARRANTY, TORT (INCLUDING NEGLIGENCE), PRODUCT
 * LIABILITY OR OTHERWISE, ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGES.
 *
 * Author:	Vasilis.Vlachoudis@cern.ch
 * Date:	31-Mar-2010
 */

#include <string>
#include <ostream>

#include <fenv.h>
//#include <ieee754.h>

#include "bmath.h"

#define FLT_SMALL	1.0E-30
#define FLT_HUGE	1.0E+30

using namespace std;

/** "fast" integer sqrt routine
 * @param n	number to find sqrt
 * @return integer square root
 *
 * WARNING: normal sqrt() from co-processor is x10 faster!
 */
unsigned isqrt(unsigned n)
{
	register unsigned root, remainder, place;

	root = 0;
	remainder = n;
	place = 0x40000000; // OR place = 0x4000; OR place = 0x40; - respectively

	while (place > remainder) place >>= 2;
	while (place) {
		if (remainder >= root + place) {
			remainder -= root + place;
			root += place << 1;
		}
		root  >>= 1;
		place >>= 2;
	}
	return root;
} // isqrt

/**
 * Cubic equation: y^3 + ay^2 + by + c = 0
 * @param a,b,c		polynomial coefficients
 * @param y[]		array filled with solution
 * @param acc		required accuracy
 * @param iter		newton iteration to improve accuracy
 * @return number of real solution 1 or 3
 */
int cubic(const double a, const double b, const double c, double x[],
	  const double acc, const int iter)
{
	double Q = (a*a - 3.*b) / 9.;
	double R = ((2.*a*a - 9.*b)*a + 27.*c)/54.;
	double R2 = R*R;
	double Q3 = Q*Q*Q;
	int num;

	if (R2 < Q3) {	// the cubic has 3 real solutions
		double theta = acos(R/sqrt(Q3));
		double sqrt_Q = sqrt(Q);
		x[0] = -2. * sqrt_Q * cos(theta/3.) - a/3.;
		x[1] = -2. * sqrt_Q * cos((theta+2.*PI)/3.) - a/3.;
		x[2] = -2. * sqrt_Q * cos((theta-2.*PI)/3.) - a/3.;
		num = 3;
	} else {
		double A = -Sign(R) * cbrt(Abs(R) + sqrt(R2 - Q3));
		double B;
		if (Abs(A)>acc)
			B = Q / A;
		else
			B = 0.0;

		x[0] = (A+B) - a/3.;
		num = 1;
	}

	// Newton refining of solutions
	if (iter>0)
		for (int i=0; i<num; i++) {
			double z = x[i];
			double zprev = z;
			double fprev = FLT_HUGE;
			for (int j=0; j<iter; j++) {
				double fx = ((z + a)*z + b)*z + c;
				if (Eq0(fx, FLT_SMALL) || Abs(fx)>=Abs(fprev)) {
					// diverging solution
					x[i] = zprev;
					break;
				}
				zprev = z;
				fprev = fx;
				double dfx = (3.0*z + 2.0*a)*z + b;
				if (Eq0(dfx, FLT_SMALL)) {
					x[i] = z;
					break;
				}
				double e = fx/dfx;
				if (Eq0(e, acc)) {
					x[i] = z;
					break;
				}
				z -= e;
			}
		}

	return num;
} // cubic

/** Quartic equation	x^4 + ax^3 + bx^2 + cx + d = 0
 * @param a,b,c,d	coefficients
 * @param x[]		vector with real solutions
 * @param acc		required accuracy
 * @param iter		newton iteration to improve accuracy
 * @return number of real solutions
 */
int quartic(const double a, const double b, const double c, const double d,
		double x[], const double acc, const int iter)
{
	int num;
	double mbd = acc*Min(Abs(b), Abs(d));
	if (Abs(a) < mbd && Abs(c)<mbd) {
		// x^4 + bx^2 + d = 0
		num = quadratic(b, d, x, x+1, acc);
		if (num==0) return 0;

		if (x[0]<0.0 && x[1]<0.0) return 0;
		if (num>1 && x[1]<0.0) num--;
		if (x[0]<0.0) {
			x[0] = x[1];
			num--;
		}

		x[0] = sqrt(x[0]);
		x[num] = -x[0];
		if (num>1) {
			x[1] = sqrt(x[1]);
			x[3] = -x[1];
		} else {
			x[2] = x[0];
			x[3] = x[1];
		}
		return 2*num;
	}

	// solution based on Ferrari's lines
	// substitute x = y - a/4 to eliminate cubic term
	// x^4 + px^2 + qx + r = 0
	double sqa = Sqr(a);
	double p = -3.0*(sqa/8.0) + b;
	double q =  (sqa/8.0 - 0.5*b)*a + c;
	double r =  (-3.0*(sqa/256) + b/16.0)*sqa - (a*c)/4.0 + d;

	// no absolute term: y(y^3 + py + q) = 0
	if (Abs(r) <= acc*(Abs(p)+Abs(q))) {
		num = cubic(0.0, p, q, x, acc, iter);
		// add also the 0 solution
		x[num++] = 0.0;
	} else {
		// solve the resolvent cubic ...
		cubic(-p/2.0, -r, r*p/2.0 - Sqr(q)/8.0, x, acc, iter);

		// take the one real solution
		double z = x[0];

		// and build two quadric equations
		double u = z*z - r;
		double v = 2.0*z - p;

		if (u<=0.0 && u>=-acc*(z*z+Abs(r)))
			u = 0.0;
		else if (u > 0.0)
			u = sqrt(u);
		else
			return 0;

		if (v<=0.0 && v>=-acc*(2.0*Abs(z)+Abs(p)))
			v = 0.0;
		else if (v > 0.0)
			v = sqrt(v);
		else
			return 0;

		num = quadratic((q<0.0? -v:v), z-u, x, x+1, acc);
		if (num==1) {
			num=2;	// always add both solution
/*
			if (iter) {
				// try to differentiate a bit the solutions
				// Most probably we are in a local min or max
				// having derivative=0 and the Newton method
				// will not work
				x[0] *= 1.0-acc;
				x[1] *= 1.0+acc;
			}
*/
		}
		int m = quadratic((q<0.0? v:-v), z+u, x+num, x+num+1, acc);
		if (m==1) {
/*
			if (iter) {
				// same as above
				x[num]   *= 1.0-acc;
				x[num+1] *= 1.0+acc;
			}
*/
			num += 2;	// always add both solutions
		} else
		if (m) num += 2;
	}

	// resubstitute
	double sub = a/4.0;
	for (int i=0; i<num; i++)
		x[i] -= sub;

	// Newton refining of solutions
	if (iter>0)
		for (int i=0; i<num; i++) {
			double z = x[i];
			double zprev = z;
			double fprev = FLT_HUGE;
			for (int j=0; j<iter; j++) {
				double fx = (((z + a)*z + b)*z + c)*z + d;
				if (Eq0(fx, FLT_SMALL) || Abs(fx)>=Abs(fprev)) {
					// diverging solution
					x[i] = zprev;
					break;
				}
				zprev = z;
				fprev = fx;
				double dfx = ((4.0*z+3.0*a)*z + 2.0*b)*z + c;
				if (Eq0(dfx, FLT_SMALL)) {
					x[i] = z;
					break;
				}
				double e = fx/dfx;
				if (Eq0(e, acc)) {
					x[i] = z;
					break;
				}
				z -= e;
			}
		}

	return num;
} // quartic

/** polpol POLynomial root POLishing
 * @param C	polynomial coefficients
 * @param n	number of coefficients
 * @param x	first guess of the solution (I/O)
 * @param eps	required accuracy
 * @return	return 0 on success, 1 on failure, 2 could not achieve
 *		required accuracy
 */
int polpol(const double C[], const int n, double *x, const double eps)
{
	const int MAXC = 100;	// maximum number of iterations
	double xorg = *x;

	double xabs = Abs(*x)*eps;
	for (int cycle=0; cycle<MAXC; cycle++) {
		double x0 = *x;
		double p  = C[n] * *x + C[n-1];
		double p1 = C[n];
		for (int i=n-2; i>=1; i--) {
			p1 = p    + p1 * *x;
			p  = C[i] + p  * *x;
		}
		if (Abs(p1) < eps) {
			// polpol: derivative should not vanish
			*x = xorg;
			return 1;
		}
		*x -=  p / p1;
		// print *,'cycle=',icycle,' x=',x,abs(x-x0),eps
		if ( Abs(*x - x0) < xabs) return 0;
	}
	return 2;
} // polpol

/** kahanSum	(from wikipedia)
 * In numerical analysis, the Kahan summation algorithm (also known as
 * compensated summation) significantly reduces the numerical error in the
 * total obtained by adding a sequence  of finite precision floating point
 * numbers, compared to the obvious approach. This is done by keeping a
 * separate running compensation (a variable to accumulate small errors).
 *
 * In particular, simply summing n numbers in sequence has a worst-case error
 * that grows proportional to n, and a root mean square error that grows as
 * sqrt(n) for random inputs (the roundoff errors form a random walk).[1] With
 * compensated summation, the worst-case error bound is independent of n, so a
 * large number of values can be summed with an error that only depends on the
 * floating-point precision.[1]

 * The algorithm is attributed to William Kahan[2]. Similar, earlier techniques
 * are, for example, Bresenham's line algorithm, keeping tab of the accumulated
 * error in integer operations (although first documented around the same
 * time[3]) and the Delta-sigma modulation[4] (integrating, not just summing
 * the error).
 *
 * @param n	number of items to add
 * @param input	array of double precision numbers to add
 * @return the sum of input
 */
double kahanSum(const int n, const double *input)
{
	double sum = input[0];
	double c = 0.0;		//A running compensation for lost low-order bits.
	for (int i=1; i<n; i++) {
		double y = input[i] - c;	// So far, so good: c is zero.
		double t = sum + y;	// Alas, sum is big, y small,
					// so low-order digits of y are lost.
		c = (t - sum) - y;	// (t - sum) recovers the
					// high-order part of y;
					// subtracting y recovers -(low part of y)
		sum = t;		// Algebraically, c should always be zero.
					// Beware eagerly optimising compilers!
		// Next time around, the lost low part will be added to y in a
		// fresh attempt.
	}
	return sum;
} // kahanSum

#if 0
// From http://www.cprogramming.com/tutorial/floating_point/understanding_floating_point_representation.html
// It's worthwhile to take a look at ieee754.h if it exists on your system.
// How does this function work? First, if the numbers are totally equal, then we
// know the answer right away. Next, we check the exponents and signs; if the
// numbers have different exponents then they cannot be equal (because of the 1.m
// representation). In the rest of the function, we use the old fabs-less-than
// idiom, but validating its assumption by actually setting the exponents of the
// two numbers to zero. We compare against 0.5*10^-sigfigs, which provides for
// rounding. Imagine sigfigs==2 then we would want 0.01 to equal 0.008 and 0.012,
// but not 0.018. Hence we compare the difference to .005, which is 0.5*10^-2.
//
// In production code, you would want to avoid computing sig_mag every time,
// probably by doing so in another function that sets the desired precision. Also
// note that I can't make any guarantees regarding the correctness or performance
// of this code; like any code I present it is intended only for instructional
// purposes.
//
// You might have noticed that this function does not do a great job comparing
// numbers to zero (since zero's representation is a special case, zero is a
// special case here as well). Experience has shown that numbers which "should be"
// zero do not tend to have a zero exponent¿they are not those special extra-small
// numbers. Instead they are just some random value around epsilon, like
// 6.12303e-17, which is what I get for the double-precision cosine of pi/2 (this
// is probably because algorithm designers often aim to get their results within
// epsilon). In this case it is better to compare the absolute value of a number
// to sig_mag/2 directly, without manipulating exponents.
//
// This function is also deficient for comparing those extra-small numbers. We
// would have to handle their different representation, in which significant
// digits begin wherever the first set bit is. The function will typically
// conclude that two such numbers are equal (since, effectively, it will just be
// comparing leading zeros), unless you ask for enough digits to reach all the way
// out to whatever bits are set.

int flt_equals(float a, float b, int sigfigs)
{
	union ieee754_float *pa, *pb;
	unsigned int aexp, bexp;
	float sig_mag;

	if (a == b)
		return 1;
	pa = (union ieee754_float*)&a;
	pb = (union ieee754_float*)&b;
	aexp = pa->ieee.exponent;
	bexp = pb->ieee.exponent;
	if (aexp != bexp || pa->ieee.negative != pb->ieee.negative)
		return 0;
	pa->ieee.exponent = pb->ieee.exponent = IEEE754_FLOAT_BIAS;
	sig_mag = pow(10, -(float)sigfigs);
	if (fabs(a-b) < sig_mag/2)
		return 1;
	return 0;
}

// The standard C library never seems to do quite what you want for printing
// floats. If you want scientific notation, you can use "%e", but then 0 prints as
// 0.000000e+00. Or you can use %f, but then large numbers yield long strings of
// digits rather than the scientific notation you'd prefer.
//
// As a parting gift, here's a routine that prints real numbers a little more
// nicely, automatically adjusting format codes depending on what kind of number
// you give it. You can specify how big or small a number can get before moving to
// scientific notation, and you can still specify field widths as in the usual
// "%n.nf" format.

#define LOG2_10 3.321928095

#define flt_zero(x) (fabs(x) < EPSILON)

int max_digs_rt = 3;  /* maximum # of 0's right of decimal before using
			scientific notation */
int max_digs_lf = 5;  /* max # of digits left of decimal */

void print_real(double r, int width, int dec)
{
	int mag;
	double fpart, temp;
	char format[8];
	char num_format[3] = {'l',0,0};
	union ieee754_double *dl;

	dl = (union ieee754_double*)&r;
	mag = (dl->ieee.exponent - IEEE754_DOUBLE_BIAS) / LOG2_10;
	if (r == 0)
		mag = 0;
	if ((mag > max_digs_lf-1) || (mag < -max_digs_rt)) {
		num_format[1] = 'e';
		temp = r/pow(10, mag);      /* see if number will have a decimal */
		fpart = temp - floor(temp); /* when written in scientific notation */
	}
	else {
		num_format[1] = 'f';
		fpart = r - floor(r);
	}
	if (flt_zero(fpart))
		dec = 0;
	if (width == 0) {
		snprintf(format, 8, "%%.%d%s", dec, num_format);
	}
	else {
		snprintf(format, 8, "%%%d.%d%s", width, dec, num_format);
	}
	printf(format, r);
}

#endif

/**
bool AlmostEqual2sComplement(float A, float B, int maxUlps)
{
	// Make sure maxUlps is non-negative and small enough that the
	// default NAN won't compare as equal to anything.
	assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
	int aInt = *(int*)&A;
	// Make aInt lexicographically ordered as a twos-complement int
	if (aInt < 0) aInt = 0x80000000 - aInt;
	// Make bInt lexicographically ordered as a twos-complement int
	int bInt = *(int*)&B;
	if (bInt < 0) bInt = 0x80000000 - bInt;
	int intDiff = abs(aInt - bInt);
	if (intDiff <= maxUlps) return true;
	return false;
}
*/

/** Sampling functions - van der Corput series
 * Chris Theis [1/10/2005]
 */
float VanDerCorput(unsigned int n, unsigned int Scramble)
{
	n = (n << 16) | (n >> 16);
	n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
	n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
	n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
	n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
	n ^= Scramble;
	return (float)n / (float)0x100000000LL;
} // VanDerCorput

/** Sampling functions - Sobol series
 * Chris Theis [1/10/2005]
 */
float Sobol2(unsigned int n, unsigned int Scramble)
{
	for (unsigned int v = 1 << 31; n != 0; n >>= 1, v ^= v >> 1)
		if (n & 0x1) Scramble ^= v;
	return (float)Scramble / (float)0x100000000LL;
} // Sobol2

/** fpetrap - trap floating point exceptions if needed */
void fpetrap()
{
#if !defined(Darwin) && !defined(__CYGWIN__)
	::feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW|FE_UNDERFLOW);
#endif
} // fpetrap
