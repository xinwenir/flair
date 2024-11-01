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

#ifndef __BMATH_H
#define __BMATH_H

#include <math.h>
#include <iostream>
#include <stdint.h>

#include "os.h"

#define PRINTDBL(val) do \
{\
	double _my_val = (val);\
	printf("%s -> %30.30f (%g) %llx\n", #val, _my_val, _my_val, *(unsigned long long *)&_my_val);\
} while(0)

#if defined(Darwin) || defined(ANDROID) || defined(CYGWIN)
// Mod by VB - 11/01/11
//   MacOS lacks few math functions
/** pow10 */
inline double pow10(const double x)                 { return pow(10.0,x);  }
/** pow10f */
inline float  pow10f(const float x)                 { return powf(10.0,x); }
inline void   sincos(const double x, double *s, double *c){ *s = sin(x); *c = cos(x);}
#endif

/** byte multiplication assuming both bytes represent a range of 0.0 .. 1.0 decimal
 * @param a, b
 * @return a * b / 256
 */
inline uint8_t byteMul(uint8_t a, uint8_t b)
{
	  return (((int)a + 1) * (int)b) >> 8;
} /* byteMul */

/** Check whether the passed integer is a power of 2 */
inline bool isPowerOf2(int v)
{
	return (v & (v - 1)) == 0;
} /* isPowerOf2 */

/** Round up to the next power of 2 */
inline unsigned int roundUpPow2(unsigned int v)
{
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return v+1;
} /* roundUpPow2 */

#if 0
/** Fast inverse square root (https://en.wikipedia.org/wiki/Fast_inverse_square_root) 
 * @param number	numner to return 1/sqrt(number)
 * @return a reasonable approximation of the inverse square root
 */
inline float Q_rsqrt( float number )
{
	int i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y  = number;
	i  = *(int *)&y;			// evil floating point bit level hacking
	i  = 0x5f3759df - ( i >> 1 );		// what the fuck?
	y  = *(float*)(&i);
	y  = y * (threehalfs - x2*y*y);		// 1st iteration
	y  = y * (threehalfs - x2*y*y);		// 2nd iteration, this can be removed
	return y;
} // Q_rsqrt
#endif

/**
 * Numerically robust quadratic equation: x^2 + bx + c = 0.
 * @param a,b,c	coefficients
 * @param x1,x2	two real solutions.
 *              In case of no real solution then -b/2a is returned
 * @param acc	accuracy of comparison operations
 * @return the number of real solutions
 */
inline int quadratic(const double b, const double c,
		double *x1, double *x2, const double acc)
{
#if 0
	double D = b*b - 4.0*c;
	if (D <= 0.0) {
		*x1 = *x2 = -0.5*b;	// Always return this as a solution!!!
		if (D >= -acc*(b*b+Abs(c)))
			return 1;
		else
			return 0;
	} else {
		double bD;
		if (b>0.0)
			bD = -b - sqrt(D);
		else
			bD = -b + sqrt(D);
		*x1 = 0.5 * bD;
		*x2 = 2.0 * c / bD;	// Vieta's Formula
		return 2;
	}
#endif
	double b2 = b*b;
	double c4 = 4.0*c;
	double D = b2 - c4;
	double prec = acc*Max(b2, Abs(c4));
	if (D <= prec) {
		*x1 = *x2 = -0.5*b;	// Always return this as a solution!!!
		if (D < -prec) return 0;
		return 1;
	} else {
		double bD;
		if (b>0.0)
			bD = -b - sqrt(D);
		else
			bD = -b + sqrt(D);
		*x1 = 0.5 * bD;
		*x2 = 2.0 * c / bD;	// Vieta's Formula
		return 2;
	}
} /* quadratic */

/** calculate the zeros of the quadratic a*z**2+b*z+c.
 * the quadratic formula, modified to avoid
 * overflow, is used to find the larger zero if the
 * zeros are real and both zeros are complex.
 * the smaller real zero is found directly from the
 * product of the zeros c/a.
 */
inline void iquadratic(double a, double b, double c,
				double* xr, double* xi,
				double* yr, double* yi)
{
	if (a==0.0) {
		*xr = (b!=0.0)? -c/b : 0.0;
		*yr = *xi = *yi = 0.0;
		return;
	}

	if (c==0.0) {
		*yr = -b/a;
		*xr = *xi = *yi = 0.0;
		return;
	}

	// compute discriminant avoiding overflow
	double b2 = b/2.;
	double e, d;
	if (Abs(b2)<Abs(c)) {
		e = (c<0.0)? -a : a;
		e = b2*(b2/Abs(c)) - e;
		d = sqrt(Abs(e))*sqrt(Abs(c));
	} else {
		e = 1.0 - (a/b2)*(c/b2);
		d = sqrt(Abs(e))*Abs(b2);
	}

	if (e>=0.0) {
		// real zeros
		if (b2>=0.0) d = -d;
		*yr = (-b2+d)/a;
		*xr = 0.0;
		if (*yr!=0.0) *xr = (c / *yr)/a;
		*xi = *yi = 0.0;
		return;
	} else {
		// complex conjugate zeros
		*xr = -b2/a;
		*yr = *xr;
		*xi = Abs(d/a);
		*yi = -*xi;
	}
} // iquadratic

/**
 * sincos but round according to eps
 */
inline void bsincos(const double x, double *s, double *c)
{
	sincos(x, s, c);
	if (*s==1.0 || *s==-1.0) *c = 0.0;
	else
	if (*c==1.0 || *c==-1.0) *s = 0.0;
} /* bsincos */

unsigned isqrt(unsigned n);
int	cubic(	const double a, const double b, const double c, double x[],
		const double acc, const int iter);
int	quartic(const double a, const double b, const double c, const double d,
		double x[], const double acc, const int iter);
int	polpol(const double C[], const int n, double *x, const double eps);
double	kahanSum(const int n, const double *input);
float	VanDerCorput(unsigned int n, unsigned int Scramble);
float	Sobol2(unsigned int n, unsigned int Scramble);
void	fpetrap();

/* ================================= KahanSum =============================== *
 * Kahan sum as a class for easier implementation
 */
class KahanSum {
public:
	double sum;
	double c;
public:
	KahanSum(const double n=0.0) : sum(n), c(0.0) {}
	void	add(const double n) {
			double y = n - c;	// So far, so good: c is zero.
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
	void	operator()(const double n)	{ add(n); }
	double	operator()()		const	{ return sum; }
}; // KahanSum

#endif
