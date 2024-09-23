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
 *
 * This algorithm is an adaptation in C++ of RPOLY from Jenkins
 * algorithm for root finding.

 * https://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm
 * http://dl.acm.org/citation.cfm?id=355643&coll=ACM&dl=ACM
 *
 */

#include <float.h>
#include <string.h>

#include "os.h"
#include "bmath.h"
#include "polynomial.h"

using namespace std;

const double PolynomialSolver::base    = 2;
const double PolynomialSolver::eta     = DBL_EPSILON;
const double PolynomialSolver::infin   = DBL_MAX;
const double PolynomialSolver::smalno  = DBL_MIN;
// are and mre refer to the unit error in + and *
// respectively. they are assumed to be the same as
// eta.
const double PolynomialSolver::are     = DBL_EPSILON;
const double PolynomialSolver::mre     = DBL_EPSILON;
const double PolynomialSolver::lo      = DBL_MIN / DBL_EPSILON;

/** Finds the zeros of a real polynomial
 * @param op      vector of coefficients in
 *                order of decreasing powers.
 * @param degree  degree of polynomial.
 * @param zeror   output vector of real parts of the zeros.
 * @param zeroi   output vector of imaginary parts of zeros
 * @return        logical parameter, true only if
 *                leading coefficient is zero or if it
 *                has found fewer than degree zeros.
 *                In the latter case degree is reset to
 *                the number of zeros found.
 *                To change the size of polynomials which can be
 *                solved, reset the dimensions of the arrays in the
 *                common area and in the following declarations.
 * The subroutine uses single precision calculations
 * for scaling, bounds and error calculations. All
 * calculations for the iterations are done in double
 * precision.
 */
int PolynomialSolver::roots(double* op, int degree, double* zeror, double* zeroi)
{
	double	pt[101], temp[101];

	// initialization of constants for shift rotation
	double xx   = sqrt(0.5);
	double yy   = -xx;
	const static double rot  = 94.0;
	const static double cosr = COSD(rot);
	const static double sinr = SIND(rot);
	int nz;

	n = degree;

	// algorithm fails if the leading coefficient is zero.
	if (op[0] == 0.0) return -1;

	// remove the zeros at the origin if any
	while (op[n] == 0.0) {
		int j = degree - n;
		zeror[j] = 0.0;
		zeroi[j] = 0.0;
		n--;
	}
	if (n<1) return -1;

	// make a copy of the coefficients
	for (int i=0; i<=n; i++) p[i] = op[i];		// @BNV memcpy

	do {	// Label 40
		// start the algorithm for one zero
		if (n==1) {
			zeror[degree-1] = -p[1]/p[0];
			zeroi[degree-1] = 0.0;
			return degree-n+1;
		} else
		if (n==2) {
			// calculate the final zero or pair of zeros
			iquadratic(p[0], p[1], p[2],
				&zeror[degree-2], &zeroi[degree-2],
				&zeror[degree-1], &zeroi[degree-1]);
			return degree-n+2;
		}

		// find largest and smallest moduli of coefficients.
		double max = 0.0;
		double min = infin;
		for (int i=0; i<n; i++) {
			double x = Abs(p[i]);
			if (x>max) max = x;
			if (x!=0.0 && x<min) min = x;
		}

		// scale if there are large or very small coefficients
		// computes a scale factor to multiply the
		// coefficients of the polynomial. the scaling is done
		// to avoid overflow and to avoid undetected underflow
		// interfering with the convergence criterion.
		// the factor is a power of the base
		double sc = lo/min;

		if ((sc<=1.0 && max>=10.0) ||
		    (sc> 1.0 && infin/sc >= max) ||
		    (infin/sc >= max && max >= 10.0)) {
			if (sc==0.0) sc = smalno;
			int l = (int)(log(sc)/log(base) + 0.5);
			double factor = pow(base,(double)l);
			if (factor != 1.0) {
				for (int i=0; i<=n; i++)
					p[i] = factor*p[i];
			}
		}

		// compute lower bound on moduli of zeros.
		for (int i=0; i<=n; i++)
			pt[i] = Abs(p[i]);
		pt[n] = -pt[n];

		// compute upper estimate of bound
		double x = exp((log(-pt[n])-log(pt[0]))/double(n));
		double xm;
		if (pt[n]!=0.0) {
			// if newton step at the origin is better, use it.
			xm = -pt[n]/pt[n-1];
			if (xm<x) x = xm;
		}

		// chop the interval (0,x) until ff <= 0
		while (true) {
			xm = x*0.1;
			double ff = pt[0];
			for (int i=1; i<=n; i++)
				ff = ff*xm + pt[i];
			if (ff<=0.0) break;
			x = xm;
		}

		double dx = x;
		// do newton iteration until x converges to two
		// decimal places
		while (Abs(dx/x)>.005) {
			double ff = pt[0];
			double df = ff;
			for (int i=1; i<n; i++) {
				ff = ff*x + pt[i];
				df = df*x + ff;
			}
			ff = ff*x + pt[n];
			dx = ff/df;
			x -= dx;
		}
		double bnd = x;

		// compute the derivative as the intial k polynomial
		// and do 5 steps with no shift
		int nm1 = n-1;
		for (int i=1; i<n; i++)
			k[i] = double(n-i)*p[i]/double(n);

		k[0] = p[0];
		double aa = p[n];
		double bb = p[n-1];
		bool zerok = (k[n-1]==0.0);

		for (int jj=0; jj<5; jj++) {
			double cc = k[n-1];
			if (!zerok) {
				// use scaled form of recurrence if value of k at 0 is
				// nonzero
				double t = -aa/cc;
				for (int i=0; i<nm1; i++) {
					int j = n-i-1;
					k[j] = t*k[j-1] + p[j];
				}

				k[0] = p[0];
				zerok = Abs(k[n-1])<=Abs(bb)*eta*10.0;
			} else {
				// use unscaled form of recurrence
				for (int i=0; i<nm1; i++) {	// BNV memcpy to shift?
					int j = n-i-1;
					k[j] = k[j-1];
				}
				k[1] = 0.0;
				zerok = (k[n-1]==0.0);
			}
		}

		// save k for restarts with new shifts
		for (int i=0; i<n; i++)		// BNV memcpy
			temp[i] = k[i];

		// loop to select the quadratic  corresponding to each
		// new shift
		for (int cnt=1; cnt<=20; cnt++) {
			// quadratic corresponds to a double shift to a
			// non-real point and its complex conjugate. the point
			// has modulus bnd and amplitude rotated by 94 degrees
			// from the previous shift
			double	xxx = cosr*xx - sinr*yy;
				yy  = sinr*xx + cosr*yy;
			xx = xxx;
			sr = bnd*xx;
			si = bnd*yy;
			u = -2.0*sr;
			v = bnd;

			// second stage calculation, fixed quadratic
			nz = fxshfr(20*cnt);

//cout << endl << "Cnt=" << cnt << " nz=" << nz << endl;
//for (int iii=0; iii<=n; iii++)
//	cout << "   " << iii+1 << ": " << p[iii] << " " << qp[iii]
//	     << " " << k[iii] << " " << qk[iii]
//	     << " " << svk[iii] << endl;

			if (nz!=0) { //goto 260
				// the second stage jumps directly to one of the third
				// stage iterations and returns here if successful.
				// deflate the polynomial, store the zero or zeros and
				// return to the main algorithm.
				int j = degree - n;
				zeror[j] = szr;
				zeroi[j] = szi;
				n -= nz;
				for (int i=0; i<=n; i++)	// @BNV memcpy? WARNING on limit
					p[i] = qp[i];

				if (nz!=1) {
					zeror[j+1] = lzr;
					zeroi[j+1] = lzi;
				}
				break;
			} else {
				// if the iteration is unsuccessful another quadratic
				// is chosen after restoring k
				for (int i=0; i<n; i++)		// BNV memcpy
					k[i] = temp[i];
			}
		}
	} while (nz);

	// return with failure if no convergence with 20
	// shifts
	return degree - n;
} // jtpolynomialsolver

/** computes up to l2 fixed shift k-polynomials,
 * testing for convergence in the linear or quadratic
 * case. initiates one of the variable shift
 * iterations and returns with the number of zeros
 * found.
 * @param l2	limit of fixed shift steps
 * @return (nz) number of zeros found
 */
int PolynomialSolver::fxshfr(int l2)
{
	double betav = .25;
	double betas = .25;
	double ots = 0.0, otv = 0.0;
	double ovv = v;
	double oss = sr;
	double ui, vi;

	int nz    = 0;
	int type  = 0;
	int iflag = 0;

	// evaluate polynomial by synthetic division
	quadsd(n, u, v, p, qp, &a, &b);
	calcsc(&type);
	for (int j=0; j<l2; j++) {
		// calculate next k polynomial and estimate v
		nextk(type);
		calcsc(&type);
		newest(type, &ui, &vi);
		double vv = vi;

		// estimate s
		double ss = 0.0;
		if (k[n-1]!=0.0) ss = -p[n]/k[n-1];
		double tv = 1.0;
		double ts = 1.0;
		if (j==0 || type==3) {
			ovv = vv;
			oss = ss;
			otv = tv;
			ots = ts;
			continue;
		}

		// compute relative measures of convergence of s and v
		// sequences
		if (vv!=0.0) tv = Abs((vv-ovv)/vv);
		if (ss!=0.0) ts = Abs((ss-oss)/ss);

		// if decreasing, multiply two most recent
		// convergence measures
		double tvv = 1.0;
		if (tv<otv) tvv = tv*otv;
		double tss = 1.0;
		if (ts<ots) tss = ts*ots;

		// compare with convergence criteria
		bool vpass = tvv<betav;
		bool spass = tss<betas;
		if (!(spass || vpass)) {
			ovv = vv;
			oss = ss;
			otv = tv;
			ots = ts;
			continue;
		}

		// at least one sequence has passed the convergence
		// test. store variables before iterating
		double svu = u;
		double svv = v;
		for (int i=0; i<n; i++)		// @BNV memcpy
			svk[i] = k[i];

		double s = ss;

		// choose iteration according to the fastest
		// converging sequence
		bool vtry = false;
		bool stry = false;
		if (spass && (!vpass || tss<tvv)) goto L40;

L20:		nz = quadit(&ui, &vi);
		if (nz>0) return nz;

		// quadratic iteration has failed. flag that it has
		// been tried and decrease the convergence criterion.
		vtry   = true;
		betav /= 4.0;

		// try linear iteration if it has not been tried and
		// the s sequence is converging
		if (stry || !spass) goto L50;

		for (int i=0; i<n; i++)		// @BNV convert to memcpy
			k[i] = svk[i];

L40:		nz = realit(&s, &iflag);
		if (nz>0) return nz;

		// linear iteration has failed. flag that it has been
		// tried and decrease the convergence criterion
		stry   = true;
		betas /= 4.0;
		if (iflag==0) goto L50;

		// if linear iteration signals an almost double real
		// zero attempt quadratic interation
		ui = -(s+s);
		vi = s*s;
		goto L20;

// restore variables
L50:		u = svu;
		v = svv;
		for (int i=0; i<n; i++)		// @BNV memcpy
			k[i] = svk[i];

		// try quadratic iteration if it has not been tried
		// and the v sequence is converging
		if (vpass && !vtry) goto L20;

		// recompute qp and scalar values to continue the
		// second stage
		quadsd(n, u, v, p, qp, &a, &b);
		calcsc(&type);
	}
	return nz;
} // fxshfr

/** variable-shift k-polynomial iteration for a
 * quadratic factor converges only if the zeros are
 * equimodular or nearly so.
 * @param uu,vv coefficients of starting quadratic
 * @return (nz) number of zero found
 */
int PolynomialSolver::quadit(double* uu, double* vv)
{
	int nz = 0;
	bool tried = false;
	u = *uu;
	v = *vv;
	int j = 0;
	int type = 0;
	double relstp = 0.0;
	double omp = 0.0;

	// main loop
	while (true) {
		iquadratic(1.0, u, v, &szr, &szi, &lzr, &lzi);

		// return if roots of the quadratic are real and not
		// close to multiple or nearly equal and  of opposite sign
		if (Abs(Abs(szr)-Abs(lzr)) > .01*Abs(lzr)) return nz;

		// evaluate polynomial by quadratic synthetic division
		quadsd(n, u, v, p, qp, &a, &b);
		double mp = Abs(a-szr*b) + Abs(szi*b);

		// compute a rigorous  bound on the rounding error in
		// evaluting p
		double zm = sqrt(Abs(v));
		double ee = 2.*Abs(qp[0]);
		double t = -szr*b;
		for (int i=1; i<n; i++)
			ee = ee*zm + Abs(qp[i]);

		ee = ee*zm + Abs(a+t);
		ee = (5.0*mre+4.0*are)*ee - (5.0*mre+2.0*are)
			* (Abs(a+t)+Abs(b)*zm)
			+ 2.0*are*Abs(t);

		// iteration has converged sufficiently if the
		// polynomial value is less than 20 times this bound
		if (mp<=20.*ee)
			return 2;
		j++;

		// stop iteration after 20 steps
		if (j>20) return nz;
		if (j>=2 && !(relstp>.01 || mp<omp || tried)) {
			// a cluster appears to be stalling the convergence.
			// five fixed shift steps are taken with a u,v close
			// to the cluster
			if (relstp<eta) relstp = eta;
			relstp = sqrt(relstp);
			u -= u*relstp;
			v += v*relstp;
			quadsd(n, u, v, p, qp, &a, &b);
			for (int i=0; i<5; i++) {
				calcsc(&type);
				nextk(type);
			}
			tried = true;
			j = 0;
		}
		omp = mp;

		// calculate next k polynomial and new u and v
		calcsc(&type);
		nextk(type);
		calcsc(&type);
		double ui,vi;
		newest(type, &ui, &vi);

		// if vi is zero the iteration is not converging
		if (vi==0.0) return nz;
		relstp = Abs((vi-v)/vi);
		u = ui;
		v = vi;
	}
} // quadit

/** realit
 * variable-shift h polynomial iteration for a real zero.
 * @param sss	starting iterate
 * @param iflag	flag to indicate a pair of zeros near real axis.
 * @return (nz) number of zero found
 */
int PolynomialSolver::realit(double* sss, int* iflag)
{
	double t = 0.0;
	double s = *sss;
	double omp = 0.0;

	int nz = 0;
	int j  = 0;

	*iflag = 0;

	// main loop
	while (true) {
		double pv = p[0];
		// evaluate p at s
		qp[0] = pv;

		for (int i=1; i<=n; i++) {
			pv = pv*s + p[i];
			qp[i] = pv;
		}

		double mp = Abs(pv);

		// compute a rigorous bound on the error in evaluating p
		double ms = Abs(s);
		double ee = (mre/(are+mre))*Abs(qp[0]);
		for (int i=1; i<=n; i++)
			ee = ee*ms + Abs(qp[i]);

		// iteration has converged sufficiently if the
		// polynomial value is less than 20 times this bound
		if (mp <= 20.*((are+mre)*ee-mre*mp)) {
			nz = 1;
			szr = s;
			szi = 0.;
			return nz;
		}
		j++;

		// stop iteration after 10 steps
		if (j>10) return nz;
		if (j>=2 && (Abs(t)<=.001*Abs(s-t) && mp>omp)) {
			// a cluster of zeros near the real axis has been
			// encountered return with iflag set to initiate a
			// quadratic iteration
			*iflag = 1;
			*sss   = s;
			return nz;
		}

		// return if the polynomial value has increased significantly
		omp = mp;

		// compute t, the next polynomial, and the new iterate
		double kv = k[0];
		qk[0] = kv;
		for (int i=1; i<n; i++) {
			kv = kv*s + k[i];
			qk[i] = kv;
		}

		if (Abs(kv)>Abs(k[n-1])*10.*eta) {
			// use the scaled form of the recurrence if the value
			// of k at s is nonzero
			t    = -pv/kv;
			k[0] = qp[0];
			for (int i=1; i<n; i++)
				k[i] = t*qk[i-1] + qp[i];
		} else {
			// use unscaled form
			k[0] = 0.0;
			for (int i=1; i<n; i++)
				k[i] = qk[i-1];
		}

		kv = k[0];
		for (int i=1; i<n; i++)
			kv = kv*s + k[i];

		t = 0.;
		if (Abs(kv) > Abs(k[n-1])*10.*eta) t = -pv/kv;
		s += t;
	}
} // realit

/** this routine calculates scalar quantities used to
 * compute the next k polynomial and new estimates of
 * the quadratic coefficients.
 * @param type integer variable set here indicating how the
 *             calculations are normalized to avoid overflow
 */
void PolynomialSolver::calcsc(int* type)
{
	// synthetic division of k by the quadratic 1,u,v
	quadsd(n-1, u, v, k, qk, &c, &d);
	if (Abs(c) <= Abs(k[n-1])*100.*eta) {
		if (Abs(d) <= Abs(k[n-2])*100.*eta) {
			// type=3 indicates the quadratic is almost a factor of k
			*type = 3;
			return;
		}
	}
	if (Abs(d)<Abs(c)) {
		// type=1 indicates that all formulas are divided by c
		*type = 1;
		e  = a/c;
		f  = d/c;
		g  = u*e;
		h  = v*b;
		a3 = a*e + (h/c+g)*b;
		a1 = b - a*(d/c);
		a7 = a + g*d + h*f;
      } else {
		// type=2 indicates that all formulas are divided by d
		*type = 2;
		e  = a/d;
		f  = c/d;
		g  = u*b;
		h  = v*b;
		a3 = (a+g)*e + h*(b/d);
		a1 = b*f - a;
		a7 = (f+u)*a + h;
	}
} // calcsc

/* computes the next k polynomials using scalars
 * computed in calcsc
 */
void PolynomialSolver::nextk(int type)
{
	if (type==3) {
		// use unscaled form of the recurrence if type is 3
		k[0] = 0.;
		k[1] = 0.;
		for (int i=2; i<n; i++)
			k[i] = qk[i-2];
	} else {
		double temp = (type==1)? b : a;
		if (Abs(a1) <= Abs(temp)*eta*10.) {
			// if a1 is nearly zero then use a special form of the
			// recurrence
			k[0] = 0.;
			k[1] = -a7*qp[0];
			for (int i=2; i<n; i++)
				k[i] = a3*qk[i-2] - a7*qp[i-1];
		} else {
			// use scaled form of the recurrence
			a7 = a7/a1;
			a3 = a3/a1;
			k[0] = qp[0];
			k[1] = qp[1] - a7*qp[0];
			for (int i=2; i<n; i++)
				k[i] = a3*qk[i-2] - a7*qp[i-1] + qp[i];
		}
	}
} // nextk

/** compute new estimates of the quadratic coefficients
 * using the scalars computed in calcsc.
 */
void PolynomialSolver::newest(int type, double *uu, double *vv)
{
	// use formulas appropriate to setting of type.
	if (type==3) {
		// if type=3 the quadratic is zeroed
		*uu = *vv = 0.0;
		return;
	}

	double a4, a5;
	if (type==2) {
		a4 = (a+g)*f + h;
		a5 = (f+u)*c + v*d;
	} else {
		a4 = a + u*b + h*f;
		a5 = c + (u+v*f)*d;
	}

	// evaluate new quadratic coefficients.
	double b1 = -k[n-1]/p[n];
	double b2 = -(k[n-2]+b1*p[n-1])/p[n];
	double c1 = v*b2*a1;
	double c2 = b1*a7;
	double c3 = b1*b1*a3;
	double c4 = c1 - c2 - c3;
	double temp = a5 + b1*a4 - c4;
	if (temp==0.0) {
		*uu = *vv = 0.0;
		return;
	}
	*uu = u - (u*(c3+c2)+v*(b1*a1+b2*a7))/temp;
	*vv = v*(1.+c4/temp);
} // newest

/** divides p by the quadratic  1,u,v  placing the
 *  quotient in q and the remainder in a,b
 */
void PolynomialSolver::quadsd(int nn, const double uu, const double vv, double pp[], double qq[], double* aa, double* bb)
{
	*bb   = pp[0];
	qq[0] = *bb;
	*aa   = pp[1] - uu*(*bb);
	qq[1] = *aa;
	for (int i=2; i<=nn; i++) {
		double cc = pp[i] - uu*(*aa) - vv*(*bb);
		qq[i] = cc;
		*bb   = *aa;
		*aa   = cc;
	}
} // quadsd

/////////////////////////////////////////////////////////////////////////////////////////
#if 0
# include <iomanip>
# include <stdlib.h>
# include <timer.h>
# include "G4JTPolynomialSolver.hh"
/** main */
int main(int ac, char* av[])
{
	PolynomialSolver poly;
	double sr,si,lr,li;

	double a=1.0;
	double b=-5.0;
	double c=1.0;
	int times = 0;
	if (ac>=2) times = atoi(av[1]);
	if (times<=0) times=1;

	double zeror[101];
	double zeroi[101];
//	double op[] = { 1.0, -6.0, 11.0, -6.0};
	double op[] = { 1.0, -2.0, -13.0, 14.0, 24.0};
//	double op[] = { 1.0, 1.0, -13.0, -25.0, -12.0};
//	double op[] = { 1.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 0.9604000000000001, -5.997600000000001, 13.951750054511718, -14.326264455924333, 5.474214401412618, -4.0, 5.0, 1.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 0.9604000000000001, -5.997600000000001, 13.951750054511718, -14.326264455924333, 5.474214401412618, -4.0, 5.0 };
	int nz;
	int degree = sizeof(op) / sizeof(double) - 1;
	cout << "degree=" << degree << endl;

#define SOL

#if 1
	Timer timer;
	timer.start();
	for (int i=0; i<times; i++)
		nz = poly.roots(op, degree, zeror, zeroi);
	timer.stop();
	cout << timer << endl;

	//cout << setprecision(16);
#ifdef SOL
	cout << "Solutions= " << nz << endl;
	for (int i=0; i<nz; i++)
		cout << i << ": " << zeror[i] << " + " << zeroi[i] << " i" << endl;
#endif
#endif

#if 1
	timer.reset();
	timer.start();
	JTPolynomialSolver solver;
	for (int i=0; i<times; i++)
		nz = solver.FindRoots(op, degree, zeror, zeroi);
#ifdef SOL
	cout << "G4 Solutions= " << nz << endl;
	for (int i=0; i<nz; i++)
		cout << i << ": " << zeror[i] << " + " << zeroi[i] << " i" << endl;
#endif
	timer.stop();
	cout << timer << endl;
#endif

#if 1
	timer.reset();
	timer.start();
	for (int i=0; i<times; i++)
		nz = quartic(op[1]/op[0], op[2]/op[0], op[3]/op[0], op[4]/op[0], zeror, 1e-10, 4);
	timer.stop();
	cout << timer << endl;
#ifdef SOL
	cout << "Quartic Solutions= " << nz << endl;
	for (int i=0; i<nz; i++)
		cout << i << ": " << zeror[i] << endl;
#endif
#endif

#if 0
	iquadratic(a,b,c,&sr,&si,&lr,&li);
	cout << "x1= " << sr << " + " << si << "i" << endl;
	cout << "x2= " << lr << " + " << li << "i" << endl;
	double x1=0.0, x2=0.0;
	cout << "n=" << quadratic(b/a,c/a,&x1,&x2,1e-10) << endl;
	cout << "x1= " << x1 << endl;
	cout << "x2= " << x2 << endl;
#endif

	return 0;
} // main
#endif
