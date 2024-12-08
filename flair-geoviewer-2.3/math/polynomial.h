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
 * Date:	10-Oct-2016
 */

#ifndef __POLYNOMIAL_H
#define __POLYNOMIAL_H

// ONLY FOR DEBUGGING
#ifdef _DEBUG
#	include <iostream>
#	include <strings.h>
#endif

/* =============================== Polynomial =============================== */
class Polynomial {
private:
	int	n;
	double*	c;

public:
	Polynomial() : n(0), c(NULL) {}
	~Polynomial() { if (c) delete [] c; }

//const	double	operator[](int i)	const	{ return c[i]; }
//	double&	operator[](int i)		{ return c[i]; }

	double	operator()(const double x) const {
			double s = c[0]*x;
			for (int i=0; i<n; i++)
				s = s*x + c[i];
			return s;
		}
	void	set(int i, double x)	{ c[i] = x; }
//	double	derivative(const double x) const { return 0.0; }
//	double	integrate(const double a, const double b);

//	int	roots();
}; // class Polynomial

/* ============================ PolynomialSolver ============================ *
 * Jenkins-Traub alogirthm
 */
class PolynomialSolver {
private:
	double	p[101];
	double	qp[101];
	double	k[101];
	double	qk[101];
	double	svk[101];

	double	sr, si;
	double  u, v;
	double	a, b, c, d;
	double	a1, a2, a3, a6, a7;
	double	e, f, g;
	double	h, szr, szi, lzr, lzi;
	int	n;

public:
	PolynomialSolver() :
			sr(0.), si(0.),
			u(0.), v(0.),
			a(0.), b(0.), c(0.), d(0.),
			a1(0.), a2(0.), a3(0.), a6(0.), a7(0.),
			e(0.), f(0.), g(0.),
			h(0.), szr(0.), szi(0.), lzr(0.), lzi(0.),
			n(0)
			{
#ifdef _DEBUG
				bzero(p,   sizeof(p));
				bzero(qp,  sizeof(qp));
				bzero(k,   sizeof(k));
				bzero(qk,  sizeof(qk));
				bzero(svk, sizeof(svk));
#endif
			}

	int	roots(double* op, int degree, double* zeror, double* zeroi);

private:
	int	fxshfr(int l2);
	int	quadit(double* uu, double* vv);
	int	realit(double* sss, int* iflag);
	void	calcsc(int* type);
	void	nextk(int type);
	void	newest(int type, double* uu, double* vv);
	void	quadsd(int nn, double uu, double vv, double pp[], double qq[], double* aa, double* bb);

public:
	// the following statements set machine constants used
	// in various parts of the program. the meaning of the
	// four constants are...
	static	const double	base;	// the base of the floating-point number
	static	const double	eta;	// the maximum relative representation error
					// which can be described as the smallest
					// positive floating point number such that
					// 1+eta is greater than 1.
	static const double	infin;	// the largest floating-point number.
	static const double	smalno; // the smallest positive floating-point number
					// if the exponent range differs in single and
					// double precision then smalno and infin
					// should indicate the smaller range.
	static const double	are;
	static const double	mre;
	static const double	lo;
}; // class PolynomialSolver

#endif
