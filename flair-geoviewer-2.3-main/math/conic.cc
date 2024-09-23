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
 * Date:	04-Feb-2010
 */

#include <math.h>
#include <ostream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <iostream>

#include "os.h"
#include "eps.h"
#include "bmath.h"
#include "conic.h"
#include "format.h"
#include "matrix3.h"
#include "bstring.h"

using namespace std;

const double Conic::CONICPREC = 1e-14;
const double CONICSMALL = 10.0*Conic::CONICPREC;
const double TOLERANCE  = 1.0E-6;
#define ISZEROH(a)	Eq0(a,1E-11)

const char *Conic::_typeStr[] = {
	"DEGENERATE",
	"POINT",
	"LINE",
	"LINES",
	"ELLIPSE",
	"PARABOLA",
	"HYPERBOLA"
};

/** set conic parameters */
void Conic::set(const double aa, const double ah, const double ab,
		const double ag, const double af, const double ac)
{
	/* copy parameters */
	a = aa;
	h = ah;
	b = ab;
	g = ag;
	f = af;
	c = ac;

	/* reset parametric coefficients */
	c1 = c2 = c3 = 0.0;
	c4 = c5 = c6 = 0.0;

	/* normalized coefficients */
	double xa, xb, xg, xf;

	/* make zero small values... */
	if (ISZEROH(h)) {
		/* calculate invariants */
		D  = a*(b*c - f*f) - b*g*g;
		I  = a + b;
		J  = a*b;
		K  = g*g + f*f;

		xa = a;
		xb = b;
		xg = g;
		xf = f;
	} else {
		// rotate to make zero the h
		double st, ct;
		double th = theta();

		bsincos(th, &st, &ct);
		//assert((b-a)*st*ct + h*(ct+st)*(ct-st) < 1e-15);

		// Try to normalize the conic before checking the parameters for the type
		xa = a*ct*ct + 2.0*h*st*ct + b*st*st;
		xb = a*st*st - 2.0*h*st*ct + b*ct*ct;
		xg =  g*ct + f*st;
		xf = -g*st + f*ct;

		// calculate the *new* invariants
		D  = xa*(xb*c - xf*xf) - xb*xg*xg;
		I  = xa + xb;
		J  = xa*xb;
		K  = xg*xg + xf*xf;
		//if (!ISZERO(xa) && !ISZERO(xb))
		if (xa!=0.0 && xb!=0.0)
			// calculate xc if zero, then it passes from center
			// therefore only intersecting lines
			if (Abs(c-xg*xg/xa-xf*xf/xb)<=CONICSMALL) D = 0.0;
	}

	/* determine and return conic type up to resolution
	 *                                                             ca-g^2
	 *   Conic                                    D     J    D/I  +bc-f^2
	 *   -------------------------------------- ----- ----- ----- --------
	 *   Real ellipse                            !0    +ve   -ve      -
	 *   Virtual ellipse                         !0    +ve   +ve      -
	 *   Hyperbola                               !0    -ve    -       -
	 *   Parabola                                !0     0     -       -
	 *   Real intersecting lines                  0    -ve    -       -
	 *   Conjugate complex intersecting lines     0    +ve    -       -
	 *   Real distinct parallel lines             0     0     -      -ve
	 *   Conjugate complex parallel lines         0     0     -      +ve
	 *   Coincident lines                         0     0     -       0
	 */
	double t = Abs(g) + Abs(f);

	_type = CONIC_NONE;
	if (Abs(a) + Abs(h) + Abs(b) <= CONICPREC*t) {
		// Single line
		if (t < CONICSMALL)
			_type = CONIC_DEGENERATE;
		else
			_type = CONIC_LINE;
	} else
	if (Abs(J) <= CONICSMALL*Sqr(I)) {
		// Parabola or Lines if J is small enough
		if (Abs(D) <= CONICSMALL*(1.0+Abs(I*K) + Abs(I*I*c))) {
			// parallel lines if invariant D is small enough
			if (Eq0(J, CONICSMALL*(1.0+Sqr(I)))) {
				if (c*(a+b) - g*g - f*f > CONICSMALL)
					_type = CONIC_DEGENERATE;
				else
					_type = CONIC_LINES;	// Parallel or Coincident lines
			} else
			if (J < 0.0)
				_type = CONIC_LINES;		// Coincident lines
			else
				_type = CONIC_DEGENERATE;	// Conjugate complex parallel lines
		} else
		// parabola if invariant J is small enough
		if (Abs(xa) < Abs(xb)) {
			if (Abs(xg)>CONICSMALL)
				_type = CONIC_PARABOLA;
			else
				_type = CONIC_DEGENERATE;
		} else {
			if (Abs(xf)>CONICSMALL)
				_type = CONIC_PARABOLA;
			else
				_type = CONIC_DEGENERATE;
		}
	} else
	if (J > 0.0) {
		// ellipse if invariant J>0
		//if (D / I < 0.0)
		if (D * I < 0.0)	// we only care for the sign
			_type = CONIC_ELLIPSE;
		else
			_type = CONIC_DEGENERATE;
	} else
	if (J < 0.0) {
		// last chance to check for lines
		// intersecting lines if constant term is small enough
		if (Abs(D)*Abs(J) < CONICSMALL)
			_type = CONIC_LINES;
		else
			_type = CONIC_HYPERBOLA;
	} else
		_type = CONIC_LINES;

	assert(_type != CONIC_NONE);
} // set

/** copy conic from source*/
void Conic::copy(const Conic& src)
{
	set(src.a, src.h, src.b, src.g, src.f, src.c);
} // copy

/**
 * @return matrix of conic
 *        /x\T  / a h g \   /x\
 *    S = |y| * | h b f | * |y| = 0
 *        \1/   \ g f c /   \1/
 */
void Conic::matrix(Matrix3 *m) const
{
	assert(m->rows()==m->columns() && m->rows()==3);
	(*m)(0,0) = a;
	(*m)(1,1) = b;
	(*m)(2,2) = c;
	(*m)(0,1) = (*m)(1,0) = h;
	(*m)(0,2) = (*m)(2,0) = g;
	(*m)(1,2) = (*m)(2,1) = f;
} // matrix

/**
 * rotate axes of conic around an angle theta
 * WARNING to rotate conic use -theta
 */
void Conic::rotate(const double th)
{
	double st, ct;
	bsincos(th, &st, &ct);

	double xa = a;
	double xb = b;
	double xh = h;
	double xg = g;
	double xf = f;

	a = xa*ct*ct + 2.0*xh*st*ct + xb*st*st;
	h = (xb-xa)*st*ct + xh*(ct+st)*(ct-st);
	b = xa*st*st - 2.0*xh*st*ct + xb*ct*ct;
	g =  xg*ct + xf*st;
	f = -xg*st + xf*ct;
} // rotate

/** @return dx,dy to eliminate the g, f components */
void Conic::translate2Origin(double* dx, double* dy) const
{
	*dx = ISZERO(a)? 0.0 : -g/a;
	*dy = ISZERO(b)? 0.0 : -f/b;
} // translate2Origin

/**
 * Translate axes of conic by dx, dy
 * WARNING: to translate the conic use -dx, -dy
 */
void Conic::translate(const double dx, const double dy)
{
	c += a*dx*dx + b*dy*dy + 2.0*(h*dx*dy +  g*dx + f*dy);
	g += a*dx + h*dy;
	f += h*dx + b*dy;
} // translate

/**
 * Return two single line conics from a degenerate one
 * XXX WRONG the sign of the second conic
 *	y^2 - 25 = 0
 *	returns:
 *		 x+5 = 0
 *		 x-5 = 0
 *	while should be
 *		 x+5 = 0
 *		-x+5 = 0
 */
void Conic::splitLines(Conic *conic1, Conic *conic2)
{
	assert(_type == CONIC_LINES);
//	double B, G, F, C;	// normalized coef
//	double x1, x2, y1, y2;

	// Check coefficients of conic
	if (ISZERO(a)) {
		//assert(!ISZERO(b));
		if (ISZEROH(h)){
			/* no rotation
			 *	(y+c1)(y+c2) = 0
			 * 0 + 0 + y^2 + 0 + (c1+c2)*y + c1*c2 = 0
			 *
			 * 2*F = c1+c2 => c2 = 2*F-c1
			 *   C = c1*c2 = c1*(2*F-c1)
			 * <=> c1^2 - 2*F*c1 + C = 0
			 */
			//assert( ISZERO(g));
			double F = f / b;
			double C = c / b;
			double x1, x2;
			quadratic(-2.0*F, C, &x1, &x2, CONICPREC);
			conic1->set(0.0, 0.0, 0.0, 0.0, 0.5, x1);
			conic2->set(0.0, 0.0, 0.0, 0.0, 0.5, x2);
		} else {
			/* with rotation
			 *	(y+c1)(x+f2*y+c2) = 0
			 * 0 + xy + f2*y^2 + c1*x + (c2+c1*f2)*y + c1*c2 = 0
			 * A=0, H=0.5
			 * B=f2
			 * 2*G = c1
			 * 2*F = c2+c1*f2
			 * C   = c1*c2
			 */
			double B = 0.5 * b / h;
			double G = g / h;
			//F = f / h;
			double C = 0.5 * c / h;
			conic1->set(0.0, 0.0, 0.0, 0.0, 0.5, G);
			conic2->set(0.0, 0.0, 0.0, 0.5, 0.5*B, C/G);
		}
	} else
	if (ISZERO(b)) {
		if (ISZEROH(h)) {
			/* no rotation, like above but with y
			 *	(x+c1)(x+c2) = 0
			 */
			//assert( ISZERO(f));
			double G = g / a;
			double C = c / a;
			double x1, x2;
			quadratic(-2.0*G, C, &x1, &x2, CONICPREC);
			conic1->set(0.0, 0.0, 0.0, 0.5, 0.0, x1);
			conic2->set(0.0, 0.0, 0.0, 0.5, 0.0, x2);
		} else {
			/* with rotation, like above with y
			 * (x+c1)(g1*x+y+c2) = 0
			 * g1*x^2 + xy + 0 + (c1*g1+c2)*x + c1*y + c1*c2 = 0
			 */
			double A = 0.5 * a / h;
			//double G = g / h;
			double F = f / h;
			double C = 0.5 * c / h;
			conic1->set(0.0, 0.0, 0.0, 0.5, 0.0, F);
			conic2->set(0.0, 0.0, 0.0, 0.5*A, 0.5, C/F);
		}
	} else
	if (ISZEROH(h)) {
		/*
		 * (x+f*y+c1)*(x-f*y+c2) = 0
		 * x^2 + 0 - f^2*y^2 + (c1+c2)*x + f*(c2-c1)*y + c1*c2 = 0
		 */
		//assert(a*b < 0.0);
		double B = b / a;
		double G = g / a;
		double F = f / a;
		double C = c / a;

		double xf;
		if (B>0.0)
			xf = 0.0;
		else
			xf = sqrt(-B);

		double x1, x2, y1, y2;
		// assume that c=0 due to the transformation applied
		if (ISZERO(F)) {
			// (x+f*y+c1)*(x-f*y+c1) = 0   (c1=c2=G)
			x1 = x2 = y1 = y2 = G;
		} else
		if (ISZERO(G)) {
			// (x+f*y-c1)*(x-f*y+c1) = 0   (-c1=c2=G)
			// c1<0 and c2>0
			y2 = x1 = - F/xf;
			x2 = y1 = -x1;

		} else {
			/*
			 * B    = -f*f
			 * 2*G  = c1+c2   =>  c2 = 2*G - c1
			 * C    = c1*c2   = c1*(2*G-c1)
			 *               => c1^2 - 2*G*c1 + C = 0
			 */
			quadratic(-2.0*G, C, &x1, &x2, CONICPREC);
			y1 = 2.0*G - x1;
			y2 = 2.0*G - x2;
		}
		// check which combination is closer to 2*F
		if (Abs(xf*(y1-x1) - 2.0*F) <= Abs(xf*(y2-x2) - 2.0*F)) {
			// normal
			conic1->set(0.0, 0.0, 0.0, 0.5, 0.5*xf, x1);
			conic2->set(0.0, 0.0, 0.0, 0.5,-0.5*xf, y1);
		} else {
			// swap
			conic1->set(0.0, 0.0, 0.0, 0.5, 0.5*xf, x2);
			conic2->set(0.0, 0.0, 0.0, 0.5,-0.5*xf, y2);
		}
	} else {
		/* generic case
		 * (x+f1*y+c1)*(x+f2*y+c2) = 0
		 * x^2 + (f1+f2)*xy + f1*f2*y^2 + (c1+c2)*x + (f1*c2+f2*c1)*y + c1*c2 = 0
		 * A   = 1
		 * 2*H = f1+f2   => f2 = 2*H-f1
		 * B   = f1*f2   = f1*(2*H-f1) => f1^2 -2*H*f1 + B = 0
		 * 2*G = c1+c2   => c2 = 2*G-c1
		 * C   = c1*c2   => c1^2 - 2*G*c1 + C = 0
		 *
		 * used for checking
		 * 2*F = f1*c2 + f2*c1
		 */
		double H = h / a;	// make A=1
		double B = b / a;
		double G = g / a;
		double F = f / a;
		double C = c / a;
		double x1, x2, y1, y2;

		quadratic(-2.0*H, B, &x1, &x2, CONICPREC);	// f1,f2

		// Substitute to 2*F to find c1,c2
		//    2*G = c1 + c2  =>  c2 = 2*G - c1
		//    2*F = f1*c2 + f2*c1 =>  c1 = (2*F - 2*G*f1) / (f2-f1)
		//                        =>  c1 = 2*(F-G*f1) / (2*H -2*f1)
		//                        =>  c1 = (F-G*f1) / (H-f1)
		if (ISZERO(G) && ISZERO(F)) {
			if (C<=0.0)
				y1 = sqrt(-C);
			else
				y1 = 0.0;
			y2 = -y1;
		} else {
			if (Eq0(J, CONICSMALL*(1.0+Abs(a*b)+Abs(h*h)))) {
				// Parallel lines
				x1 = x2 = H;
				quadratic(-2.0*G, C, &y1, &y2, CONICPREC);	// c1,2
			} else {
				// Coinciding lines
				y1 = (F - G*x1) / (H-x1);	// c1
				y2 = (F - G*x2) / (H-x2);	// c2
			}
		}
		conic1->set(0.0, 0.0, 0.0, 0.5, 0.5*x1, y1);
		conic2->set(0.0, 0.0, 0.0, 0.5, 0.5*x2, y2);

#if 0
		// check which combination f1,2, c1,2 is closer to 2*F
		if (Abs(x1*y2+x2*y1 -2.0*F) <= Abs(x1*y1+x2*y2 - 2.0*F))
#endif
	}
} // splitLines

/** addPoint to the list of points
 * @return true if the point is added
 */
bool Conic::addPoint(const Conic& conic, int *n, Vector2D pts[4], double x, double y) const
{
	if (Abs(x)>INFINITE || Abs(y)>INFINITE) return true;

	// Optimize point by performing a 2D Newton-Raphson
	double S1  = (*this)(x,y);
	double S2  = conic(x,y);
	double max = 1.0+Abs(x)+Abs(y);
	double SS  = Sqr(S1) + Sqr(S2);
	if (Abs(S1) > TOLERANCE*max || Abs(S2) > TOLERANCE*max) {
		// Improve by Newton-Raphson 2D
		for (int j=0; j<4; j++) {
			/* we try to minimize the vector
			 *	r = [x,y]T
			 *	S = [S1(r), S2(r)]T
			 *
			 * Taylor-2D expansion we have
			 *	S(r+dr) = S(r) + J*dr + O(dr^2)
			 *
			 * with J being the Jacobian
			 *	     / dS1/dx   dS1/dy \
			 *	 J = |                 |
			 *	     \ dS2/dx   dS2/dy /
			 *
			 * neglecting the O(dr^2) and setting S(r+dr)=0
			 *	J*dr = -S
			 * with
			 *	dr = inv(J) * (-S)
			 * the next approximation will be
			 *	r' = r + dr
			 */
			double S1x = 2.0*(a*x + h*y + g);
			double S1y = 2.0*(h*x + b*y + f);

			double S2x = 2.0*(conic.a*x + conic.h*y + conic.g);
			double S2y = 2.0*(conic.h*x + conic.b*y + conic.f);

			double detJ = S1x*S2y - S2x*S1y;
			if (Abs(detJ)<CONICPREC) break;
			/*
			 * Inverse of the Jacobian
			 *           1   /  S2y   -S1y \
			 * invJ = ------ |             |
			 *         detJ  \ -S2x    S1x /
			 */

			double dx = -( S2y*S1 - S1y*S2) / detJ;
			double dy = -(-S2x*S1 + S1x*S2) / detJ;
			if (Abs(dx)<CONICPREC && Abs(dy)<CONICPREC) break;

			double xn = x + dx;
			double yn = y + dy;
			S1 = (*this)(xn,yn);
			S2 = conic(xn,yn);

			double SSprev  = SS;
			SS = Sqr(S1) + Sqr(S2);
			if (SS > SSprev) break;

			x = xn;
			y = yn;
		}
#if 0
		cout << "S1: " << *this << endl;
		cout << "S2: " << conic << endl;
		cout << "CONIC-CONIC (x,y)="<<x<<',' <<y << endl;
		cout << "S1(x,y)="<< (*this)(x, y) << "\tS2(x,y)="<< conic(x, y) << endl;
#endif

	}

	// Search if already exists
	for (int i=0; i<*n; i++) {
		if (Sqr(pts[i].x-x) + Sqr(pts[i].y-y) <=
		    CONICPREC*(1.0 + Sqr(pts[i].x) + Sqr(x)+ Sqr(pts[i].y) + Sqr(y)))
			return false;
	}
	//assert(*n<sizeof(Pts)/sizeof(Vector2D));
	//if (*n>=(int)(sizeof(pts)/sizeof(Vector2D))) {
	if (*n>=4) {
		// Fail safe check. Find the smallest distance between the various
		// pairs of points and if the min-distance between the new point
		// is smaller replace one
#if _DEBUG>0
		int ii, jj, in;
#endif
		double d2old = INFINITE;
		double d2new = INFINITE;
		for (int i=0; i<*n; i++) {
			for (int j=i+1; j<*n; j++) {
				double d2 = Sqr(pts[i].x-pts[j].x) + Sqr(pts[i].y-pts[j].y);
				if (d2<d2old) {
#if _DEBUG>0
					ii = i;
					jj = j;
#endif
					d2old = d2;
				}
			}
			double d2 = Sqr(pts[i].x-x) + Sqr(pts[i].y-y);
			if (d2<d2new) {
#if _DEBUG>0
				in = i;
#endif
				d2new = d2;
			}
		}
#if _DEBUG>0
		cerr << "ERROR: conic::addPoint n>sizeof(pts)" << endl;
		cerr << "   add= " << setprecision(22) << x << ",\t" << y << endl;
		for (int i=0; i<*n; i++)
			cerr << "pts[" << i << "]= " << setprecision(22) << pts[i].x << ",\t" << pts[i].y << endl;
		cerr << " d2old= " << setprecision(22) << d2old << "\t(i,j)= " << ii << "," << jj<< endl;
		cerr << " d2new= " << d2new << ",\t(in)=" << in << endl;
		if (d2new > d2old) {
			// New solution to be accepted
			pts[in].x = x;
			pts[in].y = y;
		}
#endif
		return false;
	}
	pts[*n].x = x;
	pts[*n].y = y;
	*n += 1;
	return true;
} // addPoint

/**
 * Intersect conic with a conic-line
 * @return number of intersection points
 */
inline int Conic::_intersectLine(const Conic& line, Vector2D pts[4]) const
{
	double x0, y0, vx, vy;

	// Find a point
	if (Abs(line.g) < Abs(line.f)) {
		x0 = 0.0;
		y0 = -line.c / (2.0*line.f);
	} else {
		x0 = -line.c / (2.0*line.g);
		y0 = 0.0;
	}

	// and the slope
	vx =  line.f;
	vy = -line.g;

	// line is in the form
	//    x = x0 + vx*t
	//    y = y0 + vy*t
	// substituting we end up with an equation like
	//    A*t^2 + 2B*t + C = 0
	double A = a*vx*vx + b*vy*vy + 2.0*h*vx*vy;
	double B = a*x0*vx + h*(x0*vy + y0*vx) + b*y0*vy + g*vx + f*vy;
	double C = a*x0*x0 + b*y0*y0 + c + 2.0*(h*x0*y0 + g*x0 + f*y0);
	int n;
	double t1,t2;

	if (Abs(A) < TOOSMALL) {
		// without this the following QUA fails
		// QUA MQI_PLS1   0.0 0.0 0.0 0.8009739843649878 0.0 0.0 0.0 0.0 0.0 -1.0
		if (Eq0(B,TOOSMALL)) return 0;
		t1 = -0.5*C / B;
		n = 1;
	} else
		n = quadratic(2.0*B/A,C/A,&t1,&t2,CONICPREC);

	if (n==1) {
		pts[0].x = x0+vx*t1;
		pts[0].y = y0+vy*t1;
		return 1;
	} else
	if (n==2) {
		pts[0].x = x0+vx*t1;
		pts[0].y = y0+vy*t1;
		pts[1].x = x0+vx*t2;
		pts[1].y = y0+vy*t2;
		return 2;
	} else
		return 0;
} // _intersectLine

/**
 * Intersect conic with another one and return intersection points
 * @return number of intersection points from 0 to 4
 *         fill array pts with coordinates
 */
int Conic::intersect(const Conic& conic, Vector2D pts[4]) const
{
	// If conics are equal return 0
	if (ISEQ(a,conic.a) && ISEQ(h,conic.h) && ISEQ(b,conic.b) &&
	    ISEQ(g,conic.g) && ISEQ(f,conic.f))
		return 0;
	// FIXME could be that conic = s * (this)
	//double s;

	if (_type == CONIC_LINE) {
		if (conic._type == CONIC_LINE) {
			double d = g * conic.f - f * conic.g;
			if (Eq0(d,CONICPREC))	// FIXME or check against TOOSMALL?
				return 0;	// Parallel lines
			pts[0].x = ( f*conic.c - c*conic.f ) / d / 2.0;
			pts[0].y = ( c*conic.g - g*conic.c ) / d / 2.0;
			return 1;
		} else
			return conic._intersectLine(*this, pts);
	} else
	if (conic._type == CONIC_LINE)
		return _intersectLine(conic, pts);
	else {
#define SUBSTITUTION
#ifdef SUBSTITUTION
		// S1(x,y) = 0
		// solve for y1,2 = f(x)
		// substitute y1,2 to S2(x,y)=0
		// move all sqrt() terms on right side
		// square equation and get quartic
		int n;
		double x[4];
		double AA, BB, CC, DD, EE;

		{
			double Ab = conic.a * b;
			double aB = a * conic.b;

			double Af = conic.a * f;
			double aF = a * conic.f;

			double Bc = conic.b * c;
			double bC = b * conic.c;

			double Bf = conic.b * f;
			double bF = b * conic.f;

			double Bg = conic.b * g;
			double bG = b * conic.g;

			double Bh = conic.b * h;
			double bH = b * conic.h;

			double Ch = conic.c * h;
			double cH = c * conic.h;

			double Fh = conic.f * h;
			double fH = f * conic.h;

			double Fg = conic.f * g;
			double fG = f * conic.g;

#define KAHAN
#ifndef KAHAN
			AA = Sqr(Ab-aB) + 4.0 * (Bh-bH) * (conic.a*h - a*conic.h);

			BB = 4.0 * ( Ab * (bG-Bg - fH-Fh)
					   + aB * (Bg-bG - fH-Fh)
					+ 2.0 * (Af*Bh + aF*bH + Bh*h*conic.g + bH*conic.h*g
						 - (bG + Bg ) * h*conic.h));

			CC = 2.0 * ( (Bc-bC) * (aB-Ab)
					 + 2.0 * (bF-Bf) * (aF-Af)
					 + 2.0 * (Ch-cH) * (Bh-bH)
					 + 2.0 * (Sqr(b)*Sqr(conic.g) + Sqr(conic.b)*Sqr(g))
					 - 4.0 * ((Fh+fH) * (bG+Bg) + bG*Bg)
					 + 8.0 * (fG*Bh + Fg*bH));

			DD = -4.0 * (Bc * (bG-Bg + fH+Fh)
					  + bC * (Bg-bG + fH+Fh)
					  + 2.0 * bF * (fG-Fg-cH)
					  + 2.0 * Bf * (Fg-fG-Ch));

			EE = Sqr(Bc-bC) + 4.0 * (conic.c*f - c*conic.f) * (Bf-bF);
#else
			KahanSum kAA(Sqr(Ab-aB));
			kAA( 4.0 * (Bh-bH)*conic.a*h);
			kAA(-4.0 * (Bh-bH)*a*conic.h);
			//if (AA != kAA()) cout << AA - kAA() << endl;
			AA = kAA();

			KahanSum kBB( Ab * (bG-Bg - fH-Fh));
			kBB(aB * (Bg-bG - fH-Fh));
			kBB( 2.0 * Af*Bh);
			kBB( 2.0 * aF*bH);
			kBB( 2.0 * Bh*h*conic.g);
			kBB( 2.0 * bH*conic.h*g);
			kBB(-2.0 * (bG+Bg)*h*conic.h);
			//if (Abs(BB - 4.0*kBB())>1e-10) cout << BB - 4.0*kBB() << endl;
			BB = 4.0 * kBB();

			KahanSum kCC((Bc-bC) * (aB-Ab));
			kCC( 2.0 * (bF-Bf) * (aF-Af));
			kCC( 2.0 * (Ch-cH) * (Bh-bH));
			kCC( 2.0 * Sqr(b)*Sqr(conic.g));
			kCC( 2.0 * Sqr(conic.b)*Sqr(g));
			kCC(-4.0 * (Fh+fH)*(bG+Bg));
			kCC(-4.0 * bG*Bg);
			kCC( 8.0 * (fG*Bh + Fg*bH));
			//if (Abs(CC - 2.0*kCC())>1e-7) cout << CC - 2.0*kCC() << endl;
			CC = 2.0*kCC();

			KahanSum kDD(Bc*(bG-Bg + fH+Fh));
			kDD(bC * (Bg-bG + fH+Fh));
			kDD(2.0 * bF * (fG-Fg-cH));
			kDD(2.0 * Bf * (Fg-fG-Ch));
			//cout << DD + 4.0*kDD() << endl;
			DD = -4.0*kDD();

			KahanSum kEE(Sqr(Bc-bC));
			kEE( 4.0*conic.c*f * (Bf-bF));
			kEE(-4.0*c*conic.f * (Bf-bF));
			//cout << EE - kEE() << endl;
			EE = kEE();
#endif
		}

		// Convert largest to be of the order of 1..100
		double aAA = Abs(AA);
		double aBB = Abs(BB);
		double aCC = Abs(CC);
		double aDD = Abs(DD);
		double aEE = Abs(EE);
		double max = Max(aAA,aBB,Max(aCC,aDD));

		if (Eq0(max,TOOSMALL)) return 0;

		if (max<1.0 || max>1000.0) {
			double scale = 10.0 / max;
			aAA *= scale;
			aBB *= scale;
			aCC *= scale;
			aDD *= scale;
			aEE *= scale;
		}

		if (aAA > TOOSMALL &&
		    aAA > CONICPREC*(aBB + CONICPREC*(aCC + CONICPREC*(aDD + CONICPREC*aEE)))) {
			BB /= AA;
			CC /= AA;
			DD /= AA;
			EE /= AA;
			//AA = 1.0;
			n = quartic(BB, CC, DD, EE, x, 100.0*CONICSMALL, 4);
			if (n==0) return 0;
		} else
		if (aBB > TOOSMALL &&
		     aBB > CONICPREC*(aCC + CONICPREC*(aDD + CONICPREC*aEE))) {
			CC /= BB;
			DD /= BB;
			EE /= BB;
			//BB = 1.0;
			//AA = 0.0;
			n = cubic(CC,DD,EE, x, 10.0*CONICSMALL, 4);
		} else
		if (aCC > TOOSMALL &&
		    aCC > TOOSMALL*(aDD + aEE)) {
			DD /= CC;
			EE /= CC;
			//CC = 1.0;
			//AA = 0.0;
			//BB = 0.0;
			n = quadratic(DD,EE, x, x+1, CONICSMALL);
			if (n==0) return 0;
		} else
		if (!Eq0(DD, TOOSMALL)) {
			n = 1;
			x[0] = -EE/DD;
		} else
			return 0;

		// Check solutions for y
		int j = 0;
		for (int i=0; i<n; i++) {
			double y1, y2, y3, y4;
			if (getY(x[i], &y1, &y2) == 0) continue;
			if (conic.getY(x[i], &y3, &y4) == 0) continue;

			// add the smallest distance
			// conic1:   1  2
			// conic2:   3  4
			double d13 = Abs(y1-y3);
			double d14 = Abs(y1-y4);
			double d23 = Abs(y2-y3);
			double d24 = Abs(y2-y4);

			// Find minimum solution (warning if = 0.0)
			double acc = Min(d13,d14,d23,d24)+CONICPREC;
			// Small is too big, ignore the solution
			if (acc > 1.0e-5*max) continue;
			acc *= 1.1;	// allow 10% margin

			if (d13<=acc) { // is d13 the smallest?
				addPoint(conic, &j, pts, x[i], 0.5*(y1+y3));
				// the remaining d24 can be also a candidate
				// if it is of the same order of d13
				if (d24 <= acc*10.0) addPoint(conic, &j, pts, x[i], 0.5*(y2+y4));
			} else
			if (d14<=acc) { // is d14 the smallest?
				addPoint(conic, &j, pts, x[i], 0.5*(y1+y4));
				// check remaining d23
				if (d23 <= acc*10.0) addPoint(conic, &j, pts, x[i], 0.5*(y2+y3));
			} else
			if (d23<=acc) {	// is d23 the smallest?
				addPoint(conic, &j, pts, x[i], 0.5*(y2+y3));
				// check remaining d14
				if (d14 <= acc*10.0) addPoint(conic, &j, pts, x[i], 0.5*(y1+y4));
			} else
			if (d24<=acc) {	 // d24 should be the smallest
				addPoint(conic, &j, pts, x[i], 0.5*(y2+y4));
				// check remaining d13
				if (d13 <= acc*10.0) addPoint(conic, &j, pts, x[i], 0.5*(y1+y3));
			}
		}

		return j;

#else
		// Pencil of conics

		Matrix3 E1, E2;
		matrix(&E1);
		conic.matrix(&E2);
#ifdef _DUMP
		cout << "E1=" << endl << E1 << endl;
		cout << "E2=" << endl << E2 << endl;
#endif
		E2.negate();
//		E2i.inverse(E2);
//		E2.multiply(E1,E2i);
		Matrix3 E2i(E2);
		E2i.inverse();
		E2.multiply(E1,E2i);
#ifdef _DUMP
		cout << "inv(-E2)=" << endl << E2i << endl;
		cout << "E1*inv(-E2)=" << endl << E2 << endl;
#endif

		double aa = -E2.trace();
		double bb =  E2(0,0)*E2(1,1) - E2(1,0)*E2(0,1) +
			     E2(1,1)*E2(2,2) - E2(2,1)*E2(1,2) +
			     E2(0,0)*E2(2,2) - E2(2,0)*E2(0,2);
		double cc = -E2.det();
		if (Abs(cc)<CONICPREC) return 0;
		DUMP(cout << "aa="<< aa << " bb="<< bb << " cc=" << cc << endl);

		// Find the two conic lines
		double y[3];
		Conic M, L;
		int n = cubic(aa,bb,cc,y,CONICPREC,4);
		for (int i=0; i<n; i++) {
			DUMP(cout << endl << "y"<<i<<"= "<<y[i]<<endl);
			double r = y[i];
			if (Eq0(r,CONICPREC)) continue;
			Conic C0(a+r*conic.a, h+r*conic.h, b+r*conic.b,
				 g+r*conic.g, f+r*conic.f, c+r*conic.c);
			DUMP(cout << "C0: " << C0 << endl);
			if (C0._type == CONIC_LINES) {
				C0.splitLines(&M, &L);
				break;
			} else
			if (C0._type == CONIC_LINE) {
				M.copy(C0);
				break;
			}
		}
		if (M._type != CONIC_LINE) return 0;
#ifdef _DUMP
//		M.parametric();
//		if (L._type != CONIC_NOTHING) L.parametric();
//		cout <<"M: "<< M << endl;
//		cout <<"L: "<< L << endl;
#endif

		n = _intersectLine(M, pts);
		int k;
		if (L._type == CONIC_LINE)
			k = _intersectLine(L, pts+n);
		else
			k = 0;

		return n+k;
#endif
	}
} // intersect

/** calculate parametric coefficients for conic */
void Conic::parametric()
{
	if (_type == CONIC_LINE) {
		// Convert a conic to a parametric equation of a line
		//  x = c1 + c2*t
		//  y = c4 + c5*t
		// Find starting point
		if (Abs(g) < Abs(f)) {
			c1 = 0.0;
			c4 = -c / (2.0*f);
		} else {
			c1 = -c / (2.0*g);
			c4 = 0.0;
		}

		c2 = -f;
		c5 =  g;

		assert(c2!=0.0 || c5!=0.0);
	} else {
		double xa, xb, xg, xf;
		double dx, dy;
		double ct, st;

		// convert conic to standard one
		if (ISZEROH(h)) {
			if (a<b) {
				// theta = pi
				ct = -1.0;
				st =  0.0;
				xa =  a;
				xb =  b;
				xg = -g;
				xf = -f;
			} else {
				// theta = 0
				ct = 1.0;
				st = 0.0;
				xa = a;
				xb = b;
				xg = g;
				xf = f;
			}
		} else {
			// remove h (rotation) component
			double th = theta();
			bsincos(th, &st, &ct);
			xa = a*ct*ct + 2.0*h*st*ct + b*st*st;
			xb = a*st*st - 2.0*h*st*ct + b*ct*ct;
			xg =  g*ct + f*st;
			xf = -g*st + f*ct;
		}

		if (_type == CONIC_PARABOLA) {
			// x = c1 + c2*t + c3*t**2
			// y = c4 + c5*t + c6*t**2
			// Assume that b=ZERO otherwise interchange axes
			if (Abs(xa) < Abs(xb)) {
				assert(xb != 0.0);
				assert(xg != 0.0);
				// with the following translation
				dx = 0.5*(xf*xf/xb - c)/xg;
				dy = - xf/xb;

				// should reduce to: b*y^2 + 2*g*x = 0
				// or y^2 = 2*sqrt(-D/I^3)*x
				// or y^2 = -2*g/b*x
				double k = -0.5*xb/xg;

				c2 = -st;
				c3 =  ct*k;

				c5 =  ct;
				c6 =  st*k;
			} else {
				assert(xa != 0.0);
				assert(xf != 0.0);
				// with the following translation
				dx = - xg/xa;
				dy = 0.5*(xg*xg/xa - c)/xf;

				// should reduce to: a*x^2 + 2*f*y = 0
				// or y^2 = 2*sqrt(-D/I^3)*x
				// or x^2 = -2*f/a * y
				double k = -0.5*xa/xf;

				c2 =  ct;
				c3 = -st*k;

				c5 =  st;
				c6 =  ct*k;
			}
			c1 =  dx*ct - dy*st;
			c4 =  dx*st + dy*ct;

		} else {
			// ELLIPSE & HYPERBOLA
			//
			// Convert a conic to a parametric equation of a hyperbola
			//	x = c1 + c2*sec(t) + c3*tan(t)		[sect = 1/cost]
			//	y = c4 + c5*sec(t) + c6*tan(t)
			//	x = c1 + c2/cos(t) + c3*tan(t)
			//	y = c4 + c5/cos(t) + c6*tan(t)
			//
			// or with hyperbolic functions
			//	x = c1 + c2*cosh(t) + c3*sinh(t)
			//	y = c4 + c5*cosh(t) + c6*sinh(t)
			// Note: hyperbolic will only represent the positive side
			//       so it is better the trig functions
			dx = -xg/xa;
			dy = -xf/xb;
			double xc = c - xg*xg/xa - xf*xf/xb;

			bool inv = xb*xc<0.0;

			xa = 1.0 / sqrt(Abs(xa/xc));
			xb = 1.0 / sqrt(Abs(xb/xc));

			c1 =  dx*ct - dy*st;
			c4 =  dx*st + dy*ct;

			// For hyperbola consider two cases
			// 1st:  ax2 - by2 -c = 0 <=> -ax2 + by2 +c = 0
			// 2nd: -ax2 + by2 -c = 0 <=>  ax2 - by2 +c = 0
			// therefore 2nd rotated: b*c < 0.0
			if (_type == CONIC_HYPERBOLA && inv) {
				c3 =  xa*ct;
				c2 = -xb*st;
				c6 =  xa*st;
				c5 =  xb*ct;
			} else {
				c2 =  xa*ct;
				c3 = -xb*st;
				c5 =  xa*st;
				c6 =  xb*ct;
			}
		}
	}

#if _DEBUG>1
#if 0
	// Test parametric equation
	double R=3.14;
	if (_type==CONIC_LINE || _type==CONIC_PARABOLA) R = 100.0;
	double sR = 2.0*R/100.0;

	for (double t=-R; t<R; t += sR) {
		double x,y;
		getXY(t,&x,&y);
		double t2=getT(x,y);
		double x2,y2;
		getXY(t2,&x2,&y2);
		if (Abs(t-t2)>TOLERANCE||
		   (Abs(x)<1.0 && Abs(x-x2)>TOLERANCE) ||
		   (Abs(x)>1.0 && Abs(x-x2)/Abs(x+x2)>2.0*TOLERANCE) ||
		   (Abs(y)<1.0 && Abs(y-y2)>TOLERANCE) ||
		   (Abs(y)>1.0 && Abs(y-y2)/Abs(y+y2)>2.0*TOLERANCE)) {
			cerr << *this << endl;
			cerr	<< "\tt="   << t
				<< " t2=" << t2
				<< " Dt=" << Abs(t-t2) << endl;
			cerr	<< "\tx="   << x
				<< " x2=" << x2
				<< " Dx=" << Abs(x-x2)
				<< " Dx/x=" << Abs(x-x2)/Abs(x+x2)/2.0 << endl;
			cerr	<< "\ty="   << y
				<< " y2=" << y2
				<< " Dy=" << Abs(y-y2)
				<< " Dy/y=" << Abs(y-y2)/Abs(y+y2)/2.0 << endl;
			throw 0;
		}
	}
#endif
#endif
} // parametric

/** getT parametric t
 * @param x,y	Cartesian coordinates on the conic
 * @return parametric t
 * */
double Conic::getT(double x, double y) const
{
	double A, B, C, DD;
	double sign1, sign2;
	double cost1, sint1, cost2, sint2;
	double sect1, tant1, sect2, tant2;
	double dx1, dx2, dx3, dx4;
	double dy1, dy2, dy3, dy4;
	double tmp;

	switch (_type) {
		case CONIC_LINE:
			if (Abs(c2) < Abs(c5))
				return (y-c4) / c5;
			else
				return (x-c1) / c2;
			break;

		case CONIC_ELLIPSE:
			x -= c1;
			y -= c4;

			// Assume z=cos(t), sin(t)=+/-sqrt(1-z^2)
			A  = c2*c2 + c3*c3;
			B  = x*c2;		// -2.0*...
			C  = x*x - c3*c3;
			//n = quadratic(2.0*B/A, C/A, &cost1, &cost2, CONICPREC);
			DD = B*B - A*C;

			if (DD > CONICSMALL*(B*B+Abs(A*C))) {
				if (B>0.0)
					DD = B + sqrt(DD);
				else
					DD = B - sqrt(DD);
				cost1 = DD / A;
				tmp = (1.0-cost1)*(1.0+cost1);
				if (tmp>0.0)
					sint1 = sqrt(tmp);
				else {
					cost1 = cost1>0.0? 1.0 : -1.0;
					sint1 = 0.0;
				}

				// Evaluate the parametric equation and compare against y
				dx1 = Abs(c2*cost1 + c3*sint1 - x);
				dy1 = Abs(c5*cost1 + c6*sint1 - y);
				if (dx1<CONICSMALL && dy1<CONICSMALL) return acos(cost1);

				dx2 = Abs(c2*cost1 - c3*sint1 - x);
				dy2 = Abs(c5*cost1 - c6*sint1 - y);
				if (dx2<CONICSMALL && dy2<CONICSMALL) return -acos(cost1);

				cost2 = C / DD;
				tmp = (1.0-cost2)*(1.0+cost2);
				if (tmp>0.0)
					sint2 = sqrt(tmp);
				else {
					cost2 = cost2>0.0? 1.0 : -1.0;
					sint2 = 0.0;
				}

				dx3 = Abs(c2*cost2 + c3*sint2 - x);
				dy3 = Abs(c5*cost2 + c6*sint2 - y);
				if (dx3<CONICSMALL && dy3<CONICSMALL) return acos(cost2);

				dx4 = Abs(c2*cost2 - c3*sint2 - x);
				dy4 = Abs(c5*cost2 - c6*sint2 - y);
				if (dx4<CONICSMALL && dy4<CONICSMALL) return -acos(cost2);

				// find the smallest distance
				dx1 = dx1*dx1 + dy1*dy1;
				dx2 = dx2*dx2 + dy2*dy2;
				dx3 = dx3*dx3 + dy3*dy3;
				dx4 = dx4*dx4 + dy4*dy4;
				// else return the smallest one
				if (dx1<=dx2 && dx1<=dx3 && dx1<=dx4) return  acos(cost1);
				if (dx2<=dx1 && dx2<=dx3 && dx2<=dx4) return -acos(cost1);
				if (dx3<=dx1 && dx3<=dx2 && dx3<=dx4) return  acos(cost2);
				return -acos(cost2);
			} else {
				cost1 = B / A;
				tmp = (1.0-cost1)*(1.0+cost1);
				if (tmp>0.0)
					sint1 = sqrt(tmp);
				else {
					cost1 = cost1>0.0? 1.0 : -1.0;
					sint1 = 0.0;
				}

				dx1 = Abs(c2*cost1 + c3*sint1 - x);
				dy1 = Abs(c5*cost1 + c6*sint1 - y);
				if (dx1<CONICSMALL && dy1<CONICSMALL) return acos(cost1);

				dx2 = Abs(c2*cost1 - c3*sint1 - x);
				dy2 = Abs(c5*cost1 - c6*sint1 - y);

				if (dx1*dx1 + dy1*dy1 < dx2*dx2 + dy2*dy2)
					return acos(cost1);
				else
					return -acos(cost1);
			}
			break;

		case CONIC_PARABOLA:
			if (Abs(c3)<CONICSMALL)
				return (x - c1) / c2;

			//           A        B       C
			// equation c3*t^2 + c2*t + c1-x = 0
			//A  = c3;
			//B  = c2;
			C  = c1 - x;
			DD = c2*c2 - 4.0*c3*C;

			if (DD > CONICSMALL*(c2*c2+Abs(c3*C))) {
				if (c2>0.0)
					DD = -c2 - sqrt(DD);
				else
					DD = -c2 + sqrt(DD);
				double t1 = 0.5 * DD / c3;

				//DD = sqrt(DD);
				//if (c2>0.0)
				//	t1 = 2.0*C / (-c2 - DD);
				//else
				//	t1 = (-c2 + DD) / (2.0*c3);

				// Evaluate the parametric equation and compare against x,y
				dx1 = Abs(c1 + c2*t1 + c3*t1*t1 - x);
				dy1 = Abs(c4 + c5*t1 + c6*t1*t1 - y);
				if (dx1<CONICSMALL && dy1<CONICSMALL) return t1;

				double t2 = 2.0 * C / DD;
				//if (c2<0.0)
				//	t2 = 2.0*C / (-c2 + DD);
				//else
				//	t2 = (-c2 - DD) / (2.0*c3);

				dx2 = Abs(c1 + c2*t2 + c3*t2*t2 - x);
				dy2 = Abs(c4 + c5*t2 + c6*t2*t2 - y);

				// find the smallest distance
				if (dx1*dx1 + dy1*dy1 < dx2*dx2 + dy2*dy2)
					return t1;
				else
					return t2;
			} else
				return -c2 / (2.0*c3);
			break;

		case CONIC_HYPERBOLA:
			x -= c1;
			y -= c4;

			// Assume z=sec(t), tan(t)=+/-sqrt(z^2-1)
			A  = (c2 - c3)*(c2 + c3);
			B  = x*c2;	// 2.0*...
			C  = x*x + c3*c3;
			DD = B*B - A*C;

			if (DD > CONICSMALL*(B*B+Abs(A*C))) {
				if (B>0.0)
					DD = B + sqrt(DD);
				else
					DD = B - sqrt(DD);
				sect1 = DD / A;

				if (sect1 > 0.0)
					sign1 =  1.0;
				else
					sign1 = -1.0;

				tmp = (sect1-1.0)*(sect1+1.0);
				if (tmp>0.0)
					tant1 = sqrt(tmp);
				else {
					sect1 = sign1;
					tant1 = 0.0;
				}

				// Evaluate the parametric equation and compare against y
				dx1 = Abs(c2*sect1 + c3*tant1 - x);
				dy1 = Abs(c5*sect1 + c6*tant1 - y);
				if (dx1<CONICSMALL && dy1<CONICSMALL) return sign1*acos(1.0/sect1);

				dx2 = Abs(c2*sect1 - c3*tant1 - x);
				dy2 = Abs(c5*sect1 - c6*tant1 - y);
				if (dx2<CONICSMALL && dy2<CONICSMALL) return -sign1*acos(1.0/sect1);

				sect2 = C / DD;
				if (sect2 > 0.0)
					sign2 =  1.0;
				else
					sign2 = -1.0;

				tmp = (sect2-1.0)*(sect2+1.0);
				if (tmp > 0.0)
					tant2 = sqrt(tmp);
				else {
					sect2 = sign2;
					tant2 = 0.0;
				}

				dx3 = Abs(c2*sect2 + c3*tant2 - x);
				dy3 = Abs(c5*sect2 + c6*tant2 - y);
				if (dx3<CONICSMALL && dy3<CONICSMALL) return sign2*acos(1.0/sect2);

				dx4 = Abs(c2*sect2 - c3*tant2 - x);
				dy4 = Abs(c5*sect2 - c6*tant2 - y);
				if (dx4<CONICSMALL && dy4<CONICSMALL) return -sign2*acos(1.0/sect2);

				// find the smallest distance
				dx1 = dx1*dx1 + dy1*dy1;
				dx2 = dx2*dx2 + dy2*dy2;
				dx3 = dx3*dx3 + dy3*dy3;
				dx4 = dx4*dx4 + dy4*dy4;

				// else return the smallest one
				if (dx1<=dx2 && dx1<=dx3 && dx1<=dx4) return  sign1*acos(1.0/sect1);
				if (dx2<=dx1 && dx2<=dx3 && dx2<=dx4) return -sign1*acos(1.0/sect1);
				if (dx3<=dx1 && dx3<=dx2 && dx3<=dx4) return  sign2*acos(1.0/sect2);
				return -sign2*acos(1.0/sect2);

			} else {
				sect1 = B / A;
				if (sect1 > 0.0)
					sign1 = 1.0;
				else
					sign1 = -1.0;

				tmp = (sect1-1.0)*(sect1+1.0);
				if (tmp>0.0)
					tant1 = sqrt(tmp);
				else {
					sect1 = sign1;
					tant1 = 0.0;
				}

				dx1 = Abs(c2*sect1 + c3*tant1 - x);
				dy1 = Abs(c5*sect1 + c6*tant1 - y);
				if (dx1<CONICSMALL && dy1<CONICSMALL)
					return sign1*acos(1.0/sect1);

				dx2 = Abs(c2*sect1 - c3*tant1 - x);
				dy2 = Abs(c5*sect1 - c6*tant1 - y);

				if (dx1*dx1 + dy1*dy1 < dx2*dx2 + dy2*dy2)
					return  sign1*acos(1.0/sect1);
				else
					return -sign1*acos(1.0/sect1);
			}
			break;

		default:
			assert(0);
	}
	return 0.0;
} // _getT

/** Get X,Y coordinates from parametric t
 * @param t	parametric t
 * @param x,y	return coordinates
 * */
void Conic::getXY(const double t, double *x, double *y) const
{
	double ct, st;

	switch (_type) {
		case CONIC_LINE:
			*x = c1 + c2*t;
			*y = c4 + c5*t;
			break;

		case CONIC_ELLIPSE:
			bsincos(t, &st, &ct);
			*x = c1 + c2*ct + c3*st;
			*y = c4 + c5*ct + c6*st;
			break;

		case CONIC_PARABOLA:
			*x = c1 + (c2 + c3*t)*t;
			*y = c4 + (c5 + c6*t)*t;
			break;

		case CONIC_HYPERBOLA:
			ct = 1.0/cos(t);
			st = tan(t);
			*x = c1 + c2*ct + c3*st;
			*y = c4 + c5*ct + c6*st;
			break;

		default:
			*x = 0.0;
			*y = 0.0;
			assert(0);
	}
} // getXY

/** Get X,Y coordinates from parametric t */
int Conic::getY(const double x, double *y1, double *y2) const
{
	if (Eq0(b, TOOSMALL)) {
		double deno = 2.0*(h*x + f);
		if (Abs(deno) > CONICPREC) {
			*y1 = *y2 = -((a*x+2.0*g)*x+c)/deno;
			return 1;
		}
		return 0;
	}

	double C = (a*x + 2.0*g)*x + c;
	return quadratic(2.0*(f+h*x)/b, C/b, y1, y2,CONICSMALL);
} // getY

/**
 * return gradient vector at position (x,y).
 * It assumes that (x,y) is ON the surface.
 * @param x,y	point to find gradient
 * @param dx,dy	gradient vector
 */
void Conic::grad(const double x, const double y, double *dx, double *dy) const
{
	*dx = a*x + h*y + g;
	*dy = h*x + b*y + f;
	double invlen  = sqrt(Sqr(*dx) + Sqr(*dy));
	if (Eq0(invlen, TOOSMALL)) {
		*dx = *dy = 0.0;
		return;
	} else
	invlen  = 1.0/invlen;
	*dx *= invlen;
	*dy *= invlen;
} // grad

/** approximation to distance of point to conic.
 * My algorithm works well when the point is close to conic.
 * using the normal on the surface constructs a line and finds
 * the intersection of the line with the conic
 * @param x,y	point to find distance
 * @return distance (negative inside, positive outside)
 */
double Conic::adist(const double x, const double y) const
{
	double vx, vy;
	grad(x,y,&vx,&vy);

	double A = a*vx*vx + 2.0*h*vx*vy + b*vy*vy;
	double B = -2.0*(a*vx*x + h*(vx*y+x*vy) + b*vy*y + g*vx + f*vy);
	double C = (*this)(x,y);

	if (Eq0(A,CONICPREC)) {
		if (Eq0(B,CONICPREC))
			return INFINITE;
		return -C/B;
	} else {
		double t1,t2;
		int n = quadratic(B/A, C/A, &t1, &t2, CONICPREC);
		if (n) {
			if (Abs(t1)<Abs(t2))
				return t1;
			else
				return t2;
		} else
			return INFINITE;
	}
} // adist

/** equal
 * @param conic	to check equality with
 * @param acc	accuracy of operation
 * @return true if conics are equal
 */
bool Conic::equal(const Conic& conic, const double acc) const
{
	// Initially use the sum as scaling
	double m1 = a + h + b + g + f + c;
	double m2 = conic.a + conic.h + conic.b + conic.g + conic.f + conic.c;

	// if any of the two if zero
	if (Eq0(m1,acc) || Eq0(m2,acc)) {
		// Use the absolute maximum
		m1 = Max(Max(Abs(a), Abs(h), Abs(b)),
				Max(Abs(g), Abs(f), Abs(c)));
		m2 = Max(Max(Abs(conic.a), Abs(conic.h), Abs(conic.b)),
				Max(Abs(conic.g), Abs(conic.f), Abs(conic.c)));

		// Check again for zero
		if (Eq0(m1,acc) || Eq0(m2,acc)) {
			// compare only the maximum
			// return true if both are zero!
			return Eq0(m1,acc) && Eq0(m2,acc);
		}
	}
	m1 = 1.0 / m1;
	m2 = 1.0 / m2;

	// if we are using the maximum but this = -conic
	// it will fail
	if (!Eq0(a*m1 - conic.a*m2, acc)) return false;
	if (!Eq0(h*m1 - conic.h*m2, acc)) return false;
	if (!Eq0(b*m1 - conic.b*m2, acc)) return false;
	if (!Eq0(g*m1 - conic.g*m2, acc)) return false;
	if (!Eq0(f*m1 - conic.f*m2, acc)) return false;
	if (!Eq0(c*m1 - conic.c*m2, acc)) return false;

	return true;
} // equal

/**
 * operator <<
 * Return a string representation of the conic
 */
ostream& operator << (ostream& s, const Conic& conic)
{
	s << "S=" << setprecision(22);
	fmt(s, conic.a, "x^2");
	fmt(s, 2.0*conic.h, "x*y");
	fmt(s, conic.b, "y^2");
	fmt(s, 2.0*conic.g, "x");
	fmt(s, 2.0*conic.f, "y");
	fmt(s, conic.c);
	s << "=0" << endl;

	s << "\tD="<<conic.D << endl;
	s << "\tI="<<conic.I<<", J="<<conic.J << ", K=" << conic.K << endl;
	s << "\ttype="<< conic.typeStr();

	if (!ISZERO(conic.c1) || !ISZERO(conic.c2) || !ISZERO(conic.c3) ||
	    !ISZERO(conic.c4) || !ISZERO(conic.c5) || !ISZERO(conic.c6)) {
		s << endl << "\tx(t)=";
		switch (conic._type) {
			case CONIC_LINE:
				fmt(s, conic.c1);
				fmt(s, conic.c2, "t");
				s << endl << "\ty(t)=";
				fmt(s, conic.c4);
				fmt(s, conic.c5, "t");
				break;

			case CONIC_ELLIPSE:
				fmt(s, conic.c1);
				fmt(s, conic.c2, "cos(t)");
				fmt(s, conic.c3, "sin(t)");
				s << endl << "\ty(t)=";
				fmt(s, conic.c4);
				fmt(s, conic.c5, "cos(t)");
				fmt(s, conic.c6, "sin(t)");
				break;

			case CONIC_PARABOLA:
				fmt(s, conic.c1);
				fmt(s, conic.c2, "t");
				fmt(s, conic.c3, "t**2");
				s << endl << "\ty(t)=";
				fmt(s, conic.c4);
				fmt(s, conic.c5, "t");
				fmt(s, conic.c6, "t**2");
				break;

			case CONIC_HYPERBOLA:
				fmt(s, conic.c1);
				fmt(s, conic.c2, "1.0/cos(t)");
				fmt(s, conic.c3, "tan(t)");
				s << endl << "\ty(t)=";
				fmt(s, conic.c4);
				fmt(s, conic.c5, "1.0/cos(t)");
				fmt(s, conic.c6, "tan(t)");
				break;

			default:
				;
				// Do nothing;
		}
	}
#if 0
#if _DEBUG>1
#   ifdef _DUMP
	printf("\n");
	printf("a=%16lX ", *(long*)&(conic.a));
	printf("h=%16lX ", *(long*)&(conic.h));
	printf("b=%16lX\n",*(long*)&(conic.b));
	printf("g=%16lX ", *(long*)&(conic.g));
	printf("f=%16lX ", *(long*)&(conic.f));
	printf("c=%16lX\n",*(long*)&(conic.c));
#   endif
#endif
#endif
	return s;
} /* operator << */
