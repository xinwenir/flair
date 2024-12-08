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
#include <assert.h>
#include <string.h>
#include <iostream>

#include "os.h"
#include "eps.h"
#include "bbox.h"
#include "quad.h"
#include "bmath.h"
#include "obbox.h"
#include "format.h"
#include "bstring.h"
#include "matrix3.h"
#include "matrix4.h"

using namespace std;

/** set quad parameters */
void Quad::set( const double cxx, const double cyy, const double czz,
		const double cxy, const double cxz, const double cyz,
		const double cx,  const double cy,  const double cz,
		const double c)
{
	plane = false;
	Cxx = cxx;
	Cyy = cyy;
	Czz = czz;
	Cxy = cxy;
	Cxz = cxz;
	Cyz = cyz;
	Cx  = cx;
	Cy  = cy;
	Cz  = cz;
	C   = c;
	computeAbs();
	chk4Plane();
} // set

/** set quad parameters */
void Quad::set( const double cx,  const double cy,  const double cz,
		const double c)
{
	plane = true;
	Cxx = 0.0;
	Cyy = 0.0;
	Czz = 0.0;
	Cxy = 0.0;
	Cxz = 0.0;
	Cyz = 0.0;
	Cx  = cx;
	Cy  = cy;
	Cz  = cz;
	C   = c;
	computeAbs();
} // set

/** set quad parameters */
void Quad::set(const Quad& q)
{
	plane = q.plane;
	Cxx = q.Cxx;
	Cyy = q.Cyy;
	Czz = q.Czz;
	Cxy = q.Cxy;
	Cxz = q.Cxz;
	Cyz = q.Cyz;
	Cx  = q.Cx;
	Cy  = q.Cy;
	Cz  = q.Cz;
	C   = q.C;
	computeAbs();
} // set

/** set quad from matrix */
void Quad::set(const Matrix4& m)
{
	Cxx = m(0,0);
	Cyy = m(1,1);
	Czz = m(2,2);

	Cxy = 2.0 * m(1,0);
	Cxz = 2.0 * m(2,0);
	Cyz = 2.0 * m(2,1);

	Cx  = 2.0 * m(3,0);
	Cy  = 2.0 * m(3,1);
	Cz  = 2.0 * m(3,2);

	C = m(3,3);

	computeAbs();
	chk4Plane();
} // set

/** normalize coefficients and calculate accuracy parameter
 * for checking operations
 */
void Quad::normalize()
{
	double norm;
	if (plane)
		norm = Max(Sqr(Cx),  Sqr(Cy),  Sqr(Cz));
	else
		norm =	Max(Max(Abs(Cxx), Abs(Cyy), Abs(Czz)),
			    Max(Abs(Cxy), Abs(Cxz), Abs(Cyz)));
	if (Abs(norm) > SMALL) {
		norm = 1.0/norm;
		Cxx *= norm;
		Cyy *= norm;
		Czz *= norm;
		Cxy *= norm;
		Cxz *= norm;
		Cyz *= norm;
		Cx  *= norm;
		Cy  *= norm;
		Cz  *= norm;
		C   *= norm;
	}
/*
	if (Eq0(Cxx,TOOSMALL)) Cxx = 0.0;
	if (Eq0(Cyy,TOOSMALL)) Cyy = 0.0;
	if (Eq0(Czz,TOOSMALL)) Czz = 0.0;

	if (Eq0(Cxy,TOOSMALL)) Cxy = 0.0;
	if (Eq0(Cxz,TOOSMALL)) Cxz = 0.0;
	if (Eq0(Cyz,TOOSMALL)) Cyz = 0.0;

	if (Eq0(Cx,TOOSMALL)) Cx = 0.0;
	if (Eq0(Cy,TOOSMALL)) Cy = 0.0;
	if (Eq0(Cz,TOOSMALL)) Cz = 0.0;

	if (Eq0(C,TOOSMALL)) C = 0.0;
*/
} // normalize

/** bbox
 * Return a 3-dimensional axis-aligned bounding box

 * Find x-limits, grad(Q) = [x,0,0]
 * 1. Q(x,y,z)  = 0
 * 2. diff(Q,y) = 0
 * 3. diff(Q,z) = 0
 * 4. diff(Q,x) <> 0	If 0 then the gradient is zero and is not define
 *                      e.g. the tip of a cone
 *
 * solve([diff(Q(x,y,z),y)=0,diff(Q(x,y,z),z)=0],[y,z])
 *
 * y = 1/(4*Czz*Cyy - Cyz^2) * (-2*Cxy*Czz*x - 2*Czz*Cy +   x*Cxz*Cyz +   Cyz*Cz)
 * z = 1/(4*Cyy*Czz - Cyz^2) * (   Cxy*Cyz*x +   Cyz*Cy - 2*x*Cxz*Cyy - 2*Cyy*Cz)
 *
 * setting D = 4*Cyy*Czz - Cyz**2
 * y = [(Cxz*Cyz - 2*Cxy*Czz)*x + (Cyz*Cz- 2*Czz*Cy)] / D
 *   = K * x + L
 *     K = (Cxz*Cyz - 2*Cxy*Czz) / D
 *     L = (Cyz*Cz- 2*Czz*Cy) / D
 *
 * z = [(Cxy*Cyz - 2*Cxz*Cyy)*x + (Cyz*Cy- 2*Cyy*Cz)] / D
 *   = M * x + N
 *     M = (Cxy*Cyz - 2*Cxz*Cyy) / D
 *     N = (Cyz*Cy- 2*Cyy*Cz) / D
 *
 * substitute in #1 we have

 + Cxx        *x^2
 + Cxy*K      *x^2
 + Cxz*M      *x^2
 + Cyy*K^2    *x^2
 + Czz*M^2    *x^2
 + Cyz*K*M    *x^2

 + Cx         *x
 + Cxy*L      *x
 + Cxz*N      *x
 + 2*Cyy*L*K  *x
 + 2*Czz*N*M  *x
 + Cyz*L*M    *x
 + Cyz*K*N    *x
 + Cy*K       *x
 + Cz*M       *x

 + Cyy*L^2
 + Czz*N^2
 + Cyz*L*N
 + Cy*L
 + Cz*N
 + C
 */
static void _bboxLimits(double cxx, double cyy, double czz,
			double cxy, double cxz, double cyz,
			double cx,  double cy,  double cz,
			double c,   double norm,
			double *xmin, double *xmax)
{
	*xmin =  INFINITE;
	*xmax = -INFINITE;

	double D = 4.0*cyy*czz - Sqr(cyz);
	/* check if D is zero */
	if (!Eq0(D,SMALL)) {
		double K = (cxz*cyz - 2.0*cxy*czz) / D;
		double L = (cyz*cz  - 2.0*czz*cy ) / D;
		double M = (cxy*cyz - 2.0*cxz*cyy) / D;
		double N = (cyz*cy  - 2.0*cyy*cz ) / D;

		double AA = cxx
			+ cxy*K
			+ cxz*M
			+ cyy*Sqr(K)
			+ czz*Sqr(M)
			+ cyz*K*M;

		double BB = cx
			+ cxy*L
			+ cxz*N
			+ 2.0*cyy*L*K
			+ 2.0*czz*N*M
			+ cyz*L*M
			+ cyz*K*N
			+ cy*K
			+ cz*M;

		double CC = cyy*Sqr(L)
			+ czz*Sqr(N)
			+ cyz*L*N
			+ cy*L
			+ cz*N
			+ c;

		if (!Eq0(AA,SMALL)) {
			BB /= AA;
			CC /= AA;
			double x1, x2;
			if (quadratic(BB,CC,&x1,&x2,SMALL)) {
				// Check if dQ/dx
				double y = K*x1 + L;
				double z = M*x1 + N;
				double dx = 2*cxx*x1 + cxy*y + cxz*z + cx;
				//cout << "K=" << K << "L=" << L << endl;
				//cout << "M=" << M << "N=" << N << endl;
				//cout << "norm=" << norm << endl;
				//cout << "x1=" << x1 << endl;
				//cout << "dQ/dx(" << x1 << ") = " << dx << endl;
				//cout << "err =" << (Abs(2*cxx*x1) + Abs(cxy*y) + Abs(cxz*z) + Abs(cx))<< endl;
				if (!Eq0(dx, norm*1e-5))
					*xmin = x1;

				y = K*x2 + L;
				z = M*x2 + N;
				dx = 2*cxx*x2 + cxy*y + cxz*z + cx;
				//cout << "x2=" << x2 << endl;
				//cout << "dQ/dx(" << x2 << ") = " << dx << endl;
				//cout << "err =" << (Abs(2*cxx*x2) + Abs(cxy*y) + Abs(cxz*z) + Abs(cx))<< endl;
				if (!Eq0(dx, norm*1e-5))
					*xmax = x2;
			}
		}
	}
} // _bboxLimits

/** @return bounding box of quadric (in case of ellipsoid/cylinder/...) */
BBox Quad::bbox() const
{
	BBox bb;
	double xmin, xmax;
	double ymin, ymax;
	double zmin, zmax;

	double norm =	Max(Max(aCxx, aCyy, aCzz),
			    Max(Abs(Cxy), Abs(Cxz), Abs(Cyz)),
			    Max(Sqr(Cx),  Sqr(Cy),  Sqr(Cz)),
			    Abs(C));

	//cout << "Q=" << *this <<endl;

	/* limits on X */
	_bboxLimits(	Cxx, Cyy, Czz,
			Cxy, Cxz, Cyz,
			Cx,  Cy,  Cz,
			C, norm,
			&xmin, &xmax);

	/* limits on Y */
	_bboxLimits(	Cyy, Czz, Cxx,
			Cyz, Cxy, Cxz,
			Cy,  Cz,  Cx,
			C, norm,
			&ymin, &ymax);

	/* limits on Z */
	_bboxLimits(	Czz, Cxx, Cyy,
			Cxz, Cyz, Cxy,
			Cz,  Cx,  Cy,
			C, norm,
			&zmin, &zmax);

	bb.add(xmin, ymin, zmin);
	bb.add(xmax, ymax, zmax);

	if (!bb.isValid())
		bb.infinite();

	return bb;
} // bbox

/** @return oriented bounding box of quadric (in case of ellipsoid/cylinder/...) */
OBBox Quad::obbox() const
{
	OBBox bb;
	double xmin, xmax;
	double ymin, ymax;
	double zmin, zmax;

	double norm =	Max(Max(aCxx, aCyy, aCzz),
			    Max(Abs(Cxy), Abs(Cxz), Abs(Cyz)),
			    Max(Sqr(Cx),  Sqr(Cy),  Sqr(Cz)),
			    Abs(C));

	/* limits on X */
	_bboxLimits(	Cxx, Cyy, Czz,
			Cxy, Cxz, Cyz,
			Cx,  Cy,  Cz,
			C, norm,
			&xmin, &xmax);

	/* limits on Y */
	_bboxLimits(	Cyy, Czz, Cxx,
			Cyz, Cxy, Cxz,
			Cy,  Cz,  Cx,
			C, norm,
			&ymin, &ymax);

	/* limits on Z */
	_bboxLimits(	Czz, Cxx, Cyy,
			Cxz, Cyz, Cxy,
			Cz,  Cx,  Cy,
			C, norm,
			&zmin, &zmax);

	bb.add(xmin, ymin, zmin);
	bb.add(xmax, ymax, zmax);

	if (!bb.isValid())
		bb.infinite();

	return bb;
} // obbox

/** compute absolute values needed for accuracy */
void Quad::computeAbs()
{
	aCxx = Abs(Cxx);
	aCyy = Abs(Cyy);
	aCzz = Abs(Czz);
	aCxy = Abs(Cxy);
	aCxz = Abs(Cxz);
	aCyz = Abs(Cyz);
	aCxyz0 = Abs(Cx) + Abs(Cy) + Abs(Cz) + Abs(C);
} // computeAbs

/** check if quad represents a plane */
void Quad::chk4Plane(const double eps)
{
	double norm =	Max(Max(aCxx, aCyy, aCzz),
			    Max(Abs(Cxy), Abs(Cxz), Abs(Cyz)),
			    Max(Sqr(Cx),  Sqr(Cy),  Sqr(Cz)));
	norm = sqrt(norm);
	norm = (Abs(norm)>SMALL? 1.0/norm : 1.0);

	plane = Eq0(Cxx*norm,eps) && Eq0(Cyy*norm,eps) && Eq0(Czz*norm,eps) &&
	        Eq0(Cxy*norm,eps) && Eq0(Cxz*norm,eps) && Eq0(Cyz*norm,eps);
} // chk4Plane

/** First order approximation of the distance |f(x,y,z)| / ||grad(x,y,z)||
 * @param x,y,z	point to find distance
 * @return distance (negative inside, positive outside)
 */
double Quad::dist1(const double x, const double y, const double z) const
{
	double dx, dy, dz;
	grad(x,y,z,&dx,&dy,&dz);
	return f(x,y,z) / sqrt(dx*dx + dy*dy + dz*dz);
} // dist1

/** Second order approximation of the distance
 * @param x,y,z	point to find distance
 * @return distance (negative inside, positive outside)
 */
double Quad::dist2(const double x, const double y, const double z) const
{
	double x1, x2;
	double a = -sqrt(0.5*(Sqr(Cxy) + Sqr(Cxz) + Sqr(Cyz)) + Sqr(Cxx) + Sqr(Cyy) + Sqr(Czz));
	double b = -sqrt(Sqr(Cx + 2.0*Cxx*x +     Cxy*y +     Cxz*z) +
			 Sqr(Cy +     Cxy*x + 2.0*Cyy*y +     Cyz*z) +
			 Sqr(Cz +     Cyy*x +     Cyz*y + 2.0*Czz*z) );
	double c = Abs(f(x,y,z));
	quadratic(b/a, c/a, &x1, &x2, SMALL);
	if (x1<0.0 && x2>=0.0) return x2;
	if (x2<0.0 && x1>=0.0) return x1;
	return Min(x1,x2);
} // dist2

/** My approximation to distance of point to quadric.
 * The algorithm works well when the point is close to quad.
 * using the normal on the surface constructs a line and finds
 * the intersection of the line with the conic
 * @param x,y,z	point to find distance
 * @return distance (negative inside, positive outside)
 */
double Quad::adist(const double x, const double y, const double z) const
{
	double dx, dy, dz;
	grad(x,y,z,&dx,&dy,&dz);
	double dlen = 1.0 / sqrt(dx*dx + dy*dy + dz*dz);
	dx *= dlen;
	dy *= dlen;
	dz *= dlen;

	double t1,t2;
	int n = intersect(x,y,z, dx,dy,dz, &t1, &t2);
	if (n) {
		if (Abs(t1) < Abs(t2))
			return t1;
		else
			return t2;
	} else
		return INFINITE;
} // adist

/**
 * intersect a line with the quad and return the two intersection distances
 * @param p	position vector
 * @param v	direction vector (has to be normalized)
 * @param tmin, tmax	distance of minimum, maximum intersection points
 * @return	number of intersection points,
 *		if positive then [tmin..tmax]
 *		if negative then (-inf..tmin] | [tmax..+inf)
 */
int Quad::intersect(const double  x, const double  y, const double  z,
		    const double dx, const double dy, const double dz,
		    double *tmin, double *tmax) const
{
	if (plane) {
		double b = Cx*dx + Cy*dy + Cz*dz;
		double c = f(x, y, z);

		if (Eq0(b, SMALL)) {
			*tmin = *tmax = INFINITE;
			return 0;
		} else {
			if (b > 0.0) {
				*tmin = -INFINITE;
				*tmax = -c/b;
				return 1;
			} else {
				*tmin = -c/b;
				*tmax = INFINITE;
				return 1;
			}
		}
	} else {
		double a = (Cxx*dx + Cxy*dy + Cxz*dz)*dx +
			   (Cyy*dy + Cyz*dz)*dy +
			   Czz*dz*dz;
		double b = 2.0*(Cxx*dx*x + Cyy*dy*y + Czz*dz*z)
			 + Cxy*(dx*y + x*dy)
			 + Cxz*(dx*z + x*dz)
			 + Cyz*(dy*z + y*dz)
			 + Cx*dx + Cy*dy + Cz*dz;
		double c = f(x, y, z);

		if (Eq0(a, SMALL)) {
			if (Eq0(b, SMALL)) {
				*tmin = *tmax = INFINITE;
				return 0;
			} else {
				if (b > 0.0) {
					*tmin = -INFINITE;
					*tmax = -c/b;
					return 1;
				} else {
					*tmin = -c/b;
					*tmax = INFINITE;
					return 1;
				}
			}
		} else {
			b /= a;
			c /= a;
			int n = quadratic(b, c, tmin, tmax, SMALL);
			if (*tmin>*tmax) Swap(*tmin, *tmax);

			// Check for a negative object
			//  (-inf..tmin] | [tmax..+inf)
			// if middle point is inside or outside
			if (n>1) {
				double t = 0.5*(*tmin+*tmax);
				//if (f(x+t*dx, y+t*dy, z+t*dz) > 0.0)
				double ymid = t*t+b*t+c;	// times the sign-of-a
				if ((a>0.0 && ymid>0.0) || (a<0.0 && ymid<0.0))
					n = -n;
			}
			return n;
		}
	}
} // intersect

/** return quad matrix */
void Quad::matrix3(Matrix3 *m) const
{
	(*m)(0,0) = Cxx;
	(*m)(1,1) = Cyy;
	(*m)(2,2) = Czz;

	(*m)(0,1) = Cxy/2.0; (*m)(1,0) = Cxy/2.0;
	(*m)(0,2) = Cxz/2.0; (*m)(2,0) = Cxz/2.0;
	(*m)(1,2) = Cyz/2.0; (*m)(2,1) = Cyz/2.0;
} // matrix

/** return quad matrix */
void Quad::matrix(Matrix4 *m) const
{
	(*m)(0,0) = Cxx;
	(*m)(1,1) = Cyy;
	(*m)(2,2) = Czz;

	(*m)(0,1) = Cxy/2.0; (*m)(1,0) = Cxy/2.0;
	(*m)(0,2) = Cxz/2.0; (*m)(2,0) = Cxz/2.0;
	(*m)(1,2) = Cyz/2.0; (*m)(2,1) = Cyz/2.0;

	(*m)(0,3) = Cx/2.0;  (*m)(3,0) = Cx/2.0;
	(*m)(1,3) = Cy/2.0;  (*m)(3,1) = Cy/2.0;
	(*m)(2,3) = Cz/2.0;  (*m)(3,2) = Cz/2.0;

	(*m)(3,3) = C;
} // matrix

/**
 * Translate axes of conic by dx, dy
 * WARNING: to translate the conic use -dx, -dy
 */
void Quad::translate(const double dx, const double dy, const double dz)
{
	C += (Cxx*dx + Cxy*dy + Cxz*dz + Cx)*dx +
	     (Cyy*dy + Cyz*dz + Cy)*dy +
	     (Czz*dz + Cz)*dz;
	Cx = 2.0*Cxx*dx + Cxy*dy + Cxz*dz + Cx;
	Cy = 2.0*Cyy*dy + Cxy*dx + Cyz*dz + Cy;
	Cz = 2.0*Czz*dz + Cxz*dx + Cyz*dy + Cz;
} // translate

/** transform the quadric using matrix M     (M*Q*Mt)
 * @param M	matrix to transform quadric
 */
void Quad::transform(const Matrix4& M)
{
	if (plane) {
		double a = Cx;
		double b = Cy;
		double c = Cz;

		Cx = a*M(0,0) + b*M(1,0) + c*M(2,0);
		Cy = a*M(0,1) + b*M(1,1) + c*M(2,1);
		Cz = a*M(0,2) + b*M(1,2) + c*M(2,2);
		C += a*M(0,3) + b*M(1,3) + c*M(2,3);
	} else {
		double m00 = Cxx*M(0,0) + 0.5*(Cxy*M(1,0) + Cxz*M(2,0) + Cx*M(3,0));
		double m01 = Cxx*M(0,1) + 0.5*(Cxy*M(1,1) + Cxz*M(2,1) + Cx*M(3,1));
		double m02 = Cxx*M(0,2) + 0.5*(Cxy*M(1,2) + Cxz*M(2,2) + Cx*M(3,2));
		double m03 = Cxx*M(0,3) + 0.5*(Cxy*M(1,3) + Cxz*M(2,3) + Cx*M(3,3));

		double m10 = Cyy*M(1,0) + 0.5*(Cxy*M(0,0) + Cyz*M(2,0) + Cy*M(3,0));
		double m11 = Cyy*M(1,1) + 0.5*(Cxy*M(0,1) + Cyz*M(2,1) + Cy*M(3,1));
		double m12 = Cyy*M(1,2) + 0.5*(Cxy*M(0,2) + Cyz*M(2,2) + Cy*M(3,2));
		double m13 = Cyy*M(1,3) + 0.5*(Cxy*M(0,3) + Cyz*M(2,3) + Cy*M(3,3));

		double m20 = Czz*M(2,0) + 0.5*(Cxz*M(0,0) + Cyz*M(1,0) + Cz*M(3,0));
		double m21 = Czz*M(2,1) + 0.5*(Cxz*M(0,1) + Cyz*M(1,1) + Cz*M(3,1));
		double m22 = Czz*M(2,2) + 0.5*(Cxz*M(0,2) + Cyz*M(1,2) + Cz*M(3,2));
		double m23 = Czz*M(2,3) + 0.5*(Cxz*M(0,3) + Cyz*M(1,3) + Cz*M(3,3));

		double m30 = 0.5*(Cx*M(0,0) + Cy*M(1,0) + Cz*M(2,0)) + C*M(3,0);
		double m31 = 0.5*(Cx*M(0,1) + Cy*M(1,1) + Cz*M(2,1)) + C*M(3,1);
		double m32 = 0.5*(Cx*M(0,2) + Cy*M(1,2) + Cz*M(2,2)) + C*M(3,2);
		double m33 = 0.5*(Cx*M(0,3) + Cy*M(1,3) + Cz*M(2,3)) + C*M(3,3);

		Cxx =      M(0,0)*m00 + M(1,0)*m10 + M(2,0)*m20 + M(3,0)*m30;
		Cxy = 2.0*(M(0,0)*m01 + M(1,0)*m11 + M(2,0)*m21 + M(3,0)*m31);
		Cxz = 2.0*(M(0,0)*m02 + M(1,0)*m12 + M(2,0)*m22 + M(3,0)*m32);
		Cx  = 2.0*(M(0,0)*m03 + M(1,0)*m13 + M(2,0)*m23 + M(3,0)*m33);

		Cyy =      M(0,1)*m01 + M(1,1)*m11 + M(2,1)*m21 + M(3,1)*m31;
		Cyz = 2.0*(M(0,1)*m02 + M(1,1)*m12 + M(2,1)*m22 + M(3,1)*m32);
		Cy  = 2.0*(M(0,1)*m03 + M(1,1)*m13 + M(2,1)*m23 + M(3,1)*m33);

		Czz =      M(0,2)*m02 + M(1,2)*m12 + M(2,2)*m22 + M(3,2)*m32;
		Cz  = 2.0*(M(0,2)*m03 + M(1,2)*m13 + M(2,2)*m23 + M(3,2)*m33);

		double tmp = M(0,3)*m03 + M(1,3)*m13 + M(2,3)*m23 + M(3,3)*m33;
		C = (Abs(tmp) < Abs(C)*SMALL)? 0.0 : tmp;
	}
} // transform

/** negate quad coefficients */
void Quad::negate()
{
	Cxx = -Cxx;
	Cyy = -Cyy;
	Czz = -Czz;
	Cxy = -Cxy;
	Cxz = -Cxz;
	Cyz = -Cyz;
	Cx  = -Cx;
	Cy  = -Cy;
	Cz  = -Cz;
	C   = -C;
} // negate

/** operator << */
std::ostream& operator << (std::ostream& s, const Quad& q)
{
	if (!q.plane) {
		fmt(s, q.Cxx, "x^2"); s << ' ';
		fmt(s, q.Cyy, "y^2"); s << ' ';
		fmt(s, q.Czz, "z^2"); s << ' ';
		fmt(s, q.Cxy, "xy"); s << ' ';
		fmt(s, q.Cxz, "xz"); s << ' ';
		fmt(s, q.Cyz, "yz"); s << ' ';
	}
	fmt(s, q.Cx, "x"); s << ' ';
	fmt(s, q.Cy, "y"); s << ' ';
	fmt(s, q.Cz, "z"); s << ' ';
	fmt(s, q.C); s << ' ';
	s << "=0";
	return s;
} /* operator << */
