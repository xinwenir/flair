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

#ifndef __QUAD_H
#define __QUAD_H

#include <ostream>

#include "os.h"
#include "eps.h"
#include "vector.h"

class Matrix3;
class Matrix4;
class BBox;
class OBBox;

/**
 * Generalized quadratic surface
 * S(x,y,z) =   Cxx*x^2 + Cyy*y^2 + Czz*z^2
 *            + Cxy*xy  + Cxz*xz  + Cyz*yz
 *            + Cx*x    + Cy*y    + Cz*z    + C  = 0
 *
 * Represented in matrix format
 *          0       1      2      3
 *     0 / Cxx     Cxy/2  Cxz/2  Cx/2 \
 *     1 | Cxy/2   Cyy    Cyz/2  Cy/2 |
 *     2 | Cxz/2   Cyz/2  Czz    Cz/2 |
 *     3 \ Cx/2    Cy/2   Cz/2   C    /
 */
class Quad {
public:
	double	Cxx, Cyy, Czz;
	double	Cxy, Cxz, Cyz;
	double	Cx,  Cy,  Cz;
	double	C;
	bool	plane;			// all high orders are zero

protected:
	double	aCxx, aCyy, aCzz;
	double	aCxy, aCxz, aCyz;
	double	aCxyz0;

public:
	Quad() : Cxx(0.0), Cyy(0.0), Czz(0.0),
		 Cxy(0.0), Cxz(0.0), Cyz(0.0),
		 Cx(0.0),  Cy(0.0),  Cz(0.0), C(0.0), plane(false) {}

	Quad(const double cxx, const double cyy, const double czz,
	     const double cxy, const double cxz, const double cyz,
	     const double cx,  const double cy,  const double cz,
	     const double c)
		{ set(cxx, cyy, czz, cxy, cxz, cyz, cx,  cy,  cz, c); }

	void	set(const double cxx, const double cyy, const double czz,
		    const double cxy, const double cxz, const double cyz,
		    const double cx,  const double cy,  const double cz,
		    const double c);
	void	set(const double cx,  const double cy,  const double cz,
		    const double c);
	void	set(const Quad& q);
	void	set(const Matrix4& m);
	Quad& operator = (const Quad& q)	{ set(q); return *this; }

	/** evaluate at x,y with z=0 */
	double	operator ()(const double x, const double y) const {
			if (plane) return Cx*x + Cy*y + C;
			  else return	(Cxx*x + Cxy*y + Cx)*x + (Cyy*y + Cy)*y + C;
		}

	/** evaluate at x,y,z */
	double	operator ()(const double x, const double y, const double z) const {
			if (plane) return Cx*x + Cy*y + Cz*z + C;
			  else return	(Cxx*x + Cxy*y + Cxz*z + Cx)*x +
					(Cyy*y + Cyz*z + Cy)*y +
					(Czz*z + Cz)*z + C;
		}

	double	operator ()(const Vector& v) const { return (*this)(v.x, v.y, v.z); }

	double	f(const double x, const double y, const double z) const { return (*this)(x,y,z); }
	double	f(const Vector& v)	const	{ return (*this)(v); }

	/** evaluate accuracy at x,y with z=0 */
	double	acc(const double x, const double y, const double eps) const {
			if (plane) return eps*(Abs(x) + Abs(y) + aCxyz0);
			else {
				double ax = Abs(x);
				double ay = Abs(y);
				return eps*(
					x*x + y*y +
					/*2.0*(*/aCxx*ax + aCyy*ay/*)*/ +
					ax*ay +
					aCxy*(ax + ay) +
					ax + ay +
					aCxyz0);
			}
		}

	/** evaluate accuracy at x,y,z */
	double	acc(const double x, const double y, const double z, const double eps) const {
			if (plane) return eps*(Abs(x) + Abs(y) + Abs(z) + aCxyz0);
			else {
				double ax = Abs(x);
				double ay = Abs(y);
				double az = Abs(z);
				return eps*(
					x*x + y*y + z*z +
					/*2.0*(*/aCxx*ax + aCyy*ay + aCzz*az/*)*/ +
					ax*ay + ax*az + ay*az +
					aCxy*(ax + ay) +
					aCxz*(ax + az) +
					aCyz*(ay + az) +
					ax + ay + az +
					aCxyz0);
			}
		}

	/** evaluate gradient at x,y with z=0 */
	void	grad(const double x, const double y, double *gx, double *gy, double *gz) const {
			if (plane) {
				*gx = Cx;
				*gy = Cy;
				*gz = Cz;
			} else {
				*gx = 2.0*Cxx*x + Cxy*y + Cx;
				*gy = 2.0*Cyy*y + Cxy*x + Cy;
				*gz = Cxz*x + Cyz*y + Cz;
			}
		}
	/** evaluate gradient at x,y,z */
	void	grad(const double x, const double y, const double z,
		         double *gx,     double *gy,     double *gz) const {
			if (plane) {
				*gx = Cx;
				*gy = Cy;
				*gz = Cz;
			} else {
				*gx = 2.0*Cxx*x + Cxy*y + Cxz*z + Cx;
				*gy = 2.0*Cyy*y + Cxy*x + Cyz*z + Cy;
				*gz = 2.0*Czz*z + Cxz*x + Cyz*y + Cz;
			}
		}
	Vector	grad(const Vector& r) const {
			double gx, gy, gz;
			grad(r.x,r.y,r.z, &gx,&gy,&gz);
			return Vector(gx,gy,gz);
		}

	double	dist1(const double x, const double y, const double z) const;
	double	dist2(const double x, const double y, const double z) const;
	double	adist(const double x, const double y, const double z) const;
	int	intersect(const double  x, const double  y, const double  z,
			  const double dx, const double dy, const double dz,
			  double *tmin, double *tmax) const;
	void	translate(const double dx, const double dy, const double dz);
	void	translate(const Vector& v)	{ translate(v.x, v.y, v.z); }
	void	matrix3(Matrix3 *m) const;
	void	matrix(Matrix4 *m) const;
	void	transform(const Matrix4& m);
	void	negate();

	void	normalize();

	// bounding box
	BBox	bbox() const;
	OBBox	obbox() const;
private:
	void	computeAbs();
	void	chk4Plane(const double eps=TOOSMALL);
}; // Quad

std::ostream& operator << (std::ostream&, const Quad&);

#endif
