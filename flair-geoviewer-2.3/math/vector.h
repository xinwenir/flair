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

#ifndef __VECTOR_H
#define __VECTOR_H

#include <math.h>
#include <assert.h>
#include <ostream>

#include "os.h"
#include "bmath.h"

class Matrix4;

/* ================================= Vector2D =============================== */
class Vector2D {
public:
	double	x, y;			/** the components	*/
public:
	Vector2D(const double ax=0.0, const double ay=0.0) : x(ax), y(ay) {}
	Vector2D(const Vector2D& p) : x(p.x), y(p.y) {}
}; // Vector2D

/* ================================== Vector ================================ *
 * 3D vector
 * WARNING: all transformations for Point are different from a Vector
 *          in the vector they include only the rotation
 *          while for a point includes the translation as well
 */
class Vector {
public:
	double	x, y, z;		/** the components	*/

public:
static	const Vector	O;		/** Center of axis	*/
static	const Vector	Xo;		/** X-axis unit vector	*/
static	const Vector	Yo;		/** Y-axis unit vector	*/
static	const Vector	Zo;		/** Z-axis unit vector	*/
static	double		_epsilon;	/** comparison operations accuracy */

public:
	/** Constructors */
	Vector(const double ax=0.0, const double ay=0.0, const double az=0.0)
		: x(ax), y(ay), z(az) {}
	Vector(const Vector& p)
		: x(p.x), y(p.y), z(p.z) {}

	/** @return components by index */
	double	operator [] (int i) const {
			assert(i>=0 && i<3);
			switch (i) {
				case 0:  return x;
				case 1:  return y;
				default: return z;
			}
		}

	/** @return components by index */
	double&	operator [] (int i) {
			assert(i>=0 && i<3);
			switch (i) {
				case 0:  return x;
				case 1:  return y;
				default: return z;
			}
		}

	/** Assignment.
	 * @param p vector to assign values from
	 */
	Vector& operator = (const Vector& p) {
			x=p.x; y=p.y; z=p.z;
			return *this;
		}
	Vector& operator () (const double ax, const double ay, const double az) {
			x=ax; y=ay; z=az;
			return *this;
		}

	/** Set the components in cartesian coordinate system. */
	void	set(const double ax, const double ay, const double az)
			{ x=ax; y=ay; z=az; }
	void	set(int i, double a) {
			switch (i) {
				case 0: x = a; break;
				case 1: y = a; break;
				case 2: z = a; break;
			}
		}

	/** Set the vector directly in polar coordinates
	 * @param ma magnitude of vector
	 * @param ph azimuthal angle in radians
	 * @param th polar angle in radians
	 */
	void	polar(const double ma, const double ph, const double th) {
			double sf,cf,st,ct;
			sincos(ph, &sf, &cf);
			sincos(th, &st, &ct);
			x = ma*st*cf;
			y = ma*st*sf;
			z = ma*ct;
		}

	/** @return the azimuth angle. */
	double	phi() const {
			return Cmp0(x,_epsilon) && Cmp0(y,_epsilon)?
					0.0 : atan2(y,x);
		}

	/** @return the polar angle. */
	double	theta() const {
			return Cmp0(x,_epsilon) && Cmp0(y, _epsilon) && Cmp0(z, _epsilon)?
				0.0 : atan2(perp(),z);
		}

	/** @return cosine of the polar angle. */
	double	cosTheta() const {
			double ptot = length();
			return Cmp0(ptot, _epsilon)? 1.0 : z/ptot;
		}

	/** @return the magnitude squared
	 * (rho^2 in spherical coordinate system)
	 */
	double	length2()	const { return Sqr(x) + Sqr(y) + Sqr(z); }

	/** @return the magnitude (rho in spherical coordinate system) */
	double	length()	const { return sqrt(Sqr(x) + Sqr(y) + Sqr(z)); }

	/** Set phi keeping mag and theta constant
	 * @param ph azimuthal angle in radians
	 */
	void	phi(const double ph) {
			double ma = length();
			double th = theta();
			polar(ma, ph, th);
		}

	/** Set theta keeping mag and phi constant
	 * @param th polar angle in radians
	 */
	void	theta(const double th) {
			double ma = length();
			double ph = phi();
			polar(ma, ph, th);
		}

	/** Set magnitude keeping theta and phi constant
	 * @param ma magnitude
	 */
	void	mag(const double ma) {
			normalize();
			x *= ma;
			y *= ma;
			z *= ma;
		}

	/** @return the transverse component squared
	 * (R^2 in cylindrical coordinate system).
	 */
	double	perp2()	const { return Sqr(x) + Sqr(y); }

	/** @return the transverse component
	 * (R in cylindrical coordinate system).
	 */
	double	perp()	const { return sqrt(perp2()); }

	/** @return the transverse component
	 * w.r.t. given axis squared.
	 * @param p axis to project
	 */
	double	perp2(const Vector& p) const {
			double	tot = p.length2();
			return tot > 0.0 ? length2()-Sqr(dot(p))/tot : length2();
		}

	/** @return the transverse component w.r.t. given axis.
	 * @param p axis to project
	 */
	double	perp(const Vector& p) const { return sqrt(perp2(p)); }

	/** Compare for proximity 2 vectors
	 * @param p vector to compare to
	 * @param e proximity distance or accuracy
	 * @return true if vectors are closer than e
	 * @see eps()
	 */
	bool	isClose(const Vector& p, const double e=_epsilon) const {
			double dx = p.x - x;
			double dy = p.y - y;
			double dz = p.z - z;
			if (Abs(dx)>e*Abs(p.x+x) ||
			    Abs(dy)>e*Abs(p.y+y) ||
			    Abs(dz)>e*Abs(p.z+z)) return false;
			return Sqr(dx) + Sqr(dy) + Sqr(dz) <= Sqr(e*max());
		}
	bool operator == (const Vector& p)	const { return isClose(p);    }
	bool operator != (const Vector& p)	const { return !isClose(p);   }
	bool operator <  (const Vector& p)	const { return compare(this,&p)<0;  }
	bool operator <= (const Vector& p)	const { return compare(this,&p)<=0; }
	bool operator >  (const Vector& p)	const { return compare(this,&p)>0;  }
	bool operator >= (const Vector& p)	const { return compare(this,&p)>0;  }

	/** Addition
	 * @param p vector to add
	 */
	Vector& operator += (const Vector& p) {
			x += p.x;
			y += p.y;
			z += p.z;
			return *this;
		}

	/** Subtraction.
	 * @param p vector to add
	 */
	Vector& operator -= (const Vector& p) {
			x -= p.x;
			y -= p.y;
			z -= p.z;
			return *this;
		}

	/** Unary minus. */
	Vector operator - () const { return Vector(-x, -y, -z); }

	/** Negate vector */
	void	negate() {
			x = -x;
			y = -y;
			z = -z;
		}

	/** Scaling with real number
	 * @param a number to scale with
	 */
	Vector&	operator *= (const double a) {
			x *= a;
			y *= a;
			z *= a;
			return *this;
		}

	/** @return unit vector parallel to this. */
	Vector	unit() const {
			double	tot = length2();
			Vector	p(*this);
			return tot > 0.0 ? p *= (1.0/sqrt(tot)) : p;
		}

	/** Normalize Vector
	 * @return magnitude of vector
	 */
	double	normalize() {
			double	tot = length2();
			if (tot > 0.0) {
				tot = sqrt(tot);
				*this *= 1.0/tot;
			}
			return tot;
		}

	/** Vector orthogonal to this
	 * @return a best guess for an orthogonal vector to this
	 */
	Vector	orthogonal() const {
			double xx = Abs(x), yy = Abs(y), zz = Abs(z);
			if (xx < yy) {
				if (xx < zz)
					return Vector(0.0, z, -y);
				else
					return Vector(y, -x, 0.0);
			} else {
				if (yy < zz)
					return Vector(-z, 0.0, x);
				else
					return Vector(y, -x, 0.0);
			}
		}

	/* @return scalar product of this with p
	 * @param p vector to multiply
	 */
	double	dot(const Vector& p) const { return x*p.x + y*p.y + z*p.z; }

	/* @return cross product of this with p
	 * @param p vector to multiply
	 */
	Vector	cross(const Vector& p) const {
			return Vector(	 y*p.z - p.y*z,
					 z*p.x - p.z*x,
					 x*p.y - p.x*y);
		}

	/* @return The angle w.r.t. another 3-vector.
	 * @param q vector to find the angle with
	 */
	double	angle(const Vector& q) const {
			double ptot2 = length2()*q.length2();
			return ptot2 <= 0.0 ? 0.0 : acos(dot(q)/sqrt(ptot2));
		}

	/** Rotates the Vector around the x-axis. */
	void	rotateX(double);

	/** Rotates the Vector around the y-axis. */
	void	rotateY(double);

	/** Rotates the Vector around the z-axis. */
	void	rotateZ(double);

	/** Rotates reference frame from Uz to newUz (unit vector) */
	void	rotateUz(Vector&);

	/** Rotates around the axis specified by another Vector. */
	void	rotate(double, const Vector&);

	/** Transformation with a Rotation matrix. */
	Vector&	operator *= (const Matrix4&);
	Vector& transform(const Matrix4&);

	int	direction(const double acc) const;

	/** max for error comparison */
	double	max()	const	{ return Abs(x)+Abs(y)+Abs(z); }

	/** Set the accuracy of comparison operations */
static	void	eps(double e)		{ _epsilon = e;  }
static	double	eps()			{ return _epsilon; }
	/** Compare 2 vectors. The comparison is done for sorting purposes
	 * using the following order z->y->x.
	 * @param p vector to compare to
	 * @retval -1 if this is smaller than v
	 * @retval 0 if are equal within accuracy eps()
	 * @retval 1 if this is greater than v
	 *
	 * @see eps()
	 */
static	int	compare(const Vector* a, const Vector* b);
}; // Vector

#if 0
#ifndef __MATRIX4_H
#	include "matrix4.h"
#endif

inline Vector& Vector::operator *=(const Matrix4& m) {
	return *this = m * (*this);
} // operator *=

inline Vector& Vector::transform(const Matrix4 &m) {
	return *this = m * (*this);
} // transform
#endif

// Output to a stream.
std::ostream& operator << (std::ostream&, const Vector&);

// Addition of 3-vectors.
inline Vector operator + (const Vector& a, const Vector& b) {
	return Vector(a.x+b.x, a.y+b.y, a.z+b.z);
}

// Subtraction of 3-vectors.
inline Vector operator - (const Vector& a, const Vector& b) {
	return Vector(a.x-b.x, a.y-b.y, a.z-b.z);
}

// Cross product of 3-vectors.
inline Vector operator ^ (const Vector& a, const Vector& b) {
	return a.cross(b);
}

// Scalar dot product of 3-vectors.
inline double operator * (const Vector& a, const Vector& b) {
	return a.dot(b);
}

// Scaling of 3-vectors with a real number
inline Vector operator * (const Vector& p, double a) {
	return Vector(a*p.x, a*p.y, a*p.z);
}

inline Vector operator * (double a, const Vector& p) {
	return Vector(a*p.x, a*p.y, a*p.z);
}

Vector safePointInEdge(const Vector& a, const Vector& b);

#endif
