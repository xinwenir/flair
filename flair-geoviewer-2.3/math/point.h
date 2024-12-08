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


#ifndef __POINT_H
#define __POINT_H

#include <cstdint>
#include "vector.h"

class Matrix4;

/* ================================= Point2D ================================ */
class Point2D : public Vector2D {
public:
	Point2D(const double ax=0.0, const double ay=0.0) : Vector2D(ax,ay) {}
	Point2D(const Point2D& p) : Vector2D(p) {}
}; // Point2D

/* =================================== Point ================================ *
 * 3D Point: Overloaded class of Vector defining a point in space
 *           in contrast of a vector that typically is used to
 *           define a direction and magnitude
 *
 * WARNING: all transformations for Point are different from a Vector
 *          in the vector they include only the rotation
 *          while for a point includes the translation as well
 */
class Point : public Vector {
public:
	// Constructor
	Point(const double ax=0.0, const double ay=0.0, const double az=0.0)
		: Vector(ax, ay, az) {}

        //PS
        Point(const Point& p) { x=p.x; y=p.y; z=p.z;}
  
	// Copy constructor
	Point(const Vector& p) : Vector(p) {}

	// Assignment
	Point&	operator =(const Point& p) { x=p.x; y=p.y; z=p.z; return *this; }

	// Assignment of Vector and classes derived from it (Vector, Normal)
	Point&	operator =(const Vector &p) { x=p.x; y=p.y; z=p.z; return *this; }

	// Distance squared to the origin
	double	distance2() const { return length2(); }

	// Distance to the origin
	double	distance() const { return sqrt(distance2()); }

	// Distance squared to the point
	double	distance2(const Point &p) const {
			double delx = p.x-x;
			double dely = p.y-y;
			double delz = p.z-z;
			return delx*delx + dely*dely + delz*delz;
		}

	// Distance to the point
	double	distance(const Point &p) const { return sqrt(distance2(p)); }

	// Transformation
//	Point& transform(const Matrix4 &mat) { return *this = mat * (*this); }
}; // class Point

//#ifndef __MATRIX4_H
//#	include "matrix4.h"
//#endif
//inline Point& Point::transform(const Matrix4 &mat) {
//	return *this = mat * (*this);
//}

// Overloaded operations on Points

// point + point = point	(A bit brain twisted it should not exist)
inline Point operator + (const Point& a, const Point& b) {
	return Point(a.x+b.x, a.y+b.y, a.z+b.z);
}

// point + vector = point
inline Point operator + (const Point& a, const Vector& b) {
	return Point(a.x+b.x, a.y+b.y, a.z+b.z);
}

// vector + point = point
inline Point operator + (const Vector& a, const Point& b) {
	return Point(a.x+b.x, a.y+b.y, a.z+b.z);
}

// point - point = vector
inline Vector operator - (const Point& a, const Point& b) {
	return Vector(a.x-b.x, a.y-b.y, a.z-b.z);
}

// point - vector = point
inline Point operator - (const Point& a, const Vector& b) {
	return Point(a.x-b.x, a.y-b.y, a.z-b.z);
}

// vector - point = vector	(A bit brain twisted it should not exist)
inline Vector operator - (const Vector& a, const Point& b) {
	return Vector(a.x-b.x, a.y-b.y, a.z-b.z);
}

// f * point = point
inline Point operator * (const Point& p, double a) {
	return Point(a*p.x, a*p.y, a*p.z);
}

// point * f = point
inline Point operator * (double a, const Point& p) {
	return Point(a*p.x, a*p.y, a*p.z);
}

#endif
