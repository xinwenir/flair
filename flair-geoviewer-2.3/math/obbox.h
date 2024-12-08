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
 * Date:	09-Jan-2008
 */

#ifndef __OBBOX_H
#define __OBBOX_H

#include <iostream>

#include "os.h"
#include "eps.h"
#include "array.h"
#include "point.h"
#include "matrix3.h"

#define OBB_INF_BOUND (INFINITE/10.0)

class	Matrix4;
class   Quad;

/** OBBox - oriented bounding box class */
class OBBox {
public:
	Point	P; // Position in space
	Vector	X;
	Vector	Y;
	Vector	Z;
	Point	lowPt;
	Point	highPt;

public:
	/** Empty bounding box */
	OBBox()	{ reset(); } // Create a Non valid bounding box

	/** Copy constructor */
	OBBox(const OBBox& aBoundingBox)
		: P(aBoundingBox.P),
		  X(aBoundingBox.X),
		  Y(aBoundingBox.Y),
		  Z(aBoundingBox.Z),
		  lowPt( aBoundingBox.lowPt),
		  highPt(aBoundingBox.highPt)	{}

	/** Create an approximated bounding box from a list of vertices
	 */
	OBBox(const Array<Point>& vertices);

	/** reset bounding box to invalid */
	void	reset() {
			P = Vector::O;
			X = Vector::Xo;
			Y = Vector::Yo;
			Z = Vector::Zo;
			invalidate();
		}

	/** check if bounding box is valid
	 * @return true if lowPt < highPt
	 */
	bool	isValid() const {
			return ((lowPt.x < highPt.x) &&
				(lowPt.y < highPt.y) &&
				(lowPt.z < highPt.z));
		}

	/** return true if any of the limits is larger than OBB_INF_BOUND */
	bool	isInfinite() const {
			return ((lowPt.x < -OBB_INF_BOUND) &&
				(lowPt.y < -OBB_INF_BOUND) &&
				(lowPt.z < -OBB_INF_BOUND) &&
				(highPt.x > OBB_INF_BOUND) &&
				(highPt.y > OBB_INF_BOUND) &&
				(highPt.z > OBB_INF_BOUND));
		}

	/** reset limits */
	void invalidate() {
		lowPt.set(   INFINITE,  INFINITE,  INFINITE);
		highPt.set( -INFINITE, -INFINITE, -INFINITE);
	}

	/** set bounding box to infinite */
	void	infinite() {
			reset();
			lowPt.set(  -INFINITE, -INFINITE, -INFINITE);
			highPt.set(  INFINITE,  INFINITE,  INFINITE);
		}

	/** @return low point */
	Point	low()	const { return lowPt; }

	/** @return high point */
	Point	high()	const { return highPt; }

	/** @return low point */
	Point	&low() { return lowPt; }

	/** @return high point */
	Point	&high() { return highPt; }

	/** @return center */
	Point	center() const {
			Point c = 0.5*(lowPt + highPt);
			return P + X*c.x + Y*c.y + Z*c.z;
		}


	/** set low point
	 * @param point to be set as low
	 */
	void	low(const Point& point)		{ lowPt = point; }

	/** set high point
	 * @param point to be set as high
	 */
	void	high(const Point& point)	{ highPt = point; }

	/** Add point
	 * @param x,y,z to be added
	 */
	void	add(const double x, const double y, const double z);

	/** Add point
	 * @param p	point to be added
	 */
	void	add(const Point& p)		{ add(p.x, p.y, p.z); }

	/** Add point relative to the box system
	 * @param x,y,z to be added
	 */
	void	addRelative(const double x, const double y, const double z);
	void	addRelative(const Point& p)	{ add(p.x, p.y, p.z); }

	/** inside
	 * @param x,y,z to be checked
	 * @returns true if inside
	 */
	bool	inside(const double x, const double y, const double z) const
			{ return inside(Point(x,y,z)); }

	/** inside
	 * @param point to be checked
	 * @returns true if inside
	 */
	bool	inside(const Point& r) const;

	/** insideRay
	 * @param p position of the ray
	 * @param d direction
	 * @param t current displacement so that the point is at p+d*t
	 * @returns true if inside
	 */
	bool insideRay(const Point& p, const Vector& d, const double t) const;
	bool insideRay(const double u, const double v, const double w,
			const double du, const double dv, const double dw,
			const double t) const
	{
		return insideRay(Point(u,v,w), Vector(du,dv,dw), t);
	}

	/** insideFast
	 * It is less precise with infinite coordinates
	 * @param x,y,z to be checked
	 * @returns true if inside
	 */
	bool insideFast(const double u, const double v, const double w) const;

	/** insideRelative
	 * @param x,y,z to be checked
	 * @returns true if inside
	 */
	bool	insideRelative(const double x, const double y, const double z) const;

	/** vertices
	 * @return number of vertices
	 */
	int	vertices()	const { return 8; }

	/** vertex(i)
	 * @return one of the numbered vertices
	 */
	Point	vertex(const int i) const;
	Point	lowVertex() const { return vertex(0); }
	Point	highVertex() const { return vertex(6); }

	/** relativeVertex(i)
	 * @return one of the numbered vertices relative to P
	 */
	Point	relativeVertex(const int i) const;

	/** edges
	 * @return number of edges
	 */
	int	edges()		const { return 12; }

	/** edge(i)
	 * @return one of the numbered edges
	 */
	void	edge(const int i, int* a, int* b) const;
	void	edge(const int i, Point* p, Vector* d, double *tmin, double *tmax) const;

	/** @return parameters of numbered face i, 0 <= i <= 5 */
	bool face(const int i, Point *p, Vector *n) const;

	/** Transform the bounding box
	 * @param matrix transformation matrix to be used
	 */
	void	transform(const Matrix4& matrix);

	/** Test for overlapping 2 bounding boxes
	 * @param b second bounding box
	 * @param eps accuracy for testing
	 * @return true if the overlap within accuracy eps
	 */
	bool	overlap(const OBBox& b, double eps=TOOSMALL);

	/** Create the union with another bounding box
	 */
	void	Union(const OBBox& b);

	/** Create the intersection with another bounding box
	 */
	void	Intersect(const OBBox& b);

	/** Create the difference with another FULL bounding box
	 * normally it has to be the inside bounding box
	 */
	void	Difference(const OBBox& b);

	/** Intersect with a plane (positive or negative)
	 * @return true if intersected
	 */
	void	IntersectPlane(const Quad& plane, bool negative);

	/** @return volume of bounding box */
	double	volume() const {
			if (isValid())
				return  (high().x-low().x) *
					(high().y-low().y) *
					(high().z-low().z);
			else
				return 0.0;
		}

	// bounding box operations
	OBBox& operator |= (const OBBox& b) {
			Union(b);
			return *this;
		}

	OBBox& operator = (const OBBox& b) {
			if (this != &b) {
				this->P = b.P;
				this->X = b.X;
				this->Y = b.Y;
				this->Z = b.Z;
				this->lowPt = b.lowPt;
				this->highPt = b.highPt;
			}
			return *this;
		}

	OBBox& operator += (const OBBox& b) {
			Intersect(b);
			return *this;
		}

	OBBox& operator -= (const OBBox& b) {
			Difference(b);
			return *this;
		}

	/** @return intersection points of a vector segment with the faces of
	 * the bounding box, it will return 0,1 or 2 Vectors */
	int intersectWithSegment(const Point &pos, const Vector &dir, double tmin, double tmax, Point p[2]) const;
	int intersectWithPlane(const Matrix4& plane, double minX, double maxX, double minY, double maxY) const;

}; // OBBox

/* --- inside --- */
inline bool OBBox::inside(const Point& r) const
{
	//if (isInfinite()) return true;

	Point rt = r - P;
	double x = rt*X;
	if (lowPt.x > -OBB_INF_BOUND && ((x-lowPt.x) < -SMALL2)) return false;
	if (highPt.x < OBB_INF_BOUND && ((x-highPt.x) > SMALL2)) return false;

	double y = rt*Y;
	if (lowPt.y > -OBB_INF_BOUND && ((y-lowPt.y) < -SMALL2)) return false;
	if (highPt.y < OBB_INF_BOUND && ((y-highPt.y) > SMALL2)) return false;

	double z = rt*Z;
	if (lowPt.z > -OBB_INF_BOUND && ((z-lowPt.z) < -SMALL2)) return false;
	if (highPt.z < OBB_INF_BOUND && ((z-highPt.z) > SMALL2)) return false;

	return true;
} /* inside */

/* --- insideRay --- */
inline bool OBBox::insideRay(const Point& p, const Vector& d, const double t) const
{
	Point rpos = p + (d*t);
	Vector rdir(d*X, d*Y, d*Z);
	double x = (rpos.x-P.x)*X.x + rpos.y*X.y + rpos.z*X.z;
	if (lowPt.x > -OBB_INF_BOUND &&
	   (x <= lowPt.x) &&
	   (((x-lowPt.x) < -SMALL2) || rdir.x <= 0))
		return false;
	if (highPt.x < OBB_INF_BOUND &&
	   (x >= highPt.x) &&
	   (((x-highPt.x) > SMALL2) || rdir.x >= 0))
		return false;

	double y = rpos.x*Y.x + (rpos.y-P.y)*Y.y + rpos.z*Y.z;
	if (lowPt.y > -OBB_INF_BOUND &&
	   (y <= lowPt.y) &&
	   (((y-lowPt.y) < -SMALL2) || rdir.y <= 0))
		return false;
	if (highPt.y < OBB_INF_BOUND &&
	   (y >= highPt.y) &&
	   (((y-highPt.y) > SMALL2) || rdir.y >= 0))
		return false;

	double z = rpos.x*Z.x + rpos.y*Z.y + (rpos.z-P.z)*Z.z;
	if (lowPt.z > -OBB_INF_BOUND &&
	   (z <= lowPt.z) &&
	   (((z-lowPt.z) < -SMALL2) || rdir.z <= 0))
		return false;
	if (highPt.z < OBB_INF_BOUND &&
	   (z >= highPt.z) &&
	   (((z-highPt.z) > SMALL2) || rdir.z >= 0))
		return false;

	return true;
} /* inside */

inline bool OBBox::insideFast(const double u, const double v, const double w) const
{
#if 1
	double x = (u-P.x)*X.x + v*X.y + w*X.z;
	if (x <  lowPt.x) return false;
	if (x > highPt.x) return false;

	double y = u*Y.x + (v-P.y)*Y.y + w*Y.z;
	if (y <  lowPt.y) return false;
	if (y > highPt.y) return false;

	double z = u*Z.x + v*Z.y + (w-P.z)*Z.z;
	if (z <  lowPt.z) return false;
	if (z > highPt.z) return false;
#endif
#if 0
	double x = (u-P.x)*X.x + v*X.y + w*X.z;
	if ((x-lowPt.x) < -SMALL2) return false;
	if ((x-highPt.x) > SMALL2) return false;

	double y = u*Y.x + (v-P.y)*Y.y + w*Y.z;
	if ((y-lowPt.y) < -SMALL2) return false;
	if ((y-highPt.y) > SMALL2) return false;

	double z = u*Z.x + v*Z.y + (w-P.z)*Z.z;
	if ((z-lowPt.z) < -SMALL2) return false;
	if ((z-highPt.z) > SMALL2) return false;
#endif

	return true;
}

/** print on a stream a bounding box */
std::ostream& operator << (std::ostream& out, const OBBox& obj);

#endif
