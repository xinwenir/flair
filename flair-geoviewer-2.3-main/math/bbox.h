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

#ifndef __BBOX_H
#define __BBOX_H

#include <ostream>
#include "eps.h"
#include "point.h"

class	Quad;
class	Matrix4;

/** BBox - bounding box class */
class BBox {
public:
	Point	lowPt;
	Point	highPt;

public:
	/** Empty bounding box */
	BBox()	{ reset(); } // Create a Non valid bounding box

	BBox(const double xmin, const double xmax,
	     const double ymin, const double ymax,
	     const double zmin, const double zmax)
		: lowPt(xmin,ymin,zmin), highPt(xmax,ymax,zmax) {}

	/** Create a bounding box from low&high point
	 * @param aLow low point
	 * @param aHigh low point
	 */
	BBox(const Point& aLow, const Point& aHigh)
		: lowPt(aLow), highPt(aHigh) {}

	/** Copy constuctor */
	BBox(const BBox& aBoundingBox)
		: lowPt(aBoundingBox.lowPt),
		  highPt(aBoundingBox.highPt) {}

	/** check if bounding box is valid
	 * @return true if lowPt < highPt
	 */
	bool	isValid() const {
			return ((lowPt.x <= highPt.x) &&
				(lowPt.y <= highPt.y) &&
				(lowPt.z <= highPt.z));
		}

	/** reset bounding box to invalid */
	void	reset() {
			lowPt.set(  INFINITE,  INFINITE,  INFINITE);
			highPt.set(-INFINITE, -INFINITE, -INFINITE);
		}

	/** set bounding box to infinite */
	void	infinite() {
			lowPt.set(-INFINITE, -INFINITE, -INFINITE);
			highPt.set(INFINITE,  INFINITE,  INFINITE);
		}

	/** @return low point */
	Point	low()	const{ return lowPt; }

	/** @return high point */
	Point	high()	const	{ return highPt; }

	/** @return low point */
	Point&	low()		{ return lowPt; }

	/** @return high point */
	Point&	high()		{ return highPt; }

	/** set low point
	 * @param point to be set as low
	 */
	void	low(const Point& point)		{ lowPt = point; }

	/** set high point
	 * @param point to be set as high
	 */
	void	high(const Point& point)	{ highPt = point; }

	/** @return size of bbox */
	Vector	size() const			{ return highPt - lowPt; }

	/** @return volume of bounding box */
	double	volume() const		{ Vector s=size(); return s.x*s.y*s.z; }

	/** Add point
	 * @param x,y,z to be added
	 */
	void	add(const double x, const double y, const double z);

	/** Add point
	 * @param p	point to be added
	 */
	void	add(const Point& p)	{ add(p.x, p.y, p.z); }

	/** vertices
	 * @return number of vertices
	 */
	int	vertices()	const { return 8; }

	/** vertex(i)
	 * @return one of the numbered vertices
	 */
	Point	vertex(const int i) const;

	/** edges
	 * @return number of edges
	 */
	int	edges()		const { return 12; }

	/** edge(i)
	 * @return one of the numbered edges
	 */
	void	edge(const int i, int* a, int* b) const;

	/** Transform the bounding box
	 * @param trans transformation to be used
	 */
	void	transform(const Matrix4& matrix);

	/** Test for overlapping 2 bounding boxes
	 * @param b second bounding box
	 * @param eps accuracy for testing
	 * @return true if the overlap within accuracy eps
	 */
	bool	overlap(const BBox& b, double eps=TOOSMALL);

	/** Test is a specified point lies inside a bounding box
	 * @param point to test
	 * @param eps accuracy for testing
	 * @return true is point is inside within accuracy eps
	 */
	bool	contains(const double x, const double y, const double z, const double eps=TOOSMALL) const;
	bool	contains(const Point& point, const double eps=TOOSMALL) const
			{ return contains(point.x, point.y, point.z, eps); }

	/** Create the union with another bounding box */
	void	Union(const BBox& b);

	/** Create the intersection with another bounding box */
	void	Intersect(const BBox& b);

	/** Create the difference with another bounding box
	 * WARNING: this works only if the bbox is identical to the body!
	 *          e.g untransformed RPP, XYP, YZP, XZP
	 */
	void	Difference(const BBox& b);

	/** Intersect with a plane (positive or negative)
	 * @return true if intersected
	 */
	bool	intersectPlane(const Quad& plane, bool negative);

	/** Check for intersection against a ray
	 * @param ray to check for intersection
	 * @return true is ray intersects the bounding box
	 */
	bool	intersect(const double  x, const double  y, const double  z,
			  const double dx, const double dy, const double dz,
			  const double tmin, const double tmax) const;

	// bounding box operations
	BBox& operator |= (const BBox& b) {
			Union(b);
			return *this;
		}

	BBox& operator += (const BBox& b) {
			Intersect(b);
			return *this;
		}

	BBox& operator -= (const BBox& b) {
			Difference(b);
			return *this;
		}
}; // BBox

/** transform a bounding box */
BBox operator * (const Matrix4&, const BBox&);

/** print on a stream a bounding box */
std::ostream& operator << (std::ostream& out, const BBox& obj);

#endif
