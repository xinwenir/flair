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
 */

#include <math.h>

#include "obbox.h"
#include "quad.h"
#include "matrix4.h"

using namespace std;

OBBox::OBBox(const Array<Point>& samples)
{
	reset();
	ArrayIterator<Point> samples_iter(samples);

#if 0
	Point mean(0.0,0.0,0.0);
	Point m2(0.0,0.0,0.0);
	int n = 0;
#endif
	while (samples_iter) {
		const Point &v = samples_iter++;
#if 0
		n++;
		Point delta = v - mean;
		mean += (delta * (1.0/n));
		Point vm = v-mean;
		m2 += Point(	delta(0) * vm(0),
				delta(1) * vm(1),
				delta(2) * vm(2)
				);
#else
		add(v);
#endif
	}

}

/* --- add --- */
void OBBox::add(const double x, const double y, const double z)
{
	Point p = Point(x,y,z) - P;
	Point r(X*p, Y*p, Z*p);

	if (r.x<lowPt.x)  lowPt.x  = r.x;
	if (r.y<lowPt.y)  lowPt.y  = r.y;
	if (r.z<lowPt.z)  lowPt.z  = r.z;
	if (r.x>highPt.x) highPt.x = r.x;
	if (r.y>highPt.y) highPt.y = r.y;
	if (r.z>highPt.z) highPt.z = r.z;
} // add

/* --- addRelative --- */
void OBBox::addRelative(const double x, const double y, const double z)
{
	if (x<lowPt.x)  lowPt.x = x;
	if (y<lowPt.y)  lowPt.y = y;
	if (z<lowPt.z)  lowPt.z = z;
	if (x>highPt.x) highPt.x = x;
	if (y>highPt.y) highPt.y = y;
	if (z>highPt.z) highPt.z = z;
} // addRelative

/* --- insideRelative --- */
bool OBBox::insideRelative(const double x, const double y, const double z) const
{
	if (x<lowPt.x)  return false;
	if (y<lowPt.y)  return false;
	if (z<lowPt.z)  return false;
	if (x>highPt.x) return false;
	if (y>highPt.y) return false;
	if (z>highPt.z) return false;
	return true;
} // insideRelative

/** --- vertex ---
 * FIXME: can be cached
 */
Point OBBox::vertex(const int i) const
{
	switch (i) {
		case 0:
			return P + X*low().x  + Y*low().y  +  Z*low().z;
		case 1:
			return P + X*high().x + Y*low().y  +  Z*low().z;
		case 2:
			return P + X*high().x + Y*high().y +  Z*low().z;
		case 3:
			return P + X*low().x  + Y*high().y +  Z*low().z;
		case 4:
			return P + X*low().x  + Y*low().y  +  Z*high().z;
		case 5:
			return P + X*high().x + Y*low().y  +  Z*high().z;
		case 6:
			return P + X*high().x + Y*high().y +  Z*high().z;
		case 7:
			return P + X*low().x  + Y*high().y +  Z*high().z;
		default:
			assert(0);
	}
	return Point();
} // vertex

/** --- relativeVertex ---
 * FIXME: can be cached
 */
Point OBBox::relativeVertex(const int i) const
{
	switch (i) {
		case 0:
			return Point( low().x,  low().y,  low().z);
		case 1:
			return Point(high().x,  low().y,  low().z);
		case 2:
			return Point(high().x, high().y,  low().z);
		case 3:
			return Point( low().x, high().y,  low().z);
		case 4:
			return Point( low().x,  low().y, high().z);
		case 5:
			return Point(high().x,  low().y, high().z);
		case 6:
			return Point(high().x, high().y, high().z);
		case 7:
			return Point( low().x, high().y, high().z);
		default:
			assert(0);
	}
	return Point();
} // relativeVertex

/** --- edge ---
 * FIXME: can be cached
 */
void OBBox::edge(const int i, int* a, int* b) const
{
	const static int A[] = { 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7 };
	const static int B[] = { 1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 7, 4 };
	assert(0<=i && i<12);
	*a = A[i];
	*b = B[i];
} // edge

static inline double lowerPoint(
		double low, double high,
		double* tmin, double* tmax)
{
	double mid = 0.5*(high+low);
	double len = high-low;
	double t;

	if (Abs(mid) < Abs(low) && Abs(mid) < Abs(high)) {
		t = mid;
		*tmin = -len/2.0;
		*tmax = len/2.0;
	} else
	if (Abs(low) < Abs(high)) {
		t = low;
		*tmin = 0.0;
		*tmax = len;
	} else {
		t = high;
		*tmin = -len;
		*tmax = 0.0;
	}
	return t;
} // lowerPoint

/* setSmallestVolume
 * set result to the bb with the smallest volume among a and b
 */
static inline void setSmallestVolume(OBBox *result, const OBBox& a, const OBBox& b) {
	double va = a.volume();
	double vb = b.volume();
	if (va>SMALL && vb>SMALL) {
		if (va < vb)
			*result = a;
		else
			*result = b;
	} else
	if (va>SMALL)
		*result = a;
	else
		*result = b;
} // setSmallestVolume

/* safePointInFace
 * helper function to find a point with smallest coordinates in an axis-aligned
 * face
 */
static inline void safePointInFace(
		double minx, double maxx,
		double miny, double maxy,
		double* x, double* y)
{
	double midx = 0.5*(minx+maxx);
	double midy = 0.5*(miny+maxy);

	if (Abs(midx) < Abs(minx) && Abs(midx) < Abs(maxx))
		*x = midx;
	else if (Abs(minx) < Abs(maxx))
		*x = minx;
	else
		*x = maxx;

	if (Abs(midy) < Abs(miny) && Abs(midy) < Abs(maxy))
		*y = midy;
	else if (Abs(miny) < Abs(maxy))
		*y = miny;
	else
		*y = maxy;
} // safePointInFace


/** --- edge ---
 * FIXME: can be cached
 */
void OBBox::edge(const int i, Point* p, Vector* d, double *tmin, double *tmax) const
{
	double t;

	switch (i) {
		case 0:
			// Vertices 0 -> 1
			*d = X;
			t = lowerPoint(low().x, high().x, tmin, tmax);
			p->set(t, low().y, low().z);
			break;
		case 1:
			// Vertices 1 -> 2
			*d = Y;
			t = lowerPoint(low().y, high().y, tmin, tmax);
			p->set(high().x, t, low().z);
			break;
		case 2:
			// Vertices 2 -> 3 (reversed)
			*d = X;
			t = lowerPoint(low().x, high().x, tmin, tmax);
			p->set(t, high().y, low().z);
			break;
		case 3:
			// Vertices 3 -> 0 (reversed)
			*d = Y;
			t = lowerPoint(low().y, high().y, tmin, tmax);
			p->set(low().x, t, low().z);
			break;
		case 4:
			// Vertices 0 -> 4
			*d = Z;
			t = lowerPoint(low().z, high().z, tmin, tmax);
			p->set(low().x, low().y, t);
			break;
		case 5:
			// Vertices 1 -> 5
			*d = Z;
			t = lowerPoint(low().z, high().z, tmin, tmax);
			p->set(high().x, low().y, t);
			break;
		case 6:
			// Vertices 2 -> 6
			*d = Z;
			t = lowerPoint(low().z, high().z, tmin, tmax);
			p->set(high().x, high().y, t);
			break;
		case 7:
			// Vertices 3 -> 7
			*d = Z;
			t = lowerPoint(low().z, high().z, tmin, tmax);
			p->set(low().x, high().y, t);
			break;
		case 8:
			// Vertices 4 -> 5
			*d = X;
			t = lowerPoint(low().x, high().x, tmin, tmax);
			p->set(t, low().y, high().z);
			break;
		case 9:
			// Vertices 5 -> 6
			*d = Y;
			t = lowerPoint(low().y, high().y, tmin, tmax);
			p->set(high().x, t, high().z);
			break;
		case 10:
			// Vertices 6 -> 7 (reversed)
			*d = X;
			t = lowerPoint(low().x, high().x, tmin, tmax);
			p->set(t, high().y, high().z);
			break;
		case 11:
			// Vertices 7 -> 4 (reversed)
			*d = Y;
			t = lowerPoint(low().y, high().y, tmin, tmax);
			p->set(low().x, t, high().z);
			break;
		default:
			assert(0);
	}
	*p = P + X*p->x + Y*p->y + Z*p->z;
	return;
} // edge

/* face
 * get normal and a point in face number i, 0 <= i <= 5
 * @return true if the face should be used
 *	false if the face lies on the infinite and should not be used in
 *	other computations */
bool OBBox::face(const int i, Point* p, Vector* n) const
{
	double x, y;
	switch (i) {
		case 0:
			if (low().x < -OBB_INF_BOUND) return false;
			*n = -X;
			safePointInFace(
					low().y, high().y,
					low().z, high().z,
					&x, &y);
			*p = Point(low().x, x, y);
			//p = lowVertex();
			break;
		case 1:
			if (low().y < -OBB_INF_BOUND) return false;
			*n = -Y;
			safePointInFace(
					low().x, high().x,
					low().z, high().z,
					&x, &y);
			*p = Point(x, low().y, y);
			//p = lowVertex();
			break;
		case 2:
			if (low().z < -OBB_INF_BOUND) return false;
			*n = -Z;
			safePointInFace(
					low().x, high().x,
					low().y, high().y,
					&x, &y);
			*p = Point(x, y, low().z);
			//p = lowVertex();
			break;
		case 3:
			if (high().x > OBB_INF_BOUND) return false;
			*n = X;
			safePointInFace(
					low().y, high().y,
					low().z, high().z,
					&x, &y);
			*p = Point(high().x, x, y);
			//p = highVertex();
			break;
		case 4:
			if (high().y > OBB_INF_BOUND) return false;
			*n = Y;
			safePointInFace(
					low().x, high().x,
					low().z, high().z,
					&x, &y);
			*p = Point(x, high().y, y);
			//p = highVertex();
			break;
		case 5:
			if (high().z > OBB_INF_BOUND) return false;
			*n = Z;
			safePointInFace(
					low().x, high().x,
					low().y, high().y,
					&x, &y);
			*p = Point(x, y, high().z);
			//p = highVertex();
			break;
		default:
			assert(0);
	}
	*p = P + X*p->x + Y*p->y + Z*p->z;
	return true;
} // face

/* --- transform --- */
void OBBox::transform(const Matrix4& matrix)
{
	Point v[8];

	for(int i=0; i<vertices(); i++) v[i] = vertex(i);

//	X = matrix.multVector(X); X.normalize();
//	Y = matrix.multVector(Y); Y.normalize();
//	Z = matrix.multVector(Z); Z.normalize();
	X = matrix * X; X.normalize();
	Y = matrix * Y; Y.normalize();
	Z = matrix * Z; Z.normalize();

	// Translation
	P = Point(matrix(0,3), matrix(1,3), matrix(2,3));

	invalidate();

	for(int i=0; i<vertices(); i++) add(v[i]);
} // transform

/* --- Union --- */
void OBBox::Union(const OBBox& b)
{
	if (!b.isValid()) {
		return;
	} else if (b.isInfinite()) {
		infinite();
		return;
	}

	OBBox result1(*this);
	OBBox result2(b);
	result1.invalidate();
	result2.invalidate();

	for (int i=0; i < vertices(); i++) {
		result1.add(vertex(i));
		result2.add(vertex(i));
	}
	for (int i=0; i < b.vertices(); i++) {
		result1.add(b.vertex(i));
		result2.add(b.vertex(i));
	}

	setSmallestVolume(this, result1, result2);
} // Union

/* --- Intersect --- */
void OBBox::Intersect(const OBBox& b)
{
	if (!b.isValid()) {
		invalidate();
		return;
	} else if(b.isInfinite()) {
		return;
	}

	// Init coordinate systems
	OBBox result1(*this);
	OBBox result2(b);

	// reset limits
	result1.invalidate();
	result2.invalidate();

	// Add vertices of b if inside *this
	int nPoints = 0;
	for (int i=0; i<vertices(); i++) {
		const Point p = b.vertex(i);
		if (inside(p)) {
			nPoints++;
			result1.add(p);
			result2.add(p);
		}
	}
	// If all the vertices of b were inside the bb is correct now
	if (nPoints == 8) {
		setSmallestVolume(this, result1, result2);
		return;
	}

	// Add vertices of (*this) if inside b
	nPoints = 0;
	for (int i=0; i < vertices(); i++) {
		const Point p = vertex(i);
		if (b.inside(p)) {
			nPoints++;
			result1.add(p);
			result2.add(p);
		}
	}
	// If all the vertices of b were inside the bb is correct now
	if (nPoints == 8) {
		setSmallestVolume(this, result1, result2);
		return;
	}

	// Add intersections of the edges of b with the faces of (*this)
	for (int i=0; i < b.edges(); i++) {
		Point  pos;
		Vector dir;
		double tmin, tmax;
		b.edge(i, &pos, &dir, &tmin, &tmax);

		Point p[2];

		int n = intersectWithSegment(pos, dir, tmin, tmax, p);
		for (int j=0; j<n; j++) {
			result1.add(p[j]);
			result2.add(p[j]);
		}
	}

	// Add intersections of the edges of (*this) with the faces of b
	for (int i=0; i < edges(); i++) {
		Point  pos;
		Vector dir;
		double tmin, tmax;
		edge(i, &pos, &dir, &tmin, &tmax);

		Point p[2];

		int n = b.intersectWithSegment(pos, dir, tmin, tmax, p);
		for (int j=0; j<n; j++) {
			result1.add(p[j]);
			result2.add(p[j]);
		}
	}

	setSmallestVolume(this, result1, result2);
	return;
} // Intersect

/** --- Difference ---
 * WARNING: This works only if we assume that b is a FULL or inside bounding box
 */
void OBBox::Difference(const OBBox& b)
{
	if (!b.isValid()) {
		return;
	} else if (b.isInfinite()) {
		invalidate();
		return;
	}

	// Init coordinate systems
	OBBox result(*this);
	int nPoints;

	// reset limits
	result.invalidate();

	// Add vertices of (*this) if outside b
	nPoints = 0;
	for (int i=0; i < vertices(); i++) {
		const Point p = vertex(i);
		if (!b.inside(p)) {
			nPoints++;
			result.add(p);
		}
	}
	// If all the vertices of b were inside there will be no intersection
	// with the edges and the bb is correct now
	if (nPoints == 8) {
		*this = result;
		return;
	}

	// Add intersections of the edges of (*this) with the faces of b
	for (int i=0; i < edges(); i++) {
		Point  pos;
		Vector dir;
		double tmin, tmax;
		edge(i, &pos, &dir, &tmin, &tmax);

		Point p[2];

		int n = b.intersectWithSegment(pos, dir, tmin, tmax, p);
		for (int j=0; j<n; j++) {
			result.add(p[j]);
		}
	}
	*this = result;
	return;
} // Difference

/* --- overlap --- */
bool OBBox::overlap(const OBBox& /*b*/, double /*eps*/)
{
	assert(0);

	// Should intersect the two bodies and find out if the volume is > 0.0

	return true;
} // overlap

/* --- IntersectPlane ---
 * @return true if intersected
 */
void OBBox::IntersectPlane(const Quad& plane, bool negative) {
	OBBox bb;
	// Find out if the bounding box vertices
	// fall outside of the plane
	int inside_vertices = 0;
	for (int i=0; i < vertices(); i++) {
		const Point &v = vertex(i);
		double II = plane(v);
		if ((!negative && II <= 0.0) || (negative && II >= 0.0)) {
			inside_vertices++;
			bb.add(v);
		}
	}

	if (inside_vertices == 0) {
		// assert(0);
		reset();
		return;
	}
	else
	if (inside_vertices == 8)
		return;

	if (!Eq0(plane.Cx, SMALL)) {
		double t;

		t = - (plane.Cy * low().y + plane.Cz * low().z + plane.C) / plane.Cx;
		if (InRange(low().x, t, high().x))
			bb.add(t, low().y, low().z);

		t = - (plane.Cy * high().y + plane.Cz * low().z + plane.C) / plane.Cx;
		if (InRange(low().x, t, high().x))
			bb.add(t, high().y, low().z);

		t = - (plane.Cy * low().y + plane.Cz * high().z + plane.C) / plane.Cx;
		if (InRange(low().x, t, high().x))
			bb.add(t, low().y, high().z);

		t = - (plane.Cy * high().y + plane.Cz * high().z + plane.C) / plane.Cx;
		if (InRange(low().x, t, high().x))
			bb.add(t, high().y, high().z);
	}

	if (!Eq0(plane.Cy, SMALL)) {
		double t;

		t = - (plane.Cx * low().x + plane.Cz * low().z + plane.C) / plane.Cy;
		if (InRange(low().y, t, high().y))
			bb.add(low().x, t, low().z);

		t = - (plane.Cx * high().x + plane.Cz * low().z + plane.C) / plane.Cy;
		if (InRange(low().y, t, high().y))
			bb.add(high().x, t, low().z);

		t = - (plane.Cx * low().x + plane.Cz * high().z + plane.C) / plane.Cy;
		if (InRange(low().y, t, high().y))
			bb.add(low().x, t, high().z);

		t = - (plane.Cx * high().x + plane.Cz * high().z + plane.C) / plane.Cy;
		if (InRange(low().y, t, high().y))
			bb.add(high().x, t, high().z);
	}

	if (!Eq0(plane.Cz, SMALL)) {
		double t;

		t = - (plane.Cx * low().x + plane.Cy * low().y + plane.C) / plane.Cz;
		if (InRange(low().z, t, high().z))
			bb.add(low().x, low().y, t);

		t = - (plane.Cx * high().x + plane.Cy * low().y + plane.C) / plane.Cz;
		if (InRange(low().z, t, high().z))
			bb.add(high().x, low().y, t);

		t = - (plane.Cx * low().x + plane.Cy * high().y + plane.C) / plane.Cz;
		if (InRange(low().z, t, high().z))
			bb.add(low().x, high().y, t);

		t = - (plane.Cx * high().x + plane.Cy * high().y + plane.C) / plane.Cz;
		if (InRange(low().z, t, high().z))
			bb.add(high().x, high().y, t);
	}

	low(bb.low());
	high(bb.high());
} // IntersectPlane

/** planeSegmentIntersection
 * Intersect a segment with a plane, return the intersection point in p
 */
static bool planeLineIntersection(
		const Point& pPos, const Vector& pNormal,
		const Point& sPos, const Vector& sDir, double *t)
{
	double B = pNormal*sDir;
	double D = pNormal * (pPos - sPos);
	if (Eq0(B, SMALL)) return false;
	*t = D / B;
	return true;
}

/* intersectWithPlane
 * @param ori: origin of the plane
 * @param dirX: first axis of the plane
 * @param dirY: second axis of the plane
 * @param maxX: maximum coordinate on X, plane goes form [0..maxX]
 * @param maxY: maximum coordinate on Y, plane goes form [0..maxY]
 * @returns 0 if no intersection, 1 if the whole plane is inside, 2 if overlap
 * */
int OBBox::intersectWithPlane(const Matrix4& plane, double minX, double maxX, double minY, double maxY) const
{
	int n = 0; // Number of segments that intersect the infinite plane
	int onplane = 0; // Number of those intersections that lie within the limits of the plane

	//Point pPos = plane * Vector::O;
	Point  pPos(plane(0,3), plane(1,3), plane(2,3));
	Vector pNormal(plane(0,2), plane(1,2), plane(2,2));

	for (int i=0; i<6; i++) {
		Vector eN; // Normal of plane
		Point  eP; // Point on plane
		double tmin, tmax;

		// FIXME: Edges at infinite coordinates??
		edge(i, &eP, &eN, &tmin, &tmax);

		double t;
		if (planeLineIntersection(pPos, pNormal, eP, eN, &t)) {
			// Is on edge?
			if (InRange(tmin, t, tmax)) {
				n++;
				Point point = eP + eN*t;
				// Check if point on the viewport
				Vector pDirX(plane(0,0), plane(0,1), plane(0,2));
				Vector dir = point - pPos;
				double d = dir * pDirX;
				if (InRange(minX, d, maxX)) {
					Vector pDirY(plane(1,0), plane(1,1), plane(1,2));
					d = dir * pDirY;
					if (InRange(minY, d, maxY)) {
						onplane++;
					}
				}
			}
		}
	}
	if (n == 0) return 0;
#if 0
	// FIXME: make smarter
	if (onplane == n) {
		// Check if inside or outside
		if (inside(pPos)) return 1;
		else return 0;
	}
#endif
	return 2;
}

/** @return intersection points of a vector segment with the faces of
 * the bounding box, it will return 0,1 or 2 Vectors */
int OBBox::intersectWithSegment(
		const Point& pos, const Vector& dir,
		double tmin, double tmax,
		Point p[2]) const
{
	int n = 0;

	for (int i=0; i<6; i++) {
		Vector fN; // Normal of plane
		Point  fP; // Point on plane
		if (!face(i, &fP, &fN)) continue;

		double t;
		if (planeLineIntersection(fP, fN, pos, dir, &t)) {
			if (InRange(tmin*(1.0+SMALL2), t, tmax*(1.0+SMALL2))) {
				Point point = (dir*t) + pos;
				if (inside(point)) {
					if (n==0) {
						p[n] = point;
						n++;
					} else
					if ((p[0]-point).length2() > TOOSMALL) {
						p[n] = point;
						n++;
						break;
					}
				}
			}
		}
	}
	return n;
} // intersectWithSegment

/* ----------- operator << ------------- */
ostream& operator << (ostream& out, const OBBox& obj)
{
	out << "OBB X: " << obj.X << ", Y: " << obj.Y << ", Z: " << obj.Z << ", P: " << obj.P << ", Interval: " << obj.low() << " - " << obj.high() ;
	return out;
} /* operator << */
