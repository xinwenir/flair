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

#include <ostream>
#include <assert.h>

#include "os.h"
#include "bbox.h"
#include "quad.h"
#include "matrix4.h"

using namespace std;

/* --- add --- */
void BBox::add(const double x, const double y, const double z)
{
	if (x<lowPt.x)  lowPt.x = x;
	if (y<lowPt.y)  lowPt.y = y;
	if (z<lowPt.z)  lowPt.z = z;
	if (x>highPt.x) highPt.x = x;
	if (y>highPt.y) highPt.y = y;
	if (z>highPt.z) highPt.z = z;
} // add

/* --- vertex ---
 * FIXME: can be cached
 */
Point BBox::vertex(const int i) const
{
	switch (i) {
		case 0:
			return low();
		case 1:
			return Point(high().x,  low().y,  low().z);
		case 2:
			return Point(high().x, high().y,  low().z);
		case 3:
			return Point( low().x, high().y,  low().z);
		case 4:
			return Point( low().x,  low().y,  high().z);
		case 5:
			return Point(high().x,  low().y,  high().z);
		case 6:
			return high();
		case 7:
			return Point( low().x, high().y,  high().z);
		default:
			assert(0);
	}
	return Point();
} // vertex

/* --- vertex ---
 * FIXME: can be cached
 */
void BBox::edge(const int i, int* a, int* b) const
{
	const static int A[] = { 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7 };
	const static int B[] = { 1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 7, 4 };
	assert(0<=i && i<12);
	*a = A[i];
	*b = B[i];
} // vertex

/* --- transform --- */
void BBox::transform(const Matrix4& matrix)
{
	Point	pt[2];
	Point	lp(MAX_REAL,MAX_REAL,MAX_REAL),
		hp(-MAX_REAL,-MAX_REAL,-MAX_REAL);

	pt[0] = low();
	pt[1] = high();
	reset();

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++)
			for (int k=0; k<2; k++)
				add(matrix.mult(pt[i].x, pt[j].y, pt[k].z));
} // transform

/* --- Union --- */
void BBox::Union(const BBox& b)
{
	if (!b.isValid()) return;
	lowPt.set( Min(low().x, b.low().x),
		Min(low().y, b.low().y),
		Min(low().z, b.low().z));
	highPt.set( Max(high().x, b.high().x),
		Max(high().y, b.high().y),
		Max(high().z, b.high().z));
} // Union

/* --- Intersect --- */
void BBox::Intersect(const BBox& b)
{
	if (!b.isValid()) {
		reset();
		return;
	}
	lowPt.set( Max(low().x, b.low().x),
		   Max(low().y, b.low().y),
		   Max(low().z, b.low().z));
	highPt.set( Min(high().x, b.high().x),
		    Min(high().y, b.high().y),
		    Min(high().z, b.high().z));
} // Intersect

/* --- Difference --- */
void BBox::Difference(const BBox& b)
{
	if (!b.isValid()) return;
	Point& alow  = low();
	Point& ahigh = high();
	const Point& blow  = b.low();
	const Point& bhigh = b.high();

	// X
	{
		double& amin = low().x;
		double& amax = high().x;
		const double& bmin = b.low().x;
		const double& bmax = b.high().x;
		if ( blow.y <= alow.y && bhigh.y >= ahigh.y &&
		     blow.z <= alow.z && bhigh.z >= ahigh.z) {
			if (bmin <= amin && bmax > amin)
				amin = bmax;
			if (bmax >= amax && bmin < amax)
				amax = bmin;
			if (amin >= amax) {
				reset();
				return;
			}
		}
	}

	// Y
	{
		double& amin = low().y;
		double& amax = high().y;
		const double& bmin = b.low().y;
		const double& bmax = b.high().y;
		if ( blow.x <= alow.x && bhigh.x >= ahigh.x &&
		     blow.z <= alow.z && bhigh.z >= ahigh.z) {
			if (bmin <= amin && bmax > amin)
				amin = bmax;
			if (bmax >= amax && bmin < amax)
				amax = bmin;
			if (amin >= amax) {
				reset();
				return;
			}
		}
	}

	// Z
	{
		double& amin = low().z;
		double& amax = high().z;
		const double& bmin = b.low().z;
		const double& bmax = b.high().z;
		if ( blow.x <= alow.x && bhigh.x >= ahigh.x &&
		     blow.y <= alow.y && bhigh.y >= ahigh.y) {
			if (bmin <= amin && bmax > amin)
				amin = bmax;
			if (bmax >= amax && bmin < amax)
				amax = bmin;
			if (amin >= amax) {
				reset();
				return;
			}
		}
	}
} // Difference

/* --- contains --- */
bool BBox::contains(const double x, const double y, const double z, const double eps) const
{
	double px = Abs(x);
	double py = Abs(y);
	double pz = Abs(z);

	double e = Max(px, Abs(low().x), (double)1.0) * eps;
	if (x < low().x-e )
		return false;

	e = Max(px, Abs(high().x), 1.0) * eps;
	if (x > high().x+e )
		return false;

	e = Max(py, Abs(low().y), 1.0) * eps;
	if (y < low().y-e )
		return false;

	e = Max(py, Abs(high().y), 1.0) * eps;
	if (y > high().y+e )
		return false;

	e = Max(pz, Abs(low().z), 1.0) * eps;
	if (z < low().z-e )
		return false;

	e = Max(pz, Abs(high().z), 1.0) * eps;
	if (z > high().z+e )
		return false;

	return true;
} // contains

/* --- overlap --- */
bool BBox::overlap(const BBox& b, double eps)
{
	if (low().x   > b.high().x+eps)
		return false;
	if (b.low().x > high().x+eps)
		return false;

	if (low().y   > b.high().y+eps)
		return false;
	if (b.low().y > high().y+eps)
		return false;

	if (low().z   > b.high().z+eps)
		return false;
	if (b.low().z > high().z+eps)
		return false;

	return true;
} // overlap

/* truncate a bbox with a plane
 * @return true if intersected
 */
bool BBox::intersectPlane(const Quad& plane, bool negative)
{
	BBox bb;

	assert(plane.plane);

	// Find out if the bounding box vertices
	// fall outside of the plane
	int inside_vertices = 0;
	for (int i=0; i < 8; i++) {
		const Point &v = vertex(i);
		double II = plane(v);
		if ((!negative && II<=0.0) || (negative && II>=0.0)) {
			inside_vertices++;
			bb.add(v);
		}
	}

	if (inside_vertices == 0) {
		reset();
		return false;
	} else
	if (inside_vertices == 8)
		return false;

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

	return true;
} // intersectPlane

/* --- intersect --- */
bool BBox::intersect(const double  x, const double  y, const double  z,
		  const double dx, const double dy, const double dz,
		  const double tmin, const double tmax) const
{
	double	t, xt, yt, zt;

	if (!isValid()) return true;
	if (contains(x,y,z)) return true;

	// First search for the minimum distance
	if (dz>0.0) {
		t = (low().z - z) / dz;
		if (InRange(tmin, t, tmax )) {
			xt = dx*t + x;
			if (InRange(low().x, xt, high().x)) {
				yt = dy*t + y;
				if (InRange(low().y, yt, high().y))
					return true;
			}
		}
	} else
	if (dz<0.0) {
		t = (high().z - z) / dz;
		if (InRange(tmin, t, tmax )) {
			xt = dx*t + x;
			if (InRange(low().x, xt, high().x)) {
				yt = dy*t + y;
				if (InRange(low().y, yt, high().y))
					return true;
			}
		}
	}

	if (dy>0.0) {
		t = (low().y - y) / dy;
		if (InRange(tmin, t, tmax )) {
			xt = dx*t + x;
			if (InRange(low().x, xt, high().x)) {
				zt = dz*t + z;
				if (InRange(low().z, zt, high().z))
					return true;
			}
		}
	} else
	if (dy<0.0) {
		t = (high().y - y) / dy;
		if (InRange(tmin, t, tmax )) {
			xt = dx*t + x;
			if (InRange(low().x, xt, high().x)) {
				zt = dz*t + z;
				if (InRange(low().z, zt, high().z))
					return true;
			}
		}
	}

	if (dx>0.0) {
		t = (low().x - x) / dx;
		if (InRange(tmin, t, tmax )) {
			yt = dy*t + y;
			if (InRange(low().y, yt, high().y)) {
				zt = dz*t + z;
				if (InRange(low().z, zt, high().z))
					return true;
			}
		}
	} else
	if (dx<0.0) {
		t = (high().x - x) / dx;
		if (InRange(tmin, t, tmax )) {
			yt = dy*t + y;
			if (InRange(low().y, yt, high().y)) {
				zt = dz*t + z;
				if (InRange(low().z, zt, high().z))
					return true;
			}
		}
	}

	return false;
} // intersect

/* --- operator * --- */
BBox operator * (const Matrix4& matrix, const BBox& bounds)
{
	Point	pt[2];
	Point	low(  MAX_REAL,  MAX_REAL,  MAX_REAL),
		high(-MAX_REAL, -MAX_REAL, -MAX_REAL);

	pt[0] = bounds.low();
	pt[1] = bounds.high();

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++)
			for (int k=0; k<2; k++) {
				Point p = matrix.mult(pt[i].x, pt[j].y, pt[k].z);
				low.set(Min(low.x, p.x),
					Min(low.y, p.y),
					Min(low.z, p.z));
				high.set(Max(high.x, p.x),
					 Max(high.y, p.y),
					 Max(high.z, p.z));
			}
	return BBox(low,high);
} /* operator* */

/* ----------- operator << ------------- */
ostream& operator << (ostream& out, const BBox& obj)
{
	out << "BBox: " << obj.low() << " - " << obj.high();
	return out;
} /* operator << */
