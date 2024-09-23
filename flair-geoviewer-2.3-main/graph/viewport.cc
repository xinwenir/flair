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

#include "viewport.h"

#include "os.h"
#include "painter.h"

using namespace std;

#define MINSIZE		10
#define MAXSIZE		3000
#define MINSCALE	1e-7
#define MAXSCALE	1e14

/** ViewPort */
ViewPort::ViewPort(int w, int h)
{
	_zoom    = 1.0;
	aspect   = 1.0;
	_minu = _minv = -100.0;
	_maxu = _maxv =  100.0;
	chkminu = chkminv = _minu;
	chkmaxu = chkmaxv = _maxu;
	Uo = Vo = 0.0;
	scaleU = scaleV = 1.0;
	_fishEye = false;

	_matrix.identity();
	_invMatrix.identity();

	projection = Projection_Orthographic;
	fov(RAD(60.0));

	init(w,h);
} // ViewPort

/** init */
void ViewPort::init(int w, int h)
{
	_width  = w;
	_height = h;
	width2  = w/2;
	height2 = h/2;
	zoom(_zoom);
} // init

/**
 * Viewport window
 * @param ex	extend in X
 * @param ey	extend in Y
 */
void ViewPort::window(const double minx, const double miny,
		      const double maxx, const double maxy)
{
	_minu = Min(minx, maxx);
	_maxu = Max(minx, maxx);
	_minv = Min(miny, maxy);
	_maxv = Max(miny, maxy);

	double xm = 0.5*(_minu+_maxu);
	double ym = 0.5*(_minv+_maxv);
	double dx = maxx - minx;
	double dy = maxy - miny;

	if (dx<MINWINDOW) dx = MINWINDOW;
	else
	if (dx>MAXWINDOW) dx = MAXWINDOW;

	if (dy<MINWINDOW) dy = MINWINDOW;
	else
	if (dy>MAXWINDOW) dy = MAXWINDOW;

	_minu = xm - dx/2.0;
	_maxu = xm + dx/2.0;

	_minv = ym - dy/2.0;
	_maxv = ym + dy/2.0;

	zoom(_zoom);
	offset(xm,ym);

	ulow.set( 0.0, 0.0, 0.0,-0.5, 0.0,  _minu);
	uhigh.set(0.0, 0.0, 0.0, 0.5, 0.0, -_maxu);
	vlow.set( 0.0, 0.0, 0.0, 0.0,-0.5,  _minv);
	vhigh.set(0.0, 0.0, 0.0, 0.0, 0.5, -_maxv);
	ulow.parametric();
	uhigh.parametric();
	vlow.parametric();
	vhigh.parametric();

	if (Abs(_minu) < 1.0)
		chkminu = _minu - 1e-9;
	else
		chkminu = _minu - Abs(_minu)*1e-8;

	if (Abs(_maxu) < 1.0)
		chkmaxu = _maxu + 1e-9;
	else
		chkmaxu = _maxu + Abs(_maxu)*1e-8;

	if (Abs(_minv) < 1.0)
		chkminv = _minv - 1e-9;
	else
		chkminv = _minv - Abs(_minv)*1e-8;

	if (Abs(_maxv) < 1.0)
		chkmaxv = _maxv + 1e-9;
	else
		chkmaxv = _maxv + Abs(_maxv)*1e-8;
} // window

/** find new limits for the window
 * @return true if the window doesn't contain the previous one
 */
bool ViewPort::invalidWindow()
{
	double r = i2u(0);
	if (r < chkminu || r > chkmaxu) return true;
	r = i2u(width());
	if (r < chkminu || r > chkmaxu) return true;
	r = j2v(height());
	if (r < chkminv || r > chkmaxv) return true;
	r = j2v(0);
	if (r < chkminv || r > chkmaxv) return true;
	return false;
} // invalidWindow

/** set new limits for the window for z zoom
 * @param z	desired zoom level for the window
 */
void ViewPort::calcWindow(const double z)
{
	double minx = i2u(0);
	double maxx = i2u(width());
	double miny = j2v(height());
	double maxy = j2v(0);

	double w = maxx - minx;
	double h = maxy - miny;

	minx = Uo - w/2.0 * z;
	maxx = Uo + w/2.0 * z;

	miny = Vo - h/2.0 * z;
	maxy = Vo + h/2.0 * z;

	window(minx, miny, maxx, maxy);
	zoom(z);
} // calcWindow

/** set zoom
 * @param z	zoom factor
 */
void ViewPort::zoom(const double z)
{
	_zoom = z;
	scaleU = _zoom * (double)width()  / relWidth();
	if (scaleU < MINSCALE) {
		scaleU = MINSCALE;
		_zoom = scaleU / ((double)width() / relWidth());
	} else
	if (scaleU > MAXSCALE) {
		scaleU = MAXSCALE;
		_zoom = scaleU / ((double)width() / relWidth());
	}
	scaleV = scaleU / aspect;

	if (_fixFOV)
		calculateFocalLength();
	else
		calculateFOV();
	calculateFocalLength();
} // zoom

/** is point (x,y) inside window
 * @param u	u coordinate
 * @param v	v coordinate
 * @return true for inside, false for outside
 */
bool ViewPort::inside(double u, double v) const
{
	if (u < chkminu || u > chkmaxu) return false;
	if (v < chkminv || v > chkmaxv) return false;
	return true;
} // inside

/** matrix */
void ViewPort::matrix(const Matrix4& m)
{
	_matrix.copy(m);
	computeMatrices();
} // matrix

/** transform current matrix my multiplying with the new one
 */
void ViewPort::transform(const Matrix4& m)
{
	_matrix *= m;
	computeMatrices();
} // transform

/** set transformation matrix origin */
void ViewPort::origin(const double x, const double y, const double z)
{
	_matrix(0,3) = x;
	_matrix(1,3) = y;
	_matrix(2,3) = z;
	computeMatrices();
} // origin

/** get transformation origin
 * return origin of the center of the view screen!
 */
void ViewPort::origin(double *x, double *y, double *z) const
{
	*x = uv2x(Uo,Vo);
	*y = uv2y(Uo,Vo);
	*z = uv2z(Uo,Vo);
} // origin

/** prepare the inverse matrix and other matrices */
void ViewPort::computeMatrices()
{
	_matrix.fix();
	_invMatrix.inverse(_matrix);
	_invMatrix.fix();
} // computeMatrices

/** moveOriginTo0
 * Move view origin to center of "view"
 */
void ViewPort::moveOriginTo0()
{
	double u = Uo;
	double v = Vo;

	_matrix(0,3) += uv2dx(u,v);
	_matrix(1,3) += uv2dy(u,v);
	_matrix(2,3) += uv2dz(u,v);

	computeMatrices();

	// center view extends
	double eu = 0.5*(_maxu - _minu);
	double ev = 0.5*(_maxv - _minv);
	window(-eu,-ev,eu,ev);

	// move view origin to 0
	offset();
} // moveOriginTo0

/** @return true if conic is inside viewport inside
 */
bool ViewPort::inside(const Conic& conic) const
{
	int i,n;
	double x,y;
	Point2D pts[4];

	// Check if an intersection is inside the extends
	n = conic.intersect(ulow, pts);
	for (i=0; i<n; i++)
		if (insideV(pts[i].y))
			return true;

	n = conic.intersect(uhigh, pts);
	for (i=0; i<n; i++)
		if (insideV(pts[i].y))
			return true;

	n = conic.intersect(vlow, pts);
	for (i=0; i<n; i++)
		if (insideU(pts[i].x))
			return true;

	n = conic.intersect(vhigh, pts);
	for (i=0; i<n; i++)
		if (insideU(pts[i].x))
			return true;

	if (conic.type() == CONIC_ELLIPSE) {
		// Could lie inside without touching the borders
		// Check one point of the ellipse
		conic.getXY(0.0, &x, &y);
		if (inside(x, y))
			return true;
	}

	return false;
} // isInside

/** convert pixel coordinates to absolute xyz
 * @param i,j	pixel
 * @param x,y,z	return x,y,z coordinates
 */
void ViewPort::ij2xyz(const int i, const int j, double *x, double *y, double *z) const
{
	double u = i2u(i);
	double v = j2v(j);
	*x = uv2x(u,v);
	*y = uv2y(u,v);
	*z = uv2z(u,v);
} // ij2xyz

/** set camera fov
 * @param a	fov in radians
 */
void ViewPort::fov(const double a)
{
	if (a<=0) return;
	_fov = a;
	_fixFOV = true;
	_fishEye = _fov > PI;
	calculateFocalLength();
} // fov

/** set camera focal length
 * @param f	focal length
 */
void ViewPort::focalLength(const double f)
{
	if (f<=0) return;
	_focal = f;
	_fixFOV = false;
	calculateFOV();
} // focalLength

/** calculate the camera focal length based on fov and the window width */
void ViewPort::calculateFocalLength()
{
	if (_fishEye)
		_focal = 0.0;
	else
		_focal = relWidth() / 2.0 / _zoom / tan(_fov/2.0);
} // calculateFocalLength

/** calculate the camera fov based on focal length and the window width */
void ViewPort::calculateFOV()
{
	_fov = 2.0 * atan2(relWidth() / 2.0 / _zoom, _focal);
//	_fov = 2.0 * atan(relWidth() / 2.0 / _zoom / _focal);
} // calculateFOV

/** clip line coordinates to the appropriate clipping region
 * @param x1,y1	starting point of line
 * @param x2,y2	end point of line
 * @return true if a fraction of the line falls inside th clipping region
 */
bool ViewPort::clipLine(double *x1, double *y1, double *x2, double *y2) const
{
	union _OutCodeUnion ocu1, ocu2;
	bool   in;
	bool   out;
	bool   swap = false;

	/* Initialize 4-bit codes */
	ocu1.outcodes = 0;
	ocu1.ocs.code0 = (*x1 < chkminu);
	ocu1.ocs.code1 = (*y1 < chkminv);
	ocu1.ocs.code2 = (*x1 > chkmaxu);
	ocu1.ocs.code3 = (*y1 > chkmaxv);

	ocu2.outcodes = 0;
	ocu2.ocs.code0 = (*x2 < chkminu);
	ocu2.ocs.code1 = (*y2 < chkminv);
	ocu2.ocs.code2 = (*x2 > chkmaxu);
	ocu2.ocs.code3 = (*y2 > chkmaxv);

	in  = ((ocu1.outcodes | ocu2.outcodes) == 0);
	out = ((ocu1.outcodes & ocu2.outcodes) != 0);

	int i = 0;	// to avoid numerical precision loop few times
	while (!out && !in && i<5) {
		if (ocu1.outcodes==0) {		/* Swap endpoints if necessary so */
			Swap(*x1,*x2);
			Swap(*y1,*y2);
			Swap(ocu1.outcodes,ocu2.outcodes);
			swap = !swap;
		}

		if (ocu1.ocs.code0) {		/* Clip left */
			*y1 += ((*y2-*y1) * (_minu-*x1))/(*x2-*x1);
			*x1 = _minu;
		}
		else
		if (ocu1.ocs.code1) {		/* Clip below */
			*x1 += ((*x2-*x1) * (_minv-*y1))/(*y2-*y1);
			*y1 = _minv;
		}
		else
		if (ocu1.ocs.code2) {		/* Clip right */
			*y1 += ((*y2-*y1) * (_maxu-*x1))/(*x2-*x1);
			*x1 = _maxu;
		}
		else
		if (ocu1.ocs.code3) {		/* Clip above */
			*x1 += ((*x2-*x1) * (_maxv-*y1))/(*y2-*y1);
			*y1 = _maxv;
		}

		/* update for (*x1, *y1) */
		ocu1.outcodes = 0;
		ocu1.ocs.code0 = (*x1 < chkminu);
		ocu1.ocs.code1 = (*y1 < chkminv);
		ocu1.ocs.code2 = (*x1 > chkmaxu);
		ocu1.ocs.code3 = (*y1 > chkmaxv);


		in  = ((ocu1.outcodes | ocu2.outcodes) == 0); /* update */
		out = ((ocu1.outcodes & ocu2.outcodes) != 0); /* 4-bit codes */
		i++;
	}
	if (swap) {
		// Swap back to original order
		Swap(*x1,*x2);
		Swap(*y1,*y2);
	}
	return in;
}  /* clipLine */

/** clipLine3D */
bool ViewPort::clipLine3D(const Point& a, const Point& b, Point* ra, Point* rb) const
{
	if (projection == Projection_Orthographic) {
		// No w clipping
		xyz2uvw3D(a, ra);
		xyz2uvw3D(b, rb);
	} else {
		// Clip along w
		//    Perspective: below the viewer
		//    Combo:       below the screen plane
		xyz2uvw(a, ra);
		xyz2uvw(b, rb);

		if (projection != Projection_Orthographic) {
			// Apply offset for the perspective calculation
			ra->x -= Uo;
			ra->y -= Vo;
			rb->x -= Uo;
			rb->y -= Vo;
		}

		// The viewport is centered on the projection plane
		// ... move to viewer
		ra->z -= _focal;
		rb->z -= _focal;

		if (ra->z >= -Vector::eps()) {
			if (rb->z >= -Vector::eps()) return false;
			// clip wa
			double f = -_focal/10.0;
			ra->x = rb->x + (ra->x - rb->x)/(ra->z - rb->z)*(f - rb->z);
			ra->y = rb->y + (ra->y - rb->y)/(ra->z - rb->z)*(f - rb->z);
			ra->z = f;
		} else
		if (rb->z >= -Vector::eps()) {
			// clip wb
			double f = -_focal/10.0;
			rb->x = ra->x + (rb->x - ra->x)/(rb->z - ra->z)*(f - ra->z);
			rb->y = ra->y + (rb->y - ra->y)/(rb->z - ra->z)*(f - ra->z);
			rb->z = f;
		}

		if (Eq0(ra->z, Vector::eps()))
			ra->x = ra->y = 1e10;	// FIXME should correct for sign
		else {
			ra->x *= -_focal / ra->z;
			ra->y *= -_focal / ra->z;
		}
		if (Eq0(rb->z, Vector::eps()))
			rb->x = rb->y = 1e10;	// FIXME should correct for sign
		else {
			rb->x *= -_focal / rb->z;
			rb->y *= -_focal / rb->z;
		}
	}

	if (projection != Projection_Orthographic) {
		// Restore offset
		ra->x += Uo;
		ra->y += Vo;
		rb->x += Uo;
		rb->y += Vo;
	}

	return clipLine(&(ra->x), &(ra->y), &(rb->x), &(rb->y));
} // clipLine3D

/** operator << */
ostream& operator << (ostream& s, const ViewPort& view)
{
	s << "Viewport" << endl;
	s << "\tSize   :" << view.relWidth() << " x " << view.relHeight() << endl;
	s << "\tImage  :" << view.width() << " x " << view.height() << endl;
	s << "\tZoom   : = " << view.zoom() << endl;
	s << "\tAspect : = " << view.aspect << endl;
	s << "\tOffset : = " << view.Uofs() << ", " << view.Vofs() << endl;
	return s;
} /* operator << */
