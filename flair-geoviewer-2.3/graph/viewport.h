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

#ifndef __VIEWPORT_H
#define __VIEWPORT_H

#include <ostream>

#include "os.h"
#include "conic.h"
#include "point.h"
#include "matrix4.h"

// Viewport limits
#define MINWINDOW	1.0e-8
#define MAXWINDOW	1.0e13

enum ProjectionType {
	Projection_Orthographic = 0,	// Orthographic projection
	Projection_Perspective  = 1,	// Perspective projection starting from the camera
	Projection_Combo        = 2	// Perspective starting from the viewing plane,
					// therefore 2D projections can be overlaid
};

/** ViewPort class
 *
 *  Contains the all the information necessary to convert from
 *  the Absolute to the Viewport to the Image coordinates
 *  either for orthographic or perspective projections
 *
 *  System     Coord  Description
 *  =========  =====  ==============================
 *  Absolute   x,y,z  absolute coordinate system
 *
 *  Relative   u,v,w  relative coordinate system (plane of viewport w=0.0)
 *
 *  Image      i,j    Screen/Image coordinate system only a fraction of the viewport
 *
 *  Projected  up,vp  Projected 3D system (orthographic or perspective wp=w)
 *
 *
 *  Conversions
 *  ===========
 *  [x,y,z] = matrix * [u,v,w]		relative to absolute
 *  [u,v,w] = invMatrix * [x,y,z]	absolute to relative
 *
 *
 *  Extends
 *  =======
 *  Relative   [minu .. maxv] x [minu .. maxv]
 *
 *
 *  Convections:
 *  ============
 *  Suffix   Meaning
 *  ~~~~~~   ~~~~~~~
 *  ..d      double floating precision for image coordinates
 *  ..c      center of the pixel instead of the top-left corner
 *  ..p      projected coordinates
 *
 *  Prefix   Meaning
 *  ~~~~~~   ~~~~~~~
 *  d..      transformation of vector and not of position
 *           (translation not included)
 */
class ViewPort {
public:
	double	aspect;		/** aspect ration X/Y		*/
	enum ProjectionType	projection;	/** orthographic of perspective	*/

	Conic	ulow, uhigh;	/** Conics bounding the window	*/
	Conic	vlow, vhigh;

private:
	/* Relative 2D viewport */
	double	_minu, _maxu;	/** limits in u			*/
	double	_minv, _maxv;	/** limits in v			*/

	double	chkminu;	/** limit to check against in X	*/
	double	chkmaxu;
	double	chkminv;
	double	chkmaxv;	// taking into account numerical precision

	/* image painter coordinates */
	double	Uo, Vo;		/** offset of display		*/
	double	_zoom;		/** zoom factor			*/
	double	scaleU;		/** width/(2*extU)		*/
	double	scaleV;		/** height/(2*extV)		*/
	int	_width;		/** pixel width			*/
	int	_height;
	int	width2;
	int	height2;

	/* 3d projection information */
	bool	_fixFOV;	/** keep fixed FOV or focal	*/
	bool	_fishEye;	/** FOV>180			*/
	double	_fov;		/** field of view in radians (U)*/
	double	_focal;		/** focal length -w-location	*/

	/* transformation matrix to real coordinates */
	Matrix4	_matrix;	/** transformation matrix	*/
				/*  (u,v,w) -> (x,y,z)		*/
				/*  screen to absolute		*/
	Matrix4 _invMatrix;	/** inverse matrix		*/
				/*  (x,y,z) -> (u,v,w)		*/
				/*  absolute to screen		*/

public:
	ViewPort(int w, int h);
	~ViewPort()	{}

	/* window functions */
	void	init(int w, int h);
	void	window( const double minx, const double miny,
			const double maxx, const double maxy);
	bool	invalidWindow();
	void	calcWindow(const double z=1.0);

	void	matrix(const Matrix4& m);
const	Matrix4& matrix()	const	{ return _matrix; }
const	Matrix4& invMatrix()	const	{ return _invMatrix; }
	double	matrix(const int row, const int col)	const { return _matrix(row, col); }
	double	invMatrix(const int row, const int col)	const { return _invMatrix(row, col); }
	void	transform(const Matrix4& m);

	void	origin(const double x, const double y, const double z);
	void	origin(double *x, double *y, double *z) const;
	Point	origin()	const	{ return Point(uv2x(Uo,Vo), uv2y(Uo,Vo), uv2z(Uo,Vo)); }

	void	moveOriginTo0();

	bool	inside(double u, double v) const;
	bool	insideU(double u) const	{ return (chkminu<=u && u<=chkmaxu); }
	bool	insideV(double v) const	{ return (chkminv<=v && v<=chkmaxv); }
	bool	inside(const Conic& conic) const;

	double	relWidth()	const	{ return _maxu - _minu; }
	double	relHeight()	const	{ return _maxv - _minv; }
	double	minu()		const	{ return _minu; }
	double	minv()		const	{ return _minv; }
	double	maxu()		const	{ return _maxu; }
	double	maxv()		const	{ return _maxv; }

	double	Uofs()		const	{ return Uo; }
	double	Vofs()		const	{ return Vo; }

	/* width/height in pixels of image */
	int	width()		const	{ return _width; }
	int	height()	const	{ return _height; }

	/* width/height in real coordinates of image */
	double	imageWidth()	const	{ return (double)_width  / scaleU; }
	double	imageHeight()	const	{ return (double)_height / scaleV; }

	void	zoom(const double z);
	double	zoom()		const	{ return _zoom; }

	void	offset(const double u=0.0, const double v=0.0)	 {Uo=u; Vo=v;}
	double	Sx()		const	{ return scaleU; }
	double	Sy()		const	{ return scaleV; }
	double	Dx()		const	{ return 1.0/scaleU; }
	double	Dy()		const	{ return 1.0/scaleV; }

	/* --- conversion from image to relative coordinates --- */

	/* --- Relative to Image --- */
	/** @return rounded pixel x position of image */
	int	u2i(const double x) const
				{ return width2 + Round((x-Uo)*scaleU); }
	/** @return rounded pixel y position of image */
	int	v2j(const double y) const
				{ return height2 - Round((y-Vo)*scaleV); }

	/** @return rounded center of pixel x position of image */
	int	u2ic(const double x) const
				{ return width2 + Round((x-Uo)*scaleU-0.5); }
	/** @return rounded center of pixel y position of image */
	int	v2jc(const double y) const
				{ return height2 - Round((y-Vo)*scaleV-0.5); }

	/** @return pixel x position of image as double */
	double	u2id(const double x) const
				{ return (double)width2 + (x-Uo)*scaleU; }

	/** @return pixel y position of image as double  */
	double	v2jd(const double y) const
				{ return (double)height2 - (y-Vo)*scaleV; }

	/* --- Image to Relative --- */
	/** @return image position of pixel x */
	double	i2u(const int x) const { return (double)(x-width2)/scaleU+Uo; }
	/** @return image position of pixel y */
	double	j2v(const int y) const { return (double)(height2-y)/scaleV+Vo; }

	/** @return image position of pixel x */
	double	i2u(const double x) const { return (double)(x-width2)/scaleU+Uo; }
	/** @return image position of pixel y */
	double	j2v(const double y) const { return (double)(height2-y)/scaleV+Vo; }

	/** @return image position of center of pixel x */
	double	ic2u(const int x) const { return ((double)(x-width2)+0.5)/scaleU+Uo; }
	/** @return image position of center of pixel y */
	double	jc2v(const int y) const { return ((double)(height2-y)-0.5)/scaleV+Vo; }

	/** @return image position of center of pixel x */
	double	ic2u(const double x) const { return ((x-(double)width2)+0.5)/scaleU+Uo; }
	/** @return image position of center of pixel y */
	double	jc2v(const double y) const { return (((double)height2-y)-0.5)/scaleV+Vo; }

	/* --- Image to Absolute --- */
	void	ij2xyz(const int i, const int j, double *x, double *y, double *z) const;

	/* --- Relative to absolute --- */
	double	uv2x(const double u, const double v) const
			{ return _matrix(0,0)*u + _matrix(0,1)*v + _matrix(0,3); }
	double	uv2y(const double u, const double v) const
			{ return _matrix(1,0)*u + _matrix(1,1)*v + _matrix(1,3); }
	double	uv2z(const double u, const double v) const
			{ return _matrix(2,0)*u + _matrix(2,1)*v + _matrix(2,3); }
	Point	uv2xyz(const double u, const double v) const
			{ return Point(uv2x(u,v), uv2y(u,v), uv2z(u,v)); }

	double	uvw2x(const double u, const double v, const double w) const
			{ return _matrix(0,0)*u + _matrix(0,1)*v + _matrix(0,2)*w + _matrix(0,3); }
	double	uvw2y(const double u, const double v, const double w) const
			{ return _matrix(1,0)*u + _matrix(1,1)*v + _matrix(1,2)*w + _matrix(1,3); }
	double	uvw2z(const double u, const double v, const double w) const
			{ return _matrix(2,0)*u + _matrix(2,1)*v + _matrix(2,2)*w + _matrix(2,3); }
	Point	uvw2xyz(const double u, const double v, const double w) const
			{ return Point(uvw2x(u,v,w), uvw2y(u,v,w), uvw2z(u,v,w)); }

	double	uv2dx(const double u, const double v) const
			{ return _matrix(0,0)*u + _matrix(0,1)*v; }
	double	uv2dy(const double u, const double v) const
			{ return _matrix(1,0)*u + _matrix(1,1)*v; }
	double	uv2dz(const double u, const double v) const
			{ return _matrix(2,0)*u + _matrix(2,1)*v; }
	Point	uv2dxyz(const double u, const double v) const
			{ return Point(uv2dx(u,v), uv2dy(u,v), uv2dz(u,v)); }

	double	uvw2dx(const double u, const double v, const double w) const
			{ return _matrix(0,0)*u + _matrix(0,1)*v + _matrix(0,2)*w; }
	double	uvw2dy(const double u, const double v, const double w) const
			{ return _matrix(1,0)*u + _matrix(1,1)*v + _matrix(1,2)*w; }
	double	uvw2dz(const double u, const double v, const double w) const
			{ return _matrix(2,0)*u + _matrix(2,1)*v + _matrix(2,2)*w; }
	Point	uvw2dxyz(const double u, const double v, const double w) const
			{ return Point(uvw2dx(u,v,w), uvw2dy(u,v,w), uvw2dz(u,v,w)); }

	/* --- Absolute to relative --- */
	double	xyz2u(const double x, const double y, const double z) const
			{ return _invMatrix(0,0)*x + _invMatrix(0,1)*y + _invMatrix(0,2)*z + _invMatrix(0,3); }
	double	xyz2v(const double x, const double y, const double z) const
			{ return _invMatrix(1,0)*x + _invMatrix(1,1)*y + _invMatrix(1,2)*z + _invMatrix(1,3); }
	double	xyz2w(const double x, const double y, const double z) const
			{ return _invMatrix(2,0)*x + _invMatrix(2,1)*y + _invMatrix(2,2)*z + _invMatrix(2,3); }
	void	xyz2uv(const double x, const double y, const double z, double *u, double *v) const {
				*u = xyz2u(x,y,z);
				*v = xyz2v(x,y,z);
			}
	void	xyz2uvw(const double x, const double y, const double z, double *u, double *v, double *w) const {
				*u = xyz2u(x,y,z);
				*v = xyz2v(x,y,z);
				*w = xyz2w(x,y,z);
			}

	double	xyz2u(const Point& p) const { return xyz2u(p.x, p.y, p.z); }
	double	xyz2v(const Point& p) const { return xyz2v(p.x, p.y, p.z); }
	double	xyz2w(const Point& p) const { return xyz2w(p.x, p.y, p.z); }
	void	xyz2uv(const Point& p, double *u, double *v) const {
				*u = xyz2u(p);
				*v = xyz2v(p);
			}
	void	xyz2uvw(const Point& p, double *u, double *v, double *w) const {
				*u = xyz2u(p);
				*v = xyz2v(p);
				*w = xyz2w(p);
			}
	void	xyz2uvw(const Point& p, Point* r) const {
				r->x = xyz2u(p);
				r->y = xyz2v(p);
				r->z = xyz2w(p);
			}
	Point	xyz2uvw(const Point& p) const { return Point(xyz2u(p), xyz2v(p), xyz2w(p));}
	Point	xyz2uvw(const double x, const double y, const double z) const
			{ return Point(xyz2u(x,y,z), xyz2v(x,y,z), xyz2w(x,y,z)); }

	// absolute to relative for movements
	void	dxyz2duvw(const double dx, const double dy, const double dz,
				double *du, double *dv, double *dw) const {
				*du = _invMatrix(0,0)*dx + _invMatrix(0,1)*dy + _invMatrix(0,2)*dz;
				*dv = _invMatrix(1,0)*dx + _invMatrix(1,1)*dy + _invMatrix(1,2)*dz;
				*dw = _invMatrix(2,0)*dx + _invMatrix(2,1)*dy + _invMatrix(2,2)*dz;
			}

	// Relative to Projected (3D)
	/** @return fov half angle of camera */
	double	fov()			const { return projection==Projection_Orthographic? 0.0 :_fov; }
	void	fov(double a);
	/** @return focal length distance of camera from viewing plane */
	double	focalLength()		const { return projection==Projection_Orthographic? 0.0 : _focal; }
	void	focalLength(double f);

	/** @return camera position in absolute space */
	void	cameraPosition(double *x, double *y, double *z) const {
				// there is no meaning in Orthographic projection
				*x = uv2x(Uo,Vo) + uvw2dx(0.,0.,_focal),
				*y = uv2y(Uo,Vo) + uvw2dy(0.,0.,_focal),
				*z = uv2z(Uo,Vo) + uvw2dz(0.,0.,_focal);
			}
	Point	cameraPosition() const {
				Point p;
				cameraPosition(&p.x, &p.y, &p.z);
				return p;
			}

	/** @return x,y,z as starting position for a ray */
	void	rayPosition(double u, double v, double *x, double *y, double *z) const {
				if (projection==Projection_Perspective) {
					*x = uv2x(Uo,Vo) + uvw2dx(0.,0.,_focal),
					*y = uv2y(Uo,Vo) + uvw2dy(0.,0.,_focal),
					*z = uv2z(Uo,Vo) + uvw2dz(0.,0.,_focal);
				} else {
					*x = uv2x(u,v);
					*y = uv2y(u,v);
					*z = uv2z(u,v);
				}
			}

	// rename to uv2dxyz3D?
	/** @return normalized perspective direction taking into account the FOV */
	void	rayDirection(double u, double v, double *dx, double *dy, double *dz) const {
				if (projection == Projection_Orthographic) {
					*dx = -matrix(0,2);
					*dy = -matrix(1,2);
					*dz = -matrix(2,2);
#if 0
				} else
				if (_fishEye) {
					// Projection on a sphere
					assert(_focal == 0.0);
					double fa = (double)width()/_fov;
					double fu = PI - (u2id(u)-width2)/fa;
					double fv = (height2-v2jd(v))/fa;
					double du = cos(fu)*sin(fv);
					double dv = sin(fu)*sin(fv);
					double dw = cos(fv);
					*dx = uvw2dx(du,dv,dw);
					*dy = uvw2dy(du,dv,dw);
					*dz = uvw2dz(du,dv,dw);
#endif
				} else {
					// Projection on a plane
					u -= Uo;
					v -= Vo;
					double d = 1.0 / sqrt(u*u + v*v + _focal*_focal);
					double du =      u * d;
					double dv =      v * d;
					double dw =-_focal * d;
					*dx = uvw2dx(du,dv,dw);
					*dy = uvw2dy(du,dv,dw);
					*dz = uvw2dz(du,dv,dw);
				}
			}

	/** convert absolute coordinates to relative projected in 3D */
	void	xyz2uvw3D(const Point& p, Point* r) const {
				xyz2uvw(p, r);
				uvw2uv3D(r);
			}

	void	xyz2ij3D(const Point& p, int *i, int *j) {
				Point r;
				xyz2uvw3D(p, &r);
				*i = u2i(r.x);
				*j = v2j(r.y);
			}

	// Axes
	Vector	axisu() const { return Vector(_matrix(0,0), _matrix(1,0), _matrix(2,0)); }
	Vector	axisv() const { return Vector(_matrix(0,1), _matrix(1,1), _matrix(2,1)); }
	Vector	axisw() const { return Vector(_matrix(0,2), _matrix(1,2), _matrix(2,2)); }

	/* clipping */
	bool	clipLine(double *x1, double *y1, double *x2, double *y2) const;
	bool	clipLine3D(const Point& a, const Point& b, Point* va, Point* vb) const;

private:
	void	calculateFOV();
	void	calculateFocalLength();
	void	computeMatrices();

	/** convert relative coordinates to relative projected in 3D
	 * @param u,v,w		relative coordinates
	 * @param up		return projected u
	 * @param vp		return projected v
	 */
	void	uvw2uv3D(Point *r) const {
				if (projection == Projection_Orthographic)
					return;
				else {
					// The viewport is centered on the projection plane
					// ... move to viewer
					double w = r->z - _focal;
					if (Eq0(w,Vector::eps())) {
						r->x = r->y = 1e10;	// FIXME
						return;
					}
					w = -_focal / w;
					r->x *= w;
					r->y *= w;
				}
			}

}; // ViewPort

std::ostream& operator << (std::ostream&, const ViewPort&);

#endif
