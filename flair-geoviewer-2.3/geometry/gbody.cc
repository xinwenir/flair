/*
 *
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
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <ostream>

#include "line.h"
#include "bmath.h"
#include "conic.h"
#include "gbody.h"
#include "obbox.h"
#include "matrix3.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;

#define N_CYLINDER	16
#define N_SPHERE_LON	16
#define N_SPHERE_LAT	 8
#define N_TORUS_LON	16
#define N_TORUS_LAT	32

#define ROUND2INT(x,a)		Eq0((double)(int)x-x, a)? (double)(int)x : x

GOPRBody GBody::tnull    ("$", NULLbody);
GOPRBody GBody::tplus    ("+", PLUSbody);
GOPRBody GBody::tminus   ("-", MINUSbody);
GOPRBody GBody::tunion   ("|", UNIONbody);
GOPRBody GBody::tuniverse("@", UNIVERSEbody);
GOPRBody GBody::tleft    ("(", LEFTbody);
GOPRBody GBody::tright   (")", RIGHTbody);

const char *GBody::_typeStr[] = {
	"YZP", "XZP",
	"XYP", "PLA",

	"RPP", "BOX",
	"WED", "RAW",
	"TET", //For TET, added by zxw

	"SPH", "ELL",

	"RCC", "REC", "TRC",

	"XCC", "YCC", "ZCC",
	"XEC", "YEC", "ZEC",

	"ARB", "QUA",

	"null",
	"+","-","|",
	"@","(",")",

	"ERR"};

// locationWrt matrix
//    X = completed
//    - = to be done (if possible)
// row is the first item checked against the column. Always the smallest body index
// check against the bigger ones
//
//	 bbx  YZP XZP XYP PLA   RPP BOX WED RAW   SPH ELL   RCC REC TRC   XCC YCC ZCC XEC YEC ZEC   ARB QUA
//  bbx   X    -   -   -   -     -   -   -   -     -   -     -   -   -     -   -   -   -   -   -     -   -
//
//  YZP   X    X   X   X   X     X   X   X   X     X   -     -   -   -     X   X   X   -   -   -     X   -
//  XZP   X        X   X   X     X   X   X   X     X   -     -   -   -     X   X   X   -   -   -     X   -
//  XYP   X            X   X     X   X   X   X     X   -     -   -   -     X   X   X   -   -   -     X   -
//  PLA   X                X     X   X   X   X     X   -     -   -   -     X   X   X   -   -   -     X   -
//
//  RPP   -                      -   -   -   -     -   -     -   -   -     -   -   -   -   -   -     -   -
//  BOX   -                          -   -   -     -   -     -   -   -     -   -   -   -   -   -     -   -
//  WED   -                              -   -     -   -     -   -   -     -   -   -   -   -   -     -   -
//  RAW   -                                  -     -   -     -   -   -     -   -   -   -   -   -     -   -
//
//  SPH   -                                        X   -     -   -   -     X   X   X   -   -   -     -   -
//  ELL   -                                            -     -   -   -     -   -   -   -   -   -     -   -
//
//  RCC   -                                                  -   -   -     -   -   -   -   -   -     -   -
//  REC   -                                                      -   -     -   -   -   -   -   -     -   -
//  TRC   -                                                          -     -   -   -   -   -   -     -   -
//
//  XCC   -                                                                X   X   X   -   -   -     -   -
//  YCC   -                                                                    X   X   -   -   -     -   -
//  ZCC   -                                                                        X   -   -   -     -   -
//  XEC   -                                                                            -   -   -     -   -
//  YEC   -                                                                                -   -     -   -
//  ZEC   -                                                                                    -     -   -
//
//  ARB   -                                                                                          -   -
//  QUA   -                                                                                              -
//
// Reference for various intersections: http://www.realtimerendering.com/intersections.html

/** bbox_addRotatedEllipse
 * add to the bounding box an arbitrary rotated ellipse
 */
static void bbox_addRotatedEllipse(BBox& bb,
		const Vector& rotX, const Vector& rotY,
		const Point& pos, double ra, double rb)
{
	for (int i = 0; i < 3; i++) {
		double a = ra * rotX[i];
		double b = rb * rotY[i];
		double a2b2 = a*a + b*b;

		if (a2b2 < SMALL) continue;

		double sint = sqrt(b*b/a2b2);
		double cost = 0.0;
		if (Abs(sint) > SMALL)
			cost = a / b * sint;
		else {
			cost = 1.0;
			sint = 0.0;
		}

		Point e;
		e.x = rotX.x * ra * cost + rotY.x * rb * sint;
		e.y = rotX.y * ra * cost + rotY.y * rb * sint;
		e.z = rotX.z * ra * cost + rotY.z * rb * sint;

		bb.add(pos.x+e.x, pos.y+e.y, pos.z+e.z);
		bb.add(pos.x-e.x, pos.y-e.y, pos.z-e.z);
	}
} // bbox_addRotatedEllipse

double GBody::infinite = 100000.0;

static int zoneCompare(GZone* const& a, GZone* const& b)
{
	if ((size_t)a < (size_t)b)
		return -1;
	else
	if ((size_t)a > (size_t)b)
		return  1;
	else
		return  0;
}; // zoneCompare

static inline void roundFloat(double& f)
{
	if (Eq0(f,SMALL3)) f = 0.0;
} // roundFloat

static inline void roundVector(Vector& v)
{
	if (Eq0(v.x,SMALL3)) v.x = 0.0;
	if (Eq0(v.y,SMALL3)) v.y = 0.0;
	if (Eq0(v.z,SMALL3)) v.z = 0.0;
} // roundVector

/* ============================== GBody ================================ */
/** Geometrical Body */
GBody::GBody(const char *aname, const BodyType atype) :
	_id(-1),
	_type(atype),
	_color(0xFF00FF),
	_width(0),
	_generation(0),
	_nQ(0),
	_hasMatrix(false),
	show(0),
	X(Vector::Xo),
	Y(Vector::Yo),
	Z(Vector::Zo),
	xlen(0.0),
	ylen(0.0),
	zlen(0.0),
	_userBboxFlag(false)
{
	name(aname);
	_cached_obbox[0] = NULL;
	_cached_obbox[1] = NULL;
	_matrix.identity();
	_invMatrix.identity();
	zones.compare(zoneCompare);
	save();
} // GBody

/** ~GBody */
GBody::~GBody()
{
	clearOBB();
} /* ~GBody */

/** clearOBB */
void GBody::clearOBB() {
	if (_cached_obbox[0]) {
		delete _cached_obbox[0];
		_cached_obbox[0] = NULL;
	}
	if (_cached_obbox[1]) {
		delete _cached_obbox[1];
		_cached_obbox[1] = NULL;
	}
} // clearOBB

/** name - set name of body */
void GBody::name(const char *aname)
{
	strncpy(_name, aname, sizeof(_name));
	_name[sizeof(_name)-1] = 0;
} // name

/** @return memory used by body */
size_t GBody::memory() const
{
	return  sizeof(*this)
	      + mesh.memory()
	      + zones.memory();
} // memory

/** save */
void GBody::save()
{
	SP  = P;	// save variables
	SPo = Po;
	SX  = X;
	SY  = Y;
	SZ  = Z;
	sxlen = xlen;
	sylen = ylen;
	szlen = zlen;
	sshow = show;
} // save

/** restore */
void GBody::restore()
{
	P  = SP;
	Po = SPo;
	X  = SX;
	Y  = SY;
	Z  = SZ;
	xlen = sxlen;
	ylen = sylen;
	zlen = szlen;
	show = sshow;
} // restore

/** set position
 * @param x,y,z	new position to set
 */
void GBody::position(const Point& r)
{
	P = _hasMatrix? _invMatrix*r : r;
	Po = P + zlen*Z;
} // position

/** @return transformed position */
Point GBody::position() const
{
	return _hasMatrix? _matrix*P : P;
} // position

/** @return transformed X vector */
Vector GBody::vectorX() const
{
//	return _hasMatrix? _matrix.multVector(X) : X;
	return _hasMatrix? _matrix*X : X;
} // vectorX

/** @return transformed Y vector */
Vector GBody::vectorY() const
{
//	return _hasMatrix? _matrix.multVector(Y) : Y;
	return _hasMatrix? _matrix*Y : Y;
} // vectorY

/** @return transformed Z vector */
Vector GBody::vectorZ() const
{
//	return _hasMatrix? _matrix.multVector(Z) : Z;
	return _hasMatrix? _matrix*Z : Z;
} // vectorZ

/** @return transformed position */
Point GBody::node(const int n) const
{
	Point r;
	switch (n) {
		case 0: r = P;
			break;
		case 1: r = P+zlen*Z;
			break;
		case 2: r = P+xlen*X;
			break;
		case 3: r = P+ylen*Y;
			break;
	}
	return _hasMatrix? _matrix*r : r;
} // node

/** @return saved position of body */
Point GBody::savedPosition() const
{
	return _hasMatrix? _matrix*SP : SP;
} // savedPosition

/**  @return closest point, edge, face
 * Values:
 *	0	position of body
 *	1-6	closest quad of body
 *	10-19	closest handler of body
 *	20-	rotation handler of body
 */
int GBody::closest(const Point& r, const double d, const Vector& w) const
{
	double d2 = d*d;

	Point tr = _hasMatrix? _invMatrix*r : r;
	if ((P-tr).length2() <= d2) return 0;

	// Check possible nodes
	if (nodes()>1 && zlen>0.0 && (P+zlen*Z-tr).length2() <= d2) return 10;
	if (nodes()>2 && xlen>0.0 && (P+xlen*X-tr).length2() <= d2) return 11;
	if (nodes()>3 && ylen>0.0 && (P+ylen*Y-tr).length2() <= d2) return 12;

	int	close = -1;
	double	dmin = INFINITE;
	for (int i=0; i<nQ(); i++) {
		double dq = Q(i)(r);
		Vector g = Q(i).grad(r);
		double gg = g.normalize();
		dq = Abs(dq/gg);
		// Check if smaller but also normal is not aligned with the w axis
		if (dq<dmin && !Eq0(Abs(g*w)-1.0,0.000001)) {
			close = i;
			dmin = dq;
		}
	}

	return close+1;
} // closest

/** fixWhat
 * convert what to zero all small values
 */
void GBody::fixWhat(double *what, int n, const double eps)
{
	for (int i=0; i<n; i++)
		if (Eq0(what[i],eps))
			what[i] = 0.0;
} // fixWhat

/** matrix - set transformation matrix */
void GBody::matrix(const Matrix4& M)
{
	_hasMatrix = true;
	_matrix.copy(M);
	_matrix.fix();
	_invMatrix.inverse(M);		// Invert matrix
	_invMatrix.fix();
} // matrix

/** transform - transform quadrics with matrix */
void GBody::transform()
{
	for (int i=0; i<_nQ; i++) {
#if defined(_DUMP) && _DEBUG>1
		if (show) {
			cout << "*-* GBody::transform " << name() << ":" << i << endl;
			cout << "*-*\t Q[" << i << "] = " << Q(i) << endl;
			cout << "*-*\t M=" << endl << _invMatrix << endl;
		}
#endif
		_Q[i].transform(_invMatrix);
		_Q[i].normalize();
		DUMP(if (show) cout << "*-*\ttQ["<< i << "] = " << Q(i) << endl);
	}

	// FIXME
	// for planes and infinite bodies we should have a special transform
	// or ask to re-create the mesh!

	mesh.transform(_matrix);
} // transform

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GBody::move(int item, const Point& r, const Vector&)
{
	double len, angle;

	Point tr = _hasMatrix? _invMatrix*r : r;
	tr -= SP;

	switch (item) {
		case 10:	// Move Z freely
			len = tr.normalize();
			if (Eq0(len,SMALL3)) return;
			angle = acos(SZ * tr);
			if (!Eq0(angle,SMALL3) && !Eq0(angle-PI,SMALL3)) {
				Vector axis = SZ ^ tr;
				axis.normalize();
				Matrix4	rot;
				rot.rotate(angle, axis);
				//rot.fix();
				X = rot * SX;
				Y = rot * SY;
				Z = rot * SZ;
			} else {
				X = SX;
				Y = SY;
				Z = Eq0(angle,SMALL3)?SZ : -SZ;
			}
			zlen = len;
			Po = P + zlen*Z;
			break;

		case 11:	// Move X freely
			len = tr.normalize();
			if (Eq0(len,SMALL3)) return;
			angle = acos(SX * tr);
			if (!Eq0(angle,SMALL3) && !Eq0(angle-PI,SMALL3)) {
				Vector axis = SX ^ tr;
				axis.normalize();
				Matrix4	rot;
				rot.rotate(angle, axis);
				//rot.fix();
				X = rot * SX;
				Y = rot * SY;
				Z = rot * SZ;
			} else {
				X = Eq0(angle,SMALL3)?SX : -SX;
				Y = SY;
				Z = SZ;
			}
			xlen = len;
			Po   = P + zlen*Z;
			break;

		case 12:	// Move Y freely
			len = tr.normalize();
			if (Eq0(len,SMALL3)) return;
			angle = acos(SY * tr);
			if (!Eq0(angle,SMALL3) && !Eq0(angle-PI,SMALL3)) {
				Vector axis = SY ^ tr;
				axis.normalize();
				Matrix4	rot;
				rot.rotate(angle, axis);
				//rot.fix();
				X = rot * SX;
				Y = rot * SY;
				Z = rot * SZ;
			} else {
				X = SX;
				Y = Eq0(angle,SMALL3)?SY : -SY;
				Z = SZ;
			}
			ylen = len;
			Po   = P + zlen*Z;
			break;

		case 13:	// Move Z along SZ axis
			zlen = tr * SZ;
			if (zlen<0.0) {
				zlen = -zlen;
				Z = -SZ;
			} else
				Z = SZ;
			Po = P + zlen*Z;
			break;

		case 14:	// Move X along SX axis
			xlen = tr * SX;
			if (xlen<0.0) {
				xlen = -xlen;
				X = -SX;
			} else
				X = SX;
			Po = P + zlen*Z;
			break;

		case 15:	// Move Y along SY axis
			ylen = tr * SY;
			if (ylen<0.0) {
				ylen = -ylen;
				Y = -SY;
			} else
				Y = SY;
			Po = P + zlen*Z;
			break;

		default:
			position(r);
	}
	roundVector(P);
	roundVector(Po);
	roundVector(X);
	roundVector(Y);
	roundVector(Z);
	roundFloat(xlen);
	roundFloat(ylen);
	roundFloat(zlen);
} // move

/** rotate body around an axis and current body position
 * @param angle	to rotate
 * @param axis	to rotate
 */
void GBody::rotate(const double angle, const Vector& axis)
{
	Matrix4	rot;
	rot.rotate(angle,axis);

//	X = rot.multVector(X);
//	Y = rot.multVector(Y);
//	Z = rot.multVector(Z);
	X = rot * X;
	Y = rot * Y;
	Z = rot * Z;

	Po = P + zlen*Z;

	roundVector(Po);
	roundVector(X);
	roundVector(Y);
	roundVector(Z);
} // rotate

/** bbox */
void GBody::bbox(const BBox& bb)
{
	_userBboxFlag = true;
	_userBbox = bb;
} // bbox

/** @return bounding box of body */
BBox GBody::_bbox() const
{
	BBox bb;
	if (mesh.nvertices()) {
		// use the vertices as a guess
		bb.add(position());
		bb.Union(mesh.bbox());
	} else
		bb.infinite();

	return bb;
} // _bbox

/** obbox */
OBBox* GBody::obbox(bool in) const
{
	if (!in) {
		int idx = 0; // out obb
		if (!_cached_obbox[idx])
			_cached_obbox[idx] = updateOBB(in);
		return _cached_obbox[idx];
	}
	int idx = 1; // in obb
	if (!_cached_obbox[idx])
		_cached_obbox[idx] = updateOBB(in);
	return _cached_obbox[idx];
} // obbox

/** @return oriented bounding box of body */
OBBox* GBody::updateOBB(bool in) const
{
	OBBox* bb = new OBBox();

	if (!in) {
		if (mesh.nvertices()) {
			// use the vertices as a guess
			bb->add(position());
			for (int i=0; i<mesh.nvertices(); i++)
				bb->add(mesh.vertex(i));
		} else
			bb->infinite();
	} else {
		// TODO: Implement inside obb?
	}

	return bb;
} // obbox

/**
 * @param body	other body to check location against
 * @return location of this body with respect to another body
 *
 * WARNING: Internally it calls the _locationWrt of the body
 * with lower type() to the one with higher
 */
Location GBody::locationWrt(const GBody *body) const
{
	if (type() > body->type()) {
		// Invert the location
		Location loc = body->_locationWrt(this);
		if (loc == LOCATION_AinB)
			return LOCATION_BinA;
		else
		if (loc == LOCATION_BinA)
			return LOCATION_AinB;
		return loc;
	} else
		return _locationWrt(body);
} // locationWrt

/** check location of this bounding box against the other body's bounding box
 * @param body	other body to check location against
 * @return location of this body bbox with respect to another body bbox
 */
Location GBody::bbLocationWrt(const GBody* body) const
{
	return GBody::_locationWrt(body);
} // bbLocationWrt

/**
 * @param body	other body to check location against
 * @return location of this body with respect to another body
 */
Location GBody::_locationWrt(const GBody *body) const
{
	BBox Abb = bbox();
	BBox Bbb = body->bbox();
	if (Abb.overlap(Bbb))
		return LOCATION_OVERLAP;
	return LOCATION_OUTSIDE;
} // _locationWrt

/**
 * Add a quad: Cxx*x^2 + Cyy*y^2 + Czz*z^2 +
 *             Cxy*x*y + Cxz*x*z + Cyz*y*z +
 *              Cx*x   +  Cy*y   +  Cz*z   + C = 0
 *
 *          0       1      2      3
 *     0 / Cxx     Cxy/2  Cxz/2  Cx/2 \
 *     1 | Cxy/2   Cyy    Cyz/2  Cy/2 |
 *     2 | Cxz/2   Cyz/2  Czz    Cz/2 |
 *     3 \ Cx/2    Cy/2   Cz/2   C    /
 */
void GBody::addQuad(double Cxx, double Cyy, double Czz,
		    double Cxy, double Cxz, double Cyz,
		    double Cx,  double Cy,  double Cz,
		    double CC)
{
	assert(_nQ<BODYQUADS);
	_Q[_nQ].set(Cxx,Cyy,Czz,  Cxy,Cxz,Cyz, Cx,Cy,Cz, CC);
	_Q[_nQ].normalize();
	_nQ++;
} // body_addQuad

void GBody::addQuad(double Cx,  double Cy,  double Cz, double CC)
{
	assert(_nQ<BODYQUADS);
	_Q[_nQ].set(Cx,Cy,Cz, CC);
	_Q[_nQ].normalize();
	_nQ++;
} // body_addQuad

/** make cone quadrics from 3 vectors and 2 positions */
void GBody::makeConeQuads()
{
	Matrix4 M, Minv, T;

	Minv.make(X,Y,Z);
	Minv.inverse();
	T.translate(P);
	M.multiply(T, Minv);
	Minv.inverse(M);

	// transform first quad
	_Q[0].transform(Minv);
	_Q[0].normalize();

	// two planes for the caps
	//        Cx     Cy     Cz      C
	addQuad(-Z.x,-Z.y,-Z.z, Z*P);	// base
	addQuad( Z.x, Z.y, Z.z,-Z*Po);	// apex
} // makeConeQuads

/** check for orthogonality in X Y Z */
void GBody::checkOrthogonal(char *err)
{
	// check for orthogonality
	if (err==NULL) return;
	if (!Eq0(X*Y,SMALL2)) sprintf(err,"X-Y vectors not perpendicular X.Y=%g",X*Y);
	if (!Eq0(X*Z,SMALL2)) sprintf(err,"X-Z vectors not perpendicular X.Z=%g",X*Z);
	if (!Eq0(Y*Z,SMALL2)) sprintf(err,"Y-Z vectors not perpendicular Y.Z=%g",Y*Z);
} // checkOrthogonal

/** createEllConeMesh - create mesh of an elliptical cone
 * @param Rbx	base radius along X vector
 * @param Rby	base radius along Y vector
 * @param Rax	apex radius along X vector
 * @param Ray	apex radius along Y vector
 */
void GBody::createEllConeMesh(	const double Rbx, const double Rby,
				const double Rax, const double Ray,
				const bool isinfinite)
{
	bool first = mesh.isEmpty();
	const double step = PI2 / (double)N_CYLINDER;
	double ang = 0.0;

	Point A, B, C;
	if (isinfinite) {
		B = P - infinite*Z;
		A = P + infinite*Z;
		C = P;		// additional ring!
	} else {
		B = P;
		A = Po;
	}

	if (first) mesh.allocateVertices(2+2*N_CYLINDER);

	mesh.vertex(0) = B;	// base
	mesh.vertex(1) = A;	// apex

	// Build base
	//if (!Eq0(Rbx,SMALL) && !Eq0(Rby,SMALL))	// ignore Cones with tip
	//if (!Eq0(Rax,SMALL) && !Eq0(Ray,SMALL))	// ignore Cones with tip

	const int base = 2;
	const int apex = 2 + N_CYLINDER;

	for (int i=0; i<N_CYLINDER; i++, ang+=step) {
		double c, s;
		bsincos(ang, &s, &c);
		mesh.vertex(i+base) = B + (Rbx*c)*X + (Rby*s)*Y;
		mesh.vertex(i+apex) = A + (Rax*c)*X + (Ray*s)*Y;
	}
	mesh.calcBbox();

	if (first) {
		// XXX I could do in a single loop
		// but this way I could potentially use the fan opengl commands

		// create base faces
		for (int i=0; i<N_CYLINDER; i++) {
			int in = (i==N_CYLINDER-1? 0 : i+1); // next i
			mesh.add(0, base+in, base+i, false, true, false);
		//}

		// create apex faces
		//for (int i=0; i<N_CYLINDER; i++) {
		//	int in = (i==N_CYLINDER-1? 0 : i+1); // next i
			mesh.add(1, apex+i, apex+in, false, true, false);
		//}

		// create side faces
		//for (int i=0; i<N_CYLINDER; i++) {
		//	int in = (i==N_CYLINDER-1? 0 : i+1); // next i
			mesh.add(base+i, apex+in, apex+i,  false, true, true);
			mesh.add(base+i, base+in, apex+in, true,  true, false);
		}

		mesh.process();
		assert(mesh.isClosed());
		assert(mesh.isOrientable());
#if _DEBUG>2
		cout << endl;
		cout << "CYL Mesh ";
		cout << " isClosed=" << mesh.isClosed();
		cout << " isOrientable=" << mesh.isOrientable() << endl;
		cout << "CYL volume=" << mesh.volume() << endl;
#endif
	}
} // createEllConeMesh

/** createEllipsoidMesh - create an ellipsoidal mesh for SPH and ELL */
void GBody::createEllipsoidMesh()
{
	bool first = mesh.isEmpty();

	if (first) mesh.allocateVertices(N_SPHERE_LON*(N_SPHERE_LAT-1)+2);

//	vertices(N_SPHERE_LON*(N_SPHERE_LAT-1)+2);
//	edges(N_SPHERE_LON * 2*(N_SPHERE_LAT-1));

	// Create Mesh
	const double stepf = PI2 / (double)N_SPHERE_LON;
	const double step0 = PI  / (double)N_SPHERE_LAT;

	int k = 0;
	mesh.vertex(k++) = P + zlen*Z;
	mesh.vertex(k++) = P - zlen*Z;

	double ang0 = step0;
	for (int j=0; j<N_SPHERE_LAT-1; j++, ang0 += step0) {
		double c0, s0;
		bsincos(ang0, &s0, &c0);
		double ang = 0.0;
		for (int i=0; i<N_SPHERE_LON; i++, ang += stepf) {
			double c, s;
			bsincos(ang, &s, &c);
			mesh.vertex(k++) = P +
					(xlen*c*s0) * X +
					(ylen*s*s0) * Y +
					(zlen*c0)   * Z;
			//mesh.vertex(k++).set(what[0] + R*c*s0,
			//		what[1] + R*s*s0,
			//		what[2] + R*c0);
		}
	}
	mesh.calcBbox();

	if (first) {
		// Create faces
		const int base = 2;
		const int apex = 2 + N_SPHERE_LON*(N_SPHERE_LAT-2);
		for (int i=0; i<N_SPHERE_LON; i++) {
			int in = (i==N_SPHERE_LON-1? 0 : i+1); // next i
			mesh.add(0, base+i,  base+in);
			mesh.add(1, apex+in, apex+i);
		}

		int start = base;
		for (int j=0; j<N_SPHERE_LAT-2; j++, start += N_SPHERE_LON)
			for (int i=0; i<N_SPHERE_LON; i++) {
				int in = (i==N_SPHERE_LON-1? 0 : i+1); // next i
				mesh.add(start+i, start+i+N_SPHERE_LON, start+in, true, false, true);
				mesh.add(start+in, start+i+N_SPHERE_LON, start+in+N_SPHERE_LON, false, true, true);
			}

		mesh.process();
		assert(mesh.isClosed());
		assert(mesh.isOrientable());
#if _DEBUG>2
		cout << endl;
		cout << "SPH Mesh ";
		cout << " isClosed=" << mesh.isClosed();
		cout << " isOrientable=" << mesh.isOrientable() << endl;
		cout << "SPH volume=" << mesh.volume() << endl;
#endif
	}
} // createEllipsoidMesh

/** findXYZ - find an arbitrary XYZ orthogonal system based only
 * on Z information
 */
void GBody::findXYZ(char *err)
{
	if (Eq0(Z.length2(),SMALL)) {
		if (err) strcpy(err,"Invalid height |Z|=0");
		X = Vector::Xo;
		Y = Vector::Yo;
		Z = Vector::Zo;
	} else {
		Y = Z.orthogonal();
		// create the rotation matrix
		Y.normalize();
		X = Y ^ Z;
		X.normalize();
	}
} // findXYZ

/** newBody
 * @param aname		name of body
 * @param atype		string type of body
 * @return a new allocated GBody()
 */
GBody* GBody::newBody(const char *name, const char *type)
{
	try {
		     if (!strcmp(type, "SPH")) return new GSPHBody(name);
		else if (!strcmp(type, "RPP")) return new GBOXBody(name, RPPbody);
		else if (!strcmp(type, "BOX")) return new GBOXBody(name, BOXbody);
		else if (!strcmp(type, "RCC")) return new GRCCBody(name);
		else if (!strcmp(type, "REC")) return new GRECBody(name);
		else if (!strcmp(type, "TRC")) return new GTRCBody(name);
		else if (!strcmp(type, "ELL")) return new GELLBody(name);
		else if (!strcmp(type, "ARB")) return new GARBBody(name);
		else if (!strcmp(type, "WED") ||
			 !strcmp(type, "RAW")) return new GWEDBody(name);
		else if (!strcmp(type, "TET")) return new GTETBody(name); //For TET,added by zxw
		else if (!strcmp(type, "XYP")) return new GPLABody(name, XYPbody);
		else if (!strcmp(type, "XZP")) return new GPLABody(name, XZPbody);
		else if (!strcmp(type, "YZP")) return new GPLABody(name, YZPbody);
		else if (!strcmp(type, "PLA")) return new GPLABody(name, PLAbody);
		else if (!strcmp(type, "XCC")) return new GInfEllCylBody(name, XCCbody);
		else if (!strcmp(type, "YCC")) return new GInfEllCylBody(name, YCCbody);
		else if (!strcmp(type, "ZCC")) return new GInfEllCylBody(name, ZCCbody);
		else if (!strcmp(type, "XEC")) return new GInfEllCylBody(name, XECbody);
		else if (!strcmp(type, "YEC")) return new GInfEllCylBody(name, YECbody);
		else if (!strcmp(type, "ZEC")) return new GInfEllCylBody(name, ZECbody);
		else if (!strcmp(type, "QUA")) return new GQUABody(name);

		else if (!strcmp(type, "TRX")) return new GTorusBody(name, TRXbody);
		else if (!strcmp(type, "TRY")) return new GTorusBody(name, TRYbody);
		else if (!strcmp(type, "TRZ")) return new GTorusBody(name, TRZbody);
		else
			return new GERRBody(name,ERRbody);
	} catch ( ... ) {
		return new GERRBody(name,ERRbody);
	}
} // newBody

/** newBody
 * @param aname		name of body
 * @param atype		body type
 * @return a new allocated GBody()
 */
GBody* GBody::newBody(const char *name, const BodyType type)
{
	try {
		switch (type) {
			case SPHbody: return new GSPHBody(name);
			case RPPbody: return new GBOXBody(name, type);
			case BOXbody: return new GBOXBody(name, type);
			case RCCbody: return new GRCCBody(name);
			case RECbody: return new GRECBody(name);
			case TRCbody: return new GTRCBody(name);
			case ELLbody: return new GELLBody(name);
			case ARBbody: return new GARBBody(name);
			case WEDbody:
			case RAWbody: return new GWEDBody(name);
			case TETbody: return new GTETBody(name); //For TET, added by zxw
			case XYPbody: return new GPLABody(name, type);
			case XZPbody: return new GPLABody(name, type);
			case YZPbody: return new GPLABody(name, type);
			case PLAbody: return new GPLABody(name, type);
			case XCCbody: return new GInfEllCylBody(name, type);
			case YCCbody: return new GInfEllCylBody(name, type);
			case ZCCbody: return new GInfEllCylBody(name, type);
			case XECbody: return new GInfEllCylBody(name, type);
			case YECbody: return new GInfEllCylBody(name, type);
			case ZECbody: return new GInfEllCylBody(name, type);
			case QUAbody: return new GQUABody(name);
			default:
				return new GERRBody(name,ERRbody);
		}
	} catch ( ... ) {
		return new GERRBody(name,ERRbody);
	}
} // newBody

/**
 * Find if a point (x,y,z) is inside the body or not, using the
 * quadratic surface equations and the direction the point is moving
 * Returns: true  inside
 *          false outside
 * XXX I've tried with const Point&p and d but in the end it was slower
 *     than with the plain references x,y,z...
 */
bool GBody::inside(const double  x, const double  y, const double  z,
		   const double dx, const double dy, const double dz,
		   const int ignore1, const int ignore2, const int ignore3) const
{
	for (int i=0; i<nQ(); i++) {
		if (i==ignore1 || i==ignore2 || i==ignore3) continue;
		double q = Q(i)(x, y, z);

		// The accuracy factor contains also the error from the conversion
		// of the Quadric to Conics
		double a = Q(i).acc(x, y, z, SMALL1);

		// FIXME we have to keep the information if the point is ON the surface
		// then to go directly to the gradient check!!!!
		// could result in numerical precision errors for fancy quadrics

		if (q > a)
			return false;
		else
		if (q >= -a) {
			// check dot product with gradient
			double gx, gy, gz;
			Q(i).grad(x, y, z, &gx, &gy, &gz);
			double g = sqrt(gx*gx + gy*gy + gz*gz);
			if (dx*gx + dy*gy + dz*gz > SMALL*g)
				return false;
		}
	}
	return true;
} // inside

/**
 * Find if a point (x,y,z) is inside the body or not, using the
 * quadratic surface equations and the direction the point is moving
 * @param withinAcc	return any of the quads was checked on the boundary
 * @return: true  inside
 *          false outside
 * XXX I've tried with const Point&p and d but in the end it was slower
 *     than with the plain references x,y,z...
 */
bool GBody::inside2D(const double  x, const double  y, const double  z,
		  const double dx, const double dy, const double dz,
		  const double acc[],
#ifdef EXPERIMENTAL
		  bool *withinAcc,
#endif
		  const int ignore1, const int ignore2) const
{
#ifdef EXPERIMENTAL
	if (withinAcc) *withinAcc = false;
#endif
	for (int i=nQ()-1; i>=0; i--) {
		if (i==ignore1 || i==ignore2) continue;
		double q = Q(i)(x, y, z);

		// The accuracy factor contains also the error from the conversion
		// of the Quadric to Conics
		double a = acc[i] * Q(i).acc(x,y,z, SMALL1);

		// FIXME we have to keep the information if the point is ON the surface
		// then to go directly to the gradient check!!!!
		// could result in numerical precision errors for fancy quadrics
#ifdef EXPERIMENTAL
		if (withinAcc && InRange(-1000.0*a, q, 1000.0*a))
			*withinAcc = true;
#endif

		if (q > a)
			return false;
		else
		if (q >= -a) {
			// Do not search direction on the boundary while we are
			// in intersectBody. This might add additional vertices
			// but is better than to ignore vertices
			if (ignore1<0) {
				// we are on the boundary
				// check dot product with gradient
				double gx, gy, gz;
				Q(i).grad(x, y, z, &gx, &gy, &gz);
				double g = sqrt(gx*gx + gy*gy + gz*gz);
				if (dx*gx + dy*gy + dz*gz > SMALL*g)
					return false;
			}
		}
	}
	return true;
} // inside2D

/**
 * Find if a point (x,y,z) is inside the body or not, using the
 * quadratic surface equations.
 * @param x,y,z	location to search
 * @param ignore_a, ignore_b, ignore_c quads to ignore in the computation
 * @return:
 *	true if the point is inside or all the quads were ignored
 *	false if outside
 * WARNING is ignore1 is not supplied then if the point is on the boundary
 *         then it is checked against the normal
 */
bool GBody::inside( const double x, const double y, const double z,
		const Quad *ignore_a, const Quad *ignore_b, const Quad *ignore_c
		) const
{
	for (int i=0; i<nQ(); i++) {
		const Quad &quad = Q(i);
		if (	&quad == ignore_a ||
			&quad == ignore_b ||
			&quad == ignore_c) {
			continue;
		}
		double q = quad(x,y,z);

		// Divide by gradient
		double gx, gy, gz, g;
		quad.grad(x, y, z, &gx, &gy, &gz);
		g = sqrt(gx*gx + gy*gy + gz*gz);
		if (g > SMALL) q /= g;
		double a = quad.acc(x,y,z, SMALL2);

		if (q > a)
			return false;
	}
	return true;
} // inside

/**
 * Find if a point (x,y,z) is outside the body or not, using the
 * quadratic surface equations.
 * @param x,y,z	location to search
 * @param ignore_a, ignore_b, ignore_c quads to ignore in the computation
 * @return:
 *	true if the point is outside or all the quads were ignored
 *	false if inside
 */
bool GBody::outside(
		const double x, const double y, const double z,
		const Quad *ignore_a, const Quad *ignore_b, const Quad *ignore_c
		) const
{
	for (int i=0; i<nQ(); i++) {
		const Quad &quad = Q(i);
		if (	&quad == ignore_a ||
			&quad == ignore_b ||
			&quad == ignore_c) {
			continue;
		}
		double q = quad(x,y,z);
		double gx, gy, gz, g;
		quad.grad(x, y, z, &gx, &gy, &gz);
		g = sqrt(gx*gx + gy*gy + gz*gz);
		if (g > SMALL) q /= g;
		double a = quad.acc(x,y,z, SMALL2);

		if (q < -a)
			return false;
	}
	return true;
} // outside


/** intersectRay
 * calculates the tmin, tmax and inverse from the intersection
 * of the ray with the body
 * @param x,y,z		position vector
 * @param dx,dy,dz	direction vector (has to be normalized)
 * @param tmin, tmax	entry and exit point to body
 * @return		true/false inverse should be used
 */
bool GBody::intersectRay(const double  x, const double  y, const double  z,
			 const double dx, const double dy, const double dz,
			 double *tmin, double *tmax) const
{
	double t2min, t2max;
	bool two = false;
	bool inv = false;

	*tmin = t2min = -INFINITE;
	*tmax = t2max =  INFINITE;

	for (int i=0; i<nQ(); i++) {
		int n;
		double ta, tb;
		if ((n=Q(i).intersect(x,y,z, dx,dy,dz, &ta, &tb))) {
			if (n<0) {
				assert(n==-2 && i==0);
				// Only for QUA surface
				if (nQ()==1) {
					// only the QUA are allowed to be concave
					*tmin = ta;
					*tmax = tb;
					return true;
				} else {
					// e.g. TRC
					// two ranges [-inf..ta] and [tb..+inf]
					two = true;
					*tmax = ta;
					t2min = tb;
				}
			} else {
				if (ta > *tmin) *tmin = ta;
				if (tb < *tmax) *tmax = tb;
				if (two) {
					if (ta > t2min) t2min = ta;
					if (tb < t2max) t2max = tb;
					if (t2min>t2max)
						two = false;	// invalid second region
				}
			}
			if (*tmin>*tmax) {	// first region is invalid
				if (two) {
					// move second region to first
					*tmin = t2min;
					*tmax = t2max;
					two  = false;
				}
				if (*tmin>*tmax) {
					// if still first is invalid then return false
					*tmin =  INFINITE;
					*tmax = -INFINITE;
					return false;
				}
			}
		} else
		if (Q(i)(x,y,z)>SMALL) {
			// no intersection and we are outside
			*tmin =  INFINITE;
			*tmax = -INFINITE;
			return false;
		}
	}
	assert(!two);

	return inv;
} // intersectRay

/** normal
 * @param r	position vector
 * @return	normal from the closest surface
 */
Vector GBody::normal(const Point& r) const
{
	double g, dx, dy, dz;
	int small = -1;
	double qsmall = INFINITE;
	for (int i=0; i<nQ(); i++) {
		double q = Q(i)(r.x,r.y,r.z);
		Q(i).grad(r.x, r.y, r.z, &dx, &dy, &dz);
		g = sqrt(Sqr(dx) + Sqr(dy) + Sqr(dz));
		if (g > SMALL) {
			g = 1.0 / g;
			q *= g;
		}

		double a = Q(i).acc(r.x,r.y,r.z, SMALL4);
		q = Abs(q);
		if (q<=a)
			return Vector(dx*g, dy*g, dz*g);
		if (q < qsmall) {
			qsmall = q;
			small = i;
		}
	}

	// normalize and return the smallest value
	Q(small).grad(r.x, r.y, r.z, &dx, &dy, &dz);
	g = sqrt(Sqr(dx) + Sqr(dy) + Sqr(dz));
	if (g > SMALL) {
		g = 1.0 / g;
		return Vector(dx*g, dy*g, dz*g);
	} else
		return Vector::Zo;
} // normal

/* ---------------------------- GPLABody ------------------------------- */
void GPLABody::setWhat(double *what, char *err)
{
	fixWhat(what, nWhat(), SMALL);
	xlen = ylen = infinite;
	zlen = 1.0;
	switch (type()) {
		case PLAbody:
			Z.set(what[0], what[1], what[2]);
			zlen = Z.normalize();
			if (Eq0(zlen,SMALL)) {
				if (err) strcpy(err,"PLA normal is null");
				Z = Vector::Zo;
				zlen = 1.0;
			}
			P.set(what[3], what[4], what[5]);
			findXYZ(err);
			break;

		case XYPbody:
			Z = Vector::Zo;
			X = Vector::Xo;
			Y = Vector::Yo;
			P.set(0.0, 0.0, what[0]);
			break;

		case XZPbody:
			Z = Vector::Yo;
			X = Vector::Zo;
			Y = Vector::Xo;
			P.set(0.0, what[0], 0.0);
			break;

		case YZPbody:
			Z = Vector::Xo;
			X = Vector::Yo;
			Y = Vector::Zo;
			P.set(what[0], 0.0, 0.0);
			break;

		default: ;
	}
} // setWhat

/* --- getWhat --- */
int GPLABody::getWhat(double *what) const
{
	switch (type()) {
		case PLAbody:
			what[0] = zlen*Z.x;
			what[1] = zlen*Z.y;
			what[2] = zlen*Z.z;

			what[3] = P.x;
			what[4] = P.y;
			what[5] = P.z;
			return 6;

		case XYPbody:
			what[0] = P.z;
			return 1;

		case XZPbody:
			what[0] = P.y;
			return 1;

		case YZPbody:
			what[0] = P.x;
			return 1;

		default: return -1;
	}
} // getWhat

/** position
 * @param x,y,z	new position to set
 */
void GPLABody::position(const Point& r)
{
	GBody::position(r);

	switch (type()) {
		case XYPbody:
			P.x = P.y = 0.0;
			break;

		case XZPbody:
			P.x = P.z = 0.0;
			break;

		case YZPbody:
			P.y = P.z = 0.0;
			break;

		default: ;
	}
	Po = P + zlen*Z;
} // position

/** return sign of the w axis that has to be used for move=20 type
 * see below on move
 */
void GPLABody::signMove20(const Point& r, const Vector& w)
{
	Point tr = _hasMatrix? _invMatrix*r : r;
	tr -= P;
	Vector dw = tr ^ w;
	_sign20 = Sign(dw*Z);
} // signMove20

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 * WARNING: the w is used for the rotation (20)
 *          and should be either viewport w or -w such as
 *                    ((r-p)^w)*Z > 0
 *          The reason is that the system below gives two solutions
 *          depending on the side that the user start to drag the extremity
 *          and the program has to remember th side
 */
void GPLABody::move(int item, const Point& r, const Vector& w)
{
	Point tr = _hasMatrix? _invMatrix*r : r;
	tr -= SP;
	if (item==20) {
		// Rotate plane using r as a Second node lying on the plane
		// Calculate new normal
		//Vector w = SZ ^ tr;
		//Z = tr ^ w;

		// We have to rotate the plane around w of the viewport
		// keeping the dot(W,N')=const
		// Equations to fulfill
		// 1. N'.r - N'.p = 0 <=> N'*tr = 0
		// 2. W.N' = W.N
		// 3. |N'| = |N|

		Point&  d = tr;
		Vector dw = d ^ w;
		double  c = SZ * w;

		double E = d.x*dw.z;
		if (!Eq0(E,SMALL)) {
			double A = c*d.x*d.y;
			double B = d.y*dw.y + d.z*dw.z;
			double C = c*d.x*d.x;
			double D = d.x*dw.y;

			double AA = B*B + D*D + E*E;
			double BB = 2.0 * (A*B + C*D);
			double CC = A*A + C*C - E*E;

			if (Eq0(AA,SMALL)) return;
			double x1, x2;
			int n = quadratic(BB/AA, CC/AA, &x1, &x2, SMALL);
			if (n==0) return;
			Z.z = x1;
			Z.y =  (C + D*Z.z) / E;
			Z.x = -(A + B*Z.z) / E;

			if (Sign(dw * Z)*_sign20 < 0) {
				Z.z = x2;
				Z.y =  (C + D*Z.z) / E;
				Z.x = -(A + B*Z.z) / E;
			}
		} else {
			// rotate x->y->z->x
			E = d.y*dw.x;
			if (!Eq0(E,SMALL)) {
				double A = c*d.y*d.z;
				double B = d.z*dw.z + d.x*dw.x;
				double C = c*d.y*d.y;
				double D = d.y*dw.z;

				double AA = B*B + D*D + E*E;
				double BB = 2.0 * (A*B + C*D);
				double CC = A*A + C*C - E*E;

				if (Eq0(AA,SMALL)) return;
				double x1, x2;
				int n = quadratic(BB/AA, CC/AA, &x1, &x2, SMALL);
				if (n==0) return;
				Z.x = x1;
				Z.z =  (C + D*Z.x) / E;
				Z.y = -(A + B*Z.x) / E;

				if (Sign(dw * Z)*_sign20 < 0) {
					Z.x = x2;
					Z.z =  (C + D*Z.x) / E;
					Z.y = -(A + B*Z.x) / E;
				}
			} else {
				// rotate x->z, y->x, z->y
				E = d.z*dw.y;
				if (Eq0(E,SMALL)) return;
				double A = c*d.z*d.x;
				double B = d.x*dw.x + d.y*dw.y;
				double C = c*d.z*d.z;
				double D = d.z*dw.x;

				double AA = B*B + D*D + E*E;
				double BB = 2.0 * (A*B + C*D);
				double CC = A*A + C*C - E*E;

				if (Eq0(AA,SMALL)) return;
				double x1, x2;
				int n = quadratic(BB/AA, CC/AA, &x1, &x2, SMALL);
				if (n==0) return;
				Z.y = x1;
				Z.x =  (C + D*Z.y) / E;
				Z.z = -(A + B*Z.y) / E;

				if (Sign(dw * Z)*_sign20 < 0) {
					Z.y = x2;
					Z.x =  (C + D*Z.y) / E;
					Z.z = -(A + B*Z.y) / E;
				}
			}
		}
		assert(Eq(Z*w,c,SMALL));

		Z.normalize();
		findXYZ(NULL);
	} else
	if (item==1)
		P = SP + (tr*SZ)*SZ;
	else
		GBody::move(item, r, w);
	checkType();
} // move

/** rotate body around an axis */
void GPLABody::rotate(const double angle, const Vector& axis)
{
	GBody::rotate(angle,axis);
	checkType();
} // rotate

/** check plane type */
void GPLABody::checkType()
{
	// Check type
	switch (Z.direction(SMALL)) {
		case 1: _type = YZPbody;
			break;
		case 2: _type = XZPbody;
			break;
		case 3: _type = XYPbody;
			break;
		default:
			_type = PLAbody;
	}
} // checkType

/** createQuads */
void GPLABody::createQuads()
{
	//      Cx     Cy     Cz      C
	addQuad(Z.x, Z.y, Z.z, -Z*P);
} // createQuads

/* createMesh */
void GPLABody::createMesh()
{
	bool first = mesh.isEmpty();
	if (first) mesh.allocateVertices(10);

	//       ^ Y
	//       |
	//   6---7---8
	//   |   |   |
	//   3---4---5  ---> X
	//   |   |   |
	//   0---1---2

	mesh.vertex(0) = P - xlen*X - ylen*Y;
	mesh.vertex(1) = P          - ylen*Y;
	mesh.vertex(2) = P + xlen*X - ylen*Y;

	mesh.vertex(3) = P - xlen*X;
	mesh.vertex(4) = P;
	mesh.vertex(5) = P + xlen*X;

	mesh.vertex(6) = P - xlen*X + ylen*Y;
	mesh.vertex(7) = P          + ylen*Y;
	mesh.vertex(8) = P + xlen*X + ylen*Y;

	mesh.vertex(9) = P + zlen*Z;
	mesh.calcBbox();

	if (first) {
		mesh.add(0, 1, 3, true,  false, true );
		mesh.add(1, 4, 3, true,  true,  false);

		mesh.add(1, 2, 4, true,  false, true );
		mesh.add(2, 5, 4, true,  true,  false);

		mesh.add(3, 4, 6, true,  false, true );
		mesh.add(4, 7, 6, true,  true,  false);

		mesh.add(4, 5, 7, true,  false, true );
		mesh.add(5, 8, 7, true,  true,  false);

		mesh.add(4, 9, true);	// add edge for normal

		mesh.process();
	}
} // createMesh

/** @return bounding box of body */
BBox GPLABody::_bbox() const
{
	BBox bb;
	bb.infinite();

	Point p = GBody::position();

	switch (vectorZ().direction(SMALL)) {
		case -1:
			bb.lowPt.x = p.x;
			break;

		case 1:
			bb.highPt.x = p.x;
			break;

		case -2:
			bb.lowPt.y = p.y;
			break;

		case 2:
			bb.highPt.y = p.y;
			break;

		case -3:
			bb.lowPt.z = p.z;
			break;

		case 3:
			bb.highPt.z = p.z;
			break;
	}

	return bb;
} // bbox

/** @return oriented bounding box of body */
OBBox* GPLABody::updateOBB(bool) const
{
	OBBox* bb = new OBBox();
	bb->infinite();

	bb->P = GBody::position();
	bb->X = vectorX();
	bb->Y = vectorY();
	bb->Z = vectorZ();

	bb->highPt.z = 0.0;

	return bb;
} // obbox

/**
 * @param body	other body to check location against
 * @return location of this body with respect to another body
 *
 * FIXME TODO add bounding box limits!!!!
 */
Location GPLABody::_locationWrt(const GBody *body) const
{
	if (body->type() == PLAbody || body->type() == XYPbody ||
	    body->type() == XZPbody || body->type() == YZPbody) {
		// Check if normals are co-linear
		Vector thisN = vectorZ();
		Vector bodyN = body->vectorZ();

		double ptot2 = thisN.length2()*bodyN.length2();
		if (ptot2<=0.0) return LOCATION_OVERLAP;	// error
		double dot = (thisN*bodyN)/sqrt(ptot2);

		// if |dot|!=1 they will overlap somewhere
		if (Eq(dot,  1.0, SMALL)) {
			// Normals are parallel
			// Two situations AinB or BinA
			// Check if this is inside body or vice versa
			double pos = Q(0)(body->position());
			if (pos < 0.0)
				return LOCATION_BinA;
			else
				return LOCATION_AinB;
		} else
		if (Eq(dot, -1.0, SMALL)) {
			// Normals are anti-parallel
			// Two situations Overlap or Outside
			double pos = Q(0)(body->position());
			if (pos < 0.0)
				// FIXME could be optimized in the level of Zone::optimize
				// but here we can only return an OVERLAP
				return LOCATION_OVERLAP;
			else
				return LOCATION_OUTSIDE;
		} else
			return LOCATION_OVERLAP;

	} else
	if (body->type() == SPHbody) {
		// Check distance of SPH center to plane vs radius
		const Quad& q = Q(0);
		double d = q(body->position()) / sqrt(Sqr(q.Cx) + Sqr(q.Cy) + Sqr(q.Cz));
		double r = ((const GSPHBody*)body)->R();
		if (d > r)
			return LOCATION_OUTSIDE;
		else
		if (d < -r)
			return LOCATION_BinA;
		else
			return LOCATION_OVERLAP;
	} else
	if (body->type() == BOXbody || body->type() == RPPbody ||
	    body->type() == WEDbody || body->type() == RAWbody ||
	    body->type() == ARBbody || body->type() == TETbody { //For TET, added by zxw
		// FIXME normally this check could be done for every object
		// assuming a safe sagita distance
		double n = sqrt(Sqr(Q(0).Cx) + Sqr(Q(0).Cy) + Sqr(Q(0).Cz));
		double d = Q(0)(body->mesh.vertex(0)) / n;
		if (Eq0(d, SMALL)) return LOCATION_OVERLAP;
		bool ins = d < 0.0;
		for (int i=1; i<body->mesh.nvertices(); i++) {
			d = Q(0)(body->mesh.vertex(i)) / n;
			if (Eq0(d, SMALL)) return LOCATION_OVERLAP;
			if (ins) {
				if (d > 0.0) return LOCATION_OVERLAP;
			} else
				if (d < 0.0) return LOCATION_OVERLAP;
		}
		if (ins)
			return LOCATION_BinA;
		else
			return LOCATION_OUTSIDE;
	} else
	if (body->type() == XCCbody || body->type() == YCCbody ||
	    body->type() == ZCCbody) {
		// Verify if axis of cylinder is perpendicular to plane
		Vector thisN = vectorZ();
		Vector bodyZ = body->vectorZ();
		if (Eq0(thisN * bodyZ, SMALL)) {
			// Check distance of cylinder axis to plane
			const Quad& q = Q(0);
			double d = q(body->position()) / sqrt(Sqr(q.Cx) + Sqr(q.Cy) + Sqr(q.Cz));
			double r = ((const GInfEllCylBody*)body)->R();
			if (d > r)
				return LOCATION_OUTSIDE;
			else
			if (d < -r)
				return LOCATION_BinA;
			else
				return LOCATION_OVERLAP;
		} else
			return LOCATION_OVERLAP;
//	} else
//	if (body->type() == QUAbody) {
//		Plane: ax + by + cz + d = 0
//		Substituting z = - (ax+by+d)/c  (if c!=0) to QUA
//		we get a Conic equation on QUA(x,y)
//		Check if Conic is degenerate or not!
	} else {
		BBox bb = body->bbox();
		double n = sqrt(Sqr(Q(0).Cx) + Sqr(Q(0).Cy) + Sqr(Q(0).Cz));

		double d = Q(0)(bb.low()) / n;
		if (Eq0(d, SMALL)) return LOCATION_OVERLAP;
		bool ins = d < 0.0;

		Point p(bb.lowPt.x, bb.lowPt.y, bb.highPt.z);
		d = Q(0)(p) / n;
		if (Eq0(d, SMALL)) return LOCATION_OVERLAP;
		if (ins) {
			if (d > 0.0) return LOCATION_OVERLAP;
		} else
			if (d < 0.0) return LOCATION_OVERLAP;

		p(bb.lowPt.x, bb.highPt.y, bb.highPt.z);
		d = Q(0)(p) / n;
		if (Eq0(d, SMALL)) return LOCATION_OVERLAP;
		if (ins) {
			if (d > 0.0) return LOCATION_OVERLAP;
		} else
			if (d < 0.0) return LOCATION_OVERLAP;

		p(bb.lowPt.x, bb.highPt.y, bb.lowPt.z);
		d = Q(0)(p) / n;
		if (Eq0(d, SMALL)) return LOCATION_OVERLAP;
		if (ins) {
			if (d > 0.0) return LOCATION_OVERLAP;
		} else
			if (d < 0.0) return LOCATION_OVERLAP;

		p(bb.highPt.x, bb.lowPt.y, bb.lowPt.z);
		d = Q(0)(p) / n;
		if (Eq0(d, SMALL)) return LOCATION_OVERLAP;
		if (ins) {
			if (d > 0.0) return LOCATION_OVERLAP;
		} else
			if (d < 0.0) return LOCATION_OVERLAP;

		p(bb.highPt.x, bb.lowPt.y, bb.highPt.z);
		d = Q(0)(p) / n;
		if (Eq0(d, SMALL)) return LOCATION_OVERLAP;
		if (ins) {
			if (d > 0.0) return LOCATION_OVERLAP;
		} else
			if (d < 0.0) return LOCATION_OVERLAP;

		p(bb.highPt.x, bb.highPt.y, bb.highPt.z);
		d = Q(0)(p) / n;
		if (Eq0(d, SMALL)) return LOCATION_OVERLAP;
		if (ins) {
			if (d > 0.0) return LOCATION_OVERLAP;
		} else
			if (d < 0.0) return LOCATION_OVERLAP;

		p(bb.highPt.x, bb.highPt.y, bb.lowPt.z);
		d = Q(0)(p) / n;
		if (Eq0(d, SMALL)) return LOCATION_OVERLAP;
		if (ins) {
			if (d > 0.0) return LOCATION_OVERLAP;
		} else
			if (d < 0.0) return LOCATION_OVERLAP;

		if (ins)
			return LOCATION_BinA;
		else
			return LOCATION_OUTSIDE;
	}
	return GBody::_locationWrt(body);
} // _locationWrt

/* ---------------------------- GSPHBody ------------------------------- */
void GSPHBody::setWhat(double *what, char *err)
{
	fixWhat(what, nWhat(), SMALL);
	P.set(what[0], what[1], what[2]);
	if (err && what[3]<TOOSMALL) strcpy(err,"Negative or zero radius");
	R(Max(what[3],0.0));
} // setWhat

/* --- getWhat --- */
int GSPHBody::getWhat(double *what) const
{
	what[0] = P.x;
	what[1] = P.y;
	what[2] = P.z;
	what[3] = R();
	return 4;
} // getWhat

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GSPHBody::move(int item, const Point& r, const Vector& w)
{
	if (item==1) {
		Point tr = _hasMatrix? _invMatrix*r : r;
		R((tr-SP).length());
	} else
		GBody::move(item, r, w);
} // move

/** createQuads */
void GSPHBody::createQuads()
{
	//  Cxx  Cyy  Czz  Cxy  Cxz  Cyz
	addQuad(1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
	//   Cx       Cy       Cz
	   -2.0*P.x, -2.0*P.y, -2.0*P.z, P.length2() - SQR(R()));
} // createQuads

/** @return bounding box of body */
BBox GSPHBody::_bbox() const
{
	BBox bb;
	Point p = GBody::position();
	bb.add(p.x-R(), p.y-R(), p.z-R());
	bb.add(p.x+R(), p.y+R(), p.z+R());
	return bb;
} // bbox

/** @return oriented bounding box of body */
OBBox* GSPHBody::updateOBB(bool in) const
{
	OBBox* bb = new OBBox();

	bb->P = GBody::position();

	if (!in) {
		bb->low()  = Point(-R(), -R(), -R());
		bb->high() = Point(R(),   R(),  R());
	}
	else {
		// Add inner points
		double sc  = R()/Sqrt3;
		bb->low()  = Point(-sc, -sc, -sc);
		bb->high() = Point(+sc, +sc, +sc);
	}

	return bb;
} // obbox

/**
 * @param body	other body to check location against
 * @return location of this body with respect to another body
 */
Location GSPHBody::_locationWrt(const GBody *body) const
{
	if (body->type() == SPHbody) {
		// Minimum distance between two centers
		double dist = (position() - body->position()).length();

		if (dist + ((const GSPHBody*)body)->R() <  R())
			return LOCATION_BinA;
		else
		if (dist - ((const GSPHBody*)body)->R() < -R())
			return LOCATION_AinB;
		else
		if (dist - ((const GSPHBody*)body)->R() >  R())
			return LOCATION_OUTSIDE;
		else
			return LOCATION_OVERLAP;
	} else
	if (body->type() == XCCbody || body->type() == YCCbody ||
	    body->type() == ZCCbody) {
		// Minimum distance between center and axes
		double dist = pointLineDistance(position(),
			body->position(), body->vectorZ());

		if (dist - ((const GInfEllCylBody*)body)->R() >  R())
			return LOCATION_OUTSIDE;
		else
			return LOCATION_OVERLAP;
	}

	return GBody::_locationWrt(body);
} // _locationWrt

/* ---------------------------- GELLBody ------------------------------- */
void GELLBody::setWhat(double *what, char *err)
{
	fixWhat(what, nWhat(), SMALL);
	Point F1(what[0], what[1], what[2]);		// focus 1
	Point F2(what[3], what[4], what[5]);		// focus 2
	zlen = what[6]/2.0;				// major axis
	if (err && zlen<=TOOSMALL) strcpy(err,"Zero or negative major axis");
	P  = 0.5*(F1+F2);				// center
	Z  = F2 - F1;
	double f = Z.normalize()/2.0;			// foci half-distance
	if (err && what[6]<=f) strcpy(err,"Full length of major axis smaller than focii distance");
	double x2 = (zlen-f)*(zlen+f);
	xlen = ylen = (x2>0.0?Sqrt(x2):0.0);		// minor axes

	findXYZ(err);
} // setWhat

/* --- getWhat --- */
int GELLBody::getWhat(double *what) const
{
	double f = Sqrt((zlen+xlen)*(zlen-xlen));
	Point F1 = P - f*Z;
	Point F2 = P + f*Z;
	what[0] = ROUND2INT(F1.x, SMALL1);
	what[1] = ROUND2INT(F1.y, SMALL1);
	what[2] = ROUND2INT(F1.z, SMALL1);
	what[3] = ROUND2INT(F2.x, SMALL1);
	what[4] = ROUND2INT(F2.y, SMALL1);
	what[5] = ROUND2INT(F2.z, SMALL1);
	what[6] = 2.0*zlen;
	return 7;
} // getWhat

/** @return transformed position */
Point GELLBody::node(const int n) const
{
	double f = Sqrt((zlen+xlen)*(zlen-xlen));
	Point r;
	switch (n) {
		case 0: r = P;
			break;
		case 1: r = P - f*Z;
			break;
		case 2: r = P + f*Z;
			break;
	}
	return _hasMatrix? _matrix*r : r;
} // node

/**  @return closest point, edge, face */
int GELLBody::closest(const Point& r, const double d, const Vector& w) const
{
	Point tr = _hasMatrix? _invMatrix*r : r;

	double f = Sqrt((zlen+xlen)*(zlen-xlen));
	if ((P-f*Z - tr).length2() < Sqr(d)) return 16;
	if ((P+f*Z - tr).length2() < Sqr(d)) return 17;

	int close = GBody::closest(r,d,w);
	if (close == 1) {
		double z2 = Sqr(tr*SZ);
		double xy2 = tr.length2() - z2;
		if (xy2>z2) close = 2;
	}
	return close;
} // closest

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GELLBody::move(int item, const Point& r, const Vector& w)
{
	Point tr = _hasMatrix? _invMatrix*r : r;
	if (item==1) {
		tr -= SP;
		double z2 = Sqr(tr*SZ);
		double xy2 = tr.length2() - z2;
		double a = z2 / (1.0-xy2/Sqr(sxlen));
		zlen = a>0.0? Sqrt(a) : szlen;
		// Check if we have a valid ELL
		if (zlen < xlen) zlen = szlen;
	} else
	if (item==2) {
		tr -= SP;
		double z2 = Sqr(tr*SZ);
		double xy2 = tr.length2() - z2;
		double a = xy2 / (1.0-z2/Sqr(szlen));
		xlen = ylen = a>0.0? Sqrt(a) : sxlen;
		if (zlen < xlen) xlen = ylen = sxlen;
	} else
	if (item==16) {
		double what[7];
		getWhat(what);
		what[0] = tr.x;
		what[1] = tr.y;
		what[2] = tr.z;
		setWhat(what,NULL);
	} else
	if (item==17) {
		double what[7];
		getWhat(what);
		what[3] = tr.x;
		what[4] = tr.y;
		what[5] = tr.z;
		setWhat(what,NULL);
	} else
		GBody::move(item, r, w);
} // move

/** createQuads */
void GELLBody::createQuads()
{
	if (ISZERO(xlen) || ISZERO(ylen) || ISZERO(zlen)) return;

	// add quad centered to origin
	//      Cxx        Cyy     Czz     Cxy  Cxz  Cyz  Cx   Cy   Cz     C
	addQuad(1.0/Sqr(xlen), 1.0/Sqr(ylen), 1.0/Sqr(zlen), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 );

	// rotate to align to F1-F2 direction
	Matrix4 T, M, Minv;
	Minv.make(X,Y,Z);		// rotation matrix
	Minv.inverse();
	T.translate(P);			// translation to P
	M.multiply(T, Minv);
	Minv.inverse(M);

	// transform first quad
	_Q[0].transform(Minv);
	_Q[0].normalize();
} // createQuads

/** @return bounding box of body */
BBox GELLBody::_bbox() const
{
	return Q(0).bbox();
} // bbox

/** @return oriented bounding box of body */
OBBox* GELLBody::updateOBB(bool in) const
{
	OBBox* bb = new OBBox();

	bb->P = GBody::position();
	bb->X = vectorX();
	bb->Y = vectorY();
	bb->Z = vectorZ();

	if (!in) {
		bb->low()  = Point(-xlen, -ylen, -zlen);
		bb->high() = Point( xlen,  ylen,  zlen);
	} else {
		double sc = 1.0/Sqrt3;

		bb->low()  = Point(-xlen*sc, -ylen*sc, -zlen*sc);
		bb->high() = Point( xlen*sc,  ylen*sc,  zlen*sc);
	}

	return bb;
} // bbox

/* ---------------------------- GBOXBody ------------------------------- */
void GBOXBody::setWhat(double *what, char *err)
{
	fixWhat(what, nWhat(), SMALL);
	if (type() == RPPbody) {
		P.set( what[0], what[2], what[4]);
		Po.set(what[1], what[3], what[5]);

		X = Vector::Xo;
		Y = Vector::Yo;
		Z = Vector::Zo;

		xlen = Po.x - P.x;
		ylen = Po.y - P.y;
		zlen = Po.z - P.z;
	} else {
		P.set(what[0], what[ 1], what[ 2]);
		Z.set(what[3], what[ 4], what[ 5]);
		X.set(what[6], what[ 7], what[ 8]);
		Y.set(what[9], what[10], what[11]);

		// outer corner
		Po = P + X + Y + Z;

		xlen = X.normalize();
		ylen = Y.normalize();
		zlen = Z.normalize();

		checkOrthogonal(err);
	}
} // setWhat

/* --- getWhat --- */
int GBOXBody::getWhat(double *what) const
{
	if (type() == RPPbody) {
		what[0] = Min(P.x, Po.x);
		what[1] = Max(P.x, Po.x);
		what[2] = Min(P.y, Po.y);
		what[3] = Max(P.y, Po.y);
		what[4] = Min(P.z, Po.z);
		what[5] = Max(P.z, Po.z);
		return 6;
	} else {
		what[ 0] = P.x;
		what[ 1] = P.y;
		what[ 2] = P.z;
		what[ 3] = zlen*Z.x;
		what[ 4] = zlen*Z.y;
		what[ 5] = zlen*Z.z;
		what[ 6] = xlen*X.x;
		what[ 7] = xlen*X.y;
		what[ 8] = xlen*X.z;
		what[ 9] = ylen*Y.x;
		what[10] = ylen*Y.y;
		what[11] = ylen*Y.z;
		return 12;
	}
} // getWhat

/** position
 * @param x,y,z	new position to set
 */
void GBOXBody::position(const Point& r)
{
	GBody::position(r);
	Po = P + xlen*X + ylen*Y + zlen*Z;
} // position

/** rotate body around an axis */
void GBOXBody::rotate(const double angle, const Vector& axis)
{
	GBody::rotate(angle,axis);
	Po = P + xlen*X + ylen*Y + zlen*Z;

	// Check type
	int ax = X.direction(SMALL);
	int ay = Y.direction(SMALL);
	int az = Z.direction(SMALL);

	if (ax && ay && az) {
		_type = RPPbody;

		X = Vector::Xo;
		Y = Vector::Yo;
		Z = Vector::Zo;

		if (P.x > Po.x) Swap(P.x, Po.x);
		if (P.y > Po.y) Swap(P.y, Po.y);
		if (P.z > Po.z) Swap(P.z, Po.z);
		xlen = Po.x - P.x;
		ylen = Po.y - P.y;
		zlen = Po.z - P.z;
	} else
		_type = BOXbody;
} // rotate

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GBOXBody::move(int item, const Point& r, const Vector& w)
{
	Point tr = _hasMatrix? _invMatrix*r : r;
	tr -= SP;
	if (item==1) { // Lower X-plane
		double xp = tr*SX;
		P = SP + xp*SX;
		xlen = sxlen - xp;
		if (xlen<0.0) {
			xlen = -xlen;
			P = SP + sxlen*SX;
		}
	} else
	if (item==2) {	// Lower Y-plane
		double yp = tr*SY;
		P = SP + yp*SY;
		ylen = sylen - yp;
		if (ylen<0.0) {
			ylen = -ylen;
			P = SP + sylen*SY;
		}
	} else
	if (item==3) { // Lower Z-plane
		double zp = tr*SZ;
		P = SP + zp*SZ;
		zlen = szlen - zp;
		if (zlen<0.0) {
			zlen = -zlen;
			P = SP + szlen*SZ;
		}
	} else {
		// Combined moves
		//    16 = YZ plane (X-perpendicular)
		//    17 = XZ plane (Y-perpendicular)
		//    18 = XY plane (Z-perpendicular)

		bool found = false;
		P = SP;
		if (item==4 || item==17 || item==18) { // Upper X-plane
			found = true;
			xlen = tr*SX;
			if (xlen<0.0) {
				xlen = -xlen;
				P -= xlen*SX;
			}
		}

		if (item==5 || item==16 || item==18) { // Upper Y-plane
			found = true;
			ylen = tr*SY;
			if (ylen<0.0) {
				ylen = -ylen;
				P -= ylen*SY;
			}
		}

		if (item==6 || item==16 || item==17) { // Upper Z-plane
			found = true;
			zlen = tr*SZ;
			if (zlen<0.0) {
				zlen = -zlen;
				P -= zlen*SZ;
			}
		}

		if (!found)
			GBody::move(item, r, w);
	}
	Po = P + xlen*X + ylen*Y + zlen*Z;
} // move

/** createQuads */
void GBOXBody::createQuads()
{
	if (type() == RPPbody) {
		//        Cx   Cy   Cz   C
		addQuad(-1.0, 0.0, 0.0, P.x);
		addQuad( 0.0,-1.0, 0.0, P.y);
		addQuad( 0.0, 0.0,-1.0, P.z);

		addQuad( 1.0, 0.0, 0.0,-Po.x);
		addQuad( 0.0, 1.0, 0.0,-Po.y);
		addQuad( 0.0, 0.0, 1.0,-Po.z);
	} else {
		// first 3 planes will be
		//        Cx     Cy     Cz      C
		addQuad(-X.x,-X.y,-X.z, X*P);
		addQuad(-Y.x,-Y.y,-Y.z, Y*P);
		addQuad(-Z.x,-Z.y,-Z.z, Z*P);

		// last 3 planes will be
		addQuad( X.x, X.y, X.z,-X*Po);
		addQuad( Y.x, Y.y, Y.z,-Y*Po);
		addQuad( Z.x, Z.y, Z.z,-Z*Po);
	}
} // createQuads

/* createMesh */
void GBOXBody::createMesh()
{
	if (mesh.isEmpty())
		Mesh::createParallelepiped(mesh, P, xlen*X, ylen*Y, zlen*Z);
	else {
		// update vertices
		mesh.vertex(0) = P;			// 0
		mesh.vertex(1) = P + xlen*X;		// 1
		mesh.vertex(2) = P + xlen*X + ylen*Y;	// 2
		mesh.vertex(3) = P          + ylen*Y;	// 3

		mesh.vertex(4) = P                   + zlen*Z; // 4
		mesh.vertex(5) = P + xlen*X          + zlen*Z; // 5
		mesh.vertex(6) = P + xlen*X + ylen*Y + zlen*Z; // 6
		mesh.vertex(7) = P + ylen*Y          + zlen*Z; // 7
		mesh.calcBbox();
	}
} // GBOXBody

/** @return oriented bounding box of body */
OBBox* GBOXBody::updateOBB(bool) const
{
	OBBox* bb = new OBBox();
	bb->P = GBody::position();
	bb->X = vectorX();
	bb->Y = vectorY();
	bb->Z = vectorZ();
	bb->low()  = Point(0.0,  0.0,  0.0);
	bb->high() = Point(xlen, ylen, zlen);
	return bb;
} // obbox

/**
 * @param body	other body to check location against
 * @return location of this body with respect to another body
 */
Location GBOXBody::_locationWrt(const GBody *body) const
{
#if 0
	if (body->type() == SPHbody) {
		Location bodyloc = LOCATION_INSIDE;
		//bodysize holds half of the size of the BOX in each of its dimensions
		Vector bodysize = Vector(vectorXlen()/2, vectorYlen()/2, vectorZlen()/2);
		// Rsph locates the SPHERE center w.r.t the center of the BOX
		Vector Rsph = GBody::position() - (body->position() + bodysize);
		//d holds the closest boundary (vertex, edge) in the direction of the SPHERE position
		Vector d;
		d.x = Rsph.x ? bodysize.x * Rsph.x / Abs(Rsph.x) : 0.0;
		d.y = Rsph.y ? bodysize.y * Rsph.y / Abs(Rsph.y) : 0.0;
		d.z = Rsph.z ? bodysize.z * Rsph.z / Abs(Rsph.z) : 0.0;
		//if the SPHERE is outside the BOX...
		if (Rsph.length() >  d.length()) {
			//dist1 holds distance from SPHERE center to BOX's closest boundary (vertex, edge)
			double dist1 = (Rsph - d).length();
			//dist2 holds distance from SPHERE center to BOX's farthest boundary (vertex, edge)
			double dist2 = (Rsph + d).length();

			if (((const GSPHBody*)body)->R() < dist1)
				bodyloc = LOCATION_OUTSIDE;
			else
			if (((const GSPHBody*)body)->R() < dist2)
				bodyloc = LOCATION_OVERLAP;
			else
				bodyloc = LOCATION_AinB; //check, should be BOX in SPHERE
		} else { //...else if the sphere is inside the BOX
			//not implemented
			bodyloc = LOCATION_OVERLAP;
			float Radio = ((const GSPHBody*)body)->R();
			if (Rsph.x + Radio < xlen)
				if (Rsph.y + Radio < ylen)
					if (Rsph.z + Radio < zlen)
						bodyloc = LOCATION_BinA; //check, should be SPHERE in BOX
		}
	} //else
	// if (body->type() == XCCbody || body->type() == YCCbody ||
	//     body->type() == ZCCbody) {
	//	// Minimum distance between center and axes
	//	double dist = PointLineDistance(position(),
	//		body->position(), body->vectorZ());

	//	if (dist - ((const GInfEllCylBody*)body)->R() >  R())
	//		return LOCATION_OUTSIDE;
	//	else
	//		return LOCATION_OVERLAP;
	// } else
	// if (body->type() == BOXbody || body->type() == RPPbody) {

	//	}
		// return bodyloc;
	// }//--------------------------------------------added by Leonel---------
#endif
	return GBody::_locationWrt(body);
} // _locationWrt


/* ---------------------------- GWEDBody ------------------------------- */
void GWEDBody::setWhat(double *what, char *err)
{
	fixWhat(what, nWhat(), SMALL);
	P.set(what[0], what[ 1], what[ 2]);
	X.set(what[3], what[ 4], what[ 5]);
	Y.set(what[6], what[ 7], what[ 8]);
	Z.set(what[9], what[10], what[11]);

	// Upper point
	Po = P + Z;

	// normalize vectors
	xlen = X.normalize();
	ylen = Y.normalize();
	zlen = Z.normalize();

	checkOrthogonal(err);
} // setWhat

/* --- getWhat --- */
int GWEDBody::getWhat(double *what) const
{
	what[ 0] = P.x;
	what[ 1] = P.y;
	what[ 2] = P.z;
	what[ 3] = xlen*X.x;
	what[ 4] = xlen*X.y;
	what[ 5] = xlen*X.z;
	what[ 6] = ylen*Y.x;
	what[ 7] = ylen*Y.y;
	what[ 8] = ylen*Y.z;
	what[ 9] = zlen*Z.x;
	what[10] = zlen*Z.y;
	what[11] = zlen*Z.z;
	return 12;
} // getWhat

/** @return normal of the cutting plane */
Vector	GWEDBody::N() const
{
	// normal of the cutting plane
	Vector n = (ylen*Y - xlen*X) ^ Z;
	n.normalize();
	return n;
} // N

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GWEDBody::move(int item, const Point& r, const Vector& w)
{
	Point tr = _hasMatrix? _invMatrix*r : r;
	tr -= SP;
	switch (item) {
		case 1: {	// Lower X-plane
				double xp = tr*SX;
				P = SP + xp*SX;
				xlen = sxlen - xp;
				if (xlen<0.0) {
					X = -SX;
					xlen = -xlen;
				} else
					X = SX;
			}
			break;

		case 2: {	// Lower Y-plane
				double yp = tr*SY;
				P = SP + yp*SY;
				ylen = sylen - yp;
				if (ylen<0.0) {
					Y = -SY;
					ylen = -ylen;
				} else
					Y = SY;
			}
			break;

		case 3: {	// Lower Z-plane
				double zp = tr*SZ;
				P = SP + zp*SZ;
				zlen = szlen - zp;
				if (zlen<0.0) {
					Z = -SZ;
					zlen = -zlen;
				} else
					Z = SZ;
			}
			break;

		case 4:		// Upper Z-plane
			zlen = tr*SZ;
			if (zlen<0.0) {
				zlen = -zlen;
				P = SP - zlen*SZ;
			}
			break;

		case 5: {	// cutting plane
				// normal from the saved variables
				Vector n = (sylen*SY - sxlen*SX) ^ SZ;
				n.normalize();
				double d = tr*n;
				Point R = SP + d*n;
				if (d<0.0) {
					X = -SX;
					Y = -SY;
					d = -d;
				} else {
					X = SX;
					Y = SY;
				}
				double xr = X*(R-P);
				if (!Eq0(xr,SMALL3)) {
					xlen = Sqr(d) / xr;
					ylen = Sqr(d) / (Y*(R-P));
				}
			}
			break;

		default:
			GBody::move(item, r, w);
	}
	Po = P + zlen*Z;
} // move

/** createQuads */
void GWEDBody::createQuads()
{
	// first 3 planes will be
	//        Cx     Cy     Cz      C
	addQuad(-X.x,-X.y,-X.z, X*P);
	addQuad(-Y.x,-Y.y,-Y.z, Y*P);
	addQuad(-Z.x,-Z.y,-Z.z, Z*P);

	// upper Z plane
	addQuad( Z.x, Z.y, Z.z,-Z*Po);

	// normal of the cutting plane
	Vector n = N();

	// cutting plane
	double d = -n * (P+xlen*X);
	if (n*X < 0.0) {
		n = -n;
		d = -d;
	}
	addQuad( n.x, n.y, n.z, d);
} // createQuads

/* createMesh */
void GWEDBody::createMesh()
{
	bool first = mesh.isEmpty();
	if (first) mesh.allocateVertices(6);

	mesh.vertex(0) = P;
	mesh.vertex(1) = P  + xlen*X;
	mesh.vertex(2) = P  + ylen*Y;
	mesh.vertex(3) = Po;
	mesh.vertex(4) = Po + xlen*X;
	mesh.vertex(5) = Po + ylen*Y;
	mesh.calcBbox();

	if (first) {
		mesh.add(0, 2, 1, true,  true,  true );
		mesh.add(3, 4, 5, true,  true,  true );

		mesh.add(0, 4, 3, false, true,  true );
		mesh.add(0, 1, 4, true,  true,  false);

		mesh.add(0, 3, 5, true,  true,  false);
		mesh.add(0, 5, 2, false, true,  true );

		mesh.add(1, 2, 5, true,  true,  false);
		mesh.add(1, 5, 4, false, true,  true );

		mesh.process();
		assert(mesh.isClosed());
		assert(mesh.isOrientable());
#if _DEBUG>2
		cout << endl;
		cout << "WED Mesh ";
		cout << " isClosed=" << mesh.isClosed();
		cout << " isOrientable=" << mesh.isOrientable() << endl;
		cout << "WED volume=" << mesh.volume() << endl;
#endif
	}
} // createMesh

/** @return oriented bounding box of body */
OBBox* GWEDBody::updateOBB(bool in) const
{
	OBBox* bb = new OBBox();
	bb->P = GBody::position();
	bb->X = vectorX();
	bb->Y = vectorY();
	bb->Z = vectorZ();

	bb->low() = Point(0.0, 0.0, 0.0);

	if (!in)
		bb->high() = Point(xlen, ylen, zlen);
	else
		bb->high() = Point(xlen/2.0, ylen/2.0, zlen);

	return bb;
} // obbox

/* ---------------------------- GTETBody ------------------------------- */ //For TET, added by zxw
void GTETBody::setWhat(double* what, char* err)
{
	fixWhat(what, nWhat(), SMALL);

	mesh.allocateVertices(4);

	// add vertices
	P.set(what[3], what[4], what[5]); // P.set(what[0], what[1], what[2]);

	P1.set(what[0], what[1], what[2]);
	P2.set(what[6], what[7], what[8]);
	P3.set(what[9], what[10], what[11]);

	Vector u, v, w;
	u = P1 - P;
	v = P2 - P;
	w = P3 - P;

	X.set(u.x, u.y, u.z);
	Y.set(v.x, v.y, v.z);
	Z.set(w.x, w.y, w.z);

	// normalize vectors
	xlen = X.normalize();
	ylen = Y.normalize();
	zlen = Z.normalize();

	for (int i = 0; i < 4; i++)
		mesh.vertex(i) = Point(what[i * 3], what[i * 3 + 1], what[i * 3 + 2]);

	// checkOrthogonal(err); 
} // setWhat

/* --- getWhat --- */
int GTETBody::getWhat(double* what) const
{
	int w = 0;
	for (int i = 0; i < 4; i++)
	{
		what[w++] = mesh.vertex(i).x;
		what[w++] = mesh.vertex(i).y;
		what[w++] = mesh.vertex(i).z;
	}
	return 12;
} // getWhat

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GTETBody::move(int item, const Point& r, const Vector& w)
{
	Point dR = r - P;
	for (int i = 0; i < mesh.nvertices(); i++)
		mesh.vertex(i) += dR;

	GBody::move(item, r, w);
} // move
/** rotate body around an axis and current body position
 * @param angle	to rotate
 * @param axis	to rotate
 */
void GTETBody::rotate(const double angle, const Vector& axis)
{
	GBody::rotate(angle, axis);
	Matrix4 rot; // duplicate rot matrix with the GBody::rotate
	rot.rotate(angle, axis);
	for (int i = 1; i < 4; i++)
	{
		Vector v = mesh.vertex(i) - P;
		//		mesh.vertex(i) = rot.multVector(v) + P;
		mesh.vertex(i) = rot * v + P;
	}
} // rotate

/** create Quads */
void GTETBody::createQuads()
{
	Point V[4];
	if (mesh.isEmpty())
	{
		mesh.allocateVertices(4);
		return;
	}
	else
	{
		for (int j = 0; j < 4; j++)
		{
			V[j] = mesh.vertex(j);
		}
	}
	// create quads
	Vector N1 = (V[0] - V[1]) ^ (V[2] - V[1]);
	Vector N2 = (V[3] - V[1]) ^ (V[2] - V[1]);
	Vector N3 = (V[0] - V[1]) ^ (V[3] - V[1]);
	Vector N4 = (V[3] - V[2]) ^ (V[0] - V[2]);

	N1.normalize();
	N2.normalize();
	N3.normalize();
	N4.normalize();

	// 4 planes will be
	//        Cx     Cy     Cz      C
	// addQuad(N1.x, N1.y, N1.z, -N1 * V[1]);
	// addQuad(-N2.x, -N2.y, -N2.z, N2 * V[1]);
	// addQuad(-N3.x, -N3.y, -N3.z, N3 * V[1]);
	// addQuad(-N4.x, -N4.y, -N4.z, N4 * V[2]);
	addQuad(-N1.x, -N1.y, -N1.z, N1 * V[1]);
	addQuad(N2.x, N2.y, N2.z, -N2 * V[1]);
	addQuad(N3.x, N3.y, N3.z, -N3 * V[1]);
	addQuad(N4.x, N4.y, N4.z, -N4 * V[2]);
} // createQuads

/* createMesh */
void GTETBody::createMesh()
{
	// if (mesh.nedges())
	// 	return;
	Point V[4];
	for (int j = 0; j < 4; j++)
	{
		V[j] = mesh.vertex(j);
	}

	mesh.add(0, 1, 2, true, true, true);
	mesh.add(0, 2, 3, true, true, true);
	mesh.add(0, 3, 1, true, true, true);
	mesh.add(1, 3, 2, true, true, true);

	mesh.calcBbox();
	mesh.process();
	assert(mesh.isClosed());
	assert(mesh.isOrientable());
#if _DEBUG > 2
	cout << endl;
	cout << "TET Mesh ";
	cout << "isClosed=" << mesh.isClosed();
	cout << "isOrientable=" << mesh.isOrientable() << endl;
	cout << "TET volume=" << mesh.volume() << endl;
#endif
} // createMesh

/** @return bounding box of body */
BBox GTETBody::_bbox() const
{
	BBox bb;
	Point V[4];
	for (int j = 0; j < 4; j++)
	{
		V[j] = mesh.vertex(j);
	}
	Point p;
	p = V[1];

	bb.add(V[0].x - p.x, V[0].y - p.y, V[0].z - p.z);
	bb.add(V[2].x - p.x, V[2].y - p.y, V[2].z - p.z);
	bb.add(V[3].x - p.x, V[3].y - p.y, V[3].z - p.z);

	return bb;
} // bbox

/** @return oriented bounding box of body */
OBBox* GTETBody::updateOBB(bool in) const
{
	OBBox* bb = new OBBox();
	bb->P = GBody::position();
	Point V[4];
	if (mesh.isEmpty())
	{
		bb->infinite();
		return bb;
	}
	else
	{
		if (mesh.nvertices())
		{
			// use the vertices as a guess
			for (int i = 0; i < mesh.nvertices(); i++)
				bb->add(mesh.vertex(i));
		}
		else
			bb->infinite();
	}
	return bb;
} // obbox
/*---------------------------------------------------------------------------*/ // zxw20240926

/**
 * @param x,y	position
 * @return minimum distance to bodies conics
 */
double GBody::distance(const double x, const double y, const double z) const
{
	double min = INFINITE;
	for (int i=0; i<nQ(); i++) {
		double q = Abs(Q(i).adist(x,y,z));
		if (q < min)
			min = q;
	}
	return min;
} // distance

/* ---------------------------- GARBBody ------------------------------- */
// FIXME position, getWhat, Mesh...
void GARBBody::setWhat(double* what, char* err)
{
	fixWhat(what, nWhat(), SMALL);

	mesh.allocateVertices(8);

	// add vertices
	P.set(what[0], what[1], what[2]);
	for (int i=0; i<8; i++)
		mesh.vertex(i) = Point(what[i*3], what[i*3+1], what[i*3+2]);

	// copy faces
	bool same = true;
	for (int j=0,i=3*8; i<3*8+6; i++,j++) {
		int face = (int)what[i];
		if (face != faces[j]) {
			same = false;
			faces[j] = face;
		}

		int v = face%10 - 1;
		int zero = 0;
		if (v>=8) {
			if (err) sprintf(err,"Face %d contains an invalid vertex %d",(int)what[i],v);
			v = 0;
		}
		if (v==0) zero++;
		face /= 10;

		v = face%10 - 1;
		if (v>=8) {
			if (err) sprintf(err,"Face %d contains an invalid vertex %d",(int)what[i],v);
			v = 0;
		}
		if (v==0) zero++;
		face /= 10;

		v = face%10 - 1;
		if (v>=8) {
			if (err) sprintf(err,"Face %d contains an invalid vertex %d",(int)what[i],v);
			v = 0;
		}
		if (v==0) zero++;
		face /= 10;

		v = face%10 - 1;
		if (v>=8) {
			if (err) sprintf(err,"Face %d contains an invalid vertex %d",(int)what[i],v);
			v = 0;
		}
		if (v==0) zero++;

		if (zero>1) {
			if (err) sprintf(err,"Face %d contains too few vertices",(int)what[i]);
			continue;
		}
	}
	if (!same) {
		mesh.freeFaces();
		mesh.freeEdges();
	}
} // setWhat

/* --- getWhat --- */
int GARBBody::getWhat(double *what) const
{
	int w = 0;
	for (int i=0; i<8; i++) {
		what[w++] = mesh.vertex(i).x;
		what[w++] = mesh.vertex(i).y;
		what[w++] = mesh.vertex(i).z;
	}
	for (int i=0; i<6; i++)
		what[w++] = (double)faces[i];
	return 30;
} // getWhat

/** @return faceVertices */
bool GARBBody::faceVertices(int f, int v[4], Point* V[4])
{
	int face = faces[f];
	if (!face) {
		// The sizeof(v/V[0])*4 is needed since
		// sizeof(v/V) will return the size of int*/Point** not the size of the array
		memset(v,0,sizeof(v[0])*4);
		memset(V,0,sizeof(V[0])*4);
		return false;
	}
	for (int i=0; i<4; i++) {
		v[i] = face%10 - 1;
		if (v[i]>=8) v[i] = 0;
		if (v[i]>=0)
			V[i] = &mesh.vertex(v[i]);
		else
			V[i] = NULL;
		face /= 10;
	}
	return true;
} // faceVertices

/** createQuads */
void GARBBody::createQuads()
{
	// create faces and quadratics
	for (int i=0; i<6; i++) {
		int v[4];
		Point* V[4];
		faceVertices(i, v, V);
		if (!V[0] || !V[1] || !V[2]) continue;

		// create quads
		Vector N = (*V[1] - *V[0]) ^ (*V[2] - *V[0]);
		N.normalize();

		addQuad(-N.x, -N.y, -N.z, N*(*V[0]));

		if (V[3]!=NULL) {
			// check for co-planarity!
			if (!Eq0(_Q[_nQ-1](V[3]->x, V[3]->y, V[3]->z), SMALL))
				fprintf(stderr,"Face %d not planar\n",(int)faces[i]);
		}

		// Check remaining vertices, from which side they are
		bool first = true;
		bool positive = false;
		for (int j=0; j<8; j++) {
			if (j==v[0] || j==v[1] || j==v[2] || j==v[3]) continue;
			double q = _Q[_nQ-1](mesh.vertex(j).x, mesh.vertex(j).y, mesh.vertex(j).z);
			if (Eq0(q, SMALL))
				fprintf(stderr,"Vertex %d is coplanar with face %d\n",j+1,(int)faces[i]);
			if (first) {
				positive = (bool)(q>0.0);
				first = false;
			} else {
				if (positive != (bool)(q>0.0))
					fprintf(stderr,"Cannot resolve ARB orientation\n");
			}
		}
		// if remaining vertices are from same side with the normal
		if (positive)
			_Q[_nQ-1].negate();		// swap the sign of the quad
	}
} // createQuads

/* createMesh */
void GARBBody::createMesh()
{
	if (mesh.nedges()) return;

	for (int i=0; i<6; i++) {
		int v[4];
		Point* V[4];
		if (!faceVertices(i, v, V)) continue;
		if (!V[0] || !V[1] || !V[2]) continue;
		bool last = (bool)(v[3]>=0);
		mesh.add(v[0],v[1],v[2], true, true, !last);
		if (last)
			mesh.add(v[2],v[3],v[0], true, true, false);
	}
	mesh.calcBbox();
	mesh.process();
	if (!mesh.isClosed()) {
		fprintf(stderr,"ERROR open ARB a face is missing\n");
		return;
	}
	if (!mesh.isOrientable())
		mesh.makeOrientable();

	if (mesh.volume()<0.) mesh.flip();

	assert(mesh.isClosed());
	assert(mesh.isOrientable());
#if _DEBUG>2
	cout << endl;
	cout << "ARB Mesh " << endl;
	cout << mesh << endl;
	cout << " isClosed=" << mesh.isClosed();
	cout << " isOrientable=" << mesh.isOrientable() << endl;
	cout << "ARB volume=" << mesh.volume() << endl;
#endif
} // createMesh

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GARBBody::move(int item, const Point& r, const Vector& w)
{
	Point dR = r - P;
	for (int i=0; i<mesh.nvertices(); i++)
		mesh.vertex(i) += dR;

	GBody::move(item, r, w);
} // move

/** rotate body around an axis and current body position
 * @param angle	to rotate
 * @param axis	to rotate
 */
void GARBBody::rotate(const double angle, const Vector& axis)
{
	GBody::rotate(angle, axis);
	Matrix4	rot;			// duplicate rot matrix with the GBody::rotate
	rot.rotate(angle, axis);
	for (int i=1; i<8; i++) {
		Vector v = mesh.vertex(i) - P;
//		mesh.vertex(i) = rot.multVector(v) + P;
		mesh.vertex(i) = rot*v + P;
	}
} // rotate

/* ---------------------------- GQUABody ------------------------------- */
// FIXME position, getWhat, Mesh...
void GQUABody::setWhat(double *what, char *)
{
	fixWhat(what, nWhat(), SMALL1);
	q.set(what[0], what[1], what[2], what[3], what[4],
	      what[5], what[6], what[7], what[8], what[9]);

	Matrix3 M;
	q.matrix3(&M);
	Point v(-q.Cx/2.0, -q.Cy/2.0, -q.Cz/2.0);
	M.inverse();
	P = M * v;
} // setWhat

/* --- getWhat --- */
int GQUABody::getWhat(double *what) const
{
	what[0] = q.Cxx;
	what[1] = q.Cyy;
	what[2] = q.Czz;
	what[3] = q.Cxy;
	what[4] = q.Cxz;
	what[5] = q.Cyz;
	what[6] = q.Cx;
	what[7] = q.Cy;
	what[8] = q.Cz;
	what[9] = q.C;
	return 10;
} // getWhat

/** position
 * @param x,y,z	new position to set
 */
void GQUABody::position(const Point& r)
{
	GBody::position(r);

	Point tr = _hasMatrix? _invMatrix*r : r;
	tr -= SP;
	q = sq;
	q.translate(-tr);
} // position

/** save */
void GQUABody::save()
{
	GBody::save();
	sq = q;
} // save

/** restore */
void GQUABody::restore()
{
	q = sq;
	GBody::restore();
} // restore

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GQUABody::move(int item, const Point& r, const Vector& w)
{
	if (item==1) {
		Point tr = _hasMatrix? _invMatrix*r : r;
		tr -= SP;

		double a = sq.Cxx*Sqr(tr.x)
			 + sq.Cyy*Sqr(tr.y)
			 + sq.Czz*Sqr(tr.z);

		q = sq;
		if (!Eq0(a,TOOSMALL)) {
			double b = sq.Cxy*tr.x*tr.y
				 + sq.Cxz*tr.x*tr.z
				 + sq.Cyz*tr.y*tr.z
				 + sq.Cx*tr.x
				 + sq.Cy*tr.y
				 + sq.Cz*tr.z
				 + sq.C;
			double t = -b/a;
			q.Cxx *= t;
			q.Cyy *= t;
			q.Czz *= t;
		}
	} else
		GBody::move(item, r, w);
} // move

/** rotate body around an axis */
void GQUABody::rotate(const double angle, const Vector& axis)
{
	Matrix4	rot;
	rot.rotate(-angle, axis);
	//rot.fix();

	q = sq;
	q.translate(SP);
	q.transform(rot);
	q.translate(-SP);
} // rotate

/** @return bounding box of body */
BBox GQUABody::_bbox() const
{
	return Q(0).bbox();
} // bbox

/** @return oriented bounding box of body */
OBBox* GQUABody::updateOBB(bool in) const
{
	if (!in)
		return new OBBox(Q(0).obbox());
	else
		return new OBBox();
} // obbox

/* ---------------------------- GRCCBody ------------------------------- */
void GRCCBody::setWhat(double *what, char *err)
{
	fixWhat(what, nWhat(), SMALL2);
	P.set(what[0], what[1], what[2]);
	Z.set(what[3], what[4], what[5]);
	if (err && what[6]<TOOSMALL) strcpy(err,"Negative or Zero radius");
	R(Max(what[6],0.0));

	Po = P + Z;
	zlen = Z.normalize();

	findXYZ(err);
} // setWhat

/* --- getWhat --- */
int GRCCBody::getWhat(double *what) const
{
	what[0] = P.x;
	what[1] = P.y;
	what[2] = P.z;
	what[3] = zlen*Z.x;
	what[4] = zlen*Z.y;
	what[5] = zlen*Z.z;
	what[6] = R();
	return 7;
} // getWhat

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GRCCBody::move(int item, const Point& r, const Vector& w)
{
	Point tr = _hasMatrix? _invMatrix*r : r;
	tr -= SP;

	switch (item) {
		case 1:	// Find distance from Z-axis
			R(Sqrt(tr.length2() - Sqr(tr*Z)));
			break;

		case 2: {	// Move base
				double zp = tr*SZ;
				P = SP + zp*SZ;
				zlen = szlen - zp;
				if (zlen<0.0) {
					Z = -SZ;
					zlen = -zlen;
				} else
					Z = SZ;
			}
			break;

		case 3:		// Move apex
			zlen = tr*SZ;
			if (zlen<0.0) {
				Z = -SZ;
				zlen = -zlen;
			} else
				Z = SZ;
			Po = SP + zlen*Z;
			break;

		default:
			GBody::move(item, r, w);
	}
} // move

/** createQuads */
void GRCCBody::createQuads()
{
	// Create the cylinder around Z then transform
	//  Cxx Cyy Czz  Cxy  Cxz  Cyz  Cx   Cy   Cz    C
	addQuad(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-Sqr(R()));
	makeConeQuads();
} // createQuads

/** @return bounding box of body */
BBox GRCCBody::_bbox() const
{
	BBox bb;
	Point p = position();
	Point po = p + (vectorZ()*zlen);

	bb.add(p);
	bb.add(po);

	bbox_addRotatedEllipse(bb, vectorX(), vectorY(), p,  R(), R());
	bbox_addRotatedEllipse(bb, vectorX(), vectorY(), po, R(), R());

	return bb;
} // bbox

/** @return oriented bounding box of body */
OBBox* GRCCBody::updateOBB(bool in) const
{
	OBBox* bb = new OBBox;

	bb->P = GBody::position();
	bb->X = vectorX();
	bb->Y = vectorY();
	bb->Z = vectorZ();

	if (!in) {
		bb->low()  = Point(-R(), -R(), 0.0);
		bb->high() = Point( R(),  R(), zlen);
	} else {
		double sc = R()/Sqrt2;

		bb->low()  = Point(-sc, -sc, 0.0);
		bb->high() = Point( sc,  sc, zlen);
	}

	return bb;
} // obbox

/* ---------------------------- GRECBody ------------------------------- */
void GRECBody::setWhat(double *what, char *err)
{
	fixWhat(what, nWhat(), SMALL2);
	P.set(what[0], what[ 1], what[ 2]);
	Z.set(what[3], what[ 4], what[ 5]);
	X.set(what[6], what[ 7], what[ 8]);
	Y.set(what[9], what[10], what[11]);

	Po = P + Z;

	// correct for perpendicularity
	xlen = X.normalize();
	ylen = Y.normalize();
	zlen = Z.normalize();

	if (err && xlen<TOOSMALL) strcpy(err,"X-length is zero");
	if (err && ylen<TOOSMALL) strcpy(err,"Y-length is zero");

	checkOrthogonal(err);
} // setWhat

/* --- getWhat --- */
int GRECBody::getWhat(double *what) const
{
	what[ 0] = P.x;
	what[ 1] = P.y;
	what[ 2] = P.z;
	what[ 3] = zlen*Z.x;
	what[ 4] = zlen*Z.y;
	what[ 5] = zlen*Z.z;
	what[ 6] = xlen*X.x;
	what[ 7] = xlen*X.y;
	what[ 8] = xlen*X.z;
	what[ 9] = ylen*Y.x;
	what[10] = ylen*Y.y;
	what[11] = ylen*Y.z;
	return 12;
} // getWhat

/**  @return closest point, edge, face */
int GRECBody::closest(const Point& r, const double d, const Vector& w) const
{
	int close = GBody::closest(r,d,w);
	if (close == 1) {
		double px = Abs(Y*(r-P));
		double py = Abs(X*(r-P));
		if (px>py) return 4;
	}
	return close;
} // closest

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GRECBody::move(int item, const Point& r, const Vector& w)
{
	Point tr = _hasMatrix? _invMatrix*r : r;
	tr -= SP;
	switch (item) {
		case 1:	{ // X-radius
				// Find a solution where the surface cross the r
				Point rp = P + (tr*Z)*Z;	// projection on Z axis
				double x = (r-rp)*X;		// x,y components on X,Y axis
				double y = (r-rp)*Y;
				// equation
				//  x^2    y^2
				//  ---  + ---  = 1
				//  xl^2   yl^2
				double d = 1.0 - Sqr(y)/Sqr(sylen);
				if (d > SMALL)
					xlen = sqrt(Sqr(x) / d);
				else	// On failure find distance from Y-axis
					xlen = Sqrt(tr.length2() - Sqr(tr*Y));
			}
			break;

		case 2: {	// Move base
				double zp = tr*SZ;
				P = SP + zp*SZ;
				zlen = szlen - zp;
				if (zlen<0.0) {
					Z = -SZ;
					zlen = -zlen;
				} else
					Z = SZ;
			}
			break;

		case 3:		// Move apex
			zlen = tr*SZ;
			if (zlen<0.0) {
				Z = -SZ;
				zlen = -zlen;
			} else
				Z = SZ;
			Po = P + zlen*Z;
			break;

		case 4:	{ // X-radius
				// Find distance from Z-axis
				Point rp = P + (tr*Z)*Z;	// projection on Z axis
				double x = (r-rp)*X;	// x,y components on X,Y axis
				double y = (r-rp)*Y;
				double d = 1.0 - Sqr(x)/Sqr(sxlen);
				if (d > SMALL)
					ylen = Sqrt(Sqr(y) / d);
				else	// On failure find distance from X-axis
					ylen = Sqrt(tr.length2() - Sqr(tr*X));
			}
			break;

		default:
			GBody::move(item, r, w);
	}
} // move

/** createQuads */
void GRECBody::createQuads()
{
	if (ISZERO(xlen) || ISZERO(ylen)) return;

	// Create the cylinder around Z then transform
	//  Cxx Cyy Czz  Cxy  Cxz  Cyz  Cx   Cy   Cz    C
	addQuad( 1.0/Sqr(xlen),  1.0/Sqr(ylen), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-1.0);
	makeConeQuads();
} // createQuads

/** @return bounding box of body */
BBox GRECBody::_bbox() const
{
	BBox bb;
	Point p = position();
	Point po = p + (vectorZ()*zlen);

	bb.add(p);
	bb.add(po);

	bbox_addRotatedEllipse(bb, vectorX(), vectorY(), p,  xlen, ylen);
	bbox_addRotatedEllipse(bb, vectorX(), vectorY(), po, xlen, ylen);

	return bb;
} // bbox

/** @return oriented bounding box of body */
OBBox* GRECBody::updateOBB(bool in) const
{
	OBBox* bb = new OBBox();

	bb->P = GBody::position();
	bb->X = vectorX();
	bb->Y = vectorY();
	bb->Z = vectorZ();

	if (!in) {
		bb->low()  = Point(-xlen, -ylen, 0.0);
		bb->high() = Point( xlen,  ylen, zlen);
	} else {
		double sc = 1.0/Sqrt2;
		bb->low()  = Point(-xlen*sc, -ylen*sc, 0.0);
		bb->high() = Point( xlen*sc,  ylen*sc, zlen);
	}

	return bb;
} // obbox

/* ---------------------------- GTRCBody ------------------------------- */
void GTRCBody::setWhat(double *what, char *err)
{
	fixWhat(what, nWhat(), SMALL2);
	P.set(what[0], what[1], what[2]);
	Z.set(what[3], what[4], what[5]);
	if (err && what[6]<-TOOSMALL) strcpy(err,"Negative base radius");
	if (err && what[7]<-TOOSMALL) strcpy(err,"Negative apex radius");
	Rb(Max(what[6],0.0));	// Rb base radius
	Ra(Max(what[7],0.0));	// Ra apex radius

	Po   = P + Z;
	h2   = Z.length2();	// to keep precision
	zlen = Z.normalize();
	if (err && zlen<TOOSMALL) strcpy(err,"Zero height");

	findXYZ(err);
} // setWhat

/* --- getWhat --- */
int GTRCBody::getWhat(double *what) const
{
	what[0] = P.x;
	what[1] = P.y;
	what[2] = P.z;
	what[3] = zlen*Z.x;
	what[4] = zlen*Z.y;
	what[5] = zlen*Z.z;
	what[6] = Rb();
	what[7] = Ra();
	return 8;
} // getWhat

/**  @return closest point, edge, face */
int GTRCBody::closest(const Point& r, const double d, const Vector& w) const
{
	int close = GBody::closest(r,d,w);
	if (close == 1) {
		// Check which is closer, base or apex?
		if (Abs(Q(1)(r)) > Abs(Q(2)(r)))
			return 4;
	}
	return close;
} // closest

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GTRCBody::move(int item, const Point& r, const Vector& w)
{
	Point tr = _hasMatrix? _invMatrix*r : r;
	tr -= SP;
	switch (item) {
		case 1:	// Find distance from Z-axis
			Rb(Sqrt(tr.length2() - Sqr(tr*Z)));	// base radius
			break;

		case 2: {	// Move base
				double zp = tr*SZ;
				P = SP + zp*SZ;
				zlen = szlen - zp;
				if (zlen<0.0) {
					Z = -SZ;
					zlen = -zlen;
				} else
					Z = SZ;
				h2 = Sqr(zlen);
			}
			break;

		case 3:		// Move apex
			zlen = tr*SZ;
			if (zlen<0.0) {
				Z = -SZ;
				zlen = -zlen;
			} else
				Z = SZ;
			h2 = Sqr(zlen);
			Po = P + zlen*Z;
			break;

		case 4:	// Change appex radius
			Ra(Sqrt(tr.length2() - Sqr(tr*Z)));
			break;

		default:
			GBody::move(item, r, w);
			h2 = Sqr(zlen);
	}
} // move

/** createQuads */
void GTRCBody::createQuads()
{
	double s = Ra()-Rb();	// slope=(Ra-Rb)/zlen but avoid division
				// that generates roundup errors
				//
	if (ISZERO(zlen)) return;

	// Create the cone around Z then transform
	//  Cxx Cyy Czz  Cxy  Cxz  Cyz  Cx   Cy   Cz    C
	addQuad(1.0, 1.0, -s*s/h2, 0.0, 0.0, 0.0, 0.0, 0.0,-2.0*s*Rb()/zlen, -Rb()*Rb());
	makeConeQuads();
} // createQuads

/** @return bounding box of body */
BBox GTRCBody::_bbox() const
{
	BBox bb;
	Point p  = position();
	Point po = p + (vectorZ()*zlen);

	bb.add(p);
/** @return bounding box of quadric (in case of ellipsoid) */
	bb.add(po);

	bbox_addRotatedEllipse(bb, vectorX(), vectorY(), p,  Rb(), Rb());
	bbox_addRotatedEllipse(bb, vectorX(), vectorY(), po, Ra(), Ra());

	return bb;
} // bbox

/** @return oriented bounding box of body */
OBBox* GTRCBody::updateOBB(bool in) const
{
	OBBox* bb = new OBBox();

	bb->P = GBody::position();
	bb->X = vectorX();
	bb->Y = vectorY();
	bb->Z = vectorZ();

	if (!in) {
		double R = Max(Ra(), Rb());
		bb->low()  = Point(-R, -R, 0.0);
		bb->high() = Point( R,  R, zlen);
	} else {
		double sc, z;
		if (Rb() > Ra()) {
			sc = Rb() / 2.0;

			if (Ra() >= sc) {
				sc = Ra();
				z  = zlen;
			} else
				z  = sc*zlen/(Rb()-Ra());
			sc /= Sqrt2;

			bb->low()  = Point(-sc, -sc, 0.0);
			bb->high() = Point( sc,  sc, z);
		} else {
			sc = Ra() / 2.0;
			if (Rb() >= sc) {
				sc = Rb();
				z  = 0.0;
			} else
				z = zlen - sc*zlen/(Ra()-Rb());

			sc /= Sqrt2;
			bb->low()  = Point(-sc, -sc, z);
			bb->high() = Point( sc,  sc, zlen);
		}
	}

	return bb;
} // obbox

/* -------------------------- GInfEllCylBody --------------------------- */
void GInfEllCylBody::setWhat(double *what, char *err)
{
	fixWhat(what, nWhat(), SMALL2);
	if (err && what[2]<TOOSMALL)   strcpy(err,"Negative radius");
	if (err && what[3]<-TOOSMALL)  strcpy(err,"Negative radius");

	Rx(Max(what[2],0.0));
	Ry(Max(what[3],0.0));
	zlen = infinite;

	if (type()==XCCbody || type()==XECbody) {
		P.set(0.0, what[0], what[1]);
		if (type()==XCCbody) ylen = xlen;
		Z = Vector::Xo;
		X = Vector::Yo;
		Y = Vector::Zo;
	} else
	if (type()==YCCbody || type()==YECbody) {
		P.set(what[1], 0.0, what[0]);
		if (type()==YCCbody) ylen = xlen;
		Z = Vector::Yo;
		X = Vector::Zo;
		Y = Vector::Xo;
	} else

	if (type()==ZCCbody || type()==ZECbody) {
		P.set(what[0], what[1], 0.0);
		if (type()==ZCCbody) ylen = xlen;
		X = Vector::Xo;
		Y = Vector::Yo;
		Z = Vector::Zo;
	}
} // setWhat

/* --- getWhat --- */
int GInfEllCylBody::getWhat(double *what) const
{
	what[2] = Rx();
	what[3] = Ry();		// For [XYZ]EC only

	switch (type()) {
		case XCCbody:
		case XECbody:
			what[0] = P.y;
			what[1] = P.z;
			if (type()==XECbody) return 4;
			return 3;

		case YCCbody:
		case YECbody:
			what[0] = P.z;
			what[1] = P.x;
			if (type()==YECbody) return 4;
			return 3;

		case ZCCbody:
		case ZECbody:
			what[0] = P.x;
			what[1] = P.y;
			if (type()==ZECbody) return 4;
			return 3;

		default: return -1;
	}
} // getWhat

/**  @return closest point, edge, face */
int GInfEllCylBody::closest(const Point& r, const double d, const Vector& w) const
{
	int close = GBody::closest(r,d,w);

	switch (type()) {
		case XECbody:
			if (Abs(r.y-P.y) < Abs(r.z-P.z))
				return 2;
			return close;
		case YECbody:
			if (Abs(r.z-P.z) < Abs(r.x-P.x))
				return 2;
			return close;
		case ZECbody:
			if (Abs(r.x-P.x) < Abs(r.y-P.y))
				return 2;
			return close;
		default:
			return close;
	}
} // closest

/** move body item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 * @param w	vertical axis of the viewport
 */
void GInfEllCylBody::move(int item, const Point& r, const Vector& w)
{
	Point tr = _hasMatrix? _invMatrix*r : r;
	tr -= SP;

	double x, y;
	bool ell = true;

	switch (type()) {
		case XCCbody:
			ell = false;
		case XECbody:
			x = tr.y;
			y = tr.z;
			break;

		case YCCbody:
			ell = false;
		case YECbody:
			x = tr.z;
			y = tr.x;
			break;

		case ZCCbody:
			ell = false;
		case ZECbody:
			x = tr.x;
			y = tr.y;
			break;

		default:
			x = y = 0.0;	// Never comes here but keeps compiler happy
			assert(0);
	}

	if (item==1) {
		if (ell) {
			double d = 1.0 - Sqr(y)/Sqr(ylen);
			if (d > SMALL)
				Rx(Sqrt(Sqr(x) / d));
			else	// On failure find distance from Y-axis
				Rx(x);
		} else
			R(Sqrt(Sqr(x) + Sqr(y)));
	} else
	if (item==2) {
		double d = 1.0 - Sqr(x)/Sqr(Rx());
		if (d > SMALL)
			ylen = Sqrt(Sqr(y) / d);
		else
			ylen = y;
	} else
		GBody::move(item, r, w);
} // move

/** rotate body around an axis */
void GInfEllCylBody::rotate(const double angle, const Vector& axis)
{
	GBody::rotate(angle, axis);

	/* find closest axis to align to */
	double ax = Abs(Z.x);
	double ay = Abs(Z.y);
	double az = Abs(Z.z);

	if (az>ax && az>ay) {
		Z = Vector::Zo;
		if (InRange(XECbody, type(), ZECbody))
			_type = ZECbody;
		else
			_type = ZCCbody;
	} else
	if (ay>ax && ay>az) {
		Z = Vector::Yo;
		if (InRange(XECbody, type(), ZECbody))
			_type = YECbody;
		else
			_type = YCCbody;
	} else {
		Z = Vector::Zo;
		if (InRange(XECbody, type(), ZECbody))
			_type = XECbody;
		else
			_type = XCCbody;
	}

	findXYZ(NULL);
} // rotate

/** createQuads */
void GInfEllCylBody::createQuads()
{
	double a = 1.0, b = 1.0, c = 1.0;
	double r2 = Sqr(xlen);

	if (type()==XCCbody || type()==XECbody) {
		a = 0.0;
		if (type()==XECbody) {
			c   = r2;
			b   = Sqr(ylen);
			r2 *= b;
		}
	} else
	if (type()==YCCbody || type()==YECbody) {
		b = 0.0;
		if (type()==YECbody) {
			a   = r2;
			c   = Sqr(ylen);
			r2 *= c;
		}
	} else
	if (type()==ZCCbody || type()==ZECbody) {
		c = 0.0;
		if (type()==ZECbody) {
			b   = r2;
			a   = Sqr(ylen);
			r2 *= a;
		}
	}

	//  Cxx Cyy Czz Cxy  Cxz  Cyz
	addQuad( a,  b,  c, 0.0, 0.0, 0.0,
	//            Cx            Cy            Cz
		-2.0*a*P.x,  -2.0*b*P.y,  -2.0*c*P.z,
		 a*Sqr(P.x) + b*Sqr(P.y) + c*Sqr(P.z) - r2);

} // createQuads

/**
 * @param body	other body to check location against
 * @return location of this body with respect to another body
 */
Location GInfEllCylBody::_locationWrt(const GBody *body) const
{
//	if (body->type() == RCCbody) {
//	} else
	if (type() == XCCbody || type() == YCCbody || type() == ZCCbody) {
		if (body->type() == XCCbody || body->type() == YCCbody ||
		    body->type() == ZCCbody) {
			Vector thisZ = vectorZ();
			Vector bodyZ = body->vectorZ();

			// Minimum distance between the two axes
			double dist = lineLineDistance(position(), thisZ, body->position(), bodyZ);

			// If they are parallel.
			// FIXME could be that after a rotation thisZ = +/- bodyZ
			if (type() == body->type()) {
				if (dist + ((const GInfEllCylBody*)body)->R() <  R())
					return LOCATION_BinA;
				else
				if (dist - ((const GInfEllCylBody*)body)->R() < -R())
					return LOCATION_AinB;
				else
				if (dist - ((const GInfEllCylBody*)body)->R() >  R())
					return LOCATION_OUTSIDE;
				else
					return LOCATION_OVERLAP;
			} else {
				// any other skew orientation (after rotation)
				if (dist - ((const GInfEllCylBody*)body)->R() >  R())
					return LOCATION_OUTSIDE;
				else
					return LOCATION_OVERLAP;
			}
		}
	}

	return GBody::_locationWrt(body);
} // _locationWrt

/** @return bounding box of body */
BBox GInfEllCylBody::_bbox() const
{
	BBox bb;
	switch(type()) {
		case XCCbody:
			bb.low().set(-INFINITE, P.y-R(), P.z-R());
			bb.high().set(INFINITE, P.y+R(), P.z+R());
			break;
		case YCCbody:
			bb.low().set( P.x-R(), -INFINITE, P.z-R());
			bb.high().set(P.x+R(),  INFINITE, P.z+R());
			break;
		case ZCCbody:
			bb.low().set( P.x-R(), P.y-R(), -INFINITE);
			bb.high().set(P.x+R(), P.y+R(),  INFINITE);
			break;

		case XECbody:
			bb.low().set(-INFINITE, P.y-Rx(), P.z-Ry());
			bb.high().set(INFINITE, P.y+Rx(), P.z+Ry());
			break;
		case YECbody:
			bb.low().set( P.x-Ry(), -INFINITE, P.z-Rx());
			bb.high().set(P.x+Ry(),  INFINITE, P.z+Rx());
			break;
		case ZECbody:
			bb.low().set( P.x-Rx(), P.y-Ry(), -INFINITE);
			bb.high().set(P.x+Rx(), P.y+Ry(),  INFINITE);
			break;

		default:
			assert(0);
	}
	if (_hasMatrix) bb.transform(_matrix);

	return bb;
} // bbox

/** @return oriented bounding box of body */
OBBox* GInfEllCylBody::updateOBB(bool in) const
{
	OBBox* bb = new OBBox();
	bb->infinite();

	bb->P = GBody::position();
	bb->X = vectorX();
	bb->Y = vectorY();
	bb->Z = vectorZ();

	if (!in) {
		bb->low().x  = -Rx();
		bb->low().y  = -Ry();
		bb->high().x =  Rx();
		bb->high().y =  Ry();
	} else {
		double sc = 1.0/Sqrt2;

		bb->low().x  = -Rx()*sc;
		bb->low().y  = -Ry()*sc;
		bb->high().x =  Rx()*sc;
		bb->high().y =  Ry()*sc;
	}

	return bb;
} // obbox

/* -------------------------- GTorusBody --------------------------- */
void GTorusBody::setWhat(double *what, char *err)
{
	fixWhat(what, nWhat(), SMALL2);

	P.set(what[0], what[1], what[2]);
	if (err && what[3]<TOOSMALL)  strcpy(err,"Too small or negative radius");
	if (err && what[4]<TOOSMALL)  strcpy(err,"Too small or negative radius");
	if (err && what[5]<TOOSMALL)  strcpy(err,"Too small or negative radius");
	a = what[3];
	b = what[4];
	c = what[5];

	switch (type()) {
		case TRXbody:
			Z = Vector::Xo;
			X = Vector::Yo;
			Y = Vector::Zo;
			break;

		case TRYbody:
			Z = Vector::Yo;
			X = Vector::Zo;
			Y = Vector::Xo;
			break;

		case TRZbody:
			Z = Vector::Zo;
			X = Vector::Xo;
			Y = Vector::Yo;
			break;
		default:
			;
	}
} // setWhat

/* --- getWhat --- */
int GTorusBody::getWhat(double *what) const
{
	what[0] = P.x;
	what[1] = P.y;
	what[2] = P.z;
	what[3] = a;
	what[4] = b;
	what[5] = c;
	return 6;
} // getWhat

/** @return bounding box of body */
BBox GTorusBody::_bbox() const
{
	BBox bb;
	Point p = GBody::position();
	bb.add(Point(p - (a+c)*(vectorX()+vectorY()) - b*vectorZ()));
	bb.add(Point(p + (a+c)*(vectorX()+vectorY()) + b*vectorZ()));
	return bb;
} // bbox

/** createQuads */
void GTorusBody::createQuads()
{
	// XXX XXX XXX XXX XXX XXX
	// FIXME FOR TESTING create two cylinders and two planes
	// XXX XXX XXX XXX XXX XXX
	//
	// Create the cylinder around Z then transform
	//  Cxx Cyy Czz  Cxy  Cxz  Cyz  Cx   Cy   Cz    C
	addQuad(-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Sqr(a-c));
	addQuad( 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-Sqr(a+c));

	Matrix4 M, Minv, T;
	Minv.make(X,Y,Z);
	Minv.inverse();
	T.translate(P);
	M.multiply(T, Minv);
	Minv.inverse(M);

	// transform first quad
	_Q[0].transform(Minv);
	_Q[0].normalize();

	_Q[1].transform(Minv);
	_Q[1].normalize();

	// two planes for the caps
	//        Cx     Cy     Cz      C
	addQuad(-Z.x,-Z.y,-Z.z,  Z*P - b);	// base
	addQuad( Z.x, Z.y, Z.z,-(Z*P + b));	// apex
} // createQuads

/* createMesh */
void GTorusBody::createMesh()
{
	bool first = mesh.isEmpty();
	if (first) mesh.allocateVertices(N_TORUS_LAT * N_TORUS_LON);

	const double stepf = PI2 / (double)N_TORUS_LAT;
	const double step0 = PI2 / (double)N_TORUS_LON;

	int k = 0;
	double angf = 0.0;
	for (int j=0; j<N_TORUS_LAT; j++, angf += stepf) {
		double sf, cf;
		bsincos(angf, &sf, &cf);
		double ang0 = 0.0;
		for (int i=0; i<N_TORUS_LON; i++, ang0 += step0) {
			double s0, c0;
			bsincos(ang0, &s0, &c0);
			mesh.vertex(k++) = P +
					(a + c*c0)*cf * X +
					(a + c*c0)*sf * Y +
					b*s0*Z;
		}
	}

	mesh.calcBbox();

	if (first) {
		k = 0;
		for (int j=0; j<N_TORUS_LAT; j++) {
			for (int i=0; i<N_TORUS_LON; i++, k++) {
				int n = (j<N_TORUS_LAT-1)? k+N_TORUS_LON : i;

				int m = (i<N_TORUS_LON-1)? n+1 : n+1-N_TORUS_LON;
				mesh.add(k, n, m, true, true, false);

				n = (i<N_TORUS_LON-1)? k+1 : k+1-N_TORUS_LON;
				mesh.add(k, m, n, false, true, true);
			}
		}

		mesh.process();
		assert(mesh.isClosed());
		assert(mesh.isOrientable());
//#if _DEBUG>2
		cout << endl;
		cout << "TRXYZ Mesh ";
		cout << " isClosed=" << mesh.isClosed();
		cout << " isOrientable=" << mesh.isOrientable() << endl;
		cout << "TRZ volume=" << mesh.volume() << endl;
//#endif
	}
} // createMesh

/** operator << */
ostream& operator << (ostream& s, const GBody& body)
{
	double what[30];
	int n = body.getWhat(what);
	s << body.typeStr() << " " << body.name();
	for (int i=0; i<n; i++)
		s << " " << what[i];
	return s;
} /* operator << */
