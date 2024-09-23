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

#ifndef __GBODY_H
#define __GBODY_H

#include <math.h>
#include <assert.h>
#include <iostream>

#include "os.h"
#include "geo.h"
#include "bbox.h"
#include "edge.h"
#include "mesh.h"
#include "quad.h"
#include "array.h"
#include "bstring.h"
#include "matrix4.h"

class OBBox;
class GZone;

enum Location {
	LOCATION_OUTSIDE = 0,	// The viewport is outside the body
	LOCATION_INSIDE  = 1,	// The viewport is inside the body
	LOCATION_OVERLAP = 2,	// The viewport is overlapping the body

	LOCATION_AinB    = 3,	// A body is inside B	(only from location Wrt)
	LOCATION_BinA    = 4,	// B body is inside A	(only from location Wrt)
};

// WARNING DO NOT change the order the following depend on the order
// @see _bodyTypeStr, _locationWrt
enum BodyType {
	YZPbody,
	XZPbody,
	XYPbody,
	PLAbody,

	RPPbody,
	BOXbody,
	WEDbody,
	RAWbody,

	SPHbody,
	ELLbody,

	RCCbody,
	RECbody,
	TRCbody,

	XCCbody,
	YCCbody,
	ZCCbody,

	XECbody,
	YECbody,
	ZECbody,

	ARBbody,
	QUAbody,

	TRXbody,
	TRYbody,
	TRZbody,
	/* end of normal bodies */

	NULLbody,
	PLUSbody,
	MINUSbody,
	UNIONbody,
	UNIVERSEbody,
	LEFTbody,
	RIGHTbody,

	ERRbody
};

#define BODYQUADS	 6
#define BODYCONICS	 6
#define MAXWHATS	30

class GOPRBody;

/* ============================== GBody =============================== */
/** GBody geometrical body class */
class GBody {
protected:
	int	_id;		// index id in the bodies array
	char16	_name;		// name of body
	BodyType _type;		// FLUKA type of body
	dword	_color;		// body color
	int	_width;		// line width to use
	int	_generation;	// Generation number

	int	_nQ;		// Number of quads
	Quad	_Q[BODYQUADS];	// Surface quadratic equations

	bool    _hasMatrix;	// if rotation matrix is set
	Matrix4	_matrix;	// rotation matrix
	Matrix4 _invMatrix;	// inverse rotation matrix

public:
	int	 show;		// display body (projection segments)
	Mesh	mesh;		// mesh of body

// FIXME could be sorted list with pointer value for faster searching!!!
	Array<GZone*>	zones;  // zones that refer to this body

protected:
	Point	P, Po;		// position, outer position
	Vector	X, Y, Z;	// unit vectors along main axes
	double	xlen, ylen, zlen;// dimensions of X Y Z or radius

	bool	_userBboxFlag;	// user has defined a bbox
	BBox	_userBbox;	// user defined bbox
mutable OBBox	*_cached_obbox[2];

	Point	SP, SPo;	// Save variables
	Vector	SX, SY, SZ;	// to be restored later
	double	sxlen, sylen, szlen;
	int	sshow;		// show saved flag

protected:
static	double	infinite;	// infinite dimension to use

public:
static	const char*	_typeStr[];
static	GOPRBody	tnull;		// empty
static	GOPRBody	tplus;		// +
static	GOPRBody	tminus;		// -
static	GOPRBody	tunion;		// |
static	GOPRBody	tuniverse;	// @
static	GOPRBody	tleft;		// (
static	GOPRBody	tright;		// )

public:
	GBody(const char *aname, const BodyType atype);
	virtual ~GBody();

	// name
	void	name(const char *aname);
const	char*	name()		const	{ return _name; }

	// Id
	void	id(int i)		{ _id = i; }
	int	id()		const	{ return _id; }

	// Color
	void	color(dword c)		{ _color = c & FLAG_COLORMASK; }
	dword	color()		const	{ return _color; }

	// line width
	void	lineWidth(int w)	{ _width = w; }
	int	lineWidth()	const	{ return _width; }

	// return name hash
	dword	hash()		const	{ return hash_djb2(name()); }

	// Type
	BodyType type()		const	{ return _type; }
const	char*	typeStr()	const	{ return _typeStr[_type]; }
	bool	isOperator()	const	{ return _type >= NULLbody; }

	// Generation
	int	generation()	const	{ return _generation; }
	void	nextGeneration(int newGen) {
			_generation = newGen;
			clearOBB();
		}

	// Whats
virtual	void	setWhat(double *, char *err) = 0;
virtual	int	getWhat(double *) const { return 0; }
virtual	int	nWhat() const { return 0; }

	// Save/Restore position
virtual	void	save();
virtual	void	restore();
	Point	savedPosition()	const;

	// Zone reference
	void	addZone(GZone *zone)	{ zones.add(zone);   }
	void	delZone(GZone *zone)	{ zones.erase(zone); }
	bool	hasZone(GZone *zone)	{ return zones.find(zone)>=0; }

	// Position
	Point	position() const;
virtual	void	position(const Point& r);

	Vector	vectorX() const;
	Vector	vectorY() const;
	Vector	vectorZ() const;

	double	vectorXlen()	const	{ return xlen; }
	double	vectorYlen()	const	{ return ylen; }
	double	vectorZlen()	const	{ return zlen; }

	// Nodes
virtual	int	nodes()		const	{ return 1; }
virtual Point	node(const int)	const;

	// Selection
virtual	int	closest(const Point& r, const double d, const Vector& w) const;

	// Matrix and transformations
	void	clearMatrix()		{ _hasMatrix = false; }
	void	matrix(const Matrix4& M);
const	Matrix4& matrix()	const	{ return _matrix; }

	void	transform();
	bool	hasMatrix()	 const	{ return _hasMatrix; }
virtual void	move(int item, const Point& r, const Vector& w);
virtual	void	rotate(const double angle, const Vector& axis);

	// Quadrics
	int	nQ()		const	{ return _nQ; }
const	Quad&	Q(int n)	const	{ assert(n>=0 && n<_nQ);return _Q[n]; }
virtual	void	resetQuads()		{ _nQ = 0; }
virtual	void	createQuads()		{ }

	// bounding box
	void	bbox(const BBox& bb);
virtual	BBox	bbox() const {
			if (_userBboxFlag) return _userBbox;
			return _bbox();
		}
virtual	BBox	_bbox() const;

private:
virtual	OBBox* updateOBB(bool in = false) const;
	void clearOBB();
public:
	OBBox*	obbox(bool in) const;

	// Location
virtual	Location locationWrt(const GBody* body) const;
	Location bbLocationWrt(const GBody* body) const;

	// Mesh
virtual	void	createMesh() {};

	void	create() {
			resetQuads();
			createQuads();
			createMesh();
		}

	// Dimension labels and show flags
	// @return NULL not to show anything otherwise label
virtual	char*	showX()	const { return "X"; }
virtual char*	showY()	const { return "Y"; }
virtual char*	showZ()	const { return "Z"; }

	double	distance(const double x, const double y, const double z) const;
	// ray intersection (non-cached)
	bool	intersectRay(const double x, const double y, const double z,
			     const double dx, const double dy, const double dz,
			     double *tmin, double *tmax) const;
	bool	intersectRay(const Point& pos, const Vector& dir,
			     double *tmin, double *tmax) const {
			return intersectRay(pos.x, pos.y, pos.z,
					dir.x, dir.y, dir.z,
					tmin, tmax);
		}
	Vector	normal(const Point& r) const;
	bool	inside(	const double  x, const double  y, const double  z,
			const double dx, const double dy, const double dz,
			const int ignore1=-1, const int ignore2=-1, const int ignore3=-1) const;

	bool	inside( const double x, const double y, const double z,
			const Quad *ignore_a, const Quad *ignore_b, const Quad *ignore_c) const;
	bool	inside2D(const double  x, const double  y, const double  z,
			 const double dx, const double dy, const double dz,
			 const double acc[],
#ifdef EXPERIMENTAL
			 bool *withinAcc=NULL,
#endif
			 const int ignore1=-1, const int ignore2=-1) const;
	bool	outside(const double x, const double y, const double z,
			const Quad *ignore_a, const Quad *ignore_b, const Quad *ignore_c) const;

	size_t	memory() const;
protected:
	void	checkOrthogonal(char *err);

	// add a new quad
	void	addQuad(double Cxx, double Cyy, double Czz,
			double Cxy, double Cxz, double Cyz,
			double Cx,  double Cy,  double Cz,
			double C);
	void	addQuad(double Cx,  double Cy,  double Cz, double C);
	void	addQuad(const Quad& q)
			{ addQuad(q.Cxx, q.Cyy, q.Czz,
				  q.Cxy, q.Cxz, q.Cyz,
				  q.Cx,  q.Cy,  q.Cz, q.C); }
	void	makeConeQuads();
	void	findXYZ(char *err);

	// mesh functions
	void	createEllConeMesh(const double Rbx, const double Rby,
				  const double Rax, const double Ray,
				  const bool isinfinite=false);
	void	createEllipsoidMesh();

virtual	Location _locationWrt(const GBody *body) const;

	void	fixWhat(double *what, int n, const double eps);
//	void	setEdge(const int i, const int a, const int b) { edge(i).set(&vertex(a), &vertex(b)); }

public:
static	GBody*	newBody(const char *aname, const char *atype);
static	GBody*	newBody(const char *name, const BodyType type);
//static	int	compare(const GBody* a, const GBody* b)
//			{ return strcmp(a->name(),b->name()); }
static	const char* typeStr(BodyType t)		{ return _typeStr[t]; }
}; // GBody

/* --- GOPRBody ---
 * Operator body null, +, -, |, (, )
 */
class GOPRBody : public GBody {
public:
	GOPRBody(const char *aname, const BodyType atype)
		: GBody(aname, atype) {}

	virtual	void	setWhat(double *, char *) {};
	virtual	void	position(const Point&) {};
	virtual	int	nodes()		const	{ throw(0); }
	virtual void	move(int, const Point&, const Vector&) {};
	virtual	void	rotate(const double, const Vector&) {};
	virtual	void	createQuads() {};

	// bounding box
	virtual	BBox	_bbox()		const	{ throw(0); };
	virtual	OBBox * updateOBB(bool)	const	{ throw(0); };
protected:
	virtual	Location _locationWrt(const GBody *) const { throw(0); };
}; // GOPRBody

/* --- GERRBody ---
 * Error body
 */
class GERRBody : public GBody {
public:
	GERRBody(const char *aname, const BodyType atype)
		: GBody(aname, atype) {}

	virtual	void	setWhat(double *, char *) {};
	virtual	void	position(const Point&) {};
	virtual	int	nodes()		const	{ throw(0); }
	virtual void	move(int, const Point&, const Vector&) {};
	virtual	void	rotate(const double, const Vector&) {};
	virtual	void	createQuads() {};

	// bounding box
	virtual	BBox	_bbox()		const	{ throw(0); };
	virtual	OBBox * updateOBB(bool)	const	{ throw(0); };
protected:
	virtual	Location _locationWrt(const GBody *) const { throw(0); };
}; // GERRBody

/* --- GPLABody --- */
class GPLABody : public GBody {
private:
	int	_sign20;
public:
	GPLABody(const char *aname, const BodyType atype)
		: GBody(aname, atype) {_sign20 = 1;}

	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return type()==PLAbody?6:1; }
	virtual	void	position(const Point& r);
	virtual	int	nodes()		const	{ return type()==PLAbody?2:1; }
	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	rotate(const double, const Vector&);
	virtual	void	createQuads();
	virtual void	createMesh();

	// Dimension labels and show flags
	// @return NULL not to show anything otherwise label
	virtual	char*	showX()	const { return NULL; }
	virtual char*	showY()	const { return NULL; }
	virtual char*	showZ()	const { return NULL; }

	// bounding box
	virtual	BBox	_bbox() const;
	virtual	OBBox*	updateOBB(bool in = false) const;

	void		signMove20(const Point& r, const Vector& w);
protected:
	void	checkType();
	virtual	Location _locationWrt(const GBody *body) const;
}; // GPLABody

/* --- GSPHBody --- */
class GSPHBody : public GBody {
public:
	GSPHBody(const char *aname) : GBody(aname, SPHbody) {}

	double	R()	const		{ return xlen; }
	void	R(const double r)	{ xlen = ylen = zlen = r; }

	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return 4; }
	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	rotate(const double, const Vector&) {}
	virtual	void	createQuads();
	virtual void	createMesh() { createEllipsoidMesh(); }

	// Dimension labels and show flags
	// @return NULL not to show anything otherwise label
	virtual	char*	showX()	const { return "R"; }
	virtual char*	showY()	const { return NULL; }
	virtual char*	showZ()	const { return NULL; }

	// bounding box
	virtual	BBox	_bbox() const;
	virtual	OBBox*	updateOBB(bool in = false) const;
protected:
	virtual	Location _locationWrt(const GBody *body) const;
}; // GSPHBody

/* --- GELLBody --- */
class GELLBody : public GBody {
public:
	GELLBody(const char *aname) : GBody(aname, ELLbody) {}
	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return 7; }
	virtual	int	closest(const Point& r, const double d, const Vector& w) const;
	virtual	int	nodes()		const	{ return 3; }
	virtual Point	node(const int)	const;
	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	createQuads();
	virtual void	createMesh() { createEllipsoidMesh(); }

	// Dimension labels and show flags
	// @return NULL not to show anything otherwise label
	virtual	char*	showX()	const { return "Rm"; }
	virtual char*	showY()	const { return NULL; }
	virtual char*	showZ()	const { return "RM"; }

	// bounding box
	virtual	BBox	_bbox() const;
	virtual	OBBox*	updateOBB(bool in = false) const;
}; // GELLBody

/* --- GBOXBody --- */
class GBOXBody : public GBody {
public:
	GBOXBody(const char *aname, const BodyType atype)
		: GBody(aname, atype) {}
	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return type()==RPPbody?6:12; }
	virtual	void	position(const Point& r);
	virtual	int	nodes()		const	{ return type()==BOXbody?4:1; }
	virtual	void	rotate(const double angle, const Vector& axis);
	virtual	void	createQuads();
	virtual void	createMesh();
	virtual void	move(int item, const Point& r, const Vector& w);

	// bounding box
	virtual	OBBox*	updateOBB(bool in = false) const;
protected:
	virtual	Location _locationWrt(const GBody *body) const;
}; // GBOXBody

/* --- GWEDBody --- */
class GWEDBody : public GBody {
public:
	GWEDBody(const char *aname) : GBody(aname, WEDbody) {}

	Vector	N()	const;

	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return 12; }
	virtual	int	nodes()		const	{ return 4; }
	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	createQuads();
	virtual void	createMesh();

	// bounding box
	virtual	OBBox*	updateOBB(bool in = false) const;
}; // GWEDBody

/* --- GARBBody --- */
class GARBBody : public GBody {
private:
	int	faces[6];	// faces
public:
	GARBBody(const char *aname) : GBody(aname, ARBbody) {}
	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return 30; }
	//virtual	void	position(const Point& r);
	//virtual	void	save();
	//virtual	void	restore();
	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	rotate(const double angle, const Vector& axis);
	virtual	void	createQuads();
	virtual	void	createMesh();
private:
		bool	faceVertices(int f, int v[4], Point* V[4]);
}; // GARBBody

/* --- GQUABody --- */
class GQUABody : public GBody {
private:
	Quad	q;
	Quad	sq;
public:
	GQUABody(const char *aname) : GBody(aname, QUAbody) {}

	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return 10; }
	virtual	void	position(const Point& r);
	virtual	void	save();
	virtual	void	restore();
	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	rotate(const double angle, const Vector& axis);
	virtual	void	createQuads()	{ addQuad(q); }

	// Dimension labels and show flags
	// @return NULL not to show anything otherwise label
	virtual	char*	showX()	const { return NULL; }
	virtual char*	showY()	const { return NULL; }
	virtual char*	showZ()	const { return NULL; }

	// bounding box
	virtual	BBox	_bbox() const;
	virtual	OBBox*	updateOBB(bool in = false) const;
}; // GQUABody

/* --- GRCCBody --- */
class GRCCBody : public GBody {
public:
	GRCCBody(const char *aname) : GBody(aname, RCCbody) {}

	double	R()	const		{ return xlen; }
	void	R(const double r)	{ xlen = ylen = r; }

	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return 7; }
	virtual	int	nodes()		const	{ return 2; }
	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	createQuads();
	virtual void	createMesh()	{ createEllConeMesh(xlen,ylen,xlen,ylen); }

	// Dimension labels and show flags
	// @return NULL not to show anything otherwise label
	virtual	char*	showX()	const { return "R"; }
	virtual char*	showY()	const { return NULL; }
	virtual char*	showZ()	const { return "H"; }

	// bounding box
	virtual	BBox	_bbox() const;
	virtual	OBBox*	updateOBB(bool in = false) const;
}; // GRCCBody

/* --- GRECBody --- */
class GRECBody : public GBody {
public:
	GRECBody(const char *aname) : GBody(aname, RECbody) {}
	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return 12; }
	virtual	int	nodes()		const	{ return 4; }
	virtual	int	closest(const Point& r, const double d, const Vector& w) const;
	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	createQuads();
	virtual void	createMesh()	{ createEllConeMesh(xlen,ylen,xlen,ylen); }

	// Dimension labels and show flags
	// @return NULL not to show anything otherwise label
	virtual	char*	showX()	const { return "Rx"; }
	virtual char*	showY()	const { return "Ry"; }
	virtual char*	showZ()	const { return "H"; }

	// bounding box
	virtual	BBox	_bbox() const;
	virtual	OBBox*	updateOBB(bool in = false) const;
}; // GRECBody

/* --- GTRCBody --- */
class GTRCBody : public GBody {
private:
	double	h2;
public:
	GTRCBody(const char *aname) : GBody(aname, TRCbody) {}

	double	Rb()	const		{ return xlen; }
	double	Ra()	const		{ return ylen; }
	void	Rb(const double r)	{ xlen = r; }
	void	Ra(const double r)	{ ylen = r; }

	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return 8; }
	virtual	int	nodes()		const	{ return 2; }
	virtual	int	closest(const Point& r, const double d, const Vector& w) const;
	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	createQuads();
	virtual void	createMesh()	{ createEllConeMesh(Rb(),Rb(),Ra(),Ra()); }

	// Dimension labels and show flags
	// @return NULL not to show anything otherwise label
	virtual	char*	showX()	const { return "Rb"; }
	virtual char*	showY()	const { return "Ra"; }
	virtual char*	showZ()	const { return "H"; }

	// bounding box
	virtual	BBox	_bbox() const;
	virtual	OBBox*	updateOBB(bool in = false) const;
}; // GTRCBody

/* --- GInfEllCylBody --- */
class GInfEllCylBody : public GBody {
public:
	GInfEllCylBody(const char *aname, const BodyType atype)
		: GBody(aname, atype) {}

	Vector	direction() const;

	double	R()	const		{ return xlen; };
	double	Rx()	const		{ return xlen; };
	double	Ry()	const		{ return ylen; };
	void	R(const double r)	{ xlen = ylen = r; }
	void	Rx(const double r)	{ xlen = r; }
	void	Ry(const double r)	{ ylen = r; }

	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return type()==XCCbody || type()==YCCbody || type()==ZCCbody? 3 : 4; }
	virtual	int	closest(const Point& r, const double d, const Vector& w) const;
	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	rotate(const double angle, const Vector& axis);
	virtual	void	createQuads();
	virtual void	createMesh() { createEllConeMesh(xlen,ylen,xlen,ylen,true); }

	// Dimension labels and show flags
	// @return NULL not to show anything otherwise label
	virtual	char*	showX()	const {
				switch (type()) {
					case XCCbody:
					case YCCbody:
					case ZCCbody: return "R";
					case XECbody: return "Ry";
					case YECbody: return "Rz";
					case ZECbody: return "Rx";
					default: return NULL;
				}
			}
	virtual char*	showY()	const {
				switch (type()) {
					case XCCbody:
					case YCCbody:
					case ZCCbody: return NULL;
					case XECbody: return "Rz";
					case YECbody: return "Rx";
					case ZECbody: return "Ry";
					default: return NULL;
				}
			}
	virtual char*	showZ()	const { return NULL; }

	// bounding box
	virtual	BBox	_bbox() const;
	virtual	OBBox*	updateOBB(bool in = false) const;
protected:
	virtual	Location _locationWrt(const GBody *body) const;
}; // GInfEllCylBody

/* --- GTorusBody --- */
class GTorusBody : public GBody {
private:
	double	a, b, c;
public:
	GTorusBody(const char *aname, const BodyType atype) : GBody(aname, atype) { a = b = c = 0.0; }
	virtual	void	setWhat(double *what, char *err);
	virtual	int	getWhat(double *what) const;
	virtual int	nWhat()		const	{ return 6; }
	virtual	int	nodes()		const	{ return 3; }
//	virtual	int	closest(const Point& r, const double d, const Vector& w) const;
//	virtual void	move(int item, const Point& r, const Vector& w);
	virtual	void	createQuads();
	virtual void	createMesh();

	// Dimension labels and show flags
	// @return NULL not to show anything otherwise label
//	virtual	char*	showX()	const { return "Rx"; }
//	virtual char*	showY()	const { return "Ry"; }
//	virtual char*	showZ()	const { return "H"; }

	// bounding box
	virtual	BBox	_bbox() const;
//	virtual	OBBox*	updateOBB(bool in = false) const;
}; // GTorusBody

std::ostream& operator << (std::ostream&, const GBody&);
#endif
