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
 * Date:	05-Apr-2011
 *
 * TODO
 * GObject has to be split in two GObject and VObject
 * otherwise we will have problems when parallel drawings are done
 */

#ifndef __GOBJECT_H
#define __GOBJECT_H

#include <Python.h>

#include <X11/X.h>
#include <X11/Xutil.h>

#include "os.h"
#include "geo.h"
#include "bbox.h"
#include "viewer.h"
#include "geometry.h"

#define MAX_NODES	20	// maximum spline nodes

struct ViewerObject;

enum ClassType {
	GObjectClass  = 0,   // option offset also
	GPointClass   = 10,
	GArrowClass   = 20,
	GRulerClass   = 30,
	GRotdefiClass = 40,
	GLightClass   = 50,
	GBeamClass    = 60,
	GSplineClass  = 100, // allow for up to 10 points
	GMeshClass    = 200
};

// Text location
enum AnchorType {
	Anchor_none,	// Do not display text
	Anchor_C,	// Center
	Anchor_N,	// North
	Anchor_NE,	// North-East
	Anchor_E,	// East
	Anchor_SE,	// South-East
	Anchor_S,	// South
	Anchor_SW,	// South-West
	Anchor_W,	// West
	Anchor_NW	// North-West
};

// Object type
enum ObjectType {
	Object_Default =   0,

	Point_Cross    =  GPointClass,
	Point_Dot      =  GPointClass+1,
	Point_Square   =  GPointClass+2,
	Point_X        =  GPointClass+3,
	Point_Diamond  =  GPointClass+4,
	Point_Circle   =  GPointClass+5,

	Arrow_Head     =  GArrowClass,
	Arrow_Tail     =  GArrowClass+1,
	Arrow_HeadTail =  GArrowClass+2,
	Arrow_Line     =  GArrowClass+3,

	Ruler_Simple   =  GRulerClass,
	Ruler_Angle    =  GRulerClass+1,

	Rotdefi_XYZ    =  GRotdefiClass,

	Light_Sun      =  GLightClass,
	Light_Omni     =  GLightClass+1,
	Light_Spot     =  GLightClass+2,

	Beam_Arrow     =  GBeamClass,
        Beam_Divergence=  GBeamClass+1,

	Mesh_Default   =  GMeshClass,
};

/** GObject graphics object class */
class GObject {
public:
	ObjectType	option;		// object type to draw
	AnchorType	anchor;		// text location
	dword		color;		// drawing color
	dword		select;		// if !=0 then is selected and contains the selection color
	int		selnode;	// selected node
					// -1 = nothing, 0=object is selected, 1=first node, 2=second...
	bool		show;		// show object
	double		drawDw;		// visibility limit
	int		size;		// object dimension
	int		lineWidth;	// line width to use

protected:
	char32		_name;		// name of object
	int		_id;		// id of object
	bool		_visible;	// if object is visible in viewport
	Point		P;		// absolute position

	Point		V, Vc;		// viewer position, and clipped
	int		x1, y1;		// object position

	Point		SP;		// save position

	XGCValues	gcValues;	// context values
//	Font		font;
static	dword		nodeColor;	// selection color for nodes

public:
	GObject(const char *aname=NULL, ObjectType o=Object_Default);
	virtual ~GObject() {}

	// set/get
	virtual ClassType isA() const = 0;
	virtual const char *type() const = 0;

const	char*	name()			const	{ return _name; }
	void	name(const char *aname);
	void	id(const int i)			{ _id = i; }
	int	id()			const	{ return _id; }

	void	position(const double x, const double y, const double z) { P.x=x; P.y=y; P.z=z; }
	void	position(const Point& r)	{ P = r; }
	Point	position()		const	{ return P; }
	Point	savedPosition()		const	{ return SP; }
	Point	viewerPosition()	const	{ return V; }

	int	config(PyObject*);
	virtual PyObject* config(const char *key, PyObject *value=NULL);

	bool	visible() const { return _visible && Abs(V.z)<=drawDw; }

	// save&recall
		void	clearSave()		{ selnode = -1; }
	virtual	void	save()			{ SP = P; }
	virtual void	restore()		{ P = SP; }
	virtual Vector	savedNode(int)	const	{ return SP; }

	// transform
	virtual int	nodes() const = 0;
	virtual Vector	node(int)	const	{ return P; }
	virtual	void	node(int item, const Vector& r);
	virtual	void	rotate(const double, const Vector&) {}
	virtual void	transform(ViewerObject *self);

	// bounding box
	virtual BBox	bbox() const;
	virtual BBox	bboxView(ViewerObject *self);

	// drawing routines
	virtual void	draw(ViewerObject *self, Drawable drawable);
	virtual void	drawText(ViewerObject *self, Drawable drawable);

	// closer routines
	virtual int	closest(ViewerObject *self, const int i, const int j, const int d);
	virtual bool	enclosed(ViewerObject *self, const int l, const int t, const int r, const int b);

protected:
	void setForeground(ViewerObject* self);
	void drawSelectedPoint(ViewerObject* self, Drawable drawable, int i, int j);
}; // GObject

/** GPoint */
class GPoint : public GObject {
public:
	GPoint(const char *aname=NULL, ObjectType o=Point_Cross) : GObject(aname, o) {};
	virtual ~GPoint() {}

	virtual ClassType isA()		const	{ return GPointClass; }
	virtual const char *type()	const	{ return "point"; }

	// nodes
	virtual int nodes() const { return 1; }

	virtual void draw(ViewerObject *self, Drawable drawable);
}; // GPoint

/** GArrow */
class GArrow : public GObject {
protected:
	Vector	D;			// vector length
	Vector	SD;			// save variables

	Point	Ve, Vec;		// Viewer end position (clipped)
	int	x2, y2;

	bool	drawHead;		// draw a head arrow
	bool	drawTail;		// draw a tail arrow

public:
	GArrow(const char *aname=NULL, ObjectType o=Arrow_Head)
		: GObject(aname, o), D(10.,10.,10.), x2(0), y2(0), drawHead(true), drawTail(false) {}
	virtual ~GArrow() {}

	virtual ClassType isA()		const	{ return GArrowClass; }
	virtual const char *type()	const	{ return "arrow"; }
	virtual	PyObject* config(const char *key, PyObject *value=NULL);

		void	direction(const double x, const double y, const double z) { D.x=x; D.y=y; D.z=z; }
		void	direction(const Vector& d)	{ D = d; }
		Vector	direction()	const	{ return D; }
		Point	viewerEnd()	const	{ return Ve; }

	// save & recall
	virtual	void	save();
	virtual void	restore();
	virtual Vector	savedNode(int)	const;

	// transform
	virtual	int	nodes()		const	{ return 2; }
	virtual Vector	node(int)	const;
	virtual	void	node(int item, const Vector& r);
	virtual	void	rotate(const double angle, const Vector& axis);
	virtual void	transform(ViewerObject *self);

	// bounding box
	virtual BBox	bbox() const;
	virtual BBox	bboxView(ViewerObject *self);

	// drawing routines
	virtual void	draw(ViewerObject *self, Drawable drawable);
	virtual void	drawText(ViewerObject *self, Drawable drawable);

	// closer routines
	virtual int	closest(ViewerObject *self, const int i, const int j, const int d);
	virtual bool	enclosed(ViewerObject *self, const int left, const int top,
				const int right, const int bottom);
}; // GArrow

/** GSpline */
class GSpline : public GObject {
protected:
	int	n;			// number of nodes
	Point	point[MAX_NODES];	// up to MAX_NODES nodes!!!! node[0] == P
	Point	Vpoint[MAX_NODES];

public:
	GSpline(const char *aname=NULL, ObjectType o=Point_Cross) : GObject(aname, o), n(1) { point[0] = Vector::O; };
	virtual ~GSpline() {}

	virtual ClassType isA()		const	{ return GSplineClass; }
	virtual const char *type()	const	{ return "spline"; }
	virtual	PyObject* config(const char *key, PyObject *value=NULL);

	// nodes
	virtual int nodes() const { return n; }
	virtual void	transform(ViewerObject *self);

	virtual void draw(ViewerObject *self, Drawable drawable);
}; // GSpline

/** GRuler */
class GRuler : public GArrow {
protected:
	Vector	Da;		// angle point
	Vector	SDa;		// save

	Point	Va, Vac;	// viewer angle position + clipped
	Point	Vc2;		// viewer angle-starting position clipped

	int	x3, y3;

public:
	GRuler(const char *aname=NULL, ObjectType o=Ruler_Simple)
		: GArrow(aname,o), Da(10.,0.,0.) {}
	virtual ~GRuler() {}

	virtual ClassType isA()		const	{ return GRulerClass; }
	virtual const char *type()	const	{ return "ruler"; }

	virtual PyObject* config(const char *key, PyObject *value=NULL);

		Vector	viewerAngle()	const	{ return Va; }

	// save & recall
	virtual	void	save();
	virtual void	restore();
	virtual Vector	savedNode(int)	const;

	// transform
	virtual int	nodes()		const	{ return (option==Ruler_Angle?3:2); }
	virtual Vector	node(int)	const;
	virtual	void	node(int item, const Vector& r);
	virtual	void	rotate(const double angle, const Vector& axis);
	virtual void	transform(ViewerObject *self);

	// bounding box
	virtual BBox	bbox() const;
	virtual BBox	bboxView(ViewerObject *self);

	// drawing routines
	virtual void	draw(ViewerObject *self, Drawable drawable);
	virtual void	drawText(ViewerObject *self, Drawable drawable);

	// closer routines
	virtual int	closest(ViewerObject *self, const int i, const int j, const int d);
	virtual bool	enclosed(ViewerObject *self, const int left, const int top,
				const int right, const int bottom);
}; // GRuler

/** GRotdefi */
class GRotdefi : public GArrow {
public:
	Matrix4	orient;			// Initial orientation/position of base
	Matrix4	matrix;			// Rotdefi matrix for projection
	int	axisSize;		// axis size
	int	axisWidth;		// axis dimension

public:
	GRotdefi(const char *aname=NULL) : GArrow(aname,Arrow_Head), orient(1), matrix(1), axisSize(20)
			{ axisWidth = 3; }
	virtual ~GRotdefi() {}

	virtual ClassType isA()		const	{ return GRotdefiClass; }
	virtual const char *type()	const	{ return "axis"; }

	virtual PyObject* config(const char *key, PyObject *value=NULL);

	// nodes
//	virtual int nodes() const { return 4; }

	// bounding box
//	virtual BBox	bbox() const;
//	virtual BBox	bboxView(ViewerObject *self);

	// drawing routines
	virtual void	draw(ViewerObject *self, Drawable drawable);
//	virtual void	drawText(ViewerObject *self, Drawable drawable);

	// closer routines
//	virtual int	closest(ViewerObject *self, const int i, const int j, const int d);
	virtual bool	enclosed(ViewerObject *self, const int left, const int top,
				const int right, const int bottom);
}; // GRotdefi

/** GLight */
class GLight : public GArrow {
public:
	double	power, spec;		// light power and specular
	bool	relative;		// relative to origin or absolute
	bool	shadow;			// shadow casting
	int	falloff;		// falloff type 0,1,2

public:
	GLight(const char *aname=NULL, ObjectType o=Light_Sun);
	virtual ~GLight() {}

	virtual ClassType isA()		const	{ return GLightClass; }
	virtual const char *type()	const	{ return "light"; }

	virtual PyObject* config(const char *key, PyObject *value=NULL);

	// transform
	virtual int	nodes()		const	{ return 2; }
	virtual	void	node(const int item, const Vector& r);
	//virtual void	transform(ViewerObject *self);

	virtual void	draw(ViewerObject *self, Drawable drawable);
	virtual void	drawText(ViewerObject *self, Drawable drawable);

	// closer routines
	virtual int	closest(ViewerObject *self, const int i, const int j, const int d);

	// convert to light structure
	void toLight(Light *L) const;
}; // GLight

/** GBeam */
class GBeam : public GArrow {
public:
	double	energy;         // beam energy
	double	scale;		// scaling of arrow
	double	divergence;     // beam divergence
	double	Rin, Rout;	// ring radius

public:
	GBeam(const char *aname=NULL, ObjectType o=Beam_Arrow);
	virtual ~GBeam() {}

	virtual ClassType isA()		const	{ return GBeamClass; }
	virtual const char *type()	const	{ return "beam"; }

	virtual PyObject* config(const char *key, PyObject *value=NULL);

	// transform
//	virtual int	nodes()		const	{ return 2; }
//	virtual	void	node(const int item, const Vector& r);
//	//virtual void	transform(ViewerObject *self);

	virtual void	draw(ViewerObject *self, Drawable drawable);
//	virtual void	drawText(ViewerObject *self, Drawable drawable);

	// closer routines
//	virtual int	closest(ViewerObject *self, const int i, const int j, const int d);
}; // GBeam

#endif
