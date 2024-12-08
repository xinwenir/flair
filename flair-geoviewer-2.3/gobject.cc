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

#include <Python.h>

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "xdraw.h"
#include "pyutils.h"
#include "gobject.h"
#include "viewerobject.h"

using namespace std;

static const char *anchorStr[] = {
			"none", "C", "N",  "NE", "E",
			"SE",   "S", "SW", "W",  "NW" };

//static const char *objectTypeStr[]  = {
//		"default",	// Object_Default
//		"square",	// Point_Square
//		"cross",	// Point_Cross
//		"x",		// Point_X
//		"diamond",	// Point_Diamond
//		"circle",	// Point_Circle
//		"line",		// Arrow_Line
//		"head",		// Arrow_Head
//		"tail",		// Arrow_Tail
//		"headtail",	// Arrow_HeadTail
//		"simple",	// Ruler_Simple
//		"angle",	// Ruler_Angle
//		"xyz",		// Rotdefi_XYZ
//		"uvw",		// Rotdefi_UVW
//		"sun",		// Light_Sun
//		"omni",		// Light_Omni
//		"spot",		// Light_Spot
//		};

static const char dottedPattern[2] = {1,1};
static const char dashedPattern[2] = {3,3};

// Selection color of items
dword GObject::nodeColor = 0xFFB0B0;

/* ================================ GObject ================================= */
/** Basic class for all graphics objects */
GObject::GObject(const char *aname, ObjectType o)
{
	name(aname);
	option    = o;
	anchor    = Anchor_none;
	color     = 0x000000;    // Color
	drawDw    = INFINITE;    // limit
	size      = 3;           // pixels
	lineWidth = 1;           // default line width
	show      = true;
	select    = 0x000000;    // Selection Color
	selnode   = -1;          // None node is selected

	x1 = y1   = 0;

	gcValues.function   = GXcopy;
	gcValues.line_width = lineWidth;
	gcValues.line_style = LineSolid;
} // GObject

/** name - set the name of the object */
void GObject::name(const char *aname)
{
	if (aname) {
		strncpy(_name, aname, sizeof(_name)-1);
		_name[sizeof(_name)-1] = 0;
	} else
		_name[0] = 0;
} // name

/** setForeground */
inline void GObject::setForeground(ViewerObject *self)
{
	XSetForeground(self->display, self->gc, select? select:color);
} // setForeground

/** drawSelectedPoint */
inline void GObject::drawSelectedPoint(ViewerObject *self, Drawable drawable, int i, int j)
{
	XSetForeground(self->display, self->gc, nodeColor);
	XDrawRectangle(self->display, drawable, self->gc, i-2, j-2, 4, 4);
} // drawSelectedPoint

/** set parameters from dictionary
 * @param dict	dictionary with parameters to set
 * @param return true on error and set the python error
 * @return number of error parameters but not the names neither sets an exception
 */
int GObject::config(PyObject *dict)
{
	if (!PyDict_Check(dict)) {
		PyErr_SetString(PyExc_TypeError, "Invalid type, dictionary expected");
		return 1;
	}

	PyObject *key, *value;
#if PY_VERSION_HEX >= 0x02050300
	Py_ssize_t pos = 0;
#else
	int pos = 0;
#endif

	// Firstly set position if any
	value = PyDict_GetItemString(dict, "x");
	if (value) config("x", value);
	value = PyDict_GetItemString(dict, "y");
	if (value) config("y", value);
	value = PyDict_GetItemString(dict, "z");
	if (value) config("z", value);

	// loop over parameters
	int errcount = 0;
	while (PyDict_Next(dict, &pos, &key, &value))
		if (config(PyUnicode_AsUTF8(key), value) == NULL)
			errcount++;

	if (errcount) PyErr_Clear();

	return errcount;
} // config

/** set parameter
 * @param key	value to set
 * @param value	python object to set value from
 * @param return true on error and set the python error
 */
PyObject* GObject::config(const char *key, PyObject *value)
{
	if (!strcmp(key, "pos")) {
		if (value==NULL)
			return Py_Vector(P);
		else
			P = Py_GetVector(value);
	} else
	if (!strcmp(key, "move") && selnode<=1) {
		if (value==NULL) {
			PyErr_SetString(PyExc_SyntaxError, "'move' do not return any value");
			return NULL;
		} else {
			Vector M = Py_GetVector(value);
			P += M;
		}
	} else
	if (!strcmp(key, "x")) {
		if (value==NULL)
			return PyFloat_FromDouble(P.x);
		else
			P.x = Py_GetFloat(value);
	} else
	if (!strcmp(key, "y")) {
		if (value==NULL)
			return PyFloat_FromDouble(P.y);
		else
			P.y = Py_GetFloat(value);
	} else
	if (!strcmp(key, "z")) {
		if (value==NULL)
			return PyFloat_FromDouble(P.z);
		else
			P.z = Py_GetFloat(value);
	} else
	if (!strcmp(key, "save"))
		save();
	else
	if (!strcmp(key, "restore"))
		restore();
	else
	if (!strcmp(key, "clearsave"))
		selnode = -1;
	else
	if (!strcmp(key, "name")) {
		if (value==NULL)
			return PyString_FromString(name());
		else
			name(PyUnicode_AsUTF8(value));
	} else
	if (!strcmp(key, "move")) {
		Vector M = Py_GetVector(value);
		P += M;
	} else
	if (!strcmp(key, "type")) {
		if (value==NULL)
			return PyString_FromString(type());
		// No error is generated from this, since it is
		// set automatically during the initial config
	} else
	if (!strcmp(key, "id")) {
		if (value==NULL)
			return PyInt_FromLong(id());
		else {
			PyErr_SetString(PyExc_SyntaxError, "Cannot set the 'id'");
			return NULL;
		}
	} else
	if (!strcmp(key, "anchor")) {
		if (value==NULL)
			return PyString_FromString(anchorStr[anchor]);
		else
		if (PyInt_Check(value))
			anchor = (AnchorType)PyInt_AsLong(value);
		else {
			const char *s = PyUnicode_AsUTF8(value);
			for (unsigned i=0; i<sizeof(anchorStr)/sizeof(char*); i++)
				if (!strcmp(s, anchorStr[i])) {
					anchor = (AnchorType)i;
					Py_RETURN_NONE;
				}
			PyErr_Format(PyExc_TypeError, "Invalid anchor \"%s\"",s);
			return NULL;
		}
	} else
	if (!strcmp(key, "color")) {
		if (value==NULL)
			return PyInt_FromLong((dword)color);
		else {
			if (PyInt_Check(value))
				color = (dword)PyInt_AsLong(value);
			else
				color = (dword)strtol(PyUnicode_AsUTF8(value)+1, NULL, 16);
		}
	} else
	if (!strcmp(key, "select")) {
		if (value==NULL)
			return PyInt_FromLong((dword)select);
		else {
			if (PyInt_Check(value))
				select = (dword)PyInt_AsLong(value);
			else
				select = (dword)strtol(PyUnicode_AsUTF8(value)+1, NULL, 16);
			selnode = (select==0?-1:0); // clear or set node
		}
	} else
//	if (!strcmp(key, "selnode")) {
//		if (value==NULL)
//			return PyInt_FromLong(selnode);
//		else
//			selnode = PyInt_AsLong(value);
//	} else
	if (!strcmp(key, "nodecolor")) {
		if (value==NULL)
			return PyInt_FromLong((dword)nodeColor);
		else {
			if (PyInt_Check(value))
				nodeColor = (dword)PyInt_AsLong(value);
			else
				nodeColor = (dword)strtol(PyUnicode_AsUTF8(value)+1, NULL, 16);
		}
	} else
	if (!strcmp(key, "nodes")) {
		if (value==NULL)
			return PyInt_FromLong(nodes());
		else {
			PyErr_SetString(PyExc_SyntaxError, "Cannot set the number of nodes");
			return NULL;
		}
	} else
	if (!strcmp(key, "show"))
		if (value==NULL)
			return PyInt_FromLong((dword)show);
		else
			show = (bool)PyInt_AsLong(value);
	else
	if (!strcmp(key, "size")) {
		if (value==NULL)
			return PyInt_FromLong((long)size);
		else {
			size = Py_GetInt(value);
			if (size<=0) size = 5;
		}
	} else
	if (!strcmp(key, "linewidth")) {
		if (value==NULL)
			return PyInt_FromLong((long)lineWidth);
		else {
			lineWidth = Py_GetInt(value);
			if (lineWidth<=0) lineWidth = 0;
		}
	} else
	if (!strcmp(key, "drawDw")) {
		if (value==NULL)
			return PyFloat_FromDouble(drawDw);
		else
			drawDw = Py_GetFloat(value);
	} else
	if (!strcmp(key, "option")) {
		if (value==NULL)
			return PyInt_FromLong(option - (int)isA());
//			return PyString_FromString(objectTypeStr[option]);
		else
		if (PyInt_Check(value))
			option = (ObjectType)(PyInt_AsLong(value) + (int)isA());
//		else {
//			const char *s = PyString_AsString(value);
//			for (unsigned i=0; i<sizeof(objectTypeStr)/sizeof(char*); i++)
//				if (!strcmp(s, objectTypeStr[i])) {
//					option = (ObjectType)i;
//					Py_RETURN_NONE;
//				}
//			PyErr_Format(PyExc_TypeError, "Invalid option \"%s\"",s);
//			return NULL;
//		}
	} else
	if (!strcmp(key, "dict")) {
		if (value==NULL) {
			PyErr_SetString(PyExc_SyntaxError, "'dict' do not return any value");
			return NULL;
		} else
			config(value);
	} else {
		PyErr_Format(PyExc_KeyError,"Object: Invalid option \"%s\"", key);
		return NULL;
	}

	if (PyErr_Occurred()) return NULL;
	Py_RETURN_NONE;
} // config

/** move object item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 */
void GObject::node(int, const Vector& r)
{
	P = r;
} // node

/** transform coordinates to screen */
void GObject::transform(ViewerObject *self)
{
	self->kernel->view.xyz2uvw3D(P, &V);
	Vc = V;
	_visible = self->kernel->view.inside(V.x,V.y);
} // transform

/** @return bounding box of object */
BBox GObject::bbox() const
{
	BBox bb;
	bb.add(P);
	return bb;
} // bbox

/** @return bounding box of object */
BBox GObject::bboxView(ViewerObject *self)
{
	BBox bb;
	transform(self);
	bb.add(V);
	return bb;
} // bboxView

/** draw */
void GObject::draw(ViewerObject *self, Drawable drawable)
{
	gcValues.line_width = lineWidth;
	if (Eq0(V.z,SMALL)) {
		gcValues.line_style = LineSolid;
	} else
	if (V.z<0.0) {
		gcValues.line_style = LineOnOffDash;
		XSetDashes(self->display, self->gc, 0,
			dottedPattern, SIZE(dottedPattern));
	} else
		gcValues.line_style = LineSolid;

	XChangeGC(self->display, self->gc,
		GCFunction|GCLineWidth|GCLineStyle, &gcValues);

	if (select && selnode==1 && _visible && self->kernel->view.inside(V.x,V.y)) {
		x1 = self->kernel->view.u2i(V.x);
		y1 = self->kernel->view.v2j(V.y);
		drawSelectedPoint(self, drawable,
			self->kernel->view.u2i(V.x),
			self->kernel->view.v2j(V.y));
	}
	x1 = self->kernel->view.u2i(Vc.x);
	y1 = self->kernel->view.v2j(Vc.y);
	setForeground(self);
	if (_name[0] && anchor!=Anchor_none) drawText(self, drawable);
} // draw

/** drawText */
void GObject::drawText(ViewerObject *self, Drawable drawable)
{
	int i = x1;
	int j = y1;
	int l = (int)strlen(_name);
//	int V.z = XTextWidth(font, name, l);

	switch (anchor) {
		case Anchor_N:
			j -= size+1;
			break;
		case Anchor_NE:
			j -= size+1;
			i += size+1;
			break;
		case Anchor_E:
			i += size+1;
			break;
		case Anchor_SE:
			j += size+1;
			i += size+1;
			break;
		case Anchor_S:
			j += size+1;
			break;
		case Anchor_SW:
			j += size+1;
			i -= size+1;
			break;
		case Anchor_W:
			i -= size+1;
			break;
		case Anchor_NW:
			j -= size+1;
			i -= size+1;
			break;
		default: ;	// just to keep compiler happy
	}
	XDrawString(self->display, drawable, self->gc, i, j, _name, l);
} // drawText

/** closest
 * Check if a point lies close to the object assuming also the size of the object
 * @param i,j	check point
 * @param d	check distance
 * @return -1 on failure,
 *          0 on object
 *          # index of item which is close to hit point i,j
 */
int GObject::closest(ViewerObject *self, const int i, const int j, const int d)
{
	x1 = self->kernel->view.u2i(V.x);
	y1 = self->kernel->view.v2j(V.y);

	if (Sqr((double)(i-x1)) + Sqr((double)(j-y1)) <= (double)Sqr(size+d))
		return  1;	// Main point
	else
		return -1;	// Nowhere
} // closest

/** enclosed */
bool GObject::enclosed(ViewerObject *self, const int left, const int top, const int right, const int bottom)
{
	x1 = self->kernel->view.u2i(Vc.x);
	y1 = self->kernel->view.v2j(Vc.y);

	int xl = Max(x1-size, left);
	int yl = Max(y1-size, top);
	int xh = Min(x1+size, right);
	int yh = Min(y1+size, bottom);

	return (xl<=xh && yl<=yh);
} // enclosed

/* ================================= GPoint ================================= */
/** draw */
void GPoint::draw(ViewerObject *self, Drawable drawable)
{
	int size2;
	XPoint	pts[5];

	GObject::draw(self, drawable);

	switch (option) {
		case Point_Dot:	// Draw a dot
			XDrawPoint(self->display, drawable, self->gc, x1, y1);
			break;

		case Point_Square:
			size2 = size*2 + 1;
			XDrawRectangle(self->display, drawable, self->gc, x1-size, y1-size, size2, size2);
			XDrawPoint(self->display, drawable, self->gc, x1, y1);
			break;

		case Point_X:
			XDrawLine(self->display, drawable, self->gc, x1-size, y1-size, x1+size+1, y1+size+1);
			XDrawLine(self->display, drawable, self->gc, x1-size, y1+size, x1+size+1, y1-size-1);
			break;

		case Point_Diamond:
			pts[0].x = pts[4].x = (short)(x1);
			pts[0].y = pts[4].y = (short)(y1 - size);
			pts[1].x = (short)(x1 + size);
			pts[1].y = (short)(y1);
			pts[2].x = (short)(x1);
			pts[2].y = (short)(y1 + size);
			pts[3].x = (short)(x1 - size);
			pts[3].y = (short)(y1);
			XDrawLines(self->display, drawable, self->gc, pts, 5, CoordModeOrigin);
			XDrawPoint(self->display, drawable, self->gc, x1, y1);
			break;

		case Point_Circle:
			size2 = size*2 + 1;
			XDrawPoint(self->display, drawable, self->gc, x1, y1);
			XDrawArc(self->display, drawable, self->gc, x1-size, y1-size, size2, size2, 0, 360*64);
			break;

		default:
			XDrawLine(self->display, drawable, self->gc, x1-size, y1, x1+size+1, y1);
			XDrawLine(self->display, drawable, self->gc, x1, y1-size, x1, y1+size+1);
	}
} // draw

/* ================================= GArrow ================================== */
/** set parameter
 * @param key	value to set
 * @param value	python object to set value from
 * @param return true on error and set the python error
 */
PyObject* GArrow::config(const char *key, PyObject *value)
{
	if (!strcmp(key, "pos")) {	// Override default behaviour of GObject
		switch (selnode) {
			case 0:	// Whole arrow
				if (value==NULL)
					return Py_Vector(P);
				else
					P = Py_GetPoint(value);
				break;
			case 1: // Start point leaving end point untouched
				if (value==NULL)
					return Py_Vector(P);
				else {
					Point O = P;
					P = Py_GetPoint(value);
					D += O-P;
				}
				break;
			case 2:	// End point
				if (value==NULL)
					return Py_Vector(P+D);
				else {
					D = Py_GetVector(value);
					D -= P;
				}
				break;
		}
	} else
	if (!strcmp(key, "dx")) {
		if (value==NULL)
			return PyFloat_FromDouble(D.x);
		else
			D.x = Py_GetFloat(value);
	} else
	if (!strcmp(key, "dy")) {
		if (value==NULL)
			return PyFloat_FromDouble(D.y);
		else
			D.y = Py_GetFloat(value);
	} else
	if (!strcmp(key, "dz")) {
		if (value==NULL)
			return PyFloat_FromDouble(D.z);
		else
			D.z = Py_GetFloat(value);
	} else
	if (!strcmp(key, "dir")) {
		if (value==NULL)
			return Py_Vector(D);
		else
			D = Py_GetVector(value);
	} else
	if (!strcmp(key, "xe")) {
		if (value==NULL)
			return PyFloat_FromDouble(P.x+D.x);
		else
			D.x = Py_GetFloat(value) - P.x;
	} else
	if (!strcmp(key, "ye")) {
		if (value==NULL)
			return PyFloat_FromDouble(P.y+D.y);
		else
			D.y = Py_GetFloat(value) - P.y;
	} else
	if (!strcmp(key, "ze")) {
		if (value==NULL)
			return PyFloat_FromDouble(P.z+D.z);
		else
			D.z = Py_GetFloat(value) - P.z;
	} else
	if (!strcmp(key, "end")) {
		if (value==NULL)
			return Py_Vector(P+D);
		else {
			D = Py_GetVector(value);
			D -= P;
		}
	} else
	if (!strcmp(key, "option")) {
		PyObject* obj = GObject::config(key, value);
		switch (option) {
			case Arrow_HeadTail:
				drawHead = true;
				drawTail = true;
				break;
			case Arrow_Head:
				drawHead = true;
				drawTail = false;
				break;
			case Arrow_Tail:
				drawHead = false;
				drawTail = true;
				break;
			default:
				drawHead = false;
				drawTail = false;
				break;
		}
		return obj;
	} else
		return GObject::config(key, value);

	if (PyErr_Occurred()) return NULL;
	Py_RETURN_NONE;
} // config

/** move object item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 */
void GArrow::node(int item, const Vector& r)
{
	if (item==1) {
		P = r;
		D = SP + SD - r;
	} else
	if (item==2) {
		D  = r;
		D -= P;
	} else
		GObject::node(item,r);
} // node

/** @return node position */
Vector GArrow::node(int item) const
{
	if (item==2)
		return P + D;
	else
		return GObject::node(item);
} // node

/** @return savedNode position */
Vector GArrow::savedNode(int item) const
{
	if (item==2)
		return SP + SD;
	else
		return GObject::savedNode(item);
} // savedNode

/** transform */
void GArrow::transform(ViewerObject *self)
{
//	GObject::transform(self);
//	self->kernel->view.xyz2uvw3D(P+D, &Ve);
//	Vec = Ve;
//	_visible = self->kernel->view.clipLine(&Vc.x,&Vc.y, &Vec.x,&Vec.y);
	// FIXME WRONG Vc = V
	_visible = self->kernel->view.clipLine3D(P, P+D, &V, &Ve);
	Vc  = V;
	Vec = Ve;
} // transform

/** rotate object around an axis
 * @param angle	to rotate
 * @param axis	to rotate
 */
void GArrow::rotate(const double angle, const Vector& axis)
{
	Matrix4	rot;
	rot.rotate(angle, axis);
	D = rot * SD;
} // rotate

/** save */
void GArrow::save()
{
	GObject::save();
	SD = D;
} // save

/** restore */
void GArrow::restore()
{
	GObject::restore();
	D = SD;
} // restore

/** @return bounding box of object */
BBox GArrow::bbox() const
{
	BBox bb;
	bb.add(P);
	bb.add(P+D);
	return bb;
} // bbox

/** @return bounding box of object */
BBox GArrow::bboxView(ViewerObject *self)
{
	BBox bb;
	transform(self);
	bb.add(V);
	bb.add(Ve);
	return bb;
} // bboxView

/** draw */
void GArrow::draw(ViewerObject *self, Drawable drawable)
{
	GObject::draw(self, drawable);

	// Using the clipped relative-coordinates
	double x1d = self->kernel->view.u2id(Vc.x);
	double y1d = self->kernel->view.v2jd(Vc.y);
	double x2d = self->kernel->view.u2id(Vec.x);
	double y2d = self->kernel->view.v2jd(Vec.y);

	x1 = Round(x1d);
	y1 = Round(y1d);
	x2 = Round(x2d);
	y2 = Round(y2d);

	XDrawLine3D(self->display, drawable, self->gc, self->kernel->view, Vc, Vec);

	// Draw mid-point
	XSetForeground(self->display, self->gc, self->geometry->geometry->vertexColor);
	XDrawPoint(self->display, drawable, self->gc, (x1+x2)/2, (y1+y2)/2);
	setForeground(self);

	//XDrawLine(self->display, drawable, self->gc, x1, y1, x2, y2);
	if (option==Arrow_Line) return;

	// Switch back to solid
	gcValues.line_style = LineSolid;
	XChangeGC(self->display, self->gc, GCLineStyle, &gcValues);

	// Draw heads
	XPoint	pts[3];
	double ddx = x2d - x1d;
	double ddy = y2d - y1d;
	double length = hypot(ddx,ddy);
	if (length<0.00001) return;
	ddx /= length;
	ddy /= length;

	// Draw head  ----->
	if (self->kernel->view.inside(Ve.x,Ve.y)) {
		if (drawHead) {
			pts[0].x = (short)Round(x1d + ddx*(length-size) - ddy*size);
			pts[0].y = (short)Round(y1d + ddy*(length-size) + ddx*size);
			pts[1].x = (short)x2;
			pts[1].y = (short)y2;
			pts[2].x = (short)Round(x1d + ddx*(length-size) + ddy*size);
			pts[2].y = (short)Round(y1d + ddy*(length-size) - ddx*size);
			XDrawLines(self->display, drawable, self->gc, pts, 3, CoordModeOrigin);
		}
		if (select && selnode==2 && self->kernel->view.inside(Ve.x,Ve.y)) {
			drawSelectedPoint(self, drawable,
				self->kernel->view.u2i(Ve.x),
				self->kernel->view.v2j(Ve.y));
			setForeground(self);
		}
	}

	// Draw tail <-----
	if (drawTail && self->kernel->view.inside(V.x,V.y)) {
		pts[0].x = (short)Round(x1d + ddx*size - ddy*size);
		pts[0].y = (short)Round(y1d + ddy*size + ddx*size);
		pts[1].x = (short)x1;
		pts[1].y = (short)y1;
		pts[2].x = (short)Round(x1d + ddx*size + ddy*size);
		pts[2].y = (short)Round(y1d + ddy*size - ddx*size);
		XDrawLines(self->display, drawable, self->gc, pts, 3, CoordModeOrigin);
	}
} // draw

/** drawText */
void GArrow::drawText(ViewerObject *self, Drawable drawable)
{
	x1 = self->kernel->view.u2i(Vc.x);
	y1 = self->kernel->view.v2j(Vc.y);
	x2 = self->kernel->view.u2i(Vec.x);
	y2 = self->kernel->view.v2j(Vec.y);

	int i = (x1+x2)/2;
	int j = (y1+y2)/2;

	int l = (int)strlen(_name);
//	int V.z = XTextWidth(font, name, l);

	switch (anchor) {
		case Anchor_N:
			j -= size+1;
			break;
		case Anchor_NE:
			j -= size+1;
			i += size+1;
			break;
		case Anchor_E:
			if (x1>x2) {
				i = x1 + size+1;
				j = y1;
			} else {
				i = x2 + size+1;
				j = y2;
			}
			break;
		case Anchor_SE:
			j += size+1;
			i += size+1;
			break;
		case Anchor_S:
			j += size+1;
			break;
		case Anchor_SW:
			j += size+1;
			i -= size+1;
			break;
		case Anchor_W:
			if (x1<x2) {
				i = x1 - size+1;
				j = y1;
			} else {
				i = x2 - size+1;
				j = y2;
			}
			break;
		case Anchor_NW:
			j -= size+1;
			i -= size+1;
			break;
		default: ;	// just to keep compiler happy
	}
	XDrawString(self->display, drawable, self->gc, i, j, _name, l);
} // drawText

/** closest */
int GArrow::closest(ViewerObject *self, const int i, const int j, const int d)
{
	// Check main point
	int close = GObject::closest(self, i, j, d);
	if (close>=0) return close;

	// Check other end
	x2 = self->kernel->view.u2i(Ve.x);
	y2 = self->kernel->view.v2j(Ve.y);
	if (Sqr((double)(i-x2)) + Sqr((double)(j-y2)) <= (double)Sqr(size+d))
		return 2;	// End point

	// Check on the line using the clipped coordinates
	x2 = self->kernel->view.u2i(Vec.x);
	y2 = self->kernel->view.v2j(Vec.y);

	double Dx = (double)(x2-x1);
	double Dy = (double)(y2-y1);
	double Dlen2 = Sqr(Dx) + Sqr(Dy);
	if (Dlen2 < (double)Sqr(size)) return -1;

	double xx = (double)(i-x1);
	double yy = (double)(j-y1);

	// Check distance, calculate the cross product (to the square)
	// of the two vs length
	if (Sqr(Dy*xx - Dx*yy) > Sqr((double)(size+d)) * Dlen2)
		return -1;

	// Calculate the projection distance using the dot product
	if (InRange(0.0, Dx*xx + Dy*yy, Dlen2))
		return 0;
	else
		return -1;
} // closest

/** enclosed */
bool GArrow::enclosed(ViewerObject *self, const int left, const int top, const int right, const int bottom)
{
	x1 = self->kernel->view.u2i(Vc.x);
	y1 = self->kernel->view.v2j(Vc.y);
	x2 = self->kernel->view.u2i(Vec.x);
	y2 = self->kernel->view.v2j(Vec.y);

	return clipSegment(&x1, &y1, &x2, &y2, left, top, right, bottom);
} // enclosed

/* ================================= GSpline ================================= */
/** set parameter
 * @param key	value to set
 * @param value	python object to set value from
 * @param return true on error and set the python error
 */
PyObject* GSpline::config(const char *key, PyObject *value)
{
	if (!strcmp(key, "pos")) {
		// Override default behaviour of GObject
	} else
	if (strchr("nxyz",key[0]) && strchr("0123456789",key[1])) {
		int idx = atoi(key+1);
		if (idx < MAX_NODES) {
			n = Max(n, idx+1);
			switch (key[0]) {
				case 'x':
					if (value)
						point[idx].x = Py_GetFloat(value);
					else
						return PyFloat_FromDouble(point[idx].x);
					break;
				case 'y':
					if (value)
						point[idx].y = Py_GetFloat(value);
					else
						return PyFloat_FromDouble(point[idx].y);
					break;
				case 'z':
					if (value)
						point[idx].z = Py_GetFloat(value);
					else
						return PyFloat_FromDouble(point[idx].z);
					break;
				case 'n':
					if (value)
						point[idx] = Py_GetVector(value);
					else
						return Py_Vector(point[idx]);
					break;
			}
		}
	} else
		return GObject::config(key, value);

	if (PyErr_Occurred()) return NULL;
	Py_RETURN_NONE;
} // config

/** transform */
void GSpline::transform(ViewerObject *self)
{
	GObject::transform(self);
	for (int i=0; i<n; i++) {
		self->kernel->view.xyz2uvw3D(point[i], &Vpoint[i]);
		// FIXME very bad!!!
		if (self->kernel->view.inside(Vpoint[i].x,Vpoint[i].y))
			_visible = true;
	}
} // transform

/** draw */
void GSpline::draw(ViewerObject *self, Drawable drawable)
{
	GObject::draw(self, drawable);

	for (int i=1; i<n; i++) {
		XDrawLine3D(self->display, drawable, self->gc, self->kernel->view, Vpoint[i-1], Vpoint[i]);
	}
} // draw

/* ================================= GRuler ================================== */
/** set parameter
 * @param key	value to set
 * @param value	python object to set value from
 * @param return true on error and set the python error
 */
PyObject* GRuler::config(const char *key, PyObject *value)
{
	if (!strcmp(key, "xa")) {
		if (value==NULL)
			return PyFloat_FromDouble(P.x+Da.x);
		else
			Da.x = Py_GetFloat(value) - P.x;
	} else
	if (!strcmp(key, "ya")) {
		if (value==NULL)
			return PyFloat_FromDouble(P.y+Da.y);
		else
			Da.y = Py_GetFloat(value) - P.y;
	} else
	if (!strcmp(key, "za")) {
		if (value==NULL)
			return PyFloat_FromDouble(P.z+Da.z);
		else
			Da.z = Py_GetFloat(value) - P.z;
	} else
	if (!strcmp(key, "ang") || (!strcmp(key, "pos") && selnode==3)) {
		// Override default behaviour of GObject
		if (value==NULL)
			return Py_Vector(P+Da);
		else {
			Da = Py_GetVector(value);
			Da -= P;
		}
	} else
	if (!strcmp(key, "dxa")) {
		if (value==NULL)
			return PyFloat_FromDouble(Da.x);
		else
			Da.x = Py_GetFloat(value);
	} else
	if (!strcmp(key, "dya")) {
		if (value==NULL)
			return PyFloat_FromDouble(Da.y);
		else
			Da.y = Py_GetFloat(value);
	} else
	if (!strcmp(key, "dza")) {
		if (value==NULL)
			return PyFloat_FromDouble(Da.z);
		else
			Da.z = Py_GetFloat(value);
	} else
	if (!strcmp(key, "dang")) {
		if (value==NULL)
			return Py_Vector(Da);
		else
			Da = Py_GetVector(value);
	} else
	if (!strcmp(key, "option")) {
		PyObject* obj = GObject::config(key, value);
		drawHead = true;
		drawTail = (option==Ruler_Simple);
		return obj;
	} else
		return GArrow::config(key, value);

	Py_RETURN_NONE;
} // config

/** move object item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 */
void GRuler::node(int item, const Vector& r)
{
	if (item==3)
		Da  = r - P;
	else
		GArrow::node(item,r);
} // node

/** @return node position */
Vector GRuler::node(int item) const
{
	if (item==3)
		return P + Da;
	else
		return GArrow::node(item);
} // node

/** @return savedNode position */
Vector GRuler::savedNode(int item) const
{
	if (item==3)
		return SP + SDa;
	else
		return GArrow::savedNode(item);
} // savedNode

/** transform */
void GRuler::transform(ViewerObject *self)
{
	GArrow::transform(self);
	if (option == Ruler_Angle) {
		self->kernel->view.xyz2uvw3D(P+Da, &Va);
		Vac = Va;
		Vc2 = V;
		_visible |= self->kernel->view.clipLine(&Vc2.x,&Vc2.y, &Vac.x,&Vac.y);
	}
} // transform

/** rotate object around an axis
 * @param angle	to rotate
 * @param axis	to rotate
 */
void GRuler::rotate(const double angle, const Vector& axis)
{
	Matrix4	rot;
	rot.rotate(angle, axis);
	D  = rot * SD;
	Da = rot * SDa;
} // rotate

/** save */
void GRuler::save()
{
	GArrow::save();
	SDa = Da;
} // save

/** restore */
void GRuler::restore()
{
	GArrow::restore();
	Da = SDa;
} // restore

/** @return bounding box of object */
BBox GRuler::bbox() const
{
	BBox bb;
	bb.add(P);
	bb.add(P+D);
	if (option == Ruler_Angle)
		bb.add(P+Da);
	return bb;
} // bbox

/** @return bounding box of object */
BBox GRuler::bboxView(ViewerObject *self)
{
	BBox bb;
	transform(self);
	bb.add(V);
	bb.add(Ve);
	if (option == Ruler_Angle)
		bb.add(Va);
	return bb;
} // bboxView

/** draw */
void GRuler::draw(ViewerObject *self, Drawable drawable)
{
	gcValues.line_style = LineSolid;
	GArrow::draw(self, drawable);

	if (option==Ruler_Angle) {
		// Using the clipped relative-coordinates
		double x1d = self->kernel->view.u2id(Vc2.x);
		double y1d = self->kernel->view.v2jd(Vc2.y);
		x1 = Round(x1d);
		y1 = Round(y1d);

		double x3d = self->kernel->view.u2id(Vac.x);
		double y3d = self->kernel->view.v2jd(Vac.y);
		x3 = Round(x3d);
		y3 = Round(y3d);

		double x2d = self->kernel->view.u2id(Vec.x);
		double y2d = self->kernel->view.v2jd(Vec.y);

		// Draw end with solid line
		XDrawRectangle(self->display, drawable, self->gc,
			x3-1, y3-1, 2, 2);

		// Calculate the angles to the reference axis and the end point
		double aa = atan2(y1d-y3d, x3d-x1d);
		double ae = atan2(y1d-y2d, x2d-x1d);
		double ang = ae-aa;
		//if (ang<0.0) ang += PI2;
		if (ang<0.0) {
			Swap(aa,ae);
			ang = -ang;
		}

		// Draw the angle arc
#define RADIUS	10
		XDrawArc(self->display, drawable, self->gc,
			x1-RADIUS, y1-RADIUS, 2*RADIUS, 2*RADIUS,
			Round(DEG(aa)) << 6, Round(DEG(ang)) << 6);

		// Draw selected point
		if (select && selnode==3 && self->kernel->view.inside(Va.x,Va.y)) {
			drawSelectedPoint(self, drawable,
				self->kernel->view.u2i(Va.x),
				self->kernel->view.v2j(Va.y));
			setForeground(self);
		}

		// And reference with dashed
		gcValues.line_width = 0;
		gcValues.line_style = LineOnOffDash;
		XChangeGC(self->display, self->gc, GCLineWidth|GCLineStyle, &gcValues);
		XSetDashes(self->display, self->gc, 0,
			dashedPattern, SIZE(dashedPattern));
		XDrawLine(self->display, drawable, self->gc, x1, y1, x3, y3);
	}
} // draw

/** drawText */
void GRuler::drawText(ViewerObject *self, Drawable drawable)
{
	char   txt[256];
	double Dlen = D.length();
	sprintf(txt, "%s: %.7g", _name, Dlen);

	// Using the clipped relative-coordinates
	x1 = self->kernel->view.u2i(Vc.x);
	y1 = self->kernel->view.v2j(Vc.y);
	x2 = self->kernel->view.u2i(Vec.x);
	y2 = self->kernel->view.v2j(Vec.y);
	int i = (x1+x2)/2;
	int j = (y1+y2)/2;

	int l = (int)strlen(txt);
//	int V.z = XTextWidth(font, txt, l);

	switch (anchor) {
		case Anchor_N:
			j -= size+1;
			break;
		case Anchor_NE:
			j -= size+1;
			i += size+1;
			break;
		case Anchor_E:
			if (x1>x2) {
				i = x1 + size+1;
				j = y1;
			} else {
				i = x2 + size+1;
				j = y2;
			}
			break;
		case Anchor_SE:
			j += size+1;
			i += size+1;
			break;
		case Anchor_S:
			j += size+1;
			break;
		case Anchor_SW:
			j += size+1;
			i -= size+1;
			break;
		case Anchor_W:
			if (x1<x2) {
				i = x1 - size+1;
				j = y1;
			} else {
				i = x2 - size+1;
				j = y2;
			}
			break;
		case Anchor_NW:
			j -= size+1;
			i -= size+1;
			break;
		default: ;	// just to keep compiler happy
	}
	XDrawString(self->display, drawable, self->gc, i, j, txt, l);

	if (option==Ruler_Angle) {
		double Dalen = Da.length();
		double ang;
		if (Eq0(Dlen,SMALL) || Eq0(Dalen,SMALL))
			ang = 0.0;
		else
			ang = DEG(acos((D*Da)/Dlen/Dalen));
		sprintf(txt,"%.5g deg",ang);
		i = x1;
		j = y1;
		l = (int)strlen(txt);
		switch (anchor) {
			case Anchor_N:
				j -= size+1;
				break;
			case Anchor_NE:
				i += size+1;
				j -= size+1;
				break;
			case Anchor_E:
				i = x1 + size+1;
				//j = y1;
				break;
			case Anchor_SE:
				i += size+1;
				j += size+1;
				break;
			case Anchor_S:
				j += size+1;
				break;
			case Anchor_SW:
				i -= size+1;
				j += size+1;
				break;
			case Anchor_W:
				i = x1 - size+1;
				//j = y1;
				break;
			case Anchor_NW:
				i -= size+1;
				j -= size+1;
				break;
			default: ;	// just to keep compiler happy
		}
		XDrawString(self->display, drawable, self->gc, i, j, txt, l);
	}
} // drawText

/** closest */
int GRuler::closest(ViewerObject *self, const int i, const int j, const int d)
{
	// First check angle point just in case is overlapped with the main one
	if (option==Ruler_Angle) {
		// Check angle point
		x3 = self->kernel->view.u2i(Va.x);
		y3 = self->kernel->view.v2j(Va.y);
		if (Sqr((double)(i-x3)) + Sqr((double)(j-y3)) <= (double)Sqr(size+d))
			return 3;	// Angle point
	}

	// Check the arrow
	int close = GArrow::closest(self, i, j, d);
	if (close>=0) return close;

	// Check the arrow line
	if (option==Ruler_Angle) {
		// Check angle line
		x1 = self->kernel->view.u2i(Vc2.x);
		y1 = self->kernel->view.v2j(Vc2.y);
		x2 = self->kernel->view.u2i(Vac.x);
		y2 = self->kernel->view.v2j(Vac.y);

		double Dx = (double)(x2-x1);
		double Dy = (double)(y2-y1);
		double Dlen2 = Sqr(Dx) + Sqr(Dy);
		if (Dlen2 < (double)Sqr(size)) return -1;
		double xx = (double)(i-x1);
		double yy = (double)(j-y1);

		// Check distance, calculate the cross product (to the square)
		// of the two vs length
		if (Sqr(Dy*xx - Dx*yy) > Sqr((double)(size+d)) * Dlen2)
			return -1;

		// Calculate the projection distance using the dot product
		if (InRange(0.0, Dx*xx + Dy*yy, Dlen2))
			return 0;
	}
	return -1;
} // closest

/** enclosed */
bool GRuler::enclosed(ViewerObject *self, const int left, const int top, const int right, const int bottom)
{
	if (GArrow::enclosed(self, left, top, right, bottom))
		return true;

	if (option==Ruler_Angle) {
		x1 = self->kernel->view.u2i(Vc2.x);
		y1 = self->kernel->view.v2j(Vc2.y);
		x2 = self->kernel->view.u2i(Vac.x);
		y2 = self->kernel->view.v2j(Vac.y);

		return clipSegment(&x1, &y1, &x2, &y2, left, top, right, bottom);
	} else
		return false;
} // enclosed

/* ================================= GRotdefi ================================== */
/** set parameter
 * @param key	value to set
 * @param value	python object to set value from
 * @param return true on error and set the python error
 */
PyObject* GRotdefi::config(const char *key, PyObject *value)
{
	if (!strcmp(key, "axissize")) {
		if (value==NULL)
			return PyInt_FromLong((long)axisSize);
		else {
			axisSize = Py_GetInt(value);
			if (axisSize<=0) axisSize = 20;
		}
	} else
	if (!strcmp(key, "axiswidth")) {
		if (value==NULL)
			return PyInt_FromLong((long)axisWidth);
		else {
			axisWidth = Py_GetInt(value);
			if (axisWidth<=0) axisWidth = 2;
		}
	} else
	if (!strcmp(key, "orient")) {
		if (value!=NULL) {
			if (!PyList_AsMatrix4(value, orient)) return NULL;
		} else
			return PyList_FromMatrix4(orient);
	} else
	if (!strcmp(key, "matrix")) {
		if (value!=NULL) {
			if (!PyList_AsMatrix4(value, matrix)) return NULL;
			D = matrix*P - P;
		} else
			return PyList_FromMatrix4(matrix);
	} else
		return GArrow::config(key, value);

	if (PyErr_Occurred()) return NULL;
	Py_RETURN_NONE;
} // config

/** draw */
void GRotdefi::draw(ViewerObject *self, Drawable drawable)
{
	D = matrix*P - P;
	GArrow::draw(self, drawable);

	if (axisWidth != lineWidth) {
		gcValues.line_width = axisWidth;
		XChangeGC(self->display, self->gc, GCLineWidth, &gcValues);
	}

	if (self->kernel->view.inside(Ve.x,Ve.y)) {
		Matrix4 A;
		A.multiply(self->kernel->view.matrix(), orient);
		XDrawAxes(self->display, drawable, self->gc, x1, y1, axisSize, A);
	}

	if (self->kernel->view.inside(V.x,V.y)) {
		Matrix4 A;
		A.multiply(self->kernel->view.matrix(), matrix);
		XDrawAxes(self->display, drawable, self->gc, x2, y2, axisSize, A, false);
	}
} // draw

/** enclosed */
bool GRotdefi::enclosed(ViewerObject * /*self*/, const int /*left*/, const int /*top*/, const int /*right*/, const int /*bottom*/)
{
/*
	int x1 = self->kernel->view.u2i(V.x);
	int y1 = self->kernel->view.v2j(V.y);

	int xl = Max(x1-size, left);
	int yl = Max(y1-size, top);
	int xh = Min(x1+size, right);
	int yh = Min(y1+size, bottom);

	return (xl<=xh && yl<=yh);
*/
	return false;
} // enclosed

/* ================================= GLight ================================= */
GLight::GLight(const char *aname, ObjectType o)
	: GArrow(aname, o)
{
	power    =  1.0;
	spec     =  0.0;
	relative = false;
	shadow   = false;
	falloff  = 0;
} // GLight

PyObject* GLight::config(const char *key, PyObject *value)
{
	if (!strcmp(key, "pos")) {	// Override default behaviour of GArrow
		if (selnode<=1)
			return GObject::config(key,value);
		else {
			if (value==NULL)
				return Py_Vector(P+D);
			else {
				D = Py_GetVector(value);
				D -= P;
			}
		}
	} else
	if (!strcmp(key, "dx")) {
		if (value==NULL)
			return PyFloat_FromDouble(D.x);
		else
			D.x = Py_GetFloat(value);
	} else
	if (!strcmp(key, "dy")) {
		if (value==NULL)
			return PyFloat_FromDouble(D.y);
		else
			D.y = Py_GetFloat(value);
	} else
	if (!strcmp(key, "dz")) {
		if (value==NULL)
			return PyFloat_FromDouble(D.z);
		else
			D.z = Py_GetFloat(value);
	} else
	if (!strcmp(key, "relative")) {
		if (value==NULL)
			if (PyBool_FromLong(relative))
				return PyString_FromString("on");
			else
				return PyString_FromString("off");
		else
			relative = Py_GetBool(value);
	} else
	if (!strcmp(key, "power")) {
		if (value==NULL)
			return PyFloat_FromDouble(power);
		else
			power = Range(0.0, Py_GetFloat(value), 100000.0);
	} else
	if (!strcmp(key, "falloff")) {
		if (value==NULL)
			return PyInt_FromLong(falloff);
		else
		if (PyInt_Check(value)) {
			falloff = PyInt_AsLong(value);
			if (falloff > 2) {
				falloff -= 3;
				shadow   = true;
			} else
				shadow   = false;
		}
	} else
	if (!strcmp(key, "specular")) {
		if (value==NULL)
			return PyFloat_FromDouble(spec);
		else
			spec = Range(0.0, Py_GetFloat(value), 100000.0);
	} else
	if (!strcmp(key, "shadow")) {
		if (value==NULL)
			if (PyBool_FromLong(shadow))
				return PyString_FromString("on");
			else
				return PyString_FromString("off");
		else
			shadow = Py_GetBool(value);
	} else
		return GArrow::config(key, value);

	Py_RETURN_NONE;
} // config

/** move object item to location
 * @param item	item to move as returned by closest
 * @param r	position vector
 */
void GLight::node(int item, const Vector& r)
{
	if (item <= 1)
		P = r;
	else
		GArrow::node(item,r);
} // node

/** transform *
void GLight::transform(ViewerObject *self)
{
	if (relative) {
		V.x  = x;
		V.y  = y;
		V.z  = z;
		du = D.x;
		dv = D.y;
		dw = D.z;
	} else {
		GObject::transform(self);
		self->kernel->view.dxyz2duvw3D(D.x, D.y, D.z, &du, &dv, &dw);
	}

	_visible = self->kernel->view.inside(V.x,V.y);
} * transform */

/** draw */
void GLight::draw(ViewerObject *self, Drawable drawable)
{
	int size2;

	if (option != Light_Omni)
		GArrow::draw(self, drawable);
	else {
		GObject::draw(self, drawable);
		// unfortunatelly the point is on different plane from the maximum circle
		//if (select && selnode==2 && self->kernel->view.inside(Ve.x,Ve.y)) {
		//	drawSelectedPoint(self, drawable,
		//		self->kernel->view.u2i(Ve.x),
		//		self->kernel->view.v2j(Ve.y));
		//	setForeground(self);
		//}
	}

	size2 = size*2 + 1;
	XDrawArc(self->display, drawable, self->gc, x1-size, y1-size, size2, size2, 0, 360*64);
	gcValues.line_style = LineOnOffDash;
	XChangeGC(self->display, self->gc, GCLineStyle, &gcValues);
	XSetDashes(self->display, self->gc, 0, dashedPattern, SIZE(dashedPattern));
	XDrawArc(self->display, drawable, self->gc, x1-size-2, y1-size-2, size2+4, size2+4, 0, 360*64);

	if (option==Light_Omni && select) {
		double R = D.length() * self->kernel->view.Sx();
		if (R<1e4) {
			int r = (int)R;
			gcValues.line_style = LineSolid;
			XChangeGC(self->display, self->gc, GCLineStyle, &gcValues);
			XDrawArc(self->display, drawable, self->gc, x1-r, y1-r, 2*r, 2*r, 0, 360*64);
		}
	}
} // draw

/** drawText */
void GLight::drawText(ViewerObject *self, Drawable drawable)
{
	if (option != Light_Omni)
		GArrow::drawText(self, drawable);
	else
		GObject::drawText(self, drawable);
} // drawText

/** closest */
int GLight::closest(ViewerObject *self, const int i, const int j, const int d)
{
	int isclose = GArrow::closest(self,i,j,d);
	if (option == Light_Omni) {
		double R = D.length() * self->kernel->view.Sx();
		if (R<1e9 && Abs(sqrt(Sqr(i-x1) + Sqr(j-y1)) - (int)R) < d)
				return 2;
	}
	return isclose;
} // closest

/** complete light structure */
void GLight::toLight(Light *L) const
{
	switch (option) {
		case Light_Sun:
			L->type = LIGHT_SUN;
			break;
		case Light_Omni:
			L->type = LIGHT_OMNI;
			break;
		case Light_Spot:
			L->type = LIGHT_SPOT;
			break;
		default:
			L->type = LIGHT_SUN;
	}
	L->pos      = P;
	L->dir      = D;
	L->relative = relative;
	L->power    = power;
	L->falloff  = falloff;
	L->spec     = spec;
	L->shadow   = shadow;
} // GLight

/* ================================= GBeam ================================= */
GBeam::GBeam(const char *aname, ObjectType o)
	: GArrow(aname, o)
{
	energy     = 1.0;
	scale      = 1.0;
	divergence = 0.0;
	Rin = Rout = 0.0;
	size = 10;
} // GBeam

PyObject* GBeam::config(const char *key, PyObject *value)
{
	if (!strcmp(key, "energy")) {
		if (value==NULL)
			return PyFloat_FromDouble((double)energy);
		else {
			energy = Abs(Py_GetFloat(value));
			if (energy<=0.0) energy=1.0;
		}
	} else
	if (!strcmp(key, "scale")) {
		if (value==NULL)
			return PyFloat_FromDouble((double)scale);
		else {
			scale = Abs(Py_GetFloat(value));
			if (scale<=0.0) scale=1.0;
		}
	} else
	if (!strcmp(key, "divergence")) {
		if (value==NULL)
			return PyFloat_FromDouble((double)divergence);
		else
			divergence = Abs(Py_GetFloat(value));
	} else
	if (!strcmp(key, "dx")) {
		if (value==NULL)
			return PyFloat_FromDouble(D.x/(energy*scale));
		else
			D.x = Py_GetFloat(value)*energy*scale;
	} else
	if (!strcmp(key, "dy")) {
		if (value==NULL)
			return PyFloat_FromDouble(D.y/(energy*scale));
		else
			D.y = Py_GetFloat(value)*energy*scale;
	} else
	if (!strcmp(key, "dz")) {
		if (value==NULL)
			return PyFloat_FromDouble(D.z/(energy*scale));
		else
			D.z = Py_GetFloat(value)*energy*scale;
	} else
		return GArrow::config(key, value);

	Py_RETURN_NONE;
} // config

/** draw */
void GBeam::draw(ViewerObject *self, Drawable drawable)
{
	GArrow::draw(self, drawable);
	//GObject::drawText(self, drawable);

	//if (option==Beam_Divergence && divergence>1E-6) {
	if (divergence>1E-6) {
		int dx = x2-x1;
		int dy = y2-y1;	// up-side down screen
		int length = (int)(sqrt(dx*dx + dy*dy)*0.8);
		double angle = atan2(dy, dx);

		if (divergence<PI2) {
			XPoint	pts[3];
			double s, c;
			sincos(angle-divergence, &s, &c);
			pts[0].x = x1+(int)(c*length);
			pts[0].y = y1+(int)(s*length);

			pts[1].x = x1;
			pts[1].y = y1;

			sincos(angle+divergence, &s, &c);
			pts[2].x = x1+(int)(c*length);
			pts[2].y = y1+(int)(s*length);

			XDrawLines(self->display, drawable, self->gc,
					pts, SIZE(pts), CoordModeOrigin);
		}
		int angle1 = (int)(DEG(-angle-divergence)*64.0);
		int angle2 = (int)(DEG(divergence)*2.0*64.0);

		XDrawArc(self->display, drawable, self->gc,
				x1-length, y1-length, 2*length, 2*length,
				angle1, angle2);
	}
	GObject::drawText(self, drawable);
} // draw
