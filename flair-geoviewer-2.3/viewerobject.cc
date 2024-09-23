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
 * Date:	01-Apr-2010
 * Last change: 27-May-2022 by Paola Sala
 */

#include <Python.h>	// it has to be the first include!

#include <tk.h>
#include <tcl.h>
#include <math.h>
#include <iosfwd>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <tkDecls.h>
#include <fnmatch.h>
#include <X11/Xutil.h>

#if 0 //def OPENGL
#	include <GL/gl.h>
#	include <GL/glx.h>
#	include <GL/glu.h>
#endif

#include "os.h"
#include "geo.h"
#include "array.h"
#include "gzone.h"
#include "timer.h"
#include "vzone.h"
#include "xdraw.h"
#include "gobject.h"
#include "pyutils.h"
#include "viewerobject.h"

using namespace std;

#define VPMARGIN	1	// Margin from the border of viewport lines

static const char dashedPattern[2] = {5,5};

Py_ssize_t *lllen;

#ifdef OPENGL
static int glBuffer[] = {
		GLX_RGBA,
		GLX_RED_SIZE,    8,
		GLX_GREEN_SIZE,  8,
		GLX_BLUE_SIZE,   8,
		GLX_DEPTH_SIZE, 12,
//		GLX_DOUBLEBUFFER,
		None };
#endif

/* --- Clear viewport errors --- */
static void _clearErrors(ViewerObject *self)
{
	if (self->errors->empty()) return;
	ArrayIterator<GObject*> iter(*self->errors);
	while (iter) {
		GObject *object = iter++;
		delete object;
	}
	self->errors->clear();
} // _clearErrors

/** sendMessage: send message to tcl interpreter */
static void _sendMessage(ViewerObject* self, const char* msg)
{
	union {XEvent gen; XVirtualEvent virt;} event;

	if (!self->threaded) return;
	if (!Tk_IsMapped(self->tkWin)) return;

	memset(&event, 0, sizeof(event));
	event.gen.xany.type       = VirtualEvent;
	event.gen.xany.serial     = NextRequest(self->display);
	event.gen.xany.send_event = False;
	event.gen.xany.window     = self->window;
	event.gen.xany.display    = self->display;

	event.virt.name           = Tk_GetUid(msg);
	Tk_QueueWindowEvent(&event.gen, TCL_QUEUE_TAIL);
} // _sendMessage

/* --- loadFont, tga or system --- */
static bool loadFont(BFont& font, const char *fontpath)
{
	if (font.load(fontpath)) return true;
	PyErr_Format(PyExc_TypeError,
		"Invalid font \'%s\' specified", fontpath);
	return false;

#if 0
	// Cause a seg-fault after few font reloads with bad malloc!

	int screen = DefaultScreen(display);

	// Load bitmap system font
	XFontStruct* fontinfo = XLoadQueryFont(display, fontpath);
	if (fontinfo == NULL) {
		PyErr_Format(PyExc_TypeError,
			"Invalid system font \'%s\' specified", fontpath);
		return false;
	}
	int width  = fontinfo->max_bounds.width;
	int height = fontinfo->max_bounds.ascent + fontinfo->max_bounds.descent;

	Pixmap pixmap = XCreatePixmap(display,
				XRootWindow(display,screen),
				16*width, 16*height,
				DefaultDepth(display,screen));

	XGCValues values;
	values.font       = fontinfo->fid;
	values.function   = GXcopy;
	values.plane_mask = AllPlanes;
	values.foreground = BlackPixel(display,screen);
	values.background = BlackPixel(display,screen);

	GC gc = XCreateGC(display,
		pixmap,
		GCFont |
		GCFunction | GCPlaneMask | GCForeground | GCBackground,
		&values);

	XFillRectangle(display,pixmap,gc,0,0,16*width,16*height);
	XSetForeground(display, gc, WhitePixel(display,screen));

	for (int row=0; row<16; row++) {
		int y = row*height;
		for (int col=0; col<16; col++) {
			int x = col*width;
			char str[2];
			str[0] = row<<4 | col;
			str[1] = 0;
			XDrawImageString(display,pixmap,gc,x,y,str,1);
		}
	}
	XFreeFont(display, fontinfo);

	XImage* image = XGetImage(display, pixmap,
		0, 0, 16*width, 16*height, AllPlanes,ZPixmap);

	if (image) {
		font.setMap(fontpath, 16*width, 16*height, (dword*)image->data);
		XDestroyImage(image);
	} else {
		PyErr_Format(PyExc_TypeError,
			"Invalid system font \'%s\' specified", fontpath);
		return false;
	}
	XFreeGC(display, gc);
	XFreePixmap(display, pixmap);
	return true;
#endif
} // loadFont

/* Py_VBody
 * @param self	Viewer object
 * @param obj	object to convert to body can be string or id
 * @return body
 */
static VBody *Py_VBody(ViewerObject *self, PyObject *obj)
{
	VBody *body;
	if (obj==NULL) return NULL;
	//	if (PyString_Check(obj)) {
	if (PyString_Check(obj)) {
	  //		char *name = PyString_AsString(obj);
	  const char *name = PyUnicode_AsUTF8AndSize(obj, lllen);
		body = self->kernel->getBody(name);
		if (body==NULL) {
			PyErr_Format(PyExc_KeyError, "Body \'%s\' not found", name);
			return NULL;
		}
		return body;
	} else
	if (PyInt_Check(obj)) {
		int id = (int)PyInt_AsLong(obj);
		if (id>=0 && id<self->kernel->bodies.count()) {
			body = self->kernel->getBody(id);
			if (body==NULL) {
				PyErr_Format(PyExc_IndexError, "Body #%d not found",id);
				return NULL;
			}
			return body;
		}
		PyErr_Format(PyExc_IndexError, "Body #%d not found",id);
		return NULL;
	} else {
		PyErr_SetString(PyExc_TypeError, "Invalid body type, string or integer expected");
		return NULL;
	}
} // Py_VBody

/* Py_VRegion
 * @param self	Viewer object
 * @param obj	object to convert to region can be string or id
 * @return region
 */
static VRegion *Py_VRegion(ViewerObject *self, PyObject *obj)
{
	VRegion *region;
	if (obj==NULL) return NULL;
	//	if (PyString_Check(obj)) {
	  //		char *name = PyString_AsString(obj);
	if (PyUnicode_Check(obj)) {
	const char *name = PyUnicode_AsUTF8AndSize(obj, lllen);
		// FIXME: potential source of troubles
		region = self->kernel->getRegion(name);
		if (region==NULL) {
			PyErr_Format(PyExc_KeyError, "Region \'%s\' not found", name);
			return NULL;
		}
		return region;
	} else
	if (PyInt_Check(obj)) {
		int id = (int)PyInt_AsLong(obj);
		// FIXME: potential source of troubles
		region = self->kernel->getRegion(id);
		if (region==NULL) {
			PyErr_Format(PyExc_IndexError, "Region #%d not found",id);
			return NULL;
		}
		return region;
	} else {
		PyErr_SetString(PyExc_TypeError, "Invalid object type, string or integer expected");
		return NULL;
	}
} // Py_VRegion

//////////////////////// Viewer /////////////////////////////
// Static prototypes
static void Viewer_calcViewport(ViewerObject *self, int id);

/* --- Viewer_new --- */
static PyObject* Viewer_new(PyTypeObject *type, PyObject* /*args*/, PyObject* /*kwds*/)
{
	ViewerObject *self;

	self = (ViewerObject*)type->tp_alloc(type, 0);
	if (self == NULL) return NULL;

	self->display  = NULL;
	self->tkWin    = NULL;
	self->window   = 0;
	self->depth    = 0;
	self->gc       = NULL;
	self->ximage   = NULL;
	self->pixmap   = 0;
#ifdef OPENGL
	self->visinfo  = NULL;
#endif
	self->viewer   = NULL;
	self->geometry = NULL;
	self->showViewport  = true;
	self->crosshair     = 0;
	self->showTrackball = false;
	self->projectionChanged = true;
	self->errors   = new Array<GObject*>(16);
	self->pen      = new Array<IPoint>(16);

	for (int i=0; i<NVIEWS; i++) {
		self->viewport[i].viewer      = NULL;
		self->viewport[i].visible     = false;
		self->viewport[i].lineWidth   = 1;
		self->viewport[i].originWidth = 2;
		self->viewport[i].wLength     = 15.0;
		self->viewport[i].color       = 0xFFFF00;
	}

	return (PyObject*)self;
} // Viewer_new

/* --- Viewer_init --- */
static int Viewer_init(ViewerObject *self, PyObject *args, PyObject* /*kwds*/)
{
	PyObject *interpaddr;
	char  *winname;
	GeometryObject* geometry;

	if (!PyArg_ParseTuple(args, "OsO", &geometry, &winname, &interpaddr)) return -1;
	if (Py_TYPE(geometry) != &GeometryType) {
		PyErr_SetString(PyExc_TypeError, "Invalid type, Geometry type expected");
		return -1;
	}

	self->geometry = geometry;
	self->kernel = new GeometryKernel(*geometry->geometry);
	self->viewer = new GeometryViewer(*geometry->geometry, *self->kernel);

	Tcl_Interp *interp = (Tcl_Interp*)PyLong_AsVoidPtr(interpaddr);

	// Find if tcl is compiled with threads enabled
	int rc = Tcl_Eval(interp, "::tcl::pkgconfig get threaded");
	if (rc==0)
		self->threaded = !strcmp(Tcl_GetStringResult(interp), "1");
	else
		self->threaded = false;

	// Find Window
	self->tkWin = Tk_NameToWindow(interp, winname, Tk_MainWindow(interp));
	if (!self->tkWin) {
		PyErr_SetString(PyExc_ValueError, Tcl_GetStringResult(interp));
		return -1;
	}

	// setup display
	if (self->display == NULL) {
		self->display = Tk_Display(self->tkWin);
		self->depth   = DefaultDepth(self->display,
				DefaultScreen(self->display));
	}

	// WARNING??? window at this location is not set correctly!!!!
	assert(self->gc==NULL);
//	self->window = Tk_WindowId(self->tkWin);
//	if (self->gc) {
//		XFreeGC(self->display, self->gc);
//		self->gc = NULL;
//	}

	assert(self->ximage==NULL);
	//if (!self->ximage)
	self->ximage = XCreateImage(self->display,
		DefaultVisual(self->display, DefaultScreen(self->display)),
		self->depth, ZPixmap, 0,
		(char*)(self->viewer->painter.data()),
		self->viewer->width(),
		self->viewer->height(),
		32, 0);

	// Initialize variables
	self->startU = self->startV = 0.0;
	self->pivotU = self->pivotV = 0.0;
	self->pivotAngle = INFINITE;	// a large number not to show rotation

	self->rectX1 = -1;	// Don't draw
	self->rectY1 = -1;
	self->rectX2 = -1;
	self->rectY2 = -1;

	return 0;
} // Viewer_init

/* --- Viewer_dealloc --- */
static void Viewer_dealloc(ViewerObject *self)
{
	if (self->viewer) {
		self->viewer->stopThread();
		// will be deleted from XDestroyImage
		self->viewer->painter.dataNull();
		delete self->viewer;
		delete self->kernel;
	}
	_clearErrors(self);
	delete self->errors;
	delete self->pen;
#ifdef OPENGL
	glXDestroyGLXPixmap(self->display, self->glxpixmap);
	glXDestroyContext(self->display, self->context);
#endif
	XFreePixmap(self->display, self->pixmap);
	if (self->ximage) XDestroyImage(self->ximage);
	if (self->gc) XFreeGC(self->display, self->gc);
	//	self->ob_type->tp_free((PyObject*)self);
	Py_TYPE(self)->tp_free((PyObject*)self);
} // Viewer_dealloc

#ifdef MEM
/* --- Viewer_destroy --- */
static PyObject* Viewer_destroy(ViewerObject *self)
{
	// unfortunately we cannot call _dealloc due to all X, Tk deallocations
	delete self->viewer;
	delete self->kernel;
	_clearErrors(self);
	delete self->errors;
	delete self->pen;
	Py_RETURN_NONE;
} // Geometry_destroy
#endif

/* --- Viewer_configure --- */
static PyObject* Viewer_configure(ViewerObject *self, PyObject *args)
{
	int width, height;
	if (!PyArg_ParseTuple(args, "ii", &width, &height)) return NULL;
	if (width!=self->viewer->width() || height!=self->viewer->height()) {
		if (self->pixmap) {
			XFreePixmap(self->display, self->pixmap);
			self->pixmap = XCreatePixmap(self->display, self->window, width, height, self->depth);
#ifdef OPENGL
			glXDestroyGLXPixmap(self->display, self->glxpixmap);
			self->glxpixmap = glXCreateGLXPixmap(self->display, self->visinfo, self->pixmap);
#endif
		}

		self->viewer->resize(width, height);
		self->ximage->width  = width;
		self->ximage->height = height;
		self->ximage->data   = (char*)self->viewer->painter.data();
		self->ximage->bytes_per_line = 0;
		XInitImage(self->ximage);
	}
	Py_RETURN_NONE;
} // Viewer_configure

/* --- Viewer_expose ---
 * X11 main drawing routine.
 * Combines in a temporary pixmap
 * 1. the pixel painter that has to be already present using the GeometryViewer::draw() method
 * 2. superimpose x11 objects
 *
 */
static PyObject* Viewer_expose(ViewerObject *self, PyObject * /*args*/)
{
	XGCValues gcValues;

	int width   = Tk_Width(self->tkWin);
	int height  = Tk_Height(self->tkWin);
	int width2  = width / 2;
	int height2 = height / 2;

	// Initialize the first time
	if (self->gc == NULL) {
		if (!Tk_IsMapped(self->tkWin)) {
			Tk_MakeWindowExist(self->tkWin);
			if (!Tk_IsMapped(self->tkWin))
				Py_RETURN_NONE;
		}
		self->window = Tk_WindowId(self->tkWin);
		gcValues.function = GXcopy;
		gcValues.plane_mask = AllPlanes;
		//self->gc =  Tk_GetGC(self->tkWin, GCFunction|GCGraphicsExposures, &gcValues);
		self->gc = XCreateGC(self->display, self->window,
			GCFunction | GCPlaneMask | GCForeground | GCBackground,
			&gcValues);

		// Create initial pixmap
		self->pixmap = XCreatePixmap(self->display, self->window, width, height, self->depth);

#if OPENGL
		self->visinfo = glXChooseVisual(self->display, 0 /*self->screen*/, glBuffer);
		// Last argument is Direct rendering through False=through X11
		// Newer X-servers disable the indirect rendering!
		self->context   = glXCreateContext(self->display, self->visinfo, None, True);
		self->glxpixmap = glXCreateGLXPixmap(self->display, self->visinfo, self->pixmap);
#endif
	}

	// ---------- Draw image ----------
	//Pixmap pixmap = XCreatePixmap(self->display, self->window, width, height, self->depth);

	// -------- Initialize OpenGL pixmap --------
#if OPENGL
	glXMakeCurrent(self->display, self->glxpixmap, self->context);
	glViewport(0, 0, width, height);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
#endif

	XPutImage(self->display, self->pixmap, self->gc, self->ximage,
		0, 0, 0, 0, self->viewer->width(), self->viewer->height());

	// ---------- Do openGL stuff ---
#if OPENGL
	glEnable(GL_DEPTH_TEST);

	glShadeModel(GL_SMOOTH);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	/* light */
	static GLfloat light_pos[4] = {5.0, 20.0, 10.0, 1.0 };
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
	glEnable(GL_LIGHT0);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	/* Setup projection */
	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();

	float ratio = (float)width / (float)height;
	glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, -1.f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glBegin(GL_TRIANGLES);
	glColor3f(1.f, 0.f, 0.f);
	glVertex3f(-0.6f, -0.4f, 0.f);
	glColor3f(0.f, 1.f, 0.f);
	glVertex3f(0.6f, -0.4f, 0.f);
	glColor3f(0.f, 0.f, 1.f);
	glVertex3f(0.f, 0.6f, 0.f);
	glEnd();

//	float dw = self->viewer->view().imageWidth();
//	float dh = self->viewer->view().imageHeight();
//	glOrtho(-dw/2., dw/2., -dh/2, dh/2, 0.f, 1000.f);
//
//	double matrix[4][4];
//	memset(matrix, 0, sizeof(matrix));
//	for (int j=0; j<3; j++)
//		for (int i=0; i<3; i++)
//			matrix[j][i] = view.matrix(j,i);
//	matrix[3][3] = 1.0;
//	glMultMatrixd((GLdouble*)&matrix);
//	double x,y,z;
//	view.origin(&x, &y, &z);
//	glTranslated(-x,-y,-z);

//	/* Setup geometry */
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//
//	// Scan geometry and plot all meshes
//	for (int i=2; i<nox.geometry.bodies.size(); i++) {
//		GBody* body = nox.geometry.bodies[i];
//		cout << *body << endl;
//		GLMesh glmesh(&body->mesh);
//		glmesh.draw(GL_RENDER);
//	}

//	glFlush();	// flush contents. Doesn't ensure that the graph is completed
	glFinish();	// wait to finish and then flush
//	glXWaitGL();	// wait X11 to finish
//	glXSwapBuffers(display, glxpixmap);

#endif

	// ---------- Draw objects ----------
	if (self->geometry->objectList->count()) {
		ArrayIterator<GObject*> iter(*self->geometry->objectList);
		while (iter) {
			GObject *object = iter++;
			if (!object->show) continue;
			object->transform(self);
			if (object->visible())
				object->draw(self, self->pixmap);
		}
	}

	// ---------- Draw errors ----------
	if (self->errors->count() && self->viewer->d2.showBorders) {
		ArrayIterator<GObject*> iter(*self->errors);
		while (iter) {
			GObject *object = iter++;
			object->transform(self);
			if (object->visible())
				object->draw(self, self->pixmap);
		}
	}

	// ---------- Show crosshair (could be an object) FIXME XXX ----------
	if (self->crosshair) {
		gcValues.function   = GXxor;
		gcValues.line_width = 0;
		XSetForeground(self->display, self->gc, 0xFFFFFF);
		gcValues.line_style = LineSolid;
		XChangeGC(self->display, self->gc,
			GCFunction|GCLineWidth|GCLineStyle, &gcValues);

		XDrawLine(self->display, self->pixmap, self->gc,
			width2-self->crosshair, height2, width2+self->crosshair+1, height2);
		XDrawLine(self->display, self->pixmap, self->gc,
			width2, height2-self->crosshair, width2, height2+self->crosshair+1);
	}

	// ---------- Show trackball ----------
	if (self->showTrackball) {
		int radius = Round(Min(width, height) * self->geometry->trackballSize)/2;
		gcValues.function   = GXcopy;
		gcValues.line_width = 2;
		gcValues.line_style = LineSolid;
		XChangeGC(self->display, self->gc, GCFunction|GCLineWidth|GCLineStyle, &gcValues);
		XDrawTrackball(self->display, self->pixmap, self->gc, width2, height2,
				radius, self->kernel->view.matrix());
	}

	// ---------- Show editing body ----------
	const GBody *editBody = self->geometry->geometry->editBody();
	if (self->geometry->cursorShow) {
		int x,y;
		self->kernel->view.xyz2ij3D(self->geometry->cursor, &x, &y);
		int m = self->geometry->cursorMoveSize;
		int r = self->geometry->cursorSize;
		if (x>=-r && x<width+r && y>=-r && y<height+r) {
			// Draw the handlers
			gcValues.function   = GXxor;
			gcValues.line_width = 2;
			gcValues.line_style = LineSolid;
			XChangeGC(self->display, self->gc,
				GCFunction|GCLineWidth|GCLineStyle, &gcValues);

			XSetForeground(self->display, self->gc, 0xFFFFFF);
			XDrawArc(self->display, self->pixmap, self->gc, x-m, y-m, 2*m, 2*m, 0, 360*64);

			XSetForeground(self->display, self->gc, 0xAF00AF);
			XDrawArc(self->display, self->pixmap, self->gc, x-r, y-r, 2*r, 2*r, 0, 360*64);

			gcValues.function   = GXcopy;
			XChangeGC(self->display, self->gc, GCFunction, &gcValues);
			XDrawAxes(self->display, self->pixmap, self->gc, x, y, r, self->kernel->view.matrix());
		}
	}

	if (editBody && self->kernel->getBody(editBody)->visible) {
		// draw editing rulers if any
		if (editBody->show & BIT_SELECT && self->geometry->rulers) {
			// FIXME set position???? should be done when editing the body!
			for (int i=0; i<self->geometry->rulers; i++) {
				GRuler* ruler = self->geometry->ruler[i];
				ruler->position(editBody->position());
				char *name;
				switch (i) {
					case 0:
						name = editBody->showX();
						if (!name) continue;
						ruler->name(name);
						ruler->direction(editBody->vectorXlen() * editBody->vectorX());
						break;
					case 1:
						name = editBody->showY();
						if (!name) continue;
						ruler->name(name);
						ruler->direction(editBody->vectorYlen() * editBody->vectorY());
						break;
					case 2:
						name = editBody->showZ();
						if (!name) continue;
						ruler->name(name);
						ruler->direction(editBody->vectorZlen() * editBody->vectorZ());
						break;
					default:
						continue;
				}
				ruler->transform(self);
				if (ruler->visible())
					ruler->draw(self, self->pixmap);
			}
		}
	}

	// ---------- Show viewports ----------
	if (self->showViewport) {
		XSetDashes(self->display, self->gc, 0, dashedPattern, SIZE(dashedPattern));

		for (int i=0; i<NVIEWS; i++) {
			ViewportLine& vp = self->viewport[i];
			Viewer_calcViewport(self, i);

			gcValues.function   = GXcopy;
			gcValues.line_width = vp.lineWidth;
			gcValues.line_style = LineOnOffDash;
			XChangeGC(self->display, self->gc, GCFunction|GCLineStyle|GCLineWidth, &gcValues);

			if (vp.visible) {
				XSetForeground(self->display, self->gc, vp.color);
				// Draw line
				XDrawLine(self->display, self->pixmap, self->gc,
					vp.x[0], vp.y[0], vp.x[1], vp.y[1]);

				gcValues.line_style = LineSolid;
				XChangeGC(self->display, self->gc, GCLineStyle, &gcValues);

				// Draw Cursor
				XDrawRectangle(self->display, self->pixmap, self->gc,
					vp.xc-vp.originWidth, vp.yc-vp.originWidth,
					2*vp.originWidth, 2*vp.originWidth);

				// Draw w-vector
				XDrawLine(self->display, self->pixmap, self->gc,
					vp.xc, vp.yc, vp.xw, vp.yw);
			}
		}
	}

	// ---------- Show selection rectangle ----------
	if (self->rectX1 >= 0) {
		gcValues.function   = GXxor;
		gcValues.line_width = 0;
		gcValues.line_style = LineSolid;
		XChangeGC(self->display, self->gc, GCFunction | GCLineWidth | GCLineStyle, &gcValues);
		XSetForeground(self->display, self->gc, 0xFFFFFF);
		int x = Min(self->rectX1, self->rectX2);
		int y = Min(self->rectY1, self->rectY2);
		int w = Abs(self->rectX2 - self->rectX1);
		int h = Abs(self->rectY2 - self->rectY1);
		XDrawRectangle(self->display, self->pixmap, self->gc,
			x, y, w, h);
	}

	// ---------- Show hand drawing ----------
	if (self->pen->size() > 0) {
		gcValues.function   = GXxor;
		gcValues.line_width = 2;
		gcValues.line_style = LineSolid;
		XChangeGC(self->display, self->gc, GCFunction | GCLineWidth | GCLineStyle, &gcValues);
		XSetForeground(self->display, self->gc, 0xFFFFFF);

		for (int i=0; i < self->pen->size()-1; i++) {
			IPoint& a = self->pen->get(i);
			IPoint& b = self->pen->get(i+1);
			XDrawLine(self->display, self->pixmap, self->gc, a.x, a.y, b.x, b.y);
		}
	}

	// ---------- Display the self->pixmap ----------
	gcValues.function = GXcopy;
	XChangeGC(self->display, self->gc, GCFunction, &gcValues);
	XCopyArea(self->display, self->pixmap, self->window, self->gc,
		0, 0, width, height, 0, 0);
	//XFreePixmap(self->display, pixmap);

	Py_RETURN_NONE;
} // Viewer_expose

/* --- Viewer_export --- */
static PyObject* Viewer_export(ViewerObject *self, PyObject *args)
{
	char *filename;

	if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;
	self->viewer->exporter(filename);
	Py_RETURN_NONE;
} /* Viewer_export*/

/* _snapUV - snap the (u,v) position with d as snapping distance
 * @param u,v	position on the viewer plane
 * @param d	snapping distance
 * @return	Point(x,y,z) with absolute coordinates of the snapped position
 */
static Point _snapUV(ViewerObject *self, double u, double v, double d)
{
	double ds = d / self->kernel->view.Sx();

	// 1st Scan near-by objects
	ArrayIterator<GObject*> iter(*self->geometry->objectList);
	while (iter) {
		GObject *object = iter++;
		if (!object->show) continue;
		if (object->selnode>=0) continue;
		object->transform(self);
		if (!object->visible()) continue;

		// check the w?
		Point O = object->viewerPosition();

		if (Abs(u-O.x)<=ds && Abs(v-O.y)<=ds)
			// Use only the ou,ov coordinates of object!
			return self->kernel->view.uv2xyz(O.x, O.y);

		// Check END of GArrow and GRuler
		if (object->isA()==GArrowClass || object->isA()==GRulerClass) {
			double um = O.x;	// Remember for mid position checking
			double vm = O.y;

			// check the end position
			O = ((GArrow*)object)->viewerEnd();

			if (Abs(u-O.x)<=ds && Abs(v-O.y)<=ds)
				// Use only the ou,ov coordinates of object!
				return self->kernel->view.uv2xyz(O.x, O.y);

			um += O.x; um /= 2.0;
			vm += O.y; vm /= 2.0;
			if (Abs(u-um)<=ds && Abs(v-vm)<=ds)
				// Use only the ou,ov coordinates of object!
				return self->kernel->view.uv2xyz(um, vm);

			// check then angle handle
			if (object->isA()==GRulerClass && object->option==Ruler_Angle) {
				O = ((GRuler*)object)->viewerAngle();
				if (Abs(u-O.x)<=ds && Abs(v-O.y)<=ds)
					// Use only the ou,ov coordinates of object!
					return self->kernel->view.uv2xyz(O.x, O.y);
			}
		}
	}

	// 2nd scan viewports
	if (self->showViewport) {
		for (int ii=0; ii<NVIEWS; ii++) {
			ViewportLine& vp = self->viewport[ii];
			// FIXME (VISIBLE CENTER) Check center point
			double vx, vy, vz;
			double vu, vv;
			vp.viewer->origin(&vx, &vy, &vz);
			// Project on our viewport
			self->kernel->view.xyz2uv(vx, vy, vz, &vu, &vv);
			if (Abs(u-vu)<=ds && Abs(v-vv)<=ds)
				// Get x,y,z from projected center on our viewport!
				return self->kernel->view.uv2xyz(vu, vv);
		}
	}

	// 3rd scan for a close intersection if any
	double vu,vv;
	if (self->viewer->d2.showVertex && self->viewer->d2.closestVertex(u,v,ds, &vu,&vv))
		return self->kernel->view.uv2xyz(vu, vv);

	// 4th snap to 1/20th of the grid
	double dx = self->viewer->decoration.grid_dx / 20.0;
	double dy = self->viewer->decoration.grid_dy / 20.0;
	double dz = self->viewer->decoration.grid_dz / 20.0;
	double du = self->viewer->decoration.grid_du / 20.0;
	double dv = self->viewer->decoration.grid_dv / 20.0;

	// make it min least 2 pixels
	if (!Eq0(du,SMALL))
		while (!InRange(d, du*self->kernel->view.Sx(), d*10.0)) {
			double f = du*self->kernel->view.Sx()<d? 2.0 : 0.5;
			du *= f;
			switch (self->viewer->decoration.gridU) {
				case 'X':
					dx *= f;
					break;
				case 'Y':
					dy *= f;
					break;
				case 'Z':
					dz *= f;
					break;
			}
		}

	if (!Eq0(dv,SMALL))
		while (!InRange(d, dv*self->kernel->view.Sy(), d*10.0)) {
			double f = dv*self->kernel->view.Sy()<d? 2.0 : 0.5;
			dv *= f;
			switch (self->viewer->decoration.gridV) {
				case 'X':
					dx *= f;
					break;
				case 'Y':
					dy *= f;
					break;
				case 'Z':
					dz *= f;
					break;
			}
		}

	// Round to the closest of the grid
	double x = self->kernel->view.uv2x(u,v);
	double y = self->kernel->view.uv2y(u,v);
	double z = self->kernel->view.uv2z(u,v);

	bool calcXYZ = false;
	switch (self->viewer->decoration.gridU) {
		case 'X':
			x = (double)Round(x / dx) * dx;
			break;
		case 'Y':
			y = (double)Round(y / dy) * dy;
			break;
		case 'Z':
			z = (double)Round(z / dz) * dz;
			break;
		default:
			u = (double)Round(u / du) * du;
			calcXYZ = true;
	}
	switch (self->viewer->decoration.gridV) {
		case 'X':
			x = (double)Round(x / dx) * dx;
			break;
		case 'Y':
			y = (double)Round(y / dy) * dy;
			break;
		case 'Z':
			z = (double)Round(z / dz) * dz;
			break;
		default:
			v = (double)Round(v / dv) * dv;
			calcXYZ = true;
	}
	if (calcXYZ) {
		x = self->kernel->view.uv2x(u,v);
		y = self->kernel->view.uv2y(u,v);
		z = self->kernel->view.uv2z(u,v);
	}
	return Point(x,y,z);
} // _snapUV

/* --- Viewer_aspect --- */
static PyObject* Viewer_aspect(ViewerObject *self, PyObject *args)
{
	double aspect=-1.0;

	if (!PyArg_ParseTuple(args, "|d", &aspect)) return NULL;
	if (aspect<0.0)
		return PyFloat_FromDouble(self->kernel->view.aspect);
	else
		self->kernel->view.aspect = aspect;
	Py_RETURN_NONE;
} // Viewer_aspect

/* --- Viewer_basis ---
 * @return	return basis
 */
static PyObject* Viewer_basis(ViewerObject *self, PyObject *args)
{
	char axis;
	PyObject *basis=NULL;
	if (!PyArg_ParseTuple(args, "c|O", &axis, &basis)) return NULL;

	if (basis != NULL) {
		Py_RETURN_NONE;
	} else
		// return the normal matrix
		switch (axis) {
			case 'u':
			case 'U':
				return Py_Vector(self->kernel->view.axisu());
			case 'v':
			case 'V':
				return Py_Vector(self->kernel->view.axisv());
			case 'w':
			case 'W':
				return Py_Vector(self->kernel->view.axisw());
			default:
				Py_RETURN_NONE;
		}
} // Viewer_basis

/* --- Viewer_bbox ---
 * The same routine as the Geometry_bbox, but then transforms the bounding box
 * to the viewport coordinate system
 */
static PyObject* Viewer_bbox(ViewerObject *self, PyObject *args)
{
	char *type;
	PyObject *obj = NULL;
	BBox bbox;

	if (!PyArg_ParseTuple(args, "sO", &type, &obj)) return NULL;
	if (type[0]==0) {
		PyErr_SetString(PyExc_TypeError, "Invalid object type body, zone, region expected");
		return NULL;
	}

	if (type[0]=='b' || type[0]=='B') {
		GBody *body = Py_GBody(self->geometry,obj);
		if (body==NULL) return NULL;
		bbox = body->bbox();
	} else
	if (type[0]=='r' || type[0]=='R') {
		VRegion *region = Py_VRegion(self,obj);
		if (region==NULL) return NULL;
		bbox = region->region()->bbox();
	} else
//	if (type[0]=='z' || type[0]=='Z') {
//		// FIXME not implemented yet
//	} else
	if (type[0]=='o' || type[0]=='O') {
		GObject *object = Py_Object(self->geometry,obj);
		bbox = object->bbox();
	}

	if (bbox.isValid()) {
		bbox.transform(self->kernel->view.invMatrix());
		return Py_BuildValue("[dddddd]",
				bbox.low().x,  bbox.low().y,  bbox.low().z,
				bbox.high().x, bbox.high().y, bbox.high().z);
	} else
		Py_RETURN_NONE;
} // Viewer_bbox

/* --- Viewer_bbox2D --- */
static PyObject* Viewer_bbox2D(ViewerObject *self, PyObject *args)
{
	char *type;
	PyObject *obj = NULL;
	BBox bbox;

	if (!PyArg_ParseTuple(args, "sO", &type, &obj)) return NULL;
	if (type[0]==0) {
		PyErr_SetString(PyExc_TypeError, "Invalid object type body, zone, region expected");
		return NULL;
	}

	// FIXME: removed in refactor: there is no bbox2D in GBody
	if (type[0]=='b' || type[0]=='B') {
		VBody *body = Py_VBody(self,obj);
		if (body==NULL) return NULL;
		bbox = body->bbox2D();

	} else
	if (type[0]=='r' || type[0]=='R') {
		VRegion *region = Py_VRegion(self,obj);
		if (region==NULL) return NULL;
		bbox = region->region()->bbox2D();
	} else
	if (type[0]=='z' || type[0]=='Z') {
		// FIXME not implemented yet
	} else
	if (type[0]=='o' || type[0]=='O') {
		GObject *object = Py_Object(self->geometry,obj);
		bbox = object->bboxView(self);
	}

	if (bbox.isValid())
		return Py_BuildValue("[dddd]",
				bbox.low().x,  bbox.low().y,
				bbox.high().x, bbox.high().y);
	else
		Py_RETURN_NONE;
} // Viewer_bbox2D

/* --- Viewer_bodyVar --- */
static PyObject* Viewer_bodyVar(ViewerObject *self, GBody *body, const char *var, PyObject *value)
{
	if (!strcmp(var, "move")) {	// Move relative from the last save position
		if (value==NULL) {
			PyErr_SetString(PyExc_TypeError, "body move doesn't return anything.");
			return NULL;
		}
		// What to move
		int opt = PyInt_AsLong(value);
		self->geometry->geometry->lockWrite();
		body->move(opt, body->savedPosition()+self->move, self->kernel->view.axisw());
		body->create();
		if (body->hasMatrix()) body->transform();
		self->geometry->geometry->invalidateBody(body);
		self->geometry->geometry->unlockWrite();
	} else
	if (!strcmp(var, "rotate")) {	// rotate around an axis (axis,angle)
		// Calculate new position
		Point pos = body->savedPosition();
		// convert to relative
		double du = self->kernel->view.xyz2u(pos) - self->pivotU;
		double dv = self->kernel->view.xyz2v(pos) - self->pivotV;
		// Rotate around pivot point
		double u = self->rotateCos*du - self->rotateSin*dv + self->pivotU;
		double v = self->rotateSin*du + self->rotateCos*dv + self->pivotV;
		double w = self->kernel->view.xyz2w(pos);
		// back to absolute
		pos = self->kernel->view.uvw2xyz(u,v,w);

		self->geometry->geometry->lockWrite();

		// restore everything apart show
		int show = body->show;
		body->restore();
		body->show = show;

		// 1st move. RPP can change position if aligned to axis
		body->move(0, pos, self->kernel->view.axisw());
		// 2nd Rotate in place
		body->rotate(self->rotateAngle, self->rotateAxis);
		body->create();
		if (body->hasMatrix()) body->transform();
		self->geometry->geometry->invalidateBody(body);
		self->geometry->geometry->unlockWrite();
	} else {
		PyErr_Format(PyExc_TypeError, "Invalid type \'%s\' specified", var);
		return NULL;
	}
	Py_RETURN_NONE;
} // Viewer_bodyVar

/* --- Viewer_body --- */
static PyObject* Viewer_body(ViewerObject *self, PyObject *args)
{
	PyObject *obj;
	char *var;
	GBody *body = NULL;
	PyObject *value=NULL;
	if (!PyArg_ParseTuple(args, "Os|O", &obj, &var, &value)) return NULL;

	if (Py_Check4Pattern(obj)) {
	  //		char *pattern = PyString_AsString(obj);
	  const char *pattern = PyUnicode_AsUTF8AndSize(obj, lllen);
		ArrayIterator<GBody*> iter(self->geometry->geometry->bodies);
		while (iter) {
			body = iter++;
			if (!fnmatch(pattern, body->name(), 0)) {
				PyObject *ret = Viewer_bodyVar(self, body, var, value);
				if (ret==NULL) return NULL;
				Py_XDECREF(ret);	// ret can be NULL use XDECREF
			}
		}
		Py_RETURN_NONE;
	} else
	if (PyList_Check(obj)) {
		for (ssize_t i=0; i<PyList_GET_SIZE(obj); i++) {
			body = Py_GBody(self->geometry, PyList_GetItem(obj,i));
			if (body==NULL) return NULL;
			PyObject *ret = Viewer_bodyVar(self, body, var, value);
			if (ret==NULL) return NULL;
			Py_XDECREF(ret);	// ret can be NULL use XDECREF
		}
	} else
	if (PyTuple_Check(obj)) {
		for (ssize_t i=0; i<PyTuple_GET_SIZE(obj); i++) {
			body = Py_GBody(self->geometry, PyTuple_GetItem(obj,i));
			if (body==NULL) return NULL;
			PyObject *ret = Viewer_bodyVar(self, body, var, value);
			if (ret==NULL) return NULL;
			Py_XDECREF(ret);	// ret can be NULL use XDECREF
		}
	} else  {
		body = Py_GBody(self->geometry, obj);
		if (body==NULL) return NULL;
		return Viewer_bodyVar(self, body, var, value);
	}
	Py_RETURN_NONE;		// To keep compiler happy
} // Viewer_body

/* --- Viewer_camera --- */
static PyObject* Viewer_camera(ViewerObject *self, PyObject *args)
{
	char	*name;
	double	u=0.0, v=0.0;

	if (!PyArg_ParseTuple(args, "s|dd", &name, &u, &v)) return NULL;

	if (!strcmp(name,"position")) {
		double x,y,z;
		self->kernel->view.cameraPosition(&x,&y,&z);
		return Py_BuildValue("ddd", x, y, z);
	} else
	if (!strcmp(name,"direction")) {
		double dx,dy,dz;
		self->kernel->view.rayDirection(u,v, &dx,&dy,&dz);
		return Py_BuildValue("ddd", dx, dy, dz);
	} else
	if (!strcmp(name,"focal")) {
		if (u>0.0)
			self->kernel->view.focalLength(u);
		else
			return PyFloat_FromDouble(self->kernel->view.focalLength());
	} else
	if (!strcmp(name,"fov")) {
		if (u>0.0)
			self->kernel->view.fov(RAD(u));
		else
			return PyFloat_FromDouble(DEG(self->kernel->view.fov()));
	} else {
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid type option", name);
		return NULL;
	}
	Py_RETURN_NONE;
} // Viewer_camera

/* --- Viewer_calcViewport --- */
static void Viewer_calcViewport(ViewerObject *self, int id)
{
	ViewportLine& vp = self->viewport[id];
	const Matrix4& vpmatrix = vp.viewer->view().matrix();

	double w = (double)(self->viewer->width()  -VPMARGIN-4);
	double h = (double)(self->viewer->height() -VPMARGIN-2);

	// Find center of the viewport screen
	double vx, vy, vz;
	double vu, vv;
	vp.viewer->origin(&vx, &vy, &vz);
	self->kernel->view.xyz2uv(vx, vy, vz, &vu, &vv);

	// convert to screen coordinates
	vp.uc = self->kernel->view.u2id(vu);
	vp.vc = self->kernel->view.v2jd(vv);

	// w vector of other viewport
	double wx = vpmatrix(0,2);
	double wy = vpmatrix(1,2);
	double wz = vpmatrix(2,2);

	// The projected viewport LINE is NOT the intersection of
	// the two planes, rather a line parallel to this but
	// passing from the PROJECTED CENTER of the other viewport

	// project on current viewport (screen coordinates)
	// calculate the projection of the w-vector and get the
	// perpendicular (g,f) -> (f,-g) swapped
	const Matrix4& matrix = self->kernel->view.matrix();
	vp.g =   matrix(0,0)*wx + matrix(1,0)*wy + matrix(2,0)*wz;
	vp.f = -(matrix(0,1)*wx + matrix(1,1)*wy + matrix(2,1)*wz);
	// create perpendicular conic
	if (ISSMALL(vp.g) && ISSMALL(vp.f)) {
		vp.visible = false;
		return;
	}

	// The conic is projected on a viewport parallel to the current
	// one and passing from the center of the other viewport
	// line equation is g*x + f*y + c = 0
	vp.c = -(vp.uc*vp.g + vp.vc*vp.f);

	// normalize to be g^2+f^2=1.0
	double s = Sqr(vp.g) + Sqr(vp.f);
	if (!ISZERO(s)) {
		s = 1.0 / sqrt(s);
		vp.g *= s;
		vp.f *= s;
		vp.c *= s;
	}

	// find intersections with viewport visible window
	int n = 0;	// intersection counter
	if (!ISZERO(vp.f)) {
		// x = 0
		double y = -(vp.c + vp.g*VPMARGIN) / vp.f;
		if (InRange((double)VPMARGIN, y, h)) {
			vp.x[n] = VPMARGIN;
			vp.y[n] = Round(y);
			n++;
		}
		// x = w
		y = -(vp.c + vp.g*w) / vp.f;
		if (InRange((double)VPMARGIN, y, h)) {
			vp.x[n] = (int)w;
			vp.y[n] = Round(y);
			n++;
		}
	}
	if (!ISZERO(vp.g)) {
		// y = 0
		double x = -(vp.c + vp.f*VPMARGIN) / vp.g;
		if (n<2 && InRange((double)VPMARGIN, x, w)) {
			vp.x[n] = Round(x);
			vp.y[n] = VPMARGIN;
			n++;

		}
		// y = h
		x = -(vp.c + vp.f*h) / vp.g;
		if (n<2 && InRange((double)VPMARGIN, x, w)) {
			vp.x[n] = Round(x);
			vp.y[n] = (int)h;
			n++;
		}
	}

	// after all calculation move the center in the visible area
	if (vp.uc<VPMARGIN || vp.uc>=w || vp.vc<VPMARGIN || vp.vc>=h) {
		if (ISZERO(vp.g) || ISZERO(vp.f)) {
			// vertical or horizontal line
			vp.xc = Range(VPMARGIN, Round(vp.uc), (int)w);
			vp.yc = Range(VPMARGIN, Round(vp.vc), (int)h);
		} else
		if (n==0) {
			if (vp.g * vp.f > 0.0) {	// same sign
				// can only be the 0,0 or the w,h
				if (Abs(vp.c) < Abs(vp.g*w + vp.f*h + vp.c)) {
					vp.xc = VPMARGIN;
					vp.yc = VPMARGIN;
				} else {
					vp.xc = (int)w;
					vp.yc = (int)h;
				}
			} else {
				// can only be the w,0 or the 0,h
				if (Abs(vp.g*w+vp.c) < Abs(vp.f*h + vp.c)) {
					vp.xc = (int)w;
					vp.yc = VPMARGIN;
				} else {
					vp.xc = VPMARGIN;
					vp.yc = (int)h;
				}
			}
		} else
		if (n==1) {
			vp.xc = Round(vp.x[0]);
			vp.yc = Round(vp.y[0]);
		} else {
			// diagonal line
			double d1 = Sqr(vp.x[0]-vp.uc) + Sqr(vp.y[0]-vp.vc);
			double d2 = Sqr(vp.x[1]-vp.uc) + Sqr(vp.y[1]-vp.vc);
			if (d1<d2) {
				vp.xc = Round(vp.x[0]);
				vp.yc = Round(vp.y[0]);
			} else {
				vp.xc = Round(vp.x[1]);
				vp.yc = Round(vp.y[1]);
			}
		}
	} else {
		vp.xc = Round(vp.uc);
		vp.yc = Round(vp.vc);
	}

	// End position of the w vector
	if (vp.viewer->view().projection==Projection_Perspective) {
		double x,y,z;
		vp.viewer->view().cameraPosition(&x,&y,&z);
		double u = self->kernel->view.xyz2u(x,y,z);
		double v = self->kernel->view.xyz2v(x,y,z);
		if (self->kernel->view.inside(u,v)) {
			vp.xw = self->kernel->view.u2i(u);
			vp.yw = self->kernel->view.v2j(v);
		} else {
			vp.xw = vp.xc + Round(vp.g*vp.wLength);
			vp.yw = vp.yc + Round(vp.f*vp.wLength);
		}
	} else {
		vp.xw = vp.xc + Round(vp.g*vp.wLength);
		vp.yw = vp.yc + Round(vp.f*vp.wLength);
	}

	// corrected c for visual center
	vp.cc = -((double)vp.xc*vp.g + (double)vp.yc*vp.f);

	// ViewportLine is outside
	if (n==0) {
		// find closest line to draw
		if (ISZERO(vp.g)) {
			// horizontal
			vp.x[0] = VPMARGIN;
			vp.x[1] = (int)w;
			vp.y[0] = vp.y[1] = vp.yc;
		} else
		if (ISZERO(vp.f)) {
			// vertical
			vp.x[0] = vp.x[1] = vp.xc;
			vp.y[0] = VPMARGIN;
			vp.y[1] = (int)h;
		} else {
			vp.x[0] = vp.x[1] = vp.xc;
			vp.y[0] = vp.y[1] = vp.yc;
		}
	}

	vp.visible = true;
} // Viewer_calcViewport

/* --- Viewer_closest --- */
static PyObject* Viewer_closest(ViewerObject *self, PyObject *args)
{
	int i,j,k;
	double pi, pj, d=1.0;
	PyObject *opt=NULL, *obj=NULL;
	char *options = "vCVOBZR";

	if (!PyArg_ParseTuple(args, "ii|dOO", &i, &j, &d, &opt, &obj)) return NULL;
	pi = (double)i;
	pj = (double)j;
	double u = self->kernel->view.i2u((int)pi);
	double v = self->kernel->view.j2v((int)pj);

	//	if (opt && PyString_Check(opt)) options = PyString_AsString(opt);
	if (opt && PyUnicode_Check(opt)) options = (char* ) PyUnicode_AsUTF8AndSize(opt, lllen);
	// Loop over the items to search
	for (char *o=options; *o; o++) {
		switch (*o) {
			// ------ Search Viewports (ORIGIN ONLY) ------
			case 'v':
				if (self->showViewport) {
					int min = -1;
					double mind = INFINITE;
					double d1=d+(double)self->viewport[0].originWidth;
					// check the cursor
					for (k=0; k<NVIEWS; k++) {
						ViewportLine& vp = self->viewport[k];
						if (vp.visible) {
							double c = Sqr(pi-(double)vp.xc) + Sqr(pj-(double)vp.yc);
							if (c <= mind) {
								mind = c;
								min  = k;
							}
						}
					}
					if (min>=0 && mind<=Sqr(d1))
						return Py_BuildValue("sis","V",min,NULL);
				}
				break;

			// ------ Search Viewports (2nd) (Line) ------
			case 'V':
				if (self->showViewport) {
					int min = -1;
					double mind = INFINITE;
					double d1=d+(double)self->viewport[0].originWidth;
					// check the line
					for (k=0; k<NVIEWS; k++) {
						ViewportLine& vp = self->viewport[k];
						if (vp.visible) {
							double c = Abs(vp.g*pi + vp.f*pj + vp.cc);
							if (c <= mind) {
								mind = c;
								min  = k;
							}
						}
					}
					if (min>=0 && mind<=d1)
						return Py_BuildValue("sis","V",min,NULL);
				}
				break;

			// ------ Search Objects ------
			case 'O': {
				ArrayIterator<GObject*> iter(*self->geometry->objectList);
				while (iter) {
					GObject *object = iter++;
					if (!object->show) continue;
					object->transform(self);
					if (!object->visible()) continue;
					int close = object->closest(self,i,j,Round(d));
					if (close>=0)
						return Py_BuildValue("sii","O",object->id(),close);
				}
				} break;

			// ------ Search Editing body ------
			case 'C':
				if (self->geometry->cursorShow) {
					int ci, cj;
					self->kernel->view.xyz2ij3D(self->geometry->cursor, &ci, &cj);
					double ici = (double)(i-ci);
					double jcj = (double)(j-cj);
					double dist = sqrt(Sqr(ici) + Sqr(jcj));

					// Free move
					if (dist <= (double)self->geometry->cursorMoveSize + d)
						return Py_BuildValue("sii","C", 0, 0);

					double r = (double)self->geometry->cursorSize;
					int   d2 = Round(Sqr(d));

					// Axis X
					const Matrix4& matrix = self->kernel->view.matrix();
					double Dx =  r * matrix(0,0);
					double Dy = -r * matrix(0,1);
					double Dlen2 = Sqr(Dx) + Sqr(Dy);
					//if (Dlen2 < 2) return -1;
					// Check distance, calculate the dot product (to the square)
					// of the two vs length
					// Calculate the projection distance using the dot product
					if ( Dlen2>3.0 &&
					    (Sqr(Dy*ici - Dx*jcj) < d2*Dlen2) &&
					    (InRange(0.0, Dx*ici + Dy*jcj, Dlen2)))
						return Py_BuildValue("sii","C", 0, -1);

					// Axis Y
					Dx =  r * matrix(1,0);
					Dy = -r * matrix(1,1);
					Dlen2 = Sqr(Dx) + Sqr(Dy);
					if ( Dlen2>3.0 &&
					    (Sqr(Dy*ici - Dx*jcj) < d2*Dlen2) &&
					    (InRange(0.0, Dx*ici + Dy*jcj, Dlen2)))
						return Py_BuildValue("sii","C", 0, -2);

					// Axis Z
					Dx =  r * matrix(2,0);
					Dy = -r * matrix(2,1);
					Dlen2 = Sqr(Dx) + Sqr(Dy);
					if ( Dlen2>3.0 &&
					    (Sqr(Dy*ici - Dx*jcj) <= d2 * Dlen2) &&
					    (InRange(0.0, Dx*ici + Dy*jcj, Dlen2)))
						return Py_BuildValue("sii","C", 0, -3);

					// Rotate around w
					if (Abs(dist-r) <= d)
						return Py_BuildValue("sii","C", 0, -4);
				}
				break;

			// ------ Search Bodies ------
			case 'B': {
				int dd = Round(d);
				// FIXME: removed in refactor: cannot access kernel->getBody in all
				// cases, the body may not exist
				bool ebvisi = false;	// edit body visibility
				if (self->geometry->geometry->editBody())	// remember visibility
					ebvisi = self->kernel->getBody(self->geometry->geometry->editBody())->visible;

				// Draw in the clipping region to mark the visible ones
				self->viewer->painter.clip(Max(0,i-dd),Max(0,j-dd), i+dd,j+dd);
				self->kernel->clearVisibleBodies();

				self->viewer->draw(DRAW_SEGMENTS|DRAW_WIREFRAME);

				self->viewer->painter.resetClip();

				// First check editBody if any
				if (self->geometry->geometry->editBody() &&
				   (self->geometry->geometry->editBody()->show & BIT_SELECT) &&
				    self->kernel->getBody(self->geometry->geometry->editBody())->visible) {
					Point r = self->kernel->view.uv2xyz(u,v);
					int close = self->geometry->geometry->editBody()->closest(r,
									d/self->kernel->view.Sx(),
									self->kernel->view.axisw());

					// Special case for PLA and position is visible
					assert(YZPbody<XZPbody && XZPbody<XYPbody && XYPbody< PLAbody);
					if (InRange(YZPbody, self->geometry->geometry->editBody()->type(), PLAbody)) {
						// Check if we are close to the borders like in the viewports
						int ww = self->viewer->width()  / 20;
						int hh = self->viewer->height() / 20;
						Point p = self->geometry->geometry->editBody()->position();
						int ei = self->kernel->view.u2i(self->kernel->view.xyz2u(p));
						int ej = self->kernel->view.v2j(self->kernel->view.xyz2v(p));
						if ((i<=ww && ei>ww) ||
						    (i>=self->viewer->width()-ww && ei<self->viewer->width()-ww) ||
						    (j<=hh && ej>hh) ||
						    (j>=self->viewer->height()-hh && ej<self->viewer->height()-hh)) {
							if (close == 1) {
								// Set also the sign of rotation
								((GPLABody*)self->geometry->geometry->editBody())->signMove20(r,
										self->kernel->view.axisw());
								return Py_BuildValue("sii","B",self->geometry->geometry->editBody()->id(),20);
							}
						}
					}
					return Py_BuildValue("sii", "B", self->geometry->geometry->editBody()->id(), close);
				}

				// scan visible bodies to find the closest one
				// FIXME: potential source of troubles
				self->kernel->lock();
				double min = INFINITE;
				VBody *closestBody = NULL;
				for (int ib = 0; ib < self->kernel->nGeometryBodies(); ib++) {
					VBody *body = self->kernel->bodies[ib];
					if (!body->visible) continue;
					double x = self->kernel->view.uv2x(u,v);
					double y = self->kernel->view.uv2y(u,v);
					double z = self->kernel->view.uv2z(u,v);
					double m = body->distance(x,y,z);
					if (m<min) {
						min = m;
						closestBody = body;
					}
				}
				self->kernel->unlock();

				// FIXME: removed in refactor: cannot access kernel->getBody in all
				// cases, the body may not exist
				// restore visibility of editBody
				if (self->geometry->geometry->editBody())
					self->kernel->getBody(self->geometry->geometry->editBody())->visible = ebvisi;

				if (closestBody) {
					int close = 0;
					if (closestBody->body()==self->geometry->geometry->editBody() &&
					   (self->geometry->geometry->editBody()->show & BIT_SELECT)) {
						// Find closest node/edge on the body
						Point r = self->kernel->view.uv2xyz(u,v);
						close = closestBody->body()->closest(r,
									d/self->kernel->view.Sx(),
									self->kernel->view.axisw());
					}
					return Py_BuildValue("sii", "B", closestBody->id(), close);
				}
				} break;

			// ------ Search Zones ------
			case 'Z':
				if (obj!=NULL) {
					VRegion *region = Py_VRegion(self, obj);

					// Check in 2D mode
					if (self->viewer->d2.fillRegions) {
						if (region != NULL) {
							// Iterate over the zones composing the region
							// find find if we are inside
							int zi=0;
							// FIXME: potential source of troubles
							self->kernel->lock();
							self->kernel->engine()->incBodyCheckId();
							PyObject *zones = NULL;
							const ViewPort& V = self->kernel->view;
							const Matrix4& M = V.matrix();
							for (int iz=0; iz<region->zones().size(); iz++) {
								const VZone *zone = region->zones()[iz];
								double x = V.uv2x(u,v);
								double y = V.uv2y(u,v);
								double z = V.uv2z(u,v);
								if (zone->inside2D(self->kernel->engine(),
								    x,y,z,-M(0,2),-M(1,2),-M(2,2))) {
									if (zones==NULL) zones = PyList_New(0);
									PyList_Append(zones, PyInt_FromLong(zi));
								}
								zi++;
							}
							self->kernel->unlock();
							if (zones) return Py_BuildValue("sOs","Z",zones,NULL);
						} else
							// avoid the message, especially when voxel
							// region is involved
							PyErr_Clear();
					}

					// if failed check in 3D mode
					if (self->viewer->d3.show) {
						self->kernel->lock();
						self->kernel->engine()->incBodyCheckId();
						const ViewPort& V = self->kernel->view;
//						const Matrix4& M = V.matrix();
						double x = V.uv2x(u,v);
						double y = V.uv2y(u,v);
						double z = V.uv2z(u,v);
//						VZone *zone= region->inside2D(self->kernel->engine(), x,y,z, -M(0,2),-M(1,2),-M(2,2));
//						if (zone->region()->transparent()) {
//
						Ray ray;
						ray.use_clip      = self->kernel->clipBodyCount() > 0;
						ray.skip_1stblack = self->viewer->d3.skip_1stblack;
						self->kernel->view.rayPosition(u,v, &x,&y,&z);
						double dx,dy,dz;
						self->kernel->view.rayDirection(u,v, &dx,&dy,&dz);
						// Find the new zone of the camera
						self->kernel->engine()->incBodyCheckId();
						VZone* zone = self->kernel->engine()->whereRay(x,y,z, dx,dy,dz, SMALL3D, NULL);
						if (zone) {
							if (zone->region() != region && zone->region()->transparent()) {
								ray.push(RaySegment(x,y,z, dx,dy,dz, zone));
								self->kernel->engine()->incBodyCheckId();
								self->viewer->d3.nextIntersection(self->kernel->engine(), &ray);
								zone = ray.hitZone();
							}
						} else
							zone = NULL;
						self->kernel->unlock();
						if (zone && zone->region() == region) {
							PyObject* zones = PyList_New(0);
							PyList_Append(zones, PyInt_FromLong(zone->id()));
							return Py_BuildValue("sOs","Z",zones,NULL);
						}
					}
				}
				break;

			// ------ Search Regions ------
			case 'R':
				if ((self->viewer->d2.fillRegions||self->viewer->d3.show)) {
					// Find region under point
					// FIXME: potential source of troubles
					self->kernel->lock();
					self->kernel->engine()->incBodyCheckId();
					PyObject *regions = NULL;
					const ViewPort& V = self->kernel->view;
					const Matrix4& M = V.matrix();
					for (int ir=0; ir<self->kernel->nGeometryRegions(); ir++) {
						VRegion *region = self->kernel->regions[ir];
						double x = V.uv2x(u,v);
						double y = V.uv2y(u,v);
						double z = V.uv2z(u,v);
						VZone *zone= region->inside2D(self->kernel->engine(), x,y,z, -M(0,2),-M(1,2),-M(2,2));
						// Check for 3D if needed
						bool d3used;
						if (self->viewer->d3.show && (zone==NULL || zone->region()->transparent())) {
							//self->kernel->engine()->incBodyCheckId();
							//zone = self->kernel->engine()->whereRay(x,y,z, dx,dy,dz, SMALL3D, zone);
							Ray ray;
							ray.use_clip      = self->kernel->clipBodyCount() > 0;
							ray.skip_1stblack = self->viewer->d3.skip_1stblack;

							self->kernel->view.rayPosition(u,v, &x,&y,&z);

							double dx,dy,dz;
							self->kernel->view.rayDirection(u,v, &dx,&dy,&dz);

							// Find the new zone of the camera
							self->kernel->engine()->incBodyCheckId();
							zone = self->kernel->engine()->whereRay(x,y,z, dx,dy,dz, SMALL3D, NULL);

							ray.push(RaySegment(x,y,z, dx,dy,dz, zone));
							self->kernel->engine()->incBodyCheckId();
							self->viewer->d3.nextIntersection(self->kernel->engine(), &ray);
							zone = ray.hitZone();
							d3used = true;
						} else
							d3used = false;
						if (zone) {
							if (regions==NULL) regions = PyList_New(0);
							PyList_Append(regions, PyInt_FromLong(zone->region()->id()));
							if (d3used) break;
						}
					}
					self->kernel->unlock();
					if (regions) return Py_BuildValue("sOs","R",regions,NULL);
				}
				break;
			default:
				assert(0);
		}
	}
	Py_RETURN_NONE;
} // Viewer_closest

/* --- Viewer_derive --- */
static PyObject* Viewer_derive(ViewerObject *self)
{
	self->kernel->derive();
	Py_RETURN_NONE;
} // Viewer_derive

/* _endDraw: notifier function for end of projection */
static void _endDraw(ViewerObject *self)
{
	DUMP(cout << "_endDraw("<<self->viewer->title()<<")" << endl);
	_sendMessage(self, "EndDraw");
} // _endDraw

/* --- Viewer_draw --- */
static PyObject* Viewer_draw(ViewerObject *self, PyObject *args)
{
	int asThread = false;
	int mask = -1;
	if (!PyArg_ParseTuple(args, "|ii", &asThread, &mask)) return NULL;
	DUMP(cout << "Viewer_draw("<<self->viewer->title()<<", " << asThread << ")" << endl);
	if (asThread) {
		self->viewer->spawnDraw((NotifyFunc)_endDraw, self);
		Py_RETURN_NONE;
	} else {
		// Stop any pending drawing thread
		if (self->viewer->state() == PROJECTION_DRAW) {
			DUMP(cout << "Viewer_draw(" << self->viewer->title() << ").stopThread()" << endl);
			self->viewer->stopThread();
		}
		// Run synchronously
		return PyInt_FromLong(self->viewer->draw(mask));
	}
} // Viewer_draw

/* --- Viewer_calcWindow --- */
static PyObject* Viewer_calcWindow(ViewerObject *self, PyObject *args)
{
	double zoom;

	if (!PyArg_ParseTuple(args, "d", &zoom)) return NULL;

	self->kernel->view.calcWindow(zoom);
	Py_RETURN_NONE;
} // Viewer_calcWindow

/* --- Viewer_edit --- */
static PyObject* Viewer_edit(ViewerObject *self, PyObject *args)
{
	char *name;
	PyObject *value=NULL;
	PyObject *value2=NULL;

	if (!PyArg_ParseTuple(args, "s|OO", &name, &value, &value2)) return NULL;

	if (!strcmp(name, "move")) {	// Calculate move distance
		if (value) {
			// Expect: (uv, axis, snap, relative)
			if (!PyTuple_Check(value) || PyTuple_GET_SIZE(value)!=4) {
				PyErr_SetString(PyExc_TypeError, "tuple of size 4 expected");
				return NULL;
			}

			// What to move
			double u,v;
			if (!Py_GetUV(PyTuple_GetItem(value, 0), &u, &v)) return NULL;
			int axis = PyInt_AsLong(PyTuple_GetItem(value, 1));
			int snap = PyInt_AsLong(PyTuple_GetItem(value, 2));
			bool relative = (bool)PyInt_AsLong(PyTuple_GetItem(value, 3));
			if (PyErr_Occurred()) return NULL;

			// Real or aligned coordinates
			Point D = snap?
					_snapUV(self, u, v, self->geometry->snapDistance) :
					self->kernel->view.uv2xyz(u,v);
			self->move = self->start;

			// Find difference of self->move from absolute coordinates of projection on viewport
			Vector dS = self->move - self->kernel->view.uv2xyz(
						 self->kernel->view.xyz2u(self->move),
						 self->kernel->view.xyz2v(self->move));

			// Axis=value[2] Find new position
			switch (axis) {
				case 0:	// X-axis
					self->move.x = D.x;
					if (relative) self->move.x += dS.x;
					break;
				case 1: // Y-axis
					self->move.y = D.y;
					if (relative) self->move.y += dS.y;
					break;
				case 2: // Z-axis
					self->move.z = D.z;
					if (relative) self->move.z += dS.z;
					break;
				default:
					self->move = D;
					if (relative) self->move += dS;
			}
			self->move -= self->start;
		} else
			return Py_Vector(self->move);
	} else
	if (!strcmp(name, "rotate")) {	// Calculate rotate angle
		if (value) {
			// Expect: (uv,  axis, snap))
			if (!PyTuple_Check(value) || PyTuple_GET_SIZE(value)!=3) {
				PyErr_SetString(PyExc_TypeError, "tuple of size 3 expected");
				return NULL;
			}
			// What to move
			double u,v;
			if (!Py_GetUV(PyTuple_GetItem(value, 0), &u, &v)) return NULL;
			self->rotateAxis = Py_GetVector(PyTuple_GetItem(value, 1));
			int snap = PyInt_AsLong(PyTuple_GetItem(value, 2));
			if (PyErr_Occurred()) return NULL;

			// Calculate angle
			self->rotateAngle = atan2(v-self->pivotV, u-self->pivotU) - self->pivotAngle;
			if (snap) // FIXME angle should be past from python
				self->rotateAngle = (double)Round(
						self->rotateAngle / self->geometry->snapAngle)
					* self->geometry->snapAngle;
			bsincos(self->rotateAngle, &self->rotateSin, &self->rotateCos);
		} else
			return PyFloat_FromDouble(DEG(self->rotateAngle));
	} else
	if (!strcmp(name, "value")) {	// Calculate move by typing
		if (value) {
			// Move list(dx[,dy[,dz]]), int(option)
			if (PyList_Check(value)) {
				double a = PyFloat_AsDouble(PyList_GetItem(value, 0));
				double b=0.0, c=0.0;
				if (PyList_Size(value)>1) {
					b = PyFloat_AsDouble(PyList_GetItem(value, 1));
					if (PyList_Size(value)>2)
						c = PyFloat_AsDouble(PyList_GetItem(value, 2));
				}

				self->move.x = 0.0;
				self->move.y = 0.0;
				self->move.z = 0.0;
				if (value2==NULL) {
					self->move.x = a;
					self->move.y = b;
					self->move.z = c;
				} else
					switch (PyInt_AsLong(value2)) {
						case 0:	// X-axis
							self->move.x = a;
							break;
						case 1: // Y-axis
							self->move.y = a;
							break;
						case 2: // Z-axis
							self->move.z = a;
							break;
						default:
							self->move.x = a;
							self->move.y = b;
							self->move.z = c;
					}
			} else {
				// Angle   float, list(axis)
				self->rotateAngle = RAD(PyFloat_AsDouble(value));
				self->rotateAxis  = Py_GetVector(value2);
			}
		} else
			return Py_Vector(self->move);
	} else
	if (!strcmp(name, "start")) {	// Set start point
		if (value) {
			// Expect ( (u,v)|(x,y,z), snap )
			if (!PyTuple_Check(value) || PyTuple_GET_SIZE(value)!=2) {
				PyErr_SetString(PyExc_TypeError, "tuple of size 2 (uv|xyz, snap) expected");
				return NULL;
			}
			PyObject *pos = PyTuple_GetItem(value,0);
			int snap = PyInt_AsLong(PyTuple_GetItem(value, 1));

			// Relative if (u,v)
			if (PyTuple_Check(pos) && PyTuple_GET_SIZE(pos)==2) {
				self->startU = PyFloat_AsDouble(PyTuple_GetItem(pos,0));
				self->startV = PyFloat_AsDouble(PyTuple_GetItem(pos,1));
				if (snap) {
					Point D = _snapUV(self, self->startU, self->startV,
							self->geometry->snapDistance);
					self->startU = self->kernel->view.xyz2u(D);
					self->startV = self->kernel->view.xyz2v(D);
				}
				if (PyErr_Occurred()) return NULL;
				self->start = self->kernel->view.uv2xyz(self->startU, self->startV);
			} else {	// Otherwise absolute x,y,z
				self->start = Py_GetVector(pos);
				if (PyErr_Occurred()) return NULL;
				self->startU = self->kernel->view.xyz2u(self->start);
				self->startV = self->kernel->view.xyz2v(self->start);
			}
		} else
			return Py_BuildValue("dd",self->startU,self->startV);
	} else
	if (!strcmp(name, "pivot")) {	// Set pivot point
		if (value) {
			if (PyTuple_Check(value) && PyTuple_GET_SIZE(value)==2) {
				self->pivotU = PyFloat_AsDouble(PyTuple_GetItem(value,0));
				self->pivotV = PyFloat_AsDouble(PyTuple_GetItem(value,1));
				self->pivotAngle = atan2(self->startV-self->pivotV, self->startU-self->pivotU);
			}
		} else
			return Py_BuildValue("dd",self->pivotU,self->pivotV);
	} else
	if (!strcmp(name, "cursor")) {	// Set start position from cursor
			self->start  = self->geometry->cursor;
			self->startU = self->kernel->view.xyz2u(self->start);
			self->startV = self->kernel->view.xyz2v(self->start);
	} else
	if (!strcmp(name, "object")) {	// Set start position from object
		if (value) {
			GObject *object = Py_Object(self->geometry,value);
			int opt = PyInt_AsLong(value2);	// What to move
			if (PyErr_Occurred()) return NULL;
			self->start  = object->node(opt);
			self->startU = self->kernel->view.xyz2u(self->start);
			self->startV = self->kernel->view.xyz2v(self->start);
		} else
			return Py_Vector(self->start);
	} else
	if (!strcmp(name, "body")) {	// Set start position from body,item pair
		if (value) {
			GBody* body  = Py_GBody(self->geometry, value);
			self->start  = body->position();
			self->startU = self->kernel->view.xyz2u(self->start);
			self->startV = self->kernel->view.xyz2v(self->start);
		} else
			return Py_Vector(self->start);
	} else {
		PyErr_Format(PyExc_TypeError, "Invalid type \'%s\' specified", name);
		return NULL;
	}
	Py_RETURN_NONE;
} // Viewer_edit

/* --- Viewer_enclosed --- */
static PyObject* Viewer_enclosed(ViewerObject *self, PyObject *args)
{
	int l,t,r,b;
	char *options=NULL;

	if (!PyArg_ParseTuple(args, "iiii|s", &l, &t, &r, &b, &options)) return NULL;

	if (r<l) Swap(l,r);
	if (b<t) Swap(t,b);

//	double u = self->kernel->view.realX((int)pi);
//	double v = self->kernel->view.realY((int)pj);
//	double x = self->kernel->view.uv2x(u,v);
//	double y = self->kernel->view.uv2y(u,v);
//	double z = self->kernel->view.uv2z(u,v);

	self->viewer->painter.clip(l, t, r, b);
	self->kernel->clearVisibleBodies();
#ifdef EXPERIMENTAL
#if 0
	self->viewer->clearVisibleRegions();
#endif
#endif
	self->viewer->draw(DRAW_SEGMENTS|DRAW_WIREFRAME);
	self->viewer->painter.resetClip();

	if (options==NULL || options[0]=='B') {
		PyObject *bodies=NULL;
		// scan bodies to find the visible ones
		// FIXME: potential source of troubles
		self->kernel->lock();
		bodies = PyList_New(0);
		for (int ib = 0; ib < self->kernel->nGeometryBodies(); ib++) {
			VBody *body = self->kernel->bodies[ib];
			if (!body->visible) continue;
			PyList_Append(bodies, PyInt_FromLong(body->id()));
		}
		self->kernel->unlock();

		return bodies;
	}

	if (options==NULL || strchr(options,'O')) {
		PyObject *objects=NULL;
		// scan regions to find the visible ones
		ArrayIterator<GObject*> iter(*self->geometry->objectList);
		objects = PyList_New(0);
		while (iter) {
			GObject *object = iter++;
			if (!object->show) continue;
			object->transform(self);
			if (!object->visible()) continue;
			if (object->enclosed(self, l,t,r,b))
				PyList_Append(objects, PyInt_FromLong(object->id()));
		}
		return objects;
	}

	Py_RETURN_NONE;
} // Viewer_enclosed

/* --- Viewer_error --- */
static PyObject* Viewer_error(ViewerObject *self, PyObject *args)
{
	char *cmd;
	int err = 0;
	if (!PyArg_ParseTuple(args, "s|i", &cmd, &err)) return NULL;

	if (!strcmp(cmd, "n"))		// return number of errors
		return PyInt_FromLong(self->kernel->errors());
	else
	if (!strcmp(cmd, "clear")) {	// clear error points
		_clearErrors(self);
		self->kernel->clearErrors();
	} else
	if (!strcmp(cmd, "show")) {
		_clearErrors(self);
		// create temporary point objects
		if (err==0)
			err = self->kernel->errors();
		else
			err = Min(err, self->kernel->errors());
		for (int i=1; i<=err; i++) {
			char name[10];
			sprintf(name,"%d",i);
			GObject *obj = new GPoint(name, Point_X);
			obj->color = self->geometry->geometry->errorColor;
			obj->anchor = Anchor_E;

			ZoneOfPoint pIn, pOut;
			double u,v, xmin, xmax, ymin, ymax;
			self->kernel->error(i, &u, &v, &pIn, &pOut, &xmin, &xmax, &ymin, &ymax);
			obj->position(self->kernel->view.uv2xyz(u,v));
			self->errors->add(obj);
		}
	} else
	if (!strcmp(cmd, "get")) {
		ZoneOfPoint pIn, pOut;
		double u,v;
		double xmin, xmax, ymin, ymax;
		VBody *body = self->kernel->error(err, &u, &v, &pIn, &pOut, &xmin, &xmax, &ymin, &ymax);
		if (body==NULL) Py_RETURN_NONE;

		PyObject *inList  = PyList_New(pIn.n);
		for (int i=0; i<pIn.n; i++)
			PyList_SetItem(inList, i,
				PyString_FromFormat("%s:%d",
						pIn.zone[i]->region()->name(),
						pIn.zone[i]->id()+1));
				//PyString_FromString(pIn.zone[i]->region()->name()));

		PyObject *outList = PyList_New(pOut.n);
		for (int i=0; i<pOut.n; i++)
			PyList_SetItem(outList, i,
				PyString_FromFormat("%s:%d",
						pOut.zone[i]->region()->name(),
						pOut.zone[i]->id()+1));
				//PyString_FromString(pOut.zone[i]->region()->name()));

		return Py_BuildValue("dddddddsOO",
			self->kernel->view.uv2x(u,v),
			self->kernel->view.uv2y(u,v),
			self->kernel->view.uv2z(u,v),
			xmin, xmax, ymin, ymax,
			body->name(),
			inList, outList);
	} else
	if (!strcmp(cmd, "hash")) {
		return PyInt_FromLong(self->kernel->errorHash()*self->viewer->state());
	} else
	if (!strcmp(cmd, "kernel")) {
		if (self->kernel->error()[0])
			return PyString_FromString(self->kernel->error());
	} else
	if (!strcmp(cmd, "reset")) {
		self->kernel->errorReset();
	} else {
		PyErr_Format(PyExc_TypeError, "Invalid type \'%s\' specified", cmd);
		return NULL;
	}
	Py_RETURN_NONE;
} // Viewer_error

/* --- Viewer_extends --- */
static PyObject* Viewer_extends(ViewerObject *self, PyObject *args)
{
	double xmin, xmax, ymin, ymax;
	if (PyTuple_Size(args)!=4)
		return Py_BuildValue("dddd",
				self->kernel->view.minu(),
				self->kernel->view.minv(),
				self->kernel->view.maxu(),
				self->kernel->view.maxv());
	else {
		if (!PyArg_ParseTuple(args, "dddd", &xmin, &ymin, &xmax, &ymax)) return NULL;
		self->viewer->window(xmin, ymin, xmax, ymax);
	}
	Py_RETURN_NONE;
} // Viewer_extends

/* --- Viewer_fpe --- */
static PyObject* Viewer_fpe(ViewerObject *)
{
	fpetrap();
	Py_RETURN_NONE;
} // Viewer_fpe

/* --- Viewer_font --- */
static PyObject *Viewer_font(ViewerObject *self, PyObject *args)
{
	char *fontname = NULL;
	char *fontpath = NULL;
	if (!PyArg_ParseTuple(args, "s|s", &fontname, &fontpath)) return NULL;

	if (!strcmp(fontname,"general")) {
		if (fontpath==NULL) {
			if (self->viewer->font.name().c_str())
				return PyString_FromString(self->viewer->font.name().c_str());
		} else
			loadFont(self->viewer->font, fontpath);
	} else
	if (!strcmp(fontname,"grid")) {
		if (fontpath==NULL) {
			if (self->viewer->decoration.gridFont.name().c_str())
				return PyString_FromString(self->viewer->decoration.gridFont.name().c_str());
		} else
			loadFont(self->viewer->decoration.gridFont, fontpath);
	} else
	if (!strcmp(fontname,"palette")) {
		if (fontpath==NULL) {
			if (self->viewer->palette.font.name().c_str())
				return PyString_FromString(self->viewer->palette.font.name().c_str());
		} else
			loadFont(self->viewer->palette.font, fontpath);
	} else {
		PyErr_Format(PyExc_TypeError, "Invalid font \'%s\' specified", fontname);
		return NULL;
	}

	Py_RETURN_NONE;
} // Viewer_font

/* --- Viewer_grid --- */
static PyObject* Viewer_grid(ViewerObject *self, PyObject *args)
{
	char *name;
	PyObject *value=NULL;
	if (!PyArg_ParseTuple(args, "s|O", &name, &value)) return NULL;

	if (!strcmp(name,"label")) {
		if (value==NULL) {
			char ul[2], vl[2];
			ul[0] = self->viewer->decoration.gridU;
			ul[1] = 0;
			vl[0] = self->viewer->decoration.gridV;
			vl[1] = 0;
			return Py_BuildValue("ss", ul, vl);
		} else {
			if (PyTuple_Check(value) && PyTuple_GET_SIZE(value)==2) {
			  //				char *gu = PyString_AsString(PyTuple_GetItem(value,0));
			  //				char *gv = PyString_AsString(PyTuple_GetItem(value,1));
			  const char *gu = PyUnicode_AsUTF8AndSize(PyTuple_GetItem(value,0), lllen);
			  const char *gv = PyUnicode_AsUTF8AndSize(PyTuple_GetItem(value,1), lllen);				           char ua = (gu[0]=='-'? gu[1]:gu[0]);
				char va = (gv[0]=='-'? gv[1]:gv[0]);
				self->viewer->decoration.gridU = ua;
				self->viewer->decoration.gridV = va;
			} else {
				PyErr_SetString(PyExc_TypeError, "tuple expected of size 2");
				return NULL;
			}
		}
	} else
	if (!strcmp(name,"size")) {
		if (value==NULL)
			return Py_BuildValue("ff",
				self->viewer->decoration.grid_du,
				self->viewer->decoration.grid_dv);
		else {
			PyErr_SetString(PyExc_SyntaxError, "cannot set grid size");
			return NULL;
		}
	} else
	if (!strcmp(name,"low")) {
		if (value==NULL)
			return Py_BuildValue("ff",
				self->viewer->decoration.grid_u,
				self->viewer->decoration.grid_v);
		else {
			PyErr_SetString(PyExc_SyntaxError, "cannot set grid low value");
			return NULL;
		}
	} else
	if (!strcmp(name,"axes")) {
		if (value) {
			PyErr_SetString(PyExc_SyntaxError, "cannot set grid axes value");
			return NULL;
		}

		// Grid axes direction as tuple of strings ([-][XYZUV], [-][XYZUV])
		const Matrix4& matrix = self->kernel->view.matrix();
		Vector u(matrix(0,0), matrix(1,0), matrix(2,0));
		Vector v(matrix(0,1), matrix(1,1), matrix(2,1));
		const char *ustr = "U";
		const double accuracy = 1E-6;
		if (Eq(u.x, 1.0, accuracy))
			ustr = "X";
		else
		if (Eq(u.x, -1.0, accuracy))
			ustr = "-X";

		if (Eq(u.y, 1.0, accuracy))
			ustr = "Y";
		else
		if (Eq(u.y, -1.0, accuracy))
			ustr = "-Y";

		if (Eq(u.z, 1.0, accuracy))
			ustr = "Z";
		else
		if (Eq(u.z, -1.0, accuracy))
			ustr = "-Z";

		const char *vstr = "V";
		if (Eq(v.x, 1.0, accuracy))
			vstr = "X";
		else
		if (Eq(v.x, -1.0, accuracy))
			vstr = "-X";

		if (Eq(v.y, 1.0, accuracy))
			vstr = "Y";
		else
		if (Eq(v.y, -1.0, accuracy))
			vstr = "-Y";

		if (Eq(v.z, 1.0, accuracy))
			vstr = "Z";
		else
		if (Eq(v.z, -1.0, accuracy))
			vstr = "-Z";

		return Py_BuildValue("ss", ustr, vstr);
	} else {
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid type option", name);
		return NULL;
	}
	Py_RETURN_NONE;
} // Viewer_grid

/* --- Viewer_hit3D ---
 * @param u,v	of ray direction
 * @param options	optional B(ody), R(egion), T(ype)
 * @returns hit position or a tuple with hit and options
 */
static PyObject* Viewer_hit3D(ViewerObject *self, PyObject *args)
{
	static char *regType[] = {"NORMAL", "BLCKHOLE", "LATTICE", "VOXEL"};

	double u, v;
	char *options=NULL;
	if (!PyArg_ParseTuple(args, "dd|s", &u, &v, &options)) return NULL;

	double x, y, z;
	self->kernel->view.rayPosition(u,v, &x,&y,&z);

	double dx, dy, dz;
	self->kernel->view.rayDirection(u,v, &dx,&dy,&dz);

	self->kernel->lock();
	self->kernel->engine()->incBodyCheckId();
	VZone* zone = self->kernel->engine()->whereRay(x,y,z, dx,dy,dz, SMALL3D, NULL);

	// prepare the ray
	Ray ray;
	ray.use_clip      = self->kernel->clipBodyCount() > 0;
	ray.skip_1stblack = self->viewer->d3.skip_1stblack;
	ray.push(RaySegment(x,y,z, dx,dy,dz, zone));

	self->kernel->engine()->incBodyCheckId();
	self->viewer->d3.nextIntersection(self->kernel->engine(), &ray);
	self->kernel->unlock();

//	if (ray.hitZone() && ray.hitRegion()->type() == REGION_BLACKHOLE)
//		Py_RETURN_NONE;

	Point hit = ray.hit();
	PyObject* pyhit = Py_BuildValue("ddd", hit.x, hit.y, hit.z);
	if (options==NULL) return pyhit;

	PyObject* ret = PyTuple_New(1+strlen(options));
	PyTuple_SetItem(ret, 0, pyhit);

	// B, R, T
	for (size_t i=0; i<strlen(options); i++) {
		PyObject* item = NULL;
		switch (options[i]) {
			case 'b':
			case 'B':
				if (ray.hitBody())
					item = PyString_FromString(ray.hitBody()->name());
				break;
			case 'r':
			case 'R':
				if (ray.hitZone())
					item = PyString_FromString(ray.hitRegion()->name());
				break;
			case 't':
			case 'T':
				if (ray.hitZone())
					item = PyString_FromString(regType[ray.hitRegion()->type()]);
				break;
		}
		if (item)
			PyTuple_SetItem(ret, i+1, item);
		else {
			Py_INCREF(Py_None);
			PyTuple_SetItem(ret, i+1, Py_None);
		}
	}
	return ret;
} // Viewer_hit3D

/* --- Viewer_image --- */
static PyObject* Viewer_image(ViewerObject *self, PyObject *args)
{
	char *name;
	PyObject *value=NULL;

	if (!PyArg_ParseTuple(args, "s|O", &name, &value)) return NULL;

	if (!strcmp(name,"get")) { // Return Viewport image data
		if (value!=NULL) {
			PyErr_SetString(PyExc_SyntaxError, "option 'view' only returns data");
			return NULL;
		}
		dword *pixel = self->viewer->painter.data();
		int size = Tk_Width(self->tkWin)*Tk_Height(self->tkWin);
		for (int i=0; i<size; i++, pixel++) {
			// Convert from ?RGB to ABGR
			Color32 p;
			p.val = *pixel;
			// PIL format is ABGR
			*pixel = (*pixel&FLAG_TRANSPARENT?0:0xFF000000)	// Alpha
				| ((dword)p.rgb.blue<<16)
				| ((dword)p.rgb.green<<8)
				| (dword)p.rgb.red;
		}
		return PyString_FromStringAndSize((const char *)self->viewer->painter.data(),
				self->viewer->painter.width()*
				self->viewer->painter.height()*4);
	} else
	if (!strcmp(name,"restore")) {	// restore viewport image from previous byte swapping
		dword *pixel = self->viewer->painter.data();
		int size = Tk_Width(self->tkWin)*Tk_Height(self->tkWin);
		for (int i=0; i<size; i++, pixel++) {
			// Convert from ?RGB to ABGR
			Color32 p;
			p.val = *pixel;
			// PIL format is ABGR
			*pixel =   ((dword)(p.rgb.blue)<<16)
				 | ((dword)(p.rgb.green)<<8)
				 | (p.rgb.red);
		}
	} else
	if (!strcmp(name,"data")) {	// set image data

	  /*   Aug 2024 correction
		if (!PyString_Check(value)) {
			PyErr_SetString(PyExc_SyntaxError, "string type was expected for image data");
			return NULL;
		}
		const char *str = PyUnicode_AsUTF8AndSize(value, lllen);
		dword *data = reinterpret_cast< dword*>(const_cast< char*>(str));
		int size =  strlen(str);
	  */	
	  if (!PyBytes_Check(value)) {
		PyErr_SetString(PyExc_SyntaxError, "bytes type was expected for image data");
		return NULL;
	  }
	  dword* data = reinterpret_cast<dword*>(PyBytes_AsString(value));
	  size_t size  = PyBytes_Size(value);
	  // end Aug 2024 correction
		int s = (int)(size/sizeof(dword));
		dword *pixel = data;
		for (int i=0; i<s; i++, pixel++) {
			// Convert from ?RGB to ABGR
			Color32 p;
			p.val = *pixel;
			// PIL format is ABGR
			*pixel =   ((dword)(p.rgb.blue)<<16)
				 | ((dword)(p.rgb.green)<<8)
				 | (p.rgb.red);
		}
		//if (!self->viewer->image.data(data, size)) {
		if (!self->viewer->image.data(data, 0)) {
			PyErr_SetString(PyExc_MemoryError, "not enough memory to load image");
			return NULL;
		}
	} else
	if (!strcmp(name,"size")) {	// set image size
		if (!PyTuple_Check(value) || PyTuple_Size(value)!=2) {
			PyErr_SetString(PyExc_SyntaxError, "tuple (w,h) was expected for image size");
			return NULL;
		}
		int w = (int)PyInt_AsLong(PyTuple_GetItem(value,0));
		int h = (int)PyInt_AsLong(PyTuple_GetItem(value,1));
		self->viewer->image.size(w,h);
	} else
	if (!strcmp(name,"matrix")) {	// set transformation matrix
		Matrix3 R,M;
		if (!PyTuple_Check(value) || PyTuple_Size(value)!=2) {
			PyErr_SetString(PyExc_SyntaxError, "tuple (R,M) was expected for transformation matrix");
			return NULL;
		}
		PyList_AsMatrix3(PyTuple_GetItem(value,0), R);
		PyList_AsMatrix3(PyTuple_GetItem(value,1), M);
		self->viewer->image.matrix(&R, &M);
	} else
	if (!strcmp(name,"alpha")) {	// set image transparency
		self->viewer->image.alpha((int)PyInt_AsLong(value));
	} else
	if (!strcmp(name,"level")) {	// adjust color level
		dword low=0x000000;
		dword high=0xffffff;

		if (!PyTuple_Check(value) || PyTuple_GET_SIZE(value)!=2) {
			PyErr_SetString(PyExc_TypeError, "Tuple expected with color limits");
			return NULL;
		}
		low  = (dword)PyInt_AsUnsignedLongMask(PyTuple_GetItem(value,0));
		high = (dword)PyInt_AsUnsignedLongMask(PyTuple_GetItem(value,1));
		self->viewer->image.colorRange(low, high);
	} else {
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid option", name);
		return NULL;
	}

	Py_RETURN_NONE;
} // Viewer_image

/* --- Viewer_invalid --- */
static PyObject* Viewer_invalid(ViewerObject *self)
{
	if (self->kernel->view.invalidWindow() || !self->viewer->d2.isValid())
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
} // Viewer_invalid

#if 0
/* --- Viewer_lock --- */
static PyObject* Viewer_lock(ViewerObject *self)
{
	self->viewer->lockState();
	Py_RETURN_NONE;
} // Viewer_lock

/* --- Viewer_unlock --- */
static PyObject* Viewer_unlock(ViewerObject *self)
{
	self->viewer->unlockState();
	Py_RETURN_NONE;
} // Viewer_unlock
#endif

/* --- Viewer_matrix ---
 * @param  -1   to return the inverse matrix
 * @return	return transformation martix or the inverse if @param=-1
 */
static PyObject* Viewer_matrix(ViewerObject *self, PyObject *args)
{
	PyObject *matrixObj=NULL;
	if (!PyArg_ParseTuple(args, "|O", &matrixObj)) return NULL;

	if (matrixObj != NULL) {
		if (PyInt_Check(matrixObj) && PyInt_AsLong(matrixObj)==-1) {
			// Return the inverse matrix
			return PyList_FromMatrix4(self->kernel->view.invMatrix());
		} else {
			// Set matrix
			Matrix4 matrix;
			PyList_AsMatrix4(matrixObj, matrix);
			self->viewer->matrix(matrix);
			self->projectionChanged = true;
			Py_RETURN_NONE;
		}
	} else
		// return the normal matrix
		return PyList_FromMatrix4(self->kernel->view.matrix());
} // Viewer_matrix

/* --- Viewer_memory --- */
static PyObject *Viewer_memory(ViewerObject *self, PyObject *args)
{
	char *dump=NULL;
	if (!PyArg_ParseTuple(args, "|s", &dump)) return NULL;

	if (dump==NULL)
		return PyInt_FromLong(self->viewer->memory());
	else {
		self->viewer->printMemory();
		self->kernel->engine()->printMemory();
	}
	Py_RETURN_NONE;
} // Viewer_memory

/* --- Viewer_moveUV --- *
 * using an old coordinate (x,y,z) move to the new (xn,yn,zn) but keeping the w
 * the same. Used when moving objects in the viewports
 */
static PyObject* Viewer_moveUV(ViewerObject *self, PyObject *args)
{
	double x,y,z;
	double xn,yn,zn;

	if (!PyArg_ParseTuple(args, "(ddd)(ddd)", &x, &y, &z, &xn, &yn, &zn)) return NULL;

	double u = self->kernel->view.xyz2u(xn,yn,zn);
	double v = self->kernel->view.xyz2v(xn,yn,zn);
	double w = self->kernel->view.xyz2w(x,y,z);

	return Py_BuildValue("ddd",	self->kernel->view.uvw2x(u,v,w),
					self->kernel->view.uvw2y(u,v,w),
					self->kernel->view.uvw2z(u,v,w));
} // Viewer_moveUV

/* --- Viewer_origin --- */
static PyObject* Viewer_origin(ViewerObject *self, PyObject *args)
{
	if (PyTuple_Size(args)==0) {
		double x, y, z;
		self->viewer->origin(&x,&y,&z);
		return Py_BuildValue("ddd", x, y, z);
	}

	// set center accept a tuple/list of length 3 or 3 arguments
	if (PyTuple_Size(args)==3)
		// 3 arguments passed
		self->viewer->origin(
			PyFloat_AS_DOUBLE(PyTuple_GetItem(args,0)),
			PyFloat_AS_DOUBLE(PyTuple_GetItem(args,1)),
			PyFloat_AS_DOUBLE(PyTuple_GetItem(args,2)));
	else
	if (PyTuple_Size(args)==1) {
		// 1 argument passed
		PyObject* xyz=PyTuple_GetItem(args,0);
		if (PyTuple_Check(xyz) && PyTuple_Size(xyz)==3)
			self->viewer->origin(
				PyFloat_AS_DOUBLE(PyTuple_GetItem(xyz,0)),
				PyFloat_AS_DOUBLE(PyTuple_GetItem(xyz,1)),
				PyFloat_AS_DOUBLE(PyTuple_GetItem(xyz,2)));
		else
		if (PyList_Check(xyz) && PyList_Size(xyz)==3) {
			self->viewer->origin(
				PyFloat_AS_DOUBLE(PyList_GetItem(xyz,0)),
				PyFloat_AS_DOUBLE(PyList_GetItem(xyz,1)),
				PyFloat_AS_DOUBLE(PyList_GetItem(xyz,2)));
		} else {
			PyErr_SetString(PyExc_TypeError,"function takes exactly 1 tuple of size 3 or 3 arguments");
			return NULL;
		}
	} else {
		PyErr_SetString(PyExc_TypeError,"function takes exactly 1 tuple of size 3 or 3 arguments");
		return NULL;
	}
	self->projectionChanged = true;
	Py_RETURN_NONE;
} // Viewer_origin

/* --- Viewer_offset --- */
static PyObject* Viewer_offset(ViewerObject *self, PyObject *args)
{
	if (PyTuple_Size(args)==0) {
		return Py_BuildValue("dd",
			self->kernel->view.Uofs(),
			self->kernel->view.Vofs());
	}

	// set center accept a tuple of length 3 or 3 arguments
	PyObject *uv;

	if (PyTuple_Size(args)==2) {
		self->kernel->view.offset(
				PyFloat_AS_DOUBLE(PyTuple_GetItem(args,0)),
				PyFloat_AS_DOUBLE(PyTuple_GetItem(args,1)));
	} else
	if (PyTuple_Size(args)==1 && PyTuple_Check(uv=PyTuple_GetItem(args,0))
				  && PyTuple_Size(uv)==2) {
		self->kernel->view.offset(
				PyFloat_AS_DOUBLE(PyTuple_GetItem(uv,0)),
				PyFloat_AS_DOUBLE(PyTuple_GetItem(uv,1)));
	} else {
		PyErr_SetString(PyExc_TypeError,"function takes exactly 1 tuple of size 2 or 2 arguments");
		return NULL;
	}
	self->projectionChanged = true;
	Py_RETURN_NONE;
} // Viewer_offset

/* --- Viewer_objectVar --- */
static PyObject* Viewer_objectVar(ViewerObject *self, GObject *object, const char *var, PyObject *value)
{
	if (value==NULL) {
		PyErr_SetString(PyExc_TypeError, "object move/rotate doesn't return anything.");
		return NULL;
	} else
	if (!strcmp(var, "move")) {	// Move relative from the last save position
		int opt = PyInt_AsLong(value);	// What to move
		//self->viewer->lock();
		object->node(opt, object->savedNode(opt) + self->move);
		//self->viewer->unlock();
	} else
	if (!strcmp(var, "rotate")) {	// rotate around an axis (axis,angle)
		//self._viewer.body(bodies, "rotate", (uv,  self._axis, snap))
		if (!PyTuple_Check(value) || PyTuple_GET_SIZE(value)!=3) {
			PyErr_SetString(PyExc_TypeError, "tuple of size 3 expected");
			return NULL;
		}

		// What to move
		double u,v;
		if (!Py_GetUV(PyTuple_GetItem(value, 0), &u, &v)) return NULL;
		Vector axis = Py_GetVector(PyTuple_GetItem(value, 1));
		int snap = PyInt_AsLong(PyTuple_GetItem(value, 2));
		if (PyErr_Occurred()) return NULL;

		// Calculate angle
		double angle = atan2(v-self->pivotV, u-self->pivotU) - self->pivotAngle;
		if (snap) // FIXME dangle should be past from python
			angle = (double)Round(angle / self->geometry->snapAngle)
				* self->geometry->snapAngle;

		object->rotate(angle,axis);
	} else {
		PyErr_Format(PyExc_TypeError, "Invalid type \'%s\' specified", var);
		return NULL;
	}
	Py_RETURN_NONE;
} // Viewer_objectVar

/* --- Viewer_object --- */
static PyObject* Viewer_object(ViewerObject *self, PyObject *args)
{
	PyObject *obj;
	GObject *object = NULL;
	char const *var = (char const *)"id";
	PyObject *value=NULL;
	if (!PyArg_ParseTuple(args, "O|sO", &obj, &var, &value)) return NULL;

	if (Py_Check4Pattern(obj)) {
	  //		char *pattern = PyString_AsString(obj);
	  const char *pattern = PyUnicode_AsUTF8AndSize(obj, lllen);

		ArrayIterator<GObject*> iter(*self->geometry->objectList);
		while (iter) {
			object = iter++;
			if (!fnmatch(pattern, object->name(), 0)) {
				PyObject *ret = Viewer_objectVar(self, object, var, value);
				Py_XDECREF(ret);	// ret can be NULL use XDECREF
			}
		}
		Py_RETURN_NONE;
	} else
	if (PyList_Check(obj)) {
		for (ssize_t i=0; i<PyList_GET_SIZE(obj); i++) {
			object = Py_Object(self->geometry, PyList_GetItem(obj,i));
			if (object==NULL) return NULL;
			PyObject *ret = Viewer_objectVar(self, object, var, value);
			Py_XDECREF(ret);	// ret can be NULL use XDECREF
		}
		Py_RETURN_NONE;
	} else
	if (PyTuple_Check(obj)) {
		for (ssize_t i=0; i<PyTuple_GET_SIZE(obj); i++) {
			object = Py_Object(self->geometry, PyTuple_GetItem(obj,i));
			if (object==NULL) return NULL;
			PyObject *ret = Viewer_objectVar(self, object, var, value);
			Py_XDECREF(ret);	// ret can be NULL use XDECREF
		}
		Py_RETURN_NONE;
	} else {
		object = Py_Object(self->geometry, obj);
		if (object)
			return Viewer_objectVar(self, object, var, value);
	}
	return NULL;
} // Viewer_object

/* --- Viewer_pixel2uv --- */
static PyObject* Viewer_pixel2uv(ViewerObject *self, PyObject *args)
{
	int i,j;
	int centered=0;

	if (!PyArg_ParseTuple(args, "ii|i", &i, &j, &centered)) return NULL;

	if (centered)
		// WARNING: the uv of the viewer is always centered on the window
		// while on the viewer could be shifted...
		return Py_BuildValue("dd",	self->kernel->view.i2u(i)-self->kernel->view.Uofs(),
						self->kernel->view.j2v(j)-self->kernel->view.Vofs());
	else
		return Py_BuildValue("dd",	self->kernel->view.i2u(i),
						self->kernel->view.j2v(j));
} // Viewer_pixel2uv

/* --- Viewer_pixel2xyz --- */
static PyObject* Viewer_pixel2xyz(ViewerObject *self, PyObject *args)
{
	int i,j;

	if (!PyArg_ParseTuple(args, "ii", &i, &j)) return NULL;

	double x, y, z;
	self->kernel->view.ij2xyz(i,j,&x,&y,&z);

	return Py_BuildValue("ddd", x,y,z);
} // Viewer_pixel2xyz

/* --- Viewer_duv2dxyz --- */
static PyObject* Viewer_duv2dxyz(ViewerObject *self, PyObject *args)
{
	double du,dv;

	if (!PyArg_ParseTuple(args, "dd", &du, &dv)) return NULL;

	double dx = self->kernel->view.uv2dx(du,dv);
	double dy = self->kernel->view.uv2dy(du,dv);
	double dz = self->kernel->view.uv2dz(du,dv);

	return Py_BuildValue("ddd", dx,dy,dz);
} // Viewer_duv2dxyz

/* --- Viewer_palette --- */
static PyObject* Viewer_palette(ViewerObject *self, PyObject *args)
{
	int	 id;
	char	*name;
	PyObject *value=NULL;

	if (!PyArg_ParseTuple(args, "is|O", &id, &name, &value)) return NULL;

	if (id<0 || id>= MAXPALETTE) {
		PyErr_SetString(PyExc_ValueError,"Invalid palette index");
		return NULL;
	}

	if (!strcmp(name,"alphamin")) {
		if (value)
			self->viewer->palette[id].alphamin((bool)PyInt_AsLong(value));
		else
			return PyBool_FromLong(self->viewer->palette[id].alphamin());
	} else
	if (!strcmp(name,"alphamax")) {
		if (value)
			self->viewer->palette[id].alphamax((bool)PyInt_AsLong(value));
		else
			return PyBool_FromLong(self->viewer->palette[id].alphamax());
	} else
	if (!strcmp(name,"default")) {
		if (value)
			self->viewer->palette.def(id);
		else
			return PyBool_FromLong(self->viewer->palette.def());
	} else
	if (!strcmp(name,"interpolate")) {
		if (value)
			self->viewer->palette[id].interpolate((bool)PyInt_AsLong(value));
		else
			return PyBool_FromLong(self->viewer->palette[id].interpolate());
	} else
	if (!strcmp(name,"invert")) {
		if (value)
			self->viewer->palette[id].invert((bool)PyInt_AsLong(value));
		else
			return PyBool_FromLong(self->viewer->palette[id].invert());
	} else
	if (!strcmp(name,"label")) {
		if (value)
		  //			self->viewer->palette.label(id, PyString_AsString(value));
		  self->viewer->palette.label(id, PyUnicode_AsUTF8AndSize(value, lllen));
		else
			return PyString_FromString(self->viewer->palette.label(id).c_str());
	} else
	if (!strcmp(name,"log")) {
		if (value)
			self->viewer->palette[id].log((bool)PyInt_AsLong(value));
		else
			return PyBool_FromLong(self->viewer->palette[id].log());
	} else
	if (!strcmp(name,"max")) {
		if (value)
			self->viewer->palette[id].max(PyFloat_AS_DOUBLE(value));
		else
			return PyFloat_FromDouble(self->viewer->palette[id].max());
	} else
	if (!strcmp(name,"min")) {
		if (value)
			self->viewer->palette[id].min(PyFloat_AS_DOUBLE(value));
		else
			return PyFloat_FromDouble(self->viewer->palette[id].min());
	} else
	if (!strcmp(name,"palette")) {
		if (value) {
			if (Py_TYPE(value) != &PyList_Type) {
				PyErr_SetString(PyExc_TypeError, "Invalid type, list expected");
				return NULL;
			}
			size_t n = PyList_Size(value);
			if (n>_MAXCOLORS) {
				PyErr_SetString(PyExc_TypeError, "Maximum number of colors accepted is 256");
				return NULL;
			}
			self->viewer->palette[id].size((int)n);
			for (size_t i=0; i<n; i++)
				self->viewer->palette[id][i] = \
					(dword)PyInt_AsUnsignedLongMask(PyList_GetItem(value, i));
		}
	} else
	if (!strcmp(name,"reset")) {
		self->viewer->palette[id].reset();
		self->viewer->palette.display(id, false);
	} else
	if (!strcmp(name,"show")) {
		if (value)
			self->viewer->palette.display(id, (bool)PyInt_AsLong(value));
		else
			return PyBool_FromLong(self->viewer->palette.display(id));
	} else
	if (!strcmp(name,"smooth")) {
		if (value)
			self->viewer->palette[id].interpolate((bool)PyInt_AsLong(value));
		else
			return PyBool_FromLong(self->viewer->palette[id].interpolate());
	} else {
		// Error
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid type option", name);
		return NULL;
	}
	self->viewer->palette[id].init();	// calc palette

	Py_RETURN_NONE;
} // Viewer_palette

/** Viewer_pen
 * process handrawn array 'pen'
 */
static PyObject* Viewer_pen(ViewerObject *self, PyObject *args)
{
	char	*name;
	PyObject *value=NULL;

	if (!PyArg_ParseTuple(args, "s|O", &name, &value)) return NULL;

	if (!strcmp(name,"add")) {
		// expect pixel (i,j) coordinates
		if (!PyTuple_Check(value) || PyTuple_GET_SIZE(value)!=2) {
			PyErr_SetString(PyExc_TypeError, "pixel tuple of size 2 expected");
			return NULL;
		}
		self->pen->add(IPoint(PyInt_AsLong(PyTuple_GetItem(value,0)),
					PyInt_AsLong(PyTuple_GetItem(value,1))));
	} else
	if (!strcmp(name,"clear")) {
		self->pen->clear();
	} else {
		// Error
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid type option", name);
		return NULL;
	}

	Py_RETURN_NONE;
} // Viewer_pen

/* _endProjection: notifier function for end of projection */
static void _endProjection(ViewerObject *self)
{
	DUMP(cout << "_endProjection("<<self->viewer->title()<<")" << endl);
	_sendMessage(self, "EndProjection");
} // _endProjection

/* --- Viewer_project --- */
static PyObject* Viewer_project(ViewerObject *self, PyObject *args)
{
	int asThread = false;
	int all = true;
	if (!PyArg_ParseTuple(args, "|ii", &asThread, &all)) return NULL;
	self->viewer->d2.projectAll(all);
	_clearErrors(self);
	DUMP(cout << "Viewer_project("<<self->viewer->title()<<", " << asThread << ")" << endl);
	if (asThread) {
		self->viewer->spawnProject((NotifyFunc)_endProjection, self);
	} else {
		// Stop any pending thread
		DUMP(cout << "Viewer_project((" << self->viewer->title() << ").stopThread()" << endl);
		self->viewer->stopThread();
		// Run synchronously
		self->viewer->d2.project();
	}
	Py_RETURN_NONE;
} // Viewer_project

/* --- Viewer_projectionChanged --- */
static PyObject* Viewer_projectionChanged(ViewerObject *self, PyObject *args)
{
	int flag=-1;

	if (!PyArg_ParseTuple(args, "|i", &flag)) return NULL;
	if (flag == -1)
		return PyBool_FromLong(self->projectionChanged);
	else
		self->projectionChanged = (bool)flag;
	Py_RETURN_NONE;
} // Viewer_projectionChanged

/* --- Viewer_rectangle --- */
static PyObject* Viewer_rectangle(ViewerObject *self, PyObject *args)
{
	if (!PyArg_ParseTuple(args, "iiii",
			&(self->rectX1), &(self->rectY1),
			&(self->rectX2), &(self->rectY2))) return NULL;
	Py_RETURN_NONE;
} // Viewer_rectangle

/* --- Viewer_Sx --- */
static PyObject* Viewer_Sx(ViewerObject *self)
{
	return PyFloat_FromDouble(self->kernel->view.Sx());
} // Viewer_Sx

/* --- Viewer_Sy --- */
static PyObject* Viewer_Sy(ViewerObject *self)
{
	return PyFloat_FromDouble(self->kernel->view.Sy());
} // Viewer_Sy

/* --- Viewer_snap --- */
static PyObject* Viewer_snap(ViewerObject *self, PyObject *args)
{
	double u, v, d=self->geometry->snapDistance;

	if (!PyArg_ParseTuple(args, "dd|d", &u, &v, &d)) return NULL;

	return Py_Vector(_snapUV(self, u, v, d));
} // Viewer_snap

/* --- Viewer_set --- */
static PyObject* Viewer_set(ViewerObject *self, PyObject *args)
{
#define NOTDEFINED	-999999999
	char	*name;
	int	 value = NOTDEFINED;

	if (!PyArg_ParseTuple(args, "s|i", &name, &value)) return NULL;

	if (!strcmp(name,"3D")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->d3.show);
		else
			self->viewer->d3.show = (bool)value;
	} else
	if (!strcmp(name,"axes")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->decoration.showAxes);
		else
			self->viewer->decoration.showAxes = value;
	} else
	if (!strcmp(name,"antialias")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->viewer->d3.antialias);
		else
			self->viewer->d3.antialias = Range(1,value,16);
	} else
	if (!strcmp(name,"ambient")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->viewer->d3.ambient);
		else
			self->viewer->d3.ambient= value & 0xFF;
	} else
	if (!strcmp(name,"backgroundcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->geometry->_backgroundColor & 0xFFFFFF);
		else {
			self->viewer->backgroundColor(value);
		}
	} else
	if (!strcmp(name,"borders") || !strcmp(name,"2D")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->d2.showBorders);
		else
			self->viewer->d2.showBorders = (bool)value;
	} else
	if (!strcmp(name,"clipbody")) {
		if (value==NOTDEFINED) {
//			if (self->kernel->clipbody())
//				return PyInt_FromLong(self->kernel->clipbody());
//			else
//				return PyInt_FromLong(-1);
		} else
		if (value<0)
			self->kernel->clipBodyClear();
		else
			self->kernel->clipBodyAdd(value);
	} else
	if (!strcmp(name,"clipnegative")) {
//		if (value==NOTDEFINED)
//			return PyBool_FromLong(self->kernel->clipNegative());
//		else
			self->kernel->clipBodyNegative(value);
	} else
	if (!strcmp(name,"crosshair")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->crosshair);
		else
			self->crosshair = value;
	} else
	if (!strcmp(name,"cores")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->kernel->nthreads());
		else
			self->kernel->initThreads(value);
	}  else
	if (!strcmp(name,"deflights")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->d3.deflights);
		else
			self->viewer->d3.deflights = (bool)value;
	} else
	if (!strcmp(name,"depth"))
		return PyInt_FromLong(self->depth);
	else
	if (!strcmp(name,"drawtime")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->viewer->d3.drawTime());
		else
			self->viewer->d3.drawTime(value);
	} else
	if (!strcmp(name,"edgedetect")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->d3.drawEdges);
		else
			self->viewer->d3.drawEdges = (bool)value;
	} else
	if (!strcmp(name,"errors")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->showErrors);
		else
			self->viewer->showErrors = value;
	} else
	if (!strcmp(name,"fill")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->d2.fillRegions);
		else {
			self->viewer->d2.fillRegions = (bool)value;
			if ((bool)value) self->viewer->d2.showBorders = true;
		}
	} else
	if (!strcmp(name,"image")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->image.show);
		else
			self->viewer->image.show = value;
	} else
	if (!strcmp(name,"grid")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->decoration.showGrid);
		else
			self->viewer->decoration.showGrid = value;
	} else
	if (!strcmp(name,"gridlevel")) {
		// Old values was from 0..32 (5bits) now is 8bits
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->viewer->decoration.gridLevel>>3);
		else
			self->viewer->decoration.gridLevel = value<<3;
	} else
	if (!strcmp(name,"latticelevel")) {
		// Old values was from 0..32 (5bits) now is 8bits
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->viewer->lattice.level>>3);
		else {
			self->viewer->lattice.level = value<<3;
			self->viewer->lattice.makeHashColor();
		}
	} else
	if (!strcmp(name,"labels")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->d2.showLabels);
		else
			self->viewer->d2.showLabels = (bool)value;
	} else
	if (!strcmp(name,"lattice")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->lattice.show);
		else
			self->viewer->lattice.show = (bool)value;
	} else
	if (!strcmp(name,"palette")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->palette.show);
		else
			self->viewer->palette.show = (bool)value;
	} else
	if (!strcmp(name,"projection")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->kernel->view.projection);
		else
			self->kernel->view.projection = (ProjectionType)value;
	} else
	if (!strcmp(name,"projectbody")) {
		if (value==NOTDEFINED) {
			// return a list of projection bodies
			PyObject* obj = PyTuple_New(self->kernel->projectBodyCount());
			for (int i=0; i<self->kernel->projectBodyCount(); i++)
				PyTuple_SetItem(obj, i, PyInt_FromLong(self->kernel->projectBody(i)));
			return obj;
		} else
		if (value<0)
			self->kernel->projectBodyClear();
		else
			self->kernel->projectBodyAdd(value);
	} else
	if (!strcmp(name,"reflections")) {
		if (value==NOTDEFINED) {
			if (self->viewer->d3.maxDepth())
				return PyInt_FromLong(self->viewer->d3.maxDepth());
			else
				return PyInt_FromLong(-1);
		} else
			self->viewer->d3.maxDepth(value);
	} else
	if (!strcmp(name,"shadows")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->d3.shadows);
		else
			self->viewer->d3.shadows = (bool)value;
	} else
	if (!strcmp(name,"skip1stblack")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->d3.skip_1stblack);
		else
			self->viewer->d3.skip_1stblack = (bool)value;
	} else
	if (!strcmp(name,"textbackground")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->textBackgroundLevel);
		else
			self->viewer->textBackgroundLevel = value;
	} else
	if (!strcmp(name,"title")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->decoration.showTitle);
		else
			self->viewer->decoration.showTitle = value;
	} else
	if (!strcmp(name,"trackball")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->showTrackball);
		else
			self->showTrackball = (bool)value;
	} else
	if (!strcmp(name,"usrbinastexture")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->d3.usrbinAsTexture());
		else
			self->viewer->d3.usrbinAsTexture((bool)value);
	} else
	if (!strcmp(name,"userdump")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->userdump.show);
		else
			self->viewer->userdump.show= (bool)value;
	} else
	if (!strcmp(name,"usrbin")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->usrbin.show);
		else
			self->viewer->usrbin.show = (bool)value;
	} else
	if (!strcmp(name,"vertex")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->d2.showVertex);
		else
			self->viewer->d2.showVertex = (bool)value;
	} else
	if (!strcmp(name,"viewport")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->showViewport);
		else
			self->showViewport = (bool)value;
	} else
	if (!strcmp(name,"voxel")) {
		if (value==NOTDEFINED)
			return PyBool_FromLong(self->viewer->voxel.show);
		else
			self->viewer->voxel.show = (bool)value;
	} else
	if (!strcmp(name,"voxellevel")) {
		// Old values was from 0..32 (5bits) now is 8bits
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->viewer->voxel.level>>3);
		else {
			self->viewer->voxel.level = value<<3;
			self->viewer->voxel.makeHashColor();
		}
	} else
	if (!strcmp(name,"xray")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->viewer->d3.xray);
		else
			self->viewer->d3.xray= Range(0,value,255);
	} else {
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid type option", name);
		return NULL;
	}
	Py_RETURN_NONE;
} // Viewer_set

/* --- Viewer_regionVar --- */
static PyObject* Viewer_regionVar(ViewerObject *self, VRegion *region, const char *var, PyObject *value)
{
	if (!strcmp(var, "id")) {
		if (value==NULL)
			return PyInt_FromLong(region->id());
		else {
			PyErr_SetString(PyExc_SyntaxError, "Cannot set region id");
			return NULL;
		}
	} else
	if (!strcmp(var, "color")) {
		if (value==NULL)
			return PyInt_FromLong(region->color()&0xFFFFFF);
		else {
			region->value(0.0);	// FIXME Do I need it here ????
			region->color((dword)PyInt_AsUnsignedLongMask(value));
		}
	} else
	if (!strcmp(var, "value")) {
		if (value==NULL)
			return PyFloat_FromDouble(region->value());
		else {
			int pal = self->viewer->palette.def();
			// set color from the palette using the value
			region->value(PyFloat_AsDouble(value));
			region->color(self->viewer->palette[pal].color(region->value()));
		}
	} else
	if (!strcmp(var, "label")) {
		if (value==NULL)
			return PyString_FromString(region->label);
		else {
		  //			strncpy(region->label, PyString_AsString(value),
		  strncpy(region->label, PyUnicode_AsUTF8AndSize(value, lllen),			sizeof(region->label));
			region->label[sizeof(region->label)-1] = 0;
		}
	} else
	if (!strcmp(var, "alpha")) {
		if (value==NULL)
			return PyInt_FromLong(region->alpha()&0xFF);
		else
			region->alpha(PyInt_AsLong(value));
	} else

	if (!strcmp(var, "correct")) {
		if (value==NULL)
			return PyBool_FromLong(self->kernel->correctOverlaps(region));
		else if (PyList_Check(value)) {
			if (PyList_GET_SIZE(value)==0)
				return PyBool_FromLong(self->kernel->correctOverlaps(region));

			bool modified = false;
			for (ssize_t i=0; i<PyList_GET_SIZE(value); i++) {
				int zoneid =  PyInt_AsLong(PyList_GetItem(value,i));
				modified |= self->kernel->correctOverlaps(region, zoneid);
			}
			return PyBool_FromLong(modified);
		}
	} else {
		PyErr_Format(PyExc_TypeError, "Invalid type \'%s\' specified", var);
		return NULL;
	}
	Py_RETURN_NONE;
} // Viewer_regionVar

/* --- Viewer_region --- */
static PyObject* Viewer_region(ViewerObject *self, PyObject *args)
{
	PyObject *obj;
	char *var;
	PyObject *value=NULL;
	if (!PyArg_ParseTuple(args, "Os|O", &obj, &var, &value)) return NULL;

	if (Py_Check4Pattern(obj)) {
	  //		char *pattern = PyString_AsString(obj);
              const char *pattern = PyUnicode_AsUTF8AndSize(obj, lllen);
		for (int ir=0; ir<self->kernel->nGeometryRegions(); ir++) {
			VRegion *region = self->kernel->regions[ir];
			if (!fnmatch(pattern, region->name(), 0)) {
				PyObject *ret = Viewer_regionVar(self, region, var, value);
				Py_XDECREF(ret);	// ret can be NULL use XDECREF
			}
		}
		Py_RETURN_NONE;
	} else  {
		VRegion *region = Py_VRegion(self, obj);
		if (region==NULL) return NULL;
		return Viewer_regionVar(self, region, var, value);
	}
	Py_RETURN_NONE;		// To keep compiler happy
} // Viewer_region

/* --- Viewer_size --- */
static PyObject* Viewer_size(ViewerObject *self)
{
	return Py_BuildValue("ii", self->viewer->width(), self->viewer->height());
} // Viewer_size

/* --- Viewer_stopThread --- */
static PyObject* Viewer_stopThread(ViewerObject *self)
{
	DUMP(cout << "Viewer_stopThread("<<self->viewer->title()<<")" << endl);
	self->viewer->stopThread();
	Py_RETURN_NONE;
} // Viewer_stopThread

/* --- Viewer_title --- */
static PyObject* Viewer_title(ViewerObject *self, PyObject *args)
{
	char *title = NULL;
	if (!PyArg_ParseTuple(args, "|s", &title)) return NULL;
	if (title==NULL)
		return PyString_FromString(self->viewer->title());
	else
		self->viewer->title(title);
	Py_RETURN_NONE;
} // Viewer_title

/* --- Viewer_state --- */
static PyObject* Viewer_state(ViewerObject *self)
{
	return PyInt_FromLong(self->viewer->state());
} // Viewer_state

/* --- Viewer_userdump --- */
static PyObject* Viewer_userdump(ViewerObject *self, PyObject *args)
{
	char	 *name;
	PyObject *value=NULL;
	PyObject *opt=NULL;

	if (!PyArg_ParseTuple(args, "s|OO", &name, &value, &opt))
		return NULL;

	if (!strcmp(name,"file")) {	// load usrbin detector
		if (value!=NULL) {
		  //			self->viewer->userdump.open(PyString_AsString(value));
		  self->viewer->userdump.open(PyUnicode_AsUTF8AndSize(value, lllen));
		}
	} else
	if (!strcmp(name,"start")) {	// start event
		if (value!=NULL)
			self->viewer->userdump.start = PyInt_AsLong(value);
		else
			return PyInt_FromLong(self->viewer->userdump.start);
	} else
	if (!strcmp(name,"n")) {	// # events
		if (value!=NULL)
			self->viewer->userdump.n = PyInt_AsLong(value);
		else
			return PyInt_FromLong(self->viewer->userdump.n);
	} else
	if (!strcmp(name,"reset")) {	// # reset all information on particles
		self->viewer->userdump.reset();
	} else
	if (!strcmp(name,"show")) {	// particle to hide
		int n = PyInt_AsLong(value);
		self->viewer->userdump.display(n);
	} else
	if (!strcmp(name,"hide")) {	// particle to hide
		int n = PyInt_AsLong(value);
		self->viewer->userdump.hide(n);
	} else
	if (!strcmp(name,"color")) {	// color of particle
		int n = PyInt_AsLong(value);
		self->viewer->userdump.color(n, PyInt_AsLong(opt));
	} else {
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid option", name);
		return NULL;
	}

	Py_RETURN_NONE;
} // Viewer_userdump

/* --- Viewer_usrbin --- */
static PyObject* Viewer_usrbin(ViewerObject *self, PyObject *args)
{
	char *name;
	int index;
	PyObject *value=NULL, *opt=NULL;

	if (!PyArg_ParseTuple(args, "is|OO", &index, &name, &value, &opt))
		return NULL;

	if (!InRange(0,index,MAXUSRBIN)) {
		PyErr_SetString(PyExc_ValueError, "Invalid usrbin index");
		return NULL;
	}

	Usrbin &usrbin = self->viewer->usrbin[index];

	if (!strcmp(name,"load")) {	// load usrbin detector
		if (value != NULL && opt != NULL) {
		  //			char *filename = PyString_AsString(value);
		  const char *filename = PyUnicode_AsUTF8AndSize(value,lllen);
			int det = (int)PyInt_AsLong(opt);
			if (!usrbin.load(filename, det)) {
				PyErr_Format(PyExc_IOError, "Unable to load usrbin detector \'%s\' %d",
						filename, det);
				return NULL;
			}

//			return PyString_FromFormat(
//				"Usrbin: %s\n"
//				"Title:  %s\n"
//				"Time:   %s\n"
//				"Detector: %d %s\n"
//				"\tScore Type=%d Particle=%d\n"
//				"\tX: %.15g\t%.15g\t%d\t%.15g\n"
//				"\tY: %.15g\t%.15g\t%d\t%.15g\n"
//				"\tZ: %.15g\t%.15g\t%d\t%.15g\n",
//				filename,
//				usrbin.runtitle,
//				usrbin.runtime,
//				usrbin.detector(), usrbin.titusb,
//				usrbin.type(), usrbin.score(),
//				usrbin.xlow, usrbin.xhigh, usrbin.nx, usrbin.dx,
//				usrbin.ylow, usrbin.yhigh, usrbin.ny, usrbin.dy,
//				usrbin.zlow, usrbin.zhigh, usrbin.nz, usrbin.dz);
		}
	} else
	if (!strcmp(name,"alpha")) {	// set usrbin transparency
		if (value==NULL)
			return PyInt_FromLong(self->viewer->usrbin.alpha(index));
		else
			self->viewer->usrbin.alpha(index, (int)PyInt_AsLong(value));
	} else
	if (!strcmp(name,"palette")) {
		if (value==NULL)
			return PyInt_FromLong(self->viewer->usrbin.palette(index));
		else
			self->viewer->usrbin.palette(index, (int)PyInt_AsLong(value));
	} else
	if (!strcmp(name,"cleanup"))
		usrbin.cleanup();
	else
	if (!strcmp(name,"norm")) {	// set normalization factor
		if (value==NULL)
			return PyFloat_FromDouble((double)usrbin.norm());
		else
			usrbin.norm((float)PyFloat_AsDouble(value));
	} else
	if (!strcmp(name,"log")) {	// set normalization factor
		if (value==NULL)
			return PyBool_FromLong(usrbin.logscale());
		else
			usrbin.convert((bool)PyInt_AsLong(value));
	} else
	if (!strcmp(name,"matrix")) {	// set transformation matrix
		if (value) {
			if (PyList_Check(value)) {
				Matrix4	matrix;
				PyList_AsMatrix4(value, matrix);
				usrbin.matrix(matrix);
			} else {
				usrbin.clearMatrix();
			}
		} //else return matrix
	} else
	if (!strcmp(name,"offset")) {	// set offset
		if (value)
			self->viewer->usrbin[index].offset(Py_GetVector(value));
	} else
	if (!strcmp(name,"regioncolor")) {
		self->viewer->usrbin.regionColorFromUsrbin(index);
	} else
	if (!strcmp(name,"checker")) {	// set checker usrbin
		if (value!=NULL) {
			if (PyTuple_Check(value) && PyTuple_GET_SIZE(value)==10) {
				UsrbinType it = (UsrbinType)PyInt_AsLong(PyTuple_GetItem(value, 0));
				double xl = PyFloat_AsDouble(PyTuple_GetItem(value, 1));
				double xh = PyFloat_AsDouble(PyTuple_GetItem(value, 2));
				int    nx = PyInt_AsLong(PyTuple_GetItem(value, 3));

				double yl = PyFloat_AsDouble(PyTuple_GetItem(value, 4));
				double yh = PyFloat_AsDouble(PyTuple_GetItem(value, 5));
				int    ny = PyInt_AsLong(PyTuple_GetItem(value, 6));

				double zl = PyFloat_AsDouble(PyTuple_GetItem(value, 7));
				double zh = PyFloat_AsDouble(PyTuple_GetItem(value, 8));
				int    nz = PyInt_AsLong(PyTuple_GetItem(value, 9));
				usrbin.checker(it, xl, xh, nx, yl, yh, ny, zl, zh, nz);
			} else
				PyErr_SetString(PyExc_SyntaxError, "for 'checker' a tuple of length 10 expected");
		}
	} else
	if (!strcmp(name,"type")) {
		if (value==NULL)
			return PyInt_FromLong(usrbin.type());
	} else
	if (!strcmp(name,"score")) {
		if (value==NULL)
			return PyInt_FromLong(usrbin.score());
	} else
	if (!strcmp(name,"bins")) {
		if (value==NULL)
			return Py_BuildValue("ii",
					usrbin.nx,
					usrbin.ny,
					usrbin.nz);
	} else
	if (!strcmp(name,"low")) {
		if (value==NULL)
			return Py_BuildValue("ddd",
					usrbin.xlow,
					usrbin.ylow,
					usrbin.zlow);
	} else
	if (!strcmp(name,"high")) {
		if (value==NULL)
			return Py_BuildValue("ddd",
					usrbin.xhigh,
					usrbin.yhigh,
					usrbin.zhigh);
	} else
	if (!strcmp(name,"min")) {
		if (value==NULL)
			return PyFloat_FromDouble((double)usrbin.min);
	} else
	if (!strcmp(name,"max")) {
		if (value==NULL)
			return PyFloat_FromDouble((double)usrbin.max);
	} else
	if (!strcmp(name,"value")) {
		if (value!=NULL &&
		    usrbin.hasData() &&
		    PyTuple_Check(value) &&
		    PyTuple_GET_SIZE(value)==3) {
			double x = PyFloat_AsDouble(PyTuple_GetItem(value, 0));
			double y = PyFloat_AsDouble(PyTuple_GetItem(value, 1));
			double z = PyFloat_AsDouble(PyTuple_GetItem(value, 2));
			bool ok;
			double val = usrbin.get(x, y, z, &ok);
			if (ok) return PyFloat_FromDouble(val);
		}
	} else {
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid option", name);
		return NULL;
	}

	Py_RETURN_NONE;
} // Viewer_usrbin

/* --- Viewer_uv2pixel --- */
static PyObject* Viewer_uv2pixel(ViewerObject *self, PyObject *args)
{
	double u,v;

	if (!PyArg_ParseTuple(args, "dd", &u, &v)) return NULL;
	return Py_BuildValue("ii",	self->kernel->view.u2i(u),
					self->kernel->view.v2j(v));
} // Viewer_uv2pixel

/* --- Viewer_uv2xyz --- */
static PyObject* Viewer_uv2xyz(ViewerObject *self, PyObject *args)
{
	double u,v;

	if (!PyArg_ParseTuple(args, "dd", &u, &v)) return NULL;

// FIXME u,v should be relative to the center of the screen!!!
	double x = self->kernel->view.uv2x(u,v);
	double y = self->kernel->view.uv2y(u,v);
	double z = self->kernel->view.uv2z(u,v);

	return Py_BuildValue("ddd", x,y,z);
} // Viewer_uv2xyz

/* --- Viewer_vertex --- */
static PyObject* Viewer_vertex(ViewerObject *self, PyObject *args)
{
	char *axis;
	if (!PyArg_ParseTuple(args, "s", &axis)) return NULL;
	Array<double> vertices;

	if (('u'<=axis[0] && axis[0]<='z') || ('U'<=axis[0] && axis[0]<='Z')) {
		if (!self->viewer->d2.projectVertices(axis[0], vertices)) {
			PyErr_SetString(PyExc_SyntaxError, "Error getting vertices");
			return NULL;
		}
		double prev = -MAX_REAL;
		PyObject* vertexList = PyList_New(0);
		for (int i=0; i<vertices.size(); i++) {
			double& v = vertices[i];
			if (Abs(v-prev)>SMALL5)
				PyList_Append(vertexList, PyFloat_FromDouble(v));
			prev = v;
		}
		return vertexList;
	} else {
		PyErr_SetString(PyExc_ValueError, "Axis string expected as argument");
		return NULL;
	}
} // Viewer_vertex

/* --- Viewer_viewOto0 --- *
 * Move view origin to center of "view"
 */
static PyObject* Viewer_viewOto0(ViewerObject *self)
{
	self->viewer->moveViewOriginTo0();
	Py_RETURN_NONE;
} // Viewer_viewOto0

/* --- Viewer_viewport --- */
static PyObject* Viewer_viewport(ViewerObject *self, PyObject *args)
{
	int id;
	char *name;
	PyObject *value=NULL;
	if (!PyArg_ParseTuple(args, "is|O", &id, &name, &value)) return NULL;
	if (id<0 || id>=NVIEWS) {
		PyErr_SetString(PyExc_KeyError, "Viewport not found");
		return NULL;
	}
	if (!strcmp(name, "center")) {
		if (value != NULL) {
			if (!PyTuple_Check(value)) {
				PyErr_SetString(PyExc_TypeError, "Tuple expected");
				return NULL;
			}
			double x = PyFloat_AS_DOUBLE(PyTuple_GetItem(value, 0));
			double y = PyFloat_AS_DOUBLE(PyTuple_GetItem(value, 1));
			double z = PyFloat_AS_DOUBLE(PyTuple_GetItem(value, 2));

			self->viewport[id].viewer->origin(x,y,z);
		} else {
			double x, y, z;
			self->viewport[id].viewer->origin(&x,&y,&z);
			return Py_BuildValue("ddd",x,y,z);
		}
	} else
	if (!strcmp(name, "centerview")) {
		if (value != NULL) {
			PyErr_SetString(PyExc_TypeError, "centerview is readonly");
			return NULL;
		} else
			return Py_BuildValue("dd",
					self->viewport[id].uc,
					self->viewport[id].vc);
	} else
	if (!strcmp(name, "centerpixel")) {
		if (value != NULL) {
			PyErr_SetString(PyExc_TypeError, "centerpixel is readonly");
			return NULL;
		} else
			return Py_BuildValue("ii",
					Round(self->viewport[id].xc),
					Round(self->viewport[id].yc));
	} else
	if (!strcmp(name, "viewer")) {
		if (value != NULL) {
			if (Py_TYPE(value) != &ViewerType) {
				PyErr_SetString(PyExc_TypeError, "Invalid type, Viewer expected");
				return NULL;
			}
			self->viewport[id].viewer = ((ViewerObject*)value)->viewer;
		} else
			return PyList_FromMatrix4(self->viewport[id].viewer->view().matrix());
	} else if (!strcmp(name, "linewidth")) {
		if (value != NULL)
			self->viewport[id].lineWidth = (int)PyInt_AsLong(value);
		else
			return PyInt_FromLong(self->viewport[id].lineWidth);

	} else if (!strcmp(name, "originwidth")) {
		if (value != NULL)
			self->viewport[id].originWidth = (int)PyInt_AsLong(value);
		else
			return PyInt_FromLong(self->viewport[id].originWidth);

	} else if (!strcmp(name, "wlength")) {
		if (value != NULL)
			self->viewport[id].wLength = PyFloat_AsDouble(value);
		else
			return PyFloat_FromDouble(self->viewport[id].wLength);

	} else if (!strcmp(name, "color")) {
		if (value != NULL)
			self->viewport[id].color = (int)PyInt_AsLong(value);
		else
			return PyInt_FromLong(self->viewport[id].color);
	} else {
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid type option", name);
		return NULL;
	}
	Py_RETURN_NONE;
} // Viewer_viewport

/* --- Viewer_volume --- */
static PyObject *Viewer_volume(ViewerObject *self, PyObject *args)
{
	PyObject *obj;
	int samples;
	if (!PyArg_ParseTuple(args, "Oi", &obj, &samples)) return NULL;
	if (Py_Check4Pattern(obj)) {
/*
		char *pattern = PyString_AsString(obj);
		ArrayIterator<GRegion*> iter(self->geometry->regions);
		while (iter) {
			region = iter++;
			if (!fnmatch(pattern, region->name(), 0)) {
				PyObject ret = Geometry_regionVar(self, region, var, value);
				Py_XDECREF(ret);	// ret can be NULL use XDECREF
			}
		}
*/
	} else  {
//		GeometryKernel kernel(self->geometry);
//		kernel.initThreads(0);
		GRegion* region = Py_VRegion(self, obj)->region();
		if (region==NULL) return NULL;
		double vol, err;
		self->kernel->volume(region, samples, &vol, &err);
		return Py_BuildValue("dd",vol, err);
//		return PyFloat_FromDouble(kernel.volume(region, samples, &err));
	}
	Py_RETURN_NONE;
} // Viewer_volume

/* --- Viewer_voxel --- */
static PyObject* Viewer_voxel(ViewerObject *self, PyObject *args)
{
	char *name;
	PyObject *value=NULL, *opt=NULL;

	if (!PyArg_ParseTuple(args, "s|OO", &name, &value, &opt)) return NULL;

	// Prepare Viewer voxel if needed
	self->viewer->voxel().allocate();

	if (!strcmp(name,"color")) {	// set color of voxel
		if (value != NULL) {
			if (opt != NULL) {
				self->viewer->voxel().color((int)PyInt_AsLong(value)-1,
						(dword)PyInt_AsUnsignedLongMask(opt));
			} else
				return PyInt_FromLong(self->viewer->voxel().color((int)PyInt_AsLong(value)-1));
		}
	} else
	if (!strcmp(name, "value")) {
		if (value != NULL) {
			if (opt==NULL) {
				PyErr_SetString(PyExc_SyntaxError, "Cannot get voxel color value");
				return NULL;
			} else {
				// set color from the palette using the value
				int pal = self->viewer->palette.def();
				self->viewer->voxel().color((int)PyInt_AsLong(value)-1,
						self->viewer->palette[pal].color(PyFloat_AsDouble(opt)));
			}
		}
	} else
	if (!strcmp(name, "roi")) {
		if (value != NULL) {
			if (PyInt_Check(value)) {
				int roi = PyInt_AsLong(value);
				self->viewer->voxel().roiShow(roi, PyInt_AsLong(opt));
			}
		}
	} else
	if (!strcmp(name, "roiclear")) {
		self->viewer->voxel().roiShowClear();
	} else
	if (!strcmp(name, "roialpha")) {
		if (value != NULL && PyInt_Check(value))
			self->viewer->voxel().roiAlpha(PyInt_AsLong(value));
		else
			return PyInt_FromLong(self->viewer->voxel().roiAlpha());
	} else {
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid option", name);
		return NULL;
	}

	Py_RETURN_NONE;
} // Viewer_voxel

/* --- Viewer_where --- */
static PyObject* Viewer_where(ViewerObject *self, PyObject *args)
{
	double x, y, z;

	if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) return NULL;

	Point pos(x,y,z);
	Vector dir(-self->kernel->view.matrix()(0,2),
		   -self->kernel->view.matrix()(1,2),
		   -self->kernel->view.matrix()(2,2));

	// FIXME MAX lattices
	VZone *zone = NULL;
	for (int i=0; i<10; i++) {
		self->kernel->lock();
		self->kernel->engine()->incBodyCheckId();
		zone = self->kernel->engine()->where(pos.x,pos.y,pos.z, dir.x,dir.y,dir.z, NULL);
		self->kernel->unlock();

		if (zone==NULL)
			Py_RETURN_NONE;

		//		if (zone->gregion()->rotdefi==0)
		if (!zone->gregion()->hasMatrix())
			return PyString_FromString(zone->gregion()->name());
		else {
			pos = zone->gregion()->matrix() * pos;
			dir = zone->gregion()->matrix() * dir;
		}

		//		pos = self->geometry->geometry->rotdefi(zone->gregion()->rotdefi) * pos;
//		dir = self->geometry->geometry->rotdefi(zone->gregion()->rotdefi)->multVector(dir);
//		dir = self->geometry->geometry->rotdefi(zone->gregion()->rotdefi) * dir;
	}
	return PyInt_FromLong(zone->gregion()->id());
} // Viewer_where

/* --- Viewer_xyz2uv --- */
static PyObject* Viewer_xyz2uv(ViewerObject *self, PyObject *args)
{
	double x,y,z;

	if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) return NULL;

	double u,v;
	// FIXME should be relative to the center of the screen!!!
	self->kernel->view.xyz2uv(x,y,z,&u,&v);

	return Py_BuildValue("dd", u,v);
} // Viewer_xyz2uv

/* --- Viewer_xyz2uvw --- */
static PyObject* Viewer_xyz2uvw(ViewerObject *self, PyObject *args)
{
	double x,y,z;

	if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) return NULL;

	double u,v,w;

	self->kernel->view.xyz2uvw(x,y,z,&u,&v,&w);

	return Py_BuildValue("ddd", u,v,w);
} // Viewer_xyz2uvw

/* --- Viewer_zoom --- */
static PyObject* Viewer_zoom(ViewerObject *self, PyObject *args)
{
	double zoom=-1.0;

	if (!PyArg_ParseTuple(args, "|d", &zoom)) return NULL;
	if (zoom<0.0)
		return PyFloat_FromDouble(self->kernel->view.zoom());
	else
		self->kernel->view.zoom(zoom);
	Py_RETURN_NONE;
} // Viewer_zoom

/** _selectZone
 * @return false on error
 */
static void _selectZone(ViewerObject *self, PyObject* obj)
{
	if (PyTuple_Check(obj) && PyTuple_GET_SIZE(obj)==2) {
		// Set zone
		bool rpn = false;
		PyObject *mode = PyTuple_GetItem(obj,0);
		if (PyString_Check(mode))
		  //			rpn = (bool)!strcmp(PyString_AsString(mode),"RPN");
		  rpn = (bool)!strcmp(PyUnicode_AsUTF8AndSize(mode, lllen),"RPN");
		else
		if (PyInt_Check(mode))
			rpn = (bool)PyInt_AsLong(mode);
		PyObject *zone = PyTuple_GetItem(obj,1);

		GZone* selzone = self->geometry->geometry->editRegion().addZone(rpn,
				(int)PyList_GET_SIZE(zone));
		for (ssize_t j=0; j<PyList_GET_SIZE(zone); j++) {
		  //			char *token = PyString_AsString(PyList_GetItem(zone, j));
		  //		  const char *token = PyUnicode_AsUTF8AndSize(PyList_GetItem(zone, j), lllen);
		  const char *token=PyBytes_AsString(PyList_GetItem(zone, j));
			if (PyErr_Occurred()) {
				PyErr_SetString(PyExc_ValueError,
					"Invalid region expression, string expected");
				return;
			}
			GBody *body = self->geometry->geometry->getBody(token);
			if (!selzone->add(token, body)) {
				PyErr_Format(PyExc_ValueError,"invalid body '%s'", token);
				return;
			}
		}
	} else
		PyErr_SetString(PyExc_ValueError, "Invalid zone expression");
} // _selectZone

/* --- Viewer_zone --- */
static PyObject* Viewer_zone(ViewerObject *self, PyObject *args)
{
	char *cmd;
	PyObject *value=NULL;

	if (!PyArg_ParseTuple(args, "s|O", &cmd, &value)) return NULL;

	if (!strcmp(cmd, "find")) {
		if (value==NULL) {
			PyErr_SetString(PyExc_SyntaxError, "tuple expected with (u,v,mask)");
			return NULL;
		}
		if (PyTuple_Check(value) && PyTuple_GET_SIZE(value)==3) {
			double u = Py_GetFloat(PyTuple_GetItem(value,0));
			double v = Py_GetFloat(PyTuple_GetItem(value,1));
			dword mask = (dword)Py_GetInt(PyTuple_GetItem(value,2));

			GZone zone;
			// Do I need a lockRead/Write?
			self->kernel->engine()->incBodyCheckId();

			const ViewPort& V = self->kernel->view;
			const Matrix4& M = V.matrix();

			// First pass add only the positive bodies
			ArrayIterator<GBody*> iter(self->geometry->geometry->bodies);
			while (iter) {
				GBody *body = iter++;
				if (!(body->show & mask)) continue;
				VBody *vbody = self->kernel->getBody(body);
				double x  = V.uv2x(u,v);
				double y  = V.uv2y(u,v);
				double z  = V.uv2z(u,v);
				if (!vbody->inside2D(x,y,z,-M(0,2),-M(1,2),-M(2,2))) continue;
				zone.add(body->name(), body, false);
			}
			zone.add("$", &GBody::tnull);
			iter.restart();
			// Second pass add the negative terms
			while (iter) {
				GBody *body = iter++;
				if (!(body->show & mask)) continue;
				VBody *vbody = self->kernel->getBody(body);
				double x  = V.uv2x(u,v);
				double y  = V.uv2y(u,v);
				double z  = V.uv2z(u,v);
				if (vbody->inside2D(x,y,z, -M(0,2),-M(1,2),-M(2,2))) continue;
				zone.add(body->name(), body, false);
			}

			// Optimize the zones expression
			zone.optimize();

			// empty? only $?
			if (zone.size()==1) Py_RETURN_NONE;
			return Py_ZoneExpr(&zone);
		} else {
			PyErr_SetString(PyExc_SyntaxError, "tuple expected with (u,v,mask)");
			return NULL;
		}
	} else
	// check if it is already selected
	if (!strcmp(cmd, "has") && PyTuple_Check(value) && PyTuple_GET_SIZE(value)==2) {
		// construct the hash as in GZone::()
		dword _hash = 0;
		PyObject *zonelist = PyTuple_GetItem(value,1);	// 0 is the type
		for (ssize_t i=0; i<PyList_GET_SIZE(zonelist); i++) {
			PyObject *obj = PyList_GetItem(zonelist,i);
			_hash += (_hash<<5) + i;
			//			_hash += (_hash<<5) + hash_djb2(PyString_AsString(obj));
			_hash += (_hash<<5) + hash_djb2(PyUnicode_AsUTF8(obj));
		}

		// check all zones in editRegion
		if (&self->geometry->geometry->editRegion()){
		  ArrayIterator<GZone*> iter(self->geometry->geometry->editRegion().zones());
		  while (iter) {
		    GZone* zone = iter++;
		    if (zone->hash() == _hash) {
		      Py_RETURN_TRUE;
		    }
		  }
		}
		Py_RETURN_FALSE;
	} else
	if (!strcmp(cmd, "select")) {
		if (value==NULL)
			return Py_RegionExpr(&self->geometry->geometry->editRegion());
		else
		if (PyList_Check(value)) {
			self->geometry->geometry->lockEdit();

			for (ssize_t i=0; i<PyList_GET_SIZE(value); i++) {
				PyObject *obj = PyList_GetItem(value,i);
				_selectZone(self, obj);
				if (PyErr_Occurred()) {
					self->geometry->geometry->editRegion().clear();
					self->geometry->geometry->unlockEdit();
					return NULL;
				}
			}
			self->geometry->geometry->unlockEdit();
		} else
		if (PyTuple_Check(value)) {
			self->geometry->geometry->lockEdit();
			_selectZone(self, value);
			self->geometry->geometry->unlockEdit();
			if (PyErr_Occurred()) return NULL;
		}
	} else
	if (!strcmp(cmd, "show")) {
		if (value==NULL)
			return PyInt_FromLong(self->geometry->geometry->editRegion().show);
		else
			self->geometry->geometry->editRegion().show = PyInt_AsLong(value);
	} else
	if (!strcmp(cmd, "clear")) {
		self->geometry->geometry->lockEdit();
		self->geometry->geometry->editRegion().clear();
		self->geometry->geometry->unlockEdit();
	} else {
		PyErr_Format(PyExc_TypeError, "Invalid type \'%s\' specified", cmd);
		return NULL;
	}

	Py_RETURN_NONE;
} // Viewer_zone

/* ----------------- Viewer type ------------------- */
static PyMethodDef Viewer_methods[] = {
	{"Sx",		(PyCFunction)Viewer_Sx,		METH_NOARGS,	"Return Sx"},
	{"Sy",		(PyCFunction)Viewer_Sy,		METH_NOARGS,	"Return Sy"},
	{"aspect",	(PyCFunction)Viewer_aspect,	METH_VARARGS,	"Return/Set aspect"},
	{"basis",	(PyCFunction)Viewer_basis,	METH_VARARGS,	"Return/Set basis"},
	{"bbox",	(PyCFunction)Viewer_bbox,	METH_VARARGS,	"Return the transformed(Viewer) 3D bounding box of body/zone/region in absolute coordinates"},
	{"bbox2D",	(PyCFunction)Viewer_bbox2D,	METH_VARARGS,	"Return the 2D bounding box of body/zone/region"},
	{"body",	(PyCFunction)Viewer_body,	METH_VARARGS,	"Modify body on viewer plane"},
	{"calcWindow",	(PyCFunction)Viewer_calcWindow,	METH_VARARGS,	"Calculate and set correct window limits"},
	{"camera",	(PyCFunction)Viewer_camera,	METH_VARARGS,	"Return 3D camera parameters"},
	{"closest",	(PyCFunction)Viewer_closest,	METH_VARARGS,	"Return closest objects"},
	{"configure",	(PyCFunction)Viewer_configure,	METH_VARARGS,	"Configure"},
	{"derive",	(PyCFunction)Viewer_derive,	METH_NOARGS,	"Derive viewer geometry from Geometry"},
	{"draw",	(PyCFunction)Viewer_draw,	METH_VARARGS,	"Start a drawing thread of viewer"},
	{"duv2dxyz",	(PyCFunction)Viewer_duv2dxyz,	METH_VARARGS,	"Convert d(relative) to d(absolute)"},
	{"edit",	(PyCFunction)Viewer_edit,	METH_VARARGS,	"Get/Set editing body"},
	{"enclosed",	(PyCFunction)Viewer_enclosed,	METH_VARARGS,	"Return objects enclosed inside region"},
	{"error",	(PyCFunction)Viewer_error,	METH_VARARGS,	"Return number of error"},
	{"export",	(PyCFunction)Viewer_export,	METH_VARARGS,	"Export viewport to vector file"},
	{"expose",	(PyCFunction)Viewer_expose,	METH_VARARGS,	"Expose image to X11"},
	{"extends",	(PyCFunction)Viewer_extends,	METH_VARARGS,	"Get/Set window/extend dimensions"},
	{"font",	(PyCFunction)Viewer_font,	METH_VARARGS,	"load font file"},
	{"fpe",		(PyCFunction)Viewer_fpe,	METH_NOARGS,	"Activate FPE exceptions"},
	{"get",		(PyCFunction)Viewer_set,	METH_VARARGS,	"Return flag/parameter"},
	{"grid",	(PyCFunction)Viewer_grid,	METH_VARARGS,	"Return/set grid settings"},
	{"hit3D",	(PyCFunction)Viewer_hit3D,	METH_VARARGS,	"Return hit position in 3D"},
	{"image",	(PyCFunction)Viewer_image,	METH_VARARGS,	"Get/Set image data"},
	{"invalid",	(PyCFunction)Viewer_invalid,	METH_NOARGS,	"Return true if window is invalid"},
	{"matrix",	(PyCFunction)Viewer_matrix,	METH_VARARGS,	"Return/set transformation matrix"},
	{"memory",	(PyCFunction)Viewer_memory,	METH_VARARGS,	"Get or print memory"},
	{"moveUV",	(PyCFunction)Viewer_moveUV,	METH_VARARGS,	"Move an object on uv plane by preserving w"},
	{"object",	(PyCFunction)Viewer_object,	METH_VARARGS,	"Modify object on viewer plane"},
	{"offset",	(PyCFunction)Viewer_offset,	METH_VARARGS,	"Return/Set projection offset (in UV)"},
	{"origin",	(PyCFunction)Viewer_origin,	METH_VARARGS,	"Return/Set projection origin"},
	{"palette",	(PyCFunction)Viewer_palette,	METH_VARARGS,	"Return/set palette"},
	{"pen",         (PyCFunction)Viewer_pen,	METH_VARARGS,	"Return/set pen drawing information"},
	{"pixel2uv",	(PyCFunction)Viewer_pixel2uv,	METH_VARARGS,	"Convert image pixel coordinates (i,j) to viewport (u,v)"},
	{"pixel2xyz",	(PyCFunction)Viewer_pixel2xyz,	METH_VARARGS,	"Convert image pixel coordinates (i,j) to real (x,y,z)"},
	{"project",	(PyCFunction)Viewer_project,	METH_VARARGS,	"Perform projection of viewer"},
	{"projectionChanged",(PyCFunction)Viewer_projectionChanged,METH_VARARGS,"Return/set projection change status"},
	{"rectangle",	(PyCFunction)Viewer_rectangle,	METH_VARARGS,	"Draw rectangle"},
	{"region",	(PyCFunction)Viewer_region,	METH_VARARGS,	"Set/Get region properties"},
	{"set",		(PyCFunction)Viewer_set,	METH_VARARGS,	"Return/set flag/parameter"},
	{"size",	(PyCFunction)Viewer_size,	METH_NOARGS,	"Return viewer size"},
	{"snap",	(PyCFunction)Viewer_snap,	METH_VARARGS,	"Return closest rounded position aligned to grid"},
	{"state",	(PyCFunction)Viewer_state,	METH_NOARGS,	"Return thread state"},
	{"stopThread",	(PyCFunction)Viewer_stopThread,	METH_NOARGS,	"Stop projection/drawing thread"},
	{"title",	(PyCFunction)Viewer_title,	METH_VARARGS,	"Return/Set viewport title"},
	{"userdump",	(PyCFunction)Viewer_userdump,	METH_VARARGS,	"Set/Get userdump information"},
	{"usrbin",	(PyCFunction)Viewer_usrbin,	METH_VARARGS,	"Set/Get usrbin information"},
	{"uv2pixel",	(PyCFunction)Viewer_uv2pixel,	METH_VARARGS,	"Convert relative to image coordinates"},
	{"uv2xyz",	(PyCFunction)Viewer_uv2xyz,	METH_VARARGS,	"Convert relative to absolute"},
	{"viewOto0",	(PyCFunction)Viewer_viewOto0,	METH_NOARGS,	"Move view origin to 0"},
	{"viewport",	(PyCFunction)Viewer_viewport,	METH_VARARGS,	"Set/Get other viewport matrices"},
	{"vertex",	(PyCFunction)Viewer_vertex,	METH_VARARGS,	"Return all vertices projected on an axis"},
	{"volume",	(PyCFunction)Viewer_volume,	METH_VARARGS,	"Calculate volume of region"},
	{"voxel",	(PyCFunction)Viewer_voxel,	METH_VARARGS,	"Set/Get voxel information"},
	{"where",	(PyCFunction)Viewer_where,	METH_VARARGS,	"Get zone at location"},
	{"xyz2uv",	(PyCFunction)Viewer_xyz2uv,	METH_VARARGS,	"Convert absolute to relative"},
	{"xyz2uvw",	(PyCFunction)Viewer_xyz2uvw,	METH_VARARGS,	"Convert absolute to relative"},
	{"zone",	(PyCFunction)Viewer_zone,	METH_VARARGS,	"Manipulate selected zone"},
	{"zoom",	(PyCFunction)Viewer_zoom,	METH_VARARGS,	"Return/Set zoom level"},
//	{"lock",	(PyCFunction)Viewer_lock,	METH_NOARGS,	"Lock mutex"},
//	{"unlock",	(PyCFunction)Viewer_unlock,	METH_NOARGS,	"UnLock mutex"},
#ifdef MEM
	{"_destroy",	(PyCFunction)Viewer_destroy,	METH_NOARGS,	"Deallocate"},
#endif
	{NULL, NULL, 0, NULL}	// Sentinel
};

PyTypeObject ViewerType = {
  //	PyObject_HEAD_INIT(&PyType_Type)
  //	0,				// ob_size
        PyVarObject_HEAD_INIT(&PyType_Type,0)
	"geoviewer.Viewer",		// tp_name
	sizeof(ViewerObject),		// tp_basicsize
	0,				// tp_itemsize
	(destructor)Viewer_dealloc,	// tp_dealloc
	0,				// tp_print
	0,				// tp_getattr
	0,				// tp_setattr
	0,				// tp_compare
	0,				// tp_repr
	0,				// tp_as_number
	0,				// tp_as_sequence
	0,				// tp_as_mapping
	0,				// tp_hash
	0,				// tp_call
	0,				// tp_str
	0,				// tp_getattro
	0,				// tp_setattro
	0,				// tp_as_buffer
	Py_TPFLAGS_DEFAULT,		// tp_flags
	"Geometry Viewer object",	// tp_doc
	0,				// tp_traverse
	0,				// tp_clear
	0,				// tp_richcompare
	0,				// tp_weaklistoffset
	0,				// tp_iter
	0,				// tp_iternext
	Viewer_methods,			// tp_methods
	0,//Viewer_members,		// tp_members
	0,				// tp_getset
	0,				// tp_base
	0,				// tp_dict
	0,				// tp_descr_get
	0,				// tp_descr_set
	0,				// tp_dictoffset
	(initproc)Viewer_init,		// tp_init
	0,				// tp_alloc
	Viewer_new,			// tp_new
	0,				// tp_free
	0,				// tp_is_gc
	0,				// tp_bases
	0,				// tp_mro
	0,				// tp_cache
	0,				// tp_subclasses
	0,				// tp_weaklist
	0,				// tp_del
	// From Python 2.5.x the PyTypeObject needs one more initializer
#if PY_VERSION_HEX >= 0x02050300
	0				// tp_version_tag
#endif
};
