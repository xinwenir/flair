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
 */

#ifndef __VIEWEROBJECT_H
#define __VIEWEROBJECT_H

#include <Python.h>
#include <structmember.h>

#include <tk.h>
#include <tcl.h>
#include <X11/X.h>
#include <X11/Xlib.h>

#ifdef OPENGL
#	include <X11/Xutil.h>
#	include <GL/gl.h>
#	include <GL/glx.h>
//#	include <GL/glext.h>
#endif

#include "point.h"
#include "vector.h"
#include "painter.h"
#include "geometryobject.h"

class GeometryKernel;
class GeometryViewer;
class GObject;
template <class T> class Array;

#define NVIEWS		3	// Other views

struct ViewportLine {
	GeometryViewer	*viewer;		// viewer object
	bool		 visible;		// inside window
	double		 uc, vc;		// real center of projection
	int		 xc, yc;		// visible center of projection
	int		 xw, yw;		// direction of the w-vector from xc,yc
	int		 x[2], y[2];		// intersections with the window
	int		 lineWidth;		// line width
	int		 originWidth;		// cursor size
	int		 color;			// color of the line
	double		 g, f, c;		// conic representation of the line
	double		 cc;			// c normalized to the xc,yc center
	double		 wLength;		// length of w
};

extern PyTypeObject ViewerType;

struct ViewerObject {
	PyObject_HEAD

	Display		*display;		// X11 display
	Tk_Window	 tkWin;			// tk window handler
	Window		 window;		// X11 drawable window
	GC		 gc;			// Graphics context
	XImage		*ximage;		// XImage structure containing
	Pixmap		 pixmap;		// Pixmap used for all drawing operations
	int		 depth;			// display depth
	bool		 threaded;		// tcl is threaded or not

#ifdef OPENGL
	XVisualInfo	*visinfo;		// visual info for OpenGL
	GLXContext	 context;		// GL context
	GLXPixmap	 glxpixmap;		// GL pixmap
#endif

	GeometryViewer	*viewer;		// geometry viewer class
	GeometryKernel	*kernel;		// geometry kernel object
	GeometryObject	*geometry;		// geometry object python class

	Array<GObject*> *errors;		// drawable errors

	//XXX for the moment stay like this, it should be moved to GeometryObject
	ViewportLine	 viewport[NVIEWS];	// other viewport cross sections

	bool	 showTrackball;			// show trackball
	bool	 projectionChanged;		// true if projection has changed
	bool	 showViewport;			// show other viewport cross sections
	int	 crosshair;			// cross hair size at the center

	Point	 start;				// start absolute position
	Vector	 move;				// move vector
	double	 startU, startV;		// start point
	double	 pivotU, pivotV;		// pivot point
	double	 pivotAngle;			// pivot angle
	Vector	 rotateAxis;			// rotation axis
	double	 rotateAngle;			// angle to rotate
	double	 rotateCos, rotateSin;		// rotate Cos and Sin

	int	 rectX1, rectY1;		// Zoom rectangle
	int	 rectX2, rectY2;

	Array<IPoint>	 *pen;			// handdrawn points
};

/* type */
extern PyTypeObject ViewerType;
#endif
