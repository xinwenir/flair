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

#ifndef __GEOMETRYOBJECT_H
#define __GEOMETRYOBJECT_H

#include <structmember.h>
#include <Python.h>

#include "os.h"
#include "point.h"

class GBody;
class GZone;
class GRegion;
class GObject;
class GRuler;
class BinTree;
class Geometry;
class Material;
class GeometryViewer;
template <class T> class Array;

extern PyTypeObject GeometryType;

enum CursorType {
	CursorVector,
	CursorBody,
	CursorRegion,
	CursorZone,
	CursorObject
};

struct GeometryObject {
	PyObject_HEAD
	Geometry	*geometry;		// geometry class

	bool		 cursorShow;		// show cursor
	int		 cursorId;		// id of cursor object
	int		 cursorSize;		// size of rotate and axis locked move
	int		 cursorMoveSize;	// size of move handler (central circle)
	CursorType	 cursorType;		// type of object that set the cursor
	Point		 cursor;		// cursor position

	pthread_t	thread;			// job thread

	double		 trackballSize;		// percent of window
	double		 snapDistance;		// minimum distance in pixel to align
	double		 snapAngle;		// minimum angle in radians to align

	BinTree		*objects;		// drawable objects on screen
	Array<GObject*>	*objectList;		// drawable objects on screen

	int		 rulers;		// Rulers to show only during body editing
	GRuler*		 ruler[3];		// rulers for displaying bodies dimensions
};

/* function prototypes */
GBody*	  Py_GBody(   GeometryObject *self, PyObject *obj);
GRegion*  Py_GRegion( GeometryObject *self, PyObject *obj);
GObject*  Py_Object(  GeometryObject *self, PyObject *obj);
Material* Py_Material(GeometryObject *self, PyObject *obj);
PyObject* Py_ZoneExpr(GZone* zone);
PyObject* Py_RegionExpr(GRegion* region);

/* type */
extern PyTypeObject GeometryType;

#endif
