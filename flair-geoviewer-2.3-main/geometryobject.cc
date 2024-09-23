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
 * Last change: 25-May-2022 by Paola Sala
 */

#include <Python.h>

#include <stdio.h>
#include <fnmatch.h>

#include "os.h"
#include "geo.h"
#include "gbody.h"
#include "gmesh.h"
#include "obbox.h"
#include "voxel.h"
#include "gobject.h"
#include "gregion.h"
#include "matrix4.h"
#include "pyutils.h"
#include "geometry.h"

using namespace std;

#include "geometryobject.h"	// it has to be the first include!

#ifdef MEM
// activate them only for memory leaks debugging
#include "array.h"
#include "bbox.h"
#include "bintree.h"
#include "bmath.h"
#include "bstring.h"
#include "cbody.h"
#include "cell_line.h"
#include "color.h"
#include "conic.h"
#include "csg.h"
#include "d2layer.h"
#include "d3layer.h"
#include "decorationlayer.h"
#include "dxfexport.h"
#include "edge.h"
#include "engine.h"
#include "eps.h"
#include "eventbin.h"
#include "exportbase.h"
#include "exportlayer.h"
#include "face.h"
#include "format.h"
#include "fortran.h"
#include "gbody.h"
#include "geo.h"
#include "geometry.h"
#include "geometryobject.h"
#include "glmesh.h"
#include "gmesh.h"
#include "gobject.h"
#include "gregion.h"
#include "gzone.h"
#include "histogram.h"
#include "imagelayer.h"
#include "kernel.h"
#include "latticelayer.h"
#include "layer.h"
#include "light.h"
#include "line.h"
#include "list.h"
#include "material.h"
#include "matrix.h"
#include "matrix2.h"
#include "matrix3.h"
#include "matrix4.h"
#include "memory.h"
#include "mesh.h"
#include "nox.h"
#include "obbox.h"
#include "optimizer.h"
#include "optimizerobject.h"
#include "os.h"
#include "painter.h"
#include "palette.h"
#include "palettelayer.h"
//#include "particle.h"
#include "pbmatrix.h"
#include "pencilbeam.h"
#include "point.h"
#include "polyline.h"
#include "polynomial.h"
#include "pyutils.h"
#include "quad.h"
#include "random.h"
#include "ray.h"
#include "roi.h"
#include "scatter.h"
#include "spline.h"
#include "stl.h"
#include "stream.h"
#include "svgexport.h"
#include "tetra.h"
#include "tetramesh.h"
#include "tetrameshio.h"
#include "threadpool.h"
#include "timer.h"
#include "token.h"
//#include "units.h"
#include "userdump.h"
#include "userdumplayer.h"
#include "usrbin.h"
#include "usrbinlayer.h"
#include "vbody.h"
#include "vector.h"
#include "vertex.h"
#include "vertex2d.h"
#include "viewer.h"
#include "viewerobject.h"
#include "viewport.h"
#include "voxel.h"
#include "voxellayer.h"
#include "vregion.h"
#include "vzone.h"
#include "xdraw.h"
#endif
//
Py_ssize_t *llen;
//
/* _deleteObjects */
static void _deleteObjects(GeometryObject *self)
{
	// delete all objects
	for (int i=0; i<self->objectList->size(); i++)
		delete (*self->objectList)[i];
	self->objectList->clear();
	self->objects->destroyTree();
} /* _deleteObjects */

/* Py_GBody
 * @param self	Geometry object
 * @param obj	object to convert to body can be string or id
 * @return body
 */
GBody* Py_GBody(GeometryObject *self, PyObject *obj)
{
	GBody *body;
	if (obj==NULL) return NULL;
	//	if (PyString_Check(obj) || PyUnicode_Check(obj)) {
	//		char *name = PyString_AsString(obj);		// XXX Unicode as String?
	if (PyUnicode_Check(obj)) {
	  const char *name = PyUnicode_AsUTF8AndSize(obj,llen);
	  body = self->geometry->getBody(name);
		if (body==NULL) {
			PyErr_Format(PyExc_KeyError, "Body \'%s\' not found", name);
			return NULL;
		}
		return body;
	} else
	if (PyInt_Check(obj)) {
		int id = (int)PyInt_AsLong(obj);
		if (id>=0 && id<self->geometry->bodies.count()) {
			body = self->geometry->getBody(id);
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
} /* Py_GBody */

/* Py_GRegion
 * @param self	Geometry object
 * @param obj	object to convert to region can be string or id
 * @return region
 */
GRegion* Py_GRegion(GeometryObject *self, PyObject *obj)
{
	GRegion *region;
	if (obj==NULL) return NULL;
	if (PyString_Check(obj) || PyUnicode_Check(obj)) {
	  //		char *name = PyString_AsString(obj);		// XXX Unicode as String?
	  const char *name = PyUnicode_AsUTF8AndSize(obj, llen);
		region = self->geometry->getRegion(name);
		if (region==NULL) {
			PyErr_Format(PyExc_KeyError, "Region \'%s\' not found", name);
			return NULL;
		}
		return region;
	} else
	if (PyInt_Check(obj)) {
		int id = (int)PyInt_AsLong(obj);
		region = self->geometry->getRegion(id);
		if (region==NULL) {
			PyErr_Format(PyExc_IndexError, "Region #%d not found",id);
			return NULL;
		}
		return region;
	} else {
		PyErr_SetString(PyExc_TypeError, "Invalid region type, string or integer expected");
		return NULL;
	}
} /* Py_GRegion */

/* Py_Material
 * @param self	Geometry object
 * @param obj	object to convert to material can be string or id
 * @return material
 */
Material* Py_Material(GeometryObject *self, PyObject *obj)
{
	Material *material;
	if (obj==NULL) return NULL;
	if (PyString_Check(obj) || PyUnicode_Check(obj)) {
	  //		char *name = PyString_AsString(obj);
	  const char *name = PyUnicode_AsUTF8AndSize(obj, llen);
		material = self->geometry->getMaterial(name);
		if (material==NULL) {
			PyErr_Format(PyExc_KeyError, "Material \'%s\' not found", name);
			return NULL;
		}
		return material;
	} else
	if (PyInt_Check(obj)) {
		int id = (int)PyInt_AsLong(obj);
		material = self->geometry->getMaterial(id);
		if (material==NULL) {
			PyErr_Format(PyExc_IndexError, "Material #%d not found",id);
			return NULL;
		}
		return material;
	} else {
		PyErr_SetString(PyExc_TypeError, "Invalid material type, string or integer expected");
		return NULL;
	}
} /* Py_Material */

/* Py_Object
 * @param self	Geometry object
 * @param obj	object to convert to body can be string or id
 * @return object
 */
GObject* Py_Object(GeometryObject *self, PyObject *obj)
{
	GObject *object=NULL;
	if (obj==NULL) return NULL;
	if (PyString_Check(obj) || PyUnicode_Check(obj)) {
	  //		char *name = PyString_AsString(obj);
	  const char *name = PyUnicode_AsUTF8AndSize(obj, llen);
	  object = (GObject*)self->objects->find(name);
		if (object==NULL) {
			PyErr_Format(PyExc_KeyError, "Object \'%s\' not found", name);
			return NULL;
		}
		return object;
	} else
	if (PyInt_Check(obj)) {
		int id = (int)PyInt_AsLong(obj);
		if (id>=0 && id<self->objectList->count()) {
			object = self->objectList->get(id);
			if (object==NULL) {
				PyErr_Format(PyExc_IndexError, "Object #%d not found",id);
				return NULL;
			}
			return object;
		}
		PyErr_Format(PyExc_IndexError, "Object #%d not found",id);
		return NULL;
	} else {
		PyErr_SetString(PyExc_TypeError, "Invalid object type, string or integer expected");
		return NULL;
	}
	return NULL;
} /* Py_Object */

/** Py_ZoneExpr */
PyObject* Py_ZoneExpr(GZone* zone)
{
	if (!zone->size()) return NULL;

	PyObject *expr = PyList_New(zone->size());
	for (ssize_t j=0; j<zone->size(); j++) {
		const GBody *body = (*zone)[j];
		if (!body)
			PyList_SetItem(expr, j, PyString_FromString("$"));
		else
			PyList_SetItem(expr, j, PyString_FromString(body->name()));
	}
	return Py_BuildValue("iO",(int)zone->rpn(),expr);
} /* Py_ZoneExpr */

/** Py_RegionExpr */
PyObject* Py_RegionExpr(GRegion* region)
{
	// return selected zone
	if (!region->nzones()) Py_RETURN_NONE;

	PyObject *zones = PyList_New(0);
	ArrayIterator<GZone*> iter(region->zones());
	while (iter) {
		GZone* zone = iter++;
		PyObject* z = Py_ZoneExpr(zone);
		if (z) PyList_Append(zones, z);
	}
	return zones;
} /* Py_RegionExpr */

//////////////////////// Geometry ///////////////////////////
/* --- Geometry_new --- */
static PyObject *Geometry_new(PyTypeObject *type, PyObject* /*args*/, PyObject* /*kwds*/)
{
	GeometryObject *self;
	self = (GeometryObject*)type->tp_alloc(type, 0);
	if (self == NULL) return NULL;

	self->geometry       = new Geometry();
	self->cursor         = Vector::O;
	self->cursorType     = CursorVector;
	self->cursorId       = -1;
	self->cursorSize     = 40;
	self->cursorMoveSize = 10;
	self->trackballSize  = 0.80;
	self->snapDistance   = 4.0;
	self->snapAngle     = RAD(5.0);
	self->objects        = new BinTree();
	self->objectList     = new Array<GObject*>(16);

	self->thread  = 0;
	self->rulers  = 3;

	self->ruler[0] = new GRuler("X");
	self->ruler[0]->color  = 0xFF00FF; //COLOR_X_DARK;
	self->ruler[0]->anchor = Anchor_C;

	self->ruler[1] = new GRuler("Y");
	self->ruler[1]->color  = 0xFF00FF; //COLOR_Y_DARK;
	self->ruler[1]->anchor = Anchor_C;

	self->ruler[2] = new GRuler("Z");
	self->ruler[2]->color  = 0xFF00FF; //COLOR_Z_DARK;
	self->ruler[2]->anchor = Anchor_C;

	return (PyObject*)self;
} /* Geometry_new */

/* --- Geometry_dealloc --- */
static void Geometry_dealloc(GeometryObject *self)
{
	for (int i=0; i<3; i++)
		delete self->ruler[i];
	_deleteObjects(self);
	delete self->objectList;
	delete self->objects;
	delete self->geometry;
} /* Geometry_dealloc */

#ifdef MEM
/* --- Geometry_destroy --- */
static PyObject* Geometry_destroy(GeometryObject*self)
{
#define PRINT_SIZE(a)	printf("  %-16s %5ld\n",#a,sizeof(a));
	printf("class size for reference:\n");

	PRINT_SIZE(ArrayEmpty);
	PRINT_SIZE(ArrayEmpty);
	PRINT_SIZE(BBox);
	PRINT_SIZE(BFont);
	PRINT_SIZE(BaseSplineNode);
	PRINT_SIZE(BinLeaf);
	PRINT_SIZE(BinTree);
	PRINT_SIZE(Body3DFeeder);
	PRINT_SIZE(Body3DWorker);
	PRINT_SIZE(BodyFeeder);
	PRINT_SIZE(BodyWorker);
	PRINT_SIZE(CBody);
	PRINT_SIZE(CBodyOrderAccel);
	PRINT_SIZE(CSG);
//	PRINT_SIZE(Caloron);
	PRINT_SIZE(CardinalSpline);
	PRINT_SIZE(Cell_Line);
	PRINT_SIZE(Color);
	PRINT_SIZE(Color3D);
	PRINT_SIZE(Conic);
	PRINT_SIZE(D2Layer);
	PRINT_SIZE(D3Layer);
	PRINT_SIZE(DXFExport);
	PRINT_SIZE(DecorationLayer);
	PRINT_SIZE(Edge);
	PRINT_SIZE(Eventbin);
	PRINT_SIZE(ExportBase);
	PRINT_SIZE(ExportLayer);
	PRINT_SIZE(Face);
	PRINT_SIZE(FortranFile);
	PRINT_SIZE(FortranParser);
	PRINT_SIZE(GARBBody);
	PRINT_SIZE(GArrow);
	PRINT_SIZE(GBOXBody);
	PRINT_SIZE(GBeam);
	PRINT_SIZE(GBody);
	PRINT_SIZE(GELLBody);
	PRINT_SIZE(GERRBody);
	PRINT_SIZE(GInfEllCylBody);
	PRINT_SIZE(GLMesh);
	PRINT_SIZE(GLight);
	PRINT_SIZE(GMesh);
	PRINT_SIZE(GOPRBody);
	PRINT_SIZE(GObject);
	PRINT_SIZE(GPLABody);
	PRINT_SIZE(GPoint);
	PRINT_SIZE(GQUABody);
	PRINT_SIZE(GRCCBody);
	PRINT_SIZE(GRECBody);
	PRINT_SIZE(GRegion);
	PRINT_SIZE(GRotdefi);
	PRINT_SIZE(GRuler);
	PRINT_SIZE(GSPHBody);
	PRINT_SIZE(GSpline);
	PRINT_SIZE(GTRCBody);
	PRINT_SIZE(GTorusBody);
	PRINT_SIZE(GVoxel);
	PRINT_SIZE(GWEDBody);
	PRINT_SIZE(GZone);
	PRINT_SIZE(Geometry);
	PRINT_SIZE(GeometryEngine);
	PRINT_SIZE(GeometryKernel);
	PRINT_SIZE(GeometryViewer);
	PRINT_SIZE(H1D);
	PRINT_SIZE(Histogram);
	PRINT_SIZE(IPoint);
	PRINT_SIZE(ImageLayer);
	PRINT_SIZE(KahanSum);
	PRINT_SIZE(LatticeLayer);
	PRINT_SIZE(Layer);
	PRINT_SIZE(Material);
	PRINT_SIZE(Matrix);
	PRINT_SIZE(Matrix2);
	PRINT_SIZE(Matrix3);
	PRINT_SIZE(Matrix4);
	PRINT_SIZE(Mesh);
//	PRINT_SIZE(Neutron);
	PRINT_SIZE(OBBox);
	PRINT_SIZE(Optimizer);
	PRINT_SIZE(PBMatrix);
	PRINT_SIZE(Painter);
	PRINT_SIZE(Palette);
	PRINT_SIZE(PaletteLayer);
//	PRINT_SIZE(Particle);
	PRINT_SIZE(PencilBeam);
	PRINT_SIZE(Plane);
	PRINT_SIZE(Point);
	PRINT_SIZE(Point2D);
	PRINT_SIZE(Polyline);
	PRINT_SIZE(Polynomial);
	PRINT_SIZE(PolynomialSolver);
	PRINT_SIZE(Quad);
	PRINT_SIZE(ROICombination);
	PRINT_SIZE(ROIPlanar);
	PRINT_SIZE(ROIPlanarSlice);
	PRINT_SIZE(Random);
	PRINT_SIZE(Ray);
	PRINT_SIZE(RaySegment);
	PRINT_SIZE(Roi);
	PRINT_SIZE(STL);
	PRINT_SIZE(SVGExport);
	PRINT_SIZE(Scatter);
	PRINT_SIZE(SplineNode);
	PRINT_SIZE(Stream);
	PRINT_SIZE(Tetra);
	PRINT_SIZE(TetraMesh);
	PRINT_SIZE(TetraMeshIO);
	PRINT_SIZE(ThreadPool);
	PRINT_SIZE(ThreadPoolFeeder);
	PRINT_SIZE(ThreadPoolWorker);
	PRINT_SIZE(Timer);
	PRINT_SIZE(Token);
	PRINT_SIZE(UserDump);
	PRINT_SIZE(UserDumpLayer);
	PRINT_SIZE(Usrbin);
	PRINT_SIZE(UsrbinLayer);
	PRINT_SIZE(VBody);
	PRINT_SIZE(VRegion);
	PRINT_SIZE(VVoxel);
	PRINT_SIZE(VZone);
	PRINT_SIZE(Vector);
	PRINT_SIZE(Vector2D);
	PRINT_SIZE(VectorSplineNode);
	PRINT_SIZE(Vertex);
	PRINT_SIZE(Vertex2D);
	PRINT_SIZE(ViewPort);
	PRINT_SIZE(VolumeFeeder);
	PRINT_SIZE(VolumeWorker);
	PRINT_SIZE(VoxelLayer);
	PRINT_SIZE(XYPair);
	PRINT_SIZE(ZoneOfPoint);
// struct
	printf("\nstruct size for reference:\n");
	PRINT_SIZE(GeometryObject);
	PRINT_SIZE(GeometryViewer);
	PRINT_SIZE(ViewportLine);

	PRINT_SIZE(UserDumpTrackPos);
	PRINT_SIZE(UserDumpSourceParticle);
	PRINT_SIZE(Light);
	PRINT_SIZE(VLight);
	PRINT_SIZE(ClipRegion);
	PRINT_SIZE(drand48_data);
//	PRINT_SIZE(Memory);
	PRINT_SIZE(_ThreadData);
	PRINT_SIZE(Color32);

//	optimizer/cell_line.h:typedef struct {
//	optimizer/pbmatrix.h:typedef struct {
//	optimizer/pbmatrix.h:typedef struct {
//	optimizer/pencilbeam.h:typedef struct {
//	optimizer/roi.h:typedef struct {
//	optimizer/roi.h:typedef struct {

	fflush(stdout);
	Memory::dump();
	Geometry_dealloc(self);

	if (!Memory::fini())
		printf("*** No memory leak\n");
	exit(0);
	Py_RETURN_NONE;
} /* Geometry_destroy */
#endif

/* --- Geometry_addBody --- */
static PyObject *Geometry_addBody(GeometryObject *self, PyObject *args)
{
	char *name, *type;

	if (!PyArg_ParseTuple(args, "ss", &name, &type)) return NULL;

	// Create body
	GBody *body = self->geometry->addBody(name, type);
	if (body)
		return PyInt_FromLong(body->id());
	else {
		// setup exception
		PyErr_SetString(PyExc_SyntaxError, "invalid body");
		return NULL;
	}
} /* Geometry_addBody */

/* --- Geometry_addMaterial --- */
static PyObject *Geometry_addMaterial(GeometryObject *self, PyObject *args)
{
	char *name;
	Material *material;

	if (!PyArg_ParseTuple(args, "s", &name)) return NULL;
	material = self->geometry->addMaterial(name);
	return PyInt_FromLong(material->id());
} /* Geometry_addMaterial */

/* --- Geometry_addObject --- */
static PyObject *Geometry_addObject(GeometryObject *self, PyObject *args)
{
	char *name, *type;
	GObject *object=NULL;

	if (!PyArg_ParseTuple(args, "ss", &name, &type)) return NULL;

	if (!strcmp(type,"point") || !strcmp(type,"!point"))
		object = new GPoint(name);
	else
	if (!strcmp(type,"arrow") || !strcmp(type,"!arrow"))
		object = new GArrow(name);
	else
	if (!strcmp(type,"spline") || !strcmp(type,"!spline"))
		object = new GSpline(name);
	else
	if (!strcmp(type,"mesh") || !strcmp(type,"!mesh"))
		object = new GMesh(name);
	else
	if (!strcmp(type,"ruler") || !strcmp(type,"!ruler"))
		object = new GRuler(name);
	else
	if (!strcmp(type,"rotdefi") || !strcmp(type,"!rotdefi"))
		object = new GRotdefi(name);
	else
	if (!strcmp(type,"light") || !strcmp(type,"!light"))
		object = new GLight(name);
	else
	if (!strcmp(type,"beam") || !strcmp(type,"!beam"))
		object = new GBeam(name);

	if (object) {
		object->id(self->objectList->count());
		self->objects->add(name, object);
		self->objectList->add(object);
		return PyInt_FromLong(object->id());
	}

	Py_RETURN_NONE;
} /* Geometry_addObject */

/* --- Geometry_addRegion --- */
static PyObject *Geometry_addRegion(GeometryObject *self, PyObject *args)
{
	char *name;
	GRegion *region;

	if (!PyArg_ParseTuple(args, "s", &name)) return NULL;
	region = self->geometry->addRegion(name);
	return PyInt_FromLong(region->id());
} /* Geometry_addRegion */

/* --- _region_addZone ---
 * add a zone to existing region
 * Zone should be in the format (["RPN"]|"STD"], [zonelist...])
 * @return true on failure, false otherwise
 */
static bool _region_addZone(GeometryObject *self, GRegion *region, PyObject *zoneobj)
{
	if (!PyTuple_Check(zoneobj)) {
		PyErr_SetString(PyExc_TypeError, "Invalid region expression, tuple expected");
		return true;
	}

	PyObject *mode = PyTuple_GetItem(zoneobj, 0);
	PyObject *zone = PyTuple_GetItem(zoneobj, 1);

	if (!PyString_Check(mode) || !PyList_Check(zone)) {
		PyErr_SetString(PyExc_TypeError, "Invalid region expression, tuple (mode,zone) expected");
		return true;
	}

	//	bool rpn = (bool)!strcmp(PyString_AsString(mode),"RPN");
        bool rpn = (bool)!strcmp(PyUnicode_AsUTF8AndSize(mode, llen),"RPN");
	self->geometry->addZone(region, rpn, (int)PyList_GET_SIZE(zone));
	for (ssize_t j=0; j<PyList_GET_SIZE(zone); j++) {
		char errmsg[128];
		PyObject *token = PyList_GetItem(zone, j);
		if (!PyString_Check(token))
			PyErr_SetString(PyExc_TypeError, "Invalid region expression, string expected");
		else
		if (self->geometry->add2exp(region, PyUnicode_AsUTF8AndSize(token, llen), errmsg))
			PyErr_SetString(PyExc_SyntaxError, errmsg);
	}
	if (PyErr_Occurred()) return true;
	return false;
} /* _region_addZone */

/* --- Geometry_bbox --- */
static PyObject* Geometry_bbox(GeometryObject *self, PyObject *args)
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
		GBody *body = Py_GBody(self,obj);
		if (body==NULL) return NULL;
		bbox = body->bbox();
	} else
	if (type[0]=='r' || type[0]=='R') {
		GRegion *region = Py_GRegion(self,obj);
		if (region==NULL) return NULL;
		bbox = region->bbox();
	} else
//	if (type[0]=='z' || type[0]=='Z') {
//		// FIXME not implemented yet
//	} else
	if (type[0]=='o' || type[0]=='O') {
		GObject *object = Py_Object(self,obj);
		bbox = object->bbox();
	}

	if (bbox.isValid())
		return Py_BuildValue("[dddddd]",
				bbox.low().x,  bbox.low().y,  bbox.low().z,
				bbox.high().x, bbox.high().y, bbox.high().z);
	else
		Py_RETURN_NONE;
} /* Geometry_bbox */

/* --- Geometry_bodyVar --- */
static PyObject *Geometry_bodyVar(GeometryObject *self, GBody *body, const char *var, PyObject *value)
{
	if (!strcmp(var, "name")) {
		if (value==NULL)
			return PyString_FromString(body->name());
		else {
			PyErr_SetString(PyExc_SyntaxError, "Cannot set body name");
			return NULL;
		}
	} else
	if (!strcmp(var, "id")) {
		if (value==NULL)
			return PyInt_FromLong(body->id());
		else {
			PyErr_SetString(PyExc_SyntaxError, "Cannot set body id");
			return NULL;
		}
	} else
	if (!strcmp(var, "show")) {
		if (value==NULL)
			return PyInt_FromLong(body->show);
		else
			body->show = (int)PyInt_AsLong(value);
	} else
	if (!strcmp(var, "type")) {
		if (value==NULL)
			return PyString_FromString(body->typeStr());
		else {
			PyErr_Format(PyExc_TypeError, "\'%s\' cannot be set", var);
			return NULL;
		}
	} else
	if (!strcmp(var, "color")) {
		if (value==NULL)
			return PyInt_FromLong(body->color());
		else
			body->color((dword)PyInt_AsUnsignedLongMask(value));
	} else
	if (!strcmp(var, "linewidth")) {
		if (value==NULL)
			return PyInt_FromLong(body->lineWidth());
		else
			body->lineWidth((dword)PyInt_AsLong(value));
	} else
	if (!strcmp(var, "what")) {
		double what[30];
		if (value==NULL) {
			int n = body->getWhat(what);
			assert(n>=0);
			PyObject *obj = PyList_New(n);
			for (int i=0; i<n; i++)
				PyList_SetItem(obj, i, PyFloat_FromDouble(what[i]));
			return obj;
		} else {
			if (Py_TYPE(value) != &PyList_Type) {
				PyErr_SetString(PyExc_TypeError, "Invalid type, list expected");
				return NULL;
			}

			if (PyList_GET_SIZE(value)>30) {
				PyErr_SetString(PyExc_IndexError, "Too large what list");
				return NULL;
			}
			ssize_t i;
			for (i=0; i<PyList_GET_SIZE(value); i++)
				what[i] = PyFloat_AS_DOUBLE(PyList_GetItem(value, i));
			for (;i<30; i++)
				what[i] = 0.0;

			char err[128];
			err[0] = 0;
			try {
				body->setWhat(what, err);
				body->create();
				if (body->hasMatrix()) body->transform();
			} catch ( ... ) {
				PyErr_SetString(PyExc_SyntaxError, "Error in body parameters");
				return NULL;
			}
			if (err[0])
				return PyString_FromString(err);
		}
	} else
	if (!strcmp(var, "matrix")) {
		if (value) {
			if (PyList_Check(value)) {
				Matrix4	matrix;
				PyList_AsMatrix4(value, matrix);
				try {
					body->matrix(matrix);
					body->transform();
				} catch ( ... ) {
					PyErr_SetString(PyExc_SyntaxError, "Invalid body rotation matrix");
					return NULL;
				}
			} else {
				body->clearMatrix();
				body->transform();
			}
		} //else return matrix
	} else
	if (!strcmp(var, "bbox")) {	// Get or set user defined bbox
		if (value==NULL) {
			BBox bb = body->bbox();
			Point lw = bb.low();
			Point hg = bb.high();

			PyObject *bbox  = PyList_New(6);
			PyList_SetItem(bbox, 0, PyFloat_FromDouble(lw.x));
			PyList_SetItem(bbox, 1, PyFloat_FromDouble(hg.x));
			PyList_SetItem(bbox, 2, PyFloat_FromDouble(lw.y));
			PyList_SetItem(bbox, 3, PyFloat_FromDouble(hg.y));
			PyList_SetItem(bbox, 4, PyFloat_FromDouble(lw.z));
			PyList_SetItem(bbox, 5, PyFloat_FromDouble(hg.z));

			return bbox;
		} else {
			if (PyList_Check(value)) {
				double xmin = PyFloat_AsDouble(PyList_GetItem(value,0));
				double xmax = PyFloat_AsDouble(PyList_GetItem(value,1));
				double ymin = PyFloat_AsDouble(PyList_GetItem(value,2));
				double ymax = PyFloat_AsDouble(PyList_GetItem(value,3));
				double zmin = PyFloat_AsDouble(PyList_GetItem(value,4));
				double zmax = PyFloat_AsDouble(PyList_GetItem(value,5));
				body->bbox(BBox(xmin,xmax, ymin,ymax, zmin,zmax));
			}
		}
	} else
	if (!strcmp(var, "pos")) {	// Move absolute
		if (value==NULL)
			return Py_Vector(body->position());
		else {
			Point r = Py_GetVector(value);
			if (PyErr_Occurred()) return NULL;
			body->position(r);
			body->create();
			if (body->hasMatrix()) body->transform();
			self->geometry->invalidateBody(body);
		}
	} else
	if (!strcmp(var, "save"))
		body->save();
	else
	if (!strcmp(var, "restore")) {
		body->restore();
		body->create();
		if (body->hasMatrix()) body->transform();
		self->geometry->invalidateBody(body);
	} else
	if (!strcmp(var, "savedpos"))	// return saved position
		return Py_Vector(body->savedPosition());
	else {
		PyErr_Format(PyExc_TypeError, "Invalid type \'%s\' specified", var);
		return NULL;
	}
	Py_RETURN_NONE;
} /* Geometry_bodyVar */

/* --- Geometry_body --- */
static PyObject *Geometry_body(GeometryObject *self, PyObject *args)
{
	PyObject *obj;
	char *var;
	GBody *body = NULL;
	PyObject *value=NULL;
	if (!PyArg_ParseTuple(args, "Os|O", &obj, &var, &value)) return NULL;

	if (Py_Check4Pattern(obj)) {
	  //		char *pattern = PyString_AsString(obj);
               const char *pattern = PyUnicode_AsUTF8AndSize(obj, llen);		ArrayIterator<GBody*> iter(self->geometry->bodies);
		while (iter) {
			body = iter++;
			if (!fnmatch(pattern, body->name(), 0)) {
				PyObject *ret = Geometry_bodyVar(self, body, var, value);
				Py_XDECREF(ret);	// ret can be NULL use XDECREF
			}
		}
		Py_RETURN_NONE;
	} else
	if (PyList_Check(obj)) {
		for (ssize_t i=0; i<PyList_GET_SIZE(obj); i++) {
			body = Py_GBody(self, PyList_GetItem(obj,i));
			if (body==NULL) return NULL;
			PyObject *ret = Geometry_bodyVar(self, body, var, value);
			Py_XDECREF(ret);	// ret can be NULL use XDECREF
		}
	} else
	if (PyTuple_Check(obj)) {
		for (ssize_t i=0; i<PyTuple_GET_SIZE(obj); i++) {
			body = Py_GBody(self, PyTuple_GetItem(obj,i));
			if (body==NULL) return NULL;
			PyObject *ret = Geometry_bodyVar(self, body, var, value);
			Py_XDECREF(ret);	// ret can be NULL use XDECREF
		}
	} else  {
		body = Py_GBody(self, obj);
		if (body==NULL) return NULL;
		return Geometry_bodyVar(self, body, var, value);
	}
	Py_RETURN_NONE;		// To keep compiler happy
} /* Geometry_body */

/* --- Geometry_cleanup --- */
static PyObject *Geometry_cleanup(GeometryObject *self)
{
	_deleteObjects(self);	// delete all objects

	// clean up geometry
	//self->geometry->voxel.cleanup();	Don't clean up heavy objects, needs time
	self->geometry->cleanup();
	Py_RETURN_NONE;
} /* Geometry_cleanup */

/* --- Geometry_cursor --- */
static PyObject *Geometry_cursor(GeometryObject *self, PyObject *args)
{
	char *name;
	PyObject *value=NULL;

	if (!PyArg_ParseTuple(args, "s|O", &name, &value)) return NULL;
	if (!strcmp(name, "body")) {	// Move relative from the last save position
		if (value) {
			if (value==Py_None) {
				self->geometry->editBody(NULL);
				self->cursorShow = false;
			} else {
				self->geometry->editBody(Py_GBody(self, value));
				if (self->geometry->editBody()==NULL) {
					self->cursorShow = false;
					return NULL;
				}
				self->cursorShow = true;
				self->cursorType = CursorBody;
				self->cursor     = self->geometry->editBody()->position();
				self->cursorId   = self->geometry->editBody()->id();
			}
		} else {
			if (self->geometry->editBody())
				return PyInt_FromLong(self->geometry->editBody()->id());
		}
	} else
	if (!strcmp(name, "region")) {	// Move relative from the last save position
		if (value) {
			if (value==Py_None)
				self->cursorShow = false;
			else {
				GRegion* region  = Py_GRegion(self, value);
				if (region==NULL) {
					self->cursorShow = false;
					return NULL;
				}
				self->cursorShow = true;
				self->cursorType = CursorRegion;
				self->cursor     = region->obbox()->center();
				self->cursorId   = region->id();
#ifdef _DUMP
				cout << "Region:" << region->name() << endl;
				cout << "\tObbox:" << *region->obbox() << endl;
				cout << "\tCursor:" << self->cursor << endl;
#endif
			}
		}
	} else
	if (!strcmp(name, "cursor")) {
		if (value) {
			if (value==Py_None || value==Py_False)
				self->cursorShow = false;
			else
			if (value==Py_True)
				self->cursorShow = true;
			else {
				self->cursor     = Py_GetVector(value);
				self->cursorType = CursorVector;
				if (PyErr_Occurred()) {
					self->cursorShow = false;
					return NULL;
				}
				self->cursorShow = true;
				self->cursorId   = -1;
			}
		} else
			if (self->cursorShow)
				return Py_Vector(self->cursor);
	} else
	if (!strcmp(name, "type")) {
		switch (self->cursorType) {
			case CursorVector:
				return PyString_FromString("vector");
			case CursorBody:
				return PyString_FromString("body");
			case CursorRegion:
				return PyString_FromString("region");
			case CursorZone:
				return PyString_FromString("zone");
			case CursorObject:
				return PyString_FromString("object");
		}
	} else
	if (!strcmp(name, "id")) {
		return PyInt_FromLong(self->cursorId);
	} else
	if (!strcmp(name, "clear")) {
		self->cursorShow = false;
		self->geometry->editBody(NULL);
	} else {
		PyErr_Format(PyExc_TypeError, "Invalid type \'%s\' specified", name);
		return NULL;
	}
	Py_RETURN_NONE;
} /* Geometry_cursor */

/* --- Geometry_materialVar --- */
static PyObject *Geometry_materialVar(GeometryObject* /*self*/, Material *material,
		const char *var, PyObject *value)
{
	if (!strcmp(var, "name")) {
		if (value==NULL)
			return PyString_FromString(material->name());
		else {
			PyErr_SetString(PyExc_SyntaxError, "Cannot set material name");
			return NULL;
		}
	} else
	if (!strcmp(var, "id")) {
		if (value==NULL)
			return PyInt_FromLong(material->id());
		else {
			PyErr_SetString(PyExc_SyntaxError, "Cannot set material id");
			return NULL;
		}
	} else
	if (!strcmp(var, "density")) {
		if (value==NULL)
			return PyFloat_FromDouble(material->density());
		else
			material->density(PyFloat_AsDouble(value));
	} else
	if (!strcmp(var, "Z")) {
		if (value==NULL)
			return PyInt_FromLong(material->Z());
		else
			material->Z(PyInt_AsLong(value));
	} else
	if (!strcmp(var, "A")) {
		if (value==NULL)
			return PyInt_FromLong(material->A());
		else
			material->A(PyInt_AsLong(value));
	} else
	if (!strcmp(var, "weight")) {
		if (value==NULL)
			return PyFloat_FromDouble(material->weight());
		else
			material->weight(PyFloat_AsDouble(value));
	} else
	if (!strcmp(var, "color")) {
		if (value==NULL)
			return PyInt_FromLong((int)material->color());
		else
			material->diffuse(PyInt_AsLong(value));
	} else
	if (!strcmp(var, "diffuse")) {
		if (value==NULL)
			return PyInt_FromLong((int)material->diffuse());
		else
			material->diffuse(PyInt_AsLong(value));
	} else
	if (!strcmp(var, "specular")) {
		if (value==NULL)
			return PyFloat_FromDouble(material->specular());
		else
			material->specular(PyFloat_AsDouble(value));
	} else
//	if (!strcmp(var, "transparent")) {
//		if (value==NULL)
//			return PyInt_FromLong((int)material->transparent());
//		else
//			material->transparent(PyInt_AsLong(value));
//	} else
	if (!strcmp(var, "shine")) {
		if (value==NULL)
			return PyFloat_FromDouble(material->shine());
		else
			material->shine(PyFloat_AsDouble(value));
	} else
//	if (!strcmp(var, "shineexp")) {
//		if (value==NULL)
//			return PyFloat_FromDouble(material->shineExp());
//		else
//			material->shineExp(PyFloat_AsDouble(value));
//	} else
	if (!strcmp(var, "ior")) {
		if (value==NULL)
			return PyFloat_FromDouble(material->ior());
		else
			material->ior(PyFloat_AsDouble(value));
	} else
	if (!strcmp(var, "fuzz")) {
		if (value==NULL)
			return PyFloat_FromDouble(material->fuzz());
		else
			material->fuzz(PyFloat_AsDouble(value));
	} else {
		PyErr_Format(PyExc_TypeError, "Invalid type \'%s\' specified", var);
		return NULL;
	}
	Py_RETURN_NONE;
} /* Geometry_materialVar */

/* --- Geometry_material --- */
static PyObject *Geometry_material(GeometryObject *self, PyObject *args)
{
	PyObject *obj;
	Material *material = NULL;
	char *var;
	PyObject *value=NULL;
	if (!PyArg_ParseTuple(args, "Os|O", &obj, &var, &value)) return NULL;

	if (Py_Check4Pattern(obj)) {
	  //		char *pattern = PyString_AsString(obj);
	  const char *pattern = PyUnicode_AsUTF8AndSize(obj, llen);
	  ArrayIterator<Material*> iter(self->geometry->materials);
		while (iter) {
			material = iter++;
			if (!fnmatch(pattern, material->name(), 0)) {
				PyObject *ret = Geometry_materialVar(self, material, var, value);
				Py_XDECREF(ret);	// ret can be NULL use XDECREF
			}
		}
		Py_RETURN_NONE;
	} else  {
		material = Py_Material(self, obj);
		if (material==NULL) return NULL;
		return Geometry_materialVar(self, material, var, value);
	}
	Py_RETURN_NONE;		// To keep compiler happy
} /* Geometry_material */

/* --- Geometry_memory --- */
static PyObject *Geometry_memory(GeometryObject *self, PyObject *args)
{
	char *dump=NULL;
	if (!PyArg_ParseTuple(args, "|s", &dump)) return NULL;

#ifdef MEM
	cout << "resident="  << Memory::resident() << endl;
	cout << "items="     << Memory::items() << endl;
	cout << "allocated=" << Memory::allocated() << endl;
#endif

	if (dump==NULL)
		return PyInt_FromLong(self->geometry->memory());
	else
		self->geometry->printMemory();

	Py_RETURN_NONE;
} /* Geometry_memory */

/* --- Geometry_message --- */
static PyObject *Geometry_message(GeometryObject *self, PyObject *args)
{
	char *msg=NULL;
	int  color=-1;
	if (!PyArg_ParseTuple(args, "|si", &msg, &color)) return NULL;

	if (msg==NULL)
		return PyString_FromString(self->geometry->message());
	else
		self->geometry->message(msg);

	if (color>=0)
		self->geometry->messageColor = color & 0xFFFFFF;
	Py_RETURN_NONE;
} /* Geometry_message */

/* --- Geometry_object --- */
static PyObject *Geometry_object(GeometryObject *self, PyObject *args)
{
	PyObject *obj;
	GObject *object = NULL;
	char const *var = (char const *)"id";
	PyObject *value=NULL;
	if (!PyArg_ParseTuple(args, "O|sO", &obj, &var, &value)) return NULL;

	if (Py_Check4Pattern(obj)) {
	  //		char *pattern = PyString_AsS
               const char *pattern = PyUnicode_AsUTF8AndSize(obj, llen);
	       ArrayIterator<GObject*> iter(*self->objectList);
		while (iter) {
			object = iter++;
			if (!fnmatch(pattern, object->name(), 0)) {
				PyObject *ret = object->config(var, value);
				Py_XDECREF(ret);	// ret can be NULL use XDECREF
			}
		}
		Py_RETURN_NONE;
	} else
	if (PyList_Check(obj)) {
		for (ssize_t i=0; i<PyList_GET_SIZE(obj); i++) {
			object = Py_Object(self, PyList_GetItem(obj,i));
			if (object==NULL) return NULL;
			PyObject *ret = object->config(var, value);
			Py_XDECREF(ret);	// ret can be NULL use XDECREF
		}
		Py_RETURN_NONE;
	} else
	if (PyTuple_Check(obj)) {
		for (ssize_t i=0; i<PyTuple_GET_SIZE(obj); i++) {
			object = Py_Object(self, PyTuple_GetItem(obj,i));
			if (object==NULL) return NULL;
			PyObject *ret = object->config(var, value);
			Py_XDECREF(ret);	// ret can be NULL use XDECREF
		}
		Py_RETURN_NONE;
	} else {
		object = Py_Object(self, obj);
		if (object) {
			PyObject* ret = object->config(var, value);
			return ret;
		}
	}
	return NULL;
} /* Geometry_object */

/* --- Geometry_regionVar --- */
static PyObject *Geometry_regionVar(GeometryObject *self, GRegion *region,
		const char *var, PyObject *value)
{
	if (!strcmp(var, "name")) {
		if (value==NULL)
			return PyString_FromString(region->name());
		else {
			PyErr_SetString(PyExc_SyntaxError, "Cannot set region name");
			return NULL;
		}
	} else
	if (!strcmp(var, "id")) {
		if (value==NULL)
			return PyInt_FromLong(region->id());
		else {
			PyErr_SetString(PyExc_SyntaxError, "Cannot set region id");
			return NULL;
		}
	} else
	  //	if (!strcmp(var, "rotdefi")) {
	  //		if (value==NULL)
	  //			return PyInt_FromLong(region->rotdefi);
	  //		else
	  //			region->rotdefi = (int)PyInt_AsLong(value);
	if (!strcmp(var, "matrix")) {
		if (value) {
			if (PyList_Check(value)) {
				Matrix4	matrix;
				PyList_AsMatrix4(value, matrix);
				try {
					region->matrix(matrix);
				} catch ( ... ) {
					PyErr_SetString(PyExc_SyntaxError, "Invalid region rotation matrix");
					return NULL;
				}
			} else
				region->clearMatrix();
		} //else return matrix FIXME
	} else
	if (!strcmp(var, "hasmatrix")) {
		return PyBool_FromLong(region->hasMatrix());
	} else
	if (!strcmp(var, "show")) {
		if (value==NULL)
			return PyInt_FromLong(region->show);
		else
			region->show = PyInt_AsLong(value);
	} else
	if (!strcmp(var, "type")) {
		if (value==NULL)
			return PyInt_FromLong((int)region->type());
		else
			region->type((RegionType)PyInt_AsLong(value));
	} else
	if (!strcmp(var, "material")) {
		if (value==NULL) {
			const Material* mat = region->material();
			if (mat)
				return PyString_FromString(mat->name());
			else
				Py_RETURN_NONE;
		} else
			region->material(Py_Material(self, value));
	} else
	if (!strcmp(var, "expr")) {
		if (value==NULL)
			return Py_RegionExpr(region);
		else {
			if (!PyList_Check(value)) {
				PyErr_SetString(PyExc_KeyError,
					"Invalid region expression, list of tuples expected");
				return NULL;
			}
			if (PyList_GET_SIZE(value)>0) region->nextGeneration(self->geometry->nextGeneration());
			for (ssize_t i=0; i<PyList_GET_SIZE(value); i++)
				if (_region_addZone(self, region, PyList_GetItem(value, i)))
					return NULL;
		}
	} else {
		PyErr_Format(PyExc_TypeError, "Invalid type \'%s\' specified", var);
		return NULL;
	}
	Py_RETURN_NONE;
} /* Geometry_regionVar */

/* --- Geometry_region --- */
static PyObject *Geometry_region(GeometryObject *self, PyObject *args)
{
	PyObject *obj;
	GRegion *region = NULL;
	char *var;
	PyObject *value=NULL;
	if (!PyArg_ParseTuple(args, "Os|O", &obj, &var, &value)) return NULL;

	if (Py_Check4Pattern(obj)) {
	  //		char *pattern = PyString_AsString(obj);
               const char *pattern = PyUnicode_AsUTF8AndSize(obj, llen);
		ArrayIterator<GRegion*> iter(self->geometry->regions);
		while (iter) {
			region = iter++;
			if (!fnmatch(pattern, region->name(), 0)) {
				PyObject *ret = Geometry_regionVar(self, region, var, value);
				Py_XDECREF(ret);	// ret can be NULL use XDECREF
			}
		}
		Py_RETURN_NONE;
	} else  {
		region = Py_GRegion(self, obj);
		if (region==NULL) return NULL;
		return Geometry_regionVar(self, region, var, value);
	}
	Py_RETURN_NONE;		// To keep compiler happy
} /* Geometry_region */

/* --- Geometry_rotdefi --- */
static PyObject *Geometry_rotdefi(GeometryObject *self, PyObject *args)
{
	int id;
	PyObject *value=NULL;

	if (!PyArg_ParseTuple(args, "i|O", &id, &value)) return NULL;

	if (value == NULL) {
		// return transformation
	} else {
		if (id<=0) {
			PyErr_SetString(PyExc_ValueError,"Index larger than zero expected");
			return NULL;
		}
		// set transformation
		Matrix4 matrix;
		PyList_AsMatrix4(value, matrix);
		self->geometry->rotdefi(id, matrix);
	}

	Py_RETURN_NONE;
} /* Geometry_rotdefi */

/* --- Geometry_set --- */
static PyObject *Geometry_set(GeometryObject *self, PyObject *args)
{
#define NOTDEFINED	-999999999
	char	*name;
	int	 value=NOTDEFINED;

	if (!PyArg_ParseTuple(args, "s|i", &name, &value)) return NULL;

	if (!strcmp(name,"axislen")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->axisLen);
		else
			self->geometry->axisLen = value;
	} else
	if (!strcmp(name,"backgroundcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->backgroundColor() & FLAG_COLORMASK);
		else {
			self->geometry->backgroundColor(value);
		}
	} else
	if (!strcmp(name,"bodybboxincolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->bodyBBoxInColor());
		else
			self->geometry->bodyBBoxInColor(value);
	} else
	if (!strcmp(name,"bodybboxoutcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->bodyBBoxOutColor());
		else
			self->geometry->bodyBBoxOutColor(value);
	} else
	if (!strcmp(name,"cursor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->cursorSize);
		else {
			self->cursorSize = value;
			if (self->cursorSize < 2*self->cursorMoveSize)
				self->cursorSize = 2*self->cursorMoveSize;
		}
	} else
	if (!strcmp(name,"developer")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->developer);
		else
			self->geometry->developer = value;
	} else
	if (!strcmp(name,"errorcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->errorColor & 0xFFFFFF);
		else
			self->geometry->errorColor = value;
	} else
	if (!strcmp(name,"gridtextcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->gridTextColor & 0xFFFFFF);
		else
			self->geometry->gridTextColor = value;
	} else
	if (!strcmp(name,"lighterlevel")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->lighterLevel);
		else
			self->geometry->lighterLevel = value;
	} else
	if (!strcmp(name,"labelcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->labelColor & 0xFFFFFF);
		else
			self->geometry->labelColor = value;
	} else
	if (!strcmp(name,"latticecolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->latticeColor & 0xFFFFFF);
		else {
			self->geometry->latticeColor = value;
			self->geometry->makeLatticeColor();
		}
	} else
	if (!strcmp(name,"messagecolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->messageColor & 0xFFFFFF);
		else
			self->geometry->messageColor = value;
	} else
	if (!strcmp(name,"regionbboxcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->regionBBoxColor());
		else
			self->geometry->regionBBoxColor(value);
	} else
	if (!strcmp(name,"regionerrorcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->regionErrorColor & 0xFFFFFF);
		else
			self->geometry->regionErrorColor = value;
	} else
	if (!strcmp(name,"regioncolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->regionColor & 0xFFFFFF);
		else
			self->geometry->regionColor = value;
	} else
	if (!strcmp(name,"snapdistance")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(Round(self->snapDistance));
		else
			self->snapDistance = value;
	} else
	if (!strcmp(name,"selectcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->selectColor & 0xFFFFFF);
		else {
			self->geometry->selectColor = value;
			self->geometry->select3DColor(value);
			for (int i=0; i<self->rulers; i++)
				self->ruler[i]->color = value;
		}
	} else
	if (!strcmp(name,"snapangle")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(Round(DEG(self->snapAngle)));
		else
			self->snapAngle = RAD(value);
	} else
	if (!strcmp(name,"trackballsize")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(Round(self->trackballSize*100.0));
		else
			self->trackballSize = (double)value/100.0;
	} else
	if (!strcmp(name,"titlecolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->titleColor & 0xFFFFFF);
		else
			self->geometry->titleColor = value;
	} else
	if (!strcmp(name,"vertexcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->vertexColor & 0xFFFFFF);
		else
			self->geometry->vertexColor = value;
	} else
	if (!strcmp(name,"visiblecolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->visibleColor & 0xFFFFFF);
		else
			self->geometry->visibleColor = value;
	} else
	if (!strcmp(name,"voxelcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->voxelColor & 0xFFFFFF);
		else {
			self->geometry->voxelColor = value;
			self->geometry->makeLatticeColor();
		}
	} else
	if (!strcmp(name,"wireframecolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->wireframeColor());
		else
			self->geometry->wireframeColor(value);
	} else
	if (!strcmp(name,"zonecolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->zoneColor & 0xFFFFFF);
		else
			self->geometry->zoneColor = value;
	} else
	if (!strcmp(name,"zonebboxcolor")) {
		if (value==NOTDEFINED)
			return PyInt_FromLong(self->geometry->zoneBBoxColor());
		else
			self->geometry->zoneBBoxColor(value);
	} else {
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid type option", name);
		return NULL;
	}

	Py_RETURN_NONE;
} /* Geometry_set */

/* --- Geometry_select --- */
static PyObject *Geometry_select(GeometryObject *self, PyObject *args)
{
	char	*name=NULL;
	dword	mask=BIT_SELECT | BIT_FREEZE;

	if (!PyArg_ParseTuple(args, "s|i", &name, &mask)) return NULL;

	if (name[0]==0) Py_RETURN_NONE;

	if (name[0]=='b' || name[0]=='B') {
		int count = 0;
		for (int i=0; i<self->geometry->bodies.size(); i++) {
			if (self->geometry->bodies[i]->show & mask)
				count++;
		}
		return PyInt_FromLong(count);
	} else
	if (name[0]=='r' || name[0]=='R') {
		int count = 0;
		for (int i=0; i<self->geometry->regions.size(); i++) {
			if (self->geometry->regions[i]->show & mask)
				count++;
		}
		return PyInt_FromLong(count);
	}
	Py_RETURN_NONE;
} /* Geometry_select */

/* --- Geometry_setLights --- */
static PyObject *Geometry_setLights(GeometryObject *self)
{
	// Scan objects for the first MAXLIGHT active lights
	ArrayIterator<GObject*> iter(*self->objectList);

	self->geometry->delLights();
	while (iter) {
		Light light;
		GObject *object = iter++;
		if (object->isA() != GLightClass) continue;
		((GLight *)object)->toLight(&light);
		self->geometry->addLight(light);
	}
	if (self->geometry->lights==0)
		self->geometry->defaultLights();
	Py_RETURN_NONE;
} /* Geometry_setLights */

/** Geometry_volumeJob */
struct _VolumeArg {
	GeometryObject* self;
	int		samples;
	PyObject*	obj;

// temporary
	double		volume;
	double		error;
};
static void Geometry_volumeJob(_VolumeArg* arg)
{
	GeometryKernel kernel(*arg->self->geometry, 3);
	kernel.derive();

	arg->volume = 0.0;
	arg->error  = 0.0;

	if (Py_Check4Pattern(arg->obj)) {
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
		GRegion* region = Py_GRegion(arg->self, arg->obj);
		if (region==NULL) {
			arg->self->thread = 0;
			return;
		}
		arg->volume = kernel.volume(region, arg->samples, &arg->volume, &arg->error);
//		return Py_BuildValue("dd",vol, err);
//		return PyFloat_FromDouble(kernel.volume(region, samples, &err));
	}
	arg->self->thread = 0;
} // Geometry_volumeJob

/* --- Geometry_volume --- */
static PyObject *Geometry_volume(GeometryObject *self, PyObject *args)
{
	// FIXME not static move to self
	static _VolumeArg arg;
	PyObject *obj = NULL;
	int samples = -1;

	arg.self = self;
	if (!PyArg_ParseTuple(args, "|Oi", &obj, &samples)) return NULL;

	if (obj==NULL)
		return Py_BuildValue("dd", arg.volume, arg.error);
	else
	if (self->thread) {
		PyErr_SetString(PyExc_ValueError,"Cannot submit a volume job while the previous hasn't finished");
		return NULL;
	} else {
		arg.obj = obj;
		arg.samples = samples;
		return Py_BuildValue("d",pthread_create(&self->thread, NULL, (void*(*)(void*))Geometry_volumeJob, &arg));
	}
} /* Geometry_volume */

/* --- Geometry_voxel --- */
static PyObject *Geometry_voxel(GeometryObject *self, PyObject *args)
{
	char *name;
	PyObject *value=NULL, *opt=NULL;

	if (!PyArg_ParseTuple(args, "s|OO", &name, &value, &opt)) return NULL;
	
	if (!strcmp(name,"load")) {	// load voxel detector
		if (value != NULL) {
		  //			char *filename = PyString_AsString(value);
		  const char *filename = PyUnicode_AsUTF8AndSize(value,llen);
					   cout << " filename " << filename << endl;
		  if (!self->geometry->voxel.load(filename)) {
				PyErr_Format(PyExc_IOError, "Unable to read voxel file \'%s\'",
						filename);
				return NULL;
			}
		}
	} else
	if (!strcmp(name,"xlow")) {	// set x-low position
		if (value != NULL)
			self->geometry->voxel.xlow = PyFloat_AsDouble(value);
		else
			return PyFloat_FromDouble(self->geometry->voxel.xlow);
	} else
	if (!strcmp(name,"ylow")) {	// set y-low position
		if (value != NULL)
			self->geometry->voxel.ylow = PyFloat_AsDouble(value);
		else
			return PyFloat_FromDouble(self->geometry->voxel.ylow);
	} else
	if (!strcmp(name,"zlow")) {	// set z-low position
		if (value != NULL)
			self->geometry->voxel.zlow = PyFloat_AsDouble(value);
		else
			return PyFloat_FromDouble(self->geometry->voxel.zlow);
	} else
	if (!strcmp(name,"xhigh")) {	// set x-high position
		if (value != NULL)
			self->geometry->voxel.xhigh = PyFloat_AsDouble(value);
		else
			return PyFloat_FromDouble(self->geometry->voxel.xhigh);
	} else
	if (!strcmp(name,"yhigh")) {	// set y-high position
		if (value != NULL)
			self->geometry->voxel.yhigh = PyFloat_AsDouble(value);
		else
			return PyFloat_FromDouble(self->geometry->voxel.yhigh);
	} else
	if (!strcmp(name,"zhigh")) {	// set z-high position
		if (value != NULL)
			self->geometry->voxel.zhigh = PyFloat_AsDouble(value);
		else
			return PyFloat_FromDouble(self->geometry->voxel.zhigh);
	} else
	if (!strcmp(name,"matrix")) {	// set transformation matrix for voxel
		if (value) {
			if (PyList_Check(value)) {
				Matrix4	matrix;
				PyList_AsMatrix4(value, matrix);
				self->geometry->voxel.matrix(matrix);
			} else
				self->geometry->voxel.clearMatrix();
		} //else return matrix
	} else
	if (!strcmp(name,"get")) {	// get kreg[data[x,y,z]] value (region number)
		if (value != NULL) {
			if (!PyTuple_Check(value) || PyTuple_GET_SIZE(value)!=3) {
				PyErr_SetString(PyExc_TypeError,"Tuple with voxel coordinates was expected");
				return NULL;
			}
			double x = PyFloat_AsDouble(PyTuple_GetItem(value,0));
			double y = PyFloat_AsDouble(PyTuple_GetItem(value,1));
			double z = PyFloat_AsDouble(PyTuple_GetItem(value,2));
			return PyInt_FromLong(self->geometry->voxel.get(x,y,z));
		}
		// no setting is allowed
	} else
	if (!strcmp(name,"data")) {	// get data[x,y,z] value (HU number) WARNING!
		if (value != NULL) {
			if (!PyTuple_Check(value) || PyTuple_GET_SIZE(value)!=3) {
				PyErr_SetString(PyExc_TypeError,"Tuple with voxel coordinates was expected");
				return NULL;
			}
			double x = PyFloat_AsDouble(PyTuple_GetItem(value,0));
			double y = PyFloat_AsDouble(PyTuple_GetItem(value,1));
			double z = PyFloat_AsDouble(PyTuple_GetItem(value,2));
			return PyInt_FromLong(self->geometry->voxel.data(x,y,z));
		}
		// no setting is allowed
	} else
	if (!strcmp(name,"index")) {	// return (i,j,k) position of (x,y,z)
		if (value != NULL) {
			if (!PyTuple_Check(value) || PyTuple_GET_SIZE(value)!=3) {
				PyErr_SetString(PyExc_TypeError,"Tuple with voxel coordinates was expected");
				return NULL;
			}
			double x = PyFloat_AsDouble(PyTuple_GetItem(value,0));
			double y = PyFloat_AsDouble(PyTuple_GetItem(value,1));
			double z = PyFloat_AsDouble(PyTuple_GetItem(value,2));
			return Py_BuildValue("iii",
					self->geometry->voxel.voxeli(x),
					self->geometry->voxel.voxelj(y),
					self->geometry->voxel.voxelk(z));
		}
		// no setting is allowed
	} else
	if (!strcmp(name,"no")) {	// return no (number of regions)
		return PyInt_FromLong(self->geometry->voxel.no);
	} else
	if (!strcmp(name,"mo")) {	// return mo (number of unique voxels)
		return PyInt_FromLong(self->geometry->voxel.mo);
	} else
	if (!strcmp(name,"free")) {
		self->geometry->voxel.cleanup();
	} else {
		PyErr_Format(PyExc_SyntaxError, "'%s' is not a valid option", name);
		return NULL;
	}

	self->geometry->voxel.calcLimits();

	Py_RETURN_NONE;
} /* Geometry_voxel */

/* --- Geometry_lock --- */
static PyObject* Geometry_lock(GeometryObject *self)
{
	self->geometry->lockWrite();
	Py_RETURN_NONE;
} /* Geometry_lock */

/* --- Geometry_unlock --- */
static PyObject* Geometry_unlock(GeometryObject *self)
{
	self->geometry->unlockWrite();
	Py_RETURN_NONE;
} /* Geometry_unlock */

/* ----------------- Geometry type ------------------- */
static PyMethodDef Geometry_methods[] = {
	{"addBody",	(PyCFunction)Geometry_addBody,	METH_VARARGS,	"Add body"},
	{"addObject",	(PyCFunction)Geometry_addObject,METH_VARARGS,	"add a graphics object (point,line,...)"},
	{"addRegion",	(PyCFunction)Geometry_addRegion,METH_VARARGS,	"Add region"},
	{"addMaterial",	(PyCFunction)Geometry_addMaterial,METH_VARARGS,	"Add material"},
	{"bbox",	(PyCFunction)Geometry_bbox,	METH_VARARGS,	"Return bounding box"},
	{"body",	(PyCFunction)Geometry_body,	METH_VARARGS,	"Get/Set body properties"},
	{"cleanup",	(PyCFunction)Geometry_cleanup,	METH_NOARGS,	"cleanup everything"},
	{"cursor",	(PyCFunction)Geometry_cursor,	METH_VARARGS,	"Get/Set cursor object"},
	{"material",	(PyCFunction)Geometry_material,	METH_VARARGS,	"Set/Get material properties"},
	{"message",	(PyCFunction)Geometry_message,	METH_VARARGS,	"Set/Get message information"},
	{"memory",	(PyCFunction)Geometry_memory,	METH_VARARGS,	"Get or print memory"},
	{"object",	(PyCFunction)Geometry_object,	METH_VARARGS,	"Return/set object parameters"},
	{"region",	(PyCFunction)Geometry_region,	METH_VARARGS,	"Set/Get region properties"},
	//	{"rotdefi",	(PyCFunction)Geometry_rotdefi,	METH_VARARGS,	"Set rotdefi"},
	{"set",		(PyCFunction)Geometry_set,	METH_VARARGS,	"Return/set flag/parameter"},
	{"selection",	(PyCFunction)Geometry_select,	METH_VARARGS,	"Get number of selected items"},
	{"setLights",	(PyCFunction)Geometry_setLights,METH_NOARGS,	"set lights from objects"},
	{"volume",	(PyCFunction)Geometry_volume,	METH_VARARGS,	"Calculate volume of region"},
	{"voxel",	(PyCFunction)Geometry_voxel,	METH_VARARGS,	"Set/Get voxel information"},

	{"lock",	(PyCFunction)Geometry_lock,	METH_NOARGS,	"Lock mutex"},
	{"unlock",	(PyCFunction)Geometry_unlock,	METH_NOARGS,	"UnLock mutex"},
#ifdef MEM
	{"_destroy",	(PyCFunction)Geometry_destroy,	METH_NOARGS,	"Deallocate"},
#endif
	{NULL, NULL, 0, NULL}	// Sentinel
};

PyTypeObject GeometryType = {
  //	PyObject_HEAD_INIT(&PyType_Type)
  //	0,				// ob_size
        PyVarObject_HEAD_INIT(&PyType_Type,0)
        "geoviewer.Geometry",		// tp_name
	sizeof(GeometryObject),		// tp_basicsize
	0,				// tp_itemsize
	(destructor)Geometry_dealloc,	// tp_dealloc
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
	"Geometry object",		// tp_doc
	0,				// tp_traverse
	0,				// tp_clear
	0,				// tp_richcompare
	0,				// tp_weaklistoffset
	0,				// tp_iter
	0,				// tp_iternext
	Geometry_methods,		// tp_methods
	0,//Geometry_members,		// tp_members
	0,				// tp_getset
	0,				// tp_base
	0,				// tp_dict
	0,				// tp_descr_get
	0,				// tp_descr_set
	0,				// tp_dictoffset
	0,				// tp_init
	0,				// tp_alloc
	Geometry_new,			// tp_new
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
