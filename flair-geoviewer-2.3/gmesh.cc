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

#include "stl.h"
#include "gmesh.h"
#include "timer.h"
#include "xdraw.h"
#include "pyutils.h"
#include "viewerobject.h"

using namespace std;

/* ================================= GMesh ================================== */
/** set parameter
 * @param key	value to set
 * @param value	python object to set value from
 * @param return true on error and set the python error
 */
PyObject* GMesh::config(const char *key, PyObject *value)
{
	if (!strcmp(key, "file")) {	// load an STL file
		if (value && PyString_Check(value)) {
			STL stl;
			//			stl.read(PyString_AsString(value), mesh);
			stl.read(PyUnicode_AsUTF8(value), mesh);
		}
	} else
	if (!strcmp(key, "bodies")) {	// return list of bodies found
		List<GBody*> list;
		mesh.fit(list,0.2);

		PyObject* ret = PyList_New(0);

		ListIterator<GBody*> iter(list);
		while (iter) {
			GBody *body = iter++;
			int n = body->nWhat();
			PyObject* obj = PyList_New(1+n);
			PyList_SetItem(obj, 0, PyString_FromString(body->typeStr()));
			double what[30];
			body->getWhat(what);
			for (int i=0; i<n; i++)
				PyList_SetItem(obj, i+1, PyFloat_FromDouble(what[i]));
			PyList_Append(ret, obj);
			delete body;
		}
		return ret;
	} else
	if (!strcmp(key, "print")) {	// load an STL file
		cout << mesh << endl;
	} else
		return GObject::config(key, value);

	if (PyErr_Occurred()) return NULL;
	Py_RETURN_NONE;
} // config

/** draw */
void GMesh::draw(ViewerObject* self, Drawable drawable)
{
	GObject::draw(self, drawable);

	// FIXME VERY PRIMITIVE
	// DOUBLE draws lines
	//mesh.clearProcessedFlag();
	for (int i=0; i<mesh.nedges(); i++) {
		const Edge* edge = mesh.edge(i);
		if (edge->show) {
			Vector A = self->kernel->view.xyz2uvw(edge->A() + P);
			Vector B = self->kernel->view.xyz2uvw(edge->B() + P);
			XDrawLine3D(self->display, drawable, self->gc, self->kernel->view, A, B);
		}
		//f->processed(true);
	}
} // draw
