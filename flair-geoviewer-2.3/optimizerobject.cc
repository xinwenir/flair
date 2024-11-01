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
 * Author:	Wioletta.Kozlowska@cern.ch
 * Date:	12-Dec-2016
 */

#include <Python.h>

#include <new>
#include <string>

#include "optimizerobject.h"	// it has to be the first include!

using namespace std;

//////////////////////// Optimizer ///////////////////////////
/* --- Optimizer_new --- */
static PyObject *Optimizer_new(PyTypeObject *type, PyObject* /*args*/, PyObject* /*kwds*/)
{
	OptimizerObject* self;
	self = (OptimizerObject*)type->tp_alloc(type, 0);
	if (self == NULL) return NULL;

	self->optimizer = new Optimizer();

	return (PyObject*)self;
} /* Optimizer_new */

/* --- Optimizer_dealloc --- */
static void Optimizer_dealloc(OptimizerObject *self)
{
	delete self->optimizer;
} /* Optimizer_dealloc */

/* --- Optimizer_init --- */
static PyObject *Optimizer_init(OptimizerObject *self, PyObject *args)
{
	char *in  = NULL;
	char *out = NULL;
	bool isBioOpt;
	int  ncells;
	char *cellfile = NULL;

	if (!PyArg_ParseTuple(args, "ssiis", &in, &out, &isBioOpt, &ncells, &cellfile )) return NULL;

	string sin   = in;
	string sout  = out;
	string scell = cellfile;

	try {
		self->optimizer->init(sin, sout, isBioOpt, ncells, scell);
	} catch ( ... ) {
		PyErr_Format(PyExc_TypeError,
			"Very important \'%s\' specified", in);
		return NULL;
	}

	Py_RETURN_NONE;
} /* Optimizer_init */

/* ----------------- Optimizer type ------------------- */
static PyMethodDef Optimizer_methods[] = {
	{"init",	(PyCFunction)Optimizer_init,	METH_VARARGS,	"Initialize optimizer blah blah"},
//	{"lock",	(PyCFunction)Optimizer_lock,	METH_NOARGS,	"Lock mutex"},
	{NULL, NULL, 0, NULL}	// Sentinel
};

PyTypeObject OptimizerType = {
  //	PyObject_HEAD_INIT(&PyType_Type)
  //	0,				// ob_size
        PyVarObject_HEAD_INIT(&PyType_Type,0)
        "geoviewer.Optimizer",		// tp_name
	sizeof(OptimizerObject),	// tp_basicsize
	0,				// tp_itemsize
	(destructor)Optimizer_dealloc,	// tp_dealloc
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
	"Optimizer object",		// tp_doc
	0,				// tp_traverse
	0,				// tp_clear
	0,				// tp_richcompare
	0,				// tp_weaklistoffset
	0,				// tp_iter
	0,				// tp_iternext
	Optimizer_methods,		// tp_methods
	0,//Optimizer_members,		// tp_members
	0,				// tp_getset
	0,				// tp_base
	0,				// tp_dict
	0,				// tp_descr_get
	0,				// tp_descr_set
	0,				// tp_dictoffset
	0,				// tp_init
	0,				// tp_alloc
	Optimizer_new,			// tp_new
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
