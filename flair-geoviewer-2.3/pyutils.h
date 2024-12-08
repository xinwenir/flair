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

#ifndef __PYUTILS_H
#define __PYUTILS_H

#include <Python.h>
#include "matrix3.h"
#include "matrix4.h"

// Definitions for older pythons
#ifndef Py_RETURN_NONE
#	define Py_RETURN_NONE  return Py_INCREF(Py_None), Py_None
#	define Py_RETURN_TRUE  return Py_INCREF(Py_True), Py_True
#	define Py_RETURN_FALSE return Py_INCREF(Py_False), Py_False
#endif
#ifndef Py_TYPE
#	define Py_TYPE(ob)	(((PyObject*)(ob))->ob_type)
#endif
#if PY_MAJOR_VERSION >= 3
#define PyInt_FromLong PyLong_FromLong
#define PyInt_AsLong PyLong_AsLong
#define PyInt_FromString PyLong_FromString
#define PyInt_Check PyLong_Check
#define PyString_Check PyUnicode_Check
#define PyString_Size PyBytes_Size
#define PyUnicode_AsUTF8 PyBytes_AsString
#define PyString_FromString PyBytes_FromString
#define PyInt_AsUnsignedLongMask PyLong_AsUnsignedLongMask
#define PyString_FromFormat PyBytes_FromFormat
#define PyString_FromStringAndSize PyBytes_FromStringAndSize
#endif

// function prototypes
bool	Py_GetBool(PyObject *obj);
int	Py_GetInt(PyObject *obj);
double	Py_GetFloat(PyObject *obj);
bool	Py_GetUV(PyObject *obj, double *u, double *v);

Vector	Py_GetVector(PyObject *obj);
Point	Py_GetPoint(PyObject *obj);
// FIXME temporary for back-compatibility
inline bool Py_GetVectorOld(PyObject *obj, double *x, double *y, double *z) {
	Vector v = Py_GetVector(obj);
	*x = v.x;
	*y = v.y;
	*z = v.z;
	return PyErr_Occurred()==NULL;
}
inline	PyObject* Py_Vector(const Vector& v)	{ return Py_BuildValue("ddd", v.x, v.y, v.z); }

bool	PyList_AsMatrix3(PyObject *obj, Matrix3& matrix);
bool	PyList_AsMatrix4(PyObject *obj, Matrix4& matrix);
PyObject *PyList_FromMatrix3(const Matrix3& matrix);
PyObject *PyList_FromMatrix4(const Matrix4& matrix);
bool	Py_Check4Pattern(PyObject *obj);

#endif
