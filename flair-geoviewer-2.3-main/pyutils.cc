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

#include "pyutils.h"

using namespace std;

static const char *trueStr[] = {
	"on", "ON", "On",
	"true", "TRUE", "True" };
//
Py_ssize_t *llenn;
//

/** Py_GetBool
 * Get boolean value from string or integer
 */
bool Py_GetBool(PyObject *obj)
{
	if (PyInt_Check(obj))
		return (bool)PyInt_AsLong(obj);
	else
	if (PyString_Check(obj)) {
	  //		const char *s = PyString_AsString(obj);
		const char *s = PyUnicode_AsUTF8AndSize(obj, llenn);
		for (unsigned i=0; i<sizeof(trueStr)/sizeof(char*); i++)
			if (!strcmp(s, trueStr[i]))
				return true;
	}
	return false;
} /* Py_GetBool */

/** Py_GetInt
 * @param obj	python object to return integer accepts float/int/string
 * @return	integer point number
 */
int Py_GetInt(PyObject *obj)
{
	if (PyFloat_Check(obj))
		return (int)PyFloat_AsDouble(obj);
	else
	if (PyInt_Check(obj))
		return PyInt_AsLong(obj);
	else
	  return atoi((PyUnicode_AsUTF8AndSize(obj,llenn)));
} /* Py_GetInt */

/** Py_GetFloat
 * @param obj	python object to return float accepts float/int/string
 * @return	floating point number
 */
double Py_GetFloat(PyObject *obj)
{
	if (PyFloat_Check(obj))
		return PyFloat_AsDouble(obj);
	else
	if (PyInt_Check(obj))
		return (double)PyInt_AsLong(obj);
	else
	  return atof(PyUnicode_AsUTF8AndSize(obj,llenn));
} /* Py_GetFloat */

/** Py_GetUV
 * @param obj	python object tuple/list of size=2
 * @param u,v	return variables
 * @return	false on error
 */
bool Py_GetUV(PyObject *obj, double *u, double *v)
{
	if (PyTuple_Check(obj) && PyTuple_GET_SIZE(obj)==2) {
		*u = PyFloat_AsDouble(PyTuple_GetItem(obj,0));
		*v = PyFloat_AsDouble(PyTuple_GetItem(obj,1));
		return PyErr_Occurred()==NULL;
	} else
	if (PyList_Check(obj) && PyList_GET_SIZE(obj)==2) {
		*u = PyFloat_AsDouble(PyList_GetItem(obj,0));
		*v = PyFloat_AsDouble(PyList_GetItem(obj,1));
		return PyErr_Occurred()==NULL;
	} else
		PyErr_SetString(PyExc_TypeError, "tuple or list of size 2 expected");
	return true;
} /* Py_GetUV */

/** Py_GetVector
 * @param obj	python object tuple/list of size=3
 * @return	vector
 */
Vector Py_GetVector(PyObject *obj)
{
	Vector v;
	if (PyTuple_Check(obj) && PyTuple_GET_SIZE(obj)==3)
		v.set( PyFloat_AsDouble(PyTuple_GetItem(obj,0)),
			PyFloat_AsDouble(PyTuple_GetItem(obj,1)),
			PyFloat_AsDouble(PyTuple_GetItem(obj,2)));
	else
	if (PyList_Check(obj) && PyList_GET_SIZE(obj)==3)
		v.set( PyFloat_AsDouble(PyList_GetItem(obj,0)),
			PyFloat_AsDouble(PyList_GetItem(obj,1)),
			PyFloat_AsDouble(PyList_GetItem(obj,2)));
	else
		PyErr_SetString(PyExc_TypeError, "tuple or list of size 3 expected");
	return v;
} /* Py_GetVector */


/** Py_GetPoint
 * @param obj	python object tuple/list of size=3
 * @return	point
 */
Point Py_GetPoint(PyObject *obj)
{
	Point p;
	if (PyTuple_Check(obj) && PyTuple_GET_SIZE(obj)==3)
		p.set( PyFloat_AsDouble(PyTuple_GetItem(obj,0)),
			PyFloat_AsDouble(PyTuple_GetItem(obj,1)),
			PyFloat_AsDouble(PyTuple_GetItem(obj,2)));
	else
	if (PyList_Check(obj) && PyList_GET_SIZE(obj)==3)
		p.set( PyFloat_AsDouble(PyList_GetItem(obj,0)),
			PyFloat_AsDouble(PyList_GetItem(obj,1)),
			PyFloat_AsDouble(PyList_GetItem(obj,2)));
	else
		PyErr_SetString(PyExc_TypeError, "tuple or list of size 3 expected");
	return p;
} /* Py_GetPoint */

/* --- PyList_AsMatrix3 --- */
bool PyList_AsMatrix3(PyObject *obj, Matrix3& matrix)
{
	if (!PyList_Check(obj) || PyList_GET_SIZE(obj)!=3) goto ERROR;
	for (int j=0; j<3; j++) {
		PyObject *row = PyList_GetItem(obj, j);
		if (!PyList_Check(row) || PyList_GET_SIZE(row)!=3) goto ERROR;
		for (int i=0; i<3; i++)
			matrix(j,i) = PyFloat_AS_DOUBLE(PyList_GetItem(row,i));
	}
	return true;
ERROR:
	PyErr_SetString(PyExc_TypeError,
			"Invalid Matrix4 list of lists [3x3] expected");
	return false;
} /* PyList_AsMatrix3 */


/* --- PyList_AsMatrix4 --- */
bool PyList_AsMatrix4(PyObject *obj, Matrix4& matrix)
{
	if (!PyList_Check(obj) || PyList_GET_SIZE(obj)!=4) goto ERROR;
	for (int j=0; j<4; j++) {
		PyObject *row = PyList_GetItem(obj, j);
		if (!PyList_Check(row) || PyList_GET_SIZE(row)!=4) goto ERROR;
		for (int i=0; i<4; i++)
			matrix(j,i) = PyFloat_AS_DOUBLE(PyList_GetItem(row,i));
	}
	return true;
ERROR:
	PyErr_SetString(PyExc_TypeError,
			"Invalid Matrix4 list of lists [4x4] expected");
	return false;
} /* PyList_AsMatrix4 */

/* --- PyList_FromMatrix3 --- */
PyObject *PyList_FromMatrix3(const Matrix3& matrix)
{
	PyObject *obj = PyList_New(3);
	for (int j=0; j<3; j++) {
		PyObject *row = Py_BuildValue("[ddd]",
						matrix(j,0),
						matrix(j,1),
						matrix(j,2));
		PyList_SET_ITEM(obj, j, row);
	}
	return obj;
} /* PyList_FromMatrix3 */

/* --- PyList_FromMatrix4 --- */
PyObject *PyList_FromMatrix4(const Matrix4& matrix)
{
	PyObject *obj = PyList_New(4);
	for (int j=0; j<4; j++) {
		PyObject *row = Py_BuildValue("[dddd]",
						matrix(j,0),
						matrix(j,1),
						matrix(j,2),
						matrix(j,3));
		PyList_SET_ITEM(obj, j, row);
	}
	return obj;
} /* PyList_FromMatrix4 */

/* Py_Check4Pattern
 * @param obj	object to check if it holds a string pattern
 * @return true if it has punctuation characters inside
 */
bool Py_Check4Pattern(PyObject *obj)
{
	if (!PyString_Check(obj)) return false;
	//	char *str = PyString_AsString(obj);
	const char *str = PyUnicode_AsUTF8AndSize(obj,llenn);
        int len = strlen(str);
	//	int   len = (int)PyString_Size(obj);
	if (len==0) return false;

	if (memchr(str, '*', len)) return true;
	if (memchr(str, '?', len)) return true;
	if (memchr(str, '[', len)) return true;

	return false;
} /* Py_Check4Pattern */
