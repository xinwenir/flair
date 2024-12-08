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

#include "viewerobject.h"
#include "geometryobject.h"
#include "optimizerobject.h"

#define VERSION "2.3-0e"

#ifndef PyMODINIT_FUNC		// declarations for DLL import/export
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
//initgeoviewer(void)
PyInit_geoviewer(void)
{
#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "geoviewer",     /* m_name */
        "Geometry viewer extension type",  /* m_doc */
        -1,                  /* m_size */
        NULL,                /* m_methods */
        NULL,                /* m_reload */
        NULL,                /* m_traverse */
        NULL,                /* m_clear */
        NULL,                /* m_free */
    };
#endif
    PyObject* m;

	//GeometryType.tp_new = PyType_GenericNew;
	if (PyType_Ready(&GeometryType) < 0)  return NULL;
	if (PyType_Ready(&ViewerType) < 0)    return NULL;
	if (PyType_Ready(&OptimizerType) < 0) return NULL;
#if PY_MAJOR_VERSION >= 3
	m = PyModule_Create(&moduledef);
#else

	m = Py_InitModule3("geoviewer", NULL /*module methods*/,
		"Geometry viewer extension type.");
#endif
	Py_INCREF(&GeometryType);
	Py_INCREF(&ViewerType);
	Py_INCREF(&OptimizerType);

	PyModule_AddStringConstant(m, "__author__",  "Vasilis Vlachoudis");
	PyModule_AddStringConstant(m, "__email__",   "Vasilis.Vlachoudis@cern.ch");
	PyModule_AddStringConstant(m, "__version__",  VERSION);

	PyModule_AddObject(m, "Geometry",  (PyObject *)&GeometryType);
	PyModule_AddObject(m, "Viewer",    (PyObject *)&ViewerType);
	PyModule_AddObject(m, "Optimizer", (PyObject *)&OptimizerType);
        return m ;
} // initgeoviewer
