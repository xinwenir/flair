# -*- coding: latin1 -*-
# $Id$
#
# return the include directory for python/tcl
#
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	15-Jun-2010
#
# Mod: 14-Dec-2010 by V.Boccone to support different architectures


import os.path
import sys

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

if sys.argv[1]=="tcl":
	import tkinter
	ver = str(tkinter.TclVersion)
	for inc in ("/usr/include/tcl"+ver,
		    "/usr/include/tcl",
		    "/sw/include/tcl"+ver,
		    "/usr/include"):
		if os.path.isdir(inc):
			print (inc)
			break
elif sys.argv[1]=="tkver":
	import tkinter; print (tkinter.TclVersion)

elif sys.argv[1]=="python":
	print ("%d.%d"%(sys.version_info[0], sys.version_info[1]))

elif sys.argv[1]=="python_lib":
	print (sys.prefix + "/lib/python%d.%d"%(sys.version_info[0], sys.version_info[1]))

elif sys.argv[1]=="python_inc":
	print (sys.prefix + "/include/python%d.%d"%(sys.version_info[0], sys.version_info[1]))

# "System/Library/Frameworks/Tcl.framework/" deprecated
