#!/usr/bin/python
# -*- coding: latin1 -*-
#
# Copyright and User License
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright 2006-2019 CERN and INFN
# 
#
# Please consult the LICENSE file for the license 
#
# DISCLAIMER
# ~~~~~~~~~~
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY, OF
# SATISFACTORY QUALITY, AND FITNESS FOR A PARTICULAR PURPOSE
# OR USE ARE DISCLAIMED. THE COPYRIGHT HOLDERS AND THE
# AUTHORS MAKE NO REPRESENTATION THAT THE SOFTWARE AND
# MODIFICATIONS THEREOF, WILL NOT INFRINGE ANY PATENT,
# COPYRIGHT, TRADE SECRET OR OTHER PROPRIETARY RIGHT.
#
# LIMITATION OF LIABILITY
# ~~~~~~~~~~~~~~~~~~~~~~~
# THE COPYRIGHT HOLDERS AND THE AUTHORS SHALL HAVE NO
# LIABILITY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
# CONSEQUENTIAL, EXEMPLARY, OR PUNITIVE DAMAGES OF ANY
# CHARACTER INCLUDING, WITHOUT LIMITATION, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES, LOSS OF USE, DATA OR PROFITS,
# OR BUSINESS INTERRUPTION, HOWEVER CAUSED AND ON ANY THEORY
# OF CONTRACT, WARRANTY, TORT (INCLUDING NEGLIGENCE), PRODUCT
# LIABILITY OR OTHERWISE, ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGES.
#
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	14-Oct-2011

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import re
import sys
import math
import string
from math import *
from log import say

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import bmath
import undo
import Quad
import tkFlair
import tkExtra

import Input
import Manual
import Project

_NONE   = "-none-"
_GROUPS = ["all", "plane", "box", "sphere", "cylinder", "quadric"]

# I have to move variables in another module and import them from Viewer and Body
_VIEW   = ["all", "3D", "Red", "Green", "Blue", "Magenta"]

# FIXME should be defined in tkFlair
SMALL = 1E-10
BIG   = 1E10

_NAMEPAT = re.compile(r"^(.+?)\d+$")

#-------------------------------------------------------------------------------
# Round number to closest integer if needed
#-------------------------------------------------------------------------------
def fmt(x):
	if abs(x)<SMALL:
		return 0.0
	else:
		return bmath.format(x,15)

#-------------------------------------------------------------------------------
# Create a new body aligned to a certain viewport
#-------------------------------------------------------------------------------
def newBody(tag, xyz, size, matrix):
	# Find direction if needed
	u = bmath.Vector(matrix[0][0], matrix[1][0], matrix[2][0])
	v = bmath.Vector(matrix[0][1], matrix[1][1], matrix[2][1])
	w = bmath.Vector(matrix[0][2], matrix[1][2], matrix[2][2])

	# Create basic whats
	what = ["body"]
	if   tag == "XYP": what.append(xyz[2])
	elif tag == "XZP": what.append(xyz[1])
	elif tag == "YZP": what.append(xyz[0])
	elif tag == "PLA":
		what.extend(size*u)
		what.extend(list(xyz))

	elif tag == "RPP":
		if abs(abs(w[0])-1.0) < 1.0E-5:
			what.append(xyz[0]-size/2.0)
			what.append(xyz[0]+size/2.0)
		else:
			what.append(xyz[0])
			what.append(xyz[0]+size)

		if abs(abs(w[1])-1.0) < 1.0E-5:
			what.append(xyz[1]-size/2.0)
			what.append(xyz[1]+size/2.0)
		else:
			what.append(xyz[1])
			what.append(xyz[1]+size)

		if abs(abs(w[2])-1.0) < 1.0E-5:
			what.append(xyz[2]-size/2.0)
			what.append(xyz[2]+size/2.0)
		else:
			what.append(xyz[2])
			what.append(xyz[2]+size)

	elif tag in ("BOX", "REC"):
		what.extend(list(xyz))
		what.extend( size*u)
		what.extend( size*v)
		what.extend(-size*w)

	elif tag in ("WED", "RAW"):
		what.extend(list(xyz))
		what.extend( size*u)
		what.extend( size*v)
		what.extend(-size*w)

	elif tag == "SPH":
		what.extend(list(xyz))
		what.append(size/2.0)

	elif tag == "ELL":
		what.append(xyz[0] - size*u[0])
		what.append(xyz[1] - size*u[1])
		what.append(xyz[2] - size*u[2])
		what.append(xyz[0] + size*u[0])
		what.append(xyz[1] + size*u[1])
		what.append(xyz[2] + size*u[2])
		what.append(3.0*size)

	elif tag in ("RCC", "TRC"):
		what.extend(list(xyz))
		what.extend(size*u)
		what.append(size)
		if tag == "TRC":
			what.append(size/2.0)

	elif tag in ("XCC", "XEC"):
		what.append(xyz[1])
		what.append(xyz[2])
		what.append(size)
		if tag=="XEC": what.append(size*2.0)

	elif tag in ("YCC", "YEC"):
		what.append(xyz[2])
		what.append(xyz[0])
		what.append(size)
		if tag=="YEC": what.append(size*2.0)

	elif tag in ("ZCC", "ZEC"):
		what.append(xyz[0])
		what.append(xyz[1])
		what.append(size)
		if tag=="ZEC": what.append(size*2.0)

	elif tag == "QUA":
		what.append(1.0/size**2)
		what.append(1.0/size**2)
		what.append(1.0/size**2)
		what.append(0.0)
		what.append(0.0)
		what.append(0.0)
		what.append(-2.0*xyz[0]/size**2)
		what.append(-2.0*xyz[1]/size**2)
		what.append(-2.0*xyz[2]/size**2)
		what.append((xyz[0]**2 + xyz[1]**2 + xyz[2]**2)/size**2-1.0)

	return Input.Card(tag, what)

#-------------------------------------------------------------------------------
# return a PLA on the viewing plane, assuming w pointing inside
# @param pts	a list of points (u,v) on the viewing plane
# @return None or (tag, [whats], err)
#-------------------------------------------------------------------------------
def findPlaneUV(pts, matrix=None):
	N = len(pts)

	err = 0.0
	if N==1:
		# with one point we return
		a = b = d = 0.0
		c = 1.0
	else:
		A = bmath.Matrix(pts)
		B = bmath.Matrix([-1.0]*N)

		# default plane
		a = b = c = 0.0; d = 1.0
		try:
			# Check for line equation a*u + b*v + 1 = 0
			# and the plane equation a*u + b*v + 0*w + 1 = 0
			a,b = bmath.solveOverDetermined(A,B)

			if abs(a)>BIG or abs(b)>BIG: raise
			if abs(a)<SMALL and abs(b)<SMALL: raise
		except:
			# search for a solution a*u + b*v = 0 <=> a*u = -v
			ATA = ATB = 0.0
			for u,v in pts:
				ATA += u*u
				ATB -= u*v
			if ATA>SMALL:
				a = ATB/ATA
				b = 1.0
				d = 0.0
			else:
				# check if u=0 or v=0
				su = sv = 0.0
				for u,v in pts:
					su += abs(u)
					sv += abs(v)

				if abs(su)<SMALL and abs(sv)<SMALL: return None

				d = 0.0
				if abs(su)<SMALL:
					a = 1.0
				elif abs(sv)<SMALL:
					b = 1.0
				else:
					return None

		# find error
		for u,v in pts:
			err += (a*u+b*v+d)**2
		err = sqrt(err)/float(N)

	# Rotate if matrix is present
	if matrix is not None:
		q = Quad.Quad(Cx=a,Cy=b,Cz=c,C=d)
		q.transform(matrix)
		a = q.Cx
		b = q.Cy
		c = q.Cz
		d = q.C

	# Avoid rounding
	if abs(a) < SMALL: a = 0.0
	if abs(b) < SMALL: b = 0.0
	if abs(c) < SMALL: c = 0.0
	if abs(d) < SMALL: d = 0.0

	return "PLA", (a,b,c,d), err

#-------------------------------------------------------------------------------
# @param pts		a list of points that form the plane
#-------------------------------------------------------------------------------
def findPlane(pts):
	N = len(pts)

	if N < 3: return None

	A = bmath.Matrix(pts)
	B = bmath.Matrix([-1.0]*N)

	# default plane
	a = b = c = 0.0; d = 1.0
	try:	# try generic plane
		a, b, c = bmath.solveOverDetermined(A,B)
		if abs(a)>BIG or abs(b)>BIG or abs(c)>BIG: raise
		if abs(a)<SMALL and abs(b)<SMALL and abs(c)<SMALL: raise
	except:
		# search for a solution a*x + b*y + c*z = 0 <=> a*x + b*y = -z
		A = []
		B = []
		for x,y,z in pts:
			A.append([x,y])
			B.append(-z)

		A = bmath.Matrix(A)
		B = bmath.Matrix(B)

		d = 0.0
		try:
			a, b = bmath.solveOverDetermined(A,B)
			if abs(a)>BIG or abs(b)>BIG: raise
			if abs(a)<SMALL and abs(b)<SMALL: raise
			c = 1.0
		except:
			# Check for a solution of x=0 or y=0 or z=0
			sx = sy = sz = 0.0
			for x,y,z in pts:
				sx += abs(x)
				sy += abs(y)
				sz += abs(z)

			if abs(sx)<SMALL and abs(sy)<SMALL and abs(sz)<SMALL:
				return None

			# check if all points fall on a line
			#X = [p[0] for p in pts]
			#Y = [p[1] for p in pts]
			#Z = [p[2] for p in pts]
			#r1 = bmath.linear(X,Y)
			#r2 = bmath.linear(X,Z)
			#if r1 is not None and r2 is not None and \
			#   abs(r1[2])>0.999999 and abs(r2[2])>0.999999: return None
			#say("r1=",r1)
			#say("r2=",r2)

			if abs(sx)<SMALL:
				a = 1.0
			elif abs(sy)<SMALL:
				b = 1.0
			elif abs(sz)<SMALL:
				c = 1.0
			else:
				# Check for 1) a*x + b*y = 0 <=> a*x = -y
				#        or 2) a*x + c*z = 0 <=> a*x = -z
				#        or 3) b*y + c*z = 0 <=> b*y = -z
				# search for a solution a*u + b*v = 0 <=> a*u = -v
				ATA1 = ATB1 = 0.0
				ATA2 = ATB2 = 0.0
				ATA3 = ATB3 = 0.0

				for x,y,z in pts:
					ATA1 += x*x
					ATB1 -= x*y

					ATB2 -= x*z

					ATA3 += y*y
					ATB3 -= y*z
				ATA2 = ATA1

				if ATA1 > SMALL:
					a = ATB1/ATA1
					b = 1.0
				elif ATA3 > SMALL:
					b = ATB3/ATA3
					c = 1.0
				else:
					return None

			# FIXME Check if all points fall on a line

	# find error
	err = 0.0
	for x,y,z in pts:
		err += (a*x+b*y+c*z+d)**2
	err = sqrt(err)/float(N)

	return "PLA", (a,b,c,d), err

#-------------------------------------------------------------------------------
# Convert a generic plane to a more appropriate form
# @param tagwhat	("PLA", [a,b,c,d], err)
# @return a ("PLA", "XYP", "XZP", "YZP", [x|y|z|(nx,ny,nz,px,py,pz)], err)
#-------------------------------------------------------------------------------
def findBestPlane(tagwhat):
	if tagwhat is None: return None
	tag,what,err = tagwhat
	a,b,c,d = what

	# prefer the positive solution
	if a<=SMALL and b<=SMALL and c<=SMALL:
		a = -a
		b = -b
		c = -c
		d = -d

	if abs(a)<SMALL and abs(b)<SMALL:
		w = -d/c
		return "XYP", [fmt(w)], err
	elif abs(a)<SMALL and abs(c)<SMALL:
		w = -d/b
		return "XZP", [fmt(w)], err
	elif abs(b)<SMALL and abs(c)<SMALL:
		w = -d/a
		return "YZP", [fmt(w)], err

	# searh an intersection of the plane along  various intersection
	# for the most simpler
	solutions = []
	try: solutions.append((-d/(a+b+c), -d/(a+b+c), -d/(a+b+c)))
	except: pass
	try: solutions.append((-d/a, 0.0, 0.0))
	except: pass
	try: solutions.append((0.0, -d/b, 0.0))
	except: pass
	try: solutions.append((0.0, 0.0, -d/c))
	except: pass

	m = -1
	l = 10000
	for i,(x,y,z) in enumerate(solutions):
		s = "%s %s %s"%(bmath.format(x,22), bmath.format(y,22), bmath.format(z,22))
		if len(s)<l:
			m = i
			l = len(s)

	# normalize normal
	s = sqrt(a*a + b*b + c*c)
	a /= s
	b /= s
	c /= s

	return "PLA", list(map(fmt,[a, b, c, solutions[m][0], solutions[m][1], solutions[m][2]])), err

#-------------------------------------------------------------------------------
# Return a SPHere UV
# @param pts	a list of points (u,v) on the viewing plane
# @return None or (tag, [whats], err)
#-------------------------------------------------------------------------------
def findSphereUV(pts, matrix=None):
	N = len(pts)
	if N < 3: return None

	A = []
	B = []
	for x,y in pts:
		A.append([2.0*x, 2.0*y, -1.0])
		B.append([x*x+y*y])
	A = bmath.Matrix(A)
	B = bmath.Matrix(B)

	try:
		x,y,C = bmath.solveOverDetermined(A,B)
	except:
		return None

	R2 = C + x*x + y*y
	if R2 < SMALL: return None

	# find error
	#err = 0.0
	#for x,y,z in pts:
	#	err += ()**2
	#err = sqrt(err)/float(N)

	# Rotate if matrix is present
	p = bmath.Vector(x,y,0.0)
	if matrix is not None:
		p = matrix * p

	return "SPH", list(map(fmt,[ p[0], p[1], p[2], sqrt(R2)])), 0.0

#-------------------------------------------------------------------------------
# Return a SPHere
#-------------------------------------------------------------------------------
def findSphere(pts):
	N = len(pts)
	if N < 4: return None

	A = []
	B = []
	for x,y,z in pts:
		A.append([2.0*x, 2.0*y, 2.0*z, -1.0])
		B.append([x*x+y*y+z*z])
	A = bmath.Matrix(A)
	B = bmath.Matrix(B)

	try:
		x,y,z,C = bmath.solveOverDetermined(A,B)
	except:
		return None

	R2 = C + x*x + y*y + z*z
	if R2 < SMALL: return None

	return "SPH", list(map(fmt, [x, y, z, sqrt(R2)])), 0.0

#-------------------------------------------------------------------------------
# Return a Cylinder [XYZ]CC UV
# @param pts	a list of points (u,v) on the viewing plane
# @return None or (tag, [whats], err)
#-------------------------------------------------------------------------------
def findCylinderUV(pts, matrix=None):
	N = len(pts)
	if N < 3: return None

	A = []
	B = []
	for x,y in pts:
		A.append([2.0*x, 2.0*y, -1.0])
		B.append([x*x+y*y])
	A = bmath.Matrix(A)
	B = bmath.Matrix(B)

	try:
		x,y,C = bmath.solveOverDetermined(A,B)
	except:
		return None

	R2 = C + x*x + y*y
	if R2 < SMALL: return None

	# find error
	#err = 0.0
	#for x,y,z in pts:
	#	err += ()**2
	#err = sqrt(err)/float(N)

	# Rotate if matrix is present
	p = bmath.Vector(x,y,0.0)
	z = bmath.Vector(0.0,0.0,1.0)
	if matrix is not None:
		p = matrix * p
		z = matrix * z

	R = sqrt(R2)

	if   abs(z[0])<=SMALL and abs(z[1])<=SMALL and abs(abs(z[2])-1.0)<=SMALL:
		return "ZCC", list(map(fmt, [p[0], p[1], R])), 0.0
	elif abs(z[0])<=SMALL and abs(abs(z[1])-1.0)<=SMALL and abs(z[2])<=SMALL:
		return "YCC", list(map(fmt, [p[2], p[0], R])), 0.0
	elif abs(abs(z[0])-1.0)<=SMALL and abs(z[1])<=SMALL and abs(z[2])<=SMALL:
		return "XCC", list(map(fmt, [p[1], p[2], R])), 0.0
	else:
		q = Quad.Quad(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -2.0*x, -2.0*y, 0.0, -C)
		q.transform(matrix)
		return "QUA", list(map(fmt, [q.Cxx, q.Cyy, q.Czz, q.Cxy, q.Cxz, q.Cyz, q.Cx, q.Cy, q.Cz, q.C])), 0.0

#-------------------------------------------------------------------------------
# Return a Cylinder [XYZ]CC
#-------------------------------------------------------------------------------
def findCylinder(pts):
	N = len(pts)
	if N < 3: return None

	AX = []
	AY = []
	AZ = []
	B = bmath.Matrix([1.0]*N)
	for x,y,z in pts:
		AX.append([y*y+z*z, -2.0*y, -2.0*z])
		AY.append([x*x+z*z, -2.0*x, -2.0*z])
		AZ.append([x*x+y*y, -2.0*x, -2.0*y])
	AX = bmath.Matrix(AX)
	AY = bmath.Matrix(AY)
	AZ = bmath.Matrix(AZ)

	try:
		C,y,z = bmath.solveOverDetermined(AX,B)
		y /= C
		z /= C
		R2 = 1.0/C+(y*y+z*z)
		if R2<0.0:
			return "XCC", [-y,-z,sqrt(-R2)]
		else:
			return "XCC", [ y, z,sqrt( R2)]
	except:
		pass

	try:
		C,x,z = bmath.solveOverDetermined(AY,B)
		x /= C
		z /= C
		R2 = 1.0/C+(x*x+z*z)
		if R2<0.0:
			return "YCC", [-z,-x,sqrt(-R2)]
		else:
			return "YCC", [ z, x,sqrt( R2)]
	except:
		pass

	try:
		C,x,y = bmath.solveOverDetermined(AZ,B)
		x /= C
		y /= C
		R2 = 1.0/C+(x*x+y*y)
		if R2<0.0:
			return "ZCC", [-x,-y,sqrt(-R2)]
		else:
			return "ZCC", [ x, y,sqrt( R2)]
	except:
		pass

	return None

#-------------------------------------------------------------------------------
# Return an Elliptical Cylinder [XYZ]EC
#-------------------------------------------------------------------------------
def findEllipticalCylinder(pts):
	N = len(pts)
	if N < 4: return None

	AX = []
	AY = []
	AZ = []
	B = bmath.Matrix([1.0]*N)
	for x,y,z in pts:
		AX.append([y*y, z*z, -2.0*y, -2.0*z])
		AY.append([x*x, z*z, -2.0*x, -2.0*z])
		AZ.append([x*x, y*y, -2.0*x, -2.0*y])
	AX = bmath.Matrix(AX)
	AY = bmath.Matrix(AY)
	AZ = bmath.Matrix(AZ)

	try:
		a,b,y,z = bmath.solveOverDetermined(AX,B)
		y /= a
		z /= b
		Ry = sqrt(y*y + b/a * z*z + 1.0/a)
		Rz = sqrt(a/b)*Rx
		return "YEC", [y,z,Ry,Rz]
	except:
		pass

	try:
		a,b,x,z = bmath.solveOverDetermined(AY,B)
		x /= a
		z /= b
		Rx = sqrt(x*x + b/a * z*z + 1.0/a)
		Rz = sqrt(a/b)*Rx
		return "YEC", [x,z,Rx,Rz]
	except:
		pass

	try:
		a,b,x,y = bmath.solveOverDetermined(AZ,B)
		x /= a
		y /= b
		Rx = sqrt(x*x + b/a * y*y + 1.0/a)
		Ry = sqrt(a/b)*Rx
		return "ZEC", [x,y,Rx,Ry]
	except:
		pass

	return None

#-------------------------------------------------------------------------------
# Return a QUAdratic
#-------------------------------------------------------------------------------
def findQUA(pts):
	N = len(pts)
	if N < 5: return None

	B = bmath.Matrix([1.0]*N)

	A  = []
	AX = []		# ignore x
	AY = []		# ignore y
	AZ = []		# ignore z
	for x,y,z in pts:
		A.append([x*x, y*y, z*z, x*y, x*z, y*z, x, y, z])
		AX.append([y*y, z*z, y*z, y, z])
		AY.append([x*x, z*z, x*z, x, z])
		AZ.append([x*x, y*y, x*y, x, y])
	A = bmath.Matrix(A)

	try:
		X = bmath.solveOverDetermined(A,B)
		return "QUA", X
	except:
		pass

	try:
		X = bmath.solveOverDetermined(AX,B)
		return "QUX", X
	except:
		pass

	try:
		X = bmath.solveOverDetermined(AY,B)
		return "QUY", X
	except:
		pass

	try:
		Z = bmath.solveOverDetermined(AZ,B)
		return "QUZ", X
	except:
		pass

	return None

#===============================================================================
# Geometry Auto Body
# Detect all possible body types and propose to the user
#===============================================================================
class GeometryAutoBody(Toplevel):
	activepage = None
	width  = -1
	height = -1
	_PAGES = None

	def __init__(self, master, editor, project):
		Toplevel.__init__(self, master)
		self.editor  = editor
		self.project = project
		self.title("Geometry Auto Body Creation")
		self.layer = None

		# --- Layer ---
		row = 0
		Label(self, text="Name:").grid(row=row, column=0, sticky=E)
		self.name = Entry(self, background="White", width=9)
		self.name.insert(0, self.editor.findUniqueName(Input.Card.BODY,"body"))
		self.name.grid(row=row, column=1,sticky=EW)
		tkExtra.Balloon.set(self.name,
			"Proposed name for the newly created body")
		b = Button(self, text="A", command=self.autoname,
				padx=2,pady=1)
		b.grid(row=row, column=2,sticky=EW)
		tkExtra.Balloon.set(b, "Assign automatically a name to the body")

		# ---
		Label(self, text="# Points:").grid(row=row, column=3, sticky=E)
		self.npoints = Label(self, text=_NONE, background="White", anchor=W)
		self.npoints.grid(row=row, column=4, sticky=EW)
		tkExtra.Balloon.set(self.npoints, "Number of selected points")

		# ---
		row += 1
		Label(self, text="Comment:").grid(row=row, column=0, sticky=NE)
		self.comment = Text(self, background="White", width=30, height=3)
		self.comment.grid(row=row, column=1,columnspan=4,sticky=NSEW)
		tkExtra.Balloon.set(self.comment,
			"Comment describing the new body")

		# ---
		row += 1
		Label(self, text="View:").grid(row=row, column=0, sticky=E)
		self.view = tkExtra.Combobox(self, width=10)
		self.view.grid(row=row, column=1, columnspan=2, sticky=EW)
		self.view.fill(_VIEW)
		self.view.set(_VIEW[0])
		# Hack
		for i in range(2,len(_VIEW)):
			self.view._listbox.itemconfig(i,
				foreground="Dark%s"%(_VIEW[i]))
		tkExtra.Balloon.set(self.view,
			"Show bodies that are parallel or perpendicular to a viewport")

		Label(self, text="Group:").grid(row=row, column=3, sticky=E)
		self.group = tkExtra.Combobox(self, width=10)
		self.group.grid(row=row, column=4, sticky=EW)
		self.group.fill(_GROUPS)
		self.group.set(_GROUPS[0])
		tkExtra.Balloon.set(self.group,
			"Show only bodies belonging to a certain group")

		# ---
		row += 1
		b = Button(self, text="Bodies:",
				command = self.forceRefresh,
				image=tkFlair.icons["refresh"],
				compound=LEFT, anchor=W)
		b.grid(row=row, column=0, columnspan=5, sticky=EW)
		tkExtra.Balloon.set(b, "Click to refresh")
		row += 1
		self.text = Text(self, background="White", takefocus=1, width=40)
		self.text.grid(row=row, column=0, columnspan=5, sticky=NSEW)
		sbv = tkExtra.AutoScrollbar(self, takefocus=0)
		sbv.grid(row=row, column=2, sticky='ns')
		sbh = tkExtra.AutoScrollbar(self, takefocus=0,
				orient=HORIZONTAL)
		sbh.grid(row=row+1, column=0, columnspan=2, sticky='ew')

		self.text.config(yscrollcommand=sbv.set, xscrollcommand=sbh.set)
		sbv.config(command=self.text.yview)
		sbh.config(command=self.text.xview)
		self.text["wrap"]  = NONE
		self.text["state"] = DISABLED

		self.grid_rowconfigure(row, weight=1)
		self.grid_columnconfigure(1, weight=1)
		self.grid_columnconfigure(4, weight=1)

		# ---
		row += 1
		frame = Frame(self)
		frame.grid(row=row, column=0, columnspan=5, sticky=EW)

		# buttons
		Button(frame, text='Help', image=tkFlair.icons["info"],
			compound=LEFT, takefocus=FALSE,
			command=self.help).pack(side=LEFT, padx=5, pady=5)
		Button(frame, text='Close', image=tkFlair.icons["x"],
			compound=LEFT, takefocus=FALSE,
			command=self.close).pack(side=RIGHT, padx=5, pady=5)
		Button(frame, text='Insert', image=tkFlair.icons["add"],
			compound=LEFT, takefocus=FALSE,
			command=self.insert).pack(side=RIGHT, padx=5, pady=5)

		self.deiconify()
		self.wait_visibility()
		if GeometryAutoBody.width > 0:
			self.geometry("%dx%d" % \
				(GeometryAutoBody.width, GeometryAutoBody.height))
		else:
			self.geometry("%dx%d" % \
				(self.winfo_width(), self.winfo_height()))

		self.protocol("WM_DELETE_WINDOW", self.close)
		self.bind("<Escape>"            , self.close)
		self.bind("<Control-Key-q>"     , self.close)
		self.bind("<Control-Key-w>"     , self.close)
		self.bind("<F1>"                , self.help)
		self.bind("<FocusIn>"           , self.refresh)
		self.bind("<Control-Key-r>"     , self.refresh)

		self.text.bind("<Button-1>",	self.textFocus)
		self.text.bind("<Home>",	self.selectHome)
		self.text.bind("<End>",		self.selectEnd)
		self.text.bind("<Prior>",	self.selectPageUp)
		self.text.bind("<Next>",	self.selectPageDown)
		self.text.bind("<Up>",		self.selectUp)
		self.text.bind("<Down>",	self.selectDown)

		# Set commands only after initializing everything
		self.view.command = self.viewChanged
		self.group.command      = self.refresh
		self._bodies = []
		self._select = -1
		self._inrefresh = False
		self._hash   = -1
		self.refresh()

	# --------------------------------------------------------------------
	# Close button pressed
	# --------------------------------------------------------------------
	def close(self, event=None):
		self.update()
		GeometryAutoBody.width  = self.winfo_width()
		GeometryAutoBody.height = self.winfo_height()
		self.destroy()

	# --------------------------------------------------------------------
	# Help me please
	# --------------------------------------------------------------------
	def help(self, event=None):
		Manual.show("F4.11")

	# --------------------------------------------------------------------
	# Insert body to list
	# --------------------------------------------------------------------
	def insert(self, event=None):
		if self._select>=0:
			color, (tag, what, err) = self._bodies[self._select]
			what.insert(0, self.name.get())
			self.editor.addBody(Input.Card(tag, what))
			self.autoname()

	# --------------------------------------------------------------------
	def autoname(self):
		name = self.name.get()
		match = _NAMEPAT.match(name)
		if match:
			prefix = match.group(1)
		elif name:
			prefix = name
		else:
			prefix = "body"

		self.name.delete(0,END)
		self.name.insert(0, self.editor.findUniqueName(Input.Card.BODY,prefix))

	# --------------------------------------------------------------------
	# Projection Changed
	# --------------------------------------------------------------------
	def viewChanged(self, event=None):
		# Correct color
		view = self.view.get()
		if _VIEW.index(view) < 2:
			color = "Black"
		else:
			color = "Dark%s"%(view)
		self.view._text["foreground"] = color
		#self.refresh()

	# --------------------------------------------------------------------
	# @return a tuple with the uv view information
	# --------------------------------------------------------------------
	def _createUV(self, idx, xyz):
		# Convert to U,V (ignore W)
		viewer = self.master.viewers[idx]
		#matrix = bmath.Matrix(viewer.matrix(-1))
		matrix = bmath.Matrix(viewer.matrix())
		uv = [list(viewer.xyz2uv(*r)) for r in xyz]
		return (_VIEW[idx+2], uv, matrix)

	# --------------------------------------------------------------------
	# Refresh body list
	# --------------------------------------------------------------------
	def refresh(self, event=None):
		if self._inrefresh: return
		self._inrefresh = True
		self._refresh()
		self._inrefresh = False

	# --------------------------------------------------------------------
	# Force recalculation of bodies
	# --------------------------------------------------------------------
	def forceRefresh(self, event=None):
		self._hash = -1
		self.refresh()

	# --------------------------------------------------------------------
	def _refresh(self):
		selection = self.editor.getSelection(Input.Card.OBJECT)
		if len(selection)==0:
			self.npoints["text"] = _NONE
			return

		self.npoints["text"] = str(len(selection))
		group = self.group.get()

		#say("*-* Selection:",selection)
		#say("*-* Projection:",self.view.get())
		#say("*-* Group:",self.group.get())

		# List of points
		# Get X,Y,Z coordinates
		cards = []
		for tag in Input.FLAIR_TAGS:
			cards.extend([x for x in self.project.input[tag] if x.notIgnore()])
		cards.sort(Input.Card.cmpPos)

		xyz = []
		# Calculate hash to see if something has changed
		h = hash(self.view.get()) + hash(self.group.get())
		for card in selection:
			p = tuple(card.bodyP())
			xyz.append(p)
			h += hash(p[0]) + hash(p[1]) + hash(p[2])

		# -------------------------------
		# Check if something has changed
		if self._hash == h: return
		self._hash = h

		view = []	# a list of view to check
				# should be in the format (Name, xyz|uv, matrix|None)
		del self._bodies[:]	# Bodies found and categories
				# either: str: Category
				#	or (view-name, body-information)

		# Projection group index requested 0=all, 1=3d, 2,3,4,5 red,green,blue,magenta
		proj = _VIEW.index(self.view.get())

		# 3d view requested
		if proj in (0,1):
			view.append(("3D",xyz,None))
		# UV view requested:
		if proj==0:
			for i in range(4):
				view.append(self._createUV(i,xyz))
		elif proj>=2:
			view.append(self._createUV(proj-2, xyz))

		# ----------------
		# Check for planes
		# ----------------
		if group in ("all", "plane"):
			header = False
			for name, pts, matrix in view:
				if matrix is None:
					body = findBestPlane(findPlane(pts))
				else:
					body = findBestPlane(findPlaneUV(pts, matrix))

				if body is not None:
					if not header:
						self._bodies.append("Plane")
						header = True
					self._bodies.append((name,body))

		# ----------------
		# Check for sphere
		# ----------------
		if group in ("all", "sphere"):
			header = False
			for name, pts, matrix in view:
				if matrix is None:
					body = findSphere(pts)
				else:
					body = findSphereUV(pts, matrix)

				if body is not None:
					if not header:
						self._bodies.append("Sphere")
						header = True
					self._bodies.append((name,body))

		# ------------------
		# Check for cylinder
		# ------------------
		if group in ("all", "cylinder"):
			header = False
			for name, pts, matrix in view:
				if matrix is None:
					body = None
					#body = findCylinder(pts)
				else:
					body = findCylinderUV(pts, matrix)

				if body is not None:
					if not header:
						self._bodies.append("Cylinder")
						header = True
					self._bodies.append((name,body))

		self.fillBodies()

	# --------------------------------------------------------------------
	# Populate bodies textbox
	# --------------------------------------------------------------------
	def fillBodies(self):
		text = self.text
		text["state"] = NORMAL
		text.delete("0.0",END)
		text.update_idletasks()
		text.tag_delete("*")

		for idx,body in enumerate(self._bodies):
			if isinstance(body,str):
				text.insert(END, "%s\n"%(body), "title")
			else:
				name, (card, whats, err) = body
				if name == "3D":
					text.insert(END, name, name)
				else:
					text.insert(END, "2D", name)
				text.insert(END, " ")
				text.image_create(END, image=tkFlair.icons[card])
				text.insert(END, " ")

				tag = "B%d"%(idx)
				text.insert(END, "%s %s" % \
					(card, " ".join([bmath.format(x) for x in whats],"\t")),
					tag)
				text.tag_bind(tag, "<Button-1>", lambda e,s=self,i=idx : s.select(i))
				text.insert(END, "\t[%g]\n" % (err))

		if len(self._bodies)==0:
			self.select(-1)
			text.insert(END,"None found")
		else:
			text.tag_configure("title",   foreground="DarkBlue", font="Helvetica 14")
			text.tag_configure("3D",      background="#FFFF70")
			text.tag_configure("Red",     background="#FF7070")
			text.tag_configure("Green",   background="#70FF70")
			text.tag_configure("Blue",    background="#7070FF")
			text.tag_configure("Magenta", background="#FF70FF")
			self.select()

		text["state"] = DISABLED
		text.focus_set()

	# --------------------------------------------------------------------
	# Highlight and change select
	# --------------------------------------------------------------------
	def select(self, idx=None):
		if idx is not None:
			if self._select>=0:
				self.text.tag_configure("B%d"%(self._select),
					background="White")
			self._select = idx
		if self._select>=0:
			tag = "B%d"%(self._select)
			self.text.tag_configure(tag, background="LightBlue")
			#index = self.text.tag_nextrange(tag,"1.0")
			#if index: self.text.see(index[0])
			try: self.text.see("%s.first"%(tag))
			except TclError: pass

	# --------------------------------------------------------------------
	def selectHome(self, event=None):
		self.select(0)

	def selectEnd(self, event=None):
		self.select(len(self._bodies)-1)

	def selectUp(self, event=None):
		self.select(max(0,self._select-1))

	def selectDown(self, event=None):
		self.select(min(len(self._bodies)-1,self._select+1))

	selectPageUp   = selectUp
	selectPageDown = selectDown

	# ----------------------------------------------------------------------
	def textFocus(self, event):
		self.text.focus_set()

#===============================================================================
# General geometry transformation dialog
#===============================================================================
_TRANSFORM_PAGES = [ " Bodies ", "Repeat", "Concentric" ]
class TransformGeometryDialog(Toplevel):
	lastPage =  _TRANSFORM_PAGES[0]
	width    = -1
	height   = -1

	def __init__(self, master, flair, cardlist):
		Toplevel.__init__(self, master)
		self.transient(master)
		self.title("Transform geometry")

		self.flair    = flair
		self.input    = flair.project.input
		self.cardlist = cardlist

		# Tab page set
		self.tabPage = tkExtra.TabPageSet(self,
					pageNames=_TRANSFORM_PAGES,
					top=False)
		self.tabPage.pack(expand=YES, fill=BOTH)
		self.tabPage.changePage(TransformGeometryDialog.lastPage)

		# --------------------------------------------------------------
		# Bodies
		# --------------------------------------------------------------
		frame = self.tabPage[_TRANSFORM_PAGES[0]]
		self.bodies = TransformBodies(frame, flair, cardlist)
		self.bodies.pack(expand=YES, fill=BOTH)

		# --------------------------------------------------------------
		# Repeat
		# --------------------------------------------------------------
		frame = self.tabPage[_TRANSFORM_PAGES[1]]
		self.repeat = RepeatBodies(frame, flair, cardlist)
		self.repeat.pack(expand=YES, fill=BOTH)

		# --------------------------------------------------------------
		# Buttons
		# --------------------------------------------------------------
		frame = Frame(self)
		frame.pack(side=BOTTOM, fill=X)

		Button(frame, text="Close", command=self.close,
				image=tkFlair.icons["x"], compound=LEFT).pack(side=RIGHT)
		Button(frame, text="Apply", command=self.apply,
				image=tkFlair.icons["ok"], compound=LEFT,
				foreground="DarkRed", background="LightYellow",
				activeforeground="Red",
				activebackground="LightYellow").pack(side=RIGHT)

		# --------------------------------------------------------------
		# Bindings
		# --------------------------------------------------------------
		self.bind("<Escape>",		  self.close)
		self.protocol("WM_DELETE_WINDOW", self.close)

		# --- Show ---
		self.deiconify()
		if TransformGeometryDialog.width > 0:
			self.geometry("%dx%d" \
				%(TransformGeometryDialog.width, TransformGeometryDialog.height))
		self.lift()
		self.wait_visibility()
		self.grab_set()
		self.wait_window()

	# ----------------------------------------------------------------------
	def apply(self, event=None):
		page = self.tabPage.getActivePage()
		if page == _TRANSFORM_PAGES[0]:
			self.bodies.transform()
		elif page == _TRANSFORM_PAGES[1]:
			self.repeat.transform()
		else:
			pass

	# ----------------------------------------------------------------------
	def close(self, event=None):
		# remember last position and size
		TransformGeometryDialog.width    = self.winfo_width()
		TransformGeometryDialog.height   = self.winfo_height()
		TransformGeometryDialog.lastPage = self.tabPage.getActivePage()

		if TransformGeometryDialog.lastPage == _TRANSFORM_PAGES[0]:
			self.bodies.getTransformations()

		elif TransformGeometryDialog.lastPage == _TRANSFORM_PAGES[1]:
			self.repeat.get()

		self.destroy()

#===============================================================================
# Transform body dialog
#===============================================================================
_TRANSFORM_TYPES = [ "O", "T", "TX", "TY", "TZ", "RX" ,"RY", "RZ", "MX", "MY", "MZ", "S", "ignore" ]
class TransformBodies(Frame):
	transformList = [("RZ","90.0")]

	def __init__(self, master, flair, cardlist):
		Frame.__init__(self, master)
		self.flair    = flair
		self.input    = flair.project.input
		self.cardlist = [x for x in cardlist if x.type()==Input.Card.BODY]

		# --- ROT-DEFI ----
		frame = Frame(self)
		frame.pack(side=TOP, fill=X)
		l = Label(frame, text="ROT-DEFini:")
		l.pack(side=LEFT)
		self.combo = tkExtra.Combobox(frame,False,width=15)
		self.combo.pack(side=LEFT)
		tkExtra.Balloon.set(self.combo,
			"Import or create ROT-DEFini cards in input")
		self.combo._text.config(background="White")
		b = Button(frame, text="Add to Input", command=self.add2Input,
				pady=0, padx=3)
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b,
			"Convert transformation to ROT-DEFI cards in input")
		b = Button(frame, text="Get from Input", command=self.getFromInput,
				pady=0, padx=3)
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b,
			"Get transformation from input ROT-DEFI cards")

		# --- Transformations ----
		frame = Frame(self)
		frame.pack(side=TOP, expand=YES, fill=BOTH)

		self.listbox = tkExtra.MultiListbox(frame,
				(("Type",   10, None),
				 ("Value",  40, None)),
				 height=5)
		self.listbox.grid(row=0, rowspan=7, column=0, sticky=NSEW)
		for lst in self.listbox.lists: lst.config(takefocus=1)
		self.listbox.sortAssist = None

		b = Button(frame, image=tkFlair.icons["add"], takefocus=False,
				command=self.add)
		tkExtra.Balloon.set(b, "Add transformation")
		b.grid(row=0, column=1, sticky=EW)

		b = Button(frame, image=tkFlair.icons["del"], takefocus=False,
				command=self.delete)
		tkExtra.Balloon.set(b, "Delete transformation")
		b.grid(row=1, column=1, sticky=EW)

		b = Button(frame, image=tkFlair.icons["clone"], takefocus=False,
				command=self.clone)
		tkExtra.Balloon.set(b, "Clone Transformation")
		b.grid(row=2, column=1, sticky=EW)

		b = Button(frame, image=tkFlair.icons["up"], takefocus=False,
				command=self.moveUp)
		tkExtra.Balloon.set(b, "Move Up")
		b.grid(row=3, column=1, sticky=EW)

		b = Button(frame, image=tkFlair.icons["down"], takefocus=False,
				command=self.moveDown)
		tkExtra.Balloon.set(b, "Move Down")
		b.grid(row=4, column=1, sticky=EW)

		b = Button(frame, image=tkFlair.icons["toggle"], padx=0, takefocus=False,
				command=self.invert)
		tkExtra.Balloon.set(b, "Invert Transformation")
		b.grid(row=5, column=1, sticky=EW)

		frame.grid_columnconfigure(0, weight=1)
		frame.grid_rowconfigure(6, weight=1)

		# --- Global Variables ----
		frame = Frame(self)
		frame.pack(fill=X)

		Label(frame, text="Zero:").grid(row=0, column=0, sticky=W)
		self.zero = tkExtra.FloatEntry(frame, background="white")
		self.zero.grid(row=0, column=1, sticky=EW)
		tkExtra.Balloon.set(self.zero,
			"Threshold value to be considered zero or equal in comparison. ")

		Label(frame, text="Infinite:").grid(row=1, column=0, sticky=W)
		self.infinite = tkExtra.FloatEntry(frame, background="white")
		self.infinite.grid(row=1, column=1, sticky=EW)
		tkExtra.Balloon.set(self.infinite,
			"Length of cylinders when from infinite [XYZ]CC are converted to finite RCC")

		self.useQua = BooleanVar(self)
		self.useQua.set(Input._useQUA)

		b = Checkbutton(frame, text="Use QUA (instead of RCC/REC)",
				variable=self.useQua,
				onvalue=1, offvalue=0,
				anchor=W)
		b.grid(row=1, column=2, columnspan=2, sticky=W)
		tkExtra.Balloon.set(b,
			"If set infinite cylinders will be converted to QUA when rotated otherwise to RCC/REC")

		#frame.grid_columnconfigure(1, weight=1)

		# --- Bindings ---
		self.listbox.bindList('<ButtonRelease-1>', self.release)
		self.listbox.bindList('<Control-Key-Up>', self.moveUp)
		self.listbox.bindList('<Control-Key-Down>', self.moveDown)
		#self.listbox.bindList('<Double-1>', self.edit)
		self.listbox.bindList('<Return>',   self.release)
		self.listbox.bindList('<KP_Enter>', self.release)

		self.listbox.setPopupMenu([
				 ('Add',    0, self.add,    tkFlair.icons["add"]),
				 ('Delete', 0, self.delete, tkFlair.icons["x"]),
				 ('Clone',  0, self.clone,  tkFlair.icons["clone"]),
				 ('Invert', 0, self.invert, tkFlair.icons["invert"])] )
		self.listbox.bindLeftRight()
		self.listbox.focus_set()

		self.fill()
		self.zero.delete(0,END)
		self.infinite.delete(0,END)
		self.zero.insert(0, str(Input.zero))
		self.infinite.insert(0, str(Input.infinite))

	# ----------------------------------------------------------------------
	def fill(self):
		# Find all ROT-DEFI cards in input and fill combobox
		self.combo.fill(self.input.rotdefiList())

		# Fill transformations
		for entry in TransformBodies.transformList:
			self.listbox.insert(END, entry)

	# ----------------------------------------------------------------------
	def release(self, event=None):
		listbox = event.widget
		sel = self.listbox.curselection()
		if len(sel) != 1: return
		sel = int(sel[0])

		try: active = listbox.index(ACTIVE)
		except: return

		if active != sel: return

		lid = self.listbox.lists.index(listbox)

		if lid == 0:
			# Edit the transform type
			edit = tkExtra.InPlaceList(listbox,
					values=_TRANSFORM_TYPES,
					height=10)
			#if edit.value is None: return
		else:
			# Edit the value
			edit = tkExtra.InPlaceEdit(listbox)
			#if edit.value is None: return

	# ----------------------------------------------------------------------
	def add(self, event=None):
		self.listbox.insert(END, (_TRANSFORM_TYPES[0], "0.0 0.0 0.0"))
		self.listbox.selection_clear(0,END)
		self.listbox.selection_set(END)
		self.listbox.activate(END)
		self.listbox.see(END)
		self.listbox.focus_set()

	# ----------------------------------------------------------------------
	def delete(self, event=None):
		self.listbox.focus_set()
		lst = list(map(int, self.listbox.curselection()))
		if not lst: return
		lst.reverse()
		act = lst[0]
		for i in lst:
			self.listbox.delete(i)
		if act >= self.listbox.size(): act = self.listbox.size()-1
		if act >= 0:
			self.listbox.activate(act)
			self.listbox.selection_set(act)

	# ----------------------------------------------------------------------
	def clone(self, event=None):
		self.listbox.focus_set()
		sel = list(map(int,self.listbox.curselection()))
		if len(sel) == 0: return
		oldSize = self.listbox.size()
		for i in sel:
			entry = self.listbox.get(i)
			self.listbox.insert(END, entry)
		self.listbox.selection_clear(0,END)
		self.listbox.selection_set(oldSize,END)
		self.listbox.activate(END)

	# ----------------------------------------------------------------------
	def invert(self, event=None):
		self.getTransformations()
		self.listbox.delete(0, END)
		for type, value in TransformBodies.transformList:
			value = value.split()
			for i in range(len(value)):
				try:
					value[i] = str(-float(value[i]))
				except:
					value[i] = "0.0"

			self.listbox.insert(0, (type, " ".join(value)))

	# ----------------------------------------------------------------------
	def moveUp(self, event=None):
		for lst in self.listbox.lists:
			lst.moveUp()
		return "break"

	# ----------------------------------------------------------------------
	def moveDown(self, event=None):
		for lst in self.listbox.lists:
			lst.moveDown()
		return "break"

	# ----------------------------------------------------------------------
	def transform(self):
		# Threshold
		try:
			zero = abs(float(self.zero.get()))
			if 1e-16 < zero < 0.01:
				Input.zero = zero
		except ValueError:
			pass

		# Infinite
		try:
			infinite = float(self.infinite.get())
			if 1.0 < infinite < 1E10:
				Input.infinite = infinite
		except ValueError:
			pass

		# QUA or RCC/REC
		Input._useQUA = self.useQua.get()

		# Transformation
		self.getTransformations()
		matrix = bmath.Matrix(4, type=1)

		origin = None

		for trans, value in TransformBodies.transformList:
			if trans == "O":
				value = value.split()
				origin = bmath.Vector()
				for i in range(3):
					try:
						origin[i] = float(value[i])
					except:
						self.flair.notify("WARNING",
							"invalid origin vector",
							tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix.translate(-origin)
			elif trans == "T":
				value = value.split()
				t = bmath.Vector()
				for i in range(3):
					try:
						t[i] = float(value[i])
					except:
						self.flair.notify("WARNING",
							"invalid translate vector",
							tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix.translate(t)
			elif trans == "TX":
				t = bmath.Vector()
				try:
					t[0] = float(value)
				except:
					self.flair.notify("WARNING",
						"invalid translate-X value",
						tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix.translate(t)
			elif trans == "TY":
				t = bmath.Vector()
				try:
					t[1] = float(value)
				except:
					self.flair.notify("WARNING",
						"invalid translate-Y value",
						tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix.translate(t)
			elif trans == "TZ":
				t = bmath.Vector()
				try:
					t[2] = float(value)
				except:
					self.flair.notify("WARNING",
						"invalid translate-Z value",
						tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix.translate(t)
			elif trans == "RX":
				try:
					a = math.radians(float(value))
				except:
					a = 0.0
					self.flair.notify("WARNING",
						"invalid rotate-X value",
						tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix.rotX(-a)
			elif trans == "RY":
				try:
					a = math.radians(float(value))
				except:
					a = 0.0
					self.flair.notify("WARNING",
						"invalid rotate-Y value",
						tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix.rotY(-a)
			elif trans == "RZ":
				try:
					a = math.radians(float(value))
				except:
					a = 0.0
					self.flair.notify("WARNING",
						"invalid rotate-Z value",
						tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix.rotZ(-a)
			elif trans == "MX":
				if value!="":
					self.flair.notify("WARNING",
						"no value required for mirror-X",
						tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix(4,type=1)
				m[0][0] = -1.0
			elif trans == "MY":
				if value!="":
					self.flair.notify("WARNING",
						"no value required for mirror-Y",
						tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix(4,type=1)
				m[1][1] = -1.0
			elif trans == "MZ":
				if value!="":
					self.flair.notify("WARNING",
						"no value required for mirror-Z",
						tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix(4,type=1)
				m[2][2] = -1.0
			elif trans == "S":
				try:
					s = float(value)
				except:
					s = 1.0
					self.flair.notify("WARNING",
						"invalid Scale value",
						tkFlair.NOTIFY_WARNING)
				m = bmath.Matrix.scale(s)
			else:
				continue

			matrix = m * matrix

		if origin:
			matrix = bmath.Matrix.translate(origin) * matrix

		self.flair.log("\nTransformation matrix=")
		self.flair.log(str(matrix))
		undoinfo = [self.flair.refreshUndo()]
		for card in self.cardlist:
			undoinfo.append(self.flair.saveCardUndoInfo(card.pos()))
			self.flair.log("\n>>> Card=")
			self.flair.log(str(card))
			if self.input.transformBody(card, matrix):
				self.flair.log("*** Invalid card for transformation...")
			else:
				self.flair.log("<<< Card=")
				self.flair.log(str(card))
		undoinfo.append(self.flair.refreshUndo())
		self.flair.addUndo(undo.createListUndo(undoinfo,"Transform cards"))

		self.flair.setInputModified()
		self.flair.refresh("card")

	# ----------------------------------------------------------------------
	def getTransformations(self):
		# Remember transformations
		del TransformBodies.transformList[:]
		for i in range(self.listbox.size()):
			TransformBodies.transformList.append(self.listbox.get(i))

	# ----------------------------------------------------------------------
	def add2Input(self, event=None):
		# Get the name from the combo box
		name = self.combo.get().strip()
		hash = name.find('#')
		id   = 0
		if hash >= 0:
			try:
				id = int(name[hash+1:])
				name = ""
			except:
				pass
		else:
			try:
				id = int(name)
			except:
				id = 0
		if id==0 and name=="": return

		# Check if ROT-DEFI already exists
		undoinfo = [(self.flair.refreshUndo,True), self.flair.renumberUndo()]
		found = False
		for card in self.input["ROT-DEFI"]:
			w1 = card.intWhat(1)
			if w1 < 1000:
				(j, i) = divmod(w1, 100)
			else:
				(i, j) = divmod(w1, 1000)
			if (id!=0 and i==id) or (name!="" and name==card.sdum()):
				found = True
				break
		if found:
			ans = tkFlair.askyesno("ROT-DEFI %s already exists"%(name),
				"ROT-DEFI cards with id=\"%d\" name=\"%s\" already"
				"exists in input. Do you want to override?"%(id,name),
				parent=self)
			if not ans: return
			# Delete existing ROT-DEFI cards...
			for card in self.input["ROT-DEFI"]:
				w1 = card.intWhat(1)
				if w1 < 1000:
					(j, i) = divmod(w1, 100)
				else:
					(i, j) = divmod(w1, 1000)
				if (id!=0 and i==id) or (name!="" and name==card.sdum()):
					undoinfo.append(
						self.flair.delCardUndo(
							card._pos, False))

		# Create the necessary ROT-DEFI cards
		pos = self.input.bestPosition("ROT-DEFI")
		self.getTransformations()
		for trans, value in TransformBodies.transformList:
			what = [""]*7
			what[0] = name
			if trans == "T":
				axis  = 0
				value = value.split()
				for i in range(3):
					try:
						what[4+i] = float(value[i])
					except:
						self.flair.notify("WARNING"
							"Invalid T vector",
							tkFlair.NOTIFY_WARNING)
						break
			elif trans in ("TX", "TY", "TZ"):
				axis = ord(trans[1])-ord('X')+1
				try:
					what[3+axis] = float(value)
				except:
					self.flair.notify("WARNING",
						"Invalid %s value"%(trans),
						tkFlair.NOTIFY_WARNING)
					break
			elif trans in ("RX", "RY", "RZ"):
				axis = ord(trans[1])-ord('X')+1
				try:
					what[3] = float(value)
				except:
					self.flair.notify("WARNING",
						"Invalid %s value"%(trans),
						tkFlair.NOTIFY_WARNING)
					break
			else:
				self.flair.notify("WARNING",
					"No Scaling is allowed in ROT-DEFI",
					tkFlair.NOTIFY_WARNING)
				break
			if id < 100:
				what[1] = axis*100 + id
			else:
				what[1] = id*1000 + axis

			card = Input.Card("ROT-DEFI", what,
				"Added by flair Transformation Dialog")
			undoinfo.append(self.flair.addCardUndo(card, pos, False))
			pos += 1

		undoinfo.append(self.flair.renumberUndo())
		undoinfo.append(undoinfo[0])
		self.flair.addUndo(undo.createListUndo(undoinfo,"Add ROT-DEFI cards"))
		self.combo.fill(self.input.rotdefiList())
		self.combo._text.delete(0,END)
		self.combo._text.insert(0,name)
		self.flair.redraw()

	# ----------------------------------------------------------------------
	def getFromInput(self, event=None):
		name = self.combo.get().strip()
		hash = name.find('#')
		id   = 0
		if hash >= 0:
			try:
				id = int(name[hash+1:])
				name = ""
			except:
				pass
		else:
			try:
				id = int(name)
			except:
				id = 0
		if id==0 and name=="": return

		self.listbox.delete(0, END)
		for card in self.input["ROT-DEFI"]:
			if card.ignore(): continue
			w1 = card.intWhat(1)
			if w1 < 1000:
				(j, i) = divmod(w1, 100)
			else:
				(i, j) = divmod(w1, 1000)
			if i==id and card.sdum()==name:
				if j==0: j=3
				theta = card.numWhat(2)
				phi   = card.numWhat(3)
				x = card.numWhat(4)
				y = card.numWhat(5)
				z = card.numWhat(6)
				if x!=0.0 or y!=0.0 or z!=0.0:
					self.listbox.insert(END,
						('T',"%g %g %g"%(x,y,z)))
				if phi!=0.0:
					self.listbox.insert(END,
						('R%s'%("XYZ"[j-1]), str(phi)))
				if theta!=0.0:
					self.listbox.insert(END,
						('R%s'%("ZXY"[j-1]), str(theta)))
		self.listbox.selection_set(0)
		self.listbox.activate(0)

#===============================================================================
# Repeat bodies
#===============================================================================
class RepeatBodies(Frame):
	Nu = 5
	Nv = 1
	Nw = 1
	U  = bmath.Vector.X
	V  = bmath.Vector.Y
	W  = bmath.Vector.Z

	def __init__(self, master, flair, cardlist):
		Frame.__init__(self, master)
		self.flair    = flair
		self.input    = flair.project.input
		self.cardlist = [x for x in cardlist if x.type()==Input.Card.BODY]

		self.nu = IntVar()
		self.nv = IntVar()
		self.nw = IntVar()

		self.nu.set(RepeatBodies.Nu)
		self.nv.set(RepeatBodies.Nv)
		self.nw.set(RepeatBodies.Nw)

		# ---- Labels ----
		row = 0
		Label(self, text="Axis").grid(row=row, column=0)
		Label(self, text="N (#)").grid(row=row, column=1)
		Label(self, text="x").grid(row=row, column=2)
		Label(self, text="y").grid(row=row, column=3)
		Label(self, text="z").grid(row=row, column=4)

		# ---------- U ---------
		row += 1
		Label(self, text="U:", foreground="DarkRed").grid(row=row, column=0, sticky=E)

		nu = Spinbox(self, text=self.nu, background="White", from_=1, to=10000, width=4)
		nu.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(nu, "Number of replicas along U axis")

		self.ux = tkExtra.FloatEntry(self, background="White", width=8)
		self.ux.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(self.ux, "Ux direction")

		self.uy = tkExtra.FloatEntry(self, background="White", width=8)
		self.uy.grid(row=row, column=3, sticky=EW)
		tkExtra.Balloon.set(self.uy, "Uy direction")

		self.uz = tkExtra.FloatEntry(self, background="White", width=8)
		self.uz.grid(row=row, column=4, sticky=EW)
		tkExtra.Balloon.set(self.uz, "Uz direction")

		Button(self, text="x", command=self.setUx, takefocus=False).grid(row=row, column=5)
		Button(self, text="y", command=self.setUy, takefocus=False).grid(row=row, column=6)
		Button(self, text="z", command=self.setUz, takefocus=False).grid(row=row, column=7)

		# ---------- V ---------
		row += 1
		Label(self, text="V:", foreground="DarkGreen").grid(row=row, column=0, sticky=E)

		nv = Spinbox(self, text=self.nv, background="White", from_=1, to=10000, width=4)
		nv.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(nv, "Number of replicas along V axis")

		self.vx = tkExtra.FloatEntry(self, background="White", width=8)
		self.vx.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(self.vx, "Vx direction")

		self.vy = tkExtra.FloatEntry(self, background="White", width=8)
		self.vy.grid(row=row, column=3, sticky=EW)
		tkExtra.Balloon.set(self.vy, "Vy direction")

		self.vz = tkExtra.FloatEntry(self, background="White", width=8)
		self.vz.grid(row=row, column=4, sticky=EW)
		tkExtra.Balloon.set(self.vz, "Vz direction")

		Button(self, text="x", command=self.setVx, takefocus=False).grid(row=row, column=5)
		Button(self, text="y", command=self.setVy, takefocus=False).grid(row=row, column=6)
		Button(self, text="z", command=self.setVz, takefocus=False).grid(row=row, column=7)

		# ---------- W ---------
		row += 1
		Label(self, text="W:", foreground="DarkBlue").grid(row=row, column=0, sticky=E)

		nw = Spinbox(self, text=self.nw, background="White", from_=1, to=10000, width=4)
		nw.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(nw, "Number of replicas along W axis")

		self.wx = tkExtra.FloatEntry(self, background="White", width=8)
		self.wx.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(self.wx, "Wx direction")

		self.wy = tkExtra.FloatEntry(self, background="White", width=8)
		self.wy.grid(row=row, column=3, sticky=EW)
		tkExtra.Balloon.set(self.wy, "Wy direction")

		self.wz = tkExtra.FloatEntry(self, background="White", width=8)
		self.wz.grid(row=row, column=4, sticky=EW)
		tkExtra.Balloon.set(self.wz, "Wz direction")

		Button(self, text="x", command=self.setWx, takefocus=False).grid(row=row, column=5)
		Button(self, text="y", command=self.setWy, takefocus=False).grid(row=row, column=6)
		Button(self, text="z", command=self.setWz, takefocus=False).grid(row=row, column=7)

		self.grid_columnconfigure(1, weight=1)
		self.grid_columnconfigure(2, weight=2)
		self.grid_columnconfigure(3, weight=2)
		self.grid_columnconfigure(4, weight=2)

		# ---- Fill values ----
		self.ux.insert(0,str(RepeatBodies.U[0]))
		self.uy.insert(0,str(RepeatBodies.U[1]))
		self.uz.insert(0,str(RepeatBodies.U[2]))

		self.vx.insert(0,str(RepeatBodies.V[0]))
		self.vy.insert(0,str(RepeatBodies.V[1]))
		self.vz.insert(0,str(RepeatBodies.V[2]))

		self.wx.insert(0,str(RepeatBodies.W[0]))
		self.wy.insert(0,str(RepeatBodies.W[1]))
		self.wz.insert(0,str(RepeatBodies.W[2]))

	# ----------------------------------------------------------------------
	def setU(self, x,y,z):
		self.ux.delete(0,END)
		self.uy.delete(0,END)
		self.uz.delete(0,END)

		self.ux.insert(0, str(x))
		self.uy.insert(0, str(y))
		self.uz.insert(0, str(z))

	# ----------------------------------------------------------------------
	def setV(self, x,y,z):
		self.vx.delete(0,END)
		self.vy.delete(0,END)
		self.vz.delete(0,END)

		self.vx.insert(0, str(x))
		self.vy.insert(0, str(y))
		self.vz.insert(0, str(z))

	# ----------------------------------------------------------------------
	def setW(self, x,y,z):
		self.wx.delete(0,END)
		self.wy.delete(0,END)
		self.wz.delete(0,END)

		self.wx.insert(0, str(x))
		self.wy.insert(0, str(y))
		self.wz.insert(0, str(z))

	# ----------------------------------------------------------------------
	def setUx(self):	self.setU(1.0,0.0,0.0)
	def setUy(self):	self.setU(0.0,1.0,0.0)
	def setUz(self):	self.setU(0.0,0.0,1.0)

	def setVx(self):	self.setV(1.0,0.0,0.0)
	def setVy(self):	self.setV(0.0,1.0,0.0)
	def setVz(self):	self.setV(0.0,0.0,1.0)

	def setWx(self):	self.setW(1.0,0.0,0.0)
	def setWy(self):	self.setW(0.0,1.0,0.0)
	def setWz(self):	self.setW(0.0,0.0,1.0)

	# ----------------------------------------------------------------------
	def get(self):
		RepeatBodies.Nu = int(self.nu.get())
		RepeatBodies.Nv = int(self.nv.get())
		RepeatBodies.Nw = int(self.nw.get())

		RepeatBodies.U[0] = float(self.ux.get())
		RepeatBodies.U[1] = float(self.uy.get())
		RepeatBodies.U[2] = float(self.uz.get())

		RepeatBodies.V[0] = float(self.vx.get())
		RepeatBodies.V[1] = float(self.vy.get())
		RepeatBodies.V[2] = float(self.vz.get())

		RepeatBodies.W[0] = float(self.wx.get())
		RepeatBodies.W[1] = float(self.wy.get())
		RepeatBodies.W[2] = float(self.wz.get())

	# ----------------------------------------------------------------------
	def transform(self):
		self.get()

		bodynames = {}
		for card in self.input.cardsCache("bodies"):
			bodynames[card] = True

		pos = self.cardlist[-1].pos()+1

		undoinfo = [self.flair.refreshUndo()]
		for w in range(RepeatBodies.Nw):
			v1 = w*bmath.Vector(RepeatBodies.W)
			for v in range(RepeatBodies.Nv):
				v2 = v1 + v*bmath.Vector(RepeatBodies.V)
				for u in range(RepeatBodies.Nu):
					if u==0 and v==0 and w==0: continue
					vector = v2 + u*bmath.Vector(RepeatBodies.U)
					matrix = bmath.Matrix.translate(vector)

					for card in self.cardlist:
						clone = card.clone()
						# Find unique name
						# find number in the end of the name
						pat = Project._LAST_NUM.match(clone.sdum())
						if pat:
							name = pat.group(1)
							n = int(pat.group(2))+1
						else:
							name = clone.sdum()
							n = 1
						while True:
							sn = str(n)
							guess = "%s%s"%(name[:8-len(sn)],sn)
							if guess not in bodynames:
								clone.setSdum(guess)
								break
							n += 1
						bodynames[clone.sdum()] = True

						undoinfo.append(self.flair.addCardUndo(clone, pos, True))
						self.input.transformBody(clone, matrix)
						pos += 1

		undoinfo.append(self.flair.refreshUndo())
		self.flair.addUndo(undo.createListUndo(undoinfo,"Repeat bodies cards"))

		self.flair.setInputModified()
		self.flair.refresh("card")
