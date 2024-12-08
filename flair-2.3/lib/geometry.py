#!/usr/bin/python
# -*- coding: latin1 -*-
#
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	14-Aug-2009

import math
import bmath

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

#===============================================================================
# Plane
#===============================================================================
class Plane:
	"""Plane class
		defined by equation a*x + b*y + c*z + d = 0
	"""
	def __init__(self):
		self.a = 0.0
		self.b = 0.0
		self.c = 1.0
		self.d = 0.0

	# ----------------------------------------------------------------------
	# Fit plane to a set of points, using least square minimization
	# Points is a list of lists[x,y,z]
	# Distance of a points from a plane is
	#	di = a*xi + b*yi + c*zi + d (or to be easier -d)
	# We want to minimize the function
	#	f(a,b,c,d) = Sum(di^2) = Sum[(a*xi+b*yi+c*zi-d)^2]
	# Therefore we have to solve the system
	#	df/da = Sum[(a*xi+b*yi+c*zi-d)*xi] = 0
	#	df/db = Sum[(a*xi+b*yi+c*zi-d)*yi] = 0
	#	df/dc = Sum[(a*xi+b*yi+c*zi-d)*zi] = 0
	#	df/dd = Sum[ a*xi+b*yi+c*zi-d] = 0
	# From last equation we have
	#	df/dd = a*Sum(xi) + b*Sum(yi) + c*Sum(zi) - d = 0
	# if we divide by n
	#	a*xc + b*yc + c*zc -d 0
	# where xc,yc,zc = Sum(xi)/n, Sum(yi)/n, Sum(zi)/n	the centroid
	# The plane is passing from the centroid
	# Therefore we can modify the equation of the plane to
	#	df/da = Sum[(a*(xi-xc) + b*(yi-yc) +c*(zi-zc))*(xi-xc)] = 0
	#	df/db = Sum[(a*(xi-xc) + b*(yi-yc) +c*(zi-zc))*(yi-yc)] = 0
	#	df/dc = Sum[(a*(xi-xc) + b*(yi-yc) +c*(zi-zc))*(zi-zc)] = 0
	# which is equivalent to the over-determined linear system
	#	a00 = Sum(xi-xc)^2
	#	a11 = Sum(yi-yc)^2
	#	a22 = Sum(zi-zc)^2
	#	a01 = a10 = Sum[(xi-xc)*(yi-yc)]
	#	a02 = a20 = Sum[(xi-xc)*(zi-zc)]
	#	a12 = a21 = Sum[(yi-yc)*(zi-zc)]
	# linear system
	#	   / a00  a01  a02 \ / a \
	#	   | a10  a11  a12 | | b | = 0
	#	   \ a20  a21  a22 / \ c /
	# The system has a trivial solution a=b=c=0.
	# For the non-trivial we may require a^2+b^2+c^2 = 1
	# The eigenvectors of the matrix are mutually orthogonal and represent
	# the best, intermediate and worst plane.
	# we will choose the best with the minimum eigenvalue(=sum of square distances)
	# ----------------------------------------------------------------------
	def fromPoints(self, points):
		if len(points) < 3: return True

		# assume a!=0.0, set a=1
		A = []
		B = []
		found = False
		for x,y,z in points:
			A.append([y,z,1.0])
			B.append([-x])
		try:
			X = bmath.solveOverDetermined(bmath.Matrix(A),
					bmath.Matrix(B))
			self.a = 1.0
			self.b = X[0]
			self.c = X[1]
			self.d = X[2]
			self.normalize()
			# Find rms distance
			adist = 0.0
			asolution = (self.a, self.b, self.c, self.d)
			for x,y,z in points:
				adist += abs(self.a*x+self.b*y+self.c*z+self.d)
			found = True
		except:
			adist = 1E30

		# assume b!=0.0, set b=1
		A = []
		B = []
		for x,y,z in points:
			A.append([x,z,1.0])
			B.append([-y])
		try:
			X = bmath.solveOverDetermined(bmath.Matrix(A),
					bmath.Matrix(B))
			self.a = X[0]
			self.b = 1.0
			self.c = X[1]
			self.d = X[2]
			self.normalize()
			# Find rms distance
			bdist = 0.0
			bsolution = (self.a, self.b, self.c, self.d)
			for x,y,z in points:
				bdist += abs(self.a*x+self.b*y+self.c*z+self.d)
			found = True
		except:
			bdist = 1E30

		# assume c!=0.0, set c=1
		A = []
		B = []
		for x,y,z in points:
			A.append([x,y,1.0])
			B.append([-z])
		try:
			X = bmath.solveOverDetermined(bmath.Matrix(A),
					bmath.Matrix(B))
			self.a = X[0]
			self.b = X[1]
			self.c = 1.0
			self.d = X[2]
			self.normalize()
			# Find rms distance
			cdist = 0.0
			csolution = (self.a, self.b, self.c, self.d)
			for x,y,z in points:
				cdist += abs(self.a*x+self.b*y+self.c*z+self.d)
			found = True
		except:
			cdist = 1E30

		if not found: return True

		if   adist < bdist and adist < cdist:
			self.a, self.b, self.c, self.d = asolution
		elif bdist < adist and bdist < cdist:
			self.a, self.b, self.c, self.d = bsolution
		else:
			self.a, self.b, self.c, self.d = csolution
		return False

#		# plane should pass from the centroid
#		xc, yc, zc = 0.0, 0.0, 0.0
#		for x,y,z in points:
#			xc += x
#			yc += y
#			zc += z
#		n = float(len(points))
#		xc /= n
#		yc /= n
#		xc /= n
#		say(xc,yc,zc)
#
#		# find coefficients
#		a00, a11, a22 = 0.0, 0.0, 0.0
#		a01, a02, a12 = 0.0, 0.0, 0.0
#		for x,y,z in points:
#			X = x-xc
#			Y = y-yc
#			Z = z-zc
#			a00 += X*X
#			a11 += Y*Y
#			a22 += Z*Z
#			a01 += X*Y
#			a02 += X*Z
#			a12 += Y*Z
#
#		A = bmath.Matrix(3)
#		A[0] = [a00, a01, a02]
#		A[1] = [a01, a11, a12]
#		A[2] = [a02, a12, a22]
#		v,u = bmath.eigenvalues(A)
#		say("Eigenvalues=",v)
#		say("Eigenvectors=\n",u)
#		v = map(abs, v)
#		if   v[0]<v[1] and v[0]<v[2]: best = 0
#		elif v[1]<v[0] and v[1]<v[2]: best = 1
#		else:                         best = 2
#		ex = bmath.Vector(u[0])
#		ey = bmath.Vector(u[1])
#		ez = bmath.Vector(u[2])
#		say("ex.cross.ey=",ex.cross(ey))
#		say("ez=",ez)
#
#		self.a = u[best][0]
#		self.b = u[best][1]
#		self.c = u[best][2]
#		self.d = -(self.a*xc + self.b*yc + self.c*zc)
#		say(self)
#		return False

	# ----------------------------------------------------------------------
	def normal(self):
		"""return normal vector to plane"""
		n = bmath.Vector(self.a, self.b, self.c)
		n.norm()
		return n

	# ----------------------------------------------------------------------
	def normalize(self):
		"""normalize a^2+b^2+c^2=1"""
		f = math.sqrt(self.a**2 + self.b**2 + self.c**2)
		self.a /= f
		self.b /= f
		self.c /= f
		self.d /= f

	# ----------------------------------------------------------------------
	def distance(self, p):
		"""distance of vector p from plane"""
		return (self.a*p[0] + self.b*p[1] + self.c*p[2] + self.d) / \
			math.sqrt(self.a**2 + self.b**2 + self.c**2)

	# ----------------------------------------------------------------------
	def __str__(self):
		def fmt(s,a,b):
			if abs(a) < 1e-10:
				return s
			if abs(a-1.0) < 1e-10:
				if s!="":
					if b!="":
						return "%s+%s"%(s,b)
					else:
						return "%s+1"%(s)
				else:
					return b
			elif abs(a+1.0) < 1e-10:
				if b!="":
					return "%s-%s"%(s,b)
				else:
					return "%s-1"%(s)
			elif a > 0.0:
				if s!="":
					return "%s+%g%s"%(s,a,b)
				else:
					return "%g%s"%(a,b)
			else:
				return "%s%g%s"%(s,a,b)

		s = fmt("", self.a, "x")
		s = fmt(s, self.b, "y")
		s = fmt(s, self.c, "z")
		s = fmt(s, self.d, "")
		s += "=0"
		return s

# ------------------------------------------------------------------------------
def fitplane(points, tol=1.0E-5):
	"""
	Check if a set of points are coplanar within a certain tolerance tol.
	return the plane of the points otherwise None
	"""
	plane = Plane()
	if plane.fromPoints(points): return None

	# check distance of points from the plane and compare with
	# the tolerance
	for p in points:
		d = plane.distance(p)
		if abs(d)>tol: return None
	return plane

# ------------------------------------------------------------------------------
if __name__ == "__main__":
	import random
	from log import say

	random.seed(0.5)
	points = []
	u = bmath.Vector(1.0, 2.0,  3.0)
	v = bmath.Vector(-10.0, 1.0, 1.0)
	p = bmath.Vector(1.0, 2.5,  3.0)
	for i in range(10):
		x = random.random()
		y = random.random()
		z = random.random()
		#points.append([x,y,z])
		points.append(x*u + y*v + p)
		say(points[-1])
	plane = fitplane(points)
	say("Coplanar plane=",plane)

	# find transformation matrix so n=z' and pass from zero
	#plane.normalize()
	zaxis = plane.normal()
	p = plane.d * zaxis
	say("p=",p)
	# try xaxis
	say("zaxis=",zaxis)
	yaxis = zaxis.cross(bmath.Vector.X)
	say("yaxis=",yaxis)
	if yaxis.length2() < 1e-5:
		yaxis = zaxis.cross(bmath.Vector.Y)
		say("	yaxis=",yaxis)
	yaxis.norm()
	xaxis = yaxis.cross(zaxis)
	xaxis.norm()
	say("xaxis=",xaxis)
	R = bmath.Matrix()
	R.make(xaxis, yaxis, zaxis)
	T = bmath.Matrix.translate(p)
	M = R*T
	say("R:\n",R)
	say("T:\n",T)
	say("M:\n",M)
	points2 = []
	for x,y,z in points:
		points2.append(M * bmath.Vector(x,y,z))
		say(points2[-1])
	plane2 = fitplane(points2)
	say(plane2)
