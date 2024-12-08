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
# Date:	04-Feb-2010

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import sys
import math
import string
from math import *
from bmath import *

# Conic types
CONIC_NOTHING   = 0
#CONIC_POINT    = 1
CONIC_LINE      = 2
CONIC_LINES     = 3
CONIC_ELLIPSE   = 4
CONIC_PARABOLA  = 5
CONIC_HYPERBOLA = 6
CONIC_TYPES = ( "NOTHING", "POINT", "LINE", "LINES", "ELLIPSE", "PARABOLA", "HYPERBOLA" )

TOOSMALL  = 1.0e-15
SMALL     = 1.0e-10
CONICSTEP = 0.1
INFINITE  = 1.0e+10

# Processing of regions
_REGTRANS = string.maketrans("+-|()","     ")

#-------------------------------------------------------------------------------
def _Sx(a,suffix=""):
	if a < 0.0:
		if a==-1.0:
			if suffix!="":
				return "-%s"%(suffix)
			else:
				return "-1"
		else:
			if suffix!="":
				return "-%g*%s"%(-a,suffix)
			else:
				return "-%g"%(-a)
	elif a > 0.0:
		if a==1.0:
			if suffix!="":
				return "+%s"%(suffix)
			else:
				return "+1"
		else:
			if suffix!="":
				return "+%g*%s"%(a,suffix)
			else:
				return "+%g"%(a)
	else:
		return ""

#-------------------------------------------------------------------------------
def _cmpV(a,b):
	return cmp(a[0], b[0])

#===============================================================================
# The following parametric classes convert a conic to a parametric equation
# FIXME The order is important, in the sense that increasing the 't' to the
#       left has to be considered inside (-), and the right side outside (+)
# OR NOT?
#===============================================================================
class ConicParametric:
	def __init__(self, conic):
		self.conic = conic

	# --------------------------------------------------------------------
	def toStr(self, idx):
		return ""

	# --------------------------------------------------------------------
	def __str__(self):
		return self.toStr()

#===============================================================================
# Convert a conic to a parametric equation of a line
#	x = c1 + c2*t
#	y = c3 + c4*t
#===============================================================================
class LineParametric(ConicParametric):
	def __init__(self, conic):
		ConicParametric.__init__(self, conic)

		# Find starting point
		if abs(conic.g) < abs(conic.f):
			self.c1 = 0.0
			self.c3 = -conic.c / (2.0*conic.f)
		else:
			self.c1 = -conic.c / (2.0*conic.g)
			self.c3 = 0.0

		self.c2 = -conic.f
		self.c4 =  conic.g

	# --------------------------------------------------------------------
	# Evaluate a position t
	# --------------------------------------------------------------------
	def __call__(self, t):
		"""Return (x,y) position of t"""
		return (self.c1 + self.c2*t, self.c3 + self.c4*t)

	# --------------------------------------------------------------------
	# Return t from a point
	# --------------------------------------------------------------------
	def getT(self, x, y):
		"""Return t value from a point (x,y) beloning to the conic"""
		if abs(self.c2) < abs(self.c4):
			return (y-self.c3) / self.c4
		else:
			return (x-self.c1) / self.c2

	# --------------------------------------------------------------------
	def toStr(self, idx=""):
		s = "x%s(t)="%(idx)
		s += _Sx(self.c1)
		s += _Sx(self.c2,"t")
		s += "; y%s(t)="%(idx)
		s += _Sx(self.c3)
		s += _Sx(self.c4,"t")
		return s

#===============================================================================
# Common base class of Ellipse and Hyperbola parametric
#===============================================================================
class _EllHyperParametric(ConicParametric):
	def __init__(self, conic):
		ConicParametric.__init__(self, conic)

		# convert conic to standard one
		theta = conic.theta()
		if abs(theta)<TOOSMALL:
			ct = 1.0
			st = 0.0
			a = conic.a
			b = conic.b
			g = conic.g
			f = conic.f
		else:
			ct = cos(theta)
			st = sin(theta)
			stct = st*ct
			a = conic.a*ct**2 + 2.0*conic.h*stct + conic.b*st**2
			b = conic.a*st**2 - 2.0*conic.h*stct + conic.b*ct**2
			g =  conic.g*ct + conic.f*st
			f = -conic.g*st + conic.f*ct

		# For hyperbola consider two cases
		inverse = (conic.type == CONIC_HYPERBOLA and conic.D<0.0)

		dx = -g/a
		dy = -f/b
		c  = conic.c - g**2/a - f**2/b

		a  = 1.0 / sqrt(abs(a/c))
		b  = 1.0 / sqrt(abs(b/c))

		self.c1 =  dx*ct - dy*st
		self.c4 =  dx*st + dy*ct

		if inverse:
			self.c3 =  a*ct
			self.c2 = -b*st
			self.c6 =  a*st
			self.c5 =  b*ct
		else:
			self.c2 =  a*ct
			self.c3 = -b*st
			self.c5 =  a*st
			self.c6 =  b*ct

#===============================================================================
# Convert a conic to a parametric equation of an ellipse
#	x = c1 + c2*cos(t) + c3*sin(t)
#	y = c4 + c5*cos(t) + c6*sin(t)
#===============================================================================
class EllipseParametric(_EllHyperParametric):
	# --------------------------------------------------------------------
	# Evaluate a position t
	# --------------------------------------------------------------------
	def __call__(self, t):
		"""Return (x,y) position of t"""
		ct = cos(t)
		st = sin(t)
		return (self.c1 + self.c2*ct + self.c3*st, \
			self.c4 + self.c5*ct + self.c6*st)

	# --------------------------------------------------------------------
	# Return t from point
	# --------------------------------------------------------------------
	def getT(self, x, y):
		"""Return t value from a point (x,y) beloning to the conic"""
		x -= self.c1
		y -= self.c4

		# Assume z=cos(t), sin(t)=+/-sqrt(1-z^2)
		A = self.c2**2 + self.c3**2
		B = x*self.c2	# 2.0*...
		C = x**2 - self.c3**2
		D = B**2 - A*C

		if D > SMALL:
			D = sqrt(D)
			cost1 = (B + D) / A
			try:
				sint1 = sqrt((1.0-cost1)*(1.0+cost1))
			except ValueError:
				cost1 = 1.0
				sint1 = 0.0

			# Evaluate the parametric equation and compare against y
			dx1 = abs(self.c2*cost1 + self.c3*sint1 - x)
			dy1 = abs(self.c5*cost1 + self.c6*sint1 - y)
			if dx1<SMALL and dy1<SMALL: return acos(cost1)

			dx2 = abs(self.c2*cost1 - self.c3*sint1 - x)
			dy2 = abs(self.c5*cost1 - self.c6*sint1 - y)
			if dx2<SMALL and dy2<SMALL: return -acos(cost1)

			cost2 = (B - D) / A
			try:
				sint2 = sqrt((1.0-cost2)*(1.0+cost2))
			except ValueError:
				cost2 = 1.0
				sint2 = 0.0

			dx3 = abs(self.c2*cost2 + self.c3*sint2 - x)
			dy3 = abs(self.c5*cost2 + self.c6*sint2 - y)
			if dx3<SMALL and dy3<SMALL: return acos(cost2)

			dx4 = abs(self.c2*cost2 - self.c3*sint2 - x)
			dy4 = abs(self.c5*cost2 - self.c6*sint2 - y)
			if dx4<SMALL and dy4<SMALL: return -acos(cost2)

			# find the smallest distance
			a1 = dx1**2 + dy1**2
			a2 = dx2**2 + dy2**2
			a3 = dx3**2 + dy3**2
			a4 = dx4**2 + dy4**2
			# else return the smallest one
			if a1<=a2 and a1<=a3 and a1<=a4: return  acos(cost1)
			if a2<=a1 and a2<=a3 and a2<=a4: return -acos(cost1)
			if a3<=a1 and a3<=a2 and a3<=a4: return  acos(cost2)
			return -acos(cost2)

		else:
			#if abs(D)<SMALL:
			cost1 = B / A
			try:
				sint1 = sqrt((1.0-cost1)*(1.0+cost1))
			except:
				cost1 = 1.0
				sint1 = 0.0

			dx1 = abs(self.c2*cost1 + self.c3*sint1 - x)
			dy1 = abs(self.c5*cost1 + self.c6*sint1 - y)
			if dx1<SMALL and dy1<SMALL: return acos(cost1)

			dx2 = abs(self.c2*cost1 - self.c3*sint1 - x)
			dy2 = abs(self.c5*cost1 - self.c6*sint1 - y)

			if dx1**2 + dy1**2 < dx2**2 + dy2**2:
				return acos(cost1)
			else:
				return -acos(cost1)
		return None

	# --------------------------------------------------------------------
	def toStr(self, idx=""):
		s = "x%s(t)="%(idx)
		s += _Sx(self.c1)
		s += _Sx(self.c2,"cos(t)")
		s += _Sx(self.c3,"sin(t)")
		s += "; y%s(t)="%(idx)
		s += _Sx(self.c4)
		s += _Sx(self.c5,"cos(t)")
		s += _Sx(self.c6,"sin(t)")
		return s

#===============================================================================
# Convert a conic to a parametric equation of a hyperbola
#	x = c1 + c2*sec(t) + c3*tan(t)		[sect = 1/cost]
#	y = c4 + c5*sec(t) + c6*tan(t)
#
# or with hyperbolic functions
#	x = c1 + c2*cosh(t) + c3*sinh(t)
#	y = c4 + c5*cosh(t) + c6*sinh(t)
# however they have the drawback that they only describe the positive part(right)
# and a second instance is needed for the left part
#===============================================================================
class HyperbolaParametric(_EllHyperParametric):
	# --------------------------------------------------------------------
	# Evaluate a position t
	# --------------------------------------------------------------------
	def __call__(self, t):
		"""Return (x,y) position of t"""
		ct = 1.0/cos(t)
		st = tan(t)
		return (self.c1 + self.c2*ct + self.c3*st, \
			self.c4 + self.c5*ct + self.c6*st)

	# --------------------------------------------------------------------
	# Return t from point
	# --------------------------------------------------------------------
	def getT(self, x, y):
		"""Return t value from a point (x,y) beloning to the conic"""
		x -= self.c1
		y -= self.c4

		# Assume z=sec(t), tan(t)=+/-sqrt(z^2-1)
		A = (self.c2 - self.c3)*(self.c2 + self.c3)
		B = x*self.c2	# 2.0*...
		C = x**2 + self.c3**2
		D = B**2 - A*C

		if D > SMALL:
			D = sqrt(D)
			sect1 = (B + D) / A
			if sect1 > 0.0:
				sign1 = 1.0
			else:
				sign1 = -1.0
			try: tant1 = sqrt((sect1-1.0)*(sect1+1.0))
			except:
				sect1 = sign1
				tant1 = 0.0

			# Evaluate the parametric equation and compare against y
			dx1 = abs(self.c2*sect1 + self.c3*tant1 - x)
			dy1 = abs(self.c5*sect1 + self.c6*tant1 - y)
			if dx1<SMALL and dy1<SMALL: return sign1*acos(1.0/sect1)

			dx2 = abs(self.c2*sect1 - self.c3*tant1 - x)
			dy2 = abs(self.c5*sect1 - self.c6*tant1 - y)
			if dx2<SMALL and dy2<SMALL: return -sign1*acos(1.0/sect1)

			sect2 = (B - D) / A
			if sect2 > 0.0:
				sign2 = 1.0
			else:
				sign2 = -1.0
			try: tant2 = sqrt((sect2-1.0)*(sect2+1.0))
			except:
				sect2 = sign2
				tant2 = 0.0

			dx3 = abs(self.c2*sect2 + self.c3*tant2 - x)
			dy3 = abs(self.c5*sect2 + self.c6*tant2 - y)
			if dx3<SMALL and dy3<SMALL: return sign2*acos(1.0/sect2)

			dx4 = abs(self.c2*sect2 - self.c3*tant2 - x)
			dy4 = abs(self.c5*sect2 - self.c6*tant2 - y)
			if dx4<SMALL and dy4<SMALL: return -sign2*acos(1.0/sect2)

			# find the smallest distance
			a1 = dx1**2 + dy1**2
			a2 = dx2**2 + dy2**2
			a3 = dx3**2 + dy3**2
			a4 = dx4**2 + dy4**2
			# else return the smallest one
			if a1<=a2 and a1<=a3 and a1<=a4: return  sign1*acos(1.0/sect1)
			if a2<=a1 and a2<=a3 and a2<=a4: return -sign1*acos(1.0/sect1)
			if a3<=a1 and a3<=a2 and a3<=a4: return  sign2*acos(1.0/sect2)
			return -sign2*acos(1.0/sect2)

		else:
			#if abs(D)<SMALL:
			sect1 = B / A
			if sect1 > 0.0:
				sign = 1.0
			else:
				sign = -1.0
			try: tant1 = sqrt((sect1-1.0)*(sect1+1.0))
			except:
				sect1 = sign
				tant1 = 0.0

			dx1 = abs(self.c2*sect1 + self.c3*tant1 - x)
			dy1 = abs(self.c5*sect1 + self.c6*tant1 - y)
			if dx1<SMALL and dy1<SMALL:
				return sign*acos(1.0/sect1)

			dx2 = abs(self.c2*sect1 - self.c3*tant1 - x)
			dy2 = abs(self.c5*sect1 - self.c6*tant1 - y)

			if dx1**2 + dy1**2 < dx2**2 + dy2**2:
				return sign*acos(1.0/sect1)
			else:
				return -sign*acos(1.0/sect1)

		return None

	# --------------------------------------------------------------------
	def toStr(self, idx=""):
		s = "x%s(t)="%(idx)
		s += _Sx(self.c1)
		s += _Sx(self.c2,"1.0/cos(t)")
		s += _Sx(self.c3,"tan(t)")
		s += "; y%s(t)="%(idx)
		s += _Sx(self.c4)
		s += _Sx(self.c5,"1.0/cos(t)")
		s += _Sx(self.c6,"tan(t)")
		return s

#===============================================================================
# Convert a conic to a parametric equation of a parabola
#	x = c1 + c2*t + c3*t**2
#	y = c4 + c5*t + c6*t**2
#===============================================================================
class ParabolaParametric(ConicParametric):
	def __init__(self, conic):
		ConicParametric.__init__(self, conic)
		# convert conic to standard one
		theta = conic.theta()
		if abs(theta)<SMALL:
			ct = 1.0
			st = 0.0
			a = conic.a
			b = conic.b
			g = conic.g
			f = conic.f
		else:
			ct = cos(theta)
			st = sin(theta)
			stct = st*ct
			a = conic.a*ct**2 + 2.0*conic.h*stct + conic.b*st**2
			b = conic.a*st**2 - 2.0*conic.h*stct + conic.b*ct**2
			g =  conic.g*ct + conic.f*st
			f = -conic.g*st + conic.f*ct

		# Assume that b=ZERO otherwise interchange axes
		if abs(a)<SMALL:
			# with the following translation
			dx = 0.5*(f**2/b - conic.c)/g
			dy = - f/b
			# should reduce to: b*y^2 + 2*g*x = 0
			# or y^2 = 2*sqrt(-D/I^3)*x
			# or y^2 = -2*g/b*x
			k = -0.5*b/g

			self.c2 = -st
			self.c3 =  ct*k

			self.c5 =  ct
			self.c6 =  st*k
		else:
			# with the following translation
			dx = - g/a
			dy = 0.5*(g**2/a - conic.c)/f
			# should reduce to: a*x^2 + 2*f*x = 0
			# or y^2 = 2*sqrt(-D/I^3)*x
			# or x^2 = -2*f/a * y
			k = -0.5*a/f

			self.c2 =  ct
			self.c3 = -st*k

			self.c5 =  st
			self.c6 =  ct*k

		self.c1 =  dx*ct - dy*st
		self.c4 =  dx*st + dy*ct

	# --------------------------------------------------------------------
	# Evaluate a position t
	# --------------------------------------------------------------------
	def __call__(self, t):
		"""Return (x,y) position of t"""
		return (self.c1 + (self.c2 + self.c3*t)*t, \
			self.c4 + (self.c5 + self.c6*t)*t)

	# --------------------------------------------------------------------
	# Return t from point
	# --------------------------------------------------------------------
	def getT(self, x, y):
		"""Return t value from a point (x,y) beloning to the conic"""
		if abs(self.c3)<SMALL:
			return (x - self.c1) / self.c2

		#           A        B       C
		# equation c3*t^2 + c2*t + c1-x = 0
		A = self.c3
		B = self.c2
		C = self.c1 - x
		D = self.c2**2 - 4.0*self.c3*C

		if D > SMALL:
			D = sqrt(D)
			t1 = (-self.c2 + D) / (2.0*self.c3)

			# Evaluate the parametric equation and compare against x,y
			dx1 = abs(self.c1 + self.c2*t1 + self.c3*t1**2 - x)
			dy1 = abs(self.c4 + self.c5*t1 + self.c6*t1**2 - y)
			if dx1<SMALL and dy1<SMALL: return t1

			t2 = (-self.c2 - D) / (2.0*self.c3)
			dx2 = abs(self.c1 + self.c2*t2 + self.c3*t2**2 - x)
			dy2 = abs(self.c4 + self.c5*t2 + self.c6*t2**2 - y)

			# find the smallest distance
			if dx1**2 + dy1**2 < dx2**2 + dy2**2:
				return t1
			else:
				return t2

		else:
			#if abs(D)<SMALL:
			return -self.c2 / (2.0*self.c3)

		return None

	# --------------------------------------------------------------------
	def toStr(self, idx=""):
		s = "x%s(t)="%(idx)
		s += _Sx(self.c1)
		s += _Sx(self.c2,"t")
		s += _Sx(self.c3,"t**2")
		s += "; y%s(t)="%(idx)
		s += _Sx(self.c4)
		s += _Sx(self.c5,"t")
		s += _Sx(self.c6,"t**2")
		return s

#===============================================================================
# Conic - # Gneralized conic surface
#===============================================================================
class Conic:
	# --------------------------------------------------------------------
	# Generalized conic surface
	# S = ax^2 + 2hxy + by^2 + 2gx + 2fy + c = 0
	#        /x\T  / a h g \   /x\
	#    S = |y| * | h b f | * |y| = 0
	#        \1/   \ g f c /   \1/
	#
	#        | a h g |
	#    D = | h b f |
	#        | g f c |
	# --------------------------------------------------------------------
	def __init__(self, a, h, b, g, f, c):
		"""
		Initialize a conic in the form
		S = ax^2 + 2hxy + by^2 + 2gx + 2fy + c = 0

		WARNING: the factor 2 in the coefficients
		"""
		self.a = a
		self.h = h
		self.b = b
		self.g = g
		self.f = f
		self.c = c
		self.invariants()
		self.param = None	# Parametric representation of conic

	# --------------------------------------------------------------------
	def invariants(self):
		"""Calculate invariants and type of conic"""

		# check for zero
		if abs(self.a)<TOOSMALL: self.a = 0.0
		if abs(self.h)<TOOSMALL: self.h = 0.0
		if abs(self.b)<TOOSMALL: self.b = 0.0
		if abs(self.g)<TOOSMALL: self.g = 0.0
		if abs(self.f)<TOOSMALL: self.f = 0.0
		if abs(self.c)<TOOSMALL: self.c = 0.0

		# invariants
		self.D = self.a*(self.b*self.c-self.f**2) \
		       - self.h*(self.h*self.c-self.f*self.g) \
		       + self.g*(self.h*self.f-self.b*self.g)
		if abs(self.D) < TOOSMALL: self.D = 0.0

		self.I = self.a + self.b
		if abs(self.I) < TOOSMALL: self.I = 0.0

		self.J = self.a*self.b - self.h**2
		if abs(self.J) < TOOSMALL: self.J = 0.0

		# type of conic
		self.type = self.conicType()

	# --------------------------------------------------------------------
	def clone(self):
		"""Return a copy of the conic"""
		return Conic(self.a, self.h, self.b, self.g, self.f, self.c)

	# --------------------------------------------------------------------
	def __neg__(self):
		"""Negate conic"""
		return Conic(-self.a, -self.h, -self.b, -self.g, -self.f, -self.c)

	# --------------------------------------------------------------------
	def __add__(self, b):
		"""Add two conics"""
		c = Conic(self.a, self.h, self.b, self.g, self.f, self.c)
		c += b
		return c

	# --------------------------------------------------------------------
	def __iadd__(self, c):
		"""Add inplace conic"""
		self.a += c.a
		self.h += c.h
		self.b += c.b
		self.g += c.g
		self.f += c.f
		self.c += c.c
		self.invariants()
		return self

	# --------------------------------------------------------------------
	def __mul__(self, r):
		"""Multiply conic with r"""
		c = self.clone()
		c *= r
		return c
	__rmul__ = __mul__

	# --------------------------------------------------------------------
	def __imul__(self, r):
		"""Multiply inplace conic with r"""
		self.a *= r
		self.h *= r
		self.b *= r
		self.g *= r
		self.f *= r
		self.c *= r
		self.invariants()
		return self

	# --------------------------------------------------------------------
	def conicType(self):
		"""Determine conic type"""
		if self.D == 0.0:
			# lines?
			if self.J < 0.0:
				# Real intersecting lines
				return CONIC_LINES
			elif self.J == 0.0:
				if abs(self.a) < TOOSMALL and abs(self.b) < TOOSMALL:
					if abs(self.g) < TOOSMALL and abs(self.f) < TOOSMALL:
						# Invalid conic
						return CONIC_NOTHING
					else:
						# Single line a=h=b = 0.0
						return CONIC_LINE
				q = self.c*self.a - self.g**2 + \
				    self.b*self.c - self.f**2
				if abs(q) < TOOSMALL:
					# Coincident lines
					return CONIC_LINES
				elif q < 0.0:
					# Parallel lines
					return CONIC_LINES
				else:
					# Conjugate complex parallel lines
					return CONIC_NOTHING
			else:
				# Conjugate complex intersecting lines
				return CONIC_NOTHING
		elif self.J > 0.0:
			if self.D / self.I < 0.0:
				# Real ellipse
				return CONIC_ELLIPSE
			else:
				# Virtual ellipse
				return CONIC_NOTHING
		elif self.J < 0.0:
			# Hyperbola
			return CONIC_HYPERBOLA
		elif self.J == 0.0:
			# Parabola
			return CONIC_PARABOLA

		# Degenerate conic
		return CONIC_NOTHING

	# --------------------------------------------------------------------
	# rotate axes of conic around an angle theta
	# WARNING to rotate conic use -theta
	# --------------------------------------------------------------------
	def rotate(self, theta):
		"""Rotate conic axes by theta"""
		ct = cos(theta)
		st = sin(theta)
		stct = st*ct

		a = self.a
		b = self.b
		h = self.h
		g = self.g
		f = self.f

		self.a = a*ct**2 + 2.0*h*stct + b*st**2
		self.h = (b-a)*stct + h*(ct**2-st**2)
		self.b = a*st**2 - 2.0*h*stct + b*ct**2
		self.g =  g*ct + f*st
		self.f = -g*st + f*ct

		# check for zero
		if abs(self.a)<TOOSMALL: self.a = 0.0
		if abs(self.h)<TOOSMALL: self.h = 0.0
		if abs(self.b)<TOOSMALL: self.b = 0.0
		if abs(self.g)<TOOSMALL: self.g = 0.0
		if abs(self.f)<TOOSMALL: self.f = 0.0

	# --------------------------------------------------------------------
	# Translate axes of conic by dx, dy
	# WARNING: to translate the conic use -dx, -dy
	# --------------------------------------------------------------------
	def translate(self, dx, dy):
		"""Translate conic axes by (dx,dy)"""
		self.c += self.a*dx**2 + self.b*dy**2 + \
			 2.0*(self.h*dx*dy +  self.g*dx + self.f*dy)
		self.g += self.a*dx + self.h*dy
		self.f += self.h*dx + self.b*dy

		# check for zero
		if abs(self.g)<TOOSMALL: self.g = 0.0
		if abs(self.f)<TOOSMALL: self.f = 0.0
		if abs(self.c)<TOOSMALL: self.c = 0.0

	# --------------------------------------------------------------------
	# Rotation of axes to eliminate the h component
	# --------------------------------------------------------------------
	def theta(self):
		"""Return rotation of axes to eliminate the h component"""
		return 0.5 * atan2(2.0*self.h, self.a - self.b)

	# --------------------------------------------------------------------
	# Return a conic in the standard form
	# FIXME not correct!
	# --------------------------------------------------------------------
	def standard(self):
		"""Return a conic in the standard form"""
		# Solve the equation
		#     l^2 - I*l + J = 0
		D = self.I**2 - 4.0*self.J
		if D >= 0:
			D = sqrt(D)
			l1 = 0.5 * (-self.I - D)
			l2 = 0.5 * (-self.I + D)
			if self.D == 0.0:
				# XXX not correct!
				# D = 0, but if J=0 then D/J = anything!
				return Conic(l1, l2, 0.0, 0.0, 0.0, 0.0)
			else:
				return Conic(l1, l2, 0.0, 0.0, 0.0, self.D/self.J)
		return None

	# --------------------------------------------------------------------
	def __str__(self):
		"""Return a string represetation of the conic"""
		s = "S="
		s += _Sx(self.a, "x^2")
		s += _Sx(2*self.h, "x*y")
		s += _Sx(self.b, "y^2")
		s += _Sx(2*self.g, "x")
		s += _Sx(2*self.f, "y")
		s += _Sx(self.c)
		s += "=0"

		return s + "\t\tD=%g, I=%g, J=%g type=%s" % \
			(self.D, self.I, self.J, CONIC_TYPES[self.type])

	# --------------------------------------------------------------------
	def __call__(self, x, y):
		"""Evaluate conic at position (x,y)"""
		return	self.a*x**2  + self.b*y**2 + self.c + \
			+ 2.0 * (self.h*x*y + self.g*x + self.f*y)

	# --------------------------------------------------------------------
	# Normalized gradient
	# --------------------------------------------------------------------
	def grad(self, x, y):
		dx = self.a*x + self.h*y + self.g
		dy = self.h*x + self.b*y + self.f
		f  = hypot(dx,dy)
		return dx/f,dy/f

	# --------------------------------------------------------------------
	def matrix(self):
		"""Return conic matrix"""
		return Matrix([ [self.a, self.h, self.g],
				[self.h, self.b, self.f],
				[self.g, self.f, self.c] ])

	# --------------------------------------------------------------------
	# Return two single line conics from a degenerate one
	# --------------------------------------------------------------------
	def splitLines(self):
		"""Split a degenerate conic into two lines"""
		if self.type != CONIC_LINES:
			raise Exception("Not a lines conic")

		# convert conic to standard one
		theta = self.theta()

		# do the rotation in place
		if abs(theta) < TOOSMALL:
			ct = 1.0
			st = 0.0
			a = self.a
			b = self.b
			g = self.g
			f = self.f
		else:
			ct = cos(theta)
			st = sin(theta)
			stct = st*ct
			a = self.a*ct**2 + 2.0*self.h*stct + self.b*st**2
			b = self.a*st**2 - 2.0*self.h*stct + self.b*ct**2
			g =  self.g*ct + self.f*st
			f = -self.g*st + self.f*ct

		# translation needed
		if abs(a*b) < TOOSMALL:
			dx = dy = 0.0
			if abs(a)<TOOSMALL:
				dy = - f/b
			else:
				dx = - g/a
		else:
			dx = - g/a
			dy = - f/b
		c  = self.c + g*dx + f*dy

		# Create two new conics
		if abs(c) < TOOSMALL:
			sa = sqrt(abs(a))/2.0
			sb = sqrt(abs(b))/2.0
			conic1 = Conic(0.0,0.0,0.0, sa, sb, 0.0)
			conic2 = Conic(0.0,0.0,0.0, sa,-sb, 0.0)
		elif abs(a) < TOOSMALL:
			# b*y^2 + c = 0
			sb = sqrt(abs(b))/2.0
			sc = sqrt(abs(c))
			conic1 = Conic(0.0,0.0,0.0, 0.0, sb, sc)
			conic2 = Conic(0.0,0.0,0.0, 0.0, sb,-sc)
		else:
			#elif abs(b) < TOOSMALL:
			# a*x^2 + c = 0
			sa = sqrt(abs(a))/2.0
			sc = sqrt(abs(c))
			conic1 = Conic(0.0,0.0,0.0, sa, 0.0, sc)
			conic2 = Conic(0.0,0.0,0.0, sa, 0.0,-sc)

		conic1.translate(-dx,-dy)
		conic1.rotate(-theta)
		conic2.translate(-dx,-dy)
		conic2.rotate(-theta)
		return conic1, conic2

	# --------------------------------------------------------------------
	def parametric(self):
		"""Return a parametric class of the conic"""
		if self.type == CONIC_LINE:
			self.param = LineParametric(self)
		elif self.type == CONIC_ELLIPSE:
			self.param = EllipseParametric(self)
		elif self.type == CONIC_PARABOLA:
			self.param = ParabolaParametric(self)
		elif self.type == CONIC_HYPERBOLA:
			self.param = HyperbolaParametric(self)
		else:
			self.param = None
		return self.param

	# --------------------------------------------------------------------
	# Intersect conic with another one
	# --------------------------------------------------------------------
	def intersect(self, conic):
		"""Intersect self with conic and return intersection points"""
		#print "\nINTERSECT:"
		#print "  a:", self
		#print "  b:", conic

		if self.type == CONIC_LINE:
			if conic.type == CONIC_LINE:
				#print "+++ INTERSECT LINE-LINE"
				d = self.g*conic.f - self.f*conic.g
				if abs(d) < SMALL: return []		# Parallel lines
				x = ( self.f*conic.c - self.c*conic.f ) / d / 2.0
				y = ( self.c*conic.g - self.g*conic.c ) / d / 2.0
				return [(x,y)]
			else:
				return conic._intersectLine(self)

		elif conic.type == CONIC_LINE:
			return self._intersectLine(conic)

		else:
			#print "+++ INTERSECT CONIC-CONIC"
			E1 = self.matrix()
			#print "E1=\n",E1
			E2 = conic.matrix()
			#print "E2=\n",E2
			E2i = -E2
			#print "E2i.det=",E2i.det()
			E2i.inv()
			#print "E2i=\n", E2i
			EE = E1 * E2i
			#print "EE=\n", EE

			a = -EE.trace()
			b = (EE[0][0]*EE[1][1] - EE[1][0]*EE[0][1] +
			     EE[1][1]*EE[2][2] - EE[2][1]*EE[1][2] +
			     EE[0][0]*EE[2][2] - EE[2][0]*EE[0][2])
			c = -EE.det()
			#print a, b, c

			# Find the two conic lines
			m = None
			l = None
			for r in cubic(a,b,c):
				if isinstance(r,float):
					C0 = self + r*conic
					#print "C0: ",C0
					if C0.type == CONIC_LINES:
						m, l = C0.splitLines()
						break;
					elif C0.type == CONIC_LINE:
						m = C0;
						break

			if m is None: return []

			p = self._intersectLine(m)
			if l is not None: p.extend(self._intersectLine(l))
			#print p

			return p

	# --------------------------------------------------------------------
	# Intersect conic with a conic-line
	# --------------------------------------------------------------------
	def _intersectLine(self, line):
		#print "+++ INTERSECT CONIC-LINE"

		# Find a point
		if abs(line.g) < abs(line.f):
			x0 = 0.0
			y0 = -line.c / (2.0*line.f)
		else:
			x0 = -line.c / (2.0*line.g)
			y0 = 0.0

		# end the slope
		vx =  line.f
		vy = -line.g

		#print "  x(t) = %g + %g*t"%(x0,vx)
		#print "  y(t) = %g + %g*t"%(y0,vy)

		# line is in the form
		#    x = x0 + vx*t
		#    y = y0 + vy*t
		# substituting we end up with an equation like
		#    A*t^2 + 2B*t + C = 0
		A = self.a*vx**2 + self.b*vy**2 + 2.0*self.h*vx*vy
		B = self.a*x0*vx + self.h*(x0*vy + y0*vx) + self.b*y0*vy + self.g*vx + self.f*vy
		C = self.a*x0**2  + self.b*y0**2 + self.c + \
			+ 2.0 * (self.h*x0*y0 + self.g*x0 + self.f*y0)
		D = B**2 - A*C

		#print "A=",A
		#print "B=",B
		#print "C=",C
		#print "D=",D
		if abs(D) < SMALL:
			# only at one point
			t = -B / A
			return [(x0+vx*t, y0+vy*t)]
		elif D > 0.0:
			# two points
			D = sqrt(D)
			t1 = (-B + D) / A
			t2 = (-B - D) / A
			return [(x0+vx*t1, y0+vy*t1),
				(x0+vx*t2, y0+vy*t2) ]
		else:
			return []

	# --------------------------------------------------------------------
	# return a list of points beloning to conic
	# --------------------------------------------------------------------
	def draw(self, xlow, xhigh, xstep):
		points = []
		if abs(self.b) < SMALL:
			if     abs(self.h) < SMALL \
			   and abs(self.a) < SMALL \
			   and abs(self.g) < SMALL:
				#vertical line
				pass
			else:
				for x in frange(xlow, xhigh, xstep):
					deno = 2.0*(self.h*x + self.f)
					if abs(deno) > SMALL:
						points.append((x,
							-(self.a*x**2+2.0*self.g*x+self.c)/deno))
		else:
			A = self.b
			for x in frange(xlow, xhigh, xstep):
				B = self.f + self.h * x
				C = self.c + self.a * x**2 + 2.0*self.g*x

				D = B**2 - A*C
				if abs(D) < SMALL:
					points.append((x, -B/A))
				elif D > 0.0:
					D = sqrt(D)
					points.append((x, (-B+D)/A))
					points.append((x, (-B-D)/A))
		return points
