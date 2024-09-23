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

import bmath

#===============================================================================
# Quad- quadric equation
# Generalized quadratic surface
# S(x,y,z) =   Cxx*x^2 + Cyy*y^2 + Czz*z^2
#            + Cxy*xy  + Cxz*xz  + Cyz*yz
#            + Cx*x    + Cy*y    + Cz*z    + C  = 0
#
# Represented in matrix format
#          0       1      2      3
#     0 / Cxx     Cxy/2  Cxz/2  Cx/2 \
#     1 | Cxy/2   Cyy    Cyz/2  Cy/2 |
#     2 | Cxz/2   Cyz/2  Czz    Cz/2 |
#     3 \ Cx/2    Cy/2   Cz/2   C    /
#
#===============================================================================
class Quad:
	def __init__(self, Cxx=0.0, Cyy=0.0, Czz=0.0,
			Cxy=0.0, Cxz=0.0, Cyz=0.0,
			Cx=0.0, Cy=0.0, Cz=0.0, C=0.0):
		self.set(Cxx, Cyy, Czz, Cxy, Cxz, Cyz, Cx, Cy, Cz, C)

	# --------------------------------------------------------------------
	# Set quad parameters
	# --------------------------------------------------------------------
	def set(self,	Cxx=0.0, Cyy=0.0, Czz=0.0,
			Cxy=0.0, Cxz=0.0, Cyz=0.0,
			Cx=0.0, Cy=0.0, Cz=0.0, C=0.0):
		self.Cxx = Cxx
		self.Cyy = Cyy
		self.Czz = Czz
		self.Cxy = Cxy
		self.Cxz = Cxz
		self.Cyz = Cyz
		self.Cx  = Cx
		self.Cy  = Cy
		self.Cz  = Cz
		self.C   = C
		self._check4plane()

	# --------------------------------------------------------------------
	def _check4plane(self,small=1E-15):
		self.plane =	abs(self.Cxx)<=small and \
				abs(self.Cyy)<=small and \
				abs(self.Czz)<=small and \
				abs(self.Cxy)<=small and \
				abs(self.Cxz)<=small and \
				abs(self.Cyz)<=small

	# --------------------------------------------------------------------
	# Evaluate quad a (x,y,z)
	# --------------------------------------------------------------------
	def __call__(self, x, y=0.0, z=0.0):
		if self.plane:
			return	self.Cx*x + self.Cy*y + self.Cz*z + self.C
		else:
			return	self.Cxx*x*x + self.Cyy*y*y + self.Czz*z*z + \
				self.Cxy*x*y + self.Cxz*x*z + self.Cyz*y*z + \
				self.Cx*x + self.Cy*y + self.Cz*z + self.C

	# --------------------------------------------------------------------
	# @return quad matrix
	# --------------------------------------------------------------------
	def matrix(self):
		m = bmath.Matrix(4,4)
		m[0][0] = self.Cxx
		m[1][1] = self.Cyy
		m[2][2] = self.Czz

		m[0][1] = self.Cxy/2.0; m[1][0] = self.Cxy/2.0
		m[0][2] = self.Cxz/2.0; m[2][0] = self.Cxz/2.0
		m[1][2] = self.Cyz/2.0; m[2][1] = self.Cyz/2.0

		m[0][3] = self.Cx/2.0;  m[3][0] = self.Cx/2.0
		m[1][3] = self.Cy/2.0;  m[3][1] = self.Cy/2.0
		m[2][3] = self.Cz/2.0;  m[3][2] = self.Cz/2.0

		m[3][3] = self.C

	# --------------------------------------------------------------------
	# transform the quadric using matrix m     (m*Q*mT)
	# @param m	matrix to transform quadric
	# --------------------------------------------------------------------
	def transform(self, m):
		if self.plane:
			a = self.Cx
			b = self.Cy
			c = self.Cz

			self.Cx = a*m[0][0] + b*m[1][0] + c*m[2][0]
			self.Cy = a*m[0][1] + b*m[1][1] + c*m[2][1]
			self.Cz = a*m[0][2] + b*m[1][2] + c*m[2][2]
			self.C += a*m[0][3] + b*m[1][3] + c*m[2][3]
		else:
			m00 = self.Cxx*m[0][0] + 0.5*(self.Cxy*m[1][0] + self.Cxz*m[2][0] + self.Cx*m[3][0])
			m01 = self.Cxx*m[0][1] + 0.5*(self.Cxy*m[1][1] + self.Cxz*m[2][1] + self.Cx*m[3][1])
			m02 = self.Cxx*m[0][2] + 0.5*(self.Cxy*m[1][2] + self.Cxz*m[2][2] + self.Cx*m[3][2])
			m03 = self.Cxx*m[0][3] + 0.5*(self.Cxy*m[1][3] + self.Cxz*m[2][3] + self.Cx*m[3][3])

			m10 = self.Cyy*m[1][0] + 0.5*(self.Cxy*m[0][0] + self.Cyz*m[2][0] + self.Cy*m[3][0])
			m11 = self.Cyy*m[1][1] + 0.5*(self.Cxy*m[0][1] + self.Cyz*m[2][1] + self.Cy*m[3][1])
			m12 = self.Cyy*m[1][2] + 0.5*(self.Cxy*m[0][2] + self.Cyz*m[2][2] + self.Cy*m[3][2])
			m13 = self.Cyy*m[1][3] + 0.5*(self.Cxy*m[0][3] + self.Cyz*m[2][3] + self.Cy*m[3][3])

			m20 = self.Czz*m[2][0] + 0.5*(self.Cxz*m[0][0] + self.Cyz*m[1][0] + self.Cz*m[3][0])
			m21 = self.Czz*m[2][1] + 0.5*(self.Cxz*m[0][1] + self.Cyz*m[1][1] + self.Cz*m[3][1])
			m22 = self.Czz*m[2][2] + 0.5*(self.Cxz*m[0][2] + self.Cyz*m[1][2] + self.Cz*m[3][2])
			m23 = self.Czz*m[2][3] + 0.5*(self.Cxz*m[0][3] + self.Cyz*m[1][3] + self.Cz*m[3][3])

			m30 = 0.5*(self.Cx*m[0][0] + self.Cy*m[1][0] + self.Cz*m[2][0]) + self.C*m[3][0]
			m31 = 0.5*(self.Cx*m[0][1] + self.Cy*m[1][1] + self.Cz*m[2][1]) + self.C*m[3][1]
			m32 = 0.5*(self.Cx*m[0][2] + self.Cy*m[1][2] + self.Cz*m[2][2]) + self.C*m[3][2]
			m33 = 0.5*(self.Cx*m[0][3] + self.Cy*m[1][3] + self.Cz*m[2][3]) + self.C*m[3][3]

			self.Cxx =      m[0][0]*m00 + m[1][0]*m10 + m[2][0]*m20 + m[3][0]*m30
			self.Cxy = 2.0*(m[0][0]*m01 + m[1][0]*m11 + m[2][0]*m21 + m[3][0]*m31)
			self.Cxz = 2.0*(m[0][0]*m02 + m[1][0]*m12 + m[2][0]*m22 + m[3][0]*m32)
			self.Cx  = 2.0*(m[0][0]*m03 + m[1][0]*m13 + m[2][0]*m23 + m[3][0]*m33)

			self.Cyy =      m[0][1]*m01 + m[1][1]*m11 + m[2][1]*m21 + m[3][1]*m31
			self.Cyz = 2.0*(m[0][1]*m02 + m[1][1]*m12 + m[2][1]*m22 + m[3][1]*m32)
			self.Cy  = 2.0*(m[0][1]*m03 + m[1][1]*m13 + m[2][1]*m23 + m[3][1]*m33)

			self.Czz =      m[0][2]*m02 + m[1][2]*m12 + m[2][2]*m22 + m[3][2]*m32
			self.Cz  = 2.0*(m[0][2]*m03 + m[1][2]*m13 + m[2][2]*m23 + m[3][2]*m33)

			self.C   =      m[0][3]*m03 + m[1][3]*m13 + m[2][3]*m23 + m[3][3]*m33
		#self.normalize()

	# --------------------------------------------------------------------
	# negate quad coefficients
	# --------------------------------------------------------------------
	def negate(self):
		self.Cxx = -self.Cxx
		self.Cyy = -self.Cyy
		self.Czz = -self.Czz
		self.Cxy = -self.Cxy
		self.Cxz = -self.Cxz
		self.Cyz = -self.Cyz
		self.Cx  = -self.Cx
		self.Cy  = -self.Cy
		self.Cz  = -self.Cz
		self.C   = -self.C
