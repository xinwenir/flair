#!/usr/bin/python
# -*- coding: latin1 -*-
#
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	15-Nov-2007

import sys
import math
import bmath
import geometry

from log import say

_SQRT3  = 1.732050807568877293527446341505872366943

#===============================================================================
# Camera Calibration class
#===============================================================================
class Camera:
	"""Camera Calibration class"""
	def __init__(self):
		self.marker     = []		# Marker list
		self.marker_err = []		# ... marker error
		self.Umarker    = self.marker	# transformed markers
		self.err        = 0.0		# error of calibration
		self.initParameters()		# Initialize camera parameters
		self.initMatrices()		# Initialize transformation matrices

	# ----------------------------------------------------------------------
	# Override routine to set camera & lens parameters
	# ----------------------------------------------------------------------
	def initParameters(self):
		"""Set camera parameters"""
		self.calibrated = False		# Camera calibrated
		self.nx = 0			# Number of horizontal pixels
		self.ny = 0			# Number of vertical pixels
		self.dx = 1.0			# Pixel width in user_units
		self.dy = 1.0			# -//- height    -//-
		self.cx = 0.0			# center of screen
		self.cy = 0.0			#  -//-
		self.f  = 1.0			# Camera focal length in user_units
		self.s  = 1.0			# Scaling factor horizontal/vertical
		self.k1 = 0.0			# Radial distortion coefficient in 1/user_units^2

	# ----------------------------------------------------------------------
	# Initialize transformation matrices
	# ----------------------------------------------------------------------
	def initMatrices(self):
		"""Initialize transformation matrices"""
		# Transformation parameters
		#    rc = R * r + t
		self.t = bmath.Vector(3)	# Translation vector
		self.R = bmath.Matrix(type=1)	# Rotation matrix 4x4
		self.T = bmath.Matrix(type=1)	# Translation matrix 4x4
						# Scaling matrix
		self.S = bmath.Matrix.scale(self.s, 1.0, 1.0)
		self.M  = self.T * (self.S * self.R)	# Total transformation
		self.IM = self.M.clone()	# Inverse Transformation
		self.IM.inv()
		self.rot = (0.0, 0.0, 0.0)	# Euler rotation angles Rx,Ry,Rz
		self.loc = bmath.Vector(3)	# location of camera
		self.U   = None			# transformation matrix

	# ----------------------------------------------------------------------
	def setSize(self, nx, ny):
		"""Set image width, height"""
		self.nx = nx
		self.ny = ny
		self.cx = nx / 2.0
		self.cy = ny / 2.0

	# ----------------------------------------------------------------------
	# Clear marker list
	# ----------------------------------------------------------------------
	def delMarkers(self):
		del self.marker[:]
		del self.marker_err[:]
		self.Umarker = self.marker
		self.err = 0.0

	# ----------------------------------------------------------------------
	# Add marker
	# ----------------------------------------------------------------------
	def addMarker(self, x, y, z, xi, yi):
		try:
			x  = float(x)
			y  = float(y)
			z  = float(z)
			xi = float(xi)
			yi = float(yi)
		except:
			return True
		xd, yd = self.image2distorted(xi, yi)
		r2 = xd**2 + yd**2
		self.marker.append((x,y,z,xd,yd,r2))
		return False

	# ----------------------------------------------------------------------
	def transformMarkers(self):
		if self.U is None: return
		self.Umarker = []
		for x,y,z,xd,yd,r2 in self.marker:
			xn,yn,zn = self.U * bmath.Vector(x,y,z)
			self.Umarker.append((xn,yn,zn,xd,yd,r2))

	# ----------------------------------------------------------------------
	# Convert image point to distorted image coordinates
	# ----------------------------------------------------------------------
	def image2distorted(self, xi, yi):
		xd = (xi + 0.5 - self.cx) * self.dx
		yd = (self.cy - yi - 0.5) * self.dy
		return xd, yd

	# ----------------------------------------------------------------------
	# Convert distorted image point to image
	# ----------------------------------------------------------------------
	def distorted2image(self, xd, yd):
		xi = xd / self.dx + self.cx - 0.5
		yi = self.cy - yd / self.dy - 0.5
		return xi, yi

	# ----------------------------------------------------------------------
	# convert from distorted to undistorted sensor plane coordinates
	# ----------------------------------------------------------------------
	def distorted2undistorted(self, xd, yd):
		distortion_factor = 1 + self.k1 * (xd**2 + yd**2)
		return (xd*distortion_factor, yd*distortion_factor)

	# ----------------------------------------------------------------------
	# Convert from undistorted to distorted sensor coordinates.
	# The technique involves algebraically solving the cubic polynomial
	#		ru = rd * (1 + kappa1 * rd**2)
	# using the Cardan method.
	# Note: for kappa1 < 0 the distorted sensor plane extends out to
	#       a maximum barrel distortion radius of
	#                    rd = sqrt (-1/(3 * kappa1)).
	#
	#	To see the logic used in this routine try graphing the
	#       above polynomial for positive and negative kappa1's
	# ----------------------------------------------------------------------
	def undistorted2distorted(self, xu, yu):
		if (xu==0 and yu==0) or self.k1==0:
			return (xu, yu)

		ru = math.hypot(xu, yu)
		c = 1.0 / self.k1
		d = -c * ru

		Q = c / 3.0
		R = -d / 2.0
		D = Q**3 + R**2

		if D >= 0:	# one real root
			D = math.sqrt(D)
			S = bmath.cbrt(R + D)
			T = bmath.cbrt(R - D)
			rd = S + T

			if (rd < 0):
				rd = math.sqrt(abs(-1.0/(3.0*self.k1)))
				say("\nWarning: undistorted image point to distorted image point mapping limited by")
				say("\tmaximum barrel distortion radius of %f"%(rd))
				say("\t(xu=%g, yu=%g) -> (xd=%g, yd=%g)\n"%(xu, yu, xu*rd/ru, yu*rd/ru))
		else:	# three real roots *
			D = math.sqrt(-D)
			S = bmath.cbrt(math.hypot(R, D))
			T = math.atan2(D, R) / 3
			sinT = math.sin(T)
			cosT = math.cos(T)

			# the larger positive root is    2*S*cos(T)
			# the smaller positive root is   -S*cos(T) + SQRT(3)*S*sin(T)
			# the negative root is           -S*cos(T) - SQRT(3)*S*sin(T)

			rd = -S * cosT + _SQRT3 * S * sinT	# use the smaller positive root

		l = rd / ru
		return xu*l, yu*l

	# ----------------------------------------------------------------------
	# Compute the rotation matrix elements and the x,y components from the
	# translation vector by solving an overdetermined linear system.
	# The camera distortion cancels out since the solution is
	# based on the xc/yc = xu/yu = xd*dist/(yd*dist) = xd/yd
	# ----------------------------------------------------------------------
	def __ncc_compute_R_tx_ty_s(self):
		A = []
		B = []
		for x,y,z,xd,yd,r2 in self.Umarker:
			A.append( [yd*x, yd*y, yd*z, -xd*x, -xd*y, -xd*z, yd] )
			B.append( xd )

		# Solve over-determined system
		say(A)
		say(B)
		A = bmath.Matrix(A)
		B = bmath.Matrix(B)
		try:
			RT = bmath.solveOverDetermined(A, B)
		except:
			say("ERROR: Calculating R,tx,ty, Singular matrix")
			return True
		A.writeOctave("a.matrix","A")
		B.writeOctave("b.matrix","B")
		bmath.Matrix(RT).writeOctave("rt.matrix","RT")

		sR1 = bmath.Vector(RT[0], RT[1], RT[2])
		R2  = bmath.Vector(RT[3], RT[4], RT[5])
		lenR2 = R2.length()

		# change equations use x!
		if lenR2 < 0.0001:
			say("Error |ty|~=0")
			pass

		# ty could be also negative!! Check later sign of f
		sty = 1.0 / sR1.length()
		ty  = 1.0 / R2.length()

		self.t = bmath.Vector(RT[6]*sty, ty, 0.0)
		self.s = abs(ty / sty)

		R1n = sty * sR1
		R2n = ty  *  R2
		c   = -0.5 * R1n.dot(R2n)
		#say("Ortho-normal Correction c=",c)
		R1nn = R1n + c * R2n
		R2nn = R2n + c * R1n
		R3nn = R1nn.cross(R2nn)

		self.R.make(R1nn, R2nn, R3nn)
		return False

	# ----------------------------------------------------------------------
	# The minimization of (f, tz, k1) is a non-linear problem, that could
	# be solved iteratively by multiple solutions of linear problems (f,tz)
	# and (f,k1)
	# Solve the overdetermined linear system for the solution of the lens f
	# and the Tz translation. At this stage the camera distortion is ignored
	# ----------------------------------------------------------------------
	def __minimize_f_tz(self, k1):
		R = self.R
		t = self.t

		A = []
		B = []
		for x,y,z,xd,yd,r2 in self.Umarker:
			distortion_factor = 1 + k1 * r2
			zz = R[2][0]*x + R[2][1]*y + R[2][2]*z
			A.append([
				(R[0][0]*x + R[0][1]*y + R[0][2]*z + t[0])*self.s,
				-xd*distortion_factor ])
			B.append([ zz*xd*distortion_factor ])

			A.append([
				(R[1][0]*x + R[1][1]*y + R[1][2]*z + t[1]),
				-yd*distortion_factor ])
			B.append([ zz*yd*distortion_factor ])

		# Solve over-determined system
		try:
			PT = bmath.solveOverDetermined(
					bmath.Matrix(A),
					bmath.Matrix(B))
		except:
			say("ERROR: Calculating f,tz, singular matrix!")
			return None

		f  = PT[0]
		tz = PT[1]
		return (f, tz, k1)

	# ----------------------------------------------------------------------
	def __minimize_f_k1(self, tz):
		R = self.R
		t = self.t

		A = []
		B = []
		for x,y,z,xd,yd,r2 in self.Umarker:
			zc = R[2][0]*x + R[2][1]*y + R[2][2]*z + tz
			A.append([
				(R[0][0]*x + R[0][1]*y + R[0][2]*z + t[0])*self.s,
				-xd*zc*r2 ])
			B.append([ zc*xd ])

			A.append([
				(R[1][0]*x + R[1][1]*y + R[1][2]*z + t[1]),
				-yd*zc*r2 ])
			B.append([ zc*yd ])

		# Solve over-determined system
		try:
			PT = bmath.solveOverDetermined(
					bmath.Matrix(A),
					bmath.Matrix(B))
		except:
			say("ERROR: Calculating f,tz, singular matrix!")
			return None

		f  = PT[0]
		k1 = PT[1]
		return (f, tz, k1)

	# ----------------------------------------------------------------------
	# Function to minimize
	# ----------------------------------------------------------------------
	def __ftzk_error(self, f, tz, k1):
		S = 0.0
		R = self.R
		t = self.t
		for x,y,z,xd,yd,r2 in self.Umarker:
			xc = R[0][0]*x + R[0][1]*y + R[0][2]*z + t[0]
			yc = R[1][0]*x + R[1][1]*y + R[1][2]*z + t[1]
			zc = R[2][0]*x + R[2][1]*y + R[2][2]*z + tz
			xu = f * xc
			yu = f * yc
			if abs(zc) < 1e-10: continue
			xu /= zc
			yu /= zc
			distortion_factor = 1 + k1 * r2
			xu2 = xd * distortion_factor
			yu2 = yd * distortion_factor
			S += (xu-xu2)**2 + (yu-yu2)**2
		return S

	# ----------------------------------------------------------------------
	# Compute the tz, and f, k1, using as initial guess the previous
	# solution
	# ----------------------------------------------------------------------
	def __ncc_compute_f_tz_k1(self):
		f  = 0.0	# Initial values
		tz = 0.0
		k1 = 0.0

		ei = self.__ftzk_error(f, tz, k1)

		for i in range(1000):	# abort after 1000 repetitions
			# Minimize f, tz
			ftzk = self.__minimize_f_tz(k1)
			if ftzk:
				f, tz, k1 = ftzk

			ftzk = self.__minimize_f_k1(tz)
			if ftzk:
				f, tz, k1 = ftzk

			en = self.__ftzk_error(f, tz, k1)
			if abs(en-ei)<1e-10: break
			ei = en

		# correct for negative f
		if f > 0.0:
			self.f    = f
		else:
			self.f    = -f
			self.t[0] = -self.t[0]
			self.t[1] = -self.t[1]
			for i in range(3):
				self.R[0][i] = -self.R[0][i]
				self.R[1][i] = -self.R[1][i]
		self.t[2] = tz
		self.k1   = k1

	# ----------------------------------------------------------------------
	# Compute the rotation matrix elements
	# ----------------------------------------------------------------------
	def __cc_compute_R(self):
		A = []
		B = []
		# Assume that z=0.0 therefore divide the xd/yd to get rid
		# of the distortion corrections
		# xd/yd = (r11*x + r12*y + tx) / (r21*x + r22*y + ty)
		for x,y,z,xd,yd,r2 in self.Umarker:
			# assume ty=1
			#          r11   r12    r21    r22   tx
			A.append( [yd*x, yd*y, -xd*x, -xd*y, yd] )
			B.append( xd )

		try:
			RT = bmath.solveOverDetermined(
					bmath.Matrix(A),
					bmath.Matrix(B))
		except:
			say("ERROR: Calculating R,tx Singular matrix")
			return True

		# recover the full matrix with the following equations
		#	r11^2 + r12^2 + r13^2 = 1
		#	r21^2 + r22^2 + r23^2 = 1
		#	r1.dot(r2) = 0.0
		# assume that we have
		#	r11'^2 + r12'^2 + r13'^2 = k^2
		#	r21'^2 + r22'^2 + r23'^2 = k^2
		#	r1'.dot(r2') = 0.0
		# we get the quadratic eq
		#	k^4 - k^2 (r11'^2+r12'^2+r21'^2+r22'^2) + (r11'*r22'-r12'*r21')^2 = 0
		# therefore
		r11 = RT[0]	# assume r11=r11'
		r12 = RT[1]	# -//-
		r21 = RT[2]
		r22 = RT[3]
		k2 = 0.5 * ( (r11**2 + r12**2 + r21**2 + r22**2) \
			+ math.sqrt( ((r11-r22)**2 + (r12+r21)**2) * \
				((r11+r22)**2 + (r12-r21)**2) ) )
		k = math.sqrt(k2)

		# calculate the r13, r23
		# the remaining rij should be divided by k
		r11 /= k
		r12 /= k
		r21 /= k
		r22 /= k

		# check the signs of far away marker to verify
		# the sign of ty
		mi  = -1
		mr2 = 0.0
		i   = 0
		for x,y,z,xd,yd,r2 in self.Umarker:
			if r2>mr2:
				mi = i
				mr2 = r2
			i += 1
		x, y, z, xd, yd, r2 = self.Umarker[mi]

		tx = RT[4] / k
		ty = 1.0 / k
		XD = r11*x + r12*y + tx
		YD = r21*x + r22*y + ty

		def sign(x): return int(math.copysign(1,x))

		if sign(xd)!=sign(XD) or sign(yd)!=sign(YD):
			k = -k
			# reverse the sign
			r11 = -r11
			r12 = -r12
			r21 = -r21
			r22 = -r22

		tx = RT[4] / k
		ty = 1.0 / k

		r13 = math.sqrt(1.0 - (r11**2+r12**2))
		r23 = math.sqrt(1.0 - (r21**2+r22**2))

		# find the sign ambiguity by checking the dot product
		#	r13*r23 = -(r11*r21 + r12*r22)
		# same or opposite sign?
		same = (-(r11*r21 + r12*r22) > 0.0)

		# try switching the sign of r23
		if not same:
			r23 = - r23

		R1 = bmath.Vector(r11, r12, r13)
		R2 = bmath.Vector(r21, r22, r23)
		R3 = R1.cross(R2)

		# translation vector
		self.t = bmath.Vector(tx, ty, 0.0)

		# rotation matrix
		# check later focal length to define the sign
		self.R.make(R1, R2, R3)

	# ----------------------------------------------------------------------
	# Compute the approximate f and tz by ignoring lens distortion
	# ----------------------------------------------------------------------
	def __cc_compute_approx_f_tz(self):
		R = self.R
		t = self.t

		A = []
		B = []
		for x,y,z,xd,yd,r2 in self.Umarker:
			# z = 0.0
			zz = R[2][0]*x + R[2][1]*y
			A.append([ (R[0][0]*x + R[0][1]*y + t[0])*self.s, -xd ])
			B.append([ zz*xd ])

			A.append([ (R[1][0]*x + R[1][1]*y + t[1]), -yd ])
			B.append([ zz*yd ])

		# Solve over-determined system
		try:
			PT = bmath.solveOverDetermined(
					bmath.Matrix(A),
					bmath.Matrix(B))
		except:
			say("ERROR: Calculating f,tz, singular matrix!")
			return True

		self.f = PT[0]
		t[2]   = PT[1]

		# swap the signs on the rotation matrix
		if self.f < 0.0:
			self.f    = -self.f
			self.t[2] = -self.t[2]
			R[0][2] = -R[0][2]
			R[1][2] = -R[1][2]
			R[2][0] = -R[2][0]
			R[2][1] = -R[2][1]

	# ----------------------------------------------------------------------
	# Create transformation matrices and find the Euler rotation
	# Rx, Ry, Rz angles
	# ----------------------------------------------------------------------
	def _compute_transformation(self):
		R = self.R

		# Create transformation matrices
		self.T = bmath.Matrix.translate(self.t)
		self.S = bmath.Matrix.scale(self.s, 1.0, 1.0)
		self.M  = self.T * (self.S * self.R)
		self.IM = self.M.clone()
		self.IM.inv()
		IR = self.R.clone()
		IR.inv()
		self.loc = IR * self.t

		# ROTX(x) * ROTY(y) * ROTZ(z)
		#  cos(z)*cos(y)
		#			sin(z)*cos(y)
		#						-sin(y)
		# -sin(z)*cos(x)+cos(z)*sin(y)*sin(x)
		#			cos(z)*cos(x)+sin(z)*sin(y)*sin(x)
		#						cos(y)*sin(x)
		#  sin(z)*sin(x)+cos(z)*sin(y)*cos(x)
		#			-cos(z)*sin(x)+sin(z)*sin(y)*cos(x)
		#						cos(y)*cos(x)
		rx = -math.atan2(R[1][2], R[2][2])
		ry =  math.asin(R[0][2])
		rz = -math.atan2(R[0][1], R[0][0])
		self.rot = (rx, ry, rz)

	# ----------------------------------------------------------------------
	# dump variables
	# ----------------------------------------------------------------------
	def dump(self):
		# say(variables)
		say("Camera.size: %d x %d"%(self.nx, self.ny))
		say("R - rotation matrix")
		say(self.R)
		say("T - translation matrix")
		say(self.T)
		say("S - pixel scaling matrix")
		say(self.S)
		say("M - total transformation")
		say(self.M)
		say("inv(M)")
		say(self.IM)
		say("rot Rx=%g, Ry=%g, Rz=%g\n" % tuple(map(math.degrees, self.rot)))
		say("loc = ", self.loc)
		say("t   = ", self.t)
		say("s   = ", self.s)
		say("f   = ", self.f)
		say("k1  = ", self.k1)
		say("err = ", self.err)

	# ----------------------------------------------------------------------
	# Compute the total error of the camera calibration
	# ----------------------------------------------------------------------
	def calibrationError(self):
		total_dist2 = 0
		self.marker_err = []

		for x,y,z,xd,yd,r2 in self.marker:
			Rc = self.M * bmath.Vector(x,y,z)
			xu = Rc[0]/Rc[2]*self.f
			yu = Rc[1]/Rc[2]*self.f
			xdc, ydc = self.undistorted2distorted(xu, yu)
			dist2 = (xdc-xd)**2 + (ydc-yd)**2
			total_dist2 = total_dist2 + dist2
			self.marker_err.append(math.sqrt(dist2))

		self.err = math.sqrt(total_dist2 / len(self.marker))
		return self.err

	# ----------------------------------------------------------------------
	# Coplanar camera calibration
	# therefore we keep the default as it is given by the user
	# WARNING: 1) coplanar calibration requires that all points lie on Z=0
	#             plane
	#          2) we cannot determine the aspect ratio
	# ----------------------------------------------------------------------
	def coplanar_calibration(self):
		self.calibrated = False
		say("Coplanar camera calibration")

		if self.nx*self.ny==0:
			raise Exception("Camera size is not set")

		if len(self.marker)<7:
			raise Exception("Coplanar calibration with less than 5 markers")

		for x,y,z,xd,yd,r2 in self.Umarker:
			if abs(z)>1e-10:
				raise Exception("Coplanar calibration Z<>0.0")

		# 1) Compute the rotation matrix and translation vector
		if self.__cc_compute_R(): return True

		# 2) Compute approx f&tx ignoring lens distortion
		if self.__cc_compute_approx_f_tz(): return True

		# 3) Minimize error on f, tz and k1
		#if self.__ncc_compute_f_tz_k1(): return True

		#self.__cc_compute_exact_f_tz()
		self._compute_transformation()

		# Find calibration error
		self.calibrationError()

		self.calibrated = True

	# ----------------------------------------------------------------------
	# Non-coplanar camera calibration
	# ----------------------------------------------------------------------
	def noncoplanar_calibration(self):
		self.calibrated = False
		say("Non-coplanar camera calibration")

		if self.nx*self.ny==0:
			raise Exception("Camera size is not set")

		if len(self.marker)<7:
			raise Exception("Non-coplanar calibration with less than 7 markers")

		# 1) Compute the r11 to r23 and tx, ty
		if self.__ncc_compute_R_tx_ty_s(): return True

		# 2) Compute the tz, f, k1
		if self.__ncc_compute_f_tz_k1(): return True

		# Find transformation matrix and rotation angles
		self._compute_transformation()

		# Find calibration error
		self.calibrationError()

		self.calibrated = True
		return False

	# ----------------------------------------------------------------------
	# Check markers if they are coplanar or non-coplanar and call the
	# appropriate routine
	# ----------------------------------------------------------------------
	def calibration(self):
		# check markers
		points = []
		for x,y,z,xd,yd,r2 in self.marker:
			points.append((x,y,z))
		plane = geometry.fitplane(points)
		if plane is None:
			return self.noncoplanar_calibration()
		else:
			if abs(plane.a)<1e-10 and \
			   abs(plane.b)<1e-10 and \
			   abs(plane.d)<1e-10:
				return self.coplanar_calibration()

			# find transformation matrix so n=z' and pass from zero
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
			self.U = R*T
			say("R:\n",R)
			say("T:\n",T)
			self.transformMarkers()
			# find transformation matrix to
			# convert to z'=0
			return self.coplanar_calibration()

#-------------------------------------------------------------------------------
if __name__ == "__main__":
	import sys
	import random

	bmath._format = "%10.3g"

	rx = bmath.Matrix.rotX(math.radians(-30.0))
	ry = bmath.Matrix.rotY(math.radians(30.0))
	rz = bmath.Matrix.rotZ(math.radians(30.0))

	R = rx*ry*rz
	t  = bmath.Vector(-10.0,-5.0, 100)
	f  = 1000
	s  = 1.0
	say("Rotation Matrix")
	say(R)

	camera = Camera()
	camera.setSize(3000, 2000)
	random.seed(1.0)

	u = bmath.Vector(1.0, 2.0,  3.0)
	v = bmath.Vector(-10.0, 1.0, 1.0)
	o = bmath.Vector(1.0, 2.5,  3.0)

	u = bmath.Vector(1.0, 1.0, 0.0)
	v = bmath.Vector(0.0, 1.0, 0.0)
	o = bmath.Vector(0.0, 0.0,10.0)

	for i in range(10):
		x = (random.random()-0.5)*2.0*100.0
		y = (random.random()-0.5)*2.0*100.0
		z = (random.random()-0.5)*2.0*100.0
		#z = 100
		#p  = bmath.Vector(x, y, z)
		p = x*u + y*v + o
		say(p)
		pc = R*p + t
		xf = f * pc.x() / pc.z()
		yf = f * pc.y() / pc.z()
		xi,yi = camera.distorted2image(xf, yf)
		camera.addMarker(p.x(),p.y(),p.z(),xi,yi)
	camera.calibration()
