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
# Date:	01-Apr-2013

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import sys
import math
import bmath
from log import say

EPS  = 1E-10
EPS2 = EPS**2

#===============================================================================
# Polygon class
# Points are described as vectors, but they have to lie on the X-Y plane
#===============================================================================
class Polygon:
	def __init__(self,P=None):
		if P is not None:
			self.P = P
		else:
			self.P = []

	# ----------------------------------------------------------------------
	def __str__(self):
		if self.orientation() > 0.0:
			orient = "CCW"
		else:
			orient = "CW"
		return "Polygon%s:%s"%(str(self.P),orient)

	# ----------------------------------------------------------------------
	# Add point to polygon list
	# ----------------------------------------------------------------------
	def add(self, p):
		self.P.append(p)

	# ----------------------------------------------------------------------
	def __getitem__(self, i):
		return self.P[i]

	# ----------------------------------------------------------------------
	def __len__(self):
		return len(self.P)

	# ----------------------------------------------------------------------
	# Return previous point in the loop start..end
	# Used for the iterations
	# ----------------------------------------------------------------------
	def prev(self, idx, start=0, end=-1):
		if end>=0 and idx==start:
			return end
		idx -= 1
		if idx < 0:
			return len(self.P)-1
		return idx

	# ----------------------------------------------------------------------
	# Return next point in the loop start..end
	# Used for the iterations
	# ----------------------------------------------------------------------
	def next(self, idx, start=0, end=-1):
		if end>=0 and idx==end:
			return start
		idx += 1
		if idx == len(self.P):
			return 0
		return idx

	# ----------------------------------------------------------------------
	# @return orientation of polygon CW=-1.0, CCW=1.0 (positive or negative) z
	# ----------------------------------------------------------------------
	def orientation(self):
		# Find total angle
		ang = 0.0
		for i,B in enumerate(self.P):
			A = self.P[i-1]
			if i+1 == len(self.P):
				C = self.P[0]
			else:
				C = self.P[i+1]
			#say("A=",A)
			#say("B=",B)
			#say("C=",C)
			AB = B-A
			if AB.normalize() == 0.0: continue
			BC = C-B
			if BC.normalize() == 0.0: continue
			cross = AB^BC
			dot = AB*BC
			# WARNING Don't use the angle from the asin(cross)
			# since it can fail when ang > 90deg then it will return
			# the ang-90deg
			#phi += asin(cross)
			if   dot < -1.0: dot = -1.0
			elif dot >= 1.0: dot =  1.0
			ang += math.copysign(math.acos(dot), cross[2])

		if ang > 0.0:
			return 1.0
		else:
			return -1.0

	# ----------------------------------------------------------------------
	# @param start,end	starting ending points of sub-polygon
	# @param z		z-coordinate defining orientation/normal of polygon
	#			CW=-1, CCW=1
	# @return true if the sub-polygon[start,end] is convex
	# ----------------------------------------------------------------------
	def isConvex(self, start=0, end=-1, z=None):
		#say("# isConvex",start,end,z)
		if z is None: z = self.orientation()
		j = self.next(start,start,end)
		while j != start:
			jp = self.prev(j,start,end)
			jn = self.next(j,start,end)
			AB = self.P[j] - self.P[jp]
			BC = self.P[jn] - self.P[j]
			Y = AB ^ BC
			#say("#\t",jp,j,jn,Y,z*Y[2])
			if z * Y[2] < 0.0: return False
			j = jn

		# Check last closing corner
		jp = self.prev(j,start,end)
		jn = self.next(j,start,end)
		AB = self.P[j] - self.P[jp]
		BC = self.P[jn] - self.P[j]
		Y = AB ^ BC
		#say("#\t",jp,j,jn,Y,z*Y[2])
		if z * Y[2] < 0.0: return False
		return True

	# ----------------------------------------------------------------------
	# @return true if they intersect, false otherwise
	# ----------------------------------------------------------------------
	def segmentIntersect(self, a,b, c,d):
		# Segment1:   A + t*AB   with 0<=t<=1
		# Segment2:   C + l*CD   with 0<=l<=1
		A  = P[a]
		AB = P[b] - A
		C  = P[c]
		CD = P[d] - C

		#      | ABx     CDx |
		#  D = |             |
		#      | ABy     CDy |
		D = AB[0]*CD[1] - AB[1]*CD[0]
		if abs(D)<1e-10: return False	# Lines parallel

		#       | Cx-Ax     CDx |
		#  Dt = |               |
		#       | Cy-Ay     CDy |
		Dt = (C[0]-A[0])*CD[1] - (C[1]-A[1])*CD[0]
		t = Dt/D
		#say(t)
		if t < 1e-10 or t > 1.0-1e-10: return False	# outside segment

		#       | ABx  Cx-Ax |
		#  Dl = |            |
		#       | ABy  Cy-Ay |
		Dl = AB[0]*(C[1]-A[1]) - AB[1]*(C[0]-A[0])
		l = Dl/D
		#say(l)
		if l <= 1e-10 or l > 1.0-1e-10: return False	# outside segment
		return True

	# ----------------------------------------------------------------------
	# Check if a line segment doesn't intersect any other in the polygon
	# ----------------------------------------------------------------------
	def intersect(self, A, B):
		AB = B - A
		C = self.P[-1]
#		#print
#		#print "intersect A=",A,"B=",B
		for i,D in enumerate(self.P):
#			#print i,"C=",C,"D=",D
			CD = D-C

			if A.eq(C, EPS) or A.eq(D,EPS) or \
			   B.eq(C, EPS) or B.eq(D,EPS):
				C = D
				continue

			# check for intersection
			DD = -AB[0]*CD[1] + AB[1]*CD[0]
#			#print "DD=",DD
			if abs(DD)<EPS2:	# parallel
				C = D
				continue

			# In the first segment
			Dm = -(C[0]-A[0])*CD[1] + (C[1]-A[1])*CD[0]
			m = Dm/DD

			# In the second segment on the polygon
			Dn = AB[0]*(C[1]-A[1]) - AB[1]*(C[0]-A[0])
			n = Dn/DD
#			#print "m=",m,"n=",n
			if m>0.0 and m<1.0 and n>0.0 and n<1.0: return True

			C = D
		return False

	# ----------------------------------------------------------------------
	# @return point index that might contain a sub-polygon with the correct
	# orientation z
	# ----------------------------------------------------------------------
	def _nextConvex(self, start, z):
#		#print "="*50
#		#print "_nextConvex",start,z
#		#print "POLY=",self.P
		i = start
		j = self.next(i)
		k = self.next(j)

		A = self.P[i]
		B = self.P[j]
		C = self.P[k]

		# scan backwards until we find something
		while True:
#			#print
#			#print "i,j,k=",i,j,k
#			#print "A=",A
#			#print "B=",B
#			#print "C=",C
			AB = B - A
			BC = C - B
#			#print "AB=",AB
#			#print "BC=",BC
			Y = AB ^ BC
#			#print "Y=",Y,z*Y[2]
#			#print "intersect",i,k,"=",self.intersect(A,C)

			# Same orientation and not intersecting
			if z*Y[2] > 0.0 and not self.intersect(A,C):
				break

			k = j
			j = i
			i = self.prev(i)
			C = B
			B = A
			A = self.P[i]
			if i==start: return start
#		#print "-------------> found at",i
		return i

#		# continue scanning until it fails (to find the start)
#		D = C	# last point (keep it the same always)
#		while True:
#			AB = B - A
#			BC = C - B
#			Y = AB ^ BC
#
#			#print "intersect",i,k,"=",self.intersect(A,D)
#
#			# Wrong orientation or intersecting
#			if z*Y[2] < 0.0 or self.intersect(A,D):
#				break
#
#			k = j
#			j = i
#			i = self.prev(i)
#			C = B
#			B = A
#			A = self.P[i]
#			if i==start: break
#		#print "=============> stop at i=",i
#		return j

	# ----------------------------------------------------------------------
	# @return a list of sub polygon points that are all convex
	# WARNING destroys original polygon
	# ----------------------------------------------------------------------
	def split2Convex(self):
		polygones = []
		z = self.orientation()
		#import pdb; pdb.set_trace()
		#for end in range(len(self.P)):
		#	say(end,"isConvex",self.isConvex(0,end,z))

		# Find first point breaking convexity
		while len(self.P) > 2:
			#print
			start = self._nextConvex(0,z)
			end = self.next(self.next(start))
			#print "#1. convex=",start,"..",end

			# Scan backwards until we break convexity
			while start != end:
				start = self.prev(start)
				# Check if it is convex and
				# the start,end segment doesn't intersect polygon
				if not self.isConvex(start,end) or \
				   self.intersect(self.P[start],self.P[end]):
					start = self.next(start)
					break

			# Scan forward until we break convexity
			while start != end:
				end = self.next(end)
				# Check if it is convex and
				# the start,end segment doesn't intersect polygon
				if not self.isConvex(start,end) or \
				   self.intersect(self.P[start],self.P[end]):
					end = self.prev(end)
					break
			#print "#2. convex=",start,"..",end

			poly = Polygon()
			polygones.append(poly)
			poly.add(self.P[start])	# WARNING Copy pointers
			i = self.next(start)

			# Find intermediate (idx) points to delete
			idx  = []

			while i!=end:
				poly.add(self.P[i])
				idx.append(i)
				i = self.next(i)
			if end != start: poly.add(self.P[end])
			#print "#", start, idx, end
			for i in reversed(sorted(idx)): del self.P[i]
			#print "poly",poly
		self.mergeConvex(polygones)
		return polygones

	# ----------------------------------------------------------------------
	# Merge a list of convex polygons to a bigger one
	# ----------------------------------------------------------------------
	def mergeConvex(self, polygons):
		for i,A in enumerate(polygons[:-1]):
			j = i+1
			#print
			#print i,"A=",A
			while j<len(polygons):
				B = polygons[j]
				#print j,"B=",B
				edge = A._commonEdge(B)
				if edge is not None:
					#print
					#print "Can merge", edge
					# Test if they can merge

					newpoly = Polygon()

					# A: From [ 0 .. edge[0] )
					#print "add A: 0 ..", edge[0]
					newpoly.P.extend(A[:edge[0]])
					#print "P1=",newpoly

					# B: From [edge[0]=edge[1]) to 0)
					#print "add B:",edge[1],".. last"
					newpoly.P.extend(B[edge[1]:])
					#print "P2=",newpoly

					# B: From 0 .. edge[1].previous
					#print "add B: 0.. ",B.prev(edge[1])
					newpoly.P.extend(B[:B.prev(edge[1])])
					#print "P3=",newpoly

					# A: From (edge[0].next=edge[1]).next to last
					n = A.next(edge[0])
					if n>0:
						#print "add A:",n,".. last"
						newpoly.P.extend(A[n:])
					#print "P4=",newpoly

					if newpoly.isConvex():
						#print "new poly is convex", edge
						# replace A
						A = polygons[i] = newpoly
						#print "A=",A
						del polygons[j]
						j = i+1
						#print "continue"
						continue
				j += 1

	# ----------------------------------------------------------------------
	# Check if self and polygon have a common edge
	# WARNING: Both self,polygon should be convex with the same orientation
	# ----------------------------------------------------------------------
	def _commonEdge(self, polygon):
		# First find if they share a common edge in reverse order
		for i,A in enumerate(self.P):
			B = self.P[self.next(i)]
			# Find if A is part of polygon
			for j,C in enumerate(polygon):
				if A is C:
					D = polygon.P[polygon.prev(j)]
					if B is D:
						return i,j
		return None

# ---------------------------------------------------------------
if __name__ == "__main__":
	import sys
	import string

	P = Polygon()
	f = open(sys.argv[1],"r")
	for line in f:
		if len(line.strip())<2: break
		x,y = list(map(float,line.split()))[:2]
		P.add(bmath.Vector(x,y,0.0))
	f.close()

	say(P)
	z = P.orientation()
	say("orientation z=",z)
	say(P.isConvex())
	convex = P.split2Convex()
	for poly in convex:
		for p in poly:
			say(" ".join(map(str,p)))
		say(" ".join(map(str,poly[0])))
