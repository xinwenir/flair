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
# Date:	07-Dec-2007

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import sys
import csg
import math
import bmath
import string
import random
import Input
from log import say

_DENSITY_LIMIT = 0.1	# Ligther than that is ignored
#_FAR_AWAY      = 1E10

#-------------------------------------------------------------------------------
# Create rotation matrix using the X,Y,Z vectors
#-------------------------------------------------------------------------------
def write_matrix(fout, X, Y, Z, T):
	M = bmath.Matrix()
	M.make(X,Y,Z)
	m = M.transpose()
	for i in range(3): m[i][3] = T[i]
#	m.make(X,Y,Z,T)
	fout.write("\tmultmatrix(m=[\n")
	for r in range(4):
		if r==3: comma=""
		else: comma=","
		fout.write("\t\t\t%s%s\n"%(m[r],comma))
	fout.write("\t\t])\n")

#-------------------------------------------------------------------------------
# Create rotation matrix using the Z vector
#-------------------------------------------------------------------------------
def write_matrixZ(fout, Z, T):
	X = Z.orthogonal()
	Y = Z ^ X
	write_matrix(fout,X,Y,Z,T)

##-------------------------------------------------------------------------------
def writePlane0(fout):
	fout.write("\t\tpolyhedron( [[0,INF,0], [INF,-INF/2,0], [-INF,-INF/2,0], [0,0,-INF]],\n")
	fout.write("\t\t\t[[0,1,2], [0,3,1], [1,3,2], [2,3,0]] );");

##-------------------------------------------------------------------------------
## Return a valid povray name
##-------------------------------------------------------------------------------
#def name(s):
#	# First convert all invalid characters to underscore (_)
#	# Minimize the occurrences of _
#	s = s.encode("ascii","replace")
#	new = s.translate(_TRANS)
#	while True:
#		old = new
#		new = new.replace("__","_")
#		if old == new: break
#	return old

#-------------------------------------------------------------------------------
def header(fout, inp):
	title = inp.cards["TITLE"]
	fout.write("// flair converted to openscad\n")
	fout.write("// Title: %s\n\n"%(title[0].extra()))
	fout.write("INF=1000000;\n")
	fout.write("SPH=32;\n")
	fout.write("CYL=32;\n")
	fout.write("\n")

#-------------------------------------------------------------------------------
def exportBodies(fout, inp):
	fout.write("\n// ---- Bodies ----\n")

	bodies = []

	# Loop all cards
	for card in inp.cardlist:
		if not card.isGeo(): continue
		if len(card.tag) != 3: continue
		if card.ignore(): continue
		if card.tag == "END": continue

		name = card.sdum()
		if name == "": continue
		name = "b_"+name
		bodies.append(name)

		for line in str(card).splitlines():
			fout.write("\n// %s\n"%(line))
		fout.write("module %s()\n{\n"%(name))

		if   card.tag == "RPP":
			fout.write("\ttranslate(%s)\n" % (card.bodyP()))
			fout.write("\t\tcube([%g,%g,%g], center=false);\n"% \
				( card.numWhat(2)-card.numWhat(1),
				  card.numWhat(4)-card.numWhat(3),
				  card.numWhat(6)-card.numWhat(5)))

		elif card.tag == "BOX":
			P = card.bodyP()
			X = card.bodyX()
			Y = card.bodyY()
			Z = card.bodyZ()
			lx = X.norm()
			ly = Y.norm()
			lz = Z.norm()
			write_matrix(fout, X,Y,Z,P)
			fout.write("\t\tcube([%g,%g,%g], center=false);\n"%(lx,ly,lz))

		elif card.tag == "SPH":
			fout.write("\ttranslate(%s)\n" % (card.bodyP()))
			fout.write("\t\tsphere(%g,$fn=SPH);\n"%(card.bodyR()))

		elif card.tag == "RCC":
			P = card.bodyP()
			Z = card.bodyZ()
			h = Z.norm()
			write_matrixZ(fout,Z,P)
			fout.write("\tcylinder(h=%g, r=%g, center=false, $fn=CYL);\n" % (h, card.bodyR()))

		elif card.tag == "TRC":
			P = card.bodyP()
			Z = card.bodyZ()
			h = Z.norm()
			write_matrixZ(fout,Z,P)
			fout.write("\tcylinder(h=%g,r1=%g, r2=%g, center=false, $fn=CYL);\n" % \
				(h, card.numWhat(7), card.numWhat(8)))
		#----------------------------------------------------zxw20240827-----For TET, added by zxw	
		elif   card.tag == "TET":
			V0 = card.bodyP1()
			V1 = card.bodyP2()
			V2 = card.bodyPn(3)
			V3 = card.bodyPn(4)
			P  = V1
			u = V0 - V1
			v = V2 - V1
			w = V3 - V1

			lu = u.norm()
			lv = v.norm()
			lw = w.norm()
			write_matrix(fout, u,v,w,P)
			fout.write("\t\tcube([%g,%g,%g], center=false);\n"%(lu,lv,lw))
		#elif card.tag == "ELL":
		#	WriteEllipsoid

		#elif card.tag == "WED":
		#	WriteAngleWedge

		#elif card.tag == "ARB"
		#	WritePolyhedron

		elif card.tag == "XYP":
			fout.write("\ttranslate(%s)\n" % (card.bodyP()))
			writePlane0(fout);

		elif card.tag == "XZP":
			fout.write("\ttranslate(%s)\n" % (card.bodyP()))
			fout.write("\trotate(90,[1,0,0])\n")
			writePlane0(fout);

		elif card.tag == "YZP":
			fout.write("\ttranslate(%s)\n" % (card.bodyP()))
			fout.write("\trotate(90,[0,1,0])\n")
			writePlane0(fout);

		elif card.tag == "PLA":
			Z = card.bodyN()
			Z.norm()
			write_matrixZ(fout, Z, card.bodyP())
			writePlane0(fout);

		elif card.tag in ("XCC", "YCC", "ZCC"):
			fout.write("\ttranslate(%s)\n" % (card.bodyP()))
			if   card.tag == "XCC":
			    fout.write("\t\trotate(90,[0,1,0])\n")
			elif card.tag == "YCC":
			    fout.write("\t\trotate(90,[1,0,0])\n")

			fout.write("\t\t\tcylinder(h=INF,r=%g,center=true,$fn=CYL);"%(card.bodyR()))

		elif card.tag == "REC":
			P = card.bodyP()
			X = card.bodyX()
			Y = card.bodyY()
			Z = card.bodyZ()
			lx = X.norm()
			ly = Y.norm()
			lz = Z.norm()
			write_matrix(fout,X,Y,Z,P)
			fout.write("\tscale([%g,%g,1.0])\n"%(lx,ly))
			fout.write("\t\tcylinder(h=%g, r=1, center=false, $fn=CYL);\n"%(lz))

		elif card.tag in ("XEC", "YEC", "ZEC"):
			fout.write("\ttranslate(%s)\n" % (card.bodyP()))
			if   card.tag == "XCC":
			    fout.write("\t\trotate(90,[0,1,0])\n")
			elif card.tag == "YCC":
			    fout.write("\t\trotate(90,[1,0,0])\n")

			fout.write("\t\tscale([%g,%g,1.0])\n"%(card.numWhat(3),card.numWhat(4)))
			fout.write("\t\t\tcylinder(h=INF, r=1, center=true, $fn=CYL);\n")

		else:
			raise Exception("Unknown body %s"%(card.tag))

		fout.write("\n}\n")
	return bodies

#-------------------------------------------------------------------------------
def _exportRegion(fout, region, expr):
	# Find material
	try:
		material = region["@material"]
	except:
		sys.stderr.write("\nERROR: region %s is not assigned material\n"%(region.sdum()))
		return True

	try: density = float(material.what(3))
	except: density = 0.0
	if density <= _DENSITY_LIMIT: return False

	# Write card
	fout.write("\n// Region: %s\tMaterial: %s\n\n"%(region.sdum(), material.sdum()))
	for line in str(region).splitlines():
		fout.write("// %s\n"%(line))

	fout.write("module %s()\n{\n"%(region.name()))

	# Check expression
	exp = csg.tokenize(expr)
	try:
		csg.exp2rpn(exp)
		if exp and not csg.check(exp):
			pb.append(-1)
			pb.append("Invalid expression")
		else:
			out = ""
			stack = []
			for token in exp:
				if token in ("-", "+", "|"): #, "@"):
					if   token == "+": op = "intersection() {"
					elif token == "-": op = "difference() {"
					elif token == "|": op = "union() {"
					b = stack.pop()
					a = stack.pop()

					if isinstance(a,list) and a[0]==op:
						if isinstance(b,str):
							a.insert(len(a)-1, "\t"+b)
						else:
							a[len(a)-1:len(a)-1] = ["\t"+x for x in b]
						lst = a

					elif isinstance(b,list) and b[0]==op and token in ("+","|"):
						if isinstance(a,str):
							b.insert(1, "\t"+a)
						else:
							b[1:1] = ["\t"+x for x in a]
						lst = b

					else:
						lst = [op]
						if isinstance(a,str):
							lst.append("\t"+a)
						else:
							for i in a: lst.append("\t"+i)
						if isinstance(b,str):
							lst.append("\t"+b)
						else:
							for i in b: lst.append("\t"+i)
						lst.append("}")

					stack.append(lst)

				else:
					stack.append("b_%s();"%(token))

	except csg.CSGException:
		say("Unbalanced parenthesis")

	lst = stack.pop()
	if isinstance(lst,str):
		fout.write("\t%s\n"%(lst))
	else:
		for i in lst:
			fout.write("\t%s\n"%(i))
	fout.write("}\n")

	return True

#-------------------------------------------------------------------------------
def exportRegions(fout, input):
	fout.write("\n\n// ------ Regions ------\n")
	input.regionProperties()

	region = None
	expr   = ""
	regions = []
	for card in input.cards["REGION"]:
		if card.ignore(): continue
		if card.sdum()=="&":
			# append to previous expression
			expr += card.extra()
			continue

		if region is None:
			# Remember starting region
			region = card
			expr   = region.extra()
			continue

		if _exportRegion(fout, region, expr):
			regions.append(region.name())

		# Remember starting region
		region = card
		expr   = region.extra()

	# Last region
	if region:
		if _exportRegion(fout, region, expr):
			regions.append(region.name())

	return regions

#-------------------------------------------------------------------------------
def export(input, filename):
	fout = open(filename, "w")
	input.preprocess()
	header(fout, input)
	bodies = exportBodies(fout, input)

	for body in bodies:
		fout.write("//%s();\n"%(body))

	regions = exportRegions(fout, input)

	fout.write("\n")
	for region in regions:
		fout.write("%s();\n"%(region))
	fout.close()
