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

import re
import sys
import csg
import math
import bmath
import string
import Input

from log import say
from Input import Card
from functools import reduce

#-------------------------------------------------------------------------------
# Convert FLUKA natural material to isotopic composition
#-------------------------------------------------------------------------------
_fluka2mcnp = {
	 1:  1001,
	 2:  2004,
	 3:((3007,0.925), (3006,0.075)),
	 4:  4009,
	 5:((5011,0.801), (5010,0.199)),
	 6:((6012,0.989), (6013,0.011)),
	 7:((7014,0.996), (7015,0.004)),
	 8:  8016,
	 9:  9019,
	11: 11023,
	13: 13027,
	15: 15031,
	25: 25055,
	27: 27059,
	33: 33075,
	35:((35079,0.507), (35081,0.493)),
	36:((36084,0.570), (36078,0.003), (36080,0.023),
	    (36082,0.116), (36083,0.115), (36086,0.173)),
	38: 38090,
	41: 41093,
	43: 43099,
	53: 53127,
	55: 55133,
	56: 56138,
	73: 73181,
	74: 74186,
	75:((75187,0.626), (75185,0.374)),
	79: 79197,
	82:((82208,0.524), (82207,0.221), (82206,0.241), (82204, 0.014)),
	83: 83209,
	90: 90232,
	93: 93237,
}

PAT_CELL     = re.compile(r"(.*?)#(\d+)(.*)", re.DOTALL | re.IGNORECASE)
PAT_RMCELL   = re.compile(r"#\d+", re.DOTALL | re.IGNORECASE)
PAT_IMP      = re.compile(r"IMP_(.)=(.*)", re.IGNORECASE)
PAT_TRANS    = re.compile(r"^(\*?)TR(\d+)(\*?) (.*)", re.DOTALL | re.IGNORECASE)
PAT_MATERIAL = re.compile(r"^(M\d+) (.*)", re.DOTALL | re.IGNORECASE)

PAT_nR       = re.compile(r"^(.*?\s)(\S+)\s*\b(\d+)R\b(.*)", re.DOTALL | re.IGNORECASE)
PAT_nI       = re.compile(r"^(.*?\s)(\S+)\s*\b(\d+)I\s*(\S*)\s*(.*)", re.DOTALL | re.IGNORECASE)
PAT_xM       = re.compile(r"^(.*?\s)(\S+)\s*\b(\d+)M\b(.*)", re.DOTALL | re.IGNORECASE)
PAT_J        = re.compile(r"^(.*?)\bJ\b(.*)", re.DOTALL | re.IGNORECASE)
PAT_nJ       = re.compile(r"^(.*?)\b(\d+)J\b(.*)", re.DOTALL | re.IGNORECASE)

EPS   = 1.0e-15

CELL_LIKE    = 5000

#-------------------------------------------------------------------------------
# Convert a vector to string
#-------------------------------------------------------------------------------
def vec(v): return "%.16g %.16g %.16g"%(v[0], v[1], v[2])

#===============================================================================
# Mcnp Importer/Exporter class
#===============================================================================
class Mcnp:
	#-----------------------------------------------------------------------
	# Mcnp constructor
	#-----------------------------------------------------------------------
	def __init__(self):
		self.title      = "-Untitled-"
		self.input      = None
		self.macroBodies= True
		self.surfaces   = []	# List of Mcnp surfaces
		self.bodies     = {}	# Dictionary of FLUKA bodies
		self.transform  = []
		self.cells      = []
		self.materials  = []
		self.importance = []

	#-----------------------------------------------------------------------
	def setMacroBodies(self, flag=True):
		self.macroBodies = flag

	#-----------------------------------------------------------------------
	# Import FLUKA input
	#-----------------------------------------------------------------------
	def fromFluka(self, input):
		# Find title card
		self.input = input
		try:
			card = input.cards["TITLE"][0]
			self.title = card.extra()
		except:
			pass

		self.input.preprocess()
		self._importSurfaces()
		self._importCells()
		self._importMaterials()
		#self._importImportance()

	#-----------------------------------------------------------------------
	# Read and convert MCNP file to FLUKA
	#-----------------------------------------------------------------------
	def toFluka(self, filename, input):
		self.input = input
		self.read(filename)

	#-----------------------------------------------------------------------
	# Import Surfaces
	#-----------------------------------------------------------------------
	def _importSurfaces(self):
		# Loop all cards
		for card in self.input.cardlist:
			if not card.isGeo(): continue
			if len(card.tag) != 3: continue
			if card.ignore(): continue

			# Add a new body [##, in, out, card]
			nsurf = len(self.surfaces)
			body = [0, "", "", card]
			self.bodies[card.sdum()] = body

			tag = card.tag

			if   tag == "RPP":
				self._importRPP(card, body)
			elif tag == "BOX":
				self._importBOX(card, body)
			elif tag == "SPH":
				self._importSPH(card, body)
			elif tag == "RCC":
				self._importRCC(card, body)
			elif tag == "REC":
				self._importREC(card, body)
			elif tag == "TRC":
				self._importTRC(card, body)
			elif tag == "ELL":
				self._importELL(card, body)
			elif tag == "ARB":
				self._importARB(card, body)
			elif tag == "WED":
				self._importWED(card, body)
			elif tag == "XYP":
				self._importPlane(card, body, 'Z')
			elif tag == "XZP":
				self._importPlane(card, body, 'Y')
			elif tag == "YZP":
				self._importPlane(card, body, 'X')
			elif tag == "PLA":
				self._importPlane(card, body, 'A')
			elif tag in ("XCC", "YCC", "ZCC"):
				self._importInfCylinder(card, body)
			elif tag in ("XEC", "YEC", "ZEC"):
				self._importInfEllipticalCylinder(card, body)
			elif tag == "QUA":
				self._importQUA(card, body)
			elif tag == "END":
				pass	# ignore END card
			else:
				raise Exception("MCNP: unknown FLUKA body type %s"%(tag))

			body[0] = len(self.surfaces) - nsurf

	#-----------------------------------------------------------------------
	def _importRPP(self, card, body):
		if self.macroBodies:
			self._addSurf(body,'-',"RPP %.16g %.16g %.16g %.16g %.16g %.16g" \
				%(card.numWhat(1), card.numWhat(2), card.numWhat(3),
				  card.numWhat(4), card.numWhat(5), card.numWhat(6)))
		else:
			self._addSurf(body,'+',"PX %.16g"%(card.numWhat(1)), str(card))
			self._addSurf(body,'-',"PX %.16g"%(card.numWhat(2)))
			self._addSurf(body,'+',"PY %.16g"%(card.numWhat(3)))
			self._addSurf(body,'-',"PY %.16g"%(card.numWhat(4)))
			self._addSurf(body,'+',"PZ %.16g"%(card.numWhat(5)))
			self._addSurf(body,'-',"PZ %.16g"%(card.numWhat(6)))

	#-----------------------------------------------------------------------
	def _importBOX(self, card, body):
		P = card.bodyP()
		X = card.bodyX()
		Y = card.bodyY()
		Z = card.bodyZ()

		if self.macroBodies:
			self._addSurf(body,'-',"BOX %s %s %s %s" \
				% (vec(P), vec(X), vec(Y), vec(Z)))
		else:
			# Correct the order
			if (X.cross(Y)).dot(Z)<0:
				X, Y = Y, X		# inverse X with Y

			# outer corner
			PO = P + X + Y + Z
			X.norm()
			Y.norm()
			Z.norm()

			# first 3 planes will be
			self._addSurf(body,'+',"P %s %.16g"%(vec(X), P.dot(X)), str(card))
			self._addSurf(body,'+',"P %s %.16g"%(vec(Y), P.dot(Y)))
			self._addSurf(body,'+',"P %s %.16g"%(vec(Z), P.dot(Z)))

			# last 3 planes will be
			self._addSurf(body,'-',"P %s %.16g"%(vec(X), PO.dot(X)))
			self._addSurf(body,'-',"P %s %.16g"%(vec(Y), PO.dot(Y)))
			self._addSurf(body,'-',"P %s %.16g"%(vec(Z), PO.dot(Z)))

	#-----------------------------------------------------------------------
	def _importSPH(self, card, body):
		x = card.numWhat(1)
		y = card.numWhat(2)
		z = card.numWhat(3)
		r = card.numWhat(4)
		s = str(card)
		if x==0.0 and y==0.0 and z==0.0:
			self._addSurf(body,'-',"SO %.16g"%(r),s)
		elif y==0.0 and z==0.0:
			self._addSurf(body,'-',"SX %.16g %.16g"%(x,r),s)
		elif x==0.0 and z==0.0:
			self._addSurf(body,'-',"SY %.16g %.16g"%(y,r),s)
		elif x==0.0 and y==0.0:
			self._addSurf(body,'-',"SZ %.16g %.16g"%(z,r),s)
		else:
			self._addSurf(body,'-',"S %.16g %.16g %.16g %.16g"%(x,y,z,r),s)

	#-----------------------------------------------------------------------
	def _importRCC(self, card, body):
		P = card.bodyP()
		Z = card.bodyZ()
		r = card.bodyR()
		if self.macroBodies:
			self._addSurf(body,'-',"RCC %s %s %.16g"%(vec(P), vec(Z), r))
		else:
			Ptop = P + Z
			Z.norm()

			self._addSurf(body,'-',"CZ %.16g"%(r),str(card),"%s 6j %s"%(vec(P),vec(Z)))
			self._addSurf(body,'+',"P %s %.16g"%(vec(Z), Z.dot(P)))
			self._addSurf(body,'-',"P %s %.16g"%(vec(Z), Z.dot(Ptop)))

	#-----------------------------------------------------------------------
	def _importREC(self, card, body):
		P = card.bodyP()
		X = card.bodyX()
		Y = card.bodyY()
		Z = card.bodyZ()

		if self.macroBodies:
			self._addSurf(body,'-',"REC %s %s %s %s" \
				% (vec(P), vec(X), vec(Y), vec(Z)))
		else:
			Ptop = P + Z

			radiusX = X.length()
			radiusY = Y.length()

			# correct system
			if (X.cross(Y)).dot(Z) < 0.0:
				X = -X	# from left-hand to right-hand system

			# normalize the coordinates
			Z = X.cross(Y)
			Y = Z.cross(X)

			X.norm()
			Y.norm()
			Z.norm()

			a = 1.0/(radiusX**2)
			b = 1.0/(radiusY**2)

			self._addSurf(body,'-',"SQ %.16g %.16g 0 0 0 0 -1 0 0 0"%(a,b),
					str(card),"%s %s %s 3j"%(vec(P),vec(X),vec(Y)))
			self._addSurf(body,'+',"P %s %.16g"%(vec(Z),P.dot(Z)))
			self._addSurf(body,'-',"P %s %.16g"%(vec(Z),Ptop.dot(Z)))

	#-----------------------------------------------------------------------
	def _importTRC(self, card, body):
		P = card.bodyP()
		Z = card.bodyZ()
		baseRadius = card.numWhat(7)
		apexRadius = card.numWhat(8)

		if self.macroBodies:
			self._addSurf(body,'-',"TRC %s %s %.16g %.16g" \
				% (vec(P), vec(Z), baseRadius, apexRadius))
		else:
			Zlen = Z.length()
			Ptop = P + Z
			Z.norm()
			s = (apexRadius-baseRadius)/Zlen

			if abs(s<1e-10):
				# Normal cylinder
				self._addSurf(body,'-',"CZ %.16g"%(baseRadius),
						str(card),"%s 6j %s"%(vec(P),vec(Z)))
				self._addSurf(body,'+',"P %s %.16g"%(vec(Z), Z.dot(P)))
				self._addSurf(body,'-',"P %s %.16g"%(vec(Z), Z.dot(Ptop)))
			else:
				self._addSurf(body,'-',"KZ %.16g %.16g"%((-baseRadius/s),s**2),
						str(card),"%s 6j %s"%(vec(P), vec(Z)))
				self._addSurf(body,'+',"P %s %.16g"%(vec(Z),P.dot(Z)))
				self._addSurf(body,'-',"P %s %.16g"%(vec(Z),Ptop.dot(Z)))

	#-----------------------------------------------------------------------
	def _importELL(self, card, body):
		P1 = card.bodyP1()
		P2 = card.bodyP2()
		C  = 0.5 * (P1 + P2)
		a  = card.numWhat(7) / 2.0
		c  = (P1-C).length()
		b  = math.sqrt(a**2 - c**2)

		# Create axis system
		Z = (P1-P2)
		Z.norm()

		self._addSurf(body,'-',"SQ %.16g %.16g %.16g 0 0 0 -1 0 0 0" % \
			(1.0/b**2, 1.0/b**2, 1.0/a**2),
			str(card),
			"%s 6j %s"%(vec(C),vec(Z)))

	#-----------------------------------------------------------------------
	def _importQUA(self, card, body):
		self._addSurf(body,'-',"GQ %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g"\
			% (card.numWhat(1), card.numWhat(2), card.numWhat(3),
			   card.numWhat(4), card.numWhat(5), card.numWhat(6),
			   card.numWhat(7), card.numWhat(8), card.numWhat(9),
			   card.numWhat(10)))

	#-----------------------------------------------------------------------
	def _importARB(self, card, body):
		if self.macroBodies:
			self._addSurf(body,'-',"ARB %s %s %s %s %s %s %s %s %.16g %.16g %.16g %.16g %.16g %.16g" \
				% (card.bodyPn(1), card.bodyPn(2), card.bodyPn(3),
				   card.bodyPn(4), card.bodyPn(5), card.bodyPn(6),
				   card.bodyPn(7), card.bodyPn(8),
				   card.intWhat(25), card.intWhat(26),
				   card.intWhat(27), card.intWhat(28),
				   card.intWhat(29), card.intWhat(30)))
		else:
			V = [None]*8
			for i in range(8):
				V[i] = card.bodyPn(i+1)
			F = [None]*6
			s = str(card)
			for i in range(6):
				face = card.intWhat(i+25)
				if face == 0: continue
				v1,face = divmod(face, 1000)
				v2,face = divmod(face,  100)
				v3,v4   = divmod(face,   10)

				# use the first three to define the plane
				E1 = V[v2-1] - V[v1-1]
				E2 = V[v3-1] - V[v1-1]
				N = E1.cross(E2)
				N.norm()
				self._addSurf(body,'-',"P %s %.16g"%(vec(N),N.dot(V[v1-1])), s)
				s = None

	#-----------------------------------------------------------------------
	def _importWED(self, card, body):
		P = card.bodyP()
		X = card.bodyX()
		Y = card.bodyY()
		Z = card.bodyZ()

		if self.macroBodies:
			self._addSurf(body,'-',"WED %s %s %s %s" \
				% (vec(P), vec(X), vec(Y), vec(Z)))
		else:
			# Correct the order
			if (X.cross(Y)).dot(Z)<0:
				X, Y = Y, X		# inverse X with Y

			XC = P + X
			YC = P + Y
			ZC = P + Z

			# normal of the plane
			N  = (YC-XC).cross(ZC-XC)
			N.norm()

			# first 3 planes will be
			self._addSurf(body,'+',"P %s %.16g"%(vec(X), P.dot(X)), str(card))
			self._addSurf(body,'+',"P %s %.16g"%(vec(Y), P.dot(Y)))
			self._addSurf(body,'+',"P %s %.16g"%(vec(Z), P.dot(Z)))

			# cutting plane
			self._addSurf(body,'-',"P %s %.16g"%(vec(N), XC.dot(N)))

	#-----------------------------------------------------------------------
	def _importPlane(self, card, body, axis):
		if axis in ('X', 'Y', 'Z'):
			self._addSurf(body,'-',"P%s %.16g"%(axis,card.numWhat(1)),
					str(card))

		else:
			N = card.bodyN()
			P = card.bodyP()
			N.norm()
			self._addSurf(body,'-',"P %s %.16g"%(vec(N),N.dot(P)),
					str(card))

	#-----------------------------------------------------------------------
	def _importInfCylinder(self, card, body):
		u = card.numWhat(1)
		v = card.numWhat(2)
		r = card.numWhat(3)
		axis = card.tag[0]

		s = str(card)
		if u==0.0 and v==0.0:
			self._addSurf(body,'-',"C%s %.16g"%(axis,r),s)
		else:
			if axis=="Y": u,v = v,u		# Reverse for Y
			self._addSurf(body,'-',"C/%s %.16g %.16g %.16g" % \
				(axis,u,v,r),s)

	#-----------------------------------------------------------------------
	def _importInfEllipticalCylinder(self, card, body):
		P     = card.bodyP()
		a,b,c = card.bodyR()
		axis  = card.tag[0]
		if a!=0.0: a = 1.0/a**2
		if b!=0.0: b = 1.0/b**2
		if c!=0.0: c = 1.0/c**2

		# FIXME it should not be correct!!!
		self._addSurf(body,'-',"SQ %.16g %.16g %.16g 0 0 0 -1 %s" % \
				(a,b,c,vec(P)), str(card))

	#-----------------------------------------------------------------------
	# Add a surface to MCNP
	#-----------------------------------------------------------------------
	def _addSurf(self, body, op, surface, comment="", tr=None):
		if op == "+":
			nop = "-"
		else:
			nop = "+"

		n = len(self.surfaces)+1

		inside  = body[1]
		outside = body[2]
		if inside == "":
			inside  = "%s%d" % ( op, n)
			outside = "%s%d" % (nop, n)
		else:
			inside  += " %s%d" % ( op, n)
			outside += ":%s%d" % (nop, n)

		body[1] = inside
		body[2] = outside

		if tr is not None:
			tn = len(self.transform)+1
			self.transform.append("TR%d %s"%(tn,tr))
			surface = "%5d %d %s"%(n, tn, surface)
		else:
			surface = "%5d %s"%(n, surface)

		self.surfaces.append((surface, comment))

	#-----------------------------------------------------------------------
	# Import Cells
	#-----------------------------------------------------------------------
	def _importCells(self):
		# Loop all region
		self.input.regionProperties()
		self.input.addCardProperty("REGION","BIASING","importance",3,4,5,6)

		region = None
		expr   = ""
		for card in self.input["REGION"]:
			if card.ignore(): continue
			if card.sdum() == "&":
				# append to previous expression
				expr += card.extra()
				continue

			if region is None:
				# Remember starting region
				region = card
				expr   = region.extra()
				continue

			self.__importCell(region, expr)

			# Remember starting region
			region = card
			expr   = region.extra()

		# Last one
		if region: self.__importCell(region, expr)

		# Create universe
		# ... XXX assume is composed by only one body
		#self.cells.append(("%d 0 %s"%(len(self.cells)+1,body[2])
		#			,"Universe"))
	#-----------------------------------------------------------------------
	def __importCell(self, region, expr):
		# Find material
		try:
			material = region["@material"]
		except:
			say("\nError: Region %s is not assigned material"%(region.sdum()))
			return True

		density = material.numWhat(3)
		if density <= 0.0:
			cell_expr = "%d 0"%(region["@n"])	# create a void volume
		else:
			n = material["@n"]
			if n<0:
				say("Warning: Predefined material %s is substituted by blkhole"%(material.sdum()))
				n=1
			cell_expr = "%d %d %g"%(region["@n"], n, -density)

			# Add material to the list
			n -= 1
			if n >= len(self.materials):
				self.materials.extend([None] * \
					(n-len(self.materials)+1))
			self.materials[n] = material

		region_expr = expr.split()	# remove spaces
		if len(region_expr)==0: return True
		region_expr = reduce(lambda x,y:x+y, region_expr)

		if region_expr.find("#")>=0 or region_expr.find("(")>=0:
			# Expand parenthesis
			say("Warning: Region \"%s\" contains parentheses. Expanded!"%(region.sdum()))
			exp = csg.tokenize(region_expr)
			csg.exp2rpn(exp)
			csg.rpnorm(exp)
			region_expr = csg.rpn2exp(exp)
		else:
			region_expr = region_expr.replace("+", " + ")
			region_expr = region_expr.replace("-", " - ")
			region_expr = region_expr.replace("|", " | ")
			region_expr = region_expr.split()

		op = "+"
		first = True
		# construct the MCNP region line
		for item in region_expr:
			if item in ("+", "-"):
				op = item
			elif item == "|":
				op = "+"
				if first:
					first = False
				else:
					cell_expr += ":"
			else:
				first = False
				body = self.bodies[item]
				n, inside, outside, bcard = body
				if op=="+":
					cell_expr += " %s"%(inside)
				else:
					if n>1:
						cell_expr += " (%s)"%(outside)
					else:
						cell_expr += " %s"%(outside)
		self.cells.append((cell_expr, str(region)))

		if region["@material"].intWhat(4)==1:
			self.importance.append(0.0)
		else:
			biasing = region["@biasingN"]
			if biasing is None:
				biasing = region["@biasingA"]
			if biasing is not None:
				imp = biasing.numWhat(3)
			else:
				imp = 1.0
			self.importance.append(imp)
		return False

	#-----------------------------------------------------------------------
	# Import Materials that are in use
	#-----------------------------------------------------------------------
	def _importMaterials(self):
		for mat in self.materials:
			if mat is None: continue
			mat.zaid = []
			for za,f in self.input.materialZAID(mat):
				z = za // 1000
				zaid = _fluka2mcnp.get(z)
				if zaid is not None:
					if isinstance(zaid, int):
						mat.zaid.append((zaid,f))
					else:
						for zz,ff in zaid:
							mat.zaid.append((zz,ff*f))
				else:
					mat.zaid.append((za,f))

	#-----------------------------------------------------------------------
	# Read an mcnp file and convert it to FLUKA
	#-----------------------------------------------------------------------
	def read(self, filename):
		fin = open(filename,"r")

		if self.input is None:
			self.input = Input.Input()

		comment = ""
		cmdline = ""
		self.section = 0

		def appendComment(comment, line):
			if comment:
				return comment + "\n" + line
			else:
				return line

		# Read title as first line
		self.__parseCommand("", fin.readline().rstrip())

		# Read the remaining lines
		cont = False
		while True:
			if cont:
				line += fin.readline()
			else:
				line = fin.readline()
			if len(line)==0: break
			line = line.rstrip()
			dollar = line.find("$")
			if dollar>=0:
				comment = appendComment(comment, line[dollar+1:])
				line = line[:dollar].rstrip()

			if line and line[-1] == "&":
				line = line[:-1].rstrip()
				cont = True
				continue
			else:
				cont = False

			if line=="":
				if dollar<0:
					self.__parseCommand(comment, cmdline)
					self.__parseCommand("", None)	# Advance section
					comment = ""
					cmdline = ""

			elif line[0] in ("c", "C"):
				comment = appendComment(comment, line[1:])

			elif line[:5] == "     " or line[0]=="\t":
				# Continuation
				cmdline += line

			else:
				# Execute command
				if cmdline != "":
					self.__parseCommand(comment, cmdline)
					comment = ""
				cmdline = line

		# execute last command
		self.__parseCommand(comment, cmdline)

	#-----------------------------------------------------------------------
	# parse Mcnp command line
	#-----------------------------------------------------------------------
	def __parseCommand(self, comment, command):
		if command:
			command = self.__parseFeatures(command)

		if command is None:
			self.section += 1
			if self.section == 3:
				self.input.addCard(Card("END"))
				regions = self.__parseCells()
				self.input.addCard(Card("END"))
				self.input.addCard(Card("GEOEND"))

				# Create assignmat and optionally biasing
				for region in regions:
					name = region.sdum()
					# FIXME unique material if multiple densities are supplied!
					mat = region["@mat"]
					# biasing
					impn = region.get("@IMP:N","1")
					impp = region.get("@IMP:P","1")
					impe = region.get("@IMP:E","1")

					if mat == "0":
						if impn=="0" and impp=="0" and impe=="0":
							mat = "BLCKHOLE"
						else:
							mat = "VACUUM"
					else:
						mat = "M"+mat

					self.input.addCard(Card("ASSIGNMA", ["", mat, name]))

					if impn!="1" and impn!="0":
						self.input.addCard(Card("BIASING", ["", "3.0", "", impn, name]))
					if impe!="1" and impe!="0":
						self.input.addCard(Card("BIASING", ["", "2.0", "", impe, name]))
					if impp!="1" and impp!="0":
						self.input.addCard(Card("BIASING", ["", "2.0", "", impp, name]))
			return

		# Cell parsing
		if self.section == 0:
			self.input.addCard(Card("TITLE", [],comment,command))
			self.section += 1
			self.input.addCard(Card("GEOBEGIN",
				["COMBNAME"],"",command))

		elif self.section == 1:
			# add cells to be parsed latter
			self.cells.append((command.upper(), comment))

		elif self.section == 2:
			# parse surfaces
			self.__parseSurface(command.upper(), comment)

		else:
			# Generic commands
			m = PAT_MATERIAL.match(command.upper())
			if m:
				self.__parseMaterial(comment, m.group(1), m.group(2))
				return
			m = PAT_TRANS.match(command.upper())
			if m:
				deg = m.group(1)=="*" or m.group(3)=="*"
				self.__parseTransform(comment, m.group(2), m.group(4), deg)
				return

		#if comment: say("C",comment)
		#say(self.section, command)

	#-----------------------------------------------------------------------
	# nR: repeat the proceeding entry n times
	# nI: insert n linear interpolate entries
	# xM: multiply proceeding entry by x
	# nJ: replace by the default value (or jump) n entries (J alone is accepted)
	#-----------------------------------------------------------------------
	def __parseFeatures(self, command):
		# Parse nR
		while True:
			m = PAT_nR.match(command)
			if m is None: break
			command = m.group(1) + \
				  " ".join([m.group(2)]*(int(m.group(3))+1)) + \
				  m.group(4)

		# Parse nI
		while True:
			m = PAT_nI.match(command)
			if m is None: break
			try:
				command = m.group(1) + m.group(2)
				a = float(m.group(2))
				n = int(m.group(3))
				b = float(m.group(4))
				d = (b-a)/(n+1)
				x = a
				for i in range(n):
					x += d
					command += " %g"%(x)
				command += " " + m.group(4) + m.group(5)
			except ValueError:
				break

		# Parse xM
		while True:
			m = PAT_xM.match(command)
			if m is None: break
			try:
				command = m.group(1) + m.group(2) + " " + \
					  str(float(m.group(2))*float(m.group(3))) + \
					  m.group(4)
			except ValueError:
				break

		# Parse single J
		while True:
			m = PAT_J.match(command)
			if m is None: break
			command = m.group(1) + "0" + m.group(2)

		# Parse nJ
		while True:
			m = PAT_nJ.match(command)
			if m is None: break
			command = m.group(1) + " ".join(["0"]*int(m.group(2))) + m.group(3)

		return command

	#-----------------------------------------------------------------------
	# Parse Surfaces to Bodies
	#
	# TODO
	#   * If a surface is preceded by an asterisk * then an reflective
	#     surface is created
	#
	#   * Cones K/XYZ the last variable, single +/-1 is not taken into account
	#
	#   * Surface cards of type X, Y, and Z may be used to describe surfaces
	#     by coordinate points rather than by equation coefficients as in the previous
	#     section. The surfaces described by these cards must be symmetric about the
	#     x-, y- or z-axis, respectively, and if the surface consists of more than one
	#     sheet, the specified coordinate points must all be on the same sheet.
	#     Each of the coordinate pairs defines a geometrical point on the surface.
	#
	#     On the Y card, for example, the entries may be
	#
	#          j  Y  y1 r1  y2 r2
	#
	#     where ri = sqrt(xi^2 + zi^2) and yi is the coordinate of point i.
	#
	#   - If one coordinate pair is used, a plane (PX, PY, or PZ) is defined.
	#
	#   - If two coordinate pairs are used, a linear surface
	#       (PX, PY, PZ, CX, CY, CZ, KX, KY or KZ) is defined.
	#
	#   - If three coordinate pairs are used, a quadratic surface
	#	(PX, PY, PZ, SO, SX, SY, SZ, CX, CY, CZ, KX, KY, KZ or SQ) is defined.
	#
	#   - When a cone is specified by two points, a cone of only one sheet is
	#     generated.
	#-----------------------------------------------------------------------
	def __parseSurface(self, command, comment):
		#if Input._developer: say(">S>",command)
		token   = command.split()
		try:
			sid = token.pop(0)
			sid  = int(sid)
		except ValueError:
			say("Invalid surface id %s"%(sid))
			return False
		name    = "S%d"%(sid)
		surface = token.pop(0)
		try:
			tr = int(surface)
			surface = token.pop(0)
		except ValueError:
			tr = 0
		surface = surface
		arg = list(map(float,token))

		if surface=="P":
			# Plane
			tag  = "PLA"
			what = [name, arg[0], arg[1], arg[2],
			        arg[0]*arg[3], arg[1]*arg[3], arg[2]*arg[3]]
		elif surface=="PX":
			# Plane normal to X-axis
			tag  = "YZP"
			what = [name, arg[0]]
		elif surface=="PY":
			# Plane normal to Y-axis
			tag  = "XZP"
			what = [name, arg[0]]
		elif surface=="PZ":
			# Plane normal to Z-axis
			tag  = "XYP"
			what = [name, arg[0]]

#		elif surface=="X":
#			pass	# Not implemented yet see notes

#		elif surface=="Y":
#			pass	# Not implemented yet see notes

#		elif surface=="Z":
#			pass	# Not implemented yet see notes

		elif surface=="SO":
			# Sphere at Origin
			tag  = "SPH"
			what = [name, 0., 0., 0., arg[0]]
		elif surface=="S":
			# Sphere general
			tag  = "SPH"
			what = [name, arg[0], arg[1], arg[2], arg[3]]
		elif surface=="SX":
			# Sphere centered on X-axis
			tag  = "SPH"
			what = [name, arg[0], 0., 0., arg[1]]
		elif surface=="SY":
			# Sphere centered on Y-axis
			tag  = "SPH"
			what = [name, 0., arg[0], 0., arg[1]]
		elif surface=="SZ":
			# Sphere centered on Z-axis
			tag  = "SPH"
			what = [name, 0., 0., arg[0], arg[1]]

		elif surface=="C/X":
			# Cylinder parallel to X-axis
			tag  = "XCC"
			what = [name, arg[0], arg[1], arg[2]]
		elif surface=="C/Y":
			# Cylinder parallel to Y-axis
			tag  = "YCC"
			what = [name, arg[1], arg[0], arg[2]]
		elif surface=="C/Z":
			# Cylinder parallel to Z-axis
			tag  = "ZCC"
			what = [name, arg[0], arg[1], arg[2]]

		elif surface=="CX":
			# Cylinder on X-axis
			tag  = "XCC"
			what = [name, 0., 0., arg[0]]
		elif surface=="CY":
			# Cylinder on Y-axis
			tag  = "YCC"
			what = [name, 0., 0., arg[0]]
		elif surface=="CZ":
			# Cylinder on Z-axis
			tag  = "ZCC"
			what = [name, 0., 0., arg[0]]

		elif surface=="K/X":
			# Cone parallel to X-axis
			# FIXME the last variable, single is not taken into account
			tag = "QUA"
			x  = arg[0]
			y  = arg[1]
			z  = arg[2]
			t2 = arg[3]
			try:
				single = arg[4]
				say("Warning:", command)
				say("\tsingle side cone is not implemented yet!")
			except: single = None
			what = [name,	-t2, 1.0, 1.0,			# Cxx Cyy Czz
					0.0, 0.0, 0.0,			# Cxy Cxz Cyz
					2.0*t2*x, -2.0*y, -2.0*z,	# Cx  Cy  Cz
					y*y + z*z - t2*x*x]

		elif surface=="K/Y":
			# Cone parallel to Y-axis
			# FIXME the last variable, single is not taken into account
			tag = "QUA"
			x  = arg[0]
			y  = arg[1]
			z  = arg[2]
			t2 = arg[3]
			try:
				single = arg[4]
				say("Warning:", command)
				say("\tsingle side cone is not implemented yet!")
			except: single = None
			what = [name,	1.0, -t2, 1.0,			# Cxx Cyy Czz
					0.0, 0.0, 0.0,			# Cxy Cxz Cyz
					-2.0*x, 2.0*t2*y, -2.0*z,	# Cx  Cy  Cz
					x*x + z*z - t2*y*y]

		elif surface=="K/Z":
			# Cone parallel to Z-axis
			# FIXME the last variable, single is not taken into account
			tag = "QUA"
			x  = arg[0]
			y  = arg[1]
			z  = arg[2]
			t2 = arg[3]
			try:
				single = arg[4]
				say("Warning:", command)
				say("\tsingle side cone is not implemented yet!")
			except: single = None
			what = [name,	1.0, 1.0, -t2,			# Cxx Cyy Czz
					0.0, 0.0, 0.0,			# Cxy Cxz Cyz
					-2.0*x, -2.0*y, 2.0*t2*z,	# Cx  Cy  Cz
					x*x + y*y - t2*z*z]

		elif surface=="KX":
			# Cone on X-axis
			tag  = "QUA"
			x  = arg[0]
			t2 = arg[1]
			try:
				single = arg[2]
				say("Warning:", command)
				say("\tsingle side cone is not implemented yet!")
			except: single = None
			what = [name,	-t2, 1.0, 1.0,		# Cxx Cyy Czz
					0.0, 0.0, 0.0,		# Cxy Cxz Cyz
					2.0*t2*x, 0.0, 0.0,	# Cx  Cy  Cz
					-t2*x*x]

		elif surface=="KY":
			# Cone on Y-axis
			tag = "QUA"
			y  = arg[0]
			t2 = arg[1]
			try:
				single = arg[2]
				say("Warning:", command)
				say("\tsingle side cone is not implemented yet!")
			except: single = None
			what = [name,	1.0, -t2, 1.0,		# Cxx Cyy Czz
					0.0, 0.0, 0.0,		# Cxy Cxz Cyz
					0.0, 2.0*t2*y, 0.0,	# Cx  Cy  Cz
					-t2*y*y]

		elif surface=="KZ":
			# Cone on Z-axis
			tag = "QUA"
			tag = "QUA"
			z  = arg[0]
			t2 = arg[1]
			try:
				single = arg[2]
				say("Warning:", command)
				say("\tsingle side cone is not implemented yet!")
			except: single = None
			what = [name,	1.0, 1.0, -t2,		# Cxx Cyy Czz
					0.0, 0.0, 0.0,		# Cxy Cxz Cyz
					0.0, 0.0, 2.0*t2*z,	# Cx  Cy  Cz
					-t2*z*z]

		elif surface=="SQ":
			# Ellipsoid, Hyperboloid, Paraboloid axes parallel to X,Y,Z-axis
			tag  = "QUA"
			a,b,c, d,e,f, g, x,y,z = arg
			what = [name, a, b, c,
			              2.0*(d-a*x), 2.0*(e-b*y), 2.0*(f-c*z),
				      g + a*x**2 + b*y**2 + c*z**2
					- 2.0*(d*x + e*y + f*z)]

		elif surface=="GQ":
			# Cylinder, Cone, Ellipsoid, Hyperboloid, Paraboloid generic
			tag  = "QUA"
			what = [name, arg[0], arg[1], arg[2],
			              arg[3], arg[4], arg[5],
			              arg[6], arg[7], arg[8],
			              arg[9]]

		elif surface in ("BOX", "REC", "WED"):
			tag  = surface
			what = [name, arg[0], arg[1], arg[2],
			              arg[3], arg[4], arg[5],
			              arg[6], arg[7], arg[8],
			              arg[9], arg[10], arg[11]]

		elif surface=="RPP":
			tag  = surface
			what = [name, arg[0], arg[1], arg[2],
			              arg[3], arg[4], arg[5]]

		elif surface=="SPH":
			tag  = surface
			what = [name, arg[0], arg[1], arg[2], arg[3]]

		elif surface=="RCC":
			tag  = surface
			what = [name, arg[0], arg[1], arg[2],
			              arg[3], arg[4], arg[5],
			              arg[6]]

		elif surface=="TRC":
			tag  = surface
			what = [name, arg[0], arg[1], arg[2],
			              arg[3], arg[4], arg[5],
			              arg[6], arg[7]]

		else:
			raise Exception(command)

		if comment:
			comment += "\n%s"%(command)
		else:
			comment = command
		card = Card(tag, what, comment)

		if tr != 0: self.input.addCard(Card("$start_transform",["","rot%d"%(tr)]))
		self.input.addCard(card)
		if tr != 0: self.input.addCard(Card("$end_transform"))

		return True

	#-----------------------------------------------------------------------
	# Parse MCNP cells and convert them to FLUKA regions
	#-----------------------------------------------------------------------
	def __parseCells(self):
		regions = []
		cells   = {}
		trsurf  = {}

		# Parse cells
		likeBase = 0
		for expr,comment in self.cells:
			#say(">C>", expr)
			expr = expr.replace("IMP:","IMP_")
			expr = expr.replace("("," +( ")
			expr = expr.replace("# +("," -( ")
			expr = expr.replace(")"," ) ")
			expr = expr.replace(":"," : ")
			tokens = expr.split()
			#say("tokens=",tokens)
			idx  = tokens.pop(0)
			name = "C"+idx

			material = tokens.pop(0)
			if material == "LIKE":
				ref = tokens.pop(0)
				if tokens.pop(0) != "BUT": continue
				trcl = tokens.pop(0)
				if trcl[:5] != "TRCL=": continue
				rot = "rot"+trcl[5:]
				region   = cells[ref]
				material = region["@mat"]
				density  = 0
				likeBase += CELL_LIKE

				csg = []
				for token in region.extra().split():
					if len(token)>2 and token[1] == "S":
						bodyid = likeBase+int(token[2:])
						trsurf[bodyid] = rot
						surface = "S%d"%(bodyid)
						csg.append(token[0] + surface)
					else:
						csg.append(token)

				card = Card("REGION", [name, "5"], comment, " ".join(csg))
				card["@mat"]     = region["@mat"]
				card["@density"] = region["@density"]
				self.input.addCard(card)
				cells[idx] = card
				continue

			elif material != "0":
				density = tokens.pop(0)
			else:
				density = 0

			# Create card
			card = Card("REGION", [name, "5"], comment)
			card["@mat"]     = material
			card["@density"] = density

			csg = []
			for token in tokens:
				if token==":":
					csg.append("\n|")

				elif token in ("+(", "-(", ")"):
					csg.append(token)

				elif token[0] == "#":
					csg.append(token)

				elif token.startswith("IMP_"):
					m = PAT_IMP.match(token)
					if m: card["@IMP:%s"%(m.group(1))] = m.group(2)

				else:
					try:
						num = int(token)
						if num<0:
							op = "+"
							num = -num
						else:
							op = "-"
						csg.append("%sS%d"%(op,num))
					except ValueError:
						break

			expr = " ".join(csg)
			card.setExtra(expr)
			self.input.addCard(card)
			cells[idx] = card
			regions.append(card)

		# Correct cell references # as -()
		for card in list(cells.values()):
			expr = card.extra()
			if "#" not in card.extra(): continue
			while True:
				m = PAT_CELL.match(expr)
				if m is None: break
				try:
					ref = cells[m.group(2)]
					# Expand only first level by removing
					refexpr = re.sub(PAT_RMCELL, "", ref.extra())
					expr = "%s -(%s) %s"%(m.group(1),refexpr,m.group(3))
				except KeyError:
					expr = "%s ?C%s %s"%(m.group(1),m.group(2),m.group(3))
			card.setExtra(expr)

		# Create clones of surfaces if needed
		for trid, rot in list(trsurf.items()):
			sid = "S%d"%(trid%CELL_LIKE)
			# Find body and clone it
			for card in self.input.cardlist:
				if len(card.tag) != 3 or card.sdum() != sid: continue
				pos = card.pos()
				clone = card.clone()
				clone.setSdum("S%d"%(trid))

				self.input.addCard(Card("$start_transform",["",rot]),pos+1)
				self.input.addCard(clone,pos+2)
				self.input.addCard(Card("$end_transform"),pos+3)
				self.input.renumber(pos)
				break

		return regions

	#-----------------------------------------------------------------------
	# Parse transformations
	#-----------------------------------------------------------------------
	def __parseTransform(self, comment, idx, transform, deg):
		try:
			transform = list(map(float,transform.split()))
		except:
			say("Error: parsing transformation:", idx, transform)
		if deg:
			for i in range(3,min(len(transform),12)):
				x = math.cos(math.radians(transform[i]))
				if abs(x)<EPS: x = 0.0
				transform[i] = x

		rotidx = "rot%d"%(int(idx))

		pos = bmath.Vector(transform[0], transform[ 1], transform[ 2])
		if len(transform) > 5:
			X = bmath.Vector(transform[3], transform[ 4], transform[ 5])
			X.normalize()
			if len(transform) > 8:
				Y = bmath.Vector(transform[6], transform[ 7], transform[ 8])
				Y.normalize()
				if len(transform) > 11:
					Z = bmath.Vector(transform[9], transform[10], transform[11])
					Z.normalize()
					if len(transform)==13:
						M = int(transform[12])
					else:
						M = 1
				else:
					Z = X ^ Y

			else:
				Y = X.orthogonal()
				Z = X ^ Y

			# Get Euler rotation
			# Rx * Ry * Rz
			#       xx                   xy                    xz
			#  cos(z)*cos(y)
			#			sin(z)*cos(y)
			#						-sin(y)
			#
			#       yx                   yy                    yz
			# -sin(z)*cos(x)+cos(z)*sin(y)*sin(x)
			#			cos(z)*cos(x)+sin(z)*sin(y)*sin(x)
			#						cos(y)*sin(x)
			#
			#       zx                   zy                    zz
			#  sin(z)*sin(x)+cos(z)*sin(y)*cos(x)
			#			-cos(z)*sin(x)+sin(z)*sin(y)*cos(x)
			#						cos(y)*cos(x)
			Rx = -math.degrees(math.atan2(Y.z(), Z.z()))
			Ry =  math.degrees(math.asin(X.z()))
			Rz = -math.degrees(math.atan2(X.y(), X.x()))
		else:
			Rx = Ry = Rz = 0.0

		# Start with the translation
		def addTranslation(what, pos):
			if abs(pos.x())>EPS or abs(pos.y())>EPS or abs(pos.z())>EPS:
				what.extend(pos)

		first = True
		if abs(Rz) > EPS:
			card = Card("ROT-DEFI", [rotidx, "", "", bmath.format(Rz,10)], comment)
			first = False
			self.input.addCard(card)

		if abs(Ry) > EPS:
			card = Card("ROT-DEFI", [rotidx, 200, "", bmath.format(Ry,10)])
			if first:
				card.setComment(comment)
				first = False
			self.input.addCard(card)

		if abs(Rx) > EPS:
			card = Card("ROT-DEFI", [rotidx, 100, "", bmath.format(Rx,10)])
			if first:
				card.setComment(comment)
				first = False
			self.input.addCard(card)

		if first or abs(pos.x())>EPS or abs(pos.y())>EPS or abs(pos.z())>EPS:
			card = Card("ROT-DEFI", [rotidx, "", "", "", pos.x(),pos.y(),pos.z()], comment)
			self.input.addCard(card)

	#-----------------------------------------------------------------------
	# FIXME densities I need a dictionary of material/density from the regions
	#-----------------------------------------------------------------------
	def __parseMaterial(self, comment, name, compound): # , regions to get densities....
		tokens = compound.split()

		def za(x): return divmod(int(x.split(".")[0]), 1000)

		elems = {}

		if len(tokens)/2 == 1:
			# isotope
			z,a = za(tokens[0])
			elems[z*1000+a] = name
			if a==0: a=""
			self.input.addCard(Card("MATERIAL", [name, z, "", "1.0", "", "", a],comment))
		else:
			# compound
			self.input.addCard(Card("MATERIAL", [name, "", "", "1.0", "", "", ""],comment))
			what = [name]
			for i in range(0,len(tokens),2):
				try:
					z,a = za(tokens[i])
				except:
					break
				frac = float(tokens[i+1])

				try:
					n = elems[z*1000+a]
				except:
					# Add a new isotope
					n = "M%d"%(z*1000+a)
					elems[z*1000+a] = n
					if a==0: a=""
					self.input.addCard(Card("MATERIAL", [n, z, "", "1.0", "", "", a]))

				what.append(frac)
				what.append(n)
			# make what multiples of 6+1
			while (len(what)-1)%6 != 0: what.append("")
			self.input.addCard(Card("COMPOUND", what))

	#-----------------------------------------------------------------------
	# Write mcnp file
	#-----------------------------------------------------------------------
	def write(self, filename):
		fout = open(filename,"w")
		fout.write("%s\n"%(self.title))
		self._writeCells(fout)
		self._writeSurfaces(fout)
		self._writeTransformations(fout)
		self._writeMaterials(fout)
		self._writeBiasing(fout)
		fout.close()
		return False

	#-----------------------------------------------------------------------
	# Write cells
	#-----------------------------------------------------------------------
	def _writeCells(self, fout):
		fout.write("C ---------- Cells ----------\n")
		for cell, comment in self.cells:
			self._writeLine(fout, cell, comment)
		fout.write("\n")

	#-----------------------------------------------------------------------
	# Write surfaces
	#-----------------------------------------------------------------------
	def _writeSurfaces(self, fout):
		fout.write("C ---------- Surfaces ----------\n")
		for surface, comment in self.surfaces:
			self._writeLine(fout, surface, comment)
		fout.write("\n")

	#-----------------------------------------------------------------------
	# Write transformations
	#-----------------------------------------------------------------------
	def _writeTransformations(self, fout):
		if len(self.transform)==0: return
		fout.write("C ---------- Transformations ----------\n")
		for trans in self.transform:
			self._writeLine(fout, trans)

	#-----------------------------------------------------------------------
	# Write materials
	#-----------------------------------------------------------------------
	def _writeMaterials(self, fout):
		fout.write("C ---------- Materials ----------\n")
		for mat in self.materials:
			if mat is None: continue
			self._writeLine(fout, "M%d %s" % (mat["@n"],
				" ".join(["%d %g"%(z,f) for z,f in mat.zaid])),
				str(mat))

	#-----------------------------------------------------------------------
	# Write Biasing
	#-----------------------------------------------------------------------
	def _writeBiasing(self, fout):
		fout.write("C ---------- Biasing ----------\n")
		line   = []
		repeat = 0
		prev   = -1.0
		for imp in self.importance:
			if imp == prev:
				repeat += 1
			else:
				if prev<0:
					pass
				else:
					line.append("%.16g"%(prev))
					if repeat == 1:
						line.append(line[-1])
					elif repeat > 1:
						line.append("%dr"%(repeat))
				prev = imp
				repeat = 0
		line.append("%.16g"%(prev))
		if repeat == 1:
			line.append(line[-1])
		else:
			line.append("%dr"%(repeat))

		self._writeLine(fout, "IMP:n %s" % (" ".join(line)))

	#-----------------------------------------------------------------------
	# Write an mcnp line, with comments
	#-----------------------------------------------------------------------
	def _writeLine(self, fout, str, comment=""):
		# write first the comments
		if comment != "":
			self._writeMultiLine(fout, comment.expandtabs(), "C --- ")

		# parse the line not to exceed 72 chars
		line = None
		for a in str.split():
			if line is None:
				line = "%-5s"%(a)
			elif len(line)+len(a) > 72:
				fout.write("%s\n"%(line))
				line = "      %s"%(a)
			else:
				line += " %s"%(a)
		fout.write("%s\n"%(line))

	#-----------------------------------------------------------------------
	def _writeMultiLine(self, fout, lines, prefix=""):
		if lines:
			for line in lines.split("\n"):
				fout.write("%s%s\n"%(prefix,line))

#===============================================================================
if __name__ == "__main__":
	Input.init()
	mcnp = Mcnp()
	mcnp.read(sys.argv[1])
	mcnp.input.write(sys.argv[2])
	sys.exit(0)

	inp = Input.Input()
	inp.verbose = True
	inp.read(sys.argv[1])
	inp.convert2Names()
	mcnp = Mcnp()
	mcnp.importFluka(inp)
	mcnp.write(sys.argv[2])
