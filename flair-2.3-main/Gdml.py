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
# Date:	21-Mar-2013

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import re
import sys
import csg
import math
import bmath
import Input
import string
import Polygon
from log import say

from xml.dom import minidom

NAMEHEX = re.compile(r"^(.*?)(0x[0-9a-f]+)$")
try:
	TRANS = string.maketrans(":","_")
except AttributeError:
	TRANS = str.maketrans(":","_")

SMALL = 1e-10

PI2   = math.pi*2.0
UNITS = {	"cm"    : 1.0,
		"mm"    : 0.1,

		"g/cm3" : 1.0,

		"deg"   : math.pi/180.0,
		"rad"   : 1.0,

		"K"     : 1.0,

		"g/mole": 1.0,

		"meV"   : 1e-12,
		"eV"    : 1e-9,
		"keV"   : 1e-6,
		"MeV"   : 1e-3,
		"GeV"   : 1.0,
		"TeV"   : 1e3,
		"PeV"   : 1e6,

		"pascal"     : 1./101325.,
		"bar"        : 1./1.01325,
		"atmosphere" : 1.0,
	}

#def unit(x): return UNITS.get(x,1.0)
def unit(x): return UNITS[x]

#-------------------------------------------------------------------------------
def uniqueName(name, dictionary, prefix):
	name = name.translate(TRANS)
	pat = NAMEHEX.match(name)
	if pat is not None:
		name8 = oname8 = pat.group(1)[:8]
	else:
		name8 = oname8 = name[:8]
	if not Input._NAMEPAT.match(name8):
		name8 = prefix
	i = 0

	while name8 in dictionary:
		i += 1
		si = str(i)
		if len(oname8) + len(si) <= 8:
			name8 = oname8 + si
		else:
			name8 = oname8[:8-len(si)]+si
	dictionary[name8] = name
	if name != name8:
		dictionary[name] = name8
	return name8

#===============================================================================
# Gdml Importer/Exporter class
#===============================================================================
class Gdml:
	#-----------------------------------------------------------------------
	# Gdml constructor
	#-----------------------------------------------------------------------
	def __init__(self, filename=None):
		self.init()
		if filename is not None:
			self.read(filename)

	#-----------------------------------------------------------------------
	def init(self):
		self.constant = Input._globalDict.copy()
		self.position = {}
		self.rotation = {}
		self.isotope  = {}
		self.element  = {}
		self.material = {}
		self.materialList = []
		self.solids   = {}
		self.volumes  = {}
		self.variable = {}
		self.world    = None

	#-----------------------------------------------------------------------
	def eval(self, x, u=1.0):
		#return eval(x, Input._globalDict, self.constant)*u
		return eval(x, self.constant, self.variable)*u

	#-----------------------------------------------------------------------
	# read Gdml input
	#-----------------------------------------------------------------------
	def read(self, filename, setup=True):
		self.filename = filename

		if Input._developer: say(">>> Gdml::read",filename)

		f = open(self.filename,"r")
		xml = minidom.parse(f)
		f.close()

		if Input._developer: say(">>> Gdml::parseDefine")
		for node in xml.getElementsByTagName("define"):
			self._parseDefinitions(node)
		if Input._developer: say(">>> Gdml::parseMaterials")
		for node in xml.getElementsByTagName("materials"):
			self._parseMaterials(node)
		if Input._developer: say(">>> Gdml::parseSolids")
		for node in xml.getElementsByTagName("solids"):
			self._parseSolids(node)
		if Input._developer: say(">>> Gdml::parseStructure")
		for node in xml.getElementsByTagName("structure"):
			self._parseStructure(node)
		if Input._developer: say(">>> Gdml::parseSetup")
		if setup: self._parseSetup(xml)

	#-----------------------------------------------------------------------
	# parse constant
	#-----------------------------------------------------------------------
	def _parseNameValue(self, node):
		name = node.attributes["name"].value
		try: value = self.eval(node.attributes["value"].value)
		except KeyError: value = 0.0
		return name, value

	#-----------------------------------------------------------------------
	# parse a value with a unit
	#-----------------------------------------------------------------------
	def _parseValue(self, node):
		try: u = unit(node.attributes["unit"].value)
		except KeyError: u = 1.0
		try:
			return self.eval(node.attributes["value"].value, u)
		except KeyError:
			return 0.0

	#-----------------------------------------------------------------------
	# parse expression
	#-----------------------------------------------------------------------
	def _parseExpression(self, node):
		name = node.attributes["name"].value
		for child in node.childNodes:
			if child.nodeType != child.TEXT_NODE: continue
			expr = child.data
			break
		try: value = self.eval(expr)
		except KeyError: value = 0.0
		return name, value

	#-----------------------------------------------------------------------
	# Parse vector
	#-----------------------------------------------------------------------
	def _parseVector(self, node, u):
		name = str(node.attributes["name"].value)

		try: u = unit(node.attributes["unit"].value)
		except KeyError: pass

		try: x = self.eval(node.attributes["x"].value, u)
		except KeyError: x = 0.0

		try: y = self.eval(node.attributes["y"].value, u)
		except KeyError: y = 0.0

		try: z = self.eval(node.attributes["z"].value, u)
		except KeyError: z = 0.0

		return name,bmath.Vector(x,y,z)

	#-----------------------------------------------------------------------
	# Parse position
	#-----------------------------------------------------------------------
	def _parsePosition(self, node):
		name, pos = self._parseVector(node, 0.1)
		self.position[name] = bmath.Matrix.translate(pos)
		return name

	#-----------------------------------------------------------------------
	# Parse rotation
	#-----------------------------------------------------------------------
	def _parseRotation(self, node):
		name, rot = self._parseVector(node, 1.0)

		if rot[2] != 0.0:
			matrix = bmath.Matrix.rotZ(-rot[2])
		else:
			matrix = None

		if rot[1] != 0.0:
			m = bmath.Matrix.rotY(-rot[1])
			if matrix is None:
				matrix = m
			else:
				matrix = m * matrix

		if rot[0] != 0.0:
			m = bmath.Matrix.rotX(-rot[0])
			if matrix is None:
				matrix = m
			else:
				matrix = m * matrix

		if matrix is None:
			matrix = bmath.Matrix(4,type=1)

		self.rotation[name] = matrix
		return name

	#-----------------------------------------------------------------------
	# Parse definitions
	#-----------------------------------------------------------------------
	def _parseDefinitions(self, node):
		for child in node.childNodes:
			if child.localName is None:
				continue

			elif child.localName == "constant":
				name, value = self._parseNameValue(child)
				self.constant[name] = value

			elif child.localName == "quantity":
				name = child.attributes["name"].value
				value = self._parseValue(child)
				self.constant[name] = value

			elif child.localName == "expression":
				name, value = self._parseExpression(child)
				self.constant[name] = value

			elif child.localName == "position":
				self._parsePosition(child)

			elif child.localName == "rotation":
				self._parseRotation(child)

			elif child.localName == "variable":
				name, value = self._parseNameValue(child)
				self.variable[name] = value

			#elif child.localName == "scale":
			#	pass
			#elif child.localName == "matrix":
			#	pass
			else:
				say("define: %s is not implemented"%(child.localName))

	#-----------------------------------------------------------------------
	# Parse materials
	#-----------------------------------------------------------------------
	def _parseMaterials(self, node):
		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "define":
				self._parseDefinitions(child)
			elif child.localName == "element":
				self._parseElement(child)
			elif child.localName == "isotope":
				self._parseIsotope(child)
			elif child.localName == "material":
				self._parseMaterial(child)
			else:
				say("materials: %s is not implemented"%(child.localName))

	#-----------------------------------------------------------------------
	def _parseIsotope(self, node):
		isotope = {}
		isotope["name"] = name = str(node.attributes["name"].value)
		isotope["Z"] = node.attributes["Z"].value
		isotope["N"] = node.attributes["N"].value

		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "atom":
				isotope["atom"] = self._parseValue(child)

		self.isotope[name] = isotope

	#-----------------------------------------------------------------------
	# Parse an element
	#-----------------------------------------------------------------------
	def _parseElement(self, node):
		element = {}
		element["name"] = name = str(node.attributes["name"].value)
		try: element["formula"] = str(node.attributes["formula"].value)
		except KeyError: element["formula"] = ""
		try: element["Z"] = int(float(node.attributes["Z"].value))
		except KeyError: element["Z"] = ""

		fraction = []

		for child in node.childNodes:
			if child.localName is None:
				continue

			elif child.localName == "atom":
				if "mass" in element: raise
				# XXX FIXME WRONG to be verified!!!
				element["mass"] = self._parseValue(child)

			elif child.localName == "fraction":
				try: n = -eval(child.attributes["n"].value)
				except KeyError: raise
				try: isotope = str(child.attributes["ref"].value)
				except KeyError: continue
				fraction.append((n,isotope))

			else:
				say("element %s: %s is not implemented"%(name, child.localName))
		element["fraction"] = fraction
		self.element[name] = element
		return element

	#-----------------------------------------------------------------------
	# Parse material
	#-----------------------------------------------------------------------
	def _parseMaterial(self, node):
		material = {}
		material["name"] = name = str(node.attributes["name"].value)
		try: material["formula"] = str(node.attributes["formula"].value)
		except KeyError: pass
		try: material["Z"] = int(float(node.attributes["Z"].value))
		except KeyError: material["Z"] = ""
		material["mass"] = ""

		fraction = []

		for child in node.childNodes:
			if child.localName is None:
				continue

			elif child.localName == "D":
				material["D"] = self._parseValue(child)

			elif child.localName == "atom":
				# XXX FIXME WRONG to be verified!!!
				material["mass"] = self._parseValue(child)

			elif child.localName == "fraction":
				try: n = -eval(child.attributes["n"].value)
				except KeyError: raise
				try: isotope = str(child.attributes["ref"].value)
				except KeyError: continue
				fraction.append((n,isotope))

			elif child.localName == "composite":
				try: n = eval(child.attributes["n"].value)
				except KeyError: raise
				try: isotope = str(child.attributes["ref"].value)
				except KeyError: continue
				fraction.append((n,isotope))

			elif child.localName == "MEE":
				material["MEE"] = self._parseValue(child)

			elif child.localName == "P":
				material["P"] = self._parseValue(child)

			elif child.localName == "T":
				material["T"] = self._parseValue(child)

			else:
				say("material %s: %s is not implemented"%(name, child.localName))

		material["fraction"] = fraction
		self.material[name] = material
		self.materialList.append(material)

	#-----------------------------------------------------------------------
	# Parse solids
	#-----------------------------------------------------------------------
	def _parseSolids(self, node):
		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "define":
				self._parseDefinitions(child)
			elif child.localName in \
				("arb8",
				 "box",
				 "cone",
				 "cutTube",
				 "ellipsoid",
				 "eltube",
				 "orb",
				 "polycone",
				 "trd",
				 "tube",
				 "sphere"):
				self._parseSolid(child)
			elif child.localName == "xtru":
				self._parseXtru(child)
			elif child.localName in ("union", "intersection", "subtraction"):
				self._parseBoolean(child)
			else:
				say("solids %s: %s is not implemented.\nPlease contact the author to add it." % \
					(child.attributes["name"].value, child.localName))

	#-----------------------------------------------------------------------
	# Parse all parameters of solid and add it to dictionary
	#-----------------------------------------------------------------------
	def _parseSolid(self, node):
		solid = {}
		solid["name"] = name = str(node.attributes["name"].value)
		solid["type"] = node.localName
		try: aunit = unit(node.attributes["aunit"].value)
		except KeyError: aunit = 1.0
		try: lunit = unit(node.attributes["lunit"].value)
		except KeyError: lunit = 0.1	# mm default unit

		for a,val in list(node.attributes.items()):
			if a in ("name","aunit","lunit"):
				continue
			elif a in ("startphi","deltaphi","starttheta","deltatheta"):
				val = self.eval(val, aunit)
			else:
				val = self.eval(val, lunit)
			solid[a] = val

		# Additional information e.g. zplane
		zplanes = []
		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "zplane":
				zplane = {}
				for a,val in list(child.attributes.items()):
					if a in ("name","aunit","lunit"):
						continue
					zplane[a] = self.eval(val, lunit)
				zplanes.append(zplane)

		if zplanes: solid["zplane"] = zplanes

		self.solids[name] = solid

	#-----------------------------------------------------------------------
	# Parse extruded objects
	#-----------------------------------------------------------------------
	def _parseXtru(self, node):
		solid = {}
		solid["name"] = name = str(node.attributes["name"].value)
		solid["type"] = node.localName
		try: lunit = unit(node.attributes["lunit"].value)
		except KeyError: lunit = 0.1	# mm default unit

		polygon = Polygon.Polygon()
		solid["sections"] = sections = []

		for child in node.childNodes:
			if child.localName is None:
				continue

			elif child.localName == "twoDimVertex":
				try: x = self.eval(child.attributes["x"].value)
				except KeyError: x = 0.0
				try: y = self.eval(child.attributes["y"].value)
				except KeyError: y = 0.0
				polygon.add(bmath.Vector(x,y,0.0))

			elif child.localName == "section":
				try: order = self.eval(child.attributes["zOrder"].value)
				except KeyError: order = 0.0
				try: z = self.eval(child.attributes["zPosition"].value)
				except KeyError: z = 0.0
				try: dx = self.eval(child.attributes["xOffset"].value)
				except KeyError: dx = 0.0
				try: dy = self.eval(child.attributes["yOffset"].value)
				except KeyError: dy = 0.0
				try: scale = self.eval(child.attributes["scalingFactor"].value)
				except KeyError: scale = 1.0
				sections.append((z, dx, dy, scale))
				if dx!=0.0 or dy!=0.0 or scale!=1.0:
					say("xtru %s: dx,dy,scale are not supported"%(node.childname))

			else:
				say("xtru %s: %s is not implemented"%(node.localName, child.localName))

#		print ">>> xtru",name
#		f = open("%s.xtru"%(name),"w")
#		for x,y,z in polygon.P:
#				print>>f, x,y,z
#		print>>f, polygon.P[0][0], polygon.P[0][1], polygon.P[0][2]
#		print>>f
#		print>>f
		sections.sort(lambda x,y: cmp(x[0],y[0]))
		solid["orientation"] = polygon.orientation()
		solid["polygon"]     = polygon.split2Convex()
		self.solids[name] = solid
#		for poly in solid["polygon"]:
#			for x,y,z in poly:
#				print>>f, x,y,z
#			print>>f, poly[0][0], poly[0][1], poly[0][2]
#			print>>f
#		f.close()

	#-----------------------------------------------------------------------
	# Parse boolean solid operation
	#-----------------------------------------------------------------------
	def _parseBoolean(self, node):
		solid = {}
		solid["name"] = name = str(node.attributes["name"].value)
		solid["type"] = node.localName
		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "first":
				solid["first"] = str(child.attributes["ref"].value)
			elif child.localName == "second":
				solid["second"] = str(child.attributes["ref"].value)

			elif child.localName == "positionref":
				solid["position"] = str(child.attributes["ref"].value)
			elif child.localName == "position":
				solid["position"] = self._parsePosition(child)
			elif child.localName == "rotationref":
				solid["rotation"] = str(child.attributes["ref"].value)
			elif child.localName == "rotation":
				solid["rotation"] = self._parseRotation(child)

			elif child.localName == "firstpositionref":
				solid["firstposition"] = str(child.attributes["ref"].value)
			elif child.localName == "firstposition":
				solid["firstposition"] = self._parsePosition(child)
			elif child.localName == "firstrotationref":
				solid["firstrotation"] = str(child.attributes["ref"].value)
			elif child.localName == "firstrotation":
				solid["firstrotation"] = self._parseRotation(child)
			else:
				say("Boolean %s: %s is not implemented"%(node.localName, child.localName))

		self.solids[name] = solid

	#-----------------------------------------------------------------------
	# Parse Structure
	#-----------------------------------------------------------------------
	def _parseStructure(self, node):
		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "assembly":
				self._parseAssembly(child)
			elif child.localName == "volume":
				self._parseVolume(child)
			else:
				say("structure %s: %s is not implemented" % \
					(child.attributes["name"].value, child.localName))

	#-----------------------------------------------------------------------
	# Parse volume and add it to dictionary
	#-----------------------------------------------------------------------
	def _parseVolume(self, node):
		volume = {}
		volume["name"] = name = str(node.attributes["name"].value)
		volume["physvol"] = []

		for child in node.childNodes:
			if child.localName is None: continue

			elif child.localName == "materialref":
				volume["material"] = str(child.attributes["ref"].value)

			elif child.localName == "solidref":
				volume["solid"] = str(child.attributes["ref"].value)

			elif child.localName == "physvol":
				volume["physvol"].append(self._parsePhysvol(child))

			elif child.localName == "divisionvol":
				volume["divisionvol"] = self._parseDivisionvol(child)

			elif child.localName == "replicavol":
				volume["replicavol"] = self._parseReplicavol(child)

			elif child.localName == "paramvol":
				volume["paramvol"] = self._parseParamvol(child)

			else:
				say("volume %s: %s is not implemented"%(node.localName, child.localName))

#		# Count immediate childs
#		volume["nchlids"] = len(volume["physvol"])
#		# ... and all family
#		for x in volume["physvol"]: print x
#		volume["nfamily"] = sum([x["nfamily"] for x in volume["physvol"]])
#		for i in ("divisionvol","replicavol","paramvol"):
#			# FIXME
#			if i in volume:
#				k = volume[n].get("number",0)
#				volume["nchlids"] += k
#				volume["nfamily"] *= k

		self.volumes[name] = volume

	#-----------------------------------------------------------------------
	def _parsePhysvol(self, node):
		# Parse child nodes
		cname  = None
		cpos   = None
		crot   = None

		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "file":
				self.read(str(child.attributes["name"].value), False)
				cname = str(child.attributes["volname"].value)
			elif child.localName == "volumeref":
				cname = str(child.attributes["ref"].value)
			elif child.localName == "positionref":
				cpos = str(child.attributes["ref"].value)
			elif child.localName == "position":
				cpos = self._parsePosition(child)
			elif child.localName == "rotationref":
				crot = str(child.attributes["ref"].value)
			elif child.localName == "rotation":
				crot = self._parseRotation(child)
			else:
				say("%s: %s is not implemented"%(node.localName, child.localName))

		if node.hasAttribute("name"):
			name = str(node.getAttribute("name"))
		else:
			name = cname

		return (name, cname, cpos, crot, None)

	#-----------------------------------------------------------------------
	def _parseReplicavol(self, node):
		try: number = int(eval(node.attributes["number"].value))
		except: number = 0
		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "volumeref":
				cname = str(child.attributes["ref"].value)
			elif child.localName == "replicate_along_axis":
				axis,width,offset = self._replicateAlongAxis(child)
			else:
				say("%s: %s is not implemented"%(node.localName, child.localName))
		return (cname, axis, number, width, offset)

	#-----------------------------------------------------------------------
	# Parse replicate_along_axis of replicavol
	#-----------------------------------------------------------------------
	def _replicateAlongAxis(self, node):
		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "direction":
				try: x = int(float(child.attributes["x"].value))
				except KeyError: x = 0
				try: y = int(float(child.attributes["y"].value))
				except KeyError: y = 0
				try: z = int(float(child.attributes["z"].value))
				except KeyError: z = 0
				if   x!=0: axis = "x"
				elif y!=0: axis = "y"
				elif z!=0: axis = "z"
				else:
					say("replicate_along_axis: axis not found")
					axis = "p"
			elif child.localName == "width":
				width = self._parseValue(child)
			elif child.localName == "offset":
				offset = self._parseValue(child)
		return (axis, width, offset)

	#-----------------------------------------------------------------------
	# Parse a parametrized volume
	#-----------------------------------------------------------------------
	def _parseParamvol(self, node):
		try: ncopies = int(float(node.attributes["ncopies"].value))
		except KeyError: ncopies = 1

		params = []
		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "volumeref":
				cname = str(child.attributes["ref"].value)
			elif child.localName == "parameterised_position_size":
				for child2 in child.childNodes:
					if child2.localName is None:
						continue
					elif child2.localName == "parameters":
						params.append(self._parseParametrizedPositionSize(child2))
					else:
						say("%s: %s is not implemented"%(child.localName))
			else:
				say("%s: %s is not implemented"%(node.localName, child.localName))

		if len(params) != ncopies:
			say("%s: ncopies=%d do not match the number of parameters"%(node.localName, ncopies))

		return (cname, params)

	#-----------------------------------------------------------------------
	def _parseParametrizedPositionSize(self, node):
		param = {}

		try: aunit = unit(node.attributes["aunit"].value)
		except KeyError: aunit = 1.0
		try: lunit = unit(node.attributes["lunit"].value)
		except KeyError: lunit = 0.1	# mm default unit

		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "position":
				for a,val in list(child.attributes.items()):
					if a in ("aunit","lunit","name"):
						continue
					elif a in ("startphi","deltaphi","starttheta","deltatheta"):
						val = self.eval(val, aunit)
					else:
						val = self.eval(val, lunit)
					param["pos.%s"%(a)] = val
			else:
				for a,val in list(child.attributes.items()):
					if a in ("aunit","lunit"):
						continue
					elif a in ("startphi","deltaphi","starttheta","deltatheta"):
						val = self.eval(val, aunit)
					else:
						val = self.eval(val, lunit)
					param[a] = val
		return param

	#-----------------------------------------------------------------------
	# Parse a division volume and convert it to similar parametrized volumes
	#-----------------------------------------------------------------------
	def _parseDivisionvol(self, node):
		try: u = unit(node.attributes["unit"].value)
		except KeyError: u = 1.0
		try: axis = node.attributes["axis"].value[1].lower()
		except: axis = "z"
		try: number = int(eval(node.attributes["number"].value))
		except: number = 0
		try: width = eval(node.attributes["width"].value, u)
		except: width = 0.0
		try: offset = eval(node.attributes["offset"].value, u)
		except: offset = 0.0
		for child in node.childNodes:
			if child.localName is None:
				continue
			elif child.localName == "volumeref":
				cname = str(child.attributes["ref"].value)
			else:
				say("%s: %s is not implemented"%(node.localName, child.localName))
		return (cname,axis,number,width,offset)

	#-----------------------------------------------------------------------
	def _parseAssembly(self, node):
		volume = {}
		volume["name"] = name = str(node.attributes["name"].value)
		volume["assembly"] = []

		for child in node.childNodes:
			if child.localName is None: continue

			elif child.localName == "physvol":
				volume["assembly"].append(self._parsePhysvol(child))

			else:
				say("%s %s: %s is not implemented"%(node.localName, name, child.localName))
		self.volumes[name] = volume

	#-----------------------------------------------------------------------
	# Parse Setup
	#-----------------------------------------------------------------------
	def _parseSetup(self, xml):
		for p in xml.getElementsByTagName("setup"):
			for child in p.childNodes:
				if child.localName is None:
					continue
				elif child.localName == "world":
					self.world = str(child.attributes["ref"].value)
				else:
					say("setup: %s is not implemented"%(child.localName))

	#-----------------------------------------------------------------------
	# Create a unique fluka name out of a solid name
	#-----------------------------------------------------------------------
	def addBodyName(self, name=""):
		self.bodyPos += 1
		self.regionPos += 1
		self.assignmatPos += 1
		name = uniqueName(name, self._bodyNames, "B")
		return name

	#-----------------------------------------------------------------------
	# Return the body name from a gdml solid
	#-----------------------------------------------------------------------
	def bodyName(self, name):
		return self._bodyNames[name]

	#-----------------------------------------------------------------------
	# Create a unique fluka name out of a volume name
	#-----------------------------------------------------------------------
	def addRegionName(self, name):
		self.regionPos += 1
		self.assignmatPos += 1
		name = uniqueName(name, self._regionNames, "R")
		return name

	#-----------------------------------------------------------------------
	# Return the region name from a gdml volume
	#-----------------------------------------------------------------------
	def regionName(self, name):
		return self._regionNames[name]

	#-----------------------------------------------------------------------
	# Add a body to input with comment and matrix transformation
	#-----------------------------------------------------------------------
	def _addBody(self, inp, tag, name, what, matrix, comment=""):
		bname = self.addBodyName(name)
		what.insert(0, bname)
		card = Input.Card(tag,what,comment)
		inp.addCard(card, self.bodyPos, False)
		if matrix: inp.transformBody(card, matrix)
#		what.pop(0)
		return bname

	#-----------------------------------------------------------------------
	def _writeFlukaZCutPlanes(self, inp, solid, matrix, comment=""):
		z = solid.get("z",0.0)
		# Low cutting planes
		lowX = solid.get("lowX",0.0)
		lowY = solid.get("lowY",0.0)
		lowZ = solid.get("lowZ",0.0)
		bname1 = self._addBody(inp, "PLA", solid["name"],
				[lowX, lowY, lowZ, 0.0, 0.0, -z], matrix, comment)
		# High cutting planes
		highX = solid.get("highX",0.0)
		highY = solid.get("highY",0.0)
		highZ = solid.get("highZ",0.0)
		bname2 = self._addBody(inp, "PLA",solid["name"],
				[highX, highY, highZ, 0.0, 0.0, z], matrix, comment)
		return "+%s +%s"%(bname1,bname2)

	#-----------------------------------------------------------------------
	def _writeFlukaPhiCutPlanes(self, inp, solid, matrix, comment=""):
		# Cutting planes
		sphi = solid.get("startphi",0.0)
		dphi = solid.get("deltaphi",PI2)
		if sphi>0.0 or dphi<PI2:
			c = math.cos(sphi)
			s = math.sin(sphi)
			bname1 = self._addBody(inp, "PLA", solid["name"],
					[-s, c, 0.0, 0.0, 0.0, 0.0], matrix, comment)

			f = sphi+dphi
			c = math.cos(f)
			s = math.sin(f)
			bname2 = self._addBody(inp, "PLA", solid["name"],
					[-s, c, 0.0, 0.0, 0.0, 0.0],matrix, comment)
			return "+%s -%s"%(bname2,bname1)
		else:
			return ""

	#-----------------------------------------------------------------------
	def _writeFlukaArb8(self, inp, solid, matrix, comment):
		# Find points
		p = []
		dz = solid["dz"]
		for i in range(1,9):
			x = solid.get("v%dx"%(i),0.0)
			y = solid.get("v%dy"%(i),0.0)
			if i<=4:
				z= -dz
			else:
				z = dz
			p.append(bmath.Vector(x,y,z))

		# Top plane
		same = 0
		for i in range(4,7):
			for j in range (i+1,8):
				if (p[i]-p[j]).length2() < SMALL:
					same += 1
		if same<2:
			bname = self._addBody(inp, "XYP", solid["name"],
				[dz], matrix, comment)
			expr = "+%s"%(bname)
		else:
			expr = ""

		# Bottom plane
		same = 0
		for i in range(0,3):
			for j in range (i+1,4):
				if (p[i]-p[j]).length2() < SMALL:
					same += 1

		if same<2:
			bname = self._addBody(inp, "XYP", solid["name"],
					[-dz], matrix, comment)
			expr += "-%s"%(bname)

		# Add a face with two planes abc and acd
		# FIXME order! + or - or I have to find runtime?
		def addFace(a,b,c,d,sign):
			x,y,z = p[a]

			n1 = (p[b]-p[a]) ^ (p[c]-p[a])
			n1.normalize()

			if abs(n1[0]) >= SMALL or \
			   abs(n1[1]) >= SMALL or \
			   abs(n1[2]) >= SMALL:
				bname = self._addBody(inp, "PLA", solid["name"],
					[n1[0],n1[1],n1[2],x,y,z], matrix)
				e = "%s%s "%(sign,bname)
			else:
				e = ""

			n2 = (p[c]-p[a]) ^ (p[d]-p[a])
			n2.normalize()

			if (abs(n2[0]) >= SMALL or \
			    abs(n2[1]) >= SMALL or \
			    abs(n2[2]) >= SMALL) and \
			    n1*n2 < 0.999999:
				say("WARNING: twisted arb8 %s is approximated (%g)" \
					%(solid["name"],n1*n2))
				bname = self._addBody(inp, "PLA", solid["name"],
					[n2[0],n2[1],n2[2],x,y,z], matrix)
				e += "%s%s "%(sign,bname)
			return e

		# Check sign of operation
		A = p[1]-p[0]
		B = p[3]-p[0]
		if abs(A*B) < SMALL:
			A = p[5]-p[4]
			B = p[7]-p[4]

		if (A^B)[2] > 0.0:
			sign = "+"
		else:
			sign = "-"

		expr += addFace(0,1,5,4,sign)
		expr += addFace(1,2,6,5,sign)
		expr += addFace(2,3,7,6,sign)
		expr += addFace(3,0,4,7,sign)

		return expr

	#-----------------------------------------------------------------------
	def _writeFlukaBox(self, inp, solid, matrix, comment):
		return "+%s"%(self._addBody(inp, "RPP", solid["name"],
				[-solid["x"]/2.0, solid["x"]/2.0,
				 -solid["y"]/2.0, solid["y"]/2.0,
				 -solid["z"]/2.0, solid["z"]/2.0], matrix, comment))

	#-----------------------------------------------------------------------
	def _writeFlukaCone(self, inp, solid, matrix, comment):
		# outer cone
		bname = self._addBody(inp, "TRC", solid["name"],
				[0.0, 0.0, -solid["z"]/2.0,
				 0.0, 0.0,  solid["z"],
				 solid["rmax1"],solid["rmax2"]], matrix, comment)
		expr = "+%s "%(bname)

		if solid.get("rmin1",0.0)>0.0 or solid.get("rmin2",0.0)>0.0:
			# Inner cone
			bname = self._addBody(inp, "TRC", solid["name"],
					[0.0, 0.0, -solid["z"]/2.0,
					 0.0, 0.0,  solid["z"],
					 solid["rmin1"],solid["rmin2"]], matrix, comment)
			expr += "-%s "%(bname)
		expr += self._writeFlukaPhiCutPlanes(inp, solid, matrix)
		return expr

	#-----------------------------------------------------------------------
	def _writeFlukaCutTube(self, inp, solid, matrix, comment):
		# outer cylinder
		bname =self._addBody(inp,  "ZCC", solid["name"],
				[0.0, 0.0, solid["rmax"]], matrix, comment)
		expr = "+%s "%(bname)

		if solid.get("rmin",0.0) > 0.0:
			# Inner cylinder
			bname = self._addBody(inp, "ZCC", solid["name"],
					[0.0, 0.0, solid["rmin"]], matrix, comment)
			expr += "-%s "%(bname)

		expr += self._writeFlukaPhiCutPlanes(inp, solid, matrix, comment)
		expr += self._writeFlukaZCutPlanes(inp, solid, matrix, comment)
		return expr

	#-----------------------------------------------------------------------
	def _writeFlukaElTube(self, inp, solid, matrix, comment):
		bname = self._addBody(inp, "REC", solid["name"],
				[0.0, 0.0, -solid["dz"]/2.0,
				 0.0, 0.0,  solid["dz"],
				 solid["dx"], 0.0, 0.0,
				 0.0, solid["dy"], 0.0], matrix, comment)
		return "+%s"%(bname)

	#-----------------------------------------------------------------------
	def _writeFlukaEllipsoid(self, inp, solid, matrix, comment):
		bname = self._addBody(inp, "QUA", solid["name"],
				[1.0/solid["ax"]**2, 1.0/solid["by"]**2, 1.0/solid["cz"]**2,
				 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, -1.0],matrix,comment)
		expr = "+%s"%(bname)

		if "zcut1" in solid:
			bname = self._addBody(inp, "XYP", solid["name"],
				[solid["zcut1"]], matrix, comment)
			expr += "-%s"%(bname)

		if "zcut2" in solid:
			bname = self._addBody(inp, "XYP", solid["name"],
				[solid["zcut2"]], matrix, comment)
			expr += "+%s"%(bname)

		return expr

	#-----------------------------------------------------------------------
	def _writeFlukaOrb(self, inp, solid, matrix, comment):
		bname = self._addBody(inp, "SPH", solid["name"],
				[0.0, 0.0, 0.0, solid["r"]], matrix, comment)
		return "+%s"%(bname)

	#-----------------------------------------------------------------------
	def _writeFlukaPolycone(self, inp, solid, matrix, comment):
		# remember first z-section
		zplane = solid["zplane"][0]
		z    = zplane.get("z",0.0)	# z-location
		rmin = zplane.get("rmin",0.0)	# inner radius
		rmax = zplane.get("rmax",0.0)	# outer radius
		if rmin > rmax:
			say("ERROR: polycone %s plane z=%g rmin>rmax  rmin=%g rmax=%g"% \
						(solid["name"], zplane["z"], rmin, rmax))
			rmin, rmax = rmax, rmin

		# First plane
		zname = self._addBody(inp, "XYP", solid["name"],
				[z], matrix, comment)
		expr = ""
		phicut = self._writeFlukaPhiCutPlanes(inp, solid, matrix, comment)

		for zplane in solid["zplane"][1:]:
			# write something if there is a change in z
			zn    = zplane.get("z",0.0)	# next z plane
			rnmin = zplane.get("rmin",0.0)	# next inner radius
			rnmax = zplane.get("rmax",0.0)	# next outer radius
			if rnmin > rnmax:
				say("ERROR: polycone %s plane z=%g rmin>rmax  rmin=%g rmax=%g"%\
						(solid["name"], zn, rnmin, rnmax))
				rnmin, rnmax = rnmax, rnmin

			# WARNING polycone is the G4GenericPolycone which can decrease as well
			if zn < z:
				swap = True
				zn, z = z, zn
				rnmin, rmin = rmin, rnmin
				rnmax, rmax = rmax, rnmax
			else:
				swap = False

			if abs(zn-z) >= SMALL:
				if expr!="": expr += " | "
				if swap:
					bname = self._addBody(inp, "XYP", solid["name"],
							[z], matrix, comment)
					expr += "+%s -%s"%(zname,bname)
				else:
					bname = self._addBody(inp, "XYP", solid["name"],
							[zn], matrix, comment)
					expr += "+%s -%s"%(bname,zname)
				zname = bname	# remember last z-plane

				# Check inner radius
				if abs(rnmin-rmin) < SMALL:
					# Cylinder
					if rnmin > SMALL:
						# with non-zero radius
						bname = self._addBody(inp, "ZCC", solid["name"],
								[0.,0.,rmin], matrix, comment)
						expr += "-%s "%(bname)
				else:
					# Cone
					Ra = rnmin
					Rb = rmin
					za = zn
					zb = z

					zc = (Rb*za-Ra*zb)/(Rb-Ra)	# zero point
					s  = (Rb-Ra)/(zb-za)		# slope
					s2 = s**2
					bname = self._addBody(inp, "QUA", solid["name"],
							[1., 1., -s2,
							 0., 0., 0.,
							 0., 0., 2.*s2*zc,
							 -s2*zc*zc], matrix, comment)
					expr += "-%s "%(bname)

				# Check inner radius
				if abs(rnmax-rmax) < SMALL:
					# Cylinder
					bname = self._addBody(inp, "ZCC", solid["name"],
							[0.,0.,rmax], matrix, comment)
				else:
					# Cone
					Ra = rnmax
					Rb = rmax
					za = zn
					zb = z

					zc = (Rb*za-Ra*zb)/(Rb-Ra)	# zero point
					s  = (Rb-Ra)/(zb-za)		# slope
					s2 = s**2
					bname = self._addBody(inp, "QUA", solid["name"],
							[1., 1., -s2,
							 0., 0., 0.,
							 0., 0., 2.*s2*zc,
							 -s2*zc*zc], matrix, comment)
					#bname = self._addBody(inp, "TRC", solid["name"],
					#		[0., 0., zb,
					#		 0., 0., za-zb,
					#		 Rb, Ra], matrix, comment)
				expr += "+%s%s"%(bname,phicut)

			# remember last z-section
			if not swap:
				z    = zn	# z-location
				rmin = rnmin	# inner radius
				rmax = rnmax	# outer radius
		return expr

	#-----------------------------------------------------------------------
	def _writeFlukaSphere(self, inp, solid, matrix, comment):
		bname = self._addBody(inp, "SPH", solid["name"],
				[0.0, 0.0, 0.0, solid["rmax"]], matrix, comment)
		expr = "+%s"%(bname)

		if "rmin" in solid:
			bname = self._addBody(inp, "SPH", solid["name"],
					[0.0, 0.0, 0.0, solid["rmin"]], matrix, comment)
			expr += "-%s"%(bname)

		expr += self._writeFlukaPhiCutPlanes(inp, solid, matrix, comment)

		# Theta Cutting planes
		stheta = solid.get("starttheta",0.0)
		dtheta = solid.get("deltatheta",PI2)
		if stheta > 0.0 or dtheta<PI2:
			c = math.cos(stheta)
			s = math.sin(stheta)
			bname1 = self._addBody(inp, "PLA", solid["name"],
					[-c, s, 0.0, 0.0, 0.0, 0.0], matrix, comment)

			f = stheta+dtheta
			c = math.cos(f)
			s = math.sin(f)
			bname2 = self._addBody(inp, "PLA", solid["name"],
					[-c, s, 0.0, 0.0, 0.0, 0.0], matrix, comment)
			expr += "+%s -%s"%(bname2,bname1)
		else:
			return ""
		return expr

	#-----------------------------------------------------------------------
	def _writeFlukaTube(self, inp, solid, matrix, comment):
		rmin = solid.get("rmin",0.0)
		rmax = solid.get("rmax",0.0)
		if rmin > rmax:
			say("ERROR: tube %s rmin>rmax  rmin=%g rmax=%g"%(solid["name"], rmin, rmax))
			rmin, rmax = rmax, rmin

		# outer cylinder
		bname = self._addBody(inp, "RCC", solid["name"],
				[0.0, 0.0, -solid["z"]/2.0,
				 0.0, 0.0,  solid["z"],
				 rmax], matrix, comment)
		expr = "+%s "%(bname)

		if rmin > 0.0:
			# Inner cylinder
			bname = self._addBody(inp, "ZCC", solid["name"],
					[0.0, 0.0, rmin], matrix, comment)
			expr += "-%s "%(bname)

		expr += self._writeFlukaPhiCutPlanes(inp, solid, matrix, comment)
		return expr

	#-----------------------------------------------------------------------
	def _writeFlukaTrd(self, inp, solid, matrix, comment):
		z  = solid["z"] /2
		x1 = solid["x1"]/2
		y1 = solid["y1"]/2
		x2 = solid["x2"]/2
		y2 = solid["y2"]/2

		# cutting z-planes
		bname = self._addBody(inp, "XYP", solid["name"],
				[ z ], matrix, comment)
		expr = "+%s "%(bname)

		bname = self._addBody(inp, "XYP", solid["name"],
				[-z], matrix, comment)
		expr += "-%s "%(bname)

		# -X
		bname = self._addBody(inp, "PLA", solid["name"],
				[-2.0*z, 0., (x1-x2), -x1, -y1, -z],
				matrix, comment)
		expr += "+%s "%(bname)

		# -Y
		bname = self._addBody(inp, "PLA", solid["name"],
				[0., -2.0*z, (y1-y2), -x1, -y1, -z],
				matrix, comment)
		expr += "+%s "%(bname)

		# +X
		bname = self._addBody(inp, "PLA", solid["name"],
				[2.0*z, 0., (x1-x2), x1, y1, -z],
				matrix, comment)
		expr += "+%s "%(bname)

		# +Y
		bname = self._addBody(inp, "PLA", solid["name"],
				[0., 2.0*z, (y1-y2), x1, y1, -z],
				matrix, comment)
		expr += "+%s "%(bname)

		return expr

	#-----------------------------------------------------------------------
	# FIXME VERY PRIMITIVE!!!!
	#-----------------------------------------------------------------------
	def _writeFlukaXtru(self, inp, solid, matrix, comment):
		# Write Z plane
		s = "-"
		zexpr = ""
		for z,dx,dy,scale in solid["sections"]:
			bname = self._addBody(inp, "XYP", solid["name"],
				[z], matrix, comment)
			zexpr += "%s%s "%(s,bname)
			s = "+"

		# Write polygones
		Z = bmath.Vector(0.0,0.0,solid["orientation"])
		bname = None
		expr = ""
		lastPoint = None
		for P in solid["polygon"]:
			zone = zexpr
			for i in range(len(P)):
				j = i+1
				if j>=len(P): j=0

				# Check last polygon-plane with current one
				if lastPoint and (lastPoint[0]-P[j]).length()<SMALL:
					zone += "-%s "%(lastPoint[1])
					lastPoint = None
					continue

				N = (P[j]-P[i]) ^ Z
				N.normalize()

				bname = self._addBody(inp, "PLA", solid["name"],
					[N[0],N[1],N[2], P[i][0],P[i][1],P[i][2]], matrix, comment)
				last = (P[i],bname)
				zone += "+%s "%(bname)
			lastPoint = last

			if expr:
				expr += " | %s"%(zone)
			else:
				expr = zone

		return expr

	#-----------------------------------------------------------------------
	def _writeFlukaBoolean(self, inp, solid, matrix, comment):
		comment = "%s [%s]"%(solid["name"],comment)

		# first body
		if "firstposition" in solid or "firstrotation" in solid:
			cmatrix = self.transformMatrix(matrix,
						solid.get("firstposition"),
						solid.get("firstrotation"),
						True)
		else:
			cmatrix = matrix
		first = self.writeFlukaSolid(inp, solid["first"], cmatrix, comment+" Boolean: First")

		# second body
		cmatrix = self.transformMatrix(matrix,
					solid.get("position"),
					solid.get("rotation"),
					True)
		second = self.writeFlukaSolid(inp, solid["second"], cmatrix, comment+" Boolean: Second")

		if solid["type"] == "union":
			return "%s | %s"%(self.parenthesis(first),self.parenthesis(second))
		elif solid["type"] == "intersection":
			return "%s + %s"%(self.parenthesis(first),self.parenthesis(second))
		elif solid["type"] == "subtraction":
			return "%s - %s"%(self.parenthesis(first),self.parenthesis(second))

	#-----------------------------------------------------------------------
	# Write a SOLID into FLUKA format
	#-----------------------------------------------------------------------
	def writeFlukaSolid(self, inp, name, matrix=None, path="", param=None):
		try:
			solid = self.solids[name]
		except KeyError:
			say("ERROR: Solid %s not found"%(name))
			return ""

		# Parametrized solid, create a new dictionary
		if param is not None:
			solid = solid.copy()
			solid.update(param)

		t = solid["type"]
		comment = "%s [%s]"%(name,path)

		if t == "arb8":
			expr = self._writeFlukaArb8(inp, solid, matrix, comment)

		elif t == "box":
			expr = self._writeFlukaBox(inp, solid, matrix, comment)

		elif t == "cone":
			expr = self._writeFlukaCone(inp, solid, matrix, comment)

		elif t == "cutTube":
			expr = self._writeFlukaCutTube(inp, solid, matrix, comment)

		elif t == "eltube":
			expr = self._writeFlukaElTube(inp, solid, matrix, comment)

		elif t == "ellipsoid":
			expr = self._writeFlukaEllipsoid(inp, solid, matrix, comment)

		elif t == "orb":
			expr = self._writeFlukaOrb(inp, solid, matrix, comment)

		elif t == "polycone":
			expr = self._writeFlukaPolycone(inp, solid, matrix, comment)

		elif t == "sphere":
			expr = self._writeFlukaSphere(inp, solid, matrix, comment)

		elif t == "trd":
			expr = self._writeFlukaTrd(inp, solid, matrix, comment)

		elif t == "tube":
			expr = self._writeFlukaTube(inp, solid, matrix, comment)

		elif t == "xtru":
			expr = self._writeFlukaXtru(inp, solid, matrix, comment)

		elif t in ("union", "intersection", "subtraction"):
			expr = self._writeFlukaBoolean(inp, solid, matrix, comment)

		else:
			return ""

		return expr.strip()

	#-----------------------------------------------------------------------
	# Transform matrix with rot,pos
	# WARNING: be careful on operations to return a NEW matrix and not
	#          to modify the original
	#-----------------------------------------------------------------------
	def transformMatrix(self, matrix, pos, rot, inverse=False):
		cmatrix = matrix

		if pos is not None:
			if isinstance(pos,str):
				m = self.position[pos]
			else:
				m = pos
			if cmatrix is None:
				cmatrix = m.clone()
			else:
				cmatrix = cmatrix * m

		if rot is not None:
			if isinstance(rot,str):
				m = self.rotation[rot]
			else:
				m = rot

			if inverse:
				m = m.clone()
				m.inverse()		# WARNING inverse in place
				if cmatrix is None:
					cmatrix = m
				else:
					cmatrix = cmatrix * m

			elif cmatrix is None:
				cmatrix = m.clone()
			else:
				cmatrix = cmatrix * m

		return cmatrix

	#-----------------------------------------------------------------------
	# Use or strip the parentheses of a sub-expression
	#-----------------------------------------------------------------------
	def parenthesis(self, expr):
		if expr.find(" ") >= 0:
			return "(%s)"%(expr)
		else:
			return expr[1:]

	#-----------------------------------------------------------------------
	# Convert replicavol to physvol
	#-----------------------------------------------------------------------
	def replicavol2physvol(self, name, volume):
		solid  = self.solids[volume["solid"]]
		t = solid["type"]

		cname, axis, number, width, offset = volume["replicavol"]

		def calcDimensions(low,high,number,width,offset):
			# Calculate dimensions
			low += offset
			if width>0.0:
				n = int((high-low)/width)
				if number == 0: number = n
			else:
				width = (high-low)/number
			return low,high,number,width

		if axis in ("x","y","z"):
			h = solid[axis]/2.0
			low,high,number,width = calcDimensions(-h,h,number,width,offset)
			pos = bmath.Vector()
			idx = ord(axis) - ord('x')
			for i in range(number):
				pos[idx] = low + (i+0.5)*width
				m = bmath.Matrix.translate(pos)
				volume["physvol"].append((name, cname, m, None, None))
		elif axis=="p":
			l = solid.get("startphi",0.0)
			h = l + solid.get("deltaphi",PI2)
			low,high,number,width = calcDimensions(l,h,number,width,offset)
			for i in range(number):
				rot = low + (i+0.5)*width
				m = bmath.Matrix.rotZ(rot)
				volume["physvol"].append((name, cname, None, m, None))

		elif axis=="r":
			rmin = solid["rmin"]
			rmax = solid["rmax"]
			low,high,number,width = calcDimensions(rmin,rmax,number,width,offset)
			for i in range(number):
				param = {"rmin" : low + i*width,
					 "rmax" : low + (i+1)*width}
				volume["physvol"].append((name, cname, None, None, param.copy()))

		else:
			say("ERROR replicavol %s: solid %s type=%s axis=%s not implemented"%(name,solid,t,axis))

		del volume["replicavol"]	# delete not to reprocess next time

	#-----------------------------------------------------------------------
	# Convert paramvol to physvol
	#-----------------------------------------------------------------------
	def paramvol2physvol(self, name, volume):
		solid  = self.solids[volume["solid"]]
		t = solid["type"]
		cname, params = volume["paramvol"]
		for p in params:
			tx = p.get("pos.x",0.0)
			ty = p.get("pos.y",0.0)
			tz = p.get("pos.z",0.0)
			if tx!=0.0 or ty!=0.0 or tz!=0.0:
				m = bmath.Matrix.translate(tx,ty,tz)
			else:
				m = None
			volume["physvol"].append((name, cname, m, None, p))
		del volume["paramvol"]	# delete it so not to reprocess next time

	#-----------------------------------------------------------------------
	# Convert divisionvol to physvol
	#
	# List of volumes that support divisioning
	# G4Box		kXAxis, kYAxis, kZAxis
	# G4Tubs	kRho,   kPhi,   kZAxis
	# G4Cons	kRho,   kPhi,   kZAxis
	# G4Trd		kXAxis, kYAxis, kZAxis
	# G4Para	kXAxis, kYAxis, kZAxis
	# G4Polycone	kRho,   kPhi,   kZAxis	FIXME not implemented
	# G4Polyhedra	kRho,   kPhi,   kZAxis (*)
	#-----------------------------------------------------------------------
	def divisionvol2physvol(self, name, volume):
		solid  = self.solids[volume["solid"]]
		t = solid["type"]
		param = solid.copy()
		del param["name"]
		del param["type"]

		if t in ("polycone","polyhedra"):
			say("ERROR divisionvol %s: solid %s type=%s axis=%s not implemented"%(name,solid,t,axis))

		cname, axis, number, width, offset = volume["divisionvol"]

		def calcDimensions(low,high,number,width,offset):
			# Calculate dimensions
			low += offset
			if width>0.0:
				n = int((high-low)/width)
				if number == 0: number = n
			else:
				width = (high-low)/number
			return low,high,number,width

		if axis in ("x","y","z"):
			h = solid[axis]/2.0
			low,high,number,width = calcDimensions(-h,h,number,width,offset)
			param[axis] = width
			pos = bmath.Vector()
			idx = ord(axis) - ord('x')
			for i in range(number):
				pos[idx] = low + (i+0.5)*width
				m = bmath.Matrix.translate(pos)
				volume["physvol"].append((name, cname, m, None, param.copy()))
		elif axis=="p":
			l = solid.get("startphi",0.0)
			h = l + solid.get("deltaphi",PI2)
			low,high,number,width = calcDimensions(l,h,number,width,offset)
			param["deltaphi"] = width
			for i in range(number):
				param["startphi"] = low + i*width
				volume["physvol"].append((name, cname, None, None, param.copy()))
		elif axis=="r":
			rmin = solid.get("rmin",0.0)
			rmax = solid.get("rmax",0.0)
			low,high,number,width = calcDimensions(rmin,rmax,number,width,offset)
			for i in range(number):
				param["rmin"] = low + i*width
				param["rmax"] = low + (i+1)*width
				volume["physvol"].append((name, cname, None, None, param.copy()))
		else:
			say("ERROR divisionvol %s: solid %s type=%s axis=%s not implemented"%(name,solid,t,axis))

		del volume["divisionvol"]	# delete not to reprocess next time

	#-----------------------------------------------------------------------
	# Write FLUKA volume
	#-----------------------------------------------------------------------
	def writeFlukaVolume(self, inp, name, vname, depth, matrix=None, path="", param=None):
#		if Input._developer: say("*-* writeFlukaVolume (depth=%d) %s %s"%(depth,name,str(param)))
		try:
			volume = self.volumes[vname]
		except KeyError:
			say("Volume %s not found"%(vname))
			return ""

		newpath = path + "/" + vname

		# Write necessary bodies
		mother = self.writeFlukaSolid(inp, volume["solid"], matrix, newpath, param)
		if mother == "": return ""

		# Create division volumes if any
		if "replicavol" in volume:
			self.replicavol2physvol(vname, volume)
		elif "divisionvol" in volume:
			self.divisionvol2physvol(vname, volume)
		elif "paramvol" in volume:
			self.paramvol2physvol(vname, volume)

		# Write child bodies
		child = ""
		for pname, cname, cpos, crot, cparam in volume["physvol"]:
			cmatrix = self.transformMatrix(matrix, cpos, crot)
			cvol = self.volumes.get(cname,{})
			if "assembly" in cvol:
				# Expand assembly
				apath = newpath + "/" + cname
				for aname, cname, apos, arot, cparam in cvol["assembly"]:
					amatrix = self.transformMatrix(cmatrix, apos, arot)
					m = self.writeFlukaVolume(inp, aname, aname, depth+1, amatrix, apath, cparam)
					if m: child += "-%s "%(self.parenthesis(m))
			else:
				# Insert physvol
				m = self.writeFlukaVolume(inp, pname, cname, depth+1, cmatrix, newpath, cparam)
				if m: child += "-%s "%(self.parenthesis(m))

		# Write expression
		rname = self.addRegionName(name)
		if child!="" and "|" in mother: mother = "(%s)"%(mother)
		expr = Input.Card.splitExpr("%s %s"%(mother,child),80)
		card = Input.Card("REGION",[rname],"%s [%s]\nref: %s"%(name,path,vname),expr)
		inp.addCard(card, self.regionPos, False)

		# Write material...
		try:
			mat = self.material[volume["material"]]["fluka"]
		except:
			mat = volume.get("material","VACUUM")
		inp.addCard(Input.Card("ASSIGNMA",["",mat,rname]),
				self.assignmatPos, False)
		self.assignmatPos += 1

		return mother

	#-----------------------------------------------------------------------
	# Write everything outside world as the blackhole
	#-----------------------------------------------------------------------
	def writeBlackVolume(self, inp, world):
		expr = "-(%s)"%(world)
		name = self.addRegionName("BLKBODY")
		card = Input.Card("REGION", [name],
				"Blackbody everything outside world: %s"%(name),
				expr)
		inp.addCard(card, self.regionPos, False)
		# Write material...
		inp.addCard(Input.Card("ASSIGNMA",["","BLCKHOLE",name]),
				self.assignmatPos, False)
		self.assignmatPos += 1

	#-----------------------------------------------------------------------
	# Write Materials in FLUKA format
	#-----------------------------------------------------------------------
	def writeFlukaMaterials(self, inp):
		flukaNames = {}
		dictionary = {}

		pos = inp["GEOEND"][0].pos()

		# Write isotopes
		for name,isotope in list(self.isotope.items()):
			name8 = uniqueName(name, dictionary, "MAT")
			Z = int(float(isotope["Z"]))
			A = Z + int(float(isotope["N"]))
			inp.addCard(Input.Card("MATERIAL",
					[name8, Z, "", 1.0, "", "", A]),
					pos, False)
			pos += 1

		# Write elements
		for name,element in list(self.element.items()):
			name8 = uniqueName(name, dictionary, "MAT")
			inp.addCard(Input.Card("MATERIAL",
					[name8, element.get("Z",""),
						element.get("mass",""),
						1.0],
					element.get("formula","")),
				pos, False)
			pos += 1

			whats = [name8]
			for n,iso in element["fraction"]:
				whats.append(bmath.format(n))
				whats.append(dictionary[iso])
			if len(whats)>1:
				inp.addCard(Input.Card("COMPOUND",whats), pos, False)
			pos += 1

		# Write materials
		for mat in self.materialList:
		#for name,mat in self.material.items():
			name8 = uniqueName(mat["name"], dictionary, "MAT")
			mat["fluka"] = name8

			inp.addCard(Input.Card("MATERIAL",
					[name8, mat.get("Z",""),
						mat.get("mass",""),
						mat.get("D","")]),
					pos, False)
			pos += 1

			P   = mat.get("P","")
			mee = mat.get("MEE","")
			if P or mee:
				try: mee *= 1e9	# convert to eV
				except: pass
				inp.addCard(Input.Card("MAT-PROP",
						["", P, "", mee, name8]),
						pos, False)
				pos += 1

			whats = [name8]
			for n,iso in mat["fraction"]:
				whats.append(bmath.format(n))
				try:
					# find component name
					m = dictionary[iso]
				except KeyError:
					# if it is defined later then add it first
					m = uniqueName(iso, dictionary, "MAT")
				whats.append(m)
			card = inp.addCard(Input.Card("COMPOUND",whats), pos, False)
			card.padWhats()
			pos += 1

	#-----------------------------------------------------------------------
	# Convert to FLUKA format
	#-----------------------------------------------------------------------
	def toFluka(self, inp):
		if not inp.cardlist: inp.minimumInput()

		# Initialize fluka pos
		self._bodyNames   = {}	# Dictionary of used named
		self._regionNames = {}	# Dictionary of used named
		self.bodyPos      = inp["END"][0].pos()-1	# Find first end position
		self.regionPos    = inp["END"][1].pos()-1
		self.assignmatPos = inp["GEOEND"][0].pos()+1

		# Write first the materials
		self.writeFlukaMaterials(inp)
		# Start with world
		world = self.writeFlukaVolume(inp, self.world, self.world, 0)
		self.writeBlackVolume(inp, world)
		inp.renumber()

		del self._regionNames
		del self._bodyNames

		return inp

	#-----------------------------------------------------------------------
	# Write Gdml file
	#-----------------------------------------------------------------------
	def write(self):
		pass

#===============================================================================
if __name__ == "__main__":
	Input.init()
	gdml = Gdml(sys.argv[1])

	inp = Input.Input()
	gdml.toFluka(inp)
	inp.write(sys.argv[2])

	sys.exit(0)
