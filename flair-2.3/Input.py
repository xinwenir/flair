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
# Date:	17-May-2004

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import io
import os
import re
import csg
import sys
import time
import math
import struct
import string
from log import say
from types import *
from legacy import *
from operator import attrgetter

try:
	import configparser
except ImportError:
	import configparser as ConfigParser

try:
	import pickle as pickle
except ImportError:
	import pickle

import bmath
import fortran

__first = True

#-------------------------------------------------------------------------------
# Activate untested developer code
_developer = False
_useQUA    = True
_useBOX    = False	# enable the use of BOX/ARB/WED
_database  = "db/card.ini"

_NAMEPAT   = re.compile(r"^[A-Za-z_][A-Za-z0-9_.:!\$]*$")
_REGIONPAT = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_.:!\$]*)\s*(-?\d+)\s*(.*)$")
_VOXELPAT  = re.compile(r"^VOX[E]?[L]?(\d+)$")
_CWHAT     = re.compile(r"^[*!]@what\.(-?\d+)(=.*)$")
_OPT       = re.compile(r"^!@(\w+)=(.*)$")
_wPAT      = re.compile(r"\bw\(")
_WPAT      = re.compile(r"\bW\(")
_WNPAT     = re.compile(r"\b([wW])\((-?\d+)\)")
_ROTPAT    = re.compile(r"^ROT?#(\d+)$",re.I)
#_PRPPAT    = re.compile(r"^(#\w+)\b *(\w*)[ ,/;\\:] *(.*)$")
_PRPPAT    = re.compile(r"^(#\w+)\b *([A-Za-z0-9_.@-]*) *(.*)$")

#-------------------------------------------------------------------------------
# Configuration
commentedCards = True		# Parse comments for disabled cards
warnMultMat    = False		# Warn user on multiple materials

# transformation values
zero     = 1.0E-10
infinite = 1.0E10

#-------------------------------------------------------------------------------
# Global variables
_defaultMaterials = []		# Default material list (cards)
_defaultMatDict   = {}
_icruMaterials    = []		# ICRU special compounds (cards)
_icruMatDict      = {}
_neutronGroups    = {}		# Neutron energy groups as numbers
_neutronGroupsS   = {}		# Neutron energy groups as strings with units
_lowMaterials     = {}		# LOW-NEUT material cross sections
_usermvax         = {}		# Description of user routines

# File formats
FORMAT_SINGLE  = 0
FORMAT_DOUBLE  = 1
FORMAT_FREE    = 2
FORMAT_VOXEL   = 3

BODY_PREFIX    = "B"
REGION_PREFIX  = "R"

SCALE="*...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....\n"
#SCALE="*tag-----><----1---><---2----><---3----><---4----><---5----><---6----><--sdum--\n"
REGSCALE="*_..._____.._____.._____.._____.._____.._____.._____.._____.._____.._____\n"

BODY_TAGS      = None	# Including VOXELS
BODY_NOVXL_TAGS= None	# Without VOXELS
FLAIR_TAGS     = None
OBJECT_TAGS    = None
TRANSFORM_TAGS = None
DETECTOR_TAGS  = ("DETECT", "EVENTBIN","RESNUCLE","USRBDX", \
		  "USRBIN", "USRCOLL", "USRTRACK","USRYIELD")
NGROUPS        = []

# Preprocessor cards that increase or decrease the indent
_INDENT_INC    = ("#if", "#ifdef", "#ifndef", "#elif", "#else", "#include",
		  "$start_expansion", "$start_translat", "$start_transform")
_INDENT_DEC    = ("#elif", "#else", "#endif", "#endinclude",
		  "$end_expansion", "$end_translat", "$end_transform")
_PREPRO_BLOCK  = ("#if", "#ifdef", "#ifndef", "#elif", "#else", "#endif")

# Region types
REGION_NORMAL    = 0
REGION_BLACKHOLE = 1
REGION_LATTICE   = 2
REGION_VOXEL     = 3

ERROR   = "error"
_VOLUME = "@volume"

#-------------------------------------------------------------------------------
# flair vector returning a string with {} representation
class Vector(bmath.Vector):
	def __str__(self): return "{%s}"%(list.__str__(self)[1:-1])

#-------------------------------------------------------------------------------
_globalDict = {
	# Mathematical constants
	"fwhm"	: 2.0*math.sqrt(2.0*math.log(2.0)),

	# Physics constants
	"a"	: 1.0 / 137.035989561,
	"h"	: 4.135667516e-15,	# eV s
	"c"	: 299792458e2,
	"Na"	: 6.0221367e23,
	"qe"	: 1.60217646e-19,
	"re"	: 2.8179409183694872e-13,

	# Masses
	"amu"	: 0.93149432,
	"amuC12": 0.93123890,
	"amugr"	: 1.6605402E-24,
	"Mp"	: 0.93827231,
	"Mn"	: 0.93956563,
	"Me"	: 0.510998910e-3,

	# Distances
	"nm"	: 0.1E-6,
	"um"	: 0.1E-3,
	"mm"	: 0.1,
	"cm"	: 1.0,
	"dm"	: 10.0,
	"m"	: 100.0,
	"km"	: 100.0e3,

	# Imperial system distances
	"inch"	: 2.54,
	"feet"	: 30.48,
	"ft"	: 30.48,
	"mile"	: 160934.4,
	"mi"	: 160934.4,

	# Time
	"ns"	: 1.0e-9,
	"us"	: 1.0e-6,
	"ms"	: 0.001,
	"s"	: 1.0,
	"min"	: 60.0,
	"hour"	: 3600.0,
	"day"	: 86400.0,
	"week"	: 7.0*86400.0,
	"month"	: 365.25/12.0*86400.0,
	"year"	: 365.25*86400.0,

	# Prefixes
	"yotta" : 1e24,
	"zetta" : 1e21,
	"exa"   : 1e18,
	"peta"  : 1e15,
	"tera"  : 1e12,
	"giga"  : 1e9,
	"mega"  : 1e6,
	"kilo"  : 1e3,
	"hecto" : 1e2,
	"deca"  : 1e1,

	"deci"  : 1e-1,
	"centi" : 1e-2,
	"milli" : 1e-3,
	"micro" : 1e-6,
	"nano"  : 1e-9,
	"pico"  : 1e-12,
	"femto" : 1e-15,
	"atto"  : 1e-18,
	"zepto" : 1e-21,
	"yocto" : 1e-24,

	# Energy
	"eV"	: 1e-9,
	"keV"	: 1e-6,
	"MeV"	: 1e-3,
	"GeV"	: 1.0,
	"TeV"	: 1e3,
	"PeV"	: 1e6,
	"J"     : 1.0/1.60217646e-10,

	# Angle
	"deg"	: math.pi/180.0,
	"rad"	: 1.0,
	"mrad"	: 0.001,

	# Trig functions in degrees
	"sind"	: lambda x : math.sin(math.radians(x)),
	"cosd"	: lambda x : math.cos(math.radians(x)),
	"tand"	: lambda x : math.tan(math.radians(x)),
	"asind"	: lambda x : 180.*math.asin(x)/math.pi,
	"acosd"	: lambda x : 180.*math.acos(x)/math.pi,
	"atand"	: lambda x : 180.*math.atan(x)/math.pi,

	# Generic functions
	"abs"	: abs,
	"float"	: float,
	"int"	: int,
	"len"	: len,
	"str"	: str,

	# Physics functions
	"T2g"	: lambda T,m : 1. + T/m,
	"g2b"	: lambda g   : math.sqrt(1.-1./g**2),
	"b2g"	: lambda b   : 1./math.sqrt(1.-b**2),
	"p2T"	: lambda p,m : math.sqrt(p**2 + m**2) - m,
	"T2p"	: lambda T,m : math.sqrt(T**2 + 2.*m*T),
	"dT2dp"	: lambda T,dT,m : (T+m)/math.sqrt(T**2 + 2.*m*T)*dT,
	"omega" : lambda x : 2.*math.pi*(1.-math.cos(x)),
	"omegad": lambda x : 2.*math.pi*(1.-math.cos(math.radians(x))),

	"X0"    : lambda Z,A: 716.4*float(A)/(float(Z)*(Z+1.)*math.log(287./math.sqrt(float(Z)))),

	# Vector and matrices
	"Vector": Vector,
	"Matrix": bmath.Matrix,
}

# Populate math functions on the local dict
for name in dir(math):
	if name[0]=="_": continue
	_globalDict[name] = getattr(math, name)

#===============================================================================
# Local Dictionary
#===============================================================================
class LocalDict(dict):
	"""Implement the local dictionary"""
	def __init__(self, input):
		self.input = input
		self.card  = None
		self.clear()

	# ----------------------------------------------------------------------
	# Clear dictionary and set default variables
	# ----------------------------------------------------------------------
	def clear(self):
		dict.clear(self)
		self["w"] = lambda w,s=self:     s.card.numWhat(w)
		self["W"] = lambda w,s=self:     s.card.evalWhat(w)
		self["b"] = lambda n,w,s=self:   s.bodyWhat(n,w)
		self["C"] = lambda t,n,w,s=self: s.cardWhat(t,n,w)

	# ----------------------------------------------------------------------
	# Default getitem
	# ----------------------------------------------------------------------
	def __getitem__(self, item):
		try:
			return dict.__getitem__(self,item)
		except:
			# check defines
			defines = self.input.cardsCache("#define")
			# dummy linear search
			for define in defines:
				if define.ignore(): continue
				if define.sdum().strip()==item:
					# Convert to int or float if possible
					val = define.evalWhat(1)
					try:
						f = float(val)
						if float(int(f)) == f:
							return int(f)
						return f
					except:
						return val
			if item not in _globalDict: return str(item)
			# FIXME check for a card like BEAM(0,1)
			#if cards: return lambda n,w,t=item,s=self: s.cardWhat(t,n,w)
			#return str(item)
			raise

	# ----------------------------------------------------------------------
	def bodyWhat(self, name, what):
		bodies = self.input.cardsCache("bodies")
		body = bodies[name]
		return body.numWhat(what)

	# ----------------------------------------------------------------------
	def cardWhat(self, tag, name, what):
		cards = self.input.cardsCache(tag)
		if cards is None:
			raise Exception("Card %s do not exist"%(tag))
		if isinstance(name,str):
			for card in cards:
				if card.ignore(): continue
				if card.sdum()==name:
				#	return card.numWhat(what)
					val = card.what(what)
			raise Exception("No card %s with sdum=%s found\n"%(tag,name))
		else:
			#return cards[name].numWhat(what)
			val = cards[name].what(what)

		# return a float or int otherwise as is
		try:
			f = float(val)
			if float(int(f)) == f:
				return int(f)
			return f
		except:
			return val

#-------------------------------------------------------------------------------
# Write in utf format if necessary
#-------------------------------------------------------------------------------
def utfWrite(f, s):
#	try:
#		f.write(s)
#	except (UnicodeDecodeError, UnicodeEncodeError):
#		f.write(s.encode("utf-8"))
		if isinstance(s,str) :
			f.write(s)
		else:                        
			f.write(s.decode())

#-------------------------------------------------------------------------------
def utfWriteln(f, s):
#	try:
#		f.write(s+"\n")
#	except (UnicodeDecodeError, UnicodeEncodeError):
		if isinstance(s,str) :
			f.write(s+"\n")
		else:                        
			f.write(s.decode()+"\n")

#-------------------------------------------------------------------------------
def writeln(f, s):
	f.write(s+"\n")

#-------------------------------------------------------------------------------
# Readline and convert it to unicode if necessary
#-------------------------------------------------------------------------------
def utfReadline(f):
	line = f.readline()
	try:
		uline = line.decode("utf-8")
	except (AttributeError, UnicodeDecodeError, UnicodeEncodeError):
		return line
	try:
		sline = str(uline)
		return line
	except (AttributeError, UnicodeDecodeError, UnicodeEncodeError):
		return uline

#-------------------------------------------------------------------------------
def tagCmp(a, b):
	return cmp(a.tag, b.tag)

#-------------------------------------------------------------------------------
def isEvalStr(w):
	"""return is the what string is an evaluated function"""
	return isinstance(w,str) and w!="" and (w[0] in ("=","$") or w[:2]=="-$" or w[:2]=="+$")

#===============================================================================
# Structures as placeholders
#===============================================================================
class LowNeutMaterial:
	pass

#===============================================================================
# CardInfo class
#===============================================================================
class CardInfo:
	_db     = {}		# Card information dictionary
	_funcs  = {}		# Caching compiled functions

	def __init__(self, name, group, range_, extra, assert_, default, meaning):
		if name[0] not in ("$","#"):
			self.tag = name[0:8]
		else:
			self.tag = name
		self.name    = name
		self.group   = group
		self.range   = range_
		self.assertTrue = assert_
		self.extra   = extra
		self.default = default
		self.meaning = meaning
		self.layout  = None
		self.order   = 99999
		self.disableComment = False	# Comment disabled cards instead
						# of #if 0..#endif
		self.underline = -1

		# process range for units, particles, materials and regions
		self.rangeCode  = []
		self.assertCode = []
		self.whats      = []
		self.useUnits   = False		# True if any case is using units
		if range_ is None: return
		nwhats = 0
		for r in range_:
			if "lu" in r or "slu" in r: self.useUnits = True
			self.rangeCode.append(self.compileConditions(r))
			w = len(r)
			for j in range(len(r)-1,-1,-1):
				if r[j] != "-": break
				w -= 1
			nwhats = max(nwhats, w)
		self.nwhats = nwhats

		# compile assertion code
		for assertList in self.assertTrue:
			ass = []
			for a in assertList:
				ass.append(self.compileAssert(a))
			self.assertCode.append(ass)

	# ----------------------------------------------------------------------
	def setDisableComment(self, c=False):
		self.disableComment = c

	# ----------------------------------------------------------------------
	# Compile conditions
	# The code which is generated is in the form i.e.:
	#	lambda x: x>0
	# and should be executed with:
	#	eval(code)(x)
	# ----------------------------------------------------------------------
	def compileConditions(self, rng):
		# the switch condition is given by the f() fields and
		# by the sdum if is defined

		# Check sdums
		code = []

		if rng[0]!="-" and rng[0][0].isupper():
			sdums = rng[0].split(":")
			maxlen = max([len(x) for x in sdums])
			if maxlen>8:
				raise Exception("Error %d %s"%(maxlen,rng[0]))
			if len(sdums) == 1:
				s = "lambda x:x==\"%s\"" % (rng[0])
				code.append((0, 'f', CardInfo._compileFunc(s)))
			elif len(sdums) > 0:
				s = "lambda x:x[:8] in ("
				for i in sdums:
					if i == "-": i=""
					s += "\"%s\"," % (i)
				s += ")"
				code.append((0, 'f', CardInfo._compileFunc(s)))

		# Check if all [fri]() conditions are met
		i = 0
		for w in rng:
			if len(w)>4 and w[1] == "(":
				s = "lambda x:%s" \
					% (w[2:-1].replace(',', ' or '))
				code.append((i, w[0], CardInfo._compileFunc(s)))
			elif w[0] == 'c':
				s = "lambda x:len(str(x))<=%s"%(w[1:])
				code.append((i, w, CardInfo._compileFunc(s)))
			else:
				code.append((i, w, None))
			i += 1

		return code

	#-----------------------------------------------------------------------
	# Compile assert statement
	# First detect all what references w(#) and replace them with:
	#	w(#) -> c.numWhat(#)
	#	W(#) -> c.what(#)
	# identify all what-indexes to report as problems
	# and request in the code that are not zero!
	#
	# finally compile the code
	# example: w(8)>w(7)
	#	lambda c: c.numWhat(7)==0 or c.numWhat(8)==0 or c.numWhat(8)>c.numWhat(7)
	#-----------------------------------------------------------------------
	def compileAssert(self, assertStatement):
		whatUsed = [False]*21
		for match in _WNPAT.finditer(assertStatement):
			try: whatUsed[int(match.group(2))] = True
			except: pass
		ass = _wPAT.sub('c.numWhat(', assertStatement)
		ass = _WPAT.sub('c.what(', ass)

		whats = []
		expr = "lambda c:"
		for i,w in enumerate(whatUsed):
			if w:
				if i<20:
					whats.append(i)
					expr += "c.what(%d)=='' or "%(i)
				else:
					whats.append(-1)
		expr += "("+ass+")"
		return (whats, CardInfo._compileFunc(expr))

	#-----------------------------------------------------------------------
	# Compile or use pre-defined function
	#-----------------------------------------------------------------------
	@staticmethod
	def _compileFunc(s):
		try:
			return CardInfo._funcs[s]
		except:
			c = compile(s, "<string>", "eval")
			CardInfo._funcs[s] = c
			return c

	# ----------------------------------------------------------------------
	# Return name of compiled function (used for debugging)
	# ----------------------------------------------------------------------
	@staticmethod
	def _compileFuncName(code):
		for n,c in list(CardInfo._funcs.items()):
			if code == c:
				return n

	# ----------------------------------------------------------------------
	# Find card's range(case)
	# ----------------------------------------------------------------------
	def findCase(self, card):
		if len(self.rangeCode)==1: return 0
		for r, code in enumerate(self.rangeCode):
			# Check all 'f'-flags
			for w,f,c in code:
				if f=='f' and c:
					v = card.evalWhat(w)
					if w>0 and v=="": v=0
					try: v= int(float(v))
					except: pass
					try:
						if not eval(c)(v): break
					except:
						pass
			else:
				# All conditions met
				return r
		return 0 # return default case

	# ----------------------------------------------------------------------
	# Validate card: return a list of failing variables.
	# ----------------------------------------------------------------------
	def validate(self, card, case=None):
		if self is CardInfo._db[None]:
			return ["Error card"]

		if card.input is None: return [None]

		if case is None:
			case = self.findCase(card)
			if case is None: return [None]	# No case found!

		# FIXME check special cases
		# 1. Regions check body names inside
		# 2. All references of names if they exist mi, ri, bi, di, vi, ti
		# 3. Check step for ranges to be >= 0
		# 4. Check limits in USRxxx cards
		# 5. Check logical units usage
		# 6. Check duplicate name definition of bodies, regions...

		# Check all fields
		pb = []
		for w,f,c in self.rangeCode[case]:
			x = card.evalWhat(w)
			if x=="": continue	# Use default

			if f=='i':
				try:
					xi = int(float(x))
					if abs(float(x)-float(xi))>1e-9:
						pb.append(w)
						continue
					x = xi
				except:
					#x = 0
					pb.append(w)
					continue

			elif f=='l':
				try:
					xi = int(float(x))
					if abs(float(x)-float(xi))>1e-9:
						pb.append(w)
						continue
					x = xi
					# Special patch for RADDECAY to restore
					# the long type instead of float
					card.setWhat(w, xi)
				except:
					#x = 0
					pb.append(w)
					continue

			elif f=='r':
				try: x = float(x)
				except:
					# Check if it is a fortran number with D or internal spaces
					#try: x = float(re.sub("[dD]","e",x.replace(" ","")))
					try: x = float(re.sub("[dD]","e",x))
					except:
						pb.append(w)
						continue

			elif f=='f':
				if x=="" and w>0: x = 0
				try: x = int(float(x))
				except: pass

			elif f=="pi" or f=="spi":
				try:
					x = float(x)
					if -100.0<x<1000.0 and abs(int(x)) not in Particle._db:
						pb.append(w)
				except:
					if len(x) > 0:
						if x[0]=="-": x=x[1:]
						if x not in Particle._db:
							pb.append(w)

			# FIXME ri, mi, bi, di, vi
			try:
				if c and not eval(c)(x):
					pb.append(w)
			except:
				pb.append(w)

		# Check all assert statements
		if self.assertTrue:
			for i,(whats,code) in enumerate(self.assertCode[case]):
				if whats and not eval(code,_globalDict)(card):
					pb.extend(whats)
					pb.append("Assertion failed: %s"%(self.assertTrue[case][i]))

		# Special check for regions
		if card.tag == "REGION":
			sdum = card.sdum()
			if sdum!="&" and not _NAMEPAT.match(sdum):
				pb.append(0)
			exp = csg.tokenize(card.extra())
			# Check expression
			try:
				csg.exp2rpn(exp)
				if exp and not csg.check(exp):
					pb.append(-1)
					pb.append("Invalid expression")
				else:
					bodies = card.input.cardsCache("bodies")
					notadded = True
					for token in exp:
						if token in ("-", "+", "|", "@"):
							continue
						elif _NAMEPAT.match(token):
							# XXX Check in body list
							if token not in bodies:
								if notadded:
									pb.append(-1)
									notadded = False
								pb.append("Invalid body %r"%(token))
						else:
							if notadded:
								pb.append(-1)
								notadded = False
							pb.append("Invalid token %r"%(token))
			except csg.CSGException:
				pb.append(-1)
				pb.append("Unbalanced parenthesis")

		elif card.isGeo() and len(card.tag) == 3 and card.tag!="END":
			if not _NAMEPAT.match(card.sdum()):
				pb.append(0)

		elif card.tag == "COMPOUND":
			# Check for nesting compounds
			# XXX for them moment on level checking
			for w in range(2,7,2):
				if card.evalWhat(w) == card.sdum():
					pb.append(w)

		return pb

	# ----------------------------------------------------------------------
	# convert to names
	# ----------------------------------------------------------------------
	def toNames(self, card, case=None):
		if card.input is None: return

		if case is None:
			case = self.findCase(card)

		if case is None: return
		range_ = self.range[case]

		# Special cases
		if card.tag == "LATTICE":
			# Check the sdum and convert it to something that flair
			# and FREE format of FLUKA likes
			sdum = card.sdum().upper()
			if sdum[:4] == "ROT#" or sdum[:3] == "RO#":
				try:
					id = int(sdum[sdum.find('#')+1:])
				except:
					id = 999
				card.setSdum("rot#%03d"%(id))
		elif card.tag == "AUXSCORE":
			try:
				w = int(card.evalWhat(1))
				card.setWhat(1, ("USRBDX", "USRBIN", "USRTRACK", "USRCOLL", "USRYIELD")[w-1])
			except:
				pass

		# Go through the various whats
		for i in range(min(card.nwhats(), len(range_))):
			r = range_[i]
			w = card.what(i)
			if w == "": continue
			try: aw = int(float(w))
			except: continue
			if   r in ("mi", "smi"):
				aw = abs(aw)
				if aw==0:
					card.setWhat(i, "")
					continue
				# FIXME to correct treatment of defines...
				lst = card.input.materialList()
				aw -= 1
				if aw >= len(lst):
					card.setAbsWhat(i, "@LASTMAT")
				else:
					card.setAbsWhat(i, lst[aw])

			elif r in ("ri", "sri"):
				if aw==0:
					card.setWhat(i, "")
					continue
				if r == "ri":
					if aw==-1:
						card.setWhat(i, "@ALLREGS")
						continue
					elif aw<0:
						continue

				aw = abs(aw)

				# FIXME to correct treatment of defines...
				try:
					regions = card.input.cardsCache("REGION")
					aw -= 1

					# Normally should not enter here unless if the voxel is not found
					if aw >= len(regions):
						if "VOXELS" in card.input.cards:
							# FIXME no way to know about @LASTREG!
							if aw == len(regions):
								card.setAbsWhat(i, "VOXEL")
							else:
								if (aw-len(regions)) < 1000:
									card.setAbsWhat(i,"VOXEL%03d"%(aw-len(regions)))
								elif (aw-len(regions)) < 10000:
									card.setAbsWhat(i,"VOXE%04d"%(aw-len(regions)))
								else:
									card.setAbsWhat(i,"VOX%05d"%(aw-len(regions)))
						else:
							card.setAbsWhat(i, "@LASTREG")
					else:
						card.setAbsWhat(i, regions[aw].name())
				except KeyError:
					pass

			elif r in ("pi", "spi"):
				if aw==0:
					card.setWhat(i, "RAY")
					continue
				if r=="spi":
					# Special treatment for DISCARD and negative particle id's
					if aw > 980:
						if aw<1000: card.setSign(i,True)
						aw = -abs(1000-aw)
						card.setAbsWhat(i, Particle.convert(aw, True))
					else:
						card.setAbsWhat(i,
							Particle.convert(abs(aw), True))
				else:
					card.setWhat(i, Particle.convert(aw, True))

			elif r=="di":
				aw = abs(aw)
				if aw==0:
					card.setWhat(i, "")
					continue
				# Type of detector in sdum
				detectors = card.input.cardsCache(card.sdum(), 0)
				if detectors:
					aw -= 1
					if aw>=0 and aw<len(detectors):
						name = detectors[aw]
						# Ensure that there is only one reference
						if detectors.index(name) == aw:
							card.setWhat(i, name)

			elif r=="bi":
				aw = abs(aw)
				if aw==0:
					card.setWhat(i, "")
					continue
				# FIXME detectors EVENTBIN & USRBIN
				detectors = card.input.cardsCache("USRBIN", 0)
				if detectors:
					aw -= 1
					if aw>=0 and aw<len(detectors):
						name = detectors[aw]
						# Ensure that there is only one reference
						if detectors.index(name) == aw:
							card.setWhat(i, name)

	# ----------------------------------------------------------------------
	# find
	# @return list  - of references
	#	  tuple - for a range
	# ----------------------------------------------------------------------
	def find(self, target, case):
		starget = "s"+target

		# Go through the various whats
		found = []
		if case is None: return found
		range_ = self.range[case]
		for i in range(len(range_)):
			r = range_[i]
			if r==target or r==starget:
				found.append(i)
			if r=="i" and len(found)==2:
				found.append(i)
				return tuple(found)
		return found

	# ----------------------------------------------------------------------
	# get card information
	# ----------------------------------------------------------------------
	@staticmethod
	def get(card):
		if isinstance(card, Card):
			return CardInfo._db.get(card.tag, CardInfo._db[None])
		else:
			if card and card[0] not in ("$","#"):
				card = card[0:8]
			return CardInfo._db.get(card, CardInfo._db[None])

	# ----------------------------------------------------------------------
	@staticmethod
	def none():
		return CardInfo._db[None]

#===============================================================================
# Logical Units
#===============================================================================
class LogicalUnits:
	def __init__(self):
		self.reset()

	# ----------------------------------------------------------------------
	# Reset to default units
	# ----------------------------------------------------------------------
	def reset(self):
		self.list = [None]*100

		# Add some default units:
		self.list[ 1] = (False, "RANDOMIZ",	None)
		self.list[ 2] = (False, "system",	None)
		self.list[ 3] = (False, "system",	None)
		self.list[ 4] = (False, "system",	None)
		self.list[ 5] = (False, "GEOBEGIN",	None)	# stdin
		self.list[ 6] = (False, "stderr",	None)
		self.list[ 7] = (False, "system",	None)
		self.list[ 8] = (False, "system",	None)
		self.list[ 9] = (False, "LOW-NEUT",	None)
		self.list[10] = (False, "system",	None)
		self.list[11] = (False, "GEOBEGIN",	None)
		self.list[12] = (False, "system",	None)
		self.list[13] = (False, "EMFFLUO",	None)
		self.list[14] = (False, "data",		None)
		self.list[15] = (False, "GEOBEGIN",	None)
		self.list[16] = (False, "system",	None)
		self.list[17] = (False, "DETECT",	None)
		self.list[18] = (False, "system",	None)
		self.list[19] = (False, "DPMJET",	None)
		self.list[20] = (False, "system",	None)
		#self.list[49] = (True,  "USERDUMP",	None)

	# ----------------------------------------------------------------------
	# Assign a new unit type
	# ----------------------------------------------------------------------
	def assign(self, unit, tag, fn=None):
		"""Assign a new unit to tag"""
		try:
			unit   = int(unit)
			binary = unit<0
			unit   = abs(unit)
			self.list[unit] = (binary, tag, fn)
		except:
			pass

	# ----------------------------------------------------------------------
	def clear(self, unit):
		"""Clear a unit from list"""
		unit = abs(unit)
		self.list[unit] = None

	# ----------------------------------------------------------------------
	# Return unit type
	# ----------------------------------------------------------------------
	def get(self, u):
		"""Return unit information"""
		try:
			u = abs(int(u))
			return self.list[u]
		except:
			return None

	# ----------------------------------------------------------------------
	# Find first free
	# ----------------------------------------------------------------------
	def free(self):
		"""Return first free unit"""
		for i in range(21, len(self.list)):
			u = self.list[i]
			if u is None: return i

	#-----------------------------------------------------------------------
	# FIXME needs optimisation, we should find all cards using units in the
	# CardInfo and search only those cards
	#-----------------------------------------------------------------------
	def scan(self, input):
		self.reset()
		for card in input.cardlist:
			# Ignore cards that use only as readonly SPECSOUR, OPEN is special
			if card.tag in ("OPEN","SPECSOUR"): continue
			for u in card.units():
				if u == 0: continue
				unit = self.get(u)
				# check for a possible conflict!
				if unit is not None:
					(b, c, f) = unit
					if c!=card.tag or (b ^ (int(u)<0)):
						# Find what with units
						case = card.case()
						card.invalid = card.info.find("lu", case)
						say("WARNING: Possible conflict with unit: %s"%(str(u)))
						say("is used by: %s  #%d"%(str(self.get(u)[1]),card.pos()+1))
						say(str(card))
						break
				self.assign(u, card.tag)

		# For all open cards assign the correct name
		for card in input["OPEN"]:
			unit = card.intWhat(1)
			u    = self.list[abs(unit)]
			if u is None:
				self.assign(unit, "OPEN", card.extra())
			else:
				self.assign(unit, u[1], card.extra())

	# ----------------------------------------------------------------------
	# Return a list of all units used
	# ----------------------------------------------------------------------
	def filterList(self, tag=None, bin=True, onlymatches=False):
		"""Return a list of all units, filtered for a card"""
		lst = [("", "")]
		for i in range(1, 100):
			u = self.list[i]

			if u is not None and tag is not None:
				if u[1]!=tag: continue
			elif u is None and onlymatches:
				continue

			if bin:
				lst.append(("%02d ASC" % (i), i))
				lst.append(("%02d BIN" % (i), -i))
			else:
				lst.append(("%02d" % (i), i))
		lst.append(("0", 0))
		return lst

#===============================================================================
# Particle class
#===============================================================================
class Particle:
	# Particles dictionary
	_db           = {}
	list          = []
	hadron        = []
	lepton        = []
	beam          = []
	listAll       = []
	signedList    = []
	signedListAll = []
	particles     = []

	# ----------------------------------------------------------------------
	def __init__(self, id, name, mass=0.0, comment=None, pdg=0):
		self.id      = id
		self.name    = name
		self.mass    = mass
		self.comment = comment
		self.pdg     = pdg

	# ----------------------------------------------------------------------
	# add particle
	# ----------------------------------------------------------------------
	@staticmethod
	def add(id, name, mass=0.0, comment=None, pdg=0):
		p = Particle(id, name, mass, comment, pdg)
		Particle._db[id]   = p
		Particle._db[name] = p

		# Populate also _baseLocalDict
		if mass>0.0:
			n = name
			if name[-1] == "-":
				n = n[:-1]+"n"
			elif name[-1] == "+":
				n = n[:-1]+"p"
			if n.find("-")>=0:
				n = n.replace("-","")
			_globalDict["m%s"%(n)] = mass
		return p

	# ----------------------------------------------------------------------
	def __getitem__(self, idx):
		return Particle.particles[idx]

	# ----------------------------------------------------------------------
	@staticmethod
	def get(idx):
		return Particle._db.get(idx)

	# ----------------------------------------------------------------------
	# create lists
	# ----------------------------------------------------------------------
	@staticmethod
	def makeLists():
		if len(Particle.list)>0: return
		lastpar = Particle.convert("@LASTPAR",False)
		# All particles
		Particle.list.append("")
		Particle.hadron.append("")
		Particle.lepton.append("")
		Particle.listAll.append("")
#		for p in Particle._db.keys():
#			if isinstance(p, str):
		# Go through in ordered by id list
		for i in sorted([Particle._db[x].id for x in list(Particle._db.keys()) if isinstance(x,str)]):
			part = Particle._db[i]
			if part.name != "@LASTPAR":
				Particle.listAll.append(part.name)
			if part.id < lastpar:
				Particle.list.append(part.name)
			if abs(part.pdg)<100:
				Particle.lepton.append(part.name)
			else:
				Particle.hadron.append(part.name)

		# Beam particles
		Particle.beam = Particle.list[:]
		Particle.beam.append("ISOTOPE")
		Particle.beam.append("AMBE")
		Particle.beam.append("AMB")
		Particle.beam.append("252CF")
		Particle.beam.append("D-D")
		Particle.beam.append("D-T")
		Particle.hadron.append("ISOTOPE")
		Particle.beam.sort()
		Particle.hadron.sort()
		Particle.lepton.sort()

		# Only real particles + generalized ones
		Particle.particles = Particle.list[:]
		for p in range(201, 228):
			if p in (208, 211, 219, 220, 221, 222): continue
			Particle.particles.append(Particle._db[p].name)

		# Sort
		Particle.list.append("@LASTPAR")
		Particle.listAll.append("@LASTPAR")

		# Create signed list
		Particle.signedList += Particle.list
		for i in range(len(Particle.signedList)-1, 0, -1):
			p = Particle.signedList[i]
			if len(p)>0 and p[0] != '@':
				Particle.signedList.insert(i+1, '-'+p)

		Particle.signedListAll += Particle.listAll
		for i in range(len(Particle.signedListAll)-1, 0, -1):
			p = Particle.signedListAll[i]
			if len(p)>0 and p[0] != '@':
				Particle.signedListAll.insert(i+1, '-'+p)

	# ----------------------------------------------------------------------
	# Convert to integer or string
	# ----------------------------------------------------------------------
	@staticmethod
	def convert(pid, tonames=True):
		"""Convert particle from name to number
		pid particle type
		tonames True/False"""

		try:
			idx = int(pid)
			isInt = True
		except:
			if pid == "":
				idx = 0
				isInt = True
			else:
				isInt = False

		if tonames:
			if isInt:
				if idx<=-100 or idx>1000: return str(pid)
#				if 988 <= idx <= 1012:
#					# Only for DISCARD card
#					idx = -abs(1000 - idx)
				try:
					return Particle._db[idx].name
				except:
					say("Particle %s not found!"%str(pid))
					return "UNKNOWN"
			else:
				return pid
		else:
			if isInt:
				return idx
			try:
				return Particle._db[pid].id
			except:
				raise Exception("Particle %s not found!"%(pid))

	# ----------------------------------------------------------------------
	# scan input for particles
	# ----------------------------------------------------------------------
	@staticmethod
	def scan(input):
		# FIXME for the moment return everything
		return Particle.listAll

#===============================================================================
# Voxel class
#	parse the voxel file and create fake cards for the
#		materials, assignmats and corrfact
#	that might be embedded in the voxel file
#===============================================================================
class Voxel:
	def __init__(self, filename=None):
		self.cleanup()
		if filename is not None: self.read(filename)

	# ----------------------------------------------------------------------
	# Clean up voxel structure
	# ----------------------------------------------------------------------
	def cleanup(self):
		self.filename = ""
		self.title  = ""
		self.nx     = 0
		self.ny     = 0
		self.nz     = 0
		self.no     = 0
		self.mo     = 0
		self.dx     = 0.0
		self.dy     = 0.0
		self.dz     = 0.0
		self.kreg   = None
		self.roiN    = 0
		self.roiName  = {}	# roi structure names
		self.roiColor = {}	# colors of rois
		#self.roiComb = {}	# roi dictionary of combination of structures
		#self.roiData = None	# roi structure data
		self.input  = Input()	# Input cards embedded in voxel file
		self.input.format = FORMAT_VOXEL	# Mark as voxel

	# ----------------------------------------------------------------------
	# Read voxel information
	# ----------------------------------------------------------------------
	def read(self, filename):
		self.filename = filename
		try:
			f = open(filename,"rb")
		except:
			return False

		# Title CHAR*80
		record = fortran.read(f)
		if record is None or len(record) != 80:
			say("Invalid voxel file. No title found")
			f.close()
			return False
		self.title = record.strip()

		# Memory Size nx,ny,nz,no,mo
		record = fortran.read(f)
		if record is None or len(record) not in (20,32,36):
			say("Invalid voxel file. Wrong memory size")
			f.close()
			return False
		if len(record)==20:
			self.nx, self.ny, self.nz, self.no, self.mo = \
				struct.unpack("=5i", record)
			self.roiN = 0
			roiPlanar = 0
		else:
			self.nx, self.ny, self.nz, self.no, self.mo, self.roiN, roiCombN, roiMaxLen, roiPlanar = \
				struct.unpack("=9i", record)

		# Create voxel regions
		card = Card("REGION",["VOXEL"])
		#card["@id"] = 0
		card["@type"] = REGION_VOXEL
		card.setExtra("+VOXEL")
		self.input.addCard(card)

		for i in range(self.no+1):	# +1 to include organ 0
			card = Card("REGION",[Voxel.regionName(i+1)])
			#card["@id"] = i+1
			self.input.addCard(card)
			card["@type"] = REGION_VOXEL

		# Voxel dimension
		record = fortran.read(f)
		if record is None or len(record) != 3*8:
			say("Invalid voxel file. Wrong dimensions")
			f.close()
			return False
		self.dx, self.dy, self.dz = struct.unpack("=3d", record)

		# Voxel data (skip)
		skip = fortran.skip(f)
		if skip != self.nx*self.ny*self.nz*2:
			say("Invalid voxel file. Data size do not match")
			f.close()
			return False

		# kreg
		record = fortran.read(f)
		if record is None or len(record) != 2*self.mo:
			say("Invalid voxel file. Wrong kreg array")
			f.close()
			return False
		self.kreg = struct.unpack("=%dH"%(len(record)/2), record)

		# roiStruct
		if self.roiN>0:
			# structures
			#try:
				for i in range(1,self.roiN+1):
					record = fortran.read(f)
					color, name = struct.unpack("=i64s",record)
					name = name.strip()
					if name:
						self.roiColor[i] = color
						self.roiName[i]   = name

				# for the moment skip everything
				for i in range(roiCombN):
					record = fortran.read(f)	# skip index
				record = fortran.read(f)		# skip data

				# Read planar structures
				for i in range(roiPlanar):
					record = fortran.read(f)	# skip index
					roi,n = struct.unpack("=ii", record)
					for j in range(n):
						record = fortran.read(f)# skip data
			#except:
			#	pass

		# Read embedded cards
#		print (" input voxel  start reading cards ")
		while True:
			try: record = fortran.read(f)
			except IOError: break
			if record is None: break
#			print (" record " ,record)
			if len(record) != 80:
				say("Invalid voxel file. Wrong card found '%s'"%(record))
				f.close()
				return False
			record = record.decode()
			tag  = record[:8].strip()
			sdum = record[-10:].strip()
			what = [sdum,			# sdum
				record[10:20].strip(),	# what(1)
				record[20:30].strip(),	# what(2)
				record[30:40].strip(),	# what(3)
				record[40:50].strip(),	# what(4)
				record[50:60].strip(),	# what(5)
				record[60:70].strip()]	# what(6)
			if tag=="COMPOUND" and card.tag==tag and card.sdum()==sdum:
				# append whats to previous card
				card.appendWhats(what[1:])
			else:
				# if COMPOUND check for existing to append
				card = Card(tag, what)
				self.input.addCard(card)

		f.close()
		return True

	# ----------------------------------------------------------------------
	def __str__(self):
		return "Voxel\n" \
			"File:\t%s\n" \
			"Title:\t%s\n" \
			"nxyz:\t%d x %d x %d\n" \
			"no,mo:\t%d\t%d\n" \
			"dxyz:\t%g x %g %g\n" \
			%(self.filename, self.title,
			 self.nx, self.ny, self.nz,
			 self.no, self.mo,
			 self.dx, self.dy, self.dz)

	# ----------------------------------------------------------------------
	# return a valid region name
	# ----------------------------------------------------------------------
	@staticmethod
	def regionName(idx):
		if idx==0: return "VOXEL"
		num = "%03d"%(idx)
		return "%s%s"%("VOXEL"[:(8-len(num))], num)

#===============================================================================
# Card
#===============================================================================
class Card:
	"""Base class for a fluka card"""

	GENERIC = 0		# card types
	BODY    = 1
	REGION  = 2
	OBJECT  = 3

	# ----------------------------------------------------------------------
	def __init__(self, tag, what=None, comment="", extra=""):
		"""Initialise a fluka card
		Each card contains a tag, whats, one multi-line comment and sometimes
		some extra information
		sdum=what[0]
		continuation cards have multiple whats of 1+6*lines.
		extra is used for TITLE, GEOBEGIN and REGION
		"""

		if what is None:
			what = [""]

		self.enable  = True	# Enable or Disabled
		self.active  = True	# Preprocessor activate/deactivate
		#self.expand  = True	# Show or hide card and subcontents
		self.invalid = None	# Valid or invalid card
		self.prop    = None	# User properties
		self.input   = None	# input class holding card
		self.tag     = tag	# card tag
		self._what   = what	# what list
		self._owhat  = what	# original what list as read from input
		self._sign   = []
		self._extra  = extra
		self._pos    = -1	# if -1 the card is most probably deleted
		self._indent = 0	# preprocessor indent level
		self._geo    = False
		self._type   = Card.GENERIC
		self._userInvalid = None	# User errors for the card
		self.setComment(comment)
		self.changeTag(tag, False)

	# ----------------------------------------------------------------------
	def setModified(self):
		"""set last time modified"""
		self._modified = time.time()

	# ----------------------------------------------------------------------
	# Change tag of card, reduce the number of whats to the new card
	# XXX Do not forget to update the Input.cards... before
	# ----------------------------------------------------------------------
	def changeTag(self, tag, truncate=True):
		"""change tag of card"""
		self.tag = tag
		# find card information from _cardInfo dictionary
		self.info = CardInfo.get(tag)
		if self.info.name == ERROR: return

		if len(tag)==3 and "Geometry" in self.info.group:
			self._geo = True
		else:
			self._geo = tag in ("REGION", "LATTICE", "LATTSNGL", "VOXELS") or \
				    tag[0]=="$"	# transform types

		if tag in BODY_TAGS:
			self._type = Card.BODY
		elif tag in FLAIR_TAGS:
			self._type = Card.OBJECT
		elif tag == "REGION":
			self._type = Card.REGION
		else:
			self._type = Card.GENERIC

		# Cut extra whats
		if truncate and len(self._what) > self.info.nwhats:
			self.setNWhats(self.info.nwhats)
		self.setModified()

	# ----------------------------------------------------------------------
	# Clone card
	# ----------------------------------------------------------------------
	def clone(self):
		"""clone card object"""
		card = Card(	self.tag,
				self._what[:],
				self._comment,
				self._extra)
		card.enable = self.enable
		card._pos   = self._pos
		if self.prop:
			card.prop = self.prop.copy()	# Swallow copy
		return card

	# ----------------------------------------------------------------------
	# Copy properties from card
	# ----------------------------------------------------------------------
	def copyProperties(self, card):
		if card.prop:
			self.prop = card.prop.copy()	# Swallow copy

	# ----------------------------------------------------------------------
	# Append whats
	# ----------------------------------------------------------------------
	def appendWhats(self, what, pos=None):
		"""append what's"""
		if pos is None:
			self._what.extend(what)
		else:
			for w in what:
				if self.what(pos)=="":
					self.setWhat(pos,w)
				pos += 1
		self.setModified()

	# ----------------------------------------------------------------------
	# Pad whats to make them multiple of n*6+1
	# ----------------------------------------------------------------------
	def padWhats(self):
		rem = 6-(len(self._what)-1)%6
		if rem>0: self._what.extend([""]*rem)

	# ----------------------------------------------------------------------
	# Append comment
	# ----------------------------------------------------------------------
	def appendComment(self, comment):
		"""append comment"""
		if len(comment)==0: return
		if self._comment=="":
			self.setComment(comment)
		else:
			self.setComment(self._comment + "\n" + comment)
		self.setModified()

	# ----------------------------------------------------------------------
	# test card against a filter-function
	# Filter format  type [function] \n what \n value
	#
	# func: [type][function]
	#	type:	i=int
	#		f=float
	#		s=string
	#	function:
	#		a=abs
	# var:
	#	0..9 = whats
	#         *  = anywhere
	#	  b  = body
	#	  d  = detector
	#	  s  = sdum
	#         m  = material
	#         p  = particle	(range if tuple)
	#         r  = region		(range if tuple)
	#	  t  = transformation
	#         u  = unit
	# val:
	#	anything int, float or string
	#
	# case:
	#        True for matching case
	# ----------------------------------------------------------------------
	def test(self, func, var, val, case=False):
		"""check card against a filter function func"""
		if not case:
			try: val = val.upper()
			except: pass

		# Search everywhere
		if var=="*":
			if case:
				if self.tag.find(val)>=0: return True
				if self._comment.find(val)>=0: return True
				if self._extra.find(val)>=0: return True
				for w in self._what:
					if str(w).find(val)>=0: return True
			else:
				if self.tag.upper().find(val)>=0: return True
				if self._comment.upper().find(val)>=0: return True
				if self._extra.upper().find(val)>=0: return True
				for w in self._what:
					if str(w).upper().find(val)>=0: return True
			return False

		# Units
		elif var=="u":
			uval = abs(int(val))
			try:
				self.units().index(uval)
				return True
			except:
				return False

		# Particle
		elif var=="p":
			varlist = list(map(str,self.findType("pi")))
			if isinstance(varlist, list):
				# FIXME check for signed...
				if not case: varlist = list(map(str.upper, varlist))
				return val in varlist
			else:
				part         = Particle.convert(val, False)
				part_from    = Particle.convert(varlist[0], False)
				if part_from == part: return True
				part_to      = Particle.convert(varlist[1], False)
				if part_to == part: return True
				part_step    = varlist[2]
				if part_from < part < part_to:
					return ((part - part_from) % part_step) == 0
				return False

		# Region
		elif var=="r":
			if self.tag=="REGION":
				if case:
					return str(self.what(0)) == val
				else:
					return str(self.what(0)).upper() == val
			varlist = list(map(str, self.findType("ri")))
			if isinstance(varlist, list):
				# FIXME check for signed...
				if not case: varlist = list(map(str.upper, varlist))
				return val in varlist
			else:
				# Check for a range
				reg         = self.input.region(val, False)
				reg_from    = self.input.region(varlist[0], False)
				if reg_from == reg: return True
				reg_to      = self.input.region(varlist[1], False)
				if reg_to   == reg: return True
				reg_step    = varlist[2]
				if reg_from < reg < reg_to:
					return ((reg - reg_from) % reg_step) == 0
				return False

		# Material
		elif var=="m":
			if self.tag=="MATERIAL":
				if case:
					return str(self.what(0)) == val
				else:
					return str(self.what(0)).upper() == val
			varlist = self.findType("mi")
			if isinstance(varlist, list):
				# FIXME check for signed...
				if not case: varlist = list(map(str.upper, varlist))
				return val in varlist
			else:
				mat         = self.input.material(val, False)
				mat_from    = self.input.material(varlist[0], False)
				if mat_from == mat: return True
				mat_to      = self.input.material(varlist[1], False)
				if mat_to   == mat: return True
				mat_step    = varlist[2]
				if mat_from < mat < mat_to:
					return ((mat - mat_from) % mat_step) == 0
				return False

		# Transformation
		elif var=="t":
			# FIXME check also what(1) and Rot#\d+
			if self.tag=="ROT-DEFI":
				if case:
					return str(self.what(0)) == val
				else:
					return str(self.what(0)).upper() == val
			varlist = list(map(str,self.findType("ti")))
			if not case: varlist = list(map(str.upper, varlist))
			if val in varlist: return True
			# Check also -val
			if isinstance(val,int) or isinstance(val,float):
				if -val in varlist: return True
			elif isinstance(val,str):
				if "-%s"%(val) in varlist: return True
			return False

		# Detector
		elif var=="d":
			if self.tag in DETECTOR_TAGS:
				if case:
					return str(self.what(0)) == val
				else:
					return str(self.what(0)).upper() == val

			# FIXME ranges are not supported yet!
			varlist = list(map(str,self.findType("di")))
			varlist.extend(list(map(str,self.findType("bi"))))
			if not case: varlist = list(map(str.upper, varlist))
			return val in varlist

		# Body
		elif var=="b":
			if not self._geo: return
			if self.tag == "REGION":
				opt = re.MULTILINE|re.DOTALL
				if not case: opt |= re.I
				return re.match(r".*\b%s\b.*" % (val), self._extra,opt) is not None
			if len(self.tag)!=3: return False
			if case:
				return str(self.what(0))==val
			else:
				return str(self.what(0)).upper()==val

		#if var is None: return False
		if func=="ia":		# Absolute integer
			return abs(self.intWhat(var)) == int(val)

		elif func=="f":		# Floating point
			return self.numWhat(var) == float(val)

		elif func=="l":		# Long integer
			return self.longWhat(var) == int(val)

		elif func=="i":		# Integer
			return self.intWhat(var) == int(val)

		elif func=='sa':	# Absolute string (remove the minus '-')
			if var=="s":
				w = self.sdum()
			else:
				w = str(self.what(var))
			if len(w)>0 and w[0]=='-': w=w[1:]
			if not case: w = w.upper()
			return w == str(val)

		elif func=="s":
			if var=="s":
				return self.sdum() == str(val)
			else:
				return str(self.what(var)) == str(val)

		return False

	# ----------------------------------------------------------------------
	# Return card's cardInfo case (range)
	# ----------------------------------------------------------------------
	def case(self):
		return self.info.findCase(self)

	# ----------------------------------------------------------------------
	# Validate card values
	# ----------------------------------------------------------------------
	def validate(self, case=None):
		self.invalid = self.info.validate(self, case)
		if self._userInvalid:
			for item in self._userInvalid:
				if item not in self.invalid:
					self.invalid.append(item)
		return self.invalid

	# ----------------------------------------------------------------------
	# FIXME this check should be done in the invalid
	# ----------------------------------------------------------------------
	def clearUserInvalid(self):
		self._userInvalid = None

	# ----------------------------------------------------------------------
	# FIXME this check should be done in the invalid
	# ----------------------------------------------------------------------
	def addUserInvalid(self, item):
		if self._userInvalid:
			if item not in self._userInvalid:
				self._userInvalid.append(item)
		else:
			self._userInvalid = [item]

		# Append to the invalid cards
		if self.invalid:
			if item not in self.invalid:
				self.invalid.append(item)
		else:
			self.invalid = [item]

	# ----------------------------------------------------------------------
	# @return error message string
	# ----------------------------------------------------------------------
	def errorMessage(self, labels=[]):
		case = self.case()
		extra  = self.info.extra[case]
		range_ = self.info.range[case]
		msg = ""
		displayed = []
		for e in self.invalid:
			if e is None: continue
			if e in displayed: continue	# skip multiple messages
			displayed.append(e)

			try: lab = labels[e]
			except: lab=""

			if msg: msg += "\n"

			if isinstance(e,str):
				msg += "%s"%(e)
			elif lab:
				msg += "What=%d  Label=%s  Type=%s: %s" \
					% (e, lab, range_[e], extra[e])
			else:
				msg += "What=%d  Type=%s: %s" \
					% (e, range_[e], extra[e])
		return msg

	# ----------------------------------------------------------------------
	# Convert card to names or numbers
	# ----------------------------------------------------------------------
	def convert(self, tonames=True):
		if tonames:
			self.info.toNames(self)
		else:
			self.info.toNumbers(self)
		self.setModified()

	# ----------------------------------------------------------------------
	# I/O
	# ----------------------------------------------------------------------
	def whats(self):	return self._what
	def nwhats(self):	return len(self._what)
	def sdum(self):		return str(self.what(0))	# Force always a string
	def extra(self):	return self._extra
	def comment(self):	return self._comment
	def isGeo(self):	return self._geo
	def type(self):		return self._type
	def pos(self):		return self._pos
	def indent(self):	return max(0,self._indent)
	name = sdum

	# ----------------------------------------------------------------------
	# Compare position
	# ----------------------------------------------------------------------
	@staticmethod
	def cmpPos(a, b):	return cmp(a._pos, b._pos)

	# ----------------------------------------------------------------------
	# Return what value -1=extra, 0=sdum
	# ----------------------------------------------------------------------
	def what(self, n):
		if n<0:
			return self._extra
		else:
			try:	return self._what[n]
			except: return ""

	# ----------------------------------------------------------------------
	# Return absolute value of what
	# ----------------------------------------------------------------------
	def absWhat(self, w):
		"""return the absolute value of what[w]"""
		n = self.what(w)
		if n=="" or n==0.0: return n
		try:
			return abs(float(n))
		except:
			n = n.strip()
			if n[0] in ("+", "-"):
				return n[1:]
			return n

	# ----------------------------------------------------------------------
	# Return what as integer
	# ----------------------------------------------------------------------
	def intWhat(self, n, default=0):
		"""try to return an integer from what"""
		return int(self.numWhat(n, default))

	# ----------------------------------------------------------------------
	# Return what as long integer
	# ----------------------------------------------------------------------
	def longWhat(self, n, default=0):
		"""try to return a long integer from what"""
		return int(self.numWhat(n, default))

	# ----------------------------------------------------------------------
	# Return what as floating point number
	# ----------------------------------------------------------------------
	def numWhat(self, n, default=0.0):
		"""try to return a number from what"""
		try:
			return float(self.evalWhat(n))
		except:
			return default

	# ----------------------------------------------------------------------
	# Is evaluated what
	# ----------------------------------------------------------------------
	def isEvalWhat(self, w):
		return isEvalStr(self.what(w))

	# ----------------------------------------------------------------------
	# Return evaluated what
	# ----------------------------------------------------------------------
	def evalWhat(self, n, dollar=True):
		"""return evaluated what if needed"""
		w = self.what(n)
		if isinstance(w,str) and w!="":
			if w[0] in ("=","(","[","{"):
				if self.input is None: return w
				self.input.localDict.card = self
				if w[0]=="=":
					expr = w[1:]
				else:
					expr = w
				# Check for possible vector definitions
				# replace any { to Vector( and } to )
				expr = expr.replace("{","Vector(")
				expr = expr.replace("}",")")
				#expr = expr.replace("[[","Matrix(")
				#expr = expr.replace("]]",")")
				try:
					value = eval(expr, _globalDict, self.input.localDict)
				except:
					say("\nERROR: "+self.rawStr())
					say("ERROR: %s what(%d):%s"%(sys.exc_info()[1],n,str(w)))
					try:
						return self._owhat[n]
					except:
						return "?"

				try:
					f = float(value)
					if float(int(f)) == f:
						return int(f)
					return f
				except:
					return value

			elif dollar:
				if w[0] == "$":
					var  = w[1:]
					sign = 1.0
				elif w[:2] == "-$":
					var  = w[2:]
					sign = -1.0
				elif w[:2] == "+$":
					var  = w[2:]
					sign = +1.0
				else:
					return w

				# check defines
				defines = self.input.cardsCache("#define")
				# dummy linear search
				for define in defines:
					if define.ignore(): continue
					if define.sdum()==var:
						val = define.evalWhat(1)
						# Convert to int or float if possible
						try:
							f = sign*float(val)
							if float(int(f)) == f:
								return int(f)
							return f
						except:
							if sign<0.0:
								return "-"+str(val)
							else:
								return val
		return w

	# ----------------------------------------------------------------------
	# Return evaluated what,
	# but NOT empty space (used on $xxx cards)
	# ----------------------------------------------------------------------
	def evalWhat0(self, n, dollar=True):
		w = self.evalWhat(n, dollar)
		if w=="": return 0.
		else: return w

	# ----------------------------------------------------------------------
	def __delitem__(self, item):
		try:
			del self.prop[item]
		except:
			pass

	# ----------------------------------------------------------------------
	# get card information
	# ----------------------------------------------------------------------
	def __getitem__(self, t=None, default=None):
		"""get card information"""
		if t is None:
			return self.tag
		elif isinstance(t, int):
			return self.what(t)
		elif t in ("sdum","name"):
			return self.what(0)
		elif t=="comment":
			return self._comment
		elif t=="extra":
			return self._extra
		elif t=="info":
			return self.info
		else:
			return self.getProperty(t,default)
	get = __getitem__

	# ----------------------------------------------------------------------
	# set information to card
	# ----------------------------------------------------------------------
	def __setitem__(self, t, v):
		"""set value to what(#, sdum, comment, extra)"""
		if isinstance(t, int):
			self.setWhat(t, v)
		elif t=="sdum":
			self.setSdum(v)
		elif t=="comment":
			self.setComment(v)
		elif t=="extra":
			self.setExtra(v)
		else:
			self.setProperty(t,v)
	set = __setitem__

	# ----------------------------------------------------------------------
	# set card state
	# ----------------------------------------------------------------------
	def setEnable(self, e):
		"""set state enable(True)/disable(False)"""
		self.enable = e
		self.setModified()

	# ----------------------------------------------------------------------
	# Disable card
	# ----------------------------------------------------------------------
	def disable(self):
		"""disable card"""
		self.enable = False
		self.setModified()

	# ----------------------------------------------------------------------
	# Set Active state
	# ----------------------------------------------------------------------
	def setActive(self, a):
		"""set active state of card"""
		self.active = a

	# ----------------------------------------------------------------------
	# Ignore card if not active or not enabled
	# ----------------------------------------------------------------------
	def ignore(self):
		return not self.enable or not self.active

	# ----------------------------------------------------------------------
	def notIgnore(self):
		return self.enable and self.active

	# ----------------------------------------------------------------------
	# set comment
	# ----------------------------------------------------------------------
	def setComment(self, comment=""):
		"""set comment string"""
		self._comment = comment
		self.setModified()

	# ----------------------------------------------------------------------
	# Set the whats list
	# ----------------------------------------------------------------------
	def setWhats(self, whats):
		self._what = whats
		self.setModified()

	# ----------------------------------------------------------------------
	# Set SDUM to card (what[0])
	# ----------------------------------------------------------------------
	def setSdum(self, s):
		"""set sdum"""
		self.setWhat(0, s)
	setName = setSdum

	# ----------------------------------------------------------------------
	# Set extra information to card
	# ----------------------------------------------------------------------
	def setExtra(self, e):
		"""set extra information"""
		try:
			self._extra = str(e)
		except:
			self._extra = e.encode("ascii","replace")
		self.setModified()

	# ----------------------------------------------------------------------
	# Set a what value
	# ----------------------------------------------------------------------
	def setWhat(self, w, value):
		"""set what value"""
		if w<0:
			self._extra = value
		else:
			# Unicode characters are not accepted
			if isinstance(value, str):
				value = value
			elif isinstance(value, int):
				value = str(value)
			elif isinstance(value, float):
				value = str(value)
			else:
				value = value.decode()
			# Extend list if needed
			if w>=len(self._what):
				self._what.extend([""] * (w-len(self._what)+1))
			self._what[w] = value
		self.setModified()

	# ----------------------------------------------------------------------
	# Set the total number of whats
	# Usually to remove whats
	# ----------------------------------------------------------------------
	def setNWhats(self, n):
		"""set the total number of whats"""
		l = len(self._what)
		if n < l:
			del self._what[n:]
		elif n > l:
			self._what += [""] * (n-l)
		self.setModified()

	# ----------------------------------------------------------------------
	# Set absolute information to card. The "what" will keep the same sign
	# ----------------------------------------------------------------------
	def setAbsWhat(self, w, v):
		"""change the absolute value of what keeping the same sign"""
		s = self.sign(w)

		if v=="" or v==0.0:
			self.setWhat(w, v)
			return

		if not isEvalStr(v):
			s = self.sign(w)
			try:
				v = abs(float(v))
				if s: v = -v

			except:
				if v[0]=="+" or v[0]=="-":
					if s: v = "-" + v[1:]
					else: v = v[1:]
				else:
					if s: v = "-" + v

		self.setWhat(w, v)

	# ----------------------------------------------------------------------
	# Signs are used to avoid cases when the what[w] is ZERO
	# Should not be used by the user... only in Layout.py
	# ----------------------------------------------------------------------
	def setSign(self, w, s):
		if isEvalStr(self.what(w)): return

		if w >= len(self._sign):
			self._sign += [None] * (w-len(self._sign)+1)

		self._sign[w] = s
		n = self.what(w)

		if n=="" or n==0.0: return
		try:
			n = abs(float(n))
			if s: n = -n
		except:
			n = n.strip()
			if n[0] in ("+", "-"):
				if s: n = "-" + n[1:]
				else: n = n[1:]
			else:
				if s: n = "-" + n
		self.setWhat(w, n)

	# ----------------------------------------------------------------------
	# return (minus) sign of a what
	# ----------------------------------------------------------------------
	def sign(self, w):
		"""return (minus) sign of a number: False=+, True=-"""
		if isEvalStr(w):
			# special case
			return False

		# Default numeric value
		try:
			n = self._what[w]
		except:
			try: return self._sign[w]
			except: return False

		if n=="": n=0.0
		try:
			n = float(n)
			if n<0.0: return True
			if n>0.0: return False

			if w<len(self._sign):
				s = self._sign[w]
				if s is not None: return s
			return False
		except:
			n = n.strip()
			if n[0] == "-": return True
			if n[0] == "+": return False

			if w<len(self._sign):
				s = self._sign[w]
				if s is not None: return s
			return False

	# ----------------------------------------------------------------------
	# set/get user property
	# ----------------------------------------------------------------------
	def setProperty(self, t, v):
		if self.prop is None:
			self.prop = {}
		self.prop[t] = v

	# ----------------------------------------------------------------------
	def getProperty(self, t, default=None):
		if self.prop is None: return default
		return self.prop.get(t, default)

	# ----------------------------------------------------------------------
	def clearProperties(self, all=True):
		if self.prop is None: return
		if all:
			self.prop = None
		else:
			# Clear only user properties,
			# ignore the system ones starting with @
			for n in list(self.prop.keys()):
				if n[0]!="@":
					del self.prop[n]

	# ----------------------------------------------------------------------
	# Return a specific type from the card, see database from more info
	# ----------------------------------------------------------------------
	def findType(self, target):
		"""return whats with a specific type from the card, pi, mi, ..."""
		case = self.info.findCase(self)
		lst  = self.info.find(target, case)
		if isinstance(lst, list):
			return [self.what(x) for x in lst]
		else:
			return (self.what(lst[0]), self.what(lst[1]), lst[2])

	# ----------------------------------------------------------------------
	# Return units used by card
	# ----------------------------------------------------------------------
	def units(self, absolute=True):
		"""return unit from card"""
		case = self.info.findCase(self)
		lst  = self.info.find("lu", case)
		if absolute:
			if isinstance(lst, list):
				return [abs(self.intWhat(x)) for x in lst]
			else:
				return (abs(self.intWhat(lst[0])),
					abs(self.intWhat(lst[1])), lst[2])
		else:
			if isinstance(lst, list):
				return [self.intWhat(x) for x in lst]
			else:
				return (self.intWhat(lst[0]),
					self.intWhat(lst[1]), lst[2])

	# ----------------------------------------------------------------------
	# Return a formated comment string
	# ----------------------------------------------------------------------
	def commentStr(self):
		"""return a comment string"""
		line = ""
		for l in self._comment.splitlines():
			if len(l)>0: l = " "+l
			line += "*%s\n"%(l.rstrip())
		if self._comment!="" and self._comment[-1]=="\n": line += "*\n"
		return line

	# ----------------------------------------------------------------------
	# return bodies in single/double/free format
	# ----------------------------------------------------------------------
	def _bodyStr(self, fmt=None, prefix=""):
		i = 1
		if fmt is None or fmt==FORMAT_FREE:
			name = str(self.what(0)).strip()
			if name == "":
				name = "<name>"
			line = "%s%s %-10s" % (prefix, self.tag, name)
			prevlen = 0
			while i<self.info.nwhats:
				# use what to keep it as string
				w = self.evalWhat(i,False)
				try:
					num = float(w)	# Force to number
					w = bmath.format(w, 22)
				except:
					if not w or (w[0]!="$" and w[:2]!="-$" and w[:2]!="+$"):
						w = "0.0"
				if len(line)+len(w)+1 > prevlen+80:
					prevlen = len(line)
					line += "\n%s              %s"%(prefix,w)
				else:
					line += " " + w
				i += 1
			return line

		elif fmt==FORMAT_SINGLE:
			width = 10
			step  = 6

		else:
			width = 21
			step  = 3

		line = "%s  %3s%5s" % (prefix, self.tag, str(self.what(0)))
		while i<len(self._what):
			inext = i+step
			imin  = min(inext, len(self._what))
			while i<imin:
				w = self._what[i]
				try:
					num = float(w)
					w = bmath.format(w, width)
				except:
					w = "0.0"
				line += w.rjust(width)
				i += 1
			if i<len(self._what): line += "\n%s          "%(prefix)
		return line

	# ----------------------------------------------------------------------
	# return end string
	# ----------------------------------------------------------------------
	def _endStr(self, fmt=None):
		if fmt is None or fmt==FORMAT_FREE:
			return "END"
		else:
			return "  END"

	# ----------------------------------------------------------------------
	# Split a region expression into multi lines
	# ----------------------------------------------------------------------
	@staticmethod
	def splitExpr(expr, maxlength=100):
		line = ""
		first = True
		for l in expr.splitlines():
			l = l.rstrip()	# remove trailing spaces
			if len(l)==0: continue
			while len(l)>maxlength:
				for i in range(maxlength-1, 1, -1):
					if l[i] not in ('+','-','|','(',')','#'):
						continue
					# split line on char before
					if first:
						line += l[0:i].rstrip()
						first = False
					else:
						line += "\n%s" % (l[0:i])
					l = l[i:]
					break
				if i<=2: break
			if first:
				line += l
				first = False
			else:
				line += "\n%s" % (l)
		return line

	# ----------------------------------------------------------------------
	# return region string
	# ----------------------------------------------------------------------
	def _regionStr(self, fmt=None):
		try: neighbors = int(self._what[1])
		except: neighbors = 5

		if fmt is None or fmt==FORMAT_FREE:
			name = str(self.what(0)).strip()
			if name == "":
				name = "<name>"

			if name == "&":	# Continuation of a region
				line = "              "
			else:
				line = "%-10s%4d " % (name, neighbors)

			first = True
			for l in self._extra.splitlines():
				l = l.rstrip()	# remove trailing spaces
				if len(l)==0: continue
				while len(l)>100:
					for i in range(99, 1, -1):
						if l[i] not in ('+','-','|','(',')','#'):
							continue
						# split line on char before
						if first:
							line += l[0:i].rstrip()
							first = False
						else:
							line += "\n               %s" \
								% (l[0:i])
						l = l[i:]
						break
					if i<=2: break
				if first:
					line += l
					first = False
				else:
					line += "\n               %s" % (l)
			return line
		else:
			expr = ""
			first = True
			if self._extra.find("(")>=0 or self._extra.find(")")>=0:
				raise TypeError("No parenthesis are allowed in fixed format")

			nop = 0
			for l in self._extra.splitlines():
				# Remove spaces
				l = " ".join(l, '')
				l = l.replace("+", " +")
				l = l.replace("-", " -")
				l = l.replace("|", " | ")

				if nop>0 and expr!="":
					expr += "\n          "
					nop = 0

				orbefore = False
				for op in l.split():
					if op=="|":
						expr += "OR"
						orbefore = True
						continue
					else:
						if not orbefore: expr += "  "
						orbefore = False
					sop = str(op).rjust(5)
					if len(sop)>5: sop = sop[0:5]
					expr += sop
					nop  += 1
					if nop==9:
						expr += "\n          "
						nop = 0

			return "%5s%5d%s" % (str(self.what(0)), neighbors, expr)

	# ----------------------------------------------------------------------
	def _toStr(self, fmt=None, prefix=""):
		"""return a card string"""

		# special cards
		ctag = prefix + self.tag
		try:
			if self.tag[0] == "#":
				if self.nwhats()>0:
					return "%s %s" % \
						(ctag, " ".join(
							[str(self.evalWhat(i,False))
								for i in range(0,self.nwhats())]))
				else:
					return ctag
			elif self.tag[0] == "$":
				if self.nwhats()>0:
					return "%s %s" % \
						(ctag, " ".join(
							[str(self.evalWhat0(i,False))
								for i in range(1,self.nwhats())]))
				else:
					return ctag
		except IndexError:
			pass

		if self._geo:
			if self.tag == "END":
				return self._endStr()
			elif self._type == Card.BODY and self.tag!="VOXELS":
				return self._bodyStr(fmt,prefix)
			elif self._type == Card.REGION:
				return self._regionStr(fmt)

		if fmt==FORMAT_FREE:
			width = 22
		else:
			width = 10

		# Number of whats
		nwhats = len(self._what)
		if self.tag == "PLOTGEOM":
			#Another card with special treatment
			nwhats = min(6,nwhats)
		elif fmt==FORMAT_FREE:
			# Count real whats
			while (nwhats>0):
				i = nwhats-1
				if self._what[i] == "":
					nwhats = i
				else:
					break

		# Add whats in blocks of 6
		line   = ""
		first  = True
		i      = 1

		if prefix!="":
			ctag = "%s%-10s"%(prefix,self.tag)
		else:
			ctag = self.tag.ljust(10)
		# Cards with no whats
		if nwhats <= 1:
			if len(line)>0: line += "\n"
			line += ctag
			if nwhats == 1:
				if fmt==FORMAT_FREE:
					line += (" , "*7) + str(self.evalWhat(0,False))
				else:
					line += (" "*(width*6)) + str(self.evalWhat(0,False))
			line = line.rstrip()
		else:	# Cards with many whats
			if self.tag == "COMPOUND":
				cont_sdum = str(self.evalWhat(0,False))
			else:
				cont_sdum = " &"

			while i < nwhats:
				if len(line)>0: line += "\n"
				line += ctag
				if fmt==FORMAT_FREE:
					line += ", "

				inext = i+6
				imin  = min(inext, nwhats)
				while i<imin:
					w = self.evalWhat(i,False)
					if isinstance(w, int):	# for RADDECAY
						w = str(w)
					else:
						try:
							num = float(w)
							w = bmath.format(w, width)
						except:
							w = str(w)
					if fmt==FORMAT_FREE:
						line += w + ", "
					else:
						line += w.rjust(width)
					i += 1

				# Check for continuation lines and sdum
				if first:
					sdum = str(self.evalWhat(0,False))
				elif i>6:
					sdum = cont_sdum
					if sdum.strip()=="&": cont_sdum += "&"
				first = False

				# Do we need to pad with the line
				if sdum != "":
					if i<inext:
						if fmt==FORMAT_FREE:
							line += " ,"*(inext-i)
						else:
							line += " "*(width*(inext-i))
					line += sdum
				elif fmt==FORMAT_FREE:
					# Chop the last comma
					line = line.rstrip()
					if line[-1] == ',':
						line = line[:-1]
				line = line.rstrip()

		# XXX Double format doesn't exist for GEOBEGIN
		# simply triggers the writing of the header inline the input
		# or in an external file
		if self.tag in ("TITLE", "OPEN", "!mesh"):
			line = "%s\n%s" % (line, self._extra.rstrip())

		elif self.tag == "GEOBEGIN":
			try: iv = self.ivopt
			except: self.ivopt = 0
			try: id = self.idbg
			except: self.idbg = 0
			if fmt==FORMAT_FREE:
				line += "\n%d,%d,%s" \
					% (self.ivopt, self.idbg, self._extra.rstrip())
			else:
				line += "\n%5d%5d          %s" \
					% (self.ivopt, self.idbg, self._extra.rstrip())

		elif self.tag == "PLOTGEOM":
			w6 = self.intWhat(6)
			if w6 == 0 or w6 == 5:
				line += "\n%s" % (self._extra)
				line += "\n%10s%10s%10s%10s%10s%10s" \
					% (bmath.format(self.numWhat( 7), 10),
					   bmath.format(self.numWhat( 8), 10),
					   bmath.format(self.numWhat( 9), 10),
					   bmath.format(self.numWhat(10), 10),
					   bmath.format(self.numWhat(11), 10),
					   bmath.format(self.numWhat(12), 10) )
				line += "\n%10s%10s%10s%10s%10s%10s" \
					% (bmath.format(self.numWhat(13), 10),
					   bmath.format(self.numWhat(14), 10),
					   bmath.format(self.numWhat(15), 10),
					   bmath.format(self.numWhat(16), 10),
					   bmath.format(self.numWhat(17), 10),
					   bmath.format(self.numWhat(18), 10) )
				line += "\n%10s%10s%10s%10s" \
					% (bmath.format(self.numWhat(19), 10),
					   bmath.format(self.numWhat(20), 10),
					   bmath.format(self.numWhat(21), 10),
					   bmath.format(self.numWhat(22), 10) )
		return line

	# ---------------------------------------------------------------------
	# return string with only the cards
	# ----------------------------------------------------------------------
	def toStr(self, fmt=None):
		if not self.enable and (self.info.disableComment and commentedCards):
			prefix = "*"
		else:
			prefix = ""

		s = self._toStr(fmt, prefix)
		#if not self.enable and self.info.disableComment: s = "*"+s
		if isinstance(s,str):
			return s.encode('utf-8')
		else:
			return s

	# ----------------------------------------------------------------------
	def evalWhatStr(self):
		"""return string for special evaluated whats"""
		es = ""
		for i,w in enumerate(self._what):
			# Same as isEvalStr but with the "$" search
			if isinstance(w,str) and w!="" and w[0] == "=":
				es += "!@what.%d%s\n"%(i,w)

		if self.prop is not None:
			for n,v in list(self.prop.items()):
				if n[0] != '@' and v:	# skip system values
					es += "!@%s=%s\n"%(n,str(v))
		return es

	# ----------------------------------------------------------------------
	# Return card without evaluating whats
	# ----------------------------------------------------------------------
	def rawStr(self):
		w = [str(x) for x in self._what]
		w.insert(0,self.tag)
		w.append(self._extra)
		s = " ".join(w)
		if self._comment!="":
			return self.commentStr() + s
		else:
			return s

	# ----------------------------------------------------------------------
	# Convert to string using default format
	# ----------------------------------------------------------------------
	def __str__(self):
		"""return a card string with the comment"""
		if self._comment!="":
#			return self.commentStr() + self.toStr()
			return self.commentStr()+"".join( chr(x) for x in self.toStr())
		else:
			return "".join( chr(x) for x in self.toStr())
#			return self.toStr()

	# ----------------------------------------------------------------------
	# Dump card to pickler
	# ----------------------------------------------------------------------
	def dump(self, pickler):
		pickler.dump(self.tag)
		pickler.dump(self._what)
		pickler.dump(self._extra)
		pickler.dump(self._comment)
		if self.prop:	# Skip system variables
			pickler.dump([x
				for x in list(self.prop.items())
#					if type(x[1]) in (IntType, LongType, FloatType,
#						BooleanType, StringType, UnicodeType)])
					if type(x[1]) in (int, float,
						          bool, str, bytes )])
		else:
			pickler.dump(None)
		pickler.dump(self.enable)

	# ----------------------------------------------------------------------
	# Load card from unpickler
	# ----------------------------------------------------------------------
	def load(self, unpickler):
		self.changeTag( unpickler.load())
		self._what    = unpickler.load()
		self._extra   = unpickler.load()
		self._comment = unpickler.load()
		prop          = unpickler.load()
		self.enable   = unpickler.load()
		if prop:
			for n,v in prop:
				self[n] = v

	# ----------------------------------------------------------------------
	# Return differences of the two cards as a list of integers
	#       -3 = enable/disable
	#	-2 = Comment
	#	-1 = Extra
	#	 0 = Sdum
	#        1..n = What's
	# ----------------------------------------------------------------------
	def diff(self, card):
		df = []
		if self.enable != card.enable:
			df.append(-3)
		if self.commentStr() != card.commentStr():
			df.append(-2)
		if self._extra != card._extra:
			df.append(-1)
		minwhats = min(self.nwhats(), card.nwhats())
		maxwhats = max(self.nwhats(), card.nwhats())
		for w in range(0, minwhats):
			if w > card.nwhats():
				df.append(w)
			elif self.what(w) != card.what(w):
				df.append(w)
		if minwhats != maxwhats:
			df.extend(list(range(minwhats, maxwhats)))
		return df

	# ----------------------------------------------------------------------
	# Special body method returning various geometrical information
	# ... if they are defined
	# ----------------------------------------------------------------------
	def bodyP(self):
		"""Position of body"""
		if   self.tag == "RPP":
			return bmath.Vector(	self.numWhat(1),
						self.numWhat(3),
						self.numWhat(5))
		elif self.tag == "PLA":
			return bmath.Vector(	self.numWhat(4),
						self.numWhat(5),
						self.numWhat(6))
		elif self.tag == "XYP":
			return bmath.Vector(0.0, 0.0, self.numWhat(1))
		elif self.tag == "XZP":
			return bmath.Vector(0.0, self.numWhat(1), 0.0)
		elif self.tag == "YZP":
			return bmath.Vector(self.numWhat(1), 0.0, 0.0)
		elif self.tag in ("XCC", "XEC"):
			return bmath.Vector(0.0, self.numWhat(1), self.numWhat(2))
		elif self.tag in ("YCC", "YEC"):
			return bmath.Vector(self.numWhat(2), 0.0, self.numWhat(1))
		elif self.tag in ("ZCC", "ZEC"):
			return bmath.Vector(self.numWhat(1), self.numWhat(2), 0.0)
		elif self.tag in ("BOX", "SPH", "RCC", "TRC", "WED", "REC", "ELL", "ARB") or \
		     self.tag[0] == "!":
			return bmath.Vector(	self.numWhat(1),
						self.numWhat(2),
						self.numWhat(3))
		else:
			return None

	bodyP1 = bodyP

	# ----------------------------------------------------------------------
	def bodyP2(self):
		if   self.tag in ("ELL", "ARB"):
			return bmath.Vector(	self.numWhat(4),
						self.numWhat(5),
						self.numWhat(6))
		else:
			return None

	# ----------------------------------------------------------------------
	# WARNING index starting from 1
	# ----------------------------------------------------------------------
	def bodyPn(self, n):
		if   self.tag == "ARB":
			i = 3*n-2
			return bmath.Vector(	self.numWhat(i),
						self.numWhat(i+1),
						self.numWhat(i+2))
		else:
			return None

	# ----------------------------------------------------------------------
	def bodyN(self):
		"""Normal of body"""
		if   self.tag == "PLA":
			return bmath.Vector(	self.numWhat(1),
						self.numWhat(2),
						self.numWhat(3))
		elif self.tag in ("XYP", "ZCC", "ZEC"):
			return bmath.Vector.Z
		elif self.tag in ("XZP", "YCC", "YEC"):
			return bmath.Vector.Y
		elif self.tag in ("YZP", "XCC", "XEC"):
			return bmath.Vector.X
		else:
			return None

	# ----------------------------------------------------------------------
	def bodyX(self):
		"""X-vector of body"""
		if   self.tag == "RPP":
			return bmath.Vector(self.numWhat(2)-self.numWhat(1), 0.0, 0.0)
		elif self.tag in ("BOX", "WED"):
			return bmath.Vector(	self.numWhat(4),
						self.numWhat(5),
						self.numWhat(6))
		elif self.tag == "REC":
			return bmath.Vector(	self.numWhat(7),
						self.numWhat(8),
						self.numWhat(9))
		elif self.tag == "XEC":
			return bmath.Vector(0.0, self.numWhat(3), 0.0)
		elif self.tag == "YEC":
			return bmath.Vector(0.0, 0.0, self.numWhat(3))
		elif self.tag == "ZEC":
			return bmath.Vector(self.numWhat(3), 0.0, 0.0)
		else:
			return None

	# ----------------------------------------------------------------------
	def bodyY(self):
		"""Y-vector of body"""
		if   self.tag == "RPP":
			return bmath.Vector(0.0, self.numWhat(4)-self.numWhat(3), 0.0)
		elif self.tag in ("BOX", "WED"):
			return bmath.Vector(	self.numWhat(7),
						self.numWhat(8),
						self.numWhat(9))
		elif self.tag == "REC":
			return bmath.Vector(	self.numWhat(10),
						self.numWhat(11),
						self.numWhat(12))
		elif self.tag == "XEC":
			return bmath.Vector(0.0, 0.0, self.numWhat(4))
		elif self.tag == "YEC":
			return bmath.Vector(self.numWhat(4), 0.0, 0.0)
		elif self.tag == "ZEC":
			return bmath.Vector(0.0, self.numWhat(4), 0.0)
		else:
			return None

	# ----------------------------------------------------------------------
	def bodyZ(self):
		"""Z-vector of body"""
		if   self.tag == "RPP":
			return bmath.Vector(0.0, 0.0, self.numWhat(6)-self.numWhat(5))
		elif self.tag in ("BOX", "WED"):
			return bmath.Vector(	self.numWhat(10),
						self.numWhat(11),
						self.numWhat(12))
		elif self.tag in ("RCC", "TRC", "REC"):
			return bmath.Vector(	self.numWhat(4),
						self.numWhat(5),
						self.numWhat(6))
		elif self.tag[0] == "!":
			return bmath.Vector(	self.numWhat(7),
						self.numWhat(8),
						self.numWhat(9))
		else:
			return None

	# ----------------------------------------------------------------------
	def bodyR(self):
		"""Radius (or x-radius) of body"""
		if   self.tag == "SPH":
			return self.numWhat(4)
		elif self.tag == "RCC":
			return self.numWhat(7)
		elif self.tag in ("XCC", "YCC", "ZCC"):
			return self.numWhat(3)
		elif self.tag == "XEC":
			return bmath.Vector(0.0, self.numWhat(3), self.numWhat(4))
		elif self.tag == "YEC":
			return bmath.Vector(self.numWhat(4), 0.0, self.numWhat(3))
		elif self.tag == "ZEC":
			return bmath.Vector(self.numWhat(3), self.numWhat(4), 0.0)
		else:
			return None

	# ----------------------------------------------------------------------
	# Add a Zone to REGION
	# ----------------------------------------------------------------------
	def addZone(self, zone):
		"""Add zone to card"""
		try:
			zone = str(zone)
		except:
			zone = zone.encode("ascii","replace")

		if self._extra.strip():
			self._extra += "\n| %s"%(zone)
		else:
			self._extra = zone
		self.setModified()

	# ----------------------------------------------------------------------
	# Set a list of zones as expression
	# ... removing empty zones
	# ----------------------------------------------------------------------
	def setZones(self, zones):
		if zones:
			self.setExtra("\n| ".join(
				[x for x in map(str.strip,zones) if x!=""]))
		else:
			self.setExtra("")

	# ----------------------------------------------------------------------
	# Load voxel information
	# ----------------------------------------------------------------------
	def loadVoxel(self):
		if self.tag != "VOXELS": return
		self["@voxel"] = Voxel("%s.vxl"%(self.sdum().strip()))

#===============================================================================
# Input class
#===============================================================================
class Input:
	def __init__(self, filename=None):
		# Input file information
		self.filename    = ""		# Input Filename
		self.geoFile     = ""
		self.geoOutFile  = ""
		self.format      = FORMAT_SINGLE
		self.geoFormat   = FORMAT_SINGLE

		# Logical Units
		self.units       = LogicalUnits()

		# Cards
		self.cards       = {}
		self.cardlist    = []

		# Verbosity level (False|True)
		self.verbose     = False

		# parsing initialization
		self._cardEnable     = True
		self._commentDisable = False
		self._files          = []	# File handles used during I/O
		self._lineNo         = []	# file line number
		self._filesType      = []	# card opened file type
		self._cwhat          = []	# Commented special whats
		self._prop           = []	# Commented card properties
		self._comment        = ""	# Comment preceding card
		self._nwhats         = 0	# Number of whats of last added card

		# cache card lists
		self.cache = {}
		self.setModified()
		self.localDict = LocalDict(self)
		if filename is not None:
			self.read(filename)

	# ----------------------------------------------------------------------
	def setModified(self):
		"""set last time modified"""
		self._modified = time.time()

	# ----------------------------------------------------------------------
	def clone(self):
		"""return a clone of current input"""
		inp = Input()

		inp.filename   = self.filename
		inp.geoFile    = self.geoFile
		inp.geoOutFile = self.geoOutFile

		for card in self.cardlist:
			inp.addCard(card.clone())

		inp.renumber()
		return inp

	# ----------------------------------------------------------------------
	def __getitem__(self, item):
		return self.cards.get(item,[])

	# ----------------------------------------------------------------------
	# Create bare minimum input file
	# ----------------------------------------------------------------------
	def minimumInput(self):
		self.addCard(Card("GEOBEGIN",["COMBNAME"]))
		self.addCard(Card("END"))
		self.addCard(Card("END"))
		self.addCard(Card("GEOEND"))
		self.addCard(Card("STOP"),None)
		self.renumber()

	# ----------------------------------------------------------------------
	# Set as input/modified time the latest time from files
	# ----------------------------------------------------------------------
	def setFileTime(self):
		for filename in self.filenames():
			try:
				mtime = os.stat(filename).st_mtime
				self._modified = max(self._modified, mtime)
			except OSError:
				pass

	# ----------------------------------------------------------------------
	# Check if a new input file exists
	# ----------------------------------------------------------------------
	def checkInputFile(self):
		# FIXME Scan all #include cards
		for filename in self.filenames():
			try:
				if os.stat(filename).st_mtime > self._modified:
					return True
			except:
				pass
		return False

	#-----------------------------------------------------------------------
	# Return input filename
	#-----------------------------------------------------------------------
	def filenames(self):
		files = []
		if self.filename:	files.append(self.filename)
		if self.geoFile:	files.append(self.geoFile)
		files.extend([card.sdum() for card in self["#include"]])
		return files

	#-----------------------------------------------------------------------
	# open file and append it to the file handle list
	#-----------------------------------------------------------------------
	def _openFile(self, filename, mode="r", tag=None, backup=True):
		# check for infinite loop!
		if filename is None:
			f = None	# don't open anything, useful for #include which are disabled
		elif isinstance(filename, str) or isinstance(filename, str):
			for fn in self._files:
				if fn.name==filename:
					say("ERROR: File '%s' is recursively loaded"%(filename))
					return None

			if mode=="w" and backup:
				backupname = filename+"~"
				try: os.remove(backupname)
				except: pass
				try: os.rename(filename, backupname)
				except: pass

			try:	# open file
				f = open(filename, mode)
			except:
				say("ERROR: Cannot open file '%s'"%(filename))
				f = open("/dev/null",mode)
		else:
			f = filename

		self._files.append(f)
		self._lineNo.append(0)
		self._filesType.append(tag)
		return f

	#-----------------------------------------------------------------------
	def _closeFile(self):
		if self._files:
			f = self._files.pop()
			if f: f.close()
			self._lineNo.pop()
			self._filesType.pop()

	#-----------------------------------------------------------------------
	# read line from the last opened file
	#-----------------------------------------------------------------------
	def _readLine(self, closeinclude=True):
		while True:
			try:
				self._lineNo[-1] += 1
				line = self._files[-1].readline()
			except IndexError:
				return None
			except AttributeError:
				line = ""	# from None file

			if len(line)==0:
				# end of file?
				if self._filesType[-1] == "#include":
					if closeinclude:
						card = Card("#endinclude", comment=self._comment)
						self._addCard(card)
						if self._files[-1] is None:
							card.disable()
					else:
						return "#endinclude"
				self._closeFile()
				if len(self._files)==0: return None
				continue

			line = line.rstrip()
#			try:
#				uline = line.decode("utf-8")
#			except (UnicodeDecodeError, UnicodeEncodeError):
			return line
#			try:
#				sline = str(uline)
#				return line
#			except (UnicodeDecodeError, UnicodeEncodeError):
#				return uline

	#-----------------------------------------------------------------------
	# Parse a comment line
	#-----------------------------------------------------------------------
	def _parseComment(self, line, check4Card=True):
		line = line.rstrip()
		if len(line)==0: return None
		if line[0] == "*":
			# if next character is not space
			# check if it is a valid card
			cstart = 1
			if line!="*":
				if line[1]==" ": cstart = 2

				if line[1:4] in ("...", "234", "===", "---",
							"***", "___", "_..", "+++", "_AA" ):
					cstart = 2

				elif line[1:7] == "@what.":
					# treat special commented whats
					m = _CWHAT.match(line)
					if m:
						try: w = int(m.group(1))
						except: w = -1
						if w>=0:
							if w==len(self._cwhat):
								self._cwhat.append(m.group(2))
							else:
								if w>len(self._cwhat):
									self._cwhat.extend([""]*(w-len(self._cwhat)+1))
								self._cwhat[w] = m.group(2)
							return None

				elif commentedCards and (check4Card or line[1]=="$"):
					# Check for commented cards
					try:
						if line[cstart]==" ":
							lline = line[cstart+1:]
						else:
							lline = line[cstart:]
						(tag, what) = self._parseLine(lline, False)
						# Check tag if it exists in CardInfo
						ci = CardInfo.get(tag)
						if ci.tag == tag:
							self._cardEnable = False
							self._commentDisable = True
							return lline
					except:
						pass

			cline = line[cstart:]
			if self._comment!="":
				self._comment += "\n"+cline
			else:
				if cline=="":
					self._comment = " "
				else:
					self._comment = cline
			return None

		else:
			if self._commentDisable:
				self._cardEnable = True
				self._commentDisable = False

		if line[0] == "#":
			if line=="#if 0":
				self._cardEnable = False
			elif line=="#endif" and not self._cardEnable:
				self._cardEnable = True
			else:
				return line
			return None
		else:
			# Remove inline comments
			p = line.find("!")
			if p==0:
				if line[:8].strip() in FLAIR_TAGS:
					return line
				elif line[1] == "@":
					# treat special commented whats
					m = _CWHAT.match(line)
					if m:
						try: w = int(m.group(1))
						except: w = -1
						if w>=0:
							if w==len(self._cwhat):
								self._cwhat.append(m.group(2))
							else:
								if w>len(self._cwhat):
									self._cwhat.extend([""]*(w-len(self._cwhat)+1))
								self._cwhat[w] = m.group(2)
							return None

					m = _OPT.match(line)
					if m:
						self._prop.append((m.group(1),m.group(2)))
						return None

					# if everything fails add it as comment
					if self._comment: self._comment += "\n"
					self._comment += line
					return None

			elif p>=0:
				comment = line[p+1:]
				if self._comment: self._comment += "\n"
				self._comment += comment
				if p == 0: return None
				line = line[0:p]

			return line

	#-----------------------------------------------------------------------
	# Parse a single line to card, whats and sdum
	#-----------------------------------------------------------------------
	def _parseLine(self, line, enable=True):
		if line[0] == "$":
			what = self._parseFreeFormat(line)
			tag = what.pop(0).lower()
			what.insert(0, "")

		elif line[0] == "#":
			m = _PRPPAT.match(line)
			if m:
				tag = m.group(1)
				if tag == "#include":
					# replace the whats
					what = [line[9:].strip()]
					if self._cardEnable and enable:
						self._openFile(what[0],"r",tag)
					else:
						self._openFile(None,"r",tag)
				else:
					what = [m.group(2), m.group(3)]
			else:
				tag = line.strip()
				what = [""]

		elif self.format == FORMAT_SINGLE:
			tag = line[0:8].rstrip()
			what = [line[70:].strip()]	# sdum
			for i in range(10, 61, 10):
				w = line[i:i+10].strip()
				what.append(_str2num(w))

		else:
			what = self._parseFreeFormat(line)
			tag = what.pop(0)
			if len(tag)>8: tag=tag[0:8]
			if len(what)==7:
				sdum = what.pop()
			else:
				sdum = ""
			what.insert(0, sdum)

		return (tag, what)

	#-----------------------------------------------------------------------
	# Parse a free format line
	# Remove unwanted characters from free format line and
	# separate everything with comma ,
	#-----------------------------------------------------------------------
	def _parseFreeFormat(self, line):
		"""Remove unwanted characters from a free format line
		and separate everything with comma ,"""
		line = re.sub(" *[;:/,] *",",",line.strip())
		what = re.sub(" +",",",line).split(',')
		return list(map(_str2num, what))

	#-----------------------------------------------------------------------
	# Parse needed whats
	#-----------------------------------------------------------------------
	def _whatNeeded(self, fgeo, what, needwhat):
		# Read whats
		while len(what)<needwhat:
			self._lineNo += 1
			line = utfReadline(fgeo)
			if len(line)==0:
				raise Exception("Error paring geometry")

			line = self._parseComment(line, False)
			if line is None: continue

			self._parseBodyLine(line, what)

		# delete extra whats if any
		if len(what)>needwhat: del what[needwhat:]

	#-----------------------------------------------------------------------
	# Import/parse a file including geometry
	#-----------------------------------------------------------------------
	def parse(self):
		"""Import/parse a fluka file"""
		tag  = ""
		sdum = ""
		self._comment = ""		# current comment
		self._cwhat   = []		# Commented special whats
		self._nwhats  = 0		# Number of whats of last added card
		self._cardEnable = True
		while True:
			line = self._readLine()
			if line is None: return
			line = self._parseComment(line)
			if line is None: continue

			prev_tag  = tag
			prev_sdum = sdum
			tag, what = self._parseLine(line)
			sdum  = what[0]
			extra = ""

			# Handle card exceptions
			# 1. Continuation
			if isinstance(sdum,str) and sdum.find("&")>=0 and tag==prev_tag:
				# Find last card and append whats
				try:
					card = self.cards[tag][-1]
					card.appendComment(self._comment)
					if card.nwhats() == self._nwhats:
						# append cards at the end
						card.appendWhats(what[1:])
					else:
						# append cards in the middle in case
						# that functions are used
						card.appendWhats(what[1:], self._nwhats)
					self._nwhats += len(what)-1
					self._comment = ""
					continue
				except KeyError:
					pass

			# 2. Compound continuation (special case)
			elif tag == "COMPOUND":
				# check if previous line is COMPOUND same state
				# and same sdum
				if prev_tag == tag and prev_sdum == sdum:
					# find last COMPOUND card with that sdum
					card = self.cards[tag][-1]
					# check state
					if card.enable == self._cardEnable:
						card.appendComment(self._comment)
						if card.nwhats() == self._nwhats:
							card.appendWhats(what[1:])
						else:
							card.appendWhats(what[1:], self._nwhats)
						self._nwhats += len(what)-1
						self._comment = ""

						self.checkCompound(card)
						continue	# skip the rest

			# 3. Extra line readout
			if tag in ("TITLE", "OPEN", "!mesh"):
				extra = self._readLine()
				if extra is None: return

			# Handling of special cases needed for parsing
			# the rest of the file
			elif tag == "FREE":
				if self._cardEnable:
					self.format = FORMAT_FREE
				#else:
				#	self.commentFormat = FORMAT_FREE
				if self.verbose: say("FREE: Input format changed to", self.format)

			elif tag == "FIXED":
				if self._cardEnable:
					self.format = FORMAT_SINGLE
				#else:
				#	self.commentFormat = FORMAT_SINGLE
				if self.verbose: say("FIXED: Input format changed to", self.format)

			elif tag=="GLOBAL" and self._cardEnable:
				w4 = _intWhat(what, 4)
				if w4==2 or w4==3:
					self.format= FORMAT_FREE
				if _intWhat(what, 5)>0:
					self.geoFormat = FORMAT_FREE

			# Handle special card PLOTGEOM
			elif tag == "PLOTGEOM":
				w6 = _intWhat(what, 6)
				if self._openFile("PLG.GFXINDAT","r",tag) is None:
					say("ERROR: cannot open file PLG.GFXINDAT")
					continue

				extra = self._readLine()
				line  = self._readLine()
				what.extend([_str2num(line[x:x+10]) \
					for x in range(0, 60, 10)])

				line = self._readLine()
				what.extend([_str2num(line[x:x+10]) \
					for x in range(0, 60, 10)])

				line = self._readLine()
				what.extend([_str2num(line[x:x+10]) \
					for x in range(0, 40, 10)])

			elif tag=="COMMENT":
				# Add to comments the following lines
				for i in range(max(1, _intWhat(what, 1))):
					line = self._parseComment("*"+self._readLine())
				continue

			elif tag == "GEOBEGIN" and sdum!="FLUGG":
				# Starting position of GEO
				if self.geoFormat == FORMAT_SINGLE:
					if sdum[:8]=="COMBNAME":
						self.geoFormat = FORMAT_FREE
				card = Card(tag, what, self._comment)
				self._addCard(card)

				w3 = _intWhat(what, 3)
				w4 = _intWhat(what, 4)

				if w3 not in (0,5):
					geoFile = self._readLine()
					if geoFile is None:
						raise Exception("Error parsing geometry")
					self.geoFile = geoFile.strip()	# XXX???

				if w4 not in (0,11):
					geoOutFile = self._readLine()
					if geoOutFile is None:
						raise Exception("Error parsing geometry")
					self.geoOutFile = geoOutFile.strip()

				if self.geoFile != "":
					fgeo = self._openFile(self.geoFile,"r",tag)
				self._parseGeometry(card)
				continue

			elif tag == "LATTICE":
				self._parseLattice(what)

			elif tag == "LATTSNGL":
				try:
					what[2] = self._parseTransformation(what[2])
					what[3] = self._parseTransformation(what[3])
				except IndexError:
					pass

			# Add current card
			self._addCard(Card(tag, what, self._comment, extra))
		return False

	#-----------------------------------------------------------------------
	# Parse Geometry
	#-----------------------------------------------------------------------
	def _parseGeometry(self, geobegin):
		voxel = False
		#self._cardEnable = True
		while True:
			line = self._readLine()
			if line is None:
				raise Exception("Error parsing geometry")

			line = self._parseComment(line)
			if line is None: continue

			#line = line.rstrip()

			#if not self._cardEnable or line[0]=='#':
			if line[0]=='#':
				tag, what = self._parseLine(line)
				self._addCard(Card(tag, what, self._comment))

			elif line[0:6] == "GEOEND":
				tag, what = self._parseLine(line)
				self._addCard(Card(tag, what, self._comment))
				return

			elif line[0:6] == "VOXELS":
				tag, what = self._parseLine(line)
				card = Card(tag, what)
				self._addCard(card)
				try:
					card.loadVoxel()
				except:
					say("ERROR: loading voxel file %s.vxl"%(card.sdum()))
					say(sys.exc_info()[1])
				# Read next line
				#line = self._readLine()
				#if line is None:
				#	raise Exception("Error parsing geometry")
				voxel = True

			else:
				break

		# Except the title of the GEOBEGIN
		try: ivopt = int(line[0:5].strip())
		except: ivopt = 0
		try: idbg  = int(line[5:10].strip())
		except: idbg = 0

		# Special treatment for GEOBEGIN
		geobegin.setExtra(line[20:])
		# Save vars for use in parsing
		self.ivopt = geobegin.ivopt = ivopt
		self.idbg  = geobegin.idbg  = idbg

		if idbg<0 and self.geoFormat==FORMAT_SINGLE:
			self.geoFormat = FORMAT_DOUBLE

		if self.verbose: say("GEOBEGIN: Geometry format:", self.geoFormat)
		try:
			self._parseBodies()
			self._parseRegions()
			self._parseVolume()
		except:
			say(sys.exc_info()[0])

		# Parse Lattices only on external file
#		if external:
#			while True:
#				self._lineNo += 1
#				line = utfReadline(fgeo)
#				if len(line)==0: break
#				line = self._parseComment(line)
#				if line is None: continue
#
#				tag, what = self._parseLine(line)
#				if tag == "LATTICE":
#					self._parseLattice(what)
#					self._addCard(Card(tag, what, self._comment))
#				else:
#					say("ERROR: Unexpected card %s found in geometry"%(tag))
#					say(line)

		# Correct last body name to VOXEL
		if voxel and self.geoFormat != FORMAT_FREE:
			# Count all bodies
			nbodies = 1	# +1
			for body in BODY_TAGS:
				try:
					nbodies += len(self.cards[body])
				except KeyError:
					pass
			self.changeName("bn", _body2name(str(nbodies)), "VOXEL")

	#-----------------------------------------------------------------------
	# Parse bodies
	#-----------------------------------------------------------------------
	def _parseBodies(self):
		if self.verbose: say("   Parsing Bodies")
		while True:
			line = self._readLine()
			if line is None: return
			line = self._parseComment(line, False)
			if line is None: continue

			if line[0] == "#":
				(tag, what) = self._parseLine(line)
				self._addCard(Card(tag, what, self._comment))
				continue

			elif line[0] == "$":
				(tag, what) = self._parseLine(line)
				self._addCard(Card(tag, what, self._comment))
				continue

			what = []
			(tag, name) = self._parseBodyLine(line, what)
			name = _body2name(name)

			if tag=="END" or line=="END":
				self._addEND()
				return

			info = CardInfo.get(tag)
			if info is None:
				say("ERROR: line %d: Unknown body type %s for %s" \
					% (self._lineNo[-1], tag, name))
				say(line)
				return
			needwhat = info.nwhats-1
			#self._whatNeeded(what, needwhat)
			# Read and parse needed whats
			while len(what)<needwhat:
				line = self._readLine()
				if line is None:
					raise Exception("Error parsing geometry")
				line = self._parseComment(line, False)
				if line is None: continue
				self._parseBodyLine(line, what)
			# delete extra whats if any
			if len(what)>needwhat: del what[needwhat:]
			what.insert(0, name)
			self._addCard(Card(tag, what, self._comment))

	#-----------------------------------------------------------------------
	# Parse a body line
	#-----------------------------------------------------------------------
	def _parseBodyLine(self, line, what):
		startWhat = len(what)
		if self.geoFormat == FORMAT_SINGLE:
			body = line[2:5].strip()
			name = line[5:10].strip()
			if startWhat>0 and (len(body)>0 or len(name)>0):
				raise Exception("Invalid body definition at line: " \
					+self._lineNo[-1])
			for i in range(10, 61, 10):
				w = line[i:i+10].strip()
				what.append(_str2num(w))

		elif self.geoFormat == FORMAT_DOUBLE:
			body = line[2:5].strip()
			name = line[5:10].strip()
			if startWhat>0 and (len(body)>0 or len(name)>0):
				raise Exception("Invalid body definition at line: "\
					+self._lineNo[-1])
			for i in range(10, 55, 22):
				w = line[i:i+22].strip()
				what.append(_str2num(w))

		else:
			wlist = self._parseFreeFormat(line)
			body = name = ""
			if startWhat==0:
				body = str(wlist[0])
				if len(wlist)>1: name = wlist[1]
				i = 2
			else:
				i = 0

			for i in range(i, len(wlist)):
				what.append(_str2num(wlist[i]))

		return (body, name)

	#-----------------------------------------------------------------------
	# Parse Regions
	#-----------------------------------------------------------------------
	def _parseRegions(self):
		if self.verbose: say("   Parsing Regions")

		name          = ""
		neighbors     = ""
		expstr        = ""
		comment       = ""
		self._comment = ""
		card          = None	# active card
		_cardEnable   = True
		prop          = []
		regions       = set()	# check for duplicate names in FIXED format

		while True:
			line = self._readLine(False)
			if line is None:
				if name!="":
					self._addRegion(name, neighbors, comment, expstr, _cardEnable)
				return
			line = self._parseComment(line, False)
			if line is None: continue

			if line[:6] == "GEOEND":
				self._addRegion(name, neighbors, comment, expstr, _cardEnable)
				self._addEND()
				say("Error: Parsing geometry END is missing")
				return

			# Preprocessor
			if line[0] == "#":	# Parse defines
				# Add current region if any
				oldcomment = self._comment	# Save the comment
				card = self._addRegion(name, neighbors, comment, expstr, _cardEnable)
				# Possible region continuation
				name      = "&"
				neighbors = 0
				expstr    = ""
				comment   = ""
				self._comment = oldcomment	# Restore it

				# Add current card
				(tag, what) = self._parseLine(line)
				if tag == "#endinclude":
					self._cardEnable = True
					self._closeFile()
				self._addCard(Card(tag, what, self._comment))

			# Free format
			elif self.geoFormat == FORMAT_FREE:
				# Free format
				line = line.strip()

				if line=="END":
					oldcomment = self._comment	# remember comment
					self._addRegion(name, neighbors, comment, expstr, _cardEnable)
					self._comment = oldcomment
					self._addEND()
					return

				# Check for region continuation
#				elif line and string.find("+-|()", line[0]) >= 0:
				elif line and "+-|()".find(line[0]) >= 0:
					if self._comment:
						oldcomment = self._comment	# remember comment
						# add region
						card = self._addRegion(name, neighbors, comment, expstr, _cardEnable)

						name      = "&"			# continuation card
						neighbors = 0
						expstr    = line
						comment   = oldcomment		# restore comment
						self._comment = ""		# just for precaution
					else:
						# append to region expression
						expstr += "\n%s"%(line)

				# New region definition
				else:
					match = _REGIONPAT.match(line)
					if match is None:
						say("Error in region: %s"%(line))
						continue

					oldcomment  = self._comment	# remember comment
					thisProp    = self._prop	# and properties
					self._prop  = prop		# restore old properties
					# New region definition
					card = self._addRegion(name, neighbors, comment, expstr, _cardEnable)

					name      = match.group(1)
					neighbors = match.group(2)
					expstr    = match.group(3)
					comment   = oldcomment	# remember comment
					prop      = thisProp		# and properties
					self._comment = ""		# just for precaution
					_cardEnable = self._cardEnable

			# Single precision fixed format
			else:
				n = line[1:6].strip()	# Name

				# New region definition
				if n!="":
					thisComment = self._comment
					# Save previous region
					self._addRegion(name, neighbors, comment,
							expstr, _cardEnable)
					level = len(self._files)-1	# remember level
					if n=="END":
						self._addEND()
						return

					# New region
					name = _region2name(n)
					if name in regions:
						name = "R%d"%(len(regions))
					regions.add(name)

					neighbors = line[6:10].strip()
					expstr    = ""
					comment   = thisComment
					_cardEnable = self._cardEnable
				elif self._comment:
					if comment == "":
						comment = self._comment
					else:
						comment += "\n"+self._comment
				self._comment = ""

				# get expression
				e = line[10:]

				# Convert to free format expression string
				efree = ""

				if abs(self.idbg) == 100:
					width = 9
				else:
					width = 7
				if len(e)%width != 0: e += " "*(width-len(e)%width)

				for i in range(0, len(e), width):
					orop = e[i:i+2]
					if orop=="OR":
						efree += " |"

					body = e[i+2:i+width].strip()
					if body=="":
						say("ERROR: Bad alignment in region definition: '%s'"%(name))
						raise Exception("Error parsing region %s: " \
							"Invalid bodyname found \"%s\" " \
							"column=%d"%(name, e[i+2:i+width], i+2))
					sign = body[0]
					if sign=="+" or sign=="-":
						efree += " %s%s" % (sign,
							_body2name(body[1:]))
					else:
						efree += " +%s" % (_body2name(body))

				if expstr=="":
					expstr = efree
				else:
					expstr += "\n" + efree	# Continue region

	#-----------------------------------------------------------------------
	def _addEND(self):
		self._addCard(Card("END", [], self._comment))

	#-----------------------------------------------------------------------
	def _addRegion(self, name, neighbors, comment, expstr, _cardEnable):
		if name=="" or (name=="&" and expstr==""): return None
		saveEnable = self._cardEnable
		self._cardEnable = _cardEnable

		try: neighbors = int(neighbors)
		except: neighbors = 5

		what = [ name, neighbors ]
		# XXX should not come where with expstr=""
		if len(expstr)>1 and expstr[0]=="\n": expstr = expstr[1:]
		card = Card("REGION", what, comment, expstr)
		card[_VOLUME] = ""
		self._addCard(card)
		self._cardEnable = saveEnable
		return card

	#-----------------------------------------------------------------------
	def _parseLattice(self, what):
		sdum = what[0]
		if sdum != "":
			usdum = sdum.upper().split('#')
			if usdum[0] == "ROT" or usdum[0] == "RO":
				try: id = int(usdum[-1])
				except: id = 1
				what[0] = "Rot#%3d" % (id)

	#-----------------------------------------------------------------------
	def _parseTransformation(self, what):
		if "#" in what:
			uwhat = what.upper()
			if uwhat[:4] == "ROT#" or uwhat[:3] == "RO#":
				try: idx = int(what.split("#")[-1])
				except: idx = 1
				what = "rot#%03d" % (idx)
		return what

	#-----------------------------------------------------------------------
	# Parse Volume data if any
	#-----------------------------------------------------------------------
	def _parseVolume(self):
		if self.ivopt != 3: return
		# Expect data in the format 7E10.5
		if self.verbose: say("   Parsing Volumes")
		self.preprocess()
		lst = self["REGION"]
		lst.sort(Card.cmpPos)
		regionList = [x for x in lst
			if not x.ignore() and x.sdum()!="&"]

		nlines, rem = divmod(len(regionList), 7)
		if rem > 0: nlines += 1
		nr = 0
		while nlines>0:
			# peek next line
			pos  = self._files[-1].tell()
			line = self._readLine()
			if line is None:
				raise Exception("Error parsing geometry")
			if len(line)==0: continue
			if line[0] in ("*", "G"):
				# go back and exit
				self._files[-1].seek(pos)
				break
			for i in range(0, 61 ,10):
				vol = _str2num(line[i:i+10].strip())
				if nr < len(regionList):
					regionList[nr][_VOLUME] = vol
				nr += 1
			nlines -= 1

	#-----------------------------------------------------------------------
	# Used internally while parsing to add a card
	#-----------------------------------------------------------------------
	def _addCard(self, card):
		self._nwhats = card.nwhats()
		self._owhat  = card._what[:]	# keep a copy of original whats

		# Replace special whats with card whats
#		old = []
		for i,w in enumerate(self._cwhat):
			if len(w)>0:
#				old.append((i,card.numWhat(i)))
				card.setWhat(i,w)

		# add properties if any
		for n,v in self._prop:
			card[n] = v

		self.addCard(card)
		if not self._cardEnable: card.disable()

		# Check for possible problems during loading
#		if card.enable:
#			for i,o in old:
#				try:
#					new = card.numWhat(i)
# FIXME
#					if new!="?" and abs(new-o)>1e-10:
#						card.addUserInvalid(i)
#						# Houston we have a problem
#						say("ERROR: Function evaluation doesn't match the value entered in card")
#						say(card.tag, str(card.whats()))
#						say("\tInput WHAT(%d) = %s"%(i,str(o)))
#						say("\tFunction %s = %s"%(w, str(new)))
#				except:
#					pass

		# Reset comment and special whats
		self._comment = ""
		del self._cwhat[:]
		del self._prop[:]

	#-----------------------------------------------------------------------
	# Read an input file
	#-----------------------------------------------------------------------
	def read(self, filename):
		"""Open and read a fluka input file, return True if error"""
		self.filename = filename
		f = self._openFile(filename, "r")
		if f is None: return True
		if self.verbose: say("Reading input file",filename)
		if self.parse(): return True
		assert len(self._files)==0

		# Process input
		self.scanUnits()
		self.renumber()
		self.setFileTime()
		return False

	#-----------------------------------------------------------------------
	# #include card has been modified take appropriate action
	# @return	True if #include has modified the input
	#		False if nothing happened
	#-----------------------------------------------------------------------
	def include(self, includecard):
		if not includecard.enable: return False

		try:	# Check if file exist
			os.stat(includecard.sdum())
		except OSError:
			# File not exist do nothing
			return False

		# Check if #include contains already imported cards
		fromid = self.cardlist.index(includecard)
		# Finding ending of #include
		toid   = fromid+1
		level  = 1
		while toid<len(self.cardlist):
			card = self.cardlist[toid]
			if card.tag == "#include":
				level += 1
			elif card.tag == "#endinclude":
				level -= 1
				if level==0: break
			toid += 1

		# If range is not empty then return
		if toid-fromid>1: return False

		# Create a temporary input to load the included lines
		inp = Input()
		location = 0	# 0=input, 1=bodies, 2=regions
		for i in range(fromid):
			card = self.cardlist[i]
			if not card.enable: continue
			if card.tag == "FREE":
				inp.format = FORMAT_FREE

			elif card.tag == "FIXED":
				inp.format = FORMAT_SINGLE

			elif card.tag == "GLOBAL":
				w4 = card.intWhat(4)
				if w4==2 or w4==3:
					inp.format= FORMAT_FREE
				if card.intWhat(5)>0:
					inp.geoFormat = FORMAT_FREE

			elif card.tag == "GEOBEGIN":
				if card.sdum()[:8]=="COMBNAME":
					inp.geoFormat = FORMAT_FREE
				elif card.sdum()[:8]=="COMBINAT":
					inp.geoFormat = FORMAT_SINGLE
				location += 1

			elif card.tag == "END":
				location += 1

			elif card.tag == "GEOEND":
				location = 0

		f = inp._openFile(includecard.sdum(),"r")
		if location==0:
			inp.parse()
		elif location==1:
			inp._parseBodies()
		elif location==2:
			inp._parseRegions()
		elif location==3:
			inp._parseVolumes()
		else:
			assert False
		assert len(self._files)==0

		pos = fromid+1
		for card in inp.cardlist:
			card.input = self
			self.addCard(card, pos)
			pos += 1

		self.scanUnits()
		self.renumber(fromid)
		self.setModified()

		return True

	#-----------------------------------------------------------------------
	def writeWithInclude(self, filename, backup=True):
		"""Open and write the input to file (specified either by filename,
		or by file pointer). Optionally use a different geometry
		file and a different format"""

		# WARNING FIXME no testing is done for the correctness of the file
		if self.verbose: say("Writing input file",filename)
		if len(self.cardlist) == 0: return True

		del self._files[:]
		del self._lineNo[:]
#		if isinstance(filename, file):
		if isinstance(filename, io.IOBase):
			self._files.append(filename)
			self._lineNo.append(0)
		else:
			self._openFile(filename,"w",backup=backup)

		self.preprocess()
		# Start with single format as unknown
		self.format    = FORMAT_SINGLE
		self.geoFormat = FORMAT_SINGLE

		# Check for volume data
		self.ivopt = 0
		for region in self["REGION"]:
			vol = region.getProperty(_VOLUME,"")
			if vol != "":
				self.ivopt = 3
				break

		# main writing loop
		geolevel = 0	# 0 before GEOBEGIN
				# 1 bodies and VOXELS
				# 2 REGIONs
				# 3 LATTICEs
				# 4 GEOEND
				# 5 after
		geoExternal = False
		includeoff  = 0	# counter of how many includes to exclude when the #include is disabled
		for card in self.cardlist:
			tag = card.tag

			# if we are under an include off do not write the file
			if includeoff > 0:
				if tag=="#include":
					includeoff += 1
				elif tag=="#endinclude":
					includeoff -= 1

			elif tag=="#include":
				# Write the card even if disabled
				self.writeCard(self._files[-1], card, self.format)
				if card.notIgnore():
					self._openFile(card.what(0),"w",backup=backup)
				else:
					includeoff = 1

			elif tag=="#endinclude":
				c = card.commentStr()
				if c: utfWrite(self._files[-1], c)
				if card.notIgnore(): self._closeFile()

			elif tag=="GEOBEGIN":
				card.ivopt = self.ivopt
				self.checkFormat(card)
				flugg = card.sdum()=="FLUGG"
				geoline = card.toStr(self.format).splitlines()
				geolevel = 1
				utfWrite(self._files[-1], card.commentStr())
				utfWriteln(self._files[-1], geoline[0])

				if flugg:
					# Do nothing, skip the geometry
					pass
				else:
					if self.geoFile != "" and card.intWhat(3) not in (0,5):
						utfWriteln(self._files[-1], self.geoFile)
						geoExternal = True
					# Write geo out file before the title
					if self.geoOutFile!="" and card.intWhat(4) not in (0,11):
						utfWriteln(self._files[-1], self.geoOutFile)
					if geoExternal:
						self._openFile(self.geoFile,"w",backup=backup)

				# If VOXELS are present
				# FIXME not good. Needs to write the cards as is
				if "VOXELS" in self.cards:
					self.writeCards(self._files[-1], lambda x:x.tag=="VOXELS", self.format, 1)

				# Write title
				if not flugg:
					utfWriteln(self._files[-1], geoline[1])

			elif geolevel>0:
				if geolevel==1 and tag=="REGION":
					# XXX Could be problematic with the #if's
					geolevel = 2

				elif geolevel==2 and (tag!="REGION" and tag[0]!='#'):
					# XXX could be problematic if #include is before
					# then the END is written in the #include
					geolevel = 3

				# XXX insert region volumes? LATTICE or GEOEND
				if tag in ("LATTICE","LATTSNGL"):
                                        #tag=="LATTICE":
					self.writeCard(self._files[-1], card, self.format)
				elif tag == "VOXELS":
					pass
				elif card._geo:
					self.writeCard(self._files[-1], card, self.geoFormat)
				else:
					if tag=="GEOEND": geolevel = 0
					if geoExternal:
						self._closeFile()
					self.writeCard(self._files[-1], card, self.format)

			else:
				self.writeCard(self._files[-1], card, self.format)
				self.checkFormat(card)

		self._closeFile()
		self.setFileTime()
		assert(len(self._files)==0)
		return False

	#-----------------------------------------------------------------------
	# Write the input to a file
	#-----------------------------------------------------------------------
	def write(self, filename, backup=True):
		"""Open and write the input to file (specified either by filename,
		or by file pointer). Optionally use a different geometry
		file and a different format"""

		if self.verbose: say("Writing input file",filename)
		if len(self.cardlist) == 0: return True

		if "#include" in self.cards:
			self.writeWithInclude(filename, backup)
			return

#		if isinstance(filename, file):
		if isinstance(filename, io.IOBase):
			finp = filename
		else:
			finp = self._openFile(filename, "w", backup=backup)

		# Expressions
		exprNotGeo = lambda x: not x._geo

		# Start with single format as unknown
		self.format    = FORMAT_SINGLE
		self.geoFormat = FORMAT_SINGLE

		# 1st. Write all input cards before GEOBEGIN
		#      skipping geometry cards
		try:
			geobegin = self["GEOBEGIN"][0]
			begin_ = self.cardlist.index(geobegin)
		except:
			# add a GEOBEGIN card
			if self.geoFile!="":
				w3 = 20
			else:
				w3 = ""
			if self.geoOutFile!="":
				w4 = 21
			else:
				w4 = ""
			geobegin = Card("GEOBEGIN", ["COMBNAME","","",w3,w4])
			# Make a guess of the position before the first geometry card
			for pos,card in enumerate(self.cardlist):
				if card.isGeo():
					begin_ = pos
					break
			else:
				begin_ = len(self.cardlist)

		# Check for volume data
		self.ivopt = 0
		for region in self["REGION"]:
			vol = region.getProperty(_VOLUME,"")
			if vol != "":
				self.ivopt = 3
				break

		geobegin.ivopt = self.ivopt

		# Do not write any defines!
		# geolevel defines the position in the input file look the writeCards()
		self.writeCards(finp, exprNotGeo, None, 0, 0, begin_)

		# if there is no geometry... then exit
		if begin_ == len(self.cardlist):
			if finp is not filename:
				self._closeFile()

			self.setFileTime()
			return False

		# 2nd. Write all non-geometry cards which are inside GEOBEGIN .. GEOEND
		try:
			# Find first GEOEND
			geoend = self["GEOEND"][0]
			end_ = self.cardlist.index(geoend)
			self.writeCards(finp, exprNotGeo, None, 1, begin_+1, end_)
		except:
			geoend = None
			end_ = len(self.cardlist)

		# 3th. Write the geometry cards
		flugg = geobegin.sdum()=="FLUGG"
		self.checkFormat(geobegin)
		geoline = geobegin.toStr(self.format).splitlines()
		utfWrite(finp, geobegin.commentStr())
		utfWriteln(finp, geoline[0])

		if flugg:
			# Do nothing, skip the geometry
			fgeo = finp
		if self.geoFile != "" and geobegin.intWhat(3) not in (0,5):
			# Should be GEOBEGIN
			# XXX no defines!
			utfWriteln(finp, self.geoFile)
			try:
				fgeo = self._openFile(self.geoFile, "w", backup=backup)
			except:
				say("Cannot open output geometry file \"%s\""%(filename))
				say("Write geometry to input")
				fgeo = finp

		else:	# Inline geometry
			fgeo = finp

		# Write geo out file before the title
		if self.geoOutFile != "":
			utfWriteln(finp, self.geoOutFile)

		if not flugg:
			# Write VOXELS are present
			if "VOXELS" in self.cards:
				self.writeCards(fgeo, lambda x:x.tag=="VOXELS", self.format, 1)

			# Write title
			utfWriteln(fgeo, geoline[1])

			# 4th Write geometry
			self.writeGeometry(fgeo, self.geoFormat)

		if fgeo is not finp: self._closeFile()

		# 6th. Write remaining cards including the GEOEND if any
		if geoend is None:
			geoend = Card("GEOEND")
			self.writeCard(finp, geoend, self.format)
			self.writeCards(finp, exprNotGeo, None, 4, begin_+1)
		else:
			self.writeCards(finp, exprNotGeo, None, 4, end_)

		# Check for STOP at the end
		if self.cardlist[-1].tag != "STOP":
			self.writeCard(finp, Card("STOP"), self.format)

		if finp is not filename: self._closeFile()

		self.setFileTime()
		return False

	# ----------------------------------------------------------------------
	# Find current format
	# FIXME #defines and state!
	# ----------------------------------------------------------------------
	def checkFormat(self, card):
		if not card.enable: return
		tag = card.tag
		if   tag == "FREE":
			self.format = FORMAT_FREE
		elif tag == "FIXED":
			self.format = FORMAT_SINGLE
		elif tag == "GLOBAL":
			w = card.intWhat(4)
			if w==2 or w==3:
				self.format= FORMAT_FREE
			w = card.intWhat(5)
			if w > 0:
				self.geoFormat = FORMAT_FREE
			elif w < 0:
				# XXX XXX Accept only FREE format geometry XXX XXX
				say("WARNING: Changing GLOBAL what(5)=1")
				card.setWhat(5,1)
				self.geoFormat = FORMAT_FREE
		elif tag == "GEOBEGIN":
			# Force free geometry format
			if card.sdum()[:8]!="COMBNAME" and card.sdum()!="FLUGG":
				# XXX XXX Accept only COMBNAME! XXX XXX
				# To avoid stupid user errors
				say("WARNING: Changing GEOBEGIN sdum=%s to COMBNAME"%(card.sdum()))
				card.setSdum("COMBNAME")
			self.geoFormat = FORMAT_FREE
			#if self.geoFormat == FORMAT_SINGLE:
			#	if card.sdum()=="COMBNAME":
			#		self.geoFormat = FORMAT_FREE
			#	elif card.idbg < 0:
			#		self.geoFormat = FORMAT_DOUBLE

	#-----------------------------------------------------------------------
	# Check COMPOUND cards for correct mixture
	#-----------------------------------------------------------------------
	def checkCompound(self, card, all=False):
		# check sign of whats
		s1 = card.sign(1)
		s2 = card.sign(2)
		n  = card.nwhats()
		if all:
			start = 1
		else:
			start = n-6
		for i in range(start, n, 2):
			j = i+1
			if card.numWhat(i)!=0.0:
				if card.sign(i)!=s1 or \
				   card.sign(j)!=s2:
					say("ERROR: Inconsistent compound mixture (Check signs)")
					say("	What(%d) = %s"%(i,str(card.what(i))))
					say("	What(%d) = %s"%(j,str(card.what(j))))
					say("%s"%(str(card)))
			else:
				card.setSign(i, s1)
				card.setSign(j, s2)

	#-----------------------------------------------------------------------
	# Write geometry only
	#      preprocessor cards are written only when cards
	#      are in the correct position
	#-----------------------------------------------------------------------
	def writeGeometry(self, fgeo, geoFmt=None, fmt=None):
		if geoFmt is None: geoFmt = self.geoFormat

		geolevel = 1

		# Use the formatting of the input
		if fmt is None: fmt = self.format

		# 1st. Write bodies
		self.writeCards(fgeo, lambda x:x.tag in BODY_NOVXL_TAGS or x.tag[0]=="$", geoFmt, 1)
		ends = self.cards.get("END",None)
		if ends:
			# write all END until the first active
			for lastEnd,card in enumerate(ends):
				self.writeCard(fgeo, card, geoFmt)
				if card.notIgnore():
					break
		else:
			if geoFmt == FORMAT_FREE:
				fgeo.write("END\n")
			else:
				fgeo.write("  END\n")

		# 2nd. Write regions
		self.writeCards(fgeo, lambda x:x.tag=="REGION", geoFmt, 2)
		if ends and len(ends)>1:
			# write all END after the last one
			for card in ends[lastEnd+1:]:
				self.writeCard(fgeo, card, geoFmt)
		else:
			if geoFmt == FORMAT_FREE:
				fgeo.write("END\n")
			else:
				fgeo.write("  END\n")

		# Write volume data if needed
		try:
			ivopt = self.ivopt
		except AttributeError:
			ivopt = 0

		if ivopt == 3:
			line = ""
			for region in self.cardsSorted("REGION"):
				if region.name() == "&": continue
				vol = region.getProperty(_VOLUME,"")
				if vol=="": vol="1.0"
				line += "%10s"%(bmath.format(vol, 10))
				if len(line)==70:
					writeln(fgeo, line)
					line = ""
			if len(line)>0: writeln(fgeo, line)

		# 3rd. Write lattices
		if "LATTICE" in self.cards or \
		   "LATTSNGL" in self.cards:
			self.writeCards(fgeo,
				        lambda x:x.tag in ("LATTICE","LATTSNGL"),
                                        #=="LATTICE",
				fmt, 3)

	#-----------------------------------------------------------------------
	# write a single card
	#-----------------------------------------------------------------------
	def writeCard(self, f, card, fmt):
		if not card.enable and (not card.info.disableComment or not commentedCards):
			if0 = True
			f.write("#if 0\n")
		else:
			if0 = False

		if card.comment()!="":
			utfWrite(f, card.commentStr())
		es = card.evalWhatStr()
		if es!="": utfWrite(f, es)
		s = card.toStr(fmt)
#		try:
#			f.write("%s\n"%(s))
#		except (UnicodeDecodeError, UnicodeEncodeError):
#			f.write("%s\n"%(s.encode("utf-8")))
		utfWriteln(f, s)

		if if0: f.write("#endif\n")

	#-----------------------------------------------------------------------
	# Write cards of that fulfill a given condition(lambda)
	# with all the preprocessor cards around them
	# FIXME Write also non empty blocks: #if .. #define/#undef .. #endif
	#       Might end up in duplicate blocks when only #define/#undef
	#       exist inside the geometry
	# geolevel
	#   -1 before defines
	#    0 before geobegin
	#    1 bodies and voxels
	#    2 regions
	#    3 lattices
	#    4 after geoend
	# from_ .. to_   (to_ is not included)
	#-----------------------------------------------------------------------
	def writeCards(self, f, condition, format=None, geolevel=0, from_=0, to_=-1):
		nest    =  0
		start   = -1
		level   =  0
		writing = False

		# Force a fixed format
		if format is not None: fmt = format

		# Loop over all cards
		if to_ == -1: to_ = len(self.cardlist)
		for cid in range(from_, to_):
			card = self.cardlist[cid]

			tag  = card.tag
			if tag == "GEOBEGIN":
				level = 1
			elif tag == "GEOEND":
				level = 4

			if format is None:	# find format
				fmt = self.format

			# Preprocessor nesting
			if card.enable and tag[0]=="#":
				if writing: self.writeCard(f, card, fmt)

				if tag in ("#if", "#ifdef", "#ifndef"):
					if start < 0 and nest == 0:
						start = cid
					nest += 1

				elif tag == "#endif":
					nest -= 1
					if nest == 0:
						start  = -1
						writing = False

					elif nest < 0:
						# misplaced #elif or #else
						start = cid
						nest = 0

				elif tag in ("#elif", "#else") and nest == 0:
					# misplaced #elif or #else
					start = cid
					nest = 1

				elif tag in ("#define", "#undef"):
					if level == geolevel:
						if nest == 0:
							# Write cards on top nesting correct level
							self.writeCard(f, card, fmt)
						elif start>=0 and not writing:
							# Write everything up to now
							for c in self.cardlist[start:cid+1]:
								if c.tag[0]=='#' or condition(c):
									self.writeCard(f, c, fmt)
							start   = -1
							writing = True

					elif geolevel==-1:
						# For anything with high nesting write them anyway
						# when defines are not required
						if nest != 0 and start>=0 and not writing:
							# Write everything up to now
							for c in self.cardlist[start:cid+1]:
								if c.tag[0]=='#' or condition(c):
									self.writeCard(f, c, fmt)
							start   = -1
							writing = True

			# Target cards
			elif condition(card):
				if nest == 0:
					self.writeCard(f, card, fmt)
				else:
					if not writing:
						# Write everything up to now
						for c in self.cardlist[start:cid]:
							if c.tag[0]=='#' or condition(c):
								self.writeCard(f, c, fmt)
						start = -1
						writing = True
					self.writeCard(f, card, fmt)

			# Change the format after the current card if needed
			if card.enable:
				if format is None:	# find format
					self.checkFormat(card)
				if card.tag=="END":
					level += 1

		# In the error case where a #if..#endif is not closed
		if start >= 0:
			for c in self.cardlist[start:]:
				if c.tag[0]=='#' or condition(c):
					self.writeCard(f, c, fmt)

	# ----------------------------------------------------------------------
	# Pickle/UnPickle
	# ----------------------------------------------------------------------
	def dump(self, pickler):
		"""dump input to pickler"""

		# Input file information
		version = 1
		pickler.dump(version)
		pickler.dump(self.filename)
		pickler.dump(self.geoFile)
		pickler.dump(self.geoOutFile)

		# Cards
		for card in self.cardlist:
			card.dump(pickler)

	# ----------------------------------------------------------------------
	def load(self, unpickler):
		"""load input from unpickle"""

		# Input file information
		version         = unpickler.load()
		self.filename   = unpickler.load()
		self.geoFile    = unpickler.load()
		self.geoOutFile = unpickler.load()

		# Cards
		while True:
			try:
				card = Card(ERROR)
				card.load(unpickler)
			except EOFError:
				break
			self.addCard(card)

	#-----------------------------------------------------------------------
	# Cards
	#-----------------------------------------------------------------------
	# Return cards of a specific tag sorted with _pos
	# Possibilities on which:
	#	0: enabled and active	(requires a call to preprocess before)
	#	1: enabled
	#	2: all
	#
	# Warning on REGION maybe the "&" you want to get rid off
	#-----------------------------------------------------------------------
	def cardsSorted(self, tag, which=0):
		if tag == "bodies":
			# Special case returns a dictionary...
			lst = {}
			for tag in BODY_TAGS:
				for x in self[tag]:
					if x.notIgnore(): lst[x.name()] = x
			if "VOXELS" in self.cards:
				for card in self.cards["VOXELS"]:
					if card.notIgnore():
						lst["VOXEL"] = card
						break
		else:
			# sort the list for faster processing
			if isinstance(tag,str):
				cards = self[tag]
			else:
				cards = []
				for t in tag:
					cards.extend(self[t])
#			cards.sort(Card.cmpPos)
			cards.sort(key=attrgetter("_pos"))

			#cards.sort(key=attrgetter("_pos"))
			if which == 0:
				lst = [x for x in cards if x.notIgnore()]
			elif which == 1:
				lst = [x for x in cards if x.enable]
			else:
				lst = cards

			# Special treatment
			if tag == "MATERIAL":
				mats = {}
				for card in lst: mats[card.name()] = True
				pos = 0
				for card in _defaultMaterials:
					if card.name() not in mats:
						lst.insert(pos, card)
						pos += 1

				# FIXME Must append at the end ONLY the icru that are assigned
				for card in self.cardsSorted("ASSIGNMA", which):
					if card.name() not in mats:
						mat = _icruMatDict.get(card.name(), None)
						if mat: lst.append(mat)

			# Append or Insert voxel cards
			if tag in ("REGION", "ASSIGNMA", "MATERIAL", "COMPOUND", "CORRFACT"):
				for voxel in self["VOXELS"]:
					if voxel.ignore() or voxel["@voxel"] is None: continue
					if tag == "REGION":
						# add to the end
						lst.extend(voxel["@voxel"].input["REGION"])
					elif tag == "ASSIGNMA":
						# insert to the beginning
						lst[:0] = voxel["@voxel"].input["ASSIGNMA"]
					elif tag == "MATERIAL":
						# insert to the beginning
						lst[:0] = voxel["@voxel"].input["MATERIAL"]
					elif tag == "COMPOUND":
						# insert to the beginning
						lst[:0] = voxel["@voxel"].input["COMPOUND"]
					elif tag == "CORRFACT":
						# insert to the beginning
						lst[:0] = voxel["@voxel"].input["CORRFACT"]
		return lst

	# ----------------------------------------------------------------------
	# Return an extend list with all cards including the ones from the voxel
	# ----------------------------------------------------------------------
	def allcards(self):
		cl = self.cardlist[:]
		for voxel in self["VOXELS"]:
			if voxel.ignore() or voxel["@voxel"] is None: continue
			cl.extend(voxel["@voxel"].input.cardlist)
		return cl

	#-----------------------------------------------------------------------
	# Cache lists when they are requested
	#-----------------------------------------------------------------------
	def cardsCache(self, tag, what=None):
		try:
			return self.cache[tag]
		except KeyError:
			if what is None:
				lst = self.cardsSorted(tag)
			else:
				lst = [x.what(what) \
					for x in self.cardsSorted(tag)]
			self.cache[tag] = lst
			return lst

	#-----------------------------------------------------------------------
	def clearCache(self, tag=None):
		if tag is not None:
			try: del self.cache[tag]
			except KeyError: pass
		else:
			self.cache.clear()

	#-----------------------------------------------------------------------
	# Find best position from the sorting order in the file
	# search which card is closer to one we've asked
	#-----------------------------------------------------------------------
	def bestPosition(self, tag):
		ci  = CardInfo.get(tag)
		if ci.tag[0] == "#": return None
		cards = [x for x in self[tag] if x.notIgnore()]
		if cards:
			# Similar cards exists then add it after the last one
			pos = max([x._pos for x in cards]) + 1
			# but with indent=0
			while pos<len(self.cardlist) and \
			      (self.cardlist[pos].indent()>0 or self.cardlist[pos].tag == "#endif"):
				pos += 1
			return pos

		cur = ci.order
		mindist = 10000
		mintag  = ""
		mincl   = []
		for c in list(CardInfo._db.values()):
			cards = [x for x in self[c.tag] if x.notIgnore()]
			if cards:
				d = cur - c.order
				if abs(d)<abs(mindist):
					mindist = d
					mintag  = c.tag
					mincl   = cards

		if mintag=="" or len(mincl)==0:
			pos = 0
		elif mindist>=0:
			pos = max([x._pos for x in mincl])+1
		else:
			pos = min([x._pos for x in mincl])

		# but with indent=0
		while pos<len(self.cardlist) and \
		      (self.cardlist[pos].indent()>0 or self.cardlist[pos].tag == "#endif"):
			pos += 1

		# skip backwards any START, STOP card
		if pos<len(self.cardlist)-1:
			while pos>0 and self.cardlist[pos].tag in ("START","STOP"):
				pos -=1

		return pos

	#-----------------------------------------------------------------------
	# Add a card to list at position pos
	# WARNING: renumber is False (on contrast with delCard, moveCard)
	#-----------------------------------------------------------------------
	def addCard(self, card, pos=None, renumber=False):
		card.input = self		# Keep input file link

		try:
			taglist = self.cards[card.tag]
		except KeyError:
			taglist = []
			self.cards[card.tag] = taglist

		taglist.append(card)		# Local (tag list)

		prevCard = None
		if pos is None or pos>=len(self.cardlist):
			if self.cardlist:
				prevCard = self.cardlist[-1]	# Find prev card
			self.cardlist.append(card)	# Global cardlist
			card._pos = len(self.cardlist)-1
		else:
			if pos>0: prevCard = self.cardlist[pos-1]
			self.cardlist.insert(pos, card)
			card._pos = pos

		# Check indent level of previous card
		if prevCard and prevCard.enable:
			if card.enable and prevCard.tag in _INDENT_INC:
				card._indent = prevCard._indent + 1
			else:
				card._indent = prevCard._indent
		else:
			card._indent = 0

		if card.enable and card.tag in _INDENT_DEC:
			card._indent = max(0,card._indent-1)

		if renumber:
			p = pos
			if p is None: p = 0
			self.renumber(p)
		self.setModified()
		return card

	#-----------------------------------------------------------------------
	# Delete card by position
	#-----------------------------------------------------------------------
	def delCard(self, pos, renumber=True):
		card = self.cardlist[pos]
		card._pos  = -1
		card.input = None
		tag = card.tag

		# Delete from the main list
		del self.cardlist[pos]
		if renumber: self.renumber(pos)

		# Find tag list
		taglist = self.cards[tag]
		taglist.remove(card)
		if len(taglist)==0: del self.cards[tag]
		self.setModified()

	#-----------------------------------------------------------------------
	# Delete all cards with a specific tag
	#-----------------------------------------------------------------------
	def delTag(self, tag, renumber=True):
		try:
			for card in self.cards[tag]:
				card._pos  = -1
				card.input = None
				self.cardlist.remove(card)
			del self.cards[tag]
		except KeyError:
			return

		if renumber: self.renumber()
		self.setModified()

	#-----------------------------------------------------------------------
	# Delete geometry cards
	#-----------------------------------------------------------------------
	def delGeometryCards(self):
		self.delTag("REGION",  False)
		self.delTag("LATTICE", False)
		self.delTag("LATTSNGL", False)
		for tag in BODY_TAGS:
			self.delTag(tag, False)
		for tag in TRANSFORM_TAGS:
			self.delTag(tag, False)
		self.delTag("END", False)
		self.renumber()
		self.setModified()

	#-----------------------------------------------------------------------
	# Change card tag
	#-----------------------------------------------------------------------
	def changeTag(self, card, newtag):
		# Remove from previous list
		cl = self.cards[card.tag]
		cl.remove(card)

		# change the tag
		card.changeTag(newtag)

		# Add to new list
		try:
			taglist = self.cards[card.tag]
		except:
			taglist = []
			self.cards[card.tag] = taglist
		taglist.append(card)
		self.setModified()

	#-----------------------------------------------------------------------
	# Move card (src) at position (dest)
	#-----------------------------------------------------------------------
	def moveCard(self, src, dest, renumber=True):
		card = self.cardlist[src]
		del self.cardlist[src]
		if dest >= src: dest -= 1
		self.cardlist.insert(dest, card)
		if renumber:
			self.renumber(min(dest, src), max(dest,src)+1)
		self.setModified()

	#-----------------------------------------------------------------------
	# Replace card at position (pos) with card
	#-----------------------------------------------------------------------
	def replaceCard(self, pos, card):
		"""replace card at position (pos) with card"""
		# remember old card
		old = self.cardlist[pos]

		# replace in the cardlist
		self.cardlist[pos] = card

		# remove old from the taglist
		taglist = self.cards[old.tag]
		taglist.remove(old)
		if len(taglist)==0: del self.cards[old.tag]

		# add the new to the taglist
		try:
			taglist = self.cards[card.tag]
			taglist.append(card)
		except KeyError:
			taglist = [card]
			self.cards[card.tag] = taglist
		card._pos  = pos
		card.input = self
		self.setModified()

		return old

	#=======================================================================
	# Transform a body card, using a 4x4 transformation matrix
	#=======================================================================
	def _format(self, num):
		if abs(num)<zero: return 0.0
		return num

	#-----------------------------------------------------------------------
	def transformBody(self, card, matrix):
		tag = card.tag

		if tag=="SPH":
			center = matrix * card.bodyP()
			card.setWhat(1, self._format(center[0]))
			card.setWhat(2, self._format(center[1]))
			card.setWhat(3, self._format(center[2]))

		elif tag in ("XCC", "YCC", "ZCC"):
			self._transformInfCyl(card, matrix)

		elif tag in ("RCC", "TRC"):
			# XXX MM comment: add a check to switch to
			# infinite cylinders if it's the case?
			center = matrix * card.bodyP()
			card.setWhat(1, self._format(center[0]))
			card.setWhat(2, self._format(center[1]))
			card.setWhat(3, self._format(center[2]))

			length = matrix.multNoTranslation(card.bodyZ())
			card.setWhat(4, self._format(length[0]))
			card.setWhat(5, self._format(length[1]))
			card.setWhat(6, self._format(length[2]))

		elif tag in ("XYP", "XZP", "YZP", "PLA"):
			self._transformPlane(card, matrix)

		elif tag in ("RPP","BOX"):
			self._transformBox(card, matrix)

		elif tag=="REC":
			self._transformREC(card, matrix)

		elif tag in ("WED", "RAW"):
			self._transformWED(card, matrix)

		elif tag=="ELL":
			focus1 = matrix * card.bodyP()
			focus2 = matrix * card.bodyP2()
			card.setWhat(1, self._format(focus1[0]))
			card.setWhat(2, self._format(focus1[1]))
			card.setWhat(3, self._format(focus1[2]))
			card.setWhat(4, self._format(focus2[0]))
			card.setWhat(5, self._format(focus2[1]))
			card.setWhat(6, self._format(focus2[2]))

		elif tag in ("XEC", "YEC", "ZEC"):
			self._transformInfEllCyl(card, matrix)

		elif tag=="ARB":
			for i in range(1,9):
				point = matrix * card.bodyPn(i)
				n = 3*i-2
				card.setWhat(n  , self._format(point[0]))
				card.setWhat(n+1, self._format(point[1]))
				card.setWhat(n+2, self._format(point[2]))

		elif tag=="QUA":
			self._transformQUA(card, matrix)

		else:
			say("Invalid card %s for transformation"%(card.tag))
			return True

		self.setModified()
		return False

	# ----------------------------------------------------------------------
	def _transformPlane(self, card, matrix):
		point  = matrix * card.bodyP()
		normal = matrix.multNoTranslation(card.bodyN())
		dire   = normal.direction(zero)
		if dire == "X":
			self.changeTag(card, "YZP")
			card.setWhat(1, self._format(point[0]))

		elif dire == "Y":
			self.changeTag(card, "XZP")
			card.setWhat(1, self._format(point[1]))

		elif dire == "Z":
			self.changeTag(card, "XYP")
			card.setWhat(1, self._format(point[2]))

		else:
			self.changeTag(card, "PLA")
			card.setWhat(1, self._format(normal[0]))
			card.setWhat(2, self._format(normal[1]))
			card.setWhat(3, self._format(normal[2]))
			card.setWhat(4, self._format(point[0]))
			card.setWhat(5, self._format(point[1]))
			card.setWhat(6, self._format(point[2]))

	# ----------------------------------------------------------------------
	def _transformInfCyl(self, card, matrix):
		point = matrix * card.bodyP()
		dir   = matrix.multNoTranslation(card.bodyN())
		# XXX no scaling for the moment
		radius = card.numWhat(3)
		dirCheck = dir.direction(zero)

		if   dirCheck == "X" or dirCheck == "-X":
			self.changeTag(card, "XCC")
			card.setWhat(1, self._format(point[1]))
			card.setWhat(2, self._format(point[2]))

		elif dirCheck == "Y" or dirCheck == "-Y":
			self.changeTag(card, "YCC")
			card.setWhat(1, self._format(point[2]))
			card.setWhat(2, self._format(point[0]))

		elif dirCheck == "Z" or dirCheck == "-Z":
			self.changeTag(card, "ZCC")
			card.setWhat(1, self._format(point[0]))
			card.setWhat(2, self._format(point[1]))

		else:
			# XXX transform to old version
			if _useQUA:
				# convert card to "QUA" and then transform it
				tag = card.tag
				p = card.bodyP()
				self.changeTag(card, "QUA")
				card.setWhat(1, (1.0, 0.0)[tag=="XCC"])
				card.setWhat(2, (1.0, 0.0)[tag=="YCC"])
				card.setWhat(3, (1.0, 0.0)[tag=="ZCC"])
				card.setWhat(4, 0.0)
				card.setWhat(5, 0.0)
				card.setWhat(6, 0.0)
				card.setWhat(7, -2.0*p.x())
				card.setWhat(8, -2.0*p.y())
				card.setWhat(9, -2.0*p.z())
				card.setWhat(10, -radius**2 + p.length2())
				self._transformQUA(card, matrix)
			else:
				point -= dir  * infinite
				dir    = (2.0 * infinite) * dir
				self.changeTag(card, "RCC")
				card.setWhat(1, self._format(point[0]))
				card.setWhat(2, self._format(point[1]))
				card.setWhat(3, self._format(point[2]))
				card.setWhat(4, self._format(dir[0]))
				card.setWhat(5, self._format(dir[1]))
				card.setWhat(6, self._format(dir[2]))
				card.setWhat(7, self._format(radius))

	# ----------------------------------------------------------------------
	def _transformWED(self, card, matrix):
		point = matrix * card.bodyP()
		X = matrix.multNoTranslation(card.bodyX())
		Y = matrix.multNoTranslation(card.bodyY())
		Z = matrix.multNoTranslation(card.bodyZ())
		card.setWhat( 1, self._format(point[0]))
		card.setWhat( 2, self._format(point[1]))
		card.setWhat( 3, self._format(point[2]))
		card.setWhat( 4, self._format(X[0]))
		card.setWhat( 5, self._format(X[1]))
		card.setWhat( 6, self._format(X[2]))
		card.setWhat( 7, self._format(Y[0]))
		card.setWhat( 8, self._format(Y[1]))
		card.setWhat( 9, self._format(Y[2]))
		card.setWhat(10, self._format(Z[0]))
		card.setWhat(11, self._format(Z[1]))
		card.setWhat(12, self._format(Z[2]))

	# ----------------------------------------------------------------------
	def _transformREC(self, card, matrix):
		point = matrix * card.bodyP()
		X = matrix.multNoTranslation(card.bodyX())
		Y = matrix.multNoTranslation(card.bodyY())
		Z = matrix.multNoTranslation(card.bodyZ())
		card.setWhat( 1, self._format(point[0]))
		card.setWhat( 2, self._format(point[1]))
		card.setWhat( 3, self._format(point[2]))
		card.setWhat( 4, self._format(Z[0]))
		card.setWhat( 5, self._format(Z[1]))
		card.setWhat( 6, self._format(Z[2]))
		card.setWhat( 7, self._format(X[0]))
		card.setWhat( 8, self._format(X[1]))
		card.setWhat( 9, self._format(X[2]))
		card.setWhat(10, self._format(Y[0]))
		card.setWhat(11, self._format(Y[1]))
		card.setWhat(12, self._format(Y[2]))

	# ----------------------------------------------------------------------
	def _transformInfEllCyl(self, card, matrix):
		center = matrix * card.bodyP()
		X      = matrix.multNoTranslation(card.bodyX())
		Y      = matrix.multNoTranslation(card.bodyY())
		Z      = matrix.multNoTranslation(card.bodyN())
		dirX   = X.direction(zero)
		dirY   = Y.direction(zero)
		dirZ   = Z.direction(zero)

		if dirX=="N" or dirY=="N" or dirZ=="N":
			if _useQUA:
				# convert card to "QUA" and then transform it
				tag = card.tag
				p   = card.bodyP()
				r   = card.bodyX() + card.bodyY()
				rx2 = r.x()**2
				ry2 = r.x()**2
				rz2 = r.x()**2
				self.changeTag(card, "QUA")
				card.setWhat(1, (1.0/rx2, 0.0)[tag=="XEC"])
				card.setWhat(2, (1.0/ry2, 0.0)[tag=="YEC"])
				card.setWhat(3, (1.0/rz2, 0.0)[tag=="ZEC"])
				card.setWhat(4, 0.0)
				card.setWhat(5, 0.0)
				card.setWhat(6, 0.0)
				card.setWhat(7, -2.0*p.x()/rx2)
				card.setWhat(8, -2.0*p.y()/ry2)
				card.setWhat(9, -2.0*p.z()/rz2)
				card.setWhat(10, -1.0 + p.x()**2/rx2 \
					+ p.y()**2/ry2 + p.z()**2/rz2)
				self._transformQUA(card, matrix)
			else:
				# Transforming in REC
				self.changeTag(card, "REC")

				center -= Z   * infinite
				Z       = 2.0 * infinite * Z

				card.setWhat( 1, self._format(center[0]))
				card.setWhat( 2, self._format(center[1]))
				card.setWhat( 3, self._format(center[2]))
				card.setWhat( 4, self._format(Z[0]))
				card.setWhat( 5, self._format(Z[1]))
				card.setWhat( 6, self._format(Z[2]))
				card.setWhat( 7, self._format(X[0]))
				card.setWhat( 8, self._format(X[1]))
				card.setWhat( 9, self._format(X[2]))
				card.setWhat(10, self._format(Y[0]))
				card.setWhat(11, self._format(Y[1]))
				card.setWhat(12, self._format(Y[2]))

		elif dirZ=="X" or dirZ=="-X":
			self.changeTag(card, "XEC")
			card.setWhat(1, self._format(center[1]))
			card.setWhat(2, self._format(center[2]))
			XY = X+Y
			card.setWhat(3, self._format(abs(XY.y())))
			card.setWhat(4, self._format(abs(XY.z())))

		elif dirZ=="Y" or dirZ=="-Y":
			self.changeTag(card, "YEC")
			card.setWhat(1, self._format(center[2]))
			card.setWhat(2, self._format(center[0]))
			XY = X+Y
			card.setWhat(3, self._format(abs(XY.z())))
			card.setWhat(4, self._format(abs(XY.x())))

		else:
			self.changeTag(card, "ZEC")
			card.setWhat(1, self._format(center[0]))
			card.setWhat(2, self._format(center[1]))
			XY = X+Y
			card.setWhat(3, self._format(abs(XY.x())))
			card.setWhat(4, self._format(abs(XY.y())))

	# ----------------------------------------------------------------------
	def _transformBox(self, card, matrix):
		point = matrix * card.bodyP()
		X = matrix.multNoTranslation(card.bodyX())
		Y = matrix.multNoTranslation(card.bodyY())
		Z = matrix.multNoTranslation(card.bodyZ())

		XD = X.direction(zero)
		YD = Y.direction(zero)
		ZD = Z.direction(zero)

		if XD!="N" and YD!="N" and ZD!="N":
			self.changeTag(card, "RPP")
			diagonal = X+Y+Z
			card.setWhat(1, self._format(min(point[0], point[0]+diagonal[0])))
			card.setWhat(2, self._format(max(point[0], point[0]+diagonal[0])))

			card.setWhat(3, self._format(min(point[1], point[1]+diagonal[1])))
			card.setWhat(4, self._format(max(point[1], point[1]+diagonal[1])))

			card.setWhat(5, self._format(min(point[2], point[2]+diagonal[2])))
			card.setWhat(6, self._format(max(point[2], point[2]+diagonal[2])))
		elif _useBOX or card.tag == "BOX":
			self.changeTag(card, "BOX")
			card.setWhat( 1, self._format(point[0]))
			card.setWhat( 2, self._format(point[1]))
			card.setWhat( 3, self._format(point[2]))
			card.setWhat( 4, self._format(X[0]))
			card.setWhat( 5, self._format(X[1]))
			card.setWhat( 6, self._format(X[2]))
			card.setWhat( 7, self._format(Y[0]))
			card.setWhat( 8, self._format(Y[1]))
			card.setWhat( 9, self._format(Y[2]))
			card.setWhat(10, self._format(Z[0]))
			card.setWhat(11, self._format(Z[1]))
			card.setWhat(12, self._format(Z[2]))

	# ----------------------------------------------------------------------
	def _transformQUA(self, card, matrix):
		A = bmath.Matrix(4)
		A[0][0] = card.numWhat(1)		# Cxx
		A[1][1] = card.numWhat(2)		# Cyy
		A[2][2] = card.numWhat(3)		# Czz

		A[0][1] = A[1][0] = card.numWhat(4)/2.0	# Cxy
		A[0][2] = A[2][0] = card.numWhat(5)/2.0	# Cxz
		A[1][2] = A[2][1] = card.numWhat(6)/2.0	# Cyz

		A[0][3] = A[3][0] = card.numWhat(7)/2.0	# Cx
		A[1][3] = A[3][1] = card.numWhat(8)/2.0	# Cy
		A[2][3] = A[3][2] = card.numWhat(9)/2.0	# Cz

		A[3][3] = card.numWhat(10)		# C

		minv = matrix.clone()
		minv.inv()

		# New = matrix^T * A * matrix
		AN = minv.T() * (A * minv)

		# XXX optimize (convert to other objects if needed)
		card.setWhat(1, self._format(AN[0][0]))
		card.setWhat(2, self._format(AN[1][1]))
		card.setWhat(3, self._format(AN[2][2]))

		card.setWhat(4, self._format(AN[0][1]*2.0))
		card.setWhat(5, self._format(AN[0][2]*2.0))
		card.setWhat(6, self._format(AN[1][2]*2.0))

		card.setWhat(7, self._format(AN[0][3]*2.0))
		card.setWhat(8, self._format(AN[1][3]*2.0))
		card.setWhat(9, self._format(AN[2][3]*2.0))

		card.setWhat(10, self._format(AN[3][3]))

	#=======================================================================
	# Utilities
	#=======================================================================
	#-----------------------------------------------------------------------
	# Scan for units
	#-----------------------------------------------------------------------
	def scanUnits(self):
		self.units.scan(self)

	#-----------------------------------------------------------------------
	# Convert input to names and check for obsolete and/or non-valid cards
	#-----------------------------------------------------------------------
	def convert2Names(self):
		# Find materials for duplicate checking
		self.clearCache()

		# Check geometry for FLUGG
		try:
			geobegin = self.cards["GEOBEGIN"][0]
			notflugg = geobegin.sdum()!="FLUGG"
		except:
			notflugg = True

		matDict = {}
		for i,n in enumerate(self.materialList()):
			if n not in matDict:
				matDict[n] = [i+1]
			else:
				matDict[n].append(i+1)

		for card in self.cardlist:
			if not card.enable: continue

			# Check for obsolete cards
			tag = card.tag
			if tag in ("ACCURACY", "OUTLEVEL", "EGSCUT",
				   "EGSFIX",   "EGSFLUO",  "EGSRAY",
				   "EGS",      "LANDAU"):
				if card.enable:
					say("WARNING: Obsolete card: %s"%(tag))
					#card.disable()
				continue
			elif not _developer and tag in ("EVENTYPE"):
				if card.enable:
					say("WARNING: Obsolete card: %s"%(tag))
					#card.disable()
				continue
			elif tag=="EXTRAWEI":
				say("WARNING: Changing: %s to USERWEIG"%(tag))
				self.changeAllTags(tag, "USERWEIG")
			elif tag=="DUMPTHEM":
				say("WARNING: Changing: %s to USERDUMP"%(tag))
				self.changeAllTags(tag, "USERDUMP")
			#elif tag=="#ifdef":
			#	self.changeAllTags(tag, "#if", False)
			elif tag=="MATERIAL":
				index = matDict.get(card.sdum(),[])
				if len(index)>1:
					say("WARNING: Material %s already exists with index=%s"%(card.sdum(),str(index)))
					say(str(card))

			# Check input
			try: case  = card.case()
			except:
				say("ERROR: Unknown card.")
				say(str(card))
				#card.disable()
				continue

			if card.info is CardInfo.none():
				say("ERROR: Invalid tag name on card.")
				say(str(card))
				#card.appendComment("Old card tag was: %s"%(card.tag))
				#self.changeTag(card, ERROR)
				if not card.tag:
					self.changeTag(card, ERROR)
				card.invalid = ["Invalid tag name"]
				#card.disable()
				say(str(card))
				continue

			invalid = card.validate(case)
			if len(invalid)>0:
				if invalid[0] is None:
					say("WARNING: No matching case found for card")
				else:
					say("WARNING: Invalid what(s)=%s"%(str(invalid)))
				say(str(card))

			# Convert card to names
			if notflugg: card.convert(True)
		self.setModified()

	#-----------------------------------------------------------------------
	# Validate all active cards
	#-----------------------------------------------------------------------
	def validate(self):
		for card in self.cardlist:
			if card.ignore(): continue
			card.validate()

	#-----------------------------------------------------------------------
	# Changing all cards tag from old to new
	#-----------------------------------------------------------------------
	def changeAllTags(self, old, new, show=True):
		if show:
			say("WARNING: Changing card(s) %s to %s"%(old, new))
		try:
			cl = self.cards[old]
			self.cards[new] = cl
			del self.cards[old]
			for card in cl:
				card.changeTag(new)
		except KeyError:
			pass
		self.setModified()

	#-----------------------------------------------------------------------
	# Change name to all cards
	#-----------------------------------------------------------------------
	def changeName(self, type, old, new):
		# For bodies use special treatment
		if type == "bn":
			pat = re.compile(r"(\b%s\b)" % (re.escape(old)))
			for card in self["REGION"]:
				newextra = pat.sub(new, card.extra())
				if newextra != card.extra():
					card.setExtra(newextra)
		else:
			for card in self.cardlist:
				case = card.info.findCase(card)
				lst  = card.info.find(type, case)
				if isinstance(lst, tuple):
					lst = list(lst)
					lst.pop()	# skip step
				for w in lst:
					if card.what(w) == old:
						card.setWhat(w, new)
		self.setModified()

	#-----------------------------------------------------------------------
	# Put the correct index to cards in the range given by fromPos:toPos
	#-----------------------------------------------------------------------
	def renumber(self, fromPos=0, toPos=-1):
		cardlist = self.cardlist
		if toPos<0 or toPos>len(cardlist):
			toPos = len(cardlist)
		if fromPos == 0:
			indent = 0
		else:
			try:
				card = cardlist[fromPos-1]
				indent = cardlist[fromPos-1]._indent
				if card.enable and card.tag in _INDENT_INC:
					indent += 1
			except IndexError:
				return

		for i in range(fromPos, toPos):
			card = cardlist[i]
			card._pos = i
			tag = card.tag
			if card.enable and tag in _INDENT_DEC: indent = max(0,indent-1)
			card._indent = indent
			if card.enable and tag in _INDENT_INC: indent += 1

		self.setModified()
		self.checkNumbering()		# FIXME XXX (delete this)

	#-----------------------------------------------------------------------
	# Verify the internal order of the cards
	#-----------------------------------------------------------------------
	def checkNumbering(self):
		err = 0
		indent = 0
		for i, card in enumerate(self.cardlist):
			if card._pos != i:
				err += 1
				say("ERROR: Cards out of order pos=#%d it should be #%d"%(card.pos()+1, i))
				say(str(card))
			if card.enable and card.tag in _INDENT_DEC:
				indent = max(0,indent-1)
			if card.enable and card.indent() != indent:
				say("ERROR: Indent out of order pos=#%d card._indent=%d it should be %d"%\
					(i,card._indent,indent))
				say(str(card))
			if card.enable and card.tag in _INDENT_INC:
				indent += 1
		if err>0:
			say("ERROR: %d cards found out of order. Correcting error"%(err))
			say("*** Please contact %s ***"%(__email__))
			import traceback
			traceback.print_stack()
			self.renumber()

	#-----------------------------------------------------------------------
	# Pre-process all cards and set correspondingly the card.active flag
	# nest[] contains the evaluation of the nesting
	#	-1: inactive - do not include any subsequent #if..
	#	 0: false & active - include substitute #if..
	#	 1: true  & active
	#        2: false & active after else
	#        3: true  & active after else
	# @return list of error cards
	# FIXME values do not WORK!
	#-----------------------------------------------------------------------
	def preprocess(self, activeDefines=None):
		self.clearCache()
		define = {}
		errors = []
		self.localDict.clear()
		#if activeDefines exists and has something
		defs = {}
		if activeDefines:
			useInputDefines = False
			for var,val in activeDefines:
				define[var] = 1
				if isinstance(val,str) and val!="" and val[0]=="=":
					try: val = float(eval(val[1:], _globalDict, self.localDict))
					except: pass
				try: self.localDict[var] = float(val)
				except: self.localDict[var] = val
				defs[var] = val
		else:
			# if none or empty list e.g. default run = []
			useInputDefines = True

		nest = [1]
		active13 = True
		for i,card in enumerate(self.cardlist):
			tag    = card.tag
			var    = card.what(0)
			active = nest[-1]

			if not card.enable:
				card.setActive(active13)
				continue

			active13 = active in (1,3)
			card.setActive(active13)

			if   active13 and tag == "#define":
				if useInputDefines:
					define[var] = 1
				else:
					if self.localDict.get(var,"") == "":
						self.localDict[var] = card.numWhat(1)
					elif var not in defs:
						err = "Duplicate #define with same name is not permitted"
						errors.append((card,err))
						card.invalid.append(err)
						say(str(card))
						say(err)

			elif active13 and tag == "#undef" and useInputDefines:
				define[var] = 0

			elif tag == "#if" or tag == "#ifdef":
				if active13:
					nest.append(define.get(var, 0))
				else:
					nest.append(-1)	# ignore this nesting

			elif tag == "#ifndef":
				if active13:
					nest.append(1-define.get(var, 0))
				else:
					nest.append(-1)	# ignore this nesting

			elif tag == "#elif":
				if len(nest)<=1:
					err = "Error unexpected #elif, card: %d"%(i+1)
					errors.append((card,err))
					card.invalid.append(err)
					say(err)
				elif active==1:
					nest[-1] = -1	# terminate this nesting
				elif active==0:
					nest[-1] = define.get(var, 0)
				elif active>1:
					err = "Misplaced #elif after #else, card: %d"%(i+1)
					errors.append((card,err))
					card.invalid.append(err)
					say(err)
				active = nest[-1]
				active13 = active in (1,3)
				card.setActive(active13)

			elif tag == "#else":
				if len(nest)<=1:
					err = "Error unexpected #else, card: %d"%(i+1)
					errors.append((card,err))
					card.invalid.append(err)
					say(err)
				elif active>1:
					err = "Misplaced #else, card: %d"%(i+1)
					errors.append((card,err))
					card.invalid.append(err)
					say(err)
				elif active>=0:
					nest[-1] = 3-nest[-1]

				active = nest[-1]
				active13 = active in (1,3)
				card.setActive(active13)

			elif tag == "#endif":
				if len(nest)<=1:
					err = "Error unexpected #endif, card: %d"%(i+1)
					errors.append((card,err))
					card.invalid.append(err)
					say(err)
				else:
					nest.pop()
				active = nest[-1]
				active13 = active in (1,3)
				card.setActive(active13)

		return errors

	#-----------------------------------------------------------------------
	# Return the body properties as a dictionary of the active/enabled
	# bodies
	#-----------------------------------------------------------------------
	def bodyProperties(self):
		translat   = None
		transform  = None
		transform2 = None
		expansion  = None

		for card in self.cardlist:
			if not card.isGeo(): continue
			if card.ignore(): continue
			del card["@dx"]
			del card["@dy"]
			del card["@dz"]
			del card["@expansion"]
			del card["@transform"]
			del card["@transform2"]

			if len(card.tag)==3:
				if card.tag=="END": break
				if translat:
					card["@dx"] = translat.numWhat(1)
					card["@dy"] = translat.numWhat(2)
					card["@dz"] = translat.numWhat(3)
				if expansion:
					card["@expansion"] = expansion.numWhat(1)
				if transform:
					card["@transform"] = transform.what(1)
				if transform2:
					card["@transform2"] = transform2.what(1)

			elif card.tag[0]=="$":
				if card.tag=="$start_translat":
					if translat is not None:
						say("ERROR: Double %s found #%d"%(card.tag,card.pos()+1))
						say(str(card))
					else:
						translat = card
#
				elif card.tag=="$start_transform":
					if translat is not None:
						say("ERROR: $start_translat must be the toplevel transformation #%d"%\
							(translat.pos()+1))
					if transform is not None and transform2 is not None:
						say("ERROR: Triple %s found #%d"%(card.tag,card.pos()+1))
						say(str(card))
					elif transform is None:
						transform = card
					else:
						transform2 = card
#				elif card.tag=="$start_transform":
#					if transform is not None:
#						say("ERROR: Double %s found #%d"%(card.tag,card.pos()+1))
#						say(str(card))
#					else:
#						transform = card
				elif card.tag=="$start_expansion":
					if expansion is not None:
						say("ERROR: Double %s found #%d"%(card.tag,card.pos()+1))
						say(str(card))
					else:
						expansion = card

				elif card.tag=="$end_translat":
					if translat is None:
						say("ERROR: %s without $start_xxx found #%d"%(card.tag,card.pos()+1))
						say(str(card))
					translat = None
				elif card.tag=="$end_transform":
					if transform is None:
						say("ERROR: %s without $start_xxx found #%d"%(card.tag,card.pos()+1))
						say(str(card))
					if transform2 is not None:
						transform2 = None
					else:
						transform = None
#					transform = None
				elif card.tag=="$end_expansion":
					if expansion is None:
						say("ERROR: %s without $start_xxx found #%d"%(card.tag,card.pos()+1))
						say(str(card))
					expansion = None

	#-----------------------------------------------------------------------
	# Return the region properties as a dictionary of the active/enabled
	# regions
	#-----------------------------------------------------------------------
	def regionProperties(self):
		regionDict = {}
		regionList = []
		matList    = self.materialList(icru=True)
		matDict    = {}
		rotDefi    = {}

		# Find Materials
		for card in _defaultMaterials:	matDict[card.sdum()] = card
		for card in _icruMaterials:	matDict[card.sdum()] = card

		# finally materials in voxels and input
		for card in self.cardsSorted("MATERIAL"):
			try:
				card["@n"] = matList.index(card.sdum())+1
			except ValueError:
				card["@n"] = 1
			matDict[card.sdum()] = card

		# Find regions
		n = 0
		blckhole = matDict["BLCKHOLE"]
		vacuum   = matDict["VACUUM"]

		# Add a new region to dictionaries
		def addRegion(region,name):
			# Assign default values to region
			region["@n"]         = n
			region["@type"]      = REGION_BLACKHOLE
			region["@material"]  = blckhole
			region["@biasingA"]  = None	# All particles biasing
			region["@biasingH"]  = None	# Hadrons
			region["@biasingE"]  = None	# emf
			region["@biasingN"]  = None	# low-neutrons
			region["@assignmat"] = None	# Card that assigned the material
			region["@cont"]      = False	# has continuation cards

			regionDict[name] = region
			regionList.append(region)

		# Input and voxel regions
		for region in self.cardsSorted("REGION"):
			if region.name() == "&": continue
			n += 1
			addRegion(region, region.sdum())

		# add last region
		if regionList: regionDict["@LASTREG"] = regionList[-1]

		# Find transformations
		for card in self.cardsSorted("ROT-DEFI"):
			if card.ignore(): continue
			w1 = card.intWhat(1)
			if w1 < 1000:
				(j, i) = divmod(w1, 100)
			else:
				(i, j) = divmod(w1, 1000)

			rotDefi[i] = card.sdum()
			rotDefi[card.sdum()] = i

		# Voxel assignmat if any
		hasVoxels = False
		for voxel in self["VOXELS"]:
			if voxel.ignore() or voxel["@voxel"] is None: continue
			hasVoxels = True
			for card in voxel["@voxel"].input["ASSIGNMA"]:
				material = matDict.get(card.what(1))
				if material is None:
					say("ERROR: Material %s is not defined #%d" \
						% (card.evalWhat(1),card.pos()+1))
					say(str(card))
					material = _defaultMaterials[0]	# Blackhole
				region   = regionDict.get(card.what(2))
				region["@type"]            = REGION_VOXEL
				region["@assignmat"]       = card
				region["@material"]        = material
				region["@materialDefined"] = True
			break

		# Check all assignmat
		for card in self.cardsSorted("ASSIGNMA"):
			if card.ignore(): continue
			material = matDict.get(card.evalWhat(1))
			if material is None:
				say("ERROR: Material %s is not defined #%d" \
					% (card.evalWhat(1), card.pos()+1))
				say(str(card))
				material = _defaultMaterials[0]	# Blackhole

			matDecay = matDict.get(card.what(6), material)

			fromReg = regionDict.get(card.evalWhat(2))
			if fromReg is None:
				pat = _VOXELPAT.match(card.what(2))
				if pat is None:
					say("ERROR: Region %s is not defined #%d"\
						% (card.what(2), card.pos()+1))
					say(str(card))
				continue
			fromReg = fromReg["@n"]

			try:
				#w = card.what(3)
				#if w == "@LASTREG":
				#	toReg = len(regionList)
				#else:
				toReg = regionDict.get(card.evalWhat(3))
				toReg = toReg["@n"]
			except:
				toReg = fromReg

			step = card.intWhat(4)
			if step == 0: step = 1

			for rg in range(fromReg, toReg+1, step):
				region = regionList[rg-1]
				if hasVoxels and (region.name()=="VOXEL" or _VOXELPAT.match(region.name())):
					region["@type"] = REGION_VOXEL
				else:
					region["@type"] = int(material is blckhole)
				if region["@assignmat"] is not None and warnMultMat:
					say("WARNING: Multiple ASSIGNMAT for region %s:"%(region.name()))
					say("  Previous: #%d\n%s"%(region["@assignmat"].pos()+1, region["@assignmat"]))
					say("  Current : #%d\n%s"%(card.pos()+1, str(card)))
				region["@assignmat"]       = card
				region["@material"]        = material
				region["@matDecay"]        = matDecay
				region["@materialDefined"] = True

		# Check Biasing cards
		for card in self.cardsSorted("BIASING"):
			if card.ignore(): continue

			what1 = card.intWhat(1)
			if what1<0: continue

			elif what1 == 1:
				tag = "@biasingH"
			elif what1 == 2:
				tag = "@biasingE"
			elif what1 == 3:
				tag = "@biasingN"
			else:
				tag = "@biasingA"

			fromReg  = regionDict.get(card.what(4))
			if fromReg is None: continue
			fromReg  = fromReg["@n"]
			try:
				toReg = regionDict.get(card.what(5))
				toReg = toReg["@n"]
			except:
				toReg = fromReg
			try: step = int(card.what(6))
			except: step = 1

			for rg in range(fromReg, toReg+1, step):
				region = regionList[rg-1]
				region[tag] = card

		# Check Lattices
		for card in self.cardsSorted(("LATTICE","LATTSNGL")):
			if card.ignore(): continue

			fromReg  = regionDict.get(card.what(1))
			if fromReg is None: continue
			fromReg  = fromReg["@n"]
			if card.tag == "LATTICE":
				try:
					toReg = regionDict.get(card.what(2))
					toReg = toReg["@n"]
				except:
					toReg = fromReg
				try: step = int(card.what(3))
				except: step = 1
				rotdefi  = card.sdum()
				rotdefi2 = ""
			else:
				toReg = fromReg
				rotdefi  = card.what(2)
				rotdefi2 = card.what(3)
#			try:
#				toReg = regionDict.get(card.what(2))
#				toReg = toReg["@n"]
#			except:
#				toReg = fromReg
#			try: step = int(card.what(3))
#			except: step = 1
#
#			rotdefi = card.sdum()

			for rg in range(fromReg, toReg+1, step):
				region = regionList[rg-1]
				region["@type"]            = REGION_LATTICE
				region["@material"]        = vacuum
				region["@materialDefined"] = True
				region["@lattice"]         = card
				region["@rotdefi"]         = rotdefi
				region["@rotdefi2"]        = rotdefi2

		# Check if all regions have material definition
		#for region in regionList:
		#	if region.get("@materialDefined") is None:
		#		say("WARNING: Region %r is not assigned any material"\
		#			%(region.sdum()))
		return regionList

	#-----------------------------------------------------------------------
	# Append to card's properties the additional value from tag's
	# @param tag		cards to append the extra property
	# @param proptag	cards from where to get the property
	# @param prop		property name to add in cards
	# @param valueWhat	which value to add as property
	# @param fromWhat,toWhat,stepWhat	range of cards
	# @param sdum		filter cards with the specific sdum
	# @param func		apply a function to value
	#-----------------------------------------------------------------------
	def addCardProperty(self, tag, proptag, prop, valueWhat,
			fromWhat, toWhat=-1, stepWhat=-1, sdum=None, func=None):

		# create a dictionary of cards for faster accessing
		cardlist = self.cardsSorted(tag)
		if tag=="REGION":
			cardlist = [x for x in cardlist if x.name()!="&"]	# Remove continuations of regions

		dictionary = {}
		for card in cardlist:
			dictionary[card.sdum()] = card
			card[prop] = 0.0

		for card in self.cardsSorted(proptag):
			if sdum is not None and card.sdum() != sdum: continue
			try: from_ = dictionary[card.what(fromWhat)]["@n"]
			except: continue
			if from_ is None: continue
			to_  = from_
			step = 1
			if toWhat>=0:
				try:
					toW = card.what(toWhat)
					if toW[0]=='@':
						to_ = len(cardlist)
					else:
						to_  = dictionary[toW]["@n"]
					step = card.intWhat(stepWhat)
				except:
					pass
			if step<=0: step=1

			value = card.numWhat(valueWhat)
			if func: value = func(value)

			for c in cardlist[from_-1:to_:step]:
				c[prop] = value

	#-----------------------------------------------------------------------
	# Convert material to name/number
	#-----------------------------------------------------------------------
	def material(self, mat, toName):
		lst = self.materialList(0,False,True)
		if toName:
			mat -= 1		# 1 based
			if mat<0 or mat>=len(lst):
				if mat != -1:	# Accept -1 == 0 index
					raise Exception("Invalid material index requested idx=%d"%(mat+1))
				return ""
			return lst[mat]
		else:
			try:
				return lst.index(mat)+1
			except:
				return 0

	#-----------------------------------------------------------------------
	# Return material list
	#-----------------------------------------------------------------------
	def materialList(self, which=0, icru=False, assigned=False):
		lst = [m.sdum() for m in _defaultMaterials]

		# Add first voxel materials
		# FIXME maybe apply it on the card and not on the Input!!!!
		for voxel in self["VOXELS"]:
			if voxel.ignore() or voxel["@voxel"] is None: continue
			lst.extend([m.sdum() for m in voxel["@voxel"].input["MATERIAL"]])

		# Add user defined materials
		for card in self.cardsSorted("MATERIAL", which):
			name  = card.sdum()
			index = card.intWhat(4)-1	# 0 based!
			if index<0:
				try:    index = lst.index(name)
				except: index = len(lst)
			# Append or replace
			if index == len(lst):
				lst.append(name)
			elif index < len(lst):
				if lst[index] != "" and lst[index] != name:
					say("Warning: overriding material index %d=%s by %s" % \
						(index, lst[index], name))
				lst[index] = name
			else:
				lst.extend([""] * (index-len(lst)+1))
				lst[index] = name

		if icru:
			# Add all icru materials regardless of their order
			lst.extend([m.sdum() for m in _icruMaterials])
		elif assigned:
			# Add all materials assigned, especially for the ICRU ones
			# correct index assigned by FLUKA
			# FIXME Not very intelligent.
			for card in self.cardsSorted("ASSIGNMA", which):
				mat = card.what(1)
				if mat not in lst:
					lst.append(mat)

		return lst

	#-----------------------------------------------------------------------
	# Return material dictionary
	#-----------------------------------------------------------------------
	def materialDict(self):
		matDict = {}

		# Default materials
		for card in _defaultMaterials:
			matDict[card.sdum()] = card

		# ICRU materials
		for card in _icruMaterials:
			matDict[card.sdum()] = card

		# Voxel materials
		for voxel in self["VOXELS"]:
			if voxel.ignore() or voxel["@voxel"] is None: continue
			for card in voxel["@voxel"].input["MATERIAL"]:
				matDict[card.sdum()] = card

		# ... and finally from the input
		for card in self.cardsSorted("MATERIAL"):
			matDict[card.sdum()] = card

		return matDict

	#-----------------------------------------------------------------------
	# material ZAID description
	#-----------------------------------------------------------------------
	def materialZAID(self, material):
		# Check if it is a compound or a single material

		# Single isotope?
		Z = material.intWhat(1)
		if Z>0:
			A = material.intWhat(6)
			return [(Z*1000+A, 1.0)]

		rho = material.numWhat(3)		# Density

		# Compound search
		ZAID = []
		for compound in self["COMPOUND"]:
			if compound.ignore(): continue
			if compound.sdum() != material.sdum(): continue

			for i in range(1,len(compound.whats()),2):
				frac = compound.numWhat(i)
				mat  = compound.what(i+1)
				if mat == "": continue
				if mat[0]=="-":
					negative = True
					mat = mat[1:]
				else:
					negative = False

				# Find material
				for card in self["MATERIAL"]:
					if card.ignore(): continue
					if card.sdum() == mat: break
				else:
					# Search default materials
					for card in _defaultMaterials:
						if card.sdum() == mat: break

				# Scan recursively
				for za,f in self.materialZAID(card):
					# Scale fractions
					if negative:
						# Density
						f = abs(f)*frac*rho
					else:
						# Mass or atom fraction
						f = abs(f)*frac
					ZAID.append((za,f))
		return ZAID

	#-----------------------------------------------------------------------
	# Convert region to name/number
	#-----------------------------------------------------------------------
	def region(self, reg, toName):
		# Regions removing continuation cards
		regions = [x for x in self.cardsSorted("REGION") if x.name()!="&"]
		if len(regions) == 0: return None

		if toName:
			if reg<=0 or reg>len(regions):
				return ""
			return regions[reg-1].sdum()
		else:
			try:
				return [x.what(0) for x in regions].index(reg)+1
			except:
				return 0

	#-----------------------------------------------------------------------
	# @param idx		rotation to return
	#			Can be prefixed with "-"
	#			None to return a dictionary
	# @param rotdefi	cached dictionary
	# @param inv		return the inverse of the idx matrix.
	#			Accepted even if idx is prefixed with "-"
	# @return rotation matrix from a rot-defi cards or a dictionary with all
	#         transformations
	# FIXME I should combine the rotdefi with the self.cache!
	#-----------------------------------------------------------------------
	def getTransformation(self, idx=None, rotdefi=None, inv=False):
		neg = False
		if isinstance(idx,str):
			if len(idx)>0 and idx[0]=='-':
				neg = True
				idx = idx[1:]
			match = _ROTPAT.match(idx)
			if match:
				idx = int(match.group(1))
		elif isinstance(idx,int):
			if idx<0:
				neg = True
				idx = -idx
		else:
			# Return dictionary of transformations
			rotdefi = {}

		if inv: neg = not neg

		# Check if is cached
		if idx is not None:
			if rotdefi is not None:
				matrix = rotdefi.get(idx)
				if isinstance(matrix,int): matrix = rotdefi.get(matrix)
				if matrix is None:
					return bmath.Matrix(4, type=1)
				if neg:
					matrix = matrix.clone()
					matrix.inv()
				return matrix
			else:
				matrix = bmath.Matrix(4, type=1)

		last = 0
		# Scan transformations
		for card in self.cardsSorted("ROT-DEFI"):
			if card.ignore(): continue

			w1 = card.intWhat(1)
			if w1 < 1000:
				(j, i) = divmod(w1, 100)
			else:
				(i, j) = divmod(w1, 1000)

			if idx is None:
				if i==0: # try with the name
					i = rotdefi.get(card.sdum())
					matrix = rotdefi.get(i)
				else:
					matrix = rotdefi.get(i)
					if matrix is not None and rotdefi.get(card.sdum())!=i:
						# FIXME store the error somewhere...
						say("WARNING: ROT-DEFI Index and name do not match #%d"%(card.pos()+1))
						say(str(card))

				if matrix is None:
					matrix = bmath.Matrix(4, type=1)
					if i is None:
						last += 1
						i = last
					else:
						last = max(last,i)
					rotdefi[i] = None

					if card.sdum()!="":rotdefi[card.sdum()] = i
					rotdefi["rot#%03d"%(i)] = i

			else:
				if i!=idx and card.sdum()!=idx: continue

			if j==0: j=3
			theta = card.numWhat(2)
			phi   = card.numWhat(3)
			x = card.numWhat(4)
			y = card.numWhat(5)
			z = card.numWhat(6)

			if x!=0.0 or y!=0.0 or z!=0.0:
				matrix = bmath.Matrix.translate(x,y,z)*matrix

			if phi!=0.0:
				m = bmath.Matrix(4,type=1)
				m.rotate(math.radians(-phi), j-1)
				matrix = m*matrix

			if theta!=0.0:
				m = bmath.Matrix(4,type=1)
				m.rotate(math.radians(-theta), (j+1)%3)
				matrix = m*matrix

			if idx is None:
				rotdefi[i] = matrix

		if neg:
			matrix.inv()

		if idx is None:
			return rotdefi
		return matrix

	# ----------------------------------------------------------------------
	# Rot-Defi list
	# ----------------------------------------------------------------------
	def rotdefiList(self, negative=False):
		rot = {}
		try:
			for card in self.cards["ROT-DEFI"]:
				w1 = card.intWhat(1)
				if w1 < 1000:
					(j, i) = divmod(w1, 100)
				else:
					(i, j) = divmod(w1, 1000)
				if i>0:
					rot["rot#%03d"%(i)] = True
				if len(card.sdum()) > 0:
					rot[card.sdum()] = True
		except:
			return [""]

		l = list(rot.keys())
		l.sort()
		if negative:
			lst = [""]
			for i in l:
				lst.append(i)
				lst.append("-%s"%(i))
		else:
			lst = l
			lst.insert(0,"")

		del rot
		return lst

	# ----------------------------------------------------------------------
	# @return a list of linked ROT-DEFI cards with the same idx
	# ----------------------------------------------------------------------
	def rotdefiCards(self, idx):
		return [x for x in self.cardsSorted("ROT-DEFI") if x.sdum()==idx]

	# ----------------------------------------------------------------------
	# Refresh cards information
	# 1. ?? (renumber??)
	# 2. Reload voxels
	# 3. Properties?
	# ----------------------------------------------------------------------
	def refresh(self):
		for voxel in self["VOXELS"]:
			voxel.loadVoxel()

#===============================================================================
# Internal Utilities
#===============================================================================
#-------------------------------------------------------------------------------
def _str2num(w):
	"""Check if argument w is a valid fortran number and
	return argument converted to a valid python number"""
	try:
		if len(w)>1:
			p = string.index("+-.0123456789",   w[0])
			p = string.index(".0123456789dDeE", w[1])
			p = string.index(".0123456789",     w[-1])
			we = w.replace( "d", "E")
			we = we.replace("D", "E")
			we = we.replace(" ", "")
			wn = float(we)
			wi = int(wn)
			if abs(wi)<10000000000 and float(wi) == wn: return wi
			return we
		return int(w)
	except:
		return w

#-------------------------------------------------------------------------------
def _intWhat(what, n):
	try: return int(float(_str2num(what[n])))
	except: return 0

#-------------------------------------------------------------------------------
def _body2name(name):
	try:
		return BODY_PREFIX + str(int(name))
	except:
		return str(name)

#-------------------------------------------------------------------------------
def _region2name(name):
	try:
		return REGION_PREFIX + str(int(name))
	except:
		pass
	name = str(name)
	if _NAMEPAT.match(name): return name
	return REGION_PREFIX + name

#-------------------------------------------------------------------------------
# Initialize classes:
#	CardInfo
#	Particle
#-------------------------------------------------------------------------------
def init(filename=None):
	global _defaultMaterials, _defaultMatDict, _icruMaterials, _icruMatDict
	global _neutronGroups, _lowMaterials, _usermvax
	global NGROUPS, BODY_TAGS, BODY_NOVXL_TAGS, FLAIR_TAGS, OBJECT_TAGS, TRANSFORM_TAGS

	if filename is None:
		filename = os.path.join(os.path.dirname(__file__), _database)

	# read card database and prepare CardInfo classes
	cardini = configparser.RawConfigParser()
	cardini.read(filename)

	CardInfo._db[None] = CardInfo(ERROR, [None],
			[["-"]*7], [["?"]*7],
			[], None, "Error in tag...")

	# Go through all cards
	for name in cardini.sections():
		# Ignore sections with lower case letters
		# it does not include the '#xxx' sections
		if name[0].islower(): continue

		# name is complete, tag is only 8 characters max
		if len(name)>8 and name[0] not in ("$","#"):
			tag = name[0:8]
		else:
			tag = name

		meaning = cardini.get(name, "meaning")

		group = cardini.get(name, "group").split()

		# Capitalize each group
		for i in range(len(group)):
			group[i] = group[i].capitalize()

		# Ranges
		range_  = []
		default = []
		extra   = []
		assert_ = []
		for i in range(20):
			n = "range.%d" % (i)
			try: rg = cardini.get(name, n).split()
			except: break
			range_.append(rg)

			n = "defaults.%d" % (i)
			df = cardini.get(name, n).split()
			default.append(df)

			# Get extra info
			ex = []
			for j in range(100):
				n = "extra.%d.%d" % (i, j)
				try:
					e = cardini.get(name, n)
				except:
					# Correct sdum
					if j<7:
						e = ""
					else:
						break
				ex.append(e)
			extra.append(ex)

			# Get assert info
			ass = []
			for j in range(100):
				n = "assert.%d.%d" % (i, j)
				try:
					a = cardini.get(name, n)
				except:
					break
				else:
					ass.append(a)
			assert_.append(ass)

		# Cleanup a bit, check for empty range_, assert_
		empty = True
		for r in range_:
			if len(r)>0:
				empty = False
				break
		if empty: range_=[]

		empty = True
		for a in assert_:
			if len(a)>0:
				empty = False
				break
		if empty: assert_=[]

		# append to CardInfo database
		cinfo = CardInfo(name, group, range_,
				extra, assert_, default, meaning)
		try:
			disable = cardini.get(name, "disable")
			cinfo.setDisableComment(disable=="comment")
		except:
			pass
		CardInfo._db[tag] = cinfo

	# Create bodies cards
	BODY_NOVXL_TAGS = [ x.tag for x in list(CardInfo._db.values())
			if "Geometry" in x.group and len(x.tag)==3 and x.tag!="END" ]
	BODY_TAGS = BODY_NOVXL_TAGS[:]
	BODY_TAGS.append("VOXELS")

	BODY_NOVXL_TAGS.sort()
	BODY_TAGS.sort()

	# Geometry transformation tags
	TRANSFORM_TAGS = [ x.tag for x in list(CardInfo._db.values()) if x.tag[0]=="$" ]
	TRANSFORM_TAGS.sort()

	# Create flair tags
	FLAIR_TAGS = [ x.tag for x in list(CardInfo._db.values()) if "Flair" in x.group ]
	FLAIR_TAGS.sort()
	FLAIR_TAGS.remove("!coffee")

	# Object tags = Flair tags + ROT-DEFI
	OBJECT_TAGS = FLAIR_TAGS[:]
	OBJECT_TAGS.append("ROT-DEFI")
	OBJECT_TAGS.append("BEAM")

	# Particles
	section = "particles"
	for pid, name in cardini.items(section):
		if pid[:3] == "pid":
			tag, id = pid.split(".")
			mass    = cardini.getfloat(section,"mass."+id)
			comment = cardini.get(section,"comment."+id)
			try:
				pdg = int(cardini.get(section,"pdg."+id))
			except:
				pdg = 0
			Particle.add(int(id), name, mass, comment, pdg)
	Particle.makeLists()

	# Read default materials and create fake MATERIAL card's
	section = "materials"
	i = 1
	while True:
		n = "mat.%d" % (i)
		try: name = cardini.get(section, n)
		except: break

		desc  = cardini.get(     section, "desc.%d"  % (i))
		Amass = cardini.getfloat(section, "Amass.%d" % (i))
		Z     = cardini.getfloat(section, "Z.%d"     % (i))
		rho   = cardini.getfloat(section, "rho.%d"   % (i))
		mat   = Card("MATERIAL",[name, Z, Amass, rho, i],desc)
		mat["@n"] = i
		_defaultMaterials.append(mat)
		_defaultMatDict[name] = mat
		i += 1

	# Read special ICRU materials
	i = -1
	while True:
		n = "mat.%d" % (i)
		try: name = cardini.get(section, n)
		except: break

		desc  = cardini.get(     section, "desc.%d"  % (i))
		Amass = cardini.getfloat(section, "Amass.%d" % (i))
		Z     = cardini.getfloat(section, "Z.%d"     % (i))
		rho   = cardini.getfloat(section, "rho.%d"   % (i))
		mat   = Card("MATERIAL",[name, Z, Amass, rho],desc)
		mat["@n"] = i
		_icruMaterials.append(mat)
		_icruMatDict[name] = mat
		i -= 1

	# Low neutrons energy groups
	section = "n-groups"
	# Read number of groups
	try:
		groups = cardini.get(section,"groups")
		NGROUPS = list(map(int, groups.split()))
	except:
		say("ERROR! No neutron energy groups defined")

	pat = re.compile(r"^(\S*) *(.*)$")
	for g in NGROUPS:
		groupsList  = []
		groupsListS = []
		for i in range(1,g+2):
			n = "ene.%d.%d" % (g,i)
			try: energy = cardini.get(section, n)
			except: break
			m = pat.match(energy)
			groupsList.append(float(m.group(1)))
			groupsListS.append(m.group(2))
		_neutronGroups[g]  = groupsList
		_neutronGroupsS[g] = groupsListS

	# Low neutrons energy groups
	for g in NGROUPS:
		section = "low-neut"
		i = 1
		materials = []
		while True:
			n = "elem.%d.%d" % (g,i)
			try: elem = cardini.get(section, n)
			except: break

			mat = LowNeutMaterial()
			mat.elem = elem
			mat.desc = cardini.get(section, "mat.%d.%d"  % (g,i))
			mat.temp = cardini.get(section, "temp.%d.%d" % (g,i))
			mat.db   = cardini.get(section, "db.%d.%d" % (g,i))
			mat.name = cardini.get(section, "name.%d.%d" % (g,i))
			ids      = cardini.get(section, "ids.%d.%d" % (g,i)).split()
			mat.id1  = int(ids[0])
			mat.id2  = int(ids[1])
			mat.id3  = int(ids[2])
			mat.g    = cardini.get(section, "g.%d.%d" % (g,i))
			materials.append(mat)
			i += 1
		_lowMaterials[g] = materials

	# Read user routines
	section = "usermvax"
	for routine, desc in cardini.items(section):
		_usermvax[routine] = desc

	del cardini

#===============================================================================
if __first:
	__first = False
	init()

if __name__ == "__main__":
	inp = Input()
	inp.verbose = True
	inp.read(sys.argv[1])
	inp.convert2Names()
	inp.regionProperties()
#	say(inp.materialList())
	if len(sys.argv)>=3:
		#f = open(sys.argv[3],"w")
		#inp.writeGeometry(f, FORMAT_FREE, FORMAT_FREE)
		#f.close()
		inp.write(sys.argv[2])

#	import bz2
#	f = bz2.BZ2File(sys.argv[2],'w')
#	inp.dump(pickle.Pickler(f))
#	f.close()

#	f = bz2.BZ2File(sys.argv[2],'r')
#	inp2 = Input()
#	inp2.load(pickle.Unpickler(f))
#	f.close()
	#say(inp.getTransformation(1))
	#say(inp.getTransformation("rot#001"))
	#print inp.getTransformation("CDQ4R1")
	#inp.write(sys.argv[2])
