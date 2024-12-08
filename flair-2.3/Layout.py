 
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
# Date:	16-Jun-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import sys
import math
import string

import Input
import bmath
import undo
import bFileDialog
from log import say

#===============================================================================
# Groups order
#===============================================================================
_BODY_PARENTHESIS   = "( )"

#===============================================================================
# Layout class for the specific input/canvas display
#===============================================================================
class LayoutManager:
	"""Layout Manager class"""
	def __init__(self, flair):
		self.flair = flair
		self.updateInput()

		# lists that are dynamically linked to an Input File
		self.__regionCache   = None
		self.__materialCache = None
		self.__defineCache   = None
		self.__rotdefiCache  = None
		self.__rotdefiCache2 = None
		self.__binningCache  = None
		self.__bodyCardCache = None

	# ----------------------------------------------------------------------
	# set input file to use
	# ----------------------------------------------------------------------
	def updateInput(self):
		self.input = self.flair.project.input

	# ----------------------------------------------------------------------
	# Enable/disable caching of lists
	# ----------------------------------------------------------------------
	def cache(self, on):
		if on:
			self.__regionCache   = self.regionList()
			self.__materialCache = self.materialList()
			self.__defineCache   = self.defineList()
			self.__rotdefiCache  = self.rotdefiList()
			self.__rotdefiCache2 = self.rotdefiList2()
			self.__binningCache  = self.binningList()
		else:
			self.__bodyCardCache = None
			self.__regionCache   = None
			self.__materialCache = None
			self.__defineCache   = None
			self.__rotdefiCache  = None
			self.__rotdefiCache2 = None
			self.__binningCache  = None

	# ----------------------------------------------------------------------
	# Return body list
	# ----------------------------------------------------------------------
	def bodyList(self, inctag=True):
		# scan all cards belonging to bodies group
		if self.__bodyCardCache is not None:
			return self.__bodyCardCache

		if inctag:
			lst = ["", _BODY_PARENTHESIS]
			for body in Input.BODY_TAGS:
				try: lst.extend(["%s [%s]"%(x.sdum(),x.tag)
					for x in self.input.cards[body]])
				except: pass
		else:
			lst = []
			for body in Input.BODY_TAGS:
				try: lst.extend([x.sdum()
					for x in self.input.cards[body]])
				except: pass
		if "VOXELS" in self.input.cards:
			lst.append("VOXEL")
		lst.sort(key=str.upper)
		return lst

	# ----------------------------------------------------------------------
	# Return region list
	# ----------------------------------------------------------------------
	def regionList(self):
		if self.__regionCache is not None:
			return self.__regionCache

		lst = [""]

		# Add defined regions
		try:
			lst.extend([x.sdum()
				for x in self.input.cards["REGION"]
					if x.sdum()!="&"])
		except: pass

		# Add voxel regions
		for voxel in self.input["VOXELS"]:
			if voxel.ignore(): continue
			if voxel["@voxel"] is None:
				lst.append("VOXEL")
			else:
				lst.extend([x.sdum()
					for x in voxel["@voxel"].input["REGION"]])

		lst.append("@LASTREG")
		return lst

	# ----------------------------------------------------------------------
	def allRegionList(self):
		lst = self.regionList()[:]
		lst.append("@ALLREGS")
		return lst

	# ----------------------------------------------------------------------
	def signRegionList(self):
		lst = self.regionList()[:]
		for i in range(len(lst)-2,0,-1):
			reg = lst[i]
			lst.insert(i,"-%s"%(reg))
		return lst

	# ----------------------------------------------------------------------
	# Return define variable list
	# ----------------------------------------------------------------------
	def defineList(self):
		if self.__defineCache is not None:
			return self.__defineCache

		lst = [""]
		try: lst.extend([x.sdum() for x in self.input.cards["#define"]])
		except: pass
		lst.sort(key=str.upper)
		return lst

	# ----------------------------------------------------------------------
	# Rot-Defi list for LATTICE
	# ----------------------------------------------------------------------
	def rotdefiList(self):
		if self.__rotdefiCache is not None:
			return self.__rotdefiCache
		return self.input.rotdefiList(True)

	# ----------------------------------------------------------------------
	# Rot-Defi list for ROTPRBIN
	# ----------------------------------------------------------------------
	def rotdefiList2(self):
		if self.__rotdefiCache2 is not None:
			return self.__rotdefiCache2

		lst = self.input.rotdefiList(True)
		for i,item in enumerate(lst):
			if item.startswith("rot#"):
				lst[i] = (item,int(item[4:]))
		lst.append(("Ignore",0))
		lst.append(("Ignore",-1))
		return lst

	# ----------------------------------------------------------------------
	# Binning list
	# ----------------------------------------------------------------------
	def binningList(self):
		if self.__binningCache is not None: return self.__binningCache

		lst = [""]
		try: lst.extend([x.sdum() for x in self.input.cards["USRBIN"]])
		except: pass
		try: lst.extend([x.sdum() for x in self.input.cards["EVENTBIN"]])
		except: pass
		lst.extend(list(range(1,len(lst))))
		return lst

		# Try #1
		#cardlist = []
		#try: cardlist.extend(self.input.cards["USRBIN"])
		#except: pass
		#try: cardlist.extend(self.input.cards["EVENTBIN"])
		#except: pass
		#cardlist.sort(Input.Card.cmpPos)
		#lst = [x.sdum() for x in cardlist]
		#lst.append("")
		#lst.sort(key=str.upper)
		#lst.extend([("%d %s"%(i+1,c.sdum()),i+1)
		#		for i,c in enumerate(cardlist)])
		#return lst

		# Try #2
		#cardlist = []
		#try: cardlist.extend(self.input.cards["USRBIN"])
		#except: pass
		#try: cardlist.extend(self.input.cards["EVENTBIN"])
		#except: pass
		#cardlist.sort(Input.Card.cmpPos)
		#lst = [_empty_what]
		#for i in range(len(cardlist)):
		#	card = cardlist[i]
		#	lst.append((card.sdum(), i+1))
		#lst.sort()
		#return lst

	# ----------------------------------------------------------------------
	# Detector list
	# ----------------------------------------------------------------------
	def detectorList(self, tag):
		tag = tag[:8]
		if tag in ("USRBIN","EVENTBIN"): tag = ("USRBIN","EVENTBIN")
		elif tag in ("USRTRACK","USRCOLL"): tag = ("USRTRACK","USRCOLL")
		lst = [""]
		lst.extend([x.sdum() for x in self.input.cardsSorted(tag[:8])])
		lst.extend(list(range(1,len(lst))))
		return lst

	# ----------------------------------------------------------------------
	# Return material list, excluding material e if defined
	#-----------------------------------------------------------------------
	def materialList(self, e=None):
		if self.__materialCache is not None:
			return self.__materialCache

		lst = [""] + self.input.materialList(1,True)
		if e:
			try: lst.remove(e)
			except ValueError: pass
		lst.append("@LASTMAT")
		return lst

	# ----------------------------------------------------------------------
	# DcyTimes list
	# ----------------------------------------------------------------------
	def dcytimesList(self):
		cnt = 1
		lst = []
		dcy = {}
		for card in self.input.cardsSorted("DCYTIMES"):
			for i in range(1,7):
				w = card.numWhat(i)
				if w not in dcy:
					dcy[w] = cnt
					w = card.what(i)
					if w == "": w = "0.0"
					lst.append((str(w),cnt))
					cnt += 1
		lst.append(("Semi-Analogue",-1))
		lst.append("")
		return lst

	# ----------------------------------------------------------------------
	# Return unit list
	# ----------------------------------------------------------------------
	def unitList(self, tag=None, bin=True, onlymatches=False):
		return self.input.units.filterList(tag, bin, onlymatches)

	# ----------------------------------------------------------------------
	def neutronGroupsUsed(self):
		groups = []	# List of groups to display
		try:
			# Check all LOW-BIAS cards
			for card in self.input.cards["LOW-NEUT"]:
				w1 = card.intWhat(1)
				if w1 == 0:
					g = 260
				else:
					g = w1
				if g not in groups:
					groups.append(g)
		except KeyError:
			return [260]
		return groups

	# ----------------------------------------------------------------------
	# Low energy neutron group list
	# ----------------------------------------------------------------------
	def neutronGroupList(self):
		lst = [ ("No cut-off",      0) ]
		groups = self.neutronGroupsUsed()
		for g in groups:
			i = 1
			for e in Input._neutronGroupsS.get(g, []):
				if len(groups)==1:
					lst.append((e, i))
				else:
					lst.append(("%d: %s"%(g,e), i))
				i += 1
		return lst

	# ----------------------------------------------------------------------
	# Low energy neutron materials
	# ----------------------------------------------------------------------
	def lowmatList(self):
		groups = self.neutronGroupsUsed()
		lst = [("Unknown",0)]
		for g in groups:
			i = 1
			for mat in Input._lowMaterials.get(g,[]):
				if len(groups)==1:
					lst.append(("%s. %s, %s"  \
						% (mat.elem, mat.desc, mat.temp), i))
				else:
					lst.append(("%d: %s. %s, %s"  \
						% (g, mat.elem, mat.desc, mat.temp), i))
				i += 1
			return lst

	# ----------------------------------------------------------------------
	def setWhat(self, card, what, val):
		self.flair.addUndo(self.flair.setWhatUndo(card, what, val))

	# ----------------------------------------------------------------------
	def setAbsWhat(self, card, what, val):
		self.flair.addUndo(self.flair.setAbsWhatUndo(card, what, val))

	# ----------------------------------------------------------------------
	def setSign(self, card, what, sign):
		self.flair.addUndo(self.flair.setSignUndo(card, what, sign))

	# ----------------------------------------------------------------------
	# Absolute WHATS
	# ----------------------------------------------------------------------
	def _what_abs_type(self, card, what, val):
		if val is None:
			if card.isEvalWhat(what):
				e = card.numWhat(what)
				if e > 0.0:
					return 0
				elif e < 0.0:
					return 1
				else:
					return 2

				return 2
			elif card.sign(what):
				return 1
			else:
				return 0
		else:
			if val==0:
				self.setSign(card, what, False)
			else:
				self.setSign(card, what, True)

	# ----------------------------------------------------------------------
	def _what_abs_val(self, card, what, val):
		if val is None:
			return card.absWhat(what)
		else:
			self.setAbsWhat(card, what, val)

	# ----------------------------------------------------------------------
	def _what_bit_val(self, card, what, bit, val, max_=1022):
		w = min(card.intWhat(what), max_)
		mask = 1<<bit
		if val is None:
			if w<0:
				return w
			elif w & mask > 0:
				return 1
			else:
				return 0
		else:
			if val<0:
				self.setWhat(card, what, val)
				return
			if w<0: w = 0
			if val:
				self.setWhat(card, what, w | mask)
			else:
				self.setWhat(card, what, w & ~mask)

	# ----------------------------------------------------------------------
	# WARNING reverse order!
	# ----------------------------------------------------------------------
	def _what_digit_val(self, card, what, digit, val, length=9):
		w = "%0*d"%(length, card.longWhat(what))
		if val is None:
			return int(w[digit])
		else:
			self.setWhat(card, what,
				"%s%d%s" % (w[:digit], val, w[digit+1:]))

	# ----------------------------------------------------------------------
	# What type/val Lambdas
	# ----------------------------------------------------------------------
	def _what1_type(self,c,v=None): return self._what_abs_type(c,1,v)
	def _what1_val( self,c,v=None): return self._what_abs_val( c,1,v)
	def _what2_type(self,c,v=None): return self._what_abs_type(c,2,v)
	def _what2_val( self,c,v=None): return self._what_abs_val( c,2,v)
	def _what3_type(self,c,v=None): return self._what_abs_type(c,3,v)
	def _what3_val( self,c,v=None): return self._what_abs_val( c,3,v)
	def _what4_type(self,c,v=None): return self._what_abs_type(c,4,v)
	def _what4_val( self,c,v=None): return self._what_abs_val( c,4,v)
	def _what5_type(self,c,v=None): return self._what_abs_type(c,5,v)
	def _what5_val( self,c,v=None): return self._what_abs_val( c,5,v)
	def _what6_type(self,c,v=None): return self._what_abs_type(c,6,v)
	def _what6_val( self,c,v=None): return self._what_abs_val( c,6,v)

	# ----------------------------------------------------------------------
	# Papro (Particle properties or heavy ion)
	# Used in AUXSCORE and SPECSOU
	def _papro_what(self, card, what, val=None):
		w = card.intWhat(what)
		if val is None:
			if w <= -100:
				return ""
			else:
				return card.what(what)
		else:
			if val!="" or w>-100:
				self.setWhat(card, what, val)

	# ----------------------------------------------------------------------
	def _papro_Z(self, card, what, val=None):
		w = card.intWhat(what)
		Z = 0
		A = 0
		M = 0
		if w<=-100:
			i = -w // 100
			i, Z = divmod(i, 1000)
			M, A = divmod(i, 1000)
		if val is None:
			return Z
		else:
			try: Z = int(val)
			except: pass
			self.setWhat(card, what, -(Z*100 + A*100000 + M*100000000))

	# ----------------------------------------------------------------------
	def _papro_A(self, card, what, val=None):
		w = card.intWhat(what)
		Z = 0
		A = 0
		M = 0
		if w<=-100:
			i = -w // 100
			i, Z = divmod(i, 1000)
			M, A = divmod(i, 1000)
		if val is None:
			return A
		else:
			try: A = int(val)
			except: pass
			self.setWhat(card, what, -(Z*100 + A*100000 + M*100000000))

	# ----------------------------------------------------------------------
	def _papro_M(self, card, what, val=None):
		w = card.intWhat(what)
		Z = 0
		A = 0
		M = 0
		if w<=-100:
			i = -w // 100
			i, Z = divmod(i, 1000)
			M, A = divmod(i, 1000)
		if val is None:
			return M
		else:
			try: M = int(val)
			except: pass
			self.setWhat(card, what, -(Z*100 + A*100000 + M*100000000))

	# ----------------------------------------------------------------------
	@staticmethod
	def _papro_isotope(card, what):
		w = card.what(what)
		if w != "":
			try: return (int(w)<=-100)
			except: return False
		return True

	# ----------------------------------------------------------------------
	# Papro (Particle properties or heavy ion)
	# Used in SPOTBEAM
	def _paprospot_what(self, card, val=None):
		w = card.what(1)
		if val is None:
			if w == "RAY":
				w="BEAMPART";
			if w != "":
				try:
					w = int(w)
					if w==0 or w>999: return ""
				except:
					return w
			else:
				return w
		else:
			self.setWhat(card, 1, val)

	# ----------------------------------------------------------------------
	@staticmethod
	def _paprospot_isotope(card):
		w = card.what(1)
		if w != "":
			try:
				w = int(w)
				return w==0 or w>999
			except:
				return False
		return True

	# ----------------------------------------------------------------------
	def _paprospot_Z(self, card, val=None):
		w = card.intWhat(1)
		Z = 0
		A = 0
		M = 0
		if w>=999:
			i, A = divmod(w,1000)
			M, Z = divmod(i, 100)
		if val is None:
			return Z
		else:
			try: Z = int(val)
			except: pass
			self.setWhat(card, 1, A + 1000*Z + 100000*M)

	# ----------------------------------------------------------------------
	def _paprospot_A(self, card, val=None):
		w = card.intWhat(1)
		Z = 0
		A = 0
		M = 0
		if w>=1000:
			i, A = divmod(w,1000)
			M, Z = divmod(i, 100)
		if val is None:
			return A
		else:
			try: A = int(val)
			except: pass
			self.setWhat(card, 1, A + 1000*Z + 100000*M)

	# ----------------------------------------------------------------------
	def _paprospot_M(self, card, val=None):
		w = card.intWhat(1)
		Z = 0
		A = 0
		M = 0
		if w>=1000:
			i, A = divmod(w,1000)
			M, Z = divmod(i, 100)
		if val is None:
			return M
		else:
			try: M = int(val)
			except: pass
			self.setWhat(card, 1, A + 1000*Z + 100000*M)

	# ----------------------------------------------------------------------
	# What type/val Lambdas
	# ----------------------------------------------------------------------
	def _what1_papro_what(self, c, v=None):	return self._papro_what(c,1,v)
	def _what1_papro_Z(self, c, v=None):	return self._papro_Z(c,1,v)
	def _what1_papro_A(self, c, v=None):	return self._papro_A(c,1,v)
	def _what1_papro_M(self, c, v=None):	return self._papro_M(c,1,v)

	def _what2_papro_what(self, c, v=None):	return self._papro_what(c,2,v)
	def _what2_papro_Z(self, c, v=None):	return self._papro_Z(c,2,v)
	def _what2_papro_A(self, c, v=None):	return self._papro_A(c,2,v)
	def _what2_papro_M(self, c, v=None):	return self._papro_M(c,2,v)

	# ----------------------------------------------------------------------
	# Return what to labels conversion
	# ----------------------------------------------------------------------
	def what2labels(self, layout):
		# What to labels conversion
		# Find fields to highlight
		labels = [""]*6
		for i in range(len(layout)):
			try: lbl, dat = layout[i]
			except: lbl, dat, desc = layout[i]
			if isinstance(dat,tuple):
				dat = dat[0]

			if isinstance(dat,int):
				while dat>=len(labels): labels.append("")
				labels[dat] = lbl
		return labels

#-------------------------------------------------------------------------------
# Common types
#-------------------------------------------------------------------------------
_None2		= (None, None)
_empty_what	= ("", "")
_default_what	= ("default",	 0)
_ignore_what	= ("ignore",	 0)
_reset_what	= ("reset",	-1)
_on_what	= ("On",         1)
_off_what	= ("off",        0)
_onoff_what	= [_on_what, _off_what]
_onoffempty_what= [_empty_what, _on_what, _off_what]
_onoffreset_what= [_on_what, _off_what, _reset_what]
_onoff3_what	= [_empty_what, _on_what, ("off", -1), _ignore_what]
_log_what	= [("Linear", 1), ("Log",  -1), _empty_what]
_point_what	= [("Pointwise", 1), ("Groupwise",  0), _empty_what]

__mat456_range = [ ("Mat:",	(4, LayoutManager.materialList)),
		("to Mat:",	(5, LayoutManager.materialList)),
		("Step:",	 6)]
__reg456_range = [ ("Reg:",	(4, LayoutManager.regionList)),
		("to Reg:",	(5, LayoutManager.regionList)),
		("Step:",	 6)]
_score_type = [	"",
			"EVENTBIN",
			"RESNUCLE",
			"USRBDX",
			"USRBIN",
			"USRCOLL",
			"USRTRACK",
			"USRYIELD"
               ]
_auxscore_type = [	"",
			"EVENTBIN",
			"RESNUCLE",
			"USRBDX",
			"USRBIN",
			"USRCOLL",
			"USRTRACK",
			"USRYIELD", 1,2,3,4,5
               ]

#===============================================================================
# Layouts
#===============================================================================
_generic_layout4_1 = [
		_None2,
		("#1:", 1),
		("#2:", 2),
		("#3:", 3),
	("sdum:", 0),
		("#4:", 4),
		("#5:", 5),
		("#6:", 6)]

def _generic_layoutN(nwhats):
	n,m = divmod(nwhats-1,6)
	if m==0: n -= 1
	layout = _generic_layout4_1[:]
	i = 7
	while n>0:
		n -= 1

		layout.append(_None2)
		for j in range(3):
			layout.append(("#%d:"%(i),i))
			i += 1

		layout.append(_None2)
		for j in range(3):
			layout.append(("#%d:"%(i),i))
			i += 1
	return layout

def _generic_layout4(card):
	if card.nwhats()<=7:
		return _generic_layout4_1
	else:
		return _generic_layoutN(card.nwhats())

_generic_layout = _generic_layout4
_empty_layout4  = [ _None2, _None2, _None2, _None2 ]

#-------------------------------------------------------------------------------
_BAMJET_layout4		= _generic_layout4	# Developers
_DPM_PARA_layout4	= _generic_layout4	# Developers
_RQMD_layout4		= _generic_layout4	# Developers

_BME_layout4		= _generic_layout4	# not implemented yet
_EXPTRANS_layout4	= _generic_layout4	# not implemented yet
_MYRQMD_layout4		= _generic_layout4	# not implemented yet

#================================ Special Cards ================================
_P_DEFINE_layout4 = [ (None,0), (":",1), _None2, _None2 ]
_P_UNDEF_layout4  = [ (None,(0, LayoutManager.defineList)), _None2, _None2, _None2 ]
_P_IF_layout4     = _P_UNDEF_layout4
_P_IFDEF_layout4  = _P_IF_layout4
_P_IFNDEF_layout4 = _P_IF_layout4
_P_ELIF_layout4   = _P_UNDEF_layout4
_P_ELSE_layout4   = _empty_layout4
_P_ENDIF_layout4  = _empty_layout4

_P_INCLUDE_layout4 = [	#(None, 0),
			(None, (0, (("Fluka input",("*.inp","*.fluka")),("Include","*.inc"),("All","*")))),
			_None2,
			_None2,
			_None2 ]

#-------------------------------- default layout -------------------------------
_SOURCE_layout4		= _generic_layoutN(18)
_USRGCALL_layout4	= _generic_layout4
_USRICALL_layout4	= _generic_layout4
_USROCALL_layout4	= _generic_layout4

#================================ Geometry Cards ===============================
#-------------------------------------------------------------------------------
#-------------------------------- RPP ------------------------------------------
_RPP_layout4 = [
	(None,			0),
		("Xmin:",	1),
		("Xmax:",	2),
		_None2,
	_None2,
		("Ymin:",	3),
		("Ymax:",	4),
		_None2,
	_None2,
		("Zmin:",	5),
		("Zmax:",	6),
		_None2 ]

#-------------------------------- BOX ------------------------------------------
_BOX_layout4 = [
	(None,			0),
		("x:",		1),
		("y:",		2),
		("z:",		3),
	_None2,
		("H1x:",	4),
		("H1y:",	5),
		("H1z:",	6),
	_None2,
		("H2x:",	7),
		("H2y:",	8),
		("H2z:",	9),
	_None2,
		("H3x:",	10),
		("H3y:",	11),
		("H3z:",	12) ]

#-------------------------------- SPH ------------------------------------------
_SPH_layout4 = [
	(None,			0),
		("x:",		1),
		("y:",		2),
		("z:",		3),
	_None2,
		("R:",		4),
		_None2,
		_None2 ]

#-------------------------------- RCC ------------------------------------------
_RCC_layout4 = [
	(None,			0),
		("x:",		1),
		("y:",		2),
		("z:",		3),
	_None2,
		("Hx:",		4),
		("Hy:",		5),
		("Hz:",		6),
	_None2,
		("R:",		7),
		_None2,
		_None2 ]

#-------------------------------- REC ------------------------------------------
_REC_layout4 = [
	(None,			0),
		("x:",		1),
		("y:",		2),
		("z:",		3),
	_None2,
		("Hx:",		4),
		("Hy:",		5),
		("Hz:",		6),
	_None2,
		("R1x:",	7),
		("R1y:",	8),
		("R1z:",	9),
	_None2,
		("R2x:",	10),
		("R2y:",	11),
		("R2z:",	12) ]

#-------------------------------- TRC ------------------------------------------
_TRC_layout4 = [
	(None,			0),
		("x:",		1),
		("y:",		2),
		("z:",		3),
	_None2,
		("Hx:",		4),
		("Hy:",		5),
		("Hz:",		6),
	_None2,
		("Rbase:",	7),
		("Rapex:",	8),
		_None2 ]

#-------------------------------- TRXYZ ----------------------------------------
_TRXYZ_layout4 = [
	(None,			0),
		("x:",		1),
		("y:",		2),
		("z:",		3),
	_None2,
		("a:",		4),
		("b:",		5),
		("c:",		6) ]

#-------------------------------- ELL ------------------------------------------
_ELL_layout4 = [
	(None,		0),
		("F1x:",	1),
		("F1y:",	2),
		("F1z:",	3),
	_None2,
		("F2x:",	4),
		("F2y:",	5),
		("F2z:",	6),
	_None2,
		("L:",		7),
		_None2,
		_None2 ]
#-------------------------------- TET ------------------zxw20240827--For TET, added by zxw
_TET_layout4 = [
	(None,			0),
		("V1x:",	1),
		("V1y:",	2),
		("V1z:",	3),
	_None2,
		("V2x:",	4),
		("V2y:",	5),
		("V2z:",	6),
	_None2,
		("V3x:",	7),
		("V3y:",	8),
		("V3z:",	9),
	_None2,
		("V4x:",	10),
		("V4y:",	11),
		("V4z:",	12) ]

#-------------------------------- WED ------------------------------------------
_WED_layout4 = [
	(None,			0),
		("x:",		1),
		("y:",		2),
		("z:",		3),
	_None2,
		("H1x:",	4),
		("H1y:",	5),
		("H1z:",	6),
	_None2,
		("H2x:",	7),
		("H2y:",	8),
		("H2z:",	9),
	_None2,
		("H3x:",	10),
		("H3y:",	11),
		("H3z:",	12) ]

_RAW_layout4 = _WED_layout4

#-------------------------------- ARB ------------------------------------------
_ARB_layout4 = [
	(None,			0),
		("V1x:",	1),
		("V1y:",	2),
		("V1z:",	3),
	_None2,
		("V2x:",	4),
		("V2y:",	5),
		("V2z:",	6),
	_None2,
		("V3x:",	7),
		("V3y:",	8),
		("V3z:",	9),
	_None2,
		("V4x:",	10),
		("V4y:",	11),
		("V4z:",	12),
	_None2,
		("V5x:",	13),
		("V5y:",	14),
		("V5z:",	15),
	_None2,
		("V6x:",	16),
		("V6y:",	17),
		("V6z:",	18),
	_None2,
		("V7x:",	19),
		("V7y:",	20),
		("V7z:",	21),
	_None2,
		("V8x:",	22),
		("V8y:",	23),
		("V8z:",	24),
	_None2,
		("F1:",		25),
		("F2:",		26),
		("F3:",		27),
	_None2,
		("F4:",		28),
		("F5:",		29),
		("F6:",		30) ]

#-------------------------------- XYP ------------------------------------------
_XYP_layout4 = [
	(None,			0),
		("z:",		1),
		_None2,
		_None2 ]

#-------------------------------- YZP ------------------------------------------
_YZP_layout4 = [
	(None,			0),
		("x:",		1),
		_None2,
		_None2 ]

#-------------------------------- XZP ------------------------------------------
_XZP_layout4 = [
	(None,			0),
		("y:",		1),
		_None2,
		_None2 ]

#-------------------------------- PLA ------------------------------------------
_PLA_layout4 = [
	(None,			0),
		("Nx:",		1),
		("Ny:",		2),
		("Nz:",		3),
	_None2,
		("x:",		4),
		("y:",		5),
		("z:",		6) ]

#-------------------------------- QUA ------------------------------------------
_QUA_layout4 = [
	(None,			0),
		("Cxx:",	1),
		("Cyy:",	2),
		("Czz:",	3),
	_None2,
		("Cxy:",	4),
		("Cxz:",	5),
		("Cyz:",	6),
	_None2,
		("Cx:",		7),
		("Cy:",		8),
		("Cz:",		9),
	_None2,
		("C:",		10),
		_None2,
		_None2 ]

#-------------------------------- XCC ------------------------------------------
_XCC_layout4 = [
	(None,			0),
		("y:",		1),
		("z:",		2),
		("R:",		3) ]

#-------------------------------- YCC ------------------------------------------
_YCC_layout4 = [
	(None,			0),
		("z:",		1),
		("x:",		2),
		("R:",		3) ]

#-------------------------------- ZCC ------------------------------------------
_ZCC_layout4 = [
	(None,			0),
		("x:",		1),
		("y:",		2),
		("R:",		3) ]

#-------------------------------- XEC ------------------------------------------
_XEC_layout4 = [
	(None,			0),
		("y:",		1),
		("z:",		2),
		_None2,
	_None2,
		("Ly:",		3),
		("Lz:",		4),
		_None2 ]

#-------------------------------- YEC ------------------------------------------
_YEC_layout4 = [
	(None,			0),
		("z:",		1),
		("x:",		2),
		_None2,
	_None2,
		("Lz:",		3),
		("Lx:",		4),
		_None2 ]

#-------------------------------- ZEC ------------------------------------------
_ZEC_layout4 = [
	(None,			0),
		("x:",		1),
		("y:",		2),
		_None2,
	_None2,
		("Lx:",		3),
		("Ly:",		4),
		_None2 ]

#-------------------------------- $Start_expansion -----------------------------
_START_EXPANSION_layout4 = [
	_None2,
		("f:",		1),
		_None2,
		_None2 ]

#-------------------------------- $End_expansion -------------------------------
_END_EXPANSION_layout4 = _empty_layout4

#-------------------------------- $Start_translat ------------------------------
_START_TRANSLAT_layout4 = [
	_None2,
		("dx:",		1),
		("dy:",		2),
		("dz:",		3) ]

#-------------------------------- $End_translat --------------------------------
_END_TRANSLAT_layout4 = _empty_layout4

#-------------------------------- $Start_transform -----------------------------
_START_TRANSFORM_layout4 = [
	_None2,
		( "Trans:",	(1, LayoutManager.rotdefiList)),
		_None2,
		_None2 ]

#-------------------------------- $End_transform -------------------------------
_END_TRANSFORM_layout4 = _empty_layout4

#-------------------------------- REGION ---------------------------------------
def __REGION_volume(layout, card, val=None):
	if val is None:
		try: return card[Input._VOLUME]
		except: return ""
	else:
		card[Input._VOLUME] = val

__REGION_layout4_definition = [
		(None,		0),
		_None2,
		("Neigh:",	1),
		_None2,
#		("Volume:",	__REGION_volume, "Region volume cm3"),
	("expr:",		-1),
		_None2,
		_None2,
		_None2 ]

__REGION_layout4_cont = [ (None, -1), _None2, _None2, _None2 ]

def _REGION_layout4(card):
	if card.sdum() == "&":
		return __REGION_layout4_cont
	else:
		return __REGION_layout4_definition

#===================== Normal Cards ===================
#------------------------------- ASSIGNMAt -------------------------------------
_ASSIGNMA_field = [	_empty_what,
			("Magnetic", 1),
			("Electric", 2),
			("Magnetic+Electric", 3),
			("Prompt Magnetic", 4),
			("Prompt Electric", 5),
			("Prompt Magnetic+Electric", 6),
			("Decay Magnetic", 7),
			("Decay Electric", 8),
			("Decay Magnetic+Electric", 9),
			_off_what,
			_reset_what ]

_ASSIGNMA_layout4 = [ _None2,
		("Mat:",	(1, LayoutManager.materialList)),
		("Reg:",	(2, LayoutManager.regionList)),
		("to Reg:",	(3, LayoutManager.regionList)),
	_None2,
		("Mat(Decay):",	(6, LayoutManager.materialList)),
		("Step:",	 4),
		("Field:",	(5, _ASSIGNMA_field))]

#------------------------------- AUXSCORE --------------------------------------
def _AUXSCORE_layout4(card):
	layout = [ _None2,
		("Type:",	(1, _auxscore_type)),
		("Part:",       (LayoutManager._what2_papro_what, Input.Particle.particles),
				"Particle type or heavy ion"),
		("Set:",	(0, ["", "AMB74", "AMBDS",
				"EAP116", "EIS116", "EPA116",
				"EAP74", "EAPMP", "ERT74",
				"ERTMP", "EWT74", "EWTMP"])),
		("Delta Ray:",	(3, [_empty_what, ("assign to self",0), ("assign to parent", 1), _reset_what])) ]

	if LayoutManager._papro_isotope(card, 2):
		layout.extend([	("Z:",	LayoutManager._what2_papro_Z, "Heavy ion charge"),
				("A:",	LayoutManager._what2_papro_A, "Heavy ion mass"),
				("Isomer:",	LayoutManager._what2_papro_M, "Heavy ion isomer"),
				_None2])

	layout.extend([("Det:",	(4, lambda l,c=card:l.detectorList(c.what(1)))),
		    ("to Det:",	(5, lambda l,c=card:l.detectorList(c.what(1)))),
		    ("Step:",	 6)])
	return layout

#------------------------------- BEAM ------------------------------------------
__BEAM_type0 = [ ("Momentum",0), ("Energy",1), ("Function",2) ]
__BEAM_type1 = [ ("Flat",0), ("Gauss",1), ("Function",2) ]
__BEAM_type2 = [ ("Flat",0), ("Gauss",1), ("Isotropic",2), ("Function",3) ]
__BEAM_type3 = [ ("Rectangular",0), ("Annular",1), ("Gauss",2), ("Function", 3) ]
__BEAM_type4 = [ ("Rectangular",0), ("Gauss",2), ("Function", 3) ]
__BEAM_sdum  = [" ", "AMBE", "AMB", "252CF", "D-D", "D-T"]

# Generate a dynamic layout
def _BEAM_layout4_nsource(card):
	layout =  [ _None2,
			("Type:",	(0, __BEAM_sdum)),
#			("Type:",	(0, Input.Particle.beam)),
			("Max energy:", 1),_None2]
	# 3rd line
	__BEAM_shape_layout(card, layout)
#			("#2:", 2),
#		("#3:", 3),
#			("#4:", 4),
#			("#5:", 5),
#			("#6:", 6)]
	return layout

# Generate a dynamic layout
def _BEAM_layout4_ndeut(card):
	layout =  [ _None2,
			("Type:",	(0, __BEAM_sdum)),
#			("Type:",	(0, Input.Particle.beam)),
			("Deuteron energy:", 1),_None2]
	# 3rd line
	__BEAM_shape_layout(card, layout)
#			("#2:", 2),
#		("#3:", 3),
#			("#4:", 4),
#			("#5:", 5),
#			("#6:", 6)]
	return layout

# Beam divergence
def __BEAM_div_type(layout, card, val=None, w=3):
	if val is None:
		if card.isEvalWhat(w):
			e = card.numWhat(w)
			if e > 0.0:
				#if e >= 2000.0*math.pi:
				#	return 2
				return 0
			elif e < 0.0:
				return 1
			else:
				return 3
		elif card.sign(w):
			return 1
		elif abs(card.numWhat(w)) >= 2000.0*math.pi:
			return 2
		return 0
	else:
		n = card.absWhat(w)
		try: n = float(n)
		except: n = 0.0
		if val==2:
			if n <= 2000.0*math.pi:
				layout.setWhat(card, w, 10000.0+n)
		else:
			if n >= 2000.0*math.pi:
				n -= 10000.0
				if n >= 2000.0*math.pi: n = 0.0
				layout.setWhat(card, w, n)

			if val==0:
				layout.setSign(card, w, False)
			elif val==1:
				layout.setSign(card, w, True)

def __BEAM_shape_type3(layout, card, val=None):
	if val is None:
		if card.isEvalWhat(4):
			e = card.numWhat(4)
			if e < 0.0:
				return 2
			elif e > 0.0:
				if card.numWhat(6)<0.0:
					return 1
				return 0
			else:
				return 3
		elif card.sign(4):
			return 2
		else:
			if card.sign(6):
				return 1
			return 0
	else:
		if val==0:	# Flat
			undoinfo = [
				layout.flair.setSignUndo(card, 4,  False),
				layout.flair.setSignUndo(card, 6,  False)]
		elif val==1:	# Annular
			if card.intWhat(6) == 0:
				undoinfo = [layout.flair.setWhatUndo(card, 6, 1.0)]
			else:
				undoinfo = []
			undoinfo.extend([
				layout.flair.setSignUndo(card, 4,  False),
				layout.flair.setSignUndo(card, 5,  False),
				layout.flair.setSignUndo(card, 6,  True) ])
		else:		# Gauss
			undoinfo = [
				layout.flair.setSignUndo(card, 4,  True),
				layout.flair.setSignUndo(card, 6,  False)]
		layout.flair.addUndo(undo.createListUndo(undoinfo))

def __BEAM_shape_type4(layout, card, val=None):
	if val is None:
		if card.isEvalWhat(5):
			e = card.numWhat(5)
			if e < 0.0:
				return 2
			elif e > 0.0:
				return 0
			else:
				return 3
		elif card.sign(5):
			return 2
		else:
			return 0
	else:
		if val==0:	# Flat
			undoinfo = [
				layout.flair.setSignUndo(card, 5,  False),
				layout.flair.setSignUndo(card, 6,  False)]
		else:
			undoinfo = [
				layout.flair.setSignUndo(card, 5,  True),
				layout.flair.setSignUndo(card, 6,  False)]
		layout.flair.addUndo(undo.createListUndo(undoinfo))

# Extend layout with the beam shape
def __BEAM_shape_layout(card, layout):
	layout.append(("Shape(X):", (__BEAM_shape_type3, __BEAM_type3),
			"Beam shape in X axis"))
	if card.sign(6):
		layout.append(("Rmin:",	 LayoutManager._what5_val,
				"WHAT(5): Minimum radius of annular beam [cm]"))
		layout.append(("Rmax:",	 LayoutManager._what4_val,
				"WHAT(4): Maximum radius of annular beam [cm]"))
		layout.append(_None2)
	else:
		if card.sign(4):
			layout.append(("x(FWHM):",  LayoutManager._what4_val,
					"WHAT(4): Gaussian FWHM in x-direction [cm]"))
		else:
			layout.append(("\u0394x:",	 LayoutManager._what4_val,
					"WHAT(4): Beam width in x-direction [cm]"))

		layout.append(("Shape(Y):", (__BEAM_shape_type4, __BEAM_type4),
				"Beam shape in Y axis"))
		if card.sign(5):
			layout.append(("y(FWHM):",	 LayoutManager._what5_val,
					"WHAT(5): Gaussian FWHM in y-direction [cm]"))
		else:
			layout.append(("\u0394y:",	 LayoutManager._what5_val,
					"WHAT(5): Beam width in y-direction [cm]"))

# Generate a dynamic layout
def _BEAM_layout4_generic(card):
	# 1st line
	layout = [ _None2,
		("Beam:", (LayoutManager._what1_type, __BEAM_type0),
			"WHAT(1) sign: Energy or Momentum") ]

	if card.sign(1):
		layout.append(("E:",	 LayoutManager._what1_val,
				"WHAT(1): Energy [GeV or GeV/n for HI]"))
	else:
		layout.append(("p:",	 LayoutManager._what1_val,
				"WHAT(1): Momentum [GeV/c]"))

	layout.append(("Part:",	(0, Input.Particle.beam)))

	# 2nd line
	layout.append(("\u0394p:",	(LayoutManager._what2_type, __BEAM_type1),
			"WHAT(2): Momentum spread type"))
	if card.sign(2):
		layout.append(("\u0394p(FWHM):", LayoutManager._what2_val,
			"WHAT(2): Momentum spread Gaussian [GeV/c] FWHM=2.355*sigma"))
	else:
		layout.append(("\u0394p:",	LayoutManager._what2_val,
			"WHAT(2): Momentum spread FLAT [GeV/c]"))

	layout.append(("\u0394\u03d5:",	(__BEAM_div_type, __BEAM_type2),
			"WHAT(3): Divergence type"))
	if card.numWhat(3) >= 2000.0*math.pi:
		layout.append(_None2)
	elif card.sign(3):
		layout.append(("\u0394\u03d5 (FWHM):", LayoutManager._what3_val,
				"WHAT(3): Divergence [mrad]"))
	else:
		layout.append(("\u0394\u03d5:", LayoutManager._what3_val,
				"WHAT(3): Divergence [mrad]"))

	# 3rd line
	__BEAM_shape_layout(card, layout)

	# 3rd line
#	layout.append(("Weight:", LayoutManager._what6_val,
#			"WHAT(6): Particle weight. WARNING is important for annular beam"))

	return layout

def _BEAM_layout4(card):
	sdum = card.sdum()[:8]
	if   sdum == "AMB" or sdum == "AMBE" or sdum == "252CF" :  
		return _BEAM_layout4_nsource(card)
	elif sdum == "D-D" or sdum == "D-T":
		return _BEAM_layout4_ndeut(card)
	else:
		return _BEAM_layout4_generic(card)


#------------------------------- BEAMAXES --------------------------------------
_BEAMAXES_layout4 = [ _None2,
		("cosBxx:",	1),
		("cosBxy:",	2),
		("cosBxz:",	3),
	_None2,
		("cosBzx:",	4),
		("cosBzy:",	5),
		("cosBzz:",	6) ]

#------------------------------- BEAMPOS ---------------------------------------
__BEAMPOS_sdum = [("POSITIVE",""),
		"NEGATIVE",
		"FLOOD",
		"CART-VOL",
		"CYLI-VOL",
		"SPHE-VOL" ]
__BEAMPOS_layout4 = [ _None2,
		("x:",		1),
		("y:",		2),
		("z:",		3),
	_None2,
		("cosx:",	4),
		("cosy:",	5),
		("Type:",	(0, __BEAMPOS_sdum))]

__BEAMPOS_flood_layout4 = [ _None2,
		("R:",		1),
		_None2,
		("Type:",	(0, __BEAMPOS_sdum))]

__BEAMPOS_cart_layout4 = [ _None2,
		("Xin:",	1),
		("Xout:",	2),
		("Type:",	(0, __BEAMPOS_sdum)),
	_None2,
		("Yin:",	3),
		("Yout:",	4),
		_None2,
	_None2,
		("Zin:",	5),
		("Zout:",	6),
		_None2 ]

__BEAMPOS_cyli_layout4 = [ _None2,
		("Rin:",	1),
		("Rout:",	2),
		("Type:",	(0, __BEAMPOS_sdum)),
	_None2,
		("Hin:",	3),
		("Hout:",	4),
		_None2 ]

__BEAMPOS_sphe_layout4 = [ _None2,
		("Rin:",	1),
		("Rout:",	2),
		("Type:",	(0, __BEAMPOS_sdum))]

def _BEAMPOS_layout4(card):
	if card.sdum() == "FLOOD":
		return __BEAMPOS_flood_layout4
	elif card.sdum() == "CART-VOL":
		return __BEAMPOS_cart_layout4
	elif card.sdum() == "CYLI-VOL":
		return __BEAMPOS_cyli_layout4
	elif card.sdum() == "SPHE-VOL":
		return __BEAMPOS_sphe_layout4
	else:
		return __BEAMPOS_layout4

#------------------------------- BIASING ---------------------------------------
__BIASING_types = [	_empty_what,
			("All regions",    -1),
			("All particles",   0),
			("Hadrons & Muons", 1),
			("e-e+,\u03b3",   2),
			("Low neutrons",   3) ]
__BIASING_sdum = [ "", "PRINT", "NOPRINT", "USER", "NOUSER", "RRPRONLY" ]
__BIASING_prim = [ "", "PRIMARY", "NOPRIMAR" ]

__BIASING_layout4_regions = [ _None2,
		("Type:",	(1, __BIASING_types)),
		("RR:",		2),
		("Imp:",	3),
	("Opt:",		(0, __BIASING_sdum)) ] \
		+ __reg456_range

__BIASING_layout4_particles   = [ _None2,
		("Type:",	(1, __BIASING_types)),
		("Mod. M",	2),
		_None2,
	("Opt:",		(0, __BIASING_prim)),
		("Part:",	(3, Input.Particle.list)),
		("to Part:",	(4, Input.Particle.list)),
		("Step:",	 5) ]

def _BIASING_layout4(card):
	if card.intWhat(1) < 0:
		return __BIASING_layout4_particles
	else:
		return __BIASING_layout4_regions

#------------------------------- COMPOUND --------------------------------------
__COMPOUND_list = [ ("Atom", 0), ("Mass", 1), ("Volume", 2) ]

def __COMPOUND_components(layout, card, val=None):
	n = (card.nwhats()-1)//2
	if val is None:
		return n
	elif val != n:
		nn = int(val)*2+1

		s1 = card.sign(1)
		s2 = card.sign(2)

		undoinfo = [ layout.flair.setNWhatsUndo(card, nn) ]
		for i in range(3,card.nwhats(),2):
			undoinfo.append( layout.flair.setSignUndo(card, i,   s1))
			undoinfo.append( layout.flair.setSignUndo(card, i+1, s2))
		layout.flair.addUndo(undo.createListUndo(undoinfo))

def __COMPOUND_type(layout, card, val=None):
	if val is None:
		s1 = card.sign(1)
		s2 = card.sign(2)
		if not s1 and not s2:
			return 0
		elif s1 and not s2:
			return 1
		return 2
	else:
		if val==0:
			s1 = False;	n1 =  1.0
			s2 = False;	n2 =  1.0
		elif val==1:
			s1 = True;	n1 = -1.0
			s2 = False;	n2 =  1.0
		else:
			s1 = True;	n1 = -1.0
			s2 = True;	n2 = -1.0

		undoinfo = []
		for i in range(1,card.nwhats(),2):
			undoinfo.append(layout.flair.setSignUndo(card, i,   s1))
			undoinfo.append(layout.flair.setSignUndo(card, i+1, s2))
		layout.flair.addUndo(undo.createListUndo(undoinfo))

__COMPOUND_elements = [ ("%d..%d"%(x-2,x),x) for x in range(3,82,3) ]

def _COMPOUND_layout4(card):
	if card.nwhats()<7:
		card.setNWhats(7)

	layout = [ _None2,
		("Name:",	(0, LayoutManager.materialList)),
		("Mix:",	(__COMPOUND_type, __COMPOUND_list),
					"Type of fraction: Atom/Mass/Volume"),
		("Elements:",	(__COMPOUND_components, __COMPOUND_elements),
					"Number of elements defining compound") ]

	n = card.nwhats()-1
	i = 0
	what = 0
	while what < n:
		i += 1
		what += 1
		layout.append(("f%d:"%(i),
				lambda l,c,v=None,w=what:l._what_abs_val(c,w,v),
				"Fraction of material #%d"%(i)))
		what += 1
		layout.append(
			("M%d:"%(i),
				(lambda l,c,v=None,w=what:l._what_abs_val(c,w,v),
				lambda l,c=card: l.materialList(c.sdum())),
				"Material #%d name"%(i)))
	if (n//2)%2:
		layout.extend([_None2, _None2])

	return layout

#------------------------------- CORRFACT --------------------------------------
_CORRFACT_layout4 = [ _None2,
		("x\u03c1 (dE/dx):",	1),
		("x\u03c1 other:",	2),
		_None2,
	_None2,
		("Reg:",	(4, LayoutManager.regionList)),
		("to Reg:",	(5, LayoutManager.regionList)),
		("Step:",	 6)]

#------------------------------- DCYTIMES --------------------------------------
_DCYTIMES_layout4 = [ _None2,
		("t1:",		1),
		("t2:",		2),
		("t3:",		3),
		_None2,
		("t4:",		4),
		("t5:",		5),
		("t6:",		6) ]

#------------------------------- DCYSCORE --------------------------------------
def _DCYSCORE_layout4(card):
	return [ _None2,
		("Cooling t:",	(1, LayoutManager.dcytimesList)),
		("Kind:",	(0, _score_type)),
		_None2,
		_None2,
		("Det:",	(4, lambda l,c=card:l.detectorList(c.sdum()))),
		("to Det:",	(5, lambda l,c=card:l.detectorList(c.sdum()))),
		("Step:",	6) ]

#------------------------------- DETECT - ---------------------------------------
# FIXME help
__DETECT_type = [ ("Detector", 0), ("Combination", 1) ]
__DETECT_coin = [ _empty_what, ("Coincidence", 1), ("AntiCoincidence", -1) ]

def __DETECT_regions(layout, card, val=None):
	#                 sdum whats first regs firstreg
	n = ((card.nwhats()-1) //6    -1)  * 5  +1	# Ouf!
	if val is None:
		return max(n, 1)
	else:
		val = int(val)
		if val != n:
			nn = ((val-1)//5+1)*6+1
			undoinfo = [ layout.flair.setNWhatsUndo(card, nn) ]
			w1 = card.what(1)
			for i in range(7, nn, 6):
				undoinfo.append(layout.flair.setWhatUndo(card, i, w1))
			layout.flair.addUndo(undo.createListUndo(undoinfo))

def __DETECT_regions2(layout, card, val=None):
	#                sdum whats regs
	n = (card.nwhats()-1)//6    * 5
	if val is None:
		return n
	else:
		val = int(val)
		if val != n:
			nn = (val//5)*6+1
			undoinfo = [ layout.flair.setNWhatsUndo(card, nn) ]
			w1 = card.what(1)
			for i in range(7, nn, 6):
				undoinfo.append(layout.flair.setWhatUndo(card, i, w1))
			layout.flair.addUndo(undo.createListUndo(undoinfo))

def _DETECT_layout4(card):
	if card.intWhat(1)==0:
		layout = [ _None2,
			("Type:",	(1, __DETECT_type)),
			("Regions:",	(__DETECT_regions, list(range(1,12,5)))),
			("Name:",	0),

			("Emin:",	2),
			("Emax:",	3),
			("Ecut:",	4),
			("Trigger:",	(5, __DETECT_coin)),

			("Reg1:",	(6, LayoutManager.signRegionList))]
		j = 8
		for i in range(2, ((card.nwhats()-1)//6-1)* 5+2):
			layout.append(("Reg%d:"%(i), (j, LayoutManager.signRegionList)))
			j += 1
			if (j-1)%6 == 0: j += 1
		if len(layout)%4 != 0:
			layout.extend([_None2] * (4-len(layout)%4))
		return layout

	else:
		layout = [ _None2,
			("Type:",	(1, __DETECT_type)),
			("Regions:",	(__DETECT_regions2, list(range(5,11,5)))),
			("Name:",	0)]
		j = 2
		for i in range(1, (card.nwhats()-1)//6*5+1):
			layout.append(("Reg%d:"%(i), (j, LayoutManager.signRegionList)))
			j += 1
			if (j-1)%6 == 0: j += 1
		if len(layout)%4 != 0:
			layout.extend([_None2] * (4-len(layout)%4))
		return layout

#------------------------------- DISCARD ---------------------------------------
# FIXME Heavy ions
_DISCARD_layout4 = [ _None2,
		("p1:",	(1, Input.Particle.signedListAll)),
		("p2:",	(2, Input.Particle.signedListAll)),
		("p3:",	(3, Input.Particle.signedListAll)),
	_None2,
		("p4:",	(4, Input.Particle.signedListAll)),
		("p5:",	(5, Input.Particle.signedListAll)),
		("p6:",	(6, Input.Particle.signedListAll)) ]

#------------------------------- DEFAULTS --------------------------------------
__DEFAULTS_sdum = [ "CALORIME", "DAMAGE", "EM-CASCA", "EET/TRAN", "HADROTHE",
		"ICARUS", "NEUTRONS", "NEW-DEFA", "PRECISIO",
		"SHIELDIN" ]

_DEFAULTS_layout4 = [	_None2,
			("",	(0, __DEFAULTS_sdum)),
			_None2,	_None2 ]

#------------------------------- DELTARAY --------------------------------------
_DELTARAY_layout4 = [ _None2,
		("E thres:",	1),
		("# Log dp/dx:", 2),
		("Log width dp/dx:", 3),
	("Print",	(0, ["NOPRINT", "PRINT"])) ] \
		+ __mat456_range

#------------------------------- DPMJET ----------------------------------------
__DPMJET_verb = [ _empty_what,
		("Minimum", 0),
		("Minor Errors", 1),
		("Common Block", 2),
		("Final Residual", 3),
		("Final State Obj", 4)]
_DPMJET_layout4	= [
		_None2,
		("Format:", 1),
		("Unit:", (2, lambda l: l.unitList("DPMJET",False))),
		("Verbose:", (3, __DPMJET_verb)),
		_None2]

#------------------------------- ELCFIELD --------------------------------------
_ELCFIELD_layout4 = [ _None2,
		("Max Ang (deg):",	1),
		("Bound Acc. (cm):",	2),
		("Min step (cm):",	3),
	_None2,
		("Ex:",		4),
		("Ey:",		5),
		("Ez:",		6)]

#------------------------------- EMF -------------------------------------------
__EMF_switch = [ ("ON", ""), ("OFF", "EMF-OFF") ]
_EMF_layout4 = [ _None2,
		("",		(0, __EMF_switch)),
		_None2,
		_None2 ]

#------------------------------- EMF-BIAS --------------------------------------
__EMF_BIAS_sdum = [ "", "LPBEMF",  "LAMBEMF",  "LAMBCOMP", "LAMBBREM",
		"LBRREMF", "LBRRCOMP", "LBRRBREM" ]

__EMF_BIAS_layout4_header = [ _None2, ("Type:",	(0, __EMF_BIAS_sdum)) ]

__EMF_BIAS_layout4_LPB = __EMF_BIAS_layout4_header \
		+ [("Ethr e-e+:", 2),
		   ("Ethr \u03b3:", 3),
		("Old bremss.:", (lambda l,c,v=None,b=0:l._what_bit_val(c,1,b,v),
					_onoffreset_what),
				"LPB activated for bremsstrahlung and pair production (old default)"),
		("Bremsstrahlung:",	(lambda l,c,v=None,b=1:l._what_bit_val(c,1,b,v),
						_onoffreset_what),
				"LPB activated for bremsstrahlung"),
		("Pair Prod.:",		(lambda l,c,v=None,b=2:l._what_bit_val(c,1,b,v),
						_onoffreset_what),
				"LPB activated for pair production"),
		("e+ ann @rest:",	(lambda l,c,v=None,b=3:l._what_bit_val(c,1,b,v),
						_onoffreset_what),
				"LPB activated for positron annihilation at rest"),
		("Compton:",		(lambda l,c,v=None,b=4:l._what_bit_val(c,1,b,v),
						_onoffreset_what),
				"LPB activated for Compton scattering"),
		("Bhabha&Moller:",	(lambda l,c,v=None,b=5:l._what_bit_val(c,1,b,v),
						_onoffreset_what),
				"LPB activated for Bhabha & Moller scattering"),
		("Photo-electric:",	(lambda l,c,v=None,b=6:l._what_bit_val(c,1,b,v),
						_onoffreset_what),
				"LPB activated for photoelectric effect"),
		("e+ ann @flight:",	(lambda l,c,v=None,b=7:l._what_bit_val(c,1,b,v),
						_onoffreset_what),
				"LPB activated for positron annihilation in flight"),
		 _None2 ] + __reg456_range

__EMF_BIAS_layout4_other = __EMF_BIAS_layout4_header + [
		("e <Path>:",	 1),
		("\u03B3 <Path>:", 2),
	("Generation:",		 3) ] \
		+ __mat456_range

def _EMF_BIAS_layout4(card):
	sdum = card.sdum()
	if sdum == "LPBEMF" or sdum == "":
		return __EMF_BIAS_layout4_LPB
	else:
		return __EMF_BIAS_layout4_other

#------------------------------- EMFCUT ----------------------------------------
__EMFCUT_sdum = [ ("transport",""),
			"ANNH-THR", "DELAYED", "ELPO-THR", "PHO2-THR",
			"PHOT-THR", "PROD-CUT", "PROMPT"]

__EMFCUT_ee = [ ("Total",0), ("Kinetic",1), ("Function",2) ]

__EMFCUT_layout4_elpothr = [ _None2,
		("Type:", (0, __EMFCUT_sdum)),
		("e-e+ Brem:",	1),
		("\u03b3 Bhabha/Moller:",	2),
	("Photonuc:",	3) ] \
		+ __mat456_range

__EMFCUT_layout4_annhthr = [ _None2,
		("Type:", (0, __EMFCUT_sdum)),
		("Ekin e+ ann:",	1),
		_None2,
		_None2 ] \
		+ __mat456_range

__EMFCUT_layout4_photthr = [ _None2,
		("Type:", (0, __EMFCUT_sdum)),
		("Compton:",	1),
		("Photo-elec:",	2),
	("Pair prod:",	3) ] \
		+ __mat456_range

__EMFCUT_layout4_pho2thr = [ _None2,
		("Type:", (0, __EMFCUT_sdum)),
		("Rayleigh:",	1),
		("Photonuc:",	2),
	_None2 ] \
		+ __mat456_range

def _EMFCUT_layout4(card):
	sdum = card.sdum()
	if sdum=="ELPO-THR":
		return __EMFCUT_layout4_elpothr
	elif sdum=="ANNH-THR":
		return __EMFCUT_layout4_annhthr
	elif sdum=="PHOT-THR":
		return __EMFCUT_layout4_photthr
	elif sdum=="PHO2-THR":
		return __EMFCUT_layout4_pho2thr
	else:
		layout = [ _None2,
			("Type:",	(0, __EMFCUT_sdum)),
			_None2,
			_None2,
			_None2,
			("e-e+ Threshold:",	(LayoutManager._what1_type, __EMFCUT_ee))]

		if card.sign(1):
			layout.append(("e-e+ Ekin:",	 LayoutManager._what1_val,
					"WHAT(1): Electrons Kinetic Energy [GeV]"))
		else:
			layout.append(("e-e+ E:",	 LayoutManager._what1_val,
					"WHAT(1): Electrons Total Energy [GeV]"))
		layout.append(("\u03b3:",	2))

		if sdum=="PROD-CUT":
			layout.append(("Fudgem:", 3))
			layout.extend(__mat456_range)
		else:
			layout.append(_None2)
			layout.extend(__reg456_range)
		return layout

#------------------------------- EMFFIX ----------------------------------------
_EMFFIX_layout4 = [_None2,
		("Mat1:",	(1, LayoutManager.materialList)),
		("Max Frac.1:",	2),
		("Print:",	(0, ["", "PRINT", "NOPRINT" ])),
	_None2,
		("Mat2:",	(3, LayoutManager.materialList)),
		("Max Frac.2:",	4),
		_None2,
	_None2,
		("Mat3:",	(5, LayoutManager.materialList)),
		("Max Frac.3:",	6),
		_None2 ]

#------------------------------- EMFFLUO ---------------------------------------
_EMFFLUO_layout4 = [ _None2,
		("Fluorescence:",	(1, _onoff3_what)),
		_None2,
		_None2,
	_None2,
		("Mat:",	(2, LayoutManager.materialList)),
		("to Mat:",	(3, LayoutManager.materialList)),
		("Step:",	4) ]

#------------------------------- EMFRAY ----------------------------------------
_EMFRAY_corr = [ _empty_what,
		_ignore_what,
		("Rayleigh+Compton(bind)",	1),
		("Rayleight",		2),
		("Compton(bind)",	3),
		("R&C (bind+prof)",	4),
		("Rayleight (prof)",	5),
		("Compton (bind+prof)",	6),
		("Nothing",		-1) ]
_EMFRAY_layout4 = [ _None2,
		("Activate:",	(1, _EMFRAY_corr)),
		_None2,
		_None2,
	_None2,
		("Reg:",	(2, LayoutManager.regionList)),
		("to Reg:",	(3, LayoutManager.regionList)),
		("Step:",	4) ]

#------------------------------- EVENTDAT --------------------------------------
_EVENTDAT_layout4 = [ _None2,
		("Unit:", (1, lambda l: l.unitList("EVENTDAT"))),
		("File:", 0),
		_None2]

#------------------------------- EVENTYPE --------------------------------------
__EVENTYPE_list = [ _empty_what,
		("Reset", -2),
		("dE/dx only", 1),
		("All recoils", 2),
		("Recoils pid<=d (3)", 3),
		("Recoils pid<=t (4)", 4),
		("Recoils pid<=3-He (5)", 5),
		("Recoils pid<=4-He (6)", 6) ]
__EVENTYPE_model = [ "", "EVAP", "DPMJET" ]

_EVENTYPE_layout4 = [ _None2,
		("Ions:",	(3, __EVENTYPE_list)),
		("Model:",	(0, __EVENTYPE_model)),
		_None2 ]

#------------------------------- EVXTEST ---------------------------------------
__EVXTEST_what5 = [ _empty_what,
		("Include cascade+excitation", 0),
		("Exclude cascade+excitation", 1),
		("Evaporation", -1),
		("gamma Deexcitation", -2)]
__EVXTEST_what6 = [ _empty_what,
		("Included", 0),
		("Excluded", 1) ]
__EVXTEST_sdum = ["PRINT", "" ]
_EVXTEST_layout4 = [ _None2,
		("Events:",	1),
		("Projectile:",	(2, Input.Particle.list)),
		("Momentum:",	3),
	("Mat:",		(4, LayoutManager.materialList)),
		("Cascade:",	(5, __EVXTEST_what5)),
		("Diffractive:",(6, __EVXTEST_what6)),
		("Print",	(0, __EVXTEST_sdum)) ]

#------------------------------- FIXED -----------------------------------------
_FIXED_layout4 = [_None2,
		_None2,
		_None2,
		_None2 ]

#------------------------------- FLUKAFIX --------------------------------------
_FLUKAFIX_layout4 = [_None2,
		("Ekin frac:",	1),
		_None2, _None2, _None2 ] + __mat456_range

#------------------------------- FREE ------------------------------------------
_FREE_layout4 = [_None2,
		_None2,
		_None2,
		_None2 ]

#------------------------------- IONFLUCT --------------------------------------
__IONFLUCT_sdum = [ "", "PRIM-ION" ]

__IONFLUCT_accuracy = [ _empty_what,
			("Minimal", 1), ("Low",     2),
			("High",    3), ("Maximum", 4),
			_reset_what ]

__IONFLUCT_model = [ _empty_what,
			("Model-1", 1), ("Model-2", 2),
			("Model-3", 3), ("Model-4", 4),
			_reset_what ]

__IONFLUCT_other_layout4 = [ _None2,
		("Type:",	(0, __IONFLUCT_sdum)),
		("hadron:",	(1, _onoff3_what)),
		("e+e-:",	(2, _onoff3_what)),
	("Accuracy:",	(3, __IONFLUCT_accuracy))] \
		+ __mat456_range

__IONFLUCT_primion_layout4 = [ _None2,
		("Type:",	(0, __IONFLUCT_sdum)),
		("1st Vion:",	1),
		("Prim Ion/cm:", 2),
		("Model:",	(3, __IONFLUCT_model))] \
		+ __mat456_range

def _IONFLUCT_layout4(card):
	if card.sdum()=="PRIM-ION":
		return __IONFLUCT_primion_layout4
	else:
		return __IONFLUCT_other_layout4

#------------------------------- IONTRANS --------------------------------------
__IONTRANS_what1 = [ _empty_what,
		_ignore_what,
		("No transport of ions", 1),
		("Approximate", -1),
		("Full transport", -2),
		"HEAVYION",
		"DEUTERON",
		"TRITON",
		"3-HELIUM",
		"4-HELIUM"
		]

_IONTRANS_layout4 = [ _None2,
	("Transport:",	(1, __IONTRANS_what1)),
	_None2,
	_None2 ]

#------------------------------- IRRPROFI --------------------------------------
_IRRPROFI_layout4 = [ _None2,
		("\u0394t:",	1),
		("p/s:",	2),
		_None2,
	_None2,
		("\u0394t:",	3),
		("p/s:",	4),
		_None2,
	_None2,
		("\u0394t:",	5),
		("p/s:",	6),
		_None2 ]

#------------------------------- GEOBEGIN --------------------------------------
def __GEOBEGIN_geofile(layout, card, val=None):
	inp = layout.flair.project.input
	if val is None:
		return inp.geoFile
	else:
		inp.geoFile = layout.flair.project.relativePath(val)
		# if geometry is empty then load the file
		try:
			empty = inp["GEOEND"][0].pos()-card.pos() <= 3
		except:
			return
		if not empty: return

		# delete everything between GEOBEGIN GEOEND
		pos = card.pos()
		pos += 1
		while True:
			if inp.cardlist[pos].tag == "END":
				inp.delCard(pos)
			else:
				break

		geoinp = Input.Input()
		# Copy the old setting... correct it before!
		geoinp.geoFormat = layout.flair.project.input.geoFormat
		# prepare to load file
		try:
			geoinp._openFile(val,"r")
			geoinp._parseGeometry(card)
			for card in geoinp.cardlist:
				#say(">>>",card)
				inp.addCard(card, pos)
				pos += 1
		except:
			pass
		inp.renumber()
		layout.flair.refresh("card")

# -----------------------------------------------------------------------------
def __GEOBEGIN_geoout(layout, card, val=None):
	if val is None:
		return layout.flair.project.input.geoOutFile
	else:
		layout.flair.project.input.geoOutFile = layout.flair.project.relativePath(val)

__GEOBEGIN_layout4 = [_None2,
		("Log:", (1, [_empty_what, ("Errors",0), ("Nothing",1)])),
		("Acc:",  2),
		("Opt:", (5, [_empty_what, ("Default",0),
			("Logical", 1),
			("L+Plane", 2),
			("L+P+Bbox", 3)])),
	_None2,
		("Geometry:", (3, lambda l: l.unitList("GEOBEGIN",False))),
		("Out:",      (4, lambda l: l.unitList("GEOBEGIN",False))),
		("Fmt:",      (0, ["", "COMBNAME", "COMBINAT", "FLUGG"]))]

__GEOBEGIN_title = [("Title:",	 -1), _None2, _None2, _None2 ]

def _GEOBEGIN_layout4(card):
	layout = __GEOBEGIN_layout4[:]
	if card.intWhat(3) not in (0,5) or card.intWhat(4) not in (0,11):
		if card.intWhat(3) not in (0,5):
			layout.append(("File:",
				(__GEOBEGIN_geofile, (("Geometry file","*.geo"),("All","*")))))
		else:
			layout.append(_None2)
		layout.append(_None2)

		if card.intWhat(4) not in (0,11):
			layout.append(("Output:",
				(__GEOBEGIN_geoout,  (("Geometry file","*.out"),("All","*")))))
		else:
			layout.append(_None2)
		layout.append(_None2)
	layout.extend(__GEOBEGIN_title)
	return layout

#------------------------------- GEOEND ----------------------------------------
__GEOEND_type = ["", "DEBUG"]

__GEOEND_layout_normal = [ _None2,
		("",		(0, __GEOEND_type)),
		_None2,
		_None2 ]

__GEOEND_layout_debug = [ _None2,
		("",		(0, __GEOEND_type)),
		_None2,
		_None2,
	_None2,
		("Xmin:",	 4),
		("Xmax:",	 1),
		("NX:",		 7),
	_None2,
		("Ymin:",	 5),
		("Ymax:",	 2),
		("NY:",		 8),
	_None2,
		("Zmin:",	 6),
		("Zmax:",	 3),
		("NZ:",		 9) ]

def _GEOEND_layout4(card):
	# depending on the value of what(1) there are various layouts
	if card.sdum() == "DEBUG":
		return __GEOEND_layout_debug
	else:
		return __GEOEND_layout_normal

#------------------------------- _GLOBAL ---------------------------------------
__GLOBAL_analogue = [ _empty_what,	("Input", 0),
		("Analogue", -1),	("Biased", 2) ]
__GLOBAL_dnear    = [ _empty_what,	_default_what,
		("Always use",          1),
		("Don't use on steps", -1),
		("Do not use at all",  -2) ]
__GLOBAL_input    = [ _empty_what,	_default_what,
		("Names", 1),	("Free", 2),
		("Names+Free", 3),	("Format", 4),	_reset_what]
__GLOBAL_geometry = [ _empty_what,	_default_what,
		("Free", 1),	("Fixed", -1) ]

_GLOBAL_layout4	= [ _None2,
		("Max #reg:",	 1),
		("Analogue:",	(2, __GLOBAL_analogue)),
		("DNear:",	(3, __GLOBAL_dnear)),
	_None2,
		("Input:",	(4, __GLOBAL_input)),
		("Geometry:",	(5, __GLOBAL_geometry)),
		_None2 ]

#------------------------------- GCR-SPE ---------------------------------------
__GCR_sdum = ["Spectra", "DIPCOORD", "NO-NORM", "NORMALIZ"]
__GCR_field = [ ("naive dipole",	0),
		("exact multipole",	1),
		("MAGFLD routine",	2),
		("dipole @ infinite",	3),
		("undefined",		4) ]

def __GCR_SPE_sdum(layout, card, val=None):
	if val is None:
		if card.sdum() in ("DIPCOORD","NO-NORM","NORMALIZ"):
			return card.sdum()
		else:
			return "Spectra"
	else:
		if val in ("DIPCOORD","NO-NORM","NORMALIZ"):
			layout.setWhat(card, 0, val)
		else:
			layout.setWhat(card, 0, "")

def __GCR_what1_i0(layout, card, val=None):
	w1 = card.intWhat(1)
	i1, i0 = divmod(w1, 1000)
	if val is None:
		return i0
	else:
		try: ival = int(val)
		except: ival = 0
		layout.setWhat(card, 1, i1*1000 + ival)

def __GCR_what1_i1(layout, card, val=None):
	w1 = card.intWhat(1)
	i1, i0 = divmod(w1, 1000)
	if val is None:
		return i1
	else:
		try: ival = int(val)
		except: ival = 0
		layout.setWhat(card, 1, ival*1000 + i0)

def __GCR_what3_j0(layout, card, val=None):
	w3 = card.intWhat(3)
	j1, j0 = divmod(w3, 1000)
	if val is None:
		return j0
	else:
		try: jval = int(val)
		except: jval = 0
		layout.setWhat(card, 3, j1*1000 + jval)

def __GCR_what6_k0(layout, card, val=None):
	w6 = card.intWhat(6)
	k1, k0 = divmod(w6, 10000)
	if val is None:
		return k0
	else:
		try: kval = int(val)
		except: kval = 0
		layout.setWhat(card, 6, k1*1000 + kval)

def __GCR_what6_k1(layout, card, val=None):
	w6 = card.intWhat(6)
	k1, k0 = divmod(w6, 10000)
	if val is None:
		return k1
	else:
		try: kval = int(val)
		except: kval = 0
		layout.setWhat(card, 6, kval*1000 + k0)

__GCR_common = [_None2, ("Type:", (__GCR_SPE_sdum,__GCR_sdum)), _None2, _None2 ]

__GCR_SPE_layout4_dipcoord = __GCR_common[:] + \
	[ _None2,
		("Rdip:",	1),
		("\u03B8dip:",	2),
		("\u03C6dip:",	3),
	_None2,
		("Date:", 6),
		("\u03B8north:", 4),
		("\u03C6north:", 5)]

__GCR_SPE_layout4_default = __GCR_common[:] + \
	[ _None2,
		("Field:",	(__GCR_what1_i0, __GCR_field), "Magnetic field type"),
		("Shells:",	__GCR_what1_i1, "Number of atmospheric shells"),
		("Radius:",	4),
	("Equatorial Field:",	5),
		("Dump shell:",	__GCR_what6_k0),
		("Unit:",	__GCR_what6_k1),
		("DateFile:",	0)]

def _GCR_SPE_layout4(card):
	sdum = card.sdum()
	if sdum == "DIPCOORD":
		return __GCR_SPE_layout4_dipcoord
	elif sdum in ("NO-NORM", "NORMALIZ"):
		return __GCR_common
	else:
		return __GCR_SPE_layout4_default

#------------------------------- HI-PROPE --------------------------------------
_HI_PROPE_layout4 = [ _None2,
		("Z:",		1),
		("A:",		2),
		("Isom:",	3) ]

#------------------------------- LAM-BIAS --------------------------------------
__LAM_BIAS_sdum = ["",
		"DCDRBIAS",
		"DCY-DIRE",
		"DECALL",
		"DECPRI",
		"GDECAY",
		"INEALL",
		"INEPRI",
		"N1HSCBS" ]

__LAM_BIAS_layout4_dcydire = [ _None2,
		("Type:",	(0, __LAM_BIAS_sdum)),
		("Distribution:", (5,
				[_empty_what,
				 ("Point-like  exp(1-cos\u03B8)", 0),
				 ("Cylindrical exp(1-cos\u03B8)", 1),
				 ("Point-like  1/(1+(1-cos\u03B8)^2)", 100),
				 ("Cylindrical 1/(1+(1-cos\u03B8)^2) ", 101) ])),
		_None2,
	_None2,
		("\u03bb decay:", 4),
		("Background:", 6),
		_None2,
	_None2,
		("cosX:",	1),
		("cosY:",	2),
		("cosZ:",	3) ]

__LAM_BIAS_layout4_dcdrbias = [ _None2,
		("Type:",	(0, __LAM_BIAS_sdum)),
		("Activate:",	(1, _onoff3_what)),
		_None2,
	_None2,
		("Part:",	(4, Input.Particle.list)),
		("to Part:",	(5, Input.Particle.list)),
		("Step:",	6) ]

__LAM_BIAS_layout4_gdecay = [ _None2,
		("Type:",	(0, __LAM_BIAS_sdum)),
		("\u00d7 <\u03bb>:",	1),
		("\u00d7 \u03bb inelastic:",	2),
	("Mat:",		(3, LayoutManager.materialList)),
		("Part:",	(4, Input.Particle.list)),
		("to Part:",	(5, Input.Particle.list)),
		("Step:",	6) ]

__LAM_BIAS_layout4_n1hscbs = [ _None2,
		("Type:",	(0, __LAM_BIAS_sdum)),
		("Activate:",	(1, _onoff3_what)),
		_None2,
		_None2 ]

__LAM_BIAS_layout4_other = [ _None2,
		("Type:",	(0, __LAM_BIAS_sdum)),
		("\u00d7 mean life:",	1),
		("\u00d7 \u03bb inelastic:",	2),
	("Mat:",		(3, LayoutManager.materialList)),
		("Part:",	(4, Input.Particle.list)),
		("to Part:",	(5, Input.Particle.list)),
		("Step:",	6) ]

def _LAM_BIAS_layout4(card):
	sdum = card.sdum()
	if   sdum == "DCY-DIRE":
		return __LAM_BIAS_layout4_dcydire
	elif sdum == "DCDRBIAS":
		return __LAM_BIAS_layout4_dcdrbias
	elif sdum == "N1HSCBS":
		return __LAM_BIAS_layout4_n1hscbs
	elif sdum == "GDECAY":
		return __LAM_BIAS_layout4_gdecay
	else:
		return __LAM_BIAS_layout4_other

#------------------------------- LATTICE ---------------------------------------
_LATTICE_layout4 = [ _None2,
		("Reg:",	(1, LayoutManager.regionList)),
		("to Reg:",	(2, LayoutManager.regionList)),
		("Step:",	 3),
	( "Trans:",		(0, LayoutManager.rotdefiList)),
		("Lat:",	 4),
		("to Lat:",	 5),
		("Step:",	 6) ]

#------------------------------- LATTSNGL --------------------------------------
_LATTSNGL_layout4 = [ (None,     4),
		("Reg:",	(1, LayoutManager.regionList)),
		("R1:",	        (2, LayoutManager.rotdefiList)),
		("R2:",	        (3, LayoutManager.rotdefiList)),
		("Lat:",	 4)]

#------------------------------- LOW-BIAS --------------------------------------
_LOW_BIAS_layout4 = [ _None2,
		("Ecut:",	(1, LayoutManager.neutronGroupList)),
		("NonAnalogue:",(2, LayoutManager.neutronGroupList)),
		("Survival:",	 3),
		_None2 ] + __reg456_range

#------------------------------- LOW-DOWN --------------------------------------
_LOW_DOWN_layout4 = [ _None2,
		("Limit:",	(1, LayoutManager.neutronGroupList)),
		("Bias:",	 2),
		_None2, _None2 ] + __reg456_range

#------------------------------- LOW-MAT ---------------------------------------
def __lowmatFunc(layout, card, val=None):
	groups = layout.neutronGroupsUsed()
	if val is None:
		# Find value
		w1 = card.what(1)
		w2 = card.intWhat(2)
		w3 = card.intWhat(3)
		w4 = card.intWhat(4)
		s  = card.sdum()

		i = 1
		for g in groups:
			for mat in Input._lowMaterials[g]:
				if w2 == mat.id1 and \
				   w3 == mat.id2 and \
				   w4 == mat.id3 and \
				   s  == mat.name:
					return i
				i += 1

		if w2==0 and w3==0 and w4==0 and s=="":
			i = 1
			for g in groups:
				for mat in Input._lowMaterials[g]:
					if w1 == mat.name:
						# Force setting default value if needed
						__lowmatFunc(layout, card, i)
						return i
					i += 1
		return 0
	else:
		# Set value
		try:
			val = int(val)
		except:
			return

		if val==0:
			undoinfo = [
				layout.flair.setWhatUndo(card, 2, ""),
				layout.flair.setWhatUndo(card, 3, ""),
				layout.flair.setWhatUndo(card, 4, ""),
				layout.flair.setWhatUndo(card, 0, "") ]
		else:
			val -= 1
			group = -1
			for g in groups:
				l = len(Input._lowMaterials[g])
				if val < l:
					group = g
					break
			mat = Input._lowMaterials[group][val]
			undoinfo = [
				layout.flair.setWhatUndo(card, 2, mat.id1),
				layout.flair.setWhatUndo(card, 3, mat.id2),
				layout.flair.setWhatUndo(card, 4, mat.id3),
				layout.flair.setWhatUndo(card, 0, mat.name) ]
		layout.flair.addUndo(undo.createListUndo(undoinfo))

_LOW_MAT_layout4 = [ _None2,
		("Mat:",	(1, LayoutManager.materialList)),
		("LowMat:",	(__lowmatFunc, LayoutManager.lowmatList),
			"Combined identifier WHAT(1,2,3) of the corresponding low-energy neutron material"),
		_None2,
		_None2 ]

#------------------------------- LOW-NEUT --------------------------------------
__LOW_NEUT_ng   = [ ("260", 260), ("72", 72) ]
__LOW_NEUT_gg   = [ ("42",   42), ("22", 22), ]
__LOW_NEUT_maxe = [ "0.02", "0.0196" ]

__LOW_NEUT_print = [ _empty_what,
		("Minimum",      0),
		("Integral XS",  1),
		("+Downscatter", 2),
		("+Angles",      3),
		("+Residual",    4)]

__LOW_NEUT_pt    = [ _empty_what,
		_default_what,
		("Enable",       1),
		(">3eV S(a,b)",  2),
		(">3eV all",     3),
		("Everywhere",   4),
		("Disable",     -1)
		]
__LOW_NEUT_fiss = [ _empty_what,
		_default_what,
		("Forced 1",  1),
		("Normal",   -1)
		]

def __LOW_NEUT_what6_pt(layout, card, val=None):
	w = card.intWhat(6)
	i1,i0 = divmod(abs(w),10)
	if w<0:
		i0 = -i0
		i1 = -i1
	if val is None:
		return i0
	else:
		try:
			i0 = int(val)
			if i0>=0:
				if i1<0: i1=0
				nw = i1*10 + i0
			else:
				nw = i0
			if nw != w:
				layout.setWhat(card, 6, nw)
		except ValueError:
			pass

def __LOW_NEUT_what6_fiss(layout, card, val=None):
	w = card.intWhat(6)
	i1,i0 = divmod(abs(w),10)
	if w<0:
		i0 = -i0
		i1 = -i1
	if val is None:
		return i1
	else:
		try:
			i1 = int(val)
			if i1>=0:
				if i0<0: i0=0
				nw = i1*10 + i0
			else:
				nw = i1*10
			if nw != w:
				layout.setWhat(card, 6, nw)
		except ValueError:
			pass

_LOW_NEUT_layout4 = [ _None2,
		("n-groups:",	   (1, __LOW_NEUT_ng)),
		("\u03b3-groups:",(2, __LOW_NEUT_gg)),
		("Emax:",	   (3, __LOW_NEUT_maxe)),
	_None2,
		("Print:",	   (4, __LOW_NEUT_print)),
		("Point Wise:",	   (__LOW_NEUT_what6_pt,   __LOW_NEUT_pt),
				"Point wise treatment"),
		("(n,f) Multiplicity:", (__LOW_NEUT_what6_fiss, __LOW_NEUT_fiss),
					"Neutron induced fission multiplicity")]

#------------------------------- LOW-PWXS --------------------------------------
#__LOW_PWXS_sdum = ["1-H PWXS","2-H PWXS",]
_LOW_PWXS_layout4 = [ _None2,
		("Isotope:",	1),
		("Grid:",	2),
		("Temp.:",	3),
		("Mat:",	(4, LayoutManager.materialList)),
		("to Mat:",	(5, LayoutManager.materialList)),
		      ("Step:",	 6),
                      ("Special:",  (0,["","1-H PWXS","2-H PWXS","3-HE PWXS","4-HE PWXS","6-LI PWXS","12-C PWXS"]))]
#------------------------------- MATERIAL --------------------------------------
_MATERIAL_layout4 = [ _None2,
		("Name:",	0),
		("#",		4),
		("\u03c1:",	3),
	("Z:",			1),
		("Am:",		2),
		("A:",		6),
		("dE/dx:",	(5, LayoutManager.materialList)) ]

#------------------------------- MAT-PROP --------------------------------------
if Input._developer:
	__MAT_PROP_sdum = [ "", "DPA-ENER", "LOWNTEMP", "USERDIRE", "X-REFLEC" ]
else:
	__MAT_PROP_sdum = [ "", "DPA-ENER", "USERDIRE", "X-REFLEC" ]
__MAT_PROP_call = [ _empty_what, ("USRMED", 1), _reset_what ]

__MAT_PROP_layout4_default = [ _None2,
		("Type:",	(0, __MAT_PROP_sdum)),
		("Gas pressure:", 1),
		("RHOR:",	2),
		("Ionization:",	3) ] + __mat456_range

__MAT_PROP_layout4_dpaener = [ _None2,
		("Type:",	(0, __MAT_PROP_sdum)),
		("DPA Eth:",    1),
		_None2, _None2 ] + __mat456_range

__MAT_PROP_layout4_lowntemp = [ _None2,
		("Type:",	(0, __MAT_PROP_sdum)),
		("T_act/T_xs",	1),
		("Thermal groups:",	2),
		("Tnominal:",	3) ] + __mat456_range

__MAT_PROP_layout4_userdire = [ _None2,
		("Type:",	(0, __MAT_PROP_sdum)),
		("Call:",	(1, __MAT_PROP_call)),
		_None2, _None2 ] + __mat456_range

__MAT_PROP_layout4_xreflec = [ _None2,
		("Type:",	(0, __MAT_PROP_sdum)),
		("Data set:",	(1, [-1, 0, 998, 1000])),
		("x-ray coming from:",	(2,LayoutManager.materialList)),
		_None2, ] + __mat456_range

def _MAT_PROP_layout4(card):
	sdum = card.sdum()[:8]
	if sdum=="DPA-ENER":
		return __MAT_PROP_layout4_dpaener
	elif sdum=="LOWNTEMP":
		return __MAT_PROP_layout4_lowntemp
	elif sdum=="USERDIRE":
		return __MAT_PROP_layout4_userdire
	elif sdum=="X-REFLEC":
		return __MAT_PROP_layout4_xreflec
	else:
		return __MAT_PROP_layout4_default

#------------------------------- MCSTHRESh -------------------------------------
_MCSTHRES_layout4 = [_None2,
		("Primary:",   1),
		("Secondary:", 2),
		_None2 ]

#------------------------------- MGNFIELD --------------------------------------
_MGNFIELD_layout4 = [ _None2,
		("Max Ang (deg):",	1),
		("Bound Acc. (cm):",	2),
		("Min step (cm):",	3),
	_None2,
		("Bx:",		4),
		("By:",		5),
		("Bz:",		6)]

#------------------------------- MULSOPT ---------------------------------------
__MULSOPT_sdum = [ "", "FANO-ON", "FANO-OFF", "MLSH-ON", "MLSH-OFF",
		"GLOBAL", "GLOBEMF", "GLOBHAD" ]
__MULSOPT_strech = ["",   "1",  "2",  "3",  "4",
		        "11", "12", "13", "14",
		       "101","102","103","104",
		       "111","112","113","114"]

__MULSOPT_what1 = [ _empty_what, _on_what, _off_what, _reset_what ]

def __MULSOPT_opt(layout, card, val=None):
	w1 = card.intWhat(1)
	if w1<0:
		i0   = -1
		i2i1 =  0
	else:
		i2i1, i0 = divmod(w1, 10)
	if val is None:
		if card.what(1)=="":
			return ""
		else:
			return i0
	else:
		if val=="":
			layout.setWhat(card, 1, "")
		else:
			try:
				i0 = int(val)
				if i0<0:
					w1 = -1
				else:
					w1 = i0 + i2i1*10
				layout.setWhat(card, 1, w1)
			except ValueError:
				pass

def __MULSOPT_steps_i1(layout, card, val=None):
	w1 = card.intWhat(1)
	w1, i0 = divmod(w1, 10)
	i2, i1 = divmod(w1, 10000)
	if val is None:
		return i1
	else:
		try:
			layout.setWhat(card, 1, i0 + int(val)*10 + i2*100000)
		except ValueError:
			pass

def __MULSOPT_steps_i2(layout, card, val=None):
	i2,i1i0 = divmod(card.intWhat(1), 100000)
	if val is None:
		return i2
	else:
		try:
			layout.setWhat(card, 1, i1i0 + int(val)*100000)
		except ValueError:
			pass

__MULSOPT_what2 = [ ("No corrections",	0),
		   ("1st Born",		1),
		   ("2nd Born",		2),
		   ("1st+Finite size", -1),
		   ("2nd+Finite size", -2),
		   ("no MCS",		3),
		   ("Form factors",    -3) ]

__MULSOPT_what3 = [ ("No corrections",	0),
		   ("1st Born",		1),
		   ("2nd Born",		2),
		   ("1st+Finite size", -1),
		   ("2nd+Finite size", -2),
		   ("no MCS",		3) ]

__MULSOPT_fano_off_layout4 = [ _None2,
		("Type:",      (0, __MULSOPT_sdum)),
		("Optimize:",  (__MULSOPT_opt, __MULSOPT_what1),
			"Optimization switch for single scattering"),
		_None2,

		_None2,
		_None2,
		("h/\u03bc Corr:", (2, __MULSOPT_what2)),
		("e-e+ Corr:",      (3, __MULSOPT_what3)),
		_None2 ] + __mat456_range

__MULSOPT_fano_on_layout4 = [ _None2,
		("Type:",      (0, __MULSOPT_sdum)),
		("Optimize:",  (__MULSOPT_opt, __MULSOPT_what1),
			"WHAT(1) I0: Optimization switch for single scattering"),
		_None2,
	("h/\u03bc steps:",	__MULSOPT_steps_i1,
				"WHAT(1) I1: Number of single scattering steps for hadrons and muons"),
		("e-e+ steps:", __MULSOPT_steps_i2,
				"WHAT(1) I2: Number of single scattering steps for electrons and positrons"),
		("h/\u03bc Corr:",	(2, __MULSOPT_what2)),
		("e-e+ Corr:",		(3, __MULSOPT_what3)),
		_None2 ] + __mat456_range

__MULSOPT_global_layout4 = [ _None2,
		("Type:",      (0, __MULSOPT_sdum)),
		("Min step:",   1),
		("Streching:", (2, __MULSOPT_strech)),
	("Optimal:",    3),
		("Single scat:", (4, _onoff3_what)),
		("E<Moliere:", (5, _onoff3_what)),
		("# scatterings:", 6)]

def _MULSOPT_layout4(card):
	if card.sdum() in __MULSOPT_sdum[5:]:
		return __MULSOPT_global_layout4
	else:
		if card.what(1)=="" or card.intWhat(1)<0:
			return __MULSOPT_fano_off_layout4
		else:
			return __MULSOPT_fano_on_layout4

#------------------------------- MUPHOTON --------------------------------------
__MUPHOTON_inter = [	_empty_what,
			_ignore_what,
			("Full", 1),
			("No secondaries", 2),
			_reset_what ]
if Input._developer:
	_MUPHOTON_layout4 = [ _None2,
			("\u03bc Inter", (1, __MUPHOTON_inter)),
			("\u03c3 long/trans", 2),
			("\u03c1 inter", 3),
			_None2 ] + __mat456_range
else:
	_MUPHOTON_layout4 = [ _None2,
			("\u03bc Inter", (1, __MUPHOTON_inter)),
			_None2,
			_None2,
			_None2] + __mat456_range

#------------------------------- OPEN ------------------------------------------
__OPEN_status  = [ "", "OLD", "NEW", "UNKNOWN", "READONLY", "SCRATCH" ]
__OPEN_layout4 = [ _None2,
		("Unit:",	(1, lambda l: l.unitList(None,True))),
		("Status:",	(0, __OPEN_status)),
		_None2,
	_None2,
		("File:",	(-1,("All","*"))),
		_None2,
		_None2 ]
__OPEN_scratch_layout4 = [ _None2,
		("Unit:",	(1, lambda l: l.unitList(None,True))),
		("Status:",	(0, __OPEN_status)),
		_None2 ]

def _OPEN_layout4(card):
	if card.sdum() == "SCRATCH":
		return __OPEN_scratch_layout4
	else:
		return __OPEN_layout4

#------------------------------- OPT-PROD --------------------------------------
__OPT_PROD_sdum = [ "", "CERE-OFF", "TRD-OFF", "SCIN-OFF",
			"CERENKOV", "CEREN-WV", "CEREN-OM", "TR-RADIA",
			"SCINTILL", "SCINT-WV", "SCINT-OM"]

__OPT_PROD_off_layout4 = [ _None2,
		("Type:",	(0, __OPT_PROD_sdum)),
		_None2, _None2, _None2 ] + __mat456_range

__OPT_PROD_cerenkov_layout4 = [ _None2,
		("Type:",	(0, __OPT_PROD_sdum)),
		("Emin:",	1),
		("Emax:",	2),
	_None2 ] + __mat456_range

__OPT_PROD_cerenwv_layout4 = [ _None2,
		("Type:",	(0, __OPT_PROD_sdum)),
		("\u03bbmin:",	1),
		("\u03bbmax:",	2),
	_None2 ] + __mat456_range

__OPT_PROD_cerenom_layout4 = [ _None2,
		("Type:",	(0, __OPT_PROD_sdum)),
		("\u03a9min:",	1),
		("\u03a9max:",	2),
	_None2 ] + __mat456_range

__OPT_PROD_scintill_layout4 = [ _None2,
		("Type:",	(0, __OPT_PROD_sdum)),
		("Ei:",		1),
		("fraction:",	2),
		("Time:",	3)] + __mat456_range

__OPT_PROD_scintwv_layout4 = [ _None2,
		("Type:",	(0, __OPT_PROD_sdum)),
		("\u03bbi:",	1),
		("fraction:",	2),
		("Time:",	3)] + __mat456_range

__OPT_PROD_scintom_layout4 = [ _None2,
		("Type:",	(0, __OPT_PROD_sdum)),
		("\u03a9i:",	1),
		("fraction:",	2),
		("Time:",	3)] + __mat456_range

def _OPT_PROD_layout4(card):
	s = card.sdum()
	if   s == "CERENKOV" or s == "TR-RADIA":
		return __OPT_PROD_cerenkov_layout4
	elif s == "CEREN-WV":
		return __OPT_PROD_cerenwv_layout4
	elif s == "CEREN-OM":
		return __OPT_PROD_cerenom_layout4
	elif s == "SCINTILL":
		return __OPT_PROD_scintill_layout4
	elif s == "SCINT-WV":
		return __OPT_PROD_scintwv_layout4
	elif s == "SCINT-OM":
		return __OPT_PROD_scintom_layout4
	else:
		return __OPT_PROD_off_layout4

#------------------------------- OPT-PROP --------------------------------------
__OPT_PROP_sdum = [ "", "WV-LIMIT", "OM-LIMIT", "METAL",
			"SENSITIV", "WV-SENSI", "OM-SENSI", "SPEC-BDX", "RESET"]

__OPT_PROP_reset_layout4 = [ _None2,
		("Type:", (0, __OPT_PROP_sdum)),
		_None2, _None2, _None2] + __mat456_range

__OPT_PROP_empty_layout4 = [ _None2,
		("Type:", (0, __OPT_PROP_sdum)),
		_None2,
		_None2,
	_None2,
		("Refraction:",	1),
		("Absorption:",	2),
		("Diffusion:",	3),
	_None2 ] + __mat456_range

__OPT_PROP_wvlimit_layout4 = [ _None2,
		("Type:", (0, __OPT_PROP_sdum)),
		_None2,
		_None2,
	_None2,
		("\u03bbmin:",	1),
		("\u03bb:",	2),
		("\u03bbmax:",	3),
	_None2 ] + __mat456_range

__OPT_PROP_omlimit_layout4 = [ _None2,
		("Type:", (0, __OPT_PROP_sdum)),
		_None2,
		_None2,
	_None2,
		("\u03a9min:",	1),
		("\u03a9:",	2),
		("\u03a9max:",	3),
	_None2 ] + __mat456_range

__OPT_PROP_metal_layout4 = [ _None2,
		("Type:", (0, __OPT_PROP_sdum)),
		_None2,
		_None2,
	_None2,
		("1st:",	1),
		("2nd:",	2),
		("1-Reflectivity:", 3),
	_None2 ] + __mat456_range

__OPT_PROP_sensitiv_layout4 = [ _None2,
		("Type:", (0, __OPT_PROP_sdum)),
		("0th:",	1),
		("1st:",	2),
	_None2,
		("2nd:",	3),
		("3rd:",	4),
		("Max sensitivity:", 5)]

__OPT_PROP_wvsensi_layout4 = [ _None2,
		("Type:", (0, __OPT_PROP_sdum)),
		_None2,
		_None2,
	_None2,
		("\u03bbmin:",	1),
		("\u03bb:",	2),
		("\u03bbmax:",	3)]

__OPT_PROP_omsensi_layout4 = [ _None2,
		("Type:", (0, __OPT_PROP_sdum)),
		_None2,
		_None2,
	_None2,
		("\u03a9min:",	1),
		("\u03a9:",	2),
		("\u03a9max:",	3) ]

__OPT_PROP_specbdx_layout4 = [ _None2,
		("Type:", (0, __OPT_PROP_sdum)),
		_None2, _None2,
	_None2,
		("Activate:",	(1, _onoff3_what)),
		("Reg:",	(2, LayoutManager.regionList)),
		("to Reg:",	(3, LayoutManager.regionList)),
	_None2,
		("Activate:",	(4, _onoff3_what)),
		("Reg:",	(5, LayoutManager.regionList)),
		("to Reg:",	(6, LayoutManager.regionList)) ]

def _OPT_PROP_layout4(card):
	s = card.sdum()
	if   s == "RESET":
		return __OPT_PROP_reset_layout4
	elif s == "WV-LIMIT":
		return __OPT_PROP_wvlimit_layout4
	elif s == "OM-LIMIT":
		return __OPT_PROP_omlimit_layout4
	elif s == "METAL":
		return __OPT_PROP_metal_layout4
	elif s == "SENSITIV":
		return __OPT_PROP_sensitiv_layout4
	elif s == "WV-SENSI":
		return __OPT_PROP_wvsensi_layout4
	elif s == "OM-SENSI":
		return __OPT_PROP_omsensi_layout4
	elif s == "SPEC-BDX":
		return __OPT_PROP_specbdx_layout4
	else:
		return __OPT_PROP_empty_layout4

#------------------------------- PAIRBREM --------------------------------------
__PAIRBREM_activate = [
		_empty_what,
		_ignore_what,
		("Pair Prod",		 1),
		("Bremsstrahlung",	 2),
		("Both",		 3),
		("No Pair",		-1),
		("No Bremsstrahlung",	-2),
		("Inhibit both",	-3) ]

_PAIRBREM_layout4 = [ _None2,
		("Act:",	(1, __PAIRBREM_activate)),
		("e-e+ Thr:",	 2),
		("\u03b3 Thr:", 3),
		_None2 ] + __mat456_range

#------------------------------- PART-THR --------------------------------------
__PART_THR_types = [
		("Momentum",  0),
		("Energy",    1),
		("\u03B3 gamma (E/mc^2)",   2),
		("\u03B7 eta (\u03b2*\u03B3)",   3) ]

def __PART_THR_type(layout, card, val=None):
	if val is None:
		s = 0
		if card.sign(1): s = 1
		if card.intWhat(5)>0: s += 2
		return s
	else:
		(w5,w1) = divmod(int(val),2)
		undoinfo = [
			layout.flair.setSignUndo(card, 1, bool(w1)),
			layout.flair.setWhatUndo(card, 5, w5) ]
		layout.flair.addUndo(undo.createListUndo(undoinfo))

def _PART_THR_layout4(card):
	layout = [ _None2,
		("Type:",	(LayoutManager._what1_type, __PART_THR_types),
				"WHAT(1) sign: Energy or Momentum") ]

	if card.intWhat(5)>0:
		labels = ( "\u03B3:", "\u03B7:" )
	else:
		labels = ( "p:", "E:" )
	if card.sign(1):
		label = labels[1]
	else:
		label = labels[0]

	layout += [ (label,	 LayoutManager._what1_val),
		_None2,
		_None2,
		("Part:",	(2, Input.Particle.list)),
		("to Part:",	(3, Input.Particle.list)),
		("Step:",	 4) ]

	return layout

#------------------------------- PHOTONUC --------------------------------------
__PHOTONUC_sdum = [ "", "ELECTNUC", "MUMUPAIR", "MUMUPRIM" ]

def __PHOTONUC_int_what1(layout, card, val=None, sub=1):
	ia, ih, id, iq, ig = 0, 0, 0, 0, 0
	w1 = card.intWhat(1)
	if   w1 == -1: ia = -1
	elif w1 ==  1: ia = 1
	elif w1 ==  2: ih = 1
	elif w1 ==  3: id = 1
	elif w1 ==  4: iq = 1
	elif w1 ==  5: ig = 1
	else:
		ia = 0
		ig,w1 = divmod(w1, 1000)
		iq,w1 = divmod(w1,  100)
		id,ih = divmod(w1,   10)

	if val is None:
		if   sub==1: return ia
		elif sub==2: return ih
		elif sub==3: return id
		elif sub==4: return iq
		else:        return ig
	else:
		try: val = int(val)
		except: val = 0
		if val == -1:
			layout.setWhat(card, 1, val)
			return

		if ia<0: ia = 0
		if   sub==1:
			layout.setWhat(card, 1, val)
			return
		elif sub==2: ih=val
		elif sub==3: id=val
		elif sub==4: iq=val
		else:        ig=val
		val = ih + id*10 + iq*100 + ig*1000
		if   val ==    1: val=2
		elif val ==   10: val=3
		elif val ==  100: val=4
		elif val == 1000: val=5
		layout.setWhat(card, 1, val)

__PHOTONUC_layout4_default = [ _None2,
		("Type:", (0, __PHOTONUC_sdum)),
		_None2,
		("All E:",	(lambda l,c,v=None,s=1:
				__PHOTONUC_int_what1(l,c,v,s), _onoffreset_what),
				"WHAT(1)=1.0: photonuclear interactions are activated at all energies"),
		("E>0.7GeV",	(lambda l,c,v=None,s=2:
				__PHOTONUC_int_what1(l,c,v,s), _onoffreset_what),
				"WHAT(1) ih=1: photonuclear interactions are activated only in the high energy range (above 0.7 GeV)"),
		("\u0394 resonance", (lambda l,c,v=None,s=3:
				__PHOTONUC_int_what1(l,c,v,s), _onoffreset_what),
				"WHAT(1) id=1: photonuclear interactions are activated only in the region of the Delta resonance"),
		("Quasi D",	(lambda l,c,v=None,s=4:
				__PHOTONUC_int_what1(l,c,v,s), _onoffreset_what),
				"WHAT(1) iq=1: only quasi-deuteron photonuclear interactions are activated"),
		("Giant Dipole", (lambda l,c,v=None,s=5:
				__PHOTONUC_int_what1(l,c,v,s), _onoffreset_what),
				"WHAT(1) ig=1: only Giant Dipole Resonance photonuclear interactions are activated"),
		_None2] + __mat456_range

def __PHOTONUC_mu_what1(layout, card, val=None, sub=1):
	w1 = card.intWhat(1)
	if w1 == -1:
		ich, iqe, iin, ids = -1, 0, 0, 0
	else:
		ids,w1  = divmod(w1, 1000)
		iin,w1  = divmod(w1,  100)
		iqe,ich = divmod(w1,   10)

	if val is None:
		if   sub==1: return ich
		elif sub==2: return iqe
		elif sub==3: return iin
		else:        return ids
	else:
		val = int(val)
		if val == -1:
			layout.setWhat(card, 1, val)
		else:
			if ich<0: ich = 0
			if   sub==1: ich = val
			elif sub==2: iqe = val
			elif sub==3: iin = val
			else:        ids = val
			layout.setWhat(card, 1, ich + iqe*10 + iin*100 + ids*1000)

__PHOTONUC_layout4_mumupair = [ _None2,
		("Type:", (0, __PHOTONUC_sdum)),
		#(u"\u03bc+\u03bc- prod:", (1, __PHOTONUC_mupair_prod)),
		_None2,
		_None2,
		("Coherent:",		(lambda l,c,v=None,s=1:
				__PHOTONUC_mu_what1(l,c,v,s), _onoffreset_what),
				"activate muon pair coherent production"),
		("InC. quasielastic:",	(lambda l,c,v=None,s=2:
				__PHOTONUC_mu_what1(l,c,v,s), _onoffreset_what),
				"activate muon pair incoherent quasielastic production"),
		("Inc. Inelastic:",	(lambda l,c,v=None,s=3:
				__PHOTONUC_mu_what1(l,c,v,s), _onoffreset_what),
				"activate muon pair incoherent inelastic production"),
		("Deep Inelastic:",	(lambda l,c,v=None,s=4:
				__PHOTONUC_mu_what1(l,c,v,s), _onoffreset_what),
				"activate muon pair deep inelastic production"),
		("Bias inter-\u03bb:", 2) ] + __mat456_range

def _PHOTONUC_layout4(card):
	sdum = card.sdum()
	if sdum in __PHOTONUC_sdum[2:]:
		return __PHOTONUC_layout4_mumupair
	else:
		return __PHOTONUC_layout4_default

#------------------------------- PHYSICS ---------------------------------------
__PHYSICS_sdum = ["CHARMDEC",
		  "COALESCE",
		  "DECAYS",
		  "DPMTHRES",
		  "EM-DISSO",
		  "EVAPORAT",
		  "INFLDCAY",
		  "IONBRPAI",
		  "IONSPLIT",
                  "ISOMERS",
		  "LIMITS",
		  "NEUTRINO",
		  "PEATHRES",
		  "QMDTHRES" ]
if Input._developer:
	__PHYSICS_sdum.extend([
		  "CAPTPROB",
		  "COAPARAM",
		  "DEDX-SET",
		  "RQMDFLAG",
		  "EXPNUCLV",
		  "FISSION",
		  "FL2005.6",
		  "FL2006.3",
		  "FL2008.3",
		  "FL2011.2",
		  "IASLEVEL",
		  "LEVEL-DE",
		  "PISIGMAS",
		  "PREEQUIL",
		  "PREFTCMP",
		  "QUASI-EL",
		  "RESO-REI",
		  "RESTODPM",
		  "TAUFOR"
		  ])
__PHYSICS_sdum.sort()

__PHYSICS_layout4_generic = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			("#1:", 1),
			("#2:", 2),
		("#3:", 3),
			("#4:", 4),
			("#5:", 5),
			("#6:", 6)]

__PHYSICS_coalesce_what1 = [_empty_what,
			_on_what,
			("+Deuteron", 1001),
			("off", -1),
			_ignore_what]

__PHYSICS_layout4_captprob = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("Law:",	(1, [_empty_what, _ignore_what, ("Z-Law", 1), ("Daniel's law", 2)])),
		("MNTHCP:",	2)]

__PHYSICS_layout4_charmdec = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("Charmed:",	(1, _onoff3_what)),
		_None2 ]

if Input._developer:
	__PHYSICS_layout4_coalesce = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			("Activate:",	(1,_onoff3_what)),
			_None2,
			("Model 1:",	2),
			("Model 2:",	3),
			("Model 3:",	5),
			("Model 4:",	6)]
else:
	__PHYSICS_layout4_coalesce = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			("Activate:",	(1,_onoff3_what)),
			_None2]

def __PHYSICS_decays_what1_flag1(layout, card, val=None):
	w = card.intWhat(1)
	f100,f1 = divmod(w,100)
	if val is None:
		if w<=-1:
			return -1
		else:
			return f1
	else:
		try: val = int(val)
		except: val = 0
		if w<0 or val<0:
			layout.setWhat(card, 1, val)
			return
		layout.setWhat(card, 1, f100*100 + val)

def __PHYSICS_decays_what1_flag100(layout, card, val=None):
	w = card.intWhat(1)
	f100,f1 = divmod(w,100)
	if val is None:
		if w<100:
			return 0
		else:
			return f100
	else:
		try: val = int(val)
		except: val = 0
		if f1<0: f1=0
		layout.setWhat(card, 1, val*100 + f1)

__PHYSICS_decays1 = [ _empty_what, _ignore_what,
		("Max Acc. with Polarization", 1),
		("Max Acc. no Polarization", 2),
		("Phase space like", 3),
		_reset_what]

__PHYSICS_decays100 = [ ("", 0),
		("Leptonic only", 1),
		("Hadronic only", 2)]

__PHYSICS_layout4_decays = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("Decays:",	(__PHYSICS_decays_what1_flag1, __PHYSICS_decays1), "Flag for particle decay"),
		("Allowed:",    (__PHYSICS_decays_what1_flag100, __PHYSICS_decays100), "Allowed decays"),
	_None2,
		("Part:",	(4, Input.Particle.list)),
		("to Part:",	(5, Input.Particle.list)),
		("Step:",	 6) ]

__PHYSICS_layout4_dedxset = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		_None2,
		_None2,
	("dE/dx Nuc:",		(1, _onoff3_what)),
		("Pred.Mat.",	(2, _onoff3_what)),
		("Barkas:",	(3, _onoff3_what)),
		("Bloch:",	(4, _onoff3_what))]

__PHYSICS_layout4_dpmthres = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("Emin h:",			1),
		("Emin ions DPM:",		2),
	("Emin Ions RQMD:",			3),
		("Smearing +/-\u0394E:",	4),
		("Smearing +/-\u0394E hA:",	5),
		("h-A primaries:",		(6, _onoff3_what))]

__PHYSICS_emdisso1 = [ _empty_what, _ignore_what,
		("No EM-Dissociation",	1),
		("Proj&Target EM-Disso",2),
		("Proj EM-Disso",	3),
		("Target EM-Disso",	4),
		_reset_what]

__PHYSICS_emdisso2 = [ _empty_what, _ignore_what, _on_what, ("off", 2), _reset_what]

if Input._developer:
	__PHYSICS_layout4_emdisso = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			("EM Disso:",	(1, __PHYSICS_emdisso1)),
			_None2,
		_None2,
			("Muon+/-:",    (2, __PHYSICS_emdisso2)),
			("Deuteron:",   (3, __PHYSICS_emdisso2)),
			("Pre-Tabulation:",   (3, __PHYSICS_emdisso2))]
else:
	__PHYSICS_layout4_emdisso = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			("EM Disso:",	(1, __PHYSICS_emdisso1)),
			_None2 ]

__PHYSICS_layout4_peathres = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("Nucleons:",	1),
		("Pions:",	2),

		("Kaons:",	3),
		("Kaonbars:",	4),
		("AntiNucleon:",5),
		("(Anti)Hyperons:",6)]

if Input._developer:
	__PHYSICS_layout4_ionbrpai = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			("Direct Pair:",(1, _onoff3_what)),
			("Bremss.:",	(2, _onoff3_what)),
		_None2,
			("Form fact:",	(3, _onoff3_what)),
			_None2,
			_None2 ]
else:
	__PHYSICS_layout4_ionbrpai = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			("Direct Pair:",(1, _onoff3_what)),
			_None2,
		_None2,
			("Form fact:",	(3, _onoff3_what)),
			_None2,
			_None2 ]

def __PHYSICS_evaporat_what1(card):
	w = card.intWhat(1)
	if w<0:
		return 0, 0, -1
	AZ,i0 = divmod(w, 100)
	Amax, Zmax = divmod(AZ, 100)
	return Amax, Zmax, i0

def __PHYSICS_evaporat_what1_i0(layout, card, val=None):
	Amax, Zmax,i0 = __PHYSICS_evaporat_what1(card)
	if val is None:
		return i0
	else:
		try: val = int(val)
		except: val = 0
		if val<0:
			layout.setWhat(card, 1, -1)
		else:
			layout.setWhat(card, 1, Amax*10000 + Zmax*100 + val)

def __PHYSICS_evaporat_what1_Zmax(layout, card, val=None):
	Amax, Zmax,i0 = __PHYSICS_evaporat_what1(card)
	if val is None:
		return Zmax
	else:
		try: val = int(val)
		except: val = 0
		layout.setWhat(card, 1, Amax*10000 + val*100 + i0)

def __PHYSICS_evaporat_what1_Amax(layout, card, val=None):
	Amax, Zmax,i0 = __PHYSICS_evaporat_what1(card)
	if val is None:
		return Amax
	else:
		try: val = int(val)
		except: val = 0
		layout.setWhat(card, 1, val*10000 + Zmax*100 + i0)

__PHYSICS_evaporat = [ _ignore_what,
			("Old Evaporation",          1),
			("New Evap no heavy frag",   2),
			("New Evap with heavy frag", 3),
			_reset_what ]
if Input._developer:
	__PHYSICS_layout4_evaporat = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			("Model:",	(__PHYSICS_evaporat_what1_i0, __PHYSICS_evaporat), "FLUKA evaporation model"),
			_None2,
		_None2,
			("Zmax:",       __PHYSICS_evaporat_what1_Zmax, "Maximum fragment Z"),
			("Amax:",       __PHYSICS_evaporat_what1_Amax, "Maximum fragment A"),
			("#2:", 2),
		("#3:", 3),
			("#4:", 4),
			("#5:", 5),
			("#6:", 6)]
else:
	__PHYSICS_layout4_evaporat = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			("Model:",	(__PHYSICS_evaporat_what1_i0, __PHYSICS_evaporat), "FLUKA evaporation model"),
			_None2,
		_None2,
			("Zmax:",       __PHYSICS_evaporat_what1_Zmax, "Maximum fragment Z"),
			("Amax:",       __PHYSICS_evaporat_what1_Amax, "Maximum fragment A"),
			_None2 ]

__PHYSICS_layout4_ionsplit = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("Ion Split:",	(1, _onoff3_what)),
		("Splitting:",	(6, [("Threshold",0), ("Probabilistic",1), ("On nonelastic",2)])),
	("Emin:",	2),
		("Emax:",	3),
		("Amin:",	4),
		("Amax:",	5)]

__PHYSICS_layout4_isomers = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("Explicit isomers:",	(1, _onoff3_what)),
		_None2 ]

__PHYSICS_layout4_limits = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("Max.Pcms (pp):",	1),
		_None2 ]

__PHYSICS_neutrino_flag = [ ("ignored", 0),
			    ("NC", 1),
			    ("CC", 2),
			    ("NC+CC", 3),
			    ("No interaction", -1) ]

__PHYSICS_neutrino_flag2 = [_empty_what,
			_ignore_what,
			("Performed", 1),
			("Not performed", -1)]

__PHYSICS_layout4_neutrino = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("Quasielastic:",	(1,__PHYSICS_neutrino_flag)),
		("Resonant:",		(2,__PHYSICS_neutrino_flag)),
	_None2,
		("Deep Inelastic:",	(3,__PHYSICS_neutrino_flag)),
		("Charm production:",	(4,__PHYSICS_neutrino_flag)),
		("Forced Interactions:",(6,__PHYSICS_neutrino_flag2)),
		]

if Input._developer:

        __PHYSICS_layout4_quasiel = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("Explicit quasi-elastic:",	(1, _onoff3_what)),
		_None2 ]

        __PHYSICS_layout4_qmdthres = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			("Emin(QMD):",		1),
			("Emin(BME):",		2),
			("QMD-RQMD:",		3),
			("BME-QMD:",		4),
			("Emax(fussion):",	5),
			("BME-RQMD:",		6)]
else:
	__PHYSICS_layout4_qmdthres = [ _None2,
			("Type:",	(0, __PHYSICS_sdum)),
			_None2,
			("Emin(BME):",		2),
		_None2,
			_None2,
			("Emax(fusion):",	5),
			("BME-RQMD:",		6)]

	__PHYSICS_layout4_quasiel = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
			("#1:", 1),
			("#2:", 2),
		("#3:", 3),
			("#4:", 4),
			("#5:", 5),
			("#6:", 6)]

__PHYSICS_layout4_rqmdflag_what1 = [ ("Fast cascade", 1), _ignore_what, ("QMD mode", -1) ]
__PHYSICS_layout4_rqmdflag = [ _None2,
		("Type:",	(0, __PHYSICS_sdum)),
		("QMD flag:", (1,__PHYSICS_layout4_rqmdflag_what1)),
		("Preequilibrium step:", (2,_onoffreset_what))]

def _PHYSICS_layout4(card):
	sdum = card.sdum()[:8]
	if   sdum == "CAPTPROB":
		return __PHYSICS_layout4_captprob
	elif sdum == "CHARMDEC":
		return __PHYSICS_layout4_charmdec
	elif sdum == "COALESCE":
		return __PHYSICS_layout4_coalesce
	elif sdum == "DECAYS":
		return __PHYSICS_layout4_decays
	elif sdum == "DEDX-SET":
		return __PHYSICS_layout4_dedxset
	elif sdum == "DPMTHRES":
		return __PHYSICS_layout4_dpmthres
	elif sdum == "EM-DISSO":
		return __PHYSICS_layout4_emdisso
	elif sdum == "IONBRPAI":
		return __PHYSICS_layout4_ionbrpai
	elif sdum == "INFLDCAY":
		return __PHYSICS_layout4_generic
	elif sdum == "IONSPLIT":
		return __PHYSICS_layout4_ionsplit
	elif sdum == "ISOMERS":
		return __PHYSICS_layout4_isomers
	elif sdum == "EVAPORAT":
		return __PHYSICS_layout4_evaporat
	elif sdum == "LIMITS":
		return __PHYSICS_layout4_limits
	elif sdum == "PEATHRES":
		return __PHYSICS_layout4_peathres
	elif sdum == "NEUTRINO":
		return __PHYSICS_layout4_neutrino
	elif sdum == "QUASI-EL":
		return __PHYSICS_layout4_quasiel
	elif sdum == "QMDTHRES":
		return __PHYSICS_layout4_qmdthres
	elif sdum == "RQMDFLAG":
		return __PHYSICS_layout4_rqmdflag
	else:
		return __PHYSICS_layout4_generic

#------------------------------- PLOTGEOM --------------------------------------
_PLOTGEOM_layout4 = [ _None2,
		("Axes:",	(1, [_empty_what, ("Yes", 0), ("No", 1)])),
		("Boundaries:",	(2, [_empty_what, ("Yes", 0), ("No", 1)])),
		("Numbering:",	(3, [_empty_what, ("No", 0), ("Yes", 1)])),
		("Worm:",	 4),
		("Verbose:",	(5, [_empty_what, ("No", 0), ("Yes", 1)])),
		("Format:",	(0, [("Ascii","FORMAT"), ("Binary","")])),
		_None2,
	("Title:",	-1),
		_None2,
		_None2,
		_None2,
	_None2,
		("Xmin:", 7),
		("Ymin:", 8),
		("Zmin:", 9),
	_None2,
		("Xmax:", 10),
		("Ymax:", 11),
		("Zmax:", 12),
	_None2,
		("Yx:", 13),
		("Yy:", 14),
		("Yz:", 15),
	_None2,
		("Xx:", 16),
		("Xy:", 17),
		("Xz:", 18),
	("ExpY:",	19),
		("ExpX:", 20),
		("NX:", 21),
		("NY:", 22) ]

#------------------------------- POLARIZA --------------------------------------
__POLARIZA_flag = [	_empty_what,
			("Not polarized", 0),
			("Orthogonal",    1)
		]
_POLARIZA_layout4 = [ _None2,
		("cosx:",	1),
		("cosy:",	2),
		("cosz:",	3),
	_None2,
		("Relative:",	(4, _onoffempty_what)),
		("Fraction:",	5),
		("1-Fraction:", (6,__POLARIZA_flag))]

#------------------------------- RAD-BIOL --------------------------------------
_RAD_BIOL_layout4 = [ _None2,
		("File:",	0),
		("default",	1),
	_None2,
	_None2,
	_None2,
	_None2,
	               _None2 ]

#------------------------------- RADDECAY --------------------------------------
__RADDECAY_layout4_what1 = [ _empty_what, ("Active", 1),
		("Semi-Analogue", 2), _reset_what ]
__RADDECAY_layout4_digit = [ _ignore_what, ("Prompt", 1), ("Decay", 2), ("Both",3)]

def __RADDECAY_decay_cut(layout, card, val=None):
	wh,wl = divmod(card.longWhat(5), 100000)
	if val is None:
		return float(wh)
	else:
		try: val = bmath.limit(0, int(float(val)), 99999)
		except ValueError: val = 0
		layout.setWhat(card, 5, val*100000 + wl)

def __RADDECAY_prompt_cut(layout, card, val=None):
	wh,wl = divmod(card.longWhat(5), 100000)
	if val is None:
		return float(wl)
	else:
		try: val = bmath.limit(0, int(float(val)), 99999)
		except ValueError: val = 0
		layout.setWhat(card, 5, wh*100000 + val)

_RADDECAY_layout4 = [ _None2,
		("Decays:",	(1, __RADDECAY_layout4_what1)),
		("Patch Isom:",	(2, _onoff3_what)),
		("Replicas:",	3),
	# a
	("h/\u03bc Int:",	(lambda l,c,v=None:l._what_digit_val(c,4,0,v),
					__RADDECAY_layout4_digit),
			"WHAT(4) a: hadron/muon interaction length or decay biasing, as defined by command LAM-BIAS"),
	# b
	("h/\u03bc LPB:",	(lambda l,c,v=None:l._what_digit_val(c,4,1,v),
					__RADDECAY_layout4_digit),
			"WHAT(4) b: hadron/muon leading particle  biasing"),
	# c
	("h/\u03bc WW:",	(lambda l,c,v=None:l._what_digit_val(c,4,2,v),
					__RADDECAY_layout4_digit),
			"WHAT(4) c: hadron/muon importance and Weight Window biasing, as defined by commands BIASING and WW-FACTOr"),
	# d
	("e-e+ Int:",		(lambda l,c,v=None:l._what_digit_val(c,4,3,v),
					__RADDECAY_layout4_digit),
			"WHAT(4) d: e+/e-/gamma interaction length biasing, as defined by command EMF-BIAS"),
	# e
	("e-e+ LPB:",		(lambda l,c,v=None:l._what_digit_val(c,4,4,v),
					__RADDECAY_layout4_digit),
			"WHAT(4) e: e+/e-/gamma leading particle biasing, as defined by command EMF-BIAS"),
	# f
	("e-e+ WW:",		(lambda l,c,v=None:l._what_digit_val(c,4,5,v),
					__RADDECAY_layout4_digit),
			"WHAT(4) f: e+/e-/gamma importance and Weight Window biasing, as defined by commands BIASING and WW-FACTOr"),
	# g
	("Low-n Bias:",		(lambda l,c,v=None:l._what_digit_val(c,4,6,v),
					__RADDECAY_layout4_digit),
			"WHAT(4) g: low-energy neutron biased downscattering, as defined by command LOW-DOWN, and non-analogue absorption, as defined by LOW-BIAS"),
	# h...
	# i
	("Low-n WW:",		(lambda l,c,v=None:l._what_digit_val(c,4,8,v),
					__RADDECAY_layout4_digit),
			"WHAT(4) i: low-energy neutron importance and Weight Window biasing, as defined by commands BIASING and WW-FACTOr"),
		_None2,
		("decay cut:",	__RADDECAY_decay_cut,
			"WHAT(5) xxxxx*0.1: transport energy cut-off multiplication factor for beta+, beta- and gamma decay radiation"),
		("prompt cut:",	__RADDECAY_prompt_cut,
			"WHAT(5) yyyyy*0.1: transport energy cut-off multiplication factor for prompt e+, e- and gamma radiation"),
		("Coulomb corr:",	(6, _onoff3_what))]

#------------------------------- RANDOMIZE -------------------------------------
_RANDOMIZ_layout4 = [ _None2,
		("Unit", (1, lambda l: l.unitList("RANDOMIZ",False))),
		("Seed:", 2),
		_None2 ]

#------------------------------- RESNUCLEi -------------------------------------
_RESNUCLE_type = [("Spallation",  1),
		   ("Low Energy",  2),
		   ("All",  3)]

_RESNUCLE_layout4 = [ _None2,
		("Type:",	(1, _RESNUCLE_type)),
		("Unit:",	(2, lambda l: l.unitList("RESNUCLE"))),
		("Name:",	 0),
	("Max Z:",		 3),
		("Max M:",	 4),
		("Reg:",	(5, LayoutManager.allRegionList)),
		("Vol:",	 6) ]

#------------------------------- ROT-DEFI --------------------------------------
__ROT_DEFI_axis = [ ("X", 1), ("Y", 2), ("Z", 3), ("Z", 0) ]
def __ROT_DEFI_id(layout, card, val=None):
	w1 = card.intWhat(1)
	if w1 < 1000:
		(j, i) = divmod(w1, 100)
	else:
		(i, j) = divmod(w1, 1000)
	if val is None:
		return i
	else:
		try: ival = int(val)
		except: ival = 0
		if ival<100 and w1<1000:
			layout.setWhat(card, 1, j*100 + ival)
		else:
			layout.setWhat(card, 1, ival*1000 + j)

def __ROT_DEFI_rot(layout, card, val=None):
	w1 = card.intWhat(1)
	if w1 < 1000:
		(j, i) = divmod(w1, 100)
	else:
		(i, j) = divmod(w1, 1000)
	if val is None:
		return j
	else:
		try: ival = int(val)
		except: ival = 0
		if w1<1000:
			layout.setWhat(card, 1, ival*100 + i)
		else:
			layout.setWhat(card, 1, i*1000 + ival)

_ROT_DEFI_layout4 = [ _None2,
		("Axis:",	(__ROT_DEFI_rot, __ROT_DEFI_axis),
				"WHAT(1) j: Transformation axis"),
		("Id:",		__ROT_DEFI_id, "WHAT(1) i: Transformation index"),
		("Name:",	0),
	_None2,
		("Polar:",	2),
		("Azm:",	3),
		_None2,
	_None2,
		("\u0394x:",	4),
		("\u0394y:",	5),
		("\u0394z:",	6) ]

#------------------------------- ROTPRBIN --------------------------------------
_ROTPRBIN_layout4 = [ _None2,
		("Storage:",	 1),
		_None2,
		("Rot:",	(2, LayoutManager.rotdefiList2)),
	_None2,
		("Bin:",	(4, LayoutManager.binningList)),
		("to Bin:",	(5, LayoutManager.binningList)),
		("Step:",	 6) ]

#------------------------------- SCORE -----------------------------------------
_SCORE_layout4	= [
	_None2,
		("Part1:", (1, Input.Particle.listAll)),
		("Part2:", (2, Input.Particle.listAll)),
		_None2,
	_None2,
		("Part3:", (3, Input.Particle.listAll)),
		("Part4:", (4, Input.Particle.listAll)),
		_None2 ]

#------------------------------- START -----------------------------------------
_START_layout4 = [ _None2,
		("No.:", 1),
		("Core:", (4,[_empty_what, _on_what, _off_what])),
		_None2,
	_None2,
		("Time:", 6),
		("Report:", (5, [_default_what, ("Every history",1)])),
		_None2]

#------------------------------- STEPSIZE --------------------------------------
_STEPSIZE_layout4 = [ _None2,
		("Min (cm):",	1),
		("Max (cm):",	2),
		_None2,
	_None2,
		("Reg:",	(3, LayoutManager.regionList)),
		("to Reg:",	(4, LayoutManager.regionList)),
		("Step:",	5) ]

#------------------------------- STERNHEIme ------------------------------------
_STERNHEI_layout4 = [ _None2,
		("Cbar:",	1),
		("X0:",		2),
		("X1:",		3),
	("Mat:",		(0, LayoutManager.materialList)),
		("a:",		4),
		("m:",		5),
		("\u03b40:",	6) ]

#------------------------------- SPECSOUR --------------------------------------
__SPECSOUR_types = [	"",
			"BEAMSPOT",
			"BIN-SOUR",
			"CROSSASY",
			"CROSSSYM",
			"GCR-ALLF",
			"GCR-AMS",
			"GCR-IONF",
			"GCR-SPEC",
			"PPSOURCE",
			"SPE-2003",
			"SPE-2005",
			"SPE-SPEC",
			"SYNC-RAD",
			"SYNC-RDN",
		]

__SPECSOUR_line1 = [ _None2,
		("Type:",	(0, __SPECSOUR_types)),
		_None2,
		_None2 ]

__SPECSOUR_lastlines = [
	_None2,
		("\u03c3x:",	7),
		("\u03c3y:",	8),
		("\u03c3z:",	9),
	("limit-\u03c3:",	10),
		("Part:",	(11,Input.Particle.hadron)),
		("A:",		12),
		("Z:",		13),
	("\u03c3\u03B8C_1:",	14),
		("\u03c30_1:",	15),
		("\u03c3\u03B8C_2:",	16),
		("\u03c30_2:",	17),
	_None2,
		("NonElastic:", (lambda l,c,v=None:l._what_digit_val(c,18,2,v,3), _onoff_what),
			"WHAT(18) : flag for nuclear nonelastic interactions"),
		("Elastic:",	(lambda l,c,v=None:l._what_digit_val(c,18,1,v,3), _onoff_what),
			"WHAT(18) : flag for nuclear elastic interactions"),
		("EM dissociation:",(lambda l,c,v=None:l._what_digit_val(c,18,0,v,3), _onoff_what),
			"WHAT(18) : flag for electromagnetic dissociation interactions")]

__SPECSOUR_ppsource_layout4 = __SPECSOUR_line1[:] + \
	[_None2,
		("P1x:",	1),
		("P1y:",	2),
		("P1z:",	3),
	_None2,
		("P2x:",	4),
		("P2y:",	5),
		("P2z:",	6) ] + __SPECSOUR_lastlines

__SPECSOUR_crossasy_layout4 = __SPECSOUR_line1[:] + \
	[ _None2,
		("P1lab:",	1),
		("Polar1:",	2),
		("Azimuthal:",	3),
	_None2,
		("P2-lab:",	4),
		("Polar2:",	5),
		_None2 ] + __SPECSOUR_lastlines

__SPECSOUR_crosssym_layout4 = __SPECSOUR_line1[:] + \
	[ _None2,
		("Plab:",	1),
		("CrossAng/2:",	2),
		("Azimuthal:",	3)] + __SPECSOUR_lastlines

__SPECSOUR_other_layout4 = __SPECSOUR_line1[:] + \
		[ ("#%d"%(i), i) for i in range(1,19) ] + \
		[_None2, _None2]

__SPECSOUR_type0 = [ ("Radius",0), ("MagField",1), ("Function",2) ]

__SPECSOUR_beamspot_layout4 =  __SPECSOUR_line1[:2] + [("# spots:",    1), _None2 ]

__SPECSOUR_binsour_layout4 =  __SPECSOUR_line1[:2] + \
		[("Unit:",    (1, lambda l: l.unitList("OPEN", False, True))),
		 ("Det Id:",   2),

		 _None2,
		("Rot-before:",	(3, LayoutManager.rotdefiList2)),
		("Rot-after:",	(4, LayoutManager.rotdefiList2)),
		 _None2 ]

def _SPECSOUR_layout4(card):
	# depending on the value of sdum
	sdum = card.sdum().upper()

	if sdum == "PPSOURCE":
		return __SPECSOUR_ppsource_layout4

	elif sdum == "CROSSASY":
		return __SPECSOUR_crossasy_layout4

	elif sdum == "CROSSSYM":
		return __SPECSOUR_crosssym_layout4

	elif sdum == "BIN-SOUR":
		return __SPECSOUR_binsour_layout4

	elif sdum == "BEAMSPOT":
		return __SPECSOUR_beamspot_layout4

	elif sdum in ("SYNC-RAD", "SYNC-RDN"):
		# First line
		layout = [ _None2,
			("Type:",	(0, __SPECSOUR_types)),
			("Part:",       (LayoutManager._what1_papro_what, Input.Particle.particles),
					"Particle type or heavy ion"),
			_None2]

		# Heavy ion selection
		if LayoutManager._papro_isotope(card, 1):
			layout.extend([_None2,
				    ("Z:",	LayoutManager._what1_papro_Z, "Heavy ion charge"),
				    ("A:",	LayoutManager._what1_papro_A, "Heavy ion mass"),
				    ("Isomer:",	LayoutManager._what1_papro_M, "Heavy ion isomer")])

		# what(2) Energy or momentum
		layout.append(("E/p:", (LayoutManager._what2_type, __BEAM_type0),
					"WHAT(2) sign: Energy or Momentum"))
		if card.sign(2):
			layout.append(("E:",	 LayoutManager._what2_val,
					"WHAT(2): Energy [GeV]"))
		else:
			layout.append(("p:",	 LayoutManager._what2_val,
					"WHAT(2): Momentum [GeV/c]"))

		# what(3) Radius or Magnetic field
		layout.append(("R/B:", (LayoutManager._what3_type, __SPECSOUR_type0),
					"WHAT(3) sign: Radius or Magnetic Field"))
		if card.sign(3):
			layout.append(("B:",	 LayoutManager._what3_val,
					"WHAT(3): Magnetic Field [T]"))
		else:
			layout.append(("R:",	 LayoutManager._what3_val,
					"WHAT(3): Radius [cm]"))

		# what(4) = minimum photon energy for sampling (GeV)
		# what(5) = x-component of the magnetic field versor
		# what(6) = y-component of the magnetic field versor
		# what(7) = length (cm) of the emission segment (arc)
		layout.extend([
			("E\u03B3min:", 4),
			("Bx:",     5),
			("By:",     6),
			("Length:", 7),
			_None2,
			("x2:",     8),
			("y2:",     9),
			("z2:",    10),
			_None2,
			("cosx2:", 11),
			("cosy2:", 12),
			_None2
			])
		return layout
	else:
		return __SPECSOUR_other_layout4

#------------------------------- SPOTBEAM --------------------------------------
def __SPOTBEAM_div_type(layout, card, val=None):
	return __BEAM_div_type(layout, card, val, w=5)

__SPOTBEAM_type = [("PRINT", -99999), ("BEAM", 0)]

def _SPOTBEAM_layout4_beam(card):
	# 1st line
	layout = [ ("Spot n.",	    0),	# spot number
		("Beam:", (LayoutManager._what2_type, __BEAM_type0),
			"WHAT(2) sign: Energy or Momentum") ]

	if card.sign(2):
		layout.append(("E:",	 LayoutManager._what2_val,
				"WHAT(2): Energy [GeV or GeV/n for HI]"))
	else:
		layout.append(("p:",	 LayoutManager._what2_val,
				"WHAT(2): Momentum [GeV/c]"))

	layout.append(("Part:",       (LayoutManager._paprospot_what, Input.Particle.particles),
			"Particle type or heavy ion"))
	if card.what(1)== "":
		layout.append(("Print",(1, __SPOTBEAM_type),"what(1)=-99999:print"))
# 1bis line
	if LayoutManager._paprospot_isotope(card):
		layout.extend([("Z:",	LayoutManager._paprospot_Z, "Heavy ion charge"),
			    ("A:",	LayoutManager._paprospot_A, "Heavy ion mass"),
			    ("Isomer:",	LayoutManager._paprospot_M, "Heavy ion isomer"),
			_None2])

	# 2nd line
	layout.append(("\u0394p:",	(LayoutManager._what4_type, __BEAM_type1),
			"WHAT(4): Momentum spread type"))
	if card.sign(2):
		layout.append(("\u0394p(FWHM):", LayoutManager._what4_val,
			"WHAT(4): Momentum spread Gaussian [GeV/c] FWHM=2.355*sigma"))
	else:
		layout.append(("\u0394p:",	LayoutManager._what4_val,
			"WHAT(4): Momentum spread FLAT [GeV/c]"))

	layout.append(("X \u0394\u03d5:",	(__SPOTBEAM_div_type, __BEAM_type2),
			"WHAT(5): x Divergence type"))
	if card.numWhat(5) >= 2000.0*math.pi:
		layout.append(_None2)
	else:
		layout.append(("x \u0394\u03d5:", LayoutManager._what5_val,
				"WHAT(5): x Divergence [mrad]"))

	# 3rd line
	layout.extend([_None2,
			("Spot Weight:", 3),
                       ("y \u0394\u03d5:",	(__SPOTBEAM_div_type, __BEAM_type2),
			"WHAT(6): y Divergence type"),
                       ("y \u0394\u03d5:", LayoutManager._what6_val,
				"WHAT(6): y Divergence [mrad]"),
#			("Printout:",   (6, [_empty_what, _on_what, _default_what, _reset_what])),
			_None2])
	return layout

def _SPOTBEAM_layout4_print (card):
	layout =[ _None2,
			("print?",	(1, __SPOTBEAM_type)),
			     ("Printout:",   (6, [_empty_what, _on_what, _default_what, _reset_what])),
			_None2]
	return layout
def _SPOTBEAM_layout4(card):
	what1 = card.what(1);
	if   what1 == "":
#		return __SPOTBEAM_layout4_generic
		return _SPOTBEAM_layout4_beam (card)
	if   card.numWhat(1) < -90000:
		return _SPOTBEAM_layout4_print(card)
	else:
		return _SPOTBEAM_layout4_beam(card)


#------------------------------- STOPPOS  --------------------------------------
__SPOTPOS_layout4 = [
		(None,	    0),	# spot number
		("x:",      1),
		("y:",      2),
		("z:",      3),
	]

def _SPOTPOS_layout4(card):
	# 1st line
	layout = __SPOTPOS_layout4[:]
	# 2nd line (use the same from BEAM card)
	__BEAM_shape_layout(card, layout)
	return layout

#------------------------------- SPOTDIR  --------------------------------------
_SPOTDIR_layout4 = [
		(None,	    0),	# spot number
		("cosx:",   1),
		("cosy:",   2),
		("cosz:",   3),
	_None2,
		("cosBxx:", 4),
		("cosBxy:", 5),
		("cosBxz:", 6),
	]

#------------------------------- SPOTTRAN --------------------------------------
_SPOTTRAN_layout4 = [
		("New Rot:",   (1, LayoutManager.rotdefiList2)),
		("Old Rot:",   (2, LayoutManager.rotdefiList2)),
	_None2,
		("From spot:", 4),
		("to spot:", 5),
		("Step:", 6),
		("Transf. flag:",       (0, ["","DIRECT","INVERSE"])) ]

#------------------------------- SPOTTIME --------------------------------------
_SPOTTIME_layout4 = [ _None2,
		("Start time (s):",   1),
		("Stop time (s):",   2),
	_None2,
	_None2,
		("From spot:", 4),
		("to spot:", 5),
		("Step:", 6)]

#------------------------------- STOP  -----------------------------------------
_STOP_layout4 = [ _None2, _None2, _None2, _None2 ]

#------------------------------- TCQUENCH --------------------------------------
_TCQUENCH_layout4 = [ _None2,
		("t cut-off:",	1),
		("Birks c1:",	2),
		("Birks c2:",	3),
	_None2,
		("Bin:",	(4, LayoutManager.binningList)),
		("to Bin:",	(5, LayoutManager.binningList)),
		("Step:",	6) ]

#------------------------------- TIME-CUT --------------------------------------
_TIME_CUT_layout4 = [ _None2,
		("Tcut:",	1),
		("Tdelay:",	2),
		("Start Mat:",	(3, LayoutManager.materialList)),
	_None2,
		("Part:",	(4, Input.Particle.list)),
		("to Part:",	(5, Input.Particle.list)),
		("Step:",	 6) ]

#------------------------------- TITLE -----------------------------------------
_TITLE_layout4 = [ (None,	-1),
		_None2,
		_None2,
		_None2]

#------------------------------- TSL-PWXS --------------------------------------
_TSL_PWXS_layout4 = [ _None2,
		("Isotope:",	1),
		("Bound in",	2),
		("Temp.:",	3),
		("Mat:",	(4, LayoutManager.materialList)),
		("to Mat:",	(5, LayoutManager.materialList)),
                ("Step:",	 6),
                 ("Release",  0)    ]
#------------------------------- THRESHOLD -------------------------------------
_THRESHOL_layout4 = [ _None2,
		("Elastic:",	3),
		("Inelastic:",	4),
		("Star:",	6) ]

#------------------------------- TPSSCORE --------------------------------------
__TPSSCORE_sdum = ["BIOTOBIN"]
__TPSSCORE_layout4_generic = [ _None2,
			("Type:",	(0, __TPSSCORE_sdum)),
			("#1:", 1),
			("#2:", 2),
		("#3:", 3),
			("#4:", 4),
			("#5:", 5),
			("#6:", 6)]

__TPSSCORE_biotobin_what2 = [_ignore_what]
__TPSSCORE_biotobin_what3 = [_ignore_what]
def _TPSSCORE_layout4_biotobin (card): 
	return [ _None2,
		("Type:",	(0, __TPSSCORE_sdum)),
		("Bio_set:",	1),
		_None2,
		_None2,
		("Det:",	(4, lambda l,c=card:l.detectorList("USRBIN"))),
		("to Det:",	(5, lambda l,c=card:l.detectorList("USRBIN"))),
		("Step:",	6)]

def _TPSSCORE_layout4(card):
	sdum = card.sdum()[:8]
	if   sdum == "BIOTOBIN":
		return _TPSSCORE_layout4_biotobin(card)
	else:
		return __TPSSCORE_layout4_generic

#------------------------------- USERWEIGHT ------------------------------------
__USERWEIG_what3 = [	_empty_what,
			_ignore_what,
			("No weight", -1),
			("FLUSCW", 1),
			("FLUSCW,FLDSCP", 2),
			("FLUSCW+",  3),
			("FLUSCW,FLDSCP+",  4) ]
__USERWEIG_what6 = [	_empty_what,
			_ignore_what,
			("No weight", -1),
			("COMSCW", 1),
			("COMSCW,ENDSCP", 2),
			("COMSCW+", 3),
			("COMSCW, ENDSCP+", 4) ]

_USERWEIG_layout4 = [ _None2,
		("\u03d5 Weight:",	(3, __USERWEIG_what3)),
		("Detect:",	(4, [("No weight",0,),("USRRNC",1)])),
                      _None2,
                      _None2,
		("Resnuclei:",	(5, [("No weight",0,),("DETSCW",1)])),
		      ("Energy:",	(6, __USERWEIG_what6)),
                     _None2]

#------------------------------- USRBIN ----------------------------------------
_USRBIN_type = [ _empty_what,
		("X-Y-Z"             , 10),
		("R-\u03A6-Z"       , 11),
		("Region"            , 12),
		("|X|-Y-Z"           , 13),
		("X-|Y|-Z"           , 14),
		("X-Y-|Z|"           , 15),
		("|X|-|Y|-|Z|"       , 16),
		("R-\u03A6-|Z|"     , 17),
		("Special"           , 18),
		("X-Y-Z point"       ,  0),
		("R-\u03A6-Z point" ,  1),
		("Region point"      ,  2),
		("|X|-Y-Z point"     ,  3),
		("X-|Y|-Z point"     ,  4),
		("X-Y-|Z| point"     ,  5),
		("|X|-|Y|-|Z| point" ,  6),
		("R-\u03A6-|Z| point", 7),
		("Special point"      , 8) ]

__USRBIN_layout4_XYZ = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("USRBIN"))),
		("Name:",	0),
	("Type:",		(1, _USRBIN_type)),
		("Xmin:",	 7),
		("Xmax:",	 4),
		("NX:",		10),
	("Part:",		(2, Input.Particle.listAll)),
		("Ymin:",	 8),
		("Ymax:",	 5),
		("NY:",		11),
	_None2,
		("Zmin:",	 9),
		("Zmax:",	 6),
		("NZ:",		12) ]

__USRBIN_layout4_RZPhi = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("USRBIN"))),
		("Name:",	0),
	("Type:",		(1, _USRBIN_type)),
		("Rmin:",	 7),
		("Rmax:",	 4),
		("NR:",		10),
	("Part:",		(2, Input.Particle.listAll)),
		("X:",		 8),
		("Y:",		 5),
		("N\u03A6:",	11),
	_None2,
		("Zmin:",	 9),
		("Zmax:",	 6),
		("NZ:",		12) ]

__USRBIN_layout4_Region = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("USRBIN"))),
		("Name:",	0),
	("Type:",		(1, _USRBIN_type)),
		("R1from:",	(7, LayoutManager.regionList)),
		("R1to:",	(4, LayoutManager.regionList)),
		("Step1:",	10),
	("Part:",		(2, Input.Particle.listAll)),
		("R2from:",	(8, LayoutManager.regionList)),
		("R2to:",	(5, LayoutManager.regionList)),
		("Step2:",	11),
	_None2,
		("R3from:",	(9, LayoutManager.regionList)),
		("R3to:",	(6, LayoutManager.regionList)),
		("Step3:",	12) ]

__USRBIN_layout4_Special = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("USRBIN"))),
		("Name:",	0),
	("Type:",		(1, _USRBIN_type)),
		("I1min:",	 7),
		("I1max:",	 4),
		("Step1:",	10),
	("Part:",		(2, Input.Particle.listAll)),
		("I2min:",	 8),
		("I2max:",	 5),
		("Step2:",	11),
	_None2,
		("F3min:",	 9),
		("F3max:",	 6),
		("N3:",		12) ]

def _USRBIN_layout4(card):
	# depending on the value of what(1) there are various layouts
	w1 = card.intWhat(1) % 10
	if w1 == 2:
		return __USRBIN_layout4_Region
	elif w1 == 8:
		return __USRBIN_layout4_Special
	elif w1 == 1 or w1 == 7:
		return __USRBIN_layout4_RZPhi
	else:
		return __USRBIN_layout4_XYZ

#------------------------------- EVENTBIN --------------------------------------
__EVENTBIN_type  = ("Type:", (lambda l,c,v=None,w=1:l._what_abs_val(c,w,v), _USRBIN_type))
__EVENTBIN_print_choice = [("All Cells", 0), ("Non-Zero Cells", 1) ]
__EVENTBIN_print	= ("Print:", (lambda l,c,v=None,w=1:l._what_abs_type(c,w,v),
			__EVENTBIN_print_choice))

__EVENTBIN_layout4_XYZ = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("EVENTBIN"))),
		("Name:",	0),
	__EVENTBIN_type,
		("Xmin:",	 7),
		("Xmax:",	 4),
		("NX:",		10),
	("Part:",		(2, Input.Particle.listAll)),
		("Ymin:",	 8),
		("Ymax:",	 5),
		("NY:",		11),
	__EVENTBIN_print,
		("Zmin:",	 9),
		("Zmax:",	 6),
		("NZ:",		12) ]

__EVENTBIN_layout4_RZPhi = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("EVENTBIN"))),
		("Name:",	0),
	__EVENTBIN_type,
		("Rmin:",	 7),
		("Rmax:",	 4),
		("NR:",		10),
	("Part:",		(2, Input.Particle.listAll)),
		("X:",		 8),
		("Y:",		 5),
		("N\u03A6:",	11),
	__EVENTBIN_print,
		("Zmin:",	 9),
		("Zmax:",	 6),
		("NZ:",		12) ]

__EVENTBIN_layout4_Region = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("EVENTBIN"))),
		("Name:",	0),
	__EVENTBIN_type,
		("R1from:",	(7, LayoutManager.regionList)),
		("R1to:",	(4, LayoutManager.regionList)),
		("Step1:",	10),
	("Part:",		(2, Input.Particle.listAll)),
		("R2from:",	(8, LayoutManager.regionList)),
		("R2to:",	(5, LayoutManager.regionList)),
		("Step2:",	11),
	__EVENTBIN_print,
		("R3from:",	(9, LayoutManager.regionList)),
		("R3to:",	(6, LayoutManager.regionList)),
		("Step3:",	12) ]

__EVENTBIN_layout4_Special = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("EVENTBIN"))),
		("Name:",	0),
	__EVENTBIN_type,
		("I1min:",	 7),
		("I1max:",	 4),
		("Step1:",	10),
	("Part:",		(2, Input.Particle.listAll)),
		("I2min:",	 8),
		("I2max:",	 5),
		("Step2:",	11),
	__EVENTBIN_print,
		("F3min:",	 9),
		("F3max:",	 6),
		("N3:",		12) ]

def _EVENTBIN_layout4(card):
	# depending on the value of what(1) there are various layouts
	w1 = abs(card.intWhat(1)) % 10
	if w1 == 2:
		return __EVENTBIN_layout4_Region
	elif w1 == 8:
		return __EVENTBIN_layout4_Special
	elif w1 == 1 or w1 == 7:
		return __EVENTBIN_layout4_RZPhi
	else:
		return __EVENTBIN_layout4_XYZ

#------------------------------- USRCOLL ---------------------------------------
_USRCOLL_layout4 = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("USRCOLL"))),
		("Name:",	0),
	("Type:",		(1, _log_what)),
		("Reg:",	(4, LayoutManager.allRegionList)),
		_None2,
		("Vol:",	 5),
	("Part:",		(2, Input.Particle.listAll)),
		("Emin:",	 8),
		("Emax:",	 7),
		("Bins:",	 6) ]

#------------------------------- USRTRACK --------------------------------------
_USRTRACK_type = [ _empty_what,
		("LinE,Pointwise   ", 10000 + 1),
		("LinE,Groupwise ", +1 ),
		("LogE,Pointwise   ", 10000 - 1),
		("LogE,Groupwise ", -1 )]
                 

_USRTRACK_layout4 = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("USRTRACK"))),
		("Name:",	0),
	("Type:",		(1, _USRTRACK_type)),
		("Reg:",	(4, LayoutManager.allRegionList)),
		_None2,
		("Vol:",	 5),
	("Part:",		(2, Input.Particle.listAll)),
		("Emin:",	 8),
		("Emax:",	 7),
		("Bins:",	 6) ]

#------------------------------- USRBDX ----------------------------------------
__linElinW =  1
__logElinW = -1
__linElogW =  2
__logElogW = -2

_USRBDX_type = [ _empty_what,
		("\u03A62,LogE,Lin\u03A9", 110 + __logElinW),
		("\u03A62,LogE,Log\u03A9", 110 + __logElogW),
		("\u03A62,LinE,Lin\u03A9", 110 + __linElinW),
		("\u03A62,LinE,Log\u03A9", 110 + __linElogW),

		("\u03A62,LogE,Lin\u03A9, point", 10110 + __logElinW),
		("\u03A62,LogE,Log\u03A9, point", 10110 + __logElogW),
		("\u03A62,LinE,Lin\u03A9, point", 10110 + __linElinW),
		("\u03A62,LinE,Log\u03A9, point", 10110 + __linElogW),

		("\u03A61,LogE,Lin\u03A9", 100 + __logElinW),
		("\u03A61,LogE,Log\u03A9", 100 + __logElogW),
		("\u03A61,LinE,Lin\u03A9", 100 + __linElinW),
		("\u03A61,LinE,Log\u03A9", 100 + __linElogW),

		("\u03A61,LogE,Lin\u03A9, point", 10100 + __logElinW),
		("\u03A61,LogE,Log\u03A9, point", 10100 + __logElogW),
		("\u03A61,LinE,Lin\u03A9, point", 10100 + __linElinW),
		("\u03A61,LinE,Log\u03A9, point", 10100 + __linElogW),

		("I2,LogE,Lin\u03A9", 10 + __logElinW),
		("I2,LogE,Log\u03A9", 10 + __logElogW),
		("I2,LinE,Lin\u03A9", 10 + __linElinW),
		("I2,LinE,Log\u03A9", 10 + __linElogW),

		("I2,LogE,Lin\u03A9, point", 10010 + __logElinW),
		("I2,LogE,Log\u03A9, point", 10010 + __logElogW),
		("I2,LinE,Lin\u03A9, point", 10010 + __linElinW),
		("I2,LinE,Log\u03A9, point", 10010 + __linElogW),

		("I1,LogE,Lin\u03A9", __logElinW),
		("I1,LogE,Log\u03A9", __logElogW),
		("I1,LinE,Lin\u03A9", __linElinW),
		("I1,LinE,Log\u03A9", __linElogW),

		("I1,LogE,Lin\u03A9, point", 10000 + __logElinW),
		("I1,LogE,Log\u03A9, point", 10000 + __logElogW),
		("I1,LinE,Lin\u03A9, point", 10000 + __linElinW),
		("I1,LinE,Log\u03A9, point", 10000 + __linElogW) ]

_USRBDX_layout4 = [ _None2,
		_None2,
		("Unit:",	(3, lambda l: l.unitList("USRBDX"))),
		("Name:",	0),
	("Type:",		(1, _USRBDX_type)),
		("Reg:",	(4, LayoutManager.regionList)),
		("to Reg:",	(5, LayoutManager.regionList)),
		("Area:",	 6),
	("Part:",		(2, Input.Particle.listAll)),
		("Emin:",	 8),
		("Emax:",	 7),
		("Ebins:",	 9),
	_None2,
		("\u03A9min:",	11),
		("\u03A9max:",	10),
		("\u03A9bins:",12) ]

#------------------------------- USRYIELD --------------------------------------
__USRYIELD_type = [ "Yield", "BEAMDEF" ]
__USRYIELD_iea_type  = [
		(""                      , 0) ,
		("Ekin GeV"              , 1) ,
		("Ptot GeV/c"            , 2) ,
		("Rapidity lab"          , 3) ,
		("rapidity CMS"          , 4) ,
		("Pseudorap. lab"        , 5) ,
		("Pseudorap. CMS"        , 6) ,
		("Feynman-x lab"         , 7) ,
		("Feynman-x CMS"         , 8) ,
		("Ptr GeV/c"             , 9) ,
		("Mtr GeV"               , 10),
		("Plong lab"             , 11),
		("Plong CMS"             , 12),
		("Etot GeV"              , 13),
		("Polar \u03b8 lab"     , 14),
		("polar \u03b8 CMS"     , 15),
		("Ptr^2 (GeV/c)**2"      , 16),
		("Weighted Angle lab"    , 17),
		("Weighted Ptr GeV/c"    , 18),
		("Plab/Pbeam"            , 19),
		("Etr"                   , 20),
		("Eexit"                 , 21),
		("Particle charge"       , 22),
		("Particle LET"          , 23),
		("Polar \u03b8 lab deg" , 24),
		("Polar \u03b8 CMS deg" , 25),
		("Ekin/nucl lab"         , 26),
		("P/nucl lab"            , 27),
		("Part.Bar.Charge"       , 28),
		("4 Mom.Tran. -t"        , 29),
		("X_F_long CMS"          , 30),
		("M_x squared"           , 31),
		("M_x squared / s"       , 32),
		("Time (s)"              , 33),
		("Sin avg lab"           , 34),
		("Ptot (GeV/c) cms"      , 35),
		("Etot (GeV/c) cms"      , 36),
		("Wei.Lab.Ang deg"	 , 37),
		("Sin avg lab deg"       , 38) ]

__USRYIELD_yield = [ _empty_what,
		#("Plain Double"                 ,  0),
		("d2\u03C3/dx1dx2"             ,  0),
		("d2\u03C3/dx1dx2"             ,  1),
		("E d3\u03C3/dp3"              ,  2),
		("d2N/dx1dx2"                   ,  3),
		("d2(x2*N)/dx1dx2"              ,  4),
		("d2(x1*N)/dx1dx2"              ,  5),
		("d2N/dx1dx2/cos(\u03b8)"      ,  6),
		("d2(x2^2*N)/dx1dx2"           ,  7),
		("d2(x1^2*N)/dx1dx2"           ,  8),
		("d2(x2*N)/dx1dx2/cos(\u03b8)" , 16),
		("d2(x1*N)/dx1dx2/cos(\u03b8)" , 26),
		("d2\u03C3/dx1dx2"             , 53),
		("d2(x2*\u03C3)/dx1dx2"        , 54),
		("d2(x1*\u03C3)/dx1dx2"        , 55),
		("d2\u03C3/dx1dx2/cos(\u03b8)" , 56),
		("d2(x2^2*\u03C3)/dx1dx2"      , 57),
		("d2(x1^2*\u03C3)/dx1dx2"      , 58),
		("d2(x2*\u03C3)/dx1dx2/cos(\u03b8)" , 66),
		("d2(x1*\u03C3)/dx1dx2/cos(\u03b8)" , 76) ]

__USRYIELD_yield2 = [ _empty_what,
		("Plain Double", 0),
		("d2\u03C3/dx1dx2", 1),
		("E d3\u03C3/dp3", 2),
		("d2N/dx1dx2", 3),
		("d2(x2*N)/dx1dx2", 4),
		("d2(x1*N)/dx1dx2", 5),
		("d2(x2^2*N)/dx1dx2", 7),
		("d2(x1^2*N)/dx1dx2", 8) ]

# Will only allow the change from BEAMDEF to something else
def __USRYIELD_sdum(layout, card, val=None):
	if val is None:
		sdum = card.sdum()[:8]
		if sdum=="BEAMDEF":
			return sdum
		else:
			return "Yield"
	else:
		if val=="BEAMDEF":
			# Remove extra whats
			undoinfo = [ layout.flair.setNWhatsUndo(card, 7) ]
			undoinfo.append( layout.flair.setWhatUndo(card, 0, val))
			layout.flair.addUndo(undo.createListUndo(undoinfo))
		elif card.sdum()=="BEAMDEF" and val=="Yield":
			layout.setWhat(card, 0, val)

__USRYIELD_layout4_beamdef = [
	_None2,
		("Type:",	(__USRYIELD_sdum, __USRYIELD_type)),
		("Proj:",	(1, Input.Particle.listAll)),
		("Target:",	(2, Input.Particle.listAll)),
	("p:",			3),
		("cosx:",	4),
		("cosy:",	5),
		("cosz:",	6) ]

# Set the WHAT(1) = ie + ia*100
def __iea(card):
	w1 = card.intWhat(1)
	if w1==0: w1=101
	(ip, w1) = divmod(w1, 10000)
	(ia, ie) = divmod(w1, 100)
	if ie>50:
		ia += 1
		ie = ie-100
	return (ip,ia,ie)

def __USRYIELD_ie(layout, card, val=None):
	(ip,ia,ie) = __iea(card)
	if val is None:
		return abs(ie)
	else:
		if ie<0: val = -val
		layout.setWhat(card, 1, ip*10000 + ia*100 + val)

def __USRYIELD_ia(layout, card, val=None):
	(ip, ia, ie) = __iea(card)
	if val is None:
		return ia
	else:
		layout.setWhat(card, 1, ip*10000 + val*100 + ie)

def __USRYIELD_ip(layout, card, val=None):
	(ip, ia, ie) = __iea(card)
	if val is None:
		return ip
	else:
		layout.setWhat(card, 1, val*10000 + ia*100 + ie)

def __USRYIELD_log(layout, card, val=None):
	(ip, ia,ie) = __iea(card)
	if val is None:
		if ie<0: return -1
		return 1
	else:
		layout.setWhat(card, 1, ip*10000 + ia*100 + val*abs(ie))

def __USRYIELD_kind(layout, card, val=None):
	(ixm, ixa) = divmod(card.intWhat(12), 100)
	if val is None:
		return ixa
	else:
		if val == "":
			layout.setWhat(card, 12, val)
		elif val == 0:
			layout.setWhat(card, 12, 0)
		else:
			layout.setWhat(card, 12, ixm*100 + val)

def __USRYIELD_kindMat(layout, card, val=None):
	(ixm, ixa) = divmod(card.intWhat(12), 100)
	if val is None:
		try:
			mat = layout.input.material(ixm, True)
		except:
			say("\nERROR:", sys.exc_info()[1])
			say(card)
			return ""
		return mat
	else:
		if ixa != 0:
			ixm = layout.input.material(val, False)
			layout.setWhat(card, 12, ixm*100 + ixa)

def __USRYIELD_particle(card):
	try:
		w = int(card.what(2))
		if w<-800:
			# w = -1000-ij <=> ij = -1000-w
			entering = True
			ij = -1000-w
			part = Input.Particle.get(ij)
		else:
			part = Input.Particle.get(w)
	except ValueError:
		entering = False
		part = Input.Particle.get(card.what(2))
	return part, entering

def __USRYIELD_what2(layout, card, val=None):
	part,entering = __USRYIELD_particle(card)

	if val is None:
		if part is not None:
			return part.name
		else:
			return ""
	else:
		newpart = Input.Particle.get(val)
		if newpart is None:
			layout.setWhat(card, 2, "")
		elif card.intWhat(4)==-1 and card.intWhat(5)==-2 and entering:
			layout.setWhat(card, 2, -1000-newpart.id)
		else:
			layout.setWhat(card, 2, val)

def __USRYIELD_emerging(layout, card, val=None):
	part,entering = __USRYIELD_particle(card)
	if val is None:
		if card.intWhat(4)==-1 and card.intWhat(5)==-2:
			if entering:
				return "ENTERING"
			else:
				return "EMERGING"
		else:
			return ""
	else:
		undoinfo = []
		if val in ("EMERGING", "ENTERING"):
			undoinfo.append(layout.flair.setWhatUndo(card, 4, -1))
			undoinfo.append(layout.flair.setWhatUndo(card, 5, -2))
			if part is not None:
				if val=="EMERGING":
					undoinfo.append(layout.flair.setWhatUndo(card, 2, part.name))
				else:
					undoinfo.append(layout.flair.setWhatUndo(card, 2, -1000-part.id))
		else:
			undoinfo.append(layout.flair.setWhatUndo(card, 4, ""))
			undoinfo.append(layout.flair.setWhatUndo(card, 5, ""))
		layout.flair.addUndo(undo.createListUndo(undoinfo))

__USRYIELD_layout4_score = [
	_None2,
		("Type:",	(__USRYIELD_sdum, __USRYIELD_type),
				"SDUM: Yield for user yield or BEAMDEF to define beam axes"),
		("Unit:",	(3, lambda l: l.unitList("USRYIELD"))),
		("Name:",	 0),
	("ie:",			(__USRYIELD_ie, __USRYIELD_iea_type),
				"WHAT(1) ie: quantity to score, continuous"),
		("ia:",		(__USRYIELD_ia, __USRYIELD_iea_type),
				"WHAT(1) ia: quantity to score, one interval"),
		("ip:",		(__USRYIELD_ip, _point_what),
				"WHAT(1) ip: neutron energy binning"),
		("Log:",	(__USRYIELD_log, _log_what),
				"WHAT(1) sign: logarithmic or linear scoring"),
		("Norm:",	 6),
	("Part:",		(__USRYIELD_what2, Input.Particle.listAll)),	# FIXME special <80 means ENTERING in interaction
		("Yield:",	(__USRYIELD_emerging, ["", "EMERGING", "ENTERING"]),
				"WHAT(4): Emerging from inelastic hadronic interactions is scored"),
		("Reg:",	(4, LayoutManager.regionList)),
		("to Reg:",	(5, LayoutManager.regionList)),
	("Min1:",		 8),
		("Max1:",	 7),
		("Nbins1:",	 9),
		_None2,
	("Min2:",		11),
		("Max2:",	10),
		("Kind:",	(__USRYIELD_kind, __USRYIELD_yield),
				"WHAT(12) ixa: kind of yield or cross section desired"),
		("Mat:",	(__USRYIELD_kindMat, LayoutManager.materialList),
				"WHAT(12) ixm: material of the target for cross section or LET calculations") ]

__USRYIELD_layout4_score_EMERGING = [
	_None2,
		("Type:",	(__USRYIELD_sdum, __USRYIELD_type)),
		("Unit:",	(3, lambda l: l.unitList("USRYIELD"))),
		("Name:",	 0),
	("ie:",			(__USRYIELD_ie, __USRYIELD_iea_type),
				"WHAT(1) ie: quantity to score, continuous"),
		("ia:",		(__USRYIELD_ia, __USRYIELD_iea_type),
				"WHAT(1) ia: quantity to score, one interval"),
		("ip:",		(__USRYIELD_ip, _point_what),
				"WHAT(1) ip: neutron energy binning"),
		("Log:",	(__USRYIELD_log, _log_what),
				"WHAT(1) sign: logarithmic or linear scoring"),
		("Norm:",	 6),
	("Part:",		(__USRYIELD_what2, Input.Particle.listAll)),
		("Yield:",	(__USRYIELD_emerging, ["", "EMERGING", "ENTERING"]),
				"WHAT(4): Emerging from inelastic hadronic interactions is scored"),
		_None2,
		_None2,
	("Min1:",		 8),
		("Max1:",	 7),
		("Nbins1:",	 9),
		_None2,
	("Min2:",		11),
		("Max2:",	10),
		("Kind:",	(__USRYIELD_kind, __USRYIELD_yield2),
				"WHAT(12) ixa: kind of yield or cross section desired"),
		("Mat:",	(__USRYIELD_kindMat, LayoutManager.materialList),
				"WHAT(12) ixm: material of the target for cross section or LET calculations") ]

def _USRYIELD_layout4(card):
	# depending on the value of sdum
	if card.sdum() == "BEAMDEF":
		return __USRYIELD_layout4_beamdef
	elif card.intWhat(4)==-1 and card.intWhat(5)==-2:
		return __USRYIELD_layout4_score_EMERGING
	else:
		return __USRYIELD_layout4_score

#------------------------------- USERDUMP --------------------------------------
__USERDUMP_type = [ "Dump", "UDQUENCH" ]

__USERDUMP_what1 = [ _empty_what, _ignore_what, ("Complete",100), _reset_what ]
__USERDUMP_what3 = [ _empty_what,
			("All",			0),
			("Source",		1),
			("Traj&Cont Losses",	2),
			("Local Losses",	3),
			("Source, Traj&Cont",	4),
			("Source, Local",	5),
			("Traj & All Losses",	6),
			("No Source, Traj, Losses",	7) ]
__USERDUMP_what4 = [	_empty_what,
			_ignore_what,
			("User Defined", 1),
			_reset_what]

def __USERDUMP_sdum(layout, card, val=None):
	if val is None:
		if card.sdum()=="UDQUENCH":
			return "UDQUENCH"
		else:
			return "Dump"
	else:
		layout.setWhat(card, 0, val)

__USERDUMP_layout4_dump = [ _None2,
		("Type:",	(__USERDUMP_sdum, __USERDUMP_type), "Dump: Creates a collission tape\nUDQUENCH to specify Birks parameter for quenching"),
		("Unit:",	(2, lambda l: l.unitList("USERDUMP",False))),
		("File:",	0),
	_None2,
		("What:",	(1, __USERDUMP_what1)),
		("Score:",	(3, __USERDUMP_what3)),
		("Dump:",	(4, __USERDUMP_what4)) ]

__USERDUMP_layout4_udquench = [ _None2,
		("Type:",	(__USERDUMP_sdum, __USERDUMP_type)),
		("1st Birks m1:",	1),
		("2 Birks m1:",	2),
	("1st Birks m2:",	3),
		("2nd Birks m2:",4),
		("1st Birks m3:",5),
		("2nd Birks m3:",6)]

def _USERDUMP_layout4(card):
	if card.sdum() == "UDQUENCH":
		return __USERDUMP_layout4_udquench
	else:
		return __USERDUMP_layout4_dump

#------------------------------- VOXELS ----------------------------------------
_VOXELS_layout4 = [ _None2,
		("x:",	1),
		("y:",	2),
		("z:",	3),
	_None2,
		("Trans:", (4, LayoutManager.rotdefiList)),
		("Filename:", (0, ("Voxel file","*.vxl"))),
		_None2]

#------------------------------- WW-FACTOr -------------------------------------
__WW_FACTO_profile = [ "", 1, 2, 3, 4, 5, -1 ]
_WW_FACTO_layout4 = [ _None2,
		("RR:",		1),
		("Split:",	2),
		("mult f:",	3),
	("LowE n:", (0, __WW_FACTO_profile)) ] + __reg456_range

#------------------------------- WW-PROFIle ------------------------------------
_WW_PROFI_layout4 = [ _None2,
		("WW extra f:",	1),
		("Importa f:",	2),
		("Profile:",	(5, list(range(1,6)))),
	_None2,
		("En low:",	(3, LayoutManager.neutronGroupList)),
		("En high:",	(4, LayoutManager.neutronGroupList)),
		("Step:",	5) ]

#------------------------------- WW-THRESh -------------------------------------
__WW_THRES_layout4_normal = [ _None2,
		("E upper:",	1),
		("E lower:",	2),
		("amp f:",	3),
	("Opt:", (0, __BIASING_prim)),
		("Part:",	(4, Input.Particle.list)),
		("to Part:",	(5, Input.Particle.list)),
		("Step:",	 6) ]

__WW_THRES_layout4_lowneut = [ _None2,
		("En upper:",	(1, LayoutManager.neutronGroupList)),
		("En lower:",	(2, LayoutManager.neutronGroupList)),
		("amp f:",	3),
	("Opt:", (0, __BIASING_prim)),
		("Part:",	(4, Input.Particle.list)),
		("to Part:",	(5, Input.Particle.list)),
		("Step:",	 6) ]

def _WW_THRES_layout4(card):
	if card.what(4)=="WWLOWNEU" or card.what(5)=="WWLOWNEU":
		return __WW_THRES_layout4_lowneut
	else:
		return __WW_THRES_layout4_normal

#------------------------------- POINT -----------------------------------------
__ANCHOR_opt  = [ ("none",0), ("C",1), ("N",2), ("NE",3), ("E",4), ("SE",5), ("S",6), ("SW",7), ("W",8), ("NW",9) ]
__POINT_opt   = [ ("cross",0), ("dot",1), ("square",2), ("x",3), ("diamond",4), ("circle",5) ]
__ARROW_opt   = [ ("head",0), ("tail",1), ("headtail",2), ("line",3) ]
__SPLINE_opt  = [ ("line", 0), ("cardinal",1) ]
__RULER_opt   = [ ("simple",0), ("angle",1) ]
__LIGHT_opt   = [ ("sun",0), ("omni",1), ("spot",2) ]
__FALLOFF_opt = [ ("constant",0), ("1/r",1), ("1/r2",2),
		  ("const+shadow", 3), ("1/r+shadow",4), ("1/r2+shadow",5) ]

def __ANCHOR_type(layout, card, val=None):
	s,a = divmod(card.intWhat(5),100)
	if val is None:
		return a
	else:
		try: aval = int(val)
		except: aval = 0
		layout.setWhat(card, 5, s*100 + aval)

def __SIZE_type(layout, card, val=None):
	s,a = divmod(card.intWhat(5),100)
	if val is None:
		return s
	else:
		try: sval = int(val)
		except: sval = 0
		layout.setWhat(card, 5, sval*100 + a)

_POINT_layout4 = [ (None, 0),
		("x:",	1),
		("y:",	2),
		("z:",	3),
	("Option",		(4, __POINT_opt)),
		("Size",	__SIZE_type),
		("Anchor:",	(__ANCHOR_type, __ANCHOR_opt)),
		("Color:",	6)]

#------------------------------- ARROW ----------------------------------------
_ARROW_layout4 = [ (None, 0),
		("x:",	1),
		("y:",	2),
		("z:",	3),
	_None2,
		("dx",	7),
		("dy",	8),
		("dz",	9),
	("Option",		(4, __ARROW_opt)),
		("Size",	__SIZE_type),
		("Anchor:",	(__ANCHOR_type, __ANCHOR_opt)),
		("Color:",	6)]

#------------------------------- RULER ----------------------------------------
_RULER_layout4 = [ (None, 0),
		("x:",	1),
		("y:",	2),
		("z:",	3),
	_None2,
		("xe",	7),
		("ye",	8),
		("ze",	9),
	_None2,
		("xa",	10),
		("ya",	11),
		("za",	12),
	("Option",		(4, __RULER_opt)),
		("Size",	__SIZE_type),
		("Anchor:",	(__ANCHOR_type, __ANCHOR_opt)),
		("Color:",	6)]

#------------------------------ SPLINE ----------------------------------------
def _SPLINE_spline4(card):
	layout = [ (None, 0),
			("x1:",	1),
			("y1:",	2),
			("z1:",	3),
		("Option",		(4, __SPLINE_opt)),
			("Size",	__SIZE_type),
			("Anchor:",	(__ANCHOR_type, __ANCHOR_opt)),
			("Color:",	6)]

	n = (card.nwhats()-7)//3
	if n&1==1: n += 1
	for i in range(n):
		j = i+2
		layout.append( _None2 )
		layout.append( ("x%d"%(j), i*3+7) )
		layout.append( ("y%d"%(j), i*3+8) )
		layout.append( ("z%d"%(j), i*3+9) )
	return layout

#------------------------------- LIGHT ----------------------------------------
_LIGHT_layout4 = [ (None, 0),
		("x:",	1),
		("y:",	2),
		("z:",	3),
	_None2,
		("dx",	7),
		("dy",	8),
		("dz",	9),
	("Option",		(4, __LIGHT_opt)),
		("Size",	__SIZE_type),
		("Anchor:",	(__ANCHOR_type, __ANCHOR_opt)),
		("Color:",	6),
	_None2,
		("Power",	10),
		("Falloff",	(11,__FALLOFF_opt)),
		("Specular",	12)]

#------------------------------- MESH ------------------------------------------
_MESH_layout4 = [ (None, 0),
		("x:",	1),
		("y:",	2),
		("z:",	3),
	("Option",		(4, __POINT_opt)),
		("Size",	__SIZE_type),
		("Anchor:",	(__ANCHOR_type, __ANCHOR_opt)),
		("Color:",	6),
	("File:",		(-1, ("STL files","*.stl"))),
		_None2,
		_None2,
		_None2 ]

#------------------------------- COFFEE ---------------------------------------
_COFFEE_layout4 = [ (None, None),
		("Type:",	(1,[	("Italian Long",0),
					("Italian Espresso",1),
					("Greek",2),
					("Frappe",3),
					("Swiss",4),
					("Turkish",5)])),
		("# Cups:",	2),
		("Milk:",	(3,[	("No milk",0),
					("Small Milk",1),
					("Large Milk",2),
					("Renverse", 3)])),
		("Sugar:",	4),
		("Cocoa:",	5),
		("Temperature:", 6),
		("Blend",	(0,["COLOMBIAN", "ARABICA", "MIXED"]))]

#===============================================================================
# get layout for card
#-------------------------------------------------------------------------------
def getLayout(card):
	ci = Input.CardInfo.get(card)
	if ci:
		layout = ci.layout
	else:
		layout = _generic_layout4(card)

	if callable(layout):
		return layout(card)
	elif layout is None:
		return _generic_layout4(card)
	else:
		return layout

#-------------------------------------------------------------------------------
# order cards...
#-------------------------------------------------------------------------------
def init():
	order = [("#&define"  , _P_DEFINE_layout4 ),
		( "#&undef"   , _P_UNDEF_layout4  ),
		( "#inc&lude" , _P_INCLUDE_layout4),
		( "#&if"      , _P_IF_layout4     ),
		( "#ifdef"    , _P_IFDEF_layout4  ),
		( "#ifndef"   , _P_IFNDEF_layout4 ),
		( "#eli&f"    , _P_ELIF_layout4   ),
		( "#&else"    , _P_ELSE_layout4   ),
		( "#e&ndif"   , _P_ENDIF_layout4  ),
		( "#endinclude", _empty_layout4   ),
		( "TITLE"     , _TITLE_layout4    ),
		( "FIXED"     , _empty_layout4    ),
		( "FREE"      , _empty_layout4    ),
		( "GLOBAL"    , _GLOBAL_layout4   ),
		( "DEFAULTS"  , _DEFAULTS_layout4 ),
		( "&BEAM"     , _BEAM_layout4     ),
		( "BEAM&POS"  , _BEAMPOS_layout4  ),
		( "BEAM&AXES" , _BEAMAXES_layout4 ),
		( "HI-PROPE"  , _HI_PROPE_layout4 ),
		( "SO&URCE"   , _SOURCE_layout4   ),
		( "SPECSOUR"  , _SPECSOUR_layout4 ),
		( "SPOTBEAM"  , _SPOTBEAM_layout4 ),
		( "SPOTPOS"   , _SPOTPOS_layout4  ),
		( "SPOTDIR"   , _SPOTDIR_layout4  ),
		( "SPOTTRAN"  , _SPOTTRAN_layout4 ),
		( "SPOTTIME"  , _SPOTTIME_layout4 ),
		( "GCR-SPE"   , _GCR_SPE_layout4  ),
		( "POLARIZA"  , _POLARIZA_layout4 ),
		( "DISCARD"   , _DISCARD_layout4  ),
		( "BME"       , _BME_layout4      ),
		( "DELTARAY"  , _DELTARAY_layout4 ),
		( "DPMJET"    , _DPMJET_layout4   ),
		( "EMFFLUO"   , _EMFFLUO_layout4  ),
		( "EMFRAY"    , _EMFRAY_layout4   ),
		( "EVENTYPE"  , _EVENTYPE_layout4 ),
		( "IONFLUCT"  , _IONFLUCT_layout4 ),
		( "IONTRANS"  , _IONTRANS_layout4 ),
		( "MULSOPT"   , _MULSOPT_layout4  ),
		( "MUPHOTON"  , _MUPHOTON_layout4 ),
		( "MYRQMD"    , _MYRQMD_layout4   ),
		( "OPT-PROD"  , _OPT_PROD_layout4 ),
		( "PAIRBREM"  , _PAIRBREM_layout4 ),
		( "PHOTONUC"  , _PHOTONUC_layout4 ),
		( "PHYSICS"   , _PHYSICS_layout4  ),
		( "RQMD"      , _RQMD_layout4     ),

		( "!point"    , _POINT_layout4    ),
		( "!arrow"    , _ARROW_layout4    ),
		( "!spline"   , _SPLINE_spline4   ),
		( "!ruler"    , _RULER_layout4    ),
		( "!mesh"     , _MESH_layout4     ),
		( "!coffee"   , _COFFEE_layout4   ),
		( "!light"    , _LIGHT_layout4    ),

		( "GEOBEGIN"  , _GEOBEGIN_layout4 ),
		( "VOXELS"    , _VOXELS_layout4   ),
		( "&RPP"      , _RPP_layout4      ),
		( "&TET"      , _TET_layout4      ), #ZXW20240830-----For TET, added by zxw
		( "&BOX"      , _BOX_layout4      ),
		( "&SPH"      , _SPH_layout4      ),
		( "R&CC"      , _RCC_layout4      ),
		( "R&EC"      , _REC_layout4      ),
		( "&TRC"      , _TRC_layout4      ),
		( "&TRX"      , _TRXYZ_layout4    ),
		( "&TRY"      , _TRXYZ_layout4    ),
		( "&TRZ"      , _TRXYZ_layout4    ),
		( "E&LL"      , _ELL_layout4      ),
		( "&WED"      , _WED_layout4      ),
		( "RAW"       , _RAW_layout4      ),
		( "&ARB"      , _ARB_layout4      ),
		( "&PLA"      , _PLA_layout4      ),
		( "&QUA"      , _QUA_layout4      ),
		( "&YZP"      , _YZP_layout4      ),
		( "X&ZP"      , _XZP_layout4      ),
		( "&XYP"      , _XYP_layout4      ),
		( "XCC"       , _XCC_layout4      ),
		( "YCC"       , _YCC_layout4      ),
		( "ZCC"       , _ZCC_layout4      ),
		( "XEC"       , _XEC_layout4      ),
		( "YEC"       , _YEC_layout4      ),
		( "ZEC"       , _ZEC_layout4      ),
		( "END"       , _empty_layout4	  ),
		( "$start_expansion", _START_EXPANSION_layout4 ),
		( "$start_translat",  _START_TRANSLAT_layout4  ),
		( "$start_transform", _START_TRANSFORM_layout4 ),
		( "$end_expansion",   _END_EXPANSION_layout4   ),
		( "$end_translat",    _END_TRANSLAT_layout4    ),
		( "$end_transform",   _END_TRANSFORM_layout4   ),
		( "REGION"    , _REGION_layout4   ),
		( "LATTICE"   , _LATTICE_layout4  ),
		( "LATTSNGL"  , _LATTSNGL_layout4 ),
		( "GEO&END"   , _GEOEND_layout4   ),
		( "&MATERIAL" , _MATERIAL_layout4 ),
		( "COMPOUND"  , _COMPOUND_layout4 ),
		( "LOW-MAT"   , _LOW_MAT_layout4  ),
		( "&OPT-PROP" , _OPT_PROP_layout4 ),
		( "MAT-&PROP" , _MAT_PROP_layout4 ),
		( "CORRFACT"  , _CORRFACT_layout4 ),
		( "STERNHEI"  , _STERNHEI_layout4 ),
		( "ASSIGNMA"  , _ASSIGNMA_layout4 ),
		( "EMF"       , _EMF_layout4      ),
		( "EMF&CUT"   , _EMFCUT_layout4   ),
		( "EMFFIX"    , _EMFFIX_layout4   ),
		( "ELCFIELD"  , _ELCFIELD_layout4 ),
		( "EXPTRANS"  , _EXPTRANS_layout4 ),
		( "&FLUKAFIX" , _FLUKAFIX_layout4 ),
		( "&LOW-NEUT" , _LOW_NEUT_layout4 ),
		( "MCSTHRES"  , _MCSTHRES_layout4 ),
		( "MGNFIELD"  , _MGNFIELD_layout4 ),
		( "&PART-THR" , _PART_THR_layout4 ),
		( "STEPSIZE"  , _STEPSIZE_layout4 ),
		( "THRESHOL"  , _THRESHOL_layout4 ),
		( "TIME-&CUT" , _TIME_CUT_layout4 ),
		( "BIASING"   , _BIASING_layout4  ),
		( "EMF-BIAS"  , _EMF_BIAS_layout4 ),
		( "LAM-BIAS"  , _LAM_BIAS_layout4 ),
		( "LOW-BIAS"  , _LOW_BIAS_layout4 ),
		( "LOW-DOWN"  , _LOW_DOWN_layout4 ),
		( "LOW-PWXS"  , _LOW_PWXS_layout4 ),
		( "WW-FACTO"  , _WW_FACTO_layout4 ),
		( "WW-PROFI"  , _WW_PROFI_layout4 ),
		( "WW-THRES"  , _WW_THRES_layout4 ),
		( "BAMJET"    , _BAMJET_layout4   ),
		( "DPM-PARA"  , _DPM_PARA_layout4 ),
		( "EVXTEST"   , _EVXTEST_layout4  ),
		( "DCYSCORE"  , _DCYSCORE_layout4 ),
		( "DCYTIMES"  , _DCYTIMES_layout4 ),
		( "IRRPROFI"  , _IRRPROFI_layout4 ),
		( "RADDECAY"  , _RADDECAY_layout4 ),
		( "RAD-BIOL"  , _RAD_BIOL_layout4 ),
		( "ROT-&DEFI" , _ROT_DEFI_layout4 ),
		( "ROTPRBIN"  , _ROTPRBIN_layout4 ),
		( "AUXSCORE"  , _AUXSCORE_layout4 ),
		( "TCQUENCH"  , _TCQUENCH_layout4 ),
		( "TPSSCORE"  , _TPSSCORE_layout4 ),
		( "TSL-PWXS"  , _TSL_PWXS_layout4 ),
		( "DETECT"    , _DETECT_layout4   ),
		( "USR&BIN"   , _USRBIN_layout4   ),
		( "USRBD&X"   , _USRBDX_layout4   ),
		( "USR&COLL"  , _USRCOLL_layout4  ),
		( "USR&TRACK" , _USRTRACK_layout4 ),
		( "USR&YIELD" , _USRYIELD_layout4 ),
		( "&EVENTBIN" , _EVENTBIN_layout4 ),
		( "EVENTDAT"  , _EVENTDAT_layout4 ),
		( "&USERDUMP" , _USERDUMP_layout4 ),
		( "USER&WEIG" , _USERWEIG_layout4 ),
		( "&RESNUCLE" , _RESNUCLE_layout4 ),
		( "&SCORE"    , _SCORE_layout4    ),
		( "OPEN"      , _OPEN_layout4     ),
		( "USRGCALL"  , _USRGCALL_layout4 ),
		( "USRICALL"  , _USRICALL_layout4 ),
		( "USROCALL"  , _USROCALL_layout4 ),
		( "PLOTGEOM"  , _PLOTGEOM_layout4 ),
		( "&RANDOMIZ" , _RANDOMIZ_layout4 ),
		( "&START"    , _START_layout4    ),
		( "STOP"      , _STOP_layout4     )]

	# Set ordering and underline
	i = 0
	for (t,l) in order:
		underline = t.find('&')
		if underline>=0:
			tag = t.replace('&','')
		else:
			tag = t
		ci = Input.CardInfo.get(tag)
		ci.order     = i
		ci.layout    = l
		ci.underline = underline
		i += 1
