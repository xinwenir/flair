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
# Date:	13-Jun-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import io
import os
import re
import sys
import math
import time
import random
import string
import binascii
from log import say
from operator import attrgetter
try:
	import pickle as pickle
except ImportError:
	import pickle

try:
	from io import StringIO
except ImportError:
	from io import StringIO

try:
	from tkinter import *
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	import tkinter.messagebox as messagebox

import Unicode

import csg
import Ribbon
import tkFlair
import tkExtra
import FlairRibbon
import GeometryBody
import GeometryExtra
import GeometryViewer
import GeometryLayers

import bmath
import undo
import Input
import Manual
import Project
import Palette

# For debugging, deactivate otherwise
import pdb
import traceback

from Constants import *

# Body colors for Visible, Locked similar to the ones in GeometryViewer
# BUT for Tkinter notation (Geometry Viewer = int) (Tkinter = string #xxxxxx)
VISIBLE_COLOR   = "DarkMagenta"
LOCK_COLOR      = "Orange"
FREEZE_COLOR    = "#00A0A0"
ERROR_COLOR     = "Red"
WIREFRAME_COLOR = "#007070"

ADDZONE         = "+zone"
UNDOTIMER       = 200

#FONT = tkFont.Font(font=("Helvetica",8))
#FONT_LIST   = "Helvetica 9"

try:
	OPER_TRANS      = string.maketrans("+-|()\n\t","       ")
except AttributeError:
	OPER_TRANS      =    str.maketrans("+-|()\n\t","       ")

# Warning on the trailing space! Is need to avoid accidentally the user typing in
# a ON/OFF/MULTI pattern
_ON    = "[X] "
_OFF   = "[  ] "
_MULTI = "[-] "

_PROPTABS = ["Properties","Attributes"]

# Tuple of card tags that contain gobject information to be drawn
_CARD_OBJECTS = ("BEAM", "ROT-DEFI", "SPOTBEAM")

#-------------------------- Bodies Cards ---------------------------------------
_RPP_prop = [	("Xmin",  1),
		("Xmax",  2),
		("Ymin",  3),
		("Ymax",  4),
		("Zmin",  5),
		("Zmax",  6),
		("@Xmid",-1),
		("@Ymid",-2),
		("@Zmid",-3),
		("@Dx",  -4),
		("@Dy",  -5),
		("@Dz",  -6)]

_BOX_prop = [	("x",     1),
		("y",     2),
		("z",     3),
		("H1x",   4),
		("H1y",   5),
		("H1z",   6),
		("H2x",   7),
		("H2y",   8),
		("H2z",   9),
		("H3x",  10),
		("H3y",  11),
		("H3z",  12),
		("@Xmid",-1),
		("@Ymid",-2),
		("@Zmid",-3),
		("@Dx",  -4),
		("@Dy",  -5),
		("@Dz",  -6)]

_SPH_prop = [	("x",     1),
		("y",     2),
		("z",     3),
		("R",     4),
		("@Diam", -1)]

_RCC_prop = [	("x",      1),
		("y",      2),
		("z",      3),
		("Hx",     4),
		("Hy",     5),
		("Hz",     6),
		("R",      7),
		("@Xmid", -1),
		("@Ymid", -2),
		("@Zmid", -3),
		("@H",    -4),
		("@Diam", -5),
		("@azm",  -6),
		("@polar",-7)]

_REC_prop = [	("x",      1),
		("y",      2),
		("z",      3),
		("Hx",     4),
		("Hy",     5),
		("Hz",     6),
		("R1x",    7),
		("R1y",    8),
		("R1z",    9),
		("R2x",   10),
		("R2y",   11),
		("R2z",   12),
		("@Xmid", -1),
		("@Ymid", -2),
		("@Zmid", -3),
		("@H",    -4),
		("@Rx",   -5),
		("@Ry",   -6),
		("@azm",  -7),
		("@polar",-8)]

_TRC_prop = [	("x",      1),
		("y",      2),
		("z",      3),
		("Hx",     4),
		("Hy",     5),
		("Hz",     6),
		("Rbase",  7),
		("Rapex",  8),
		("@Xmid", -1),
		("@Ymid", -2),
		("@Zmid", -3),
		("@H",    -4),
		("@azm",  -5),
		("@polar",-6)]

_ELL_prop = [	("F1x",   1),
		("F1y",   2),
		("F1z",   3),
		("F2x",   4),
		("F2y",   5),
		("F2z",   6),
		("L",     7),

		("@Cx",    -1),
		("@Cy",    -2),
		("@Cz",    -3),
		("@Rmajor",-4),
		("@Rminor",-5)]

_WED_prop = [	("x",     1),
		("y",     2),
		("z",     3),
		("H1x",   4),
		("H1y",   5),
		("H1z",   6),
		("H2x",   7),
		("H2y",   8),
		("H2z",   9),
		("H3x",  10),
		("H3y",  11),
		("H3z",  12),
		("@Dx",  -1),
		("@Dy",  -2),
		("@Dz",  -3)]

_RAW_prop = _WED_prop

_ARB_prop = [	("V1x",	1),
		("V1y",	2),
		("V1z",	3),
		("V2x",	4),
		("V2y",	5),
		("V2z",	6),
		("V3x",	7),
		("V3y",	8),
		("V3z",	9),
		("V4x",	10),
		("V4y",	11),
		("V4z",	12),
		("V5x",	13),
		("V5y",	14),
		("V5z",	15),
		("V6x",	16),
		("V6y",	17),
		("V6z",	18),
		("V7x",	19),
		("V7y",	20),
		("V7z",	21),
		("V8x",	22),
		("V8y",	23),
		("V8z",	24),
		("F1",	25),
		("F2",	26),
		("F3",	27),
		("F4",	28),
		("F5",	29),
		("F6",	30) ]

_XYP_prop = [	("z",	1) ]

_YZP_prop = [	("x",	1) ]

_XZP_prop = [	("y",	1) ]

_PLA_prop = [	("Nx",	   1),
		("Ny",	   2),
		("Nz",	   3),
		("x",	   4),
		("y",	   5),
		("z",	   6),
		("@azm",  -1),
		("@polar",-2)]

_XCC_prop = [	("y",	1),
		("z",	2),
		("R",	3),
		("@Diam", -1)]

_YCC_prop = [	("z",	1),
		("x",	2),
		("R",	3),
		("@Diam", -1)]

_ZCC_prop = [	("x",	1),
		("y",	2),
		("R",	3),
		("@Diam", -1)]

_XEC_prop = [	("y",	1),
		("z",	2),
		("Ly",	3),
		("Lz",	4) ]

_YEC_prop = [	("z",	1),
		("x",	2),
		("Lz",	3),
		("Lx",	4) ]

_ZEC_prop = [	("x",	1),
		("y",	2),
		("Lx",	3),
		("Ly",	4) ]

_QUA_prop = [	("Cxx",	1),
		("Cyy",	2),
		("Czz",	3),
		("Cxy",	4),
		("Cxz",	5),
		("Cyz",	6),
		("Cx",	7),
		("Cy",	8),
		("Cz",	9),
		("C",	10) ]

_TRXYZ_prop = [	("x",	1),
		("y",	2),
		("z",	3),
		("a",	4),
		("b",	5),
		("c",	6) ]

_VOXEL_prop = [	("x",	1),
		("y",	2),
		("z",	3),
		("rotdefi", 4)]

_BODIES = {	"ARB": _ARB_prop,
		"BOX": _BOX_prop,
		"ELL": _ELL_prop,
		"PLA": _PLA_prop,
		"QUA": _QUA_prop,
		"RAW": _RAW_prop,
		"RCC": _RCC_prop,
		"REC": _REC_prop,
		"RPP": _RPP_prop,
		"SPH": _SPH_prop,
		"TRC": _TRC_prop,
		"WED": _WED_prop,
		"XCC": _XCC_prop,
		"XEC": _XEC_prop,
		"XYP": _XYP_prop,
		"XZP": _XZP_prop,
		"YCC": _YCC_prop,
		"YEC": _YEC_prop,
		"YZP": _YZP_prop,
		"ZCC": _ZCC_prop,
		"ZEC": _ZEC_prop,
		"TRX": _TRXYZ_prop,
		"TRY": _TRXYZ_prop,
		"TRZ": _TRXYZ_prop,
		"VOXELS" : _VOXEL_prop }

_POINT_prop = [ ("option",	4),
		("anchor",	5),
		("size",	5),
		("color",	6),
		("x",		1),
		("y",		2),
		("z",		3)]

_ARROW_prop = _POINT_prop[:]
_ARROW_prop.extend([	("dx",	7),
			("dy",	8),
			("dz",	9),
			("@length", -1)])

_SPLINE_prop = _POINT_prop[:]
for i in range(0, 19):
	_SPLINE_prop.append(("x%d"%(i+2), i*3+7))
	_SPLINE_prop.append(("y%d"%(i+2), i*3+8))
	_SPLINE_prop.append(("z%d"%(i+2), i*3+9))

_RULER_prop = _POINT_prop[:]
_RULER_prop.extend([	("xe",	7),
			("ye",	8),
			("ze",	9),
			("xa",	10),
			("ya",	11),
			("za",	12),
			("@length",-1),
			("@angle", -2)])

_LIGHT_prop = _POINT_prop[:]
_LIGHT_prop.extend([	("dx",	7),
			("dy",	8),
			("dz",	9),
			("power",	10),
			("falloff",	11),
			("spec",	12)])

_MESH_prop = _POINT_prop[:]
_MESH_prop.append(("file", -1))

_MESH_prop.append(("file", -1))

_POINT_opt  = [ "cross", "dot", "square", "x", "diamond", "circle" ]
_ARROW_opt  = [ "head", "tail", "headtail", "line" ]
_SPLINE_opt = [ "line", "cardinal" ]
_RULER_opt  = [ "simple", "angle" ]
_LIGHT_opt  = [ "sun", "omni", "spot" ]
_MESH_opt   = [ "stl" ]

_OBJECTS = {	"!point"  : (_POINT_prop,  _POINT_opt),
		"!arrow"  : (_ARROW_prop,  _ARROW_opt),
		"!spline" : (_SPLINE_prop, _SPLINE_opt),
		"!mesh"   : (_MESH_prop,   _MESH_opt),
		"!ruler"  : (_RULER_prop,  _RULER_opt),
		"!light"  : (_LIGHT_prop,  _LIGHT_opt)
	}

_ANCHOR_list = [ "none", "C", "N", "NE", "E", "SE", "S", "SW", "W", "NW" ]

_ONOFF_list  = [ "off", "on"]
_FALLOFF_list = [ "constant", "1/r", "1/r2", "const+shadow", "1/r+shadow", "1/r2+shadow" ]

# Patterns
_FIRST_WORD = re.compile(r"^\s*(\S+)\s+(.*)$",re.DOTALL)
_LAST_NUM   = re.compile(r"^(.*?)(\d+)$")

_ASSIGNMA_field = [ "",
		"Magnetic",
		"Electric",
		"Magnetic+Electric",
		"Prompt Magnetic",
		"Prompt Electric",
		"Prompt Magnetic+Electric",
		"Decay Magnetic",
		"Decay Electric",
		"Decay Magnetic+Electric" ]

#===============================================================================
class InPlaceText(tkExtra.InPlaceText):
	def defaultBinds(self):
		tkExtra.InPlaceText.defaultBinds(self)
		self.edit.bind("<Escape>", self.ok)

#===============================================================================
class InPlaceRegion(InPlaceText):
	# ----------------------------------------------------------------------
	def __init__(self, listbox, flair, item=ACTIVE, value=None, **kw):
		self.flair = flair
		self._inBodies = False
		InPlaceText.__init__(self, listbox, item, value, **kw)

	# ----------------------------------------------------------------------
	def defaultBinds(self):
		InPlaceText.defaultBinds(self)
		self.edit.bind("<Key>",  self._insertBody)

	# ----------------------------------------------------------------------
	def focusOut(self, event=None):
		if self._inBodies:
			return "break"
		else:
			self.ok()

	# ----------------------------------------------------------------------
	def _insertBody(self, event):
		if event.char in ("+","-") or event.keysym == "Insert":
			self.edit.insert(INSERT, event.char)
			try:
				x,y,x2,y2 = self.edit.bbox(INSERT)
			except TypeError:
				x = y = 0
			x += self.edit.winfo_rootx()
			y += self.edit.winfo_rooty()
			self.flair._layout.updateInput()
			bodies = self.flair._layout.bodyList()
			self._inBodies = True
			self.toplevel.grab_release()
			selected = tkExtra.PopupList(self.toplevel, bodies).show(x,y)
			self._inBodies = False
			if selected:
				try: selected = selected.split('[')[0]
				except: pass
				self.edit.insert(INSERT,selected)

			# After closing a toplevel with overrideredirect(1)
			# the whole application is loosing the focus so a
			# simple focus_set() is not sufficient. I have to use
			# focus_force() on the toplevel
			self.toplevel.focus_force()
			self.edit.focus_set()
			self.toplevel.grab_set()
			return "break"

#===============================================================================
class Properties:
	_lastActive  = 0	# last element active in listbox(s)
	_lastActive2 = 0	# remember them in case for a moment everything is
				# unselected, when an object is selected back
				# to restore the active element

	# ----------------------------------------------------------------------
	def __init__(self, editor, cards, freeze=False, visible=False):
		self.editor   = editor
		self.flair    = editor.flair
		self.listbox  = editor.propList
		self.listbox2 = editor.attrList
		self.cards    = cards
		self.many     = len(cards)>1
		self.input    = self.editor.project.input
		self.items    = ["name", "comment"]
		self.items2   = ["name"]

		for card in cards:
			if card.input is not None and card.input.format==Input.FORMAT_VOXEL:
				self.voxel = True
				break
		else:
			self.voxel = False

		if freeze:	self.items2.append("freeze")
		if visible:	self.items2.append("visible")

	# ----------------------------------------------------------------------
	def type(self): return self.cards[0].type()
	def name(self): return self.cards[0].name()

	# ----------------------------------------------------------------------
	# Return common tag of all cards if any otherwise None
	# ----------------------------------------------------------------------
	def tag(self):
		tag = self.cards[0].tag
		for c in self.cards:
			if tag != c.tag:
				return None
		return tag

	# ----------------------------------------------------------------------
	def populateList(self):
		if self.listbox.size()  > 0: Properties._lastActive  = self.listbox.index(ACTIVE)
		if self.listbox2.size() > 0: Properties._lastActive2 = self.listbox2.index(ACTIVE)
		self.listbox.delete(0,END)
		self.listbox2.delete(0,END)
		for item in self.items:
			self.listbox.insert(END, (item, self.value(item)))
		for item in self.items2:
			self.listbox2.insert(END, (item, self.value(item)))
		self.listbox.activate(Properties._lastActive)
		self.listbox2.activate(Properties._lastActive2)

	# ----------------------------------------------------------------------
	def value(self, var=None, val=None):
		if var == "name":
			if self.many:
				return "-"
			elif val is None:
				return self.cards[0].name()
			else:
				self.editor.changeName(self.cards[0].name(), val)
				return None

		elif var == "comment":
			if val is None:
				v = self.cards[0].comment()
				for c in self.cards:
					if v != c.comment():
						return "-"
				return v
			else:
				undoinfo = []
				for c in self.cards:
					undoinfo.append(self.flair.setCommentUndo(c,val))
				self.addUndo(undoinfo)

		elif var=="visible":
			if val is None:
				v = self.cards[0].get(SELECT,0)&BIT_VISIBLE
				for c in self.cards:
					if v != c.get(SELECT,0)&BIT_VISIBLE:
						return _MULTI
				if v: return _ON
				else: return _OFF
			else:
				for c in self.cards:
					if val: c[SELECT] = c.get(SELECT,0) | BIT_VISIBLE
					else:   c[SELECT] = c.get(SELECT,0) & BIT_UNVISIBLE

		elif var=="freeze":
			if val is None:
				v = self.cards[0].get(SELECT,0)&BIT_FREEZE
				for c in self.cards:
					if v != c.get(SELECT,0)&BIT_FREEZE:
						return _MULTI
				if v: return _ON
				else: return _OFF
			else:
				for c in self.cards:
					if val: c[SELECT] = c.get(SELECT,0) | BIT_FREEZE
					else:   c[SELECT] = c.get(SELECT,0) & BIT_UNFREEZE

		return ""

	# ----------------------------------------------------------------------
	def edit(self, listbox):
		var, value = listbox.get(ACTIVE)
		if var == "name" and self.many:
			return None
		elif var == "comment":
			return InPlaceText(listbox.lists[1])
		else:
			return tkExtra.InPlaceEdit(listbox.lists[1], exportselection=False)

	# ----------------------------------------------------------------------
	# return true if fields accepts only one number
	# then multiple values to be split over the next fields
	# e.g. a copy paste of a vector ##x## ##y## ##z##
	# will be split in to the next fields
	# ----------------------------------------------------------------------
	def split(self):
		return False

	def split2(self):
		return False

	# ----------------------------------------------------------------------
	def addUndo(self, undoinfo):
		self.flair.addUndo(undoinfo)

#===============================================================================
class BodyProperties(Properties):
	def __init__(self, editor, cards):
		Properties.__init__(self, editor, cards, freeze=True, visible=True)
		self.items.append("type")
		t = self.tag()
		if t:
			self.bodytype = _BODIES[t]
			self.items.extend([x for x,y in self.bodytype])
		else:
			self.bodytype = None
		self.items2.extend(["bbox","wireframe",
				"dx","dy","dz","expansion","transform","transform2"])

	# ----------------------------------------------------------------------
	def value(self, var=None, val=None):
		if var=="type":
			if val is None:
				if self.bodytype: return self.cards[0].tag
				else: return "-"

		elif var=="wireframe":
			if val is None:
				v = self.cards[0].get(SELECT,0)&BIT_WIREFRAME
				for c in self.cards:
					if v != c.get(SELECT,0)&BIT_WIREFRAME:
						return _MULTI
				if v: return _ON
				else: return _OFF
			else:
				for c in self.cards:
					if val: c[SELECT] = c.get(SELECT,0) | BIT_WIREFRAME
					else:   c[SELECT] = c.get(SELECT,0) & BIT_UNWIREFRAME

		elif var=="bbox":
			if val is None:
				v = self.cards[0].get(SELECT,0)&BIT_BBOX
				for c in self.cards:
					if v != c.get(SELECT,0)&BIT_BBOX:
						return _MULTI
				if v: return _ON
				else: return _OFF
			else:
				for c in self.cards:
					if val: c[SELECT] = c.get(SELECT,0) | BIT_BBOX
					else:   c[SELECT] = c.get(SELECT,0) & BIT_UNBBOX

		elif var=="dx":
			if val is None:
				v = self.cards[0]["@dx"]
				for c in self.cards:
					if v != c["@dx"]:
						return ""
				if v is None:
					return ""
				else:
					return v
			else:
				for c in self.cards: c["@dx"] = val
				self.flair.addUndo(self.flair.correctBodyPropertiesUndo())

		elif var=="dy":
			if val is None:
				v = self.cards[0]["@dy"]
				for c in self.cards:
					if v != c["@dy"]:
						return ""
				if v is None:
					return ""
				else:
					return v
			else:
				for c in self.cards: c["@dy"] = val
				self.flair.addUndo(self.flair.correctBodyPropertiesUndo())

		elif var=="dz":
			if val is None:
				v = self.cards[0]["@dz"]
				for c in self.cards:
					if v != c["@dz"]:
						return ""
				if v is None:
					return ""
				else:
					return v
			else:
				for c in self.cards: c["@dz"] = val
				self.flair.addUndo(self.flair.correctBodyPropertiesUndo())

		elif var=="expansion":
			if val is None:
				v = self.cards[0]["@expansion"]
				for c in self.cards:
					if v != c["@expansion"]:
						return ""
				if v is None:
					return ""
				else:
					return v
			else:
				for c in self.cards: c["@expansion"] = val
				self.flair.addUndo(self.flair.correctBodyPropertiesUndo())

		elif var=="transform":
			if val is None:
				v = self.cards[0]["@transform"]
				for c in self.cards:
					if v != c["@transform"]:
						return ""
				if v is None:
					return ""
				else:
					return v
			else:
				for c in self.cards: c["@transform"] = val
				self.flair.addUndo(self.flair.correctBodyPropertiesUndo())
		elif var=="transform2":
			if val is None:
				v = self.cards[0]["@transform2"]
				for c in self.cards:
					if v != c["@transform2"]:
						return ""
				if v is None:
					return ""
				else:
					return v
			else:
				for c in self.cards: c["@transform2"] = val
				self.flair.addUndo(self.flair.correctBodyPropertiesUndo())


		elif var in ("@Dx", "@H"):
			if val is None:
				v = None
				for c in self.cards:
					if c.tag=="RPP":
						d = str(c.numWhat(2) - c.numWhat(1))
					else:
						d = str(math.sqrt(c.numWhat(4)**2 + c.numWhat(5)**2 + c.numWhat(6)**2))
					if v is None: v = d
					elif v != d: return ""
				return bmath.format(v,10)

			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					if c.tag == "RPP":
						undoinfo.append(
							self.flair.setWhatUndo(c, 2, val+c.numWhat(1)))
					else:
						if c.tag in ("BOX", "WED", "RAW"):
							V = c.bodyX()
						else:	# RCC, TRC, REC
							V = c.bodyZ()
						V.normalize()
						V = V*val
						undoinfo.append(self.flair.setWhatUndo(c, 4, V[0]))
						undoinfo.append(self.flair.setWhatUndo(c, 5, V[1]))
						undoinfo.append(self.flair.setWhatUndo(c, 6, V[2]))
				self.addUndo(undoinfo)

		elif var in ("@Dy", "@Rx"):
			if val is None:
				v = None
				for c in self.cards:
					if c.tag=="RPP":
						d = str(c.numWhat(4) - c.numWhat(3))
					else:
						d = str(math.sqrt(c.numWhat(7)**2 + c.numWhat(8)**2 + c.numWhat(9)**2))
					if v is None: v = d
					elif v != d: return ""
				return bmath.format(v,10)

			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					if c.tag == "RPP":
						undoinfo.append(
							self.flair.setWhatUndo(
								c, 4, val+c.numWhat(3)))
					else:
						if c.tag in ("BOX", "WED", "RAW"):
							V = c.bodyY()
						else:
							V = c.bodyX()
						V.normalize()
						V = V*val
						undoinfo.append(self.flair.setWhatUndo(c, 7, V[0]))
						undoinfo.append(self.flair.setWhatUndo(c, 8, V[1]))
						undoinfo.append(self.flair.setWhatUndo(c, 9, V[2]))
				self.addUndo(undoinfo)

		elif var in ("@Dz", "@Ry"):
			if val is None:
				v = None
				for c in self.cards:
					if c.tag=="RPP":
						d = str(c.numWhat(6) - c.numWhat(5))
					else:
						d = str(math.sqrt(c.numWhat(10)**2 + c.numWhat(11)**2 + c.numWhat(12)**2))
					if v is None: v = d
					elif v != d: return ""
				return bmath.format(v,10)

			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					if c.tag == "RPP":
						undoinfo.append(
							self.flair.setWhatUndo(
								c, 6, val+c.numWhat(5)))
					else:
						if c.tag in ("BOX", "WED", "RAW"):
							V = c.bodyZ()
						else:
							V = c.bodyY()
						V.normalize()
						V = V*val
						undoinfo.append(self.flair.setWhatUndo(c, 10, V[0]))
						undoinfo.append(self.flair.setWhatUndo(c, 11, V[1]))
						undoinfo.append(self.flair.setWhatUndo(c, 12, V[2]))
				self.addUndo(undoinfo)

		elif var == "@Diam":
			if val is None:
				v = None
				for c in self.cards:
					if c.tag == "RCC":
						d = str(2*c.numWhat(7))
					elif c.tag == "SPH":
						d = str(2*c.numWhat(4))
					elif c.tag in ("XCC", "YCC", "ZCC"):
						d = str(2*c.numWhat(3))
					if v is None: v = d
					elif v != d: return ""
				return bmath.format(v,10)

			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					if c.tag == "RCC":
						undoinfo.append(self.flair.setWhatUndo(c, 7, val/2))
					elif c.tag == "SPH":
						undoinfo.append(self.flair.setWhatUndo(c, 4, val/2))
					elif c.tag in ("XCC", "YCC", "ZCC"):
						undoinfo.append(self.flair.setWhatUndo(c, 3, val/2))
				self.addUndo(undoinfo)

		elif var == "@Xmid":
			if val is None:
				v = None
				for c in self.cards:
					if c.tag == "RPP":
						d = str((c.numWhat(1) + c.numWhat(2))/2.0)
					elif c.tag in ("RCC","TRC","REC"):
						h = c.bodyZ()
						d = str(c.numWhat(1) + h[0]/2.0)
					elif c.tag == "BOX":
						h = c.bodyX() + c.bodyY() + c.bodyZ()
						d = str(c.numWhat(1) + h[0]/2.0)
					if v is None: v = d
					elif v != d: return ""
				return bmath.format(v,10)

			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					if c.tag == "RPP":
						d = (c.numWhat(2) - c.numWhat(1))/2.0
						undoinfo.append(self.flair.setWhatUndo(c, 1, val-d))
						undoinfo.append(self.flair.setWhatUndo(c, 2, val+d))
					elif c.tag in ("RCC","TRC","REC"):
						h = c.bodyZ()
						undoinfo.append(self.flair.setWhatUndo(c, 1, val-h[0]/2.0))
					elif c.tag == "BOX":
						h = c.bodyX() + c.bodyY() + c.bodyZ()
						undoinfo.append(self.flair.setWhatUndo(c, 1, val-h[0]/2.0))
				self.addUndo(undoinfo)

		elif var == "@Ymid":
			if val is None:
				v = None
				for c in self.cards:
					if c.tag == "RPP":
						d = str((c.numWhat(3) + c.numWhat(4))/2.0)
					elif c.tag in ("RCC","TRC","REC"):
						h = c.bodyZ()
						d = str(c.numWhat(2) + h[1]/2.0)
					elif c.tag == "BOX":
						h = c.bodyX() + c.bodyY() + c.bodyZ()
						d = str(c.numWhat(2) + h[1]/2.0)
					if v is None: v = d
					elif v != d: return ""
				return bmath.format(v,10)

			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					if c.tag == "RPP":
						d = (c.numWhat(4) - c.numWhat(3))/2.0
						undoinfo.append(self.flair.setWhatUndo(c, 3, val-d))
						undoinfo.append(self.flair.setWhatUndo(c, 4, val+d))
					elif c.tag in ("RCC","TRC","REC"):
						h = c.bodyZ()
						undoinfo.append(self.flair.setWhatUndo(c, 2, val-h[1]/2.0))
					elif c.tag == "BOX":
						h = c.bodyX() + c.bodyY() + c.bodyZ()
						undoinfo.append(self.flair.setWhatUndo(c, 2, val-h[1]/2.0))
				self.addUndo(undoinfo)

		elif var == "@Zmid":
			if val is None:
				v = None
				for c in self.cards:
					if c.tag == "RPP":
						d = str((c.numWhat(5) + c.numWhat(6))/2.0)
					elif c.tag in ("RCC","TRC","REC"):
						h = c.bodyZ()
						d = str(c.numWhat(3) + h[2]/2.0)
					elif c.tag == "BOX":
						h = c.bodyX() + c.bodyY() + c.bodyZ()
						d = str(c.numWhat(3) + h[2]/2.0)
					if v is None: v = d
					elif v != d: return ""
				return bmath.format(v,10)

			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					if c.tag == "RPP":
						d = (c.numWhat(6) - c.numWhat(5))/2.0
						undoinfo.append(self.flair.setWhatUndo(c, 5, val-d))
						undoinfo.append(self.flair.setWhatUndo(c, 6, val+d))
					elif c.tag in ("RCC","TRC","REC"):
						h = c.bodyZ()
						undoinfo.append(self.flair.setWhatUndo(c, 3, val-h[2]/2.0))
					elif c.tag == "BOX":
						h = c.bodyX() + c.bodyY() + c.bodyZ()
						undoinfo.append(self.flair.setWhatUndo(c, 3, val-h[2]/2.0))
				self.addUndo(undoinfo)

		elif var == "@azm":
			if val is None:
				v = None
				for c in self.cards:
					if c.tag == "PLA":
						d = str(math.degrees(c.bodyN().phi()))
					elif c.tag in ("RCC","TRC","REC"):
						d = str(math.degrees(c.bodyZ().phi()))
					if v is None: v = d
					elif v != d: return ""
				return bmath.format(v,10)

			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					if c.tag == "PLA":
						n      = c.bodyN()
						ma    = n.length()
						phi   = n.phi()
						theta = n.theta()
						n.setPolar(ma, math.radians(val), theta)
						undoinfo.append(self.flair.setWhatUndo(c, 1, n[0]))
						undoinfo.append(self.flair.setWhatUndo(c, 2, n[1]))
						undoinfo.append(self.flair.setWhatUndo(c, 3, n[2]))

					elif c.tag in ("RCC","TRC","REC"):
						n     = c.bodyZ()
						ma    = n.length()
						phi   = n.phi()
						theta = n.theta()
						n.setPolar(ma, math.radians(val), theta)
						undoinfo.append(self.flair.setWhatUndo(c, 4, n[0]))
						undoinfo.append(self.flair.setWhatUndo(c, 5, n[1]))
						undoinfo.append(self.flair.setWhatUndo(c, 6, n[2]))
				self.addUndo(undoinfo)

		elif var == "@polar":
			if val is None:
				v = None
				for c in self.cards:
					if c.tag == "PLA":
						d = str(math.degrees(c.bodyN().theta()))
					elif c.tag in ("RCC","TRC","REC"):
						d = str(math.degrees(c.bodyZ().theta()))
					if v is None: v = d
					elif v != d: return ""
				return bmath.format(v,10)

			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					if c.tag == "PLA":
						n = c.bodyN()
						ma    = n.length()
						phi   = n.phi()
						theta = n.theta()
						n.setPolar(ma, phi, math.radians(val))
						undoinfo.append(self.flair.setWhatUndo(c, 1, n[0]))
						undoinfo.append(self.flair.setWhatUndo(c, 2, n[1]))
						undoinfo.append(self.flair.setWhatUndo(c, 3, n[2]))
					elif c.tag in ("RCC","TRC","REC"):
						n     = c.bodyZ()
						ma    = n.length()
						phi   = n.phi()
						theta = n.theta()
						n.setPolar(ma, phi, math.radians(val))
						undoinfo.append(self.flair.setWhatUndo(c, 4, n[0]))
						undoinfo.append(self.flair.setWhatUndo(c, 5, n[1]))
						undoinfo.append(self.flair.setWhatUndo(c, 6, n[2]))
				self.addUndo(undoinfo)

		# Special for ELL!
		elif var in ("@Cx", "@Cy", "@Cz"):
			idx = ord(var[-1])-ord('x')
			if val is None:
				v = None
				for c in self.cards:
					F1 = c.bodyP()
					F2 = c.bodyP2()
					C  = 0.5*(F1 + F2)
					d  = C[idx]
					if abs(d)<1e-10: d=0.0
					d  = bmath.format(d,10)
					if v is None: v = d
					elif v != d: return ""
				return v
			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					F1 = c.bodyP()
					F2 = c.bodyP2()
					C  = 0.5*(F1 + F2)
					d  = val - C[idx]
					undoinfo.append(self.flair.setWhatUndo(c, 1+idx, F1[idx]+d))
					undoinfo.append(self.flair.setWhatUndo(c, 4+idx, F2[idx]+d))
				self.addUndo(undoinfo)

		# Special for ELL!
		elif var == "@Rmajor":
			if val is None:
				v = None
				for c in self.cards:
					d = str(c.numWhat(7)/2.0)
					if v is None: v = d
					elif v != d: return ""
				return bmath.format(v,10)
			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					# Current
					F1 = c.bodyP()
					F2 = c.bodyP2()
					C  = 0.5*(F1 + F2)
					Z  = F2 - F1
					f  = Z.normalize()/2.0
					RM = c.numWhat(7)/2.0
					rm2 = (RM-f)*(RM+f)
					# rm2 = RM2 - f2

					# New values
					# RM2 = rm2 - f2
					try: f  = math.sqrt(val**2 - rm2)
					except: f = 0.0
					F1 = C + f * Z
					F2 = C - f * Z

					for i in range(3):
						undoinfo.append(self.flair.setWhatUndo(c, 1+i, F1[i]))
						undoinfo.append(self.flair.setWhatUndo(c, 4+i, F2[i]))
					undoinfo.append(self.flair.setWhatUndo(c, 7, 2.0*val))

				self.addUndo(undoinfo)

		# Special for ELL!
		elif var == "@Rminor":
			if val is None:
				v = None
				for c in self.cards:
					F1 = c.bodyP()
					F2 = c.bodyP2()
					Z  = F2 - F1
					f  = Z.normalize()/2.0
					RM = c.numWhat(7)/2.0
					rm = (RM-f)*(RM+f)
					# rm2 = RM2 - f2
					if rm<0.0:
						rm = str(0.0)
					else:
						rm = str(math.sqrt(rm))
					if v is None: v = rm
					elif v != rm: return ""
				return v
			else:
				if val=="": return
				try: val = float(val)
				except ValueError: return
				undoinfo = []
				for c in self.cards:
					# Current
					F1 = c.bodyP()
					F2 = c.bodyP2()
					C  = 0.5*(F1 + F2)
					Z  = F2 - F1
					f  = Z.normalize()/2.0
					RM = c.numWhat(7)/2.0
					rm2 = (RM-f)*(RM+f)

					# New values
					# RM2 = rm2 - f2
					try: f  = math.sqrt(RM**2 - val**2)
					except: f = 0.0
					F1 = C + f * Z
					F2 = C - f * Z
					for i in range(3):
						undoinfo.append(self.flair.setWhatUndo(c, 1+i, F1[i]))
						undoinfo.append(self.flair.setWhatUndo(c, 4+i, F2[i]))
					undoinfo.append(self.flair.setWhatUndo(c, 7, 2.0*RM))

				self.addUndo(undoinfo)

		else:
			if not self.bodytype:
				return Properties.value(self, var, val)

			for lbl,what in self.bodytype:
				if var==lbl: break
			else:
				return Properties.value(self, var, val)

			if val is None:
				v = self.cards[0].what(what)
				for c in self.cards:
					if v != c.what(what): return ""
				return v
			else:
				if val=="" and self.cards[0].tag != "VOXELS": return
				undoinfo = []
				for c in self.cards:
					undoinfo.append(self.flair.setWhatUndo(c,what,val))
				self.addUndo(undoinfo)
				return

	# ----------------------------------------------------------------------
	def edit(self, listbox):
		var, value = listbox.get(ACTIVE)
		if var == "type":
			return None
		elif var == "name":
			if self.many: return None
			return tkExtra.InPlaceMaxLength(listbox.lists[1], maxlength=8, exportselection=True)
		elif var in ("transform", "rotdefi"):
			return tkExtra.InPlaceList(listbox.lists[1],
				values=self.input.rotdefiList(True))
		else:
			return Properties.edit(self,listbox)

	# ----------------------------------------------------------------------
	def split(self):
		return int(self.listbox.index(ACTIVE)) >= 3

	# ----------------------------------------------------------------------
	def split2(self):
		return int(self.listbox2.index(ACTIVE)) >= 4

#===============================================================================
# FIXME
# Allow linking bodies. Bodies refer to other bodies and the only thing that
# changes is the transformation of the body
#===============================================================================
class LinkBodyProperties(Properties):
	def __init__(self, editor, cards):
		Properties.__init(self, editor, cards)

#===============================================================================
class RegionProperties(Properties):
	def __init__(self, editor, cards):
		Properties.__init__(self, editor, cards, freeze=True)

		if self.voxel:
			self.items.extend(["material"])
		else:
			self.items.extend(["material", ADDZONE])
		self.items2.extend(["bbox", "alpha", "naz", "rotdefi", "field",
				"imp-All", "imp-H", "imp-E", "imp-N"])
		# FIXME NAZ
		# move material, rotdefi to attributes
		self.nzones = 0

	# ----------------------------------------------------------------------
	def populateList(self):
		# Dynamically adjust number of zones
		if not self.many and not self.voxel:
			zones = self.cards[0]["@zone"]
			if zones:
				if self.nzones != len(zones):
					del self.items[-self.nzones-1:]
					self.nzones = len(zones)
					for i in range(self.nzones):
						self.items.append("zone%02d"%(i+1))
					self.items.append(ADDZONE)
			else:
				del self.items[3:]
				self.items.append(ADDZONE)
				self.nzones = 0
		Properties.populateList(self)
		self.setColor()

	# ----------------------------------------------------------------------
	def setColor(self):
		color = self.editor.viewers.materialColor(self.cards[0]["@material"])
		self.listbox.lists[1].itemconfig(2,background="#%06X"%(color))

	# ----------------------------------------------------------------------
	def value(self, var=None, val=None):
		if var=="bbox":
			if val is None:
				v = self.cards[0].get(SELECT,0)&BIT_BBOX
				for c in self.cards:
					if v != c.get(SELECT,0)&BIT_BBOX:
						return _MULTI
				if v: return _ON
				else: return _OFF
			else:
				for c in self.cards:
					if val: c[SELECT] = c.get(SELECT,0) | BIT_BBOX
					else:   c[SELECT] = c.get(SELECT,0) & BIT_UNBBOX

		elif var == "material":
			if val is None:
				try:
					v = self.cards[0]["@material"].name()
					for c in self.cards:
						if v != c["@material"].name():
							return "-"
					return v
				except:
					return "error"
			else:
				# Set last material
				self.editor._lastMaterial = val

				# find the new material card to set the color
				for card in self.input["MATERIAL"]:
					if card.ignore(): continue
					if card.name()==val:
						for c in self.cards:
							c["@material"] = card
						self.setColor()
						break
				else:
					for card in Input._defaultMaterials[:] + Input._icruMaterials:
						if card.name()==val:
							for c in self.cards:
								c["@material"] = card
							self.setColor()
							break

				# find which ASSIGNMAt card set the material
				# FIXME collect all undo in a list
				undoinfo = []
				for c in self.cards:
					assignmat = c["@assignmat"]
					if assignmat is None:
						# Create a new assignmat card in input
						pos = self.input.bestPosition("ASSIGNMA")
						card = Input.Card("ASSIGNMA", ["",val,c.name()])
						undoinfo.append(self.flair.addCardUndo(card, pos, True))
						c["@assignmat"] = card

					elif assignmat.what(3)=="" or assignmat.what(2)==assignmat.what(3):
						# Set region assignmat
						undoinfo.append(self.flair.setWhatUndo(assignmat, 1, val))

					else:
						# Leave range untouched and save it below
						card = Input.Card("ASSIGNMA", ["",val,c.name()])
						undoinfo.append(self.flair.addCardUndo(card, assignmat.pos()+1, True))
						c["@assignmat"] = card

				if undoinfo:
					self.addUndo(undoinfo)
					self.flair.refresh("card")
					self.input.regionProperties()
					self.editor.updateGeometry()
					self.editor.redraw()	# FIXME maybe unnecessary call

		elif var == "field":
			if val is None:
				try:
					v = self.cards[0]["@assignmat"].intWhat(5)
					for c in self.cards:
						if v != c["@assignmat"].intWhat(5):
							return "-"
					return _ASSIGNMA_field[v]
				except:
					return "error"
			else:
				try:
					ival = _ASSIGNMA_field.index(val)
				except IndexError:
					return
				# find which ASSIGNMAt card set the material
				# FIXME collect all undo in a list
				undoinfo = []
				for c in self.cards:
					assignmat = c["@assignmat"]
					if assignmat is None:
						# Create a new assignmat card in input
						pos = self.input.bestPosition("ASSIGNMA")
						card = Input.Card("ASSIGNMA", ["","",c.name(),"","",ival])
						undoinfo.append(self.flair.addCardUndo(card, pos, True))
						c["@assignmat"] = card

					elif assignmat.what(3)=="" or assignmat.what(2)==assignmat.what(3):
						# Set region assignmat
						undoinfo.append(self.flair.setWhatUndo(assignmat, 5, ival))

					else:
						# Leave range untouched and save it below
						card = Input.Card("ASSIGNMA", ["",c.what(1),c.name(),c.what(3),
								c.what(4),ival,c.what(6)])
						undoinfo.append(self.flair.addCardUndo(card, assignmat.pos()+1, True))
						c["@assignmat"] = card

				if undoinfo:
					self.addUndo(undoinfo)
					self.flair.refresh("card")
					self.input.regionProperties()
					self.editor.updateGeometry()
					self.editor.redraw()	# FIXME maybe unnecessary call

		elif var == "rotdefi":
			if val is None:
				v = self.cards[0].get("@rotdefi","")
				for c in self.cards:
					if v != c.get("@rotdefi",""):
						return "-"
				return v
			else:
				undoinfo = []
				for c in self.cards:
					lattice = c["@lattice"]

					if val=="":
						if lattice is None: continue
						if lattice.what(2)=="" or lattice.what(1)==lattice.what(2):
							undoinfo.append(
								self.flair.delCardUndo(lattice.pos(), True))
							c["@rotdefi"] = ""
							c["@lattice"] = None
						else:
							# cannot  rotdefi use flair
							say("Cannot remove LATTICE card for REGION ",c.name())
					else:
						try:	# for integer
							rotdefi = "Rot%03d"%(int(val))
						except:
							rotdefi = val

						c["@rotdefi"] = rotdefi
						if lattice is None:
							# Create a new lattice card in input
							pos = self.input.bestPosition("LATTICE")
							card = Input.Card("LATTICE", [rotdefi,c.name()])
							undoinfo.append(self.flair.addCardUndo(card, pos, True))
							c["@lattice"] = card

						elif lattice.what(2)=="" or lattice.what(1)==lattice.what(2):
							# Set region lattice
							undoinfo.append(
								self.flair.setWhatUndo(lattice, 0, val))

						else:
							# Leave range untouched and save it below
							card = Input.Card("LATTICE",[rotdefi,c.name()])
							undoinfo.append(
								self.flair.addCardUndo(card, lattice.pos()+1, True))
							c["@lattice"] = card

				if undoinfo:
					self.addUndo(undoinfo)
					self.flair.refresh("card")
					self.input.regionProperties()
					self.editor.updateGeometry()
					self.editor.redraw()	# FIXME maybe unnecessary call

		elif var in ("imp-All", "imp-H", "imp-E", "imp-N"):
			if var == "imp-All":
				tag   = "@biasingA"
				what1 = 0
			elif var == "imp-H":
				tag = "@biasingH"
				what1 = 1
			elif var == "imp-E":
				tag = "@biasingE"
				what1 = 2
			elif var == "imp-N":
				tag = "@biasingN"
				what1 = 3

			if val is None:
				biasing = self.cards[0][tag]
				if biasing is not None:
					v = biasing.what(3)
				else:
					v = ""
				for c in self.cards:
					biasing = c[tag]
					if biasing is not None:
						vi = biasing.what(3)
					else:
						vi = ""
					if v != vi: return "-"
				return v
			else:
				undoinfo = []

				for c in self.cards:
					biasing = c[tag]

					if val=="":
						if biasing is None: continue
						if biasing.what(2)=="":
							undoinfo.append(
								self.flair.delCardUndo(biasing.pos(), True))
							c[tag] = None
						else:
							# cannot delete biasing use flair
							say("Cannot remove BIASING card for REGION", c.name())
					else:
						if biasing is None:
							# Create a new biasing card in input
							pos = self.input.bestPosition("BIASING")
							card = Input.Card("BIASING", ["",what1,"",val,c.name()])
							undoinfo.append(self.flair.addCardUndo(card, pos, True))
							c[tag] = card

						elif biasing.what(5)=="" or biasing.what(4)==biasing.what(5):
							# Set region importance
							undoinfo.append( self.flair.setWhatUndo(biasing, 3, val))

						else:
							# Leave range untouched and save it below
							card = Input.Card("BIASING", ["",what1,"",val,c.name()])
							undoinfo.append(
								self.flair.addCardUndo(card, biasing.pos()+1, True))
							c[tag] = card

				if undoinfo:
					self.addUndo(undoinfo)
					self.flair.refresh("card")
					self.editor.loadProject()	# FIXME maybe unnecessary call

		elif var.startswith("zone"):
			if self.many:
				return "-"
			i = int(var[4:])-1
			zones = self.cards[0]["@zone"]
			if i >= len(zones): return ""
			if val is None:
				return zones[i]
			else:
				if val != zones[i]:
					# create new expression by substituting the zone
					zones[i] = val
					self.flair.addUndo(self.flair.setZonesUndo(self.cards[0],zones))

		elif var == ADDZONE:
			if self.many:
				return "-"
			elif val is None:
				return ""
			else:
				if val: self.flair.addUndo(
					   self.flair.addZoneUndo(self.cards[0], val))

		elif var == "alpha":
			if val is None:
				v = self.cards[0].get(ALPHA,0)
				for c in self.cards:
					if v != c.get(ALPHA,0):
						return ""
				return v
			else:
				try: val = int(val)
				except: val = 0
				if val<0: val=0
				elif val>100: val=100

				undoinfo = []
				for c in self.cards:
					undoinfo.append(self.flair.setCardUndo(c, ALPHA, val, 0))

				if undoinfo:
					self.addUndo(undoinfo)
					self.editor.redraw()
					self.editor.setInputModified()

		else:
			return Properties.value(self, var, val)

	# ----------------------------------------------------------------------
	def edit(self, listbox):
		var, value = listbox.get(ACTIVE)
		if var == "material":
			return tkExtra.InPlaceList(listbox.lists[1],
				values=self.input.materialList(icru=True))
		elif var == "rotdefi":
			return tkExtra.InPlaceList(listbox.lists[1],
				values=self.input.rotdefiList(True))
		elif var.startswith("zone") or var==ADDZONE:
			if self.many: return None
			return InPlaceRegion(listbox.lists[1], self.flair)
		elif var == "name":
			if self.many: return None
			return tkExtra.InPlaceMaxLength(listbox.lists[1], maxlength=8, exportselection=True)
		elif var == "field":
			return tkExtra.InPlaceList(listbox.lists[1], values=_ASSIGNMA_field)
		elif var == "alpha":
			return tkExtra.InPlaceInteger(listbox.lists[1], exportselection=True)
		else:
			return Properties.edit(self, listbox)

#===============================================================================
class BeamProperties(Properties):
	def __init__(self, editor, cards):
		Properties.__init__(self, editor, cards, visible=True)
		self.items.extend(["energy","scale"])

	# ----------------------------------------------------------------------
	def populateList(self):
		Properties.populateList(self)

	# ----------------------------------------------------------------------
	def value(self, var=None, val=None):
		if var == "name":
			return Properties.value(self,var,val)

		elif var == "energy":
			if val is None:
				e = None
				for c in self.cards:
					if c.tag=="BEAM":
						energy = c.numWhat(1)
					else:
						energy = c.numWhat(2)
					if e is None:
						e = energy
					elif e != energy: return ""
				return e
			else:
				undoinfo = []
				for c in self.cards:
					if c.tag=="BEAM":
						w = 1
					else:
						w = 2
					undoinfo.append(self.flair.setWhatUndo(c, w, val))

				if undoinfo:
					self.addUndo(undoinfo)
					self.flair.refresh("card")
					self.editor.reload()	# Very nasty!!!

		elif var == "scale":
			if val is None:
				v = self.cards[0].get("scale",1.0)
				for c in self.cards:
					if v != c.get("scale",1.0): return ""
				return v
			else:
				undoinfo = []
				for c in self.cards:
					undoinfo.append(self.flair.setCardUndo(c, "scale", val, 1.0))

				if undoinfo:
					self.addUndo(undoinfo)
					self.editor.reload()	# Very nasty!!!

		else:
			return Properties.value(self, var, val)

	# ----------------------------------------------------------------------
	def edit(self, listbox):
		var, value = listbox.get(ACTIVE)
		if var == "color":
			edit = tkExtra.InPlaceColor(listbox.lists[1])
			return edit
		else:
			return Properties.edit(self, listbox)

#===============================================================================
class RotdefiProperties(Properties):
	def __init__(self, editor, cards):
		Properties.__init__(self, editor, cards, visible=True)
		self._rotdefi = self.input.rotdefiCards(self.cards[0].name())

	# ----------------------------------------------------------------------
	def populateList(self):
		# Dynamically adjust number of cards
		if not self.many:
			n = (len(self.items)-2) // 6
			if n != len(self._rotdefi):
				del self.items[2:]
				for i in range(1,len(self._rotdefi)+1):
					self.items.append("dx.%d"%(i))
					self.items.append("dy.%d"%(i))
					self.items.append("dz.%d"%(i))
					self.items.append("axis.%d"%(i))
					self.items.append("azm.%d"%(i))
					self.items.append("polar.%d"%(i))
				self.items.append("<add>")
		Properties.populateList(self)

	# ----------------------------------------------------------------------
	def value(self, var=None, val=None):
		# Search all patterns with prefix and what#
		for prefix,what in (("dx.",4), ("dy.",5), ("dz.",6), ("azm.",3), ("polar.",2)):
			if not var.startswith(prefix): continue
			if self.many: return "-"
			try:
				idx = int(var[len(prefix):])-1
				card = self._rotdefi[idx]
			except:
				if val is None: return ""
				# Add a new card!

			if val is None:
				return card.what(what)
			else:
				self.flair.addUndo(
					self.flair.setWhatUndo(
						card, what, val))
			return

		if var.startswith("axis."):
			if self.many: return "-"
			idx = int(var[5:])-1
			card = self._rotdefi[idx]

			w1 = card.intWhat(1)
			if w1 < 1000:
				(j, i) = divmod(w1, 100)
			else:
				(i, j) = divmod(w1, 1000)

			if val is None:
				try: return _ROTDEFI_axes[j]
				except: return "?"

			else:
				try: j = _ROTDEFI_axes.index(val)
				except: return
				if w1 < 1000:
					w1 = i + j*100
				else:
					w1 = j + i*1000
				self.flair.addUndo(self.flair.setWhatUndo(card, 1, w1))

		elif var == "<add>":
			if val is None:
				return ""
			else:
				# add a new ROTDEFI card
				pos = self._rotdefi[-1].pos()+1
				card = Input.Card("ROT-DEFI",[self._rotdefi[0].name()])
				self.flair.addUndo(self.flair.addCardUndo(card, pos, True))
				self.flair.refresh("card")
				self._rotdefi.append(card)
				self.populateList()

		else:
			return Properties.value(self, var, val)

	# ----------------------------------------------------------------------
	def edit(self, listbox):
		var, value = listbox.get(ACTIVE)
		if var.startswith("axis."):
			return tkExtra.InPlaceList(listbox.lists[1],
				values=_ROTDEFI_axes)
		elif var == "<add>":
			return tkExtra.InPlaceList(listbox.lists[1],
				values=["","add"])
		return Properties.edit(self, listbox)

	# ----------------------------------------------------------------------
	def split(self):
		return int(self.listbox.index(ACTIVE)) >= 2

#===============================================================================
class ObjectProperties(Properties):
	def __init__(self, editor, cards):
		Properties.__init__(self, editor, cards, visible=True)
		self.items.append("type")
		self.objtype = cards[0].tag
		for c in cards:
			if self.objtype != c.tag:
				self.objtype = None
				break
		if self.objtype:
			self.objprop, self.option = _OBJECTS[self.objtype]
			self.items.extend([x for x,y in self.objprop])
		else:
			self.objprop = self.option = None

	# ----------------------------------------------------------------------
	def populateList(self):
		Properties.populateList(self)
		self.setColor()

	# ----------------------------------------------------------------------
	def setColor(self):
		if self.objprop is not None:
			self.listbox.lists[1].itemconfig(6,
				background="#%06X"%(self.cards[0].intWhat(6)))

	# ----------------------------------------------------------------------
	def value(self, var=None, val=None):
		if var == "name":
			return Properties.value(self,var,val)

		elif var == "type":
			if val is None:
				if self.objprop: return self.cards[0].tag
				else: return "-"

		elif var == "option":
			if val is None:
				v = self.cards[0].intWhat(4)
				for c in self.cards:
					if v != c.intWhat(4): return ""
				try: return self.option[v]
				except: return self.option[0]
			else:
				try:
					opt = self.option.index(val)
					# FIXME undoinfo!!!
					undoinfo = []
					for c in self.cards:
						undoinfo.append(self.flair.setWhatUndo(c, 4, opt))
					self.addUndo(undoinfo)
				except IndexError:
					pass

		elif var == "size":
			if val is None:
				s,a = divmod(self.cards[0].intWhat(5),100)
				for c in self.cards:
					if s != c.intWhat(5) // 100: return ""
				return s
			else:
				s,a = divmod(self.cards[0].intWhat(5),100)
				try:
					ival = int(val)
					for c in self.cards:
						self.flair.addUndo(
							self.flair.setWhatUndo(
								c, 5, ival*100+a))
				except:
					pass

		elif var == "anchor":
			if val is None:
				s,a = divmod(self.cards[0].intWhat(5),100)
				for c in self.cards:
					if a != c.intWhat(5) % 100: return ""
				try:
					return _ANCHOR_list[a]
				except:
					return _ANCHOR_list[0]
			else:
				s,a = divmod(self.cards[0].intWhat(5),100)
				try:
					ival = _ANCHOR_list.index(val)
					for c in self.cards:
						self.flair.addUndo(
							self.flair.setWhatUndo(
								c, 5, s*100+ival))
				except IndexError:
					pass

		elif var == "color":
			if val is None:
				v = self.cards[0].intWhat(6)
				for c in self.cards:
					if v != c.intWhat(6): return ""
				return "#%06X"%(v)
			else:
				try:
					color = int(str(val)[1:], 16)
					for c in self.cards:
						self.flair.addUndo(
							self.flair.setWhatUndo(c, 6, val))
				except:
					pass

		elif var == "falloff":
			if val is None:
				v = self.cards[0].intWhat(11)
				for c in self.cards:
					if v != c.intWhat(11): return ""
				try:
					return _FALLOFF_list[v]
				except IndexError:
					return _FALLOFF_list[0]
			else:
				try:
					for c in self.cards:
						self.flair.addUndo(
							self.flair.setWhatUndo(
								c, 11, _FALLOFF_list.index(val)))
				except:
					pass

		elif var == "@length":
			if self.many: return None
			card = self.cards[0]
			if card.tag not in ("!arrow", "!ruler"): return None
			P = bmath.Vector(card.numWhat(1), card.numWhat(2), card.numWhat(3))
			D = bmath.Vector(card.numWhat(7), card.numWhat(8), card.numWhat(9))
			if card.tag == "!ruler": D -= P
			if val is None:
				return bmath.format(D.length(),10)
			else:
				D.normalize()
				try:
					undoinfo = [self.flair.refreshUndo(),
						self.flair.saveCardUndoInfo(card.pos())]
					D = float(val)*D
					if card.tag == "!ruler": D += P
					card.setWhat(7, bmath.format(D[0],10))
					card.setWhat(8, bmath.format(D[1],10))
					card.setWhat(9, bmath.format(D[2],10))
					undoinfo.append(self.flair.refreshUndo())
					self.flair.addUndo(undo.createListUndo(undoinfo,"object length"))
				except:
					pass

		elif var == "@angle":
			if self.many: return None
			card = self.cards[0]
			P = bmath.Vector(card.numWhat( 1), card.numWhat( 2), card.numWhat( 3))
			E = bmath.Vector(card.numWhat( 7), card.numWhat( 8), card.numWhat( 9))
			A = bmath.Vector(card.numWhat(10), card.numWhat(11), card.numWhat(12))
			D  = E-P
			D2 = A-P
			l = D.normalize()
			D2.normalize()
			if val is None:
				return bmath.format(math.degrees(math.acos(D*D2)),10)
			else:
				try: ang = math.radians(float(val))
				except: return
				W = D2 ^ D		# Find perpendicular vector
				l = W.normalize()
				if l>0.001:
					M = bmath.Matrix(4)	# Create rotation matrix
					M.rotate(ang,W)
					E = l*(M*D2) + P	# Find new ending position
					undoinfo = [self.flair.refreshUndo(),
						self.flair.saveCardUndoInfo(card.pos())]
					card.setWhat(7, bmath.format(E[0],10))
					card.setWhat(8, bmath.format(E[1],10))
					card.setWhat(9, bmath.format(E[2],10))
					undoinfo.append(self.flair.refreshUndo())
					self.flair.addUndo(undo.createListUndo(undoinfo,"object length"))

		else:
			if self.objprop is None:
				return Properties.value(self, var, val)

			for lbl,what in self.objprop:
				if var==lbl: break
			else:
				return Properties.value(self, var, val)

			if val is None:
				v = self.cards[0].what(what)
				for c in self.cards:
					if v != c.what(what): return ""
				return v
			else:
				if val=="": return
				for c in self.cards:
					self.flair.addUndo(
						self.flair.setWhatUndo(c, what, val))
				return

	# ----------------------------------------------------------------------
	def edit(self, listbox):
		var, value = listbox.get(ACTIVE)
		if var == "type":
			return None
		elif var == "option":
			return tkExtra.InPlaceList(listbox.lists[1], values=self.option)
		elif var == "anchor":
			return tkExtra.InPlaceList(listbox.lists[1], values=_ANCHOR_list)
		elif var == "falloff":
			return tkExtra.InPlaceList(listbox.lists[1], values=_FALLOFF_list)
		elif var == "color":
			edit = tkExtra.InPlaceColor(listbox.lists[1])
			if edit.value is not None:	# convert to int
				try:
					edit.value = int(str(edit.value)[1:], 16)
				except:
					pass
			return edit
		elif var == "file":
			edit = tkExtra.InPlaceFile(listbox.lists[1],
				title="Load file",
				filetypes=(("STL files","*.stl"),("All","*")),
				save=False)
			if edit.value is not None:	# convert to int
				edit.value = self.flair.project.relativePath(edit.value)
			return edit
		else:
			return Properties.edit(self, listbox)

	# ----------------------------------------------------------------------
	def split(self):
		return int(self.listbox.index(ACTIVE)) >= 7

#===============================================================================
class MaterialProperties(Properties):
	def __init__(self, editor, cards):
		Properties.__init__(self, editor, cards, visible=True)
		self.items.extend(["color","alpha","specular","shine","fuzz","ior"])
		self.items2.extend(["global"])

	# ----------------------------------------------------------------------
	def populateList(self):
		Properties.populateList(self)
		self.setColor()

	# ----------------------------------------------------------------------
	def setColor(self):
		color = self.flair.project.getMaterial(self.cards[0].name())["color"]
		self.listbox.lists[1].itemconfig(2, background=color)

	# ----------------------------------------------------------------------
	def value(self, var=None, val=None):
		if var == "name":
			return Properties.value(self,var,val)

		elif var in ("color", "alpha", "specular", "shine", "fuzz", "ior"):
			project = self.flair.project
			if val is None:
				v = project.getMaterial(self.cards[0].name())[var]
				for c in self.cards:
					if v != project.getMaterial(c.name())[var]: return ""
				return v
			else:
				for c in self.cards:
					project.getMaterial(c.name())[var] = val
				self.editor.loadProjectColors()
				self.editor.setInputModified()

		elif var=="global":
			project = self.flair.project
			if val is None:
				v = project.getMaterial(self.cards[0].name()).isGlobal()
				for c in self.cards:
					if v != project.getMaterial(c.name()).isGlobal(): return ""
				if v: return _ON
				else: return _OFF
			else:
				for c in self.cards:
					if val: project.material2Global(c.name())
					else:   project.material2Local(c.name())

		else:
			return Properties.value(self, var, val)

	# ----------------------------------------------------------------------
	def edit(self, listbox):
		var, value = listbox.get(ACTIVE)
		if var == "color":
			edit = tkExtra.InPlaceColor(listbox.lists[1])
			return edit
		else:
			return Properties.edit(self, listbox)

#===============================================================================
# Body Menu
#===============================================================================
class _BodyMenuButton(FlairRibbon.FlairMenu):
	def createMenu(self):
		if self.page.project is None: return
		menu = Menu(self.page.ribbon, tearoff=1)
		tkFlair.bodiesMenu(menu, self.page.viewers.newBody)
		return menu

#===============================================================================
class _ObjectMenuButton(FlairRibbon.FlairMenu):
	def createMenu(self):
		if self.page.project is None: return
		menu = Menu(self.page.ribbon, tearoff=1)
		self.page.viewers._addObjectMenu(menu)
		return menu

#===============================================================================
# Layer Menu
#===============================================================================
class _LayerMenuButton(Ribbon.MenuButton):
	def __init__(self, master, page, **kw):
		Ribbon.MenuButton.__init__(self, master, None, **kw)
		self.page = page

	#----------------------------------------------------------------------
	def showMenu(self):
		if self.page.project is None: return
		menu = Menu(self, tearoff=0)
		for i,name in enumerate(self.page.project.geometry.layersList()):
			menu.add_command(label=name,
				command=lambda v=self.page.viewers,i=i: v.setLayer(i))
		self._showMenu(menu)

#===============================================================================
# Geometry Editor Empty
#===============================================================================
class GeometryEditorEmpty(FlairRibbon.FlairPage):

	_name_ = "Geometry"
	_icon_ = "geoedit"

	"""Graphical geometry editor"""
	def createPage(self, scroll=False):
		FlairRibbon.FlairPage.createPage(self, True)
		fn=os.path.join(tkFlair.prgDir,"icons","grab-geometry-page.gif")
		self.img = PhotoImage(file=fn)
		l = Label(self.frame,
			image=self.img,
			text="Please install the geoviewer module",
			foreground="Darkred",
			background="LightYellow",
			compound=BOTTOM)
		l.pack(expand=YES, fill=BOTH)

	# ----------------------------------------------------------------------
	# Mockup methods
	# ----------------------------------------------------------------------
	def _empty_(self, *args): pass

	loadProject = _empty_
	refreshUndo = _empty_
	loadState   = _empty_
	saveState   = _empty_
	activateIfNeeded = _empty_

#===============================================================================
# Tool group menu
#===============================================================================
class _ToolGroup(FlairRibbon.FlairMenuGroup):
	#----------------------------------------------------------------------
	def createMenu(self):
		menu = Menu(self, tearoff=0)

		menu.add_radiobutton(
				label="Correct",
				image=tkFlair.icons["edit"], compound=LEFT,
				accelerator="E",
				variable=self.page._mouseAction,
				value=ACTION_EDIT)
#		tkExtra.Balloon.set(b, "Correct region by removing overlaps")

		menu.add_radiobutton(
				label="Paint",
				image=tkFlair.icons["paint"], compound=LEFT,
				state=DISABLED,
#				accelerator="P",
				variable=self.page._mouseAction,
				value=ACTION_PAINT)
		#tkExtra.Balloon.set(b, "Set region properties (materials, biasing, cuts) by painting [p]")

		menu.add_radiobutton(
				label="Center",
				image=tkFlair.icons["center"], compound=LEFT,
				variable=self.page._mouseAction,
				accelerator="C",
				value=ACTION_VIEW_CENTER)

#		menu.add_radiobutton(
#				label="Pencil",
#				image=tkFlair.icons["pencil"], compound=LEFT,
#				state=DISABLED,
#				accelerator="P",
#				variable=self.page._mouseAction,
#				value=ACTION_PEN)

		menu.add_checkbutton(
				label="Snap",
				image=tkFlair.icons["snap"], compound=LEFT,
				variable = self.page.snap2Grid)
		return menu

#===============================================================================
# View group menu
#===============================================================================
class _ViewGroup(FlairRibbon.FlairMenuGroup):
	#----------------------------------------------------------------------
	def createMenu(self):
		menu = Menu(self, tearoff=0)

		submenu = Menu(menu)
		menu.add_cascade(label="Export", menu=submenu,
				compound=LEFT, image=tkFlair.icons["save"])
		for i in range(len(GeometryViewer.VIEWNAME)):
			submenu.add_command(label=GeometryViewer.VIEWNAME[i],
					foreground="Dark"+(GeometryViewer.VIEWNAME[i]),
					activeforeground="Dark"+(GeometryViewer.VIEWNAME[i]),
					underline=0,
					command=lambda p=self.page,v=i: p.export(v),
					compound=LEFT, image=tkFlair.icons["save"],
					accelerator="Ctrl-P")

		submenu = Menu(menu)
		menu.add_cascade(label="Notes", menu=submenu,
				compound=LEFT, image=tkFlair.icons["document"])
		for i in range(len(GeometryViewer.VIEWNAME)):
			submenu.add_command(label=GeometryViewer.VIEWNAME[i],
					foreground="Dark"+(GeometryViewer.VIEWNAME[i]),
					activeforeground="Dark"+(GeometryViewer.VIEWNAME[i]),
					underline=0,
					command=lambda p=self.page,v=i: p.toNotes(v),
					compound=LEFT, image=tkFlair.icons["document"],
					accelerator="Ctrl-N")

		return menu

#===============================================================================
# Geometry Editor Window
#===============================================================================
class GeometryEditor(FlairRibbon.FlairPage):
	"""Graphical geometry editor"""

	_name_ = "Geometry"
	_icon_ = "geoedit"

	#----------------------------------------------------------------------
	def init(self):
		# project variables
		self._inFocus     = False
		self._projectName = None	# XXX FIXME should be eliminated
						# check the project id
						# since now we have a new project
						# every time we load one
		self.editObject   = None	# Object currently being edited
		self._undoAfter   = None	# Undo after

		self._loadStatePending = True	# Make a refresh on first load
		self._loadingProject   = False	# Do not send any update to flair while loading project
		self._refreshOn        = ("card","voxel","run","plot")

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Tools ==============
		group = _ToolGroup(self.ribbon, self, "Tools")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ======== Select-Tool ===========
		col,row = 0,0
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["select32"],
				text="Select",
				compound=TOP,
				variable=self._mouseAction,
				value=ACTION_SELECT,
				padx=5)
		b.grid(row=row, rowspan=3, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Select [s]")

		# ======== Pan/Rotate/Select ========
		col,row=1,0
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["pan"],
				text="Pan",
				compound=LEFT,
				anchor=W,
				variable=self._mouseAction,
				value=ACTION_PAN)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Pan viewport [x]")

		# ---
		col,row=1,1
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["rotate"],
				text="Orbit",
				compound=LEFT,
				anchor=W,
				variable=self._mouseAction,
				value=ACTION_ORBIT)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Orbit viewport [t]")

		# ---
		col,row=1,2
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["target"],
				text="Info",
				compound=LEFT,
				anchor=W,
				variable=self._mouseAction,
				value=ACTION_INFO)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Display body and region information on cursor location [i]")

		# ======== Zoom ========
		col,row = 2,0
		b = Ribbon.LabelRadiobutton(group.frame, image=tkFlair.icons["zoom_in"],
				#text="Zoom In",
				#compound=LEFT,
				#anchor=W,
				variable=self._mouseAction,
				value=ACTION_ZOOM_IN)
		b.grid(row=row, column=col, padx=5, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Drag a rectangle to zoom in")

		# ----
		col,row = 2,1
		b = Ribbon.LabelRadiobutton(group.frame, image=tkFlair.icons["zoom_out"],
				#text="Zoom Out",
				#compound=LEFT,
				#anchor=W,
				variable=self._mouseAction,
				value=ACTION_ZOOM_OUT)
		b.grid(row=row, column=col, padx=5, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Drag a rectangle to zoom out")

		# ----
		col,row = 2,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["zoom_on"],
				#text="Zoom On",
				#compound=LEFT,
				#anchor=W,
				command=self.zoomOnSelected,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=5, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Zoom on selected objects [f]")

		# ========== Add ==============
		group = Ribbon.LabelGroup(self.ribbon, "Add")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ----
		col,row = 0,0
		b = _BodyMenuButton(group.frame, self,
				image=tkFlair.icons["cone32"],
				text="Body",
				compound=TOP,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add Body [b]. Repeat body [B]")

		# ----
		col,row = 1,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["zone32"],
				text="Zone",
				compound=TOP,
				command=self.viewers.setActionZone,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add Zone to select region, else create a new region to add zone [d]")

		# ---
		col,row=2,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["REGION"],
				text="Region",
				compound=LEFT,
				anchor=W,
				command=self.addRegion,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Add REGION card [R]")

		# ----
		col,row=2,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["ROT-DEFI"],
				text="ROT-DEFI",
				compound=LEFT,
				anchor=W,
				command=self.addRotdefi,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add a ROT-DEFI card")

		# ----
		col,row=2,2
		b = _ObjectMenuButton(group.frame, self,
				image=tkFlair.icons["pin"],
				text="Object " + Unicode.BLACK_DOWN_POINTING_TRIANGLE,
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add a flair object [space]")

		# ----
		col,row=4,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Delete",
				compound=LEFT,
				anchor=W,
				command=self.delete,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Delete cards and remove their references [Del] or [BackSpace]")

		# ----
		col,row=4,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clone"],
				text="Clone",
				compound=LEFT,
				anchor=W,
				command=self.clone,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Clone cards and all associated cards to them [Ctrl-D]")

		# ========== Selection ==============
		group = Ribbon.LabelGroup(self.ribbon, "Selection")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ----
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["toggle32"],
				command=self.toggleVisibility,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, rowspan=2, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Toggle Visibility [v]")

		col,row = 0,2
		menulist = [	("Show",      "enable",   self.setVisibility),
				("Hide",      "disable",  self.unsetVisibility),
				None,
				("Show All",  "enable",   self.setAllVisibility),
				("Hide All",  "disable",  self.unsetAllVisibility)]
		b = Ribbon.MenuButton(group.frame, menulist,
				text="Visibility",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Set/Unset visibility")

		# ---
		col,row=1,1
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["wireframe"],
				command=self.toggle3D,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Toggle wireframe [#]")

		menulist = [	("Show",     "wireframe",     self.set3D),
				("Hide",     "wireframe_off", self.unset3D),
				None,
				("Hide All", "wireframe_off", self.unset3DAll)]
		b = Ribbon.MenuButton(group.frame, menulist,
				image=tkFlair.icons["triangle_down"],
				text="Wireframe",
				compound=RIGHT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+1, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Set/Unset wireframe")

		# ---
		col,row=1,0
		b = Ribbon.LabelCheckbutton(group.frame,
				image=tkFlair.icons["lock"],
				text="Lock",
				compound=LEFT,
				anchor=W,
				variable=self.lockEditing)
				#value=ACTION_SELECT,
		b.grid(row=row, column=col, columnspan=3, padx=0, pady=0, sticky=NSEW)
#		b = Ribbon.LabelButton(group.frame,
#				image=tkFlair.icons["lock"],
#				command=self.toggleLock,
#				background=Ribbon._BACKGROUND)
#		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
#		tkExtra.Balloon.set(b, "Toggle editing lock [l]")
#		menulist = [	("Lock Bodies",   "lock",     self.lock),
#				("Unlock bodies", "unlock",   self.unlock),
#				None,
#				("Lock All bodies",  "lock",   self.lockAll),
#				("Unlock All bodies","unlock", self.unlockAll)]
#		b = Ribbon.MenuButton(group.frame, menulist,
#				image=tkFlair.icons["triangle_down"],
#				text="Lock",
#				compound=RIGHT,
#				anchor=W,
#				background=Ribbon._BACKGROUND)
#		b.grid(row=row, column=col+1, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Lock/Unlock graphical editing")

		# ---
		col,row=1,2
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["freeze"],
				command=self.toggleFreeze,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Freeze selection [F]")

		menulist = [	("Freeze",     "freeze",     self.freeze),
				("Unfreeze",   "heat",       self.unfreeze),
				None,
				("Unfreeze All","heat",      self.unfreezeAll)]
		b = Ribbon.MenuButton(group.frame, menulist,
				image=tkFlair.icons["triangle_down"],
				text="Selection",
				compound=RIGHT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+1, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Freeze/Unfreeze selection")

		# ========== View ===========
		group = _ViewGroup(self.ribbon, self, "View")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		group.grid3rows()
		group.frame.grid_columnconfigure(1, weight=1)
		group.frame.grid_columnconfigure(3, weight=1)

		# ===
		col,row=0,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["run"],
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=3, pady=0, sticky=E)
		tkExtra.Balloon.set(b, "Select active run settings")

		self.runCombo = Ribbon.LabelCombobox(group.frame,
					command=self.runChanged,
					width=30)
		self.runCombo.grid(row=row, column=col+1, columnspan=4, padx=0, pady=0, sticky=EW)

		# ---
		col,row=0,1
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["layers"],
				command=self.viewers.setLayer,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Rotate layer on all viewports")

		b = _LayerMenuButton(group.frame, self,
				image=tkFlair.icons["triangle_down"],
				text="Layer",
				compound=RIGHT,
				anchor=W,
#				width=8,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+1, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Select layer for all viewports")

		# ---
		col,row=2,1
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["layout"],
				command=self.viewLayout,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Select windows layout [w,W,#]")

		menulist = []
		for name in GeometryViewer.VIEWORDER:
			icon = "layout_"+name
			menulist.append((name, icon,
				lambda s=self,n=name : s.viewLayout(n)))

		b = Ribbon.MenuButton(group.frame, menulist,
				image=tkFlair.icons["triangle_down"],
				text="Layout",
				compound=RIGHT,
				anchor=W,
#				width=8,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+1, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Select windows layout")

		# ---
		col,row=0,2
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["refresh"],
				text="Reload",
				compound=LEFT,
				anchor=W,
				command=self.reload,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Reload input file")

		# ---
		col,row=2,2
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["reset"],
				text="Reset",
				compound=LEFT,
				anchor=W,
				command=self.resetEditor,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Reset all viewports to system default")

		# ---
#		menulist = []
#		for i in range(len(GeometryViewer.VIEWNAME)):
#			name = GeometryViewer.VIEWNAME[i]
#			cmd  = lambda p=self.page,v=i: p.export(v)
#			menulist.append((name,"save",cmd))
#
#		col,row=2,2
#		b = Ribbon.MenuButton(group.frame, menulist,
#				image=tkFlair.icons["save"],
#				text="Viewport",
#				compound=LEFT,
#				anchor=W,
#				background=Ribbon._BACKGROUND)
#		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=EW)
#		tkExtra.Balloon.set(b, "Save view port into image")

	#----------------------------------------------------------------------
	# Ribbon method to create context groups
	#----------------------------------------------------------------------
	def createContextGroups(self):
		# ========== Body ===========
		group = Ribbon.LabelGroup(self.ribbon, "Body")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		self.addContextGroup("body", group)
		group.grid3rows()

		# ---
		col,row=0,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["axes"],
				text="Transform",
				compound=LEFT,
				anchor=W,
				command=self.transformGeometry,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Transform bodies (Ctrl-T)")

		# ---
		col,row=0,1
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["expand_parenthesis"],
				text="Expand",
				compound=LEFT,
				anchor=W,
				command=self.expandBody,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Expand macro bodies to infinite ones")

		# ---
		col,row=0,2
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["ROT-DEFI"],
				text="ROT-DEFI",
				compound=LEFT,
				anchor=W,
				command=self.gotoRotdefi,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Show body ROT-DEFI transformation card")

		# ========== Region ===========
		group = Ribbon.LabelGroup(self.ribbon, "Region")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		self.addContextGroup("region", group)
		group.grid3rows()

		# ---
		col,row=0,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["MATERIAL"],
				text="Material",
				compound=LEFT,
				anchor=W,
				command=self.gotoMaterial,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Edit properties of selected material")

		# ---
		col,row=0,1
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["RPP"],
				text="Volume",
				compound=LEFT,
				state=DISABLED,
				anchor=W,
				command=self.regionVolume,
				background=Ribbon._BACKGROUND)
		if Input._developer: b["state"] = NORMAL
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Calculate volume of selected regions")

		# ---
		col,row=0,2
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["ROT-DEFI"],
				text="ROT-DEFI",
				compound=LEFT,
				anchor=W,
				command=self.gotoRotdefi,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Show lattice region ROT-DEFI transformation card")

		# ---
		col,row=1,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["expand_parenthesis"],
				text="Expand",
				compound=LEFT,
				anchor=W,
				command=self.expandRegions,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Expand parenthesis if any. WARNING can take long")

		# ---
		col,row=1,1
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["light"],
				command=self.toggleAlpha,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Toggle region transparency [T]")

		menulist = []
		for a in range(0,101,10):
			menulist.append(("%02d%%"%(a), None,
				lambda s=self,a=a : s.toggleAlpha(a)))

		b = Ribbon.MenuButton(group.frame, menulist,
				image=tkFlair.icons["triangle_down"],
				text="Transparency",
				compound=RIGHT,
				anchor=W,
#				width=8,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+1, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Set region transparency")

		# ========== Mesh ===========
		group = Ribbon.LabelGroup(self.ribbon, "Mesh")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		self.addContextGroup("mesh", group)
		group.grid3rows()

		# ---
		col,row=0,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["TRC"],
				text="Bodies",
				compound=LEFT,
				anchor=W,
				command=self.mesh2Bodies,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Create bodies out of the mesh")

	#----------------------------------------------------------------------
	# GeometryEditor toplevel class
	#----------------------------------------------------------------------
	def createPage(self, scroll=False):
		FlairRibbon.FlairPage.createPage(self, False)

		# check buttons
		self.showBodies   = BooleanVar();	self.showBodies.set(True)
		self.showRegions  = BooleanVar();	self.showRegions.set(True)
		self.showRotdefi  = BooleanVar();	self.showRotdefi.set(True)
		self.showMaterial = BooleanVar();	self.showMaterial.set(False)
		self.showObjects  = BooleanVar();	self.showObjects.set(True)
		self.showVisible  = BooleanVar()
		self.snap2Grid    = BooleanVar();	self.snap2Grid.set(True)
		self.lockEditing  = BooleanVar();	self.lockEditing.set(False)

		self._lastBodyName   = "body"
		self._lastRegionName = "REG"
		self._lastMaterial   = "VACUUM"

		self._mouseAction = IntVar()
		self._mouseAction.set(ACTION_SELECT)
		self._mouseAction.trace("w",self.mouseActionChanged)

		# ---------- create layout -------------
		self.splitter = tkExtra.HSplitter(self.frame,
					tkFlair.getInt(tkFlair._PAGE_SECTION,
						self.name+".split", 250),
						True)
		self.splitter.pack(side=BOTTOM, expand=YES, fill=BOTH)

		# =============== TabSet Geometry/Layers/Debug=================
		self._tabSet = tkExtra.TabPageSet(self.splitter.leftFrame(),
				pageNames=[("Geometry", tkFlair.icons["REGION"]),
					   ("Layers",   tkFlair.icons["layers"]),
					   ("Errors",   tkFlair.icons["debug" ])],
				hidetext=False)
		self._tabSet.pack(expand=YES, fill=BOTH)
		self._tabSet.bind("<<ChangePage>>", self.changePage)

		# Fill Layers
		tabframe = self._tabSet['Layers']
		self._layers = GeometryLayers.GeometryLayersFrame(tabframe, self.flair)
		self._layers.bind("<<LayersUpdate>>", self.updateLayers)
		self._layers.pack(expand=YES, fill=BOTH)

		# Fill Errors
		tabframe = self._tabSet['Errors']
		self._errors = GeometryExtra.ErrorsFrame(tabframe, self)
		self._errors.pack(expand=YES, fill=BOTH)

		# Fill Geometry
		tabframe = self._tabSet['Geometry']
		self._tabSet.changePage()

		# --------------------------------------------------------------
		frame = Frame(tabframe)
		frame.pack(fill=X)

		# 3state buttons!
		b = Checkbutton(frame,image=tkFlair.icons["TRC"],
				command=self.objListFill,
				variable=self.showBodies,
				indicatoron=False,
				highlightthickness=0,
				selectcolor="White",
				padx=0, pady=0)
		b.bind("<Control-Button-1>", self.onlyBodies)
		b.pack(side=LEFT, padx=1)
		tkExtra.Balloon.set(b, "Filter Bodies. <Ctrl-Click> to select only Bodies")

		b = Checkbutton(frame,image=tkFlair.icons["REGION"],
				command=self.objListFill,
				variable=self.showRegions,
				indicatoron=False,
				highlightthickness=0,
				selectcolor="White",
				padx=0, pady=0)
		b.bind("<Control-Button-1>", self.onlyRegions)
		b.pack(side=LEFT, padx=1)
		tkExtra.Balloon.set(b, "Filter Regions. <Ctrl-Click> to select only Regions")

		b = Checkbutton(frame,image=tkFlair.icons["ROT-DEFI"],
				command=self.objListFill,
				variable=self.showRotdefi ,
				indicatoron=False,
				highlightthickness=0,
				selectcolor="White",
				padx=2, pady=0)
		b.bind("<Control-Button-1>", self.onlyRotdefi)
		b.pack(side=LEFT, padx=1)
		tkExtra.Balloon.set(b, "Filter Transformations. <Ctrl-Click> to select only Transformations")

		b = Checkbutton(frame,image=tkFlair.icons["MATERIAL"],
				command=self.objListFill,
				variable=self.showMaterial,
				indicatoron=False,
				highlightthickness=0,
				selectcolor="White",
				padx=0, pady=0)
		b.bind("<Control-Button-1>", self.onlyMaterials)
		b.pack(side=LEFT, padx=1)
		tkExtra.Balloon.set(b, "Filter Materials. <Ctrl-Click> to select only Materials")

		b = Checkbutton(frame,image=tkFlair.icons["pin"],
				command=self.objListFill,
				variable=self.showObjects,
				indicatoron=False,
				highlightthickness=0,
				selectcolor="White",
				padx=0, pady=0)
		b.bind("<Control-Button-1>", self.onlyObjects)
		b.pack(side=LEFT, padx=1)
		tkExtra.Balloon.set(b, "Filter Objects. <Ctrl-Click> to select only Objects")

		b = Checkbutton(frame,image=tkFlair.icons["select"],
				command=self.objListFill,
				variable=self.showVisible,
				indicatoron=False,
				highlightthickness=0,
				selectcolor="White",
				padx=2, pady=0)
		b.pack(side=LEFT, padx=1)
		tkExtra.Balloon.set(b, "Show Selected & Visible items")

		b = Button(frame,image=tkFlair.icons["clean"],
				command=self.selectionClear,
				padx=2, pady=0)
		b.pack(side=RIGHT, padx=1)
		tkExtra.Balloon.set(b, "Clear selected items")

		# --------------------------------------------------------------
		frame = Frame(tabframe)
		frame.pack(fill=X)
		self.filter = tkExtra.LabelEntry(frame, label="Filter", labelcolor="DarkGray",
					background="White", width=7)
		self.filter.pack(side=LEFT, expand=YES, fill=X)
		self.filter.bind("<Return>",   self.objListFill)
		self.filter.bind("<KP_Enter>", self.objListFill)
		tkExtra.Balloon.set(self.filter, "Filter items")
		b = Button(frame,text="X",command=self.filterClear, padx=1,pady=1)
		b.pack(side=RIGHT)
		tkExtra.Balloon.set(b, "Clear filter")

		# --------------------------------------------------------------
		vsplit = tkExtra.VSplitter(tabframe, 0.5)
		vsplit.pack(expand=YES,fill=BOTH)

		# ------------------- Objects ----------------------------------
		self.objList = tkExtra.MultiListbox(vsplit.topFrame(),
					(("Type",  8, None),
					 ("Name", 10, None)),
					 header=True,
					 #font=FONT_LIST,
					 stretch="last",
					 background="White")
		self.objListCard = []	# Contains the card of every item in the list
		self.objList.pack(expand=YES, fill=BOTH)
		self.objList.bind("<<ListboxSelect>>",	     self.objListSelect)
		self.objList.bind("<<ListboxSort>>",	     self.objListSort)
		self.objList.bindList("<Delete>",	     self.objListDelete)
		self.objList.bindList("<Button-1>",          self.objListClick)
		self.objList.bindList("<ButtonRelease-1>",   self.objListRelease)
		self.objList.bindList("<Double-1>",	     self.objListName)
		self.objList.bindList("<Return>",	     self.objListName)
		self.objList.bindList("<KP_Enter>",	     self.objListName)
		self.objList.bindList("<Button-3>",	     self.objListPopup)
		self.objList.bindList("<F11>",		     self.toggleSplitter)

		self.objList.bindList("<Control-Key-space>", self.objListPopup)
		self.objList.bindList("<Control-Key-c>",     self.objListCopy)
		self.objList.bindList("<Control-Key-C>",     self.objListCopy)
#		self.objList.bindList("<Control-Key-d>",     self.objListClone)
#		self.objList.bindList("<Control-Key-D>",     self.objListClone)
		self.objList.bindList("<F2>",                self.objListName)
		self.objList.bindList("<Control-Key-e>",     self.objListName)
		self.objList.bindList("<Control-Key-v>",     self.objListPaste)
		self.objList.bindList("<Control-Key-V>",     self.objListPasteSpecial)
		self.objList.bindList("<Control-Key-x>",     self.objListCut)

		# -------------------- Properties/Attributes -------------------
		# -------- TabFrame ----
		self.propTab = tkExtra.TabPageSet(vsplit.bottomFrame(),
			pageNames=[("Properties",tkFlair.icons["properties"]),
				   ("Attributes",tkFlair.icons["attributes"])])
		self.propTab.pack(side=BOTTOM, expand=YES, fill=BOTH)

		tabframe = self.propTab[_PROPTABS[0]]
		self.propList = tkExtra.MultiListbox(tabframe,
					(("Property", 8, None),
					 ("Value",   10, None)),
					 #font=FONT_LIST,
					 header=False,
					 stretch="last",
					 background="White")
		self.propList.sortAssist = None
		self.propList.pack(side=BOTTOM, expand=YES, fill=BOTH)
		self.propTab.changePage()

		self.propList.bind("<<ListboxSelect>>",		self.propListSelect)
		self.propList.bindList("<Delete>",		self.propListDelete)
		self.propList.bindList("<BackSpace>",		self.propListDelete)
		self.propList.bindList("<F2>",			self.propListEditEvent)
		self.propList.bindList("<Return>",		self.propListEditEvent)
		self.propList.bindList("<KP_Enter>",		self.propListEditEvent)
		self.propList.bindList("<F11>",		        self.toggleSplitter)
		self.propList.bindList("<ButtonRelease-1>",	self.propListRelease)
		self.propList.bindList("<Button-3>",		self.propListPopup)
		self.propList.bindList("<Control-Key-space>",	self.propListPopup)
		self.propList.bindList("<Control-Key-c>",	self.propListCopy)
		self.propList.bindList("<Control-Key-x>",	self.propListCut)
		self.propList.bindList("<Control-Key-v>",	self.propListPaste)
		self.propList.bindList("<Control-Key-Up>",	self.propListMoveUp)
		self.propList.bindList("<Control-Key-Down>",	self.propListMoveDown)

		tabframe = self.propTab[_PROPTABS[1]]
		self.attrList = tkExtra.MultiListbox(tabframe,
					(("Attribute",  8, None),
					 ("Value",     15, None)),
					 #font=FONT_LIST,
					 header=False,
					 stretch="last",
					 background="White")
		self.attrList.sortAssist = None
		self.attrList.pack(side=BOTTOM, expand=YES, fill=BOTH)

		#self.attrList.bind("<<ListboxSelect>>",		self.attrListSelect)
		self.attrList.bindList("<Delete>",		self.attrListDelete)
		self.attrList.bindList("<BackSpace>",		self.attrListDelete)
		self.attrList.bindList("<F2>",			self.attrListEdit)
		self.attrList.bindList("<Return>",		self.attrListEdit)
		self.attrList.bindList("<KP_Enter>",		self.attrListEdit)
		self.attrList.bindList("<F11>",		        self.toggleSplitter)
		self.attrList.bindList("<ButtonRelease-1>",	self.attrListRelease)
		self.attrList.bindList("<Button-3>",		self.attrListPopup)
		self.attrList.bindList("<Control-Key-space>",	self.attrListPopup)
		#self.attrList.bindList("<Control-Key-c>",	self.attrListCopy)
		#self.attrList.bindList("<Control-Key-x>",	self.attrListCut)
		#self.attrList.bindList("<Control-Key-x>",	self.attrListCut)
		#self.attrList.bindList("<Control-Key-v>",	self.attrListPaste)
		#self.attrList.bindList("<Control-Key-Up>",	self.attrListMoveUp)
		#self.attrList.bindList("<Control-Key-Down>",	self.attrListMoveDown)

		# --------------------------------------------------------------
		self.viewers = GeometryViewer.Geometry4Frame(
					self.splitter.rightFrame(), self)
		self.viewers.pack(side=RIGHT, expand=YES, fill=BOTH)
		self.viewers.focus_set()
		self.viewers.bind("<<Errors>>", self.updateErrors)

	#----------------------------------------------------------------------
	def createStatus(self):
		FlairRibbon.FlairPage.createStatus(self)
		Label(self.status, text="x:", foreground="DarkGreen", pady=1).pack(side=LEFT)
		self.xlabel = Label(self.status, width=14, anchor=W, pady=1)
		self.xlabel.pack(side=LEFT, fill=BOTH, padx=1, pady=1)

		Label(self.status, text=" y:", foreground="DarkGreen", pady=1).pack(side=LEFT)
		self.ylabel = Label(self.status, width=14, anchor=W, pady=1)
		self.ylabel.pack(side=LEFT, fill=BOTH, padx=1, pady=1)

		Label(self.status, text=" z:", foreground="DarkGreen", pady=1).pack(side=LEFT)
		self.zlabel = Label(self.status, width=14, anchor=W, pady=1)
		self.zlabel.pack(side=LEFT, fill=BOTH, padx=1, pady=1)

		return self.status

	# ----------------------------------------------------------------------
	def vsplit(self):
		if not self.created(): return None
		return self.viewers.vsplit.split

	# ----------------------------------------------------------------------
	def tsplit(self):
		if not self.created(): return None
		return self.viewers.tsplit.split

	# ----------------------------------------------------------------------
	def bsplit(self):
		if not self.created(): return None
		return self.viewers.bsplit.split

	# ----------------------------------------------------------------------
	def resetEditor(self):
		self.viewers.resetViewports()
		self.viewers.resetOrigin()
		self.viewers.layout("2x2")
		self.reload()

	# ----------------------------------------------------------------------
	def reset(self):
		self._loadStatePending = True	# Force refreshing after showing the window

	# ----------------------------------------------------------------------
	def beforeDock(self):
		self.reset()

	# ----------------------------------------------------------------------
	def refresh(self, event=None):
		FlairRibbon.FlairPage.refresh(self)
		if not self.created(): return
		self.project.input.refresh()
		self.show()
		self.redraw()

	# ----------------------------------------------------------------------
	def reload(self, event=None):
		self.saveState()
		self.flair.checkInputFile(self.page)
		self.refresh()

	# ----------------------------------------------------------------------
	def show(self):
		"""Show project on viewports"""
		if self.project is None: return
		# FIXME should be eliminated
		if self._projectName != self.project.name:
			for card in self.project.input.allcards():
				if card.ignore(): continue
				card[SELECT] = card.get(SELECT,0) & BIT_UNSELECT
			self._projectName = self.project.name

		self.loadProject()
		if self._loadStatePending:
			self._loadStatePending = False
			self.loadState()
		self.selectShow(False)

	# ----------------------------------------------------------------------
	# Load project and refresh information
	# ----------------------------------------------------------------------
	def loadProject(self, all=False):
		if GeometryViewer._PROFILE:
			say("\n*** GeometryEditor.loadProject started")
			start = time.time()
		self._loadingProject = True
		cursor = self.page["cursor"]
		self.page["cursor"] = "watch"
		self.fillRuns()
		self.page["cursor"] = cursor
		self._loadingProject = False
		if GeometryViewer._PROFILE:
			say("*** GeometryEditor.loadProject ended", time.time()-start,"s")

	# ----------------------------------------------------------------------
	# Fill the runCombo with active runs from the flair project
	# ----------------------------------------------------------------------
	def fillRuns(self):
		activeRun = self.runCombo.get()
		self.runCombo.clear()
		found = False
		for run in self.project.runs:
			if run.parent != "": continue
			runname = run.getInputName2()
			self.runCombo.insert(END, runname)
			if activeRun == runname:
				found = True

		# Set the run in the combo to actually load the project
		if found:
			self.runCombo.set(activeRun)
		else:
			self.runCombo.set(self.runCombo.get(0))

	# ----------------------------------------------------------------------
	def runChanged(self):
		"""Called by runCombo when run is changed"""
		run = self.project.getRunByName(self.runCombo.get())
		tkExtra.Balloon.set(self.runCombo, "Active Run: %s"%(self.runCombo.get()))
		self.viewers.loadProject(self.project, run)
		self.updateErrors()
		self.objListFill()
		self.redraw()

	# ----------------------------------------------------------------------
	def export(self, view):
		"""Export viewport to image"""
		self.viewers.frames[view].export()

	# ----------------------------------------------------------------------
	def toNotes(self, view):
		"""Export viewport to notes"""
		self.viewers.frames[view].toNotes()

	# ----------------------------------------------------------------------
	# Save&Load current state and refresh flair window if any
	# ----------------------------------------------------------------------
	def loadState(self):
		if not self.created(): return
		self.viewers.loadState()

	# ----------------------------------------------------------------------
	def saveState(self):
		if not self.created(): return
		self.viewers.saveState()

	# ----------------------------------------------------------------------
	def get(self):
		self.saveState()

	# ----------------------------------------------------------------------
	# Focus In event
	# ----------------------------------------------------------------------
	def _focusIn(self, event=None):
		"""Handle the <FocusIn> event. Reload project if necessary"""
		if self._inFocus: return
		FlairRibbon.FlairPage._focusIn(self,event)

		# FocusIn is generated for all sub-windows
		if event and event.widget is not self or self.project is None: return
		self._inFocus = True

		try:	# FIXME should become absolute paths in the future!
			# Always change the directory to the project one
			os.chdir(self.project.dir)
		except:
			pass

		self.flair.checkInputFile(self)

		self._inFocus = False

	# ----------------------------------------------------------------------
#	def _focusOut(self, event):
#		"""Handle the <FocusOut> event"""
#		if self is not event.widget: return
#		self.saveState()

	# ----------------------------------------------------------------------
	def cut(self, event=None):
		if event is None:
			self.objListCut(event)
		elif event.widget is self.propList or \
		     self.editType() == Input.Card.REGION and \
		     self.propList.curselection():
			self.propListCut(event)
#		elif event.widget is self.attrList:
#			self.attrListCut(event)
		elif isinstance(event.widget, Entry) or \
		     isinstance(event.widget, Text):
			FlairRibbon.FlairPage.cut(self, event)
		else:
			self.objListCut(event)

	# ----------------------------------------------------------------------
	def copy(self, event=None):
		if event is None:
			self.objListCopy(event)
		elif event.widget is self.propList:
			self.propListCopy(event)
#		elif event.widget is self.attrList:
#			self.attrListCopy(event)
		elif isinstance(event.widget, Entry) or \
		     isinstance(event.widget, Text):
			FlairRibbon.FlairPage.copy(self, event)
		else:
			self.objListCopy(event)

	# ----------------------------------------------------------------------
	def paste(self, event=None):
		if event is None:
			self.objListPaste(event)
		elif event.widget is self.propList:
			self.propListPaste(event)
#		elif event.widget is self.attrList:
#			self.attrListPaste(event)
		elif isinstance(event.widget, Entry) or \
		     isinstance(event.widget, Text):
			FlairRibbon.FlairPage.paste(self, event)
		else:
			self.objListPaste(event)

	# ----------------------------------------------------------------------
	def pasteSpecial(self, event=None):
		return self.objListPasteSpecial(event)

	# ----------------------------------------------------------------------
	def clone(self, event=None):
		self.copy(event)
		self.pasteSpecial(event)
		return "break"

	# ----------------------------------------------------------------------
	def openProject(self, event=None):
		self.flair.lift(self)
		return self.flair.openProject()

	# ----------------------------------------------------------------------
	def saveProject(self, event=None, filename=None):
		"""Save current project"""
		self.saveState()
		self.flair.lift(self)
		self.flair.saveProject()
		try: self.lift(self.flair)
		except TclError: pass

	# ----------------------------------------------------------------------
	def saveProjectAs(self, event=None):
		"""Save current project with different filename"""
		self.saveState()
		self.flair.lift(self)
		return self.flair.saveProjectAs()

	# ----------------------------------------------------------------------
	def setInputModified(self, m=True):
		self.flair.setInputModified(m)

	# ----------------------------------------------------------------------
	def filterClear(self):
		"""Clear the active filter"""
		self.filter.delete(0,END)
		self.objListFill()

	# ----------------------------------------------------------------------
	def onlyBodies(self, event=None):
		self.showBodies.set(True)
		self.showRegions.set(False)
		self.showRotdefi.set(False)
		self.showMaterial.set(False)
		self.showObjects.set(False)
		self.objListFill()
		return "break"

	# ----------------------------------------------------------------------
	def onlyRegions(self, event=None):
		self.showBodies.set(False)
		self.showRegions.set(True)
		self.showRotdefi.set(False)
		self.showMaterial.set(False)
		self.showObjects.set(False)
		self.objListFill()
		return "break"

	# ----------------------------------------------------------------------
	def onlyRotdefi(self, event=None):
		self.showBodies.set(False)
		self.showRegions.set(False)
		self.showRotdefi.set(True)
		self.showMaterial.set(False)
		self.showObjects.set(False)
		self.objListFill()
		return "break"

	# ----------------------------------------------------------------------
	def onlyMaterials(self, event=None):
		self.showBodies.set(False)
		self.showRegions.set(False)
		self.showRotdefi.set(False)
		self.showMaterial.set(True)
		self.showObjects.set(False)
		self.objListFill()
		return "break"

	# ----------------------------------------------------------------------
	def onlyObjects(self, event=None):
		self.showBodies.set(False)
		self.showRegions.set(False)
		self.showRotdefi.set(False)
		self.showMaterial.set(False)
		self.showObjects.set(True)
		self.objListFill()
		return "break"

	# ----------------------------------------------------------------------
	def objListFill(self, event=None):
		"""Populate object list with requested items"""
		search = self.filter.get().upper()

		# empty the list
		self.objList.delete(0,END)
		del self.objListCard[:]

		inp = self.project.input
		visible = self.showVisible.get()

		# FIXME items with error highlight RED

		def searchCard(card, search):
			if card.tag.find(search)>=0: return True
			if card.name().upper().find(search)>=0: return True
			if card.comment().upper().find(search)>=0: return True
			if card.extra().upper().find(search)>=0: return True
			if card.tag=="REGION" and card["@material"].name().upper().find(search)>=0: return True
			for w in card.whats():
				if str(w).upper().find(search)>=0:
					return True
			return False

		# Bodies
		if self.showBodies.get():
			for card in inp.cardlist:
				if not card.isGeo(): continue
				if card.ignore(): continue
				if len(card.tag)!=3 and card.tag!="VOXELS": continue
				if card.tag=="END": continue
				if card.tag=="VOXELS":
					name = "VOXEL"
				else:
					name = card.name()
				if not search or searchCard(card, search):
					name = (card.tag,name)
					if visible:
						if card.get(SELECT,0)>0:
							self.objList.insert(END, name)
							self.objListCard.append(card)
					else:
						self.objList.insert(END, name)
						self.objListCard.append(card)

		# Regions
		if self.showRegions.get():
			for card in inp.cardsSorted("REGION"):
				if card.ignore() or card.name()=="&": continue
				if not search or searchCard(card, search):
					name = (card.tag,card.name())
					if visible:
						if card.get(SELECT,0)>0:
							self.objList.insert(END, name)
							self.objListCard.append(card)
					else:
						self.objList.insert(END, name)
						self.objListCard.append(card)

		# Transformations
		if self.showRotdefi.get():
			rotdefi = {}
			for card in inp.cardsSorted("ROT-DEFI"):
				if card.ignore(): continue
				if not search or searchCard(card, search):
					# FIXME for old inputs without names!!!
					name = card.name()
					name = (card.tag,name)
					if name in rotdefi: continue
					rotdefi[name] = card
					if visible:
						if card.get(SELECT,0)>0:
							self.objList.insert(END, name)
							self.objListCard.append(card)
					else:
						self.objList.insert(END, name)
						self.objListCard.append(card)

		# Materials
		if self.showMaterial.get():
			for card in inp.cardsSorted("MATERIAL"):
				#if card.ignore() or card.name()=="&": continue
				if not search or searchCard(card, search):
					name = (card.tag,card.name())
					if visible:
						if card.get(SELECT,0)>0:
							self.objList.insert(END, name)
							self.objListCard.append(card)
					else:
						self.objList.insert(END, name)
						self.objListCard.append(card)

		# Objects
		if self.showObjects.get():
			cards = []
			for tag in Input.FLAIR_TAGS:
				cards.extend(inp[tag])
			cards.extend(inp["BEAM"])
			cards.extend(inp["SPOTBEAM"])
			cards.sort(key=attrgetter("_pos"))

			for card in cards:
				if card.ignore(): continue
				if not search or searchCard(card, search):
					name = (card.tag,card.name())
					if visible:
						if card.get(SELECT,0)>0:
							self.objList.insert(END, name)
							self.objListCard.append(card)
					else:
						self.objList.insert(END, name)
						self.objListCard.append(card)

		# scan items to select
		last = None
		for i in range(self.objList.size()):
			card = self.objListCard[i]
			flag = card.get(SELECT,0)

			if card.invalid:
				color = ERROR_COLOR
			elif flag & BIT_LOCK:
				color = LOCK_COLOR
			elif flag & BIT_FREEZE:
				color = FREEZE_COLOR
			elif flag & BIT_VISIBLE:
				color = VISIBLE_COLOR
			elif flag & BIT_WIREFRAME:
				color = WIREFRAME_COLOR
			else:
				color = None
			if color is not None and color != self.objList.itemcget(i, "foreground"):
				self.objList.itemconfig(i,foreground=color)

			if flag & BIT_SELECT:
				self.objList.selection_set(i)
				last = i

		if last is not None:
			self.objList.activate(last)
			self.objList.see(last)

		self.selectShow(False)

	# ----------------------------------------------------------------------
	# Object list popup menu
	# ----------------------------------------------------------------------
	def objListPopup(self, event):
		"""Create pop-up menu with default actions"""
		self.page.focus_set()
		menu=Menu(self.page, tearoff=0)
		# --- Std menu ---
		menu.add_command(label="Cut", underline=2,
				compound=LEFT, image=tkFlair.icons["cut"],
				command=self.cut,
				accelerator="Ctrl-X")
		menu.add_command(label="Copy", underline=0,
				compound=LEFT, image=tkFlair.icons["copy"],
				command=self.copy,
				accelerator="Ctrl-C")
		menu.add_command(label="Paste", underline=0,
				compound=LEFT, image=tkFlair.icons["paste"],
				command=self.paste,
				accelerator="Ctrl-V")
		menu.add_command(label="Clone", underline=2,
				compound=LEFT, image=tkFlair.icons["clone"],
				command=self.clone,
				accelerator="Ctrl-D")

		# --- Add menu ---
		menu.add_separator()
		self.viewers._addMenus(menu)
		menu.add_command(label="Delete", underline=1,
				command=self.delete,
				accelerator="Del",
				compound=LEFT, image=tkFlair.icons["x"])

		# --- Visibility and Selection ---
		menu.add_separator()
		self.viewers._visibilityMenu(menu)
		self.viewers._wireframeMenu(menu)
		self.viewers._freezeMenu(menu)

		# -----------------------------------
		menu.add_separator()
		menu.add_command(label="All", underline=0,
				command=self.objList.selectAll,
				compound=LEFT, image=tkFlair.icons["all"])
		menu.add_command(label="Clear", underline=0,
				command=self.selectionClear,
				compound=LEFT, image=tkFlair.icons["none"])
		menu.add_command(label="Invert", underline=0,
				command=self.objList.selectInvert,
				compound=LEFT, image=tkFlair.icons["toggle"])

		menu.tk_popup(event.x_root, event.y_root)
		return "break"

	# ----------------------------------------------------------------------
	# Clear the selection except when the Control is pressed
	# ----------------------------------------------------------------------
	def objListClick(self, event=None):
		if not event.state&GeometryViewer.CONTROL_MASK:
			# unselect everything
			self.selectionClear(False)

	# ----------------------------------------------------------------------
	# Button-1 released
	# ----------------------------------------------------------------------
	def objListRelease(self, event=None):
		# Find activated item
		active = int(event.widget.nearest(event.y))
		if 0 <= active < len(self.objListCard):
			self.viewers.updateEdit(self.objListCard[active])
			self.viewers.expose()
#			if card.type() == Input.Card.BODY:
#				self.viewers.geometry().edit("body",card[ID])
#				self.viewers.expose()

	# ----------------------------------------------------------------------
	# Called only if the selection on object listbox is changed
	# ----------------------------------------------------------------------
	def objListSelect(self, event=None):
		"""Handle the <<ListboxSelect>> message on the object list"""
		if self.objList.size()==0: return
		# (Un)select objects
		for i in range(self.objList.size()):
			card = self.objListCard[i]
			if self.objList.selection_includes(i):
				card[SELECT] = card.get(SELECT,0) | BIT_SELECT
			else:
				card[SELECT] = card.get(SELECT,0) & BIT_UNSELECT
		if not self.showMaterial.get():
			for c in Input._defaultMaterials:
				if c.get(SELECT,0) & BIT_SELFREEZE:
					c[SELECT] = c.get(SELECT,0) & BIT_UNSELECT
			for c in Input._icruMaterials:
				if c.get(SELECT,0) & BIT_SELFREEZE:
					c[SELECT] = c.get(SELECT,0) & BIT_UNSELECT

		# Find active entry and show cursor on it
		active = self.objList.index(ACTIVE)
		card = self.objListCard[active]
		self.viewers.updateEdit(card)
		self.selectShow(True)

	# ----------------------------------------------------------------------
	# Set selected items in objList box
	# ----------------------------------------------------------------------
	def objListSelectSet(self):
		if self.showVisible.get():
			self.objListFill()
			return
		self.objList.selection_clear(0,END)
		last = None
		for i in range(self.objList.size()):
			if self.objListCard[i].get(SELECT, 0) & BIT_SELECT:
				self.objList.selection_set(i)
				last = i

		if last is not None:
			self.objList.activate(last)
			self.objList.see(last)

	# ----------------------------------------------------------------------
	# List has been sorted correct the objListCard and set Colors
	# ----------------------------------------------------------------------
	def objListSort(self, event=None):
		# reorder the objListCard
		oc = self.objListCard[:]
		self.objListCard = []
		for i in self.objList._sortOrder:
			self.objListCard.append(oc[i])
		self.objListColors()

	# ----------------------------------------------------------------------
	# Double click on object list to change name
	# ----------------------------------------------------------------------
	def objListName(self, event=None):
		active = self.objList.lists[1].index(ACTIVE)
		if active is not None and 0 <= int(active) < len(self.objListCard):
			card = self.objListCard[int(active)]
			edit = tkExtra.InPlaceEdit(self.objList.lists[1], exportselection=True)
			if edit.value is not None:
				self.changeName(card.name(), edit.value, card)

	# ----------------------------------------------------------------------
	# Restore/Set objList colors
	# Should be called after sort and after any change to
	#	selection concerning Visibility and Lock
	# FIXME should highlight the errors!
	# ----------------------------------------------------------------------
	def objListColors(self):
		for i in range(self.objList.size()):
			card = self.objListCard[i]
			flag = card.get(SELECT,0)
			if card.invalid:
				color = ERROR_COLOR
			elif flag & BIT_LOCK:
				color = LOCK_COLOR
			elif flag & BIT_FREEZE:
				color = FREEZE_COLOR
			elif flag & BIT_VISIBLE:
				color = VISIBLE_COLOR
			elif flag & BIT_WIREFRAME:
				color = WIREFRAME_COLOR
			else:
				color = ""
			if color != self.objList.itemcget(i, "foreground"):
				self.objList.itemconfig(i,foreground=color)

	# ----------------------------------------------------------------------
	# Delete action on object list
	# delete only highlighted items in objList
	# ----------------------------------------------------------------------
	def objListDelete(self, event=None):
		undoinfo = []
		bodiesRemoved = []
		regionsRemoved = []
		for item in reversed(list(map(int,self.objList.curselection()))):
			card = self.objListCard[item]
			if not undoinfo:
				# if we don't renumber cards everything gets screwed up
				undoinfo.append(self.flair.refreshUndo())
				#undoinfo.append(self.flair.renumberUndo())

			card[SELECT] = card.get(SELECT,0) & BIT_UNSELECT
			if card.type() == Input.Card.BODY:
				bodiesRemoved.append(card.name())
			elif card.type() == Input.Card.REGION:
				regionsRemoved.append(card.name())
			undoinfo.append(self.flair.delCardUndo(card.pos(),True))

		# Correct references to bodies
		for body in bodiesRemoved:
			undoinfo.extend(self.flair.removeBodyRefUndo(body))

		# Correct references to regions
		for region in regionsRemoved:
			undoinfo.extend(self.flair.removeRegionRefUndo(region))

		if undoinfo:
			#undoinfo.append(self.flair.renumberUndo())
			undoinfo.append(self.flair.refreshUndo())
			self.addUndo(undo.createListUndo(undoinfo, "Delete card(s)"))
		self.refreshUndoButton()

	# ----------------------------------------------------------------------
	# Perform the selection of objects
	# action < 0	un-select
	# action = 0	select
	# action > 0	toggle
	# ----------------------------------------------------------------------
	def select(self, card, clear=False, action=0):
		#say("GeometryEditor.select card=",card,"clear=",clear,"action=",action)
		if isinstance(card,tuple):
			obj, idx = card
			if card[0] == "Z":
				# Special treatment select zones
				if clear:
					self.propList.selection_clear(0,END)
				zz  = 0
				for z in card[1]:
					zz = z+3	# XXX WARNING should be in sync
							# with RegionProperties
					if action > 0:	# toggle
						if self.propList.selection_includes(zz):
							self.propList.selection_clear(zz)
						else:
							self.propList.selection_set(zz)
					else:
						self.propList.selection_set(zz)
				if zz: self.propList.see(zz)
				self.propListSelect(clear=clear)
				return None

		if clear:
			for c in self.project.input.allcards():
				c[SELECT] = c.get(SELECT,0) & BIT_UNSELECT
			for c in Input._defaultMaterials:
				if c.get(SELECT,0) & BIT_SELFREEZE:
					c[SELECT] = c.get(SELECT,0) & BIT_UNSELECT
			for c in Input._icruMaterials:
				if c.get(SELECT,0) & BIT_SELFREEZE:
					c[SELECT] = c.get(SELECT,0) & BIT_UNSELECT

		if isinstance(card, Input.Card):
			if action > 0:		# toggle
				card[SELECT] = card.get(SELECT,0) ^ BIT_SELECT
			elif action == 0:	# select
				card[SELECT] = card.get(SELECT,0) | BIT_SELECT
			else:			# un-select
				card[SELECT] = card.get(SELECT,0) & BIT_UNSELECT

		elif isinstance(card, list):
			for c in card:
				if action > 0:
					if c.get(SELECT,0) & BIT_SELECT:
						action = -1	# un-select
					else:
						action =  0	# select

				if action == 0:		# select
					c[SELECT] = c.get(SELECT,0) | BIT_SELECT
				else:			# un-select
					c[SELECT] = c.get(SELECT,0) & BIT_UNSELECT
			card = c # remember last card

		else:
			card = None
			if isinstance(idx, list):
				for i in idx:
					card = self.findCard(obj,i)
					if card is None: continue
					if action > 0:	# toggle
						if card.get(SELECT,0) & BIT_SELECT:
							action = -1	# un-select
						else:
							action =  0	# select

					if action == 0:	# select
						card[SELECT] = card.get(SELECT,0) | BIT_SELECT
					else:		# un-select
						card[SELECT] = card.get(SELECT,0) & BIT_UNSELECT
			else:
				card = self.findCard(obj,idx)
				if card is not None:
					if action > 0:		# toggle
						card[SELECT] = card.get(SELECT,0) ^ BIT_SELECT
					elif action == 0:	# select
						card[SELECT] = card.get(SELECT,0) | BIT_SELECT
					else:
						card[SELECT] = card.get(SELECT,0) & BIT_UNSELECT

		# Populate listbox
		self.objListSelectSet()

		# Update last edit
		self.viewers.updateEdit(card)

		# show selected items
		self.selectShow()

	# ----------------------------------------------------------------------
	# Show on all viewports the selected and visible objects
	# @return number of selected objects
	# ----------------------------------------------------------------------
	def selectShow(self, redraw=True, updateprop=True):
		"""Show on viewports the selected and visible objects"""
		self.viewers.body(  "*", "show",   0)
		self.viewers.region("*", "show",   0)
		self.viewers.object("*", "show",   0)
		self.viewers.object("*", "select", 0)	# select 0 clears node also
		self.viewers.zone("clear")

		selcards = []

		# FIXME check if a selection flag has changed from bodies/regions/zones
		# then only redraw!
#		needredraw = False
		# Scan on materials cards if any
		#if self.showMaterial.get():
		for card in Input._defaultMaterials:
			if card.get(SELECT,0) & BIT_SELFREEZE:
				selcards.append(card)
		for card in Input._icruMaterials:
			if card.get(SELECT,0) & BIT_SELFREEZE:
				selcards.append(card)

		# Scan over all cards (input+voxel)
		for card in self.project.input.allcards():
			if card.ignore(): continue
			flag = card.get(SELECT,0)
#			needredraw = True
			if flag == 0: continue
			if card.type() == Input.Card.BODY:
				self.viewers.body(card[ID], "show", flag)
				if flag & BIT_SELECT:
					self.viewers.body(card[ID], "color",
						GeometryViewer._SELECT_COLOR)
#					self.viewers.body(card[ID], "linewidth", 2)
				elif flag & BIT_FREEZE:
					self.viewers.body(card[ID], "color",
						GeometryViewer._FREEZE_COLOR)
#					self.viewers.body(card[ID], "linewidth", 0)
				elif flag & BIT_LOCK:
					self.viewers.body(card[ID], "color",
						GeometryViewer._LOCK_COLOR)
#					self.viewers.body(card[ID], "linewidth", 0)
				elif flag & BIT_VISIBLE:
					self.viewers.body(card[ID], "color",
						GeometryViewer._VISIBLE_COLOR)
#					self.viewers.body(card[ID], "linewidth", 0)

			elif card.type()==Input.Card.REGION:
#				needredraw = True
				if card.name() == "&": continue
				try:
					self.viewers.region(card[ID], "show", flag)
				except IndexError:
					pass

			elif card.type()==Input.Card.OBJECT or card.tag in _CARD_OBJECTS:
				try:
					self.viewers.object(card[ID], "show", 1)
					if flag & BIT_UNVISIBLE:
						self.viewers.object(card[ID], "select",
								GeometryViewer._SELECT_COLOR)
				except IndexError:
					pass

			# Append to the list of selected cards
			if flag & BIT_SELFREEZE: selcards.append(card)
			#if flag & BIT_SELECT: selcards.append(card)

		# Redraw viewports if required
		if redraw: self.viewers.draw()

		# Do not update properties
		if not updateprop: return

		# Update property list only if edit button is not pressed
		if not selcards:
			self.showNothing()
		else:
			self.showProperties(selcards)

		self.viewers.updateSelection()

	# ----------------------------------------------------------------------
	# @return true if any item is selected
	# ----------------------------------------------------------------------
	def isAnythingSelected(self, cardtype=None):
		for card in self.project.input.allcards():
			if card.ignore(): continue
			if cardtype is not None and card.type() != cardtype:
				continue
			if card.get(SELECT,0) & BIT_SELECT:
				return True
		return False

	# ----------------------------------------------------------------------
	# Selection clear (ignore editing object)
	# ----------------------------------------------------------------------
	def selectionClear(self, show=True, cardtype=None):
		"""Clear the selection list"""
		for card in self.project.input.allcards():
			if card.ignore(): continue
			if cardtype is not None and card.type() != cardtype: continue
			card[SELECT] = card.get(SELECT,0) & BIT_UNSELECT
		if show:
			self.objListSelectSet()
			self.selectShow()

	# ----------------------------------------------------------------------
	# @return true if object is selected
	# ----------------------------------------------------------------------
	def isSelected(self, card):
		"""Return true if card is selected"""
		if isinstance(card,tuple):
			card = self.findCard(*card)
		if card:
			return card.get(SELECT,0) & BIT_SELECT
		else:
			return 0

	# ----------------------------------------------------------------------
	def getSelection(self, mask=BIT_SELECT, objtype=None):
		"""Return the list of selected objects with mask"""
		objects = []
		for card in self.project.input.allcards():
			if card.ignore(): continue
			if card.get(SELECT,0) & mask:
				if objtype is None or card.type() == objtype:
					objects.append(card)
		return objects

	# ----------------------------------------------------------------------
	def selectAllBodies(self):
		"""Select all bodies from all selected regions"""
		found = False
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.type() != Input.Card.REGION: continue
			if card.get(SELECT,0) & BIT_SELECT == 0: continue
			found = True
			for name in csg.tokenize(card.extra()):
				if name in ("(",")","-","+","|"): continue
				body = self.findCard("B", name)
				body[SELECT] = body.get(SELECT,0) | BIT_SELECT
		self.objListSelectSet()
		self.selectShow()
		return found

	# ----------------------------------------------------------------------
	def getVisible(self):
		"""Return the list of visible objects"""
		return self.getSelection(BIT_VISIBLE)

	# ----------------------------------------------------------------------
	def setBit(self, bit, types, msg, tags=()):
		"""set bit to all selected objects"""
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.type() not in types and card.tag not in tags: continue
			value = card.get(SELECT,0)
			if value & BIT_SELECT:
				card[SELECT] = value | bit
		self.setStatus(msg)

	# ----------------------------------------------------------------------
	def unsetBit(self, bit, types, msg, tags=()):
		"""unset the selection bit to all selected objects"""
		unbit = ~bit
		mask  = bit | BIT_SELECT
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.type() not in types and card.tag not in tags: continue
			value = card.get(SELECT,0)
			if value & mask == mask:
				card[SELECT] = value & unbit
		self.setStatus("Un-"+msg)

	# ----------------------------------------------------------------------
	def toggleBit(self, bit, types, msg, tags=()):
		"""Change bit on selected items"""
		for card in self.project.input.allcards():
			if card.ignore(): continue
			if card.type() not in types and card.tag not in tags: continue
			value = card.get(SELECT,0)
			if value & BIT_SELECT:
				if value & bit:
					self.unsetBit(bit, types, msg, tags)
				else:
					self.setBit(bit, types, msg, tags)
				return

	# ----------------------------------------------------------------------
	def setBitAll(self, bit, types, msg, tags=()):
		"""unset the selection bit to all objects"""
		unbit = ~bit
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.type() not in types and card.tag not in tags: continue
			card[SELECT] = card.get(SELECT,0) | bit
		self.setStatus(msg)

	# ----------------------------------------------------------------------
	def unsetBitAll(self, bit, types, msg, tags=()):
		"""unset the selection bit to all objects"""
		unbit = ~bit
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.type() not in types and card.tag not in tags: continue
			value = card.get(SELECT,0)
			if value & bit == bit:
				card[SELECT] = value & unbit
		self.setStatus("Un-"+msg)

	# ----------------------------------------------------------------------
	def toggleBitAll(self, bit, types, msg, tags=()):
		"""Change bit on all items"""
		for card in self.project.input.allcards():
			if card.ignore(): continue
			if card.type() not in types and card.tag not in tags: continue
			value = card.get(SELECT,0)
			if value & bit:
				self.unsetBitAll(bit, types, msg, tags)
			else:
				self.setBitAll(bit, types, msg, tags)
			return

	# ----------------------------------------------------------------------
	def toggleLock(self):
		self.lockEditing.set(not self.lockEditing.get())
		if self.lockEditing.get():
			self.setStatus("Graphical editing is LOCKED")
		else:
			self.setStatus("Graphical editing is UN-LOCKED")

	# ----------------------------------------------------------------------
	def _freezeEnd(self):
		self.viewers.check4Frozen()
		self.objListColors()
		self.updateProperties()
		self.selectShow()

	def toggleFreeze(self):
		self.toggleBit(BIT_FREEZE, (Input.Card.BODY, Input.Card.REGION), "Freeze selected objects")
		self._freezeEnd()

	def freeze(self):
		self.setBit(BIT_FREEZE, (Input.Card.BODY, Input.Card.REGION), "Freeze selected objects")
		self._freezeEnd()

	def unfreeze(self):
		self.unsetBit(BIT_FREEZE, (Input.Card.BODY, Input.Card.REGION), "Freeze selected objects")
		self._freezeEnd()

	def unfreezeAll(self):
		self.unsetBitAll(BIT_FREEZE, (Input.Card.BODY, Input.Card.REGION), "Freeze all objects")
		self._freezeEnd()

	# ----------------------------------------------------------------------
	def _3dEnd(self):
		self.objListColors()
		self.selectShow(True, False)
		self.updateProperties()

	def toggle3D(self):
		self.toggleBit(BIT_WIREFRAME, (Input.Card.BODY,), "Set 3D wireframe")
		self._3dEnd()

	def set3D(self):
		self.setBit(BIT_WIREFRAME, (Input.Card.BODY,), "Set 3D wireframe")
		self._3dEnd()

	def unset3D(self):
		self.unsetBit(BIT_WIREFRAME, (Input.Card.BODY,), "Set 3D wireframe")
		self._3dEnd()

	def unset3DAll(self):
		self.unsetBitAll(BIT_WIREFRAME, (Input.Card.BODY,), "Set 3D wireframe")
		self._3dEnd()

	# ----------------------------------------------------------------------
	def _visibilityEnd(self):
		self.objListColors()
		self.updateProperties()
		self.selectShow()

	def toggleVisibility(self):
		self.toggleBit(BIT_VISIBLE, (Input.Card.BODY,Input.Card.OBJECT),
				"Set visibility", _CARD_OBJECTS)
		self._visibilityEnd()

	def setVisibility(self):
		self.setBit(BIT_VISIBLE, (Input.Card.BODY,Input.Card.OBJECT),
				"Set visibility", _CARD_OBJECTS)
		self._visibilityEnd()

	def unsetVisibility(self):
		self.unsetBit(BIT_VISIBLE, (Input.Card.BODY,Input.Card.OBJECT),
				"Set visibility", _CARD_OBJECTS)
		self._visibilityEnd()

	def toggleVisibilityAll(self):
		self.toggleBitAll(BIT_VISIBLE, (Input.Card.BODY,Input.Card.OBJECT),
				"Set visibility", _CARD_OBJECTS)
		self._visibilityEnd()

	def setAllVisibility(self):
		self.setBitAll(BIT_VISIBLE, (Input.Card.BODY,Input.Card.OBJECT),
				"Set visibility", _CARD_OBJECTS)
		self._visibilityEnd()

	def unsetAllVisibility(self):
		self.unsetBitAll(BIT_VISIBLE, (Input.Card.BODY,Input.Card.OBJECT),
				"Set visibility", _CARD_OBJECTS)
		self._visibilityEnd()

	# ----------------------------------------------------------------------
	# Toggle the alpha channel for all selected regions
	# ----------------------------------------------------------------------
	def toggleAlpha(self, a=None):
		undoinfo = []
		for card in self.project.input["REGION"]:
			if card.ignore(): continue
			value = card.get(SELECT,0)
			if value & BIT_SELECT:
				if a is None:
					alpha = 100 - card.get(ALPHA,0)
				else:
					alpha = a
				undoinfo.append(self.flair.setCardUndo(card, ALPHA, alpha, 0))
		if undoinfo:
			self.flair.addUndo(undo.createListUndo(undoinfo,
					undoinfo[-1][0]))
			self.redraw()
			self.updateProperties()

	# ----------------------------------------------------------------------
	# Return all movable selected objects id
	# Movable are: bodies and objects
	# ----------------------------------------------------------------------
	def getMoveableSelection(self, mask=BIT_SELECT):
		"""Return selection list of moveable objects (bodies and objects)"""
		objects = []
		bodies  = {}
		bodiesSelected = False	# Avoid moving the regions if bodies&regions are selected

		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.get(SELECT,0) & mask == 0: continue

			if card.type()==Input.Card.BODY:
				bodies[card.name()] = card
				bodiesSelected = True

			elif card.type()==Input.Card.OBJECT or card.tag=="ROT-DEFI":
				objects.append(card)

			elif card.type()==Input.Card.REGION:
				# Return bodies of region if only
				# regions are selected! No bodies
				if bodiesSelected: continue
				for name in csg.tokenize(card.extra()):
					if name in ("(",")","-","+","|"): continue
					if name in bodies: continue
					body = self.findCard("B", name)
					bodies[name] = body

		return (objects, list(bodies.values()))

	# ----------------------------------------------------------------------
	# Find Input Card from obj, idx
	# ----------------------------------------------------------------------
	def findCard(self, obj, idx):
		if obj == "B":
			if isinstance(idx,int):
				# Search by index
				for card in self.project.input.cardlist:
					if not card.isGeo(): continue
					if card.ignore(): continue
					if len(card.tag) != 3 and card.tag!="VOXELS": continue
					if idx == card[ID]:
						return card
			else:
				# Search by name
				for card in self.project.input.cardlist:
					if not card.isGeo(): continue
					if card.ignore(): continue
					if len(card.tag) != 3 and card.tag!="VOXELS": continue
					if card.tag=="END": continue
					if card.tag=="VOXELS":
						if idx=="VOXEL":
							return card
					else:
						if idx==card.name():
							return card

		elif obj == "R":
			if isinstance(idx,int):
				if idx >= 0:
					# Scan in normal input
					for card in self.project.input["REGION"]:
						if card.ignore(): continue
						if idx == card[ID]:
							return card
				else:
					# Scan in voxels by name!!!
					idx = -idx
					if idx < 1000:
						name = "VOXEL%03d"%(idx)
					elif idx < 10000:
						name = "VOXE%04d"%(idx)
					else:
						name = "VOX%05d"%(idx)
					for voxel in self.project.input["VOXELS"]:
						if voxel.ignore() or voxel["@voxel"] is None: continue
						for card in voxel["@voxel"].input["REGION"]:
							if card.name() == name:
							#if idx == card[ID]:
								return card
			else:
				for card in self.project.input["REGION"]:
					if card.ignore(): continue
					if idx == card.name(): return card

		elif obj=="O":
			if isinstance(idx,int):
				for tag in Input.OBJECT_TAGS:
					for card in self.project.input[tag]:
						if card.ignore(): continue
						if idx == card[ID]:
							return card
			else:
				for tag in Input.OBJECT_TAGS:
					for card in self.project.input[tag]:
						if card.ignore(): continue
						if idx == card.name():
							return card

		elif obj=="T":
			if isinstance(idx,int):
				for card in self.project.input["ROT-DEFI"]:
					if card.ignore(): continue
					if idx == card[ID]:
						return card
			else:
				for card in self.project.input["ROT-DEFI"]:
					if card.ignore(): continue
					if idx == card.name(): return card
		return None

	# ----------------------------------------------------------------------
	def zoomOnSelected(self, event=None):
		"""Zoom on selected objects"""
		objects = self.getSelection()
		if objects is None: return

		for frame in self.viewers.frames:
			frame.zoomOnSelected(objects)

	# ----------------------------------------------------------------------
	def showNothing(self):
		self.propList.delete(0,END)
		self.attrList.delete(0,END)
		self.showContextGroup()
		self.editObject = None

	# ----------------------------------------------------------------------
	# Show properties of selected object
	# ----------------------------------------------------------------------
	def showProperties(self, cards):
		"""Show properties of the selected object"""
		select = [0,0,0,0]
		freeze = [0,0,0,0]
		tags  = {}	# FIXME replace with set
		for c in cards:
			flag = c.get(SELECT,0)
			if flag & BIT_SELECT:
				select[c.type()] += 1
			elif flag & BIT_FREEZE:
				freeze[c.type()] += 1
			tags[c.tag] = True
		tags  = list(tags.keys())

		lbody = select[Input.Card.BODY]
		lreg  = select[Input.Card.REGION]
		lobj  = select[Input.Card.OBJECT]
		lgen  = select[Input.Card.GENERIC]
		#print ">>> B:%d  R:%d  O:%d  G:%d"%(lbody,lreg,lobj,lgen)

		# XXX Check before if editObject contain the same cards
		if lreg==1 or \
		   (lreg==0 and freeze[Input.Card.REGION]==1 and lobj==0 and lgen==0):
			self.showContextGroup("region")
			# Find first region
			if lreg==1: # Which is selected
				for c in cards:
					if c.tag == "REGION" and c.get(SELECT,0) & BIT_SELECT:
						region = c
						break
			else: # or locked
				for c in cards:
					if c.tag == "REGION" and c.get(SELECT,0) & BIT_FREEZE:
						region = c
						break
			if self.editObject is not None and \
			   len(self.editObject.cards) == 1 and \
			   self.editObject.cards[0] is region:
				# Do nothing, most probably zone selection
				return
			self.editObject = RegionProperties(self, [region])

		elif lbody!=0 and lobj==0 and lreg==0 and lgen==0:
			selcards = [x for x in cards if x[SELECT]&BIT_SELECT != 0]
			if selcards:
				self.editObject = BodyProperties(self, selcards)
				self.showContextGroup("body")
			else:
				self.showNothing()
				return

		elif lbody==0 and lobj==0 and lreg!=0 and lgen==0:
			self.editObject = RegionProperties(self, cards)
			self.showContextGroup("region")

		elif lbody==0 and lobj!=0 and lreg==0 and lgen==0:
			self.editObject = ObjectProperties(self, cards)
			if self.editObject.objtype == "!mesh":
				self.showContextGroup("mesh")
			else:
				self.showContextGroup("object")

		#elif lbody==0 and lobj==0 and lreg==0 and lgen!=0 and len(tags)==1:
		elif lgen!=0:
			if tags[0] == "MATERIAL":
				self.editObject = MaterialProperties(self, cards)
				self.showContextGroup("material")
			elif tags[0] == "ROT-DEFI":
				self.editObject = RotdefiProperties(self, cards)
				self.showContextGroup("rotdefi")
			elif tags[0] in ("BEAM","SPOTBEAM"):
				self.editObject = BeamProperties(self, cards)
				self.showContextGroup("beam")
			else:
				self.showNothing()
				return

		else:
			self.showNothing()
			return

		# FIXME highlight sub if any
		self.updateProperties()
		self.propListSelect()

	# ----------------------------------------------------------------------
	# Update the properties and attribute lists of the current object
	# ----------------------------------------------------------------------
	def updateProperties(self):
		#print "\nupdateProperties"
		#import traceback; traceback.print_stack()
		if self.editObject is None: return
		if not self._loadingProject:
			self.flair.activateIfNeededInputFrame(self.editObject.cards[0])
		self.editObject.populateList()

	# ----------------------------------------------------------------------
	def editType(self):
		"""Return type of editing objects if any else None"""
		if self.editObject is None or self.editObject.many:
			return None
		return self.editObject.type()

	# ----------------------------------------------------------------------
	# Property List popup menu
	# ----------------------------------------------------------------------
	def propListPopup(self, event):
		"""Create popup menu with default actions"""
		if self.editObject is None: return
		self.page.focus_set()

		menu=Menu(self.page, tearoff=0)

		menu.add_command(label="Cut", underline=2,
				command=self.propListCut,
				compound=LEFT, image=tkFlair.icons["cut"],
				accelerator="Ctrl-X")
		menu.add_command(label="Copy", underline=0,
				command=self.propListCopy,
				compound=LEFT, image=tkFlair.icons["copy"],
				accelerator="Ctrl-C")
		menu.add_command(label="Paste", underline=0,
				command=self.propListPaste,
				compound=LEFT, image=tkFlair.icons["paste"],
				accelerator="Ctrl-V")
		menu.add_command(label="Clone", underline=0,
				state=DISABLED,
				compound=LEFT, image=tkFlair.icons["clone"])
		if self.editObject.type() == Input.Card.REGION:
			menu.add_separator()
			menu.add_command(label="Edit Material", underline=0,
					command=self.gotoMaterial,
					compound=LEFT, image=tkFlair.icons["MATERIAL"],
					accelerator="Ctrl-M")
		menu.add_separator()
		menu.add_command(label="Delete", underline=0,
				command=self.propListDelete,
				compound=LEFT, image=tkFlair.icons["x"])
		menu.add_separator()

		menu.add_command(label="All", underline=0,
				command=self.propList.selectAll,
				compound=LEFT, image=tkFlair.icons["all"])
		menu.add_command(label="Clear", underline=0,
				command=self.propList.selectClear,
				compound=LEFT, image=tkFlair.icons["none"])
		menu.add_command(label="Invert", underline=0,
				command=self.propList.selectInvert,
				compound=LEFT, image=tkFlair.icons["toggle"])

		menu.tk_popup(event.x_root, event.y_root)
		return "break"

	# ----------------------------------------------------------------------
	# Copy properties
	# ----------------------------------------------------------------------
	def propListCopy(self, event=None):
		if self.editObject is None: return
		lst = self.propList.lists[1]
		if self.editObject.type() == Input.Card.REGION:
			txt = " ".join([str(lst.get(i)) for i in self.propList.curselection()],"\n|")
		else:
			txt = " ".join([str(lst.get(i)) for i in self.propList.curselection()])
		self.page.clipboard_clear()
		self.page.clipboard_append(txt)
		return "break"

	# ----------------------------------------------------------------------
	# Cut properties
	# ----------------------------------------------------------------------
	def propListCut(self, event=None):
		if self.editObject is None: return
		self.propListCopy()
		self.propListDelete()
		return "break"

	# ----------------------------------------------------------------------
	# Paste properties
	# ----------------------------------------------------------------------
	def propListPaste(self, event=None):
		if self.editObject is None: return
		if self.editObject.type() == Input.Card.REGION:
			clipboard = self.flair.clipboard()
			if clipboard.startswith(tkFlair._FLAIRCLIP):
				self.flair.notify("Paste error",
					"Clipboard contains cards and not text.",
					tkFlair.NOTIFY_WARNING)
				return
			self.addZone(clipboard)
		return "break"

	# ----------------------------------------------------------------------
	def addZone(self, zone):
		"""Add or update zone information from mouse selection ACTION_ZONE"""
		if self.editObject.type() != Input.Card.REGION: return
		if self.editObject.many: return

		# Find if any zone is selected
		sel = self.propList.curselection()
		if len(sel)==1:
			name = self.propList.get(sel[0])[0]
			if name.startswith("zone"):
				i = int(name[4:])-1
				zones = self.editObject.cards[0]["@zone"]
				zones[i] = zone
				self.flair.addUndo(
					self.flair.setZonesUndo(self.editObject.cards[0], zones))
			else:
				self.flair.addUndo(
					self.flair.addZoneUndo(self.editObject.cards[0], zone))
		elif len(sel)>1:
			self.flair.notify("Too many selected zones",
				"Too many zones are selected. Please select " \
				"only one zone to update, or none to add a new one",
				tkFlair.NOTIFY_WARNING)
			return
		else:
			self.flair.addUndo(
				self.flair.addZoneUndo(self.editObject.cards[0], zone))

		# FIXME Show properties after re-parsing the zones...
		self.updateProperties()
#		self.showProperties(self.editObject.cards)
#		self.propListSelect()

	# ----------------------------------------------------------------------
	# Return # of zones selected:
	# ----------------------------------------------------------------------
	def propListZoneSelected(self):
		# Find if any zone is selected
		count = 0
		for sel in self.propList.curselection():
			if self.propList.get(sel)[0].startswith("zone"):
				count += 1
		return count

	# ----------------------------------------------------------------------
	# Return selected zones
	# ----------------------------------------------------------------------
	def propListSelectedZones(self):
		if self.editObject.type() == Input.Card.REGION:
			zones = []
			for sel in self.propList.curselection():
				name,zone = self.propList.get(sel)
				if not name.startswith("zone"): continue
				zones.append(int(name[4:])-1)
			return zones
		return None

	# ----------------------------------------------------------------------
	def propListClear(self):
		self.propList.selection_clear(0,END)
		self.propListSelect()

	# ----------------------------------------------------------------------
	# Selection changed in propList
	# For Regions:
	#	- Select all bodies appearing in selected zones
	# ----------------------------------------------------------------------
	def propListSelect(self, event=None, clear=True):
		if self.editObject is None: return

		if self.editObject.type() == Input.Card.REGION:
			# Select all bodies from zone
			bodies = set()
			selzones = []
			for sel in self.propList.curselection():
				name,zone = self.propList.get(sel)
				if not name.startswith("zone"): continue
				# find bodies
				try:
					for body in zone.translate(OPER_TRANS).split():
						bodies.add(body)
				except TypeError:	# Maybe some unicode will not work
					pass

				# Select zone in viewers
				try:
					exp = csg.splitZones(csg.tokenize(zone))
				except csg.CSGException:
					messagebox.showerror("Zone error",
							sys.exc_info()[1],
							parent=self.page)
				else:
					try:
						selzones.append(GeometryViewer.regionExp(exp)[0])
					except:
						selzones.append(zone)

			if bodies:
				# Clear any body selection
				if clear:
					for card in self.project.input.cardlist:
						if card.ignore(): continue
						if card.type() == Input.Card.BODY:
							card[SELECT] = card.get(SELECT,0) & BIT_UNSELECT

				# Select bodies, but don't use the self.select() to avoid
				# the filling of the propList and de-selection of zone
				for body in bodies:
					card = self.findCard("B",body)
					if card is not None:
						card[SELECT] = card.get(SELECT,0) | BIT_SELECT
				self.objListFill()

				if selzones:
					self.viewers.stopThreads()
					self.viewers.zone("select", selzones)
					self.viewers.zone("show",
						self.viewers.region(
							self.editObject.cards[0][ID],"show"))
				self.viewers.draw()

			else:
				self.viewers.zone("clear")
				self.viewers.draw()

			self.viewers.updateEdit()
			del bodies

		else:
			# FIXME
			pass
		self.viewers.updateSelection()

	# ----------------------------------------------------------------------
	def propListRelease(self, event):
		"""Button Release on the properties listbox"""
		act = event.widget.nearest(event.y)
		self.propList.activate(act)
		# Find list box
		if self.propList.lists.index(event.widget) == 0: return
		self.propListEdit()
		return "break"

	# ----------------------------------------------------------------------
	def editName(self):
		if self.editObject is None: return
		self._tabSet.changePage("Geometry")
		self.propTab.changePage(_PROPTABS[0])
		self.propList.selection_clear(0,END)
		self.propList.selection_set(0)
		self.propList.activate(0)
		self.propListEdit(False)

	# ----------------------------------------------------------------------
	def editMaterial(self):
		if self.editObject is None or self.editObject.type() != Input.Card.REGION: return
		self._tabSet.changePage("Geometry")
		self.propTab.changePage(_PROPTABS[0])
		self.propList.selection_clear(0,END)
		self.propList.selection_set(2)
		self.propList.activate(2)
		self.propListEdit(False)

	# ----------------------------------------------------------------------
	def editProperty(self):
		if self.editObject is None: return
		self._tabSet.changePage("Geometry")
		if self.propTab.getActivePage() == _PROPTABS[0]:
			self.propList.selection_clear(0,END)
			self.propList.selection_set(ACTIVE)
			self.propListEdit(False)
		else:
			self.attrList.selection_clear(0,END)
			self.attrList.selection_set(ACTIVE)
			self.attrListEdit(False)

	# ----------------------------------------------------------------------
	def gotoMaterial(self, event=None):
		if self.editObject and self.editObject.type() == Input.Card.REGION:
			self._tabSet.changePage("Geometry")
			materials = [card["@material"] for card in self.editObject.cards]
			self.select(materials, True)

	# ----------------------------------------------------------------------
	# ----------------------------------------------------------------------
	# Show body transformation ROT-DEFI card
	# ----------------------------------------------------------------------
	def gotoRotdefi(self):
		if not self.editObject: return
		if self.editObject.type() == Input.Card.BODY:
			key = "@transform"
		elif self.editObject.type() == Input.Card.REGION:
			key = "@rotdefi"
		lst = [card[key] for card in self.editObject.cards]
		if lst[0] == None: return
		# Find cards to display
		rotdefi = []
		for r in lst:
			if r and r[0]=="-": r=r[1:]
			rot = self.findCard("T",r)
			if rot: rotdefi.append(rot)

		self._tabSet.changePage("Geometry")
		self.select(rotdefi, True)

	# -----------------------------------------------------------------------
	# Delete action on property list
	# Regions: Delete zones
	# ----------------------------------------------------------------------
	def propListDelete(self, event=None):
		if self.editObject is None: return
		if self.editObject.type() == Input.Card.REGION:
			card = self.editObject.cards[0]
			zones = card["@zone"]
			lst = list(map(int,self.propList.curselection()))
			lst.reverse()
			for sel in lst:
				name = self.propList.lists[0].get(sel)
				if not name.startswith("zone"): continue
				i = int(name[4:])-1
				del zones[i]
			self.flair.addUndo(self.flair.setZonesUndo(card,zones))
			self.setInputModified()

			# FIXME temporarily reload all
			self.loadProject()

			# FIXME Show properties after re-parsing the zones...
			self.updateProperties()

		else:
			cursel = self.propList.curselection()
			for sel in cursel:
				name = self.propList.lists[0].get(sel)
				self.editObject.value(name,"")

			# Special check for ROT-DEFI
			if isinstance(self.editObject, RotdefiProperties):
				# if the selected card has all elements removed
				# remove also the card
				# Find index of cards
				idx = []
				for sel in cursel:
					name = self.propList.lists[0].get(sel)
					try: i = int(name.split(".")[-1])
					except: continue
					if i not in idx: idx.append(i)
				idx.reverse()
				# Check if card is empty
				undoinfo = []
				for i in idx:
					card = self.editObject._rotdefi[i-1]
					say(i,">>>",card)
					for j in range(1,7):
						if card.what(j)!="": break
					else:
						if not undoinfo:
							undoinfo = [self.flair.refreshUndo(),
								self.flair.renumberUndo()]
						undoinfo.append(self.flair.delCardUndo(card.pos(),False))
				if undoinfo:
					undoinfo.append(self.flair.renumberUndo())
					undoinfo.append(self.flair.refreshUndo())

		self.refreshUndoButton()

	# ----------------------------------------------------------------------
	# Clone action on property list
	# Regions: Clone zone
	# ----------------------------------------------------------------------
	def propListClone(self, event=None):
		if self.editObject is None: return
		if self.editObject.type() == Input.Card.REGION:
			# FIXME
			pass

	# ----------------------------------------------------------------------
	# Move Up action in properties
	# Region: Move up zone
	# ----------------------------------------------------------------------
	def propListMoveUp(self, event=None):
		if self.editObject is None: return

		if self.editObject.type() == Input.Card.REGION:
			# FIXME
			pass

	# ----------------------------------------------------------------------
	# Move Down action in properties
	# Region: Move down zone
	# ----------------------------------------------------------------------
	def propListMoveDown(self, event=None):
		if self.editObject is None: return

		if self.editObject.type() == Input.Card.REGION:
			# FIXME
			pass

	# ----------------------------------------------------------------------
	def propListEditEvent(self, event=None):
		if len(self.propList.curselection())==0:
			active = self.propList.index(ACTIVE)
			self.propList.selection_set(active)
		self.propListEdit()

	# ----------------------------------------------------------------------
	def propListEdit(self, next=True):
		if len(self.propList.curselection())!=1: return

		active = int(self.propList.curselection()[0])
		self.propList.activate(active)

		lbl, value = self.propList.get(active)
		edit = self.editObject.edit(self.propList)

		if edit is not None:
			if edit.value is None:
				return

			elif edit.value != value:
				if self.editObject.split() and edit.value.find(" ")>=0:
					# Split value and propagate change to next fields
					value = edit.value
					while True:
						m = _FIRST_WORD.match(value)
						if m is None:
							self.propList.set(active,(lbl, value))
							self.editObject.value(lbl, value)
							break
						self.propList.set(active,(lbl, m.group(1)))
						self.editObject.value(lbl, m.group(1))
						value = m.group(2)
						active += 1
						if active >= self.propList.size(): break
						self.propList.activate(active)
						lbl = self.propList.get(active)[0]
				else:
					self.editObject.value(lbl, edit.value)
					if self.editObject.type() == Input.Card.REGION:
						self.updateProperties()

			if edit.lastkey == "Up":
				active -= 2
			elif edit.lastkey not in ("Return", "KP_Enter", "Down"):
				return

		# edit next item
		if not next: return
		active += 1
		if active < self.propList.size():
			self.propList.selection_clear(0,END)
			self.propList.selection_set(active)
			self.propList.activate(active)
			self.propList.see(active)
			self.page.event_generate("<Return>")

	# ----------------------------------------------------------------------
	# Attribute List popup menu
	# ----------------------------------------------------------------------
	def attrListPopup(self, event):
		"""Create popup menu with default actions"""
		if self.editObject is None: return
		self.page.focus_set()

		menu=Menu(self.page, tearoff=0)

#		menu.add_command(label="Cut", underline=2,
#				command=self.propListCut,
#				compound=LEFT, image=tkFlair.icons["cut"],
#				accelerator="Ctrl-X")
#		menu.add_command(label="Copy", underline=0,
#				command=self.propListCopy,
#				compound=LEFT, image=tkFlair.icons["copy"],
#				accelerator="Ctrl-C")
#		menu.add_command(label="Paste", underline=0,
#				command=self.propListPaste,
#				compound=LEFT, image=tkFlair.icons["paste"],
#				accelerator="Ctrl-V")
#		menu.add_command(label="Clone", underline=0,
#				state=DISABLED,
#				compound=LEFT, image=tkFlair.icons["clone"])
#		menu.add_separator()
		menu.add_command(label="Delete", underline=0,
				command=self.attrListDelete,
				compound=LEFT, image=tkFlair.icons["x"])
		menu.add_separator()

		menu.add_command(label="All", underline=0,
				command=self.attrList.selectAll,
				compound=LEFT, image=tkFlair.icons["all"])
		menu.add_command(label="Clear", underline=0,
				command=self.attrList.selectClear,
				compound=LEFT, image=tkFlair.icons["none"])
		menu.add_command(label="Invert", underline=0,
				command=self.attrList.selectInvert,
				compound=LEFT, image=tkFlair.icons["toggle"])

		menu.tk_popup(event.x_root, event.y_root)
		return "break"

	# ----------------------------------------------------------------------
	def attrListRelease(self, event):
		"""Button Release on the attribute listbox"""
		act = event.widget.nearest(event.y)
		self.attrList.activate(act)
		# Find list box
		if self.attrList.lists.index(event.widget) == 0: return

		var, value = self.attrList.get(act)
		if value in (_ON, _OFF, _MULTI):
			if value == _OFF:
				value = _ON
			else:
				value = _OFF
			self.editObject.value(var, value==_ON)
			self.attrList.lists[1].set(act, value)
			self.attrList.lists[1].selection_set(act)
			self.viewers.check4Frozen()
			self.selectShow(True, False)
			self.objListColors()
		else:
			self.attrListEdit()
		return "break"

	# ----------------------------------------------------------------------
	def attrListDelete(self, event=None):
		if self.editObject is None: return
		for sel in self.attrList.curselection():
			name, value = self.attrList.get(sel)
			if value in (_ON, _OFF, _MULTI):
				value = False
			else:
				value = ""
			self.editObject.value(name,value)

		self.refreshUndoButton()

	# ----------------------------------------------------------------------
	def attrListEdit(self, next=True):
		if len(self.attrList.curselection())!=1: return

		active = int(self.attrList.curselection()[0])
		self.attrList.activate(active)
		lbl, value = self.attrList.get(active)
		edit = self.editObject.edit(self.attrList)

		if edit is not None:
			if edit.value is None:
				return

			elif edit.value != value:
				if self.editObject.split2() and edit.value.find(" ")>=0:
					# Split value and propagate change to next fields
					value = edit.value
					while True:
						m = _FIRST_WORD.match(value)
						if m is None:
							self.attrList.set(active,(lbl, value))
							self.editObject.value(lbl, value)
							break
						self.attrList.set(active,(lbl, m.group(1)))
						self.editObject.value(lbl, m.group(1))
						value = m.group(2)
						active += 1
						if active >= self.attrList.size(): break
						self.attrList.activate(active)
						lbl = self.attrList.get(active)[0]
				else:
					self.editObject.value(lbl, edit.value)
					if self.editObject.type() == Input.Card.REGION:
						self.updateProperties()

			if edit.lastkey == "Up":
				active -= 2
			elif edit.lastkey not in ("Return", "KP_Enter", "Down"):
				return

		# edit next item
		if not next: return
		active += 1
		if active < self.attrList.size():
			self.attrList.selection_clear(0,END)
			self.attrList.selection_set(active)
			self.attrList.activate(active)
			self.page.event_generate("<Return>")

	# ----------------------------------------------------------------------
	def findUniqueName(self, obj, name=None):
		last = False
		if obj == Input.Card.BODY:
			func = self.viewers.body
			if name is None:
				last = True
				name = self._lastBodyName
		elif obj == Input.Card.REGION:
			func = self.viewers.region
			if name is None:
				last = True
				name = self._lastRegionName
		elif obj == Input.Card.OBJECT:
			func = self.viewers.object
		else:
			return None

		# Find unique name
		if last:
			guess = "%s1"%(name[:7])
		else:
			guess = name[:8]
		try:
			func(guess,"id")
			n = 1
			# find number in the end of the name
			pat = _LAST_NUM.match(name)
			if pat:
				name = pat.group(1)
				n = int(pat.group(2))+1
			while True:
				# Constraint in 8 characters
				sn = str(n)
				guess = "%s%s"%(name[:8-len(sn)],sn)
				try:
					func(guess,"id")
					n += 1
				except KeyError:
					break
		except KeyError:
			pass
		return guess

	# ----------------------------------------------------------------------
	def changeName(self, oldname, newname, card=None):
		if newname == "": return
		if card is None:
			card = self.editObject.cards[0]

		if card.type() in (Input.Card.BODY, Input.Card.REGION):
			pat = _LAST_NUM.match(newname)
			if pat: prefix = pat.group(1)
			else: prefix = newname

			count = 0
			if card.type() == Input.Card.REGION:
				ctype = "ri"
				stype = "Region"
				# Store last name
				self._lastRegionName = prefix
				# Count duplicates
				for c in self.project.input.cards["REGION"]:
					if c.name() == oldname:
						count += 1
			else:
				ctype = "bn"
				stype = "Body"
				# Store last name
				self._lastBodyName = prefix
				# Count duplicates
				for body in Input.BODY_TAGS:
					for c in self.project.input[body]:
						if c.name() == oldname:
							count += 1

			if oldname!="" and count>1:
				change = tkFlair.askyesno("Name changed",
					    "%s name has changed\n"
					    "%s \u2192 %s\n" \
					    "Change all references in input?" \
					    % (stype, oldname, newname), parent=self.page)
			else:
				change = (oldname!="")

			if change:
				undoinfo = undo.createListUndo(
					[ (self.flair.refreshUndo,),
					   self.flair.changeNameUndo(ctype, oldname, newname),
					   self.flair.setWhatUndo(card,0,newname),
					  (self.flair.refreshUndo,) ])
				self.objListFill()
			else:
				undoinfo = self.flair.setWhatUndo(card,0,newname)

			self.flair.addUndo(undoinfo)

		else:
			self.flair.addUndo(self.flair.setWhatUndo(card,0,newname))
		self.project.input.clearCache()

	# ----------------------------------------------------------------------
	def _findEnd(self):
		if "END" not in self.project.input.cards:
			if not self.project.input.cardlist:
				# create the bare minimum input
				self.project.input.minimumInput()
				self.flair.redraw()
			else:
				self.flair.notify("Incomplete GEOMETRY card structure",
					"The basic geometry cards are not found in your input\n" \
					"GEOBEGIN, END, END, GEOEND",
					tkFlair.NOTIFY_ERROR)
		ends = self.project.input["END"]
		if len(ends) != 2:
			self.flair.notify("Incorrect GEOMETRY card structure",
				"The basic geometry cards are not found in your input\n" \
				"GEOBEGIN, END, END, GEOEND",
				tkFlair.NOTIFY_ERROR)
		return ends

	# ----------------------------------------------------------------------
	# Add a body at initial location, with size and z-direction
	# ----------------------------------------------------------------------
	def addBody(self, card):
		"""Add a new body of type card"""
		self.propTab.changePage(_PROPTABS[0])

		end = self._findEnd()
		try:
			pos = end[0].pos()
		except:
			pos = None

		undoinfo = [(self.flair.refreshUndo,)]
		undoinfo.append(self.flair.addCardUndo(card,pos,True))
		undoinfo.append((self.flair.refreshUndo,))
		self.addUndo(undo.createListUndo(undoinfo, "Add body %s"%(card.tag)))

		self.objListFill()

		self.select(card, True)		# Select
		card[SELECT] = card.get(SELECT,0) | BIT_VISIBLE	# and make visible

		self.objListColors()

		return card

	# ----------------------------------------------------------------------
	def addRegion(self, edit=True, mat=None, expr=None):
		"""Add a new region"""
		self.propTab.changePage(_PROPTABS[0])

		end = self._findEnd()
		try:
			pos = end[1].pos()
		except:
			pos = None
		name = self.findUniqueName(Input.Card.REGION)

		undoinfo = [(self.flair.refreshUndo,)]

		# Create REGION card
		card = Input.Card("REGION", [name])
		if expr: card.setExtra(expr)
		undoinfo.append(self.flair.addCardUndo(card, pos, True))

		# Create ASSIGNMA card
		pos = self.project.input.bestPosition("ASSIGNMA")

		if mat is None: mat = self._lastMaterial
		undoinfo.append(self.flair.addCardUndo(
			Input.Card("ASSIGNMA", ["",mat, name]),
				pos, True))

		undoinfo.append((self.flair.refreshUndo,))
		self.addUndo(undo.createListUndo(undoinfo, "Add region"))

		self.loadProject() # FIXME
		self.flair.refresh("card")

		self.select(card, edit)
		if edit: self.editName()
		return card

	# ----------------------------------------------------------------------
	def addRotdefi(self):
		"""Add a ROT-DEFI card"""
		self.propTab.changePage(_PROPTABS[0])

		card = Input.Card("ROT-DEFI")
		pos = self.project.input.bestPosition(card.tag)
		undoinfo = [(self.flair.refreshUndo,)]
		undoinfo.append(self.flair.addCardUndo(card, pos, True))
		undoinfo.append((self.flair.refreshUndo,))
		self.addUndo(undo.createListUndo(undoinfo, "Add ROT-DEFI"))

		self.loadProject() # FIXME
		self.objListFill()

		self.select(card, True)					# Select
		card[SELECT] = card.get(SELECT,0) | BIT_VISIBLE	# and make visible
		self.objListColors()

		self.propListEdit()
		self.editName()
		return card

	# ----------------------------------------------------------------------
	def addObject(self, card):
		"""Add a new object"""
		self.propTab.changePage(_PROPTABS[0])

		pos = self.project.input.bestPosition(card.tag)
		# Add only card do not refresh the display
		self.addUndo(self.flair.addCardUndo(card,pos,True))

		self.objListFill()

		self.select(card, True)					# Select
		card[SELECT] = card.get(SELECT,0) | BIT_VISIBLE	# and make visible
		self.objListColors()

		self.propListEdit()
		return card

	# ----------------------------------------------------------------------
	# Update object variables through a dictionary
	# ----------------------------------------------------------------------
	def updateObject(self, object):
		self.viewers.object(object[ID], "dict", object.var)
		self.draw()

	# ----------------------------------------------------------------------
	# create bodies from mesh
	# ----------------------------------------------------------------------
	def mesh2Bodies(self, event=None):
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if not card.get(SELECT,0) & BIT_SELECT: continue
			if card.type() != Input.Card.OBJECT: continue
			if card.tag != "!mesh": continue

			oid = card[ID]
			say("Object:", card.name(), oid)
			bodies = self.viewers._geometry.object(oid,"bodies")
			say("Bodies:", bodies)
			for i,what in enumerate(bodies):
				tag = what[0]
				what[0] = "%s%d"%(card.name(),i+1)
				body = Input.Card(tag, what)
				self.project.input.addCard(body)
				say(">>>",body)
				#self.addBody(body)
		self.loadProject()
		self.redraw()

	# ======================================================================
	# Undo
	# ======================================================================
	def addUndo(self, undoinfo): self.flair.addUndo(undoinfo)
	def canUndo(self):	return self.flair.canUndo()
	def canRedo(self):	return self.flair.canRedo()

	def undo(self, event=None):
		self.flair.undoInput()

	# ----------------------------------------------------------------------
	def redo(self, event=None):
		self.flair.redoInput()

	# ----------------------------------------------------------------------
	def undolist(self, event=None):
		if tkFlair.UndoListbox.posted is not None:
			tkFlair.UndoListbox.posted.cancel()
		else:
			b = self.tool["undolist"]
			x = b.winfo_rootx()
			y = b.winfo_rooty() + b.winfo_height()
			dlg = tkFlair.UndoListbox(self, self.flair, x, y)

	# ----------------------------------------------------------------------
	def refreshUndo(self):
		"""Called from flair when undo is requested"""
		if not self.created(): return
		if self._undoAfter is not None:
			self.page.after_cancel(self._undoAfter)
		self._undoAfter = self.page.after(UNDOTIMER, self.undoLate)

	# ----------------------------------------------------------------------
	# Force refresh if any pending
	# ----------------------------------------------------------------------
	def forceRefresh(self):
		if self._undoAfter is not None:
			self.page.after_cancel(self._undoAfter)
			self._undoAfter = None
		self.undoLate()

	# ----------------------------------------------------------------------
	# Reload/draw waiting for all undo commands to finalize resizes the action handle
	#----------------------------------------------------------------------
	def undoLate(self):
		self._undoAfter = None
		self.loadProject()	# Load project calls fillRun -> runChange -> redraw
		if self.editObject is not None:
			self.editObject.populateList()
			self.viewers.updateEdit()
		self.redraw()

	# ----------------------------------------------------------------------
	def activateIfNeeded(self, card):
		# Active only cards handled by GeometryEditor
		if card.type() not in (Input.Card.BODY, Input.Card.REGION, Input.Card.OBJECT):
			return

		if self.editObject and card in self.editObject.cards:
			return

		self.select(card, True)

	# ----------------------------------------------------------------------
	# Standard cut, copy and paste
	# ----------------------------------------------------------------------
	def objListCut(self, event=None):
		self.objListCopy()
		self.delete()
		return "break"

	# ----------------------------------------------------------------------
	# Copy selected items to clipboard
	# ----------------------------------------------------------------------
	def objListCopy(self, event=None):
		# Copy also as pickle format
#             sio = StringIO()
		sio = io.BytesIO()
#		sio.write(tkFlair._FLAIRCLIP)
#		sio.write("<card>")
		pickler = pickle.Pickler(sio)
		pickler.dump(tkFlair._FLAIRCLIP)
		pickler.dump("<card>")

		something = False
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.get(SELECT,0) & BIT_SELECT:
				#clipboard.append(card.clone())
				something = True
				card.dump(pickler)
		#clipboard.sort(Input.Card.cmpPos)

		if something:
			self.page.clipboard_clear()
#			self.page.clipboard_append(sio.getvalue())
			self.page.clipboard_append(binascii.b2a_hex(sio.getvalue()), type='STRING')
		return "break"

	#-----------------------------------------------------------------------
	# Deep copy of a region card
	#-----------------------------------------------------------------------
	def regionDeepCopy(self, region, origname, bodies, bodyPos, undoinfo):
		expr = region.extra()

		# 1st clone the bodies
		for item in csg.tokenize(expr):
			if item in ("(",")","-","+","|"): continue
			if item in bodies: continue

			# clone the body
			body = self.findCard("B", item)
			# cloning body body.name()
			card = body.clone()
			name = self.findUniqueName(card.type(),card.name())
			card.setSdum(name)
			self.viewers.geometry().addBody(name, card.tag)
			bodies[item] = card
			undoinfo.append(self.flair.addCardUndo(card, bodyPos, False))
			bodyPos += 1

		# 2nd correct expression
		for name,card in list(bodies.items()):
			expr = re.sub(r"\b%s\b"%(name), card.name(), expr)
		region.setExtra(expr)

		# 3rd add any associated card, ASSIGNMA, BIASING?
		# FIXME will not work properly as we are not renumbering before!!!!
		oldRegion = self.findCard("R",origname)
		if oldRegion:
			pos = self.project.input.bestPosition("ASSIGNMA")
			mat = oldRegion["@material"].name()
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ASSIGNMA", ["",mat, region.name()]),
					pos, False))

		return len(bodies)

	# ----------------------------------------------------------------------
	# Paste cards and return undo information (with no refresh)
	# ----------------------------------------------------------------------
	def _cardsPaste(self, cards, undoinfo, deep=False):
		# find positions
		end = self._findEnd()
		try:
			bodyPos = end[0].pos()
		except:
			bodyPos = self.project.input.bestPosition("REGION")-1
		try:
			regPos  = end[1].pos()
		except:
			regPos  = self.project.input.bestPosition("GEOEND")-1
		objPos  = self.project.input.bestPosition("!object")+1

		pasted   = []

		bodies = {}	# Duplicated bodies
		renamed = ""
		undoinfo.append(self.flair.renumberUndo())
		for card in cards:
			if card.type() == Input.Card.BODY:
				pos = bodyPos

			elif card.type() == Input.Card.REGION:
				pos = regPos

			elif card.type() == Input.Card.OBJECT:
				pos = objPos

			else:
				continue

			# Increase pointers
			if bodyPos >= pos: bodyPos += 1
			if regPos  >= pos: regPos  += 1
			if objPos  >= pos: objPos  += 1

			# Set visibility flag for bodies/objects when duplicate
			visible = False

			origname = card.name()
			name = self.findUniqueName(card.type(),card.name())
			if name != origname:
				renamed += "   %s:\t%s\t-> %s\n"%(card.tag, card.name(), name)
				card.setName(name)
				visible = card.type() in (Input.Card.BODY, Input.Card.OBJECT)

			# Update geometry to register the name, to avoid
			# duplicate names from findUniqueName
			if card.type() == Input.Card.BODY:
				try:
					whats = [card.numWhat(i) for i in range(1,card.nwhats()) ]
				except ValueError:
					whats = []
				card[ID] = self.viewers.addBody(name, card.tag, whats)

			elif card.type() == Input.Card.REGION:
				if deep:
					# Create a deep clone of the card!
					bodiesAdded = self.regionDeepCopy(card, origname, bodies,
									bodyPos, undoinfo)
					bodyPos += bodiesAdded
					if regPos  >= bodyPos: regPos  += bodiesAdded
					if objPos  >= bodyPos: objPos  += bodiesAdded
					if pos     >= bodyPos: pos += bodiesAdded
				card[ID] = self.viewers.geometry().addRegion(name)
			if visible: card[SELECT] = card.get(SELECT,0) | BIT_VISIBLE

			pasted.append(card)
			undoinfo.append(self.flair.addCardUndo(card, pos, False))
		undoinfo.append(self.flair.renumberUndo())

		self.viewers.derive()

		return pasted, renamed

	# ----------------------------------------------------------------------
	# Paste clipboard items to input
	# ----------------------------------------------------------------------
	def objListPaste(self, event=None):
		cards = self.flair.clipboard2Cards()
		if not cards:
			if self.editType() == Input.Card.REGION:
				return self.propListPaste(event)
			return "break"

		self.selectionClear(False)
		self.viewers.stopThreads()

		undoinfo = [self.flair.refreshUndo()]
		pasted,renamed = self._cardsPaste(cards,undoinfo)
		undoinfo.append(self.flair.refreshUndo())

		self.addUndo(undo.createListUndo(undoinfo, "Paste card(s)"))

		if renamed:
			self.flair.notify("Cards renamed",
				"The following cards have been renamed:\n%s"%(renamed),
				tkFlair.NOTIFY_WARNING)

		self.project.input.checkNumbering()	# XXX
		self.project.input.clearCache()

		if pasted:
			self.flair.activateIfNeededInputFrame(pasted[-1])
		return "break"

	# ----------------------------------------------------------------------
	# Paste special makes a deep copy of the objects
	# ----------------------------------------------------------------------
	def objListPasteSpecial(self, event=None):
		cards = self.flair.clipboard2Cards()
		if not cards: return

		# find position to insert the cards
		#unknownpos = True
		#active = self.objList.index(ACTIVE)
		#if active:
		#	activecard = self.objListCard[active]
		#	if activecard:
		#		unknownpos = False
		#		pos = activecard.card.pos()+1

		self.selectionClear(False)
		self.viewers.stopThreads()

		undoinfo = [self.flair.refreshUndo()]
		pasted,renamed = self._cardsPaste(cards, undoinfo, True)
		undoinfo.append(self.flair.refreshUndo())

		self.addUndo(undo.createListUndo(undoinfo, "Paste card(s)"))

# Testing possibility to reload the page immediately and select
#		if self._undoAfter is not None:
#			self.page.after_cancel(self._undoAfter)
		if renamed:
			self.flair.notify("Cards renamed",
				"The following cards have been renamed:\n%s"%(renamed),
				tkFlair.NOTIFY_WARNING)

		self.project.input.checkNumbering()	# XXX
		self.project.input.clearCache()

		# Unfortunately loading of the geometry is done at later state
		# from the refreshUndo() undoLate
		# no possibility to select now

		# if we force the loading we can select
		# but we end up with double loading...
		if pasted:
			self.flair.activateIfNeededInputFrame(pasted[-1])
		return "break"

	# ----------------------------------------------------------------------
	# Delete key was pressed
	# Location     Action
	# --------     --------------------------
	# objList      Delete all selected items (call objListDelete)
	# propList     Delete zones if any (call propListDelete)
	# Viewports    1. If zone selected call propListDelete
	#              2. If bodies are selected delete only the bodies
	#              3. Delete all selected cards
	# ----------------------------------------------------------------------
	def delete(self, event=None):
		focus = self.page.focus_get()

		if focus is self.objList:
			self.objListDelete()

		elif focus is self.propList:
			self.propListDelete()

		else:
			t = False	# type of objects to delete
			# 1st if a region is selected
			if self.editType() == Input.Card.REGION and \
			   self.propList.curselection():
				self.propListDelete()
				return

			# 2nd if objects are selected
			elif self.isAnythingSelected(Input.Card.OBJECT):
				t = Input.Card.OBJECT

			# 3rd if bodies are selected
			elif self.isAnythingSelected(Input.Card.BODY):
				t = Input.Card.BODY

			# bodies and regions removed, to correct references!
			bodiesRemoved = []
			regionsRemoved = []

			undoinfo = []
			for card in reversed(self.project.input.cardlist):
				if card.ignore(): continue
				if card.get(SELECT,0) & BIT_SELECT == 0: continue
				if t and card.type() != t: continue
				if not undoinfo:
					# if we don't renumber cards everything gets screwed up
					#undoinfo.append(self.flair.renumberUndo())
					undoinfo.append(self.flair.refreshUndo())

				card[SELECT] = card.get(SELECT,0) & BIT_UNSELECT
				if card.type() == Input.Card.BODY:
					bodiesRemoved.append(card.name())
				elif card.type() == Input.Card.REGION:
					regionsRemoved.append(card.name())

				undoinfo.append(self.flair.delCardUndo(card.pos(),True))

			# Correct references to bodies
			for body in bodiesRemoved:
				undoinfo.extend(self.flair.removeBodyRefUndo(body))

			# Correct references to regions
			for region in regionsRemoved:
				undoinfo.extend(self.flair.removeRegionRefUndo(region))

			if undoinfo:
				#undoinfo.append(self.flair.renumberUndo())
				undoinfo.append(self.flair.refreshUndo())
				self.addUndo(undo.createListUndo(undoinfo, "Delete card(s)"))
			self.refreshUndoButton()

	# ----------------------------------------------------------------------
	# Update region and body properties in geometry
	# ----------------------------------------------------------------------
	def updateGeometry(self):
		self.viewers.updateGeometry()

	# ----------------------------------------------------------------------
	def updateLayers(self, event=None):
		self.setModified()
		self.viewers.updateLayers(self._layers.layer.name)

	# ----------------------------------------------------------------------
	def updateErrors(self, event=None):
		if self._tabSet.getActivePage() == "Errors":
			self._errors.fill()

	# ----------------------------------------------------------------------
	def changePage(self, event):
		page = self._tabSet.getActivePage()
		if page == "Layers":
			self._layers.show(self.project)
		if page == "Errors":
			self._errors.refresh()

	# ----------------------------------------------------------------------
	def showLayer(self, layer):
		self._layers.show(self.project, layer)
		self._tabSet.changePage("Layers")

	# ----------------------------------------------------------------------
	def draw(self):
		"""Force drawing of viewports"""
		if self.page is not None:
			self.viewers.draw()

	# ----------------------------------------------------------------------
	def redraw(self, event=None):
		"""Force redrawing of viewports"""
		if not self.created(): return
		if self._undoAfter: return	# If a reload/redraw is pending quit
		if self._tabSet.getActivePage() == "Layers":
			self._layers.refresh()
		try:
			start = time.time()
			self.viewers.redraw()
			if GeometryViewer._PROFILE:
				#traceback.print_stack()
				say("GeometryEditor.redraw=",time.time()-start)
			self.refreshUndoButton()
		except TclError:
			# XXX if we close the geoviewer then the
			# flair project is not informed containing
			# a self.geoedit pointer to a non Tk object
			pass

	# ----------------------------------------------------------------------
	def setMouseAction(self, action):
		"""Set mouseAction IntVar with the appropriate action"""
		self._mouseAction.set(action)

	# ----------------------------------------------------------------------
	def mouseActionChanged(self, a=False, b=False, c=False):
		"""Inform viewers that mouse Action has changed"""
		self.viewers._setMouseAction(self._mouseAction.get())

	# ----------------------------------------------------------------------
	def toggleSnap2Grid(self, event=None):
		self.snap2Grid.set(not self.snap2Grid.get())

	# ----------------------------------------------------------------------
	def viewLayout(self, name=None):
		if name is None: name = "2x2"
		self.viewers.layout(name)

	# ----------------------------------------------------------------------
	def transformGeometry(self, event=None):
		cardlist = []
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.type() == Input.Card.BODY and card.get(SELECT,0)&BIT_SELECT:
				cardlist.append(card)

		GeometryBody.TransformGeometryDialog(self.page, self.flair, cardlist)

	# ----------------------------------------------------------------------
	def regionVolume(self, event=None):
		self.viewers.stopThreads()

		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.type() == Input.Card.REGION and card.get(SELECT,0)&BIT_SELECT:
				say("\nRegion:", card.name())
				n = 32
				n = 1<<30
				print(n)
				try:
					#for i in range(18):
					for i in range(1):
						vol = self.viewers._geometry.volume(card[ID], n)
						say("Geometry Volume:", vol, n)
						#vol = self.viewers[0].volume(card[ID], n)
						#say("Viewer Volume:", vol, n)
						n <<= 1
				except KeyboardInterrupt:
					pass

	# ----------------------------------------------------------------------
	def loadProjectColors(self):
		self.viewers.loadProjectColors()
		self.viewers.redraw()

	# ----------------------------------------------------------------------
	# Expand macro body to infinite ones
	# FIXME should be moved to GeometryBody
	# ----------------------------------------------------------------------
	def expandBody(self):
		if self.editObject is None: return
		cards = self.getSelection(objtype=Input.Card.BODY)
		if not cards: return
		self.busy()

		def addUniqueBody(card, tag, whats, pos, comment=""):
			n = self.findUniqueName(Input.Card.BODY, card.name())
			self.viewers.geometry().addBody(n, tag)
			whats[0] = n
			body = Input.Card(tag, whats, comment)
			body.copyProperties(card)
			body[SELECT] = body.get(SELECT,0) | BIT_VISIBLE
			undoinfo.append(self.flair.addCardUndo(body, pos, True))

		undoinfo = [(self.flair.refreshUndo,)]
		for body in cards:
			pos  = body.pos()+1

			#if body.tag not in ("RPP","BOX","RCC","TRC","WED","REC","ARB"): continue
			if body.tag == "RPP":
				addUniqueBody(body, "YZP", [None, body.numWhat(1)], pos, body.comment())
				addUniqueBody(body, "YZP", [None, body.numWhat(2)], pos+1)
				addUniqueBody(body, "XZP", [None, body.numWhat(3)], pos+2)
				addUniqueBody(body, "XZP", [None, body.numWhat(4)], pos+3)
				addUniqueBody(body, "XYP", [None, body.numWhat(5)], pos+4)
				addUniqueBody(body, "XYP", [None, body.numWhat(6)], pos+5)

			elif body.tag == "BOX":
				P = body.bodyP()
				X = body.bodyX()
				Y = body.bodyY()
				Z = body.bodyZ()
				PO = P+X+Y+Z

				addUniqueBody(body, "PLA", ([None]+X)+P,  pos, body.comment())
				addUniqueBody(body, "PLA", ([None]+X)+PO, pos+1)

				addUniqueBody(body, "PLA", ([None]+Y)+P,  pos+2)
				addUniqueBody(body, "PLA", ([None]+Y)+PO, pos+3)

				addUniqueBody(body, "PLA", ([None]+Z)+P,  pos+4)
				addUniqueBody(body, "PLA", ([None]+Z)+PO, pos+5)

			elif body.tag in ("RCC", "REC", "TRC"):
				P = body.bodyP()
				Z = body.bodyZ()
				PZ = P + Z

				dire = Z.direction(Input.zero)
				if dire=="X":
					addUniqueBody(body, "YZP", [None, P.x()], pos, body.comment())
					addUniqueBody(body, "YZP", [None, PZ.x()], pos+1)

				elif dire=="Y":
					addUniqueBody(body, "XZP", [None, P.y()], pos, body.comment())
					addUniqueBody(body, "XZP", [None, PZ.y()], pos+1)

				elif dire=="Z":
					addUniqueBody(body, "XYP", [None, P.z()], pos, body.comment())
					addUniqueBody(body, "XYP", [None, PZ.z()], pos+1)

				else:
					addUniqueBody(body, "PLA", ([None]+Z)+P, pos, body.comment())
					addUniqueBody(body, "PLA", ([None]+Z)+PZ, pos+1)

				if body.tag == "RCC":
					if dire=="X":
						addUniqueBody(body, "XCC", [None, P.y(),P.z(),body.what(7)], pos+2)
					elif dire=="Y":
						addUniqueBody(body, "YCC", [None, P.z(),P.x(),body.what(7)], pos+2)
					elif dire=="Z":
						addUniqueBody(body, "ZCC", [None, P.x(),P.y(),body.what(7)], pos+2)
					#else:
						#addUniqueBody(body, "QUA", [None, P.x(),P.y(),body.what(7)], pos+2)
				elif body.tag == "REC":
					if dire=="X":
						# Not correct
						addUniqueBody(body, "XEC",
							[None, P.y(),P.z(),body.bodyX().length(), body.bodyY().length()],
							pos+2)
					elif dire=="Y":
						# Not correct
						addUniqueBody(body, "YEC",
							[None, P.z(),P.y(),body.bodyX().length(), body.bodyY().length()],
							pos+2)
					elif dire=="Z":
						# Not correct
						addUniqueBody(body, "ZEC",
							[None, P.x(),P.y(),body.bodyX().length(), body.bodyY().length()],
							pos+2)
				#elif body.tag == "TRC":

		if len(undoinfo)>1:
			undoinfo.append(undoinfo[0])
			self.flair.addUndo(undo.createListUndo(undoinfo,"Expand bodies"))

		self.refresh()
		self.notBusy()

	# ----------------------------------------------------------------------
	# expand parenthesis from regions
	# FIXME should be replaced with the C++ version
	# ----------------------------------------------------------------------
	def expandRegions(self):
		if self.editObject is None: return
		cards = self.getSelection(objtype=Input.Card.REGION)
		if not cards: return

		undoinfo = [(self.flair.refreshUndo,)]
		self.busy()
		for region in cards:
			undoinfo.append(self.flair.saveCardUndoInfo(region.pos()))
			try:
				exp = csg.tokenize(region.extra())
				csg.exp2rpn(exp)
				csg.rpnorm(exp)
				expnorm = csg.rpn2exp(exp)
			except csg.CSGException:
				self.flair.notify("Expansion Error",
					"Region %s expansion error:\n%s"% \
						(region.sdum(),sys.exc_info()[1]),
					tkFlair.NOTIFY_ERROR)
				continue
			region["@zone"] = list(map(csg.toString, expnorm))
			region_expr = " ".join(expnorm,"")
			region_expr = region_expr.replace("+", " +")
			region_expr = region_expr.replace("-", " -")
			region_expr = region_expr.replace("|", "\n| ")
			region.setExtra(region_expr)

		undoinfo.append(undoinfo[0])
		self.flair.addUndo(undo.createListUndo(undoinfo,"Expand regions"))

		self.setInputModified()
		self.loadProject()
		self.updateProperties()
		self.notBusy()

	# ----------------------------------------------------------------------
	def toggleSplitter(self, event=None):
		self.splitter.toggle()

	# ----------------------------------------------------------------------
	def help(self, event=None):
		Manual.show("F4")

#===============================================================================
def installed():
	"""Return true if geoviewer is install and loaded"""
	return GeometryViewer.geoviewer is not None

def PILInstalled():
	"""Return true if PIL is installed"""
	return GeometryViewer.Image is not None

#-------------------------------------------------------------------------------
def geoviewer_error():
	"""Return geometry viewer errors during import"""
	return GeometryViewer.geoviewer_error

# -----------------------------------------------------------------------------
# Return geoviewer.so version
# -----------------------------------------------------------------------------
def version():
	if GeometryViewer.geoviewer is not None:
		return GeometryViewer.geoviewer.__version__
	else:
		return None

# -----------------------------------------------------------------------------
# Check version against flair
# -----------------------------------------------------------------------------
def checkVersion(flair):
	v = version()
	if v and v != tkFlair.__version__:
		flair.write("WARNING: flair - geoviewer. Version mismatch.\n")
		flair.write("         flair     Version: %s\n"%(tkFlair.__version__))
		flair.write("         geoviewer Version: %s\n"%(v))
		flair.write("         Please use the same version for both packages\n")
		messagebox.showwarning("Version mismatch.",
			"flair Version: %s\n" \
			"geoviewer Version: %s\n" \
			"Please use the same version for both packages" % \
				(tkFlair.__version__, GeometryViewer.geoviewer.__version__),
			parent=flair)

#-------------------------------------------------------------------------------
def configGet():
	"""Read configuration file"""
	global VISIBLE_COLOR, LOCK_COLOR, FREEZE_COLOR, ERROR_COLOR
	if GeometryViewer.geoviewer_error is None: return

	GeometryViewer._BACKGROUND_COLOR= tkFlair.getIColor("geometry.background",
						GeometryViewer._BACKGROUND_COLOR)
	GeometryViewer._SELECT_COLOR    = tkFlair.getIColor("geometry.select",
						GeometryViewer._SELECT_COLOR)
	GeometryViewer._VISIBLE_COLOR   = tkFlair.getIColor("geometry.visible",
						GeometryViewer._VISIBLE_COLOR)
	GeometryViewer._LOCK_COLOR      = tkFlair.getIColor("geometry.lock",
						GeometryViewer._LOCK_COLOR)
	GeometryViewer._FREEZE_COLOR    = tkFlair.getIColor("geometry.freeze",
						GeometryViewer._FREEZE_COLOR)
	GeometryViewer._WIREFRAME_COLOR = tkFlair.getIColor("geometry.wireframe",
						GeometryViewer._WIREFRAME_COLOR)
	GeometryViewer._BODY_BBIN_COLOR = tkFlair.getIColor("geometry.bodybboxincolor",
						GeometryViewer._BODY_BBIN_COLOR)
	GeometryViewer._BODY_BBOUT_COLOR= tkFlair.getIColor("geometry.bodybboxoutcolor",
						GeometryViewer._BODY_BBOUT_COLOR)
	GeometryViewer._ZONE_BB_COLOR   = tkFlair.getIColor("geometry.zonebboxcolor",
						GeometryViewer._ZONE_BB_COLOR)
	GeometryViewer._REGION_BB_COLOR = tkFlair.getIColor("geometry.regionbboxcolor",
						GeometryViewer._REGION_BB_COLOR)
	GeometryViewer._REGION_COLOR    = tkFlair.getIColor("geometry.region",
						GeometryViewer._REGION_COLOR)
	GeometryViewer._ZONE_COLOR      = tkFlair.getIColor("geometry.zone",
						GeometryViewer._ZONE_COLOR)
	GeometryViewer._ERROR_COLOR     = tkFlair.getIColor("geometry.error",
						GeometryViewer._ERROR_COLOR)
	GeometryViewer._VERTEX_COLOR    = tkFlair.getIColor("geometry.vertex",
						GeometryViewer._VERTEX_COLOR)
	GeometryViewer._LATTICE_COLOR   = tkFlair.getIColor("geometry.lattice",
						GeometryViewer._LATTICE_COLOR)
#	GeometryViewer._REGERROR_COLOR  = tkFlair.getIColor("geometry.region.error",
#						GeometryViewer._REGERROR_COLOR)
	GeometryViewer._TITLE_COLOR     = tkFlair.getIColor("geometry.title",
						GeometryViewer._TITLE_COLOR)
	GeometryViewer._GRIDTEXT_COLOR  = tkFlair.getIColor("geometry.gridtext",
						GeometryViewer._GRIDTEXT_COLOR)

	VISIBLE_COLOR   = "#%06X"%(GeometryViewer._VISIBLE_COLOR)
	LOCK_COLOR      = "#%06X"%(GeometryViewer._LOCK_COLOR)
	FREEZE_COLOR    = "#%06X"%(GeometryViewer._FREEZE_COLOR)
	WIREFRAME_COLOR = "#%06X"%(GeometryViewer._WIREFRAME_COLOR)
	ERROR_COLOR     = "#%06X"%(GeometryViewer._ERROR_COLOR)

	section = tkFlair._GEOMETRY_SECTION
	GeometryViewer.LAPTOP       = tkFlair.getBool(section, "laptop",      GeometryViewer.LAPTOP)
	GeometryViewer.INVERTWHEEL  = tkFlair.getBool(section, "invertwheel", GeometryViewer.INVERTWHEEL)
	GeometryViewer.PINZOOM      = tkFlair.getBool(section, "pinzoom",     GeometryViewer.PINZOOM)
	GeometryViewer.CORES        = tkFlair.getInt( section, "cores",	      GeometryViewer.CORES)
	GeometryViewer.SNAPDISTANCE = tkFlair.getInt( section, "snap",	      GeometryViewer.SNAPDISTANCE)
	GeometryViewer.CURSORSIZE   = tkFlair.getInt( section, "cursorsize",  GeometryViewer.CURSORSIZE)
