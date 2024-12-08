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

import os
import re
import pdb
import sys
import math
import time
import string
import platform
import traceback
import binascii
import io

from operator import attrgetter
#from log import say
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
	import tkinter.dnd as Dnd
	import tkinter.font
except ImportError:
	from tkinter import *
	import tkinter.dnd as Dnd
	import tkinter.font as tkFont

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkinter as tk
import undo
import tkTree
import tkFlair
import tkExtra
import Unicode
import tkDialogs

import awari
import Input
import Layout
import Manual
import Ribbon
import Project
import Validate
import Materials
import InputExtra
import FlairRibbon
import GeometryBody

SHIFT_MASK	= 1
CONTROL_MASK	= 4
# Find ALT key based on OS
_ALT_KEY = 8
# check for Vista
if sys.platform == "cygwin":
	s = platform.system().split("-")
	if s[1][0] == "6":
		_ALT_KEY = 0x20000
elif sys.platform == "win32":
	if getattr(sys, "getwindowsversion", None) is not None:
		if sys.getwindowsversion()[0] == 6:
			_ALT_KEY = 0x20000

# Global variables
_ALL            = "*all*"
_PAD            = 1
_INDENT         = 20
_COLUMNS        = 4
_TABPOS         = [0.10, 0.32, 0.56, 0.81, 1.0]
_ICONWIDTH      = 17		# width+1pixel
_MINTAB         = 0.05
_NORMAL_BORDER  = 0
_SEL_BORDER     = 2
_MIN_WIDTH      = 500
_TRIANGLE_DOWN  = " "+Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE

#
_SHOW_INTER        = True
_SHOW_SCALE        = True
_SHOW_PREPROCESSOR = True

_INSERT_COMMENT    = False
_BODY_INSERT       = True

_KEY_TIMEOUT   = 5000		# after how much time to commit the change
_CARDS_TIMEOUT = 1.0		# maximum time to show cards

# Colors
_ACTIVE_COLOR  = "LightYellow"
_COMMENT_COLOR = "Blue"
_DEFINE_COLOR  = "Orange"
_DIFF_COLOR    = "Yellow"
_DISABLE_COLOR = "DarkGray"
_ERROR_COLOR   = "Red"
_HIDDEN_COLOR  = "DarkBlue"
_INVALID_COLOR = "MistyRose"
_LABEL_COLOR   = "#30A030"
_NORMAL_COLOR  = "White"
_SDUM_COLOR    = "Magenta"
_TAG_COLOR     = "DarkRed"
_VALUE_COLOR   = "Black"

# Fonts
_COMMENT_FONT = ("helvetica", "-12")
_TAG_FONT     = ("helvetica", "-12", "bold")
_LABEL_FONT   = ("helvetica", "-10")
_VALUE_FONT   = ("helvetica", "-12")
_FIXED_FONT   = ("courier",   "-12")
_HIDDEN_FONT  = ("helvetica", "-10", "italic")

# ------------------------------------------------------------------------------
# Find word/body limits in the expression
# ------------------------------------------------------------------------------
def findWord(expr, cursor):
	lexpr = len(expr)
	if lexpr == 0: return None, None

	if cursor<0: cursor = 0
	if cursor>=lexpr: cursor=lexpr-1
	begin = cursor
	accept="-+_."
	while begin>=0:
		ch = expr[begin]
		if ch.isalnum() or accept.find(ch)>=0:
			begin -= 1
		else:
			break

	end = cursor
	while end<lexpr:
		ch = expr[end]
		if ch.isalnum() or accept.find(ch)>=0:
			end += 1
		else:
			break
	return (begin+1,end-1)

#===============================================================================
class CardWidgetBase:
	dialogShown = False		# Keep track of dialogs shown
	width       = _MIN_WIDTH

	# ----------------------------------------------------------------------
	def page(self):		return self.canvas.page
	def input(self):	return self.canvas.input
	def flair(self):	return self.canvas.flair

	# ----------------------------------------------------------------------
	# Initialize the canvas metrics for all cards
	# ----------------------------------------------------------------------
	@staticmethod
	def initMetrics(canvas):
		width = canvas.winfo_width()
		if width<_MIN_WIDTH:
			width = _MIN_WIDTH
		CardWidgetBase.width = width
		return width

	# ----------------------------------------------------------------------
	# Last draw time
	# ----------------------------------------------------------------------
	def setModified(self):
		"""set last time modified"""
		self._modified = time.time()

	# ----------------------------------------------------------------------
	# Return true if card is newer than the widget
	# ----------------------------------------------------------------------
	def isModified(self):
		return True

	# ----------------------------------------------------------------------
	# lower end of card
	# ----------------------------------------------------------------------
	def yend(self):
		if self.canvas.active is self:
			self.canvas.itemconfig(self.box, width=_NORMAL_BORDER)
			x1,y1,x2,y2 = self.canvas.bbox(self.box)
			self.canvas.itemconfig(self.box, width=_SEL_BORDER)
		else:
			x1,y1,x2,y2 = self.canvas.bbox(self.box)
		return y2 + _PAD

	# ----------------------------------------------------------------------
	# Update color state
	# ----------------------------------------------------------------------
	def updateState(self, validate=False):
		pass

	# ----------------------------------------------------------------------
	# edit Item
	# ----------------------------------------------------------------------
	def editItem(self, n=None, x=0, y=0, select=True, focus=True):
		self.canvas.focus("")

	# ----------------------------------------------------------------------
	def indent(self):
		pass

	# ----------------------------------------------------------------------
	# update label
	# ----------------------------------------------------------------------
	def updateLabel(self):
		self.page().updateLabel(self.card)

	# ----------------------------------------------------------------------
	# make visible
	# ----------------------------------------------------------------------
	def see(self, item=None):
		canvas = self.canvas

		if item is None:	# Card
			try:
				x1, y1, x2, y2 = canvas.bbox(self.box)
			except TypeError:
				return
		else:
			x,y = canvas.cursorPosition()
			bbox = canvas.bbox(item)
			if y is not None:
				x1 = x
				x2 = x + 16
				y1 = y
				y2 = y + 16
			else:
				if bbox is None: return
				x1, y1, x2, y2 = bbox

		while x2 > canvas.canvasx(0) + canvas.winfo_width():
			old=canvas.canvasx(0)
			canvas.xview('scroll', 1, 'units')
			# avoid endless loop if we can't scroll
			if old == canvas.canvasx(0):
				break

		while y2 > canvas.canvasy(0) + canvas.winfo_height():
			old=canvas.canvasy(0)
			canvas.yview('scroll', 1, 'units')
			if old == canvas.canvasy(0):
				break

		# done in this order to ensure upper-left of object is visible
		while x1 < canvas.canvasx(0):
			old=canvas.canvasx(0)
			canvas.xview('scroll', -1, 'units')
			if old == canvas.canvasx(0):
				break

		while y1 < canvas.canvasy(0):
			old=canvas.canvasy(0)
			canvas.yview('scroll', -1, 'units')
			if old == canvas.canvasy(0):
				break

	# ----------------------------------------------------------------------
	def findItem(self, x, y):
		return None

	# ----------------------------------------------------------------------
	# Position of starting card
	# ----------------------------------------------------------------------
	def fromPos(self):
		return None

	# ----------------------------------------------------------------------
	# Position of ending card
	# ----------------------------------------------------------------------
	toPos = fromPos

	# ----------------------------------------------------------------------
	# Tag of starting card
	# ----------------------------------------------------------------------
	def fromTag(self):
		return None

	# ----------------------------------------------------------------------
	# Tag of ending card
	# ----------------------------------------------------------------------
	toTag = fromTag

	# ----------------------------------------------------------------------
	# Commit changes
	# WARNING: Return True if LAYOUT has changed, else False
	# ----------------------------------------------------------------------
	def commit(self, timeout=False):
		# Nothing to commit
		return False

	# ----------------------------------------------------------------------
	# editNext
	# ----------------------------------------------------------------------
	def editNext(self, select=True):
		pass

	# ----------------------------------------------------------------------
	# editPrev
	# ----------------------------------------------------------------------
	def editPrev(self, select=True):
		pass

#===============================================================================
# Hidden card widget
#===============================================================================
class HiddenCardWidget(CardWidgetBase):
	def __init__(self, canvas, fromCard, toCard, ypos=0):
		self.canvas   = canvas
		self.card     = None
		self.fromCard = fromCard
		self.toCard   = toCard
		self.create(ypos)

	# ----------------------------------------------------------------------
	# create canvas items
	# ----------------------------------------------------------------------
	def create(self, ypos):
		canvas = self.canvas

		X = 2
		Y = ypos
		M = X + CardWidgetBase.width/2

		if self.fromCard is self.toCard:
			txt = "----------- %s : 1 card hidden-----------" % \
				(self.fromCard.tag)
		else:
			txt = "----------- %s ... %s : %d cards hidden -----------" % \
				(self.fromCard.tag, self.toCard.tag,
				 self.toCard.pos() - self.fromCard.pos() + 1 )

		self.label = canvas.create_text(M,Y,text=txt, font=_HIDDEN_FONT,
					fill=_HIDDEN_COLOR, anchor=N)
		(x1,y1,x2,y2) = canvas.bbox(self.label)
		w = CardWidgetBase.width - 2
		self.box = canvas.create_rectangle(
				(_PAD, ypos, max(_MIN_WIDTH,w), y2),
				tags="box",
				fill=_NORMAL_COLOR, width=_NORMAL_BORDER)
		canvas.tag_lower(self.box)

		# Speedup access to CardWidget from boxes in click()
		canvas.addBox(self)
		self.setModified()

	# ----------------------------------------------------------------------
	# refresh hidden card
	# ----------------------------------------------------------------------
	def refresh(self):
		if self.fromCard is self.toCard:
			txt = "----------- %s : 1 card hidden-----------" % \
				(self.fromCard.tag)
		else:
			txt = "----------- %s ... %s : %d cards hidden -----------" % \
				(self.fromCard.tag, self.toCard.tag,
				 self.toCard.pos() - self.fromCard.pos() + 1 )

		self.canvas.itemconfig(self.label, text=txt)

	# ----------------------------------------------------------------------
	# Search for a string inside fields
	# ----------------------------------------------------------------------
	def search(self, findStr, matchCase, edit, pos):
		return False

	# ----------------------------------------------------------------------
	def delComment(self):
		pass

	# ----------------------------------------------------------------------
	# Position of starting card
	# ----------------------------------------------------------------------
	def fromPos(self):
		return self.fromCard.pos()

	# ----------------------------------------------------------------------
	# Position of ending card
	# ----------------------------------------------------------------------
	def toPos(self):
		return self.toCard.pos()

	# ----------------------------------------------------------------------
	# Tag of starting card
	# ----------------------------------------------------------------------
	def fromTag(self):
		return self.fromCard.tag

	# ----------------------------------------------------------------------
	# Tag of ending card
	# ----------------------------------------------------------------------
	def toTag(self):
		return self.toCard.tag

	# ----------------------------------------------------------------------
	# Return a list of cards
	# ----------------------------------------------------------------------
	def cards(self):
		return self.canvas.input.cardlist[self.fromCard.pos() : self.toCard.pos()+1]

	# ----------------------------------------------------------------------
	# Change card state
	# ----------------------------------------------------------------------
	def setEnable(self, enable):
		undoinfo = undo.UndoInfo()
		for card in self.cards():
			undoinfo.append(self.flair().setEnableUndo(card, enable))
		return undoinfo

#===============================================================================
# Basic Card Widget
#===============================================================================
class CardWidget(CardWidgetBase):
	save_edit  = -1
	save_start =  0

	# ----------------------------------------------------------------------
	def __init__(self, canvas, card, ypos=0):
		self.canvas  = canvas
		self.card    = card
		self.layout  = Layout.getLayout(self.card)
		self.save    = (None, None)
		self.create(ypos)

	# ----------------------------------------------------------------------
	def createComment(self,x,y,txt):
		self.comment = self.canvas.create_text(x,y, text=txt,
				fill=_COMMENT_COLOR, font=_COMMENT_FONT,
				anchor=NW)
		(cx,cy,cw,Y) = self.canvas.bbox(self.comment)
		return Y

	# ----------------------------------------------------------------------
	# create canvas items
	# ----------------------------------------------------------------------
	def create(self, ypos):
		canvas = self.canvas
		card   = self.card
		indent = self._indent = card.indent()

		X = 2 + indent * _INDENT
		Y = ypos
		width = float(CardWidgetBase.width)
		DY = 0

		# comment in blue
		if card.comment()!="":
			Y = self.createComment(X,Y,card.comment())
		else:
			self.comment = None
		y = Y

		nr = self.rows()
		self.cell = [(None,None)] * (nr*_COLUMNS)

		# Create the items
		self.card.clearUserInvalid()
		for i in range(nr*_COLUMNS):
			try:
				lbl, dat = self.layout[i]
				desc = ""
			except:
				lbl, dat, desc = self.layout[i]

			(r,c) = divmod(i,_COLUMNS)
			x = X + int(width*canvas._tabs[c])
			if r>0:
				if DY==0: DY = tagy2-tagy1
				y = Y + DY*r

			if i==0:
				tag = card.tag
				if tag and tag[0] in ("#", "$"):
					color = _DEFINE_COLOR
				else:
					color = _TAG_COLOR

				# Add the image
				try: icon = tkFlair.icons[tag]
				except: icon = tkFlair.icons["tag"]
				self.image = canvas.create_image(X,y,
					anchor=NW, image=icon)
				imgy2 = canvas.bbox(self.image)[-1]

				lblWidget = canvas.create_text(X+_ICONWIDTH, y,
					text=card.tag, font=_TAG_FONT,
					fill=color, anchor=NW)
				x1,tagy1,x2,tagy2 = canvas.bbox(lblWidget)
				self.lastWidget = lblWidget

			elif lbl is not None:
				lblWidget = canvas.create_text(x-_PAD-1, y,
					text=lbl, font=_LABEL_FONT,
					tags="Label%d.%d"%(c,indent),
					fill=_LABEL_COLOR, anchor=NE)

			else:
				continue

			# Create data entry
			anchor = NW
			font   = _VALUE_FONT
			color  = _VALUE_COLOR

			dataerror = False
			if isinstance(dat, tuple):
				datatxt = str(self.listValue(dat)) \
					+ _TRIANGLE_DOWN
				if datatxt[0] == '?':
					dataerror = True
			elif callable(dat):
				datatxt = dat(canvas.layout, card)
				if datatxt is None: datatxt = ""
				try:
					datatxt = str(datatxt)
				except:
					datatxt = str(datatxt)
			elif isinstance(dat, int):
				if dat==0: color = _SDUM_COLOR
				datatxt = card.what(dat)
				try:
					datatxt = str(datatxt)
				except:
					datatxt = str(datatxt)
			else:
				datatxt = None

			if datatxt is None:
				datWidget = None
			else:
				datWidget = canvas.create_text(x, y, text=datatxt,
						font=font, fill=color,
						tags="Data%d.%d"%(c,indent),
						anchor=anchor)
				if DY==0:
					x1,y1,x2,y2 = canvas.bbox(datWidget)
					DY = y2-y1

				# Highlight error
				if dataerror and card.enable:
					canvas.itemconfig(datWidget, fill=_ERROR_COLOR)

				# Find last widget
				self.lastWidget = datWidget

			self.cell[i] = (lblWidget, datWidget)

		# background
		w = CardWidgetBase.width - 2
		x1,y1,x2,y2 = canvas.bbox(self.lastWidget)
		self.box = canvas.create_rectangle(
				(_PAD, ypos, max(_MIN_WIDTH,w), max(y2,tagy2,imgy2)),
				tags="box",
				fill=_NORMAL_COLOR, width=_NORMAL_BORDER)
		canvas.tag_lower(self.box)

		# Speedup access to CardWidget from boxes in click()
		canvas.addBox(self)
		if not card.enable or card.invalid:
			self.updateState()

		self.setModified()

	# ----------------------------------------------------------------------
	# find list value
	# ----------------------------------------------------------------------
	def listValue(self, dat):
		(w, lst) = dat
		if callable(w):
			# Call function to return result
			v = vn = w(self.canvas.layout, self.card)
		else:
			# Get it from card directly
			vn = self.card.what(w)			# Actual what
			v  = self.card.evalWhat(w,True)		# Evaluated

		# Dynamic lists
#		print ("what w ",w," v ",v," vn ",vn);
		if callable(lst):
			lst = lst(self.canvas.layout)

		if not isinstance(lst, list): return vn

		# Convert to integer if necessary
		try:
			vf = float(v)
			vi = int(vf)
			if vf == float(vi): v = vi
		except ValueError:
			pass

		# Check also length 8
		try:
			v8 = v[:8]
		except TypeError:
			v8 = None

		# Make first a fast check
		if v in lst or (v8 is not None and v8 in lst):
			return vn
		else:
			# check for tuples in list
			for item in lst:
				if isinstance(item,tuple) and \
					(v in item or (v8 is not None and v8 in item)):
						return item[0]

		if v!="" or len(lst)==0:
			if self.card.enable:
				msg = "%s do not exist" % (v)
				try:
					int(w)
					self.card.addUserInvalid(w)
				except:
					self.card.addUserInvalid(msg)
#				self.flair().write("\n%s\n%s" % (msg,str(self.card).encode("ascii","replace")))
				self.flair().write(f"{msg}{self.card}\n")
			return "?%s"%(vn)

		item = lst[0]
		if isinstance(item, tuple):
			return item[0]
		else:
			return item

	# ----------------------------------------------------------------------
	# redraw card
	# force re-generation of the card
	# ----------------------------------------------------------------------
	def redraw(self):
		canvas = self.canvas

		fillcolor = canvas.itemcget(self.box, "fill")
		canvas.itemconfig(self.box, width=_NORMAL_BORDER)
		x1,y1,x2,y2 = canvas.bbox(self.box)
		oldHeight = y2-y1
		ax1, ay1, ax2, ay2 = canvas.bbox("all")

		# Delete everything
		if self.comment: canvas.delete(self.comment)
		for (lblWidget, datWidget) in self.cell:
			if lblWidget is not None: canvas.delete(lblWidget)
			if datWidget is not None: canvas.delete(datWidget)
		canvas.delete(self.image)
		canvas.delBox(self)
		canvas.delete(self.box)

		# Mark below as moveable
		canvas.dtag("move")
		canvas.addtag_overlapping("move", ax1, y2, ax2, ay2)

		# Create new box and move everything
		# For CYGWIN use y1+1
		self.create(y1)

		x1,y1,x2,y2 = canvas.bbox(self.box)
		canvas.move("move", 0, (y2-y1)-oldHeight)
		canvas.itemconfig(self.box, fill=fillcolor)

		# Activate box
		if canvas.active is self:
			canvas.focus_set()
			canvas.itemconfig(self.box,
					fill=_ACTIVE_COLOR,
					width=_SEL_BORDER)
					# FIXME check adjustTextBox

			canvas.dtag("active","active")
			canvas.addtag_withtag("active", self.box)

	# ----------------------------------------------------------------------
	# refresh card widget
	# ----------------------------------------------------------------------
	def refresh(self):
		canvas = self.canvas

		# Comment
		if self.card.comment()!="":
			if self.comment is None:
				self.insComment()
			self.adjustTextBox(False, -1)
			canvas.itemconfig(self.comment,
				text=self.card.comment())
			self.adjustTextBox(True, -1)
		else:
			self.delComment()

		# Go through all fields
#		self.card.clearUserInvalid()
		for i in range(len(self.cell)):
			lblWidget, datWidget = self.cell[i]
			if datWidget is None: continue
			dat = self.layout[i][1]

			if isinstance(dat, tuple):
				datatxt = str(self.listValue(dat)) \
					+ _TRIANGLE_DOWN
				# Highlight errors
				if datatxt[0]=='?' and self.card.enable:
					canvas.itemconfig(datWidget,
						fill=_ERROR_COLOR)
			elif callable(dat):
				datatxt = dat(canvas.layout, self.card)
			else:
				if dat == -1: self.adjustTextBox(False, i)
				datatxt = str(self.card.what(dat))

			if datatxt is None: datatxt = ""
			canvas.itemconfig(datWidget, text=datatxt)
			if dat == -1: self.adjustTextBox(True, i)

		self.updateState()
		self.setModified()

	# ----------------------------------------------------------------------
	# correct indent level of the card if needed
	# ----------------------------------------------------------------------
	def indent(self):
		if self._indent != self.card.indent():
			self.redraw()

	# ----------------------------------------------------------------------
	# rows needed (excluding the comment line)
	# ----------------------------------------------------------------------
	def rows(self):
		if self.layout is Layout._generic_layout:
			if len(self.card.whats())>7:
				return 4
			else:
				return 2
		else:
			return len(self.layout)//_COLUMNS

	# ----------------------------------------------------------------------
	# change the layout
	# ----------------------------------------------------------------------
	def changeLayout(self, layout):
		if self.layout is layout: return False

		# Compare layouts to see if they are identical
		if len(self.layout) == len(layout):
			equal = True
			for i in range(len(self.layout)):
				if self.layout[i][0] != layout[i][0]:
					equal = False
					break
			if equal: return False

		self.layout = layout
		self.redraw()
		self.canvas._updateScrollBars()
		return True

	# ----------------------------------------------------------------------
	# highlight label
	# ----------------------------------------------------------------------
	def highlightLabel(self, what=None):
		self.page().highlightLabel(what)

	# ----------------------------------------------------------------------
	# Return true if card is newer than the widget
	# ----------------------------------------------------------------------
	def isModified(self):
		return self.card._modified > self._modified

	# ----------------------------------------------------------------------
	# Insert comment
	# ----------------------------------------------------------------------
	def insComment(self):
		if self.comment: return
		self.adjustTextBox(False)
		x1,y1,x2,y2 = self.canvas.bbox(self.box)
		self.createComment(2+self.card.indent()*_INDENT, y1+_PAD, "Comment")
		self.adjustTextBox(True)
		self.setModified()

	# ----------------------------------------------------------------------
	# Delete comment
	# ----------------------------------------------------------------------
	def delComment(self):
		if self.comment is None: return
		self.adjustTextBox(False)
		self.canvas.delete(self.comment)
		self.comment = None
		self.adjustTextBox(True)
		self.setModified()

	# ----------------------------------------------------------------------
	# move everything below comment to adjust properly the
	# size necessary for the comment or extra field
	# ----------------------------------------------------------------------
	def adjustTextBox(self, adjust, edit=None):
		canvas = self.canvas
		if edit is None:
			edit = canvas.edit

		if edit<0:
			h = 0
			if self.comment is not None:
				try:
					x1,y1,x2,y2 = canvas.bbox(self.comment)
					h = y2-y1
				except:
					pass
		else:
			bb = canvas.bbox(self.cell[edit][1])
			if bb:
				x1,y1,x2,y2 = bb
			else:
				return False
			h = y2-y1

		if not adjust:
			self._before = h
			return False

		# some times _before is not set FIXME
		try: dist = h - self._before
		except: return False

		if dist == 0: return False

		ax1, ay1, ax2, ay2 = canvas.bbox("all")
		bx1, by1, bx2, by2 = canvas.bbox(self.box)

		canvas.dtag("move")
		canvas.addtag_overlapping("move", ax1, by1+1, ax2, ay2-1)

		# untag current box,image and comment
		canvas.dtag(self.box, "move")
		if edit>=0: canvas.dtag(self.image, "move")
		if self.comment is not None:
			canvas.dtag(self.comment, "move")
		if edit>=0:
			for i in range(len(self.cell)):
				(l,d) = self.cell[i]
				if i<=edit:
					if l: canvas.dtag(l, "move")
					if d: canvas.dtag(d, "move")

		canvas.move('move', 0, dist)
		# XXX why +/- 1 it comes from the width!!!
		canvas.coords(self.box, bx1+(_SEL_BORDER-1),
					by1+(_SEL_BORDER-1),
					bx2-(_SEL_BORDER-1),
					by2+dist-(_SEL_BORDER-1))
		if edit<0:
			self.see(self.comment)
		else:
			self.see(self.cell[edit][1])
		canvas._updateScrollBars()
		return True

	# ----------------------------------------------------------------------
	# check for a possible name change in the card
	# ----------------------------------------------------------------------
	def changeName(self, old, new):
		if old=="" or old is None: return []
		tag = self.card.tag

		# Check card type
		if len(tag)==3 and self.card.isGeo():
			stype = "Body"
			ctype = "bn"
			lst   = self.canvas.layout.bodyList(False)
		elif tag == "#define":
			stype = "#define"
			ctype = "vi"
			lst   = self.canvas.layout.defineList()
		elif tag == "REGION":
			stype = "Region"
			ctype = "ri"
			lst   = self.canvas.layout.regionList()
		elif tag == "MATERIAL":
			stype = "Material"
			ctype = "mi"
			lst   = self.canvas.layout.materialList()
		elif tag in ("EVENTBIN", "RESNUCLE", "USRBDX", "USRBIN",
			     "USRCOLL", "USRTRACK", "USRYIELD"):
			stype = "Detector"
			ctype = "di"
			lst   = self.canvas.layout.detectorList(tag)
		elif tag == "ROT-DEFI":
			stype = "Transformation"
			ctype = "ti"
			lst   = self.canvas.layout.rotdefiList()
		else:
			return []

		# Check for duplicate (apart from ROT-DEFI where duplicated are permitted!)
		if tag=="ROT-DEFI" or len([x for x in lst if x==old])>1:
			change = tkFlair.askyesno("Name changed",
				    "%s name has changed\n"
				    "%s -> %s\n" \
				    "Change all references in input?" \
				    % (stype, old, new), parent=self.canvas)
		else:
			change = True

		CardWidgetBase.dialogShown = True
		if change:
			self.input().clearCache()
			self.canvas.layout.cache(False)

			undoinfo = [self.flair().changeNameUndo(ctype, old, new)]
			if tag in ("EVENTBIN", "USRBIN"):
				undoinfo.append( self.flair().changeNameUndo("bi", old, new) )
			self.flair().notify("Name change",
				"Changed all references of %s\n%s to %s"%(stype, old, new))

			self.canvas.layout.cache(True)
			self.canvas.refresh()
			self.canvas.layout.cache(False)
			return undoinfo
		else:
			return []

	# ----------------------------------------------------------------------
	# Commit changes
	# WARNING: Return True if LAYOUT has changed, else False
	# timeout = True is the committing is done due to timeout (see _handleKey)
	#           False when committing is due to a user action
	# ----------------------------------------------------------------------
	def commit(self, timeout=False):
		canvas = self.canvas

		# active can be different than self
		active  = canvas.active
		edit    = canvas.edit
		changed = False

		# Comment?
		if edit is None or active is None:
			return False
		elif edit < 0:
			if active.comment is None: return False
			val = canvas.itemcget(active.comment,"text")
			if val != active.save[1]:
				canvas.page.setInputModified()
				try:
					sval = str(val)
				except UnicodeEncodeError:
					sval = val.encode("utf-8")
				canvas.addUndo(canvas.flair.setCommentUndo(active.card, sval))

		elif edit >= len(self.cell):	# layout changed
			pass
		else:
			try: lbl, dat = active.layout[edit]
			except: lbl, dat, desc = active.layout[edit]

			if isinstance(dat, tuple): # List or filename
				self.card.clearUserInvalid()
				widget = active.cell[edit][1]
				ptxt = canvas.itemcget(widget, "text")
				if not self.card.enable:
					canvas.itemconfig(widget,fill=_DISABLE_COLOR)

				if isinstance(dat[1], tuple):
					txt = val = canvas.project().relativePath(InputExtra.filename)
					# if user hits cancel then do nothing
					if txt == "": return
					if active.card.tag == "VOXELS": # XXX strip extension
						# Special case
						fn,ext = os.path.splitext(val)
						if ext == ".vxl":
							txt = val = fn
				else:
					txt,val = canvas.popupList.selection()

				if txt is not None: # Text, Value
					txt = str(txt)
					txt += _TRIANGLE_DOWN
					canvas.itemconfig(widget, text=txt)
					w = dat[0]
					if callable(w):
						if w(canvas.layout, active.card) != val:
							changed = True
							w(canvas.layout, active.card, val)
					else:
						try:
							sval = str(val)
						except UnicodeEncodeError:
							sval = val.encode("ascii","replace")
						if active.card.what(w) != sval:
							changed = True
							canvas.addUndo(
								canvas.flair.setWhatUndo(
									active.card,
									w, sval))

				# Check for a possible change in the layout
				if txt != ptxt:
					# Check for inclusion of a new file
					if active.card.tag == "#include":
# XXX XXX XXX
# FIXME include() should return the cards to include or delete
# and pass it through the addCardUndo, delCardUndo
# XXX XXX XXX
						if self.input().include(active.card):
							self.canvas.redraw()
					elif active.card.tag == "VOXELS":
						if dat[0]==0:	# Special case!!!
							# Load new voxel file
							active.card.loadVoxel()
					else:
						layout = Layout.getLayout(self.card)
						canvas.page.setInputModified()
						if not self.changeLayout(layout):
							self.refresh()
				canvas.highlight()

			else:	# Normal edit text
				(lblWidget, datWidget) = active.cell[edit]
				if datWidget is not None:
					val = canvas.itemcget(datWidget, "text")
					if not self.card.enable:
						canvas.itemconfig(datWidget,
							fill=_DISABLE_COLOR)
					if val != active.save[1]:
						changed = True
						if callable(dat):
							dat(canvas.layout, active.card, val)
						else:
							try:
								sval = str(val)
							except UnicodeEncodeError:
								sval = val.encode("ascii","replace")

							# Check if value didn't change
							if not timeout and dat==0 and val is not None:
								undoinfo = [canvas.flair.setWhatUndo(
											active.card,
											dat, sval)]
								undoinfo.extend(active.changeName(
											active.save[1],
											val))
								canvas.addUndo(undo.createListUndo(undoinfo, undoinfo[0][0]))
							elif active.card.what(dat) != sval:
									canvas.addUndo(canvas.flair.setWhatUndo(
										active.card,
										dat, sval))
						canvas.page.setInputModified()
						if not timeout:
							active.save = (active.save[0], val)

						# if a body name has changed
						if dat==0 and active.card.tag in Input.BODY_TAGS:
							self.input().clearCache("bodies")
						self.input().clearCache(active.card.tag)

			active.updateState(True)
			active.updateLabel()
		return changed

	# ----------------------------------------------------------------------
	# Find selected item from list and inserted as text to the active item
	# at the cursor position. Mark it also as selected. Move the cursor
	# either to the right or the left(unchanged) of the insertion point
	# ----------------------------------------------------------------------
	def insertSelection(self, moveLeft=False):
		canvas = self.canvas
		txt, val = canvas.popupList.selection()
		try: body = txt.split('[')[0]
		except: return
		item   = canvas.focus()
		insert = canvas.index(item, INSERT)
		#canvas.select_from(item, INSERT)
		canvas.pasteText(body)
		#canvas.select_to(item, canvas.index(item, INSERT)-1)
		if moveLeft:
			canvas.icursor(item, insert)
		else:
			if body == Layout._BODY_PARENTHESIS:
				canvas.icursor(item, insert+2)
		# Print commit changes
		self.commit()

	# ----------------------------------------------------------------------
	# Is present item a list
	# ----------------------------------------------------------------------
	#def isEditList(self):
	#	if self.canvas.edit<1:
	#		return False
	#	(lbl,dat) = self.canvas.active.layout[self.canvas.edit]	# FIXME
	#	return isinstance(dat, tuple)

	# ----------------------------------------------------------------------
	# restore value before editing
	# ----------------------------------------------------------------------
	def restore(self):
		canvas = self.canvas
		edit   = canvas.edit
		self.adjustTextBox(False)
		if edit < 0:
			item = self.comment
		else:
			(lbl, item) = self.cell[edit]

		canvas.itemconfig(item, text=self.save[1])
		self.adjustTextBox(True)

	# ----------------------------------------------------------------------
	# Return a list of cards
	# ----------------------------------------------------------------------
	def cards(self):
		return [self.card]

	# ----------------------------------------------------------------------
	# Change card state
	# ----------------------------------------------------------------------
	def setEnable(self, enable):
		undoinfo = self.flair().setEnableUndo(self.card, enable)
		self.updateState()
		return undoinfo

	# ----------------------------------------------------------------------
	# change colors to display card state
	# ----------------------------------------------------------------------
	def updateState(self, validate=False):
		canvas = self.canvas

		err = None
		if self.card.enable:	# Enabled
			comColor = _COMMENT_COLOR
			if self.card.tag[0] in ("#", "$"):
				tagColor = _DEFINE_COLOR
			else:
				tagColor = _TAG_COLOR
			lblColor  = _LABEL_COLOR
			datColor  = _VALUE_COLOR
			sdumColor = _SDUM_COLOR

			# Check for errors in case they exist
			if self.card.invalid or validate:
				case = self.card.case()
				err  = self.card.validate(case)

		else:	# Disabled
			comColor  = _DISABLE_COLOR
			tagColor  = _DISABLE_COLOR
			lblColor  = _DISABLE_COLOR
			datColor  = _DISABLE_COLOR
			sdumColor = _DISABLE_COLOR

		if self.comment is not None:
			canvas.itemconfig(self.comment, fill=comColor)

		for i in range(len(self.cell)):
			try: lbl, dat = self.layout[i]
			except: lbl, dat, desc = self.layout[i]
			(lblWidget,datWidget) = self.cell[i]

			if lblWidget is not None:
				if i==0:
					canvas.itemconfig(lblWidget, fill=tagColor)
				else:
					canvas.itemconfig(lblWidget, fill=lblColor)

			if datWidget is not None:
				if isinstance(dat,int) and dat==0:
					color = sdumColor
				else:
					color = datColor
				if err:
					if isinstance(dat,tuple):
						w = dat[0]
					elif isinstance(dat,int):
						w = dat
					else:
						w = None
					if err[0] is None or w in err:
						color = _ERROR_COLOR
				canvas.itemconfig(datWidget, fill=color)

		if self is not self.canvas.active and self.card.invalid:
			canvas.itemconfig(self.box, fill=_INVALID_COLOR)

	# ----------------------------------------------------------------------
	# find item by position x,y in canvas
	# ----------------------------------------------------------------------
	def findItem(self, x, y):
		# Find region to edit
		bx1,by1,bx2,by2 = self.canvas.bbox(self.box)

		if self.comment:
			cx1,cy1,cx2,cy2 = self.canvas.bbox(self.comment)
			if y<cy2:
				return -1
		else:
			cy2 = by1

		# Find column
		frac = float(x - bx1) / float(bx2 - bx1)
		for col in range(_COLUMNS):
			if self.canvas._tabs[col] <= frac < self.canvas._tabs[col+1]: break
		else:
			col = 0

		if self.lastWidget:
			lx1,ly1,lx2,ly2 = self.canvas.bbox(self.lastWidget)
			nr1  = self.rows()-1
			if nr1==0 or y >= ly1:
				row = nr1
			else:
				h   = (ly1 - cy2) // nr1
				row = (y   - cy2) // h
		else:
			nr  = self.rows()
			h   = (by2-cy2) // nr
			row = (y - cy2) // h

		if row>=0 or col>=0:
			i = row*_COLUMNS + col
			if i<0:
				i=0
			elif i>=len(self.cell):
				i = len(self.cell)-1

			closest = self.canvas.find_closest(x,y)[0]
			if i == 0:	# Special treatment on tag
				if closest in (self.cell[0][0], self.image):
					return None
			for ii,cell in enumerate(self.cell):
				if cell and closest in cell:
					return ii
			#say("Empty cell",i)
			if i<len(self.cell) and self.cell[i] is not None: return i
			return None
		return None

	# ----------------------------------------------------------------------
	# edit Item
	# ----------------------------------------------------------------------
	def editItem(self, n=None, x=0, y=0, select=True, focus=True):
		if self.canvas.edit != n:
			if self.canvas.edit is not None:
				self.canvas.selectClear()

		lst = None
		color = _VALUE_COLOR
		if n is None:
			datWidget = ""
		elif n<0:
			dat = 0
			datWidget = self.comment
		elif n>=len(self.cell):	# if layout changed
			dat = 0
			datWidget = ""
		else:
			(lblWidget, datWidget) = self.cell[n]
			if datWidget is not None:
				dat  = self.layout[n][1]
				what = None
				if isinstance(dat, int):
					what = dat
					if dat==0: color = _SDUM_COLOR
				elif isinstance(dat, tuple):
					# create an option menu
					(var, lst) = dat
					if isinstance(var, int): what = var
					if callable(lst):
						lst = lst(self.canvas.layout)
					w = dat[0]
					if callable(w):
						sel = w(self.canvas.layout, self.card)
					else:
						sel = self.card.what(w)
				self.highlightLabel(what)
			else:
				datWidget = ""
				n = None

		if datWidget is None: datWidget = ""

		self.canvas.edit = n
		if datWidget != "":
			self.canvas.highlight(datWidget)

		if n is not None:
			if n>=0 and datWidget != "":
				self.canvas.itemconfig(datWidget, fill=color)

		self.canvas.multiline = False
		if lst:
			# save value
			self.save = (self.canvas.edit, self.canvas.itemcget(datWidget, "text")[:-2])
			if focus:
				if isinstance(lst, list):
					self.canvas.focus("")
					self.see(datWidget)
					self.canvas.postList(lst, sel)
				elif isinstance(lst, tuple):
					# Opening bFiledialog FocusOut will send a commit
					# and reset the canvas.edit to None
					InputExtra.openFile(self.canvas,self.save[1],lst)
					# Reset canvas.edit and commit the change
					self.canvas.edit = n
					self.commit()
				else:
					raise Exception("Invalid layout defined "+self.cell[n])
		else:
			# save value
			self.save = (self.canvas.edit, self.canvas.itemcget(datWidget, "text"))
			if focus:
				self.canvas.focus_set()
				self.canvas.focus(datWidget)
			if self.canvas.edit is not None:
#				if n<0 or dat<0:
				if n<0:
					self.canvas.multiline = True
				if isinstance(dat, int):
                                        mydum = dat
				elif isinstance(dat, tuple):
                                        mydum=dat[0]        
				else:
                                        mydum=999
				if mydum<0:
					self.canvas.multiline = True

				self.canvas.icursor(datWidget, "@%d,%d" % (x, y))
				self.see(datWidget)
				if select:
					self.canvas.selectAll()
				else:
					self.canvas.selectClear()

	# ----------------------------------------------------------------------
	# editNext
	# ----------------------------------------------------------------------
	def editNext(self, select=True):
		oldEdit = self.canvas.edit
		n = oldEdit
		if n is None: n = -2
		lupe = 0
		while 1:
			n += 1
			if n >= len(self.cell): # or n>len(self.card.whats())+1:
				if self.comment:
					n = -1
					break
				else:
					n = 0
				lupe += 1
				if lupe>2: return

			if n>=0:
				if self.cell[n][1]: break
			else:
				if self.comment: break

		if oldEdit != n:
			self.editItem(n, select)
		else:
			self.editItem()

	# ----------------------------------------------------------------------
	# editPrev
	# ----------------------------------------------------------------------
	def editPrev(self, select=True):
		oldEdit = self.canvas.edit
		n = oldEdit
		if n is None: n = -1
		if self.comment:
			low = -1
		else:
			low = 0
		lupe = 0
		while 1:
			n -= 1
			if n < low:
				n = len(self.cell)-1
				lupe += 1
				if lupe>2: return
			if self.cell[n][1]: break

		if oldEdit != n:
			self.canvas.edit += 1
			self.editItem(n, select)
		else:
			self.editItem()

	# ----------------------------------------------------------------------
	# Search for a string inside fields
	# if found activate the card
	# ----------------------------------------------------------------------
	def search(self, findStr, matchCase, edit, start):
		#if self.card.name() == "BA18R__d":
		#	import pdb; pdb.set_trace()
		if self is not self.canvas.active:
			edit = None
			pos  = 0
		else:
			# remember last found...
			if edit is None:
				edit  = CardWidget.save_edit
				start = CardWidget.save_start
			else:
				CardWidget.save_edit  = edit
				CardWidget.save_start = start

		# Start search in a circular order
		# comment -> cells [label -> data] ...
		if edit is None or edit<0:
			txt = self.card.comment()
			if not matchCase: txt = txt.upper()
			pos = txt.find(findStr, start)
			if pos>=0:
				self.canvas.activate(self)
				self.see()
				self.canvas.highlight(self.comment)
				self.canvas.select_from(self.comment, pos)
				self.canvas.select_to(self.comment, pos+len(findStr)-1)
				CardWidget.save_edit  = -1
				CardWidget.save_start = pos+len(findStr)
				return True

		if edit is None:
			startItem = 0
		else:
			startItem = edit

		for i in range(startItem, len(self.cell)):
			(lblWidget, datWidget) = self.cell[i]
			if lblWidget:
				txt = self.canvas.itemcget(lblWidget,"text")
				if not matchCase: txt = txt.upper()
				pos = txt.find(findStr, start)
				if pos>=0:
					self.canvas.activate(self)
					self.see()
					self.canvas.highlight(lblWidget)
					self.canvas.select_from(lblWidget, pos)
					self.canvas.select_to(lblWidget, pos+len(findStr)-1)
					CardWidget.save_edit = i
					CardWidget.save_start = pos+len(findStr)
					return True

			if datWidget:
				txt = self.canvas.itemcget(datWidget,"text")
				if not matchCase: txt = txt.upper()
				pos = txt.find(findStr, start)
				if pos>=0:
					self.canvas.activate(self)
					self.see()
					self.canvas.highlight(datWidget)
					self.canvas.select_clear()
					self.canvas.select_from(datWidget, pos)
					self.canvas.select_to(datWidget, pos+len(findStr)-1)
					CardWidget.save_edit = i
					CardWidget.save_start  = pos+len(findStr)
					return True
		return False

	# ----------------------------------------------------------------------
	# Position from starting card
	# ----------------------------------------------------------------------
	def fromPos(self):
		return self.card.pos()

	# ----------------------------------------------------------------------
	# Position from ending card
	# ----------------------------------------------------------------------
	toPos = fromPos

	# ----------------------------------------------------------------------
	# Tag of starting card
	# ----------------------------------------------------------------------
	def fromTag(self):
		return self.card.tag

	# ----------------------------------------------------------------------
	# Tag of ending card
	# ----------------------------------------------------------------------
	toTag = fromTag

#===============================================================================
# Input Canvas filled with CardWidget
#===============================================================================
class InputCanvas(Canvas):
	def __init__(self, master, page, *arg, **kw_args):
		Canvas.__init__(self, master, *arg, **kw_args)
		# Master link
		self.page     = page
		self.flair    = page.flair

		# Interrupt long process
		self._interrupt = False

		# Undo/Redo list
		self.lastKey  = (-1,0)

		# Layout attached to canvas
		self.input  = None
		#self.layout = Layout.LayoutManager(self.flair)
		self.layout = self.flair._layout

		# Card list
		self.cardlist = []
		self.boxes    = {}	# Dictionary of boxes
					# for faster access to CardWidgets
					# in click()
		self._resizeWidth = 0
		self.setTabs()
		self._highlightItem = None

		# Active field
		self.active     = None	# Active card widget
		self.edit       = None	# editing field in active card
		self.multiline  = False	# allow Return to open a new line
		self.selAnchor  = None	# Card selection range
		self.selFrom    = -1
		self.selTo      = -1
		self.selEditAnchor = -1	# Edit item anchor
		self.shift      = False	# Shift Key
		self.dnd_item   = None	# Drag item
		self.dnd_time   = None	# Initialize drag time
		self.dnd_active = None	# Old active
		self.scandrag   = False	# Scan drag
		self.delay      = tkExtra.Balloon.delay	# Delay from Balloon

		# Popup list
		self.popupList  = InputExtra.SelectionList(self)

		# Canvas bindings
		self.bind('<Button-1>',  self.click)
		self.bind('<Double-1>',  self.double)
		#self.bind('<Button-2>',  self.paste)
		self.bind('<Button-3>',  self.popupMenu)
		self.bind('<Control-Button-3>',  self.popupInsertMenu)
		self.bind('<Button-4>',   lambda e,s=self: s.yview(SCROLL, -1, UNITS))
		self.bind('<Button-5>',   lambda e,s=self: s.yview(SCROLL,  1, UNITS))
		self.bind('<ButtonRelease-2>', self.scanRelease)
		self.bind('<MouseWheel>', lambda e,s=self: s.yview(SCROLL,
						-e.delta//120, UNITS))

		self.bind('<Key>',	self._handleKey)
		self.bind('<Tab>',	self.tabKey)
		self.bind('<Shift-Tab>',self.tabShiftKey)
		if sys.platform == "linux2":
			self.bind('<ISO_Left_Tab>', self.tabShiftKey)

		self.bind('<Control-Key-slash>',self.selectAll)

		# Ignore events are handled by the class_bind in flair.py
		self.bind('<Control-Key-Return>',self.popupInsertMenu)
		self.bind('<Control-Key-Up>',	self.moveUp)
		self.bind('<Control-Key-Down>',	self.moveDown)
		self.bind('<Control-Key-space>',self.popupMenu)

		self.bind("<Insert>",		self.page.add)

		self.bind('<Control-Key-a>',	self.selectAll)
		self.bind('<Control-Key-A>',	self.unSelectAll)
#		self.bind('<Control-Key-d>',	self.clone)
#		self.bind('<Control-Key-D>',	self.clone)
		self.bind("<Control-Key-f>",	lambda e,s=self : s.page.findString.focus_set())
		self.bind("<Control-Key-F>",	lambda e,s=self : s.page.findString.focus_set())
		self.bind("<Control-Key-e>",	self.page.editCard)
		self.bind("<Control-Key-g>",	self.page.find)
		self.bind("<Control-Key-G>",	self.page.find)
#		self.bind("<Control-Key-h>",	self.replace)
#		self.bind("<Control-Key-H>",	self.replace)
		self.bind('<Control-Key-i>',	self.toggleShow)
		self.bind('<Control-Key-l>',	self.toggleState)
		self.bind('<Control-Key-L>',	self.toggleState)
		self.bind("<Control-Key-p>",	self.page.hardcopy)
		self.bind("<Control-Key-P>",	self.page.hardcopy)
		self.bind("<Control-Key-t>",	self.page.transformGeometry)
		self.bind("<Control-Key-T>",	self.page.transformGeometry)

		self.bind("<F11>",		self.page.toggleSplitter)

		self.bind('<Configure>',	self.resize)
		self.bind('<FocusOut>',		self._focusOut)
		self.bind('<Leave>',		self.hideBalloon)
		self.motionBind()

	# ----------------------------------------------------------------------
	def project(self): return self.flair.project

	# ----------------------------------------------------------------------
	def setTabs(self):
		if self.project() is None:
			self._tabs = _TABPOS[:]
		else:
			self._tabs = self.project().tabs

	# ----------------------------------------------------------------------
	# There is a bug/conflict with dnd motion bind!
	# ----------------------------------------------------------------------
	def motionBind(self):
		self.motionId   = -1
		self.motionTime = 0.0
		self.bind('<Motion>',	 self.motion)
		self.bind('<B1-Motion>', self.mouseSelect)
		self.bind('<B2-Motion>', self.scanDrag)

	# ----------------------------------------------------------------------
	# update input
	# ----------------------------------------------------------------------
	def updateInput(self):
		self.active = None
		self.edit   = None
		self.input  = self.project().input
		self.layout.updateInput()

	# ----------------------------------------------------------------------
	def ignore(self, event):
		pass

	# ----------------------------------------------------------------------
	# CANVAS DRAWING
	#
	# Return item position in widget coordinates
	# ----------------------------------------------------------------------
	# ----------------------------------------------------------------------
	# Move the columns inside the canvas
	# ----------------------------------------------------------------------
	def resize(self, event=None):
		self.popupList.unPost()
		width = CardWidgetBase.initMetrics(self)

		if event is not None and self._resizeWidth == width: return
		self._resizeWidth = width

		if len(self.cardlist)==0: return

		# Scale box
		bbox = self.bbox(self.cardlist[0].box)
		xscale = float(width) / float(bbox[2]-bbox[0])
		self.scale("box", _PAD, 0, xscale, 1.0)

		# Move columns
		for indent in range(20):
			xind = indent*_INDENT
			found = False
			xind += 2
			for c in range(_COLUMNS):
				x = xind + int(float(width)*self._tabs[c])

				tag = "Label%d.%d" % (c,indent)
				bbox = self.bbox(tag)
				if bbox:
					self.move(tag,(x-1)-bbox[2], 0)
					found = True

				tag = "Data%d.%d" % (c,indent)
				bbox = self.bbox(tag)
				if bbox:
					self.move(tag,x-bbox[0], 0)
					found = True
			if not found: break

		if self.edit is not None:
			item = self.focus()
			self.highlight(item)

		self.highlight()

	# ----------------------------------------------------------------------
	def yview(self, *args):
		self.popupList.unPost()
		Canvas.yview(self, *args)

	# ----------------------------------------------------------------------
	# Find cursor position inside the focused item
	# ----------------------------------------------------------------------
	def cursorPosition(self):
		item = self.focus()
		if item:
			font = tkinter.font.Font(name=self.itemcget(item,"font"), exists=True)
			x1,y1,x2,y2 = self.bbox(item)
			fontheight = font.metrics().get("linespace")
			index = self.index(item,INSERT)
			if index==0:
				return x1,y1+fontheight//2
			txt = self.itemcget(item,"text")
			# count the lines before
			subtxt = txt[:index]
			lines = subtxt.count("\n")
			y = y1 + lines*fontheight + fontheight//2
			if lines>0:
				x = x1 + font.measure(subtxt[subtxt.rfind("\n")+1:])
			else:
				x = x1 + font.measure(subtxt)
			return (x,y)
		else:
			return (None,None)

	# ----------------------------------------------------------------------
	# Determine drag location in canvas coordinates. event.x & event.y
	# don't seem to be what we want.
	# ----------------------------------------------------------------------
	def where(self, event):
		# where the corner of the canvas is relative to the screen:
		x_org = self.winfo_rootx()
		y_org = self.winfo_rooty()
		# where the pointer is relative to the canvas widget,
		# including scrolling
		x = self.canvasx(event.x_root - x_org)
		y = self.canvasy(event.y_root - y_org)
		return x, y

	# ----------------------------------------------------------------------
	# Tag selected card range
	# ----------------------------------------------------------------------
	def _tagRange(self, tag, card_from, card_to=-1):
		# Find bounding box of cards to be deleted
		fx1, fy1, fx2, fy2 = self.bbox(self.cardlist[card_from].box)
		if card_to < 0:
			tx2, ty2 = fx2, fy2
		else:
			tx1, ty1, tx2, ty2 = self.bbox(self.cardlist[card_to-1].box)

		# Mark cards to be deleted
		self.dtag(tag)
		self.addtag_overlapping(tag, fx1, fy1, tx2, ty2-1)
		return fx1, fy1, tx2, ty2

	# ----------------------------------------------------------------------
	# Focus
	# ----------------------------------------------------------------------
	def _focusOut(self, event=None):
		try:
			if self.focus_get() is None:
				return
		except: return

		self.commit()
		if not self.popupList._isPosted and self.active and self.edit:
			self.active.editItem(focus=False)
			self.focus("")

	# ----------------------------------------------------------------------
	# Update scrollbars
	# ----------------------------------------------------------------------
	def _updateScrollBars(self):
		"""Update scroll region for new size"""
		bb = self.bbox('all')
		if bb is None: return
		x1,y1,x2,y2 = bb
		self.configure(scrollregion=(x1+1,y1+1,x2,y2))

	# ----------------------------------------------------------------------
	# Find CardWidget from mouse position and activated if not active
	# ----------------------------------------------------------------------
	def findCardWidget(self, y, dy=1):
		try:
			y    = self.canvasy(y)
			item = self.find_overlapping(0, y-dy, _PAD+1, y+dy)
			cw   = self.boxes[item[0]]
			return cw
		except:
			return None

	# ----------------------------------------------------------------------
	def focusIn(self):
		self.page._focusIn()
		self.focus_set()

	# ----------------------------------------------------------------------
	# Mouse click
	# ----------------------------------------------------------------------
	def click(self, event):
		self.page._focusIn()
		CardWidgetBase.dialogShown = False
		self.dnd_time = None
		self.shift = (event.state & 1)
		self.focusIn()
		if self.shift: self.popupList.unPost()

		# Get coordinates
		x = int(self.canvasx(event.x))
		y = int(self.canvasy(event.y))
		self.mouseXY = (x,y)

		# Remember selected range (if any)
		card_from, card_to = self.getSelectRange()
		if card_from is None: return

		# Test if we clicked a CardWidget
		try:
			item = self.find_overlapping(0, y, _PAD+1, y+1)
			cw   = self.boxes[item[0]]

		except (KeyError, IndexError):
			# Unpost list if posted! but only when event is not handled
			# by the CanvasWidget.click() method. Otherwise the list
			# cannot be posted
			self.focus("")
			if self.active:
				self.commit()
				self.active.editItem()
			self.popupList.unPost()
			self.cardSelectClear()

		else:
			#if self.active is not cw:
			self.activate(cw)
			self.selAnchor = cw
			#else:
			item = cw.findItem(x,y)
			if item is None:
				cw.editItem()
			else:
				cw.editItem(item, x, y, False)

			if self.edit is not None: return
			if len(self.cardlist)<=1: return
			if self.selFrom==0 and self.selTo==len(self.cardlist):
				return

			if CardWidgetBase.dialogShown:
				CardWidgetBase.dialogShown = False
				return

			# Check if recent click change range
			cid = self.cardlist.index(self.active)
			if self.selFrom<0 and card_from>=0 \
					and cid>=card_from \
					and cid<card_to:
				self.dragFrom = card_from
				self.dragTo   = card_to
			else:
				self.dragFrom = self.selFrom
				self.dragTo   = self.selTo

			# Mark the time of the button press
			self.dnd_time  = event.time

	# ----------------------------------------------------------------------
	# Select a range of cards enclosed with a start_tag .. stop_tag
	# ----------------------------------------------------------------------
	def _selectCardRange(self, start_tag, stop_tag):
		cid   = self.cardlist.index(self.active)
		level = 1
		if self.active.card.tag == start_tag:
			fromid = cid
			toid   = fromid+1
			while toid<len(self.cardlist):
				cw = self.cardlist[toid]
				if cw.card is None:
					continue
				elif cw.card.tag == start_tag:
					level += 1
				elif cw.card.tag == stop_tag:
					level -= 1
					if level==0: break
				toid += 1
			self.cardSelect(fromid, toid)

		else:
			toid   = cid
			level  = 1
			fromid = toid-1
			while fromid>0:
				cw = self.cardlist[fromid]
				if cw.card is None:
					continue
				elif cw.card.tag == start_tag:
					level -= 1
					if level==0: break
				elif cw.card.tag == stop_tag:
					level += 1
				fromid -= 1
			self.cardSelect(fromid, toid)

	# ----------------------------------------------------------------------
	# Double click
	# ----------------------------------------------------------------------
	def double(self, event):
		if self.edit is not None:
			item   = self.focus()
			try: cursor = self.index(item, INSERT)
			except TclError: return

			begin, end = findWord(self.itemcget(item,"text"), cursor)
			if begin is None: return

			self.select_clear()
			self.select_from(item, begin)
			self.select_to(item, end)

		# Out in the blue...
		elif self.active is None:
			self.cardSelectClear()

		# On hidden cards
		elif self.active.card is None:
			self.expandCard()

		elif event.state&CONTROL_MASK:
			self.hideCard()

		# On preprocessor block #if...#endif

		# XXX if double click on preprocessor #if..#else..#endif
		# select all block
		elif self.active.card.tag in Input._PREPRO_BLOCK:
			cid = self.cardlist.index(self.active)
			indent = self.active._indent

			# Scan for the starting card with the same ident level
			fromid = cid
			while fromid>0:
				cw = self.cardlist[fromid]
				if cw.card is not None and \
				   cw.card.tag in ("#if", "#ifdef", "#ifndef") and \
				   cw._indent == indent: break
				fromid -= 1

			# Scan for the #endif card with the same indent level
			toid = cid
			while toid<len(self.cardlist):
				cw = self.cardlist[toid]
				if cw.card is not None and \
				   cw.card.tag == "#endif" and \
				   cw._indent == indent: break
				toid += 1
			self.cardSelect(fromid, toid)

		# if on #include..#endinclude
		elif self.active.card.tag in ("#include", "#endinclude"):
			self._selectCardRange("#include", "#endinclude")

		# start..end transform
		elif self.active.card.tag in ("$start_transform", "$end_transform"):
			self._selectCardRange("$start_transform", "$end_transform")

		# start..end translat
		elif self.active.card.tag in ("$start_translat", "$end_translat"):
			self._selectCardRange("$start_translat", "$end_translat")

		# start..end expansion
		elif self.active.card.tag in ("$start_expansion", "$end_expansion"):
			self._selectCardRange("$start_expansion", "$end_expansion")

		# Select range of similar cards
		else:
			self.cardSelectClear()
			if self.active is None: return
			if isinstance(self.active, HiddenCardWidget):
				pass
			else:
				pos = self.cardlist.index(self.active)
				tag = self.cardlist[pos].card.tag
				begin = pos
				while begin>=0:
					card = self.cardlist[begin].card
					if card and card.tag == tag:
						begin -= 1
					else:
						break
				end = pos
				while end < len(self.cardlist):
					card = self.cardlist[end].card
					if card and card.tag == tag:
						end += 1
					else:
						break
				self.cardSelect(begin+1, end-1)

	# ----------------------------------------------------------------------
	# Accept dnd messages, i.e. we're a legit drop target, and we do
	# implement d&d functions.
	# ----------------------------------------------------------------------
	def dnd_accept(self, source, event):
		return self

	# ----------------------------------------------------------------------
	# Create drag object
	# ----------------------------------------------------------------------
	def dnd_enter(self, source, event):
		self.focusIn()
		if source.active:
			if source.dragFrom >= 0:
				source.cardSelect(self.dragFrom, self.dragTo-1)
				if source.dragFrom != source.dragTo-1:
					source.dnd_label = "%s ... %s" % (
						source.cardlist[source.dragFrom].fromTag(),
						source.cardlist[source.dragTo-1].toTag())
				else:
					source.dnd_label = \
						source.cardlist[source.dragFrom].fromTag()
			elif isinstance(source.active, HiddenCardWidget):
				if source.active.fromCard is not source.active.toCard:
					source.dnd_label = "%s ... %s" % \
						(source.active.fromTag(), source.active.toTag())
				else:
					source.dnd_label = source.active.fromTag()
			else:
				source.dnd_label = source.active.fromTag()
			source.cut()
			source.dnd_active = source.active
			source.deactivate()

		self.dnd_item = self.create_text(event.x, event.y,
					text = source.dnd_label,
					font = _TAG_FONT,
					fill = _TAG_COLOR)
		self.dnd_motion(source, event)

	# ----------------------------------------------------------------------
	# Move drag icon
	# ----------------------------------------------------------------------
	def dnd_motion(self, source, event):
		if self.page.sbv.winfo_ismapped():
			# Scroll window if necessary
			# Vertical scroll
			h   = self.winfo_height()
			p10 = max(h//20, 20)
			if event.y <= p10 and self.canvasy(event.y) > p10:
				self.yview(SCROLL, -1, UNITS)
			elif event.y >= h - p10:
				self.yview(SCROLL, 1, UNITS)

			# Horizontal scroll
			w   = self.winfo_height()
			p10 = max(w//20, 20)
			if event.x <= p10:
				self.xview(SCROLL, -1, UNITS)
			elif event.x >= w - p10:
				self.xview(SCROLL, 1, UNITS)

		# Move widget
		x, y = self.where(event)
		try:
			x1, y1, x2, y2 = self.bbox(self.dnd_item)
		except TypeError:
			return
		xm = (x1+x2)//2		# From the middle
		ym = (y1+y2)//2
		self.move(self.dnd_item, x-xm, y-ym)

	# ----------------------------------------------------------------------
	def dnd_leave(self, source, event):
		self.delete(self.dnd_item)
		self.dnd_item = None

	# ----------------------------------------------------------------------
	# Notification that dnd processing has been ended. It DOES NOT imply
	# that we've been dropped somewhere useful, we could have just been
	# dropped into deep space and nothing happened to any data structures,
	# or it could have been just a plain mouse-click w/o any dragging.
	# ----------------------------------------------------------------------
	def dnd_end(self, target, event):
		if self.dnd_item:
			self.delete(self.dnd_item)
			self.dnd_item = None

		if target is None and self.dnd_label is not None:
			self.page.undo()	# Restore cards
			self.activate(self.dnd_active.card)

		self.motionBind()	# Restore motion binding
		self.dnd_label = None

	# ----------------------------------------------------------------------
	# Object has been dropped here
	# ----------------------------------------------------------------------
	def dnd_commit(self, source, event):
		if self.dnd_item:
			self.delete(self.dnd_item)
			self.dnd_item = None
		if source.dnd_label is None: return

		x, y = self.where(event)
		cw = self.findCardWidget(event.y, _PAD)
		if cw is None and len(self.cardlist)>0:
			x1, y1, x2, y2 = self.bbox(self.cardlist[0].box)
			if y < y1:
				cw = self.cardlist[0]
			else:
				cw = self.cardlist[-1]

		up = False
		if cw is not None:
			x1, y1, x2, y2 = self.bbox(cw.box)
			if y < (y1+y2)/2:
				cid = self.cardlist.index(cw)
				if cid>0:
					cw = self.cardlist[cid-1]
				else:
					up = True
			self.activate(cw)

		if source is not self and source.dnd_active is not None:
			source.activate(source.dnd_active)

		self.paste()
		self.clipboard_clear()
		if up: self.moveUp()

	# ----------------------------------------------------------------------
	# ADD/DEL CARDs
	# ----------------------------------------------------------------------
	# ----------------------------------------------------------------------
	# Add a card (Normal or Hidden) or range of cards [card..to_card] to InputCanvas
	# @param card	card to be added
	# @param to_card	range as hidden cards
	# @param pos	position where to add the card
	# @param first	when multiple cards are inserted (in the middle) mark what is
	#		is to be moved when first is True
	# @param last	when multiple cards are inserted (in the middle) perform
	#		the move when last is True
	# @return card widget
	# ----------------------------------------------------------------------
	def addCard(self, card, to_card=None, pos=None, first=True, last=True):
		if len(self.cardlist) == 0:
			pos  = 0
			ypos = _PAD
		elif pos is None or pos>=len(self.cardlist):
			# Add after the last one
			ypos = self.cardlist[-1].yend()
			pos  = len(self.cardlist)
		else:
			if pos<=0:
				pos  = 0
				ypos = _PAD
			else:
				ypos = self.cardlist[pos-1].yend()
			if first:
				# mark for move all cards that follow
				self.itemconfig("active", width=_NORMAL_BORDER)
				ax1, ay1, ax2, ay2 = self.bbox("all")
				self.dtag("move")
				self.addtag_overlapping("move", ax1, ypos, ax2, ay2)
				self._ypos = ypos

		# Create new card widget
		if isinstance(card, int):
			card = self.input.cardlist[card]

		if to_card is not None:
			if isinstance(to_card, int):
				to_card = self.input.cardlist[to_card]
			cw = HiddenCardWidget(self, card, to_card, ypos)
		else:
			cw = CardWidget(self, card, ypos)

		# Append card to list
		if pos>=len(self.cardlist):
			self.cardlist.append(cw)
		else:
			self.cardlist.insert(pos, cw)
			if last:
				by2 = self.bbox(cw.box)[3]
				self.move('move', 0, by2 - self._ypos + _PAD)
				self.itemconfig("active", width=_SEL_BORDER)

		# Activate card if necessary
		if self.active:
			if self.active.card:
				if self.active.card is card:
					self.activate(cw)
			elif self.active.fromCard is card:
				self.activate(cw)
		return cw

	# ----------------------------------------------------------------------
	# Add card but in correct position
	# ----------------------------------------------------------------------
	def addCardSorted(self, card):
		cardpos = card.pos()
		i = 0
		for i, cw in enumerate(self.cardlist):
			if cw.toPos() < cardpos:
				continue

			elif cardpos == cw.fromPos():
				# Same starting position
				if isinstance(cw,HiddenCardWidget) and cw.fromPos() == cw.toPos():
					# Is hidden with only one card
					# replace the hidden with a normal cardwidget
					self.delCard(i)
					return self.addCard(card,None,i)
				else:
					# Add a card at the start
					newcw = self.addCard(card,None,i)
					# and update the widget
					cw.fromCard = self.input.cardlist[cardpos+1]
					cw.refresh()
					return newcw

			elif cardpos == cw.toPos():
				# at the end of the range
				cw.toCard = self.input.cardlist[cardpos-1]
				cw.refresh()
				return self.addCard(card,None,i+1)

			elif cardpos < cw.fromPos():
				# has to be added before this HiddenCardWidget
				return self.addCard(card,None,i)

			else:
				# Insert in the HiddenCardWidget
				# by splitting into two
				tocard = cw.toCard
				cw.toCard = self.input.cardlist[cardpos-1]
				cw.refresh()
				newcw = self.addCard(card,None,i+1)
				self.addCard(self.input.cardlist[cardpos+1], tocard, i+2)
				return newcw
		return self.addCard(card,None,i)

	# ----------------------------------------------------------------------
	# Delete a card(s) widgets from the canvas box
	# ----------------------------------------------------------------------
	def delCard(self, card_from, card_to=-1):
		if len(self.cardlist)==0:
			return undo.NullUndo
		elif card_to < 0:
			card_to = card_from+1
		elif card_to > len(self.cardlist):
			card_to = len(self.cardlist)

		if self.active is not None:
			actid = self.cardlist.index(self.active)
		else:
			actid = -1

		# Find bounding box of cards to be deleted
		self.itemconfig("active", width=_NORMAL_BORDER)
		x1,y1,x2,y2 = self._tagRange("del", card_from, card_to)

		# activate next card
		if card_from <= actid < card_to:
			if card_to < len(self.cardlist):
				self.activate(self.cardlist[card_to])
			elif card_from>0:
				self.activate(self.cardlist[card_from-1])
			else:
				# Delete all cards
				self.active = None
				self.page.updateLabel()

		# delete the widgets
		self.delete("del")

		if card_to < len(self.cardlist):
			ax1, ay1, ax2, ay2 = self.bbox("all")
			self.dtag("move")
			self.addtag_overlapping("move", ax1, y2, ax2, ay2)
			self.move('move', 0, y1 - y2 - _PAD)

		# Delete cards and prepare undo list
		for cw in self.cardlist[card_from:card_to]:
			self.delBox(cw)
			if isinstance(cw, HiddenCardWidget):
				to_card = cw.toPos()
			else:
				to_card = None

		del self.cardlist[card_from:card_to]
		self.selFrom = -1

	# ----------------------------------------------------------------------
	# Delete all cards widgets
	# ----------------------------------------------------------------------
	def deleteAll(self):
		#undoinfo = []
		for cw in self.cardlist:
			if isinstance(cw, HiddenCardWidget):
				to_card = cw.toPos()
			else:
				to_card = None
			#undoinfo.append((InputCanvas.addCard, self,
			#		cw.fromPos(), to_card, 0))

		self.delete('all')
		del self.cardlist[:]
		del self.boxes
		self.boxes = {}
		#return undoinfo

	# ----------------------------------------------------------------------
	# Add/Del Box
	# Speedup access to CardWidget from boxes in click()
	# ----------------------------------------------------------------------
	def addBox(self, card_widget):
		self.boxes[card_widget.box] = card_widget

	# ----------------------------------------------------------------------
	def delBox(self, card_widget):
		del self.boxes[card_widget.box]

	# ----------------------------------------------------------------------
	# Create the CardWidgets for the cardlist provided
	# ----------------------------------------------------------------------
	def displayCards(self, cardlist):
		self.deleteAll()

		# Remember active
		if self.active is not None:
			actpos = self.active.fromPos()
		else:
			actpos = 0

		self.active = None

		# Initialize metrics for faster drawing
		CardWidgetBase.initMetrics(self)

		# insert to canvas
		if len(self.input.cardlist) > 0:
			fromCard = self.input.cardlist[0]
			startTime = self.resetInterrupt()
			for card in cardlist:
				pos = card.pos()
				if fromCard and fromCard.pos() != pos:
					toCard = self.input.cardlist[pos-1]
					self.addCard(fromCard, toCard)
				self.addCard(card)
				if pos+1 < len(self.input.cardlist):
					fromCard = self.input.cardlist[pos+1]
				else:
					fromCard = None
				if time.time()-startTime > _CARDS_TIMEOUT:
					break
			if fromCard:
				self.addCard(fromCard, self.input.cardlist[-1])
			del cardlist

		# Try to find which card to activate
		for cw in self.cardlist:
			if cw.fromPos() <= actpos <= cw.toPos():
				self.activate(cw)
				break

		if self.active is None and self.cardlist:
			self.activate(self.cardlist[0])

		self._updateScrollBars()
		self.page.refreshStatus()

		if self.active is not None:
			self.active.see()
		#print ">>>>>>>>>>>>>>>> displayCards t=",time.time()-startTime

	# ----------------------------------------------------------------------
	# Check card widget list and input card list for non-existing cards
	# ----------------------------------------------------------------------
	def insertAdditionalCards(self):
		lastpos       = len(self.input.cardlist)
		startTime     = self.resetInterrupt()
		longNotify    = False
		updateCounter = 0
		for i in range(len(self.cardlist)-1,0,-1):	# In reverse order
			cw = self.cardlist[i]
			if cw.card is not None:
				if cw.card.input is None:	# Card deleted
					self.delCard(i)
				elif cw.card.pos() < 0:
					# delete card
					self.delCard(i)
				elif cw.card.pos() < lastpos-1:
					# insert additional cards
					for card in self.input.cardlist[cw.card.pos()+1 : lastpos]:
						i += 1
						# timeout add them as HiddenCardWidget
						if time.time()-startTime > _CARDS_TIMEOUT:
							if not self.flair.notifyIsVisible():
								self.flair.notify("Too long",
									"Process takes too long",
									tkFlair.NOTIFY_NORMAL,
									"Interrupt?", self.interrupt)
								longNotify = True
							updateCounter -= 1
							if updateCounter<=0:
								updateCounter = 100
								self.flair.update()

							if self._interrupt:
								self.addCard(card, self.input.cardlist[lastpos],i)
								startTime = self.resetInterrupt()
								break
						else:
							self.addCard(card, None, i)
					lastpos = cw.card.pos()
				else:
					lastpos = cw.card.pos()
			else:
				if cw.fromPos() < 0 or cw.toPos() < 0:
					self.delCard(i)
				elif cw.toPos() < lastpos-1:
					# insert additional cards
					for card in self.input.cardlist[cw.toPos()+1 : lastpos]:
						i += 1
						# timeout add them as HiddenCardWidget
						if time.time()-startTime > _CARDS_TIMEOUT:
							if not self.flair.notifyIsVisible():
								self.flair.notify("Too long",
									"Process takes too long",
									tkFlair.NOTIFY_NORMAL,
									"Interrupt?", self.interrupt)
								longNotify = True
							updateCounter -= 1
							if updateCounter<=0:
								updateCounter = 100
								self.flair.update()
							if self._interrupt:
								self.addCard(card, self.input.cardlist[lastpos],i)
								startTime = self.resetInterrupt()
								break
						self.addCard(card, None, i)
					lastpos = cw.fromPos()
				else:
					lastpos = cw.fromPos()
		if longNotify: self.flair.notifyHide()

	# ----------------------------------------------------------------------
	# Refresh canvas, refresh only cardwidgets with new values
	# ----------------------------------------------------------------------
	def refresh(self):
		# XXX not sure if we need that!
		self.insertAdditionalCards()

		for c in self.cardlist:
			if c.isModified():
				c.refresh()
		self._updateScrollBars()

		if self.selFrom>=0:
			self.cardSelect(self.selFrom, self.selTo-1)

		self.page.refreshStatus()
		self.input.checkNumbering()	# XXX

	# ----------------------------------------------------------------------
	# Force redrawing of the canvas
	# ----------------------------------------------------------------------
	def redraw(self):
		# FIXME Could well be adding only to the self.cardlist!
		self.insertAdditionalCards()

		self.boxes.clear()
		self.delete('all')
		cards = self.cardlist	# Save card list
		self.cardlist = []
		self.edit = None
		for c in cards:
			if c.card:
				self.addCard(c.card)
			else:
				self.addCard(c.fromCard, c.toCard)
		self._updateScrollBars()

		if self.selFrom>=0:
			self.cardSelect(self.selFrom, self.selTo-1)

		self.page.refreshStatus()
		self.input.checkNumbering()	# XXX

	# ----------------------------------------------------------------------
	# Activate the CardWidget containing the card
	# ----------------------------------------------------------------------
	def activate(self, card):
		# If normal Input card then search for card widget
		if isinstance(card, Input.Card):
			for cw in self.cardlist:
				if cw.card is card:
					break
			else:
				cw = self.addCardSorted(card)
				self._updateScrollBars()
				#raise Exception("Card to be activated not found")
			card = cw

		# If card is not found do nothing
		if card is None: return

		# Card Widget
		if self.active is not None:
			if self.edit is not None:
				self.commit()

			if card is not self.active:
				card.updateState(True)

			if self.active is card:
				return undo.NullUndo

			#undoinfo = (self.activate, card)
		#else:
			#undoinfo = undo.NullUndo

		self.highlight("")
		if self.active is not None:
			self.itemconfig(self.active.box, width=_NORMAL_BORDER)

			try:
				cid = self.cardlist.index(card)
			except ValueError:
				cid = -1

			if self.shift:
				if self.selAnchor is None:
					self.selAnchor = self.active
				# To avoid problem when cardSelectClear is called
				self.active = card
				self.cardSelect(self.selAnchor, card)
			else:
				if self.selFrom>=0:
					if cid<self.selFrom or cid>=self.selTo:
						self.selAnchor = None
						self.cardSelectClear()
						if self.active.card and self.active.card.invalid:
							self.itemconfig(self.active.box,
								fill=_INVALID_COLOR)
						else:
							self.itemconfig(self.active.box,
								fill=_NORMAL_COLOR)
				else:
					self.selAnchor = None
					self.cardSelectClear()
					if self.active.card and self.active.card.invalid:
						self.itemconfig(self.active.box,
							fill=_INVALID_COLOR)
					else:
						self.itemconfig(self.active.box,
							fill=_NORMAL_COLOR)
		self.active = card

		# FIXME adjustTextBox
		self.itemconfig(card.box, fill=_ACTIVE_COLOR, width=_SEL_BORDER)

		self.dtag("active","active")
		self.addtag_withtag("active", card.box)

		card.updateLabel()
		# Avoid unnecessary scrolling for undo/redo on the second canvas
		# If canvas doesn't have the focus do not scroll to see it
		try:
			if self is self.focus_get():
				card.see()
		except KeyError:
			# In case of "tearoff" for a submenu!
			pass
		#self._updateScrollBars()
		self.page.refreshStatus()

		#return undoinfo

	# ----------------------------------------------------------------------
	def activateIfNeeded(self, card):
		if self.active is None or self.active.card is not card:
			self.commit()
			self.activate(card)
			self.active.see()

		if self.active is not None and self.active.card is not None :
			self.page.updateLabel(self.active.card)
			layout = Layout.getLayout(card)
			layoutChanged = self.active.changeLayout(layout)
			if not layoutChanged: self.active.refresh()
			self.highlight(self._highlightItem)

	# ----------------------------------------------------------------------
	# Activate next card
	# ----------------------------------------------------------------------
	def nextActive(self):
		if len(self.cardlist) == 0:
			return

		# activate next widget
		try:
			a = self.cardlist.index(self.active)
		except ValueError:
			self.activate(self.cardlist[0])
		else:
			a += 1
			if a >= len(self.cardlist): return
			self.activate(self.cardlist[a])

	# ----------------------------------------------------------------------
	# Activate previous card
	# ----------------------------------------------------------------------
	def prevActive(self):
		if len(self.cardlist) == 0:
			return

		# activate previous widget
		try:
			a = self.cardlist.index(self.active)
		except:
			self.activate(self.cardlist[0])
		else:
			a -= 1
			if a < 0: return
			self.activate(self.cardlist[a])

	# ----------------------------------------------------------------------
	# Deactivate card
	# ----------------------------------------------------------------------
	def deactivate(self):
		if self.active is None: return
		if self.active.card and self.active.card.invalid:
			self.itemconfig("active",
				fill=_INVALID_COLOR,
				width=_NORMAL_BORDER)
		else:
			self.itemconfig("active",
				fill=_NORMAL_COLOR,
				width=_NORMAL_BORDER)
		self.dtag("active","active")
		self.commit()
		self.active.editItem()
		self.active = None

	# ----------------------------------------------------------------------
	# Commit changes to active card and to selected range with similar cards
	# ----------------------------------------------------------------------
	def commit(self):
		if self.active is None: return

		# Remember active (to be restored in the end)
		active = self.active

		try:
			if not active.commit(): return
		except:
			# Recover from possibly fatal errors
			self.focus("")
			self.edit   = None
			active = None
			raise

		if  self.selFrom < 0 \
		   or active.card is None \
		   or self.edit is None:
			#self.page._refreshOtherFrames()
			return

		if self.edit < 0:	# Comment change
			val = active.card.comment()
		else:
			try: lbl, dat = active.layout[self.edit]
			except: lbl, dat, desc = active.layout[self.edit]
			# Get the value!
			if isinstance(dat, tuple): #List
				dat = dat[0]
			if callable(dat):
				val = dat(self.layout, active.card)
			else:
				val = active.card.what(dat)

		card_from, card_to = self.getSelectRange()
		if card_from is None: return

		# Skip expression change in regions
		if active.card.tag == "REGION" and dat==-1:
			#self.page._refreshOtherFrames()
			return

		undoinfo = []
		for cw in self.cardlist[card_from:card_to]:
			# Check for the same card type
			card = cw.card
			if card is None: continue

			self.active = cw	# fool the system

			if self.edit < 0:
				undoinfo.append(self.flair.setCommentUndo(card, val))
				cw.refresh()
				continue

			if card.tag != active.card.tag: continue
			# Modify the value
			if callable(dat):
				dat(self.layout, card, val)
				cw.refresh()
			else:
				undoinfo.append(self.flair.setWhatUndo(card, dat, val))
				cw.refresh()
		if undoinfo:
			self.addUndo(undo.createListUndo(undoinfo,
						undoinfo[0][0]))
		self.active = active
		#self.master._refreshOtherFrames()

	# ----------------------------------------------------------------------
	# CARD SELECTION
	# ----------------------------------------------------------------------
	# Return selected range (from, to+1)
	# ----------------------------------------------------------------------
	def getSelectRange(self):
		if self.active is None: return (None, None)
		if self.selFrom>=0:
			return (min(self.selFrom,len(self.cardlist)), \
			        min(self.selTo,  len(self.cardlist)))
		else:
			try:
				cid = self.cardlist.index(self.active)
			except ValueError:
				return 0,1
			return cid, cid+1

	# ----------------------------------------------------------------------
	# Clear Card selection
	# ----------------------------------------------------------------------
	def cardSelectClear(self):
		if self.selFrom >= 0:
			# Remove colors
			for cw in self.cardlist[self.selFrom:self.selTo]:
				if cw is not self.active:
					if cw.card and cw.card.invalid:
						self.itemconfig(cw.box,
							fill=_INVALID_COLOR)
					else:
						self.itemconfig(cw.box,
							fill=_NORMAL_COLOR)
		self.selFrom = -1

	# ----------------------------------------------------------------------
	# Select range
	# ----------------------------------------------------------------------
	def cardSelect(self, from_card, to_card):
		self.cardSelectClear()
		if isinstance(from_card, int):
			f = from_card
		else:
			try: f = self.cardlist.index(from_card)
			except: return

		if isinstance(to_card, int):
			t = to_card
		elif to_card==END:
			t = len(self.cardlist)-1
		else:
			t = self.cardlist.index(to_card)

		# swap values
		if f>t: f,t = t,f

		# Remove colors
		t += 1
		for cw in self.cardlist[f:t]:
			self.itemconfig(cw.box,fill=_ACTIVE_COLOR)

		self.selFrom = f
		self.selTo   = t

	# ----------------------------------------------------------------------
	# CARD/TEXT
	#
	# Select All, text or cards
	# ----------------------------------------------------------------------
	def selectAll(self, event=None):
		item = self.focus()
		if item and self.edit is not None:
			self.select_from(item, 0)
			self.select_to(item, END)
		else:
			self.cardSelect(0, END)

	# ----------------------------------------------------------------------
	def unSelectAll(self, event=None):
		#item = self.focus()
		#if item and self.edit is not None:
		#	self.select_from(item, 0)
		#	self.select_to(item, END)
		#else:
		self.cardSelectClear()

	# ----------------------------------------------------------------------
	# delete selection, text or cards
	# ----------------------------------------------------------------------
	def selectClear(self, event=None):
		self.select_clear()
		self.selEditAnchor = -1

	# ----------------------------------------------------------------------
	# Copy selected range to clipboard
	# ----------------------------------------------------------------------
	def _copy2Clipboard(self):
		hidden = 0		# Hidden group
		card_from, card_to = self.getSelectRange()
		if card_from is None: return

		# Copy also as pickle format
#		sio = StringIO()
		sio = io.BytesIO()
#		sio.write(tkFlair._FLAIRCLIP)
#		sio.write("<card>")
		pickler = pickle.Pickler(sio)
		pickler.dump(tkFlair._FLAIRCLIP)
		pickler.dump("<card>")

		for cw in self.cardlist[card_from : card_to]:
			if cw.card is None:
				hidden += 1
				for c in cw.cards():
					c["@hidden"] = hidden
					c.dump(pickler)
			else:
				cw.card["@hidden"] = 0
				cw.card.dump(pickler)
		self.clipboard_clear()
#		self.clipboard_append(sio.getvalue())
		self.clipboard_append(binascii.b2a_hex(sio.getvalue()), type='STRING')

	# ----------------------------------------------------------------------
	# Insert clipboard at active position
	# ----------------------------------------------------------------------
	def _pasteBuffer(self, cards):
		pos   = 0
		cwpos = 0
		if self.active is not None:
			if isinstance(self.active, HiddenCardWidget):
				pos = self.active.toPos()+1
			else:
				pos = self.active.fromPos()+1
			cwpos = self.cardlist.index(self.active)+1

		frompos   = pos	# Remember position
		fromcwpos = cwpos
		hidden    = 0
		fromcard  = None
		fromfirst = True
		newcard   = None
		undoinfo  = [(self.flair.refreshUndo,), self.flair.renumberUndo()]
		self.layout.cache(None)	# Clean caching of lists

		renamed = ""
		for card in cards:
			h = card.get("@hidden",0)

			tag = card.tag
			if card.type() == Input.Card.BODY:
				lst = self.layout.bodyList(False)
			elif tag == "#define":
				lst = self.layout.defineList()
			elif tag == "REGION":
				lst = self.layout.regionList()
			elif tag == "MATERIAL":
				lst = self.layout.materialList()
			elif tag in ("EVENTBIN", "RESNUCLE", "USRBDX", "USRBIN",
				     "USRCOLL", "USRTRACK", "USRYIELD"):
				lst = self.layout.binningList()
			#elif tag == "ROT-DEFI":
			#	lst = self.layout.rotdefiList()
			else:
				lst = None

			lastcard = newcard
			newcard  = card.clone()

			# Check for duplicate
			if card.sdum():
				if lst is not None and card.sdum() in lst:
					# Find unique name
					# find number in the end of the name
					pat = Project._LAST_NUM.match(card.sdum())
					if pat:
						name = pat.group(1)
						n = int(pat.group(2))+1
					else:
						name = card.sdum()
						n = 1
					while True:
						sn = str(n)
						guess = "%s%s"%(name[:8-len(sn)],sn)
						if guess not in lst:
							renamed += "   %s:\t%s\t-> %s\n"%(tag, name, guess)
							newcard.setName(guess)
							break
						n += 1

			# Add to undo
			undoinfo.append(self.flair.addCardUndo(newcard, pos, False))
			pos += 1
			if h != 0:
				if h != hidden:
					if fromcard:
						# Create widget
						self.addCard(fromcard, lastcard, cwpos, fromfirst, False)
						cwpos += 1
					hidden = h
					fromcard = newcard
					fromfirst = card is cards[0]
			else:
				if fromcard:
					# Create widget
					self.addCard(fromcard, lastcard, cwpos, fromfirst, False)
					cwpos += 1
					fromcard = None
					fromfirst = card is cards[0]
				self.addCard(newcard, None, cwpos, card is cards[0], card is cards[-1])
				cwpos += 1
		if fromcard:
			self.addCard(fromcard, newcard, cwpos, fromfirst, True)

		undoinfo.append(self.flair.renumberUndo())
		undoinfo.append((self.flair.refreshUndo,))
		self.addUndo(undo.createListUndo(undoinfo, "Paste card(s)"))

		if renamed:
			self.flair.notify("Cards renamed",
				"The following cards have been renamed:\n%s"%(renamed),
				tkFlair.NOTIFY_WARNING)

		self._updateScrollBars()
		self.layout.cache(False)	# Disable caching of lists
		self.input.checkNumbering()	# XXX
		self.input.clearCache()
		self.page._refreshOtherFrames()
		return fromcwpos, cwpos

	# ----------------------------------------------------------------------
	def addUndo(self, undoinfo):
		self.flair.addUndo(undoinfo)

	# ======================================================================
	# Copy selected content to clipboard
	# ======================================================================
	def copy(self, event=None):
		item = self.focus()
		if self.edit is not None:
			if not self.has_selection(): return
			text = self.itemcget(item,"text")
			sel_text = text[self.index(item, SEL_FIRST):
					self.index(item, SEL_LAST)+1]
			self.clipboard_clear()
			self.clipboard_append(sel_text.encode("ASCII"))

		elif self.active:
			self._copy2Clipboard()

	# ----------------------------------------------------------------------
	# Cut selected content to clipboard
	# ----------------------------------------------------------------------
	def cut(self, event=None):
		item = self.focus()
		if self.edit is not None:
			if not self.has_selection(): return
			self.copy()

			if self.multiline: self.active.adjustTextBox(False)
			self.dchars(item, SEL_FIRST, SEL_LAST)
			self.selectClear()
			if self.multiline: self.active.adjustTextBox(True)
			self.highlight(item)

		elif self.active:
			self._copy2Clipboard()
			t0 = time.time()
			self.delSelectedCards()
			print("Del time",time.time()-t0)

	# ----------------------------------------------------------------------
	# Paste clipboard content
	# ----------------------------------------------------------------------
	def paste(self, event=None):
		try: clipboard = self.selection_get(selection='CLIPBOARD')
		except: return

		if self.edit is not None:
			try: clipboard = self.selection_get(selection='CLIPBOARD')
			except: return
			self.pasteText(clipboard)
		else:
			cards = self.flair.clipboard2Cards()
			if not cards: return
			self.cardSelectClear()
			card_from, card_to = self._pasteBuffer(cards)
			self.activate(self.cardlist[card_from])
			if card_to > card_from+1:
				self.cardSelect(card_from, card_to-1)
			self.focusIn()
			self.page.refreshStatus()

	# ----------------------------------------------------------------
	# Paste text
	# ----------------------------------------------------------------
	def pasteText(self, txt):
		item = self.focus()
		if self.multiline: self.active.adjustTextBox(False)
		if self.has_selection():
			self.dchars(item, SEL_FIRST, SEL_LAST)
			self.selectClear()
		if not self.multiline:
			if txt.find("\n")>=0 or len(txt)>80:
				self.flair.notify("Cannot paste",
					"String too long to paste",
					tkFlair.NOTIFY_ERROR)
				return
		self.insert(item, INSERT, txt)
		if self.multiline: self.active.adjustTextBox(True)
		self.highlight(item)
		self.active.see(item)
		self._updateScrollBars()

	# ----------------------------------------------------------------------
	def clone(self, event=None):
		# During editing it pastes the old text....
		if self.active is not None:
			self.commit()
			self.active.editItem()
			# Find last item to active
			# Remember selected range (if any)
			card_from, card_to = self.getSelectRange()
			if card_from is None: return
			self.activate(self.cardlist[card_to-1])
			self.cardSelect(card_from, card_to-1)
			self.copy()
			self.paste()

	# ----------------------------------------------------------------
	# POPUP LIST
	# ----------------------------------------------------------------
	# listbox handling methods
	# ----------------------------------------------------------------------
	def postList(self, lst=[], sel=None):
		self.hideBalloon()

		# Find geometry
		x1,y1,x2,y2 = self.bbox(self.active.box)
		w1 = (3*(x2-x1))/16
		x1,y1,x2,y2 = self.bbox(self.active.cell[self.edit][1])

		wx = x1 - int(self.canvasx(0))
		x  = wx + self.winfo_rootx()
		w  = max(w1, x2-x1)
		if w>w1+5: w += (w1 - w%w1)

		wy = y2-int(self.canvasy(0))
		y = wy + self.winfo_rooty()
		h = self.popupList.height()
		y = min(y, self.winfo_screenheight()-h)

		#self.popup.configure(width=x2-x1)
		self.popupList._commit = True
		#self.popupList._selectKeys = ""
		self.popupList.post("%dx%d+%d+%d" % (w,h,x,y), lst, sel)

	# ----------------------------------------------------------------------
	# Display a popup list with current bodies
	# @see CardWidget.insertSelection()
	# ----------------------------------------------------------------------
	def insertBody(self):
		x,y = self.cursorPosition()
		x  += self.winfo_rootx() - self.canvasx(0)
		y  += self.winfo_rooty() - self.canvasy(0)
		w   = self.winfo_width()//3
		h   = self.popupList.height()

		if self.has_selection():
			self.dchars(self.focus(), SEL_FIRST, SEL_LAST)
		self.popupList._commit = False
		#self.popupList._selectKeys = "-+|# "
		self.popupList.post("%dx%d+%d+%d"%(w,h,x,y),
				self.layout.bodyList(), "")

	# ----------------------------------------------------------------
	# TEXT EDITING
	#
	# text editing methods
	# ----------------------------------------------------------------------
	def highlight(self, item=None):
		# mark focused item.  note that this code recreates the
		# rectangle for each update, but that's fast enough for
		# this case.
		self.delete("highlight")
		if item is None:
			item = self._highlightItem
		else:
			self._highlightItem = item
		try:
			bbox = self.bbox(item)
		except TclError:
			return
		if bbox:
			i = self.create_rectangle(bbox, tag="highlight",
				outline="Gray")
			self.lower(i, item)

	# ----------------------------------------------------------------------
	def has_selection(self):
		# hack to work around bug in Tkinter 1.101 (Python 1.5.1)
		return self.tk.call(self._w, 'select', 'item')

	# ----------------------------------------------------------------------
	# Popup balloon help, binded to motion
	# ----------------------------------------------------------------------
	def motion(self, event):
		tkExtra.Balloon.hide()
		if self.motionId>=0:
			if event.time - self.motionTime < self.delay:
				self.after_cancel(self.motionId)
				self.motionId = -1
		if self.motionId<0:
			self.motionX = event.x
			self.motionY = event.y
			self.motionId = -1
#			self.motionId = self.after(self.delay, self.showBalloon)
			self.after(self.delay, self.showBalloon)
		self.motionTime = event.time

	# ----------------------------------------------------------------------
	def showBalloon(self):
		self.motionId = -1
		self.delay = 500	# reduce delay to 500 ms for next popup
		# Find the CardWidget
		x = int(self.canvasx(self.motionX))
		y = int(self.canvasy(self.motionY))
		try:
			item = self.find_overlapping(0, y, _PAD+1, y+1)
			cw   = self.boxes[item[0]]
		except (KeyError, IndexError):
			return

		# Prepare the text
		self._help = "Sorry: no info available..."	# default message

		if cw.card is None:	# Hidden Cards
			self._help = "%d hidden cards" \
				% (cw.toCard.pos() - cw.fromCard.pos()+1)
		else:
			item = cw.findItem(x,y)
			if item is None:
				self._help = "Card: %s\n%s" \
					% (cw.card.tag, cw.card.info.meaning)
			elif item == -1:
				self._help = "Comment text"
			else:
				try: lbl,dat,desc = cw.layout[item]
				except:
					lbl,dat = cw.layout[item]
					desc = None

				# We know "what"
				if isinstance(dat, tuple): dat = dat[0]
				if isinstance(dat,int):
					if dat==0:
						what = "SDUM"
					else:
						what = "WHAT(%d)"%(dat)
				else:
					what = None

				# if description is present
				if desc is not None:
					if what is None:
						self._help = desc
					else:
						self._help = "%s: %s"%(what,desc)
				else:
					if what is not None:
						case = cw.card.case()
						if case is not None:
							try:
								self._help = "%s: %s\nDefault: %s" \
									%(what,
									cw.card.info.extra[case][dat],
									cw.card.info.default[case][dat])
							except (IndexError, TypeError):
								self._help = what
						else:
							self._help = what

				# REGION expression, check for bodies
				if cw.card.tag == "REGION" and dat==-1:
					# find index position
					try:
						lbl,widget = cw.cell[item]
						cursor = self.index(widget, "@%d,%d"%(x,y))
						expr = self.itemcget(widget,"text")
						begin, end = findWord(expr, cursor)

						if begin is not None:
							body = expr[begin+1:end+1]
							# Find body card
							found = False
							if body=="VOXEL" and "VOXELS" in self.input.cards:
								found = True
								card = self.input.cards["VOXELS"][0]
							else:
								for tag in Input.BODY_TAGS:
									try:
										for card in self.input.cards[tag]:
											if card.sdum() == body:
												found = True
												break
									except KeyError:
										pass
									if found: break
							if found:
								self._help = str(card)
					except TclError:
						pass
		# Show the balloon
		tkExtra.Balloon.setWidget(self,
			self.motionX + self.winfo_rootx(),
			self.motionY + self.winfo_rooty())
		tkExtra.Balloon.show()

	# ----------------------------------------------------------------------
	def hideBalloon(self, event=None):
		if self.motionId>=0:
			self.after_cancel(self.motionId)
			self.motionId = -1
		else:
			tkExtra.Balloon.hide()

	# ----------------------------------------------------------------------
	# Select using the mouse
	# ----------------------------------------------------------------------
	def mouseSelect(self, event=None):
		item = self.focus()
		if item and self.edit is not None and not self.popupList._isPosted:
			# translate to the canvas coordinate system
			x = self.canvasx(event.x)
			y = self.canvasy(event.y)
			try:
				pos = self.index(item, "@%d,%d" % (x,y))
			except:
				return
			if self.selEditAnchor<0:
				self.selEditAnchor = self.index(item, "@%d,%d"
						% self.mouseXY)
				self.select_from(item, self.selEditAnchor)
				self.select_to(item, pos)
			elif pos == self.selEditAnchor:
				self.select_clear()
				self.selEditAnchor = -1
			else:
				self.select_to(item, pos)

		elif self.dnd_time:
			dt  = event.time - self.dnd_time
			self.dnd_time = None
			if dt<30: return
			event.num = 1	# Fake B1 press
			if Dnd.dnd_start(self, event):
				self.popupList.unPost()
				# Remove conflict with Drag Motion bind
				# we will fix it at the end of the dnd event
				self.unbind("<B1-Motion>")
				self.dnd_label = None

	# ----------------------------------------------------------------------
	def scanRelease(self, event):
		self.scandrag = False
		self.config(cursor="")

	# ----------------------------------------------------------------------
	def scanDrag(self, event):
		if self.scandrag:
			self.scan_dragto(event.x, event.y, 3)
		else:
			self.config(cursor="hand2")
			self.scan_mark(event.x, event.y)
			self.scandrag = True

	# ----------------------------------------------------------------------
	def tabKey(self, event):
		item = self.focus()
		if item:
			self.commit()
			self.active.editNext()
			return "break"

	# ----------------------------------------------------------------------
	def tabShiftKey(self, event):
		item = self.focus()
		if item:
			self.commit()
			self.active.editPrev()
			return "break"

	# ----------------------------------------------------------------------
	# Move cursor and select text when moving to left, home or up
	# ----------------------------------------------------------------------
	def editLeft(self, item, curr, pos):
		if self.shift:
			if self.selEditAnchor<0:
				self.selEditAnchor = curr
				self.select_from(item, curr)
				self.select_to(item, pos)
			elif pos == self.selEditAnchor:
				self.select_clear()
				self.selEditAnchor = -1
			elif pos > self.selEditAnchor:
				self.select_to(item, pos-1)
			else:
				self.select_to(item, pos)
		else:
			self.select_clear()
			self.selEditAnchor = -1

	# ----------------------------------------------------------------------
	# Move cursor and select text when moving to right, end or down
	# ----------------------------------------------------------------------
	def editRight(self, item, curr, pos):
		if self.shift:
			if self.selEditAnchor<0:
				self.selEditAnchor = curr
				self.select_from(item, curr)
				self.select_to(item, pos-1)
			elif pos == self.selEditAnchor:
				self.select_clear()
				self.selEditAnchor = -1
			elif pos < self.selEditAnchor:
				self.select_to(item, pos)
			else:
				self.select_to(item, pos-1)
		else:
			self.select_clear()
			self.selEditAnchor = -1

	# ----------------------------------------------------------------------
	def editUp(self):
		item  = self.focus()
		x, y = self.cursorPosition()
		bbox = self.bbox(item)
		y = bbox[1] - 12	# XXX XXX XXX should be properly calculated
		bbox = self.bbox(self.active.box)

		self.commit()
		self.active.editItem()

		if y < bbox[1]:
			self.prevActive()
			y = self.bbox(self.active.box)[3] - 12
		item = self.active.findItem(x,y)
		if item==0 and self.active.cell[0][1] is None:	# Find next item to highlight if any
			inext = item + 1
			while inext<len(self.active.cell) and self.active.cell[inext][0] is None:
				inext += 1
			if inext<len(self.active.cell): item = inext
		self.active.editItem(item, x, y, False)

	# ----------------------------------------------------------------------
	def editDown(self):
		item  = self.focus()
		x, y = self.cursorPosition()
		bbox = self.bbox(item)
		y = bbox[3] + 12	# XXX XXX XXX should be properly calculated
		bbox = self.bbox(self.active.box)

		self.commit()
		self.active.editItem()

		if y > bbox[3]:
			self.nextActive()
			y = self.bbox(self.active.box)[1] + 12
		item = self.active.findItem(x,y)
		if item==0 and self.active.cell[0][1] is None:	# Find next item to highlight if any
			inext = item + 1
			while inext<len(self.active.cell) and self.active.cell[inext][0] is None:
				inext += 1
			if inext<len(self.active.cell): item = inext
		self.active.editItem(item, x, y, False)

	# ----------------------------------------------------------------------
	# Find highlighted item and edit-it
	# ----------------------------------------------------------------------
	def editHighlighted(self):
		if self.active is None: return

		if self.active.card is None:
			self.expandCard()
			return

		# check for highlight
		bbox = self.bbox("highlight")
		self.commit()
		if bbox is not None:
			x1, y1, x2, y2 = bbox
			hitems = self.find_overlapping(x1+1,y1+1,x2-1,y2-1)
			if self.active.comment in hitems:
				self.active.editItem(-1)
				return
			else:
				i = 0
				editnext = False
				for (l,d) in self.active.cell:
					if d is not None and editnext:
						self.active.editItem(i)
						return
					if l in hitems:
						if d is not None:
							self.active.editItem(i)
							return
						else:
							editnext = True
					if d in hitems:
						self.active.editItem(i)
						return
					i += 1
		self.active.editNext()

	# ----------------------------------------------------------------------
	# Keyboard event handler
	# ----------------------------------------------------------------------
	def _handleKey(self, event):
		# widget-wide key dispatcher
		if event.state & _ALT_KEY == _ALT_KEY: return	# Ignore ALT key
		self.shift = (event.state & 1)
		self.hideBalloon()
		self.delay = tkExtra.Balloon.delay	# reset delay when key pressed

		item   = self.focus()
		active = self.active

		if not item:
			# operate on canvas
			if event.keysym == "Prior":
				self.yview(SCROLL, -1, PAGES)
				cw = self.findCardWidget(self.winfo_height()-_PAD)
				if cw and active is not cw:
					self.activate(cw)
			elif event.keysym == "Next":
				self.yview(SCROLL,  1, PAGES)
				cw = self.findCardWidget(_PAD)
				if cw and active is not cw:
					self.activate(cw)
			elif event.keysym == "Home":
				if len(self.cardlist)>0:
					self.activate(self.cardlist[0])
			elif event.keysym == "End":
				if len(self.cardlist)>0:
					self.activate(self.cardlist[-1])
			elif event.keysym == "Up":
				self.prevActive()
			elif event.keysym == "Down":
				self.nextActive()
			elif event.keysym == "Left":
				self.xview(SCROLL, -1, UNITS)
			elif event.keysym == "Right":
				self.xview(SCROLL, 1, UNITS)
			elif event.keysym in ("Return", "KP_Enter", "F2"):
				self.editHighlighted()
			elif event.keysym == "Escape":
				self.cardSelectClear()
			elif event.keysym in ("BackSpace", "Delete"):
				self.delSelectedCards()

			# Delete/Insert keys are handled by the global binding
			self.shift = False
			return

		insert = self.index(item, INSERT)
		text   = self.itemcget(item,"text")

		# Save undo info in case of too long wait
		if self.lastKey[0]==item and \
			event.time-self.lastKey[1] >= _KEY_TIMEOUT:	# 5 seconds?
				self.active.commit(True)

		if self.multiline: active.adjustTextBox(False)

		# editing
		if event.keysym == "BackSpace" or event.keysym=="Terminate_Server":
			if self.has_selection():
				self.dchars(item, SEL_FIRST, SEL_LAST)
				self.selectClear()
			else:
				if insert > 0:
					self.dchars(item, insert-1, insert-1)
			self.highlight(item)
			self.page.setInputModified()

		elif event.keysym == "Delete":
			if self.has_selection():
				self.dchars(item, SEL_FIRST, SEL_LAST)
				self.selectClear()
			else:
				self.dchars(item, insert, insert)
			self.page.setInputModified()

		# navigation
		elif event.keysym == "Home":
			pos = text.rfind('\n',0,insert)+1
			self.icursor(item, pos)
			self.editLeft(item, insert, pos)

		elif event.keysym == "End":
			try:
				pos = text.index('\n',insert)
				self.icursor(item, pos)
			except:
				self.icursor(item, END)
				pos = len(text)

			if self.shift:
				self.editRight(item, insert, pos)
			else:
				self.selectClear()
				self.selEditAnchor = -1

		elif event.keysym == "Right":
			pos = insert+1
			self.icursor(item, pos)
			self.editRight(item, insert, pos)

		elif event.keysym == "Left":
			pos = insert-1
			self.icursor(item, pos)
			self.editLeft(item, insert, pos)

		elif event.keysym == "Up":
			pl = text.rfind('\n',0,insert)
			if pl>=0:
				ppl = text.rfind('\n',0,pl-1)
				pos = insert - pl + ppl
				if pos>pl: pos = pl
				self.icursor(item, pos)
				self.editLeft(item, insert, pos)
			elif not self.shift:
				if self.multiline:
					if active.adjustTextBox(True):
						self.highlight(item)
				self.editUp()

		elif event.keysym == "Down":
			pl = text.rfind('\n',0,insert)
			nl  = text.find('\n',insert)
			if nl>0:
				nnl = text.find('\n',nl+1)
				pos = nl + insert - pl
				if nnl>0 and pos>nnl: pos = nnl
				self.icursor(item, pos)
				self.editRight(item, insert, pos)
			elif not self.shift:
				if self.multiline:
					if active.adjustTextBox(True):
						self.highlight(item)
				self.editDown()

		elif event.keysym in ("Return", "KP_Enter"):
			if self.multiline:	# Comment?
				if self.has_selection():
					self.dchars(item, SEL_FIRST, SEL_LAST)
					self.selectClear()
				self.insert(item, "insert", '\n')
				self.highlight(item)
			else:
				self.commit()
				active.editNext()

		elif event.keysym == "Escape":
			self.commit()
			active.editItem()

		# anything else editing
		elif event.char >= " ":
			# printable character
			if self.has_selection():
				self.dchars(item, SEL_FIRST, SEL_LAST)
				self.selectClear()

			self.insert(item, INSERT, event.char)
			self.highlight(item)
			self.page.setInputModified()

			if _BODY_INSERT:
				if event.char in ("+","-"):
					if active.card.tag=="REGION" and self.edit > 0:
						# Check next char
						try: next = text[insert]
						except: next = " "
						if " +-|)\n#".find(next)>=0:
							self.insertBody()

		if self.multiline and item==self.focus():
			if active.adjustTextBox(True):
				self.highlight(item)

		self._updateScrollBars()
		active.see(item)

		self.shift = False
		self.lastKey = (item, event.time)

	# ----------------------------------------------------------------------
	# POPUP MENU
	# ----------------------------------------------------------------------
	def _addCardLambda(self, tag):
		if isinstance(tag, str):
			ci = Input.CardInfo.get(tag)
		else:
			ci = tag
		return lambda t=ci.tag,s=self:s.addNewCardUndo(t)

	# ----------------------------------------------------------------------
	def _addTransformMenu(self, menu, cards):
		menu.add_command(label="ROT-DEFINI",
				underline='0',
				command=self._addCardLambda("ROT-DEFINI"),
				compound=LEFT, image=tkFlair.icons["ROT-DEFI"])
		menu.add_command(label=tkFlair._EXPANSION,
				underline='0',
				command=self.startend_expansion,
				compound=LEFT, image=tkFlair.icons["tags"])
		menu.add_command(label=tkFlair._TRANSLAT,
				underline='0',
				command=self.startend_translat,
				compound=LEFT, image=tkFlair.icons["tags"])
		menu.add_command(label=tkFlair._TRANSFORM,
				underline='1',
				command=self.startend_transform,
				compound=LEFT, image=tkFlair.icons["tags"])
		underline = tkFlair.initUnderline(cards)
		for c in cards:
			tkFlair.addCardCommand(menu, c,
				self._addCardLambda(c),
				underline)

	# ----------------------------------------------------------------------
	def _addGroupMenu(self, menu, group):
		# scan all cards belonging to a group
		cards = [x for x in list(Input.CardInfo._db.values()) if group in x.group]

		# Sort cards
#		cards.sort(Input.tagCmp)
		cards.sort(key=attrgetter("tag"))

		# Check for bodies
		if group == tkFlair._GEOMETRY_GROUP:
			bodiesmenu = Menu(menu, tearoff=1, activebackground=Ribbon._ACTIVE_COLOR)
			transmenu  = Menu(menu, activebackground=Ribbon._ACTIVE_COLOR)

			# remove transformations from the cards list
			trans = [x for x in cards if x.tag[0]=="$"]
			cards = [x for x in cards if x.tag[0]!="$" and len(x.tag)!=3 or x.tag=="END"]
#			trans.sort(Input.tagCmp)
			trans.sort(key=attrgetter("tag"))

			# Bodies list
			tkFlair.bodiesMenu(bodiesmenu, self.addNewCardUndo)

			# Transformation list
			self._addTransformMenu(transmenu, trans)

		elif group == "Bodies":
			tkFlair.bodiesMenu(menu, self.addNewCardUndo)
			return

		underline = tkFlair.initUnderline(cards)
		for c in cards:
			if c.tag=="#endinclude": continue
			if group == tkFlair._GEOMETRY_GROUP and c.tag=="END":
				menu.add_cascade(label="Bodies",
						menu=bodiesmenu,
						compound=LEFT, image=tkFlair.icons["TRC"],
						underline=0)
				menu.add_cascade(label="Transformations",
						menu=transmenu,
						compound=LEFT, image=tkFlair.icons["ROT-DEFI"],
						underline=0)
			tkFlair.addCardCommand(menu, c,
				self._addCardLambda(c),
				underline)

		if group == tkFlair._PREPROCESSOR_GROUP:
			menu.add_command(label=tkFlair._IFENDIF,
					underline='0',
					command=self.ifendif,
					compound=LEFT,
					image=tkFlair.icons[tkFlair._IFENDIF])

	# ----------------------------------------------------------------------
	# Create add a new card menu
	# ----------------------------------------------------------------------
	def _addCardMenu(self, menu, tags, submenu=True, prepro=True):
		if tags:
			# Specific cards
			tags.sort()

			# Check if all bodies are present remove them and
			# add a new bodies menu
			bodies = True
			for b in Input.BODY_TAGS:
				if b == "VOXELS": continue
				if b not in tags:
					bodies = False
					break
			if bodies:
				tags = [x for x in tags if len(x)!=3]
				tags.append("END")
				tags.sort()
				bodiesmenu = Menu(menu, tearoff=1, activebackground=Ribbon._ACTIVE_COLOR)
				tkFlair.bodiesMenu(bodiesmenu, self.addNewCardUndo)

			# remove transformations
			trans = [x for x in tags if x[0]=="$"]
			if trans:
				tags = [x for x in tags if x[0]!="$"]

				transmenu  = Menu(menu, activebackground=Ribbon._ACTIVE_COLOR)
				# Transformation list
				self._addTransformMenu(transmenu, trans)

			# Check remaining tags
			if (submenu or bodies) and len(tags)>4:
				cardmenu = Menu(menu, activebackground=Ribbon._ACTIVE_COLOR)
				menu.add_cascade(label="Add", underline=0,
					compound=LEFT, image=tkFlair.icons["add"],
					menu=cardmenu)
				add = True
			else:
				cardmenu = menu
				add = False

			underline = tkFlair.initUnderline(tags)
			for t in tags:
				tkFlair.addCardCommand(cardmenu, t,
					self._addCardLambda(t),
					underline)

			if bodies:
				menu.add_cascade(label="Bodies",
						menu=bodiesmenu,
						compound=LEFT,
						image=tkFlair.icons["TRC"],
						underline=0)
			if trans:
				menu.add_cascade(label="Transformations",
						menu=transmenu,
						compound=LEFT,
						image=tkFlair.icons["ROT-DEFI"],
						underline=0)

			# Add also the group menu that the cards belongs to
			if not add and tags:
				ci = Input.CardInfo.get(tags[0])
				for group in ci.group:
					groupmenu = Menu(menu, activebackground=Ribbon._ACTIVE_COLOR)
					menu.add_cascade(label=group, underline=0,
						compound=LEFT, image=tkFlair.icons["add"],
						menu=groupmenu)
					self._addGroupMenu(groupmenu, group)

		else:	# Groups
			if submenu:
				cardmenu = Menu(menu, activebackground=Ribbon._ACTIVE_COLOR)
				menu.add_cascade(label="Add", underline=0,
					compound=LEFT, image=tkFlair.icons["add"],
					menu=cardmenu)
			else:
				cardmenu = menu

			for group in tkFlair._groups_order:
				if group == "Developers" and not Input._developer: continue
				# Group Sub Menu
				groupmenu = Menu(cardmenu, activebackground=Ribbon._ACTIVE_COLOR)
				# Find underline
				underline = tkFlair._groups_underline[
							tkFlair._groups_order.index(group)]
				cardmenu.add_cascade(label=group, menu=groupmenu,
						compound=LEFT, image=tkFlair.icons["empty"],
						underline=underline)
				self._addGroupMenu(groupmenu, group)

		if prepro:
			premenu = Menu(menu, activebackground=Ribbon._ACTIVE_COLOR)
			menu.add_cascade(label=tkFlair._PREPROCESSOR_GROUP,
					compound=LEFT, image=tkFlair.icons["hash"],
					menu=premenu, underline=0)
			self._addGroupMenu(premenu, tkFlair._PREPROCESSOR_GROUP)

	# ----------------------------------------------------------------------
	# Filter menu
	# ----------------------------------------------------------------------
	def _addFilterMenu(self, menu):
		filtermenu = Menu(menu, activebackground=Ribbon._ACTIVE_COLOR)
		if InputExtra.filterMenu(self.page,
					filtermenu,
					self.active.card):
			menu.add_cascade(label="Filter on", menu=filtermenu,
					compound=LEFT, image=tkFlair.icons["filter"],
					underline=0)
		return

	# ----------------------------------------------------------------------
	# Create the popup to add a new card
	# ----------------------------------------------------------------------
	def popupInsertMenu(self, event=None):
		if self.active:
			self.commit()
			self.active.editItem()
		tags = self.page._tags
		if tags is not None and	len(tags)==1:
			self.addNewCardUndo(tags[0])
		else:
			menu = Menu(self, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
			self._addCardMenu(menu, tags, False)
			menu.tk_popup(self.winfo_pointerx(), self.winfo_pointery())
			return "break"	# is important otherwise the handlekey
					# sends the focus to the canvas

	# ----------------------------------------------------------------------
	# Create the popup menu
	# ----------------------------------------------------------------------
	def popupMenu(self, event):
		self.shift = False
		self.popupList.unPost()
		if int(event.type) == 4:	# Mouse event
			cw = self.findCardWidget(event.y)
			if cw is not None:
				# Find if is already selected
				if self.itemcget(cw.box, "fill") != _ACTIVE_COLOR:
					if self.active is not None:
						self.commit()
					self.click(event)
					cw.editItem()

		menu = Menu(self, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		active = self.active
		master = self.master

		menu.add_command(label="Cut", underline=2, command=self.cut,
				compound=LEFT, image=tkFlair.icons["cut"],
				accelerator="Ctrl-X")
		menu.add_command(label="Copy", underline=0, command=self.copy,
				compound=LEFT, image=tkFlair.icons["copy"],
				accelerator="Ctrl-C")
		menu.add_command(label="Paste", underline=0, command=self.paste,
				compound=LEFT, image=tkFlair.icons["paste"],
				accelerator="Ctrl-V")
		menu.add_command(label="Clone", underline=0, command=self.clone,
				compound=LEFT, image=tkFlair.icons["clone"],
				accelerator="Ctrl-D")

		menu.add_separator()

		self._addCardMenu(menu, self.page._tags)
		menu.add_command(label="Delete", underline=0,
				command=self.page.delete,
				compound=LEFT, image=tkFlair.icons["x"],
				accelerator="Del")

		if active is not None:
			if active.card is not None:
				menu.add_separator()
				self._addFilterMenu(menu)

				if not active.comment or self.selFrom>=0:
					menu.add_command(label='Insert Comment',
							compound=LEFT, image=tkFlair.icons["empty"],
							underline=0, command=self.insComment)
				if active.comment or self.selFrom>=0:
					menu.add_command(label='Delete Comment',
							compound=LEFT, image=tkFlair.icons["empty"],
							underline=2, command=self.delComment)
				menu.add_command(label="Edit Card", underline=0,
						compound=LEFT, image=tkFlair.icons["edit"],
						command=self.editCard)

			menu.add_separator()
			if active.card is None or (not active.card.enable) or self.selFrom>=0:
				menu.add_command(label='Enable Card(s)', underline=0,
					compound=LEFT, image=tkFlair.icons["empty"],
					command=lambda s=self:s.changeState(True))
			if active.card is None or (active.card.enable or self.selFrom>=0):
				menu.add_command(label='Disable Card(s)', underline=0,
					compound=LEFT, image=tkFlair.icons["empty"],
					command=lambda s=self:s.changeState(False))

			menu.add_separator()
			menu.add_command(label='Show Card(s)', underline=0,
				compound=LEFT, image=tkFlair.icons["enable"],
				command=lambda s=self:s.expandCard())
			menu.add_command(label='Hide Card(s)', underline=0,
				compound=LEFT, image=tkFlair.icons["disable"],
				command=lambda s=self:s.hideCard())

		menu.add_separator()
		menu.add_command(label='Help', command=self.help,
				compound=LEFT, image=tkFlair.icons["info"])

		menu.tk_popup(event.x_root, event.y_root)
		return "break"	# is important otherwise the handlekey
				# sends the focus to the canvas

	# ----------------------------------------------------------------------
	# Insert comment
	# ----------------------------------------------------------------------
	def insComment(self):
		if self.active is None or \
		   isinstance(self.active, HiddenCardWidget): return
		self.commit()
		self.active.editItem(None)
		self.active.insComment()
		self.active.see()
		self.active.editItem(-1)
		self.page.setInputModified()

	# ----------------------------------------------------------------------
	# Delete comment
	# ----------------------------------------------------------------------
	def delComment(self):
		if self.active is None: return
		self.commit()
		self.active.editItem(None)
		self.highlight("")

		undoinfo = []
		card_from, card_to = self.getSelectRange()
		if card_from is None: return
		for cw in self.cardlist[card_from : card_to]:
			for card in cw.cards():
				undoinfo.append(self.flair.setCommentUndo(card,""))
		self.addUndo(undo.createListUndo(undoinfo, "Delete Comment"))

	# ----------------------------------------------------------------------
	def toggleComment(self):
		if self.active is None or \
		   isinstance(self.active, HiddenCardWidget): return
		if self.active.comment:
			self.delComment()
		else:
			self.insComment()

	# ----------------------------------------------------------------------
	# Add start..stop cards around the selected ones
	# ----------------------------------------------------------------------
	def addStartStop(self, startTag, stopTag):
		self.shift = False
		card_from, card_to = self.getSelectRange()
		if card_from is None: return

		undoinfo = [(self.flair.refreshUndo,)]
		self.activate(self.cardlist[card_to-1])
		undoinfo.append(self.addNewCard(stopTag, edit=False, refresh=False))
		self.activate(self.cardlist[max(0,card_from-1)])
		undoinfo.append(self.addNewCard(startTag, edit=False, refresh=False))
		undoinfo.append((self.flair.refreshUndo,))
		self.addUndo(undo.createListUndo(undoinfo,"Add %s..%s"%(startTag,stopTag)))

		# re-indent cards
		for cw in self.cardlist[card_from+1:card_to+2]:
			cw.indent()

	# ----------------------------------------------------------------------
	# Add #include..#endinclude only if a range of cards is selected
	# otherwise add an empty #include..#endinclude
	# ----------------------------------------------------------------------
	def addIncludeCard(self):
		self.shift = False
		card_from, card_to = self.getSelectRange()
		if card_from is None: return

		empty = card_to-card_from==1
		undoinfo = [(self.flair.refreshUndo,)]
		self.activate(self.cardlist[max(0,card_from-1)])
		# on purpose in capital see below
		undoinfo.append(self.addNewCard("#INCLUDE", edit=False, refresh=False))
		if not empty: self.activate(self.cardlist[card_to])
		undoinfo.append(self.addNewCard("#endinclude", edit=False, refresh=False))
		undoinfo.append((self.flair.refreshUndo,))
		self.addUndo(undo.createListUndo(undoinfo,"Add #include"))
		self.activate(self.cardlist[max(0,card_from)])

		if not empty:
			# re-indent cards
			for cw in self.cardlist[card_from+1:card_to+2]:
				cw.indent()

	# ----------------------------------------------------------------------
	# add cards around the selected ones
	# #if ... #endif
	# ----------------------------------------------------------------------
	def ifendif(self):		self.addStartStop("#if", "#endif")
	def startend_expansion(self):	self.addStartStop("$start_expansion","$end_expansion")
	def startend_translat(self):	self.addStartStop("$start_translat", "$end_translat")
	def startend_transform(self):	self.addStartStop("$start_transform","$end_transform")

	# ----------------------------------------------------------------------
	# Add a new card and return undo information
	# ----------------------------------------------------------------------
	def addNewCard(self, tag, edit=True, refresh=True):
		"""add a new card given by tagname or by Card class"""

		if self.active is not None:
			self.commit()
		else:
			if len(self.cardlist)>0:
				self.activate(self.cardlist[0])
		active = self.active

		if isinstance(tag,str):
			# Special treatment
			if tag == tkFlair._IFENDIF:
				self.ifendif()
				return
			elif tag == tkFlair._EXPANSION:
				self.startend_expansion()
				return
			elif tag == tkFlair._TRANSLAT:
				self.startend_translat()
				return
			elif tag == tkFlair._TRANSFORM:
				self.startend_transform()
				return
			elif tag == "#INCLUDE":
				tag = "#include"

		self.cardSelectClear()
		project = self.project()

		# Find location to add card
		inside  = False
		if active is None or active.card is None:
			pos = project.input.bestPosition(tag)
			if active is None:
				cid = len(self.cardlist)
			else:
				if pos is None: pos = active.toPos()+1

				cid = self.cardlist.index(active)
				if pos <= active.fromPos():
					cid -= 1
				elif pos <= active.toPos():
					inside = True
				elif pos > active.toPos():
					pos = active.toPos()+1
		else:
			cid = self.cardlist.index(active)
			pos = active.card.pos()+1

		if isinstance(tag, Input.Card):
			card = tag.clone()
		else:
			card = Input.Card(tag)
			if _INSERT_COMMENT:
				card.setComment("Comment")

			# Special cases
			if card.tag == "TITLE" or \
			   card.tag=="GEOBEGIN" or \
			   card.tag=="PLOTGEOM":
				card.setExtra(project.title)

		undoinfo = self.flair.addCardUndo(card, pos, True)
		self.deactivate()
		self.addCard(card, None, cid+1)
		if len(self.cardlist)==1:
			cw = self.cardlist[0]
		else:
			cw = self.cardlist[cid+1]
		self.activate(cw)

		# Add card inside a range of hidden cards
		if inside:
			# Modify Hidden Card before...
			to_card = active.toCard
			active.toCard = project.input.cardlist[pos-1]
			active.refresh()
			self.addCard(project.input.cardlist[pos+1], to_card, cid+2)

		if refresh:
			self._updateScrollBars()
			self.update_idletasks()
			cw.see()
			self.focusIn()
			self.page.refreshStatus()
			self.page._refreshOtherFrames()

		if edit:
			# Call when system is idle
			self.after_idle(self.active.editNext)

		return undoinfo

	# ----------------------------------------------------------------------
	# Add a new card and add also the undo + refresh information
	# ----------------------------------------------------------------------
	def addNewCardUndo(self, tag, edit=True, refresh=True):
		if tag == "#include":
			self.addIncludeCard()
		else:
			undoinfo = [	(self.flair.refreshUndo,),
					self.addNewCard(tag, edit, refresh),
					(self.flair.refreshUndo,)]
			if undoinfo[1]:	# if addNewCard returned undo info
					# if..endif, start..stop do not return anything
				self.addUndo(
					undo.createListUndo(undoinfo, "New card %s"%(tag)))

	# ----------------------------------------------------------------------
	# Change Card type
	# ----------------------------------------------------------------------
	def changeTag(self, tag):
		if self.active is None: return
		self.commit()
		self.active.editItem()
		self.highlight("")

		card_from, card_to = self.getSelectRange()
		if card_from is None: return
		undoinfo = [self.flair.refreshUndo()]
		for cw in self.cardlist[card_from : card_to]:
			if cw.card:
				undoinfo.append(self.flair.changeTagUndo(cw.card, tag))
				cw.layout = Layout.getLayout(cw.card)
				cw.redraw()
			else:
				for c2 in cw.cards():
					undoinfo.append(self.flair.changeTagUndo(c2, tag))
				cw.refresh()
		undoinfo.append(self.flair.refreshUndo())
		self.active.updateLabel()
		self.cardSelect(card_from, card_to-1)
		self.addUndo(undo.createListUndo(undoinfo, "Change card type"))
		self.page.setInputModified()

	# ----------------------------------------------------------------------
	# Delete selected cards
	# ----------------------------------------------------------------------
	def delSelectedCards(self):
		if self.active is None: return
		card_from, card_to = self.getSelectRange()
		if card_from is None: return

		frompos = self.cardlist[card_from].fromPos()
		topos   = self.cardlist[card_to-1].toPos()

		# Delete the widgets
		self.delCard(card_from, card_to)

		# Prepare unfo info
		undoinfo = [(self.flair.refreshUndo,), self.flair.renumberUndo()]

		for c in range(topos, frompos-1, -1):
			undoinfo.append(self.flair.delCardUndo(c, False))

		undoinfo.append(self.flair.renumberUndo())
		undoinfo.append((self.flair.refreshUndo,))

		self.addUndo(undo.createListUndo(undoinfo, "Delete card(s)"))

		self._updateScrollBars()
		self.page.setInputModified()
		self.page.refreshStatus()
		self.page._refreshOtherFrames()
		self.input.clearCache("bodies")
		self.input.checkNumbering()	# XXX

	# ----------------------------------------------------------------------
	# Change card state
	# ----------------------------------------------------------------------
	def changeState(self, enable):
		if self.active is None: return
		card_from, card_to = self.getSelectRange()
		if card_from is None: return
		undoinfo = undo.UndoInfo()
		for cw in self.cardlist[card_from : card_to]:
			undoinfo.append(cw.setEnable(enable))
		self.addUndo(undoinfo.create("Change state"))
		self.page._refreshOtherFrames()

	# ----------------------------------------------------------------------
	def toggleState(self, event=None):
		if self.active is None: return
		card = self.input.cardlist[self.active.fromPos()]
		self.changeState(not card.enable)

	# ----------------------------------------------------------------------
	# Reset interrupt
	# ----------------------------------------------------------------------
	def resetInterrupt(self):
		self._interrupt = False
		return time.time()

	# ----------------------------------------------------------------------
	# Interrupt a long process from notify
	# ----------------------------------------------------------------------
	def interrupt(self):
		self.flair.notifyHide()
		self._interrupt = True

	# ----------------------------------------------------------------------
	# Show/Expand hidden cards
	# ----------------------------------------------------------------------
	def expandCard(self):
		if self.active is None: return
		card_from, card_to = self.getSelectRange()
		if card_from is None: return
		selCards = list(range(card_from, card_to))
		selCards.reverse()
		selFrom  = None
		selTo    = None
		first    = True
		startTime = self.resetInterrupt()
		updateCounter = 0
		self.deactivate()
		for i in selCards:
			cw = self.cardlist[i]
			if isinstance(cw, HiddenCardWidget):
				cards = cw.cards()
				self.delCard(i)
				pos    = i
				first2 = True
				for card in cards:
					if time.time()-startTime > _CARDS_TIMEOUT:
						if not self.flair.notifyIsVisible():
							self.flair.notify("Too long",
								"Process takes too long",
								tkFlair.NOTIFY_NORMAL,
								"Interrupt?", self.interrupt)
						updateCounter -= 1
						if updateCounter<=0:
							updateCounter = 100
							self.flair.update()
						if self._interrupt:
							cw = self.addCard(card, cards[-1], pos, first2, True)
							selfTo = cw
							break
					cw = self.addCard(card, None, pos, first2, card is cards[-1])
					if first2:
						selFrom = cw
						first2  = False
					selTo = cw
					pos += 1
			else:
				selFrom = cw
				if first: selTo = cw
			if self._interrupt: break
			first = False
		self.flair.notifyHide()
		self.activate(selFrom)
		if selFrom != selTo:
			self.cardSelect(selFrom, selTo)
		self._updateScrollBars()
#		print ">>>>>>>>>>>>>>>> expandCard t=",time.time()-startTime

	# ----------------------------------------------------------------------
	# Hide selected cards
	# ----------------------------------------------------------------------
	def hideCard(self):
		# remove selected cards from card list and replace them with
		# a HiddenCardWidget
		if self.active is None: return

		#undoinfo = []
		card_from, card_to = self.getSelectRange()
		if card_from is None: return

		if isinstance(self.cardlist[card_from], HiddenCardWidget):
			fromCard = self.cardlist[card_from].fromCard
		else:
			fromCard = self.cardlist[card_from].card

		if isinstance(self.cardlist[card_to-1], HiddenCardWidget):
			toCard = self.cardlist[card_to-1].toCard
		else:
			toCard = self.cardlist[card_to-1].card

		#undoinfo = [self.delCard(card_from, card_to)]
		self.delCard(card_from, card_to)
		#undoinfo.append(self.addCard(fromCard, toCard, card_from))
		self.addCard(fromCard, toCard, card_from)
		#self.addUndo(undo.createListUndo(undoinfo, "Hide Card(s)"))
		self.activate(self.cardlist[card_from])

		if self.active: self.active.see()
		self._updateScrollBars()

	# ----------------------------------------------------------------------
	# Show/Hide selected cards
	# ----------------------------------------------------------------------
	def toggleShow(self, event=None):
		card_from, card_to = self.getSelectRange()
		if card_from is None: return
		if isinstance(self.cardlist[card_from], HiddenCardWidget):
			self.expandCard()
		else:
			self.hideCard()

	# ----------------------------------------------------------------------
	# Edit card in dialog
	# ----------------------------------------------------------------------
	def editCard(self):
		if self.active is None or self.active.card is None: return
		self.commit()
		dialog = InputExtra.EditCardDialog(self, self.flair, self.active.card)
		if dialog.show():
			self.active.layout = Layout.getLayout(self.active.card)
			self.active.redraw()
			self.highlight()
			self.active.updateLabel()
			self.page._refreshOtherFrames()

	# ----------------------------------------------------------------------
	# Move selected cards up
	# ----------------------------------------------------------------------
	def moveUp(self, event=None):
		if self.active is None: return
		if self.edit is not None:
			self.commit()
			self.active.editItem()
		self.addUndo(self.moveUpCards(self.getSelectRange()))
		self.page.setInputModified()

	# ----------------------------------------------------------------------
	# Move selected cards down
	# ----------------------------------------------------------------------
	def moveDown(self,event=None):
		if self.active is None: return
		if self.edit is not None:
			self.commit()
			self.active.editItem()
		self.addUndo(self.moveDownCards(self.getSelectRange()))
		self.page.setInputModified()

	# ----------------------------------------------------------------------
	# Move range of cards(fromid,toid) up one card and return undo info
	# ----------------------------------------------------------------------
	def moveUpCards(self, cards):
		fromid, toid = cards
		if fromid is None: return
		if fromid==0: return undo.NullUndo

		self.itemconfig("active", width=_NORMAL_BORDER)

		topid = fromid - 1
		top = self._tagRange("top", topid)
		sel = self._tagRange("sel", fromid, toid)

		if self.selFrom >= 0:
			self.selFrom -= 1
			self.selTo   -= 1

		# Move canvas objects
		self.move("sel", 0, top[1]-sel[1])
		self.move("top", 0, sel[3]-top[3])

		# update cardlist
		cw = self.cardlist[topid]
		if toid >= len(self.cardlist):
			pos = self.cardlist[-1].toPos()+1
		else:
			pos = self.cardlist[toid].toPos()
		del self.cardlist[topid]
		self.cardlist.insert(toid-1, cw)

		# update input
		frompos = cw.fromPos()
		for card in cw.cards():
			self.input.moveCard(frompos, pos, False)
		self.input.renumber(frompos, pos+1)

		# re-indent cards
		for cw in self.cardlist[topid:toid]:
			cw.indent()

		self.itemconfig("active", width=_SEL_BORDER)
		if self.active is not None: self.active.see()

		self.input.checkNumbering() # XXX
		self.page._refreshOtherFrames()
		return ("Move up cards", self.moveDownCards, (fromid-1, toid-1))

	# ----------------------------------------------------------------------
	# Move range of cards(fromid,toid) down one card and return undo info
	# ----------------------------------------------------------------------
	def moveDownCards(self, cards):
		fromid, toid = cards
		if fromid is None: return
		if toid >= len(self.cardlist): return undo.NullUndo

		self.itemconfig("active", width=_NORMAL_BORDER)

		sel = self._tagRange("sel", fromid, toid)
		bot = self._tagRange("bot", toid)

		if self.selFrom>=0:
			self.selFrom += 1
			self.selTo   += 1

		# Move canvas objects
		self.move("sel", 0, bot[3]-sel[3])
		self.move("bot", 0, sel[1]-bot[1])

		# update cardlist
		cw  = self.cardlist[toid]
		pos = self.cardlist[fromid].fromPos()
		del self.cardlist[toid]
		self.cardlist.insert(fromid, cw)

		# update input
		pos2 = pos
		frompos = cw.fromPos()
		for card in cw.cards():
			self.input.moveCard(frompos, pos, False)
			frompos += 1
			pos += 1
		self.input.renumber(pos2, frompos+1)

		# re-indent cards
		for cw in self.cardlist[fromid:toid+1]:
			cw.indent()

		self.itemconfig("active", width=_SEL_BORDER)
		if self.active is not None: self.active.see()

		# XXX
		self.input.checkNumbering()
		self.page._refreshOtherFrames()
		return ("Move down cards", self.moveUpCards, (fromid+1, toid+1))

	# ----------------------------------------------------------------------
	# Find item
	# ----------------------------------------------------------------------
	def search(self, findStr, matchCase=False):
		if not matchCase: findStr = findStr.upper()

		self.selectClear()
		if self.active is not None:
			edit = self.edit
			item = self.focus()
			if item:
				pos = self.index(item,INSERT)
			else:
				pos  = 0
			self.commit()
			self.active.editItem(focus=False)
			cid = self.cardlist.index(self.active)
		else:
			edit = None
			pos  = 0
			cid  = 0

		for cw in self.cardlist[cid:]:
			if cw.search(findStr, matchCase, edit, pos):
				return
			edit = None
			pos  = 0
		for cw in self.cardlist[:cid]:
			if cw.search(findStr, matchCase, edit, pos):
				return
			edit = None
			pos  = 0

		self.flair.notify("Not Found",
			"Target \"%s\" not found" % (findStr),
			tkFlair.NOTIFY_WARNING)

	# ----------------------------------------------------------------------
	def replace(self, findStr, replaceStr, matchCase):
		#say("replace string...")
		pass

	# ----------------------------------------------------------------------
	def replaceAll(self, findStr, replaceStr, matchCase):
		#say("replaceAll string...")
		pass

	# ----------------------------------------------------------------------
	def help(self):
		if self.active is not None and self.active.card is not None:
			tag = self.active.card.tag
			if tag == "REGION":
				tag = ["8.2.7", "8.6"]
			elif tag == "LATTICE":
				tag = ["8.2.10", "8.7"]
			elif tag == "GEOEND":
				tag = ["8.2.11","8.8"]
			elif tag == "VOXELS":
				tag = ["8.9", "8.3"]
			elif tag in ("$start_expansion","$end_expansion"):
				tag = ["8.2.5.1", "8.4"]
			elif tag in ("$start_translat","$end_translat"):
				tag = ["8.2.5.2", "8.4"]
			elif tag in ("$start_transform","$end_transform"):
				tag = ["8.2.5.3", "8.4"]
			elif tag[0] == "$":
				tag = ["8.2.5", "8.4"]
			elif tag == "END":
				tag = ["8.2.6"]
			elif tag[0] == "!":
				tag = "F3"
			Manual.show(tag)
		else:
			Manual.show("F3")
		return "break"

	# ----------------------------------------------------------------------
	# For debugging dump cardlist
	# ----------------------------------------------------------------------
	#def printCardList(self):
	#	for i, cw in enumerate(self.cardlist):
	#		if isinstance(cw,HiddenCardWidget):
	#			say(i,cw.fromCard.tag,":",cw.fromPos(),"..",cw.toCard.tag,":",cw.toPos())
	#		else:
	#			say(i,cw.card.tag,":",cw.fromPos())

#===============================================================================
# Input Node
#===============================================================================
class _CardNode(tkTree.Node):
	def __init__(self, *args, **kw_args):
		tkTree.Node.__init__(self, *args, **kw_args)
		self.widget.tag_bind(self.label,  '<Button-1>', self.select)
		self.widget.tag_bind(self.symbol, '<Button-1>', self.select_image)

	# ----------------------------------------------------------------------
	# select node and show frame
	# ----------------------------------------------------------------------
	def select(self, event=None):
		tkTree.Node.select(self)
		self.PVT_enter()
		self.widget.move_cursor(self)
		self.findTags()

	# ----------------------------------------------------------------------
	def select_image(self, event=None):
		self.PVT_click(event)
		self.select(event)

	# ----------------------------------------------------------------------
	# show frame of node
	# ----------------------------------------------------------------------
	def findTags(self):
		nid = self.full_id()

		tag = None
		for idx in nid:
			if idx[0].isupper() or idx[0] in ('#',"!"):
				tag = idx

		filters = []
		for idx in nid:
			try: # Prepare filter
				# Filter format:
				#	type [function] \t what \t value
				# type:
				#	i=int
				#	f=float
				#	s=string
				# function:
				#	a=abs
				# what:
				#	[0-#]=whats
				#	s=sdum
				# value: anything int, float or string

				[f, w, v] = idx.split("\t")

				# just to speed up
				if f[0]=='i':
					v = int(v)
				elif f[0]=='f':
					v = float(v)
				try: w = int(w)
				except: pass

				filters.append((f, w, v, True))
			except ValueError:
				pass

		self.widget.page.filterCards(tag, filters)

#===============================================================================
# Add Card Group Menu
#===============================================================================
class _AddCardGroupMenu(FlairRibbon.FlairMenu):
	def createMenu(self):
		menu = Menu(self.page.ribbon, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		self.page.canvas._addCardMenu(menu, self.page._tags, False, False)
		return menu

#===============================================================================
class _AddCardAllMenu(FlairRibbon.FlairMenu):
	def createMenu(self):
		menu = Menu(self.page.ribbon, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		tkFlair.cardsMenu(menu, self.page.canvas.addNewCardUndo, tearoff=1, activebackground=Ribbon._ACTIVE_COLOR)
		menu.add_separator()
		tkFlair.cardsABCMenu(menu, self.page.canvas.addNewCardUndo, activebackground=Ribbon._ACTIVE_COLOR)
		return menu

#===============================================================================
class _AddCardPreprocessorMenu(FlairRibbon.FlairMenu):
	def createMenu(self):
		menu = Menu(self.page.ribbon, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		self.page.canvas._addGroupMenu(menu, tkFlair._PREPROCESSOR_GROUP)
		return menu

#===============================================================================
class _ChangeCardMenu(FlairRibbon.FlairMenu):
	def createMenu(self):
		menu = Menu(self.page.ribbon, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		tkFlair.cardsMenu(menu, self.page.canvas.changeTag, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		return menu

#===============================================================================
# Filter group menu
#===============================================================================
class _FilterGroup(FlairRibbon.FlairMenuGroup):
	#----------------------------------------------------------------------
	def createMenu(self):
		if self.page.canvas.active is None: return
		menu = Menu(self, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		InputExtra.filterMenu(  self.page,
					menu,
					self.page.canvas.active.card,
					self.page.findString.get(),
					self.page.findCase.get())
		return menu

#===============================================================================
# Input Tab Page
#===============================================================================
class InputPage(FlairRibbon.FlairPage):
	"""FLUKA input editor."""

	_name_ = "Input"
	_icon_ = "input"

	#-------------------------------------------------------------------------------
	def init(self):
		self.hash    = 0
		self.tags    = None
		self.filters = None
		# Remember for diff
		self.labelCard  = None
		self.labelWhat  = None
		self.labelText  = None
		self.labelError = None
		self._refreshOn = ("run","card")

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Input ==============
		group = Ribbon.LabelGroup(self.ribbon, "Input")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row=0,0
		b = Ribbon.LabelButton(group.frame,
				command=self.newInput,
				image=tkFlair.icons["new"],
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "New Input from the basic template")

		col,row=1,0
		b = FlairRibbon._TemplateMenuButton(group.frame, self,
				text="New",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "New input from template")

		# ---
		col,row=0,1
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["OPEN"],
				anchor=W,
				command=self.flair.loadInput,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Load existing input")

		# ---
		col,row=1,1
		b = Ribbon.LabelButton(group.frame,
				text="Load",
				anchor=W,
				command=self.flair.loadInput,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Load existing input")

		# ---
		col,row=0,2
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["save"],
				anchor=W,
				command=self.flair.saveInput,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Save input")

		col,row=1,2
		b = Ribbon.LabelButton(group.frame,
				text="Save",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				anchor=W,
				command=self.flair.saveInputAs,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Save input as")

		# ---
		col,row=2,0
		menulist = [	("Mcnp",     "mcnp",    self.flair.exportMcnp),
				("Povray",   "povray",  self.flair.exportPovray),
				("OpenSCAD", "openscad",self.flair.exportOpenscad) ]
		b = Ribbon.MenuButton(group.frame, menulist,
				text="Export "+Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				image=tkFlair.icons["export"],
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Export input to other formats")

		# ---
		col,row=2,2
		menulist = [	("Fluka", "input",   self.flair.importFluka),
				("Gdml",  "geant",   self.flair.importGdml),
				("Mcnp",  "mcnp",    self.flair.importMcnp) ]
		b = Ribbon.MenuButton(group.frame, menulist,
				text="Import "+Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				image=tkFlair.icons["import"],
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Import input from other formats")

		# ========== Cards ==============
		group = Ribbon.LabelGroup(self.ribbon, "Card")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row=0,0
		b = _AddCardGroupMenu(group.frame, self,
					image=tkFlair.icons["add32"],
					background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=3, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add Card from currently selected group [Ins] [Ctrl-Enter]")

		col,row=0,2
		b = _AddCardAllMenu(group.frame, self,
				text="Add",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add FLUKA Card")

		# ---
		col,row=1,0
		b = _AddCardPreprocessorMenu(group.frame, self,
				image=tkFlair.icons["define"],
				text="Preprocessor " + Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Add a preprocessor card")

		# ---
		col,row=1,1
		self._materialBtn = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["MATERIAL"],
				text="Material " + Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				compound=LEFT,
				anchor=W,
				command=self.addMaterial,
				background=Ribbon._BACKGROUND)
		self._materialBtn.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(self._materialBtn, "Add Material from database")

		# ---
		col,row=1,2
		b = _ChangeCardMenu(group.frame, self,
				image=tkFlair.icons["change"],
				text="Change " + Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Change card to another type")

		# ---
		col,row=2,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Delete",
				command=self.delete,
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Delete cards [Del]")

		# ---
		col,row=2,2
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["clone"],
				text="Clone",
				compound=LEFT,
				anchor=W,
				command=self.canvas.clone,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Clone cards [Ctrl-D]")

		# ========== Edit ==============
		group = Ribbon.LabelGroup(self.ribbon, "Edit")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row=0,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["toggle32"],
				command=self.toggle,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Change card state (Active/Disable) [Ctrl-L]")

		col,row=0,2
		menulist = [	("Enable",  "enable",  self.enableCard),
				("Disable", "disable", self.disableCard)]
		b = Ribbon.MenuButton(group.frame, menulist,
				text="State",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Change card state (Active/Disable) [Ctrl-L]")

		# ---
		col,row=1,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["expand"],
				command=self.toggleShow,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Expand or collapse cards [Ctrl-Double-click] [Ctrl-I]")

		col,row=2,0
		menulist = [	("Show",  "expand", self.expandCard),
				("Hide",  "shrink", self.hideCard)]
		b = Ribbon.MenuButton(group.frame, menulist,
				text="Show",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Expand or collapse cards [Ctrl-Double-click] [Ctrl-I]")

		# ---
		col,row=1,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["comment"],
				command=self.toggleComment,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Insert cards comment")

		col,row=2,1
		menulist = [	("Add Comment",    "comment_add",    self.insComment),
				("Remove Comment", "comment_delete", self.delComment)]
		b = Ribbon.MenuButton(group.frame, menulist,
				text="Comment",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add/Remove cards comment")

		# ---
		col,row=3,0
		f = Frame(group.frame)
		f.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NSEW)

		b = Ribbon.LabelButton(f, image=tkFlair.icons["up"],
				text="Move Up",
				compound=LEFT,
				command=self.canvas.moveUp,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.pack(side=TOP, expand=YES, fill=BOTH)
#		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Move up card [Ctrl-Up]")

		# ---
		col,row=3,1
		b = Ribbon.LabelButton(f, image=tkFlair.icons["down"],
				text="Move Down",
				compound=LEFT,
				command=self.canvas.moveDown,
				anchor=W,
				background=Ribbon._BACKGROUND)
		#b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		b.pack(side=BOTTOM, expand=YES, fill=BOTH)
		tkExtra.Balloon.set(b, "Move down card [Ctrl-Down]")

		# ---
		col,row=1,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["refresh"],
				text="Refresh",
				compound=LEFT,
				command=self.redraw,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Refresh Cards [Ctrl-R]")

		# ========== Filter ==============
		group = _FilterGroup(self.ribbon, self, "Filter")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()
		group.frame.grid_columnconfigure(0, weight=1)

		# ===
		col,row = 0,0
		self.runCombo = Ribbon.LabelCombobox(group.frame,
					width=24,
					command=self.runChanged)
		self.runCombo.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(self.runCombo, "Show only active cards for the selected run")

		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["run"],
				command=self.showActiveCards,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+2, padx=3, pady=0, sticky=E)
		tkExtra.Balloon.set(b, "Show only active cards for the selected run")

		# ===
		col,row = 0,1
		self.findString = tkExtra.LabelEntry(group.frame,
					"Search",
					"DarkGray",
					background="White", #Ribbon._BACKGROUND,
					width=20)
		self.findString.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(self.findString , "Find cards with string")
		self.findString.bind("<Return>",   self.find)
		self.findString.bind("<KP_Enter>", self.find)

		self.findCase = BooleanVar()
		self.findCase.set(False)
		b = Checkbutton(group.frame,
				variable=self.findCase,
				padx=0,
				pady=0,
				borderwidth=0,
				highlightthickness=0,
				font=Ribbon._FONT,
				activebackground="LightYellow",
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+1, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Case sensitive")

		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["find"],
				command=self.find,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+2, padx=3, pady=0, sticky=E)
		tkExtra.Balloon.set(b, "Find cards with string")

		# ===
		col,row = 0,2
		self.replaceString = tkExtra.LabelEntry(group.frame,
					"Replace",
					"DarkGray",
					state=DISABLED,
					background=Ribbon._BACKGROUND,
					width=20)
		self.replaceString.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(self.replaceString, "Replace with string")
		self.replaceString.bind("<Return>",   self.replace)
		self.replaceString.bind("<KP_Enter>", self.replace)

		b = Checkbutton(group.frame,
				state=DISABLED,
				padx=0,
				pady=0,
				borderwidth=0,
				highlightthickness=0,
				font=Ribbon._FONT,
				activebackground="LightYellow",
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+1, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Replace all matches")

		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["replace"],
				state=DISABLED,
				command=self.replace,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+2, padx=3, pady=0, sticky=E)
		tkExtra.Balloon.set(b, "Replace with string")

		# ========== Validation ===========
		if Input._developer:
			group = Ribbon.LabelGroup(self.ribbon, "Validation")
			group.pack(side=LEFT, fill=Y, padx=0, pady=0)

			# ---
			b = Ribbon.LabelButton(group.frame, text="Viewer",
					image=tkFlair.icons["viewer32"],
					compound=TOP,
					command=self.validate,
					anchor=W,
					background=Ribbon._BACKGROUND)
			b.pack(side=TOP, expand=YES, fill=BOTH)
			tkExtra.Balloon.set(b, "Validate input")

		# ========== View ===========
		group = Ribbon.LabelGroup(self.ribbon, "View")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		# ---
		b = Ribbon.LabelButton(group.frame, text="Viewer",
				image=tkFlair.icons["view"],
				compound=LEFT,
				command=self.viewer,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.pack(side=TOP, expand=YES, fill=BOTH)
		tkExtra.Balloon.set(b, "Open in internal viewer")

		# ---
		b = Ribbon.LabelButton(group.frame, text="Editor",
				image=tkFlair.icons["edit"],
				compound=LEFT,
				command=self.editor,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.pack(side=TOP, expand=YES, fill=BOTH)
		tkExtra.Balloon.set(b, "Open in external editor")

		# ---
		b = Ribbon.LabelButton(group.frame, text="Print",
				image=tkFlair.icons["print"],
				compound=LEFT,
				command=self.hardcopy,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.pack(side=TOP, expand=YES, fill=BOTH)
		tkExtra.Balloon.set(b, "Print input file [Ctrl-P]")

	#----------------------------------------------------------------------
	# Create Input Editor page
	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self, False)
		frame = self.frame

		self.splitter = tkExtra.HSplitter(frame,
					tkFlair.getInt(tkFlair._PAGE_SECTION,
						self.name+".split", 130),
						True)
		self.splitter.pack(expand=YES, fill=BOTH)
		treeFrame = self.splitter.leftFrame()
		cardFrame = self.splitter.rightFrame()

		# Tree
		self.createInputTreeList()
		self.tree = tkTree.Tree(master=treeFrame,
				root_id="/",
				root_label="Input",
				get_contents_callback=self.treeContents,
				width=300,
				font=tkFlair._TREE_FONT,
				node_class=_CardNode)
		self.tree.grid(row=0, column=0, sticky=NSEW)
		self.tree.page = self

		icon = tkFlair.icons["input"]
		self.tree.root.set_expanded_icon(icon)
		self.tree.root.set_collapsed_icon(icon)

		treeFrame.grid_rowconfigure(0, weight=1)
		treeFrame.grid_columnconfigure(0, weight=1)

		# Scroll bars
		sb=tkExtra.AutoScrollbar(treeFrame, takefocus=0)
		sb.grid(row=0, column=1, sticky=NS)
		self.tree.configure(yscrollcommand=sb.set)
		sb.configure(command=self.tree.yview)

		sb=tkExtra.AutoScrollbar(treeFrame,
				orient=HORIZONTAL, takefocus=0)
		sb.grid(row=1, column=0, sticky=EW)
		self.tree.configure(xscrollcommand=sb.set)
		sb.configure(command=self.tree.xview)

		self.tree.root.expand()

		# Canvas
		self.canvas = InputCanvas(cardFrame, self, cursor="")
		self.canvas.grid(row=0, column=0, sticky=NSEW)

		self.sbv = Scrollbar(cardFrame, takefocus=0)
		self.sbv.grid(row=0, column=1, sticky=NS)
		self.sbh = tkExtra.AutoScrollbar(cardFrame,
				orient=HORIZONTAL, takefocus=0)
		self.sbh.grid(row=1, column=0, sticky=EW)

		self.canvas.config( yscrollcommand=self.sbv.set,
				    xscrollcommand=self.sbh.set )
		self.sbv.config(command=self.canvas.yview)
		self.sbh.config(command=self.canvas.xview)

		# Place resize tab buttons
		self.tabs = []
		for c in range(_COLUMNS):
			self.tabs.append(
				Label(self.canvas,
					text=Unicode.MATH_UP_TACK,
					font=_COMMENT_FONT,
					foreground=_LABEL_COLOR,
					background="White",
					cursor="arrow",
					padx=0, pady=1,
					takefocus=False))
			self.tabs[-1].bind("<Button1-Motion>",lambda e,c=c:self.tabsMotion(e,c))
		self.placeTabs()

		cardFrame.grid_rowconfigure(0, weight=1)
		cardFrame.grid_columnconfigure(0, weight=1)

		# Card Interpretation
		self.labelFrame = Frame(cardFrame)
		self.label = Text(self.labelFrame, font=_FIXED_FONT,
				relief=SUNKEN,
				takefocus=1,
				#state=DISABLED,
				height=4,
				wrap=NONE)
		self.label.bind("<ButtonRelease-1>",self.labelFocus)
		#self.label.bind("<Double-1>",      self.editCard)
		self.label.bind("<Double-3>",       self.awari)
		sbv2 = Scrollbar(self.labelFrame, takefocus=0)
		sbh2 = tkExtra.AutoScrollbar(self.labelFrame, orient=HORIZONTAL, takefocus=0)

		# To initialize the grid
		self.label.grid(row=0,column=0, sticky=NSEW)
		sbv2.grid(row=0, column=1, sticky=NS)
		sbh2.grid(row=1, column=0, sticky=EW)

		self.label.config( yscrollcommand=sbv2.set,
				   xscrollcommand=sbh2.set )
		sbv2.config(command=self.label.yview)
		sbh2.config(command=self.label.xview)

		self.labelFrame.grid_columnconfigure(0, weight=1)

		if _SHOW_INTER:
			txt = Unicode.BLACK_UP_POINTING_TRIANGLE
		else:
			txt = Unicode.BLACK_DOWN_POINTING_TRIANGLE
		self.labelButton = Button(cardFrame, text=txt,
				command=self.labelButtonShow,
				font=tkFlair._TITLE_FONT,
				border=0,
				highlightthickness=0,
				width=2,
				activeforeground=tkFlair._TITLE_FOREGROUND_COLOR,
				takefocus=False, padx=1, pady=1)
		self.labelButtonShow()
		self.labelButton.bind("<B1-Motion>", self.labelButtonResize)
		tkExtra.Balloon.set(self.labelButton, "Hide/Show/Drag to Resize")

		b = Button(self.label, image=tkFlair.icons["edit"],
				command=self.editCard,
				cursor="arrow",
				takefocus=False)
		b.place(relx=1.0, x=-20, rely=0.0)
		tkExtra.Balloon.set(b, "Edit contents of current card")

	# ----------------------------------------------------------------------
	# Input I/O
	# ----------------------------------------------------------------------
	def newInput(self):
		self.flair.openInputTemplate(tkFlair.inputTemplate())

	# ----------------------------------------------------------------------
	def openTemplate(self, template):
		self.flair.openInputTemplate(template)

	#-------------------------------------------------------------------------------
	# generate Input Tree List from CardInfo
	#-------------------------------------------------------------------------------
	def createInputTreeList(self):
		#-----------------------------------------------------------------------
		# special treatment
		#-----------------------------------------------------------------------
		_inputSpecial = {
			"USRBIN"   : self.treeUSRXXX ,
			"USRBDX"   : self.treeUSRXXX ,
			"USRCOLL"  : self.treeUSRXXX ,
			"EVENTBIN" : self.treeUSRXXX ,
			"RESNUCLE" : self.treeUSRXXX ,
			"USRTRACK" : self.treeUSRXXX ,
			"USRYIELD" : self.treeUSRXXX }

		self.inputTreeList = []
		for group in tkFlair._groups_order:
			if group == "Developers" and not Input._developer: continue
			# scan all cards belonging to a group
			cards=list([c for c in list(Input.CardInfo._db.values()) if group in c.group])
			tags = " ".join([c.tag for c in cards])

			subList = []

			# special treatment for Geometry
			if group==tkFlair._GEOMETRY_GROUP:
				# remove them from the cards list
				bodies = list([x for x in cards if len(x.tag)==3 or x.tag[0]=="$"])
				#bodies.sort(lambda x,y:cmp(x.order, y.order))
				bodies.sort(key=attrgetter("order"))
				cards  = list([x for x in cards if len(x.tag)!=3 and x.tag[0]!="$"])

				# Append a special sub/category Bodies after GEOBEGIN
				bodyTags = " ".join([x.tag for x in bodies])

				# Bodies list
				bodyList = []
				for c in bodies:
					cname = string.capwords(c.name.lower(),"-")
					bodyList.append((cname, c.tag, None))

			#cards.sort(lambda x,y:cmp(x.order, y.order))
			cards.sort(key=attrgetter("order"))

			for c in cards:
				cname = string.capwords(c.name.lower(),"-")
				p = _inputSpecial.get(c.tag, None)
				subList.append((cname, c.tag, p))
				if c.tag=="GEOBEGIN":
					# Correct bodies
					subList.append(("Bodies", bodyTags, bodyList))

			self.inputTreeList.append((group, tags, subList))

	#===============================================================================
	# Tree build up
	#-------------------------------------------------------------------------------
	# default USRxxx cards
	#-------------------------------------------------------------------------------
	def treeUSRXXX(self, node):
		if self.project.input is None: return

		nid = node.full_id()
		tag = nid[2]
		try:
			cardlist = self.project.input.cards[tag]
		except KeyError:
			return

		# scan for the available units and particles
		if tag == "RESNUCLE":
			pwhat = None
			uwhat = 2
		else:
			pwhat = 2
			uwhat = 3

		units = []
		parts = []
		for card in cardlist:
			# Check if it exists
			u = abs(card.intWhat(uwhat))
			if u not in units: units.append(u)

			if pwhat is not None:
				p = card.what(pwhat)
				if p not in parts: parts.append(p)

		# sort units
		icon = tkFlair.icons["new"]
		units.sort()
		for i in units:
			idname = "ia\t%d\t%d" % (uwhat, i)
			self.tree.add_node(name=str(i),
					id=idname, flag=False,
					expanded_icon=icon, collapsed_icon=icon)

		# particles
		if pwhat is not None:
			icon = tkFlair.icons["particle"]
			parts.sort()
			for i in parts:
				idname = "s\t%d\t%s" % (pwhat, i)
				self.tree.add_node(name=i,
						id=idname, flag=False,
						expanded_icon=icon, collapsed_icon=icon)

	#-------------------------------------------------------------------------------
	def treeContents(self, node):
		treelist = None
		for nid in node.full_id():
			if nid == "/":
				treelist = self.inputTreeList
			else:
				if callable(treelist): break

				# scan list for target
				for name,idx,sublist in treelist:
					if nid==idx:
						treelist = sublist
						break
				#else:
				#	raise Exception("flukaTreeContents: no sublist found! for "+nid)

		if callable(treelist):
			treelist(node)
		else:
			nid = node.full_id()
			for name,idx,sublist in treelist:
				if len(nid)>=2 and self.project.input is not None:
					if name not in ("Bodies", "Editor"):
						if idx not in self.project.input.cards:
							continue

				if len(nid)>1 and sublist is None:
					icon = tkFlair.icons["tag"]
				else:
					icon = None
				try:
					icon = tkFlair.icons[idx]
				except KeyError:
					node.widget.add_node(name=name,
						id=idx, flag=(sublist!=None),
						expanded_icon=icon, collapsed_icon=icon)
				else:
					node.widget.add_node(name=name,
						id=idx, flag=(sublist!=None),
						expanded_icon=icon, collapsed_icon=icon)

	# ----------------------------------------------------------------------
	def placeTabs(self):
		w = float(self.canvas.winfo_width())
		for c,tab in enumerate(self.tabs):
			#tab.place(x=int(w*self.canvas._tabs[c]), y=0)
			tab.place(relx=self.canvas._tabs[c], y=0, anchor=CENTER)

	# ----------------------------------------------------------------------
	def tabsMotion(self, event, c):
		x = event.x_root - self.canvas.winfo_rootx()	# XXX bbox size
		x1,y1,x2,y2 = self.canvas.master.bbox('all')
		rel = float(x) / float(x2-x1)

		if c>0:
			rel = max(rel, self.canvas._tabs[c-1] + _MINTAB)
		else:
			rel = max(rel, _MINTAB)

		diff = rel - self.canvas._tabs[c]

		self.canvas._tabs[c] = rel
		for i in range(c+1,len(self.tabs)):
			self.canvas._tabs[i] += diff

		# Correct for right limits
		limit = 1.0 - _MINTAB
		for i in range(len(self.tabs)-1,c-1,-1):
			if self.canvas._tabs[i] > limit:
				self.canvas._tabs[i] = limit
			limit -= _MINTAB

		# Place tabs
		for i in range(c,len(self.tabs)):
			self.tabs[i].place(relx=self.canvas._tabs[i])

		self.canvas.resize()

	# ----------------------------------------------------------------------
	# label
	# ----------------------------------------------------------------------
	def labelButtonShow(self):
		global _SHOW_INTER
		if self.labelButton["text"] == Unicode.BLACK_DOWN_POINTING_TRIANGLE:
			_SHOW_INTER = False	# Don't be based on that if we have more InputPages open
			self.labelButton["text"] = Unicode.BLACK_UP_POINTING_TRIANGLE
			self.labelButton.place(in_=self.splitter.rightFrame(),
					relx=0.5, rely=1.0,
					anchor=S)
			self.labelFrame.grid_forget()
		else:
			_SHOW_INTER = True
			self.labelFrame.grid(row=2, column=0, columnspan=2, sticky=NSEW)
			self.labelButton["text"] = Unicode.BLACK_DOWN_POINTING_TRIANGLE
			self.labelButton.place(in_=self.labelFrame,
					relx=0.5, rely=0.0,
					anchor=N)

	# ----------------------------------------------------------------------
	def labelButtonResize(self, event):
		if event.y < -3:
			h = self.label["height"] + 1
			self.label["height"] = min(h,20)
		elif event.y > self.labelButton.winfo_height() + 3:
			h = self.label["height"] - 1
			self.label["height"] = max(h,1)

	#-------------------------------------------------------------------------------
	def updateLabel(self, card=None):
		isGeo = False
		if self.project is None: return
		if card is not None:
			if not card.isGeo() or card.tag in ("LATTICE","VOXELS"):
				fmt = self.project.input.format
			else:
				fmt = Input.FORMAT_FREE
				isGeo = True

			# if card is the same then do nothing
			cardText  = card.toStr(fmt)
			if card.invalid:
				cardError = card.errorMessage(
						self.canvas.layout.what2labels(
							self.canvas.active.layout))
			else:
				cardError = None
			if card is self.labelCard and \
			   self.labelText == cardText and \
			   self.labelError == cardError:
				return

		self.label.config(state=NORMAL)
		self.label.delete("0.0", END)

		if card is None:
			self.labelText  = None
			self.labelError = None
			self.labelCard  = card
			return

		# Show errors if needed
		if card.invalid:
			# what to labels conversion
			self.label.insert(END, cardError+"\n", "error")
			self.label.tag_configure("error",foreground="DarkRed")

		# Show alignment scale
		if _SHOW_SCALE and self.label["height"]>1:
			self.label.tag_delete("scale")
			if card.tag=="REGION":
				txt = Input.REGSCALE
			else:
				txt = Input.SCALE
			index = str(self.label.index(CURRENT))
			self.label.insert(END, txt)
			self.label.tag_add("scale",index,END)
			self.label.tag_configure("scale", foreground=_COMMENT_COLOR)

		# Insert text and highlight
		self.label.tag_delete("card")
		index = str(self.label.index(CURRENT))
		for line in cardText.splitlines():
			self.label.insert(END, line.ljust(80).decode()+"\n")
		self.label.tag_add("card",index,END)
		self.label.tag_configure("card", foreground="Black")

		if not isGeo:
			self.label.tag_delete("tag")
			self.label.tag_add("tag", index, "%s + %s chars"%(index,len(card.tag)+1))
			for l in range(1,(card.nwhats()-2)//6+1):
				#if card.tag != "TITLE" and card.tag != "GEOBEGIN":
				self.label.tag_add("tag", "%s + %d lines"%(index,l),
					"%s + %d lines + %s chars"%(index,l,len(card.tag)+1))
			self.label.tag_configure("tag", foreground=_TAG_COLOR)

		# find differences and mark them
		if card is self.labelCard:
			self.label.tag_delete("diff")
			start = None
			stop  = None
			self.labelText = self.labelText.ljust(len(cardText))
			for i in range(min(len(cardText), len(self.labelText))):
				new = cardText[i]
				old = self.labelText[i]
				if new!=old:
					if start is None: start = i
					stop  = i
				else:
					if start is not None:
						istart = '%s + %d chars' % (index,start)
						istop  = '%s + %d chars' % (index,stop+1)
						self.label.tag_add("diff",istart,istop)
						start = None
				if start is not None:
					istart = '%s + %d chars' % (index,start)
					istop  = '%s + %d chars' % (index,stop+1)
					self.label.tag_add("diff",istart,istop)
			self.label.tag_configure("diff",background=_DIFF_COLOR)
		else:
			self.label.tag_delete("what")
			self.labelWhat = None

		self.label.config(state=DISABLED)

		# Save card
		self.labelText  = cardText
		self.labelError = cardError
		self.labelCard  = card
		self.labelIndex = index

	# ----------------------------------------------------------------------
	# highlight label
	# ----------------------------------------------------------------------
	def highlightLabel(self, what):
		if self.labelCard is None: return
		if what == self.labelWhat: return
		if self.labelCard.isGeo() and \
			self.labelCard.tag not in ("LATTICE","VOXELS"): return

		self.label.tag_delete("what")
		if what is not None:
			if what == 0:
				self.label.tag_add("what",
					self.labelIndex+" + 70 chars",
					self.labelIndex+" + 1 lines")
			else:
				lines = (what-1) // 6
				if lines>0:
					linesstr = "+ %d lines " % (lines)
					what = (what-1) % 6 + 1
				else:
					linesstr = ""

				self.label.tag_add("what",
					"%s %s+ %d chars" % (self.labelIndex,linesstr,10*what),
					"%s %s+ %d chars" % (self.labelIndex,linesstr,10*(what+1)))
			self.label.tag_configure("what", background=_ACTIVE_COLOR)
		self.label.update_idletasks()
		self.labelWhat = what

	# ----------------------------------------------------------------------
	# canvas hardcopy
	# ----------------------------------------------------------------------
	def hardcopy(self, event=None):
		self.get()
		printer = tkFlair.printer(self.page)
		if printer is None: return

		canvas = self.canvas

		# Check canvas
		bbox = canvas.bbox("all")
		if bbox is None: return
		bbox = list(map(int, bbox))

		# Paper dimensions in cm
		margin = 1.0	# 1 cm margin from all sides
		cm2pt  = 72.0/2.54
		page_width, page_height = printer.PAPER_FORMAT[printer.paper]
		if printer.landscape:
			page_width, page_height = page_height, page_width

		# Find cm dimensions
		cmx  = canvas.create_rectangle(0, 0, "100c", "100c")
		cmbb = list(map(int, canvas.bbox(cmx)))
		canvas.delete(cmx)
		cmx  = (cmbb[2] - cmbb[0]) / 100
		cmy  = (cmbb[3] - cmbb[1]) / 100

		# Find number of pages
		width_cm  = float(bbox[2]-bbox[0]) / cmx
		height_cm = float(bbox[3]-bbox[1]) / cmy
		nop = math.ceil(height_cm / (page_height-4*margin))

		old_width = None

		if width_cm > (page_width-2.0*margin) * 0.9:
			old_width = bbox[2]-bbox[0]
			# Fake a resize event
			event = Event()
			event.width = int(0.9 * (page_width-2.0*margin) * cmx)
			# FIXME resize wont do anything!!!!
			# it gets the width from canvas.winfo_width()
			canvas.resize(event)
			canvas.update_idletasks()
			#bbox = canvas.bbox("all")
			#bbox = map(int, bbox)

		# Open file or a pip
		printer.open()
		printer.write("%!PS-Adobe-3.0\n")
		printer.write("%%%%Creator: %s-%s %s\n" \
			% (tkFlair.__name__, tkFlair.__version__, tkFlair.__author__))
		printer.write("%%%%Title: %s.inp: %s\n"%(self.project.inputName,
					self.project.title))
		printer.write("%%%%CreationDate: %s\n"%(time.ctime()))
		printer.write("%%%%DocumentPaperSizes: %s\n"%(printer.paper))
		printer.write("%%%%Pages: %d\n"%(nop))
		printer.write("%%EndComments\n")
		printer.write("/EncapDict 200 dict def EncapDict begin\n")
		printer.write("/showpage {} def /erasepage {} def /copypage {} def end\n")
		printer.write("/BeginInclude {0 setgray 0 setlinecap 1 setlinewidth\n")
		printer.write("0 setlinejoin 10 setmiterlimit [] 0 setdash\n")
		printer.write("/languagelevel where {\n")
		printer.write("  pop\n")
		printer.write("  languagelevel 2 ge {\n")
		printer.write("    false setoverprint\n")
		printer.write("    false setstrokeadjust\n")
		printer.write("  } if\n")
		printer.write("} if\n")
		printer.write("newpath\n")
		printer.write("save EncapDict begin} def\n")
		printer.write("/EndInclude {restore end} def\n")
		printer.write("%%EndProlog\n")

		PrologTarget = "%%BeginProlog\n"

		pageno = 1
		yoff   = 0

		while yoff < height_cm:
			printer.write("BeginInclude\n")

			if printer.landscape:
				# Header
				printer.write("%%%%Page: %d %d\n" % (pageno, pageno))
				printer.write("gsave\n")
				printer.write("/Times-Bold findfont 12 scalefont setfont\n")
				printer.write("90 rotate %d %d moveto ( %s:   %s ) show\n" % \
					(margin*cm2pt, -1.8*margin*cm2pt,
					 self.project.inputFile, self.project.title))
				printer.write("newpath 0.25 setlinewidth "
						"%d %d moveto %d %d lineto stroke\n" % \
					(margin*cm2pt, -2*margin*cm2pt,
					(page_width-margin)*cm2pt, -2*margin*cm2pt))

				# Page number
				printer.write("/Times findfont 12 scalefont setfont\n")
				printer.write("%d %d moveto ( - %d - ) show\n" % \
					((page_width/2.0-margin)*cm2pt,
					 -(page_height-margin)*cm2pt, pageno))
				printer.write("newpath 0.25 setlinewidth "
					"%d %d moveto %d %d lineto stroke\n" % \
					(margin*cm2pt, -(page_height-margin)*cm2pt+15,
					(page_width-margin)*cm2pt, -(page_height-margin)*cm2pt+15))
				printer.write("grestore\n")

				ps = canvas.postscript(
					colormode="color",
					rotate=1,
					x=0, y="%dc"%(yoff),
					pageanchor=NW,
					pagex="%dc"%(2*margin),
					pagey="%dc"%(margin),
					width="%dc"%(page_width-margin),
					height="%dc"%(page_height-4*margin),
					pagewidth="%dc"%(page_width-margin),
					pageheight="%dc"%(page_height-margin))
			else:
				# Header
				printer.write("%%%%Page: %d %d\n" % (pageno, pageno))
				printer.write("gsave\n")
				printer.write("/Times-Bold findfont 12 scalefont setfont\n")
				printer.write("%d %d moveto ( %s:   %s ) show\n" % \
					(margin*cm2pt, (page_height-1.8*margin)*cm2pt,
					 self.project.inputFile, self.project.title))
				printer.write("newpath 0.25 setlinewidth "
						"%d %d moveto %d %d lineto stroke\n" % \
					(margin*cm2pt, (page_height-2*margin)*cm2pt,
					(page_width-margin)*cm2pt, (page_height-2*margin)*cm2pt))

				# Page number
				printer.write("/Times findfont 12 scalefont setfont\n")
				printer.write("%d %d moveto ( - %d - ) show\n" % \
					((page_width/2.0-margin)*cm2pt, (margin)*cm2pt,
					 pageno))
				printer.write("newpath 0.25 setlinewidth "
					"%d %d moveto %d %d lineto stroke\n" % \
					(margin*cm2pt, margin*cm2pt+15,
					(page_width-margin)*cm2pt, margin*cm2pt+15))
				printer.write("grestore\n")

				ps = canvas.postscript(
					colormode="color",
					x=0, y="%dc"%(yoff),
					pageanchor=NW,
					pagex="%dc"%(margin),
					pagey="%dc"%(page_height-2*margin),
					width="%dc"%(page_width-margin),
					height="%dc"%(page_height-4*margin),
					pagewidth="%dc"%(page_width-margin),
					pageheight="%dc"%(page_height-margin))
			EOC = ps.find(PrologTarget) + len(PrologTarget)
			EOF = ps.find("%%EOF\n")
			for line in ps[EOC:EOF].split("\n"):
				if len(line)>0 and line[0] == "%": continue
				printer.write(line+"\n")
			printer.write("EndInclude showpage\n")
			yoff += page_height - 4*margin
			pageno += 1
		printer.write("%%EOF\n")

		# Restore width if needed
		if old_width:
			event.width = old_width
			canvas.resize(event)

		printer.close()
		self.setStatus("Input printed",  "DarkBlue")

	# ----------------------------------------------------------------------
	# delete all
	# ----------------------------------------------------------------------
	def deleteAll(self):
		if self.page is None: return
		undoinfo = self.canvas.deleteAll()
#		if len(undoinfo)>0:
#			self.canvas.addUndo(undo.createListUndo(undoinfo))
		self.updateLabel()

	# ----------------------------------------------------------------------
	# Activate and set focus
	# ----------------------------------------------------------------------
	def activate(self):
		FlairRibbon.FlairPage.activate(self)
		self.canvas.focus_set()

	# ----------------------------------------------------------------------
	# redraw all widgets
	# ----------------------------------------------------------------------
	def redraw(self):
		self.hash = 0
		self.showCards(self.tags, self.filters, False)

	# ----------------------------------------------------------------------
	# Refresh canvas
	# ----------------------------------------------------------------------
	def refresh(self):
		FlairRibbon.FlairPage.refresh(self)
		self.fillRuns()
		self.tree.refresh()
		self.canvas.setTabs()
		if self.canvas.active is not None and self.canvas.edit is not None:
			self.canvas.commit()
		self.busy()
		self.canvas.layout.cache(True)	# Enable caching of lists
		self.canvas.refresh()
		self.canvas.layout.cache(False)
		self.notBusy()

	# ----------------------------------------------------------------------
	# Fill the runCombo with active runs from the flair project
	# ----------------------------------------------------------------------
	def fillRuns(self):
		activeRun = self.runCombo.get()
		self.runCombo.clear()
		self.runCombo.insert(END, _ALL)
		found = activeRun == _ALL
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
	# refresh status
	# ----------------------------------------------------------------------
	def refreshStatus(self):
		FlairRibbon.FlairPage.refreshStatus(self)
		active = self.canvas.active
		if active:
			if isinstance(active,HiddenCardWidget):
				num = "%d-%d"%(active.fromCard.pos()+1,
					active.toCard.pos()+1)
			else:
				num = str(active.card.pos()+1)
		else:
			num = "-"
		# Count cards
		ncards = 0
		for cw in self.canvas.cardlist:
			if isinstance(cw,CardWidget): ncards += 1

		card_from, card_to = self.canvas.getSelectRange()
		if card_from is None or card_to - card_from <= 1:
			selected = ""
		else:
			selected = "Selected:%d"%(card_to-card_from)

		if ncards != len(self.project.input.cardlist):
			self.setStatus("Active:%s %s Displayed:%d Total:%d" \
				% (num, selected, ncards, len(self.project.input.cardlist)))
		else:
			self.setStatus("Active:%s %s Total:%d" \
				% (num, selected, len(self.project.input.cardlist)))

	# ----------------------------------------------------------------------
	def addMaterial(self):
		x = self._materialBtn.winfo_rootx()
		y = self._materialBtn.winfo_rooty() + self._materialBtn.winfo_height()
		if not InputExtra.MaterialList.shown:
			mat = InputExtra.MaterialList(self.ribbon).show(x,y)
			if mat:
				Materials.insert2Input(self.flair, mat)
		else:
			InputExtra.MaterialList.shown.destroy()

	# ----------------------------------------------------------------------
	# Refresh other frames of flair
	# ----------------------------------------------------------------------
	def _refreshOtherFrames(self):
		self.flair.refreshGeoedit()

	# ----------------------------------------------------------------------
	# filter card and add card to canvas
	# the filters are defined as lists of lists of tuples
	# filters = [ [list-of-OR] AND [list-of-OR] AND ...]
	# [list-of-OR] = (f,w,v)...
	# ----------------------------------------------------------------------
	def filterCard(self, cardlist, card, filters):
		if filters and len(filters)>0:
			inp = self.project.input
			for flst in filters: # list of AND's
				if isinstance(flst,tuple):
					f,w,v,c = flst
					if not card.test(f,w,v,c): return
				else: # list of OR's
					res_or = False	# result of OR
					for (f,w,v,c) in flst:
						if card.test(f,w,v,c):
							res_or = True
							break
					if not res_or: return
			cardlist.append(card)
		else:
			cardlist.append(card)

	# ----------------------------------------------------------------------
	# Filter and show cards
	# ----------------------------------------------------------------------
	def filterCards(self, tags=None, filters=None):
		self.canvas.setTabs()
		self.showCards(tags, filters)

	# ----------------------------------------------------------------------
	def filterOn(self, f):
		case = self.findCase.get()
		self.showCards(None,[('as',f[0],f[1],case)])

	# ----------------------------------------------------------------------
	# Show project cards
	# ----------------------------------------------------------------------
	def showCards(self, tags=None, filters=None, addUndo=True):
		"""tags - tags as words to show
		filters - commands see main program"""
#		print
#		print "="*80
#		import traceback; traceback.print_stack()
#		print "="*80

		if self.project.input is not self.canvas.input:
			# Set the input file to use
			self.canvas.updateInput()

		# Measure time to show
		startTime = time.time()

		# Validate last item
		if self.canvas.active is not None:
			if self.canvas.edit is not None:
				self.canvas.commit()
				self.canvas.active.editItem()
			self.canvas.cardSelectClear()

		inp = self.project.input
		if inp is None:
			self.deleteAll()
			self.canvas.active = None
			self.label.config(state=NORMAL)
			self.label.delete('0.0', END)
			self.label.config(state=DISABLED)
			return

		# Prepare hash from parameters
		h = hash(self.project.inputFile) ^ hash(self.project.input) \
			^ hash(tags) ^ len(self.project.input.cardlist)
		if filters:
			for lst in filters:
				if isinstance(lst,list):
					for (f,w,v) in lst:
						h ^= hash(f) \
						  ^  hash(w) \
						  ^  hash(v)
				else:
					h ^= hash(lst[0]) \
					  ^  hash(lst[1]) \
					  ^  hash(lst[2])

		# Is the same page?
		if self.hash == h: return
		self.hash    = h
		self.tags    = tags
		self.filters = filters

		# Exclude cards that are not active
		run = self.runCombo.get()
		if run != _ALL:
			run = self.project.getRunByName(run)
			if run:
				defines = run.defines
			else:
				defines = None
			errors = self.project.input.preprocess(defines)
			if errors:
				for card,err in errors:
					self.flair.write("\nERROR: %s\n"%(card))
					self.flair.write(err)

				self.flair.notify("Pre-process errors",
					"There are several preprocessor errors",
					tkFlair.NOTIFY_WARNING,
					"Output",self.flair.showOutput)
			ignore = True
		else:
			ignore = False

		# Fill with visible cards
		self.canvas.layout.cache(True)	# Enable caching of lists
		self.busy()
		cardlist = []
		if tags is None:
			self._tags = None
			for card in inp.cardlist:
				if ignore and not card.active: continue
				self.filterCard(cardlist, card, filters)
		else:
			self._tags = tags.split()
			for t in self._tags:
				try:
					for card in inp.cards[t]:
						if ignore and not card.active: continue
						self.filterCard(cardlist, card, filters)
				except KeyError:
					pass

			# Add preprocessor cards unfiltered if needed
			if _SHOW_PREPROCESSOR:
				for t in [x.tag for x in list(Input.CardInfo._db.values()) \
						if tkFlair._PREPROCESSOR_GROUP in x.group ]:
					# Is already added?
					if t in tags: continue
					try:
						for card in inp.cards[t]:
							cardlist.append(card)
					except KeyError:
						pass
		# Sort list
#		cardlist.sort(Input.Card.cmpPos)
		cardlist.sort(key=lambda card: card.pos())

		# Remove preprocessor cards and subcontent when everything is hidden
		rmUnwanted = False
		if _SHOW_PREPROCESSOR:
			if tags and tags[0][0]=='#' and "#include" not in self.project.input.cards:
				rmUnwanted = False
			else:
				rmUnwanted = True
		preerror = False
		if rmUnwanted:
			start   = [-1]	# starting position of block
			count   = [0]	# line count of disabled cards
			include = [-1]	# expanded or hidden include lines
			i = 0
			while i < len(cardlist):
				card = cardlist[i]
				tag  = card.tag
				if not card.enable:
					# Increase count
					count[-1] += 1
				elif tag in ("#if", "#ifdef", "#ifndef", "#include"):
					# Start new branch
					start.append(i)
					count.append(0)
					if tag == "#include":
						include.append(i)
				elif tag in ("#elif", "#else"):
					pass
				elif tag in ("#endif", "#endinclude"):
					if len(start)==1:
						# Error in input TOO many endif
						if not preerror:
							self.flair.write("ERROR Too many closing #endif Card:%s\n"%(i+1))
							self.flair.notify("Preprocessing error",
								"Too many closing #endif Card:%s"%(i+1),
								tkFlair.NOTIFY_ERROR,
								"Output",self.flair.showOutput)
							preerror = True
					else:
						# check last count
						s = start.pop()
						c = count.pop()
						if tag=="#endinclude":
							inc = include.pop()
							# FIXME this should be treated for only the EXPANDED cards
							# Maybe I have to put a special icon
							if len(include)==1:
								# remove cards
								s += 1
								del cardlist[s:i]
								i = s + 1
								continue
						elif c == 0:
							# remove cards
							del cardlist[s:i+1]
							i = s
							continue
						else:
							# keep the card and move count
							count[-1] += c
				else:
					# Increase count
					count[-1] += 1
				i += 1

			if len(start)>1:
				if not preerror:
					self.flair.write("ERROR Too many #if\n")
					self.flair.notify("Preprocessing error",
						"Too many #if",
						tkFlair.NOTIFY_ERROR,
						"Output",self.flair.showOutput)

		# Display cards
		self.canvas.displayCards(cardlist)
		del cardlist
		self.notBusy()
		self.canvas.layout.cache(False)	# Disable caching of lists

		# Measure time to show
		if Input._developer:
			dt = time.time() - startTime
			if dt > 1.0:
				self.flair.write("InputPage: draw time= %5.2gs\n"%(dt))

	# ----------------------------------------------------------------------
	# See active card
	# ----------------------------------------------------------------------
	def seeActive(self):
		if self.canvas.active:
			self.canvas.active.see()

	# ----------------------------------------------------------------------
	def showInvalid(self):
		cardlist = [ x.card for x in self.canvas.cardlist
				if x.card is not None ]
		if len(cardlist) == 0: self.showCards()

		# Add a new view
		self.hash = 0

		# Remove invalid cards
		for i in range(len(cardlist)-1,-1,-1):
			card = cardlist[i]
			if not card.invalid:
				if cardlist[i] == self.canvas.active:
					self.canvas.active = None
				del cardlist[i]
		self.canvas.displayCards(cardlist)

	# ----------------------------------------------------------------------
	def afterDock(self):
		self.hash = 0

	# ----------------------------------------------------------------------
	def togglePreprocessor(self):
		global _SHOW_PREPROCESSOR
		_SHOW_PREPROCESSOR = not _SHOW_PREPROCESSOR
		tkFlair.setBool(tkFlair._INPUT_SECTION, "showpreprocessor", _SHOW_PREPROCESSOR)
		self.hash = 0
		self.showCards(self.tags, self.filters)

	# ----------------------------------------------------------------------
	# Run has changed show active cards
	# ----------------------------------------------------------------------
	def runChanged(self):
		self.showCards(self.tags, self.filters)
		if self.runCombo.get() != _ALL: self.showActiveCards()

	# ----------------------------------------------------------------------
	# Filter only active cards
	# ----------------------------------------------------------------------
	def showActiveCards(self):
		cardlist = [ x.card for x in self.canvas.cardlist
				if x.card is not None ]

		if len(cardlist) == 0: self.showCards()
		# Add a new view
		self.hash = 0

		# Exclude cards that are not active
		run = self.runCombo.get()
		if run != _ALL:
			run = self.project.getRunByName(run)
			if run:
				defines = run.defines
			else:
				defines = None
			errors = self.project.input.preprocess(defines)
			if errors:
				for card,err in errors:
					self.flair.write("\nERROR: %s\n"%(card))
					self.flair.write(err)

				self.flair.notify("Pre-process errors",
					"There are several preprocessor errors",
					tkFlair.NOTIFY_WARNING,
					"Output",self.flair.showOutput)

			for i in range(len(cardlist)-1,-1,-1):
				card = cardlist[i]
				#if card.tag[0]!='#' and card.ignore():
				if card.tag[0]!='#' and not card.active:
					if cardlist[i] == self.canvas.active:
						self.canvas.active = None
					del cardlist[i]
			self.canvas.displayCards(cardlist)
		else:
			self.showCards(self.tags, self.filters)

	# ======================================================================
	# Undo
	# ======================================================================
	def canUndo(self): return self.flair.canUndo()
	def canRedo(self): return self.flair.canRedo()
	def resetUndo(self):
		self.deleteAll()
		self.flair.undoredo.reset()

	# ----------------------------------------------------------------------
	def undo(self, event=None): self.flair.undoInput()
	def redo(self, event=None): self.flair.redoInput()

	# ----------------------------------------------------------------------
	def commit(self):
		"""Commit changes in canvas"""
		if self.page is not None and self.canvas.active is not None:
			self.canvas.commit()
			self.canvas.active.editItem()

	# ----------------------------------------------------------------------
	# Add a new card or body if we are in the region expression
	# ----------------------------------------------------------------------
	def add(self, event=None):
		canvas = self.canvas
		active = canvas.active
		if active and active.card and \
				active.card.tag=="REGION" and canvas.edit > 0:
			canvas.insertBody()
		else:
			self.canvas.popupInsertMenu()

	# ----------------------------------------------------------------------
	def delete(self, event=None):
		item = self.canvas.focus()
		if not item: self.canvas.delSelectedCards()
		# else is handled by _handleKey

	# ----------------------------------------------------------------------
	def toggle(self):                 self.canvas.toggleState()
	def toggleShow(self, event=None): self.canvas.toggleShow()
	def toggleComment(self):	  self.canvas.toggleComment()

	def ifendif(self):	self.canvas.ifendif()
	def selectAll(self):	self.canvas.selectAll()
	def unSelectAll(self):	self.canvas.unSelectAll()
	def insComment(self):	self.canvas.insComment()
	def delComment(self):	self.canvas.delComment()
	def enableCard(self):	self.canvas.changeState(True)
	def disableCard(self):	self.canvas.changeState(False)
	def expandCard(self):	self.canvas.expandCard()
	def hideCard(self):	self.canvas.hideCard()

	# ----------------------------------------------------------------------
	def cut(self, event=None):	self.canvas.cut(event)

	# ----------------------------------------------------------------------
	def copy(self, event=None):
		focus = self.page.focus_get()
		if focus is self.canvas:
			self.canvas.copy(event)
		elif focus is self.label:
			self.labelCopy()

	def paste(self, event=None):	self.canvas.paste(event)
	def clone(self, event=None):	self.canvas.clone(event)
	def help(self):			self.canvas.help()
	def editCard(self,event=None):	self.canvas.editCard()

	# ----------------------------------------------------------------------
	def labelFocus(self, event=None):
		self.label.focus_set()
		return "break"

	# ----------------------------------------------------------------------
	def labelCopy(self, event=None):
		self.page.clipboard_clear()
		try:
			self.page.clipboard_append(self.label.get("sel.first", "sel.last"))
		except TclError:
			pass
		return "break"

	# ----------------------------------------------------------------------
	def find(self, event=None):
		txt = self.findString.get()
		if not txt: return
		case = self.findCase.get()
		self.canvas.search(txt, case)

	# ----------------------------------------------------------------------
	def replace(self, event=None):
		src = self.findString.get()
		if not src: return
		dst = self.replaceString.get()
		if not dst: return
		self.canvas.replace(src, dst)

	# ----------------------------------------------------------------------
	def editor(self):
		if self.project.isModified():
			if tkFlair.askyesno("Save Project",
					"Current Project is modified\n" \
					"Save before opening?", parent=self.page):
				self.flair.saveProject()
		tkFlair.editor([self.project.inputFile])

	# ----------------------------------------------------------------------
	def viewer(self):
		if self.project.isModified():
			if tkFlair.askyesno("Save Project",
					"Current Project is modified\n" \
					"Save before opening?", parent=self.page):
				self.flair.saveProject()
		self.flair.viewer([self.project.inputFile])

	# ----------------------------------------------------------------------
	def get(self):
		if self.page is None: return
		if self.canvas.active is not None:
			self.canvas.commit()

	# ----------------------------------------------------------------------
	def setInputModified(self):
		self.flair.setInputModified()

	# ----------------------------------------------------------------------
	def transformGeometry(self, event=None):
		if self.canvas.active is None: return
		self.canvas.commit()
		self.canvas.active.editItem()

		cardlist = []
		card_from, card_to = self.canvas.getSelectRange()
		if card_from is None: return
		for c in self.canvas.cardlist[card_from : card_to]:
			if c.card is not None:
				cardlist.append(c.card)
			else:
				cardlist.extend(c.cards())

		GeometryBody.TransformGeometryDialog(self.page, self.flair, cardlist)

	# ----------------------------------------------------------------------
	# Validate input
	# ----------------------------------------------------------------------
	def validate(self):
		validate = Validate.Validate(self.project, self.flair.log)
		run = self.runCombo.get()
		if run != _ALL:
			run = self.project.getRunByName(run)
		else:
			run = None
		validate.check(run)

	# ----------------------------------------------------------------------
	def toggleSplitter(self, event=None):
		self.splitter.toggle()

	# ----------------------------------------------------------------------
	def awari(self, event):
		iniFil = os.path.join(tkFlair.iniDir,"awari.var")
		imgDir = os.path.join(tkFlair.prgDir, "icons")
		awari.startGame(self.master, iniFil, imgDir)

	# ----------------------------------------------------------------------
	def configSave(self):
		section = tkFlair._INPUT_SECTION
		tkFlair.config.set(section, "showcard",	_SHOW_INTER)

#-------------------------------------------------------------------------------
def configGet():
	global _ACTIVE_COLOR, _NORMAL_COLOR, _TAG_COLOR, _DEFINE_COLOR, \
	       _LABEL_COLOR, _VALUE_COLOR, _COMMENT_COLOR, _DISABLE_COLOR, \
	       _DIFF_COLOR, _ERROR_COLOR, _HIDDEN_COLOR, _SDUM_COLOR, \
	       _LABEL_FONT, _COMMENT_FONT, _VALUE_FONT, _TAG_FONT, \
	       _FIXED_FONT, _HIDDEN_FONT, _INVALID_COLOR, \
	       _PAD, _COLUMNS, _KEY_TIME_THRESHOLD, \
	       _INSERT_COMMENT, _SHOW_SCALE, _SHOW_PREPROCESSOR, \
	       _BODY_INSERT, _SHOW_INTER, _CARDS_TIMEOUT

	_ACTIVE_COLOR  = tkFlair.getColor("input.active",   _ACTIVE_COLOR)
	_COMMENT_COLOR = tkFlair.getColor("input.comment",  _COMMENT_COLOR)
	_DEFINE_COLOR  = tkFlair.getColor("input.define",   _DEFINE_COLOR)
	_DIFF_COLOR    = tkFlair.getColor("input.diff",     _DIFF_COLOR)
	_DISABLE_COLOR = tkFlair.getColor("input.disable",  _DISABLE_COLOR)
	_ERROR_COLOR   = tkFlair.getColor("input.error",    _ERROR_COLOR)
	_HIDDEN_COLOR  = tkFlair.getColor("input.hidden",   _HIDDEN_COLOR)
	_INVALID_COLOR = tkFlair.getColor("input.invalid",  _INVALID_COLOR)
	_LABEL_COLOR   = tkFlair.getColor("input.label",    _LABEL_COLOR)
	_NORMAL_COLOR  = tkFlair.getColor("input.normal",   _NORMAL_COLOR)
	_SDUM_COLOR    = tkFlair.getColor("input.sdum",     _SDUM_COLOR)
	_TAG_COLOR     = tkFlair.getColor("input.tag",      _TAG_COLOR)
	_VALUE_COLOR   = tkFlair.getColor("input.value",    _VALUE_COLOR)

	_COMMENT_FONT  = tkFlair.getFont("input.comment",   _COMMENT_FONT)
	_FIXED_FONT    = tkFlair.getFont("input.fixed",     _FIXED_FONT)
	_HIDDEN_FONT   = tkFlair.getFont("input.hidden",    _HIDDEN_FONT)
	_LABEL_FONT    = tkFlair.getFont("input.label",     _LABEL_FONT)
	_TAG_FONT      = tkFlair.getFont("input.tag",       _TAG_FONT)
	_VALUE_FONT    = tkFlair.getFont("input.value",     _VALUE_FONT)

	section = tkFlair._INPUT_SECTION
	# Only configurable by hand in the .ini file
	_SHOW_PREPROCESSOR= tkFlair.getBool(section, "showpreprocessor",_SHOW_PREPROCESSOR)
	_INSERT_COMMENT   = tkFlair.getBool(section, "insertcomment",   _INSERT_COMMENT)
	_SHOW_SCALE       = tkFlair.getBool(section, "showscale",       _SHOW_SCALE)
	_BODY_INSERT      = tkFlair.getBool(section, "bodyinsert",	_BODY_INSERT)
	_SHOW_INTER       = tkFlair.getBool(section, "showcard",	_SHOW_INTER)
	_CARDS_TIMEOUT    = float(tkFlair.getInt(section,  "cardstimeout",	int(_CARDS_TIMEOUT)))
