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
# Date:	24-Sep-2006

# This file implements the Element and PeriodicTable
# classes.
# The first displayes a dialog with the element properties.
# The creates the Medeleyv table

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

try:
	from tkinter import *
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	import tkinter.messagebox as messagebox
	from legacy import cmp
from legacy import cmp_to_key

try:
	import configparser
except ImportError:
	import configparser as ConfigParser

import tkExtra
import tkFlair
import Unicode

import Ribbon
import FlairRibbon

from legacy import cmp

#-------------------------------------------------------------------------------
def _cmpA(x, y):
	x = x[0]
	y = y[0]
	try:
		ix = int(x)
		nx = ""
	except:
		ix = int(x[:-1])
		nx = x[-1:]
	try:
		iy = int(y)
		ny = ""
	except:
		iy = int(y[:-1])
		ny = y[-1:]

	c = cmp(ix, iy)
	if c == 0:
		return cmp(nx, ny)
	return c

#===============================================================================
# Database
#-------------------------------------------------------------------------------
_elements = []
_footnote = {}

#-------------------------------------------------------------------------------
# Static-local Variables
#-------------------------------------------------------------------------------
_PTBL_ROWS = 12
_PTBL_COLS = 20

#-------------------------------------------------------------------------------
# Labels of periodic table
#-------------------------------------------------------------------------------
_ptableLbl = [	"1",			# -1
		"2",			# -2
		"3",			# -3
		"4",			# -4
		"5",			# -5
		"6",			# -6
		"7",			# -7
		"8",			# -8
		"9",			# -9
		"10",			# -10
		"11",			# -11
		"12",			# -12
		"13",			# -13
		"14",			# -14
		"15",			# -15
		"16",			# -16
		"17",			# -17
		"18",			# -18
		"Group",		# -19
		"Period",		# -20
		"*",			# -21
		"**",			# -22
		"*\nLanthanides",	# -23
		"**\nActinides"]	# -24

#-------------------------------------------------------------------------------
# Positive is given the Z of the elements
# Zero is empty space
# Negative is a string indexed from the above array ptable[abs(num)]
#-------------------------------------------------------------------------------
_ptable = [
   [-19, -1, -2,  0, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,-16,-17,-18],
   [-20,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
   [ -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2],
   [ -2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  6,  7,  8,  9, 10],
   [ -3, 11, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 13, 14, 15, 16, 17, 18],
   [ -4, 19, 20,  0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
   [ -5, 37, 38,  0, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54],
   [ -6, 55, 56,-21, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86],
   [ -7, 87, 88,-22,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118],
   [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
   [-23,  0,  0,-21, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70,  0,  0],
   [-24,  0,  0,-22, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,  0,  0] ]


#-------------------------------------------------------------------------------
#            A:0      B:1        C:2        D:3        E:4         F:5       G:6
#-------------------------------------------------------------------------------
#               w         A         B           K          O
_bcolor  = [ "White", "#CCCCCC", "#CCCCFF", "#FFFFCC", "#FFCCCC",
#               G         M         C          P           T
	    "#CCFFCC", "#FFCC00", "#CCFFFF", "#FFCCFF", "#AADDFF" ]
#0 1 2 3 4 5 6 7 8 9
(w,A,B,K,O,G,M,C,P,T) = list(range(10))
_ecolor = [	[ A,A,A,w,A, A,A,A,A,A, A,A,A,A,A, A,A,A,A,A ],
		[ A,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,w ],
		[ A,G,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,O ],
		[ A,M,C,w,w, w,w,w,w,w, w,w,w,w,G, G,G,G,K,O ],
		[ A,M,C,w,w, w,w,w,w,w, w,w,w,w,P, G,G,G,K,O ],
		[ A,M,C,w,B, B,B,B,B,B, B,B,B,B,P, P,G,G,K,O ],
		[ A,M,C,w,B, B,B,B,B,B, B,B,B,B,P, P,P,G,K,O ],
		[ A,M,C,w,B, B,B,B,B,B, B,B,B,B,P, P,P,P,K,O ],
		[ A,M,C,w,B, B,B,B,B,B, B,B,B,B,A, A,A,A,A,A ],
		[ w,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,w ],
		[ A,A,A,w,T, T,T,T,T,T, T,T,T,T,T, T,T,T,w,w ],
		[ A,A,A,w,T, T,T,T,T,T, T,T,T,T,T, T,T,T,w,w ] ]

#             w         B      R        G
_fcolor = [ "Black", "Blue", "Red", "#707070" ]
(w,B,R,G) = list(range(4))
_tcolor = [	[ w,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,w ],
		[ w,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,w ],
		[ w,B,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,B ],
		[ w,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,B,B,B,B ],
		[ w,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,B,B ],
		[ w,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,R,B ],
		[ w,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,B ],
		[ w,R,w,w,w, w,w,w,w,w, w,w,w,R,w, w,w,w,w,B ],
		[ w,R,w,w,w, G,G,G,G,G, G,G,G,G,G, G,G,G,G,G ],
		[ w,w,w,w,w, w,w,w,w,w, w,w,w,w,w, w,w,w,w,w ],
		[ w,w,w,w,w, w,w,w,G,w, w,w,w,w,w, w,w,w,w,w ],
		[ w,w,w,w,w, w,w,w,G,G, G,G,G,G,G, G,G,G,w,w ] ]

#===============================================================================
# Element Class
#===============================================================================
class Element:
	def __init__(self, symbol, Z):
		self.Z          = Z
		self.symbol     = symbol
		self.name       = ""
		self.AmassStr   = ""
		self.densityStr = ""
		self.melting    = ""
		self.boiling    = ""
		self.oxidation  = ""
		self.isotopes   = []

	# ----------------------------------------------------------------------
	def Amass(self):
		return float(self.AmassStr.split()[0])

	# ----------------------------------------------------------------------
	def density(self):
		return float(self.densityStr.split()[0])

	# ----------------------------------------------------------------------
	def densityNote(self):
		d = self.densityStr.split()
		if len(d)>0:
			return d[-1]
		return None

#===============================================================================
# Element Dialog
#===============================================================================
class ElementDialog:
	def __init__(self, master, element):
		if master is None:
			self.toplevel = Toplevel(master)
		else:
			self.toplevel = master

		frame = self.toplevel
		self.element = element

		self.symbol = Label(frame, font="Helvetica -72 bold",
				relief=GROOVE, padx=30)
		self.symbol.grid(row=0, rowspan=5, column=0, columnspan=3, sticky=NSEW)

		b = Button(frame, text=Unicode.WHITE_LEFT_POINTING_TRIANGLE, padx=2, pady=0,
				command=self.prevZ)
		b.grid(row=5, column=0, sticky=W)
		self.name = Label(frame, foreground="DarkBlue",
				width=13, relief=GROOVE)
		self.name.grid(row=5, column=1, sticky=EW)
		b = Button(frame, text=Unicode.WHITE_RIGHT_POINTING_TRIANGLE, padx=2, pady=0,
				command=self.nextZ)
		b.grid(row=5, column=2, sticky=E)

		# ------
		row = 0
		col = 3
		l = Label(frame, text="Z:", justify=RIGHT, anchor=E)
		l.grid(row=row, column=col, sticky=EW)

		row += 1
		l = Label(frame, text="Atomic Weight:", justify=RIGHT, anchor=E)
		l.grid(row=row, column=col, sticky=EW)

		row += 1
		l = Label(frame, text="Density:", justify=RIGHT, anchor=E)
		l.grid(row=row, column=col, sticky=EW)

		row += 1
		l = Label(frame, text="Melting:", justify=RIGHT, anchor=E)
		l.grid(row=row, column=col, sticky=EW)

		row += 1
		l = Label(frame, text="Boiling:", justify=RIGHT, anchor=E)
		l.grid(row=row, column=col, sticky=EW)

		row += 1
		l = Label(frame, text="Oxidation:", justify=RIGHT, anchor=E)
		l.grid(row=row, column=col, sticky=EW)

		# ------
		row = 0
		col += 1
		self.Z = Label(frame, justify=LEFT,
			background="White", foreground="DarkBlue",
			anchor=W)
		self.Z.grid(row=row, column=col, sticky=EW, pady=1)

		row += 1
		self.Amass = Label(frame, justify=LEFT,
			background="White", foreground="DarkBlue",
			anchor=W)
		self.Amass.grid(row=row, column=col, sticky=EW, pady=1)

		row += 1
		self.density = Label(frame, justify=LEFT,
			background="White", foreground="DarkBlue",
			anchor=W)
		self.density.grid(row=row, column=col, sticky=EW, pady=1)

		row += 1
		self.melting = Label(frame, justify=LEFT,
			background="White", foreground="DarkBlue",
			anchor=W)
		self.melting.grid(row=row, column=col, sticky=EW, pady=1)

		row += 1
		self.boiling = Label(frame, justify=LEFT,
			background="White", foreground="DarkBlue",
			anchor=W)
		self.boiling.grid(row=row, column=col, sticky=EW, pady=1)

		row += 1
		self.oxidation = Label(frame, justify=LEFT,
			background="White", foreground="DarkBlue",
			anchor=W)
		self.oxidation.grid(row=row, column=col, sticky=EW, pady=1)

		# ------
		row += 1
		self.isotopeList = tkExtra.MultiListbox(frame,
				(('A',                         5, None),
				 ('J\u03c0',                  5, None),
				 ('\u0394 (MeV)',             8, None),
				 ('T1/2, \u0393, Abundance', 15, None),
				 ('Decay Mode',                8, None)))
		self.isotopeList.grid(row=row, column=0, columnspan=5, sticky=NSEW)
		self.isotopeList.bindLeftRight()

		row += 1
		self.footnote = Label(frame, justify=LEFT, anchor=W,
				relief=GROOVE, background="White")
		self.footnote.grid(row=row, column=0, columnspan=5, sticky=NSEW)

		frame.grid_rowconfigure(row-1, weight=3)
		frame.grid_columnconfigure(3, weight=1)
		frame.grid_columnconfigure(4, weight=3)

		self.fill()

		if master is None:
			self.toplevel.bind("<Escape>", self.close)
			self.toplevel.bind("<Control-Key-q>", self.close)
			self.toplevel.protocol("WM_DELETE_WINDOW", self.close)
			#self.toplevel.transient(master)
			#self.toplevel.deiconify()
			#self.toplevel.wait_visibility()
			self.toplevel.focus_set()
			#self.toplevel.grab_set()
			self.toplevel.wait_window()

	# ----------------------------------------------------------------------
	def fill(self):
		element = self.element

		try:
			self.toplevel.title("Element: %d-%s %s" \
				% (element.Z, element.symbol, element.name))
		except:
			pass

		# Find color
		color = "White"
		for j in range(_PTBL_ROWS):
			try:
				i = _ptable[j].index(element.Z)
				fcol = _fcolor[_tcolor[j][i]]
				bcol = _bcolor[_ecolor[j][i]]
				break
			except:
				pass

		self.symbol.config(text=element.symbol,
				background=bcol, foreground=fcol)
		self.name.config(text=element.name)

		self.Z.config(text=element.Z)
		self.Amass.config(text=element.AmassStr)
		self.density.config(text=element.densityStr)
		self.melting.config(text=element.melting)
		self.boiling.config(text=element.boiling)
		self.oxidation.config(text=element.oxidation)

		# Fill isotope List
		self.isotopeList.delete(0,END)
		for (a,info) in element.isotopes:
			jpi = info[0:13].strip()
			d   = info[13:23].strip()
			t12 = info[23:41].strip()
			dec = info[41:].strip()
			self.isotopeList.insert(END, (a,jpi,d,t12,dec))

		# Show footnote if needed
		f = element.densityNote()
		if f:
			try:
				note = _footnote[f]
				self.footnote.config(text="%s) %s" % (f,note))
			except:
				pass

	# ----------------------------------------------------------------------
	def prevZ(self, event=None):
		if self.element.Z>1:
			self.element = element(self.element.Z-1)
			self.fill()

	# ----------------------------------------------------------------------
	def nextZ(self, event=None):
		if self.element.Z<118:
			self.element = element(self.element.Z+1)
			self.fill()

	# ----------------------------------------------------------------------
	# Close window
	# ----------------------------------------------------------------------
	def close(self, event=None):
		self.toplevel.destroy()

#===============================================================================
# Periodic Table class
#===============================================================================
class PeriodicTable(Frame):
	def __init__(self, master, isotopeFrame=None):
		global _elements
		self.selected = None
		self.element  = None
		self.isotopeFrame = isotopeFrame

		Frame.__init__(self, master, class_="PeriodicTable", background="White")

		for j in range(_PTBL_ROWS):
			for i in range(_PTBL_COLS):
				val = _ptable[j][i]
				fcol = _fcolor[_tcolor[j][i]]
				bcol = _bcolor[_ecolor[j][i]]

				if val>0:
					elem   = _elements[val-1]
					l = Label(self,
						text="%d\n%s" % (val, elem.symbol),
						foreground=fcol,
						background=bcol,
						justify=CENTER,
						borderwidth=2,
						width = 3,
						padx = 1,
						pady = 0,
						relief=FLAT)
					l.bind('<Button-1>', self.click)
				elif val<0:
					l = Label(self,
						pady=0,
						text=_ptableLbl[-val-1])
				else:
					continue
				l.grid(row=j, column=i, sticky=NSEW, padx=1, pady=1)

		# ---
		row, col, span = _PTBL_ROWS, 2, 4
		Label(self, background="White").grid(row=row, column=col)

		row += 1
		Label(self, text="H - gas",
			foreground=_fcolor[1], background="White").grid( \
			row=row, column=col, columnspan=span, sticky=W)

		col += span
		Label(self, text="Li - solid",
			foreground=_fcolor[0], background="White").grid( \
			row=row, column=col, columnspan=span, sticky=W)

		col += span
		Label(self, text="Br - liquid",
			foreground=_fcolor[2], background="White").grid( \
			row=row, column=col, columnspan=span, sticky=W)

		col += span
		Label(self, text="Tc - synthetic",
			foreground=_fcolor[3], background="White").grid( \
			row=row, column=col, columnspan=span, sticky=W)

		# ---
		row += 1
		col  = 1
		Message(self, text=" ", background=_bcolor[5],
			relief=GROOVE).grid(row=row, column=col, pady=3)
		col += 1
		span = 3
		Label(self, text="Non-Metals", background="White").grid(
			row=row, column=col, columnspan=span, sticky=W)

		col += span
		Message(self, text=" ", background=_bcolor[2],
			relief=GROOVE).grid(row=row, column=col)
		col += 1
		Label(self, text="Transition Metals", background="White").grid(
			row=row, column=col, columnspan=span, sticky=W)

		col += span
		Message(self, text=" ", background=_bcolor[9],
			relief=GROOVE).grid(row=row, column=col)
		col += 1
		Label(self, text="Rare Earth Metals", background="White").grid(
			row=row, column=col, columnspan=span, sticky=W)

		col += span
		Message(self, text=" ", background=_bcolor[3],
			relief=GROOVE).grid(row=row, column=col)
		col += 1
		Label(self, text="Halogens", background="White").grid( \
			row=row, column=col, columnspan=span, sticky=W)

		# ---
		row += 1
		col  = 1
		Message(self, text=" ", background=_bcolor[6],
			relief=GROOVE).grid(row=row, column=col)
		col += 1
		span = 3
		Label(self, text="Alkali Metals", background="White").grid( \
			row=row, column=col, columnspan=span, sticky=W)

		col += span
		Message(self, text=" ", background=_bcolor[7],
			relief=GROOVE).grid(row=row, column=col)
		col += 1
		Label(self, text="Alkali Earth Metals", background="White").grid( \
			row=row, column=col, columnspan=span, sticky=W)

		col += span
		Message(self, text=" ", background=_bcolor[8],
			relief=GROOVE).grid(row=row, column=col)
		col += 1
		Label(self, text="Other Metals", background="White").grid( \
			row=row, column=col, columnspan=span, sticky=W)

		col += span
		Message(self, text=" ", background=_bcolor[4],
			relief=GROOVE).grid(row=row, column=col)
		col += 1
		Label(self, text="Inert Elements", background="White").grid( \
			row=row, column=col, columnspan=span, sticky=W)

	# ----------------------------------------------------------------------
	def click(self, event):
		if self.selected:
			#self.selected.config(borderwidth=1)
			self.selected.config(relief=FLAT)
		z = self.buttonZ(event.widget)
		self.selected = event.widget
		self.selected.config(relief=RAISED) #borderwidth=3)
		dlg = ElementDialog(self.isotopeFrame, element(z))
		self.element = dlg.element

	# ----------------------------------------------------------------------
	def buttonZ(self, widget):
		return int(widget["text"].split()[0])

#===============================================================================
#class PeriodicTableSelect(PeriodicTable):
#	def __init__(self, master):
#		PeriodicTable.__init__(self, master)
#
#	# ----------------------------------------------------------------------
#	def click(self, event):
#		self.selected = self.buttonZ(event.widget)
#		self.destroy()

#===============================================================================
# Periodic Elements List
#===============================================================================
class ElementsList(Frame):
	def __init__(self, master, isotopeFrame=None):
		global _elements

		self.element = None
		self.isotopeFrame = isotopeFrame
		Frame.__init__(self, master, class_="ElementsList")

		self.listbox = tkExtra.MultiListbox(self,
				(("Z"        ,  3, None),
				 ("El"       ,  3, None),
				 ("Name"     , 12, None),
				 ("Weight"   , 14, None),
				 ("Density"  , 10, None),
				 ("Melting"  ,  8, None),
				 ("Boiling"  ,  8, None),
				 ("Oxidation",  8, None)),
				 background="White",
				 height=20)
		self.listbox.pack(expand=YES, fill=BOTH)
		self.listbox.bindLeftRight()

		for elem in _elements:
			if elem is None: continue
			self.listbox.insert(END,
				(elem.Z, elem.symbol, elem.name, elem.AmassStr,
				elem.densityStr, elem.melting, elem.boiling, elem.oxidation))
		self.listbox.bind('<<ListboxSelect>>', self.show)

	# --------------------------------------------------------------------
	def show(self, event=None):
		global _elements

		sel = [ int(x) for x in self.listbox.curselection() ]
		if len(sel) == 1:
			isotope = self.listbox.get(sel[0])
			elem = _elements[isotope[0]-1]
			if elem is None: return
			dlg = ElementDialog(self.isotopeFrame, elem)
			self.element = dlg.element

#-------------------------------------------------------------------------------
def element(Z):
	return _elements[Z-1]

#-------------------------------------------------------------------------------
def elementBySymbol(symbol):
	symbol = symbol.upper()
	for element in _elements:
		if element.symbol.upper() == symbol:
			return element
	return None

#-------------------------------------------------------------------------------
def elementByName(name):
	name = name.upper()
	for element in _elements:
		if element.name.upper() == name:
			return element
	return None

#-------------------------------------------------------------------------------
def init(filename):
	global _elements, _footnote
	c = configparser.RawConfigParser()
	c.read(filename)

	# Read Elements
	_elements = [None] * len(c.sections())
	for symb in c.sections():
		if symb[0].islower(): continue
		Z         = c.getint(symb, "Z")
		elem      = Element(symb, Z)
		elem.name = c.get(symb, "name")
		try:	elem.AmassStr   = c.get(symb, "Amass")
		except: pass
		try:	elem.densityStr = c.get(symb, "density")
		except: pass
		try:	elem.melting = c.get(symb, "melting")
		except: pass
		try:	elem.boiling = c.get(symb, "boiling")
		except: pass
		try:	elem.oxidation = c.get(symb, "oxidation")
		except: pass

		elem.isotopes = []
		for item,value in c.items(symb):
			if item.startswith("isotope."):
				iso = item.split(".")[1]
				if "\n" in value:
					for line in value.splitlines():
						elem.isotopes.append((iso, line[1:]))
				else:
					elem.isotopes.append((iso, value[1:]))
		elem.isotopes.sort(key=cmp_to_key(_cmpA))
		_elements[Z-1] = elem

	# Remove all None parts
	try:
		while True: _elements.remove(None)
	except ValueError:
		pass

	for (i,n) in c.items("footnotes"):
		_footnote[i] = n #.replace("\n"," ")

	del c

#===============================================================================
# Periodic Table Frame class
#===============================================================================
class PeriodicTablePage(FlairRibbon.FlairPage):

	_name_ = "Elements"
	_icon_ = "toi"

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Input ===========
		group = Ribbon.LabelGroup(self.ribbon, "View")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		kw = {
			"compound" : TOP,
			"variable" : self.tabVar,
			"command"  : self.change,
			"padx"     : 5
		}

		# ---
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["toi32"],
				text="Table",
				value=0,
				**kw)
		b.pack(side=LEFT, fill=Y, expand=YES)
		tkExtra.Balloon.set(b, "Show as periodic table")

		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["list32"],
				text="List",
				value=1,
				**kw)
		tkExtra.Balloon.set(b, "Show as list")
		b.pack(side=LEFT, fill=Y, expand=YES)

	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self, False)

		self.frame.config(background="White")
		self.isotope = Frame(self.frame)

		self.periodic_frame = Frame(self.frame, background="White")

		self.periodic_table = PeriodicTable(self.periodic_frame, self.isotope)
		self.periodic_table.pack(expand=YES, fill=X)

		self.periodic_list = ElementsList(self.frame, self.isotope)

		self.tabVar = IntVar()
		self.tabVar.set(0)

		self.change()
		self.isotope.pack(side=RIGHT, expand=YES, fill=BOTH)

	#----------------------------------------------------------------------
	def change(self):
		if self.tabVar.get() == 0:
			self.periodic_frame.pack(side=LEFT, anchor=NW)
			self.periodic_list.pack_forget()
		else:
			self.periodic_frame.pack_forget()
			self.periodic_list.pack(side=LEFT, anchor=NW, fill=Y)

#===============================================================================
# TOIFrame
#===============================================================================
class TOIFrame(Frame):
	def __init__(self, master, isotopeFrame=None):
		Frame.__init__(self, master, class_="TOIFrame")

		self.tabPage = tkExtra.TabPageSet(self, pageNames=['Table','List'])
		self.tabPage.pack(expand=YES, fill=BOTH)

		frame = self.tabPage['Table']
		self.periodic_table = PeriodicTable(frame, isotopeFrame)
		self.periodic_table.pack(expand=YES, fill=BOTH)

		frame = self.tabPage['List']
		self.periodic_list = ElementsList(frame, isotopeFrame)
		self.periodic_list.pack(expand=YES, fill=BOTH)

		self.tabPage.changePage()

	# ----------------------------------------------------------------------
	def getLast(self):
		actPage = self.tabPage.getActivePage()
		if actPage == 'Table':
			return self.periodic_table.element
		else:
			return self.periodic_list.element

#===============================================================================
def show(master, swallow=False, select=False):
	if swallow:
		toplevel = master
	else:
		toplevel = Toplevel(master)
		toplevel.bind("<Escape>", lambda e:toplevel.destroy())
		toplevel.bind("<Control-Key-q>", lambda e:toplevel.destroy())

	toplevel.title("Table of Elements")	#Periodic Table")
	toi = TOIFrame(toplevel)
	toi.pack(side=TOP, expand=YES, fill=BOTH)
	toplevel.deiconify()
	toplevel.wm_resizable(False,False)
	toplevel.update_idletasks()
	toplevel.focus_set()
	toplevel.lift()
	if select:
		toplevel.grab_set()
		toplevel.wait_window()
		return toi.getLast()
	else:
		toplevel.wait_visibility()
		return None

#===============================================================================
if __name__ == "__main__":
	import sys, os
	init(os.path.join(os.path.dirname(sys.argv[0]),"db/isotopes.ini"))
	root = Tk()
	root.deiconify()
	show(root,True)
	root.mainloop()
