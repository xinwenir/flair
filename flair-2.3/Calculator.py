#!/usr/bin/python
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
# Date:	04-Oct-2012

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import re
try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkFlair
import tkExtra
import Unicode
import sys

import Input
import Manual
import Ribbon
import FlairRibbon

tk = None

ERROR_COLOR    = "Red"
RESULT_COLOR   = "Blue"
COMMENT_COLOR  = "DarkGreen"

_ASSIGN = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$")
_ASSRES = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*) = ([+\-]?\d*\.?\d*[eE]?[+\-]?\d*)$",re.MULTILINE)
_PATANS = re.compile(r"^ans\d+ = .*$",re.MULTILINE)
_PATCOM = re.compile(r"^ *#.*$")
_PATERR = re.compile(r"^Error: .*$",re.MULTILINE)

#===============================================================================
class Calculator(Frame):
	def __init__(self, master, *args):
		Frame.__init__(self, master, *args)
		self.text = Text(self, background="White")
		self.text.pack(expand=YES, fill=BOTH, side=LEFT)
		sb = Scrollbar(self, orient=VERTICAL, command=self.text.yview)
		self.text.config(yscrollcommand=sb.set)
		sb.pack(fill=Y, side=RIGHT)

		self.text.bind("<Return>",		self.returnKey)
		self.text.bind("<KP_Enter>",		self.returnKey)
		self.text.bind("<Shift-Key-Return>",	self.shiftReturn)
#		self.text.bind("<<Paste>>",		self.paste)
		if tk is not None:
			self.text.bind("<Escape>", self.quit)
			self.text.bind("<Control-Key-q>", self.quit)
		self.localDict = Input.Input().localDict

	# ----------------------------------------------------------------------
	def setInput(self, inp):
		self.localDict = inp.localDict

	# ----------------------------------------------------------------------
	# Evaluate all lines
	# ----------------------------------------------------------------------
#	def evalEverything(self, event=None):
#		say("evalEverything")
#		say(self.text.get(1.0,END))

	# ----------------------------------------------------------------------
	def returnKey(self, event=None):
		rowe,cole = list(map(int,self.text.index(END).split(".")))
		row, col  = list(map(int,self.text.index(INSERT).split(".")))

		# Read expression line
		expr = self.text.get("%d.0"%(row),"%d.end"%(row))
		if len(expr)==0 or re.match("^ans\d+ = .*$", expr) or expr[0]=="#":
			self.text.insert(INSERT,"\n")
			if expr and expr[0]=="#":
				self.text.tag_add("comment", "%d.0"%(row), "%d.end"%(row))
				self.text.tag_config("comment", foreground=COMMENT_COLOR)
			self.text.see(INSERT)
			return "break"

		# Get continuation lines
		r = row
		while r>0 and expr[:1]==" ":
			r -= 1
			expr = self.text.get("%d.0"%(r),"%d.end"%(r))+expr

		# Possibly check previous line for continuation
		# Maybe first char could be space
		m = _ASSIGN.match(expr)
		if m:
			variable = m.group(1)
			expr     = m.group(2)
		else:
			# check next line
			nex = self.text.get("%d.0"%(row+1),"%d.end"%(row+1))
			if nex and re.match("^ans\d+ = .*$", nex):
				variable = nex.split("=")[0].strip()
			else:
				variable = self.findUniqueAns()

		# Check for possible vector definitions
		# replace any [ to Vector( and ] to )
		expr = expr.replace("{","Vector(")
		expr = expr.replace("}",")")

		# Check if it is a simple assignment
		try:
			result = float(expr)
			self.localDict[variable] = result
			self.localDict["ans"] = result
			self.text.insert(END,"\n")
			self.text.see(INSERT)
			return "break"
		except:
			pass
		# Evaluate the expression
		try:
			result = eval(expr, Input._globalDict, self.localDict)
			tag    = "result"
			color  = RESULT_COLOR
			self.localDict[variable] = result
			self.localDict["ans"] = result
			result = "%s = %s"%(variable,str(result))
		except:
			result = "Error: "+str(sys.exc_info()[1])
			tag    = "error"
			color  = ERROR_COLOR
		row += 1
		if row == rowe:
			self.text.insert(END,"\n%s\n"%(result))
		else:
			self.text.delete("%d.0"%(row),"%d.end"%(row))
			self.text.insert("%d.0"%(row),"%s"%(result))
			self.text.mark_set(INSERT,"%d.0"%(row+1))
		self.text.tag_add(tag, "%d.0"%(row), "%d.%d"%(row,len(result)))
		self.text.tag_config(tag, foreground=color)
		self.text.see(INSERT)
		return "break"

	# ----------------------------------------------------------------------
	def shiftReturn(self, event=None):
		self.text.insert(INSERT,"\n")
		return "break"

	# ----------------------------------------------------------------------
#	def paste(self, event=None):
#		# Send an evaluate everything message
#		# to be treated after the paste
#		self.event_generate("<<Evaluate>>")

	# ----------------------------------------------------------------------
	def findUniqueAns(self):
		i = 1
		while True:
			ans="ans%d"%(i)
			if ans not in self.localDict: return ans
			i += 1

	# ----------------------------------------------------------------------
	def quit(self, event=None):
		if tk is not None:
			tk.destroy()

	# ----------------------------------------------------------------------
	# Set text and highlight answers and errors
	# ----------------------------------------------------------------------
	def set(self, txt):
		self.localDict.clear()

		self.text.delete("1.0",END)
		self.text.insert("1.0", txt)

		self.text.tag_delete("error")
		self.text.tag_delete("result")

		for m in _PATCOM.finditer(txt):
			self.text.tag_add("comment",
					"1.0 + %d chars"%(m.start()),
					"1.0 + %d chars"%(m.end()))
		for m in _PATANS.finditer(txt):
			self.text.tag_add("result",
					"1.0 + %d chars"%(m.start()),
					"1.0 + %d chars"%(m.end()))
		for m in _PATERR.finditer(txt):
			self.text.tag_add("error",
					"1.0 + %d chars"%(m.start()),
					"1.0 + %d chars"%(m.end()))
		for m in _ASSRES.finditer(txt):
			try:
				self.localDict[m.group(1)] = float(m.group(2))
			except:
				pass

		self.text.tag_config("comment", foreground=COMMENT_COLOR)
		self.text.tag_config("result", foreground=RESULT_COLOR)
		self.text.tag_config("error", foreground=ERROR_COLOR)

	# ----------------------------------------------------------------------
	def get(self):
		return self.text.get("1.0",END).strip()

	# ----------------------------------------------------------------------
	def clear(self):
		self.text.delete("1.0",END)

#===============================================================================
# Button with Label that popups a menu
#===============================================================================
class MenuButton(Ribbon.MenuButton):
	#----------------------------------------------------------------------
	def __init__(self, master, insert, menulist, **kw):
		self.insert = insert
		Ribbon.MenuButton.__init__(self, master, menulist, **kw)

	#----------------------------------------------------------------------
	def createMenuFromList(self, menulist):
		mainmenu = menu = Menu(self, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		for item in menulist:
			if item is None:
				menu.add_separator()
			elif isinstance(item,str):
				menu = Menu(mainmenu)
				mainmenu.add_cascade(label=item, menu=menu)
			else:
				name, acc = item
				cmd = lambda t=name, ins=self.insert: ins(t)
				menu.add_command(label=name,
						accelerator=acc,
						command=cmd)
		return mainmenu

#===============================================================================
# Calculator frame class
#===============================================================================
class CalculatorPage(FlairRibbon.FlairPage):
	"""Scientific calculator for FLUKA"""

	_name_ = "Calculator"
	_icon_ = "calculator"

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Plot List ===========
		group = Ribbon.LabelGroup(self.ribbon, "Functions")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		menulist = [
			("a",	"Fine constant"),
			("amu",	"Atomic Mass Unit"),
			("amuC12","Atomic Mass Unit C12"),
			("amugr","Atomic Mass Unit in gr"),
			("c",	"Speed of light"),
			("fwhm","Full Width at Half Max"),
			("Na",	"Avogadro number"),
			("pi",	"3.14159..."),
			("qe",	"Electron charge"),
			("re",	"Electron radius") ]

		col,row = 1,0
		b = MenuButton(group.frame, self.insert, menulist,
				image=tkFlair.icons["constant"],
				text="Constants "+Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Constants")

		# ---
		menulist = [
			"Angle",
			(	"deg"	, "pi/180"),
			(	"rad"	, "1"),
			(	"mrad"	, "0.001"),

			"Distance",
			(	"nm",	"0.1E-6"),
			(	"um",	"0.1E-3"),
			(	"mm",	"0.1"),
			(	"cm",	"1"),
			(	"m",	"100"),
			(	"km",	"100E3"),
			None,
			(	"inch"	, "2.54"),
			(	"feet"	, "30.48"),
			(	"ft"	, "30.48"),
			(	"mile"	, "160934.4"),
			(	"mi"	, "160934.4"),

			"Energy",
			(	"eV"	, "1E-9"),
			(	"keV"	, "1E-6"),
			(	"MeV"	, "1e-3"),
			(	"GeV"	, "1"),
			(	"TeV"	, "1000"),
			(	"PeV"	, "1E6"),
			(	"J"     , "eV/qe"),

			"Time",
			(	"ns"	, "1.0E-9"),
			(	"us"	, "1.0E-6"),
			(	"ms"	, "0.001"),
			(	"s"	, "1"),
			(	"min"	, "60"),
			(	"hour"	, "3600"),
			(	"day"	, "24*hour"),
			(	"week"	, "7*day"),
			(	"month"	, "365.25/12*day"),
			(	"year"	, "365.25*day") ]

		col,row = 1,1
		b = MenuButton(group.frame, self.insert, menulist,
				image=tkFlair.icons["units"],
				text="Units "+Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Conversion Units")

		# ---
		menulist = [
			"Vector",
			("{x,y,z}",	"vector"),
			(".x()",	"x component"),
			(".y()",	"y component"),
			(".z()",	"z component"),
			("*",		"dot product"),
			(".dot(v)",	"dot product"),
			("^",		"cross product"),
			(".cross(v)",	"cross product"),
			(".length()",	"length"),
			("abs(v)",	"length of vector"),
			(".unit()",	"unit vector"),
			(".phi()",	"azimuthal angle"),
			(".theta()",	"polar angle"),
#			"Matrix",
#
			"List",
			("[x,y,z]",	"list or array")]

		col,row = 1,2
		b = MenuButton(group.frame, self.insert, menulist,
				image=tkFlair.icons["matrix"],
				text="Vector "+Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Vector functions")

		# ---
		col,row = 2,0
		menulist = [
			("Mp",	"Proton mass"),
			("Mn",	"Neutron mass"),
			("Me",	"Electron mass") ]

		b = MenuButton(group.frame, self.insert, menulist,
				image=tkFlair.icons["particle"],
				text="Particles "+Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Particle masses")

		# ---
		col,row = 2,1
		menulist = [
			"Physics",
			("T2p(E,m)"	, "Energy to Momentum"),
			("p2T(p,m)"	, "Momentum to Energy"),
			("dT2dp(E,dE,m)", "Energy spread in momentum"),
			("T2g(E,m)"	, "Energy to gamma"),
			("g2b(g)"	, "gamma to beta"),
			("X0(Z,A)"      , "Radiation length in g/cm2"),

			"Logarithmic",
			("exp(x)",	"e**x"),
			("expm1(x)",	"e**x-1"),
			("log(x)",	"Natural logarithm of x (base e)"),
			("log(x,b)",	"logarithm of x to the given base"),
			("log1p(x)",	"Natural logarithm of 1+x (base e)"),
			("log10(x)",	"Base-10 logarithm of x"),
			("pow(x,y)",	"x raised to the power y"),
			("sqrt(x)",	"Return the square root of x"),

			"Trigonometric",
			("acos(x)",	"Arc cosine to radians"),
			("asin(x)",	"Arc sine to radians"),
			("atan(x)",	"Arc tangent to radians"),
			("atan2(y,x)",	"Atan(y/x) to radians"),
			("cos(x)",	"Cosine from radians"),
			("hypot(x,y)",	"Euclidean norm, sqrt(x*x + y*y)"),
			("sin(x)",	"Sine from radians"),
			("tan(x)",	"Tangent from radians"),
			None,
			("acosd(x)"	, "acos to degrees"),
			("asind(x)"	, "asin to degrees"),
			("atand(x)"	, "atan to degrees"),
			("cosd(x)"	, "cos from degrees"),
			("sind(x)"	, "sin from degrees"),
			("tand(x)"	, "tan from degrees"),
			("omega(x)"	, "Solid angle from radians"),

			"Hyperbolic",
			("acosh(x)",	"Inverse hyperbolic cosine"),
			("asinh(x)",	"Inverse hyperbolic sine"),
			("atanh(x)",	"Inverse hyperbolic tangent"),
			("cosh(x)",	"Hyperbolic cosine"),
			("sinh(x)",	"Hyperbolic sine"),
			("tanh(x)",	"Hyperbolic tangent"),
			("omegad(x)",	"Solid angle from degrees"),

			"Special",
			("erf(x)",	"Error function"),
			("erfc(x)",	"Complementary error function"),
			("gamma(x)",	"Gamma function"),
			("lgamma(x)",	"Natural logarithm of the absolute Gamma function"),
			]
		b = MenuButton(group.frame, self.insert, menulist,
				image=tkFlair.icons["function"],
				text="Functions "+Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Physics and mathematical functions")

		# ---
		col,row = 2,2
		menulist = [
			("w(n)",	"n'th what of current card"),
			("b(n,w)",	"w'th what of body named n"),
			("C(t,n,w)",	"n'th what of card tag t with sdum n")
			]
		b = MenuButton(group.frame, self.insert, menulist,
				image=tkFlair.icons["input"],
				text="Card "+Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Card cross reference")

		# ---
		# ========== Plot List ===========
		group = Ribbon.LabelGroup(self.ribbon, "Action")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		group.frame.grid_rowconfigure(0, weight=1)
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clean32"],
				text="Clear",
				compound=TOP,
				command=self.clear,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Clear calculator")

	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self, "Calculator")

		# --- Frame with listbox ---
		self.calculator = Calculator(self.frame)
		self.calculator.text.config(height=5)
		self.calculator.pack(side=TOP, expand=YES, fill=BOTH)

	# ----------------------------------------------------------------------
	# Refresh calculator
	# ----------------------------------------------------------------------
	def refresh(self):
		self.calculator.set(self.project.calculator)

	# ----------------------------------------------------------------------
	# Update fields
	# ----------------------------------------------------------------------
	def get(self, event=None):
		if self.page is None or self.project is None: return

		txt= self.calculator.get()
		if txt != self.project.calculator:
			self.project.calculator = txt
			self.setModified(True)

	#----------------------------------------------------------------------
	def clear(self):
		self.calculator.clear()

	# ----------------------------------------------------------------------
	def insert(self, txt):
		self.calculator.text.insert(INSERT, txt)

#-------------------------------------------------------------------------------
def quit(event):
	tk.destroy()

if __name__ == "__main__":
	Input.init()
	tk = Tk()
	tkFlair.openIni()
#	tkFlair.loadIcons()

	tk.title("Flair calculator")
	tk.iconbitmap("@%s/icons/calculator.xbm"%(tkFlair.prgDir))
#	tk.iconwindow(icon)
	calc = Calculator(tk)
	calc.pack(expand=YES, fill=BOTH)
	tk.bind("<Escape>", quit)
	tk.bind("<Control-Key-q>", quit)
	calc.text.focus_set()
	tk.mainloop()
#	tkFlair.delIcons()
