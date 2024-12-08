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
# Date:	05-Sep-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import tkinter as tk
import os
import re
import time
import math

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import bmath
import Input
import Ribbon
import Manual
import Project
import FlairRibbon

import tkFlair
import tkExtra
import Unicode
import tkDialogs
import bFileDialog

import Plot
import Gnuplot
import PlotEngine

_GEO_NO    = "-No-"
_GEO_AUTO  = "-Auto-"

_AXES = ["x","y","z","cb","x2","y2"]

_NORM_BALLOON = "Normalization factor (eg. 1e11) or a constant expression (eg. 7e12*1.6e-10) or an expression using variable 'x' for the value of the bin (eg. 7e12*1.6e-10*x/4.2)"

#-------------------------------------------------------------------------------
def setFloat(entry, value):
	value = float(value)
	if abs(value) < 1e-15: value=0.0
	s = bmath.format(value,10,useD=False)

	# remember state
	state = entry["state"]
	entry["state"] = NORMAL

	# set value
	entry.delete(0,END)
	entry.insert(0,s)

	# restore state
	entry["state"] = state

#===============================================================================
# Base of Plot Frames
#===============================================================================
class BasePlotFrame(Frame):
	# ----------------------------------------------------------------------
	def __init__(self, master, flair):
		Frame.__init__(self, master)
		self.flair    = flair
		self.project  = flair.project
		self.widget   = {}
		self.variable = {}
		self.plot     = None

	# ----------------------------------------------------------------------
	def type(self): raise

	# ----------------------------------------------------------------------
	# set fields
	# ----------------------------------------------------------------------
	def set(self, plot):
		if plot is None: return
		self.plot = plot
		# Set entries
		for name, widget in list(self.widget.items()):
			if name[0]=="@": continue
			if isinstance(widget,tuple): widget = widget[0]
			var = self.variable.get(name)
			if var is not None:
				value = plot[name]
				if value in ("True","False"): value = (value=="True")
				if isinstance(var,IntVar) or isinstance(var,BooleanVar):
					try:
						var.set(int(value))
					except ValueError:
						var.set(0)
				else:
					var.set(value)

			elif isinstance(widget, Entry):
				widget.delete(0, END)
				widget.insert(0, str(plot[name]))

			elif isinstance(widget, Text):
				widget.delete("1.0", END)
				widget.insert("1.0", str(plot[name]))

			elif isinstance(widget, tkExtra.Combobox):
				value = plot[name]
				if value:
					widget.set(value)
				else:
					widget.select(0)

	# ----------------------------------------------------------------------
	# get fields
	# ----------------------------------------------------------------------
	def get(self, plot):
		if plot is None: return
		# Get values from widgets
		for name,widget in list(self.widget.items()):
			if name[0]=="@": continue
			if isinstance(widget,tuple):
				for w in widget:
					if w.cget("state") == NORMAL:
						break
				else:
					continue
				widget = widget[0]
			elif widget.cget("state") == DISABLED:
				continue

			var = self.variable.get(name)

			if var is not None:
				try:
					value = var.get()
				except ValueError:
					value = ""

			elif isinstance(widget, Entry):
				value = str(widget.get())

			elif isinstance(widget, Text):
				value = widget.get(1.0, tk.END).strip()

			elif isinstance(widget, tkExtra.Combobox):
				value = widget.get()

			else:
				value = ""

			if value=="" and name in plot.var:
				del plot.var[name]
			else:
				plot.var[name] = value

	# ----------------------------------------------------------------------
	# update fields after a plot
	# ----------------------------------------------------------------------
	def afterPlot(self):
		pass

	# ----------------------------------------------------------------------
	# setState in all widgets
	# ----------------------------------------------------------------------
	def setState(self, state):
		#for name,widget in self.widget.items():
		#	say(">>>",name)
		for widget in list(self.widget.values()):
			if isinstance(widget,tuple):
				for w in widget: w.config(state=state)
			else:
				widget.config(state=state)

	# ----------------------------------------------------------------------
	def toggleField(self, event):
		if event.widget["state"] == NORMAL:
			state = DISABLED
		else:
			state = NORMAL

		if isinstance(event.widget.master, tkExtra.Combobox):
			event.widget.master.config(state=state)
		else:
			event.widget.config(state=state)

		# Search linked list of widgets
		for widget in list(self.widget.values()):
			if isinstance(widget,tuple) and event.widget in widget:
				for w in widget:
					w.config(state=state)

	#-----------------------------------------------------------------------
	def bindToggle(self):
		for name,widget in list(self.widget.items()):
			if isinstance(widget,tuple):
				for w in widget: w.bind("<Button-3>", self.toggleField)
			elif isinstance(widget,tkExtra.Combobox):
				widget.bindWidgets("<Button-3>", self.toggleField)
			else:
				widget.bind("<Button-3>", self.toggleField)

	# ----------------------------------------------------------------------
	def styleFrame(self, frame):
		row,col = 0,0
		Label(frame, text="Type:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_with = tkExtra.Combobox(frame, width=10)
		self.plot_with.fill(Plot.STYLE)
		self.plot_with.set(Plot.STYLE[Plot.STYLE_DEF])
		self.plot_with.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.plot_with, "Plotting style")
		self.widget["with.0"] = self.plot_with

		# ===
		row += 1
		col = 0
		Label(frame, text="Color:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_lc = tkExtra.Combobox(frame, width=15)
		self.plot_lc.fill(PlotEngine._COLOR_NAMES)
		self.plot_lc.grid(row=row,column=col,sticky=EW)
		tkExtra.Balloon.set(self.plot_lc, "Line color if gnuplot>=4.2")
		self.widget["lc.0"] = self.plot_lc

		# ---
		col += 1
		Label(frame, text=" Line width:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_lw = IntVar()
		self.plot_lw.set(1)
		self.plot_lw_spin = Spinbox(frame, text=self.plot_lw, from_=1, to=6,
				background="White",width=3)
		self.plot_lw_spin.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(self.plot_lw_spin, "Line width")
		self.widget["lw.0"] = self.plot_lw_spin
		self.variable["lw.0"] = self.plot_lw

		# ===
		row, col = row+1, 0
		Label(frame, text=" Point type:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_pt = tkExtra.Combobox(frame, width=15)
		self.plot_pt.fill(PlotEngine.POINTS)
		self.plot_pt.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(self.plot_pt, "Point type to use")
		self.widget["@pt.0"] = self.plot_pt

		# ---
		col += 1
		Label(frame, text="Point size:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_ps = IntVar()
		self.plot_ps.set(1)
		self.plot_ps_spin = Spinbox(frame, text=self.plot_ps, from_=1, to=6,
				background="White",width=3)
		self.plot_ps_spin.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(self.plot_ps_spin, "Point size")
		self.widget["ps.0"] = self.plot_ps_spin
		self.variable["ps.0"] = self.plot_ps

		return frame

#===============================================================================
# Common of Plot Frames
#===============================================================================
class CommonPlotFrame(BasePlotFrame):
	# ----------------------------------------------------------------------
	def __init__(self, master, flair):
		BasePlotFrame.__init__(self, master, flair)
#		splitter = tkExtra.VSplitter(self,0.7)
#		splitter.pack(expand=YES, fill=BOTH)

		# === Top frame with common info ===
		frame = Frame(self) #splitter.topFrame())
		frame.pack(side=TOP, fill=X)

		f = Frame(frame)
		f.pack(side=TOP, fill=X)

#		Label(f, text="Title:").pack(side=LEFT)

		self.titleOptions = BooleanVar()
		b = Checkbutton(f, text="Title: ", #+Unicode.BLACK_DOWN_POINTING_TRIANGLE,
				variable=self.titleOptions,
				selectcolor="LightYellow",
				indicatoron=FALSE,
				command=self.showTitleOptions)
		tkExtra.Balloon.set(b, "Show plot options")
		b.pack(side=LEFT)

		b = Entry(f, background="White")
		tkExtra.Balloon.set(b, "Plot Title")
		b.pack(side=LEFT, expand=YES, fill=X)
		self.widget["title"] = b

		var = IntVar()
		var.set(0)
		Label(f, text="Display:").pack(side=LEFT)
		b = Spinbox(f, text=var, from_=0,
			to=100, background="White", width=3)
		tkExtra.Balloon.set(b, "Gnuplot display window index. " \
			"Change index to display many windows")
		b.pack(side=LEFT)
		self.widget["index"] = b
		self.variable["index"] = var

		# ==== Title Options Frame ====
		self.titleOptionsFrame = LabelFrame(frame, text="Options ", foreground="DarkBlue")

		# --- 1st row ---
		f = Frame(self.titleOptionsFrame)
		f.pack(side=TOP, fill=X)
		self.optionsFrame(f, "title")

		# --- 2nd row ---
		f = Frame(self.titleOptionsFrame)
		f.pack(side=TOP, fill=X)

		# - grid -
		var = IntVar()
		var.set(0)
		b = Checkbutton(f, text="grid",
				foreground="DarkBlue",
				variable=var,
				anchor=W)
		b.pack(side=LEFT, padx=5)
		tkExtra.Balloon.set(b,"Display plot grid")
		self.widget["grid"] = b
		self.variable["grid"] = var

		# - aspect -
		Label(f, text=" aspect:", foreground="DarkBlue").pack(side=LEFT)
		b = tkExtra.FloatEntry(f, background="white", width=7)
		tkExtra.Balloon.set(b, "Plotting aspect ratio")
		b.pack(side=LEFT)
		self.widget["ratio"] = b

		# - auto -
		b = Button(f, text="Auto", foreground="DarkBlue",
				command=self.calcAspectRatio,
				pady=1)
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b,
			"Calculate aspect ratio, using the plot range information. " \
			"Make a plot before clicking")

		# line type
		Label(f, text="  lines:", foreground="DarkBlue").pack(side=LEFT)
		b = tkExtra.Combobox(f, width=6)
		b.fill(Plot.LINE_TYPES)
		b.set(Plot.LINE_TYPES[0])
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b,"Line type style dashed or solid")
		self.widget["linetype"] = b

		# === Top frame with common info ===
		topframe = Frame(self) #splitter.topFrame())
		topframe.pack(fill=X)

		groupframe = LabelFrame(topframe, text="Axes")
		groupframe.pack(side=LEFT, fill=X, expand=YES, padx=3, pady=3)

		row = 0
		col = 0

		self.axesButton = Button(groupframe,
				text=Unicode.BLACK_DOWN_POINTING_TRIANGLE,
				command=self.axesMenu,
				width=3,
				padx=0,
				pady=0)
		self.axesButton.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(self.axesButton, "Select axes to display")
		col += 1
		for label,anchor in [
			      ("Label",	W),
			      ("Log",	W),
			      ("Min",	CENTER),
			      (None,	None),
			      ("Max",	CENTER)]:
			if label:
				Label(groupframe,
					text=label,
					anchor=anchor).grid(row=row, column=col, sticky=NSEW)
			col += 1

		# Axes frame
		self.axesVar = [("z", BooleanVar()),
				("cb",BooleanVar()),
				("x2",BooleanVar()),
				("y2",BooleanVar())]

		# Options
		self.axesOptions = []

		for axis in _AXES:
			row += 1
			# name
			col = 0

			#b = Label(groupframe, text="%s:"%(axis), anchor=E)
			#b.grid(row=row, column=col, sticky=EW)
			#self.widget["%sname"%(axis)] = b

			# Options button
			optvar = BooleanVar()
			b = Checkbutton(groupframe, text="%s:"%(axis), #text=Unicode.BLACK_DOWN_POINTING_TRIANGLE,
					anchor=E,
					variable=optvar,
					selectcolor="LightYellow",
					indicatoron=FALSE,
					width=2,
					command=self.axesShowOptions)
			tkExtra.Balloon.set(b, "Show axis options")
			self.widget["%sname"%(axis)] = b
			#self.widget["%sshow"%(axis)] = b
			b.grid(row=row, column=col, sticky=EW)

			# Label
			col += 1
			b = Entry(groupframe, background="White")
			b.grid(row=row, column=col, sticky=EW)
			tkExtra.Balloon.set(b, "%s axis label"%(axis))
			self.widget["%slabel"%(axis)] = b

			# Log
			col += 1
			var = BooleanVar()
			var.set(False)
			b = Checkbutton(groupframe,
					text="",
					variable=var,
					anchor=W)
			b.grid(row=row, column=col, sticky=EW)
			tkExtra.Balloon.set(b, "Log/Linear")
			self.widget["%slog"%(axis)] = b
			self.variable["%slog"%(axis)] = var

			if axis == "cb":	# Special handling for color band
				self.cblogVar = var
				var.trace("w", self.cpdChange)

			# Limits:  min - max
			col += 1
			if axis == "cb":	# Special handling for color band
				self.cb_min_var = var = StringVar()
				self.cb_min_var.trace("w", self.cpdChange)
			else:
				var = None

			b = tkExtra.FloatEntry(groupframe,
					textvariable=var,
					width=15,
					background="White")
			b.grid(row=row, column=col, sticky=EW)
			tkExtra.Balloon.set(b, "%s axis minimum range"%(axis))
			self.widget["%smin"%(axis)] = b

			col += 1
			b = Label(groupframe, text="-")
			b.grid(row=row, column=col, sticky=EW)
			self.widget["%sdash"%(axis)] = b

			col += 1
			b = tkExtra.FloatEntry(groupframe, width=15, background="White")
			b.grid(row=row, column=col, sticky=EW)
			tkExtra.Balloon.set(b, "%s axis maximum range"%(axis))
			self.widget["%smax"%(axis)] = b

			# Special handling for color band
			row += 1
			frame = Frame(groupframe)
			# grid now now
			self.axesOptions.append((axis,optvar,row,frame))

			if axis == "cb":
				row += 1
				f = Frame(frame)
				f.pack(side=TOP, fill=X)

				# ---
				var = IntVar()
				var.set(1)
				b = Checkbutton(f, text="Round", variable=var)
				b.pack(side=RIGHT)
				tkExtra.Balloon.set(b, "Round limits so colors match always a decade boundary")
				self.widget["cbround"]   = b
				self.variable["cbround"] = var
				var.trace("w", self.cpdChange)

				# ---
				var = IntVar()
				var.set(3)
				b = Spinbox(f, text=var,
						from_=0, to=50,
						background="White",
						width=3)
				tkExtra.Balloon.set(b, "Colors per decade")
				b.pack(side=RIGHT)
				Label(f, text="  CPD:", foreground="DarkBlue").pack(side=RIGHT)
				self.widget["cpd"]   = b
				self.variable["cpd"] = var
				var.trace("w", self.cpdChange)

				# ---
				var = IntVar()
				var.set(30)
				b = Spinbox(f, text=var,
						from_=4, to=256,
						background="White",
						width=4)
				tkExtra.Balloon.set(b, "Total number of colors")
				b.pack(side=RIGHT)
				Label(f, text="  Colors:", foreground="DarkBlue").pack(side=RIGHT)
				self.widget["cbcolors"]   = b
				self.variable["cbcolors"] = var
				var.trace("w", self.cpdChange)

				# ---
				b = tkExtra.Combobox(f, width=10)
				b.fill(Plot.PALETTE)
				b.set(Plot.PALETTE[0])
				b.pack(side=RIGHT)
				Label(f, text="Palette:", foreground="DarkBlue").pack(side=RIGHT)
				self.widget["cbpalette"] = b

				# Replace frame with the axes options
				frame = Frame(frame)
				frame.pack(side=BOTTOM, fill=X)

			# option field
			self.optionsFrame(frame, axis, 0, "label")
			self.optionsFrame(frame, axis+"tics", 1, "tics")

		# -----
		groupframe.grid_columnconfigure(1, weight=1)
		self.axesShow()

		# =====
#		self.plot_keys = IntVar()
#		self.plot_keys.set(1)
#		b = Checkbutton(frame, text="legend", variable=self.plot_keys,
#				anchor=W)
#		b.pack(fill=X)
#		tkExtra.Balloon.set(b, "Display label keys")
#
#		# =====
#		Label(frame, text="Width:").grid(row=row, column=col, sticky=E)
#		col += 1
#		self.plot_sizex = tkExtra.FloatEntry(frame, background="white", width=4)
#		tkExtra.Balloon.set(self.plot_sizex, "Plotting width")
#		self.plot_sizex.grid(row=row, column=col, sticky=EW)
#
#		Label(frame, text="Height:").grid(row=row, column=col, sticky=E)
#		self.plot_sizey = tkExtra.FloatEntry(frame, background="white", width=4)
#		tkExtra.Balloon.set(self.plot_sizey, "Plotting height")
#		self.plot_sizey.grid(row=row, column=col, sticky=EW)
#
		# === Common Main Frame ===
		self.common = Frame(self) #splitter.topFrame())
		self.common.pack(fill=BOTH)

		# === Commands frame with common info ===
		frame = Frame(self)
		frame.pack(side=BOTTOM,fill=BOTH, expand=YES, pady=3)

		b = Text(frame, #splitter.bottomFrame(),
				background="white", height=5, undo=True)
		b.pack(side=LEFT,fill=BOTH, expand=YES)
		self.widget["commands"] = b

		sb = Scrollbar(frame, #splitter.bottomFrame(),
				takefocus=False, command=b.yview)
		sb.pack(side=RIGHT, fill=Y)
		b.config(yscrollcommand=sb.set)

		self.bindToggle()

	# ----------------------------------------------------------------------
	def optionsFrame(self, frame, prefix, row=0, label=""):
		col = 0
		if label:
			Label(frame, text=label, fg="DarkBlue").grid(row=row, column=col, sticky=E)

		# font
		col += 1
		Label(frame, text="font:", fg="DarkBlue").grid(row=row, column=col, sticky=E)
		col += 1
		b = tkExtra.Combobox(frame, width=12)
		b.grid(row=row, column=col, sticky=EW)
		b.fill(Plot.FONTS)
		b.set(Plot.FONTS[0])
		tkExtra.Balloon.set(b,"Font family")
		self.widget["%sfont"%(prefix)] = b

		# size
		col += 1
		b = tkExtra.Combobox(frame, width=4)
		b.grid(row=row, column=col, sticky=EW)
		b.fill(Plot.SIZE)
		b.set(Plot.SIZE[0])
		tkExtra.Balloon.set(b,"Font size")
		self.widget["%ssize"%(prefix)] = b

		# color
		col += 1
		Label(frame, text="color:", fg="DarkBlue").grid(row=row, column=col, sticky=EW)
		col += 1
		b = tkExtra.Combobox(frame, width=12)
		b.grid(row=row, column=col, sticky=EW)
		b.fill(PlotEngine._COLOR_NAMES)
		b.set(PlotEngine._COLOR_NAMES[0])
		tkExtra.Balloon.set(b,"Font color")
		self.widget["%scolor"%(prefix)] = b

		# --- 3rd row ---
		# options
		col += 1
		Label(frame, text="options:", fg="DarkBlue").grid(row=row, column=col, sticky=EW)
		col += 1
		b = Entry(frame, background="White")
		b.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(b, "Additional options")
		self.widget["%soptions"%(prefix)] = b

		frame.grid_columnconfigure(col, weight=1)

	# ----------------------------------------------------------------------
	def showTitleOptions(self):
		if self.titleOptions.get():
			self.titleOptionsFrame.pack(side=BOTTOM, fill=X)
		else:
			self.titleOptionsFrame.pack_forget()

	# ----------------------------------------------------------------------
	# Display axes menu to select from
	# ----------------------------------------------------------------------
	def axesMenu(self):
		menu = Menu(self, tearoff=0)
		for axis,var in self.axesVar:
			menu.add_checkbutton(label=axis,  variable=var,  command=self.axesShow)
		menu.tk_popup(	self.axesButton.winfo_rootx(),
				self.axesButton.winfo_rooty() + self.axesButton.winfo_height())

	# ----------------------------------------------------------------------
	# Hide/Show axes
	# ----------------------------------------------------------------------
	def axesShow(self):
		for j,(axis,var) in enumerate(self.axesVar):
			for i,label in enumerate(["name", "label", "log", "min", "dash", "max"]): #, "show"]):
				name = axis+label
				if var.get():
					# Show
					self.widget[name].grid(row=j*2+5, column=i, sticky=EW)
				else:
					# Hide
					self.widget[name].grid_forget()
			if not var.get():
				for a,v,r,f in self.axesOptions:
					if axis==a: v.set(False)
		self.axesShowOptions()

	# ----------------------------------------------------------------------
	def axesShowOptions(self):
		for axis,var,row,frame in self.axesOptions:
			if var.get():
				# Show
				frame.grid(row=row, column=1, columnspan=6,sticky=NSEW)
			else:
				# Hide
				frame.grid_forget()

	# ----------------------------------------------------------------------
	# set fields
	# ----------------------------------------------------------------------
	def set(self, plot):
		BasePlotFrame.set(self, plot)

		# Special values
		for i,(axis, var) in enumerate(self.axesVar):
			axis = plot.get("%stics"%(axis),0)
			try:
				val = int(axis)
			except:
				val = int(axis=="True")
			if val == 0:
				# find option frame to remove
				a,v,r,f = self.axesOptions[i+2]
				v.set(0)
				f.grid_forget()
			var.set(val)
		self.axesShow()

	# ----------------------------------------------------------------------
	# get fields
	# ----------------------------------------------------------------------
	def get(self, plot):
		BasePlotFrame.get(self, plot)

		# Special values
		for axis, var in self.axesVar:
			plot["%stics"%(axis)] = int(var.get())

	# ----------------------------------------------------------------------
	# Override method to calculate aspect ratio from plot extends
	# ----------------------------------------------------------------------
	def calcAspectRatio(self):
		if not self.flair.plotEngine: return
		xlow, xhigh = self.flair.plotEngine.getRange("x")
		if xhigh is None: return
		ylow, yhigh = self.flair.plotEngine.getRange("y")
		if yhigh is None: return
		try: ratio = (yhigh-ylow) / (xhigh-xlow)
		except: return
		self.widget["ratio"].set(ratio)

	# ----------------------------------------------------------------------
	# CPD Change
	# ----------------------------------------------------------------------
	def cpdChange(self, a, b, c):
		if int(self.cblogVar.get()):
#			self.widget["cbmax"].config(state=DISABLED)
			state = NORMAL
		else:
#			self.widget["cbmax"].config(state=NORMAL)
			state = DISABLED
		self.widget["cpd"].config(state=state)
		self.widget["cbround"].config(state=state)
		if state == DISABLED: return

		try:
			cpd = float(self.variable["cpd"].get())
		except ValueError:
			cpd = 0.0

		if cpd < 1.0: return

		# Calculate max
		try:
			cbmin  = float(self.widget["cbmin"].get())
		except ValueError:
			return
		try:
			colors = float(self.widget["cbcolors"].get())
		except ValueError:
			return

		if self.variable["cbround"].get():
			# round minimum to lowest round color value
			try:
				cbmin = math.pow(10.0, \
					math.floor(math.log10(cbmin)*cpd)/cpd)
			except (ValueError, OverflowError):
				pass

		if cbmin>0.0:
			try:
				self.widget["cbmax"].set(str(cbmin*math.pow(10.0, colors/cpd)))
			except (ValueError, OverflowError):
				pass
