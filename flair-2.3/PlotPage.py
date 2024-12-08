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

import os
import re
import math
import string
import io
from log import say
try:
	from io import StringIO
except ImportError:
	from io import StringIO
try:
	import pickle as pickle
except ImportError:
	import pickle

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkinter as tk
import undo
import bmath
import Input
import Ribbon
import Manual
import Project
import FlairRibbon

import Unicode
import tkFlair
import tkExtra
import tkDialogs

import Plot
import GPPlot
import MPPlot
import UsrPlot
import BasePlot
import Usr2Plot
import UsrbinPlot
import GeometryPlot
import UserdumpPlot
import ResnucleiPlot

#from flair import PLOT_ENGINES
PLOT_ENGINES = ["Gnuplot", "Matplotlib"]

_PLOT_CLIP = tkFlair._FLAIRCLIP + "<plot>"

#===============================================================================
# Action Group
#===============================================================================
class _ActionGroup(FlairRibbon.FlairMenuGroup):
	#----------------------------------------------------------------------
	def createMenu(self):
		menu = Menu(self, tearoff=0)

		for i,n in enumerate(PLOT_ENGINES):
			menu.add_radiobutton(
					label=n,
					image=tkFlair.icons["plot"], compound=LEFT,
					variable=self.page._plotEngine,
					value=n)
		return menu

#===============================================================================
# Plot List Frame
#===============================================================================
class PlotListPage(FlairRibbon.FlairPage):
	"""Create plots of geometry or from the detectors defined in the input."""

	_name_ = "Plot"
	_icon_ = "plot"

	#----------------------------------------------------------------------
	def init(self):
		self._plotEngine = StringVar()
		self._plotEngine.set(PLOT_ENGINES[0])
		self.undoredo = undo.UndoRedo()
		self.selected = None

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Plot List ===========
		group = Ribbon.LabelGroup(self.ribbon, "Plot List")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row = 0,0
		menulist = [	("Geometry", "GEOMETRY",  self.addGeometry),
				("USR-1D",   "USR-1D",    self.addUsr1D),
				("USR-2D",   "USR-2D",    self.addUsr2D),
				("USRBIN",   "USRBIN",    self.addUsrbin),
				("USERDUMP", "USERDUMP",  self.addUserdump),
				("RESNUCLE", "RESNUCLE",  self.addResnuclei)]

		b = Ribbon.MenuButton(group.frame, menulist,
				image=tkFlair.icons["add32"],
				text="Add\n"+Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				compound=TOP,
				command=self.add,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, rowspan=3, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add new plot [Ins]")

		# ---
		col,row = 1,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["wizard32"],
				text="Oz",
				compound=TOP,
				anchor=W,
				command=self.wizard,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Automatically create plots")

		# ---
		col,row = 2,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["up"],
				text="Move Up",
				compound=LEFT,
				command=self.moveUp,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Move selected plots up [Ctrl-Up]")

		# ---
		col,row = 2,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["down"],
				text="Move Down",
				compound=LEFT,
				command=self.moveDown,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Move selected plots down [Ctrl-Down]")

		# ---
		col,row = 2,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["rename"],
				text="Rename",
				compound=LEFT,
				command=self.rename,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Rename plot [F2]")

		# ---
		col,row = 3,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Delete",
				compound=LEFT,
				anchor=W,
				command=self.delete,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Delete selected plots [Del]")

		# ---
		col,row = 3,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clone"],
				text="Clone",
				compound=LEFT,
				anchor=W,
				command=self.clone,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Clone selected plots [Ctrl-D]")

		# ========== Plot ===========
		group = _ActionGroup(self.ribbon, self, "Action")
		#group = Ribbon.LabelGroup(self.ribbon, "Action")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()
		group.frame.grid_columnconfigure(2, weight=1)

		# ---
		col,row = 0,0
		menulist = []
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["save"],
				command=self.save,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Save plot")

		col,row = 1,0
		menulist = []
		for fmt in sorted(tkFlair._PLOT_FORMAT.keys()):
			cmd = lambda s=self,f=fmt: s.save(f)
			menulist.append((fmt, "save", cmd))
		b = Ribbon.MenuButton(group.frame, menulist,
				image=tkFlair.icons["triangle_down"],
				text="Save",
				compound=RIGHT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Save plot as")

		# ---
		col,row = 0,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["print"],
				text="Print",
				compound=LEFT,
				command=self.hardcopy,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Print selected plots")

		# ---
		col,row = 0,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["document"],
				text="Notes",
				compound=LEFT,
				command=self.toNotes,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Insert to notes")

		# ---
		col,row = 2,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clean32"],
				text="Clean",
				compound=TOP,
				command=self.clean,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Remove all generated files")

		# ---
		col,row = 3,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["plot32"],
				text="Plot",
				compound=TOP,
				command=self.execute,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, rowspan=3, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Do the plot. Generate files when needed")

	#----------------------------------------------------------------------
	# Create Project page
	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self)

		self.splitter = tkExtra.HSplitter(self.frame,
					tkFlair.getInt(tkFlair._PAGE_SECTION,
						self.name+".split", 240),
						True)
		self.splitter.pack(expand=YES, fill=BOTH)

		# -----------
		self.listbox = tkExtra.ImageListbox(self.splitter.leftFrame(),
					width=16)
		self.listbox.pack(fill=BOTH, expand=YES, side=LEFT)
		sb = Scrollbar(self.splitter.leftFrame(),
				orient=VERTICAL,
				command=self.listbox.yview)
		sb.pack(fill=Y, side=RIGHT)
		self.listbox.config(yscrollcommand=sb.set)

		self.listbox.bind('<ButtonRelease-1>', self.release)
		self.listbox.bind('<Control-Key-e>',   self.rename)
		self.listbox.bind('<F2>',              self.rename)
		self.listbox.bind('<Return>',          self.rename)
		self.listbox.bind('<Control-Return>',  self.execute)
		self.listbox.bind('<KP_Enter>',        self.rename)
		self.listbox.bind('<Delete>',          self.delete)
		self.listbox.bind('<Control-Key-Up>',  self.moveUp)
		self.listbox.bind('<Control-Key-Down>',self.moveDown)
		self.listbox.bind('<<ListboxSelect>>', self.select)
		self.listbox.bind('<Button-3>',        self.popup)
		self.listbox.bind('<<Copy>>',          self.copy)
		self.listbox.bind('<<Cut>>',           self.cut)
		self.listbox.bind('<<Paste>>',         self.paste)

		self.listbox.focus_set()
		self.lastrelease = (-1, 0.0)

		# -----------
		self.base = BasePlot.CommonPlotFrame(self.splitter.rightFrame(), self.flair)
		self.base.pack(fill=BOTH, expand=YES)

		# ----- Create Plot Frames ----
		self.frames = {}
		for C in [GeometryPlot.GeometryPlotFrame,
			  UsrPlot.Usr1DPlotFrame,
			  Usr2Plot.Usr2DPlotFrame,
			  UsrbinPlot.UsrbinPlotFrame,
			  UserdumpPlot.UserdumpPlotFrame,
			  ResnucleiPlot.ResnucleiPlotFrame]:
			frame = C(self.base.common, self.flair)
			self.frames[frame.type()] = frame
		self.displayed = None

	# ----------------------------------------------------------------------
	def help(self):
		Manual.show("F5")

	# ----------------------------------------------------------------------
	def canUndo(self):	return self.undoredo.canUndo()
	def canRedo(self):	return self.undoredo.canRedo()
	def resetUndo(self):	self.undoredo.reset()

	# ----------------------------------------------------------------------
	def undo(self, event=None):
		if self.canUndo():
			self.undoredo.undo()
			self.selected = None
			self.refresh()

	# ----------------------------------------------------------------------
	def redo(self, event=None):
		if self.canRedo():
			self.undoredo.redo()
			self.selected = None
			self.refresh()

	# ----------------------------------------------------------------------
	def copy(self, event=None):
		# Operate only on the listbox
		if self.page.focus_get() is not self.listbox:
			FlairRibbon.FlairPage.copy(self,event)
		else:
			self.get()

			plots = []
			for i in self.selected:
				plot = self.project.plots[i]
				plots.append((plot.name, plot.type, plot.var))

			# Write to clipboard
#			sio = StringIO()
			sio = io.BytesIO()
			sio.write(_PLOT_CLIP)
			pickler = pickle.Pickler(sio)
			pickler.dump(plots)
			self.page.clipboard_clear()
			self.page.clipboard_append(sio.getvalue())
			return "break"

	# ----------------------------------------------------------------------
	def cut(self, event=None):
		# Operate only on the listbox
		if self.page.focus_get() is not self.listbox: return
		self.copy()
		self.delete()
		return "break"

	# ----------------------------------------------------------------------
	def paste(self, event=None):
		# Operate only on the listbox
		if self.page.focus_get() is not self.listbox:
			FlairRibbon.FlairPage.paste(self,event)
		else:
			self.get()
			self.current = None

			try: clipboard = self.page.selection_get(selection='CLIPBOARD')
			except: return
			if clipboard.startswith(_PLOT_CLIP):
				# Pickler format
#				unpickler = pickle.Unpickler(StringIO(clipboard[len(_PLOT_CLIP):]))
				unpickler = pickle.Unpickler(io.BytesIO(clipboard[len(_PLOT_CLIP):]))
				try:
					plots = unpickler.load()
				except:
					return "break"
			else:
				self.flair.notify("Cannot paste","Invalid content in clipboard",tkFlair.NOTIFY_WARNING)
				return "break"

			active = self.listbox.index(ACTIVE)
			if active < 0 or active >= self.listbox.size()-1:
				active = -1
			self.listbox.selection_clear(0,END)

			new = []
			for n,t,v in plots:
				plot = Project.PlotInfo(n)
				plot.type = t
				plot.var = v
				if active < 0:
					new.append(len(self.project.plots))
					self.project.plots.append(plot)
				else:
					active += 1
					new.append(active)
					self.project.plots.insert(active,plot)

			self.refresh()
			self.listbox.selection_clear(0,END)
			for s in new: self.listbox.selection_set(s)
			self.listbox.activate(new[0])
			self.listbox.see(new[0])
			self.select()
			return "break"

	# ----------------------------------------------------------------------
	# Delete plot
	# ----------------------------------------------------------------------
	def delete(self, event=None):
		if self.page.focus_get() is not self.listbox: return

		self.listbox.focus_set()
		lst = list(map(int, self.listbox.curselection()))
		if not lst: return
		lst.reverse()
		act = lst[0]
		for i in lst:
			del self.project.plots[i]
			self.listbox.delete(i)
		if act >= self.listbox.size(): act = self.listbox.size()-1
		if act >= 0:
			self.listbox.activate(act)
			self.listbox.selection_set(act)
		del self.selected[:]
		self.select()
		self.setModified()
		self.flair.refreshGeoedit()
		return "break"

	# ----------------------------------------------------------------------
	# Clone plots
	# ----------------------------------------------------------------------
	def clone(self, event=None):
		self.listbox.focus_set()
		sel = list(map(int,self.listbox.curselection()))
		if len(sel) == 0: return
		oldSize = self.listbox.size()
		self.listbox.selection_clear(0,END)
		insertPoint = sel[-1]+1
		for i in reversed(sel):
			plot = self.project.plots[i].clone()
			self.project.findUniquePlotName(plot)
			if insertPoint >= self.listbox.size():
				self.project.plots.append(plot)
				self.listbox.insert(END, tkFlair.icons[plot.type.upper()], plot.name)
			else:
				self.project.plots.insert(insertPoint,plot)
				self.listbox.insert(insertPoint, tkFlair.icons[plot.type.upper()], plot.name)
			self.listbox.selection_set(insertPoint)
		self.setModified()
		self.listbox.activate(i+1)
		self.flair.refreshGeoedit()
		self.select()

	# ----------------------------------------------------------------------
	# Activate and set focus
	# ----------------------------------------------------------------------
	def activate(self):
		#self.refresh()
		FlairRibbon.FlairPage.activate(self)

	# ----------------------------------------------------------------------
	def refresh(self):
		FlairRibbon.FlairPage.refresh(self)
		if self.flair.geoedit:
			self.flair.geoedit.saveState()
		saveSelected = list(map(int, self.listbox.curselection()))
		self.listbox.delete(0,END)
		for plot in self.project.plots:
			self.listbox.insert(END,tkFlair.icons[plot.type.upper()], plot.name)
		self.listbox.selection_clear(0,END)
		for s in saveSelected:
			self.listbox.selection_set(s)
		if len(saveSelected)>0:
			self.listbox.activate(saveSelected[0])
		else:
			self.listbox.selection_set(0)
			self.listbox.activate(0)
		self.select()

	# ----------------------------------------------------------------------
	# Reset before loading new projects
	# ----------------------------------------------------------------------
	def reset(self):
		if self.page is None: return
		self.selected = None
		self.listbox.delete(0,END)

	# ----------------------------------------------------------------------
	# Get information from the plots
	# ----------------------------------------------------------------------
	def get(self):
		if not self.selected: return
		if self.displayed is None: return
		for i in self.selected:
			plot = self.project.plots[i]

			# Remember for undo
			oldplot = plot.clone()

			# Get information
			try:
				self.base.get(plot)
				if plot.type == self.displayed.type():
					self.displayed.get(plot)
			except TclError:
				continue

			#if Input._developer:
			#	say("----old---")
			#	oldplot.write(sys.stdout)
			#	say("\n---new---")
			#	plot.write(sys.stdout)

			if plot.hash() != oldplot.hash():
				undoinfo = (plot.setUndo, oldplot)
				self.undoredo.addUndo(undoinfo)
				self.flair.setModified()

	# ----------------------------------------------------------------------
	# Selection changed in the listbox
	# ----------------------------------------------------------------------
	def select(self, event=None):
		self.get()

		# As a precaution remove everything
		for child in self.base.common.winfo_children():
			child.pack_forget()

		self.selected = list(map(int, self.listbox.curselection()))
		if len(self.selected) == 0: return

		# Find common type
		plot = self.project.plots[self.selected[0]]
		t = plot.type
		for i in self.selected:
			if t != self.project.plots[i].type:
				t = None
				break

		# Display the new frame
		if t:
			self.displayed = self.frames[plot.type]
			self.displayed.pack(expand=YES, fill=BOTH)

		# Fill values with plot type
		self.displayed.setState(NORMAL)
		self.base.setState(NORMAL)

		# Fill in reverse order so maybe the set to modify the plot like in usrbin
		self.displayed.set(plot)
		self.base.set(plot)

		if len(self.selected)>1:
			self.base.setState(DISABLED)
			self.displayed.setState(DISABLED)
		else:
			self.base.setState(NORMAL)
			self.displayed.setState(NORMAL)

	# ----------------------------------------------------------------------
	# Add a new plot
	# ----------------------------------------------------------------------
	def add(self, typ=None):
		plot = Project.PlotInfo("")
		self.project.plots.append(plot)
		plot["title"] = "Plot #%d" % (len(self.project.plots))
		if self.project.name == "":
			plot.name   = "plot%02d" % (len(self.project.plots))
		else:
			(fn, ext) = os.path.splitext(self.project.name)
			plot.name   = "%s_plot%02d" % (fn, len(self.project.plots))
		plot.format = Plot.DEFAULT_FORMAT
		if typ:
			plot.type = typ
		else:
			plot.type = Plot.PLOT_TYPES[0]

		# Detector specific
		if plot.type == "USRBIN":
			plot["cbtics"]   = 1
			plot["cblog"]    = 1
			plot["cbcolors"] = 30

		elif plot.type == "USR-1D":
			plot["ndetectors"] = 1
			plot["name.0"] = "#Detector 1"
			plot["det.0"]  = 0

		elif plot.type == "USR-2D":
			plot["ztics"]  = 1
			plot["cbtics"] = 1

		elif plot.type == "RESNUCLE":
			plot["cbtics"] = 1

		elif plot.type == "USERDUMP":
			plot["type"] = "Histogram"

		self.listbox.insert(END,tkFlair.icons[plot.type.upper()], plot.name)
		self.setModified()
		self.listbox.selection_clear(0,END)
		self.listbox.selection_set(END)
		self.listbox.activate(END)
		self.listbox.see(END)
		self.listbox.focus_set()
		self.select()
		self.flair.refreshGeoedit()
		self.lastrelease = (self.listbox.index(END), 0)

	# ----------------------------------------------------------------------
	def addGeometry(self):	self.add("Geometry")
	def addUsr1D(self):	self.add("USR-1D")
	def addUsr2D(self):	self.add("USR-2D")
	def addUsrbin(self):	self.add("USRBIN")
	def addUserdump(self):	self.add("USERDUMP")
	def addResnuclei(self):	self.add("RESNUCLE")

	# ----------------------------------------------------------------------
	# Swap two items in the list
	# ----------------------------------------------------------------------
	def swap(self, a, b):
		lb = self.listbox
		if a>b:
			a, b = b, a

		ima, at = lb.elicit(a);	al = self.project.plots[a]
		imb, bt = lb.elicit(b);	bl = self.project.plots[b]

		lb.delete(b)
		lb.delete(a)

		lb.insert(a, imb, bt);	self.project.plots[a] = bl
		lb.insert(b, ima, at);	self.project.plots[b] = al

	# ----------------------------------------------------------------------
	# Move up select items by one
	# ----------------------------------------------------------------------
	def moveUp(self, event=None):
		self.listbox.focus_set()
		self.get()
		self.selected = None
		lb  = self.listbox
		for i in map(int,lb.curselection()):
			if i==0: continue
			prev = i-1
			if not lb.selection_includes(prev):
				act = lb.index(ACTIVE)
				self.swap(prev,i)
				lb.selection_set(prev)
				if act == i: lb.activate(prev)
		self.selected = list(map(int, self.listbox.curselection()))
		return "break"

	# ----------------------------------------------------------------------
	# Move down select items by one
	# ----------------------------------------------------------------------
	def moveDown(self, event=None):
		self.listbox.focus_set()
		self.get()
		self.selected = None
		lb  = self.listbox
		sz  = lb.size()-1
		lst = list(map(int,lb.curselection()))
		lst.reverse()
		for i in lst:
			if i >= sz: continue
			next = i+1
			if not lb.selection_includes(next):
				act = lb.index(ACTIVE)
				self.swap(i,next)
				lb.selection_set(next)
				if act == i: lb.activate(next)
		self.selected = list(map(int, self.listbox.curselection()))
		return "break"

	# ----------------------------------------------------------------------
	def wizard(self):
		"""scan input file for possible plots"""
		if self.project.name == "":
			pattern = "%s"
		else:
			(fn, ext) = os.path.splitext(self.project.name)
			pattern = "%s_%%s" % (fn)

		# Add one geometry plot if doesn't exist
		for p in self.project.plots:
			if p.type == Plot.PLOT_TYPES[0]: break
		else:
			plot = Project.PlotInfo(pattern % ("geometry"))
			plot["title"] = "Geometry %s"%(self.project.title)
			plot.type     = Plot.PLOT_TYPES[0]
			plot.format   = ".eps"
			self.project.plots.append(plot)

		# For default Run find all usrinfos
		run = self.project.runs[0]
		for usr in run.usrinfo:
			name = os.path.splitext(usr.name())[0]

			plot = Project.PlotInfo("%s_plot"%(name))
			self.project.plots.append(plot)
			if usr.type=="b":
				plot.type = "USRBIN"
				plot["datafile"] = usr.name()
				plot["geo"]      = BasePlot._GEO_AUTO
				plot["cbtics"]   = 1
				plot["cblog"]    = 1
				plot["cbcolors"] = 30

			elif usr.type=="r":
				plot.type = "RESNUCLE"
				plot["datafile"] = usr.name()

			else:
				plot.type = "USR-1D"
				plot["xlog"] = 1
				plot["ylog"] = 1
				plot["ndetectors"] = 1

				tablis = "%s_tab.lis"%(name)

				plot["name.0"] = "Detector 1"
				plot["file.0"] = tablis
				plot["det.0"]  = "0"
				plot["x.0"]    = "2"
				plot["y.0"]    = "0"
				plot["lw.0"]   = "2"
				plot["lc.0"]   = "red"
				plot["with.0"] = "histerror"
			plot["title"] = "%s %s"%(plot.type, name.replace("_"," "))
			plot.format   = ".eps"

		self.refresh()
		self.setModified()

	# ----------------------------------------------------------------------
	def release(self, event=None):
		listbox = event.widget
		act = listbox.nearest(event.y)
		if act: self.listbox.activate(act)
		sel = self.listbox.curselection()
		if len(sel) != 1: return
		sel = int(sel[0])

		try: active = listbox.index(ACTIVE)
		except: return
		if active != sel: return

		# Character position under cursor
		if self.lastrelease[0] == active:
			if event.time <= 500:
				self.lastrelease = (active, event.time)
				return
		else:
			self.lastrelease = (active, event.time)
			return


		index = self.listbox.index("@%d,%d"%(event.x,event.y))
		if index is None:
			self.lastrelease = (-1, event.time)
			return

		self.lastrelease = (active, event.time)

		if int(index.split(".")[1]) > 1:
			self.rename(listbox)
		else:
			self.changeTypePopup()

	# ----------------------------------------------------------------------
	def rename(self, event=None):
		active = self.listbox.index(ACTIVE)
		if active is None: return
		plot = self.project.plots[active]

		# Edit the filename
		edit = tkExtra.InPlaceFile(self.listbox, value=plot.name+plot.format,
			title="Save plot file as",
			filetypes=[
				("Encapsulated Postscript",   "*.eps"),
				("Portable Network Graphics", "*.png"),
				("Scalable Vector Graphics",  "*.svg"),
				("All","*")])
		if edit.value is None:
			# Break not to insert the <Return> in the Text and modify the text
			return "break"

		fn = self.project.relativePath(edit.value)
		name, ext = os.path.splitext(fn)
		if name != plot.name:
			plot.name = name
		self.setModified()
		if len(ext) > 0: plot.format = ext
		self.listbox.delete(active)
		self.listbox.insert(active,tkFlair.icons[plot.type.upper()], plot.name)
		self.listbox.selection_set(active)
		self.listbox.activate(active)
		self.listbox.see(active)
		# Break not to insert the <Return> in the Text and modify the text
		return "break"

	# ----------------------------------------------------------------------
	def popup(self, event):
		menu=Menu(self.listbox, tearoff=0)

		submenu = Menu(menu)
		menu.add_command(label="Rename", image=tkFlair.icons["rename"], compound=LEFT,
				command=self.rename)

		menu.add_cascade(label="Add", image=tkFlair.icons["add"], compound=LEFT,
				underline=0, menu=submenu)
		for t in Plot.PLOT_TYPES:
			icon = tkFlair.icons[t.upper()]
			submenu.add_command(label=t,
				compound=LEFT, image=icon,
				command=lambda s=self,t=t:s.add(t))

		submenu = Menu(menu)
		menu.add_cascade(label="Change to", image=tkFlair.icons["change"], compound=LEFT,
				underline=7, menu=submenu)

		self._changeTypeMenu(submenu)

		menu.add_command(label="Delete", image=tkFlair.icons["x"], compound=LEFT,
				underline=0, command=self.delete)
		menu.add_command(label="Clone", image=tkFlair.icons["clone"], compound=LEFT,
				underline=0, command=self.clone)

		menu.tk_popup(event.x_root, event.y_root)

	# ----------------------------------------------------------------------
	def changeTypePopup(self):
		menu=Menu(self.listbox, tearoff=0)
		self._changeTypeMenu(menu)
		bbox = self.listbox.dlineinfo(ACTIVE)
		menu.tk_popup(self.listbox.winfo_rootx(),
			      self.listbox.winfo_rooty() + bbox[1] + bbox[3])

	# ----------------------------------------------------------------------
	def _changeTypeMenu(self, menu):
		for t in Plot.PLOT_TYPES:
			icon = tkFlair.icons[t.upper()]
			menu.add_command(label=t,
				compound=LEFT, image=icon,
				command=lambda s=self,t=t:s.changePlotType(t))

	# ----------------------------------------------------------------------
	# Change the Plot type
	# ----------------------------------------------------------------------
	def changePlotType(self, newtype):
		active = self.listbox.index(ACTIVE)
		if active is None: return
		plot   = self.project.plots[active]
		if plot.type != newtype:
			self.setModified()
			plot.type = newtype

			self.listbox.delete(active)
			self.listbox.insert(active,tkFlair.icons[plot.type.upper()], plot.name)
			self.listbox.selection_set(active)
			self.listbox.activate(active)
			self.listbox.see(active)

			self.select()
			self.flair.refreshGeoedit()

	# ----------------------------------------------------------------------
	# generate the plot
	# ----------------------------------------------------------------------
	def execute(self, show=True, save=False, export=False, clean=False,
			printer=None, fmt=None):
		if not self.selected: return

		self.get()
		if self.flair.startPlotEngine(self._plotEngine.get()) is None: return

		doall    = False
		saved    = []
		printed  = []
		exported = []
		indexes  = {}	# display indexes updated
		for i, p in enumerate(self.selected):
			plot = self.project.plots[p]

			# Ask to wait before next plot
			if show and i>0 and not doall and int(plot("index",0)) in indexes:
				ans = tkDialogs.Dialog(self.page.winfo_toplevel(),
					{"title":"Multiple Plot",
					"text":"Continue with the next\nOr do all plots",
					"bitmap": "questhead",
					"default": 0,
					"strings": ("Next","All","Cancel")})

				indexes.clear()

				if ans.num==1:
					doall = True
				elif ans.num in (2,-1):
					break

			self.base.get(plot)
			if plot.type == self.displayed.type():
				self.displayed.get(plot)

			self.busy()
			#try
			log = self.flair.newLog(plot.type, plot.name)

			# FIXME to be changed with one class!!!
			if self._plotEngine.get() == PLOT_ENGINES[0]:
				plotter = GPPlot.GPPlot(
						self.flair.plotEngine,
						self.project,
						plot,
						log)
			elif self._plotEngine.get() == PLOT_ENGINES[1]:
				plotter = MPPlot.MPPlot(
						self.flair.plotEngine,
						self.project,
						plot,
						log)
			if clean: plotter.cleanup()
			if show:
				self.flair.setStatus("Processing plot...", color="DarkGreen")
				self.flair.update_idletasks()
				try:
					err = plotter.show()
				except (OSError, IOError):
					err = sys.exc_info()[1]
					self.flair.notify("OS Error",
						err,
						tkFlair.NOTIFY_ERROR)
				self.displayed.afterPlot()
				self.flair.setStatus("Plot completed")
			else:
				err = None
			self.notBusy()
			#except
			#...
			if err:
				if isinstance(err,tuple):
					err,inp,pid = err
				else:
					inp = pid = None
				self.flair.notify("Plot error",
					"%s\nCheck Output window for additional messages."%(err),
					tkFlair.NOTIFY_ERROR,
					"Output", self.flair.showOutput)
				if inp is not None:
					tkFlair.OutputFilesDialog(self.flair,
							inp, pid, 1, self.project)
					if Project.cleanup:
						self.project.deleteFiles(inp, pid)
						if isinstance(inp,str):
							os.remove(inp+".inp")
						else:
							os.remove(inp.decode()+".inp")
				break

			# Save file if needed
			if save:
				if fmt is None:
					saved.append(plotter.save(plot.format))
				else:
					plot.format = fmt
					saved.append(plotter.save(fmt))

			# Print if needed
			if printer:
				plotter.hardcopy(printer)
				printed.append(plot.name)

			# Update which display indexes have been used up to now
			indexes[int(plot("index",0))] = True

		if clean:
			self.setStatus("Plot files cleaned", "DarkGreen")

		if save:
			self.setStatus("Saved: %s"%(" ".join(saved)), "DarkBlue")

		if printer:
			self.setStatus("Printed: %s"%(" ".join(printed)), "DarkBlue")

	# ----------------------------------------------------------------------
	def replot(self):		self.execute(replot=True)
	def export(self):		self.execute(export=True)
	def save(self,fmt=None):	self.execute(save=True, fmt=fmt)
	def clean(self):		self.execute(show=False, clean=True)

	# ----------------------------------------------------------------------
	# Select plot by name
	# ----------------------------------------------------------------------
	def showPlot(self, name):
		for i,plot in enumerate(self.project.plots):
			if plot.name == name:
				self.listbox.selection_clear(0,END)
				self.listbox.selection_set(i)
				self.listbox.activate(i)
				self.listbox.see(i)
				self.listbox.focus_set()
				self.select()
				break

	# ----------------------------------------------------------------------
	def hardcopy(self):
		printer = tkFlair.printer(self.page)
		if printer is None: return
		self.execute(printer=printer)
		printer.close()

	# ----------------------------------------------------------------------
	def showNotesPage(self):
		self.flair.notifyHide()
		self.flair.tabs.changePage("Flair")

	# ----------------------------------------------------------------------
	# Insert to project notes a link to the image
	# ----------------------------------------------------------------------
	def toNotes(self):
		added = 0

		self.save(".png")

		for p in self.selected:
			plot = self.project.plots[p]
			link = "{img:%s.png}"%(plot.name)
			if link not in self.project.notes:
				self.project.notes += "\n\n%s"%(link)
				added += 1

		if added:
			self.flair.setModified()
			self.flair.refresh("image")
			self.setStatus("%d plots added to notes"%(added))
			self.flair.notify("Plots to Notes",
				"%d plots added to notes"%(added))
		else:
			self.setStatus("NO plots added to notes")
			self.flair.notify("Plots to Notes",
				"No plots added to notes",
				tkFlair.NOTIFY_NORMAL,
				"Notes",
				self.showNotesPage)
