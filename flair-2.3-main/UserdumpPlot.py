# 
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
# Date:	28-Feb-2007

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import math
import time
import string

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import Data

import tkFlair
import tkExtra
import bFileDialog

import Input

import Plot
import BasePlot

#===============================================================================
# USRBIN Data Plot Frame
#===============================================================================
class UserdumpPlotFrame(BasePlot.BasePlotFrame):
	def __init__(self, master, flair):
		BasePlot.BasePlotFrame.__init__(self, master, flair)

		# --- Variables ---
		self.case    = -1
		self.mgfile  = ""
		self.limfile = ""
		self.dmin    = None
		self.dmax    = None

		# === Fill main frame ===
		# === Axes ===
#		frame = self.axesRangeFrame(self, True)
#		frame.pack(side=TOP, expand=YES, fill=X, padx=3)

		# === File & Detector ===
		frame = LabelFrame(self, text="Collision File ")
		frame.pack(side=TOP, expand=YES, fill=BOTH, padx=3)

		row, col = 0, 0
		Label(frame, text="File:").grid(row=row,
				column=col, sticky=E)
		col += 1; span = 4
		self.datafile = Entry(frame, foreground="DarkBlue", background="White")
		self.datafile.grid(row=row, column=col, columnspan=span,
				sticky=EW)
		tkExtra.Balloon.set(self.datafile, "USERDUMP binary output file")
		self.datafile.bind('<Control-Button-1>', self.load)

		col += span
		b = Button(frame, image=tkFlair.icons["load"],
				anchor=W,padx=0,pady=0,
				command=self.load)
		b.grid(row=row, column=col, sticky=W)
		self.widget["datafile"] = self.datafile, b

		# --
		col += 1
		self.particle_list = tkExtra.ExListbox(frame,
				selectmode=MULTIPLE,
				selectborderwidth=0,
				background="White",
				takefocus=True,
				height=6,
				exportselection=FALSE)
		self.particle_list.grid(row=row, rowspan=4, column=col, sticky=NSEW)

		col += 1
		sb = Scrollbar(frame, orient=VERTICAL, takefocus=False,
				command=self.particle_list.yview)
		sb.grid(row=row, rowspan=4, column=col, sticky=NS)
		self.particle_list.config(yscrollcommand=sb.set)
		self.widget["@particles"] = self.particle_list #,sb

		# ---
		row, col = row+1, 0
		Label(frame, text="From Event:").grid(row=row,
				column=col, sticky=E)
		col += 1
		var = IntVar()
		var.set(1)
		s = Spinbox(frame, text=var, from_=1, to_=10000000,
				background="White",width=5)
		tkExtra.Balloon.set(s, "Starting event number")
		s.grid(row=row, column=col, sticky=W)
		self.widget["from"] = s
		self.variable["from"] = var

		# ---
		col += 1
		Label(frame, text="# Events:").grid(row=row,
				column=col, sticky=E)
		col += 1
		var = IntVar()
		var.set(0)
		s = Spinbox(frame, text=var, from_=0, to_=1000000,
				background="White",width=7)
		tkExtra.Balloon.set(s, "Total events to load")
		s.grid(row=row, column=col, sticky=W)
		self.widget["n"] = s
		self.variable["n"] = var

		# ===
		row, col = row+1, 0
		Label(frame, text="Emin:").grid(row=row,
				column=col, sticky=E)
		col += 1
		e = tkExtra.FloatEntry(frame, background="White")
		e.grid(row=row, column=col, sticky=W)
		self.widget["emin"] = e

		# ---
		col += 1
		Label(frame, text="Emax:").grid(row=row,
				column=col, sticky=E)
		col += 1
		e = tkExtra.FloatEntry(frame, background="White")
		e.grid(row=row, column=col, sticky=W)
		self.widget["emax"] = e

		frame.grid_columnconfigure(4, weight=2)	# <<<<-------
		frame.grid_columnconfigure(6, weight=1)

		# ============= Info =============
		infoframe = Frame(self)
		infoframe.pack(side=TOP, expand=TRUE, fill=BOTH, pady=3, padx=3)

		# Tabs
		self.tabPages = tkExtra.TabPageSet(infoframe,
				pageNames=Plot.USERDUMP_CASE)
		self.tabPages.changePage()
		self.tabPages.pack(side=LEFT, expand=TRUE, fill=BOTH)
#		self.widget["case"] = self.tabPages

		# ------- Geometry ---------
		frame = LabelFrame(infoframe, text="Geometry ")
		frame.pack(side=RIGHT, fill=Y, padx=3)

		# ---
		row, col = 0, 0
		Label(frame, text="Use:").grid(row=row, column=col, sticky=E)

		col += 1
		self.geo = tkExtra.Combobox(frame, width=10)
		self.geo.insert(END, Plot.GEO_NO)
		self.geo.set(Plot.GEO_NO)
		self.geo.grid(row=row, column=col, columnspan=2, sticky=EW)
		self.widget["geo"] = self.geo

		# ---
		row, col = row+1, 0
		Label(frame, text="Axes:").grid(row=row, column=col, sticky=E)
		col += 1
		self.geo_axes = tkExtra.Combobox(frame)
		self.geo_axes.fill(Plot.COORD)
		self.geo_axes.set(Plot.COORD[0])
		self.geo_axes.grid(row=row, column=col, sticky=EW)
		self.widget["geoaxes"] = self.geo_axes

		# ---
		#row, col = row+1, 0
		#Label(frame, text="Plot:").grid(row=row, column=col, sticky=E)
		#col += 1
		#self.plot_axes = StringVar()
		#self.plot_axes.set(Plot.GP_AXES[0])
		#self.plot_axes_menu = OptionMenu(frame, self.plot_axes,
		#			*Plot.GP_AXES)
		#self.plot_axes_menu.grid(row=row, column=col, sticky=W)

		frame.grid_columnconfigure(1, weight=1)

		# ============= Source =============
		pageframe = self.tabPages[Plot.USERDUMP_CASE[2]]

		# ----------------
		frame = LabelFrame(pageframe, text="Type ")
		frame.pack(side=LEFT,fill=Y)

		self.source_type = StringVar()
		self.source_type.set(Plot.HISTOGRAM)
		b1 = Radiobutton(frame, text=Plot.HISTOGRAM,
				variable=self.source_type, value=Plot.HISTOGRAM,
				command=self.typeChange, anchor=W)
		tkExtra.Balloon.set(b1, "Draw one variable histogram")
		b1.pack(side=TOP,fill=X)

		b2 = Radiobutton(frame, text=Plot.SCATTER,
				variable=self.source_type, value=Plot.SCATTER,
				command=self.typeChange, anchor=W)
		tkExtra.Balloon.set(b2, "Draw a 2D scatter plot")
		b2.pack(fill=X)
		self.widget["type"]   = b1,b2
		self.variable["type"] = self.source_type

		# ------- X ---------
		frame = LabelFrame(pageframe, text="X ")
		frame.pack(side=LEFT, expand=YES, fill=BOTH)

		row, col = 0, 0
		Label(frame, text="Plot:", anchor=W).grid(row=row,
				column=col, sticky=W, pady=0)
		col += 1
		s = tkExtra.Combobox(frame, width=10,
					command=self.histChange)
		s.fill(Plot.SOURCE_HIST)
		s.set(Plot.SOURCE_HIST[0])
		s.grid(row=row, column=col, columnspan=3, sticky=EW)
		self.widget["x"] = s

		# ---
		row, col = row+1, 0
		Label(frame, text="Min:", anchor=W).grid(row=row,
				column=col, sticky=W, pady=0)
		col += 1
		self.hmin = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.hmin.grid(row=row, column=col, sticky=EW)

		col += 1
		Label(frame, text="Max:", anchor=W).grid(row=row,
				column=col, sticky=W, pady=0)
		col += 1
		self.hmax = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.hmax.grid(row=row, column=col, sticky=EW)

		# ----
		row, col = row+1, 0
		Label(frame, text="Bins:", anchor=W).grid(row=row,
				column=col, sticky=W, pady=0)
		col += 1
		self.bins = tkExtra.IntegerEntry(frame, background="White", width=10)
		tkExtra.Balloon.set(self.bins, "Number of bins to create histogram")
		self.bins.grid(row=row, column=col, sticky=EW)
		self.widget["bins"] = self.bins

		# ----
		col += 1
		var = IntVar()
		var.set(1)
		b = Checkbutton(frame, text="Weighted", variable=var)
		tkExtra.Balloon.set(b, "Weight histogram with particle weight")
		b.grid(row=row, column=col, columnspan=2, sticky=W)
		self.widget["weighted"] = b
		self.variable["weighted"] = var

		# Widgets to enable/disable
		self.x_widgets = [self.bins, b]

		frame.grid_columnconfigure(1,weight=1)
		frame.grid_columnconfigure(3,weight=1)

		# -------- Y --------
		frame = LabelFrame(pageframe, text="Y ")
		frame.pack(side=LEFT, expand=YES, fill=BOTH)

		row, col = 0, 0
		Label(frame, text="Plot:", anchor=W).grid(row=row,
				column=col, columnspan=3, sticky=W, pady=0)
		col += 1
		s = tkExtra.Combobox(frame, width=10,
					command=self.histChange)
		s.fill(Plot.SOURCE_HIST)
		s.set(Plot.SOURCE_HIST[0])
		s.grid(row=row, column=col, columnspan=3, sticky=EW)
		self.widget["y"] = s

		# Widgets to enable/disable
		self.y_widgets = [s]

		# ---
		row,col = row+1, 0
		Label(frame, text="Min:", anchor=W).grid(row=row,
				column=col, sticky=W, pady=0)
		col += 1
		self.h2min = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.h2min.grid(row=row, column=col, sticky=EW)

		col += 1
		Label(frame, text="Max:", anchor=W).grid(row=row,
				column=col, sticky=W, pady=0)
		col += 1
		self.h2max = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.h2max.grid(row=row, column=col, sticky=EW)

		frame.grid_columnconfigure(1,weight=1)
		frame.grid_columnconfigure(3,weight=1)

		self.typeChange()

		# --- Particle List ---
		for p in Input.Particle.list:
			if p != "": self.particle_list.insert(END, p)
		self._datafile = None

		self.bindToggle()

	# ----------------------------------------------------------------------
	def type(self): return "USERDUMP"

	# ----------------------------------------------------------------------
	# Show Frame
	# ----------------------------------------------------------------------
	def set(self, plot):
		self.fillGeometry()
		if plot["proj"]=="": plot["proj"] = "Z"
		BasePlot.BasePlotFrame.set(self, plot)
		if self.plot is None: return

		self.tabPages.changePage(self.plot["case"])
		self._datafile = self.plot["datafile"]

		self.particle_list.selection_clear(0,END)
		for p in self.plot["particles"].split():
			self.particle_list.selection_set(
				Input.Particle.list.index(p)-1)

		# Try to read file?
		self.mgfile  = ""
		self.limfile = ""
		self.dmin    = None
		self.dmax    = None
		self.fillGeometry()
		self.typeChange()

	# ----------------------------------------------------------------------
	def get(self, plot):
		BasePlot.BasePlotFrame.get(self, plot)
		if plot is None: return

		# Enable all fields
		tkExtra.multiConfig(self.x_widgets, state=NORMAL)
		tkExtra.multiConfig(self.y_widgets, state=NORMAL)

		partlist = [self.particle_list.get(x)
				for x in self.particle_list.curselection()]
		plot["particles"] = " ".join(partlist)
		plot["case"]      = self.tabPages.getActivePage()

		self.typeChange()

	# ----------------------------------------------------------------------
	def typeChange(self):
		if self.source_type.get() == Plot.HISTOGRAM:
			tkExtra.multiConfig(self.x_widgets, state=NORMAL)
			tkExtra.multiConfig(self.y_widgets, state=DISABLED)
		else:
			tkExtra.multiConfig(self.x_widgets, state=DISABLED)
			tkExtra.multiConfig(self.y_widgets, state=NORMAL)

	# ----------------------------------------------------------------------
	def histChange(self):
		if self.dmin:
			col = Plot.SOURCE_HIST.index(self.widget["x"].get())
			self.hmin.config(text="%10g"%(self.dmin[col]))
			self.hmax.config(text="%10g"%(self.dmax[col]))

			col = Plot.SOURCE_HIST.index(self.widget["y"].get())
			self.h2min.config(text="%10g"%(self.dmin[col]))
			self.h2max.config(text="%10g"%(self.dmax[col]))

	# ----------------------------------------------------------------------
	def fillGeometry(self):
		self.geo.clear()
		self.geo.insert(END, Plot.GEO_NO)

		for plot in self.project.plots:
			if plot.type == "Geometry":
				self.geo.insert(END, plot.name)

	# ----------------------------------------------------------------------
	def load(self, event=None):
		inpname = self.project.inputName
		fn = bFileDialog.askopenfilename(master=self,
			title="Load USERDUMP file",
			initialfile=self.datafile.get(),
			filetypes=[
			        ("%s USERDUMP files"%(inpname),"%s*dump*"%(inpname)),
				("USERDUMP files","*dump*"),
				("All","*") ])
		if len(fn) > 0:
			self._datafile = self.project.relativePath(fn)
			self.datafile.delete(0,END)
			self.datafile.insert(0,self._datafile)
