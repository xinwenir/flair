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
import math
import bmath
import string
import sys
from log import say

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import Data
import Input
import Project
import Layout

import tkFlair
import tkExtra
import bFileDialog
import tkinter as tk

import Plot
import Gnuplot
import BasePlot
import PlotEngine
import GeometryPlot



_ACCURACY = 8
_format = "%.7g"

#-------------------------------------------------------------------------------
# range of floating point numbers
#-------------------------------------------------------------------------------
def binrange(start,stop,step):
	"""formatted number range"""
	if abs(start)<1e-9: start = 0.0
	if abs(stop)<1e-9:   stop = 0.0
	lst = [""]
	try:
		n = int(round((stop-start)/step))
	except (ZeroDivisionError, OverflowError):
		return lst

	# set an upper limit to avoid huge lists
	for i in range(min(n+1,10000)):
		x = start + float(i)*step
		if abs(x)<1e-9:
			lst.append("0.")
		else:
			num = _format%(x)
			if "." not in num and "e" not in num: num += "."
			lst.append(num)
	return lst

#===============================================================================
# USRBIN Data Plot Frame
#===============================================================================
class UsrbinPlotFrame(BasePlot.BasePlotFrame):
	def __init__(self, master, flair):
		BasePlot.BasePlotFrame.__init__(self, master, flair)

		# === File & Detector ===
		frame = LabelFrame(self, text="Binning Detector ")
		frame.pack(side=TOP, expand=YES, fill=X, padx=3)

		row, col = 0, 0
		Label(frame, text="File:").grid(row=row,
				column=col, sticky=E)
		col += 1; span = 3
		self.datafile = Entry(frame, foreground="DarkBlue", background="White")
		tkExtra.Balloon.set(self.datafile, "USRBIN binary output file")
		self.datafile.grid(row=row, column=col, columnspan=span,
				sticky=EW)
		self.datafile.bind('<Control-Button-1>', self.load)
		self.datafile.bind('<Return>',   self.datafileFocusOut)
		self.datafile.bind('<KP_Enter>', self.datafileFocusOut)
		self.datafile.bind('<FocusOut>', self.datafileFocusOut)

		col += span
		b = Button(frame, image=tkFlair.icons["load"],
				anchor=W,padx=0,pady=0,
				command=self.load)
		b.grid(row=row, column=col, sticky=W)
		self.widget["datafile"] = self.datafile, b

		# ---
		col += 1
		Label(frame, text="Title:", anchor=E).grid(row=row,
				column=col, sticky=E)
		col += 1
		self.usr_title = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_title.grid(row=row, column=col, columnspan=3, sticky=EW)

		# ---
		row, col = row+1, 0
		Label(frame, text="Cycles:", anchor=E).grid(row=row,
				column=col, sticky=E)
		col += 1
		self.usr_cycles = Label(frame, anchor=E,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_cycles.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Primaries:", anchor=E).grid(row=row,
				column=col, sticky=E)
		col += 1
		self.usr_ncase = Label(frame, anchor=E,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_ncase.grid(row=row, column=col, columnspan=2, sticky=EW)

		# ---
		col += 2
		Label(frame, text="Weight:", anchor=E).grid(row=row,
				column=col, sticky=E)
		col += 1
		self.usr_weight = Label(frame, anchor=E,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_weight.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Time:", anchor=E).grid(row=row,
				column=col, sticky=E)
		col += 1
		self.usr_time = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_time.grid(row=row, column=col, sticky=EW)

		# ---
		frame.grid_columnconfigure(1, weight=1)	# <<<<-------
		frame.grid_columnconfigure(3, weight=1)
		frame.grid_columnconfigure(6, weight=1)
		frame.grid_columnconfigure(8, weight=5)


		# === Info ===
		frame = LabelFrame(self, text="Binning Info ")
		frame.pack(expand=YES, fill=X, padx=3)

		# ---
		row, col = row+1, 0
		Label(frame, text="Det:", anchor=W).grid(row=row,
				column=col, sticky=W, pady=0)
		col += 1
		self.det = tkExtra.Combobox(frame, width=10, command=self.fillDetectorInfo)
		self.det.grid(row=row, column=col, sticky=EW, pady=0)
		self.widget["@det"] = self.det

		# ---
		col += 1
		self.usr_x = Label(frame, text="X:")
		self.usr_x.grid(row=row, column=col, sticky=E, pady=1)
		col += 1
		self.usr_x_val = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_x_val.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Min:").grid(row=row, column=col, sticky=E)
		col += 1
		self.usr_min = Label(frame, anchor=W,
				foreground=tkFlair._CLABEL_FOREGROUND_COLOR,
				background=tkFlair._CLABEL_BACKGROUND_COLOR)
		self.usr_min.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.usr_min, "Current projection minimum value")

		# ---
		row, col = row+1, 0
		Label(frame, text="Type:", anchor=W).grid(row=row,
				column=col, sticky=W)
		col += 1
		self.usr_type = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_type.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		self.usr_y = Label(frame, text="Y:")
		self.usr_y.grid(row=row,
				column=col, sticky=E, pady=1)
		col += 1
		self.usr_y_val = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_y_val.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Max:").grid(row=row, column=col, sticky=E)
		col += 1
		self.usr_max = Label(frame, anchor=W,
				foreground=tkFlair._CLABEL_FOREGROUND_COLOR,
				background=tkFlair._CLABEL_BACKGROUND_COLOR)
		self.usr_max.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.usr_max, "Current projection maximum value")

		# ---
		row, col = row+1, 0
		Label(frame, text="Score:", anchor=W).grid(row=row,
				column=col, sticky=W)
		col += 1
		self.usr_score = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_score.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		self.usr_z = Label(frame, text="Z:")
		self.usr_z.grid(row=row, column=col, sticky=E, pady=1)
		col += 1
		self.usr_z_val = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_z_val.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Int:").grid(row=row, column=col, sticky=E)
		col += 1
		self.usr_int = Label(frame, anchor=W,
				foreground=tkFlair._CLABEL_FOREGROUND_COLOR,
				background=tkFlair._CLABEL_BACKGROUND_COLOR)
		self.usr_int.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.usr_int, "Current projection integral value")

		# ---
		frame.grid_columnconfigure(1, weight=1)	# <<<<-------
		frame.grid_columnconfigure(3, weight=1)
		frame.grid_columnconfigure(5, weight=1)

		# === MultiFrame: Projection, plotting ... ===
		multiframe = Frame(self)
		multiframe.pack(expand=YES, fill=BOTH)

		# === 1st Projection & Limits ===
		frame = LabelFrame(multiframe, text="Projection & Limits ")
		frame.pack(side=LEFT, expand=YES, fill=X, padx=3)

		# ---
		row, col = 0, 0
		self.proj  = StringVar()
		self.proj.set("Z")
		self.proj_x = Radiobutton(frame, text="X",
				variable=self.proj, value='X')
		self.proj_x.grid(row=row, column=col, sticky=W)

		col += 1
		self.x_min = tkExtra.Combobox(frame, label=True,
				width=8, command=self.xChanged)
		self.x_min._text.config(background="White")
		self.x_min.grid(row=row, column=col, sticky=EW)
		self.widget["bxmin"] = self.x_min

		col += 1
		self.x_rebin = IntVar()
		self.x_rebin.set(1)
		self.updating = False
		self.x_rebin.trace("w", self.xChanged)
		s = Spinbox(frame, text=self.x_rebin, from_=1, to=1000,
				background="White",width=3)
		s.grid(row=row, column=col)
		tkExtra.Balloon.set(s, "Rebin factor")
		self.widget["xrebin"] = s
		self.variable["xrebin"] = self.x_rebin

		col += 1
		self.x_max = tkExtra.Combobox(frame, label=True, width=8)
		self.x_max._text.config(background="White")
		self.x_max.grid(row=row, column=col, sticky=EW)
		self.widget["bxmax"] = self.x_max

		# ---
		col += 1
		b = Button(frame, text="Get", padx=2, pady=0,
				command=self.get_limits)
		b.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(b, "Get limits from plot window")

		# ---
		row, col = row+1, 0
		self.proj_y = Radiobutton(frame, text="Y",
				variable=self.proj, value='Y')
		self.proj_y.grid(row=row, column=col, sticky=W)

		col += 1
		self.y_min = tkExtra.Combobox(frame, label=True,
				width=8, command=self.yChanged)
		self.y_min._text.config(background="White")
		self.y_min.grid(row=row, column=col, sticky=EW)
		self.widget["bymin"] = self.y_min

		col += 1
		self.y_rebin = IntVar()
		self.y_rebin.set(1)
		self.y_rebin.trace("w", self.yChanged)
		s = Spinbox(frame, text=self.y_rebin, from_=1, to=1000,
				background="White",width=3)
		s.grid(row=row, column=col)
		tkExtra.Balloon.set(s, "Rebin factor")
		self.widget["yrebin"] = s
		self.variable["yrebin"] = self.y_rebin

		col += 1
		self.y_max = tkExtra.Combobox(frame, label=True, width=8)
		self.y_max._text.config(background="White")
		self.y_max.grid(row=row, column=col, sticky=EW)
		self.widget["bymax"] = self.y_max

		# ---
		col += 1
		self.swap = IntVar()
		self.swap.set(0)
		b = Checkbutton(frame, text="swap", variable=self.swap)
		b.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(b, "Swap axes on binning data")
		self.widget["swap"] = b
		self.variable["swap"] = self.swap

		# ---
		row, col = row+1, 0

		self.proj_z = Radiobutton(frame, text="Z",
				variable=self.proj, value='Z')
		self.proj_z.grid(row=row, column=col, sticky=W)
		self.widget["proj"]   = (self.proj_x, self.proj_y, self.proj_z)
		self.variable["proj"] = self.proj

		col += 1
		self.z_min = tkExtra.Combobox(frame, label=True,
				width=8, command=self.zChanged)
		self.z_min._text.config(background="White")
		self.z_min.grid(row=row, column=col, sticky=EW)
		self.widget["bzmin"] = self.z_min

		col += 1
		self.z_rebin = IntVar()
		self.z_rebin.set(1)
		self.z_rebin.trace("w", self.zChanged)
		s = Spinbox(frame, text=self.z_rebin, from_=1, to=1000,
				background="White",width=3)
		s.grid(row=row, column=col)
		tkExtra.Balloon.set(s, "Rebin factor")
		self.widget["zrebin"] = s
		self.variable["zrebin"] = self.z_rebin

		col += 1
		self.z_max = tkExtra.Combobox(frame, label=True, width=8)
		self.z_max._text.config(background="White")
		self.z_max.grid(row=row, column=col, sticky=EW)
		self.widget["bzmax"] = self.z_max

		# ---
		col += 1
		var = IntVar()
		var.set(0)
		b = Checkbutton(frame, text="errors", variable=var)
		b.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(b, "Show relative errors")
		self.widget["errors"] = b
		self.variable["errors"] = var

		# ---
		row, col = row+1, 0
		Label(frame, text="Norm:").grid(row=row, column=col, sticky=E)
		col += 1; span=3
		b = Entry(frame, background="White")
		tkExtra.Balloon.set(b, BasePlot._NORM_BALLOON)
		b.grid(row=row, column=col, columnspan=span, sticky=EW)
		self.widget["norm"] = b

		# ---
		tkExtra.Balloon.set(self.x_min, "Project limit minimum")
		tkExtra.Balloon.set(self.y_min, "Project limit minimum")
		tkExtra.Balloon.set(self.z_min, "Project limit minimum")
		tkExtra.Balloon.set(self.x_max, "Project limit maximum")
		tkExtra.Balloon.set(self.y_max, "Project limit maximum")
		tkExtra.Balloon.set(self.z_max, "Project limit maximum")

		frame.grid_columnconfigure(1, weight=1)	# <<<<-------
		frame.grid_columnconfigure(3, weight=1)

		# === 2nd Tabs depending on type of plot ===

		frame = Frame(multiframe)
		frame.pack(side=TOP, expand=YES, fill=BOTH)
		button_frame = Frame(frame)
		button_frame.pack(side=TOP, fill=X)

		self.typeCombo = tkExtra.Combobox(button_frame,
			width=15, command=self.tabPageChange)
		self.typeCombo.pack(side=RIGHT)
		Label(button_frame, text="Type:").pack(side=RIGHT)
		for i in Plot.USRBIN_HIST: self.typeCombo.insert(END, i)
		self.widget["hist"] = self.typeCombo

		tab_frame = Frame(frame, relief=RAISED)
		tab_frame.pack(side=BOTTOM, fill=BOTH, expand=YES)

		self.frame2d = Frame(tab_frame)
		self.frame2d.pack(expand=YES, fill=BOTH)

		# === Geometry ===
		frame = LabelFrame(self.frame2d, text="Geometry ")
		frame.pack(side=LEFT, expand=YES, fill=X, padx=3)

		# ---
		row, col = 0, 0
		Label(frame, text="Use:").grid(row=row, column=col, sticky=E)

		col += 1
		b = tkExtra.Combobox(frame)
		b.fill((Plot.GEO_NO, Plot.GEO_AUTO))
		b.set(Plot.GEO_NO)
		b.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(b, "Geometry plot to superimpose to data")
		self.widget["geo"] = b

		# ---
		row, col = row+1, 0
		Label(frame, text="Pos:").grid(row=row, column=col, sticky=E)
		col += 1
		b = tkExtra.FloatEntry(frame, background="White", width=10)
		b.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(b,
			"Geometry position on the projecting plane to draw when Auto is selected")
		self.widget["geopos"] = b

		# ---
		row, col = row+1, 0
		Label(frame, text="Axes:").grid(row=row, column=col, sticky=E)
		col += 1
		b = tkExtra.Combobox(frame, width=5)
		b.fill(Plot.COORDAUTO)
		b.set(Plot.COORDAUTO[0])
		b.grid(row=row, column=col, sticky=EW)
		self.widget["axes"] = b

		frame.grid_columnconfigure(1, weight=1)

		# === 2nd tab 1D Histogram ===
		self.frame1d = LabelFrame(tab_frame, text="Options ")
		self.styleFrame(self.frame1d).pack(side=LEFT,fill=BOTH)

		# ---- variables -----
		self.inp  = []
		self.usr  = Data.Usrbin()
		self._datafile = None

		self.bindToggle()

	# ----------------------------------------------------------------------
	def type(self): return "USRBIN"

	# ----------------------------------------------------------------------
	# Show Frame
	# ----------------------------------------------------------------------
	def set(self, plot):
		# Set cbtics if needed
		try:
			cb = int(plot.get("cbtics",0))
		except:
			cb = int(plot.get("cbtics","False")=="True")
		if cb==0 and ( \
			plot["cblabel"] or \
			plot["cbmin"]!="" or \
			plot["cbmax"] or \
			plot["hist"]==Plot.USRBIN_HIST[0]):
				plot["cbtics"] = 1

		self.fillGeometry()

		if plot["proj"]=="": plot["proj"] = "Z"

		BasePlot.BasePlotFrame.set(self, plot)
		if plot is None: return

		self.updating = True

		# ---
		self.x_rebin.set(int(plot.get("xrebin",1)))
		self.y_rebin.set(int(plot.get("yrebin",1)))
		self.z_rebin.set(int(plot.get("zrebin",1)))

		# ---
		self._datafile = plot["datafile"]
		self.loadDataFile(self._datafile)
		try:
			d = int(plot("det", 1))
		except:
			d = 1
		self.detector(d)

		try: pt = int(plot["pt.0"]) % 14
		except: pt = 0
		self.plot_pt.set(PlotEngine.POINTS[pt])

		self.updating = False
		self.xChanged()
		self.yChanged()
		self.zChanged()
		self.afterPlot()

	# ----------------------------------------------------------------------
	def get(self, plot):
		BasePlot.BasePlotFrame.get(self, plot)
		if plot is None: return

		self.updating = True
		self.update_idletasks()
		if self.det.cget("state") == NORMAL:
			det = self.detector()
			if det>=0: plot["det"] = str(det)

		if self.plot_pt["state"] == NORMAL:
			plot["pt.0"] = PlotEngine.POINTS.index(self.plot_pt.get())
		self.updating = False

	# ----------------------------------------------------------------------
	def afterPlot(self):
		self.usr_min.config(text=self.plot["min"])
		self.usr_max.config(text=self.plot["max"])
		self.usr_int.config(text=self.plot["int"])

	# --------------------------------------------------------------------
	def xChanged(self, a=None, b=None, c=None): self.fillMaxLimit(0)
	def yChanged(self, a=None, b=None, c=None): self.fillMaxLimit(1)
	def zChanged(self, a=None, b=None, c=None): self.fillMaxLimit(2)

	# --------------------------------------------------------------------
	def fillMaxLimit(self, axis):
		if self.updating: return
		d = self.detector()
		if d<=0 or d>len(self.usr.detector): return
		det = self.usr.detector[d-1]

		if   axis == 0:
			cmin   = self.x_min
			crebin = self.x_rebin
			cmax   = self.x_max
			dlow   = det.xlow
			dhigh  = det.xhigh
			n      = det.nx
			dd     = det.dx
		elif axis == 1:
			cmin   = self.y_min
			crebin = self.y_rebin
			cmax   = self.y_max
			dlow   = det.ylow
			dhigh  = det.yhigh
			n      = det.ny
			dd     = det.dy
		else:
			cmin   = self.z_min
			crebin = self.z_rebin
			cmax   = self.z_max
			dlow   = det.zlow
			dhigh  = det.zhigh
			n      = det.nz
			dd     = det.dz

		try:
			vmin = float(cmin.get())
		except ValueError:
			cmax.fill(binrange(dlow, dhigh, dd))
			return

		try:
			drebin = dd * float(crebin.get())
		except ValueError:
			drebin = dd

		# Round vmin to the next bin
		minbin = round((vmin - dlow)/dd)
		vmin = dlow + (minbin*dd + drebin)

		# Fill max list
		cmax.fill(binrange(vmin, dhigh, drebin))

		# Find closest max
		try:
			vmax   = float(cmax.get())
			maxbin = round((vmax - dlow)/dd)
		except ValueError:
			return

		if maxbin < minbin:
			cmax.set("")
		else:
			# Round vmax to closest re-bin
			cmax.set(_format%(round((vmax-vmin)/drebin)*drebin+vmin))

	# ----------------------------------------------------------------------
	def clearInfo(self, all=True):
		if all:
			self.usr_title.config(text="")
			self.usr_time.config(text="")
			self.usr_weight.config(text="")
			self.usr_ncase.config(text="")
			self.usr_cycles.config(text="")

		self.usr_type.config(text="")
		self.usr_score.config(text="")
		self.usr_x_val.config(text="")
		self.usr_y_val.config(text="")
		self.usr_z_val.config(text="")
		self.usr_min.config(text="")
		self.usr_max.config(text="")
		self.usr_int.config(text="")

	# ----------------------------------------------------------------------
	def loadDataFile(self, datafile):
		# --- Load USRBIN ---
		if len(datafile)==0:
			self._datafile = None
			self.usr.reset()
			self.clearInfo()
			self.det.clear()
			return

		try:
			self.usr.readHeader(datafile)
		except (IOError, MemoryError):
			if self.winfo_ismapped():
				# FIXME shoud display a nice message
				self.flair.notify("Error loading USRBIN",
					"Error loading USRBIN file: %s"%(sys.exc_info()[1]),
					tkFlair.NOTIFY_ERROR)
				say("ERROR:",sys.exc_info()[1])
			self.clearInfo()
			self.det.clear()
			return

		self.usr_title.config(text=self.usr.title)
		self.usr_time.config(text=self.usr.time)
		self.usr_weight.config(text=self.usr.weight)
		self.usr_ncase.config(text=self.usr.ncase)
		self.usr_cycles.config(text=self.usr.nbatch)

		# Populate the detector menu
		self.det.delete(0,END)
		for i,det in enumerate(self.usr.detector):
			self.det.insert(END, "%d %s" % (i+1, det.name))

	# --------------------------------------------------------------------
	def fillDetectorInfo(self):
		# --- Display Detector ---
		d = self.detector()
		if d<=0 or d>len(self.usr.detector):
			self.clearInfo(False)
			return

		det = self.usr.detector[d-1]
		for (n,i) in Layout._USRBIN_type:
			if i==det.type:
				self.usr_type.config(
					text="%d: %s"%(i, n))
				break
		self.usr_score.config(text=\
			Input.Particle.convert(det.score))

		b10 = det.type % 10
		if   b10==1 or b10==7:
			labels = ("R:", "\u03A6:", "Z:")
		elif b10==2:
			labels = ("R1:", "R2:", "R3:")
		elif b10==8:
			labels = ("I1:", "I2:", "I3:")
		else:
			labels = ("X:", "Y:", "Z:")

		self.usr_x.config(text=labels[0])
		self.usr_y.config(text=labels[1])
		self.usr_z.config(text=labels[2])

		self.proj_x.config(text=labels[0])
		self.proj_y.config(text=labels[1])
		self.proj_z.config(text=labels[2])

		self.usr_x_val.config(text="[%g .. %g] x %d (%g)" \
			% (det.xlow, det.xhigh, det.nx, det.dx))
		self.usr_y_val.config(text="[%g .. %g] x %d (%g)" \
			% (det.ylow, det.yhigh, det.ny, det.dy))
		self.usr_z_val.config(text="[%g .. %g] x %d (%g)" \
			% (det.zlow, det.zhigh, det.nz, det.dz))

		# Fill Comboboxes
		self.x_min.fill(binrange(det.xlow, det.xhigh, det.dx))
		self.y_min.fill(binrange(det.ylow, det.yhigh, det.dy))
		self.z_min.fill(binrange(det.zlow, det.zhigh, det.dz))

		self.fillMaxLimit(0)
		self.fillMaxLimit(1)
		self.fillMaxLimit(2)

		self.selectProjection(det)

	# --------------------------------------------------------------------
	# Override the setState to correct for the correct projection
	# --------------------------------------------------------------------
	def setState(self, state):
		BasePlot.BasePlotFrame.setState(self, state)
		if state==NORMAL:
			d = self.detector()
			if d<=0 or d>len(self.usr.detector): return
			self.selectProjection(self.usr.detector[d-1])

	# --------------------------------------------------------------------
	def selectProjection(self, det):
		page = self.typeCombo.get()

		self.proj_x.config(state=NORMAL)
		self.proj_y.config(state=NORMAL)
		self.proj_z.config(state=NORMAL)

		# Check how many bins are requested
		# and active tab
		if page == Plot.USRBIN_HIST[0]:
			if det.type not in (2,12):
				if   det.nx == 1:
					self.proj.set('X')
					self.proj_y.config(state=DISABLED)
					self.proj_z.config(state=DISABLED)
				elif det.ny == 1:
					self.proj.set('Y')
					self.proj_x.config(state=DISABLED)
					self.proj_z.config(state=DISABLED)
				elif det.nz == 1:
					self.proj.set('Z')
					self.proj_x.config(state=DISABLED)
					self.proj_y.config(state=DISABLED)

				ones = (det.nx==1) + (det.ny==1) + (det.nz==1)
				if ones>1:
					self.flair.notify("USRBIN 2D Projection",
						"Two dimensions have ONLY 1 bin,\n" \
						"Therefore 2D projection IS NOT possible ONLY 1D",1)
			else:
				self.proj_x.config(state=DISABLED)
				self.proj_y.config(state=DISABLED)
				self.proj_z.config(state=DISABLED)

		#elif page == Plot.USRBIN_HIST[1]:
		else:
			if   det.ny == 1 and det.nz == 1:
				self.proj.set('X')
			elif det.nx == 1 and det.nz == 1:
				self.proj.set('Y')
			elif det.nx == 1 and det.ny == 1:
				self.proj.set('Z')

			if det.nx == 1:
				self.proj_x.config(state=DISABLED)
				if self.proj.get()=="X": self.proj.set("Y")

			if det.ny == 1:
				self.proj_y.config(state=DISABLED)
				if self.proj.get()=="Y": self.proj.set("Z")

			if det.nz == 1:
				self.proj_z.config(state=DISABLED)
				if self.proj.get()=="Z":
					if det.nx != 1:
						self.proj.set("X")
					else:
						self.proj.set("Y")

	# --------------------------------------------------------------------
	def tabPageChange(self, event=None):
		if self.typeCombo.get() == Plot.USRBIN_HIST[0]:
			self.frame2d.pack(expand=YES, fill=BOTH)
			try: self.frame1d.pack_forget()
			except: pass
		else:
			self.frame1d.pack(expand=YES, fill=BOTH)
			try: self.frame2d.pack_forget()
			except: pass
		d = self.detector()
		if d<=0 or d>len(self.usr.detector): return
		self.selectProjection(self.usr.detector[d-1])

	# ----------------------------------------------------------------------
	def detector(self, val=None):
		if val is None:
			# Get the value
			try:
				return int(self.det.get().split()[0])
			except:
				return -1
		else:
			if 0 < val <= len(self.usr.detector):
				det = self.usr.detector[val-1]
				self.det.set("%d %s" % (val, det.name))
			else:
				self.clearInfo(False)

	# ----------------------------------------------------------------------
	def fillGeometry(self):
		geo = self.widget["geo"]
		geo.clear()
		geo.insert(END, Plot.GEO_NO)
		geo.insert(END, Plot.GEO_AUTO)

		for plot in self.project.plots:
			if plot.type == "Geometry":
				geo.insert(END, plot.name)

	# ----------------------------------------------------------------------
	def load(self, event=None):
		inpname = self.project.inputName
		fn = bFileDialog.askopenfilename(master=self,
			title="Load USRBIN file",
			initialfile=self.datafile.get(),
			filetypes=[
			        ("Fluka USRBIN files",("*usrbin*","*.bnn")),
			        ("%s USRBIN files"%(inpname),"%s*usrbin*"%(inpname)),
				("Fluka fort files",("*_ftn.*","*_fort.*")),
				("All","*")])
		if len(fn) > 0:
			self._datafile = str(self.project.relativePath(fn))
			self.datafile.delete(0,END)
			self.datafile.insert(0,self._datafile)
			self.loadDataFile(self._datafile)
			self.detector(1)
			self.fillDetectorInfo()

	# ----------------------------------------------------------------------
	def datafileFocusOut(self, event=None):
		if self._datafile == self.datafile.get(): return
		self._datafile = self.project.relativePath(self.datafile.get())
		self.datafile.delete(0,END)
		self.datafile.insert(0,self._datafile)
		self.loadDataFile(self._datafile)
		self.detector(1)
		self.fillDetectorInfo()

	# ----------------------------------------------------------------------
	def replot(self):
		self.inp  = []
		self.plot.what = [""] * 24
		BasePlot.BasePlotFrame.replot(self)

	# ----------------------------------------------------------------------
	def get_limits(self):
		if not self.flair.plotEngine: return
		xlow, xhigh = self.flair.plotEngine.getRange("x")
		if xhigh is None: return
		ylow, yhigh = self.flair.plotEngine.getRange("y")
		if yhigh is None: return

		proj = self.proj.get()
		if proj == "X":
			xlow_field  = self.y_min
			xhigh_field = self.y_max
			ylow_field  = self.z_min
			yhigh_field = self.z_max

		elif proj == "Y":
			xlow_field  = self.x_min
			xhigh_field = self.x_max
			ylow_field  = self.z_min
			yhigh_field = self.z_max

		elif proj == "Z":
			xlow_field  = self.x_min
			xhigh_field = self.x_max
			ylow_field  = self.y_min
			yhigh_field = self.y_max

		if self.swap.get():
			xlow_field,  ylow_field  = ylow_field,  xlow_field
			xhigh_field, yhigh_field = yhigh_field, xhigh_field

		xlow_field.set( _format%(xlow))
		xhigh_field.set(_format%(xhigh))
		ylow_field.set( _format%(ylow))
		yhigh_field.set(_format%(yhigh))

#===============================================================================
if __name__ == "__main__":
	tk = Tk()
	f = UsrbinPlotFrame(tk)
	f.pack(expand=YES,fill=BOTH)
	tk.mainloop()
