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
# Authors:	Vasilis.Vlachoudis@cern.ch,
#		Adonai.Herrera.Martinez@cern.ch
# Date:	31-Oct-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
from log import say

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import Data
import Input
import Layout

import tkFlair
import tkExtra
import bFileDialog

import Plot
import BasePlot

#===============================================================================
# Residual Nuclei Plot Frame
#===============================================================================
class ResnucleiPlotFrame(BasePlot.BasePlotFrame):
	def __init__(self, master, flair):
		BasePlot.BasePlotFrame.__init__(self, master, flair)

		# === File & Detector ===
		frame = LabelFrame(self, text="Residuals Detector ")
		frame.pack(side=TOP, expand=YES, fill=BOTH, padx=3)

		row, col = 0, 0
		Label(frame, text="File:", anchor=E).grid(row=row,
				column=col, sticky=E)
		col += 1; span = 3
		self.datafile = Entry(frame, foreground="DarkBlue", background="White")
		tkExtra.Balloon.set(self.datafile, "RESNUCLEi binary output file")
		self.datafile.grid(row=row, column=col, columnspan=span,
				sticky=EW)
		self.datafile.bind('<Control-Button-1>', self.load)
		self.datafile.bind('<Return>',   self.datafileFocusOut)
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
		frame = LabelFrame(self, text="Detector Info ")
		frame.pack(expand=YES, fill=BOTH, padx=3)

		# ---
		row, col = 0, 0
		Label(frame, text="Det:", anchor=W).grid(row=row,
				column=col, sticky=W, pady=0)
		col += 1
		self.det = tkExtra.Combobox(frame, width=10, command=self.fillDetectorInfo)
		self.det.grid(row=row, column=col, sticky=EW, pady=0)
		self.widget["@det"] = self.det

		# ---
		col += 1
		Label(frame, text="Type:", anchor=W).grid(row=row,
				column=col, sticky=W)
		col += 1
		self.usr_type = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_type.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Region:", anchor=W).grid(row=row,
				column=col, sticky=W)
		col += 1
		self.usr_region = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_region.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Volume:", anchor=W).grid(row=row,
				column=col, sticky=W)
		col += 1
		self.usr_volume = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_volume.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Min:", anchor=E).grid(row=row,
				column=col, sticky=E)
		col += 1
		self.usr_min = Label(frame, anchor=E, width=10,
				foreground=tkFlair._CLABEL_FOREGROUND_COLOR,
				background=tkFlair._CLABEL_BACKGROUND_COLOR)
		self.usr_min.grid(row=row, column=col, sticky=EW)

		# ---
		row, col = row+1, 0
		Label(frame, text="Plot:", anchor=W).grid(row=row,
				column=col, sticky=W)

		col += 1
		self.plotType = tkExtra.Combobox(frame, width=5)
		self.plotType.fill(Plot.RESNUCLEI_TYPES)
		self.plotType.set(Plot.RESNUCLEI_TYPES[0])
		self.plotType.grid(row=row, column=col, sticky=EW)
		self.widget["plot"] = self.plotType

		# ---
		col += 1
		Label(frame, text="Zhigh:", anchor=W).grid(row=row,
				column=col, sticky=W)
		col += 1
		self.usr_zhigh = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_zhigh.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Ahigh:", anchor=W).grid(row=row,
				column=col, sticky=W)
		col += 1
		self.usr_ahigh = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_ahigh.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Mhigh:", anchor=W).grid(row=row,
				column=col, sticky=W)
		col += 1
		self.usr_mhigh = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.usr_mhigh.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(frame, text="Max:", anchor=E).grid(row=row,
				column=col, sticky=E)
		col += 1
		self.usr_max = Label(frame, anchor=E, width=10,
				foreground=tkFlair._CLABEL_FOREGROUND_COLOR,
				background=tkFlair._CLABEL_BACKGROUND_COLOR)
		self.usr_max.grid(row=row, column=col, sticky=EW)

		# ---
		row, col = row+1, 0
		Label(frame, text="Norm:", anchor=W).grid(row=row,
				column=col, sticky=W)

		col += 1; span=3
		b = Entry(frame, background="White")
		tkExtra.Balloon.set(b, BasePlot._NORM_BALLOON)
		b.grid(row=row, column=col, columnspan=span, sticky=EW)
		self.widget["norm"] = b

		# --- Options ---
		self._styleRow = row+1
		self.style = LabelFrame(frame, text="Options ")
		self.styleFrame(self.style)
		self.plotType.command = self.typeChange

		# ---
		frame.grid_columnconfigure(1, weight=1)	# <<<<-------
		frame.grid_columnconfigure(3, weight=1)
		frame.grid_columnconfigure(5, weight=1)
		frame.grid_columnconfigure(7, weight=1)

		# ---- variables -----
		self.usr = Data.Resnuclei()
		self.previous = []
		self._datafile = None

		self.bindToggle()

	# ----------------------------------------------------------------------
	def typeChange(self):
		if self.plotType.get() in (Plot.RESNUCLEI_Z, Plot.RESNUCLEI_A):
			# Show style frame
			self.style.grid(row=self._styleRow, column=0, columnspan=10, sticky=NSEW)
		else:
			# Hide style frame
			self.style.grid_forget()

	# ----------------------------------------------------------------------
	def type(self): return "RESNUCLE"

	# ----------------------------------------------------------------------
	# Set plot
	# ----------------------------------------------------------------------
	def set(self, plot):
		BasePlot.BasePlotFrame.set(self, plot)
		if plot is None: return

		self._datafile = plot["datafile"]
		self.loadDataFile(self._datafile)
		try: det = int(plot("det",1))
		except: det = 1
		self.detector(det)
		self.fillDetectorInfo()

	# ----------------------------------------------------------------------
	def get(self, plot):
		BasePlot.BasePlotFrame.get(self, plot)
		if plot is None: return

		if self.det.cget("state") == NORMAL:
			det = self.detector()
			if det>=0: plot["det"] = str(det)

	# ----------------------------------------------------------------------
	def clearInfo(self, all=True):
		if all:
			self.usr_title.config(text="")
			self.usr_time.config(text="")
			self.usr_weight.config(text="")
			self.usr_ncase.config(text="")
			self.usr_cycles.config(text="")

		self.usr_type.config(text="")
		self.usr_region.config(text="")
		self.usr_volume.config(text="")
		self.usr_zhigh.config(text="")
		self.usr_ahigh.config(text="")
		self.usr_mhigh.config(text="")

		self.usr_min.config(text="")
		self.usr_max.config(text="")

	# ----------------------------------------------------------------------
	def loadDataFile(self, datafile):
		# --- Load RESNUCLEi ---
		if len(datafile)==0:
			self._datafile = None
			self.usr.reset()
			self.clearInfo()
			self.det.clear()
			return

		#if self.usr.file!=datafile or self.usr.ncase==0:
		try:
			self.usr.readHeader(datafile)
		except IOError:
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
		for det in self.usr.detector:
			self.det.insert(END,"%d %s" % (det.nb, det.name))

	# --------------------------------------------------------------------
	def fillDetectorInfo(self):
		# --- Display Detector ---
		d = self.detector()
		if d<=0:
			self.clearInfo(False)
			return

		det = self.usr.detector[d-1]
		for (n,i) in Layout._RESNUCLE_type:
			if i==det.type:
				self.usr_type.config(
					text="%d: %s"%(i, n))
				break
		self.usr_region.config(text=det.region)
		self.usr_volume.config(text="%-11g"%(det.volume))
		self.usr_zhigh.config(text=det.zhigh)
		self.usr_ahigh.config(text=str(2*det.zhigh + det.mhigh + det.nmzmin))
		self.usr_mhigh.config(text=det.mhigh)

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
				self.det.set("%d %s" % (det.nb, det.name))
			else:
				self.clearInfo(False)

	# ----------------------------------------------------------------------
	def load(self, event=None):
		inpname = self.project.inputName
		fn = bFileDialog.askopenfilename(master=self,
			title="Load Residual Nuclei file",
			filetypes=[
				("Fluka RESNUCLEi files",("*resnucle*","*.rnc")),
				("%s RESNUCLEi files"%(inpname),"%s*resnucle*"%(inpname)),
				("All","*")])
		if len(fn) > 0:
			self._datafile = self.project.relativePath(fn)
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

