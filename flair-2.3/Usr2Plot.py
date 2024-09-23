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
# Date:	02-Jul-2008

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os

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

import BasePlot

#===============================================================================
# Double Differential Plot Frame (mainly for USRBDX)
#===============================================================================
class Usr2DPlotFrame(BasePlot.BasePlotFrame):
	def __init__(self, master, flair):
		BasePlot.BasePlotFrame.__init__(self, master, flair)

		# === File & Detector ===
		lblframe = tkExtra.ExLabelFrame(self, text="Detector ")
		lblframe.pack(side=TOP, expand=YES, fill=BOTH, padx=3)
		frame = lblframe.frame

		row, col = 0, 0
		Label(frame, text="File:", anchor=E).grid(row=row,
				column=col, sticky=E)
		col += 1
		self.datafile = Entry(frame, foreground="DarkBlue", background="White")
		tkExtra.Balloon.set(self.datafile, "USRBDX tab.lis output file")
		self.datafile.grid(row=row, column=col, sticky=EW)
		self.datafile.bind('<Control-Button-1>', self.load)
		self.datafile.bind('<Return>',   self.datafileFocusOut)
		self.datafile.bind('<KP_Enter>', self.datafileFocusOut)
		self.datafile.bind('<FocusOut>', self.datafileFocusOut)

		col += 1
		b = Button(frame, image=tkFlair.icons["load"],
				anchor=W,padx=0,pady=0,
				command=self.load)
		b.grid(row=row, column=col, sticky=W)
		self.widget["datafile"] = self.datafile, b

		col += 1
		Label(frame, text="Det:", anchor=W).grid(row=row,
				column=col, sticky=W, pady=0)
		col += 1
		self.det = tkExtra.Combobox(frame, width=10)
		self.det.grid(row=row, column=col, sticky=EW, pady=0)
		self.widget["@det"] = self.det

		frame.grid_columnconfigure(1, weight=1)
		frame.grid_columnconfigure(4, weight=1)

		self._datafile = None
		self.bindToggle()

	# ----------------------------------------------------------------------
	def type(self): return "USR-2D"

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

	# ----------------------------------------------------------------------
	def get(self, plot):
		BasePlot.BasePlotFrame.get(self, plot)
		if plot is None: return

		if self.det.cget("state") == NORMAL:
			det = self.detector()
			if det>=0: plot["det"] = str(det)

	# ----------------------------------------------------------------------
	def detector(self, val=None):
		if val is None:
			try:
				# Return the value
				return int(self.det.get().split()[0])-1
			except:
				return 0
		else:
			self.det.set("")
			if val >= 0:
				val1 = val + 1
				for item in self.det.get(0,END):
					if int(item.split()[0])==val1:
						self.det.set(item)

	# ----------------------------------------------------------------------
	def loadDataFile(self, datafile):
		# --- Load a _tab.lis file ---
		if len(datafile)==0:
			self._datafile = None
			return

		# Load detector file
		self.det.clear()
		cnt   = 1
		first = None
		half  = 0
		name  = ""
		try:
			f = open(datafile,"r")
			for line in f:
				if line=="\n":
					half += 1
					if half == 2:
						half = 0
						cnt += 1

				elif line.startswith(" # Detector n:"):
					name = line[17:]
					p = name.find("(")
					if p>0: name = name[:p]
					name = "%d %s"%(cnt+1,name.strip())

				elif line.startswith(" # double differential"):
					self.det.insert(END, name)
					if not first:
						first = name

				else:
					half = 0
			f.close()
			if first:
				self.det.set(first)
		except IOError:
			pass

	# ----------------------------------------------------------------------
	def load(self, event=None):
		inpname = self.project.inputName
		fn = bFileDialog.askopenfilename(master=self,
			title="Load USRxxx file",
			filetypes=[
				("Fluka _tab.lis files","*_tab.lis"),
			        ("%s _tab.lis files"%(inpname),"%s*_tab.lis"%(inpname)),
				("All","*")])
		if len(fn) > 0:
			self._datafile = self.project.relativePath(fn)
			self.datafile.delete(0,END)
			self.datafile.insert(0,self._datafile)
			self.loadDataFile(self._datafile)

	# ----------------------------------------------------------------------
	def datafileFocusOut(self, event=None):
		if self._datafile == self.datafile.get(): return
		self._datafile = self.project.relativePath(self.datafile.get())
		self.datafile.delete(0,END)
		self.datafile.insert(0,self._datafile)
		self.loadDataFile(self._datafile)
