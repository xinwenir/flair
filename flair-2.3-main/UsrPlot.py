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
# Authors:	Vasilis.Vlachoudis@cern.ch
# Date:	31-Oct-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import re
try:
	from cStringIO import StringIO
except ImportError:
	from io import StringIO
try:
	import cPickle as pickle
except ImportError:
	import pickle

try:
	from Tkinter import *
except ImportError:
	from tkinter import *

import Data
import Input
import Layout

import Plot
import tkFlair
import tkExtra
import bFileDialog

import Gnuplot
import BasePlot
import PlotEngine

_X_OPTION_DEF = 2
_Y_OPTION_DEF = 0

_detectorPattern = re.compile(r"^ ?# ?Detector ?n?:\s*\d*\s*(.*)\s*", re.MULTILINE)
_blockPattern    = re.compile(r"^ ?# ?Block ?n?:\s*\d*\s*(.*)\s*", re.MULTILINE)

_DETECTOR_CLIP = tkFlair._FLAIRCLIP + "<usr1d-detector>"

#===============================================================================
# Residual Nuclei Plot Frame
#===============================================================================
class Usr1DPlotFrame(BasePlot.BasePlotFrame):
	def __init__(self, master, flair):
		BasePlot.BasePlotFrame.__init__(self, master, flair)

		# === Detectors ===
		frame = LabelFrame(self, text="Detectors ")
		frame.pack(side=LEFT, expand=YES, fill=BOTH, padx=3)

		row, col = 0,0
		rowspan  = 7
		self.listbox = tkExtra.ExListbox(frame, borderwidth=0,
				selectmode=EXTENDED,
				selectborderwidth=0,
				background="White",
				takefocus=True,
				width=15,
				height=10,
				exportselection=FALSE)
		self.listbox.grid(row=row, rowspan=rowspan, column=col, sticky=NSEW)
		self.widget["listbox"] = self.listbox

		col += 1
		sbv = Scrollbar(frame, orient=VERTICAL, takefocus=False,
				command=self.listbox.yview)
		sbv.grid(row=row, rowspan=rowspan, column=col, sticky=NS)
		self.listbox.config(yscrollcommand=sbv.set)

		col += 1
		b = Button(frame, image=tkFlair.icons["rename"], takefocus=False,
				command=self.rename)
		b.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(b, "Rename detector")
		self.widget["rename"] = b

		row += 1
		b = Button(frame, image=tkFlair.icons["add"], takefocus=False,
				command=self.add)
		b.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(b, "Add detector")
		self.widget["add"] = b

		row += 1
		b = Button(frame, image=tkFlair.icons["del"], takefocus=False,
				command=self.delete)
		b.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(b, "Delete detector")
		self.widget["del"] = b

		row += 1
		b = Button(frame, image=tkFlair.icons["clone"], takefocus=False,
				command=self.clone)
		b.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(b, "Clone detector")
		self.widget["clone"] = b

		row += 1
		b = Button(frame, image=tkFlair.icons["up"], takefocus=False,
				command=self.moveUp)
		tkExtra.Balloon.set(b, "Move Up")
		b.grid(row=row, column=col, sticky=EW)
		self.widget["up"] = b

		row += 1
		b = Button(frame, image=tkFlair.icons["down"], takefocus=False,
				command=self.moveDown)
		tkExtra.Balloon.set(b, "Move Down")
		b.grid(row=row, column=col, sticky=EW)
		self.widget["down"] = b

		frame.grid_columnconfigure(0, weight=1)	# <<<<--------
		frame.grid_rowconfigure(row+1, weight=1)	# <<<<--------

		# === Detectors ===
		labelframe = LabelFrame(self, text="Detector Info ")
		labelframe.pack(side=LEFT, expand=YES, fill=X, padx=3)

		# ---
		frame = Frame(labelframe)
		frame.pack(side=TOP, fill=X)

		Label(frame, text="File:", anchor=E).pack(side=LEFT)

		self.datafile = Entry(frame, foreground="DarkBlue", background="White")
		self.datafile.pack(side=LEFT, expand=YES, fill=X)
		tkExtra.Balloon.set(self.datafile, "Tabulated _tab.lis output file")
		self.datafile.bind('<Control-Button-1>', self.load)
		self.datafile.bind('<Return>',   self.datafileFocusOut)
		self.datafile.bind('<FocusOut>', self.datafileFocusOut)

		self.datafile_button = Button(frame, image=tkFlair.icons["load"],
				anchor=W,padx=0,pady=0,
				command=self.load)
		self.datafile_button.pack(side=LEFT)
		self.widget["@datafile"] = self.datafile, self.datafile_button

		# ---
		Label(frame, text=" Det:").pack(side=LEFT)
		self.det = tkExtra.Combobox(frame, width=32)
		self.det.pack(side=LEFT, fill=X, expand=YES)
		tkExtra.Balloon.set(self.det, "Select detector from the data file to plot")
		self.widget["@det"] = self.det

		# ====================================
		frame = Frame(labelframe)
		frame.pack(fill=X)

		# ====================================
		lf = LabelFrame(frame, text="Show ")
		lf.pack(side=LEFT)

		self.show_plot = IntVar()
		self.show_plot.set(1)
		self.show_plot_button = Checkbutton(lf, text="graph",
			variable=self.show_plot, anchor=E, justify=RIGHT)
		self.show_plot_button.pack(side=TOP)
		tkExtra.Balloon.set(self.show_plot_button,
				"Include detector/graph in plot")
		self.widget["@show_plot"] = self.show_plot_button

		# ---
		self.show_key = IntVar()
		self.show_key.set(1)
		self.show_key_button = Checkbutton(lf, text="legend",
			variable=self.show_key, anchor=E, justify=RIGHT)
		self.show_key_button.pack(side=TOP)
		tkExtra.Balloon.set(self.show_key_button,
				"Show key/legend of detector")
		self.widget["@show_key"] = self.show_key_button

		# ---
		lf = LabelFrame(frame, text="Plot ")
		lf.pack(side=RIGHT, expand=YES, fill=X)

		# ===
		row,col = 0,0
		Label(lf, text="Type:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_with = tkExtra.Combobox(lf, width=10)
		self.plot_with.fill(Plot.STYLE)
		self.plot_with.set(Plot.STYLE[Plot.STYLE_DEF])
		self.plot_with.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.plot_with, "Plotting style")
		self.widget["@plot_with"] = self.plot_with

		# ---
		col += 1
		Label(lf, text=" X Norm:").grid(row=row, column=col, sticky=E)
		col += 1
		self.x_norm = tkExtra.Combobox(lf, label=False, width=10, background="White")
		self.x_norm.fill(Plot.X_SCALE)
		self.x_norm.set(Plot.X_SCALE[0])
		self.x_norm.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.x_norm, "X axis normalization factor or expression (eg. 3*x-10)")
		self.widget["@x_norm"] = self.x_norm

		# ===
		row,col = row+1,0

		Label(lf, text="Value:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_value = tkExtra.Combobox(lf, width=10)
		self.plot_value.fill(Plot.VALUE)
		self.plot_value.set(Plot.VALUE[0])
		self.plot_value.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.plot_value, "Plotting Y value")
		self.widget["@plot_value"] = self.plot_value

		# ---
		col += 1
		Label(lf, text=" Y Norm:").grid(row=row, column=col, sticky=E)
		col += 1
		self.norm = Entry(lf, background="White", justify=LEFT)
		tkExtra.Balloon.set(self.norm, "Normalization factor or expression (eg. x**2+10)")
		self.norm.grid(row=row, column=col, sticky=EW)
		self.widget["@norm"] = self.norm

		lf.grid_columnconfigure(3, weight=1)

		# ====================================
		frame = Frame(labelframe)
		frame.pack(fill=X)

		# ====================================
		lf = LabelFrame(frame, text="Options ")
		lf.pack(side=LEFT)

		row,col = 0,0
		Label(lf, text="Color:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_lc = tkExtra.Combobox(lf, width=15)
		self.plot_lc.fill(PlotEngine._COLOR_NAMES)
		self.plot_lc.grid(row=row,column=col,sticky=EW)
		tkExtra.Balloon.set(self.plot_lc, "Line color if gnuplot>=4.2")
		self.widget["@plot_lc"] = self.plot_lc

		# ---
		col += 1
		Label(lf, text=" Line width:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_lw = IntVar()
		self.plot_lw.set(1)
		self.plot_lw_spin = Spinbox(lf, text=self.plot_lw, from_=1, to=6,
				background="White",width=3)
		self.plot_lw_spin.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(self.plot_lw_spin, "Line width")
		self.widget["@plot_lw"] = self.plot_lw_spin

		# ===
		row, col = row+1, 0
		Label(lf, text=" Point type:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_pt = tkExtra.Combobox(lf, width=15)
		self.plot_pt.fill(PlotEngine.POINTS)
		self.plot_pt.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(self.plot_pt, "Point type to use")
		self.widget["@plot_pt"] = self.plot_pt

		# ---
		col += 1
		Label(lf, text="Point size:").grid(row=row, column=col, sticky=E)
		col += 1
		self.plot_ps = IntVar()
		self.plot_ps.set(1)
		self.plot_ps_spin = Spinbox(lf, text=self.plot_ps, from_=1, to=6,
				background="White",width=3)
		self.plot_ps_spin.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(self.plot_ps_spin, "Point size")
		self.widget["@plot_ps"] = self.plot_ps_spin

		# ---- bindings -----
		self.listbox.bind('<Return>',          self.rename)
		self.listbox.bind('<KP_Enter>',        self.rename)
		self.listbox.bind('<F2>',              self.rename)
		self.listbox.bind('<Control-Key-e>',   self.rename)
		self.listbox.bind('<<Cut>>',           self.cut)
		self.listbox.bind('<<Copy>>',          self.copy)
		self.listbox.bind('<<Paste>>',         self.paste)
		self.listbox.bind('<Double-1>',        self.rename)
		self.listbox.bind('<<ListboxSelect>>', self.selectDetector)
		self.listbox.bind('<Control-Key-Up>',  self.moveUp)
		self.listbox.bind('<Control-Key-Down>',self.moveDown)

		# ---- variables -----
		self.current = None		# Current detector
		self.selectDetector(DISABLED)
		self._datafile = None		# Remember last data file shown

		self.bindToggle()

	# ----------------------------------------------------------------------
	def type(self): return "USR-1D"

	# ----------------------------------------------------------------------
	# Show Frame
	# ----------------------------------------------------------------------
	def set(self, plot):
		BasePlot.BasePlotFrame.set(self, plot)

		self.current = None
		self.clearInfo()

		self.listbox.delete(0,END)
		for i in range(int(plot.var.get("ndetectors",0))):
			self.listbox.insert(END, plot.get("name.%d"%(i),""))
		self.listbox.selection_set(0)
		self.selectDetector()

	# ----------------------------------------------------------------------
	def get(self, plot):
		BasePlot.BasePlotFrame.get(self, plot)
		if not self.current: return
		plot["ndetectors"] = str(self.listbox.size())
		for cur in self.current:
			detname = self.listbox.get(cur)
			if self.listbox["state"] == NORMAL:
				plot["name.%d"%(cur)] = detname

			if self.datafile["state"] == NORMAL:
				plot["file.%d"%(cur)] = self.datafile.get()

			if self.det._text["state"] == NORMAL:
				det,block = self.detector()
				plot["det.%d"%(cur)] = str(det)
				plot["block.%d"%(cur)]= str(block)

			if self.show_plot_button["state"] == NORMAL:
				plot["show.%d"%(cur)] = str(self.show_plot.get())

			if self.show_key_button["state"] == NORMAL:
				plot["key.%d"%(cur)]  = str(self.show_key.get())

			if self.norm["state"] == NORMAL:
				plot["norm.%d"%(cur)] = self.norm.get()

			if self.plot_with["state"] == NORMAL:
				plot["with.%d"%(cur)]  = self.plot_with.get()

			if self.plot_value["state"] == NORMAL:
				plot["y.%d"%(cur)] = Plot.VALUE.index(self.plot_value.get())

			if self.x_norm["state"] == NORMAL:
				plot["xnorm.%d"%(cur)] = self.x_norm.get()

			if self.plot_lc["state"] == NORMAL:
				plot["lc.%d"%(cur)] = self.plot_lc.get()

			if self.plot_pt["state"] == NORMAL:
				try: idx = PlotEngine.POINTS.index(self.plot_pt.get())
				except: idx = 0
				plot["pt.%d"%(cur)] = idx

			if self.plot_lw_spin["state"] == NORMAL:
				plot["lw.%d"%(cur)] = str(self.plot_lw.get())

			if self.plot_ps_spin["state"] == NORMAL:
				plot["ps.%d"%(cur)] = str(self.plot_ps.get())

	# ----------------------------------------------------------------------
	def selectDetector(self, event=None):
		self.get(self.plot)

#		self.current = map(int, self.listbox.curselection())
		self.current = list(map(int, self.listbox.curselection()))
		if len(self.current)!=1:
			state = DISABLED
		else:
			state = NORMAL

		self.datafile.config(state=state)
		self.datafile_button.config(state=state)
		self.det.config(state=state)

		self.show_plot_button.config(state=state)
		self.show_key_button.config(state=state)

		self.plot_with.config(state=state)
		self.plot_value.config(state=state)
		self.norm.config(state=state)
		self.x_norm.config(state=state)
		self.plot_lc.config(state=state)
		self.plot_pt.config(state=state)
		self.plot_lw_spin.config(state=state)
		self.plot_ps_spin.config(state=state)

		self.clearInfo()
		self.fillDetectorInfo()

	# ----------------------------------------------------------------------
	def fillDetectorInfo(self):
		if not self.current: return

		cur = self.current[0]
		if len(self.current)>1:
			# initial values
			Idatafile = self.plot.get("file.%d"%(cur),"")

			for cur in self.current:
				if Idatafile != self.plot.get("file.%d"%(cur)):
					# Fill detector only with numbers!
					for i in range(1,100):
						self.det.insert(END, "%d ?"%(i))
					break
			else:
				self.loadDataFile(Idatafile)

		else:
			datafile = self.plot.get("file.%d"%(cur),"")
			if datafile != "":
				self.loadDataFile(datafile)

			self.show_plot.set(int(self.plot.get("show.%d"%(cur),1)))
			self.show_key.set(int(self.plot.get("key.%d"%(cur),1)))

			style = self.plot["with.%d"%(cur)]
			if style not in Plot.STYLE: style = "histogram"	# Switch to histogram if not existing
			self.plot_with.set(style)

			try:
				yval = int(self.plot["y.%d"%(cur)])
				if yval: yval = 1	# Isolethargic
			except:
				yval = 0
			self.plot_value.set(Plot.VALUE[yval])

			self.norm.delete(0,END)
			self.norm.insert(0,self.plot["norm.%d"%(cur)])

			self.x_norm.set(self.plot["xnorm.%d"%(cur)])

			self.plot_lc.set(self.plot["lc.%d"%(cur)])
			try:
				pt = int(self.plot["pt.%d"%(cur)]) % 14
			except:
				pt = 0
			self.plot_pt.set(PlotEngine.POINTS[pt])

			try:
				self.plot_lw.set(int(self.plot["lw.%d"%(cur)]))
			except:
				self.plot_lw.set(1)
			try:
				self.plot_ps.set(int(self.plot["ps.%d"%(cur)]))
			except:
				self.plot_ps.set(1)

		self.detector(int(self.plot.get("det.%d"%(cur),0)),
			      int(self.plot.get("block.%d"%(cur),0)))

	# ----------------------------------------------------------------------
	def clearInfo(self):
		self._datafile = None
		self.datafile.delete(0,END)
		self.det.clear()

	# ----------------------------------------------------------------------
	def loadDataFile(self, datafile):
		# --- Load a _tab.lis file ---
		if len(datafile)==0: return

		if datafile == self._datafile: return
		self._datafile = datafile
		self.datafile.delete(0,END)
		self.datafile.insert(0,datafile)

		# Load detector file
		self.det.clear()
		ind   = 1
		first = None
		half  = 0	# half 2xempty lines = 1-line :)
		block = 1
		name = "Detector"
		blockname = ""
		blockspresent = False
		lastind = 0
		try:
			f = open(datafile,"r")
			for line in f:
				if line=="\n":
					half += 1
					if half == 2:
						half = 0
						ind += 1

				elif line.find("#")>=0:
					if line.find("# Detector")>=0:
						m = _detectorPattern.match(line)
						if m:
							name = m.group(1)
						else:
							name = "-unknown-"
						p = name.find("(")
						if p>0: name = name[:p]
						name = name.strip()
						entry = "%d %s"%(ind, name)
						self.det.insert(END, entry)
						lastind = ind
						if not first:
							first = entry
						half  = 0
						block = 0
						blockspresent = False

					elif line.find("# Block")>=0:
						m = _blockPattern.match(line)
						if m:
							blockname = m.group(1)
							blockspresent = True
						else:
							blockname = ""
						half = 1

					else:
						if lastind != ind:
							self.det.insert(END, "%d %s"%(ind,line[1:].strip()))
							lastind = ind

				else:
					if half == 1:
						block += 1
						if blockname != "":
							self.det.insert(END, "%d-%d %s %s"%(ind,block,name,blockname))
							lastind = ind
						elif not blockspresent:
							self.det.insert(END, "%d-%d %s ?"%(ind,block,name))
							lastind = ind
						blockname = ""
					half = 0
			f.close()
			if first:
				self.det.set(first)
		except IOError:
			pass

	# ----------------------------------------------------------------------
	def load(self, event=None):
		if self.datafile["state"] != NORMAL: return
		inpname = self.project.inputName
		fn = bFileDialog.askopenfilename(master=self,
			title="Load USRxxx file",
			initialfile=self.datafile.get(),
			filetypes=[
				("Fluka _tab.lis files","*_tab.lis"),
			        ("%s _tab.lis files"%(inpname),"%s*_tab.lis"%(inpname)),
				("Fluka .lis files","*.lis"),
			        ("Data files","*.dat"),
				("All","*")])
		if len(fn) > 0:
			self.loadDataFile(self.project.relativePath(fn))

	# ----------------------------------------------------------------------
	def datafileFocusOut(self, event=None):
		self.loadDataFile(self.project.relativePath(self.datafile.get()))

	# ----------------------------------------------------------------------
	# Set/Get detector=d/block=b
	# ----------------------------------------------------------------------
	def detector(self, d=None, b=None):
		if d is None:
			d = b = 0
			try:
				db = self.det.get().split()[0]
				if "-" in db:
					d,b = map(int,db.split("-"))
				else:
					d = int(db)
				d -= 1
			except:
				return 0,0
			return d,b
		else:
			self.det.set("")
			if d < 0: return
			d += 1

			for item in self.det.get(0,END):
				db = item.split()[0]
				di = bi = 0
				if "-" in db:
					di,bi = map(int,db.split("-"))
				else:
					di = int(db)
				if d==di and b==bi:
					self.det.set(item)

	# ----------------------------------------------------------------------
	# Return active block of the detector
	# ----------------------------------------------------------------------
	def block(self):
		db = self.det.get()
		if "-" in db:
			return int(db.split("-")[1])
		else:
			return 0

	# ----------------------------------------------------------------------
	def add(self):
		self.get(self.plot)
		ndets = int(self.plot.get("ndetectors",0))
		self.plot["ndetectors"] = ndets+1
		name = "#Detector %d" % (ndets+1)
		self.plot["name.%d"%(ndets)] = name
		self.plot["det.%d"%(ndets)] = 0
		self.listbox.insert(END,name)
		self.listbox.selection_clear(0,END)
		self.listbox.selection_set(END)
		self.listbox.activate(END)
		self.listbox.see(END)

		self.current = None
		self.selectDetector()
		self.load()

	# ----------------------------------------------------------------------
	def rename(self, event=None):
		edit = tkExtra.InPlaceEdit(self.listbox)

	# ----------------------------------------------------------------------
	# find all vars
	# create dynamically the list of variables controlling the detector
	# ----------------------------------------------------------------------
	def _findDetectorVariables(self, det):
		lst = []
		for k in self.plot.var.keys():
			try:
				ks = k.split(".")
				if int(ks[1])==det:
					lst.append(ks[0])
			except:
				pass
		return lst

	# ----------------------------------------------------------------------
	def _moveDetectorVariables(self, src, dst):
		lst = self._findDetectorVariables(src)
		# Make a copy of all variables, and delete the original
		for k in lst:
			vsrc = "%s.%d"%(k,src)
			self.plot["%s.%d"%(k,dst)] = self.plot[vsrc]
			del self.plot.var[vsrc]

	# ----------------------------------------------------------------------
	# Delete list
	# ----------------------------------------------------------------------
	def delete(self):
		# protect against accidental delete
		if self.focus_get() is not self.listbox: return

		self.get(self.plot)
		self.current = None
		act = self.listbox.index(ACTIVE)
		sel = map(int, self.listbox.curselection())
		sel.reverse()
		ndets = int(self.plot.get("ndetectors",0))

		for det in sel:
			for k in self._findDetectorVariables(det):
				del self.plot.var["%s.%d"%(k,det)]
			if det < self.listbox.size()-1:
				for d in range(det+1,ndets):
					self._moveDetectorVariables(d, d-1)
			ndets -= 1
			self.listbox.delete(det)

		self.plot["ndetectors"] = ndets
		self.selectDetector()

	# ----------------------------------------------------------------------
	# Copy detector in clipboard
	# ----------------------------------------------------------------------
	def copy(self, event=None):
		if self.focus_get() is not self.listbox: return
		self.get(self.plot)

		clip = {}
		ndets = 0
		for det in map(int, self.listbox.curselection()):
			for k in self._findDetectorVariables(det):
				clip["%s.%d"%(k,ndets)] = self.plot["%s.%d"%(k,det)]
			ndets += 1
		clip["@n"] = ndets

		# Write to clipboard
		sio = StringIO()
		sio.write(_DETECTOR_CLIP)
		pickler = pickle.Pickler(sio)
		pickler.dump(clip)
		self.clipboard_clear()
		self.clipboard_append(sio.getvalue())
		return "break"

	# ----------------------------------------------------------------------
	def cut(self, event=None):
		self.copy()
		self.delete()

	# ----------------------------------------------------------------------
	def paste(self, event=None):
		if self.focus_get() is not self.listbox: return
		self.get(self.plot)
		self.current = None

		try: clipboard = self.selection_get(selection='CLIPBOARD')
		except: return
		if clipboard.startswith(_DETECTOR_CLIP):
			# Pickler format
			unpickler = pickle.Unpickler(StringIO(clipboard[len(_DETECTOR_CLIP):]))
			try:
				clip = unpickler.load()
			except:
				return
		else:
			self.flair.notify("Cannot paste","Invalid content in clipboard",tkFlair.NOTIFY_WARNING)
			return

		self.listbox.selection_clear(0,END)

		ndets = int(self.plot.get("ndetectors",0))

		for i in range(clip.get("@n",0)):
			for k,v in clip.items():
				try:
					ks = k.split(".")
					if int(ks[1]) == i:
						self.plot["%s.%d"%(ks[0], ndets)] = v
				except:
					pass
			self.listbox.insert(END, self.plot["name.%d"%(ndets)])
			self.listbox.selection_set(END)
			self.listbox.activate(END)
			self.listbox.see(END)
			ndets += 1

		self.plot["ndetectors"] = ndets
		self.selectDetector()

	# ----------------------------------------------------------------------
	# Clone detector
	# ----------------------------------------------------------------------
	def clone(self):
		self.get(self.plot)
		self.current = None

		sel = map(int, self.listbox.curselection())
		self.listbox.selection_clear(0,END)

		ndets = int(self.plot.get("ndetectors",0))

		for det in sel:
			for k in self._findDetectorVariables(det):
				self.plot["%s.%d"%(k,ndets)] = self.plot["%s.%d"%(k,det)]
			self.listbox.insert(END, self.plot["name.%d"%(ndets)])
			self.listbox.selection_set(END)
			self.listbox.activate(END)
			self.listbox.see(END)
			ndets += 1

		self.plot["ndetectors"] = ndets
		self.selectDetector()

	# ----------------------------------------------------------------------
	# Swap two items in the list
	# ----------------------------------------------------------------------
	def swap(self, a, b):
		self.listbox.swap(a,b)

		self._moveDetectorVariables(a, 999)
		self._moveDetectorVariables(b,   a)
		self._moveDetectorVariables(999, b)

	# ----------------------------------------------------------------------
	# Move up select items by one
	# ----------------------------------------------------------------------
	def moveUp(self, event=None):
		self.get(self.plot)
		self.current = None
		for i in map(int,self.listbox.curselection()):
			if i==0: continue
			prev = i-1
			if not self.listbox.selection_includes(prev):
				act = self.listbox.index(ACTIVE)
				self.swap(prev,i)
				self.listbox.selection_set(prev)
				self.listbox.see(prev)
				if act == i: self.listbox.activate(prev)
			self.selectDetector()
		return "break"

	# ----------------------------------------------------------------------
	# Move down select items by one
	# ----------------------------------------------------------------------
	def moveDown(self, event=None):
		self.get(self.plot)
		self.current = None
		sz  = self.listbox.size()-1
		lst = map(int,self.listbox.curselection())
		lst.reverse()
		for i in lst:
			if i >= sz: continue
			next = i+1
			if not self.listbox.selection_includes(next):
				act = self.listbox.index(ACTIVE)
				self.swap(i,next)
				self.listbox.selection_set(next)
				self.listbox.see(next)
				if act == i: self.listbox.activate(next)
			self.selectDetector()
		return "break"
