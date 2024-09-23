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
# Date:	06-Sep-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import re
import time
import string
import fnmatch
from stat import *
from log import say

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import pipe
import Input
import Utils
import Manual
import Ribbon
import Project
import RunList
import tkFlair
import tkExtra
import FlairRibbon
import bFileDialog

_INPUT_CYCLE = "input"
_DATA_CYCLE  = "data"
_PLOT_CYCLE  = "plot"
_COMP_CYCLE  = "compile"
_TEMP_CYCLE  = "temporary"
SPECIAL_CYCLES = [_INPUT_CYCLE, _DATA_CYCLE, _PLOT_CYCLE, _COMP_CYCLE, _TEMP_CYCLE]

#===============================================================================
# Files Tab in Page Frame
#===============================================================================
class FilesTab(FlairRibbon.FlairTab):
	#----------------------------------------------------------------------
	def __init__(self, master, page, **kw):
		FlairRibbon.FlairTab.__init__(self, master, page, **kw)

		# --- Run-Cycle splitter ---
		splitter = PanedWindow(self, orient=HORIZONTAL,
			showhandle=0,
			sashrelief=RAISED, sashwidth=4, opaqueresize=1)
		splitter.pack(expand=YES, fill=BOTH)

		# --- Cycles ---
		frame = LabelFrame(splitter, text="Cycles", foreground="DarkBlue")
		splitter.add(frame, minsize=40)

		self.cycleList = tkExtra.ExListbox(frame, borderwidth=0,
			selectmode=EXTENDED,
			selectborderwidth=0,
			takefocus=True,
			width=15,
			height=7,
			exportselection=FALSE)
		self.cycleList.pack(side=LEFT,expand=YES,fill=BOTH)

		sb = Scrollbar(frame, orient=VERTICAL, takefocus=False,
			command=self.cycleList.yview)
		sb.pack(side=RIGHT,fill=Y)
		self.cycleList.config(yscrollcommand=sb.set)

		# Bind events
		self.cycleList.bind('<<ListboxSelect>>', self.refresh)

		# --- Files ---
		frame = Frame(splitter)
		splitter.add(frame)

		self.fileList = tkExtra.ColorMultiListbox(frame,
				(('File', 30, None),
				 ('Type', 10, None),
				 ('Size',  8, None),
				 ('Date', 15, None)),
				 background="White")
		self.fileList.pack(expand=YES, fill=BOTH)

		# Bind events
		self.fileList.bindList('<Double-1>',     self.viewer)
		self.fileList.bindList('<Return>',       self.viewer)
		self.fileList.bindList('<KP_Enter>',     self.viewer)
		self.fileList.bindList('<Delete>',       self.delete)

		self.fileList.setPopupMenu(
				[('View',     0, self.viewer,  tkFlair.icons["view"]),
				 ('Edit',     0, self.editor,  tkFlair.icons["edit"]),
				 ('To ascii', 0, self.toAscii, tkFlair.icons["ascii"]),
				 ('Refresh',  0, self.refresh, tkFlair.icons["refresh"]),
				 None,
				 ('Delete',   0, self.delete,  tkFlair.icons["x"])] )
		self.ascending = True

		self.saveCyclesSelected = None
		self.pattern = None
		self._files = {}	# to avoid duplicates

	# ----------------------------------------------------------------------
	# create the ribbon buttons for the Run tab
	# ----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairTab.createRibbon(self)

		group = Ribbon.LabelGroup(self.ribbon, "Files")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.label["background"] = Ribbon._BACKGROUND_GROUP3

		group.grid2rows()

		# ---
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["ascii"],
				text=tkExtra.ARROW_RIGHT+"Ascii",
				compound=LEFT,
				anchor=W,
				command=self.toAscii,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Convert file to ASCII")

		# ----
		col,row = 0,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Delete",
				compound=LEFT,
				command=self.delete,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Delete files [Del]")

		# ---
		col,row = 0,2
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["refresh"],
				text="Refresh",
				compound=LEFT,
				anchor=W,
				command=self.refresh,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Refresh files")

		# ---
		col,row = 1,1
		self.patternEntry = tkExtra.LabelEntry(group.frame,
					"Filter", "DarkGray",
					background="White", #Ribbon._BACKGROUND,
					width=20)
		self.patternEntry.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(self.patternEntry,
				"Filter files with a pattern (wildchars: * and ?)")
		self.patternEntry.bind("<Return>",   self.refresh)
		self.patternEntry.bind("<KP_Enter>", self.refresh)

		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["filter"],
				command=self.refresh,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col+1, padx=0, pady=0, sticky=E)
		tkExtra.Balloon.set(b, "Filter files with pattern")

#		b = Ribbon.LabelButton(group.frame,
#				text="X",
#				command=self.patternClear,
#				background=Ribbon._BACKGROUND)
#		b.grid(row=row, column=col+2, padx=0, pady=0, sticky=E)
#		tkExtra.Balloon.set(b, "Delete filter")

		# ======================
		group = Ribbon.LabelGroup(self.ribbon, "Action")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.label["background"] = Ribbon._BACKGROUND_GROUP3

		group.grid3rows()

		# ----
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["viewer32"],
				text="Viewer",
				compound=TOP,
				command=self.viewer,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Open in viewer")

		# ----
		col,row = 1,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["edit32"],
				text="Editor",
				compound=TOP,
				command=self.editor,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Open in editor")

		return self.ribbon

	# ----------------------------------------------------------------------
	def activate(self):
		self.fillCycles()
		self.refreshStatus()

	# ----------------------------------------------------------------------
	# Lose focus
	# ----------------------------------------------------------------------
	def _focusOut(self, event=None):
		if self.project is None: return
		self.page.runList.updateSelection()
		self.saveCyclesSelected = self.cycleList.getSelected()

	# ----------------------------------------------------------------------
	# Show Frame
	# ----------------------------------------------------------------------
	def showFrame(self):
		tkFlair.TitleFrame.showFrame(self)

		self.page.runList.show(self.flair)
		self.page.runList.update_idletasks()
		self.page.runList.focus_set()
		self.fillCycles(select=False)
		self.cycleList.selectSaved(self.saveCyclesSelected, 0)
		self.refresh()

	# ----------------------------------------------------------------------
	def patternClear(self):
		self.patternEntry.delete(0,END)
		self.pattern = ""
		self.refresh()

	# ----------------------------------------------------------------------
	def selectAll(self):
		self.page.runList.selection_set(0,END)
		self.fillCycles()

	# ----------------------------------------------------------------------
	def selectClear(self):
		self.page.runList.selection_clear(0,END)
		self.fillCycles()

	# ----------------------------------------------------------------------
	# Select Run
	# ----------------------------------------------------------------------
	def runSelect(self, event):
		i = self.page.runList.nearest(event.y)
		self.page.runList.focus_set()
		self.page.runList.activate(i)
		self.fillCycles()

	# ----------------------------------------------------------------------
	# Fill Run Cycles
	# ----------------------------------------------------------------------
	def fillCycles(self, event=None, select=True):
		#runs = self.page.runList.selection()
		cycles = SPECIAL_CYCLES[:]	#copy
		for run in self.page.runList.runSelection():
			if run is None: continue
			inpname = run.getInputBaseName()
			inpfile = "%s.inp"%(inpname)
			inpname = re.escape(inpname)
			rundir  = run.getDir()

			cpat1  = re.compile(r"^%s(\d\d\d)[._].+$" % (inpname))
			cpat2  = re.compile(r"^ran%s(\d\d\d)$"    % (inpname))
			patout = re.compile(r"^%s(\d\d\d).out$"   % (inpname))
			patlog = re.compile(r"^%s(\d\d\d).log$"   % (inpname))

			try:
				os.stat(rundir)
			except OSError:
				continue

			for f in Utils.listdir(rundir):
				m = cpat1.match(f)
				if m is None: m = cpat2.match(f)
				if m:
					c = m.group(1)
					try: i = cycles.index(c)
					except ValueError:
						cycles.append(c)

				# Check sub directory if ...out file exist
				if not Project._PAT_DIR.match(f): continue

				fn = os.path.join(rundir,f)
				try: s = os.lstat(fn)
				except OSError: continue

				if f.startswith(Project._FLUKA_CYCLE) and S_ISDIR(s[ST_MODE]):
					for fsub in Utils.listdir(fn):
						if fsub==inpfile or \
						   patout.match(fsub) or \
						   patlog.match(fsub):
							cycles.append(f)
							break
		cycles.sort()

		# Populate list
		self.cycleList.delete(0,END)
		for i in cycles:
			self.cycleList.insert(END,i)

		self.cycleList.selection_set(0)
		if select: self.refresh()

	# ----------------------------------------------------------------------
	# Select Cycle
	# ----------------------------------------------------------------------
	def cycleSelect(self, event):
		i = self.page.runList.nearest(event.y)
		self.cycleList.focus_set()
		self.cycleList.activate(i)
		self.refresh()

	# ----------------------------------------------------------------------
	def _compileFiles(self):
		totalSize = 0

		# Source files
		d = self.project.dir
		for f in self.project.sourceList:
			fn,ext = os.path.splitext(f)
			if ext in (".f", ".for", ".fpp", ".c", ".cpp", ".cc"):
				fo = "%s.o"%(fn)
				try:
					os.stat(fo)
					totalSize += self.addFile(d, fo)
				except:
					pass

		# ... executable
		try:
			os.stat(self.project.exe)
			totalSize += self.addFile(d, self.project.exe)
		except:
			pass

		# ... and map
		fn,ext = os.path.splitext(self.project.exe)
		fm = "%s.map"%(fn)
		try:
			os.stat(fm)
			totalSize += self.addFile(d, fm)
		except:
			pass

		return totalSize

	# ----------------------------------------------------------------------
	def _plotFiles(self):
		totalSize = 0
		d = self.project.dir
		for f in Utils.listdir(d):
			fn,ext = os.path.splitext(f)
			if fn=="gplevh" or fn=="doslev":
				totalSize += self.addFile(d, f)
			else:
				for plot in self.project.plots:
					if fn==plot.name:
						totalSize += self.addFile(d, f)
		return totalSize

	# ----------------------------------------------------------------------
	def _tmpFiles(self):
		totalSize = 0
		tmp = Project.tmpPrefix
		rantmp = "ran%s"%(tmp)
		d = self.project.dir
		for f in Utils.listdir(d):
			fn,ext = os.path.splitext(f)
			if fn.startswith(tmp) or fn.startswith(rantmp):
				totalSize += self.addFile(d, f)
		return totalSize

	# ----------------------------------------------------------------------
	# Refresh file list
	# ----------------------------------------------------------------------
	def refresh(self, event=None):
		self.page.busy()
		savesort = self.fileList.saveSort()
		#runs   = self.page.runList.curselection()
		cycles = [self.cycleList.get(int(x))
				for x in self.cycleList.curselection()]
		self.pattern = self.patternEntry.get()

		self._files.clear()
		self.fileList.delete(0,END)
		self.totalSize = 0
		tdata = 0
		tcomp = 0
		tplot = 0
		ttemp = 0
		first = True

		tstart = time.time()
#		self.listdir.clear()
		for run in self.page.runList.runSelection():
			#run     = self.page.runList.getRun(r)
			inpname = re.escape(run.getInputBaseName())
			rundir  = run.getDir()
			mainout = False

			for c in cycles:
				if c == _INPUT_CYCLE:
					t = 0
					inpfile = "%s/%s.inp"%(rundir, run.getInputBaseName())
					try:
						os.stat(inpfile)
						t += self.addFile(rundir, inpfile)
					except OSError:
						pass
					tdata     += t
					self.totalSize += t

				elif c == _DATA_CYCLE:
					# Populate with output files
					t = 0
					try:
						for f in Utils.listdir(rundir):
							for usr in run.usrinfo:
								name = os.path.splitext(usr.name())[0]
								if f.startswith(os.path.basename(name)):
									t += self.addFile(rundir, f)
					except OSError:
						pass

					tdata     += t
					self.totalSize += t

				elif c == _COMP_CYCLE and first:
					# Show compile only for the first cycle
					t = self._compileFiles()
					tcomp += t
					self.totalSize += t

				elif c == _PLOT_CYCLE and first:
					t = self._plotFiles()
					tplot += t
					self.totalSize += t

				elif c == _TEMP_CYCLE and first:
					t = self._tmpFiles()
					ttemp += t
					self.totalSize += t

				elif c[0:6] == Project._FLUKA_CYCLE and first:
					cdir = os.path.join(rundir,c)
					try:
						s = os.stat(cdir)
						self.addFile(rundir, c)
						for f in os.listdir(cdir):
							self.totalSize += self.addFile(cdir, f,
								not tkFlair._showFlukaFiles)
					except OSError:
						pass
				else:
					# Populate with cycle files
					cpat1 = re.compile(r"^%s%s[._].+$" % (inpname,c))
					cpat2 = re.compile(r"^ran%s%s$" % (inpname,c))
					try:
						for f in Utils.listdir(rundir):
							if cpat1.match(f) or cpat2.match(f):
								self.totalSize += self.addFile(rundir, f)
					except OSError:
						self.flair.notify("OSError",sys.exc_info()[1],
								tkFlair.NOTIFY_WARNING)

					if not mainout:
						mainout = True
						outf = "%s.out"%(run.getInputBaseName())
						try:
							os.stat(outf)
							self.totalSize += self.addFile(rundir, outf)
						except OSError:
							pass

			first = False

		# Remove non existing cycles
		#if tdata == 0: self.cycleList.deleteByName(_DATA_CYCLE)
		#if tcomp == 0: self.cycleList.deleteByName(_COMP_CYCLE)
		#if tplot == 0: self.cycleList.deleteByName(_PLOT_CYCLE)
		#if ttemp == 0: self.cycleList.deleteByName(_TEMP_CYCLE)

		# Display total
		if Input._developer: say("FileFrame::refresh time=",time.time() - tstart)
		self.fileList.restoreSort(savesort)
		self.page.notBusy()

	# ----------------------------------------------------------------------
	def refreshStatus(self):
		self.page.setStatus("Files: %d Total Size: %d" \
			% (self.fileList.size(), self.totalSize))

	# ----------------------------------------------------------------------
	# Add a file to the list
	# ----------------------------------------------------------------------
	def addFile(self, path, filename, checklink=False):
		if self.pattern:
			if not fnmatch.fnmatch(filename, self.pattern):
				return 0
		f  = os.path.join(path, filename)
		if f in self._files: return 0
		self._files[f] = True
		fn = self.project.relativePath(f)
		try:
			s = os.lstat(f)
			if checklink and S_ISLNK(s[ST_MODE]): return 0
			ext,color = bFileDialog.fileTypeColor(f)
			self.fileList.insert(END, (fn, ext, s[ST_SIZE],
				 time.strftime(bFileDialog._TIME_FORMAT,
					time.localtime(s[ST_MTIME]))))
			if color is not None:
				self.fileList.itemconfig(END,foreground=color)
			return s[ST_SIZE]
		except OSError:
			self.fileList.insert(END, (fn, "", 0, "-"))
			return 0

	# ----------------------------------------------------------------------
	# Delete selected files
	# ----------------------------------------------------------------------
	def delete(self, event=None):
		log = self.flair.newLog("Files", "Files")

		files = [ self.fileList.get(int(x))[0]
				for x in self.fileList.curselection() ]
		if len(files)==0:
			self.flair.notify("No files selected",
				"Please select the files you want to delete",
				tkFlair.NOTIFY_WARNING)
			return

		elif len(files)>20:
			msg = "Are you sure you want to delete %d files?" \
				% (len(files))
		else:
			msg = "Are you sure you want to delete the files:\n%s" \
				% (" ".join(files))

		if not tkFlair.askyesno("Delete Files", msg, parent=self.ribbon): return

		# First delete files
		refreshCycles = False

		for f in files:
			try: s = os.lstat(f)
			except OSError: continue

			# Skip input file
			if f == self.project.inputFile: continue

			if S_ISDIR(s[ST_MODE]): continue
			if S_ISLNK(s[ST_MODE]):
				say("Removing link: %s"%(f))
			else:
				say("Deleting file: %s"%(f))
			try:
				os.remove(f)
			except OSError:
				say("Error deleting file: %s"%(f))

		# Try to delete subdirectories if empty
		for f in files:
			try: s = os.stat(f)
			except OSError: continue
			if not S_ISDIR(s[ST_MODE]): continue

			# Try to delete all links in directory
			if not tkFlair._showFlukaFiles:
				for fn in os.listdir(f):
					ff = os.path.join(f,fn)
					s = os.lstat(ff)
					if S_ISLNK(s[ST_MODE]):
						say("Removing link: %s"%(ff))
						try:
							os.remove(ff)
						except OSError:
							say("Error removing link: %s"%(ff))

			say("Removing directory: %s"%(f))
			try:
				os.rmdir(f)
				refreshCycles = True
			except OSError:
				say("Error removing directory: %s"%(f))

		self.refresh()
		if refreshCycles: self.fillCycles()

	# ----------------------------------------------------------------------
	# Edit files
	# ----------------------------------------------------------------------
	def editor(self):
		tkFlair.editor([ self.fileList.get(int(x))[0]
				for x in self.fileList.curselection() ])

	# ----------------------------------------------------------------------
	# Edit files
	# ----------------------------------------------------------------------
	def viewer(self, event=None):
		self.flair.viewer([ self.fileList.get(int(x))[0]
				for x in self.fileList.curselection() ])

	# ----------------------------------------------------------------------
	def toAscii(self, event=None):
		cmd = Project.command("usbrea")

		files = []
		log = self.flair.newLog("Files", "Files")
		for x in self.fileList.curselection():
			f = self.fileList.get(int(x))[0]
			files.append(f)

			say("Processing: %s cmd=%s"%(f, cmd))
			self.page.setStatus("Converting %s to ascii"%(f))

			inp = [f, "%s.lis"%(f)]
			try:
				rc, out = pipe.system(cmd, inp)
			except (KeyboardInterrupt, IOError):
				say("*** Interrupt ***")
				return True
			for line in out:
				say(line)
		self.refresh()
		self.flair.notify("To Ascii",
			"Files converted to ascii: %s\n"%(" ".join(files)))

	# ----------------------------------------------------------------------
	def help(self):
		Manual.show(":Run:Files:")
