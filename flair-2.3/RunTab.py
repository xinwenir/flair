#!bin/env python
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

import os, re
import time
import signal
import string
import subprocess

from stat  import *
from log   import say
from bmath import ZeroIntDict

try:
	from tkinter import *
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	import tkinter.messagebox as messagebox

import Input
import Utils
import Manual
import Ribbon
import Project
import RunList
import tkFlair
import FlairRibbon
import bFileDialog

import Unicode
import tkExtra
import tkDialogs

STATUS_STR = [
	"Not Running",
	"Waiting to attach",
	"Running",
	"Finished OK",
	"Finished with ERRORS",
	"*** TIMED-OUT ***" ]

_ON    = Unicode.BALLOT_BOX_WITH_X
_OFF   = Unicode.BALLOT_BOX
_MULTI = Unicode.DIE_FACE_6

#===============================================================================
# Run Tab in Page Frame
#===============================================================================
class RunTab(FlairRibbon.FlairTab):
	#----------------------------------------------------------------------
	def __init__(self, master, page, **kw):
		FlairRibbon.FlairTab.__init__(self, master, page, **kw)

		# -----
		row = 0
		Label(self, text="Title").grid(row=row, column=0, sticky=E)
		self.title = Entry(self, background="White", justify=LEFT)
		tkExtra.Balloon.set(self.title, "Input Title")
		self.title.grid(row=row, column=1, columnspan=5, sticky=EW)
		self.title.bind('<Button-3>', self.toggleField)

		# -----
		row += 1
		Label(self, text="Primaries").grid(row=row, column=0, sticky=NE)
		self.primaries = tkExtra.FloatEntry(self, background="White", width=5)
		tkExtra.Balloon.set(self.primaries,
				"Override number of starting particles")
		self.primaries.grid(row=row, column=1, sticky=EW)
		self.primaries.bind('<Button-3>', self.toggleField)

		# -----
		Label(self, text="Rnd").grid(row=row, column=2, sticky=NE)
		self.rnd = tkExtra.IntegerEntry(self, background="White", width=5)
		tkExtra.Balloon.set(self.rnd, "Override random number seed")
		self.rnd.grid(row=row, column=3, sticky=EW)
		self.rnd.bind('<Button-3>', self.toggleField)

		self.rndButton = Button(self, image=tkFlair.icons["RANDOMIZ"],
				command=self.getRandomNumberAsSeed)
		tkExtra.Balloon.set(self.rndButton, "Radnom number")
		self.rndButton.grid(row=row, column=4, sticky=NSEW)
		self.rndButton.bind('<Button-3>', self.toggleField)

		# -----
		row += 1
		Label(self, text="Time").grid(
				row=row, column=0, sticky=E)
		self.time = tkExtra.FloatEntry(self, background="White", width=5)
		tkExtra.Balloon.set(self.time,
				"Time limit of the run")
		self.time.grid(row=row, column=1, sticky=EW)
		self.time.bind('<Button-3>', self.toggleField)

		# -----
		Label(self, text="Exe").grid(row=row, column=2, sticky=E)
		self.exe = Label(self,
				foreground="DarkBlue", background="White",
				anchor=W, relief=SUNKEN, padx=1, width=10)
		tkExtra.Balloon.set(self.exe, "Override Executable file")
		self.exe.grid(row=row, column=3, sticky=EW)
		self.exe.bind('<Button-1>', self.loadExecutable)
		self.exe.bind('<Button-3>', self.toggleField)

		self.bed = Button(self, image=tkFlair.icons["x"],
				command=self.defaultExecutable)
		tkExtra.Balloon.set(self.bed, "Default Executable")
		self.bed.grid(row=row, column=4, sticky=NSEW)
		self.bed.bind('<Button-3>', self.toggleField)
		self.bel = Button(self, image=tkFlair.icons["load"],
				command=self.loadExecutable)
		tkExtra.Balloon.set(self.bel, "Load Executable")
		self.bel.grid(row=row, column=5, sticky=NSEW)
		self.bel.bind('<Button-3>', self.toggleField)

		# -----
		row += 1
		Label(self, text="Defines").grid(row=row, column=0, sticky=W)
		self.bd = Button(self, text="Default Defines", pady=0, padx=2,
				takefocus=False, command=self.selectDefaultDefines)
		tkExtra.Balloon.set(self.bd, "Reset to default defines")
		self.bd.grid(row=row, column=1, sticky=EW)
		self.bd.bind('<Button-3>', self.toggleField)

		# -----
		row += 1
		self.defineList = tkExtra.MultiListbox(self,
				((' ',     2, None),
				 ('Name', 16, None),
				 ('Value',32, None)),
				 background="White",
				 stretch="last",
				 height=5)
		self.defineList.bindList("<space>",		self.toggleDefine)
		self.defineList.bindList("<Double-1>",		self.toggleDefine)
		self.defineList.bindList("<Return>",		self.editDefine)
		self.defineList.bindList("<KP_Enter>",		self.editDefine)
		self.defineList.lists[0].bind("<Button-1>",	self.defineListClick)
		self.defineList.bindList("<ButtonRelease-1>",	self.defineListRelease)
		self.defineList.bindLeftRight()
		#self.defineList.setPopupMenu([('Default', 0,	self.selectDefaultDefines)])
		# unbind popupmenu
		for l in self.defineList.lists: l.unbind("<Button-3>")
		self.defineList.bindList('<Button-3>',		self.toggleField)
		self.defineList.grid(row=row,  column=0, columnspan=6, sticky=NSEW)

		self.grid_columnconfigure(1, weight=1)
		self.grid_columnconfigure(3, weight=1)
		self.grid_rowconfigure(4, weight=1)

		# --- Progress ---
		frame = LabelFrame(self, text="Progress",
				foreground="DarkBlue")
		frame.grid(row=6, column=0, columnspan=6, sticky=EW)

		row = 0
		col = 0
		Label(frame, text="Status:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.pstatus = Label(frame, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.pstatus.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.pstatus, "Run status")

		col += 1
		Label(frame, text="Input:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.pinput = Label(frame, text="", justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.pinput.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.pinput, "Input file name")

		col += 1
		Label(frame, text="Dir:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.pdir = Label(frame, text="", justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.pdir.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.pdir, "Temporary directory where run takes place")

		row += 1
		col = 0
		Label(frame, text="Started:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.pstarted = Label(frame, text="", justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.pstarted.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.pstarted, "Time when run started")

		col += 1
		Label(frame, text="ETA:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.peta = Label(frame, text="", justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.peta.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.peta, "Estimated time of accomplishment assuming no other process are running on CPU")

		col += 1
		Label(frame, text="Time/prim:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.ptime = Label(frame, text="", justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.ptime.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.ptime, "CPU Time per primary")

		row += 1
		col = 0
		Label(frame, text="Elapsed:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.pelapsed = Label(frame, text="", justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.pelapsed.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.pelapsed, "Time elapsed since the start of the current cycle")

		col += 1
		Label(frame, text="Cycle:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.premain = Label(frame, text="", justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.premain.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.premain, "CPU Time remaining for present cycle")

		col += 1
		Label(frame, text="Run:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.prunremain = Label(frame, text="", justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.prunremain.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.prunremain, "Total CPU Time remaining for run")

		frame.grid_columnconfigure(1, weight=1)
		frame.grid_columnconfigure(3, weight=1)
		frame.grid_columnconfigure(5, weight=1)

		# Progress bar
		row += 1
		col = 0
		Label(frame, text="Cycles:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.cyclePBar = tkExtra.ProgressBar(frame, background="DarkGray", height=20)
		self.cyclePBar.grid(row=row, column=col, columnspan=5, sticky=NSEW)
		self.cyclePBar.setShowTime(False)

		row += 1
		col = 0
		Label(frame, text="Primaries:", pady=2).grid(row=row, column=col,
				sticky=E)
		col += 1
		self.primaryPBar = tkExtra.ProgressBar(frame, background="DarkGray", height=20)
		self.primaryPBar.grid(row=row, column=col, columnspan=5, sticky=NSEW)
		self.primaryPBar.setShowTime(False)

		# Set initial values
		self.pcycle_perc  = (0.0, 0.0)
		self.pprim_perc   = (0.0, 0.0)
		self.timer	  = None
		self._updating	  = False	# Semaphore for updating the values
		self._varSetting  = True	# Semaphore for setting the cycles

	# ----------------------------------------------------------------------
	# create the ribbon buttons for the Run tab
	# ----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairTab.createRibbon(self)

		# ========== Input ===========
		group = Ribbon.LabelGroup(self.ribbon, "Input")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["add32"],
				text="Add",
				compound=TOP,
				command=self.page.add,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, rowspan=3, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add run [Ins]")

		# ---
		col,row = 1,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["up"],
				text="Move Up",
				compound=LEFT,
				command=self.page.moveUp,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Move selected run up [Ctrl-Up]")

		# ---
		col,row = 1,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["down"],
				text="Move Down",
				compound=LEFT,
				command=self.page.moveDown,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Move selected runs down [Ctrl-Down]")

		# ---
		col,row = 1,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["rename"],
				text="Rename",
				compound=LEFT,
				command=self.page.rename,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Rename run [F2]")

		# ---
		col,row = 2,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Remove",
				compound=LEFT,
				anchor=W,
				command=self.page.delete,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Remove run [Del]")

		# ---
		col,row = 2,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["loop"],
				text="Loop",
				compound=LEFT,
				anchor=W,
				command=self.page.loop,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Loop runs over a variable")

		# ---
		col,row = 2,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clone"],
				text="Clone",
				compound=LEFT,
				anchor=W,
				command=self.page.clone,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Clone run [Ctrl-D]")

		# ========== Job ===========
		group = Ribbon.LabelGroup(self.ribbon, "Job")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# Cycles
		self.prev   = IntVar()
		self.cycles = IntVar()
		self.last   = IntVar()

		self.prev.trace(  "w", self.prevChange)
		self.cycles.trace("w", self.prevChange)	# Use the same function
		self.last.trace(  "w", self.lastChange)

		# ---
		col,row = 0,0
		b = Label(group.frame,
				image=tkFlair.icons["queue"],
				font=Ribbon._FONT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=E)

		self.queueList = Ribbon.LabelCombobox(group.frame,
					width=10)
		self.queueList.grid(row=row, column=col+1, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(self.queueList, "Submitting queue")

		# ----
		col,row = 0,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["continue"],
				text="Continue",
				compound=LEFT,
				command=self.continueFromLast,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Continue from last cycle")

		# ----
		col,row = 0,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["attach"],
				text="Attach",
				compound=LEFT,
				command=self.attach2Run,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, columnspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Attach to a running process")

		# --- Cycles
		col,row = 2,0

		b = Label(group.frame, text="Prev:",
				anchor=E,
				font=Ribbon._FONT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		self.prevbox = Spinbox(group.frame,
					text=self.prev,
					from_=0, to=1000,
					font=Ribbon._FONT,
					background="White", #Ribbon._BACKGROUND,
					width=4)
		self.prevbox.grid(row=row, column=col+1, padx=0, pady=0, sticky=NSEW)
		self.prevbox.bind('<Button-3>', self.toggleField)
		tkExtra.Balloon.set(self.prevbox, "Previous run last cycle number")

		# ---
		col,row = 2,1
		b = Label(group.frame, text="No:",
				anchor=E,
				font=Ribbon._FONT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		self.cyclesbox = Spinbox(group.frame,
					text=self.cycles,
					from_=0, to=1000,
					font=Ribbon._FONT,
					background="White", #Ribbon._BACKGROUND,
					width=4)
		self.cyclesbox.grid(row=row, column=col+1, padx=0, pady=0, sticky=NSEW)
		self.cyclesbox.bind('<Button-3>', self.toggleField)
		tkExtra.Balloon.set(self.cyclesbox, "Number of cycles to run")

		# ---
		col,row = 2,2
		b = Label(group.frame, text="To:",
				anchor=E,
				font=Ribbon._FONT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		self.lastbox = Spinbox(group.frame,
					text=self.last,
					from_=0, to=1000,
					font=Ribbon._FONT,
					background="White", #Ribbon._BACKGROUND,
					width=4)
		self.lastbox.grid(row=row, column=col+1, padx=0, pady=0, sticky=NSEW)
		self.lastbox.bind('<Button-3>', self.toggleField)
		tkExtra.Balloon.set(self.lastbox, "Last cycle number")

		# ==================
		group = Ribbon.LabelGroup(self.ribbon, "Action")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clean32"],
				text="Clean",
				compound=TOP,
				command=self.clean,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Remove all generated files")

		# ---
		col,row = 1,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["STOP"],
				text="Cycle",
				compound=LEFT,
				anchor=W,
				command=self.stopCycle,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Stop current cycle by creating a fluka.stop")

		# ---
		col,row = 1,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["STOP"],
				text="Run",
				compound=LEFT,
				anchor=W,
				command=self.stopRun,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Stop run by creating a rfluka.stop")

		# ---
		col,row = 1,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Kill",
				compound=LEFT,
				anchor=W,
				command=self.killRun,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Kill run process if possible")

		# ----
		col,row = 2,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["refresh32"],
				text="Refresh",
				compound=TOP,
				command=self.refresh,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Refresh progress information")


		# ---
		col,row = 3,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["play32"],
				text="Start",
				compound=TOP,
				command=self.execute,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, rowspan=3, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Execute run [Ctrl-Enter]")

		return self.ribbon

	# ----------------------------------------------------------------------
	# select default defines
	# ----------------------------------------------------------------------
	def selectDefaultDefines(self, runs=None):
		self.defineList.delete(0,END)

		# Update defines
		lst = []
		selected = []
		for card in self.project.input.cardsSorted("#define",2):
			if card._indent != 0: continue
			sdum = card.sdum()
			for n,w in lst:
				if n==sdum: break
			else:
				lst.append((sdum, ""))
				if runs is None and card.enable:
					selected.append((sdum, "",False))

		# Highlight
		if runs:
			selected = []
			for n,v in lst:
				value = None
				count = 0
				for run in runs:
					for rn,rv in run.defines:
						if n!=rn: continue	# Skip

						if value is None:
							value = rv	# First time set the value
							count = 1
						elif value == rv:
							count += 1	# Count same occurrences
							break
				if value is not None:
					selected.append((n,value,count!=len(runs)))

		# Show defines
		for n,v in lst:
			# Search in selected ones
			for sn,sv,m in selected:
				if n==sn:
					if m:
						self.defineList.insert(END,(_MULTI,sn,sv))
					else:
						self.defineList.insert(END,(_ON,sn,sv))
					break
			else:
				self.defineList.insert(END,(_OFF,n,v))

	# ----------------------------------------------------------------------
	def defineListClick(self, event):
		if len(event.widget.curselection())>1:
			active = event.widget.nearest(event.y)
			if event.widget.selection_includes(active):
				self.toggleDefine()
				return "break"

	# ----------------------------------------------------------------------
	def defineListRelease(self, event):
		"""Button Release on the define listbox"""
		if len(self.defineList.curselection())>1: return
		active = self.defineList.nearest(event.y)
		if active is not None:
			self.defineList.activate(active)
			try:
				idx = self.defineList.lists.index(event.widget)
			except ValueError:
				return
			if idx == 0:
				self.toggleDefine()
			elif idx == 2:
				self.editDefine()

	# ----------------------------------------------------------------------
	def toggleDefine(self, event=None):
		# toggle state
		lst = self.defineList.lists[0]
		sel = lst.curselection()
		if not sel: return

		status = lst.get(sel[0])
		if status == _ON:		# [X] -> [ ]
			status = _OFF
		elif len(self.page.runList.curselection())>1:
			if status == _OFF:	# [ ] -> [-]
				status = _MULTI
			else:
				status = _ON	# [-] -> [X]
		else:
			status = _ON		# [ ] -> [X]

		for i in sel:
			lst.set(i,status)
			self.defineList.selection_set(i)
			self.defineList._updateSelect()
		self.get()

	# ----------------------------------------------------------------------
	def editDefine(self, event=None):
		sel = list(map(int,self.defineList.curselection()))
		idx = self.defineList.index(ACTIVE)
		self.defineList.activate(idx)
		oldValue = self.defineList.lists[2].get(idx)
		edit = tkExtra.InPlaceEdit(self.defineList.lists[2])
		if edit.value is not None and edit.value != oldValue:
			# Activate define
			for i in sel:
				self.defineList.lists[0].set(i,_ON)
				self.defineList.lists[2].set(i,edit.value)
			self.get()

	# ----------------------------------------------------------------------
	# set default executable to use
	# ----------------------------------------------------------------------
	def defaultExecutable(self):
		run = self.page.runList.getRun(self.page.selected[0])
		run.exe = ""
		self.page.setModified()
		self.exe.config(text=run.exe)

	# ----------------------------------------------------------------------
	# set executable to use
	# ----------------------------------------------------------------------
	def loadExecutable(self, event=None):
		if self.exe["state"] == DISABLED: return
		if not self.page.selected: return
		run = self.page.runList.getRun(self.page.selected[0])
		fn = bFileDialog.askopenfilename(master=self,
				initialfile=run.exe,
				title="Load executable")
		if len(fn) > 0:
			run.setExecutable(fn)
			self.page.setModified()
			self.exe.config(text=run.exe)
#			self.flair.refresh("exe", self)

	# ----------------------------------------------------------------------
	def fillQueueList(self):
		queues = tkFlair.config.options(tkFlair._BATCH_SECTION)
		try: queues.remove("spawnname")
		except ValueError: pass
		queues.sort()
		self.queueList.fill(queues)
		if self.queueList.get() == "":
			if queues: self.queueList.set(queues[0])

	# ----------------------------------------------------------------------
	def setQueueList(self):
		queue = self.project.submit
		if tkFlair.config.has_option(tkFlair._BATCH_SECTION, queue):
			self.queueList.set(queue)
		else:
			self.queueList.select(0)

	# ----------------------------------------------------------------------
	# set Random parameter in the Run Frame
	#   from Fluka manual:	WHAT(2) : any number < 9.E8 (see Note 5)  Default: 54217137
	# ----------------------------------------------------------------------
	def getRandomNumberAsSeed(self):
		if self.rnd["state"] != NORMAL: return
		for r in self.page.runList.curselection():
			run = self.page.runList.getRun(r)
			if run.name == Project.DEFAULT_INPUT: continue
			self.project.setRunRandomSeed(run)

		s = self.rnd["state"]
		self.rnd.config(state=NORMAL)
		self.rnd.delete(0, END)

		run = self.page.runList.getRun(ACTIVE)
		self.rnd.insert(0, run.rnd)

		self.rnd.config(state=s)

	# ----------------------------------------------------------------------
	# disable overrides
	# ----------------------------------------------------------------------
	def setState(self, state, all=False):
		# Run / Info frame
		#self.renameButton.config(state=state)
		#self.upButton.config(state=state)
		#self.downButton.config(state=state)
		#self.cloneButton.config(state=state)

		# Override option frame
		self.title.config(state=state)
		self.defineList.config(state=state)
		self.bd.config(state=state)
		self.rnd.config(state=state)
#		self.rndButton.config(state=state)
		self.primaries.config(state=state)
		self.time.config(state=state)
		self.bed.config(state=state)
		self.bel.config(state=state)
		self.exe.config(state=state)
		if state==DISABLED:
			self.exe["background"] = "LightGray"
		else:
			self.exe["background"] = "White"

		# Enable the remaining buttons if not requested
		if not all: state = NORMAL
#		self.cont.config(state=state)
		self.prevbox.config(state=state)
		self.cyclesbox.config(state=state)
		self.lastbox.config(state=state)

	# ----------------------------------------------------------------------
	def toggleField(self, event):
		if len(self.page.selected)==1 and self.page.selected[0]==0:
			return

		if event.widget["state"] == NORMAL:
			state = DISABLED
		else:
			state = NORMAL

		if event.widget in (self.rnd, self.rndButton):
			state = NORMAL
			self.rnd.config(state=state)
			self.rndButton.config(state=state)

		elif event.widget in (self.exe, self.bed, self.bel):
			self.exe.config(state=state)
			self.bed.config(state=state)
			self.bel.config(state=state)

		elif isinstance(event.widget.master, tkExtra.Combobox):
			event.widget.master.config(state=state)

		elif event.widget in self.defineList.lists:
			self.defineList.config(state=state)

		else:
			event.widget.config(state=state)

	# ----------------------------------------------------------------------
	# select item
	# ----------------------------------------------------------------------
	def select(self, event=None):
		if self.page.selected: self.get()

		selected = list(map(int, self.page.runList.curselection()))
		if selected == self.page.selected and self.page.valid(): return
		if len(selected) == 0: return
		self.page.selected = selected

		runs = [self.page.runList.getRun(x) for x in self.page.selected]
		run = runs[0]
		if run is None: return

		# Update Information information
		# enable first to allow inserting info
		self.setState(NORMAL)
		self._varSetting = True
		self.last.set(run.last)
		self.prev.set(run.prev)
		self.cycles.set(run.last-run.prev)
		self._varSetting = False

		self.title.delete(0, END)

		self.rnd.delete(0, END)
		self.rnd.insert(0, run.rnd)

		self.primaries.delete(0, END)
		self.primaries.insert(0, run.primaries)

		self.time.delete(0, END)
		self.time.insert(0, run.time)

		if self.page.selected[0] == 0:
			self.title.insert(0, self.project.title)
			self.exe.config(text=self.project.exe)
			self.selectDefaultDefines()
		else:
			self.title.insert(0, run.title)
			self.exe.config(text=run.exe)
			self.selectDefaultDefines(runs)

		self.refresh()

	# ----------------------------------------------------------------------
	# Change in the cycles
	# ----------------------------------------------------------------------
	def prevChange(self, a, b, c):
		if self._varSetting: return
		self._varSetting = True
		try: self.last.set(self.prev.get() + self.cycles.get())
		except: pass
		self._varSetting = False

	# ----------------------------------------------------------------------
	# Change in the cycles
	# ----------------------------------------------------------------------
	def lastChange(self, a, b, c):
		if self._varSetting: return
		self._varSetting = True
		try: self.cycles.set(self.last.get() - self.prev.get())
		except: pass
		self._varSetting = False

	# ----------------------------------------------------------------------
	# Update all fields
	# ----------------------------------------------------------------------
	def get(self, event=None):
		#print "Run.get", self.page.selected, self.page.invalid()
		# executable can change from Compile
		# do not save anything if invalid
		if len(self.page.selected)==0 or self.page.invalid():
			return

		try:
			s = self.queueList.get()
			if self.project.submit != s:
				self.project.submit = s
				self.page.setModified()
		except TclError:
			return

		for sel in self.page.selected:
			run = self.page.runList.getRun(sel)
			if run is None:
				continue

			oldhash = hash(run)

			family = [self.project.getRunByName(x,run.name) for x in run.family ]

			if self.prevbox["state"] == NORMAL:
				try: run.prev = self.prev.get()
				except: pass
				for r in family: r.prev = run.prev

			if self.cyclesbox["state"]==NORMAL or self.lastbox["state"]==NORMAL:
				try: run.last = self.last.get()
				except: pass
				for r in family: r.last = run.last

			if self.exe["state"] == NORMAL:
				if sel==0:
					self.project.exe = self.exe["text"]
				else:
					run.exe = self.exe["text"]
				for r in family: r.exe = run.exe


			# skip the rest for the main run
			if sel != 0:
				if self.defineList["state"] == NORMAL:
					olddefs = run.defines
					run.defines = []
					for s,n,v in self.defineList.get(0,END):
						if s==_ON:
							# add define
							run.defines.append((n,v))
						elif s==_MULTI:
							# Append the existing one if any
							for nn,vv in olddefs:
								if n==nn:
									run.defines.append((nn,vv))
									break
					for r in family: r.defines = run.defines

				if self.title["state"] == NORMAL:
					run.title = self.title.get().strip()
					for r in family: r.title = run.title

				if self.rnd["state"] == NORMAL:
					try:	run.rnd   = int(self.rnd.get())
					except: run.rnd   = 0
					# No propagation of rnd in the family

				if self.primaries["state"] == NORMAL:
					try:	run.primaries = float(self.primaries.get())
					except: run.primaries = 0
					for r in family: r.primaries = run.primaries

				if self.time["state"] == NORMAL:
					try:	run.time = float(self.time.get())
					except: run.time = 0
					for r in family: r.time = run.time

			if hash(run) != oldhash:
				#self.flair.redraw() or refresh()???
				self.page.setModified()

	# ----------------------------------------------------------------------
	# Continue from last run
	# ----------------------------------------------------------------------
	def continueFromLast(self):
		# Scan directory to find last random number seed
		for r in self.page.runList.curselection():
			run = self.page.runList.getRun(r)
			last = run.lastRun()
			run.last = last + (run.last - run.prev)
			run.prev = last
			if self.page.selected[0] == int(r): self.prev.set(last)

	# ----------------------------------------------------------------------
	# Start all runs
	# ----------------------------------------------------------------------
	def execute(self):
		self.page.busy()
		self.get()
		for r in self.page.runList.curselection():
			run = self.page.runList.getRun(r)
			if run.family:
				for child in run.family:
					r = self.project.getRunByName(child, run.name)
					try:
						rc = self._startRun(r)
					except (OSError, IOError):
						self.flair.notify("Error",
							sys.exc_info()[1],
							tkFlair.NOTIFY_ERROR)
						rc = 1
						break
			else:
				try:
					rc = self._startRun(run)
				except (OSError, IOError):
					self.flair.notify("Error",
						sys.exc_info()[1],
						tkFlair.NOTIFY_ERROR)
					rc = 1
					break

			if rc:
			#	messagebox.showerror("Run",
			#		"Run submission stopped due to errors on run: %s"%(run.getInputName()),
			#		parent=self)
				break

		self.setState(DISABLED, True)
		self.page.setModified()
		if self.timer is None: self._timer()

		# Show present information
		self.refresh()
		self.page.notBusy()

	# ----------------------------------------------------------------------
	# Start one run and try to attach
	# ----------------------------------------------------------------------
	def _startRun(self, run):
		if run.status == Project.STATUS_WAIT2ATTACH or \
		   run.status == Project.STATUS_RUNNING: return

		if self.project.inputName == "":
			messagebox.showerror("Input file not saved",
				"Please save the input file before",
				parent=self.ribbon)
			return True

		try:
			log = self.flair.newLog("Run", run.getInputBaseName())
			rc = run.start(log=log)
			if rc:
				self.flair.notify("Run Errors",
					"Errors or warnings during writing run input.",
					tkFlair.NOTIFY_ERROR, "Output", self.flair.showOutput)
			self.flair.restoreLog()
			return rc
		except Project.RunException:
			self.flair.notify("Error submitting run",
				sys.exc_info()[1],
				tkFlair.NOTIFY_ERROR)
			return True

	# ----------------------------------------------------------------------
	def activate(self):
		self.select()
		self.refreshStatus()
		if self.timer is None: self._timer()

	# ----------------------------------------------------------------------
	# Clean files
	# ----------------------------------------------------------------------
	def clean(self):
		if tkFlair.askyesno("Warning",
				"Do you want to remove all generated files?",
				parent=self.ribbon):
			log = self.flair.newLog("Run", "Clean")
			warn = False
			for r in self.page.runList.curselection():
				run = self.page.runList.getRun(r)
				if run.family:
					for child in run.family:
						r = self.project.getRunByName(child, run.name)
						if not self._cleanRun(r, log):
							warn = True
							break
				if not self._cleanRun(run, log):
					warn = True
					break

			if warn:
				self.flair.notify("Clean", "Run files removed with warnings", tkFlair.NOTIFY_WARNING)
			else:
				self.flair.notify("Clean", "Run files removed")
			self.page.setStatus("Run files cleaned","DarkBlue")

	# ----------------------------------------------------------------------
	def _cleanRun(self, run, log):
		if run.status in (Project.STATUS_WAIT2ATTACH, Project.STATUS_RUNNING):
			messagebox.showwarning("Running",
				"Input '%s' is running. Please stop before deleting all files"%(run.name),
				parent=self.ribbon)
			return False

		if run.clean():
			self.flair.notify("Errors during run cleaning",
				"Please check the output for error messages during cleaning",
				tkFlair.NOTIFY_ERROR)
		return True

	# ----------------------------------------------------------------------
	def _flukaStop(self, stopfile):
		for r in self.page.runList.curselection():
			self.page.runList.getRun(r)._stop(stopfile)

	# ----------------------------------------------------------------------
	# Stop smoothly all runs
	# ----------------------------------------------------------------------
	def stopCycle(self, msg=True):
		if not tkFlair.askyesno("Stop Running Cycle(s)",
				"Are you sure you want to stop the running "
				"cycle of the selected run(s)",
				parent=self.ribbon):
			return
		self._flukaStop("fluka.stop")

	# ----------------------------------------------------------------------
	# Stop smoothly all runs
	# ----------------------------------------------------------------------
	def stopRun(self, msg=True):
		if not tkFlair.askyesno("Stop Run(s)",
			"Are you sure you want to stop the selected run(s)",
			parent=self.ribbon):
			return
		self._flukaStop("rfluka.stop")

	# ----------------------------------------------------------------------
	# Kill all runs
	# ----------------------------------------------------------------------
	def killRun(self):
		if not tkFlair.askyesno("Stop Run(s)",
			"Are you sure you want to KILL the selected run(s)",
			parent=self.ribbon):
			return
		for r in self.page.runList.curselection():
			run = self.page.runList.getRun(r)
			# Check status, if not running start the run
			if run.status == Project.STATUS_RUNNING and run.pid>0:
				# Create a file fluka.stop
				#f = open("%s/fluka_%d/fluka.stop" \
				#	% (run.getDir(), run.pid),"w")
				#f.close()
				# and kill the rfluka command
				if Project.kill == "":
					try:
						# For Linux
						os.spawnv(os.P_WAIT, "/usr/bin/pkill",
							[ "/usr/bin/pkill", "-P", str(run.pid) ])
						os.kill(run.pid, signal.SIGKILL)
					except OSError:
						messagebox.showwarning("Kill Failed",
							"Kill Run has failed. Probably is "
							"running thru a batch system. "
							"Please kill it manually.",
							parent=self.ribbon)
				else:
					try:
						os.spawnv(os.P_NOWAIT, Project.kill,
							[ Project.kill, str(run.pid) ] )
					except:
						messagebox.showerror("Error executing command",
							"Error executing command:\n%s %d" \
							% (Project.kill, run.pid),
							parent=self.ribbon)
			run.status = Project.STATUS_FINISHED_ERR
			run.pid    = 0
		self.refresh()
		self.setState(DISABLED)

	# ----------------------------------------------------------------------
	# Update run progress information
	# and Enable or Disable the fields
	# ----------------------------------------------------------------------
	def refresh(self):
		if self._updating: return
		self._updating = True

		# Update messages
		self.pstatus["text"]	= ""
		self.pinput["text"]	= ""
		self.pdir["text"]	= ""
		self.ptime["text"]	= ""
		self.pstarted["text"]	= ""
		self.pelapsed["text"]	= ""
		self.premain["text"]	= ""
		self.peta["text"]	= ""
		self.prunremain["text"] = ""

		# Scan status of all runs
		changed = False
		for run in self.project.runs:
			s = run.status
			if run.status == Project.STATUS_WAIT2ATTACH or \
			  (run.status == Project.STATUS_RUNNING and run.pid==0):
				if run.family:
					run.status = Project.STATUS_NOT_RUNNING
					changed = True
					for child in run.family:
						rchild = self.project.getRunByName(child, run.name)
						if rchild.status in (Project.STATUS_NOT_RUNNING,
								     Project.STATUS_WAIT2ATTACH):
							self.attach(rchild)
				else:
					self.attach(run)
			elif run.status == Project.STATUS_RUNNING:
				if run.family:
					run.status = Project.STATUS_NOT_RUNNING
					changed = True
					for child in run.family:
						rchild = self.project.getRunByName(child, run.name)
						rchild.refresh()
				else:
					run.refresh()
			if run.status != s: changed = True

		def setRunStatus(run, status):
			if run.status != status:
				run.status = status
				changed = True

		# Update families colors
		for run in self.project.runs:
			if not run.family: continue

			# Operate only on families
			count = ZeroIntDict()
			for child in run.family:
				rchild = self.project.getRunByName(child, run.name)
				count[rchild.status] += 1
			n = len(run.family)

			# Check if all childs are the same
			if len(count)==1:
				# check against last child status
				setRunStatus(run, rchild.status)
			else:
				if Project.STATUS_FINISHED_ERR in count:
					setRunStatus(run, Project.STATUS_FINISHED_ERR)

				elif Project.STATUS_RUNNING in count:
					setRunStatus(run, Project.STATUS_RUNNING)

		# update colors on listbox
		if changed: self.page.runList.updateColors()

		# Show active run
		if len(self.page.selected)==0 or \
		   self.page.selected[0] >= len(self.project.runs):
			self._updating = False
			return

		run = self.page.runList.getRun(self.page.selected[0])
		if run is None:
			self._updating = False
			return

		# Find information
		if not run.family and (run.pid>0 or run.status in (Project.STATUS_RUNNING, Project.STATUS_TIMEOUT)):
			run.refresh(True)

		self.pstatus["text"] = STATUS_STR[run.status]
		self.pinput["text"]  = run.getInputName()

		# Show information for a family
		#if run.family:
		#	startTime = 1e99
		#	for x in run.family:
		#		rchild = self.project.getRunByName(x,run.name)
		#		if rchild.pid>0 or run.status in (Project.STATUS_RUNNING, Project.STATUS_TIMEOUT):
		#			rchild.refresh(True)
		#		if rchild.startTime>0:
		#			startTime = min(startTime, rchild.startTime)
		#	if startTime < 1e99:
		#		run.startTime = startTime
		#		self.pstarted["text"] = time.strftime(
		#			bFileDialog._TIME_FORMAT,
		#			time.localtime(run.startTime))
		#	else:
		#		run.startTime = 0

		# Individual run
		if not run.family and run.status in (Project.STATUS_RUNNING, Project.STATUS_TIMEOUT):
			self.setState(DISABLED, True)
			self.pdir["text"] = "fluka_%d" % (run.pid)
			if run.cycle>0:
				ncycles = run.last - run.prev
				self.cyclePBar.setLimits(run.prev, run.last)
				self.cyclePBar.setProgress(run.cycle)

			else:
				self.cyclePBar.setProgress(0.0, 0.0, "")

			nprim = run.handled + run.remaining
			if run.handled>0:
				self.ptime["text"] = Utils.friendlyTime(run.timeperprim)

				self.primaryPBar.setLimits(0.0, nprim)
				self.primaryPBar.setProgress(run.handled+1)

				self.pelapsed["text"] = Utils.friendlyTime(
					run.handled*run.timeperprim)

				cycle_remain = min(run.remaining*run.timeperprim,
						run.maxtime - run.handled*run.timeperprim)
				self.premain["text"] = Utils.friendlyTime(cycle_remain)

				run_remain = cycle_remain + (run.last-run.cycle) * \
						(min(nprim*run.timeperprim, run.maxtime) + run.initime)
				self.prunremain["text"] = Utils.friendlyTime(run_remain)
			else:
				self.primaryPBar.setProgress(0.0, 0.0, "")

			if run.startTime>0:
				self.pstarted["text"] = time.strftime(
					bFileDialog._TIME_FORMAT,
					time.localtime(run.startTime))
				self.peta["text"] = time.strftime(
					bFileDialog._TIME_FORMAT,
					time.localtime(time.time() + \
					  (run.remaining + (run.last-run.cycle)*nprim) \
					  * run.timeperprim))

		else:
			# ------------------------------------
			# Finally enable/disable fields
			# ------------------------------------
			if len(self.page.selected)==1:
				if self.page.selected[0]==0 or run.parent != "":
					# Disable if child apart from rnd
					self.setState(DISABLED)
					# Apart from exe
					self.exe["background"] = "White"
					self.exe["state"] = NORMAL
					self.bed["state"] = NORMAL
					self.bel["state"] = NORMAL
					# and maybe the run
					if run.parent != "":
						self.rnd.config(state=NORMAL)
				else:
					self.setState(NORMAL)
			elif len(self.page.selected)>1:
				self.setState(DISABLED, True)
			else:
				self.setState(NORMAL)
			self.cyclePBar.setProgress(0.0, 0.0, "")
			self.primaryPBar.setProgress(0.0, 0.0, "")

		self._updating = False

	# ----------------------------------------------------------------------
	# Refresh Status
	# ----------------------------------------------------------------------
	def refreshStatus(self):
		if self.flair is None: return
		# Count processes
		running = 0
		total   = 0
		for run in self.project.runs:
			if run.family: continue	# Skip groups
			if run.status == Project.STATUS_RUNNING: running += 1
			total += 1
		self.page.setStatus("Running %d out of %d" \
			% (running, total))

	# ----------------------------------------------------------------------
	# Attach to run
	# ----------------------------------------------------------------------
	def attach2Run(self):
		self.page.busy()
		for r in self.page.runList.curselection():
			run = self.page.runList.getRun(r)
			run.exclude = None
			if run.family:
				for child in run.family:
					r = self.project.getRunByName(child, run.name)
					r.exclude = None
					r.status = Project.STATUS_WAIT2ATTACH
			else:
				run.status = Project.STATUS_WAIT2ATTACH
			self.refresh()
		self.page.notBusy()

	# ----------------------------------------------------------------------
	# Attach to a process
	# ----------------------------------------------------------------------
	def attach(self, run):
		found = run.attach()
		if isinstance(found,list):
			# Multiple choices
			# Popup a list with the possible choices
			aid = SelectRunningDir(self,run,found).show()
			if aid<0: return

			run.status = Project.STATUS_RUNNING
			run.pid    = int(Project._PAT_DIR2.match(found[aid][0]).group(2))
			run.attachTime = time.time()

		self.page.setModified()
		#stat = os.stat("%s/%s" % (run.getDir(), found[aid][0]))

	# ----------------------------------------------------------------------
	# Timer function
	# ----------------------------------------------------------------------
	def _timer(self):
		self.timer = None
		self.refresh()
		# Set a new timer
		if not self.winfo_ismapped(): return
		self.timer = self.after(Project.refreshInterval, self._timer)

	# ----------------------------------------------------------------------
	# update progress bars
	# ----------------------------------------------------------------------
	def _resize(self, event):
		width  = event.width
		height = event.height

		p1, p2 = self.pcycle_perc
		box1, box2 = self.pcycle_jauge

		if box1:
			w = int(width*p1)
			self.pcycle_canvas.coords(box1, 0, 0, w, height)
		w = int(width*p2)
		self.pcycle_canvas.coords(box2, 0, 0, w, height)

		p1, p2 = self.pprim_perc
		box1, box2 = self.pprim_jauge
		if box1:
			w = int(width*p1)
			self.pprim_canvas.coords(box1, 0, 0, w, height)
		w = int(width*p2)
		self.pprim_canvas.coords(box2, 0, 0, w, height)

	# ----------------------------------------------------------------------
	def help(self):
		Manual.show(":Run:Run:")

#===============================================================================
# Select Running Directory
#===============================================================================
class SelectRunningDir(Toplevel):
	def __init__(self, master, run, dirlist):
		Toplevel.__init__(self, master)

		l = Label(self, text="Run: %s"%(run.name),
			justify=RIGHT, anchor=NW, pady=2)
		l.pack(side=TOP)

		l = Label(self, text="Directories:",
			justify=RIGHT, anchor=NW, pady=2)
		l.pack(side=TOP)

		self.listbox = tkExtra.ExListbox(self, borderwidth=0,
			selectmode=BROWSE,
			selectborderwidth=0,
			background="White",
			takefocus=True,
			exportselection=FALSE)
		self.listbox.pack(expand=YES,fill=BOTH)

		frame = Frame(self)
		Button(frame, text="Cancel", command=self.cancel).pack(side=RIGHT)
		Button(frame, text="Ok", command=self.ok).pack(side=RIGHT)
		self.listbox.bind("<Double-1>",self.ok)
		self.bind("<Return>",   self.ok)
		self.bind("<KP_Enter>", self.ok)
		self.bind("<Escape>",   self.cancel)
		frame.pack(side=BOTTOM)

		# Populate directories
		for d,f in dirlist:
			self.listbox.insert(END,d)

		self.listbox.activate(0)
		self.listbox.selection_set(0)

		self.sel = -1
		self.title("Select Running Dir")

	# ----------------------------------------------------------------------
	def show(self):
		self.focus_set()
		self.listbox.focus_set()
		self.transient(self.master)
		self.wait_window()
		return self.sel

	# ----------------------------------------------------------------------
	def ok(self, event=None):
		self.sel = self.listbox.index(ACTIVE)
		self.destroy()

	# ----------------------------------------------------------------------
	def cancel(self, event=None):
		self.destroy()
