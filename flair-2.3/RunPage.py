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
from stat import *

try:
	from tkinter import *
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	import tkinter.messagebox as messagebox

import Input
import Manual
import Ribbon
import Project
import RunList
import tkFlair
import bFileDialog
import FlairRibbon

import tkExtra
import tkDialogs

import RunTab
import DataTab
import FilesTab

STATUS_STR = [
	"Not Running",
	"Waiting to attach",
	"Running",
	"Finished OK",
	"Finished with ERRORS",
	"*** TIMED-OUT ***" ]

#===============================================================================
# Run Frame
#===============================================================================
class RunPage(FlairRibbon.FlairPage):
	"""Execute a simulation, inspect the output files and merge the generated data."""

	_name_ = "Run"
	_icon_ = "run"

	#----------------------------------------------------------------------
	def init(self):
		self.selected   = []
		self._refreshOn = ("run","exe","files","data")

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Tools ===========
		group = Ribbon.LabelGroup(self.ribbon, "View")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		kw = {
			"compound" : TOP,
			"variable" : self._tabVar,
			"command"  : self.changeTab,
			"width"    : 48,
		}

		# ---
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["run32"],
				text="Runs",
				value=0,
				**kw)
		b.pack(side=LEFT, fill=Y, expand=YES)
		tkExtra.Balloon.set(b, "Run setting page")

		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["file-manager32"],
				text="Files",
				value=1,
				**kw)
		tkExtra.Balloon.set(b, "Output files inspecting page")
		b.pack(side=LEFT, fill=Y, expand=YES)

		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["data32"],
				text="Data",
				value=2,
				**kw)
		tkExtra.Balloon.set(b, "Merge output data page")
		b.pack(side=LEFT, fill=Y, expand=YES)

	#----------------------------------------------------------------------
	# Create Project page
	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self, False)

		split = PanedWindow(self.frame, orient=HORIZONTAL,
			showhandle=0,
			sashrelief=RAISED, sashwidth=4, opaqueresize=1)
		split.pack(expand=YES, fill=BOTH)

		# ======= Runs ========
		self.runList = RunList.RunListbox(split, True, height=5, stretch="middle")
		split.add(self.runList, minsize=50) #, padx=3, pady=3)

		# Bind events
		self.runList.bind('<<ListboxSelect>>',		self.select)
		self.runList.bindList("<Insert>",		self.add)
		self.runList.bindList("<Delete>",		self.delete)
		self.runList.bindList("<Control-Key-Up>",	self.moveUp)
		self.runList.bindList("<Control-Key-Down>",	self.moveDown)

		# ======== Tabs =======
		tabFrame = Frame(split, padx=0, pady=0)
		split.add(tabFrame)

		# --------- Run tab ---------
		self.run = RunTab.RunTab(tabFrame, self)
		self.addTab(self.run)
		# --------- Files tab ---------
		self.addTab(FilesTab.FilesTab(tabFrame, self))
		# --------- Data tab ---------
		self.addTab(DataTab.DataTab(tabFrame, self))

		self.runList.set(self.flair)

	# ----------------------------------------------------------------------
	def select(self, event=None):
		self.activeTab.activate()
		self.setValid()

	# ----------------------------------------------------------------------
	def get(self):
		self.run.get()

	# ----------------------------------------------------------------------
	# Activate and set focus
	# ----------------------------------------------------------------------
	def activate(self):
		FlairRibbon.FlairPage.activate(self)

	# ----------------------------------------------------------------------
	# Update run progress information
	# and Enable or Disable the fields
	# ----------------------------------------------------------------------
	def refresh(self):
		FlairRibbon.FlairPage.refresh(self)

		if self.project is None: return

		# Refresh run listbox
		self.run.fillQueueList()	# If list is changed from ConfigDialog
		self.run.setQueueList()

		# Run list
		self.runList.fill()
		self.select()

	# ----------------------------------------------------------------------
	def refreshStatus(self):
		FlairRibbon.FlairPage.refreshStatus(self)
		self.activeTab.refreshStatus()

	# ----------------------------------------------------------------------
	def reset(self):
		self.selected = []

	# ----------------------------------------------------------------------
	def beforeDock(self):
		self.reset()

	# ----------------------------------------------------------------------
	def execute(self):
		self.activeTab.execute()

	# ----------------------------------------------------------------------
	def clean(self):
		self.activeTab.clean()

	# ----------------------------------------------------------------------
	# Move runs to a different position
	# expect a list of runs (idx, run)
	# ----------------------------------------------------------------------
	def _moveRunsToPosition(self, runs, pos):
		self.runList.updateSelection()

		# Delete selected runs from the list
		runs.reverse()
		for run in runs:
			self.project.runs.remove(run)

		# Insert them at pos
		for run in runs:
			append = (pos >= len(self.project.runs))
			if append:
				self.project.runs.append(run)
			else:
				self.project.runs.insert(pos, run)

		self.runList.fill()
		self.runList.activate(pos)
		self.select()
		self.flair.refresh("run")
		self.setModified()

	# ----------------------------------------------------------------------
	# Move up select items by one
	# ----------------------------------------------------------------------
	def moveUp(self, event=None):
		self.get()
		del self.selected[:]
		# FIXME broken for spawned jobs

		runs = self._selectedParentsList()
		if not runs: return

		# Find position to insert them
		pos = self.project.getRunIndex(runs[0]) - 1
		if pos <= 0: return

		# find a parent position
		while pos>0:
			if self.project.runs[pos].parent == "": break
			pos -= 1

		if pos <= 0: return
		self._moveRunsToPosition(runs, pos)
		return "break"

	# ----------------------------------------------------------------------
	# Move down select items by one
	# ----------------------------------------------------------------------
	def moveDown(self, event=None):
		self.get()
		del self.selected[:]
		runs = self._selectedParentsList()
		if not runs: return

		# Find position to insert them
		pos = self.project.getRunIndex(runs[-1]) + 1
		sz  = len(self.project.runs) - 1
		if pos > sz: return

		# find a parent position
		while pos <= sz:
			if self.project.runs[pos].parent == "":
				# move to the end of the family
				if self.project.runs[pos].family:
					# skip child's
					pos += 1
					while pos <= sz:
						if self.project.runs[pos].parent == "": break
						pos += 1
				else:
					pos += 1
				break
			pos += 1

		pos -= len(runs)

		self._moveRunsToPosition(runs, pos)
		return "break"

	# ----------------------------------------------------------------------
	# Add new input to run
	# ----------------------------------------------------------------------
	def add(self, event=None):
		self.get()
		inpname = self.project.inputName
		fn = bFileDialog.asksaveasfilename(master=self.page,
			title="Create new run-input file as",
			filetypes=[
				("Fluka %s*.inp files"%(inpname),"%s*.inp"%(inpname)),
				("Fluka Input files","*.inp"),
				("All","*")])
		if len(fn) > 0:
			f = self.project.relativePath(fn)
			(fn,ext) = os.path.splitext(f)

			# Check if run is there
			newfn = self.project.findUniqueRunName(fn)
			if newfn != fn:
				messagebox.showerror("Run already exists",
					"Run '%s' already exists. Renamed to '%s'"%(fn,newfn),
					parent=self.ribbon)
				fn = newfn

			run = Project.RunInfo(self.project, fn)
			run.exe = self.project.exe	# Default executable

			self.project.runs.append(run)
			del self.selected[:]
			self.runList.fill()
			self.flair.refresh("run")
			self.runList.selection_clear(0,END)
			self.runList.selection_set(END)
			self.runList.activate(END)
			self.runList.see(END)
			self.select()
			self.run.selectDefaultDefines()
			self.setModified(True)

	# ----------------------------------------------------------------------
	# Delete input from running
	# ----------------------------------------------------------------------
	def delete(self, event=None):
		# protect against accidental delete
		#if self.flair.event is not None:
		#	if not self.focus_get() in self.runList.lists: return

		sel = self.runList.curselection()
		if len(sel) == 0: return "break"

		if tkFlair.askyesno("Remove Runs",
		   "Are you sure you want to remove the selected runs?", parent=self.ribbon):
			del self.selected[:]
			sel = list(self.runList.curselection())
			first = int(sel[0])
			sel.reverse()
			for r in sel:
				if int(r) == 0: continue
				run = self.runList.getRun(r)
				if run is not None:
					self.project.delRun(run)

			if first > self.runList.size():
				first = END

			self.runList.fill()
			self.runList.selection_clear(0,END)
			self.runList.selection_set(first)
			self.runList.activate(first)
			self.select()
			self.flair.refresh("run")
			self.setModified()
		return "break"

	# ----------------------------------------------------------------------
	# Clone runs
	# ----------------------------------------------------------------------
	def clone(self, event=None):
		self.get()
		# protect against accidental delete
		#if self.flair.event is not None:
		#	if not self.focus_get() in self.runList.lists: return

		sel = self.runList.curselection()
		if len(sel) == 0: return "break"
		oldSize = self.runList.size()
		for r in sel:
			run = self.runList.getRun(r)

			# Cannot clone children
			if run.parent != "": continue

			child = run.clone()

			if run.name == Project.DEFAULT_INPUT:
				child.name = self.project.inputName+"-1"
				child.defaultDefines()

			child.name = self.project.findUniqueRunName(child.name)
			self.project.runs.append(child)

		self.runList.fill()
		self.flair.refresh("run")
		self.runList.selection_clear(0,END)
		self.runList.selection_set(oldSize,END)
		self.runList.activate(END)
		self.runList.see(END)
#		self.selected = map(int,self.runList.curselection())
		self.select()
		self.setModified()
		return "break"

	# ----------------------------------------------------------------------
	# Loop runs
	# ----------------------------------------------------------------------
	def loop(self):
		runs = []
		for r in self.runList.curselection():
			run = self.runList.getRun(r)
			if run.parent != "": continue
			runs.append(run)

		if not runs: return

		dlg = LoopRunsDialog(self.page, self.flair, runs)
		self.flair.refresh("run")
		self.runList.fill()

	# ----------------------------------------------------------------------
	# Rename selected runs
	# ----------------------------------------------------------------------
	def rename(self, event=None):
		self.get()
		sel = list(self.runList.curselection())
		if len(sel) == 0: return
		sel.reverse()
		for r in sel:
			# Skip default and childs
			if int(r)==0: continue
			run = self.runList.getRun(r)
			if run.parent != "": continue

			newname = bFileDialog.asksaveasfilename(master=self.page,
				title="Rename run '%s' to"%(run.name),
				initialfile=run.name,
				filetypes=[("Fluka Input files","*.inp"),
					("All","*")])
			if len(newname) == 0: continue
			f = self.project.relativePath(newname)
			(fn, ext) = os.path.splitext(f)

			self.runList.renameRun(run, fn)
			self.runList.fill()
			self.flair.refresh("run")
			self.setModified()

	# ----------------------------------------------------------------------
	# Return a list of the selected parents only
	# ----------------------------------------------------------------------
	def _selectedParentsList(self):
		runs = []
		for r in self.runList.curselection():
			if int(r)==0: continue		# ignore default run
			run = self.runList.getRun(r)
			if run.parent != "": continue	# ignore child nodes
			runs.append(run)
			if run.family:
				# Append all its family
				for childname in run.family:
					child = self.project.getRunByName(childname, run.name)
					runs.append(child)

		return runs

	# ----------------------------------------------------------------------
	def help(self):
		if str(self.page.focus_get()).startswith(str(self.runList)):
			FlairRibbon.FlairPage.help(self)
		else:
			self.activeTab.help()

#===============================================================================
# Create clone of Runs by looping on one variable
#===============================================================================
class LoopRunsDialog(Toplevel):
	def __init__(self, master, flair, runs):
		Toplevel.__init__(self, master)
		self.transient(master)
		self.title("Loop runs")

		self.flair = flair
		self.runs  = runs

		self.incparent = BooleanVar(self)

		l = Label(self, text=	"Create clone of the selected runs by\n" \
					"looping the following variable over\n" \
					"the defined range",
				justify=LEFT,
				anchor=NW,
				pady=2)
		l.pack(side=TOP)

		frame = Frame(self)
		frame.pack(fill=BOTH)

		# ---
		row = 0
		Label(frame, text="Variable:").grid(
			row=row, column=0, sticky=E)

		self.variable = tkExtra.Combobox(frame, width=16)
		self.variable.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.variable,"Select the variable to loop")

		# Update defines
		lst = []
		for card in self.flair.project.input.cardsSorted("#define",2):
			if card._indent != 0: continue
			lst.append(card.sdum())
		self.variable.fill(lst)

		# ---
		row += 1
		Label(frame, text="Start value:").grid(
			row=row, column=0, sticky=E)
		self.start = Entry(frame, background="White", justify=LEFT)
		self.start.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.start, "Enter starting value or expression")

		# ---
		row += 1
		Label(frame, text="End value:").grid(
			row=row, column=0, sticky=E)
		self.stop = Entry(frame, background="White", justify=LEFT)
		self.stop.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.stop, "Enter ending value or expression")

		# ---
		row += 1
		Label(frame, text="Step:").grid(
			row=row, column=0, sticky=E)
		self.step = Entry(frame, background="White", justify=LEFT)
		self.step.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.step, "Step value or expression")

		# ---
		row += 1
		b = Checkbutton(frame, text="Use Parent as 1st iteration",
				variable=self.incparent,
				onvalue=1, offvalue=0,
				anchor=W)
		tkExtra.Balloon.set(b,
			"Check to use the parent run as the first iteration of the loop.\n" \
			"Uncheck to create a new run (suffix -1) for the first iteration")
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		for run in runs:
			if run.name == Project.DEFAULT_INPUT:
				b.config(state=DISABLED)
				break

		# ---
		frame.grid_columnconfigure(1, weight=1)

		# ----
		frame = Frame(self)
		Button(frame, text="Cancel",
				image=tkFlair.icons["x"],
				compound=LEFT,
				command=self.cancel).pack(side=RIGHT)
		Button(frame, text="Ok",
				image=tkFlair.icons["ok"],
				compound=LEFT,
				command=self.ok).pack(side=RIGHT)
		self.bind("<Return>",  self.ok)
		self.bind("<KP_Enter>",self.ok)
		self.bind("<Escape>",  self.cancel)
		frame.pack(side=BOTTOM, fill=X)

		self.deiconify()
		self.wait_visibility()
		# Get the focus
		self.focus_set()
		self.grab_set()
		self.wait_window()

	# ----------------------------------------------------------------------
	def ok(self, event=None):
		var = self.variable.get()
		if var=="":
			messagebox.showerror("No variable selected",
				"Please select a valid variable first",
				parent=self)
			return

		inp = self.flair.project.input
		try:
			start = float(eval(self.start.get(),
					Input._globalDict,
					inp.localDict))
		except:
			messagebox.showerror("Error start value",
				"Please enter a valid starting value",
				parent=self)
			return

		try:
			stop = float(eval(self.stop.get(),
					Input._globalDict,
					inp.localDict))
		except:
			messagebox.showerror("Error stop value",
				"Please enter a valid stopping value",
				parent=self)
			return

		step = self.step.get()
		if step == "": step="1"
		try:
			step = float(eval(step,
					Input._globalDict,
					inp.localDict))
		except:
			messagebox.showerror("Error step value",
				"Please enter a valid step value",
				parent=self)
			return

		if step==0.0: step = 1.0

		if start<stop and step<0.0:
			messagebox.showerror("Wrong step value",
				"Start<Stop and step<0. Please change sign of step",
				parent=self)
			return
		elif start>stop and step>0.0:
			messagebox.showerror("Wrong step value",
				"Start>Stop and step>0. Please change sign of step",
				parent=self)
			return

		# Find run entry positions
		self.findPositions()

		# Create runs
		val   = start
		first = True
		if step > 0.0:
			while val <= stop:
				self.createRuns(var, val, first)
				val += step
				first = False

		else:
			while val >= stop:
				self.createRuns(var, val, first)
				val += step
				first = False

		self.flair.setModified()
		self.destroy()

	# ----------------------------------------------------------------------
	def cancel(self, event=None):
		self.destroy()

	# ----------------------------------------------------------------------
	def findPositions(self):
		project = self.flair.project
		sz  = len(project.runs) - 1
		self.pos = []
		for run in self.runs:
			pos = project.getRunIndex(run) + 1
			while pos<sz and project.runs[pos].parent!="":
				pos += 1
			self.pos.append(pos)

	# ----------------------------------------------------------------------
	def createRuns(self, var, val, first):
		try:
			f = float(val)
			if float(int(f)) == f:
				val = str(int(f))
			else:
				val = str(f)
		except:
			val = str(val)
		project = self.flair.project
		for i,run in enumerate(self.runs):
			if first and self.incparent.get():
				clone = run
			else:
				clone = run.clone()
				if run.name == Project.DEFAULT_INPUT:
					clone.name = run.getInputBaseName()
					clone.defaultDefines()
				clone.name = project.findUniqueRunName(clone.name+"-1")

				if clone.title == "":
					title = project.title
				else:
					title = clone.title
				clone.title = "%s [%s=%s]"%(title, var, val)

			# Check defines and append the var,value
			for j,nv in enumerate(clone.defines):
				if nv[0]==var:
					del clone.defines[j]
					break
			clone.defines.append((var,val))

			if not (first and self.incparent.get()):
				project.runs.insert(self.pos[i],clone)

				# increase all positions higher than us by 1
				for k in range(i,len(self.pos)):
					self.pos[k] += 1
