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
import Project
import tkFlair
import tkExtra

_CHILD = "."

#===============================================================================
# Run Listbox
#===============================================================================
class RunListbox(tkExtra.MultiListbox):
	def __init__(self, master=None, edit=False, **kw):
		tkExtra.MultiListbox.__init__(self, master,
				(("+",    1, None),
				 ('Run', 20, None),
				 ('Spawn',3, None)),
				 **kw)
		self.sortAssist = None

		# Bind events
		self.bindList("<space>",		self.toggle)
		self.bindList("<F2>",			self.renameInPlace)
		self.bindList("<ButtonRelease-1>",	self.buttonRelease)
		self.bindLeftRight()

		self.lists[1].bind("<Return>",		self.renameInPlace)
		self.lists[1].bind("<KP_Enter>",	self.renameInPlace)
		self.lists[1].bind("<Double-1>",	self.renameInPlace)
		self.lists[1].bind("<Control-Key-e>",	self.renameInPlace)

		self.lists[2].bind("<F2>",		self.spawnInPlace)
		self.lists[2].bind("<Return>",		self.spawnInPlace)
		self.lists[2].bind("<KP_Enter>",	self.spawnInPlace)
		self.lists[2].bind("<Button-1>",	self.spawnButton1)

		self.setPopupMenu([
#			("Expand",     0, self.expand,   tkFlair.icons["+"]),
#			("Collapse",   0, self.collapse, tkFlair.icons["-"]),
#			None,
#			("Add",        0, self.add,           tkFlair.icons["add"]),
#			("Remove",     0, self.remove,        tkFlair.icons["x"]),
#			("Clone",      0, self.clone,         tkFlair.icons["clone"]) ])
			("Rename",     0, self.renameInPlace, tkFlair.icons["rename"]) ])

		# Remember selected runs
		self.edit     = edit
		self.selected = []
		self.expanded = []
		self.flair    = None
		self.project  = None
		self.lastrelease = -1

	# ----------------------------------------------------------------------
	def set(self, flair):
		self.flair = flair
		self.project = flair.project
		self.fill()

	# ----------------------------------------------------------------------
	# fill run list
	# ----------------------------------------------------------------------
	def fill(self):
		active = self.index(ACTIVE)
		self.delete(0, END)
		for run in self.project.runs:
			# ignore child
			if run.parent != "": continue
			name = run.getInputName2()
			if run.family:
				if name in self.expanded:
					exp = "-"
				else:
					exp = "+"
				self.insert(END, (exp, name, str(len(run.family))))
				self.itemconfig(END, foreground=Project._statusColor[run.status])
				if name in self.selected: self.selection_set(END)

				if exp == "-":
					for child in run.family:
						childRun = self.project.getRunByName(child, run.name)
						self.insert(END, (_CHILD, "  %s"%(child),""))
						self.itemconfig(END,
							foreground=Project._statusColor[childRun.status])
						if child in self.selected: self.selection_set(END)
			else:
				self.insert(END, ("",name,""))
				self.itemconfig(END, foreground=Project._statusColor[run.status])
				if name in self.selected: self.selection_set(END)

		if not self.selected:
			self.selection_set(0)

		self.activate(active)
		self.see(active)
		#self.selected =[]	# To avoid updating the fields the first time
		#self.select()
		#self.focus_set()
		#self.update_idletasks()

	# ----------------------------------------------------------------------
	# Update run colors
	# ----------------------------------------------------------------------
	def updateColors(self):
		for i in range(self.size()):
			run = self.getRun(i)
			if run:
				self.itemconfig(i, foreground=Project._statusColor[run.status])

	# ----------------------------------------------------------------------
	# Expand selected items
	# ----------------------------------------------------------------------
	def expand(self):
		sel = list(map(int, self.curselection()))
		if len(sel)==0: return

		sel.reverse()
		for i in sel:
			if self.lists[0].get(i) != "+": continue
			self.expanded.append(self.lists[1].get(i))

		self.updateSelection()
		self.fill()
		self.activate(sel[-1])

	# ----------------------------------------------------------------------
	# Collapse selected items
	# ----------------------------------------------------------------------
	def collapse(self):
		sel = list(map(int, self.curselection()))
		if len(sel)==0: return
		sel.reverse()

		for i in sel:
			if self.lists[0].get(i) != "-": continue
			try: self.expanded.remove(self.lists[1].get(i))
			except ValueError: pass

		self.updateSelection()
		self.fill()
		self.activate(sel[-1])

	# ----------------------------------------------------------------------
	def toggle(self, event=None):
		sel = list(map(int, self.curselection()))
		if len(sel)==0: return
		sel.reverse()

		for i in sel:
			item = self.lists[0].get(i)
			if item not in ("+","-"): continue
			name = self.lists[1].get(i)

			if item == "+":
				self.expanded.append(name)
			else:
				try: self.expanded.remove(name)
				except ValueError: pass

		self.updateSelection()
		self.fill()
		self.activate(sel[-1])

	# ----------------------------------------------------------------------
	def buttonRelease(self, event):
		sel = list(map(int, self.curselection()))
		listbox = event.widget
		act = listbox.nearest(event.y)
		self.activate(act)
		lid = self.lists.index(listbox)

		if lid==0:
			if listbox.get(act) in ("+","-"): self.toggle()
			return

		elif lid==2:
			# Change the number of spawned jobs
			if self.lists[0].get(act) == _CHILD: return
			self.spawnInPlace()

	# ----------------------------------------------------------------------
	# Update/remember selection
	# ----------------------------------------------------------------------
	def updateSelection(self):
		del self.selected[:]
		for i in self.curselection():
			self.selected.append(self.lists[1].get(i).strip())

	# ----------------------------------------------------------------------
	def getName(self, pos):
		return self.lists[1].get(pos).strip()

	# ----------------------------------------------------------------------
	# Return run from position
	# ----------------------------------------------------------------------
	def getRun(self, pos):
		idx = int(self.index(pos))
		if idx==0: return self.project.runs[0]
		name = self.lists[1].get(idx).strip()
		if self.lists[0].get(idx) == _CHILD:
			# Child find parent
			for i in range(idx-1,0,-1):
				if self.lists[0].get(i) in ("+","-"):
					parent = self.lists[1].get(i).strip()
					break
			else:
				parent = self.project.runs[0].name
		else:
			parent = None
		return self.project.getRunByName(name, parent)

	# ----------------------------------------------------------------------
	# Return run list selection
	# ----------------------------------------------------------------------
	def runSelection(self):
		runs = []
		for i in self.curselection():
			run = self.getRun(i)
			if run is None: continue
			if run.family and run.name not in self.expanded:
				runs.append(run)
				for child in run.family:
					runs.append(self.project.getRunByName(child, run.name))
			else:
				runs.append(run)
		return runs

	# ----------------------------------------------------------------------
	def spawnButton1(self, event):
		act = self.lists[2].nearest(event.y)
		if self.lists[2].selection_includes(act):
			self.activate(act)
			return "break"

	# ----------------------------------------------------------------------
	# Hitting F2 rename in place
	# ----------------------------------------------------------------------
	def renameInPlace(self, event=None):
		if not self.edit: return

		sel = list(map(int,self.curselection()))
		if len(sel) != 1 or sel[0]==0: return "break"
		pos = sel[0]
		self.activate(pos)

		item = self.lists[0].get(pos)
		if item==_CHILD: return "break"		# Ignore child

		name = self.lists[1].get(pos)
		run = self.project.getRunByName(name.strip())
		if run is None or run.parent != "": return "break"

		# Set the correct name
		self.lists[1].set(pos, name.strip())
		self.lists[1].focus_set()
		edit = tkExtra.InPlaceFile(self.lists[1],
				title="Rename run '%s' to"%(run.name),
				filetypes=[("Fluka Input files","*.inp"),
					("All","*")])

		if edit.value is not None and edit.value != name.strip():
			f = self.project.relativePath(edit.value)
			n,e = os.path.splitext(f)
			self.renameRun(run, n)
			self.fill()
			self.flair.setModified()
		else:
			self.lists[1].set(pos, name)
		return "break"

	# ----------------------------------------------------------------------
	# Rename a run
	# @param run		to rename
	# @param name		new name
	# ----------------------------------------------------------------------
	def renameRun(self, run, name):
		# Check if run is there
		newname = self.project.findUniqueRunName(name)
		if newname != name:
			messagebox.showerror("Run already exists",
				"Run '%s' already exists. Renamed to '%s'"%(name,newname),
				parent=self)
			name = newname

		self.updateSelection()

		try: self.selected.remove(run.name)
		except ValueError: pass
		self.selected.append(name)

		try:
			self.expanded.remove(run.name)
			self.expanded.append(name)
		except ValueError:
			pass

		self.project.renameRun(run,name)

	# ----------------------------------------------------------------------
	# In place spawn runs
	# ----------------------------------------------------------------------
	def spawnInPlace(self, event=None):
		if not self.edit: return

		#self.get()
		# protect against accidental delete
		#if self.flair.event is not None:
		#	if self.focus_get() is not self.runList: return

		sel = list(map(int, self.curselection()))
		if len(sel) == 0: return
		active = self.index(ACTIVE)
		if active not in sel:
			self.activate(sel[0])
		self.lists[2].focus_set()

		try: oldSpawn = int(self.lists[2].get(sel[0]))
		except: oldSpawn = 0
		edit = tkExtra.InPlaceEdit(self.lists[2])

		if edit.value is not None:
			try:
				spawn = int(edit.value)
			except:
				spawn = 0
		else:
			return

		if spawn == oldSpawn: return

		sel.reverse()
		for pos in sel:
			run = self.getRun(pos)
			if run.parent != "": continue
			index = self.project.getRunIndex(run)

			# Fill runs
			for i in range(spawn):
				name = self.project.getSpawnRunName(run,i+1)
				if self.project.getRunByName(name, run.name) is not None:
					continue

				child = run.clone()
				child.name   = name
				child.parent = run.name
				child.index  = i+1
				if name == Project.DEFAULT_INPUT:
					# XXX Copy the default executable
					child.exe = self.project.exe
				self.project.setRunRandomSeed(child, run, True)
				self.project.runs.insert(index+i+1, child)
				run.family.append(child.name)

			# Remove additional runs
			i = spawn+1
			while len(run.family) > spawn:
				name = self.project.getSpawnRunName(run,i)
				if name in run.family:
					child = self.project.getRunByName(name, run.name)
					if child is not None: self.project.delRun(child)
				i += 1
				if i>1000: break

		self.flair.setModified()
		self.updateSelection()
		self.fill()

#===============================================================================
if __name__ == "__main__":
	import sys
	tk = Tk()
	Input.init()
	project = Project.Project(sys.argv[1])
	lst = RunListbox(tk, #borderwidth=0,
			#selectmode=EXTENDED,
			#selectborderwidth=0,
			#background="White",
			#takefocus=True,
			#width=24,
			height=20,
			#exportselection=FALSE
			)
	lst.pack(expand=YES,fill=BOTH)
	lst.show(project)
	tk.mainloop()
