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
# Date:	20-Sep-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import re
import time
import string

from stat import *

try:
	from tkinter import *
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	import tkinter.messagebox as messagebox

from Input import Input

import Manual
import Ribbon
import tkFlair
import tkExtra
import Project
import RunList
import Unicode
import FlairRibbon
import bFileDialog

import DataProcess

from log import say

USRTYPES = [x for x in list(Project._usr_type.keys()) if len(x)>1]
USRTYPES.sort()

#===============================================================================
# Commands Menu
#===============================================================================
class _CommandsMenu(FlairRibbon.FlairMenu):
	#----------------------------------------------------------------------
	def createMenu(self):
		menu = Menu(self, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		for cmd in sorted(Project._usr_type):
			if len(cmd)==1: continue
			try:
				icon = tkFlair.icons[cmd.upper()[:8]]
			except:
				icon = tkFlair.icons["empty"]
			menu.add_command(label=cmd,
				compound=LEFT, image=icon,
				command=lambda t=cmd,p=self.page:p.addUsr(t))

		if not Project.Command.commands: return menu

		menu.add_separator()
#		for cmd in sorted(Project.Command.commands.keys()): #, key=lambda s: s.lower()):
		for cmd in sorted([x for x in Project.Command.commands
						if isinstance(x,str)]):
			if isinstance(cmd,int): continue
			icon = tkFlair.icons["empty"]
			menu.add_command(label=cmd,
				compound=LEFT, image=icon,
				command=lambda t=cmd,p=self.page:p.addUsr(t))

		return menu

#===============================================================================
# Data Tab in Page Frame
#===============================================================================
class DataTab(FlairRibbon.FlairTab):
	#----------------------------------------------------------------------
	def __init__(self, master, page, **kw):
		FlairRibbon.FlairTab.__init__(self, master, page, **kw)

		# --- Run-Usrxxx splitter ---
		splitter = PanedWindow(self, orient=VERTICAL,
			showhandle=0,
			sashrelief=RAISED, sashwidth=4, opaqueresize=1)
		splitter.pack(expand=YES, fill=BOTH)

		# --- Usrxxx ---
		frame = LabelFrame(splitter, text="Usrxxx", foreground="DarkBlue")
		splitter.add(frame, minsize=50)

		self.usrList = tkExtra.MultiListbox(frame,
				(("Run",      15, None),
				 ("Command",   8, None),
				 ("Output",   30, None),
				 ("Unit",      4, None)),
				 background="White",
				 height=10)
		self.usrList.pack(expand=YES, fill=BOTH)
#		self.usrList.sortAssist = None

		self.usrList.bind('<<ListboxSelect>>',     self.fillListboxes)
		self.usrList.bindList('<F2>',              self.usrInPlaceEdit)
		self.usrList.bindList('<Control-Key-e>',   self.usrInPlaceEdit)
		self.usrList.bindList('<Return>',          self.usrInPlaceEdit)
		self.usrList.bindList('<KP_Enter>',        self.usrInPlaceEdit)
		self.usrList.bindList('<ButtonRelease-1>', self.usrInPlaceEdit)
		self.usrList.bindList('<Delete>',          self.deleteUsr)
		self.usrList.lists[3].bind('<Button-1>',   self.button1Ignore)

		self.usrList.setPopupMenu(
				[('Add',     0, self.addUsr,    tkFlair.icons["add"]),
				 ('Remove', -1, self.deleteUsr, tkFlair.icons["x"]),
				 ('Scan',    0, self.filter,    tkFlair.icons["find"]) ])

		# --- Frame with Multilistbox ---
		self.tabPage = tkExtra.TabPageSet(splitter,
				pageNames=[("Parameters", tkFlair.icons["gear"]),
					   ("Files",      tkFlair.icons["files"])])
		splitter.add(self.tabPage)
		self.tabPage.changePage("Files")

		# --- Parameters ---
		frame = self.tabPage["Parameters"]
		self.inputList = tkExtra.MultiListbox(frame,
				(("Name",   10, None),
				 ("Tag",     5, None),
				 ("Value",  32, None),
				 ("Type",   10, None),
				 ("Default",10, None)),
				 background="White")
		self.inputList.pack(expand=YES, fill=BOTH)
		self.inputList.sortAssist = None

		self.inputList.lists[2].bind('<F2>',              self.inputInPlaceEdit)
		self.inputList.lists[2].bind('<Control-Key-e>',   self.inputInPlaceEdit)
		self.inputList.lists[2].bind('<Return>',          self.inputInPlaceEdit)
		self.inputList.lists[2].bind('<KP_Enter>',        self.inputInPlaceEdit)
		self.inputList.lists[2].bind('<ButtonRelease-1>', self.inputInPlaceEdit)

		# --- Files ---
		frame = self.tabPage["Files"]
		self.fileList = tkExtra.ColorMultiListbox(frame,
				(("File", 30, None),
				 ("Type", 10, None),
				 ("Size",  8, None),
				 ("Date", 15, None)),
				 background="White")
		self.fileList.pack(expand=YES, fill=BOTH)

		self.fileList.setPopupMenu(
				[('Add',     0, self.add,     tkFlair.icons["add"]),
				 ('Remove', -1, self.delete,  tkFlair.icons["x"]),
				 ('Filter',  0, self.filter,  tkFlair.icons["filter"]),
				 ('Refresh', 0, self.refresh, tkFlair.icons["refresh"])])
		self.fileList.bindList('<Delete>',        self.delete)

		self.process = DataProcess.DataProcess(self.project)

	# ----------------------------------------------------------------------
	# create the ribbon buttons for the Run tab
	# ----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairTab.createRibbon(self)

		group = Ribbon.LabelGroup(self.ribbon, "Command")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.label["background"] = Ribbon._BACKGROUND_GROUP4

		group.grid3rows()

		# ---
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["scan32"],
				text="Scan",
				compound=TOP,
				command=self.scan,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, rowspan=3, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Scan input for detectors")

		# ---
		col,row = 1,0
		b = _CommandsMenu(group.frame, self,
				text="New " + Unicode.BLACK_DOWN_POINTING_SMALL_TRIANGLE,
				image=tkFlair.icons["add"],
				compound=LEFT,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add new command")

		# ---
		col,row = 1,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Remove",
				compound=LEFT,
				anchor=W,
				command=self.deleteUsr,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Remove command")

		# ---
		col,row = 1,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["rename"],
				text="Rename",
				compound=LEFT,
				anchor=W,
				command=self.rename,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Rename command")

		# ========== Data/Files ===========
		group = Ribbon.LabelGroup(self.ribbon, "Files")
		group.label["background"] = Ribbon._BACKGROUND_GROUP4
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		group.grid3rows()

		# ---
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["refresh32"],
				text="Refresh",
				compound=TOP,
				command=self.refresh,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, rowspan=3, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Refresh file")

		# ---
		col,row = 1,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["add"],
				text="Add",
				compound=LEFT,
				anchor=W,
				command=self.add,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add detector")

		# FIXME SHOULD  we ask for the UNIT?

		# ---
		col,row = 1,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Remove",
				compound=LEFT,
				anchor=W,
				command=self.delete,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Remove detector")

		# ---
		col,row = 1,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["filter"],
				text="Rules",
				compound=LEFT,
				anchor=W,
				command=self.filter,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Open the filter/rules dialog")

		# ======= Data =========
		group = Ribbon.LabelGroup(self.ribbon, "Action")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.label["background"] = Ribbon._BACKGROUND_GROUP4

		group.frame.grid_rowconfigure(0, weight=1)

		# ---
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clean32"],
				text="Clean",
				compound=TOP,
				command=self.clean,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, expand=YES, fill=Y)
		tkExtra.Balloon.set(b, "Remove all generated files")

		# ---
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["compress32"],
				text="Process",
				compound=TOP,
				command=self.execute,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, expand=YES, fill=Y)
		tkExtra.Balloon.set(b, "Merged data files [Ctrl-Enter]")

		return self.ribbon

	# ----------------------------------------------------------------------
	def activate(self):
		self.fillUsrList()
		self.selectAllUsr()
		self.fillListboxes()
		self.refreshStatus()

	# --------------------------------------------------------------------
	# Handle global refresh message
	# --------------------------------------------------------------------
	def refresh(self):
		# Restore selections
		self.fillUsrList()
		self.fillListboxes()

	# --------------------------------------------------------------------
	def scan(self, all=False):
		# XXX FIXME XXX when finding ASCII files issue a warning message
		for r in self.page.runList.curselection():
			run = self.page.runList.getRun(r)
			self.project.findRunUsrInfo(run)
		self.fillUsrList()
		self.selectAllUsr()
		self.fillListboxes()

	# --------------------------------------------------------------------
	# Fill USRxxx list
	# --------------------------------------------------------------------
	def fillUsrList(self, event=None):
		self.log = self.flair.newLog("Data", "Process")
		oldSize = self.log.size()	# remember size before loading
		self.page.busy()
		savesort = self.usrList.saveSort()
		self.usrList.delete(0,END)
		for r in self.page.runList.curselection():
			run = self.page.runList.getRun(r)
			if run is None: continue
			if len(run.usrinfo) == 0:
				self.project.findRunUsrInfo(run)
			runname = run.getInputName2()
			for usr in run.usrinfo:
				usrname = self.project.relativePath(usr.name())
				self.usrList.insert(END, (runname, usr.typeLong(), usrname, usr.unit))
		self.usrList.restoreSort(savesort)
		self.page.notBusy()
		if self.log.size() != oldSize:
			self.flair.notify("Errors during file scanning",
				"Error or warning messages appeared while USRxxx file scanning. " \
				"Open the Output window for more information.", tkFlair.NOTIFY_ERROR,
				"Output", self.flair.showOutput)

	# --------------------------------------------------------------------
	def selectAllUsr(self):
		self.usrList.selection_set(0, END)
		self.usrList.activate(0)
		self.usrList.see(0)
		#self.usrList.update_idletasks()

	# --------------------------------------------------------------------
	# Find UsrInfo for a specific usrname from the selected ones
	# --------------------------------------------------------------------
	def findUsrInfo(self, usrname, runname=None):
		if runname is not None:
			run = self.project.getRunByName(runname)
			for usr in run.usrinfo:
				if usrname == self.project.relativePath(usr.name()):
					return usr
		else:
			for r in self.page.runList.curselection():
				run = self.page.runList.getRun(r)
				for usr in run.usrinfo:
					if usrname == self.project.relativePath(usr.name()):
						return usr
		return None

	# --------------------------------------------------------------------
	# Return selected UsrInfo list
	# --------------------------------------------------------------------
	def getUsr(self):
		usr = []
		for u in self.usrList.curselection():
			runname = self.usrList.lists[0].get(u)
			usrname = self.usrList.lists[2].get(u)
			u = self.findUsrInfo(usrname,runname)
			if u is not None: usr.append(u)
		return usr

	# --------------------------------------------------------------------
	# Return selected UsrInfo
	# --------------------------------------------------------------------
	def selectedUsr(self):
		sel = self.usrList.curselection()
		if len(sel) != 1: return None
		return self.findUsrInfo(self.usrList.lists[2].get(sel[0]), self.usrList.lists[0].get(sel[0]))

	# --------------------------------------------------------------------
	# Set files in list
	# --------------------------------------------------------------------
	def fillListboxes(self, event=None):
		self.page.busy()
		savesort = self.fileList.saveSort()

		# Empty listboxes
		self.inputList.delete(0,END)
		self.fileList.delete(0,END)

		# Check if all usr selected are the same
		# FIXME NEEDS CHECKING FOR MULTIPLE SELECTION
		t = None
		for usr in self.getUsr():
			if t is None: t = usr.type
			if t != usr.type:
				t = None
				break

		# Fill filelist
		self.files = []
		repeat = None
		for usr in self.getUsr():
			self.files.extend(usr.depends())
			if t is not None and usr.cmd is not None:
				for i,(ctag,ctype,cdef,cname) in enumerate(usr.cmd.arg):
					try:
						u = usr.arg[i]
					except:
						u = ""
					if not u:
						if ctype=="output":
							u = usr.name()
						elif cdef:
							u = usr.expand(cdef)
					self.inputList.insert(END,(cname,ctag,u,ctype,cdef))

				i = 0	# rule index
				j = 0	# input index
				while i < len(usr.cmd.inp):
					ctype,cdef,cname = usr.cmd.inp[i]
					i += 1
					try:
						u = usr.inp[j]
						j += 1
					except:
						u = ""
					if not u:
						if ctype=="output":
							u = usr.name()
						elif cdef:
							u = usr.expand(cdef)
					if ctype=="repeat":
						repeat = int(cname)
						if u != cdef:
							# correct pointers and repeat sequence
							i  = repeat-1
							j -= 1
							continue
						else:
							ctype,cdef,cname = usr.cmd.inp[int(cname)-1]
							cname += " [repeat?]"
							ctype = "%s:%d"%(ctype,repeat)
							repeat = None

					if repeat is not None:
						ctype = "%s:%d"%(ctype,i)

					self.inputList.insert(END,(cname,"",u,ctype,cdef))

				t = None

		self.files.sort()
		self.updateFiles()
		self.fileList.restoreSort(savesort)
		self.page.notBusy()

	# --------------------------------------------------------------------
	# Update file list
	# --------------------------------------------------------------------
	def updateFiles(self):
		# Populate list
		self.fileList.delete(0,END)
		totalSize = 0
		for fn in self.files:
			try:
				s  = os.stat(fn)
				sz = s[ST_SIZE]
				ti = s[ST_MTIME]
			except OSError:
				self.log("WARNING File not found: %s"%(fn))
				continue
			totalSize += sz
			ext,color = bFileDialog.fileTypeColor(fn)
			self.fileList.insert(END,
				(fn, ext, sz, time.strftime(bFileDialog._TIME_FORMAT,
					time.localtime(ti))))
			if color is not None:
				self.fileList.itemconfig(END,foreground=color)
		#self.fileList.update_idletasks()
		if self.flair is None: return

	# --------------------------------------------------------------------
	def refreshStatus(self):
		if self.process.isAlive():
			self.page.setStatus("%2d%%  %s: %s" % \
				(self.process.percent,
				 self.process.message[0],
				 self.process.message[1]), "DarkGreen")
		else:
			self.page.setStatus("Files: %d" % (len(self.files)))

	# --------------------------------------------------------------------
	# Add a new file
	# --------------------------------------------------------------------
	def add(self):
		usr = self.selectedUsr()
		if usr is None:
			messagebox.showwarning("Please select one usrxxx",
				"Select only one detector to add a file(s) to it",
				parent=self)
			return

		inp = usr.run.getInputName()
		files = bFileDialog.askopenfilenames(master=self,
				 initialdir=self.project.dir,title='Add Files',
					 filetypes=[
						("Fluka Run fort files",
							("%s*_fort.%d"%(inp,usr.unit),
							 "%s*_ftn.%d"%(inp,usr.unit))),
						("Fluka unit %d fort files" % (usr.unit),
							("*_fort.%d"%(usr.unit),
							 "*_ftn.%d"%(usr.unit))),
						("Fluka fort files","*_fort.*"),
						("Fluka ftn files","*_ftn.*"),
						("All","*")])
		if len(files) > 0:
			for f in files:
				fn = self.project.relativePath(f)
				usr._doRule(self.files, '+', "", fn)
				usr.rules.append(('+',fn))
				self.page.setModified()
			self.updateFiles()

	# --------------------------------------------------------------------
	# Delete a file from list
	# --------------------------------------------------------------------
	def delete(self, event=None):
		usr = self.selectedUsr()
		if usr is None:
			messagebox.showwarning("Please select one usrxxx",
				"Select only one detector to delete a file(s) from it",
				parent=self)
			return
		cursel = list(map(int, self.fileList.curselection()))

		# Delete all?
		if len(cursel) == self.fileList.size():
			usr.rules  = []
			self.files = []
		else:
			# Delete selected files
			for i in cursel:
				fn = self.fileList.get(i)[0]
				# search if + rule already exists then remove it
				for i,(o,f) in enumerate(usr.rules):
					if f==fn:
						if o=="+": del usr.rules[i]
						break
				else:
					usr.rules.append(("-",fn))
				usr._doRule(self.files,'-',"",fn)
		self.updateFiles()

	# --------------------------------------------------------------------
	def filter(self):
		usr = self.getUsr()
		dlg = RulesDialog(self.master)
		rules = dlg.show(usr)
		if rules is not None:
			# Set rules for all runs
			for u in usr:
				if u is not None:
					u.rules = rules[:]
		self.fillListboxes()
		self.page.setModified()

	# --------------------------------------------------------------------
	# Set USRxxx
	# --------------------------------------------------------------------
	def addUsr(self, cmd="usrbin"):
		t = cmd
		if len(cmd)>1:
			try:
				t = Project._usr_type[cmd]
			except KeyError:
				pass

		# do for all runs
		added = []
		for r in self.page.runList.curselection():
			run = self.page.runList.getRun(r)
			# do for all selected usr
			#for u in self.getUsr():
			usr = Project.UsrInfo(run)
			usr.set(type=t)
			if usr.cmd: usr.name(usr.cmd.defaultName())
			run.usrinfo.append(usr)
			added.append((run.getInputName2(), usr.typeLong(), usr.name()))

		# Update display
		self.fillUsrList()
		self.usrList.selection_clear(0,END)

		see = END
		for r,t,n in added:
			for i,row in enumerate(self.usrList.get(0,END)):
				if row[0]==r and row[1]==t and row[2]==n:
					if see == END: see = i
					self.usrList.selection_set(i)
					self.usrList.activate(i)
		self.fillListboxes()
		self.usrList.see(see)
		self.page.setModified()

	# --------------------------------------------------------------------
	def deleteUsr(self, event=None):
		sel = list(map(int, self.usrList.curselection()))
		if len(sel) == 0: return
		if not tkFlair.askyesno("Delete USRxxx files", \
			"You are about to delete %d USRxxx summary data files. " \
			"Do you want to continue?" % (len(sel)), parent=self.ribbon):
			return

		sel.reverse()
		for u in sel:
			usrname = self.usrList.lists[2].get(u)
			usr = self.findUsrInfo(usrname)
			usr.run.usrinfo.remove(usr)
		self.fillUsrList()
		self.page.setModified()
		self.fillListboxes()

	# --------------------------------------------------------------------
	def rename(self, event=None):
		for u in map(int , self.usrList.curselection()):
			usrname = self.usrList.lists[2].get(u)
			usr = self.findUsrInfo(usrname)
			if usr is None: return
			fn = bFileDialog.asksaveasfilename(master=self,
					title="Rename %s"%(usrname),
					initialfile=usrname)
			if len(fn)>0:
				usr._name = self.project.relativePath(fn)
		self.fillUsrList()
		self.fillListboxes()
		self.page.setModified()

	# --------------------------------------------------------------------
	# Button1 hit on selected item, ignore then action so not to unselect
	# --------------------------------------------------------------------
	def button1Ignore(self, event):
		act = event.widget.nearest(event.y)
		if event.widget.selection_includes(act):
			event.widget.focus_set()
			event.widget.activate(act)
			return "break"

	# --------------------------------------------------------------------
	# USR list edit in Place
	# --------------------------------------------------------------------
	def usrInPlaceEdit(self, event=None):
		listbox = event.widget
		if event.type == "2":
			# keyboard event
			active = listbox.index(ACTIVE)
		else:
			active = listbox.nearest(event.y)
			self.usrList.activate(active)
		lid = self.usrList.lists.index(listbox)
		sel = self.usrList.curselection()

		cmd = self.usrList.lists[1].get(active)

		ypos = listbox.yview()[0]	# remember y position

		changed = 0
		if lid == 1:
			lst = USRTYPES[:]
			for cmd in sorted(Project.Command.commands.keys()): #, key=lambda s: s.lower()):
				if isinstance(cmd,int): continue
				lst.append(cmd)
			# Change the usrxxx type
			edit = tkExtra.InPlaceList(listbox, values=lst, height=8)
			if edit.value is None: return

			for usr in self.getUsr():
				if usr.typeLong() != edit.value:
					try:
						usr.type = Project._usr_type[edit.value]
						changed += 1
						self.page.setModified()
					except:
						pass

		elif lid == 2:
			# Edit the filename
			edit = tkExtra.InPlaceFile(listbox, title="Rename usrxxx")
			if edit.value is None: return
			listbox.set(active, edit.value)
			for usr in self.getUsr():
				if usr.typeLong() != cmd: continue
				fn = self.project.relativePath(edit.value, usr.run)
				usr.name(fn)
				changed += 1
				self.page.setModified()

		elif lid == 3:
			edit = tkExtra.InPlaceSpinbox(listbox, from_=21, to_=99)
			if edit.value is None: return
			try: edit.value = int(edit.value)
			except: return

			changed = 2	# Force filling of UsrList below
			for usr in self.getUsr():
				if usr.typeLong() != cmd: continue
				if usr.unit != edit.value:
					usr.unit = edit.value
					row = self.usrList.get(active)
					row[1] = usr.typeLong()
					row[2] = self.project.relativePath(usr.name())
					changed += 1
					self.page.setModified()

		if changed>1: self.fillUsrList()
		for s in sel: self.usrList.selection_set(s)
		self.usrList.selection_set(active)
		self.usrList.activate(active)
		self.usrList.yview_moveto(ypos)
		self.fillListboxes()

	# --------------------------------------------------------------------
	# Input list edit in Place
	# --------------------------------------------------------------------
	def inputInPlaceEdit(self, event=None):
		listbox = event.widget
		if event.type == "2":
			# keyboard event
			active = listbox.index(ACTIVE)
		else:
			active = listbox.nearest(event.y)
			self.inputList.activate(active)
		#lid = self.inputList.lists.index(listbox)

		usr = self.selectedUsr()
		if usr is None or len(usr.type)==1: return
		ypos = listbox.yview()[0]	# remember y position
		oldname = usr.name()

		# depends on type
		name,tag,value,vartype,default = self.inputList.get(ACTIVE)
		if ":" in vartype:
			vartype,repeat = vartype.split(":")
			repeat = int(repeat)-1
		else:
			repeat = None

		if vartype == "int":
			edit = tkExtra.InPlaceInteger(listbox)
		elif vartype == "float":
			edit = tkExtra.InPlaceFloat(listbox)
		elif vartype == "output":
			edit = tkExtra.InPlaceFile(listbox, title=name, save=True,
				filetypes=[("All","*"),
					("Data Files","*.dat"),
					("Fluka fort files","*_fort.*"),
					("Tab lis","*_tab.lis"),
					("Voxel",".vxl"),
					("Usrbin","*.bnn")])
		elif vartype == "file":
			edit = tkExtra.InPlaceFile(listbox, title=name, save=False,
				filetypes=[("All","*"),
					("Data Files","*.dat"),
					("Fluka fort files","*_fort.*"),
					("Tab lis","*_tab.lis"),
					("Voxel",".vxl"),
					("Usrbin","*.bnn")])
		elif vartype == "filelist":
			edit = tkExtra.InPlaceEdit(listbox)
		elif "," in vartype:
			choices = [""]
			choices.extend(vartype.split(","))
			edit = tkExtra.InPlaceList(listbox, values=choices)
		elif vartype == "fixed":
			# Do nothing
			return
		else:
			edit = tkExtra.InPlaceEdit(listbox)

		if edit.value is None: return

		if vartype in ("output","file"): #,"filelist")
			edit.value = self.project.relativePath(edit.value)

		# Argument?
		if active < len(usr.cmd.arg):
			idx = active
			if not edit.value:
				ctag,ctype,cdef,cname = usr.cmd.arg[idx]
				edit.value = usr.expand(cdef)
			usr.setArg(idx, edit.value)

		# or Input?
		else:
			idx = active-len(usr.cmd.arg)
			modified = False
			# Empty
			if not edit.value:
				# if before was not empty and repeat is active
				if value and repeat is not None:
					# delete all repeated terms
					usr.deleteRepeat(active, repeat)
					modified = True

				if repeat is not None:
					ctype,cdef,cname = usr.cmd.inp[repeat]
				else:
					ctype,cdef,cname = usr.cmd.inp[idx]
				edit.value = usr.expand(cdef)

			else:   # Or has a new value
				# if before was empty and repeat is active
				if not value and repeat is not None:
					# insert repeated terms
					usr.insertRepeat(active,repeat)
					modified = True
			usr.setInp(idx, edit.value)

			if modified: self.fillListboxes()

		if vartype == "output":
			usr.name(edit.value)

		listbox.set(active, edit.value)

		if usr.name() != oldname:
			r = usr.run.getInputName2()
			t = usr.typeLong()
			n = usr.name()

			# Name has changed re-fill UsrList and select the usr with the new name
			self.fillUsrList()
			for i,row in enumerate(self.usrList.get(0,END)):
				if row[0]==r and row[1]==t and row[2]==n:
					self.usrList.selection_set(i)
					self.usrList.activate(i)
					self.usrList.see(i)
					break


		self.inputList.selection_set(active)
		self.inputList.activate(active)
		self.inputList.yview_moveto(ypos)

	# --------------------------------------------------------------------
	# Process
	# --------------------------------------------------------------------
	def execute(self):
		if self.process.isAlive():
			if not tkFlair.askyesno("Data process is active",
					"A data merging process is still running.\n"
					"Do you want to kill it?",
					parent=self.ribbon):
				return
			self.process.kill()
			self.process.wait()

		self.log = self.flair.newLog("Data", "Process")
		self.process.setUsr(self.getUsr())
		self.flair.addJob(self.process, self.log)
		self.page.setStatus("Data merging started...")

	# ----------------------------------------------------------------------
	# clean files
	# ----------------------------------------------------------------------
	def clean(self):
		if not tkFlair.askyesno("Warning",
				"Do you want to remove all generated files?",
				parent=self.ribbon):
			return

		if self.process.isAlive():
			if not tkFlair.askyesno("Data process is active",
					"A data merging process is still running.\n"
					"Do you want to kill it?",
					parent=self.ribbon):
				return
			self.process.kill()
			self.process.wait()

		self.log = self.flair.newLog("Data", "Process")
		self.process.setUsr(self.getUsr())
		rc = self.process.clean()
		self.log.write(self.process.output())
		if self.process.message:
			self.flair.notify(self.process.message[0],
				self.process.message[1])
			self.page.setStatus(self.process.message[1])

	# ----------------------------------------------------------------------
	def help(self):
		Manual.show(":Run:Data:")

#=============================================================================
# Rules dialog
#=============================================================================
class RulesDialog(Toplevel):
	def __init__(self, master):
		Toplevel.__init__(self, master)
		self.title("File Selection Rules")

		frame = LabelFrame(self, text="Rules")
		frame.pack(side=TOP, expand=YES, fill=BOTH, padx=5, pady=3)

		# ---
		self.rulesList = tkExtra.ExListbox(frame,
				borderwidth=0,
				selectmode=EXTENDED,
				selectborderwidth=0,
				background="White",
				takefocus=True,
				width=50,
				exportselection=FALSE)
		self.rulesList.grid(row=0, column=0, rowspan=7, sticky=NSEW)
		sb = Scrollbar(frame, orient=VERTICAL, takefocus=False,
				command=self.rulesList.yview)
		sb.grid(row=0, column=1, rowspan=7, sticky=NSEW)
		self.rulesList.config(yscrollcommand=sb.set)

		# --- Buttons ---
		b = Button(frame, image=tkFlair.icons["edit"],
				command=self.edit)
		tkExtra.Balloon.set(b, "Edit Rule")
		b.grid(row=0, column=2, sticky=EW)

		b = Button(frame, image=tkFlair.icons["add"],
				command=self.add)
		tkExtra.Balloon.set(b, "Add new rule")
		b.grid(row=1, column=2, sticky=EW)

		b = Button(frame, image=tkFlair.icons["del"],
				command=self.delete)
		tkExtra.Balloon.set(b, "Delete rule(s)")
		b.grid(row=2, column=2, sticky=EW)

		b = Button(frame, image=tkFlair.icons["clone"],
				command=self.clone)
		tkExtra.Balloon.set(b, "Duplicate rule(s)")
		b.grid(row=3, column=2, sticky=EW)

		b = Button(frame, image=tkFlair.icons["up"],
				command=self.rulesList.moveUp)
		tkExtra.Balloon.set(b, "Move up")
		b.grid(row=4, column=2, sticky=EW)

		b = Button(frame, image=tkFlair.icons["down"],
				command=self.rulesList.moveDown)
		tkExtra.Balloon.set(b, "Move down")
		b.grid(row=5, column=2, sticky=EW)

		m = Label(frame,
			text="Syntax: [+-][filename | ^regexpr$]\n" \
				"Special Characters:\n" \
				" \\I\tInput name\t" \
				" \\U\tUnit number\n" \
				" \\T\tType (usrtrack...)\t" \
				" \\t\tShort type (t,x,...)\n" \
				" \\e\tDefault extension\n" \
				"  .\tAny character\t" \
				"  *\t0 or more char\n" \
				"  +\t1 or more char\t" \
				"  ?\t0 or 1 match of char\n" \
				"  \\d\tDigit\t\t" \
				"  \\D\tNon Digit\n" \
				"http://docs.python.org/library/re.html",
				foreground="DarkBlue", justify=LEFT,
				relief=GROOVE, anchor=W)
		m.grid(row=7, column=0, columnspan=3, sticky=NSEW)

		frame.grid_rowconfigure(6, weight=1)
		frame.grid_columnconfigure(0, weight=1)

		# --------------
		frame = Frame(self)
		frame.pack(side=BOTTOM, fill=X)

		Button(frame, text="Cancel",
				image=tkFlair.icons["x"],
				compound=LEFT,
				command=self.cancel).pack(side=RIGHT)
		Button(frame, text="Ok",
				image=tkFlair.icons["ok"],
				compound=LEFT,
				command=self.ok).pack(side=RIGHT)

		# --------------
		self.rulesList.setPopupMenu(
				[('Edit',   0, self.edit,   tkFlair.icons["edit"]),
				 ('Add',    0, self.add,    tkFlair.icons["add"]),
				 ('Clone',  0, self.clone,  tkFlair.icons["clone"]),
				 ('Remove',-1, self.delete, tkFlair.icons["x"]) ])
		# --------------
		self.rulesList.bind("<Return>",          self.edit)
		self.rulesList.bind("<KP_Enter>",        self.edit)
		self.rulesList.bind("<Double-Button-1>", self.edit)
		self.rulesList.bind("<Insert>",          self.add)
		self.rulesList.bind("<Delete>",          self.delete)
		self.rulesList.bind("<Control-Up>",      self.rulesList.moveUp)
		self.rulesList.bind("<Control-Down>",    self.rulesList.moveDown)
		self.bind("<Escape>", self.cancel)
		self.protocol("WM_DELETE_WINDOW", self.cancel)

	# ----------------------------------------------------------------------
	def show(self, usr):
		self.rules = None
		self.rulesList.delete(0,END)
		rulesstr = {}
		for u in usr:
			if u is None: continue
			rules = u.rules
			for pat,file in rules:
				rule = "%s%s"%(pat,file)
				if rule not in rulesstr:
					rulesstr[rule] = True
					self.rulesList.insert(END, rule)
		self.transient(self.master)
		self.wait_window()
		return self.rules

	# ----------------------------------------------------------------------
	def ok(self, event=None):
		self.rules = []
		for rule in self.rulesList.get(0,END):
			pat  = "+"
			file = rule
			if len(rule)>1:
				pat  = rule[0]
				if pat=="+" or pat=="-":
					file = rule[1:]
				else:
					pat = "+"
			self.rules.append([pat, file])

		self.destroy()

	# ----------------------------------------------------------------------
	def cancel(self, event=None):
		self.destroy()

	# ----------------------------------------------------------------------
	# Edit rule
	# ----------------------------------------------------------------------
	def edit(self, event=None):
		edit = tkExtra.InPlaceEdit(self.rulesList)

	# ----------------------------------------------------------------------
	# Add new rules
	# ----------------------------------------------------------------------
	def add(self, event=None):
		self.rulesList.insert(END,"")
		self.rulesList.selection_clear(0,END)
		self.rulesList.selection_set(END)
		self.rulesList.activate(END)
		self.edit()
		if self.rulesList.get(ACTIVE).strip() == "":
			self.rulesList.delete(ACTIVE)
			self.rulesList.selection_set(END)
			self.rulesList.activate(END)

	# ----------------------------------------------------------------------
	# Delete rules
	# ----------------------------------------------------------------------
	def delete(self, event=None):
		sel = list(self.rulesList.curselection())
		if len(sel) == 0: return
		sel.reverse()
		for i in sel:
			self.rulesList.delete(i)

	# ----------------------------------------------------------------------
	# Clone
	# ----------------------------------------------------------------------
	def clone(self, event=None):
		sel = list(map(int, self.rulesList.curselection()))
		self.rulesList.selection_clear(0, END)
		for i in sel:
			rule = self.rulesList.get(i)
			self.rulesList.insert(END, rule)
			self.rulesList.selection_set(END)
		self.rulesList.activate(END)
