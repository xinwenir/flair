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
# Author: Vasilis.Vlachoudis@cern.ch
# Date:	06-Sep-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import time
import shutil
from stat import *
from log import say

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import Input
import Utils
import Ribbon
import Project
import FlairRibbon
import CompileProcess

import tkFlair
import tkExtra
import bFileDialog

#===============================================================================
# Database scan
#===============================================================================
class _DatabaseMenuButton(FlairRibbon.FlairMenu):
	#----------------------------------------------------------------------
	def createMenu(self):
		user = list(UserRoutinesDialog.scanInput(self.page.project.input))
		user.sort()

		if len(user) == 0:
			self.page.usermvax()
			return

		menu = Menu(self, tearoff=0)
		for u in user:
			menu.add_command(label=u,
				compound=LEFT, image=tkFlair.icons["new"],
				command=lambda f=u+".f",p=self.page : p.addUsermvaxFile(f))
		return menu

#===============================================================================
# Compile Page
#===============================================================================
class CompilePage(FlairRibbon.FlairPage):
	"""Compile a FLUKA executable from user routines."""

	_name_ = "Compile"
	_icon_ = "compile"

	#----------------------------------------------------------------------
	def init(self):
		self._refreshOn = ("exe",)

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Executable ===========
		group = FlairRibbon.FlairMenuGroup(self.ribbon, self, "Executable",
			[("Export Makefile","compile", self.flair.exportMakefile)])
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		# ---
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["save32"],
				text="Save As",
				compound=TOP,
				command=self.saveAs,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, expand=YES, fill=Y)
		tkExtra.Balloon.set(b, "Save executable as")

		# ========== Files ===========
		group = Ribbon.LabelGroup(self.ribbon, "Files")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		# ===
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["add32"],
				command=self.add,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add fortran file or library")

		col,row = 0,2
		menulist = [	("Fortran", "new",   self.addFortran),
				("C/C++",   "new",   self.addC),
				("Object",  "new",   self.addObject),
				("Library", "new",   self.addLibrary) ]
		b = Ribbon.MenuButton(group.frame, menulist,
				text="Add",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add fortran file or library")

		# ===
		col,row = 1,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["database32"],
				command=self.usermvax,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Copy from FLUKA user routines to the current directory/project")

		col,row = 1,2
		b = _DatabaseMenuButton(group.frame, self,
				text="Database",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Copy from FLUKA user routines to the current directory/project")

		# ---
		col,row = 2,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Remove",
				compound=LEFT,
				command=self.delete,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Remove file from the compile list")

		# ---
		col,row = 2,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["up"],
				text="Move Up",
				compound=LEFT,
				command=self.moveUp,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Move up card")

		# ---
		col,row = 2,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["down"],
				text="Move Down",
				compound=LEFT,
				command=self.moveDown,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Move down card")

		# ========== View ===========
		group = Ribbon.LabelGroup(self.ribbon, "View")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		# ---
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["viewer32"],
				text="Viewer",
				compound=TOP,
				command=self.viewer,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, expand=YES, fill=Y)
		tkExtra.Balloon.set(b, "Open in internal viewer")

		# ---
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["edit32"],
				text="Editor",
				compound=TOP,
				command=self.editor,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, expand=YES, fill=Y)
		tkExtra.Balloon.set(b, "Open in external editor")

		# ========== Build ===========
		group = Ribbon.LabelGroup(self.ribbon, "Action")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		group.grid2rows()

		# ===
		col,row = 0,0
		b = Label(group.frame, image=tkFlair.icons["link"],
					text="Link",
					compound=LEFT,
					font=Ribbon._FONT,
					anchor=E,
					background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=SW)
		tkExtra.Balloon.set(b, "Linking program")

		# ===
		col,row = 0,1
		self.lfluka = tkExtra.Combobox(group.frame,
					width=10,
					font=Ribbon._FONT,
					background=Ribbon._BACKGROUND)
		Frame.config(self.lfluka, background=Ribbon._BACKGROUND, padx=0, pady=0)
		self.lfluka.grid(row=row, column=col, padx=0, pady=0, sticky=N+EW)
		tkExtra.Balloon.set(self.lfluka, "Linking program")

		self.lfluka.fill(Project.LINK_PRG)
		self.lfluka.set(Project.LINK_PRG[0])

		# ---
		col,row = 1,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clean32"],
				text="Clean",
				compound=TOP,
				command=self.clean,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Remove all generated files")

		# ---
		col,row = 2,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["build32"],
				text="Build",
				compound=TOP,
				command=self.execute,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Build executable file")

	#----------------------------------------------------------------------
	# Create Project page
	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self)

		# -----
		row = 0
		Label(self.frame, text="Executable:").grid(row=row, column=0, sticky=E)

		self.exe = Label(self.frame, foreground="DarkBlue", background="White",
				relief=SUNKEN, anchor=W)
		tkExtra.Balloon.set(self.exe, "Executable filename")
		self.exe.grid(row=row, column=1, sticky=EW)
		self.exe.bind('<Button-1>', self.saveAs)

		b = Button(self.frame, image=tkFlair.icons["x"], command=self.defaultExecutable, pady=1, padx=1)
		b.grid(row=row, column=2)
		tkExtra.Balloon.set(b, "Remove executable")

		# --- Developer options ---
		if Input._developer:
			row += 1
			frame = Frame(self.frame)
			frame.grid(row=row, column=0, columnspan=3, sticky=NSEW)

			self.defmain = IntVar()
			self.defmain.set(0)
			c = Checkbutton(frame, text="Default main:",
				variable=self.defmain, command=self.mainEdit)
			c.pack(side=LEFT)
			tkExtra.Balloon.set(c, "Default or override the main program of FLUKA")

			self.main = Entry(frame, background="White", width=15)
			self.main.pack(side=LEFT)
			tkExtra.Balloon.set(self.main, "Object to use as main program")

			Label(frame, text=" Options:").pack(side=LEFT)
			self.opt = Entry(frame, background="White")
			self.opt.pack(side=LEFT, fill=X, expand=TRUE)
			tkExtra.Balloon.set(self.opt, "Options to be passed to the compiler")

			self.dline = IntVar()
			self.dline.set(0)
			c = Checkbutton(frame, text="D Line", variable=self.dline)
			c.pack(side=LEFT)
			tkExtra.Balloon.set(c, "Enable compilation of D lines")

			self.bound = IntVar()
			self.bound.set(1)
			c = Checkbutton(frame, text="Bound Check", variable=self.bound)
			c.pack(side=LEFT)
			tkExtra.Balloon.set(c, "Enable boundary checking of arrays")

		# --- fileList
		row += 1
		self.fileList = tkExtra.ColorMultiListbox(self.frame,
				(('File', 10, None),
				 ('Type',  8, None),
				 ('Size',  6, None),
				 ('Date', 15, None)),
				 height=7)
		self.fileList.grid(row=row,
				column=0, columnspan=3, sticky=NSEW)

		# ------
		self.frame.grid_rowconfigure(row, weight=1)
		self.frame.grid_columnconfigure(1, weight=1)

		self.fileList.setPopupMenu([
				("Edit",   0, self.editor, tkFlair.icons["edit"]),
				("View",   0, self.viewer, tkFlair.icons["view"]),
				None,
				("Add",	   0, self.add,    tkFlair.icons["add"]),
				("Remove", 0, self.delete, tkFlair.icons["x"]),
				("Rename", 0, self.rename, tkFlair.icons["edit"]) ])
		self.fileList.bindList('<Double-1>',        self.editor)
		self.fileList.bindList('<F2>',              self.rename)
		self.fileList.bindList('<Return>',          self.editor)
		self.fileList.bindList('<KP_Enter>',        self.editor)
		self.fileList.bindList('<Control-Key-Up>',  self.moveUp)
		self.fileList.bindList('<Control-Key-Down>',self.moveDown)
		self.fileList.focus_set()

		self.process = CompileProcess.CompileProcess(self.project)

	# ----------------------------------------------------------------------
	# Rename in place filename
	# ----------------------------------------------------------------------
	def rename(self, event = None):
		listbox  = self.fileList.lists[0]
		selected = list(map(int,listbox.curselection()))
		if len(selected) != 1: return
		pos = selected[0]
		old = listbox.get(pos)
		edit = tkExtra.InPlaceEdit(listbox)
		if edit.value and os.path.isfile(edit.value):
			# Replace file
			fn = self.project.relativePath(edit.value)
			self.project.sourceList[pos] = fn
			ext,color = bFileDialog.fileTypeColor(fn)
			s = os.stat(edit.value)
			self.fileList.set(pos,[fn, ext, s[ST_SIZE],
				time.strftime(bFileDialog._TIME_FORMAT, time.localtime(s[ST_MTIME]))])
			if color is not None:
				self.fileList.itemconfig(pos,foreground=color)
			self.flair.notify("Success: File Changed",
				"Success: file %s successfully added"%(edit.value) ,
				tkFlair.NOTIFY_NORMAL)
			self.setModified(True)
		else:
			self.refresh()
			self.flair.notify("Error: File not found",
					"Error: file %s NOT found"%(edit.value) ,
					tkFlair.NOTIFY_ERROR)
		self.fileList.selection_set(pos)
		self.fileList.activate(pos)

	# ----------------------------------------------------------------------
	# Show Frame
	# ----------------------------------------------------------------------
	def refresh(self):
		FlairRibbon.FlairPage.refresh(self)
		if Input._developer:
			self.opt.delete(0, END)
			self.opt.insert(0, self.project.f77opts)
			self.bound.set(int(self.project.f77bound))
			self.dline.set(int(self.project.f77dline))
			self.defmain.set(self.project.defmain)
			self.main["state"] = NORMAL
			self.main.delete(0, END)
			self.main.insert(0, self.project.main)
			self.mainEdit()

		self.lfluka.set(self.project.lfluka)
		self.exe.config(text=Project.abbrevPath(self.project.exe))

		self.fileList.delete(0,END)
		for f in self.project.sourceList:
			try:
				s = os.stat(f)
				ext,color = bFileDialog.fileTypeColor(f)
				self.fileList.insert(END, (f, ext, s[ST_SIZE],
					time.strftime(bFileDialog._TIME_FORMAT,
						time.localtime(s[ST_MTIME]))))
				if color is not None:
					self.fileList.itemconfig(END,foreground=color)
			except OSError:
				self.fileList.insert(END, (f, '-', '-', '-'))

	# --------------------------------------------------------------------
	def refreshStatus(self):
		if self.process.isAlive():
			self.setStatus("%2d%%  %s: %s" % \
				(self.process.percent,
				 self.process.message[0],
				 self.process.message[1]))
		else:
			self.setStatus("Files: %d" % (self.fileList.size()))

	# ----------------------------------------------------------------------
	# Update fields
	# ----------------------------------------------------------------------
	def updateProject(self, event=None):
		if self.page is None or self.project is None: return

		if Input._developer:
			o = self.opt.get()
			if o!=self.project.f77opts:
				self.project.f77opts = o
				self.setModified(True)

			if self.defmain.get() != self.project.defmain:
				self.project.defmain = self.defmain.get()
				self.setModified(True)

			m = self.main.get()
			if m!=self.project.main:
				self.project.main = m
				self.setModified(True)

			d = bool(self.dline.get())
			if d!=self.project.f77dline:
				self.project.f77dline = d
				self.setModified(True)

			b = bool(self.bound.get())
			if b!=self.project.f77bound:
				self.project.f77bound = b
				self.setModified(True)

		l = self.lfluka.get()
		if l!=self.project.lfluka:
			self.project.lfluka = l
			self.setModified(True)

	# ----------------------------------------------------------------------
	def addFortran(self):
		self.add(filetypes=[("Fortran files", ("*.f","*.for","*.F","*.FOR")),
				   ("All","*")])

	def addC(self):
		self.add(filetypes=[("C/C++ files", ("*.c", "*.cc", "*.cpp","*.C","*.CC","*.CPP")),
				   ("All","*")])
	def addObject(self):
		self.add(filetypes=[("Object Files","*.o"),
				   ("All","*")])
	def addLibrary(self):
		self.add(filetypes=[("Libraries",("*.a","*.so")),
				   ("All","*")])

	# ----------------------------------------------------------------------
	# Add a new file
	# ----------------------------------------------------------------------
	def add(self, event=None, filetypes=None):
		if filetypes is None:
			filetypes=[("Fortran files", ("*.f","*.for","*.F","*.FOR")),
				   ("C/C++ files", ("*.c", "*.cc", "*.cpp","*.C","*.CC","*.CPP")),
				   ("Object Files","*.o"),
				   ("Libraries",("*.a","*.so")),
				   ("All","*")]

		files = bFileDialog.askopenfilenames(master=self.page,
				filetypes=filetypes)
		if files is None: return
		for f in files: self.addFile(f)

	# ----------------------------------------------------------------------
	def addFile(self, filename):
		fn = self.project.relativePath(filename)
		ext,color = bFileDialog.fileTypeColor(fn)
		s = os.stat(filename)
		self.fileList.insert(END,[fn, ext, s[ST_SIZE],
			time.strftime(bFileDialog._TIME_FORMAT, time.localtime(s[ST_MTIME]))])
		if color is not None:
			self.fileList.itemconfig(END,foreground=color)
		self.project.sourceList.append(fn)
		self.setModified(True)

	# ----------------------------------------------------------------------
	def addUsermvaxFile(self, filename):
		# Copy the file
		src = os.path.join(Project.flukaDir, "usermvax", filename)
		dst = os.path.join(self.project.dir, filename)
		try:
			os.stat(dst)
			if not tkFlair.askyesno(
				"Override existing file %s"%(filename),
				"File %s already exists in project directory. " \
				"Do you want to copy and override the existing file?" %\
				(filename), parent=self.page):
					return

			# FIXME XXX
			try: self.project.sourceList.remove(filename)
			except: pass
			for i in range(self.fileList.size()):
				fn = self.fileList.lists[0].get(i)
				if fn==filename:
					self.fileList.delete(i)
					break
		except OSError:
			pass

		try:
			shutil.copyfile(src, dst)
			# Add file to the compile list of files
			self.addFile(dst)
		except IOError:
			self.flair.notify("Error in copying file %s"%(filename),
				"Error during copying file %s"%(filename),
				tkFlair.NOTIFY_ERROR)

	# ----------------------------------------------------------------------
	# Delete region
	# ----------------------------------------------------------------------
	def delete(self, event=None):
		# protect against accidental delete
		lst = list(map(int, self.fileList.curselection()))
		lst.reverse()
		src = self.project.sourceList
		for i in lst:
			self.fileList.delete(i)
			del src[i]

	# ----------------------------------------------------------------------
	# Move up select items by one
	# ----------------------------------------------------------------------
	def moveUp(self, event=None):
		src = self.project.sourceList
		for i in map(int,self.fileList.curselection()):
			if i==0: continue
			prev = i-1
			if not self.fileList.selection_includes(prev):
				act = self.fileList.index(ACTIVE)
				self.swap(prev,i)
				self.fileList.selection_set(prev)
				if act == i: self.fileList.activate(prev)
		return "break"

	# ----------------------------------------------------------------------
	# Move down select items by one
	# ----------------------------------------------------------------------
	def moveDown(self, event=None):
		sz  = self.fileList.size()-1
		lst = list(map(int,self.fileList.curselection()))
		lst.reverse()
		for i in lst:
			if i >= sz: continue
			inext = i+1
			if not self.fileList.selection_includes(inext):
				act = self.fileList.index(ACTIVE)
				self.swap(i,inext)
				self.fileList.selection_set(inext)
				if act == i: self.fileList.activate(inext)
		return "break"

	# ----------------------------------------------------------------------
	# Swap two items in the list
	# ----------------------------------------------------------------------
	def swap(self, a, b):
		lb = self.fileList
		if a>b: a, b = b, a

		at = lb.get(a);		al = self.project.sourceList[a]
		bt = lb.get(b);		bl = self.project.sourceList[b]

		lb.delete(b)
		lb.delete(a)

		lb.insert(a, bt);	self.project.sourceList[a] = bl
		lb.insert(b, at);	self.project.sourceList[b] = al

	# ----------------------------------------------------------------------
	# Edit file
	# ----------------------------------------------------------------------
	def editor(self, event=None):
		tkFlair.editor([self.fileList.get(int(x))[0]
				for x in self.fileList.curselection() ])

	# ----------------------------------------------------------------------
	# View file
	# ----------------------------------------------------------------------
	def viewer(self, event=None):
		self.flair.viewer([self.fileList.get(int(x))[0]
				for x in self.fileList.curselection() ])

	# ----------------------------------------------------------------------
	def mainEdit(self):
		if self.defmain.get():
			self.main["state"] = DISABLED
		else:
			self.main["state"] = NORMAL

	# ----------------------------------------------------------------------
	def usermvax(self, event=None):
		UserRoutinesDialog(self.page, self, self.project)

	# ----------------------------------------------------------------------
	# Save executable as
	# ----------------------------------------------------------------------
	def saveAs(self, event=None):
		self.activate()
		fn = bFileDialog.asksaveasfilename(master=self.page,
			title="Save executable as",
			filetypes=[ ("All","*"), ("Exe","*.exe") ])
		if len(fn) > 0:
			self.project.setExecutable(fn)
			self.flair.refresh("exe")
			self.setModified()
			self.exe.config(text=self.project.exe)

	# ----------------------------------------------------------------------
	# Default executable
	# ----------------------------------------------------------------------
	def defaultExecutable(self):
		self.project.setExecutable("")
		self.exe.config(text="")
		self.flair.refresh("exe")
		self.setModified()

	# ----------------------------------------------------------------------
	# Select All
	# ----------------------------------------------------------------------
	def selectAll(self):
		self.fileList.selection_set(0,END)

	# ----------------------------------------------------------------------
	# build executable
	# ----------------------------------------------------------------------
	def execute(self):
		self.updateProject()

		if self.project.exe == "":
			self.saveAs()
			if self.project.exe == "": return

		if self.process.isAlive():
			if not tkFlair.askyesno("Compile process active",
					"A compile process is still running.\n"
					"Do you want to kill it?",
					parent=self.page):
				return
			self.process.kill()
			self.process.wait()

		log = self.flair.newLog(self.name,self.project.exe)

		# Check if libflukahp is newer than the executable
		# then rebuild everything
		lib = os.path.join(Project.flukaDir,Project.DEFAULT_LIB)
		if Utils.isNewer(lib, self.project.exe):
			say("Executable is older than %s rebuilding\n"%(Project.DEFAULT_LIB))
			self.process.clean()

		self.setStatus("Compiling executable")
		self.flair.addJob(self.process, log)

	# ----------------------------------------------------------------------
	# clean files
	# ----------------------------------------------------------------------
	def clean(self):
		if self.process.isAlive():
			if tkFlair.askyesno("Compile process active",
					"A compile process is still running.\n"
					"Do you want to kill it?",
					parent=self.page):
				return
			self.process.kill()
			self.process.wait()

		log = self.flair.newLog(self.name, self.project.exe)
		self.setStatus("Cleaning files")
		rc = self.process.clean()
		log.write(self.process.output())
		if self.process.message:
			self.flair.notify(self.process.message[0],
				self.process.message[1])

#===============================================================================
# User Routines Dialog
#===============================================================================
class UserRoutinesDialog(Toplevel):
	def __init__(self, master, page, project):
		Toplevel.__init__(self, master)
		self.page    = page
		self.project = project
		self.path = os.path.join(Project.flukaDir, "usermvax")
		self.transient(master)

		self.fileList = tkExtra.MultiListbox(self,
				(('File', 10, None),
				 ('Size',  6, None),
				 ('Date', 15, None),
				 ('Desc', 50, None)))
		self.fileList.config(height=20)
		self.fileList.bindList('<Double-1>', self.view)
		self.fileList.pack(side=TOP, expand=YES, fill=BOTH)

		frame = Frame(self)
		frame.pack(side=BOTTOM, fill=X)

		btn = Button(frame, text="Copy to Project", underline=0,
				command=self.copy)
		btn.pack(side=LEFT)

		btn = Button(frame, text="Scan Input", underline=0,
				command=self.scan)
		btn.pack(side=LEFT)

		btn = Button(frame, text="View", underline=0,
				command=self.view)
		btn.pack(side=LEFT)

		btn = Button(frame, text="Close", underline=0,
			command=self.close)
		btn.pack(side=LEFT)

		self.bind("<Escape>", self.close)
		self.protocol("WM_DELETE_WINDOW", self.close)
		self.title("FLUKA User routines")

		self.fill()
		self.scan()

		self.wait_visibility()
		self.grab_set()
		self.focus_set()
		self.fileList.focus_set()
		self.wait_window()

	# ----------------------------------------------------------------------
	def close(self, event=None):
		self.destroy()

	# ----------------------------------------------------------------------
	# Fill with files generated by fluka
	# ----------------------------------------------------------------------
	def fill(self):
		try:
			for filename in os.listdir(self.path):
				fn,ext = os.path.splitext(filename)
				if ext != ".f": continue
				try:
					s = os.lstat(os.path.join(self.path,filename))
					desc = Input._usermvax.get(fn, "")
					self.fileList.insert(END,
						(filename, s[ST_SIZE],
						 time.strftime(bFileDialog._TIME_FORMAT,
							time.localtime(s[ST_MTIME])),
						 desc))
				except OSError:
					say("ERROR accessing file", filename)
		except OSError:
			say("ERROR accessing folder", self.path)
		self.fileList.sort(0)

	# ----------------------------------------------------------------------
	# View document
	# ----------------------------------------------------------------------
	def view(self, event=None):
		self.page.flair.viewer([
			os.path.join(self.path,self.fileList.get(int(x))[0])
				for x in self.fileList.curselection() ])

	# ----------------------------------------------------------------------
	# Copy selected files to project
	# ----------------------------------------------------------------------
	def copy(self, event=None):
		for sel in self.fileList.curselection():
			filename, size, date, desc = self.fileList.get(sel)
			self.page.addUsermvaxFile(filename)

	# ----------------------------------------------------------------------
	# Scan input and highlight necessary files
	# ----------------------------------------------------------------------
	def scan(self):
		user = UserRoutinesDialog.scanInput(self.project.input)

		# Return a set of routines
		active = self.fileList.index(ACTIVE)
		self.fileList.selection_clear(0,END)
		activate = -1
		for i in range(self.fileList.size()):
			fn,ext = os.path.splitext(self.fileList.get(i)[0])
			if fn in user:
				self.fileList.selection_set(i)
				if active == i: activate = i
				if activate < 0:
					activate = i
		if activate>=0:
			self.fileList.activate(activate)
			self.fileList.see(activate)

	# ----------------------------------------------------------------------
	@staticmethod
	def scanInput(input):
		user = set()

		# Check for abscff, dffcff, frghns
		for card in input["OPT-PROP"]:
			sdum = card.sdum().strip()
			user.add("frghns")
			user.add("rflctv")
			if sdum=="":
				if card.intWhat(1)<-99: user.add("rfrndx")
				if card.intWhat(2)<-99: user.add("abscff")
				if card.intWhat(3)<-99: user.add("dffcff")
			elif sdum=="SPEC-BDX":
				user.add("ophbdx")
			elif sdum=="SENSITIV":
				if card.intWhat(1)<-99: user.add("queffc")

		if "OPT-PROD" in input.cards:
			user.add("frghns")

		# Check for comscw, fluscw, endscp, fldscp, usrrnc
		for card in input["USERWEIG"]:
			w3 = card.intWhat(3)
			if w3 > 0: user.add("fluscw")
			if w3==2 or w3==4: user.add("fldscp")
			if card.intWhat(5)>0: user.add("usrrnc")
			w6 = card.intWhat(6)
			if w6 > 0: user.add("comscw")
			if w6==2 or w6==4: user.add("endscp")

		# Check for formfu
		if "MULSOPT" in input.cards:
			user.add("formfu")

		# FUSRBV, LUSRBL, MUSRBR
		for card in input["USRBIN"]:
			if card.intWhat(1) in (8,18):
				user.add("fusrbv")
				user.add("lusrbl")
				user.add("musrbr")
				break

		# LATTIC
		for card in input["LATTICE"]:
			if card.sdum().strip()=="":
				user.add("lattic")
				break

		# MAGFLD
		for card in input["ASSIGNMA"]:
			if card.intWhat(5) == 1:
				user.add("magfld")
				break
		if "MGNFIELD" in input.cards:
			user.add("magfld")

		# MGDRAW
		if "USERDUMP" in input.cards:
			user.add("mgdraw")

		# Source
		if "SOURCE" in input.cards:
			user.add("source")

		# Biasing
		cards = input.cards.get("BIASING")
		if cards is not None:
			user.add("ubsset")
			for card in cards:
				if card.sdum()=="USER":
					user.add("usimbs")
					break
		for tag in ["EMF-BIAS", "EMFCUT", "EXPTRANS",
				"EXPTRANS", "LOW-BIAS", "LOW-DOWN", "WW-FACTO"]:
			if tag in input.cards:
				user.add("ubsset")
				break

		# UDCDRL
		for card in input["LAM-BIAS"]:
			if card.sdum() == "DCY-DIRE" and card.intWhat(4)==0:
				user.add("udcdrl")
				break
			elif card.sdum() == "N1HSCBS":
				user.add("uesdrl")
				break

		# USRINI
		if "USRICALL" in input.cards:
			user.add("usrini")

		# USROUT
		if "USROCALL" in input.cards:
			user.add("usrout")

		# USRGLO
		if "USRGCALL" in input.cards:
			user.add("usrglo")

		# USRMED
		for card in input["MAT-PROP"]:
			if card.sdum()[:8] == "USERDIRE":
				if card.intWhat(1)>0:
					user.add("usrmed")

		return user
