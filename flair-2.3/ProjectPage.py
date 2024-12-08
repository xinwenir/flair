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
# Date:	13-Jun-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import re
import sys
from stat import *
from log import say

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkinter as tk
import Input
import RichText
import tkDialogs
import bFileDialog

import tkFlair
import tkExtra
import Unicode

import Manual
import Ribbon
import Project
import FlairRibbon

try:
	import PIL.Image
	import PIL.ImageTk
	have_PIL = True
except:
	have_PIL = False

#===============================================================================
# Recent Menu button
#===============================================================================
class _RecentMenuButton(FlairRibbon.FlairMenu):
	#----------------------------------------------------------------------
	def createMenu(self):
		rf = [	self.page.loadRecent0,
			self.page.loadRecent1,
			self.page.loadRecent2,
			self.page.loadRecent3,
			self.page.loadRecent4,
			self.page.loadRecent5,
			self.page.loadRecent6,
			self.page.loadRecent7,
			self.page.loadRecent8,
			self.page.loadRecent9 ]
		menu = Menu(self, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		empty = True
		for i in range(tkFlair._maxRecent):
			try:
				filename = tkFlair.config.get(tkFlair.__name__,
						"recent.%d"%(i))
				fn = os.path.basename(filename)
				menu.add_command(label="%d %s"%(i+1, fn),
					compound=LEFT, image=tkFlair.icons["project"],
					command=rf[i])
				empty = False
			except:
				pass
		if empty:
			self.page.flair.openProject()
			return None
		else:
			return menu

#===============================================================================
# Project Tab Page
#===============================================================================
class ProjectPage(FlairRibbon.FlairPage):
	"""Create and/or edit a flair project."""

	_name_ = "Flair"
	_icon_ = "project"

	#----------------------------------------------------------------------
	def init(self):
		self._refreshOn = ("exe","image")

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Project ===========
		group = Ribbon.LabelGroup(self.ribbon, "Project")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row=0,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["newflair32"],
				#command=self.flair.newProject,
				command=self.openTemplate,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "New project with the basic template as input")

		col,row=0,2
		b = FlairRibbon._TemplateMenuButton(group.frame, self,
				text="New",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "New project + input from template")

		# ---
		col,row=1,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["openflair32"],
				command=self.flair.openProject,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Open project")

		col,row=1,2
		b = _RecentMenuButton(group.frame, self,
				text="Open",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Open recent project")

		# ---
		col,row=2,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["saveflair32"],
				command=self.flair.saveProject,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Save project")

		col,row=2,2
		b = Ribbon.LabelButton(group.frame,
				text="Save",
				image=tkFlair.icons["triangle_down"],
				compound=RIGHT,
				command=self.flair.saveProjectAs,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Save project as")

		# ========== Edit ==============
		group = FlairRibbon.FlairMenuGroup(self.ribbon, self, "Edit",
			[("Edit style", "edit", self.editStyle)])
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		col,row = 0,0
		self.styleCombo = Ribbon.LabelCombobox(group.frame,
					width=12,
					command=self.setStyleFromCombo)
		self.styleCombo.grid(row=row, column=col, columnspan=5, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(self.styleCombo, "Set editing style")
		self.styleCombo.fill(self.notes.styleList())

		# ---
		col,row = 0,1
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["alignleft"],
				variable=self._style,
				value="Left")
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Align left text")

		# ---
		col,row = 1,1
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["aligncenter"],
				variable=self._style,
				value="Center")
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Align center text")

		# ---
		col,row = 2,1
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["alignright"],
				variable=self._style,
				value="Right")
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Align right text")

		# ---
		col,row = 3,1
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["hyperlink"],
				variable=self._style,
				value="Link")
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Hyperlink")

		# ---
		col,row = 4,1
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["image"],
				command=self.insertImage,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Insert image")

		# ---
		col,row = 0,2
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["bold"],
				variable=self._style,
				value="Bold")
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Bold font [Ctrl-B]")

		# ---
		col,row = 1,2
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["italic"],
				variable=self._style,
				value="Italic")
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Italics font [Ctrl-I]")

		# ---
		col,row = 2,2
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["underline"],
				variable=self._style,
				value="Underline")
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Underline font [Ctrl-U]")

		# ---
		col,row = 3,2
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["overstrike"],
				variable=self._style,
				value="Overstrike")
		b.grid(row=row,  column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Overstrike font [Ctrl-O]")

		# ========== Publish ==============
		group = Ribbon.LabelGroup(self.ribbon, "Publish")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row=0,0
		b = Ribbon.LabelButton(group.frame, self.page, "<<Document>>",
				text="Document",
				image=tkFlair.icons["new32"],
				#anchor=S,
				compound=TOP,
				state=DISABLED,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NS)
		tkExtra.Balloon.set(b, "Create document from project")

		# ---
		col,row=1,0
		b = Ribbon.LabelButton(group.frame,
				text="Print",
				command=self.hardcopy,
				image=tkFlair.icons["print32"],
				compound=TOP,
				#anchor=S,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NS)
		tkExtra.Balloon.set(b, "Print input")

		# ---
		col,row=2,0
		b = Ribbon.LabelButton(group.frame,
				text="Refresh",
				command=self.refreshButton,
				image=tkFlair.icons["refresh32"],
				compound=TOP,
				#anchor=S,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NS)
		tkExtra.Balloon.set(b, "Refresh document")

		# ========== Tools ==============
		group = Ribbon.LabelGroup(self.ribbon, "Tools")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row=0,0
		b = Ribbon.LabelButton(group.frame, #self.page, "<<Config>>",
				text="Config",
				image=tkFlair.icons["config32"],
				command=self.flair.preferences,
				compound=TOP,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NS)
		tkExtra.Balloon.set(b, "Open configuration dialog")

		# ===
		col,row=1,0
		b = Ribbon.LabelButton(group.frame,
				text="Report",
				image=tkFlair.icons["debug"],
				compound=LEFT,
				command=tkFlair.ReportDialog.sendErrorReport,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Send Error Report")

		# ---
		col,row=1,1
		b = Ribbon.LabelButton(group.frame,
				text="Updates",
				image=tkFlair.icons["GLOBAL"],
				compound=LEFT,
				command=self.flair.checkUpdates,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Check Updates")

		col,row=1,2
		b = Ribbon.LabelButton(group.frame,
				text="About",
				image=tkFlair.icons["about"],
				compound=LEFT,
				command=lambda s=self: tkFlair.aboutDialog(s.page),
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "About flair")

		# ========== Tools ==============
		group = Ribbon.LabelGroup(self.ribbon, "Close")
		group.pack(side=RIGHT, fill=Y, padx=0, pady=0)

		# ---
		b = Ribbon.LabelButton(group.frame,
				text="Exit",
				image=tkFlair.icons["exit32"],
				compound=TOP,
				command=self.flair.quit,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y, padx=0, pady=0)
		tkExtra.Balloon.set(b, "Close program")

	#----------------------------------------------------------------------
	# Create Project page
	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self, False)

		lf = LabelFrame(self.frame, text="Notes", padx=2, pady=2)
		lf.pack(side=BOTTOM, fill=BOTH, expand=YES)

		Label(self.frame, text="Title:").pack(side=LEFT)
		self.title = Entry(self.frame, background="White")
		self.title.pack(side=RIGHT, fill=X, expand=YES)
		tkExtra.Balloon.set(self.title, "Project title")

		self._style = StringVar()
		self._setVar = False
		self._style.set("Normal")
		self._style.trace("w",self.styleChanged)

		self.notes = RichText.RichText(lf,
		#self.notes = tkTextEditor.TextEditor(lf,
			height=4, wrap=WORD, background="white", undo=True)
		self.notes.grid(row=0, column=0, sticky=NSEW)
		self.notes.bind("<ButtonRelease-1>",   self.checkStyle)
		self.notes.bind("<Key>", self.checkStyle)
		self.notes.bind('<Control-Key-f>',self.findDialog)
		self.notes.bind('<Control-Key-r>',self.replaceDialog)
		self.notes.bind('<Control-Key-b>',lambda e,s=self,t="Bold": s.setStyle(t))
		self.notes.bind('<Control-Key-i>',lambda e,s=self,t="Italic": s.setStyle(t))
		self.notes.bind('<Control-Key-u>',lambda e,s=self,t="Underline": s.setStyle(t))
		self.notes.bind('<Control-Key-o>',lambda e,s=self,t="Overstrike": s.setStyle(t))

		self.notes.bind("<<Image>>",     self.openImage)
		self.notes.bind("<<Link>>",      self.openLink)
		self.notes.bind("<<Selection>>", self.selection)
		self.notes_First = self.notes_Last  = None

		sbv = tkExtra.AutoScrollbar(lf, takefocus=False)
		sbv.grid(row=0, column=1, sticky='ns')
		sbh = tkExtra.AutoScrollbar(lf,
				takefocus=False,
				orient=HORIZONTAL)
		sbh.grid(row=1, column=0, sticky='ew')

		self.notes.config(  yscrollcommand=sbv.set,
				    xscrollcommand=sbh.set )
		sbv.config(command=self.notes.yview)
		sbh.config(command=self.notes.xview)

		lf.grid_rowconfigure(0, weight=1)
		lf.grid_columnconfigure(0, weight=1)

		self.title.bind("<Return>",   self.updateTitle)
		self.title.bind("<KP_Enter>", self.updateTitle)
		self.title.bind("<FocusOut>", self.updateTitle)

	#----------------------------------------------------------------------
	def selection(self, event=None):
		if self.notes.focus_get() is self.notes:
			try:
				self.notes_First = self.notes.index(SEL_FIRST)
				self.notes_Last  = self.notes.index(SEL_LAST)
			except TclError:
				self.notes_First = self.notes_Last = None

	#----------------------------------------------------------------------
	def createStatus(self):
		FlairRibbon.FlairPage.createStatus(self)
		self.dirLabel = Label(self.status, text="Dir:",
					justify=LEFT, anchor=W, width=10,
					relief=SUNKEN)
		self.dirLabel.pack(side=LEFT, expand=YES, fill=X)
		self.exeLabel = Label(self.status, text="Exe:",
					foreground="DarkBlue",
					justify=LEFT, anchor=W, width=10,
					relief=SUNKEN)
		self.exeLabel.pack(side=LEFT, expand=YES, fill=X)
		self.exeLabel.bind("<Button-1>", self.loadExecutable)
		return self.status

	# ----------------------------------------------------------------------
	# refresh contents of files
	# ----------------------------------------------------------------------
	def refresh(self):
		#say("ProjectPage.refresh")
		FlairRibbon.FlairPage.refresh(self)
		self.title.delete(0, END)
		self.title.insert(0, self.project.title)
		# Find styles remove the ones from the buttons
		styles = self.notes.styleList()
		for s in ("Bold", "Italic", "Underline", "Overstrike", "Left", "Center", "Right", "Link"):
			styles.remove(s)
		self.styleCombo.fill(styles)
		self.notes.setText(self.project.notes)
		self.checkStyle()

	# ----------------------------------------------------------------------
	def refreshButton(self):
		self.get()
		self.refresh()

	# ----------------------------------------------------------------------
	def refreshStatus(self):
		FlairRibbon.FlairPage.refreshStatus(self)
		self.dirLabel["text"] = "Dir: %s"%(self.project.dir)
		self.exeLabel["text"] = "Exe: %s"%(self.project.exe)

	# ----------------------------------------------------------------------
	# Load executable
	# ----------------------------------------------------------------------
	def loadExecutable(self, event=None):
		self.get()
		file = bFileDialog.askopenfilename(master=self.page,
			title="Load Executable",
			filetypes=[ ("All","*"), ("Exe","*.exe") ],
			initialfile=self.project.exe)
		if len(file) > 0:
			self.project.setExecutable(file)
			self.setModified()
			self.flair.refresh("exe")
			self.refreshStatus()

	# ----------------------------------------------------------------------
	def setStyle(self, style):
		self._style.set(style)

	# ----------------------------------------------------------------------
	def setStyleFromCombo(self):
		if self.notes_First is not None and self.notes.focus_get() is not self.notes:
			self.notes.tag_add(SEL, self.notes_First, self.notes_Last)
		self._style.set(self.styleCombo.get())
		self.notes.focus_set()

	# ----------------------------------------------------------------------
	def styleChanged(self, a=False, b=False, c=False):
		if self._setVar: return
		self.notes.setStyle(self._style.get())
		self.styleCombo.set(self._style.get())

	# ----------------------------------------------------------------------
	def checkStyle(self, event=None):
		self._setVar = True
		style = self.notes.styleAtCursor()
		self._style.set(style)
		self.styleCombo.set(style)
		self._setVar = False

	# ----------------------------------------------------------------------
	def editStyle(self, event=None):
		RichText.StyleConfigurationDialog(self.page, self.notes, self._style.get())
		self.notes.configureTags()

	# ----------------------------------------------------------------------
	def insertImage(self, event=None):
		filename = bFileDialog.askopenfilename(master=self.page,
			title="Insert image",
			filetypes=[
				("PNG - Portable Network Graphics", "*.png"),
				("GIF - Graphics Interchange Format", "*.gif"),
				("JPG - Joint Photographic Experts Group", ("*.jpg","*.jpeg")),
				("All","*")])

		if not filename: return
		filename = self.project.relativePath(filename)
		self.notes.insertImage(filename)

	# ----------------------------------------------------------------------
	def openImage(self, event=None):
		try:
			image = self.notes.image_cget(INSERT,"image")
		except TclError:
			try:
				image = self.notes.image_cget("%s - 1 char"%(INSERT),"image")
			except TclError:
				return "break"


		imageName = self.notes.imageName(image)
		if imageName is not None:
			# Switch to Plots page and activate the plot
			name,ext = os.path.splitext(imageName)
			self.flair.showPlot(name)
		return "break"

	# ----------------------------------------------------------------------
	def openLink(self, event=None):
		prev = self.notes.tag_prevrange("link", INSERT)
		if prev:
			from_,to_ = prev
			url = Text.get(self.notes, from_,to_)
			if "://" not in url: url = "http://%s"%(url)
			tkFlair.openurl(url)

	# --------------------------------------------------------------------
	def findDialog(self, event=None):
		try: txt = self.notes.selection_get()
		except:
			txt = ""
		fd = tkDialogs.FindReplaceDialog(self, replace=False)
		fd.show(self.find, None, None, txt)
		try:
			self.notes.tag_delete("find")
		except TclError:
			pass

	# --------------------------------------------------------------------
	def replaceDialog(self, event=None):
		try: txt = self.notes.selection_get()
		except: txt = ""
		fd = tkDialogs.FindReplaceDialog(self,replace=True)
		fd.show(self.find, self.replace, self.replaceAll, txt)
		self.notes.tag_delete("find")

	# --------------------------------------------------------------------
	def find(self, findStr=None, matchCase=None):
		global save_findStr, save_matchCase
		if findStr:
			save_findStr   = findStr
			save_matchCase = matchCase
		else:
			try:
				findStr   = save_findStr
				matchCase = save_matchCase
			except: return

		self.notes.tag_remove(SEL,"0.0",END)
		self.notes.tag_delete("find")
		try:
			index = self.notes.search(findStr,
				"insert + 1 chars",
				nocase = matchCase)
			self.notes.mark_set('insert', index)
			endex = '%s + %d chars' % ('insert', len(findStr))
			self.notes.tag_add(SEL, index, endex)
			self.notes.tag_add("find", index, endex)
			self.notes.tag_configure("find",
					background = '#800000',
					foreground = '#FFFFFF')
			self.notes.see(index)
			self.notes.update_idletasks()
			self.found = True
			return True
		except:
			self.flair.notify("Not Found",
				"Target \"%s\" not found" % (findStr))
			self.found = False
			return False

	# --------------------------------------------------------------------
	def findNext(self, event=None):
		self.find()

	# --------------------------------------------------------------------
	def replace(self, findStr, replaceStr, matchCase):
		if not self.found:
			if not self.find(findStr, matchCase): return
		index = self.notes.index(SEL_FIRST)
		self.notes.delete(SEL_FIRST, SEL_LAST)
		self.notes.insert(index, replaceStr)
		endex = "%s + %d chars" % (index, len(replaceStr))
		self.notes.tag_add(SEL, index, endex)
		self.notes.tag_add("find", index, endex)
		self.found = False

	# --------------------------------------------------------------------
	def replaceAll(self, findStr, replaceStr, matchCase):
		say("replaceAll", findStr, replaceStr, matchCase)

	# ----------------------------------------------------------------------
	# Update (Get) fields
	# ----------------------------------------------------------------------
	def get(self, event=None):
		if self.project is None: return
		self.updateTitle()

		#n = self.notes.txt.get(1.0, END).strip()
		n = self.notes.getText().strip()
		if n != self.project.notes:
			self.project.notes = n
			self.setModified(True)

	# ----------------------------------------------------------------------
	def updateTitle(self, event=None):
		t = self.title.get().strip()
		if t != self.project.title:
			self.project.title = t
			self.flair.refresh()
			self._updateInputTitle()
			self.setModified(True)
			self.flair.setInputModified(True)

	# ----------------------------------------------------------------------
	# update input title's
	# ----------------------------------------------------------------------
	def _updateInputTitle(self):
		if self.project.title == "": return

		try: titlecard = self.project.input.cards["TITLE"][0]
		except: return

		# Project title has been modified
		# ask user if he wants to update the titles
		if not tkFlair.askyesno("Update title?",
				"Project title is modified.\n"
				"Update TITLE card?",
				parent=self.title):
			return
		titlecard.setExtra(self.project.title)

	# ----------------------------------------------------------------------
	# Load Recent
	# ----------------------------------------------------------------------
	def loadRecent(self, r):
		if self.project.isModified():
			if not tkFlair.askyesno("Load recent project",
					"Current Project is modified\nDiscard changes",
					parent=self.page):
				return
		filename = tkFlair.config.get(tkFlair.__name__, \
					"recent.%d"%(r))
#					"recent.%d"%(r)).decode("utf-8")
		self.flair.loadProject(filename)

	def loadRecent0(self): self.loadRecent(0)
	def loadRecent1(self): self.loadRecent(1)
	def loadRecent2(self): self.loadRecent(2)
	def loadRecent3(self): self.loadRecent(3)
	def loadRecent4(self): self.loadRecent(4)
	def loadRecent5(self): self.loadRecent(5)
	def loadRecent6(self): self.loadRecent(6)
	def loadRecent7(self): self.loadRecent(7)
	def loadRecent8(self): self.loadRecent(8)
	def loadRecent9(self): self.loadRecent(9)

	# ----------------------------------------------------------------------
	def openTemplate(self, template=None):
		self.flair.newProject()
		if template is None:
			self.flair.openInputTemplate(tkFlair.inputTemplate())
		else:
			self.flair.openInputTemplate(template)

	# ----------------------------------------------------------------------
	def hardcopy(self, event=None):
		self.get()
		printer = tkFlair.printer(self.page)
		if printer:
			printer.open()
			printer.write("Project: %s\n"%(self.project.name))
			printer.write("Dir:   %s\n"%(self.project.dir))
			printer.write("Title: %s\n"%(self.project.title))
			printer.write("Input: %s\n"%(self.project.inputFile))
			printer.write("Executable: %s\n"%(self.project.exe))
			printer.write("Geometry: %s\n"%(self.project.input.geoFile))
			printer.write("Geometry Output: %s\n"%(self.project.input.geoOutFile))
			printer.write("Notes:\n%s\n"%(self.project.notes))
			printer.write("\n")
			printer.close()
