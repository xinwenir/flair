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
# Date:	12-Sep-2013

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
#from log import say

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import Manual
import Ribbon
import Unicode
import tkFlair
import tkExtra

#===============================================================================
# Tab pages inside a Page with a dedicated Ribbon
#===============================================================================
class FlairTab(Frame):
	#----------------------------------------------------------------------
	def __init__(self, master, page, **kw):
		Frame.__init__(self, master, **kw)
		self.page    = page
		self.flair   = page.flair
		self.project = page.project

	#----------------------------------------------------------------------
	def createRibbon(self):
		self.ribbon = Frame(self.page.ribbon,
					padx=0,
					pady=0,
					background=Ribbon._BACKGROUND)
		return self.ribbon

	#----------------------------------------------------------------------
	# Called when tab is activated
	#----------------------------------------------------------------------
	def activate(self):
		pass

	#----------------------------------------------------------------------
	# Called when tab is de-activated
	#----------------------------------------------------------------------
	def deactivate(self):
		pass

#===============================================================================
# Flair dynamic menu
#===============================================================================
class FlairMenu(Ribbon.MenuButton):
	def __init__(self, master, page, **kw):
		Ribbon.MenuButton.__init__(self, master, None, **kw)
		self.page = page

#===============================================================================
# Flair menu group
#===============================================================================
class FlairMenuGroup(Ribbon.MenuGroup):
	def __init__(self, master, page, name, menulist=None, **kw):
		Ribbon.MenuGroup.__init__(self, master, name, menulist, **kw)
		self.page = page

#===============================================================================
# Input template menu button
#===============================================================================
class _TemplateMenuButton(FlairMenu):
	def __init__(self, master, page, **kw):
		FlairMenu.__init__(self, master, page, **kw)
#		self.openproject = openproject

	#----------------------------------------------------------------------
	def createMenu(self):
		# Find templates
		files = {}
		for path in (os.path.join(tkFlair.prgDir, tkFlair._TEMPLATES_DIR),
			     os.path.join(tkFlair.iniDir, tkFlair._TEMPLATES_DIR)):
			try:
				for filename in os.listdir(path):
					(fn, ext) = os.path.splitext(filename)
					if ext == ".inp":
						files[fn] = path
			except OSError:
				pass
		filenames = list(files.keys())
		filenames.sort()

		menu = Menu(self, tearoff=0, activebackground=Ribbon._ACTIVE_COLOR)
		for fn in filenames:
			f = os.path.join(files[fn], fn+".inp")
#			if self.openproject:
			cmd = lambda f=f,p=self.page : p.openTemplate(f)
#			else:
#				cmd = lambda f=f,p=self.page : p.openInputTemplate(f)
			menu.add_command(label=fn,
				image=tkFlair.icons["input"],
				compound=LEFT,
				command=cmd)
		return menu

#===============================================================================
# Flair Tab Page
#===============================================================================
class FlairPage(Ribbon.Page):		# <--- should be possible to be a toplevel as well
	active = None

	_name_ = "Page"			# To be overridden by the end class
	_icon_ = "empty"

	#----------------------------------------------------------------------
	def __init__(self, flair, name, icon, show, **kw):
		Ribbon.Page.__init__(self, name, icon, show, **kw)
		self.flair   = flair
		self.project = flair.project

	#----------------------------------------------------------------------
	def createPage(self, scroll=True):
		Ribbon.Page.createPage(self)

		self._title = []

		# Do not show title when undocked
		if not self._undock:
			f = Frame(self.page, relief=GROOVE)
			f.grid(row=0, column=0, columnspan=2, sticky=NSEW)
			self._title.append(f)

			# --- Pin --- and close button ---
			self._pinbtn = Button(f, text=Unicode.BULLSEYE,
				command=self.pinButton,
				font=tkFlair._TITLE_FONT,
				border=0, highlightthickness=0,
				foreground=tkFlair._TITLE_FOREGROUND_COLOR,
				width=3,
				takefocus=False, padx=1, pady=1)
			tkExtra.Balloon.set(self._pinbtn, "Pin page to remain always visible")
			self.showPin()
			self._title.append(self._pinbtn)
			self._title[-1].pack(side=LEFT, fill=Y)

			# --- title ---
			if FlairPage.active is self:
				bg=tkFlair._TITLE_BACKACTIVE_COLOR
			else:
				bg=tkFlair._TITLE_BACKINACTIVE_COLOR

			self._title.append(Label(f, text=self.name,
					font=tkFlair._TITLE_FONT,
					foreground=tkFlair._TITLE_FOREGROUND_COLOR,
					background=bg,
					justify="center",
					padx=10,
					pady=1))
			self._title[-1].bind("<Button-1>",        self._dragStart)
			self._title[-1].bind("<B1-Motion>",       self._dragMove)
			self._title[-1].bind("<ButtonRelease-1>", self._dragEnd)
			self._title[-1].bind("<Double-1>",        self.minmaxButton)
			self._title[-1].pack(side=LEFT, expand=YES, fill=X)

			# --- close ---
			self._title.append(Button(f, text=tkExtra.CROSS,
					command=self.close,
					font=tkFlair._TITLE_FONT,
					border=0, highlightthickness=0,
					width=2,
					foreground=tkFlair._TITLE_FOREGROUND_COLOR,
					takefocus=False, padx=1, pady=1))
			self._title[-1].pack(side=RIGHT, fill=Y)
			tkExtra.Balloon.set(self._title[-1], "Close page and display the previous one")

			# --- max/min ---
			self._minmax = Button(f, text=Unicode.BLACK_UP_POINTING_TRIANGLE,
				command=self.minmaxButton,
				font=tkFlair._TITLE_FONT,
				border=0, highlightthickness=0,
				width=2,
				foreground=tkFlair._TITLE_FOREGROUND_COLOR,
				takefocus=False, padx=1, pady=1)
			self._title.append(self._minmax)
			self._title[-1].pack(side=RIGHT, fill=Y)
			tkExtra.Balloon.set(self._title[-1], "Maximize/Minimize page")

		else:
			self._pinbtn = None
			self._minmax = None

		# Do we need a scroll frame?
		if scroll:
			self._scrollFrame = tkExtra.ScrollFrame(self.page, takefocus=True)
			self._scrollFrame.grid(row=1, column=0, sticky=NSEW, padx=2, pady=2)

			self.frame = self._scrollFrame()
			self.frame.config(takefocus=True)

			sbx=tkExtra.AutoScrollbar(self.page, orient=HORIZONTAL,
					command=self._scrollFrame.xview)
			sbx.grid(row=2,column=0,sticky=NSEW)
			sby=tkExtra.AutoScrollbar(self.page, orient=VERTICAL,
					command=self._scrollFrame.yview)
			sby.grid(row=1,column=1,sticky=NSEW)
			self._scrollFrame["xscrollcommand"] = sbx.set
			self._scrollFrame["yscrollcommand"] = sby.set
			self._mark = None
			self._pos  = None

		else:	# Or fixed
			self._scrollFrame = None
			self.frame = Frame(self.page, takefocus=True)
			self.frame.grid(row=1, column=0, columnspan=2, sticky=NSEW)

		self.page.grid_columnconfigure(0, weight=1)
		self.page.grid_rowconfigure(1, weight=1)

		self.frame.bind("<FocusIn>", self._focusIn)
#		self.page.bind("<Map>", self.mapped)

		self._invalid = False
		self._tabs   = []
		self._tabVar = IntVar()
		self._tabVar.set(0)
		self.activeTab = None
#		self.log = None

	#----------------------------------------------------------------------
	def createRibbon(self):
		Ribbon.Page.createRibbon(self)
		for tab in self._tabs:
			tab.createRibbon()

	#----------------------------------------------------------------------
	def createStatus(self):
		Ribbon.Page.createStatus(self)

		self.inputLabel = Label(self.status, relief=SUNKEN,
				width=30,
				justify=LEFT, anchor=W, takefocus=0)
		self.inputLabel.pack(side=LEFT, fill=X)

		b = Button(self.status,image=tkFlair.icons["config"],
				command=self.preferences)
		b.pack(side=RIGHT)
		tkExtra.Balloon.set(b, "Open configuration dialog")
		b = Button(self.status,image=tkFlair.icons["output"],
				command=self.showOutput)
		b.pack(side=RIGHT)
		tkExtra.Balloon.set(b, "Show output")

		self.statusLabel = Label(self.status, relief=SUNKEN,
				justify=LEFT, anchor=W, takefocus=0)
		self.statusLabel.pack(side=RIGHT, expand=YES, fill=X)

		return self.status

	#----------------------------------------------------------------------
	def preferences(self):
		self.master.mastertab.event_generate("<<Preferences>>")

	#----------------------------------------------------------------------
	def showOutput(self):
		self.master.mastertab.event_generate("<<Output>>")

	#----------------------------------------------------------------------
	def addTab(self, tab):
		self._tabs.append(tab)

	#----------------------------------------------------------------------
	def changeTab(self, tab=None):
		if not self._tabs: return
		if tab is not None:
			self._tabVar.set(tab)
		if self.activeTab:
			self.activeTab.deactivate()
			self.activeTab.pack_forget()
			self.activeTab.ribbon.pack_forget()
		self.activeTab = self._tabs[self._tabVar.get()]
		self.activeTab.pack(expand=YES, fill=BOTH)
		self.activeTab.ribbon.pack(side=LEFT)
		self.activeTab.activate()

	# ----------------------------------------------------------------------
	def cut(self, event=None):
		f = self.page.focus_get()
		if event and f is event.widget: return	# Stop infinite loop
		if f: f.event_generate("<<Cut>>")

	# ----------------------------------------------------------------------
	def copy(self, event=None):
		f = self.page.focus_get()
		if event and f is event.widget: return	# Stop infinite loop
		if f: f.event_generate("<<Copy>>")

	# ----------------------------------------------------------------------
	def paste(self, event=None):
		f = self.page.focus_get()
		if event and f is event.widget: return	# Stop infinite loop
		if f: f.event_generate("<<Paste>>")

	# ----------------------------------------------------------------------
	def clone(self, event=None):
		self.copy(event)
		self.paste(event)

	# ----------------------------------------------------------------------
	def help(self, event=None):
		Manual.show(":%s:"%(self.name))

	# ----------------------------------------------------------------------
	def undo(self, event=None):
		f = self.page.focus_get()
		if event and f is event.widget: return	# Stop infinite loop
		if f: f.event_generate("<<Undo>>")

	# ----------------------------------------------------------------------
	def redo(self, event=None):
		f = self.page.focus_get()
		if event and f is event.widget: return	# Stop infinite loop
		if f: f.event_generate("<<Redo>>")

	# ----------------------------------------------------------------------
	def busy(self):
		try:
			self.flair.config(cursor="watch")
			self.flair.update_idletasks()
		except TclError:
			pass

	# ----------------------------------------------------------------------
	def notBusy(self):
		try:
			self.flair.config(cursor="")
		except TclError:
			pass

	# ----------------------------------------------------------------------
	def minmaxButton(self, event=None):
		self.activate()
		self.master.maximize(self)
		self.showminmax()

	# ----------------------------------------------------------------------
	def showminmax(self):
		if self._minmax is None: return
		if self.page is self.master._pageFrame.maxchild():
			self._minmax["text"] = Unicode.BLACK_DOWN_POINTING_TRIANGLE
		else:
			self._minmax["text"] = Unicode.BLACK_UP_POINTING_TRIANGLE

	# ----------------------------------------------------------------------
	def pinButton(self, event=None):
		self.pinned = not self.pinned
#		if not self.pinned: self.master.hide(self)
		self.showPin()

	# ----------------------------------------------------------------------
	def setpin(self, pin):
		Ribbon.Page.setpin(self, pin)
		self.showPin()

	# ----------------------------------------------------------------------
	def close(self):
		self.master.remove(self)

	# ----------------------------------------------------------------------
	def hide(self, replace):
		self.master.hide(self, replace)

	#----------------------------------------------------------------------
	def afterDock(self):
		self._invalid = True
		self.refreshTitle()

	# ----------------------------------------------------------------------
	def showPin(self):
		if self._pinbtn:
			if self.pinned:
				self._pinbtn.config(text=Unicode.FISHEYE)
			else:
				self._pinbtn.config(text=Unicode.BULLSEYE)

	# ----------------------------------------------------------------------
	# Show busy cursor and execute
	# ----------------------------------------------------------------------
	def execute(self):
		#self._setFocus()
		self.busy()
		self._execute()
		self.notBusy()

	# ----------------------------------------------------------------------
	# Override execute function
	# ----------------------------------------------------------------------
	def _execute(self):
		pass

	# ----------------------------------------------------------------------
	# Event handlers for dragging and placing
	# Bind to <Button-1>
	# ----------------------------------------------------------------------
	def _dragStart(self, event):
		self.activate()
		self.master._pageFrame.dragStart(event)

	# ----------------------------------------------------------------------
	# Bind to <B1-Motion>
	# ----------------------------------------------------------------------
	def _dragMove(self, event):
		self.master._pageFrame.dragMove(event)

	# ----------------------------------------------------------------------
	# Bind to <ButtonRelease-1>
	# ----------------------------------------------------------------------
	def _dragEnd(self, event):
		wasdragging = self.master._pageFrame._drag is not None
		if not self.master._pageFrame.dragEnd(event) and wasdragging:
			# Check if outside to undock it
			if event.x < -10 or event.y <-10 or \
			   event.x > self.master.winfo_width() + 10 or \
			   event.y > self.master.winfo_height() + 10:
				self.master.undock(self.name, event.x_root, event.y_root)

	# ----------------------------------------------------------------------
	# set invalid or refresh
	# ----------------------------------------------------------------------
	def setInvalid(self):
		if self.page and self.page.winfo_ismapped():
			self._invalid = False
			self.refresh()
		else:
			self._invalid = True

	# ----------------------------------------------------------------------
	def setValid(self):
		self._invalid = False

	# ----------------------------------------------------------------------
	def invalid(self): return self._invalid
	def valid(self): return not self._invalid

	# ----------------------------------------------------------------------
	def reset(self): pass

	# ----------------------------------------------------------------------
	# Refresh window contents
	# ----------------------------------------------------------------------
	def refresh(self):
		if self.page is None: return
		self.changeTab()
		self.refreshUndoButton()
		self._invalid = False

	#----------------------------------------------------------------------
	def refreshTitle(self):
		if self.page is None: return
		if self.project.name == "":
			title = "[untitled]"
		else:
			title = self.project.name
		try: self.master.master.title("%s - %s" % (title, self.name))
		except TclError: pass

	#----------------------------------------------------------------------
	def refreshStatus(self):
		if self.page is None: return
		if self.project.isInputModified():
			c = " +"
		else:
			c = ""
		self.inputLabel.config(text = "Inp: %s%s"%(self.project.inputFile,c))

	# ----------------------------------------------------------------------
	# Set status message
	# ----------------------------------------------------------------------
	def setStatus(self, txt, fg="#000000", bg=None):
		if self.page is None: return
		try:
			self.statusLabel.config(text=txt)
			self.statusLabel.config(fg=fg)
			if bg is None:
				self.statusLabel.config(bg=self.statusLabel.winfo_toplevel().cget("bg"))
			else:
				self.statusLabel.config(bg=bg)
		except TclError:
			pass

	# ----------------------------------------------------------------------
	# Redraw window contents
	# ----------------------------------------------------------------------
	redraw = refresh

	# ----------------------------------------------------------------------
	# Update contents
	# ----------------------------------------------------------------------
	def get(self):
		pass

	# ----------------------------------------------------------------------
	def setModified(self, m=True):
		self.flair.setModified(m)

	# ----------------------------------------------------------------------
	# Set focus, activate page if it is mapped
	# ----------------------------------------------------------------------
	def _focusIn(self, event=None):
		if self.page and self.page.winfo_ismapped():
			self.activate()

	# ----------------------------------------------------------------------
	def isactive(self):
		return FlairPage.active is self

	# ----------------------------------------------------------------------
	@staticmethod
	def clearActive():
		FlairPage.active = None

	# ----------------------------------------------------------------------
	def activate(self):
		#say("\nactivate",self.name)
		#import traceback
		#traceback.print_stack()
		if FlairPage.active is self: return

		if FlairPage.active is not None:
			FlairPage.active.deactivate()
		FlairPage.active = self

		for w in self._title:
			try: w["background"] = tkFlair._TITLE_BACKACTIVE_COLOR
			except TclError: pass
		self.showminmax()
		#if self._scrollFrame is not None:
		#	self._scrollFrame.update_scrollregion()
		if self.master.activePage.get() != self.name:
			self.master.changePage(self)

		# Get the focus if needed
		try:
			found = False
			focus = self.page.focus_get()
			while focus is not None:
				focus = focus.master
				if focus is self.page:
					found = True
					break
			if not found: self.page.focus_set()
		except TclError:
			pass

		self.refreshUndoButton()
		self.refreshStatus()

		if self._invalid: self.refresh()

		try: os.waitpid(-1,os.WNOHANG)	# Cleanup child processes
		except OSError: pass
		tkFlair.incStats(self.name)

	# ----------------------------------------------------------------------
	def deactivate(self):
		#print "\ndeactivate",self.name, self.page
		for w in self._title:
			try: w["background"] = tkFlair._TITLE_BACKINACTIVE_COLOR
			except TclError: pass
		#if self.page and self.page.winfo_ismapped(): self.get()
		if self.page: self.get()

	# ----------------------------------------------------------------------
	def configSave(self):
		pass
