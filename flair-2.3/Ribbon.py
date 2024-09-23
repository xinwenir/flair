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

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import time
import tkFlair
import tkExtra
import Unicode
from log import say

_TABFONT    = ("Sans","-14")
_FONT       = ("Sans","-11")

#_BACKGROUND_DISABLE = "#D6D2D0"
#_BACKGROUND         = "White"
#_BACKGROUND_GROUP   = "LightGray"
#_BACKGROUND_GROUP2  = "#D0E0E0"
#_FOREGROUND_GROUP   = "Black"

_BACKGROUND_DISABLE = "#A6A2A0"
_BACKGROUND         = "#E6E2E0"
_BACKGROUND_GROUP   = "#B6B2B0"

_BACKGROUND_GROUP2  = "#B0C0C0"
_BACKGROUND_GROUP3  = "#A0C0A0"
_BACKGROUND_GROUP4  = "#B0C0A0"

_FOREGROUND_GROUP   = "White"
_ACTIVE_COLOR       = "LightYellow"
_LABEL_SELECT_COLOR = "#C0FFC0"

# Ribbon show state
RIBBON_HIDDEN =  0	# Hidden
RIBBON_SHOWN  =  1	# Displayed
RIBBON_TEMP   = -1	# Show temporarily

#===============================================================================
# Frame Group with a button at bottom
#===============================================================================
class LabelGroup(Frame):
	def __init__(self, master, name, command=None, **kw):
		Frame.__init__(self, master, **kw)
		self.name = name
		self.config(	#bg="Green",
				background=_BACKGROUND,
				borderwidth=0,
				highlightthickness=0,
				pady=0)

		# right frame as a separator
		f = Frame(self, borderwidth=2, relief=GROOVE, background=_BACKGROUND_DISABLE)
		f.pack(side=RIGHT, fill=Y, padx=0, pady=0)

		# frame to insert the buttons
		self.frame = Frame(self,
				#bg="Orange",
				background=_BACKGROUND,
				padx=0,
				pady=0)
		self.frame.pack(side=TOP, expand=TRUE, fill=BOTH, padx=0, pady=0)

		if command:
			self.label = LabelButton(self, self, "<<%s>>"%(name), text=name)
			self.label.config(command=command,
				image=tkFlair.icons["triangle_down"],
				foreground=_FOREGROUND_GROUP,
				background=_BACKGROUND_GROUP,
				highlightthickness=0,
				borderwidth=0,
				pady=0,
				compound=RIGHT)
		else:
			self.label = Label(self, text=name,
					font       = _FONT,
					foreground = _FOREGROUND_GROUP,
					background = _BACKGROUND_GROUP,
					padx=2,
					pady=0)	# Button takes 1px for border width
		self.label.pack(side=BOTTOM, fill=X, pady=0)

	#-----------------------------------------------------------------------
	def grid2rows(self):
		self.frame.grid_rowconfigure(0, weight=1)
		self.frame.grid_rowconfigure(1, weight=1)

	#-----------------------------------------------------------------------
	def grid3rows(self):
		self.grid2rows()
		self.frame.grid_rowconfigure(2, weight=1)

#===============================================================================
class _KeyboardFocus:
	#-----------------------------------------------------------------------
	def _bind(self):
		self.bind("<Return>",		self._invoke)
		self.bind("<FocusIn>",		self._focusIn)
		self.bind("<FocusOut>",		self._focusOut)

	#-----------------------------------------------------------------------
	def _focusIn(self, event):
		self.__backgroundColor = self.cget("background")
		self.config(background = _ACTIVE_COLOR)

	#-----------------------------------------------------------------------
	def _focusOut(self, event):
		self.config(background = self.__backgroundColor)

	#-----------------------------------------------------------------------
	def _invoke(self, event):
		self.invoke()

#===============================================================================
# Button with Label that generates a Virtual Event or calls a command
#===============================================================================
class LabelButton(Button, _KeyboardFocus):
	def __init__(self, master, recipient=None, event=None, **kw):
		Button.__init__(self, master, **kw)
		self.config(	relief           = FLAT,
				activebackground = _ACTIVE_COLOR,
				font             = _FONT,
				borderwidth      = 1,
				highlightthickness = 0,
				padx             = 2,
				pady             = 0)
		_KeyboardFocus._bind(self)
		if recipient is not None:
			self.config(command = self.sendEvent)
			self._recipient = recipient
			self._event     = event
		else:
			self._recipient = None
			self._event     = None

	#-----------------------------------------------------------------------
	def sendEvent(self):
		self._recipient.event_generate(self._event)

#===============================================================================
class LabelCheckbutton(Checkbutton, _KeyboardFocus):
	def __init__(self, master, **kw):
		Checkbutton.__init__(self, master, **kw)
		self.config(	selectcolor        = _LABEL_SELECT_COLOR,
				activebackground   = _ACTIVE_COLOR,
				background         = _BACKGROUND,
				indicatoron        = FALSE,
				relief             = FLAT,
				borderwidth        = 0,
				highlightthickness = 0,
				padx               = 0,
				pady               = 0,
				font               = _FONT
			)
		_KeyboardFocus._bind(self)

#===============================================================================
class LabelRadiobutton(Radiobutton, _KeyboardFocus):
	def __init__(self, master, **kw):
		Radiobutton.__init__(self, master, **kw)
		self.config(
			selectcolor        = _LABEL_SELECT_COLOR,
			activebackground   = _ACTIVE_COLOR,
			background         = _BACKGROUND,
			indicatoron        = FALSE,
			borderwidth        = 0,
			highlightthickness = 0,
			pady               = 0,
			font               = _FONT
		)
		_KeyboardFocus._bind(self)

#===============================================================================
class LabelCombobox(tkExtra.Combobox, _KeyboardFocus):
	def __init__(self, master, **kw):
		tkExtra.Combobox.__init__(self, master, **kw)
		self.config(background=_BACKGROUND, font=_FONT)
		Frame.config(self, background=_BACKGROUND, padx=0, pady=0)
		_KeyboardFocus._bind(self)

	#-----------------------------------------------------------------------
	def _focusOut(self, event):
		self.config(background = _BACKGROUND) #self.__backgroundColor)
		Frame.config(self, background= _BACKGROUND) #self.__backgroundColor)

#===============================================================================
# Button with Label that popup a menu
#===============================================================================
class MenuButton(Button, _KeyboardFocus):
	def __init__(self, master, menulist, **kw):
		Button.__init__(self, master, **kw)
		self.config(	relief           = FLAT,
				activebackground = _ACTIVE_COLOR,
				font             = _FONT,
				borderwidth      = 0,
				highlightthickness= 0,
				padx             = 2,
				pady             = 0,
				command          = self.showMenu)

		_KeyboardFocus._bind(self)
		self.bind("<Return>", self.showMenu)
		if menulist is not None:
			self._menu = self.createMenuFromList(menulist)
		else:
			self._menu = None

	#-----------------------------------------------------------------------
	def showMenu(self, event=None):
		if self._menu is not None:
			self._showMenu(self._menu)
		else:
			self._showMenu(self.createMenu())

	#-----------------------------------------------------------------------
	def _showMenu(self, menu):
		if menu is not None:
			menu.tk_popup(
				self.winfo_rootx(),
				self.winfo_rooty() + self.winfo_height())

	#-----------------------------------------------------------------------
	def createMenu(self):
		return None

	#-----------------------------------------------------------------------
	def createMenuFromList(self, menulist):
		mainmenu = menu = Menu(self, tearoff=0, activebackground=_ACTIVE_COLOR)
		for item in menulist:
			if item is None:
				menu.add_separator()
			elif isinstance(item,str):
				menu = Menu(mainmenu)
				mainmenu.add_cascade(label=item, menu=menu)
			else:
				name, icon, cmd = item
				if icon is None: icon = "empty"
				menu.add_command(label=name,
						image=tkFlair.icons[icon],
						compound=LEFT,
						command=cmd)
		return menu

#===============================================================================
# A label group with a drop down menu
#===============================================================================
class MenuGroup(LabelGroup):
	def __init__(self, master, name, menulist=None, **kw):
		LabelGroup.__init__(self, master, name, command=self._showMenu, **kw)
		self._menulist = menulist

	#-----------------------------------------------------------------------
	def createMenu(self):
		if self._menulist is not None:
			return self.createMenuFromList(self._menulist)
		else:
			return None

	#-----------------------------------------------------------------------
	def createMenuFromList(self, menulist):
		mainmenu = menu = Menu(self, tearoff=0, activebackground=_ACTIVE_COLOR)
		for item in menulist:
			if item is None:
				menu.add_separator()
			elif isinstance(item,str):
				menu = Menu(mainmenu)
				mainmenu.add_cascade(label=item, menu=menu)
			else:
				name, icon, cmd = item
				if icon is None: icon = "empty"
				menu.add_command(label=name,
						image=tkFlair.icons[icon],
						compound=LEFT,
						command=cmd)
		return menu

	#-----------------------------------------------------------------------
	def _showMenu(self):
		menu = self.createMenu()
		if menu is not None:
			menu.tk_popup(
				self.winfo_rootx(),
				self.winfo_rooty() + self.winfo_height())

#===============================================================================
# Context group for a specific item in the Ribbon
#===============================================================================
#class ContextGroup(LabelGroup):
#	def __init__(self, master, name, **kw):
#		LabelGroup.__init__(self, master, name, **kw)

#===============================================================================
# Page Tab buttons
#===============================================================================
class TabButton(Radiobutton):
	def __init__(self, master, **kw):
		Radiobutton.__init__(self, master, **kw)
		self.config(	selectcolor        = _BACKGROUND,
				activebackground   = _ACTIVE_COLOR,
				indicatoron        = FALSE,
				relief             = FLAT,
				font               = _TABFONT,
				borderwidth        = 0,
				highlightthickness = 0,
				padx               = 5,
				pady               = 0,
				background         = _BACKGROUND_DISABLE
			)
		self.bind("<FocusIn>",		self._focusIn)
		self.bind("<FocusOut>",		self._focusOut)

	#-----------------------------------------------------------------------
	# Bind events on TabFrame
	#----------------------------------------------------------------------
	def bindClicks(self, tabframe):
		self.bind("<Double-1>",         tabframe.double)
		self.bind("<Button-1>",         tabframe.dragStart)
		self.bind("<B1-Motion>",        tabframe.drag)
		self.bind("<ButtonRelease-1>",  tabframe.dragStop)
		self.bind("<Control-ButtonRelease-1>", tabframe.pinActive)

		self.bind("<Left>",		tabframe._tabLeft)
		self.bind("<Down>",		tabframe._tabDown)

	#----------------------------------------------------------------------
	def _focusIn(self, evenl=None):
		self.config(selectcolor = _ACTIVE_COLOR)

	#----------------------------------------------------------------------
	def _focusOut(self, evenl=None):
		self.config(selectcolor = _BACKGROUND)

#===============================================================================
# Page
#===============================================================================
class Page:		# <--- should be possible to be a toplevel as well
	_motionClasses = (LabelButton, LabelRadiobutton, LabelCheckbutton, LabelCombobox, MenuButton)

	def __init__(self, name, icon, show, **kw):
		self.master    = None
		self.name      = name
		self._icon     = icon
#		self._packed   = False	# If button is packed due to lack of space
		self._showTab  = show
		self.pinned    = False
		self._undock   = False
#		self._shortcut = None
		self._previous = None	# Insert after previous if any
		self.width     = 640	# Window dimensions when undocked
		self.height    = 480
		self._tab      = None	# Tab button
		self.ribbon    = None	# Ribbon frame
		self.status    = None	# Status frame
		self.page      = None	# Page frame
		self._context  = {}	# dictionary of context groups
		self._refreshOn= ()	# Refresh trigger events
		self.init()

	#----------------------------------------------------------------------
	# Override initialization
	#----------------------------------------------------------------------
	def init(self):
		pass

	#----------------------------------------------------------------------
	# The tab page can change master if undocked
	#----------------------------------------------------------------------
	def create(self, master):
		if self.page is not None and self.master is master: return
		self.master = master
		self.createPage()
		self.createRibbon()
		self.createStatus()
		self.createContextGroups()
		self.ribbonBindMotion()
		self.refresh()

	#----------------------------------------------------------------------
	def createPage(self):
		self.page = Frame(self.master._pageFrame)
		return self.page

	#----------------------------------------------------------------------
	def createRibbon(self):
		self.ribbon = Frame(self.master, pady=0, background=_BACKGROUND)
		return self.ribbon

	#----------------------------------------------------------------------
	def ribbonBindMotion(self, widget=None):
		if widget is None: widget = self.ribbon
		for child in widget.winfo_children():
			for class_ in Page._motionClasses:
				if isinstance(child, class_):
					child.bind("<Left>",  self._ribbonLeft)
					child.bind("<Right>", self._ribbonRight)
					child.bind("<Up>",    self._ribbonUp)
					child.bind("<Down>",  self._ribbonDown)
			else:
				self.ribbonBindMotion(child)

	#----------------------------------------------------------------------
	def createStatus(self):
		self.status = Frame(self.master) #background=_BACKGROUND,
#					relief=SUNKEN)
		return self.status

	#----------------------------------------------------------------------
	def createContextGroups(self):
		pass

	#----------------------------------------------------------------------
	# Add context group
	#----------------------------------------------------------------------
	def addContextGroup(self, name, group):
		self._context[name] = group

	#----------------------------------------------------------------------
	# Show context group
	#----------------------------------------------------------------------
	def showContextGroup(self, name=None, append=False):
		for n, g in list(self._context.items()):
			if name==n:
				g.pack(side=LEFT, fill=Y, padx=0, pady=0)
			elif not append:
				g.pack_forget()

	#----------------------------------------------------------------------
	def setpin(self, pin):
		self.pinned = pin

	#----------------------------------------------------------------------
	def created(self):
		return self.page is not None

	#----------------------------------------------------------------------
	def ismapped(self):
		return self.page is not None and self.page.winfo_ismapped()

	#----------------------------------------------------------------------
	# Called when a page is activated
	#----------------------------------------------------------------------
	def activate(self):
		pass

	#----------------------------------------------------------------------
	# Called before the page is docked or undocked
	#----------------------------------------------------------------------
	def beforeDock(self):
		pass

	#----------------------------------------------------------------------
	# Called after the page is docked or undocked
	#----------------------------------------------------------------------
	def afterDock(self):
		pass

	#----------------------------------------------------------------------
	def dock(self):
		self.master.dock()

	# ----------------------------------------------------------------------
	def isdocked(self):	return not self._undock
	def isundocked(self):	return self._undock
	def canUndo(self):	return True
	def canRedo(self):	return True
	def resetUndo(self):	pass
	def undo(self, event=None): pass
	def redo(self, event=None): pass
	def refresh(self):	pass

	# ----------------------------------------------------------------------
	def refreshUndoButton(self):
		# Check if frame provides undo/redo
		if self.master is None: return
		if self.page is None: return

		if self.canUndo():
			state = NORMAL
		else:
			state = DISABLED
		try:
			self.master.tool["undo"].config(state=state)
			self.master.tool["undolist"].config(state=state)
		except TclError:
			pass

		if self.canRedo():
			state = NORMAL
		else:
			state = DISABLED
		try:
			self.master.tool["redo"].config(state=state)
		except TclError:
			pass

	#----------------------------------------------------------------------
	def keyboardFocus(self):
		self._tab.focus_set()

	#----------------------------------------------------------------------
	# Return the closest widget along a direction
	#----------------------------------------------------------------------
	@staticmethod
	def __compareDown(x,y,xw,yw):  return yw>y+1
	@staticmethod
	def __compareUp(x,y,xw,yw):    return yw<y-1
	@staticmethod
	def __compareRight(x,y,xw,yw): return xw>x+1
	@staticmethod
	def __compareLeft(x,y,xw,yw):  return xw<x-1

	#----------------------------------------------------------------------
	@staticmethod
	def __closest(widget, compare, x, y):
		closest = None
		dc2 = 10000000
		if widget is None: return closest, dc2
		for child in widget.winfo_children():
			for class_ in Page._motionClasses:
				if isinstance(child, class_):
					if child["state"] == DISABLED: continue
					xw = child.winfo_rootx()
					yw = child.winfo_rooty()
					if compare(x,y,xw,yw):
						d2 = (xw-x)**2 + (yw-y)**2
						if d2 < dc2:
							closest = child
							dc2 = d2
					break
			else:
				c,d2 = Page.__closest(child, compare, x, y)
				if d2 < dc2:
					closest = c
					dc2 = d2
		return closest, dc2

	#----------------------------------------------------------------------
	# Select/Focus the closest element
	#----------------------------------------------------------------------
	def _ribbonUp(self, event=None):
		x = event.widget.winfo_rootx()
		y = event.widget.winfo_rooty()
		closest,d2 = Page.__closest(self.ribbon, Page.__compareUp, x, y)
		if closest is not None:
			closest.focus_set()

	#----------------------------------------------------------------------
	def _ribbonDown(self, event=None):
		x = event.widget.winfo_rootx()
		y = event.widget.winfo_rooty()
		closest,d2 = Page.__closest(self.ribbon, Page.__compareDown, x, y)
		if closest is not None:
			closest.focus_set()

	#----------------------------------------------------------------------
	def _ribbonLeft(self, event=None):
		x = event.widget.winfo_rootx()
		y = event.widget.winfo_rooty()
		closest,d2 = Page.__closest(self.ribbon, Page.__compareLeft, x, y)
		if closest is not None:
			closest.focus_set()

	#----------------------------------------------------------------------
	def _ribbonRight(self, event=None):
		x = event.widget.winfo_rootx()
		y = event.widget.winfo_rooty()
		closest,d2 = Page.__closest(self.ribbon, Page.__compareRight, x, y)
		if closest is not None:
			closest.focus_set()

#===============================================================================
# TabRibbonFrame
#===============================================================================
class TabRibbonFrame(Frame):
	def __init__(self, master, mastertab=None, **kw):
		Frame.__init__(self, master, kw)
		self.config(background=_BACKGROUND_DISABLE)

		self.oldActive  = None	# remember old active page
		self.activePage = StringVar(self)
		self.tool       = {}
		self.pages      = {}
		self.undocked   = []
		self.history    = []
		self.hideTime   = 0

#		self.winfo_toplevel().bind("<ButtonRelease-1>", self.release)

		if mastertab is not None:
			self.mastertab = mastertab	# Master tab frame for Undocked frames
		else:
			self.mastertab = master

		self.defaultPage    = None
		self.defaultDynamic = None

		# === Top frame with buttons ===
		frame = Frame(self, background=_BACKGROUND_DISABLE)
		frame.pack(side=TOP, fill=X)

		# --- Basic buttons ---
		b = LabelButton(frame, self, "<<Save>>",
				image=tkFlair.icons["save"],
				background=_BACKGROUND_DISABLE)
		tkExtra.Balloon.set(b, "Save project [Ctrl-S]")
		b.pack(side=LEFT)

		b = LabelButton(frame, self, "<<Undo>>",
				image=tkFlair.icons["undo"],
				background=_BACKGROUND_DISABLE)
		tkExtra.Balloon.set(b, "Undo [Ctrl-Z]")
		b.pack(side=LEFT)
		self.tool["undo"] = b

		b = LabelButton(frame, image=tkFlair.icons["triangle_down"],
				command=self.undolist,
				background=_BACKGROUND_DISABLE)
		b.pack(side=LEFT)
		self.tool["undolist"] = b

		b = LabelButton(frame, self, "<<Redo>>",
				image=tkFlair.icons["redo"],
				background=_BACKGROUND_DISABLE)
		tkExtra.Balloon.set(b, "Redo [Ctrl-Y]")
		b.pack(side=LEFT)
		self.tool["redo"] = b

		Label(frame, image=tkFlair.icons["sep"],
				background=_BACKGROUND_DISABLE).pack(side=LEFT, padx=3)

		# --- Help ---
		b = LabelButton(frame, self, "<<Help>>",
				image=tkFlair.icons["info"],
				background=_BACKGROUND_DISABLE)
		tkExtra.Balloon.set(b, "Help [F1]")
		b.pack(side=RIGHT, padx=2)

		if self.mastertab is self.master:
			b = Button(frame, text=Unicode.BLACK_DOWN_POINTING_TRIANGLE,
					command          = self.dynamicMenu,
					relief           = FLAT,
					activebackground = _ACTIVE_COLOR,
					background       = _BACKGROUND_DISABLE,
					highlightbackground = _BACKGROUND_DISABLE,
					takefocus        = FALSE,
					borderwidth      = 1,
					width            = 2,
					padx             = 0,
					pady             = 0 )
			b.pack(side=RIGHT)

		self.dynamic = TabButton(frame,
				image=tkFlair.icons["empty"], text=" ", compound=LEFT,
				value    = "",
				variable = self.activePage,
				command  = self.changePage)
		if self.mastertab is self.master:
			self.dynamic.pack(side=RIGHT, fill=Y)
			self.dynamic.bindClicks(self)

		Label(frame, image=tkFlair.icons["sep"],
				background=_BACKGROUND_DISABLE).pack(side=RIGHT, padx=3)

		# --- TabBar ---
		self._tabFrame = Frame(frame, background=_BACKGROUND_DISABLE)
		self._tabFrame.pack(side=LEFT, fill=BOTH, expand=YES)

		# === Ribbon Frame with permanent and dynamic ribbons ===
		self._allRibbons = Frame(self, background=_BACKGROUND, pady=0)
		self._allRibbons.pack(side=TOP, fill=BOTH, pady=0)

		# --- Permanent Groups ---
		self.createClipboardGroup(self._allRibbons)

		# ==== Ribbon Frame ====
		self._ribbonFrame = Frame(self._allRibbons,
						#bg="Yellow",
						background=_BACKGROUND,
						relief=RAISED)
		self._ribbonFrame.pack(side=LEFT, fill=BOTH, padx=0, pady=0)
		self._ribbonFrameStatus = RIBBON_SHOWN

		# ==== Page Frame ====
		self._pageFrame = tkExtra.TreeSplitter(self) #, relief=SUNKEN, borderwidth=2)
		self._pageFrame.pack(fill=BOTH, expand=YES)

#		self.bind("<Configure>", self.configure)
		self._dragging = 0

	#----------------------------------------------------------------------
	def createClipboardGroup(self, frame):
		group = LabelGroup(frame, "Clipboard")
		group.pack(side=LEFT, fill=BOTH, padx=0, pady=0)
		group.grid2rows()

		# ---
		b = LabelButton(group.frame, self, "<<Paste>>",
				image=tkFlair.icons["paste32"],
				text="Paste",
				compound=TOP,
				takefocus=FALSE,
				background=_BACKGROUND)
		b.grid(row=0, column=0, rowspan=2, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Paste [Ctrl-V]")
		#self.tool["paste"] = b

		# ---
		b = LabelButton(group.frame, self, "<<Cut>>",
				image=tkFlair.icons["cut"],
				text="Cut",
				compound=LEFT,
				anchor=W,
				takefocus=FALSE,
				background=_BACKGROUND)
		tkExtra.Balloon.set(b, "Cut [Ctrl-X]")
		b.grid(row=0, column=1, padx=0, pady=1, sticky=NSEW)

		# ---
		b = LabelButton(group.frame, self, "<<Copy>>",
				image=tkFlair.icons["copy"],
				text="Copy",
				compound=LEFT,
				anchor=W,
				takefocus=FALSE,
				background=_BACKGROUND)
		tkExtra.Balloon.set(b, "Copy [Ctrl-C]")
		b.grid(row=1, column=1, padx=0, pady=1, sticky=NSEW)

	#----------------------------------------------------------------------
	def undolist(self, event=None): self.event_generate("<<UndoList>>")

	#----------------------------------------------------------------------
	def getActivePage(self):
		return self.pages[self.activePage.get()]

	#----------------------------------------------------------------------
	# Add page to the tabs
	#----------------------------------------------------------------------
	def addPage(self, page, toplevel=False):
		if page._showTab:
			if self.defaultPage is None:
				self.defaultPage = page
				self.activePage.set(page.name)
		else:
			if self.defaultDynamic is None:
				self.defaultDynamic = page

		if page.name not in self.pages:
			self.pages[page.name] = page
		if toplevel or not page._undock:
			if page._showTab:
				self._addTab(page)
		else:
			self._addUndocked(page)

	# ----------------------------------------------------------------------
	# Add page as tab
	# ----------------------------------------------------------------------
	def _addTab(self, page):
		page._tab = TabButton(self._tabFrame,
				image    = page._icon,
				text     = page.name,
				compound = LEFT,
				value    = page.name,
				variable = self.activePage,
				command  = self.changePage)
		tkExtra.Balloon.set(page._tab,
				"%s\nDrag to move or undock" % (page.__doc__))

		if page._previous and not page._undock:
			# insert after previous if it exists
			slaves = self._tabFrame.pack_slaves()
			for i in range(len(slaves)):
				if page._previous == slaves[i]["text"]:
					page._tab.pack(after=slaves[i],
						side=LEFT, fill=Y, padx=5)
					break
			else:
				page._tab.pack(side=LEFT, fill=Y, padx=5)
		else:
			page._tab.pack(side=LEFT, fill=Y, padx=5)
		if self.mastertab is self.master:
			page._tab.bindClicks(self)

	# ----------------------------------------------------------------------
	# Add a page as a toplevel Undocked
	# ----------------------------------------------------------------------
	def _addUndocked(self, page, rootx=None, rooty=None):
		toplevel = Toplevel(self)
		toplevel.pack_propagate(0)
		toplevel.title(page.name)
		try:
			toplevel.tk.call('wm','iconphoto',toplevel._w, page._icon)
		except TclError:
		#	#self.iconbitmap("@%s/icons/marmite.xbm"%(tkFlair.prgDir))
			pass

		tabframe = self.__class__(toplevel, self)
		tabframe.pack(side=TOP, fill=BOTH, expand=YES)
		page._showTab = True
		page.create(tabframe)
		tabframe.addPage(page, True)
		toplevel.protocol("WM_DELETE_WINDOW", tabframe.dock)

		if rootx is not None:
			toplevel.geometry("+%d+%d"%(rootx,rooty))
		if page.width > 0:
			toplevel.geometry("%dx%d" % (page.width, page.height))
		toplevel.wait_visibility()
		tabframe.changePage(page)	# only at the end
						# otherwise it steals the focus from other
						# pages like the GeometryEditor and
						# forces a get() while still incomplete
		return toplevel

	# ----------------------------------------------------------------------
	# Remove page if mapped
	# ----------------------------------------------------------------------
	def remove(self, page):
		# find another active page
		if page.page is None: return
		page.setpin(False)
		self._pageFrame.remove(page.page)
		page.ribbon.pack_forget()
		page.status.pack_forget()

		if self._pageFrame.isempty(): # Pop from history
			if self.history and self.history[-1] == page.name:
				self.history.pop()
				if self.history:
					self.changePage(self.history[-1], False)
				else:
					self.changePage(self.defaultPage, False)
				return
		else:
			# Find page that is maximized
			for page in list(self.pages.values()):
				if page.page is not None and \
				   page.isdocked() and \
				   page.page.winfo_ismapped():
					page.activate()
					return

		# If everything else fails
		self.changePage(self.defaultPage, False)

	# ----------------------------------------------------------------------
	# Hide or remove a page from the display
	# ----------------------------------------------------------------------
	def hide(self, page, replace=False):
		# find another active page
		if page.page is None: return
		#waspinned = page.pinned
		if not replace: self._pageFrame.remove(page.page)
		page.setpin(False)
		page.ribbon.pack_forget()
		page.status.pack_forget()

		if self.history and self.history[-1] == page.name:
			self.history.pop()
			if self.history:
				self.changePage(self.history[-1], False)
			else:
				self.changePage(self.defaultPage, False)
		else:
			self.changePage(self.defaultPage, False)

	# ----------------------------------------------------------------------
	# undock tab with pagename
	# ----------------------------------------------------------------------
	def undock(self, pagename=None, x_root=None, y_root=None):
		if pagename is None: pagename = self.activePage.get()
		page = self.pages[pagename]
		if page is self.defaultPage: return

		if page.created():
			self.remove(page)
		else:
			page.create(self)

		page.beforeDock()

		# find tab before
		if self.dynamic["text"] == page.name:
			self.changeDynamic()
			page._previous = None
		else:
			slaves = self._tabFrame.pack_slaves()
			for i in range(len(slaves)):
				if page.name == slaves[i]["text"]:
					page._previous = slaves[i-1]["text"]
					break
			else:
				page._previous = None

		# Remove tab
		for tab in self._tabFrame.pack_slaves():
			if tab["text"] == page.name:
				tab.pack_forget()
				break

		# create undocked window
		page._undock = True
		self._addUndocked(page, x_root, y_root)
		page.afterDock()

	# ----------------------------------------------------------------------
	# Dock window
	# ----------------------------------------------------------------------
	def dock(self):
		if self.master is self.mastertab: return
		page = list(self.pages.values())[0]
		page.beforeDock()
		page.deactivate()
		page._undock = False
		if page._previous is None:
			page._showTab = False
		page.width  = self.master.winfo_width()
		page.height = self.master.winfo_height()
		page.ribbon = None
		page.status = None
		page.page   = None
		self.mastertab.addPage(page)
		self.master.destroy()
		page.afterDock()
#		if page._previous is None:
#			self.mastertab.changeDynamic(page.name)
#		self.mastertab.changePage(page)

	# ----------------------------------------------------------------------
	# Dock all undocked windows
	# ----------------------------------------------------------------------
	def dockAll(self):
		for page in list(self.pages.values()):
			if page.page and page.isundocked():
				page.dock()

	#----------------------------------------------------------------------
	# Move from the dynamic tab to the user tabs
	#----------------------------------------------------------------------
	def dynamic2tab(self):
		count = len([x for x in list(self.pages.values()) if not x._showTab])
		if count == 0: return
		page = self.pages[self.dynamic["text"]]
		self.addPage(page)
		page._showTab = True
		self.dynamic["text"] = ""
		self.changeDynamic()

	#----------------------------------------------------------------------
	# Move from user tabs to dynamic
	#----------------------------------------------------------------------
	def tab2dynamic(self, tab):
		page = self.pages[tab["text"]]
		page._showTab = False
		tab.pack_forget()
		self.changeDynamic(page.name)

#	def removePage(self, page):
#		if page.name == self.activePage.get():
#			self._pageFrame.remove(page.page)
#			self.changePage(self.defaultPage)
#		else:
#			self._pageFrame.remove(page.page)

	# ----------------------------------------------------------------------
	def addHistory(self, page):
		if self.history and self.history[-1] == page: return
		self.history.append(page)
		if len(self.history) > 10:
			self.history.pop(0)

	# ----------------------------------------------------------------------
	# Change ribbon and page
	# ----------------------------------------------------------------------
	def changePage(self, page=None, add2History=True):
		#say("\n----> changePage",page,self.activePage.get())
		#import traceback
		#traceback.print_stack()

		if page:
			if isinstance(page, str):
				try:
					page = self.pages[page]
				except KeyError:
					return
			self.activePage.set(page.name)
		else:
			try:
				page = self.pages[self.activePage.get()]
			except KeyError:
				return

		if add2History: self.addHistory(page.name)

		# Remember top splitter orientation
		if self._pageFrame.tree:
			hori = not self._pageFrame.tree.hori
			pos  = self._pageFrame.tree.pos
		else:
			hori = True
			pos  = 0.5

		if page._undock and self.master is self.mastertab:
			toplevel = page.master.winfo_toplevel()
			toplevel.deiconify()
			toplevel.lift()
			page.master.focus_set()
			self.oldActive = page
			return

		if page.page is None:
			# Create page if needed
			page.create(self)
		else:
			# Check if current page is already maximized
			if self._pageFrame.maxchild() is page.page:
				self.oldActive = page
				return

		# Check if old active was unpinned to replace this one
		node = self._pageFrame.node(page.page)
		if node is None and self.oldActive:
			if not self.oldActive.pinned and not self.oldActive._undock:
				node = self._pageFrame.node(self.oldActive.page)

		# Remove unwanted pages
		for p in list(self.pages.values()):
			if p is page: continue
			if p.page is None: continue
			if not self._pageFrame.node(p.page): continue

			if p._undock: continue
			if self.mastertab is not self.master: continue

			if p._showTab:
				if p.ribbon: p.ribbon.pack_forget()
				if p.status: p.status.pack_forget()
			if not p.pinned and node is None:
				# Node to be replaced
				node = self._pageFrame.node(p.page)

		# Show ribbon
		if page.ribbon is None:
                        self.create(self,master)
                        
		if not page.ribbon.winfo_ismapped():
			page.ribbon.pack(in_=self._ribbonFrame, fill=BOTH, expand=YES)

		if not page.status.winfo_ismapped():
			page.status.pack(side=BOTTOM, fill=X)

		# Show page 
		if not page.page.winfo_ismapped():
			# Show page
			if node is None:
				self._pageFrame.add(None, page.page, pos, hori)
			else:
				self._pageFrame.replace(node, page.page)
			page.ribbon.pack(in_=self._ribbonFrame, fill=BOTH, expand=YES)
			page.status.pack(side=BOTTOM, fill=X)

		# show in dynamic if needed, and remove old one
		if not page._showTab:
			oldpage = self.pages.get(self.dynamic["text"])
			if oldpage:
				oldpage._showTab = False
				oldpage._tab = None
			self.dynamic.config(text=page.name, image=page._icon, value=page.name)
			tkExtra.Balloon.set(self.dynamic,
				"%s\nDrag to move or undock" % (page.__doc__))
			page._tab = self.dynamic
			page._showTab = True

		if self._pageFrame.maxchild() and self._pageFrame.maxchild() is not page.page:
			self._pageFrame.maximize(page.page)
		page.activate()
		self.event_generate("<<ChangePage>>", data=page.name)
		self.oldActive = page

	#----------------------------------------------------------------------
	# Change the dynamic to name
	#----------------------------------------------------------------------
	def changeDynamic(self, name=None):
		#say("----> changeDynamic",name)
		if self.mastertab is not self.master: return

		if name is None:
			# Prefer default
			page = self.defaultDynamic
			if page._showTab:
				# Search for the first non-displayed item
				pages = [p.name for p in list(self.pages.values())
						if not p._showTab]
				if pages:	# Choose one random
					page = self.pages[pages[0]]
				else:
					return

		else:
			# Change to new name
			page = self.pages[name]

		# old dynamic page
		oldpage = self.pages.get(self.dynamic["text"])

		if oldpage:
			# If it is the same then simple activate it
			if oldpage.name == page.name:
				self.changePage(page)
				return

			# Remove old dynamic if needed
			if oldpage.page and  not oldpage.pinned and \
			  self._pageFrame.node(page.page) is not None:
				self._pageFrame.remove(oldpage.page)
				oldpage.ribbon.pack_forget()
				oldpage.status.pack_forget()

		# Set new dynamic page
		self.dynamic.config(text=page.name, image=page._icon, value=page.name)
		tkExtra.Balloon.set(self.dynamic,
				"%s\nDrag to move or undock" % (page.__doc__))
		page._tab = self.dynamic
		page._showTab = True

		if name:
			self.changePage(page)
			self.activePage.set(page.name)

		if oldpage:
			oldpage._showTab = False
			oldpage._tab = None

	#----------------------------------------------------------------------
	def maximize(self, page):
		# Page should be visible
		self._pageFrame.maximize(page.page)

	#----------------------------------------------------------------------
	# Show the page menu with the non-shown pages
	#----------------------------------------------------------------------
	def dynamicMenu(self, event=None):
		menu = Menu(self, tearoff=0, activebackground=_ACTIVE_COLOR)
		menu.bind("<Left>", lambda e : e.widget.destroy())
		menu.bind("<Right>", lambda e : e.widget.destroy())

		pages = [p.name for p in list(self.pages.values()) if not p._showTab]
		pages.append(self.dynamic["text"])
		pages.sort()

		for name in pages:
			page = self.pages[name]
			menu.add_command(label=name,
				image=page._icon,
				compound=LEFT,
				command=lambda n=name,s=self : s.changeDynamic(n))

		menu.tk_popup(self.dynamic.winfo_rootx(),
				self.dynamic.winfo_rooty() + self.dynamic.winfo_height())

	#----------------------------------------------------------------------
	def pinActive(self, event=None):
		self._tabFrame["cursor"] = ""
		page = self.getActivePage()
		page.setpin(True)

	#----------------------------------------------------------------------
	# Mouse is clicked => hide ribbon if temporary is active
	#----------------------------------------------------------------------
#	def release(self, event):
#		if self._ribbonFrameStatus == RIBBON_HIDDEN:
#			if not isinstance(event.widget, TabButton): return
#			if time.time() - self.hideTime > 1.0:
#				self._ribbonFrameStatus = RIBBON_TEMP
#				self._allRibbons.pack(side=TOP, fill=BOTH)
#				self.hideTime = time.time()
#
#		elif self._ribbonFrameStatus == RIBBON_TEMP and \
#		     self._allRibbons.winfo_ismapped():
#			if isinstance(event.widget, TabButton): return
#			self._ribbonFrameStatus = RIBBON_HIDDEN
#			self._allRibbons.pack_forget()

	#----------------------------------------------------------------------
	# Double click hide/show ribbon
	#----------------------------------------------------------------------
	def double(self, event):
		if self._ribbonFrameStatus == RIBBON_HIDDEN:
			self._ribbonFrameStatus = RIBBON_SHOWN
			self._pageFrame.pack_forget()
			self._allRibbons.pack(side=TOP, fill=BOTH)
			self._pageFrame.pack(fill=BOTH, expand=YES)
		else:
			self._ribbonFrameStatus = RIBBON_HIDDEN
			self._allRibbons.pack_forget()
			self.hideTime = time.time()

	# ----------------------------------------------------------------------
	def tabFocus(self, event=None):
		self.getActivePage().keyboardFocus()

	#----------------------------------------------------------------------
	# Give focus to the tab on the left
	#----------------------------------------------------------------------
	def _tabLeft(self, event=None):
		slaves = self._tabFrame.pack_slaves()
		try:
			pos = slaves.index(event.widget)-1
		except ValueError:
			if event.widget is self.dynamic:
				pos = len(slaves)-1
			else:
				return
		if pos < 0: return	# Do not replace First tab
		slaves[pos].select()
		#self.changePage()
		slaves[pos].focus_set()

	#----------------------------------------------------------------------
	# Give focus to the tab on the right
	#----------------------------------------------------------------------
	def _tabRight(self, event=None):
		slaves = self._tabFrame.pack_slaves()
		try:
			pos = slaves.index(event.widget)+1
		except ValueError:
			return
		if pos < len(slaves):
			slaves[pos].select()
			#self.changePage()
			slaves[pos].focus_set()
		else:
			# Open dynamic menu
			self.dynamic.select()
			self.dynamic.focus_set()
			self.dynamicMenu()

	#----------------------------------------------------------------------
	def _tabDown(self, event=None):
		self.getActivePage()._ribbonDown(event)

	#----------------------------------------------------------------------
	def tabNext(self, event=None):
		self._tabRight(event)
		self.changePage()
		return "break"

	#----------------------------------------------------------------------
	def tabPrev(self, event=None):
		self._tabLeft(event)
		self.changePage()
		return "break"

	# ----------------------------------------------------------------------
	# Activate next window
	# ----------------------------------------------------------------------
	def nextWindow(self, event):
		active = self.getActivePage()
		found = False
		for page in list(self.pages.values()):
			if page is active:
				found = True

			elif found:
				if page.page is not None and \
				   page.isdocked() and \
				   self._pageFrame.node(page.page):
					self.changePage(page)
					return "break"

		# Start again
		for page in list(self.pages.values()):
			if page.page is not None and \
			   page.isdocked() and \
			   self._pageFrame.node(page.page):
				self.changePage(page)
				return "break"

	# ----------------------------------------------------------------------
	# Activate prev window
	# ----------------------------------------------------------------------
	def prevWindow(self, event):
		# doesn't work
		pass

	#----------------------------------------------------------------------
	# Re-order tabs
	# _dragging flag
	#     0     not dragging
	#     1     dragging
	#     2     finished
	#----------------------------------------------------------------------
	def dragStart(self, event):
		self._dragging = 0
		#self._dragtop  = None

	#----------------------------------------------------------------------
	def drag(self, event):
		if not self._dragging:
		#if not self._tabFrame["cursor"]:
			self._tabFrame.config(cursor="hand1")
			self.dynamic.config(cursor="hand1")
			self._dragging = 1
			#self._dragtop = Toplevel(self)
			#self._dragtop.overrideredirect(1)
			#Label(self._dragtop, text="Bingo", image=tkFlair.icons["save"], compound=LEFT).pack()
			#self._dragtop.geometry("+%d+%d"%(event.x_root, event.y_root))
			#self._dragtop.wait_visibility()
			#self._dragtop.lift()
		#elif self._dragtop is not None:
		#	self._dragtop.geometry("+%d+%d"%(event.x_root, event.y_root))

		if self._dragging != 1: return

		if event.y<-32 or event.y > event.widget.winfo_height()+32:
			try:
				pos = self._tabFrame.pack_slaves().index(event.widget)
			except ValueError:
				if event.widget is self.dynamic:
					pos = -1
			self.dragStop(event)
			self._dragging = 2
			if pos==0:
				# Open a new window
				self.event_generate("<<NewWindow>>")
				return
			self.undock(event.widget["text"], event.x_root, event.y_root)

		elif event.x < -32:
			if event.widget is self.dynamic:
				if self._dragging:
					self.dragStop(event)
					self._dragging = 2
					self.dynamic2tab()
				return
			slaves = self._tabFrame.pack_slaves()
			try: pos = slaves.index(event.widget)
			except ValueError: return
			if pos <= 1: return	# Do not replace First tab
			event.widget.pack(side=LEFT,before=slaves[pos-1])

		elif event.x > event.widget.winfo_width() + 32:
			if event.widget is self.dynamic: return
			slaves = self._tabFrame.pack_slaves()
			try: pos = slaves.index(event.widget)
			except ValueError: return
			if pos==0: return	# Do not replace First tab
			if pos >= len(slaves)-1:
				if self._dragging:
					self.dragStop(event)
					self._dragging = 2
					self.tab2dynamic(event.widget)
				return
			event.widget.pack(side=LEFT,after=slaves[pos+1])

	#----------------------------------------------------------------------
	def dragStop(self, event):
		self._tabFrame.config(cursor="")
		self.dynamic.config(cursor="")
		self._dragging = 0
		#if self._dragtop is not None:
		#	self._dragtop.destroy()
		#self._dragtop = None

	#----------------------------------------------------------------------
#	def configure(self, event):
#		return
#		if not self.winfo_ismapped(): return
#		try:
#			active = self.pages[self.activePage.get()]
#		except KeyError:
#			return
#
#		activex = active._tab.winfo_x()
#		width = self._tabFrame.winfo_width()
#
#		# Pack text
#		if self._tabFrame.winfo_reqwidth() > width:
#			# Find the farthest unpacked and pack it
#			maxx    = 0
#			maxpage = None
#			for page in self.pages.values():
#				if page is active: continue
#				if page.packed: continue
#				dx = abs(page._tab.winfo_x() - activex)
#				if dx>maxx:
#					maxx    = dx
#					maxpage = page
#			if maxpage:
#				maxpage.packed = True
#				maxpage._tab.config(text="")
#				maxpage._tab.config(width=maxpage._tab.winfo_reqwidth())
#
#		else:
#			# Find the closest packed and expand it
#			minx    = 10000
#			minpage = None
#			for page in self.pages.values():
#				if page is active: continue
#				if not page.packed: continue
#				dx = abs(page._tab.winfo_x() - activex)
#				if dx<minx:
#					minx    = dx
#					minpage = page
#			if minpage:
#				minpage.packed = False
#				minpage._tab.config(text=minpage.name)


	#----------------------------------------------------------------------
	# Get/Set pages according to the RPN expression of the tree
	#----------------------------------------------------------------------
	def RPN(self, rpn=None):
		if rpn is None:
			rpn = []
			for item in self._pageFrame.RPN():
				if isinstance(item,Widget):
					# Find page name
					for page in list(self.pages.values()):
						if page.page is item:
							rpn.append(page.name)
							break
				else:
					rpn.append(item)
			return rpn
		else:
			active = None
			for i,item in enumerate(rpn):
				page = self.pages.get(item)
				if page is not None:
					page.create(self)
					rpn[i] = page.page
					if active is None: active = page.name
				else:
					try:
						rpn[i] = float(item)
					except ValueError:
						pass
			self._pageFrame.RPN(rpn)
			self._pageFrame.update_idletasks()
			return active

	#----------------------------------------------------------------------
	# Get/Set tab order
	#----------------------------------------------------------------------
	def tabs(self, tabs=None):
		if tabs is None:
			# return tab order
			tabs = []
			for child in self._tabFrame.pack_slaves():
				for page in list(self.pages.values()):
					if page._tab is child:
						tabs.append(page.name)
			return tabs
		else:
			# Set tab order
			self.activePage.set(tabs[0])
			self.defaultPage = self.pages[tabs[0]]
			for item in tabs:
				page = self.pages.get(item)
				if page is not None:
					page._showTab = True
					self._addTab(page)

	#----------------------------------------------------------------------
	# Print a human readable page tree
	#----------------------------------------------------------------------
	def _printTree(self):
		self.__printNode(self._pageFrame.tree, 0)

	# ----------------------------------------------------------------------
	def __printNode(self, node, depth):
		if node is None: return

		self.__printNode(node.left, depth+1)
		if node.child:
			for p in list(self.pages.values()):
				if p.page is node.child:
					say("   "*depth, p.name, p.pinned)
					break
		else:
			say("   "*depth, " ======== H=",node.hori," pos=",node.pos)
		self.__printNode(node.right, depth+1)
