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

import pdb
import traceback

import os
import sys
import binascii
import io
import tkinter as tk

#sys.path.append(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'lib'))

#if sys.version_info < (2,4) or sys.version > (2,10):
#	sys.stderr.write("FATAL ERROR: flair only works with only with python V2 from 2.4+\n")
#	sys.exit(0)

PRGPATH=os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(PRGPATH, 'lib'))

import re
import bz2
import time
import getopt
import string
import signal
try:
	import configparser
except ImportError:
	import configparser as ConfigParser
try:
	import pickle as pickle
except ImportError:
	import pickle
try:
	from io import StringIO
except ImportError:
	from io import StringIO

try:
	import tkinter
# using the "as" the CallWrapper setting doesn't work
	from tkinter import *
	import tkinter.messagebox as messagebox
except ImportError:
	import tkinter
	from tkinter import *
	import tkinter.messagebox as messagebox
	Tkinter=tkinter

from stat import *

import log
import undo
import tkTree
import tkFlair
import tkExtra

import tkDialogs
import bFileDialog

import Input
# Fast scan for developer option and initialize it if necessary
if "-d" in sys.argv:   Input._developer = True
elif "-D" in sys.argv: Input._developer = False

import Manual
import Layout
import Ribbon
import Project
import Updates
import FlairRibbon

# Plot Engines
import GPPlot
import MPPlot
import Gnuplot

# Importers/Exporters
import Gdml
import Mcnp
import Povray
import OpenSCAD

# Frames
import Output
import RunPage
import PetPage
import PlotPage
import DicomPage
import InputPage
import Materials
import Calculator
import ViewerPage
import CompilePage
import ProjectPage
import PeriodicTable
import GeometryEditor

import EmptyPage

REFRESH_TIME   = 5000
PLOT_ENGINES = ["Gnuplot", "Matplotlib"]

#===============================================================================
root           = None
checkFlukaExe  = True

#-------------------------------------------------------------------------------
def malloc_trim():
	return
	import ctypes
	#ctypes.CDLL('libc.so.6').malloc_trim(0)
	ctypes.CDLL('libc.so').malloc_trim(0)

#-------------------------------------------------------------------------------
def garbage_collect(event=None):
	import gc
	gc.set_debug(gc.DEBUG_STATS | gc.DEBUG_COLLECTABLE | gc.DEBUG_UNCOLLECTABLE \
		| gc.DEBUG_LEAK )
	tkFlair.write("gc isenabled=%d\n"%(gc.isenabled()))
	tkFlair.write("gc threshold=%s\n"%(str(gc.get_threshold())))
	#tkFlair.write("gc count=%s\n"%(str(gc.get_count())))
	tkFlair.write("gc collect=%d"%(gc.collect()))

	try:
		import objgraph
		objgraph.show_most_common_types(limit=20)
	except:
		pass
	malloc_trim()

# XXX Maybe combined with Geometry Viewer? But if module do not exists?
##		( "Debug"     , "debug"     , Ribbon.FlairPage, False),
#===============================================================================
_PAGES =  [	#  Class
		Calculator.CalculatorPage,
		CompilePage.CompilePage  ,
		DicomPage.DicomPage      ,
		ProjectPage.ProjectPage  ,
		InputPage.InputPage      ,
		Materials.MaterialPage   ,
		Output.OutputPage        ,
		PetPage.PetPage          ,
		PlotPage.PlotListPage    ,
		RunPage.RunPage          ,
		PeriodicTable.PeriodicTablePage,
		ViewerPage.ViewerPage
#		EmptyPage.EmptyPage
	]
if GeometryEditor.installed():
	_PAGES.append(GeometryEditor.GeometryEditor)
else:
	_PAGES.append(GeometryEditor.GeometryEditorEmpty)

## FIXME add navigation between pages with Control-Arrow? (oops no InputPage needs them for Up/Down)

#===============================================================================
# FlairTabRibbonFrame
#===============================================================================
class FlairTabRibbonFrame(Ribbon.TabRibbonFrame):
	# ----------------------------------------------------------------------
	def __init__(self, master, mastertab=None, **kw):
		Ribbon.TabRibbonFrame.__init__(self, master, mastertab, **kw)

		self.createNotify()

		# Create bindings for events that can be of general interest
		# in the whole application
		# All private events are handled by their respective page
		toplevel = self.winfo_toplevel()

		# Event generators
		toplevel.event_add("<<Help>>",		"<F1>")
		toplevel.event_add("<<Rename>>",	"<F2>")

		toplevel.event_add("<<Copy>>",		"<Control-c>")
		toplevel.event_add("<<Copy>>",		"<Control-C>")
		toplevel.event_add("<<Clone>>",		"<Control-d>")
		toplevel.event_add("<<Clone>>",		"<Control-D>")

		toplevel.bind("<Control-Key-B>",	garbage_collect)
		toplevel.bind("<Control-Key-q>",	programExit)
		toplevel.bind("<Control-Key-Q>",	programExit)

		toplevel.bind("<Control-Key-w>",	self.close)
		toplevel.bind("<Control-Key-W>",	self.close)
		toplevel.bind("<Control-Key-Return>",	self.execute)

		toplevel.bind("<Control-Tab>",		self.nextWindow)
#		toplevel.bind("<Shift-Control-Tab>",	self.prevWindow)
#		toplevel.bind("<Control-Tab>",		self.tabNext)
#		toplevel.bind("<Control-Shift-Tab>",	self.tabPrev)

		toplevel.bind("<F10>",			self.tabFocus)
		toplevel.bind("<Menu>",			self.tabFocus)
		toplevel.bind("<F12>",			self.minmax)

		toplevel.bind("<Control-r>",		self.refreshKey)
		toplevel.bind("<Control-R>",		self.refreshKey)

		toplevel.event_add("<<Save>>",		"<Control-s>")
		toplevel.event_add("<<Save>>",		"<Control-S>")

		toplevel.event_add("<<Paste>>",		"<Control-v>")
		toplevel.event_add("<<Paste>>",		"<Control-V>")

		toplevel.event_add("<<Cut>>",		"<Control-x>")
		toplevel.event_add("<<Cut>>",		"<Control-X>")
		# Correct Ctrl-y behavior
		toplevel.event_delete("<<Paste>>",	"<Control-y>")
		toplevel.event_delete("<<Paste>>",	"<Control-Y>")
		toplevel.event_add("<<Redo>>",		"<Control-y>")
		toplevel.event_add("<<Redo>>",		"<Control-Y>")

		toplevel.event_add("<<Undo>>",		"<Control-z>")
		toplevel.event_add("<<Undo>>",		"<Control-Z>")

		# Event bindings
		toplevel.bind("<<Help>>",	self.help)
		toplevel.bind("<<Cut>>",	self.cut)
		toplevel.bind("<<Copy>>",	self.copy)
		toplevel.bind("<<Clone>>",	self.clone)
		toplevel.bind("<<Paste>>",	self.paste)
		toplevel.bind("<<Undo>>",	self.undo)
		toplevel.bind("<<Redo>>",	self.redo)
		toplevel.bind("<<Save>>",	self.save)
		toplevel.bind("<<NewWindow>>",	newWindow)

	# ----------------------------------------------------------------------
	# Create pages on initialization
	# ----------------------------------------------------------------------
	def create(self, flair):
		# Create pages
		#for name, iconname, PageClass in _PAGES:
		for PageClass in _PAGES:
			name = PageClass._name_
			iconname = PageClass._icon_
			if isinstance(iconname,str):
				icon = tkFlair.icons[iconname]
			else:
				icon = iconname
			page = PageClass(flair, name, icon, False)
			# page settings
			try:
				settings = tkFlair.config.get(tkFlair._PAGE_SECTION,name).split()
				page.width  = int(settings[0])
				page.height = int(settings[1])
				page.pinned = bool(int(settings[2]))
				# shortcut?
				# undocked?
			except: # ConfigParser.NoOptionError:
				pass
			self.addPage(page)

		# Set tab order
		try:
			tabs = tkFlair.config.get(tkFlair._PAGE_SECTION,"tabs").split()
		except: # ConfigParser.NoOptionError:
			tabs = []
		if not tabs:
			tabs = ["Flair", "Input", "Geometry", "Run", "Data", "Plot"]
		self.tabs(tabs)

		# Remember setup
		self.changeDynamic()	#<--- must ensure that is not shown
		active = None
		try:
			rpn = tkFlair.config.get(tkFlair._PAGE_SECTION,"pages").split()
		except configparser.NoOptionError:
			rpn = []

		if rpn:
			# Checked for errors
			for page in list(self.pages.values()):
				if page.pinned and page.name not in rpn:
					#log.say("Error: page",page.name,"was pinned and not in the rpn")
					#log.say("RPN:", rpn)
					page.pinned = False

			# Count doubles
			for item in rpn:
				try:
					i = float(item)
				except:
					if rpn.count(item) > 1:
						#log.say("Error: duplicate page",item," found in rpn")
						#log.say("RPN:", rpn)
						rpn = []
						break

		# Apply rpn
		if rpn: active = self.RPN(rpn)

		# If project file exist then go to the right page else to Flair
		self.changePage(flair.project.page)	#<--- must ensure that exists

		# Undock pages
		try:
			for name in tkFlair.config.get(tkFlair._PAGE_SECTION,"undocked").split():
				self.undock(name)
		except configparser.NoOptionError:
			pass

	# ----------------------------------------------------------------------
	# retrieve page settings
	# ----------------------------------------------------------------------
	def saveSettings(self):
		# Find Undocked
		undocked = []
		for page in list(self.pages.values()):
			if page.isundocked():
				undocked.append(page.name)
		tkFlair.config.set(tkFlair._PAGE_SECTION, "undocked"," ".join(undocked))

		# Remember page setup
		rpn = self.RPN()
		tkFlair.config.set(tkFlair._PAGE_SECTION, "pages"," ".join(map(str,rpn)))

		# dock all undocked windows
		self.dockAll()

		# Tabs
		tabs = self.tabs()
		tkFlair.config.set(tkFlair._PAGE_SECTION, "tabs"," ".join(tabs))

		# Write individual settings
		for page in list(self.pages.values()):
			tkFlair.config.set(tkFlair._PAGE_SECTION, page.name,
				"%d %d %d"%(page.width, page.height, page.pinned))
			page.configSave()

	# ----------------------------------------------------------------------
	def flair(self):
		return self.getActivePage().flair

	# ----------------------------------------------------------------------
	def close(self, event=None):
		page = self.getActivePage()
		if page.isundocked():
			page.dock()
		else:
			page.flair.close()

	# ----------------------------------------------------------------------
	def save(self, event=None):
		self.flair().saveProject()

	# ----------------------------------------------------------------------
	def cut(self,   event=None):	self.getActivePage().cut(event)
	def copy(self,  event=None):	self.getActivePage().copy(event)
	def paste(self, event=None):	self.getActivePage().paste(event)
	def clone(self, event=None):	self.getActivePage().clone(event)
	def help(self,  event=None):	self.getActivePage().help()
	def minmax(self, event=None):	self.getActivePage().minmaxButton()

	# ----------------------------------------------------------------------
	# undo command
	# ----------------------------------------------------------------------
	def undo(self, event=None):
		page = self.getActivePage()
		page.undo(event)
		page.refreshUndoButton()

	# ----------------------------------------------------------------------
	# redo command
	# ----------------------------------------------------------------------
	def redo(self, event=None):
		page = self.getActivePage()
		page.redo(event)
		page.refreshUndoButton()

	# ----------------------------------------------------------------------
	def undolist(self, event=None):
		if tkFlair.UndoListbox.posted is not None:
			tkFlair.UndoListbox.posted.cancel()
		else:
			b = self.tool["undolist"]
			x = b.winfo_rootx()
			y = b.winfo_rooty() + b.winfo_height()
			dlg = tkFlair.UndoListbox(self, self.getActivePage().flair, x, y)

	# ----------------------------------------------------------------------
	def refreshUndoButton(self):
		page = self.getActivePage()
		page.refreshUndoButton()

	# ----------------------------------------------------------------------
	# refresh a page on specific events
	# If what is None then refresh all viewports
	# otherwise check if the what is in the refreshTrigger list of
	# each page
	#			therefore only the pages dealing with the
	#			executable will be refreshed
	# ----------------------------------------------------------------------
	def refresh(self, what=None, ignore=None):
		#log.say("Tabs.refresh=",what)
		for page in list(self.pages.values()):
			if page is ignore: continue
			if what is None or what in page._refreshOn:
				#log.say("page-invalid:",page.name)
				page.setInvalid()

	# ----------------------------------------------------------------------
	def refreshKey(self, event=None):
		self.refresh()

	# ----------------------------------------------------------------------
	# redraw all frames
	# ----------------------------------------------------------------------
	def redraw(self, event=None):
		for page in list(self.pages.values()):
			if page.created():
				page.redraw()
		self.refreshUndoButton()

	# ----------------------------------------------------------------------
	# Reset settings on each page (before loading new project)
	# ----------------------------------------------------------------------
	def reset(self):
		for page in list(self.pages.values()):
			page.reset()

	# ----------------------------------------------------------------------
	# execute default command in current page
	# ----------------------------------------------------------------------
	def execute(self, event=None):
		self.getActivePage().execute()
		return "break"

	# ----------------------------------------------------------------------
	# active frame update
	# ----------------------------------------------------------------------
	def get(self):
		for page in list(self.pages.values()):
			if page.ismapped():
				page.get()

	# ----------------------------------------------------------------------
	# Dock all undocked windows
	# ----------------------------------------------------------------------
	def dockAll(self):
		for page in list(self.pages.values()):
			if page.page and page.isundocked():
				page.dock()

	#-----------------------------------------------------------------------
	# Create notification frame
	#-----------------------------------------------------------------------
	def createNotify(self):
		self.notifyFrame = Frame(self, relief=GROOVE, background="#FFA0A0")

		self._notifyTitle = Label(self.notifyFrame, relief=RAISED, font="bold")
		self._notifyTitle.grid(row=0, column=0, columnspan=2, sticky=NSEW)

		self._notifyIcon = Label(self.notifyFrame)
		self._notifyIcon.grid(row=1, column=0, sticky=NSEW)

		self._notifyMsg = Message(self.notifyFrame, aspect=500)
		self._notifyMsg.grid(row=1, column=1, sticky=NSEW)

		self._notifyAction = Button(self.notifyFrame)

		self.notifyFrame.grid_rowconfigure(1, weight=1)
		self.notifyFrame.grid_columnconfigure(1, weight=1)

		self._notifyTitle.bind("<Button-1>", self.notifyHide)
		self._notifyTitle.bind("<Enter>",    self._notifyShow)
		self._notifyMsg.bind("<Button-1>",   self.notifyHide)
		self._notifyMsg.bind("<Enter>",      self._notifyShow)
		self._notifyId = None
		self._notifyY = 0

	#-----------------------------------------------------------------------
	def notify(self, title, msg, level=0, action=None, actionCmd=None, timeout=15000):
		self._notifyTitle.config(text=title)
		self._notifyMsg.config(text=str(msg))
		if action:
			self._notifyAction.config(text=action, command=actionCmd)
			self._notifyAction.grid(row=2, column=0, columnspan=2)
		else:
			self._notifyAction.grid_forget()

		color = ("LightGreen", "#FFDD70", "#FFA0A0")[level]
		icon  = ("info32", "warning32", "error32")[level]
		self.notifyFrame.config(background=color)
		self._notifyTitle.config(background=color)
		self._notifyIcon.config(background=color, image=tkFlair.icons[icon])
		self._notifyMsg.config(background=color)
		self._notifyAction.config(background=color, activebackground=color)

		# Increase timeout by the size of the string
		timeout += (len(str(msg))*1000)//20

		if self._notifyId:
			self.notifyFrame.after_cancel(self._notifyId)

		self._notifyShow(timeout=timeout)

	#----------------------------------------------------------------------
	def _notifyShow(self, event=None, timeout=2000):
		if self._notifyId: self.notifyFrame.after_cancel(self._notifyId)
		self.update_idletasks()
		width  = max(300,self.notifyFrame.winfo_reqwidth())
		height = max(100,self.notifyFrame.winfo_reqheight())
		self.notifyFrame.place(x=-width, width=width, relx=1.0,
				y=-height-30, height=height, rely=1.0)
		self._notifyY = -height
		self._notifyId = self.notifyFrame.after(timeout, self._notifyAfter)

	#----------------------------------------------------------------------
	def _notifyAfter(self):
		if self._notifyY < 0:
			self._notifyY += 1
			self.notifyFrame.place(y=self._notifyY-30)
			self._notifyId = self.notifyFrame.after(10, self._notifyAfter)
		else:
			self._notifyId = None
			self.notifyHide()

	#----------------------------------------------------------------------
	def notifyHide(self, event=None):
		if self._notifyId: self.notifyFrame.after_cancel(self._notifyId)
		self.notifyFrame.place_forget()
		self._notifyY =  0
		self._notifyId = None

	#----------------------------------------------------------------------
	def notifyIsVisible(self):
		return self._notifyId is not None

#===============================================================================
# Flair Project Frame
#===============================================================================
class Flair(Toplevel):
	#----------------------------------------------------------------------
	# Initialize Flair
	#----------------------------------------------------------------------
	def __init__(self, master=None, project=None, view=None):
		Toplevel.__init__(self, master, class_="Flair")

		try: geom = "%sx%s" % (tkFlair.config.get(tkFlair.__name__, "width"),
			tkFlair.config.get(tkFlair.__name__, "height"))
		except: geom = "900x700"
		try: self.geometry(geom)
		except: pass

		try:
			self.tk.call('wm','iconphoto',self._w, tkFlair.icons["project"])
		except TclError:
			pass
		self.iconbitmap("@%s/icons/marmite.xbm"%(tkFlair.prgDir))

		self.log          = None
		self.geoedit      = None
		self.plotEngine   = None
		self._plotEngines = {}
		self.undoredo     = undo.UndoRedo()
		self._renumber    = True
		self._renumberPos = sys.maxsize
		self._autosave    = time.time()
		self._afterId     = None

		# Attach a project file
		if project is None:
			self.project = Project.Project()
		else:
			self.project = project
			self._scanProjectDirs()

		self._layout = Layout.LayoutManager(self)

		# Custom Tab
		self.tabs = FlairTabRibbonFrame(self)
		self.tabs.create(self)
		self.tabs.pack(side=TOP, fill=BOTH, expand=YES)

		# Create frame
		self.refreshTitle()

		self.inputPage  = self.tabs.pages["Input"]
		self.viewerPage = self.tabs.pages["Viewer"]
		self.geoedit    = self.tabs.pages["Geometry"]
		self.outputPage = self.tabs.pages["Output"]

		# Check if there are messages in the output
		err = tkFlair._output != ""
		self.log = self.newLog("flair","flair")
		log.set(self.log)

		self.jobs = []	# list of active jobs
		self._afterJob = None

		if view:
			self.tabs.changePage("Viewer")
			for fn in view:
				self.viewerPage.load(fn)

		self.flush()
		self.tabs.refresh()

		self.protocol("WM_DELETE_WINDOW", self.close)
		self.bind("<<Output>>",           self.showOutput)
		self.bind("<<Preferences>>",      self.preferences)

		# focusIn semaphore
		self._inFocus = False
		self.bind("<FocusIn>", self.focusIn)

		if err: self.errorsDuringLoading()

	# --------------------------------------------------------------------
	def notify(self, title, msg, level=0, action=None, actionCmd=None, timeout=10000):
		self.tabs.notify(title, msg, level, action, actionCmd, timeout)
		#self.write("\n%s: %s\n%s"%(("Info", "Warning", "Error")[level], title,msg))

	# --------------------------------------------------------------------
	def notifyHide(self):      self.tabs.notifyHide()
	def notifyIsVisible(self): self.tabs.notifyIsVisible()

	# --------------------------------------------------------------------
	# Add a job to the watch list
	# --------------------------------------------------------------------
	def addJob(self, process, log):
		self.jobs.append((process,log))
		process.start()
		if self._afterJob is None:
			self._afterJob = self.after(1000, self._watchJobs)

	# --------------------------------------------------------------------
	def _watchJobs(self):
		try: os.waitpid(-1,os.WNOHANG)	# Cleanup child processes
		except OSError: pass

		alive = False
		for i in range(len(self.jobs)-1, -1, -1):
			process, log = self.jobs[i]
			if process.isAlive():
				log.write(process.output())
				log.status(process.status, process.percent)
				self.refreshStatus()
				alive = True
			else:
				log.write(process.output())
				log.status(process.status, process.percent)
				self.refreshStatus()
				log.timestamp("Ended:")
				if process.message:
					log.write("\n%s\n%s\n"%tuple(process.message))
					if process.message:
						if process.rc:
							self.notify(process.message[0],
								process.message[1],
								tkFlair.NOTIFY_ERROR,
								"Output",
								self.showOutput)
						else:
							self.notify(process.message[0],
								process.message[1])
				del self.jobs[i]
		if alive:
			self._afterJob = self.after(1000, self._watchJobs)
		else:
			self._afterJob = None

	#----------------------------------------------------------------------
	# Create a new log entry
	#----------------------------------------------------------------------
	def newLog(self, typ, name, append=False):
		if self.outputPage.page is None:
			self.outputPage.create(self.tabs)
		newlog = self.outputPage.newLog(typ, name)
		log.set(newlog)
		return newlog

	#----------------------------------------------------------------------
	# Restore log to flair
	#----------------------------------------------------------------------
	def restoreLog(self):
		log.set(self.log)

	#-----------------------------------------------------------------------
	# Write a string to log entry
	#-----------------------------------------------------------------------
	def write(self, msg):
		if self.log:
			self.flush()
			self.log.write(msg)

	# ----------------------------------------------------------------------
	def flush(self):
		if self.log and tkFlair._output != "":
			self.log.write(tkFlair._output)
			tkFlair._output = ""

	#-----------------------------------------------------------------------
	def setModified(self, m=True):
		self.project.setModified(m)
		self.refreshTitle()
		self.tabs.getActivePage().refreshUndoButton()

	#-----------------------------------------------------------------------
	def setInputModified(self, m=True):
		self.project.setInputModified(m)
		self.refreshTitle()
		self.refreshStatus()
		self.tabs.getActivePage().refreshUndoButton()
		if self._afterId is None:
			self._afterId = self.after(REFRESH_TIME, self.refreshLate)

	# ----------------------------------------------------------------------
	# Late refreshing of information after the input was modified
	# ----------------------------------------------------------------------
	def refreshLate(self):
		if self.inputPage.page: self.inputPage.tree.refresh()
		if int(time.time() - self._autosave) > tkFlair._autoSaveTime:
			self._saveInputPickle("flair-autosave.pickle",False)
			self._autosave = time.time()
		self._afterId = None

	#-----------------------------------------------------------------------
	# Refresh frame title
	#-----------------------------------------------------------------------
	def refreshTitle(self):
		if self.project.isModified():
			c = "+ "
		else:
			c = ""

		if self.project.name == "":
			title = "[untitled]"
		else:
			title = self.project.name

		try: self.title("%s%s - %s" % (c, title, tkFlair.__name__))
		except TclError: pass

	# ----------------------------------------------------------------------
	# Refresh status bar
	# ----------------------------------------------------------------------
	def refreshStatus(self):
		if FlairRibbon.FlairPage.active:
			FlairRibbon.FlairPage.active.refreshStatus()

	# ----------------------------------------------------------------------
	# Set status message
	# ----------------------------------------------------------------------
	def setStatus(self, txt, color="#000000"):
		if FlairRibbon.FlairPage.active:
			FlairRibbon.FlairPage.active.setStatus(txt,color)

	# ----------------------------------------------------------------------
	def refreshGeoedit(self):
		if self.geoedit.created():
			self.geoedit.loadProject()

	# ----------------------------------------------------------------------
	def resetUndo(self):
		# Reset undo info on all frames
		for page in list(self.tabs.pages.values()):
			try:
				page.resetUndo()
				page.refreshUndoButton()
			except AttributeError:
				pass

	# ======================================================================
	# Event handlers
	# ======================================================================
	# Focus In event
	# ----------------------------------------------------------------------
	def focusIn(self, event):
		if self._inFocus: return
		# FocusIn is generated for all sub-windows, handle only the main window
		if self is not event.widget: return
		self._inFocus = True

		try:	# FIXME I should use absolute paths in the future
			# Always change the directory to the project one
			os.chdir(self.project.dir)
		except OSError:
			pass
		self.checkInputFile()
		self._inFocus = False

	# ----------------------------------------------------------------------
	def checkInputFile(self, parent=None):
		if self.project.input.checkInputFile():
			if self.project.isInputModified():
				if parent is None: parent = self
				if tkFlair.askyesno("Warning",
				    "Input file %s has changed since editing started\n"
				    "Reload new version" % (self.project.inputFile),
				    parent=parent):
					self.project.loadInput()
					self.refresh()
				else:
					self.setInputModified()
			else:
				self.project.loadInput()
				self.refresh()

	# ----------------------------------------------------------------------
	# Project Files
	# ----------------------------------------------------------------------
	def newProject(self, event=None):
		self.get()
		if self.project.isModified():
			if not tkFlair.askyesno("New Project",
					"Current Project is modified\nDiscard changes",
					parent=self):
				return

		self.reset()
		self.project.clear()
		self.project.input.minimumInput()
		self.tabs.refresh()
		self.refreshStatus()
		self.refreshTitle()
		self.resetUndo()

	# ----------------------------------------------------------------------
	# Scan all project directories and append them to the bFileDialog history
	# ----------------------------------------------------------------------
	def _scanProjectDirs(self):
		def add(name):
			if not isinstance(name,str): return
			if name == Project.DEFAULT_INPUT: return
			try:
				bFileDialog.append2History(os.path.dirname(os.path.abspath(name)))
			except:
				pass

		add(self.project.exe)
		add(self.project.defmain)

		for fn in self.project.sourceList:
			add(fn)

		for run in self.project.runs:
			add(run.name)
			for usr in run.usrinfo:
				add(usr.name())

		for plot in self.project.plots:
			add(plot.name)
			for name,value in list(plot.var.items()):
				if "file" in name:
					add(value)

		for layer in list(self.project.geometry.layers.values()):
			for name,value in list(layer.var.items()):
				if "file" in name:
					add(value)

	# ----------------------------------------------------------------------
	# Load project
	# ----------------------------------------------------------------------
	def loadProject(self, filename):
		self.logReset()
		self.reset()
		try:
			self.project.clear()
			self.project.load(filename)
		except (IOError, OSError):
			self.notify("Error loading project",
				"Can not load project: \"%s\" either the file " \
				"does not exist or you have insufficient " \
				"permissions to access it"%(filename),
				tkFlair.NOTIFY_ERROR)
			return False

		self._scanProjectDirs()

		self._autosave = time.time()
		self.refreshStatus()
		self.setInputModified(False)
		self.setModified(False)
		self.refreshTitle()
		self.resetUndo()
		addRecent(filename)
		self.tabs.refresh()
		self.geoedit.loadState()
		self.tabs.changePage(self.project.page)	#<--- must ensure that exists

		if not self.errorsDuringLoading():
			self.notify("Loaded project",
				"Project: %s\nInput: %s\nDir: %s" % \
					(os.path.basename(filename),
					 self.project.inputName,
					 self.project.dir),
				timeout=2000)

		return True

	# ----------------------------------------------------------------------
	# Open Project
	# ----------------------------------------------------------------------
	def openProject(self, event=None):
		self.get()
		if self.project.isModified():
			if not tkFlair.askyesno("Open Project",
					"Current Project is modified\nDiscard changes",
					parent=self):
				return
		filename = bFileDialog.askopenfilename(master=self,
			title="Open flair project",
			filetypes=[("Flair Project files","*.flair"),
				("All","*")])
		if len(filename) > 0:
			fn, ext = os.path.splitext(filename)
			if ext == ".inp":
				# Just as a failsafe precaution
				self.project.loadInput(filename)
				self.redraw()
			else:
				if self.loadProject(filename):
					addRecent(self.project.projFile)

	# ----------------------------------------------------------------------
	def _saveProject(self, filename=None):
		try:
			self.project.save(filename)
			return False
		except (IOError, OSError):
			messagebox.showerror("Error saving project file",
				"The following error occured during saving of the project\nfile: %s/%s\nError: %s" \
				%(self.project.dir, self.project.name, sys.exc_info()[1]),
				parent=self)
			return True

	# ----------------------------------------------------------------------
	# Save Project
	# ----------------------------------------------------------------------
	def saveProject(self, event=None):
		"""Save current project"""
		self.get()
		self.project.page = self.tabs.activePage.get()
		if self.project.name == "":
			if self.saveProjectAs(): return True
			if self.project.inputName == "":
				name,ext = os.path.splitext(self.project.name)
				self.project.inputName = name
				self.project.inputFile = name+".inp"
		else:
			if self.project.isInputModified():
				if self.saveInput(): return True
			if self._saveProject(): return True
		self.refreshTitle()
		self.refreshStatus()
		return False

	# ----------------------------------------------------------------------
	# Save Project As
	# ----------------------------------------------------------------------
	def saveProjectAs(self, event=None):
		"""Save current project with different filename"""
		self.get()
		if self.project.inputName == "" and len(self.project.input.cardlist)>0:
			self.refreshTitle()
			if self.project.isInputModified():
				if self.saveInput(): return True
				# Give the same filename to flair project
				fn, ext = os.path.splitext(self.project.inputFile)
				if fn:
					self._saveProject("%s.flair"%(fn))
					addRecent(self.project.projFile)
					return self.saveProject()
				else:
					return False

		self.project.setModified(False)
		if self.project.name == "":
			initfile = os.path.splitext(self.project.inputName)[0]
		else:
			initfile = self.project.name
		filename = bFileDialog.asksaveasfilename(master=self,
			title="Save flair project as",
			initialfile=initfile,
			filetypes=[("Flair Project files",("*.flair")),
				("All","*")])
		if len(filename) > 0:
			# Correct for accidental mistype of extension
			fn, ext = os.path.splitext(filename)
			if ext == ".inp": filename = "%s.flair"%(fn)
			if self._saveProject(filename): return True
			if self.project.isInputModified():
				if self.saveInput(): return True
			#self.tabs.refresh()
			self.refreshTitle()
			addRecent(self.project.projFile)

		self.refreshTitle()
		self.refreshStatus()
		return False

	# ----------------------------------------------------------------------
	def loadInput(self):
		self.get()
		filename = bFileDialog.askopenfilename(master=self,
			title="Load input file",
			filetypes=[("Fluka Input files",("*.inp","*.fluka")),
				("Fluka Input pickle files","*.pickle"),
				("All","*")],
			initialfile=self.project.inputFile)
		if len(filename) > 0:
			self.logReset()
			# Correct for accidental mistype of extension
			fn, ext = os.path.splitext(filename)
			if ext == ".flair": filename = "%s.inp"%(fn)
			if ext == ".pickle":
				self._loadInputPickle(filename)
			else:
				self.project.loadInput(filename)
			self.resetUndo()	# << before the refresh since it deletes all cards
			self.refresh()
			self.errorsDuringLoading()

	# ----------------------------------------------------------------------
	def _loadInputPickle(self, filename=None):
		try:
			f = bz2.BZ2File(filename,'r')
			self.project.input.load(pickle.Unpickler(f))
			f.close()
			return False
		except (IOError, OSError):
			messagebox.showerror("Error loading input pickled file",
				"The following error occured during saving of the input\nfile: %s/%s\nError: %s" \
				%(self.project.dir, self.project.inputFile, sys.exc_info()[1]),
				parent=self)
			return True

	# ----------------------------------------------------------------------
	def _saveInputPickle(self, filename=None, showerror=True):
		try:
			f = bz2.BZ2File(filename,'w')
			self.project.input.dump(pickle.Pickler(f))
			f.close()
			return False
		except (IOError, OSError):
			if showerror:
				messagebox.showerror("Error saving input pickled file",
					"The following error occurred during saving of the input\nfile: %s/%s\nError: %s" \
					%(self.project.dir, self.project.inputFile, sys.exc_info()[1]),
					parent=self)
			return True


	# ----------------------------------------------------------------------
	def _saveInput(self, filename=None):
		try:
			self.project.saveInput(filename)
			return False
		except (IOError, OSError):
			messagebox.showerror("Error saving input file",
				"The following error occurred during saving of the input\nfile: %s/%s\nError: %s" \
				%(self.project.dir, self.project.inputFile, sys.exc_info()[1]),
				parent=self)
			return True

	# ----------------------------------------------------------------------
	def saveInput(self):
		self.get()
		self._autosave = time.time()
		if self.project.inputFile == "":
			return self.saveInputAs()
		else:
			if self.project.isInputModified():
				return self._saveInput()
		return False

	# ----------------------------------------------------------------------
	def saveInputAs(self):
		self.tabs.get()
		self._autosave = time.time()

		initialfile=self.project.inputFile
		if initialfile == "":
			name,ext = os.path.splitext(self.project.name)
			initialfile = name

		filename = bFileDialog.asksaveasfilename(master=self,
			title="Save input file as",
			initialfile=initialfile,
			filetypes=[("Fluka Input files","*.inp"),
				("Fluka Input pickle files","*.pickle"),
				("All","*")])
		if len(filename)==0: return True

		# Correct for accidental mistype of extension
		fn, ext = os.path.splitext(filename)
		if ext == ".flair": filename = "%s.inp"%(fn)
		if ext == ".pickle":
			if self._saveInputPickle(filename):
				return True
		elif self._saveInput(filename): return True
		self.tabs.refresh()

		self.refreshTitle()
		self.refreshStatus()
		return False

	# ----------------------------------------------------------------------
	def openInputTemplate(self, template=None):
		if self.project.isInputModified() and len(self.project.input.cardlist)>0:
				if not tkFlair.askyesno("Exit",
					"Project is modified.\n"
					"Create a new input file?",
					parent=self):
					return
		self.project.inputFile = ""
		self.project.input = Input.Input()
#		self._autosave  = time.time()
		# Read template
		if template is not None:
			self.project.input.read(template)
			self.notify("New Input",
				"New Input file has been created based on" \
				" template: %s"%(os.path.basename(template)))
		else:
			self.project.input.minimumInput()
		#self.resetUndo()
		self.redraw()
		self.refresh()

	# ----------------------------------------------------------------------
	# Import cards from another fluka file
	# ----------------------------------------------------------------------
	def importFluka(self):
		self.notify("Import",
			"Fluka input importing is not implemented for the moment",
			tkFlair.NOTIFY_WARNING)

		#inputFrame = self.frames["input"][0]
		#inputFrame.show()
		#inputFrame.insert2Input()

	# ----------------------------------------------------------------------
	def importGdml(self):
		filename = bFileDialog.askopenfilename(master=self,
			title="Import GDML file",
			filetypes=[("Gdml files","*.gdml"),
				("All","*")])
		if len(filename) == 0: return
		oldSize = self.log.size()	# remember size before loading
		tkFlair.incStats("importGdml")
		relfn = self.project.relativePath(filename)
		gdml = Gdml.Gdml(relfn)
		gdml.toFluka(self.project.input)
		if self.log.size() != oldSize:
			self.notify("GDML Import",
				"Errors or warnings during GDML import.\nOpen the Output window for more information.",
				tkFlair.NOTIFY_ERROR, "Output", self.showOutput)
		else:
			self.notify("GDML Import",
				"GDML file \"%s\" imported"%(relfn))
		self.setInputModified()
		self.redraw()

	# ----------------------------------------------------------------------
	def importMcnp(self):
		filename = bFileDialog.askopenfilename(master=self,
			title="Import MCNP file",
			filetypes=[("Mcnp files",("*.mcnp", "*.m")),
				("All","*")])
		if len(filename) == 0: return
		oldSize = self.log.size()	# remember size before loading
		tkFlair.incStats("importMcnp")
		relfn = self.project.relativePath(filename)
		mcnp = Mcnp.Mcnp()
		mcnp.toFluka(relfn, self.project.input)
		if self.log.size() != oldSize:
			self.notify("MCNP Import",
				"Errors or warnings during MCNP import.\nOpen the Output window for more information.",
				tkFlair.NOTIFY_ERROR, "Output", self.showOutput)
		else:
			self.notify("MCNP Import",
				"MCNP file \"%s\" imported"%(relfn))
		self.setInputModified()
		self.redraw()
		#self.refreshGeoedit()

	# ----------------------------------------------------------------------
	def exportMcnp(self):
		name,ext = os.path.splitext(self.project.name)
		filename = bFileDialog.asksaveasfilename(master=self,
			title="Save Input as MCNP file",
			initialfile = name + ".mcnp",
			filetypes=[("Mcnp files",("*.mcnp", "*.m")),
				("All","*")])
		if len(filename) == 0: return
		oldSize = self.log.size()	# remember size before loading
		tkFlair.incStats("exportMcnp")
		relfn = self.project.relativePath(filename)
		mcnp = Mcnp.Mcnp()
		mcnp.setMacroBodies(tkFlair.askyesno("Export to MCNP",
			"Do you want to use macrobodies in MCNP/X like RPP, BOX, RCC...?",
			parent=self))
		mcnp.fromFluka(self.project.input)
		mcnp.write(relfn)
		if self.log.size() != oldSize:
			self.notify("MCNP Export",
				"Errors or warnings during MCNP export.\nOpen the Output window for more information.",
				tkFlair.NOTIFY_ERROR, "Output", self.showOutput)
		else:
			self.notify("MCNP Export",
				"Input exported to MCNP format as %s"%(relfn))

	# ----------------------------------------------------------------------
	def exportPovray(self):
		name,ext = os.path.splitext(self.project.name)
		filename = bFileDialog.asksaveasfilename(master=self,
			title="Save Input as Povray file",
			initialfile = name + ".pov",
			filetypes=[("Povray files",("*.pov", "*.inc")),
				("All","*")])
		if len(filename) == 0: return
		oldSize = self.log.size()	# remember size before loading
		tkFlair.incStats("writePovray")
		relfn = self.project.relativePath(filename)
		Povray.export(self.project.input, relfn)
		if self.log.size() != oldSize:
			self.notify("Povray Export",
				"Errors or warnings during Povray export.\nOpen the Output window for more information.",
				tkFlair.NOTIFY_ERROR, "Output", self.showOutput)
		else:
			self.notify("Povray Export",
				"Input exported to Povray format as %s"%(relfn))

	# ----------------------------------------------------------------------
	def exportOpenscad(self):
		name,ext = os.path.splitext(self.project.name)
		filename = bFileDialog.asksaveasfilename(master=self,
			title="Save Input as OpenSCAD file",
			initialfile = name + ".scad",
			filetypes=[("Povray files","*.scad"),
				("All","*")])
		if len(filename) == 0: return
		oldSize = self.log.size()	# remember size before loading
		tkFlair.incStats("writeOpenscad")
		relfn = self.project.relativePath(filename)
		OpenSCAD.export(self.project.input, relfn)
		if self.log.size() != oldSize:
			self.notify("OpenSCAD Export",
				"Errors or warnings during OpenSCAD export.\nOpen the Output window for more information.",
				tkFlair.NOTIFY_ERROR, "Output", self.showOutput)
		else:
			self.notify("OpenSCAD Exporting",
				"Input exported to OpenSCAD format as %s"%(relfn))

	# ----------------------------------------------------------------------
	def exportMakefile(self):
		if self.project.name == "":
			initialfile = "Makefile"
		else:
			fn,ext = os.path.splitext(self.project.name)
			initialfile = "%s.make"%(fn)
		tkFlair.incStats("writeMakefile")
		fmake = bFileDialog.asksaveasfilename(master=self,
			title="Save makefile as",
			initialfile=initialfile,
			filetypes=[ ("All","*"),
				("Make files",("*.mak","*.make"))])
		if len(fmake) > 0:
			self.project.writeMakefile(fmake)

	# ----------------------------------------------------------------------
	def logReset(self):
		if self.log: self.log.reset()

	# ----------------------------------------------------------------------
	def errorsDuringLoading(self):
		if self.log and self.log.size():
			self.notify("Warnings during loading",
				"Error or warning messages appeared while loading input. " \
				"Open the Output window for more information.", tkFlair.NOTIFY_ERROR,
				"Output", self.showOutput)
			return True
		else:
			return False

	# ----------------------------------------------------------------------
	# Change current directory
	# ----------------------------------------------------------------------
	def chdir(self, event=None):
		if self.project.inputName != "" or self.project.name != "":
			self.notify("WARNING Change directory", \
				"Changing the current directory while working " \
				"on a project can end-up in broken file links",
				tkFlair.NOTIFY_WARNING)

		d = bFileDialog.askdirectory(master=self,
			filetypes=[("Flair Project files","*.flair"),
				("All","*")])
		if len(d) > 0:
			self.project.chdir(d)
			self.setModified()
			self.setInputModified()
			self.redraw()

	# ----------------------------------------------------------------------
	# redraw all frames
	# ----------------------------------------------------------------------
	def redraw(self, event=None):
		self.tabs.redraw()

	# ----------------------------------------------------------------------
	# refresh
	# ----------------------------------------------------------------------
	def refresh(self, what=None, ignore=None):
		self.tabs.refresh(what, ignore)

	# ----------------------------------------------------------------------
	def reset(self):
		self.tabs.reset()

	# ----------------------------------------------------------------------
	def get(self): self.tabs.get()

	#-------------------------------------------------------------------------------
	# Open files in viewer
	#-------------------------------------------------------------------------------
	def viewer(self, files):
		if len(files) == 0: return
		tkFlair.incStats("viewer")
		for f in files:
			if isinstance(f,tuple):
				# (name,data)
				self.tabs.changePage(self.viewerPage)
				self.viewerPage.show(*f)
			else:
				# filename
				try:
					s = os.stat(f)
				except:
					self.write("Error: Viewer cannot open file \'%s\'\n" %(f))
					continue

				filename = os.path.basename(f)

				# Check for a directory
				if S_ISDIR(s[ST_MODE]):
					continue

				elif self.project is not None and \
						(filename=="core" or Project._PAT_CORE.match(filename)):
					# Check for core dump
					exe = self.project.executable()
					self.write("%s -e %s %s %s&"%(Project.terminal,
						Project.debugger, exe, f))
					os.system("%s -e %s %s %s&"%(Project.terminal,
						Project.debugger, exe, f))
				else:
					self.tabs.changePage(self.viewerPage)
					if self.viewerPage.hasFile(f):
						self.viewerPage.reload(f)
					else:
						self.viewerPage.load(f)

	# ----------------------------------------------------------------------
	def showPlot(self, name):
		# Search plots for the appropriate one
		for plot in self.project.plots:
			if name == plot.name:
				# Switch to Plots page and activate the plot
				plotPage = self.tabs.pages["Plot"]
				self.tabs.changePage(plotPage)
				plotPage.showPlot(name)
				break

	# ----------------------------------------------------------------------
	def showOutput(self, event=None):
		self.notifyHide()
		if self.outputPage.isdocked() and self.tabs.getActivePage() is self.outputPage:
			self.outputPage.hide(True)
		else:
			self.tabs.changePage(self.outputPage)

	# ----------------------------------------------------------------------
	def preferences(self, event=None):
		global checkFlukaExe
		self.notifyHide()
		if tkFlair.preferences(self):
			self.redraw()
			checkFlukaExe = True
			self.checkFluka()

	# ----------------------------------------------------------------------
	# Close Window
	# ----------------------------------------------------------------------
	def close(self, event=None):
		global root

		self.get()
		if self.project.isModified():
			title = self.project.name
			if title=="": title = "[untitled]"
			ans = tkDialogs.Dialog(self,
				{"title":"Close: %s"%(title),
				"text":"Project %s has been modified. "
					"Do you want to save it before "
					"exiting the application?"%(title),
				"bitmap": "questhead",
				"default": 2,
				"strings": ('Save and Close','Discard and Close','Cancel')})
			if ans.num==0:
				self.saveProject()
			elif ans.num in (2,-1):
				return

		# FIXME check if there are unfinished jobs in the queue and ask the user
		# XXX XXX XXX
		# XXX XXX XXX
		# XXX XXX XXX

		# Save page settings
		self.tabs.saveSettings()

		# special I/O
		# splitters
		for name in ("Input", "Geometry", "Plot"):
			page = self.tabs.pages[name]
			if not page.created(): continue
			try:
				tkFlair.config.set(tkFlair._PAGE_SECTION,
						"%s.split"%(name),
						page.splitter.split)
			except AttributeError:
				pass

			# Save viewer splitters
			if name == "Geometry":
				try:
					tkFlair.config.set(tkFlair._PAGE_SECTION,
							"%s.vsplit"%(name),
							page.vsplit())
					tkFlair.config.set(tkFlair._PAGE_SECTION,
							"%s.tsplit"%(name),
							page.tsplit())
					tkFlair.config.set(tkFlair._PAGE_SECTION,
							"%s.bsplit"%(name),
							page.bsplit())
				except AttributeError: pass

		# Save window info on config
		try:
			tkFlair.config.set(tkFlair.__name__, "width",
					self.winfo_width())
			tkFlair.config.set(tkFlair.__name__, "height",
					self.winfo_height())
		except TclError:
			pass

		flairCount = 0
		for child in root.winfo_children():
			if isinstance(child, Flair):
				flairCount += 1

		FlairRibbon.FlairPage.clearActive()
		self.destroy()
		if root is not None:
			if flairCount==1:
				# Before closing the program
				check4Errors()
				root.destroy()
				root = None
				Manual.destroy()

	# ----------------------------------------------------------------------
	# Quit program
	# ----------------------------------------------------------------------
	@staticmethod
	def quit():
		programExit()

	# ----------------------------------------------------------------------
	# Check for updates
	# ----------------------------------------------------------------------
	def checkUpdates(self):
		# Find fluka version
		fluka_version = Updates.flukaVersion(Project.flukaDir)
		Updates.CheckUpdateDialog(self, tkFlair.__version__, fluka_version)

	# ----------------------------------------------------------------------
	# Check fluka installation
	# ----------------------------------------------------------------------
	def checkFluka(self):
		# Verify if flukahp exists in this directory
		global checkFlukaExe
		if not checkFlukaExe: return
		checkFlukaExe = False

		if Project.flukaDir != "":
			exe = os.path.join(Project.flukaDir,Project.flukaExe)
		else:
			exe = "?"
		try:
			os.stat(exe)
			return
		except:
			pass
		self.write("Error: %s do not exist\n"%(exe))
		self.write("Please set the correct FLUKA directory in Preferences")
		self.notify("Error FLUKA executable not found",
			"Cannot find FLUKA executable \"%s\"\n" \
			"Please set the correct FLUKA directory in Preferences"%(exe),
			tkFlair.NOTIFY_ERROR,
			"Preferences", self.preferences)

	# ----------------------------------------------------------------------
	def clipboard(self):
		return self.selection_get(selection='CLIPBOARD', type='ASCII')

	# ----------------------------------------------------------------------
	def clipboardClear(self):
		self.clipboard_clear()

	# ----------------------------------------------------------------------
	# Parse text contents as flair cards normally comes from clipboard
	# ----------------------------------------------------------------------
	def clipboard2Cards(self):
		try: clipboard = self.selection_get(selection='CLIPBOARD',type='STRING')
		except: return

		cards = []

		cardtarget = tkFlair._FLAIRCLIP+"<card>"
		if clipboard.startswith(cardtarget) or True: 
			# Pickler format
#			unpickler = pickle.Unpickler(StringIO(clipboard[len(cardtarget):]))
			dummy_io = io.BytesIO() # creates arficial file
			dummy_io.write(binascii.unhexlify(clipboard)) # writes the decoded pickle files to artificial file
			dummy_io.seek(0) # sets the "reader point to the first values"
			unpickler = pickle.Unpickler(dummy_io) # creates the pickler
			print(unpickler.load()) # These two loads load cardtaget
			print(unpickler.load()) # These two loads load +<card>
			while True:
				try:
					card = Input.Card("error")
					card.load(unpickler)
					cards.append(card)
				except EOFError:
					break

		else:
			# We are in trouble, text format
			# Try to guess the format of the cards

			# Create a temporary input
			inp = Input.Input()

			inp._openFile(StringIO(clipboard))

			# find location and call the appropriate parser
			# FIXME I should scan the lines to find the appropriate parser
			inp.parse()

			# Move cards to clipboard
			#cards = inp.cardlist

			# Remove error cards
			cards = [x for x in inp.cardlist \
					if Input.CardInfo.get(x).name != Input.ERROR]
			del inp

		return cards

	# ======================================================================
	# Undo/redo input functions
	# ======================================================================
	def addUndo(self, undoinfo):
		if not undoinfo:
			return

		if isinstance(undoinfo, list):
			undoinfo = undo.createListUndo(undoinfo)

		if undoinfo is not undo.NullUndo:
			self.undoredo.addUndo(undoinfo)
			self.setInputModified()

	# ----------------------------------------------------------------------
	# Create undo list with refresh and renumber calls
	# ----------------------------------------------------------------------
	def undolistRefresh(self, undoinfo=None):
		if undoinfo is None:
			return undo.UndoInfo([ self.refreshUndo(), self.renumberUndo()])
		else:
			undoinfo.append(self.renumberUndo())
			undoinfo.append(self.refreshUndo())

	# ----------------------------------------------------------------------
	# add undolist with refreshUndo and renumberUndo
	# ----------------------------------------------------------------------
	def addUndolistRefresh(self, undoinfo, msg):
		self.undolistRefresh(undoinfo)
		self.addUndo(undo.createListUndo(undoinfo,msg))

	# ----------------------------------------------------------------------
	def canUndo(self):	return self.undoredo.canUndo()
	def canRedo(self):	return self.undoredo.canRedo()

	# ----------------------------------------------------------------------
	def undoInput(self):
		self.inputPage.commit()
		if self.canUndo():
			self.undoredo.undo()
			self.refresh()
			self.project.input.clearCache("bodies")

	# ----------------------------------------------------------------------
	def redoInput(self):
		self.inputPage.commit()
		if self.canRedo():
			self.undoredo.redo()
			self.refresh()
			self.project.input.clearCache("bodies")

	# ----------------------------------------------------------------------
	# Activate card if needed in input page
	# ----------------------------------------------------------------------
	def activateIfNeededInputFrame(self, card):
		"""activate card if needed in Input Editor"""
		if card.input is None: return
		if self.inputPage.page:
			self.inputPage.canvas.activateIfNeeded(card)

	# ----------------------------------------------------------------------
	# Activate card in InputPage editor
	# ----------------------------------------------------------------------
	def activateIfNeeded(self, card):
		"""activate card if needed in all Input editor and Geometry Viewer"""
		self.activateIfNeededInputFrame(card)
		if self.geoedit.created():
			self.geoedit.refreshUndo()
			self.geoedit.activateIfNeeded(card)
			self.geoedit.redraw()

	# ----------------------------------------------------------------------
	# set what value to card and return undo info
	# ----------------------------------------------------------------------
	def setWhatUndo(self, card, w, v):
		"""set what value and return undo info"""
		undoinfo = ("Card %s:%s what(%d)=%s"%(card.tag, card.sdum(), w, str(v)),
			self.setWhatUndo, card, w, card.what(w))
		card.setWhat(w, v)
		self.setInputModified()
		self.activateIfNeeded(card)
		if card.info.useUnits:
			# Print scan units
			self.project.input.scanUnits()
		return undoinfo

	# ----------------------------------------------------------------------
	def setAbsWhatUndo(self, card, w, v):
		"""set absolute what value and return undo info"""
		undoinfo = ("Card %s:%s abs-what(%d)=%s"%(card.tag, card.sdum(), w, str(v)),
			self.setAbsWhatUndo, card, w, card.absWhat(w))
		card.setAbsWhat(w, v)
		self.setInputModified()
		self.activateIfNeeded(card)
		return undoinfo

	# ----------------------------------------------------------------------
	def setSignUndo(self, card, w, s):
		"""Change the sign of a card.what, and return undo info"""
		undoinfo = ("Card %s:%s sign-what(%d)"%(card.tag, card.sdum(), w),
			self.setSignUndo, card, w, card.sign(w))
		card.setSign(w, s)
		self.setInputModified()
		self.activateIfNeeded(card)
		return undoinfo

	# ----------------------------------------------------------------------
	def setCommentUndo(self, card, comment):
		"""set comment and return undo information"""
		undoinfo = ("Card %s:%s comment"%(card.tag, card.sdum()),
				self.setCommentUndo, card, card.comment())
		card.setComment(comment)
		self.setInputModified()
		self.activateIfNeeded(card)
		return undoinfo

	# ----------------------------------------------------------------------
	def setCardUndo(self, card, name, value, default=None):
		"""set an associated variable var with value and return undo information"""
		undoinfo = ("Card %s:%s variable name=%s"%(card.tag, card.sdum(), name),
				self.setCardUndo, card, name, card.get(name,default), default)
		card[name] = value
		self.setInputModified()
		return undoinfo

	# ----------------------------------------------------------------------
	# Create the undo sequence for modifying a region with multiple cards
	# as well delete the cards and merge the comments
	# ----------------------------------------------------------------------
	def _modifyContRegionZoneUndo(self, region, msg):
		regions = self.project.input.cardsSorted("REGION")
		try:
			pos = regions.index(region)
		except ValueError:
			return undo.NullUndo

		undoinfo = [(self.refreshUndo,), self.renumberUndo()]
		undoinfo.append((self.setWhatUndo, region, -1, region.extra()))
		undoinfo.append((self.setCommentUndo, region, region.comment()))

		delcards = []
		expr     = region.extra()
		comment  = region.comment()
		for card in regions[pos+1:]:
			if card.sdum()=="&":
				cmt = card.comment()
				if cmt != "":
					if comment == "":
						comment = cmt
					else:
						comment += "\n"+cmt
				expr += "\n" + card.extra()
				delcards.append(card)
			else:
				break
		region.setComment(comment)
		region.setExtra(expr)

		for card in reversed(delcards):
			undoinfo.append(self.delCardUndo(card.pos(), False))

		undoinfo.append(self.renumberUndo())
		undoinfo.append((self.refreshUndo,))
		region["@cont"] = False
		return undo.createListUndo(undoinfo, msg)

	# ----------------------------------------------------------------------
	# add zone to region and return undo information
	# ----------------------------------------------------------------------
	def addZoneUndo(self, card, zone):
		msg = "Card %s:%s add zone"%(card.tag, card.sdum())
		if card["@cont"]:
			undoinfo = self._modifyContRegionZoneUndo(card, msg)
			if undoinfo is undo.NullUndo:
				return undoinfo
			self.refresh("card")
		else:
			undoinfo = (msg, self.setWhatUndo, card, -1, card.extra())

		card.addZone(zone)
		self.setInputModified()
		self.activateIfNeeded(card)
		return undoinfo

	# ----------------------------------------------------------------------
	# set zones to region and return undo information
	# ----------------------------------------------------------------------
	def setZonesUndo(self, card, zones):
		msg ="Card %s:%s set zones"%(card.tag, card.sdum())
		if card["@cont"]:
			undoinfo = self._modifyContRegionZoneUndo(card, msg)
			if undoinfo is undo.NullUndo:
				return undoinfo
			self.refresh("card")
		else:
			undoinfo = (msg, self.setWhatUndo, card, -1, card.extra())

		card.setZones(zones)
		self.setInputModified()
		self.activateIfNeeded(card)
		return undoinfo

	# ----------------------------------------------------------------------
	# set the number of card whats
	# ----------------------------------------------------------------------
	def setNWhatsUndo(self, card, n):
		if n > card.nwhats():
			undoinfo = ("Card %s:%s set nwhats=%d"%(card.tag, card.sdum(), n),
				self.setNWhatsUndo, card, card.nwhats())
			card.setNWhats(n)
			self.setInputModified()
			self.activateIfNeeded(card)
			return undoinfo
		elif n < card.nwhats():
			undoinfo = []
			for w in range(n+1, card.nwhats()):
				undoinfo.append(self.setWhatUndo(card, w, ""))
			undoinfo.append((self.setNWhatsUndo, card, card.nwhats()))
			card.setNWhats(n)
			self.setInputModified()
			self.activateIfNeeded(card)
			return undo.createListUndo(undoinfo)
		else:
			#if card is not self.active.card: self.activate(card)
			return undo.NullUndo

	# ----------------------------------------------------------------------
	# Add/Delete a card and return undo info
	# WARNING: For successive add/del cards disable renumbering and call
	#          twice renumberUndo() before and after the commands.
	#          The double call is needed since redo will change the order
	#          of the commands
	# ----------------------------------------------------------------------
	def addCardUndo(self, card, pos, renumber):
		"""add card to input as position pos and return undo info"""
		undoinfo = ("Add card %s:%s"%(card.tag, card.sdum()),
			self.delCardUndo, pos, renumber)
		self.project.input.addCard(card, pos, renumber)
		self.setInputModified()
		if not renumber:
			self._renumberPos = min(self._renumberPos, pos)
		return undoinfo

	# ----------------------------------------------------------------------
	def delCardUndo(self, pos, renumber):
		"""delete a card from input and return unfo info"""
		card  = self.project.input.cardlist[pos]
		undoinfo = ("Delete card %s:%s"%(card.tag, card.sdum()),
			self.addCardUndo, card, pos, renumber)
		self.project.input.delCard(pos, renumber)
		self.setInputModified()
		if not renumber:
			self._renumberPos = min(self._renumberPos, pos)
		return undoinfo

	# ----------------------------------------------------------------------
	# Remove all body references
	# ----------------------------------------------------------------------
	def removeBodyRefUndo(self, body):
		# Search all regions for possibly a reference to this body
		# and remove it
		pat = re.compile(r"([\+\-]\b%s\b)" % (re.escape(body)))
		undoinfo = []
		for card in self.project.input["REGION"]:
			if card.ignore(): continue
			expr = pat.sub("", card.extra())
			if expr != card.extra():
				undoinfo.append(self.setWhatUndo(card,-1,expr))
		return undoinfo

	# ----------------------------------------------------------------------
	# Remove all body references
	# Warning expects input to be correctly numbered
	# ----------------------------------------------------------------------
	def removeRegionRefUndo(self, region):
		# Search all cards referring to region and remove them
		remove = []
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			case = card.info.findCase(card)
			lst  = card.info.find("ri", case)
			if isinstance(lst, tuple):
				lst = list(lst)
				lst.pop()	# skip step
			for w in lst:
				if card.what(w) == region:
					remove.append(card.pos())
					break
		undoinfo = []
		remove.reverse()
		for pos in remove:
			undoinfo.append(self.delCardUndo(pos,True))
		return undoinfo

	# ----------------------------------------------------------------------
	# Every second call it refreshes the display
	# It needs to be called twice on each undo:
	# - in the beginning
	# - in the end
	# ----------------------------------------------------------------------
	def refreshUndo(self):
		try:
			self._refresh = not self._refresh
		except AttributeError:
			self._refresh = False

		if self._refresh:
			if self.geoedit.created():
				self.geoedit.refreshUndo()
			self.redraw()
		return (self.refreshUndo,)

	# ----------------------------------------------------------------------
	# WARNING: renumberUndo should always come in pairs
	# one before and one in the end of the operation
	#
	# Every second call it renumbers the cards
	# It needs to be called twice on each undo:
	# - in the beginning
	# - in the end
	# ----------------------------------------------------------------------
	def renumberUndo(self):
		self._renumber = not self._renumber

		if not self._renumber:
			self._renumberPos = sys.maxsize	# an invalid number

		undoinfo = (self.renumberUndo, )

		if self._renumber and self._renumberPos<sys.maxsize:
			self.project.input.renumber(self._renumberPos-1)
			self.project.input.scanUnits()	# FIXME

		return undoinfo

	# ----------------------------------------------------------------------
	def saveCardUndoInfo(self, pos):
		card    = self.project.input.cardlist[pos]
		comment = card.comment()
		tag     = card.tag
		whats   = card.whats()[:]	# clone list
		extra   = card.extra()
		return ("Save Card", self.restoreCardUndo,
				pos, comment, tag, whats, extra)

	# ----------------------------------------------------------------------
	def restoreCardUndo(self, pos, comment, tag, whats, extra):
		undoinfo = self.saveCardUndoInfo(pos)
		card = self.project.input.cardlist[pos]
		card.setComment(comment)
		if card.tag != tag:
			self.project.input.changeTag(card, tag)
		card.setWhats(whats)
		card.setExtra(extra)
		return undoinfo

	# ----------------------------------------------------------------------
	def setEnableUndo(self, card, enable):
		if card.tag == "#include":
			undoinfo = self.undolistRefresh()
			undoinfo.append((self.setEnableUndo, card, card.enable))
			card.setEnable(enable)
			self.project.input.include(card)
			undoinfo.append(self.renumberUndo())
			undoinfo.append(self.refreshUndo())
		else:
			undoinfo = (self.setEnableUndo, card, card.enable)
			card.setEnable(enable)
			if card.isGeo() and len(card.tag)==3:
				self.project.input.clearCache("bodies")
			elif card.tag in Input._INDENT_INC or card.tag in Input._INDENT_DEC:
				self.project.input.renumber(card.pos())

		self.setInputModified()
		return undoinfo

	# ----------------------------------------------------------------------
	def changeTagUndo(self, card, tag):
		undoinfo = (self.changeTagUndo, card, card.tag)
		self.project.input.changeTag(card, tag)
		self.setInputModified()
		return undoinfo

	# ----------------------------------------------------------------------
	def changeNameUndo(self, kind, old, new):
		self.project.input.changeName(kind, old, new)
		return (self.changeNameUndo, kind, new, old)

	#-----------------------------------------------------------------------
	# Insert or delete $start_xxx $end_xxx cards in order to make an input
	# that reflects body properties
	#-----------------------------------------------------------------------
	def correctBodyPropertiesUndo(self):
		translat  = None
		transform = None
		expansion = None
		idx = 0

		undoinfo = [self.refreshUndo(), self.renumberUndo()]

		def insCard(tag, whats=None):
			#nonlocal idx (only in python3)
			c = Input.Card(tag, whats)
			undoinfo.append(self.addCardUndo(c, idx-1, False))
			#idx += 1
			return c

		cardlist = self.project.input.cardlist

		while idx < len(cardlist):
			card = cardlist[idx]
			idx += 1
			prepro = card.tag[0]=="#"
			if card.ignore() or not(card.isGeo() or prepro): continue

			if len(card.tag)==3 or prepro:
				# if body has changed translation
				try: dx = float(card.get("@dx",0.0))
				except: dx = 0.0
				try: dy = float(card.get("@dy",0.0))
				except: dy = 0.0
				try: dz = float(card.get("@dz",0.0))
				except: dz = 0.0
				try: e = float(card.get("@expansion",0.0))
				except: e = 0.0
				t = card.get("@transform","")

				if translat is None:
					if dx!=0.0 or dy!=0.0 or dz!=0.0:
						translat = insCard("$start_translat", ["", dx, dy, dz])
						idx += 1
				else:
					if dx==0.0 and dy==0.0 and dz==0.0:
						# remove any previous translat card
						# Insert an end_translat to close previous one
						if cardlist[idx-2].tag == "$start_translat":
							undoinfo.append(self.delCardUndo(idx-2, False))
							idx -= 1
						else:
							insCard("$end_translat")
							idx += 1
						translat = None

					elif abs(dx-translat.numWhat(1)) > 1e-10 or \
					     abs(dy-translat.numWhat(2)) > 1e-10 or \
					     abs(dz-translat.numWhat(3)) > 1e-10:
						if cardlist[idx-2].tag == "$start_translat":
							undoinfo.append(self.delCardUndo(idx-2, False))
							idx -= 1
						else:
							insCard("$end_translat")
							idx += 1
						translat = insCard("$start_translat", ["", dx, dy, dz])
						idx += 1

				if expansion is None:
					if e!=0.0 and e!=1.0:
						expansion = insCard("$start_expansion", ["", e])
						idx += 1
				else:
					if e==0.0 or e==1.0:
						if cardlist[idx-2].tag == "$start_expansion":
							undoinfo.append(self.delCardUndo(idx-2, False))
							idx -= 1
						else:
							insCard("$end_expansion")
							idx += 1
						expansion = None

					elif abs(e-expansion.numWhat(1)) > 1e-10:
						if cardlist[idx-2].tag == "$start_expansion":
							undoinfo.append(self.delCardUndo(idx-2, False))
							idx -= 1
						else:
							insCard("$end_expansion")
							idx += 1
						expansion = insCard("$start_expansion", ["", e])
						idx += 1

				if transform is None:
					if t != "":
						transform = insCard("$start_transform", ["", t])
						idx += 1
				else:
					if t == "":
						if cardlist[idx-2].tag == "$start_transform":
							undoinfo.append(self.delCardUndo(idx-2, False))
							idx -= 1
						else:
							insCard("$end_transform")
							idx += 1
						transform = None

					elif t != transform.what(1):
						if cardlist[idx-2].tag == "$start_transform":
							undoinfo.append(self.delCardUndo(idx-2, False))
							idx -= 1
						else:
							insCard("$end_transform")
							idx += 1
						transform = insCard("$start_transform", ["", t])
						idx += 1

				if card.tag=="END": break

			elif card.tag[0]=="$":
				if card.tag=="$start_translat":
					if translat is not None:
						# Close any pending translat
						insCard("$end_translat")
						idx += 1
					translat = card
				elif card.tag=="$start_transform":
					if transform is not None:
						insCard("$end_transform")
						idx += 1
					transform = card
				elif card.tag=="$start_expansion":
					if expansion is not None:
						insCard("$end_expansion")
						idx += 1
					expansion = card

				elif card.tag=="$end_translat":
					if translat is None:
						idx -= 1
						undoinfo.append(self.delCardUndo(idx, False))
					translat = None
				elif card.tag=="$end_transform":
					if transform is None:
						idx -= 1
						undoinfo.append(self.delCardUndo(idx, False))
					transform = None
				elif card.tag=="$end_expansion":
					if expansion is None:
						idx -= 1
						undoinfo.append(self.delCardUndo(idx, False))
					expansion = None

		undoinfo.append(self.renumberUndo())
		undoinfo.append(self.refreshUndo())

		return undo.createListUndo(undoinfo)

	# ----------------------------------------------------------------------
	# Base class to start plotting (Launches gnuplot)
	# ----------------------------------------------------------------------
	def startPlotEngine(self, engine):
		self.plotEngine = self._plotEngines.get(engine)
		if self.plotEngine is None:
			if engine == PLOT_ENGINES[0]:
				log.say("Starting gnuplot")
				self.plotEngine = GPPlot.start()
				if self.plotEngine is None:
					self.notify("Gnuplot not responding",
						"Gnuplot program is not responding. Please " \
						"ensure that the path is correct or the program exist",
						tkFlair.NOTIFY_ERROR)
					return None
				log.say("Gnuplot version = %g"%(self.plotEngine.version()))
				log.say("Default terminal = %s"%(Gnuplot.defaultTerminal))

			elif engine == PLOT_ENGINES[1]:
				self.plotEngine = MPPlot.start()
				if self.plotEngine is None:
					self.notify("Matplotlib not responding",
						"Matplotlib python library is not responding. Please " \
						"ensure that the necessary files are installed",
						tkFlair.NOTIFY_ERROR)
					return None
				log.say("Matplotlib version = %g"%(self.plotEngine.version()))

			else:
				self.notify("Plot engine",
					"Cannot find plot engine %s"%(engine),
					tkFlair.NOTIFY_ERROR)
				return None

			self._plotEngines[engine] = self.plotEngine

		if self.plotEngine.version() == 0:
			return None

		self.plotEngine.cd(self.project.dir)
		return self.plotEngine

#-------------------------------------------------------------------------------
# New flair window
#-------------------------------------------------------------------------------
def newWindow(event=None):
	Flair()

#-------------------------------------------------------------------------------
# Material Database I/O
#-------------------------------------------------------------------------------
def loadMaterialIni():
	Materials.init([os.path.join(tkFlair.prgDir, "db/material.ini"),
			tkFlair.iniMaterial])

#-------------------------------------------------------------------------------
def saveMaterial():
	if not Materials.changed: return
	tkFlair.createIniDir()
	matini = Materials.makeIni()
	try:
		fini = open(tkFlair.iniMaterial, "w")
	except IOError:
		messagebox.showerror("Cannot open file",
			"Cannot open file %s for writing"%(tkFlair.iniMaterial))
		tkFlair.write("Error: %s cannot write %s file\n"%(tkFlair.__name__,
				tkFlair.iniMaterial))
		return
	matini.write(fini)
	fini.close()
	Materials.changed = False

#-------------------------------------------------------------------------------
# Add Recent
#-------------------------------------------------------------------------------
def addRecent(filename):
	try:
		sfn = str(filename)
	except UnicodeEncodeError:
		sfn = filename.encode("utf-8")
	try:
		found = tkFlair._maxRecent-1
		for i in range(tkFlair._maxRecent):
			if tkFlair.config.get(tkFlair.__name__, \
					"recent.%d"%(i)) == sfn:
				if i==0: return
				found = i-1

	except configparser.NoOptionError:
		pass

	# FIXME should modify the recent menu
	# Shift everything by one
	for i in range(found, -1, -1):
		try: fn = tkFlair.config.get(tkFlair.__name__, "recent.%d"%(i))
		except: continue
		tkFlair.config.set(tkFlair.__name__, "recent.%d"%(i+1), fn)

	tkFlair.config.set(tkFlair.__name__, "recent.0", sfn)

#-------------------------------------------------------------------------------
# check for errors
#-------------------------------------------------------------------------------
def check4Errors():
	if tkFlair.errors and tkFlair._errorReport:
		tkFlair.ReportDialog.sendErrorReport()

#-------------------------------------------------------------------------------
# Global Program Exit
#-------------------------------------------------------------------------------
def programExit(event=None):
	global root
	for child in root.winfo_children():
		if isinstance(child, Flair):
			try:
				child.close()
			except TclError:
				pass
		else:
			try:
				child.destroy()
			except TclError:
				pass

#-------------------------------------------------------------------------------
# show short help
#-------------------------------------------------------------------------------
def usage(rc=0):
	sys.stdout.write("\n")
	sys.stdout.write("Usage: %s [options] [filename]\n" \
			% (tkFlair.__name__))
	sys.stdout.write("\n")
	sys.stdout.write("Options:\n")
	sys.stdout.write("\t-1\t\tLoad the first flair file in the folder\n")
	sys.stdout.write("\t-c | --compile\tCompile executable\n")
#	sys.stdout.write("\t-d\t\tActivate the beta-development features\n")
#	sys.stdout.write("\t-D\t\tDeactivate beta-development features (default)\n")
	sys.stdout.write("\t-i file\t\tAlternative flair ini file\n")
	sys.stdout.write("\t--ini file\t... (Default to $HOME/.flair/flair.ini)\n")
	sys.stdout.write("\t-e executable\tFluka executable. (Default to flukahp)\n")
	sys.stdout.write("\t--exe executable\n")
	sys.stdout.write("\t-r | --recent\tLoad most recent project\n")
	sys.stdout.write("\t-R #\t\tLoad recent project (number 1..10 or filename)\n")
	sys.stdout.write("\t-l | --list\tList recent projects\n")
	sys.stdout.write("\t-s\t\tSkip About dialog\n")
	sys.stdout.write("\t-u\tupdate\tRecalculate and save input file variables\n")
	sys.stdout.write("\tfilename\tproject file or input file or imported files\n")
	sys.stdout.write("\t\t\textensions supported: <none|.flair>,  .inp, .fluka,\n")
	sys.stdout.write("\t\t\t                      .pickle, .mcnp, .gdml\n")
	sys.stdout.write("\n")
	sys.stdout.write("Environment variables:\n")
	sys.stdout.write("  PYTHON\tSet this variable only in the case that a different\n")
	sys.stdout.write("\t\tpython installation is needed. The variable should\n")
	sys.stdout.write("\t\tpoint the python executable\n\n")
	sys.stdout.write("Author:\t\t%s\n"      % (tkFlair.__author__))
	sys.stdout.write("Email:\t\t%s\n"       % (tkFlair.__email__))
	sys.stdout.write("Web:\t\t%s\n"         % (tkFlair.__www__))
	sys.stdout.write("Created:\t%s\n"       % (tkFlair.__created__))
	sys.stdout.write("Last Change:\t%s\n"   % (tkFlair.__lastchange__))
	sys.stdout.write("Version:\t%s [R%s]\n" % (tkFlair.__version__, tkFlair.__revision__))
	sys.exit(rc)

#-------------------------------------------------------------------------------
# parse command line arguments
#-------------------------------------------------------------------------------
def parseArgs(arglist): 
    # This defines a function named 'parseArgs' that takes a parameter 'arglist'
	project = Project.Project() 
    # Creates an instance of the 'Project.Project' class and assigns it to the variable 'project'
	try:
        # Starts a try-except block to handle potential errors
		optlist, args = getopt.getopt(arglist,
					"?hci:e:R:rlsdDg1u",
					["help", "compile", "ini=", "exe=", "recent", "list", "skip", "update"])
        # Uses the 'getopt' module to parse the 'arglist' and gets the option list and remaining arguments
	except getopt.GetoptError: 
        # Handles the exception if there's an error during the parsing
		usage(1) 
        # Calls the 'usage' function with the parameter 1

	recent = None 
    # Initializes the variable'recent' to None
	update = False 
    # Initializes the variable 'update' to False
	compileExe = False 
    # Initializes the variable 'compileExe' to False

	for opt, val in optlist: 
        # Starts a loop to iterate over the parsed option list
		if opt in ("-h", "-?", "--help"): 
            # Checks if the option is one of the help options
			usage(0) 
            # Calls the 'usage' function with the parameter 0
		elif opt in ("-i", "--ini"): 
            # Checks if the option is related to 'ini'
			pass 
            # Does nothing for this option (might be implemented later)
		elif opt in ("-c", "--compile"): 
            # Checks if the option is for compilation
			compileExe = True 
            # Sets the 'compileExe' flag to True
		elif opt in ("-s", "--skip"): 
            # Checks if the option is for skipping
			tkFlair._SKIP_INTRO = True 
            # Sets the '_SKIP_INTRO' attribute of 'tkFlair' to True
		elif opt in ("-e", "--exe"): 
            # Checks if the option is related to an executable
			project.setExecutable(val) 
            # Calls the'setExecutable' method of the 'project' object with the given value
		elif opt in ("-d","-D"): 
            # Checks for specific options
			pass  
            # Does nothing for this option (might be handled elsewhere)
		elif opt == "-1": 
            # Checks for a specific option
			for filename in os.listdir("."): 
                # Loops over the files in the current directory
				fn,ext = os.path.splitext(filename) 
                # Gets the file name and extension
				if ext == ".flair": 
                    # Checks if the extension is '.flair'
					try:
                        # Starts a try block
						project.load(filename) 
                        # Attempts to load the project
					except (IOError, OSError):
                        # Handles file-related errors
						project.clear() 
                        # Clears the project
						tkFlair.write("Error: %s\n"%(sys.exc_info()[1]))
                        # Writes the error message
            #... (similar explanations for other options and conditions follow)
	if recent: # Checks if'recent' has a value
		args.append(recent) # Appends the value to the 'args' list
		 
	if len(args) == 0:
		print("zaizaiwoaini!!")
	view = []

	for filename in args: 
        # Loops over the 'args' list
		(name, ext) = os.path.splitext(filename) 
        # Gets the file name and extension		
		if ext in ("","."): 
            # Checks for empty or dot extension
			try:  
                # Starts a try block
				filename = name+".flair"
				os.stat(filename)
				ext = ".flair"
			except OSError:
                # Handles file-related errors
				try:  
                    # Starts another try block
					filename = name+".inp"
					os.stat(filename)
					ext = ".inp"
				except OSError:
                    # Handles more file-related errors
					project.clear()
					project.name = name+".flair"
					project.inputName = name
					project.inputFile = name+".inp"

		if ext in (".inp", ".fluka"): 
            # Checks for specific extensions
			try:
                # Starts a try block
				project.loadInput(filename)
			except:
                # Handles loading errors
				project.clear()
				tkFlair.write("Error: %s\n"%(sys.exc_info()[1]))
				traceback.print_exc()

		elif ext in (".mcnp", ".m"): 
            # Checks for other specific extensions
			mcnp = Mcnp.Mcnp()
			mcnp.toFluka(filename, project.input)

		elif ext == ".gdml": 
            # Checks for a specific extension
			gdml = Gdml.Gdml(filename)
			gdml.toFluka(project.input)

		elif ext == ".pickle": 
            # Checks for a specific extension
			try:
				f = bz2.BZ2File(filename,'r')
				project.input.load(pickle.Unpickler(f))
				f.close()
			except:
                # Handles pickle-related errors
				project.clear()
				tkFlair.write("Error: %s\n"%(sys.exc_info()[1]))
				traceback.print_exc()
		elif ext == ".flair": 
            # Checks for a specific extension
			try:
				project.load(filename)
			except (IOError, OSError):
                # Handles loading errors
				project.clear()
				project.name = filename
				name,ext = os.path.splitext(project.name)
				project.inputName = name
				project.inputFile = name+".inp"
		else:
            # Handles other extensions
			view.append(filename)
	if update: 
        # Checks the 'update' flag
		project.saveInput()
        # Saves the input
		sys.exit(0)
        # Exits the program

	if compileExe: 
        # Checks the 'compileExe' flag
		import CompileProcess
		p = CompileProcess.CompileProcess(project)
		p.output = lambda x:sys.stdout.write(x)
		p.execute()
		p.output("%s\n%s\n"%(p.message))
		sys.exit(p.rc)
	return project, view 
    # Returns the 'project' object and the 'view' list
#===============================================================================
# Main subroutine
#===============================================================================
if __name__ == "__main__":
	#signal.signal(signal.SIGCHLD, signal.SIG_IGN)
	# Convert path to absolute
	sys.argv[0] = os.path.abspath(sys.argv[0])

	#Tkinter.NoDefaultRoot()
	try:
		root = Tk()
	except TclError:
		sys.stderr.write("Error: %s\n"%(str(sys.exc_info()[1])))
		sys.exit(-1)

	root.title(tkFlair.__name__)
	tkinter.CallWrapper = tkFlair.CallWrapper
	tkExtra.bindClasses(root)
	root.withdraw()

	# Check for alternative ini file
	ini = None
	for opt in ("-i","--ini"):
		if opt in sys.argv:
			try:
				ini = sys.argv[sys.argv.index(opt)+1]
				break
			except IndexError:
				pass

	# Initialise options
	PeriodicTable.init(os.path.join(tkFlair.prgDir, "db/isotopes.ini"))
	loadMaterialIni()
	Layout.init()
	tkFlair.init(ini)

	if not GeometryEditor.installed():
		tkFlair.write("Warning:  Viewer not found\n")
		tkFlair.write("Error: %s\n"%(GeometryEditor.geoviewer_error()))

	if not GeometryEditor.PILInstalled():
		tkFlair.write("Warning: PIL.Image and PIL.ImageTk not found\n")

	# Temporary set log to tkFlair.write. It will be overridden once the
	# output page is created
	log.set(tkFlair.write)

	project,view = parseArgs(sys.argv[1:])

	tkFlair.addOptions(tkFlair.__name__, root)
	tkFlair.addOptions(tkFlair._FONT_SECTION, root)
	tkFlair.addOptions(tkFlair._COLOR_SECTION, root)

	InputPage.configGet()
	if GeometryEditor.installed(): GeometryEditor.configGet()

	if project.projFile: addRecent(project.projFile)

	flair = Flair(root, project, view)

	if not tkFlair._SKIP_INTRO: tkFlair.aboutDialog(flair, 2500)
	if Project.flukaDir == "":
		messagebox.showinfo("Set Fluka Directory",
			"The environment variable %s is not found. "
			"Please set the Fluka Directory in the Preferences Dialog." \
			% (Project.flukaVar))
		flair.preferences()
	Manual.init((tkFlair.prgDir, Project.flukaDir))
	if Updates.need2Check(): flair.checkUpdates()
	flair.checkFluka()
	GeometryEditor.checkVersion(flair)

	# Main loop
	if root is not None: root.mainloop()
#		count = 0
#		while True:
#			try:
#				root.mainloop()
#				break
#			except KeyboardInterrupt:
#				count += 1
#				log.say("***Interrupt***")
	Manual.done()
	tkFlair.fini()

	# Write initialization files
	tkFlair.saveStats()
	tkFlair.writeIni()
	saveMaterial()
#	saveObject()

	import gc
	del flair
	del root
	if Input._developer: tkFlair.write("gc collect=%d"%(gc.collect()))
