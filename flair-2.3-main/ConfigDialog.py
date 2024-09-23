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

# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	24-Dec-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os

from tkinter import *
import tkinter.font
import tkFlair
import tkExtra
import tkTextEditor
import bFileDialog

import Input
import Manual
import Project
import configparser
try:
	import multiprocessing
except:
	multiprocessing = None

_FONT_SIZE = (   '6',  '7',  '8',  '9', '10', '11', '12', '13',
		'14', '16', '18', '20', '22', '24', '26', '28',
		'32', '36', '40', '48', '56', '64', '72' )

_FONT_LIST = [
		"balloon",
		"input.comment",
		"input.fixed",
		"input.hidden",
		"input.label",
		"input.tag",
		"input.value",
		"ribbon",
		"ribbontab",
		"title",
		"tree",
		"TkDefaultFont",
		"TkFixedFont",
		"TkMenuFont",
		"TkTextFont" ]

_COLOR_LIST = [ "input.active",
		"input.comment",
		"input.define",
		"input.diff",
		"input.disable",
		"input.error",
		"input.hidden",
		"input.label",
		"input.normal",
		"input.sdum",
		"input.tag",
		"input.value",

		"geometry.background",
		"geometry.bodybboxincolor",
		"geometry.bodybboxoutcolor",
		"geometry.zonebboxcolor",
		"geometry.regionbboxcolor",
		"geometry.error",
		"geometry.gridtext",
		"geometry.label",
		"geometry.lattice",
		"geometry.lock",
		"geometry.region",
		"geometry.select",
		"geometry.title",
		"geometry.vertex",
		"geometry.visible",
		"geometry.wireframe",
		"geometry.zone",
		"geometry.zonebody",

		"balloon.foreground",
		"balloon.background",
		"title.foreground",
		"title.background",
		"title.backinactive",
		"entrylabel.foreground",
		"entrylabel.background",
		"infolabel.foreground",
		"infolabel.background",
		"calclabel.foreground",
		"calclabel.background",
#		"*entry*foreground",
#		"*entry*background",
#		"*label*foreground",
#		"*label*background",
#		"*text*foreground",
#		"*text*background",
		]

#===============================================================================
# Configuration Dialog
#===============================================================================
class ConfigDialog(Toplevel):
	activepage = None
	width  = -1
	height = -1
	_PAGES = None

	def __init__(self, master, config):
		Toplevel.__init__(self, master)
		self.transient(master)
		self.config = config
		self.title("Preferences")

		# Action buttons
		frameButtons = Frame(self)
		frameButtons.pack(side=BOTTOM)
		self.buttonHelp = Button(frameButtons,
				text='Help',
				image=tkFlair.icons["info"],
				compound=LEFT,
				command=self.help, takefocus=FALSE)
		self.buttonOk = Button(frameButtons,
				text='Ok',
				image=tkFlair.icons["ok"],
				compound=LEFT,
				command=self.ok,
				takefocus=FALSE)
		self.buttonCancel = Button(frameButtons,
				text='Cancel',
				image=tkFlair.icons["x"],
				compound=LEFT,
				command=self.cancel, takefocus=FALSE)
		self.buttonHelp.pack(side=RIGHT, padx=5, pady=5)
		self.buttonOk.pack(side=LEFT, padx=5, pady=5)
		self.buttonCancel.pack(side=LEFT, padx=5, pady=5)

		# Selection listbox
		self.listbox = Listbox(self,
				width=15,
				height=20,
				selectmode=SINGLE,
				selectborderwidth=0,
				background="White",
				takefocus=True,
				exportselection=FALSE,
				borderwidth=2,
				relief=SUNKEN)
		self.listbox.pack(side=LEFT, fill=Y, padx=2, pady=5)
		self.listbox.bind("<<ListboxSelect>>", self.listboxChanged)

		# main frame
		#self.labelframe = LabelFrame(self, text="")
		#self.labelframe.pack(expand=YES, fill=BOTH, padx=2, pady=5)
		self.mainframe = Frame(self, borderwidth=2, relief=SUNKEN)
		self.mainframe.pack(expand=YES, fill=BOTH, padx=2, pady=5)

		self.initFontDictionary()
		self.initColorDictionary()

		# Initialize pages, but only after the class is defined!
		if ConfigDialog._PAGES is None:
			ConfigDialog._PAGES = (
			   ("Programs",	ConfigDialog.programsFrame),
			   ("Colors",	ConfigDialog.colorsFrame),
			   ("Data",	ConfigDialog.dataProcessingFrame),
			   ("Fonts",	ConfigDialog.fontsFrame),
			   ("Geometry",	ConfigDialog.geometryFrame),
			   ("Gnuplot",	ConfigDialog.gnuplotFrame),
			   ("Interface",ConfigDialog.interfaceFrame),
			   ("Run",	ConfigDialog.batchSystemFrame))

		# Populate pages
		self.pages = {}
		for page, method in ConfigDialog._PAGES:
			self.listbox.insert(END, page)
			self.pages[page] = method(self)

		self.bind('<F1>', self.help)
		self.bind('<Escape>', self.cancel)

		self.loadFonts()
		self.rc = False

	# --------------------------------------------------------------------
	def show(self):
		# Select last page
		self.listbox.selection_clear(0,END)
		if ConfigDialog.activepage is None:
			self.showPage(ConfigDialog._PAGES[0][0])
		else:
			self.showPage(ConfigDialog.activepage)

		if ConfigDialog.width > 0:
			self.geometry("%dx%d" % \
				(ConfigDialog.width, ConfigDialog.height))
		self.init()

		self.deiconify()
		self.wait_visibility()
		# Get the focus
		self.focus_set()
		self.grab_set()
		self.wait_window()
		return self.rc

	# --------------------------------------------------------------------
	# Display new page
	# --------------------------------------------------------------------
	def listboxChanged(self, event):
		sel = self.listbox.get(self.listbox.curselection())
		if sel: self.showPage(sel)

	# --------------------------------------------------------------------
	def showPage(self, page):
		for i in range(len(ConfigDialog._PAGES)):
			if page == ConfigDialog._PAGES[i][0]:
				if not self.listbox.selection_includes(i):
					self.listbox.selection_set(i)
				break

		if ConfigDialog.activepage:
			self.pages[ConfigDialog.activepage].pack_forget()
		ConfigDialog.activepage = page
		self.pages[page].pack(expand=YES, fill=BOTH)

	# --------------------------------------------------------------------
	def initFontDictionary(self):
		# save font dictionary...
		self.fontDict = {}

		# and populate font list with extra fonts
		for n,v in self.config.items(tkFlair._FONT_SECTION):
			if n not in _FONT_LIST:
				_FONT_LIST.append(n)

			vt = tuple(v.split(','))
			try:    size = int(vt[1])
			except: size = 12

			weight = tkinter.font.NORMAL
			slant  = tkinter.font.ROMAN

			try:
				if vt[2].find(tkinter.font.BOLD)>=0:
					weight = tkinter.font.BOLD

				# FIXME old version
				if vt[2].find(tkinter.font.ITALIC)>=0:
					slant = tkinter.font.ITALIC
			except:
				pass

			try:
				slant = vt[3]
			except:
				try:
					if vt[2].find(tkinter.font.ITALIC)>=0:
						slant = tkinter.font.ITALIC
				except:
					pass

			self.fontDict[n] = (vt[0], size, weight, slant)

	# --------------------------------------------------------------------
	def initColorDictionary(self):
		self.colorDict = {}

		# and populate color list with extra colors
		for n,v in self.config.items(tkFlair._COLOR_SECTION):
			if n not in _COLOR_LIST:
				if not n.startswith("material.") and \
				   not n.startswith("alpha."):
					_COLOR_LIST.append(n)

			self.colorDict[n] = v
		_COLOR_LIST.sort()

	# --------------------------------------------------------------------
	# Create programs frame
	# --------------------------------------------------------------------
	def programsFrame(self):
		pageframe = Frame(self.mainframe)

		# ---
		row = 0
		if Input._developer:	# Hide fluka var for normal users
			Label(pageframe, text="Fluka Var").grid(
					row=row, column=0, sticky=W)
		self.flukavar = Entry(pageframe, background="white")
		if Input._developer:
			self.flukavar.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.flukavar,
			"Environment variable pointing to FLUKA directory. Default: FLUPRO")

		# ---
		if Input._developer: row += 1

		Label(pageframe, text="Fluka Directory").grid(
				row=row, column=0, sticky=W)
		self.flukadir = Entry(pageframe, background="white")
		self.flukadir.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.flukadir,
			"Override or set FLUKA directory. Default: ${FlukaVar}")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
			command=self.chooseDirectory)
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="rfluka").grid(
				row=row, column=0, sticky=W)
		self.rfluka = Entry(pageframe, background="white")
		self.rfluka.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.rfluka,
			"rfluka executable. Default: ${FlukaVar}/flutil/rfluka")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.rfluka:
					self.chooseFile("rfluka", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="Executable").grid(
				row=row, column=0, sticky=W)
		self.flukaexe = Entry(pageframe, background="white")
		self.flukaexe.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.flukaexe,
			"flukahp executable. Default: ${FlukaVar}/flukahp")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.flukaexe:
					self.chooseFile("flukahp", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="Compiler").grid(
				row=row, column=0, sticky=W)
		self.flufor = tkExtra.Combobox(pageframe, label=True)
		self.flufor.fill(Project.FLUFOR_OPTIONS)
		self.flufor.grid(row=row, column=1, columnspan=2, sticky=EW)
		tkExtra.Balloon.set(self.flufor,
			"Override or set FLUKA compiler. Default: ${FLUFOR}")

		# ---
		row += 1
		Label(pageframe, text="fff").grid(
				row=row, column=0, sticky=W)
		self.fff = Entry(pageframe, background="white")
		self.fff.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.fff,
			"fff compiling program. Default: ${FlukaVar}/flutil/fff")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.fff:
					self.chooseFile("fff", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="USBMAX trace scan").grid(
				row=row, column=0, sticky=W)
		self.usbmax = Entry(pageframe, background="white")
		self.usbmax.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.usbmax,
			"USBMAX max trace scanning program. Default: ${FlukaVar}/flutil/usbmax")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.usbmax:
					self.chooseFile("usbmax", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="Gplevbin projection").grid(
				row=row, column=0, sticky=W)
		self.gplevbin = Entry(pageframe, background="white")
		self.gplevbin.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.gplevbin,
			"USRBIN processing program for gnuplot. Default: ${FlukaVar}/flutil/gplevbin")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.gplevbin:
					self.chooseFile("gplevbin", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="USBREA usrbin to ascii").grid(
				row=row, column=0, sticky=W)
		self.usbrea = Entry(pageframe, background="white")
		self.usbrea.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.usbrea,
			"USBREA convert USRBIN files to ascii. Default: ${FlukaVar}/flutil/usbrea")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.usbrea:
					self.chooseFile("usbrea", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="Kill Command").grid(
				row=row, column=0, sticky=W)
		self.killCmd = Entry(pageframe, background="white")
		self.killCmd.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.killCmd,
			"Kill command. No default. If empty uses the kill system command. Override with the appropriate program for batch systems")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.killCmd:
					self.chooseFile("kill", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="Editor").grid(
				row=row, column=0, sticky=W)
		self.editor = Entry(pageframe, background="white")
		self.editor.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.editor,
			"Editor program. No default")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.editor:
					self.chooseFile("editor", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="Terminal").grid(
				row=row, column=0, sticky=W)
		self.terminal = Entry(pageframe, background="white")
		self.terminal.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.terminal,
			"Terminal program. No default")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.terminal:
					self.chooseFile("terminal", e))
		b.grid(row=row, column=2)

		row += 1
		Label(pageframe, text="Debugger").grid(
				row=row, column=0, sticky=W)
		self.debugger = Entry(pageframe, background="white")
		self.debugger.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.debugger,
			"Debugging program.")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.debugger:
					self.chooseFile("debugger", e))
		b.grid(row=row, column=2)

		row += 1
		Label(pageframe, text="Gnuplot").grid(
				row=row, column=0, sticky=W)
		self.gnuplot = Entry(pageframe, background="white")
		self.gnuplot.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.gnuplot, "Gnuplot program.")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.gnuplot:
					self.chooseFile("gnuplot", e))
		b.grid(row=row, column=2)

		pageframe.grid_columnconfigure(1, weight=1)

		return pageframe

	# --------------------------------------------------------------------
	# Choose Directory
	# --------------------------------------------------------------------
	def chooseDirectory(self):
		directory = self.flukadir.get().strip()
		if len(directory)==0:
			try:
				directory = os.getenv(self.flukavar.get())
			except:
				directory = ""
		directory = bFileDialog.askdirectory(master=self,
			title="Select FLUKA directory",
			initialdir=directory)
#		if directory > 0:
		if len(directory) >0 :
			self.flukadir.delete(0,END)
			self.flukadir.insert(0, directory)

	# --------------------------------------------------------------------
	def chooseFile(self, default, entry, flukadir=True):
		if flukadir:
			directory = self.flukadir.get().strip()
			if len(directory)==0:
				directory = os.getenv(self.flukavar.get())
		else:
			directory = tkFlair.prgDir
		fn = entry.get().strip()
		if len(fn)==0:
			fn = default
		ft=[("%s files"%(default),"*%s*"%(default)),
				("All","*")]

		fn = bFileDialog.askopenfilename(master=self,
				title="Select %s"%(default),
				initialdir=directory,
				initialfile=fn,
				filetypes=ft)
		if len(fn) > 0:
			entry.delete(0, END)
			entry.insert(0, fn)

	# --------------------------------------------------------------------
	# Create batch system frame
	# --------------------------------------------------------------------
	def batchSystemFrame(self):
		pageframe = Frame(self.mainframe, padx=5, pady=5)
		self.enableBatchSystem = BooleanVar(self)

		# row = 0
		# b = Checkbutton(pageframe, text="Enable Batch System",
		#		variable=self.enableBatchSystem,
		#		onvalue=1, offvalue=0,
		#		anchor=W, command=self.batchEdit)
		# if self.enableBatchSystem : b.select()

		#tkExtra.Balloon.set(b,
		#	"Enable Batch System")
		#b.grid(row=row, column=0, columnspan=2, sticky=W)
		#b.pack(side=TOP, fill=X)

		# ---
		row = 0
		Label(pageframe, text="Spawn Run Name").grid(
				row=row, column=0, sticky=E)
		self.spawnRunName = Entry(pageframe, background="white")
		self.spawnRunName.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.spawnRunName,
			"Prefix to be used for all temporary files created by flair")

		# --- Queue Types ---
		row += 1
		frame = LabelFrame(pageframe, text="Queues / Submit command",
				relief=GROOVE)
		frame.grid(row=row, column=0, columnspan=2, sticky=NSEW, padx=5, pady=5)

		pageframe.grid_rowconfigure(1, weight=1)
		pageframe.grid_columnconfigure(1, weight=1)

		# ----
		self.queueTypes = tkExtra.MultiListbox(frame,
			(("Queue ", 16, None),
			 ("Execution command", 30, None)))
		self.queueTypes.pack(fill=BOTH, expand=YES)

		self.queueTypes.bindList('<ButtonRelease-1>',
				self.queueRelease)
		self.queueTypes.bindList('<Return>',
				self.queueEdit)
		self.queueTypes.bindLeftRight()
		self.queueTypes.setPopupMenu(
				[('Add',    0, self.queueAdd,    tkFlair.icons["add"]),
				 ('Delete', 0, self.queueDelete, tkFlair.icons["x"]),
				 ('Clone',  0, self.queueClone,  tkFlair.icons["clone"])] )
		return pageframe

	# ----------------------------------------------------------------------
	def queueEdit(self, event):
		sel = self.queueTypes.curselection()
		if len(sel) != 1: return
		sel = int(sel[0])
		try: active = event.widget.index(ACTIVE)
		except: return
		if active != sel: return
		edit = tkExtra.InPlaceEdit(event.widget)

	# ----------------------------------------------------------------------
	def queueRelease(self, event=None):
		act = event.widget.nearest(event.y)
		if act: self.queueTypes.activate(act)
		self.queueEdit(event)

	# ----------------------------------------------------------------------
	def queueAdd(self):
		self.queueTypes.insert(END, ("myQueue", ""))
		self.queueTypes.selection_clear(0,END)
		self.queueTypes.selection_set(END)
		self.queueTypes.activate(END)
		self.queueTypes.see(END)
		self.queueTypes.focus_set()

	# ----------------------------------------------------------------------
	def queueDelete(self):
		lst = list(map(int, self.queueTypes.curselection()))
		if not lst: return
		lst.reverse()
		for i in lst:
			self.queueTypes.delete(i)

	# ----------------------------------------------------------------------
	def queueClone(self):
		sel = list(map(int,self.queueTypes.curselection()))
		if len(sel) == 0: return
		oldSize = self.queueTypes.size()
		for i in sel:
			item = self.queueTypes.get(i)
			self.queueTypes.insert(END, item)
		self.queueTypes.selection_clear(0,END)
		self.queueTypes.selection_set(oldSize,END)
		self.queueTypes.activate(END)

	# ----------------------------------------------------------------------
	def batchEdit(self):
		if self.enableBatchSystem.get():
			self.queueTypes["state"] = NORMAL
		else:
			self.queueTypes["state"] = DISABLED

	# --------------------------------------------------------------------
	# Graphical Interface frame
	# --------------------------------------------------------------------
	def interfaceFrame(self):
		pageframe = Frame(self.mainframe)

		self.skipintro     = BooleanVar(self)
		self.flukafiles    = BooleanVar(self)
		self.keepbackup    = BooleanVar(self)
		self.cleanup       = BooleanVar(self)
		self.warnmultmat   = BooleanVar(self)

		row = 0
		b = Checkbutton(pageframe, text="Skip About dialog",
				variable=self.skipintro,
				onvalue=1, offvalue=0,
				anchor=W)
		tkExtra.Balloon.set(b,
			"Skip the about dialog during startup")
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		b = Checkbutton(pageframe, text="Warn on multiple material assignment to region",
				variable=self.warnmultmat,
				onvalue=1, offvalue=0,
				anchor=W)
		tkExtra.Balloon.set(b,
			"Show warnings if a region is assigned multiple materials")
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		b = Checkbutton(pageframe, text="Show fluka files in fluka_XXX dir",
				variable=self.flukafiles,
				onvalue=1, offvalue=0,
				anchor=W)
		tkExtra.Balloon.set(b,
			"Show the fluka files (usually links) in the running fluka_XXX directory listing")
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		b = Checkbutton(pageframe, text="Keep backups",
				variable=self.keepbackup,
				onvalue=1, offvalue=0,
				anchor=W)
		b.grid(row=row, column=0, columnspan=2, sticky=W)
		tkExtra.Balloon.set(b,
			"Flag whether you want to keep backups of input and project files")

		# ---
		row += 1
		b = Checkbutton(pageframe, text="Cleanup temporary files",
				variable=self.cleanup,
				onvalue=1, offvalue=0,
				anchor=W)
		b.grid(row=row, column=0, columnspan=2, sticky=W)
		tkExtra.Balloon.set(b,
			"Cleanup flair temporary files (Recommended: ON)")

		# ---
		row += 1
		Label(pageframe, text="Icon zoom").grid(
				row=row, column=0, sticky=W)
		self.iconZoom = Scale(pageframe, from_=1, to_=5,
				orient=HORIZONTAL, showvalue=1)
		self.iconZoom.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.iconZoom, "Zoom icons by an integer factor")

		# ---
		row += 1
		Label(pageframe, text="Temporary prefix").grid(
				row=row, column=0, sticky=E)
		self.tmpprefix = Entry(pageframe, background="white")
		self.tmpprefix.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.tmpprefix,
			"Prefix to be used for all temporary files created by flair")

		# ---
		row += 1
		Label(pageframe, text="Auto save time").grid(
				row=row, column=0, sticky=E)
		self.autosave = tkExtra.IntegerEntry(pageframe, background="white")
		self.autosave.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.autosave,
			"Input file auto save period. File is written in pickle format")
		Label(pageframe, text="[s]").grid(row=row, column=2, sticky=W)

		# ---
		#row += 1
		#Label(pageframe, text="Key time threshold").grid(
		#	row=row, column=0, sticky=E)
		#self.keytime = tkExtra.IntegerEntry(pageframe, width=7, background="white")
		#self.keytime.grid(row=row, column=1, sticky=EW)
		#tkExtra.Balloon.set(self.keytime,
		#	"Set the keyboard time out threshold during listbox searches")
		#Label(pageframe, text="[ms]").grid(row=row, column=2, sticky=W)

		# ---
		#row += 1
		#Label(pageframe, text="Balloon delay").grid(
		#	row=row, column=0, sticky=E)
		#self.balloondelay = tkExtra.IntegerEntry(pageframe, width=7, background="white")
		#self.balloondelay.grid(row=row, column=1, sticky=EW)
		#tkExtra.Balloon.set(self.balloondelay,
		#	"Delay time to display the help-balloon")
		#Label(pageframe, text="[ms]").grid(row=row, column=2, sticky=W)

		# ---
		row += 1
		Label(pageframe, text="Attach timeout").grid(
				row=row, column=0, sticky=E)
		self.timethreshold = tkExtra.IntegerEntry(pageframe, background="white")
		self.timethreshold.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.timethreshold,
			"Time out threshold for attaching to a running job")
		Label(pageframe, text="[s]").grid(row=row, column=2, sticky=W)

		# ---
		row += 1
		Label(pageframe, text="Refresh Interval (s)").grid(
				row=row, column=0, sticky=W)
		self.refreshinterval = tkExtra.IntegerEntry(pageframe, background="white")
		self.refreshinterval.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.refreshinterval,
			"Refresh interval for updating information on running programs")
		Label(pageframe, text="[s]").grid(row=row, column=2, sticky=W)

		# ---
		row += 1
		Label(pageframe, text="Time format").grid(
			row=row, column=0, sticky=E)
		self.timeformat = Entry(pageframe, width=20, background="white")
		self.timeformat.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.timeformat,
			"Set the format to display time in file dialogs. Check strftime python syntax")

		pageframe.grid_columnconfigure(1, weight=1)

		return pageframe

	# --------------------------------------------------------------------
	def dataProcessingFrame(self):
		pageframe = Frame(self.mainframe)

		row = 0
		self.usrTypes = tkExtra.MultiListbox(pageframe,
			(("Type", 12, None),
			 ("Default Output File Name", 20, None)),
			 height=6)
		self.usrTypes.grid(row=row, column=0, columnspan=3, sticky=NSEW)
		self.usrTypes.bindList('<ButtonRelease-1>',
				self.usrTypesRelease)
		self.usrTypes.bindLeftRight()

		# ---
		row += 1
		self.usrRules = tkExtra.ExListbox(pageframe,
					height=4,
					takefocus=TRUE,
					exportselection=TRUE)
		self.usrRules.grid(row=row, column=0, columnspan=3, sticky=NSEW)
		self.usrRules.bind('<ButtonRelease-1>', self.usrRulesRelease)
		self.usrRules.setPopupMenu(
				[('Add',    0, self.usrRulesAdd,    tkFlair.icons["add"]),
				 ('Delete', 0, self.usrRulesDelete, tkFlair.icons["x"]),
				 ('Clone',  0, self.usrRulesClone,  tkFlair.icons["clone"])] )

		# ---
		row += 1
		Label(pageframe, text="USRBIN processing").grid(
				row=row, column=0, sticky=W)
		self.usbsuw = Entry(pageframe, background="white")
		self.usbsuw.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.usbsuw,
			"USRBIN processing program. Default: ${FlukaVar}/flutil/usbsuw")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.usbsuw:
					self.chooseFile("usbsuw", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="USRTRACK processing").grid(
				row=row, column=0, sticky=W)
		self.ustsuw = Entry(pageframe, background="white")
		self.ustsuw.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.ustsuw,
			"USRTRACK/USRCOLL processing program. Default: ${FlukaVar}/flutil/ustsuw")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.ustsuw:
					self.chooseFile("ustsuw", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="USRBDX processing").grid(
				row=row, column=0, sticky=W)
		self.usxsuw = Entry(pageframe, background="white")
		self.usxsuw.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.usxsuw,
			"USRBDX processing program. Default: ${FlukaVar}/flutil/usxsuw")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.usxsuw:
					self.chooseFile("usxsuw", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="USRYIELD processing").grid(
				row=row, column=0, sticky=W)
		self.usysuw = Entry(pageframe, background="white")
		self.usysuw.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.usysuw,
			"USRYIELD processing program. Default: ${FlukaVar}/flutil/usysuw")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.usysuw:
					self.chooseFile("usysuw", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="RESNUCLEI processing").grid(
				row=row, column=0, sticky=W)
		self.usrsuw = Entry(pageframe, background="white")
		self.usrsuw.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.usrsuw,
			"RESNUCLEI processing program. Default: ${FlukaVar}/flutil/usrsuw")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.usrsuw:
					self.chooseFile("usrsuw", e))
		b.grid(row=row, column=2)

		# ---
		row += 1
		Label(pageframe, text="DETECT processing").grid(
				row=row, column=0, sticky=W)
		self.detsuw = Entry(pageframe, background="white")
		self.detsuw.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.detsuw,
			"DETECT processing program. Default: ${FlukaVar}/flutil/detsuw")

		b = Button(pageframe, image=tkFlair.icons["load"], pady=0,
				command=lambda e=self.detsuw:
					self.chooseFile("detsuw", e))
		b.grid(row=row, column=2)

		# ---
		pageframe.grid_rowconfigure(0, weight=1)
		pageframe.grid_columnconfigure(1, weight=1)

		return pageframe

	# ----------------------------------------------------------------------
	def usrTypesRelease(self, event=None):
		listbox = event.widget
		act = listbox.nearest(event.y)
		if act: self.usrTypes.activate(act)
		sel = self.usrTypes.curselection()
		if len(sel) != 1: return
		sel = int(sel[0])
		try: active = listbox.index(ACTIVE)
		except: return
		if active != sel: return
		# Find list box
		lid = self.usrTypes.lists.index(listbox)
		if lid==0: return
		edit = tkExtra.InPlaceEdit(listbox)

	# ----------------------------------------------------------------------
	def usrRulesRelease(self, event=None):
		act = self.usrRules.nearest(event.y)
		if act: self.usrRules.activate(act)
		sel = self.usrRules.curselection()
		if len(sel) != 1: return
		sel = int(sel[0])
		try: active = self.usrRules.index(ACTIVE)
		except: return
		if active != sel: return
		edit = tkExtra.InPlaceEdit(self.usrRules)
		self.usrRules.selection_set(sel)

	# ----------------------------------------------------------------------
	def usrRulesAdd(self):
		self.usrRules.insert(END, r"+\I\d\d\d_fort.\U")
		self.usrRules.selection_clear(0,END)
		self.usrRules.selection_set(END)
		self.usrRules.activate(END)
		self.usrRules.see(END)
		self.usrRules.focus_set()

	# ----------------------------------------------------------------------
	def usrRulesDelete(self):
		lst = list(map(int, self.usrRules.curselection()))
		if not lst: return
		lst.reverse()
		for i in lst:
			self.usrRules.delete(i)

	# ----------------------------------------------------------------------
	def usrRulesClone(self):
		sel = list(map(int,self.usrRules.curselection()))
		if len(sel) == 0: return
		oldSize = self.usrRules.size()
		for i in sel:
			item = self.usrRules.get(i)
			self.usrRules.insert(END, item)
		self.usrRules.selection_clear(0,END)
		self.usrRules.selection_set(oldSize,END)
		self.usrRules.activate(END)

	# --------------------------------------------------------------------
	# Gnuplot Frame
	# --------------------------------------------------------------------
	def gnuplotFrame(self):
		pageframe = Frame(self.mainframe)

		# --- Project ---
		frame = Frame(pageframe)
		frame.pack(side=TOP, padx=5, pady=5, fill=X)
		Label(frame, text="Terminal:").pack(side=LEFT)
		self.gnuplotTerminal = Entry(frame, background="White")
		self.gnuplotTerminal.pack(side=RIGHT, expand=YES, fill=X)
		tkExtra.Balloon.set(self.gnuplotTerminal,
			"Gnuplot terminal. Default=empty, let gnuplot decide")

		frame = LabelFrame(pageframe, text="Global Commands",
				relief=GROOVE)
		frame.pack(side=TOP, padx=5, pady=5, fill=BOTH, expand=YES)

		self.gnuplotCommands = Text(frame, background="white", width=40, height=7)
		self.gnuplotCommands.pack(fill=BOTH, expand=YES)

		# --- File Types ---
		frame = LabelFrame(pageframe, text="File Types",
				relief=GROOVE)
		frame.pack(padx=5, pady=5, fill=BOTH, expand=YES)

		self.gnuplotFileTypes = tkExtra.MultiListbox(frame,
			(("Type", 5, None),
			 ("Settings", 30, None)))
		self.gnuplotFileTypes.pack(fill=BOTH, expand=YES)

		self.gnuplotFileTypes.bindList('<ButtonRelease-1>',
				self.gnuplotRelease)
		self.gnuplotFileTypes.bindLeftRight()
		self.gnuplotFileTypes.setPopupMenu(
				[('Add',    0, self.gnuplotAdd,    tkFlair.icons["add"]),
				 ('Delete', 0, self.gnuplotDelete, tkFlair.icons["x"]),
				 ('Clone',  0, self.gnuplotClone,  tkFlair.icons["clone"])] )
		return pageframe

	# ----------------------------------------------------------------------
	def gnuplotRelease(self, event=None):
		listbox = event.widget
		act = listbox.nearest(event.y)
		if act: self.gnuplotFileTypes.activate(act)
		sel = self.gnuplotFileTypes.curselection()
		if len(sel) != 1: return
		sel = int(sel[0])
		try: active = listbox.index(ACTIVE)
		except: return
		if active != sel: return
		edit = tkExtra.InPlaceEdit(listbox)

	# ----------------------------------------------------------------------
	def gnuplotAdd(self):
		self.gnuplotFileTypes.insert(END, (".ext", ""))
		self.gnuplotFileTypes.selection_clear(0,END)
		self.gnuplotFileTypes.selection_set(END)
		self.gnuplotFileTypes.activate(END)
		self.gnuplotFileTypes.see(END)
		self.gnuplotFileTypes.focus_set()

	# ----------------------------------------------------------------------
	def gnuplotDelete(self):
		lst = list(map(int, self.gnuplotFileTypes.curselection()))
		if not lst: return
		lst.reverse()
		for i in lst:
			self.gnuplotFileTypes.delete(i)

	# ----------------------------------------------------------------------
	def gnuplotClone(self):
		sel = list(map(int,self.gnuplotFileTypes.curselection()))
		if len(sel) == 0: return
		oldSize = self.gnuplotFileTypes.size()
		for i in sel:
			item = self.gnuplotFileTypes.get(i)
			self.gnuplotFileTypes.insert(END, item)
		self.gnuplotFileTypes.selection_clear(0,END)
		self.gnuplotFileTypes.selection_set(oldSize,END)
		self.gnuplotFileTypes.activate(END)

	# --------------------------------------------------------------------
	# Geometry Frame
	# --------------------------------------------------------------------
	def geometryFrame(self):
		pageframe = Frame(self.mainframe)

		self.laptopMode  = BooleanVar(self)
		self.invertWheel = BooleanVar(self)
		self.pinZoom     = BooleanVar(self)
		self.useQUA      = BooleanVar(self)
		self.useBOX      = BooleanVar(self)

		# ---
		row = 0
		b = Checkbutton(pageframe, text="Laptop Mode",
				variable=self.laptopMode,
				onvalue=1, offvalue=0,
				anchor=W)
		b.grid(row=row, column=0, columnspan=2, sticky=W)
		tkExtra.Balloon.set(b,
			"Enable laptop mode for geometry editor: " \
			"Swap middle-mouse button with right-mouse button." \
			"Requires restart of the program" )

		# ---
		row += 1
		b = Checkbutton(pageframe, text="Invert Mouse Wheel",
				variable=self.invertWheel,
				onvalue=1, offvalue=0,
				anchor=W)
		b.grid(row=row, column=0, columnspan=2, sticky=W)
		tkExtra.Balloon.set(b,
			"Invert the zooming with the mouse wheel." \
			"Requires restart of the program" )

		# ---
		row += 1
		b = Checkbutton(pageframe, text="Zoom on Mouse location",
				variable=self.pinZoom,
				onvalue=1, offvalue=0,
				anchor=W)
		b.grid(row=row, column=0, columnspan=2, sticky=W)
		tkExtra.Balloon.set(b,
			"Check to zoom in on mouse location. " \
			"Uncheck to zoom on viewport center. " \
			"[Alt] inverts the action" )

		# ---
		row += 1
		b = Checkbutton(pageframe, text="Prefer QUA on rotation",
				variable=self.useQUA,
				onvalue=1, offvalue=0,
				anchor=W)
		b.grid(row=row, column=0, columnspan=2, sticky=W)
		tkExtra.Balloon.set(b,
			"Create infinite QUA rotated instead of RCC " \
			"when rotating infinite cylinders [XYZ]CC or [XYZ]EC")

		# ---
		if Input._developer:
			row += 1
			b = Checkbutton(pageframe, text="Allow the use of obsolete bodies (BOX,RAW,ARB)",
					variable=self.useBOX,
					onvalue=1, offvalue=0,
					anchor=W)
			b.grid(row=row, column=0, columnspan=2, sticky=W)
			tkExtra.Balloon.set(b,
				"Allow the use of obsolete bodies BOX/RAW/WED/ARB for back compatibility")

		# ---
		row += 1
		Label(pageframe, text="# CPU/Cores").grid(
				row=row, column=0, sticky=W)
		if multiprocessing:
			cores = multiprocessing.cpu_count()
		else:
			cores = 8
		self.cores = Scale(pageframe, from_=0, to_=cores,
				orient=HORIZONTAL, showvalue=1)
		self.cores.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.cores, "Number of cores per viewport to use.\n" \
						"0=all available (default); 1=synchronous, nn=user defined.\n"\
						"Restart of flair is required.")

		# ---
		row += 1
		Label(pageframe, text="Snapping distance").grid(
				row=row, column=0, sticky=W)
		self.snap = Scale(pageframe, from_=1, to_=30,
				orient=HORIZONTAL, showvalue=1)
		self.snap.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.snap, "Distance in pixels for snapping.")

		# ---
		row += 1
		Label(pageframe, text="Zero").grid(
				row=row, column=0, sticky=W)
		self.zero = tkExtra.FloatEntry(pageframe, background="white")
		self.zero.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.zero,
			"Threshold value to be considered zero or equal in comparison.")

		# ---
		row += 1
		Label(pageframe, text="Infinite:").grid(row=row, column=0, sticky=W)
		self.infinite = tkExtra.FloatEntry(pageframe, background="white")
		self.infinite.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(self.infinite,
			"Length of cylinders when from infinite [XYZ]CC are converted to finite RCC")

		# ---
		row += 1
		Label(pageframe, text="Accuracy:").grid(row=row, column=0, sticky=W)
		self.accuracy = IntVar()
		self.accuracy.set(15)
		s = Spinbox(pageframe, text=self.accuracy, from_=3, to=17, background="white")
		s.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(s,
			"Number of digits to keep during transformation")

		pageframe.grid_columnconfigure(1, weight=1)

		return pageframe

	# --------------------------------------------------------------------
	# Fonts Frame
	# --------------------------------------------------------------------
	def fontsFrame(self):
		fontframe = Frame(self.mainframe)

		self.fontChanged = False
		self.fontSize   = StringVar(self)
		self.fontAbs    = BooleanVar(self)
		self.fontBold   = BooleanVar(self)
		self.fontItalic = BooleanVar(self)
		self.fontFamily = StringVar(self)
		self.font       = tkinter.font.Font(self,
					('helvetica', 12, 'normal'))
		self.fontSize.set('12')
		self.fontFamily.set('helvetica')

		# --- Which Font ---
		self.fontList = tkExtra.ExListbox(fontframe,
					takefocus=TRUE,
					selectmode=EXTENDED,
					exportselection=TRUE)
		self.fontList.bind("<<ListboxSelect>>", self.fontListChange)
		self.fontList.pack(side=LEFT, padx=2, pady=2,
					fill=BOTH, expand=YES)
		for font in _FONT_LIST:
			self.fontList.insert(END, font)

		# --- Font ---
		rightframe = Frame(fontframe)
		rightframe.pack(side=RIGHT, padx=2, pady=2, fill=BOTH, expand=YES)

		frame = Frame(rightframe)
		frame.pack(side=TOP, padx=2, pady=2, fill=BOTH, expand=YES)

		#Label(frame, justify=LEFT, text='Font:').pack(side=TOP, anchor=W)
		self.listFontName = tkExtra.ExListbox(frame,
					height=5,
					takefocus=FALSE,
					exportselection=FALSE)
		self.listFontName.pack(side=LEFT, expand=TRUE, fill=BOTH)
		self.listFontName.bind('<ButtonRelease-1>',
					self.listFontButtonRelease)
		#self.listFontName.bind("<<ListboxSelect>>",
		#			self.listFontButtonRelease)
		self.listFontName.bind("<KeyRelease>",
					self.listFontKeyRelease)
		sb = Scrollbar(frame)
		sb.pack(side=LEFT, fill=Y)
		sb.config(command=self.listFontName.yview)
		self.listFontName.config(yscrollcommand=sb.set)

		# --- Param ---
		frame = Frame(rightframe)
		frame.pack(padx=5, pady=2, fill=X)
		Label(frame, text='Size:').pack(side=LEFT, anchor=W)
		self.menuFontSize = tkExtra.ExOptionMenu(frame,
					self.fontSize,
					command=self.setFontSample,
					*_FONT_SIZE)
		self.menuFontSize.pack(side=LEFT, anchor=W)

		frame = Frame(rightframe)
		frame.pack(padx=5, pady=2, fill=X)
		self.fontAbs.set(1)
		self.checkFontAbs = Checkbutton(frame, text='DPI',
					variable=self.fontAbs,
					onvalue=1, offvalue=0,
					command=self.setFontSample)
		self.checkFontAbs.pack(side=LEFT, anchor=W)

		self.checkFontBold = Checkbutton(frame, text='Bold',
					variable=self.fontBold,
					onvalue=1, offvalue=0,
					command=self.setFontSample)
		self.checkFontBold.pack(side=LEFT, anchor=W)

		self.checkFontItalic = Checkbutton(frame, text='Italic',
					variable=self.fontItalic,
					onvalue=1, offvalue=0,
					command=self.setFontSample)
		self.checkFontItalic.pack(side=LEFT, anchor=W)

		# --- Sample ---
		frame = Frame(rightframe, relief=SOLID, borderwidth=1)
		frame.pack(side=TOP, padx=5, pady=5, expand=TRUE, fill=BOTH)

		self.labelFontSample = Label(frame,
			text=   "AaBbCcDdEe\nFfGgHhIiJjK\n"
				"A\u03b1B\u03b2\u0393\u03b3\u0394\u03b4E\u03b5\n"
				"Z\u03b6H\u03b7\u0398\u03b8I\u03b9K\u03ba\n"
				"\u039b\u03bbM\u03bcN\u03bd\u039e\u03beOo\n"
				"\u03a0\u03c0P\u03c1\u03a3\u03c3T\u03c4Y\u03c5\u03a6\n"
				"\u03d5X\u03c7\u03a8\u03c8\u03a9\u03c9\n"
				"1234567890\n#:+=(){}[]",
			justify=LEFT, font=self.font)
		self.labelFontSample.pack(expand=TRUE, fill=BOTH)
		self.fontFamily.trace("w", self.fontFamilyChange)
		return fontframe

	# --------------------------------------------------------------------
	def loadFonts(self):
		fonts = list(tkinter.font.families(self))
		fonts.sort()
		# remove duplicates
		self.listFontName.insert(END,'')
		for font in fonts:
			self.listFontName.insert(END, font)
		self.listFontName.see(0)
		self.listFontName.activate(0)
		self.listFontName.selection_set(0)

	# --------------------------------------------------------------------
	def listFontButtonRelease(self, event=None):
		font = self.listFontName.get(ANCHOR)
		self.fontFamily.set(font)
		self.setFontSample()

	# --------------------------------------------------------------------
	def listFontKeyRelease(self, event=None):
		self.after_idle(self.listFontButtonRelease)

	# --------------------------------------------------------------------
	def setFontSample(self, event=None):
		self.fontChanged = True
		fontFamily = self.fontFamily.get()

		if self.fontAbs.get():
			fontSign = 1
		else:
			fontSign = -1

		if self.fontBold.get():
			fontWeight = tkinter.font.BOLD
		else:
			fontWeight = tkinter.font.NORMAL

		if self.fontItalic.get():
			fontSlant = tkinter.font.ITALIC
		else:
			fontSlant = tkinter.font.ROMAN

		try:
			self.font.config(size=fontSign*int(self.fontSize.get()),
				weight=fontWeight, slant=fontSlant,
				family=fontFamily)
		except UnicodeEncodeError:
			pass

	# --------------------------------------------------------------------
	def fontListChange(self, event=None):
		if self.fontCurrent and self.fontChanged:
			for font in self.fontCurrent:
				self.fontDict[font] = \
					(self.font.cget("family"), self.font.cget("size"),
					 self.font.cget("weight"), self.font.cget("slant"))

		sel = self.fontList.curselection()
		self.fontCurrent = list(map(self.fontList.get, sel))
		if self.fontCurrent:
			info = self.fontDict.get(self.fontCurrent[0])
			family = info[0]
			self.fontFamily.set(family)
			sz = int(info[1])
			if sz>0:
				self.fontAbs.set(1)
			else:
				self.fontAbs.set(0)
			self.fontSize.set(abs(sz))
			self.fontBold.set(info[2] == tkinter.font.BOLD)
			self.fontItalic.set(info[3] == tkinter.font.ITALIC)
		else:
			family = ""
			self.fontFamily.set("")
			self.fontAbs.set(1)
			self.fontSize.set(12)
			self.fontBold.set(0)
			self.fontItalic.set(0)

		self.listFontName.selection_clear(0, END)
		for i in range(self.listFontName.size()):
			if self.listFontName.get(i) == family:
				self.listFontName.selection_set(i)
				self.listFontName.selection_anchor(i)
				self.listFontName.activate(i)
				self.listFontName.see(i)
				break
		self.setFontSample()
		self.fontChanged = False

	# --------------------------------------------------------------------
	def fontFamilyChange(self, a=None, b=None, c=None):
		if self.fontFamily.get()=="":
			state = DISABLED
		else:
			state = NORMAL

		self.menuFontSize["state"] = state
		self.checkFontAbs["state"] = state
		self.checkFontBold["state"] = state
		self.checkFontItalic["state"] = state

	# --------------------------------------------------------------------
	def colorsFrame(self):
		pageframe = Frame(self.mainframe)

		self.colorsList = tkExtra.MultiListbox(pageframe,
					(("Item", 20, None),
					 ("Color", 15, None)))
		self.colorsList.pack(expand=YES, fill=BOTH)
		self.colorsList.setPopupMenu(
				[("Edit",   0, self.editColor,tkFlair.icons["edit"])])
#				 ("Add",    0, self.addColor, tkFlair.icons["add"]),
#				 ("Delete", 0, self.delColor, tkFlair.icons["x"])])

		self.colorsList.bindList('<ButtonRelease-1>', self.editColor)
		self.colorsList.sortAssist = None

		for name in _COLOR_LIST:
			color = str(self.colorDict.get(name,""))
			self.colorsList.insert(END, (name,color))
			try:
				self.colorsList.lists[1].itemconfig(END, background=color)
			except:
				self.colorDict[name] = ""

		return pageframe

	# ----------------------------------------------------------------------
	# Edit color
	# ----------------------------------------------------------------------
	def editColor(self, event=None):
		listbox = self.focus_get()
		if listbox is None: return

		if event is not None:
			active = self.colorsList.nearest(event.y)
			if active is not None:
				self.colorsList.activate(active)
			else:
				return
		else:
			active = self.colorsList.index(ACTIVE)

		if len(self.colorsList.curselection())>1: return

		values = self.colorsList.get(active)
		lid    = self.colorsList.lists.index(listbox)

		# Edit and change value
		if lid == 1:
			name = self.colorsList.lists[0].get(active)
			edit = tkExtra.InPlaceColor(listbox)

			if edit.value is None or edit.value=="":
				self.colorsList.lists[lid].set(active, "")
			else:
				self.colorsList.lists[lid].set(active, edit.value)
				if edit.value!="":
					self.colorsList.lists[lid].itemconfig(active,
						background=edit.value)
					self.colorDict[name] = str(edit.value)

	# --------------------------------------------------------------------
	def _setBool(self, var, section, name):
		try: var.set(self.config.getboolean(section, name))
		except configparser.NoOptionError: pass

	# --------------------------------------------------------------------
	def _setStr(self, var, section, name):
		var.delete(0,END)
		try: var.insert(0, self.config.get(section, name))
		except: pass

	# --------------------------------------------------------------------
	def _setInt(self, var, section, name):
		try: var.set(0, self.config.get(section, name))
		except: pass

	# --------------------------------------------------------------------
	def init(self):
		section = tkFlair._FLAIR_SECTION
		self._setBool(self.skipintro,	section, "skipintro")
		self._setBool(self.warnmultmat,	section, "warnmultmat")
		self._setBool(self.flukafiles,	section, "showflukafiles")
		self._setStr(self.autosave,	section, "autosave")
		#self._setStr(self.keytime,	section, "keytimethreshold")
		#self._setStr(self.balloondelay,	section, "balloondelay")
		self._setStr(self.timeformat,	section, "timeformat")
		self.iconZoom.set(self.config.get(section,"iconzoom"))

		section = tkFlair._PROJECT_SECTION
		self._setBool(self.keepbackup, section, "keepbackup")
		self._setBool(self.cleanup,    section, "cleanup")
		self._setStr(self.flukavar,    section, "flukavar")
		self._setStr(self.flukadir,    section, "flukadir")
		try: self.flufor.set(self.config.get(section, "flufor"))
		except: pass

		self._setStr(self.rfluka,      section, "rfluka")
		self._setStr(self.flukaexe,    section, "flukaexe")
		self._setStr(self.fff,         section, "fff")
		self._setStr(self.usbsuw,      section, "usbsuw")
		self._setStr(self.ustsuw,      section, "ustsuw")
		self._setStr(self.usxsuw,      section, "usxsuw")
		self._setStr(self.usrsuw,      section, "usrsuw")
		self._setStr(self.usysuw,      section, "usysuw")
		self._setStr(self.detsuw,      section, "detsuw")
		self._setStr(self.usbmax,      section, "usbmax")
		self._setStr(self.gplevbin,    section, "gplevbin")
		self._setStr(self.usbrea,      section, "usbrea")

		self._setStr(self.killCmd,     section, "kill")

		self._setStr(self.editor,      section, "editor")
		self._setStr(self.terminal,    section, "terminal")
		self._setStr(self.debugger,    section, "debugger")
		self._setStr(self.tmpprefix,   section, "tmpprefix")
		self._setStr(self.timethreshold,section,  "timethreshold")
		self._setStr(self.refreshinterval,section,"refreshinterval")

		# Populate gnuplot commands
		section = tkFlair._PLOT_SECTION
		for k,v in self.config.items(section):
			if k[0] == ".": self.gnuplotFileTypes.insert(END, (k,v))
		self.gnuplotFileTypes.sort(0)
		self.gnuplotCommands.insert("0.0",
			self.config.get(section, "commands"))
		self._setStr(self.gnuplotTerminal,    section, "terminal")
		self._setStr(self.gnuplot,    section, "program")

		# Populate Data
		section = tkFlair._DATA_SECTION
		for k,v in self.config.items(section):
			if k.startswith("rule."): continue
			type = Project._usr_type.get(k)
			self.usrTypes.insert(END, (type, v))
		self.usrTypes.sort(0)

		for op,r in Project._usr_rule:
			self.usrRules.insert(END,"%s%s"%(op,r))

		#Add custom batch system script
		section = tkFlair._BATCH_SECTION
		self._setStr(self.spawnRunName,    section, "spawnname")
		for k,v in self.config.items(section):
			if k=="spawnname": continue
			self.queueTypes.insert(END, (k,v))
		self.queueTypes.sort(0)

		# Geometry Section
		section = tkFlair._GEOMETRY_SECTION
		self._setBool(self.laptopMode,  section, "laptop")
		self._setBool(self.invertWheel, section, "invertwheel")
		self._setBool(self.pinZoom,     section, "pinzoom")
		self._setBool(self.useQUA,      section, "usequa")
		if Input._developer:
			self._setBool(self.useBOX,      section, "usebox")
		self.cores.set(self.config.get( section, "cores"))
		self.snap.set(self.config.get(  section, "snap"))
		self._setStr(self.zero,		section, "zero")
		self._setStr(self.infinite,	section, "infinite")
		self._setInt(self.accuracy,	section, "accuracy")

		# Fonts
		self.fontCurrent = []
		self.fontList.selection_clear(0,END)
		self.fontList.selection_set(0)
		self.fontListChange()

	# --------------------------------------------------------------------
	# Ok button pressed
	# --------------------------------------------------------------------
	def ok(self, event=None):
		self.fontListChange()

		# Copy all selected fonts to config
		section = tkFlair._FONT_SECTION
		for n,v in list(self.fontDict.items()):
			if len(v[0]) == 0:	# Delete font
				try:
					self.config.remove_option(section, n)
				except: pass
			else:
				tkFlair.setFont(n, v)

		# Copy all selected colors
		section = tkFlair._COLOR_SECTION
		for n,v in list(self.colorDict.items()):
			if v is not None:
				tkFlair.setColor(n, str(v))
			else:
				self.config.remove_option(section, n)

		# Get General options
		section = tkFlair._FLAIR_SECTION
		self.config.set(section, "skipintro",       str(bool(self.skipintro.get())))
		self.config.set(section, "warnmultmat",     str(bool(self.warnmultmat.get())))
		self.config.set(section, "showflukafiles",  str(bool(self.flukafiles.get())))
		self.config.set(section, "autosave",        self.autosave.get())
		#self.config.set(section, "keytimethreshold",self.keytime.get())
		#self.config.set(section, "balloondelay",    self.balloondelay.get())
		self.config.set(section, "timeformat",      self.timeformat.get())
		self.config.set(section, "iconzoom",        self.iconZoom.get())

		section = tkFlair._PROJECT_SECTION
		self.config.set(section, "keepbackup",      str(bool(self.keepbackup.get())))
		self.config.set(section, "cleanup",         str(bool(self.cleanup.get())))
		self.config.set(section, "flukavar",        self.flukavar.get())
		self.config.set(section, "flukadir",        self.flukadir.get())
		self.config.set(section, "flufor",          self.flufor.get())

		self.config.set(section, "rfluka",          self.rfluka.get())
		self.config.set(section, "flukaexe",        self.flukaexe.get())
		self.config.set(section, "fff",             self.fff.get())
		self.config.set(section, "usbsuw",          self.usbsuw.get())
		self.config.set(section, "ustsuw",          self.ustsuw.get())
		self.config.set(section, "usxsuw",          self.usxsuw.get())
		self.config.set(section, "usrsuw",          self.usrsuw.get())
		self.config.set(section, "usysuw",          self.usysuw.get())
		self.config.set(section, "detsuw",          self.detsuw.get())
		self.config.set(section, "usbmax",          self.usbmax.get())
		self.config.set(section, "gplevbin",        self.gplevbin.get())
		self.config.set(section, "usbrea",          self.usbrea.get())

		self.config.set(section, "kill",            self.killCmd.get())
		self.config.set(section, "editor",          self.editor.get())
		self.config.set(section, "terminal",        self.terminal.get())
		self.config.set(section, "debugger",        self.debugger.get())
		self.config.set(section, "tmpprefix",       self.tmpprefix.get())
		self.config.set(section, "timethreshold",   self.timethreshold.get())
		self.config.set(section, "refreshinterval", self.refreshinterval.get())

		# Gnuplot
		section = tkFlair._PLOT_SECTION
		self.config.remove_section(section)
		self.config.add_section(section)
		for name,value in self.gnuplotFileTypes.get(0,END):
			self.config.set(section, name, value)

		self.config.set(section, "commands",
			self.gnuplotCommands.get("1.0",END))
		self.config.set(section, "terminal", self.gnuplotTerminal.get())
		self.config.set(section, "program",  self.gnuplot.get())

		# Batch
		section = tkFlair._BATCH_SECTION
		self.config.remove_section(section)
		self.config.add_section(section)
		for name,value in self.queueTypes.get(0,END):
			# say("ConfigDialog:ok",t name, value)
			# queueName = Project._usr_type.get(name)
			self.config.set(section, name, value)
		self.config.set(section, "spawnname",       self.spawnRunName.get())

		# Data processing
		section = tkFlair._DATA_SECTION
		self.config.remove_section(section)
		self.config.add_section(section)
		for name,value in self.usrTypes.get(0,END):
			type = Project._usr_type.get(name)
			self.config.set(section, type, value)

		i = 1
		for rule in self.usrRules.get(0,END):
			self.config.set(section, "rule.%d"%(i), rule)
			i += 1

		# Geometry Section
		section = tkFlair._GEOMETRY_SECTION
		self.config.set(section, "laptop",	str(bool(self.laptopMode.get())))
		self.config.set(section, "invertwheel",	str(bool(self.invertWheel.get())))
		self.config.set(section, "pinzoom",	str(bool(self.pinZoom.get())))
		self.config.set(section, "usequa",	str(bool(self.useQUA.get())))
		if Input._developer:
			self.config.set(section, "usebox",	str(bool(self.useBOX.get())))
		self.config.set(section, "cores",	self.cores.get())
		self.config.set(section, "snap",	self.snap.get())
		self.config.set(section, "zero",	self.zero.get())
		self.config.set(section, "infinite",	self.infinite.get())
		self.config.set(section, "accuracy",	self.accuracy.get())

		# Destroy the dialog and exit
		self.rc = True
		self.close()

	# --------------------------------------------------------------------
	def cancel(self, event=None):
		self.rc = False
		self.close()

	# --------------------------------------------------------------------
	def close(self):
		ConfigDialog.width  = self.winfo_width()
		ConfigDialog.height = self.winfo_height()
		self.destroy()

	# --------------------------------------------------------------------
	def help(self, event=None):
		Manual.show("F6")

#-------------------------------------------------------------------------------
if __name__ == "__main__":
	tk = Tk()
	tkFlair.loadIcons()
	tkFlair.openIni()
	cd = ConfigDialog(tk, tkFlair.config)
	cd.show()
	tk.mainloop()
