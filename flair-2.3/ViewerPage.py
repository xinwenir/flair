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
#
# Author: Marco.Mauri@cern.ch
# Date:	1-Nov-2007

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.sala@mi.infn.it"

import os
import re
from stat import *

try:
	from tkinter import *
except ImportError:
	from tkinter import *
        
import tkinter as tk
import tkTree
import tkFlair
import tkExtra
import bFileDialog
import tkTextEditor
import tkFileViewer

import Ribbon
import FlairRibbon

#===============================================================================
# Generic TOC File Data
#===============================================================================
class TOCFileData(tkFileViewer.FileData):
	# ----------------------------------------------------------------------
	# Create Table of Contents
	# ----------------------------------------------------------------------
	def makeTOC(self):
		self.sections.clear()

		self.sections[self.check[0][0]] = (0, None)

		#Check each string....
		for section, start_end, offset, fingerprint, afterthis in self.check:
			pat = re.compile(fingerprint, re.MULTILINE) # | re.DOTALL)
			#By default start searching from the beginning..
			fromhere = 0
			if afterthis is not None:
				match = self.sections.get(afterthis, None)
				if match is not None:
					fromhere = match[0]

			match = pat.search(self.filedata, fromhere)
			if match is not None:
				#Ok the string matches. what shell we do now?
				startingpoint = match.start()
				#1) Advancing by a # of lines specified by "offset"
				if offset<0:
					for i in range(-offset):
						n = self.filedata.rfind("\n",startingpoint)
						if n<0: break
						startingpoint = n - 1
				elif offset>0:
					for i in range(offset):
						n = self.filedata.find("\n",startingpoint)
						if n<0: break
						startingpoint = n + 1

				if start_end:
					#2) If this is a True  create a section...
					self.sections[section] = (startingpoint, None)
				else:
					#3) If this is an False we find the proper section and close it.
					match = self.sections.get(section)
					if match is not None and match[1] is None:
						self.sections[section] = (match[0], startingpoint)

		#Now it's time to check if something went wrong
		match = self.sections.get(self.check[-1][0])
		if match is None or match[1] is None:
			self.sections["ERROR"] = (len(self.filedata)-1000, len(self.filedata))

		#Let's close all sections that are still open.
		for key,value in list(self.sections.items()):
			if value[1] is None:
				temp_end = len(self.filedata)
				for f,t in list(self.sections.values()):
					if f>value[0] and f<temp_end:
						temp_end = f
				self.sections[key] = (value[0], temp_end)

	# ----------------------------------------------------------------------
	# Return icon to display
	# ----------------------------------------------------------------------
	def icon(self, name):
		return tkFlair.icons.get(name,None)

#===============================================================================
# Output File Data
#===============================================================================
class OutputFileData(TOCFileData):
	STARTLINE = "1 *====================================================================*\n"

	# root is Linked to "Notes Page"
	OVSummaryTree = [("Totals/CPU time",		"cputime",	None),
			 ("# of stars",			"stars",	None),
			 ("# of secondaries in stars",	"2stars",	None),
			 ("# of fissions",		"fissions",	None),
			 ("# of decay products",	"decayprod",	None),
			 ("# of particles decayed",	"decpart",	None),
			 ("# of stopping particles",	"stoppart",	None),
			 ("# of part. from low en. neutrons","lowneut",	None),
			 ("Energy balance",		"enebalance",	None)]

	OVScoringTree = [("USRBIN",			"USRBIN",	None),
			 ("USRBDX",			"USRBDX",	None),
			 ("USRTRACK",			"USRTRACK",	None),
			 ("USRCOLL",			"USRCOLL",	None),
			 ("USRYIELD",			"USRYIELD",	None),
			 ("RESNUCLE",			"RESNUCLE",	None),
			 ("DETECT",			"DETECT",	None),
			 ("Unknown section",		"UNKNOWN",	None)]	 #Not implemented yet

	OVInputEchoTree = [("Body data",		"bodydata",	None),
			   ("Region data",		"REGION",	None),
			   ("Body echo",		"bodyecho",	None),
			   ("Region echo",		"regionecho",	None),
			   ("Coord. transformations",	"ROT-DEFI",	None)]

	OVTreeList = [	("License/version",		"license",	None),
			("Input Echo",			"input",	OVInputEchoTree),
			("Nuclear Data",		"nucleardata",	None),
			("Mulmix",			"mulmix",	None),
			("Products/Decays",		"isotopedata",	None),
			("Neutron",			"neutrondata",	None),
			("dp/dx",			"dpdx",		None),
			("Blank Common",		"emfmem",	None),
			("Media Parameters",		"media",	None),
			("EMF",				"emf",		None),
			("Multispot",                   "spotdata",     None),
			("Particles",			"partdata",	None),
			("Beam",			"beam",		None),
			("Particle Thresholds",		"thresholds",	None),
			("Termination Conditions",	"termcond",	None),
			("Mult. Coulomb Scattering",	"multcoulomb",	None),
			("EM Showers",			"EMShowers",	None),
			("Importances",			"importance",	None),
			("Scoring",			"scoring",	OVScoringTree),
			("Material",			"material",	None),
			("Regions"	,		"regsum",	None),
			("Initialization Time",		"inittime",	None),
			("Output During Transport",	"rundata",	None),
			("Events by Region",		"inforeg",	None),
			("Scattering Statistics",	"scatstat",	None),
			("Biasing Counters",		"biasing",	None),
			("Run Summary",			"summary",	OVSummaryTree),
			("ERROR",			"ERROR",	None)]

	# The following informations are used to parse the output ad 'split' it into different sections.
	#Fields meaning
	#		Section name	start/end, lines-offset,   fingerprint,    Start searching after this section
	#				start=True
	#				end=False

	checkStrings = [("license",	False,	7, r"^ *distributor of the code",				None),
			("input",	True,	8, r"^ *distributor of the code",				None),
			("input",	False,	0, r"^ Total time used for input reading",			None),
			("bodydata",	True,	0, r"^ {50}Body data",						"inputecho"),
			("REGION",	True,	0, r"^ {50}Input region data",					"inputecho"),
			("bodyecho",	True,	0, r"^ Interpreted body echo",					"inputecho"),
			("regionecho",	True,	0, r"^ Interpreted region echo",				"inputecho"),
			("regionecho",	False,	0, r"\*{3}",							"regionecho"),
			("nucleardata",	True,	0, r"^ *\*{3} Reading evaporation and nuclear data from unit",	None),
			("mulmix",	True,	0, r"^ *\*{4} Subroutine Mulmix",				None),
			("isotopedata",	True,	0, r"^ *\*{4} Isotope tabulation data start",			None),
			("neutrondata",	True,	0, r"^ *\*{4} Non analog absorption factors start",		None),
			("dpdx",	True,	0, r"^ *\*{3} dp/dx tab\. generated up to",			None),
			("emfmem",	True,	0, r"^ \*\*\* Blank common temporary cells from",		None),
			("media",	True,	0, r"^1 Quantities/Biasing associated with each media:",	None),
			("emf",		True,	0, r"^1 Correspondence of regions and EMF-FLUKA material numbers and names:",None),
			("spotdata",	True,	0, r"=== Spotbeam properties ===",                              None),
			("partdata",	True,	0, r"=== Output before the actual run - Particle properties: ===",None),
			("beam",	True,	0, r"^ *=== Output before the actual run - Beam properties ===",None),
			("thresholds",	True,	0, r"^ *=== Particle transport thresholds:",			None),
			("termcond",	True,	0, r"^ *=== Termination conditions: ===",			None),
			("multcoulomb",	True,	0, r"^ *=== Multiple Coulomb scattering: ===",			None),
			("EMShowers",	True,	0, r"^ *=== Electromagnetic Showers: ===",			None),
			("importance",	True,	0, r"^ *Region        Particle   importances",			None),
			("scoring",	True,	0, r'^ \*{7} "usrbin" option:',					None),
			("scoring",	False,	0, r"^  === Material compositions: ===",			"scoring"),
			("USRBIN",	True,	0, r'^ \*{7} "usrbin" option:',					"scoring"),
			("USRBDX",	True,	0, r'^ \*{7} "USRBDX" option:',					"USRBIN"),
			("USRTRACK",	True,	0, r'^ \*{7} "USRTRACK" option:',				"USRBIN"),
			("USRCOLL",	True,	0, r'^ \*{7} "USRCOLL"  option:',				"USRBIN"),
			("USRYIELD",	True,	0, r'^ \*{7} "Usryield" option:',				"USRBIN"),
			("RESNUCLE",	True,	0, r'^ \*{7} "RESNUCLEI" option:',				"USRBIN"),
			("DETECT",	True,	0, r'^ \*{7} "Detect" option:',					"USRBIN"),
			("ROT-DEFI",	True,	0, r'^ \*{3} Transformation #',					"USRBIN"),
			("material",	True,	0, r"^  === Material compositions: ===",			"DETECT"),
			("regsum",	True,	0, r"^  === Regions: materials and fields ===",			"DETECT"),
			("inittime",	True,	0, r"^  === End of the output associated with the input ===",	"DETECT"),
			("rundata",	True,	0, r"^1NUMBER OF BEAM *NUMBER OF BEAM",				"inittime"),
			("inforeg",	True,	0, r"^1Region # name *volume ",					"inittime"),
			("scatstat",	True,	0, r"^ \*{4} Total number of not-performed scatterings in FLUKA","inittime"),
			("biasing",     True,   0, r"^   Hadron importance RR/Splitting counters",		"inittime"),
			("summary",	True,	0, r"^ \* ====== Summary of FLUKA-",				"inittime"),
			("summary",	False,	4, r"^ \* ======== End of FLUKA-",				"summary"),
			("cputime",	True,	0, r"^ *Total number of primaries run:",			"summary"),
			("stars",	True,	0, r"^ Number of stars generated per beam particle",		"cputime"),
			("2stars",	True,	0, r"^ Number of secondaries generated in inelastic interactions","cputime"),
			("fissions",	True,	0, r"^ Number of high energy fissions generated",		"cputime"),
			("decayprod",	True,	0, r"^ Number of decay products produced",			"cputime"),
			("decpart",	True,	0, r"^ Number of particles decayed per beam particle",		"cputime"),
			("stoppart",	True,	0, r"^ Number of stopping particles",				"cputime"),
			("lowneut",	True,	0, r"^ Number of secondaries created by low energy neutron",	"cputime"),
			("enebalance",	True,	0, r"^ .* GeV available per beam particle divided into",	"cputime"),
			("end",		True,  -2, r"^ \* ======== End of FLUKA-",				"enebalance"),
			("end",		False,  4, r"^ \* ======== End of FLUKA-",				"enebalance")]

	# ----------------------------------------------------------------------
	def __init__(self, *kw):
		tkFileViewer.FileData.__init__(self,*kw)
		self.tree  = OutputFileData.OVTreeList
		self.check = OutputFileData.checkStrings

	# ----------------------------------------------------------------------
	# Create Table of Contents
	# ----------------------------------------------------------------------
	def makeTOC(self):
		# Check if it is a FLUKA output file
		if self.filedata[:len(OutputFileData.STARTLINE)] != OutputFileData.STARTLINE:
			return
		else:
			TOCFileData.makeTOC(self)

#===============================================================================
# Listing File Data
#===============================================================================
class ListingFileData(tkFileViewer.FileData):
	_detectorPattern = re.compile(r"^ ?#? Detector n:\s*(.*)\s*", re.MULTILINE)
	_regionPattern   = re.compile(r"^ *\w* +binning n\. +(\d+) *\"(\w*) *\" *,.*$", re.MULTILINE)
	_pattern = [	("detector", _detectorPattern, "DarkRed"),
			("detector", _regionPattern,   "DarkRed"),
			("resp",     re.compile(r"^     Tot\. resp\..*", re.MULTILINE), "Blue"),
			("resp",     re.compile(r"^     Tot\. response .*", re.MULTILINE), "Blue"),
			("resp",     re.compile(r"^     \( -->.*", re.MULTILINE), "Blue")]

	# ----------------------------------------------------------------------
	def __init__(self, *kw):
		tkFileViewer.FileData.__init__(self,*kw)

	# ----------------------------------------------------------------------
	# Create Table of Contents
	# ----------------------------------------------------------------------
	def makeTOC(self):
		self.sections.clear()

		name  = None
		start = -1
		stop  = -1
		while True:
			match = ListingFileData._detectorPattern.search(self.filedata, start+1)
			if match is None:
				match = ListingFileData._regionPattern.search(self.filedata, start+1)
				if match is None: break
				newname  = "%s %s"%(match.group(1),match.group(2))
			else:
				newname  = match.group(1)
			if name:
				self.tree.append((name,name,None))
				self.sections[name] = (start, match.start())
			name = newname
			start = match.start()

		if name is not None:
			self.tree.append((name,name,None))
			self.sections[name] = (start,len(self.filedata))

	# ----------------------------------------------------------------------
	# Return icon to display
	# ----------------------------------------------------------------------
	def icon(self, name):
		return tkFlair.icons.get(name,None)

	# ----------------------------------------------------------------------
	def highlight(self, widget):
		widget.tag_delete('*')

		text = widget.get("1.0",END)

		for tag, pat, color in ListingFileData._pattern:
			for m in pat.finditer(text):
				widget.tag_add(tag,	"1.0 + %d chars"%(m.start()),
							"1.0 + %d chars"%(m.end()))
			widget.tag_configure(tag, foreground=color)

#===============================================================================
# Input File Data
#===============================================================================
class InputFileData(tkFileViewer.FileData):
	_pattern = [	("comment",  re.compile(r"^[*!].*", re.MULTILINE), "input.comment"),
			("define",   re.compile(r"^#.*",    re.MULTILINE), "input.define"),
			("transform",re.compile(r"^\$.*",   re.MULTILINE), "input.define"),
			("tag",      re.compile(r"^[\w-]*", re.MULTILINE), "input.tag")]

	# ----------------------------------------------------------------------
	def highlight(self, widget):
		widget.tag_delete('*')

		text = widget.get("1.0",END)

		for tag, pat, color in InputFileData._pattern:
			for m in pat.finditer(text):
				widget.tag_add(tag,	"1.0 + %d chars"%(m.start()),
							"1.0 + %d chars"%(m.end()))
			color = tkFlair.getColor(color,None)
			if color: widget.tag_configure(tag, foreground=color)

#===============================================================================
# Output File Viewer
#===============================================================================
class OutFileViewer(tkFileViewer.FileViewer):
	def __init__(self, master, **kw):
		"""Open a viewer assuming that the file exists."""
		tkFileViewer.FileViewer.__init__(self, master, **kw)
#		self.iconbitmap("@%s/icons/view.xbm"%(tkFlair.prgDir))

	# ----------------------------------------------------------------------
	# Routine to be overwritten if loading has to be done with another
	# FileData class
	# ----------------------------------------------------------------------
	def loadData(self, filename, maxsize):
		fn,ext = os.path.splitext(filename)
		if ext == ".out":
			curfile = OutputFileData(filename, maxsize)
		elif ext == ".lis":
			curfile = ListingFileData(filename, maxsize)
		elif ext == ".inp":
			curfile = InputFileData(filename, maxsize)
		else:
			curfile = tkFileViewer.FileData(filename, maxsize)

		if curfile.load(): return None
		curfile.makeTOC()
		return curfile

#===============================================================================
# Viewer Tab Page
#===============================================================================
class ViewerPage(FlairRibbon.FlairPage):
	"""Flair file viewer diving in sections the output files"""

	_name_ = "Viewer"
	_icon_ = "view"

	#----------------------------------------------------------------------
	def init(self):
		self.files = []

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Find ===========
		group = Ribbon.LabelGroup(self.ribbon, "Search")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		self.findString = tkExtra.LabelEntry(group.frame,
					"Search",
					"DarkGray",
					background="White", #Ribbon._BACKGROUND,
					width=20)
		self.findString.pack(side=LEFT)
		tkExtra.Balloon.set(self.findString , "Find cards with string")
		self.findString.bind("<Return>",   self.find)
		self.findString.bind("<KP_Enter>", self.find)

		self.findCase = BooleanVar()
		self.findCase.set(False)
		b = Checkbutton(group.frame,
				variable=self.findCase,
				padx=0,
				pady=0,
				borderwidth=0,
				highlightthickness=0,
				font=Ribbon._FONT,
				activebackground="LightYellow",
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b, "Case sensitive")

		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["find"],
				command=self.find,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b, "Find cards with string")

		# ========== File ===========
		group = Ribbon.LabelGroup(self.ribbon, "File")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row=0,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["open32"],
				command=self.load,
				text="Load",
				compound=TOP,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, stick=NSEW)
		tkExtra.Balloon.set(b, "Load file")

		# ---
		col,row=1,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["reload32"],
				text="Reload",
				compound=TOP,
				command=self.reload,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, stick=NSEW)
		tkExtra.Balloon.set(b, "Reload selected file")

		# ---
		col,row=2,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["edit32"],
				text="Edit",
				compound=TOP,
				command=self.editor,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, stick=NSEW)
		tkExtra.Balloon.set(b, "Edit selected file")

		# ---
		col,row=3,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["print32"],
				text="Print",
				compound=TOP,
				command=self.hardcopy,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, stick=NSEW)
		tkExtra.Balloon.set(b, "Print selected file")

		# ---
		col,row=4,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["x32"],
				text="Close",
				compound=TOP,
				command=self.closeFile,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, stick=NSEW)
		tkExtra.Balloon.set(b, "Close selected file")

	#----------------------------------------------------------------------
	# Create Project page
	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self, False)
		self.viewer = OutFileViewer(self.frame)
		self.viewer.editor.txt["background"] = "White"
		self.viewer.pack(expand=YES, fill=BOTH)
		self.viewer.editor.txt.bind("<Escape>",        self.closeFile)
		self.viewer.editor.txt.bind("<Control-Key-f>", lambda e,s=self : s.findString.focus_set())
		self.viewer.tree.bind("<Escape>",              self.closeFile)

	#----------------------------------------------------------------------
	def load(self, filename=None):
		if filename is None:
			filename = bFileDialog.askopenfilename(master=self.page)
			if filename is None: return
		filename = self.project.relativePath(filename)
		self.viewer.load(filename)
		self.files = [x.filename for x in self.viewer.files if x]

	#----------------------------------------------------------------------
	def reload(self, filename=None):
		if filename is None:
			self.viewer.reloadSelected()
		else:
			filename = self.project.relativePath(filename)
			self.viewer.reload(filename)

	#----------------------------------------------------------------------
	def show(self, name, text):
		self.viewer.show(name, text)

	#----------------------------------------------------------------------
	def hasFile(self, filename):
		filename = self.project.relativePath(filename)
		return self.viewer.hasFile(filename)

	#----------------------------------------------------------------------
	def editor(self):
		tkFlair.editor([x.filename for x in self.viewer.files if x])

	#----------------------------------------------------------------------
	def closeFile(self, event=None):
		self.viewer.close()
		self.files = [x.filename for x in self.viewer.files if x]
		if not self.files:
			self.viewer.clear()
			self.hide(False)

	#----------------------------------------------------------------------
	def afterDock(self):
		for f in self.files:
			self.load(f)

	# ----------------------------------------------------------------------
	def hardcopy(self, event=None):
		printer = tkFlair.printer(self.page)
		if printer:
			printer.open()
			printer.write(self.viewer.editor.txt.get("1.0",END))
			printer.close()

	#----------------------------------------------------------------------
	def find(self, event=None):
		txt = self.findString.get()
		if not txt: return
		self.viewer.editor.find(txt, not self.findCase.get())

#===============================================================================
# Toplevel File Viewer
#===============================================================================
class TopFileViewer(Toplevel):
	width  = -1
	height = -1

	def __init__(self, master, **kw):
		"""Open a viewer assuming that the file exists."""
		Toplevel.__init__(self, master, class_="TopFileViewer", **kw)

		self.title("FileViewer")

		self.viewer = OutFileViewer(self)
		self.viewer.pack(expand=TRUE, fill=BOTH)

		# Menu
		menu = Menu(self)
		self.config(menu=menu)

		# File Menu
		filemenu = Menu(menu)
		menu.add_cascade(label="File", underline=0, menu=filemenu)
		filemenu.add_command(label="Open", underline=0,
			command=self.viewer.openFile, accelerator="Ctrl-O")
		filemenu.add_command(label="Reload", underline=0,
			command=self.viewer.reloadSelected, accelerator="Ctrl-R")
		filemenu.add_command(label="Close", underline=0,
			command=self.close, accelerator="Ctrl-W")
		filemenu.add_command(label="Print", underline=0,
			command=self.printer,
			accelerator="Ctrl-P")
		filemenu.add_separator()
		filemenu.add_command(label="Quit", underline=0,
			command=self.quit, accelerator="Ctrl-Q")

		# Edit Menu
		editmenu = Menu(menu)
		menu.add_cascade(label="Edit", underline=0, menu=editmenu)
		editmenu.add_command(label="Copy", underline=0,
			command=lambda v=self.viewer.editor:v.event_generate("<<Copy>>"),
			accelerator="Ctrl-C")
		editmenu.add_separator()
		editmenu.add_command(label="Find", underline=0,
			command=self.viewer.editor.findDialog,
			accelerator="Ctrl-F")
		editmenu.add_command(label="Find Next", underline=5,
			command=self.viewer.editor.findNext,
			accelerator="Ctrl-G")

		viewmenu = Menu(menu)
		menu.add_cascade(label="View", underline=0, menu=viewmenu)
		viewmenu.add_command(label="Toggle TOC", underline=7,
				accelerator="F11",
				command=self.viewer.toggleTree)

		# Bindings
		self.bind('<Control-Key-o>', self.viewer.openFile)
		self.bind('<Control-Key-p>', self.printer)
		self.bind('<Control-Key-q>', self.quit)
		self.bind('<Control-Key-r>', self.viewer.reloadSelected)
		self.bind('<Control-Key-w>', self.close)
		self.bind('<F11>', self.viewer.toggleTree)
		self.protocol("WM_DELETE_WINDOW", self.quit)

		self.bind('<Control-Key-f>', self.viewer.editor.findDialog)
		self.bind('<Key-slash>',     self.viewer.editor.findDialog)
		self.bind('<Control-Key-g>', self.viewer.editor.findNext)

		self.deiconify()
		if TopFileViewer.width > 0:
			self.geometry("%dx%d"%(TopFileViewer.width, TopFileViewer.height))
		else:
			self.geometry("800x480")
		self.wait_visibility()
		self.lift()

	# ----------------------------------------------------------------------
	def load(self, filename, maxsize=None):
		self.viewer.load(filename, maxsize)

	# ----------------------------------------------------------------------
	def close(self, event=None):
		if self.viewer.close():
			self.quit()

	# ----------------------------------------------------------------------
	def quit(self, event=None):
		# Save window info on config
		TopFileViewer.width  = self.winfo_width()
		TopFileViewer.height = self.winfo_height()
		self.destroy()
		tk.destroy()

	# ----------------------------------------------------------------------
	def printer(self, event=None):
		printer = tkDialogs.Printer(self)
		if printer.show():
			printer.open()
			printer.write(self.editor.txt.get("1.0",END))
			printer.close()

#===============================================================================
if __name__ == "__main__":
	tk = Tk()
	tk.withdraw()

	# Set fonts if needd
	tkFlair.init()

	view = TopFileViewer(tk)
	view.bind("<Control-Key-q>", lambda e: tk.destroy())
	view.protocol("WM_DELETE_WINDOW", lambda : tk.destroy())
	view.geometry("850x700")
	try:
		tk.tk.call('wm','iconphoto',view._w, tkFlair.icons["view"])
	except TclError:
		tk.iconbitmap("@%s/icons/view.xbm"%(tkFlair.prgDir))
	for fn in sys.argv[1:]:
		view.load(fn)
	tk.mainloop()
	tkFlair.fini()
