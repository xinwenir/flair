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
# Date:	4-Sep-2006

# __version__ string from here modifies all VERSION/RELEASE info in specs and geoviewer
__author__     = "Vasilis Vlachoudis"
__email__      = "Paola.Sala@mi.infn.it"
__version__    = "2.3-0"
__credits__    = "David Sinuela," \
		 "Eleftherios Skordis," \
		 "Pablo Garcia Ortega," \
		 "Wioletta Kozlowska"

__name__       = "flair"
__created__    = "Jun 2006"
__www__        = "http://www.fluka.org/flair"
__revision__   = "4cf6d2c"
__lastchange__ = "Fri Apr 28 10:14:43 2017 +0200"

import os
import re
import sys
import time
import pickle
import random
import string
import traceback
from stat import *
from log import say
try:
	import configparser
except ImportError:
	import configparser as ConfigParser
import webbrowser

import gettext
try:
	import builtins
except:
	import builtins as __builtin__
from operator import attrgetter
import tkinter as tk
import tkinter.font as tkFont

# --- Global Variables ---
prgDir      = os.path.dirname(__file__)+os.sep
iniDir      = os.path.expanduser("~/.%s" % (__name__))
iniFile     = os.path.join(iniDir, "%s.ini"%(__name__))
iniMaterial = os.path.join(iniDir, "material.ini")
STAT_FILE   = os.path.join(iniDir, "stats.pickle")

# dirty way of substituting the "_" on the builtin namespace
#__builtin__.__dict__["_"] = gettext.translation('bCNC', 'locale', fallback=True).ugettext
builtins._ = gettext.translation('bCNC', os.path.join(prgDir,'locale'), fallback=True).gettext
builtins.N_ = lambda message: message

import Input
import Project
import Gnuplot
import RichText

try:
	from tkinter import *
	import tkinter.font
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	import tkinter.font as tkFont
	import tkinter.messagebox as messagebox

#import Palette
import tkExtra
import tkDialogs
import bFileDialog

# Variables
config           = None
#_lastFrame      = None
_SKIP_INTRO      = False
_SAVE_FRAME      = False
_errorReport     = True
_showFlukaFiles  = False
_gnuplotTerminal = ""
_maxRecent       = 10
_showToolbar     = True
_showStatusbar   = True
_output          = ""			# Global output
_autoSaveTime    = 120			# auto save to input.pickle every # seconds
_iconZoom        = 1

# --- Sections ---
_BATCH_SECTION    = Project.RunInfo.section
_COLOR_SECTION    = "Color"
_DATA_SECTION     = Project.Command.section
_FLAIR_SECTION    = "flair"
_FONT_SECTION     = "Font"
_GEOMETRY_SECTION = "Geometry"
_INPUT_SECTION    = "Input"
_LAYERS_SECTION   = "Layers"
_MATERIAL_SECTION = Project.MaterialInfo.section
_PAGE_SECTION     = "Page"
_PLOT_SECTION     = "Plot"
_PROJECT_SECTION  = Project.section
_STYLE_SECTION    = "Style"

__sections = [	_BATCH_SECTION,
		_COLOR_SECTION,
		_DATA_SECTION,
		_FLAIR_SECTION,
		_FONT_SECTION,
		_GEOMETRY_SECTION,
		_INPUT_SECTION,
		_LAYERS_SECTION,
		_MATERIAL_SECTION,
		_PAGE_SECTION,
		_PLOT_SECTION,
		_PROJECT_SECTION,
		_STYLE_SECTION]

# --- Global Fonts ---
_TREE_FONT  = ("helvetica", "-12")
_TITLE_FONT = ("helvetica", "-12", "bold")

# --- Global Colors ---
_TITLE_FOREGROUND_COLOR = None
_TITLE_BACKACTIVE_COLOR = None
_TITLE_BACKINACTIVE_COLOR = None

_ELABEL_FOREGROUND_COLOR = None
_ELABEL_BACKGROUND_COLOR = None

_ILABEL_FOREGROUND_COLOR = None
_ILABEL_BACKGROUND_COLOR = None

_CLABEL_FOREGROUND_COLOR = None
_CLABEL_BACKGROUND_COLOR = None

# Special tags
_ALL_GROUP          = "*all*"
_GEOMETRY_GROUP     = "Geometry"
_PREPROCESSOR_GROUP = "Preprocessor"
_BODIES_SUBGROUP    = "...Bodies"
_groups_order = ["General",
		"Primary",
		_GEOMETRY_GROUP,
		"Media",
		"Physics",
		"Transport",
		"Biasing",
		"Scoring",
		"Flair",
		"Developers",
		_PREPROCESSOR_GROUP]
_groups_underline = []

# Special card
_EXPANSION  = "expansion"
_TRANSLAT   = "translat"
_TRANSFORM  = "transform"

_IFENDIF   = "#if..#endif"

# Arrays
icons  = {}
errors = []
stats  = {}

# Plot formats types .eps, .ps, ... read from .ini
_PLOT_FORMAT   = {}
_PLOT_COMMANDS = ""

# Flair clipboard pattern
_FLAIRCLIP     = "<flair:pickle>"

NOTIFY_NORMAL  = 0
NOTIFY_WARNING = 1
NOTIFY_ERROR   = 2

_TEMPLATES_DIR = "templates"
_DEFAULT_INP   = "basic.inp"

#-------------------------------------------------------------------------------
def write(msg):
	global _output
	_output += msg + "\n"
	#if Input._developer:
	if msg.startswith("Error:"):
		sys.stderr.write(tkExtra.ANSI_BOLD)
		sys.stderr.write(msg)
		sys.stderr.write(tkExtra.ANSI_NORMAL)
	else:
		sys.stdout.write(msg)

#-------------------------------------------------------------------------------
def loadIcons():
	global icons
	if len(icons) > 0: return

	pattern = "%sicons%s%%s.gif" % (prgDir, os.sep)

	# Cards icons
	for name in list(Input.CardInfo._db.keys()):
		try:
			icons[name] = PhotoImage(file=pattern%(name))
		except TclError:
			pass

	# icons directory
	try:
		f = open(os.path.join(prgDir,"icons","icons.dir"),"r")
	except:
		return

	for line in f:
		line = line.strip()
		if not line: continue
		if line[0] == '#': continue
		line = line.split()
		name = line[0]
		try:
			icons[name] = PhotoImage(file=pattern%(name))
			for aux in line[1:]:
				icons[aux] = icons[name]
		except TclError:
			pass
	f.close()

	if _iconZoom > 1:
		for name,icon in list(icons.items()):
			icons[name] = icon.zoom(_iconZoom)

#-------------------------------------------------------------------------------
# FIXME doesn't work?
#-------------------------------------------------------------------------------
def __getitem__(name):
	return icons[name]

#-------------------------------------------------------------------------------
def delIcons():
	global icons
	if len(icons) > 0:
		for i in list(icons.values()):
			del i
		icons = {}

#-------------------------------------------------------------------------------
# Create a font string
#-------------------------------------------------------------------------------
def fontString(font):
	name  = str(font[0])
	size  = str(font[1])
	if name.find(' ')>=0:
		s = '"%s" %s'%(name,size)
	else:
		s = '%s %s'%(name,size)

	try:
		if font[2] == tkinter.font.BOLD: s += " bold"
	except: pass
	try:
		if font[3] == tkinter.font.ITALIC: s += " italic"
	except: pass
	return s

#-------------------------------------------------------------------------------
# Return a font from a string
#-------------------------------------------------------------------------------
def makeFont(name, value=None, root=None):
	try:
		font = tkinter.font.Font(root=root, name=name, exists=True)
	except TclError:
		font = tkinter.font.Font(root=root, name=name)
		font.delete_font = False
	except AttributeError:
		return None

	if value is None: return font

	if isinstance(value, tuple):
		font.configure(family=value[0], size=value[1])
		try:
			font.configure(weight=value[2])
		except:
			pass
		try:
			font.configure(slant=value[3])
		except:
			pass
	return font

#-------------------------------------------------------------------------------
# Get font from configuration
#-------------------------------------------------------------------------------
def getFont(name, default=None, root=None):
	try:
		value = config.get(_FONT_SECTION, name)
	except:
		font = makeFont(name, default, root)
		setFont(name, font)
		return font

	if isinstance(value, str):
		value = tuple(value.split(','))

	if isinstance(value, tuple):
		font = makeFont(name, value, root)
		if font is not None: return font
	return value

#-------------------------------------------------------------------------------
# Set font in configuration
#-------------------------------------------------------------------------------
def setFont(name, font):
	if font is None: return
	if isinstance(font,str):
		config.set(_FONT_SECTION, name, font)
	elif isinstance(font,tuple):
		config.set(_FONT_SECTION, name, ",".join(map(str,font)))
	else:
		config.set(_FONT_SECTION, name, "%s,%s,%s" % \
			(font.cget("family"),font.cget("size"),font.cget("weight")))

#-------------------------------------------------------------------------------
def getColor(name, default):
	try:
		return config.get(_COLOR_SECTION, name)
	except:
		setColor(name, default)
		return default

#-------------------------------------------------------------------------------
def setColor(name, value):
	config.set(_COLOR_SECTION, name, value)

#-------------------------------------------------------------------------------
# return color as Integer
#-------------------------------------------------------------------------------
def getIColor(name, default):
	return int(str(getColor(name,"#%06X"%(default)))[1:],16)

#-------------------------------------------------------------------------------
# set color from integer
#-------------------------------------------------------------------------------
def setIColor(name, value):
	setColor(name, "#%06X"%(value))

#-------------------------------------------------------------------------------
def setBool(section, name, value):
	config.set(section, name, str(bool(value)))

#-------------------------------------------------------------------------------
def getInt(section, name, default):
	try:
		return config.getint(section, name)
	except:
		config.set(section, name, str(default))
		return default

#-------------------------------------------------------------------------------
def getFloat(section, name, default):
	try:
		return config.getfloat(section, name)
	except:
		config.set(section, name, str(default))
		return default

#-------------------------------------------------------------------------------
def getStr(section, name, default):
	try:
		return config.get(section, name)
	except:
		config.set(section, name, str(default))
		return default

#-------------------------------------------------------------------------------
def getBool(section, name, default):
	try:
		return config.getboolean(section, name)
	except:
		config.set(section, name, str(bool(default)))
		return default

#-------------------------------------------------------------------------------
# Open and Read Initialization file
#-------------------------------------------------------------------------------
def openIni(ini=None, loaduserini=True):
	global config, iniDir, iniFile

	systemIni = os.path.join(prgDir, "%s.ini"%(__name__))

	if ini is not None:
		iniFile = os.path.abspath(ini)
		iniDir  = os.path.dirname(iniFile)

	config = configparser.RawConfigParser()
	config.optionxform=str	# preserve case
	if loaduserini:
		config.read([systemIni, iniFile])
	else:
		config.read(systemIni)

	# XXX XXX XXX Check for old configuration file....
	try:
		ver = config.get(__name__, "version")
		# Check version if needed
	except:
		# Patch... XXX XXX XXX XXX
		del config
		config = configparser.RawConfigParser()
		config.optionxform=str
		config.read(systemIni)
	# XXX XXX XXX XXX

	# Add the necessary sections
	for s in __sections:
		try: config.add_section(s)
		except: pass

#-------------------------------------------------------------------------------
# Read Initialization file
#-------------------------------------------------------------------------------
def readIni():
	global _SKIP_INTRO, _SAVE_FRAME
	global _TREE_FONT, _TITLE_FONT
	global _TITLE_FOREGROUND_COLOR, _TITLE_BACKACTIVE_COLOR, _TITLE_BACKINACTIVE_COLOR
	global _ELABEL_FOREGROUND_COLOR, _ELABEL_BACKGROUND_COLOR
	global _ILABEL_FOREGROUND_COLOR, _ILABEL_BACKGROUND_COLOR
	global _CLABEL_FOREGROUND_COLOR, _CLABEL_BACKGROUND_COLOR
	global _PLOT_COMMANDS
	global _showFlukaFiles, _errorReport, _lastFrame
	global _gnuplotTerminal, _autoSaveTime
	global _showToolbar, _showStatusbar
	global _iconZoom

	# Read basic fonts from the font section
	_TREE_FONT  = getFont("tree",  _TREE_FONT)
	_TITLE_FONT = getFont("title", _TITLE_FONT)

	import Ribbon
	Ribbon._TABFONT = getFont("ribbontab", Ribbon._TABFONT)
	Ribbon._FONT    = getFont("ribbon",    Ribbon._FONT)

	# read Tk fonts to initialize them
	font = getFont("TkDefaultFont")
	font = getFont("TkFixedFont")
	font = getFont("TkMenuFont")
	font = getFont("TkTextFont")

	_TITLE_FOREGROUND_COLOR  = getColor("title.foreground", "white")
	_TITLE_BACKACTIVE_COLOR  = getColor("title.background", "dark slate blue")
	_TITLE_BACKINACTIVE_COLOR= getColor("title.backinactive", "darkgray")

	_ELABEL_FOREGROUND_COLOR = getColor("entrylabel.foreground", "darkblue")
	_ELABEL_BACKGROUND_COLOR = getColor("entrylabel.background", "white")

	_ILABEL_FOREGROUND_COLOR = getColor("infolabel.foreground", "darkblue")
	_ILABEL_BACKGROUND_COLOR = getColor("infolabel.background", "lightgray")

	_CLABEL_FOREGROUND_COLOR = getColor("calclabel.foreground", "darkred")
	_CLABEL_BACKGROUND_COLOR = getColor("calclabel.background", "lightgray")

	tkExtra.Balloon.font = getFont("balloon", tkExtra.Balloon.font)
	tkExtra.Balloon.foreground = getColor("balloon.foreground", tkExtra.Balloon.foreground)
	tkExtra.Balloon.background = getColor("balloon.background", tkExtra.Balloon.background)

	# Correct materials
	for name, color in config.items(_COLOR_SECTION):
		# Correct material location
		if name.startswith("material."):
			config.set(_MATERIAL_SECTION, name.replace("material.","color."), color)
			config.remove_option(_COLOR_SECTION, name)
		elif name.startswith("alpha."):
			config.set(_MATERIAL_SECTION, name, color)
			config.remove_option(_COLOR_SECTION, name)

	# Read global colors
	for name, color in config.items(_COLOR_SECTION):
		if name.startswith("file."):
			bFileDialog.COLORS[name[5:]] = color

	# Flair
	section = _FLAIR_SECTION
	_showFlukaFiles  = getBool(section, "showflukafiles", _showFlukaFiles)
	_errorReport     = getBool(section, "errorreport",    _errorReport)
	_showToolbar     = getBool(section, "showtoolbar",    _showToolbar)
	_showStatusbar   = getBool(section, "showstatusbar",  _showStatusbar)
	if "-d" not in sys.argv and "-D" not in sys.argv:
		Input._developer = getBool(section, "developer",      Input._developer)
	_iconZoom        = getInt(section,  "iconzoom",       _iconZoom)

	bFileDialog._TIME_FORMAT  = getStr( section, "timeformat",bFileDialog._TIME_FORMAT)

	# Project
	Project.loadConfig(config)
	# Load materials
	Project.MaterialInfo.loadConfig(config)
	# Load all commands
	Project.Command.loadConfig(config)
	Project.UsrInfo.loadConfig(config)
	Project.RunInfo.loadConfig(config)

	# Printer
	tkDialogs.Printer.cmd   = getStr(_FLAIR_SECTION, "printercmd", "lpr")
	tkDialogs.Printer.paper = getStr(_FLAIR_SECTION, "paper",      "A4")
	tkDialogs.Printer.landscape = getBool(_FLAIR_SECTION,"landscape","True")

	_SKIP_INTRO = getBool(_FLAIR_SECTION, "skipintro", _SKIP_INTRO)
	_SAVE_FRAME = getBool(_FLAIR_SECTION, "rememberframe", _SAVE_FRAME)
	#if _SAVE_FRAME:
	#	_lastFrame = getStr(_FLAIR_SECTION, "lastframe", "project")

	Input.warnMultMat  = getBool(_FLAIR_SECTION, "warnmultmat", Input.warnMultMat)

	try: _autoSaveTime = getInt(_FLAIR_SECTION, "autosave", _autoSaveTime)
	except: pass

	try: tkExtra.ExListbox._KEY_TIME_THRESHOLD = \
		getInt(_FLAIR_SECTION, "keytimethreshold",
			tkExtra.ExListbox._KEY_TIME_THRESHOLD)
	except: pass

	try: tkExtra.Balloon.delay = \
		getInt(_FLAIR_SECTION, "balloondelay", tkExtra.Balloon.delay)
	except: pass

	# Check plot formats
	section = _PLOT_SECTION
	items = config.items(section)
	if not items:
		for k,v in list(_PLOT_FORMAT.items()):
			if k[0] == ".": config.set(section, k, v)
	else:
		_PLOT_FORMAT.clear()
		for k,v in items:
			if k[0] == ".":
				_PLOT_FORMAT[k] = v
	_PLOT_COMMANDS   = getStr(section, "commands","")
	_gnuplotTerminal = getStr(section, "terminal", _gnuplotTerminal)
	Gnuplot.program  = getStr(section, "program",  Gnuplot.program)

	# Geometry Editor
	import InputPage
	import GeometryEditor
	InputPage.configGet()
	GeometryEditor.configGet()

	section = _GEOMETRY_SECTION
	Input.zero     = getFloat(section, "zero",	Input.zero)
	Input.infinite = getFloat(section, "infinite",	Input.infinite)
	Input._useQUA  = getBool(section,   "usequa",	Input._useQUA)
	Input._useBOX  = getBool(section,   "usebox",	Input._useBOX)
#	GeometryViewer.LAPTOP       = getBool(section, "laptop",      GeometryViewer.LAPTOP)
#	GeometryViewer.INVERTWHEEL  = getBool(section, "invertwheel", GeometryViewer.INVERTWHEEL)
#	GeometryViewer.PINZOOM      = getBool(section, "pinzoom",     GeometryViewer.PINZOOM)
#	GeometryViewer.CORES        = getInt( section, "cores",	      GeometryViewer.CORES)
#	GeometryViewer.SNAPDISTANCE = getInt( section, "snap",	      GeometryViewer.SNAPDISTANCE)
#	GeometryViewer.CURSORSIZE   = getInt( section, "cursorsize",  GeometryViewer.CURSORSIZE)

	# Layers section - load default layers if any
	initGlobalLayers()

	# Input (Not user configurable option)
	Input.commentedCards = getBool(_INPUT_SECTION, "parsecomment",Input.commentedCards)

	# ------------- Style Section ------------------
	section = _STYLE_SECTION
	for k,v in config.items(section):
		name, var = k.split(".")
		if "," in v:
			nv = []
			for i in v.split(","):
				if i == "None": i=None
				nv.append(i)
			v = tuple(nv)
		else:
			if v=="None": v=None
		RichText.RichText.TAGS[name][var] = v

#-------------------------------------------------------------------------------
def initGlobalLayers():
	# Layers section - load default layers if any
	Project.initGlobalLayers()
	section = _LAYERS_SECTION
	names = []
	for i in range(100):
		try: name = config.get(section, "layer.%d"%(i))
		except: break
		names.append(name)
		Project._globalLayers[name] = Project.GeometryLayerInfo(name, True)

	# Scan all variables in the layers section
	pat = re.compile(r"^(.*)\.(\d+)$")
	for option in config.options(section):
		m = pat.match(option)
		if m:
			var = m.group(1)
			idx = int(m.group(2))
			try: name = names[idx]
			except KeyError: name= 'Media'
			val = config.get(section, option)
			if val in ("True","False"): val = (val=="True")
			try: val=int(val)
			except ValueError: pass
#			Project._globalLayers[names[idx]][var] = val
			Project._globalLayers[name][var] = val

#-------------------------------------------------------------------------------
# Update the layers from Project -> config
#-------------------------------------------------------------------------------
def updateLayers():
	section = _LAYERS_SECTION
	config.remove_section(section)
	config.add_section(section)

	for i,layer in enumerate(Project._globalLayers.values()):
		config.set(section, "layer.%d"%(i), layer.name)

		for var,val in list(layer.var.items()):
			config.set(section, "%s.%d"%(var,i), str(val))

#-------------------------------------------------------------------------------
# Create the flair directory if needed
#-------------------------------------------------------------------------------
def createIniDir():
	# Create first the ini directory if doesn't exist
	try:
		os.stat(iniDir)
	except OSError:
		os.mkdir(iniDir)

#-------------------------------------------------------------------------------
# Remove from config all the unchanged items
#-------------------------------------------------------------------------------
def cleanIni():
	global config
	newconfig = config	# Remember config
	from tkinter import _default_root
	root = tk.Tk()
	root.withdraw()

	openIni(None,False)	# Open only system ini
	readIni()		# Reading will create all default sections and variables

	import InputPage
	import GeometryEditor
	InputPage.configGet()
	GeometryEditor.configGet()

	# Compare items
	for section in config.sections():
		for item, value in config.items(section):
			try:
				new = newconfig.get(section, item)
				if value==new:
					newconfig.remove_option(section, item)
			except configparser.NoOptionError:
				pass

	config = newconfig
	root.destroy()

#-------------------------------------------------------------------------------
# Write initialization file
#-------------------------------------------------------------------------------
def writeIni():
	createIniDir()
	try:
		fini = open(iniFile, "w")
	except IOError:
		write("Error: Cannot write %s file\n"%(iniFile))
		return
	config.set(__name__, "version", __version__)

	# write unsaved variables
	# 1. global colors	(XXX maybe move to Project)
	Project.MaterialInfo.saveConfig(config)

	# Cleanup
	cleanIni()

	# Write file
	config.write(fini)
	fini.close()

#-------------------------------------------------------------------------------
# Accept global font/color/setting options for all Tkinter widgets
# They have to contain with a star '*' in their name
#	*font = Times 12		Set the global for everything
#	*text*font = Courier -14	Set the Text() widget font (help/viewer..)
#	*text*foreground = White	Set the Text() widget font (help/viewer..)
#-------------------------------------------------------------------------------
def addOptions(section, root):
	if root is None: return
	for name,value in config.items(section):
		if name.find("*")>=0:
			if section == _FONT_SECTION:
				value = fontString(string.split(value,','))
			root.option_add(string.capwords(name,"*"), value)

#-------------------------------------------------------------------------------
def writePrintIni():
	config.set(_FLAIR_SECTION, "printercmd", tkDialogs.Printer.cmd)
	config.set(_FLAIR_SECTION, "paper",      tkDialogs.Printer.paper)
	config.set(_FLAIR_SECTION, "landscape",  str(tkDialogs.Printer.landscape))

#-------------------------------------------------------------------------------
# Load statistics
#-------------------------------------------------------------------------------
def loadStats():
	global stats
	try: f = open(STAT_FILE,"rb")
	except: return
	unpickler = pickle.Unpickler(f)
	try: stats = unpickler.load()
	except: pass
	f.close()

#-------------------------------------------------------------------------------
# Save stats to file
#-------------------------------------------------------------------------------
def saveStats():
	global stats
	try: f = open(STAT_FILE,"wb")
	except: return
	pickler = pickle.Pickler(f)
	pickler.dump(stats)
	f.close()

#-------------------------------------------------------------------------------
def incStats(var):
	global stats
	try: stats[var] += 1
	except: stats[var] = 1

#-------------------------------------------------------------------------------
def resetStats():
	stats.clear()

#-------------------------------------------------------------------------------
# Return an input template
#-------------------------------------------------------------------------------
def inputTemplate(fn=None):
	if fn is None: fn=_DEFAULT_INP
	template = os.path.join(prgDir, _TEMPLATES_DIR, fn)
	try:
		os.stat(template)
	except OSError:
		template = os.path.join(prgDir, _TEMPLATES_DIR, fn)
		try:
			os.stat(template)
		except OSError:
			template = None
	return template

#===============================================================================
# Menu commands
#===============================================================================
# ------------------------------------------------------------------------------
# Initialize underline string
# ------------------------------------------------------------------------------
def initUnderline(cards):
	underline = {}
	for t in cards:
		if isinstance(t, str):
			c = Input.CardInfo.get(t)
		else:
			c = t
		if c.underline>=0:
			underline[c.tag[c.underline]] = True
	return underline

# ------------------------------------------------------------------------------
def findUnderline(ci, underline):
	# find character to underline
	u = ci.underline
	if u < 0:
		u = 0
		for ch in ci.name.upper():
			if ch not in underline:
				underline[ch] = True
				break
			u += 1
		else:
			u = None
	return u

# ------------------------------------------------------------------------------
# create menu items for FLUKA cards
# ------------------------------------------------------------------------------
def addCardCommand(menu, tag, cmd, underline):
	global icons

	if isinstance(tag, str):
		ci = Input.CardInfo.get(tag)
	else:
		ci = tag

	#name = string.capwords(ci.name.lower(),"-")
	u = findUnderline(ci, underline)

	try:
		icon = icons[ci.tag]
	except KeyError:
		icon = icons["tag"]

	menu.add_command(label=ci.name,
			command=cmd,
			underline=u,
			compound=LEFT,
			image=icon)
	return u

# ------------------------------------------------------------------------------
# Create a nice two column body menu
# ------------------------------------------------------------------------------
def _bodiesMenuTearOff(menuname, top):
	# Remove the tearoff bar
	from tkinter import _default_root
	try:
		menu = _default_root.nametowidget(menuname)
	except KeyError:
		# Apparently tear of from the std menu gives a name
		# with # instead of . which generates a keyerror
		return
	if Input._developer:
		menu.delete(20)
	else:
		menu.delete(18)
	# as tkinter doesn't know about the toplevel nametowidget doesn't work
	_default_root.tk.call('wm', 'title', top,  'Bodies')

# ------------------------------------------------------------------------------
def bodiesMenu(menu, callback):
	"""Create a nice bodies menu"""
	def _addCommand(tag, underline=-1):
		try:
			icon = icons[tag]
		except KeyError:
			icon = icons["tag"]
		menu.add_command(label=tag, image=icon, compound=LEFT,
			underline=underline,
			command=lambda t=tag,c=callback: c(t))

	menu.add_command(label="Plane",state=DISABLED)
	_addCommand("YZP", 0)	# X
	_addCommand("XYP", 0)	# Z

	menu.add_command(label="Box",state=DISABLED)
	_addCommand("RPP", 0)	# R
	if Input._useBOX:
		_addCommand("WED", 0)	# W

	menu.add_command(label="Sphere",state=DISABLED)
	_addCommand("SPH", 0)	# S

	menu.add_command(label="Cylinder",state=DISABLED)
	_addCommand("RCC", 2)	# C
	_addCommand("XCC")
	_addCommand("YCC")
	_addCommand("ZCC")
	_addCommand("TRC", 0)	# T

	menu.add_command(label="Other",state=DISABLED)
	_addCommand("QUA", 0)	# Q

	if Input._developer:
		_addCommand("TRX")
		_addCommand("TRZ")

	menu.add_command(label="",state=DISABLED)

	# Second column
	if menu["tearoff"]:
		menu.add_command(label="",state=DISABLED)
		menu["tearoffcommand"] = _bodiesMenuTearOff

	_addCommand("XZP", 1)	# Y
	_addCommand("PLA", 0)	# P

	if Input._useBOX:
		menu.add_command(label="",state=DISABLED)
		_addCommand("BOX", 0)	# B
		_addCommand("RAW")
	else:
		menu.add_command(label="",state=DISABLED)
		menu.add_command(label="",state=DISABLED)

	menu.add_command(label="",state=DISABLED)
	_addCommand("ELL", 1)	# L

	menu.add_command(label="",state=DISABLED)
	_addCommand("REC", 1)	# E
	_addCommand("XEC")
	_addCommand("YEC")
	_addCommand("ZEC")
	menu.add_command(label="",state=DISABLED)
	menu.add_command(label="",state=DISABLED)

	if Input._useBOX:
		_addCommand("ARB")
	else:
		menu.add_command(label="",state=DISABLED)

	column = 15
	if Input._developer:
		_addCommand("TRY")
		column += 2

	if Input._useBOX: column += 1
	if menu["tearoff"]: column += 1
	menu.entryconfigure(column, columnbreak=1)

#-------------------------------------------------------------------------------
def __addABCgroup(menu, callback, abc, cards, **kw):
	# Group Sub Menu
	groupmenu = Menu(menu, **kw)
	menu.add_cascade(label=abc, menu=groupmenu) #, underline=underline)
	underline = initUnderline(cards)
	for c in cards:
		cmd = lambda t=c.tag,c=callback: c(t)
		addCardCommand(groupmenu, c.name, cmd, underline)

#-------------------------------------------------------------------------------
def cardsABCMenu(menu, callback, **kw):
	cards = [c for c in list(Input.CardInfo._db.values()) if "Obsolete" not in c.group]
#	cards.sort(Input.tagCmp)
	cards.sort(key=attrgetter("tag"))
	nsubmenus = len(cards) // 7
	i = 0
	n = 0
	category = ""
	crom     = 0
	lasti    = 0
	while i < len(cards):
		ch = cards[i].tag[0]
		if category.find(ch)<0:
			if n >= nsubmenus:
				if nsubmenus - lasti < n - nsubmenus:
					__addABCgroup(menu,
						callback,
						category[:-1],
						cards[crom : lasti],
						**kw)
					crom = lasti
					i    = lasti
				else:
					__addABCgroup(menu,
						callback,
						category,
						cards[crom : i],
						**kw)
					crom  = i
					lasti = i
				category = ""
				n = 0
			else:
				category += ch
				lasti = i
		n += 1
		i += 1
	__addABCgroup(menu, callback, category[:-1], cards[crom : -1], **kw)

#-------------------------------------------------------------------------------
def cardsMenu(menu, callback, tearoff=1, **kw):
	for group in _groups_order:
		if group == "Developers" and not Input._developer: continue
		# Find underline
		underline = _groups_underline[
				_groups_order.index(group)]

		# Group Sub Menu
		groupmenu = Menu(menu, tearoff=tearoff, **kw)
		menu.add_cascade(label=group,
				menu=groupmenu, underline=underline)

		# scan all cards belonging to a group
		cards = [x for x in list(Input.CardInfo._db.values()) if group in x.group]

		# special treatment for Geometry
		if group==_GEOMETRY_GROUP:

			bodiesmenu = Menu(groupmenu, tearoff=0, **kw)
			transmenu  = Menu(groupmenu, tearoff=tearoff, **kw)

			# remove bodies and transformations from the card list
			trans = [x for x in cards if x.tag[0]=="$"]
			cards = [x for x in cards if x.tag[0]!="$" and \
					len(x.tag)!=3 or x.tag=="END"]

			# Bodies list
			bodiesMenu(bodiesmenu, callback)

			# Transformations list
			transmenu.add_command(label=_EXPANSION,
					command=lambda t=_EXPANSION,
						c=callback: c(t),
					compound=LEFT, image=icons["tags"],
					underline=0)
			transmenu.add_command(label=_TRANSLAT,
					command=lambda t=_TRANSLAT,
						c=callback: c(t),
					compound=LEFT, image=icons["tags"],
					underline=0)
			transmenu.add_command(label=_TRANSFORM,
					command=lambda t=_TRANSFORM,
						c=callback: c(t),
					compound=LEFT, image=icons["tags"],
					underline=1)

			underline = initUnderline(trans)
			for c in trans:
				cmd = lambda t=c.tag,c=callback: c(t)
				addCardCommand(transmenu, c.name, cmd, underline)

		# Sort cards
#		cards.sort(Input.tagCmp)
		cards.sort(key=attrgetter("tag"))
		underline = initUnderline(cards)
		for c in cards:
			cmd = lambda t=c.tag,c=callback: c(t)
			if group==_GEOMETRY_GROUP and c.tag=="END":
				groupmenu.add_cascade(label="Bodies",
						menu=bodiesmenu,
						image=icons["TRC"],
						compound=LEFT,
						underline=0)
				groupmenu.add_cascade(label="Transformations",
						menu=transmenu,
						image=icons["ROT-DEFI"],
						compound=LEFT,
						underline=0)
			addCardCommand(groupmenu, c.name, cmd, underline)

		if group==_PREPROCESSOR_GROUP:
			groupmenu.add_command(label=_IFENDIF,
					command=lambda t=_IFENDIF,
						c=callback: c(t),
					compound=LEFT, image=icons[_IFENDIF],
					underline=0)

#===============================================================================
def addException():
	global _output
	#self.widget._report_exception()
	try:
		typ, val, tb = sys.exc_info()
		traceback.print_exception(typ, val, tb)
		if errors: errors.append("")
		exception = traceback.format_exception(typ, val, tb)
		errors.extend(exception)
		_output += "\n".join(exception)
		if len(errors) > 100:
			# If too many errors are found send the error report
			ReportDialog(self.widget)
	except:
		say(sys.exc_info())

#===============================================================================
class CallWrapper:
	"""Replaces the Tkinter.CallWrapper with extra functionality"""
	def __init__(self, func, subst, widget):
		"""Store FUNC, SUBST and WIDGET as members."""
		self.func   = func
		self.subst  = subst
		self.widget = widget

	# ----------------------------------------------------------------------
	def __call__(self, *args):
		"""Apply first function SUBST to arguments, than FUNC."""
		try:
			if self.subst:
				args = self.subst(*args)
			return self.func(*args)
		# One possible fix is to make an external file for the wrapper
		# and import depending the version
		#except SystemExit, msg:	# python2.4 syntax
		#except SystemExit as msg:	# python3 syntax
		#	raise SystemExit(msg)
		except SystemExit:		# both
			raise SystemExit(sys.exc_info()[1])
		except KeyboardInterrupt:
			pass
		except:
			addException()

#===============================================================================
# Error message reporting dialog
#===============================================================================
class ReportDialog(Toplevel):
#	_shown = False	# avoid re-entry when multiple errors are displayed
	_shown = True	# avoid re-entry when multiple errors are displayed

	def __init__(self, master):
		if ReportDialog._shown: return
		ReportDialog._shown = True

		Toplevel.__init__(self, master)
		self.title("%s Error Reporting"%(__name__))
		#self.transient(master)

		# Label Frame
		frame = LabelFrame(self, text="Report")
		frame.pack(side=TOP, expand=YES, fill=BOTH)

		l = Label(frame, text="The following report is about to be sent "\
				"to the mantainer  of this version of  %s"%(__name__), justify=LEFT, anchor=W)
		l.pack(side=TOP)

		self.text = Text(frame, background="White")
		self.text.pack(side=LEFT, expand=YES, fill=BOTH)

		sb = Scrollbar(frame, orient=VERTICAL, command=self.text.yview)
		sb.pack(side=RIGHT, fill=Y)
		self.text.config(yscrollcommand=sb.set)

		# email frame
		frame = Frame(self)
		frame.pack(side=TOP, fill=X)

		l = Label(frame, text="Your email")
		l.pack(side=LEFT)

		self.email = Entry(frame, background="White")
		self.email.pack(side=LEFT, expand=YES, fill=X)

		# Automatic error reporting
		self.err = BooleanVar()
		self.err.set(_errorReport)
		b = Checkbutton(frame, text="Automatic error reporting",
			variable=self.err, anchor=E, justify=RIGHT)
		b.pack(side=RIGHT)

		# Buttons
		frame = Frame(self)
		frame.pack(side=BOTTOM, fill=X)

		b = Button(frame, text="Close",
				image=icons["x"],
				compound=LEFT,
				command=self.cancel)
		b.pack(side=RIGHT)
		b = Button(frame, text="Send report",
				image=icons["debug"],
				compound=LEFT,
				command=self.send)
		b.pack(side=RIGHT)

		import GeometryEditor

		# FIXME check version PIL, numpy, pydicom
		try:
			import numpy
			numpy_version = numpy.__version__
		except ImportError:
			numpy_version = "?"

		try:
			import pydicom
			dicom_version = pydicom.__version__
		except ImportError:
			try:
				import dicom
				dicom_version = dicom.__version__
			except ImportError:
				dicom_version = "?"

		# Fill report
		txt = [ "Version     : %s"%(__version__),
			"geoviewer   : %s %s"%(GeometryEditor.version(),
					__version__!=GeometryEditor.version() and "WARNING"),
			"Revision    : %s"%(__revision__),
			"Last Change : %s"%(__lastchange__),
			"dicom       : %s"%(dicom_version),
			"numpy       : %s"%(numpy_version),
			"Platform    : %s"%(sys.platform),
			"Python      : %s"%(sys.version),
			"TkVersion   : %s"%(TkVersion),
			"TclVersion  : %s"%(TclVersion),
			"\nTraceback:" ]
		for e in errors:
			if e!="" and e[-1] == "\n":
				txt.append(e[:-1])
			else:
				txt.append(e)

		self.text.insert('0.0', "\n".join(txt))

		# Guess email
		user = os.getenv("USER")
		host = os.getenv("HOSTNAME")
		if user and host:
			email = "%s@%s"%(user,host)
		else:
			email = ""
		email = getStr(_FLAIR_SECTION, "email", email)
		self.email.insert(0,email)

		self.protocol("WM_DELETE_WINDOW", self.close)
		self.bind('<Escape>', self.close)

		# Wait action
		self.wait_visibility()
		self.grab_set()
		self.focus_set()
		self.wait_window()

	# ----------------------------------------------------------------------
	def close(self, event=None):
		ReportDialog._shown = False
		self.destroy()

	# ----------------------------------------------------------------------
	def send(self):
		import http.client, urllib.request, urllib.parse, urllib.error
		global errors
		email = self.email.get()
		desc  = self.text.get('1.0', END).strip()
		config.set(_FLAIR_SECTION, "email", self.email.get())

		# Send information
		self.config(cursor="watch")
		self.text.config(cursor="watch")
		self.update_idletasks()
		params = urllib.parse.urlencode({"email":email, "desc":desc})
		headers = {"Content-type": "application/x-www-form-urlencoded",
			"Accept": "text/plain"}
		conn = http.client.HTTPConnection("www.fluka.org:80")
		try:
			conn.request("POST", "/flair/send_email.php", params, headers)
			response = conn.getresponse()
		except:
			messagebox.showwarning("Error sending report",
				"There was a problem connecting to the FLUKA web site",
				parent=self)
		else:
			if response.status == 200:
				messagebox.showinfo("Report successfully send",
					"Report was successfully uploaded to %s web site"%\
					(__name__),
					parent=self)
				del errors[:]
			else:
				messagebox.showwarning("Error sending report",
					"There was an error sending the report\nCode=%d %s"%\
					(response.status, response.reason),
					parent=self)
		conn.close()
		self.config(cursor="")
		self.cancel()

	# ----------------------------------------------------------------------
	def cancel(self):
		global _errorReport, errors
		_errorReport = self.err.get()
		config.set(_FLAIR_SECTION, "errorreport", str(bool(self.err.get())))
		del errors[:]
		self.close()

	# ----------------------------------------------------------------------
	@staticmethod
	def sendErrorReport():
		ReportDialog(None)

#===============================================================================
def preferences(self):
	global config
	import ConfigDialog
	cd = ConfigDialog.ConfigDialog(self, config)
	if cd.show():
		from tkinter import _default_root
		# Load fonts...
		readIni()
		addOptions(_FONT_SECTION, self.winfo_toplevel())
		addOptions(_COLOR_SECTION, self.winfo_toplevel())
		import InputPage
		import GeometryEditor
		InputPage.configGet()
		GeometryEditor.configGet()
		return True
	return False

#===============================================================================
# About dialog box
#===============================================================================
def aboutDialog(parent, timer=None):
	toplevel = Toplevel(parent)
	toplevel.transient(parent)
#	toplevel.title("About %s" % (string.capitalize(__name__)))
	toplevel.title("About %s" % (str.capitalize(__name__)))

	bg = "#707070"
	fg = "#ffffff"

	font1 = 'Helvetica -32 bold'
	font2 = 'Helvetica -12'
	font3 = 'Helvetica -10'

	frame = Frame(toplevel, borderwidth=2,
			relief=SUNKEN, background=bg)
	frame.pack(side=TOP, expand=TRUE, fill=BOTH, padx=5, pady=5)

	l = Label(frame, image=icons["flair"],
			foreground=fg, background=bg,
			relief=RAISED,
			padx=0, pady=0)
	l.grid(row=0, column=0, columnspan=2, padx=5, pady=5)

	l = Label(frame, text=__name__,
			foreground=fg, background=bg,
			font=font1)
	l.grid(row=1, column=0, columnspan=2, sticky=W, padx=10, pady=5)

	l = Label(frame, text="FLUKA Advanced Interface",
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=2, column=0, columnspan=2, sticky=W, padx=10, pady=2)

	l = Label(frame, text="/fle\u0259(r)/\tn [U,C] natural or " \
			"instinctive ability - to do something\n" \
			"\twell, to select or recognize what is best, " \
			"more useful, etc.\n" \
			"\t\t\t[Oxford Advanced Dictionary]",
			font = font3,
			foreground=fg, background=bg, justify=LEFT)
	l.grid(row=3, column=0, columnspan=2, sticky=W, padx=10, pady=1)

	f = Frame(frame, borderwidth=1, relief=SUNKEN,
		height=2, background=bg)
	f.grid(row=4, column=0, columnspan=2, sticky=EW, padx=5, pady=5)

	# -----
	l = Label(frame, text='www:',
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=5, column=0, sticky=E, padx=10, pady=2)

	l = Label(frame, text=__www__,
			foreground=fg, background=bg, justify=LEFT,
			activeforeground="Blue",
			font=font2, cursor="hand1")
	l.grid(row=5, column=1, sticky=W, padx=2, pady=2)
	l.bind('<Button-1>', aboutWebsite)

	# -----
	l = Label(frame, text='email:',
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=6, column=0, sticky=E, padx=10, pady=2)

	l = Label(frame, text=__email__,
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=6, column=1, sticky=W, padx=2, pady=2)

	# -----
	l = Label(frame, text='author:',
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=7, column=0, sticky=NE, padx=10, pady=2)

	l = Label(frame, text=__author__.replace(",","\n"),
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=7, column=1, sticky=NW, padx=2, pady=2)

	# -----
	l = Label(frame, text='contributors:',
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=8, column=0, sticky=NE, padx=10, pady=2)

	l = Label(frame, text=__credits__.replace(",","\n"),
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=8, column=1, sticky=NW, padx=2, pady=2)

	# -----
	l = Label(frame, text='version:',
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=9, column=0, sticky=E, padx=10, pady=2)

	l = Label(frame, text="%s [R%s]" % (__version__, __revision__),
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=9, column=1, sticky=NW, padx=2, pady=2)

	# -----
	l = Label(frame, text='last change:',
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=10, column=0, sticky=E, padx=10, pady=2)

	l = Label(frame, text=__lastchange__,
			foreground=fg, background=bg, justify=LEFT,
			font=font2)
	l.grid(row=10, column=1, sticky=NW, padx=2, pady=2)

	closeFunc = lambda e=None,t=toplevel: t.destroy()
	b = Button(toplevel, text="Close", command=closeFunc)
	b.pack(pady=5)
	frame.grid_columnconfigure(1, weight=1)

	toplevel.bind('<Escape>',   closeFunc)
	toplevel.bind('<Return>',   closeFunc)
	toplevel.bind('<KP_Enter>', closeFunc)

	toplevel.deiconify()
	toplevel.wait_visibility()
	toplevel.resizable(False, False)
	try: toplevel.grab_set()
	except: pass
	b.focus_set()
	toplevel.lift()
	if timer:
		toplevel.after(timer, closeFunc)
	toplevel.wait_window()

#-------------------------------------------------------------------------------
def aboutWebsite(event=None):
	if openurl(__www__):
		messagebox.showwarning("No valid browser",
			"Cannot find any valid browser to open web page")

#===============================================================================
# FlukaOutputFiles list
#===============================================================================
class FlukaFilesList(tkExtra.ColorMultiListbox):
	def __init__(self, master, flair):
		tkExtra.ColorMultiListbox.__init__(self, master,
				(('File', 30, None),
				 ('Type', 10, None),
				 ('Size',  5, None),
				 ('Date', 13, None)))
		self.bindList('<Double-1>', self.view)
		self.flair = flair

	# ----------------------------------------------------------------------
	# Fill with files generated by fluka
	# ----------------------------------------------------------------------
	def fill(self, inpname, pid=None, cycle=None, project=None):
		self.project = project
		# Scan directory
		if isinstance(cycle, int):
			pat = re.compile(r"^%s%03d[._].+$" % (inpname,cycle))
		else:
			pat = re.compile(r"^%s\d\d\d[._].+$" % (inpname))

		for filename in os.listdir("."):
			if pat.match(filename):
				try:
					s = os.lstat(filename)
					if S_ISLNK(s[ST_MODE]):
						continue
					ext,color = bFileDialog.fileTypeColor(filename)
					self.insert(END,
						(filename, ext, s[ST_SIZE],
						 time.strftime(bFileDialog._TIME_FORMAT,
							time.localtime(s[ST_MTIME]))))
					if color is not None:
						self.itemconfig(END,foreground=color)
				except OSError:
					pass

		# Scan temporary directory
		if pid:
			try:
				for fn in os.listdir("fluka_%d"%(pid)):
					filename = "fluka_%d/%s"%(pid, fn)
					try:
						s = os.lstat(filename)
						if S_ISLNK(s[ST_MODE]):
							continue
						ext,color = bFileDialog.fileTypeColor(filename)
						self.insert(END,
							(filename, ext, s[ST_SIZE],
							 time.strftime(bFileDialog._TIME_FORMAT,
								time.localtime(s[ST_MTIME]))))
						if color is not None:
							self.itemconfig(END,foreground=color)
					except OSError:
						pass
			except OSError:
				pass

		if self.size() == 0:
			messagebox.showerror("Cannot open output",
				"Cannot find any output file.",
				parent=self)
			return

	# ----------------------------------------------------------------------
	# View document
	# ----------------------------------------------------------------------
	def view(self, event=None):
		files = []
		for sel in self.curselection():
			files.append(self.get(sel)[0])
		self.flair.viewer(files)

	# ----------------------------------------------------------------------
	# Edit document
	# ----------------------------------------------------------------------
	def edit(self, event=None):
		files = []
		for sel in self.curselection():
			files.append(self.get(sel)[0])
		editor(files)

#===============================================================================
# OutputFilesDialog
#===============================================================================
class OutputFilesDialog:
	def __init__(self, flair, inpname, pid=None, cycle=None, project=None):
		self.toplevel = Toplevel(flair)
		self.toplevel.transient(flair)
		self.flair = flair

		listbox = FlukaFilesList(self.toplevel, flair)
		listbox.pack(side=TOP, expand=YES, fill=BOTH)
		listbox.fill(inpname, pid, cycle, project)

		frame = Frame(self.toplevel)
		frame.pack(side=BOTTOM, fill=X)

		btn = Button(frame, text="View", underline=0,
				command=listbox.view)
		btn.pack(side=LEFT)

		btn = Button(frame, text="Edit", underline=0,
				command=listbox.edit)
		btn.pack(side=LEFT)

		btn = Button(frame, text="Close", underline=0,
			command=self.close)
		btn.pack(side=LEFT)

		self.toplevel.bind("<Escape>", self.close)
		self.toplevel.protocol("WM_DELETE_WINDOW", self.close)
		self.toplevel.title("Output File Viewer")

		# Wait action
		try:
			#self.toplevel.grab_set()
			self.toplevel.focus_set()
			self.toplevel.wait_window()
		except TclError:
			pass

	# ----------------------------------------------------------------------
	def close(self, event=None):
		self.toplevel.destroy()

#===============================================================================
# Undo listbox
#===============================================================================
class UndoListbox(Toplevel):
	posted = None

	def __init__(self, master, flair, x, y):
		Toplevel.__init__(self, master)
		self.overrideredirect(1)
		self.transient(master)

		self.flair = flair

		# Listbox
		frame = Frame(self)
		frame.pack(expand=YES,fill=BOTH)

		self.listbox = tkExtra.ExListbox(frame, borderwidth=0,
			selectmode=SINGLE,
			selectborderwidth=0,
			background="White",
			takefocus=True,
			exportselection=FALSE)
		self.listbox.pack(side=LEFT, expand=YES, fill=BOTH)

		sb = Scrollbar(frame, orient=VERTICAL, command=self.listbox.yview)
		sb.pack(side=RIGHT, fill=Y)
		self.listbox.config(yscrollcommand=sb.set)

		frame = Frame(self)
		frame.pack(side=BOTTOM)

		self.listbox.bind("<Motion>", self.motion)
		self.listbox.bind("<ButtonRelease-1>", self.release)
		self.bind("<Escape>", self.cancel)
		self.listbox.bind('<FocusOut>',	self.cancel)
		self.bind("<Return>",   self.ok)
		self.bind("<KP_Enter>", self.ok)

		for lbl in self.flair.undoredo.undoTextList():
			self.listbox.insert(END, lbl)

		self.listbox.activate(0)
		self.listbox.selection_set(0)
		self.listbox.focus_set()
		self.focus_set()

		self.geometry("+%d+%d"%(x,y))
		UndoListbox.posted = self

	#-----------------------------------------------------------------------
	def motion(self, event):
		y = self.listbox.nearest(event.y)
		self.listbox.selection_clear(0,END)
		self.listbox.selection_set(0,y)
		self.listbox.activate(y)

	#-----------------------------------------------------------------------
	def release(self, event):
		y = self.listbox.nearest(event.y)
		self.listbox.selection_clear(0,END)
		self.listbox.selection_set(y)
		self.listbox.activate(y)
		self.ok()

	#-----------------------------------------------------------------------
	def ok(self, event=None):
		sel = self.listbox.curselection()
		if not sel: return
		act = int(sel[0])
		action = self.listbox.get(act)
		for i in range(act+1):
			self.flair.write(self.flair.undoredo.undoText()+"\n")
			self.flair.undoredo.undo()
			self.flair.refresh()
		self.flair.project.input.clearCache("bodies")
		self.destroy()
		UndoListbox.posted = None
		self.master.refreshUndoButton()

	#-----------------------------------------------------------------------
	def cancel(self, event=None):
		self.destroy()
		UndoListbox.posted = None

#-------------------------------------------------------------------------------
# Edit files
#-------------------------------------------------------------------------------
def editor(files):
	if len(files) == 0: return
	if isinstance(files,str): files = [ files ]
	try:
		os.spawnv(os.P_NOWAIT, Project.editor,
				[ Project.editor ] + files )
	except:
		messagebox.showerror("Error executing command",
			"Error executing command:\n%s %s" \
			% (Project.editor, " ".join(files)))

#-------------------------------------------------------------------------------
def openurl(url):
	try:
		if webbrowser.get().basename in ("links", "elinks", "lynx"):
			write("Error: No valid web browser found\n")
			return True
	except:
		write("Error: No valid web browser found\n")
		return True
	webbrowser.open_new(url)
	return False

#-------------------------------------------------------------------------------
def printer(master):
	printer = tkDialogs.Printer(master)
	if printer.show():
		writePrintIni()
		return printer
	else:
		return None

#-------------------------------------------------------------------------------
# Stupid fix for python 2.3 the askyesno
# returns sometimes "yes" and sometimes "True"
#-------------------------------------------------------------------------------
def askyesno(title, message, **options):
	ans = messagebox.askquestion(title, message, **options)
	return str(ans)==messagebox.YES or ans==True

#-------------------------------------------------------------------------------
# Show info message and close dialog after few seconds
#-------------------------------------------------------------------------------
class ShowInfo(Toplevel):
	def __init__(self, master, title, message, timeout=5000):
		Toplevel.__init__(self,master)
		self.title(title)
		self.transient(master)

#		Label(self, image="info").pack(side=LEFT)
		Message(self, text=message, font=",-10,bold", aspect=300).pack(side=RIGHT)

		self.after(timeout, self.close)
		self.bind("<Button-1>",  self.close)
		self.bind("<Button-2>",  self.close)
		self.bind("<Button-3>",  self.close)
		self.bind("<Key>",	 self.close)
		self.protocol("WM_DELETE_WINDOW", self.close)
		self.wait_window()

	#-----------------------------------------------------------------------
	def close(self, event=None):
		self.destroy()

#-------------------------------------------------------------------------------
def __init_group_underline():
	# Groups underline
	gunderline = ""
	for group in _groups_order:
		# Find underline
		u = 0
		for ch in group.upper():
			if gunderline.find(ch)==-1:
				gunderline += ch
				break
			u += 1
		else:
			u = None
		_groups_underline.append(u)

#-------------------------------------------------------------------------------
def init(ini=None):
	openIni(ini)
	readIni()
	loadStats()
	loadIcons()
	__init_group_underline()

#-------------------------------------------------------------------------------
def fini():
	delIcons()
