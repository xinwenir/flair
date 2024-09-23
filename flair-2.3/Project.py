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

__author__  = "Vasilis Vlachoudis"
__email__   = "Vasilis.Vlachoudis@cern.ch"
__version__ = "2"

import os
import re
import sys
import stat
import time
import shlex
import string
import random
import shutil
import subprocess
import tempfile
from log import say
from fnmatch import fnmatch

import bmath
import Input
import Utils

__first = True
section  = "Project"

#===============================================================================
# Global constants
DEFAULT_INPUT = "<default>"
DEFAULT_EXE   = "flukahp"
DEFAULT_LIB   = "libflukahp.a"

STATUS_NOT_RUNNING  = 0
STATUS_WAIT2ATTACH  = 1
STATUS_RUNNING      = 2
STATUS_FINISHED     = 3
STATUS_FINISHED_ERR = 4
STATUS_TIMEOUT      = 5

MAX_TIME            = 1.0E100

_statusColor = {
	STATUS_NOT_RUNNING  : "Black",
	STATUS_WAIT2ATTACH  : "Orange",
	STATUS_RUNNING      : "DarkGreen",
	STATUS_FINISHED     : "DarkBlue",
	STATUS_FINISHED_ERR : "Red",
	STATUS_TIMEOUT      : "Purple"
}

BOLD   = "\033[1m"
NORMAL = "\033[m"

FLUFOR_OPTIONS = ["", "g77", "gfortran" ]

#===============================================================================
# Global variables
timeThreshold   = 120			# more that XXs the run is considered
					# ... to be finished
refreshInterval = 10 * 1000		# ms
tmpPrefix       = "flair_"		# Commands and tmp files
cleanup         = True			# Cleanup temporary files
editor          = "/usr/bin/emacs"	# Default editor
kill            = ""			# kill command
terminal        = "xterm"		# terminal command
debugger        = "gdb"			# debugger
#alphaSortOrder = True			# For tree view menus
keepBackup      = True			# Keep backups
flukaVar        = "FLUPRO"		# Default Fluka variable
flufor          = ""			# Default compiler
flukaDir        = ""			# Override fluka directory
flukaExe        = DEFAULT_EXE		# Default fluka executable
ccprogram       = "/usr/bin/cc"		# C compiler
cppprogram      = "/usr/bin/c++"	# C++ compiler

_checkedVar     = None			# Is directory already checked?

# fluka programs, the name is the default program
_fluka_programs  = {    "fff"    : "",
			"rfluka" : "",
			"detsuw" : "",
			"usbsuw" : "",
			"ustsuw" : "",
			"usxsuw" : "",
			"usbmax" : "",
			"usrsuw" : "",
			"usbrea" : "",
			"gplevbin" : "",
			"d": "detsuw",
			"r": "usrsuw",
			"b": "usbsuw",
			"x": "usxsuw",
			"t": "ustsuw",
			"y": "usysuw",
			"c": "ustsuw",
		}

LINK_PRG         = [	"lfluka", "lflukac",
			"ldpmqmd",    "ldpm2qmd", "ldpm3qmd" ]

# USRxxx processing
# default types
_usr_type = {	"d": "detect",
		"r": "resnuclei",
		"b": "usrbin",
		"x": "usrbdx",
		"t": "usrtrack",
		"y": "usryield",
		"c": "usrcoll",
		"detect"   : "d",
		"resnuclei": "r",
		"usrbin"   : "b",
		"usrbdx"   : "x",
		"usrtrack" : "t",
		"usryield" : "y",
		"usrcoll"  : "c",
	}

# USRxxx extension
_usr_ext = {	"d": "dtc",
		"r": "rnc",
		"b": "bnn",
		"x": "bnx",
		"t": "trk",
		"y": "yie",
		"c": "cll"
	}

# Old plot variable translate to new one
_plot_translate = {
	"titleopts"   : "titleoptions",
	"xlabelopts"  : "xoptions",
	"ylabelopts"  : "yoptions",
	"cblabelopts" : "cboptions",
}

# default names
_DEFAULT_USRNAME = r'\I_\U.\e'
_usr_name = {}

# default clone input name
_DEFAULT_RUNNAME = r'\I_\a'

# fluka cycle directory
_FLUKA_CYCLE = "fluka_"

# submit commands
_submit_command = { }

# default rules
_usr_rule = [	('+',r"\I\d\d\d_fort\.\U"),
		('+',r"\I\d\d\d_ftn\.\U") ]

# patterns
_LAST_NUM = re.compile(r"^(.*?) *(\d+)$")
_PAT_CORE = re.compile(r"^core\.\d+$")
_PAT_DIR  = re.compile(r"^fluka_(\d+)$")
_PAT_DIR2 = re.compile(r"^(.*/?)fluka_(\d+)$")

# list of global layers GeometryLayerInfo
_systemLayers = {}
_globalLayers = {}

__gnucolors = ("black","red", "green","blue", "magenta","cyan", "yellow")

#-------------------------------------------------------------------------------
# Initialize global layers to the 4 basic ones
#-------------------------------------------------------------------------------
def initSystemLayers():
	global _systemLayers
	_systemLayers.clear()

	layer = GeometryLayerInfo("Borders",True)
	_systemLayers[layer.name] = layer
	layer["_system"]    = True
	layer["_order"]     = 0
	layer["Show"]       = True
	layer["Show.color"] = "None"

	layer = GeometryLayerInfo("Media",True)
	_systemLayers[layer.name] = layer
	layer["_system"]    = True
	layer["_order"]     = 1
	layer["Show"]       = True
	layer["Show.color"] = "Material"
	layer["Options.shortcut"] = 'F3'

	layer = GeometryLayerInfo("3D",True)
	_systemLayers[layer.name] = layer
	layer["_system"]      = True
	layer["_order"]       = 2
	layer["Show"]         = False
	layer["Show.lattice"] = True
	layer["Show.voxel"]   = True
	layer["3D"]           = True
	layer["deflights"]    = True
	layer["Options.shortcut"] = 'F4'

	layer = GeometryLayerInfo("Lattice",True)
	_systemLayers[layer.name] = layer
	layer["_system"]      = True
	layer["_order"]       = 3
	layer["Show"]         = True
	layer["Show.color"]   = "Material"
	layer["Show.lattice"] = True
	layer["Show.voxel"]   = True

	initGlobalLayers()

#-------------------------------------------------------------------------------
def initGlobalLayers():
	global _globalLayers
	_globalLayers.clear()

	for layer in list(_systemLayers.values()):
		_globalLayers[layer.name] = layer.clone()

#-------------------------------------------------------------------------------
# Load configuration
#-------------------------------------------------------------------------------
def loadConfig(config):
	global tmpPrefix, cleanup, timeThreshold, refreshInterval, keepBackup
	global editor, terminal, debugger, kill, flukaVar, flufor, flukaExe
	from tkFlair import getStr, getBool, getInt

	tmpPrefi        = getStr(section, "tmpprefix",       tmpPrefix)
	cleanup         = getBool(section,"cleanup",         cleanup)
	timeThreshold   = getInt(section, "timethreshold",   timeThreshold)
	refreshInterval = getInt(section, "refreshinterval", refreshInterval/1000)*1000
	keepBackup      = getBool(section,"keepbackup",      keepBackup)
	editor          = getStr(section, "editor",          editor)
	terminal        = getStr(section, "terminal",        terminal)
	debugger        = getStr(section, "debugger",        debugger)
	kill            = getStr(section, "kill",            kill)
	flukaVar        = getStr(section, "flukavar",        flukaVar)
	flufor          = getStr(section, "flufor",          flufor)

	setCmd("rfluka",  getStr(section, "rfluka",  ""))
	setCmd("rfluka",  getStr(section, "rfluka",  ""))
	setCmd("fff",     getStr(section, "fff",     ""))
	setCmd("usbsuw",  getStr(section, "usbsuw",  ""))
	setCmd("ustsuw",  getStr(section, "ustsuw",  ""))
	setCmd("usxsuw",  getStr(section, "usxsuw",  ""))
	setCmd("usrsuw",  getStr(section, "usrsuw",  ""))
	setCmd("usysuw",  getStr(section, "usysuw",  ""))
	setCmd("detsuw",  getStr(section, "detsuw",  ""))
	setCmd("usbmax",  getStr(section, "usbmax",  ""))
	setCmd("usbrea",  getStr(section, "usbrea",  ""))
	setCmd("gplevbin",getStr(section, "gplevbin",""))

	setFlukaDir(getStr(section, "flukadir", ""))
	flukaExe =  getStr(section, "flukaexe", "flukahp")

#-------------------------------------------------------------------------------
# Add rule from configuration file
#-------------------------------------------------------------------------------
def addRule(rule, pos=None):
	pat  = "+"
	name = rule
	if len(rule)>1:
		pat  = rule[0]
		if pat=="+" or pat=="-":
			name = rule[1:]
		else:
			pat = "+"
	rule = (pat, name)

	if pos is None:
		_usr_rule.append(rule)
	else:
		if pos >= len(_usr_rule):
			_usr_rule.extend([("","")]*(pos-len(_usr_rule)+1))
		_usr_rule[pos] = rule

#-------------------------------------------------------------------------------
# Helper function, return hash of a dictionary
#-------------------------------------------------------------------------------
def dicthash(dictionary):
	h = 0
	for k,v in list(dictionary.items()):
		if v != "":
			try:
				try:
					v = str(v)
				except UnicodeEncodeError:
					v = str(v).encode("utf-8")
				h += hash(k)*167 + hash(v)*171
			except TypeError:
				pass
	return h

#===============================================================================
# Read line and split it to variable and value
#-------------------------------------------------------------------------------
def _readVar(f):
	while True:
		line = Input.utfReadline(f)
		if len(line)==0: return (None, None)
		line = line.strip()
		if line=="": continue
		if line[0] == "#" or len(line) == 0: continue
		try:
			split = line.index(":")
			tag   = str(line[0:split].strip())
			val   = line[split+1:].strip()
			return (tag, val)
		except:
			return (line, "")

#-------------------------------------------------------------------------------
# read a block of lines until ^L character
#-------------------------------------------------------------------------------
def _readBlock(f):
	block = ""
	while True:
		line = Input.utfReadline(f)
		if len(line)==0: break
		ep = line.find('\f')
		if ep >= 0:
			block += line[0:ep]
			break
		else:
			block += line
	return block

#-------------------------------------------------------------------------------
def setFlukaDir(path=""):
	global flukaDir
	flukaDir = path
	if flukaDir == "":
		d = os.getenv(flukaVar)
		if d is not None:
			flukaDir = d
	else:
		os.putenv(flukaVar, flukaDir)
	if flufor:
		os.putenv("FLUFOR", flufor)

#-------------------------------------------------------------------------------
def checkDir():
	global flukaDir, _checkedVar

	try:
		os.stat(flukaDir)
		return False
	except OSError:
		say("ERROR: FLUKA directory is not found or not set")

	# Avoid repeating messages
	if flukaVar==_checkedVar: return True

	_checkedVar = flukaVar
	say("WARNING: No %s variable declared."%(flukaVar))
	say("         Please set %s variable prior calling flair"%(flukaVar))
	say("         or set the FLUKA directory in the Preferences Dialog")
	return True

#-------------------------------------------------------------------------------
def setCmd(cmd, val):
	val = val.strip()
	path = os.path.dirname(val)
	if len(path)==0 or path==os.path.join(flukaDir, "flutil"):
		_fluka_programs[cmd] = os.path.basename(val)
	else:
		_fluka_programs[cmd] = val

#-------------------------------------------------------------------------------
def command(prg):
	"""Return a fluka command, default if cmd=empty"""
	path = _fluka_programs.get(prg, "")
	if checkDir(): return path
	if not path:
		cmd = os.path.join(flukaDir, "flutil", prg)
		# Try fluka directory
		try:
			os.stat(cmd)
			return cmd
		except:
			pass

		# Try flair directory
		cmd = os.path.join(os.path.dirname(__file__), prg)
		try:
			os.stat(cmd)
			return cmd
		except:
			pass
		return prg
	else:
		if path[0] == os.path.sep:
			return path
		else:
			return os.path.join(flukaDir, "flutil", path)

#-------------------------------------------------------------------------------
# Make path relative or smallest representation
#-------------------------------------------------------------------------------
def relativePath(path, cwd):
	"""Convert path to relative to the cwd current working directory"""
	if len(path)==0: return path
	dirpath = os.path.dirname(os.path.normpath(path))

	# if already relative return
	if len(dirpath) == 0 or \
		(dirpath[0] != os.sep and dirpath[0] != "."): return path

	if path==cwd: return ""

#	try:
#		# python 2.6
#		# fails to return absolute e.g.  /home/bnv  /clueet/bnv
#		return os.path.relpath(path, cwd)
#	except:
#		pass

	# Find common path
	d = cwd.split(os.sep)
	p = path.split(os.sep)
	cp = ""
	for i in range(min(len(d),len(p))):
		if d[i]!="":
			if d[i] == p[i]:
				cp += os.sep + d[i]
			else:
				break
	if len(cp)>0: cp += os.sep
	dirsuffix = cwd[len(cp):]
	if len(dirsuffix)>0:
		relpath = (os.pardir+os.sep) * len(dirsuffix.split(os.sep))
	else:
		relpath = ""
	relpath += path[len(cp):]
	if len(relpath)<len(path):
		return relpath
	return path

#-------------------------------------------------------------------------------
# Abbreviate a large path to fit in small entry
#-------------------------------------------------------------------------------
def abbrevPath(path, maxlen=40):
	if len(path)<maxlen: return path

	prefix = ""
	if path.find(os.sep)>=0:
		d = os.path.dirname(path).split(os.sep)
		# keep only the first and last
		del d[1:-2]
		d[1] = "..."
		prefix = " ".join(d, os.sep) + os.sep
		if len(prefix) > maxlen:
			del d[-1]
			prefix = " ".join(d, os.sep) + os.sep

	filename = os.path.basename(path)
	maxlen -= len(prefix)
	if len(filename) >= maxlen:
		# remove middle characters
		rm = (maxlen - 3)//2
		if rm<1: rm = 1
		filename = "%s...%s"%(filename[:rm], filename[-rm:])

	return "%s%s"%(prefix,filename)

#===============================================================================
# Rules for commands
#===============================================================================
class Command:
	section  = "Data"
	commands = {}	# available commands

	# ----------------------------------------------------------------------
	def __init__(self, idx, name, cmd):
		self.name = name
		if cmd:
			self.cmd = cmd
		else:
			self.cmd = name
		self.idx  = idx
		self.arg  = []
		self.inp  = []

	# ----------------------------------------------------------------------
	# Append an argument string separated with ":"
	#
	#  tag, type, default, name
	#
	#  tag:  is the argument tag as "-x"
	#  type: int, float, output, file, filelist, fixed, or list of values
	#        look the DataTab.inputInPlaceEdit for data types
	#  default: default value, can be any expression on the input
	#        e.g. =C(VOXELS,0,0)+".vxl"
	#  name: to appear in listbox
	# ----------------------------------------------------------------------
	def addArgument(self, arg):
		self.arg.append(arg.split(":"))

	# ----------------------------------------------------------------------
	# Append an input string separated with ":"
	#
	#  type, default, name
	#
	#  type: like in the arguments
	#        + conditional
	#	repeat : value : label	(repeat from label until value is given)
	# ----------------------------------------------------------------------
	def addInput(self, inp):
		self.inp.append(inp.split(":"))

	# ----------------------------------------------------------------------
	# Return default name
	# ----------------------------------------------------------------------
	def defaultName(self):
		for tag,ty,default,name in self.arg:
			if ty=="output": return default
		for ty,default,name in self.inp:
			if ty=="output": return default
		return _DEFAULT_USRNAME

	# ----------------------------------------------------------------------
	@staticmethod
	def addCommand(idx, name, cmd):
		command = Command(idx, name, cmd)
		Command.commands[name] = command
		Command.commands[idx]  = command
		return command

	# ----------------------------------------------------------------------
	@staticmethod
	def clearCommands():
		Command.commands.clear()

	# ----------------------------------------------------------------------
	@staticmethod
	def get(name):
		return Command.commands.get(name)

	# ----------------------------------------------------------------------
	# Load command configuration
	# From 1..99 system configuration and from 100.. user configurations
	# ----------------------------------------------------------------------
	@staticmethod
	def loadConfig(config):
		# Load all commands
		Command.clearCommands()
		i = 1
		while True:
			key = "cmd.%d"%(i)
			try: value = config.get(Command.section, key)
			except:
				if i<100:
					i=100
					continue
				else:
					break

			name,prg = value.split(":")
			cmd = Command.addCommand(i,name,prg)
			for a in range(1,100):
				key = "arg.%d.%d"%(i,a)
				try:
					cmd.addArgument(config.get(Command.section, key))
				except: break
			for j in range(1,100):
				key = "inp.%d.%d"%(i,j)
				try: cmd.addInput(config.get(Command.section, key))
				except: break
			i += 1

	# ----------------------------------------------------------------------
	@staticmethod
	def saveConfig(config):
		pass

#===============================================================================
# Info in project
#===============================================================================
class Info:
	# ----------------------------------------------------------------------
	def __init__(self, name):
		self.name   = name
		self.type   = "Info"
		self.prefix = ""
		self.var    = {}

	# ----------------------------------------------------------------------
	# Write information
	# ----------------------------------------------------------------------
	def write(self, f):
		if not self.var: return
		Input.utfWrite(f,"\n%s# %s\n"%(self.prefix, self.type))
		Input.utfWrite(f,"%s%s: %s\n"%(self.prefix, self.type, self.name))
		for k in sorted(self.var.keys()):
			v = self.var[k]
			try:
				v = str(v)
			except UnicodeEncodeError:
				v = str(v).encode("utf-8")
			if v != "": Input.utfWrite(f,"%s\t%s: %s\n" % (self.prefix, k,v))
		Input.utfWrite(f,"%sEnd\n"%(self.prefix))

	# ----------------------------------------------------------------------
	# Read geometry from file.
	# ----------------------------------------------------------------------
	def read(self, f):
		endstr = "End"+self.type
		while True:
			(tag, val) = _readVar(f)
			if tag is None: break
			if tag==endstr or tag=="End": break
			# Add in dictionary
			if val != "":
				if val in ("True","False"):
					self.var[tag] = val=="True"
				else:
					try: val = int(val)
					except: pass
					self.var[tag] = val
			else:
				self.var[tag] = _readBlock(f)

	# ----------------------------------------------------------------------
	def clone(self):
		c = self.__class__(self.name)
		c.copy(self)
		return c

	# ----------------------------------------------------------------------
	def copy(self, src):
		self.name = src.name
		self.type = src.type
		self.var  = src.var.copy()

	# ----------------------------------------------------------------------
	def hash(self):
		return hash(self.name) + \
		       hash(self.type)*3 + \
		       dicthash(self.var)
	__hash__ = hash

	# ----------------------------------------------------------------------
	def get(self, key, default=None):
		return self.var.get(key,default)
	__call__ = get

	# ----------------------------------------------------------------------
	def setdefault(self, key, value):
		return self.var.setdefault(key,value)

	# ----------------------------------------------------------------------
	def __getitem__(self, key):
		return self.var[key]

	# ----------------------------------------------------------------------
	def __setitem__(self, key, value):
		self.var[key] = value

	# ----------------------------------------------------------------------
	def __delitem__(self, key):
		del self.var[key]

	# ----------------------------------------------------------------------
	def items(self):  return list(self.var.items())
	def keys(self):   return list(self.var.keys())
	def values(self): return list(self.var.values())

	# ----------------------------------------------------------------------
	# Delete all keys starting with prefix
	# ----------------------------------------------------------------------
	def delPrefix(self, prefix):
		for name in list(self.var.keys()):
			if name.startswith(prefix):
				del self.var[name]

	# ----------------------------------------------------------------------
	# change name or set a default value
	# ----------------------------------------------------------------------
	def _change(self, old, new):
		if old in self.var:
			self.var[new] = self.var[old]
			del self.var[old]

	# ----------------------------------------------------------------------
	# set default value if doesn't exist
	# ----------------------------------------------------------------------
	def _default(self, name, default):
		if name not in self.var:
			self.var[name] = default

	# ----------------------------------------------------------------------
	def setUndo(self, undo):
		"""set all parameters and return undo as info"""
		undoinfo = (self.setUndo, self.clone())
		self.copy(undo)
		return undoinfo

	# ----------------------------------------------------------------------
	# @return a list of dependencies
	# ----------------------------------------------------------------------
	def depends(self):
		return []

	# ----------------------------------------------------------------------
	# @return a list of what files it provides
	# ----------------------------------------------------------------------
	def provides(self):
		return [self.name()]

	# ----------------------------------------------------------------------
	# Verify that all provided files are newer than the dependencies
	# ----------------------------------------------------------------------
	def verify(self, t0=0):
		for p in self.provides():
			try:
				pt = os.stat(p).st_mtime
			except:
				return False
			for d in self.depends():
				try:
					dt = os.stat(d).st_mtime
				except:
					return False
				if pt < dt or pt < t0:
					return False
		return True

	# ----------------------------------------------------------------------
	# create the 'provided' files
	# ----------------------------------------------------------------------
#	def make(self):
#		pass

	# ----------------------------------------------------------------------
	# clean all temporary and provided files
	# ----------------------------------------------------------------------
#	def clean(self):
#		pass

#===============================================================================
# UsrInfo
#===============================================================================
class UsrException(Exception): pass
class UsrInfo:
	def __init__(self, run, name=_DEFAULT_USRNAME):
		self._name = name
		self.run   = run
		self.type  = 'b'
		self.unit  =  0
		self.rules = []
		self.arg   = []		# Arguments for generic commands
		self.inp   = []		# Input parameters for generic commands
		self.cmd   = None

	# ----------------------------------------------------------------------
	# Write usr data information
	# ----------------------------------------------------------------------
	def write(self, f):
		Input.utfWrite(f,"\n\t# USRxxx data file: %s\n"%(self.name()))
		Input.utfWrite(f,"\tData: %s\n" % (self._name))
		Input.utfWrite(f,"\t\tUnit: %d\n" % (self.unit))
		Input.utfWrite(f,"\t\tType: %s\n" % (self.type))
		if len(self.type)==1:
			for i in self.rules:
				Input.utfWrite(f,"\t\tRule: %s,%s\n" % (i[0],i[1]))
		else:
			for a in self.arg:
				Input.utfWrite(f,"\t\tArg: %s\n" % (a))
			for i in self.inp:
				Input.utfWrite(f,"\t\tInp: %s\n" % (i))
		Input.utfWrite(f,"\tEnd\n")

	# ----------------------------------------------------------------------
	# Read usr data information
	# ----------------------------------------------------------------------
	def read(self, f):
		self.rules = []
		while True:
			(tag, val) = _readVar(f)
			if tag is None: break

			if   tag=="Unit": self.unit  = int(val)
			elif tag=="Type":
				self.type = val
				self.cmd = Command.get(self.type)
			elif tag=="Recursive": pass	# from V1
			elif tag=="Rule":
				# --- XXX removing old rule XXX ---
				# Added 21.03.13
				if val==r"+,_\I_\d\d\d\d\d_fort\.\U": continue
				# ----------------------------------
				self.rules.append(tuple(val.split(",")))
			elif tag=="Arg":
				self.arg.append(val)
			elif tag=="Inp":
				self.inp.append(val)
			elif tag in ("EndData","End"):
				break
			else:	# Add in dictionary
				raise UsrException("Unknown Data keyword "+tag)

	# ----------------------------------------------------------------------
	def set(self, unit=None, type=None):
		if unit is not None: self.unit = unit
		if type is not None:
			# Only if name is the default one change it
			if self._name == _usr_name.get(self.type, _DEFAULT_USRNAME):
				self._name = _usr_name.get(type,  _DEFAULT_USRNAME)
			self.type = type
			self.cmd  = Command.get(self.type)
		self.rules = _usr_rule[:]

	# ----------------------------------------------------------------------
	# return/set name
	# ----------------------------------------------------------------------
	def name(self, name=None):
		if name is not None:
			if name:
				self._name = name
			else:
				self._name = self.cmd.defaultName()
			return
		filename = self.expand(self._name)
		if "/" in filename: return filename
		d = self.run.getInputDirName()
		if d != "": filename = os.path.join(d,filename)
		if self.run.parent != "":
			parent = self.run.project.findParentRun(self.run)
			d = parent.getInputDirName()
			if d != "": filename = os.path.join(d,filename)
		return filename

	# ----------------------------------------------------------------------
	# Set argument, expand array if needed
	# ----------------------------------------------------------------------
	def setArg(self, idx, value):
		if idx >= len(self.arg):
			self.arg.extend([""]*(idx-len(self.arg)+1))
		self.arg[idx] = value

	# ----------------------------------------------------------------------
	def setInp(self, idx, value):
		if idx >= len(self.inp):
			self.inp.extend([""]*(idx-len(self.inp)+1))
		self.inp[idx] = value

	# ----------------------------------------------------------------------
	# Insert empty values for repeated statements starting from repeat
	# ----------------------------------------------------------------------
	def insertRepeat(self, pos, repeat):
		while repeat < len(self.cmd.inp):
			ctype,cdef,cname = self.cmd.inp[repeat]
			repeat += 1
			if ctype == "repeat": break
			self.inp.insert(pos, "")
			pos += 1

	# ----------------------------------------------------------------------
	# Delete all terms until "repeat"
	# ----------------------------------------------------------------------
	def deleteRepeat(self, pos, repeat):
		while repeat < len(self.cmd.inp):
			ctype,cdef,cname = self.cmd.inp[repeat]
			repeat += 1
			if ctype == "repeat": break
			del self.inp[pos]

	# ----------------------------------------------------------------------
	def __str__(self):
		return "UsrInfo{u=%d, t=%s, n=%s, files=%s}" \
			% (self.unit, self.type, self._name, str(self.rules))

	# ----------------------------------------------------------------------
	def __cmp__(self, b):
		return cmp(self.unit, b.unit)

	# ----------------------------------------------------------------------
	def clone(self):
		c = UsrInfo(self.run)
		c.copy(self)
		return c

	# ----------------------------------------------------------------------
	def copy(self, src):
		self._name = src._name
		self.run   = src.run	# will override the one from clone! and will break the RunInfo.clone()
		self.type  = src.type
		self.unit  = src.unit
		self.rules = src.rules[:]
		self.arg   = src.arg[:]
		self.inp   = src.inp[:]
		self.cmd   = src.cmd

	# ----------------------------------------------------------------------
	def hash(self):
		h = hash(self._name) ^ hash(self.type) ^ self.unit
		for o,r in self.rules:
			h ^= r
		h ^= hash(self.arg)
		h ^= hash(self.inp)
		return h
	__hash__ = hash

	# ----------------------------------------------------------------------
	# Return a list for possible file patterns based on the rule
	# By substituting all escape characters with * for possible file filters
	# ----------------------------------------------------------------------
	def fileFilters(self):
		rule = self._name.split("\\")
		items = []
		for i, r in enumerate(rule):
			if i==0:
				rest = r
			else:
				t = r[:1]
				if t=="F":
					items.append(self.run.project.name)
				elif t=="I":
					items.append(self.run.getInputBaseName())
				elif t=="U":
					items.append(str(self.unit))
				elif t=="t":
					items.append(self.type)
				elif t=="T":
					items.append(_usr_type[self.type])
				rest = r[1:]

			if rest!='':
				if rest.find("_")>=0:
					items.extend(rest.split('_'))
				elif rest.find(".")>=0:
					items.extend(rest.split('.'))
				else:
					items.append(rest)
		try:
			while True: items.remove('')
		except ValueError:
			pass

		for i, r in enumerate(items):
			if i==0:
				r += "*"
			elif i<len(items)-1:
				r = "*%s*"%(r)
			else:
				r = "*"+r
			items[i] = r

		filters = []
		for i in range(1,len(items)+1):
			for uc in bmath.xuniqueCombinations(items,i):
				f = " ".join(uc,'').replace("**","*")
				filters.append(f)

		return filters

	# ----------------------------------------------------------------------
	def typeLong(self):
		if len(self.type)==1:
			return _usr_type[self.type]
		else:
			return self.type

	# ----------------------------------------------------------------------
	# Return argument value for name or index
	# ----------------------------------------------------------------------
	def getArg(self, name):
		if isinstance(name,int):
			try:
				value = self.arg[name]
			except:
				value = ""
			ctag,ctype,cdef,cname = self.cmd.arg[name]
			if not value and cdef: value = self.expand(cdef)
			return value
		else:
			return None

	# ----------------------------------------------------------------------
	# expand escape characters in string
	# ----------------------------------------------------------------------
	def expand(self, s, esc=False, run=None):
		if run is None: run = self.run
		if esc:
			s = s.replace(r"\F", re.escape(run.project.name))
			s = s.replace(r"\I", re.escape(run.getInputBaseName()))
		else:
			s = s.replace(r"\F", run.project.name)
			s = s.replace(r"\I", run.getInputBaseName())
		s = s.replace(r"\U", str(self.unit))
		s = s.replace(r"\t", self.type)
		s = s.replace(r"\T", _usr_type.get(self.type,self.type))
		s = s.replace(r"\e", _usr_ext.get(self.type,"dat"))

		if self.cmd is not None:
			for i in range(1,10):
				target = r"\%d"%(i)
				if target in s:
					s = s.replace(target, self.getArg(i-1))


		if s and s[0]=="=":
			# evaluate expression and return result
			try:
				s = eval(s[1:], Input._globalDict, run.project.input.localDict)
			except:
				pass

		return s

	# ----------------------------------------------------------------------
	def compile(self, rule, run=None):
		rule = self.expand(rule, True, run)
		if rule[0]  != "^": rule = "^%s"%(rule)
		if rule[-1] != "$": rule = "%s$"%(rule)
		return re.compile(rule)

	# --------------------------------------------------------------------
	# Return processing program command
	# --------------------------------------------------------------------
	def command(self, cwd=None):
		if self.cmd is None:
			inp = self.depends(cwd)
			if len(inp) == 0: return None,None
			inp.append("")
			if cwd:
				inp.append(os.path.join(cwd,self.name()))
			else:
				inp.append(self.name())
			cmd = command(_fluka_programs.get(self.type, "us%ssuw"%(self.type)))
			return cmd,inp
		else:
			cmd = command(self.cmd.cmd)
			for i,(ctag,ctype,cdef,cname) in enumerate(self.cmd.arg):
				try:
					u = self.arg[i]
				except:
					u = ""
				if not u and cdef: u = self.expand(cdef)
				if u: cmd += " %s %s"%(ctag,u)

			inp = []
			i = 0	# rule index
			j = 0	# input index
			while i < len(self.cmd.inp):
				ctype,cdef,cname = self.cmd.inp[i]
				i += 1
				try:
					u = self.inp[j]
					j += 1
				except:
					u = ""
				if ctype == "repeat":
					if u != cdef:
						i  = int(cname)-1
						j -= 1
						continue
				if not u and cdef: inp.append(self.expand(cdef))
				inp.append(u)
			return cmd,inp

	# --------------------------------------------------------------------
	# Find files by executing the rules
	# --------------------------------------------------------------------
	def depends(self, cwd=None):
		files = []
		if cwd is None:
			rundir = self.run.getDir()
		else:
			rundir = cwd

		if self.cmd:
			for i,(ctag,ctype,cdef,cname) in enumerate(self.cmd.arg):
				if ctype not in ("file","filelist"): continue
				try:
					files.append(self.arg[i])
				except:
					if cdef: files.append(self.expand(cdef))

			i = 0	# rule index
			j = 0	# input index
			while i < len(self.cmd.inp):
				ctype,cdef,cname = self.cmd.inp[i]
				i += 1
				try:
					u = self.inp[j]
					j += 1
				except:
					u = ""
				if ctype == "repeat":
					if u != cdef:
						i  = int(cname)-1
						j -= 1
						continue
				if ctype not in ("file","filelist"): continue
				try:
					files.append(u)
				except:
					if cdef: files.append(self.expand(cdef))
		else:
			for op,rule in self.rules:
				if rule.find("*")>=0 or rule.find("\\")>=0:
					try:
						pattern = [(rundir, self.compile(rule))]
						for name in self.run.family:
							child = self.run.project.getRunByName(name, self.run.name)
							if child is not None:
								if cwd is None:
									childdir = child.getDir()
								else:
									childdir = cwd
								pattern.append((childdir, self.compile(rule, child)))
					except:
						say("Error:",sys.exc_info()[1])
						continue
					self._scanDir(files, op, pattern)
				else:
					self._doRule(files, op, rundir, rule)

		return files

	# --------------------------------------------------------------------
	# Scan Directory for Files
	# --------------------------------------------------------------------
	def _scanDir(self, files, op, patterns):
		msg = False
		for path, pattern in patterns:
			try: os.stat(path)
			except OSError: continue
#			if os.sep in pattern.pattern:
#				# Potentially there is a subdirectory in the pattern
#				for dirname, dirlist, filelist in os.walk(path):
#					for f in filelist:
#						fn = os.path.join(dirname[len(path)+1:],f)
#						if pattern.match(fn):
#							self._doRule(files, op, os.path.join(path,dirname), f)
#			else:
			# In the current path
			for f in Utils.listdir(path):
				if pattern.match(f):
					self._doRule(files, op, path, f)

	# --------------------------------------------------------------------
	# Add/Del a single file to the list
	# --------------------------------------------------------------------
	def _doRule(self, files, op, path, filename):
		if op == "-":	# Delete file
			# delete filename from list
			while True:
				try:
					files.remove(filename)
				except ValueError:
					return
		else:	# Add filename
			files.append(self.run.project.relativePath(
					os.path.join(path, filename)))

	# ----------------------------------------------------------------------
	# Verify that all provided files are newer than the dependencies
	# ----------------------------------------------------------------------
	def verify(self, t0=0):
		try:
			pt = os.stat(self.name()).st_mtime
		except:
			return False

		for d in self.depends():
			try:
				dt = os.stat(d).st_mtime
				if pt<dt or pt<t0:
					return False
			except:
				pass
		return True

	# ----------------------------------------------------------------------
	# Config load
	# ----------------------------------------------------------------------
	@staticmethod
	def loadConfig(config):
		del _usr_rule[:]	# clear the rules
		#items = config.items(Command.section)
		#if not items:
		#	for k,v in Project._usr_name.items():
		#		config.set(Command.section, k, v)
		#	i = 1
		#	for t,r in Project._usr_rule:
		#		config.set(Command.section, "rule.%d"%(i), "%s%s"%(t,r))
		#		i += 1
		for k,v in config.items(Command.section):
			# --- XXX  cleanup of rule: +_\I_\d\d\d\d\d_fort\.\U  XXX ---
			# To be cleaned latter
			# Added: 21.03.2013
			if v==r"+_\I_\d\d\d\d\d_fort\.\U": continue
			# --- XXX XXX XXX ---
			if k.startswith("rule."):
				try:
					idx = int(k[5:])-1
				except:
					idx = None
				addRule(v, idx)
			elif "." not in k:
				_usr_name[k] = v

	# ----------------------------------------------------------------------
	@staticmethod
	def saveConfig(config):
		pass

#===============================================================================
# Run Information Class
#===============================================================================
class RunException(Exception): pass
class RunInfo:
	section = "Batch"

	def __init__(self, project, name):
		self.project   = project
		self.name      = name	# Input name
		self.title     = ""	# Title of the run
		self.defines   = []	# override defines
		self.rnd       = 0	# -//- random seed
		self.primaries = 0	# -//- random No of primaries
		self.time      = 0      # -//- time limit
		self.exe       = ""	# -//- executable
		self.prev      = 0	# Previous cycle
		self.last      = 5	# Last cycle
		self.pid       = 0	# Process id when running
					# directory is fluka_%d % (pid)
		self.cycle     = 0
		self.status    = STATUS_NOT_RUNNING
		self.startRun  = 0	# time when run button was pressed
		self.attachTime= 0	# when the last attached was requested
		self.usrinfo   = []	# Usrinfos list

		self.handled   = 0
		self.remaining = 0
		self.startTime = 0
		self.timeperprim = 0
		self.initime   = 0
		self.maxtime   = MAX_TIME
		self.exclude   = None
		self.spawnName = ""	# Spawn pattern name (if not defined use default)

		self.parent    = ""	# parent name
		self.family    = []	# other members of the family
		self.index     = 0

		self._runDir   = None	# running directory
		self._lastOut  = None	# last output file opened
		self._lastPos  = 0	# last output file position scanned

	# ----------------------------------------------------------------------
	# Write run information
	# ----------------------------------------------------------------------
	def write(self, f):
		f.write("\n# Run information\n")
		Input.utfWrite(f,"Run: %s\n" % (self.name))
		if self.title != "":      Input.utfWrite(f,"\tTitle:    %s\n" % (self.title))
		for n,v in self.defines:
			v = v.strip()
			if v!="":
				f.write("\tDefine:   %s=%s\n"%(n,v))
			else:
				f.write("\tDefine:   %s\n"%(n))
		if self.rnd > 0:    f.write("\tRnd:      %d\n" % (self.rnd))
		if self.primaries>0:f.write("\tStart:    %d\n" % (self.primaries))
		if self.time > 0:   f.write("\tTime:     %d\n" % (self.time))
		if self.exe != "":  Input.utfWrite(f,"\tExe:      %s\n" % (self.exe))
		f.write("\tPrev:     %d\n" % (self.prev))
		f.write("\tLast:     %d\n" % (self.last))
		Input.utfWrite(f,"\tStatus:   %d\n" % (self.status))
		f.write("\tPid:      %d\n" % (self.pid))
		f.write("\tStartRun: %d\n" % (self.startRun))
		if self.spawnName:      Input.utfWrite(f,"\tSpawnName: %s\n" % (self.spawnName))
		if self.parent:		Input.utfWrite(f,"\tParent:    %s\n" % (self.parent))
		if len(self.family) > 0: f.write("\tFamily:  %s\n" \
					% (" ".join(self.family)))
		for u in self.usrinfo:
			u.write(f)
		f.write("End\n")

	# ----------------------------------------------------------------------
	# Read run from file. The file pointer should be exactly
	# after the Run: tag
	# ----------------------------------------------------------------------
	def read(self, f):
		child = False
		while True:
			(tag, val) = _readVar(f)
			if tag is None: break
			if   tag=="Title":    self.title     = val
			elif tag=="Rnd":      self.rnd       = int(val)
			elif tag=="Start":    self.primaries = int(val)
			elif tag=="Time":     self.time      = int(val)
			elif tag=="Exe":      self.exe       = val
			elif tag=="Prev":     self.prev      = int(val)
			elif tag=="Last":     self.last      = int(val)
			elif tag=="Status":   self.status    = int(val)
			elif tag=="Pid":      self.pid       = int(val)
			elif tag=="StartRun": self.startRun  = int(val)
			elif tag=="SpawnName":self.spawnName = val
			elif tag=="Parent":   self.parent    = val
			elif tag=="Child":    child          = bool(int(val))
			elif tag=="Defines":
				self.defines  = [(n,"") for n in val.split()]	# deprecated
			elif tag=="Define":
				e = val.find("=")
				if e>0:
					n = val[:e]
					v = val[e+1:].strip()
				else:
					n = val
					v = ""
				self.defines.append((n,v))
			elif tag=="Data":
				usr = UsrInfo(self, val)
				usr.read(f)
				self.usrinfo.append(usr)
			elif tag=="Family": self.family  = val.split()
			elif tag in ("EndRun","End"):
				break
			else:
				raise RunException("Unknown Run keyword "+tag)

		if child and self.parent == "":
			# Try to find parent run for old projects
			for parent in self.project.runs:
				if self.name in parent.family:
					self.parent = parent.name
					break

	# ----------------------------------------------------------------------
	def clone(self):
		c = RunInfo(self.project, self.name)
		c.copy(self)
		return c

	# ----------------------------------------------------------------------
	def copy(self, src):
		self.project     = src.project
		self.name        = src.name
		self.title       = src.title
		self.defines     = src.defines[:]
		self.rnd         = src.rnd
		self.primaries   = src.primaries
		self.time        = src.time
		if src.name == DEFAULT_INPUT:
			self.exe = src.project.exe
		else:
			self.exe = src.exe
		self.prev        = src.prev
		self.last        = src.last
		self.pid         = src.pid
		self.cycle       = src.cycle
		self.status      = src.status
		self.startRun    = src.startRun
		self.attachTime  = src.attachTime
		self.handled     = src.handled
		self.remaining   = src.remaining
		self.startTime   = src.startTime
		self.timeperprim = src.timeperprim
		self.exclude     = src.exclude
		self.usrinfo     = [x.clone() for x in src.usrinfo]
		for u in self.usrinfo: u.run = self	# update run pointer

	# ----------------------------------------------------------------------
	def getInputName(self):
		if self.name == DEFAULT_INPUT:
			return self.project.inputName
		else:
			return self.name

	# ----------------------------------------------------------------------
	def getInputName2(self):
		if self.name == DEFAULT_INPUT:
			if self.project.name == "":
				return "<untitled>"
			else:
				return "<%s>"%(self.project.inputName)
		else:
			return self.name

	# ----------------------------------------------------------------------
	def getInputBaseName(self):
		return os.path.basename(self.getInputName())

	# ----------------------------------------------------------------------
	def getInputDirName(self):
		return os.path.dirname(self.getInputName())

	# ----------------------------------------------------------------------
	# Get default directory
	# ----------------------------------------------------------------------
	def getDir(self):
		subdir = os.path.dirname(self.getInputName())
		if len(subdir)>0:
			if subdir[0] == os.path.sep:
				return subdir
			else:
				if self.parent != "":
					parent = self.project.findParentRun(self)
					if parent is not None:
						return os.path.join(parent.getDir(), subdir)

				return os.path.join(self.project.dir, subdir)

		elif self.parent != "":
			parent = self.project.findParentRun(self)
			if parent is not None:
				return parent.getDir()

		return self.project.dir

	# ----------------------------------------------------------------------
	# Get temporary running directory if any otherwise default
	# ----------------------------------------------------------------------
	def _getRunDir(self):
		if self._runDir is not None:
			return self._runDir
		else:
			return self.getDir()

	# ----------------------------------------------------------------------
	def resetRunDir(self):
		self._runDir = None

	# ----------------------------------------------------------------------
	# setExecutable: always as absolute
	# ----------------------------------------------------------------------
	def setExecutable(self, fn):
		if fn != "":
			path = os.path.dirname(fn)
			if path != self.project.dir:
				self.exe = fn		# Absolute path
			else:
				self.exe = os.path.basename(fn)
		else:
			self.exe = fn

	# ----------------------------------------------------------------------
	def defaultDefines(self):
		self.defines = []
		for card in self.project.input.cardsSorted("#define",2):
			if card._indent != 0: continue
			if card.ignore(): continue
			sdum = card.sdum()
			for n,w in self.defines:
				if n==sdum: break
			else:
				self.defines.append((sdum, card.what(1)))

	# ----------------------------------------------------------------------
	# return a tuple to spawn a fluka command
	# ----------------------------------------------------------------------
	def command(self, cwd=None):
		"""Return a tuple containing the command to spawn a fluka run"""

		cmd = []
		submit = _submit_command.get(self.project.submit,"")
		if submit: cmd.extend(shlex.split(submit))
		cmd.append(command("rfluka"))

		exe = self.exe
		if self.name == DEFAULT_INPUT:
			if self.project.exe != "": exe = self.project.exe
		else:
			if self.parent != "": # Child
				parent = self.project.findParentRun(self)
				if parent is not None:
					if parent.name == DEFAULT_INPUT: exe = self.project.exe
					elif parent.exe != "": exe = parent.exe

		if exe:
			if cwd is not None and exe.find("/") < 0:
				exe = os.path.join(self.project.dir, exe)
			cmd.append("-e")
			cmd.append(exe)

		if self.prev > 0:
			cmd.append("-N")
			cmd.append(str(self.prev))
		cmd.append("-M")
		cmd.append(str(self.last))
		cmd.append(self.getInputBaseName())

		return tuple(cmd)

	# ----------------------------------------------------------------------
	# Find possible fluka_XXXX directories where the run takes place
	# ----------------------------------------------------------------------
	def findFlukaDirs(self, check_core=True):
		# Scan all subdirs with pattern "fluka_\d+"
		# for a file named "input\d\d\d.out"
		# check for a core.\d+ file
		found  = []
		inpname= re.escape(self.getInputBaseName())
		rundir = self._getRunDir()
		patout = re.compile(r"^%s(\d\d\d).out$" % (inpname))
		patlog = re.compile(r"^%s(\d\d\d).log$" % (inpname))

		try:
			os.stat(rundir)
		except:
			return found

		for d in Utils.listdir(rundir):
			if _PAT_DIR.match(d):	# Directory found
				fullpath = os.path.join(rundir,d)
				relpath  = self.project.relativePath(fullpath)
				try:   s = os.stat(fullpath)
				except: continue
				if not stat.S_ISDIR(s.st_mode): continue
				#if check_core and \
				#	s.st_mtime < self.startRun - timeout:
				#		continue
				if self.exclude and d in self.exclude: continue
				# Check files in directory
				out  = None
				log  = None
				core = None
				for f in Utils.listdir(fullpath):
					# ignore links
					try: s = os.lstat(os.path.join(fullpath,f))
					except: continue
					if stat.S_ISLNK(s.st_mode):
						continue
					elif check_core and (f=="core" or _PAT_CORE.match(f)):
						core = (relpath,None)
					elif patout.match(f):
						out = (relpath,f)
					elif patlog.match(f):
						log = (relpath,f)
				if core and (out or log):
					found.append(core)	# core.xxx file
				elif out is not None:	# Prefer out file
					found.append(out)
				elif log is not None:	# Otherwise log file
					found.append(log)
		return found

	# ----------------------------------------------------------------------
	# Start the run and try to attach
	# ----------------------------------------------------------------------
	def start(self, cwd=None, defines=None, log=None):
		if cwd is not None:
			self._runDir = cwd
		else:
			self.resetRunDir()
			cwd = self.getDir()
		if log is None: log = say

		# Check if directory exist
		try:
			os.stat(cwd)
		except OSError:
			os.makedirs(cwd)

		# Check last random number if exist
		self.status = STATUS_NOT_RUNNING
		if self.prev > 0:
			inpname = self.getInputBaseName()
			try:
				s = os.stat("%s/ran%s%03d" \
					% (cwd, inpname, self.prev+1))
			except OSError:
				raise RunException(
					"Error running input \"%s\". Random number from "\
					"previous cycle %s/ran%s%03d not found" \
					%(inpname, cwd, inpname, self.prev+1))

		# Prepare input file
		rc = 0
		errors = self.project.writeInput(self, cwd, defines)
		if errors:
			for card,err in errors:
				log(card)
				log(err)
			rc = len(errors)

		# Check executable
		try:
			os.stat(self.project.executable(self))
		except:
			raise RunException("Fluka executable \"%s\" not found" % (self.project.executable(self)))

		# Create run command
		if self.project.submit == "null": return 0
		cmd = self.command(cwd)
		try:
			os.stat(cmd[0])
		except OSError:
			raise RunException("Cannot find executable %s\n" \
				"Please verify that the FLUKA directory\n" \
				" %s=%s\nor the submit command are correct!" \
				% (cmd[0], flukaVar, flukaDir))

		stdout = open("%s.out"%(self.getInputBaseName()), "w")
		stdout.write("Dir: %s\n"%(cwd))
		stdout.write("Cmd: %s\n"%(" ".join(list(cmd))))
		log("Dir: %s"%(cwd))
		log("Cmd: %s"%(" ".join(list(cmd))))

		# Exclude existing directories
		self.exclude = [d for d,f in self.findFlukaDirs(check_core=False)]

		# spawn in the background
		self._popen = subprocess.Popen(cmd,
					cwd=cwd,
					stdout=stdout,
					stderr=subprocess.STDOUT,
					close_fds=True,
					preexec_fn=os.setpgrp)	# start a new session, detach from the father

		# Change status to reflect changes
		self.status     = STATUS_WAIT2ATTACH	# Wait for the second update
		self.startRun   = time.time()
		self.attachTime = self.startRun
		self.initime    = 0
		self.maxtime    = MAX_TIME
		self._lastOut   = None			# last output file opened
		self._lastPos   = 0			# last output file position scanned
		return rc

	# ----------------------------------------------------------------------
	# Attach to a process
	# ----------------------------------------------------------------------
	def attach(self):
		# Set default info
		if time.time() - self.attachTime > timeThreshold:
			self.status = STATUS_TIMEOUT
		self.pid     = 0
		self.cycle   = 0
		self.handled = 0

		# Scan all subdirs with pattern "fluka_\d+"
		# for a file named "input\d\d\d.out"
		# check for a core.\d+ file
		inpname = self.getInputBaseName()
		try:
			found = self.findFlukaDirs()
		except OSError:
			found = []

		# Check if run is already finished
		if len(found) == 0:
			try:
				tm = os.stat("%s/%s%03d.out"%(self._getRunDir(),inpname,self.last)).st_mtime
				if tm >= self.startRun:
					self.status = STATUS_FINISHED
			except:
				pass
			return None # Nothing found

		# Filter out all core files
		found = [ (d,f) for d,f in found if f ]

		if len(found) == 0:
			self.status = STATUS_FINISHED_ERR
			return None
		elif len(found) > 1:
			return found	# Multiple choices ask the user

		# Else take the first one
		self.status = STATUS_RUNNING
		self.pid    = int(_PAT_DIR2.match(found[0][0]).group(2))
		pat	    = re.compile(r"^%s(\d\d\d).(out|log)$" % (re.escape(inpname)))
		m = pat.match(found[0][1])
		if m is None:
			self.status = STATUS_FINISHED_ERR
			m = pat.match(found[0][1])
		if m is not None:
			self.cycle  = int(m.group(1))
		self.attachTime = time.time()
		return found[0]

	# ----------------------------------------------------------------------
	# Remove all generated files and return True on error, False otherwise
	# ----------------------------------------------------------------------
	def clean(self, exclude=None, log=None):
		if log is None: log = say
		inpname = self.getInputBaseName()
		p  = re.compile(r"^%s\d\d\d[._].+$" % (inpname))
		rp = re.compile(r"^ran%s\d\d\d$" % (inpname))
		rundir = self._getRunDir()

		try:
			os.stat(rundir)
		except:
			return False

		errors = False
		for fn in Utils.listdir(rundir):
			if exclude:
				if isinstance(exclude,list):
					skip = False
					for pattern in exclude:
						if fnmatch(fn,pattern):
							skip = True
							break
					if skip: continue
				elif fnmatch(fn,exclude): continue
			if p.match(fn) or rp.match(fn):
				say("Remove: %s"%(fn))
				os.remove(os.path.join(rundir,fn))
			elif fn.startswith(_FLUKA_CYCLE):
				dfn = os.path.normpath(os.path.join(rundir,fn))
				try: s = os.lstat(dfn)
				except OSError: continue
				if not stat.S_ISDIR(s.st_mode): continue
				found = False
				inpfile = "%s.inp"%(inpname)
				for fsub in Utils.listdir(dfn):
					if  p.match(fsub) or \
					   rp.match(fsub) or \
					   fsub==inpfile:
						found = True
						break
				if found:
					say("Remove folder: "+fn)
					try:
						shutil.rmtree(dfn)
					except OSError:
						say(sys.exc_info()[1])
						errors = True
		return errors

	# ----------------------------------------------------------------------
	# refresh status if running
	# @param output	if True check the running output for # events handled
	# ----------------------------------------------------------------------
	def refresh(self, output=False):
		inpname = self.getInputBaseName()

		# Check directory if it exists
		rundir = "%s/fluka_%d" % (self._getRunDir(), self.pid)
		try:
			os.stat(rundir)
		except:
			self.status = STATUS_FINISHED
			self.pid    = 0
			return

		# find if cycle still exist
		found = False
		for self.cycle in range(self.cycle, self.last+1):
			try:
				fnout = "%s/%s%03d.out" % (rundir, inpname, self.cycle)
				outstat = os.stat(fnout)
				found = True
				break
			except OSError:
				pass
			try:
				logstat = os.stat("%s/%s%03d.log" \
					% (rundir, inpname, self.cycle))
				found = True
				break
			except OSError:
				pass

		# Check if a core.### file exists
		maxtime = 0
		try:
			for f in Utils.listdir(rundir):
				try: s = os.lstat(os.path.join(rundir,f))
				except: continue
				if stat.S_ISLNK(s.st_mode): continue
				maxtime = max(maxtime, s.st_mtime)
				if f=="core" or _PAT_CORE.match(f):
					self.status = STATUS_FINISHED_ERR
					self.pid    = 0
					return
		except OSError:
			self.status = STATUS_FINISHED
			self.pid    = 0
			return

		# check maximum time of all files
		if time.time() - maxtime > timeThreshold:
			# Hmm it appears dead!
			self.status = STATUS_TIMEOUT
			return

		# No output, log file and no core try to re-attach
		# FIXME Should have a counter, how many times to attach
		if not found:
			self.status    = STATUS_WAIT2ATTACH
			self.pid       = 0
			self.cycle     = 0
			self.handled   = 0
			self.startTime = 0
			return

		# Running
		self.status = STATUS_RUNNING
		if not output: return

		# Find progress
		if fnout != self._lastOut:
			# New file reset information
			self._lastOut    = fnout
			self._lastPos    = 0
			self.handled     = 0
			self.remaining   = 0
			self.timeperprim = 0.0
		try:
#			f = open(self._lastOut,"r")
			f = open(self._lastOut,"r",errors='ignore')
		except:
			# Rest information
			self._lastOut    = None
			self._lastPos    = 0
			self.handled     = 0
			self.remaining   = 0
			self.timeperprim = 0.0
			return

		if self._lastPos: f.seek(self._lastPos)
		nextseed = prev = ""
		offset=0
		for line in f:
			offset += len(line)
			if line.startswith(" NEXT SEEDS:"):
				nextseed = prev
#				self._lastPos = f.tell()
				self._lastPos = offset
			elif self.initime == 0:
				if line.startswith(" Total time used for initialization:"):
					try: self.initime = float(line.split()[-2])
					except: pass
				elif line.startswith(" Maximum cpu-time allocated for this run:"):
					try: self.maxtime = float(line.split()[-2])
					except: pass
			else:
				prev = line
		f.close()

		line = nextseed.split()
		if len(line)>0:
			try:
				self.handled     = int(line[0])
				self.remaining   = int(line[1])
				self.timeperprim = float(line[3])
			except:
				pass

		# Get the starting time from when the fort.2 file was modified
		try:
			self.startTime = os.lstat("%s/fort.1"%(rundir)).st_mtime
		except:
			self.startTime = 0

	# ----------------------------------------------------------------------
	def _stop(self, stopfile):
		def _stopRun(run):
			# Check status, if not running start the run
			if run.status != STATUS_RUNNING: return
			# Create fluka.stop or rfluka.stop files
			try:
				f = open("%s/fluka_%d/%s" \
					% (run._getRunDir(), run.pid, stopfile),"w")
				f.close()
			except:
				pass

		_stopRun(self)
		if self.family:
			for child in self.family:
				_stopRun(self.project.getRunByName(child, self.name))

	# ----------------------------------------------------------------------
	def stopRun(self): _stop("rfluka.stop")
	def stopCycle(self): _stop("fluka.stop")

	# ----------------------------------------------------------------------
	# Scan directory to find last random number seed for a possible cont.
	# ----------------------------------------------------------------------
	def lastRun(self):
		ranpat = re.compile(r"^ran%s(\d\d\d)$" % (re.escape(self.getInputBaseName())))
		# Scan run directory...
		last = 1
		for fn in Utils.listdir(self._getRunDir()):
			m = ranpat.match(fn)
			if m:
				no = int(m.group(1))
				if last < no: last = no
		return last-1

	# ----------------------------------------------------------------------
	def hash(self):
		h = 0
		for i in self.family: h ^= hash(i)
		for n,v in self.defines: h ^= hash(n)^hash(v)
		return hash(self.name) ^ hash(self.title) \
		  ^ self.rnd ^ hash(self.primaries) ^ hash(self.time) \
		  ^ hash(self.exe) ^ self.prev \
		  ^ self.status ^ self.pid ^ int(self.startRun) \
		  ^ hash(self.parent) ^ h
	__hash__ = hash

	# ----------------------------------------------------------------------
	# Load commands
	# ----------------------------------------------------------------------
	@staticmethod
	def loadConfig(config):
		global _submit_command, _DEFAULT_RUNNAME
		from tkFlair import getStr

		_DEFAULT_RUNNAME = getStr(RunInfo.section, "spawnname",	_DEFAULT_RUNNAME)

		_submit_command.clear()
		for tag, value in config.items(RunInfo.section):
			_submit_command[tag] = value

	# ----------------------------------------------------------------------
#	@staticmethod
#	def saveConfig(config):
#		for name,material in RunInfo.globalMaterials.items():
#			for var in ("color", "alpha"):
#				config.set(MaterialInfo.section, "%s.%s"%(var,name), str(material[var]))

#===============================================================================
# Plot Information Class
#===============================================================================
class PlotInfo(Info):
	def __init__(self, name):
		Info.__init__(self, name)
		self.format = ""	# Format of plot (file extension)
		self.type   = ""	# Type of the plot
		self.prefix = "\t"
		self.time   = 0		# last plotting time

	# ----------------------------------------------------------------------
	# Write run information
	# ----------------------------------------------------------------------
	def write(self, f):
		Input.utfWrite(f,"\n# %s plot \"%s\"\n"%(self.type,self["title"]))
		Input.utfWrite(f,"Plot: %s\n"  % (self.name))
		f.write("\tFormat: %s\n" % (self.format))
		f.write("\tType:   %s\n" % (self.type))
		for k in sorted(self.var.keys()):
			v = self.var[k]
			try:
				v = str(v)
			except UnicodeEncodeError:
				v = str(v).encode("utf-8")
			if v != "":
				if k == "commands":
					Input.utfWrite(f,"\t%s:\n%s\f\n" % (k,v))
				else:
					Input.utfWrite(f,"\t%s: %s\n" % (k,v))
		f.write("End\n")

	# ----------------------------------------------------------------------
	# Read run from file. The file pointer should be exactly
	# after the Run: tag
	# ----------------------------------------------------------------------
	def read(self, f):
		while True:
			(tag, val) = _readVar(f)
			if tag is None: break

			if tag=="Type":     self.type   = val
			elif tag=="Format": self.format = val
			elif tag in ("EndPlot","End"):
				break
			else:
				tag = tag.lower()
				tag = _plot_translate.get(tag, tag)
				if tag=="commands":
					self.var["commands"] = _readBlock(f)
				elif val in ("True","False"):
					self.var[tag] = val=="True"
				else:	# Add in dictionary
					self.var[tag] = val

	# ----------------------------------------------------------------------
	def copy(self, src):
		Info.copy(self, src)
		self.format = src.format

	# ----------------------------------------------------------------------
	def hash(self):
		return Info.hash(self) + hash(self.format)*5

	# ----------------------------------------------------------------------
	def __getitem__(self, key):
		return self.var.get(key,"")

	# ----------------------------------------------------------------------
	def correct(self, version):
		# XXX cleanup of material colors
		# Legacy code correction
		# FIXME we should move it to the geometry section
		# move to project var
		# since we can delete within the loop

#			for name, color in self.var.items():
#				if name.startswith("material."):
#					self.var[name] = color

		# Change from 1->2
		if version < 2:
			if self.type == "USRBIN":
				# move [xyz][min|max] -> b\0
				self._change("xmin","bxmin")
				self._change("xmax","bxmax")
				self._change("ymin","bymin")
				self._change("ymax","bymax")
				self._change("zmin","bzmin")
				self._change("zmax","bzmax")
				self._change("log", "cblog")
				self._default("cbtics", 1)
				lc = self["lc.0"]
				if lc=="":
					try:
						lt = int(self["lt.0"])
						if lt<=0: lc = __gnucolors[0]
						else: lc = __gnucolors[lt%7]
						self["lc.0"] = lc
						del self.var["lt.0"]
					except:
						pass

			elif self.type == "USR-1D":
				for i in range(int(self.get("ndetectors",0))):
					lc = self["lc.%d"%(i)]
					if lc == "":
						try:
							lt = int(self["lt.%d"%(i)])
							if lt<=0: lc = __gnucolors[0]
							else: lc = __gnucolors[lt%7]
							self["lc.%d"%(i)] = lc
							del self.var["lt.%d"%(i)]
						except:
							pass

					# It should be consistent with the Plot.VALUE
					y = int(self["y.%d"%(i)])
					if 1<=y<=4:	# Y*<X>, Y*<Xgeo>, Y*Xl, Y*Xh
						y=1	# convert to isoolethargic
					elif y==5:
						y=2	# convert to Y*DX
					else:
						y=0
					self["y.%d"%(i)] = y

			elif self.type == "USR-2D":
				self._change("log", "cblog")

#===============================================================================
class GeometryLayerInfo(Info):
	def __init__(self, name, isglobal=False):
		Info.__init__(self, name)
		self.type    = "Layer"
		self.prefix  = "\t"
		self._global = isglobal

	# ----------------------------------------------------------------------
	def read(self,f):
		Info.read(self,f)
		#############################
		# FIXME
		# Rename Colorband -> Palette
		#############################
		for key,value in list(self.var.items()):
			if key.startswith("Colorband"):
				newkey = "Palette" + key[9:]
				self.var[newkey] = value
				del self.var[key]
			elif key == "Show.label" and not value:
				self.var[key] = "Region"
		#############################

	# ----------------------------------------------------------------------
	def copy(self, src):
		Info.copy(self, src)
		self._global = src._global

#===============================================================================
class GeometryObjectInfo(Info):
	def __init__(self, name):
		Info.__init__(self, name)
		self.type   = "Object"

#===============================================================================
# Geometry Info
#===============================================================================
class GeometryInfo:
	def __init__(self):
		self.var    = {}
		self.layers = {}
		self.layers = _globalLayers.copy()	# Make a copy of global layers
		#self.objects = []

	# ----------------------------------------------------------------------
	# Write geometry information
	# ----------------------------------------------------------------------
	def write(self, f):
		Input.utfWrite(f,"\n# Geometry Information\n")
		Input.utfWrite(f,"Geometry:\n")
		for k in sorted(self.var.keys()):
			v = self.var[k]
			try:
				v = str(v)
			except UnicodeEncodeError:
				v = str(v).encode("utf-8")
			if v != "": Input.utfWrite(f,"\t%s: %s\n" % (k,v))
		for layer in list(self.layers.values()):
			if not layer._global: layer.write(f)
		#for object in self.objects:
		#	object.write(f)
		Input.utfWrite(f,"End\n")

	# ----------------------------------------------------------------------
	# Read geometry from file.
	# ----------------------------------------------------------------------
	def read(self, f):
		while True:
			(tag, val) = _readVar(f)
			if tag is None: break
			if tag in ("EndGeometry","End"): break
			elif tag == "Layer":
				# Parse layer
				layer = GeometryLayerInfo(val)
				layer.read(f)
				self.layers[layer.name] = layer

			elif tag == "Object":
				# Parse object
				object = GeometryObjectInfo(val)
				object.read(f)
				#self.objects.append(object)

			# Add in dictionary
			elif val != "":
				if val in ("True","False"):
					self.var[tag] = val=="True"
				else:
					self.var[tag] = val

			else:
				self.var[tag] = _readBlock(f)

	# ----------------------------------------------------------------------
	def hash(self):
		h = dicthash(self.var)
		for layer in list(self.layers.values()):
			h += layer.hash()
		return h
	__hash__ = hash

	# ----------------------------------------------------------------------
	def __getitem__(self, key):
		return self.var[key]

	# ----------------------------------------------------------------------
	def __setitem__(self, key, value):
		self.var[key] = value

	# ----------------------------------------------------------------------
	def get(self, key, default=None):
		return self.var.get(key,default)

	# ----------------------------------------------------------------------
	def setdefault(self, key, value):
		return self.var.setdefault(key,value)

	# ----------------------------------------------------------------------
	# Scan all layers to see if something has moved out from global
	# or the opposite
	# ----------------------------------------------------------------------
	def correctGlobalLayers(self):
		for layer in list(self.layers.values()):
			# Layer was in global before
			if layer in list(_globalLayers.values()):
				# if it was moved out from globals
				if not layer._global:
					# FIXME Warning if other flair are open
					# they will keep a duplicate!
					del _globalLayers[layer.name]

			# else if wasn't in global but now it has to move there
			elif layer._global:
				_globalLayers[layer.name] = layer

	# ----------------------------------------------------------------------
	def layersList(self):
		lst = [x.name for x in list(self.layers.values()) if x.get("_system") and int(x.get("_order",999))<2]
		lst.sort()
		l = [x.name for x in list(self.layers.values()) if x.get("_system") and int(x.get("_order",999))>=2]
		l.sort()
		lst.extend(l)
		l = [x.name for x in list(self.layers.values()) if not x.get("_system")]
		l.sort()
		lst.extend(l)
		return lst

#===============================================================================
# DicomInfo
#===============================================================================
class DicomInfo(Info):
	def __init__(self, name):
		Info.__init__(self, name)
		self.type   = "Dicom"
		self.files  = None

	# ----------------------------------------------------------------------
	def correct(self):
		# correct uid -> name
		s = self.get("uid")
		if s is not None and s != "":
			self.name = s
		if "directory" not in self.var:
			self["directory"] = "."

	# ----------------------------------------------------------------------
	def clone(self):
		c = Info.clone(self)
		c.files = self.files[:]
		return c

#===============================================================================
# PetInfo
#===============================================================================
class PetInfo(Info):
	def __init__(self, name):
		Info.__init__(self, name)
		self.type   = "Pet"

#===============================================================================
class MaterialInfo(Info):
	globalMaterials = {}
	section = "Material"

	def __init__(self, name):
		Info.__init__(self, name)
		self.type   = "Material"
		self.random()
		# default variables
		self["specular"] =  0.0
		self["shine"]    = 10.0
		self["fuzz"]     =  0.0
		self["ior"]      =  1.0
		self["alpha"]    =  0

	# ----------------------------------------------------------------------
	def random(self):
		self["color"] = "#%06X"%(bmath.rndColor(hash(self.name)))
		return self["color"]

	# ----------------------------------------------------------------------
	def color32(self):
		try:
			return int(self.var["color"][1:],16)
		except:
			return 0x707070

	# ----------------------------------------------------------------------
	def isGlobal(self):
		return self in list(MaterialInfo.globalMaterials.values())

	# ----------------------------------------------------------------------
	# Check if a color is valid or not
	# ----------------------------------------------------------------------
	@staticmethod
	def valid(color):
		if len(color)!=7: return False
		if color[0]!='#': return False
		try: int(color[1:],16)
		except: return False
		return True

	# ----------------------------------------------------------------------
	# Return global material.
	# If it doesn't exist create it
	# ----------------------------------------------------------------------
	@staticmethod
	def getGlobal(name):
		mat = MaterialInfo.globalMaterials.get(name)
		if mat is None:
			MaterialInfo.globalMaterials[name] = mat = MaterialInfo(name)
		return mat

	# ----------------------------------------------------------------------
	# Load materials
	# ----------------------------------------------------------------------
	@staticmethod
	def loadConfig(config):
		for tag, value in config.items(MaterialInfo.section):
			var, name = tag.split(".")
			material  = MaterialInfo.getGlobal(name)
			material[var] = value

	# ----------------------------------------------------------------------
	@staticmethod
	def saveConfig(config):
		for name,material in list(MaterialInfo.globalMaterials.items()):
			for var in ("color", "alpha"):
				config.set(MaterialInfo.section, "%s.%s"%(var,name), str(material[var]))

#===============================================================================
# Project class holder
#===============================================================================
class Project:
	def __init__(self, filename=None):
		"""Initialize a fluka project"""
		if filename is not None:
			self.load(filename)
		else:
			self.clear()

	# ----------------------------------------------------------------------
	# clear
	# ----------------------------------------------------------------------
	def clear(self):
		"""Set to default all parameters of fluka project"""
		# File Names
		self.version    = __version__
		self.name       = ""
		self.projFile   = ""
		self.inputFile  = ""
		self.inputName  = ""
		self.exe        = ""
		self.dir        = os.path.realpath(os.getcwd())
		self.tabs       = [0.10, 0.32, 0.56, 0.81, 1.0]
		self.page       = "Flair"	# Last page visited

		# Project information
		self.title      = ""
		self.notes      = ""
		self.calculator = ""
		self.submit     = "*Default"

		self.projectModified = False
		self.inputModified   = False

		self.input      = Input.Input()
		self.debugList  = []

		# Fortran var
		self.sourceList = []
		self.f77opts    = ""
		self.f77bound   = True
		self.f77dline   = False
		self.main       = ""
		self.defmain    = True
		self.lfluka     = LINK_PRG[0]

		# Create default run
		self.runs       = [RunInfo(self, DEFAULT_INPUT)]
		self.plots      = []
		self.materials  = {}
		self.geometry   = GeometryInfo()

		self.dicoms     = []
		self.pet        = PetInfo("pet")

		self.var        = {}	# other variables

	# ----------------------------------------------------------------------
	def hash(self):
		h = hash(self.name) ^ hash(self.projFile) ^ hash(self.inputName) \
		  ^ hash(self.exe) ^ hash(self.dir) ^ hash(self.title) \
		  ^ hash(self.notes) ^ hash(self.f77opts) ^ self.f77bound \
		  ^ self.f77dline ^ hash(self.main) ^ hash(self.defmain) \
		  ^ hash(self.calculator)
		for d in self.debugList:
			for i in d:
				h ^= hash(i)
		for f in self.sourceList:
			h ^= hash(f)
		return h
	__hash__ = hash

	# ----------------------------------------------------------------------
	# Change the directory of the project
	# ----------------------------------------------------------------------
	def chdir(self, path):
		try:	# FIXME, I should use absolute paths in the future!
			self.dir = os.path.realpath(os.path.abspath(path))
			os.chdir(self.dir)
		except OSError:
			say("No such directory",path)

		if self.projFile != "":
			fn = os.path.basename(self.projFile)
			self.projFile = os.path.join(self.dir, fn)

		if self.inputFile:
			fn = os.path.basename(self.inputFile)
			fn = os.path.join(self.dir, fn)
			self.inputFile = self.relativePath(fn)

	# ----------------------------------------------------------------------
	# load project
	# ----------------------------------------------------------------------
	def load(self, filename):
		"""Load a project filename and change directory to it.
		All files are referred relative to this directory"""
		self.clear()
		self.version = 0
		self.setModified(False)
		self.setInputModified(False)
		self.projFile = os.path.abspath(filename)
		self.name     = os.path.basename(self.projFile)
		self.dir      = os.path.realpath(os.path.dirname(self.projFile))

		# Change directory
		os.chdir(self.dir)

		# Open project filename
		f = open(self.projFile,"r") #, encoding="utf-8")
		while True:
			(tag, val) = _readVar(f)
			if tag is None: break
			if   tag=="Version":
					try: self.version = int(val)
					except: self.version = 0
			elif tag=="Title":    self.title   = val
			elif tag=="Input":    self.loadInput(val)
			# FIXME move executable to the first RunInfo
			elif tag=="Exec":     self.setExecutable(val)
			elif tag=="Submit":   self.submit = val
			elif tag=="Notes":    self.notes  = _readBlock(f)
			elif tag=="Tabs":     self.tabs   = list(map(float,val.split()))
			elif tag=="Page":     self.page   = val
			elif tag=="Calculator":self.calculator = _readBlock(f)
			elif tag=="Debug":
				vl  = val.split('"')
				# FIXME use of quotes can be problematic
				self.debugList.append([vl[1]] + vl[2].split())
			elif tag=="Main":     self.main      = val
			elif tag=="DefMain":  self.defmain   = (val=="True")
			elif tag=="Source":   self.sourceList.append(val)
			elif tag=="LinkPrg":  self.lfluka    = val
			elif tag=="F77opts":  self.f77opts   = val
			elif tag=="F77bound": self.f77bound  = (val=="True")
			elif tag=="F77dline": self.f77dline  = (val=="True")
			elif tag=="Run":
				if val == DEFAULT_INPUT:
					run = self.runs[0]
				else:
					run = RunInfo(self, val)
					self.runs.append(run)
				run.read(f)
			elif tag=="Material":
				material = self.materials.get(val)
				if material is None:
					self.materials[val] = material = MaterialInfo(val)
				material.read(f)
			elif tag=="Geometry":
				self.geometry.read(f)
			elif tag=="Dicom":
				dicom = DicomInfo(val)
				dicom.read(f)
				dicom.correct()
				if "uid" in dicom.var: # and "directory" in dicom.var:
					self.dicoms.append(dicom)
			elif tag=="Pet":
				self.pet.read(f)
			elif tag=="Plot":
				plot = PlotInfo(val)
				self.plots.append(plot)
				plot.read(f)
				if self.version != __version__:
					plot.correct(self.version)
			else:	# Add in dictionary v.0 correction
				if tag.startswith("material."):
					name = tag[9:]
					self.materials[name] = material = MaterialInfo(name)
				elif val != "":
					self.var[tag] = val
				elif tag != "End":
					self.var[tag] = _readBlock(f)

		# Verify parent/child relations
		for run in self.runs:
			if run.family:
				family = []
				for childname in run.family:
					child = self.getRunByName(childname, run.name)
					if child is not None:
						family.append(childname)
				if len(family) != len(run.family):
					run.family = family

			elif run.parent != "":
				parent = self.getRunByName(run.parent)
				if run.name not in parent.family:
					parent.family.append(run.name)

		f.close()
		return False

	# ----------------------------------------------------------------------
	# save project
	# ----------------------------------------------------------------------
	def save(self, filename=None):
		"""Save a fluka project file. If a filename is specified
		it will be saved under this name.
		WARNING it will also change the current directory"""

		if filename is not None:
			# Change of filename is required
			# Warning also change of directory!!!
			self.projFile = os.path.abspath(filename)
			self.name     = os.path.basename(self.projFile)
			self.dir      = os.path.realpath(os.path.dirname(self.projFile))
			#if olddir != self.dir:
			#	self.inputFile = self.newPath(self.inputFile)
			#	if self.inputFile != "":
			#	inp = os.abspath(self.inputFile)
			#	exe = os.abspath(self.exe)
			#	# Change directory
			os.chdir(self.dir)

			# FIXME Warning change relative filenames
			# input, exe, ...

		if keepBackup:
			backupname = self.name+"~"
			try: os.remove(backupname)
			except: pass
			try: os.rename(self.name, backupname)
			except: pass

		self.setModified(False)

		f = open(self.name,"w")
		f.write("#!%s/flair\n"%(os.path.dirname(os.path.abspath(sys.argv[0]))))
		f.write("# FLUKA Project file\n")
		f.write("Version: %s\n" % (__version__))
#		f.write("Title: %s\n" % (self.title.encode("utf-8")))
		f.write("Title: %s\n" % (self.title))
		Input.utfWrite(f,"Input: %s\n" % (self.inputFile))
		if self.exe != "": Input.utfWrite(f,"Exec: %s\n" % (self.exe))
		if self.submit != "": Input.utfWrite(f,"Submit: %s\n" % (self.submit))
#		if self.notes != "": f.write("Notes:\n%s\f\n" % (self.notes.encode("utf-8")))
		if self.notes != "": f.write("Notes:\n%s\f\n" % (self.notes))
		if self.tabs: f.write("Tabs: %s\n" % (" ".join(map(str,self.tabs))))
		if self.page: f.write("Page: %s\n" % (self.page))
		if self.calculator != "": f.write("Calculator:\n%s\f\n" % (self.calculator.encode("utf-8")))
		for dl in self.debugList:
			Input.utfWrite(f,"Debug: \"%s\" %s %s %s %s %s %s %s %s %s\n" % tuple(dl))
		for fn in self.sourceList:
			Input.utfWrite(f,"Source: %s\n" % (fn))
		if self.main != "": Input.utfWrite(f,"Main: %s\n" % (self.main))
		f.write("DefMain: %s\n"  % (str(self.defmain)))
		if self.f77opts != "": Input.utfWrite(f,"F77opts: %s\n" % (self.f77opts))
		f.write("LinkPrg: %s\n"  % (str(self.lfluka)))
		f.write("F77bound: %s\n" % (str(self.f77bound)))
		f.write("F77dline: %s\n" % (str(self.f77dline)))
		for k in sorted(self.var.keys()):
			v = self.var[k]
			try:
				v = str(v)
			except UnicodeEncodeError:
				v = str(v).encode("utf-8")
			if v != "": Input.utfWrite(f,"%s: %s\n" % (k,v))

		for run in self.runs: run.write(f)
		for material in list(self.materials.values()): material.write(f)
		self.geometry.write(f)
		for dicom in self.dicoms: dicom.write(f)
		self.pet.write(f)
		for plot in self.plots: plot.write(f)
		f.close()

		# Convert to executable
#		if Input._developer:
#			try:
#				mode = os.stat(self.name).st_mode
#				# Add executable flag only when a read flag exists
#				if mode & stat.S_IRUSR: mode |= stat.S_IXUSR
#				if mode & stat.S_IRGRP: mode |= stat.S_IXGRP
#				if mode & stat.S_IROTH: mode |= stat.S_IXOTH
#				os.chmod(self.name,mode)
#			except:
#				pass

	# ----------------------------------------------------------------------
	# Make relative or smallest representation
	# ----------------------------------------------------------------------
	def relativePath(self, filename, run=None):
		"""Convert path to relative to the current directory"""
		if run is None:
			directory = self.dir
		else:
			directory = os.path.abspath(run.getInputDirName())
			if directory == "" and run.parent != "":
				parent = self.findParentRun(run)
				directory = os.path.abspath(parent.getInputDirName())
				if directory == "":
					directory = self.dir

		# XXX from 2.6 could be replaced to relPath
		return relativePath(os.path.abspath(filename), directory)

	# ----------------------------------------------------------------------
	#def _newPath(self, filename):
	#	if len(filename) > 0:
	#		return self.relativePath(os.abspath(filename))
	#	else:
	#		return filename

	# ----------------------------------------------------------------------
	# Load input file
	# ----------------------------------------------------------------------
	def loadInput(self, filename=None):
		"""Load input file"""
		if filename is not None:
			if self.name == "":
				# Change directory
				self.dir = os.path.dirname(os.path.abspath(filename))
				os.chdir(self.dir)
				filename = os.path.basename(filename)
			self.inputFile = self.relativePath(filename)
			self.inputName, ext = os.path.splitext(self.inputFile)
		self.setInputModified(False)
		del self.input
		self.input = Input.Input()
		if self.inputFile != "":
			try:
				self.input.read(self.inputFile)
				self.input.convert2Names()
			except (IOError, OSError):
				say("ERROR:",sys.exc_info()[1])
			if len(self.title)==0:
				try:
					tcard = self.input.cards["TITLE"][0]
					self.title = tcard.extra().strip()
				except:
					pass

	# ----------------------------------------------------------------------
	# Save input file
	# ----------------------------------------------------------------------
	def saveInput(self, filename=None):
		"""Save input file"""
		if self.name == "" and filename is not None:
			# XXX if no name exists switch to the input directory
			# to avoid problems of wrong relative reference
			# of the project file to the input
			self.dir = os.path.dirname(filename)
			os.chdir(self.dir)

		if filename is None:
			filename=self.inputFile
		else:
			self.inputFile = self.relativePath(filename)
			(self.inputName, ext) = os.path.splitext(self.inputFile)
		self.input.write(self.inputFile, keepBackup)
		self.setInputModified(False)

	# ----------------------------------------------------------------------
	# Set exec file
	# ----------------------------------------------------------------------
	def setExecutable(self, file):
		"""Set the executable name"""
		oldexe = self.exe
		if file != "":
			path = os.path.dirname(file)
			if path != self.dir:
				self.exe = file		# Absolute path
			else:
				self.exe = os.path.basename(file)
		else:
			self.exe = file

		# Modify all runs with the same exe
		for run in self.runs[1:]:
			if run.exe == oldexe:
				run.exe = self.exe

	# ----------------------------------------------------------------------
	# executable name of a run
	# ----------------------------------------------------------------------
	def executable(self, run=None):
		if run is not None:
			exe = run.exe
			if run.name == DEFAULT_INPUT:
				if self.exe != "": exe = self.exe
			else:
				if run.parent != "": # Child
					parent = self.findParentRun(run)
					if parent is not None:
						if parent.name == DEFAULT_INPUT: exe = self.exe
						elif parent.exe != "": exe = parent.exe

			if exe:
				return os.path.join(self.dir, exe)

		else:
			if self.exe != "":
				return self.exe

		if flukaExe != "":
			return os.path.join(flukaDir, flukaExe)
		else:
			return os.path.join(flukaDir, DEFAULT_EXE)

	# ----------------------------------------------------------------------
	# Set project modified flag
	# ----------------------------------------------------------------------
	def setModified(self, m=True):
		"""Set the flag that project file is modified"""
		if self.projectModified != m:
			#import traceback; traceback.print_stack()
			self.projectModified = m	# Modified flag

	# ----------------------------------------------------------------------
	def setInputModified(self, m=True):
		"""Set the flag that project file is modified, and keep track
		   of modification time"""
		if self.inputModified != m:
			self.inputModified = m
		self.input.setModified()	# Modification time

	# ----------------------------------------------------------------------
	# Return True if project is modified
	# ----------------------------------------------------------------------
	def isModified(self):
		"""Return TRUE is project is modified"""
		return self.projectModified or self.inputModified

	# ----------------------------------------------------------------------
	def isInputModified(self):
		"""Return TRUE if input is modified"""
		return self.inputModified

	# ----------------------------------------------------------------------
	def inputModifiedTime(self):
		"""Return last modification time"""
		return self.input._modified

	# ----------------------------------------------------------------------
	def _writeGeometry(self, fgeo, title):
		geobegin = Input.Card("GEOBEGIN", ["COMBNAME"], extra=title)
		geobegin.ivopt = 0
		geobegin.ivdbg = 0
		self.input.ivopt = 0
		# Copy important options from real GEOBEGIN
		flugg = False
		try:
			card = self.input.cards["GEOBEGIN"][0]
			geobegin.setWhat(1, card.what(1))
			geobegin.setWhat(2, card.what(2))
			geobegin.setWhat(5, card.what(5))
			geobegin.setWhat(6, card.what(6))
			if card.sdum()=="FLUGG":
				geobegin.setSdum(card.sdum())
				flugg = True
		except (KeyError, TypeError):
			pass

		fgeo.write("* Temporary input file only for geometry checking\n")
		fgeo.write(Input.SCALE)
		glcard = Input.Card("GLOBAL", ["","","","",3.0,1.0])
		try:
			for card in self.input.cards["GLOBAL"]:
				if card.notIgnore():
					for i in range(4):
						glcard.setWhat(i, card.what(i))
		except KeyError:
			pass
		# XXX Lattice card should be saved with FIXED FORMAT, in this case!

		self.input.writeCard(fgeo, glcard, Input.FORMAT_SINGLE)
		self.input.writeCards(fgeo, lambda x:False, Input.FORMAT_FREE)

		# GEOBEGIN More tricky
		geoline = geobegin.toStr(Input.FORMAT_FREE).splitlines()
		fgeo.write(geoline[0].decode()+"\n")
		if flugg: return

		# If VOXELS are present, write the first active one
		if "VOXELS" in self.input.cards:
			for card in self.input.cards["VOXELS"]:
				if card.notIgnore():
					fgeo.write(card.toStr(Input.FORMAT_FREE)+"\n")
					break
		fgeo.write(geoline[1].decode()+"\n")		# Write title
		self.input.writeGeometry(fgeo, Input.FORMAT_FREE, Input.FORMAT_FREE)

	# ----------------------------------------------------------------------
	# Write an input file to debug geometry
	# ----------------------------------------------------------------------
	def writeDebugInput(self, debugRegion):
		"""Write a temporary file for debuging the geometry"""
		global tmpPrefix
		[name,xmin,ymin,zmin,xmax,ymax,zmax,nx,ny,nz] = debugRegion

		# Write a minimized version of the file
		# FIXME preprocessor defines are missing
		# GEOBEGIN
		# ...Bodies
		# ...Regions
		# ...Lattices
		# GEOEND with debug info
		# ROT-DEFI if any
		# STOP

		# FIXME Skip annoying RuntimeWarning
		curdir=os.getcwd()

# tempfile.NamedTemporaryFile(mode='w+b', buffering=- 1, encoding=None, newline=None, suffix=None, prefix=None, dir=None, delete=True, *, errors=None)¶
		fgeo=tempfile.NamedTemporaryFile(mode="w+t",suffix=".inp",dir=curdir,delete=False)
#		tmpname = os.path.basename(os.tempnam("",tmpPrefix))
		tmpname=fgeo.name
		x=tmpname.split(".")
		tmpname=x[-2]
		tmpname = os.path.basename(tmpname)
		# Create the input file
#		fgeo = open(tmpname+".inp", "w")
		self._writeGeometry(fgeo, name)
		geoend   = Input.Card("GEOEND",
			   ["DEBUG", xmax, ymax, zmax, xmin, ymin, zmin, nx, ny, nz],
			   "Geometry debugging",
			   name)
		self.input.writeCard( fgeo, geoend, Input.FORMAT_FREE)
		self.input.writeCards(fgeo,
				lambda x:x.tag=="ROT-DEFI",
				Input.FORMAT_FREE)
		fgeo.write("STOP\n")
		fgeo.close()

		return tmpname

	# ----------------------------------------------------------------------
	# Write an input file to be used with plotgeom
	# ----------------------------------------------------------------------
	def writePlotgeomInput(self, plotgeom):
		"""Write a temporary file for plotting the geometry"""
		global tmpPrefix

		# Write a minimized version of the file
		# FIXME preprocessor defines are missing
		# GEOBEGIN
		# ...Bodies
		# ...Regions
		# ...Lattices
		# GEOEND
		# ROT-DEFI if any
		# MATERIAL
		# ASSIGNMA
		# PLOTGEOM card
		# STOP
		curdir=os.getcwd()
#tempfile.mkstemp(suffix=None, prefix=None, dir=None, text=False)		
		fgeo=tempfile.NamedTemporaryFile(mode="w+t",suffix=".inp",dir=curdir,delete=False)
#		tmpname = os.path.basename(os.tempnam("",tmpPrefix))
		tmpname=fgeo.name
		x=tmpname.split(".")
		tmpname=x[-2]
		tmpname = os.path.basename(tmpname)
		# Create the input file
#		fgeo = open(tmpname+".inp", "w")
		self._writeGeometry(fgeo, plotgeom.extra())
		geoend = Input.Card("GEOEND")
		self.input.writeCard(fgeo, geoend, Input.FORMAT_FREE)
		self.input.writeCards(fgeo,
				lambda x:x.tag in ("ROT-DEFI", "MATERIAL", "ASSIGNMA",
						   "USRICALL", "USRGCALL", "MGNFIELD"),
				Input.FORMAT_FREE)
		self.input.writeCard( fgeo, plotgeom, Input.FORMAT_FREE)
		fgeo.write("STOP\n")
		fgeo.close()
		return tmpname

	# ----------------------------------------------------------------------
	# Run Input file
	#
	# FIXME move to RunInfo
	# ----------------------------------------------------------------------
	def writeInput(self, run=None, cwd=None, defines=None):
		"""Write the input file with the info in Run"""

		# If default input
		if run is None or run.name == DEFAULT_INPUT:
			# if we are writing the default input in the default directory
			# ---> use the standard writing routine
			if (cwd is None or cwd == self.dir) and defines is None:
				if self.isInputModified(): self.saveInput()
				return
			else:
				# else create a clone of the input
				defaultDefines = False
				run = self.runs[0]

		else:
			# Create a clone and override parameters if necessary
			defaultDefines = False
			if run.parent != "":
				parent = self.findParentRun(run)
				if parent is not None and parent.name == DEFAULT_INPUT:
					defaultDefines = True

		# Special Run override in needed
		runInput = self.input.clone()

		if run.title != "":
			# Modify title cards
			for card in runInput["TITLE"]:
				card.setExtra(run.title)

		if run.rnd != 0:
			# Modify title randomize cards
			for card in runInput["RANDOMIZ"]:
				card.setWhat(2, run.rnd)

		if run.primaries != 0 or run.time != 0:
			# Modify title start cards
			for card in runInput["START"]:
				if run.primaries != 0:
					card.setWhat(1, run.primaries)
				if run.time  != 0:
					card.setWhat(6, run.time)

		if defines is None: defines = run.defines

		# if there are defines and we unselect all
		# then the defines will be empty and don't work!!
		if not defaultDefines: #and defines:
			# Modify #define cards
			for card in runInput["#define"]:
				if card._indent != 0: continue
				# Search if it exists on defines
				name = card.sdum()
				for n,v in defines:
					if n == name:
						card.setEnable(True)
						if v != "": card.setWhat(1,v)
						break
				else:
					card.setEnable(False)

		if cwd is not None:
			inpfile = "%s/%s.inp"%(cwd, run.getInputBaseName())
		else:
			inpfile = "%s.inp"%(run.name)
		errors = runInput.preprocess()
		runInput.write(inpfile)
		del runInput
		return errors

	# ----------------------------------------------------------------------
	# return a tuple to spawn a fluka command
	# FIXME only used by Plot, should be canceled
	# ----------------------------------------------------------------------
	def runCmd(self, inpname, last, prev=-1):
		"""Return a tuple containing the command to spawn a fluka run"""
		cmd = []
		cmd.append( command("rfluka") )
		if self.exe != "":
			cmd.append("-e")
			cmd.append(self.exe)
		if prev > 0:
			cmd.append("-N")
			cmd.append(str(prev))
		cmd.append("-M")
		cmd.append(str(last))
		if inpname:
			cmd.append(inpname)
		else:
			cmd.append(self.inputName)
#		print ("command to run " , cmd)
		return tuple(cmd)

	# ----------------------------------------------------------------------
	# return a tuple with the compile command to be spawned
	# ----------------------------------------------------------------------
	def compileCmd(self, filename):
		"""Return a tuple containing the compile command to be spawned"""
		(fn, ext) = os.path.splitext(filename)

		if ext == ".c":
			return ( ccprogram, "-m32", "-g", "-c", "-I"+os.getenv(flukaVar), filename )

		elif ext==".cc" or ext==".cpp":
			return ( cppprogram, "-m32", "-g", "-c", "-I"+os.getenv(flukaVar), filename )

		else:
			c = [ command("fff") ]
			if self.f77bound:
				c.append("-b")
			else:
				c.append("-n")
			if self.f77dline:
				c.append("-D")
			else:
				c.append("-N")
			if self.f77opts != "":
				for i in self.f77opts.split():
					c.append("-o")
					c.append(i)
			c.append(filename)
			return tuple(c)

	# ----------------------------------------------------------------------
	# Link command
	# ----------------------------------------------------------------------
	def linkCmd(self):
		"""Return a tuple containing the link command to be spawned"""
		if checkDir(): return
		c = [ os.path.join(flukaDir, "flutil", self.lfluka),
			"-o", self.exe]

		if self.defmain:
			# Default main
			c.extend(["-m", "fluka"])
		else:
			if self.main != "":
				# User main function if any
				c.extend(["-m", self.main])

		dirs = []
		libs = []
		for f in self.sourceList:
			(fn, ext) = os.path.splitext(f)
			if ext in (".a", ".so"):
				dir = os.path.dirname(fn)
				if dir=="": dir = "."
				if dir not in dirs:
					dirs.append(dir)
				lib = os.path.basename(fn)
				if lib[0:3] == "lib":
					libs.append(lib[3:])
				else:
					libs.append(f)
			else:
				c.append(fn+".o")

		if len(dirs)>0:
			c.append("-L")
			c.append(flukaDir)

			for d in dirs:
				c.append("-L")
				c.append(d)

			for l in libs:
				c.append("-l")	# Before
				c.append(l)
				c.append("-a")	# and after
				c.append(l)

		return tuple(c)

	# ----------------------------------------------------------------------
	# write a makefile to build the executable
	# ----------------------------------------------------------------------
	def writeMakefile(self, filename=None):
		"""Write a makefile to build the executable"""
		if filename is None:
			filename = "Makefile"

		f = open(filename, "w")

		Input.utfWrite(f,"# Exported Makefile from project: %s\n" % (self.name))
		Input.utfWrite(f,"# Title: %s\n" % (self.title))

		f.write("\nF77  = %s\n" % (command("fff")))
		f.write("LINK = $(%s)/flutil/%s\n\n" % \
				(flukaVar, self.lfluka))
		Input.utfWrite(f,"EXE  = %s\n" % (self.exe))

		srcs = []
		dirs = []
		libs = []
		for file in self.sourceList:
			(fn, ext) = os.path.splitext(file)
			if ext == ".a":
				dir = os.path.dirname(fn)
				if dir!="": dirs.append(dir)
				lib = os.path.basename(fn)
				if lib[0:3] == "lib":
					libs.append(lib[3:])
				else:
					libs.append(file)
			else:
				srcs.append(file)

		Input.utfWrite(f,"SRCS = %s\n" % (" ".join(srcs)))
		Input.utfWrite(f,"OBJS = $(SRCS:.f=.o)\n\n")
		comp = self.compileCmd("")
		Input.utfWrite(f,"%%.o: %%.f\n\t$(F77) %s $<\n\n" % (" ".join(comp[1:-1])))
		f.write("$(EXE): $(OBJS)\n")

		dirstr = " ".join([n for n in dirs if "-L "+n])
		libstr = " ".join([n for n in libs if "-l %s -a %s" % (n,n)])

		if self.defmain:
			Input.utfWrite(f,"\t$(LINK) -o $@ -m fluka %s %s $^\n" \
				% (dirstr, libstr))
		else:
			if self.main!="":
				Input.utfWrite(f,"\t$(LINK) -o $@ -m %s %s %s $^\n" \
					% (self.main, dirstr, libstr))
			else:
				Input.utfWrite(f,"\t$(LINK) -o $@ %s %s $^\n" \
					% (dirstr, libstr))

		f.write("\nclean:\n")
		f.write("\trm -f $(EXE) $(EXE).map $(OBJS)\n")
		f.close()

	# ----------------------------------------------------------------------
	# Delete files generated by running fluka
	# [ran]name[0-9][0-9][0-9]*
	# ----------------------------------------------------------------------
	def deleteFiles(self, name, pid=None):
		"""Delete temporary files and directory"""
		p  = re.compile(r"^%s\d\d\d[._].+$" % (name))
		rp = re.compile(r"^ran%s\d\d\d$" % (name))
		for fn in Utils.listdir(self.dir):
			if p.match(fn) or rp.match(fn):
				os.remove(fn)

		if pid is not None:
			rundir = "fluka_%d" % (pid)
			try:
				os.stat(rundir)
				for f in Utils.listdir(rundir):
					os.remove("%s/%s"%(rundir,f))
				os.rmdir(rundir)
			except OSError:
				pass

	# ----------------------------------------------------------------------
	def _scanUSR(self, usr, score, tag):
		for card in self.input[tag]:
			if card.ignore(): continue
			if tag=="DETECT":
				unit = -17
			else:
				try: unit = int(card.units(False)[0])
				except: continue

			# Use only binary files
			if unit>0: continue
			unit = abs(unit)
			try:
				i = usr[unit]
			except:
				usr[unit] = score

	# ----------------------------------------------------------------------
	# Scan Input for available binary scoring cards
	# ----------------------------------------------------------------------
	def scanInput4Usr(self):
		"""Scan Input for available binary scoring cards"""
		usr = {}
		self._scanUSR(usr, 'r', "RESNUCLE")
		self._scanUSR(usr, 'x', "USRBDX"  )
		self._scanUSR(usr, 'b', "USRBIN"  )
		self._scanUSR(usr, 'c', "USRCOLL" )
		self._scanUSR(usr, 't', "USRTRACK")
		self._scanUSR(usr, 'y', "USRYIELD")
		self._scanUSR(usr, 'd', "DETECT"  )
		return usr

	# ----------------------------------------------------------------------
	# For each run fill default UsrInfo
	# ----------------------------------------------------------------------
	def findRunUsrInfo(self, run):
		if run.name == DEFAULT_INPUT:
			defines = None
		else:
			defines = run.defines
		self.input.preprocess(defines)
		usr = self.scanInput4Usr()
		# add all non existing units in run
		for u,t in list(usr.items()):
			for i in run.usrinfo:
				if i.unit == u:
					break
			else:
				ui = UsrInfo(run)
				ui.set(u, t)
				run.usrinfo.append(ui)
		run.usrinfo.sort(key=lambda x: x.unit )

	# ----------------------------------------------------------------------
	# return run by name under a certain parent if not None
	# ----------------------------------------------------------------------
	def getRunByName(self, name, parentName=None):
		if name and name[0]=="<": return self.runs[0]	# DEFAULT_INPUT
		for run in self.runs:
			if run.name == name:
				if parentName is not None:
					if parentName == run.parent: return run
				else:
					return run
		return None

	# ----------------------------------------------------------------------
	def getRunIndex(self, run):
		return self.runs.index(run)

	# ----------------------------------------------------------------------
	def findUniqueRunName(self, name, parent=""):
		while True:
			for run in self.runs:
				if run.name==name and run.parent==parent:
					pat = _LAST_NUM.match(name)
					if pat:
						prefix = pat.group(1)
						n = int(pat.group(2))+1
						name = prefix+str(n)
					else:
						name += "-1"
					break
			else:
				return name

	# ----------------------------------------------------------------------
	# Rename a run and its childs
	# ----------------------------------------------------------------------
	def renameRun(self, run, name):
		# rename all children with the same way
		oldname = run.name
		run.name = name
		if run.family:
			newfamily = []
			basename = os.path.basename(name)
			for i,childname in enumerate(run.family):
				newname = self.getSpawnRunName(run, i)
				#newname = "_%s_%s"%(basename, childname[-2:])
				child = self.getRunByName(childname, oldname)
				if child is not None:
					newfamily.append(newname)
					child.name = newname
					child.parent = name
			run.family = newfamily

	# ----------------------------------------------------------------------
	# expand escape characters in string
	# ----------------------------------------------------------------------
	def getSpawnRunName(self, run, idx, esc=False):
		name = run.spawnName or _DEFAULT_RUNNAME
		if esc:
			name = name.replace(r"\F", re.escape(self.name))
			name = name.replace(r"\I", re.escape(run.getInputBaseName()))
		else:
			name = name.replace(r"\F", self.name)
			name = name.replace(r"\I", run.getInputBaseName())

		name = name.replace(r"\0n", "%d"%(idx-1))
		name = name.replace(r"\02", "%02d"%(idx-1))
		name = name.replace(r"\03", "%03d"%(idx-1))
		name = name.replace(r"\04", "%04d"%(idx-1))
		name = name.replace(r"\n", "%d"%(idx))
		name = name.replace(r"\2", "%02d"%(idx))
		name = name.replace(r"\3", "%03d"%(idx))
		name = name.replace(r"\4", "%04d"%(idx))

		if name.find(r"\a")>=0 or name.find(r"\A")>=0:
			hi,lo = divmod(idx-1,26)
			name = name.replace(r"\a", chr(hi+ord('a')) + chr(lo+ord('a')))
			name = name.replace(r"\A", chr(hi+ord('A')) + chr(lo+ord('A')))
		name = name.replace(r"\x", "%02X"%(idx))
		if name.find(r"\X")>=0:
			name = name.replace(r"\X", bmath.int2roman(idx))
		return name

	# ----------------------------------------------------------------------
	def findUniquePlotName(self, plot):
		# Find a unique file name
		if not plot.name[-1].isdigit():
			name = plot.name
			start = 1
		else:
			m = _LAST_NUM.match(plot.name)
			name = m.group(1)
			start = int(m.group(2))+1

		for i in range(start,100):
			namei = "%s%02d"%(name,i)
			for pl in self.plots:
				if pl.name == namei:
					break
			else:
				plot.name = namei
				break

	# ----------------------------------------------------------------------
	def setRunRandomSeed(self, run, parent=None, sequential=False):
		if sequential:
			# Get father and add child index
			if parent is None or parent.rnd == 0:
				rnd = 0
				# Scan in input to find the randomize card
				for card in self.input["RANDOMIZ"]:
					if card.ignore(): continue
					rnd = card.intWhat(2) + run.index
					break
			else:
				rnd = parent.rnd + run.index
			if rnd == 0:
				rnd = run.index
		else:
			# Random random seed
			for i in range(1000):
				rnd = int(random.randint(3, 0x7fffffff))
				found = False
				for r in self.runs:
					if rnd == r.rnd:
						found = True
						break
				if not found: break

		run.rnd = rnd

	# ----------------------------------------------------------------------
	# Find parent of a child run
	# ----------------------------------------------------------------------
	def findParentRun(self, run):
		if run.parent == "": return None
		for parent in self.runs:
			if run.parent == parent.name:
			#if run.name in parent.family:
				return parent
		return None

	# ----------------------------------------------------------------------
	# Delete run
	# ----------------------------------------------------------------------
	def delRun(self, run):
		# Check for child
		if run.parent != "":
			# find parent
			for parent in self.runs:
				found = False
				if run.parent == parent.name:
				#if run.name in parent.family:
					parent.family.remove(run.name)
					found = True
					break
				if found: break
		# Check for family
		elif run.family:
			family = run.family	# save family
			run.family = None	# set it to none
			for childname in family:
				child = self.getRunByName(childname, run.name)
				# to make it faster cut the link with the parent
				child.parent = ""
				self.delRun(child)

		self.runs.remove(run)

	# ----------------------------------------------------------------------
	# Return plot by name
	# ----------------------------------------------------------------------
	def getPlot(self, name, t=None):
		for plot in self.plots:
			if name == plot.name:
				if t is None:
					return plot
				elif plot.type == t:
					return plot
		return None

	# ----------------------------------------------------------------------
	# @return local or global material
	# If it doesn't exist create it in local
	# ----------------------------------------------------------------------
	def getMaterial(self, name):
		mat = self.materials.get(name)
		if mat is None:
			mat = MaterialInfo.globalMaterials.get(name)
			if mat is None:
				self.materials[name] = mat = MaterialInfo(name)
		return mat

	# ----------------------------------------------------------------------
	# Move a local material to global
	# ----------------------------------------------------------------------
	def material2Global(self, name):
		mat = self.materials.get(name)
		if mat is None:
			MaterialInfo.globalMaterials[name] = name
			del self.materials[name]

	# ----------------------------------------------------------------------
	# Move a local material to global
	# ----------------------------------------------------------------------
	def material2Local(self, name):
		mat = MaterialInfo.globalMaterials.get(name)
		if mat:
			self.materials[name] = mat
			del MaterialInfo.globalMaterials[name]

	# ----------------------------------------------------------------------
	# Prepare the gnuplot palette
	# ----------------------------------------------------------------------
	def gnuplotPalette(self, materialList):
		i = 0
		palette = []
		for name in materialList:
			i += 1
			palette.append("%d \"%s\""%(i,self.getMaterial(name)["color"]))
		return palette

	# ----------------------------------------------------------------------
	# return dicom info structure from UID
	# ----------------------------------------------------------------------
	def dicomInfo(self, uid):
		for di in self.dicoms:
			if uid == di["uid"]:
				return di
		return None

#===============================================================================
def init():
	initSystemLayers()
	setFlukaDir()

if __first:
	__first = False
	init()

#===============================================================================
if __name__ == "__main__":
	#path = "/one/two/three/four/five_six_seven_eight"
	#for i in range(10,len(path)):
	#	say(i)
	#	say(abbrevPath(path, i))
	#sys.exit(0)
	#Input.init()
	p = Project()
	p.load(sys.argv[1])
	for run in p.runs:
		p.findRunUsrInfo(run)
		for v in run.usrinfo:
			say(v)
			items = v.fileFilters()
			say(items)
	p.writeMakefile("test.mak")
