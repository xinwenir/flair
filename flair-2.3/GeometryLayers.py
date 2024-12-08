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
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR 'AS IS'
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

__author__ = 'Vasilis Vlachoudis'
__email__  = 'Vasilis.Vlachoudis@cern.ch'

import io
import sys
import os
import re
import math
import string
import subprocess
from log import say

try:
	from io import StringIO
except ImportError:
	from io import StringIO
try:
	import pickle as pickle
except ImportError:
	import pickle

try:
	from tkinter import *
	import tkinter.messagebox as messagebox
	from tkinter.colorchooser import askcolor
except ImportError:
	from tkinter import *
	import tkinter.messagebox as messagebox
	from tkinter.colorchooser import askcolor

import undo
import tkFlair
import tkExtra

import Data
import Input
import Manual
import Palette
import Project
import bFileDialog
import CalibrateImage
from Constants import *	# V2

_LAYER_CLIP = tkFlair._FLAIRCLIP + "<layer>"

_ON  = '[X]'
_OFF = '[  ]'

_MAXUSRBIN   =  10
_MAXPALETTE  =   5
_APPLY_DELAY = 500

_ADD = '<add>'
_ORDER = [
	 'Options',
	 'Show',
	 '3D',
	 'Beam',
	 'Image',
	 'Palette',
	 'Userdump',
	 'Usrbin',
	 'Voxel',
	 ]

#	 'Options',
#	 'Show',
#	 'Image',
#	 'Beam',
#	 'Userdump',
#	 'Usrbin',
#	 '3D',
#	 'Voxel',
#	 'Palette' ]

_PROJECTIONS = ["Orthographic", "Perspective", "Combo"]

DEFLAYERS = {
	'_system'           : False,

	# Options
	'Options'           : True,
	'Options._visible'  : True,
	'Options.axes'      : True,
	'Options.title'     : True,
	'Options.errors'    : True,
	'Options.shortcut'  : '',
	'Options.grid'      : True,
	'Options.gridlevel' : 25,
	'Options.latticelevel': 25,
	'Options.crosshair' : False,
	'Options.viewport'  : True,
	'Options.font'      : 'fixed8x13',
	'Options.gridfont'  : 'fixed8x13',
	'Options.palettefont': 'fixed8x13',
	'Options.textbackground': 200,

	# Show
	'Show'              : True,
	'Show._visible'     : True,
	'Show.vertex'       : True,
	'Show.label'        : 'Region',
	'Show.color'        : 'Material',
	'Show.lattice'      : False,
	'Show.voxel'        : False,
	'Show.palette'      : 'Palette',

	# Background image
	'Image'             : False,
	'Image._visible'    : False,
	'Image.file'        : '',
	'Image.alpha'       : 127,
	'Image.black'       : 0,
	'Image.white'       : 0xffffff,
	'Image.prompt'      : False,

	# Transformations
	'Transformations'   : False,
	'Transformations._visible' : False,

	# Beam
	'Beam'              : False,
	'Beam._visible'     : False,

	# Userdump
	'Userdump'          : False,
	'Userdump._visible' : False,
	'Userdump.file'     : '',
	'Userdump.start'    : 0,
	'Userdump.n'        : 100,

	# Usrbin
	'Usrbin'            : False,
	'Usrbin._visible'   : False,
	'Usrbin.input'      : False,
	'Usrbin.file'       : '',
	'Usrbin.det'        : 1,
	'Usrbin.norm'       : 1.0,
	'Usrbin.xofs'       : '',
	'Usrbin.yofs'       : '',
	'Usrbin.zofs'       : '',
	'Usrbin.rotdefi'    : '',
	'Usrbin.palette'    : 'Palette',
	'Usrbin.alpha'      : 0,

	# 3D
	'3D'                : False,
	'3D._visible'       : True,
	'3D.ambient'        : 64,
	'3D.antialias'      : 1,
#	'3D.focal'          : 100.0,
	'3D.fov'            : 60,
	'3D.clip.0'         : '',
	'3D.clip.1'         : '',
	'3D.clip.2'         : '',
	'3D.clipnegative.0' : False,
	'3D.clipnegative.1' : False,
	'3D.clipnegative.2' : False,
	'3D.deflights'      : False,
	'3D.edgedetect'     : False,
	'3D.project.0'      : '',
	'3D.project.1'      : '',
	'3D.project.2'      : '',
	'3D.projection'     : 0,
	'3D.reflections'    : 0,
	'3D.shadows'        : True,
	'3D.skip1stblack'   : False,
	'3D.usrbinastexture': True,
	'3D.xray'           : 0,

	# Color band
	'Palette'         : False,
	'Palette._visible': False,
	'Palette.palette' : 'FLUKA',
	'Palette.log'     : True,
	'Palette.min'     : 1.0E-5,
	'Palette.max'     : 1.0E5,
	'Palette.smooth'  : False,
	'Palette.n'       : 30,
	'Palette.cpd'     : 3,
	'Palette.inv'     : False,
	'Palette.alphamin': False,	# Transparent below minimum limit
	'Palette.alphamax': False,	# Transparent above maximum limit
	'Palette.label'   : '',

	# Voxel, ROI - regions of interest in voxel
	'Voxel'           : False,
	'Voxel._visible'  : False,
	'Voxel.roialpha'  : 200,
	'Voxel.roi'       : '',
}

# Add particle lists in userdump
for particle in Input.Particle.list:
	DEFLAYERS["Userdump.%s"%(particle)] = 0

# ------------------------- Colors ---------------------------------------------
COLORVALUES = [ 'None',			# 0
		'Region',		# 1
		'Material',		# 2
		'Material-Decay',	# 3
		'Corrfactor',
		'Corrfactor-dE/dx',
		'Corrfactor-rho',
		'Deltaray',
		'Density',
		'Field',
		'e-Production',
		'e-Transport',
		'g-Production',
		'g-Transport',
		'Importance',
		'Importance-E',
		'Importance-H',
		'Importance-N',
		'Splitting',
		'Splitting-E',
		'Splitting-H',
		'Splitting-N' ]

# Labels
LABELS = [ 'None', 'Region', 'Material', 'Value' ]
LABEL_NONE     = 0
LABEL_REGION   = 1
LABEL_MATERIAL = 2
LABEL_VALUE    = 3

SHORTCUTS    = [ "",
		 "F3",
		 "F4",
		 "F5",
		 "F6",
		 "F7",
		 "F8",
		 "F9",
		 "F10",
		 "F11",
		 "F12",
#		 "Control-Key-F1",
#		 "Control-Key-F2",
#		 "Control-Key-F3",
#		 "Control-Key-F4",
#		 "Control-Key-F5",
#		 "Control-Key-F6",
#		 "Control-Key-F7",
#		 "Control-Key-F8",
#		 "Control-Key-F9",
#		 "Control-Key-F10",
#		 "Control-Key-F11",
#		 "Control-Key-F12",
#		 "Shift-Key-F1",
#		 "Shift-Key-F2",
#		 "Shift-Key-F3",
#		 "Shift-Key-F4",
#		 "Shift-Key-F5",
#		 "Shift-Key-F6",
#		 "Shift-Key-F7",
#		 "Shift-Key-F8",
#		 "Shift-Key-F9",
#		 "Shift-Key-F10",
#		 "Shift-Key-F11",
#		 "Shift-Key-F12"
		]

USRBINNAMES  = [ 'Usrbin' ]
PALETTENAMES = [ 'Palette' ]

_usrbinPat = re.compile(r'(\d+) (.+) \[(\d+)\]')

# Execute for first time
# Insert Additional USRBINS
pos  = _ORDER.index('Usrbin')+1
keys = [ x for x in DEFLAYERS if x.startswith('Usrbin') ]
for i in range(_MAXUSRBIN,1,-1):
	_ORDER.insert(pos, 'Usrbin%d'%(i))
	for k in keys:
		newkey = 'Usrbin%d%s'%(i,k[6:])
		DEFLAYERS[newkey] = DEFLAYERS[k]
for i in range(2,_MAXUSRBIN+1):
	USRBINNAMES.append('Usrbin%d'%(i))

# Insert Additional Palette
pos  = _ORDER.index('Palette')+1
keys = [ x for x in DEFLAYERS if x.startswith('Palette') ]
for i in range(_MAXPALETTE,1,-1):
	_ORDER.insert(pos, 'Palette%d'%(i))
	for k in keys:
		newkey = 'Palette%d%s'%(i,k[7:])
		DEFLAYERS[newkey] = DEFLAYERS[k]
for i in range(2,_MAXPALETTE+1):
	PALETTENAMES.append('Palette%d'%(i))
del keys, i, pos

## ------------------------------------------------------------------------------
def option(dictionary, key, default=None):
	try:
		return dictionary[key]
	except:
		if default is None:
			return dictionary.get(key, DEFLAYERS.get(key))
		else:
			return dictionary.get(key, default)

# ------------------------------------------------------------------------------
def getMarkers(layers):
	markers = []
	n = int(layers.get('Image.marker',0))
	for i in range(n):
		key = 'Image.marker.%d'%(i)
		mark = layers.get(key)
		if mark is None: return
		markers.append(list(map(float,mark.split())))
	return markers

#===============================================================================
class Combobox(tkExtra.Combobox):
	def setGeometryLayersDialog(self, dlg):
		self._geoDlg = dlg

	# --------------------------------------------------------------------
	def beforeShow(self):
		self._geoDlg.get()
		self._geoDlg.fillLayers()

#===============================================================================
# Geometry Layers Frame
#===============================================================================
class GeometryLayersFrame(Frame):
	_PAGES = None

	def __init__(self, master, flair):
		Frame.__init__(self, master)
		self.flair        = flair
		self.layer        = None
		self.activepage   = None
		self._apply_after = None
		self._init        = True	# Avoid updating

		# --- Layer ---
		frame = Frame(self)
		frame.pack(side=TOP, fill=X)
		self.layerBox = Combobox(frame, False, width=8,
				command=self.changeLayer)
		self.layerBox.setGeometryLayersDialog(self)
		self.layerBox._text.config(background='White')
		self.layerBox._text.bind('<Return>',   self.nothing)
		self.layerBox._text.bind('<KP_Enter>', self.nothing)
		self.layerBox.pack(side=LEFT, expand=YES, fill=X)
		tkExtra.Balloon.set(self.layerBox, 'Layer name currently being edited')

		# buttons
		b = Button(frame, image=tkFlair.icons['clone'], command=self.cloneLayer,
				pady=0, padx=3)
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b, 'Clone active layer')

		b = Button(frame, image=tkFlair.icons['add'], command=self.addLayer,
				pady=0, padx=3)
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b, 'Add layer to project')

		b = Button(frame, image=tkFlair.icons['del'], command=self.delLayer,
				pady=0, padx=3)
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b, 'Delete layer from project')

		self.globalLayer = BooleanVar(self)
		self.globalButton = Checkbutton(frame, text='Global',
				variable=self.globalLayer,
				onvalue=1, offvalue=0,
				anchor=W)
		tkExtra.Balloon.set(self.globalButton,
			'If selected the layer will be saved in the flair.ini, ' \
			'instead of the project and will be available to all projects')
		self.globalButton.pack(side=RIGHT)

		# --- Action buttons ---
		frame = Frame(self)
		frame.pack(side=BOTTOM, fill=X)
		Button(frame, text='Help', image=tkFlair.icons['info'],
			compound=LEFT, takefocus=FALSE,
			command=self.help).pack(side=LEFT, padx=5, pady=5)
		Button(frame, text='Apply', image=tkFlair.icons['ok'],
			compound=LEFT, takefocus=FALSE,
			command=self.apply).pack(side=RIGHT, padx=5, pady=5)

		# --- Center frame ---
		splitter = tkExtra.VSplitter(self, 0.35)
		splitter.pack(side=TOP, expand=YES, fill=BOTH)

		# --- Selection listbox ---
		self.listbox = tkExtra.MultiListbox(splitter.topFrame(),
					(("Value", 16, None),
					 ("Sel",    3, None)),
					header=False,
					selectmode=BROWSE,
					height=5,
					stretch="first",
					background="White")
		self.listbox.pack(side=TOP, expand=YES, fill=BOTH, padx=2, pady=2)
		self.listbox.bind('<<ListboxSelect>>',		self.listboxChanged)
		self.listbox.bindList('<space>',		self.selectToggle)
		self.listbox.bindList('<Return>',		self.add2Listbox)
		self.listbox.bindList('<KP_Enter>',		self.add2Listbox)
		self.listbox.bindList('<Insert>',		self.add2Listbox)
		self.listbox.bindList('<Delete>',		self.removeListbox)
		self.listbox.lists[0].bind("<ButtonRelease-1>",	self.listbox0Release)
		self.listbox.lists[1].bind("<ButtonRelease-1>",	self.listbox1Release)
		self.listbox.bindList('<<Copy>>',		self.copy)
		self.listbox.bindList('<<Cut>>',		self.cut)
		self.listbox.bindList('<<Paste>>',		self.paste)
		self.listbox.setPopupMenu([
				('Add',    0, self.add2Listbox,   tkFlair.icons["add"]),
				('Clone',  0, self.cloneListbox,  tkFlair.icons["clone"]),
				('Remove', 0, self.removeListbox, tkFlair.icons["x"]) ])

		# ---- Bottom frame ----
		# Frame title
		self.title = Label(splitter.bottomFrame(),
				foreground=tkFlair._TITLE_FOREGROUND_COLOR,
				background=tkFlair._TITLE_BACKACTIVE_COLOR,
				relief=RAISED)
		#self.title.pack(side=TOP, fill=X, padx=2, pady=2)
		self.title.grid(row=0, column=0, columnspan=2, sticky=EW)

		# main frame
		self._scrollFrame = tkExtra.ScrollFrame(splitter.bottomFrame(), borderwidth=2, relief=SUNKEN)
		self._scrollFrame.grid(row=1, column=0, sticky=NSEW, padx=2, pady=2)
		#self._scrollFrame.pack(side=LEFT, expand=YES, fill=BOTH, padx=2, pady=2)
		self.mainframe = self._scrollFrame()

		#sbx=tkExtra.AutoScrollbar(splitter.bottomFrame(), orient=HORIZONTAL,
		#		command=self._scrollFrame.xview)
		#sbx.grid(row=2,column=0,sticky=NSEW)
		sby=tkExtra.AutoScrollbar(splitter.bottomFrame(), orient=VERTICAL,
				command=self._scrollFrame.yview)
		sby.grid(row=1, column=1, sticky=NS)
		#self._scrollFrame["xscrollcommand"]=sbx.set
		self._scrollFrame["yscrollcommand"]=sby.set
		splitter.bottomFrame().grid_columnconfigure(0, weight=1)
		splitter.bottomFrame().grid_rowconfigure(1, weight=1)

		# Bindings
		self.bind("<F1>",       self.help)
		self.bind("<Return>",   self.apply)
		self.bind("<KP_Enter>", self.apply)

		# Initialize pages, but only after the class is defined!
		if GeometryLayersFrame._PAGES is None:
			GeometryLayersFrame._PAGES = {
			   'Options'   : GeometryLayersFrame.optionsFrame,
			   'Show'      : GeometryLayersFrame.showFrame,
			   'Image'     : GeometryLayersFrame.imageFrame,
			   'Beam'      : GeometryLayersFrame.emptyFrame,
			   'Userdump'  : GeometryLayersFrame.userdumpFrame,
			   'Usrbin'    : GeometryLayersFrame.usrbinFrame,
			   '3D'        : GeometryLayersFrame.d3Frame,
			   'Voxel'     : GeometryLayersFrame.voxelFrame,
			   'Palette'   : GeometryLayersFrame.paletteFrame,
			}
		# Populate fonts
		self.populateFonts()

		# Populate pages
		self.pages = {}

		# default values
		for page, method in list(GeometryLayersFrame._PAGES.items()):
			self.pages[page] = method(self)
		self._scrollFrame.defaultBinds()

		# copy pages for Usrbin# and Colorband#
		page = self.pages['Usrbin']
		for i in range(2,_MAXUSRBIN+1):
			self.pages['Usrbin%d'%(i)] = page
		page = self.pages['Palette']
		for i in range(2,_MAXPALETTE+1):
			self.pages['Palette%d'%(i)] = page
		self._init = False	# Avoid updating

	# --------------------------------------------------------------------
	def show(self, project, layer=None):
		self.project = project

		# set the list of layers
		layers = self.project.geometry.layersList()
		if layer is not None or not self.layerBox.get():
			if not layer in layers: layer = layers[0]
			self.layerBox.fill(layers)
			self.layerBox.set(layer)

		# Select last page
		self.listbox.selection_clear(0,END)
		oldPage = self.activepage
		self.activepage = None
		self.showPage(_ORDER[0])	# Show options by default
		# do not show immediately Usrbin/Colorband otherwise it
		# will save empty options
		if oldPage is not None:
			self.listbox.selection_clear(0,END)
			self.showPage(oldPage)

	# --------------------------------------------------------------------
	def refresh(self):
		self.d3Refresh()
		self.voxelRefresh()

	# --------------------------------------------------------------------
	# Apply button
	# --------------------------------------------------------------------
	def apply(self, event=None):
		if self._init: return
		if self._apply_after is not None:
			self.after_cancel(self._apply_after)
			self._apply_after = None
		oldHash = self.project.geometry.hash()
		self.get()
		self.project.geometry.correctGlobalLayers()
		tkFlair.updateLayers()
		if self.project.geometry.hash() != oldHash:
			# FIXME unfortunately data= is one of the missing parameters in Tkinter!!
			# and it is not copied in the event class
			self.event_generate('<<LayersUpdate>>') #,data=self.layer.name)

			if self.projection.get() == _PROJECTIONS[0]:
				self.fov["state"] = DISABLED
				self._fovLabel["state"] = DISABLED
			else:
				self.fov["state"] = NORMAL
				self._fovLabel["state"] = NORMAL

# For debugging multiple applies
#		import traceback
#		say('\n====================================================')
#		traceback.print_stack()

	# --------------------------------------------------------------------
	def applyLate(self, event=None):
		if self._apply_after is not None: self.after_cancel(self._apply_after)
		self._apply_after = self.after(_APPLY_DELAY, self.apply)

	# --------------------------------------------------------------------
	def nothing(self, event=None):
		return 'break'

	# --------------------------------------------------------------------
	def entry(self, key):
		return '%s.%s' % (self.activepage, key)

	# ----------------------------------------------------------------------
	def copy(self, event=None):
		self.get()
		selection = {}
		item = None
		for i in self.listbox.curselection():
			item = self.listbox.get(i)[0]
			for name,value in list(self.layer.items()):
				if name.startswith(item):
					selection[name] = value
		if item is None: return

		# Write to clipboard
#		sio = StringIO()
		sio = io.BytesIO()
		#sio.write(_LAYER_CLIP)
		pickler = pickle.Pickler(sio)
		pickler.dump(_LAYER_CLIP)
		pickler.dump(item)
		pickler.dump(selection)
		self.clipboard_clear()
		self.clipboard_append(binascii.b2a_hex(sio.getvalue()), type='STRING')
#		self.clipboard_append(sio.getvalue())
		return "break"

	# ----------------------------------------------------------------------
	def cut(self, event=None):
		return "break"

	# ----------------------------------------------------------------------
	def paste(self, event=None):
		self.get()

		try: clipboard = self.selection_get(selection='CLIPBOARD')
		except: return
		if clipboard.startswith(_LAYER_CLIP):
			# Pickle format
#			unpickler = pickle.Unpickler(StringIO(clipboard[len(_LAYER_CLIP):]))
			unpickler = pickle.Unpickler(io.BytesIO(clipboard[len(_LAYER_CLIP):]))
			try:
				item      = unpickler.load()
				selection = unpickler.load()
			except:
				return "break"
		else:
			self.flair.notify("Cannot paste","Invalid content in clipboard",tkFlair.NOTIFY_WARNING)
			return "break"

		self.layer.var.update(selection)
		self.showPage(item)
		return "break"

	# --------------------------------------------------------------------
	# Display new page
	# --------------------------------------------------------------------
	def listboxChanged(self, event):
		sel = self.listbox.get(self.listbox.curselection())
		if sel and sel[0] != _ADD:
			self.showPage(sel[0])

	# --------------------------------------------------------------------
	# Populate listbox with visible and active layers
	# --------------------------------------------------------------------
	def populateListbox(self):
		self.listbox.delete(0,END)
		self.listbox.insert(END, (_ORDER[0], _ON))
		for page in _ORDER[1:]:
			on = option(self.layer, page)
			if on or option(self.layer, page+'._visible', False):
				self.layer[page+'._visible'] = True
				if on:
					onoff = _ON
				else:
					onoff = _OFF
				self.listbox.insert(END, (page, onoff))
		self.listbox.insert(END,(_ADD,''))

		if int(self.layer.get('_system',False)):
			self.globalButton.config(state=DISABLED)
		else:
			self.globalButton.config(state=NORMAL)
		self.listbox.itemconfigure(0, fg='DarkRed')

	# --------------------------------------------------------------------
	# Add layer to listbox
	# --------------------------------------------------------------------
	def add2Listbox(self, event=None):
		self.flair.setModified()
		# Add only on the non-system layers
		self.get()

		# Layers that could be added
		layers = []
		# Add only one of these layers
		usrbin  = False
		palette = False
		for page in _ORDER[1:]:
			if not option(self.layer, page) and \
			   not option(self.layer, page+"._visible", False):
				if page.startswith("Usrbin"):
					if usrbin: continue
					usrbin = True
				elif page.startswith("Palette"):
					if palette: continue
					palette = True
				layers.append(page)
		edit = tkExtra.InPlaceList(self.listbox.lists[0], values=layers)

		if edit and edit.value is not None:
			self.layer[edit.value] = True
			self.layer[edit.value+"._visible"] = True
			self.populateListbox()
			self.showPage(edit.value)
		return 'break'

	# --------------------------------------------------------------------
	def cloneListbox(self, event=None):
		if not self.activepage.startswith("Usrbin") and \
		   not self.activepage.startswith("Palette"):
			return
		self.flair.setModified()

		# Search what to add
		if self.activepage.startswith("Usrbin"):
			for page in USRBINNAMES:
				if not option(self.layer, page) and \
				   not option(self.layer, page+"._visible", False):
					break
			else:
				return

		elif self.activepage.startswith("Palette"):
			for page in PALETTENAMES:
				if not option(self.layer, page) and \
				   not option(self.layer, page+"._visible", False):
					break
			else:
				return

		self.layer[page] = True
		self.layer[page+"._visible"] = True
		self.populateListbox()
		# Don't call showPage() but simply activate listbox
		for i,p in enumerate(self.listbox.lists[0].get(0,END)):
			if p==page:
				self.listbox.selection_set(i)
				self.listbox.activate(i)
				break
		self.activepage = page

	# --------------------------------------------------------------------
	def removeListbox(self, event=None):
		if self.activepage == _ORDER[0]: return
		self.flair.setModified()

		# find selected item
		active = self.listbox.index(ACTIVE)

		self.layer[self.activepage] = False
		self.layer[self.activepage + "._visible"] = False
		self.populateListbox()

		try:
			page = self.listbox.lists[0].get(active)
			self.showPage(page)
		except:
			pass

	# --------------------------------------------------------------------
	def showPage(self, page):
		if not page or \
		   not option(self.layer, page+"._visible"):
			page = _ORDER[0]

		for i,p in enumerate(self.listbox.lists[0].get(0,END)):
			if p==page:
				self.listbox.selection_set(i)
				self.listbox.activate(i)
				break

		if self.activepage is not None and \
		   self.activepage != page:
			self.usrbinUpdate()
			self.paletteUpdate()
			if self.activepage == "Voxel": self.voxelUpdate()
			try: self.pages[self.activepage].pack_forget()
			except KeyError: pass

		self.activepage = page
		self.title.config(text=page)
		if page == '3D': self.d3Refresh()
		elif page == "Voxel": self.voxelRefresh()
		self.usrbinInit()
		self.paletteInit()

		try: self.pages[page].pack(expand=YES, fill=BOTH)
		except KeyError: pass
		self._scrollFrame.updateScrollRegion()

	# --------------------------------------------------------------------
	def listbox0Release(self, event):
		if len(self.listbox.curselection())>1: return

		active = self.listbox.nearest(event.y)
		self.listbox.activate(active)
		if active and self.listbox.lists[0].get(active) == _ADD:
			self.add2Listbox()

	# --------------------------------------------------------------------
	def listbox1Release(self, event):
		if len(self.listbox.curselection())>1: return

		active = self.listbox.nearest(event.y)
		if active is not None:
			self.listbox.activate(active)
		else:
			return

		self.selectToggle()

	# --------------------------------------------------------------------
	def selectToggle(self, event=None):
		lst = self.listbox.lists[1]
		for i in map(int,lst.curselection()):
			if i==0: continue
			sel = lst.get(i)
			if not sel: continue
			if sel[1] == ' ':
				sel = _ON
			else:
				sel = _OFF
			lst.set(i,sel)
			lst.selection_set(i)
		self.applyLate()

	# --------------------------------------------------------------------
	# Populate the font list
	# Scan the .tga in the system/font directory then ini/font
	# --------------------------------------------------------------------
	def populateFonts(self):
		self.fonts = []

		fonts = {}
		# User defined tga fonts
		try:
			for f in os.listdir(os.path.join(tkFlair.prgDir,'fonts')):
				fn,ext = os.path.splitext(f)
				if ext=='.tga': fonts[fn] = True
		except OSError:
			messagebox.showwarning('Fonts directory not found',
				'Please check your flair installation. the "fonts" ' \
				'directory is not found!')
			return

		# flair default fonts
		try:
			for f in os.listdir(os.path.join(tkFlair.iniDir,'fonts')):
				fn,ext = os.path.splitext(f)
				if ext=='.tga': fonts[fn] = True
		except OSError:
			pass

		# Try system bitmap fonts
		#for f in subprocess.check_output('xlsfonts').splitlines():
		#	if f[0]=='-': continue
		#	fonts[f] = True

		self.fonts = list(fonts.keys())
		self.fonts.sort()

	# --------------------------------------------------------------------
	# Create options frame
	# --------------------------------------------------------------------
	def optionsFrame(self):
		pageframe = Frame(self.mainframe)

		self.axes       = BooleanVar(self)
		self.titleVar   = BooleanVar(self)
		self.showErrors = BooleanVar(self)
		self.gridVar    = BooleanVar(self)
		self.viewport   = BooleanVar(self)

		# ---
		row = 0
		b = Checkbutton(pageframe, text='Title',
				variable=self.titleVar,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Draw title of plot on top right corner')
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		b = Checkbutton(pageframe, text='Coordinate system',
				variable=self.axes,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Draw coordinate system on bottom left corner')
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		b = Checkbutton(pageframe, text='Viewport lines',
				variable=self.viewport,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Display other viewport cross sections')
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		b = Checkbutton(pageframe, text='Show Errors',
				variable=self.showErrors,
				onvalue=1, offvalue=0,
				foreground="Red",
				activeforeground="Red",
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b, 'Show Errors. Disabling this you will not be notified for geometry errors')
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		b = Checkbutton(pageframe, text='Grid',
				variable=self.gridVar,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Draw grid')
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		Label(pageframe, text='Grid Level').grid(
			row=row, column=0, sticky=W)
		self.gridlevel = Scale(pageframe, from_=0, to_=32,
				orient=HORIZONTAL, showvalue=0,
				command=self.apply)
		tkExtra.Balloon.set(self.gridlevel,
			'Set the grid lines intensity level')
		self.gridlevel.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Lattice Level').grid(
			row=row, column=0, sticky=W)
		self.latticeLevel = Scale(pageframe, from_=0, to_=32,
				orient=HORIZONTAL, showvalue=0,
				command=self.apply)
		tkExtra.Balloon.set(self.latticeLevel,
			'Set the lattice&voxel hash lines intensity level')
		self.latticeLevel.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Crosshair:').grid(
			row=row, column=0, sticky=W)
		self.crosshair = Scale(pageframe, from_=0, to_=32,
				orient=HORIZONTAL, showvalue=0,
				command=self.apply)
		tkExtra.Balloon.set(self.crosshair,
			'Draw crosshair at the center of the plot')
		self.crosshair.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Text Background').grid(
			row=row, column=0, sticky=W)
		self.background = Scale(pageframe, from_=0, to_=255,
				orient=HORIZONTAL, showvalue=0,
				command=self.apply)
		tkExtra.Balloon.set(self.background,
			'Text background level')
		self.background.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Short cut:').grid(
			row=row, column=0, sticky=W)
		self.shortCut = tkExtra.Combobox(pageframe, width=10)
		tkExtra.Balloon.set(self.shortCut,
			'Select keyboard shortcut to access this layer')
		self.shortCut.grid(row=row, column=1, sticky=EW)
		self.shortCut.fill(SHORTCUTS)

		# ---
		row += 1
		Label(pageframe, text='General Font:').grid(
			row=row, column=0, sticky=W)
		self.generalFont = tkExtra.Combobox(pageframe, width=10,
					command=self.apply)
		tkExtra.Balloon.set(self.generalFont,
			'Bitmap font used in Geometry (located in [flair]/fonts folder')
		self.generalFont.grid(row=row, column=1, sticky=EW)
		self.generalFont.fill(self.fonts)

		# ---
		row += 1
		Label(pageframe, text='Grid Font:').grid(
			row=row, column=0, sticky=W)
		self.gridFont = tkExtra.Combobox(pageframe, width=10,
					command=self.apply)
		tkExtra.Balloon.set(self.gridFont,
			'Bitmap font used in Geometry (located in [flair]/fonts folder')
		self.gridFont.grid(row=row, column=1, sticky=EW)
		self.gridFont.fill(self.fonts)

		# ---
		row += 1
		Label(pageframe, text='Palette Font:').grid(
			row=row, column=0, sticky=W)
		self.paletteFont = tkExtra.Combobox(pageframe, width=10,
					command=self.apply)
		tkExtra.Balloon.set(self.paletteFont,
			'Bitmap font used in Geometry (located in [flair]/fonts folder')
		self.paletteFont.grid(row=row, column=1, sticky=EW)
		self.paletteFont.fill(self.fonts)

		# ---
		pageframe.grid_columnconfigure(1, weight=1)

		return pageframe

	# --------------------------------------------------------------------
	# Show frame
	# --------------------------------------------------------------------
	def showFrame(self):
		pageframe = Frame(self.mainframe)

		self.showlattice = BooleanVar(self)
		self.showvertex  = BooleanVar(self)
		self.showvoxel   = BooleanVar(self)

		# ---
		row = 0
		b = Checkbutton(pageframe, text='Vertices',
				variable=self.showvertex,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Display vertices on the intersections of selected bodies')
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		l = Label(pageframe, text='Labels:')
		l.grid(row=row, column=0, sticky=W)

		self.labels = tkExtra.Combobox(pageframe, width=10,
					command=self.apply)
		self.labels.grid(row=row, column=1, sticky=EW)
		self.labels.fill(LABELS)

		# ---
		row += 1
		l = Label(pageframe, text='Color:')
		l.grid(row=row, column=0, sticky=W)

		self.colorvalue = tkExtra.Combobox(pageframe, width=10,
				command=self.check4ShowPalette)
		self.colorvalue.grid(row=row, column=1, sticky=EW)
		self.colorvalue.fill(COLORVALUES)

		# ---
		row += 1
		b = Checkbutton(pageframe, text='Lattices',
				variable=self.showlattice,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b, 'Show lattice regions')
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		b = Checkbutton(pageframe, text='Voxel',
				variable=self.showvoxel,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b, 'Show voxel regions')
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1
		l = Label(pageframe, text='Palette:')
		l.grid(row=row, column=0, sticky=W)

		self.showpalette = tkExtra.Combobox(pageframe, width=10,
				command=self.check4ShowPalette)
		self.showpalette.grid(row=row, column=1, sticky=EW)
		self.showpalette.fill(PALETTENAMES)

		# ---
		pageframe.grid_columnconfigure(1, weight=1)

		return pageframe

	# --------------------------------------------------------------------
	def check4ShowPalette(self, event=None):
		if self._init:
			self.apply()
			return

		try: idx = COLORVALUES.index(self.colorvalue.get())
		except: return
		if idx < 4:
			self.apply()
			return

		# Check if palette is visible otherwise added
		palette = self.showpalette.get()

		if not option(self.layer, palette+'._visible', False):
			self.layer[palette] = True
			self.layer[palette+'._visible'] = True
			self.populateListbox()
			self.showPage(self.activepage)
		self.apply()

	# --------------------------------------------------------------------
	# Background Image Frame
	# --------------------------------------------------------------------
	def imageFrame(self):
		pageframe = Frame(self.mainframe)

		self.imagePrompt = BooleanVar(self)

		# ---
		row = 0
		Label(pageframe, text='Image:').grid(
				row=row, column=0, sticky=W)

		self.imageFilename = ''
		self.backgroundImage = Button(pageframe,
						image=tkFlair.icons["load"],
						text='',
						compound=LEFT,
						command=self.loadImage)
		tkExtra.Balloon.set(self.backgroundImage,
			'Select background image')
		self.backgroundImage.grid(row=row, column=1, sticky=EW)

		pageframe.grid_columnconfigure(1, weight=1)

		# ---
		if CalibrateImage.Image is None:
		# XXX FIXME give a message
			state = DISABLED
		else:
			state = NORMAL

		# ---
		row += 1
		b = Button(pageframe,	text='Calibrate',
					image=tkFlair.icons['edit'],
					compound=LEFT,
					command=self.calibrateImage,
					state=state)
		b.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(b, 'Calibrate image. Requires PIL package')

		# ---
		row += 1
		Label(pageframe, text='Transparency:').grid(
			row=row, column=0, sticky=W)
		self.imageAlpha = Scale(pageframe, from_=0, to_=255,
				orient=HORIZONTAL, showvalue=0,
				command=self.apply)
		tkExtra.Balloon.set(self.imageAlpha,
			'Set the transparency level of the superimposed plot')
		self.imageAlpha.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Color Adjust:').grid(
				row=row, column=0, sticky=W)

		# ---
		row += 1
		Label(pageframe, text='Black:').grid(
				row=row, column=0, sticky=W)

		self.blackImageColor = Button(pageframe,text='  ',
				background='#000000', activebackground='#000000',
				command=self.chooseBlackColor)
		tkExtra.Balloon.set(self.blackImageColor,
			'Choose color for black pixels')
		self.blackImageColor.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='White:').grid(
				row=row, column=0, sticky=W)

		self.whiteImageColor = Button(pageframe,text='  ',
				background='#FFFFFF', activebackground='#FFFFFF',
				command=self.chooseWhiteColor)
		tkExtra.Balloon.set(self.whiteImageColor,
			'Choose color for black pixels')
		self.whiteImageColor.grid(row=row, column=1, sticky=EW)

		pageframe.grid_columnconfigure(1, weight=1)

		# ---
		row += 1
		Label(pageframe, text='Prompt draw', foreground='DarkRed').grid(
				row=row, column=0, sticky=W)

		b = Checkbutton(pageframe, text='',
				variable=self.imagePrompt,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Default (Off) is to draw image only when the display is inactive. '\
			'Enable to draw the image immediately, ' \
			'BUT slowing down the refreshing')
		b.grid(row=row, column=1, sticky=W)

		return pageframe

	# --------------------------------------------------------------------
	def loadImage(self):
		filename = bFileDialog.askopenfilename(master=self,
			title='Load background image',
			initialfile=self.imageFilename,
			filetypes=[
				('All images', ('*.png', '*.gif', '*.jpg', '*.jpeg')),
				('PNG - Portable Network Graphics', '*.png'),
				('GIF - Graphics Interchange Format', '*.gif'),
				('JPG - Joint Photographic Experts Group', ('*.jpg','*.jpeg')),
				('All','*')])

		if len(filename) > 0: self.setImagefile(filename)

	# --------------------------------------------------------------------
	def setImagefile(self, filename):
		self.imageFilename = filename
		self.backgroundImage['text'] = os.path.basename(filename)

	# --------------------------------------------------------------------
	def calibrateImage(self):
		if self.imageFilename:
			dlg = CalibrateImage.Calibrate(self)
			dlg.loadImage(self.imageFilename)

			# prepare markers
			markers = getMarkers(self.layer.var)

			dlg.setMarkers(markers)

			if dlg.show():
				# Calibration changed save to input and flair
				markers = dlg.getMarkers()

				# Save to layers and create points out of the markers
				points = self.project.input['!point'][:]
				if points: pos = points[-1].pos()+1
				else: pos = self.project.input.bestPosition('!point')

				self.layer['Image.marker'] = len(markers)

				undoinfo = [self.flair.refreshUndo(), self.flair.renumberUndo()]
				for i,mark in enumerate(markers):
					self.layer['Image.marker.%d'%(i)] = \
						" ".join(map(str,mark))
					xi,yi,x,y,z = mark
					for card in points:
						if abs(card.numWhat(1)-x)<1e-10 and \
						   abs(card.numWhat(2)-y)<1e-10 and \
						   abs(card.numWhat(3)-z)<1e-10:
							break
					else:
						# Add a new point
						j=0
						while True:	# Find unique name
							j += 1
							name = 'marker%d'%(j)
							for c in points:
								if c.sdum()==name: break
							else: break

						card = Input.Card('!point',[name,x,y,z])
						card[SELECT] = BIT_VISIBLE
						undoinfo.append(self.flair.addCardUndo(card, pos, False))
						#self.project.input.addCard(card, pos, True)
						points.append(card)
						#self.project.setInputModified()
						pos += 1
				undoinfo.append(self.flair.renumberUndo())
				undoinfo.append(self.flair.refreshUndo())
				self.flair.addUndo(undo.createListUndo(undoinfo,'Add Image markers'))
				self.flair.setInputModified()
				self.apply()

	# --------------------------------------------------------------------
	def setImageBlack(self, col):
		if isinstance(col,int):
			color = '#%06X'%(col)
		else:
			color = col

		if self.blackImageColor['background'] != color:
			self.blackImageColor['background'] = color
			self.blackImageColor['activebackground'] = color
			self.apply()

	# --------------------------------------------------------------------
	def setImageWhite(self, col):
		if isinstance(col,int):
			color = '#%06X'%(col)
		else:
			color = col
		if self.whiteImageColor['background'] != color:
			self.whiteImageColor['background'] = color
			self.whiteImageColor['activebackground'] = color
			self.apply()

	# --------------------------------------------------------------------
	def chooseBlackColor(self):
		try:
			rgb, color = askcolor(title='Select Black color',
				parent=self,
				initialcolor=self.blackImageColor['background'])
			self.setImageBlack(color)
		except TclError:
			pass

	# --------------------------------------------------------------------
	def chooseWhiteColor(self):
		try:
			rgb, color = askcolor(title='Select White color',
				parent=self,
				initialcolor=self.whiteImageColor['background'])
			self.setImageWhite(color)
		except TclError:
			pass

	# --------------------------------------------------------------------
	# Usrdump Frame
	# --------------------------------------------------------------------
	def userdumpFrame(self):
		pageframe = Frame(self.mainframe)

		# ---
		row = 0
		Label(pageframe, text='Userdump file:').grid(
				row=row, column=0, sticky=W)

		self.userdumpFilename = ''
		self.userdumpButton = Button(pageframe,text='',
				command=self.loadUserdump)
		tkExtra.Balloon.set(self.userdumpButton,
			'Select userdump file')
		self.userdumpButton.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Start:').grid(
				row=row, column=0, sticky=W)

		self.userdumpStart = IntVar()
		b = tkExtra.Spinbox(pageframe,
					text = self.userdumpStart,
					from_=1,
					to=sys.maxsize,
					command=self.applyLate,
					background='white')
		tkExtra.Balloon.set(b, 'Starting event for userdump')
		b.bind("<Return>",   self.apply)
		b.bind("<KP_Enter>", self.apply)
		b.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='# Events:').grid(
				row=row, column=0, sticky=W)

		self.userdumpEvents = IntVar()
		b = tkExtra.Spinbox(pageframe,
					text = self.userdumpEvents,
					from_=0,
					to=10000,
					command=self.applyLate,
					background='white')
		tkExtra.Balloon.set(b, 'Userdump events to load')
		b.bind("<Return>",   self.apply)
		b.bind("<KP_Enter>", self.apply)
		b.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		self.userdumpList = tkExtra.MultiListbox(pageframe,
					(("  ",       3, None),
					 ("Particle", 8, None),
					 ("Emin",     5, None),
					 ("Emax",     5, None),
					 ("Color",    8, None)),
					 height=5,
					 background="White")
		self.userdumpList.bindList("<Button-1>",        self.userdumpClick)
		self.userdumpList.bindList("<ButtonRelease-1>", self.userdumpRelease)
		self.userdumpList.grid(row=row, column=0, columnspan=2, sticky=NSEW)
		self.userdumpList.setPopupMenu([
				("Show", 0, self.userdumpShow, tkFlair.icons["enable"]),
				("Hide", 0, self.userdumpHide, tkFlair.icons["disable"])])
		pageframe.columnconfigure(1,weight=1)
		pageframe.rowconfigure(row,weight=1)

		for particle in Input.Particle.list:
			if not particle or particle[0]=='@': continue
			self.userdumpList.insert(END, (_OFF, particle, "", "", ""))

		return pageframe

	# --------------------------------------------------------------------
	def loadUserdump(self):
		filename = bFileDialog.askopenfilename(master=self,
			title='Load USERDUMP file',
			initialfile=self.userdumpFilename,
			filetypes=[
			        ('Fluka USERDUMP files',('*dump*','*draw*')),
				('Fluka fort files',('*_ftn.*','*_fort.*')),
				('All','*')])

		if len(filename) > 0:
			self.setUserdumpFile(filename)
			self.apply()

	# --------------------------------------------------------------------
	def setUserdumpFile(self, filename):
		self.userdumpFilename = filename
		self.userdumpButton['text'] = os.path.basename(filename)

	# --------------------------------------------------------------------
	def userdumpClick(self, event):
		y = event.widget.nearest(event.y)
		if event.widget is self.userdumpList.lists[0]:
			if event.widget.selection_includes(y):
				sel = event.widget.curselection()
				if len(sel)>1:
					status = event.widget.get(y)
					if status == _ON:
						status = _OFF
					else:
						status = _ON
					for i in sel:
						event.widget.set(i, status)
				self.applyLate()
				return "break"

	# --------------------------------------------------------------------
	def userdumpRelease(self, event):
		y = event.widget.nearest(event.y)
		if event.widget is self.userdumpList.lists[0]:
			sel = event.widget.curselection()
			if len(sel)==1:
				status = event.widget.get(y)
				if status == _ON:
					status = _OFF
				else:
					status = _ON
				event.widget.set(y, status)
				self.applyLate()
			return "break"

	# --------------------------------------------------------------------
	def userdumpShow(self, event=None):
		sel = self.userdumpList.curselection()
		if not sel:
			for i in range(self.userdumpList.size()):
				self.userdumpList.lists[0].set(i,_ON)
		else:
			for i in sel:
				self.userdumpList.lists[0].set(i,_ON)
		self.applyLate()

	# --------------------------------------------------------------------
	def userdumpHide(self, event=None):
		sel = self.userdumpList.curselection()
		if not sel:
			for i in range(self.userdumpList.size()):
				self.userdumpList.lists[0].set(i,_OFF)
		else:
			for i in sel:
				self.userdumpList.lists[0].set(i,_OFF)
		self.applyLate()

	# --------------------------------------------------------------------
	# Usrbin Frame
	# --------------------------------------------------------------------
	def usrbinFrame(self):
		pageframe = Frame(self.mainframe)

		self.usrbinInput = BooleanVar(self)

		# ---
		row = 0
		b = Checkbutton(pageframe, text='USRBIN from Input',
				variable=self.usrbinInput,
				command=self.usrbinInputChanged,
				onvalue=1, offvalue=0,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Unselected load a USRBIN file with data to display. ' \
			'Select to display a USRBINs grid from Input as checker.')
		b.grid(row=row, column=0, columnspan=2, sticky=W)

		# ---
		row += 1

		self.usrbinLabel = Label(pageframe, text='File:')
		self.usrbinLabel.grid(row=row, column=0, sticky=W)

		self.usrbinFile = Entry(pageframe, width=5, foreground="DarkBlue", background="White")
		tkExtra.Balloon.set(self.usrbinFile, "USRBIN binary output file")
		self.usrbinFile.grid(row=row, column=1, sticky=EW)

		self.usrbinFile.bind('<Return>',   self.usrbinFileFocusOut)
		self.usrbinFile.bind('<KP_Enter>', self.usrbinFileFocusOut)
		self.usrbinFile.bind('<FocusOut>', self.usrbinFileFocusOut)

		self.usrbinButton = Button(pageframe,
					image=tkFlair.icons["load"],
					command=self.loadUsrbin)
		tkExtra.Balloon.set(self.usrbinButton, 'Load usrbin file')
		self.usrbinButton.grid(row=row, column=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Detector:').grid(
				row=row, column=0, sticky=W)

		self.usrbinDet = tkExtra.Combobox(pageframe,
				width=10,
				command=self.usrbinDetChanged)
		tkExtra.Balloon.set(self.usrbinDet,
			'Usrbin detector to load')
		self.usrbinDet.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		self.usrbinNormLabel = Label(pageframe, text='Norm:')
		self.usrbinNormLabel.grid(row=row, column=0, sticky=W)
		self.usrbinNorm = Entry(pageframe, background='white')
		self.usrbinNorm.bind('<Return>',   self.apply)
		self.usrbinNorm.bind('<KP_Enter>', self.apply)
		tkExtra.Balloon.set(self.usrbinNorm,
			'Usrbin normalisation. Default=1.0')
		self.usrbinNorm.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		self.usrbinXofs = Label(pageframe, text='X-offset:')
		self.usrbinXofs.grid(row=row, column=0, sticky=W)
		self.usrbinXofs = tkExtra.FloatEntry(pageframe, background='white')
		self.usrbinXofs.bind('<Return>',   self.apply)
		self.usrbinXofs.bind('<KP_Enter>', self.apply)
		tkExtra.Balloon.set(self.usrbinXofs,
			'Usrbin X offset. Default=0.0')
		self.usrbinXofs.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		self.usrbinYofs = Label(pageframe, text='Y-offset:')
		self.usrbinYofs.grid(row=row, column=0, sticky=W)
		self.usrbinYofs = tkExtra.FloatEntry(pageframe, background='white')
		self.usrbinYofs.bind('<Return>',   self.apply)
		self.usrbinYofs.bind('<KP_Enter>', self.apply)
		tkExtra.Balloon.set(self.usrbinYofs,
			'Usrbin Y offset. Default=0.0')
		self.usrbinYofs.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		self.usrbinZofs = Label(pageframe, text='Z-offset:')
		self.usrbinZofs.grid(row=row, column=0, sticky=W)
		self.usrbinZofs = tkExtra.FloatEntry(pageframe, background='white')
		self.usrbinZofs.bind('<Return>',   self.apply)
		self.usrbinZofs.bind('<KP_Enter>', self.apply)
		tkExtra.Balloon.set(self.usrbinZofs,
			'Usrbin Z offset. Default=0.0')
		self.usrbinZofs.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Rotdefi:').grid(
				row=row, column=0, sticky=W)

		self.usrbinRotdefi = tkExtra.Combobox(pageframe, width=10,
					command=self.apply)
		tkExtra.Balloon.set(self.usrbinRotdefi,
			'Rotdefi to be applied to USRBIN')
		self.usrbinRotdefi.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Rotdefi2:').grid(
				row=row, column=0, sticky=W)

		self.usrbinRotdefi2 = tkExtra.Combobox(pageframe, width=10,
					command=self.apply)
		tkExtra.Balloon.set(self.usrbinRotdefi2,
			'Second Rotdefi to be applied to USRBIN')
		self.usrbinRotdefi2.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Transparency:').grid(
			row=row, column=0, sticky=W)
		self.usrbinAlpha = Scale(pageframe, from_=0, to_=255,
					orient=HORIZONTAL,
					command=self.apply,
					showvalue=0)
		tkExtra.Balloon.set(self.usrbinAlpha,
			'Set the transparency level of the superimposed usrbin plot')
		self.usrbinAlpha.grid(row=row, column=1, columnspan=2, sticky=EW)

		pageframe.grid_columnconfigure(1, weight=1)
		pageframe.columnconfigure(1,weight=1)

		# ---
		row += 1
		l = Label(pageframe, text='Palette:')
		l.grid(row=row, column=0, sticky=W)

		self.usrbinPalette = tkExtra.Combobox(pageframe, width=10,
				command=self.usrbinPaletteChange)
		self.usrbinPalette.grid(row=row, column=1, columnspan=2, sticky=EW)
		self.usrbinPalette.fill(PALETTENAMES)

		return pageframe

	# --------------------------------------------------------------------
	def usrbinPaletteChange(self, event=None):
		# Check if palette is visible otherwise added
		palette = self.usrbinPalette.get()
		self.layer[self.entry('palette')] = palette

		if not option(self.layer, palette+'._visible', False):
			self.layer[palette] = True
			self.layer[palette+'._visible'] = True
			self.populateListbox()
			self.showPage(self.activepage)
		self.apply()

	# --------------------------------------------------------------------
	# Initialize usrbin values set to fields
	# --------------------------------------------------------------------
	def usrbinInit(self):
		if self.activepage is None or \
		   not self.activepage.startswith('Usrbin'): return
#		if not option(self.layer, self.activepage): return

		oldInit = self._init
		self._init = True
		self.usrbin = Data.Usrbin()
		self._usrbinFile = ''
		self.usrbinInput.set(option(self.layer,   self.entry('input')))
		self.usrbinInputChanged()
		self.setUsrbinFile(option(self.layer,     self.entry('file')))
		self.setUsrbindet(option(self.layer,      self.entry('det')))
		self.usrbinNorm.delete(0, END)
		self.usrbinNorm.insert(0, option(self.layer, str(self.entry('norm'))))
		self.usrbinXofs.set(option(self.layer, str(self.entry('xofs'))))
		self.usrbinYofs.set(option(self.layer, str(self.entry('yofs'))))
		self.usrbinZofs.set(option(self.layer, str(self.entry('zofs'))))
		self.usrbinRotdefi.clear()
		self.usrbinRotdefi.fill(self.project.input.rotdefiList(True))
		self.usrbinRotdefi.set(option(self.layer, self.entry('rotdefi')))
		self.usrbinRotdefi2.clear()
		self.usrbinRotdefi2.fill(self.project.input.rotdefiList(True))
		self.usrbinRotdefi2.set(option(self.layer, self.entry('rotdefi2')))
		self.usrbinAlpha.set(option(self.layer,   self.entry('alpha')))
		self.usrbinPalette.set(option(self.layer, self.entry('palette')))
		self._init = oldInit

	# --------------------------------------------------------------------
	# usrbinUpdate: Copies the current values on the interface to the
	# project data (only updates current usrbin)
	# --------------------------------------------------------------------
	def usrbinUpdate(self):
		if not self.activepage.startswith('Usrbin'): return
#		if not option(self.layer, self.activepage): return

		self.layer[self.entry('input')] = self.usrbinInput.get()
		try:
			self.layer[self.entry('file')] = self._usrbinFile
		except OSError:
			pass

		try:    det = int(self.usrbinDet.get().split()[0])
		except: det = 1
		self.layer[self.entry('det')] = det

		try:    self.layer[self.entry('norm')] = self.usrbinNorm.get()
		except: self.layer[self.entry('norm')] = 1.0
		try:    self.layer[self.entry('xofs')] = self.usrbinXofs.get()
		except: self.layer[self.entry('xofs')] = ""
		try:    self.layer[self.entry('yofs')] = self.usrbinYofs.get()
		except: self.layer[self.entry('yofs')] = ""
		try:    self.layer[self.entry('zofs')] = self.usrbinZofs.get()
		except: self.layer[self.entry('zofs')] = ""
		self.layer[self.entry('rotdefi')] = self.usrbinRotdefi.get()
		self.layer[self.entry('rotdefi2')]     = self.usrbinRotdefi2.get()
		self.layer[self.entry('alpha')]   = self.usrbinAlpha.get()
		self.layer[self.entry('palette')] = self.usrbinPalette.get()

	# --------------------------------------------------------------------
	def usrbinInputChanged(self, event=None):
		#self.setUsrbinFile('')
		self.populateUsrbinDet()
		if self.usrbinInput.get():
			self.usrbinLabel.config(state=DISABLED)
			self.usrbinFile.config(state=DISABLED)
			self.usrbinButton.config(state=DISABLED)
			self.usrbinNormLabel.config(state=DISABLED)
			self.usrbinNorm.config(state=DISABLED)
			self.usrbinXofs.config(state=DISABLED)
			self.usrbinYofs.config(state=DISABLED)
			self.usrbinZofs.config(state=DISABLED)
		else:
			self.usrbinLabel.config(state=NORMAL)
			self.usrbinFile.config(state=NORMAL)
			self.usrbinButton.config(state=NORMAL)
			self.usrbinNormLabel.config(state=NORMAL)
			self.usrbinNorm.config(state=NORMAL)
			self.usrbinXofs.config(state=NORMAL)
			self.usrbinYofs.config(state=NORMAL)
			self.usrbinZofs.config(state=NORMAL)

	# --------------------------------------------------------------------
	def populateUsrbinDet(self):
		oldInit = self._init
		self._init = True
		self.usrbinDet.clear()

		if self.usrbinInput.get():
			self.usrbinDet.delete(0,END)
			# Populate with input usrbins + EVENTBIN
			for idx, card in enumerate(self.project.input.cardsSorted('USRBIN')):
				self.usrbinDet.insert(END, '%d %s [%02d]' % \
					(idx+1, card.sdum(), abs(card.intWhat(3))))
		else:
			# --- Load USRBIN ---
			if self._usrbinFile:
				say("Loading:", self._usrbinFile)
				try:
					self.usrbin.readHeader(self._usrbinFile)
				except:
					# FIXME should display a nice message
					# Displaying an error before the parent (from the __init__)
					# is displayed makes funny things!!!
					#messagebox.showerror('Error loading USRBIN',
					#	'Error loading USRBIN file: %s'%(sys.exc_info()[1]),
					#	parent=self)
					say('ERROR loading USRBIN file:',sys.exc_info()[1])
					return

				# Populate the detector menu
				self.usrbinDet.delete(0,END)
				for i,det in enumerate(self.usrbin.detector):
					self.usrbinDet.insert(END, '%d %s' % (i+1, det.name))
		self._init = oldInit

	# ----------------------------------------------------------------------
	def usrbinFileFocusOut(self, event=None):
		if self._usrbinFile == self.usrbinFile.get(): return
		self.setUsrbinFile(self.usrbinFile.get())
		self.apply()

	# --------------------------------------------------------------------
	def loadUsrbin(self):
		filename = bFileDialog.askopenfilename(master=self,
			title='Load USRBIN file',
			initialfile=self._usrbinFile,
			filetypes=[
			        ('Fluka USRBIN files',('*usrbin*','*.bnn')),
			        #('%s USRBIN files'%(inpname),'%s*usrbin*'%(inpname)),
				('Fluka fort files',('*_ftn.*','*_fort.*')),
				('All','*')])

		if len(filename) > 0:
			self.setUsrbinFile(filename)
			self.apply()

	# --------------------------------------------------------------------
	def setUsrbinFile(self, filename):
		self._usrbinFile= self.project.relativePath(filename)
		self.usrbinFile.delete(0,END)
		self.usrbinFile.insert(0, self._usrbinFile)
		self.populateUsrbinDet()

	# --------------------------------------------------------------------
	def setUsrbindet(self, det):
		try: det = int(det)
		except: det = 1

		oldInit = self._init
		self._init = True
		if self.usrbinInput.get():
			try:
				card = self.project.input.cardsSorted('USRBIN')[det-1]
				self.usrbinDet.set('%d %s [%02d]' %
					(det, card.sdum(), abs(card.intWhat(3))))
			except IndexError:
				self.usrbinDet.set('%d ? [??]'%(det))
		else:
			if 0 < det <= len(self.usrbin.detector):
				detector = self.usrbin.detector[det-1]
				newdet = '%d %s' % (det, detector.name)
				if self.usrbinDet.get() != newdet:
					# FIXME self.apply()
					self.usrbinDet.set(newdet)
			else:
				if self.usrbinDet.get():
					# FIXME self.apply()
					self.usrbinDet.clear()
		self._init = oldInit

	# --------------------------------------------------------------------
	# usrbin detector changed, check if there is rotation associated with it
	# --------------------------------------------------------------------
	def usrbinDetChanged(self):
		if self.usrbinInput.get():
			match = re.match(r'(\d+) (.+) \[(\d+)\]', self.usrbinDet.get())
		else:
			match = re.match(r'(\d+) (.+)', self.usrbinDet.get())

		oldInit = self._init
		self._init = True
		if match:
			idx  = int(match.group(1))
			name = match.group(2)

			# Check ROTPRBIN
			usrbins = self.project.input.cardsSorted('USRBIN')
			names   = [x.name() for x in usrbins]
			self.usrbinRotdefi.set("")
			self.usrbinRotdefi2.set("")
			for card in self.project.input['ROTPRBIN']:
				if card.ignore(): continue
				from_ = card.what(4)
				to_   = card.what(5)
				step  = card.intWhat(6)

				try:    ifrom = names.index(from_)
				except: ifrom = -1
				try:    ito = names.index(to_)
				except: ito = ifrom

				try:
					if ifrom <= names.index(name) <= ito:
						self.usrbinRotdefi.set(card.what(2))
						rot2 = card.sdum()
						if rot2 != "RHO-FACT":
							self.usrbinRotdefi2.set(rot2)
						break
				except:
					pass
			#else:
			#	self.usrbinRotdefi.set("")

			# Check USRBIN in input (maybe)
			if self.usrbinInput.get():
				self.usrbinXofs.config(state=NORMAL)
				self.usrbinYofs.config(state=NORMAL)
				self.usrbinZofs.config(state=NORMAL)
				self.usrbinXofs.set("")
				self.usrbinYofs.set("")
				try:
					card = usrbins[idx-1]
					if card.intWhat(1) in (1,7,11,17):
						self.usrbinXofs.set(card.numWhat(8))
						self.usrbinYofs.set(card.numWhat(5))
				except:
					pass
				self.usrbinXofs.config(state=DISABLED)
				self.usrbinYofs.config(state=DISABLED)
				self.usrbinZofs.config(state=DISABLED)
		self._init = oldInit
		self.apply()

	# --------------------------------------------------------------------
	# 3D ray tracing frame
	# --------------------------------------------------------------------
	def d3Frame(self):
		pageframe = Frame(self.mainframe)

		self.deflights       = BooleanVar(self)
		self.shadows         = BooleanVar(self)
		self.edgedetect      = BooleanVar(self)
		self.skip1stblack    = BooleanVar(self)
		self.light           = IntVar(self)
		self.clipnegative    = BooleanVar(self)
		self.clipnegative1   = BooleanVar(self)
		self.clipnegative2   = BooleanVar(self)
		self.usrbinastexture = BooleanVar(self)
		self.usrbinastexture.set(True)

		# ---
		row = 0
		Label(pageframe, text='Projection:').grid(
				row=row, column=0, sticky=W)
		self.projection = tkExtra.Combobox(pageframe, width=12,
					command=self.apply)
		self.projection.fill(_PROJECTIONS)
		tkExtra.Balloon.set(self.projection,
			"Orthographic, Perspective or Combo=Perspective starting from the Orthographic/Viewing plane")
		self.projection.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		self._fovLabel = Label(pageframe, text='FOV')
		self._fovLabel.grid(row=row, column=0, sticky=W)
		self.fov = Scale(pageframe, from_=0.1, to_=160.0,  # FIXME For fishEye 360.0,
				resolution=0.1,
				orient=HORIZONTAL, showvalue=0,
				command=self.apply)
		self.fov.grid(row=row, column=1, columnspan=2, sticky=EW)
		tkExtra.Balloon.set(self.fov,
			"Camera Field of View (Angle in U). "\
			"Enabled only when perspective or combo mode is enabled")

		# ---
		row += 1
		Label(pageframe, text='Def. Lights:').grid(
				row=row, column=0, sticky=W)
		b = Checkbutton(pageframe, text='',
				variable=self.deflights,
				onvalue=1, offvalue=0,
				command=self.apply,
				anchor=W)
		tkExtra.Balloon.set(b, 'Use default (predefined) lights')
		b.grid(row=row, column=1, sticky=W)

		# ---
		row += 1
		Label(pageframe, text='Ambient Light').grid(
			row=row, column=0, sticky=W)
		self.ambient = Scale(pageframe, from_=0, to_=255,
				orient=HORIZONTAL,
				command=self.apply,
				showvalue=1)
		tkExtra.Balloon.set(self.ambient,
			'Set the level of ambient background intensity')
		self.ambient.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Antialias:', foreground='DarkRed').grid(
			row=row, column=0, sticky=W)
		self.antialias = Scale(pageframe, from_=1, to_=4,
				orient=HORIZONTAL,
				command=self.apply,
				showvalue=0)
		self.antialias.grid(row=row, column=1, columnspan=2, sticky=EW)
		tkExtra.Balloon.set(self.antialias,
			'Enable anti-aliasing by super sampling for \n' \
			'each pixel to enhance image quality.\n' \
			'WARNING Increases exponentially the rendering time.')

		# ---
		row += 1
		Label(pageframe, text='Reflections:', foreground='DarkRed').grid(
			row=row, column=0, sticky=W)
		self.reflections = Scale(pageframe, from_=0, to_=4,
				orient=HORIZONTAL,
				command=self.apply,
				showvalue=0)
		self.reflections.grid(row=row, column=1, columnspan=2, sticky=EW)
		tkExtra.Balloon.set(self.reflections,
			'Enable reflections and refractions from objects\n' \
			'WARNING Increases the rendering time.')

		# ---
		row += 1
		Label(pageframe, text='Shadows:', foreground='DarkRed').grid(
			row=row, column=0, sticky=W)
		b = Checkbutton(pageframe, text='',
				variable=self.shadows,
				onvalue=1, offvalue=0,
				command=self.apply,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Enable casting of shadows from lights\n' \
			'WARNING Increases the rendering time.')
		b.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Edge Detection:', foreground='DarkRed').grid(
				row=row, column=0, sticky=W)
		b = Checkbutton(pageframe, text='',
				variable=self.edgedetect,
				onvalue=1, offvalue=0,
				command=self.apply,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Draw edges as black lines.\n' \
			'WARNING Increases the rendering time')
		b.grid(row=row, column=1, sticky=W)

		# ---
		row += 1
		Label(pageframe, text='Skip BLCKHOLE:').grid(
				row=row, column=0, sticky=W)
		b = Checkbutton(pageframe, text='',
				variable=self.skip1stblack,
				onvalue=1, offvalue=0,
				command=self.apply,
				anchor=W)
		tkExtra.Balloon.set(b, 'Do not stop the ray if it starts from BLACKHOLE. Can lead to errors in display.')
		b.grid(row=row, column=1, sticky=W)

		# ---
		row += 1
		Label(pageframe, text='Xray Level').grid(
			row=row, column=0, sticky=W)
		self.xray = Scale(pageframe, from_=0, to_=255,
				orient=HORIZONTAL,
				command=self.apply,
				showvalue=1)
		tkExtra.Balloon.set(self.xray,
			'Automatic transparency for each region')
		self.xray.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Clipped by:').grid(
				row=row, column=0, sticky=W)

		b = Checkbutton(pageframe, text='',
				variable=self.clipnegative,
				onvalue=1, offvalue=0,
				command=self.apply,
				anchor=W)
		tkExtra.Balloon.set(b, 'Use the inverse of the clipping body')
		b.grid(row=row, column=1, sticky=W)

		self.clipbody = tkExtra.Combobox(pageframe,
					width=9,
					command=self.apply)
		tkExtra.Balloon.set(self.clipbody,
			'Clipping body')
		self.clipbody.grid(row=row, column=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Clipped by #2:').grid(
				row=row, column=0, sticky=W)

		b = Checkbutton(pageframe, text='',
				variable=self.clipnegative1,
				onvalue=1, offvalue=0,
				command=self.apply,
				anchor=W)
		tkExtra.Balloon.set(b, 'Use the inverse of the clipping body')
		b.grid(row=row, column=1, sticky=W)

		self.clipbody1 = tkExtra.Combobox(pageframe,
					width=9,
					command=self.apply)
		tkExtra.Balloon.set(self.clipbody1,
			'Second Clipping body')
		self.clipbody1.grid(row=row, column=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Clipped by #3:').grid(
				row=row, column=0, sticky=W)

		b = Checkbutton(pageframe, text='',
				variable=self.clipnegative2,
				onvalue=1, offvalue=0,
				command=self.apply,
				anchor=W)
		tkExtra.Balloon.set(b, 'Use the inverse of the clipping body')
		b.grid(row=row, column=1, sticky=W)

		self.clipbody2 = tkExtra.Combobox(pageframe,
					width=9,
					command=self.apply)
		tkExtra.Balloon.set(self.clipbody2,
			'Third Clipping body ')
		self.clipbody2.grid(row=row, column=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Project body:').grid(
				row=row, column=0, sticky=W)

		self.projectbody = tkExtra.Combobox(pageframe,
					width=9,
					command=self.apply)
		tkExtra.Balloon.set(self.projectbody,
			'USRBIN projection body')
		self.projectbody.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Project body #2:').grid(
				row=row, column=0, sticky=W)

		self.projectbody1 = tkExtra.Combobox(pageframe,
					width=9,
					command=self.apply)
		tkExtra.Balloon.set(self.projectbody1,
			'USRBIN second projection body')
		self.projectbody1.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Project body #3:').grid(
				row=row, column=0, sticky=W)

		self.projectbody2 = tkExtra.Combobox(pageframe,
					width=9,
					command=self.apply)
		tkExtra.Balloon.set(self.projectbody2,
			'USRBIN third projection body')
		self.projectbody2.grid(row=row, column=1, columnspan=2, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Usrbin as texture:').grid(
				row=row, column=0, sticky=W)
		b = Checkbutton(pageframe, text='',
				variable=self.usrbinastexture,
				onvalue=1, offvalue=0,
				command=self.apply,
				anchor=W)
		tkExtra.Balloon.set(b, 'Usrbin color mapping as region texture')
		b.grid(row=row, column=1, sticky=W)

		# ---
		row += 1
		LabelFrame(pageframe, text='Light Information').grid(
			row=row, column=0, columnspan=2, sticky=NSEW)

		pageframe.grid_columnconfigure(2, weight=1)

		return pageframe

	# --------------------------------------------------------------------
	def d3Refresh(self):
		# Populate with clipping bodies
		bodies = list(self.project.input.cardsCache('bodies').keys())
		bodies.sort()
		bodies.insert(0,'')

		self.clipbody.fill(bodies)
		self.clipbody1.fill(bodies)
		self.clipbody2.fill(bodies)
		self.projectbody.fill(bodies)
		self.projectbody1.fill(bodies)
		self.projectbody2.fill(bodies)

	# --------------------------------------------------------------------
	# Voxel structures frame
	# --------------------------------------------------------------------
	def voxelFrame(self):
		pageframe = Frame(self.mainframe)

		# ---
		row = 0
		Label(pageframe, text='ROI Colors:').grid(
			row=row, column=0, sticky=W)
		self.roiAlpha = Scale(pageframe, from_=0, to_=255,
					orient=HORIZONTAL,
					command=self.apply,
					showvalue=0)
		tkExtra.Balloon.set(self.roiAlpha,
			'Set blending level of roi colors')
		self.roiAlpha.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		self.roiList = tkExtra.MultiListbox(pageframe,
					(('Id',    3, None),
					 ('Name', 16, None)),
					#selectmode=BROWSE,
					height=10,
					#stretch="first",
					background="White")
		tkExtra.Balloon.set(self.roiList, 'Select ROIS to be displayed')
		self.roiList.grid(row=row, column=0, columnspan=2, sticky=NSEW)
		self.roiList.bind("<<ListboxSelect>>", self.apply)
		pageframe.grid_rowconfigure(row, weight=1)

		# ---
		row += 1
		b = Button(pageframe, text="Clear", command=self.roiClear)
		b.grid(row=row, column=0, columnspan=2)
		tkExtra.Balloon.set(b, 'Clear roi selection')

		pageframe.grid_columnconfigure(1, weight=1)

		return pageframe

	# --------------------------------------------------------------------
	def voxelRefresh(self):
		self.roiList.delete(0,END)
		for card in self.project.input["VOXELS"]:
			if card.ignore() or card["@voxel"] is None: continue
			voxel = card["@voxel"]
			break
		else:
			return

		rois = list(map(int,str(option(self.layer, 'Voxel.roi')).split()))
		for roi,name in list(voxel.roiName.items()):
			self.roiList.insert(END,(roi, name))
			if roi in rois:
				self.roiList.selection_set(END)

	# --------------------------------------------------------------------
	def voxelUpdate(self):
		rois = []
		for sel in self.roiList.curselection():
			rois.append(str(self.roiList.get(sel)[0]))
		self.layer['Voxel.roi'] = ' '.join(rois)

	# --------------------------------------------------------------------
	def roiClear(self, event=None):
		self.roiList.selection_clear(0,END)
		self.apply()

	# --------------------------------------------------------------------
	# Color band frame
	# --------------------------------------------------------------------
	def paletteFrame(self):
		pageframe = Frame(self.mainframe)

		self.cblog      = BooleanVar(self)
		self.cbinv      = BooleanVar(self)
		self.cbmin      = StringVar()
		self.cbmax      = StringVar()
		self.cbsmooth   = BooleanVar(self)
		self.cbcolors   = IntVar()
		self.cpd        = IntVar()
		self.cbalphamin = BooleanVar(self)
		self.cbalphamax = BooleanVar(self)
		self.cb_min_variable = StringVar()

		# ---
		row = 0
		Label(pageframe, text='Label:').grid(
				row=row, column=0, sticky=W)
		self.cblabel = Entry(pageframe, background="White")
		self.cblabel.bind("<Return>",   self.apply)
		self.cblabel.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.cblabel, 'Palette label')
		self.cblabel.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Palette:').grid(
				row=row, column=0, sticky=W)

		self.cbpalette = tkExtra.Combobox(pageframe, width=10,
					command=self.apply)
		palettes = list(Palette.PALETTE.keys())
		palettes.sort()
		self.cbpalette.fill(palettes)
		tkExtra.Balloon.set(self.cbpalette,
			'Color palette to use')
		self.cbpalette.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Log:').grid(
				row=row, column=0, sticky=W)
		b = Checkbutton(pageframe, text='',
				variable=self.cblog,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Logarithmic or linear')
		b.grid(row=row, column=1, sticky=W)

		# ---
		row += 1
		Label(pageframe, text='Invert:').grid(
				row=row, column=0, sticky=W)
		b = Checkbutton(pageframe, text='',
				variable=self.cbinv,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Reverse order of colors')
		b.grid(row=row, column=1, sticky=W)

		# ---
		row += 1
		Label(pageframe, text='Transparent (<Min):').grid(
				row=row, column=0, sticky=W)
		b = Checkbutton(pageframe, text='',
				variable=self.cbalphamin,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Transparent below the minimum limit')
		b.grid(row=row, column=1, sticky=W)

		# ---
		row += 1
		Label(pageframe, text='Transparent (>Max):').grid(
				row=row, column=0, sticky=W)
		b = Checkbutton(pageframe, text='',
				variable=self.cbalphamax,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b,
			'Transparent above the maximum limit')
		b.grid(row=row, column=1, sticky=W)

		# ---
		row += 1
		Label(pageframe, text='Minimum:').grid(
				row=row, column=0, sticky=W)
		e = tkExtra.FloatEntry(pageframe,
				textvariable=self.cbmin,
				background='white')
		e.bind('<Return>',   self.apply)
		e.bind('<KP_Enter>', self.apply)
		tkExtra.Balloon.set(e,
			'Minimum value for color band')
		e.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Maximum:').grid(
				row=row, column=0, sticky=W)
		e = tkExtra.FloatEntry(pageframe,
				textvariable=self.cbmax,
				background='LightYellow')
		e.bind('<Return>',   self.apply)
		e.bind('<KP_Enter>', self.apply)
		tkExtra.Balloon.set(e,
			'Maximum value for color band')
		e.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(pageframe, text='Smooth:').grid(
				row=row, column=0, sticky=W)
		b = Checkbutton(pageframe, text='',
				variable=self.cbsmooth,
				onvalue=1, offvalue=0,
				command=self.applyLate,
				anchor=W)
		tkExtra.Balloon.set(b, 'Smooth palette')
		b.grid(row=row, column=1, sticky=W)

		# ---
		row += 1
		Label(pageframe, text='Colors:').grid(
				row=row, column=0, sticky=W)
		self.cbcolors.set(30)
		self.cbspin = Spinbox(pageframe, text=self.cbcolors,
					from_=2,
					to=256,
					command=self.applyLate,
					background='LightYellow',width=4)
		self.cbspin.grid(row=row, column=1, sticky=W)
		tkExtra.Balloon.set(self.cbspin, 'Total number of colors')

		# ---
		row += 1
		Label(pageframe, text='CPD:').grid(
				row=row, column=0, sticky=W)
		self.cpd_spin = Spinbox(pageframe,
					text=self.cpd,
					from_=1,
					to=50,
					command=self.applyLate,
					background='LightYellow',width=4)
		self.cpd_spin.grid(row=row, column=1, sticky=W)
		tkExtra.Balloon.set(self.cpd_spin,
			'Colors per decay for log scales')

		# set trace sentinels
		self.cbmin.trace(   'w', self.cpdChange)
		self.cblog.trace(   'w', self.cpdChange)
		self.cpd.trace(     'w', self.cpdChange)
		self.cbsmooth.trace('w', self.cpdChange)
		self.cbcolors.trace('w', self.cpdChange)

		pageframe.grid_columnconfigure(1, weight=1)

		return pageframe

	# --------------------------------------------------------------------
	def cpdChange(self, a, b, c):
		if self.cbsmooth.get():
			self.cbspin.config(state=DISABLED)
			self.cpd_spin.config(state=DISABLED)
			return
		else:
			self.cbspin.config(state=NORMAL)

		if int(self.cblog.get()):
			self.cpd_spin.config(state=NORMAL)
		else:
			self.cpd_spin.config(state=DISABLED)
			return

		try:
			cpd = float(self.cpd.get())
		except ValueError:
			cpd = 0.0

		if cpd < 1.0: return

		# Calculate max
		try: cbmin = float(self.cbmin.get())
		except: return
		try: colors = float(self.cbcolors.get())
		except: return

		try:
			cbmax = cbmin * math.pow(10.0, colors/cpd)
			self.cbmax.set('%g'%(cbmax))
		except (ValueError, OverflowError):
			pass

	# --------------------------------------------------------------------
	def paletteInit(self):
		if self.activepage is None or \
		   not self.activepage.startswith('Palette'): return

		oldInit = self._init
		self._init = True
		self.cbpalette.set(option(self.layer,      self.entry('palette')))
		self.cblog.set(option(self.layer,          self.entry('log')))
		self.cbmin.set(option(self.layer,          self.entry('min')))
		self.cbmax.set(option(self.layer,          self.entry('max')))
		self.cbsmooth.set(option(self.layer,       self.entry('smooth')))
		self.cbcolors.set(option(self.layer,       self.entry('n')))
		self.cpd.set(option(self.layer,            self.entry('cpd')))
		self.cbinv.set(option(self.layer,          self.entry('inv')))
		self.cbalphamin.set(option(self.layer,     self.entry('alphamin')))
		self.cbalphamax.set(option(self.layer,     self.entry('alphamax')))

		self.cblabel.delete(0,END)
		self.cblabel.insert(0, option(self.layer, self.entry('label')))

		self._init = oldInit

	# --------------------------------------------------------------------
	def paletteUpdate(self):
		if not self.activepage.startswith('Palette'): return
#		if not option(self.layer, self.activepage): return

		self.layer[self.entry('palette')]     = self.cbpalette.get()
		self.layer[self.entry('log')]         = self.cblog.get()
		self.layer[self.entry('inv')]         = self.cbinv.get()
		self.layer[self.entry('alphamin')]    = self.cbalphamin.get()
		self.layer[self.entry('alphamax')]    = self.cbalphamax.get()
		self.layer[self.entry('smooth')]      = self.cbsmooth.get()
		try:    self.layer[self.entry('min')] = float(self.cbmin.get())
		except: self.layer[self.entry('min')] = 1E-5
		try:    self.layer[self.entry('max')] = float(self.cbmax.get())
		except: self.layer[self.entry('max')] = 1E5
		try:    self.layer[self.entry('n')]   = int(self.cbcolors.get())
		except: self.layer[self.entry('n')]   = 30
		try:    self.layer[self.entry('cpd')] = int(self.cpd.get())
		except: self.layer[self.entry('cpd')] = 3
		self.layer[self.entry('label')]       = self.cblabel.get()

	# --------------------------------------------------------------------
	# Empty frame
	# --------------------------------------------------------------------
	def emptyFrame(self):
		pageframe = Frame(self.mainframe)

		Label(pageframe, text='Not implemented yet').pack(expand=YES,fill=BOTH)

		return pageframe

	# --------------------------------------------------------------------
	# Initialize and set defaults options
	# --------------------------------------------------------------------
	def init(self):
		self._init = True	# Avoid updating

		# Layers
		for i in range(1, self.listbox.size()):
			page, sel = self.listbox.get(i)
			if sel=='': continue
			s = self.listbox.selection_includes(i)
			if option(self.layer, page):
				sel = _ON
			else:
				sel = _OFF
			self.listbox.set(i, (page,sel))
			if s: self.listbox.selection_set(i)

		# Global layer
		self.globalLayer.set(self.layer._global)

		# Options
		self.axes.set(option(self.layer,         "Options.axes"))
		self.titleVar.set(option(self.layer,     "Options.title"))
		self.showErrors.set(option(self.layer,   "Options.errors"))
		self.gridVar.set(option(self.layer,      "Options.grid"))
		self.gridlevel.set(option(self.layer,    "Options.gridlevel"))
		self.latticeLevel.set(option(self.layer, "Options.latticelevel"))
		self.crosshair.set(option(self.layer,    "Options.crosshair"))
		self.background.set(option(self.layer,   "Options.textbackground"))
		self.viewport.set(option(self.layer,     "Options.viewport"))
		self.shortCut.set(option(self.layer,     "Options.shortcut"))
		self.generalFont.set(option(self.layer,  "Options.font"))
		self.gridFont.set(option(self.layer,     "Options.gridfont"))
		self.paletteFont.set(option(self.layer,  "Options.palettefont"))

		# Show
		self.showvertex.set(option(self.layer,  "Show.vertex"))
		self.labels.set(option(self.layer,      "Show.label"))
		self.colorvalue.set(option(self.layer,  "Show.color"))
		self.showlattice.set(option(self.layer, "Show.lattice"))
		self.showvoxel.set(option(self.layer,   "Show.voxel"))
		self.showpalette.set(option(self.layer, "Show.palette"))

		# Background image
		self.setImagefile(option(self.layer,      "Image.file"))
		self.imageAlpha.set(option(self.layer,    "Image.alpha"))
		self.setImageBlack(int(option(self.layer, "Image.black")))
		self.setImageWhite(int(option(self.layer, "Image.white")))
		self.imagePrompt.set(option(self.layer,   "Image.prompt"))

		# Userdump
		self.setUserdumpFile(option(self.layer,   "Userdump.file"))
		self.userdumpStart.set(option(self.layer, "Userdump.start"))
		self.userdumpEvents.set(option(self.layer,"Userdump.n"))
		for i,name in enumerate(self.userdumpList.lists[1].get(0,END)):
			on    = self.layer.get("Userdump.%s"%(name), 0) and _ON or _OFF
			emin  = self.layer.get("Userdump.%s.emin"%(name),  "")
			emax  = self.layer.get("Userdump.%s.emax"%(name),  "")
			color = self.layer.get("Userdump.%s.color"%(name), "")
			self.userdumpList.set(i, (on,name,emin,emax,color))
			if color!="":
				self.userdumpList.itemconfig(i, background=color)

		# 3D
		self.projection.set(_PROJECTIONS[option(self.layer,'3D.projection')])
		self.ambient.set(option(self.layer,        '3D.ambient'))
		self.antialias.set(option(self.layer,      '3D.antialias'))
		self.clipbody.set(option(self.layer,       '3D.clip.0'))
		self.clipbody1.set(option(self.layer,      '3D.clip.1'))
		self.clipbody2.set(option(self.layer,      '3D.clip.2'))
		self.clipnegative.set(option(self.layer ,  '3D.clipnegative.0'))
		self.clipnegative1.set(option(self.layer , '3D.clipnegative.1'))
		self.clipnegative2.set(option(self.layer , '3D.clipnegative.2'))
		self.deflights.set(option(self.layer ,     '3D.deflights'))
		self.edgedetect.set(option(self.layer ,    '3D.edgedetect'))
		self.fov.set(option(self.layer,            '3D.fov'))
		self.projectbody.set(option(self.layer,    '3D.project.0'))
		self.projectbody1.set(option(self.layer,   '3D.project.1'))
		self.projectbody2.set(option(self.layer,   '3D.project.2'))
		self.reflections.set(option(self.layer,    '3D.reflections'))
		self.shadows.set(option(self.layer,        '3D.shadows'))
		self.skip1stblack.set(option(self.layer ,  '3D.skip1stblack'))
		self.usrbinastexture.set(option(self.layer,'3D.usrbinastexture'))
		self.xray.set(option(self.layer,           '3D.xray'))

		# Voxel
		self.roiAlpha.set(option(self.layer,       'Voxel.roialpha'))
		self.voxelRefresh()

		# Multi init Usrbin/Palette
		self._usrbinFile = ''
		self.usrbinInit()
		self.paletteInit()

		self._init = False

	# --------------------------------------------------------------------
	# Get layers information and save it in project
	# --------------------------------------------------------------------
	def get(self):
		if self.layer is None: return

		# Get the active layers
		for page, sel in self.listbox.get(0,END):
			if sel!='': self.layer[page] = (sel[1]=='X')

		# Update layer
		oldname = self.layer.name
		self.layer.name = self.layerBox.get()
		if oldname != self.layer.name:
			if self.layer._global:
				try: del Project._globalLayers[oldname]
				except KeyError: pass
			try: del self.project.geometry.layers[oldname]
			except KeyError: pass
			self.project.geometry.layers[self.layer.name] = self.layer

		# Global layer or not
		self.layer._global = self.globalLayer.get()

		# Options
		self.layer["Options.axes"]          = self.axes.get()
		self.layer["Options.title"]         = self.titleVar.get()
		self.layer["Options.errors"]        = self.showErrors.get()
		self.layer["Options.grid"]          = self.gridVar.get()
		self.layer["Options.gridlevel"]     = self.gridlevel.get()
		self.layer["Options.latticelevel"]  = self.latticeLevel.get()
		self.layer["Options.crosshair"]     = self.crosshair.get()
		self.layer["Options.textbackground"]= self.background.get()
		self.layer["Options.viewport"]      = self.viewport.get()
		self.layer["Options.font"]          = self.generalFont.get()
		self.layer["Options.shortcut"]      = self.shortCut.get()
		self.layer["Options.gridfont"]      = self.gridFont.get()
		self.layer["Options.palettefont"]   = self.paletteFont.get()

		# Show
		self.layer["Show.vertex"]  = self.showvertex.get()
		self.layer["Show.label"]   = self.labels.get()
		self.layer["Show.color"]   = self.colorvalue.get()
		self.layer["Show.lattice"] = self.showlattice.get()
		self.layer["Show.voxel"]   = self.showvoxel.get()
		self.layer["Show.palette"] = self.showpalette.get()

		# Background image
		try:
			self.layer['Image.file']   = self.project.relativePath(self.imageFilename)
		except OSError:
			pass
		self.layer["Image.alpha"]  = self.imageAlpha.get()
		self.layer["Image.black"]  = int(self.blackImageColor["background"][1:],16)
		self.layer["Image.white"]  = int(self.whiteImageColor["background"][1:],16)
		self.layer["Image.prompt"] = self.imagePrompt.get()

		# Userdump
		try:
			self.layer["Userdump.file"] = self.project.relativePath(self.userdumpFilename)
		except OSError:
			pass
		try:    self.layer["Userdump.start"] = int(self.userdumpStart.get())
		except: self.layer["Userdump.start"] = 1
		try:    self.layer["Userdump.n"] = int(self.userdumpEvents.get())
		except: self.layer["Userdump.n"] = 100
		for on, name, emin, emax, color in self.userdumpList.get(0,END):
			self.layer["Userdump.%s"%(name)] = int(on==_ON)
			if emin != "":
				self.layer["Userdump.%s.emin"%(name)] = emin
			else:
				try: del self.layer["Userdump.%s.emin"%(name)]
				except KeyError: pass

			if emax != "":
				self.layer["Userdump.%s.emax"%(name)] = emax
			else:
				try: del self.layer["Userdump.%s.emax"%(name)]
				except KeyError: pass

			if color != "":
				self.layer["Userdump.%s.color"%(name)] = color
			else:
				try: del self.layer["Userdump.%s.color"%(name)]
				except KeyError: pass

		# 3D
		try: self.layer["3D.projection"]  = _PROJECTIONS.index(self.projection.get())
		except: self.layer["3D.projection"] = 0
		self.layer["3D.deflights"]        = self.deflights.get()
		try:    self.layer["3D.fov"]      = float(self.fov.get())
		except: self.layer["3D.fov"]      = 60.0
		self.layer["3D.ambient"]          = self.ambient.get()
		self.layer["3D.antialias"]        = self.antialias.get()
		self.layer["3D.clip.0"]           = self.clipbody.get()
		self.layer["3D.clip.1"]           = self.clipbody1.get()
		self.layer["3D.clip.2"]           = self.clipbody2.get()
		self.layer["3D.clipnegative.0"]   = self.clipnegative.get()
		self.layer["3D.clipnegative.1"]   = self.clipnegative1.get()
		self.layer["3D.clipnegative.2"]   = self.clipnegative2.get()
		self.layer["3D.edgedetect"]       = self.edgedetect.get()
		self.layer["3D.project.0"]        = self.projectbody.get()
		self.layer["3D.project.1"]        = self.projectbody1.get()
		self.layer["3D.project.2"]        = self.projectbody2.get()
		self.layer["3D.reflections"]      = self.reflections.get()
		self.layer["3D.shadows"]          = self.shadows.get()
		self.layer["3D.skip1stblack"]     = self.skip1stblack.get()
		self.layer["3D.usrbinastexture"]  = self.usrbinastexture.get()
		self.layer["3D.xray"]             = self.xray.get()

		# Voxel
		self.layer["Voxel.roialpha"]      = self.roiAlpha.get()
		self.voxelUpdate()

		# Multi updates Usrbin/Palette
		self.usrbinUpdate()
		self.paletteUpdate()

		# Correct layer
		self.removeSameWithDefaults()
		self.checkSystem()

	# --------------------------------------------------------------------
	# Remove layer options that are the same with the default ones
	# --------------------------------------------------------------------
	def removeSameWithDefaults(self):
		remove = []
		for name,value in list(self.layer.var.items()):
			if name == '_system': continue
			if "._visible" in name: continue
			if value == DEFLAYERS.get(name):
				remove.append(name)
		for name in remove:
			del self.layer.var[name]

	# --------------------------------------------------------------------
	# Check if system was changed and save it to local
	# --------------------------------------------------------------------
	def checkSystem(self):
		if not self.layer.var.get("_system",False): return

		# Find in the global
		default = Project._systemLayers.get(self.layer.name)
		if default is None:
			self.makeLocal()
			return
		else:
			# Compare with the default
			for name,value in list(self.layer.var.items()):
				# skip options
				if name.startswith("Options"): continue
				defvalue = option(default,name)
				if defvalue is not None and value != defvalue:
					self.makeLocal()
					return

	# --------------------------------------------------------------------
	# Make a system layer local
	# --------------------------------------------------------------------
	def makeLocal(self):
		# clone layer to project
		self.layer = self.layer.clone()
		try:	# remove the system flag
			del self.layer.var["_system"]
			del self.layer.var["_order"]
		except KeyError:
			pass
		self.layer._global = False

		self.project.geometry.layers[self.layer.name] = self.layer
		self.globalLayer.set(False)
		self.globalButton.config(state=NORMAL)

		# restore default layer
		try:
			Project._globalLayers[self.layer.name] = \
				Project._systemLayers[self.layer.name].clone()
		except KeyError:
			pass

	# --------------------------------------------------------------------
	def addLayer(self, event=None):
		self.flair.setModified()
		self.get()
		self.layer = None

		# find unique name
		i = 1
		while True:
			name = 'Layer %02d'%(i)
			if name not in self.project.geometry.layers: break
			i += 1
		self.project.geometry.layers[name] = Project.GeometryLayerInfo(name)
		self.fillLayers(name)
		self.layerBox._text.select_range(0,END)
		self.layerBox._text.focus_set()

	# --------------------------------------------------------------------
	def cloneLayer(self, event=None):
		self.flair.setModified()
		self.get()

		layer = self.layer.clone()
		try:	# remove the system flag
			del layer.var["_system"]
		except KeyError:
			pass
		layer._global = False

		# Find unique name
		pat = Project._LAST_NUM.match(layer.name)
		if pat:
			name = pat.group(1)
			n    = int(pat.group(2))+1
		else:
			name = layer.name
			n    = 1
		while True:
			guess = '%s %02d'%(name,n)
			if guess not in self.project.geometry.layers:
				layer.name = guess
				break
			n += 1

		self.project.geometry.layers[guess] = layer
		self.fillLayers(layer.name)
		self.layerBox._text.select_range(0,END)
		self.layerBox._text.focus_set()

	# --------------------------------------------------------------------
	def delLayer(self, event=None):
		if self.layer.get('_system'): return

		if not tkFlair.askyesno('Delete Layer',
				'Do you want to delete the active layer \'%s\''\
				%(self.layer.name),
				parent=self.listbox):
			return

		lst = self.project.geometry.layersList()

		if self.layer._global:
			try: del Project._globalLayers[self.layer.name]
			except KeyError: pass

		try: del self.project.geometry.layers[self.layer.name]
		except KeyError: pass
		tkFlair.updateLayers()

		try:
			idx = lst.index(self.layer.name)
			del lst[idx]
			while idx>0 and idx >= len(lst):
				idx -= 1
			name = lst[idx]
			self.layer = self.project.geometry.layers[name]
			self.fillLayers(name)
		except:
			self.layer = None
			self.addLayer()

	# --------------------------------------------------------------------
	def fillLayers(self, layername=None):
		layers = self.project.geometry.layersList()
		self.layerBox.fill(layers)
		if layername is not None and layername in self.project.geometry.layers:
			self.layerBox.set(layername)

	# --------------------------------------------------------------------
	def changeLayer(self, event=None):
		self.setLayer(self.layerBox.get())

	# --------------------------------------------------------------------
	def setLayer(self, layer):
		self.layerBox._text.delete(0, END)
		self.layerBox._text.insert(0, layer)
		self.layer = self.project.geometry.layers[layer]
		self.populateListbox()
		self.showPage(self.activepage)
		self.init()

	# --------------------------------------------------------------------
	def help(self, event=None):
		Manual.show('F4.14')

#-------------------------------------------------------------------------------
if __name__ == '__main__':
	import Project

	if len(sys.argv) == 1: sys.exit(0)

	Input.init()

	filename = sys.argv[1]
	fn,ext = os.path.splitext(filename)

	project = Project.Project()
	if ext == '.flair':
		project.load(filename)
	else:
		project.loadInput(filename)
	project.input.convert2Names()

	palette = Palette.Palette()
	palette.load(project)

	tk = Tk()
	tkFlair.loadIcons()
	tkFlair.openIni()
	tkFlair.readIni()

	#layers = {}
	cd = GeometryLayers(tk, project, None)
