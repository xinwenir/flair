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
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	16-Aug-2012

import os
from tkinter import *
import math
import time

try:
	import pydicom as dicom
except ImportError:
	import dicom

import tkFlair
import tkExtra
import tkDialogs
import bFileDialog

import undo
import Dicom
import Input
import Ribbon
import FlairRibbon
from Layout import _ASSIGNMA_field

NEW = "<new>"

#-------------------------------------------------------------------------------
class UnitListAssist(tkExtra.SortAssist):
	def __call__(self, x):
		try:
			return float(x[self.column])
		except:
			return x[self.column]

#===============================================================================
# Voxel Tab
#===============================================================================
class VoxelTab(FlairRibbon.FlairTab):
	#-----------------------------------------------------------------------
	def __init__(self, master, page, **kw):
		FlairRibbon.FlairTab.__init__(self, master, page, **kw)

		# -------------- Materials ---------------
		row,col = 0,0
		Label(self, text="Materials:").grid(
				row=row, column=col, pady=1, sticky=E)
		col += 1
		self.materials = Label(self,
				foreground=tkFlair._ELABEL_FOREGROUND_COLOR,
				background=tkFlair._ELABEL_BACKGROUND_COLOR,
				relief=SUNKEN, anchor=W)
		self.materials.grid(row=row, column=col, sticky=EW)
		self.materials.bind('<Button-1>', self.loadMaterials)
		tkExtra.Balloon.set(self.materials, "Materials FLUKA input file to embed to voxel")

		col += 1
		b = Button(self, image=tkFlair.icons["load"], pady=0,
				command=self.loadMaterials)
		b.grid(row=row, column=col)
		tkExtra.Balloon.set(b, "Load materials FLUKA input file")

		col += 1
		b = Button(self, image=tkFlair.icons["edit"], pady=0,
				state = DISABLED,
				command=self.editMaterials)
		b.grid(row=row, column=col)
		tkExtra.Balloon.set(b, "Edit materials")

		col += 1
		b = Button(self, image=tkFlair.icons["x"], pady=0,
				command=self.cancelMaterials)
		b.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(b, "Cancel materials")

		# ---
		col += 1
		Label(self, text="# materials:").grid(row=row, column=col, sticky=E)
		col += 1
		self.nmaterials = Label(self, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.nmaterials.grid(row=row, column=col, sticky=EW)

		# ---
		col += 1
		Label(self, text="Field:").grid(row=row, column=col, sticky=E)
		col += 1
		self.field = tkExtra.Combobox(self, width=30,
					foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
					background="White")
		self.field.grid(row=row, column=col, columnspan=2, sticky=EW)
		self.field.fill([x[0] for x in _ASSIGNMA_field])

		# -------------- Unit to Material ---------------
		row += 1
		col   = 0
		Label(self, text="Unit to Material :").grid(
				row=row, column=col, pady=1, sticky=E)
		col += 1
		self.unit2mat = Label(self,
				foreground=tkFlair._ELABEL_FOREGROUND_COLOR,
				background=tkFlair._ELABEL_BACKGROUND_COLOR,
				relief=SUNKEN, anchor=W)
		self.unit2mat.grid(row=row, column=col, sticky=EW)
		self.unit2mat.bind('<Button-1>', self.loadUnit2Material)
		tkExtra.Balloon.set(self.unit2mat, "Unit to material conversion")

		col += 1
		b = Button(self, image=tkFlair.icons["load"], pady=0,
				command=self.loadUnit2Material)
		b.grid(row=row, column=col)
		tkExtra.Balloon.set(b, "Load unit to material conversion")

		col += 1
		b = Button(self, image=tkFlair.icons["edit"], pady=0,
				command=self.editUnit2Material)
		b.grid(row=row, column=col)
		tkExtra.Balloon.set(b, "Edit hounsfield units")

		col += 1
		b = Button(self, image=tkFlair.icons["save"], pady=0,
				command=self.saveUnit2Material)
		b.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(b, "Save unit to material conversion")

		# ---
		col += 1
		Label(self, text="RTSTRUCT:").grid(row=row, column=col, sticky=E)
		col += 1
		self.rtstruct = tkExtra.Combobox(self, width=30, foreground="DarkBlue")
		self.rtstruct.grid(row=row, column=col, columnspan=4, sticky=EW)

		# --- Table editing ---
		row += 1
		col  = 0
		self.unit2matList = tkExtra.MultiListbox(self,
				(#('< Unit',        5, None),
				 ('HUmin',          5, None),
				 ('HUmax',          5, None),
				 ('Material',      10, None),
				 ('Crho_min',      10, None),
				 ('Crho_max',      10, None),
				 ('CdEdx_rel_min', 10, None),
				 ('CdEdx_rel_max', 10, None)),
				 background="White")
		self.unit2matList.sortAssist = UnitListAssist
		self.unit2matList.grid(row=row, column=0, columnspan=9, rowspan=4, sticky=NSEW)
		self.unit2matList.bindList("<ButtonRelease-1>", self.unit2matRelease)
		self.unit2matList.bindList("<Return>",          self.unit2matEdit)
		self.unit2matList.bindList("<KP_Enter>",        self.unit2matEdit)
		self.unit2matList.bindList("<F2>",              self.unit2matEdit)
		self.unit2matList.bindList("<Insert>",          self.add)
		self.unit2matList.bindList("<Control-Return>",  self.add)
		self.unit2matList.bindList("<Delete>",          self.delete)
		self.unit2matList.bindList("<Control-Key-D>",   self.clone)
		self.unit2matList.bindLeftRight()
		self.unit2matList.sortAssist = None
		self.unit2matModified = False

		# ---
		col += 9
		b = Button(self, image=tkFlair.icons["add"],
				takefocus=False, command=self.add)
		tkExtra.Balloon.set(b, "Add unit range")
		b.grid(row=row, column=col, sticky=EW)

		row += 1
		b = Button(self, image=tkFlair.icons["del"],
				takefocus=False, command=self.delete)
		tkExtra.Balloon.set(b, "Delete unit range")
		b.grid(row=row, column=col, sticky=EW)

		# ---
		row += 1
		b = Button(self, image=tkFlair.icons["clone"],
				takefocus=False, command=self.clone)
		tkExtra.Balloon.set(b, "Clone unit range")
		b.grid(row=row, column=col, sticky=EW)

		# Stretching
		self.grid_rowconfigure(row+1, weight=1)
		self.grid_columnconfigure(1, weight=1)
		self.grid_columnconfigure(6, weight=3)
		self.grid_columnconfigure(8, weight=1)

		self.lastPath = ""

	#-----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairTab.createRibbon(self)

		# ========== Files ===========
		group = Ribbon.LabelGroup(self.ribbon, "Execute")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# -------------- Output voxel file ---------------
		col, row = 0,0
		self.voxelfile = Label(group.frame,
				foreground=tkFlair._ELABEL_FOREGROUND_COLOR,
				background=tkFlair._ELABEL_BACKGROUND_COLOR,
				width=30,
				relief=SUNKEN,
				anchor=W)
		self.voxelfile.grid(row=row, column=col, columnspan=4, sticky=EW)
		self.voxelfile.bind('<Button-1>', self.saveVoxelFile)
		tkExtra.Balloon.set(self.voxelfile, "Output voxel file")

		col += 4
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["dicom"],
				text="VOXEL",
				compound=LEFT,
				command=self.makeVoxel,
				foreground="Darkred",
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(b, "Process DICOM files and generate VOXEL")

		col += 1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["input"],
				text="Card",
				compound=LEFT,
				command=self.voxelCard,
				background=Ribbon._BACKGROUND)

		b.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(b, "Add or modify VOXEL card")

		# -------------- Output USRBIN file ---------------
		col, row = 0, 1
		self.usrbinfile = Label(group.frame,
				foreground=tkFlair._ELABEL_FOREGROUND_COLOR,
				background=tkFlair._ELABEL_BACKGROUND_COLOR,
				width=30,
				relief=SUNKEN, anchor=W)
		self.usrbinfile.grid(row=row, column=col, columnspan=4, sticky=EW)
		self.usrbinfile.bind('<Button-1>', self.saveUsrbinFile)
		tkExtra.Balloon.set(self.usrbinfile, "Output Usrbin file")

		col += 4
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["USRBIN"],
				text="USRBIN",
				compound=LEFT,
				command=self.makeUSRBIN,
				foreground="Darkred",
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(b, "Process DICOM files and generate USRBIN file")

		col += 1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["input"],
				text="Card",
				compound=LEFT,
				command=self.usrbinCard,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(b, "Add or modify USRBIN card")

		# ---
		col, row = 0, 2
		Label(group.frame, text = "Gantry:",
			background=Ribbon._BACKGROUND).grid(
				row=row, column=col, sticky=E)

		col += 1
		self.gantry = tkExtra.FloatEntry(group.frame, width=5, background="White")
		self.gantry.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.gantry, "Gantry rotation angle in degrees")

		# ---
		col += 1
		Label(group.frame, text = "Patient:",
			background=Ribbon._BACKGROUND).grid(
				row=row, column=col, sticky=E)
		col += 1
		self.patient = tkExtra.FloatEntry(group.frame, width=5, background="White")
		self.patient.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.patient, "Patient rotation angle in degrees")

		# ---
		col += 2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["input"],
				text="ROT-DEFI",
				compound=LEFT,
				command=self.rotdefiCard,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(b, "Add ROT-DEFI cards in input")

		self._updated = False

	# ----------------------------------------------------------------------
	def dicom(self):	return self.page.dicom

	# ----------------------------------------------------------------------
	def activate(self):
		if self.page.loadDicom():
			self.flair.notify("No DICOM selected",
				"Please select a DICOM set first")
		else:
			self.refresh()

	# ----------------------------------------------------------------------
	def deactivate(self):
		self.get()

	# ----------------------------------------------------------------------
	# Refresh values
	# ----------------------------------------------------------------------
	def refresh(self):
		if self.dicom() is None: return
		di = self.page.dicominfo
		if di is None: return

		v = di.get("materials","")
		self._updated = True
		if v != self.materials["text"]:
			self.materials["text"] = v
			self.updateMaterials()

		v = di.get("unit2mat","")
		if v != self.unit2mat["text"]:
			self.unit2mat["text"] = v
			self.dicom().readHounsfield(v)
			self.updateUnit2Material()

		self.voxelfile["text"]  = di.get("voxel","")
		self.usrbinfile["text"] = di.get("usrbin","")

		self.field.set(_ASSIGNMA_field[int(di.get("field",0))][0])

		self.gantry.set( di.get("gantry",""))
		self.patient.set(di.get("patient",""))

		self.rtstruct.set(di.get("rtstruct","").strip())

		# Fill a list of RTSTRUCT
		self.rtstruct.delete(0,END)
		self.rtstruct.insert(END,"")
		for di in self.project.dicoms:
			try:
				if di.modality=="RTSTRUCT": self.rtstruct.insert(END,di.name)
			except AttributeError:
				pass

	# ----------------------------------------------------------------------
	def get(self):
		if self.dicom() is None: return
		if not self._updated: return
		di = self.page.dicominfo
		if di is None: return

		di["materials"] = self.materials["text"]
		di["unit2mat"]  = self.unit2mat["text"]
		di["voxel"]     = self.voxelfile["text"]
		di["usrbin"]    = self.usrbinfile["text"]
		di["rtstruct"]  = self.rtstruct.get()
		di["gantry"]    = self.gantry.get()
		di["patient"]   = self.patient.get()
		di["field"]     = self.fieldIndex()

	# ----------------------------------------------------------------------
	def fieldIndex(self):
		mag = self.field.get()
		for n,i in _ASSIGNMA_field:
			if n==mag:
				if i=="": i=0
				return i
		return 0

	# ----------------------------------------------------------------------
	# Load material definitions
	# ----------------------------------------------------------------------
	def loadMaterials(self, event=None):
		fn = bFileDialog.askopenfilename(master=self,
			title="Import Materials from FLUKA file",
			filetypes=[("Fluka Input files",("*.inp","*.fluka")),
				("All","*")],
			initialdir=self.lastPath,
			initialfile=self.materials["text"])

		if len(fn) == 0: return

		self.lastPath = os.path.abspath(os.path.dirname(fn))
		relfn = self.project.relativePath(fn)
		self.materials["text"] = relfn
		self.updateMaterials()

	# ----------------------------------------------------------------------
	def editMaterials(self, event=None):
		pass
#		import flair
#		flair.newWindow()

	# ----------------------------------------------------------------------
	def cancelMaterials(self, event=None):
		self.materials["text"] = ""
		self.updateMaterials()

	# ----------------------------------------------------------------------
	def updateMaterials(self):
		if self.dicom() is None: return
		if self.materials["text"] == "":
			self.dicom().material = None
			self.nmaterials["text"] = "none"
		else:
			self.dicom().material = Input.Input(self.materials["text"])
			n = len(self.dicom().material["MATERIAL"])
			if n>0:
				self.nmaterials["text"] = str(n)
			else:
				self.nmaterials["text"] = "none"

	# ----------------------------------------------------------------------
	def loadUnit2Material(self, event=None):
		if self.dicom() is None: return
		if self.unit2matModified:
			if not tkFlair.askyesno("Load Unit 2 Material",
					"Current Unit to material list is modified\nDiscard changes",
					parent=self):
				return
		fn = bFileDialog.askopenfilename(master=self,
			title="Unit to material conversion",
			filetypes=[("Unit to material data files","*.mat"),
				("Data files","*.dat"),
				("Text files","*.txt"),
				("All","*")],
			initialdir=self.lastPath,
			initialfile=self.unit2mat["text"])

		if len(fn) == 0: return False

		relfn = self.project.relativePath(fn)
		self.lastPath = os.path.abspath(os.path.dirname(relfn))
		self.unit2mat["text"] = relfn
		self.dicom().readHounsfield(relfn)
		self.updateUnit2Material()
		self.unit2matModified = False
		return True

	# ----------------------------------------------------------------------
	def saveUnit2Material(self, event=None):
		if self.dicom() is None: return
		#if not self.unit2matModified: return
		self.getHounsfieldFromList()
		fn = self.unit2mat["text"]
		try:
			f = open(fn,"w")
		except:
			self.flair.notify("Cannot open file",
				"Cannot open file '%s' for writing hounsfield units"%(fn),
				tkFlair.NOTIFY_ERROR)
			return

		for hu in self.dicom().hounsfield:
			f.write("%s\n"%(hu.str()))
		f.close()

		self.unit2matModified = False
		self.flair.notify("Unit to Materials file written",
			"File %s with units to materials updated"%(fn),
			tkFlair.NOTIFY_NORMAL)

	# ----------------------------------------------------------------------
	def editUnit2Material(self, event=None):
		tkFlair.editor(self.unit2mat["text"])

	# ----------------------------------------------------------------------
	def updateUnit2Material(self):
		if self.dicom() is None: return
		self.unit2matList.delete(0,END)
		humin = -3072
		for hu in self.dicom().hounsfield:
			self.unit2matList.insert(END,
				(humin, hu.humax-1, hu.mat,
				 hu.crho_min, hu.crho_max,
				 hu.cdedx_min, hu.cdedx_max))
			humin = hu.humax
		self.unit2matList.sort(0,False)

	# ----------------------------------------------------------------------
	def saveVoxelFile(self, event=None):
		fn = bFileDialog.asksaveasfilename(master=self,
			title="Voxel file",
			filetypes=[("Voxel files","*.vxl"),
				("All","*")],
			initialdir=self.lastPath,
			initialfile=self.voxelfile["text"])

		if len(fn) == 0: return False

		self.lastPath = os.path.abspath(os.path.dirname(fn))
		relfn = self.project.relativePath(fn)
		self.voxelfile["text"] = relfn
		return True

	# ----------------------------------------------------------------------
	# Create the voxel file
	# ----------------------------------------------------------------------
	def makeVoxel(self):
		if self.dicom() is None or len(self.dicom()._slices) == 0:
			self.flair.notify("No DICOM processed",
				"Please select a DICOM dataset",
				tkFlair.NOTIFY_WARNING)
			return
		self.get()

		# Do not complain if materials are not there!
		#if len(self.dicom().material["MATERIAL"])==0:
		#	self.flair.notify("No Material file",
		#		"Please select a FLUKA file containing needed materials",
		#		tkFlair.NOTIFY_WARNING)
		#	return

		if len(self.dicom().hounsfield)==0:
			self.flair.notify("No Unit to Material file",
				"Please select or create a 'Unit to Material' conversion file",
				tkFlair.NOTIFY_WARNING)
			return

		if len(self.voxelfile["text"])==0:
			if not self.saveVoxelFile():
				return

		# time keeping
		t0 = time.time()

		# Interpolate correction factors
		if self.unit2matModified:
			self.getHounsfieldFromList()
			self.dicom().interpolateCorrections()

		# Setup the assignmats
		self.dicom().assignMaterials2Regions(self.fieldIndex())

		# Do the real processing
		if not self.dicom()._writeVoxelInit(self.voxelfile["text"]):
			self.flair.notify("DICOM to Voxel",
				"ERROR: Cannot open output Voxel file",
				tkFlair.NOTIFY_ERROR)
			return

		dlg = tkDialogs.ProgressDialog(self, "Creating VOXEL file")
		if self.dicom().slice:
			# Multiple slices
			dlg.setLimits(0, len(self.dicom()._slices))
			for i,(l,fn) in enumerate(self.dicom()._slices):
				# Show progress
				if dlg.show(i, "Processing: "+os.path.basename(fn)): return
				self.dicom()._writeVoxelSlice(fn,i)
		else:
			# Unique 3D data
			fn = self.dicom()._slices[0][1]
			dlg.setLimits(0.0, self.dicom().nz)
			# Show progress
			dlg.show(0, "Processing: "+os.path.basename(fn))
			def update(s): return dlg.show(s)
			self.dicom()._writeVoxelSlice(fn, update=update)
		dlg.stop()

		# time keeping
		t1 = time.time()
		self.flair.log("CT processed in %.1f s"%(t1-t0))

		# Processing RTSTRUCT
		if self._processRTSTRUCT():
			t2 = time.time()
			self.flair.log("RTSTRUCT processed in %.1f s"%(t2-t1))

		# Finalize voxel by writing data
		if not self.dicom()._writeVoxelEnd():
			self.flair.notify("DICOM to Voxel",
				"Error writing Voxel file",
				tkFlair.NOTIFY_ERROR)
			return

		# time keeping
		t3 = time.time()
		self.flair.log("Voxel generated in %.1f s"%(t3-t0))

		self.flair.notify("DICOM to Voxel",
			"Voxel file generated successfully from DICOM in %.1fs"%(t3-t0))

	# ----------------------------------------------------------------------
	def _processRTSTRUCT(self):
		# Processing RTSTRUCT
		rtuid = self.rtstruct.get()
		if rtuid == "": return False

		# Find and open RTSTRUCT dataset
		for di in self.project.dicoms:
			if di["uid"] == rtuid:
				try:
					dataset = dicom.read_file(
							os.path.join(di["directory"], di.files[0][0]),
							force=True)
				except:
					return False
				break
		else:
			return False

		# For all(cropped) SOP of DICOM scan the RTSTRUCT
		dlg = tkDialogs.ProgressDialog(self, "Processing RTSTRUCTs")
		dlg.setLimits(0, len(self.page.dicominfo.files))
		for i,(fn,z,sop) in enumerate(self.page.dicominfo.files):
			#if i<70 or i>80: continue
			self.dicom()._roiStructSlice(dataset,sop)
			if dlg.show(i, "Processing: "+sop): return
		dlg.stop()
		del dataset
		return True

	# ----------------------------------------------------------------------
	# Create or modify the input VOXEL card
	# ----------------------------------------------------------------------
	def voxelCard(self, event=None):
		if self.dicom() is None: return
		self.page.get()

		undoinfo = self.flair.undolistRefresh()

		fn,ext = os.path.splitext(os.path.basename(self.voxelfile["text"]))
		whats = [fn,
			self.dicom().xmin()/10.0,
			self.dicom().ymin()/10.0,
			self.dicom().zmin()/10.0]

		inp = self.flair.project.input

		card = inp["VOXELS"]
		if card:
			card = card[0]
			msg  = "Modify VOXELS card"

			undoinfo.append(self.flair.saveCardUndoInfo(card.pos()))
			for i,w in enumerate(whats):
				card.setWhat(i,w)
		else:
			if not inp.cardlist: inp.minimumInput()
			geobegin = inp["GEOBEGIN"][0]
			pos = geobegin.pos() + 1

			# add card to input
			card   = Input.Card("VOXELS",whats)
			msg    = "Add VOXELS card"
			undoinfo.append(self.flair.addCardUndo(card, pos, False))

		# Load the new voxel file
		card.loadVoxel()

		# Finalize
		self.flair.addUndolistRefresh(undoinfo,msg)
		self.flair.refresh("card")

		self.flair.notify(msg,
			"A VOXELS card corresponding to the dicom file " \
			" has been added/modified to your input")

	# ----------------------------------------------------------------------
	def saveUsrbinFile(self, event=None):
		fn = bFileDialog.asksaveasfilename(master=self,
			title="Usrbin file",
			filetypes=[("Usrbin files",("*.bnn","*usrbin*")),
				("All","*")],
			initialdir=self.lastPath,
			initialfile=self.usrbinfile["text"])

		if len(fn) == 0: return False

		self.lastPath = os.path.abspath(os.path.dirname(fn))
		relfn = self.project.relativePath(fn)
		self.usrbinfile["text"] = relfn
		return True

	# ----------------------------------------------------------------------
	def makeUSRBIN(self):
		if self.dicom() is None: return
		if len(self.dicom()._slices) == 0:
			self.flair.notify("No DICOM processed",
				"Please select a DICOM dataset",
				tkFlair.NOTIFY_WARNING)
			return

		if len(self.usrbinfile["text"])==0:
			if not self.saveUsrbinFile():
				return

		# Do the real processing
		t0 = time.time()
		dlg = tkDialogs.ProgressDialog(self, "Creating USRBIN file")

		if not self.dicom()._writeUsrbinInit(self.usrbinfile["text"]):
			self.flair.notify("DICOM to USRBIN",
				"ERROR: Cannot open output Usrbin file",
				tkFlair.NOTIFY_ERROR)
			return

		if self.dicom().slice:
			# Multiple slices
			dlg.setLimits(0, len(self.dicom()._slices))
			for i,(l,fn) in enumerate(self.dicom()._slices):
				# Show progress
				if dlg.show(i, "Processing: "+os.path.basename(fn)): return
				self.dicom()._writeUsrbinSlice(fn,i)
		else:
			# Unique 3D data
			fn = self.dicom()._slices[0][1]
			# Show progress
			dlg.setLimits(0.0, self.dicom().nz)
			dlg.show(0, "Processing: "+os.path.basename(fn))
			def update(s): return dlg.show(s)
			self.dicom()._writeUsrbinSlice(fn, update=update)
		dlg.stop()

		if not self.dicom()._writeUsrbinEnd():
			self.flair.notify("DICOM to USRBIN",
				"Error writing USRBIN file",
				tkFlair.NOTIFY_ERROR)
			return

		self.flair.log("DICOM to Usrbin file conversion finish in %.1f s"%(time.time()-t0))
		self.flair.notify("DICOM to USRBIN",
			"DICOM to Usrbin file conversion finish in %.1f s"%(time.time()-t0))

	# ----------------------------------------------------------------------
	# Create a USRBIN card compatible with the dicom
	# ----------------------------------------------------------------------
	def usrbinCard(self, event=None):
		if self.dicom() is None: return
		self.page.get()
		undoinfo = self.flair.undolistRefresh()

		fn,ext = os.path.splitext(os.path.basename(self.voxelfile["text"]))
		whats = ["dicom",
			10.0,	# Cartesian mesh
			"ENERGY",
			-99.0,	# File
			self.dicom().xmax()/10.0,
			self.dicom().ymax()/10.0,
			self.dicom().zmax()/10.0,

			self.dicom().xmin()/10.0 ,
			self.dicom().ymin()/10.0 ,
			self.dicom().zmin()/10.0 ,
			self.dicom().columns(),
			self.dicom().rows(),
			self.dicom().slices()]

		inp = self.flair.project.input
		if not inp.cardlist: inp.minimumInput()
		pos = inp.bestPosition("USRBIN")

		# add card to input
		card   = Input.Card("USRBIN",whats)
		undoinfo.append(self.flair.addCardUndo(card, pos, False))

		# Finalize
		self.flair.addUndolistRefresh(undoinfo,"Add USRBIN card")
		self.flair.refresh("card")
		self.flair.notify("Added USRBIN card",
			"A USRBIN card corresponding to the dicom file has been added to your input")

	# ----------------------------------------------------------------------
	# Create ROT-DEFI cards for the gantry and patient angles
	# ----------------------------------------------------------------------
	def rotdefiCard(self):
		try:
			dGantry  = float(self.gantry.get())
		except:
			self.flair.notify("Invalid Gantry",
				"Please enter a valid Gantry angle in degrees",
				tkFlair.NOTIFY_WARNING)
			return

		try:
			dPatient = float(self.patient.get())
		except:
			self.flair.notify("Invalid Patient",
				"Please enter a valid Patient angle in degrees",
				tkFlair.NOTIFY_WARNING)
			return

		rGantry  = math.radians(dGantry)
		rPatient = math.radians(dPatient)

		GimbalLock = False

		dTheta = dPhi = dPsi = 0.
		rTheta = rPhi = rPsi = 0.

		eps = 0.0001
		if (abs(abs(dGantry)-90)<eps) and (abs(abs(dPatient)-90)<eps):
			GimbalLock = True
			if dGantry*dPatient > 0: #same sign
				dTheta = 180.
			else:
				dTheta =   0.
			if dGantry > 0:
				dPhi = -90.
			else:
				dPhi =  90.
			dPsi = 0.
		else:
			cosTheta = -math.sin(rGantry)*math.sin(rPatient)
			sinTheta =  math.sqrt(1.0-cosTheta**2)
			rTheta   =  math.acos(cosTheta)

			cosPhi = math.cos(rGantry)/sinTheta
			sinPhi = math.sin(rGantry)*math.cos(rPatient)/sinTheta
			rPhi   = math.acos(cosPhi)
			if sinPhi<0: rPhi = -rPhi

			cosPsi =  math.cos(rPatient)/math.sin(rTheta)
			sinPsi = -math.cos(rGantry) /math.sin(rPatient)/sinTheta
			rPsi   =  math.acos(cosPsi)
			if sinPsi<0: rPsi = -rPsi

			dTheta = math.degrees(rTheta)
			dPhi   = math.degrees(rPhi)
			dPsi   = math.degrees(rPsi)

		self.page.get()
		undoinfo = self.flair.undolistRefresh()
		pos = self._rotdefiCards(1,dPhi,dTheta,dPsi, undoinfo,
					self.project.input.bestPosition("ROT-DEFI"))
		pos = self._rotdefiCards(2,dPhi,dTheta,dPsi, undoinfo, pos)
		self.flair.addUndolistRefresh(undoinfo,"Add ROT-DEFI cards")
		self.flair.refresh("card")
		self.flair.notify("Gantry and Patient transformations",
				"ROT-DEFI cards have been added for the " \
				"specified Gantry and Patient angles")

	# ----------------------------------------------------------------------
	def _rotdefiCards(self, iTrasf, dPhi, dTheta, dPsi, undoinfo, pos):
		transformationNames = ('Beam2CT','CT2Beam','IECgan2CT','CT2IECgan')
		# assuming as convention that even-numbered transforms
		# are inverse of odd-numbered ones
		if iTrasf%2 == 1:
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", [ transformationNames[iTrasf-1], 300, "", dPhi]),
				pos, False))
			pos += 1

			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", [ transformationNames[iTrasf-1], 100, "", dTheta]),
				pos, False))
			pos += 1

			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", [ transformationNames[iTrasf-1], 300, "", dPsi]),
				pos, False))
			pos += 1

		else:
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", [ transformationNames[iTrasf-1], 300, "", -dPsi]),
				pos, False))
			pos += 1

			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", [ transformationNames[iTrasf-1], 100, "", -dTheta]),
				pos, False))
			pos += 1

			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", [ transformationNames[iTrasf-1], 300, "", -dPhi]),
				pos, False))
		return pos

	# ----------------------------------------------------------------------
	def unit2matRelease(self, event):
		act = event.widget.nearest(event.y)
		if act: self.unit2matList.activate(act)
		return self.unit2matEdit(event)

	# ----------------------------------------------------------------------
	def unit2matEdit(self, event):
		if self.dicom() is None: return
		sel = self.unit2matList.curselection()
		if len(sel) != 1: return
		sel = int(sel[0])

		listbox = event.widget
		try: active = listbox.index(ACTIVE)
		except: return
		if active != sel: return

		# Find list box
		lid = self.unit2matList.lists.index(listbox)

		# remember current record
		save = self.unit2matList.get(sel)

		if lid == 0:	# HUmin changed
			edit = tkExtra.InPlaceInteger(listbox)
			# Value cannot be:
			#  1. smaller than the previous minimum
			#  2. bigger than the maximum
			try:
				value = int(edit.value)
				if save[0]==NEW:
					save[0] = save[1] = edit.value
					for i,hu in enumerate(self.unit2matList.get(0,END)):
						try:
							humax = int(hu[1])
						except ValueError:
							self.unit2matList.set(sel,save)
							break
						if value<=humax:
							hu = list(hu)
							hu[1] = str(value-1)
							self.unit2matList.delete(sel)
							self.unit2matList.set(i,hu)
							self.unit2matList.insert(i+1,save)
							self.correctValues()
							active = i+1
							break
				elif value > save[1]:
					self.unit2matList.set(sel, save)	# restore record
				elif sel>0:
					prev = self.unit2matList.get(sel-1)
					if value <= int(prev[0]):
						self.unit2matList.set(sel, save) # restore record
					else:
						# Correct also previous record
						prev[1] = str(value-1)
						self.unit2matList.set(sel-1, prev)
			except (TypeError, ValueError):
				self.unit2matList.set(sel, save)	# restore record

		elif save[0]==NEW:
			return

		elif lid == 1:	# HUmax changed
			edit = tkExtra.InPlaceInteger(listbox)
			# Value cannot be:
			#  1. bigger than the maximum of the next one
			#  2. smaller than the minimum
			try:
				value = int(edit.value)
				if value < save[0]:
					self.unit2matList.set(sel, save)	# restore record
				elif sel<self.unit2matList.size()-1:
					next = self.unit2matList.get(sel+1)
					if value >= int(next[1]):
						self.unit2matList.set(sel, save) # restore record
					else:
						# Correct also previous record
						next[0] = str(value+1)
						self.unit2matList.set(sel+1, next)
			except (TypeError, ValueError):
				self.unit2matList.set(sel, save) # restore record

		elif lid == 2:
			matList = self.dicom().material.materialList(0,True)
			edit = tkExtra.InPlaceList(listbox, values=matList, height=10)
		else:
			edit = tkExtra.InPlaceFloat(listbox)
		if edit.value is None: return

		#if self.unit2matList._sortColumn == lid:
		#	self.unit2matList.sort(lid, self.unit2matList._sortReverse)

		self.unit2matList.activate(active)
		self.unit2matList.selection_set(active)
		self.unit2matModified = True
		return "break"

	# ----------------------------------------------------------------------
	# Scan the list to correct the hounsfield unit values if needed
	# ----------------------------------------------------------------------
	def correctValues(self):
		active = self.unit2matList.index(ACTIVE)
		sel = self.unit2matList.curselection()
		for i,hu in enumerate(self.unit2matList.get(0,END)):
			hu = list(hu)
			try:    humin = int(hu[0])
			except: pass
			try:    humax = int(hu[1])
			except: humax = humin

			if i==0:
				huprev = humax
				continue

			changed = False
			if humin != huprev+1:
				humin = huprev+1
				hu[0] = str(humin)
				changed = True
			if humax < humin:
				humax = humin
				hu[1] = str(humin)
				changed = True
			if changed: self.unit2matList.set(i, hu)
			huprev = humax
		self.unit2matList.activate(active)
		for s in sel: self.unit2matList.selection_set(s)

	# ----------------------------------------------------------------------
	def getHounsfieldFromList(self):
		if self.dicom() is None: return
		Dicom.Hounsfield.reset()
		del self.dicom().hounsfield[:]

		self.unit2matList.sort(0,False)

		prev = -9999
		dups = []

		for elem in self.unit2matList.get(0,END):
			try: humin = int(elem[0])
			except: continue
			try: humax = int(elem[1])+1
			except: continue
			if humax == prev: dups.append(humax)
			prev = humax
			material = elem[2]
			try:	crho_min = float(elem[3])
			except: crho_min = 1.0
			try:	crho_max = float(elem[4])
			except: crho_max = 1.0
			try:	cdedx_min = float(elem[5])
			except: cdedx_min = 1.0
			try:	cdedx_max = float(elem[6])
			except: cdedx_max = 1.0
			hu = Dicom.Hounsfield(humax, material, crho_min, crho_max, cdedx_min, cdedx_max)
			self.dicom().hounsfield.append(hu)

		if dups:
			self.flair.notify("Duplicate Hounsfield units",
				"Duplicate units found "+str(dups),
				tkFlair.NOTIFY_WARNING)
			self.unit2matList.selection_clear(0,END)
			for i,hu in enumerate(self.unit2matList.lists[0].get(0,END)):
				if hu in dups:
					self.unit2matList.selection_set(i)

	# ----------------------------------------------------------------------
	def add(self):
		# add an empty record (only once)
		for hu in self.unit2matList.lists[0].get(0,END):
			if hu == NEW: return
		self.unit2matList.insert(END,(NEW,"","","","","",""))
		self.unit2matList.selection_clear(0,END)
		self.unit2matList.activate(END)
		self.unit2matList.selection_set(END)
		self.unit2matModified = True
		self.setModified(True)

	# ----------------------------------------------------------------------
	def delete(self, event=None):
		sel = list(map(int,self.unit2matList.curselection()))
		if len(sel) == 0: return
		sel.reverse()
		for item in sel:
			self.unit2matList.delete(item)
		self.unit2matModified = True
		self.setModified(True)
		self.correctValues()

	# ----------------------------------------------------------------------
	def clone(self, event=None):
		sel = list(map(int,self.unit2matList.curselection()))
		if len(sel) == 0: return
		for item in sel:
			self.unit2matList.insert(END, self.unit2matList.get(item))
		self.unit2matModified = True
		self.setModified(True)
		self.correctValues()
