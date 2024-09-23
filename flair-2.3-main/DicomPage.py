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
import time
try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkFlair
import tkExtra
import tkDialogs

import Dicom
import Ribbon
import Project
import FlairRibbon

try:
	import pydicom as dicom
except ImportError:
	try:
		import dicom
	except ImportError:
		dicom = None
try:
	import DicomTab
	import InfoTab
	import SliceTab
	import VoxelTab
	import RTPlanTab
	import RTViewerTab
except ImportError:
	dicom = None

#===============================================================================
# Dicom Frame class
#===============================================================================
class DicomPage(FlairRibbon.FlairPage):
	"""Generate Voxel and USRBINs from DICOM files CT, PET, RT..."""

	_name_ = "Dicom"
	_icon_ = "dicom"

	#-----------------------------------------------------------------------
	def init(self):
		# Variables
		self.dicom     = Dicom.Dicom()
		self.dicominfo = None

	#-----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		if dicom is None: return

		# ========== Tools ===========
		group = Ribbon.LabelGroup(self.ribbon, "Process")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		kw = {
			"compound" : TOP,
			"variable" : self._tabVar,
			"command"  : self.changeTab,
			"padx"	   : 6
		}

		# ---
		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["dicom32"],
				text="DICOM",
				value=0,
				**kw)
		b.pack(side=LEFT, fill=Y, expand=YES)
		tkExtra.Balloon.set(b, "Select DICOM set")

		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["slice32"],
				text="Slice",
				value=1,
				**kw)
		tkExtra.Balloon.set(b, "Select/View slices")
		b.pack(side=LEFT, fill=Y, expand=YES)

		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["total32"],
				text="Information",
				value=2,
				**kw)
		tkExtra.Balloon.set(b, "View information on selected slices")
		b.pack(side=LEFT, fill=Y, expand=YES)

		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["voxel32"],
				text="Voxel",
				value=3,
				**kw)
		tkExtra.Balloon.set(b, "Create VOXEL or USRBIN from DICOM")
		b.pack(side=LEFT, fill=Y, expand=YES)

		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["plan32"],
				text="RTPlan",
				value=4,
				**kw)
		tkExtra.Balloon.set(b, "Process RTPlan Data")
		b.pack(side=LEFT, fill=Y, expand=YES)


		b = Ribbon.LabelRadiobutton(group.frame,
				image=tkFlair.icons["dicom32"],
				text="RTViewer",
				value=5,
				**kw)
		tkExtra.Balloon.set(b, "View RTPlan Dose")
		b.pack(side=LEFT, fill=Y, expand=YES)

	#-----------------------------------------------------------------------
	# Create Project page
	#-----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self)

		if dicom is None:
			t = Message(self.frame,
				text =	"pydicom library is not installed\n\n" \
					"Website: http://code.google.com/p/pydicom/\n\n" \
					"Install:\n" \
					"\tUbuntu\t# apt-get install python-dicom\n" \
					"\tall\t# easy_install pydicom" \
					"\t\t(Requires the python-setuptools)",
				aspect = 600,
				relief=RAISED)
			t.pack(expand=YES, fill=BOTH)
			b = Button(self.frame,
				text = "Help",
				command=self.help)
			b.pack()
			return

		# --------- Dicom tab ---------
		self.dicomTab = DicomTab.DicomTab(self.frame, self)
		self.addTab(self.dicomTab)

		# --------- Slice tab ---------
		self.sliceTab = SliceTab.SliceTab(self.frame, self)
		self.addTab(self.sliceTab)

		# --------- Info tab ---------
		self.infoTab = InfoTab.InfoTab(self.frame, self)
		self.addTab(self.infoTab)

		# --------- Voxel tab ---------
		self.voxelTab = VoxelTab.VoxelTab(self.frame, self)
		self.addTab(self.voxelTab)

		# -------- RTPlan tab ---------
		self.RTPlanTab = RTPlanTab.RTPlanTab(self.frame, self)
		self.addTab(self.RTPlanTab)

		# -------- RTViewer tab ---------
		self.RTViewerTab = RTViewerTab.RTViewerTab(self.frame, self)
		self.addTab(self.RTViewerTab)

	# ----------------------------------------------------------------------
	# Update fields
	# ----------------------------------------------------------------------
	def get(self, event=None):
		if self.page is None or dicom is None: return
		if self.project is None: return
		if self.dicominfo is None: return

		oldhash = self.dicominfo.hash()
		self.dicomTab.get()
		self.sliceTab.get()
#		self.infoTab.get()
		self.voxelTab.get()
#		self.RTPlanTab.get()
#		self.RTViewerTab.get()
		if self.dicominfo.hash() != oldhash:
			self.setModified(True)

	# ----------------------------------------------------------------------
	# Refresh directory listing
	# ----------------------------------------------------------------------
	def refresh(self):
		FlairRibbon.FlairPage.refresh(self)
		if self.page is None or dicom is None: return

		self.dicomTab.refresh()
		self.sliceTab.refresh()
#		self.infoTab.refresh()
		self.voxelTab.refresh()
		self.RTPlanTab.refresh()
		self.RTViewerTab.refresh()
#		self.RTOptTab.refresh()

	# ----------------------------------------------------------------------
	# Change selection
	# ----------------------------------------------------------------------
	def datasetChange(self):
		di = self.dicominfo
		self.get()
		if di == self.selectedDicom(): return

		# Refresh only the rest of tabs
		self.sliceTab.refresh()
#		self.infoTab.refresh()
		self.voxelTab.refresh()
		self.RTPlanTab.refresh()
		self.RTViewerTab.refresh()
#		self.RTOptTab.refresh()

	# ----------------------------------------------------------------------
	def selectedDicom(self):
		sel = list(map(int,self.dicomTab.dicomList.curselection()))
		if len(sel) != 1: return None
		self.dicominfo = self.project.dicomInfo(self.dicomTab.dicomList.get(sel[0])[0])
		return self.dicominfo

	# ----------------------------------------------------------------------
	# Load dicom structure
	# @return True on success (or already existing)
	#	  False otherwise
	# ----------------------------------------------------------------------
	def loadDicom(self):
		self.selectedDicom()

		if self.dicominfo is None:
			self.flair.notify("DICOM select",
				"Please select only ONE DICOM dataset")
			return True

		if self.dicom.uid == self.dicominfo["uid"] and self.dicom.scanned():
			# Check last slices loaded
			if not self.dicom.setLimits(*self.sliceTab.limits()):
				return None
		else:
			self.dicom.setLimits(*self.sliceTab.limits())

		self.log = self.flair.newLog("Dicom",self.dicominfo["uid"])
		self.dicom.setLog(self.log)

		oldSize = self.log.size()	# remember size before loading
		err = self._loadDicom()
		if isinstance(err,str) or self.log.size() != oldSize:
			if isinstance(err,str): self.log(err)
			self.flair.notify("Errors during loading",
				"Error or warning messages appeared while loading dicom slices. " \
				"Open the Output window for more information.", tkFlair.NOTIFY_ERROR,
				"Output", self.flair.showOutput)
		return err

	# ---------------------------------------------------------------------
	def _loadDicom(self):
		files = self.dicominfo.files
		path  = self.dicominfo["directory"]

		# Find slices and zmin
		self.dicom.clear()

		t = time.time()

		dlg = tkDialogs.ProgressDialog(self.page, "Analyzing slices")

		if self.dicominfo["frames"] == 0:
			# dicom with multiple slices
			z = min([x[1] for x in files])
			fn = os.path.join(path,files[0][0])
			if not self.dicom.init(fn, len(files), z):
				#self.dicom.clear()
				dlg.stop()
				return True #"Error loading slice %s"%(files[0][0])
			self.dicom.calcLimits()

			dlg.setLimits(0.0, len(files))
			for i,(fn,z,sop) in enumerate(files):
				if not self.dicom._scanSlice(os.path.join(path,fn)):
					# Error
					self.dicom.clear()
					dlg.stop()
					return "Error loading slice %s"%(fn)

				# Show progress
				if dlg.show(i,"Slice: "+fn):
					self.dicom.clear()
					return "Progress bar internal error"
		else:
			# Unique slice dicom
			# get selected dicom file from slice list
			sel = list(map(int,self.sliceTab.sliceList.curselection()))
			if len(sel) == 1:
				fn = self.sliceTab.sliceList.get(sel[0])[0]
			else:
				if self.sliceTab.sliceList.size() != 1:
					self.dicom.clear()
					dlg.stop()
					return "No slice selected"
				fn = self.sliceTab.sliceList.get(0)[0]

			if not self.dicom.init(os.path.join(path,files[sel[0]][0])):
				#self.dicom.clear()
				dlg.stop()
				return "Error loading slice %s"%(files[0][0])
			self.dicom.calcLimits()	# calculate correctly the limits now we know the dz

			dlg.setLimits(0.0, self.dicom.nz)
			dlg.show(0,"Slice: "+fn)

			def update(n): return dlg.show(n)
			if not self.dicom._scanSlice(os.path.join(path,fn),update=update):
				self.dicom.clear()
				dlg.stop()
				return "Error scanning slice %s"%(fn)

		# End scanning
		err = self.dicom._scanSliceEnd()

		# Close dialog
		dlg.stop()

		t = time.time() - t
		self.flair.log("Analyzing DICOM slices finished in %.1fs"%(t))
		self.flair.notify("Analyzing DICOM",
			"Analyzing DICOM slices finished in %.1fs"%(t))

		return None

	# ----------------------------------------------------------------------
	# Add a new dicom dataset
	# ----------------------------------------------------------------------
	def addDicomInfo(self, path, uid):
		di = Project.DicomInfo(uid)	# should change the name later
		self.project.dicoms.append(di)
		di["uid"]	= uid
		di["directory"] = path
		di["frames"]	= Dicom.Scan.frames(uid)
		di.files	= Dicom.Scan.files(uid)
		di.modality	= Dicom.Scan.modality(uid)
		return di

	# ----------------------------------------------------------------------
	# Find specific slice from SOPInstanceUID
	# ----------------------------------------------------------------------
	def findDicomFromSOPInstanceUID(self, sop):
		# First time scan all dicom info sets
		for di in self.project.dicoms:
			for fn,z,s in di.files:
				if sop == s:
					return di
		return None
