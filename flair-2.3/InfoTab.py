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

import bmath
import tkFlair
import tkExtra
import FlairRibbon

try:
	import pydicom as dicom
except ImportError:
	import dicom
import Dicom

#===============================================================================
# Info Tab
#===============================================================================
class InfoTab(FlairRibbon.FlairTab):
	#----------------------------------------------------------------------
	def __init__(self, master, page, **kw):
		FlairRibbon.FlairTab.__init__(self, master, page, **kw)

		topf = Frame(self)
		topf.pack(side=TOP, fill=X)

		# ---
		row = 0
		Label(topf, text="Global:", foreground="DarkRed").grid(row=row, column=0, sticky=W)

		# ---
		row += 1
		Label(topf, text="UID:").grid(row=row, column=0, sticky=E)
		self.uid = Label(topf, justify=LEFT, anchor=W,
				width=20,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.uid.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(topf, text="Patients Name:").grid(row=row, column=0, sticky=E)
		self.patientname = Label(topf, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.patientname.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(topf, text="Modality:").grid(row=row, column=0, sticky=E)
		self.modality = Label(topf, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.modality.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(topf, text="Date/Time:").grid(row=row, column=0, sticky=E)
		self.datetime = Label(topf, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.datetime.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(topf, text="Patients Age:").grid(row=row, column=0, sticky=E)
		self.age = Label(topf, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.age.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(topf, text="Patients Sex:").grid(row=row, column=0, sticky=E)
		self.sex = Label(topf, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.sex.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(topf, text="Institution:").grid(row=row, column=0, sticky=E)
		self.institution = Label(topf, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.institution.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(topf, text="Study:").grid(row=row, column=0, sticky=E)
		self.study = Label(topf, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.study.grid(row=row, column=1, sticky=EW)

		topf.grid_columnconfigure(1, weight=1)

		# ====== Create frames for modality based information ====
		self.frames = {}

		f = Frame(self)
		self.frames["CT"] = f
		self.frames["MR"] = f
		self.frames["PT"] = f
		self.frames["RTDOSE"] = f

		row = 0
		Label(f, text="Dataset:", foreground="DarkRed").grid(row=row, column=0, sticky=W)

		# ---
		row += 1
		Label(f, text="# slices:").grid(row=row, column=0, sticky=E)
		self.nslices = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.nslices.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(f, text="# voxels [Nx x Ny x Nz]:").grid(row=row, column=0, sticky=E)
		self.nxnynz = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.nxnynz.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(f, text="Voxel size [dx x dy x dz]:").grid(row=row, column=0, sticky=E)
		self.dxdydz = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.dxdydz.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(f, text="Dimensions [Dx x Dy x Dz]:").grid(row=row, column=0, sticky=E)
		self.DxDyDz = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.DxDyDz.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(f, text="Position:").grid(row=row, column=0, sticky=E)
		self.position = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.position.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(f, text="Rescale:").grid(row=row, column=0, sticky=E)
		self.rescale = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.rescale.grid(row=row, column=1, sticky=EW)

		# ===
		row += 1
		Label(f, text="Selected Volume:", foreground="DarkRed").grid(row=row, column=0, sticky=W)

		# ---
		row += 1
		Label(f, text="# voxels [Nx x Ny x Nz]:").grid(row=row, column=0, sticky=E)
		self.crop_nxnynz = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.crop_nxnynz.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(f, text="Dimensions [Dx x Dy x Dz]:").grid(row=row, column=0, sticky=E)
		self.crop_DxDyDz = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.crop_DxDyDz.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(f, text="Position:").grid(row=row, column=0, sticky=E)
		self.crop_position = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.crop_position.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(f, text="Unit range [Min..Max]:").grid(row=row, column=0, sticky=E)
		self.unitrange = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.unitrange.grid(row=row, column=1, sticky=EW)

		# ---
		row += 1
		Label(f, text="Unique Units [Regions]:").grid(row=row, column=0, sticky=E)
		self.uniqueunits = Label(f, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.uniqueunits.grid(row=row, column=1, sticky=EW)

		f.grid_columnconfigure(1, weight=1)

		#---------------------------------------------------------------
		self.frames["RTSTRUCT"] = f = Frame(self)
		self.rtlist = tkExtra.MultiListbox(f,
					(("Id",     3, None),
					 ("Name",   16, None),
					 ("Dicom Volume (cm3)",   8, None),
					 ("Voxel Volume (cm3)",   8, None),
					 ("Polygon Volume (cm3)", 8, None),
					 ("Dicom/Polygon (+%)",   4, None),
					 ("Voxel/Polygon (+%)",   4, None)))
		self.rtlist.pack(fill=BOTH, expand=YES)

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairTab.createRibbon(self)

	# ----------------------------------------------------------------------
	def dicom(self):	return self.page.dicom

	# ----------------------------------------------------------------------
	def showModality(self, modality):
		frame = None
		for k,f in list(self.frames.items()):
			if k==modality:
				frame = f
			else:
				f.pack_forget()

		# Show frame in the end, since multiple modalities share the
		# same frame, eg CT,MR,RTDOSE and one will cancel the other
		if frame:
			frame.pack(side=BOTTOM, fill=BOTH, expand=YES)

	# ----------------------------------------------------------------------
	def activate(self):
		err = self.page.loadDicom()
		self.showModality(self.dicom().modality)
		self.clearVoxelInformation()
		di = self.page.selectedDicom()
		if di is None:
			self.flair.notify("No DICOM selected",
				"Please select a DICOM set first")
			return
		self.headerVoxelInformation()

		if err is None:
			self.updateVoxelInformation()

		elif self.dicom().modality == "RTSTRUCT":
			fn = os.path.join(di["directory"], di.files[0][0])
			self.updateRtStruct(fn)

	# --------------------------------------------------------------------
	# Clear voxel information
	# --------------------------------------------------------------------
	def clearVoxelInformation(self):
		self.uid["text"]         = ""
		self.datetime["text"]    = ""
		self.patientname["text"] = ""
		self.institution["text"] = ""
		self.study["text"]       = ""
		self.nslices["text"]     = ""
		self.nxnynz["text"]      = ""
		self.dxdydz["text"]      = ""
		self.DxDyDz["text"]      = ""
		self.position["text"]    = ""
		self.rescale["text"]     = ""
		self.unitrange["text"]   = ""
		self.uniqueunits["text"] = ""

	# --------------------------------------------------------------------
	def headerVoxelInformation(self):
		uid = self.page.dicominfo["uid"]
		fn  = os.path.join(self.page.dicominfo["directory"],
				self.page.dicominfo.files[0][0])
		header = Dicom.dicomHeader(self.page.dicominfo["uid"], fn)

		self.uid["text"]         = header[0]
		self.modality["text"]    = "%s : %s" % \
				(header[1], Dicom.MODALITY.get(header[1],"?"))
		self.datetime["text"]    = header[2]
		self.patientname["text"] = header[3]
		self.age["text"]         = header[4]
		self.sex["text"]         = header[5]
		self.institution["text"] = header[6]
		self.study["text"]       = header[7]
		self.nslices["text"]     = "%d"%(len(self.page.dicominfo.files))

	# --------------------------------------------------------------------
	# Update RT structure information
	# --------------------------------------------------------------------
	def updateRtStruct(self, filename):
		self.rtlist.delete(0,END)

		try:
			dataset = dicom.read_file(filename, force=True)
		except:
			return
		ctinfo = None
		for cnt in dataset.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].ContourImageSequence:
			ctinfo = self.page.findDicomFromSOPInstanceUID(cnt.ReferencedSOPInstanceUID)
			if ctinfo is not None:
				break

		# Load first file
		if ctinfo:
			ctfilename = os.path.join(ctinfo["directory"], ctinfo.files[0][0])
		else:
			ctfilename = None

		rt = Dicom.RTStruct(dataset, ctfilename) #, x0, y0, dx, dy, dz)
		for idx in sorted(list(rt.idx.values())):
			name  = rt.names[idx]
			d,v,p = rt.volumes[idx]
			v /= 1000.0
			p /= 1000.0
			if d>0.0:
				dp = 100.0*(d/p-1.0)
			else:
				dp = 0.0
			try:
				vp = 100.0*(v/p-1.0)
			except ZeroDivisionError:
				vp = 0.0
			self.rtlist.insert(END, (idx, name, round(d,4), round(v,4), round(p,4), round(dp,2), round(vp,2)))

	# --------------------------------------------------------------------
	# Scan slices to update voxel information
	# --------------------------------------------------------------------
	def updateVoxelInformation(self, event=None):
		self.nxnynz["text"] = "%d x %d x %d" % \
				(self.dicom().nx, self.dicom().ny, self.dicom().nz)

		self.dxdydz["text"] = "%g x %g x %g mm^3" % \
				(self.dicom().dx,
				 self.dicom().dy,
				 self.dicom().dz)

		self.DxDyDz["text"] = "%g x %g x %g mm^3" % \
				(self.dicom().nx*self.dicom().dx,
				 self.dicom().ny*self.dicom().dy,
				 self.dicom().nz*self.dicom().dz)

		self.position["text"] = "[%g, %g, %g] mm" % \
			(self.dicom().x, self.dicom().y, self.dicom().z)

		if self.dicom().b > 0:
			self.rescale["text"] = "%g * val + %g"%(self.dicom().m, self.dicom().b)
		elif self.dicom().b < 0:
			self.rescale["text"] = "%g * val - %g"%(self.dicom().m, -self.dicom().b)
		elif self.dicom().b == 0:
			self.rescale["text"] = "%g * val"%(self.dicom().m)

		# Cropped area
		nx = self.dicom().columns()
		ny = self.dicom().rows()
		nz = self.dicom().slices()

		self.crop_nxnynz["text"] = "%d x %d x %d"%(nx, ny, nz)
		self.crop_DxDyDz["text"] = "%g x %g x %g mm^3" % \
				(nx*self.dicom().dx,
				 ny*self.dicom().dy,
				 nz*self.dicom().dz)
		self.crop_position["text"]    = "[%g, %g, %g] mm" % \
			(self.dicom().xmin(), self.dicom().ymin(), self.dicom().zmin())

		self.unitrange["text"]   = "%d .. %d" % \
			(self.dicom().minimum, self.dicom().maximum)
		self.uniqueunits["text"] = "%d"%(len(self.dicom().reg2hu))
