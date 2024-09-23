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
# Author: Wioletta.Kozlowska@cern.ch
# Date:	16-Aug-2015

from tkinter import *

import os
import sys
import math
try:
	import pydicom as dicom
except ImportError:
	import dicom

import Data
import Ribbon
import tkFlair
import tkExtra
import Palette
import FlairRibbon
import DicomBrowser
import bFileDialog

try:
	import ZoomImage
	from PIL import Image
	from PIL import ImageDraw
	from PIL import ImageOps
	have_PIL = True
except:
	have_PIL = False

try:
	import numpy
	have_numpy = True
except ImportError:
	have_numpy = False

from Unicode import RIGHTWARDS_ARROW

_PALETTE_RES   = 1
_PALETTE_MAX   = 10
_PALETTE_SIZE  = 30

#===============================================================================
# RTViewer Tab
#===============================================================================
class RTViewerTab(FlairRibbon.FlairTab):
	#----------------------------------------------------------------------
	def __init__(self, master, page, **kw):
		FlairRibbon.FlairTab.__init__(self, master, page, **kw)

		# -------------- Slices ---------------
		frame = Frame(self)
		frame.pack(fill=BOTH, expand=YES)

#		 self.roiList = tkExtra.MultiListbox(frame,
#						    (("Name", 10, None),("z",10,None)))
#		self.roiList.grid(column=4, sticky=NS)
#		self.roiList.bind("<Double-1>",        self.viewer)
#		self.roiList.bind("<<ListboxSelect>>", self.roiSelect)
#		self.roiList.bind("<Control-Key-f>", lambda e,s=self : s.findString.focus_set())
#		self.roiList.sortAssist = None

		# --------------- Canvas -----------
		self.canvas = []
		for i in ("RTDOSE", "FLUKA", "Difference"):
			canvas = ZoomImage.ZoomImage(frame)#, takefocus=True)
			canvas.pack(side=LEFT, fill=BOTH, expand=YES)
			canvas.bindMotion(wheel=False)
			#canvas.bindArrows()
			#canvas.bind('<<Zoom>>',  self.imageZoomChanged)
			#canvas.bind('<<View>>',  self.viewChanged)
			canvas.bind('<Button-4>',	self.imageZoomIn);
			canvas.bind('<Button-5>',	self.imageZoomOut);
			canvas.bind('<Key-1>',		self.imageZoomOne)
			canvas.bind('<B2-Motion>',	self.imageScanDrag)
			canvas.bind('<ButtonRelease-2>', self.imageScanRelease)
			#canvas.create_text((5,5), anchor=NW, fill="white", text=i)
			self.canvas.append(canvas)

			l = Label(canvas, text=i, anchor=N, justify=CENTER)
			l.place(relx=0.5, y=5, anchor=N)

		frame = Frame(self)
		frame.pack(fill=BOTH)

		row,col = 0,0
		Label(frame, text="Min:").grid(row=row, column=col, sticky=E)
		col += 1
		self.palmin = Scale(frame, from_=0, to_=_PALETTE_MAX,
#					showvalue=0,
					resolution=_PALETTE_RES,
					orient=HORIZONTAL,
					command=self.paletteLevel)
		self.palmin.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.palmin, "Minimum palette values to display")

		row,col = 0,col+1
		Label(frame, text="Max:").grid(row=row, column=col, sticky=E)
		col += 1
		self.palmax = Scale(frame, from_=_PALETTE_RES, to_=_PALETTE_MAX,
#					showvalue=0,
					resolution=_PALETTE_RES,
					orient=HORIZONTAL,
					command=self.paletteLevel)
		self.palmax.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.palmax, "Maximum palette values to display")

		col += 1
		Label(frame, text=" x :").grid(row=row, column=col, sticky=E)
		col += 1
		self.palmult = Scale(frame, from_=1, to_=20,
#					showvalue=0,
					orient=HORIZONTAL,
					command=self.paletteLevel)
		self.palmult.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.palmult, "Multiply palette values for difference")

		row,col = row+1,0
		Label(frame, text="Window Center:").grid(row=row, column=col, sticky=S+E)

		col += 1
		self.center = Scale(frame, from_=0, to_=5000,
					showvalue=0,
					orient=HORIZONTAL,
					command=self.levelChanged)
		self.center.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.center, "Center or Level of values to display")

		# ---
		col += 1
		Label(frame, text="Width:").grid(row=row, column=col, sticky=S+E)

		col += 1
		self.window = Scale(frame, from_=0, to_=5000,
					showvalue=0,
					orient=HORIZONTAL,
					command=self.levelChanged)
		self.window.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.window, "Window around level/center values to display")

		# -------------
		row,col = row+1,0
		Label(frame, text="Slice:").grid(row=row, column=col, sticky=S+E)

		col += 1
		self.sliceScale = Scale(frame, from_=0, to_=100,
					orient=HORIZONTAL,
					command=self.sliceChange)
		self.sliceScale.grid(row=row, column=col, columnspan=5, sticky=EW)
		tkExtra.Balloon.set(self.sliceScale, "Slice/Frame to display")

		frame.grid_columnconfigure(1, weight=1)
		frame.grid_columnconfigure(3, weight=1)
		frame.grid_columnconfigure(5, weight=1)

		self.view	= "Axial"
		self.lastPath	= ""
		self.dataset	= None	# current loaded dataset
		self.frame	= 0	# Number of frames for scale
		self.PlanFile	= False	 # Flag for existing Dose Files
		self.CalcFile	= False
		self.CTInfo	= None	# CT project DicomInfo structure related to this dataset
		self.arrayDicom = None

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairTab.createRibbon(self)

		#==
		group = Ribbon.LabelGroup(self.ribbon, "View")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		b = Ribbon.LabelRadiobutton(group.frame,
					image=tkFlair.icons["axial32"],
					text="Axial",
					compound=TOP,
					command=lambda:self.viewChange("Axial"),
					value = 0,
					padx=5,
					pady=10)
		b.pack(side=LEFT, fill=Y, expand=YES)
		tkExtra.Balloon.set(b, "Axial view")
		b.select()

		b = Ribbon.LabelRadiobutton(group.frame,
					image=tkFlair.icons["coronal32"],
					text="Coronal",
					compound=TOP,
					command=lambda:self.viewChange("Coronal"),
					value = 1,
					padx=5)
		b.pack(side=LEFT, fill=Y, expand=YES)
		tkExtra.Balloon.set(b, "Coronal View")

		b = Ribbon.LabelRadiobutton(group.frame,
					image=tkFlair.icons["sagittal32"],
					text="Sagittal",
					compound=TOP,
					command=lambda:self.viewChange("Sagittal"),
					value = 2,
					padx = 5)
		b.pack(side=LEFT, fill=Y, expand=YES)
		tkExtra.Balloon.set(b, "Sagittal View")

		# ========== View ===========
		group = Ribbon.LabelGroup(self.ribbon, "Dose")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		row,col = 0,0
		Label(group.frame, text="Planned Dose:", background=Ribbon._BACKGROUND).grid(
				row=row, column=col, pady=1, sticky=E)
		col += 1
		self.PDose = Label(group.frame,
				foreground=tkFlair._ELABEL_FOREGROUND_COLOR,
				background=tkFlair._ELABEL_BACKGROUND_COLOR,
				width=30,
				relief=SUNKEN, anchor=W)
		self.PDose.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.PDose, "Planned Dose from RTPLAN")

		col += 1
		b = Button(group.frame, image=tkFlair.icons["load"], pady=0,
				command=self.loadPDose)
		b.grid(row=row, column=col)
		tkExtra.Balloon.set(b, "Load USRBIN file")

		row,col = 1,0
		Label(group.frame, text="Calculated Dose:", background=Ribbon._BACKGROUND).grid(
				row=row, column=col, pady=1, sticky=E)
		col += 1
		self.CDose = Label(group.frame,
				foreground=tkFlair._ELABEL_FOREGROUND_COLOR,
				background=tkFlair._ELABEL_BACKGROUND_COLOR,
				width=30,
				relief=SUNKEN, anchor=W)
		self.CDose.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.CDose, "Calculated Dose from RTPLAN")

		col += 1
		b = Button(group.frame, image=tkFlair.icons["load"], pady=0,
				command=self.loadCDose)
		b.grid(row=row, column=col)
		tkExtra.Balloon.set(b, "Load USRBIN file")

		# ========== RTStruct ===========
		group = Ribbon.LabelGroup(self.ribbon, "ROI")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		row,col = 0,0
		Label(group.frame, text="RTSTRUCT:", background=Ribbon._BACKGROUND).grid(row=row, column=col, sticky=E)
		col += 1
		self.rtstruct = tkExtra.Combobox(group.frame, width=30, foreground="DarkBlue")
		self.rtstruct.grid(row=row, column=col, columnspan=4, sticky=EW)

		col = 0
		row += 1
		Label(group.frame, text="ROI:", background=Ribbon._BACKGROUND).grid(row=row, column=col, sticky=E)
		col += 1
		self.roilist = tkExtra.Combobox(group.frame, width=30, foreground="DarkBlue")
		self.roilist.grid(row=row, column=col, columnspan=4, sticky=EW)

	# ----------------------------------------------------------------------
	def dicom(self):	return self.page.dicom

	# ----------------------------------------------------------------------
	def refresh(self):
		"""Refresh the view in the beginning and set the variables"""

		di = self.page.dicominfo
		if di is None:
			self._photo = None
			self.view = "Axial"
			self.lastPath = ""
			self.frame   = 0	# Number of frames for scale
			return

		self.dataset = None

		try:
			if di.modality != "CT":# and  di.modality != "RTSTRUCT":
				self.flair.notify("Wrong DICOM selected",
				"Please select a CT or RTSTRUCT set first")
				return
		except:
			self.flair.notify("No DICOM selected",
				"Please select a CT or RTSTRUCT set first")
			return

#		for struct in self.project.dicoms:
#			if struct.modality=="RTSTRUCT":
#				self.rtstruct.insert(END,struct.name)
#				rtuid = self.rtstruct.get()
#
#		for struct in self.project.dicoms:
#			if struct["uid"] == rtuid:
#				try:
#					self.rtstruct = dicom.read_file(
#							os.path.join(struct["directory"], struct.files[0][0]),
#							force=True)
#					for roi in self.rtstruct.StructureSetROISequence:
#						self.roilist.insert(END,(str(roi.ROIName)))
#				except:
#					return False
#				break

		# TODO: Correct listing of roi
		# Changing the center/window will update the image/crop
		self.window.set(float(di.get("window",0)))
		self.center.set(float(di.get("center",0)))
		self.sliceScale.set(float(di.get("sliceScale",0)))
		self.palmin.set(float(di.get("palmin",0)))
		self.palmax.set(float(di.get("palmax",_PALETTE_MAX)))
		self.palmult.set(float(di.get("palmult",_PALETTE_MAX)))

	# ----------------------------------------------------------------------
	def get(self):
		"""Get the variables from the sliders"""

		if self.dicom() is None: return
		di = self.page.dicominfo
		if di is None: return

		di["window"] = self.window.get()
		di["center"] = self.center.get()
		di["palmin"] = self.palmin.get()
		de["palmax"] = self.palmax.get()
		de["palmult"]= self.palmult.get()
		di["rtstruct"]	= self.rtstruct.get()

	# ----------------------------------------------------------------------
	# Fill the dicom array
	# ----------------------------------------------------------------------
	def activate(self, event=None):
		"""Fill the dicom array"""

		di = self.page.selectedDicom()

		if di is None:
			self.flair.notify("No DICOM selected",
				"Please select a DICOM set first")
			return

		if (di.modality != "CT" and di.modality != "RTSTRUCT"):
			self.flair.notify("Wrong DICOM selected",
				"Please select a CT or RTSTRUCT file first")
			return

		try:
			self.dataset = dicom.read_file(
						os.path.join(self.page.dicominfo["directory"], di.files[0][0]),
						defer_size=256,
						force=True)

		except dicom.filereader.InvalidDicomError:
			self.dataset = None
			return

		self.CTInfo = di

		# Init dicom with the first slice
		self.CT_low_pos  = [0,0,0]
		self.CT_high_pos = [1,1,1]
		self.CT_size	 = [1,1,1]
		self.CT_dim	 = [1,1,1]
		self.CT_bin	 = [1,1,1]

		try:
			PixelDims    = [int(self.dataset.Rows), int(self.dataset.Columns), len(di.files[:])]
			PixelSpacing = [float(self.dataset.PixelSpacing[0]),
					float(self.dataset.PixelSpacing[1]),
					float(self.dataset.SliceThickness)]
			Position     = [float(self.dataset.ImagePositionPatient[0]),
					float(self.dataset.ImagePositionPatient[1]),
					float(self.dataset.ImagePositionPatient[2])]

		except AttributeError:
			PixelDims    = [1,1,1]
			PixelSpacing = [1,1,1]
			Position     = [0,0,0]
			self.flair.notify("CT DICOM problem",
				"Please select correct CT file")
			return

		#Calculating max values for axes, and resize factor
		self.CT_dim	 = PixelDims
		self.CT_bin	 = PixelSpacing
		for i in range(len(PixelSpacing)):
			self.CT_size[i]     = int((PixelDims[i])*PixelSpacing[i])
			self.CT_low_pos[i]  = Position[i]-PixelSpacing[i]*0.5
			self.CT_high_pos[i] = self.CT_low_pos[i]+((PixelDims[i])*PixelSpacing[i])

		self.arrayDicom = numpy.zeros(self.CT_dim, dtype=self.dataset.pixel_array.dtype)

		try:
			i=0
			for filenameDCM in sorted(di.files, key=lambda z:di.files[:][1])[:]:
				# read the file
				ds = dicom.read_file(os.path.join(self.CTInfo["directory"],filenameDCM[0]),force=True)
				self.arrayDicom[:, :,i] = ds.pixel_array
				i+=1
		except AttributeError:
			self.flair.notify("CT DICOM problem",
				"Please select correct CT file")
			return

	# ----------------------------------------------------------------------
	# For RTSTRUCT return a dictionary with sop and CT filename
	# ----------------------------------------------------------------------
	def scan4ReferencedCT(self):
		"""For RTSTRUCT return a dictionary with sop and CT filename"""
		self.CTInfo = None
		for cnt in self.dataset.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].ContourImageSequence:
			self.CTInfo = self.page.findDicomFromSOPInstanceUID(cnt.ReferencedSOPInstanceUID)
			if self.CTInfo is not None:
				return

	# ----------------------------------------------------------------------
	# Zoom all canvases
	# ----------------------------------------------------------------------
#	def imageZoomChanged(self, event):
#		z = event.widget.zoom
#		for canvas in self.canvas:
#			if canvas is event.widget: continue
#			canvas.zoomImage(z)

	# ----------------------------------------------------------------------
	def imageZoomIn(self, event):
		for c in self.canvas:
			c.zoomIn(event)
		return "break"

	# ----------------------------------------------------------------------
	def imageZoomOut(self, event):
		for c in self.canvas:
			c.zoomOut(event)
		return "break"

	# ----------------------------------------------------------------------
	def imageZoomOne(self, event):
		for c in self.canvas:
			c.zoomOne(event)
		return "break"

	# ----------------------------------------------------------------------
	def imageScanDrag(self, event):
		for c in self.canvas:
			c.scanDrag(event)
		return "break"

	# ----------------------------------------------------------------------
	def imageScanRelease(self, event):
		for c in self.canvas:
			c.scanRelease(event)
		return "break"

	# ----------------------------------------------------------------------
	# Select a slice
	# ----------------------------------------------------------------------
	def sliceChange(self, event=None):
		self.frame = self.sliceScale.get()
		self.updateImage(self.frame)

	# ----------------------------------------------------------------------
	# Select a view
	# ----------------------------------------------------------------------
	def viewChange(self,view, event=None):
		self.view = view
		self.updateImage(self.frame)

	# ----------------------------------------------------------------------
	# Convert dicom to image and show it on canvas
	# ----------------------------------------------------------------------
	def updateImage(self,event=None):
		"""Convert dicom to mage and show it on canvas"""

		if self.arrayDicom is None: return
		di = self.page.selectedDicom()

		if di is None:
			self.flair.notify("No DICOM selected",
				"Please select a DICOM set first")
			return

		try:
			if di.modality != "CT":# and  di.modality != "RTSTRUCT":
				self.flair.notify("Wrong DICOM selected",
				"Please select a CT or RTSTRUCT set first")
				return
		except:
			self.flair.notify("No DICOM selected",
				"Please select a CT or RTSTRUCT set first")
			return

		try:
			self.dataset = dicom.read_file(
						os.path.join(self.page.dicominfo["directory"], di.files[0][0]),
						defer_size=256,
						force=True)
		except dicom.filereader.InvalidDicomError:
			self.dataset = None
			return

		if self.dataset is None:
			self.canvas[0].setImage()
			return

		try:	window = int(self.window.get())
		except: window = None

		try:	center = int(self.center.get())
		except: center = None

		if window == 0: window = None

		try:	palmin = float(self.palmin.get())
		except: palmin = 0

		try:	palmax = float(self.palmax.get())
		except: palmax = _PALETTE_MAX

		if palmax<=palmin:
			palmax = palmin+_PALETTE_RES
			self.palmax.set(palmax)

		try:	palmult = int(self.palmult.get())
		except: palmult = 10

		paletteLev = (palmin,palmax,palmult)

		#-----------------------------------------------------------------------------
		#Setting the DICOM slice on 3 canvas

		images = None

		try:
			if self.view == "Axial":
				images = self.views(self.view,window,center,self.arrayDicom.shape[2]-1,paletteLev)
			elif self.view == "Coronal":
				images = self.views(self.view,window,center,self.arrayDicom.shape[0]-1,paletteLev)
			elif self.view == "Sagittal":
				images = self.views(self.view,window,center,self.arrayDicom.shape[1]-1,paletteLev)

		except:
			tkFlair.addException()
			self.flair.notify("Cannot display dicom",
					sys.exc_info()[1],
					tkFlair.NOTIFY_WARNING)

			self.canvas[0].setImage()
			images = None
			return

		self.canvas[0].setImage(images[0])
		self.canvas[1].setImage(images[1])
		self.canvas[2].setImage(images[2])

	# ----------------------------------------------------------------------
	def reset(self):
		return

	# ----------------------------------------------------------------------
	def release(self, event):
		pass

	# ----------------------------------------------------------------------
	def downZ(self, event=None):
		self.sliceScale.set(max(0,self.frame-1))

	# ----------------------------------------------------------------------
	def upZ(self, event=None):
		self.sliceScale.set(min(self.sliceScale["to"]-1,self.frame+1))

	# ----------------------------------------------------------------------
	def levelChanged(self, s):
		self.updateImage(self.frame)

	# ----------------------------------------------------------------------
	def paletteLevel(self, event=None):
		self.updateImage(self.frame)

	# ----------------------------------------------------------------------
	def roiSelect(self, event=None):
		self.updateImage(self.frame)

	#-----------------------------------------------------------------------
	# Set the views on each canvas
	#-----------------------------------------------------------------------
	def views(self,view,window,center,sliceLen,paletteLev):
		"""Set the view on each canvas"""

		# Check the slice frame
		if self.frame>sliceLen:
			self.frame = sliceLen
		self.sliceScale.config(to_=sliceLen)

		di = self.page.selectedDicom()
		size=(128,128)

		# Set the image and the size
		if view   ==   "Axial":
			im   = self.arrayDicom[:,:,self.frame]
			size = [self.CT_size[1],self.CT_size[0]]
		elif view ==   "Coronal":
			im   = self.arrayDicom[self.frame,:,:]
			im   = im.copy()
			im   = numpy.rot90(im,1)
			size = (self.CT_size[0],self.CT_size[2])
		elif view ==   "Sagittal":
			im   = self.arrayDicom[:,self.frame,:]
			im   = im.copy()
			im   = numpy.rot90(im,1)
			size = (self.CT_size[1],self.CT_size[2])
		else:
			return None

		dicom =		      dicom2image(self.dataset,im,window,center,size)
		#TODO: rtview
		rtview=		      dicom#rtstruct2image(self.rtstruct,self.dataset,dicom,sop, rois, window,center,self.view)
		Image1=Image2=Image3= rtview

		#TODO better interpolation between slices
		if self.PlanFile:
			if view   ==   "Axial":
				zaxis = ((self.frame*self.CT_bin[2])-self.Plan_low_pos[2])/self.Plan_bin[2]
				if zaxis<0:
					zaxis=0
				elif zaxis>=self.Plan_dim[2]:
					zaxis=self.Plan_dim[2]-1
				im   = self.PlanArray[:,:,int(round(zaxis))]*1.0
				bnnSize= (self.Plan_size[0],self.Plan_size[1])
				bnnPos = (self.Plan_low_pos[0], self.Plan_low_pos[1])
				imP   = im.copy()
				imP   = numpy.swapaxes(imP,0,1)

			elif view ==   "Coronal":
				zaxis = ((self.frame*self.CT_bin[1])-self.Plan_low_pos[1])/self.Plan_bin[1]
				if zaxis<0:
					zaxis=0
				elif zaxis>=self.Plan_dim[1]:
					zaxis=self.Plan_dim[1]-1
				im    = self.PlanArray[:,int(round(zaxis)),:]*1.0
				imP   = im.copy()
				imP   = numpy.rot90(imP,1)
				bnnSize= (self.Plan_size[0],self.Plan_size[2])
				bnnPos = (self.Plan_low_pos[0], self.Plan_high_pos[2])

			elif view ==   "Sagittal":
				zaxis = ((self.frame*self.CT_bin[0])-self.Plan_low_pos[0])/self.Plan_bin[0]
				if zaxis<0:
					zaxis=0
				elif zaxis>=self.Plan_dim[0]:
					zaxis=self.Plan_dim[0]-1
				im  = self.PlanArray[int(round(zaxis)),:,:]*1.0
				imP   = im.copy()
				imP   = numpy.rot90(imP,1)
				bnnSize= (self.Plan_size[1],self.Plan_size[2])
				bnnPos = (self.Plan_low_pos[1], self.Plan_high_pos[2])
			else:
				return None
#			 print bnnPos, bnnSize, self.Plan_size, self.Plan_low_pos, self.Plan_high_pos
			# As a result of bug mentioned
			# http://stackoverflow.com/questions/10854903/
			# what-is-causing-dimension-dependent-attributeerror-in-pil-fromarray-function
			pal_imP = get_palette_value(imP,paletteLev)
			Image1	= bnn2Image(rtview,pal_imP,size,bnnSize,bnnPos,paletteLev)

		if self.CalcFile:
			if view   ==   "Axial":
				zaxis = ((self.frame*self.CT_bin[2])-self.Calc_low_pos[2])/self.Calc_bin[2]
				if zaxis<0:
					zaxis=0
				elif zaxis>=self.Calc_dim[2]:
					zaxis=self.Calc_dim[2]-1
				im    = self.CalcArray[:,:,int(round(zaxis))]*1.0
				bnnSize= (self.Calc_size[0],self.Calc_size[1])
				bnnPos = (self.Calc_low_pos[0], self.Calc_low_pos[1])
				imC   = im.copy()
				imC   = numpy.swapaxes(imC,0,1)

			elif view ==   "Coronal":
				zaxis = ((self.frame*self.CT_bin[1])-self.Calc_low_pos[1])/self.Calc_bin[1]
				if zaxis<0:
					zaxis=0
				elif zaxis>=self.Calc_dim[1]:
					zaxis=self.Calc_dim[1]-1
				im    = self.CalcArray[:,int(round(zaxis)),:]*1.0
				imC   = im.copy()
				imC   = numpy.rot90(imC,1)
				bnnSize= (self.Calc_size[0],self.Calc_size[2])
				bnnPos = (self.Calc_low_pos[0], self.Calc_high_pos[2])
			elif view ==   "Sagittal":
				zaxis = ((self.frame*self.CT_bin[0])-self.Calc_low_pos[0])/self.Calc_bin[0]
				if zaxis<0:
					zaxis=0
				elif zaxis>=self.Calc_dim[0]:
					zaxis=self.Calc_dim[0]-1
				im = self.CalcArray[int(round(zaxis)),:,:]*1.0
				imC   = im.copy()
				imC   = numpy.rot90(imC,1)
				bnnSize= (self.Calc_size[1],self.Calc_size[2])
				bnnPos = (self.Calc_low_pos[1], self.Calc_high_pos[2])
			else:
				return None
			pal_imC =get_palette_value(imC,paletteLev)
			Image2	= bnn2Image(rtview,pal_imC,size,bnnSize,bnnPos,paletteLev)

		if self.PlanFile and self.CalcFile:
			im     = get_dif_palette_value((imC-imP), paletteLev)
			im     = im.copy()
			Image3 = bnn2Image(rtview,im,size,bnnSize,bnnPos,paletteLev,dif=True)

		return (Image1,Image2,Image3)

	#-----------------------------------------------------------------------
	# Reading the TPS and FLUKA dose files variables
	#-----------------------------------------------------------------------

	def _readUsrbin(self,filename):
		""" read Usrbin variables """
		usrbin=Data.Usrbin(filename)
		BinMatrix =Data.unpackArray(usrbin.readData(0)) #there should be only one detector or the first one

		Dim	=[int(usrbin.detector[0].nx),int(usrbin.detector[0].ny),int(usrbin.detector[0].nz)]
		BinSize =[float(usrbin.detector[0].dx)*10,float(usrbin.detector[0].dy)*10,float(usrbin.detector[0].dz)*10]
		Size	=[int(round(BinSize[0]*Dim[0])),int(round(BinSize[1]*Dim[1])),int(round(BinSize[2]*Dim[2]))]
		BinLowPos  =[int(round((usrbin.detector[0].xlow)*10)-self.CT_low_pos[0]),
			     int(round((usrbin.detector[0].ylow)*10)-self.CT_low_pos[1]),
			     int(round((usrbin.detector[0].zlow)*10)-self.CT_low_pos[2])]
		BinHighPos  =[-int(round((usrbin.detector[0].xlow)*10+float(BinSize[0]*Dim[0]))-self.CT_high_pos[0]),
			      -int(round((usrbin.detector[0].ylow)*10+float(BinSize[1]*Dim[1]))-self.CT_high_pos[1]),
			      -int(round((usrbin.detector[0].zlow)*10+float(BinSize[2]*Dim[2]))-self.CT_high_pos[2])]

		return (BinMatrix,Dim, BinSize,Size, BinLowPos, BinHighPos)

	def _planDose(self,filename):
		"""PlanDose variables"""

		BinMatrix,self.Plan_dim, self.Plan_bin, self.Plan_size, self.Plan_low_pos, self.Plan_high_pos = self._readUsrbin(filename)
		self.PlanArray=numpy.reshape(BinMatrix,(self.Plan_dim),order="F")
		self.PlanFile=True

	def _calcDose(self,filename):
		"" "CalcDose variables"""

		BinMatrix,self.Calc_dim, self.Calc_bin, self.Calc_size, self.Calc_low_pos, self.Calc_high_pos = self._readUsrbin(filename)
		self.CalcArray=numpy.reshape(BinMatrix,(self.Calc_dim),order="F")
		self.CalcFile=True

	#-----------------------------------------------------------------------
	# Reading the TPS and FLUKA dose files in usrbin format
	#------------------------------------------------------------------------
	#TODO: catch an errors

	def loadPDose(self, event=None):
		"""load Prescribed Dose from bnn File"""

		fn=bFileDialog.askopenfilename(master=self,
			title="Import Dose from USRBIN file",
			filetypes=[("Fluka USRBIN files",("*.bnn")),
				("All","*")],
			initialdir=self.lastPath,
			initialfile=self.PDose["text"])

		if len(fn) == 0: return

		self.lastPath = os.path.abspath(os.path.dirname(fn))
		relfn = self.project.relativePath(fn)
		self.PDose["text"] = relfn
		self._planDose(self.PDose["text"])

	#-------------------------------------------------------------------------
	def loadCDose(self):
		"""load Calculated Dose from bnn File"""

		fn=bFileDialog.askopenfilename(master=self,
			title="Import Dose from USRBIN file",
			filetypes=[("Fluka USRBIN files",("*.bnn")),
				("All","*")],
			initialdir=self.lastPath,
			initialfile=self.CDose["text"])

		if len(fn) == 0: return

		self.lastPath = os.path.abspath(os.path.dirname(fn))
		relfn = self.project.relativePath(fn)
		self.CDose["text"] = relfn
		self._calcDose(self.CDose["text"])

#------------------------------------------------------------------------------
#Rescale preparing for palette
#-----------------------------------------------------------------------------
def get_palette_value(im, paletteLev):
	"""Apply the Palette Look-Up Table for the given data and min/max  value."""
	if not have_numpy:
		raise ImportError("Numpy is not available. See http://numpy.scipy.org/ to download and install")
	pal_window=paletteLev[1]-paletteLev[0]

	# Fix for numpy 1.9.2 error: all arrays must have the same dimension
	below = (im <= paletteLev[0])
	above = (im >  paletteLev[1])
	totlist = numpy.logical_or.reduce([below,above], axis=0)

	return numpy.piecewise(im,
		[below, above, ~totlist],
		[0,255, lambda im: ((im-paletteLev[0])/(pal_window))*256])


#	return numpy.piecewise(im,[(im <= paletteLev[0]),
#					 (im > paletteLev[1])],
#				   [0,255, lambda im: ((im-paletteLev[0])/(pal_window))*256])

#------------------------------------------------------------------------------
#Rescale preparing for palette diff
#-----------------------------------------------------------------------------
def get_dif_palette_value(im, paletteLev):
	"""Apply the Palette LUT for the given diff data concerning the multiplication factor."""
	if not have_numpy:
		raise ImportError("Numpy is not available. See http://numpy.scipy.org/ to download and install")
	pal_window=(paletteLev[1]-paletteLev[0])/paletteLev[2]

	# Fix for numpy 1.9.2 error: all arrays must have the same dimension
	below = (im <= -paletteLev[1]/paletteLev[2])
	above = (im >=	paletteLev[1]/paletteLev[2])
	totlist = numpy.logical_or.reduce([below,above], axis=0)

	return numpy.piecewise(im,
		[below, above, ~totlist],
		[0,255, lambda im: ( im/(pal_window)*128+127 )])

#	return numpy.piecewise(im,[(im <= -paletteLev[1]/paletteLev[2]),
#					 (im >= paletteLev[1]/paletteLev[2])],
#				     [0,255, lambda im: (im/(pal_window)*128+127)])

#-------------------------------------------------------------------------------
# Rescale using the Window, Center from Look Up Table
#-------------------------------------------------------------------------------
def get_LUT_value(data, window, level):
	"""Apply the RGB Look-Up Table for the given data and window/level value."""
	if not have_numpy:
		raise ImportError("Numpy is not available. See http://numpy.scipy.org/ to download and install")

	window = float(window)
	level  = float(level)

	# Fix for numpy 1.9.2 error: all arrays must have the same dimension
	below = (data <= (level-0.5-(window-1.0)/2.0))
	above = (data >  (level-0.5+(window-1.0)/2.0))
	totlist = numpy.logical_or.reduce([below,above], axis=0)

	return numpy.piecewise(data,
		[below, above, ~totlist],
		[0, 255, lambda data: ((data-(level-0.5))/(window-1.0)+0.5)*255])

#-------------------------------------------------------------------------------
# Convert a DICOM slice to Tk photoimage
#-------------------------------------------------------------------------------
def dicom2image(dataset, im,  window=None, center=None, size=None):
	if 'PixelData' not in dataset:
		raise TypeError("Cannot show image -- DICOM dataset does not have pixel data")

	# Multi frame image
	if 'NumberOfFrames' in dataset:
		try:
			frame = dataset.PerFrameFunctionalGroupsSequence[0]
		except:
			frame = None
		if window:
			window = float(window)
		else:
			try: window = float(frame.FrameVOILUTSequence[0].WindowWidth)
			except: window = 0
		if center:
			center = float(center)
		else:
			try: center = float(frame.FrameVOILUTSequence[0].WindowCenter)
			except: center = 0
		try:
			image = get_LUT_value(im, window, center)
		except IndexError:
			return

		# Convert mode to L since LUT has only 256 values
		# http://www.pythonware.com/library/pil/handbook/image.htm
		# when image is not in int16 (e.g. uint16) format the PIL convert doesn't work
		image = Image.fromarray(image.astype(numpy.int16)).convert('L')
	# can only apply LUT if these values exist
	elif 'WindowWidth' not in dataset or 'WindowCenter' not in dataset:
		bits = dataset.BitsAllocated
		samples = dataset.SamplesPerPixel
		if bits == 8 and samples == 1:
			mode = "L"
		elif bits == 8 and samples == 3:
			mode = "RGB"
		elif bits == 16:
			# not sure about this -- PIL source says is 'experimental'
			# and no documentation. Also, should bytes swap depending
			# on endian of file and system??
			mode = "I;16"
		else:
			raise TypeError("Don't know PIL mode for %d BitsAllocated and %d SamplesPerPixel" \
					% (bits, samples))

		# PIL size = (width, height)
		size = (dataset.Columns, dataset.Rows)

		# Recommended to specify all details by http://www.pythonware.com/library/pil/handbook/image.htm
		image = Image.frombuffer(mode, size, dataset.PixelData, "raw", mode, 0, 1).convert('L')

	# can only apply LUT if these values exist
	else:
		window = window or dataset.WindowWidth
		center = center or dataset.WindowCenter
		try:	window = float(window)
		except: window = float(window[-1])
		try:	center = float(center)
		except: center = float(center[-1])
		image = get_LUT_value(im, window, center)
		# Convert mode to L since LUT has only 256 values
		# http://www.pythonware.com/library/pil/handbook/image.htm
		# when image is not in int16 (e.g. uint16) format the PIL convert doesn't work
		image = Image.fromarray(image.astype(numpy.int16)).convert('L')

	if size is not None:
		image= image.resize((size),resample=Image.BILINEAR)

	return image

#------------------------------------------------------------------------------
#bnn to image
#------------------------------------------------------------------------------
def bnn2Image(dicom,bnnIm,size, bnnSize, bnnPos,paletteLev, dif=False):

	pal_window  = float(paletteLev[1]-paletteLev[0])

	if dif == True:
		palette = Palette.PILPalette("Differ")
	else:
		palette = Palette.PILPalette("FLUKA")
	val_tab=[]
	dummy_tab=[]

	# PALETTE resolution
	# TODO to be corrected
	r = 1.0
	steps =0
	while(True):
		steps = round((paletteLev[1]-paletteLev[0])/r)
		if dif == True:
			steps*=2
		if (steps <10):
			r=r/10.0
		else:
			if (steps >=95):
				r=r*10.0
			elif (steps >=20):
				r = r*5.0
			else:
				break

	if dif==True:
		min_val=0
		for val in numpy.arange(min_val,paletteLev[1]/paletteLev[2],r/paletteLev[2]):
			if (abs(val)>paletteLev[1]/paletteLev[2]):
			    continue
			val_tab.append(((1,size[1]*(0.5-(round(val,2)/(2*pal_window/paletteLev[2])))),str(round(val,2))))
			if val==min_val:
			    continue
			val_tab.insert(0,((1,size[1]*(0.5+(round(val,2)/(2*pal_window/paletteLev[2])))),str(-round(val,2))))
	else:
		min_val = float(int(paletteLev[0] / r)*r)
		for val in numpy.arange(min_val,paletteLev[1],r):
			if ((val<paletteLev[0]) or (val>paletteLev[1])):
				continue
			val_tab.append(((1,size[1]*(1-(round(val,2)- paletteLev[0])/pal_window)),str(round(val,2))))

	bnnIm	 = bnnIm.copy()

	#arrange palete table on the right side
	pal    = [numpy.arange(255,-1,-1)]*_PALETTE_SIZE
	pal    = numpy.swapaxes(pal,0,1)
	pal    = pal.copy()
	pal    = Image.fromarray(pal.astype(numpy.int16)).convert('L')
	pal    = pal.resize((_PALETTE_SIZE,size[1]))
	pal.putpalette(palette)
	pal    = pal.convert("RGB")
	drawpal   = ImageDraw.Draw(pal)
	drawpal.rectangle([0,0,_PALETTE_SIZE,size[1]], outline="black")

	text   = Image.new("L",(_PALETTE_SIZE,size[1]), "white")
	drawtext   = ImageDraw.Draw(text)

	for value in val_tab:
		drawtext.text(value[0],value[1],fill="black")
		drawpal.line((value[0][0],value[0][1],value[0][0]+_PALETTE_SIZE,value[0][1]),fill="black")

	#arrange dicom and bnn results based on the palette

	bnnIm	 = Image.fromarray(bnnIm.astype(numpy.int16)).convert('L')
	bnnIm.putpalette(palette)
	bnnImage = bnnIm.resize((bnnSize),resample=Image.BILINEAR)

	im   = Image.new("L",(size))
	im.putpalette(Palette.PILPalette("FLUKA"))
	image=im.convert("RGB")
	image.paste(bnnImage,(bnnPos))
	mask = image.convert("L")

	image.paste(dicom,mask)

	#combine dicom image + palette
	new_image = Image.new("RGB",(size[0]+2*_PALETTE_SIZE,size[1]))
	new_image.paste(image,(0,0))
	new_image.paste(pal,(size[0],0))
	new_image.paste(text,(size[0]+_PALETTE_SIZE,0))

	return new_image

#-------------------------------------------------------------------------------
# TODO: Draw an RT structure on image
#-------------------------------------------------------------------------------
def rtstruct2image(dataset, ds, image,	sop, roi, window, center,view):
	try:
		image = image.convert('RGB')
		x0, y0, z0 = list(map(float,ds.ImagePositionPatient))
		dx, dy	   = list(map(float,ds.PixelSpacing))
		dz	   = 2
	except:
		image = Image.new("RGB", (512,512))
		x0, y0 = 0., 0.
		dx, dy = 1., 1.

	# Find coordinates for transformation of ROI

	draw  = ImageDraw.Draw(image)
#	for roi in rois:
		# XXX damn it ROI sequence is not incremental can be anything!!
		# Find the sequence idx from the roi
	for sequence in dataset.ROIContourSequence:
		if int(sequence.ReferencedROINumber) == int(roi):
			break
	else:
		return image

		# Find contour if any based on SOP
	for contour in sequence.ContourSequence:
		if contour.ContourImageSequence[0].ReferencedSOPInstanceUID == sop:
			break
	else:
		return image

	# Color
	try:
		color = sequence.ROIDisplayColor
		color = int(color[2])<<16 | int(color[1])<<8 | int(color[0])
	except AttributeError:
		color = 128

	# Ignore non-planar (for the moment)
	if contour.ContourGeometricType in ("CLOSED_PLANAR", "POINT"):
		npoints = contour.NumberOfContourPoints
		coords = list(map(float, contour.ContourData))
	 #	 print coords[2]
		if npoints == 1:
			iX = (coords[0]-x0)/dx
			iY = (coords[1]-y0)/dy
			draw.line((iX,iY-5,iX,iY+5), fill=color)
			draw.line((iX-5,iY,iX+5,iY), fill=color)

		else:
			X = coords[::3]
			Y = coords[1::3]
			Z = coords[2::3]
	   #		 print Z
			iX = [int((x-x0)/dx) for x in X]
			iY = [int((y-y0)/dy) for y in Y]
			iZ = [int((z-z0)/dz) for z in Z]
			if view =="Axial":
				draw.polygon(list(zip(iX,iY)), outline=color)
			if view=="Coronal":
				draw.polygon(list(zip(iX,iZ)), outline=color)
	#elif contour.ContourGeometricType == "POINT":
	#	say("Contour:",contour)

	del draw
	return image

#-------------------------------------------------------------------------------
def _cmp(column, x, y):
	if column==0:
		return cmp(x[column], y[column])
	else:
		return cmp(float(x[column]), float(y[column]))
