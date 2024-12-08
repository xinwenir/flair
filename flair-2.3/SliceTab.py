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

from tkinter import *

import os
import sys
import math
import Ribbon
import tkFlair
import tkExtra
import Palette
import ZoomImage
import FlairRibbon
import DicomBrowser

try:
	import pydicom as dicom
except ImportError:
	import dicom

try:
	from PIL import Image
	from PIL import ImageDraw
	have_PIL = True
except:
	have_PIL = False

try:
	import numpy
	have_numpy = True
except ImportError:
	have_numpy = False

from Unicode import RIGHTWARDS_ARROW

_MARGIN      = 5
_NEW         = 0
_TOPLEFT     = 1
_TOP         = 2
_TOPRIGHT    = 3
_RIGHT       = 4
_BOTTOMRIGHT = 5
_BOTTOM      = 6
_BOTTOMLEFT  = 7
_LEFT        = 9

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
		[0, 255, lambda x: ((x-(level-0.5))/(window-1.0)+0.5)*255])

#-------------------------------------------------------------------------------
# Convert a DICOM slice to Tk photoimage
#-------------------------------------------------------------------------------
def dicom2image(dataset, window=None, center=None, z=0):
	#print
	#print "="*50
	#import traceback; traceback.print_stack()

	if 'PixelData' not in dataset:
		raise TypeError("Cannot show image -- DICOM dataset does not have pixel data")

	# Multi frame image
	if 'NumberOfFrames' in dataset:
		try:
			frame = dataset.PerFrameFunctionalGroupsSequence[0]
		except:
			frame = None

		try: window = window or dataset.WindowWidth
		except: pass
		try: center = center or dataset.WindowCenter
		except: pass

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
			image = get_LUT_value(dataset.pixel_array[z], window, center)
		except IndexError:
			return

		# Convert mode to L since LUT has only 256 values
		# http://www.pythonware.com/library/pil/handbook/image.htm
		# when image is not in int16 (e.g. uint16) format the PIL convert doesn't work
		image = Image.fromarray(image.astype(numpy.int16)).convert('L')

		# Testing palette for DOSE distribution
		#image2 = Image.new("P", image.size)
		#image2.putpalette(Palette.PILPalette("FLUKA"))
		#image2.putdata(image.tostring())
		#image = image2

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

	else:
		window = window or dataset.WindowWidth
		center = center or dataset.WindowCenter
		try:    window = float(window)
		except: window = float(window[-1])
		try:    center = float(center)
		except: center = float(center[-1])
		image = get_LUT_value(dataset.pixel_array, window, center)
		# Convert mode to L since LUT has only 256 values
		# http://www.pythonware.com/library/pil/handbook/image.htm
		# when image is not in int16 (e.g. uint16) format the PIL convert doesn't work
		image = Image.fromarray(image.astype(numpy.int16)).convert('L')

	return image

#-------------------------------------------------------------------------------
# Draw an RT structure on image
#-------------------------------------------------------------------------------
def rtstruct2image(dataset, ctfilename, z, sop, rois, window, center):
	try:
		ds = dicom.read_file(ctfilename, force=True)
		image = dicom2image(ds, window, center).convert('RGB')
		x0, y0, z0 = list(map(float,ds.ImagePositionPatient))
		dx, dy     = list(map(float,ds.PixelSpacing))
	except:
		image = Image.new("RGB", (512,512))
		x0, y0 = 0., 0.
		dx, dy, dz = 1., 1., 1.

	try:
		dz = ds.SliceThickness
	except:
		dz = 1.

	# Find coordinates for transformation of ROI
	draw = ImageDraw.Draw(image)
	for roi in rois:
		# XXX damn it ROI sequence is not incremental can be anything!!
		# Find the sequence idx from the roi
		for sequence in dataset.ROIContourSequence:
			if int(sequence.ReferencedROINumber) == roi:
				break
		else:
			continue

		# Color
		try:
			color = sequence.ROIDisplayColor
			color = int(color[2])<<16 | int(color[1])<<8 | int(color[0])
		except AttributeError:
			color = 128

		# Find contour if any based on SOP
		for contour in sequence.ContourSequence:
			try:
				if contour.ContourImageSequence[0].ReferencedSOPInstanceUID == sop and \
				   contour.ContourGeometricType == "CLOSED_PLANAR":
					npoints = contour.NumberOfContourPoints
					coords = list(map(float, contour.ContourData))
					X = coords[::3]
					Y = coords[1::3]
					iX = [int((x-x0)/dx) for x in X]
					iY = [int((y-y0)/dy) for y in Y]
					draw.polygon(list(zip(iX,iY)), outline=color)
			except AttributeError:
				if contour.ContourGeometricType == "POINT" and \
				   abs(float(contour.ContourData[2]) - z) <= dz/2:
					npoints = contour.NumberOfContourPoints
					coords = list(map(float, contour.ContourData))
					if npoints == 1:
						iX = (coords[0]-x0)/dx
						iY = (coords[1]-y0)/dy
						draw.line((iX,iY-5,iX,iY+5), fill=color)
						draw.line((iX-5,iY,iX+5,iY), fill=color)
	del draw
	return image

#-------------------------------------------------------------------------------
def _cmp(column, x, y):
	if column==0:
		return cmp(x[column], y[column])
	else:
		return cmp(float(x[column]), float(y[column]))

#===============================================================================
# Slice Tab
#===============================================================================
class SliceTab(FlairRibbon.FlairTab):
	#----------------------------------------------------------------------
	def __init__(self, master, page, **kw):
		FlairRibbon.FlairTab.__init__(self, master, page, **kw)

		# -------------- Slices ---------------
		pane = PanedWindow(self, orient=HORIZONTAL, opaqueresize=1)
		pane.pack(fill=BOTH, expand=YES)

		labelframe = Frame(pane)
		pane.add(labelframe)

		self.sliceList = tkExtra.MultiListbox(labelframe,
					(("Name", 16, None),
					 ("z",    10, None)))
		self.sliceList.pack(side=LEFT, expand=YES, fill=BOTH)
		self.sliceList.bind("<Double-1>",        self.viewer)
		self.sliceList.bind("<<ListboxSelect>>", self.sliceSelect)
		self.sliceList.bind("<Control-Key-f>", lambda e,s=self : s.findString.focus_set())
		self.sliceList.sortAssist = None

		# -------------- Tabs -----
		self.tabPage = tkExtra.TabPageSet(pane, pageNames=[
					("Viewer",  tkFlair.icons["dicom"]),
					("Browser", tkFlair.icons["new"])])
		self.tabPage.bind("<<ChangePage>>", self.changePage)
		pane.add(self.tabPage)

		# -----
		frame = self.tabPage.page("Viewer")

		# --------------- Canvas -----------
		row,col = 0,0
		self.canvas = ZoomImage.ZoomImage(frame, takefocus=True)
		self.canvas.grid(row=row, column=col, columnspan=2, sticky=NSEW)

		self.canvas.bindMotion(wheel=False)
		self.canvas.bindArrows()
		self.canvas.bindZoom()
		self.canvas.bind('<<Zoom>>',          self.drawRectangle)
		self.canvas.bind('<Control-Button-4>',self.downZ)
		self.canvas.bind('<Control-Button-5>',self.upZ)
		self.canvas.bind('<Button-1>',        self.click)
		self.canvas.bind('<B1-Motion>',       self.drag)
		self.canvas.bind('<ButtonRelease-1>', self.release)
		self.canvas.bind('<Motion>',          self.motion)
		self.canvas.bind('<Next>',            self.downZ)
		self.canvas.bind('<Prior>',           self.upZ)

		col += 2
		sbv = Scrollbar(frame, orient=VERTICAL, takefocus=False,
				command=self.canvas.yview)
		sbv.grid(row=row, column=col, sticky=NS)

		row,col = row+1,0
		sbh = Scrollbar(frame, orient=HORIZONTAL, takefocus=False,
				command=self.canvas.xview)
		sbh.grid(row=row, column=col, columnspan=2, sticky=EW)

		self.canvas.config(xscrollcommand=sbh.set, yscrollcommand=sbv.set)

		# -------------
		row,col = row+1,0
		f = Frame(frame)
		f.grid(row=row, column=col, columnspan=3, sticky=NSEW)

		Label(f, text="x:").pack(side=LEFT)
		self.xpixel = Label(f, width=10, foreground="DarkBlue", anchor=W)
		self.xpixel.pack(side=LEFT)

		Label(f, text=" y:").pack(side=LEFT)
		self.ypixel = Label(f, width=10, foreground="DarkBlue", anchor=W)
		self.ypixel.pack(side=LEFT)

		Label(f, text=" z:").pack(side=LEFT)
		self.zpixel = Label(f, width=10, foreground="DarkBlue", anchor=W)
		self.zpixel.pack(side=LEFT)

		Label(f, text=" Dicom:").pack(side=LEFT)
		self.dicompixel = Label(f, width=10, foreground="DarkBlue", anchor=W)
		self.dicompixel.pack(side=LEFT)

		Label(f, text=" Value:").pack(side=LEFT)
		self.valuepixel = Label(f, width=10, foreground="DarkBlue", anchor=W)
		self.valuepixel.pack(side=LEFT)

		# -------------
		row,col = row+1,0
		lf = Frame(frame)
		lf.grid(row=row, column=col, columnspan=3, sticky=NSEW)

		r,c = 0,0
		Label(lf, text="Window Center:").grid(row=r, column=c, sticky=S+E)

		c += 1
		self.center = Scale(lf, from_=0, to_=5000,
					showvalue=0,
					orient=HORIZONTAL,
					command=self.levelChanged)
		self.center.grid(row=r, column=c, sticky=EW)
		tkExtra.Balloon.set(self.center, "Center or Level of values to display")

		# ---
		c += 1
		Label(lf, text="Width:").grid(row=r, column=c, sticky=S+E)

		c += 1
		self.window = Scale(lf, from_=0, to_=5000,
					showvalue=0,
					orient=HORIZONTAL,
					command=self.levelChanged)
		self.window.grid(row=r, column=c, sticky=EW)
		tkExtra.Balloon.set(self.center, "Window around level/center values to display")
		lf.grid_columnconfigure(c, weight=1)

		# -------------
		r,c = r+1,0
		Label(lf, text="Slice:").grid(row=r, column=c, sticky=S+E)

		c += 1
		self.sliceScale = Scale(lf, from_=0, to_=100,
					orient=HORIZONTAL,
					command=self.sliceChange)
		self.sliceScale.grid(row=r, column=c, columnspan=3, sticky=EW)
		tkExtra.Balloon.set(self.sliceScale, "Slice/Frame to display")

		lf.grid_columnconfigure(1, weight=1)
		lf.grid_columnconfigure(3, weight=1)

		# -------------
		frame.grid_rowconfigure(0, weight=1)
		frame.grid_columnconfigure(1, weight=1)
#		frame.grid_columnconfigure(3, weight=1)

		# -----
		frame = self.tabPage.page("Browser")
		self.browser = DicomBrowser.DicomBrowser(frame)
		self.browser.pack(side=LEFT, expand=YES, fill=BOTH)

		self.tabPage.changePage("Viewer")

		self.grid_columnconfigure(1, weight=1)
		self._photo  = None
		self._rect   = None

		self.dataset = None	# current loaded dataset
		self.frame   = 0	# Number of frames for scale
		self.CTInfo  = None	# CT project DicomInfo structure related to this dataset

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairTab.createRibbon(self)

		self.groups = {}

		# ========== Viewer ===========
		self.groups["Viewer:Crop"] = group = Ribbon.LabelGroup(self.ribbon, "Crop")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.grid3rows()

		# --- Xcrop ---
		row,column = 0,0
		Label(group.frame, text="Xmin:").grid(row=row, column=column, sticky=E)
		column += 1
		self.xfrom = tkExtra.FloatEntry(group.frame, background="white", width=8)
		self.xfrom.grid(row=row, column=column, sticky=EW)
		self.xfrom.bind("<Return>",   self.drawRectangle)
		self.xfrom.bind("<KP_Enter>", self.drawRectangle)
		self.xfrom.bind("<FocusOut>", self.drawRectangle)
		group.frame.grid_columnconfigure(column, weight=1)

		column += 1
		Label(group.frame, text="Xmax:").grid(row=row, column=column, sticky=E)
		column += 1
		self.xto = tkExtra.FloatEntry(group.frame, background="white", width=8)
		self.xto.grid(row=row, column=column, sticky=EW)
		self.xto.bind("<Return>",   self.drawRectangle)
		self.xto.bind("<KP_Enter>", self.drawRectangle)
		self.xto.bind("<FocusOut>", self.drawRectangle)
		group.frame.grid_columnconfigure(column, weight=1)

		# --- Ycrop ---
		row,column = row+1,0
		Label(group.frame, text="Ymin:").grid(row=row, column=column, sticky=E)
		column += 1
		self.yfrom = tkExtra.FloatEntry(group.frame, background="white", width=8)
		self.yfrom.grid(row=row, column=column, sticky=EW)
		self.yfrom.bind("<Return>",   self.drawRectangle)
		self.yfrom.bind("<KP_Enter>", self.drawRectangle)
		self.yfrom.bind("<FocusOut>", self.drawRectangle)

		column += 1
		Label(group.frame, text="Ymax").grid(row=row, column=column, sticky=E)
		column += 1
		self.yto = tkExtra.FloatEntry(group.frame, background="white", width=8)
		self.yto.grid(row=row, column=column, sticky=EW)
		self.yto.bind("<Return>",   self.drawRectangle)
		self.yto.bind("<KP_Enter>", self.drawRectangle)
		self.yto.bind("<FocusOut>", self.drawRectangle)

		# --- Zcrop ---
		row,column = row+1,0
		Label(group.frame, text="Zmin:").grid(row=row, column=column, sticky=E)
		column += 1
		self.zfrom = tkExtra.FloatEntry(group.frame, background="white", width=8)
		self.zfrom.grid(row=row, column=column, sticky=EW)
		column += 1
		Label(group.frame, text="Zmax:").grid(row=row, column=column, sticky=E)
		column += 1
		self.zto = tkExtra.FloatEntry(group.frame, background="white", width=8)
		self.zto.grid(row=row, column=column, sticky=EW)

		# --- Buttons
		row,column = 0, column+1
		Ribbon.LabelButton(group.frame,
				text="Reset",
				anchor=W,
				command=self.reset,
				background=Ribbon._BACKGROUND).grid(
					row=row, column=column, sticky=EW)
		row += 1
		Ribbon.LabelButton(group.frame,
				text="Slice "+RIGHTWARDS_ARROW+" Zmin",
				anchor=W,
				command=self.setZFrom,
				background=Ribbon._BACKGROUND).grid(
					row=row, column=column, sticky=EW)
		row += 1
		Ribbon.LabelButton(group.frame,
				text="Slice "+RIGHTWARDS_ARROW+" Zmax",
				anchor=W,
				command=self.setZTo,
				background=Ribbon._BACKGROUND).grid(
					row=row, column=column, sticky=EW)

		# ========== Browser ===========
		self.groups["Browser"] = group = Ribbon.LabelGroup(self.ribbon, "Browser")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row = 0,1
		self.findString = tkExtra.LabelEntry(group.frame,
					"Search",
					"DarkGray",
					background="White", #Ribbon._BACKGROUND,
					width=20)
		self.findString.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(self.findString , "Search string")
		self.findString.bind("<Return>",   self.find)
		self.findString.bind("<KP_Enter>", self.find)

		# ---
		col += 1
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
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Case sensitive")

		# ---
		col += 1
		row = 0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["find"],
				text = "Search",
				compound = LEFT,
				command=self.find,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=W)
		tkExtra.Balloon.set(b, "Find cards with string")

		# ---
		row += 1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["filter"],
				text = "Filter",
				compound = LEFT,
				command=self.filterBrowser,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=W)
		tkExtra.Balloon.set(b, "Filter: display only cards matching the search string")

		# ---
		row += 1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clean"],
				text = "Clear",
				compound = LEFT,
				command=self.filterClear,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=W)
		tkExtra.Balloon.set(b, "Filter clear")

		# ---
		col += 1
		row = 0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["viewer32"],
				text="Viewer",
				compound=TOP,
				command=self.viewer,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, rowspan=3, column=col, padx=0, pady=0, sticky=W)
		tkExtra.Balloon.set(b, "Open in internal viewer")

		self.changePage()

	# ----------------------------------------------------------------------
	def changePage(self, event=None):
		tab =  self.tabPage.getActivePage()
		for name,group in list(self.groups.items()):
			if name.startswith(tab):
				group.pack(side=LEFT, fill=Y, padx=0, pady=0)
			else:
				group.pack_forget()

	# ----------------------------------------------------------------------
	def dicom(self):	return self.page.dicom

	# ----------------------------------------------------------------------
	def refresh(self):
		di = self.page.dicominfo
		if di is None:
			self.canvas.setImage()
			self.browser.clear()
			self.drawRectangle()
			self._photo = None
			return

		self.dataset = None

		self.xfrom.set(di.get("xfrom",""))
		self.xto.set(  di.get("xto",""))

		self.yfrom.set(di.get("yfrom",""))
		self.yto.set(  di.get("yto",""))

		self.zfrom.set(di.get("zfrom",""))
		self.zto.set(  di.get("zto",""))

		# Changing the center/window will update the image/crop
		self.window.set(float(di.get("window",0)))
		self.center.set(float(di.get("center",0)))

	# ----------------------------------------------------------------------
	def get(self):
		di = self.page.dicominfo
		if di is None: return

		di["xfrom"] = self.xfrom.get()
		di["xto"]   = self.xto.get()

		di["yfrom"] = self.yfrom.get()
		di["yto"]   = self.yto.get()

		di["zfrom"] = self.zfrom.get()
		di["zto"]   = self.zto.get()

		di["window"] = self.window.get()
		di["center"] = self.center.get()

	# ----------------------------------------------------------------------
	# Return crop image in pixel coordinates
	# ----------------------------------------------------------------------
	def crop(self):
		if self.dataset is None: return

		try:
			xfrom = max(0,min(self.dicom().nx-1,
				int(math.floor((float(self.xfrom.get())-self.dicom().x)/self.dicom().dx))))
		except (ValueError, ZeroDivisionError):
			xfrom = 0

		try:
			xto = max(0,min(self.dicom().nx-1,
				int(math.floor((float(self.xto.get())-self.dicom().x)/self.dicom().dx))))
		except (ValueError, ZeroDivisionError):
			xto = self.dicom().nx-1

		try:
			yfrom = max(0,min(self.dicom().ny-1,
				int(math.floor((float(self.yfrom.get())-self.dicom().y)/self.dicom().dy))))
		except (ValueError, ZeroDivisionError):
			yfrom = 0

		try:
			yto = max(0,min(self.dicom().ny-1,
				int(math.floor((float(self.yto.get())-self.dicom().y)/self.dicom().dy))))
		except (ValueError, ZeroDivisionError):
			yto = self.dicom().ny-1

		if xto < xfrom:
			xfrom, xto = xto, xfrom
			tmp = self.xfrom.get()
			self.xfrom.set(self.xto.get())
			self.xto.set(tmp)

		if yto < yfrom:
			yfrom, yto = yto, yfrom
			tmp = self.yfrom.get()
			self.yfrom.set(self.yto.get())
			self.yto.set(tmp)

		return xfrom, yfrom, xto, yto

	# ----------------------------------------------------------------------
	# @return limits as entered by the user
	# ----------------------------------------------------------------------
	def limits(self):
		try:    xfrom = float(self.xfrom.get())
		except: xfrom = None

		try:    xto   = float(self.xto.get())
		except: xto   = None

		try:    yfrom = float(self.yfrom.get())
		except: yfrom = None

		try:    yto   = float(self.yto.get())
		except: yto   = None

		try:    zfrom = float(self.zfrom.get())
		except: zfrom = None

		try:    zto   = float(self.zto.get())
		except: zto   = None

		return (xfrom, xto, yfrom, yto, zfrom, zto)

	# ----------------------------------------------------------------------
	# Select a new dataset and populate the listbox
	# ----------------------------------------------------------------------
	def activate(self, event=None):
		if self.dataset is not None and self.dataset.SeriesInstanceUID == self.dicom().uid:
			return

		self.clearDatasetSelection()
		di = self.page.selectedDicom()
		if di is None:
			self.flair.notify("No DICOM selected",
				"Please select a DICOM set first")
			return

		# Init dicom with the first slice
		try:
			self.dataset = dicom.read_file(
						os.path.join(self.page.dicominfo["directory"], di.files[0][0]),
						defer_size=256,
						force=True)
		except dicom.filereader.InvalidDicomError:
			self.dataset = None
			return

		if self.dataset:
			self.dicom().init(self.dataset, z=di.files[0][1])

		self.sliceList.delete(0,END)
		self.CTInfo = None

		# Fill list according to modality
		if self.dicom().modality == "RTSTRUCT":
			self.sliceList.labels(["Name","Id"])
			for roi in self.dataset.StructureSetROISequence:
				self.sliceList.insert(END,(str(roi.ROIName), str(roi.ROINumber)))
			self.loadSlice(di.files[0][0])

			# Find the corresponding CT data set
			self.scan4ReferencedCT()
			if self.CTInfo is None:
				self.flair.notify("Referenced CT not found",
					"RT structure is referring to an DICOM image which is not loaded",
					tkFlair.NOTIFY_ERROR)
				return

			# Set scale to the number of SOP UIDs
			self.sliceScale.config(to_=len(self.CTInfo.files)-1)

		else:
			self.sliceList.labels(["Name","z"])
			for fnz in di.files:
				self.sliceList.insert(END,fnz)

		self.sliceList.sort(1)
		self.sliceList.selection_set(0)
		self.sliceSelect()

	# ----------------------------------------------------------------------
#	def deactivate(self):
#		self.get()

	# ----------------------------------------------------------------------
	# For RTSTRUCT return a dictionary with sop and CT filename
	# ----------------------------------------------------------------------
	def scan4ReferencedCT(self):
		self.CTInfo = None
		for cnt in self.dataset.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].ContourImageSequence:
			self.CTInfo = self.page.findDicomFromSOPInstanceUID(cnt.ReferencedSOPInstanceUID)
			if self.CTInfo is not None:
				return

	# ----------------------------------------------------------------------
	# Select a slice
	# ----------------------------------------------------------------------
	def sliceSelect(self, event=None):
		sel = self.sliceList.curselection()
		if not sel: return
		if self.dicom().modality == "RTSTRUCT":
			self.updateImage()
		else:
			if len(sel) == 1:
				active = sel[0]
			else:
				active = self.sliceList.index(ACTIVE)
			self.loadSlice(self.sliceList.get(active)[0])
			self.updateImage()
			if event is not None and self.dicom().slice:
#				if len(sel)>1:
#					self.updateImage()
#				else:
				if len(sel)<=1:
					self.sliceScale.set(active)

	# ----------------------------------------------------------------------
	def sliceChange(self, event=None):
		self.frame = self.sliceScale.get()
		if self.dicom().slice:
			self.sliceList.selection_clear(0,END)
			self.sliceList.activate(self.frame)
			self.sliceList.selection_set(self.frame)
			self.sliceList.see(self.frame)
			self.sliceSelect()

		# update canvas
		self.updateImage()
		#print self.frame, self.dicom().dz, self.dicom().z, self.dicom().x, self.dicom().y
		self.zpixel["text"] = "%.4g"%(float(self.frame)*self.dicom().dz + self.dicom().z)

	# ----------------------------------------------------------------------
	# Load slice
	# ----------------------------------------------------------------------
	def loadSlice(self, filename):
		# remember position
		yview = self.browser.yview()
		self.clearSlice()
		try:
			self.dataset = dicom.read_file(
						os.path.join(self.page.dicominfo["directory"], filename),
						force=True)
		except dicom.filereader.InvalidDicomError:
			return

		# update browser
		self.browser.fill(self.dataset)
		self.browser.yview_moveto(yview[0])
		txt = self.findString.get()
		if txt: self.search(txt)

		self.dicom().init(self.dataset, z=self.page.selectedDicom().files[0][1])
		# Force change to Browser tab if no image is present
		if not self.dicom().hasPixel:
			self.tabPage.changePage("Browser")
			return

		if self.dicom().slice:
			# 2D slice data
			self.sliceScale.config(to_=self.sliceList.size()-1)
		else:
			# 3D multi frame data
			self.sliceScale.config(to_=self.dicom().nz-1)

	# ----------------------------------------------------------------------
	# Convert dicom to image and show it in canvas
	# ----------------------------------------------------------------------
	def updateImage(self, event=None):
		if self.dataset is None:
			self.drawRectangle()
			self.canvas.setImage()
			self._photo = None
			return

		try:    window = int(self.window.get())
		except: window = None
		try:    center = int(self.center.get())
		except: center = None
		if window == 0: window = None

		if self.dicom().modality == "RTSTRUCT":
			sel = list(map(int,self.sliceList.curselection()))
			# XXX WARNING ROI order is not the same with the sequence index!!!
			rois = []
			for i in sel:
				rois.append(int(self.sliceList.lists[1].get(i)))
			if self.CTInfo is None:
				self.flair.notify("Referenced SOP not found",
					"RT structure is referring to an DICOM image which is not loaded",
					tkFlair.NOTIFY_ERROR)
			else:
				fn,z,sop  = self.CTInfo.files[self.frame]
				ct = os.path.join(self.CTInfo["directory"], fn)
				self._photo = rtstruct2image(self.dataset, ct, z, sop, rois, window, center)
		else:
			try:
				self._photo = dicom2image(self.dataset, window, center, self.frame)
			except:
				tkFlair.addException()
				self.flair.notify("Cannot display dicom",
					sys.exc_info()[1],
					tkFlair.NOTIFY_WARNING)

				self.tabPage.changePage("Browser")

				self.canvas.setImage()
				self._photo = None
				return

		self.canvas.setImage(self._photo)
		self.drawRectangle()

		bb = self.canvas.bbox('all')
		if bb is None: return
		self.canvas.configure(scrollregion=bb)

	# ----------------------------------------------------------------------
	def drawRectangle(self, event=None):
		if self.dataset is None: return

		xfrom, yfrom, xto, yto = self.crop()

		xfrom *= self.canvas.zoom
		yfrom *= self.canvas.zoom
		xto   *= self.canvas.zoom
		yto   *= self.canvas.zoom

		if self._rect is None:
			self._rect = self.canvas.create_rectangle((xfrom,yfrom, xto, yto), outline="Orange")
		else:
			self.canvas.coords(self._rect, (xfrom,yfrom, xto, yto))

	# ----------------------------------------------------------------------
	def reset(self):
		self.xfrom.set("")
		self.xto.set("")
		self.yfrom.set("")
		self.yto.set("")
		self.zfrom.set("")
		self.zto.set("")
		self.drawRectangle()

	# ----------------------------------------------------------------------
	def setZFrom(self):
		sel = self.sliceList.curselection()
		if len(sel)==1 and self.dataset:
			try:
				self.zfrom.set(self.dataset.SliceLocation)
			except AttributeError:
				pass
		else:
			self.setZRange()

	# ----------------------------------------------------------------------
	def setZTo(self):
		sel = self.sliceList.curselection()
		if len(sel)==1 and self.dataset:
			try:
				self.zto.set(self.dataset.SliceLocation)
			except AttributeError:
				pass
		else:
			self.setZRange()

	# ----------------------------------------------------------------------
	def setZRange(self):
		sel = self.sliceList.curselection()
		if len(sel)<=1: return

		zmin = 9999999.0
		zmax = -zmin

		for i in sel:
			z = float(self.sliceList.get(i)[1])
			zmin = min(z,zmin)
			zmax = max(z,zmax)

		self.zfrom.set(zmin)
		self.zto.set(zmax)

	# --------------------------------------------------------------------
	def motion(self, event):
		if self._photo  is None or \
		   self.dataset is None or \
		   self.dicom() is None or \
		   self.dicom().modality == "RTSTRUCT":
			return

		i = min(self.dicom().nx-1, max(0, int(self.canvas.canvasx(event.x)/self.canvas.zoom)))
		j = min(self.dicom().ny-1, max(0, int(self.canvas.canvasy(event.y)/self.canvas.zoom)))

		x = float(i)*self.dicom().dx + self.dicom().x
		y = float(j)*self.dicom().dy + self.dicom().y
		z = float(self.frame)*self.dicom().dz + self.dicom().z

		self.xpixel["text"] = "%.4g"%(x)
		self.ypixel["text"] = "%.4g"%(y)
		self.zpixel["text"] = "%.4g"%(z)

		if self.dicom().slice:
			d = self.dataset.pixel_array[j][i]
		else:
			d = self.dataset.pixel_array[self.frame][j][i]

		self.dicompixel["text"] = "%g"%(d)
		self.valuepixel["text"] = "%g"%(self.dicom().m*d + self.dicom().b)

	# --------------------------------------------------------------------
	def click(self, event):
		if self._photo is None: return
		self.canvas.focus_set()
		x = self.canvas.canvasx(event.x)
		y = self.canvas.canvasy(event.y)
		enclosed = self.canvas.find_overlapping(x-5, y-5, x+5, y+5)

		if self._rect in enclosed:
			# Find item to drag
			xfrom, yfrom, xto, yto = self.crop()

			if abs(y - yfrom) <= _MARGIN:	# On top line
				if  abs(x - xfrom) <= _MARGIN:
					self._item = _TOPLEFT
				elif abs(x - xto) <= _MARGIN:
					self._item = _TOPRIGHT
				else:
					self._item = _TOP

			elif abs(y - yto) <= _MARGIN:	# On bottom line
				if  abs(x - xfrom) <= _MARGIN:
					self._item = _BOTTOMLEFT
				elif abs(x - xto) <= _MARGIN:
					self._item = _BOTTOMRIGHT
				else:
					self._item = _BOTTOM

			elif abs(x - xfrom) <= _MARGIN:
				self._item = _LEFT

			else:
				self._item = _RIGHT

		else:
			self._item = _NEW

		self._startx = x
		self._starty = y

	# --------------------------------------------------------------------
	def drag(self, event):
		if self._photo is None: return
		i = min(self.dicom().nx-1, max(0, int(self.canvas.canvasx(event.x)/self.canvas.zoom)))
		j = min(self.dicom().ny-1, max(0, int(self.canvas.canvasy(event.y)/self.canvas.zoom)))
		x = "%.4g"%(float(i)*self.dicom().dx + self.dicom().x)
		y = "%.4g"%(float(j)*self.dicom().dy + self.dicom().y)

		if   self._item == _NEW:
			if abs(self._startx - self.canvas.canvasx(event.x))>3 or \
			   abs(self._starty - self.canvas.canvasy(event.y))>3:
				self.xfrom.set(x)
				self.yfrom.set(y)
				self._item = _BOTTOMRIGHT

		elif self._item == _TOPLEFT:
			self.xfrom.set(x)
			self.yfrom.set(y)

		elif self._item == _TOP:
			self.yfrom.set(y)

		elif self._item == _TOPRIGHT:
			self.xto.set(x)
			self.yfrom.set(y)

		elif self._item == _RIGHT:
			self.xto.set(x)

		elif self._item == _BOTTOMRIGHT:
			self.xto.set(x)
			self.yto.set(y)

		elif self._item == _BOTTOM:
			self.yto.set(y)

		elif self._item == _BOTTOMLEFT:
			self.xfrom.set(x)
			self.yto.set(y)

		elif self._item == _LEFT:
			self.xfrom.set(x)

		if self._item != _NEW: self.drawRectangle()

	# --------------------------------------------------------------------
	def release(self, event):
		pass

	# --------------------------------------------------------------------
	def zoomIn(self, event=None):
		self.canvas.zoomIn(event)
		self.drawRectangle()

	# ----------------------------------------------------------------------
	def zoomOut(self, event=None):
		self.canvas.zoomOut(event)
		self.drawRectangle()

	# ----------------------------------------------------------------------
	def zoomOne(self, event=None):
		self.canvas.zoomImage(1)
		self.drawRectangle()

	# ----------------------------------------------------------------------
	def downZ(self, event=None):
		self.sliceScale.set(max(0,self.frame-1))

	# ----------------------------------------------------------------------
	def upZ(self, event=None):
		self.sliceScale.set(min(self.sliceScale["to"]-1,self.frame+1))

	# ----------------------------------------------------------------------
	def levelChanged(self, s):
		self.updateImage()

	# ----------------------------------------------------------------------
	# Clear current selection
	# ----------------------------------------------------------------------
	def clearDatasetSelection(self):
		self.sliceList.delete(0,END)
		self.clearSlice()

	# ----------------------------------------------------------------------
	def clearSlice(self):
		self.dataset = None
		self.dicom().clear()
		if self._photo:
			self.drawRectangle()
			self.canvas.setImage()
			self._photo = None
		#self.browser.clear()

	# --------------------------------------------------------------------
	# Find dialog
	# --------------------------------------------------------------------
	def find(self, event=None):
		txt = self.findString.get()
		self.search(txt, not self.findCase.get())

	# --------------------------------------------------------------------
	# Search inside dataset information header
	# --------------------------------------------------------------------
	def search(self, needle, upper=False):
		tag = self.browser.search(needle, upper)
		if tag:
			self.browser.selectTag(tag)
		else:
			self.flair.notify("Not found",
				"'%s' not found"%(needle),
				tkFlair.NOTIFY_WARNING)

	# ----------------------------------------------------------------------
	def filterBrowser(self, event=None):
		txt = self.findString.get()
		self.browser.filter(txt, not self.findCase.get())

	# ----------------------------------------------------------------------
	def filterClear(self, event=None):
		self.findString.set("")
		self.browser.filter()

	# --------------------------------------------------------------------
	# Open file in viewer
	# --------------------------------------------------------------------
	def viewer(self, event=None):
		lst = []
		if self.dicom().modality == "RTSTRUCT":
			lst.append(("RTSTRUCT", str(self.dataset)))

		else:
			for i in self.sliceList.curselection():
				name = self.sliceList.get(i)[0]
				try:
					dataset = dicom.read_file(
							os.path.join(self.page.dicominfo["directory"], name),
							force=True,
							defer_size=256,
							stop_before_pixels=True)
					lst.append((name, str(dataset)))
				except dicom.filereader.InvalidDicomError:
					pass

		if lst:
			self.flair.viewer(lst)
