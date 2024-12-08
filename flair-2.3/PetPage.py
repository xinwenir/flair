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
# Date:	03-Oct-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import re
import sys
import math
import string
from math import *
from log import say

try:
	from tkinter import *
except ImportError:
	from tkinter import *
try:
	import configparser
except ImportError:
	import configparser as ConfigParser

import bmath
import undo
import Quad
import tkFlair
import tkExtra

import Input
import Ribbon
import Manual
import Project
import FlairRibbon

CUSTOM    = "custom"
EXTENSION = ".pet"

#===============================================================================
# Pet Frame class
#===============================================================================
class PetPage(FlairRibbon.FlairPage):
	"""Generate the geometry of PET structures"""

	_name_ = "Pet"
	_icon_ = "pet"

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Executable ===========
		group = Ribbon.LabelGroup(self.ribbon, "Templates")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		group.grid2rows()

		# ---
		self.templates = tkExtra.Combobox(group.frame,
					width=40,
					command=self.loadTemplate,
					font=Ribbon._FONT,
					background=Ribbon._BACKGROUND)
		tkExtra.Balloon.set(self.templates, "PET scanner templates")
		self.templates.pack(side=LEFT)

		self.scanTemplates()

		# ------
		group = Ribbon.LabelGroup(self.ribbon, "Action")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		# ---
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["gear32"],
				text="Create",
				compound=TOP,
				command=self.execute,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y)
		tkExtra.Balloon.set(b, "Create PET geometry")


	#----------------------------------------------------------------------
	# Create Project page
	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self)

		self.modVar = IntVar()
		self.nrvar  = IntVar()
		self.nzvar  = IntVar()
		self.navar  = IntVar()

		self.mrvar  = IntVar()
		self.mzvar  = IntVar()
		self.mavar  = IntVar()

		# ==================== Block ====================
		lframe = tkExtra.ExLabelFrame(self.frame, text="Block")
		lframe.grid(row=0, column=0, sticky=NSEW)
		frame = lframe.frame

		self.frame.grid_columnconfigure(0, weight=1)
		self.frame.grid_columnconfigure(1, weight=1)

		# Image
		Label(frame, image=tkFlair.icons["pet_block"]).grid(row=0, column=0, rowspan=11)

		# --- Radial ---
		row = 0
		WIDTH=12
		Label(frame, text="NR").grid(row=row, column=1)
		s = Spinbox(frame, text=self.nrvar,
					width=WIDTH,
					background="White",
					from_=1, to=10000)
		s.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(s, "Number of radial crystals")

		# ----
		row += 1
		Label(frame, text="\u0394R").grid(row=row, column=1)
		self.drentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.drentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.drentry, "Radial dimension of the crystals")

		# ----
		row += 1
		Label(frame, text="\u03B4R").grid(row=row, column=1)
		self.srentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.srentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.srentry, "Radial separation between crystals")

		# ----
		row += 1
		Frame(frame, relief=RAISED, height=2, background="White").grid(
			row=row, column=1, columnspan=3, sticky=EW)

		# --- Axial ---
		row += 1
		Label(frame, text="Nz").grid(row=row, column=1)
		s = Spinbox(frame, text=self.nzvar,
					width=WIDTH,
					background="White",
					from_=1, to=10000)
		s.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(s, "Number of axial crystals")

		# ----
		row += 1
		Label(frame, text="\u0394z").grid(row=row, column=1)
		self.dzentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.dzentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.dzentry, "Axial dimension of the crystals")

		# ----
		row += 1
		Label(frame, text="\u03B4z").grid(row=row, column=1)
		self.szentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.szentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.szentry, "Axial separation between crystals")

		# ----
		row += 1
		Frame(frame, relief=RAISED, height=2, background="White").grid(
			row=row, column=1, columnspan=3, sticky=EW)

		# --- Transaxial ---
		row += 1
		Label(frame, text="N\u03B8").grid(row=row, column=1)
		s = Spinbox(frame, text=self.navar,
					width=WIDTH,
					background="White",
					from_=1, to=10000)
		s.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(s, "Number of transaxial (azimuthal) crystals")

		# ----
		row += 1
		Label(frame, text="\u0394\u03B8").grid(row=row, column=1)
		self.daentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.daentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.daentry, "Transaxial (azimuthal) dimension of the crystals")

		# ----
		row += 1
		Label(frame, text="\u03B4\u03B8").grid(row=row, column=1)
		self.saentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.saentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.saentry, "Transaxial (azimuthal) separation between crystals")

		# ==================== Module ====================
		lframe = tkExtra.ExLabelFrame(self.frame, text="Module")
		lframe.grid(row=1, column=0, sticky=NSEW)
		frame = lframe.frame

		# Image
		Label(frame, image=tkFlair.icons["pet_module"]).grid(row=0, column=0, rowspan=9)

		# --- Radial ---
		row = 0
		Label(frame, text="MR").grid(row=row, column=1)
		s = Spinbox(frame, text=self.mrvar,
					width=WIDTH,
					background="White",
					from_=1, to=10000)
		s.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(s, "Number of radial crystal blocks")

		# ----
		row += 1
		Label(frame, text="SR").grid(row=row, column=1)
		self.lrentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.lrentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.lrentry, "Radial separation between crystal blocks")

		# ----
		row += 1
		Frame(frame, relief=RAISED, height=2, background="White").grid(
			row=row, column=1, columnspan=3, sticky=EW)

		# --- Axial ---
		row += 1
		Label(frame, text="Mz").grid(row=row, column=1)
		s = Spinbox(frame, text=self.mzvar,
					width=WIDTH,
					background="White",
					from_=1, to=10000)
		s.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(s, "Number of axial crystal blocks")

		# ----
		row += 1
		Label(frame, text="Sz").grid(row=row, column=1)
		self.lzentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.lzentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.lzentry, "Axial separation between crystal blocks")

		# ----
		row += 1
		Frame(frame, relief=RAISED, height=2, background="White").grid(
			row=row, column=1, columnspan=3, sticky=EW)

		# --- Transaxial ---
		row += 1
		Label(frame, text="M\u03B8").grid(row=row, column=1)
		s = Spinbox(frame, text=self.mavar,
					width=WIDTH,
					background="White",
					from_=1, to=10000)
		s.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(s, "Number of transaxial (azimuthal) crystal blocks")

		# ----
		row += 1
		Label(frame, text="S\u03B8").grid(row=row, column=1)
		self.laentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.laentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.laentry, "Transaxial (azimuthal) separation between crystal blocks")

		# ==================== Pet ====================
		lframe = tkExtra.ExLabelFrame(self.frame, text="Pet")
		lframe.grid(row=0, column=1, sticky=NSEW)
		frame = lframe.frame

		# Image
		Label(frame, image=tkFlair.icons["pet_array"]).grid(row=0, column=0, rowspan=11)

		# --- Radius ---
		row = 0
		Label(frame, text="Modules").grid(row=row, column=1, sticky=E)
		s = Spinbox(frame, text=self.modVar,
					width=WIDTH,
					background="White",
					from_=0, to=10000)
		s.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(s, "Number of Modules")

		# ----
		row += 1
		Label(frame, text="Radius").grid(row=row, column=1, sticky=E)
		self.rentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.rentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.rentry, "PET radius")

		# ----
		row += 1
		Label(frame, text="Angle").grid(row=row, column=1, sticky=E)
		self.angentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.angentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="deg").grid(row=row, column=3)
		tkExtra.Balloon.set(self.angentry, "Opening angle of PET (0 .. 180deg]")

		# ---- Axis direction cosine x
		row += 1
		Label(frame, text="Axis Ux").grid(row=row, column=1, sticky=E)
		self.auxentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.auxentry.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(self.auxentry, "PET axis x direction cosines as a vector")

		# ---- Axis direction cosine y
		row += 1
		Label(frame, text="Axis Uy").grid(row=row, column=1, sticky=E)
		self.auyentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.auyentry.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(self.auyentry, "PET axis y direction cosines as a vector")

		# ---- Axis direction cosine z
		row += 1
		Label(frame, text="Axis Uz").grid(row=row, column=1, sticky=E)
		self.auzentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.auzentry.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(self.auzentry, "PET axis z direction cosines as a vector")

		# ---- PET scanner center position in x axis
		row += 1
		Label(frame, text="Center X").grid(row=row, column=1, sticky=E)
		self.cnxentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.cnxentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.cnxentry, "PET center position in X axis")

		# ---- PET scanner center position in y axis
		row += 1
		Label(frame, text="Center Y").grid(row=row, column=1, sticky=E)
		self.cnyentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.cnyentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.cnyentry, "PET center position in Y axis")

		# ---- PET scanner center position in z axis
		row += 1
		Label(frame, text="Center Z").grid(row=row, column=1, sticky=E)
		self.cnzentry = tkExtra.FloatEntry(frame,
					width=WIDTH,
					background="White")
		self.cnzentry.grid(row=row, column=2, sticky=EW)
		Label(frame, text="cm").grid(row=row, column=3)
		tkExtra.Balloon.set(self.cnzentry, "PET center position in Z axis")

		# --- Material ---
		row += 1
		Label(frame, text="Material").grid(row=row, column=1, sticky=E)
		self.mat = tkExtra.Combobox(frame, width=10)
		self.mat.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(self.mat, "Select crystal material")

		# --- Region ---
		row += 1
		Label(frame, text="Region").grid(row=row, column=1, sticky=E)
		self.reg = tkExtra.Combobox(frame, width=10)
		self.reg.grid(row=row, column=2, sticky=EW)
		tkExtra.Balloon.set(self.mat, "Select surrounding region to correct")

	# ----------------------------------------------------------------------
	# Scan for templates and fill combo box
	# ----------------------------------------------------------------------
	def scanTemplates(self):
		parser = configparser.RawConfigParser()

		paths = (("D:", os.path.join(tkFlair.prgDir, "templates")),
			 ("U:", os.path.join(tkFlair.iniDir, "templates")))

		lst = []
		for label, path in paths:
			try:
				for filename in os.listdir(path):
					(fn, ext) = os.path.splitext(filename)
					if ext == EXTENSION:
						parser.read(os.path.join(path,filename))
						try: name = parser.get("pet","name")
						except: name = "?"
						lst.append("%s %s - %s"%(label,fn,name))
			except OSError:
				pass
		lst.sort()
		lst.insert(0,CUSTOM)

		self.templates.fill(lst)
		self.templates.set(lst[0])

	# ----------------------------------------------------------------------
	def loadTemplate(self, event=None):
		scanner = self.templates.get()

		self.mat.fill(self.project.input.materialList())

		regions = [x.sdum() for x in self.project.input.cardsSorted("REGION") if x.name()!="&"]
		self.reg.fill(regions)
		if regions:
			self.reg.set(self.project.pet.get("region",regions[-1]))

		if scanner != CUSTOM:	# From project
			self.get()
			scanner = scanner.split()

			if scanner[0] == "D:":
				path = os.path.join(tkFlair.prgDir, "templates")
			else:
				path = os.path.join(tkFlair.iniDir, "templates")
			filename = os.path.join(path, scanner[1]+EXTENSION)
			parser = configparser.RawConfigParser()
			parser.read(filename)

			self.modVar.set(parser.get("pet","modules"))

			# Block
			self.nrvar.set(parser.get("block","nr"))
			self.nzvar.set(parser.get("block","nz"))
			self.navar.set(parser.get("block","na"))

			self.drentry.set(str(parser.get("block","dr")))
			self.srentry.set(str(parser.get("block","sr")))

			self.dzentry.set(str(parser.get("block","dz")))
			self.szentry.set(str(parser.get("block","sz")))

			self.daentry.set(str(parser.get("block","da")))
			self.saentry.set(str(parser.get("block","sa")))

			# Module
			self.mrvar.set(parser.get("module","mr"))
			self.mzvar.set(parser.get("module","mz"))
			self.mavar.set(parser.get("module","ma"))

			self.lrentry.set(str(parser.get("module","lr")))
			self.lzentry.set(str(parser.get("module","lz")))
			self.laentry.set(str(parser.get("module","la")))

			# Pet
			self.rentry.set(str(parser.get("pet","radius")))
			self.angentry.set(str(parser.get("pet","angle")))
			self.auxentry.set(str(parser.get("pet","axisux")))
			self.auyentry.set(str(parser.get("pet","axisuy")))
			self.auzentry.set(str(parser.get("pet","axisuz")))
			self.cnxentry.set(str(parser.get("pet","centrx")))
			self.cnyentry.set(str(parser.get("pet","centry")))
			self.cnzentry.set(str(parser.get("pet","centrz")))

			self.mat.set(parser.get("pet","material"))

	# ----------------------------------------------------------------------
	def get(self):
		if self.project is None: return
		oldhash = self.project.pet.hash()

		var = self.project.pet.var

		# Block
		try: var["nr"] = self.nrvar.get()
		except: pass
		try: var["dr"] = self.drentry.get()
		except: pass
		try: var["sr"] = self.srentry.get()
		except: pass

		try: var["nz"] = self.nzvar.get()
		except: pass
		try: var["dz"] = self.dzentry.get()
		except: pass
		try: var["sz"] = self.szentry.get()
		except: pass

		try: var["na"] = self.navar.get()
		except: pass
		try: var["da"] = self.daentry.get()
		except: pass
		try: var["sa"] = self.saentry.get()
		except: pass

		# Module
		try: var["mr"] = self.mrvar.get()
		except: pass
		try: var["lr"] = self.lrvar.get()
		except: pass

		try: var["mz"] = self.mzvar.get()
		except: pass
		try: var["lz"] = self.lzvar.get()
		except: pass

		try: var["ma"] = self.mavar.get()
		except: pass
		try: var["la"] = self.lavar.get()
		except: pass

		# Pet
		try: var["modules"]  = self.modVar.get()
		except: pass
		try: var["radius"]   = self.rentry.get()
		except: pass
		try: var["angle"]    = self.angentry.get()
		except: pass
		try: var["axisux"]   = self.auxentry.get()
		except: pass
		try: var["axisuy"]   = self.auyentry.get()
		except: pass
		try: var["axisuz"]   = self.auzentry.get()
		except: pass
		try: var["centrx"]   = self.cnxentry.get()
		except: pass
		try: var["centry"]   = self.cnyentry.get()
		except: pass
		try: var["centrz"]   = self.cnzentry.get()
		except: pass
		try: var["material"] = self.mat.get()
		except: pass
		try: var["region"]   = self.reg.get()
		except: pass

		if oldhash != self.project.pet.hash():
			self.setModified(True)

	# ----------------------------------------------------------------------
	def execute(self):
		err = ""

		# Block
		try: Nr = int(self.nrvar.get())
		except: err += "Invalid NR\n"
		try: Dr = float(self.drentry.get())
		except: err += "Invalid DR\n"
		try: Sr = float(self.srentry.get())
		except: err += "Invalid dR\n"

		try: Nz = int(self.nzvar.get())
		except: err += "Invalid Nz\n"
		try: Dz = float(self.dzentry.get())
		except: err += "Invalid Dz\n"
		try: Sz = float(self.szentry.get())
		except: err += "Invalid dz\n"

		try: Na = int(self.navar.get())
		except: err += "Invalid Ntheta\n"
		try: Da = float(self.daentry.get())
		except: err += "Invalid Dtheta\n"
		try: Sa = float(self.saentry.get())
		except: err += "Invalid dtheta\n"

		# Module
		try: Mr = int(self.mrvar.get())
		except: err += "Invalid MR\n"
		try: Lr = float(self.lrentry.get())
		except: err += "Invalid SR\n"

		try: Mz = int(self.mzvar.get())
		except: err += "Invalid Mz\n"
		try: Lz = float(self.lzentry.get())
		except: err += "Invalid Sz\n"

		try: Ma = int(self.mavar.get())
		except: err += "Invalid Mtheta\n"
		try: La = float(self.laentry.get())
		except: err += "Invalid Stheta\n"

		# Pet
		try: modules  = int(self.modVar.get())
		except: err += "Invalid Modules\n"
		try: radius   = float(self.rentry.get())
		except: err += "Invalid Radius\n"
		try: angle    = float(self.angentry.get())
		except: err += "Invalid Angle\n"
		try: axisux  = float(self.auxentry.get())
		except: err += "Invalid Axis X Component\n"
		try: axisuy  = float(self.auyentry.get())
		except: err += "Invalid Axis Y Component\n"
		try: axisuz  = float(self.auzentry.get())
		except: err += "Invalid Axis Z Component\n"
		try: centrx  = float(self.cnxentry.get())
		except: err += "Invalid Scanner X Center\n"
		try: centry  = float(self.cnyentry.get())
		except: err += "Invalid Scanner Y Center\n"
		try: centrz  = float(self.cnzentry.get())
		except: err += "Invalid Scanner Z Center\n"
		try: material = self.mat.get()
		except: err += "Invalid Material\n"
		try: region   = self.reg.get()
		except: err += "Invalid Region\n"

		if err:
			self.flair.notify("Error",
				"Invalid variables %s"%(err),
				tkFlair.NOTIFY_ERROR)
			return

		def fmt(x): return bmath.format(x)

		# Check consistency
		if modules % 2 == 1:
			self.flair.notify("Odd number of modules",
				"Only even number of modules are supported",
				tkFlair.NOTIFY_ERROR)
			return False

		# Find location and dimensions
		DDx = Ma*(Na*Da + (Na-1)*Sa) + (Ma-1)*La
		DDy = Mr*(Nr*Dr + (Nr-1)*Sr) + (Mr-1)*Lr
		DDz = Mz*(Nz*Dz + (Nz-1)*Sz) + (Mz-1)*Lz

		MNa = Ma * Na
		MNr = Mr * Nr
		MNz = Mz * Nz

		DDx2 =  DDx / 2.0
		x    = -DDx2
		y    =  radius
		z    = -DDz / 2.0

		# Normalize axis vector and translate it into polar and azimuthal angles
		Axnr = math.sqrt( axisux**2. + axisuy**2. + axisuz**2. )
		if Axnr == 0.0: self.flair.notify("Wrong axis vector",
				"Please introduce a valid axis vector",
				tkFlair.NOTIFY_ERROR)
		Axpl = math.acos( axisuz / Axnr )*180./math.pi
		try:
			Axat = math.atan( axisuy / axisux )*180./math.pi
		except ZeroDivisionError:
			Axat = 0.0
			if axisuy != 0.0: Axat = 90.0

		# If modules is not supplied calculate the roundest number of modules
		phi = math.degrees(2.0*math.atan2(DDx2, radius))
		if modules == 0:
			if angle < 180:
				n = int(2.0*angle/phi)
			else:
				n = int(360.0/phi)
			if n%2 == 1: n -= 1
			modules = n

		initialAngle = 0.0
		if angle >= 180.0:
			angstep = 360.0 / float(modules)
			hang    = 0.0
		else:
			angstep = angle / (modules/2 - 1)
			hang    = angle/2+0.001

			# Check if we have odd or even modules in half group
			if (modules/2)%2 == 0:
				# Rotate everything by -angstep/2
				initialAngle = angstep/2.0

		if Input._developer:
			say("modules=",modules)
			say("angle=",angle)
			say("DDx=",DDx)
			say("DDy=",DDy)
			say("DDz=",DDz)
			say("angstep=",angstep)
			say("phi=",phi)
			say("initialAngle=",initialAngle)
			say("hang=",hang)

		if angstep < phi:
			self.flair.notify("Too many modules",
				"Too many modules are selected resulting in overlapping regions",
				tkFlair.NOTIFY_ERROR)
			return False

		transfo = []

		# Rotate axis basic block
		if Axpl != 0.0:
			transfo.append(["pet000",200,"",fmt(-Axpl)])
		if Axat != 0.0:
			transfo.append(["pet000","","",fmt(-Axat)])
		if centrx != 0.0 or centry != 0.0 or centrz != 0.0:
			transfo.append(["pet000","","","",centrx,centry,centrz])

		# Rotate basic block
		if initialAngle != 0.0:
			transfo.append(["pet001","","",fmt(initialAngle)])
		if Axpl != 0.0:
			transfo.append(["pet001",200,"",fmt(-Axpl)])
		if Axat != 0.0:
			transfo.append(["pet001","","",fmt(-Axat)])
		if centrx != 0.0 or centry != 0.0 or centrz != 0.0:
			transfo.append(["pet001","","","",centrx,centry,centrz])

		# --- CHECK INPUT ---
		input = self.project.input
		if not input.cardlist:
			input.minimumInput()

		# --- BODIES ---
		# Card generation variables
		if "END" not in input.cards:
			# create the bare minimum input
			input.minimumInput()
		pos = input["END"][0].pos()

		# Undo
		undoinfo = [(self.flair.refreshUndo,), self.flair.renumberUndo()]

		# name mapping
		_name = {	"PETx0001"         : "PETl001",  # Left
				"PETx%04d"%(2*MNa) : "PETr001",  # Right
				"PETy0001"         : "PETb001",  # Bottom
				"PETy%04d"%(2*MNr) : "PETt001",  # Top
				"PETz0001"         : "PETback",  # Back
				"PETz%04d"%(2*MNz) : "PETfront", # Front
			}
		def name(x): return _name.get(x,x)
		def rname(x,s="+"):
			n = _name.get(x,x)
			if n[0] != "-":
				return s+n
			elif s=="-":
				return "+"+n[1:]
			else:
				return n

		# First two concentric ZCC
		if Axat != 0.0 or Axpl != 0.0 or centrx != 0.0 or centry != 0.0 or centrz != 0.0:
			# Add transformation
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("$start_transform",["","pet000"]),pos,False))
			pos += 1

		undoinfo.append(
			self.flair.addCardUndo(Input.Card("ZCC",
				["PETco", "0.0", "0.0", fmt(sqrt((radius+DDy)**2+DDx2**2)+0.0001)]),pos,False))
		pos += 1
		undoinfo.append(
			self.flair.addCardUndo(
				Input.Card("ZCC",["PETci", "0.0", "0.0", fmt(radius-0.0001)]),pos,False))
		pos += 1

		# Rotate if needed basic axis block
		if Axat != 0.0 or Axpl != 0.0 or centrx != 0.0 or centry != 0.0 or centrz != 0.0:
			# Add transformation
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("$end_transform"),pos,False))
			pos += 1

		# Rotate if needed basic block
		if initialAngle != 0.0 or Axpl != 0.0 or centrx != 0.0 or centry != 0.0 or centrz != 0.0:
			# Add transformation
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("$start_transform",["","pet001"]),pos,False))
			pos += 1

		# Along X
		idx = 1
		for i in range(MNa):
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("YZP",[name("PETx%04d"%(idx)), fmt(x)]),pos,False))
			idx += 1; pos += 1

			x += Da

			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("YZP",[name("PETx%04d"%(idx)), fmt(x)]),pos,False))
			idx += 1; pos += 1

			if (i+1)%Na == 0:
				x += La
			else:
				x += Sa

		# Along Y
		idx = 1
		for i in range(MNr):
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("XZP",[name("PETy%04d"%(idx)), fmt(y)]),pos,False))
			idx += 1; pos += 1

			y += Dr

			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("XZP",[name("PETy%04d"%(idx)), fmt(y)]),pos,False))
			idx += 1; pos += 1

			if (i+1)%Nr == 0:
				y += Lr
			else:
				y += Sr

		# Along z
		idx = 1
		for i in range(MNz):
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("XYP",[name("PETz%04d"%(idx)), fmt(z)]),pos,False))
			idx += 1; pos += 1

			z += Dz

			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("XYP",[name("PETz%04d"%(idx)), fmt(z)]),pos,False))
			idx += 1; pos += 1

			if (i+1)%Nz == 0:
				z += Lz
			else:
				z += Sz

		# Rotate if needed basic block
		if initialAngle != 0.0 or Axpl != 0.0 or centrx != 0.0 or centry != 0.0 or centrz != 0.0:
			# Add transformation
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("$end_transform"),pos,False))
			pos += 1

		# Create additional planes for replicas and transformations
		idx = 2

		# Break points
		if hang>0:
			break1 = int(round(float(modules)/4.0)) - 1
			break2 = break1 + modules/2
		else:
			break1 = break2 = -1

		ang = initialAngle + angstep
		for i in range(1,modules):
			if initialAngle != 0.0 or Axat != 0.0 or Axpl != 0.0 or centrx != 0.0 or centry != 0.0 or centrz != 0.0:
				trname = "petb%03d"%(idx)
				if centrx != 0.0 or centry != 0.0 or centrz != 0.0:
					transfo.append([trname,"","","",-centrx,-centry,-centrz])
				if Axat != 0.0:
					transfo.append([trname,"","",fmt(Axat)])
				if Axpl != 0.0:
					transfo.append([trname,200,"",fmt(Axpl)])
				transfo.append([trname,"","",fmt(-ang)])

				if centrx != 0.0 or centry != 0.0 or centrz != 0.0:
					transfo.append(["pet%03d"%(idx),"","","",-centrx,-centry,-centrz])
				if Axat != 0.0:
					transfo.append(["pet%03d"%(idx),"","",fmt(Axat)])
				if Axpl != 0.0:
					transfo.append(["pet%03d"%(idx),200,"",fmt(Axpl)])
				transfo.append(["pet%03d"%(idx),"","",fmt(-ang+initialAngle)])
				if Axpl != 0.0:
					transfo.append(["pet%03d"%(idx),200,"",fmt(-Axpl)])
				if Axat != 0.0:
					transfo.append(["pet%03d"%(idx),"","",fmt(-Axat)])
				if centrx != 0.0 or centry != 0.0 or centrz != 0.0:
					transfo.append(["pet%03d"%(idx),"","","",centrx,centry,centrz])
			else:
				trname = "pet%03d"%(idx)
				transfo.append([trname,"","",fmt(-ang)])

			# Add transformation
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("$start_transform",["","-"+trname]),pos,False))
			pos += 1

			if ang < 179.9999:
				undoinfo.append(
					self.flair.addCardUndo(
						Input.Card("YZP",["PETl%03d"%(idx),-DDx2]),pos,False))
				pos += 1

				undoinfo.append(
					self.flair.addCardUndo(
						Input.Card("YZP",["PETr%03d"%(idx),DDx2]),pos,False))
				pos += 1
			else:
				_name["PETl%03d"%(idx)] = "-PETr%03d"%(idx-modules/2)
				_name["PETr%03d"%(idx)] = "-PETl%03d"%(idx-modules/2)

			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("XZP",["PETb%03d"%(idx),radius]),pos,False))
			pos += 1

			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("XZP",["PETt%03d"%(idx),radius+DDy]),pos,False))
			pos += 1

			# End transformation
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("$end_transform"),pos,False))
			pos += 1
			idx += 1

			if i == break1:
				ang = 180.0 - angle / 2.0
			elif i == break2:
				ang = 360.0 - angle / 2.0
			else:
				ang += angstep

		# XXX
		input.renumber()

		# --- REGIONS ---
		# Correct surrounding region
		if region != "":
			for card in input["REGION"]:
				if card.ignore(): continue
				if card.sdum() == region:
					undoinfo.append(
						self.flair.setWhatUndo(card, -1,
							card.extra()+" -(+PETco -PETci +PETfront -PETback)"))

		# Where to add cards
		pos = input["END"][1].pos()

		# Gaps between modules
		expr = ""
		for mod in range(1,modules+1):
			nxt = mod % modules + 1
			right = rname("PETl%03d"%(nxt),"+")
			left  = rname("PETr%03d"%(mod),"-")
			expr += "| %s %s +PETfront -PETback +PETco -PETci\n"%(right,left)

		# Gaps between modules and cylinders
		for mod in range(1,modules+1):
			right = rname("PETr%03d"%(mod),"+")
			left  = rname("PETl%03d"%(mod),"-")

			if mod <= modules/2:
				bottom  = rname("PETb%03d"%(mod),"+")
				bottom2 = rname("PETb%03d"%(mod+modules/2),"+")
				expr += "| %s %s %s %s +PETfront -PETback -PETci\n"%(right,left,bottom,bottom2)

			top   = rname("PETt%03d"%(mod),"-")
			expr += "| %s %s %s +PETfront -PETback +PETco\n"%(right,left,top)

		undoinfo.append(
			self.flair.addCardUndo(
				Input.Card("REGION",["PETgaps"],"",expr[2:-1]),pos,False))
		pos += 1

		# Crystal Gaps

		expr = ""
		lim = "-PETb001 +PETt001 -PETback +PETfront\n"
		for i in range(2,2*MNa,2):
			expr += "| -PETx%04d +PETx%04d %s"%(i, i+1, lim)

		lim = "-PETl001 +PETr001 -PETback +PETfront\n"
		for i in range(2,2*MNr,2):
			expr += "| -PETy%04d +PETy%04d %s"%(i, i+1, lim)

		lim = "-PETl001 +PETr001 -PETb001 +PETt001\n"
		for i in range(2,2*MNz,2):
			expr += "| -PETz%04d +PETz%04d %s"%(i, i+1, lim)

		undoinfo.append(
			self.flair.addCardUndo(
				Input.Card("REGION",["PETmgap"],"",expr[2:-1]),pos,False))
		pos += 1

		# Create basic module
		n = 0
		for k in range(1,2*MNz,2):
			for j in range(1,2*MNr,2):
					for i in range(1,2*MNa,2):
						right  = rname("PETx%04d"%(i+1),"+")
						left   = rname("PETx%04d"%(i), "-")
						top    = rname("PETy%04d"%(j+1),"+")
						bottom = rname("PETy%04d"%(j), "-")
						front  = rname("PETz%04d"%(k+1),"+")
						back   = rname("PETz%04d"%(k), "-")
						expr   = "%s %s %s %s %s %s" % \
							(right,left, top,bottom, front, back)
						undoinfo.append(
							self.flair.addCardUndo(
								Input.Card("REGION",
									["PET%05d"%(n)],
									"",
									expr),
									pos,False))
						n += 1
						pos += 1

		# Create empty lattices
		for mod in range(2,modules+1):
			right  = rname("PETr%03d"%(mod),"+")
			left   = rname("PETl%03d"%(mod),"-")
			top    = rname("PETt%03d"%(mod),"+")
			bottom = rname("PETb%03d"%(mod),"-")
			expr = "%s %s %s %s +PETfront -PETback" % (right,left, top,bottom)
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("REGION",["PETm%03d"%(mod)],"",expr),pos,False))
			pos += 1

		# --- LATTICEs ---
		pos += 1	# Skip END
		for mod in range(2,modules+1):
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("LATTICE",["pet%03d"%(mod),"PETm%03d"%(mod)]),pos,False))
			pos += 1

		# XXX
		input.renumber()

		# --- ASSIGNMAT ---
		if "ASSIGNMA" not in input.cards:
			pos = input["GEOEND"][0].pos()
		else:
			pos = input["ASSIGNMA"][-1].pos()

		undoinfo.append(
			self.flair.addCardUndo(
				Input.Card("ASSIGNMA",["","AIR","PETgaps"]),pos,False))
		pos += 1
		undoinfo.append(
			self.flair.addCardUndo(
				Input.Card("ASSIGNMA",["","AIR","PETmgap"]),pos,False))
		pos += 1

		# XXX
		undoinfo.append(
			self.flair.addCardUndo(
				Input.Card("ASSIGNMA",["",material,"PET00000","PET%05d"%(n-1)]),pos,False))
		undoinfo.append(
			self.flair.addCardUndo(
				Input.Card("ASSIGNMA",["","VACUUM","PETm002","PETm%03d"%(modules)]),pos,False))
		pos += 1

		# XXX
		input.renumber()

		# --- ROT-DEFI ---
		if "ROT-DEFI" in input.cards:
			pos = input["ROT-DEFI"][-1].pos()

		for tr in transfo:
			undoinfo.append(
				self.flair.addCardUndo(
					Input.Card("ROT-DEFI",tr),pos,False))
			pos += 1

		# Finalize
		undoinfo.append(self.flair.renumberUndo())
		undoinfo.append(undoinfo[0])
		self.flair.addUndo(undo.createListUndo(undoinfo,"Add PET geometry"))
		self.flair.setInputModified()
		self.flair.refresh("card")

		self.flair.notify("PET",
			"Pet inserted to input\n" \
			"Don't forget to add the necessary materials if any.")
