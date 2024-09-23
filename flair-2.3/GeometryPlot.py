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
# Date:	05-Sep-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import time
import math
import bmath
import string
import signal

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import Input
import tkFlair
import tkExtra
import Unicode
import Project

import Plot
import BasePlot

#===============================================================================
# Geometry Plot Frame
#===============================================================================
class GeometryPlotFrame(BasePlot.BasePlotFrame):
	def __init__(self, master, flair):
		BasePlot.BasePlotFrame.__init__(self, master, flair)

		# Enable aspect ratio calculation
		#self.calcAspect.config(state=NORMAL)

		# === Fill main frame ===

		rowframe = Frame(self)
		rowframe.pack(side=TOP,fill=X)

		# --- Origin/Center ---
		frame = LabelFrame(rowframe, text="Center ")
		frame.pack(side=LEFT, expand=YES, fill=X, padx=3)

		Label(frame, text="x:", anchor=W).grid(row=0,column=0)
		self.x_origin =  tkExtra.VectorEntry(frame, background="White")
		self.x_origin.grid(row=0,column=1,sticky=EW)
		self.x_origin.bind("<Return>",   self.getOrigin)
		self.x_origin.bind("<KP_Enter>", self.getOrigin)

		Label(frame, text="y:", anchor=W).grid(row=1,column=0)
		self.y_origin =  tkExtra.FloatEntry(frame, background="White")
		self.y_origin.grid(row=1,column=1,sticky=EW)

		Label(frame, text="z:", anchor=W).grid(row=2,column=0)
		self.z_origin =  tkExtra.FloatEntry(frame, background="White")
		self.z_origin.grid(row=2,column=1,sticky=EW)

		tkExtra.Balloon.set(self.x_origin, "X coordinate of the center")
		tkExtra.Balloon.set(self.y_origin, "Y coordinate of the center")
		tkExtra.Balloon.set(self.z_origin, "Z coordinate of the center")

		frame.grid_columnconfigure(1, weight=1)

		self.widget["@x"] = self.x_origin
		self.widget["@y"] = self.y_origin
		self.widget["@z"] = self.z_origin

#		b = Button(frame, text="Reset", padx=2, pady=0,
#				command=self.reset_origin)
#		b.pack(side=LEFT)
#		tkExtra.Balloon.set(b, "Reset center to (0,0,0)")

		# --- Basis ---
		frame = LabelFrame(rowframe, text="Basis ")
		frame.pack(side=LEFT, expand=YES, fill=X, padx=3)

		# ---------- Axes label --------
		f = Frame(frame)
		f.pack(side=TOP, fill=X)
		self.basisOptions = BooleanVar()
		b = Checkbutton(f, text="Axes "+Unicode.BLACK_DOWN_POINTING_TRIANGLE,
				variable=self.basisOptions,
				selectcolor="LightYellow",
				indicatoron=FALSE,
				command=self.showBasis)
		tkExtra.Balloon.set(b, "Show basis options")
		b.pack(side=LEFT)

		self.basisAxes = Label(f, width=6, background="White", foreground="DarkBlue")
		self.basisAxes.pack(side=LEFT, expand=YES, fill=X)

		# -----
		self.basisButtons = Frame(frame)
		self.basisButtons.pack(side=BOTTOM, expand=YES, fill=BOTH, anchor=CENTER)

		bxy = Button(self.basisButtons, text="x-y",
				padx=2, pady=0, command=self.basis_xy)
		bxy.grid(row=1, column=0, sticky=EW)
		tkExtra.Balloon.set(bxy, "Use X-Y as basis")

		bxz = Button(self.basisButtons, text="x-z",
				padx=2, pady=0, command=self.basis_xz)
		bxz.grid(row=2, column=0, sticky=EW)
		tkExtra.Balloon.set(bxz, "Use X-Z as basis")

		# -----
		byz = Button(self.basisButtons, text="y-z",
				padx=2, pady=0, command=self.basis_yz)
		byz.grid(row=1, column=1, sticky=EW)
		tkExtra.Balloon.set(byz, "Use Y-Z as basis")

		bswap = Button(self.basisButtons, text="swap",
				padx=2, pady=0, command=self.basis_swap)
		bswap.grid(row=2, column=1, sticky=EW)
		tkExtra.Balloon.set(bswap, "Swap basis from A-B to B-A")

		# -----
		bmu = Button(self.basisButtons, text="-u",
				width=3, padx=2, pady=0, command=self.basis_negu)
		bmu.grid(row=1, column=2, sticky=EW)
		tkExtra.Balloon.set(bmu, "Negate vector U")
		bmv = Button(self.basisButtons, text="-v",
				width=3, padx=2, pady=0, command=self.basis_negv)
		bmv.grid(row=2, column=2, sticky=EW)
		tkExtra.Balloon.set(bmv, "Negate vector V")

		# ----- Custom basis -----
		self.basisFrame = Frame(frame)
#		self.basisFrame.pack(side=BOTTOM, expand=YES, fill=X)

		Label(self.basisFrame, text="u:", anchor=W).grid(row=0, column=0, sticky=W)
		Label(self.basisFrame, text="v:", anchor=W).grid(row=1, column=0, sticky=W)

		self.ux_basis =  tkExtra.VectorEntry(self.basisFrame,
					background="White", width=6)
		self.ux_basis.grid(row=0, column=1, sticky=EW)
		self.uy_basis =  tkExtra.FloatEntry(self.basisFrame,
					background="White", width=6)
		self.uy_basis.grid(row=0, column=2, sticky=EW)
		self.uz_basis =  tkExtra.FloatEntry(self.basisFrame,
					background="White", width=6)
		self.uz_basis.grid(row=0, column=3, sticky=EW)
		self.ux_basis.bind("<Return>",   self.getBasis)
		self.ux_basis.bind("<KP_Enter>", self.getBasis)

		self.vx_basis =  tkExtra.VectorEntry(self.basisFrame,
					background="White", width=6)
		self.vx_basis.grid(row=1, column=1, sticky=EW)
		self.vy_basis =  tkExtra.FloatEntry(self.basisFrame,
					background="White", width=6)
		self.vy_basis.grid(row=1, column=2, sticky=EW)
		self.vz_basis =  tkExtra.FloatEntry(self.basisFrame,
					background="White", width=6)
		self.vz_basis.grid(row=1, column=3, sticky=EW)
		self.vx_basis.bind("<Return>",   self.getBasis)
		self.vx_basis.bind("<KP_Enter>", self.getBasis)

		tkExtra.Balloon.set(self.ux_basis, "Define basis vector U")
		tkExtra.Balloon.set(self.uy_basis, "Define basis vector U")
		tkExtra.Balloon.set(self.uz_basis, "Define basis vector U")
		tkExtra.Balloon.set(self.vx_basis, "Define basis vector V")
		tkExtra.Balloon.set(self.vy_basis, "Define basis vector V")
		tkExtra.Balloon.set(self.vz_basis, "Define basis vector V")

		self.basisFrame.grid_columnconfigure(1, weight=1)
		self.basisFrame.grid_columnconfigure(2, weight=1)
		self.basisFrame.grid_columnconfigure(3, weight=1)

		self.widget["@axes"] = (self.basisAxes,
			self.ux_basis, self.uy_basis, self.uz_basis,
			self.vx_basis, self.vy_basis, self.vz_basis,
			bxy, bxz, byz, bswap, bmu, bmv)

		# --- Style ---
#		frame = LabelFrame(mframe, text="Style ")
#		frame.pack(side=LEFT, fill=X)

#		b = Button(frame, text="Palette", padx=2, pady=0,
#				compound=LEFT, image=tkFlair.icons["palette"],
#				command=self.editPalette)
#		b.grid(row=0, column=0, columnspan=2, sticky=EW)
#		tkExtra.Balloon.set(b, "Edit Material Palette")
#		Label(frame, text="Font:", anchor=W).grid(row=1, column=0, sticky=W)
#		self.label_font = Entry(frame, background="White", width=10)
#		self.label_font.grid(row=1, column=1, sticky=EW)
#		tkExtra.Balloon.set(self.label_font, "Label font")

		# --- Extends ---
		frame = LabelFrame(rowframe, text="Extends ")
		frame.pack(side=LEFT, expand=YES, fill=BOTH, padx=3)

		Label(frame, text="\u0394u:", anchor=W).grid(
				row=0,column=0, sticky=W)
		Label(frame, text="\u0394v:", anchor=W).grid(
				row=1, column=0, sticky=W)

		self.u_extend = tkExtra.FloatEntry(frame, background="White", width=7)
		self.u_extend.grid(row=0, column=1, sticky=EW)
		self.v_extend = tkExtra.FloatEntry(frame, background="White", width=7)
		self.v_extend.grid(row=1, column=1, sticky=EW)
		tkExtra.Balloon.set(self.u_extend, "Extend in U axis. Full width is twice")
		tkExtra.Balloon.set(self.v_extend, "Extend in V axis. Full height is twice")
		frame.grid_columnconfigure(1, weight=1)

		b = Button(frame, text="Get", padx=2, pady=0, command=self.get_extend)
		b.grid(row=2, column=1, sticky=E)
		tkExtra.Balloon.set(b, "Get extends and origin from plot window")

		self.widget["@du"] = self.u_extend
		self.widget["@dv"] = self.v_extend

		# --- Plot Options ---
		frame = LabelFrame(rowframe, text="Plot ")
		frame.pack(side=LEFT, fill=BOTH, padx=3)

		# --- type ---
		Label(frame,text="Type:").grid(row=0, column=0, sticky=E)
		b = tkExtra.Combobox(frame, width=10)
		b.fill(Plot.SUBTYPES)
		b.set(Plot.SUBTYPES[0])
		b.grid(row=0, column=1, sticky=EW)
		tkExtra.Balloon.set(b, "Select information to plot")
		self.widget["subtype"] = b

		# ---run ---
		Label(frame,text="Run:").grid(row=1, column=0, sticky=E)
		b = tkExtra.Combobox(frame, width=12)
		b.grid(row=1, column=1, sticky=EW)
		tkExtra.Balloon.set(b, "Select run configuration to use")
		self.widget["run"] = b

		# ---options---
		self.advancedOptions = BooleanVar()
		b = Checkbutton(frame, text="Advanced "+Unicode.BLACK_DOWN_POINTING_TRIANGLE,
				variable=self.advancedOptions,
				indicatoron=FALSE,
				selectcolor="LightYellow",
				command=self.showAdvancedOptions)
		tkExtra.Balloon.set(b, "Show additional plot options")
		b.grid(row=2, column=1, sticky=E)

		frame.grid_rowconfigure(1, weight=1)

		# ===========================================================
		self.optionsFrame = Frame(self)
		#self.optionsFrame.pack(side=TOP,fill=X)

		# --- Options ---
		frame = LabelFrame(self.optionsFrame, text="Options ", foreground="DarkBlue")
		frame.pack(side=RIGHT, fill=X)

		# ---
		Label(frame, text="Vector Scale:", anchor=W).grid(row=0, column=0, sticky=W)
		b = tkExtra.FloatEntry(frame, background="White", width=5)
		b.grid(row=0, column=1, sticky=EW)
		tkExtra.Balloon.set(b, "Scaling factor for magnetic field plot vectors")
		self.widget["vectscale"] = b

		# ---
		Label(frame, text="Plot Coordinates:").grid(row=1, column=0, sticky=E)
		b = tkExtra.Combobox(frame, width=5)
		b.fill(Plot.COORD)
		b.set(Plot.COORD[2])
		b.grid(row=1, column=1, sticky=EW)
		tkExtra.Balloon.set(b, "Select plot coordinates to use")
		frame.grid_columnconfigure(1, weight=1)
		self.widget["coord"] = b

		# ---
		var = IntVar()
		var.set(0)
		b = Checkbutton(frame, text="boundaries", variable=var, anchor=W)
		b.grid(row=0, column=2, sticky=W)
		tkExtra.Balloon.set(b, "Plot boundaries between materials or regions")
		self.widget["boundaries"] = b
		self.variable["boundaries"] = var

		# ---
		var = IntVar()
		var.set(0)
		b = Checkbutton(frame, text="labels", variable=var, anchor=W)
		b.grid(row=1, column=2, sticky=W)
		tkExtra.Balloon.set(b, "Plot labels of materials, regions or lattices")
		self.widget["labels"] = b
		self.variable["labels"] = var

		# --- Grid ---
		frame = LabelFrame(self.optionsFrame, text="Grid ", foreground="DarkBlue")
		frame.pack(side=RIGHT, fill=X)

		Label(frame, text="Nu:", anchor=W).grid(row=0, column=0, sticky=W)
		b = tkExtra.IntegerEntry(frame, background="White", width=5)
		b.grid(row=0, column=1, sticky=EW)
		tkExtra.Balloon.set(b, "Number of bins horizontal in U axis")
		self.widget["nu"] = b

		Label(frame, text="Nv:", anchor=W).grid(row=1, column=0, sticky=W)
		b = tkExtra.IntegerEntry(frame, background="White", width=5)
		b.grid(row=1, column=1, sticky=EW)
		tkExtra.Balloon.set(b, "Number of bins horizontal in V axis")
		self.widget["nv"] = b

		self.bindToggle()

		# --- Reset all values ---
		self.reset_origin()
		self.reset_basis()
		self.reset_extend()

	# ----------------------------------------------------------------------
	def type(self): return "Geometry"

	# ----------------------------------------------------------------------
	def showBasis(self):
		if self.basisOptions.get():
			self.basisButtons.pack_forget()
			self.basisFrame.pack(side=BOTTOM, expand=YES, fill=X)
		else:
			self.basisButtons.pack(side=BOTTOM, expand=YES, fill=X)
			self.basisFrame.pack_forget()

	# ----------------------------------------------------------------------
	def showAdvancedOptions(self):
		if self.advancedOptions.get():
			self.optionsFrame.pack(side=TOP, fill=X)
		else:
			self.optionsFrame.pack_forget()

	# ----------------------------------------------------------------------
	# Show Frame
	# ----------------------------------------------------------------------
	def set(self, plot):
		BasePlot.BasePlotFrame.set(self, plot)
		if plot is None: return

		# Special handling
		try:
			origin = bmath.Vector(list(map(float,plot.var["origin"].split())))
			self.setOrigin(origin)
		except: pass
		try:
			u = bmath.Vector(list(map(float,plot.var["basisu"].split())))
			v = bmath.Vector(list(map(float,plot.var["basisv"].split())))
			self.setBasis(u,v)
		except: pass
		try:
			extU,extV = list(map(float, plot.var["extends"].split()))
			BasePlot.setFloat(self.u_extend, extU)
			BasePlot.setFloat(self.v_extend, extV)
		except: pass

		# Fill runs
		combo = self.widget["run"]
		combo.clear()
		for run in self.project.runs:
			if run.parent != "": continue
			runname = run.getInputName2()
			combo.insert(END, runname)

	# ----------------------------------------------------------------------
	def get(self, plot):
		BasePlot.BasePlotFrame.get(self, plot)
		if plot is None: return

		# Special treatment for ENABLE/DISABLE
		# Origin
		if self.x_origin["state"]==NORMAL or \
		   self.y_origin["state"]==NORMAL or \
		   self.z_origin["state"]==NORMAL:
			O = self.getOrigin()
			x,y,z = list(map(float,plot.get("origin","0 0 0").split()))
			if self.x_origin["state"]==NORMAL: x = O.x()
			if self.y_origin["state"]==NORMAL: y = O.y()
			if self.z_origin["state"]==NORMAL: z = O.z()
			plot.var["origin"] = "%g %g %g" % (x, y, z)

		# Axes
		axes = self.basisAxes["state"]==NORMAL
		U,V = None,None
		if axes or self.ux_basis["state"]==NORMAL \
			or self.uy_basis["state"]==NORMAL \
			or self.uz_basis["state"]==NORMAL:
				U,V  = self.getBasis()
				x,y,z = list(map(float,plot.get("basisu","1 0 0").split()))
				if axes or self.ux_basis["state"]==NORMAL: x = U.x()
				if axes or self.uy_basis["state"]==NORMAL: y = U.y()
				if axes or self.uz_basis["state"]==NORMAL: z = U.z()
				plot.var["basisu"] = "%g %g %g" % (x, y, z)

		if axes or self.vx_basis["state"]==NORMAL \
			or self.vy_basis["state"]==NORMAL \
			or self.vz_basis["state"]==NORMAL:
				if U is None: U,V  = self.getBasis()
				x,y,z = list(map(float,plot.get("basisv","0 1 0").split()))
				if axes or self.vx_basis["state"]==NORMAL: x = V.x()
				if axes or self.vy_basis["state"]==NORMAL: y = V.y()
				if axes or self.vz_basis["state"]==NORMAL: z = V.z()
				plot.var["basisv"] = "%g %g %g" % (x, y, z)

		# Extends
		if self.u_extend["state"]==NORMAL or \
		   self.v_extend["state"]==NORMAL:
			extU = self.u_extend.getfloat()
			extV = self.v_extend.getfloat()
			extModified = False
			if extU<1e-16:
				extModified = True
				extU = 100.0
				BasePlot.setFloat(self.u_extend, extU)
				extModified = True
				extV = 100.0
				BasePlot.setFloat(self.v_extend, extV)
			if extModified:
				# nonsense extends
				self.flair.notify("Nonsense extends",
					"Nonsense value on plot extends:\n" \
					"Value is forced to 100.0",
					tkFlair.NOTIFY_ERROR)
			x,y = list(map(float,plot.get("extends","100 100").split()))
			if self.u_extend["state"]==NORMAL: x = extU
			if self.v_extend["state"]==NORMAL: y = extV
			plot.var["extends"] = "%g %g" % (x, y)

	# ----------------------------------------------------------------------
	def getOrigin(self, event=None):
		self.x_origin.split((self.y_origin, self.z_origin))
		return bmath.Vector(
				self.x_origin.getfloat(),
				self.y_origin.getfloat(),
				self.z_origin.getfloat())

	# ----------------------------------------------------------------------
	def setOrigin(self, x, y=None, z=None):
		if y is None:
			z = x.z()
			y = x.y()
			x = x.x()
		BasePlot.setFloat(self.x_origin, x)
		BasePlot.setFloat(self.y_origin, y)
		BasePlot.setFloat(self.z_origin, z)

	# ----------------------------------------------------------------------
	def getBasis(self, event=None):
		self.ux_basis.split((self.uy_basis, self.uz_basis))
		self.vx_basis.split((self.vy_basis, self.vz_basis))

		U  = bmath.Vector(
				self.ux_basis.getfloat(),
				self.uy_basis.getfloat(),
				self.uz_basis.getfloat())
		V  = bmath.Vector(
				self.vx_basis.getfloat(),
				self.vy_basis.getfloat(),
				self.vz_basis.getfloat())
		U.norm()
		V.norm()
		self.setBasisAxes(U,V)
		return (U, V)

	# ----------------------------------------------------------------------
	def setBasis(self, u, v):
		BasePlot.setFloat(self.ux_basis, u.x())
		BasePlot.setFloat(self.uy_basis, u.y())
		BasePlot.setFloat(self.uz_basis, u.z())

		BasePlot.setFloat(self.vx_basis, v.x())
		BasePlot.setFloat(self.vy_basis, v.y())
		BasePlot.setFloat(self.vz_basis, v.z())
		self.setBasisAxes(u,v)

	# ----------------------------------------------------------------------
	def setBasisAxes(self, U, V):
		Udir = U.direction()
		if Udir == "N": Udir="U"
		Vdir = V.direction()
		if Vdir == "N": Vdir="V"
		self.basisAxes.config(text="%s:%s"%(Udir,Vdir))

	# ----------------------------------------------------------------------
	def basis_xy(self):
		if self.basisAxes["state"]==DISABLED: return
		if len(self.widget["coord"].get())==3: self.widget["coord"].set("X-Y")
		self.setBasis(bmath.Vector.X, bmath.Vector.Y)

	# ----------------------------------------------------------------------
	def basis_xz(self):
		if self.basisAxes["state"]==DISABLED: return
		if len(self.widget["coord"].get())==3: self.widget["coord"].set("X-Z")
		self.setBasis(bmath.Vector.X, bmath.Vector.Z)

	# ----------------------------------------------------------------------
	def basis_yz(self):
		if self.basisAxes["state"]==DISABLED: return
		if len(self.widget["coord"].get())==3: self.widget["coord"].set("Y-Z")
		self.setBasis(bmath.Vector.Y, bmath.Vector.Z)

	# ----------------------------------------------------------------------
	def basis_swap(self):
		if self.basisAxes["state"]==DISABLED: return
		pc = self.widget["coord"].get()
		if len(pc)==3: self.widget["coord"].set("%s-%s"%(pc[-1],pc[0]))
		U, V = self.getBasis()
		self.setBasis(V, U)

	# ----------------------------------------------------------------------
	def basis_negu(self):
		if self.basisAxes["state"]==DISABLED: return
		U, V = self.getBasis()
		self.setBasis(-U, V)

	# ----------------------------------------------------------------------
	def basis_negv(self):
		if self.basisAxes["state"]==DISABLED: return
		U, V = self.getBasis()
		self.setBasis(U, -V)

	# ----------------------------------------------------------------------
	def reset_origin(self):
		self.setOrigin(0.0, 0.0, 0.0)

	# ----------------------------------------------------------------------
	def reset_basis(self): self.basis_xy()

	# ----------------------------------------------------------------------
	def reset_extend(self):
		BasePlot.setFloat(self.u_extend, 100.0)
		BasePlot.setFloat(self.v_extend, 100.0)

	# ----------------------------------------------------------------------
	def get_extend(self):
		if not self.flair.plotEngine:
			self.flair.notify("Get extends",
				"You need to plot first before getting the extends from plot window",
				tkFlair.NOTIFY_ERROR)
			return
		if self.u_extend["state"]==DISABLED and self.v_extend["state"]==DISABLED: return
		xlow, xhigh = self.flair.plotEngine.getRange("x")
		ylow, yhigh = self.flair.plotEngine.getRange("y")

		if xhigh is None: return
		if self.u_extend["state"]!=DISABLED:
			BasePlot.setFloat(self.u_extend, (xhigh-xlow)/2.0)
		if self.v_extend["state"]!=DISABLED:
			BasePlot.setFloat(self.v_extend, (yhigh-ylow)/2.0)

		# find the origin of axis
		x = (xhigh+xlow)/2.0
		y = (yhigh+ylow)/2.0

		# depending the plot configuration...
		O  = self.getOrigin()
		cd = self.widget["coord"].get()
		if cd == Plot.COORD[0] or cd == Plot.COORD[1]:
			# relative
			U, V = self.getBasis()
			if cd == Plot.COORD[0]:
				self.setOrigin(O + x*U + y*V)
			else:
				self.setOrigin(O + x*V + y*U)
		else:
			# absolute
			O[ ord(cd[ 0])-ord('X') ] = x
			O[ ord(cd[-1])-ord('X') ] = y
			self.setOrigin(O)

	# ----------------------------------------------------------------------
	# Calculate aspect ratio from extends
	# ----------------------------------------------------------------------
	def calcAspectRatio(self):
		extU = self.u_extend.getfloat()
		extV = self.v_extend.getfloat()
		try: ratio = extV / extU
		except: return
		self.plot_ratio.delete(0,END)
		self.plot_ratio.insert(0,ratio)

	# ----------------------------------------------------------------------
	def editPalette(self):
		self.flair.paletteDialog()

#===============================================================================
if __name__ == "__main__":
	tk = Tk()
	tkFlair.loadIcons()
	f = GeometryPlotFrame(tk)
	f.pack(expand=YES,fill=BOTH)
	tk.mainloop()
	tkFlair.delIcons()
