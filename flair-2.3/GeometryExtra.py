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
# Date:	27-Jan-2012

import string

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkExtra
import Input
import Manual
import tkFlair
import Unicode

from bmath import *
from Constants import *

FMT = "%.10g"

#===============================================================================
# Base frame for origin manipulation
#===============================================================================
class _BaseFrame(Frame):
	def __init__(self, master, geometry, *kw):
		Frame.__init__(self, master, *kw)
		self.geometry = geometry
		self.viewer   = geometry.viewer
		self._viewer  = geometry.viewer._viewer

	# ----------------------------------------------------------------------
	# Return the apply button
	# ----------------------------------------------------------------------
	def applyButton(self, msg, compound=LEFT):
		b = Button(self, text="Apply",
				image=tkFlair.icons["ok"],
				compound=compound,
				command=self.apply,
				padx=1,
				pady=0)
		tkExtra.Balloon.set(b, msg)
		return b

	# ----------------------------------------------------------------------
	# fill information from viewer on every draw
	# ----------------------------------------------------------------------
	def fill(self):  pass

	# ----------------------------------------------------------------------
	# load information from project/input
	# ----------------------------------------------------------------------
	def load(self):	 pass

	# ----------------------------------------------------------------------
	# apply changes to viewer
	# ----------------------------------------------------------------------
	def apply(self, event=None):
		self.viewer.focus_set()

	# ----------------------------------------------------------------------
	# update viewer
	# ----------------------------------------------------------------------
	def draw(self, invalid, expose):
		self.viewer.draw(invalid, expose)

#===============================================================================
class _NavigationFrame(_BaseFrame):
	name = "Navigation"
	icon = "navigation"

	def __init__(self, master, geometry, *kw):
		_BaseFrame.__init__(self, master, geometry, *kw)

		# back
		self.backwardButton = Button(self, image=tkFlair.icons["left"],
				takefocus=0,
				command=self.viewer.historyBackward,
				state=DISABLED,
				padx=0, pady=0)
		self.backwardButton.pack(side=LEFT)
		tkExtra.Balloon.set(self.backwardButton, "Move back one view [<]")

		# forward
		self.forwardButton = Button(self, image=tkFlair.icons["right"],
				takefocus=0,
				command=self.viewer.historyForward,
				state=DISABLED,
				padx=0, pady=0)
		self.forwardButton.pack(side=LEFT)
		tkExtra.Balloon.set(self.forwardButton, "Go forward one view [>]")

		# reload
		b = Button(self, image=tkFlair.icons["refresh"],
				takefocus=0,
				command=geometry.redraw,
				padx=0, pady=0)
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b, "Refresh display")

		# sep
		Label(self, image=tkFlair.icons["sep"]).pack(side=LEFT, padx=3)

		# layers
		b = Button(self, image=tkFlair.icons["layers"],
				takefocus=0,
				command=self.geometry.layersDialog,
				padx=0, pady=0)
		tkExtra.Balloon.set(b, "Configure the layers [Ctrl-L]")
		b.pack(side=LEFT)

		# layer Combo
		self.layerCombo = tkExtra.Combobox(self, label=True, width=10,
				takefocus=0,
				command=self.geometry.changeLayer)
		self.layerCombo.pack(side=LEFT, padx=1, pady=0, expand=YES, fill=X)
		self.layerCombo.select(1)

		# Bookmarks
		b = Button(self, image=tkFlair.icons["star"], command=self.addPlot)
		b.pack(side=LEFT)

		tkExtra.Balloon.set(b, "Bookmark current view as a new geometry plot")

		self.plotCombo = tkExtra.Combobox(self, False, width=9,
				command=self.changePlot)
		self.plotCombo.pack(side=LEFT, fill=X, expand=YES)
		self.plotCombo._text.config(background="White")
		self.plotCombo._text.bind("<FocusOut>", self.renamePlot)
		self.plotCombo._text.bind("<Return>",   self.returnHit)
		self.plotCombo._text.bind("<KP_Enter>", self.returnHit)

		b = Button(self, image=tkFlair.icons["del"], command=self.delPlot)
		b.pack(side=LEFT)
		tkExtra.Balloon.set(b, "Delete selected bookmark")

		self.lastPlot = self.geometry.name


	# ----------------------------------------------------------------------
	def historyButton(self, left, right):
		if left:
			self.backwardButton["state"]=NORMAL
		else:
			self.backwardButton["state"]=DISABLED

		if right:
			self.forwardButton["state"]=NORMAL
		else:
			self.forwardButton["state"]=DISABLED

	# ----------------------------------------------------------------------
	def _loadListbox(self):
		lst = []
		for plot in self.geometry.project().plots:
			if plot.type == "Geometry":
				lst.append(plot.name)
		self.plotCombo.fill(lst)

	# ----------------------------------------------------------------------
	def load(self):
		self._loadListbox()
		if self.plotCombo.get() != self.lastPlot:
			self.plotCombo.set(self.geometry.name)

	# ----------------------------------------------------------------------
	def apply(self, event=None):
		_BaseFrame.apply(self)
		self.geometry.loadPlot(self.plotCombo.get())
		self.draw(True, True)

	# ----------------------------------------------------------------------
	def changePlot(self):
		if self.lastPlot == self.plotCombo.get(): return
		self.lastPlot = self.plotCombo.get()
		self.apply()

	# ----------------------------------------------------------------------
	def addPlot(self, event=None):
		for i in range(1,1000):
			name = "Geometry%02d"%(i)
			for plot in self.geometry.project().plots:
				if plot.type == "Geometry" and plot.name == name:
					break
			else:
				break
		self.geometry.savePlot(name)
		self.plotCombo.insert(END,name)
		self.lastPlot = name
		self.plotCombo.set(name)
		self.geometry.parent.reloadOriginFrame()

	# ----------------------------------------------------------------------
	def delPlot(self, event=None):
		name = self.plotCombo.get()
		for i, plot in enumerate(self.geometry.project().plots):
			if plot.name == name and plot.type == "Geometry":
				del self.geometry.project().plots[i]
				self.plotCombo.delete(i)
				break
		self._loadListbox()
		self.lastPlot = self.geometry.name
		self.plotCombo.set(self.lastPlot)
		self.geometry.parent.reloadOriginFrame()

	# ----------------------------------------------------------------------
	def renamePlot(self, event=None):
		if self.lastPlot == self.plotCombo.get(): return
		if self.focus_get() == self.plotCombo._text: return
		for plot in self.geometry.project().plots:
			if plot.name == self.lastPlot and plot.type == "Geometry":
				plot.name = self.plotCombo.get()
				self._loadListbox()
				self.lastPlot = plot.name
				self.plotCombo.set(self.lastPlot)
				self.geometry.parent.reloadOriginFrame()
				return

	# ----------------------------------------------------------------------
	def returnHit(self, event=None):
		self.viewer.focus_set()


#===============================================================================
class _OriginFrame(_BaseFrame):
	name = "Origin"
	icon = "center"
	def __init__(self, master, geometry, *kw):
		_BaseFrame.__init__(self, master, geometry, *kw)

		Label(self,text="x:").pack(side=LEFT)
		self.xentry = tkExtra.VectorEntry(self, width=6, background="White")
		self.xentry.pack(side=LEFT, fill=X, expand=YES)
		self.xentry.bind("<Return>", self.apply)
		self.xentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.xentry, "Origin X-coordinate or all coordinates (x y z) as a vector")

		Label(self,text="y:").pack(side=LEFT, )
		self.yentry = tkExtra.FloatEntry(self, width=6, background="White")
		self.yentry.pack(side=LEFT, fill=X, expand=YES)
		self.yentry.bind("<Return>", self.apply)
		self.yentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.yentry, "Origin Y-coordinate")

		Label(self,text="z:").pack(side=LEFT, )
		self.zentry = tkExtra.FloatEntry(self, width=6, background="White")
		self.zentry.pack(side=LEFT, fill=X, expand=YES)
		self.zentry.bind("<Return>", self.apply)
		self.zentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.zentry, "Origin Z-coordinate")

		self.applyButton("Set Origin").pack(side=RIGHT)

	# ----------------------------------------------------------------------
	def fill(self):
		xyz = self._viewer.origin()
		self.xentry.set(FMT%(xyz[0]))
		self.yentry.set(FMT%(xyz[1]))
		self.zentry.set(FMT%(xyz[2]))

	# ----------------------------------------------------------------------
	def apply(self, event=None):
		_BaseFrame.apply(self)
		x,y,z = self._viewer.origin()
		self.xentry.split((self.yentry, self.zentry))
		x = self.xentry.getfloat(x)
		y = self.yentry.getfloat(y)
		z = self.zentry.getfloat(z)
		self._viewer.viewOto0()
		self._viewer.origin(x,y,z)
		self.draw(True, True)

#===============================================================================
class _ExtendsFrame(_BaseFrame):
	name = "Extends"
	icon = "expand"
	def __init__(self, master, geometry, *kw):
		_BaseFrame.__init__(self, master, geometry, *kw)

		Label(self,text="\u0394u:").pack(side=LEFT)
		self.duentry = tkExtra.FloatEntry(self, width=6, background="White")
		self.duentry.pack(side=LEFT, fill=X, expand=YES)
		self.duentry.bind("<Return>", self.apply)
		self.duentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.duentry, "Extends on U")

		Label(self,text="\u0394v:").pack(side=LEFT)
		self.dventry = Label(self, background="White", anchor=W, relief=SUNKEN)
		self.dventry.pack(side=LEFT, fill=X, expand=YES)
		tkExtra.Balloon.set(self.dventry, "Extends on V")

		Label(self,text="Aspect(X/Y):").pack(side=LEFT)
		self.aspectentry = tkExtra.FloatEntry(self, width=6, background="White")
		self.aspectentry.pack(side=LEFT, fill=X, expand=YES)
		self.aspectentry.bind("<Return>", self.apply)
		self.aspectentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.aspectentry, "Aspect Ratio")

		self.applyButton("Set extends and aspect ratio").pack(side=RIGHT)

	# ----------------------------------------------------------------------
	def fill(self):
		x1,y1,x2,y2 = self._viewer.extends()
		z = self._viewer.zoom()
		du = (x2-x1)/2.0/z
		w,h = self._viewer.size()
		dv = float(h)/float(w) * du
		self.duentry.set(FMT%(du))
		self.dventry.config(text=FMT%(dv))

	# ----------------------------------------------------------------------
	def apply(self, event=None):
		_BaseFrame.apply(self)
		try: du = float(self.duentry.get())
		except: return
		try:
			aspect = float(self.aspectentry.get())
		except:
			aspect = 1.0
		du = 2.0*du
		self._viewer.extends(-du, -du, du, du)
		self._viewer.zoom(2.0)
		self._viewer.aspect(aspect)
		self.draw(True, True)

#===============================================================================
class _MoveFrame(_BaseFrame):
	name = "Move"
	icon = "pan"
	def __init__(self, master, geometry, *kw):
		_BaseFrame.__init__(self, master, geometry, *kw)

		Label(self,text="+u:").pack(side=LEFT)
		self.uentry = tkExtra.VectorEntry(self, width=6, background="White")
		self.uentry.pack(side=LEFT, fill=X, expand=YES)
		self.uentry.bind("<Return>", self.apply)
		self.uentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.uentry, "Shift Origin on U or by vector (ux, uy, uz)")

		Label(self,text="+v:").pack(side=LEFT)
		self.ventry = tkExtra.FloatEntry(self, width=6, background="White")
		self.ventry.pack(side=LEFT, fill=X, expand=YES)
		self.ventry.bind("<Return>", self.apply)
		self.ventry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.ventry, "Shift Origin on V")

		Label(self,text="+w:").pack(side=LEFT)
		self.wentry = tkExtra.FloatEntry(self, width=6, background="White")
		self.wentry.pack(side=LEFT, fill=X, expand=YES)
		self.wentry.bind("<Return>", self.apply)
		self.wentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.wentry, "Shift Origin on W")

		self.applyButton("Move origin by").pack(side=RIGHT)

	# ----------------------------------------------------------------------
	def apply(self, event=None):
		_BaseFrame.apply(self)
		u,v = self._viewer.offset()
		self.uentry.split((self.ventry,self.wentry))
		u += self.uentry.getfloat()
		v += self.ventry.getfloat()
		self._viewer.offset(u,v)

		w = self.wentry.getfloat()
		if abs(w)>1e-10:
			matrix = Matrix(self._viewer.matrix())
			self._viewer.matrix(matrix*Matrix.translate(0.0, 0.0, w))
			self.draw(True, True)
		else:
			self.draw(False, True)

#===============================================================================
class _BasisFrame(_BaseFrame):
	name = "Basis"
	icon = "axes"
	def __init__(self, master, geometry, *kw):
		_BaseFrame.__init__(self, master, geometry, *kw)

		row = 0
		Label(self,text="ux:").grid(row=row, column=0)
		self.uxentry = tkExtra.VectorEntry(self, width=5, background="White")
		self.uxentry.grid(row=row, column=1, sticky=EW)
		self.uxentry.bind("<Return>", self.apply)
		self.uxentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.uxentry, "Enter U-axis x-direction cosine or all (ux uy uz) as a vector")

		Label(self,text="uy:").grid(row=row, column=2)
		self.uyentry = tkExtra.FloatEntry(self, width=5, background="White")
		self.uyentry.grid(row=row, column=3, sticky=EW)
		self.uyentry.bind("<Return>", self.apply)
		self.uyentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.uyentry, "Enter U-axis y-direction cosine")

		Label(self,text="uz:").grid(row=row, column=4)
		self.uzentry = tkExtra.FloatEntry(self, width=5, background="White")
		self.uzentry.grid(row=row, column=5, sticky=EW)
		self.uzentry.bind("<Return>", self.apply)
		self.uzentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.uzentry, "Enter U-axis z-direction cosine")

		f = Frame(self)
		f.grid(row=row, column=6, sticky=EW)
		Button(f, text="x-y", padx=3, pady=1, command=self.orientXY).pack(side=LEFT)
		Button(f, text="x-z", padx=2, pady=1, command=self.orientXZ).pack(side=LEFT)
		Button(f, text="y-z", padx=2, pady=1, command=self.orientYZ).pack(side=LEFT)

		f = Frame(self)
		f.grid(row=row, column=7, sticky=EW)
		Button(f, image=tkFlair.icons["rotate_left"], command=self.rotWp, pady=1).pack(side=LEFT)
		Button(f, image=tkFlair.icons["rotate_vleft"],command=self.rotVp, pady=1).pack(side=LEFT)
		Button(f, image=tkFlair.icons["rotate_up"],   command=self.rotUp, pady=1).pack(side=LEFT)

		row = 1
		Label(self,text="vx:").grid(row=row, column=0)
		self.vxentry = tkExtra.VectorEntry(self, width=5, background="White")
		self.vxentry.grid(row=row, column=1, sticky=EW)
		self.vxentry.bind("<Return>", self.apply)
		self.vxentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.vxentry, "Enter V-axis x-direction cosine or all (vx vy vz) as a vector")

		Label(self,text="vy:").grid(row=row, column=2)
		self.vyentry = tkExtra.FloatEntry(self, width=5, background="White")
		self.vyentry.grid(row=row, column=3, sticky=EW)
		self.vyentry.bind("<Return>", self.apply)
		self.vyentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.vyentry, "Enter V-axis y-direction cosine")

		Label(self,text="vz:").grid(row=row, column=4)
		self.vzentry = tkExtra.FloatEntry(self, width=5, background="White")
		self.vzentry.grid(row=row, column=5, sticky=EW)
		self.vzentry.bind("<Return>", self.apply)
		self.vzentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.vyentry, "Enter V-axis z-direction cosine")

		f = Frame(self)
		f.grid(row=row, column=6, sticky=EW)
		Button(f, text="-v",   padx=2, pady=1, command=self.orient_V).pack(side=LEFT)
		Button(f, text="-u",   padx=2, pady=1, command=self.orient_U).pack(side=LEFT)
		Button(f, text="swap", padx=2, pady=1, command=self.orientSWAP).pack(side=LEFT)

		f = Frame(self)
		f.grid(row=row, column=7, sticky=EW)
		Button(f, image=tkFlair.icons["rotate_right"], command=self.rotWm, pady=1).pack(side=LEFT)
		Button(f, image=tkFlair.icons["rotate_vright"],command=self.rotVm, pady=1).pack(side=LEFT)
		Button(f, image=tkFlair.icons["rotate_down"],  command=self.rotUm, pady=1).pack(side=LEFT)

		self.applyButton("Set basis orientation", compound=TOP).grid(row=0, column=13, rowspan=2)

		self.grid_columnconfigure(1,weight=1)
		self.grid_columnconfigure(3,weight=1)
		self.grid_columnconfigure(5,weight=1)

	# ----------------------------------------------------------------------
	def fill(self):
		matrix = Matrix(self._viewer.matrix())
		self.orientSet(	matrix[0][0], matrix[1][0], matrix[2][0],
				matrix[0][1], matrix[1][1], matrix[2][1])

	# ----------------------------------------------------------------------
	def rotUp(self): self.rotate( 90.0,  0.0,  0.0)
	def rotUm(self): self.rotate(-90.0,  0.0,  0.0)
	def rotVp(self): self.rotate(  0.0, 90.0,  0.0)
	def rotVm(self): self.rotate(  0.0,-90.0,  0.0)
	def rotWp(self): self.rotate(  0.0,  0.0,-90.0)
	def rotWm(self): self.rotate(  0.0,  0.0, 90.0)

	# ----------------------------------------------------------------------
	def rotate(self, ru, rv, rw):
		_BaseFrame.apply(self)
		self._viewer.viewOto0()
		rotated = False
		ang = math.radians(ru)
		if abs(ang)>1e-10:
			self._viewer.matrix(Matrix(self._viewer.matrix())*Matrix.rotX(ang))
			rotated = True

		ang = math.radians(rv)
		if abs(ang)>1e-10:
			self._viewer.matrix(Matrix(self._viewer.matrix())*Matrix.rotY(ang))
			rotated = True

		ang = math.radians(rw)
		if abs(ang)>1e-10:
			self._viewer.matrix(Matrix(self._viewer.matrix())*Matrix.rotZ(ang))
			rotated = True
		if rotated: self.draw(True, True)

	# ----------------------------------------------------------------------
	def apply(self, event=None):
		_BaseFrame.apply(self)
		try: u,v,w = self.orientNorm()
		except: return

		xyz    = self._viewer.origin()
		matrix = self._viewer.matrix()
		for i in range(3):
			matrix[i][0] = u[i]
			matrix[i][1] = v[i]
			matrix[i][2] = w[i]
		self._viewer.matrix(matrix)
		self._viewer.origin(xyz)
		self.draw(True, True)

	# ----------------------------------------------------------------------
	def orientSet(self, ux=None, uy=None, uz=None, vx=None, vy=None, vz=None):
		if ux is not None: self.uxentry.set(str(ux))
		if uy is not None: self.uyentry.set(str(uy))
		if uz is not None: self.uzentry.set(str(uz))
		if vx is not None: self.vxentry.set(str(vx))
		if vy is not None: self.vyentry.set(str(vy))
		if vz is not None: self.vzentry.set(str(vz))

	# ----------------------------------------------------------------------
	def orientXY(self, event=None):
		self.orientSet(1.0,0.0,0.0, 0.0,1.0,0.0)
		self.apply()

	# ----------------------------------------------------------------------
	def orientXZ(self, event=None):
		self.orientSet(1.0,0.0,0.0, 0.0,0.0,1.0)
		self.apply()

	# ----------------------------------------------------------------------
	def orientYZ(self, event=None):
		self.orientSet(0.0,1.0,0.0, 0.0,0.0,1.0)
		self.apply()

	# ----------------------------------------------------------------------
	def orientSWAP(self, event=None):
		self.orientSet(	self.vxentry.get(), self.vyentry.get(), self.vzentry.get(),
				self.uxentry.get(), self.uyentry.get(), self.uzentry.get())
		self.apply()

	# ----------------------------------------------------------------------
	# Negate U vector
	# ----------------------------------------------------------------------
	def orient_U(self, event=None):
		try: ux = -float(self.uxentry.get())
		except: ux = None
		try: uy = -float(self.uyentry.get())
		except: uy = None
		try: uz = -float(self.uzentry.get())
		except: uz = None
		self.orientSet(	ux, uy, uz )
		self.apply()

	# ----------------------------------------------------------------------
	# Negate V vector
	# ----------------------------------------------------------------------
	def orient_V(self, event=None):
		try: vx = -float(self.vxentry.get())
		except: vx = None
		try: vy = -float(self.vyentry.get())
		except: vy = None
		try: vz = -float(self.vzentry.get())
		except: vz = None
		self.orientSet(	None, None, None, vx, vy, vz )
		self.apply()

	# ----------------------------------------------------------------------
	def orientNorm(self, event=None):
		# try it as a vector
		self.uxentry.split((self.uyentry, self.uzentry))
		self.vxentry.split((self.vyentry, self.vzentry))

		# Get vectors
		u = Vector(self.uxentry.getfloat(),
			   self.uyentry.getfloat(),
			   self.uzentry.getfloat())

		v = Vector(self.vxentry.getfloat(0.0),
			   self.vyentry.getfloat(1.0),
			   self.vzentry.getfloat(0.0))

		# normalize vectors
		u.norm()
		v.norm()

		# find w
		w = u.cross(v)
		if w.length() < 0.01:
			# try another vector
			w = u.cross(Vector.Z)
			if w.length() < 0.01:
				w = u.cross(Vector.X)
		w.norm()
		v = w.cross(u)
		v.norm()

		self.orientSet(u[0], u[1], u[2], v[0], v[1], v[2])
		return u,v,w

#===============================================================================
class _EulerFrame(_BaseFrame):
	name = "Euler"
	icon = "rotate"
	def __init__(self, master, geometry, *kw):
		_BaseFrame.__init__(self, master, geometry, *kw)

		Label(self,text="Rx:").pack(side=LEFT)
		self.rxentry = tkExtra.VectorEntry(self, width=6, background="White")
		self.rxentry.pack(side=LEFT, fill=X, expand=YES)
		self.rxentry.bind("<Return>", self.apply)
		self.rxentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.rxentry, "Enter Rx Euler angle or all (Rx Ry Rz)")

		Label(self,text="Ry:").pack(side=LEFT)
		self.ryentry = tkExtra.FloatEntry(self, width=6, background="White")
		self.ryentry.pack(side=LEFT, fill=X, expand=YES)
		self.ryentry.bind("<Return>", self.apply)
		self.ryentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.ryentry, "Enter Ry Euler angle")

		Label(self,text="Rz:").pack(side=LEFT)
		self.rzentry = tkExtra.FloatEntry(self, width=6, background="White")
		self.rzentry.pack(side=LEFT, fill=X, expand=YES)
		self.rzentry.bind("<Return>", self.apply)
		self.rzentry.bind("<KP_Enter>", self.apply)
		tkExtra.Balloon.set(self.rzentry, "Enter Rz Euler angle")

		self.applyButton("Set orientation using Euler angles").pack(side=RIGHT)

	# ----------------------------------------------------------------------
	def fill(self):
		matrix = Matrix(self._viewer.matrix())
		rxyz = matrix.getEulerRotation()
		self.rxentry.set(FMT%(math.degrees(rxyz[0])))
		self.ryentry.set(FMT%(math.degrees(rxyz[1])))
		self.rzentry.set(FMT%(math.degrees(rxyz[2])))

	# ----------------------------------------------------------------------
	def apply(self, event=None):
		_BaseFrame.apply(self)
		self._viewer.viewOto0()

		self.rxentry.split((self.ryentry, self.rzentry))
		rx = math.radians(self.rxentry.getfloat())
		ry = math.radians(self.ryentry.getfloat())
		rz = math.radians(self.rzentry.getfloat())

		xyz = self._viewer.origin()
		self._viewer.matrix(Matrix.eulerRotation(rx,ry,rz))
		self._viewer.origin(xyz)
		self.draw(True, True)

#===============================================================================
class _TransformFrame(_BaseFrame):
	name = "Transform"
	icon = "ROT-DEFI"
	def __init__(self, master, geometry, *kw):
		_BaseFrame.__init__(self, master, geometry, *kw)
		Label(self, text="Rotdefi:").pack(side=LEFT)
		self.transform = tkExtra.Combobox(self, width=10)
		self.transform.pack(side=LEFT,fill=X, expand=YES)
		self.applyButton("Transform using ROT-DEFI").pack(side=RIGHT)

	# ----------------------------------------------------------------------
	def fill(self):
		if self.geometry.parent.editor.editType() == Input.Card.REGION:
			rotdefi = self.geometry.parent.editor.editObject.cards[0]["@rotdefi"]
			if rotdefi: self.transform.set(rotdefi)

	# ----------------------------------------------------------------------
	def load(self):
		trans = self.geometry.input().rotdefiList(True)
		self.transform.fill(trans)

	# ----------------------------------------------------------------------
	def apply(self, event=None):
		_BaseFrame.apply(self)
		trans = self.geometry.input().getTransformation(str(self.transform.get()))
		self._viewer.viewOto0()
		matrix = Matrix(self._viewer.matrix())
		self._viewer.matrix(trans * matrix)
		self.draw(True, True)

#===============================================================================
# Projection Frame
#===============================================================================
class ProjectionFrame(Frame):
	def __init__(self, master, geometry, *kw):
		Frame.__init__(self, master, *kw)
		self.comboButton = Button(self,image=tkFlair.icons["axes"],
				text=" "+Unicode.BLACK_DOWN_POINTING_TRIANGLE,
				compound=LEFT,
				command=self.change,
				padx=0, pady=0)
		self.comboButton.grid(row=0, column=0, sticky=N+EW)
		tkExtra.Balloon.set(self.comboButton, "Origin and projection [o]")

		#b = Button(self, text="Apply",
		#		image=tkFlair.icons["ok"],
		#		compound=LEFT,
		#		command=self.apply)
		#b.grid(row=0, column=2, padx=2, sticky=EW)
		#tkExtra.Balloon.set(b, "Apply operation")

		self._frames = []
		self._frames.append(_NavigationFrame(self,geometry))
		self._frames.append(_OriginFrame(self,    geometry))
		self._frames.append(_ExtendsFrame(self,   geometry))
		self._frames.append(_MoveFrame(self,      geometry))
		self._frames.append(_BasisFrame(self,     geometry))
		self._frames.append(_EulerFrame(self,     geometry))
		self._frames.append(_TransformFrame(self, geometry))
		self.lastFrame = None
		self.show(self._frames[0])

		self.grid_columnconfigure(1, weight=1)
		return

	# ----------------------------------------------------------------------
	def change(self):
		menu = Menu(self, tearoff=0)
		for frame in self._frames:
			menu.add_command(label=frame.name,
					image=tkFlair.icons[frame.icon],
					compound=LEFT,
					command=lambda s=self,f=frame:s.show(f))
		menu.tk_popup(self.comboButton.winfo_rootx(),
			      self.comboButton.winfo_rooty() + self.comboButton.winfo_height())

	# ----------------------------------------------------------------------
	def show(self, frame):
		if self.lastFrame: self.lastFrame.grid_forget()
		self.comboButton["image"] = tkFlair.icons[frame.icon]
		frame.grid(row=0, column=1, sticky=EW)
		frame.fill()
		self.lastFrame = frame

	# ----------------------------------------------------------------------
	def load(self):
		for frame in self._frames:
			frame.load()

	# ----------------------------------------------------------------------
	def apply(self, event=None):
		self.lastFrame.apply()

#===============================================================================
# Errors Frame
#===============================================================================
class ErrorsFrame(Frame):
	def __init__(self, master, editor):
		Frame.__init__(self, master)
		self.editor = editor

		# --------------------------------------------------------------
		Button(self, text="Refresh", image=tkFlair.icons["refresh"],
			compound=LEFT, takefocus=FALSE,
			command=self.refresh).pack(side=TOP, fill=X)

		# --------------------------------------------------------------
		sbh = Scrollbar(self, orient=HORIZONTAL)
		sbh.pack(side=BOTTOM, fill=X)
		sbv = Scrollbar(self, orient=VERTICAL)
		sbv.pack(side=RIGHT, fill=Y)

		self.text = Text(self, width=80, cursor="arrow", wrap="none")
		self.text.pack(side=LEFT, expand=YES, fill=BOTH)
		self.text.bind("<Button-1>",  self.button1)
		self.text.bind("<Escape>",    self.escape)
		self.text.bind("<B1-Motion>", self.donothing)
		self.text.bind("<Double-1>",  self.donothing)
		self.text.bind("<F1>",        self.help)

		self.text.config(yscrollcommand=sbv.set, xscrollcommand=sbh.set, state=DISABLED)
		sbh.config(command=self.text.xview)
		sbv.config(command=self.text.yview)

		# Remember error-list
		self._hash   = [0, 0, 0, 0]
		self._errors = [0, 0, 0, 0]
		self._cache  = [[], [], [], []]
		self._show   = [True, True, True, True, True]
		self._select = ""
		self._expand = {}

	# ----------------------------------------------------------------------
	# Refresh all errors
	# ----------------------------------------------------------------------
	def refresh(self, event=None):
		for i,name in enumerate(VIEWNAME):
			viewer = self.editor.viewers[i]
			self.readErrors(i,viewer)
		self.fill()

	# ----------------------------------------------------------------------
	# Read errors from viewport
	# ----------------------------------------------------------------------
	def readErrors(self, idx, viewer):
		self._hash[idx]   = viewer.error("hash")
		self._errors[idx] = viewer.error("n")
		lst = self._cache[idx]
		del lst[:]
		for i in range(1,min(MAXERRORS, self._errors[idx]+1)):
			err = viewer.error("get",i)
			if err is None: break
			lst.append(err)

	# ----------------------------------------------------------------------
	# Fill text with error messages
	# ----------------------------------------------------------------------
	def fill(self, event=None):
		self.text.config(state=NORMAL)
		self.text.delete("0.0",END)
		self.text.tag_delete("*")

		for i,name in enumerate(VIEWNAME):
			viewer = self.editor.viewers[i]
			if viewer.error("hash") != self._hash[i]:
				self.readErrors(i,viewer)
			self.insertErrors(i, name)

		ne = len(self.editor.viewers.errors)
		nw = len(self.editor.viewers.warnings)

		if self._show[4]:
			header = "%s Input [Errors:%d, Warnings:%d]\n"%(Unicode.BLACK_DOWN_POINTING_TRIANGLE, ne, nw)
		else:
			header = "%s Input [Errors:%d, Warnings:%d]\n"%(Unicode.BLACK_RIGHT_POINTING_TRIANGLE, ne, nw)

		self.text.insert(END, header)
		self.text.tag_add("Input", "end - 2 lines", "end -1 lines")
		self.text.tag_configure("Input",   background="#FFFF70")

		if self._show[4]:
			if self.editor.viewers.errors:
				tag = "Errors"
				self.text.insert(END, "Errors:\n",tag)
				self.text.tag_configure(tag, foreground="DarkRed")

				for i,err in enumerate(self.editor.viewers.errors):
					tag = "Error.%d"%(i)
					self.text.insert(END, "%2d: %s"%(i+1,err), tag)
#					if not i%2: self.text.tag_config(tag, background="DarkGray")

			if self.editor.viewers.warnings:
				tag = "Warnings"
				if self.editor.viewers.errors: self.text.insert(END, "\n")
				self.text.insert(END, "Warnings:\n",tag)
				self.text.tag_configure(tag, foreground="DarkBlue")
				for i,err in enumerate(self.editor.viewers.warnings):
					tag = "Warning.%d"%(i)
					self.text.insert(END, "%2d: %s"%(i+1,err), tag)
#					if not i%2: self.text.tag_config(tag, background="DarkGray")

		self.text.config(state=DISABLED)

	# ----------------------------------------------------------------------
	# Insert errors in Text from a specific viewport
	# ----------------------------------------------------------------------
	def insertErrors(self, idx, name):
		if self._show[idx]:
			header = "%s %s [%d]\n"%(Unicode.BLACK_DOWN_POINTING_TRIANGLE, name, self._errors[idx])
		else:
			header = "%s %s [%d]\n"%(Unicode.BLACK_RIGHT_POINTING_TRIANGLE, name, self._errors[idx])

		self.text.insert(END, header)
		self.text.tag_add(name,"end - 2 lines", "end -1 lines")
		self.text.tag_configure(name, background="#%06X"%(VIEWCOLOR[idx]))

		if self._errors[idx] and self._show[idx]:
			viewer = self.editor.viewers[idx]
			viewer.error("show", MAXERRORS)
			viewer.expose()
			i = 0
			for err in self._cache[idx]:
				i += 1
				x,y,z,xmin,xmax,ymin,ymax,body,regIn,regOut = err
				tag = "%s.%d"%(name,i)
				if tag in self._expand:
					self.text.insert(END, "-%2d:  %8s %8s %8s\n" \
								 "   body: %s\n" \
								 "   Reg+: %s\n" \
								 "   Reg-: %s\n" \
							% (i, format(x,8),format(y,8),format(z,8), body,
							   " ".join(regOut,","),
							   " ".join(regIn,",")),
							tag)
				else:
					self.text.insert(END, "+%2d:  %8s %8s %8s\n" \
							% (i, format(x,8),format(y,8),format(z,8)),
							tag)
#				if not i%2: self.text.tag_config(tag, background="DarkGray")

		self.text.insert(END, "\n")

	# ----------------------------------------------------------------------
	def button1(self, event=None):
		index = self.text.index("@%d,%d"%(event.x,event.y))
		row,col = list(map(int,index.split(".")))
		tags = self.text.tag_names(index)
		if not tags: return

		for tag in tags:
			if tag == SEL: continue

			try:
				name,err = tag.split(".")
				err = int(err)
			except ValueError:
				name = tag
				err = -1

			try:
				idx = VIEWNAME.index(name)
			except ValueError:
				if name == "Input":
					idx = 4
				else:
					idx = -1

			if err>=0:
				if col<3:
					if tag in self._expand:
						del self._expand[tag]
					else:
						self._expand[tag] = True
					self.fill()
				self.select(tag)
				if col>=3 and idx>=0:
					# Zoom on selected error
					self.editor.viewers.frames[idx].viewer.zoomOnError(err)
			else:
				self._show[idx] = not self._show[idx]
				self.fill()
			break

		self.text.see(index)
		return "break"

	# ----------------------------------------------------------------------
	def escape(self, event=None):
		self.select("")

	# ----------------------------------------------------------------------
	def select(self, tag):
		self.text.tag_config(self._select,
					background="",
					foreground="")
		self._select = tag
		if self._select:
			self.text.tag_config(self._select,
					background=self.text["selectbackground"],
					foreground=self.text["selectforeground"])

	# ----------------------------------------------------------------------
	def donothing(self, event=None):
		return "break"

	# ----------------------------------------------------------------------
	def help(self, event=None):
		Manual.show("F4.8")

	# ----------------------------------------------------------------------
#	def write(self, event=None):
#		filename = bFileDialog.asksaveasfilename(master=self,
#			title="Save errors",
#			filetypes=[("Errors file",("*.err")),
#				("Text file",("*.txt")),
#				("All","*")])

#		if len(filename) > 0:
#			f = open(filename,"w")
#			for i in range(self.listbox.size()):
#				err,x,y,z,body,regIn,regOut = self.listbox.get(i)
#				f.write("Error: %d\n"%(err))
#				f.write("\tPosition:   %g,%g,%g\n"%(x,y,z))
#				f.write("\tBody:       %s\n"%(body)
#				f.write("\tRegion-In:  %s\n"%(regIn)
#				f.write("\tRegion-Out: %s\n"%(regOut)
#			f.close()
#			messagebox.showinfo("Write Errors",
#				"Errors are writen in file %s"%(filename),
#				parent=self)
#

#===============================================================================
# Volume calculation frame/dialog
#===============================================================================
class VolumeDialog(Toplevel):
	width  = -1
	height = -1

	def __init__(self, master, geometry, project, regions):
		Toplevel.__init__(self, master)
		self.transient(master)
		self.title("Volume calculation")
		self.project  = project
		self.geometry = geometry
		self.regions  = regions

		self.listbox = tkExtra.MultiListbox(self,
					(("Region", 10, None),
					 ("Volume", 15, None),
					 ("Error",  15, None),
					 ("Iteration", 8, None)),
					height=10)
		self.listbox.pack(side=TOP, fill=BOTH, expand=YES)

		f = Frame(self)
		f.pack(side=BOTTOM, fill=X)
		b = Button(f, text="Next", command=self.__next__)
		b.pack(side=LEFT)
		b = Button(f, text="Stop", command=self.stop)
		b.pack(side=LEFT)

		self.startCalculation()
		self.monitor()

		self.lift()
		self.wait_visibility()
		self.grab_set()
		self.wait_window()

	#-----------------------------------------------------------------------
	def __next__(self):
		pass

	#-----------------------------------------------------------------------
	def stop(self):
		self.destroy()

	#-----------------------------------------------------------------------
	def monitor(self):
		self.after(VolumeDialog.TIMER, self.monitor)

	#-----------------------------------------------------------------------
	# Start calculation in the background
	#-----------------------------------------------------------------------
	def startCalculation(self):
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.type() == Input.Card.REGION and card.get(SELECT,0)&BIT_SELECT:
				say("\nRegion:", card.name())
				n = 32
				try:
					for i in range(18):
						vol = self.viewers._geometry.volume(card[ID], n)
						say("Geometry Volume:", vol, n)
						#vol = self.viewers[0].volume(card[ID], n)
						#say("Viewer Volume:", vol, n)
						n <<= 1
				except KeyboardInterrupt:
					pass
