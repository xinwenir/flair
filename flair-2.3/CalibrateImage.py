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
# Date:	13-Jun-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import sys
import math
from log import say

try:
	from tkinter import *
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	import tkinter.messagebox as messagebox

import bmath
import tkFlair
import tkExtra
import Ribbon
try:
	import PIL.Image as Image
	import ZoomImage
except ImportError:
	Image = None

CURSOR = 3

MXi = 0
MYi = 1
MX  = 2		# Position in markers list
MY  = 3
MZ  = 4
L1  = 5
L2  = 6
TXT = 7

ACTION_NONE   = 0
ACTION_ADD    = 1
ACTION_EDIT   = 2
ACTION_DELETE = 3

_FORMATXYZ = "[%g, %g, %g]"

root = None

#===============================================================================
class Marker:
	idx   = 0
	color = "Magenta"
	font  = "Helvetica,14"
	width = 5.0

	#-----------------------------------------------------------------------
	def __init__(self, xi, yi, x, y, z):
		self.xi = float(xi)
		self.yi = float(yi)
		self.x  = x
		self.y  = y
		self.z  = z
		self.tag = "M%d"%(Marker.idx)
		Marker.idx += 1

	#-----------------------------------------------------------------------
	def create(self, canvas):
		self.canvas = canvas
		zoom = canvas.zoom
		xi = round(self.xi*zoom)
		yi = round(self.yi*zoom)
		self.line1 = canvas.create_line(
				xi - Marker.width,
				yi - Marker.width,
				xi + Marker.width + 1,
				yi + Marker.width + 1,
				fill=Marker.color,
				tags=self.tag)

		self.line2 = canvas.create_line(
				xi - Marker.width,
				yi + Marker.width,
				xi + Marker.width + 1,
				yi - Marker.width - 1,
				fill=Marker.color,
				tags=self.tag)

#		canvas.create_rectangle(( xi, yi, xi+2, yi+2), fill="#FFFF00")

		self.txt = canvas.create_text(
				xi, yi + Marker.width,
				text=str(self),
				fill=Marker.color,
				font=Marker.font,
				tags=self.tag,
				justify=CENTER, anchor=N)

	#-----------------------------------------------------------------------
	def zoom(self):
		zoom = self.canvas.zoom
		self.canvas.coords(self.line1,
				self.xi*zoom - Marker.width,
				self.yi*zoom - Marker.width,
				self.xi*zoom + Marker.width + 1,
				self.yi*zoom + Marker.width + 1)
		self.canvas.coords(self.line2,
				self.xi*zoom - Marker.width,
				self.yi*zoom + Marker.width + 1,
				self.xi*zoom + Marker.width + 1,
				self.yi*zoom - Marker.width)
		self.canvas.coords(self.txt,
				self.xi*zoom,
				self.yi*zoom + Marker.width)

	#-----------------------------------------------------------------------
	def delete(self):
		self.canvas.delete(self.tag)

	#-----------------------------------------------------------------------
	def __str__(self):
		return _FORMATXYZ%(self.x, self.y, self.z)

	#-----------------------------------------------------------------------
	def edit(self, master):
		self.toplevel = Toplevel(master)
		self.toplevel.title("Coordinates")

		Label(self.toplevel, text="X:").grid(row=0,column=0)
		Label(self.toplevel, text="Y:").grid(row=1,column=0)
		Label(self.toplevel, text="Z:").grid(row=2,column=0)

		self.xentry = tkExtra.FloatEntry(self.toplevel, background="White")
		self.yentry = tkExtra.FloatEntry(self.toplevel, background="White")
		self.zentry = tkExtra.FloatEntry(self.toplevel, background="White")

		self.xentry.grid(row=0, column=1, columnspan=3, sticky=EW)
		self.yentry.grid(row=1, column=1, columnspan=3, sticky=EW)
		self.zentry.grid(row=2, column=1, columnspan=3, sticky=EW)

		b = Button(self.toplevel, text="Ok",
				image=tkFlair.icons["ok"],
				compound=LEFT,
				command=self.ok,
				padx=2, pady=2)
		b.grid(row=3, column=2)
		b = Button(self.toplevel, text="Cancel",
				image=tkFlair.icons["x"],
				compound=LEFT,
				command=self.close,
				padx=2, pady=2)
		b.grid(row=3, column=3)
		self.toplevel.grid_columnconfigure(1, weight=1)

		self.xentry.insert(0, str(self.x))
		self.yentry.insert(0, str(self.y))
		self.zentry.insert(0, str(self.z))
		self.xentry.focus_set()

		self.toplevel.bind("<Return>",   self.ok)
		self.toplevel.bind("<KP_Enter>", self.ok)
		self.toplevel.bind("<Escape>",   self.close)
		self.toplevel.protocol("WM_DELETE_WINDOW", self.close)

		self.toplevel.transient(master)
		self.toplevel.wait_visibility()
		self.toplevel.grab_set()
		self.toplevel.wait_window()

	#-----------------------------------------------------------------------
	def ok(self, event=None):
		try: self.x = float(self.xentry.get())
		except: pass
		try: self.y = float(self.yentry.get())
		except: pass
		try: self.z = float(self.zentry.get())
		except: pass

		self.canvas.itemconfig(self.txt, text=str(self))
		self.close()

	#-----------------------------------------------------------------------
	def close(self, event=None):
		self.toplevel.destroy()
		del self.toplevel

#===============================================================================
# Image Calibrate class
#===============================================================================
class Calibrate(Toplevel):
	width  = -1
	height = -1

	#-----------------------------------------------------------------------
	def __init__(self, master, *options):
		Toplevel.__init__(self, master, *options)
		self.title("Calibrate Image")

		# ---------- create layout -------------
		#self.splitter = tkExtra.HSplitter(self, 250, True)
		#self.splitter.pack(expand=YES, fill=BOTH)

		# --- working canvas ---
		frame = Frame(self)
		frame.pack(side=TOP, expand=YES, fill=BOTH)

		self.canvas = ZoomImage.ZoomImage(frame,
				#width=640, height=480,
				takefocus=True)
		self.canvas.grid(row=0, column=0, sticky=NSEW)

		sbv = Scrollbar(frame, orient=VERTICAL)
		sbv.grid(row=0, column=1, sticky=NSEW)
		sbv.config(command=self.canvas.yview)

		sbh = Scrollbar(frame, orient=HORIZONTAL)
		sbh.grid(row=1, column=0, sticky=NSEW)
		sbh.config(command=self.canvas.xview)

		self.canvas.config( yscrollcommand=sbv.set,
				    xscrollcommand=sbh.set)
		frame.grid_columnconfigure(0, weight=1)
		frame.grid_rowconfigure(0, weight=1)

		# Override zooming
		self.canvas.bindMotion()
		self.canvas.bindArrows()
		self.canvas.bind('<Button-4>',  self.zoomIn)
		self.canvas.bind('<Button-5>',  self.zoomOut)
		self.canvas.bind('<Key-minus>', self.zoomOut)
		self.canvas.bind('<Key-equal>', self.zoomIn)
		self.canvas.bind('<Key-plus>',  self.zoomIn)
		self.canvas.bind('<Key-1>',     self.zoomOne)

		self.canvas.bind('<Button-1>',  self.click)
		self.canvas.bind('<Control-Button-1>', self.addMarkerAction)
		self.canvas.bind('<Double-1>',  self.editMarkerAction)
		self.canvas.bind('<Button-3>',  self.popup)
		self.canvas.bind('<Key-space>', self.popup)

		self.canvas.bind('<B1-Motion>', self.motion)
		self.canvas.bind('<ButtonRelease-1>', self.release)
		self.canvas.bind('<Enter>',	self.enter)
		self.canvas.bind('<Insert>',	self.addMarkerAction)
		self.canvas.focus_set()

		# --- canvas ---
		ribbon = Frame(self)
		ribbon.place(x=0, y=0)

		# ===
		group = Ribbon.LabelGroup(ribbon, "Edit")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		b = Ribbon.LabelButton(group.frame,
				text="Add",
				image=tkFlair.icons["add32"],
				compound=TOP,
				command=self.addMarker,
				padx=3,
				pady=0,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y)

		b = Ribbon.LabelButton(group.frame,
				text="Edit",
				image=tkFlair.icons["edit32"],
				compound=TOP,
				command=self.editMarker,
				padx=3,
				pady=0,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y)

		b = Ribbon.LabelButton(group.frame,
				text="Delete",
				image=tkFlair.icons["x32"],
				compound=TOP,
				command=self.delMarker,
				padx=3,
				pady=0,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y)

		# ===
		group = Ribbon.LabelGroup(ribbon, "View")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid2rows()

		b = Ribbon.LabelButton(group.frame,
				text="Zoom In",
				image=tkFlair.icons["zoom_in"],
				compound=LEFT,
				command=self.zoomIn,
				padx=3,
				pady=0,
				background=Ribbon._BACKGROUND)
		b.grid(row=0, column=0, stick=W)

		b = Ribbon.LabelButton(group.frame,
				text="Zoom Out",
				image=tkFlair.icons["zoom_out"],
				compound=LEFT,
				command=self.zoomOut,
				padx=3,
				pady=0,
				background=Ribbon._BACKGROUND)
		b.grid(row=1, column=0, stick=W)

		# ---
		b = Ribbon.LabelButton(group.frame,
				text="Info",
				image=tkFlair.icons["info32"],
				compound=TOP,
				command=self.showInfo,
				padx=3,
				pady=0,
				background=Ribbon._BACKGROUND)
		b.grid(row=0, rowspan=2, column=1, stick=NSEW)

		# ===
		group = Ribbon.LabelGroup(ribbon, "Confirm")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		b = Ribbon.LabelButton(group.frame,
				text="Ok",
				image=tkFlair.icons["ok32"],
				compound=TOP,
				command=self.ok,
				padx=3,
				pady=0,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y)

		b = Ribbon.LabelButton(group.frame,
				text="Cancel",
				image=tkFlair.icons["close32"],
				compound=TOP,
				command=self.close,
				padx=3,
				pady=0,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y)

		# Global variables
		self.imagefile = ""
		self.filename  = None

		#self.info = Frame(self.splitter.right())
		#self.info.pack(fill=BOTH)

		self.action    = ACTION_NONE
		self.w         = 5
		self.marker    = []
		self.point     = []
		self.closest   = None
		self.closest_x = 0.0
		self.closest_y = 0.0
		self.oldZoom   = 1.0
		self.rc        = False

#		self.protocol("WM_DELETE_WINDOW", self.close)

#		self.cursor = (self.canvas.create_line(-2, 0, 2, 0, fill="Pink"),
#			       self.canvas.create_line( 0,-2, 0, 2, fill="Pink"))

		# ---
#		viewmenu = Menu(self)
#		viewmenu.add_command(label="Zoom 1:1", underline=0,
#					accelerator="1", command=self.zoomOne)
#		viewmenu.add_command(label="Zoom In", underline=0,
#					accelerator="=", command=self.zoomIn)
#		viewmenu.add_command(label="Zoom Out", underline=0,
#					accelerator="-", command=self.zoomOut)

	#-----------------------------------------------------------------------
	def show(self):
		self.transient(self.master)
		self.deiconify()
		self.wait_visibility()
		if Calibrate.width > 0:
			self.geometry("%dx%d" % \
				(Calibrate.width, Calibrate.height))
		else:
			self.geometry("%dx%d" % \
				(self.winfo_width(), self.winfo_height()))
		self.wait_window()
		return self.rc

	#-----------------------------------------------------------------------
	def ok(self):
		self.rc = True
		self.close()

	#-----------------------------------------------------------------------
	def close(self):
		Calibrate.width  = self.winfo_width()
		Calibrate.height = self.winfo_height()
		self.destroy()
		if root: root.destroy()

	#-----------------------------------------------------------------------
	def showInfo(self):
		InfoDialog(self)

	#-----------------------------------------------------------------------
#	def printInfo(self):
#		markers = self.getMarkers()
#		self.matrix = compute(markers)
#		i = 1
#		say("##     x          y          z         xi       yi      err    dx     dy")
#		say("-- ---------- ---------- ---------- -------- -------- ------ ------ ------")
#		for xi,yi,x,y,z in markers:
#			pixel = self.matrix * bmath.Vector(x,y,z)
#			e = math.sqrt((pixel[0]-xi)**2 + (pixel[1]-yi)**2)
#			say("%2d %10g %10g %10g %8.4g %8.4g %6.3g %6.3g %6.3g" % \)
#				(i, x,y,z, xi,yi, e,(pixel[0]-xi),(pixel[1]-yi))
#			i += 1
#		say("matrix=\n",self.matrix)
#		say("Error=",error(markers, self.matrix))

	#-----------------------------------------------------------------------
	def loadImage(self, fn):
		self.config(cursor="watch")
		try:
			image = Image.open(fn)
			image.load()
		except:
			self.config(cursor="")
			messagebox.showerror("Error opening file",
				"Error opening image file %s"%(fn),
				parent=self)
			return
		self.imagefile = fn
		self.image     = image
		self.canvas.setImage(image)
		self.config(cursor="")

	# ----------------------------------------------------------------------
	def setImage(self, img):
		self.canvas.setImage(img)

	#-----------------------------------------------------------------------
	def addMarker(self):
		self.action = ACTION_ADD
		self.config(cursor="tcross")

	#-----------------------------------------------------------------------
	def editMarker(self):
		self.action = ACTION_EDIT
		self.config(cursor="question_arrow")

	#-----------------------------------------------------------------------
	def delMarker(self):
		self.action = ACTION_DELETE
		self.config(cursor="pirate")

	# ----------------------------------------------------------------------
	def setMarkers(self, markers):
		self.deleteAll()
		for xi,yi,x,y,z in markers:
			mark = Marker(xi,yi,x,y,z)
			self.marker.append(mark)
			mark.create(self.canvas)

	# ----------------------------------------------------------------------
	def getMarkers(self):
		markers = []
		for mark in self.marker:
			markers.append([mark.xi, mark.yi, mark.x, mark.y, mark.z])
		return markers

	# ----------------------------------------------------------------------
	def zoomMarkers(self):
		for mark in self.marker:
			mark.zoom()

	# ----------------------------------------------------------------------
	def find(self, x, y):
		"""find closest marker"""
		#items = self.canvas.find_closest(x, y,
		#		self.w*self.canvas.zoom,
		#		self.canvas.image)
		items = self.canvas.find_overlapping(x-1, y-1, x+2, y+2)

		for item in items:
			if item == self.canvas.image: continue

			for mark in self.marker:
				if item==mark.line1 or item==mark.line2 or item==mark.txt:
					tag = self.canvas.gettags(mark.line1)[0]
					return tag

		return None

	# ----------------------------------------------------------------------
	def _click(self, event):
		self.cx = self.canvas.canvasx(event.x)
		self.cy = self.canvas.canvasy(event.y)
		self.x = float(self.cx)/self.canvas.zoom
		self.y = float(self.cy)/self.canvas.zoom

		self.closest   = self.find(self.cx, self.cy)
		self.start_x   = event.x
		self.start_y   = event.y
		self.closest_x = event.x
		self.closest_y = event.y

		pos = "%6.1f %6.1f"%(self.x,self.y)
		self.clipboard_clear()
		self.clipboard_append(pos)

	# ----------------------------------------------------------------------
	def addMarkerAction(self, event):
		self.action = ACTION_NONE
		self._click(event)
		if self.marker:
			mark = self.marker[-1]
			x = mark.x
			y = mark.y
			z = mark.z
		else:
			x = y = z = 0.0
		mark = Marker(self.x,self.y,x,y,z)
		self.marker.append(mark)
		mark.create(self.canvas)
		mark.edit(self.canvas)

	# ----------------------------------------------------------------------
	def delMarkerAction(self, event):
		self.action = ACTION_NONE
		self._click(event)
		self.deleteMarker()

	# ----------------------------------------------------------------------
	def click(self, event):
		if self.canvas.tkimage is None: return
		self._click(event)

		if self.action == ACTION_ADD:
			self.addMarkerAction(event)
		elif self.action == ACTION_EDIT:
			self.editMarkerAction(event)
		elif self.action == ACTION_DELETE:
			self.delMarkerAction(event)
		#if self.closest is None:
		self.config(cursor="")

	# ----------------------------------------------------------------------
	def motion(self, event):
		if self.canvas.tkimage is None: return
		if self.closest is None: return
		if event.x == self.closest_x and event.y == self.closest_y:
			return

		# Moving a marker or an object
		self.canvas.move(self.closest,
			event.x-self.closest_x,
			event.y-self.closest_y)

		self.closest_x = event.x
		self.closest_y = event.y
		self.canvas.update_idletasks()

	# ----------------------------------------------------------------------
	def release(self, event):
		if self.canvas.tkimage is None: return
		if self.closest is None: return

		# Modify position of marker
		x = float(self.canvas.canvasx(event.x))/self.canvas.zoom
		y = float(self.canvas.canvasy(event.y))/self.canvas.zoom

		for mark in self.marker:
			if mark.tag == self.closest:
				mark.xi += (event.x-self.start_x)/self.canvas.zoom
				mark.yi += (event.y-self.start_y)/self.canvas.zoom
				return

	# ----------------------------------------------------------------------
#	def setCursor(self, x, y):
#		self.canvas.coords(self.cursor[0], x-CURSOR, y, x+CURSOR+1, y)
#		self.canvas.coords(self.cursor[1], x, y-CURSOR, x, y+CURSOR+1)

	# ----------------------------------------------------------------------
	def editMarkerAction(self, event=None):
		if self.closest is None: return
		# Edit marker
		for mark in self.marker:
			if mark.tag == self.closest:
				mark.edit(self.canvas)
				return

	# ----------------------------------------------------------------------
	def deleteMarker(self, event=None):
		if self.closest is None: return

		# Delete marker
		for mark in self.marker:
			if mark.tag == self.closest:
				mark.delete()
				self.marker.remove(mark)
				return

	# ----------------------------------------------------------------------
	def deleteAll(self):
		for mark in self.marker:
			mark.delete()
		del self.marker[:]

	# ----------------------------------------------------------------------
	def popup(self, event=None):
		if self.canvas.tkimage is None: return
		self._click(event)
#		self.setCursor(self.cx, self.cy)

		menu = Menu(self, tearoff=0)
		menu.add_command(label="Add",    command=self.addMarker)
		menu.add_command(label="Edit",   command=self.editMarker)
		menu.add_command(label="Delete", command=self.delMarker)
		#menu.add_separator()
		menu.tk_popup(event.x_root+1, event.y_root+1)

	# ----------------------------------------------------------------------
	def zoomIn(self, event=None):
		self.oldZoom = self.canvas.zoom
		self.canvas.zoomIn(event)
		self.zoomMarkers()

	# ----------------------------------------------------------------------
	def zoomOut(self, event=None):
		self.oldZoom = self.canvas.zoom
		self.canvas.zoomOut(event)
		self.zoomMarkers()

	# ----------------------------------------------------------------------
	def zoomOne(self, event=None):
		self.oldZoom = self.canvas.zoom
		self.canvas.zoomImage(1)
		self.zoomMarkers()

	# ----------------------------------------------------------------------
	def enter(self, event):
		self.canvas.focus_set()

#-------------------------------------------------------------------------------
class InfoDialog(Toplevel):
	def __init__(self, master, **kw):
		Toplevel.__init__(self, master, class_="InfoDialog", **kw)
		self.transient(master)

		# calibration data
		self.markers = master.getMarkers()

		self.listbox = tkExtra.MultiListbox(self,
					(("#",   2, None),
					 ("xi",  7, None),
					 ("yi",  7, None),
					 ("x",  10, None),
					 ("y",  10, None),
					 ("z",  10, None),
					 ("\u0394x",  7, None),
					 ("\u0394y",  7, None),
					 ("err", 7, None)),
					 background='White')
		self.listbox.pack(side=TOP, expand=YES, fill=BOTH)
		self.listbox.bindList("<ButtonRelease-1>", self.edit)
		self.listbox.bindList("<Return>",          self.edit)
		self.listbox.bindList("<KP_Enter>",        self.edit)
		self.listbox.bindLeftRight()
		self.listbox.setPopupMenu([
				("Edit",   0, self.edit,   tkFlair.icons["edit"]),
				("Add",	   0, self.add,    tkFlair.icons["add"]),
				("Delete", 0, self.delete, tkFlair.icons["x"]),
				("Clone",  0, self.clone,  tkFlair.icons["clone"])])

		# Total error
		frame = Frame(self)
		frame.pack(fill=X)
		Button(frame, image=tkFlair.icons["add"], command=self.add).pack(side=LEFT)
		Button(frame, image=tkFlair.icons["del"], command=self.delete).pack(side=LEFT)
		Button(frame, image=tkFlair.icons["clone"], command=self.clone).pack(side=LEFT)

		self.pixelerror = Label(frame, width=10, background="White")
		self.pixelerror.pack(side=RIGHT)
		Label(frame, text="Error Pixel:").pack(side=RIGHT)

		# M matrix
		frame = LabelFrame(self, text="Transformation: Real to Image plane")
		frame.pack(fill=X)
		self.mLabel = []
		for j in range(3):
			Label(frame, text="#%d"%(j+1)).grid(
				row=j+1, column=0, sticky=EW)
			Label(frame, text="#%d"%(j+1)).grid(
				row=0, column=j+1, sticky=EW)
			frame.grid_columnconfigure(j+1, weight=1)
			row = []
			for i in range(3):
				l = Label(frame, width=10, background="White", anchor=E)
				l.grid(row=j+1, column=i+1, padx=2, pady=2, sticky=EW)
				row.append(l)
			self.mLabel.append(row)

		# R matrix
		frame = LabelFrame(self, text="Image Plane to Pixel")
		frame.pack(fill=X)
		self.rLabel = []
		for j in range(2):
			Label(frame, text="#%d"%(j+1)).grid(row=j+1, column=0, sticky=EW)
			Label(frame, text="#%d"%(j+1)).grid(row=0, column=j+1, sticky=EW)
			frame.grid_columnconfigure(j+1, weight=1)
			row = []
			for i in range(3):
				l = Label(frame, width=10, background="White", anchor=E)
				l.grid(row=j+1, column=i+1, padx=2, pady=2, sticky=EW)
				row.append(l)
			self.rLabel.append(row)
		Label(frame, text="#3").grid(row=0, column=3, sticky=EW)
		frame.grid_columnconfigure(3, weight=1)

		# Buttons
		frame = Frame(self)
		frame.pack(side=BOTTOM, fill=X)
		Button(frame, text="Close", command=self.quit).pack(side=RIGHT)

		self.fill()

		# Window
		self.protocol("WM_DELETE_WINDOW", self.quit)
		try:
			self.wait_visibility()
			self.grab_set()
			self.wait_window()
		except TclError:
			pass

	#-------------------------------------------------------------------------------
	def fill(self):
		# Matrix
		R,M = compute(self.markers)
		for j in range(3):
			for i in range(3):
				if M is None:
					value = ""
				else:
					value = str(M[j][i])
				self.mLabel[j][i].config(text=value)
		for j in range(2):
			for i in range(3):
				if R is None:
					value = ""
				else:
					value = str(R[j][i])
				self.rLabel[j][i].config(text=value)

		# Error
		if M is not None:
			value = "%g"%(error(self.markers,R,M))
			self.pixelerror.config(text=value)
		else:
			self.pixelerror.config(text="")

		# Markers
		sel = self.listbox.curselection()
		act = self.listbox.index(ACTIVE)
		self.listbox.delete(0,END)
		i = 1
		for xi,yi,x,y,z in self.markers:
			if R is not None:
				uvw = M * bmath.Vector(x,y,z)
				uvw[2] = 1.0
				pixel = R * uvw
				dx = (pixel[0]-xi)
				dy = (pixel[1]-yi)
				e  = math.sqrt(dx**2 + dx**2)
				if abs(dx<1e-04): dx = 0.0
				if abs(dy<1e-04): dy = 0.0
				if abs(e <1e-04): e  = 0.0
				dx = "%-.4g"%(dx)
				dy = "%-.4g"%(dy)
				e  = "%-.4g"%(e)
			else:
				dx = ""
				dy = ""
				e  = ""
			try:   x  = "%g"%(x)
			except: pass
			try: y  = "%g"%(y)
			except: pass
			try: z  = "%g"%(z)
			except: pass
			try: xi = "%g"%(xi)
			except: pass
			try: yi = "%g"%(yi)
			except: pass
			self.listbox.insert(END, (i, xi,yi, x,y,z, dx, dy, e))
			i += 1

		for s in sel: self.listbox.selection_set(s)
		self.listbox.activate(act)

	#-------------------------------------------------------------------------------
	def edit(self, event=None):
		listbox = self.listbox.focus_get()
		nl = self.listbox.lists.index(listbox)
		if nl<1 or nl>5: return
		sel = listbox.curselection()
		if len(sel) != 1: return
		self.listbox.activate(sel[0])
		edit = tkExtra.InPlaceEdit(listbox)
		if edit:
			try: value = float(edit.value)
			except: return
			self.markers[int(sel[0])][nl-1] = value
			self.fill()
			self.master.setMarkers(self.markers)

	# ----------------------------------------------------------------------
	def add(self, event=None):
		self.markers.append([0,0,0,0,0])
		self.listbox.insert(END, (len(self.markers), "0", "0", "0", "0", "0", "", "", ""))
		self.listbox.selection_clear(0,END)
		self.listbox.activate(END)
		self.listbox.selection_set(END)

	# ----------------------------------------------------------------------
	def delete(self, event=None):
		lst = list(map(int,self.listbox.curselection()))
		lst.reverse()
		for i in lst:
			self.listbox.delete(i)
			del self.markers[i]
		self.fill()
		self.master.setMarkers(self.markers)

	# ----------------------------------------------------------------------
	def clone(self, event=None):
		sel = list(map(int, self.listbox.curselection()))
		self.listbox.selection_clear(0,END)
		for i in sel:
			row = self.listbox.get(i)
			row[0] = len(self.markers)
			self.markers.append(self.markers[i][:])
			self.listbox.insert(END, row)
			self.listbox.selection_set(END)
		self.fill()
		self.master.setMarkers(self.markers)

	# ----------------------------------------------------------------------
	def quit(self):
		self.destroy()

#-------------------------------------------------------------------------------
def compute(markers):
	if markers is None or len(markers)<3: return None,None

	# Find perpendicular (plot or projection) plane
	abcd = bmath.fitPlane([x[2:] for x in markers])
	if abcd:
		w = bmath.Vector(abcd[0], abcd[1], abcd[2])
		norm = w.normalize()
	else:
		return None,None

	# Create a 3D system u,v,w   with w=plane
	u = w.orthogonal()
	u.normalize()
	v = w ^ u

	M = bmath.Matrix(3)
	M.make(u,v,w)

	# Transformation from real(x,y,z) to plane(u,v,w) coordinates
	# [u,v,w] = M * [x,y,z]
	A = []
	B = []
	C = []
	# Find transformation from plane(u,v) coord to image(i,j) coordinates
	# ignoring u
	# Solve first for Xi,Yi = k*u + l*v + m
	for i,j,x,y,z in markers:
		uvw = M*bmath.Vector(x,y,z)
		uvw[2]=1.0	# ignore w set it to 1.0
		A.append(uvw)
		B.append([i])
		C.append([j])
	A = bmath.Matrix(A)
	B = bmath.Matrix(B)
	C = bmath.Matrix(C)
	try:
		RTi = bmath.solveOverDetermined(A, B)
	except:
		return None,None
	try:
		RTj = bmath.solveOverDetermined(A, C)
	except:
		return None,None

	# plane to image. Normally should be 2x3
	R = bmath.Matrix(3,type=1)
	for i in range(3):
		R[0][i] = RTi[i]
		R[1][i] = RTj[i]
	return R,M

#-------------------------------------------------------------------------------
def error(markers, R,M):
	# Find pixel error
	err = 0.0
	for i,j,x,y,z in markers:
		uvw = M * bmath.Vector(x,y,z)
		uvw[2] = 1.0
		ijk = R * uvw
		err += (ijk[0]-i)**2 + (ijk[1]-j)**2
	return math.sqrt(err/len(markers))

#-------------------------------------------------------------------------------
def realerror(markers, matrix):
	# Find real error
	err = 0.0
	invMatrix = matrix.clone()
	invMatrix.inv()
	for xi,yi,x,y,z in markers:
		r = invMatrix * bmath.Vector(xi,yi,0.0)
		err += (r[0]-x)**2 + (r[1]-y)**2 + (r[2]-z)**2
	return math.sqrt(err/len(markers))

if __name__ == "__main__":
	markers = [	(2410.5, 1566.5,  112.183, 1804.7, -170.583),
			(1400.5, 1569.5, -228.707, 1804.7,  323.172),
			(1400.5,  641.5, -228.707, 2354.7,  323.172),
			(2412.5,  641.5,  112.183, 2354.7, -170.583),
			(2042.0, 1569.0,  -12.81,  1804.7,   10.46) ]
	R,M = compute(markers)
	say("Error=",error(markers, R,M))
#	say(planeFit(markers))
#	root = Tk()
#	tkFlair.loadIcons()

#	calibrate = Calibrate(None)
#	root.withdraw()

#	if len(sys.argv)>1:
#		calibrate.loadImage(sys.argv[1])
#		calibrate.setMarkers(MARKERS)	# TESTING...
#	root.geometry("1024x640")
#	calibrate.show()
#	root.mainloop()
#	tkFlair.delIcons()
#	say(calibrate.getMarkers())
