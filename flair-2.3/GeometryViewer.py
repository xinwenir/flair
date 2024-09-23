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
import re
import math
import time
import struct
import string
import sys
import traceback
import pdb
#import platform

from stat import *
from log import say
from operator import attrgetter

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkExtra
import Unicode

from bmath import *
import bFileDialog

import csg
import undo
import Input
import Project
import tkFlair
import platform

import Palette
import CalibrateImage

import GeometryBody
import GeometryTool
import GeometryExtra
import GeometryLayers

from Constants import *
from Layout import _ASSIGNMA_field

geoviewer = None
try:
	import geoviewer
	geoviewer_error = ""
except:
	geoviewer_error = sys.exc_info()[1]

try:
	import PIL.Image
	import PIL.ImageTk	# Check if ImageTk exists to be used by ZoomImage
	have_PIL = True
except:
	have_PIL = False

if geoviewer is None:
	if platform.machine()=="x86_64":
		try:
			import geoviewer64 as geoviewer
			geoviewer_error = ""
		except:
			pass
	else:
		try:
			import geoviewer32 as geoviewer
			geoviewer_error = ""
		except:
			pass

# General variables
MOVESTEP   = 10.0
MOVESTEP2  = 100.0
ROTSTEP    = math.pi / 36.0	#  5deg
ROTSTEP2   = math.pi / 2.0	# 90deg
MAXORBIT   = 100.0		# maximum orbit depth * window dimensions
ZOOM       = 1.20
SZOOM      = 1.02
DEFZOOM    = 4.0
MAXZOOM    = 50.0
TIMER      = 250
TIMERDRAW  = 50
TIMERPAN   = 25
FRICTION   = math.pow(1.0-0.0035, float(TIMERPAN))

# Masks
SHIFT_MASK   = 1
CONTROL_MASK = 4
ALT_MASK     = 8

PROJECTION_NOTHING    = 0	# No projection information
PROJECTION_SPAWN      = 1
PROJECTION_START      = 2	# start projection (after stop)
PROJECTION_CONICS     = 3	# Transform bodies and prepare conics
PROJECTION_LOCATION   = 4	# Calculate the bodies&region location
PROJECTION_INTERSECT  = 5	# Intersect all bodies
PROJECTION_SCAN       = 6	# Scan segments
PROJECTION_DRAW       = 7	# Late drawing
PROJECTION_FINISHED   = 10	# Projection ended

UVTITLE	= {	( "X","Y")	: "Front",
		("-X","Y")	: "Back",

		( "Z","Y")	: "Left",
		("-Z","Y")	: "Right",

		( "Z","X")	: "Top",
		("-Z","X")	: "Bottom"
	}

UNIT_MATRIX  = Matrix(4,type=1)
FRONT_MATRIX = UNIT_MATRIX
BACK_MATRIX  = Matrix.rotY(math.pi)
TOP_MATRIX   = Matrix( [[0.0, 1.0, 0.0, 0.0],
			[0.0, 0.0, 1.0, 0.0],
			[1.0, 0.0, 0.0, 0.0],
			[0.0, 0.0, 0.0, 1.0]])
BOTTOM_MATRIX = TOP_MATRIX * Matrix.rotY(math.pi)
RIGHT_MATRIX  = Matrix.rotY(math.pi/2.0)
LEFT_MATRIX   = Matrix.rotY(-math.pi/2.0)
FMT = "%.10g"

DEF_CURSOR = "left_ptr"
MOUSE_CURSOR = {
	ACTION_SELECT        : DEF_CURSOR,
	ACTION_SELECT_REGION : DEF_CURSOR,
	ACTION_SELECT_ZONE   : DEF_CURSOR,
	ACTION_SELECT_BODIES : DEF_CURSOR,

	ACTION_AREA_SELECT   : "right_ptr",
	ACTION_ZONE          : "center_ptr",
	ACTION_ZONEPAINT     : "spraycan",	# "pencil"
	ACTION_PEN           : "pencil",
	ACTION_PAINT         : "spraycan",
	ACTION_INFO          : "tcross",	# "target"
	ACTION_VIEW_CENTER   : "target",	# "target"

	ACTION_PAN           : "fleur",
	ACTION_ORBIT         : "exchange",
	ACTION_ZOOM_IN       : "sizing",
	ACTION_ZOOM_OUT      : "sizing",
	ACTION_ZOOM_ON       : "sizing",

	ACTION_VIEW_CENTER   : "cross",
	ACTION_VIEW_MOVE     : "fleur",
	ACTION_VIEW_ROTATE   : "exchange",

	ACTION_ADD           : "tcross",
	ACTION_ADD_NEXT      : "tcross",

	ACTION_MOVE          : "hand1",
	ACTION_ROTATE        : "exchange",

	ACTION_EDIT          : "pencil",
}

# Viewers
VIEWFORMAT = {	VIEWNAME[0] : (1.0, 1.0, 0.5),
		VIEWNAME[1] : (1.0, 0.0, 0.5),
		VIEWNAME[2] : (0.0, 0.5, 1.0),
		VIEWNAME[3] : (0.0, 0.5, 0.0),
		"2x2"       : (0.5, 0.5, 0.5),
		"2x1"       : (1.0, 0.5, 0.5),
		"1x2"       : (0.5, 1.0, 1.0),
		"1+2"       : (0.5, 1.0, 0.5),
		"2+1"       : (0.5, 0.5, 1.0) }
VIEWORDER  =  VIEWNAME[:]
VIEWORDER.extend(["2x2", "2x1", "1x2", "1+2", "2+1", "custom"])

# Color definitions
_BACKGROUND_COLOR = 0x707070	# geometry.background
_SELECT_COLOR     = 0xFF00FF	# geometry.select
_LOCK_COLOR       = 0xFF8033	# geometry.lock
_FREEZE_COLOR     = 0x00FFFF	# geometry.freeze
_WIREFRAME_COLOR  = 0x007070	# geometry.wireframe
_VISIBLE_COLOR    = 0x7F007F	# geometry.visible
_REGION_COLOR     = 0x000000	# geometry.region
_ZONE_COLOR       = 0x606060	# geometry.zone
_ERROR_COLOR      = 0xFF0000	# geometry.error
_VERTEX_COLOR     = 0xFFFF00	# geometry.vertex
_LABEL_COLOR      = 0xA05010	# geometry.label
_LATTICE_COLOR    = 0xFFFFFF	# geometry.lattice
_BODY_BBIN_COLOR  = 0xE6E600	# geometry.bodyBBoxInColor
_BODY_BBOUT_COLOR = 0x00FF00	# geometry.bodyBBoxOutColor
_ZONE_BB_COLOR    = 0x00FFFF	# geometry.zoneBBoxColor
_REGION_BB_COLOR  = 0xFF8000	# geometry.regionBBoxColor

#_REGERROR_COLOR  = 0xFFB0B0	# geometry.region.error

_TITLE_COLOR      = 0xD0D000	# geometry.title
_GRIDTEXT_COLOR   = 0x505050	# geometry.gridlabel
#_COORD_COLOR      = # geometry.coordinate

# Configuration variables
TRACKBALLSIZE  = 0.8
CORES          = 0		# CPU/cores to use per viewport
LAPTOP         = False		# Laptop mode swap button-2 and button-3
INVERTWHEEL    = False		# Invert zooming with wheel
PINZOOM        = False		# Zoom in as google, center on mouse location
SNAPDISTANCE   = 4		# Snapping distance in pixels
CURSORSIZE     = 40		# Grab/Rotation cursor size
CROSSHAIR      = 10		# Cross hair size
_PROFILE       = False

# ------------------------------------------------------------------------------
# Variables to retrieve when manipulating an object
_OBJECT_VARS = {
	"!point" : (("x",  1), ("y",  2), ("z",  3)),
	"!arrow" : (("x",  1), ("y",  2), ("z",  3),
		    ("dx", 7), ("dy", 8), ("dz", 9)),
	"!light" : (("x",  1), ("y",  2), ("z",  3),
		    ("dx", 7), ("dy", 8), ("dz", 9)),
	"!mesh"  : (("x",  1), ("y",  2), ("z",  3)),
	"!ruler" : (("x",  1), ("y",  2), ("z",  3),
		    ("xe", 7), ("ye", 8), ("ze", 9),
		    ("xa",10), ("ya",11), ("za",12))
}

# Sequence of items to manipulate on ACTION_ADD_NEXT
_ADD_SEQUENCE = {	# use a list as a sequence since python 2.4 doesn't have tuple.index()
	"!arrow" : [ 2, ],
	"!ruler" : [ 2, ],
	"!light" : [ 2, ],
	"BOX"    : [ 10, 4 ],
	"ELL"    : [ 10, 2 ],
	"PLA"    : [ 20, ],
	"RAW"    : [ 11, 15 ],
	"RCC"    : [ 10, 1 ],
	"RPP"    : [ 18, 6 ],
	"REC"    : [ 10, 1 ],
	"SPH"    : [  1, ],
	"TRC"    : [ 10, 4, 1 ],
	"WED"    : [ 11, 15 ],
	"XCC"    : [  1, ],
	"YCC"    : [  1, ],
	"ZCC"    : [  1, ],
	"XEC"    : [  1, 2 ],
	"YEC"    : [  1, 2 ],
	"ZEC"    : [  1, 2 ],
}

# ------------------------------------------------------------------------------
def mouseCursor(action):
	return MOUSE_CURSOR.get(action, DEF_CURSOR)

# ------------------------------------------------------------------------------
# Convert EMFCUT what1 to kinetic energy
# ------------------------------------------------------------------------------
def _emfwhat1(x):
	if x<0.0: return abs(x)
	return max(0.0,x-0.51099906000000E-03)

# ------------------------------------------------------------------------------
# return a signed rotdefi index
# ------------------------------------------------------------------------------
def signedRotdefi(rotdefi, idx):
	if idx and idx[0]=='-':
		return -rotdefi.get(idx[1:],0)
	else:
		return rotdefi.get(idx,0)

# ------------------------------------------------------------------------------
# Prepare the region expression
# ------------------------------------------------------------------------------
def regionExp(expr):
	zones = []

	for zone in expr:
		# Check mode
		if "(" in zone:
			mode = "RPN"
			csg.exp2rpn(zone)
			if csg.check(zone):
				zones.append(("RPN", zone))
			else:
				raise Exception("Invalid RPN region expression")
		else:
			plus  = []
			minus = []
			prev  = "|"	# Just to mark the start
			for token in zone:
				if token not in ("+", "-"):
					if prev in ("+", "|"):
						plus.append(token)
					else:
						minus.append(token)
				else:
					if prev in ("+", "-"):
						raise Exception("Invalid region expression")
				prev = token

			plus.append("$")
			plus.extend(minus)
			zones.append(("STD",plus))
	return zones

# ------------------------------------------------------------------------------
# Convert a list of STD expression to a zone string
# ------------------------------------------------------------------------------
def list2Zone(zone):
	szone = ""
	op = " +"
	for b in zone:
#		if b=="$":
		if b.decode("UTF-8")=="$":
			op=" -"
		else:
			szone += op
			szone += b.decode("UTF-8")
	if szone: szone = szone[1:]
	return szone

# ------------------------------------------------------------------------------
# Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
# if we are away from the center of the sphere.
# ------------------------------------------------------------------------------
def projection2Sphere(r, x, y):
	d = sqrt(x**2 + y**2)
	if d < r:				# Inside sphere
		return sqrt(r**2 - d**2)
	else:					# On hyperbola
#		t = r / 1.41421356237309504880
#		return t*t / d
		return 0.0

# ------------------------------------------------------------------------------
# Ok, simulate a track-ball.  Project the points onto the virtual
# trackball, then figure out the axis of rotation, which is the cross
# product of P1 P2 and O P1 (O is the center of the ball, 0,0,0)
# Note:  This is a deformed trackball-- is a trackball in the center,
# but is deformed into a hyperbolic sheet of rotation away from the
# center.  This particular function was chosen after trying out
# several variations.
#
# It is assumed that the arguments to this routine are in the range (-1.0 ... 1.0)
# ------------------------------------------------------------------------------
# using quaternions
def trackball(p1x, p1y, p2x, p2y):
	if (abs(p1x-p2x)<=1.0e-5 and abs(p1y-p2y)<=1.0e-5):
		# Zero rotation
		return Quaternion([0.0, 0.0, 0.0, 1.0])

	# First, figure out z-coordinates for projection of P1 and P2 to
	# deformed sphere
	p1 = Vector(p1x, p1y, projection2Sphere(TRACKBALLSIZE, p1x, p1y))
	p2 = Vector(p2x, p2y, projection2Sphere(TRACKBALLSIZE, p2x, p2y))

	a = p2.cross(p1)		# rotation axis
	a.norm()

	# Figure out how much to rotate around that axis.
	d = p1-p2
	t = d.length() / (2.0*TRACKBALLSIZE)

	# Avoid problems with out-of-control values...
	if   t >  1.0:	t = 1.0
	elif t < -1.0:	t = -1.0
	phi = 2.0 * asin(t)		# how much we rotate on this axis

	return Quaternion(a, phi)

# using matrices
def trackball2(p1x, p1y, p2x, p2y, axis):
	if (abs(p1x-p2x)<=1.0e-5 and abs(p1y-p2y)<=1.0e-5):
		return None, None

	# First, figure out z-coordinates for projection of P1 and P2 to
	# deformed sphere
	p1 = Vector(p1x, p1y, projection2Sphere(TRACKBALLSIZE, p1x, p1y))
	p2 = Vector(p2x, p2y, projection2Sphere(TRACKBALLSIZE, p2x, p2y))

	if axis is not None:
		# find projection of p1,p2 on the perpendicular plane to axis (passing from 0)
		pp1 = p1 - (axis*p1)*axis
		pp2 = p2 - (axis*p2)*axis

		# Then the cross product should be parallel to axis
		# and the length to the sin()
		pp1.norm()
		pp2.norm()
		d = pp1 - pp2

		# Figure out how much to rotate around that axis.
		axis = pp2 ^ pp1
		l = axis.norm()
		if abs(l)<0.0001: axis = Vector.Z
		t = d.length() / 2.0

	else:
		axis = p2 ^ p1		# rotation axis
		l = axis.norm()
		if abs(l)<0.0001: axis = Vector.Z

		# Figure out how much to rotate around that axis.
		d = p1-p2
		t = d.length() / (2.0*TRACKBALLSIZE)

	# Avoid problems with out-of-control values...
	if   t >  1.0:	t = 1.0
	elif t < -1.0:	t = -1.0
	phi = 2.0 * asin(t)		# how much we rotate on this axis

	return axis, phi

#===============================================================================
# Geometry viewer frame
#===============================================================================
class GeometryViewer(Frame):
	def __init__(self, master, frame, parent, **kw):
		"""
		Initialize GeometryViewer
			master is the GeometryFrame
			parent is the Geometry4Frame
		"""
		Frame.__init__(self, master, **kw)
		self._viewer   = geoviewer.Viewer(parent._geometry, self._w, self.tk.interpaddr())
		self._viewer.set("cores", CORES)
		self.parent    = parent
		self.frame     = frame
		self.editor    = parent.editor
		self._cursor   = DEF_CURSOR
		self._geometry = parent._geometry
		self.config(background="", cursor=self._cursor)

		self.bindings()
		self.viewports  = []

		self.action     = ACTION_SELECT	# current action
		self.selected   = None		# item selected
		self.lateDraw   = False		# enable late drawing
		self.curHistory = -1
		self.history    = []
		self._save      = None		# General purpose variable
						# typically holds the mouse position
		self._closest   = None		# Tuple containing the closest mouse hit

		# FIXME maybe move to parent!
		self._objects     = []		# Selected moveable objects
		self._objectsid   = []		# @see _select()
		self._objectLight = False	# if object has a light selected @see _select
		self._bodies      = []		# Selected moveable bodies
		self._bodiesid    = []
		self.__opt        = 0		# During adding option to modify
		self.x            = 0		# starting position of the last action or segment
		self.y            = 0		#             -//-
		self._x           = 0		# Mouse position when clicked
		self._y           = 0		#             -//-
		self.__x          = 0		# starting position during adding
		self.__y          = 0		#             -//-
		self._orbitDepth  = 0.		# distance of the last centered <B2> object for orbit
		self._msg         = ""
		self._clone       = False	# Clone or Move (Control was pressed?)

		self.resetAnchor()

		if self._viewer.get("depth")<24:
			self.flair().notify("Insufficient Display color depth",
				"Display color depth is %d bits while 24 bits required" \
				% (self._viewer.get("depth")),
				tkFlair.NOTIFY_ERROR)

	# ----------------------------------------------------------------------
	def bindings(self):
		self.bind("<Enter>",		self.enter)
		self.bind("<Expose>",		self._viewer.expose)
		self.bind("<Configure>",	self.configure)
		#self.bind("<Visibility>",	self.visibility)
		#self.bind("<Property>",	self.property)
		self.bind("<<Project>>",	self._project)
		self.bind("<<EndProjection>>",	self.endProjection)
		self.bind("<<EndDraw>>",	self.endDraw)

		self.bind("<Motion>",		self.motion)

		self.bind("<Button-1>",		self.button1)
		self.bind("<Button-2>",		self.button2)
		self.bind("<Button-3>",		self.button3)

#		self.bind("<Double-1>",		self.double1)

		self.bind("<B1-Motion>",	self.buttonMotion)
		self.bind("<B2-Motion>",	self.buttonMotion)
		self.bind("<B3-Motion>",	self.buttonMotion)

		self.bind("<ButtonRelease-1>",	self.buttonRelease)
		self.bind("<ButtonRelease-2>",	self.buttonRelease)
		self.bind("<ButtonRelease-3>",	self.buttonRelease)

		if INVERTWHEEL:
			self.bind("<Button-4>",		self.cursorOrZoomOut)
			self.bind("<Button-5>",		self.cursorOrZoomIn)
		else:
			self.bind("<Button-4>",		self.cursorOrZoomIn)
			self.bind("<Button-5>",		self.cursorOrZoomOut)
		self.bind("<MouseWheel>",	self.wheel)

		self.bind("<Control-Button-4>",	self.panBack)
		self.bind("<Control-Button-5>",	self.panFront)

		#self.bind("<Shift-Button-4>",	self.panLeft)
		#self.bind("<Shift-Button-5>",	self.panRight)

		#self.bind("<Shift-Control-Button-4>",	self.panTop)
		#self.bind("<Shift-Control-Button-5>",	self.panBottom)

		self.bind("<KeyPress-Control_R>",   self.controlPress )
		self.bind("<KeyPress-Control_L>",   self.controlPress )
		self.bind("<KeyRelease-Control_R>", self.controlRelease )
		self.bind("<KeyRelease-Control_L>", self.controlRelease )

		self.bind("<Key>",                  self.handleKey)

		# Some global bindings
		self.bind("<Control-Key-l>",	self.frame.layersDialog)
		self.bind("<Control-Key-m>",	self.editor.gotoMaterial)
		self.bind("<Control-Key-n>",	self.toNotes)
		self.bind("<Control-Key-p>",	self.export)
		self.bind("<Control-Key-t>",	self.editor.transformGeometry)
		self.bind("<Control-Key-V>",	self.editor.pasteSpecial)
		self.bind("<F11>",		self.editor.toggleSplitter)
		if Input._developer:
			self.bind("<Control-Key-b>",	self.parent.autoBodyDialog)
			self.bind("<Control-Key-m>",	self.editor.mesh2Bodies)

		# Timers
		self._projId = None	# Projection
		self._drawId = None	# Late drawing
		self._panId  = None	# Panning

	# ----------------------------------------------------------------------
	def flair(self):	return self.editor.flair
	def project(self):	return self.parent.project
	def input(self):	return self.parent.project.input

	# ----------------------------------------------------------------------
	def addViewport(self, viewer, color):
		id = len(self.viewports)
		self._viewer.viewport(id, "viewer", viewer._viewer)
		self._viewer.viewport(id, "color",  color)
		self.viewports.append(viewer)

	# ----------------------------------------------------------------------
	def set(self, var, value=None):
		if value is None:
			return self._viewer.set(var)
		else:
			self._viewer.set(var,value)

	# ----------------------------------------------------------------------
	def _select(self, objects=None, bodies=None):
		if objects:
			self._objects = objects
			self._objectsid = [x[ID] for x in self._objects]
			for card in objects:
				if card.tag == "!light":
					self._objectLight = True
					break
			else:
				self._objectLight = False
		else:
			del self._objects[:]
			del self._objectsid[:]

		if bodies:
			self._bodies = bodies
			self._bodiesid = [x[ID] for x in self._bodies]
		else:
			del self._bodies[:]
			del self._bodiesid[:]

		return self._objects or self._bodies

	# ----------------------------------------------------------------------
	def enter(self, event):
		if not self.editor.isactive(): return
		# Steal keyboard focus when mouse enters only if window has the
		try: focus = self.focus_get()
		except: return
		if focus is None: return

		# Avoid stealing the focus from the Projection Dialog
		if isinstance(focus, Entry):
			widget = focus
			while True:
				widget = widget.master
				if widget is None: break
				if widget is self.master:
					return

		# Steal the focus for all other widgets withing the same toplevel
		if focus.winfo_toplevel() == self.winfo_toplevel():	# for V2
			self.focus_set()

	# ----------------------------------------------------------------------
	def expose(self, event=None):
		self._viewer.expose()

	# ----------------------------------------------------------------------
	def configure(self, event):
		w,h = self._viewer.size()
		expose = (event.width<=w) and (event.height<=h)
		self._viewer.configure(event.width, event.height)
		self.draw(False, expose)

	# ----------------------------------------------------------------------
	#def visibility(self, event):
		#say("*-*", self._viewer.title(), "visibility")

	# ----------------------------------------------------------------------
	#def property(self, event):
		#say("*-*", self._viewer.title(), "property")

	# ----------------------------------------------------------------------
	def _project(self, event=None):
		self.draw(True)

	# ----------------------------------------------------------------------
	def draw(self, invalid=False, expose=True, history=True):
		#if self.frame.name == "Magenta":
		#say("*-*", self._viewer.title(), "draw", invalid,
		#		self._viewer.invalid(), self._viewer.zoom())
		#traceback.print_stack()
		# Kill any previous request to draw
		if self._drawId is not None:
			self.after_cancel(self._drawId)
			self._drawId = None

		if invalid or self._viewer.invalid() or self._viewer.zoom()>MAXZOOM:
			# VERY VERY DANGEROUS!!
			# especially in the case that there is an EndProjection as
			# a pending task
			#if expose: self.update_idletasks()
			self._viewer.calcWindow(DEFZOOM)
			self._viewer.viewOto0()
			if history and self._panId is None: self.historyAdd()

			if self._projId is not None: self.after_cancel(self._projId)
			self._viewer.project(True)
			self._projId = self.after(TIMER, self.projectAfter)

			# Set as grid label the closest axis match
			axes = self._viewer.grid("axes")
			self._viewer.grid("label", axes)
			saxes = UVTITLE.get(axes)
			if saxes is None: saxes = "%s:%s"%axes
			self._viewer.title(saxes)
			self._viewer.draw()

		else:
			self._viewer.draw()
			if expose: self._viewer.expose()
			if self.lateDraw and \
			   self._viewer.state()==PROJECTION_FINISHED and \
			   self._panId is None:
				self._drawId = self.after(TIMERDRAW, self.startDrawThread)

		# expose other viewports
		if self._viewer.projectionChanged():
			for v in self.viewports:
				v.expose()
			self._viewer.projectionChanged(False)

		self.master.fillOrigin()

	# ----------------------------------------------------------------------
	def drawFast3D(self):
		if self._viewer.get("3D"):
			self._viewer.draw(False, DRAW_FAST)
			self._viewer.expose()

	# ----------------------------------------------------------------------
	def drawFast3DOtherViewports(self):
		for view in self.viewports:
			if view._viewer.get("3D"):
				view._viewer.draw(False, DRAW_FAST)

	# ----------------------------------------------------------------------
	def projectAfter(self):
		s = self._viewer.draw()
		#say("*-*", self._viewer.title(), "projectAfter",s)
		self._viewer.expose()
		if s == PROJECTION_FINISHED:
			if self._viewer.error("n"):
				# FIX ME we should clear all previous errors
				self.parent.event_generate("<<Errors>>")

			if self.lateDraw:
				self._drawId = self.after(TIMERDRAW, self.startDrawThread)
			else:
				self._projId = None
		else:
			self._projId = self.after(TIMER, self.projectAfter)

	# ----------------------------------------------------------------------
	def endProjection(self, event=None):
		#say("*-*", self._viewer.title(), "endProjection")
		if self._projId is not None:
			self.after_cancel(self._projId)
		self.projectAfter()

	# ----------------------------------------------------------------------
	def startDrawThread(self, event=None):
		#say("*-*", self._viewer.title(), "startDrawThread")
		self._viewer.draw(True)
		self._drawId = self.after(TIMER, self.drawAfter)

	# ----------------------------------------------------------------------
	def drawAfter(self):
		#say("*-*", self._viewer.title(), "drawAfter", self._viewer.state())
		self._viewer.expose()
		if self._viewer.state() == PROJECTION_DRAW:
			self._drawId = self.after(TIMER, self.drawAfter)
		else:
			self._drawId = None
			self.check4Errors()

	# ----------------------------------------------------------------------
	def endDraw(self, event=None):
		#say("*-*", self._viewer.title(), "endDraw")
		if self._drawId is not None:
			self.after_cancel(self._drawId)
		self._viewer.expose()
		self._drawId = None
		self.check4Errors()

	# ----------------------------------------------------------------------
	def check4Errors(self):
		if self._viewer.error("kernel"):
			self.flair().notify("Kernel Error",
				self._viewer.error("kernel"),
				tkFlair.NOTIFY_ERROR)

	# ----------------------------------------------------------------------
	# Stop all pending threads
	# ----------------------------------------------------------------------
	def stop(self):
		if self._projId is not None: self.after_cancel(self._projId)
		if self._drawId is not None: self.after_cancel(self._drawId)
		self._projId = None
		self._afterd = None
		self._viewer.stopThread()

	# ----------------------------------------------------------------------
	# Return  snap to grid status or closest handle
	# ----------------------------------------------------------------------
	def snap(self, event):
		snap = self.editor.snap2Grid.get()
		if event.state & SHIFT_MASK:
			if event.state & CONTROL_MASK:
				self._geometry.set("snapdistance",  5*SNAPDISTANCE)
			else:
				self._geometry.set("snapdistance",  SNAPDISTANCE)
				snap = not snap
		return snap

	# ----------------------------------------------------------------------
	# Mouse button1 pressed
	# ----------------------------------------------------------------------
	def button1(self, event):
		self.focus_set()
		self.panStop()
		self.x = self._x = event.x
		self.y = self._y = event.y

		if self.action in (ACTION_MOVE, ACTION_ROTATE):
			self.actionMoveStop(event)
			return

		self._cursor = self.cget("cursor")

		if self.parent.mouseAction in (ACTION_ZOOM_IN, ACTION_ZOOM_OUT):
			self.action = self.parent.mouseAction
			self.config(cursor=mouseCursor(self.action))

		elif self.parent.mouseAction == ACTION_PAN:
			self.actionPan(event)

		elif self.parent.mouseAction == ACTION_ORBIT:
			self.actionOrbit(event, True)

		elif self.parent.mouseAction == ACTION_ADD_ZONE:
			self.actionZone(event.x, event.y)
			self.editor.setMouseAction(ACTION_SELECT)

		elif self.parent.mouseAction == ACTION_ADD:
			self.actionAdd(event)

		elif self.parent.mouseAction == ACTION_ADD_NEXT:
			self._viewer.edit("start", (self._viewer.pixel2uv(event.x, event.y), self.snap(event)))
			self.resetAnchor()

#		elif event.state & SHIFT_MASK:
#			if Input._developer:
#				self.action = ACTION_PEN
#				self._viewer.pen("clear")
#				self._viewer.pen("add",(event.x,event.y))
#				self._viewer.expose()
#			else:
#				self.action = ACTION_AREA_SELECT
#			self.config(cursor=mouseCursor(self.action))

		elif self.parent.mouseAction == ACTION_ZONEPAINT:
			self.action = self.parent.mouseAction
			self.parent.zone("clear")
			self.actionZoneShow(event)

		elif self.parent.mouseAction == ACTION_VIEW_CENTER:
			self.centerViewports(event)

		elif self.parent.mouseAction == ACTION_INFO:
			self.actionInfo(event)

		else:
			self.actionSelect(event)

	# ----------------------------------------------------------------------
	def button2(self, event):
		if self.action == ACTION_ADD_NEXT: return
		self.focus_set()
		self.panStop()
		self._cursor = self.cget("cursor")
		self.x = self._x = event.x
		self.y = self._y = event.y
		if LAPTOP:
			self.buttonMenu(event)
		else:
			self.buttonPan(event)

	# ----------------------------------------------------------------------
	def button3(self, event):
		if self.action == ACTION_ADD_NEXT: return
		self.focus_set()
		self.panStop()
		self.x = self._x = event.x
		self.y = self._y = event.y

		if self.action == ACTION_PAN:	# button-2 already pressed
			self.actionOrbit(event, False)
			return

		self._cursor = self.cget("cursor")
		if LAPTOP:
			self.buttonPan(event)
		else:
			if event.state & CONTROL_MASK:
				self.actionOrbit(event, False)
			elif event.state & SHIFT_MASK:
				self.actionPan(event)
			else:
				self.buttonMenu(event)

	# ----------------------------------------------------------------------
	def double1(self, event):
		pass

	# ----------------------------------------------------------------------
	# Event: move mouse without any button pressed
	# Display the position of the mouse, and/or perform object
	# transformations like move,rotate,scale
	# ----------------------------------------------------------------------
	def motion(self, event):
		self.parent.setStatus(self._viewer.pixel2xyz(event.x, event.y))
		if self.parent.balloon and \
		   self.parent.mouseAction not in (ACTION_ADD_ZONE, ACTION_INFO):
			self.parent.hideBalloon()

		if self.action in (ACTION_MOVE, ACTION_ADD_NEXT):
			self.actionMove(event)

		elif self.action == ACTION_ROTATE:
			self.actionRotate(event)

		elif self.parent.mouseAction in (ACTION_ADD_ZONE, ACTION_ZONEPAINT):
			self.actionZoneShow(event.x, event.y)

	# ----------------------------------------------------------------------
	# Event: move mouse while pressing any button
	# ----------------------------------------------------------------------
	def buttonMotion(self, event):
		if self.action == ACTION_ZONEPAINT:
			self.actionZoneShow(event)
			return

		if self.snap(event):
			uv  = self._viewer.pixel2uv(event.x, event.y)
			xyz = self._viewer.snap(uv[0], uv[1])
		else:
			xyz = self._viewer.pixel2xyz(event.x, event.y)
		self.parent.setStatus(xyz)

		if self.action == ACTION_SEL_OR_MOVE:
			# Check for the type of action
			if abs(event.x-self._x)>1 or abs(event.y-self._y)>1 and \
			   self._closest[0] in ('B','O','C'):
				self._actionMoveOrRotate(event)

		elif self.action == ACTION_SEL_OR_AREA:
			if abs(event.x-self._x)>1 or abs(event.y-self._y)>1:
				if self.parent.mouseAction == ACTION_EDIT:
					self.action = ACTION_PEN
					self._viewer.pen("add", (event.x, event.y))
					self._viewer.expose()
				else:
					self.action = ACTION_AREA_SELECT
					self.config(cursor=mouseCursor(self.action))
					self._viewer.rectangle(self.x, self.y, event.x, event.y)
					self._viewer.expose()

		elif self.action == ACTION_PEN:
			self._viewer.pen("add", (event.x, event.y))
			self._viewer.expose()

		elif self.action in (ACTION_MOVE, ACTION_ADD, ACTION_ADD_NEXT):
			self.actionMove(event)

		elif self.action == ACTION_ROTATE:
			self.actionRotate(event)

		elif self.action==ACTION_VIEW_MOVE: # or \
		   #(self.action==ACTION_VIEW_CENTER and event.state&CONTROL_MASK):
			matrix  = self._viewer.matrix()
			matrixv = Matrix(self.viewports[self.selected]._viewer.matrix())
			w  = Vector(matrix[0][2], matrix[1][2], matrix[2][2])
			wv = Vector(matrixv[0][2],matrixv[1][2],matrixv[2][2])

			# Move along _axis (xyz-save) * axis
			d = (xyz[0]-self._save[0])*self._axis[0] + \
			    (xyz[1]-self._save[1])*self._axis[1] +\
			    (xyz[2]-self._save[2])*self._axis[2]

			# The center calculated as
			#     center = self._save + d*self._axis
			#        C   = S + d*A
			# is along the projection of the real center
			# Therefore we need to find the intersection of the two
			# lines:
			#    1) perpendicular to the screen passing from 'center'
			#        S' = C + w*t = S + d*A + w*t
			#    2) along the viewport center towards the wv vector
			#        S' = S + wv*tv
			# From 1,2
			#       wv*tv = d*A + w*t
			# Taking the cross product with w
			#       w^wv*tv = d*w^A + t*w^w = d*w^A
			# Therefore we can solve for the tv and projected it
			# along the line 2

			A = w ^ wv
			B = d * (w^self._axis)
			try: tv =B[0]/A[0]
			except: pass
			try: tv = B[1]/A[1]
			except: pass
			try: tv = B[2]/A[2]
			except: pass

			center = self._save + wv*tv

			self._viewer.viewport(self.selected, "center", tuple(center))
			self.drawFast3DOtherViewports()
			self.parent.expose()

		elif self.action == ACTION_VIEW_CENTER:
			# Move center of viewport
			if self.moveAbsolute(event):
				self._viewer.viewport(self.selected, "center", xyz)
			else:
				self._viewer.viewport(self.selected, "center",
					self._viewer.moveUV(self._save, xyz))
			self.drawFast3DOtherViewports()
			self.parent.expose()

		elif self.action == ACTION_VIEW_ROTATE:
			# Rotate viewport
			x,y = self._viewer.viewport(self.selected, "centerview")
			phi = atan2(event.y-y, event.x-x)
			#if event.state & SHIFT_MASK:
			#	phi = round(phi/ROTSTEP)*ROTSTEP

			phi -= self._phi
			m = Matrix(4)
			m.rotate(-phi, self._axis)
			m = m*self._save
			if event.state & SHIFT_MASK:
				rx, ry, rz = m.getEulerRotation()
				rx = round(rx/ROTSTEP)*ROTSTEP
				ry = round(ry/ROTSTEP)*ROTSTEP
				rz = round(rz/ROTSTEP)*ROTSTEP
				m = Matrix.eulerRotation(rx,ry,rz)
			xyz = self.viewports[self.selected]._viewer.origin()
			self.viewports[self.selected]._viewer.matrix(m)
			self.viewports[self.selected]._viewer.origin(xyz)
			self.drawFast3DOtherViewports()
			self.parent.expose()

		elif self.action in (ACTION_AREA_SELECT, ACTION_ZOOM_IN, ACTION_ZOOM_OUT):
			self._viewer.rectangle(self.x, self.y, event.x, event.y)
			self._viewer.expose()

		elif self.action == ACTION_PAN:
			dx = event.x - self.x
			dy = event.y - self.y
			dt = float(event.time - self.t)
			if dt>0.0:
				self.vx = float(dx) / dt
				self.vy = float(dy) / dt

			ox, oy = self._viewer.offset()
			if self._viewer.get("3D") and self._viewer.get("projection"):	# Perspective or Combo
			#if self.frame.getInt("3D") and self.frame.getInt("3D.projection"):
				z = self._viewer.zoom()
			else:
				z = 1.0
			self._viewer.offset(ox-dx/self._viewer.Sx()*z, oy+dy/self._viewer.Sy()*z)
			self.x = event.x
			self.y = event.y
			self.t = event.time
			if self._viewer.get("3D"):
				self._viewer.draw(False, DRAW_FAST)
				self._viewer.expose()
			else:
				self.draw()

		elif self.action == ACTION_ORBIT:
			width  = float(self.winfo_width())
			height = float(self.winfo_height())
			m = min(width, height)

			if self._anchorAxis == 0:
				axis = Vector(	self._save[0][0],
						self._save[0][1],
						self._save[0][2])
			elif self._anchorAxis == 1:
				axis = Vector(	self._save[1][0],
						self._save[1][1],
						self._save[1][2])
			elif self._anchorAxis == 2:
				axis = Vector(	self._save[2][0],
						self._save[2][1],
						self._save[2][2])
			else:
				axis = None

			#self.quat = trackball((2.0*self.x - width)  / m,
			#		  (height - 2.0*self.y) / m,
			#		  (2.0*event.x - width)  / m,
			#		  (height - 2.0*event.y) / m) + self.quat
			axis,phi = trackball2((2.0*self.x - width) / m,
					  (height - 2.0*self.y)    / m,
					  (2.0*event.x - width)    / m,
					  (height - 2.0*event.y)   / m, axis)
			if axis is None: return

			m = Matrix(4)
			m.rotate(phi,axis)
			self._matrix *= m	# XXX could have many errors creep in!!!
			if event.state & SHIFT_MASK:
				#rx,ry,rz = self.quat.matrix().getEulerRotation()
				rx,ry,rz = self._matrix.getEulerRotation()
				# round with ROTSTEP
				rrx = round(rx/ROTSTEP)*ROTSTEP
				rry = round(ry/ROTSTEP)*ROTSTEP
				rrz = round(rz/ROTSTEP)*ROTSTEP
				if self._anchorAxis == 0:
					rry = ry
					rrz = rz
				elif self._anchorAxis == 1:
					rrx = rx
					rrz = rz
				elif self._anchorAxis == 2:
					rrx = rx
					rry = ry
				m = Matrix.eulerRotation(rrx,rry,rrz)
				self._viewer.matrix(m)
			else:
				self._viewer.matrix(self._matrix)
			#	self._viewer.matrix(self.quat.matrix())

			# Restore origin
			if self._viewer.get("3D"):
				# Pan backward
				w  = self._viewer.basis(b'w')
				x  = self._origin[0] + self._orbitDepth*w[0]
				y  = self._origin[1] + self._orbitDepth*w[1]
				z  = self._origin[2] + self._orbitDepth*w[2]
				self._viewer.origin((x,y,z))
			else:
				self._viewer.origin(self._origin)

			self._viewer.draw(False, DRAW_FAST)
			self._viewer.expose()
			self.x = event.x
			self.y = event.y

		elif self.parent.mouseAction == ACTION_INFO:
			self.actionInfo(event)

	# ----------------------------------------------------------------------
	def buttonRelease(self, event):
		self.config(cursor=self._cursor)

		if self.action in (ACTION_SEL_OR_MOVE, ACTION_SEL_OR_AREA):
			if self._closest is not None:
				if event.state&CONTROL_MASK or \
				   self._closest[0] not in ('B','O','C') or \
				   self._closest[2]>0:
					self.actionDoSelection(event)
				else:
					self._actionMoveOrRotate(event)
					return	# Return not to change action

		elif self.action in (ACTION_MOVE, ACTION_ROTATE):
			self.actionMoveStop(event)

		elif self.action in (ACTION_VIEW_CENTER, ACTION_VIEW_MOVE, ACTION_VIEW_ROTATE):
			viewer = self.viewports[self.selected]
			viewer.draw(True,True)
			viewer.expose()
			for v in self.viewports: v.expose()
			self.parent.setStatusMsg()

		elif self.action == ACTION_AREA_SELECT:
			toggle = True
			clear = not bool(event.state&CONTROL_MASK)
			if self.editor.editType() == Input.Card.REGION and \
			   self._geometry.selection("b",BIT_SELECT)==0:
				clear  = bool(event.state&SHIFT_MASK)
				toggle = bool(event.state&CONTROL_MASK)

			objs = self._viewer.enclosed(self.x, self.y,
					event.x, event.y, "B")

			if objs:
				self.editor.select(("B",objs), clear, toggle)
				clear = False

			objs = self._viewer.enclosed(self.x, self.y,
					event.x, event.y, "O")
			if objs:
				self.editor.select(("O",objs), clear, toggle)

			self._viewer.rectangle(-1,-1,-1,-1)
			self.draw(False, True)

		elif self.action == ACTION_PAN:
			# Only click
			if abs(event.x-self._x)<2 and abs(event.y-self._y)<2:
				self._viewer.viewOto0()

				# Center current position
				u,v = self._viewer.pixel2uv(event.x, event.y)
				xyz = Vector(self._viewer.uv2xyz(u,v))
				origin = xyz

				# For 3D displays find hit position to center
				if self._viewer.get("3D"): # and self._viewer.get("projection"):
					# For Perspective or Combo
					# Calculate distance for orbit operations
					hit,regtype = self._viewer.hit3D(u,v,"T")
					if regtype != "BLCKHOLE":
						hit = Vector(hit)
						# Is not correct, it should be calculated after the move but
						# should be a reasonable approximation
						xyz0 = Vector(self._viewer.uv2xyz(0.,0.))
						cpos = Vector(self._viewer.camera("position"))
						dhit = (hit-cpos).length()
						width = float(self.winfo_width())/self._viewer.Sx()
						if self._viewer.get("projection")==1:	# Perspective
							# Find the location of the hit object above/below
							# viewing plane
							#if (hit-xyz) * Vector(self._viewer.basis('W')) > 0.0:
							#	self._orbitDepth = -self._orbitDepth
							self._orbitDepth = 0.0;
							origin = hit
						else:
							self._orbitDepth = (hit-xyz).length()
							if self._viewer.get("projection")==2:	# Combo
								dxyz  = (xyz-cpos).length()
								origin = dhit/dxyz*xyz + (1.0-dhit/dxyz)*xyz0

				# Set new origin
				self._viewer.origin(*origin)
				self.draw(True)

			#-----------------------
			# if we are still moving
			#-----------------------
			elif event.time - self.t < 100:
				self.action = ACTION_PAN_SLOWDOWN
				self._panId = self.after(TIMERPAN, self.panSlowdown)
				return

			else:
				self.draw()

		elif self.action in (ACTION_ZOOM_IN, ACTION_ZOOM_OUT):
			x1,y1 = self._viewer.pixel2uv(self.x, self.y)
			x2,y2 = self._viewer.pixel2uv(event.x, event.y)
			self._viewer.offset(0.5*(x1+x2), 0.5*(y1+y2))

			zoom = self.action==ACTION_ZOOM_IN
			if event.state & CONTROL_MASK: zoom = not zoom
			if abs(event.x-self._x)>1 or (event.y-self._y)>1:
				if zoom:
					w = abs(x2-x1)
					h = abs(y2-y1)
					if w>0 and h>0:
						x1,y1,x2,y2 = self._viewer.extends()
						zx = (x2-x1) / w
						zy = (y2-y1) / h
						self._viewer.zoom(min(zx,zy))
				else:
					z = self._viewer.zoom()
					zx = float(abs(event.x-self._x)) / \
					     float(self.winfo_width())  * z
					zy = float(abs(event.y-self._y)) / \
					     float(self.winfo_height()) * z
					self._viewer.zoom(max(zx,zy))
			else:
				if zoom:
					self.zoomIn(event)
				else:
					self.zoomOut(event)
			self._viewer.rectangle(-1,-1,-1,-1)
			self.draw(False, True)

		elif self.action == ACTION_ORBIT:
			if self.parent.mouseAction != ACTION_ORBIT:
				self._viewer.set("trackball", 0)
			self.draw(True, True)

		elif self.action in (ACTION_ADD, ACTION_ADD_NEXT):
			self.actionAddNext(event)
			return

		self.action = ACTION_SELECT

	# ----------------------------------------------------------------------
	def buttonMenu(self, event):
		if self.action not in (ACTION_SELECT, ACTION_PAN, ACTION_ORBIT): return
		self.popupMenu(event)

	# ----------------------------------------------------------------------
	def buttonPan(self, event):
		if event.state & CONTROL_MASK:
			self.actionOrbit(event, False)
		elif event.state & SHIFT_MASK:
			self.actionZoomRect(event)
		else:
			self.actionPan(event)

	# ----------------------------------------------------------------------
	def panSlowdown(self):
		self.vx *= FRICTION
		self.vy *= FRICTION
		dx = self.vx * TIMERPAN
		dy = self.vy * TIMERPAN

		if self.action== ACTION_PAN_SLOWDOWN and (abs(dx)>=2.0 or abs(dy)>=2.0):
			ox, oy = self._viewer.offset()
			self._viewer.offset(ox-dx/self._viewer.Sx(), oy+dy/self._viewer.Sy())
			self._panId = self.after(TIMERPAN, self.panSlowdown)
			if self._viewer.get("3D"):
				self._viewer.draw(False, DRAW_FAST)
				self._viewer.expose()
			else:
				self.draw()
		else:
			self.action = ACTION_SELECT
			self.historyAdd()
			self._panId = None
			self.draw()

	# ----------------------------------------------------------------------
	def panStop(self):
		if self.action == ACTION_PAN_SLOWDOWN:
			self.action = ACTION_SELECT
			if self._panId: self.after_cancel(self._panId)
			self._panId = None
			self.historyAdd()
			self.draw()

	# ----------------------------------------------------------------------
	def actionStop(self, event=None):
		if self.action == ACTION_ORBIT:
			self._viewer.matrix(self._save)
			self._viewer.set("trackball", 0)

		elif self.action in (ACTION_VIEW_CENTER, ACTION_VIEW_MOVE):
			self._viewer.viewport(self.selected, "center", tuple(self._save))

		elif self.action == ACTION_VIEW_ROTATE:
			self.viewports[self.selected]._viewer.matrix(self._save)

		elif self.action in (ACTION_MOVE, ACTION_ROTATE):
			self.actionMoveCancel()

		elif self.action == ACTION_ADD_NEXT:
			self.actionMoveCancel()
			self.parent.loadProject()
			self.parent.redraw()

		else:
			self._viewer.rectangle(-1,-1,-1,-1)

		self.editor.setMouseAction(ACTION_SELECT)
		self.action = ACTION_SELECT
		self._viewer.expose()

	# ----------------------------------------------------------------------
	# Default mouse action depending on the proximity of the mouse to
	# objects as well the key-combinations Ctrl/Shift
	# ----------------------------------------------------------------------
	def actionSelect(self, event):
		if self.parent.mouseAction < 0:
			opt = name = None
			if self.parent.mouseAction == ACTION_SELECT_BODIES:
				opt = "B"
			elif self.parent.mouseAction == ACTION_SELECT_REGION:
				opt = "R"
			elif self.parent.mouseAction == ACTION_SELECT_ZONE:
				if self.editor.editType() == Input.Card.REGION:
					opt = "Z"
					name = self.editor.editObject.name()
				else:
					opt = "RZ"
			self._closest = self._viewer.closest(event.x, event.y, 4.0, opt, name)
		elif self.editor.editType() == Input.Card.REGION:
			self._closest = self._viewer.closest(event.x, event.y, 4.0,
						None,
						self.editor.editObject.name())
		else:
			self._closest = self._viewer.closest(event.x, event.y, 4.0)

		if Input._developer: say("GeometryViewer.actionSelect: self._closest=", self._closest)

		if self._closest is None:
			self.action = ACTION_SEL_OR_AREA
			self.parent.setStatusMsg()
			return

		obj,idx,opt = self._closest

		if obj=="V" and idx>=0:
			self.viewports[idx]._viewer.viewOto0()
			self.selected = idx

			x,y = self._viewer.viewport(idx, "centerpixel")

			width  = self.winfo_width()
			height = self.winfo_height()
			ww = width  / 20
			hh = height / 20

			# default action
			self.action = ACTION_VIEW_MOVE

			# check for rotation, if we click far from center
			# and close to borders
			if (event.x-x)**2 + (event.y-y)**2 > ww*hh:
				if (event.x<=ww and x>ww) or (event.x>=width-ww and x<width-ww) or \
				   (event.y<=hh and y>hh) or (event.y>=height-hh and y<height-hh):
					self.action = ACTION_VIEW_ROTATE

			elif abs(event.x-x) < 5 and abs(event.y-y) < 5:
				self.action = ACTION_VIEW_CENTER

			# remember the position of the viewport
			self._save = Vector(self._viewer.viewport(idx, "center"))
			matrix  = self._viewer.matrix()
			matrixv = Matrix(self.viewports[idx]._viewer.matrix())

			w  = Vector(matrix[0][2], matrix[1][2], matrix[2][2])
			wv = Vector(matrixv[0][2],matrixv[1][2],matrixv[2][2])

			if self.action == ACTION_VIEW_CENTER:
				self._axis = wv ^ w

			elif self.action == ACTION_VIEW_MOVE:
				# Tricky we have to move along the wv
				# however the viewport line is NOT the intersection
				# of the two planes rather a vertical line at the CENTER position!
				ww = w ^ wv
				ww.norm()
				self._axis = ww ^ w

			else:
				# remember the matrix of the viewport
				self._save = matrixv
				self._phi  = atan2(event.y-y, event.x-x)
				self._axis = Vector(matrix[0][2],matrix[1][2],matrix[2][2])

			self.config(cursor=mouseCursor(self.action))

		# Initial hit on Zone or Region?
		elif obj in ("R","Z"):
			if self.parent.mouseAction == ACTION_EDIT:
				self._viewer.pen("clear")
				self._viewer.pen("add",(event.x,event.y))
				self._viewer.expose()
			# Postpone the selection after checking if
			# the cursor moved in the mean time
			self.action = ACTION_SEL_OR_AREA

		# Cursor
		elif obj == "C":
			self.action = ACTION_SEL_OR_MOVE

		# Check if we hit on a selected item
		elif self.editor.isSelected((obj,idx)):
			if obj=="B" and self._geometry.cursor("body") != idx:
				self._geometry.cursor("body", idx)
				self._closest = self._viewer.closest(
					event.x, event.y, 4.0, "B")
			# Postpone the selection after checking if
			# the cursor moved in the mean time
			self.action = ACTION_SEL_OR_MOVE

		else:
			self.actionDoSelection(event)

	# ----------------------------------------------------------------------
	# Perform the actual selection
	# ----------------------------------------------------------------------
	def actionDoSelection(self, event):
		t, idx, opt = self._closest
		if t=="R":
			msg = "Region: "
			for i in idx:
				name = self._geometry.region(i,"name")
				if name == "VOXEL":
					if self._viewer.get("3D"):
						uv = self._viewer.pixel2uv(event.x, event.y)
						xyz = self._viewer.hit3D(*uv)
					else:
						xyz = self._viewer.pixel2xyz(event.x, event.y)

					voxelReg = self._geometry.voxel("get",xyz)+1
					idx  = -voxelReg
					msg += "VOXEL%03d "%(voxelReg)
				else:
					msg += "%s "%(name.decode())
			self.parent.setStatusMsg(msg)

		elif t=="Z":
			self.parent.setStatusMsg("Zones: %s"%\
				(" ".join(["zone%02d"%(i+1) for i in idx])))

		elif t=="B":
			self.parent.setStatusMsg("Body: %s [%s]" % \
				(self._geometry.body(idx, "name").decode(),
				 self._geometry.body(idx, "type").decode()))

		elif t=="O":
			self.parent.setStatusMsg("Object: %s" % \
				(self._geometry.object(idx, "name")))
			self.parent.expose()

		# Do the selection
		if self.parent.mouseAction < 0:
			clear = False
		else:
			clear = not bool(event.state&CONTROL_MASK)
		self.editor.select((t,idx), clear, True)

		if self.parent.mouseAction == ACTION_EDIT:
			self.actionOverlaps()

	# ----------------------------------------------------------------------
	# Correct overlaps in 2D
	# ----------------------------------------------------------------------
	def actionOverlaps(self):
		if self.editor.editType() != Input.Card.REGION: return
		region = self.editor.editObject.cards[0]
		zones  = self.editor.propListSelectedZones()
		modified = self._viewer.region(region[ID],"correct", zones)
		if modified:
			expr = self._geometry.region(region[ID],"expr")
			zones = []
			for z in expr:
				zones.append(list2Zone(z[1]))
			self.flair().addUndo(self.flair().setZonesUndo(region,zones))

	# ----------------------------------------------------------------------
	# Display object information at location
	# ----------------------------------------------------------------------
	def actionInfo(self, event):
		closest = self._viewer.closest(event.x, event.y, 4.0, "BR")
		if self.frame.getInt("3D"):
			uv = self._viewer.pixel2uv(event.x, event.y)
			xyz = self._viewer.hit3D(*uv)
		else:
			xyz = self._viewer.pixel2xyz(event.x, event.y)
		message = "Location:\n  x: %s\n  y: %s\n  z: %s"%(tuple(map(format,xyz)))

		# Specific information
		if closest is None:
			pass

		elif closest[0] == "B":
			# Body information
			idx  = closest[1]
			message += "\nBody: %s [%s]"%(
				self._geometry.body(idx,"name"),
				self._geometry.body(idx,"type"))
			what = self._geometry.body(idx,"what")
			message += "\nWhat: %s"%(" ".join(map(format,what)))
			card = self.editor.findCard("B",idx)

			dx = card["@dx"]
			dy = card["@dy"]
			dz = card["@dz"]
			exp = card["@expansion"]
			tra = card["@transform"]
			if dx or dy or dz:
				message += "\nTranslat: %s, %s, %s"%(str(dx), str(dy), str(dz))
			if exp:
				message += "\nExpansion: %s"%(str(exp))
			if tra:
				message += "\nTransform: %s"%(str(tra))

		elif closest[0] == "R":
			# Region information
			idx = closest[1][0]

			name = self._geometry.region(idx,"name")
			name2 = None

			if name == "VOXEL" and self.frame.getInt("Show.voxel"):
				voxelReg = self._geometry.voxel("get",xyz)+1
				idx = -voxelReg
				name2 = "VOXEL%03d"%(voxelReg)

			elif self._geometry.region(idx,"hasmatrix")!=0 and self.frame.getInt("Show.lattice"):
				idx = self._viewer.where(*xyz)
				name2 = self._geometry.region(idx,"name")

			if name2:
				message += "\nRegion: %s [%s]"%(name2.encode(), name.encode())
			else:
#				message += "\nRegion: %s"%(name.encode())
				message += "\nRegion: %s"%(name.decode())

			card = self.editor.findCard("R",idx)
			if card:
				mat = card["@material"]
				if mat: message += "\nMaterial: %s"%(mat.sdum())
				mat2 = card["@matDecay"]
				if mat2 and mat2 is not mat:
					message += "\nmatDecay: %s"%(mat2.sdum())

				#if name2:
				#	message += "\nValue: %g"%(self._viewer.voxel(card[ID],"value"))
				#else:
				#	message += "\nValue: %g"%(self._viewer.region(card[ID],"value"))

				rot = card["@rotdefi"]
				if rot: message += "\nrotdefi: %s"%(rot)

				for prop in card.prop:
					if prop in (	# skip the following
						"@assignmat",
						"@color",
						"@cont",
						"@id",
						"@lattice",
						"@material",
						"@materialDefined",
						"@matDecay",
						"@n",
						"@rotdefi",
						"@select",
						"@type",
						"@zone"): continue
					value = card[prop]
					if prop[0] == "@": prop = prop[1:]
					if isinstance(value,Input.Card):
						if prop in ("biasingA","biasingH","biasingE","biasingN"):
							message += "\n%s"%(prop)
							message += "\n   Imp: %g"%(value.numWhat(3))
							message += "\n   Split: %g"%(value.numWhat(2))
						else:
							message += "\n%s: %s"%(prop,value.sdum())
					elif value:
						message += "\n%s: %s"%(prop,str(value))

				# Check for magnetic field
				if "@assignmat" in card.prop:
					field = card["@assignmat"].intWhat(5)
					if field > 0:
						message += "\nfield: %s"%(_ASSIGNMA_field[field][0])

			# Show information on voxel
			if name == "VOXEL" and self.frame.getInt("Show.voxel"):
				for card in self.input()["VOXELS"]:
					if card.ignore() or card["@voxel"] is None: continue
					voxel = card["@voxel"]
					message += "\nVoxel:"
					message += "\n   Data: %s = %d"%(
								str(self._geometry.voxel("index",xyz)),
								self._geometry.voxel("get",xyz))
					message += "\n   File: %s"%(voxel.filename)
					message += "\n   Title: %s"%(voxel.title)
					message += "\n   Nxyz: %d x %d x %d"%(voxel.nx, voxel.ny, voxel.nz)
					message += "\n   dxyz: %g x %g x %g"%(voxel.dx, voxel.dy, voxel.dz)
					message += "\n   No: %d"%(voxel.no)
					message += "\n   Mo: %d"%(voxel.mo)
					message += "\n   Rtstruct: %d"%(voxel.roiN)

			# Any usrbin visible?
			for i,tag in enumerate(GeometryLayers.USRBINNAMES):
				if self.frame.getInt(tag):
					value = self._viewer.usrbin(i,"value",xyz)
					if value is not None:
						message += "\n%s: %s"%(tag, format(value))

		self.clipboard_clear()
		self.clipboard_append(message)
		self.parent.setBalloon(message, 60000)

	# ----------------------------------------------------------------------
	def actionOrbit(self, event, lockAxis):
		if self.action in (ACTION_MOVE, ACTION_ROTATE, ACTION_ADD_NEXT): return

		self.action = ACTION_ORBIT
		self.config(cursor=mouseCursor(self.action))
		self._viewer.viewOto0()

		# Save matrix before any transformation
		self._save   = self._viewer.matrix()
		self._matrix = Matrix(self._save)
		self._origin = Vector(self._viewer.origin())

		width  = float(self.winfo_width())
		height = float(self.winfo_height())
		m = min(width,height)
		x = (2.0*event.x - width)  / m
		y = (height - 2.0*event.y) / m
		v = Vector(x, y, projection2Sphere(TRACKBALLSIZE, x, y))

		if lockAxis:
			THRESHOLD = 0.05
			self._anchorValue = ""

			vv = list(map(abs, Matrix(self._viewer.matrix()).multNoTranslation(v)))
			if vv[0] <= vv[1] and vv[0] <= vv[2] and abs(vv[0])<THRESHOLD:
				self._anchorAxis = 0
			elif vv[1] <= vv[0] and vv[1] <= vv[2] and abs(vv[1])<THRESHOLD:
				self._anchorAxis = 1
			elif vv[2] <= vv[0] and vv[2] <= vv[1] and abs(vv[2])<THRESHOLD:
				self._anchorAxis = 2
			else:
				self._anchorAxis = -1
		else:
			self.resetAnchor()

		if self.frame.getInt("3D"):
			# Pan backward
			w  = Vector(self._viewer.basis(b'W'))

			# Window width
			width = float(self.winfo_width())/self._viewer.Sx()
			if self._viewer.get("projection") == 1:
				# Perspective
				minw = -self._viewer.camera("focal")
			else:
				# Orthographic and Combo
				minw = 0.0

			# Limit orbit distance to -minw
			if self._orbitDepth < minw:
				self._orbitDepth = minw/2.0

			self._origin -= self._orbitDepth*w

		#self.quat   = Quaternion(Matrix(self._viewer.matrix()))
		self._viewer.set("trackball", 1)

		self._viewer.expose()

	# ----------------------------------------------------------------------
	def actionZoomRect(self, event):
		self.action = ACTION_ZOOM_IN
		self.config(cursor=mouseCursor(self.action))

	# ----------------------------------------------------------------------
	def actionPan(self, event):
		self.action = ACTION_PAN
		self.config(cursor=mouseCursor(self.action))
		self.vx = 0.0
		self.vy = 0.0
		self.t  = event.time

	# ----------------------------------------------------------------------
	# Add a new object/body
	# ----------------------------------------------------------------------
	def actionAdd(self, event):
		if self.action in (ACTION_MOVE, ACTION_ROTATE, ACTION_ADD_NEXT): return

		u,v = self._viewer.pixel2uv(event.x, event.y)
		if self.snap(event):
			xyz = self._viewer.snap(u, v)
		else:
			xyz = self._viewer.pixel2xyz(event.x, event.y)

		self._viewer.edit("start", ((u,v), self.snap(event)))
		self.parent.stopThreads()

		grid_size = self._viewer.grid("size")[0]
		matrix    = self._viewer.matrix()
		self.__opt = 0
		self.__x   = self.x
		self.__y   = self.y

		# Adding flair object
		if self.parent.actionOption[0] == "!":
			tag = self.parent.actionOption
			what = [self.editor.findUniqueName(Input.Card.OBJECT,tag[1:].lower())]
			what.extend(list(xyz))
			if tag not in ("!point", "!mesh"):
				what.append(0.)		# Options
				if tag == "!ruler":
					what.append(1.)	# simple ruler (double head arrow)
				else:
					what.append(0.)
				what.append(0.)

				if tag == "!ruler":
					what.extend(list(xyz))	# other end
					what.extend(list(xyz))	# angle
				else:
					what.append(0.)		# vector length
					what.append(0.)
					what.append(0.)

				if tag == "!light":
					what.append(1.0)	# power
					what.append(0.0)	# falloff
					what.append(0.0)	# specular

			card = Input.Card(tag, what)
			card[ID] = self.parent.addObject(card)
			self._select([card],None)
			self._geometry.object(self._objectsid, "save", True)

			self._msg = "Add object %s"%(tag)

		else:	# Adding body
			card = GeometryBody.newBody(
					self.parent.actionOption,
					xyz,
					grid_size,
					matrix)
			card.setSdum(self.editor.findUniqueName(Input.Card.BODY))
			whats = [card.numWhat(i) for i in range(1,card.nwhats()) ]
			card[ID] = self.parent.addBody(card.sdum(), card.tag, whats)

			# FIXME a bit overshoot!
			self.parent.derive()

			self.parent.viewportProject(False, False)
			self.parent.draw()
			self.parent.expose()	# all viewports

			self._select(None,[card])
			self._geometry.cursor("body", self._bodiesid[0])
			self._geometry.body(self._bodiesid, "save", True)
			self._geometry.body(self._bodiesid, "show", BIT_SELECT|BIT_WIREFRAME|BIT_MOVE)

			self._msg = "Add body %s"%(card.tag)

			# Special treatment for the add sequence
			matrix = self._viewer.matrix()
			w = Vector(matrix[0][2], matrix[1][2], matrix[2][2])
			if card.tag == "RPP":
				if   abs(abs(w[0])-1.0) < 1.0E-5:
					_ADD_SEQUENCE[card.tag] = [16]
				elif abs(abs(w[1])-1.0) < 1.0E-5:
					_ADD_SEQUENCE[card.tag] = [17]
				elif abs(abs(w[2])-1.0) < 1.0E-5:
					_ADD_SEQUENCE[card.tag] = [18]
				else:
					_ADD_SEQUENCE[card.tag] = [18, 6]

			elif card.tag == "XEC":
				if   abs(abs(w[1])-1.0) < 1.0E-5:
					_ADD_SEQUENCE[card.tag] = [2]
				elif abs(abs(w[2])-1.0) < 1.0E-5:
					_ADD_SEQUENCE[card.tag] = [1]
				else:
					_ADD_SEQUENCE[card.tag] = [1, 2]

			elif card.tag == "YEC":
				if   abs(abs(w[0])-1.0) < 1.0E-5:
					_ADD_SEQUENCE[card.tag] = [1]
				elif abs(abs(w[2])-1.0) < 1.0E-5:
					_ADD_SEQUENCE[card.tag] = [2]
				else:
					_ADD_SEQUENCE[card.tag] = [1, 2]

			elif card.tag == "ZEC":
				if   abs(abs(w[0])-1.0) < 1.0E-5:
					_ADD_SEQUENCE[card.tag] = [2]
				elif abs(abs(w[1])-1.0) < 1.0E-5:
					_ADD_SEQUENCE[card.tag] = [1]
				else:
					_ADD_SEQUENCE[card.tag] = [1, 2]
		self.grab_set()
		self.actionAddNext(event)

	# ----------------------------------------------------------------------
	def actionAddNext(self, event=None):
		self.editor.setMouseAction(ACTION_ADD_NEXT)
		self.action = ACTION_ADD_NEXT

#		self._viewer.edit("start", (self._viewer.pixel2uv(event.x, event.y), self.snap(event)))
		self.resetAnchor()

		if self.parent.actionOption[0] == "!":
			tag = self._objects[0].tag
			self._geometry.object(self._objectsid, "save", True)
			#self._viewer.edt("object", self._objectsid[0])
		else:
			tag = self._bodies[0].tag
			self._geometry.body(self._bodiesid, "save", True)
			self._viewer.edit("body", self._bodiesid[0])

		seq = _ADD_SEQUENCE.get(tag)
		if seq is not None:
			try:	# Find current item
				idx = seq.index(self.__opt)
			except ValueError:	# On error start with first one
				self.__opt = seq[0]
				return
			else:
				# If button1 released at the same position
				# set parameter without click-dragging, but simply dragging
				if idx==0 and abs(event.x-self.__x)<=1 and abs(event.y-self.__y)<=1:
					return
				else:
					try:	# Get the next item
						self.__opt = seq[idx+1]
						return
					except:	# On end stop
						pass
		self.actionMoveStop()
		if self._bodies:
			self.editor.addBody(self._bodies[0])
			self._geometry.cursor("body", self._bodiesid[0])

			self.parent.viewportProject(False, False)
			self.parent.draw()
			self.parent.expose()	# all viewports
		else:
			self.editor.addObject(self._objects[0])
		self.editor.setMouseAction(ACTION_SELECT)

	# ----------------------------------------------------------------------
	# Decide whether to rotate or move next
	# ----------------------------------------------------------------------
	def _actionMoveOrRotate(self, event):
		# Clone if Control is pressed
		self._clone = bool(event.state&SHIFT_MASK)
		if self._clone:
			self.editor.clone()
			self.editor.forceRefresh()

		self.resetAnchor()
		if self._closest[2] == 0:
			self.actionMoveStart(event)
		elif -4<self._closest[2]<0:
			self.resetAnchor(-self._closest[2]-1)
			self.actionMoveStart(event)
		elif self._closest[2] == -4:
			self.actionRotateStart(event)
		else:
			self.actionMoveStart(event, self._closest[2])

		# Correct the cursor after cloning
		if self._clone:
			if self._bodiesid:
				self._geometry.cursor("body", self._bodiesid[0])
			self._geometry.cursor("cursor", False)

	# ----------------------------------------------------------------------
	# Start moving action on objects or bodies
	# ----------------------------------------------------------------------
	def actionMoveStart(self, event, opt=0):
		if self.action in (ACTION_ADD_NEXT, ACTION_MOVE, ACTION_ROTATE): return
		if self.editor.lockEditing.get(): return
		if not self._select(*self.editor.getMoveableSelection()): return

		self.parent.stopThreads()

		self.action = ACTION_MOVE
		self.parent.setStatusMsg("Move mouse or Type coordinates e.g. 5 10 20 or x50 /5 or x10 *5",
				"Blue","White")
		self._cursor = self.cget("cursor")
		if opt>=20:
			self.config(cursor=mouseCursor(ACTION_ROTATE))
		else:
			self.config(cursor=mouseCursor(self.action))

		self.x = event.x
		self.y = event.y
		if self._closest is None:
			self._viewer.edit("start",
				(self._viewer.pixel2uv(event.x, event.y),
				 self.snap(event)))
		elif self._closest[0] == 'B':
			self._viewer.edit("body", self._closest[1])
		elif self._closest[0] == 'O':
			self._viewer.edit("object", self._closest[1], self._closest[2])
		elif self._closest[0] == 'C':
			self._viewer.edit("cursor", None)
		else:
			self._viewer.edit("start",
				(self._viewer.pixel2uv(event.x, event.y),
				 self.snap(event)))

		if self._objectsid:
			self._geometry.object(self._objectsid, "save", True)
		if self._bodiesid:
			self._geometry.body(self._bodiesid, "save", True)
			self._geometry.body(self._bodiesid, "show",
					BIT_SELECT|BIT_WIREFRAME|BIT_MOVE)
		# XXX It should be saved for each object!!!
		self.__opt = opt
		if opt<=0:
			self._msg  = "Move selected objects"
		else:
			self._msg  = "Resize selected objects"

		self._geometry.cursor("cursor", False)
		self.grab_set()

	# ----------------------------------------------------------------------
	# FIXME merge _anchorAxis with _axis
	# ----------------------------------------------------------------------
	def actionMove(self, event):
		if self.editor.lockEditing.get(): return
		if not self._anchorValue:
			# Calculate move distance
			if self._closest and self._closest[0] == 'B':
				# Bodies
				# Position absolute everything else relative
				if self.__opt>0:
					relative = False
				else:
					relative = not bool(event.state & CONTROL_MASK)
			else:
				# Objects
				# Everything absolute unless CONTROL is pressed
				relative = not bool(event.state & CONTROL_MASK)
			self._viewer.edit("move", (
					self._viewer.pixel2uv(event.x, event.y),	# u,v
					self._anchorAxis,				# axis
					self.snap(event),				# snap
					relative))					# relative
#					not bool(event.state & CONTROL_MASK)))		# relative
			move = self._viewer.edit("move")
			self.parent.setStatusMsg("Move: dx=%g  dy=%g dz=%g"%tuple(move), timeout=0)
		else:
			# Parse anchorValue
			val = []
			for x in self._anchorValue.split():
				try: val.append(float(x))
				except ValueError: pass
			if val: self._viewer.edit("value", val, self._anchorAxis)

		# Move objects
		self._actionDoMove()

		# Redraw if needed
		if self._objectLight: self._geometry.setLights()
		if self._bodiesid:
			self.parent.viewportProject(False,False)
			self.parent.draw()
		else:
			self.parent.drawFast3D()
		self.parent.expose()	# all viewports

	# ----------------------------------------------------------------------
	# Move selected objects in the dx.dy direction
	# ----------------------------------------------------------------------
	def actionMoveBy(self, event, dx, dy):
		# Find cursor position
		cursor = self._geometry.cursor("cursor")
		u,v = self._viewer.xyz2uv(*cursor)
		event.x, event.y = self._viewer.uv2pixel(u,v)
		self.actionMoveStart(event)

		gx, gy = self._viewer.grid("size")
		u += gx*dx/20.
		v += gy*dy/20.
		self._viewer.edit("move", ((u,v), -1, False, True))
		move = self._viewer.edit("move")

		# Move objects
		self._actionDoMove()

		# Redraw if needed
		if self._bodiesid:
			self.parent.viewportProject(False,False)
			self.parent.draw()
		self.drawFast3DOtherViewports()
		self.parent.expose()	# all viewports

		self.actionMoveStop(event)

	# ----------------------------------------------------------------------
	def _actionDoMove(self):
		self.parent.stopThreads()
		# Move objects
		if self._objectsid:
			if self._closest is not None and self._closest[0] == 'O':
				opt = self.__opt
			else:
				opt = 0
			self._viewer.object(self._objectsid, "move", self.__opt)

		# Move bodies
		if self._bodiesid:
			if self._closest is not None and self._closest[0] == 'B':
				opt = self.__opt
			else:
				opt = 0
			self._viewer.body(self._bodiesid, "move", self.__opt)

	# ----------------------------------------------------------------------
	# Start rotating action on objects or bodies
	# ----------------------------------------------------------------------
	def actionRotateStart(self, event):
		if self.action in (ACTION_ADD_NEXT, ACTION_MOVE, ACTION_ROTATE): return
		if self.editor.lockEditing.get(): return
		if not self._select(*self.editor.getMoveableSelection()): return

		cursor = self._geometry.cursor("cursor")
		if cursor is None: return

		self.parent.stopThreads()

		self.action  = ACTION_ROTATE
		self.parent.setStatusMsg("Move mouse or Type coordinates e.g. 180", "Blue","White")
		self._cursor = self.cget("cursor")
		self.config(cursor=mouseCursor(self.action))

		self._viewer.edit("start", (self._viewer.pixel2uv(event.x, event.y), self.snap(event)))

		if self._bodiesid:
			self._geometry.body(self._bodiesid, "save", True)
			self._geometry.body(self._bodiesid, "show", BIT_SELECT|BIT_WIREFRAME|BIT_MOVE)

			self._viewer.draw(False)
			self.parent.expose()

		# Set pivot point
		self._viewer.edit("pivot", self._viewer.xyz2uv(*cursor))

		matrix = self._viewer.matrix()
		self._axis = Vector(matrix[0][2], matrix[1][2], matrix[2][2])

		self._msg = "Rotate selected items"
		self.grab_set()

	# ----------------------------------------------------------------------
	def actionRotate(self, event):
		if self.editor.lockEditing.get(): return
		if not self._anchorValue:
			self._viewer.edit("rotate", (
					self._viewer.pixel2uv(event.x, event.y),
					self._axis,
					self.snap(event)))
			self.parent.setStatusMsg("Angle: %g"%(self._viewer.edit("rotate")),timeout=0)
		else:
			try: val = float(self._anchorValue)
			except: return
			self._viewer.edit("value", val, self._axis)

		if self._objectsid:
			self._viewer.object(self._objectsid, "rotate")
		if self._bodiesid:
			self._viewer.body(self._bodiesid, "rotate")

		if self._bodiesid:
			self._viewer.project(False, False)
			self._viewer.draw(False)
		self.parent.expose()	# all viewports

	# ----------------------------------------------------------------------
	def actionMoveStop(self, event=None):
		if self.editor.lockEditing.get(): return
		self.grab_release()

		undoinfo = []

		refresh = False
		if self._bodies:
			refresh = True
		else:
			# If a 3D viewport is enabled!
			for view in self.viewports:
				if view._viewer.get("3D"):
					# and there is a light selected force a refresh
					if self._objectLight:
						refresh = True
						break
					break

		if refresh:
			undoinfo.append(self.flair().refreshUndo())

		if self._clone:
			for x in self._anchorValue.split():
				if "*" in x:
					try:
						x = int(x.replace("*",""))
						if self.action == ACTION_MOVE:
							self._actionMoveClone(x, undoinfo)
						elif self.action == ACTION_ROTATE:
							self._actionRotateClone(x, undoinfo)
					except ValueError:
						pass
					break
				elif "/" in x:
					try:
						x = -int(x.replace("/",""))
						if self.action == ACTION_MOVE:
							self._actionMoveClone(x, undoinfo)
						elif self.action == ACTION_ROTATE:
							self._actionRotateClone(x, undoinfo)
					except ValueError:
						pass
					break
			else:
				self.actionMoveUpdate(undoinfo)
		else:
			self.actionMoveUpdate(undoinfo)

		if undoinfo:
			if refresh:
				undoinfo.append(self.flair().refreshUndo())
			else:
				self.editor.updateProperties()
			self.flair().addUndo(undo.createListUndo(undoinfo,self._msg))
		self.actionMoveFinalize()

	# ----------------------------------------------------------------------
	# Create clones of the moved objects
	# Multiply length if times > 0
	# Divide   length if times < 0
	# ----------------------------------------------------------------------
	def _actionMoveClone(self, times, undoinfo):
		self.parent.stopThreads()

		mx,my,mz = self._viewer.edit("move")

		if self._closest[2] > 0:
			s = self._viewer.uv2xyz(*self._viewer.edit("start"))
			#say("s=",s)
#			val = []
#			for x in self._anchorValue.split():
#				try: val.append(float(x))
#				except ValueError: pass
			#if not self._anchorValue:
			uv  = self._viewer.pixel2uv(self.x, self.y)
			xyz = self._viewer.snap(uv[0], uv[1])
			#say("xyz=",xyz)
			sx  = xyz[0] - s[0]
			sy  = xyz[1] - s[1]
			sz  = xyz[2] - s[2]
			mx -= sx
			my -= sy
			mz -= sz
			#say("sxyz=",sx,sy,sz)
			#say("mxyz=",mx,my,mz)
		else:
			sx = sy = sz = 0.0

		if times > 0:
			# Multiply
			start = 0	# Skip first copy (is already there)
			stop  = times	# Do the rest
					# step is equal to move
		else:
			times = -times
			start = 0	# start from the first copy
			stop  = times	# until one before the last
			mx   /= float(times)	# divide step
			my   /= float(times)
			mz   /= float(times)

		for i in range(start, stop):
			dx = sx + mx * (i+1)
			dy = sy + my * (i+1)
			dz = sz + mz * (i+1)
			#say("dxyz=",dx,dy,dz)

			pasted = []
			if i>0:
				# paste the original bodies still in buffer
				self.editor.selectionClear(False)
				cards = self.flair().clipboard2Cards()
				pasted, renamed = self.editor._cardsPaste(cards,undoinfo)
				self._select(*self.editor.getMoveableSelection())
				if self._objectsid:
					self._geometry.object(self._objectsid, "save", True)
				if self._bodiesid:
					self._geometry.body(self._bodiesid, "save", True)
			self._viewer.edit("value",[dx,dy,dz])
			self._actionDoMove()
			self.actionMoveUpdate(undoinfo)

	# ----------------------------------------------------------------------
	# Update the variables back to the object
	# ----------------------------------------------------------------------
	def actionMoveUpdate(self, undoinfo=[]):
		cursorType = self._geometry.cursor("type")
		cursorId   = self._geometry.cursor("id")

		for card in self._objects:
			if card.tag == "ROT-DEFI":
				x,y,z = self._geometry.object(card[ID], "pos")
				card["x"] = x
				card["y"] = y
				card["z"] = z

			else:
				if self.action != ACTION_ADD_NEXT:
					undoinfo.append(self.flair().saveCardUndoInfo(card.pos()))
				for var,w in _OBJECT_VARS[card.tag]:
					card.setWhat(w, format(self._geometry.object(card[ID], var),15))

		for card in self._bodies:
			tag = self._geometry.body(card[ID], "type")
			tag = tag.decode()
			if self.action != ACTION_ADD_NEXT:
				undoinfo.append(self.flair().saveCardUndoInfo(card.pos()))
				if card.tag != "VOXELS" and card.tag != tag:
					if Input._useBOX:
						card.input.changeTag(card,tag)
					elif tag not in ("BOX","WED","RAW","ARB"):
						card.input.changeTag(card,tag)

			what = [format(x,15) for x in self._geometry.body(card[ID], "what")]
			if card.tag == "VOXELS":
				# Get it from the RPP
				if tag == "RPP":
					card.setWhat(1, what[0])
					card.setWhat(2, what[2])
					card.setWhat(3, what[4])
				else:	# BOX	we have to assign a rotation
					card.setWhat(1, what[0])
					card.setWhat(2, what[1])
					card.setWhat(3, what[2])
			else:
				what.insert(0,card.sdum())
				if not Input._useBOX:
					if tag not in ("BOX","WED","RAW","ARB"):
						card.setWhats(what)
				else:
					card.setWhats(what)

		self.flair().setInputModified()
		if self.action != ACTION_ADD_NEXT:
			self._msg = ""

		# FIXME should check for lights!
		if self._objectLight: self._geometry.setLights()
		#self.editor.updateProperties()

		# Restore edit body
		try: cursor = self._geometry.cursor(cursorType, cursorId)
		except: pass

	# ----------------------------------------------------------------------
	# Create clones of the rotated objects
	# Multiply angle if times > 0
	# Divide   angle if times < 0
	# ----------------------------------------------------------------------
	def _actionRotateClone(self, times, undoinfo):
		pass

	# ----------------------------------------------------------------------
	def actionMoveFinalize(self):
		self.resetAnchor()
		self.config(cursor=self._cursor)
		self.action = ACTION_SELECT
		self._clone = False

	# ----------------------------------------------------------------------
	def actionMoveCancel(self, event=None):
		self.grab_release()
		if self._clone:
			# Call undo for the last cloning
			say("FIXME actionMoveCancel: mustUndo")
		self._clone = False
		self._geometry.cursor("cursor", True)
		if self._objectsid:
			self._geometry.object(self._objectsid, "restore")
		if self._bodiesid:
			self._geometry.body(self._bodiesid, "restore")
		self.actionMoveFinalize()
		self._viewer.project(False, False)
		self._viewer.draw(False)
		self.parent.expose()	# all viewports

	# ----------------------------------------------------------------------
	# Reset anchor axis and value
	# ----------------------------------------------------------------------
	def resetAnchor(self, axis=-1):
		self._anchorAxis  = axis
		self._anchorValue = ""

	# ----------------------------------------------------------------------
	# Set and toggle anchor axis
	# ----------------------------------------------------------------------
	def setAnchorAxis(self, axis):
		a = ord(axis) - ord('x')
		if self._anchorAxis == a:
			self._anchorAxis  = -1
			self.parent.setStatusMsg("Clear axis")
		else:
			self.parent.setStatusMsg("Anchor on %r axis"%(axis))
			self._anchorAxis = a
		self._anchorValue = ""

	# ----------------------------------------------------------------------
	# Return True is move has to be absolute on 2D plane or relative on 3D
	# ----------------------------------------------------------------------
	def moveAbsolute(self, event):
		if event.state & CONTROL_MASK:
			moveAbs = False
		else:
			moveAbs = True
		if self.parent.move3D:
			moveAbs = 1 - moveAbs
		return moveAbs

	# ----------------------------------------------------------------------
	# Show zone as the mouse moves
	# ----------------------------------------------------------------------
	def actionZoneShow(self, x, y, clear=True):
		u,v = self._viewer.pixel2uv(x, y)
		zone = self._viewer.zone("find",(u,v,BIT_SELFREEZE))
		if zone and not self._viewer.zone("has",zone):
			#self.parent.stopThreads()
			if clear: self.parent.zone("clear")
			self.parent.zone("select", zone)
			self.parent.draw()
			self.parent.expose()
			self.parent.setStatusMsg(list2Zone(zone[1]),timeout=10000)

	# ----------------------------------------------------------------------
	# create a new zone definition on the location we clicked
	# ----------------------------------------------------------------------
	def actionZone(self, x, y):
		u,v = self._viewer.pixel2uv(x, y)
		self.parent.stopThreads()
#		self._viewer.zone("select",(u,v,BIT_SELFREEZE))

		zone = self._viewer.zone("select")
		for view in self.viewports:
			view._viewer.zone("select",zone)
		self.parent.draw()

		if zone is None or len(zone)==0:
			self.parent.setStatusMsg("-empty-")
			return
		szone = list2Zone(zone[0][1])

		xx,yy,zz = self._viewer.pixel2xyz(x, y)

		if Input._developer:
			say("\nMatrix   =\n", str(Matrix(self._viewer.matrix())))
			say("Origin   =",self._viewer.origin())
			say("OriginUV =",self._viewer.offset())
			say("Zoom     =",self._viewer.zoom())
			say("Extends  =",self._viewer.extends())
			say("\nxi=",x)
			say("yi=",y)
			say("u=",repr(u))
			say("v=",repr(v))
			say("x=",repr(xx))
			say("y=",repr(yy))
			say("z=",repr(zz))
			say("Zone=", szone)

		self.parent.setStatusMsg(szone)

		if self.editor.editType() == Input.Card.REGION:
			self.editor.selectionClear(False, Input.Card.BODY)
			self.editor.addZone(szone)
		else:
			self.editor.addRegion(False, expr=szone)

	# ----------------------------------------------------------------------
	# Zoom on selected object
	# ----------------------------------------------------------------------
	def zoomOn(self, xmin, ymin, xmax, ymax, uv=None):
		# Origin
		if uv is None:
			self._viewer.offset(0.5*(xmin+xmax), 0.5*(ymin+ymax))
		else:
			self._viewer.offset(*uv)

		# Zoom and Window
		x1,y1,x2,y2 = self._viewer.extends()
		if abs(xmax-xmin) < 1.0e-8:
			zx = MAXZOOM
		else:
			zx = 0.5*(x2-x1)/(xmax-xmin)

		if abs(ymax-ymin) < 1.0e-8:
			zy = MAXZOOM
		else:
			zy = 0.5*(y2-y1)/(ymax-ymin)
		self._viewer.zoom(min(zx,zy))
		self.draw()
		self._viewer.expose()

	# ----------------------------------------------------------------------
	def zoomOnError(self, err):
		if not (1<=err<=self._viewer.error("n")): return

		error = self._viewer.error("get",err)
		if err is None:
			return

		x,y,z,xmin,xmax,ymin,ymax,body,regIn,regOut = error
		self.zoomOn(xmin,ymin, xmax,ymax, self._viewer.xyz2uv(x,y,z))

	# ----------------------------------------------------------------------
	def panLeft(self, event=None):
		x,y = self._viewer.offset()
		if event is not None and event.state & SHIFT_MASK:
			step = MOVESTEP2
		else:
			step = MOVESTEP
		x -= float(self.winfo_width())/step / self._viewer.Sx()
		self._viewer.offset(x,y)
		self.draw()

	# ----------------------------------------------------------------------
	def panRight(self, event=None):
		x,y = self._viewer.offset()
		if event is not None and event.state & SHIFT_MASK:
			step = MOVESTEP2
		else:
			step = MOVESTEP
		x += float(self.winfo_width())/step / self._viewer.Sx()
		self._viewer.offset(x,y)
		self.draw()

	# ----------------------------------------------------------------------
	def panTop(self, event=None):
		x,y = self._viewer.offset()
		if event is not None and event.state & SHIFT_MASK:
			step = MOVESTEP2
		else:
			step = MOVESTEP
		y += float(self.winfo_height())/step / self._viewer.Sy()
		self._viewer.offset((x,y))
		self.draw()

	# ----------------------------------------------------------------------
	def panBottom(self, event=None):
		x,y = self._viewer.offset()
		if event is not None and event.state & SHIFT_MASK:
			step = MOVESTEP2
		else:
			step = MOVESTEP
		y -= float(self.winfo_height())/step / self._viewer.Sy()
		self._viewer.offset((x,y))
		self.draw()

	# ----------------------------------------------------------------------
	def panFront(self, event=None):
		matrix = Matrix(self._viewer.matrix())
		matrix = Matrix(self._viewer.matrix())
		width = float(self.winfo_width())/self._viewer.Sx()
		if event is not None and event.state & SHIFT_MASK:
			step = MOVESTEP2
		else:
			step = MOVESTEP
		if self._viewer.get("3D") and self._viewer.get("projection"):	# Perspective or Combo
			width = max(abs(self._orbitDepth)/10.0,width)
		self._viewer.matrix(matrix*Matrix.translate(0.0, 0.0, width/2.0/step))
		self.draw(True)

	# ----------------------------------------------------------------------
	def panBack(self, event=None):
		matrix = Matrix(self._viewer.matrix())
		width = float(self.winfo_width())/self._viewer.Sx()
		if event is not None and event.state & SHIFT_MASK:
			step = MOVESTEP2
		else:
			step = MOVESTEP
		if self._viewer.get("3D") and self._viewer.get("projection"):	# Perspective or Combo
			width = max(abs(self._orbitDepth)/10.0,width)
		self._viewer.matrix(matrix*Matrix.translate(0.0, 0.0, -width/2.0/step))
#			-float(self.winfo_width())/step/2.0/self._viewer.Sx()))
		self.draw(True)

	# ----------------------------------------------------------------------
	def rotateX(self, ang):
		self._viewer.viewOto0()
		self._viewer.matrix(Matrix(self._viewer.matrix())*Matrix.rotX(ang))
		self.draw(True)

	# ----------------------------------------------------------------------
	def rotateY(self, ang):
		self._viewer.viewOto0()
		self._viewer.matrix(Matrix(self._viewer.matrix())*Matrix.rotY(ang))
		self.draw(True)

	# ----------------------------------------------------------------------
	def rotateZ(self, ang):
		self._viewer.viewOto0()
		self._viewer.matrix(Matrix(self._viewer.matrix())*Matrix.rotZ(ang))
		self.draw(True)

	# ----------------------------------------------------------------------
	# ZoomIn or increase cursor size
	# ----------------------------------------------------------------------
	def cursorOrZoomIn(self, event):
		if event.state & SHIFT_MASK and self._geometry.cursor("cursor"):
			global CURSORSIZE
			CURSORSIZE = self._geometry.set("cursor")-5
			self._geometry.set("cursor", CURSORSIZE)
			self.expose()
		else:
			self.zoomIn(event)

	# ----------------------------------------------------------------------
	# ZoomOut or decrease cursor size
	# ----------------------------------------------------------------------
	def cursorOrZoomOut(self, event):
		if event.state & SHIFT_MASK and self._geometry.cursor("cursor"):
			global CURSORSIZE
			CURSORSIZE = self._geometry.set("cursor")+5
			self._geometry.set("cursor", CURSORSIZE)
			self.expose()
		else:
			self.zoomOut(event)

	# ----------------------------------------------------------------------
	def zoomIn(self, event):
		if event.state & SHIFT_MASK:
			zoom = SZOOM
		else:
			zoom = ZOOM
		# Hold Alt to zoom on center!
		self.zoomPin(zoom, event.x, event.y, not bool(event.state&ALT_MASK))
		self.draw()
		self.parent.expose()	# all viewports

	# ----------------------------------------------------------------------
	def zoomOut(self, event):
		if event.state & SHIFT_MASK:
			zoom = 1.0 / SZOOM
		else:
			zoom = 1.0 / ZOOM
		# Hold Alt to zoom on center!
		self.zoomPin(zoom, event.x, event.y, not bool(event.state&ALT_MASK))
		self.draw()
		self.parent.expose()	# all viewports

	# ----------------------------------------------------------------------
	def wheel(self, event):
		zoom = math.pow(ZOOM*(event.delta//120))
		self.zoomPin(zoom, event.x, event.y, not bool(event.state&ALT_MASK))
		self.draw()

	# ----------------------------------------------------------------------
	# Zoom around a pin point
	# ----------------------------------------------------------------------
	def zoomPin(self, zoomFactor, x, y, pin=True):
		if (PINZOOM and pin) or (not PINZOOM and not pin):
			pu, pv = self._viewer.pixel2uv(x, y)
			ou, ov = self._viewer.offset()
			du = (ou - pu) / zoomFactor
			dv = (ov - pv) / zoomFactor
			self._viewer.offset(pu + du, pv + dv)
		self._viewer.zoom(self._viewer.zoom() * zoomFactor)
		self.parent.expose()	# all viewports

	# ----------------------------------------------------------------------
	# Center all viewports on mouse location
	# ----------------------------------------------------------------------
	def centerViewports(self, event):
		if self._viewer.get("3D"):
			u,v = self._viewer.pixel2uv(event.x, event.y)
			xyz = self._viewer.hit3D(u,v)
		elif self.snap(event):
			uv  = self._viewer.pixel2uv(event.x, event.y)
			xyz = self._viewer.snap(uv[0], uv[1])
		else:
			xyz = self._viewer.pixel2xyz(event.x, event.y)
		for v in range(3):
			if not self.viewports[v].set("3D"):
				self._viewer.viewport(v, "center", xyz)
		self.parent.expose()	# all viewports
		for v in self.viewports:
			v.draw(True,True)

	# ----------------------------------------------------------------------
	# Transform viewport to go to the lattice/prototype position
	# ----------------------------------------------------------------------
	def go2Lattice(self):
		# It has to be a REGION with a ROT-DEFI associated
		if not self.editor.editType() == Input.Card.REGION: return
		rotdefi = self.editor.editObject.cards[0]["@rotdefi"]
		if not rotdefi: return
		trans = self.input().getTransformation(rotdefi)
		self._viewer.viewOto0()
		matrix = Matrix(self._viewer.matrix())
		self._viewer.matrix(trans * matrix)
		self.draw(True, True)

	# ----------------------------------------------------------------------
	def zoomOnSelected(self, objects=None):
		if objects is None:
			objects = self.editor.getSelection()
			if objects is None: return

		cardTypes = {
			Input.Card.BODY:   "B",
			Input.Card.REGION: "R",
			Input.Card.OBJECT: "O" }

		# Search in 2 and 3D
		bbox  = None
		bbox2 = None
		bbox3 = None

		for card in objects:
			t = cardTypes.get(card.type())
			if t is None: continue
			if card[ID] is None: continue

			# First in 2D
			bb = self._viewer.bbox2D(t,card[ID])
			if bbox2 is None:
				bbox2 = bb
			elif bb is not None:
				# make the union with bbox2
				bbox2[0] = min(bbox2[0], bb[0])
				bbox2[1] = min(bbox2[1], bb[1])
				bbox2[2] = max(bbox2[2], bb[2])
				bbox2[3] = max(bbox2[3], bb[3])

			# Then in Viewer3D
			bb = self._viewer.bbox(t,card[ID])
			if bbox3 is None:
				bbox3 = bb
			elif bb is not None:
				# make the union with bbox
				bbox3[0] = min(bbox3[0], bb[0])
				bbox3[1] = min(bbox3[1], bb[1])
				bbox3[2] = min(bbox3[2], bb[2])
				bbox3[3] = max(bbox3[3], bb[3])
				bbox3[4] = max(bbox3[4], bb[4])
				bbox3[5] = max(bbox3[5], bb[5])

			# Then in Geometry 3D
			bb = self._geometry.bbox(t,card[ID])
			if bbox is None:
				bbox = bb
			elif bb is not None:
				# make the union with bbox
				bbox[0] = min(bbox[0], bb[0])
				bbox[1] = min(bbox[1], bb[1])
				bbox[2] = min(bbox[2], bb[2])
				bbox[3] = max(bbox[3], bb[3])
				bbox[4] = max(bbox[4], bb[4])
				bbox[5] = max(bbox[5], bb[5])

		if bbox2 is not None:
			self.zoomOn(*bbox2)

		elif bbox3 is not None:
			x = (bbox[0] + bbox[3])/2.0
			y = (bbox[1] + bbox[4])/2.0
			z = (bbox[2] + bbox[5])/2.0
			self._viewer.viewOto0()
			self._viewer.origin(x,y,z)
			#self.draw(True,True)

			# get limits
			umin = bbox3[0]
			vmin = bbox3[1]
			umax = bbox3[3]
			vmax = bbox3[4]

			# Check for infinite
			x1,y1,x2,y2 = self._viewer.extends()
			if abs(umin)>1e10: umin = x1
			if abs(vmin)>1e10: vmin = y1
			if abs(umax)>1e10: umax = x2
			if abs(vmax)>1e10: vmax = y2

			self.zoomOn(umin, vmin, umax, vmax)

	# ----------------------------------------------------------------------
	# Views History
	# ----------------------------------------------------------------------
	def historyAdd(self):
		t = time.time()
		if self.curHistory>0 and t-self.history[-1][0]<2.0:
			del self.history[-1]
			self.curHistory -= 1

		self.curHistory += 1

		if self.curHistory == MAXHISTORY:
			del self.history[0]
			self.curHistory -= 1

		if 0 <= self.curHistory < len(self.history):
			del self.history[self.curHistory:]


		m = self._viewer.matrix()
		e = self._viewer.extends()
		o = self._viewer.offset()
		z = self._viewer.zoom()

		self.history.append((t,m,e,o,z))
		self.frame.historyButton(self.curHistory>0, False)

	# ----------------------------------------------------------------------
	def historyShow(self):
		t,m,e,o,z = self.history[self.curHistory]
		self._viewer.matrix(m)
		self._viewer.extends(e)
		self._viewer.offset(o)
		self._viewer.zoom(z)
		self.draw(True,True,False)
		self.frame.historyButton(self.curHistory>0,
			self.curHistory<len(self.history)-1)

	# ----------------------------------------------------------------------
	def historyBackward(self):
		if self.curHistory > 0:
			self.curHistory -= 1
			self.historyShow()

	# ----------------------------------------------------------------------
	def historyForward(self):
		if self.curHistory < len(self.history)-1:
			self.curHistory += 1
			self.historyShow()

	# ----------------------------------------------------------------------
	# Menus
	# ----------------------------------------------------------------------
	def popupMenu(self, event):
		menu = Menu(self, tearoff=0)

		menu.add_command(label="Repeat Last Action", underline=0,
				command=self.parent.repeatLastAction,
				compound=LEFT, image=tkFlair.icons["empty"])

		# --- Std cut & copy & paste menus ---
		menu.add_command(label="Cut", underline=2,
				command=self.editor.cut,
				accelerator="Ctrl-X",
				compound=LEFT, image=tkFlair.icons["cut"])
		menu.add_command(label="Copy", underline=0,
				command=self.editor.copy,
				accelerator="Ctrl-C",
				compound=LEFT, image=tkFlair.icons["copy"])
		menu.add_command(label="Paste", underline=0,
				command=self.editor.paste,
				accelerator="Ctrl-V",
				compound=LEFT, image=tkFlair.icons["paste"])
		menu.add_command(label="Clone", underline=1,
				command=self.editor.clone,
				accelerator="Ctrl-D",
				compound=LEFT, image=tkFlair.icons["clone"])
		menu.add_separator()

		# --- Add Menus ---
		submenu = Menu(menu)
		menu.add_cascade(label="Add", menu=submenu, underline=0,
				compound=LEFT, image=tkFlair.icons["add"],
				accelerator="Space")
		self.parent._addMenus(submenu)

		menu.add_command(label="Delete", underline=1,
				command=self.editor.delete,
				accelerator="Del",
				compound=LEFT, image=tkFlair.icons["x"])

		# --- General ---
		menu.add_separator()
		if self.editor.editType() == Input.Card.REGION:
			menu.add_command(label="Go to Lattice", underline=0,
					command=self.go2Lattice,
					accelerator="l",
					compound=LEFT, image=tkFlair.icons["LATTICE"])

		menu.add_command(label="Zoom On Selected", underline=0,
				command=self.zoomOnSelected,
				accelerator="f",
				compound=LEFT, image=tkFlair.icons["zoom_on"])

		menu.add_separator()
		self.parent._visibilityMenu(menu)
		self.parent._wireframeMenu(menu)
		self.parent._freezeMenu(menu)

		menu.tk_popup(event.x_root, event.y_root)

	# ----------------------------------------------------------------------
	def _exportImage(self, filename):
		if not have_PIL:
			self.flair().notify("Missing PIL module",
				"PIL libraries are missing.\nCannot export the image.",
				tkFlair.NOTIFY_ERROR)
			return False

		data = self._viewer.image("get")
		image = PIL.Image.frombuffer("RGBA",
			self._viewer.size(),
			data, "raw", "RGBA", 0, 1)
		try:
			image.save(filename)
		except:
			self.flair().notify("Python Imaging Error",
				"Unknown type: %s"%(sys.exc_info()[1]),
				tkFlair.NOTIFY_ERROR)
			error = True
		del image
		del data
		# restore the image data
		self._viewer.image("restore")

	# ----------------------------------------------------------------------
	# Export view to file
	# ----------------------------------------------------------------------
	def export(self, event=None):
		filename = bFileDialog.asksaveasfilename(master=self,
			title="Export viewport as",
			filetypes=[
				("PNG - Portable Network Graphics", "*.png"),
				("GIF - Graphics Interchange Format", "*.gif"),
				("JPG - Joint Photographic Experts Group", ("*.jpg","*.jpeg")),
				("DXF - Data Exchange Format","*.dxf"),
				("SVG - Scalable Vector Graphics","*.svg"),
				("NOX - No-X11 Format","*.nox"),
				#("PS  - Postscript","*.ps"),
				#("EPS - Encapsulated Postscript","*.eps"),
				("All","*")])

		if not filename: return

		fn,ext = os.path.splitext(filename)
		error = False

		if ext in (".dxf", ".svg"):
			self._viewer.export(filename)

		elif ext == ".nox":
			self._saveNox(filename)

		else:
			error = self._exportImage(filename)

		if not error:
			self.flair().notify("Export to Image",
				"Viewport exported to %r"%(filename))

	# ----------------------------------------------------------------------
	# Save view port into notes
	# ----------------------------------------------------------------------
	def toNotes(self, event=None):
		if not have_PIL:
			self.flair().notify("Missing PIL module",
				"PIL libraries are missing.\nCannot export the image.",
				tkFlair.NOTIFY_ERROR)
			return
		filename = "%s-%s.png"%(self.frame.name, self.frame.layerCombo.get())
		error = self._exportImage(filename)

		link = "{img:%s}"%(filename)
		if link not in self.project().notes:
			self.project().notes += "\n\n%s"%(link)

		self.flair().setModified()
		self.flair().refresh("image")
		self.flair().notify("Viewport to Notes",
			"Viewport %s added to notes\nwith the name: %s"%(self.frame.name, filename))

	# ----------------------------------------------------------------------
	# Read a file with camera coordinates and move the viewport accordingly
	# exporting one image per location
	# ----------------------------------------------------------------------
	def exportMovie(self, event=None):
		filename = bFileDialog.askopenfilename(master=self.parent,
			title="Movie: camera locations",
			filetypes=[
				("Camera positions", ("*.dat","*.txt","*.csv")),
				("All","*")])

		if not filename: return

		fn,ext = os.path.splitext(filename)
		error = False

		f = open(filename,"r")
		if not f:
			self.flair().notify("Error opening file",
				"Cannot open file %r"%(filename),
				tkFlair.NOTIFY_ERROR)
			return

		pos = Vector()
		w   = Vector(0.,0.,1.)
		v   = Vector(0.,1.,0.)
		u   = Vector(1.,0.,0.)
		extu = 100.
		extv = 100.
		M = Matrix(4)
		mask = ~(DRAW_SEGMENTS|DRAW_REGIONS)	# Draw everything apart segments
		drawtime = self._viewer.get("drawtime")
		self._viewer.set("drawtime", 0)
		say("Creating frames:")
		try:
			frame = 0
			for i,line in enumerate(f):
				line = line.strip()
				if not line: continue
				if line[0] in ("#","*"): continue
				try:
					value = list(map(float,line.replace(","," ").split()))
				except ValueError:
					say("Error line #",i,":",line)
					continue
				if len(value) < 3:
					say("Error line #",i,":",line)
					continue

				pos[0] = value[0]
				pos[1] = value[1]
				pos[2] = value[2]

				if len(value) >= 6:
					w[0] = -value[3]
					w[1] = -value[4]
					w[2] = -value[5]

					if len(value) >= 9:
						v[0] = value[6]
						v[1] = value[7]
						v[2] = value[8]

						if len(value) >= 11:
							extu = value[9]
							extv = value[10]

				# Create vectors
				w.normalize()
				v.normalize()

				u = v ^ w
				u.normalize()

				v = w ^ u
				v.normalize()

				# Set position
				M.make(u,v,w)
				self._viewer.matrix(M.transpose())
				self._viewer.origin(pos[0], pos[1], pos[2])
				self._viewer.extends(-2.0*extu, -2.0*extv, 2.0*extu, 2.0*extv)
				self._viewer.zoom(2.0)

				# Draw
				start = time.time()

				self._viewer.draw(False,mask)	# draw as thread
				self._viewer.expose()
				data = self._viewer.image("get")

				# filename
				framename = "%s%04d.png"%(fn,frame)
				frame += 1

				say(">>>",frame,os.path.basename(framename),time.time()-start,"s")

				image = PIL.Image.frombuffer("RGBA",
					self._viewer.size(),
					data, "raw", "RGBA", 0, 1)
				try:
					image.save(framename)
				except:
					self.flair().notify("Python Imaging Error",
						sys.exc_info()[1],
						tkFlair.NOTIFY_ERROR)
					break
				# restore the image data
				self._viewer.image("restore")
		except KeyboardInterrupt:
			pass
		self._viewer.set("drawtime",drawtime)
		f.close()

	# ---------------------------------------------------------------------
	def controlPress(self, event):
		if self.action in (ACTION_VIEW_CENTER, \
				ACTION_MOVE, ACTION_ADD, ACTION_ADD_NEXT):
			# When Control is pressed check for the appropriate action
			event.state |= CONTROL_MASK
			self.buttonMotion(event)

	# ---------------------------------------------------------------------
	def controlRelease(self, event):
		if self.action in (ACTION_VIEW_CENTER, \
				ACTION_MOVE, ACTION_ADD, ACTION_ADD_NEXT):
			# When Control is pressed check for the appropriate action
			event.state &= ~CONTROL_MASK
			self.buttonMotion(event)

	# ---------------------------------------------------------------------
	# Scan other viewports for errors along one axis, stop on first error
	#
	# FIXME should be done in parallel and stop on errors
	# ---------------------------------------------------------------------
	def scan4Errors(self):
		# FIXME somehow the user should enter the requested axis!
		axis = 'x'
		idx = ord(axis)-ord('x')
		origins = [list(x._viewer.origin()) for x in self.viewports]
		vertices = self._viewer.vertex(axis)
		if len(vertices)<2: return
		a = vertices[0]
		for b in vertices[1:]:
			m = 0.5*(a+b)
			for o,v in zip(origins, self.viewports):
				o[idx] = m
				v._viewer.origin(o)
				v.draw(True)
			a = b
		for o,v in zip(origins, self.viewports):
			v._viewer.origin(o)
			v.draw(True)

	# ---------------------------------------------------------------------
#	def _scan4ErrorsOnAxis(self, axis, vertices):
#		origin = list(self._viewer.origin())
#		a = ord(axis)-ord('x')
#		for v in vertices:
#			origin[a] = v
#			self._viewer.origin(origin)
#			self.draw(True)
#			# wait to finish and check for errors

	# ---------------------------------------------------------------------
	# Handle keyboard commands
	# ---------------------------------------------------------------------
	def handleKey(self, event):
		if event.keysym in ("Shift_L", "Shift_R", "Alt_L", "Alt_R"):
			return

		# Cancel action
		elif event.keysym == "Escape":
			if self.action in (ACTION_MOVE, ACTION_ROTATE) and self._anchorValue != "":
				self._anchorValue = ""
				# Call action Move or Rotate to update
				if self.action == ACTION_ROTATE:
					self.actionRotate(event)
				else:
					self.actionMove(event)
				return

			self.grab_release()
			# 1st stop current Action
			if self.action != ACTION_SELECT or self.parent.mouseAction < 0:
				if self.parent.mouseAction == ACTION_ADD_ZONE:
					self.parent.zone("clear")
				self.actionStop()
				self.parent.setStatusMsg("Cancel action")

			# 2nd stop any pending action
			elif self.parent.mouseAction in (ACTION_ADD, ACTION_ADD_NEXT, ACTION_EDIT):
				if self.parent.mouseAction == ACTION_ADD_ZONE:
					self.parent.zone("clear")
				self.editor.setMouseAction(ACTION_SELECT)
				self.parent.setStatusMsg("Cancel action")

			# 3rd if a region is selected
			elif self.editor.editType() == Input.Card.REGION:
				# 3.1 clear property list selection if any
				if self.editor.propList.curselection():
					self.editor.propList.selection_clear(0,END)
					self.editor.propListSelect()
					self.parent.zone("clear")
					self.parent.setStatusMsg("Properties cleared")

				# 3.2 clear objects
				elif self.editor.isAnythingSelected(Input.Card.OBJECT):
					self.editor.selectionClear(True,Input.Card.OBJECT)
					self.parent.setStatusMsg("Objects unselected")

				# 3.3 clear bodies
				elif self.editor.isAnythingSelected(Input.Card.BODY):
					self.editor.selectionClear(True,Input.Card.BODY)
					self.parent.setStatusMsg("Bodies unselected")

				# 3.4 clear everything else
				elif self.editor.isAnythingSelected():
					self.editor.selectionClear()
					self.parent.setStatusMsg("Unselected all")

			# 4rd unselect anything selected
			elif self.editor.isAnythingSelected():
				self.editor.selectionClear()
				self.parent.setStatusMsg("Unselected all")

			self.parent.updateEdit()
			self.parent.draw()
			self.parent.expose()

		# Move
		elif event.keysym == "Left":
			if event.state & CONTROL_MASK:
				if event.state & SHIFT_MASK:
					ang = ROTSTEP2
				else:
					ang = ROTSTEP
				self.rotateY(ang)
			elif event.state & ALT_MASK:
				if event.state & SHIFT_MASK:
					self.actionMoveBy(event, -0.1, 0)
				else:
					self.actionMoveBy(event, -1, 0)
			else:
				self.panLeft(event)

		elif event.keysym == "Right":
			if event.state & CONTROL_MASK:
				if event.state & SHIFT_MASK:
					ang = ROTSTEP2
				else:
					ang = ROTSTEP
				self.rotateY(-ang)
			elif event.state & ALT_MASK:
				if event.state & SHIFT_MASK:
					self.actionMoveBy(event, 0.1, 0)
				else:
					self.actionMoveBy(event, 1, 0)
			else:
				self.panRight(event)

		elif event.keysym == "Up":
			if event.state & CONTROL_MASK:
				if event.state & SHIFT_MASK:
					ang = ROTSTEP2
				else:
					ang = ROTSTEP
				self.rotateX(ang)
			elif event.state & ALT_MASK:
				if event.state & SHIFT_MASK:
					self.actionMoveBy(event, 0, 0.1)
				else:
					self.actionMoveBy(event, 0, 1)
			else:
				self.panTop(event)

		elif event.keysym == "Down":
			if event.state & CONTROL_MASK:
				if event.state & SHIFT_MASK:
					ang = ROTSTEP2
				else:
					ang = ROTSTEP
				self.rotateX(-ang)
			elif event.state & ALT_MASK:
				if event.state & SHIFT_MASK:
					self.actionMoveBy(event, 0, -0.1)
				else:
					self.actionMoveBy(event, 0, -1)
			else:
				self.panBottom(event)

		elif event.keysym == "Prior":
			if event.state & CONTROL_MASK:
				if event.state & SHIFT_MASK:
					ang = ROTSTEP2
				else:
					ang = ROTSTEP
				self.rotateZ(-ang)
			else:
				self.panFront(event)

		elif event.keysym == "Next":
			if event.state & CONTROL_MASK:
				if event.state & SHIFT_MASK:
					ang = ROTSTEP2
				else:
					ang = ROTSTEP
				self.rotateZ(ang)
			else:
				self.panBack(event)

		elif self.action in (ACTION_MOVE, ACTION_ROTATE, ACTION_ADD_NEXT, ACTION_ORBIT):
			self._handleKeyMove(event)

		else:
			self._handleKeyNormal(event)

	# ---------------------------------------------------------------------
	# Handle key method for all other actions
	# ---------------------------------------------------------------------
	def _handleKeyNormal(self, event):
		invalid = False

		if event.keysym == "space":
			if event.state & CONTROL_MASK:
				self.popupMenu(event)
			else:
				self.parent.insertMenu()
			return

		elif event.keysym == "Menu":
			self.popupMenu(event)
			return

		elif event.char == 'a':
			# with or without control
			if self.editor.editObject is not None and \
			   self.editor.editObject.type() == Input.Card.REGION:
				# Select all bodies from all regions
				if not self.editor.selectAllBodies():
					self.editor.objList.selectAll()
			else:
				self.editor.objList.selectAll()

		elif event.char == 'A':
			# with or without control
			self.editor.objList.selectClear()

		elif event.char == 'b':
			self.parent.insertBodyMenu()

		elif event.char == 'B':
			if self.parent.actionOption:
				self.parent.newBody()
			else:
				self.parent.insertBodyMenu()

		elif event.char == 'd':
			if self.parent.actionTool is self.parent._zoneTool:
				self.parent.actionTool.setAction(ACTION_ADD_ZONE)
			else:
				self.editor.setMouseAction(ACTION_ZONE)
			self.actionZoneShow(event.x, event.y)
			return

#		elif event.char == 'D':
#			if Input._developer:
#				self.editor.setMouseAction(ACTION_ZONEPAINT)
#			elif self.parent.setActionNewZone():
#				self.actionZoneShow(event.x, event.y)
#			return

		elif event.char == 'G':
			self.editor.toggleSnap2Grid()

		elif event.char == 'c':
			self.centerViewports(event)
			return

		elif event.char == 'C':
			self.parent.toggleAction(ACTION_VIEW_CENTER)
			return

		elif event.char == 'e':
			self.editor.editProperty()
			return

		elif event.char == 'E':
			self.parent.toggleAction(ACTION_EDIT)
			return

		elif event.char == 'f':	# find
			self.zoomOnSelected()
			return

		elif event.char == 'F':
			self.editor.toggleFreeze()

		elif event.char == 'g':
			self._closest = None
			self.resetAnchor()
			self.actionMoveStart(event, 0)
			self.actionMove(event)
			return

		elif event.char == 'i':
			self.parent.toggleAction(ACTION_INFO)
			return

		elif event.char == 'K':
			if Input._developer:
				say("\nMemory Information:",self._viewer.title())
				self._geometry.memory("dump")
				self._viewer.memory("dump")

				# Self destroy for MEMORY leak testing
				# needs to be compiled with MEM=yes
				try:
					for view in self.parent.viewers:
						view._destroy()
					self._geometry._destroy()
				except AttributeError:
					pass
				return

		elif event.char == 'l':
			self.go2Lattice()

		elif event.char == 'L':
			self.editor.toggleLock()

		elif event.char == 'm':
			self.editor.editMaterial()

		elif event.char == 'M':
			self.exportMovie()

		elif event.char == 'n' or event.keysym == "F2":
			self.editor.editName()
			return

		elif event.char == 'o':
			self.master.toggleOrigin()

		elif event.char == 'O':
			self.master.nextOrigin()

		elif event.char == 'p':
			self.parent.newPoint()
			return

		elif event.char == 'r':		# rotate object
			self.resetAnchor()
			self.actionRotateStart(event)
			return

		elif event.char == 'R':
			self.editor.addRegion()
			return

		elif event.char == 's':
			self.parent.toggleAction(ACTION_SELECT)
			return

		elif event.char == 'S':
			if Input._developer: self.scan4Errors()
			return

		elif event.char == 't':
			self.parent.toggleAction(ACTION_ORBIT)
			return

		elif event.char == 'v':
			self.editor.toggleVisibility()

		elif event.char == 'V':
			self.editor.toggleAlpha()

		elif event.char == 'w':
			self.parent.rotateLayout()

		elif event.char == 'W':
			self.parent.customLayout()

		elif event.char == 'x':
			self.parent.toggleAction(ACTION_PAN)
			return

		elif event.char == 'X':
			self.frame.rotateLayer()

		elif event.char == 'z':
			self.parent.toggleAction(ACTION_ZOOM_IN)
			return

		elif event.char == 'Z':
			self.parent.toggleAction(ACTION_ZOOM_OUT)
			return

		elif event.keysym in ("Insert", "Return", "KP_Enter"):
			self.parent.insertMenu()
			return

		elif event.keysym in ("Delete", "BackSpace"):
			self.editor.delete()

		elif event.char == '.':
			self.parent.repeatLastAction()
			return

		elif event.char == "#":
			self.parent.layout("2x2")

		elif event.char == "%":
			self.editor.toggle3D()
			return

		# Zoom
		elif event.char in ('=','+'):
			self.zoomIn(event)
			return

		elif event.char in ('-','_'):
			self.zoomOut(event)
			return

		elif event.char == '<':
			self.historyBackward()

		elif event.char == '>':
			self.historyForward()

		# Views
		elif event.keysym == '1' and event.state & CONTROL_MASK:
			xyz = self._viewer.origin()
			self._viewer.matrix(FRONT_MATRIX)
			self._viewer.origin(xyz)
			invalid = True

		elif event.keysym == '2' and event.state & CONTROL_MASK:
			xyz = self._viewer.origin()
			self._viewer.matrix(BACK_MATRIX)
			self._viewer.origin(xyz)
			invalid = True

		elif event.keysym == '3' and event.state & CONTROL_MASK:
			xyz = self._viewer.origin()
			self._viewer.matrix(TOP_MATRIX)
			self._viewer.origin(xyz)
			invalid = True

		elif event.keysym == '4' and event.state & CONTROL_MASK:
			xyz = self._viewer.origin()
			self._viewer.matrix(BOTTOM_MATRIX)
			self._viewer.origin(xyz)
			invalid = True

		elif event.keysym == '5' and event.state & CONTROL_MASK:
			xyz = self._viewer.origin()
			self._viewer.matrix(LEFT_MATRIX)
			self._viewer.origin(xyz)
			invalid = True

		elif event.keysym == '6' and event.state & CONTROL_MASK:
			xyz = self._viewer.origin()
			self._viewer.matrix(RIGHT_MATRIX)
			self._viewer.origin(xyz)
			invalid = True

		elif event.keysym == '0'  and event.state & CONTROL_MASK or \
		     event.keysym == "Home":		# center screen
			self._viewer.viewOto0()
			self._viewer.origin(0.0, 0.0, 0.0)
			invalid = True

		# Layer
		elif event.keysym in GeometryLayers.SHORTCUTS:
			for layer in list(self.project().geometry.layers.values()):
				try:
					if layer["Options.shortcut"] == event.keysym:
						self.frame.layerCombo.set(layer.name)
						return
				except KeyError:
					pass

		else:
			if Input._developer:
				say("GeometryViewer::handleKey: char=%s keysym=%s (%X) state=%X" % \
					(event.char, event.keysym, event.keycode, event.state))

			return

		self.draw(invalid)

	# ---------------------------------------------------------------------
	# Handle key method for move actions
	# ---------------------------------------------------------------------
	def _handleKeyMove(self, event):
		if event.char == 'G':
			self.editor.toggleSnap2Grid()

		#elif event.char == 'u':
			#self.setAnchorAxis(event.char)

		#elif event.char == 'v':
			#self.setAnchorAxis(event.char)

		#elif event.char == 'w':
			#self.setAnchorAxis(event.char)

		elif event.char == 'x':
			self.setAnchorAxis(event.char)

		elif event.char == 'y':
			self.setAnchorAxis(event.char)

		elif event.char == 'z':
			self.setAnchorAxis(event.char)

		# Zoom
		elif event.char in ('0','1','2','3','4','5','6','7','8','9','.','+','-','e','E',' '):
			self._anchorValue += event.char
			self.parent.setStatusMsg("Value: %s"%(self._anchorValue),"Blue","White")

		elif event.char in ("*","/"):
			if self._anchorValue and self._anchorValue[-1] != ' ':
				self._anchorValue += ' '	# add a space
			self._anchorValue += event.char
			self.parent.setStatusMsg("Value: %s"%(self._anchorValue),"Blue","White")

		elif event.char == ",":
			self._anchorValue += ' '	# replace to space
			self.parent.setStatusMsg("Value: %s"%(self._anchorValue),"Blue","White")

		elif event.keysym  == "BackSpace":
			if self._anchorValue:
				self._anchorValue = self._anchorValue[:-1]
				self.parent.setStatusMsg("Value: %s"%(self._anchorValue),"Blue","White")

		elif event.keysym in ("Return", "KP_Enter"):
			# Equivalent to release button
			self.buttonRelease(event)
			return

		elif Input._developer:
			say("GeometryViewer::handleKeyMove: char=%s keysym=%s state=%X" % \
				(event.char, event.keysym, event.state))

		# Call action Move or Rotate to update
		if self.action == ACTION_ROTATE:
			self.actionRotate(event)
		else:
			self.actionMove(event)

	# ---------------------------------------------------------------------
	# Save input&view to a special file to be used by the test program
	# ---------------------------------------------------------------------
	def _saveNox(self, filename):
		fout = open(filename,"w")
		inp = self.input()
		inp.preprocess()
		inp.regionProperties()
		rotdefi = inp.getTransformation()

		# Write matrix sub-function
		def writeMatrix(matrix):
			for row in range(4):
				for col in range(4):
					fout.write(" %.15g"%(matrix[row][col]))
			fout.write("\n")

		# Write Transformations
		fout.write("# Transformations\n")
		for key,matrix in list(rotdefi.items()):
			if isinstance(key,int):
				fout.write("rotdefi %d"%(key))
				writeMatrix(matrix)
		fout.write("\n")

		# --------------------
		#       Bodies
		# --------------------
		fout.write("# Bodies\n")
		hasVoxel = False
		matrix = [Matrix(4,type=1)]	# list of transformation matrices
		for card in inp.cardlist:
			if not card.isGeo(): continue
			if card.ignore(): continue
			if len(card.tag)==3 and card.isGeo():
				if card.tag=="END": break

				# Add the body
				if len(matrix)>1:
					fout.write("matrix ")
					writeMatrix(matrix[-1])
				fout.write("body %s %s %s\n"%(card.tag, card.sdum(),
					" ".join(map(str,
						[card.numWhat(x) for x in range(1, card.nwhats())]))))

			elif card.tag[0]=="$":
				if card.tag=="$start_translat":
					matrix.append(matrix[-1]*Matrix.translate(
						card.numWhat(1), card.numWhat(2), card.numWhat(3)))
				elif card.tag=="$start_transform":
					# To be compatible with FLUKA
#					matrix.append(matrix[-1] *
#						inp.getTransformation(card.evalWhat(1),rotdefi))
					matrix.append(matrix[-1]*input.getTransformation(card.evalWhat(1),rotdefi))

				elif card.tag=="$start_expansion":
					matrix.append(matrix[-1]*Matrix.scale(card.numWhat(1)))
				elif card.tag[:5]=="$end_":
					if len(matrix)>1:
						del matrix[-1]
					else:
						self.addError("Too many $end_xxx cards")

			elif card.tag == "VOXELS":
				if card.what(4):
					mat = inp.getTransformation(card.evalWhat(4),rotdefi)
					fout.write("matrix ")
					writeMatrix(mat)
				fout.write("body %s %s.vxl %g %g %g\n"%(card.tag, card.sdum().strip(),
					card.numWhat(1), card.numWhat(2), card.numWhat(3)))
				hasVoxel = True
		fout.write("\n")

		# --------------------
		#       Regions
		# --------------------
		def writeRegion(region,expr):
			rot = signedRotdefi(rotdefi, region["@rotdefi"])

			color = self._viewer.region(region.sdum(),"color")
			alpha = self._viewer.region(region.sdum(),"alpha")

			sex = " ".join(["%s %s"%(m," ".join(e))
				for m,e in regionExp(csg.splitZones(csg.tokenize(expr)))]," : ")

			fout.write("region %s %d %d %d %d %s\n" %\
				(region.sdum(), region["@type"], rot, color, alpha, sex))

		fout.write("# Regions\n")
		i = 0
		regionCards = inp.cardsSorted("REGION")
		n = len(regionCards)
		expr  = None
		color = None
		while i < n:
			card = regionCards[i]
			i += 1
			if card.ignore(): continue

			if card.sdum() == "&": # append to previous expression
				expr += card.extra()
				continue

			if expr is not None:
				writeRegion(region,expr)

			# move old card to region and process new card
			region = card
			expr   = card.extra()


		# Add last region
		writeRegion(region,expr)
		fout.write("\n")

		# Write voxel information
		if hasVoxel:
			fout.write("# Voxel information\n")
			for i in range(1,self._geometry.voxel("no")+1):
				color = self._viewer.voxel("color",i)
				fout.write("voxel color %d %d\n"%(i,color))
			fout.write("\n")

		# Write view
		fout.write("# View\n")
		fout.write("viewmatrix ")
		writeMatrix(self._viewer.matrix())
		fout.write("extends %g %g %g %g\n"%self._viewer.extends())
		fout.write("zoom %g\n"%self._viewer.zoom())
		fout.write("origin %g %g %g\n"%self._viewer.origin())

		# Write view
		fout.close()

#===============================================================================
# Geometry frame
#===============================================================================
class GeometryFrame(Frame):
	# ----------------------------------------------------------------------
	def __init__(self, master, parent, name, **kw):
		Frame.__init__(self, master, **kw)
		self.parent = parent
		self.name   = name

		self.viewer = GeometryViewer(self, self, parent, **kw)
		self.viewer.pack(side=BOTTOM, expand=YES, fill=BOTH)
		self.createToolbar()

		# Remember loaded image
		self._clearImage()

	# ----------------------------------------------------------------------
	def flair(self):	return self.parent.editor.flair
	def project(self):	return self.parent.project
	def input(self):	return self.parent.project.input

	# ----------------------------------------------------------------------
	def createToolbar(self):
		toolbar = Frame(self, relief=RAISED, takefocus=0, borderwidth=2)
		toolbar.pack(side=TOP, fill=X)

		self.maxWindow = IntVar()

		self.comboButton = Button(toolbar,image=tkFlair.icons["axes"],
				text=" "+Unicode.BLACK_DOWN_POINTING_TRIANGLE,
				compound=LEFT,
				command=self.changeOrigin,
				padx=0, pady=0)
		self.comboButton.grid(row=0, column=0, sticky=N+EW)
		self._lastOrigin = 1
		tkExtra.Balloon.set(self.comboButton, "Origin and projection [o]")

		self._frames = []
		self._frames.append(GeometryExtra._NavigationFrame(toolbar,self))
		self._frames.append(GeometryExtra._OriginFrame(toolbar,    self))
		self._frames.append(GeometryExtra._ExtendsFrame(toolbar,   self))
		self._frames.append(GeometryExtra._MoveFrame(toolbar,      self))
		self._frames.append(GeometryExtra._BasisFrame(toolbar,     self))
		self._frames.append(GeometryExtra._EulerFrame(toolbar,     self))
		self._frames.append(GeometryExtra._TransformFrame(toolbar, self))
		self.lastFrame = None
		self.showOrigin(self._frames[0])

		# FIXME dirty hack
		self.layerCombo = self._frames[0].layerCombo

		# max button
		self.maxButton = Checkbutton(toolbar, text=Unicode.BLACK_UP_POINTING_TRIANGLE,
				command=self.maximizeWindow,
				variable=self.maxWindow,
				indicatoron=FALSE,
				highlightthickness=0,
				padx=2, pady=0)
		self.maxButton.grid(column=2, row=0, pady=2, sticky=N+E)
		tkExtra.Balloon.set(self.maxButton, "Maximize viewport [W]")

		toolbar.grid_columnconfigure(1, weight=1)

	# ----------------------------------------------------------------------
	def changeOrigin(self):
		menu = Menu(self, tearoff=0)
		for frame in self._frames:
			if frame is self.lastFrame:
				state = DISABLED
			else:
				state = NORMAL
			menu.add_command(label=frame.name,
					image=tkFlair.icons[frame.icon],
					state=state,
					compound=LEFT,
					command=lambda s=self,f=frame:s.showOrigin(f))
		menu.tk_popup(self.comboButton.winfo_rootx(),
			      self.comboButton.winfo_rooty() + self.comboButton.winfo_height())

	# ----------------------------------------------------------------------
	def toggleOrigin(self):
		idx = self._frames.index(self.lastFrame)
		if idx == 0:
			idx = self._lastOrigin
		else:
			idx = 0
		self.showOrigin(self._frames[idx])

	# ----------------------------------------------------------------------
	def nextOrigin(self):
		idx = self._frames.index(self.lastFrame) + 1
		if idx >= len(self._frames): idx = 0
		self.showOrigin(self._frames[idx])

	# ----------------------------------------------------------------------
	def showOrigin(self, frame):
		if self.lastFrame: self.lastFrame.grid_forget()
		self.comboButton["image"] = tkFlair.icons[frame.icon]
		frame.grid(row=0, column=1, sticky=EW)
		frame.fill()
		self.lastFrame = frame
		idx = self._frames.index(self.lastFrame)
		if idx>0: self._lastOrigin = idx

	# ----------------------------------------------------------------------
	def loadOrigin(self):
		for frame in self._frames:
			frame.load()

	# ----------------------------------------------------------------------
	def fillOrigin(self):
		self.lastFrame.fill()

	# ----------------------------------------------------------------------
	def setColor(self, color):
		hexcolor = "#%06X"%(color)
		self.maxButton.config(background=hexcolor,
				selectcolor=hexcolor,
				activebackground=hexcolor)

	# ----------------------------------------------------------------------
	def draw(self):
		"""Redraw screen"""
		self.viewer.draw()

	# ----------------------------------------------------------------------
	def redraw(self):
		"""Force regeneration of projections & plots"""
		# save current layer
		self.layer()
		self.viewer.draw(True)

	# ----------------------------------------------------------------------
	def changeLayer(self):
		self.viewer.stop()
		self.layer()
		try:
			if self.viewer.winfo_ismapped():
				self.viewer.draw()
		except TclError:
			pass

	# ----------------------------------------------------------------------
	def layersDialog(self, event=None):
		self.parent.showLayer(self.layerCombo.get())

	# ----------------------------------------------------------------------
	def layer(self):
		if self.project() is None: return

		# FIXME
		# Maybe we should avoid reloading if the last layer was
		# the same. Probably with a hash on the layer
		# Problem. On large voxel files it takes too long
		# to update the regions etc..

		self.activeLayer = self.layerCombo.get()
		if self.activeLayer not in self.project().geometry.layers:
			self.activeLayer = "Media"
			if self.layerCombo.get() != self.activeLayer:
				self.layerCombo.set(self.activeLayer)

		self.viewer.lateDraw = True
		_viewer   = self.viewer._viewer
		_geometry = self.viewer._geometry

		_viewer.set("borders", 0)
		_viewer.set("vertex",  0)
		_viewer.set("fill",    0)
		_viewer.set("lattice", 0)
		_viewer.set("voxel",   0)
		_viewer.set("palette", 0)
		_viewer.set("userdump",0)
		_viewer.set("usrbin",  0)
		_viewer.set("image",   0)
		_viewer.set("3D",      0)
		_viewer.set("projection",0)
		_viewer.set("labels", self.getVar("Show.label")!="None")

		_viewer.set("axes",		self.getInt("Options.axes"))
		_viewer.set("title",		self.getInt("Options.title"))
		_viewer.set("viewport",		self.getInt("Options.viewport"))
		_viewer.set("grid",		self.getInt("Options.grid"))
		_viewer.set("errors",		self.getInt("Options.errors"))
		_viewer.set("gridlevel",	self.getInt("Options.gridlevel"))
		_viewer.set("latticelevel",	self.getInt("Options.latticelevel"))
		_viewer.set("voxellevel",	self.getInt("Options.latticelevel"))	# <-- same as lattice
		_viewer.set("crosshair",	self.getInt("Options.crosshair"))
		_viewer.set("textbackground",	self.getInt("Options.textbackground"))

		# Check fonts
		self.setFont("general", self.getVar("Options.font"))
		self.setFont("grid",    self.getVar("Options.gridfont"))
		self.setFont("palette", self.getVar("Options.palettefont"))

		# First set palette colors
		for tag in GeometryLayers.PALETTENAMES:
			if self.getInt(tag+"._visible"):
				self._paletteLayer()
				break
		else:
			# reset all palettes
			for i in range(len(GeometryLayers.PALETTENAMES)):
				_viewer.palette(i,"reset")

		self._showLayer()	# setup information to display

		if self.getInt("Show"):
			if self.getVar("Show.color") != "None":
				_viewer.set("fill",    1)
			_viewer.set("borders", 1)
			_viewer.set("vertex",  self.getInt("Show.vertex"))

		if self.getInt("Image"):
			self._imageLayer()
		elif self._imageData:
			# Free memory
			self._clearImage()

		if self.getInt("3D"):		self._3dLayer()	# check before the usrbin
		if self.getInt("Userdump"):	self._userdumpLayer()

		self.viewer._viewer.voxel("roiclear")
		if self.getInt("Voxel"):	self._voxelLayer()

		for tag in GeometryLayers.USRBINNAMES:
			if self.getInt(tag):
				self._usrbinLayer()
				break

		# Common
		try:
			if self.getInt("Show.lattice"):
				_viewer.set("lattice",1)
		except TypeError:
			pass

		try:
			if self.getInt("Show.voxel") and self.parent.voxels is not None:
				_viewer.set("voxel", 1)
		except TypeError:
			pass

	# ----------------------------------------------------------------------
	# Load fonts
	# ----------------------------------------------------------------------
	def setFont(self, name, font):
		vfont = self.viewer._viewer.font(name)
		if vfont is None or font != os.path.splitext(os.path.basename(vfont))[0]:
			# Find font and load it
			fn = os.path.join(tkFlair.iniDir,"fonts/%s.tga"%(font))
			try: os.stat(fn)
			except: fn = os.path.join(tkFlair.prgDir,"fonts/%s.tga"%(font))

			try: os.stat(fn)
			except: return

			self.viewer._viewer.font(name,fn)

	# ----------------------------------------------------------------------
	def _regionsLayer(self):
		_viewer = self.viewer._viewer
		showlabel = self.getVar("Show.label")
		for card in self.input()["REGION"]:
			if card.ignore(): continue
			if card.sdum()=="&": continue
			try:
				color = self.parent.regionColor(card)
				rid = _viewer.region(card.sdum(), "id")
				_viewer.region(rid, "color", color)
				if showlabel == "Material":
					try: _viewer.region(rid, "label", card["@material"].sdum())
					except KeyError: pass
				else:
					_viewer.region(rid, "label", card.sdum())
			except KeyError:
				pass

	# ----------------------------------------------------------------------
	def _materialsLayer(self, key="@material"):
		_viewer   = self.viewer._viewer
		geometry  = self.parent._geometry
		showlabel = self.getVar("Show.label")

		# Scan regions
		for card in self.input()["REGION"]:
			if card.ignore(): continue
			if card.sdum()=="&": continue
			try:
				mat   = card[key]
				color = self.parent.materialColor(mat)
				rid   = card[ID]
				_viewer.region(rid, "color", color)

				if showlabel in ("Material","Value"):
					_viewer.region(rid, "label", mat.sdum())
				else:
					_viewer.region(rid, "label", card.sdum())

				alpha = int(card.get(ALPHA,0))
				if alpha:
					alpha = int(round((alpha*255)/100))
					_viewer.region(rid, "alpha", alpha)
				elif card["@type"]==Input.REGION_NORMAL:
					alpha = self.parent.materialAlpha(mat)
					alpha = int(round((alpha*255)/100))
					_viewer.region(rid, "alpha", alpha)
				else:
					_viewer.region(rid, "alpha", 0)
			except (KeyError, TypeError):
				pass

		# Scan for voxel assignments if needed
		if self.parent.voxels and (self.getInt("Show.voxel") or self.getInt("3D")):
			matDict = self.input().materialDict()
			_viewer.region("VOXEL", "alpha", 0)

			# set the voxel and user assignmats
			for card in self.input().cardsSorted("ASSIGNMA"):
				if card.ignore(): continue

				if len(card.what(3)) > 0:
					toVoxel = Input._VOXELPAT.match(card.what(3))
					if toVoxel is None: continue
				else:
					toVoxel = None

				fromVoxel = Input._VOXELPAT.match(card.what(2))
				if fromVoxel is None: continue

				fromVoxel = int(fromVoxel.group(1))
				if toVoxel:
					toVoxel = int(toVoxel.group(1))
				else:
					toVoxel = fromVoxel

				step = card.intWhat(4)
				if step<=0: step=1

				for voxel in range(fromVoxel, toVoxel+1, step):
					try:
						material = matDict.get(card.what(1))
						if material:
							alpha = self.parent.materialAlpha(material)
							if alpha == 100:
								# transparent
								color = 0xFFFFFFFF
							else:
								color = self.parent.materialColor(material)
						else:
							color = _LATTICE_COLOR
							if alpha>0:
								color |= int((255*alpha)/100)<<24
						_viewer.voxel("color", voxel, color)
					except KeyError:
						pass

	# ----------------------------------------------------------------------
	def _showLayer(self):
		_viewer   = self.viewer._viewer
		colorValue = self.getVar("Show.color")
		# None will get the default value of 1.0
		if colorValue == "Zone":
			self._regionsLayer()
			return
		elif colorValue == "Region":
			self._regionsLayer()
			return
		elif colorValue == "Material":
			self._materialsLayer()
			return
		elif colorValue == "Material-Decay":
			self._materialsLayer("@matDecay")
			return

		# Setup viewer
		palette = self.getVar("Show.palette")
		try:
			_viewer.palette(GeometryLayers.PALETTENAMES.index(palette), "default", 0)
		except ValueError:
			pass
		cblog = self.getInt(palette+".log")

		propValue = "@%s"%(colorValue.lower())
		if colorValue=="Corrfactor":
			self.input().addCardProperty("REGION",
						"CORRFACT",
						propValue,
						2,4,5,6)
		elif colorValue=="Corrfactor-dE/dx":
			self.input().addCardProperty("REGION",
						"CORRFACT",
						propValue,
						1,4,5,6)
		elif colorValue=="Corrfactor-rho":
			self.input().addCardProperty("REGION",
						"CORRFACT",
						"@corrfactor",
						2,4,5,6)
		elif colorValue=="e-Production":
			self.input().addCardProperty("MATERIAL",
						"EMFCUT",
						propValue,
						1,4,5,6,
						"PROD-CUT")
		elif colorValue=="e-Transport":
			self.input().addCardProperty("REGION",
						"EMFCUT",
						propValue,
						1,4,5,6,"",
						_emfwhat1)
		elif colorValue=="g-Production":
			self.input().addCardProperty("MATERIAL",
						"EMFCUT",
						propValue,
						1,4,5,6,"PROD-CUT",
						_emfwhat1)
		elif colorValue=="g-Transport":
			self.input().addCardProperty("REGION",
						"EMFCUT",
						propValue,
						2,4,5,6,"")
		elif colorValue=="Deltaray":
			self.input().addCardProperty("MATERIAL",
						"DELTARAY",
						propValue,
						1,4,5,6)
		elif colorValue=="Field":
			self.input().addCardProperty("REGION",
						"ASSIGNMA",
						propValue,
						5,2,3,4)

		# ---- @return color value ----
		def findColorValue(card):
			try:
				if colorValue=="Density":
					value = card["@material"].numWhat(3)

				elif colorValue.startswith("Importance"):	# Suffix -H, -E, -N
					part = colorValue[11:]	# particle suffix
					biasing = card["@biasing"+part]
					if biasing is None: biasing = card["@biasingA"]
					if biasing is not None:
						value = biasing.numWhat(3)
					else:
						value = 1.0

				elif colorValue.startswith("Splitting"):
					part = colorValue[10:]	# particle suffix
					biasing = card["@biasingH"]
					if biasing is None: biasing = card["@biasingA"]
					if biasing is not None:
						value = biasing.numWhat(2)
					else:
						value = 1.0

				elif colorValue=="Corrfactor-rho":
					try:
						value = float(card["@corrfactor"])
						if value == 0.0: value = 1.0
					except KeyError: value = 1.0
					value *= card["@material"].numWhat(3)
					#say(card.sdum(), card["@corrfactor"], card["@material"].numWhat(3), value)

				elif colorValue=="e-Production":
					value = float(card["@material"][propValue])

				elif colorValue=="g-Production":
					value = float(card["@material"][propValue])

				elif colorValue=="Deltaray":
					value = float(card["@material"][propValue])

				else:
					value = float(card[propValue])

			except (TypeError, AttributeError):
				value = 1.0

			label = str(value)

			if cblog:
				try: value = math.log10(value)
				except: value = 0.0
			else:
				try: value = float(value)
				except: value = 0.0

			return label, value

		# Loop over regions
		showlabel = self.getVar("Show.label")
		for card in self.input()["REGION"]:
			if card.ignore(): continue
			if card.sdum()=="&": continue
			label, value = findColorValue(card)
			try:
				rid = _viewer.region(card.sdum(), "id")
				_viewer.region(rid, "value", value)
				if showlabel == "Material":
					try: _viewer.region(rid, "label", card["@material"].sdum())
					except KeyError: pass
				elif showlabel == "Value":
					_viewer.region(rid, "label", label)
				else:
					_viewer.region(rid, "label", card.sdum())

				alpha = int(card.get(ALPHA,0))
				if alpha:
					alpha = int(round((alpha*255)/100))
					_viewer.region(rid, "alpha", alpha)
				elif card["@type"]==Input.REGION_NORMAL:
					alpha = self.parent.materialAlpha(card["@material"])
					alpha = int(round((alpha*255)/100))
					_viewer.region(rid, "alpha", alpha)
				else:
					_viewer.region(rid, "alpha", 0)
			except KeyError:
				pass

		# Loop voxel regions
		if self.parent.voxels and (self.getInt("Show.voxel") or self.getInt("3D")):
			for voxel in self.input()["VOXELS"]:
				if voxel.ignore() or voxel["@voxel"] is None: continue
				for card in voxel["@voxel"].input["REGION"]:
					label, value = findColorValue(card)
					pat = Input._VOXELPAT.match(card.name())
					if pat is None: continue
					vid = int(pat.group(1))
					_viewer.voxel("value", vid, value)

	# ----------------------------------------------------------------------
	def _clearImage(self):
		self._imageName  = None
		self._imageTime  = 0
		self._imageBlack = 0
		self._imageWhite = 0
		self._imageData  = None

	# ----------------------------------------------------------------------
	def _imageLayer(self):
		_viewer = self.viewer._viewer

		self.viewer.lateDraw = True
		if self.getInt("Image.prompt"):
			_viewer.set("image", 2)	# 1=Late, 2=Prompt
		else:
			_viewer.set("image", 1)

		# Check for background image
		markers = GeometryLayers.getMarkers(
				self.project().geometry.layers[self.activeLayer].var)

		R,M = CalibrateImage.compute(markers)
		if M is None:
			self.flair().notify("Image Calibration",
				"Error in calibrating image",
				tkFlair.NOTIFY_ERROR)
			return

		iw = Vector(M[2][0], M[2][1], M[2][2])
		invmatrix = _viewer.matrix(-1)
		u  = Vector(invmatrix[0][0], invmatrix[0][1], invmatrix[0][2])
		v  = Vector(invmatrix[1][0], invmatrix[1][1], invmatrix[1][2])
		w  = Vector(invmatrix[2][0], invmatrix[2][1], invmatrix[2][2])
		dot = w * iw

		if abs(dot)<=0.9999:
			_viewer.viewOto0()

			# Not visible rotate the viewport
			if abs(dot)<0.01:
				u = iw.orthogonal()
			else:
				if dot<0.0: iw = -iw
				u = v^iw
			u.normalize()
			v = iw^u
			v.normalize()
			matrix = Matrix()
			matrix.make(u,v,iw)
			try:
				matrix.inverse()
				xyz = _viewer.origin()
				_viewer.matrix(matrix)
				_viewer.origin(xyz)
				self.viewer.draw(True,True)
			except:
				self.flair().write("\nViewport: %s\n"%(self.name))
				self.flair().write("Image matrix:\n")
				self.flair().write(str(M))
				self.flair().write("Viewport matrix:\n")
				self.flair().write(str(matrix))
				self.flair().notify("Invert Image matrix",
					"Singular image matrix found",
					tkFlair.NOTIFY_ERROR,
					"Output",
					self.flair().showOutput)

		alpha    = self.getInt("Image.alpha")
		black    = self.getInt("Image.black")
		white    = self.getInt("Image.white")
		filename = self.getVar("Image.file")

		load = False	# Flag whether to load image

		# Check if color range has changed
		if black!=self._imageBlack or \
		   white!=self._imageWhite:
			load = True
			self._imageBlack = black
			self._imageWhite = white

		# Check if file is different
		if not load:
			load = (filename != self._imageName)
			self._imageName = filename

		# Check if image has been modified
		if not load:
			try:	# Check time
				mtime = os.stat(filename)[ST_MTIME]
				load = (mtime > self._imageTime)
				self._imageTime = mtime
			except OSError:
				return

		if load:
			try:
				image = PIL.Image.open(filename)
				self._imageName = filename
			except:
				self.flair().notify("Error opening background image",
					"Cannot open file %r"%(filename),
					tkFlair.NOTIFY_ERROR)
				return

			image = image.convert("RGBX")
			image.load()
			try:
				self._imageData = image.tobytes()
			except AttributeError:
				self._imageData = image.tostring()
			size = image.size
			del image
			try:
				_viewer.image("data", self._imageData)
			except MemoryError:
				del image
				self.flair().notify("Not enough memory",
					sys.exc_info()[1],
					tkFlair.NOTIFY_ERROR)
				return
			_viewer.image("size", size)
			if black!=0 or white!=0xffffff:
				_viewer.image("level", (black,white))

		_viewer.image("alpha", alpha)
		_viewer.image("matrix", (R,M))

	# ----------------------------------------------------------------------
	def _usrbinLayer(self):
		_viewer = self.viewer._viewer

		self.viewer.lateDraw = True
		_viewer.set("usrbin", 1)

		for i,tag in enumerate(GeometryLayers.USRBINNAMES):
			if not self.getInt(tag):
				_viewer.usrbin(i, "cleanup")
				continue

			# FIXME remember and do not reload if image is not modified
			inp = self.getInt(tag+".input")
			try: det = self.getInt(tag+".det")
			except: det = 1

			palette = self.getVar(tag+".palette")
			try:
				_viewer.usrbin(i, "palette", GeometryLayers.PALETTENAMES.index(palette))
			except ValueError:
				palette = "Palette"
				_viewer.usrbin(i, "palette", 0)

			if inp:
				# Find usrbin from input
				try:
					usrbins = self.input().cardsSorted('USRBIN')
					card    = usrbins[det-1]
					case    = card.case()
					default = card.info.default[case]
					_viewer.usrbin(i, "checker",
						(card.intWhat(1),
						 card.numWhat(7), card.numWhat(4), card.intWhat(10, int(default[10])),
						 card.numWhat(8), card.numWhat(5), card.intWhat(11, int(default[11])),
						 card.numWhat(9), card.numWhat(6), card.intWhat(12, int(default[12]))))
					_viewer.usrbin(i, "offset", (0.0, 0.0, 0.0))
				except (KeyError, IndexError):
					continue
			else:
				# Load usrbin from file
				filename = self.getVar(tag+".file")
				if not filename:
					_viewer.usrbin(i, "cleanup")
					continue
				try:
					_viewer.usrbin(i, "load", filename, det)
				except IOError:
					self.flair().notify("Error loading USRBIN",
						sys.exc_info()[1],
						tkFlair.NOTIFY_ERROR)
					continue

				norm = self.getVar(tag+".norm")
				if not isinstance(norm,float):
					if norm and isinstance(norm,str):
						if norm[0]=="=": norm=norm[1:]
					try:
						norm = float(norm)
					except ValueError:
						try:
							norm = float(eval(str(norm), Input._globalDict))
						except:
							norm = 1.0
				else:
					norm = 1.0

				_viewer.usrbin(i, "norm", norm)
				_viewer.usrbin(i, "log",  self.getInt(palette+".log"))
				if _viewer.usrbin(i, "type") in (2,8,12,18):
					_viewer.set("fill",1)
					_viewer.usrbin(i, "regioncolor")

				try: xofs = float(self.getVar(tag+".xofs"))
				except: xofs = 0.0
				try: yofs = float(self.getVar(tag+".yofs"))
				except: yofs = 0.0
				try: zofs = float(self.getVar(tag+".zofs"))
				except: zofs = 0.0
				_viewer.usrbin(i, "offset", (xofs, yofs, zofs))

			rot = self.getVar(tag+".rotdefi")
			if rot is not None and rot!="":
				_viewer.usrbin(i, "matrix", self.input().getTransformation(rot))
			else:
				_viewer.usrbin(i, "matrix", None)
			rot = self.getVar(tag+".rotdefi2")
			if rot is not None and rot!="":
				_viewer.usrbin(i, "matrix", self.input().getTransformation(rot))

			_viewer.usrbin(i, "alpha", self.getInt(tag+".alpha"))

	# ----------------------------------------------------------------------
	def _3dLayer(self):
		self.viewer.lateDraw = True

		self.viewer._geometry.setLights()
		self.viewer._viewer.set("3D", 1)
		self.viewer._viewer.set("drawtime",       15)
		self.viewer._viewer.set("ambient",        self.getInt("3D.ambient"))
		self.viewer._viewer.set("antialias",      self.getInt("3D.antialias"))
		self.viewer._viewer.set("deflights",      self.getInt("3D.deflights"))
		self.viewer._viewer.set("edgedetect",     self.getInt("3D.edgedetect"))
		self.viewer._viewer.set("projection",     self.getInt("3D.projection"))
		self.viewer._viewer.set("reflections",    self.getInt("3D.reflections"))
		self.viewer._viewer.set("shadows",        self.getInt("3D.shadows"))
		self.viewer._viewer.set("skip1stblack",   self.getInt("3D.skip1stblack"))
		self.viewer._viewer.set("usrbinastexture",self.getInt("3D.usrbinastexture"))
		self.viewer._viewer.set("xray",           self.getInt("3D.xray"))
		self.viewer._viewer.camera("fov",         self.getFloat("3D.fov"))

		self.viewer._viewer.error("clear")

		self.viewer._viewer.set("clipbody",-1)	# Clear clipping bodies
		for i in range(3):
			try:
				clipbody = self.viewer._geometry.body(self.getVar("3D.clip.%d"%(i)),"id")
				self.viewer._viewer.set("clipbody", clipbody)
				if self.getInt("3D.clipnegative.%d"%(i)):
					self.viewer._viewer.set("clipnegative", i)
			except (KeyError, TypeError):
				pass

		self.viewer._viewer.set("projectbody",-1)	# Clear projection bodies
		for i in range(1000):
			try:
				projbody = self.viewer._geometry.body(self.getVar("3D.project.%d"%(i)),"id")
				self.viewer._viewer.set("projectbody", projbody)
			except (KeyError, TypeError):
				if i>3: break

	# ----------------------------------------------------------------------
	def _userdumpLayer(self):
		self.viewer.lateDraw = True
		self.viewer._viewer.set("userdump", 1)
		self.viewer._viewer.userdump("file",	self.getVar("Userdump.file"))
		self.viewer._viewer.userdump("start",	self.getInt("Userdump.start")-1)
		self.viewer._viewer.userdump("n",	self.getInt("Userdump.n"))

		for particle in Input.Particle.list:
			if not particle or particle[0]=='@': continue
			idx = Input.Particle.get(particle).id
			try:
				show = self.getInt("Userdump.%s"%(particle))
			except TypeError:
				show = False
			if show:
				self.viewer._viewer.userdump("show", idx)
				try:
					emin  = self.getFloat("Userdump.%s.emin"%(particle))
				except TypeError:
					pass
				try:
					emax  = self.getFloat("Userdump.%s.emin"%(particle))
				except TypeError:
					pass
				try:
					color = self.getInt("Userdump.%s.color"%(particle))
					self.viewer._viewer.userdump("color", idx, color)
				except TypeError:
					pass
			else:
				self.viewer._viewer.userdump("hide", idx)

	# ----------------------------------------------------------------------
	def _voxelLayer(self):
		for card in self.input()["VOXELS"]:
			if card.ignore() or card["@voxel"] is None: continue
			voxel = card["@voxel"]
			break
		else:
			return
		for roi in map(int,self.getStr('Voxel.roi').split()):
			self.viewer._viewer.voxel("roi", roi, 0)
		self.viewer._viewer.voxel("roialpha", self.getInt('Voxel.roialpha'))

	# ----------------------------------------------------------------------
	def _paletteLayer(self):
		_viewer = self.viewer._viewer
		_viewer.set("palette", 1)
		for idx,tag in enumerate(GeometryLayers.PALETTENAMES):
			_viewer.palette(idx, "reset")
			if not self.getInt(tag+"._visible"):
				continue
			cbmin = self.getFloat(tag+".min")
			cbmax = self.getFloat(tag+".max")
			cbn   = self.getInt(tag+".n")
			cblog = self.getInt(tag+".log")
			if cbn<=0: cbn = 32
			elif cbn>256: cbn=256
			palette = Palette.resizePalette(Palette.PALETTE[ \
					self.getVar(tag+".palette")], cbn)
			if cblog:
				try: cbmin = math.log10(cbmin)
				except: cbmin = -4.0
				try: cbmax = math.log10(cbmax)
				except: cbmax = 4.0

			_viewer.palette(idx, "show", self.getInt(tag))
			_viewer.palette(idx, "min", cbmin)
			_viewer.palette(idx, "max", cbmax)
			_viewer.palette(idx, "log", cblog)
			_viewer.palette(idx, "smooth",  self.getInt(tag+".smooth"))
			_viewer.palette(idx, "invert",  self.getInt(tag+".inv"))
			_viewer.palette(idx, "alphamin",self.getInt(tag+".alphamin"))
			_viewer.palette(idx, "alphamax",self.getInt(tag+".alphamax"))
			_viewer.palette(idx, "label",   self.getStr(tag+".label"))
			_viewer.palette(idx, "palette", palette)

	# ----------------------------------------------------------------------
	def rotateLayer(self):
		self.layerCombo.select(self.nextLayer())

	# ----------------------------------------------------------------------
	def nextLayer(self):
		return (self.layerCombo.select()+1)%self.layerCombo.size()

	# ----------------------------------------------------------------------
	def getOption(self, name):
		return GeometryLayers.option(self.project().geometry.var, name)

	# ----------------------------------------------------------------------
	def getVar(self, name, default=None):
		try:
			layer = self.project().geometry.layers[self.activeLayer].var
		except:
			return GeometryLayers.DEFLAYERS[name]

		if default is not None:
			return layer.setdefault(name, default)
		else:
			return GeometryLayers.option(layer, name)

	# ----------------------------------------------------------------------
	def getStr(self, name, default=None):
		return str(self.getVar(name,default))

	# ----------------------------------------------------------------------
	def getInt(self, name, default=None):
		return int(self.getVar(name,default))

	# ----------------------------------------------------------------------
	def getFloat(self, name, default=None):
		return float(self.getVar(name,default))

	# ----------------------------------------------------------------------
	def maximizeWindow(self, callLayout=True):
		if self.maxWindow.get():
			self.maxButton.config(text=Unicode.BLACK_DOWN_POINTING_TRIANGLE)
			if callLayout:
				self.parent.layout(self.name)
		else:
			self.maxButton.config(text=Unicode.BLACK_UP_POINTING_TRIANGLE)
			if callLayout:
				self.parent.layout("custom")

	# ----------------------------------------------------------------------
	def matrix(self, matrix=None):
		if matrix is None:
			return self.viewer._viewer.matrix()
		else:
			self.viewer._viewer.matrix(matrix)

	# ----------------------------------------------------------------------
	def export(self): self.viewer.export()
	def toNotes(self): self.viewer.toNotes()

	# ----------------------------------------------------------------------
	def historyButton(self, left, right):
		self._frames[0].historyButton(left, right)

	# ----------------------------------------------------------------------
	def zoomOnSelected(self, objects=None):
		self.viewer.zoomOnSelected(objects)

	# ----------------------------------------------------------------------
	def loadPlot(self, name=None):
		"""get projection information for flair PlotInfo"""
		if name is None: name = self.name
		for plot in self.project().plots:
			if plot.name==name and plot.type=="Geometry": break
		else:
			return

		_viewer = self.viewer._viewer
		var = plot.var
		try:
			u = Vector(list(map(float,var["basisu"].split())))
			v = Vector(list(map(float,var["basisv"].split())))
			w = u ^ v
			m = Matrix(4)
			m.make(u,v,w)
			_viewer.matrix(m.transpose())
		except KeyError: pass

		try:
			origin = list(map(float,var["origin"].split()))
			_viewer.origin(origin[0], origin[1], origin[2])
		except KeyError: pass

		try:
			extU,extV = list(map(float, var["extends"].split()))
			_viewer.extends(-2.0*extU, -2.0*extV, 2.0*extU, 2.0*extV)
			_viewer.zoom(2.0)
		except KeyError: pass

		if "layer" in var: self.layerCombo.set(var["layer"])

	# ----------------------------------------------------------------------
	def savePlot(self, name=None):
		"""set projection information to flair PlotInfo"""
		if name is None: name = self.name
		for plot in self.project().plots:
			if plot.name == name and plot.type == "Geometry": break
		else:
			plot = Project.PlotInfo(name)
			plot["title"] = "GeoViewer %s plot" % (name)	# for V2
			plot.format   = ".eps"
			plot.type     = "Geometry"
			self.project().plots.append(plot)

		_viewer = self.viewer._viewer
		var = plot.var

		var["origin"]  = "%g %g %g" % _viewer.origin()
		var["basisu"]  = "%g %g %g" % _viewer.basis(b'u')
		var["basisv"]  = "%g %g %g" % _viewer.basis(b'v')
		x1,y1,x2,y2 = _viewer.extends()
		z = _viewer.zoom()
		var["extends"] = "%g %g" % ((x2-x1)/2.0/z, (y2-y1)/2.0/z)

		u,v = _viewer.grid("axes")	# will return {-}[UXYZ]:{-}[VXYZ]
		u = u[-1]	# get rid of sign if any
		v = v[-1]
		if u=="U" or v=="V":
			var["coord"]   = "U-V"
		else:
			var["coord"]   = "%s-%s"%(u, v)
		var["boundaries"] = "1"
		try:	# when a second main window is opened and we close the
			# the window through the "protocol" I get a TclError
			dt = self.layerCombo.get()
			var["layer"]   = dt
			var["subtype"] = (dt=="Borders") and "Geometry" or "Material"
		except TclError:
			pass
		var["labels"]  = "1"

#===============================================================================
# Geometry4frame
#===============================================================================
class Geometry4Frame(Frame):
	# ----------------------------------------------------------------------
	def __init__(self, master, editor, **kw):
		Frame.__init__(self, master, **kw)
		self._geometry = geoviewer.Geometry()
		self._geometry.set("trackballsize", int(TRACKBALLSIZE*100))
		self._geometry.set("snapdistance",  SNAPDISTANCE)
		self._geometry.set("cursor",        CURSORSIZE)
		self._geometry.set("developer",     Input._developer)

		self.project        = None
		self.editor         = editor
		self.oldMouseAction = ACTION_SELECT
		self.action         = ACTION_SELECT
		self.mouseAction    = ACTION_SELECT	# Current mouse action when mouse is clicked
		self.actionTool     = None		# Tool dialog of action
		self.actionOption   = None		# option of mouse action e.g. what to add
		self.errors         = []
		self.warnings       = []
		self.voxels         = None		# voxels filename
		self.colorValue     = 0			# which color value to plot
		self.move3D         = False		# moving 2D or 3D
		self._regionColor   = {}
		self._autoBody      = None		# auto body dialog
		self.balloon        = False

		# Create frames
		self.vsplit = tkExtra.VSplitter(self,
					tkFlair.getFloat(tkFlair._PAGE_SECTION,
						"Geometry.vsplit", 0.5))
		self.vsplit.pack(side=TOP, expand=YES, fill=BOTH)

		self.tsplit = tkExtra.HSplitter(self.vsplit.topFrame(),
					tkFlair.getFloat(tkFlair._PAGE_SECTION,
						"Geometry.tsplit", 0.5))
		self.tsplit.pack(side=TOP, expand=YES, fill=BOTH)
		self.bsplit = tkExtra.HSplitter(self.vsplit.bottomFrame(),
					tkFlair.getFloat(tkFlair._PAGE_SECTION,
						"Geometry.bsplit", 0.5))
		self.bsplit.pack(side=BOTTOM, expand=YES, fill=BOTH)

		# The 4 geometry viewer frames
		view1 = GeometryFrame(self.tsplit.leftFrame(), self,
				VIEWNAME[0], takefocus=False)
		view1.pack(expand=YES, fill=BOTH)
		view2 = GeometryFrame(self.tsplit.rightFrame(), self,
				VIEWNAME[1], takefocus=False)
		view2.pack(expand=YES, fill=BOTH)
		view3 = GeometryFrame(self.bsplit.leftFrame(), self,
				VIEWNAME[2], takefocus=False)
		view3.pack(expand=YES, fill=BOTH)
		view4 = GeometryFrame(self.bsplit.rightFrame(), self,
				VIEWNAME[3], takefocus=False)
		view4.pack(expand=YES, fill=BOTH)

		# Save frames
		self.frames  = (view1, view2, view3, view4)
		self.viewers = [x.viewer._viewer for x in self.frames]

		self.resetViewports()

		for i,f in enumerate(self.frames):
			f.setColor(VIEWCOLOR[i])
			for j,f2 in enumerate(self.frames):
				if f is f2: continue
				f.viewer.addViewport(f2.viewer, VIEWCOLOR[j])

		self.custom = (0.5, 0.5, 0.5)

		self._layout = "custom"

		# Create tools
		self._zoneTool = GeometryTool.ZoneTool(self)
		self._zoneTool.hide()

	# ----------------------------------------------------------------------
	# Access
	def geometry(self):	return self._geometry
	def flair(self):	return self.editor.flair

	# ----------------------------------------------------------------------
	def __getitem__(self, item): return self.viewers[item]

	# ----------------------------------------------------------------------
	def resetViewports(self):
		self.frames[0].matrix(FRONT_MATRIX)
		self.frames[1].matrix(TOP_MATRIX)
		self.frames[2].matrix(LEFT_MATRIX)
		self.frames[3].matrix(BACK_MATRIX)

	# ----------------------------------------------------------------------
	def resetOrigin(self):
		for v in self.viewers:
			v.viewOto0()
			v.origin(0.0, 0.0, 0.0)

	# ----------------------------------------------------------------------
	def reloadOriginFrame(self):
		for frame in self.frames:
			frame.loadOrigin()

	# ----------------------------------------------------------------------
	def loadProjectColors(self):
		for card in self.project.input.cardsSorted("MATERIAL",True):
			# FIXME do I really need that?
			card[COLOR] = None
			self.materialColor(card)
			self.updateMaterial(card)

	# ----------------------------------------------------------------------
	def materialColor(self, card):
		if card is None: return _BACKGROUND_COLOR
		color = card[COLOR]
		if color is not None: return color
		mat = self.project.getMaterial(card.sdum())
		card[COLOR] = color = mat.color32()
		return color

	# ----------------------------------------------------------------------
	def materialAlpha(self, card):
		if card is None: return 0
		return int(self.project.getMaterial(card.sdum())["alpha"])

	# ----------------------------------------------------------------------
	def regionColor(self, card):
		if card is None: return _BACKGROUND_COLOR
		color = card[COLOR]
		if color is not None: return color
		try:
			card[COLOR] = self._regionColor[card.sdum()]
		except KeyError:
			card[COLOR] \
				= self._regionColor[card.sdum()] \
				= rndColor(len(self._regionColor))
		return card[COLOR]

	# ----------------------------------------------------------------------
	def setStatus(self, xyz):
		self.editor.xlabel.config(text=FMT%(xyz[0]))
		self.editor.ylabel.config(text=FMT%(xyz[1]))
		self.editor.zlabel.config(text=FMT%(xyz[2]))

	# ----------------------------------------------------------------------
	# Set status message
	# In order not to show the balloon set timeout=0
	# ----------------------------------------------------------------------
	def setStatusMsg(self, msg="", fg="DarkBlue", bg=None, timeout=4000):
		self.editor.setStatus(msg, fg, bg)
		if timeout and msg:
			self.setBalloon(msg, timeout)

	# ---------------------------------------------------------------------
	def setBalloon(self, msg, timeout=4000):
		self._help = msg
		tkExtra.Balloon.setWidget(self,
			self.winfo_pointerx(),
			self.winfo_pointery())
		tkExtra.Balloon.show()
		self.balloon = True
		if timeout: self.after(timeout, self.hideBalloon)

	# ----------------------------------------------------------------------
	def hideBalloon(self):
		self.balloon = False
		tkExtra.Balloon.hide()

	# ----------------------------------------------------------------------
	def updateSelection(self):
		if self.actionTool:
			self.actionTool.updateSelection()

	# ----------------------------------------------------------------------
	# Perform the mouse action in all frames
	# ----------------------------------------------------------------------
	def _setMouseAction(self, action=ACTION_SELECT):
		if self.mouseAction == action: return

		if self.mouseAction != ACTION_ADD_ZONE and self.mouseAction>=0:
			self.oldMouseAction = self.mouseAction
		self.mouseAction = action

		# Hide and Show action tool dialog
		if action>=0:
			if self.actionTool: self.actionTool.hide()
			if action == ACTION_ZONE:
				self.actionTool = self._zoneTool
			else:
				self.actionTool = None
			if self.actionTool: self.actionTool.show()

		cursor = mouseCursor(action)
#		import traceback; traceback.print_stack()
#		print "Action=",action, cursor
		for f in self.frames:
			f.viewer._cursor = cursor
			f.viewer.config(cursor=cursor)
			f.viewer._viewer.zone("clear")

		# Color level
		ll = self.geometry().set("lighterlevel")
		if action == ACTION_ADD_ZONE:
			newll = 200
		else:
			newll = 0
		if newll != ll:
			self.geometry().set("lighterlevel",newll)
			self.draw()

		# Update viewers if needed
		trackball = (action == ACTION_ORBIT)
		for v in self.viewers:
			v.set("trackball", trackball)
			v.expose()

		if self.actionTool and self.mouseAction>0: self.actionTool.nextAction()

	# ----------------------------------------------------------------------
	# Toggle action to new or previous
	# ----------------------------------------------------------------------
	def toggleAction(self, action):
		if self.mouseAction != action:
			self.editor.setMouseAction(action)
		else:
			self.editor.setMouseAction(self.oldMouseAction)

	# ----------------------------------------------------------------------
	# Repeat last action
	# ----------------------------------------------------------------------
	def repeatLastAction(self):
		self.editor.setMouseAction(self.oldMouseAction)

	# ----------------------------------------------------------------------
	def setLayer(self, dt=None):
		if dt is None:
			focus = self.focus_get()
			if focus.master in self.frames:
				dt = focus.master.nextLayer()
			else:
				dt = self.frames[0].nextLayer()
		for frame in self.frames:
			frame.layerCombo.select(dt)

	# ----------------------------------------------------------------------
	def showLayer(self, layer):
		self.editor.showLayer(layer)

	# ----------------------------------------------------------------------
	def rotateLayout(self):
		n = VIEWORDER.index(self._layout) + 1
		if n>=len(VIEWORDER): n = 0
		self.layout(VIEWORDER[n])

	# ----------------------------------------------------------------------
	def customLayout(self):
		focus = self.focus_get()
		if focus.master in self.frames:
			fmt = focus.master.name
		else:
			fmt = ""
		if self._layout==fmt:
			self.layout("custom")
		else:
			self.layout(fmt)

	# ----------------------------------------------------------------------
	def compareLayout(self, v, t, b):
		if abs(self.vsplit.split-v) > 0.02: return False
		if v < 0.02:
			# Compare only the bottom
			return abs(self.bsplit.split-b) < 0.02
		elif v > 0.98:
			# Compare only the top
			return abs(self.tsplit.split-t) < 0.02
		else:
			# Compare both top and bottom
			return	abs(self.bsplit.split-b) < 0.02 and \
				abs(self.tsplit.split-t) < 0.02

	# ----------------------------------------------------------------------
	def layout(self, fmt):
		v,t,b = VIEWFORMAT.get(fmt, self.custom)
		if self.compareLayout(v,t,b): return	# same

		self._layout = fmt
		self.custom = ( self.vsplit.split,
				self.tsplit.split,
				self.bsplit.split)

		self.vsplit.setSplit(v)
		self.tsplit.setSplit(t)
		self.bsplit.setSplit(b)

		for f in self.frames: f.maxWindow.set(False)
		if v>0.98:
			if t>0.98:	# top left is maxized
				self.frames[0].maxWindow.set(True)
			elif t<0.02:	# top right is maximized
				self.frames[1].maxWindow.set(True)
		elif v<0.02:
			if b>0.98:	# bottom left is maxized
				self.frames[2].maxWindow.set(True)
			elif b<0.02:	# bottom right is maximized
				self.frames[3].maxWindow.set(True)
		for f in self.frames: f.maximizeWindow(False)

	# ----------------------------------------------------------------------
	def saveState(self):
		"""Save views as flair plots"""
		if self.project is None: return
		for frame in self.frames:
			frame.savePlot()

		self.setOption("Frame.vsplit", self.vsplit.split)
		self.setOption("Frame.tsplit", self.tsplit.split)
		self.setOption("Frame.bsplit", self.bsplit.split)

	# ----------------------------------------------------------------------
	def loadState(self):
		"""Load views from flair plots"""
		for frame in self.frames:
			frame.loadPlot()

	# ----------------------------------------------------------------------
	def log(self, msg):
		# FIXME Maybe replace with log() from log.py
		if not self.flair(): return
		if not self.flair().log: return
		self.flair().log(msg)

	# ----------------------------------------------------------------------
	def addError(self, msg):
		# FIXME Maybe replace with log() from log.py
		self.log(msg+"\n")
		self.errors.append(msg+"\n")

	# ----------------------------------------------------------------------
	def addWarning(self, msg):
		# FIXME Maybe replace with log() from log.py
		self.log(msg+"\n")
		self.warnings.append(msg+"\n")

	# ----------------------------------------------------------------------
	# set a variable in all views
	# ----------------------------------------------------------------------
	def set(self, var, value=None):
		if value is None:
			return self._geometry.set(var)
		else:
			self._geometry.set(var, value)

	# ----------------------------------------------------------------------
	# Add material in geometry
	# ----------------------------------------------------------------------
	def addMaterial(self, card):
		try:
			return self._geometry.addMaterial(card.name())
		except:
			self.addError("Material %r %s"%(card.name(),sys.exc_info()[1]))
			return None

	# ----------------------------------------------------------------------
	def updateMaterial(self, card):
		mid = card[ID]
		if mid is None: return
		self._geometry.material(mid, "Z",       card.intWhat(1))
		self._geometry.material(mid, "weight",  card.numWhat(2))
		self._geometry.material(mid, "density", card.numWhat(3))
		self._geometry.material(mid, "A",       card.intWhat(6))

		mat = self.project.getMaterial(card.name())
		self._geometry.material(mid, "specular", float(mat["specular"]))
		self._geometry.material(mid, "shine",    float(mat["shine"]))
		self._geometry.material(mid, "fuzz",     float(mat["fuzz"]))
		self._geometry.material(mid, "ior",      float(mat["ior"]))
#		self._geometry.material(mid, "alpha",    int(mat["alpha"]))

	# ----------------------------------------------------------------------
	# add a body in all views
	# ----------------------------------------------------------------------
	def addBody(self, name, tag, what, matrix=None):
		if not Input._NAMEPAT.match(name):
			self.addError("Body %r invalid name"%(name))
		try:
			bid = self._geometry.addBody(name, tag)
			err = self._geometry.body(bid, "what", what)
			if err:
				self.addError("Body %r %s"%(name,err))
			if matrix is not None:
				self._geometry.body(bid, "matrix", matrix)
			return bid
		except:
			self.addError("Body %r %s"%(name,sys.exc_info()[1]))
			return None

	# ----------------------------------------------------------------------
	# modify body in all views
	# ----------------------------------------------------------------------
	def body(self, idx, var, value=None):
		if idx is None: return
		if value is None:
			return self._geometry.body(idx, var)
		else:
			self._geometry.body(idx, var, value)

	# ----------------------------------------------------------------------
	# add a region in all views
	# ----------------------------------------------------------------------
	def addRegion(self, name, expr):
		if not Input._NAMEPAT.match(name):
			self.addError("Region %r invalid name"%(name))

		rid = None
		try:
			rid = self._geometry.addRegion(name)
			self._geometry.region(rid, "expr", expr)
		except:
			self.addError("Region %r %s"%(name,sys.exc_info()[1]))
		return rid

	# ----------------------------------------------------------------------
	# modify region in all views
	# ----------------------------------------------------------------------
	def region(self, idx, var, value=None):
		if idx is None: return
		if value is None:
			return self._geometry.region(idx, var)
		else:
			self._geometry.region(idx, var, value)

	# ----------------------------------------------------------------------
	# Manipulate zone selection
	# ----------------------------------------------------------------------
	def zone(self, cmd, value=None):
		if value is None:
			for view in self.viewers:
				view.zone(cmd)
		else:
			for view in self.viewers:
				# Ignore invalid bodies from zones
				try: view.zone(cmd, value)
				except ValueError: pass

	# ----------------------------------------------------------------------
	# Add object
	# ----------------------------------------------------------------------
	def addObject(self, card):
		oid = card[ID] = self._geometry.addObject(card.sdum(), card.tag)
		self._geometry.object(oid, "x", card.numWhat(1))
		self._geometry.object(oid, "y", card.numWhat(2))
		self._geometry.object(oid, "z", card.numWhat(3))
		self._geometry.object(oid, "option", card.intWhat(4))
		s,a = divmod(card.intWhat(5),100)
		self._geometry.object(oid, "size",   s)
		self._geometry.object(oid, "anchor", a)
		self._geometry.object(oid, "color", card.intWhat(6))

		# Additional variables
		if   card.tag == "!arrow":
			self._geometry.object(oid, "dx", card.numWhat(7))
			self._geometry.object(oid, "dy", card.numWhat(8))
			self._geometry.object(oid, "dz", card.numWhat(9))
		elif card.tag == "!light":
			self._geometry.object(oid, "dx",      card.numWhat(7))
			self._geometry.object(oid, "dy",      card.numWhat(8))
			self._geometry.object(oid, "dz",      card.numWhat(9))
			self._geometry.object(oid, "power",   card.numWhat(10))
			self._geometry.object(oid, "falloff", card.intWhat(11))
			self._geometry.object(oid, "specular",card.numWhat(12))
		elif card.tag == "!ruler":
			self._geometry.object(oid, "xe", card.numWhat(7))
			self._geometry.object(oid, "ye", card.numWhat(8))
			self._geometry.object(oid, "ze", card.numWhat(9))
			self._geometry.object(oid, "xa", card.numWhat(10))
			self._geometry.object(oid, "ya", card.numWhat(11))
			self._geometry.object(oid, "za", card.numWhat(12))
		elif card.tag == "!mesh":
			self._geometry.object(oid, "file", card.extra())
		elif card.tag == "!spline":
			n = (card.nwhats()-7)//3
			for i in range(n):
				w = card.what(i*3+7)
				if w=="": break
				self._geometry.object(oid, "x%d"%(i+1), card.numWhat(i*3+7))
				self._geometry.object(oid, "y%d"%(i+1), card.numWhat(i*3+8))
				self._geometry.object(oid, "z%d"%(i+1), card.numWhat(i*3+9))

		return oid

	# ----------------------------------------------------------------------
	def object(self, idx, var, value=None):
		if value is None:
			return self._geometry.object(idx, var)
		else:
			self._geometry.object(idx, var, value)

	# ----------------------------------------------------------------------
	# Add a rotdefi as object to geometry
	# ----------------------------------------------------------------------
	def addRotdefi(self, card, name, matrix):
		oid = card[ID] = self._geometry.addObject(name, "!rotdefi")
		self._geometry.object(oid, "axissize", 40)
		self._geometry.object(oid, "linewidth", 2)
		self._geometry.object(oid, "axiswidth", 3)
		self._geometry.object(oid, "size", 7)
		x = card["x"]
		if x: self._geometry.object(oid, "x", float(x))
		y = card["y"]
		if y: self._geometry.object(oid, "y", float(y))
		z = card["z"]
		if z: self._geometry.object(oid, "z", float(z))
		self._geometry.object(oid, "matrix", matrix)
		return oid

	# ----------------------------------------------------------------------
	# Add a BEAM card as object to geometry
	# ----------------------------------------------------------------------
	def addBeam(self, card, name):
		oid = card[ID] = self._geometry.addObject(name, "!beam")
		if card.tag == "BEAM":
			ew = 1
			dw = 3
		else:
			ew = 2
			dw = 5
		self._geometry.object(oid, "energy",     abs(card.numWhat(ew)))
		try:
			self._geometry.object(oid, "scale",      abs(float(card.get("scale",1.0))))
		except ValueError:
			pass
		self._geometry.object(oid, "divergence", card.numWhat(dw)/1000.0)
		return oid

	# ----------------------------------------------------------------------
	def setBeamPos(self, oid, card):
		sdum = card.sdum().strip()
		if sdum=="" or sdum=="NEGATIVE":
			card[ID] = oid
			self._geometry.object(oid, "x", card.numWhat(1))
			self._geometry.object(oid, "y", card.numWhat(2))
			self._geometry.object(oid, "z", card.numWhat(3))
			dx = card.numWhat(4)
			dy = card.numWhat(5)
			dz = sqrt(1.0 - dx*dx - dy*dy)
			if sdum=="NEGATIVE": dz = -dz
			self._geometry.object(oid, "dx", dx)
			self._geometry.object(oid, "dy", dy)
			self._geometry.object(oid, "dz", dz)
			return True
		else:
			return False

	# ----------------------------------------------------------------------
	# Update the edit handle of the geometry
	# ----------------------------------------------------------------------
	def updateEdit(self, card=None):
		# Don't show handler within an action
		# FIXME needs better thinking
		if self.mouseAction<0:
			self._geometry.cursor("clear")
			return

		# Check last card
		if card and card.get(SELECT,0) & BIT_SELECT:
			if card.type() == Input.Card.BODY:
				self._geometry.cursor("body", card[ID])
		else:
			lastBody = lastRegion = lastObject = None
			for card in self.project.input.allcards():
				if card.ignore(): continue
				if card.get(SELECT,0) & BIT_SELECT == 0: continue
				t = card.type()
				if t == Input.Card.BODY:
					lastBody = card
				elif t == Input.Card.REGION:
					lastRegion = card
				elif t == Input.Card.OBJECT:
					lastObject = card

			# Highest priority bodies
			if lastBody is not None:
				self._geometry.cursor("body", lastBody[ID])

			elif lastRegion is not None:
				try:
					self._geometry.cursor("region", lastRegion[ID])
				except IndexError:
					# Voxels will fail
					pass

#			elif lastObject is not None:
#				self._geometry.cursor("object", lastObject[ID])

			else:
				self._geometry.cursor("clear")

	# ----------------------------------------------------------------------
	# Update region and body properties in geometry
	# ----------------------------------------------------------------------
	def updateGeometry(self):
		rotdefi = self.project.input.getTransformation()
		for region in self.project.input["REGION"]:
			if region.ignore(): continue
			if region.sdum()=="&": continue
			try:
				rid = region[ID]
				self._geometry.region(rid, "type", region["@type"])
				if region["@rotdefi"] is not None:
					rot = signedRotdefi(rotdefi, region["@rotdefi"])
					self._geometry.region(rid, "matrix", self.project.input.getTransformation(rot,rotdefi))
			except (KeyError, TypeError):
				pass

	# ----------------------------------------------------------------------
	# Load project
	# ----------------------------------------------------------------------
	def loadProject(self, project=None, run=None):
		self.flair().logReset()
		self.setStatusMsg()
		if _PROFILE:
			say("\n*** GeometryViewer.loadProject started")
			start = time.time()
			#traceback.print_stack()

		if project: self.project = project
		self.checkLayers()

		# Remember cursor
		cursor     = self._geometry.cursor("cursor")
		cursorType = self._geometry.cursor("type")
		cursorId   = self._geometry.cursor("id")

		# Stop threads
		self.stopThreads()

		# Cleanup everything
		self._geometry.cleanup()

		# Load the new project
		self._geometry.lock()
		try:	self._loadInput(run)
		except: traceback.print_exc()
		self._geometry.unlock()

		# Restore selection cursor if any
		if cursor:
			try: self._geometry.cursor(cursorType.decode('utf-8'), cursorId)
			except IndexError: pass

		# Derive...
		for viewer in self.viewers:
			viewer.derive()

		# Load information
		for frame in self.frames:
			frame.loadOrigin()

		# Check for possible errors
		if self.errors or self.warnings:
			self.setStatusMsg("Errors found while loading geometry.", "Red")

		if _PROFILE:
			say("*** GeometryViewer.loadProject ended", time.time()-start,"s")

		self.flair().errorsDuringLoading()

	# ----------------------------------------------------------------------
	# Load Input file to geoviewer
	# ----------------------------------------------------------------------
	def _loadInput(self, run):
		if _PROFILE: start = time.time()
		del self.errors[:]
		del self.warnings[:]

		# prepare geometry
		self.voxels = None

		input = self.project.input
		if run:
			defines = run.defines
		else:
			defines = None
		try:
			input.preprocess(defines)
		except:
			self.addError("Preprocessing file")

		#input.validate()
		#lap=time.time(); say("GeometryViewer._loadInput:validate=",lap-start); start=lap
#		self.palette.loadInput(input)
		if _PROFILE:
			lap=time.time(); say("GeometryViewer._loadInput:paletteInput=",lap-start); start=lap

		# Scan input for body and region properties
		input.bodyProperties()
		if _PROFILE:
			lap=time.time(); say("GeometryViewer._loadInput:bodyProperties=",lap-start); start=lap
		regionList = input.regionProperties()
		if _PROFILE:
			lap=time.time(); say("GeometryViewer._loadInput:regionProperties=",lap-start); start=lap

		# transformations
		rotdefi = input.getTransformation()
		for key,matrix in list(rotdefi.items()):
			if isinstance(key,int):
				try: self._geometry.rotdefi(key, matrix)
				except: pass

		# --------------------
		#       Objects
		# --------------------
		cards = []
		for tag in Input.FLAIR_TAGS:
			cards.extend(input[tag])
		cards.sort(key=attrgetter("_pos"))

		for card in cards:
			if card.ignore(): continue
			self.addObject(card)

		# --------------------
		#       BEAM
		# --------------------
		for card in input["BEAM"]:
			if card.ignore(): continue
			oid = self.addBeam(card, "beam")
			break
		else:
			oid = None
		if oid is not None:
			for card in input["BEAMPOS"]:
				if card.ignore(): continue
				if self.setBeamPos(oid, card): break

		# --------------------
		#  All SPOTxxx cards
		# --------------------
		if "SPOTBEAM" in input.cards:
			nspot = 0
			spots = {}
			for card in input.cardlist:
				if card.ignore(): continue
				if card.tag == "SPOTBEAM":
					oid = self.addBeam(card, "spot")
					if card.intWhat(0)==0:
						nspot += 1
					else:
						nspot = card.intWhat(0)
					spots[nspot] = oid
				elif card.tag == "SPOTPOS":
					if card.intWhat(0)!=0:
						try:
							oid = spots[card.intWhat(0)]
						except:
							#print "ERROR mixup in spots id"
							pass
					card[ID] = oid
					self._geometry.object(oid, "x", card.numWhat(1))
					self._geometry.object(oid, "y", card.numWhat(2))
					self._geometry.object(oid, "z", card.numWhat(3))
				elif card.tag == "SPOTDIR":
					if card.intWhat(0)!=0:
						try:
							oid = spots[card.intWhat(0)]
						except:
							#print "ERROR mixup in spots id"
							pass
					card[ID] = oid
					self._geometry.object(oid, "dx", card.numWhat(1))
					self._geometry.object(oid, "dy", card.numWhat(2))
					self._geometry.object(oid, "dz", card.numWhat(3))

		# --------------------
		#       ROT-DEFI
		# --------------------
		rot = {}
		for card in input["ROT-DEFI"]:
			if card.ignore(): continue
			name = card.sdum().strip()
			#if not name: ... find name from what(1)
			if name in rot: continue
			matrix = input.getTransformation(name, rotdefi)
			if matrix: self.addRotdefi(card, name, matrix)

		self._geometry.setLights()

		if _PROFILE:
			lap=time.time(); say("GeometryViewer._loadInput:transformations=",lap-start); start=lap

		# --------------------
		#      MATERIAL
		# -------------------#-
		for card in input.cardsSorted("MATERIAL",True):
			card[ID] = self.addMaterial(card)
			self.updateMaterial(card)

		# --------------------
		#       Bodies
		# --------------------
		rid = -1
		rname = ""
		matrix = [Matrix(4,type=1)]	# list of transformation matrices
		istranslat=0
		for card in input.cardlist:
			if not card.isGeo(): continue
			if card.ignore(): continue
			if len(card.tag)==3:
				if card.tag=="END": continue
				if card.invalid:
					self.addError("Body %r %s"%(card.sdum(),
						card.errorMessage()))

				try:	# We have to check if body already exists
					self._geometry.body(card.sdum(),"id")
					self.addError("Body %r Duplicate body"%(card.sdum()))
					#card.addInvalid("Duplicate body")
				except KeyError:
					pass
				except:
					say(sys.exc_info()[0])
					say(card)
					say(card.sdum(),type(card.sdum()))

				try:
					whats = [card.numWhat(i) for i in range(1,card.nwhats()) ]
					#whats = map(float,card.whats()[1:])
				except ValueError:
					self.addError("Body %r invalid whats"%(card.sdum()))
					continue

				# Add the body
				if len(matrix)==1:
					card[ID] = self.addBody(card.sdum(), card.tag, whats)
				else:
					card[ID] = self.addBody(card.sdum(), card.tag, whats, matrix[-1])
				bbox = card["bbox"]
				if bbox:
					try:
						self._geometry.body(card[ID], "bbox", list(map(float,bbox.split())))
					except ValueError:
						say("ERROR: Invalid number in body's bbox")
						say(card)
						say(card.sdum(),type(card.sdum()))

			elif card.tag[0]=="$":
				if card.tag=="$start_translat":
					istranslat=1
					matrix.append(matrix[-1]*Matrix.translate(
						card.numWhat(1), card.numWhat(2), card.numWhat(3)))
#					print ("transl",istranslat,len(matrix))
				elif card.tag=="$start_transform":
					# FIXME accept transformation ONLY if matrix is identity
					# otherwise issue an error
					if istranslat ==1 :
						self.addError("$start_transform must be before $start_translat\n"+str(card))
					# Fluka way applying $start_translat before $start_transform
						matrix.append(input.getTransformation(card.evalWhat(1),rotdefi)*matrix[-1])
						istranslat = -1
#						print ("transf",istranslat,len(matrix))
					elif istranslat == -1 :
						self.addError("$start_transform and $start_translat messed up \n"+str(card))
						self.addError("please correct, no transformation is applied \n")
						del matrix[-1]

					# normal way
					elif istranslat == 0 :

						matrix.append(matrix[-1]*input.getTransformation(card.evalWhat(1),rotdefi))

				elif card.tag=="$start_expansion":
					matrix.append(matrix[-1]*Matrix.scale(card.numWhat(1)))
				elif card.tag[:5]=="$end_":
					istranslat =0 
#					print ("end",len(matrix))
					if len(matrix)>1:
						del matrix[-1]
#						print ("end-dele",len(matrix))
					else:
						self.addError("Too many $end_xxx cards")

			elif card.tag == "VOXELS":
				if self.voxels is not None:
					self.addError("More than one voxel cards found")
					continue

				self.voxels = "%s.vxl"%(card.sdum().strip())
#				print ("voxels", type(self.voxels), self.voxels)
#				pdb.set_trace()
				try:
					self._geometry.voxel("load", self.voxels)
				except:
					self.addError("Cannot open voxel file %s"%(self.voxels))
					self.voxels = None
					continue

				try:
					x = card.numWhat(1)
					y = card.numWhat(2)
					z = card.numWhat(3)
					rot = card.what(4)

					self._geometry.voxel("xlow", x)
					self._geometry.voxel("ylow", y)
					self._geometry.voxel("zlow", z)

					card[ID] = self.addBody("VOXEL", "RPP",
						[ x, self._geometry.voxel("xhigh"),
						  y, self._geometry.voxel("yhigh"),
						  z, self._geometry.voxel("zhigh") ])

					if rot:
						irot = input.getTransformation(rot,rotdefi,True)
						# Unfortunately bodies and regions/voxels
						# use opposite transformation
						self._geometry.voxel("matrix",irot) # None will clear the matrix
						if irot is not None:
							self._geometry.body(card[ID], "matrix",
								input.getTransformation(rot,rotdefi))
					else:
						self._geometry.voxel("matrix",None)
						self._geometry.body(card[ID], "matrix", None)
				except:
					self.addError("loading voxel %s"%(sys.exc_info()[1]))

		# Cleanup maybe old voxel
		if self.voxels is None:
			self._geometry.voxel("free")

		if len(matrix)!=1:
			self.addError("Unbalanced body transformations $start .. $end")

		if _PROFILE:
			lap=time.time(); say("GeometryViewer._loadInput:bodies=",lap-start); start=lap

		lattice_error = False
		# --------------------
		#       Regions
		# --------------------
		def _addRegion():
			# add last region
			if region.invalid:
				self.addError("Region %r  %s"%(region.sdum(),
					region.errorMessage()))
			try:	# We have to check if region already exists
				self._geometry.region(region.sdum(),"id")
				self.addError("Duplicate region %r"%(region.sdum()))
			except KeyError:
				pass
			except:
				say(sys.exc_info()[0])
				say(card)
				say(card.sdum(),type(card.sdum()))

			if region["@materialDefined"] is None:
				self.addWarning("Region %r not assigned any material"\
						%(region.sdum()))
			#if region["@type"] == Input.REGION_VOXEL: return
			try:
				expList = csg.splitZones(csg.tokenize(expr))
				if len(expList)==0 and region["@type"] != Input.REGION_VOXEL:
					self.addWarning("Region %r empty expression"%(region.sdum()))
					expList = []
				else:
					region["@zone"] = list(map(csg.toString, expList))
					expList = regionExp(expList)
			except:
				self.addError("Region %r %s"%\
					(region.sdum(),sys.exc_info()[1]))
				expList = []

			rid = self.addRegion(region.name(), expList)
			if rid is not None:
				region[ID] = rid
				self.region(rid, "type", region["@type"])
				if region["@rotdefi"] is not None  or region["@rotdefi2"] is not None:
					rot = signedRotdefi(rotdefi, region["@rotdefi"])
					rot2 = signedRotdefi(rotdefi, region["@rotdefi2"])

					if rot==0 and rot2==0:
						self.addWarning("Cannot display LATTICE region %r, no ROTDEFI is assigned, therefore it requires an external lattic.f routine"\
								%(region.sdum()))
					if rot is not None:
						self.region(rid, "matrix", input.getTransformation(rot,  rotdefi))
					if rot2 is not None:
						self.region(rid, "matrix", input.getTransformation(rot2, rotdefi))
#					self.region(rid, "rotdefi", rot)
				mat = region["@material"]
				if mat and mat[ID]:
					self._geometry.region(rid, "material", mat[ID])

		regionCards = input.cardsSorted("REGION")
		if regionCards:
			expr   = None
			region = None
			color  = None

			i = 0
			n = len(regionCards)
			while i < n:
				card = regionCards[i]
				i += 1
				if card.ignore(): continue

				if card.name() == "&" and expr is not None: # append to previous expression
					expr += card.extra()
					# FIXME should moved to regionProperties
					region["@cont"] = True
					continue

				if expr is not None: _addRegion()

				# move old card to region and process new card
				region = card
				# create a region color if needed
				self.regionColor(card)
				# use by default the materialColor, we change it latter
				color = self.materialColor(card["@material"])
				expr = card.extra()

			if expr is not None: _addRegion()

			#if self.voxels:
			#	rid = self.addRegion("VOXEL", [("STD", ["VOXEL","$"])])
			#	if rid is not None:
			#		self.region(rid, "type", Input.REGION_VOXEL)

		if _PROFILE:
			lap=time.time(); say("GeometryViewer._loadInput:regions=",lap-start); start=lap
		# --------------------
		#       Usrbins
		# --------------------
		#idx = 0
		#for card in input["USRBIN"]:
		#	if card.ignore(): continue
		#	idx += 1
		#	name = "USRBIN.%d"%(idx)

		#	if card.intWhat(1) in (0,10):
		#		tag   = "RPP"
		#		whats = [card.numWhat(7), card.numWhat(4),
		#		         card.numWhat(8), card.numWhat(5),
		#		         card.numWhat(9), card.numWhat(6)]
		#		self.addBody(name, tag, whats)

		#	elif card.intWhat(1) in (1,11):
		#		tag   = "RCC"
		#		whats = [card.numWhat(7), card.numWhat(4),
		#		         card.numWhat(8), card.numWhat(5),
		#		         card.numWhat(9), card.numWhat(6)]
		#		self.addBody(name, tag, whats)

		#	else:
		#		continue

	# ----------------------------------------------------------------------
	# Display message if objects are locked
	# ----------------------------------------------------------------------
	def check4Frozen(self):
		frozen = False
		for card in self.project.input.cardlist:
			if card.ignore(): continue
			if card.get(SELECT,0) & BIT_FREEZE:
				frozen = True
				break
		if frozen:
			self._geometry.message("Objects frozen", _FREEZE_COLOR)
		else:
			self._geometry.message("")
		self.draw()

	# ----------------------------------------------------------------------
	def viewportProject(self, asthread=False, all=False):
		"""project all viewports"""
		for view in self.viewers:
			view.project(asthread, all)

	# ----------------------------------------------------------------------
	def derive(self):
		"""Derive generation in all viewers"""
		for view in self.viewers:
			view.derive()

	# ----------------------------------------------------------------------
	def expose(self, event=None):
		"""Expose all viewports"""
		for view in self.viewers:
			view.expose()

	# ----------------------------------------------------------------------
	def stopThreads(self):
		"""Stop all threads"""
		for frame in self.frames:
			frame.viewer.stop()

	# ----------------------------------------------------------------------
	def reload(self):
		self.editor.reload()

	# ----------------------------------------------------------------------
	def draw(self):
		"""Redraw screen"""
		for view in self.frames:
			view.draw()

	# ----------------------------------------------------------------------
	def drawFast3D(self):
		for view in self.viewers:
			if view.get("3D"):
				view.draw(False, DRAW_FAST)

	# ----------------------------------------------------------------------
	def redraw(self):
		"""Force regeneration of projections & plots"""
		# Update variables (just in case)
		self.set("backgroundcolor",	_BACKGROUND_COLOR)
		self.set("selectcolor",		_SELECT_COLOR)
		self.set("visiblecolor",	_VISIBLE_COLOR)
		self.set("regioncolor",		_REGION_COLOR)
		self.set("zonecolor",		_ZONE_COLOR)
		self.set("errorcolor",		_ERROR_COLOR)
		self.set("wireframecolor",	_WIREFRAME_COLOR)
		self.set("bodybboxincolor",	_BODY_BBIN_COLOR)
		self.set("bodybboxoutcolor",	_BODY_BBOUT_COLOR)
		self.set("zonebboxcolor",	_ZONE_BB_COLOR)
		self.set("regionbboxcolor",	_REGION_BB_COLOR)
		self.set("titlecolor",		_TITLE_COLOR)
		self.set("gridtextcolor",	_GRIDTEXT_COLOR)
		self.set("labelcolor",		_LABEL_COLOR)
		self.set("latticecolor",	_LATTICE_COLOR)

		# refresh viewports
		for frame in self.frames:
			frame.viewer.set("backgroundcolor", _BACKGROUND_COLOR)
			frame.redraw()

	# ----------------------------------------------------------------------
	def checkLayers(self):
		# Create list of layers
		for frame in self.frames:
			frame.layerCombo.fill(self.project.geometry.layersList())

	# ----------------------------------------------------------------------
	# Update all viewports with the new layer configuration
	# ----------------------------------------------------------------------
	def updateLayers(self, layer):
		self.checkLayers()
		for frame in self.frames:
			if layer == frame.layerCombo.get():
				frame.viewer.stop()
				frame.layer()
				frame.viewer.drawFast3D()
				frame.draw()

	# ----------------------------------------------------------------------
	def autoBodyDialog(self, event=None):
		if self._autoBody is None:
			self._autoBody = GeometryBody.GeometryAutoBody(self, self.editor, self.project)
		else:
			try:
				self._autoBody.deiconify()
				self._autoBody.lift()
			except:
				self._autoBody = GeometryBody.GeometryAutoBody(self,
					self.editor, self.project)

	# ----------------------------------------------------------------------
	def getOption(self, name, default=None):
		return GeometryLayers.option(self.project.geometry, name, default)

	# ----------------------------------------------------------------------
	def setOption(self, name, value):
		self.project.geometry.var[name] = value

	# ----------------------------------------------------------------------
	def getVar(self, name, default=None):
		layer = self.project.geometry.layers[0].var
		return GeometryLayers.option(layer, name, default)

	# ----------------------------------------------------------------------
	def setVar(self, name, value):
		self.project.geometry.layers[0].var[name] = value

	# ----------------------------------------------------------------------
	# Add menus
	# ----------------------------------------------------------------------
	def _addObjectMenu(self, menu):
		menu.add_command(label="Point", underline=0, command=self.newPoint,
				compound=LEFT, image=tkFlair.icons["!point"],
				accelerator="p")
		menu.add_command(label="Arrow", underline=0, command=self.newArrow,
				compound=LEFT, image=tkFlair.icons["!arrow"])
		#if Input._developer:
		menu.add_command(label="Mesh", underline=0, command=self.newMesh,
				compound=LEFT, image=tkFlair.icons["!mesh"])
		menu.add_command(label="Ruler", underline=0, command=self.newRuler,
				compound=LEFT, image=tkFlair.icons["!ruler"])
		#menu.add_command(label="Rotdefi", underline=0, command=self.newRotdefi,
		#		compound=LEFT, image=tkFlair.icons["!axes"])
		menu.add_command(label="Light", underline=0, command=self.newLight,
				compound=LEFT, image=tkFlair.icons["!light"])

	# ----------------------------------------------------------------------
	def _visibilityMenu(self, menu):
		submenu = Menu(menu)
		menu.add_cascade(label="Visibility", menu=submenu, underline=0,
				compound=LEFT, image=tkFlair.icons["view"])
		submenu.add_command(label="Set",     underline=0,
				compound=LEFT, image=tkFlair.icons["enable"],
				command=self.editor.setVisibility)
		submenu.add_command(label="Clear",   underline=0,
				compound=LEFT, image=tkFlair.icons["disable"],
				command=self.editor.unsetVisibility)
		submenu.add_command(label="Toggle",  underline=0,
				compound=LEFT, image=tkFlair.icons["toggle"],
				accelerator="v",
				command=self.editor.toggleVisibility)

	# ----------------------------------------------------------------------
	def _wireframeMenu(self, menu):
		submenu = Menu(menu)
		menu.add_cascade(label="Wireframe", menu=submenu, underline=0,
				compound=LEFT, image=tkFlair.icons["wireframe"])
		submenu.add_command(label="Show",     underline=0,
				compound=LEFT, image=tkFlair.icons["wireframe"],
				command=self.editor.set3D)
		submenu.add_command(label="Hide",   underline=0,
				compound=LEFT, image=tkFlair.icons["wireframe_off"],
				command=self.editor.unset3D)
		submenu.add_command(label="Toggle",  underline=0,
				compound=LEFT, image=tkFlair.icons["toggle"],
				accelerator="#",
				command=self.editor.toggle3D)

	# ----------------------------------------------------------------------
	def _freezeMenu(self, menu):
		submenu = Menu(menu)
		menu.add_cascade(label="Selection", menu=submenu, underline=0,
				compound=LEFT, image=tkFlair.icons["freeze"])
		submenu.add_command(label="Freeze",     underline=0,
				compound=LEFT, image=tkFlair.icons["freeze"],
				command=self.editor.freeze)
		submenu.add_command(label="UnFreeze",   underline=0,
				compound=LEFT, image=tkFlair.icons["heat"],
				command=self.editor.unfreeze)
		submenu.add_command(label="Toggle",  underline=0,
				compound=LEFT, image=tkFlair.icons["toggle"],
				accelerator="F",
				command=self.editor.toggleFreeze)

	# ----------------------------------------------------------------------
	def _addMenus(self, menu):
		bodymenu = Menu(menu, tearoff=1)
		menu.add_cascade(label="Bodies", menu=bodymenu, underline=0,
				compound=LEFT, image=tkFlair.icons["TRC"],
				accelerator="b")
		tkFlair.bodiesMenu(bodymenu, self.newBody)

		menu.add_command(label="Zone", underline=0,
				compound=LEFT, image=tkFlair.icons["zone"],
				command=self.setActionZone,
				accelerator="d")

		menu.add_command(label="REGION", underline=0,
				compound=LEFT, image=tkFlair.icons["REGION"],
				command=self.editor.addRegion,
				accelerator="R")

		menu.add_command(label="ROT-DEFI", underline=2,
				compound=LEFT, image=tkFlair.icons["ROT-DEFI"],
				command=self.editor.addRotdefi,
				accelerator="T")

		objmenu = Menu(menu)
		menu.add_cascade(label="Object", underline=0, menu=objmenu,
				compound=LEFT, image=tkFlair.icons["pin"])
		self._addObjectMenu(objmenu)

	# ----------------------------------------------------------------------
	def insertMenu(self, event=None):
		menu = Menu(self, tearoff=0)
		self._addMenus(menu)
		menu.tk_popup(self.winfo_pointerx(), self.winfo_pointery())

	# ----------------------------------------------------------------------
	def insertBodyMenu(self):
		menu = Menu(self, tearoff=1)
		tkFlair.bodiesMenu(menu, self.newBody)
		menu.tk_popup(self.winfo_pointerx(), self.winfo_pointery())

	# ----------------------------------------------------------------------
	def insertObjectMenu(self):
		menu = Menu(self, tearoff=0)
		self._addObjectMenu(menu)
		menu.tk_popup(self.winfo_pointerx(), self.winfo_pointery())

	# ----------------------------------------------------------------------
	# Add objects
	# ----------------------------------------------------------------------
	def newPoint(self, event=None):
		self.setStatusMsg("Adding Point object")
		self.actionOption = "!point"
		self.editor.setMouseAction(ACTION_ADD)

	# ----------------------------------------------------------------------
	def newArrow(self, event=None):
		self.setStatusMsg("Adding Arrow object")
		self.actionOption = "!arrow"
		self.editor.setMouseAction(ACTION_ADD)

	# ----------------------------------------------------------------------
	def newRuler(self, event=None):
		self.setStatusMsg("Adding Ruler object")
		self.actionOption = "!ruler"
		self.editor.setMouseAction(ACTION_ADD)

	# ----------------------------------------------------------------------
	def newMesh(self, event=None):
		self.setStatusMsg("Adding Mesh object")
		self.actionOption = "!mesh"
		self.editor.setMouseAction(ACTION_ADD)

	# ----------------------------------------------------------------------
	def newRotdefi(self, event=None):
		self.setStatusMsg("Adding Rotdefi object")
		self.actionOption = "!rotdefi"
		self.editor.setMouseAction(ACTION_ADD)

	# ----------------------------------------------------------------------
	def newLight(self, event=None):
		self.setStatusMsg("Adding Light object")
		self.actionOption = "!light"
		self.editor.setMouseAction(ACTION_ADD)

	# ----------------------------------------------------------------------
	def newBody(self, card=None):
		if card is not None: self.actionOption = card
		self.setStatusMsg("Adding body \'%s\'"%(self.actionOption))
		self.editor.setMouseAction(ACTION_ADD)

	# ----------------------------------------------------------------------
	# Action zone called from popup menu
	# ----------------------------------------------------------------------
#	def setActionNewZone(self):
#		self.editor.propList.selectClear()
#		self.editor.setMouseAction(ACTION_ZONE)

	# ----------------------------------------------------------------------
	def setActionZone(self):
		self.editor.setMouseAction(ACTION_ZONE)
