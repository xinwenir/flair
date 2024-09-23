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

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkinter
from Constants import *
import Input
import tkFlair
import tkExtra
import Unicode

#===============================================================================
# ToolDialog base class for a floating dialog inside the geometry editor
#===============================================================================
class ToolDialog(Frame):
	TITLE_BG = "DarkBlue"
	TITLE_FG = "White"

	TEXT_NO  = "Black"
	TEXT_OK  = "DarkGreen"
	TEXT_ERR = "Red"

	EXPAND   = Unicode.BLACK_RIGHT_POINTING_TRIANGLE
	COLLAPSE = Unicode.BLACK_DOWN_POINTING_TRIANGLE

	INITIAL  = -1000

	def __init__(self, master, title="", ok=False, *kw, **kwargs):
		Frame.__init__(self, master, *kw, **kwargs)

		self.bind("<Configure>", self._configure)
		self.master.bind("<Configure>", self._configure)

		self.action = IntVar()
		self.action.set(ACTION_SELECT)
		self.action.trace("w", self.actionChange)

		frame = Frame(self, bg=ToolDialog.TITLE_BG)
		frame.pack(side=TOP, fill=X)
		frame.bind("<B1-Motion>", self._motion)
		frame.bind("<ButtonRelease-1>", self._motionEnd)

		self._minimize = Label(frame, text=ToolDialog.COLLAPSE, bg=ToolDialog.TITLE_BG, fg=ToolDialog.TITLE_FG)
		self._minimize.pack(side=LEFT)
		self._minimize.bind("<Button-1>", self._toggle)

		self._title = Label(frame, text=title, bg=ToolDialog.TITLE_BG, fg=ToolDialog.TITLE_FG)
		self._title.pack(side=LEFT, fill=X) #, expand=YES)
		self._title.bind("<B1-Motion>", self._motion)
		self._title.bind("<ButtonRelease-1>", self._motionEnd)

		self._frame = Frame(self, relief=GROOVE, padx=2, pady=2)
		self._frame.pack(side=TOP, fill=BOTH, expand=YES, padx=2, pady=2)

		self._buttonFrame = Frame(self)
		self._buttonFrame.pack(side=TOP, fill=X)

		self._cancelButton = Button(self._buttonFrame, text="Cancel", command=self.cancel, padx=1, pady=1)
		self._cancelButton.pack(side=RIGHT)

		if ok:
			self._okButton = Button(self._buttonFrame, text="Ok", command=self.ok, padx=1, pady=1)
			self._okButton.pack(side=RIGHT)

		b = Checkbutton(self._buttonFrame, image=tkFlair.icons["info"],
						command=self.toggleInfo,
						indicatoron=0,
						padx=1, pady=1)
		b.pack(side=LEFT)

		self._info = Message(self, relief=SUNKEN)
		self.fillInfo()

		self._oldx = ToolDialog.INITIAL
		self._oldy = 50
		self._x0 = None
		self._y0 = None

	#-----------------------------------------------------------------------
	def geometry(self):	return self.master.geometry()
	def editor(self):	return self.master.editor

	#-----------------------------------------------------------------------
	def show(self, side=None):
		self.place(relx=0.0, rely=0.0, x=self._oldx, y=self._oldy, anchor=NW, bordermode="outside")
		self.wait_visibility()
		self.placeAt(side)

	#-----------------------------------------------------------------------
	def hide(self):
		if not self.winfo_ismapped(): return
		self._oldx = self.winfo_x()
		self._oldy = self.winfo_y()
		self.place_forget()

	#-----------------------------------------------------------------------
	def _configure(self, event=None):
		if not self.winfo_ismapped(): return
		self.placeAt()

	#-----------------------------------------------------------------------
	def placeAt(self, side=None):
		W = self.master.winfo_width()
		H = self.master.winfo_height()
		w = self.winfo_width()
		h = self.winfo_height()

		if side is None:
			if self._oldx == ToolDialog.INITIAL:
				x = W-w
				y = self._oldy
			else:
				x = self.winfo_x()
				y = self.winfo_y()
			if x<0:     x = 0
			elif x>W-w: x = W-w
			if y<0:	    y = 0
			elif y>H-h: y = H-h
		else:
			if side in (N, S, CENTER):
				x = (W-w)/2
			elif W in side:
				x = 0
			else:
				x = W-w

			if side in (E, W, CENTER):
				y = (H-h)/2
			elif N in side:
				y = 0
			else:
				y = H-h

		self.place(x=x, y=y)

	#-----------------------------------------------------------------------
	def _toggle(self, event=None):
		if self._minimize["text"] == ToolDialog.COLLAPSE:
			self._frame.pack_forget()
			self._buttonFrame.pack_forget()
			self._minimize["text"] = ToolDialog.EXPAND
		else:
			self._frame.pack(side=TOP, fill=BOTH, expand=YES)
			self._buttonFrame.pack(side=BOTTOM, fill=X)
			self._minimize["text"] = ToolDialog.COLLAPSE

	#-----------------------------------------------------------------------
	def _motion(self, event):
		if self._x0 is None:
			self._x0 = event.x
			self._y0 = event.y
		x = event.x_root - self.master.winfo_rootx() - self._x0
		y = event.y_root - self.master.winfo_rooty() - self._y0

		W = self.master.winfo_width()
		H = self.master.winfo_height()
		w = self.winfo_width()
		h = self.winfo_height()
		if x<0:     x = 0
		elif x>W-w: x = W-w
		if y<0:	    y = 0
		elif y>H-h: y = H-h

		self.place(x=x, y=y)

	#-----------------------------------------------------------------------
	def _motionEnd(self, event):
		self._x0 = None
		self._y0 = None

	#-----------------------------------------------------------------------
	def cancel(self):
		self.action.set(ACTION_SELECT)
		self.hide()

	#-----------------------------------------------------------------------
	def ok(self):
		self.hide()

	#-----------------------------------------------------------------------
	def toggleInfo(self):
		if self._info.winfo_ismapped():
			self._info.pack_forget()
		else:
			self._info.pack(side=BOTTOM, fill=BOTH, expand=YES)

	#-----------------------------------------------------------------------
	def fillInfo(self):
		#self._info.config(state=NORMAL)
		#self._info.delete("0.0",END)
		#self._info.insert("0.0",self.__doc__)
		#self._info.config(state=DISABLED)
		self._info["text"] = self.__doc__

	#-----------------------------------------------------------------------
	def nextAction(self):
		pass

	#-----------------------------------------------------------------------
	def setAction(self, act):
		self.action.set(act)

	#-----------------------------------------------------------------------
	def actionChange(self, a=None, b=None, c=None):
		self.master._setMouseAction(self.action.get())

	#-----------------------------------------------------------------------
	def addRadiobuttonRow(self, row, label, hlp, image, value, clear=None):
		l = Label(self._frame, text=label)
		l.grid(row=row, column=0, sticky=W)
		btn = Radiobutton(self._frame,
					image=image, compound=LEFT,
					variable=self.action,
					value=value,
					selectcolor="LightGreen",
					indicatoron=0)
		btn.grid(row=row, column=1, sticky=EW)
		tkExtra.Balloon.set(btn, hlp)
		if clear:
			b = Button(self._frame, text=Unicode.MATH_HEAVY_MULTIPLICATION,
					padx=1, pady=1,
					command=clear,
					relief=FLAT)
			b.grid(row=row, column=2, sticky=EW)
			tkExtra.Balloon.set(b, "Clear %s"%(label))
		return btn

#===============================================================================
# ZoneTool
#===============================================================================
class ZoneTool(ToolDialog):
	"""
To enable this action you need to SELECT:
- None or ONE Region
- ONE or MORE Bodies.
And to:
- Add a new zone: NO Zone selected
- Edit a zone: only ONE Zone selected
	"""
	def __init__(self, master,*kw, **kwargs):
		ToolDialog.__init__(self, master, "Zone Edit", *kw, **kwargs)

		# ===
		row = 0
		self._regionBtn = self.addRadiobuttonRow(row,
					"Region",
					"Region selection:\nNone - to add a new one\nOne - To append a new zone",
					tkFlair.icons["select"],
					ACTION_SELECT_REGION,
					self.clearRegion)
		row += 1
		self._zoneBtn   = self.addRadiobuttonRow(row,
					"Zone",
					"Zone selection:\nNone - to add a new zone\nOne - to modify the selected one",
					tkFlair.icons["select"],
					ACTION_SELECT_ZONE,
					self.clearZone)
		row += 1
		self._bodiesBtn = self.addRadiobuttonRow(row,
					"Bodies",
					"Select the bodies enclosing the zone you want to describe",
					tkFlair.icons["select"],
					ACTION_SELECT_BODIES,
					self.clearBodies)
		row += 1
		self._opBtn     = self.addRadiobuttonRow(row,
					"Operation",
					"Create a new region, new zone or modify an existing one [d]",
					tkFlair.icons["add"],
					ACTION_ADD_ZONE)

	#-----------------------------------------------------------------------
	def show(self, side=None):
		ToolDialog.show(self, side)
		self.updateSelection()

	#-----------------------------------------------------------------------
	def clearRegion(self):
		self.editor().selectionClear(True, Input.Card.REGION)
		self.action.set(ACTION_SELECT_REGION)

	#-----------------------------------------------------------------------
	def clearZone(self):
		self.editor().propListClear()
		if self.action.get() != ACTION_SELECT_BODIES:
			self.action.set(ACTION_SELECT_ZONE)

	#-----------------------------------------------------------------------
	def clearBodies(self):
		self.editor().selectionClear(True, Input.Card.BODY)
		self.action.set(ACTION_SELECT_BODIES)

	# ----------------------------------------------------------------------
	def updateSelection(self):
		nb = self.geometry().selection("b")
		nr = self.geometry().selection("r")
		nz = self.editor().propListZoneSelected()

		error = False
		if nb==0:
			self._bodiesBtn.config(text="Select", fg=ToolDialog.TEXT_NO)
			error = True
		else:
			self._bodiesBtn.config(text="%d selected"%(nb), fg=ToolDialog.TEXT_OK)
			self._opBtn.config(text="Add Zone", fg=ToolDialog.TEXT_NO, state=NORMAL)

		if nz==0:
			self._zoneBtn.config(text="Select", fg=ToolDialog.TEXT_NO)
			self._opBtn.config(text="Add Zone", fg=ToolDialog.TEXT_NO, state=NORMAL)
		elif nz==1:
			self._zoneBtn.config(text="1 selected", fg=ToolDialog.TEXT_OK)
			self._opBtn.config(text="Overwrite Zone", fg=ToolDialog.TEXT_NO, state=NORMAL)
			if self.action.get() == ACTION_SELECT_ZONE:
				self.action.set(ACTION_SELECT_BODIES)
		else:
			self._zoneBtn.config(text="%d selected"%(nz), fg=ToolDialog.TEXT_ERR)
			error = True

		if nr==0:
			self._regionBtn.config(text="Select", fg=ToolDialog.TEXT_NO)
			self._opBtn.config(text="Add Region", fg=ToolDialog.TEXT_NO, state=NORMAL)
		elif nr==1:
			self._regionBtn.config(text="1 selected", fg=ToolDialog.TEXT_OK)
			if self.action.get() == ACTION_SELECT_REGION:
				self.action.set(ACTION_SELECT_ZONE)
		else:
			self._regionBtn.config(text="%d selected"%(nr), fg=ToolDialog.TEXT_ERR)
			error = True

		if error:
			self._opBtn.config(text="Error", fg=ToolDialog.TEXT_ERR, state=DISABLED)

	# ----------------------------------------------------------------------
	def nextAction(self):
		# Bodies should be selected
		nb = self.geometry().selection("b")
		nr = self.geometry().selection("r")
		nz = self.editor().propListZoneSelected()

#		print
#		import traceback; traceback.print_stack()
#		print "Next Action", "B:",nb, "Z:", nz, "R:", nr

		if nr<=1 and nz<=1 and nb>0:
			self.action.set(ACTION_ADD_ZONE)
		elif nr>1:
			self.action.set(ACTION_SELECT_REGION)
		elif nr==1:
			self.action.set(ACTION_SELECT_BODIES)
		elif nz>1:
			self.action.set(ACTION_SELECT_ZONE)
		elif nb==0:
			self.action.set(ACTION_SELECT_BODIES)
		else:
			self.action.set(ACTION_ADD_ZONE)
#
