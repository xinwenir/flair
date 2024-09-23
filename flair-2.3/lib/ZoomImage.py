# 
# -*- coding: latin1 -*-
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

from tkinter import *
import tkinter as tk
import sys
from math import ceil
from PIL import ImageTk

#===============================================================================
# Zoom image canvas
#===============================================================================
class ZoomImage(Canvas):
	def __init__(self, master, *arg, **kw_args):
		Canvas.__init__(self, master, *arg, **kw_args)
		self.image = self.create_image(0,0, anchor='nw')

		self.bind("<Configure>", self.resize)

		self.tkimage  = None
		self.photo    = None
		self.zoom     = 1.0
		self.scandrag = False

	# ----------------------------------------------------------------------
	def bindZoom(self):
		self.bind('<Button-4>', self.zoomIn);
		self.bind('<Button-5>', self.zoomOut);

		self.bind('<Key-minus>', self.zoomOut)
		self.bind('<Key-equal>', self.zoomIn)
		self.bind('<Key-plus>',  self.zoomIn)
		self.bind('<Key-1>',     self.zoomOne)

	# ----------------------------------------------------------------------
	def bindMotion(self, wheel=True):
		if wheel:
			self.bind('<Shift-Button-4>',     lambda e,s=self: s.xview(SCROLL, -1, UNITS))
			self.bind('<Shift-Button-5>',     lambda e,s=self: s.xview(SCROLL,  1, UNITS))
			self.bind('<Control-Button-4>',   lambda e,s=self: s.yview(SCROLL, -1, UNITS))
			self.bind('<Control-Button-5>',   lambda e,s=self: s.yview(SCROLL,  1, UNITS))

		self.bind('<B2-Motion>',       self.scanDrag)
		self.bind('<ButtonRelease-2>', self.scanRelease)

	# ----------------------------------------------------------------------
	def bindArrows(self):
		self.bind('<Key-Up>',    lambda e,s=self: s.yview(SCROLL, -1, UNITS))
		self.bind('<Key-Down>',  lambda e,s=self: s.yview(SCROLL,  1, UNITS))
		self.bind('<Key-Left>',  lambda e,s=self: s.xview(SCROLL, -1, UNITS))
		self.bind('<Key-Right>', lambda e,s=self: s.xview(SCROLL,  1, UNITS))

	# ----------------------------------------------------------------------
	def setImage(self, img=None):
		self.photo = img
		if img is None:
			self.itemconfig(self.image, image=None)
		else:
			self.zoomImage(self.zoom)
		self.event_generate("<<Zoom>>")

	# ----------------------------------------------------------------------
	def zoomIn(self, event=None):
		if self.zoom > 20.0: return
		if event is not None:
			x0,y0 = event.x, event.y
		else:
			x0 = y0 = None

		if self.zoom < 1.0:
			self.zoomImage(self.zoom+0.125,x0,y0)
		else:
			self.zoomImage(self.zoom+1.0,x0,y0)
		self.event_generate("<<Zoom>>")

	# ----------------------------------------------------------------------
	def zoomOut(self, event=None):
		if self.zoom <= 0.25: return
		if event is not None:
			x0,y0 = event.x, event.y
		else:
			x0 = y0 = None
		if self.zoom <= 1.0:
			self.zoomImage(self.zoom-0.125,x0,y0)
		else:
			self.zoomImage(self.zoom-1.0,x0,y0)
		self.event_generate("<<Zoom>>")

	# ----------------------------------------------------------------------
	def zoomOne(self, event=None):
		if event is not None:
			x0,y0 = event.x, event.y
		else:
			x0 = y0 = None
		self.zoomImage(1,x0,y0)
		self.event_generate("<<Zoom>>")

	# ----------------------------------------------------------------------
	def zoomImage(self, zoom, x0=None, y0=None):
		if self.photo is None: return
		if abs(zoom-1.0)<1e-10: zoom = 1.0
		oldZoom = self.zoom
		self.zoom = zoom

		if abs(zoom - 1.0) < 1E-9:
			self.tkimage = ImageTk.PhotoImage(self.photo)
			self.itemconfig(self.image, image=self.tkimage)
			self.coords(self.image, 0, 0)
			self.configure(scrollregion=(
				(0,0,	int(self.photo.size[0]*zoom),
					int(self.photo.size[1]*zoom))))

		elif zoom < 1.0:
			w, h = self.photo.size
			img = self.photo.copy()
			self.tkimage = ImageTk.PhotoImage(
				img.resize((int(w*zoom), int(h*zoom))))
			self.coords(self.image, 0, 0)
			cx, cy = 0.0, 0.0
			self.configure(scrollregion=(
				(0,0,	int(self.photo.size[0]*zoom),
					int(self.photo.size[1]*zoom))))
			self.itemconfig(self.image, image=self.tkimage)

		else:
			self.configure(scrollregion=(
				(0,0,	int(self.photo.size[0]*zoom),
					int(self.photo.size[1]*zoom))))
			if x0 is not None:
				cx = (zoom*(self.canvasx(0) + x0))/oldZoom - x0
				cy = (zoom*(self.canvasy(0) + y0))/oldZoom - y0
				cx = max(0.0, cx)
				cy = max(0.0, cy)
				self.xview_moveto(float(cx)/(self.photo.size[0]*zoom))
				self.yview_moveto(float(cy)/(self.photo.size[1]*zoom))
			self.updateImage()

	# ----------------------------------------------------------------------
	def updateImage(self, nx=None, ny=None):
		if self.photo is None: return
		x = int(self.canvasx(0) / self.zoom)
		y = int(self.canvasy(0) / self.zoom)
		w = int(ceil(self.winfo_width() / self.zoom))+1
		h = int(ceil(self.winfo_height() / self.zoom))+1
		img = self.photo.crop((x,y,x+w,y+h))

		self.tkimage = ImageTk.PhotoImage(
			img.resize((int(w*self.zoom), int(h*self.zoom))))

		self.coords(self.image, x*self.zoom, y*self.zoom)
		self.itemconfig(self.image, image=self.tkimage)

	# ----------------------------------------------------------------------
	def xview(self, *args):
		Canvas.xview(self, *args)
		if self.zoom > 1.0:
			self.updateImage()
			if args: self.event_generate("<<View>>")

	# ----------------------------------------------------------------------
	def yview(self, *args):
		Canvas.yview(self, *args)
		if self.zoom > 1.0:
			self.updateImage()
			if args: self.event_generate("<<View>>")

	# ----------------------------------------------------------------------
	def scan_dragto(self, x, y):
		Canvas.scan_dragto(self, x, y, int(max(1,self.zoom/2.0)))
		if self.zoom > 1.0:
			self.updateImage()

	# ----------------------------------------------------------------------
	def resize(self, event):
		if self.zoom > 1.0:
			self.updateImage()

	# ----------------------------------------------------------------------
	def scanDrag(self, event):
		if self.scandrag:
			self.scan_dragto(event.x, event.y)
			self.event_generate("<<View>>")
		else:
			self.config(cursor="hand2")
			self.scan_mark(event.x, event.y)
			self.scandrag = True

	# ----------------------------------------------------------------------
	def scanRelease(self, event):
		self.scandrag = False
		self.config(cursor="")

#-------------------------------------------------------------------------------
if __name__ == "__main__":
	import PIL.Image

	root = Tk()
	root.geometry("800x600")
	root.title("ZoomImage")

	zi = ZoomImage(root)
	zi.pack(expand=YES, fill=BOTH)
	zi.config(cursor="watch")
	try:
		image = PIL.Image.open(sys.argv[1])
	except:
		#print >>sys.stderr, "Error opening image file %s"%(sys.argv[1])
		sys.exit(1)
	zi.setImage(image)
	zi.config(cursor="")
	root.mainloop()
