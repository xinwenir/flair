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
# Date:	03-Apr-2014

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import re
import sys
from stat import *

try:
	from tkinter import *
	from tkinter.colorchooser import askcolor
	import tkinter.font
except ImportError:
	from tkinter import *
	import tkinter.font as tkFont
	from tkinter.colorchooser import askcolor

import tkExtra

#import pdb

_FONT_SIZE = (   '6',  '7',  '8',  '9', '10', '11', '12', '13',
		'14', '16', '18', '20', '22', '24', '26', '28',
		'32', '36', '40', '48', '56', '64', '72' )

_FONT_JUSTIFY = (LEFT, CENTER, RIGHT)

try:
	import PIL.Image
	import PIL.ImageTk
	have_PIL = True
except:
	have_PIL = False

# User readable styles
ABBREV = {
	"Heading 1"  : "h1",
	"Heading 2"  : "h2",
	"Heading 3"  : "h3",
	"Bold"       : "b",
	"Italic"     : "i",
	"BoldItalic" : "bi",
	"Underline"  : "u",
	"Overstrike" : "o",
	"Normal"     : "n",
	"Link"       : "link",
	"Caption"    : "caption",

	"Left"       : "l",
	"Right"      : "r",
	"Center"     : "c",
	}
# and inverted ones
INV_ABBREV = dict([(y,x) for (x,y) in list(ABBREV.items())])

#===============================================================================
# RichText widget
#===============================================================================
class RichText(Text):
	TAGPAT = re.compile(r"{(\w+):\s*(.*?)}", re.MULTILINE|re.DOTALL)
	TAGS = {
		"n"   : {"font":("Sans",10,""), "justify":LEFT, "foreground":"#000000", "background":"#FFFFFF"},
		"h1"  : {"font":(None,16,"bold"), "foreground":"DarkBlue", "justify":CENTER},
		"h2"  : {"font":(None,14,"bold"), "foreground":"DarkBlue", "justify":CENTER},
		"h3"  : {"font":(None,14,"italic"), "foreground":"DarkBlue", "justify":CENTER},
		"b"   : {"font":(None,None,"bold")},
		"i"   : {"font":(None,None,"italic")},
		"bi"  : {"font":(None,None,"bold italic")},
		"u"   : {"font":(None,None,"underline")},
		"o"   : {"font":(None,None,"overstrike")},
		"l"   : {"justify":LEFT},
		"r"   : {"justify":RIGHT},
		"c"   : {"justify":CENTER},
		"img" : {},
		"link": {"font":(None,None,"underline"), "foreground":"Blue"},
		"caption" : {"foreground":"DarkBlue"}
		}
	LINETAGS = ["h1", "h2", "h3", "l", "c", "r"]

	# --------------------------------------------------------------------
	def __init__(self, master, *kw, **args):
		Text.__init__(self, master, *kw, **args)
		self.clear()

	# --------------------------------------------------------------------
	def clear(self):
		self.images = {}
		self._img   = {}

	# --------------------------------------------------------------------
	# Open file and highlight
	# @return list of errors
	# --------------------------------------------------------------------
	def open(self, filename):
		try:
			f = open(filename, "r")
		except IOError:
			# Error file not found!!!
			self.clear()
			return "File %s not found"%(filename)

		data = f.read()
		f.close()
		return self.set(data)

	# --------------------------------------------------------------------
	# Return style with inheritance
	# --------------------------------------------------------------------
	@staticmethod
	def getStyle(tag):
		style = RichText.TAGS.get(tag)
		if style is None: return None

		based = style.get("based")
		if based is None:
			based = RichText.TAGS.get("n")
		else:
			based = RichText.getStyle(based)
		if based is None: return None

		newstyle = based.copy()
		try: del newstyle["background"]
		except: pass
		for tag,value in list(style.items()):
			if tag == "based": continue
			if value is None: continue

			if isinstance(value,tuple):
				newvalue = []
				for i,v in enumerate(value):
					if v is None:
						v = based[tag][i]
					newvalue.append(v)
				newstyle[tag] = tuple(newvalue)
			else:
				# Replace with the new value
				newstyle[tag] = value

		return newstyle

	# --------------------------------------------------------------------
	def configureTags(self):
		# Configure tags
		for tag in list(RichText.TAGS.keys()):
			style = RichText.getStyle(tag)
			if style: self.tag_config(tag, style)

		# Basic colors of text window
		style = RichText.TAGS["n"]
		for n in ("font", "background", "foreground"):
			self[n] = style[n]

		self.tag_bind("image","<Double-1>", self.imageDouble)
		self.tag_bind("link", "<Double-1>", self.linkDouble)

	# --------------------------------------------------------------------
	def imageDouble(self, event):
		self.event_generate("<<Image>>")
		return "break"

	# --------------------------------------------------------------------
	def linkDouble(self, event):
		self.event_generate("<<Link>>")
		return "break"

	# --------------------------------------------------------------------
	def insertImage(self, filename):
		image = self.image(filename)
		if image:
			img = self.image_create(INSERT, image=image)
			self.tag_add("image",INSERT,INSERT)
			self._img[img] = filename

	# --------------------------------------------------------------------
	def imageName(self, img): return self._img.get(img)

	# --------------------------------------------------------------------
	# Set data to widget and highlight it
	# --------------------------------------------------------------------
	def setText(self, data):
		self.clear()
		last = 0
		self.delete("0.0",END)
		for pat in re.finditer(RichText.TAGPAT, data):
			# text outside tag
			self.insert(END, data[last:pat.start()])
			tag = pat.group(1)
			value = pat.group(2)
			if tag == "img":
				image = self.image(value)
				if image:
					img = self.image_create(END, image=image)
					self.tag_add("image",self.index(img))
					self._img[img] = value
				else:
					self.insert(END,pat.group(0))
			else:
				self.insert(END, value, tag)
			last = pat.end()
		self.insert(END, data[last:])

		self.configureTags()
		return None

	# --------------------------------------------------------------------
	# Get rich text from widget
	# --------------------------------------------------------------------
	def getText(self):
		data = ""
		tag = None
		for key,value,index in self.dump("0.0",END):
			if key=="tagon":
				if value not in ("sel","image"):
					if tag is not None: data += "}"
					tag = value
					data += "{%s:"%(tag)
			elif key=="text":
				data += value
			elif key=="tagoff":
				if tag is not None: data += "}"
				tag  = None
			elif key=="image":
				if tag is not None: data += "}"
				name = self._img.get(value)
				data += "{img:%s}"%(name)
				if tag is not None: data += "{%s: "%(tag)
			elif key=="mark":
				pass
			elif key=="window":
				pass
		return data

	# --------------------------------------------------------------------
	# Find image
	# --------------------------------------------------------------------
	def image(self, filename):
		if not have_PIL: return None
		try:
			s = os.stat(filename)
			mtime = s[ST_MTIME]
		except:
			return None

		image = self.images.get(filename)
		if image is not None:
			image, mt = image
			if mtime > mt: image = None

		if image is None:
			try:
				image = PIL.ImageTk.PhotoImage(PIL.Image.open(filename))
				self.images[filename] = image, mtime
			except:
				image = None
		return image

	# --------------------------------------------------------------------
	def setStyle(self, style):
		insert = self.index(INSERT)
		style = ABBREV.get(style,style)
		try:
			sel_first = self.index(SEL_FIRST)
			sel_last  = self.index(SEL_LAST)

			if style in RichText.LINETAGS:
				sel_first = "%s linestart"%(sel_first)
				sel_last  = "%s lineend"%(sel_last)
		except TclError:
			cur = self.index(INSERT)
			if style in RichText.LINETAGS:
				sel_first = "%s linestart"%(cur)
				sel_last  = "%s lineend"%(cur)
			else:
				sel_first = "%s - 1 char"%(cur)
				sel_last  = cur

		names = self.tag_names(sel_first)

		for n in names:
			if n in RichText.LINETAGS:
				first = "%s linestart"%(sel_first)
				last  = "%s lineend"%(sel_last)
				break
		else:
			first = sel_first
			last  = sel_last

		for tag in list(RichText.TAGS.keys()):
			self.tag_remove(tag, first, last)

		if style!="n" and style not in names:
			self.tag_add(style, sel_first, sel_last)

		# Why? + 1 chars (but it works)
		self.mark_set(INSERT,"%s + 1 chars"%(insert))

	# --------------------------------------------------------------------
	def styleList(self):
		styles = list(ABBREV.keys())
		# Add custom styles here
		styles.sort()
		return styles

	# --------------------------------------------------------------------
	def styleAtCursor(self):
		for n in self.tag_names(INSERT):
			if n in RichText.TAGS:
				return INV_ABBREV.get(n,n)
		return "Normal"

#===============================================================================
# Style configuration
#===============================================================================
class StyleConfigurationDialog(Toplevel):
	def __init__(self, master, text, active=None):
		Toplevel.__init__(self, master)
		self.title("Style configuration")
		self.transient(master)
		self.text = text

		# Variables
		self.fontChanged = False
		self.fontSize   = StringVar(self)
		self.fontAbs    = BooleanVar(self)
		self.fontBold   = BooleanVar(self)
		self.fontItalic = BooleanVar(self)
		self.fontUnderline = BooleanVar(self)
		self.fontStrike = BooleanVar(self)
		self.fontFamily = StringVar(self)
		self.font       = tkinter.font.Font(self, ('helvetica', 12, 'normal'))
		self.fontSize.set('12')
		self.fontFamily.set('helvetica')
		self.fontJustify = StringVar(self)
		self.fontJustify.set(LEFT)

		# Style selection frame
		frame = Frame(self, relief=RAISED)
		frame.pack(side=TOP, fill=X, pady=2)

		Label(frame, text="Style:").pack(side=LEFT)

		self.styleCombo = tkExtra.Combobox(frame, #False,
					width=12,
					background="White",
					command=self.setStyleFromCombo)
		self.styleCombo.pack(side=LEFT)
		self.styleCombo.bindWidgets("<Button-3>", self.toggleField)
		tkExtra.Balloon.set(self.styleCombo, "Select style to edit")

		# Style editing frame
		frame = Frame(self, relief=SUNKEN, pady=2, padx=1)
		frame.pack(expand=YES, fill=BOTH)

		f = Frame(frame)
		f.pack(side=TOP, padx=2, pady=2, fill=BOTH, expand=YES)

		#Label(f, justify=LEFT, text='Font:').pack(side=TOP, anchor=W)
		self.listFontName = tkExtra.ExListbox(f,
					height=5,
					takefocus=FALSE,
					exportselection=FALSE)
		self.listFontName.pack(side=LEFT, expand=TRUE, fill=BOTH)
		self.listFontName.bind('<ButtonRelease-1>',
					self.listFontButtonRelease)
		#self.listFontName.bind("<<ListboxSelect>>",
		#			self.listFontButtonRelease)
		self.listFontName.bind("<KeyRelease>",
					self.listFontKeyRelease)
		self.listFontName.bind("<Button-3>", self.toggleField)
		sb = Scrollbar(f)
		sb.pack(side=LEFT, fill=Y)
		sb.config(command=self.listFontName.yview)
		self.listFontName.config(yscrollcommand=sb.set)

		# --- Param ---
		f = Frame(frame)
		f.pack(padx=5, pady=2, fill=X)
		Label(f, text='Size:').pack(side=LEFT, anchor=W)
		self.menuFontSize = tkExtra.ExOptionMenu(f,
					self.fontSize,
					command=self.setFontSample,
					*_FONT_SIZE)
		self.menuFontSize.pack(side=LEFT, anchor=W)
		self.menuFontSize.bind("<Button-3>", self.toggleField)

		self.fontAbs.set(1)
		self.checkFontAbs = Checkbutton(f, text='DPI',
					variable=self.fontAbs,
					onvalue=1, offvalue=0,
					command=self.setFontSample)
		self.checkFontAbs.pack(side=LEFT, anchor=W)
		self.checkFontAbs.bind("<Button-3>", self.toggleField)

		self.checkFontBold = Checkbutton(f, text='Bold',
					variable=self.fontBold,
					onvalue=1, offvalue=0,
					command=self.setFontSample)
		self.checkFontBold.pack(side=LEFT, anchor=W)
		self.checkFontBold.bind("<Button-3>", self.toggleField)

		self.checkFontItalic = Checkbutton(f, text='Italic',
					variable=self.fontItalic,
					onvalue=1, offvalue=0,
					command=self.setFontSample)
		self.checkFontItalic.pack(side=LEFT, anchor=W)
		self.checkFontItalic.bind("<Button-3>", self.toggleField)

		self.checkFontUnderline = Checkbutton(f, text='Underline',
					variable=self.fontUnderline,
					onvalue=1, offvalue=0,
					command=self.setFontSample)
		self.checkFontUnderline.pack(side=LEFT, anchor=W)
		self.checkFontUnderline.bind("<Button-3>", self.toggleField)

		self.checkFontStrike = Checkbutton(f, text='Strike',
					variable=self.fontStrike,
					onvalue=1, offvalue=0,
					command=self.setFontSample)
		self.checkFontStrike.pack(side=LEFT, anchor=W)
		self.checkFontStrike.bind("<Button-3>", self.toggleField)

		# ---
		f = Frame(frame)
		f.pack(padx=5, pady=2, fill=X)
		Label(f, text="Color: Foreground").pack(side=LEFT)
		self.fontForeground = Button(f, text="         ", command=self.foregroundColor,
					background="Black", activebackground="Black")
		self.fontForeground.pack(side=LEFT)
		self.fontForeground.bind("<Button-3>", self.toggleField)

		Label(f, text="Background").pack(side=LEFT)
		self.fontBackground = Button(f, text="         ", command=self.backgroundColor,
					background="White", activebackground="White")
		self.fontBackground.pack(side=LEFT)
		self.fontBackground.bind("<Button-3>", self.toggleField)

		# ---
		f = Frame(frame)
		f.pack(padx=5, pady=2, fill=X)
		Label(f, text="Justify:").pack(side=LEFT)
		self.menuJustify = tkExtra.ExOptionMenu(f,
					self.fontJustify,
					command=self.setFontSample,
					*_FONT_JUSTIFY)
		self.menuJustify.pack(side=LEFT, anchor=W)

		# --- Sample ---
		f = Frame(frame, relief=SOLID, borderwidth=1)
		f.pack(side=TOP, padx=5, pady=5, expand=TRUE, fill=BOTH)

		self.labelFontSample = Label(f,
			text=   "AaBbCcDdEe\nFfGgHhIiJjK\n"
				"A\u03b1B\u03b2\u0393\u03b3\u0394\u03b4E\u03b5\n"
				"Z\u03b6H\u03b7\u0398\u03b8I\u03b9K\u03ba\n"
				"\u039b\u03bbM\u03bcN\u03bd\u039e\u03beOo\n"
				"\u03a0\u03c0P\u03c1\u03a3\u03c3T\u03c4Y\u03c5\u03a6\n"
				"\u03d5X\u03c7\u03a8\u03c8\u03a9\u03c9\n"
				"1234567890\n#:+=(){}[]",
			justify=LEFT, font=self.font)
		self.labelFontSample.pack(expand=TRUE, fill=BOTH)

		self.protocol("WM_DELETE_WINDOW", self.close)
		self.bind('<Escape>', self.close)

		# Button frame
		frame = Frame(self)
		frame.pack(side=BOTTOM, fill=X)
		Button(frame, text="Help", state=DISABLED, command=self.help).pack(side=RIGHT)
		Button(frame, text="Close", command=self.close).pack(side=RIGHT)
		Button(frame, text="Apply", command=self.apply).pack(side=RIGHT)

		# Load fonts
		self.loadFonts()
		self.styleCombo.fill(self.text.styleList())
		self._current = None
		if active:
			self.styleCombo.set(active)
		else:
			self.styleCombo.set("Normal")

		# Wait action
		self.wait_visibility()
		self.grab_set()
		self.focus_set()
		self.wait_window()

	# --------------------------------------------------------------------
	def apply(self, event=None):
		if self.fontChanged: self.saveStyle()
		if self.text is not None:
			self.text.configureTags()

	# --------------------------------------------------------------------
	def close(self, event=None):
		self.apply()
		self.destroy()

	# --------------------------------------------------------------------
	def saveStyle(self):
		name  = self.styleCombo.get()
		short = ABBREV.get(name, name)
		style = RichText.getStyle(short)

		if self.listFontName["state"] == NORMAL:
			family = self.fontFamily.get()
		else:
			family = None

		if self.menuFontSize["state"] == NORMAL:
			size = int(self.fontSize.get())
			if not self.fontAbs.get():
				size = -size
		else:
			size = None

		options = ""
		if self.checkFontBold["state"] == NORMAL or \
		   self.checkFontItalic["state"] == NORMAL or \
		   self.checkFontUnderline["state"] == NORMAL or \
		   self.checkFontStrike["state"] == NORMAL:
			if self.fontBold.get(): options += " bold"
			if self.fontItalic.get(): options += " italic"
			if self.fontUnderline.get(): options += " underline"
			if self.fontStrike.get(): options += " overstrike"

		style["font"] = (family, size, options.strip())

		if self.menuJustify["state"] == NORMAL:
			style["justify"] = self.fontJustify.get()

		if self.fontForeground["state"] == NORMAL:
			style["foreground"] = self.fontForeground["background"]

		if self.fontBackground["state"] == NORMAL:
			b = self.fontBackground["background"]
			if short != "n" and RichText.TAGS["n"]["background"] == b:
				try: del style["background"]
				except: pass
			else:
				style["background"] = b
		RichText.TAGS[short] = style
		self.fontChanged = False

	# --------------------------------------------------------------------
	def setStyleFromCombo(self):
		# Save current!
		if self._current and self.fontChanged:
			self.saveStyle()

		# New style
		name  = self.styleCombo.get()
		short = ABBREV.get(name, name)
		if short == self._current: return
		style = RichText.getStyle(short)

		self._current = short

		# Font
		family, size, options = style["font"]
		self.fontFamily.set(family)
		self.fontSize.set(size)
		self.fontBold.set("bold" in options)
		self.fontItalic.set("italic" in options)
		self.fontUnderline.set("underline" in options)
		self.fontStrike.set("overstrike" in options)

		self.fontJustify.set(style.get("justify",LEFT))

		c = style.get("foreground","#000000")
		self.fontForeground.config(background=c, activebackground=c)
		c = style.get("background","#FFFFFF")
		self.fontBackground.config(background=c, activebackground=c)

		self.listFontName.selection_clear(0, END)
		for i in range(self.listFontName.size()):
			if self.listFontName.get(i) == family:
				self.listFontName.selection_set(i)
				self.listFontName.selection_anchor(i)
				self.listFontName.activate(i)
				self.listFontName.see(i)
				break
		self.setFontSample()
		self.fontChanged = False

	# --------------------------------------------------------------------
	def loadFonts(self):
		fonts = list(tkinter.font.families(self))
		fonts.sort()
		# remove duplicates
		self.listFontName.insert(END,'')
		for font in fonts:
			self.listFontName.insert(END, font)
		self.listFontName.see(0)
		self.listFontName.activate(0)
		self.listFontName.selection_set(0)

	# --------------------------------------------------------------------
	def listFontButtonRelease(self, event=None):
		font = self.listFontName.get(ANCHOR)
		self.fontFamily.set(font)
		self.setFontSample()

	# --------------------------------------------------------------------
	def listFontKeyRelease(self, event=None):
		self.after_idle(self.listFontButtonRelease)

	# --------------------------------------------------------------------
	def setFontSample(self, event=None):
		self.fontChanged = True
		fontFamily = self.fontFamily.get()

		if self.fontAbs.get():
			fontSign = 1
		else:
			fontSign = -1

		if self.fontBold.get():
			fontWeight = tkinter.font.BOLD
		else:
			fontWeight = tkinter.font.NORMAL

		if self.fontItalic.get():
			fontSlant = tkinter.font.ITALIC
		else:
			fontSlant = tkinter.font.ROMAN

		fontUnder  = bool(self.fontUnderline.get())
		fontStrike = bool(self.fontStrike.get())

		try:
			self.font.config(
				size=fontSign*int(self.fontSize.get()),
				weight=fontWeight,
				slant=fontSlant,
				underline=fontUnder,
				overstrike=fontStrike,
				family=fontFamily)
		except UnicodeEncodeError:
			pass

		j = self.fontJustify.get()
		self.labelFontSample["justify"] = j
		if   j==LEFT:  self.labelFontSample["anchor"] = W
		elif j==RIGHT: self.labelFontSample["anchor"] = E
		else:          self.labelFontSample["anchor"] = CENTER

		self.labelFontSample["foreground"] = self.fontForeground["background"]
		self.labelFontSample["background"] = self.fontBackground["background"]

	# --------------------------------------------------------------------
	def chooseColorFromButton(self, button, title):
		try:
			rgb, colorStr = askcolor(title=title,
						initialcolor=button["background"],
						parent=self)
		except TclError:
			colorStr = None
		if colorStr is not None:
			colorStr = str(colorStr)
			button.config(background=colorStr, activebackground=colorStr)
			self.setFontSample()

	# --------------------------------------------------------------------
	def foregroundColor(self, event=None):
		self.chooseColorFromButton(self.fontForeground, "Foreground color")

	# --------------------------------------------------------------------
	def backgroundColor(self, event=None):
		self.chooseColorFromButton(self.fontBackground, "Background color")

	# ----------------------------------------------------------------------
	def toggleField(self, event):
		if event.widget["state"] == NORMAL:
			state = DISABLED
		else:
			state = NORMAL

		if isinstance(event.widget.master, tkExtra.Combobox):
			event.widget.master.config(state=state)
		else:
			event.widget.config(state=state)

	# ----------------------------------------------------------------------
	def help(self):
		pass

# -----------------------------------------------------------------------------
if __name__ == "__main__":
	def close(event=None):
		global root
		root.destroy()

	root = Tk()
	root.protocol("WM_DELETE_WINDOW", root.destroy)
	root.bind("<Escape>", close)
	root.bind("<Control-Key-q>", close)
	root.bind("<Control-Key-w>", close)

	toolbar = Frame(root)
	toolbar.pack(side=TOP, fill=X)

	text = RichText(root)
	text.pack(side=TOP, expand=YES, fill=BOTH)

	for name in list(RichText.TAGS.keys()):
		cmd = lambda w=text,n=name: w.setStyle(n)
		b = Button(toolbar, text=name, command=cmd)
		b.pack(side=LEFT)

	text.open(sys.argv[1])
	#print text.getText()

	root.mainloop()
