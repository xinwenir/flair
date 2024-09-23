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

# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	21-Nov-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import re
import sys
import glob
import string
from stat import *
from log import say

try:
	from tkinter import *
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	import tkinter.messagebox as messagebox
import tkinter as tk
import tkFlair
import tkExtra
import Unicode
import tkDialogs

_directories = None
current     = None
images      = {}

#-------------------------------- Names ----------------------------------------
_change = [	("\\`e",  "\u00e8"),
		("\\'e",  "\u00e9"),
		("\\`o",  "\u00f2"),
		("\\\"u", "\u00fc"),
		("\\o",   "\u00f6"),
		("\\ie\\", "i.e.") ]
_patterns = [	re.compile(r"^(\d+)} +[A-Z].*$"),
		re.compile(r"^([A-Za-z\-]\d+)} +[A-Za-z].*$"),
		re.compile(r"^{([A-Za-z0-9\-:]+)}[a-z]*$"),
		re.compile(r"^{([A-Za-z0-9\-:]+)} +[A-Za-z].*$") ]
_subpatterns = [ re.compile(r"^ *(\d+\.\d+)} +[A-Za-z0-9].*$"),
		 re.compile(r"^ *(\d+\.\d+.\d+)} +[A-Za-z].*$"),
		 re.compile(r"^ *(\d+\.\d+.\d+\.\d+)} +[A-Za-z].*$"),
		 re.compile(r"^ *([A-Za-z]\d+\.\d+)} +[A-Za-z].*$"),
		 re.compile(r"^ *([A-Za-z]\d+\.\d+\.\d+)} +[A-Za-z].*$"),
		 re.compile(r"^ *([A-Za-z]\d+\.\d+\.\d+\.\d+)} +[A-Za-z].*$") ]
_patbody1   = re.compile(r"^ +\d+\. [A-Z].*\. Code: ([A-Z].*)$")
_patbody2   = re.compile(r"^ +\d+\.\d+\.\d+\.\d+} [A-Z].*\. Code: ([A-Z].*)$")
_patwhat    = re.compile(r"\bWHAT\([1-9.]+\)")
_patsdum    = re.compile(r"\bSDUM\b")
_patcomment = re.compile(r"^ \*.*$", re.MULTILINE)
_pathtml    = re.compile(r"<([a-z]+)[ >].*?</\1>")

_MAX_VIEWS  = 20

#===============================================================================
class Manual:
	# Manual is stored as a dictionary of sections
	manual     = {}
	section    = []
	title      = []
	pattern    = []
	topSection = None

	# ----------------------------------------------------------------------
	def __init__(self, path):
		if len(Manual.manual) == 0:
			for dir in path:
				filename = Manual.findManual(dir)
				if filename:
					Manual.load(filename)
		self.views   = []
		self.curview = 0
		self._move2End = False

	# ----------------------------------------------------------------------
	# show manual window
	# ----------------------------------------------------------------------
	def show(self, master=None):
		global current

		if master is not None:
			self.toplevel = master
		else:
			# Start a new interpreter (not Toplevel),
			# to have the manual page responsive at all times
			self.toplevel = Tk()

		self.toplevel.title("fluka/flair manual")

		self.defaultFont = tkFlair.getFont("manual.default",("Sans",-12), self.toplevel)
		self.fixedFont   = tkFlair.getFont("manual.fixed",  ("Mono",-12), self.toplevel)
		self.titleFont   = tkFlair.getFont("manual.title",  ("Sans",-18), self.toplevel)

		# Paned Window
		self.splitter = tkExtra.HSplitter(self.toplevel, 250, True)
		self.splitter.pack(side=TOP, fill=BOTH, expand=YES)
		tocFrame  = self.splitter.leftFrame()
		textFrame = self.splitter.rightFrame()

		# Search
		self.backwardBtn = Button(tocFrame, text=Unicode.LEFTWARDS_ARROW,
					font=self.defaultFont,
					padx=2, pady=0,
					command=self.backward)
		self.backwardBtn.grid(row=0, column=0, sticky=NSEW)
		self.forwardBtn = Button(tocFrame, text=Unicode.RIGHTWARDS_ARROW,
					font=self.defaultFont,
					padx=2, pady=0,
					command=self.forward)
		self.forwardBtn.grid(row=0, column=1, sticky=NSEW)
		self.entry = Entry(tocFrame,
					font=self.defaultFont,
					background='White')
		self.entry.grid(row=0, column=2, sticky=NSEW)
		b = Button(tocFrame, text=Unicode.MATH_HEAVY_MULTIPLICATION,
				font=self.defaultFont,
				padx=0, pady=0,
				command=self.clear)
		b.grid(row=0, column=3, sticky=NSEW)

		self.tocList = tkExtra.ExListbox(tocFrame, borderwidth=0,
				selectmode=BROWSE,
				selectborderwidth=0,
				background="White",
				font=self.defaultFont,
				takefocus=True,
				exportselection=FALSE)
		self.tocList.grid(row=1, column=0, columnspan=3, sticky=NSEW)

		sb = Scrollbar(tocFrame, orient=VERTICAL, takefocus=0,
				command=self.tocList.yview)
		sb.grid(row=1, column=3, sticky=NSEW)
		self.tocList.config(yscrollcommand=sb.set)
		tocFrame.grid_rowconfigure(1, weight=1)
		tocFrame.grid_columnconfigure(2, weight=1)

		# Text Viewer
		self.txt = Text(textFrame, background="White", font=self.fixedFont, takefocus=1)
		self.txt.config(wrap=NONE, cursor="top_left_arrow")
		self.txt.grid(row=0, column=0, sticky='nsew')
		sbv = tkExtra.AutoScrollbar(textFrame, takefocus=0)
		sbv.grid(row=0, column=1, sticky='ns')
		sbh = tkExtra.AutoScrollbar(textFrame, takefocus=0,
				orient=HORIZONTAL)
		sbh.grid(row=1, column=0, sticky='ew')

		self.txt.config(yscrollcommand=sbv.set, xscrollcommand=sbh.set)
		sbv.config(command=self.txt.yview)
		sbh.config(command=self.txt.xview)

		textFrame.grid_rowconfigure(0, weight=1)
		textFrame.grid_columnconfigure(0, weight=1)

		# Bindings
		self.txt.bind('<Button-1>',		self.textFocus)
		self.tocList.bind('<ButtonRelease-1>',	self.searchClick)
		self.tocList.bind('<Return>',		self.searchReturn)
		self.tocList.bind('<KP_Enter>',		self.searchReturn)
		self.toplevel.bind('<Control-Key-f>',	self.findDialog)
		self.toplevel.bind('<Control-Key-g>',	self.findNext)
		self.entry.bind('<Return>',		self.entryReturn)
		self.entry.bind('<KP_Enter>',		self.entryReturn)
		self.txt.bind('<Up>',	lambda e, s=self.txt:
						s.yview(SCROLL, -1, UNITS))
		self.txt.bind('<Down>',	lambda e, s=self.txt:
						s.yview(SCROLL,  1, UNITS))
		self.toplevel.bind("<F11>",		self.toggleTree)
		self.toplevel.bind("<F12>",		self.toggleHeight)
		self.toplevel.bind("<Escape>",		self.close)
		self.toplevel.bind('<Control-Key-w>',	self.close)
		self.toplevel.bind('<Control-Key-q>',	self.close)
		self.toplevel.protocol("WM_DELETE_WINDOW", self.close)

		self.tocList.bind('<Button-3>',		self.popupMenu)
		self.tocList.bind('<Control-Key-space>',self.popupMenu)
		self.tocList.bind('<<ListboxSelect>>',	self.searchAfter)
		self.txt.bind('<Key-slash>',		self.findDialog)
		self.txt.bind('<Button-3>',		self.popupMenu)
		self.txt.bind('<Control-Key-space>',	self.popupMenu)
		self.txt.bind('<Prior>',		self.prior)
		self.txt.bind('<Next>',			self.next)
		self.txt.bind('<Control-Prior>',	self.pageUp)
		self.txt.bind('<Control-Next>',		self.pageDown)

		# Fill TOC
		first = True
		for title in Manual.title:
			if not first and title[0:3] == "---":
				self.tocList.insert(END,"")
			self.tocList.insert(END, title)
			if title[0]!=" ":
				if title[0]=="-":
					self.tocList.itemconfig(END,foreground="DarkRed")
				else:
					self.tocList.itemconfig(END,foreground="DarkBlue")
			first = False

		# Initial values
		self.current_section = None
		self.searchStr       = None
		current              = self
		self.showSection(Manual.section[0])

		# Show window
		self.toplevel.update_idletasks()
		self.toplevel.deiconify()
		self.toplevel.focus_set()
		self.tocList.focus_set()
		self._oldHeight = 700
		self.toplevel.wait_visibility()
		self.toplevel.geometry("1024x800")
		self.splitter_pos = self.splitter.split
		self.inafter = False

	# ----------------------------------------------------------------------
	def close(self, event=None):
		global current
		if isinstance(self.toplevel, Tk):
			# Destroy local icons...
			done()
		self.toplevel.destroy()
		current = None

	# ----------------------------------------------------------------------
	def textFocus(self, event):
		self.txt.focus_set()

	# ----------------------------------------------------------------------
	def toggleTree(self, event):
		if self.splitter.split > self.splitter.min:
			self.splitter.doubleClick()
			self.splitter.rightFrame().focus_set()
			self.toplevel.focus_set()
		else:
			self.splitter.doubleClick()

	# ----------------------------------------------------------------------
	def toggleHeight(self, event):
		self._oldHeight = tkExtra.toggleHeight(self.toplevel,
					self._oldHeight)

	# ----------------------------------------------------------------------
	def clear(self, event=None):
		self.entry.delete(0,END)
		self.search("")

	# ----------------------------------------------------------------------
	def backward(self, event=None):
		if self.curview>=len(self.views):
			self.curview -= 1
		if self.curview>0:
			self.curview -= 1
			section, searchstr = self.views[self.curview]
			if searchstr != self.entry.get():
				self.entry.delete(0,END)
				self.entry.insert(0, searchstr)
				self.search(searchstr, False)
			self.showSection(section)

	# ----------------------------------------------------------------------
	def forward(self, event=None):
		if self.curview < len(self.views)-1:
			self.curview += 1
			section, searchstr = self.views[self.curview]
			if searchstr != self.entry.get():
				self.entry.delete(0,END)
				self.entry.insert(0, searchstr)
				self.search(searchstr, False)
			self.showSection(section)

	# ----------------------------------------------------------------------
	# Create the popup menu
	# ----------------------------------------------------------------------
	def popupMenu(self, event):
		menu=Menu(self.toplevel, tearoff=0)
		menu.add_command(label='Backward', underline=0,
				command=self.backward)
		menu.add_command(label='Forward', underline=0,
				command=self.forward)
		menu.add_command(label='Clear', underline=0,
				command=self.clear)
		menu.add_separator()
		menu.add_command(label='Find', underline=0,
				command=self.findDialog)
		menu.add_command(label='Find Next', underline=5,
				command=self.findNext)
		menu.add_separator()
		menu.add_command(label='Close',
				command=self.close)
		menu.tk_popup(event.x_root, event.y_root)
		return "break"

	# ----------------------------------------------------------------------
	# Return hit on entry box
	# ----------------------------------------------------------------------
	def entryReturn(self, event=None):
		self.search(self.entry.get())

	# ----------------------------------------------------------------------
	# Search for an expression in help pages
	# ----------------------------------------------------------------------
	def search(self, expr, select=True):
		self.searchStr = expr.upper().split()
		self.current_section = None
		self.tocList.delete(0,END)
		for i in range(len(Manual.section)):
			section = Manual.section[i]
			txt = Manual.manual[section].upper()
			found = True
			for w in self.searchStr:
				if txt.find(w) == -1:
					found = False
					break
			if found:
				title = Manual.title[i]
				self.tocList.insert(END,title)
				if title[0]!=" ":
					if title[0]=="-":
						self.tocList.itemconfig(END,foreground="DarkRed")
					else:
						self.tocList.itemconfig(END,foreground="DarkBlue")
		if select and self.tocList.size() > 0:
			self.tocList.selection_clear(0,END)
			self.tocList.selection_set(0)
			self.tocList.activate(0)
			self.tocList.focus_set()
			self.searchReturn()

	# ----------------------------------------------------------------------
	# button click
	# ----------------------------------------------------------------------
	def searchClick(self, event):
		self.tocList.focus_set()
		idx = self.tocList.nearest(event.y)
		if idx<0: return
		title = self.tocList.get(idx)
		if len(title)==0:
			idx += 1
			title   = self.tocList.get(idx)
		item    = Manual.title.index(title)
		section = Manual.section[item]
		self.showSection(section)

	# ----------------------------------------------------------------------
	def searchReturn(self, event=None):
		if self.tocList.size()==0: return
		sel = self.tocList.index(ACTIVE)
		if sel<0: return
		self.tocList.selection_clear(0,END)
		self.tocList.selection_set(sel)
		title = self.tocList.get(sel)
		if len(title) == 0:
			self.tocList.selection_clear(sel)
			sel = int(sel) + 1
			self.tocList.selection_set(sel)
			title = self.tocList.get(sel)
		item    = Manual.title.index(title)
		section = Manual.section[item]
		self.showSection(section)
		if self._move2End:
			self.txt.see(END)
			self._move2End = False

	# ----------------------------------------------------------------------
	def searchAfter(self, event=None):
		if self.inafter: return
		self.inafter = True
		self.toplevel.after(1000, self.searchReturn)

	# ----------------------------------------------------------------------
	def prior(self, event=None):
		top, bot = self.txt.yview()
		diff = bot-top
		if top < 0.001:
			self.pageUp()
			return "break"

	# ----------------------------------------------------------------------
	def next(self, event=None):
		top, bot = self.txt.yview()
		diff = bot-top
		if bot > 1.0 - diff/5.0:
			self.pageDown()
			return "break"

	# ----------------------------------------------------------------------
	def pageUp(self, event=None):
		active = self.tocList.index(ACTIVE)
		if active==0: return
		try: active = int(active)-1
		except: pass
		self.tocList.activate(active)
		self._move2End = True
		self.searchReturn()

	# ----------------------------------------------------------------------
	def pageDown(self, event=None):
		active = self.tocList.index(ACTIVE)
		try: active = int(active)+1
		except: pass
		self.tocList.activate(active)
		self.searchReturn()

	# ----------------------------------------------------------------------
	# show section
	# ----------------------------------------------------------------------
	def showSection(self, section=None, err=True):
		self.inafter = False
		if section is None:
			section = Manual.topSection

		try:
			idx = Manual.section.index(section)
		except:
			# Try to find the best match
			idx = 0
			su = section[:8].upper()
			for s in Manual.section:
				if su == s[:8].upper():
					section = s
					break
				idx += 1

		# Search item in list to highlight
		try: title=Manual.title[idx]
		except IndexError:
			if err:
				try:
					messagebox.showwarning("Section not found",
							"Section %s not found"%(section),
							parent=self.toplevel)
				except TclError:
					pass
			return False
		for i in range(0, self.tocList.size()):
			item = self.tocList.get(i)
			if item == title or item.strip() == section:
				self.tocList.selection_clear(0,END)
				self.tocList.selection_set(i)
				self.tocList.activate(i)
				self.tocList.see(i)
				break

		if self.current_section == section: return True
		self.current_section = section

		# Add to view history
		if self.curview == 0:
			self.backwardBtn.config(state=DISABLED)
		else:
			self.backwardBtn.config(state=NORMAL)
		if self.curview >= len(self.views):
			self.curview += 1
			self.views.append((section, self.entry.get()))
			if len(self.views) > _MAX_VIEWS:
				self.views.pop(0)
				self.curview -= 1
		else:
			if self.views[self.curview][0] != section:
				# delete the rest and append the new view
				self.curview += 1
				del self.views[self.curview:]
				self.views.append((section, self.entry.get()))
		if self.curview >= len(self.views)-1:
			self.forwardBtn.config(state=DISABLED)
		else:
			self.forwardBtn.config(state=NORMAL)

		txt = self.txt
		txt.config(state=NORMAL)
		txt.delete(0.0, END)
		text = Manual.manual.get(section,"")
		txt.insert(0.0, text)

		# Highlight items
		txt.tag_delete("title")
		txt.tag_add("title", "1.0", "1.end")
		txt.tag_configure("title", foreground="DarkRed", font=self.titleFont)
		txt.update_idletasks()

		# WHAT's & SDUM
		txt.tag_delete("what")
		for m in _patwhat.finditer(text):
			txt.tag_add("what",
					"1.0 + %d chars"%(m.start()),
					"1.0 + %d chars"%(m.end()))
		for m in _patsdum.finditer(text):
			txt.tag_add("what",
					"1.0 + %d chars"%(m.start()),
					"1.0 + %d chars"%(m.end()))
		txt.tag_configure("what", foreground="DarkMagenta")

		# Hyper links
		txt.tag_delete("card")
		tag = "card"
		for i in range(len(Manual.section)-1,-1,-1):
			s = Manual.section[i]
			if s == section: continue
			if s[0]>="A" and s[0]<="Z":
				prefix = 0
			else:
				prefix = 1
			pat1, pat2 = Manual.pattern[i]
			for m in pat1.finditer(text):
				txt.tag_add(tag,
						"1.0 + %d chars"%(m.start()+prefix),
						"1.0 + %d chars"%(m.end()))

			if pat2 is not None:
				for m in pat2.finditer(text):
					txt.tag_add(tag,
							"1.0 + %d chars"%(m.start()+prefix),
							"1.0 + %d chars"%(m.end()))
		txt.tag_remove(tag, "1.0", "2.0")

		# Highlight current section
		try:
			i = Manual.section.index(section)
			tag = "me"
			pat1, pat2 = Manual.pattern[i]
			for m in pat1.finditer(text):
				txt.tag_add(tag,
						"1.0 + %d chars"%(m.start()),
						"1.0 + %d chars"%(m.end()))

			if pat2 is not None:
				for m in pat2.finditer(text):
					txt.tag_add(tag,
							"1.0 + %d chars"%(m.start()),
							"1.0 + %d chars"%(m.end()))
			txt.tag_configure("me", foreground="DarkBlue")
		except ValueError:
			pass
			#if section == Manual.section[i]:
			#	tag = "me"
			#else:
		txt.tag_configure("card", foreground="Blue", underline=1)
		txt.tag_bind("card", "<Button-1>", self.hyperlink)

		# Highlight search string
		if self.searchStr is not None:
			txt.tag_delete("find")
			firstIndex = None
			for target in self.searchStr:
				index = "1.0"
				while True:
					index = str(txt.search(target, index,
							nocase=1, stopindex=END))
					if not index: break
					endex = '%s + %d chars' % (index, len(target))
					txt.tag_add("find", index, endex)
					txt.update_idletasks()
					if firstIndex is None:
						firstIndex = index
					index = endex
			if firstIndex is not None:
				txt.see(firstIndex)
			txt.tag_configure("find",
					background = '#800000',
					foreground = '#FFFFFF')

		# Find&insert all image tags
		index = "0.0"
		while True:
			index = str(txt.search("{img:", index))
			if len(index)==0: break
			indexEnd = txt.search("}", index)
			imageName = txt.get(index, indexEnd)[5:]
			txt.delete(index, "%s + 1 char"%(indexEnd))
			image = images.get(imageName)
			if image is None:
				# Try to load from icons directory
				if imageName.find('/') != -1:
					path =  imageName
				else :
					path = os.path.join(tkFlair.prgDir, "icons", imageName)

				try:
					image = PhotoImage(master=self.toplevel,
								file=path)
				except TclError:
					image = None
					say("WARNING Image '%s' not found"%(imageName))
				images[imageName] = image

			if image is not None:
				try:
					txt.image_create(index, image=image)
				except TclError:
					say("ERROR displaying image '%s'"%(imageName))

		# No more editing is allowed
		txt.config(state=DISABLED)
		return True

	# --------------------------------------------------------------------
	# Write HTML
	# --------------------------------------------------------------------
	def writeHTML(self, section=None):
		if section is None:
			section = Manual.topSection

		try:
			idx = Manual.section.index(section)
		except:
			# Try to find the best match
			idx = 0
			su = section[:8].upper()
			for s in Manual.section:
				if su == s[:8].upper():
					section = s
					break
				idx += 1

		# Search item in list to highlight
		try:
			title=Manual.title[idx]
		except IndexError:
			say("Section %s not found\n"%(section))
			return False

		# Add to view history
		if idx==0:
			backward = ''
		else:
			backward = '<a href="%s.html">Previous</a>'%(Manual.section[idx-1])

		if idx>=len(Manual.section)-1:
			foreward = ''
		else:
			foreward = '<a href="%s.html">Next</a>'%(Manual.section[idx+1])

		top = '<a href="%s.html">Index</a>'%(Manual.section[1])

#		text = Manual.manual.get(section,"").encode("utf-8") #.splitlines()
#		lf = text.find("\n".encode("utf-8"))
		text = Manual.manual.get(section,"") #.splitlines()
		lf = text.find("\n")
		if lf>0:
			firstline = "<h1>%s</h1>\n"%(text[:lf])
			text = text[lf+1:]
		else:
			firstline = ""

		# Special characters
		text = re.sub(r"&",'&amp;', text)
		text = re.sub(r"<",'&lt;',  text)
		text = re.sub(r">",'&gt;',  text)

		# WHAT's & SDUM
		text = _patwhat.sub(r'<div id="what">\g<0></div>', text)
		text = _patsdum.sub(r'<div id="what">\g<0></div>', text)
		text = _patcomment.sub(r'<div id="comment">\g<0></div>', text)

		text = re.sub(r"\bDefault\b",r'<div id="default">\g<0></div>', text)
		text = re.sub(r"\bHint:",   r'<div id="warning">\g<0></div>', text)
		text = re.sub(r"\bAdvice:", r'<div id="warning">\g<0></div>', text)
		text = re.sub(r"\bNote:",   r'<div id="default">\g<0></div>', text)
		text = re.sub(r"\bNotes:",  r'<div id="default">\g<0></div>', text)
		text = re.sub(r"\bExample:",r'<div id="default">\g<0></div>', text)
		text = re.sub(r"\bExample \d:",r'<div id="default">\g<0></div>', text)

		text = re.sub(r"\bWARNING:",r'<div id="warning">\g<0></div>', text)

		text = re.sub("1st", "1<sup>st</sup>", text)
		text = re.sub("2nd", "2<sup>nd</sup>", text)
		text = re.sub("3rd", "3<sup>rd</sup>", text)
		text = re.sub("4th", "4<sup>th</sup>", text)

		# Hyper links
		for i in range(len(Manual.section)-1,-1,-1):
			s = Manual.section[i]
			#if s == section: continue
			if s[0]>="A" and s[0]<="Z":
				prefix = 0
			else:
				prefix = 1
			pat1, pat2 = Manual.pattern[i]

			# loop both iterators
			iter1 = pat1.finditer(text)
			if pat2:
				iter2 = pat2.finditer(text)
			else:
				iter2 = None

			# Loop over all matches in the string
			out = ""
			last = 0

			# skip html syntax
			iter0 = _pathtml.finditer(text)
			try:    m0 = next(iter0)
			except: m0 = None
			try:    m1 = next(iter1)
			except: m1 = None
			try:    m2 = next(iter2)
			except: m2 = None

			while m1 is not None or m2 is not None:
				if (m1 and m2 and m1.start()<m2.start()) or m2 is None:
					if m0 is None or m1.start() < m0.start():
						out += '%s<a href="%s.html">%s</a>' % \
							(text[last:m1.start()],
							 m1.group(1),
							 m1.group(0))
						last = m1.end()
						try:    m1 = next(iter1)
						except: m1 = None
					else:
						out += text[last:m0.end()]
						last = m0.end()
						while m1 and m1.start()<m0.end():
							try:    m1 = next(iter1)
							except: m1 = None
						try:    m0 = next(iter0)
						except: m0 = None

				else:
					if m0 is None or m2.start() < m0.start():
						out += '%s<a href="%s.html">%s</a>' % \
							(text[last:m2.start()],
							 m2.group(1),
							 m2.group(0))
						last = m2.end()
						try:    m2 = next(iter2)
						except: m2 = None
					else:
						out += text[last:m0.end()]
						last = m0.end()
						while m2 and m2.start()<m0.end():
							try:    m2 = next(iter2)
							except: m2 = None
						try:    m0 = next(iter0)
						except: m0 = None

			text = out + text[last:]

		# write html
		say("Writing: "+section+".html")
		f = open(section+".html", "w")
		f.write("<!DOCTYPE html>\n")
		f.write("<html>\n")
		f.write("<head>\n")
		f.write("<style>\n")
		f.write("#what    { display:inline; color:maroon;     font-weight:bold;}\n")
		f.write("#default { display:inline; color:darkgreen;  font-weight:bold;}\n")
		f.write("#comment { display:inline; color:blue;}\n")
		f.write("</style>\n")
		f.write("<title>FLUKA: %s</title>\n"%(title))
		f.write("</head>\n")

		f.write("<body>\n")
		f.write("%s %s %s\n"%(backward, top, foreward))
		f.write("<p>\n")
		f.write(firstline)
		f.write("<pre>\n")
		f.write(text)
		f.write("</pre>\n")
		if idx==1:
			f.write("<h2>Complete Index</h2>\n")
			f.write("<ul>\n");
			prev = 0
			for i,s in enumerate(Manual.section):
				title = Manual.title[i]
				# Count front spaces
				for sp,ch in enumerate(title):
					if ch!=' ': break
				sp /= 3
				while sp>prev:
					f.write("<ul>\n");
					prev += 1
				while sp<prev:
					f.write("</ul>\n");
					prev -= 1
				f.write('<li><a href="%s.html">%s</a></li>\n'%(s,Manual.title[i]))
			f.write("</ul>\n");
			f.write("<p>\n")

		f.write("<p>\n")
		f.write("%s %s %s\n"%(backward, top, foreward))
		f.write("</body>\n")
		f.write("</html>\n")
		f.close()
		return True

	# --------------------------------------------------------------------
	# FindDialog
	# --------------------------------------------------------------------
	def findDialog(self, event=None):
		txt = self.entry.get()
		fd = tkDialogs.FindReplaceDialog(self.toplevel, replace=False)
		fd.show(self.find, None, None, txt)
		try:
			self.txt.tag_delete("find")
		except TclError:
			pass

	# ----------------------------------------------------------------------
	# Find target(s)
	# ----------------------------------------------------------------------
	def find(self, findStr=None, matchCase=None):
		if findStr is not None:
			self.searchStr = [findStr]
		else:
			if self.searchStr is None: return

		minrow = 100000
		mincol = 100000
		minlen = 0

		self.txt.tag_delete("find")
		for target in self.searchStr:
			index = str(self.txt.search(target,
					"insert + 1 chars",
					nocase = matchCase))
			if not index: continue
			row,col = list(map(int,index.split('.')))
			if row<minrow or (row==minrow and col<mincol):
				minrow = row
				mincol = col
				minlen = len(target)

		if minrow < 100000:
			index = "%d.%d" % (minrow,mincol)
			self.txt.mark_set('insert', index)
			self.txt.tag_add("find", 'insert', "insert + %d chars"%(minlen))
			self.txt.see(index)
			self.txt.tag_configure("find",
					background = '#800000',
					foreground = '#FFFFFF')

	# ----------------------------------------------------------------------
	def findNext(self, event=None):
		self.find()

	# ----------------------------------------------------------------------
	# hyper link
	# ----------------------------------------------------------------------
	def hyperlink(self, event):
#		row,col = list(map(int, string.split(
#			str(self.txt.index("@%d,%d"%(event.x,event.y))),".")))
		row,col = list(map(int, (str(self.txt.index("@%d,%d"%(event.x,event.y)))).split(".")))

		# Find line being clicked
		lines = Manual.manual[self.current_section].splitlines()
		line = lines[row-1]
		del lines

		# Find word
		if line[col]=='}': col -= 1
		begin = col
		while begin>=0:
			ch = line[begin]
			if ch.isalnum() or ch in ('-','.',':'):
				begin -= 1
			else:
				break
		begin += 1

		end = col
		while end<len(line):
			ch = line[end]
			if ch.isalnum() or ch in ('-','.',':'):
				end += 1
			else:
				break

		section = line[begin:end]
		if section[-1]==".": section = section[:-1]

		try:
			item = Manual.section.index(section)
		except:
			# Try to find the best match
			item = 0
			su = section[:8].upper()
			for s in Manual.section:
				if su == s[:8].upper():
					section = s
					break
				item += 1

		self.tocList.selection_clear(0,END)
		self.tocList.selection_set(item)
		self.tocList.activate(item)
		self.searchStr = None
		self.showSection(section)

	# ----------------------------------------------------------------------
	# load manual and split it into sections
	# ----------------------------------------------------------------------
	@staticmethod
	def load(file):
		try:
			f = open(file, "r")
		except IOError:
			messagebox.showerror("Error",
				sys.exc_info()[1])
			return
		(fn, ext) = os.path.splitext(os.path.basename(file))
		index     = True
		txt       = ""
		section   = "ABOUT-"+fn
		if Manual.topSection is None:
			Manual.topSection = section
		firstline = "----- %s -----" % (fn)
		inGeo     = False
		depth     = 1
		for line in f:
			if section=="":
				if index:
					section   = "INDEX-%s"%(fn)
					firstline = "Index %s"%(fn)
					index     = False
					txt       = ""
					continue

				line = line.strip()
				if line=="": continue
				for pat in _patterns:
					m = pat.match(line)
					if m is not None:
						section = m.group(1)
						if line[0]=="{":
							line = line.replace("{","")
							line = line.replace("}","")
							firstline = "      "+line
						else:
							firstline = line
						txt = line+"\n"
						inGeo = (section=="8")
						break

			elif line.startswith("1******************"):
				Manual._addSection(section, firstline, txt)
				section = ""

			elif line == "1\n":
				txt += "\n"
			else:
				# Check for a sub section
				for pat in _subpatterns:
					m = pat.match(line)
					if m is not None:
						if inGeo and _patbody2.match(line):
							depth = m.group(1).count(".")
							break
						Manual._addSection(section, firstline, txt)
						section   = m.group(1)
						depth     = section.count(".")
						firstline = ("   "*depth)+line.strip()
						txt = ""
						break

				if inGeo:
					m = _patbody1.match(line)
					if m is None:
						m = _patbody2.match(line)
					if m:
						Manual._addSection(section, firstline, txt)
						section   = m.group(1).replace(","," ")
						section   = section.replace("or","").split()
						firstline = ("   "*depth)+m.group(1)
						txt       = ""

					elif line==" Region data\n":
						Manual._addSection(section, firstline, txt)
						section   = "REGION"
						firstline = line.rstrip()
						txt       = ""

					elif line==" Region Volumes\n":
						Manual._addSection(section, firstline, txt)
						section   = line.rstrip()
						firstline = section
						txt       = ""

					elif line==" LATTICE card\n":
						Manual._addSection(section, firstline, txt)
						section   = "LATTICE"
						firstline = line.rstrip()
						txt       = ""

					elif line==" GEOEND card\n":
						Manual._addSection(section, firstline, txt)
						section   = "GEOEND card"
						firstline = "  "+section
						txt       = ""

					elif line=="1                             Voxel Geometry\n":
						Manual._addSection(section, firstline, txt)
						section   = "Voxel"
						firstline = " Voxel Geometry"
						txt       = ""
						line      = firstline+"\n"

				# Convert characters
				for src, dst in _change:
					try:
						line = line.replace(src,dst)
					except UnicodeDecodeError:
						pass
				try:
					txt += line
				except (UnicodeDecodeError, UnicodeEncodeError):
					txt += line.decode('utf-8','ignore')
		if section!="": Manual._addSection(section, firstline, txt)
		f.close()

	# ----------------------------------------------------------------------
	@staticmethod
	def _addSection(section, firstline, txt):
		if isinstance(section, list):
			# Special case for bodies
			for s in section:
				Manual.manual[s] = txt
				Manual.section.append(s)
				Manual.title.append("         "+s)
				pat1 = re.compile(r"\b(%s)\b"%(s[:8]))
				Manual.pattern.append((pat1, None))
			#section = section[0]
		else:
			section = section.strip()
			if len(Manual.manual)==0:
				Manual.manual["ABOUT"] = txt
			Manual.manual[section] = txt

			Manual.section.append(section)
			Manual.title.append(firstline)
			ch = section[0]
			if ch>='A' and ch<='Z':
				pat1 = re.compile(r"\b(%s)\b"%(section[:8]))
				f = firstline.split()[0]
				if f !=  section[:8]:
					pat2 = re.compile(r"\b(%s)%s\b"%(f[:8],f[8:]))
				else:
					pat2 = None
			elif ch==":":
				pat1 = re.compile(r"(%s)"%(section.replace(".","\.")))
				pat2 = None
			else:
				pat1 = re.compile(r" \b(%s)\}"%(section.replace(".","\.")))
				pat2 = None
			Manual.pattern.append((pat1, pat2))

	# ----------------------------------------------------------------------
	# Find latest manual file in the FLUKA directory
	# ----------------------------------------------------------------------
	@staticmethod
	def findManual(path):
		maxt   = 0
		manual_file = None
		if path is None or path=="": return
		for f in glob.glob("%s/*.manual"%(path)):
			try:
				t = os.stat(os.path.join(f))[ST_MTIME]
				if t > maxt: # get the newer only
					manual_file = f
					maxt = t
			except:
				pass
		return manual_file

#-------------------------------------------------------------------------------
# Show keyword
#-------------------------------------------------------------------------------
def show(keyword=None):
	global _directories, current

	if current is not None:
		current.clear()
#		current.showSection(keyword)
		current.toplevel.deiconify()
		current.toplevel.lift()
		current.toplevel.focus_set()
		manual = current
	else:
		if _directories is None:
			_directories = (".", os.getenv("FLUPRO"))
		manual = Manual(_directories)
		try:
			manual.show()
		except IndexError:
			messagebox.showwarning("Manual error",
				"Error opening manual")

	if isinstance(keyword,list):
		for k in keyword:
			if manual.showSection(k,False): break
	else:
		manual.showSection(keyword)

#-------------------------------------------------------------------------------
def destroy():
	if current is not None:
		current.close()

#-------------------------------------------------------------------------------
# Initialize Manual by settings the default search path
#-------------------------------------------------------------------------------
def init(dirs):
	global _directories
	_directories = dirs

#-------------------------------------------------------------------------------
def done():
	global images
	for i in list(images.values()):
		del i
	images = {}

#-------------------------------------------------------------------------------
def writeHTML(path):
	manual = Manual((path,))
	for s in manual.section:
		manual.writeHTML(s)

#-------------------------------------------------------------------------------
if __name__ == "__main__":
	import configparser
	tkFlair.openIni()

	var = "FLUPRO"
	flukadir = os.getenv(var)

	if len(sys.argv)>1 and sys.argv[1]=='html':
		writeHTML(flukadir)
		sys.exit()

	if flukadir == "" or flukadir is None:
		var = "FLUKA"
		flukadir = os.getenv(var)
		if flukadir == "" or flukadir is None:
			try:
				flukadir = tkFlair.config.get("Project","flukaDir")
			except:
				flukadir = "-"
	manual = Manual((".", flukadir))
	root = Tk()
	tkFlair.loadIcons()

	root.geometry("1024x800")
	manual.show(root)
	root.mainloop()
	tkFlair.delIcons()
	done()
