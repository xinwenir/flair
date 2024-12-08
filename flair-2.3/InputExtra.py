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
# Date:	13-Feb-2007

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import sys
import undo

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkFlair
import tkExtra
import InputPage
import Materials
import bFileDialog
import tkTextEditor

import Input
import Layout
import Manual

#===============================================================================
# Open filename dialog
#===============================================================================
filename = ""

# ------------------------------------------------------------------------------
def openOldFile(master, fn="", ft=None):
	global filename
	if ft is None: ft=("All","*")
	filename = bFileDialog.askopenfilename(master=master,
			title="Open file ReadOnly",
			initialfile=fn, filetypes=ft)
	return filename

# ------------------------------------------------------------------------------
def openNewFile(master, fn="", ft=None):
	global filename
	if ft is None: ft=("All","*")
	filename = bFileDialog.asksaveasfilename(master=master,
			title="Open file Read/Write",
			initialfile=fn, filetypes=ft)
	return filename

# ------------------------------------------------------------------------------
def openFile(master, fn="", ft=None):
	global filename
	if ft is None: ft=("All","*")
	filename = bFileDialog.askfilename(master=master,
			title="Open file Read/Write",
			initialfile=fn,
			filetypes=ft)
	return filename

#===============================================================================
# Error Dialog
#===============================================================================
class ErrorDialog(Toplevel):
	def __init__(self, master, card, errmsg):
		Toplevel.__init__(self, master)
		self.transient(master)
		self.title("Card Validation Errors")

		from InputPage import _FIXED_FONT

		if card.tag=="REGION":
			txt = Input.REGSCALE
		else:
			txt = Input.SCALE
		txt += str(card)

		lbl = Label(self, text=txt,
				font=_FIXED_FONT,
				justify=LEFT,
				anchor=NW,
				relief=SUNKEN)
		lbl.pack(side=TOP, fill=X, pady=3)

		msg = Label(self, text=errmsg,
				justify=LEFT,
				anchor=NW,
				relief=SUNKEN)
		msg.pack(side=TOP, expand=YES, fill=BOTH, pady=3)

		b = Button(self, text="Close", command=self.close)
		b.pack(side=BOTTOM)
		self.bind("<Escape>",   self.close)
		self.bind("<Return>",   self.close)
		self.bind("<KP_Enter>", self.close)

		self.deiconify()
		self.grab_set()
		self.wait_window()

	# ----------------------------------------------------------------------
	def close(self, event=None):
		self.destroy()

#===============================================================================
# Select List Popup
#===============================================================================
class SelectionList(Toplevel):
	def __init__(self, master):
		Toplevel.__init__(self, master, relief=RAISED)
		self.overrideredirect(1)
		self.transient(master)
		self.withdraw()

		self._isPosted = False
		self._commit   = True		# commit changes and move to next field
		#self._selectKeys = ""		# Keys that trigger the select&unpost

		# Create the listbox inside the dropdown window
		sb = Scrollbar(self, takefocus=0)
		self._listbox = tkExtra.SearchListbox(self,
					width=10,
					yscrollcommand=sb.set,
					selectmode=BROWSE,
					foreground=InputPage._VALUE_COLOR,
					font=InputPage._VALUE_FONT)
		sb.config(command=self._listbox.yview)

		self._listbox.pack(side=LEFT, expand=YES, fill=BOTH)
		sb.pack(side=RIGHT, fill=Y)

		# Bindings
		self.bind('<Escape>',	self._unselectUnpost)
		self.bind('<FocusOut>',	self._focusOut)
		self.bind('<Shift-Tab>',self._selectUnpostPrev)
		self.bind('<Tab>',	self._selectUnpostNext)
		if sys.platform == "linux2":
			self.bind('<ISO_Left_Tab>', self._selectUnpostPrev)
		self.bind('<Left>',	self._selectUnpostPrev)
		self.bind('<Right>',	self._selectUnpostNext)
		self.bind('<Return>',	self._selectUnpostNext)
		self.bind('<KP_Enter>',	self._selectUnpostNext)
		self._listbox.bind('<ButtonRelease-1>',self._buttonRelease)
		self._listbox.bind('<Button-3>',self._unselectUnpost)

	# ----------------------------------------------------------------------
	# listbox handling methods
	# ----------------------------------------------------------------------
	def post(self, geometry, lst=[], sel=None):
		if self._isPosted: return

		self._listbox.config(	font=InputPage._VALUE_FONT,
					foreground=InputPage._VALUE_COLOR)
		self._isPosted = True

		# Populate listbox
		self._list = lst
		self._listbox.delete(0, END)
		tkExtra.ExListbox.resetSearch()

		# Find best match
		selected = False

		# Convert selection to int if possible
		try:
			sf = float(sel)
			si = int(sf)
			if sf == float(si): sel = si
		except ValueError:
			pass

		try:
			sel8 = sel[:8]
		except TypeError:
			sel8 = None

		# Fill list
		for i,item in enumerate(lst):
			if isinstance(item, tuple):
				n, v = item
			else:
				n = v = item

			self._listbox.insert(END, str(n))
			if v==sel or v==sel8:
				self._listbox.selection_set(i)
				self._listbox.activate(i)
				self._listbox.see(i)
				selected = True

		if not selected:
#			self._listbox.selection_set(0)
#			self._listbox.activate(0)
			self._listbox.see(0)

		self.deiconify()
		self.wm_geometry(geometry)
		self.update_idletasks()
		self.lift()
		self.focus_set()
		self._listbox.focus_set()
		self.wait_visibility()
		# Check position
		i = 0
		while self.winfo_rootx()==0 and self.winfo_rooty()==0:
			self.wm_geometry(geometry)
			self.update_idletasks()
			i += 1
			if i>10: break

	# ----------------------------------------------------------------------
	def unPost(self, event=None):
		if not self._isPosted: return
		self._isPosted = False
		tkExtra.ExListbox.resetSearch()
		self.withdraw()

	# ----------------------------------------------------------------------
	def _selectUnpost(self, event=None):
		cur = self._listbox.index(ACTIVE)
		self._listbox.selection_clear(0, END)
		self._listbox.selection_set(cur)
		self.unPost()

	# ----------------------------------------------------------------------
	def _buttonRelease(self, event):
		self._listbox.selection_clear(0, END)
		i = self._listbox.nearest(event.y)
		self._listbox.activate(i)
		self._listbox.see(i)
		self._selectUnpost(event)
		if self._commit:
			self.master.commit()
			self.master.active.editItem()
		else:
			# Insert selection
			self.master.active.insertSelection()
			self.master.focus_set()

	# ----------------------------------------------------------------------
	def _selectUnpostPrev(self, event=None):
		self._selectUnpost(event)
		if self._commit:
			self.master.commit()
			self.master.active.editPrev()
		else:
			self.master.active.insertSelection(True)
			self.master.focus_set()

	# ----------------------------------------------------------------------
	def _selectUnpostNext(self, event=None):
		self._selectUnpost(event)
		if self._commit:
			self.master.commit()
			self.master.active.editNext()
		else:
			self.master.active.insertSelection()
			self.master.focus_set()

	# ----------------------------------------------------------------------
	def _unselectUnpost(self, event=None):
		self._listbox.selection_clear(0, END)
		self.unPost()
		#self.master.focus_set()
		if self._commit:
			self.master.commit()
			self.master.active.editItem()
		else:
			self.master.focus_set()

	# ----------------------------------------------------------------------
	def _focusOut(self, event):
		try:
			f = self.focus_get()
		except KeyError:
			pass
		else:
			if f == self or f == self._listbox: return
		self._listbox.selection_clear(0, END)
		self.unPost()

	# ----------------------------------------------------------------------
	# return selected item
	# ----------------------------------------------------------------------
	def selection(self):
		sel = self._listbox.curselection()
		if len(sel)==0: return (None, None)
		item = self._list[int(sel[0])]
		if isinstance(item, tuple):
			return item
		else:
			return (item, item)

	# ----------------------------------------------------------------------
	# return height of the frame
	# ----------------------------------------------------------------------
	def height(self):
		h = self.winfo_height()
		if h<10:
			self.deiconify()
			self.update_idletasks()
			h = self.winfo_height()
		return h

#===============================================================================
# Select groups dialog
#===============================================================================
class SelectGroupsDialog(Toplevel):
	# ----------------------------------------------------------------------
	def __init__(self, master, input):
		Toplevel.__init__(self, master)
		self.transient(master)
		self.title("Select cards to append")

		# --- Categories ---
		frame = LabelFrame(self, text="Categories to Append")
		frame.pack(side=TOP, expand=YES, fill=BOTH)

		self.categoryList = tkExtra.ExListbox(frame, borderwidth=0,
			selectmode=EXTENDED,
			selectborderwidth=0,
			background="White",
			takefocus=True,
			exportselection=FALSE)
		self.categoryList.pack(side=LEFT, expand=YES, fill=BOTH)

		sb = Scrollbar(frame, orient=VERTICAL,
				command=self.categoryList.yview)
		sb.pack(side=RIGHT, fill=Y)
		self.categoryList.config(yscrollcommand=sb.set)

		# --- Buttons ---
		frame = Frame(self)
		frame.pack(side=BOTTOM)
		b = Button(frame, text="Ok",
				image=tkFlair.icons["ok"],
				compound=LEFT,
				command=self.ok)
		b.pack(side=LEFT)
		b = Button(frame, text="Cancel",
				image=tkFlair.icons["x"],
				compound=LEFT,
				command=self.cancel)
		b.pack(side=RIGHT)

		# Fill categories
		self.categoryList.insert(END, Layout._ALL_GROUP)
		for group in Layout._groups_order:
			# Check if category is not empty
			cards = [ x for x in list(Input.CardInfo._db.values())
					if group in x.group ]
			empty = True
			for c in cards:
				if c.tag in input.cards:
					empty = False
					break
			if not empty:
				self.categoryList.insert(END, group)

		self.categoryList.selection_set(0)
		self.categoryList.activate(0)

		self.groups = None	# Empty list is all, None is none

	# ----------------------------------------------------------------------
	def show(self):
		self.wait_window()
		return self.groups

	# ----------------------------------------------------------------------
	def ok(self):
		items = list(map(int, self.categoryList.curselection()))
		if len(items)==0 or items[0] == Layout._ALL_GROUP:
			self.groups = []
		else:
			self.groups = [self.categoryList.get(i) for i in items]
		self.destroy()

	# ----------------------------------------------------------------------
	def cancel(self):
		self.groups = None
		self.destroy()

#===============================================================================
# Manual entry of card
#===============================================================================
class EditCardDialog(Toplevel):
	# ----------------------------------------------------------------------
	def __init__(self, canvas, flair, card):
		Toplevel.__init__(self, canvas)
		self.transient(canvas)
		self.card   = card
		self.canvas = canvas
		self.flair  = flair

		self.title("Edit Card")

		# --- Card Frame ---
		frame = Frame(self, relief=GROOVE, padx=3, pady=3)
		frame.pack(side=TOP, expand=YES, fill=BOTH, padx=3, pady=3)

		# ---
		row,col = 0, 0
		l = Label(frame, text="Comment:")
		l.grid(row=row, column=col, sticky=NE)

		col += 1
		self.comment = tkTextEditor.TextEditor(frame, background="white")
		self.comment.txt.config(height=5, foreground=InputPage._COMMENT_COLOR)
		self.comment.grid(row=row, column=col, columnspan=8, sticky=NSEW)

		# ---
		row,col = row+1, 0
		l = Label(frame, text="Card:")
		l.grid(row=row, column=col, sticky=E)

		col += 1
		self.tag = Label(frame, text=self.card.tag, anchor=W,
				relief=SUNKEN,
				foreground=InputPage._TAG_COLOR,
				background=tkFlair._ELABEL_BACKGROUND_COLOR)
		self.tag.bind('<Button-1>', self.tagButton)
		self.tag.bind('<Button-3>', self.tagButton)
		self.tag.grid(row=row, column=col, sticky=EW)

		col += 1
		b = Button(frame, text="\u25BC", padx=1, pady=0,
				command=self.tagButton)
		b.grid(row=row, column=col)

		# -
		col += 1
		l = Label(frame, text="Lines:")
		l.grid(row=row, column=col, sticky=E)

		col += 1
		self.linesVar = IntVar()
		lines = [x for x in range(1, 11) ]
		self.linesVar.trace("w", self.linesChange)
		l = OptionMenu(frame, self.linesVar, *lines)
		l.grid(row=row, column=col, columnspan=2, sticky=W)
		col += 1

		# -
		self.what = []

		col += 1
		l = Label(frame, text="sdum:")
		l.grid(row=row, column=col, sticky=E)

		col += 1
		cmd = lambda e=None,s=self,w=0:s.whatButton(w)
		e = Entry(frame,
				foreground=InputPage._VALUE_COLOR,
				background="White")
		e.bind('<Button-3>', cmd)
		e.grid(row=row, column=col, sticky=EW)

		col += 1
		b = Button(frame, text="\u25BC", padx=1, pady=0,
				command=cmd)
		b.grid(row=row, column=col)

		self.what.append((l, e, b, row, col-2))

		# --- Whats ---
		w = 0
		ln = 0
		for i in range(10):
			ln += 1
			for j in range(6):
				if j % 3 == 0:
					row,col = row+1, 0
				else:
					col += 1
				w += 1
				if w < 7:
					lbl = "%d:"%(j+1)
				else:
					lbl = "%d-%d:"%(ln, j+1)
				l = Label(frame, text=lbl)
				l.grid(row=row, column=col, sticky=E)

				col += 1
				cmd = lambda e=None,s=self,w=w:s.whatButton(w)

				e = Entry(frame,
						foreground=InputPage._VALUE_COLOR,
						background="White")
				e.bind('<Button-3>', cmd)
				e.grid(row=row, column=col, sticky=EW)

				col += 1
				b = Button(frame, text="\u25BC", padx=1, pady=0,
						command=cmd)
				b.grid(row=row, column=col)

				self.what.append((l, e, b, row, col-2))

		# ---
		row, col = row+1, 0
		l = Label(frame, text="Extra:")
		l.grid(row=row, column=col, sticky=NE)

		col += 1
		self.extra = tkTextEditor.TextEditor(frame, background="white")
		self.extra.txt.config(height=6, foreground=InputPage._VALUE_COLOR)
		self.extra.grid(row=row, column=col, columnspan=8, sticky=NSEW)
		tkExtra.Balloon.set(self.extra,"Extra information like title or region description of card")

		# ---
		row, col = row+1, 0
		l = Label(frame, text="Properties:")
		l.grid(row=row, column=col, sticky=NE)

		col += 1
		self.properties = tkTextEditor.TextEditor(frame, background="white")
		self.properties.txt.config(height=6, foreground=InputPage._VALUE_COLOR)
		self.properties.grid(row=row, column=col, columnspan=8, sticky=NSEW)
		tkExtra.Balloon.set(self.properties,"Additional variables/properties associated with card")

		frame.grid_rowconfigure(0, weight=1)
		frame.grid_rowconfigure(row, weight=1)
		frame.grid_columnconfigure(1, weight=1)
		frame.grid_columnconfigure(4, weight=1)
		frame.grid_columnconfigure(7, weight=1)

		# --- Buttons Frame ---
		frame = Frame(self)
		frame.pack(side=BOTTOM, fill=X)

		b = Button(frame, text="Help",
				image=tkFlair.icons["info"],
				compound=LEFT,
				command=self.help)
		b.pack(side=RIGHT)
		b = Button(frame, text="Cancel",
				image=tkFlair.icons["x"],
				compound=LEFT,
				command=self.cancel)
		b.pack(side=RIGHT)
		b = Button(frame, text="Ok",
				image=tkFlair.icons["ok"],
				compound=LEFT,
				command=self.ok)
		b.pack(side=RIGHT)
		self.rc = False

		self.bind('<Escape>', self.cancel)
		self.bind('<F1>', self.help)

	# ----------------------------------------------------------------------
	def fill(self):
		lines, rem = divmod(self.card.nwhats()-1, 6)
		if rem>0: lines += 1
		lines = max(1, lines)

		# fill comment
		self.comment.txt.delete(0.0, END)
		self.comment.txt.insert(0.0, self.card.comment())

		# set tag
		self.tag.config(text=self.card.tag)

		# line number
		self.linesVar.set(lines)

		# set whats
		w = 0
		for l,e,b,r,c in self.what[:lines*6+1]:
			e.delete(0, END)
			e.insert(0, self.card.what(w))
			w += 1

		# fill the extra
		self.extra.txt.delete(0.0, END)
		self.extra.txt.insert(0.0, self.card.extra())

		# fill the properties
		self.properties.txt.delete(0.0, END)
		if self.card.prop:
			props = ["%s=%s"%(n,str(v))
					for n,v in list(self.card.prop.items())
						if n[0]!="@"]
			self.properties.txt.insert(0.0, "\n".join(props))

	# ----------------------------------------------------------------------
	def linesChange(self, a, b, c):
		lines = self.linesVar.get()
		# display whats
		for l,e,b,r,c in self.what[:lines*6+1]:
			l.grid(row=r, column=c, sticky=E)
			e.grid(row=r, column=c+1, sticky=NSEW)
			b.grid(row=r, column=c+2)

		# remove the rest
		for l,e,b,r,c in self.what[lines*6+1:]:
			l.grid_forget()
			e.grid_forget()
			b.grid_forget()

	# ----------------------------------------------------------------------
	def show(self):
		self.fill()
		self.wait_visibility()
		self.grab_set()
		self.transient(self.canvas)
		self.wait_window()
		return self.rc

	# ----------------------------------------------------------------------
	# Add the information to the card and exit
	# ----------------------------------------------------------------------
	def ok(self, event=None):
		undoinfo = [(self.flair.refreshUndo,),
				self.flair.saveCardUndoInfo(self.card._pos),
				(self.flair.refreshUndo,)]
		# Collect the data

		# Change the tag if needed!
		if self.card.tag != self.tag["text"]:
			self.canvas.input.changeTag(self.card, self.tag["text"])

		# set the comment
		val = self.comment.txt.get(1.0, END).rstrip()
		try:	# Convert to unicode if needed
			val = str(val)
		except UnicodeEncodeError:
			val = val.encode("utf-8")
		self.card.setComment(val)

		# Set number of whats
		lines  = self.linesVar.get()

		# Find number of whats
		wend = lines*6
		while wend >= 0:
			l,e,b,r,c = self.what[wend]
			if e.get().strip() != "":
				break
			wend -= 1
		wend += 1

		self.card.setNWhats(wend)

		w = 0
		for l,e,b,r,c in self.what[:wend]:
			val = e.get().strip()
			try:	# Convert to ASCII
				val = str(val)
			except UnicodeEncodeError:
				val = val.encode("ascii","replace")
			self.card.setWhat(w, val)
			w += 1

		# Set extra
		val = self.extra.txt.get("1.0", END).rstrip()
		try:	# Convert to ASCII
			val = str(val)
		except UnicodeEncodeError:
			val = val.encode("ascii","replace")
		self.card.setExtra(val)

		# Set properties
		prop = self.properties.txt.get("1.0", END).splitlines()
		if prop:
			self.card.clearProperties(False)
			var = value = None
			for line in prop:
				line = line.strip()
				if line:
					try:
						equal = line.index("=")
						if var is not None: self.card[var] = value
						var   = line[:equal].strip()
						value = line[equal+1:].strip()
					except:
						value += line.strip()
			if var is not None: self.card[var] = value

		self.flair.addUndo(undo.createListUndo(undoinfo,"Edit card"))
		self.rc = True
		self.destroy()
		self.flair.setInputModified()

	# ----------------------------------------------------------------------
	def cancel(self, event=None):
		self.rc = False
		self.destroy()

	# ----------------------------------------------------------------------
	def help(self, event=None):
		Manual.show(self.card.tag)

	# ----------------------------------------------------------------------
	def tagButton(self, event=None):
		l,e,b,c,r = self.what[1]
		x = self.tag.winfo_rootx()
		y = self.tag.winfo_rooty() + self.tag.winfo_height()
		width = self.tag.winfo_width() + b.winfo_width()

		cards = list(Input.CardInfo._db.keys())
		cards.sort()
		tag = self.tag["text"]
		sel = SelectionListDialog(self).show(x, y, width, cards, tag)
		if sel:
			self.tag.config(text=sel)

	# ----------------------------------------------------------------------
	def whatButton(self, w):
		l,e,b,c,r = self.what[w]
		self.x = e.winfo_rootx()
		self.y = e.winfo_rooty() + e.winfo_height()
		self.w = w
		self.width = e.winfo_width() + b.winfo_width()

		# Show menu first
		menu = Menu(self, tearoff=0)
		menu.add_command(label='#define',  underline=0, command=self.defines)
		menu.add_command(label='Binning',  underline=1, command=self.binnings)
		menu.add_command(label='Body',     underline=0, command=self.bodies)
		menu.add_command(label='Card',     underline=0, command=self.cards)
		menu.add_command(label='Material', underline=0, command=self.materials)
		menu.add_command(label='Particle', underline=0, command=self.particles)
		menu.add_command(label='Region',   underline=0, command=self.regions)
		menu.add_command(label='Rot-defi', underline=1, command=self.rotdefis)
		menu.tk_popup(self.x, self.y)

	# ----------------------------------------------------------------------
	def defines(self):
		lst = self.canvas.layout.defineList()
		self.showList(lst)

	# ----------------------------------------------------------------------
	def binnings(self):
		lst = self.canvas.layout.binningList()
		self.showList(lst)

	# ----------------------------------------------------------------------
	def cards(self):
		cards = list(Input.CardInfo._db.keys())
		cards.sort()
		self.showList(cards)

	# ----------------------------------------------------------------------
	def bodies(self):
		lst = self.canvas.layout.bodyList()
		self.showList(lst)

	# ----------------------------------------------------------------------
	def materials(self):
		lst = self.canvas.layout.materialList()
		self.showList(lst)

	# ----------------------------------------------------------------------
	def particles(self):
		lst = Input.Particle.list
		self.showList(lst)

	# ----------------------------------------------------------------------
	def regions(self):
		lst = self.canvas.layout.regionList()
		self.showList(lst)

	# ----------------------------------------------------------------------
	def rotdefis(self):
		lst = self.canvas.layout.rotdefiList2()
		self.showList(lst)

	# ----------------------------------------------------------------------
	def showList(self, lst):
		sel = SelectionListDialog(self).show(
				self.x, self.y, self.width, lst,
				self.what[self.w][1].get())
		if sel:
			l,e,b,c,r = self.what[self.w]
			e.delete(0, END)
			e.insert(0, sel)

#===============================================================================
# Selection List Dialog
#===============================================================================
class SelectionListDialog(Toplevel):
	# ----------------------------------------------------------------------
	def __init__(self, master):
		Toplevel.__init__(self, master)
		self.transient(master)
		self.overrideredirect(1)
		self.title("")

		sb = Scrollbar(self, takefocus=0)
		sb.pack(side=RIGHT, fill=Y)
		self.listbox = tkExtra.ExListbox(
					self,
					width=10,
					yscrollcommand=sb.set,
					selectmode=BROWSE,
					takefocus=1,
					#background="White",
					font=InputPage._VALUE_FONT)
		self.listbox.pack(side=LEFT, expand=YES, fill=BOTH)
		sb.config(command=self.listbox.yview)

		self.bind('<Escape>',	self.cancel)
		self.bind('<FocusOut>', self.focusOut)
		self.bind('<Return>',	self.ok)
		self.bind('<KP_Enter>',	self.ok)
		self.listbox.bind('<ButtonRelease-1>',	self.ok)
		self.listbox.bind('<Double-1>',	self.ok)
		self.listbox.bind('<Button-3>',	self.cancel)

	# ----------------------------------------------------------------------
	def show(self, x, y, w, lst, val):
		self.ret = None
		self._grab_window = self.grab_current()
		if self._grab_window is not None:
			self._grab_window.grab_release()

		self.listbox.delete(0, END)
		for item in lst:
			self.listbox.insert(END, item)
			if val == item:
				self.listbox.selection_set(END)
				self.listbox.activate(END)
				self.listbox.see(END)

		h = 200
		self.wm_geometry("%dx%d+%d+%d" % (w,h,x,y))
		self.deiconify()
		#self.lift()
		#self.wait_visibility()
		self.listbox.focus_set()
		self.wait_window()
		try:
			if self._grab_window is not None:
				self._grab_window.grab_set()
			self.master.focus_set()
		except TclError:
			pass
		return self.ret

	# ----------------------------------------------------------------------
	# Add the information to the card and exit
	# ----------------------------------------------------------------------
	def ok(self, event=None):
		try:
			sel = int(self.listbox.curselection()[0])
			self.listbox.activate(sel)
		except:
			pass
		self.ret = self.listbox.get(ACTIVE)
		self.destroy()

	# ----------------------------------------------------------------------
	def focusOut(self, event=None):
		focus = self.focus_get()
		if focus is None: return
		self.cancel()

	# ----------------------------------------------------------------------
	def cancel(self, event=None):
		self.destroy()

# ------------------------------------------------------------------------------
# Filter menu
# ------------------------------------------------------------------------------
def filterMenu(page, menu, card, search=None, case=False):
	if search:
		if case:
			target = search
		else:
			target = search.upper()

		rc = True
		cmd=lambda f=('*',target),p=page:p.filterOn(f)
		menu.add_command(label="Any: "+search,
				image=tkFlair.icons["GLOBAL"], compound=LEFT,
				command=cmd)
		cmd=lambda f=('b',target),p=page:p.filterOn(f)
		menu.add_command(label="Body: "+search,
				image=tkFlair.icons["TRC"], compound=LEFT,
				command=cmd)
		cmd=lambda f=('r',target),p=page:p.filterOn(f)
		menu.add_command(label="Region: "+search,
				image=tkFlair.icons["REGION"], compound=LEFT,
				command=cmd)
		cmd=lambda f=('m',target),p=page:p.filterOn(f)
		menu.add_command(label="Material: "+search,
				image=tkFlair.icons["MATERIAL"], compound=LEFT,
				command=cmd)
		cmd=lambda f=('t',target),p=page:p.filterOn(f)
		menu.add_command(label="Transformation: "+search,
				image=tkFlair.icons["ROT-DEFI"], compound=LEFT,
				command=cmd)
		cmd=lambda f=('d',target),p=page:p.filterOn(f)
		menu.add_command(label="Detector: "+search,
				image=tkFlair.icons["USRBIN"], compound=LEFT,
			command=cmd)

		if target in Input.Particle.particles:
			cmd=lambda f=('p',target),p=page:p.filterOn(f)
			menu.add_command(
				label="Particle: "+search,
				image=tkFlair.icons["particle"], compound=LEFT,
				command=cmd)

		try:
			target = str(int(target))
			cmd=lambda f=('u',target),p=page:p.filterOn(f)
			menu.add_command(
				label="Unit: "+search,
				image=tkFlair.icons["new"], compound=LEFT,
				command=cmd)
		except ValueError:
			pass
	else:
		rc = False

	if card is None:
		pass

	elif card.tag[0] == "#":
		for w in card.whats():
			if w != "":
				cmd=lambda f=('*',w),p=page:p.filterOn(f)
				menu.add_command(label="Any: "+str(w),
						image=tkFlair.icons["GLOBAL"], compound=LEFT,
						command=cmd)

	elif card.isGeo() and len(card.tag)==3:
		w = card.sdum()
		cmd=lambda f=('b',w),p=page:p.filterOn(f)
		if rc: menu.add_separator()
		menu.add_command(label="Body: "+str(w),
				image=tkFlair.icons["TRC"], compound=LEFT,
				command=cmd)
		rc = True

	elif card.tag == "REGION":
		w = card.sdum()
		cmd=lambda f=('r',w),p=page:p.filterOn(f)
		if rc: menu.add_separator()
		menu.add_command(label="Region: "+str(w),
				image=tkFlair.icons["REGION"], compound=LEFT,
				command=cmd)
		rc = True
	elif card.tag == "MATERIAL":
		w = card.sdum()
		cmd=lambda f=('m',w),p=page:p.filterOn(f)
		if rc: menu.add_separator()
		menu.add_command(label="Material: "+str(w),
				image=tkFlair.icons["MATERIAL"], compound=LEFT,
				command=cmd)
		rc = True
	elif card.tag == "ROT-DEFI":
		w = card.sdum()
		cmd=lambda f=('t',w),p=page:p.filterOn(f)
		if rc: menu.add_separator()
		menu.add_command(label="Transformation: "+str(w),
				image=tkFlair.icons["ROT-DEFI"], compound=LEFT,
				command=cmd)
		rc = True
	else:
		case = card.case()
		if case is None: return
		range_ = card.info.range[case]
		if card.tag in Input.DETECTOR_TAGS:
			w = card.sdum()
			if w!="":
				cmd=lambda f=('d',w),p=page:p.filterOn(f)
				if rc: menu.add_separator()
				menu.add_command(label="Detector: "+str(w),
					image=tkFlair.icons["USRBIN"], compound=LEFT,
					command=cmd)
				rc = True

		prev = ""
		for i in range(len(range_)):
			r = range_[i]
			w = card.what(i)
			if w=="": continue
			# strip leading minus if any
			if isinstance(w,str) and w[0]=='-': w = w[1:]
			if   r=="mi" or r=="smi":
				cmd=lambda f=('m',w),p=page:p.filterOn(f)
				if rc and prev!=r:
					menu.add_separator()
					prev = r
				menu.add_command(
					label="Material: "+str(w),
					image=tkFlair.icons["MATERIAL"], compound=LEFT,
					command=cmd)
				rc = True
			elif r=="ri" or r=="sri":
				cmd=lambda f=('r',w),p=page:p.filterOn(f)
				if rc and prev!=r:
					menu.add_separator()
					prev = r
				menu.add_command(
					label="Region: "+str(w),
					image=tkFlair.icons["REGION"], compound=LEFT,
					command=cmd)
				rc = True
			elif r=="pi" or r=="spi":
				cmd=lambda f=('p',w),p=page:p.filterOn(f)
				if rc and prev!=r:
					menu.add_separator()
					prev = r
				menu.add_command(
					label="Particle: "+str(w),
					image=tkFlair.icons["particle"], compound=LEFT,
					command=cmd)
				rc = True
			elif r=="di" or r=="bi":
				cmd=lambda f=('d',w),p=page:p.filterOn(f)
				if rc and prev!=r:
					menu.add_separator()
					prev = r
				menu.add_command(
					label="Detector: "+str(w),
					image=tkFlair.icons["USRBIN"], compound=LEFT,
					command=cmd)
				rc = True
			elif r=="ti":
				if isinstance(w,str) and w[0]=="-": w=w[1:]
				cmd=lambda f=('t',w),p=page:p.filterOn(f)
				if rc and prev!=r:
					menu.add_separator()
					prev = r
				menu.add_command(
					label="Transformation: "+str(w),
					image=tkFlair.icons["ROT-DEFI"], compound=LEFT,
					command=cmd)
				rc = True
			elif r=="lu" or r=="slu":
				w = str(int(card.absWhat(i)))
				cmd=lambda f=('u',w),p=page:p.filterOn(f)
				if rc and prev!=r:
					menu.add_separator()
					prev = r
				menu.add_command(
					label="Unit: "+str(w),
					image=tkFlair.icons["new"], compound=LEFT,
					command=cmd)
				rc = True

	# Filter invalid cards
	if rc: menu.add_separator()
	menu.add_command(label="Invalid cards",
			image=tkFlair.icons["debug"], compound=LEFT,
			command=page.showInvalid)

	if InputPage._SHOW_PREPROCESSOR:
		icon = tkFlair.icons["undef"]
		txt = "Hide preprocessor"
	else:
		icon = tkFlair.icons["define"]
		txt = "Show preprocessor"

	menu.add_command(label=txt, image=icon, compound=LEFT,
			command=page.togglePreprocessor)
	return rc

# ------------------------------------------------------------------------------
# Material Selection Popup
# ------------------------------------------------------------------------------
class MaterialList(Toplevel):
	selected = None
	shown    = None
	def __init__(self, master, **kw):
		if MaterialList.shown: return
		Toplevel.__init__(self, master, **kw)
		self.overrideredirect(1)
		self.transient(master)

		sb = Scrollbar(self, takefocus=0)
		self.listbox = tkExtra.SearchListbox(self,
					selectmode=BROWSE,
					yscrollcommand=sb.set,
					width=30,
					height=20,
					exportselection=FALSE,
					takefocus=FALSE)
		self.listbox.pack(side=LEFT, fill=BOTH)
		sb.config(command=self.listbox.yview)
		sb.pack(side=RIGHT, fill=Y)

		lst = list(Materials.materials.keys())
		lst.sort()
		sel = 0
		for m in lst:
			self.listbox.insert(END, m)
			if MaterialList.selected == m:
				sel = self.listbox.size()
		MaterialList.selected = None
		self.listbox.selection_set(sel)

		# Bind events to the dropdown window.
		self.listbox.bind('<Escape>',	self.cancel)
		self.listbox.bind('<Return>',	self.ok)
		self.listbox.bind('<KP_Enter>',	self.ok)
#		self.listbox.bind('<Button-1>',self._click)
		self.listbox.bind('<ButtonRelease-1>', self.ok)
		self.listbox.bind('<FocusOut>',	self.focusOut)

	# ----------------------------------------------------------------------
	def show(self, x, y):
		if MaterialList.shown:
			MaterialList.shown.destroy()
			MaterialList.shown = None
			return
		try:
			self._focus = self.focus_get()
		except TclError:
			self._focus = None
		MaterialList.shown = self
		self.deiconify()
		self.geometry('+%d+%d' % (x,y))
		self.lift()
		self.focus_set()
		self.listbox.focus_set()
		self.wait_window()
		MaterialList.shown = None
		return MaterialList.selected

	# ----------------------------------------------------------------------
	def ok(self, event=None):
		sel = self.listbox.curselection()
		if sel:
			MaterialList.selected = self.listbox.get(sel[0])
		tkExtra.ExListbox.resetSearch()
		self.destroy()

	# ----------------------------------------------------------------------
	def focusOut(self, event=None):
		focus = self.focus_get()
		self.cancel()

	# ----------------------------------------------------------------------
	def cancel(self, event=None):
		MaterialList.selected = None
		tkExtra.ExListbox.resetSearch()
		if self._focus: self._focus.focus_set()
		self.destroy()
