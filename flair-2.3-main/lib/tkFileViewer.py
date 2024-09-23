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
#
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	14-Sep-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import time
from stat import *

try:
	from tkinter import *
	from tkinter.colorchooser import askcolor
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	import tkinter.messagebox

import tkTree
import tkExtra
import tkDialogs
import bFileDialog
import tkTextEditor

_DEFAULT_SPLIT = 150

#===============================================================================
# OutputViewer Node
#===============================================================================
class OVNode(tkTree.Node):
	def __init__(self, *args, **kw_args):
		tkTree.Node.__init__(self, *args, **kw_args)
		self.widget.tag_bind(self.label,  '<1>', self.select)
		self.widget.tag_bind(self.symbol, '<1>', self.selectImage)
		self.widget.tag_bind(self.label,  '<3>', self.popupMenu)

	# ----------------------------------------------------------------------
	def select(self, event=None):
		tkTree.Node.select(self)
		self.PVT_enter(event)
		self.widget.move_cursor(self)
		self.showSection()

	# ----------------------------------------------------------------------
	def selectImage(self, event=None):
		self.PVT_click(event)
		self.select(event)

	# ----------------------------------------------------------------------
	def showSection(self):
		self.widget.viewer.showSection(self.full_id())
		self.widget.viewer.splitter.rightFrame().focus_set()

	# ----------------------------------------------------------------------
	def popupMenu(self, event=None):
		self.select()
		menu=Menu(self.widget, tearoff=0)
		menu.add_command(label='Reload', command=self.widget.viewer.reloadSelected)
		menu.add_command(label='Close', command=self.widget.viewer.close)
		menu.tk_popup(event.x_root, event.y_root)

#===============================================================================
# File Data Viewer
#===============================================================================
class FileData:
	def __init__(self, filename, maxsize=None):
		self.filename = filename
		self.fn, self.ext = os.path.splitext(filename)
		self.maxsize  = maxsize
		self.filedata = ""
		self.mtime    = 0
		self.ypos     = 0.0

		# Sections contains a dictionary with the range in the file
		# to display: self.sections[idname] = (from, to)
		self.sections = {}
		# Tree contains a nested list of lists to construct the
		# displayed tree
		# [ (title, idname, None|another-tree),... ]
		self.tree     = []

	# ----------------------------------------------------------------------
	# Returns True on failure
	# ----------------------------------------------------------------------
	def load(self):
		try:
			stat = os.stat(self.filename)
			size = stat[ST_SIZE]
			self.mtime = stat[ST_MTIME]
			if self.maxsize is not None:
				if size >= self.maxsize: return True
			f = open(self.filename, "r")
		except (OSError, IOError):
			return True
		try:
			self.filedata = f.read()
		except (OverflowError, MemoryError):
			return True
		f.close()
		return False

	# ----------------------------------------------------------------------
	def set(self, text):
		self.filedata = text

	# ----------------------------------------------------------------------
	# Override routine for creating the sections table
	# Sections dictionary should contain tuples with the range start:end
	# ----------------------------------------------------------------------
	def makeTOC(self):
		pass

	# ----------------------------------------------------------------------
	# Return icon to display
	# ----------------------------------------------------------------------
	def icon(self, name):
		return None

	# ----------------------------------------------------------------------
	# Override routine to highlight messages
	# @param widget the text widget to highlight
	# ----------------------------------------------------------------------
	def highlight(self, widget):
		pass

#===============================================================================
# File Viewer frames
#===============================================================================
class FileViewer(Frame):
	split  = _DEFAULT_SPLIT

	def __init__(self, master, **kw):
		"""Open a viewer assuming that the file exists."""
		Frame.__init__(self, master, **kw)

		self._inFocus = False
		self.splitter  = tkExtra.HSplitter(self, FileViewer.split, True)
		self.splitter.pack(expand=YES, fill=BOTH)

		self.treeFrame = self.splitter.leftFrame()
		self.editorFrame = self.splitter.rightFrame()

		# Tree
		self.tree = tkTree.Tree(master=self.treeFrame,
			    root_id="/",
			    root_label="Viewer",
			    get_contents_callback=self.treeContents,
			    width=150,
			    node_class=OVNode
			    )
		self.tree.grid(row=0, column=0, sticky=NSEW)
		self.tree.viewer = self

		self.treeFrame.grid_rowconfigure(0, weight=1)
		self.treeFrame.grid_columnconfigure(0, weight=1)

		# Scroll bars
		sb=tkExtra.AutoScrollbar(self.treeFrame, takefocus=0)
		sb.grid(row=0, column=1, sticky=NS)
		self.tree.configure(yscrollcommand=sb.set)
		sb.configure(command=self.tree.yview)

		sb=tkExtra.AutoScrollbar(self.treeFrame,
				orient=HORIZONTAL, takefocus=0)
		sb.grid(row=1, column=0, sticky=EW)
		self.tree.configure(xscrollcommand=sb.set)
		sb.configure(command=self.tree.xview)

		self.tree.focus_set()

		self.editor = tkTextEditor.TextEditor(self.editorFrame)
		self.editor.txt.config(state=DISABLED)
		self.editor.txt.focus_set()
		self.editor.pack(side=LEFT, expand=YES, fill=BOTH)

		self.editor.txt.bind('<Up>',   lambda e, s=self.editor.txt:
					s.yview(SCROLL, -1, UNITS))
		self.editor.txt.bind('<Down>', lambda e, s=self.editor.txt:
					s.yview(SCROLL,  1, UNITS))
		self.editor.txt.bind('<Button-1>', lambda e, s=self.editor.txt:
					s.focus_set())
		self.editor.txt.bind("<F11>",	self.toggleSplitter)
		self.tree.bind("<F11>",		self.toggleSplitter)
		self.bind("<F11>",	self.toggleSplitter)
		self.bind('<FocusIn>', self._focusIn)

		self.files = []
		self.curid = None

	# ----------------------------------------------------------------------
	def _append(self, curfile):
		self.files.append(curfile)

		#self.tree.root.set_label(curfile.filename)
		self.tree.root.set_label("Files")
		if self.tree.root.expanded():
			self.tree.root.collapse()
		self.tree.root.expand()

		node = self.tree.root.child_nodes[-1]
		if node is not None:
			node.select()
			if node.expandable(): node.expand()

	# ----------------------------------------------------------------------
	# Routine to be overwritten
	# @return the FileData structure
	# ----------------------------------------------------------------------
	def loadData(self, filename, maxsize):
		curfile = FileData(filename, maxsize)
		if curfile.load():
			return None
		curfile.makeTOC()
		return curfile

	# ----------------------------------------------------------------------
	# Load file, returns True on failure
	# ----------------------------------------------------------------------
	def load(self, filename, maxsize=None):
		curfile = self.loadData(filename, maxsize)
		if curfile is None: return True
		self._append(curfile)

		return False

	# ----------------------------------------------------------------------
	# Show a text(file) from memory
	# ----------------------------------------------------------------------
	def show(self, name, text):
		curfile = FileData(name)
		curfile.set(text)
		self._append(curfile)

	# ----------------------------------------------------------------------
	def hasFile(self, filename):
		return filename in [x.filename for x in self.files if x]

	# ----------------------------------------------------------------------
	def reload(self, filename):
		for i,f in enumerate(self.files):
			if f is not None and f.filename == filename:
				node = self.tree.find_full_id_expand(["/",i])
				break
		if node is not None:
			node.select()
			self.reloadSelected()

	# ----------------------------------------------------------------------
	def reloadSelected(self, event=None):
		selid = self.tree.cursor_node().full_id()
		if len(selid)==1: return
		curfile = self.files[selid[1]]
		curfile.load()
		self.showSection(selid)

	# ----------------------------------------------------------------------
	def openFile(self, event=None):
		fn = bFileDialog.askopenfilename(master=self)
		if fn:
			self.load(fn)

	# ----------------------------------------------------------------------
	def close(self, event=None):
		if len(self.files)<1:
			return True
		else:
			selid = self.tree.cursor_node().full_id()
			if len(selid)==1: return
			if len(selid)>2:
				node = self.tree.find_full_id_expand(selid[:2])
			else:
				node = self.tree.find_full_id_expand(selid)
			del self.files[selid[1]]
			self.curid = None

			if self.tree.root.expanded():
				self.tree.root.collapse()
			self.tree.root.expand()

			return False

	# ----------------------------------------------------------------------
	def _focusIn(self, event=None):
		if self._inFocus: return
		# FocusIn is generated for all sub-windows, handle only the main window
		if self is not event.widget: return
		self._inFocus = True

		# Check time stamp for all files
		changed = []
		fn      = []
		for i,f in enumerate(self.files):
			if f is None: continue
			try:
				mtime = os.stat(f.filename)[ST_MTIME]
				if mtime > f.mtime:
					changed.append(i)
					fn.append(f.filename)
			except OSError:
				pass

		if changed:
			if messagebox.askyesno("File has changed",
					"Files:\n%s\nhave changed. Reload?"%("\n".join(fn)),
					parent=self):
				for i in changed:
					try:
						node = self.tree.root.child_nodes[i]
						self.tree.move_cursor(node)
					except IndexError:
						pass
					self.reloadSelected()
			else:
				for i in changed:
					f = self.files[i]
					try:
						f.mtime = os.stat(f.filename)[ST_MTIME]
					except OSError:
						pass

		self._inFocus = False

	# ----------------------------------------------------------------------
	def toggleTree(self, event=None):
		if self.splitter.split > self.splitter.min:
			self.splitter.toggle()
			self.splitter.rightFrame().focus_set()
			self.focus_set()
		else:
			self.splitter.toggle()
			if self.splitter.split <= self.splitter.min:
				self.splitter.split = _DEFAULT_SPLIT
				self.splitter.place()

	# ----------------------------------------------------------------------
	# Callback routine to fill tree contents
	# ----------------------------------------------------------------------
	def treeContents(self,node):
		nid = node.full_id()
		if len(nid)==1:
			for idx,f in enumerate(self.files):
				if f is not None:
					self.tree.add_node(name=f.filename,
						id=idx,
						flag=len(self.files[idx].sections)>0)
			return

		treelist = None
		for nid in node.full_id():
			if nid == "/":
				treelist = None

			elif isinstance(nid, int):
				curfile  = self.files[nid]
				treelist = curfile.tree

			else:
				# scan list for target
				for name,id,sublist in treelist:
					if nid==id:
						treelist = sublist
						break
				else:
					treelist = None

		if treelist is None: return

		nid = node.full_id()
		for name,idx,sublist in treelist:
			icon = curfile.icon(idx)
			if sublist != None or idx in curfile.sections:
				 self.tree.add_node(name=name,
					id=idx,
					flag=(sublist!=None),
					expanded_icon=icon,
					collapsed_icon=icon)

	# ----------------------------------------------------------------------
	def showSection(self, fullid):
		self.editor.txt.config(state=NORMAL)

		# Remember position
		if isinstance(self.curid, int):
			self.files[self.curid].ypos = self.editor.txt.yview()[0]

		# Delete everything
		self.editor.txt.delete(1.0, END)

		# find section infos
		found = 0
		if len(fullid)==1:
			# Root node don't display anything
			self.curid = None
			if isinstance(self.master,Toplevel):
				self.title("FileViewer")
			pass

		else:
			self.curid = fullid[1]
			curfile = self.files[self.curid]
			if curfile is None: return
			if isinstance(self.master,Toplevel):
				self.master.title("FileViewer: %s"%(curfile.filename))

			if len(fullid)==2:
				self.editor.txt.insert(1.0, curfile.filedata)
				self.editor.txt.yview_moveto(curfile.ypos)
				curfile.highlight(self.editor.txt)
			else:
				section = curfile.sections.get(fullid[-1])
				if section is not None:
					self.editor.txt.insert(INSERT,
						curfile.filedata[section[0]:section[1]])
					curfile.highlight(self.editor.txt)
				self.curid = None	# forget it

		self.editor.txt.config(state=DISABLED)

	# ----------------------------------------------------------------------
	def clear(self):
		self.editor.txt.config(state=NORMAL)
		self.editor.txt.delete(1.0, END)
		self.editor.txt.config(state=DISABLED)

	# ----------------------------------------------------------------------
	def toggleSplitter(self, event=None):
		self.splitter.toggle()
