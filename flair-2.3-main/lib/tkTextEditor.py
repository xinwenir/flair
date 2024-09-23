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
# Date:	24-Jun-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

try:
	from tkinter import *
	from tkinter import _cnfmerge
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	from tkinter import _cnfmerge
	import tkinter.messagebox

import tkExtra
import tkDialogs
from log import say

#=============================================================================
# TextEditor
#=============================================================================
class TextEditor(Frame):
	def __init__(self, master, *arg, **kw_args):
		Frame.__init__(self, master)
		self.txt = Text(self, *arg, **kw_args)
		self.txt.grid(row=0, column=0, sticky='nsew')
		sbv = tkExtra.AutoScrollbar(self,
				takefocus=False)
		sbv.grid(row=0, column=1, sticky='ns')
		sbh = tkExtra.AutoScrollbar(self,
				takefocus=False,
				orient=HORIZONTAL)
		sbh.grid(row=1, column=0, sticky='ew')

		self.txt.config(  yscrollcommand=sbv.set,
				  xscrollcommand=sbh.set )
		sbv.config(command=self.txt.yview)
		sbh.config(command=self.txt.xview)

		self.grid_rowconfigure(0, weight=1)
		self.grid_columnconfigure(0, weight=1)

	# ----------------------------------------------------------------------
	def cget(self,item):
		if not hasattr(self.txt,item):
			return Frame.cget(self.txt,item)
		else:
			getattr(self.txt,item)
	__getitem__ = cget

	# ----------------------------------------------------------------------
	def configure(self,cnf=None,**kw):
		if kw: cnf=_cnfmerge((cnf,kw))
		for key in cnf.keys():
			if not hasattr(self.txt,key):
				self.txt.configure(cnf)
			else:
				setattr(self.txt,key,cnf[key])
	config=configure

	# ----------------------------------------------------------------------
	# Use this method to get the parent widget of the frame
	# ----------------------------------------------------------------------
	def __call__(self): return self.txt

	# --------------------------------------------------------------------
	def findDialog(self, event=None):
		try: txt = self.txt.selection_get()
		except:
			txt = ""
		fd = tkDialogs.FindReplaceDialog(self, replace=False)
		fd.show(self.find, None, None, txt)
		try:
			self.txt.tag_delete("find")
		except TclError:
			pass

	# --------------------------------------------------------------------
	def replaceDialog(self, event=None):
		try: txt = self.txt.selection_get()
		except: txt = ""
		fd = tkDialogs.FindReplaceDialog(self,replace=True)
		fd.show(self.find, self.replace, self.replaceAll, txt)
		self.txt.tag_delete("find")

	# --------------------------------------------------------------------
	def find(self, findStr=None, matchCase=None):
		global save_findStr, save_matchCase
		if findStr:
			save_findStr   = findStr
			save_matchCase = matchCase
		else:
			try:
				findStr   = save_findStr
				matchCase = save_matchCase
			except: return

		self.txt.tag_remove(SEL,"0.0",END)
		self.txt.tag_delete("find")
		try:
			index = self.txt.search(findStr,
				"insert + 1 chars",
				nocase = matchCase)
			self.txt.mark_set('insert', index)
			endex = '%s + %d chars' % ('insert', len(findStr))
			self.txt.tag_add(SEL, index, endex)
			self.txt.tag_add("find", index, endex)
			self.txt.tag_configure("find",
					background = '#800000',
					foreground = '#FFFFFF')
			self.txt.see(index)
			self.txt.update_idletasks()
			self.found = True
			return True
		except:
			messagebox.showinfo("Not Found",
				"Target \"%s\" not found" % (findStr))
			self.found = False
			return False

	# --------------------------------------------------------------------
	def findNext(self, event=None):
		self.find()

	# --------------------------------------------------------------------
	def replace(self, findStr, replaceStr, matchCase):
		if not self.found:
			if not self.find(findStr, matchCase): return
		index = self.txt.index(SEL_FIRST)
		self.txt.delete(SEL_FIRST, SEL_LAST)
		self.txt.insert(index, replaceStr)
		endex = "%s + %d chars" % (index, len(replaceStr))
		self.txt.tag_add(SEL, index, endex)
		self.txt.tag_add("find", index, endex)
		self.found = False

	# --------------------------------------------------------------------
	def replaceAll(self, findStr, replaceStr, matchCase):
		say("TextEditor::replaceAll", findStr, replaceStr, matchCase)

	# ----------------------------------------------------------------------
	def delete(self, index1, index2=None): return self.txt.delete(index1, index2)
	def insert(self, index, chars, *args): return self.txt.insert(index, chars, *args)

#=============================================================================
# Testing...
if __name__ == "__main__":
	root = Tk()
	editor = TextEditor(root)
	editor.pack(side=LEFT,expand=1,fill=BOTH)
	root.bind_class(root,'<Control-Key-q>',lambda e,r=root: r.destroy())
	root.bind_class(root,'<Control-Key-f>',editor.findDialog)
	root.bind_class(root,'<Control-Key-r>',editor.replaceDialog)
	root.bind_class(root,'<F3>',editor.findNext)

	root.mainloop()
