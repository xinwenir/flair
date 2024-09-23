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
# Date:	21-Aug-2014

try:
	import pydicom as dicom
except ImportError:
	import dicom

from tkinter import *
import tkExtra
import sys

#===============================================================================
# Display the DICOM information in a browsable
#===============================================================================
class DicomBrowser(tkExtra.MultiListbox):
	_NONE  = "\u25CC"	# dotted circle
	_OPEN  = "\u25BC"	# triangle down
	_CLOSE = "\u25B6"	# triangle right
	_NEST  = "   "		# nesting spaces

	# ----------------------------------------------------------------------
	def __init__(self, master, dataset=None, **options):
		tkExtra.MultiListbox.__init__(self, master,
				(('Name',  32, None),
				 ('Value', 48, None),
				 ('Tag',   12, None),
				 ('VM',     4, None),
				 ('VR',     4, None)), **options)
		self.sortAssist = None
		self.lists[0].bind("<ButtonRelease-1>", self.buttonRelease)
		self.lists[0].bind("<Right>",           self.expand)
		self.lists[0].bind("<Left>",            self.collapse)
		#self.bindLeftRight()

		self._tags = []		# tags per line
		self._open = {}		# open tags
		self._dataset = dataset	# currently displayed dataset
		self._filter     = None
		self._filterCase = False

	# ----------------------------------------------------------------------
	def dataset(self, ds=None):
		if ds is not None: self._dataset = ds
		return self._dataset

	# ----------------------------------------------------------------------
	def clear(self):
		del self._tags[:]
		self._open.clear()
		self._dataset = None
		self.delete(0,END)

	# ----------------------------------------------------------------------
	# Recursively fill a DICOM file
	# ----------------------------------------------------------------------
	def _fill(self, root, nesttag, depth=0, prefix=""):
		for i,item in enumerate(root):
			if isinstance(item, dicom.dataelem.RawDataElement):
				name  = "RAW"
				value = item.value
				tag   = item.tag
				VM    = "-"
				VR    = item.VR
			elif isinstance(item, dicom.dataelem.DataElement):
				name  = item.description()
				value = item.value
				tag   = item.tag
				VM    = item.VM
				VR    = item.VR
			elif isinstance(item, dicom.dataset.Dataset):
				ntag = "%s:%s"%(nesttag,i)
				self._fill(item, ntag, depth, "%d. "%(i+1))
				continue
			else:
				raise Exception("Unknown item type: "+str(type(item)))

			ntag = nesttag + ":" + str(tag.numerator)
			if isinstance(value, dicom.sequence.Sequence):
				if self._open.get(ntag,0):
					handle = DicomBrowser._OPEN
				else:
					handle = DicomBrowser._CLOSE
				svalue = ""
			else:
				handle = DicomBrowser._NONE
				if VR=="OB":
					svalue = "Binary data length=%d"%(len(value))
				else:
					try:
						svalue = str(value)
					except UnicodeDecodeError:
						svalue = "Binary data length=%d"%(len(value))

			sname = "%s%s %s%s"%("   "*depth, handle, prefix, name)
			tag = str(tag)
			if self._filter:
				if self._filterCase:
					sn = sname.upper()
					sv = svalue.upper()
				else:
					sn = sname
					sv = svalue
				if self._filter in sn or \
				   self._filter in sv or \
				   self._filter in tag:
					self.insert(END, (sname, svalue, tag, VM, VR))
			else:
				self.insert(END, (sname, svalue, tag, VM, VR))
			self._tags.append(ntag)

			# If open sequence
			if handle == DicomBrowser._OPEN:
				self._fill(value, ntag, depth+1)

	# ----------------------------------------------------------------------
	# Fill a dicom
	# ----------------------------------------------------------------------
	def fill(self, dataset=None):
		yview = self.yview()
		del self._tags[:]
		self.delete(0,END)
		if dataset is not None: self._dataset = dataset
		if self._dataset is None: return
		self._fill(self._dataset,'')
		self.yview_moveto(yview[0])

	# ----------------------------------------------------------------------
	def _changeTag(self, tag):
		self.fill()
		try:
			idx = self._tags.index(tag)
			self.selection_set(idx)
			self.activate(idx)
			self.see(idx)
		except:
			pass

	# ----------------------------------------------------------------------
	# Open close a subelement when the button-1 is released
	# ----------------------------------------------------------------------
	def buttonRelease(self, event):
		if len(event.widget.curselection()) != 1: return
		if event.x > event.widget.winfo_width()/5: return
		i = event.widget.nearest(event.y)
		item = event.widget.get(i)
		handle = item.strip()[0]
		if handle == DicomBrowser._NONE: return
		tag = self._tags[i]
		if handle == DicomBrowser._CLOSE:
			self._open[tag] = True
		else:
			self._open[tag] = False
		self._changeTag(tag)

	# ----------------------------------------------------------------------
	# Expand current node
	# ----------------------------------------------------------------------
	def expand(self, event=None):
		i = self.lists[0].index(ACTIVE)
		item = event.widget.get(i)
		handle = item.strip()[0]
		if handle == DicomBrowser._NONE: return
		tag = self._tags[i]
		if handle == DicomBrowser._CLOSE:
			self._open[tag] = True
		else:
			return
		self._changeTag(tag)

	# ----------------------------------------------------------------------
	# Collapse current or master node
	# ----------------------------------------------------------------------
	def collapse(self, event=None):
		i = self.lists[0].index(ACTIVE)
		item = event.widget.get(i)
		handle = item.strip()[0]
		tag = self._tags[i]
		if handle==DicomBrowser._NONE or handle==DicomBrowser._CLOSE:
			if item[0] != ' ': return
			# We are in sub nesting find the head node and close it
			last = tag.rfind(":")
			if last==-1: return
			prev = tag.rfind(":",0,last)
			if prev==-1: return
			tag = tag[:prev]
			self._open[tag] = False
		else:
			if handle == DicomBrowser._CLOSE:
				return
			else:
				self._open[tag] = False
		self._changeTag(tag)

	# ----------------------------------------------------------------------
	# Search for a tag
	# @param txt	txt to search
	# @param upper	case insensitive search
	# @return None or tag nesting
	# ----------------------------------------------------------------------
	def filter(self, txt=None, case=False):
		if txt is None:
			self._filter = None
		else:
			if case:
				self._filter = txt.upper()
			else:
				self._filter = txt
			self._filterCase = case
		self.fill()

	# ----------------------------------------------------------------------
	# Search for a tag
	# @param txt	txt to search
	# @param upper	case insensitive search
	# @return None or tag nesting
	# ----------------------------------------------------------------------
	def search(self, txt, case=False):
		cur = self.curselection()
		if len(cur)==1:
			self._skip = self._tags[int(cur[0])]
		else:
			self._skip = None
		if case: txt = txt.upper()
		return self._search(txt, self._dataset, '', case)

	# ----------------------------------------------------------------------
	# Recursive function to perform the search
	# ----------------------------------------------------------------------
	def _search(self, needle, root, nesttag, upper):
		for i,item in enumerate(root):
			if isinstance(item, dicom.dataelem.RawDataElement):
				name  = "RAW"
				value = item.value
				tag   = item.tag
			elif isinstance(item, dicom.dataelem.DataElement):
				name  = item.description()
				value = item.value
				tag   = item.tag
			elif isinstance(item, dicom.dataset.Dataset):
				ntag = "%s:%s"%(nesttag,i)
				pos = self._search(needle, item, ntag, upper)
				if pos is not None: return pos
				continue
			else:
				raise Exception("Unknown item type: "+str(type(item)))

			ntag = nesttag + ":" + str(tag.numerator)
			if isinstance(value, dicom.sequence.Sequence) or \
			   item.VR=="OB":
				svalue = ""
			else:
				try:
					svalue = str(value)
				except UnicodeDecodeError:
					svalue = ""

			if self._skip is None:
				if upper:
					name   = name.upper()
					svalue = svalue.upper()
				if needle in name or \
				   needle in svalue:
					return ntag

			elif self._skip == ntag:
				self._skip = None

			# If open sequence
			if isinstance(value, dicom.sequence.Sequence):
				pos = self._search(needle, item, ntag, upper)
				if pos is not None: return pos
		return None

	# ----------------------------------------------------------------------
	# Open all tags and select the current one
	# ----------------------------------------------------------------------
	def selectTag(self, tag):
		# select to open all parent tags
		parent = ""
		refill = False
		for t in tag.split(":")[1:-1]:
			parent = "%s:%s"%(parent, t)
			if not refill and not self._open.get(parent,False):
				refill = True
			self._open[parent] = True

		if refill: self.fill()

		self.selection_clear(0,END)
		try:
			idx = self._tags.index(tag)
			self.selection_set(idx)
			self.activate(idx)
			self.see(idx)
		except:
			pass

#===============================================================================
if __name__ == "__main__":
	from log import say
	if len(sys.argv) != 2:
		say("syntax: viewer.py <dicom_file>")
		sys.exit()

	try:
		ds = dicom.read_file(sys.argv[1],
				force=True,
				defer_size=256,
				stop_before_pixels=True)
		ds.decode()
	except:
		say("Error opening file:", sys.argv[1])
		sys.exit()

	tk = Tk()
	tk.title("DICOM Viewer: "+sys.argv[1])
	tk.bind("<Control-q>", lambda e : tk.destroy())

	tree = DicomBrowser(tk, ds, height=30)
	tree.pack(expand=TRUE, fill=BOTH)
	entry = Entry(tk)
	entry.pack(side=BOTTOM, fill=X)

	def searchTree(event):
		upper = True
		needle = event.widget.get()
		tag = tree.search(needle, upper)
		if tag:
			say("Found:", tag)
			tree.selectTag(tag)

		else:
			say("Not found")
	entry.bind("<Return>", searchTree)

	tree.fill()
	tk.mainloop()
