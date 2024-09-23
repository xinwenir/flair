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
# Date:	16-Aug-2012

import os
from tkinter import *
import time

import tkFlair
import tkExtra
import tkDialogs
import bFileDialog

import Dicom

import Ribbon
import FlairRibbon

#===============================================================================
# Dicom Tab
#===============================================================================
class DicomTab(FlairRibbon.FlairTab):
	"""Generate Voxel and USRBINs from DICOM files CT, PET, RT..."""
	#----------------------------------------------------------------------
	def __init__(self, master, page, **kw):
		FlairRibbon.FlairTab.__init__(self, master, page, **kw)

		# ---
		pane = PanedWindow(self, orient=VERTICAL, opaqueresize=1)
		pane.pack(expand=YES, fill=BOTH)

		frame = Frame(pane, border=0)
		try: pane.add(frame, minsize=16, stretch="always")
		except: pane.add(frame, minsize=16)	# tk8.4

		# -------------- DataSets ---------------
		lblframe = LabelFrame(frame, text="Data sets", foreground="DarkBlue")
		lblframe.grid(row=1, column=0, columnspan=3, sticky=NSEW)

		self.dicomList = tkExtra.MultiListbox(lblframe,
				(('Series Instance UID', 30, None),
				 ('Modality',             8, None),
				 ('Date',                12, None),
				 ('Patient',             20, None),
				 ('Age',                  4, None),
				 ('Sex',                  3, None),
				 ('Institution',         10, None),
				 ('Study',               15, None)),
				height=5)
		self.dicomList.pack(expand=TRUE, fill=BOTH)
		self.dicomList.bind("<<ListboxSelect>>", self.datasetChange)
		self.dicomList.bindList("<Double-1>", self.datasetProcess)

		frame.grid_columnconfigure(1, weight=1)
		frame.grid_rowconfigure(1, weight=1)

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairTab.createRibbon(self)

		# ========== Set ===========
		group = Ribbon.LabelGroup(self.ribbon, "Set?")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		group.grid3rows()

		# ---
		col,row=0,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["add32"],
					text="Add",
					command=self.add,
					compound=TOP,
					background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=3, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Add a dicom set to project")

		# ---
		col,row = 1,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Delete",
				compound=LEFT,
				command=self.delete,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Delete dicom set")

		# ---
		col,row = 1,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["change"],
				text="Change",
				compound=LEFT,
				command=self.change,
				anchor=W,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Change dicom data set")

		# ---
		col,row = 1,2
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["clone"],
				text="Clone",
				compound=LEFT,
				anchor=W,
				command=self.clone,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=EW)
		tkExtra.Balloon.set(b, "Clone dicom set")

		self._lastHash = 0

	# ----------------------------------------------------------------------
	def activate(self):
		pass

	# ----------------------------------------------------------------------
	def get(self):
		pass

	# ----------------------------------------------------------------------
	def refresh(self):
		# Refresh the display only if the dicoms have changed
		if self._lastHash == self._dicomsHash(): return
		self.fill()
		self._lastHash = self._dicomsHash()

	# ----------------------------------------------------------------------
	def _dicomsHash(self):
		h = 0
		for di in self.project.dicoms:
			h += hash(di["uid"])
		return h

	# ----------------------------------------------------------------------
	# Called every time a dataset is changed
	# ----------------------------------------------------------------------
	def datasetChange(self, event=None):
		self.page.datasetChange()

	# ----------------------------------------------------------------------
	def datasetProcess(self, event=None):
		self.page.changeTab(1)

	# ----------------------------------------------------------------------
	# Add a new dicom set
	# ----------------------------------------------------------------------
	def add(self, event=None):
		result = DicomSelectDialog(self, self.project).show()
		if result is None: return

		path = result.pop(0)
		first = None
		for uid in result:
			s = self._addDicomInfo(path, uid)
			if first is None: first = s
		if first:
			self.dicomList.activate(first)
			self.dicomList.see(first)
		self.page.refresh()

	# ----------------------------------------------------------------------
	# Add a new dicom dataset
	# ----------------------------------------------------------------------
	def _addDicomInfo(self, path, uid):
		di = self.page.addDicomInfo(path, uid)
		s = self.dicomList.size()
		self.insertDicom(END, di)
		self.dicomList.selection_set(s)
		self._lastHash += hash(uid)
		self.page.refresh()
		return s

	# ----------------------------------------------------------------------
	def change(self, event=None):
		for i in self.dicomList.curselection():
			uid = self.dicomList.get(i)[0]
			di  = self.project.dicomInfo(uid)
			result = DicomSelectDialog(self, self.project, di["directory"], uid).show()

			if result is None: continue

			path = result.pop(0)
			for uid in result:
				if di is None:
					self._addDicomInfo(path, uid)
				else:
					# remove old hash
					self._lastHash -= hash(di["uid"])

					di["uid"]       = uid
					di["directory"] = path
					di["frames"]    = Dicom.Scan.frames(uid)
					di.files        = Dicom.Scan.files(uid)
					di.modality     = Dicom.Scan.modality(uid)

					filename = os.path.join(path, di.files[0][0])
					header = Dicom.dicomHeader(uid, filename)
					self.dicomList.set(i, header)
					self._lastHash += hash(uid)

	# ----------------------------------------------------------------------
	# Delete selected data sets
	# ----------------------------------------------------------------------
	def delete(self, event=None):
		sel = list(self.dicomList.curselection())
		sel.reverse()
		for i in sel:
			uid = self.dicomList.get(i)[0]
			self.dicomList.delete(i)
			for j,di in enumerate(self.project.dicoms):
				if di["uid"] == uid:
					del self.project.dicoms[j]
					break
		self.page.refresh()

	# ----------------------------------------------------------------------
	# Clone selected data sets
	# ----------------------------------------------------------------------
	def clone(self, event=None):
		sel = self.dicomList.curselection()
		self.dicomList.selection_clear(0,END)
		for i in sel:
			header = self.dicomList.get(i)
			di = self.project.dicomInfo(header[0])
			self.project.dicoms.append(di.clone())
			self.dicomList.insert(END, header)
			self.dicomList.selection_set(END)
			self._lastHash += hash(header[0])

	# ----------------------------------------------------------------------
	# Fill list with dicom sets
	# ----------------------------------------------------------------------
	def fill(self):
		self.dicomList.delete(0,END)
		i = 0
		while i < len(self.project.dicoms):
			di = self.project.dicoms[i]
			if self.insertDicom(END, di):
				i += 1
			else:
				del self.project.dicoms[i]

	# ----------------------------------------------------------------------
	def insertDicom(self, pos, di):
		uid = di["uid"]
		if uid is None: return False
		if di.files is None:
			try:
				Dicom.Scan.scan(self, di["directory"])
				di["frames"] = Dicom.Scan.frames(uid)
			except KeyError:
				return False
			di.files     = Dicom.Scan.files(uid)
			di.modality  = Dicom.Scan.modality(uid)
		filename = os.path.join(di["directory"], di.files[0][0])
		header = Dicom.dicomHeader(di["uid"], filename)
		self.dicomList.insert(pos, header)
		return True

#===============================================================================
# Dicom Select Dialog
#===============================================================================
class DicomSelectDialog(Toplevel):
	def __init__(self, master, project, path="", uid=None):
		Toplevel.__init__(self, master)
		self.title("Select DICOM data set")
		self.transient(master)

		# -------------- Directory ---------------
		Label(self, text="Directory:").grid(
				row=0, column=0, pady=1, sticky=E)
		self.path = Label(self,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				#background=tkFlair._ILABEL_BACKGROUND_COLOR,
				background="LightYellow",
				relief=SUNKEN, anchor=W)
		self.path.grid(row=0, column=1, sticky=EW)
		self.path.bind('<Button-1>', self.chdir)
		tkExtra.Balloon.set(self.path, "Dicom directory")

		b = Button(self, image=tkFlair.icons["load"], pady=0,
				command=self.chdir)
		b.grid(row=0, column=2)
		tkExtra.Balloon.set(b, "Change Dicom directory")

		# -------------- DataSets ---------------
		lblframe = LabelFrame(self, text="Data sets", foreground="DarkBlue")
		lblframe.grid(row=1, column=0, columnspan=3, sticky=NSEW)

		self.dicomList = tkExtra.MultiListbox(lblframe,
				(('Series Instance UID', 30, None),
				 ('Modality',             8, None),
				 ('Date',                12, None),
				 ('Patient',             20, None),
				 ('Age',                  4, None),
				 ('Sex',                  3, None),
				 ('Institution',         10, None),
				 ('Study',               15, None)),
				height=10)
		self.dicomList.pack(expand=TRUE, fill=BOTH)
		self.dicomList.bindList("<Double-1>", self.ok)
		self.dicomList.bindList("<Return>",   self.ok)

		frame = Frame(self)
		frame.grid(row=2, column=0, columnspan=3, sticky=NSEW)
		b = Button(frame, text="Cancel",
				image=tkFlair.icons["x"],
				compound=LEFT,
				command=self.close)
		b.pack(side=RIGHT)
		b = Button(frame, text="Ok",
				image=tkFlair.icons["ok"],
				compound=LEFT,
				command=self.ok)
		b.pack(side=RIGHT)

		self.grid_columnconfigure(1, weight=1)
		self.grid_rowconfigure(1, weight=1)

		self.bind("<Escape>", self.close)
		self.protocol("WM_DELETE_WINDOW", self.close)

		self.project   = project

		if os.path.abspath(path) == Dicom.Scan.directory:
			self.path["text"] = path
			self.fill(uid)

		elif path != "":
			self.path["text"] = path
			self.scan()

		elif Dicom.Scan.directory != "":
			self.path.config(text=self.project.relativePath(Dicom.Scan.directory))
			self.fill(uid)

		self._result = None

	# ----------------------------------------------------------------------
	def show(self):
		self.deiconify()
		self.grab_set()
		if self.dicomList.size()==0:
			self.chdir()
		self.wait_window()
		return self._result

	# ----------------------------------------------------------------------
	def ok(self, event=None):
		sel = self.dicomList.curselection()
		self._result = [self.path["text"]]
		for i in sel:
			self._result.append(self.dicomList.get(i)[0])
		self.close()

	# ----------------------------------------------------------------------
	def close(self, event=None):
		self.grab_release()
		self.destroy()

	# ----------------------------------------------------------------------
	# Change dicom directory
	# ----------------------------------------------------------------------
	def chdir(self, event=None):
		self.grab_release()
		d = bFileDialog.askdirectory(master=self,
			title="Choose DICOM directory",
			initialdir=self.path["text"],
			filetypes=[("Dicom files","*.dcm"),
				("All","*")])
		if len(d) > 0:
			self.path["text"] = self.project.relativePath(d)
			self.config(cursor="watch")
			self.scan()
			self.config(cursor="")
		self.grab_set()

	# ----------------------------------------------------------------------
	# Scan directory for unique dicom files
	# ----------------------------------------------------------------------
	def scan(self):
		Dicom.Scan.scan(self, self.path["text"])
		self.fill()

	# ----------------------------------------------------------------------
	# Fill listbox
	# ----------------------------------------------------------------------
	def fill(self, select=None):
		self.dicomList.delete(0,END)
		for uid in list(Dicom.Scan.datasets.keys()):
			self.dicomList.insert(END, Dicom.Scan.headerUID(uid))
			if uid == select:
				self.dicomList.selection_set(END)
				self.dicomList.see(END)
				self.dicomList.activate(END)
