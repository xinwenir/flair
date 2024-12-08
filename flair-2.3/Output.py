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
# Date:	03-Oct-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import sys
import time
import Ribbon
import Project
import tkFlair
import tkExtra
import Process
import bFileDialog
import FlairRibbon

import Input
try:
	from tkinter import *
except ImportError:
	from tkinter import *

stdout_orig = sys.stdout
stderr_orig = sys.stderr

_ORDER = {"flair"   : 0,
	  "Compile" : 1,
	  "Run"     : 2,
	  "Data"    : 3,
	  "Plot"    : 4,
	  "Dicom"   : 5
	 }

#===============================================================================
# NullLog
# Empty wrapper to display messages on screen
#===============================================================================
class NullLog:
	def __init__(self, process):
		sys.stdout.write("%s%s %s%s\n"% \
			(tkExtra.ANSI_BOLD, time.ctime(), process,
			 tkExtra.ANSI_NORMAL))

	# ----------------------------------------------------------------------
	def write(self, text):
		sys.stdout.write(text)

	# ----------------------------------------------------------------------
	def writeln(self, text):
		sys.stdout.write("%s\n"%(text))

	# ----------------------------------------------------------------------
	def __call__(self, text):
		self.writeln(text)

#===============================================================================
# Log
#===============================================================================
class Log:
	def __init__(self, page, typ, name):
		self.page = page
		self.type = typ
		self.idx  = 0
		self.name = name
		self._status  = 0
		self._percent = 0
		self.reset()

	# ----------------------------------------------------------------------
	def reset(self):
		self.messages  = ""
		self._status   = 0

	# ----------------------------------------------------------------------
	def status(self, st, p):
		if self._status!=st or self._percent!=p:
			self._status  = st
			self._percent = p
			self.page.updateStatus(self)

	# ----------------------------------------------------------------------
	def statusText(self):
		if self._status == Process.STATUS_RUNNING:
			return "%s %d%%"%(Process.STATUS[self._status], self._percent)
		else:
			return Process.STATUS[self._status]

	# ----------------------------------------------------------------------
	def statusColor(self):
		return Process.COLOR[self._status]

	# ----------------------------------------------------------------------
	def timestamp(self, msg):
		self.messages  += "%s %s\n" % \
			(msg, time.strftime(bFileDialog._TIME_FORMAT, time.localtime(time.time())))

	# ----------------------------------------------------------------------
	def write(self, text):
		if Input._developer:
			sys.stdout.write(str(text))
		if not self.messages: self.timestamp("Started")
		self.messages += str(text)
		self.page.refreshLate()		# FIXME not very brilliant
						# it needs to insert just the last string
						# to insert(END, text\n)
						# if the proper item is selected
						# otherwise highlight the item
						# that new info is printed

	# ----------------------------------------------------------------------
	def writeln(self, text):
		self.write(str(text)+"\n")

	# ----------------------------------------------------------------------
	def __call__(self, text):
		self.writeln(text)

	# ----------------------------------------------------------------------
	def size(self):
		return len(self.messages)

#===============================================================================
# Output Text Stream
#===============================================================================
class OutputBaseStream:
	txt = None

	def __init__(self, master):
		if OutputBaseStream.txt is None:
			OutputBaseStream.txt = Text(master)

	# ----------------------------------------------------------------------
	def write(self, text):
		pass

#===============================================================================
class OutputStream(OutputBaseStream):
	def __init__(self, master):
		OutputBaseStream.__init__(self, master)

	# ----------------------------------------------------------------------
	def write(self, text):
		OutputBaseStream.txt.insert(END, text)

#===============================================================================
class ErrorStream(OutputBaseStream):
	def __init__(self, master):
		OutputBaseStream.__init__(self, master)

	# ----------------------------------------------------------------------
	def write(self, text):
		OutputBaseStream.txt.insert(END, text)

#===============================================================================
# Output frame class
#===============================================================================
class OutputPage(FlairRibbon.FlairPage):
	"""Monitor the progress of a process and its output"""

	_name_ = "Output"
	_icon_ = "output"

	# ----------------------------------------------------------------------
	def init(self):
		self._log = []	# list of log structures
		self._refreshAfter = None

	# ----------------------------------------------------------------------
	# create the ribbon buttons
	# ----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		group = Ribbon.LabelGroup(self.ribbon, "Process")
		group.pack(side=LEFT, fill=BOTH, expand=TRUE, padx=0, pady=0)

		# ---
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["play32"],
				text="Run",
				compound=TOP,
				state=DISABLED,
#				command=self.execute,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y, padx=0, pady=0)
		tkExtra.Balloon.set(b, "Execute selected processes [Ctrl-Enter]")

		# ---
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["stop32"],
				text="Stop",
				compound=TOP,
				state=DISABLED,
#				command=self.stopRun,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y, padx=0, pady=0)
		tkExtra.Balloon.set(b, "Stop selected processes")

		# ---
#		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x32"],
#				text="Files",
#				compound=TOP,
#				state=DISABLED,
##				command=self.cleanFiles,
#				background=Ribbon._BACKGROUND)
#		b.pack(side=LEFT, fill=Y, padx=0, pady=0)
#		tkExtra.Balloon.set(b, "Delete all generated files from this process.")

		# ---
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clean32"],
				text="Clean",
				compound=TOP,
				command=self.deleteOutput,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y, padx=0, pady=0)
		tkExtra.Balloon.set(b, "Clean output of selected processes")

		# ---
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x32"],
				text="Close",
				compound=TOP,
				command=self.close,
				background=Ribbon._BACKGROUND)
		b.pack(side=LEFT, fill=Y, padx=0, pady=0)
		tkExtra.Balloon.set(b, "Close output")
		return self.ribbon

	#----------------------------------------------------------------------
	# Create Project page
	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self)

		# --- Run-Usrxxx splitter ---
		splitter = PanedWindow(self.frame, orient=VERTICAL,
			showhandle=0,
			sashrelief=RAISED, sashwidth=4, opaqueresize=1)
		splitter.pack(expand=YES, fill=BOTH)

		# List of commands
		self.listbox = tkExtra.MultiListbox(splitter,
				(("Type",     8, None),
				 ("Process", 20, None),
				 ("Status",  10, None)),
				 height=6)
		self.listbox.sortAssist = None
		splitter.add(self.listbox)

		# Text with output messages
		frame = Frame(splitter)
		splitter.add(frame)

		self.msg = Text(frame, state=DISABLED, takefocus=TRUE, background="White", height=6)
		self.msg.grid(row=0, column=0, sticky=NSEW)
		sbv = Scrollbar(frame, orient=VERTICAL, takefocus=False,
			command=self.msg.yview)
		sbv.grid(row=0, column=1, sticky=NS)
		sbh = Scrollbar(frame, orient=HORIZONTAL, takefocus=False,
			command=self.msg.xview)
		sbh.grid(row=1, column=0, sticky=EW)
		self.msg.config(yscrollcommand=sbv.set,
				xscrollcommand=sbh.set)

		self.msg.bind("<Button-1>", self._focusIn)

		frame.grid_columnconfigure(0, weight=1)
		frame.grid_rowconfigure(0, weight=1)

		self.populateLog()

		self.listbox.bind("<<ListboxSelect>>", self.refresh)

	# ----------------------------------------------------------------------
	def populateLog(self):
		self.listbox.delete(0,END)
		for log in self._log:
			self.listbox.insert(END, (log.type, log.name, log.statusText()))
		self.listbox.selection_set(END)
		self.listbox.see(END)

	# ----------------------------------------------------------------------
	# Return log with name
	# ----------------------------------------------------------------------
	def newLog(self, typ, name, append=False):
		self.listbox.selection_clear(0,END)
		for i, log in enumerate(self._log):
			if log.name == name and log.type==typ:
				if not append: log.reset()
				self.listbox.selection_set(i)
				self.refresh()
				return log
		else:
			log = Log(self,typ,name)
			self.insert(log)
			self.refresh()
			return log

	# ----------------------------------------------------------------------
	# Find position to insert log
	# ----------------------------------------------------------------------
	def insert(self, log):
		order = _ORDER.get(log.type, len(_ORDER))

		# Find the last
		for i,l in enumerate(self._log):
			if _ORDER.get(l.type, len(_ORDER)) > order:
				self._log.insert(i,log)
				self.listbox.insert(i, (log.type, log.name, log.statusText()))
				self.listbox.selection_set(i)
				self.listbox.see(i)
				break
		else:
			self._log.append(log)
			self.listbox.insert(END, (log.type, log.name, log.statusText()))
			self.listbox.selection_set(END)
			self.listbox.see(END)

		for i,l in enumerate(self._log):
			l.idx = i

	# ----------------------------------------------------------------------
	def refreshLate(self):
		# FIXME Normally it should be protected with semaphors!!!!
		if self._refreshAfter is not None:
			self.page.after_cancel(self._refreshAfter)
		if self.page:
			try:
				self._refreshAfter = self.page.after(1000, self.refresh)
			except TclError:
				# catch possible infinite loop error when multiple writes happen
				pass

	# ----------------------------------------------------------------------
	# Update the display
	# ----------------------------------------------------------------------
	def refresh(self, event=None):
		self._refreshAfter = None
		if self.page is None: return
		FlairRibbon.FlairPage.refresh(self)
		try:
			sel = self.listbox.curselection()
		except TclError:
			self.page = None
			return

		self.msg.config(state=NORMAL)
		self.msg.delete("0.0",END)
		if len(sel)==1:
			idx = int(sel[0])
			try:
				log = self._log[idx]
				self.msg.insert("0.0", log.messages)
				self.msg.see(END)
			except IndexError:
				pass

		self.msg.config(state=DISABLED)

	# ----------------------------------------------------------------------
	def deleteOutput(self, event=None):
		for sel in self.listbox.curselection():
			idx = int(sel)
			try:
				log = self._log[idx]
				log.reset()
			except IndexError:
				pass
		self.refresh()

	# ----------------------------------------------------------------------
	def updateStatus(self, log):
		sel = self.listbox.selection_includes(log.idx)
		self.listbox.set(log.idx, (log.type, log.name, log.statusText()))
		if log._status != 0:
			self.listbox.itemconfig(log.idx, foreground=log.statusColor())

		if sel: self.listbox.selection_set(log.idx)

	# ----------------------------------------------------------------------
	def capture(self):
		sys.stdout = self.stdout
		sys.stderr = self.stderr

	# ----------------------------------------------------------------------
	def release(self):
		sys.stdout = stdout_orig
		sys.stderr = stderr_orig

	# ----------------------------------------------------------------------
	def copy(self, event=None):
		if event.widget is self.msg: return
		self.msg.event_generate("<<Copy>>")
	cut = copy

	# ----------------------------------------------------------------------
	def close(self, event=None):
		self.hide(True)
