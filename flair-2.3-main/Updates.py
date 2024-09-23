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
# Date:	5-Apr-2007

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import re
import time

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkFlair
import tkExtra

flair_pat = re.compile(r".*^flair-(\d+\.\d+.\d+)\.tgz (.*?)$", re.MULTILINE | re.DOTALL)
fluka_pat = re.compile(r".*^The latest version is.*?: (\d+..*?)$.*?: (.*?)$", \
			re.MULTILINE | re.DOTALL)

#===============================================================================
# Check for updates in flair and fluka
#===============================================================================
class CheckUpdateDialog(Toplevel):
	def __init__(self, master, flair_version, fluka_version):
		Toplevel.__init__(self, master)
		self.transient(master)
		self.title("Check for updates")

		# Variables
		self.flair_version  = flair_version
		self.fluka_version  = fluka_version
		self.download_flair = False
		self.download_fluka = False

		# ----
		frame = LabelFrame(self, text="Flair", padx=3, pady=5)
		frame.pack(side=TOP, fill=BOTH)

		l = Label(frame, text="Installed Version:")
		l.grid(row=0, column=0, sticky=E, pady=1)

		l = Label(frame, text=flair_version, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		l.grid(row=0, column=1, sticky=EW)
		tkExtra.Balloon.set(l, "Running version of flair")

		l = Label(frame, text="Web Version:")
		l.grid(row=1, column=0, sticky=E, pady=1)

		self.flair_webversion = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.flair_webversion.grid(row=1, column=1, sticky=EW)
		tkExtra.Balloon.set(self.flair_webversion,
			"Version on the web site of flair")
		l = Label(frame, text="Released on:")
		l.grid(row=2, column=0, sticky=E, pady=1)

		self.flair_released = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.flair_released.grid(row=2, column=1, sticky=EW)
		tkExtra.Balloon.set(self.flair_released,
			"Release date of the flair web version")

		frame.grid_columnconfigure(1, weight=1)

		# ----
		frame = LabelFrame(self, text="Fluka", padx=3, pady=5)
		frame.pack(fill=BOTH)

		l = Label(frame, text="Installed Version:")
		l.grid(row=0, column=0, sticky=E, pady=1)

		l = Label(frame, text=fluka_version, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		l.grid(row=0, column=1, sticky=EW)
		tkExtra.Balloon.set(l, "FLUKA version on the $FLUPRO directory")

		l = Label(frame, text="Web Version:")
		l.grid(row=1, column=0, sticky=E, pady=1)

		self.fluka_webversion = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.fluka_webversion.grid(row=1, column=1, sticky=EW)
		tkExtra.Balloon.set(self.fluka_webversion,
			"Version on the web site of fluka")

		l = Label(frame, text="Released on:")
		l.grid(row=2, column=0, sticky=E, pady=1)

		self.fluka_released = Label(frame, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.fluka_released.grid(row=2, column=1, sticky=EW)
		frame.grid_columnconfigure(1, weight=1)
		tkExtra.Balloon.set(self.flair_released,
			"Release date of the FLUKA web version")

		# ----
		frame = LabelFrame(self, text="Check Interval", padx=3, pady=5)
		frame.pack(fill=BOTH)

		l = Label(frame, text="Last Check:")
		l.grid(row=0, column=0, sticky=E, pady=1)

		# Last check
		lastCheck = tkFlair.getInt(tkFlair._FLAIR_SECTION,"lastcheck",0)
		if lastCheck == 0:
			lastCheckStr = "unknown"
		else:
			lastCheckStr = time.asctime(time.localtime(lastCheck))

		l = Label(frame, text=lastCheckStr, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		l.grid(row=0, column=1, sticky=EW)
		tkExtra.Balloon.set(l, "Date of last checking")

		l = Label(frame, text="Interval (days):")
		l.grid(row=1, column=0, sticky=E, pady=1)

		checkInt = tkFlair.getInt(tkFlair._FLAIR_SECTION,"checkinterval",30)
		self.checkInterval = IntVar()
		self.checkInterval.set(checkInt)

		s = Spinbox(frame, text=self.checkInterval, from_=0, to_=365,
				background="White")
		s.grid(row=1, column=1, sticky=EW)
		frame.grid_columnconfigure(1, weight=1)
		tkExtra.Balloon.set(s, "Days-interval to remind again for checking")

		# ----
		frame = Frame(self)
		frame.pack(side=BOTTOM,fill=X)
		b = Button(frame,text="Close",
				image=tkFlair.icons["x"],
				compound=LEFT,
				command=self.later)
		b.pack(side=RIGHT)

		self.checkButton = Button(frame,text="Check Now",
				image=tkFlair.icons["GLOBAL"],
				compound=LEFT,
				command=self.check)
		self.checkButton.pack(side=RIGHT)
		tkExtra.Balloon.set(self.checkButton,
				"Check the web site for new versions of flair and FLUKA")

		self.bind('<Escape>', self.close)

		#x = master.winfo_rootx() + 200
		#y = master.winfo_rooty() + 50
		#self.geometry("+%d+%d" % (x,y))
		#self.wait_visibility()
		self.wait_window()

	# ----------------------------------------------------------------------
	def check(self):
		import socket

		# Find latest version of flair
		self.flair_webversion.config(text="...waiting...")
		self.fluka_webversion.config(text="...waiting...")
		self.flair_webversion.update_idletasks()
		self.fluka_webversion.update_idletasks()

		version  = None
		released = ""
		self.flair_webversion.config(
			text="Establishing connection\nwww.fluka.org/flair")
		self.flair_webversion.update_idletasks()
		try:
			s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
			s.settimeout(10.0)
			try:
				s.connect(("www.fluka.org",80))
				s.send("GET /flair/version.tag HTTP/1.0\n\n")
				receive = s.recv(1024)
			except:
				self.flair_webversion.config(text="Error in connection")
				self.fluka_webversion.config(text="... Interrupted ...")
				return
			m = flair_pat.match(receive)
			if m:
				version, released = m.groups()
			else:
				self.flair_webversion.config(text="Error in connection")
				self.fluka_webversion.config(text="... Interrupted ...")
				return
			s.close()
			self.flair_webversion.config(text=version)
			self.flair_released.config(text=released)
			if self.flair_version != version:
				self.flair_webversion.config(
					foreground=tkFlair._CLABEL_FOREGROUND_COLOR)
				self.download_flair = True
		except KeyboardInterrupt:
			self.flair_webversion.config(text="... Interrupted ...")
			self.fluka_webversion.config(text="... Interrupted ...")
			return

		version  = None
		released = ""
		self.fluka_webversion.config(text="Establishing connection\nwww.fluka.org/flair")
		self.fluka_webversion.update_idletasks()
		try:
			s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
			s.settimeout(10.0)
			try:
				s.connect(("www.fluka.org",80))
				s.send("GET /Version.tag HTTP/1.0\n\n")
				receive = s.recv(1024)
			except:
				self.fluka_webversion.config(text="Error in connection")
				return
			m = fluka_pat.match(receive)
			if m:
				version, released = m.groups()
			else:
				self.fluka_webversion.config(text="Error in connection")
				return
			s.close()
			self.fluka_webversion.config(text=version)
			self.fluka_released.config(text=released)
			if self.fluka_version != version:
				self.fluka_webversion.config(
					foreground=tkFlair._CLABEL_FOREGROUND_COLOR)
				self.download_fluka = True
		except KeyboardInterrupt:
			pass

		if self.download_flair or self.download_fluka:
			self.checkButton.config(text="Download", command=self.download)
			tkExtra.Balloon.set(self.checkButton,
					"Open web browser to download flair/FLUKA")
		else:
			self.checkButton.config(state=DISABLED)
		#self.laterButton.config(state=DISABLED)

		# Save today as lastcheck date
		tkFlair.config.set(tkFlair._FLAIR_SECTION,
			"lastcheck", str(int(time.time())))

	# ----------------------------------------------------------------------
	def later(self):
		# Save today as lastcheck date
		tkFlair.config.set(tkFlair._FLAIR_SECTION,
			"lastcheck", str(int(time.time())))
		self.close()

	# ----------------------------------------------------------------------
	def download(self):
		if self.download_flair:
			tkFlair.openurl("http://www.fluka.org/flair/download.html")
		if self.download_fluka:
			tkFlair.openurl("http://www.fluka.org/download.html")

	# ----------------------------------------------------------------------
	def close(self, event=None):
		try:
			tkFlair.config.set(tkFlair._FLAIR_SECTION,
				"checkinterval", str(int(self.checkInterval.get())))
		except TypeError:
			pass
		self.destroy()

#-------------------------------------------------------------------------------
# Find FLUKA version
#-------------------------------------------------------------------------------
def flukaVersion(dir):
	try:
		f = open(os.path.join(dir,"Version.tag"),"r")
		str = f.read()
		m = fluka_pat.match(str)
		if m:
			version, released = m.groups()
		f.close()
		return version
	except:
		return "Unknown"

#-------------------------------------------------------------------------------
# Check if interval has passed from last check
#-------------------------------------------------------------------------------
def need2Check():
	lastCheck = tkFlair.getInt(tkFlair._FLAIR_SECTION,"lastcheck",0)
	if lastCheck == 0:	# Unknown
		return True

	checkInt = tkFlair.getInt(tkFlair._FLAIR_SECTION,"checkinterval",30)
	if checkInt == 0:	# Check never
		return False

	return lastCheck + checkInt*86400 < int(time.time())

#===============================================================================
if __name__ == "__main__":
	tk = Tk()
	dlg = CheckUpdateDialog(tk,0,0)
	tk.mainloop()
