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
# Date:	20-Apr-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import re
import sys
import select
import string
import histogram
import subprocess
import time
import os
import PlotEngine

_range_pat  = re.compile(r"set ([xy]2?)range \[(.*?):(.*?)]")
_range_pat2 = re.compile(r".*currently \[(.*?):(.*?)]")

STYLE = [	"lines",
		"points",
		"linespoints",
		"impulses",
		"dots",
		"steps",
		"fsteps",
		"histeps",
		"errorbars",
		"yerrorbars",
		"xerrorbars",
		"xyerrorbars",
		"errorlines",
		"yerrorlines",
		"xerrorlines",
		"xyerrorlines",
		"filledcurves",
		"boxes",
		"boxerrorbars",
		"boxxyerrorbars" ]

SPEC = [	2,	#"lines",	To be checked!!!
		2,	#"points",
		2,	#"linespoints",
		2,	#"impulses",
		2,	#"dots",
		2,	#"steps",
		2,	#"fsteps",
		2,	#"histeps",
		3,	#"errorbars",
		3,	#"yerrorbars",
		3,	#"xerrorbars",
		4,	#"xyerrorbars",
		3,	#"errorlines",
		3,	#"yerrorlines",
		3,	#"xerrorlines",
		4,	#"xyerrorlines",
		2,	#"filledcurves",
		2,	#"boxes",
		4,	#"boxerrorbars",
		4 ]	#"boxxyerrorbars"

SMOOTH = [	"",
		"unique",
		"frequency",
		"cumulative",
		"csplines",
		"acsplines",
		"bezier" ]

defaultTerminal = ""
program = "gnuplot"

#===============================================================================
# Gnuplot class
#===============================================================================
class Gnuplot(PlotEngine.PlotEngine):
	"""
	Interface class to call gnuplot
	"""
	def __init__(self, args=None):
		PlotEngine.PlotEngine.__init__(self)
		self.init(args)

	# ----------------------------------------------------------------------
	def init(self, args=None):
		if args:
			prg = "%s %s"%(program, args)
		else:
			prg = program
			if not prg: prg = "gnuplot"
		pipe = subprocess.Popen(prg,
				shell=True, bufsize=0,
				stdin=subprocess.PIPE, stdout=subprocess.PIPE,
				stderr=subprocess.STDOUT, close_fds=True)
		self.pipe_in  = pipe.stdin
		self.pipe_out = pipe.stdout

	# ----------------------------------------------------------------------
	# Close the pipe
	# ----------------------------------------------------------------------
	def close(self):
		self.pipe_in.close()
		self.pipe_out.close()

	# ----------------------------------------------------------------------
	def __del__(self):
		self.quit()
		self.close()

	# ----------------------------------------------------------------------
	# Read from the pipe
	# ----------------------------------------------------------------------
	def readline(self):
		return self.pipe_out.readline().decode()

	# ----------------------------------------------------------------------
	# Send a string to command, followed by newline
	# ----------------------------------------------------------------------
	def writeline(self, s):
		self.pipe_in.write((s+'\n').encode())
		self.pipe_in.flush()

	# ----------------------------------------------------------------------
	# Send a string to command, followed by newline.
	# ----------------------------------------------------------------------
	def write(self, s):
		self.logMessage(s)
		try:
			self.writeline(s)
		except UnicodeEncodeError:
			sys.stderr.write("ERROR: Unicode characters are not accepted by gnuplot. Replaced with ?\n")
			self.writeline(s.encode("ascii","replace"))
		except IOError:
			# Try to restart gnuplot
			self.init()
			# Re-execute commands...
			self.writeline(s)

	# ----------------------------------------------------------------------
	def __call__(self, s):
		self.write(s)

	# ----------------------------------------------------------------------
	# Flush the buffer to the piped command
	# ----------------------------------------------------------------------
	def flush(self):
		self.pipe_in.flush()

	# ----------------------------------------------------------------------
	def check(self, timeout=None):
		#ready = select.select([self.pipe_out], [], [], 0.0001)
		#if self.pipe_out in ready[0]: return True

		# give a little time for buffers to fill
		try:
			ready = select.select([self.pipe_out], [], [], timeout)
		except KeyboardInterrupt:
			return False
		return self.pipe_out in ready[0]

	# ----------------------------------------------------------------------
	def reset(self):
		self("reset")

	# ----------------------------------------------------------------------
	def quit(self):
		self("\nquit")

	# ----------------------------------------------------------------------
	def title(self, title, options=""):
		self("set title \'%s\' %s" % (self.escape(title), options))

	# ----------------------------------------------------------------------
	def label(self, axis, label, options=""):
		self("set %slabel \'%s\' %s" % (axis, self.escape(label), options))

	# ----------------------------------------------------------------------
	def grid(self, on=True):
		if on:
			self("set grid")
		else:
			self("unset grid")

	# ----------------------------------------------------------------------
	def tics(self, axis, options=""):
		self("set %stics %s" % (axis, options))

	# ----------------------------------------------------------------------
	def log(self, axis="y", log=True):
		if log:
			self("set logscale %s"%(axis))
		else:
			self("unset logscale %s"%(axis))

	# ----------------------------------------------------------------------
	def range(self, axis, amin="", amax=""):
		if amin=="" and amax=="": return
		self("set %srange [%s:%s]" % (axis, amin, amax))

	# ----------------------------------------------------------------------
	def plot(self, cmd, data=None):
		if data is None:
			self("plot %s" % (cmd))
		else:
			self("plot \"-\" %s"%(cmd))
			# accept a list-of-lists of data and plot it interactively
			if isinstance(data,list):
				for line in data:
					self(" ".join(map(str,line)))
			else:
				self(data)
			self("e")

	# ----------------------------------------------------------------------
	def splot(self, cmd, data=None):
		if data is None:
			self("splot %s" % (cmd))
		else:
			self("splot \"-\" %s"%(cmd))
			# accept a list-of-lists of data and plot it interactively
			if isinstance(data,list):
				for line in data:
					self(" ".join(map(str,line)))
			else:
				self(data)
			self("e")

	# ----------------------------------------------------------------------
	def replot(self):
		self("replot")

	# ----------------------------------------------------------------------
	def x11(self):
		self("set terminal x11 raise")

	# ----------------------------------------------------------------------
	def escape(self, str):
		"""Escape characters in string"""
		return str.replace("'","''")

	# ----------------------------------------------------------------------
	def output(self, output=None):
		if output is None:
			self("unset output")
		else:
			self("set output \"%s\"" % (output))

	# ----------------------------------------------------------------------
	def set(self, cmd=""):
		self("set %s"%(cmd))

	# ----------------------------------------------------------------------
	# Send end data marker
	# ----------------------------------------------------------------------
	def end(self):
		self("e")

	# ----------------------------------------------------------------------
	def empty(self):
		while self.check(1.0):
#			sys.stdout.write("\033[1m<e<%s\033[m\n"%(self.readline()[:-1]))
			print("\033[1m<e<%s\033[m\n"%(self.readline()[:-1]))

	# ----------------------------------------------------------------------
	def version(self):
		"""return version number of gnuplot"""
		if Gnuplot._version is None:
			self("show version")
			cnt = 0
			while self.check(1.0):
				words = self.readline().strip().split()
				if len(words)>0 and words[0] == "Version":
					Gnuplot._version = float(words[1])
				cnt += 1
				if cnt > 1000:
					Gnuplot._version = 0
					raise Exception("Gnuplot is not responding")
		return Gnuplot._version

	# ----------------------------------------------------------------------
	def getTerminal(self):
		"""return the default terminal of gnuplot"""
		self("show terminal")
		while self.check(1.0):
			words = self.readline().strip().split()
			if len(words)>0 and words[0] == "terminal":
				return words[3]
		return None

	# ----------------------------------------------------------------------
	def getRange(self, axis="x"):
		self("show %srange"%(axis))
		low, high  = None, None
		timeout = 1.0
		while self.check(timeout):
			line = self.readline().strip()
			m = _range_pat.match(line)
			if m:
				a = m.group(1)
				if a != axis: continue
				l = m.group(2)
				h = m.group(3)
				try:
					low  = float(l)
					high = float(h)
					if low > high: high,low = low,high
				except:
					pass

			if low is None or high is None:
				m = _range_pat2.match(line)
				if m:
					l = m.group(1)
					h = m.group(2)
					try:
						low  = float(l)
						high = float(h)
						if low > high: high,low = low,high
					except:
						pass
			timeout = 0.1
		return low, high

	# ----------------------------------------------------------------------
	def getRange2(self, axes):
		self("show xrange")
		self("show yrange")
		xlow, xhigh  = None, None
		ylow, yhigh  = None, None
		timeout = 1.0
		while self.check(timeout):
			line = self.readline()
			m = _range_pat.match(line.strip())
			if m:
				axis = m.group(1)
				low  = m.group(2)
				high = m.group(3)
				try: low = float(low)
				except: continue
				try: high = float(high)
				except: continue
				if low>high: high,low = low,high
				if axis=="x":
					xlow, xhigh = low, high
				else:
					ylow, yhigh = low, high
			timeout = 0.1
		return xlow, xhigh, ylow, yhigh

	# ----------------------------------------------------------------------
	# Change plot engine directory to path
	# ----------------------------------------------------------------------
	def cd(self, path):
		self("cd '%s'"%(path))

	# ----------------------------------------------------------------------
	def test(self):
		self("test")

#-------------------------------------------------------------------------------
if __name__ == "__main__":
	import time
	gp = Gnuplot()
	mx = 10.0
	#print("Version=",gp.version())
	while mx<10000.0:
		#print(mx)
		gp.reset()
		gp("plot [0:%d] sin(x)" % (mx))
		gp("show xrange")
		#gp("show version")
		while gp.check(0.05):
			print(gp.readline())
#			sys.stdout.write("%s\n"%(gp.readline()))
		#	sys.stdout.write("%...time.time(), gp.readline()
		j = 0
		#for i in range(0,10):
		#	j += 1
		#	sys.stdin.readline()
		mx += 1
