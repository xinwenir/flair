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
# Date:	12-Dec-2014

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import ntuple
import histogram

_COLOR_NAMES = ["",
		"antiquewhite",
		"aquamarine",
		"beige",
		"bisque",
		"black",
		"blue",
		"brown",
		"chartreuse",
		"coral",
		"cyan",
		"dark-blue",
		"dark-cyan",
		"dark-goldenrod",
		"dark-gray",
		"dark-green",
		"dark-khaki",
		"dark-magenta",
		"dark-olivegreen",
		"dark-orange",
		"dark-red",
		"dark-salmon",
		"dark-spring-green",
		"dark-turquoise",
		"gold",
		"goldenrod",
		"gray",
		"green",
		"greenyellow",
		"honeydew",
		"khaki",
		"lemonchiffon",
		"light-blue",
		"light-coral",
		"light-cyan",
		"light-goldenrod",
		"light-gray",
		"light-green",
		"light-pink",
		"light-salmon",
		"magenta",
		"navy",
		"olive",
		"orange",
		"orange-red",
		"orchid",
		"pink",
		"plum",
		"purple",
		"red",
		"royalblue",
		"salmon",
		"sandybrown",
		"seagreen",
		"skyblue",
		"slategray",
		"steelblue",
		"turquoise",
		"violet",
		"white",
		"yellow"
]

POINTS = [	"dot",
		"+",
		"x",
		"*",
		"square",
		"square-filled",
		"circle",
		"circle-filled",
		"triangle-up",
		"triangle-up-filled",
		"triangle-down",
		"triangle-down-filled",
		"diamond",
		"diamond-filled"
	]

#===============================================================================
# Abstract Plot Engine class
#===============================================================================
class PlotEngine:
	_version = None

	"""
	Interface class to any plot egine
	"""
	def __init__(self, logger=None):
		self._logger = logger

	# ----------------------------------------------------------------------
	def setLogger(self, logger=None):
		self._logger = logger

	# ----------------------------------------------------------------------
	def logMessage(self, msg):
		if self._logger: self._logger(msg)

	# ----------------------------------------------------------------------
	def close(self):
		pass

	# ----------------------------------------------------------------------
	def reset(self):
		pass

	# ----------------------------------------------------------------------
	def quit(self):
		pass

	# ----------------------------------------------------------------------
	def title(self, title, opts=''):
		pass

	# ----------------------------------------------------------------------
	def label(self, axis, label, opts=''):
		pass

	# ----------------------------------------------------------------------
	def xlabel(self, label, opts=""):  self.label("x", label,opts)
	def ylabel(self, label, opts=""):  self.label("y", label,opts)
	def zlabel(self, label, opts=""):  self.label("z", label,opts)
	def cblabel(self, label, opts=""): self.label("cb",label,opts)

	# ----------------------------------------------------------------------
	def grid(self, on=True):
		pass

	# ----------------------------------------------------------------------
	def tics(self, axis, opts=''):
		pass

	# ----------------------------------------------------------------------
	def range(self, axis, amin="", amax=""):
		pass
	def xrange(self, amin="", amax=""):  self.range("x",  amin, amax)
	def yrange(self, amin="", amax=""):  self.range("y",  amin, amax)
	def zrange(self, amin="", amax=""):  self.range("z",  amin, amax)
	def cbrange(self, amin="", amax=""): self.range("cb", amin, amax)

	# ----------------------------------------------------------------------
	def plot(self, cmd, data):
		pass

	# ----------------------------------------------------------------------
	def splot(self, cmd, data):
		pass

	# ----------------------------------------------------------------------
	def draw(self, data, options=None):
		"""Generic drawing of data with options"""
		if isinstance(data,histogram.H1) or \
		   isinstance(data,histogram.H2):
			data.draw(self, options)
		else:
			raise Exception("Don't know how to handle datatype %s"%(str(type(data))))

	# ----------------------------------------------------------------------
	def replot(self):
		pass

	# ----------------------------------------------------------------------
	def escape(self, s):
		"""Escape characters in string"""
		return s.replace("'","''")

	# ----------------------------------------------------------------------
	def output(self, output=None):
		pass

	# ----------------------------------------------------------------------
	def set(self, cmd=""):
		pass

	# ----------------------------------------------------------------------
	def title(self, title):
		pass

	# ----------------------------------------------------------------------
	def version(self):
		"""return version number of plotting engine"""
		return PlotEngine._version

	# ----------------------------------------------------------------------
	def log(self, axis="y", log=True):
		pass

	# ----------------------------------------------------------------------
	def xlog(self):  self.log("x")
	def ylog(self):  self.log("y")
	def zlog(self):  self.log("z")
	def cblog(self): self.log("cb")

	# ----------------------------------------------------------------------
	def linear(self, axis):
		self.log(axis,False)

	def xlin(self):  self.linear("x")
	def ylin(self):  self.linear("y")
	def zlin(self):  self.linear("z")
	def cblin(self): self.linear("cb")

	# ----------------------------------------------------------------------
	def getRange(self, axis="x"):
		return 0.,1.

	# ----------------------------------------------------------------------
	def getRange2(self, axes):
		return 0., 1., 0., 1.

	# ----------------------------------------------------------------------
	# Change plot engine directory to path
	# ----------------------------------------------------------------------
	def cd(self, path):
		pass

	# ----------------------------------------------------------------------
	# Convert color name to something acceptable by the engine
	# To be overridden if needed
	# ----------------------------------------------------------------------
	def color(self, name):
		return name
