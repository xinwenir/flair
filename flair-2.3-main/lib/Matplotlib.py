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

try:
	import matplotlib
	import matplotlib.pyplot as pyplot
except ImportError:
	matplotlib = None

import PlotEngine

# Colors with different name in Matplotlib than gnuplot
_COLOR_TRANSLATE = {
	"dark-blue"         : "darkblue",
	"dark-cyan"         : "darkcyan",
	"dark-goldenrod"    : "darkgoldenrod",
	"dark-gray"         : "darkgray",
	"dark-green"        : "darkgreen",
	"dark-khaki"        : "darkkhaki",
	"dark-magenta"      : "darkmagenta",
	"dark-olivegreen"   : "darkolivegreen",
	"dark-orange"       : "darkorange",
	"dark-red"          : "darkred",
	"dark-salmon"       : "darksalmon",
	"dark-spring-green" : "darkseagreen",
	"dark-turquoise"    : "darkturquoise",
	"light-blue"        : "lightblue",
	"light-coral"       : "lightcoral",
	"light-cyan"        : "lightcyan",
	"light-goldenrod"   : "lightgoldenrodyellow",
	"light-gray"        : "lightgray",
	"light-green"       : "lightgreen",
	"light-pink"        : "lightpink",
	"light-salmon"      : "lightsalmon",
	"orange-red"        : "orangered",
}

#===============================================================================
# Matplotlib class
#===============================================================================
class Matplotlib(PlotEngine.PlotEngine):
	"""
	Interface class to call matplotlib
	"""
	def __init__(self):
		PlotEngine.PlotEngine.__init__(self)

	# ----------------------------------------------------------------------
	def version(self):
		if Matplotlib._version is None:
			a,b,c = matplotlib.__version__.split(".")
			Matplotlib._version = float(a) + float(b)/10.0
		return Matplotlib._version

	# ----------------------------------------------------------------------
	def title(self, title, *options):
		pyplot.title(title, *options)

	# ----------------------------------------------------------------------
	def label(self, axis, label, *options):
		if axis=="x":
			pyplot.xlabel(label, *options)
		elif axis=="y":
			pyplot.ylabel(label, *options)
		elif axis=="z":
			pyplot.zlabel(label, *options)

	# ----------------------------------------------------------------------
	def grid(self, on):
		pyplot.grid(on)

	# ----------------------------------------------------------------------
	def tics(self, axis, options=''):
		pass

	# ----------------------------------------------------------------------
	def log(self, axis="y", log=True):
		if axis=="x":
			pyplot.gca().set_xscale(log and 'log' or 'linear')	# symlog?
		elif axis=="y":
			pyplot.gca().set_yscale(log and 'log' or 'linear')	# symlog?
#		elif axis=="z":
#			pyplot.gca().set_zscale(log and 'log' or 'linear')	# symlog?

	# ----------------------------------------------------------------------
	def range(self, axis, amin, amax):
		axes = pyplot.gca()
		args = {}
		if axis=="x":
			if amin:
				try: axes.set_xlim(left=float(amin))
				except: pass
			if amax:
				try: axes.set_xlim(right=float(amax))
				except: pass
			#axes.set_xlim(**args)
			#axes.set_xlim((amin, amax))
		elif axis=="y":
			if amin:
				try: axes.set_ylim(bottom=float(amin))
				except: pass
			if amax:
				try: axes.set_ylim(top=float(amax))
				except: pass
			#axes.set_ylim((amin, amax))
#		elif axis=="z":
#			axes.set_zlim((amin, amax))
#		print "args=",args

	# ----------------------------------------------------------------------
	def show(self, block=False):
		pyplot.show(block)

	#-----------------------------------------------------------------------
	def color(self, name):
		if name in _COLOR_TRANSLATE:
			return _COLOR_TRANSLATE[name]
		return name

#-------------------------------------------------------------------------------
def hasMatPlotLib():
	return matplotlib is not None
