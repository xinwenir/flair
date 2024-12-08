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
# Date:	05-Sep-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import time
import math
#import string
#import random

try:
	import tkinter.messagebox as messagebox
except ImportError:
	import tkinter.messagebox as messagebox

import Data
import bmath
import tkFlair
#import tkExtra

import Input
import Project
try:
	import Matplotlib
except ImportError:
	Matplotlib = None
try:
	import numpy
except ImportError:
	numpy = None

from Plot import *

try:
	import matplotlib.pyplot as pyplot
except ImportError:
	pyplot = None

#===============================================================================
# ----------------------------------------------------------------------
# Base class to start plotting (Launches gnuplot)
# ----------------------------------------------------------------------
def start():
	if Matplotlib is None: return None
	engine = Matplotlib.Matplotlib()
	while 20:	# wait a while to get version
		try:
			ver = engine.version()
		except:
			return None
		if ver is not None: break

	# Initialize mouse and functions
	#engine("bin(x,m,s)=floor((x-m)/s)*s+m")
	#engine("logbin(x,m,s)=10**(floor((log10(x)-m)/s)*s+m)")
	#engine("rnd(x)=abs(int((x*1103515245+12345)/65536)%32768)")
	return engine

#===============================================================================
# MatPlotlib Plot class
#===============================================================================
class MPPlot(Plot):
	# ----------------------------------------------------------------------
	def __call__(self, s):
		self.write(s)
		self.plt.write(s)

	# ----------------------------------------------------------------------
	# Save plot to file
	# ----------------------------------------------------------------------
	def save(self, ext):
#		filename = Plot.save(self, ext)
#		if ext in (".gnu", ".gp"):
#			self("save \"%s\""%(filename))
#		else:
#			opt = tkFlair._PLOT_FORMAT.get(ext,None)
#			if opt:
#				self("set terminal %s %s"%(opt, self.plot["linetype"]))
#			else:
#				self("set terminal %s"%(ext[1:]))
#			self.engine.output(filename)
#			self.engine.replot()
#			self("set terminal pop")
#			self.engine.output()
#		self.checkOutput()
		return filename

	# ----------------------------------------------------------------------
	# Send plot on printer
	# ----------------------------------------------------------------------
	def hardcopy(self, printer):
#		Plot.hardcopy(self, printer)
#		if printer.landscape:
#			orient = "landscape"
#		else:
#			orient = "portrait"
#
#		self("set terminal postscript %s %s enhanced color" \
#			% (orient, self.plot["linetype"]))
#
#		if printer.printTo:
#			self.engine.output("|"+printer.command())
#		else:
#			self.engine.output(printer.filename)
#		self.engine.replot()
#		self("set terminal pop")
#		self.engine.output()
		pass

	# ----------------------------------------------------------------------
	def show(self):
		Plot.show(self)
		self.engine.show(False)

	# ----------------------------------------------------------------------
	# @return plot options for prefix: title, x, y, ...
	# ----------------------------------------------------------------------
	def _plotOptions(self, prefix):
		# Find the corresponding one in Matplotlib
		font = self.plot["%sfont" % (prefix)]
		size = self.plot["%ssize" % (prefix)]
		color = self.engine.color(self.plot["%scolor" % (prefix)])
		font = {
			# "fontname" : font,
			# "fontsize" : size,
			# "color"    : color
		}
		return font

	# ----------------------------------------------------------------------
	# Set color band
	# ----------------------------------------------------------------------
	def plotColorBand(self, pm3dopts=None):
		# Color band
		try: cbmin = str(float(self.plot["cbmin"]))
		except: cbmin = ""

		try: cbmax = str(float(self.plot["cbmax"]))
		except: cbmax = ""

		try: cpd = float(self.plot("cbcpd",3.0))
		except: cpd = 3.0

		cblog = int(self.plot("cblog",1))

		# If logscale... (not for balance)
		if cbmin!="" or cbmax!="":
			if cblog and int(self.plot("cbround",1)) and cbmin!="" and cpd>=1.0:
				# round limits that colors boundaries falls on decades
				try:
					nmin = math.floor(math.log10(float(cbmin))*cpd)/cpd
					cbmin = "%g"%(math.pow(10.0, nmin))
				except:
					cbmin = "*"
			self("set cbrange [%s:%s]" % (cbmin,cbmax))

		# Set the pm3d
		try:
			colors = int(self.plot("cbcolors", 30))
		except ValueError:
			colors = 30
		self.plotPm3D(colors, self.plot("cbpalette", "FLUKA"), pm3dopts)

		if cblog:
			self("set logscale cb")
			if 4.4 < self.engine.version() < 5.0:
				self("set logscale z")	# Necessary for 4.4.4+
								# to treat log(0) as 0 and not as undefined

	# ----------------------------------------------------------------------
	# Get axes range from engine
	# ----------------------------------------------------------------------
#	def getAxesRangeFromGnuplot(self):
#		if not self.engine: return
#		xlow, xhigh = self.engine.getRange("x")
#		ylow, yhigh = self.engine.getRange("y")
#		try:
#			zlow, zhigh = self.engine.getRange("z")
#		except:
#			pass
#
#		self.resetAxesRange(False)
#		if xlow is not None:  self.x_min.insert(0, xlow)
#		if ylow is not None:  self.y_min.insert(0, ylow)
#		if xhigh is not None: self.x_max.insert(0, xhigh)
#		if yhigh is not None: self.y_max.insert(0, yhigh)
#
#		try:
#			if zlow is not None:  self.z_min.insert(0, zlow)
#			if zhigh is not None: self.z_max.insert(0, zhigh)
#		except:
#			pass
#
#		try:
#			if self.x2_tics.get():
#				self.x2_min.delete(0,END)
#				self.x2_max.delete(0,END)
#				xlow, xhigh = self.engine.getRange("x2")
#
#				if xlow is not None:  self.x2_min.insert(0, xlow)
#				if xhigh is not None: self.x2_max.insert(0, xhigh)
#
#			if self.y2_tics.get():
#				self.y2_min.delete(0,END)
#				self.y2_max.delete(0,END)
#				ylow, yhigh = self.engine.getRange("y2")
#
#				if ylow is not None:  self.y2_min.insert(0, ylow)
#				if yhigh is not None: self.y2_max.insert(0, yhigh)
#		except:
#			pass

	# ----------------------------------------------------------------------
	# Plot Using
	# return using string for line plots depending on parameters requested
	# ----------------------------------------------------------------------
	def _plotUsing(self, style, x, y, ex, ey):
		try:
			specs = Gnuplot.SPEC[Gnuplot.STYLE.index(style)]
		except ValueError:
			specs = 2
		if   specs == 2:
			return "%s:%s" % (x, y)
		elif specs == 3:
			if style in ["xerrorbars", "xerrorlines"]:
				return "%s:%s:%s" % (x, y, ex)
			else:
				return "%s:%s:%s" % (x, y, ey)
		elif specs == 4:
			return "%s:%s:%s:%s" % (x, y, ex, ey)

	# ----------------------------------------------------------------------
	# Plot Style
	# @return string with plot options
	# ----------------------------------------------------------------------
	def _plotStyle(self, det, style=None):
		# axes
		axes = self.plot("axes.%d"%(det),"x1y1")
		if axes != "x1y1":
			cmd = "axes %s " % (axes)
		else:
			cmd = ""

		# with
		if style is None:
			style = self.plot("with.%d"%(det),STYLE_DEF)
			if style in ("histogram","histerror"): style = "steps"
		cmd += "w %s " % (style)

		# smooth
		smooth = self.plot("smooth.%d"%(det),"")
		if smooth != "":
			cmd += "smooth %s "%(smooth)

		# line style
		#if style not in ["points", "dots"]:
		x = int(self.plot("ls.%d"%(det),0))
		if x:
			cmd += "ls %s "%(x)
		else:
			x = self.plot["lt.%d"%(det)]
			if x: cmd += "lt %s "%(x)

			x = self.plot["lw.%d"%(det)]
			if x: cmd += "lw %s "%(x)


			if self.engine.version() >= 4.2:
				x = self.plot["lc.%d"%(det)]
				if x: cmd += "lc rgb '%s' "%(x)

			# point style
			if style not in ["lines", "steps", "impulses", "fsteps", "dots",\
					"filledcurves", "histeps", "boxes", \
					"boxerrorbars", "boxxyerrorbars"]:
				cmd += "pt %s ps %s " \
					% (self.plot("pt.%d"%(det),"1"),
					   self.plot("ps.%d"%(det),"1"))

		#if style in ["filledcurves", "boxes", \
		#		"boxerrorbars", "boxxyerrorbars"]:
		#	cmd += "fs %s " % (self.plot("fs.%d"%(det),"1"))
		return cmd

	# ----------------------------------------------------------------------
	# Perform plot
	# ----------------------------------------------------------------------
	def _plot(self, cmd):
		if self.engine.version() == 0: return None
		self.plotCommands()
		if self.plotappend != "":
			if cmd:
				cmd += ",%s"%(self.plotappend)
			else:
				cmd = self.plotappend
		if cmd: self("plot %s"%(cmd))
		return self._plotEnd()

	# ----------------------------------------------------------------------
	# Perform splot
	# ----------------------------------------------------------------------
	def _splot(self, cmd):
		if self.engine.version() == 0: return None
		self.plotCommands()
		if self.plotappend != "":
			if cmd:
				cmd += ",%s"%(self.plotappend)
			else:
				cmd = self.plotappend
		if cmd: self("splot %s"%(cmd))
		return self._plotEnd()

	# ----------------------------------------------------------------------
	# Base class to start plotting (Launches engine)
	# ----------------------------------------------------------------------
	def _plotEnd(self):
		self.plotCommands(False)
		return self.checkOutput()

	# ----------------------------------------------------------------------
	def checkOutput(self):
		self.engine.writeline('print "end"')
		# wait until we get the line end
		line  = "\n"
		END   = "end\n"
		err   = None
		start = time.time()
		while line!=END and time.time()-start < 10.0:	# 10s
			while self.engine.check(1.0):
				line = self.engine.readline()
				if line!="\n" and line!=END:
					self.log("<e< %s"%(line))
					if err is None:
						err = "Gnuplot error"
		return err

	# ----------------------------------------------------------------------
	# Enable pm3d mode for gnuplot
	# ----------------------------------------------------------------------
	def plotPm3D(self, colors=None, palette="FLUKA", pm3dopts=None):
		self("set colorbox vertical")
		if pm3dopts is None:
			self("set pm3d map explicit corners2color c1")
		else:
			self("set pm3d %s"%(pm3dopts))

		if not isinstance(palette, str):
			palette.gnuplot_setPalette(self.engine)
			return
		elif palette == PALETTE[0]:
			self('set palette defined (' \
				' 1. 1.0 1.0 1.0,  2. 0.9 0.6 0.9,  3. 1.0 0.4 1.0,' \
				' 4. 0.9 0.0 1.0,  5. 0.7 0.0 1.0,  6. 0.5 0.0 0.8,' \
				' 7. 0.0 0.0 0.8,  8. 0.0 0.0 1.0,  9. 0.0 0.6 1.0,' \
				'10. 0.0 0.8 1.0, 11. 0.0 0.7 0.5, 12. 0.0 0.9 0.2,' \
				'13. 0.5 1.0 0.0, 14. 0.8 1.0 0.0, 15. 1.0 1.0 0.0,' \
				'16. 1.0 0.8 0.0, 17. 1.0 0.5 0.0, 18. 1.0 0.0 0.0,' \
				'19. 0.8 0.0 0.0, 20. 0.6 0.0 0.0, 21. 0.0 0.0 0.0 )' )
		elif palette == PALETTE[1]:
			self('set palette defined ( 1 "grey",' \
				'2 "magenta", 3 "blue",4 "green",' \
				'5 "yellow", 6 "red", 7 "black")')
		elif palette == PALETTE[2]:
			self('set palette')
		elif palette == "rainbow":
			self('set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )')
		elif palette == "rgb":
			self('set palette defined ( 0 0 0 0, 0.25 0 0 1, 0.50 0 1 0, 0.75 1 0 0, 1 1 1 1 )')
		elif palette == "red-yellow":
			self("set palette model XYZ functions gray**0.35, gray**0.5, gray**0.8")
		elif palette == "geometry":
			self('set palette defined ( 1 "grey",' \
				'2 "magenta", 3 "blue",4 "green",' \
				'5 "yellow", 6 "red", 7 "white")')
		elif palette == "field":
			# Fluka colors where COLOR = (1+COLOR)/2
			self('set palette defined (' \
				' 1.  1.00  1.00  1.00,  2.  0.90  0.90  0.90,' \
				' 3.  0.95  0.80  0.95,  4.  1.00  0.70  1.00,' \
				' 5.  0.95  0.50  1.00,  6.  0.85  0.50  1.00,' \
				' 7.  0.75  0.50  0.90,  8.  0.50  0.50  0.90,' \
				' 9.  0.50  0.50  1.00, 10.  0.50  0.80  1.00,' \
				'11.  0.50  0.90  1.00, 12.  0.50  0.85  0.75,' \
				'13.  0.50  0.95  0.60, 14.  0.75  1.00  0.50,' \
				'15.  0.90  1.00  0.50, 16.  1.00  1.00  0.50,' \
				'17.  1.00  0.90  0.50, 18.  1.00  0.75  0.50,' \
				'19.  1.00  0.50  0.50, 20.  0.90  0.50  0.50,' \
				'21.  0.80  0.50  0.50, 22.  0.70  0.50  0.50,' \
				'23.  0.50  0.50  0.50 )' )
		elif palette == "Matlab":
			self('set palette defined (' \
				'0  0.0 0.0 0.5,' \
				'1  0.0 0.0 1.0,' \
				'2  0.0 0.5 1.0,' \
				'3  0.0 1.0 1.0,' \
				'4  0.5 1.0 0.5,' \
				'5  1.0 1.0 0.0,' \
				'6  1.0 0.5 0.0,' \
				'7  1.0 0.0 0.0,' \
				'8  0.5 0.0 0.0 )')
		else:
			self('set palette %s'%(palette))

		if colors:
			# Stupid bug fix of gnuplot
			# doesn't work when printing!!!!
			if self.engine.version() < 4.2: colors += 1
			self('set palette maxcolors %d' % (colors))

	# ----------------------------------------------------------------------
	# Set normalization factor
	# ----------------------------------------------------------------------
	def plotNorm(self, norm, column, axis="y"):
		# Normalization and default using
		if norm=="":
			return str(column)
		try:
			int(column)
			column="($%s)"%(column)
		except:
			pass

		try:
			norm = float(eval(norm, Input._globalDict))
			if isinstance(column,int):
				return "($%d*%g)" % (column, norm)
			else:
				return "(%s*%g)" % (column, norm)
		except:
			func=axis+"norm"
			if norm.find("$")>=0:
				return "(%s)"%(norm)
			else:
				self("%s(x)=%s"%(func,norm))
				if isinstance(column,int):
					return "(%s($%d))"%(func,column)
				else:
					return "(%s(%s))"%(func,column)

	# ----------------------------------------------------------------------
	# Plot commands
	# ----------------------------------------------------------------------
	def plotCommands(self, prefix=True):
		first   = True
		pline   = ""
		self.plotappend = ""
		for line in self.plot["commands"].splitlines():
			line = line.strip()
			if pline != "":
				line = pline + line
				pline = ""
			if line!="" and line[-1] == "\\":
				pline = line[:-1]
				continue

			words = line.split()
			if len(words)==0: continue
			# stop on the first plot command
			if words[0] in ["plot", "replot", "splot"]:
				if words[0] == "plot":
					self.plotappend = line[5:]
				if prefix: break
				if words[0] != "plot":
					self(line+"\n")
				first = False
			elif prefix == first:
				self(line+"\n")

	# ----------------------------------------------------------------------
	# Generate the plot
	# ----------------------------------------------------------------------
	def _geometryPlot(self):
		return
		err = Plot._geometryPlot(self)
		if err is not None: return err

		# Important variables
		st   = self.plot("subtype", SUBTYPES[0])
		cd   = self.plot("coord", COORD[0])
		O    = bmath.Vector(list(map(float,self.plot["origin"].split())))
		U    = bmath.Vector(list(map(float,self.plot["basisu"].split())))
		V    = bmath.Vector(list(map(float,self.plot["basisv"].split())))
		extU,extV = list(map(float, self.plot["extends"].split()))
		loP  = -extU*U -extV*V + O
		hiP  =  extU*U +extV*V + O
		Us = Vs = ""	# Sign of field operation
		if cd == COORD[0]:
			center = "0,0"
			using  = "us 1:2"
			self("set xrange [%g:%g]"%(-extU, extU))
			self("set yrange [%g:%g]"%(-extV, extV))
		elif cd == COORD[1]:
			center   = "0,0"
			using  = "us 2:1"
			self("set xrange [%g:%g]"%(-extV, extV))
			self("set yrange [%g:%g]"%(-extU, extU))
		else:
			u = ord(cd[ 0])-ord('X')
			v = ord(cd[-1])-ord('X')
			center = "%g,%g"%(O[u],O[v])
			self("set xrange [%g:%g]"%(loP[u],hiP[u]))
			self("set yrange [%g:%g]"%(loP[v],hiP[v]))
			u += 3
			v += 3
			using  = "us %d:%d" % (u,v)
			if U.direction()[0]=="-": Us = "-"
			if V.direction()[0]=="-": Vs = "-"

		# Styles
		self("set style line 1 lt -1 lw 1")
		self("set style line 2 lt 6 lw 1")

		# Vector field scaling
		try:
			self("vscale=%g"%(float(self.plot("vectscale", 0.1))))
		except ValueError:
			self("vscale=0.1")

		# Do the plot
		if st == SUBTYPES[0]:
			return self._plot("'%s.dat' ind 0 %s w l ls 1 notitle" \
					% (self.plot.name, using))

		elif st == SUBTYPES[-2]:		# Field Vector
			every = 5
			return self._plot("'%s.dat' ind 0 %s w l ls 1 notitle," \
				"'' ind 1 every %d:%d %s:(%s$9*vscale):(%s$10*vscale) w vectors ls 2 notitle" \
				% (self.plot.name, using, every, every, using, Us, Vs))

		else:
			extra = ""
			colors = 32
			if st == SUBTYPES[2]:		# Materials
				try:
					colors = self.project.gnuplotPalette(self.plot._materialList)
					self.plotPm3D(len(colors),
						"defined (%s)" % " ".join(colors, ","))
					# FIXME correct bug that might set cblog by mistake to 1!
					if int(self.plot("cblog",0)):
						self("unset logscale cb")
					self("set cbrange [1:%d]"%(len(colors)))
				except AttributeError:
					# Could be that due to some error the _materialList is not defined
					pass

			elif st == SUBTYPES[-1]: # Field Intensity + Vectors
				every = 5
				extra=",'' ind 1 every %d:%d %s:(0):(%s$9*vscale):(%s$10*vscale):(0) w vectors ls 2 notitle" \
					%(every, every, using, Us, Vs)
				self.plotPm3D(colors, "field")

			elif st in (SUBTYPES[-3], SUBTYPES[-1]):
				self.plotPm3D(colors, "FLUKA")
			else:
				self.plotPm3D(colors, "geometry")

			if st in (SUBTYPES[-3], SUBTYPES[-1]):
				col=11
				u3 ="11"
			else:
				col = SUBTYPES.index(st) + 5
				if st == SUBTYPES[2]:
					u3 = str(col)
				else:
					# Use random colors for region/lattice plots
					u3 ="(rnd($%d)%%%d)"%(col,colors)
				self("unset colorbox")

			if int(self.plot("labels",1)) and col<8:
				# Check version of engine
				if self.engine.version() < 4.2:
					# Convert to notify
					messagebox.showwarning( \
						"No Material labels", \
						"Gnuplot version 4.2 or higher is required " \
						"for plotting material labels. " \
						"Option will be disabled.")
					self.plot_labels.set(0)
				else:
					font = self.plot("labelfont")
					if font is not None:
						font = "font '%s' "%(font)
					else:
						font = ""
					extra = ",'' ind 2 %s:(0):%d w labels " \
						"offset 0,.3 point ps 0.7 pt 11 " \
						"tc lt 0 %snotitle"%(using, col, font)

			# FIXME if one day I convert to with image (not working until 4.2.3)
			#if self.engine.version() > 4.2:
#			if self.engine.version() >= 5.0:
#				return self._plot("'%s.dat' ind 1 %s:%s with image notitle," \
#					"'' ind 0 %s w l lt -1 lw 1 notitle%s" \
#					% (self.plot.name, using, u3, using, extra))
#			else:
			return self._splot("'%s.dat' ind 1 %s:%s notitle," \
					"'' ind 0 %s:(0) w l ls 1 notitle%s" \
					% (self.plot.name, using, u3, using, extra))

	# ----------------------------------------------------------------------
	# Generate the plot
	# ----------------------------------------------------------------------
	def _usr1DPlot(self):
		# --- Loop over detectors to build plot string ---
		for i in range(int(self.plot("ndetectors", 0))):
			if int(self.plot("show.%d" % (i), 1)) != 1:
				continue  # Skip graph
			filename = self.plot("file.%d" % (i), "")
			det = int(self.plot("det.%d" % (i), -1))
			if filename == "" or det < 0: continue

			# block to display
			try:
				block = int(self.plot("block.%d"%(i),0))
				if block > 0:
					block -= 1
					det += " every :::%d::%d "%(block,block)
			except:
				pass

			# load data
			try:
				name, xlow, xhigh, y, yerr = Data.tabLis(filename, det, block)
			except:
				# error cannot open file...
				continue

			# ===================================
			# block to display
			# try:
			#	block = int(self.plot("block.%d"%(i),0))
			#	if block > 0:
			#		block -= 1
			#		det += " every :::%d::%d "%(block,block)
			# except:
			#	pass
			# cmd += "'%s' %s " % (filename, det)
			# ===================================

			xlog = int(self.plot("xlog", 0))
			style = self.plot["with.%d" % (i)]
			addErrors = style == "histerror"  # Additional errors

			# Default using (middle of bin)
			if xlog:
				X = numpy.sqrt(xlow * xhigh)
			else:
				X = (xlow + xhigh) / 2.0

			if style in ("histogram", "histerror"):
				kwargs = {'drawstyle': 'steps-mid'}
				X = xlow
				usingex = ""

			yval = int(self.plot["y.%d" % (i)])
			if yval == 1:  # <X>*Y
				Y = X * y
				EY = X * y * yerr / 100.0
			elif yval == 2:  # DX*Y
				Y = (xhigh - xlow) * y
				EY = (xhigh - xlow) * y * yerr / 100.0

			else:  # Leave default Y
				Y = y
				EY = y * yerr / 100.0


			pyplot.plot(X, Y, **kwargs)
			continue

			# Normalization of x-axis
			xnorm = self.plot["xnorm.%d"%(i)]
			if xnorm:
				usingx  = self.plotNorm(xnorm, usingx,  "x%d"%(i))
				if usingex:
					usingex = self.plotNorm(xnorm, usingex, "ex%d"%(i))
				if addErrors:
					usingxe = self.plotNorm(xnorm, usingxe, "xe%d"%(i))

			# Normalization and default using
			norm = self.plot["norm.%d"%(i)]
			if norm:
				usingy  = self.plotNorm(norm, usingy,  "y%d"%(i))
				usingey = self.plotNorm(norm, usingey, "ey%d"%(i))

			cmd += "us %s "% \
				(self._plotUsing(style, usingx,
					usingy, usingex, usingey))

			# Plotting style
			cmd += self._plotStyle(i, style)

			# Label
			title = self.plot("name.%d"%(i),"")
			key   = int(self.plot("key.%d"%(i),1))
			if title.strip()=="" or key==0:
				cmd += "notitle"
			else:
				cmd += "t '%s'" % (title)

			if addErrors:
				cmd += ",'' %s us %s "% \
					(det, self._plotUsing("errorbars", usingxe,
						usingy, "", usingey))
				cmd += self._plotStyle(i, "errorbars")
				cmd += "notitle"

			# return self._plot(cmd)

	# ----------------------------------------------------------------------
	# USR-2D plot
	# ----------------------------------------------------------------------
	def _usr2DPlot(self):
		return
		# Normalization and default using
		using = "us 1:2:%s"%(self.plotNorm(self.plot["norm"],"($1*$3)"))

		# Color band
		self("set hidden3d")
		self.plotColorBand("at s")

		return self._splot("'%s' ind %s %s with lines lt -1 notitle" \
			% (self.plot["datafile"], self.plot["det"], using))

	# ----------------------------------------------------------------------
	# Usrbin plot
	# ----------------------------------------------------------------------
	def _usrbinPlot(self):
		return
		err = Plot._usrbinPlot(self)
		#if err is not None: return err
		if self._det is None: return err

		regionPlot = self._det.type in (2, 8, 12, 18)
		b10  = self._det.type % 10
		proj = self.plot("proj","Z")
		swap = int(self.plot("swap",0))
		norm = self.plot["norm"]
		font = ""
		try:
			h = USRBIN_HIST.index(self.plot["hist"])
		except ValueError:
			h = 0

		# Styles
		self("set style line 1 lt -1 lw 1")

		# --- 2D Histogram ---
		if h == 0:
			# Geometry information
			geo  = self.plot["geo"]
			axes = self.plot("axes",COORDAUTO[0])
			plotgeo = ""
			if geo == GEO_NO or err:
				plotgeo = ""
			else:
				if geo == GEO_AUTO:
					name = self.plot.name + ".geo"
				else:
					name = geo

				if axes == COORDAUTO[0]:
					if   proj == "X":
						u1, u2 = 5, 4
					elif proj == "Y":
						u1, u2 = 5, 3
					else:
						u1, u2 = 3, 4

					if regionPlot:
						# Find plot information
						for p in self.project.plots:
							if p.name == geo:
								axes = p.get("coord","U-V")
								u1 = ord(axes[ 0])-ord('U')+1
								u2 = ord(axes[-1])-ord('U')+1
								if u1>3: u1 -= 1
								if u2>3: u2 -= 1
#								font = p.get("labelfont")
#								if font is not None:
#									font = "font '%s' "%(font)
								break

				else:
					u1 = ord(axes[ 0])-ord('U')+1
					u2 = ord(axes[-1])-ord('U')+1
					if u1>3: u1 -= 1
					if u2>3: u2 -= 1
				#if self.engine.version() > 4.2:
				#if self.engine.version() > 5.0:
				#	plotgeo = ",'%s.dat' ind 0 us %d:%d w l ls 1 notit" \
				#		% (name, u1, u2)
				#else:
				plotgeo = ",'%s.dat' ind 0 us %d:%d:(0.) w l ls 1 notit" \
						% (name, u1, u2)

			# Plot region or normal usrbin plot
			extra = ""
			if regionPlot:
				if geo==GEO_NO or geo==GEO_AUTO:
					# Convert to notify
					messagebox.showwarning("Warning - Plotting by region",
						"Plot of Region or Special binning has been requested,\n" \
						"however no geometry to plot-onto has been selected\n" \
						"-Auto- or -No- plot cannot not work with this type of binning.\n"\
						"Please create first a geometry plot")
					return "Plot by region..."
				ind   = "ind 1"
				using = "%d:%d:"%(u1,u2)

				# label plot
				extra = ",'' ind 2 us %s(0):6 w labels " \
					"offset 0,.3 point ps 0.7 pt 11 " \
					"tc lt 0 %snotitle"%(using,font)
				# Color band
				if int(self.plot("errors",0)) == 0:
					using += self.plotNorm(norm,7)
				else:
					using += self.plotNorm(norm,8)

			else:
				ind   = ""

				# Plotting limits
				xmin,xmax = "",""
				ymin,ymax = "",""
				if (b10==1 or b10==7) and proj=="Z":
					if swap:
						using = "($2*sin($1)):($2*cos($1)):"
					else:
						using = "($2*cos($1)):($2*sin($1)):"
				#elif b10==2:
				#	pass
				#elif b10==8:
				#	pass
				else:
					if (swap and proj!="Z") or (not swap and proj=="Z"):
						using = "2:1:"
					else:
						using = "1:2:"

#					# Find plotting limits, if defined by the user
#					if proj=="X":
#						axis1 = "z"
#						axis2 = "y"
#						xmin = self.plot["zmin"]
#						xmax = self.plot["zmax"]
#						ymin = self.plot["ymin"]
#						ymax = self.plot["ymax"]
#					elif proj=="Y":
#						axis1 = "z"
#						axis2 = "x"
#						xmin = self.plot["zmin"]
#						xmax = self.plot["zmax"]
#						ymin = self.plot["xmin"]
#						ymax = self.plot["xmax"]
#					else:
#						axis1 = "y"
#						axis2 = "x"
#
#					if using[0]=="1":
#						xaxis = axis1
#						yaxis = axis2
#					else:
#						xaxis = axis2
#						yaxis = axis1
#
#					# --------------------------------------
#					# avoid rounding errors add a small gap
#					# --------------------------------------
#					def _up(x):
#						try:
#							f = float(x)
#							if   f>0.0: return str(f*1.0000001)
#							elif f<0.0: return str(f*0.9999999)
#						except:
#							pass
#						return x
#
#					# --------------------------------------
#					def _dn(x):
#						try:
#							f = float(x)
#							if   f>0.0: return str(f*0.9999999)
#							elif f<0.0: return str(f*1.0000001)
#						except:
#							pass
#						return x
#
#					xmin = _dn(self.plot["%smin"%(xaxis)])
#					xmax = _up(self.plot["%smax"%(xaxis)])
#					ymin = _dn(self.plot["%smin"%(yaxis)])
#					ymax = _up(self.plot["%smax"%(yaxis)])

#				if xmin!="" or xmax!="":
#					self("set xrange [%s:%s]"%(xmin,xmax))
#				if ymin!="" or ymax!="":
#					self("set yrange [%s:%s]"%(ymin,ymax))

				# Color band
				if int(self.plot("errors",0)) == 0:
					using += self.plotNorm(norm,3)
				else:
					using += self.plotNorm(norm,4)

			# -- Common plotting --
			self.plotColorBand()
			#if self.engine.version() >= 5.0:
			#	err2 = self._plot("'%s.dat' %s us %s with image notitle %s%s" \
			#		% (self.plot.name, ind, using, plotgeo, extra))
			#else:
			err2 = self._splot("'%s.dat' %s us %s notitle %s%s" \
					% (self.plot.name, ind, using, plotgeo, extra))
			return err2 or err

		# --- 1D Histogram ---
		elif h == 1:
			if regionPlot:
				self("set xtics rotate by 90")
				y  = self.plotNorm(norm,4)
				ey = self.plotNorm(norm,"($4*$5/100)")
				return self._plot("'%s.dat' ind 0 us 0:%s:xticlabels(3) with boxes not" \
					%(self.plot.name, y))
			else:
				style = self.plot["with.0"]
				if style in ("histogram", "histerror"):
					x  = 1
				else:
					x  = "(($1+$2)/2.0)"
				ex = "(0)"
				y  = self.plotNorm(norm, 3)
				ey = self.plotNorm(norm, "($3*$4/100)")
				if style=="histerror":
					return self._plot("'%s.dat' us %s %s notitle,'' us %s %s notitle" \
						% (self.plot.name,
						   self._plotUsing(style,x,y,ex,ey),
						   self._plotStyle(0),
						   self._plotUsing("errorbars","(($1+$2)/2.0)",y,ex,ey),
						   self._plotStyle(0,"errorbars")))
				else:
					return self._plot("'%s.dat' us %s %s notitle" \
						% (self.plot.name,
						   self._plotUsing(style,x,y,ex,ey),
						   self._plotStyle(0)))

		# --- 1D Max + ...
		else:
			style = self.plot["with.0"]
			if style in ("histogram", "histerror"):
				x  = 1
			else:
				x  = "(($1+$2)/2.0)"
			ex = "(0)"
			ind = max(0,h-3)
			if ind==0:	# Max trace
				y  = self.plotNorm(norm, 3)
				ey = self.plotNorm(norm, "($3*$4/100)")
				if style=="histerror":
					return self._plot("'%s.dat' ind %d us %s %s notitle,'' ind %d us %s %s notitle" \
						% (self.plot.name,
						   ind,
						   self._plotUsing(style,x,y,ex,ey),
						   self._plotStyle(0),
						   ind,
						   self._plotUsing("errorbars","(($1+$2)/2.0)",y,ex,ey),
						   self._plotStyle(0,"errorbars")))
				else:
					return self._plot("'%s.dat' ind %d us %s %s notitle" \
						% (self.plot.name,
						   ind,
						   self._plotUsing(style,x,y,ex,ey),
						   self._plotStyle(0)))

			else:	# horizontal/vertical trace
				# make plotstyle
				axes = self.plot("axes.0","x1y1")
				if axes != "x1y1":
					plotstyle = "axes %s " % (axes)
				else:
					plotstyle = ""

				plotstyle += "w yerrorbars lt %s" % (self.plot("lt.0","1"))
				if self.engine.version() >= 4.2:
					lc = self.plot["lc.0"]
					if lc != "": plotstyle += "lc rgb '%s' "%(lc)

				lw = "lw %s" % (self.plot("lw.0","2"))

				return self._plot("'%s.dat' ind %d us %s:3:4:5 %s %s notitle" \
						% (self.plot.name, ind, x, plotstyle, lw))

	# ----------------------------------------------------------------------
	# RESNUCLEi plot
	# ----------------------------------------------------------------------
	def _resnucleiPlot(self):
		return
		err = Plot._resnucleiPlot(self)
		if err is not None: return err

		plotType = self.plot["plot"]
		plotId   = RESNUCLEI_TYPES.index(plotType)

		# Styles
		self("set style line 1 lt -1 lw 1")
		self("set style line 2 lt -1 lw 0.5")

		# --- histogram plotting ---
		norm = self.plot["norm"]
		if plotType in (RESNUCLEI_A, RESNUCLEI_Z):
			# --- 1D plotting ---
			index = "ind %d" % (plotType == RESNUCLEI_A)

			x  = "($1-0.5)"
			ex = "(0)"
			y  = self.plotNorm(norm,3)
			ey = self.plotNorm(norm,"($3*$4/100)")
			style = self.plot["with.0"]

			if style == "histerror":
				return self._plot("'%s.dat' %s us %s %s notitle,'' %s us %s %s notitle" % \
					(self.plot.name,
					 index,
					 self._plotUsing(style,x,y,ex,ey),
					 self._plotStyle(0),
					 index,
					 self._plotUsing("errorbars","1",y,ex,ey),
					 self._plotStyle(0,"errorbars")))
			else:
				return self._plot("'%s.dat' %s us %s %s notitle" % \
					(self.plot.name,
					 index,
					 self._plotUsing(style,x,y,ex,ey),
					 self._plotStyle(0)))

		else:
			# --- 2D plotting ---
#			self("set hidden3d")
			self.plotColorBand() #"at s")

			# Assume format x/y
			using = "us "
			for p in plotType.split("/"):
				if   p=="Z":
					using += "($1-0.5):"
				elif p=="A":
					using += "($2-0.5):"
				elif p=="N":
					using += "($2-$1)>0?($2-$1):1/0:"

			# Normalization and default using
			using += self.plotNorm(norm,3)

			return self._splot("'%s.dat' ind 2 %s notitle"% (self.plot.name, using))

	# ----------------------------------------------------------------------
	def _userdumpPlot(self):
		return
		err = Plot._userdumpPlot(self)
		if err: return err

		# --- Geometry ---
		geo  = self.plot["geo"]
		axes = self.plot["geoaxes"]
		if geo == GEO_NO:
			plotgeo = ""
		else:
			name = geo
			u1 = ord(axes[ 0])-ord('U')+1
			u2 = ord(axes[-1])-ord('U')+1
			if u1>3: u1 -= 1
			if u2>3: u2 -= 1
			plotgeo = ",'%s.dat' ind 0 us %d:%d w l lt -1 lw 1 notit" \
				% (name, u1, u2)

		case = self.plot["case"]
		if case == "Continuous":
			return self._splot("'%s.dat' us 4:5:6 w l"%(self.plot.name))
		elif case == "Point":
			pass
		elif case == "Source":
			column = SOURCE_HIST.index(self.plot["x"])

			# Limits
			xmin  = self.plot.get("xmin","")
			if xmin == "":
				xmin = self.dmin[column]
			else:
				xmin = float(xmin)

			xmax  = self.plot.get("xmax","")
			if xmax == "":
				xmax = self.dmax[column]
			else:
				xmax = float(xmax)

			# --- Histogram or Scatter
			if self.plot["type"] == HISTOGRAM:
				log   = int(self.plot["xlog"])
				try: bins = float(self.plot["bins"])
				except: bins = 50.0
				if log:
					if xmin <= 0.0: xmin = 1e-10
					if xmax <= 0.0: xmax = 1.
					step  = (math.log10(xmax) - math.log10(xmin)) / bins
					func = "logbin"
				else:
					step  = (xmax - xmin) / bins
					func = "bin"

				if int(self.plot.get("weight",1)):
					weight = "3"
				else:
					weight = "(1)"

				# Using...
				using = "(%s($%d,%g,%g)):%s" % \
					(func, column+1, self.dmin[column],
					 step, weight)
				self("set style fill solid 0.5 noborder")
				return self._plot("'%s.dat' us %s " \
						"smooth freq with boxes lt 7 notitle," \
						"'' us %s smooth freq with histeps lt -1 notitle" % \
						(self.plot.name, using, using))

			else:	# Scatter plot
				column2 = SOURCE_HIST.index(self.plot["y"])

				return self._plot("'%s.dat' us %d:%d with points %s" % \
						(self.plot.name, column+1, column2+1, plotgeo))
