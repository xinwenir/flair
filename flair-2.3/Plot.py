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
import sys
import time
import math
import string
import random
import subprocess

import pipe
import bmath
import Utils
import tkFlair

import Data
import Input
import Project
import Palette

DEFAULT_FORMAT = ".eps"

GEO_NO    = "-No-"
GEO_AUTO  = "-Auto-"

PLOT_TYPES = [ "Geometry",
		"USR-1D",
		"USR-2D",
		"USRBIN",
		"USERDUMP",
		#"Output",	<- Thresholds, WW, Energy balance, etc...
		"RESNUCLE" ]

#AXES = ["x","y","z","cb","x2","y2"]

COORD     = ["U-V", "V-U", "X-Y", "X-Z", "Y-Z", "Y-X", "Z-X", "Z-Y" ]
COORDAUTO = ["Auto"] + COORD

AXES = [	"x1y1",
		"x1y2",
		"x2y1",
		"x2y2" ]

STYLE = [
		"errorbars",		# X
		"errorlines",		# X
		"histerror",		# histogram with error
		"histogram",		# Xlow
		"lines",		# X
		"linespoints",		# X
		"points",
		"xyerrorbars"]
STYLE_DEF = 3

VALUE = [	"Y",			# 0
		"<X>*Y",		# 1
		"DX*Y" ]		# 2

X_SCALE =  [	"",
		"1./eV",
		"1./keV",
		"1./MeV",
		"1./GeV",
		"1./um",
		"1./mm",
		"1./cm",
		"1./dm",
		"1./m",
		"1./km",
#		"in",
#		"ft",
#		"mi",
		]

PALETTE = [	"FLUKA",
		"FLUKA2",
		"gnuplot",
		"defined",
		"grey",
		"rainbow",
		"rgb",
		"red-yellow",
		"geometry",
		"field",
		"Matlab"]

LINE_TYPES = ["", "solid", "dashed"]
FONTS      = ["",
		"Arial",
		"Courier",
		"Helvetica",
		"Sans",
		"Serif",
		"Symbol",
		"Times" ]
SIZE       = ["","8","9","10","11","12","14","16","18","20","22","24","26","28","36","48","72"]

#                 0                   1            2           3            4           5
USRBIN_HIST = ['2D Projection', '1D Projection', '1D Max', '1D Min', '1D Trace-H', '1D Trace-V']


SUBTYPES = ["Geometry", "Region", "Material", "Lattice", "Field Intensity", "Field Vector", "Field" ]
#ERRGEO = "Error in geometry or plotting region is empty\nTry another region."

RESNUCLEI_Z   = "Z"
RESNUCLEI_A   = "A"
RESNUCLEI_Z_A = "Z/A"
RESNUCLEI_Z_N = "Z/N"
RESNUCLEI_A_Z = "A/Z"
RESNUCLEI_N_Z = "N/Z"

RESNUCLEI_TYPES  = [ RESNUCLEI_Z, RESNUCLEI_A, RESNUCLEI_Z_A, RESNUCLEI_Z_N, RESNUCLEI_A_Z, RESNUCLEI_N_Z ]
RESNUCLEI_TITLE  = ["Atomic Number",
		"Mass Number",
		"Mass vs Atomic Number",
		"Neutron vs Atomic Number",
		"Atomic vs Mass Number",
		"Atomic vs Neutron Number" ]
#RESNUCLEI_STABLE = ["No", "Lines", "Dots"]

SOURCE_HIST = [ 'Particle',
		 'Total-Energy',
		 'Weight',
		 'x',
		 'y',
		 'z',
		 'tx',
		 'ty',
		 'tz',
		 'Kinetic-Energy',	# Are not in the file
		 'Momentum' ]

HISTOGRAM = "Histogram"
SCATTER   = "Scatter"

USERDUMP_CASE = ["Continuous", "Point", "Source"]

# -----------------------------------------------------------------------------
# perform a flood fill algorithm to determine regions
# -----------------------------------------------------------------------------
def _scan4Labels(hist):
	regions = []	# List of lists
			# Each row contains [reg/0, n/1, sumx/2, sumy/3]
	rows = len(hist)
	cols = len(hist[0])

	rnd = random.Random()
	rnd.seed(314159)
	pos = []
	pos_append = pos.append

	for row in range(rows):
		for col in range(cols):
			reg = hist[row][col]
			if reg < 0: continue

			# Create a new entry
			index = -len(regions)-1
			n, sx, sy  = 1, col, row
			stack = [(col,row)]
			stack_append = stack.append
			del pos[:]

			while len(stack)>0:
				xo, y = stack.pop()
				# Scan to the left
				x = xo+1
				yup = y-1
				ydn = y+1
				while x>0:
					x -= 1
					if hist[y][x] != reg: break
					hist[y][x] = index
					n  += 1
					sx += x
					sy += y
					pos_append((x,y))

					# Check up/down
					if yup>=0 and hist[yup][x]==reg:
						# Push up value
						# XXX Check if neighbor is pushed!
						stack_append((x, yup))
					if ydn<rows and hist[ydn][x]==reg:
						# Push down value
						# XXX Check if neighbor is pushed!
						stack_append((x, ydn))

				# Scan to the right
				x = xo
				while True:
					x += 1
					if x==cols or hist[y][x] != reg: break
					hist[y][x] = index
					n  += 1
					sx += x
					sy += y
					pos_append((x,y))

					# Check up/down
					if yup>=0 and hist[yup][x]==reg:
						# Push up value
						# XXX Check if neighbor is pushed!
						stack_append((x, yup))
					if ydn<rows and hist[ydn][x]==reg:
						# Push down value
						# XXX Check if neighbor is pushed!
						stack_append((x, ydn))

			# Ignore small regions
			if n<5: continue
			ax = int(float(sx)/float(n))
			ay = int(float(sy)/float(n))

			# Check neighboring materials
			SCAN=5
			ok = True
			for y in range(max(0, ay-SCAN), min(ay+SCAN,rows)):
				for x in range(max(0, ax-SCAN), min(ax+SCAN,cols)):
					if hist[y][x] != index:
						ok = False
						break
				if not ok: break
			if ok:
				regions.append((reg, ax+0.5, ay+0.5))
			else:
				# Try few times to find a point not close
				# to anything else
				max_dist2 = 0
				max_x = -1
				max_y = -1
				if len(regions) == 0:
					# Select randomly
					r = rnd.randint(0, len(pos)-1)
					max_x,max_y = pos[r]
				else:
					for i in range(SCAN):
						# Select randomly
						r = rnd.randint(0, len(pos)-1)
						rx,ry = pos[r]
						min_dist2 = 10000000
						min_x = -1
						min_y = -1
						for rg,x,y in regions:
							d2 = (x-rx)**2 + (y-ry)**2
							if d2 <= min_dist2:
								min_dist2 = d2
								min_x = rx
								min_y = ry

						if min_dist2 > max_dist2:
							max_dist2 = min_dist2
							max_x = min_x
							max_y = min_y

				regions.append((reg, max_x+0.5, max_y+0.5))
	return regions

#===============================================================================
# Plot class
#===============================================================================
class Plot:
	def __init__(self, engine, project, plot, log=None):
		self.project = project
		self.engine  = engine
		self.plot    = plot
		self._log    = log
		self.engine.setLogger(log)

	# ----------------------------------------------------------------------
	def log(self, *kw):
		s = " ".join(map(str,kw))
		if self._log is None:
			sys.stdout.write("%s\n"%(s))
		else:
			self._log(s)

	# ----------------------------------------------------------------------
	def write(self, s):
		self.log(s)

	# ----------------------------------------------------------------------
	# Show plot on screen
	# ----------------------------------------------------------------------
	def show(self):
		self.engine.reset()
		self._startPlot()
		if self.plot.type == "Geometry":
			err = self._geometryPlot()
		elif self.plot.type == "USR-1D":
			err = self._usr1DPlot()
		elif self.plot.type == "USR-2D":
			err = self._usr2DPlot()
		elif self.plot.type == "USRBIN":
			err = self._usrbinPlot()
		elif self.plot.type == "RESNUCLE":
			err = self._resnucleiPlot()
		elif self.plot.type == "USERDUMP":
			err = self._userdumpPlot()
		else:
			err = None

		if err: self.log(err)

		# Keep plot information
		self.plot.time = time.time()
#		self.plotId = self.id
		return err

	# ----------------------------------------------------------------------
	# Save plot to file
	# ----------------------------------------------------------------------
	def save(self, ext):
		return self.plot.name+ext

	# ----------------------------------------------------------------------
	# Send plot on printer
	# ----------------------------------------------------------------------
	def hardcopy(self, printer):
		pass

	# ----------------------------------------------------------------------
	# Clean up files and variables
	# ----------------------------------------------------------------------
	def cleanup(self):
		if self.plot.type in ("Geometry", "USRBIN"):
			self.plot._what = [""] * 24
			if self.plot.type == "USRBIN":
				self.plot._inp = None

		elif self.plot.type == "USR-1D":
			pass

		elif self.plot.type == "USR-2D":
			pass

		elif self.plot.type == "RESNUCLE":
			self.plot._previous = None

		elif self.plot.type == "USERDUMP":
			self.plot._case    = -1
			self.plot._mgfile  = ""
			self.plot._limfile = ""

		for f in [self.plot.name+".dat", "gplevh.lim", "gplevh.npo", "gplevh.poi"]:
			try: os.remove(f)
			except OSError: pass

	# ----------------------------------------------------------------------
	# Perform plot
	# ----------------------------------------------------------------------
	def _plot(self, cmd):
		return None

	# ----------------------------------------------------------------------
	# Perform splot
	# ----------------------------------------------------------------------
	def _splot(self, cmd):
		return None

	# ======================================================================
	# Generate the plot
	# ----------------------------------------------------------------------
	def _geometryPlot(self):
		if self._getGeometryPlotgeom():
			return self._runFluka(self.plot.name)
		return None

	# ----------------------------------------------------------------------
	# Initialize plot
	# ----------------------------------------------------------------------
	def _startPlot(self):
		self.engine.title(self.plot["title"], self._plotOptions("title"))
		self.engine.grid(self.plot["grid"])
		self._setupAxes()

	# ----------------------------------------------------------------------
	# Setup plotting axes
	# ----------------------------------------------------------------------
	def _setupAxes(self):
		self.engine.label('x', self.plot["xlabel"], self._plotOptions("x"))
		self.engine.tics( 'x', self._plotOptions("xtics"))

		self.engine.label('y', self.plot["ylabel"], self._plotOptions("y"))
		self.engine.tics( 'y', self._plotOptions("ytics"))

		if self.plot["ztics"]:
			self.engine.label('z',  self.plot["zlabel"],  self._plotOptions("z"))
			self.engine.tics( 'z', self._plotOptions("ztics"))

		if self.plot["cbtics"]:
			self.engine.label('cb', self.plot["cblabel"], self._plotOptions("cb"))
			self.engine.tics( 'cb', self._plotOptions("cbtics"))

		if self.plot["x2tics"]:
			self.engine.label('x2', self.plot["x2label"], self._plotOptions("x2"))
			self.engine.tics( 'x2', self._plotOptions("x2tics"))
			self("set x2tics")

		if self.plot["y2tics"]:
			self.engine.label('y2', self.plot["y2label"], self._plotOptions("y2"))
			self.engine.tics( 'y2', self._plotOptions("y2tics"))
			self("set y2tics")

		self._plotRange("x")
		self._plotRange("y")
		self._plotRange("z")
		self._plotRange("cb")
		self._plotRange("x2")
		self._plotRange("y2")

	# ----------------------------------------------------------------------
	# Set Plot range and log/linear
	# ----------------------------------------------------------------------
	def _plotRange(self, axis):
		self.engine.log(axis, int(self.plot("%slog"%(axis),0)))
		amin = self.plot["%smin"%(axis)]
		amax = self.plot["%smax"%(axis)]
		self.engine.range(axis, amin, amax)

	# ----------------------------------------------------------------------
	# Generate the plot
	# ----------------------------------------------------------------------
	def _usr1DPlot(self):
		return None

	# ----------------------------------------------------------------------
	# Usrbin plot
	# ----------------------------------------------------------------------
	def _usrbinPlot(self):
		# Read detector
		usr = Data.Usrbin()
		self._det = None

		# Read file
		try:
			usr.readHeader(self.plot["datafile"])
		except (IOError, MemoryError):
			return "Error loading USRBIN file: %s"%(sys.exc_info()[1])

		# Using...
		d = int(self.plot.get("det",1))
		try:
			self._det = usr.detector[d-1]
		except IndexError:
			return "Detector %d not found"%(d)

		# Start plotting
		regionPlot = self._det.type in (2, 8, 12, 18)
		try:
			h = USRBIN_HIST.index(self.plot["hist"])
		except ValueError:
			h = 0

		if h<=1:
			if regionPlot:
				# Region plot
				err = self._mergeGeoData(usr, d, self._det)
				if err is not None: return err

			else:	# Normal plot with GPLEVBIN
				err = self._runLevbin(usr,d)
				if err is not None: return err

			# --- 2D Histogram ---
			if h==0 and self.plot["geo"]==GEO_AUTO:
				name = self.plot.name + ".geo"
				if self._getUsrbinPlotgeom(name, self._det):
					err = self._runFluka(name)
					if err is not None:
						self.plot._what = [""] * 24
						return err

		else:
			if regionPlot: return
			return self._runUsbmax(d)

	# ======================================================================
	# Run Fluka and prepare the output file
	# ----------------------------------------------------------------------
	def _runFluka(self, name):
		# Check executable
		try:
			os.stat(self.project.executable())
		except OSError:
			return "FLUKA executable '%s' not found."%(self.project.executable()),None,None

		# Write the inputfile
		title = self.plot["title"]
		if title == "": title = "-no title-"
		plotgeom = Input.Card("PLOTGEOM", self.plot._what, title, extra=title)

		inpname = self.project.writePlotgeomInput(plotgeom)

		# Execute fluka
		stime = time.time()
		cmd = self.project.runCmd(inpname, 1)
		aaa=list(cmd)
		bbb=""
		for i in range(len(aaa)):
			if isinstance(aaa[i],str):
				bbb+=aaa[i]
			else:
				bbb+=aaa[i].decode()

#		self.log("\nCmd: %s"%(" ".join(list(cmd))))
		self.log("\nCmd: %s"%(bbb))

		process = subprocess.Popen(cmd,
				stdout=subprocess.PIPE,
				stderr=subprocess.STDOUT)
		for line in process.stdout:
			self.log(line[:-1])
		pid = process.pid
		rc  = process.returncode
		if rc is None: rc=0
		self.log("... Fluka run: %g s"%(time.time()-stime))

		# Error?
		if rc != 0:
			return "Return Code=%d"%(rc),inpname,pid

		# Process output files
		# if everything went well we should have a file called
		# $[inpname}001_PLOTGEOM.STORE
		try:
			fpg = open("%s001_PLOTGEOM.STORE" % (inpname),"r")
		except:
			return "No %s001_PLOTGEOM.STORE file"%(inpname),inpname,pid

		try:
			fgeo = open("%s.dat"%(name),"w")
		except:
			return "Error opening file %s.dat for writing"%(self.plot.name),inpname,pid

		# Work out the origin information from the what's
		O = bmath.Vector((self.plot._what[7]+self.plot._what[10])/2.0,
				 (self.plot._what[8]+self.plot._what[11])/2.0,
				 (self.plot._what[9]+self.plot._what[12])/2.0)
		V = bmath.Vector(self.plot._what[13], self.plot._what[14], self.plot._what[15])
		U = bmath.Vector(self.plot._what[16], self.plot._what[17], self.plot._what[18])

		loP  = O - bmath.Vector(self.plot._what[7], self.plot._what[8], self.plot._what[9])
		extU = U.dot(loP)
		extV = V.dot(loP)

		# Skip first 7 lines
		for i in range(7):
			fgeo.write("#"+fpg.readline())

		TARG1 = " Cold storage";		len1 = len(TARG1)
		TARG2 = " Worm no";			len2 = len(TARG2)
		TARG3 = " ***MATERIAL***FIELD***";	len3 = len(TARG3)
		TARG4 = " ***MATERIAL***NAMES***";	len4 = len(TARG4)
		TARG5 = " ***REGION***NAMES***";	len5 = len(TARG5)

		nworm      = 0
		firstblock = 0	# Check that at least one worm have different size
		blocksize  = 0	# ... otherwise gnuplot thinkgs we are dealing with isolines

		self.log("Processing PLOTGEOM.STORE file")
		stime = time.time()
		while True:
			line = fpg.readline()
			if line=="": break
			elif line[0:len3] == TARG3: break
			elif line[0:len1] == TARG1 or line[0:len2] == TARG2:
				if nworm>0:
					# Write the 2 last points
					lx -= extU
					ly -= extV
					r = O + lx*U + ly*V
					fgeo.write("%.8g %.8g %.8g %.8g %.8g\n" \
						% (lx,ly,r[0],r[1],r[2]))

					mx -= extU
					my -= extV
					r = O + mx*U + my*V
					fgeo.write("%.8g %.8g %.8g %.8g %.8g\n\n" \
						% (mx,my,r[0],r[1],r[2]))
					blocksize += 2

					if firstblock == 0:
						firstblock = blocksize
					elif firstblock != blocksize:
						firstblock = -1

				nworm = 0
				blocksize = 0
				continue

			try:
				x, y = list(map(float, line.split()))
			except:
				return  "Processing PLOTGEOM",inpname,pid

			# Check if we have to forget mid point
			nworm += 1
			if nworm >= 3:
				dx1 = x  - lx
				dy1 = y  - ly
				dx2 = mx - lx
				dy2 = my - ly
				diff = (dx1*dx2 + dy1*dy2) - \
					(math.sqrt(dx1**2 + dy1**2) * math.sqrt(dx2**2 + dy2**2))

				if diff != 0.0:
					lx -= extU
					ly -= extV
					r = O + lx*U + ly*V
					fgeo.write("%.8g %.8g %.8g %.8g %.8g\n" \
						% (lx,ly,r[0],r[1],r[2]))
					blocksize += 1
					lx, ly = mx, my
				else:
					nworm -= 1
				mx, my = x, y
			elif nworm == 2:
				mx, my = x, y
			else:
				lx, ly = x, y

		# Write the 2 last points
		if nworm >= 1:
			lx -= extU
			ly -= extV
			r = O + lx*U + ly*V
			fgeo.write("%.8g %.8g %.8g %.8g %.8g\n" \
				% (lx,ly,r[0],r[1],r[2]))
			blocksize += 1

			if nworm >= 2:
				mx -= extU
				my -= extV
				r = O + mx*U + my*V
				fgeo.write("%.8g %.8g %.8g %.8g %.8g\n" \
					% (mx,my,r[0],r[1],r[2]))
				blocksize += 1
			else:
				mx, my = lx, ly

		if firstblock == blocksize:
			if blocksize == 0:
				return "Geometry Problem or empty plot at the specified location."\
					"Try to change the location or enlarge the plot. "\
					"For USRBIN change either the Use: to -No- or the Pos:", inpname, pid
			else:
				fgeo.write("%.8g %.8g %.8g %.8g %.8g\n" \
					% (mx,my,r[0],r[1],r[2]))
		self.log("... Worms:     %g s"%(time.time()-stime))
		stime = time.time()

		# Material section
		if line!="":
			fgeo.write("\n\n")

			# Split #Grid nx ny line
			grid = fpg.readline().split()
			nu = int(grid[1])
			nv = int(grid[2])

			# New 2010 release
			region2material = [0]*(len(self.project.input["REGION"])+20)

			# Histogram
			histogram = [[0]*nu for i in range(nv)]
			us = (2.0*extU)/nu
			vs = (2.0*extV)/nv
			v  = -extV
			#v  = -extV + vs/2.0	# with image
			for i in range(nv):
				u  = -extU
				#u  = -extU + us/2.0.0	# with image
				histrow = histogram[i]
				for j in range(nu):
					line = fpg.readline()
					reg = int(line[ 0: 5])
					mat = int(line[ 5:10])
					lat = int(line[10:15])
					bu  = float(line[15:27])
					bv  = float(line[27:39])
					b   = float(line[39:51])
					r   = O + u*U + v*V
					fgeo.write("%.8g %.8g %.8g %.8g %.8g %d %d %d %g %g %g\n" \
						% (u,v,r[0],r[1],r[2],reg,mat,lat,bu,bv,b))
					histrow[j] = reg
					try:
						region2material[reg] = mat
					except IndexError:	# It should come here for VOXELS
						region2material.extend([0] * (reg-len(region2material)+20))
						region2material[reg] = mat

					u += us
				fgeo.write("\n")
				v += vs
			self.log("... Histogram: %g s"%(time.time()-stime))
			stime = time.time()

			line = fpg.readline()

			try:
				del self.plot._materialList[:]
			except AttributeError:
				self.plot._materialList = []

			if line[0:len4] == TARG4:
				# New 2010 release
				# materials
				while True:
					line = fpg.readline()
					if len(line)==0:
						break	# XXX should never arrive here ###
					if line[0:len5] == TARG5:
						break
					# should be in serial order!!!!
					try:
						self.plot._materialList.append(line.split()[1])
					except IndexError:
						self.plot._materialList.append("?")

				# regions
				regionList = []
				while True:
					line = fpg.readline()
					if len(line)==0:
						break
					# should be in serial order!!!!
					try:
						regionList.append(line.split()[1])
					except IndexError:
						regionList.append("?")
			else:
				# Old <= 2008 fluka
				# XXX XXX should be removed in the future XXX XXX
				self.project.input.preprocess()
				cardlist = self.project.input["REGION"]
#				cardlist.sort(Input.Card.cmpPos)
				cardlist.sort(key=lambda card: card.pos())
				regionList = [x.sdum() for x in cardlist if not x.ignore() and x.sdum()!="&"]
				self.plot._materialList.extend([x for x in self.project.input.materialList()])
				if "VOXELS" in self.project.input.cards:
					regionList.append("VOXEL")
					fromIdx = len(regionList)
					regionList.extend(["?"] * (len(region2material)-fromIdx+1))
					voxel = 1
					for reg in range(fromIdx, len(region2material)):
						regionList[reg] = "VOXEL%03d"%(voxel)
					voxel += 1

			# Scan histograms for labels
			first = True
			for reg,x,y in _scan4Labels(histogram):
				if first:
					fgeo.write("\n\n")
					first = False
				u = x*us-extU
				v = y*vs-extV
				r = O + u*U + v*V
				try:
					mat = str(self.plot._materialList[region2material[reg]-1]).replace("_",r"\\_")
					fgeo.write("%g %g %g %g %g %s %s\n" \
						% (u,v, r[0],r[1],r[2], regionList[reg-1], mat))
				except IndexError:
					if reg >= len(regionList):
						sreg = str(reg)
					else:
						sreg = regionList[reg-1]
					mat = region2material[reg]
					if mat >= len(self.plot._materialList):
						smat = str(mat)
					else:
						smat = self.plot._materialList[mat-1]
					smat = smat.replace("_",r"\\_")
					fgeo.write("%g %g %g %g %g %s %s\n" \
						% (u,v, r[0],r[1],r[2], sreg, smat))
			self.log("... Labels:    %g s"%(time.time()-stime))

		# Close files and do the plot
		fgeo.close()
		fpg.close()

		# Clean up
		if Project.cleanup:
			self.project.deleteFiles(inpname, pid)
			try:
				os.remove(inpname+".inp")
			except OSError:
				pass
		return None

	# ----------------------------------------------------------------------
	# Get plotgeom info and return True if it has changed
	# ----------------------------------------------------------------------
	def _getGeometryPlotgeom(self):
		try:
			what = self.plot._what
		except AttributeError:
			self.plot._what = [""]*24
			what = self.plot._what

		# Get variables
		O = bmath.Vector(list(map(float,self.plot["origin"].split())))
		U = bmath.Vector(list(map(float,self.plot["basisu"].split())))
		V = bmath.Vector(list(map(float,self.plot["basisv"].split())))
		extU,extV = list(map(float, self.plot["extends"].split()))
		loP  = -extU*U -extV*V + O
		hiP  =  extU*U +extV*V + O

		try:
			nu = float(self.plot.get("nu", 200.0))
		except ValueError:
			nu = 200
		try:
			nv = float(self.plot.get("nv", 200.0))
		except ValueError:
			nv = 200

		# Run fluka only if something has changed
		datfile = self.plot.name+".dat"
		try:
			# Check the output data file
			os.stat(datfile)
			changed = self.project.inputModifiedTime() > self.plot.time
		except OSError:
			changed = True

		if not changed: changed = what[ 2] != 1-int(self.plot.get("boundaries",0))
		if not changed: changed = what[ 7] != loP.x()
		if not changed: changed = what[ 8] != loP.y()
		if not changed: changed = what[ 9] != loP.z()
		if not changed: changed = what[10] != hiP.x()
		if not changed: changed = what[11] != hiP.y()
		if not changed: changed = what[12] != hiP.z()
		if not changed: changed = what[13] != V.x()
		if not changed: changed = what[14] != V.y()
		if not changed: changed = what[15] != V.z()
		if not changed: changed = what[16] != U.x()
		if not changed: changed = what[17] != U.y()
		if not changed: changed = what[18] != U.z()
		if not changed: changed = what[21] != nu
		if not changed: changed = what[22] != nv
		if not changed: changed = Utils.isNewer(self.project.input,
						self.plot.name+".dat")

		if not changed: return False

		what[ 0] = "FORMAT"
		what[ 2] = 1-int(self.plot.get("boundaries",0))	# Plot regions

		what[ 7] = loP.x()
		what[ 8] = loP.y()
		what[ 9] = loP.z()
		what[10] = hiP.x()
		what[11] = hiP.y()
		what[12] = hiP.z()

		what[13] = V.x()
		what[14] = V.y()
		what[15] = V.z()
		what[16] = U.x()
		what[17] = U.y()
		what[18] = U.z()

		what[19] = 1.0
		what[20] = 1.0
		what[21] = nu
		what[22] = nv

		return True

	# ----------------------------------------------------------------------
	# Get plotgeom info and return True if it has changed
	# ----------------------------------------------------------------------
	def _getUsrbinPlotgeom(self, name, det):
		b10  = det.type % 10
		proj = ord(self.plot.get("proj","Z")) - ord("X")
		pos  = self.plot["geopos"]
		try:
			what = self.plot._what
		except AttributeError:
			self.plot._what = [""]*24
			what = self.plot._what

		# Get variables
		if   b10==1 or b10==7:
			# Cylindrical x=R, y=Phi, z=z
			# ymin=X, yhigh=Y
			# XXX bug in fluka does not keep the X,Y in the binning
			X = 0
			Y = 0
			if   proj==0: # X
				return False	# No geometry plot when Z-Phi
			elif proj==1: # Y
				# More tricky pos should be the rotation
				loP = bmath.Vector(X, Y, det.zlow)
				hiP = bmath.Vector(X+det.xhigh,  Y, det.zhigh)
				U = bmath.Vector.X
				V = bmath.Vector.Z
			elif proj==2: # Z
				if pos!="":
					try:
						z = float(pos)
					except ValueError:
						z = (det.zlow+det.zhigh)/2
				else:
					z = (det.zlow+det.zhigh)/2
				loP = bmath.Vector(-det.xhigh + X,
						   -det.xhigh + Y,
						    z)
				hiP = bmath.Vector( det.xhigh + X,
						    det.xhigh + Y,
						    z)
				U = bmath.Vector.X
				V = bmath.Vector.Y
		elif b10==2:
			return False
		elif b10==8:
			return False
		else:
			loP = bmath.Vector(det.xlow,  det.ylow,  det.zlow)
			hiP = bmath.Vector(det.xhigh, det.yhigh, det.zhigh)

			# Override limits if defined by user
			try: loP[0] = float(self.plot["bxmin"])
			except: pass
			try: loP[1] = float(self.plot["bymin"])
			except: pass
			try: loP[2] = float(self.plot["bzmin"])
			except: pass

			try: hiP[0] = float(self.plot["bxmax"])
			except: pass
			try: hiP[1] = float(self.plot["bymax"])
			except: pass
			try: hiP[2] = float(self.plot["bzmax"])
			except: pass

			# Add a small margin to avoid numerical precision
			# errors from the PLOTGEOM
			for i in range(3):
				loP[i] -= abs(loP[i]+1.0)*1e-6
				hiP[i] += abs(hiP[i]+1.0)*1e-6

			# Unless there are limits set...
			if pos != "":
				loP[proj] = hiP[proj] = float(pos)
			else:
				loP[proj] = hiP[proj] = (loP[proj]+hiP[proj])/2

			if   proj == 0:	# X
				U = bmath.Vector.Y
				V = bmath.Vector.Z
			elif proj == 1: # Y
				U = bmath.Vector.X
				V = bmath.Vector.Z
			elif proj == 2: # Z
				U = bmath.Vector.X
				V = bmath.Vector.Y

		changed = self.project.inputModifiedTime() > self.plot.time
		if not changed: changed = what[ 7] != loP.x()
		if not changed: changed = what[ 8] != loP.y()
		if not changed: changed = what[ 9] != loP.z()
		if not changed: changed = what[10] != hiP.x()
		if not changed: changed = what[11] != hiP.y()
		if not changed: changed = what[12] != hiP.z()
		if not changed: changed = what[13] != V.x()
		if not changed: changed = what[14] != V.y()
		if not changed: changed = what[15] != V.z()
		if not changed: changed = what[16] != U.x()
		if not changed: changed = what[17] != U.y()
		if not changed: changed = what[18] != U.z()
		if not changed: changed = Utils.isNewer(self.project.input,
						name+".dat")

		if not changed: return False

		what[ 0] = "FORMAT"
		what[ 2] = 1.0	# Plot materials

		what[ 7] = loP.x()
		what[ 8] = loP.y()
		what[ 9] = loP.z()
		what[10] = hiP.x()
		what[11] = hiP.y()
		what[12] = hiP.z()

		what[13] = V.x()
		what[14] = V.y()
		what[15] = V.z()
		what[16] = U.x()
		what[17] = U.y()
		what[18] = U.z()

		what[19] = 1.0
		what[20] = 1.0
		what[21] = 2.0
		what[22] = 2.0

		return True

	# ----------------------------------------------------------------------
	# Run the gplevbin program to process the data
	# ----------------------------------------------------------------------
	def _runLevbin(self, usr, d):
		if d<=0 or d>len(usr.detector):
			return "Requested detector out of range"
		det = usr.detector[d-1]

		# Find limits and bins
		# --- X ---
		xmin = self.plot["bxmin"]
		xmax = self.plot["bxmax"]
		if xmin=="" and xmax=="":
			limx = ""
			xmin = det.xlow
			xmax = det.xhigh
		else:
			try: xmin = float(xmin)
			except ValueError: xmin = det.xlow
			try: xmax = float(xmax)
			except ValueError: xmax = det.xhigh
			if xmin<det.xlow:  xmin = det.xlow
			if xmax>det.xhigh: xmax = det.xhigh
			limx = "%g %g"%(xmin,xmax)
		try:
			xrebin = int(self.plot.get("xrebin",1))
			if xrebin<=0: xrebin=1
		except:
			xrebin = 1
		dx = (det.xhigh - det.xlow) / float(det.nx) * float(xrebin)
		if abs(dx)<1e-19:
			self.log("Too small bin in x",dx)
			return False
		nx = str(int(round((xmax-xmin) / dx)))

		# --- Y ---
		ymin = self.plot["bymin"]
		ymax = self.plot["bymax"]
		if ymin=="" and ymax=="":
			limy = ""
			ymin = det.ylow
			ymax = det.yhigh
		else:
			try: ymin = float(ymin)
			except ValueError: ymin = det.ylow
			try: ymax = float(ymax)
			except ValueError: ymax = det.yhigh
			if ymin<det.ylow:  ymin = det.ylow
			if ymax>det.yhigh: ymax = det.yhigh
			limy = "%g %g"%(ymin,ymax)
		try:
			yrebin = int(self.plot.get("yrebin",1))
			if yrebin<=0: yrebin=1
		except:
			yrebin = 1
		dy = (det.yhigh - det.ylow) / float(det.ny) * float(yrebin)
		if abs(dy)<1e-19:
			self.log("Too small bin in y",dy)
			return False
		ny = str(int(round((ymax-ymin) / dy)))

		# --- Z ---
		zmin = self.plot["bzmin"]
		zmax = self.plot["bzmax"]
		if zmin=="" and zmax=="":
			limz = ""
			zmin = det.zlow
			zmax = det.zhigh
		else:
			try: zmin = float(zmin)
			except ValueError: zmin = det.zlow
			try: zmax = float(zmax)
			except ValueError: zmax = det.zhigh
			if zmin<det.zlow:  zmin = det.zlow
			if zmax>det.zhigh: zmax = det.zhigh
			limz = "%g %g"%(zmin,zmax)
		try:
			zrebin = int(self.plot.get("zrebin",1))
			if zrebin<=0: zrebin=1
		except:
			zrebin = 1
		dz = (det.zhigh - det.zlow) / float(det.nz) * float(zrebin)
		if abs(dz)<1e-19:
			self.log("Too small bin in z",dz)
			return False
		nz = str(int(round((zmax-zmin) / dz)))

		p = self.plot.get("proj","Z")
		try:
			h = USRBIN_HIST.index(self.plot["hist"])
		except ValueError:
			h = 0
		if h == 0:
			if   p == "X": nx = "1"
			elif p == "Y": ny = "1"
			else:          nz = "1"
		elif h == 1:
			if   p == "X": ny,nz = "1","1"
			elif p == "Y": nx,nz = "1","1"
			else:          nx,ny = "1","1"
		else:
			# Special type... not handled now
			nz = "1"
			pass

		inp = [ "-",			# PLOTGEOM
			"",			# Swap plotgeom
			self.plot["datafile"],	# Bin file
			"",			# Normalization
			"",			# Threshold density
			d,			# Binning
			limx, nx,		# Limits and bins
			limy, ny,
			limz, nz,
			"" ]			# End

		# Compare self.inp and inp
		try:
			if self.plot._inp == inp: return None
		except AttributeError:
			pass

		cmd = Project.command("gplevbin")
		self.log("Cmd: "+cmd+" <<EOF")
		for line in inp:
			self.log(line)
		self.log("EOF")
		rc, out = pipe.system(cmd, inp)
		for line in out: self.log(line)
		if rc: self.log("Return code %d"%(rc))

		try:
			# To be propagated into the viewport
			minmax = out[-2].split()
			self.plot["min"] = minmax[0]
			self.plot["max"] = minmax[1]
			n = -3
			while True:
				foo = out[n].split()
				if foo[0] == "Binning":
					self.plot["int"] = foo[-1]
					break
				n -= 1
		except:
			self.plot["min"] = ""
			self.plot["max"] = ""
			self.plot["int"] = ""

		try:
			os.stat("gplevh.dat")
		except:
			return "Error: File gplevh.dat is not created."

		if not Utils.isNewer("gplevh.dat",self.plot["datafile"]):
			return "Error: gplevh.dat is older than %s"%(self.plot["datafile"])

		fndat = self.plot.name+".dat"
		try:
			os.stat(fndat)
			try:
				os.remove(fndat)
			except OSError:
				return "Error: Cannot remove old %s file"%(fndat)
		except OSError:
			pass

		try:
			os.rename("gplevh.dat", fndat)
		except:
			return "Error: Cannot rename gplevh.dat to %s"%(fndat)

		if Project.cleanup:
			try: os.remove("gplevh.lim")
			except: pass
			try: os.remove("gplevh.npo")
			except: pass
			try: os.remove("gplevh.poi")
			except: pass
			try: os.remove("doslev.dat")
			except: pass

		self.plot._inp = inp

	# ----------------------------------------------------------------------
	# Run the usbmax program to process the data
	# ----------------------------------------------------------------------
	def _runUsbmax(self, det):
		# Find limits
		xmin = self.plot["bxmin"]
		xmax = self.plot["bxmax"]
		ymin = self.plot["bymin"]
		ymax = self.plot["bymax"]
		zmin = self.plot["bzmin"]
		zmax = self.plot["bzmax"]
		axis = self.plot.get("proj","Z")

		cmd = [ Project.command("usbmax"),
			"-d", str(det),
			"-a", axis,
			"-o", '%s.dat'%(self.plot.name)]
		if self.plot["hist"] == "1D Min":
			cmd.append("-m")
		cmd.append(self.plot["datafile"]) # Bin file

		if xmin!="": cmd.extend(["-x", str(xmin)])
		if xmax!="": cmd.extend(["-X", str(xmax)])
		if ymin!="": cmd.extend(["-y", str(ymin)])
		if ymax!="": cmd.extend(["-Y", str(ymax)])
		if zmin!="": cmd.extend(["-z", str(zmin)])
		if zmax!="": cmd.extend(["-Z", str(zmax)])

		try:
			if self.plot._inp == cmd: return None
		except AttributeError:
			pass

		self.log("\nCmd: "+" ".join(cmd))

		try:
			#out = subprocess.check_output(cmd)
			p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			out,err = p.communicate()
			self.log(out)
			if p.returncode != 0: return out
		#except subprocess.CalledProcessError, ex:	# as ex: (python3)
		except subprocess.CalledProcessError:
			self.log(str(sys.exc_info()[1]))
			return str(sys.exc_info()[1])
		self.plot._inp = cmd

	# ----------------------------------------------------------------------
	# Get plotgeom info and return True is it has changed
	# ----------------------------------------------------------------------
	def getPlotgeom(self, name):
		var  = self.plot.var
		det  = self.usr.detector[int(self.plot("det",1))-1]
		b10  = det.type % 10
		proj = ord(self.plot("proj","Z")) - ord("X")
		pos  = self.plot("geopos","")
		try:
			what = self.plot._what
		except AttributeError:
			self.plot._what = [""]*24
			what = self.plot._what

		# Get variables
		if   b10==1 or b10==7:
			# Cylindrical x=R, y=Phi, z=z
			# ymin=X, yhigh=Y
			# XXX bug in fluka does not keep the X,Y in the binning
			X = 0
			Y = 0
			if   proj==0: # X
				return False	# No geometry plot when Z-Phi
			elif proj==1: # Y
				# More tricky pos should be the rotation
				loP = bmath.Vector(X, Y, det.zlow)
				hiP = bmath.Vector(X+det.xhigh,  Y, det.zhigh)
				U = bmath.Vector.X
				V = bmath.Vector.Z
			elif proj==2: # Z
				if pos!="":
					z = float(pos)
				else:
					z = (det.zlow+det.zhigh)/2
				loP = bmath.Vector(-det.xhigh + X,
						   -det.xhigh + Y,
						    z)
				hiP = bmath.Vector( det.xhigh + X,
						    det.xhigh + Y,
						    z)
				U = bmath.Vector.X
				V = bmath.Vector.Y
		elif b10==2:
			return False
		elif b10==8:
			return False
		else:
			loP = bmath.Vector(det.xlow,  det.ylow,  det.zlow)
			hiP = bmath.Vector(det.xhigh, det.yhigh, det.zhigh)

			# Override limits if defined by user
			try: loP[0] = float(self.plot("xmin",""))
			except: pass
			try: loP[1] = float(self.plot("ymin",""))
			except: pass
			try: loP[2] = float(self.plot("zmin",""))
			except: pass

			try: hiP[0] = float(self.plot("xmax",""))
			except: pass
			try: hiP[1] = float(self.plot("ymax",""))
			except: pass
			try: hiP[2] = float(self.plot("zmax",""))
			except: pass

			# Unless there are limits set...
			if pos != "":
				loP[proj] = hiP[proj] = float(pos)
			else:
				loP[proj] = hiP[proj] = (loP[proj]+hiP[proj])/2

			if   proj == 0:	# X
				U = bmath.Vector.Y
				V = bmath.Vector.Z
			elif proj == 1: # Y
				U = bmath.Vector.X
				V = bmath.Vector.Z
			elif proj == 2: # Z
				U = bmath.Vector.X
				V = bmath.Vector.Y

		changed = self.project.inputModifiedTime() > self.plot.time
		if not changed: changed = what[ 7] != loP.x()
		if not changed: changed = what[ 8] != loP.y()
		if not changed: changed = what[ 9] != loP.z()
		if not changed: changed = what[10] != hiP.x()
		if not changed: changed = what[11] != hiP.y()
		if not changed: changed = what[12] != hiP.z()
		if not changed: changed = what[13] != V.x()
		if not changed: changed = what[14] != V.y()
		if not changed: changed = what[15] != V.z()
		if not changed: changed = what[16] != U.x()
		if not changed: changed = what[17] != U.y()
		if not changed: changed = what[18] != U.z()
		if not changed: changed = Utils.isNewer(self.project.input,
						name+".dat")

		if not changed: return False

		what[ 0] = "FORMAT"
		what[ 2] = 1.0	# Plot materials

		what[ 7] = loP.x()
		what[ 8] = loP.y()
		what[ 9] = loP.z()
		what[10] = hiP.x()
		what[11] = hiP.y()
		what[12] = hiP.z()

		what[13] = V.x()
		what[14] = V.y()
		what[15] = V.z()
		what[16] = U.x()
		what[17] = U.y()
		what[18] = U.z()

		what[19] = 1.0
		what[20] = 1.0
		what[21] = 2.0
		what[22] = 2.0

		return True

	# ----------------------------------------------------------------------
	# Run the gplevbin program to process the data
	# ----------------------------------------------------------------------
	def _mergeGeoData(self, usr, d, det):
		datfile = self.plot["datafile"]
		outfile = "%s.dat"%(self.plot.name)
		geofile = "%s.dat"%(self.plot["geo"])

		# check if conditions have changed
		try:
			what = self.plot._what
		except AttributeError:
			self.plot._what = [""]*24
			what = self.plot._what

		changed = datfile != what[0]
		if not changed: changed = geofile != what[1]
		if not changed: changed = Utils.isNewer(self.plot["datafile"], outfile)
		if not changed: changed = Utils.isNewer(geofile, outfile)
		if not changed: return None
		what[0] = datfile
		what[1] = geofile

		# If no geometry is provided then make histogram plot
		# of the data "Region" / "Value"
		# if a geometry is provided then superimpose the data on
		# the geometry
		try:
			data = usr.readData(d-1)
			stat = usr.readStat(d-1)
		except IOError:
			return "Error reading data file %s"%(d)

		fdata  = Data.unpackArray(data)
		edata  = Data.unpackArray(stat)
		del data
		del stat

		# Find regions
		cardlist   = self.project.input["REGION"]
#		cardlist.sort(Input.Card.cmpPos)
		cardlist.sort(key=lambda card: card.pos())
		regionList = [x.sdum() for x in cardlist
				if x.notIgnore() and x.sdum()!="&"]

		# open output file
		try:
			f = open(outfile,"w")
		except IOError:
			return "Error opening output filename"

		# Dump energy/region
		rstart  = int(det.xlow)-1	# zero based
		rend    = int(det.xhigh)-1
		rstep   = int(det.dx)
		rlen    = int(det.nx)

		lstart  = int(det.ylow)
		lstep   = int(det.dy)
		llen    = int(det.ny)

		special = det.type in (8,18)
		reglen  = len(regionList)
		if special:
			# Special binning
			# accept only if the last field has 1 bin!
			if det.nz != 1:
				return "Error cannot process special binning with more than one binning in the third variable"

			idx = 0
			for j in range(llen):
				l  = j*lstep + lstart
				jl = j*rlen
				for i in range(rlen):
					r = int(i*rstep) + rstart
					if r<reglen:
						f.write("%d %d %s:%d %g %g\n" % \
							(r+1, l, regionList[r],l,
							 fdata[idx], edata[idx]))
					idx += 1

		else:	# Region binning
			for i in range(rlen):
				r = int(i*rstep) + rstart
				if r<reglen:
					f.write("%d 0 %s:0 %g %g\n" % \
						(r+1, regionList[r],
						 fdata[i], edata[i]))

		# If custom geometry plot is available, get the second part
		# and add the energy/region
		if self.plot["geo"] not in (GEO_NO, GEO_AUTO):
			try:
				fgeo = open(geofile,"r")

				# Make a new separation
				f.write("\n\n")

				# Scan fgeo for the scan part (double empty line)
				count = 0
				for line in fgeo:
					if line=="\n":
						count += 1
						if count == 2: break
					else:
						count = 0

				# Read remaining lines and append energy deposition if any
				count = 0
				for line in fgeo:
					if line=="\n":
						f.write(line)
						count += 1
						if count == 2: break
						continue
					count = 0
					words = line.split()

					# Check if it is a round integer
					r,m = divmod(int(words[5])-1-rstart, rstep)
					l = int(words[7])	# before we substitute
					# set the default 0.0
					words[6] = '0.0'
					words[7] = '0.0'
					if m==0 and 0 <= r < rlen:
						if special:
							l,m = divmod(l-lstart, lstep)
							if m==0 and 0 <= l < llen:
								idx = l*rlen+r
								words[6] = str(fdata[idx])
								words[7] = str(edata[idx])
						else:
							words[6] = str(fdata[r])
							words[7] = str(edata[r])
					f.write("%s\n"%(" ".join(words[:8])))

				# Copy the remaining lines
				for line in fgeo: f.write(line)

				fgeo.close()
			except IOError:
				return sys.exc_info()[1]

			f.close()

		del fdata
		del edata

	# ----------------------------------------------------------------------
	# Prepare the data for the resnuclei plot
	# ----------------------------------------------------------------------
	def _resnucleiPlot(self):
		# Check if we have to re-write the file
		d = int(self.plot("det",1))-1
		if d<0: return "Error: RESNUCLEi invalid detector specified"

		usr = Data.Resnuclei()
		try:
			usr.readHeader(self.plot["datafile"])
		except IOError:
			return "Error: reading resnuclei file %s"%(sys.exc_info()[1])

		if len(usr.detector)==0: return "Error: empty RESNUCLEi file"
		current = [self.plot.name, self.plot["datafile"], d]
		try:
			if current == self.plot._previous:
				return None
		except AttributeError:
			pass

		det = usr.detector[d]
		try:
			data = usr.readData(d)
			stat = usr.readStat(d)
		except IOError:
			return "Error reading data file %s"%(d)

		fdata = Data.unpackArray(data)
		del data

		if stat:
			total = Data.unpackArray(stat[0])
			A     = Data.unpackArray(stat[1])
			eA    = Data.unpackArray(stat[2])
			Z     = Data.unpackArray(stat[3])
			eZ    = Data.unpackArray(stat[4])
			edata = Data.unpackArray(stat[5])
			del stat
			stat  = True
		else:
			stat  = False

		# Fill min/max values
		try:
			fmin = min([x for x in fdata if x>0.0])
		except:
			return "Error: Empty data set"

		fmax = max(fdata)
		#self.usr_min.config(text="%-10g"%(fmin))
		#self.usr_max.config(text="%-10g"%(fmax))

		# open output file
		try:
			f = open(self.plot.name+".dat","w")
		except IOError:
			return "Error opening output filename"

		zh  = det.zhigh
		mh  = det.mhigh
		nmz = det.nmzmin

		length = len(fdata)

		# Z distribution
		f.write("# Detector n: 1 Z\n")
		if stat:
			for z in range(zh):
				f.write("%d %d %15g %g\n"%(z+1,z+2,Z[z],eZ[z]*100.0))
		else:
			for z in range(zh):
				s  = 0.0
				for m in range(det.mhigh):
					pos = z + m * zh
					if 0 <= pos < length:
						s += fdata[pos]
				f.write("%d %d %15g 0.\n"%(z+1,z+2,s))

		# A distribution
		f.write("\n\n# Detector n: 2 A\n")
		amax = 2*zh + mh + nmz
		if stat:
			for a in range(1,min(amax,len(A))+1):
				f.write("%d %d %15g %g\n"%(a,a+1, A[a-1],eA[a-1]))
		else:
			for a in range(1,amax+1):
				s = 0.0
				for z in range(zh):
					m = a - 2*z - nmz - 3
					pos = z + m * zh
					if 0 <= pos < length:
						s += fdata[pos]
				f.write("%d %d %15g 0.\n"%(a,a+1, s))

		# Z,M distribution
		#f.write("\n\n#Isotope Yield as a function of Atomic " \
		#	"Number / N-Z-NMZmin (nuclei / cmc / pr)\n")
		#pos = 0
		#for m in range(mh):
		#	for z in range(zh):
		#		f.write("%d %d %15g\n"%(z+1,m+1,fdata[pos]))
		#		pos += 1
		#	f.write("\n")

		# Z,A distribution
		f.write("\n\n# Detector n: 3 Z/A\n")
		for a in range(1,amax+1):
			f.write("# Block n: 3 A:%d\n"%(a))
			for z in range(zh):
				m = a - 2*z - nmz - 3
				if m<0 or m>=mh:
					f.write("%d %d 0. 0.\n"%(z+1,a))
				else:
					pos = z + m * zh
					if stat:
						err = edata[pos]*100
					else:
						err = 0.0
					f.write("%d %d %15g %g\n"%(z+1,a,
						fdata[pos], err))
			f.write("\n")

		f.close()
		del fdata
		if stat: del edata
		self.plot._previous = current

	# ----------------------------------------------------------------------
	# Userdump file creation
	# ----------------------------------------------------------------------
	def _userdumpPlot(self):
		try:
			case = USERDUMP_CASE.index(self.plot["case"])
		except:
			case = 0

		mgfile = self.plot.name + ".dat"
		#if self.case == case and \
		#   mgfile == self.mgfile and \
		#   Utils.isNewer(self.mgfile, self.plot["datafile"]):
		#	return True
		try:
			fout = open(os.path.join(self.project.dir,mgfile),"w")
		except IOError:
			return "Error: cannot open file %s for writing"%(mgfile)

		mg = Data.Mgdraw(self.plot["datafile"])
		if mg.hnd is None: return "Error: opening userdump file %s"%(self.plot["datafile"])

		try: from_ = int(self.plot("from",0))
		except: from_ = 0
		try: n     = int(self.plot("n",0))
		except: n = 0
		try: emin = float(self.plot["emin"])
		except: emin = None
		try: emax = float(self.plot["emax"])
		except: emax = None
		particles = []
		for p in self.plot["particles"].split():
			particles.append(Input.Particle._db[p].id)
		if len(particles)==0:
			particles = None

		if n==0:
			to_ = from_ + 10000
		else:
			to_ = from_ + n

		if case == 2:
			# Find limits
			dmin = [ 1e30]*11
			dmax = [-1e30]*11

		i = 0
		startTime = time.time()
		last = startTime
		try:
			while True:
				try:
					t = mg.readEvent(case)
				except:
					return "ERROR not a valid MGREAD file"

				if t is None: break
				elif t != case: continue

				if case == 0:		# Tracking
					n      = mg.ntrack
					ptype  = mg.jtrack
					weight = mg.wtrack
					try:
						mass = Input.Particle._db[int(ptype)].mass
					except:
						# FIXME XXX XXX
						# for ptype < -6!!!!!
						# or ptype == -2
						# only kinetic energy is possible
						mass = 0.0

					ekin   = mg.etrack - mass
					if emin and ekin < emin: continue
					if emax and ekin > emax: continue
					if particles and ptype not in particles:
						continue

					i += 1
					if i < from_: continue
					elif i > to_: break

					j = 0
					while j < n:
						x = mg.data[j];	j += 1
						y = mg.data[j];	j += 1
						z = mg.data[j];	j += 1
						fout.write("%d %g %g %g %g %g\n" % \
							(ptype, ekin, weight, x, y, z))
					fout.write("\n")

				elif case == 2:		# Source particles
					(ptype, etot, weight, x, y, z, tx, ty, tz) \
							= mg.data
					try:
						mass = Input.Particle._db[int(ptype)].mass
					except KeyError:
						return "Error: Unknown particle id=%d"%(int(ptype))
					ekin = etot - mass
					if ekin < 0.0: ekin = 0.0
					if emin and ekin < emin: continue
					if emax and ekin > emax: continue
					if particles and ptype not in particles:
						continue

					i += 1
					if   i < from_: continue
					elif i > to_: break

					pmom = math.sqrt(ekin*(etot+mass))
					fout.write("%d %g %g %g %g %g %g %g %g %g %g\n" % \
						(ptype, etot, weight, \
							x, y, z, tx, ty, tz, ekin, pmom))

					for d in range(9):
						dmin[d] = min(dmin[d], mg.data[d])
						dmax[d] = max(dmax[d], mg.data[d])
					dmin[ 9] = min(dmin[ 9], ekin)
					dmin[10] = min(dmin[10], pmom)
					dmax[ 9] = max(dmin[ 9], ekin)
					dmax[10] = max(dmin[10], pmom)

				#if time.time() - last > 2.0:
				#	last = time.time()
					#self.processed.config(text="%d"%(i))
					#self.processed.update_idletasks()

		except KeyboardInterrupt:
			self.log("*** Interrupted after %s s ***")
		else:
			self.log("Processing finished in %g s"%(time.time()-startTime))
#		self.processed.config(text="%d"%(i))
		mg.close()
		fout.close()

		if case == 2:
			#self.limfile = "%s-%d.lim" % (self.plot.name, case)
			self.limfile = self.plot.name + ".lim"
			try:
				fout = open(os.path.join(self.project.dir, self.limfile),"w")
			except IOError:
				return "Error: cannot open file %s for writing"%(self.limfile)
			for d in range(len(dmin)):
				fout.write("%g %g\n"%(dmin[d], dmax[d]))
			fout.close()
			self.dmin = dmin
			self.dmax = dmax
		else:
			self.limfile   = ""
			self.dmin = None
			self.dmax = None
