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

# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	31-Mar-2014

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import sys
import bmath
import io
from math import *

#===============================================================================
# Histogram 1D
#===============================================================================
class Histogram1D:
	"""1D histogram with variable or equidistance binning"""
	UNDER = -1
	OVER  = -2

	# ----------------------------------------------------------------------
	def __init__(self, bins=100, xlow=0.0, xhigh=1.0):
		"""
		Initialize a histogram with bins spanning from xlow to xhigh
		bins can the number of bins, or
		be a list of bins or a file for a variable binning histogram
		"""
		self.under   = 0.0
		self.over    = 0.0
		self.total   = 0.0
		self.entries = 0

		if isinstance(bins,str) or isinstance(bins,file):
			self.load(bins)
		elif isinstance(bins,list):
			self.nbins  = len(bins)-1	# number of bins
			self.xbins  = bins		# x bins list
			self._xlow  = bins[0]		# limits
			self._xhigh = bins[-1]
			self.xstep  = None		# Step is variable
			self.h      = [0.0]*self.nbins
			self.eh     = [0.0]*self.nbins
		else:
			self.nbins  = bins
			self.xbins  = None
			self._xlow  = xlow
			self._xhigh = xhigh
			self.xstep  = (self._xhigh-xlow)/float(bins)
			self.h      = [0.0]*self.nbins
			self.eh     = [0.0]*self.nbins

	# ----------------------------------------------------------------------
	def logBinning(self, bins, xlow, xhigh):
		"""Create equal log10 binnings"""
		self.nbins  = bins
		self._xlow  = xlow
		self._xhigh = xhigh
		self.xstep  = None
		self.h      = [0.0]*self.nbins
		self.eh     = [0.0]*self.nbins

		self.xbins  = []
		x = log10(xlow)
		s = (log10(xhigh)-x)/float(bins)
		for i in range(bins):
			self.xbins.append(10.0**x)
			if abs(x)<1e-15: x = 0.0
			x += s
		self.xbins.append(10.0**x)

	# ----------------------------------------------------------------------
	def clear(self):
		"""Delete data"""
		for i in range(bins):
			self.h[i] = self.eh[i] = 0.0
		self.under   = 0.0
		self.over    = 0.0
		self.total   = 0.0
		self.entries = 0

	# ----------------------------------------------------------------------
	def load(self, fin):
		"""Load a histogram from file"""
		# FIXME variable binning
		if isinstance(fin,str):
			fin = open(fin,"r")
			close = True
		else:
			close = False

		self.h     = []
		self.eh    = []
		self.xbins = []
		self.nbins = 0
		first      = True
		for line in fin:
			if line[0]=="#":
				word = line[1:].lower().split()
				if word[0] in ("bins","nbins"):
					self.nbins = int(word[1])
				elif word[0] in ("xmin", "xlow"):
					self._xlow = float(word[1])
				elif word[0] in ("xmax", "xhigh"):
					self._xhigh = float(word[1])
			else:
				word = line.split()
				if len(word)==1:
					self.h.append(float(word[0]))
					self.eh.append(0.0)
				elif len(word)==2:
					self.h.append(float(word[1]))
					self.eh.append(0.0)
				elif len(word)==3:
					self.h.append(float(word[1]))
					self.eh.append(float(word[2]))
				elif len(word)==4:
					if self.nbins==0:
						if first:
							self.xbins.append(float(word[0]))
							first = False
						self.xbins.append(float(word[1]))
					self.h.append(float(word[2]))
					self.eh.append(float(word[3]))

		if self.xbins:
			self._xlow  = self.xbins[0]
			self._xhigh = self.xbins[-1]
			self.nbins  = len(self.xbins)-1
			self.xstep  = None
		else:
			self.xbins = None
			self.xstep = (self._xhigh-self._xlow)/float(self.nbins)

		if close:
			fin.close()
	read = load

	# ----------------------------------------------------------------------
	# Save histogram into a file
	# ----------------------------------------------------------------------
	def save(self, fout=sys.stdout, cols=None):
		"""Save histogram to file"""
		# FIXME variable binning
		if isinstance(fout,str):
			fout = open(fout,"w")
			close = True
		else:
			close = False
		fout.write("# xmin    %.10g\n"%(self._xlow))
		fout.write("# xmax    %.10g\n"%(self._xhigh))
		if self.xstep is not None:
			fout.write("# nbins   %d\n"   %(self.nbins))
			fout.write("# xstep   %.10g\n"%(self.xstep))
		fout.write("# entries %d\n"   %(self.entries))
		fout.write("# weight  %.10g\n"%(self.total))
		fout.write("# under   %.10g\n"%(self.under))
		fout.write("# over    %.10g\n"%(self.over))
		x = self._xlow

		if self.xstep is None:
			for i in range(self.nbins):
				x  = self.xbins[i]
				xh = self.xbins[i+1]
				y  = self.h[i]
				e  = self.eh[i]
				fout.write("%15.10g %15.10g %15.10g %15.10g\n"%(x,xh,y,e))
		else:
			for y,e in zip(self.h,self.eh):
				xh = x + self.xstep
				if abs(xh)<1e-14: xh = 0.0
				if cols==2:
					fout.write("%15.10g %15.10g %15.10g\n"%(x,y))
				elif cols==3:
					fout.write("%15.10g %15.10g %15.10g\n"%(x,y,e))
				else:
					fout.write("%15.10g %15.10g %15.10g %15.10g\n"%(x,xh,y,e))
				x = xh

		if close:
			fout.close()
	write = save

	# ----------------------------------------------------------------------
	def dump(self):
		"""Dump histogram to stdout"""
		self.save()

	# ----------------------------------------------------------------------
	def dumpString(self):
		"""Dump histogram to StringIO"""
		sio = io.StringIO()
		self.save(sio)
		return sio

	# ----------------------------------------------------------------------
	def copy(self, hist):
		"""Copy histogram from hist"""
		self.nbins  = hist.nbins
		if hist.xbins is None:
			self.xbins = None
		else:
			self.xbins = hist.xbins[:]
		self._xlow   = hist._xlow
		self._xhigh  = hist._xhigh
		self.xstep   = hist.xstep
		self.h       = hist.h[:]
		self.eh      = hist.eh[:]
		self.under   = hist.under
		self.over    = hist.over
		self.entries = hist.entries
		self.total   = hist.total

	# ----------------------------------------------------------------------
	def clone(self):
		"""Return a clone of the current histogram"""
		hist = Histogram()
		hist.copy(self)
		return hist

	# ----------------------------------------------------------------------
	def empty(self):
		"""Zero histogram"""
		self.h       = [0.0] * self.nbins
		self.eh      = [0.0] * self.nbins
		self.under   = 0.0
		self.over    = 0.0
		self.entries = 0
		self.total   = 0.0

	# ----------------------------------------------------------------------
	def __getitem__(self, i):
		"""return ith content of histogram"""
		return self.h[i]

	def __setitem__(self, i, y):
		"""set ith content of histogram"""
		self.h[i] = y

	# ----------------------------------------------------------------------
	def bin(self, x):
		"""Return bin corresponding to value x,
		return UNDER on underflow and OVER on overflow
		"""
		if x<self._xlow:
			return Histogram.UNDER

		elif x>=self._xhigh:
			return Histogram.OVER

		if self.xstep is None:
			# Perform a binary search to find the correct interval
			low  = 0
			high = self.nbins-1
			while True:
				mid = (low+high)//2
				if mid==low: return mid
				if x > self.xbins[mid]:
					low = mid
				elif x < self.xbins[mid]:
					high = mid
				else:
					return mid
		else:
			return int((x - self._xlow) // self.xstep)

	# ----------------------------------------------------------------------
	def x(self, i=None):
		"""Return low value of ith bin"""
		if self.xstep is None:
			if i is None:
				return self.xbins
			else:
				return self.xbins[i]
		else:
			if i is None:
				return [self.xstep*i + self._xlow for i in range(self.nbins)]
			else:
				return self.xstep*i + self._xlow
	getx = x

	# ----------------------------------------------------------------------
	def xlow(self, i):
		"""return low value of ith bin"""
		if self.xstep is None:
			return self.xbins[i]
		else:
			return self.xstep*i + self._xlow

	# ----------------------------------------------------------------------
	def xhigh(self, i):
		"""return high value of ith bin"""
		if self.xstep is None:
			return self.xbins[i+1]
		else:
			return self.xstep*(i+1) + self._xlow

	# ----------------------------------------------------------------------
	def xmid(self, i):
		"""Return x-mean of bin i"""
		return 0.5*(self.xlow(i)+self.xhigh(i))

	# ----------------------------------------------------------------------
	def xgeo(self, i):
		"""Return x-geometric mean of bin i"""
		return sqrt(self.xlow(i)*self.xhigh(i))

	# ----------------------------------------------------------------------
	def dx(self, i):
		"""return size of ith bin"""
		if self.xstep is None:
			return self.xbins[i+1]-self.xbins[i]
		else:
			return self.xstep

	# ----------------------------------------------------------------------
	def error(self, i, e=None):
		"""return or set the error value for the ith bin"""
		if e is None:
			return self.eh[i]
		else:
			self.eh[i] = e

	# ----------------------------------------------------------------------
	def convertError2Relative(self):
		"""Convert error to relative in percent"""
		for i,y in enumerate(self.h):
			self.eh[i] *= 100.0/y

	# ----------------------------------------------------------------------
	def convertError2Absolute(self):
		"""Convert errot to absolute value"""
		for i,y in enumerate(self.h):
			self.eh[i] *= y/100.0

	# ----------------------------------------------------------------------
	def fill(self, x, w=1.0):
		"""Fill/add to position x weight w"""
		self.entries += 1
		self.total += w
		i = self.bin(x)
		if i==Histogram.UNDER:
			self.under += w

		elif i==Histogram.OVER:
			self.over  += w

		else:
			self.h[i] += w

	# ----------------------------------------------------------------------
	# Normalize histogram with a factor f
	# ----------------------------------------------------------------------
	def norm(self, f=None):
		"""Normalize histogram. if f is None divide by the bin width,
		   else multiply with the factor provided"""
		if f is None:
			for i in range(self.nbins):
				dx = self.dx(i)
				self.h[i]  /= dx
				self.eh[i] /= dx
		else:
			for i in range(self.nbins):
				self.h[i]  *= f
				self.eh[i] *= f
	normalize = norm

	# ----------------------------------------------------------------------
	def isolethargic(self):
		# FIXME variable binning
		"""Convert a log10(histogram) to isolethargic"""
		f = 10**(self.xstep/2) / (10**self.xstep-1.0)
		self.norm(f)

	# ----------------------------------------------------------------------
	def isSame(self, hist):
		# FIXME variable binning
		"""Return true if histograms have the same limits"""
		if self.nbins != hist.nbins: return False
		if self._xlow  != hist._xlow:  return False
		if self._xhigh != hist._xhigh: return False
		return True

	# ----------------------------------------------------------------------
	def __iadd__(self, a):
		# FIXME variable binning
		if isinstance(a,float):
			for i in range(self.nbins):
				self.h[i]  += a

		elif isinstance(a,Histogram):
			if not a.isSame(self): raise
			for i in range(self.nbins):
				self.eh[i] = sqrt(self.eh[i]**2 + a.eh[i]**2)
				self.h[i] += a[i]
		return self

	# ----------------------------------------------------------------------
	def __isub__(self, a):
		# FIXME variable binning
		if isinstance(a,float):
			for i in range(self.nbins):
				self.h[i]  -= a

		elif isinstance(a,Histogram):
			if not a.isSame(self): raise
			for i in range(self.nbins):
				self.eh[i] = sqrt(self.eh[i]**2 + a.eh[i]**2)
				self.h[i] -= a[i]
		return self

	# ----------------------------------------------------------------------
	def __imul__(self, a):
		# FIXME variable binning
		if isinstance(a,float):
			for i in range(self.nbins):
				self.h[i]  *= a
				self.eh[i] *= a

		elif isinstance(a,Histogram):
			if not a.isSame(self): raise
			for i in range(self.nbins):
				self.eh[i] = sqrt((a[i]*self.eh[i])**2 + \
						(self[i]*a.eh[i])**2)
				self.h[i] *= a[i]
		return self

	# ----------------------------------------------------------------------
	def __idiv__(self, a):
		# FIXME variable binning
		if isinstance(a,float):
			for i in range(self.nbins):
				self.h[i]  /= a
				self.eh[i] /= a

		elif isinstance(a,Histogram):
			if not a.isSame(self): raise
			for i in range(self.nbins):
				try:
					self.eh[i] = sqrt((self.eh[i]/a[i])**2 + \
						     (self[i]/a[i]**2*a.eh[i])**2)
				except:
					self.eh[i] = 0.0
				self.h[i] /= a[i]
		return self

	# ----------------------------------------------------------------------
	def __radd__(self, a):
		# FIXME variable binning
		hist = a.clone()
		hist += self
		return hist

	# ----------------------------------------------------------------------
	def __rsub__(self, a):
		# FIXME variable binning
		hist = a.clone()
		hist -= self
		return hist

	# ----------------------------------------------------------------------
	def __rmul__(self, a):
		# FIXME variable binning
		hist = a.clone()
		hist *= self
		return hist

	# ----------------------------------------------------------------------
	def __rdiv__(self, a):
		# FIXME variable binning
		hist = a.clone()
		hist /= self
		return hist

	# ----------------------------------------------------------------------
	# Convert histogram to a cumulative function
	# @return total sum
	# ----------------------------------------------------------------------
	def cumulative(self):
		"""Convert to a running sum(0,n,h(i))"""
		S = 0.0
		for i,y in enumerate(self.h):
			S += y
			self.h[i] = S
		return S

	# ----------------------------------------------------------------------
	# Convert to a running integral
	# @return total integral
	# ----------------------------------------------------------------------
	def integrate(self):
		"""Convert to a running integral int(0,n,h(i)*dx)"""
		S = 0.0
		for i,y in enumerate(self.h):
			S += y*self.dx(i)
			self.h[i] = S
		return S

	# ----------------------------------------------------------------------
	# Calculate the derivative of the histogram
	# ----------------------------------------------------------------------
	def derivative(self):
		for i in range(self.nbins-1,1,-1):
			self.h[i] = (self.h[i]-self.h[i-1]) / self.dx(i)
		self.h[0] = 0.0

	# -----------------------------------------------------------------------------
	def mean(self):
		"""return mean and rms value"""
		# FIXME variable binning
		sy = 0.0
		sxy = 0.0
		sx2y = 0.0

		x = self.xmid(0)
		for y in self.h:
			sy   += y
			sxy  += x*y
			sx2y += x*x*y
			x += self.xstep

		if sy == 0.0: return None,None

		mx  = sxy/sy
		var = sx2y/sy - mx**2
		if var<0.0: var = 0.0
		return mx, sqrt(var)

	# -----------------------------------------------------------------------------
	def rms(self):
		# FIXME variable binning
		return self.mean()[1]

	# -----------------------------------------------------------------------------
	# Draw histogram using plot engine
	# -----------------------------------------------------------------------------
	def draw(self, engine, options=None):
		"""plot histogram using engine and options"""
		if options is None:
			options = "using 1:3 w steps not",
		engine.plot(options, [(self.xlow(i), self.xhigh(i), self.h[i], self.eh[i])
						for i in range(self.nbins)])

	# -----------------------------------------------------------------------------
	@staticmethod
	def average(histograms):
		"""create a new histogram containing the mean value of all histograms
		ALL HISTOGRAMS MUST WITH THE SAME PARAMETERS
		and for the same number of events
		"""
		#FIXME variable binning
		h = histograms[0]
		new = Histogram(h.nbins, h._xlow, h._xhigh)

		n = float(len(histograms))
		sn = sqrt(n-1.0)

		for i in range(h.nbins):
			t = t2 = 0.0
			for h in histograms:
				t  += h[i]
				t2 += h[i]**2

			new[i] = t/n
			var = t2/n - new[i]**2
			if var<0.0: var = 0.0
			new.error(i, sqrt(var) / sn)

		# Under, Over
		u = 0.0
		o = 0.0
		for h in histograms:
			u  += h.under
			o  += h.over

		new.under = u/n
		new.over  = o/n

		return new

# aliases
Histogram = Histogram1D
H1D = Histogram1D
H1  = Histogram1D

#===============================================================================
# Histogram 2D
#===============================================================================
class Histogram2D:
	"""2D histogram with equi distance bins"""

	# ----------------------------------------------------------------------
	def __init__(self, xbins=10, xlow=0.0, xhigh=0.1, ybins=10, ylow=0.0, yhigh=1.0):
		"""Initialize a 2D histogram with xbins from xlow to xhigh,
		and ybins from ylow to yhigh
		"""
		#self.xunder = 0.0
		#self.yunder = 0.0
		#self.under  = [0.0]*xbins
		#self.over   = [0.0]*xbins
		self.total   = 0.0
		self.entries = 0
		if isinstance(xbins,str) or isinstance(xbins,file):
			self.load(xbins)
		else:
			self.xbins  = xbins
			self.xlow   = xlow
			self.xhigh  = xhigh
			self.xstep  = (xhigh-xlow)/float(xbins)

			self.ybins  = ybins
			self.ylow   = ylow
			self.yhigh  = yhigh
			self.ystep  = (yhigh-ylow)/float(ybins)

			# List of lists
			self.h = []
			self.eh = []
			for i in range(xbins):
				self.h.append([0.0]*ybins)
				self.eh.append([0.0]*ybins)

	# ----------------------------------------------------------------------
	def fill(self, x, y, w=1.0):
		"""Fill/add to position x,y weight w"""
		self.entries += 1
		self.total += w
		if x>=self.xlow and x<self.xhigh and \
		   y>=self.ylow and y<self.yhigh:
			i = int((x - self.xlow) // self.xstep)
			j = int((y - self.ylow) // self.ystep)
			self.h[i][j] += w
#		if x<self.xlow:
#			self.under += w
#		elif x>self.xhigh:
#			self.over  += w
#		else:

#	# ----------------------------------------------------------------------
#	def copy(self, hist):
#		self.xbins  = hist.xbins
#		self.xlow   = hist.xlow
#		self.xhigh  = hist.xhigh
#		self.h      = hist.h[:]
#		self.eh     = hist.eh[:]
#		self.under  = hist.under
#		self.over   = hist.over
#		self.xstep  = hist.xstep
#		self.entries= hist.entries
#		self.total  = hist.total

#	# ----------------------------------------------------------------------
#	def clone(self):
#		hist = Histogram(self.xbins, self.xlow, self.xhigh)
#		hist.copy(self)
#		return hist

	# ----------------------------------------------------------------------
#	def __getitem__(self, i): return self.h[i]
#	def __setitem__(self, i, y): self.h[i] = y

	# ----------------------------------------------------------------------
	def x(self, i=None):
		"""return x value of ith column"""
		if i is None:
			return [self.xstep*i + self.xlow for i in range(self.xbins)]
		else:
			return self.xstep*i + self.xlow
	getx = x

	# ----------------------------------------------------------------------
	def y(self, j=None):
		"""return y value of jth row"""
		if j is None:
			return [self.ystep*j + self.ylow for j in range(self.ybins)]
		else:
			return self.ystep*j + self.ylow
	gety = y

	# ----------------------------------------------------------------------
	def xmid(self, i):
		"""return middle x value of a column i"""
		return self.xstep*(i+0.5) + self.xlow
	def ymid(self, i):
		"""return middle y value of a row i"""
		return self.ystep*(i+0.5) + self.ylow

	# ----------------------------------------------------------------------
	def error(self, i, j, e=None):
		"""return error of [i,j]] or set if e is not None"""
		if e is None:
			return self.eh[i][j]
		else:
			self.eh[i][j] = e

	# ----------------------------------------------------------------------
	def norm(self, f=None):
		"""normalize histogram with a factor f"""
		if f is None: f = 1.0 / (self.xstep * self.ystep)
		for i in range(self.xbins):
			for j in range(self.ybins):
				self.h[i][j]  *= f
				self.eh[i][j] *= f
	normalize = norm

	# ----------------------------------------------------------------------
	def xslice(self, ifrom, ito=None):
		"""return a x-slice (sum) (Y-histogram) from [ifrom : ito)"""
		hsum = Histogram1D(self.ybins, self.ylow, self.yhigh)
		h1   = Histogram1D(self.ybins, self.ylow, self.yhigh)
		if ito is None: ito = ifrom + 1
		for i in range(ifrom, ito):
			for j,(v,e) in enumerate(zip(self.h[i], self.eh[i])):
				h1[j] = v
				h1.error(j, e)
			hsum += h1	# WARNING errors will be WRONG!!!!!!!
		return hsum

	# ----------------------------------------------------------------------
	def rebiny(self, n):
		"""rebin y every n bins"""
		for i,yh in enumerate(self.h):
			row = []
			s   = 0.0
			c   = 0.0
			for val in yh:
				s += val
				c += 1.0
				if c == n:
					row.append(s/c)
					s = 0.0
					c = 0.0
			if c > 0:
				# FIXME should I change the limit of yhigh?
				row.append(s/c)

			self.h[i]  = row
		self.ybins = len(row)

#	# ----------------------------------------------------------------------
#	def isSame(self, hist):
#		"""Return true if histograms have the same limits"""
#		if self.xbins != hist.xbins: return False
#		if self.xlow  != hist.xlow:  return False
#		if self.xhigh != hist.xhigh: return False
#		return True

#	# ----------------------------------------------------------------------
#	def __iadd__(self, a):
#		if isinstance(a,float):
#			for i in range(self.xbins):
#				self.h[i]  += a
#
#		elif isinstance(a,Histogram):
#			if not a.isSame(self): raise
#			for i in range(self.xbins):
#				self.eh[i] = sqrt(self.eh[i]**2 + a.eh[i]**2)
#				self.h[i] += a[i]
#		return self

#	# ----------------------------------------------------------------------
#	def __isub__(self, a):
#		if isinstance(a,float):
#			for i in range(self.xbins):
#				self.h[i]  -= a
#
#		elif isinstance(a,Histogram):
#			if not a.isSame(self): raise
#			for i in range(self.xbins):
#				self.eh[i] = sqrt(self.eh[i]**2 + a.eh[i]**2)
#				self.h[i] -= a[i]
#		return self

#	# ----------------------------------------------------------------------
#	def __imul__(self, a):
#		if isinstance(a,float):
#			for i in range(self.xbins):
#				self.h[i]  *= a
#				self.eh[i] *= a
#
#		elif isinstance(a,Histogram):
#			if not a.isSame(self): raise
#			for i in range(self.xbins):
#				self.eh[i] = sqrt((a[i]*self.eh[i])**2 + (self[i]*a.eh[i])**2)
#				self.h[i] *= a[i]
#		return self

#	# ----------------------------------------------------------------------
#	def __idiv__(self, a):
#		if isinstance(a,float):
#			for i in range(self.xbins):
#				self.h[i]  /= a
#				self.eh[i] /= a
#
#		elif isinstance(a,Histogram):
#			if not a.isSame(self): raise
#			for i in range(self.xbins):
#				try:
#					self.eh[i] = sqrt((self.eh[i]/a[i])**2 + (self[i]/a[i]**2*a.eh[i])**2)
#				except:
#					self.eh[i] = 0.0
#				self.h[i] /= a[i]
#		return self
#
#	# ----------------------------------------------------------------------
#	def __radd__(self, a):
#		hist = a.clone()
#		hist += self
#		return hist
#
#	# ----------------------------------------------------------------------
#	def __rsub__(self, a):
#		hist = a.clone()
#		hist -= self
#		return hist
#
#	# ----------------------------------------------------------------------
#	def __rmul__(self, a):
#		hist = a.clone()
#		hist *= self
#		return hist
#
#	# ----------------------------------------------------------------------
#	def __rdiv__(self, a):
#		hist = a.clone()
#		hist /= self
#		return hist

	# ----------------------------------------------------------------------
	def save(self, fout=sys.stdout):
		"""save to file fout"""
		if isinstance(fout,str):
			fout = open(fout,"w")
			close = True
		else:
			close = False
		fout.write("# xbins   %d\n"   %(self.xbins))
		fout.write("# xmin    %.10g\n"%(self.xlow))
		fout.write("# xmax    %.10g\n"%(self.xhigh))
		fout.write("# xstep   %.10g\n"%(self.xstep))
		fout.write("# ybins   %d\n"   %(self.ybins))
		fout.write("# ymin    %.10g\n"%(self.ylow))
		fout.write("# ymax    %.10g\n"%(self.yhigh))
		fout.write("# ystep   %.10g\n"%(self.ystep))
#		fout.write("# entries %d\n"   %(self.entries))
		fout.write("# weight  %.10g\n"%(self.total))
#		fout.write("# under   %.10g\n"%(self.under))
#		fout.write("# over    %.10g\n"%(self.over))
		x = self.xlow
		for i in range(self.xbins):
			y = self.ylow
			for v,e in zip(self.h[i],self.eh[i]):
				fout.write("%.10g %.10g %.10g %.10g\n"%(x,y,v,e))
				y += self.ystep
				if abs(y)<1e-15: y = 0.0
			fout.write("\n")
			x += self.xstep
			if abs(x)<1e-15: x = 0.0

		if close:
			fout.close()
	write = save

	# ----------------------------------------------------------------------
	def load(self, fin):
		"""load histogram from file"""
		if isinstance(fin,str):
			fin = open(fin,"r")
			close = True
		else:
			close = False

		self.h  = []
		self.eh = []
		row = []
		erow = []

		for line in fin:
			line = line.strip()
			if not line:
				if row:
					self.h.append(row)
					self.eh.append(erow)
				row = []
				erow = []
			elif line[0]=="#":
				word = line[1:].lower().split()
				if word[0] in ("bins", "xbins"):
					self.xbins = int(word[1])
				elif word[0] in ("xmin", "xlow"):
					self.xlow = float(word[1])
				elif word[0] in ("xmax", "xhigh"):
					self.xhigh = float(word[1])
				elif word[0] == "ybins":
					self.ybins = int(word[1])
				elif word[0] in ("ymin", "ylow"):
					self.ylow = float(word[1])
				elif word[0] in ("ymax", "yhigh"):
					self.yhigh = float(word[1])
			else:
				word = line.split()
				if len(word) in (1,3):
					row.append(float(word[-1]))
					erow.append(0.0)
				else:
					row.append(float(word[-2]))
					erow.append(float(word[-1]))

		if row:
			self.h.append(row)
			self.eh.append(erow)

		self.xstep  = (self.xhigh-self.xlow)/float(self.xbins)
		self.ystep  = (self.yhigh-self.ylow)/float(self.ybins)

		if close:
			fin.close()
	read = load

	# -----------------------------------------------------------------------------
	@staticmethod
	def average(histograms):
		"""
		create a new histogram containing the mean value of all histograms
		ALL HISTOGRAMS MUST WITH THE SAME PARAMETERS
		and for the same number of events
		"""
		h = histograms[0]
		new = Histogram2D(h.xbins, h.xlow, h.xhigh, h.ybins, h.ylow, h.yhigh)

		n = float(len(histograms))
		sn = sqrt(n-1.0)

		for i in range(h.xbins):
			for j in range(h.ybins):
				t = t2 = 0.0
				for h in histograms:
					t  += h[i][j]
					t2 += h[i][j]**2

				new[i][j] = t/n
				new.error(i, j, sqrt(t2/n - new[i][j]**2) / sn)

#		# Under, Over
#		u = 0.0
#		o = 0.0
#		for h in histograms:
#			u  += h.under
#			o  += h.over
#
#		new.under = u/n
#		new.over  = o/n
		return new

# Aliases
H2D = Histogram2D
H2  = Histogram2D
