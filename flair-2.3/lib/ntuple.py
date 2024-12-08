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
import math
import bmath
from math import *
from collections import namedtuple

# Populate math functions on the local dict
_globalDict = {}
for name in dir(math):
	if name[0]=="_": continue
	_globalDict[name] = getattr(math, name)

#===============================================================================
# NTuple
#===============================================================================
class NTuple:
	"""NTuple"""

	# ----------------------------------------------------------------------
	def __init__(self, variables):
		"""
		Initialize the Ntuple with the variables
		defined as a stirng:  name[/type][:name/[type]]...
		or a list of [name/type...] or a list of tuples
		types are:
			f : floating point or nothing
			i : integer
			s : string
		the self.Record contains the namedtuple definition for each row

		Example:
			t = NTuple("n/i x y z xp yp")
			t.load("psblack001_dump")
			h = t.project(H1(100,-2.0,2.0), "x")
			gnuplot = Gnuplot()
			h.draw(gnuplot)
			raw_input()
			gnuplot.plot("u 1:2 w p not", t.scatter("xp:x"))
		"""
		if isinstance(variables, str):
			if ":" in variables:
				variables = variables.split(":")
			else:
				variables = variables.split()

		self.names    = []	# names of variables
		self.types    = []	# types of variables
		self._toint   = []	# which columns should be converted to int
		self._tofloat = []	# which columns should be converted to float

		for i,var in enumerate(variables):
			if isinstance(var,tuple):
				n,t = var
			elif "/" in var:
				n,t = var.split("/")
				t = t.lower()
			else:
				n = var
				t = "f"
			self.names.append(n)
			self.types.append(t)
			if t=="f":
				self._tofloat.append(i)
			elif t=="i":
				self._toint.append(i)

		self.Record = namedtuple("Record",self.names)
		self.rows = []

	# ----------------------------------------------------------------------
	def __getitem__(self, i):
		"""return ith row of data"""
		return self.rows[i]

	# ----------------------------------------------------------------------
	def events(self):
		"""return number of events in Ntuple"""
		return len(self.rows)

	# ----------------------------------------------------------------------
	def columns(self):
		"""return number of columns (variables)"""
		return len(self.names)

	# ----------------------------------------------------------------------
	def append(self, values):
		"""append a new record"""
		self.rows.append(self.Record(*values))

	# ----------------------------------------------------------------------
	# Load an ascii file
	# ----------------------------------------------------------------------
	def load(self, filename, delimiter=None, usecols=None, skiprows=0):
		"""Load text file in tuple"""
		if isinstance(filename,str):
			f = open(filename,"r")
		else:
			f = filename
		n0 = len(self.rows)
		for line in f:
			if skiprows>0:
				skiprows -= 1
				continue
			line = line.strip()
			if not line: continue
			if line[0] in ("#","*"): continue
			val = line.split(delimiter)
			for i in self._toint:
				val[i] = int(val[i])
			for i in self._tofloat:
				val[i] = float(val[i])

			if usecols is not None:
				if isinstance(usecols,int):
					val = [val[usecols]]
				else:
					val2 = []
					for c in usecols:
						val2.append(val[c])
					val = val2
			self.rows.append(self.Record(*val))
		return len(self.rows)-n0
	readFile = load

	# ----------------------------------------------------------------------
	def clear(self):
		"""Delete data"""
		del self.rows[:]

	# ----------------------------------------------------------------------
	# Dump to stdout
	# ----------------------------------------------------------------------
	def dump(self):
		"""dump NTuple on stdout"""
		print("#","\t".join(["%d:%s/%s"%(i,n,t)
			for i,(n,t) in enumerate(zip(self.names,self.types))]))
		for i,row in enumerate(self.rows):
			print(i,"\t".join(map(str,row)))

	# ----------------------------------------------------------------------
	def min(self, column):
		"""return minimum value of column"""
		mini = 1e999
		for row in self.rows:
			mini = min(mini, row[column])
		return mini

	# ----------------------------------------------------------------------
	def max(self, column):
		"""return maximum value of column"""
		maxi = -1e999
		for row in self.rows:
			maxi = max(maxi, row[column])
		return maxi

	# ----------------------------------------------------------------------
	@staticmethod
	def _compile(expr, weight=None):
		if ":" in expr:
			# Multi dimensional
			expry, exprx = expr.split(":")
			exprx = compile("float(%s)"%(exprx),"<exprx>","eval")
			expry = compile("float(%s)"%(expry),"<expry>","eval")
		else:
			exprx = None
			expry = compile("float(%s)"%(expr),"<expry>","eval")

		if isinstance(weight,str):
			weightc = compile("float(%s)"%(weight),"<weight>","eval")
		else:
			weightc = None
		return exprx, expry, weightc

	# ----------------------------------------------------------------------
	# Project to a histogram the expression, applying a weight
	# ----------------------------------------------------------------------
	def project(self, hist, expr, weight=1.0, from_=0, to_=-1):
		"""fill 1D or2D histogram hist with the expression expr"""
		exprx, expry, weightc = NTuple._compile(expr, weight)

		for row in self.rows[from_:to_]:
			if weightc is not None:
				weight = eval(weightc, _globalDict, row._asdict())
			if weight>0.0:
				valuey = eval(expry, _globalDict, row._asdict())
				if exprx is None:
					hist.fill(valuey, weight)
				else:
					valuex = eval(exprx, _globalDict, row._asdict())
					hist.fill(valuex, valuey, weight)
		return hist

	# ----------------------------------------------------------------------
	def scatter(self, expr, weight=1.0, from_=0, to_=-1):
		"""return a list of two columns with the scatter data for plotting"""
		exprx, expry, weightc = NTuple._compile(expr, weight)

		output = []
		for row in self.rows[from_:to_]:
			if weightc is not None:
				weight = eval(weightc, _globalDict, row._asdict())
			if weight>0.0:
				valuey = eval(expry, _globalDict, row._asdict())
				if exprx is None:
					output.append((valuey, weight))
				else:
					valuex = eval(exprx, _globalDict, row._asdict())
					output.append((valuex, valuey, weight))
		return output

	# ----------------------------------------------------------------------
	def evaluate(self, expr):
		"""return a list with the evaluation of the expression expr"""
		exprc = compile("float(%s)"%(expr),"<expr>","eval")
		output = []
		for row in self.rows:
			output.append(eval(exprc, _globalDict, row._asdict()))
		return output

#-------------------------------------------------------------------------------
if __name__ == "__main__":
	from Gnuplot import Gnuplot
	from histogram import H1

	t = NTuple("n/i x y z xp yp")
	t.load("psblack001_dump")
	for i in range(t.columns()):
		print(i, t.min(i), t.max(i))
	#t.dump()
	h = t.project(H1(100,-2.0,2.0), "x")
	gnuplot = Gnuplot()
	#h.draw(gnuplot)
	#raw_input()
	gnuplot.grid()
	gnuplot.plot("u 1:2 w p not", t.scatter("xp:x"))
	input()
