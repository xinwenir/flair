#!bin/env python
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
# ~~~~~~~~~~~~~~~~~~~~~t)
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
# Date:	03-Mar-2015

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os

#===============================================================================
# Cached listdir
#===============================================================================
_cache = {}

#-------------------------------------------------------------------------------
def listdir_reset():
	global _cache
	_cache = {}

#-------------------------------------------------------------------------------
def listdir(path):
	global _cache

	"""List directory contents, using cache."""
	try:
		cached_mtime, lst = _cache[path]
		del _cache[path]
	except KeyError:
		cached_mtime, lst = -1, []

	mtime = os.stat(path).st_mtime
	if mtime != cached_mtime:
		try:
			lst = os.listdir(path)
		except OSError:
			lst = []
#		lst.sort()

	_cache[path] = mtime, lst
	return lst

#-------------------------------------------------------------------------------
# Display friendly time
#-------------------------------------------------------------------------------
def friendlyTime(t):
	if t<1.0:
		return "%g ms" % (t*1000.0)
	elif t<60.0:
		return "%g s" % (t)

	(days,  tm) = divmod(t, 86400)
	(hours, tm) = divmod(tm, 3600)
	(minu, sec) = divmod(tm,   60)

	if days>0:
		s = "%dd" % (days)
	else:
		s = ""
	if hours>0 or s!="":
		if s!="": s += " "
		s += "%dh" % (hours)
	if minu>0 or s!="":
		if s!="": s += " "
		s += "%dm" % (minu)
	if s!="": s += " "

	s += "%ds" % (int(sec))
	return s

#-------------------------------------------------------------------------------
# Search all matching words in the string
#-------------------------------------------------------------------------------
def search(s, needle):
	for w in needle:
		if s.find(w) < 0: return False
	return True

#-------------------------------------------------------------------------------
# Check if the first file is newer than the second
#-------------------------------------------------------------------------------
def isNewer(newfile, oldfile):
	try:
		n = os.stat(newfile).st_mtime
		o = os.stat(oldfile).st_mtime
		return n > o
	except:
		return False
