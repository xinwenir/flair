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
# Date:	14-Nov-2014

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

#==============================================================================
# cmp function
#==============================================================================
def cmp(a,b):
	return (a>b) - (a<b)

#==============================================================================
# Convert a compare(mycmp) to a key for list sorting
#==============================================================================
def cmp_to_key(mycmp):
	'Convert a cmp= function into a key= function'
	class K:
		def __init__(self, obj, *args):
			self.obj = obj
		def __lt__(self, other):
			return mycmp(self.obj, other.obj) < 0
		def __gt__(self, other):
			return mycmp(self.obj, other.obj) > 0
		def __eq__(self, other):
			return mycmp(self.obj, other.obj) == 0
		def __le__(self, other):
			return mycmp(self.obj, other.obj) <= 0
		def __ge__(self, other):
			return mycmp(self.obj, other.obj) >= 0
		def __ne__(self, other):
			return mycmp(self.obj, other.obj) != 0
	return K
