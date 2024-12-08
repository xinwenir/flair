#!/usr/bin/python
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
# DAMAGES.
#
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	04-Oct-2012

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkFlair
import tkExtra
import Unicode

import Input
import Manual
import Ribbon
import FlairRibbon

#===============================================================================
# Empty frame class
#===============================================================================
class EmptyPage(FlairRibbon.FlairPage):
	"""Empty page for flair"""

	_name_ = "Empty"
	_icon_ = "empty"

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Plot List ===========
		group = Ribbon.LabelGroup(self.ribbon, "Empty")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		group.frame.grid_rowconfigure(0, weight=1)
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clean32"],
				text="Clear",
				compound=TOP,
				#command=self.clear,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Clear")

	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self, "Empty")

	# ----------------------------------------------------------------------
	# Refresh
	# ----------------------------------------------------------------------
	def refresh(self):
		pass

	# ----------------------------------------------------------------------
	# Update fields
	# ----------------------------------------------------------------------
	def get(self, event=None):
		if self.page is None or self.project is None: return
		#self.setModified(True)
