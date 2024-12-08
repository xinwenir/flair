# 
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
# Date:	13-Jun-2006

__author__  = "Vasilis Vlachoudis"
__email__   = "Vasilis.Vlachoudis@cern.ch"
__version__ = "1"

import os
import sys
import time
import string
import select
import subprocess

from stat import *

import tkFlair
import Project
import Process
import FlairProcess

#===============================================================================
# Data process
#===============================================================================
class DataProcess(FlairProcess.FlairProcess):
	#----------------------------------------------------------------------
	def setUsr(self, usr):
		self.usr = usr

	#----------------------------------------------------------------------
	# Cleaning files in synchronous way
	# @return rc (and output filled)
	#----------------------------------------------------------------------
	def _clean(self):
		rc = 0
		count = 0
		for usr in self.usr:
			fn,ext = os.path.splitext(usr.name())
			files = [fn+ext, fn+"_tab.lis", fn+"_sum.lis"]
			for f in files:
				# Remove files
				if os.path.isfile(f):
					try:
						os.remove(f)
						count += 1
						self.output("Removing: "+f)
					except Exception as err:
						self.output(str(err))
						rc=1
			self.message = ("Data clean", "%d files removed" % (count))
		return rc

	#----------------------------------------------------------------------
	# Compile executable
	#----------------------------------------------------------------------
	def _execute(self):
		Process.Process.execute(self)
		files = []
		errfiles = []
		rc  = 0
		for i,usr in enumerate(self.usr):
			if usr is None:
				self.message = ("Error", "Invalid detector. usr=None")
				return 1
			self.message = ("Processing", usr.name())
			if self.status == Process.STATUS_KILLED:
				self.output("*** Killed ***\n")
				return rc
			if usr is None:
				self.message = ("Warning", "No USRxxx rules found")
				return 1
			if self._processUsrInfo(usr):
				rc = 2
				self.output("Error processing: %s\n"%(usr.name()))
				errfiles.append(usr.name())
			else:
				files.append(usr.name())
			# Update percentage
			self.percent = (i*100) // len(self.usr)

		if rc:
			self.message = ("Warning",
				"Possible errors: %s\n" \
				"Correctly processed: %s"\
					%(" ".join(errfiles)," ".join(files)))
		elif files:
			self.message = ("Data Merging",
				"Files created: %s" % (" ".join(files)))
		else:
			self.message = ("Data Merging", "No files were created")
		return rc

	# --------------------------------------------------------------------
	def _processUsrInfo(self, usr):
		startTime = time.time()-2.0
		# command, input
		cmd,inp = usr.command()
		if cmd is None:
			self.output("\nWARNING: Cannot execute command %s"%(usr.name()))
			return True

		self.output("\nProcessing: %s cmd=%s"%(usr.name(), cmd))

		# start execution
		try:
			pipe = subprocess.Popen(cmd,
					shell=True,
					stdin=subprocess.PIPE,
					stdout=subprocess.PIPE,
					stderr=subprocess.STDOUT,
					close_fds=True)
			for x in inp:
				ready = select.select([pipe.stdout], [], [], 0)
				while len(ready[0]) > 0:
					# Use the os.read instead of the pout.readline
					# to avoid blocking situations
					# readline waits for the \n or EOF
					data = os.read(pipe.stdout.fileno(), 1024)
					if len(data)==0: break
					self.output(data)
					ready = select.select([pipe.stdout], [], [], 0)

				self.output(">>> %s"%(x))
				pipe.stdin.write(("%s\n"%(x)).encode())
				pipe.stdin.flush()
			pipe.stdin.close()

			while True:
				line = pipe.stdout.readline()
				if not line: break
				self.output(line[:-1])

		except (IOError, OSError):
			self.output("ERROR executing command %s"%(cmd))
			self.message = ("Error", "ERROR executing command %s"%(cmd))
			return True

		# Assert that generated file is newer than inp
		verify = usr.verify(startTime)
		if not verify:
			self.output("Processed file is older than some of the dependencies (%g). "\
					"Maybe the run is still going on?")
		return not verify
