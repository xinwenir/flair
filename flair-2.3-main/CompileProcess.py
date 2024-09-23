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
import subprocess

from bFileDialog import _TIME_FORMAT

import Project
import Process
import FlairProcess

#===============================================================================
# Compile FLUKA project process
#===============================================================================
class CompileProcess(FlairProcess.FlairProcess):
	#----------------------------------------------------------------------
	# Cleaning files in synchronous way
	# @return rc (and output filled)
	#----------------------------------------------------------------------
	def _clean(self):
		# Remove files
		rc = 0
		count = 0
		try:
			os.remove(self.project.exe)
			self.output("Removing: "+self.project.exe)
			count += 1
### XXX PYTHON3 XXX
		except OSError:
			errno, strerror = sys.exc_info()[1]
			if errno!=2:	# Ignore No such file or directory
				self.output(strerror)
		except:
			self.output(str(sys.exc_info()[1]))
			self.message = ("Compile clean", "Cannot delete executable")
			rc = 1

		try:
			os.remove(self.project.exe+".map")
			self.output("Removing: %s.map"%(self.project.exe))
			count += 1
### XXX PYTHON3 XXX
		except OSError:
			errno, strerror = sys.exc_info()[1]
			if errno!=2:	# Ignore No such file or directory
				self.output(strerror)
		except:
			self.output(str(sys.exc_info()[1]))
			self.message = ("Compile clean", "Cannot delete map file")
			rc = 2

		for f in self.project.sourceList:
			(fn, ext) = os.path.splitext(f)
			# Skip libraries
			if ext in (".a",".so",".o"): continue
			obj = fn + ".o"
			try:
				os.remove(obj)
				count += 1
				self.output("Removing: " + obj)
### XXX PYTHON3 XXX
			except OSError:
				errno, strerror = sys.exc_info()[1]
				if errno!=2:	# Ignore No such file or directory
					self.output(strerror)
			except:
				self.output(str(sys.exc_info()[1]))
				self.message = ("Compile clean", "Cannot delete map %s"%(obj))
				rc = 3
		if rc == 0:
			if count:
				self.message = ("Compile clean", "%d files deleted"%(count))
			else:
				self.message = ("Compile clean", "No files deleted")
		return rc

	#----------------------------------------------------------------------
	# Compile executable
	#----------------------------------------------------------------------
	def _execute(self):
		# compile all sources that need compilation
		try: exe_time = os.stat(self.project.exe).st_mtime
		except: exe_time = 0

		libfluka = os.path.join(Project.flukaDir,Project.DEFAULT_LIB)
		try:
			lib_time = os.stat(libfluka).st_mtime
		except:
			lib_time = 0

		maxobj_time = 0
		tot = len(self.project.sourceList) + 1
		for i,f in enumerate(self.project.sourceList):
			self.message = ("Compiling", f)
			if self.status == Process.STATUS_KILLED:
				self.output("*** Killed ***")
				return 0
			(fn, ext) = os.path.splitext(f)
			if ext in (".a",".so",".o"):
				# Check date of libraries
				try: o_time = os.stat(f).st_mtime
				except: o_time = 0
			else:
				# or object
				try: o_time = os.stat(fn+".o").st_mtime
				except: o_time = 0

				try: f_time = os.stat(f).st_mtime
				except:
					self.message = ("File do not exist",
							"Cannot find file %s" %(f))
					return 1

				if o_time <= f_time or o_time <= lib_time:
					# Compile Fortran
					cmd = self.project.compileCmd(f)
					self.output(">>> Compiling: "+" ".join(list(cmd)))
					#self.output(subprocess.check_output(cmd,
					#		stderr=subprocess.STDOUT))
					try:
						p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
					except:
						self.message = ("Error", sys.exc_info()[1])
						return 2
					out,err = p.communicate()
					self.output(out)
					try: o_time = os.stat(fn+".o").st_mtime
					except: o_time = 0
					if o_time < f_time:
						if o_time == 0:
							self.message = ("Error compiling",
									"Error: compiling file: %s\n" \
									"No object generated" %(f))
						else:
							self.message = ("Error compiling",
									"Error: file %s time: %s\n" \
									"is newer than the object %s" \
									%(f,
									 time.strftime(_TIME_FORMAT, time.localtime(f_time)),
									 time.strftime(_TIME_FORMAT, time.localtime(o_time))))
						return 2

				if self.isKilled():
					self.message = ("Compilation Stopped",
							"User killed the compilation")
					return 3
			if o_time > maxobj_time:
				maxobj_time = o_time
			# Update percentage
			self.percent = (i*100) // tot

		# Link executable
		if exe_time <= maxobj_time:
			cmd = self.project.linkCmd()
			self.output(">>> Linking: "+" ".join(list(cmd)))
			#self.output(subprocess.check_output(cmd,
			#			stderr=subprocess.STDOUT))
			p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			out,err = p.communicate()
			self.output(out)
			try:
				exe_time = os.stat(self.project.exe).st_mtime
				if exe_time < maxobj_time:
					self.message = ("Error linking",
						"Error: executable %s time: %s\n" \
						"is older than objects %s" \
							%(self.project.exe,
							 time.strftime(_TIME_FORMAT, time.localtime(exe_time)),
							 time.strftime(_TIME_FORMAT, time.localtime(maxobj_time))))
					return 4
			except:
				self.message = ("Error linking", "Error: executable is not generated")
				return 4
			try:
				exe_time = os.stat(self.project.exe).st_mtime
				self.message = ("Successful Build",
					"Executable %s is successfully built" \
					% (self.project.exe))
				return 0
			except:
				self.message = ("Error linking",
					"Error linking executable "+self.project.exe)
				return 5
		else:
			self.message = ("Up to date", "All files are up to date")
			return 0
