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

import sys
import string
import threading
import subprocess
from log import say

# Constants defining the status
STATUS_IDLE     = 0
STATUS_WAITING  = 1
STATUS_RUNNING  = 2
STATUS_FINISHED = 3
STATUS_ERROR    = 4
STATUS_KILLED   = 5

#           0        1           2          3          4         5
STATUS = ["Idle", "Waiting", "Running", "Finished", "Error", "Killed"]
COLOR  = ["Black", "Orange", "Blue",    "DarkGreen", "Red", "DarkRed"]

#===============================================================================
# Base Process
# The base class for all offering monitoring tools for the state of a system
# process e.g. Compile, Running, Data merging, etc.
#===============================================================================
class Process:
	#----------------------------------------------------------------------
	def __init__(self, notify=None):
		self._notify = notify
		self._mutex  = threading.Semaphore()
		self.clearStatus()

	#----------------------------------------------------------------------
	def clearStatus(self):
		self._thread = None
		self.rc      = 0
		self._out    = []
		self.message = ("","")
		self.percent = 0
		self.status  = STATUS_IDLE

	#----------------------------------------------------------------------
	def statusMessage(self):
		return STATUS[self.status]

	#----------------------------------------------------------------------
	def start(self):
		if self.status != STATUS_IDLE:
			say("Process is not Idle")
			return False
		if self._thread and self._thread.isAlive():
			say("Process is still alive")
			return

		self.rc      = None
		self.percent = 0
		self.status  = STATUS_WAITING
		self._thread = threading.Thread(target=self._executeThread)
		self._thread.start()
		return True

	#----------------------------------------------------------------------
	def _executeThread(self):
		try:
			self.rc = self.execute()
		except:
			self.message = ("Exception", sys.exc_info()[1])
			self.status  = STATUS_ERROR
			raise
		else:
			if self.rc != 0:
				self.status = STATUS_ERROR
			else:
				self.percent = 100
				self.status  = STATUS_FINISHED
		if self._notify:
			self._notify(self)

	#----------------------------------------------------------------------
	def isAlive(self):
		if self.status in (STATUS_WAITING, STATUS_RUNNING):
			return True
		if self._thread is not None:
#			return self._thread.isAlive()
			return self._thread.is_alive()
		return False

	#----------------------------------------------------------------------
	def kill(self):
		self.status = STATUS_KILLED

	#----------------------------------------------------------------------
	def isKilled(self):
		return self.status == STATUS_KILLED

	#----------------------------------------------------------------------
	def wait(self, timeout=None):
		self._thread.join(timeout)

	#----------------------------------------------------------------------
	# Method to override with the execution of the task
	#----------------------------------------------------------------------
	def execute(self):
		self.status = STATUS_RUNNING
		return 0

	#----------------------------------------------------------------------
	def clean(self):
		return 0

	#----------------------------------------------------------------------
	def lock(self):
		self._mutex.acquire()

	#----------------------------------------------------------------------
	def release(self):
		self._mutex.release()

	#----------------------------------------------------------------------
	# Return/Set output
	#----------------------------------------------------------------------
	def output(self, out=None):
		self._mutex.acquire()
		if out is None:
#			out = " ".join(self._out.encode("utf-8"))
#			out+="\n"
			out="\n"
			for i in range(len(self._out)):
				aaa=isinstance( self._out[i],str)
				if aaa:                                
					out += self._out[i]
				else:                                
					out += self._out[i].decode()
				out+= "\n"
			del self._out[:]
		else:
			self._out.append(out)
			out = None
		self._mutex.release()
		return out

#------------------------------------------------------------------------------
class MyProcess(Process):
	def __init__(self, path):
		Process.__init__(self)
		self.path = path

	#----------------------------------------------------------------------
	def execute(self):
		self.p = subprocess.Popen(["/bin/ls","-lR",self.path], stdout=subprocess.PIPE)
		self.stdout = self.p.stdout
		Process.execute(self)
		return self.p.wait()

	#----------------------------------------------------------------------
	def kill(self):
		Process.kill(self)
		self.p.terminate()

#------------------------------------------------------------------------------
if __name__ == "__main__":
	proc = MyProcess("/")
	proc.start()
	while proc.status == STATUS_WAITING:
		say(proc.status)
	try:
		while proc.status == STATUS_RUNNING:
			say(">1>",proc.stdout.readline())
	except KeyboardInterrupt:
		proc.kill()
	for line in proc.stdout:
		say(">2>",proc.status,line)
	say("Status=",proc.status,proc.rc)
