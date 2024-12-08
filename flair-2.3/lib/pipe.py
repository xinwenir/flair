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
#
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	20-Apr-2006

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import string
import select
import subprocess

#-------------------------------------------------------------------------------
# Execute a system command using pipes
# FIXME to unblock the output buffer
#-------------------------------------------------------------------------------
def system(cmd, inp=[], cwd=None):
	"""
	execute a command by piping the input (list of strings) and
	returning the return code and a list of string of the output
	of the command. Newline character is not required in the input list
	and is removed from the output
	"""

	try:
		pipe = subprocess.Popen(cmd,
				shell=True,
				cwd=cwd,
				stdin=subprocess.PIPE,
				stdout=subprocess.PIPE,
				stderr=subprocess.STDOUT,
				close_fds=True)
	except OSError:
		return 1,["Cannot execute command: %s"%(cmd)]

	out = []
	try:
		rc = 0
		for x in inp:
			ready = select.select([pipe.stdout], [], [], 0)
			while len(ready[0]) > 0:
				# Use the os.read instead of the pipe.stdout.readline
				# to avoid blocking situations
				# readline waits for the \n or EOF
				data = os.read(pipe.stdout.fileno(), 1024)
				if len(data)==0: break
				out.append(data)
				ready = select.select([pipe.stdout], [], [], 0)
			pipe.stdin.write( ("%s\n"%(str(x))).encode() )
			pipe.stdin.flush()
		pipe.stdin.close()
		while True:
			line = pipe.stdout.readline()
			if len(line)==0: break
			out.append(line.decode())
	except OSError:
		rc = 1

#	out = string.join(out,"").splitlines()
	out = "".join(out).splitlines()
	return rc, out
