/*
 * $Id$
 *
 * Copyright and User License
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright 2006-2019 CERN and INFN
 * 
 *
 * Please consult the LICENSE file for the license 
 *
 * DISCLAIMER
 * ~~~~~~~~~~
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
 * NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY, OF
 * SATISFACTORY QUALITY, AND FITNESS FOR A PARTICULAR PURPOSE
 * OR USE ARE DISCLAIMED. THE COPYRIGHT HOLDERS AND THE
 * AUTHORS MAKE NO REPRESENTATION THAT THE SOFTWARE AND
 * MODIFICATIONS THEREOF, WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADE SECRET OR OTHER PROPRIETARY RIGHT.
 *
 * LIMITATION OF LIABILITY
 * ~~~~~~~~~~~~~~~~~~~~~~~
 * THE COPYRIGHT HOLDERS AND THE AUTHORS SHALL HAVE NO
 * LIABILITY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
 * CONSEQUENTIAL, EXEMPLARY, OR PUNITIVE DAMAGES OF ANY
 * CHARACTER INCLUDING, WITHOUT LIMITATION, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES, LOSS OF USE, DATA OR PROFITS,
 * OR BUSINESS INTERRUPTION, HOWEVER CAUSED AND ON ANY THEORY
 * OF CONTRACT, WARRANTY, TORT (INCLUDING NEGLIGENCE), PRODUCT
 * LIABILITY OR OTHERWISE, ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGES.
 *
 * Author:	Vasilis.Vlachoudis@cern.ch
 * Date:	11-Jan-2006
 */

#include <assert.h>
#include "stream.h"

using namespace std;

/* --- open --- */
bool Stream::open(StreamIO io, bool bin)
{
	_io    = io;
	binary = bin;
	header = false;
	footer = false;
	if (isOpen()) close();
	if (binary)
		stream.open(filename.c_str(),
			(isWrite()? ios::out : ios::in | ios::binary));
	else
		stream.open(filename.c_str(),
			(isWrite()? ios::out : ios::in));
	return isOpen();
} // open

/* --- open --- */
bool Stream::open(const char* fn, StreamIO io, bool bin)
{
	filename = fn;
	return open(io, bin);
} // open

/* --- close --- */
void Stream::close()
{
	if (isOpen()) {
		if (isWrite())
			if (header && !footer) writeFooter();
		stream.close();
	}
} // close

/** skipLine */
void Stream::skipLine()
{
	int c;
	do {
		c = stream.get();
	} while (c!='\n' && c!=-1);
} // skipLine

/* --- write --- */
void Stream::write()
{
	assert(isWrite());
	if (!header)
		writeHeader();
} // write

/* --- writeHeader --- */
void Stream::writeHeader()
{
	assert(isWrite());
	assert(isOpen());
	assert(header==false);

	header = true;
} // writeHeader

/* --- writeFooter --- */
void Stream::writeFooter()
{
	assert(isWrite());
	assert(footer==false);

	footer = true;
} // writeFooter
