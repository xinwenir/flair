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

#ifndef __STREAM_H
#define __STREAM_H

#include <fstream>

enum StreamIO {
	StreamRead,
	StreamWrite
};

/**
 * Geometry streaming I/O class
 */
class Stream {
protected:
	std::string	filename;	/** filename to open	*/
	std::fstream	stream;		/** Input stream	*/
	StreamIO	_io;		/** I/O type		*/
	bool		binary;		/** binary file		*/
	bool		header;		/** is header written	*/
	bool		footer;		/** is footer written	*/

public:
	Stream() : _io(StreamRead) {}
	Stream(const char* name) : filename(name), _io(StreamRead) {}
	Stream(const char* name, StreamIO io=StreamRead, bool bin=false) { open(name, io, bin); }
	~Stream()	{ close(); }

	/** open a new file for output */
	bool	open(StreamIO io=StreamRead, bool bin=false);
	bool	open(const char *name, StreamIO io=StreamRead, bool bin=false);

	/** close file and check headers if are needed */
	void	close();

	/** @return true is file is open for read */
	bool	isRead() const	{ return _io == StreamRead; }
	/** @return true if opened for writing */
	bool	isWrite() const	{ return _io == StreamWrite; }

	/** @return true is file is open */
	bool	isOpen()		{ return stream.is_open(); }

	/** skip entire line */
	void	skipLine();

	/** write a header before any other write statement
	 * @see close
	 * @see writeFooter
	 */
virtual	void	writeHeader();

	/** write a footer before closing the file
	 * @see close
	 * @see writeHeader
	 */
virtual	void	writeFooter();

	/** test if header is needed
	 * @see writeHeader
	 */
virtual	void	write();
}; // Stream

#endif
