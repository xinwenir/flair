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
 * Date:	11-Oct-2016
 */

#ifndef __TETRAMESHIO_H
#define __TETRAMESHIO_H

#include "stream.h"
#include "tetramesh.h"

/**
 * TetraMesh I/O helper class
 */
class TetraMeshIO : public Stream {
public:
	TetraMeshIO() : Stream() {}
	TetraMeshIO(const char *name, StreamIO io=StreamRead, bool bin=false)
		: Stream(name, io, bin) {}
	virtual ~TetraMeshIO() {close();}

	bool	read(TetraMesh& mesh);
//	bool	write(TetraMesh& mesh);

#if 0
	/** prepare a valid DXF header */
	void writeHeader();
	/** prepare a valid DXF footer */
	void writeFooter();

	/** write a mesh and check if header is needed
	 * @param mesh mesh to be written
	 * @param color color code to use
	 * @param name override name to use
	 */
	virtual void write(Mesh& mesh, int color=0, const char *name=NULL);

private:
	/** write a vector in DXF format */
	void writeVector(const double x, const double y, const double z);
	void writeVector(const Vector& v)	{ writeVector(v.x, v.y, v.z); }
#endif

private:
	bool	_readAscii(TetraMesh& mesh);
	bool	_readBinary(TetraMesh& mesh);
}; // TetraMeshIO

#endif
