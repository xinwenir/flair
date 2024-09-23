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
 * Date:	2012
 *
 */

#ifndef EXPORTBASE_H
#define EXPORTBASE_H

#include <iostream>
#include <fstream>
#include <string>

#include "os.h"
#include "array.h"
#include "vector.h"

typedef Array<Vector2D>	Path;

class ExportBase {
protected:
	std::ofstream	_file;
	std::string	_fileName;

protected:
	// prevent copies
	ExportBase( const ExportBase& rhs );
	ExportBase& operator=( const ExportBase& rhs );

public:
	ExportBase(const std::string& filename ) : _fileName(filename)
			{ _file.open(filename.c_str()); };
	virtual ~ExportBase()	{ if (_file) _file.close(); }

	operator bool () const	{ return _file.is_open(); } // overload the boolean conversion op

	// prevent instantiation of base class & enforce minimum exporter features
	virtual bool line(double x1, double y1, double x2, double y2, dword color=0, const char* layer=NULL)=0;
	virtual bool circle(double x, double y, double radius, dword color=0, const char* layer=NULL)=0;
	virtual bool arc(double x, double y, double radius, double startPhi, double endPhi, dword color=0, const char* layer=NULL)=0;
	virtual bool point(double x, double y, dword color=0, const char* layer=NULL) {
				return line(x,y,x,y,color,layer);
			}
	virtual bool polyline(Path& path, dword color=0, const char* layer=NULL) {
				bool b = true;
				for (int i=0; i<path.size()-1; i++)
					b &= line(path[i].x, path[i].y, path[i+1].x, path[i+1].y, color, layer);
				return b;
			}
	virtual bool rectangle(double left, double bottom, double right, double top, dword color=0, const char* layer=NULL) {
				return	line(left, bottom, right, bottom, color, layer) &&
					line(right, bottom, right, top, color, layer)   &&
					line(right, top, left, top, color, layer)       &&
					line(left, top, left, bottom, color, layer);
			}

}; // ExportBase

#endif
