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
 * BasedOn:	Chris.Theis@cern.ch
 * Date:	unknown
 *
 */
#ifndef DXF_EXPORT_H
#define DXF_EXPORT_H

#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>

#include "color.h"
#include "vector.h"
#include "exportbase.h"

class DXFExport : public ExportBase {
public:
	DXFExport(const std::string& filename) : ExportBase(filename) {
			if (_file.is_open()) writeHeader();
			_lastColor = 0;
			_best      = 0;
		}
	virtual ~DXFExport(void) {
			if (_file.is_open()) writeEOF();
		}
	bool line(double x1, double y1, double x2, double y2, dword color=0, const char* layer=NULL);
	bool point(double x, double y, dword color=0, const char* layer=NULL);
	bool circle(double x, double y, double radius, dword color=0, const char* layer=NULL );
	bool arc(double x, double y, double radius, double startPhi, double endPhi, dword color=0,  const char* layer=NULL);
	bool polyline(Path& path, dword color=0, const char* layer=NULL);

protected:
	bool    writeHeader();
	bool    writeEOF();

	void	write(const int id, const int value) {
			_file << std::setw(3) << id << std::endl;
			_file << std::setw(0) << value << std::endl;
		}

	void	write(const int id, const double value) {
			_file << std::setw(3) << id << std::endl;
			_file << std::setw(0) << value << std::endl;
		}

	void	write(const int id, const char* str) {
			_file << std::setw(3) << id << std::endl;
			_file << str << std::endl;
		}

	void	write(const int id, const double x, const double y) {
			write(10+id, x);
			write(20+id, y);
		}

	void	write(const int id, const double x, const double y, const double z) {
			write(10+id, x);
			write(20+id, y);
			write(30+id, z);
		}

	void	write2D(const int id, const Vector& v) {
			write(10+id, v.x);
			write(20+id, v.y);
		}

	void	write(const int id, const Vector& v) {
			write(10+id, v.x);
			write(20+id, v.y);
			write(30+id, v.z);
		}

	void	writeLayer(const char *layer=NULL) {
			if (layer) write(8,layer);
		}

	void	writeColor(const int color) {
			if (color>=0) write(62, closestColor(color));
		}

	int	closestColor(const dword color);

protected:
	dword	_lastColor;
	int	_best;
}; // DXFExport

#endif

