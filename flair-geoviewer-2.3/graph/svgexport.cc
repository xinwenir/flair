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
 * Author:	David.Sinuela.Pastor@cern.ch
 * Date:	Sep-2012
 *
 */

#include <string>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <iostream>

#include "svgexport.h"

using namespace std;

/* =============================== SVGExport =============================== */
SVGExport::SVGExport(const std::string& filename)
	: ExportBase(filename)
{
	if (_file)
		writeHeader();
} // SVGExport

SVGExport::~SVGExport(void)
{
	if (_file) {
		writeEOF();
		_file.close();
	}
} // ~SVGExport

/** attributes */
void SVGExport::attributes(dword color, float width, const char* layer)
{
	_file << " stroke=\"#" << hex << setw(6) << setfill('0') << color << "\"";
	_file << " stroke-width=\"" << width << "\"";
	_file << " fill=\"\"";
	_file << " class=\"layer" << layer << '\"';
} // attributes

/** line */
bool SVGExport::line(double x1, double y1, double x2, double y2, dword color, const char* layer)
{
	_file << "<line"
	      << " x1=\"" << x1 << '\"'
	      << " y1=\"" << y1 << '\"'
	      << " x2=\"" << x2 << '\"'
	      << " y2=\"" << y2 << '\"';
	attributes(color, 0.25, layer);
	_file << "/>" << endl;
	return true;
} // line

/** circle */
bool SVGExport::circle(double x, double y, double radius, dword color, const char* layer)
{
	_file << "<circle"
	      << " cx=\"" << x      << '\"'
	      << " cy=\"" << y      << '\"'
	      << " r=\""  << radius << '\"';
	attributes(color, 0.25, layer);
	_file << "/>" << endl;
	return true;
} // circle

/** arc */
bool SVGExport::arc(double, double, double, double, double, dword, const char*)
{
	return true;
} // arc

/** rectangle */
bool SVGExport::rectangle(double left, double bottom, double right, double top, dword color, const char* layer)
{
	_file << "<rect"
	      << " x=\"" << left            << '\"'
	      << " y=\"" << bottom          << '\"'
	      << " width=\""  << right-left << '\"'
	      << " height=\"" << top-bottom << '\"';
	attributes(color, 0.25, layer);
	_file << "/>" << endl;
	return true;
} // rectangle

/** writeHeader */
bool SVGExport::writeHeader()
{
	_file << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">" << endl;
	return true;
} // WriteHeader

/** WriteEOF */
bool SVGExport::writeEOF()
{
	_file <<  "</svg>" << endl;

	return true;
} // writeEOF
