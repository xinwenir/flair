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
#ifndef SVG_EXPORT_H
#define SVG_EXPORT_H

#include <iostream>
#include <fstream>
#include <string>

#include "exportbase.h"

class SVGExport : public ExportBase {
public:
	SVGExport(const std::string& filename);
	virtual ~SVGExport(void);

	bool line(double x1, double y1, double x2, double y2, dword color=0, const char* layer=NULL);
//	bool path(double x1, double y1, double x2, double y2, dword color=0, const char* layer=NULL);
	bool circle(double x, double y, double radius, dword color=0, const char* layer=NULL);
	bool arc(double x, double y, double radius, double startPhi, double endPhi, dword color=0,  const char* layer=NULL);
	bool rectangle(double left, double bottom, double right, double top, dword color=0, const char* layer=NULL);

protected:
	bool	writeHeader();
	bool	writeEOF();

private:
	void	attributes(dword color, float width, const char* layer);
}; // SVGExport

#endif
