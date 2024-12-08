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
 * Date:	07-Mar-2016
 */

#ifndef __SCATTER_H
#define __SCATTER_H

#include "array.h"

//-----------------------------------------------------------------------------
// Pair x,y value
//-----------------------------------------------------------------------------
class XYPair {
public:
	double	x;
	double	y;
public:
	XYPair(const double ax=0.0, const double ay=0.0) : x(ax), y(ay) {}
	XYPair(const XYPair& p) : x(p.x), y(p.y) {}
static	int compare(const XYPair& a, const XYPair& b) { return Cmp(a.x, b.x); }

	// not used but to keep compiler happy
	bool	operator==(const XYPair& p) const	{ return x == p.x; }
};

//=============================================================================
// Scatter data class
//=============================================================================
class Scatter {
private:
	Array<XYPair>	data;
public:
	Scatter() { data.compare(XYPair::compare); }
	~Scatter();

	bool	load(const char* filename);
	bool	isEmpty() const				{ return data.size()==0; }

const	XYPair&	operator[](const int i)	const		{ return data[i]; }

	double	interpolate(const double x) const;
	double	interpolateLog(const double x) const;
	double	operator()(const double x) const	{ return interpolate(x); }
}; // class Scatter

#endif
