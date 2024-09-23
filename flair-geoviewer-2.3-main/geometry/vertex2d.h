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
 * Date:	04-Feb-2010
 */

#ifndef __VERTEX2D_H
#define __VERTEX2D_H

#include "os.h"

class VBody;
class VZone;

enum SegmentType {		// bit mask
	SEGMENT_IGNORE = 0,	// Do not draw anything (ok or error)
	SEGMENT_EMPTY  = 1,	// Do not draw anything (ok or error)
	SEGMENT_ERROR  = 1<<1,	// Draw the segment as error
	SEGMENT_REGION = 1<<2,	// Draw the segment as boundary or regions
	SEGMENT_ZONE   = 1<<3,	// Draw the segment as boundary of zones
	SEGMENT_BODY   = 1<<4	// Draw the segment as body
};

/** A specialized Vertex2D class definition for dealing with 2D
 * conical segment intersections. It holds the location as well
 * the parametric position, body,zone that it are involved
 * */
class Vertex2D {
public:
	double	 t;		// parametric t position of body
	double	 x, y;		// location of vertex
	VBody	*body;		// principle body / conic
	VZone	*zone;		// principle zone
				// needed for displaying zone lines if editRegion is selected
	int	 err;		// error index id
	int	 type;		// segment type
	bool	 invalid;	// if segment needs recalculation

public:
	Vertex2D() {}

	Vertex2D(const double at, const double ax, const double ay, VBody* b) :
		t(at), x(ax), y(ay),
		body(b),
		zone(NULL),
		err(0),
		type(SEGMENT_IGNORE),
		invalid(true)
		{}

//static	bool compare(const Vertex2D& a, const Vertex2D& b) { return a.t < b.t; }
static	int compare(const Vertex2D& a, const Vertex2D& b) { return Cmp(a.t, b.t); }
}; // Vertex2D

#endif
