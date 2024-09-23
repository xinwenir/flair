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
 */

#ifndef __VZONE_H
#define __VZONE_H

#include <iostream>

#include "gzone.h"
#include "vector.h"
#include "vregion.h"

class GBody;
class CBody;
class OBBox;
class VBody;
class GeometryEngine;

class VZone {
private:
	GZone*		_zone;		// Zone definition
	VRegion*	_region;	// parent region
	int		_generation;	// Generation number, to be checked for
					// validity against the _zone generation
public:
	bool	location;		// true if in is correct
	bool	in;			// if window is inside or outside the region
					// condition to be checked ONLY if location=true

public:
	VZone(GZone *z=NULL, VRegion *r=NULL)	{ init(z,r); }
	void	init(GZone *z=NULL, VRegion *r=NULL);

	bool	valid()		const { return _generation == _zone->generation(); }
	bool	invalid()	const { return _generation != _zone->generation(); }
	void	setInvalid()          { _generation = -1;  }
	void	setValid()	      { _generation = _zone->generation(); }

	GZone*	zone()		const { return _zone; }
	VRegion* region()	const { return _region; }
const	GRegion* gregion()	const { return _region->region(); }
const	char*	name()		const { return _region->name(); }
	int	id()		const { return _zone->id(); }

	int	size()		const { return _zone->size();  }
	bool	rpn()		const { return _zone->rpn();   }
const	GBody*	gexpr(int i)	const { return _zone->expr[i]; }
inline	VBody*	vbody(const GBody* body) const;

	bool	inside2D(GeometryEngine* engine,
			 const double  x, const double  y, const double  z,
			 const double dx, const double dy, const double dz) const;
	bool	inside(GeometryEngine* engine,
			 const double  x, const double  y, const double  z,
			 const double dx, const double dy, const double dz) const;
	bool	insideRay(GeometryEngine* engine,
			 const double  x, const double  y, const double  z,
			 const double dx, const double dy, const double dz,
			 const double  t) const;
	bool	insideRay(GeometryEngine* engine,
			  const Vector& p, const Vector &d, double t) const
			{ return insideRay(engine, p.x,p.y,p.z, d.x,d.y,d.z, t); }

	CBody*	intersectRay(GeometryEngine* engine,
			     const double  x, const double  y, const double  z,
			     const double dx, const double dy, const double dz,
			     double *tmin, const double tmax) const;
	CBody*	intersectRay(GeometryEngine* engine,
			     const Vector& p, const Vector &d,
			     double *tmin, const double tmax) const
			{ return intersectRay(engine, p.x,p.y,p.z, d.x,d.y,d.z, tmin, tmax); }

	void	updateLocation();
	bool	ignore()	const { return (location && !in); }

	OBBox*	obbox() const { return _zone->obbox(); };

	size_t	memory() const { return sizeof(VZone); }

private:
	bool	_inside2D(GeometryEngine* engine,
			  const double  x, const double  y, const double  z,
			  const double dx, const double dy, const double dz) const;

inline	static CBody* cbody(GeometryEngine* engine, const GBody* body);
}; // VZone

/** inside */
inline bool VZone::inside2D(GeometryEngine* engine,
			    const double  x, const double  y, const double  z,
			    const double dx, const double dy, const double dz) const
{
	if (location) return in;
#if 0
//#if _DEBUG>1 && defined(EXPERIMENTAL)
	bool ginside = _zone->inside(x,y,z, dx,dy,dz);
	bool zinside = _inside2D(engine, x,y,z, dx,dy,dz);
	if (ginside != zinside)
		std::cout << "*** VZone::inside2D ERROR ***" << endl;
//	assert(ginside == zinside);
#endif
	return _inside2D(engine, x,y,z, dx,dy,dz);
} /* inside2D */

std::ostream& operator << (std::ostream&, const VZone&);

#endif
