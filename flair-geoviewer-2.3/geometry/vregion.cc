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

#include <assert.h>
#include <ostream>
#include <string>

#include "vzone.h"
#include "engine.h"
#include "vregion.h"

class GZone;

/** FLUKA viewport region definition */
void VRegion::init(GRegion *reg, GeometryKernel *k)
{
	_region     = reg;
	_generation = -1;
	_kernel     = k;
	location    = false;
	in          = NULL;

	clear();	// FIXME clears the zones.
			// Not really needed we could have reused the zones!
	if (_region) {
		_hash = _region->hash();
		for (int i=0; i<reg->nzones(); i++)
			_zones.add(new VZone(reg->_zones[i],this));
	}
} // init

/** clear */
void VRegion::clear()
{
	for (int i=0; i<_zones.size(); i++)
		if (_zones[i]) delete _zones[i];
	_zones.clear();
} // clear

/** updateLocation
 *
 * find the location of the Region on the Z=0 plane and within the viewport
 * coordinates
 *
 * FIXME: Improvement on STD expressions we can ignore directly the zones
 *        with a negative result, and remove the bodies
 *        Probably add a reference count to the bodies!
 *
 * We have 3 possible locations inside(1) / outside(0) / unknown(?:2)
 *		  a | b
 *	0 | 0 = 0	? | 0 = ?
 *	0 | 1 = 1	? | 1 = 1
 *	1 | 0 = 1	0 | ? = ?
 *	1 | 1 = 1	1 | ? = 1
 */
void VRegion::updateLocation()
{
	location = false;
	in       = NULL;

	bool allknown = true;

	for (int iz = 0; iz < _zones.size(); iz++) {
		VZone* z = _zones[iz];
		z->updateLocation();
		if (z->location) {
			if (z->in) {
				location = true;
				in = z;
			}
		} else
			allknown = false;
	}
	if (allknown && !location) {
		location = true;
		in = NULL;
	}
} // updateLocation

/** inside2D
 * WARNING this routine is optimized for points close to the 2D viewport ONLY!
 * @param x,y,z		location of point to search
 * @param dx,dy,dz	direction in case on boundary
 * @return zone		return zone that location is inside otherwise NULL
 */
VZone *VRegion::inside2D(GeometryEngine* engine,
			 const double  x, const double  y, const double  z,
			 const double dx, const double dy, const double dz) const
{
	if (location) return in;

	for (int iz = 0; iz < _zones.size(); iz++)
		if (_zones[iz]->inside2D(engine, x,y,z, dx,dy,dz))
			return _zones[iz];
	return NULL;
} // inside2D

/** inside
 * @param x,y,z		location of point to search
 * @param dx,dy,dz	direction in case on boundary
 * @return zone		return zone that location is inside otherwise NULL
 */
VZone *VRegion::inside(GeometryEngine* engine,
			 const double  x, const double  y, const double  z,
			 const double dx, const double dy, const double dz) const
{
	for (int iz = 0; iz < _zones.size(); iz++)
		if (_zones[iz]->inside(engine, x,y,z, dx,dy,dz))
			return _zones[iz];
	return NULL;
} // inside

/** insideRay */
VZone *VRegion::insideRay(GeometryEngine* engine,
			  const double  x, const double  y, const double  z,
			  const double dx, const double dy, const double dz,
			  const double t) const
{
	//const Vector p(x+dx*t, y+dy*t, z+dz*t);
	//if (!obbox().inside(p)) return NULL;
	for (int iz = 0; iz < _zones.size(); iz++) {
		if (_zones[iz]->insideRay(engine, x,y,z, dx,dy,dz, t))
			return _zones[iz];
	}
	return NULL;
} // insideRay

#if 0
/** intersectRay
 * @param tmin	I/O return min distance to the next intersection inside region
 * @return body where exiting intersection was found
 */
CBody* VRegion::intersectRay(const double  x, const double  y, const double  z,
			    const double dx, const double dy, const double dz,
			    double *tmin, const double tmax) const
{
	VBody *body = NULL;
	for (unsigned iz = 0; iz < _zones.size(); iz++) {
		body = _zones[iz]->intersectRay(x,y,z,dx,dy,dz,tmin,tmax);
		if (body) return body;
	}
	return NULL;
} // intersectRay
#endif

/** @return memory used by region */
size_t VRegion::memory() const
{
	size_t zmem = sizeof(VRegion);
	for (int iz = 0; iz < _zones.size(); iz++)
		zmem += _zones[iz]->memory();
	return zmem;
} // memory

/** operator << */
std::ostream& operator << (std::ostream& s, const VRegion& region)
{
	s << "VRegion: " << region.name() << std::endl;
#if 0
	for (int i=0; i<region.n; i++) s << region[i];
	//for (int i=0; i<region.expr.count(); i++) {
		VBody *body = (VBody*)region.expr[i];
		if (body == tplus)		s << "+ " << std::endl;
		else if (body == tminus)	s << "- " << std::endl;
		else if (body == tunion)	s << "| " << std::endl;
		else if (body == tuniverse)	s << "@ " << std::endl;
		else
			s << body->name() << " ";
	}
#endif
	return s;
} /* operator << */
