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

#include <assert.h>
#include <string.h>

#include "gzone.h"
#include "obbox.h"
#include "gregion.h"
class GBody;

using namespace std;

/* ================================ GREGION ================================= */
/** FLUKA region definition */
GRegion::GRegion(const char *aname) :
	_id(-1),
	_type(REGION_NORMAL),
	_generation(0),
//	_volume(0.0),
	_hash(0),
	cached_obbox(NULL),
	_material(NULL),
	rotdefi(0),
	show(0)
{
	name(aname);
} // GRegion

/** name - set name of body */
void GRegion::name(const char *aname)
{
	strncpy(_name, aname, sizeof(_name));
	_name[sizeof(_name)-1] = 0;
} // name

/** @return memory used by region */
size_t GRegion::memory() const
{
	size_t mem = sizeof(GRegion);
	ArrayIterator<GZone*> iter(_zones);
	while (iter)
		mem += (iter++)->memory();
	return mem;
} // memory

/** addZone */
GZone* GRegion::addZone(const bool mode, const int size)
{
	GZone* z = new GZone(this);

	z->id(_zones.count());
	z->rpn(mode);
	if (size) z->resize(size);

	_zones.add(z);
	return z;
} // addZone

/** addZone */
GZone* GRegion::addZone(GZone* z)
{
	GZone *nz = addZone(z->rpn(), z->size());
	for (int i=0; i<z->size(); i++) {
		GBody *body = (*z)[i];
		nz->add(body->name(), body);
	}
	return nz;
} // addZone

/** delZone */
void GRegion::delZone(int i)
{
	// FIXME correct the id() of the zone!
	GZone *z = _zones[i];
	delete z;
	_zones.erase(i);
} // delZone

/** add2exp */
bool GRegion::add2exp(const char *token, GBody *body)
{
	assert(_zones.count() > 0);
	//nextGeneration();
	return _zones.tail()->add(token, body);
} // add2exp

/** @return hash value of expression */
dword GRegion::hash() const
{
	if (_hash>0) return _hash;
	_hash = 0;
	ArrayIterator<GZone*> zones_iter(_zones);
	while(zones_iter) {
		GZone *gzone = zones_iter++;
		_hash += zones_iter.index() * gzone->hash();
	}
	return _hash;
} // hash

/** clear */
void GRegion::clear()
{
	_hash = 0;
	// Clean the zones by popping the elements in order
	for (int i=0; i<_zones.count(); i++)
		delete _zones[i];
	_zones.clear();
	clearOBB();
} // clear

/** clearOBB */
void GRegion::clearOBB()
{
	if (cached_obbox) {
		delete cached_obbox;
		cached_obbox = NULL;
	}
} // clearOBB

/** 2D bounding box on viewer plane
 * @return bounding box of region
 * @WARNING z coordinates are meaning less
 */
BBox GRegion::bbox2D() const
{
	BBox bb;
//	assert(0);
	return bb;
} // bbox2D

/** guess of 3D bounding box
 * @return bounding box of zone
 */
BBox GRegion::bbox() const
{
	BBox bb;

	ArrayIterator<GZone*> iter(_zones);
	while (iter) {
		GZone* z = iter++;
		if (bb.isValid())
			bb |= z->bbox();
		else
			bb = z->bbox();
	}
	return bb;
} // bbox

/** guess of 3D oriented bounding box
 * @return bounding box of zone
 */
OBBox* GRegion::updateOBB() const
{
	OBBox* _obbox = new OBBox();
	_obbox->reset();

	ArrayIterator<GZone*> iter(_zones);
	while (iter) {
		GZone* z = iter++;
		if (_obbox->isValid())
			_obbox->Union(*z->obbox());
		else
			*_obbox = *z->obbox();
	}
	return _obbox;
} // obbox

/** inside
 * @param x,y,z		location of point to search
 * @param dx,dy,dz	direction in case on boundary
 * @return zone		return zone that location is inside otherwise NULL
 *
 * @warning	NON-Caching function
 */
GZone *GRegion::inside(const double  x, const double  y, const double  z,
		      const double dx, const double dy, const double dz) const
{
	ArrayIterator<GZone*> iter(_zones);
	while (iter) {
		GZone* zzz = iter++;
		if (zzz->inside(x,y,z,dx,dy,dz)) return zzz;
	}
	return NULL;
} // inside

/** operator << */
ostream& operator << (ostream& s, const GRegion& region)
{
	s << region.name() << endl;

	for (int i=0; i<region.nzones(); i++)
		s << "    " << *region.zone(i) << endl;
	return s;
} /* operator << */
