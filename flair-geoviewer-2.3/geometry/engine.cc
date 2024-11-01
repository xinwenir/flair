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

#include <math.h>
#include <time.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include <iomanip>
#include <iostream>

#include "point.h"
#include "engine.h"
#include "kernel.h"
#include "viewport.h"

using namespace std;

/* ========================== GeometryEngine ========================== */
/** GeometryEngine */
GeometryEngine::GeometryEngine(const Geometry* g, GeometryKernel* k) :
	_geometry(g),
	_kernel(k),
	bodies(32),
	zonesSorted(32),
	projectBodies(4),
	irurAccel(32)
{
	gBodyCheckId    = 0;
	gBodyMaxCheckId = 0;
//	nProjectBodies  = 0;
	nClipBodies     = 0;

	// Accelerators
	irurAccel.compare(CBodyOrderAccel::compare);
	lastIRURCheckID = -1;

	pthread_mutex_init(&mutexZones, NULL);
} // GeometryEngine

/** ~GeometryEngine */
GeometryEngine::~GeometryEngine()
{
	bodies.clear();
	zonesSorted.clear();
	pthread_mutex_destroy(&mutexZones);
} /* ~GeometryEngine */

/** view */
ViewPort& GeometryEngine::view() const
{
	return _kernel->view;
} // view

/** getRegion */
inline VRegion* GeometryEngine::getRegion(const int id)
{
	return _kernel->getRegion(id);
} // getRegion

/** cleanup
 * Do not delete the data, simply destroy the links V[Body/Region] -> G[Body/Region]
 */
void GeometryEngine::cleanup()
{
	for (int i=0; i<bodies.size(); i++)
		bodies[i].init(NULL);
} // cleanup

#ifdef EXPERIMENTAL
/** clearProximity
 * increase checkid of proximity bodies and clear the proximity list
 */
void GeometryEngine::clearProximity()
{
	for (int i=0; i<proximity.size(); i++) {
		DUMP(cout << "clearProximity(" << proximity[i]->name() << ")" << endl);
		proximity[i]->resetCheckId();
	}
	proximity.clear();
} // clearProximity
#endif

// Convert the Bodies,Zones,Regions to Array(VBody,VZone,VRegion)
// do not delete anything but reassign
//
/** derive - derive the VBodies, VZones and VRegions from Geometry
 * It initializes the list of bodies, zones and regions in the Viewer,
 * note that it deletes all previous objects and creates new ones
 */
void GeometryEngine::derive()
{
	// WARNING: bodies and regions vectors can contain more or equal
	// number of VBodies vs _geometry->bodies/regions
	//
	// Synchronize the bodies
	bodies.clear();
	bodies.allocate(_geometry->bodies.count());
	for (int i=0; i<_geometry->bodies.count(); i++) {
		GBody* gbody = _geometry->bodies[i];
		VBody* vbody = _kernel->bodies[i];
		assert(gbody == vbody->body());
		bodies[i].init(gbody, vbody);
		bodies[i].checkId(&gBodyCheckId);
	}

	// Rebuild zones caching array
	zonesSorted.clear();
	for (int i=0; i<_geometry->regions.count(); i++) {
		VRegion* region = _kernel->getRegion(i);
		for (int j=0; j<region->zones().size(); j++) {
			VZone* zone = region->zones()[j];
			zonesSorted.add(zone);
			fillBodyZoneReferences(zone);
		}
	}
} // derive

/** fillBodyZoneReferences */
void GeometryEngine::fillBodyZoneReferences(VZone *zone)
{
	// Fill Body/Zone references
	for (int i=0; i<zone->size(); i++) {
		const GBody *gbody = zone->gexpr(i);
		if (!gbody->isOperator())
			getBody(gbody)->addZone(zone);
	}
} // fillBodyZoneReferences

/** find zone where point is
 * WARNING this routine is optimized for points close to the 2D viewport ONLY!
 * @param x,y,z		location to search
 * @param dx,dy,dz	direction in case we are on the boundary
 * @param zone          candidate zone to be checked first
 * @return first zone that belongs to x,y,z with direction (dx,dy,dz)
 */
VZone* GeometryEngine::where2D(const double  x, const double  y, const double  z,
			     const double dx, const double dy, const double dz,
			     VZone *zone)
{
	if (zone)
		if (zone->inside2D(this, x,y,z, dx,dy,dz)) return zone;

	for (int i=_geometry->regions.count()-1; i>=0; i--) {
		VZone *vzone = getRegion(i)->inside2D(this, x,y,z, dx,dy,dz);
		if (vzone) return vzone;
	}
	return NULL;
} // where2D

/** find zone where a 3D point is
 * @param x,y,z		location to search
 * @param dx,dy,dz	direction in case we are on the boundary
 * @param zone          candidate zone to be checked first
 * @return first zone that belongs to x,y,z with direction (dx,dy,dz)
 */
VZone* GeometryEngine::where(const double  x, const double  y, const double  z,
			     const double dx, const double dy, const double dz,
			     VZone *zone)
{
	if (zone)
		if (zone->inside(this, x,y,z, dx,dy,dz)) return zone;

	for (int i=_geometry->regions.count()-1; i>=0; i--) {
		VZone *vzone = getRegion(i)->inside(this, x,y,z, dx,dy,dz);
		if (vzone) return vzone;
	}
	return NULL;
} // where

/** find zone where point is
 * WARNING this routine is optimized for points close to the 2D viewport ONLY!
 * @param p		location structure to search
 * @return first zone that belongs to x,y,z with direction (dx,dy,dz)
 */
int GeometryEngine::where2D(ZoneOfPoint *p)
{
	// FIXME use the "valid/invalid" information
	// if body is not modified but from the previous check
	// then do not check it's zones...
	// WARNING: e.g the second body can be modified, by the check below
	//              then is important to keep it!!!
	p->init();
#ifdef EXPERIMENTAL
	clearProximity();
#endif
	if (p->nOld != 0) {
		// First check all zones referring to that body
		for (int i=p->from; i<=p->to; i++) {
			VBody *body = (*p->V)[i].body;
			if (body==NULL) continue;
			CBody *cbody = getBody(body);
			for (int j=0; j<cbody->zones.size() && p->more(); j++) {
				STATS(stats.where_body_count++);
				VZone *zone = cbody->zones[j];
				if (zone->inside2D(this, p->x,p->y,p->z, p->dx,p->dy,p->dz)) {
					STATS(stats.where_body_found++);
					p->add(zone);
				}
			}
		}

		// loop over body zones
		for (int i=0; i<getBody(p->body)->zones.size() && p->more(); i++) {
			STATS(stats.where_body_count++);
			VZone *zone = getBody(p->body)->zones[i];
			if (zone->inside2D(this, p->x,p->y,p->z, p->dx,p->dy,p->dz)) {
				STATS(stats.where_body_found++);
				p->add(zone);
			}
		}

		// Then all zones in zoneOld
		for (int i=0; i<p->nOld && p->more(); i++) {
			STATS(stats.where_zoneold_count++);
			if (p->zoneOld[i]->inside2D(this, p->x,p->y,p->z, p->dx,p->dy,p->dz)) {
				STATS(stats.where_zoneold_found++);
				p->add(p->zoneOld[i]);
			}
		}
	}

	// Nothing found check all the remaining zones
	if (p->n==0) {
		STATS(stats.where_zone_none++);
		// We should check for touching Conics and
		// scan all regions sharing a conic...
		// otherwise it wont find possible errors
		// with touching surfaces
		// Check everything...
		//
		// Check all zones sorted
		for (int i=zonesSorted.size()-1; i >= 0 && p->more(); i--) {
			pthread_mutex_lock(&mutexZones);
			VZone *zone = zonesSorted[i];
			pthread_mutex_unlock(&mutexZones);
			if (zone->inside2D(this, p->x,p->y,p->z, p->dx,p->dy,p->dz)) {
				// promote zone
				if (i<=zonesSorted.size()-1) {
					pthread_mutex_lock(&mutexZones);
					zonesSorted.erase(i);
					zonesSorted.add(zone);
					pthread_mutex_unlock(&mutexZones);
				}
				p->add(zone);
			}
		}
	}

	return p->n;
} // where2D

/** find zone where point is in the transformed system
 * @param x,y,z		location to search
 * @param dx,dy,dz	direction in case we are on the boundary
 * @param t		distance to check
 * @param zone		current zone to check
 * @param body		check zones linked to body
 * @return first region that belongs to (x+dx*t,y+dy*t,z+dz*t)
 */
VZone* GeometryEngine::whereRay(const double  x, const double  y, const double  z,
				const double dx, const double dy, const double dz,
				const double t, VZone *zone, CBody *body)
{
	STATS(stats.whereRay_call_count++);

	// Check current zone if specified
	if (zone) {
		STATS(stats.whereRay_zoneCheck_count++);
		if (zone->insideRay(this, x,y,z, dx,dy,dz, t)) {
			STATS(stats.whereRay_zoneHit_count++);
			STATS(stats.whereRay_zoneIndex_sum += 1);
			return zone;
		}
	}

	// Check zones linked to body if any
	if (body) {
		for (int i=body->zones.size()-1; i>=0; i--) {
			VZone *zzz = body->zones[i];
#if 0
				if (!zzz->region()->obbox()->insideRay(x,y,z, dx,dy,dz, t)) continue;
				if (!zzz->obbox()->insideRay(x,y,z, dx,dy,dz, t)) continue;
#endif
			STATS(stats.whereRay_bodyCheck_count++);
			if (zzz->insideRay(this, x,y,z, dx,dy,dz, t)) {
				STATS(stats.whereRay_bodyHit_count++);
				STATS(stats.whereRay_zoneIndex_sum += (body->zones.size() - i));
				// promote zone
				body->zones.erase(i);
				body->zones.add(zzz);
				return zzz;
			}
		}
	}

	// Check everything...
	for (int i=zonesSorted.size()-1; i>=0; i--) {
		pthread_mutex_lock(&mutexZones);
		VZone *zzz = zonesSorted[i];
		pthread_mutex_unlock(&mutexZones);
#if 0
			if (!zzz->region()->obbox()->insideRay(x,y,z, dx,dy,dz, t)) continue;
			if (!zzz->obbox()->insideRay(x,y,z, dx,dy,dz, t)) continue;
#endif
		STATS(stats.whereRay_otherCheck_count++);
		if (zzz->insideRay(this, x,y,z, dx,dy,dz, t)) {
			STATS(stats.whereRay_otherHit_count++);
			//int idx = iter.index()+1;
			//STATS(stats.whereRay_zoneIndex_sum += (zones.size() - idx));
			// promote zone
			if (i<=zonesSorted.size()-1) {
				pthread_mutex_lock(&mutexZones);
				zonesSorted.erase(i);
				zonesSorted.add(zzz);
				pthread_mutex_unlock(&mutexZones);
			}
			return zzz;
		}
		// No need to check the zones
	}
	STATS(stats.whereRay_error_count++);
	return NULL;
} // whereRay

/** fillIRURegionAccel
 * Fill Accelerator structure for intersectRayUndefinedRegion
 * Initializes a list of bodies ordered by the distance of their intersections
 * with the current ray.
 * Checks if last checkID is still valid and if so skips all computations
 */
void GeometryEngine::fillIRURegionAccel(Ray *ray)
{
	if (gBodyCheckId == lastIRURCheckID) return; // Previous list still valid
	lastIRURCheckID = gBodyCheckId;

	irurAccel.clear();

	RaySegment& segment = ray->segment();
	double tmin = segment.tmin;
	double tmax = segment.tmax;

	for (int ib=0; ib<_geometry->bodies.count(); ib++) {
		CBody *body = getBody(ib);

		if (body->zones.empty()) continue;

		if (body->intersectRay(segment.pos, segment.dir)) {
			CBodyOrderAccel accel;
			if (InRangeOpen(tmin, body->tmin, tmax)) { // Add entry point
				accel.t = body->tmin;
				accel.body = body;
				irurAccel.add(accel);
			}
			if (InRangeOpen(tmin, body->tmax, tmax)) { // Add exit point
				accel.t = body->tmax;
				accel.body = body;
				irurAccel.add(accel);
			}
		}
	}
} // fillIRURegionAccel

/** intersectRayUndefinedRegion
 * assume an empty region as transparent and check for next intersection
 * @param ray	ray structure filled
 * @return false if ray arrives to end without hitting anything
 *	   true  if ray hits any defined region
 *
 * NOTE: Similar logic has the method Zone::intersectRay
 */
bool GeometryEngine::intersectRayUndefinedRegion(Ray *ray)
{
	STATS(stats.intersectRayUndefined_call_count++);

	if (ray->clip && applyClip(ray)) {
		RaySegment& segment = ray->segment();

		if (ray->ended()) return false;
		if (ray->n==0) ray->prevZone(segment.zone);	// update previous/master zone
		if (segment.zone) return true;
	}

	RaySegment& segment = ray->segment();
	ray->error = true;

	fillIRURegionAccel(ray);
	ArrayIterator<CBodyOrderAccel> iter(irurAccel);
	while (iter) {
		CBodyOrderAccel &accel = iter++;
		if (accel.t < segment.tmin) continue;
		segment.tmin = accel.t;
		segment.body = accel.body;

		for (int iz=segment.body->zones.size()-1; iz>=0; iz--) {
//		for (int iz=0; iz<segment.body->zones.size(); iz++) {
			VZone *zone = segment.body->zones[iz];
			double t = segment.tmin*(1.0+SMALL3D1);
			STATS(stats.whereRay_zoneIndex_sum++);
			if (zone->insideRay(this, segment.pos, segment.dir, t)) {
				STATS(stats.whereRay_errorHit_count++);
				// Lattices are trickier than this
				//if (zone->insideRay(ray->segment(0).pos, ray->segment(0).dir, ray->getT())) {
				segment.tmin = t;
				segment.zone = zone;
				return true;
			}
		}

	}
	return false;
} // intersectRayUndefinedRegion

/** initClip */
void GeometryEngine::initClip(Ray* ray)
{
	// First find if there is a clipping body in sight...
	const RaySegment& segment = ray->segment(0);
	ray->clip     = false;
	ray->clip_hit = false;

	for (int i=0; i<nClipBodies; i++) {
		CBody* cbody = clipBodies[i];
		if (!cbody) continue;

		bool intersect = cbody->intersectRay(segment.pos, segment.dir);
		if (intersect || clipNegative[i]) {
			bool tinv;
			if (clipNegative[i])
				tinv = !cbody->tinverse;
			else
				tinv = cbody->tinverse;

			double a = Max(cbody->tmin, segment.tmin);
			double b = Min(cbody->tmax, segment.tmax);
			if (!tinv)
				ray->clip |= a<b;
			else {
				ray->clip |= segment.tmin < Min(cbody->tmin, segment.tmax);
				ray->clip |= Max(cbody->tmax, segment.tmin) < segment.tmax;
			}
		}
	}
} // initClip

/** @return true if the clip body was applied */
inline bool GeometryEngine::applyClip(Ray *ray)
{
	assert(ray->clip);
	RaySegment &segment = ray->segment();
	double tpos = ray->T();
	ray->clip_hit = false;

	// loop until we are out of all clips
	bool more;
	do {
		more = false;
		for (int i=0; i<nClipBodies; i++) {
			CBody* cbody = clipBodies[i];
			if (!cbody) break;
			bool tinv;

			if (clipNegative[i])
				tinv = !cbody->tinverse;
			else
				tinv =  cbody->tinverse;

			if (!tinv) {	// normal clip
				if (InRange(cbody->tmin, tpos, cbody->tmax)) {
					// Hit inside the clip body, move to clip end and continue
					segment.tmin = (cbody->tmax - ray->Tsum())*(1.0+SMALL3D1);
					if (segment.tmax <= segment.tmin) {
						segment.tmin = segment.tmax;
						ray->clip = ray->clip_hit = false;
						return false;
					}
					ray->clip_hit = more = true;
					tpos = ray->T();
					segment.zone = whereRay(segment.pos,
							segment.dir,
							segment.tmin);
					if (segment.zone) segment.body = cbody;
				}
			} else	// clip with concave or outside of body
			if (tpos < cbody->tmin) {
				// Hit inside cbody (first segment),
				// move to boundary and continue
				segment.tmin = (cbody->tmin - ray->Tsum())*(1.0+SMALL3D1);
				if (segment.tmax <= segment.tmin) {
					segment.tmin = segment.tmax;
					ray->clip = ray->clip_hit = false;
					return false;
				}
				ray->clip_hit = more = true;
				tpos = ray->T();
				segment.zone = whereRay(segment.pos,
							segment.dir,
							segment.tmin);
				if (segment.zone)
					segment.body = cbody;
			} else
			if (tpos >= cbody->tmax) {
				// Hit inside cbody (second segment)
				// therefore it will never leave
				segment.tmin = INFINITE;
				ray->clip = ray->clip_hit = false;
				return false;
			}
		}
	} while (more && !segment.ended());

	return ray->clip_hit;
} // applyClip

/** initProject */
void GeometryEngine::initProject(Ray *ray)
{
	// First find if there is a projecting body in sight...
	const RaySegment& segment = ray->segment(0);
	ray->project     = false;
	ray->project_hit = false;

	for (int i=0; i<projectBodyCount(); i++) {
		CBody* pbody = projectBodies[i];
		if (!pbody) continue;

		if (pbody->intersectRay(segment.pos, segment.dir)) {
			ray->project |= InRange(segment.tmin, pbody->tmin, segment.tmax) ||
				        InRange(segment.tmin, pbody->tmax, segment.tmax);
		}
	}
} // initProject

/** applyProject */
inline bool GeometryEngine::applyProject(Ray *ray, const double prevT)
{
	assert(ray->project);
	RaySegment &segment = ray->segment();
	double T = ray->T();
	ray->project_hit = false;
	for (int i=0; i<projectBodyCount(); i++) {
		CBody* pbody = projectBodies[i];
		if (!pbody) continue;
		if (InRange(prevT, pbody->tmin, T)) {
			ray->project_hit = true;
			T = pbody->tmin;
			segment.tmin = pbody->tmin - ray->Tsum();
			segment.body = pbody;
		}
		if (InRange(prevT, pbody->tmax, T)) {
			ray->project_hit = true;
			T = pbody->tmax;
			segment.tmin = pbody->tmax - ray->Tsum();
			segment.body = pbody;
		}
	}
	return ray->project_hit;
} // applyProject

/** intersectLattice
 * @return	true  if lattice was applied
 *		false if an error occurred
 */
inline bool GeometryEngine::intersectLattice(Ray *ray)
{
	RaySegment& segment = ray->segment();
	VZone *zone = segment.zone;

	VRegion *region = zone->region();
	assert(region->region()->hasMatrix());

	//	RaySegment new_segment;
	//	new_segment.rotdefi = region->region()->rotdefi;
	RaySegment new_segment(region->region());

	// Find exit distance location from lattice
	double texit = segment.tmin;

	do {
		// find exit from zone
		CBody* body = zone->intersectRay(this,
					  segment.pos,
					  segment.dir,
					 &texit,
					  segment.tmax);

		// end of the ray?
		if (body==NULL) break;

		// check if we exited the region
		zone = region->insideRay(this,
					 segment.pos,
					 segment.dir,
					 texit);
	} while (zone!=NULL);

	// find new limits in the lattice
	// from 0 to Min(texit,tmax) - tmin
	// with a small shift
	new_segment.tmin = SMALL3D2;
	new_segment.tmax = Min(texit, segment.tmax) - segment.tmin;

	// Transform ray at the location that enters the lattice
	// to the prototype and come back to screen coordinates
	// Transformation: M^(-1) * R * M
	// M = screen transformation
	// R = rotdefi transformation
	//	const Matrix4& rot = _geometry->rotdefi(new_segment.rotdefi);
	//	new_segment.pos = rot * (Point)(segment.pos
	//				+ segment.tmin * segment.dir);
	new_segment.pos = new_segment.matrix() * (Point)(segment.pos
				+ segment.tmin * segment.dir);
//	new_segment.acc = 8.0*SMALL3D * (Abs(segment.pos.x) + Abs(segment.pos.y) + Abs(segment.pos.z) +
	//				Abs(rot(0,3)) + Abs(rot(1,3)) + Abs(rot(2,3)));
	new_segment.acc = 8.0*SMALL3D * ( Abs(segment.pos.x)
					+ Abs(segment.pos.y)
					+ Abs(segment.pos.z)
					+ Abs(new_segment.matrix()(0,3))
					+ Abs(new_segment.matrix()(1,3))
					+ Abs(new_segment.matrix()(2,3)));
//	new_segment.dir = rot.multVector(segment.dir);
//	new_segment.dir = rot * segment.dir;
	new_segment.dir = new_segment.matrix() * segment.dir;

	// Find the new zone in the prototype
	segment.bodyCheckId = incBodyCheckId();
	new_segment.zone = whereRay(new_segment.pos,
				    new_segment.dir,
				    new_segment.tmin);

	if (new_segment.zone == NULL) {
		// error in geometry, move to texit
		segment.zone = NULL;
		return false;
	}

	if (ray->push(new_segment)) { // segment is no longer the last segment
		_kernel->error("GeometryEngine::intersectLattice n>=MAXLEVEL. Possibly LATTICE infinite loop");
#if _DEBUG>1
		cerr << "ERROR: GeometryEngine::intersectLattice n>=MAXLEVEL. Possibly LATTICE infinite loop" << endl;
#endif
		bodyCheckId(segment.bodyCheckId);
		segment.zone = NULL;
		return false;
	}
	return true;
} // intersectLattice

/** intersectVoxel
 * @return	 1	hit found or ended
 *		-1	skip this zone
 *		-2	exit from voxel
 */
inline int GeometryEngine::intersectVoxel(Ray* ray)
{
	// Normal tracking, find next intersection
	RaySegment& segment = ray->segment();
	VZone* zone = segment.zone;

	double tentry = segment.tmin;
VOX:
	double tmin = segment.tmin;
	double T = ray->T();
	segment.tmin -= segment.tmin * SMALL3D;	// move backwards a bit

	if (_kernel->voxel.intersectRay(
			segment.pos,
			segment.dir,
			&segment.tmin)) {

		// if previous hit was clip
		if (ray->clip_hit && segment.tmin - tmin < SMALL3D4) return 1;

		// Apply clip
		if (ray->clip && applyClip(ray)) {
			if (ray->project && applyProject(ray, T)) return 1;
			if (ray->ended()) return 1;
			if (segment.zone != zone) return -1;
			// continue to check if we hit on a filled voxel
			goto VOX;
		}

		if (ray->project && applyProject(ray, T)) return 1;
		if (segment.tmin - tmin > SMALL3D4) segment.body = NULL;
		// either we hit the clipbody or the voxel, update master zone
		if (ray->n==0) ray->prevZone(segment.zone);

		// we are outside the region
		ray->moveby(SMALL3D);	// move a bit to ensure we are in the new voxel
		return 1;
	}

	// restore tmin and exit voxel zone
	segment.tmin = tentry;
	return -2;
} // intersectVoxel

/** intersectRay
 * find the next intersection of the ray with any zone and return.
 * It takes care of multiple-lattices as well as voxels
 * @param ray	structure contains the ray/depth
 * @param step	if true, the function returns after every step otherwise only
 *		when hitting a non-transparent object
 * @return true	 if we arrived to the end of the ray or error (and
 *		 ray->segment().zone will be NULL)
 *	   false the ray hit a solid object or next step was computed and there
 *		 is still way to go (if step=true)
 */
bool GeometryEngine::intersectRay(Ray *ray, bool step)
{
	STATS(stats.intersectRay_call_count++);

	if (!step) {
		if (ray->use_clip    && nClipBodies)        initClip(ray);
		if (ray->use_project && projectBodyCount()) initProject(ray);
	}

	// Loop until we exit the region (zone->region())
	while (true) {
		// Check for ending condition
		// XXX it doesn't work for negative rays!!!
		RaySegment& segment = ray->segment();

		// Segment has ended
		if (segment.ended()) {
			// Initial segment and ended, we arrived to the end
			if (ray->n == 0) return true;

			ray->pop();
			ray->segment().tmin += segment.tmax;

			RaySegment &s = ray->segment();
			//			if (segment.rotdefi) {
				// find new zone
			if (segment.region) {
				bodyCheckId(s.bodyCheckId);
				s.zone = whereRay(s.pos,
						  s.dir,
						  s.tmin,
						  s.zone,
						  NULL);
			}

			if (step) return false;
			continue;
		}

		VZone *zone = segment.zone;

		// Check for error
		if (zone == NULL) {
			if (!intersectRayUndefinedRegion(ray))
				return true;
			else
				continue;
		}

		// Check for clipping body
		if (!step) {
			if (zone->region()->type() == REGION_BLACKHOLE &&
			    (!ray->skip_1stblack || ray->prevZone() != NULL))
				return false;

			// Inside a solid object
			// TODO: Review this condition
			if (!ray->skip_current && zone->region()->type() == REGION_NORMAL &&
					(!ray->skip_transparent || !zone->region()->transparent() )) {

				double T = ray->T();
				// Check for clipping body
				if (ray->clip && applyClip(ray)) {
					if (ray->project && applyProject(ray, T))
						return ray->ended();
					if (ray->ended()) return true;
					if (ray->n==0) ray->prevZone(zone);	// update master zone
					continue;
				} else {
					return false;
				}
			}
			ray->skip_current = false;
		}

		// Only increment stats when considering a new intersection or segment
		STATS(stats.intersectRay_retry_count++);

		// ------------------- LATTICE -------------------
		// Could be that a region is marked as LATTICE but there is
		// no transformation associated with it
		//		if (zone->region()->rotdefi()) {
		if (zone->gregion()->hasMatrix()) {
			if (intersectLattice(ray) && step) return false;
			continue;
		} else	// end of LATTICE

		// ------------------- VOXEL -------------------
		if (zone->region()->type() == REGION_VOXEL) {
			int ret = intersectVoxel(ray);
			if (ret == 1) return segment.ended(); // Intersection in voxel
			if (ret == -1) continue; // Skip this zone
			// if (ret == -2) let it continue...
		} // End of VOXEL

		// ------------------- NORMAL tracking -------------------
		if (ray->n==0) ray->prevZone(segment.zone);	// update master zone
		// Move to next intersection
		//
		double T = ray->T();

		segment.body = segment.zone->intersectRay(this,
							  segment.pos,
							  segment.dir,
							 &segment.tmin,
							  segment.tmax);
		if (ray->project && applyProject(ray, T)) return ray->ended();

		if (segment.ended()) continue;
		if (segment.body==NULL) return false;

		// Check for clipping body
		if (ray->clip && applyClip(ray)) {
			if (ray->project && applyProject(ray, T)) return ray->ended();
			if (ray->ended()) return true;
			if (ray->n==0) ray->prevZone(segment.zone);	// update master zone
		}

		// Find the zone we are moving to
		segment.zone = whereRay(segment.pos,
					segment.dir,
					segment.tmin,
					NULL,
					segment.body);

		if (step) return (segment.zone==NULL);
	}
} // intersectRay

/** bodiesMemory */
size_t GeometryEngine::bodiesMemory() const
{
	size_t bmem = bodies.memory();
	for (int i=0; i<bodies.count(); i++)
		bmem += bodies[i].memory();
	return bmem;
} // bodiesMemory

/** memory */
size_t GeometryEngine::memory() const
{
	return sizeof(GeometryEngine) + bodiesMemory() + zonesSorted.memory() + irurAccel.memory();
} // memory

/** printMemory */
void GeometryEngine::printMemory() const
{
	cout << endl << "GeometryEngine:" << endl;
	cout << "Memory:"     << endl;
	cout << "\tSelf:\t"   << sizeof(GeometryEngine) << endl;
	cout << "\tBodies:\t" << bodiesMemory() << endl;
	cout << "\tZones:\t"  << zonesSorted.memory() << endl;
	cout << "\tTotal:\t"  << memory() << endl;
} // printMemory

/* ========================== GeometryStats ========================== */
#ifdef _STAT
/** reset drawing statistics */
void GeometryEngineStats::reset()
{
	whereRay_call_count              = 0;
	whereRay_zoneCheck_count         = 0;
	whereRay_zoneHit_count           = 0;
	whereRay_bodyCheck_count         = 0;
	whereRay_bodyHit_count           = 0;
	whereRay_otherCheck_count        = 0;
	whereRay_otherHit_count          = 0;
	whereRay_error_count             = 0;
	intersectRay_call_count          = 0;
	intersectRay_retry_count         = 0;
	intersectRayUndefined_call_count = 0;
	whereRay_zoneHit_count           = 0;
	whereRay_zoneIndex_sum           = 0;
	whereRay_errorHit_count          = 0;
	where_body_count		 = 0;
	where_body_found		 = 0;
	where_zoneold_count		 = 0;
	where_zoneold_found		 = 0;
	where_zone_none			 = 0;
	where_region_count		 = 0;
	where_region_found		 = 0;
	where_all_count			 = 0;
	where_all_notskip		 = 0;
	where_all_found			 = 0;
} // reset

/** operator << */
ostream & operator<<(ostream &os, const GeometryEngineStats &gv)
{
	if (gv.intersectRay_call_count) {
		os << "intersectRay call count: "          << gv.intersectRay_call_count << endl;
		os << "intersectRay retry count: "         << gv.intersectRay_retry_count;
		os << " (" << setprecision(5) << ((double)gv.intersectRay_retry_count/gv.intersectRay_call_count)*100.0 << "%)" << endl;

		os << "intersectRayUndefined call count: " << gv.intersectRayUndefined_call_count << endl;
		os << endl;
	}

	if (gv.whereRay_call_count) {
		os << "whereRay call count: "              << gv.whereRay_call_count << endl;
		os << "whereRay zone check count: "        << gv.whereRay_zoneCheck_count << endl;
		os << "whereRay zone hit count: "          << gv.whereRay_zoneHit_count;
		os << " (" << setprecision(5)  << ((double)gv.whereRay_zoneHit_count/gv.whereRay_call_count)*100.0       << "% of whereRay calls)" << endl;
		os << "whereRay body check count: "        << gv.whereRay_bodyCheck_count << endl;
		os << "whereRay body hit count: "          << gv.whereRay_bodyHit_count;
		os << " (" << setprecision(5) << ((double)gv.whereRay_bodyHit_count/gv.whereRay_call_count)*100.0       << "% of whereRay calls)" << endl;
		os << "whereRay other check count: "       << gv.whereRay_otherCheck_count << endl;
		os << "whereRay other hit count: "         << gv.whereRay_otherHit_count;
		os << " (" << setprecision(5) << ((double)gv.whereRay_otherHit_count/gv.whereRay_call_count)*100.0      << "% of whereRay calls)" << endl;
		os << "whereRay error count: "             << gv.whereRay_error_count;
		os << " (" << setprecision(5) << ((double)gv.whereRay_error_count/gv.whereRay_call_count)*100.0         << "% of whereRay calls)" << endl;
		os << "whereRay total zone count: " << gv.whereRay_zoneCheck_count + gv.whereRay_bodyCheck_count + gv.whereRay_otherCheck_count << endl;
		os << endl;
		os << "whereRay total zone count: " << gv.whereRay_zoneIndex_sum << endl;
		os << "whereRay zone index average: " << (double)gv.whereRay_zoneIndex_sum / (gv.whereRay_zoneHit_count + gv.whereRay_bodyHit_count + gv.whereRay_otherHit_count + gv.whereRay_errorHit_count) << endl;
		os << endl;
	}

	os << " where_body_count:      " << gv.where_body_count      << endl;
	os << " where_body_found:      " << gv.where_body_found      << endl;
	os << " where_zoneold_count:   " << gv.where_zoneold_count   << endl;
	os << " where_zoneold_found:   " << gv.where_zoneold_found   << endl;
	os << " where_zone_none:       " << gv.where_zone_none       << endl;
	os << " where_region_count:    " << gv.where_region_count    << endl;
	os << " where_region_found:    " << gv.where_region_found    << endl;
	os << " where_all_count:       " << gv.where_all_count       << endl;
	os << " where_all_notskip:     " << gv.where_all_notskip     << " skipped: "
					  << gv.where_all_count - gv.where_all_notskip<< endl;
	os << " where_all_found:       " << gv.where_all_found       << endl;

	return os;
} /* operator << */
#endif
