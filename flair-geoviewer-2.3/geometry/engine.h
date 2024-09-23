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
 */

#ifndef __GEOMETRY_ENGINE_H
#define __GEOMETRY_ENGINE_H

#include <time.h>
#include <iosfwd>
#include <stdio.h>
#include <ostream>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <pthread.h>

#include "os.h"
#include "geo.h"
#include "array.h"
#include "vector.h"
#include "random.h"
#include "matrix4.h"

#include "ray.h"
#include "cbody.h"
#include "vzone.h"
#include "vregion.h"

#include "geometry.h"

class GeometryKernel;

#define MAXZONE		 5
#define MAXCLIP		 5

#ifdef _STAT
/* =========================== GeometryEngineStats ========================= */
/** Geometry Engine Statistics */
class GeometryEngineStats {
public:
	long long whereRay_call_count;
	long long whereRay_zoneCheck_count;
	long long whereRay_zoneHit_count;
	long long whereRay_bodyCheck_count;
	long long whereRay_bodyHit_count;
	long long whereRay_otherCheck_count;
	long long whereRay_otherHit_count;
	long long whereRay_error_count;
	long long whereRay_errorHit_count;
	long long intersectRay_call_count;
	long long intersectRay_retry_count;
	long long intersectRayUndefined_call_count;

	long long whereRay_zoneIndex_sum;

	long long where_body_count;
	long long where_body_found;
	long long where_zoneold_count;
	long long where_zoneold_found;
	long long where_zone_none;
	long long where_region_count;
	long long where_region_found;
	long long where_all_count;
	long long where_all_notskip;
	long long where_all_found;

	GeometryEngineStats() {reset();}

	void reset();
}; // GeometryStats
std::ostream& operator << (std::ostream&, const GeometryEngineStats&);
#endif

/* ============================= ZoneOfPoint ========================== */
// Scanning of region information
class ZoneOfPoint {
public:
	int	 nMax;			// Maximum number of regions to fill
	int	 n, nOld;		// New/Old zones found
	double	 x,  y,  z;		// Position		// FIXME convert to Vector()
	double	 dx, dy, dz;		// Direction		// FIXME convert to Vector()
	VBody	*body;
	Array<Vertex2D>	*V;		// Vertices array of body where point belongs to
	int	from, to;		// list of vertices/bodies V[fromId:toId] that
					// are involved in the search
	VZone	*zone[MAXZONE];		// new zones (only first zone from each region)
	VZone	*zoneOld[MAXZONE];	// old zones (only first zone from each region)
	int	 nRegions;		// New/Old regions found

public:
	ZoneOfPoint()		: nMax(MAXZONE), n(0), nOld(0) {}
	void zero() {
			n = nOld  = 0;
			from = to = 0;
			nRegions  = 0;
		}
	void copy2Old() {
			// copy previous search
			if (n>0) {
				memcpy(zoneOld, zone, n*sizeof(VZone*));
			}
			nOld = n;
			n = 0;
			nRegions = 0;
		}

	void init()		{ n = nRegions = 0; }
	bool more()	const	{ return n < nMax; }
	void add(VZone *zz) {
			for (int i=0; i<n; i++)
				if (zone[i]->region() == zz->region()) return;
			nRegions++;
			zone[n++] = zz;
			assert(n<=nMax);
		}

	bool hasZone(const VZone *zo) const {
			for (int i=0; i<n; i++)
				if (zone[i] == zo) return true;
			return false;
		}
}; // ZoneOfPoint

/* ========================== CBodyOrderAccel ========================= */
/** CBodyOrderAccel
 * Used by intersectRayUndefinedRegion to maintain a list of bodies ordered by its
 * tmin and tmax distance to the ray origin
 */
class CBodyOrderAccel {
public:
	CBody *body;
	double t;

	CBodyOrderAccel() : body(NULL), t(-INFINITE) {}

/** Compare distance to ray origin */
static  int compare(const CBodyOrderAccel &a, const CBodyOrderAccel &b) { return Cmp(a.t,b.t); }
}; // CBodyOrderAccel

/* ========================== GeometryEngine ========================== */
/** Geometry Engine class */
class GeometryEngine {
protected:
const	Geometry*	_geometry;	/** father geometry object	*/
	GeometryKernel*	_kernel;	/** geometry kernel object	*/

public:
	Array<CBody>	bodies;		/** body array			*/
	Array<VZone*>	zonesSorted;	/** sorted list of zones	*/
#ifdef EXPERIMENTAL
	Array<CBody*>	proximity;	/** list of bodies that last inside2D check */
					/** fell into proximity limits of accuracy  */
#endif
	Random	random;			/** random number generator	*/

protected:
	int	nClipBodies;		/** number of clip bodies	*/
	CBody*	clipBodies[MAXCLIP];	/** clipping bodies		*/
	bool	clipNegative[MAXCLIP];	/** negative side of clip body	*/
	Array<CBody*>	projectBodies;	/** projection bodies		*/

	int	lastIRURCheckID;
	Array<CBodyOrderAccel>	 irurAccel;	/** intersectRayUndefinedRegion accelerator structure */

	int	gBodyCheckId;		/** checkid for bodies		*/
	int	gBodyMaxCheckId;	/** checkid for bodies maximum previous value*/
	pthread_mutex_t	mutexZones;	/** zones Sorted mutex		*/

public:
#ifdef _STAT
	GeometryEngineStats stats;      /** Statistics                  */
#endif

public:
	GeometryEngine(const Geometry* g=NULL, GeometryKernel* k=NULL);
	~GeometryEngine();

	void	cleanup();

	// get/set
	void	geometry(const Geometry* g)	{ _geometry = g; }
	void	kernel(GeometryKernel* v)	{ _kernel = v; }
	ViewPort& view() const;

	/* bodies */
	// WARNING unprotected for limits. Check before
	CBody*	getBody(const int id)		{ return &bodies[id]; }
	CBody*	getBody(const GBody* body)	{ return &bodies[body->id()]; }
	CBody*	getBody(const VBody* body)	{ return &bodies[body->id()]; }
	CBody*	getBody(const char *name) {
			const GBody *body = _geometry->getBody(name);
			return body?getBody(body):NULL;
		}

	/* regions */
inline	VRegion* getRegion(const int id);

	/* clip body */
	void	clipBodyClear()			{ nClipBodies = 0; }
	void	clipBodyAdd(int id) {
			if (id>=0 && id<bodies.size() && nClipBodies<MAXCLIP-1) {
				clipBodies[nClipBodies] = getBody(id);
				clipNegative[nClipBodies] = false;
				nClipBodies++;
			}
		}
	void	clipBodyNegative(int id) {
			if (id>=0 && id<nClipBodies)
				clipNegative[id] = true;
		}
	int	clipBody(int id)	const	{
			return id<nClipBodies? clipBodies[id]->id() : -1;
		}
	int	clipBodyCount()		const	{ return nClipBodies; }

	/* project body */
	void	projectBodyClear()			{ projectBodies.clear(); }
	void	projectBodyAdd(int id) {
			if (id>=0 && id<bodies.size())
				projectBodies.append(getBody(id));
		}
	int	projectBody(int id)	const	{
			return id<projectBodies.count()? projectBodies[id]->id() : -1;
		}
	int	projectBodyCount()	const	{ return projectBodies.count(); }

#ifdef EXPERIMENTAL
	/** proximity */
	void	addProximity(CBody *body)	{ proximity.add(body); }
	void	clearProximity();
#endif

	/** Derive */
	void	 derive();

	/* 2D scan functions */

	/* where */
	VZone*	where2D(const double  x, const double  y, const double  z,
			const double dx, const double dy, const double dz,
			VZone *zone=NULL);
	VZone*	where2D(const Vector &pos, const Vector &dir, VZone *zone=NULL)
			{ return where2D(pos.x, pos.y, pos.z, dir.x, dir.y, dir.z, zone); };
	int	where2D(ZoneOfPoint *p);

	VZone*	where(const double  x, const double  y, const double  z,
			const double dx, const double dy, const double dz,
			VZone *zone=NULL);
	VZone*	where(const Vector &pos, const Vector &dir, VZone *zone=NULL)
			{ return where(pos.x, pos.y, pos.z, dir.x, dir.y, dir.z, zone); };

	/* ray methods */
	VZone*	whereRay(const double  x, const double  y, const double  z,
			 const double dx, const double dy, const double dz,
			 const double t,
			 VZone *zone=NULL, CBody *body=NULL);
	VZone*	whereRay(const Vector& p, const Vector &d, const double t,
			VZone *zone=NULL, CBody *body=NULL)
			{ return whereRay(p.x,p.y,p.z, d.x,d.y,d.z, t, zone, body); }

	bool	intersectRay(Ray *ray, bool step);
	bool	intersectRayUndefinedRegion(Ray *ray);

	/** Increase body check id. Needed when we want to calculate for the bodies
	 *  a new location
	 */
	int	incBodyCheckId()	{
			int old = gBodyCheckId;
			gBodyCheckId = ++gBodyMaxCheckId;
#ifdef EXPERIMENTAL
			proximity.clear();
#endif
			return old;
		}
	void	bodyCheckId(int value)	{ gBodyCheckId = value; }
	int	bodyCheckId()	const	{ return gBodyCheckId; }

	/* Info */
	size_t	bodiesMemory()	const;
	size_t	memory()	const;
	void	printMemory()	const;

private:
	void    fillBodyZoneReferences(VZone *zone);

	void	initClip(Ray *ray);
	void	initProject(Ray *ray);

	bool	applyClip(Ray *ray);
	bool	applyProject(Ray *ray, const double prevT);
	bool	intersectLattice(Ray *ray);
	int	intersectVoxel(Ray *ray);

	void	fillIRURegionAccel(Ray *ray);

friend class VZone;
friend class VRegion;
}; // GeometryEngine

#endif
