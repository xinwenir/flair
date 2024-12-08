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

#ifndef __RAY_H
#define __RAY_H

#include <string.h>

#include "os.h"
#include "geo.h"
#include "point.h"
#include "vector.h"
#include "matrix4.h"

#include "cbody.h"
#include "vzone.h"

/* ============================= RaySegment =========================== */
class RaySegment {
public:
	Point	pos;
	Vector	dir;
	VZone	*zone;			// Current location
	CBody	*body;			// entering body
const	GRegion	*region;		// lattice region
	double	tmin;			// range of segment
	double	tmax;
	double	acc;
  //	int	rotdefi;		// associated rotdefi
	int	bodyCheckId;		// bodyCheckId, used for lattices

public:
  //	RaySegment() :  pos(0.,0.,0.),
	RaySegment(GRegion* reg=NULL) :
            pos(0.,0.,0.),
			dir(0.,0.,0.),
			zone(NULL),
			body(NULL),
			region(reg),
			tmin(SMALL3D),
			tmax(INFINITE),
			acc(SMALL3D),
			//			rotdefi(0),
			bodyCheckId(0) {}
	RaySegment(const double x, const double y, const double z,
		   const double dx, const double dy, const double dz,
		   VZone* zz) :
			pos(x,y,z),
			dir(dx,dy,dz),
			zone(zz),
			body(NULL),
			region(NULL),
			tmin(SMALL3D),
			tmax(INFINITE),
			acc(SMALL3D),
			//			rotdefi(0),
			bodyCheckId(0) {}
	RaySegment(const Point& p, const Vector& d, VZone* zz) :
			pos(p),
			dir(d),
			zone(zz),
			body(NULL),
			region(NULL),
			tmin(SMALL3D),
			tmax(INFINITE),
			acc(SMALL3D),
			//	rotdefi(0),
			bodyCheckId(0) {}

	bool	ended()	const	{ return tmin*(1.0+acc) >= tmax; }

	// Relative hit position inside the segment-frame
	// LATTICE: are screen coordinates transformed by lattice rot-translation
	// VOXEL: are absolute coordinates
	Point	hit()	const	{ return pos + tmin*dir; }

	Point	hit(const double boost) const
				{ return pos + tmin*(1.0+boost)*dir; }
const	Matrix4& matrix()	{ assert(region); return region->matrix(); }
const	Matrix4& invMatrix()	{ assert(region); return region->invMatrix(); }
}; // RaySegment

/* ================================= Ray ============================== */
/** Intersect Ray structure
 * This structure contains a list of RaySegments,
 *   the last segment of the list, the one pointed by n (segments[n]) is
 *   considered the working one and the rest in the list are the previous
 *   segments.
 */
class Ray {
public:
	bool	error;			// If ray passed through an undefined region
	bool	shadow;			// calculate shadows for ray
	int	lights;			// number of lights to use
	bool	skip_current;		// Skip current zone as if it was transparent
	bool	skip_transparent;	// Skip transparent zones
	bool	skip_1stblack;		// Skip the starting blackhole if any
	int	voxelreg;		// voxel region if any

	Vector	normal;			// keep normal at hit position

	/* Clipping body */
	bool	use_clip;		// flag if we have to use clipping
	bool	clip;			// possible clip along the ray
	bool	clip_hit;		// clipping body was hit

	/* Project body */
	bool	use_project;		// flag if we have to use projection
	bool	project;		// possible projection along the ray
	bool	project_hit;		// project was hit
	int	project_alpha;		// projection transparency

	int	depth;			// segment depth
	int	max_depth;		// maximum segment depth

private:
	VZone	*_prevzone;		// current upstream zone (NOT the hit zone)
	double	tsum;			// tmin running length of previous levels
public:
	int	n;			// lattice depth level
private:
	RaySegment segments[MAXLEVEL];

public:
	Ray()	: shadow(false),
		  lights(0),
		  skip_1stblack(false),
		  use_clip(false),
		  use_project(false),
		  depth(0),
		  max_depth(0)
			{init();}

	Ray(const Ray& ray) {
			shadow        = ray.shadow;
			lights        = ray.lights;
			use_clip      = ray.use_clip;
			use_project   = ray.use_project;
			skip_1stblack = ray.skip_1stblack;
			depth         = ray.depth+1;
			max_depth     = ray.max_depth;
			init();
		}

	void	init() {
			voxelreg      = -1;

			error         = false;
			skip_current  = false;
			skip_transparent = true;

			clip          = false;
			clip_hit      = false;

			project       = false;
			project_hit   = false;
			project_alpha = 0;

			_prevzone     = NULL;
			tsum          = 0.0;
			n             = -1;
		}

	// Return true if ray segments buffer is full
	bool	isFull()	const	{ return n+1 >= MAXLEVEL; }
	bool	isEmpty()	const	{ return n == -1; }
	bool	start()		const	{ assert(n>-1); return T() <= SMALL3D3; }

	void	incDepth()		{ depth++; }

	// Consolidate the current segment and make current a new one
	bool	push(const RaySegment& _segment) {
			if (isFull()) {
#if _DEBUG>1
				std::cerr << "Error pushing raySegment, ray is full" << std::endl;
#endif
				return true;
			}
			if (n >= 0) {
				tsum += segments[n].tmin;
				if (n == 0)
					_prevzone = segments[n].zone;
			}
			n++;
			segments[n] = _segment;
			return false;
		}

	// Discards the current segment and moves to current the previous one
	bool	pop() {
			if (isEmpty()) return true;
			n--;
			if (n >= 0) {
				tsum -= segments[n].tmin;
			}
			return false;
		}

	// Return reference to any segment
	RaySegment&	segment(int idx)	{ assert(idx<=n); return segments[idx]; }
const	RaySegment&	segment(int idx) const	{ assert(idx<=n); return segments[idx]; }
	// Current segment
	RaySegment&	segment()		{ return segments[n]; }
const	RaySegment&	segment()	const	{ return segments[n]; }

const	Point&	pos()		const { return segments[0].pos; }
const	Vector&	dir()		const { return segments[0].dir; }

	// Absolute screen coordinates hit position
	Point	hit()		const { return pos() + (segments[n].tmin+tsum)*dir(); }
	Point	hit(const double boost) const {
			assert(boost!=0.0);
			if (boost>0.0)	// Forward boost
				return pos() + ((segments[n].tmin+tsum)*(1.0+boost))*dir();
			else		// Backward boost
				return pos() + ((segments[n].tmin+tsum)/(1.0-boost))*dir();
		}

	// set/get
	void	prevZone(VZone *z)    { _prevzone = z; }
	VZone*	prevZone()	const { return _prevzone; }
	CBody*	hitBody()	const { return segment().body; }
	VZone*	hitZone()	const { return segment().zone; }
	VRegion* hitRegion()	const { return segment().zone->region(); }

	void	moveby(double s)      { segment().tmin = T()*(1.0+s) - tsum; }

	double	T()		const { return segments[n].tmin + tsum; }
	double	Tsum()		const { return tsum; }
	bool	ended()		const { return T()>=segments[0].tmax; }
}; // Ray

#endif
