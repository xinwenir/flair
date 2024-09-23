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

#ifndef __CBODY_H
#define __CBODY_H

#include <iostream>

#include "array.h"
#include "gbody.h"
#include "vbody.h"

class VZone;
class GeometryEngine;

/* ============================== CBody =============================== */
/** CBody cached viewer body class
 * WARNING: Do not mix 3D and 2D operations at the same time
 */
class CBody {
private:
const	GBody*	gbody;	// reference to Gbody
const	VBody*	vbody;	// reference to Vbody

	bool	checkInside;	// Result of the last inside check
	int	_checkId;	// instance checked Id
	int*	gCheckId;	// global current check Id

public:
	bool	tinverse;	// if [tmin,tmax] solution should be
				// used in inverse e.g. for concave objects
				// or cones
	double	tmin, tmax;	// ray intersection min, max distances
	Array<VZone*>	zones;  // zones that refer to this body

#if _DEBUG>2
private:
	double	xx,yy,zz,ddx,ddy,ddz;	// STRONG debugging remember pos
#endif

public:
	CBody(GBody *gb=NULL, VBody *vb=NULL) : gCheckId(NULL) { init(gb, vb); }

	void	init(GBody *gb=NULL, VBody *vb=NULL) {
			gbody = gb;
			vbody = vb;
			_checkId = -1;
			zones.clear();
		}

const	char*	name()		const;
const	GBody*	body()		const { return gbody; }
	int	id()		const { return gbody->id(); }

	// Quadrics
	int	nC()		const { return vbody->nC; }

	Location location()	const { return vbody->location; }

	// Inside2D version before use please take care of incBodyCheckId2D()
	bool	inside2D(
#ifdef EXPERIMENTAL
			 GeometryEngine* engine,
#endif
			 const double x, const double y, const double z,
			 const double dx, const double dy, const double dz);

	bool	inside( const double x, const double y, const double z,
			const double dx, const double dy, const double dz) {
				if (*gCheckId != _checkId) {
					checkInside = gbody->inside(x,y,z, dx,dy,dz);
					checked();
				}
				return checkInside;
			}

	/**
	 * Find if a point (x,y,z)+t*(dx,dy,dz) is inside the body or not, using the
	 * information from the intersectRay
	 * Returns: true  inside
	 *          false outside
	 */
	bool	insideRay(const double x, const double y, const double z,
			  const double dx, const double dy, const double dz,
			  const double t) {
				if (*gCheckId != _checkId) {
					tinverse = gbody->intersectRay(x,y,z, dx,dy,dz, &tmin,&tmax);
					_checkId = *gCheckId;
#if _DEBUG>2
					xx = x;
					yy = y;
					zz = z;
					ddx = dx;
					ddy = dy;
					ddz = dz;
				} else {
					assert(x==xx);
					assert(y==yy);
					assert(z==zz);
					assert(dx==ddx);
					assert(dy==ddy);
					assert(dz==ddz);
#endif
				}

				if (tinverse)
					return (bool)(t<tmin || tmax<t);
				else
					return (bool)(tmin<=t && t<=tmax);
			}

	bool	insideRay(const Vector& p, const Vector &d, const double t)
			{ return insideRay(p.x,p.y,p.z, d.x,d.y,d.z, t); }

	/** intersectRay cached
	 * calculates the tmin, tmax and tinverse variables from the intersection
	 * of the ray with the body
	 * @param x,y,z		position vector
	 * @param dx,dy,dz	direction vector (has to be normalized)
	 */
	bool	intersectRay(const double x, const double y, const double z,
			     const double dx, const double dy, const double dz) {
				if (*gCheckId != _checkId) {
					tinverse  = gbody->intersectRay(x,y,z, dx,dy,dz, &tmin,&tmax);
					_checkId  = *gCheckId;
#if _DEBUG>2
					xx = x;
					yy = y;
					zz = z;
					ddx = dx;
					ddy = dy;
					ddz = dz;
				} else {
					assert(x==xx);
					assert(y==yy);
					assert(z==zz);
					assert(dx==ddx);
					assert(dy==ddy);
					assert(dz==ddz);
#endif
				}
				return (bool)(tmin<tmax);
			}
	bool	intersectRay(const Vector& p, const Vector &d)
			{ return intersectRay(p.x,p.y,p.z, d.x,d.y,d.z); }

	// Checkid
	void	checkId(int *i)		{ gCheckId = i; _checkId = -1; }
	void	resetCheckId()		{ _checkId  = -1; }
	void	checked()		{ _checkId  = *gCheckId; }
	bool	isChecked()		{ return *gCheckId == _checkId; }

	/** add zones that refers to this body
	 * @param zone	to be added to the list
	 */
	void	addZone(VZone *zone)	{ zones.add(zone); }
	size_t	memory()	const	{ return sizeof(CBody) + zones.memory(); }
}; // CBody
#endif
