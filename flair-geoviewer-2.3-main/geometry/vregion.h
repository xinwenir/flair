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

#ifndef __VREGION_H
#define __VREGION_H

#include <ostream>
#include "os.h"
#include "geo.h"
#include "array.h"
#include "vector.h"
#include "gregion.h"

class OBBox;
class GZone;
class VZone;
class VBody;
class GeometryKernel;
class GeometryEngine;

/** Viewport Region class */
class VRegion {
private:
	GRegion*	_region;	// region reference
	Array<VZone*>	_zones;		// List of zones
	GeometryKernel*	_kernel;	// Kernel where region belongs

public:
	char	label[16];	// label or value to display

private:
	int	_generation;	// Generation number, to be checked for validity
				// against the _region generation
	dword   _hash;
public:
	bool	location;	// if true viewport location vs region is known (given by in)
	VZone	*in;		// if viewport is inside or outside the region
	dword	_color;		// region color
	double	_value;		// region color value
	int	_alpha;		// Transparency level [0..255]

public:
	VRegion(GRegion *reg=NULL, GeometryKernel *k=NULL) :
		_color(0xA0A0A0),
		_value(0.0),
		_alpha(0)
		{ label[0] = 0; init(reg, k); }
	~VRegion()			{ clear(); _region=NULL; }

	void	init(GRegion *reg=NULL, GeometryKernel *k=NULL);
	void	clear();

	GRegion* region()	const	{ return _region; }
const	char*	name()		const	{ return _region->name(); }
	RegionType type()	const	{ return _region->type(); }
	int	id()		const	{ return _region->id(); }
	dword	hash()		const	{ return _hash; }

	int	show()		const	{ return _region->show; }
  //	int	rotdefi()	const	{ return _region->rotdefi; }

	void	color(dword c)		{ _color = (c & FLAG_COLORMASK) | FLAG_REGION; }
	dword	color()		const	{ return _color; }

	void	value(double v)		{ _value = v; }
	double	value()		const	{ return _value; }

  //	void	alpha(byte a)		{ _alpha = a;     }
  //	byte	alpha()		const	{ return _alpha;      }
        void    alpha(uint8_t a)                { _alpha = a;     }
        uint8_t alpha()         const   { return _alpha;      }
	bool	transparent()	const	{ return _alpha==255; }
	bool	notTransparent()const	{ return _alpha==0;   }

	// Generation/valid/invalid
	bool	valid()		const	{ return _generation == _region->generation(); }
	bool	invalid()	const	{ return _generation != _region->generation(); }
	void	setInvalid()		{ _generation = -1;  }
	void	setValid()		{ _generation = _region->generation(); }

	int	nzones()	const	{ return _zones.size(); }
const	Array<VZone*>& zones()	const	{ return _zones; }
	VZone*	zone(int i)	const	{ return _zones[i]; }

	bool	ignore()	const	{ return (location && !in); }

	void	updateLocation();
	VZone*	inside2D(GeometryEngine* engine,
			 const double  x, const double  y, const double  z,
			 const double dx, const double dy, const double dz) const;
	VZone*	inside2D(GeometryEngine* engine,
			 const Vector& p, const Vector& d) const
			{ return inside2D(engine, p.x,p.y,p.z, d.x,d.y,d.z); }
	VZone*	inside(GeometryEngine* engine,
			 const double  x, const double  y, const double  z,
			 const double dx, const double dy, const double dz) const;
	VZone*	inside(GeometryEngine* engine,
			 const Vector& p, const Vector& d) const
			{ return inside(engine, p.x,p.y,p.z, d.x,d.y,d.z); }
	VZone*	insideRay(GeometryEngine* engine,
			  const double  x, const double  y, const double  z,
			  const double dx, const double dy, const double dz,
			  const double t ) const;
	VZone*	insideRay(GeometryEngine* engine,
			  const Vector& p, const Vector &d, const double t) const
			{ return insideRay(engine, p.x,p.y,p.z, d.x,d.y,d.z, t); }

	void	showBodies(int show);

	OBBox*	obbox() const { return _region->obbox(); };

	size_t	memory() const;

friend	class	VZone;
}; // VRegion

std::ostream& operator << (std::ostream&, const VRegion&);

#endif
