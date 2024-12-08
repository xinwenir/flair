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

#ifndef __GREGION_H
#define __GREGION_H

#include <iostream>

#include "geo.h"
#include "bbox.h"
#include "array.h"
#include "vector.h"
#include "material.h"

class GBody;
class GZone;
class OBBox;

// Special region types
enum RegionType {
	REGION_NORMAL    = 0,
	REGION_BLACKHOLE = 1,
	REGION_LATTICE   = 2,
	REGION_VOXEL     = 3
};

/** Geometry Region class */
class GRegion {
protected:
	int	_id;		// index id in regions array
	char16	_name;		// name of the region
	RegionType _type;	// type of region: normal/blackhole/lattice/voxel
	int	_generation;	// Generation number

//	double	_volume;	// volume of region

	Array<GZone*>	_zones;	// zones defining the region
mutable	dword	_hash;		 // Cached hash value
mutable	OBBox	*cached_obbox;
	Material *_material;
	bool    _hasMatrix;	// if rotation matrix is set
	Matrix4	_matrix;	// rotation matrix
	Matrix4 _invMatrix;	// inverse rotation matrix

public:
  //	int	rotdefi;	// Rotdefi index if any
				//       0    nothing
				//     1.. n  1-based
				//    -1..-n  for inverse
	int	show;		// show region

public:
	GRegion(const char *aname);
	~GRegion()			{ clear(); };
	void	clear();

	void	name(const char *aname);
const	char*	name()		const	{ return _name; }

	void	type(const RegionType t) { _type = t; }
	RegionType type()	const	{ return _type; }

	void	id(const int i)		{ _id = i; }
	int	id()		const	{ return _id; }
	// Matrix and transformations
	void	clearMatrix()		{ _hasMatrix = false; }
	bool	hasMatrix()	 const	{ return _hasMatrix; }
	void	matrix(const Matrix4& M);
const	Matrix4& matrix()	const	{ return _matrix; }
const	Matrix4& invMatrix()	const	{ return _invMatrix; }


	void	material(Material* mat)	{ _material = mat; }
const	Material* material()	const	{ return _material; }

	int	nzones()	const	{ return _zones.count(); }
const	Array<GZone*>& zones()	const	{ return _zones; }
const	GZone*	zone(int i)	const	{ return _zones[i]; }
	GZone*	addZone(const bool mode=false, const int size=0);
	GZone*	addZone(GZone *z);
	void	delZone(int i);
	bool	add2exp(const char *token, GBody *body);

	// Generation
	int	generation()	const	{ return _generation; }
	void	nextGeneration(int newGen) {
			_generation = newGen;
			clearOBB();
			_hash = 0;
			//cout << "GRegion: " << name() << " NEXT GEN=" << _generation << endl;
		}

	// Hash
	dword	hash()   const;

	// bounding box
	BBox	bbox2D() const;
	BBox	bbox()   const;
	OBBox*	obbox()  const {
			if (!cached_obbox)
				cached_obbox = updateOBB();
			return cached_obbox;
		};

	// Inside (WARNING: non-caching functions)
	GZone*	inside( const double  x, const double  y, const double  z,
			const double dx, const double dy, const double dz) const;
	GZone*	inside( const Vector& p, const Vector& d) const
			{ return inside(p.x, p.y, p.z, d.x, d.y, d.z); }

	size_t	memory() const;

private:
	OBBox*	updateOBB() const;
	void	clearOBB();

friend class VRegion;
}; // GRegion

std::ostream& operator << (std::ostream&, const GRegion&);
#endif
