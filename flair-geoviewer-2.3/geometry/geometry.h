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
 * The geometry module is split into the following classes
 *
 * Classes Prefix:
 *    Gxxx: geometry information (Geometry)
 *    Cxxx: caching information, (Engine)
 *    Vxxx: viewer related information (Viewer)
 *
 * Meaning:
 *    +---	owns the element
 *    +-->      pointer or reference to element
 *
 *    GBody
 *    =====
 *    +--- BodyType type
 *    +--- char16   name
 *    +--- .... parameters
 *                Holder of the geometrical properties
 *                of each body. Provides setting
 *                intersection calculation, ...
 *                as well save/restore functions
 *
 *    GZone
 *    =====
 *    +--- GBody    expr[]
 *    +--- BBox     bbox
 *                Description of a single zone (normal)
 *                or an advanced with RPN representation.
 *                Calculation of intersection, inclusion
 *                etc..
 *
 *    GRegion
 *    +--- GZone    zones[]
 *    +--- int      rotdefi
 *    +--- RegionType type
 *    +--- Material material
 *                Region as a union of a list of zones
 *
 *    CBody
 *    +--- float    tmin,tmax
 *    +--> GBody    body
 *                 Cached body information, per viewport
 *                 (or per kernel) with the entry,exit
 *                 information
 *
 *    CZone, CRegion - no longer exist as it is faster
 *                 to recalculate than to cache the
 *                 information
 *
 *    Geometry
 *    ========
 *    +--- GBody    bodies
 *    +--- GRegion  regions
 *    +--- GVoxel   voxel
 *    +--- Matrix4  rotdefi, invRotdefi
 *    +--- Material materials
 *    +--- Light    lights
 *                Owner class of geometry properties
 *                bodies, regions, lattices, lights, materials etc.
 *                Contains all methods for manipulating
 *                the geometry structure
 *
 *    Engine
 *    ======
 *    +--- CBody bodies[]
 *    +--- CZone zonesSorted[]
 *    +--- CBody clipbody
 *    +--- CBody projectbody
 *    +--> Geometry
 *    +--> Kernel
 *                An asynchronous low level class that
 *                owns C(ached) structures calculating and
 *                caching information on the viewport/geometry
 *
 *    Kernel
 *    ======
 *    +--- Engine   engine, engines (of threadpool)
 *    +--> Geometry geometry
 *    +--- VBody    bodies
 *    +--- VRegion  regions
 *    +--- VVoxel   voxel
 *    +--- Viewport view
 *    +--- ThreadPool pool
 *                Job scheduling class owner of the V(iewer)
 *                structures
 *
 *    Viewer
 *    ======
 *    +--> Geometry
 *    +--> Kernel
 *    +--- *Layer	*layer
 *                Holder class of viewer depended properties
 *                like layers, colors, etc Contains all methods
 *                for drawing either synchronous or asynchronous
 *                Each viewer is referring to the geometry
 *                as well the drawing kernels
 *
 * Author:	Vasilis.Vlachoudis@cern.ch
 * Date:	04-Feb-2010
 */

#ifndef __GEOMETRY_H
#define __GEOMETRY_H

#include "os.h"
#include "geo.h"
#include "array.h"
#include "color.h"
#include "voxel.h"
#include "light.h"
#include "bintree.h"
#include "gregion.h"
#include "material.h"

class GBody;
class GZone;
class Matrix4;

/* ============================= Geometry ============================= */
/** Geometry class */
class Geometry {
public:
	/* drawing flags */
	int	axisLen;		/** axis length			*/

	int	gridStep;		/** pixel step for grid		*/
	int	lighterLevel;		/** lighter non-selected regions*/

	/* colors */
	dword	_backgroundColor;	/** background Color		*/
	dword	errorColor;		/** error line color		*/
	dword	gridTextColor;		/** color of grid text		*/
	dword	latticeColor;		/** lattice color		*/
	dword	regionColor;		/** region boundary colors	*/
	dword	regionErrorColor;	/** region error color		*/
	dword	selectColor;		/** selected item color		*/
	dword	titleColor;		/** color for title		*/
	dword	vertexColor;		/** vertex color		*/
	dword	visibleColor;		/** visible item color		*/
	dword	voxelColor;		/** voxel color			*/
	dword	zoneColor;		/** zone line color		*/
	dword	messageColor;		/** message color		*/
	dword	labelColor;		/** labels color		*/

	Color3D select3DColor;		/** 3d selection color		*/
	Color3D	wireframeColor;		/** wireframe color		*/
	Color3D	bodyBBoxInColor;	/** body BB inside color	*/
	Color3D	bodyBBoxOutColor;	/** body BB outside color	*/
	Color3D	zoneBBoxColor;		/** zone BB color		*/
	Color3D	regionBBoxColor;	/** region BB color		*/

	int	lights;
	Light	light[MAXLIGHT];	/** lights in the system	*/

	GVoxel	voxel;			/** voxel object		*/

	BinTree bodiesTree;		/** tree with all bodies	*/
	Array<GBody*>	bodies;		/** list of bodies		*/

	BinTree	regionsTree;		/** tree with all regions	*/
	Array<GRegion*> regions;	/** list of regions		*/

	BinTree materialsTree;		/** tree with all materials	*/
	Array<Material*> materials;	/** list of materials		*/

protected:
	int	_nrotdefi;
	Matrix4	**_rotdefi;		/** transformations		*/
	Matrix4	**_invRotdefi;		/** inverse transformations	*/

	GBody*	_editBody;		/** currently editing body	*/
	GRegion	_editRegion;		/** editing region		*/
	dword	_transparentColor;	/** transparent/ignore Color	*/
	char	_msg[128];		/** msg to display in viewports	*/

	int    _generation;

mutable pthread_rwlock_t rwlock;        /** Read/write lock,		*/
					/*  protects geometry definition*/
mutable	pthread_mutex_t	mutexEdit;	/** editRegion  mutex		*/

public:
static	bool	developer;		/** developer flag		*/

public:
	Geometry();
	~Geometry();
	void	cleanup();

	/* mutex */
	int	lockWrite()	const	{ return pthread_rwlock_wrlock(&rwlock); }
	int	unlockWrite()	const	{ return pthread_rwlock_unlock(&rwlock); }
	int	lockRead()	const	{ return pthread_rwlock_rdlock(&rwlock); }
	int	unlockRead()	const	{ return pthread_rwlock_unlock(&rwlock); }
	int	lockEdit()	const	{ return pthread_mutex_lock(&mutexEdit); }
	int	unlockEdit()	const	{ return pthread_mutex_unlock(&mutexEdit); }

	void	backgroundColor(const dword color) { _backgroundColor = color | FLAG_TRANSPARENT; }
	dword	backgroundColor()	const	{ return _backgroundColor; }

	void	makeLatticeColor();

	void	rotdefi(int id, const Matrix4& matrix);
const	Matrix4& rotdefi(int id)	const	{ assert(id); return id>0? *_rotdefi[id-1]    : *_invRotdefi[-id-1]; }
const	Matrix4& invRotdefi(int id)	const	{ assert(id); return id>0? *_invRotdefi[id-1] : *_rotdefi[-id-1]; }

	dword	transparentColor()	const	{ return _transparentColor; }

	/* generation */
	int	nextGeneration()		{ return ++_generation; }

	/* bodies */
	GBody*	addBody(const char *name, const char *type);
	GBody*	getBody(const char *name)	{ return (GBody*)bodiesTree[name]; }
	GBody*	getBody(const int  id)		{ return (GBody*)bodies[id]; }	// WARNING check limits
const	GBody*	getBody(const char *name) const	{ return (GBody*)bodiesTree[name]; }
const	GBody*	getBody(const int  id)	const	{ return (GBody*)bodies[id]; }	// WARNING check limits
	void	delBodies();
	void	invalidateBody(GBody *body);

	/* regions */
	GRegion* addRegion(const char *name);
	GRegion* getRegion(const char *name)	{ return (GRegion*)regionsTree[name]; }
	GRegion* getRegion(const int  id)	{
			if (id >=0 && id < regions.count())
				return regions[id];
			else
				return NULL;
		}
const	GRegion* getRegion(const char *name) const { return (GRegion*)regionsTree[name]; }
const	GRegion* getRegion(const int  id) const {
			if (id >=0 && id < regions.count())
				return regions[id];
			else
				return NULL;
		}
	void	delRegions();

	/* zone setup */
	GZone*	addZone(GRegion *region, const bool mode, const int size=0);
	bool	add2exp(GRegion *reg, const char *token, char *err);

	/* materials */
	Material* addMaterial(const char *name);
	Material* getMaterial(const char *name)		{ return (Material*)materialsTree[name]; }
const	Material* getMaterial(const char *name) const	{ return (Material*)materialsTree[name]; }
	Material* getMaterial(const int  id)	{
			if (id >= 0 && id < materials.count())
				return (Material*)materials[id];
			else
				return NULL;
			}
const	Material* getMaterial(const int  id) const {
			if (id >= 0 && id < materials.count())
				return (Material*)materials[id];
			else
				return NULL;
			}
	void	delMaterials();

	/* lights */
	bool	addLight(Light& l);
	void	delLights()		{ lights = 0; }
	void	defaultLights();

const	GBody*	editBody()	const	{ return _editBody;   }
	GBody*	editBody()		{ return _editBody;   }
	void	editBody(GBody* b)	{ _editBody = b;      }

	GRegion& editRegion()		{ return _editRegion; }
const   GRegion& editRegion()	const	{ return _editRegion; }

	/* message */
	void	message(const char *);
const	char*	message()	const	{ return _msg; }

	/* Info */
	size_t	bodiesMemory()	const;
	size_t	regionsMemory()	const;
	size_t	memory()	const;
	void	printInfo()	const;
	void    printMemory()	const;

friend class GeometryViewer;
friend class DecorationLayer;
friend class D2Layer;
}; // Geometry

#endif
