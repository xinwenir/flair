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
 *
 */

#include <ostream>

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

#include "geo.h"
#include "gbody.h"
#include "gzone.h"
#include "timer.h"
#include "vector.h"
#include "painter.h"
#include "bintree.h"
#include "matrix4.h"
#include "geometry.h"

using namespace std;

bool Geometry::developer = false;

/* ============================== Geometry ================================= */

/** Geometry */
Geometry::Geometry() :
	bodies(32),
	regions(32),
	materials(32),
//	_rotdefi(8),
//	_invRotdefi(8),
	_editBody(NULL),
	_editRegion("<editReg>")
{
	gridStep     =    3;
	lighterLevel =    0;
	_nrotdefi    =    0;
	_rotdefi     = NULL;
	_invRotdefi  = NULL;
	axisLen      =   25;
	_msg[0]      =    0;

	// lights
	lights       = 0;
	_generation  = 0;

	_backgroundColor = 0x707070 | FLAG_TRANSPARENT;
	errorColor       = 0xFF0000;
	gridTextColor    = 0xA0A0A0;
	regionColor      = 0x000000;
	regionErrorColor = 0xFFB0B0;
	selectColor      = 0xFF00FF;
	titleColor       = 0xD0D000;
	vertexColor      = 0xFFFF00;
	visibleColor     = 0x7F007F;
	latticeColor     = 0xFFFFFF;
	voxelColor       = 0xFFFFFF;
	zoneColor        = 0x606060;
	messageColor     = 0xFF9000;
	labelColor       = 0xA05010;
	_transparentColor= 0xFFFFFF | FLAG_TRANSPARENT;
	makeLatticeColor();

	pthread_rwlock_init(&rwlock, NULL);
	pthread_mutex_init(&mutexEdit, NULL);
} // Geometry

/** ~Geometry */
Geometry::~Geometry()
{
	// clean up local structures
	cleanup();
	if (_rotdefi) {
		delete [] _rotdefi;
		delete [] _invRotdefi;
	}
	pthread_rwlock_destroy(&rwlock);
	pthread_mutex_destroy(&mutexEdit);
} /* ~Geometry */

/** cleanup everything. All threads must be stopped */
void Geometry::cleanup()
{
	// Clean local arrays of bodies
	delRegions();
	delBodies();
	delMaterials();

	lockEdit();
	_editBody = NULL;
	_editRegion.clear();
	unlockEdit();

	// Do not cleanup the voxel maybe we needed latter

	// remove transformations
//	_rotdefi.clear();
//	_invRotdefi.clear();
	if (_rotdefi) {
		for (int i=0; i<_nrotdefi; i++) {
			if (_rotdefi[i]) {
				delete _rotdefi[i];
				delete _invRotdefi[i];
			}
		}
		memset(_rotdefi,    0, _nrotdefi*sizeof(Matrix4*));
		memset(_invRotdefi, 0, _nrotdefi*sizeof(Matrix4*));
	}
} // cleanup

/** add a rotdefi rotation matrix (+ inverse)
 * @param id		id of rotdefi (1 based)
 * @param matrix	rotation matrix of rotdefi
 * */
void Geometry::rotdefi(int id, const Matrix4& matrix)
{
	assert(id>0);
	id -= 1;
#if 0
	if (id >= _rotdefi.capacity()) {
		_rotdefi.resize(id);
		_invRotdefi.resize(id);
	}
#endif
	if (id >= _nrotdefi) {
		// find next size
		int nold = _nrotdefi;
		int sz = id+1;
		_nrotdefi = (sz%ROTDELTA)? ((sz+ROTDELTA)/ROTDELTA)*ROTDELTA : sz;
		if (_rotdefi) {
			Matrix4** old = _rotdefi;
			_rotdefi =  new Matrix4*[_nrotdefi];
			if (!_rotdefi) { // ERROR
				_rotdefi = old;
				return;
			}
			memcpy(_rotdefi, old, nold*sizeof(Matrix4*));
			delete [] old;

			old = _invRotdefi;
			_invRotdefi = new Matrix4*[_nrotdefi];
			if (!_invRotdefi) { // ERROR
				_invRotdefi = old;
				return;
			}
			memcpy(_invRotdefi, old, nold*sizeof(Matrix4*));
			delete [] old;
		} else {
			_rotdefi    = new Matrix4*[_nrotdefi];
			_invRotdefi = new Matrix4*[_nrotdefi];
		}
		// zero the rest
		memset(_rotdefi+nold,    0, (_nrotdefi-nold)*sizeof(Matrix4*));
		memset(_invRotdefi+nold, 0, (_nrotdefi-nold)*sizeof(Matrix4*));
	}
	if (_rotdefi[id] == NULL) {
		_rotdefi[id]    = new Matrix4(matrix);
		_invRotdefi[id] = new Matrix4(matrix);
	}
	_rotdefi[id]->copy(matrix);
	_invRotdefi[id]->copy(matrix);
	_invRotdefi[id]->inverse();
} // rotdefi

/** makeLatticeColor */
void Geometry::makeLatticeColor()
{
	latticeColor    |= FLAG_LATTICE;
	voxelColor      |= FLAG_VOXEL;
} // makeLatticeColor

/** addBody */
GBody* Geometry::addBody(const char *name, const char *type)
{
	GBody *body = GBody::newBody(name, type);
	body->id(bodies.count());
	body->nextGeneration(nextGeneration());

	bodies.add(body);
	bodiesTree.add(body->name(), body);

	return body;
} // addBody

/** delBodies */
void Geometry::delBodies()
{
	for (int i=0; i<bodies.count(); i++)
		delete bodies[i];
	bodies.clear();
	bodiesTree.destroyTree();
} // delBodies

/** invalidateBody */
void Geometry::invalidateBody(GBody *body)
{
	assert(body);
	body->nextGeneration(nextGeneration());
	ArrayIterator<GZone*> iter(body->zones);
	while (iter) {
		GZone *zone = iter++;
		zone->nextGeneration(nextGeneration());
		zone->region->nextGeneration(nextGeneration());
	}
} // invalidateBody

/** addRegion */
GRegion* Geometry::addRegion(const char *name)
{
	GRegion *region = new GRegion(name);
	region->id(regions.count());
	region->nextGeneration(nextGeneration());

	// XXX we could use the region->id to check if the region
	// is already defined
	regions.add(region);
	regionsTree.add(region->name(), region);

	return region;
} // addRegion

/** delRegions */
void Geometry::delRegions()
{
	for (int i=0; i<regions.count(); i++)
		delete regions[i];
	regions.clear();
	regionsTree.destroyTree();
} // delRegions

/** addZone to region */
GZone* Geometry::addZone(GRegion *region, const bool mode, const int size)
{
	GZone *zone = region->addZone(mode,size);
	zone->nextGeneration(nextGeneration());
	return zone;
} // addZone

/** add2exp add token to regions expression
 * @param region	region to add token to
 * @param token		token to be added either body name or + - | ()
 * @param err		returned error string in case of problem
 * @return true on error, false otherwise
 */
bool Geometry::add2exp(GRegion *region, const char *token, char *err)
{
	assert(region);
	assert(token);
	if (err) err[0] = 0;
	GBody *body = getBody(token);
	if (!region->add2exp(token, body)) {
		if (err) sprintf(err,"invalid token '%s'", token);
		return true;
	}
	return false;
} // add2exp

/** addMaterial */
Material* Geometry::addMaterial(const char *name)
{
	Material *material = new Material(name);
	material->id(materials.count());
	materials.add(material);
	materialsTree.add(material->name(), material);
	return material;
} // addMaterial

/** delMaterials */
void Geometry::delMaterials()
{
	for (int i=0; i<materials.count(); i++)
		delete materials[i];
	materials.clear();
	materialsTree.destroyTree();
} // delMaterials

/** addLight to the light list
 * @param l	light structure to add
 * @return true on error, false otherwise
 */
bool Geometry::addLight(Light& l)
{
	if (lights >= MAXLIGHT) return true;
	memcpy(&(light[lights]), &l, sizeof(l));
	light[lights].distance = light[lights].dir.normalize();
	lights++;
	return false;
} // addLight

/** set default lights */
void Geometry::defaultLights()
{
	// http://www.secondpicture.com/tutorials/3d/three-point_lighting_in_3ds_max_01.html
	lights = 0;
	light[lights].type     = LIGHT_SUN;			// sun light
	light[lights].relative = true;
	light[lights].dir      = Vector( 5.,-8.,-10.);
	light[lights].distance = INFINITE;
	light[lights].power    = 0.8;
	light[lights].spec     = true;
	light[lights].falloff  = 0;
	light[lights].shadow   = true;
	lights++;

	light[lights].type     = LIGHT_SUN;			// fill light
	light[lights].relative = true;
	light[lights].dir      = Vector(-1.,1.,-1.);
	light[lights].distance = INFINITE;
	light[lights].power    = 0.25;
	light[lights].spec     = false;
	light[lights].falloff  = 0;
	light[lights].shadow   = false;
	lights++;

	light[lights].type     = LIGHT_SUN;			// back light
	light[lights].relative = true;
	light[lights].dir      = Vector(-4.,-2.,6.);
	light[lights].distance = INFINITE;
	light[lights].power    = 0.8;
	light[lights].spec     = false;
	light[lights].falloff  = 0;
	light[lights].shadow   = false;
	lights++;

	for (int i=0; i<lights; i++)
		light[i].dir.normalize();
} // defaultLights

/** set message to display in all viewports */
void Geometry::message(const char *m)
{
	strncpy(_msg, m, sizeof(_msg));
	_msg[sizeof(_msg)-1] = 0;
} // message

/** printInfo */
void Geometry::printInfo() const
{
	cout << "Geometry Info" << endl;
	cout << "Number of bodies: " << bodies.count() << endl;
	cout << "Number of regions: " << regions.count() << endl;
	int nzones = 0;
	ArrayIterator<GRegion*> region_iter(regions);
	while (region_iter) {
		GRegion *region = region_iter++;
		nzones += region->zones().count();
	}
	cout << "Number of zones: " << nzones << endl;
} // printInfo

/** bodiesMemory */
size_t Geometry::bodiesMemory() const
{
	size_t bmem = bodies.memory() + bodiesTree.memory();
	for (int i=0; i<bodies.count(); i++)
		bmem += bodies[i]->memory();
	return bmem;
} // bodiesMemory

/** regionsMemory */
size_t Geometry::regionsMemory() const
{
	size_t rmem = regions.memory() + regionsTree.memory();
	for (int i=0; i<regions.count(); i++)
		rmem += regions[i]->memory();
	return rmem;
} // regionsMemory

/** memory */
size_t Geometry::memory() const
{
	return   sizeof(Geometry)
	       + bodiesMemory()
	       + regionsMemory();
//	       + voxel.memory();
} // memory

/** printMemory */
void Geometry::printMemory() const
{
	cout << endl << "Geometry:" << endl;
	cout << "Bodies:\t\t" << bodies.count()   <<  endl;
	cout << "Regions:\t"  << regions.count()  <<  endl;
	cout << "Memory:"     << endl;
	cout << "\tSelf:\t"   << sizeof(Geometry) << endl;
	cout << "\tBodies:\t" << bodiesMemory()   << endl;
	cout << "\tRegion:\t" << regionsMemory()  << endl;
	cout << "\tTotal:\t"  << memory() << endl;
} // printMemory
