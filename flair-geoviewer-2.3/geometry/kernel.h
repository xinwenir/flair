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

#ifndef __GEOMETRY_KERNEL_H
#define __GEOMETRY_KERNEL_H

//#include <iosfwd>
//#include <ostream>
//#include <iostream>

//#include <time.h>
//#include <stdio.h>
//#include <string.h>
#include <assert.h>
#include <pthread.h>

#include "os.h"
#include "geo.h"
#include "array.h"
#include "gbody.h"
#include "vbody.h"
#include "voxel.h"
#include "vector.h"
#include "vregion.h"
#include "gregion.h"
#include "matrix4.h"
#include "geometry.h"
#include "viewport.h"
#include "threadpool.h"

#include "engine.h"

//class BBox;
//class OBBox;
//class Vertex;
//class VBodyOrderAccel;

#ifdef _STAT
/* =========================== GeometryKernelStats ========================= */
/** Geometry Engine Statistics */
class GeometryKernelStats {
public:
	long long segments_total;
	long long segments_invalid;
	long long segments_skipped;

	long long valid_bodies;
	long long invalid_bodies;
	long long valid_regions;
	long long invalid_regions;

	GeometryKernelStats() {reset();}
	void reset();
}; // GeometryStats
std::ostream& operator << (std::ostream&, const GeometryKernelStats&);
#endif

#ifdef THREAD
class VolumeFeeder;
/* ============================== VolumeWorker ============================== */
class VolumeWorker : public ThreadPoolWorker {
public:
	GeometryEngine*	engine;		/** engine assigned to worker	*/
	int	idx;			/** index of event		*/
	int	to;			/** up to which event to calc	*/
	long	hits;			/** hits count			*/
public:
	VolumeWorker(): engine(NULL), idx(0), to(0), hits(0) {}
	virtual void operator()();
}; // VolumeWorker

/* ============================== VolumeFeeder ============================== */
class VolumeFeeder : public ThreadPoolFeeder {
private:
	int		nworkers;	/** number of workers allocated	*/
	VolumeWorker*	workers;	/** array of body workers	*/

	int	idx;			/** next start index		*/
	int	step;			/** step per worker		*/

public:
	BBox*	 bbox;			/** bounding box		*/
	GRegion* region;		/** region to check		*/

	int	samples;		/** maximum number of samples	*/
	int	randomNr1;		/** random scrabbler 1		*/
	int	randomNr2;		/** random scrabbler 2		*/

public:
	VolumeFeeder(ThreadPool& p) : ThreadPoolFeeder(p),
		  nworkers(0),
		  workers(NULL)
		  {}
	virtual	~VolumeFeeder()		{ if (workers) delete [] workers; }

	// Initialize workers
	void	reset(int nsamples, int nstep, int rnd1, int rnd2, BBox* bb, GRegion* reg);

virtual	ThreadPoolWorker* feed(int threadId);
	long	hits()		{ long h=0; for (int i=0; i<nworkers; i++) h += workers[i].hits; return h; }
}; // VolumeFeeder
#endif

/* ========================== GeometryKernel ========================== */
/** Geometry Kernel class */
class GeometryKernel {
protected:
const   Geometry&	geometry;	/** father geometry object	*/
	GeometryEngine	_engine;	/** main geometry engine	*/
#ifdef THREAD
	GeometryEngine*	threadengines;	/** geometry engine threads	*/
#endif

public:
	Array<VBody*>	bodies;		/** body array			*/
	Array<VRegion*> regions;	/** region array		*/
	VVoxel          voxel;		/** voxel structure		*/
	ViewPort	view;		/** Viewing window		*/

protected:
	double		_minSx2;	/** min Sx2 factor ignore points*/
	int		_errors;	/** geometry errors count	*/
	dword		_errorhash;	/** error hash variable		*/
	char		_errmsg[128];	/** error message		*/
#ifdef THREAD
	ThreadPool	threadpool;
#endif
#ifdef _STAT
	GeometryKernelStats stats;	/** Statistics                  */
#endif

	pthread_mutex_t	mutex;		/** thread mutex		*/
public:
	GeometryKernel(Geometry& g, int nt=-1);
	~GeometryKernel();

	void	clearBodies();
	void	clearRegions();

	/** window */
	void	window(double ex, double ey=0.0);
	void	window(double xmin, double ymin, double xmax, double ymax)
			{ view.window(xmin, ymin, xmax, ymax); }
	void	resize(int w, int h);
	int     width()		const	{ return view.width(); }
	int     height()	const	{ return view.height(); }

	GeometryEngine*	engine()	{ return &_engine; }

	void	initThreads(int nt=0);
#ifdef THREAD
	GeometryEngine*	engine(int i)	{ return &threadengines[i]; }
	int	nthreads()	const	{ return threadpool.nthreads(); }
#else
	int	nthreads()	const	{ return 0; }
#endif

	/* bodies */
	VBody*	getBody(const int id)	{ return bodies[id]; }	// XXX WARNING unprotected for limits
	VBody*	getBody(const char *name) {
			const GBody *body = geometry.getBody(name);
			return body?getBody(body->id()):NULL;
		}
	VBody*	getBody(const GBody* body)	{ return bodies[body->id()]; }
	void	clearVisibleBodies();
	int	nGeometryBodies()	const	{ return geometry.bodies.count(); }

	/* projection bodies */
	void	projectBodyClear();
	void	projectBodyAdd(int idx);
	int	projectBody(int i)	const	{ return _engine.projectBody(i); }
	int	projectBodyCount()	const	{ return _engine.projectBodyCount(); }

	/* clipping bodies */
	void	clipBodyClear();
	void	clipBodyAdd(int idx);
	void	clipBodyNegative(int idx);
	int	clipBody(int i)		const	{ return _engine.clipBody(i); }
	int	clipBodyCount()		const	{ return _engine.clipBodyCount(); }

	/* regions */
	VRegion* getRegion(const int id) {
			if (id>=0 && id<regions.size())
				return regions[id];
			else
				return NULL;
		}
	VRegion* getRegion(const char *name) {
			const GRegion *region = geometry.getRegion(name);
			return region? getRegion(region->id()) : NULL;
		}
	VRegion* getRegion(GRegion *region)	{ return regions[region->id()]; }
	int	 nGeometryRegions()	const	{ return geometry.regions.count(); }

	/* transformations */
	void	matrix(const Matrix4& m);
	void	origin(const double x, const double y, const double z);
	void	origin(double *x, double *y, double *z)	{ view.origin(x,y,z); }
	void	moveViewOriginTo0();

	/* errors */
	int	errors()	const	{ return _errors; }
	void	clearErrors()		{ _errors = 0; }
	void	clearErrorHash()	{ _errorhash = 0; }
	dword	errorHash()	const	{ return _errorhash; }
	void	setError()		{ if (!_errors) _errors = 1; }
	VBody*	error(const int err, double *u, double *v,
			ZoneOfPoint *pIn, ZoneOfPoint *pOut,
			double *umin, double *umax,
			double *vmin, double *vmax);

	/* error message */
	void	errorReset()		{ _errmsg[0] = 0; }
	void	error(const char *msg)	{ if (_errmsg[0]==0) strcpy(_errmsg, msg); }
const	char*	error()		const	{ return _errmsg; }

	/* body operations */
	void	makeBodyConics(bool all);
	void	intersectBody(VBody* a, bool all);
	void	scanDirection(VBody* body, const int cid,
			const double u, const double v,
			ZoneOfPoint* pIn, ZoneOfPoint* pOut);
	void	scanBodySegments(bool all);
	void	scanBodySegments(VBody* body, GeometryEngine* eng, bool all);

	/* region operations (MOVE TO KERNEL) */
	void	updateRegionLocation(bool all);
	bool	correctOverlaps(VRegion *region, int zoneid=-1);
	VBody*	zoneOverlaps(VZone* zone, int *side);
	int	bodySide(VZone* zone, VBody* body);
	int	segmentSide(VZone* zone, VBody* body, int cid, int from, int to);
	bool	isZoneExist(VZone* zone, VBody* body);

	/* geometry functions */
	double	volume(GRegion* region, int samples, double* vol, double* err);

	/* threads */
	void	lock() {
			pthread_mutex_lock(&mutex);
			geometry.lockRead();
		}
	void	unlock() {
			geometry.unlockRead();
			pthread_mutex_unlock(&mutex);
		}
	void	stop()	{
#ifdef THREAD
			threadpool.stop();
#endif
		}
	bool	isRunning() {
#ifdef THREAD
			return threadpool.isRunning();
#else
			return false;
#endif
		}

	/** Derive */
	void	derive();

	/* Info */
	size_t	bodiesMemory()	const;
	size_t	regionsMemory()	const;
	size_t	memory()	const;
	void	printMemory()	const;

friend class Layer;
friend class D2Layer;
friend class D3Layer;
}; // GeometryKernel
#endif
