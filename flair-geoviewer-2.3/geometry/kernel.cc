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
 * Changes:
 * - Added DXF export - 22.06.2010 C. Theis
 */

#include <math.h>
#include <time.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "bbox.h"
#include "edge.h"
#include "quad.h"
#include "bmath.h"
#include "conic.h"
#include "gzone.h"
#include "obbox.h"
#include "timer.h"
#include "vzone.h"
#include "engine.h"
#include "kernel.h"
#include "random.h"
#include "vertex2d.h"

// Body/Zone border type flags
#define	SIDE_NONE 0x00		// Body not found on any side (border)
#define	SIDE_IN   0x01		// Body found on side and is an overlap in the IN
#define	SIDE_OUT  0x02		// Body found on side and is an overlap in the OUT
#define	SIDE_OK   0x04		// Body found on side which is OK
#define SIDE_INOUT  (SIDE_IN|SIDE_OUT)
#define SIDE_BORDER (SIDE_IN|SIDE_OUT|SIDE_OK)

#define SIDE_STR(s)	(s&SIDE_IN? 'I':'-') << (s&SIDE_OUT?'O':'-') << (s&SIDE_OK? 'B':'-')

using namespace std;

#ifdef THREAD
/* =============================== VolumeFeeder ============================= */
/** operator() */
void VolumeWorker::operator()()
{
	VolumeFeeder* f = (VolumeFeeder*)feeder();
	Vector size = f->bbox->size();
	while (idx<to) {
		// sample coordinates using low-discrepancy series
		double x = f->bbox->low().x + size.x * VanDerCorput(idx, f->randomNr1);
		double y = f->bbox->low().y + size.y * Sobol2(idx, f->randomNr2);
		double z = f->bbox->low().z + size.z / (double)f->samples * idx;

		// check the geometry with containment tests
		if (f->region->inside(x,y,z,0.,0.,1.0)) ++hits;
		idx++;
	}
} // operator

/** reset */
void VolumeFeeder::reset(int nsamples, int nstep, int rnd1, int rnd2, BBox* bb, GRegion* reg)
{
	pool.reset();
	if (nworkers != pool.nthreads()) {
		if (workers) delete [] workers;
		nworkers = pool.nthreads();
		workers = new VolumeWorker[nworkers];
		for (int i=0; i<nworkers; i++)
			workers[i].feeder(this);
	}

	idx     = 0;
	step    = nstep;
	samples = nsamples;
	bbox    = bb;
	region  = reg;
	randomNr1 = rnd1;
	randomNr2 = rnd2;

	// initialize workers
	for (int i=0; i<nworkers; i++) {
		workers[i].idx    = 0;
		workers[i].to     = 0;
	}
} // reset

/** feed */
ThreadPoolWorker* VolumeFeeder::feed(int threadId)
{
	// Check for ending condition
	if (idx>=samples) return NULL;

	// update worker
	workers[threadId].idx = idx;
	workers[threadId].to  = idx = min(idx + step, samples);	// next value

	return (ThreadPoolWorker*)&workers[threadId];
} // feed
#endif

/* ============================== GeometryKernel ============================ */
/** GeometryKernel */
GeometryKernel::GeometryKernel(Geometry& g, int nt) :
	geometry(g),
	_engine(&g, this),
#ifdef THREAD
	threadengines(NULL),
#endif
	bodies(32),
	regions(32),
	voxel(g.voxel),
	view(640,480)
{
	_minSx2    = Sqr(1e-5); // min of a pixel to be ignored
	_errors    = 0;
	_errorhash = 0;

	// thread
	if (nt>=0) initThreads(nt);
	pthread_mutex_init(&mutex, NULL);
} // GeometryKernel

/** ~GeometryKernel */
GeometryKernel::~GeometryKernel()
{
	//stop();
	clearBodies();
	clearRegions();
#ifdef THREAD
	if (threadengines) delete [] threadengines;
#endif
	pthread_mutex_destroy(&mutex);
} /* ~GeometryKernel */

/** initialize threadpool
 * @param nt	number of threads, 0=all cpu/cores, 1=synchronous mode
 */
#ifdef THREAD
void GeometryKernel::initThreads(int nt)
#else
void GeometryKernel::initThreads(int)
#endif
{
#ifdef THREAD
	if (nt==0)	// 0 = use all available threads
		threadpool.init(ThreadPool::num_cores());
	else
	if (nt==1)	// 1 = run in sync mode, no threadpool
		threadpool.init(0);
	else	// n = up to max number of threads to use
#ifdef DEBUG
		// allow to up to what the user defined in debug mode
		threadpool.init(nt);
#else
		// limit to maximum number of cores
		threadpool.init(Min(nt, ThreadPool::num_cores()));
#endif
	if (threadengines) {
		delete [] threadengines;
		threadengines = NULL;
	}
#endif
} // initThreads

/** clearBodies */
void GeometryKernel::clearBodies() {
	for (int i=0; i<bodies.size(); i++)
		if (bodies[i])
			delete bodies[i];
	bodies.clear();
	clipBodyClear();
	projectBodyClear();
} // clearBodies

/** clearRegions */
void GeometryKernel::clearRegions() {
	for (int i=0; i<regions.size(); i++)
		if (regions[i])
			delete regions[i];
	regions.clear();
} // clearRegions

/** clear visible flag from all bodies */
void GeometryKernel::clearVisibleBodies()
{
	for (int i=0; i<nGeometryBodies(); i++) {
		VBody *vbody = bodies[i];
		vbody->visible = false;
	}
} // clearVisibleBodies

/** projectBodyClear */
void GeometryKernel::projectBodyClear()
{
	_engine.projectBodyClear();
#ifdef THREAD
	if (threadengines!=NULL)
		for (int i=0; i<nthreads(); i++)
			engine(i)->projectBodyClear();
#endif
} // projectBodyClear

/** projectBodyAdd */
void GeometryKernel::projectBodyAdd(int idx)
{
	_engine.projectBodyAdd(idx);
#ifdef THREAD
	for (int i=0; i<nthreads(); i++)
		engine(i)->projectBodyAdd(idx);
#endif
} // projectBodyAdd

/** clipBodyClear */
void GeometryKernel::clipBodyClear()
{
	_engine.clipBodyClear();
#ifdef THREAD
	if (threadengines!=NULL)
		for (int i=0; i<nthreads(); i++)
			engine(i)->clipBodyClear();
#endif
} // clipBodyClear

/** clipBodyAdd */
void GeometryKernel::clipBodyAdd(int idx)
{
	_engine.clipBodyAdd(idx);
#ifdef THREAD
	if (threadengines!=NULL)
		for (int i=0; i<nthreads(); i++)
			engine(i)->clipBodyAdd(idx);
#endif
} // clipBodyAdd

/** clipBodyNegative */
void GeometryKernel::clipBodyNegative(int idx)
{
	_engine.clipBodyNegative(idx);
#ifdef THREAD
	if (threadengines!=NULL)
		for (int i=0; i<nthreads(); i++)
			engine(i)->clipBodyNegative(idx);
#endif
} // clipBodyNegative

// Convert the Bodies,Zones,Regions to Array(VBody,VZone,VRegion)
// do not delete anything but reassign
//
/** derive - derive the VBodies, VZones and VRegions from Geometry
 * It initializes the list of bodies, zones and regions in the Viewer,
 * note that it deletes all previous objects and creates new ones
 */
void GeometryKernel::derive()
{
#ifdef _STAT
	Timer timer;
	timer.start();
#endif
	lock();

#if 0
	// WARNING: bodies and regions vectors can contain more or equal
	// number of VBodies vs geometry.bodies/regions
	//
	// Synchronize the bodies
	clearBodies();
	bodies.resize(nGeometryBodies());
	for (int i=0; i<nGeometryBodies(); i++) {
		GBody *gbody = geometry.bodies[i];
		bodies[i] = new VBody(gbody);
	}
#else
//#define DUMP(x) x
	DUMP(cout << endl);
	DUMP(cout << " ===================== DERIVE ======================" << endl);
	int oldSize = bodies.size();

	// Synchronize the bodies
	if (nGeometryBodies() < oldSize) {
		// Free the extra bodies if size if reduced
		for (int i=nGeometryBodies(); i<bodies.size(); i++)
			if (bodies[i]) {
				delete bodies[i];
				bodies[i] = NULL;
			}
		oldSize = nGeometryBodies();
	}

	bodies.allocate(nGeometryBodies());

	for (int i=0; i<oldSize; i++) {
		// Update existing elements
		GBody *gbody = geometry.bodies[i];
		VBody *vbody = bodies[i];
		DUMP(cout << "VBody " << gbody->name() << " eq=" << (vbody->body() == gbody) <<
			" invalid=" << vbody->invalid() << endl);
		if (vbody->body() != gbody || vbody->invalid()) {
			DUMP(cout << "VBody " << gbody->name() << " [UPDATE]" << endl);
			vbody->init(gbody);
		} else {
			DUMP(cout << "VBody " << gbody->name() << " [OK]" << endl);
		}
	}
	for (int i=oldSize; i<nGeometryBodies(); i++) {
		GBody *gbody = geometry.bodies[i];
		DUMP(cout << "VBody " << gbody->name() << " [NEW]" << endl);
		bodies[i] = new VBody(gbody);
	}
#endif

#if 0
	clearRegions();
	regions.resize(nGeometryRegions());
	for (int i=0; i<nGeometryRegions(); i++) {
		GRegion *gregion = geometry.regions[i];
		regions[i] = new VRegion(gregion, this);
	}
#else
	// Synchronize regions
	oldSize = regions.size();
	if (nGeometryRegions() < oldSize) {
		// Free the extra bodies if size if reduced
		for (int i=nGeometryRegions(); i<regions.size(); i++)
			if (regions[i]) {
				delete regions[i];
				regions[i] = NULL;
			}
		oldSize = nGeometryRegions();
	}

	regions.allocate(nGeometryRegions());

	for (int i=0; i<oldSize; i++) {
		GRegion *gregion = geometry.regions[i];
		VRegion *vregion = regions[i];
		DUMP(cout << "VRegion " << gregion->name() << " eq=" << (vregion->region() == gregion) <<
			" hash_eq=" << (vregion->hash() == gregion->hash()) <<
			" ggeneration=" << gregion->generation() <<
//			" vgeneration=" << vregion->generation() <<
			" invalid=" << vregion->invalid() << endl);
		if(vregion->region() != gregion || vregion->hash() != gregion->hash() || vregion->invalid()) {
			DUMP(cout << "VRegion " << gregion->name() << " [UPDATE]" << endl);
			vregion->init(gregion, this);
		} else {
			DUMP(cout << "VRegion " << gregion->name() << " [OK]" << endl);
		}
	}
	for (int i=oldSize; i<nGeometryRegions(); i++) {
		GRegion *gregion = geometry.regions[i];
		DUMP(cout << "VRegion " << gregion->name() << " [NEW]" << endl);
		regions[i] = new VRegion(gregion, this);
	}
	DUMP(cout << " ===================================================" << endl);
	DUMP(cout << endl);
#endif

	_engine.derive();
#ifdef THREAD
	if (threadengines==NULL && nthreads()>0) {
		threadengines = new GeometryEngine[nthreads()];
		for (int i=0; i<nthreads(); i++) {
			threadengines[i].geometry(&geometry);
			threadengines[i].kernel(this);
		}
	}
	for (int i=0; i<nthreads(); i++)
		threadengines[i].derive();
#endif

	unlock();
#ifdef _STAT
	timer.stop();
	cerr << "Time to derive geometry: " << timer << endl;
#endif
} // derive

/** origin */
void GeometryKernel::origin(const double x, const double y, const double z)
{
	stop();
	view.origin(x,y,z);
} // origin

/** moveViewOriginTo0
 * Move view origin to center of "view"
 */
void GeometryKernel::moveViewOriginTo0()
{
	stop();
	view.moveOriginTo0();
} // moveViewOriginTo0

/** matrix */
void GeometryKernel::matrix(const Matrix4& m)
{
	stop();		// Stop thread before modifying the matrix
	view.matrix(m);
} // matrix

/** window */
void GeometryKernel::window(double ex, double ey)
{
	if (ISZERO(ey)) ey = ex;
	view.window(-ex,-ey, ex,ey);
} // window

/** resize */
void GeometryKernel::resize(int w, int h)
{
	view.init(w, h);
} // resize

/* makeBodyConics
 * Creates the viewport conics for all the bodies
 * takes into account the _projectAll flag, otherwise it only affects invalid
 * bodies
 */
void GeometryKernel::makeBodyConics(bool all)
{
	_engine.incBodyCheckId();
	for (int i=0; i<nGeometryBodies(); i++) {
		VBody* vbody = getBody(i);
		if (all || vbody->invalid()) {
			STATS(stats.invalid_bodies++);
			vbody->makeConics(view);
			vbody->intersectSelf(view);
			vbody->intersectViewport(view);
			vbody->updateLocation(view);
			DUMP(cout << "VBody:" << *vbody << endl);
		} else
		if (!all) {
			STATS(stats.valid_bodies++);
			vbody->removeInvalidVertices();
		}
	}
} // makeBodyConics

/** intersectBody */
void GeometryKernel::intersectBody(VBody* a, bool all)
{
#if 1
	// Generate all body combinations and intersect them
	for (int j=a->id()+1; j<nGeometryBodies(); j++) {
		VBody* b = getBody(j);
		if (b->ignore()) continue;
		if (all || a->invalid() || b->invalid())
			a->intersectBody(b, view);
	}
#else
	// it goes slightly slower (I don't know why?)
	// Use a pseudo-random starting point
	int from = a->id()+1;
	int end   = _viewer->nGeometryBodies();
	if (end==from) return;

	// start from a random position
	int start = (from * 12997 + 981251) % (end-from) + from;

	// Generate all body combinations and intersect them
	// Split into two loops just to avoid workers processing the same bodies
	for (int j=start; j<end; j++) {
		VBody* b = _viewer->getBody(j);
		if (b->ignore()) continue;
		if (all || a->invalid() || b->invalid())
			a->intersectBody(b, view);
	}

	for (int j=from; j<start; j++) {
		VBody* b = _viewer->getBody(j);
		if (b->ignore()) continue;
		if (_projectAll || a->invalid() || b->invalid())
			a->intersectBody(b, view);
	}
#endif
} // intersectBody
/** compare zone list a against zone list b if regions match
 * @return true if are different
 */
static inline bool equalZones(VZone *a[], VZone *b[], int n)
{
	for (int i=0; i<n; i++) {
		// check first the same position for faster lookup
		if (a[i]->region() == b[i]->region()) continue;

		bool found = false;
		for (int j=0; j<n; j++) {
			if (i == j) continue;
			if (a[i]->region() == b[j]->region()) {
				found = true;
				break;
			}
		}

		if (!found)
			return false;
	}
	return true;
} // equalZones

/** scanDirection fill the pIn, pOut structure with the search position/direction
 * @param body	body whose conic to scan
 * @param cid	conic id
 * @param u	u location of the middle segment
 * @param v	v location of the middle segment
 * @param pIn	zone of point with inward direction
 * @param pOut	zone of point with outward direction
 */
void GeometryKernel::scanDirection(VBody* body, const int cid,
		const double u, const double v,
		ZoneOfPoint* pIn, ZoneOfPoint* pOut)
{
	// --------------------
	//   get the gradient
	// --------------------
	double tx = view.uv2x(u,v);
	double ty = view.uv2y(u,v);
	double tz = view.uv2z(u,v);

	double dx, dy, dz;
	body->Q(body->c2q[cid]).grad(tx,ty,tz,&dx,&dy,&dz);
	// normalize
	double dd = 1.0/sqrt(dx*dx+dy*dy+dz*dz);
	dx *= dd;
	dy *= dd;
	dz *= dd;

	/*
	 *             S\  ^ w
	 *               \ |
	 *                \|_-'normal
	 * ---view plane---+-------------------
	 *                / \`-_OUT
	 *               /   \
	 *              IN    \
	 * calculate the half angle from the surface
	 * with the view plane as IN and OUT
	 *
	 * IN = -normal-z
	 * AFTER the calculation of the shift
	 * adjust for the case of dz=0.0
	 */
	pIn->dx = -dx - view.matrix(0,2);
	pIn->dy = -dy - view.matrix(1,2);
	pIn->dz = -dz - view.matrix(2,2);

	//OUT = normal-z
	pOut->dx = dx - view.matrix(0,2);
	pOut->dy = dy - view.matrix(1,2);
	pOut->dz = dz - view.matrix(2,2);

	// Normalize directions
	dd = 1.0/sqrt(Sqr(pIn->dx) + Sqr(pIn->dy) + Sqr(pIn->dz));
	pIn->dx *= dd;
	pIn->dy *= dd;
	pIn->dz *= dd;
	dd = 1.0/sqrt(Sqr(pOut->dx) + Sqr(pOut->dy) + Sqr(pOut->dz));
	pOut->dx *= dd;
	pOut->dy *= dd;
	pOut->dz *= dd;
#if 1
	// Move a bit in the direction of the search
	pIn->x  = tx + pIn->dx *SMALL*Abs(tx);
	pIn->y  = ty + pIn->dy *SMALL*Abs(ty);
	pIn->z  = tz + pIn->dz *SMALL*Abs(tz);
	pOut->x = tx + pOut->dx*SMALL*Abs(tx);
	pOut->y = ty + pOut->dy*SMALL*Abs(ty);
	pOut->z = tz + pOut->dz*SMALL*Abs(tz);
#else
	pIn->x  = tx;
	pIn->y  = ty;
	pIn->z  = tz;
	pOut->x = tx;
	pOut->y = ty;
	pOut->z = tz;
#endif
} // scanDirection

/** scanSegments, a simple non-threaded version that scans
 * segments only of referenced bodies
 */
void GeometryKernel::scanBodySegments(bool all)
{
	clearErrorHash();
	for (int ib=0; ib<nGeometryBodies(); ib++) {
		VBody* b = getBody(ib);
		if (b->notReferenced()) continue;
		scanBodySegments(b, &_engine, all);
	}
} // scanSegments

/** scanBodySegments */
void GeometryKernel::scanBodySegments(VBody* body, GeometryEngine* eng, bool all)
{
	ZoneOfPoint pIn;
	ZoneOfPoint pOut;

#ifdef _DUMP
	if (body->show()) {
		cout << "********************************" << endl;
		cout << "+++ VBody-segments: " << body->name() << endl;
		for (int jj=0; jj < body->nC; jj++) {
			cout << "   Acc: " << body->acc[jj] << endl;
			for (int ii=0; ii<body->V[jj].size(); ii++) {
				Vertex2D& vv = body->V[jj][ii];
				cout << "   " << ii << " "
				     << (vv.body?vv.body->name():"NULL")
				     << '\t' << vv.t
				     << " [" << vv.x << ", " << vv.y << "]"
				     << " [" << view.uv2x(vv.x, vv.y)
				     << ", " << view.uv2y(vv.x, vv.y)
				     << ", " << view.uv2z(vv.x, vv.y) << "]"
				     << " invalid=" << vv.invalid
				     << endl;
			}
			cout << endl;
		}
		cout << "********************************" << endl;
	}
#endif
	// Check if any of the zones that refer to this body is visible
	// in the current viewport
	bool referenced = !body->notReferenced();

	for (int i=0; i < body->nC; i++) {
		DUMP(cout << "+++ " << body->name() << " conic=" << i << endl);
		Array<Vertex2D>& V = body->V[i];
		bool isHyperbola = body->C[i].type()==CONIC_HYPERBOLA;
		bool isEllipse   = body->C[i].type()==CONIC_ELLIPSE;

		if (V.size() <= 1) continue;
		pIn.body  = body;
		pOut.body = body;
		pIn.V = pOut.V = &V;
		Vertex2D* v = &V[0]; // Start point

		double tp = v->t;
		double xp = v->x;
		double yp = v->y;
		v->err = 0;
		v->type = SEGMENT_IGNORE;
		v->invalid = false;

		pIn.zero();
		pOut.zero();

		STATS(stats.segments_total += V.size());

		// Segment invalid flag
		bool invalid = all || body->invalid();

		eng->incBodyCheckId();

		// Go through all the segments and find their inside and
		// outside zones by checking at midpoint
		for (int j=1; j<V.size(); j++) {
			v = &V[j]; // next endpoint
			double t = v->t;
			double x = v->x;
			double y = v->y;

#ifdef EXPERIMENTAL
			// Reset the checkId of all objects crossed
			if (v->body!=NULL && v->body!=body)
				eng->getBody(v->body)->resetCheckId();
			else
				eng->incBodyCheckId();
#endif

			DUMP(cout << "<*< "<<body->name()<<" conic=" << i<<" vertex="<<j
				<< "\tv->body: " << (v->body? v->body->name() :"<win>") << endl);
			DUMP(cout << "-*- "<<body->name()<<" x=" <<x<<" y="<<y<<" t="<<t<<endl);

			// Check if there is a reference to invalid body
			if (v->invalid) invalid = true;

			// Skip singularities
			if (Sqr((x-xp)*view.Sx())+Sqr((y-yp)*view.Sy())<_minSx2 || t-tp<SMALL)
				// do not skip empty ellipse with t=PI and tp=-PI (only two points)
				if (!isEllipse || (isEllipse && j!=1)) {
					STATS(stats.segments_skipped++);
					// set the previous type
					tp = t;
					xp = x;
					yp = y;
					v->err = 0;
					v->type = SEGMENT_IGNORE;
					v->invalid = false;

					// reset the current body for the next search
					// increase the id
					// if goes outside the window or
					// intersect it self
#if 0
					// Reset the checkId of all objects crossed
					if (v->body != NULL && v->body != &body)
						eng->getBody(v->body)->resetCheckId();
					else
						eng->incBodyCheckId();
#endif
					continue;
				}

			// Hyperbola needs special checking
			if (isHyperbola &&
			   ((tp<-PIover2 && t>-PIover2) ||
			    (tp< PIover2 && t> PIover2))) {
				STATS(stats.segments_skipped++);
				// Hyperbola on joining branches from - to + or vice versa
				// ignore them
				tp = t;
				xp = x;
				yp = y;
				v->err = 0;
				v->type = SEGMENT_IGNORE;
				v->invalid = false;
				// Is like moving out of window
#ifdef EXPERIMENTAL
				eng->incBodyCheckId();
#endif
				pIn.init();
				pOut.init();
				continue;
			}

			// Check window location of the segment
			// e.g. an arc at the viewport boundary
			double tm = 0.5*(t+tp);
			double xm, ym;
			body->C[i].getXY(tm,&xm,&ym);
			if (!view.inside(xm,ym)) {
				STATS(stats.segments_skipped++);
				tp = t;
				xp = x;
				yp = y;
				v->err = 0;
				v->type = SEGMENT_IGNORE;
				v->invalid = false;
				// outside window increase id
#ifdef EXPERIMENTAL
				eng->incBodyCheckId();
#endif
				pIn.init();
				pOut.init();
				continue;
			}

			// If segment is valid ignore it
			// FIXME needs checking sometimes creates problems
			if (!invalid) { // reset flag
				STATS(stats.segments_invalid++);
				pIn.init();
				pOut.init();
				if (v->invalid) invalid = true;
				tp = t;
				xp = x;
				yp = y;
#ifdef EXPERIMENTAL
				eng->incBodyCheckId();
#endif
				continue;
			}

			// Find search direction
			scanDirection(body, i, xm, ym, &pIn, &pOut);

			// --- Find inside point ---
			pIn.to = j;
			pIn.copy2Old();

			// XXX XXX we should cache the information on the bodies
			// and reset the checkId only on the bodies that are crossing the Vertex2D
			DUMP(cout <<"<+< "<<body->name()<<" vertex="<<j<<" xm="<< xm << " ym="<<ym<<endl);
			// check also against the body
#ifdef EXPERIMENTAL
			eng->getBody(body)->resetCheckId();
#else
			eng->incBodyCheckId();
#endif
			bool bodyIn = eng->getBody(body)->inside2D(
#ifdef EXPERIMENTAL
					eng,
#endif
					pIn.x,pIn.y,pIn.z, pIn.dx,pIn.dy,pIn.dz);

			// and inside the geometry
			if (referenced) eng->where2D(&pIn);
			pIn.from = pIn.to;
			v->zone  = pIn.n? pIn.zone[0] : NULL;

			DUMP(cout << "<+< "<<body->name()<<" inside="<<bodyIn<<endl);

			// --- Find outside point ---
			pOut.to = j;
			pOut.copy2Old();

			DUMP(cout << "<-< "<<body->name()<<" vertex="<<j<<" xm="<<xm<<" ym="<<ym<<endl);
			// check against the body
#ifdef EXPERIMENTAL
			eng->getBody(body)->resetCheckId();
#else
			eng->incBodyCheckId();
#endif
			bool bodyOut = eng->getBody(body)->inside2D(
#ifdef EXPERIMENTAL
					eng,
#endif
					pOut.x,pOut.y,pOut.z, pOut.dx,pOut.dy,pOut.dz);
			// and inside the geometry
			if (referenced) eng->where2D(&pOut);
			pOut.from = pOut.to;

#ifdef EXPERIMENTAL
			// Invalidate last crossed bodies
			if (v->body!=NULL && v->body!=body)
				eng->getBody(v->body)->resetCheckId();
			else
				eng->incBodyCheckId();
#endif

#ifdef EXPERIMENTAL
			// Invalidate current body
			eng->getBody(body)->resetCheckId();
#endif

#ifdef _DUMP
			// STRONG DEBUGGING....
			if (body->show() && bodyIn==bodyOut) {
				cout << endl << "BREAK: " << body->name() << " j=" << j << " i="
					<< i << " in=" << bodyIn << " out=" << bodyOut << endl;
//					body->QN[body->c2q[i]].grad(xm,ym,&dx,&dy,&dz);
//					dd = sqrt(dx*dx+dy*dy+dz*dz);
				cout << "\ttm="<< tm << " xm="<< xm << " ym=" << ym << endl;
//					cout << "\t   : dx="<< dx << " dy=" << dy << " dz=" << dz << " Abs()<=" << (Abs(dz)<SMALL) << endl;
//					cout << "QN[i] = " << body->QN[body->c2q[i]] << endl;
				cout << "C[i]  = " << body->C[i] << endl;
//					cout << "\td= " << dx << ", " << dy << ", " << dz << "\t|d|=" << dd << endl;
				cout << "\tacc[i]= " << body->acc[i] << "  *SMALL=" << body->acc[i]*SMALL << endl;
//					cout << "\tQ(xm,ym,0)=" << body->QN[body->c2q[i]](xm,ym,0.) << endl;
//					cout << "\taccQ(xm,ym,0)=" << body->QN[body->c2q[i]].acc(xm,ym,0.) << endl;
//					cout << "\tQin(xi,yi,zi) =" << body->QN[body->c2q[i]](pIn.x,pIn.y,pIn.z) << endl;
//					cout << "\tQout(xo,yo,zo)=" << body->QN[body->c2q[i]](pOut.x,pOut.y,pOut.z) << endl;
//					cout << "\taccQin(xi,yi,zi) =" << body->QN[body->c2q[i]].acc(pIn.x,pIn.y,pIn.z) << endl;
//					cout << "\taccQout(xo,yo,zo)=" << body->QN[body->c2q[i]].acc(pOut.x,pOut.y,pOut.z) << endl;
				cout << "\tIn : x="<< pIn.x << " y=" << pIn.y << " z=" << pIn.z << endl;
				cout << "\t     dx="<< pIn.dx << " dy=" << pIn.dy << " dz=" << pIn.dz << endl;
				cout << "\tOut: x="<< pOut.x << " y=" << pOut.y << " z=" << pOut.z << endl;
				cout << "\t     dx="<< pOut.dx << " dy=" << pOut.dy << " dz=" << pOut.dz << endl;
			}
#endif
			DUMP(cout << "<-< "<<body->name()<<" outside="<<bodyOut<<endl);

#if 0
			if (v->body==NULL || &body==v->body) {
				pIn.zero();
				pOut.zero();
			}
#endif

#if 0
			if (body->show()) {
			cout << "*-* t=" << tm
				<< " name=" << (v->body?v->body->name():"-null-") << " In=" << endl;
			for (int ii=0; ii<pIn.n; ii++)
				cout << pIn.zone[ii]->region->name() <<":" << *(pIn.zone[ii]) << endl;
			cout << "    Out=" << endl;
			for (int ii=0; ii<pOut.n; ii++)
				cout << pOut.zone[ii]->region->name() <<":" << *(pOut.zone[ii]) << endl;
			}
#endif
			// Compare the regions...
			if (pIn.nRegions == pOut.nRegions) {
				if (pIn.nRegions == 1) {
					v->err = 0;
					if (pIn.zone[0]->region() != pOut.zone[0]->region())
						v->type = SEGMENT_REGION;
					else
					if (pIn.zone[0] != pOut.zone[0])
						v->type = SEGMENT_ZONE;
					else
						v->type = SEGMENT_EMPTY;
				} else {
					if (pIn.n == pOut.n &&
					    equalZones(pIn.zone, pOut.zone, pIn.n)) {
						// it is an error but the same on
						// both sides, ignore to show only the
						// borders below
						v->err  = 0;
						v->type = SEGMENT_EMPTY;
					} else {
						if (pIn.n  != pIn.nOld  ||
						    pOut.n != pOut.nOld ||
						    !equalZones(pIn.zone,  pIn.zoneOld,  pIn.n) ||
						    !equalZones(pOut.zone, pOut.zoneOld, pOut.n))
// FIXME  _errors and errorhash should be protected with semaphores!!!!
							_errors++;
						_errorhash += (body->id()+1)*(i+1)*(j+1);
						v->err  = _errors;
						v->type = SEGMENT_ERROR;
					}
				}
#if 0
			} else
			// Compare the zones...
			if (pIn.n == pOut.n) {
				if (pIn.n == 1) {
					v->err = 0;
					if (pIn.zone[0]->region() != pOut.zone[0]->region())
						v->type = SEGMENT_REGION;
					else
					if (pIn.zone[0] != pOut.zone[0])
						v->type = SEGMENT_ZONE;
					else
						v->type = SEGMENT_EMPTY;
				} else {
					if (equalZones(pIn.zone, pOut.zone, pIn.n)) {
						// it is an error but the same on
						// both sides, ignore to show only the
						// borders below
						v->err  = 0;
						v->type = SEGMENT_EMPTY;
					} else {
						if (pIn.n  != pIn.nOld  ||
						    pOut.n != pOut.nOld ||
						    !equalZones(pIn.zone,  pIn.zoneOld,  pIn.n) ||
						    !equalZones(pOut.zone, pOut.zoneOld, pOut.n)) {
							//cout << "error 2" << endl;
							_errors++;
						}
						_errorhash += (body->id()+1)*(i+1)*(j+1);
						v->err  = _errors;
						v->type = SEGMENT_ERROR;
					}
				}
#endif
			} else
			if (( (  pIn.n==1 &&  pIn.zone[0]->region()->type()==REGION_BLACKHOLE ) && pOut.n==0 ) ||
			    ( ( pOut.n==1 && pOut.zone[0]->region()->type()==REGION_BLACKHOLE ) &&  pIn.n==0 )) {
				// if one of the two is zero and the other is BLACKHOLE then ignore
				v->err  = 0;
				v->type = SEGMENT_REGION;
			} else {
				if (pIn.n  != pIn.nOld  ||
				    pOut.n != pOut.nOld ||
				    !equalZones(pIn.zone,  pIn.zoneOld,  pIn.n) ||
				    !equalZones(pOut.zone, pOut.zoneOld, pOut.n))
					_errors++;
				_errorhash += (body->id()+1)*(i+1)*(j+1);
				v->err  = _errors;
				v->type = SEGMENT_ERROR;
			}

			// On body boundary
			if (bodyIn != bodyOut) v->type |= SEGMENT_BODY;

#if 0
			// FIXME doesn't work properly
			// Check&Skip next vertex if pointing on notReferenced body
			for (int j1=j+1; j1<V.size(); j1++) {
				Vertex2D *vn = V[j1];
				if (vn->body && body!=vn->body && vn->body.notReferenced()) {
					// Copy information from previous segment
					t = vn->t;
					x = vn->x;
					y = vn->y;
					vn->zone  = v->zone;
					vn->err  = v->err;
					vn->type = v->type;
					vn->invalid = false;
					j = j1;
				} else
					break;
			}
#endif

			// save position as previous of the next step
			tp = t;
			xp = x;
			yp = y;

			// reset flag
			invalid = all || body->invalid();
			if (v->invalid) {
				invalid = true;
				v->invalid = false;
			}
		}
	}

	// Set valid flag
	body->setValid();
} // scanBodySegments

/** error
 * return error information on err #err, or closest to point (x,y) if err<0
 * @param err	error index (if <0 then search closest to point (x,y)
 * @param x,y	I/O point to check error, and return error position
 * @param pIn, pOut	region of point information
 * @param xmin, xmax	x min and max limits of the error
 * @param ymin, ymax	y min and max limits of the error
 * @return NULL if no error is found otherwise the body containing the error
 */
VBody* GeometryKernel::error(const int err, double *u, double *v,
			ZoneOfPoint *pIn, ZoneOfPoint *pOut,
			double *umin, double *umax,
			double *vmin, double *vmax)
{
	int i=0,j;
	VBody *body = NULL;
	Array<Vertex2D> *V = NULL;
	Vertex2D *vertex;
	int from=-1, to=-1;		// error range

	for (int ib=0; ib<nGeometryBodies(); ib++) {
		VBody* b = getBody(ib);
		for (i=0; i<b->nC; i++) {
			V = &(b->V[i]);
			if (V->size() < 2) continue;
			for (j=1; j<V->size(); j++) {
				vertex = &(*V)[j];
				if (vertex->err == err) {
					if (from < 0) {
						from = j-1;
						*umin = *umax = (*V)[from].x;
						*vmin = *vmax = (*V)[from].y;
					}
					to = j;
					*umin = Min(*umin, vertex->x);
					*umax = Max(*umax, vertex->x);
					*vmin = Min(*vmin, vertex->y);
					*vmax = Max(*vmax, vertex->y);
				} else
				if (from>=0 && vertex->err!=err) {
					body = b;
					goto FOUND;
				}
			}
			if (from>=0) {
				body = b;
				goto FOUND;
			}
		}
	}
FOUND:
	if (from<0) return NULL;

	vertex = &(*V)[from];
	Vertex2D *vn = &(*V)[to];
	double tm = 0.5 * (vertex->t + vn->t);
	body->C[i].getXY(tm, u, v);
	scanDirection(body,i,*u,*v, pIn, pOut);

	// Prepare for search
	pIn->from = pOut->from = from;
	pIn->to   = pOut->to   = to;
	pIn->nOld = pOut->nOld = 0;
	pIn->body = pOut->body = body;

	// Find inside point
	_engine.incBodyCheckId();
	_engine.where2D(pIn);

	// Find outside point
	_engine.incBodyCheckId();
	_engine.where2D(pOut);

	return body;
} // error

/* updateRegionLocation */
void GeometryKernel::updateRegionLocation(bool all)
{
	for (int i=0; i<nGeometryRegions(); i++) {
		VRegion* region = getRegion(i);
		if (all || region->invalid()) {
			STATS(stats.invalid_regions++);
			DUMP(if (!all) cout << "Found invalid region: " << region->name() << endl);
			region->updateLocation();
			region->setValid();
		}
#ifdef _STAT
		else
			STATS(stats.valid_regions++);
#endif
	}
} // updateRegionLocation

//#define DUMP(x) x
/** correctOverlaps, remove overlaps with other regions from that zone
 * @param region	region to remove overlaps from
 * @param zoneid	zone from region to remove overlaps from (-1=ALL)
 */
bool GeometryKernel::correctOverlaps(VRegion* region, int zoneid)
{
	bool regionModified = false;
	int  side;

	DUMP(cout << endl << "<<< " << *region->region() << endl);

	// Remember initial number of zones
	int nzones = region->nzones();
	while (true) {
		GZone* zone2 = NULL;	// additional zone
		VBody* body = NULL;
		side = SIDE_NONE;

		for (int i=0; i<region->nzones(); i++) {
			if (zoneid>=0 && i<nzones && i!=zoneid) continue;
			VZone *zone = region->zones()[i];
			DUMP(cout << endl << "Zone: " << *zone->zone() << endl);
			// FIXME don't check for bodies already been checked
			// keep id and start from that
			body = zoneOverlaps(zone, &side);
			DUMP(cout << "    SIDE= " << SIDE_STR(side) << endl);

			// Do we need both sides?
			if (body==NULL && (side&SIDE_INOUT)) {
				DUMP(cout << "   *********** We have to add another body..." << endl);
				// I should guess the remaining bodies from the zones
				// of the region that we have the problem
				side = SIDE_NONE;
			} else
			if ((side & SIDE_INOUT) == SIDE_INOUT) {
				DUMP(cout << "    ADD BOTH +/-" << body->name() << endl);
				zone2 = region->region()->addZone(zone->zone());

				zone->zone()->addPlus(body->body());
				zone->zone()->optimize();

				zone2->addMinus(body->body());
				zone2->optimize();
				zone2 = NULL;	// do not scan to remove it later
			} else
			// Only the inside
			if (side & SIDE_IN) {
				DUMP(cout << "    ADD +" << body->name() << endl);
				zone2 = region->region()->addZone(zone->zone());

				zone->zone()->addPlus(body->body());
				zone->zone()->optimize();

				zone2->addMinus(body->body());	// add the other side to test
				zone2->optimize();
			} else
			// Only the outside
			if (side & SIDE_OUT) {
				DUMP(cout << "    ADD -" << body->name() << endl);
				zone2 = region->region()->addZone(zone->zone());

				zone->zone()->addMinus(body->body());
				zone->zone()->optimize();

				zone2->addPlus(body->body());	// add the other side to test
				zone2->optimize();
			}

			// break
			if (side & SIDE_INOUT) break;
		}

		if (!(side & SIDE_INOUT)) break;
		regionModified = true;
		region->init(region->region(), this);
		derive();

		body->setInvalid();
		scanBodySegments(false);

		// check if we have to keep the temporary zone2
		if (zone2) {
			int side2;
			// Find vzone2
			VZone* vzone2=NULL;
			for (int i=0; i<region->nzones(); i++)
				if (region->zone(i)->zone() == zone2) {
					vzone2 = region->zone(i);
					break;
				}
			assert(vzone2!=NULL);
			DUMP(cout << "    Test: " << *zone2 << endl);
			VBody* body2 = zoneOverlaps(vzone2, &side2);
			DUMP(cout << "    SIDE= " << SIDE_STR(side2) << endl);
			if (body2==NULL || side2==SIDE_NONE) {
				DUMP(cout << "*** WRONG SIDE DELETE ***" << endl);
				// Delete last zone
				region->region()->delZone(region->nzones()-1);
				region->init(region->region(), this);
				derive();
				body->setInvalid();
				scanBodySegments(false);
				DUMP(cout << "*-* " << *region->region() << endl);
			} else {
				DUMP(cout << "*** OK KEEP ***" << endl);
			}
		}
	}

	DUMP(cout << ">>> " << *region->region() << endl);
	return regionModified;
} // correctOverlaps

/** isZoneValid *
bool GeometryKernel::isZoneValid(VZone* zone, VBody* body)
{
	return false;
} * isZoneValid */

/** zoneOverlaps
 * @param zone	zone to scan
 * @param side	return which side of body to to keep
 * @return body	that overlaps or appears on the border of the zone
 */
VBody* GeometryKernel::zoneOverlaps(VZone* zone, int* side)
{
	*side = SIDE_NONE;
	int side2 = SIDE_NONE;
	for (int ib=0; ib<nGeometryBodies(); ib++) {
		VBody* body = getBody(ib);
		if (body->notReferenced()) continue;
		//if (body->body()->hasZone(zone->zone())) continue;
		int s = bodySide(zone, body);
		DUMP(cout << "  Body: " << body->name() << "\tside= " << SIDE_STR(s) << endl);
		if (body->body()->hasZone(zone->zone()))
			side2 |= s;
		else {
			*side |= s;
			if (*side & SIDE_INOUT) return body;
		}
	}
	*side |= side2;
	return NULL;
} // zoneOverlaps

/** bodySide, remove overlaps with other regions from that zone
 * @param zone	zone to check for body if belongs to border
 * @param body	body to check
 * @return sides of body found
 */
int GeometryKernel::bodySide(VZone* zone, VBody* body)
{
	/* Scan all errors to find problematic bodies */
//	const int mask = SEGMENT_EMPTY | SEGMENT_REGION | REGION_ERROR;
	int side = SIDE_NONE;

	for (int i=0; i<body->nC; i++) {
		Array<Vertex2D>* V = &(body->V[i]);
		if (V->size() < 2) continue;
		for (int j=1; j<V->size(); j++) {
			Vertex2D* vertex = &(*V)[j];
			if (vertex->type == SEGMENT_IGNORE) continue;
			side |= segmentSide(zone, body, i, j-1, j);
		}
		// early stop
		if ((side&SIDE_BORDER) == SIDE_BORDER) return side;
	}
	return side;
} // bodySide

/** segmentSide, remove overlaps with this body
 * @param zone	zone to remove overlaps from
 * @param body	body might causing the overlap
 * @param i	conic of body
 * @param from	from vertex
 * @param to	to vertex
 * @return	SIDE_FLAGS
 */
int GeometryKernel::segmentSide(VZone* zone, VBody* body, int cid, int from, int to)
{
	ZoneOfPoint pIn;
	ZoneOfPoint pOut;

	Vertex2D* vf = &(body->V[cid][from]);
	Vertex2D* vt = &(body->V[cid][to]);
	double tm = 0.5 * (vf->t + vt->t);
	double u,v;
	body->C[cid].getXY(tm, &u, &v);

	scanDirection(body,cid,u,v, &pIn, &pOut);

	// Prepare for search
	pIn.from = pOut.from = from;
	pIn.to   = pOut.to   = to;
	pIn.nOld = pOut.nOld = 0;
	pIn.body = pOut.body = body;

	// Find inside point
	engine()->incBodyCheckId();
	engine()->where2D(&pIn);

	// Find outside point
	engine()->incBodyCheckId();
	engine()->where2D(&pOut);

	// To be an overlap:
	//   - zone should exist on both sides
	//   - one side should have only one zone (the current)
	//   - the other side should have more than one zones/regions and the current
	bool inZone  = pIn.hasZone(zone);
	bool outZone = pOut.hasZone(zone);

	// If zone is on neither side
	if (!inZone && !outZone) return SIDE_NONE;	// nothing to do with this zone

	// If zone is only on one side then...
	if ( inZone && !outZone) {
		if (pIn.nRegions ==1)
			return SIDE_OK;		// valid for this zone
		else
			return SIDE_IN;		// overlap in IN
	}
	if (!inZone &&  outZone) {
		if (pOut.nRegions==1)
			return SIDE_OK;		// valid for this zone
		else
			return SIDE_OUT;	// overlap in OUT
	}

	// Both sides contain zone
	if (pIn.nRegions==1 && pOut.nRegions==1) return SIDE_NONE;	// valid for the region of the zone
	if (pIn.nRegions>1  && pOut.nRegions>1 ) return SIDE_NONE;	// overlap from both sides

	// Overlap: one side has more than one regions
	if (pIn.nRegions==1)
		return SIDE_IN;
	else
		return SIDE_OUT;
} // segmentSide

/** Calculate the volume of the current region using (Quasi) Monte Carlo integration
 * @param region	region to find volume of
 * @param samples	number of samples to check
 * @return volume of region
 *
 * Based on the implementation of Chris Theis [1/10/2005]
 */
double GeometryKernel::volume(GRegion* region, int samples, double* vol, double* err)
{
	BBox bbox = region->bbox();
	*err = 1.0;
	if (!bbox.isValid()) return -1.0;

	// get random scrambling values
	int RandomNr1 = engine()->random.integer();
	int RandomNr2 = engine()->random.integer();

	// our sequence samplers require the number to be a power of 2
	if( !isPowerOf2(samples)) samples = roundUpPow2(samples);

	long hits = 0;

	*vol = 0.0;
	*err = 1.0;

#ifdef THREAD
	if (threadengines) {
		VolumeFeeder feeder(threadpool);
		for (int n=32; n<samples; n<<=1) {
			feeder.reset(n, 1000, RandomNr1, RandomNr2, &bbox, region);
			threadpool.execute(&feeder);
			hits = feeder.hits();
			// Calculate temporarily
			*vol = (double)hits / n * bbox.volume();
			*err = 1.0 / sqrt((double)hits);
			cout << n << ":  Volume=" << *vol<< " +/- " << *err << endl;
		}
	} else
#else
	{
		Vector size = bbox.size();
		// sample interior
		for (int i=0; i<samples; i++) {
			// sample coordinates using low-discrepancy series
			double x = bbox.low().x + size.x * VanDerCorput(i, RandomNr1);
			double y = bbox.low().y + size.y * Sobol2(i, RandomNr2);
			double z = bbox.low().z + size.z / (double)samples * i;

			// check the geometry with containment tests
			if (region->inside(x,y,z,0.,0.,1.0)) ++hits;
		}
		//cout << "Linear hits=" << hits << endl;
	}
#endif
//	assert(hits == feeder.hits());

	if (hits) {
		*vol = (double)hits / samples * bbox.volume();
		*err = 1.0 / sqrt((double)hits);
	}
	return *vol;
} // volume

/** bodiesMemory */
size_t GeometryKernel::bodiesMemory() const
{
	size_t bmem = bodies.memory();
	for (int i=0; i<bodies.count(); i++)
		bmem += bodies[i]->memory();
	return bmem;
} // bodiesMemory

/** regionsMemory */
size_t GeometryKernel::regionsMemory() const
{
	size_t rmem = regions.memory();
	for (int i=0; i<regions.count(); i++)
		rmem += regions[i]->memory();
	return rmem;
} // regionsMemory

/** memory */
size_t GeometryKernel::memory() const
{
	size_t mem = sizeof(GeometryKernel)
	       + bodiesMemory()
	       + regionsMemory();
#ifdef THREAD
	mem += (nthreads()+1)*_engine.memory();
#else
	mem += _engine.memory();
#endif
	return mem;
} // memory

/** printMemory */
void GeometryKernel::printMemory() const
{
	cout << endl << "GeometryKernel:" << endl;
	cout << "Memory:"     << endl;
	cout << "\tSelf:\t"   << sizeof(GeometryKernel) << endl;
	cout << "\tBodies:\t" << bodiesMemory() << endl;
	cout << "\tRegion:\t" << regionsMemory() << endl;

#ifdef THREAD
	cout << "\tEngine:\t"	<< (nthreads()+1)
				<< " x " << _engine.memory()
				<< " = " << (nthreads()+1)*_engine.memory() << endl;
#else
	cout << "\tEngine:\t"	<< _engine.memory() << endl;
#endif

	cout << "\tTotal:\t"   << memory() << endl;
} // printMemory

/* ========================== GeometryKernelStats ========================== */
#ifdef _STAT
/** reset projection statistics */
void GeometryKernelStats::reset()
{
	segments_total        = 0;
	segments_skipped      = 0;
	segments_invalid      = 0;

	valid_bodies          = 0;
	valid_regions         = 0;
	invalid_bodies        = 0;
	invalid_regions       = 0;
} // reset

/** operator << */
ostream & operator<<(ostream &os, const GeometryKernelStats &gv)
{
	os << " valid_bodies:          " << gv.valid_bodies       << endl;
	os << " invalid_bodies:        " << gv.invalid_bodies     << endl;
	os << " valid_regions:         " << gv.valid_regions      << endl;
	os << " invalid_regions:       " << gv.invalid_regions    << endl;

	cout << " Segments:" << endl;
	cout << "\tTotal   : " << gv.segments_total     << endl;
	cout << "\tSkipped : " << gv.segments_skipped   << endl;
	cout << "\tInvalid : " << gv.segments_invalid   << endl;
	return os;
} /* operator << */
#endif
