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

#include "quad.h"
#include "bbox.h"
#include "edge.h"
#include "bmath.h"
#include "conic.h"
#include "gzone.h"
#include "vzone.h"
#include "timer.h"
#include "obbox.h"
#include "vertex2d.h"
#include "viewer.h"
#include "kernel.h"

using namespace std;

#define ERRORMSG	"Errors found"

static const char *projectionMessage[] = {
	/* 0 */ "",
	/* 1 */ "",
	/* 2 */ "Starting",
	/* 3 */ "Creating conics",
	/* 4 */ "Finding bodies and region location",
	/* 5 */ "Intersecting bodies",
	/* 6 */ "Scanning segments",
	/* 7 */ ""
};

/** GeometryViewer */
GeometryViewer::GeometryViewer(const Geometry& g, GeometryKernel& k) :
	geometry(g),
	kernel(k),

	painter(640,480),

	decoration(g, k, *this),
	palette(g, k, *this),

	d2(g, k, *this),
	d3(g, k, *this),

	lattice(g, k, *this),
	voxel(g,   k, *this),
	userdump(g,k, *this),
	usrbin(g,  k, *this),
	image(g,   k, *this),

	exporter(g, k, *this)
{
	title("Viewer");

	showErrors  = true;

	textBackgroundLevel = 200;

	// thread
	pthread_mutex_init(&mutexSpawn,   NULL);
	pthread_mutex_init(&mutexStop,    NULL);
	thread    = 0;
	pStop     = false;
	pState    = PROJECTION_NOTHING;
	_percent  = 0;
	endThread = NULL;

	backgroundColor(geometry.backgroundColor());

	font.load(DEFAULT_FONT);
} // GeometryViewer

/** ~GeometryViewer */
GeometryViewer::~GeometryViewer()
{
	stopThread();
	pthread_mutex_destroy(&mutexSpawn);
	pthread_mutex_destroy(&mutexStop);
} /* ~GeometryViewer */

#if 0
/** makeLatticeColor */
void GeometryViewer::makeLatticeColor()
{
	latticeHashColor = Darken(geometry.latticeColor, latticeLevel);
	voxelHashColor   = Darken(geometry.voxelColor,   latticeLevel);
} // makeLatticeColor
#endif

/** backgroundColor */
void GeometryViewer::backgroundColor(const dword color)
{
	painter.background((color & FLAG_COLORMASK) | FLAG_TRANSPARENT);
} // backgroundColor

/** title */
void GeometryViewer::title(const char *t)
{
	strncpy(_title, t, sizeof(_title));
	_title[sizeof(_title)-1] = 0;
} // title

/** origin */
void GeometryViewer::origin(const double x, const double y, const double z)
{
	stopThread();
	kernel.origin(x,y,z);
	image.matrix();
} // origin

/** moveViewOriginTo0
 * Move view origin to center of "view"
 */
void GeometryViewer::moveViewOriginTo0()
{
	stopThread();
	kernel.moveViewOriginTo0();
	image.matrix();
} // moveViewOriginTo0

/** matrix */
void GeometryViewer::matrix(const Matrix4& m)
{
	stopThread();	// Stop thread before modifying the matrix
	kernel.matrix(m);
	image.matrix();
} // matrix

/** resize */
void GeometryViewer::resize(int w, int h)
{
	if (pState == PROJECTION_DRAW) stopThread();
	pthread_mutex_lock(&mutexSpawn);
	painter.init(w, h);
	kernel.resize(w,h);
	pthread_mutex_unlock(&mutexSpawn);
} // resize

/** set projection state and percentage done
 * @return if thread should stop
 */
bool GeometryViewer::state(ProjectionState s, int p)
{
	bool ps;

	ps = pStop;
	if (ps)
		_percent = 0;
	else {
		pState   = s;
		_percent = p;
	}
	return ps;
} // state

/* run project as a thread */
int GeometryViewer::spawnProject(NotifyFunc func, void *arg)
{
	if (pState == PROJECTION_SPAWN) return 0;
	pthread_mutex_lock(&mutexSpawn);
	DUMP(cout << "GeometryViewer::spawnProject(" << title() << ").stopThread()" << endl);
	stopThread();
	endThread = func;
	endArg    = arg;
	pState    = PROJECTION_SPAWN;
	assert(thread==0);
	int rc = pthread_create(&thread, NULL, (void*(*)(void*))GeometryViewer::runProject, this);
	if (rc) {
		perror("pthread_create");
		d2.project();
		kernel.error("System error spawning in background the projection.\n"
				"Switching to synchronous mode");
#if _DEBUG>1
		cerr << "System error spawning in background the projection." << endl;
		cerr << "Switching to synchronous mode" << endl;
#endif
	}
	pthread_mutex_unlock(&mutexSpawn);
	return rc;
} // spawnProject

/* static entry to run a project */
void* GeometryViewer::runProject(GeometryViewer *gthis)
{
	DUMP(cout << "*-* GeometryViewer::runProject(" << gthis->title() << ")" << endl);
	if (gthis->pStop) return NULL;
	if (gthis->pState != PROJECTION_SPAWN) return NULL;

	gthis->d2.project();

	if (gthis->endThread) gthis->endThread(gthis->endArg);
	return NULL;
} // runProject

/* run draw as a thread */
int GeometryViewer::spawnDraw(NotifyFunc func, void *arg)
{
	if (pState != PROJECTION_FINISHED) return 0;
	pthread_mutex_lock(&mutexSpawn);
	DUMP(cout << "GeometryViewer::spawnDraw(" << title() << ").stopThread()" << endl);
	stopThread();
	endThread = func;
	endArg    = arg;
	assert(thread==0);
	assert(!kernel.isRunning());
	int rc = pthread_create(&thread, NULL, (void*(*)(void*))GeometryViewer::runDraw, this);
	if (rc) {
		perror("pthread_create");
		draw(-1, true);
		kernel.error("System error spawning in background the late drawing.\n"
				"Switching to synchronous mode");
#if _DEBUG>1
		cerr << "System error spawning in background the late drawing." << endl;
		cerr << "Switching to synchronous mode" << endl;
#endif
	}
	pthread_mutex_unlock(&mutexSpawn);
	return rc;
} // spawnDraw

/* static entry to submit late drawing */
void* GeometryViewer::runDraw(GeometryViewer *gthis)
{
	assert(!gthis->kernel.isRunning());
	DUMP(cout << "*-* GeometryViewer::runDraw(" << gthis->title() << ")" << endl);
	ProjectionState oldState = gthis->state();
	gthis->state(PROJECTION_DRAW, 0);
	gthis->draw(-1, true);
	gthis->state(oldState,0);
	if (gthis->endThread) gthis->endThread(gthis->endArg);
	return NULL;
} // runDraw

/* stop projection/drawing thread */
void GeometryViewer::stopThread()
{
	DUMP(cout << "*-* GeometryViewer::stopThread(" << title() << ")" << endl);
	lockStop();
	endThread = NULL;	/* do not send end project cmd	*/

	pStop = true;
	kernel.stop();
	if (thread) {
		pthread_join(thread, NULL);
		thread = 0;
	}
	pStop = false;

	// stopThread has to finish with state either
	//  1. FINISHED
	//  2. NOTHING
	if (pState==PROJECTION_DRAW || pState==PROJECTION_FINISHED)
		pState = PROJECTION_FINISHED;
	else {
		_percent = 0;
		pState  = PROJECTION_NOTHING;
		d2.setInvalid();
	}
	unlockStop();
} // stopThread

/** draw
 * @param mask	of options to draw
 * @return	projection state that was drawn
 */
ProjectionState	GeometryViewer::draw(int mask, bool late)
{
	int  line = 0;
	ProjectionState _state = pState;		// Remember projection state

	DUMP(cout << "*-* GeometryViewer::draw(" << title() << ", " << mask << ", " << late << ")" << endl);

	kernel.errorReset();

#ifdef _STAT
	Timer timer;
	timer.start();
#endif

	// Clean display
	if (mask & DRAW_CLEAR) painter.clear();

	if (_state>=PROJECTION_SCAN && mask&DRAW_SEGMENTS && d2.showBorders)
		d2.draw(painter);
	else
		d2.showWireframe = true;

	if (_state==PROJECTION_FINISHED || _state==PROJECTION_DRAW) {
		if (mask&DRAW_REGIONS && d2.fillRegions)	 d2.drawRegions(painter);
		if (late) {
			// Warning on the priority order
			if (mask&DRAW_LATTICES   && lattice.show)	lattice.draw(painter);
			if (mask&DRAW_VOXEL      && voxel.show)		voxel.draw(painter);
			if (mask&DRAW_USRBIN     && usrbin.show)	usrbin.draw(painter);
			if (mask&DRAW_USERDUMP   && userdump.show)	userdump.draw(painter);
			if (mask&DRAW_3D         && d3.show)		d3.draw(painter);
			if (mask&DRAW_IMAGE      && image.show)		image.draw(painter);
			if (mask&DRAW_REGIONS    && d2.fillRegions)	d2.drawRegionsLate(painter);
			if (mask&DRAW_LABELS     && d2.showLabels)	d2.drawLabels(painter);
			if (mask&DRAW_PALETTE    && palette.show)       palette.draw(painter);

		} else {
			if (mask&DRAW_IMAGE && image.show==SHOW_PROMPT) image.draw(painter);
		}
	} else {
		decoration.drawProgress(painter, percent(), projectionMessage[_state]);
		line++;
	}

	// Special 3D real-time drawing
	if (!late && (_state==PROJECTION_FINISHED || _state==PROJECTION_NOTHING)) {
		// WARNING should not be called during the SPAWN projection is ON
		if (mask&DRAW_3D && d3.show)	d3.draw(painter, true);
	}

	if (_state==PROJECTION_FINISHED || _state==PROJECTION_DRAW) {
		// Draw 3D wireframe of bodies
		if (mask&DRAW_WIREFRAME && d2.showWireframe) d3.drawWireframe(painter);

		// Draw bounding boxes
		if (mask&DRAW_BBOX) {
			if (d2.showBBox) d3.drawBodiesBBox(painter);
			d3.drawRegionsBBox(painter);
		}
	}

	if (_state>=PROJECTION_SCAN && mask&DRAW_VERTICES && d2.showVertex && d2.showBorders)
		d2.drawVertices(painter);

	// Independent from the state
	decoration.draw(painter, mask);

	if (mask&DRAW_STATUS && kernel.errors() && showErrors)
		decoration.drawMessage(painter, ERRORMSG, geometry.errorColor, line++);

	if (mask&DRAW_MESSAGE && geometry.message())
		decoration.drawMessage(painter, geometry.message(), geometry.messageColor, line++);
#ifdef _STAT
	timer.stop();
	cerr << "Time to draw: " << timer << endl;
#endif
	return _state;
} // draw

/** memory */
size_t GeometryViewer::memory() const
{
	size_t mem = sizeof(GeometryViewer)
		+ kernel.memory()
//		+ voxel.memory()
		+ painter.memory();
	return mem;
} // memory

/** printMemory */
void GeometryViewer::printMemory() const
{
	cout << endl << "GeometryViewer:" << endl;
	cout << "Memory:"	<< endl;
	cout << "\tSelf:\t"	<< sizeof(GeometryViewer) << endl;
	cout << "\tVideo:\t"	<< painter.memory() << endl;

	cout << "\tKernel:\t"	<< kernel.memory() << endl;
	cout << "\tTotal:\t"	<< memory() << endl;
} // printMemory

