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

#ifndef __GEOMETRY_VIEWER_H
#define __GEOMETRY_VIEWER_H

#include <iosfwd>
#include <ostream>
#include <iostream>

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>

#include "os.h"
#include "geo.h"
#include "array.h"
#include "gbody.h"
#include "vbody.h"
#include "vector.h"
#include "vregion.h"
#include "gregion.h"
#include "matrix4.h"
#include "painter.h"
#include "geometry.h"
#include "viewport.h"

#include "kernel.h"
#include "d2layer.h"
#include "d3layer.h"
#include "imagelayer.h"
#include "voxellayer.h"
#include "exportlayer.h"
#include "userdumplayer.h"
#include "usrbinlayer.h"
#include "latticelayer.h"
#include "palettelayer.h"
#include "decorationlayer.h"

// maximum t step on curved conic
#define MAXSTEP	PI/6.0
#define MINSTEP	0.05

enum {
	BIT_SELECT	= 1,
	BIT_LOCK	= 1<<1,
	BIT_FREEZE	= 1<<2,
	BIT_VISIBLE	= 1<<3,
	BIT_WIREFRAME   = 1<<4,
	BIT_MOVE        = 1<<5,
	BIT_BBOX        = 1<<6
};

enum DrawMask {
	DRAW_CLEAR	= 1,
	DRAW_SEGMENTS	= 1<<1,
	DRAW_REGIONS	= 1<<2,
	DRAW_GRID	= 1<<3,
	DRAW_AXES	= 1<<4,
	DRAW_TITLE	= 1<<5,
	DRAW_STATUS	= 1<<6,
	DRAW_LABELS	= 1<<7,
	DRAW_IMAGE	= 1<<8,
	DRAW_USRBIN	= 1<<9,
	DRAW_3D         = 1<<10,
	DRAW_EDGE       = 1<<11,
	DRAW_VERTICES	= 1<<12,
	DRAW_USERDUMP	= 1<<13,
	DRAW_PALETTE    = 1<<14,
	DRAW_LATTICES   = 1<<15,
	DRAW_VOXEL      = 1<<16,
	DRAW_WIREFRAME  = 1<<17,
	DRAW_BBOX       = 1<<18,
	DRAW_MESSAGE    = 1<<19,

	DRAW_TEST       = 1<<30
};
#define DRAW_LATE	(DRAW_LABELS | DRAW_BACKGROUND | DRAW_USRBIN | \
			 DRAW_3D     | DRAW_VERTICES   | DRAW_USERDUMP)

// Describe current step of the Geometry::project() method
// XXX should be in sync with GeometryViewer.py
enum ProjectionState {
	PROJECTION_NOTHING    = 0,	// No projection information
	PROJECTION_SPAWN      = 1,	// spawn projection (waiting to start)
	PROJECTION_START      = 2,	// start projection
	PROJECTION_CONICS     = 3,	// Transform bodies and prepare conics
	PROJECTION_LOCATION   = 4,	// Calculate the bodies&region location
	PROJECTION_INTERSECT  = 5,	// Intersect all bodies
	PROJECTION_SCAN       = 6,	// Scan segments
	PROJECTION_DRAW       = 7,	// Late drawing
	PROJECTION_FINISHED   = 10	// Projection ended
};

typedef void (*NotifyFunc)(void *arg);

/* ========================== GeometryViewer ========================== */
/** Geometry Viewer class */
class GeometryViewer {
protected:
const   Geometry&	geometry;	/** father geometry object	*/
	GeometryKernel&	kernel;		/** geometry kernel with engines*/
	char		_title[64];	/** view title			*/

public:
	bool	showErrors;		/** flag ON/OFF for errors	*/
	int	textBackgroundLevel;	/** text background		*/

public:
	Painter		painter;	/** painter class		*/

// FIXME review (move to protected?)
public:					// Layers
	DecorationLayer decoration;	/** decoration layer		*/
	PaletteLayer	palette;	/** Color box layer		*/

	D2Layer		d2;		/** 2D layer			*/
	D3Layer		d3;		/** 3D layer			*/
	LatticeLayer	lattice;	/** lattice layer		*/
	VoxelLayer	voxel;		/** voxel layer			*/
	UserDumpLayer	userdump;	/** userdump layer		*/
	UsrbinLayer	usrbin;		/** usrbin layer		*/
	ImageLayer	image;		/** image layer			*/
	ExportLayer	exporter;	/** export layer		*/

public:
	BFont		font;		/** general font		*/

protected:				// Layers
	// Treading variables
	pthread_t	thread;		/** thread handler		*/
	pthread_mutex_t	mutexSpawn;	/** spawn mutex			*/
	pthread_mutex_t	mutexStop;	/** stop thread mutex		*/
	bool		pStop;		/** stop signal			*/

// FIXME move to d2layer
	ProjectionState	pState;		/** projection state		*/
	int		_percent;	/** projection percentage	*/

	VBody*		pBody;		/** last processed body		*/
	NotifyFunc	endThread;	/** end project notification	*/
	void*		endArg;		/** arguments to notify	func	*/
public:
	GeometryViewer(const Geometry& g, GeometryKernel& k);
	~GeometryViewer();

	void	makeLatticeColor();

	/* set/get */
	void	title(const char *t);
const	char*	title()		const	{ return _title; }

	void	percent(int p)		{ _percent = p; }
	int	percent()	const	{ return _percent; }

	/** window */
const	ViewPort& view()	const	{ return kernel.view; }
	ViewPort& view()		{ return kernel.view; }
	void	window(double ex, double ey=0.0);
	void	window(double xmin, double ymin, double xmax, double ymax)
			{ view().window(xmin, ymin, xmax, ymax); }
	void	resize(int w, int h);
	int     width()		const	{ return view().width(); }
	int     height()	const	{ return view().height(); }

	void	backgroundColor(const dword color);

	/* transformations */
	void	matrix(const Matrix4& m);
	void	origin(const double x, const double y, const double z);
	void	origin(double *x, double *y, double *z)	{ view().origin(x,y,z); }
	void	moveViewOriginTo0();

	/* threads */
	int	spawnProject(NotifyFunc func=NULL, void *arg=NULL);
	int	spawnDraw(NotifyFunc func=NULL, void *arg=NULL);
	void	stopThread();

	ProjectionState state()	const	{ return pState; }

	void	lockStop()		{ pthread_mutex_lock(&mutexStop); }
	void	unlockStop()		{ pthread_mutex_unlock(&mutexStop); }
	bool	stop()		const	{ return pStop; }

	/* drawing */
	ProjectionState	draw(int mask=-1, bool late=false);

	/* Info */
	size_t	memory()	const;
	void	printMemory()	const;

protected:
	bool	state(ProjectionState s, int p);
static	void*	runProject(GeometryViewer *gthis);
static	void*	runDraw(GeometryViewer *gthis);

friend class Layer;
friend class D2Layer;
}; // GeometryViewer
#endif
