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

#ifndef __D3_LAYER_H
#define __D3_LAYER_H

#include "geo.h"
#include "voxel.h"
#include "layer.h"
#ifdef THREAD
#	include "threadpool.h"
#endif

class GeometryEngine;
class D3Layer;

#ifdef THREAD
/* ============================== Body3DWorker ============================== */
class Body3DWorker : public ThreadPoolWorker {
public:
	D3Layer*	d3;		/** link to D3layer		*/
	GeometryEngine*	engine;		/** engine assigned to worker	*/
	Ray	ray;			/** ray object			*/
	dword*	ptr;			/** pixmap pointer position	*/
	int	W;			/** width of screen		*/
	int	H;			/** height of screen		*/
	int	ii;			/** pixmap array position	*/
	int	jj;
	int	step;			/** step size			*/
	VZone*	zone;			/** last zone information	*/

public:
	Body3DWorker(): d3(NULL), engine(NULL), zone(NULL) {}
	virtual void operator()();
}; // Body3DWorker

/* ============================== Body3DFeeder ============================== */
class Body3DFeeder : public ThreadPoolFeeder {
public:
	int		nworkers;	/** number of workers allocated	*/
	Body3DWorker*	workers;	/** array of body workers	*/
	D3Layer*	d3;		/** link to D3layer		*/
	Painter*	painter;	/** painter object		*/
	dword*		pline;		/** current line		*/
	dword*		ptr;		/** current pointer		*/
	int		W;		/** width of screen		*/
	int		H;		/** height of screen		*/
	int		step;		/** step size to scan		*/
	int		ii;		/** current i screen position	*/
	int		jj;		/** current j screen position	*/

public:
	Body3DFeeder(ThreadPool& p, D3Layer *d) : ThreadPoolFeeder(p),
		  nworkers(0),
		  workers(NULL),
		  d3(d),
		  painter(NULL)
		  {}

	virtual	~Body3DFeeder()		{ if (workers) delete [] workers; }

	// Initialize workers
	void	allocate();
virtual	void	reset(Painter* p, Ray& ray);
	void	init(int astep);
	bool	loop();

virtual	ThreadPoolWorker* feed(int threadId);
}; // Body3DFeeder
#endif

#ifdef _STAT
/* ============================ D3LayerStats ========================= */
/** Geometry Viewer Statistics */
class D3LayerStats {
public:
	long long shade_call_count;
	long long shade_correctRotdefiNormal_count;

	D3LayerStats() { reset(); }

	void reset();
}; // D3LayerStats
std::ostream& operator << (std::ostream&, const D3LayerStats&);
#endif

/* ================================== D3Layer ================================ */
class D3Layer : public Layer {
public:
	/* 3D settings */
	bool	deflights;		/** default lights		*/
	int	ambient;		/** ambient intensity [0..255]	*/
	int	xray;			/** xray transparency		*/
	int	antialias;		/** super-sampling antialias	*/
	bool	drawEdges;		/** detect and draw edges	*/
	bool	shadows;		/** show shadows		*/
	bool	skip_1stblack;		/** skip 1st blackhole		*/
	VLight	light[MAXLIGHT];	/** rel. lights in the system	*/

protected:
	bool	_usrbinastexture;	/** usrbin mapping as texture	*/
	bool	_maxDepth;		/** enable reflections/refractions*/
	int	_step;			/** pixel step of 3D drawing	*/
	bool	_showErrors;		// temporary copy of showErrors
	bool	_edgeDetect;		// temporary copy of edgeDetect
	int	_samples;		// temporary copy of antialias

	clock_t	_drawTime;		// how long to draw

#ifdef THREAD
	Body3DFeeder	feeder;		/** threadpool body feeder	*/
#endif
#ifdef _STAT
	D3LayerStats stats;             /** Statistics                  */
#endif

public:
	D3Layer(const Geometry& g, GeometryKernel& k, GeometryViewer& v);

	bool	usrbinAsTexture()	const	{ return _usrbinastexture; }
	void	usrbinAsTexture(bool b)		{ _usrbinastexture = b; }

	int	maxDepth()		const	{ return _maxDepth; }
	void	maxDepth(int b)			{ _maxDepth = b; }

	clock_t	drawTime()		const	{ return _drawTime?CLOCKS_PER_SEC/_drawTime:0; }
	void	drawTime(int div)		{ _drawTime = div? CLOCKS_PER_SEC/div : 0; }

	/* draw method */
	void	draw(Painter& painter, bool checkTime=false);

	Color	shade(GeometryEngine* engine, Ray *ray);
	bool	nextIntersection(GeometryEngine* engine, Ray *ray);
	dword	shadeXray(GeometryEngine* engine, Ray* ray,
			const double u, const double v,
			int alpha=255, VRegion* last_region=NULL);

	bool	draw3Dline(Painter& painter,
			const Point& a, const Point& b,
			const Color3D& color);

	void	drawWireframe(Painter& painter);
	void	drawWireframe(Painter& painter, VBody *body);
	void	drawBBox(Painter& painter, const BBox& bbox, const Color3D& color);
	void	drawOBBox(Painter& painter, const OBBox* bbox, const Color3D& color);

	void	drawBodiesBBox(Painter& painter);
	void	drawRegionsBBox(Painter& painter);

protected:
	bool	edgePixel(GeometryEngine* engine, Ray& ray, const double u, const double v, VZone* zone);
	VZone*	drawPixel(GeometryEngine* engine, Ray& ray,
			dword *ptr, int W, int H,
			int i, int j, int step, int sample,
			double u, double v,
			VZone* zone);

friend class Body3DWorker;
}; // D3Layer
#endif
