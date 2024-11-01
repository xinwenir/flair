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
 *
 */

#include "timer.h"
#include "viewer.h"
#include "d2layer.h"
#include "geometry.h"

// maximum line length^2 on curved conic
#define MAXLEN2 Sqr(10.0)

#define SETPIXEL(i,j,c) if (painter.insideClip(i,j)) painter.pixel(i,j,c);

// flags for debugging printouts
#ifdef _DUMP
#	define DRAWVERTEX
#endif

using namespace std;

#ifdef THREAD
/* =============================== BodyFeeder =============================== */
/** allocate */
void BodyFeeder::allocate()
{
	if (nworkers != pool.nthreads()) {
		if (workers) delete [] workers;
		nworkers = pool.nthreads();
		workers = new BodyWorker[nworkers];
		for (int i=0; i<nworkers; i++)
			workers[i].feeder(this);
	}
} // allocate

/** reset */
void BodyFeeder::reset(D2Method method, int pfrom, int pto)
{
	pool.reset();
	allocate();
	idx = 0;
	perc_from = pfrom;
	perc_span = pto - pfrom;

	// initialize workers
	for (int i=0; i<nworkers; i++) {
		workers[i].d2     = d2;
		workers[i].method = method;
		workers[i].engine = d2->kernel.engine(i);
		workers[i].body   = NULL;
	}
} // reset

/** feed */
ThreadPoolWorker* BodyFeeder::feed(int threadId)
{
	// Check for ending condition
	if (idx >= d2->kernel.nGeometryBodies() || d2->viewer.stop()) return NULL;

	// update worker
	workers[threadId].body = d2->kernel.getBody(idx++);

	// Update information done up to now
	d2->viewer.percent(perc_from + (perc_span*idx)/d2->kernel.nGeometryBodies());

	return (ThreadPoolWorker*)&workers[threadId];
} // feed
#endif

/* ============================== D2Layer ============================= */
/** D2Layer */
D2Layer::D2Layer(const Geometry& g, GeometryKernel& k, GeometryViewer& v) :
		Layer(g, k, v),
		labels(32)
#ifdef THREAD
		, feeder(k.threadpool, this)
#endif
{
	_projectAll   = true;
	showBorders   = true;
	showVertex    = true;
	showLabels    = true;
	fillRegions   = true;
	showWireframe = false;
	showBBox      = false;
	_projectionValid = false;
} // D2Layer

/** _intersectBody threaded version of intersect body, called from the threadpool */
void D2Layer::_intersectBody(VBody* a, GeometryEngine*)
{
	kernel.intersectBody(a, _projectAll);
} // _intersectBody

/** project - create the projection on the viewport
 * two main cases can happen depending on the value of _projectAll:
 * - if _projectAll is true: all the bodies have to be processed
 * - if _projectAll is false: only invalid bodies have changed and only those have
 *   to be processed
 */
void D2Layer::project()
{
	if (!showBorders) {
		viewer.state(PROJECTION_FINISHED, 100);
		_projectionValid = false;
		return;
	} else
		_projectionValid = true;
	kernel.lock();

	// XXX FIME XXX
	// I could use some adaptive technique from the previous runs to learn
	// the time needed for the plots to make a better estimate
	// e.g. keep the 2 or 3 last runs information
	// and predict using the average for the next ones

	viewer.state(PROJECTION_START, 0);
	kernel.clearErrors();

#ifdef _STAT
	Timer timer;
	timer.start();
	STATS(kernel.stats.reset());
#endif
	DUMP(cout << endl << "*-* D2Layer::project("<< viewer.title() <<")" << endl);

	/* generate and transform bodies */
	if (viewer.state(PROJECTION_CONICS, 0)) {
		kernel.unlock();
		return;
	}

	// Project conics to viewer plane
	kernel.makeBodyConics(_projectAll);

	if (viewer.state(PROJECTION_LOCATION, 1)) {
		kernel.unlock();
		return;
	}

#ifdef _STAT
	timer.stop();
	cerr << "Time to make body conics: " << timer << endl;
	timer.start();
#endif

	// update region location
	kernel.updateRegionLocation(_projectAll);

#ifdef _STAT
	timer.stop();
	cerr << "Time to update region locations: " << timer << endl;
	timer.start();
#endif

	// Compute all intersections
#ifdef THREAD
	if (kernel.nthreads()) {
		viewer.pState = PROJECTION_INTERSECT;
		feeder.reset(&D2Layer::_intersectBody, 2,20);
		kernel.threadpool.execute(&feeder);
	} else
#endif
	{
		for (int i=0; i<kernel.nGeometryBodies(); i++) {
			VBody* a = kernel.getBody(i);
			if (a->ignore()) continue;

			if (viewer.state(PROJECTION_INTERSECT, 2+(18*i)/kernel.nGeometryBodies())) {
				kernel.unlock();
				return;
			}

			// Generate all body combinations and intersect them
			for (int j=i+1; j<kernel.nGeometryBodies(); j++) {
				VBody* b = kernel.getBody(j);
				if (b->ignore()) continue;
				if (_projectAll || a->invalid() || b->invalid())
					a->intersectBody(b, view());
			}
		}
	}

#ifdef _STAT
	timer.stop();
	cerr << "Time for body intersections: " << timer << endl;
	timer.start();
#endif

	/* Clean up bodies */
	for (int i=0; i<kernel.nGeometryBodies(); i++) {
		VBody* body = kernel.getBody(i);
		body->removeWrongVertices();

		if (_projectAll || body->invalid()) {
// FIXME review: maybe calcAcc should be done for all bodies?
			body->calculateAccuracy(view());
// FIXME maybe this calculation can be moved to the scanSegments
// and remove the notref variable from the VBody!
			if (!body->body()->zones.empty()) {
				bool notused = true;

				CBody *cbody = engine()->getBody(body);

				// Check if all referring zones are not displayed
				for (int iz = 0; iz < cbody->zones.size(); iz++) {
					if (!(cbody->zones[iz])->ignore()) {
						notused = false;
						break;
					}
				}
				if (notused)
					body->notReferenced(true);
			}
		}
	}

#ifdef _STAT
	timer.stop();
	cerr << "Time body cleanup: " << timer << endl;
	timer.start();
#endif

	scanSegments();
	viewer.state(PROJECTION_FINISHED, 100);

#ifdef _STAT
	timer.stop();
	cerr << "Time to process segments: " << timer << endl;

	cerr << "project: all= " << _projectAll << endl;
	cerr << "\t# Bodies:  " << kernel.nGeometryBodies()  << endl;
	cerr << "\t# Regions: " << kernel.nGeometryRegions()  << endl;
	cerr << "Stats\n" << kernel.stats << endl;
	cerr << "Engine Stats\n" << engine()->stats << endl;
#endif

	DUMP(cout << endl << "*-* D2Layer::endproject('"<< viewer.title() <<"')" << endl);
	kernel.unlock();
} // project

/** scanSegments
 * loop through all the segments in bodies and assign them an in and out region
 */
void D2Layer::scanSegments()
{
	DUMP(cout << endl << "*-* D2Layer::scanSegments("<<viewer.title()<<")" << endl);

	viewer.pState = PROJECTION_SCAN;
	kernel.clearErrorHash();

#ifdef THREAD
	if (kernel.nthreads()) {
		if (kernel.nGeometryBodies()) {
			feeder.reset(&D2Layer::_scanBodySegments, 20,100);
			kernel.threadpool.execute(&feeder);
		}
	} else
#endif
	{
		for (int ib=0; ib<kernel.nGeometryBodies(); ib++) {
			kernel.scanBodySegments(kernel.getBody(ib), engine(), _projectAll);
			// Update information done up to now
			viewer.pState = PROJECTION_SCAN;
			if (kernel.nGeometryBodies())
				viewer.percent(20 + (90*ib)/kernel.nGeometryBodies());
			else
				viewer.percent(99);
			if (stop()) break;
		}
	}
} // scanSegments

/** scanBodySegments called from the threadpool */
void D2Layer::_scanBodySegments(VBody* body, GeometryEngine* eng)
{
	kernel.scanBodySegments(body, eng, _projectAll);
} // _scanBodySegments

/** find closest vertex
 * @param u,v	location to search
 * @param d	minimum distance
 * @param vu, vv	return vertex u,v location
 * @return true	if vertex is found
 * */
bool D2Layer::closestVertex(const double u, const double v, const double d,
				double *vu, double *vv)
{
	// Can only check if the projection exists
	if (viewer.state()!=PROJECTION_FINISHED && viewer.state()!=PROJECTION_DRAW)
		return false;

	double d2 = d*d;
	kernel.lock();	// <--- Locking Geometry?
	for (int ib=0; ib<kernel.nGeometryBodies(); ib++) {
		VBody* body = kernel.getBody(ib);
		if (!body->show() || body->show()&BIT_MOVE) continue;

		// Check body center on UVW system
		Point p = view().xyz2uvw(body->body()->position());
		if (Eq0(p.z,SMALL)) {
			if (Sqr(p.x-u) + Sqr(p.y-v) <= d2) {
				*vu = p.x;
				*vv = p.y;
				kernel.unlock();
				return true;
			}
		}

		// Check conic vertices
		for (int i=0; i < body->nC; i++) {
			if (body->V[i].size()<2) continue;

			Array<Vertex2D>& V = body->V[i];

			int kprev = -1;
			for (int k=0; k<V.size(); k++) {
				Vertex2D& vertex = V[k];
				if (vertex.body && vertex.body->show()) {
					if (Sqr(vertex.x-u) + Sqr(vertex.y-v) <= d2) {
						*vu = vertex.x;
						*vv = vertex.y;
						kernel.unlock();
						return true;
					}
					if (body == vertex.body) {
						if (kprev<0)
							kprev = k;
						else {
							double tm = 0.5 * (vertex.t + V[kprev].t);
							double xm, ym;
							body->C[i].getXY(tm,&xm,&ym);
							if (Sqr(xm - u) + Sqr(ym - v) <= d2) {
								*vu = xm;
								*vv = ym;
								kernel.unlock();
								return true;
							}
						}
					}
				}
			}
		}
	}
	kernel.unlock();
	return false;
} // closestVertex

/** projectVertices
 * project vertices along an axis an populate the array
 */
bool D2Layer::projectVertices(const char axis, Array<double>& vertices)
{
	// Can only check if the projection exists
	if (viewer.state()!=PROJECTION_FINISHED && viewer.state()!=PROJECTION_DRAW)
		return false;

	vertices.compare(Cmp);

	kernel.lock();	// <--- Locking Geometry?
	for (int ib=0; ib<kernel.nGeometryBodies(); ib++) {
		VBody* body = kernel.getBody(ib);
		// Check conic vertices
		for (int i=0; i < body->nC; i++) {
			if (body->V[i].size()<2) continue;
			Array<Vertex2D>& V = body->V[i];
			for (int k=0; k<V.size(); k++) {
				Vertex2D& vertex = V[k];
				switch (axis) {
					case 'x':
					case 'X':
						vertices.add(view().uv2x(vertex.x,vertex.y));
						break;
					case 'y':
					case 'Y':
						vertices.add(view().uv2y(vertex.x,vertex.y));
						break;
					case 'z':
					case 'Z':
						vertices.add(view().uv2z(vertex.x,vertex.y));
						break;
					case 'u':
					case 'U':
						vertices.add(vertex.x);
						break;
					case 'v':
					case 'V':
						vertices.add(vertex.y);
						break;
				}
			}
		}
	}
	kernel.unlock();
	return true;
} // projectVertices

/** drawSegment
 * @param body	body of which to draw the segment
 * @param ic	index of conic to draw
 * @param from	starting index in the vertex array
 * @param to	ending index (+1) in the vertex array
 * @param color	color to use for drawing
 * @param step	I/O step used in drawing (<0.0 to reset)
 * @return true is segment is visible (inside clipping area)
 */
bool D2Layer::drawSegment(Painter& painter, VBody *body, const int ic,
			const int from,    const int to,
			const dword color, const int width,
			double *step)
{
	double x, y;
	int x1,y1, x2,y2;
	bool visible=false;
	const Conic& conic = body->C[ic];
	const Array<Vertex2D>& V = body->V[ic];
	const Vertex2D *v = &V[from];
	const double ts = v->t;
	const double xs = v->x;
	const double ys = v->y;
	v = &V[to];
	const double te = v->t;
	const double xe = v->x;
	const double ye = v->y;

	if (te-ts < SMALL2) return false;
	double maxstep = (conic.type()==CONIC_PARABOLA)? (te-ts)/10.0 : MAXSTEP;

#ifdef DRAWVERTEX
	// Draw each segment
	for (int i=from+1; i<to-1; i++) {
		v = &V[i];
		SETPIXEL(view().u2i(v->x), view().v2j(v->y), 0xFFFFFF);
	}
	// Draw mid point
	double tm = 0.5*(ts+te);
	double xm,ym;
	conic.getXY(tm,&xm,&ym);
	x1 = view().u2i(xm);
	y1 = view().v2j(ym);
	SETPIXEL(x1,  y1,  0xFFFFFF);
	SETPIXEL(x1-1,y1-1,0xFFFFFF);
	SETPIXEL(x1+1,y1-1,0xFFFFFF);
	SETPIXEL(x1-1,y1+1,0xFFFFFF);
	SETPIXEL(x1+1,y1+1,0xFFFFFF);
#endif

	if (ISZERO(Sqr(xe-xs) + Sqr(ye-ys)) && ISEQ(ts,te))
		return false;

	// draw lines
	x1 = view().u2i(xs);
	y1 = view().v2j(ys);
	if (conic.type() == CONIC_LINE) {
		x2 = view().u2i(xe);
		y2 = view().v2j(ye);
		if (!width)
			visible |= painter.line(x1,y1,x2,y2,color);
		else
			visible |= painter.lineThick(x1,y1,x2,y2,width,color);
	} else {
		// initial guess
		if (*step<=0.0) *step = maxstep;
		if (*step > te-ts) *step = te-ts;
		double t1 = ts;
		do {
			double t2 = t1 + *step;
			if (t2>te) {
				t2 = te;
				*step = te - ts;
			}
			conic.getXY(t2, &x, &y);
			x2 = view().u2i(x);
			y2 = view().v2j(y);

			// find at half way what is the distance from a line
			while (1) {
				if (Abs(x2-x1)<=2 && Abs(y2-y1)<=2) {
					if (!width)
						visible |= painter.line(x1,y1,x2,y2,color);
					else
						visible |= painter.lineThick(x1,y1,x2,y2,width,color);
#ifdef DRAWVERTEX
					SETPIXEL(x1,y1,0xFFFFFF);
#endif
					x1 = x2;
					y1 = y2;
					t1 = t2;
					*step *= 2.0;
					if (*step > maxstep) *step = maxstep;
					break;
				}

				double xh, yh;
				double t3 = t1 + *step/2.0;
				assert(t3<=te);
				conic.getXY(t3, &xh, &yh);
				int x3 = view().u2i(xh);
				int y3 = view().v2j(yh);

				// distance^2 of line 1,2
				double l2 = Sqr((double)(x2-x1)) + Sqr((double)(y2-y1));
				// perpendicular distance of midpoint 3 wrt to line 1,2
				double d2 = Sqr((double)(x3-x1)*(double)(y2-y1)
					      - (double)(y3-y1)*(double)(x2-x1))
					      / (double)l2;
				if (d2<=0.01 && l2<MAXLEN2) {
					if (!width)
						visible |= painter.line(x1,y1,x2,y2,color);
					else
						visible |= painter.lineThick(x1,y1,x2,y2,width,color);
#ifdef DRAWVERTEX
					SETPIXEL(x1,y1,0xFFFFFF);
#endif
					x1 = x2;
					y1 = y2;
					t1 = t2;
					*step *= 2.0;
					if (*step > maxstep) *step = maxstep;
					break;
				} else
				if (d2<=0.1 || l2<MAXLEN2) {
					if (!width)
						visible |= painter.line(x1,y1,x2,y2,color);
					else
						visible |= painter.lineThick(x1,y1,x2,y2,width,color);
#ifdef DRAWVERTEX
					SETPIXEL(x1,y1,0xFFFFFF);
#endif
					x1 = x2;
					y1 = y2;
					t1 = t2;
					break;
				} else
				if (Abs(x3-x1)<=2 && Abs(y3-y1)<=2) {
					if (!width)
						visible |= painter.line(x1,y1,x2,y2,color);
					else
						visible |= painter.lineThick(x1,y1,x2,y2,width,color);
#ifdef DRAWVERTEX
					SETPIXEL(x1,y1,0xFFFFFF);
#endif
					x1 = x3;
					y1 = y3;
					t1 = t3;
					*step *= 2.0;
					break;
				} else {
					x2 = x3;	// half
					y2 = y3;
					t2 = t3;
					*step /= 2.0;
				}
			}
		} while (t1 < te - *step);
		x2 = view().u2i(xe);
		y2 = view().v2j(ye);
		if (!width)
			visible |= painter.line(x1,y1,x2,y2,color);
		else
			visible |= painter.lineThick(x1,y1,x2,y2,width,color);
	}
#ifdef DRAWVERTEX
	x1 = view().u2i(xs);
	y1 = view().v2j(ys);
	SETPIXEL(x1,y1, 0xFFFF00);
	SETPIXEL(x1-1,y1-1, 0xFFFF00);
	SETPIXEL(x1+1,y1-1, 0xFFFF00);
	SETPIXEL(x1-1,y1+1, 0xFFFF00);
	SETPIXEL(x1+1,y1+1, 0xFFFF00);

	SETPIXEL(x2,y2, 0xFFFF00);
	SETPIXEL(x2-1,y2-1, 0xFFFF00);
	SETPIXEL(x2+1,y2-1, 0xFFFF00);
	SETPIXEL(x2-1,y2+1, 0xFFFF00);
	SETPIXEL(x2+1,y2+1, 0xFFFF00);
#endif
	return visible;
} // drawSegment

/** draw geometry to view */
void D2Layer::drawSegments(Painter& painter, VBody *body, dword mask)
{
	dword	maskZone;
	maskZone = mask | SEGMENT_ZONE;

	for (int i=0; i < body->nC; i++) {
		Array<Vertex2D>& V = body->V[i];
		if (V.size()<2) continue;

		// Initial vertex
		int kp = 0;
		double step = -1.0;

		// Next vertex
		Vertex2D *v = &V[1];
		dword type;
		// WARNING the order is important to avoid segfault when v->in==NULL
		bool zone = mask!=SEGMENT_BODY &&
		            v->type&SEGMENT_ZONE &&
			    v->zone->region()->show()&BIT_SELECT;
		if (zone)
			type = v->type & maskZone;
		else
			type = v->type & mask;

		for (int k=2; k<V.size(); k++) {
			// Check next vertex
			v = &V[k];

			if (zone == ( mask!=SEGMENT_BODY &&
					  v->type&SEGMENT_ZONE &&
					  v->zone->region()->show()&BIT_SELECT) &&
				    (v->type&mask)==type) {
#ifdef DRAWVERTEX
				int x = view().u2i(v->x);
				int y = view().v2j(v->y);
				SETPIXEL(x,y,0x00FFFF);
#endif
				continue;
			}

			if (zone || (type&mask)) {
				dword color;
				int   width=0;

				if (mask==SEGMENT_BODY) {
					color = body->color();
					width = body->lineWidth();
				} else
				if (type&SEGMENT_ERROR)
					color = viewer.showErrors?
							geometry.errorColor : geometry.regionColor;
				else
				if (type&SEGMENT_ZONE)
					color = geometry.zoneColor;
				else
				if (type&SEGMENT_REGION)
					color = geometry.regionColor;
				else {
					color = body->color();
					width = body->lineWidth();
				}

				if (drawSegment(painter, body, i, kp, k-1, color, width, &step))
					body->visible = true;
			}
#ifdef DRAWVERTEX
			else {
				int x1 = view().u2i(V[kp].x);
				int y1 = view().v2j(V[kp].y);
				SETPIXEL(x1,  y1,  0x00FFFF);
				SETPIXEL(x1-1,y1-1,0x00FFFF);
				SETPIXEL(x1+1,y1-1,0x00FFFF);
				SETPIXEL(x1-1,y1+1,0x00FFFF);
				SETPIXEL(x1+1,y1+1,0x00FFFF);
				x1 = view().u2i(v->x);
				y1 = view().v2j(v->y);
				SETPIXEL(x1,  y1,  0x00FFFF);
				SETPIXEL(x1-1,y1-1,0x00FFFF);
				SETPIXEL(x1+1,y1-1,0x00FFFF);
				SETPIXEL(x1-1,y1+1,0x00FFFF);
				SETPIXEL(x1+1,y1+1,0x00FFFF);
			}
#endif

			kp = k-1;
			zone = mask!=SEGMENT_BODY &&
				v->type&SEGMENT_ZONE &&
				v->zone->region()->show()&BIT_SELECT;
			if (zone)
				type = v->type & maskZone;
			else
				type = v->type & mask;
		}

		if (zone || (type&mask)) {
			dword color;
			int   width = 0;

			if (mask==SEGMENT_BODY) {
				color = body->color();
				width = body->lineWidth();
			} else
			if (type&SEGMENT_ERROR)
				color = viewer.showErrors?
						geometry.errorColor : geometry.regionColor;
			else
			if (type&SEGMENT_ZONE)
				color = geometry.zoneColor;
			else
			if (type&SEGMENT_REGION)
				color = geometry.regionColor;
			else {
				color = body->color();
				width = body->lineWidth();
			}

			if (drawSegment(painter, body, i, kp, V.size()-1, color, width, &step))
				body->visible = true;
		}
#ifdef DRAWVERTEX
		else {
			int x1 = view().u2i(V[kp].x);
			int y1 = view().v2j(V[kp].y);
			SETPIXEL(x1,y1,0x00FFFF);
			SETPIXEL(x1-1,y1-1,0x00FFFF);
			SETPIXEL(x1+1,y1-1,0x00FFFF);
			SETPIXEL(x1-1,y1+1,0x00FFFF);
			SETPIXEL(x1+1,y1+1,0x00FFFF);
			x1 = view().u2i(v->x);
			y1 = view().v2j(v->y);
			SETPIXEL(x1,y1,0x00FFFF);
			SETPIXEL(x1-1,y1-1,0x00FFFF);
			SETPIXEL(x1+1,y1-1,0x00FFFF);
			SETPIXEL(x1-1,y1+1,0x00FFFF);
			SETPIXEL(x1+1,y1+1,0x00FFFF);
		}
#endif
	}
} // drawSegments

/** drawNodes */
void D2Layer::drawNodes(Painter& painter, VBody *body)
{
	const GBody* gbody = body->body();

	for (int n=0; n<gbody->nodes(); n++) {
		double u, v, w;
		view().xyz2uvw(gbody->node(n), &u, &v, &w);
		if (n==0 || Eq0(w,SMALL3)) {
			int x = view().u2i(u);
			int y = view().v2j(v);
			if (painter.rectangle(x-2,y-2,x+2,y+2, geometry.selectColor))
				body->visible = true;
		}
	}
} // drawNodes

/** draw vertices */
void D2Layer::drawVertices(Painter& painter)
{
	for (int ib=0; ib<kernel.nGeometryBodies(); ib++) {
		VBody* body = kernel.getBody(ib);
		if (!body->show()) continue;

		for (int i=0; i < body->nC; i++) {
			if (body->V[i].size()<2) continue;
			Array<Vertex2D>& V = body->V[i];

			int kprev = -1;
			for (int k=0; k<V.size(); k++) {
				Vertex2D& v = V[k];
				if (v.body && v.body->show()) {
					// XXX shows start/end of circles!!!
					int x = view().u2i(v.x);
					int y = view().v2j(v.y);
					SETPIXEL(x,y,  geometry.vertexColor);
					SETPIXEL(x-1,y,geometry.vertexColor);
					SETPIXEL(x+1,y,geometry.vertexColor);
					SETPIXEL(x,y-1,geometry.vertexColor);
					SETPIXEL(x,y+1,geometry.vertexColor);

					// Draw also middle of segment (between two visible vertices!
					// one pixel is enough
					if (body == v.body) {
						if (kprev<0)
							kprev = k;
						else {
							double tm = 0.5 * (v.t + V[kprev].t);
							double xm, ym;
							body->C[i].getXY(tm,&xm,&ym);
							x = view().u2i(xm);
							y = view().v2j(ym);
							SETPIXEL(x,y,geometry.vertexColor);
						}
					}
				}
			}
		}
	}
} // drawVertices

/** draw regions */
void D2Layer::drawRegions(Painter& painter)
{
	int W = painter.width();
	int widthadd = W - (painter.clip().right - painter.clip().left)-1;
	dword *pixelPtr = painter.pixelPtr(painter.clip().left, painter.clip().top);

	double dx = -view().matrix(0,2);
	double dy = -view().matrix(1,2);
	double dz = -view().matrix(2,2);

	geometry.lockRead();	// No need to lock the kernel, the segments are there
	for (int j=painter.clip().top; j<=painter.clip().bottom; j++) {
		if (stop()) break;
		//assert(pixelPtr == painter.pixelPtr(painter._clip.left, j));
		for (int i=painter.clip().left; i<=painter.clip().right; i++, pixelPtr++) {
			// Find a pixel that has at least on background pixel
			// in every direction
			//     b b b
			//     b + b
			//     b b b
			// we start with the up-left corner then check for i+1, j+1
			// 2x2 sometimes it fails
			if (*pixelPtr       != painter.background() ||
			    pixelPtr[1]     != painter.background() ||
			    pixelPtr[2]     != painter.background() ||
			    pixelPtr[W]     != painter.background() ||
			    pixelPtr[W+1]   != painter.background() ||
			    pixelPtr[W+2]   != painter.background() ||
			    pixelPtr[2*W]   != painter.background() ||
			    pixelPtr[2*W+1] != painter.background() ||
			    pixelPtr[2*W+2] != painter.background()) continue;
			double u = view().ic2u(i+1);
			double v = view().jc2v(j+1);

			double x = view().uv2x(u,v);
			double y = view().uv2y(u,v);
			double z = view().uv2z(u,v);

			engine()->incBodyCheckId();

			geometry.lockEdit();
			// Draw selected zones
			if (geometry.editRegion().nzones()) {
				GZone *zone = geometry.editRegion().inside(x,y,z, dx,dy,dz);
				if (zone) {
					painter.fill(i, j,
						0xFFFFFF,
						(((const GZone*)zone==geometry.editRegion().zones().tail())?
							geometry.zoneColor : geometry.visibleColor),
						FILL_X);
					geometry.unlockEdit();
					continue;
				}
			}
			geometry.unlockEdit();

			VZone *zone = engine()->where2D(x,y,z, dx,dy,dz);
			if (zone) {
				const VRegion *region = zone->region();
				if (region->show() & BIT_SELECT) {
					Color32 col;
					col.val = region->color();
					col.rgb.red   = col.rgb.red  >0x7F? 0:0xFF;
					col.rgb.green = col.rgb.green>0x7F? 0:0xFF;
					col.rgb.blue  = col.rgb.blue >0x7F? 0:0xFF;
					painter.fill(i, j, region->color(),
						col.val,
						//(region->color()^geometry.selectColor),
						(FillType)(region->show()&BIT_SELECT));
				} else
				if (region->type() == REGION_LATTICE)
					painter.fill(i, j, geometry.latticeColor,
						viewer.lattice.hashColor, FILL_HASHR);
				else
				if (region->type() == REGION_VOXEL)
					painter.fill(i, j, geometry.voxelColor,
						viewer.voxel.hashColor, FILL_HASHR);
				else
				if (viewer.d3.show && region->transparent())
					painter.fill(i, j, geometry.transparentColor());
				else
				if (geometry.lighterLevel)
					painter.fill(i,j, Lighten(region->color(),
						geometry.lighterLevel));
				else
					painter.fill(i,j, region->color());
			} else	{ // Error
				if (viewer.showErrors)
					painter.fill(i,j,
						(viewer.d3.show?
							geometry.transparentColor():
							0xFFFFFF | FLAG_ERROR),
						geometry.regionErrorColor,
						FILL_HASH);
				else
					painter.fill(i,j,
						viewer.d3.show?
							geometry.transparentColor():
							0xFFFFFF | FLAG_ERROR);
			}
		}
		pixelPtr += widthadd;
	}
	geometry.unlockRead();
} // drawRegions

/** late draw regions */
void D2Layer::drawRegionsLate(Painter& painter)
{
	int H = painter.height();	// save variables that might change
	int W = painter.width();	// during the scanning...
	dword *ptr = painter.data();	// data pointer

	double dx = -view().matrix(0,2);
	double dy = -view().matrix(1,2);
	double dz = -view().matrix(2,2);

	geometry.lockRead();
	for (int j=0; j<H; j++) {
		if (stop()) break;
		double v = view().jc2v(j);
		for (int i=0; i<W; i++, ptr++) {
			if (*ptr != painter.background()) continue;
			double u = view().ic2u(i);
			double x = view().uv2x(u,v);
			double y = view().uv2y(u,v);
			double z = view().uv2z(u,v);

			engine()->incBodyCheckId();
			VZone *zone = engine()->where2D(x,y,z, dx,dy,dz);

			// FIXME: maybe check for 3D, lattice or voxel
			if (zone) *ptr = zone->region()->color();
		}
	}
	geometry.unlockRead();
} // drawRegionsLate

/** flood fill scanning for label
 * @param x,y	I/O	input staring pixel to scan, output average pixel
 * @param pos	ref	return list pos points belonging to this label
 */
void D2Layer::scanLabel(Painter& painter, int *x, int *y, Array<IPoint>& pos)
{
#define STACKSIZE 500
	int tos;
	int stackX[STACKSIZE];
	int stackY[STACKSIZE];

	int H = painter.height();	// save variables that might change
	int W = painter.width();	// during the scanning...

	pos.clear();
	dword startColor = painter(*x, *y);
	long long sx=0, sy=0;

	stackX[0] = *x;
	stackY[0] = *y;
	tos = 1;
	while (tos) {
		int xx = stackX[--tos];
		int yy = stackY[  tos];

		dword *pix = painter.pixelPtr(xx,yy);		// get pointer
		if (*pix != startColor) continue;

		// Move to the left side
		while (xx>0 && pix[-1]==startColor) {
			xx--;
			pix--;
		}

		/* scan lines */
		dword *pixelup = (yy>0)? pix-W : NULL;
		dword *pixeldown = (yy<H-1)? pix+W: NULL;

		int up   = 0;
		int down = 0;
		while (*pix == startColor) {
			*pix |= FLAG_LABEL;
			sx += xx;
			sy += yy;
			pos.add(IPoint(xx,yy));

			int pup = (pixelup && *pixelup==startColor);
			if (!up && pup && tos<STACKSIZE) {
				stackX[tos  ] = xx;
				stackY[tos++] = yy-1;
			}
			up = pup;

			int pdown = (pixeldown && *pixeldown==startColor);
			if (!down && pdown && tos<STACKSIZE) {
				stackX[tos  ] = xx;
				stackY[tos++] = yy+1;
			}
			down = pdown;

			if (++xx>=W) break;

			pix++;
			if (pixelup) pixelup++;
			if (pixeldown) pixeldown++;
		}
	}
	assert(pos.size()>0);
	// Find center
	*x = Round((double)sx / double(pos.size()));
	*y = Round((double)sy / double(pos.size()));
} // scanLabel

/** drawLabel if it fits! */
bool D2Layer::drawLabel(Painter& painter, int i, int j, dword color)
{
	int H = painter.height();	// save variables that might change
	int W = painter.width();	// during the scanning...

	if (j<=viewer.font.height() || j>H-viewer.font.height()) return false;

	// Scan text area if a possible label fits (4 characters wide minimum)
	int dx = viewer.font.width() * 2 + 1;
	if (i<=dx || i>W-dx) return false;

	for (int r=Max(0,j-viewer.font.height()); r<Min(j+3,H); r++)
		for (int c=Max(0,i-dx); c<Min(i+dx,W); c++)
			if (painter(c,r) != color)
				return false;

	double u = view().ic2u(i);
	double v = view().jc2v(j);

	double x = view().uv2x(u,v);
	double y = view().uv2y(u,v);
	double z = view().uv2z(u,v);

	engine()->incBodyCheckId();
	VZone *zone = engine()->where2D(x,y,z,
				-view().matrix(0,2),
				-view().matrix(1,2),
				-view().matrix(2,2));

	if (zone) {
		color = FLAG_LABEL | geometry.labelColor;
		painter.pixel(i,j,   FLAG_LABEL | 0xFFFFFF);

		const char *name = zone->region()->label;
		i -= painter.measure(viewer.font, name) / 2;
		painter.printf(viewer.font, i, j-viewer.font.height()-2, color, name);
	}

	return true;
} // drawLabel

/** drawLabels */
void D2Layer::drawLabels(Painter& painter)
{
	srand(314159);

	int H = painter.height();	// save variables that might change
	int W = painter.width();	// during the scanning...
	dword *data = painter.data();	// data pointer
	dword *ptr = data;

	Array<IPoint> pos(1024);

	minLabelSize = 10*viewer.font.width();

	kernel.lock();
	// Check previous labels
	for (int i=(int)labels.size()-1; i>=0; i--) {
		int x = view().u2i(labels[i].x);
		int y = view().v2j(labels[i].y);
		int xx=x, yy=y;
		if (x>=0 && x<W && y>=0 && y<H)
			if ((painter(x,y)&FLAG_INFOMASK) == FLAG_REGION) {
				scanLabel(painter, &x, &y, pos);
				if ((int)pos.size() >= minLabelSize*viewer.font.height() &&
				    drawLabel(painter, xx, yy, painter(xx,yy)))
					continue;
			}
		labels.erase(i);
	}

	// Scan whole image
	for (int j=0; j<H; j++) {
		if (stop()) break;
		for (int i=0; i<W; i++, ptr++) {
			if ((*ptr&FLAG_INFOMASK) != FLAG_REGION) continue;
			int x=i, y=j;

			scanLabel(painter, &x, &y, pos);

			// Ignore small regions
			if ((int)pos.size() < minLabelSize*viewer.font.height()) continue;

			if (drawLabel(painter, x, y, *ptr)) {
				Vector2D v;
				v.x = view().ic2u(x),
				v.y = view().jc2v(y);
				labels.add(v);
				continue;
			}

			for (int k=0; k<20; k++) {
				int m = rand() % pos.size();
				if (drawLabel(painter, pos[m].x, pos[m].y, *ptr)) {
					Vector2D v;
					v.x = view().ic2u(pos[m].x),
					v.y = view().jc2v(pos[m].y);
					labels.add(v);
					break;
				}
			}
		}
	}
	kernel.unlock();
} // drawLabels

/** drawEdge */
void D2Layer::drawEdge(Painter& painter)
{
//#ifdef EXPERIMENTAL
//#if 1
	/* Edge detection if requested */
	// XXX needs speed-up!!!
	int H = painter.height() - 1;	// save variables that might change
	int W = painter.width()  - 1;	// during the scanning...

	for (int j=1; j<H; j++) {
		int jm1 = j-1;
		int jp1 = j+1;
		for (int i=1; i<W; i++) {
			Color32 pixel;
			pixel.val = painter(i,j);
			if (pixel.rgb.alpha != 0x01) continue;
			int im1 = i-1;
			int ip1 = i+1;

			// Construct the following pixels
			//    +----+----+----+
			//    | mm | 0m | pm |
			//    +----+----+----+
			//    | m0 | ij | p0 |
			//    +----+----+----+
			//    | mp | 0p | pp |
			//    +----+----+----+
			Color32 pmm, p0m, ppm;
			Color32 pm0, pp0;
			Color32 pmp, p0p, ppp;

			pmm.val = painter(im1, jm1);
			if (pmm.rgb.alpha == 0) continue;
			p0m.val = painter(i  , jm1);
			if (p0m.rgb.alpha == 0) continue;
			ppm.val = painter(ip1, jm1);
			if (ppm.rgb.alpha == 0) continue;

			pm0.val = painter(im1, j  );
			if (pm0.rgb.alpha == 0) continue;
			pp0.val = painter(ip1, j  );
			if (pp0.rgb.alpha == 0) continue;

			pmp.val = painter(im1, jp1);
			if (pmp.rgb.alpha == 0) continue;
			p0p.val = painter(i  , jp1);
			if (p0p.rgb.alpha == 0) continue;
			ppp.val = painter(ip1, jp1);
			if (ppp.rgb.alpha == 0) continue;

			//if (p0m.val&0x02000000 || pm0.val&0x02000000)
			//	continue;

			// Find the Sobel gradient in the color space
			int gxR = ( pmm.rgb.red +
			          2*p0m.rgb.red +
			            ppm.rgb.red)
			       - (  pmp.rgb.red +
			          2*p0p.rgb.red +
			            ppp.rgb.red);
			int gxG = ( pmm.rgb.green +
			          2*p0m.rgb.green +
			            ppm.rgb.green)
			       - (  pmp.rgb.green +
			          2*p0p.rgb.green +
			            ppp.rgb.green);
			int gxB = ( pmm.rgb.blue +
			          2*p0m.rgb.blue +
			            ppm.rgb.blue)
			       - (  pmp.rgb.blue +
			          2*p0p.rgb.blue +
			            ppp.rgb.blue);

			int gyR= (  ppm.rgb.red +
			          1*pp0.rgb.red +
			            ppp.rgb.red)
			       - (  pmm.rgb.red +
			          2*pm0.rgb.red +
			            pmp.rgb.red);
			int gyG= (  ppm.rgb.green +
			          2*pp0.rgb.green +
			            ppp.rgb.green)
			       - (  pmm.rgb.green +
			          2*pm0.rgb.green +
			            pmp.rgb.green);
			int gyB= (  ppm.rgb.blue +
			          2*pp0.rgb.blue +
			            ppp.rgb.blue)
			       - (  pmm.rgb.blue +
			          2*pm0.rgb.blue +
			            pmp.rgb.blue);
			int g = (int)isqrt(gxR*gxR + gxG*gxG + gxB*gxB
			                 + gyR*gyR + gyG*gyG + gyB*gyB);

#if 0
			int gx = (  Intensity(painter.pixel(im1, jm1)) +
			          2*Intensity(painter.pixel(i,   jm1)) +
			            Intensity(painter.pixel(ip1, jm1)))
			       - (  Intensity(painter.pixel(im1, jp1)) +
			          2*Intensity(painter.pixel(i,   jp1)) +
			            Intensity(painter.pixel(ip1, jp1)));

			int gy = (  Intensity(painter.pixel(ip1, jm1)) +
			          2*Intensity(painter.pixel(ip1, j  )) +
			            Intensity(painter.pixel(ip1, jp1)))
			       - (  Intensity(painter.pixel(im1, jm1)) +
			          2*Intensity(painter.pixel(im1, j  )) +
			            Intensity(painter.pixel(im1, jp1)));
			int g = (int)isqrt(gx*gx + gy*gy);
#endif

			if (g>200) {
				pixel.rgb.alpha = 0x02;
				painter.pixel(i,j,pixel.val);
			}
		}
	}

	for (int j=1; j<H; j++)
		for (int i=1; i<W; i++)
			if (painter.pixel(i,j)&0x02000000)
				painter.pixel(i,j,geometry.zoneColor);
//#endif
//#endif
} // drawEdge

/** draw bodies boundaries */
void D2Layer::draw(Painter& painter)
{
	bool bodiesSelected = 0;

	showWireframe = false;
	showBBox      = false;

	geometry.lockRead();	// No need to lock the kernel, since the segments/nodes are there
				// only the color is not known
	for (int i=0; i<kernel.nGeometryBodies(); i++) {
		VBody *body = kernel.bodies[i];
		if (body->show()==0)
			drawSegments(painter, body, SEGMENT_REGION|SEGMENT_ERROR);
		else {
			// show flag contains info for wireframe and bbox
			if (body->show() & BIT_WIREFRAME)
				showWireframe = true;
			if (body->show() & BIT_BBOX)
				showBBox = true;
			bodiesSelected = true;	// there are bodies selected
		}
	}

	// if exists selected bodies make a second pass to draw them
	if (bodiesSelected) {
		for (int i=0; i<kernel.nGeometryBodies(); i++) {
			VBody *body = kernel.bodies[i];
			if (!body->show()) continue;
			if (body->show() & (BIT_SELECT|BIT_FREEZE)) {
				drawSegments(painter, body, SEGMENT_BODY);
				if (body->show() & BIT_SELECT)
					drawNodes(painter, body);
			} else
				drawSegments(painter, body, SEGMENT_BODY|SEGMENT_REGION|SEGMENT_ERROR);
		}
	}
//	kernel.unlock();
	geometry.unlockRead();
} // draw
