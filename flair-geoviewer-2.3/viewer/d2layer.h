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

#ifndef __2D_LAYER_H
#define __2D_LAYER_H

#include "geo.h"
#include "array.h"
#include "layer.h"
#ifdef THREAD
#	include "threadpool.h"
#endif

class VBody;
class GRegion;
class GeometryEngine;
class D2Layer;

typedef void (D2Layer::*D2Method)(VBody*, GeometryEngine*);

#ifdef THREAD
/* =============================== BodyWorker =============================== */
class BodyWorker : public ThreadPoolWorker {
public:
	D2Layer*	d2;		/** link to D2layer		*/
	D2Method	method;		/** D2Layer.method to invoke	*/
	GeometryEngine*	engine;		/** engine assigned to worker	*/
	VBody*		body;		/** body to process		*/

public:
	BodyWorker() : d2(NULL), method(NULL), engine(NULL), body(NULL) {}
	virtual void operator()()	{ (*d2.*method)(body,engine); }
}; // BodyWorker

/* =============================== BodyFeeder =============================== */
class BodyFeeder : public ThreadPoolFeeder {
protected:
	int		nworkers;	/** number of workers allocated	*/
	BodyWorker*	workers;	/** array of body workers	*/
	D2Layer*	d2;		/** link to D2layer		*/

	int		idx;		/** current processed body idx	*/
	int		perc_from;	/** percentage range from	*/
	int		perc_span;	/** percentage range span	*/

public:
	BodyFeeder(ThreadPool& p, D2Layer *d) : ThreadPoolFeeder(p),
		  nworkers(0),
		  workers(NULL),
		  d2(d),
		  idx(0) {}

	virtual	~BodyFeeder()		{ if (workers) delete [] workers; }

	// Initialize workers
	void	allocate();
	void	reset(D2Method method, int pfrom, int pto);

virtual	ThreadPoolWorker* feed(int threadId);
}; // BodyFeeder
#endif

/* ================================= D2Layer ================================ */
class D2Layer : public Layer {
public:
	/* drawing flags */
	bool	showBorders;		/** display region boundaries	*/
	bool	showVertex;		/** show vertices		*/
	bool	showLabels;		/** show region labels		*/
	bool	fillRegions;		/** which information to fill	*/

protected:				// Layers
	bool	_projectAll;		/** project all or only invalid	*/
	bool	_projectionValid;	/** if projection is valid	*/

	Array<Vector2D> labels;		/** remember label positions	*/
	int	minLabelSize;		/** minimum label size (width to hold label) */

#ifdef THREAD
	BodyFeeder	feeder;		/** threadpool body feeder	*/
#endif
	bool	showWireframe;
	bool	showBBox;

public:
	D2Layer(const Geometry& g, GeometryKernel& k, GeometryViewer& v);

	/** projection */
	void	projectAll(const bool b)	{ _projectAll = b; }
	void	project();
	bool	isValid()		const	{ return _projectionValid; }
	void	setInvalid()			{ _projectionValid = false; }

	/* information */
	bool	closestVertex(const double u, const double v, const double d, double *vu, double *vv);
	bool	projectVertices(const char axis, Array<double>& vertices);

	/* draw */
	void	draw(Painter& painter);
	void	drawNodes(Painter& painter, VBody *body);
	void	drawVertices(Painter& painter);
	void	drawRegions(Painter& painter);
	void	drawRegionsLate(Painter& painter);
	void	drawLabels(Painter& painter);
	void	drawEdge(Painter& painter);

private:
	void	_intersectBody(VBody* a, GeometryEngine*);
	void	scanSegments();
	void	_scanBodySegments(VBody* body, GeometryEngine* engine);

	void	drawSegments(Painter& painter, VBody *body, dword mask);
	bool	drawSegment(Painter& painter, VBody *body, const int ic,
			const int from, const int to,
			const dword color, const int width,
			double *step);

	void	scanLabel(Painter& painter, int *i, int *j, Array<IPoint>& pos);
	bool	drawLabel(Painter& painter, int i, int j, dword color);

friend	class	GeometryViewer;
friend	class	BodyBodyFeeder;
}; // D2Layer

#endif
