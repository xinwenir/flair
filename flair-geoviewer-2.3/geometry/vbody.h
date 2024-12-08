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

#ifndef __VBODY_H
#define __VBODY_H

#include <iostream>
#ifdef THREAD
#	include <pthread.h>
#endif

#include "array.h"
#include "gbody.h"
#include "conic.h"
#include "vertex2d.h"

class VZone;

#define BODYQUADS	 6
#define BODYCONICS	 6

/* ============================== VBody =============================== */
/** VBody viewer body class */
class VBody {
private:
	GBody	*_body;		// reference to Gbody
	int	_generation;	// Generation number, to be checked for validity against the _body generation

public:
	bool	notref;		// set when NO visible zone is referring to the body
	bool	visible;	// set if body is visible in clipping area
	int	nC;		// Number of conics for the region plotting
	Location location;	// if window is inside or outside the body
				// condition to be checked ONLY if nC=0
	Conic	C[BODYCONICS];	// Conics
	double	acc[BODYQUADS];	// accuracy factor to be multiplied on each quad
	int	c2q[BODYCONICS];	// Conic to quad index conversion
	Array<Vertex2D>	V[BODYCONICS];	// Array of intersecting vertices

protected:
#ifdef THREAD
	pthread_mutex_t	mutex;	// locking mutex when running through threadpool
#endif

public:
	VBody(GBody *abody=NULL);
	~VBody();
	void	init(GBody *abody);

const	char*	name()		const { return _body->name(); }
	BodyType type()		const { return _body->type(); }
const	char*	typeStr()	const { return _body->typeStr(); }
	GBody*	body()		const { return _body; }
	int	id()		const { return _body->id(); }
	dword	color()		const { return _body->color(); }
	dword	lineWidth()	const { return _body->lineWidth(); }
	int	show()		const { return _body->show; }

	// if body is valid
	bool	valid()		const { return _generation == _body->generation(); }
	bool	invalid()	const { return _generation != _body->generation(); }
	void	setInvalid()	      { _generation = -1; }
	void	setValid()	      { _generation = _body->generation(); }

	// Conics
	void	makeConics(const ViewPort& view);
	void	intersectSelf(const ViewPort& view);
	void	intersectViewport(const ViewPort& view);
	void	calculateAccuracy(const ViewPort& view);
	void	updateLocation(const ViewPort& view);
	bool	intersectBody(VBody* abody, const ViewPort& view);
	void	intersectConic(const Conic& conic, const ViewPort& view);

#ifdef THREAD
	// Mutex
	void	lock()		{ pthread_mutex_lock(&mutex); }
	void	unlock()	{ pthread_mutex_unlock(&mutex); }
#else
	void	lock()		{}
	void	unlock()	{}
#endif

	// Vertices
	void	addVertex(const int i, const double x, const double y, VBody* b) {
			lock();
			double t = C[i].getT(x,y);
			V[i].add(Vertex2D(t,x,y,b));
			unlock();
		}
	void	delVertices();
	void	removeInvalidVertices();
	void    markInvalidVertices(const VBody *invalidBody);
	void	removeWrongVertices();

	// Quadrics
	int	nQ()		const { return _body->nQ(); }
const	Quad& Q(int n)		const { return _body->Q(n); }

	bool	inside(const double x, const double y, const double z,
			const double dx, const double dy, const double dz,
			const int ignore1=-1, const int ignore2=-1, const int ignore3=-1)
			{ return _body->inside(x,y,z, dx,dy,dz, ignore1, ignore2, ignore3); }
	// Internal inside2D NON CACHED!
	bool	inside2D(const double x, const double y, const double z,
			 const double dx, const double dy, const double dz,
			 const int ignore1=-1, const int ignore2=-1)
#ifdef EXPERIMENTAL
			{ return _body->inside2D(x,y,z, dx,dy,dz, acc, NULL, ignore1, ignore2); }
#else
			{ return _body->inside2D(x,y,z, dx,dy,dz, acc, ignore1, ignore2); }
#endif

//	bool	isOn(const double x, const double y, const double d) const;
	double	distance(const double x, const double y, const double z) const
			{ return _body->distance(x,y,z); }

	// bounding box
	BBox	bbox2D() const;

	// No visible zone is referring to this body
	void	notReferenced(const bool n) { notref = n; }
	bool	notReferenced() const { return notref || _body->zones.empty(); }
	bool	ignore()	const { return (location != LOCATION_OVERLAP); }

	size_t	memory() const;

private:
	double	accuracy(const int q, const double x, const double y, const double z) const
			{ return Abs(Q(q)(x,y,z)) / Q(q).acc(x,y,z,SMALL1); }
}; // VBody

std::ostream& operator << (std::ostream&, const VBody&);

#endif
