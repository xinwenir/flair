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
 */

#ifndef __FACE_H
#define __FACE_H

#include <assert.h>

#include "edge.h"
#include "array.h"
#include "vertex.h"

/* ============================== Face ================================ */
class Face {
protected:
	// Doesn't own the vertices or the edges
	Vertex*	_vertex[3];		/** vertex indices	*/
	Edge*	_edge[3];		/** pointers to edges	*/
	Vector	_normal;		/** normal of face	*/
	double	_errorNormal;		/** error on normal	*/
	int	_processed;		/** processed flag	*/

public:
	/** constructors */
	Face();
	Face(Edge* ab, Edge* bc, Edge* ca);

	/** copy constructor */
//	Face(const Face&);

	/** destructor */
	~Face() {};

	void	set(Edge* ab, Edge* bc, Edge* ca);
	void	calc();

	/** @return vertex */
	Vertex&	A()	const	{ return *_vertex[0]; }
	Vertex&	B()	const	{ return *_vertex[1]; }
	Vertex&	C()	const	{ return *_vertex[2]; }

	/** @return vertex index */
	Vertex*	vertex(const int i)	const	{ assert(InRange(0,i,2)); return _vertex[i]; }
	Vertex*	operator[] (const int i) const	{ return vertex(i); }

	/** @return edge */
const	Edge*	edge(const int i)	const	{ return _edge[i]; }

	/** @return edge vector */
	Vector	AB()	const	{return Vector(B()-A()); }
	Vector	BC()	const	{return Vector(C()-B()); }
	Vector	CA()	const	{return Vector(A()-C()); }

	Vector	BA()	const	{return Vector(A()-B()); }
	Vector	CB()	const	{return Vector(B()-C()); }
	Vector	AC()	const	{return Vector(C()-A()); }

	/** @return normal */
const	Vector&	normal()		const	{ return _normal; }
	/** set normal */
	void	normal(const Vector& n);

	/** @return vary-center of a face */
	Vertex	center() const	{ return (A() + B() + C()) * (1.0/3.0); }

	/** @return neighbor of edge e */
	Face*	neighbor(int e)	const {
			if (_edge[e]->fA()==this)
				return _edge[e]->fB();
			else
				return _edge[e]->fA();
		}

	int	findEdge(const Vertex* a, const Vertex* b, const bool strict=true) const;

	/** flip a face to up-side-down */
	void	flip();

	double	surface()	const	{ return Abs(AB().cross(AC()).length())/2.0; }

	/** @return point face-plane distance */
	// the plane equation is   n*(r-A) = n*r - n*A = 0
	// distance is given as n*(p-A)/|n| = (n*p-n*A)/|n|
	// normal n is already normalized |n| = 1
	double	distance(const Point& p) const { return normal()*(p - A()); }

#if 0
//	bool	show(const int i)	const	{ return _show[i]; }
	bool	interectRay(const Vertex& O, const Vector& D, double* out, double eps=1e-7);
#endif

	/** set/get processed flag */
	int	processed()		 const	{ return _processed; }
	void	processed(int p)		{ _processed = p; }

protected:
	void	clear();
	double	volume() const;	/* Bizarre function :) */

	/** Print face. */
friend	std::ostream& operator << (std::ostream&, const Face& face);
friend class Mesh;
friend class TetraMesh;
}; // Face

typedef Array<Face*>	FaceArray;
typedef int	(*FaceFunc)(Face*,void*);

#endif
