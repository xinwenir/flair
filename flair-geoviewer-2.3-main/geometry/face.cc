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

#include <string.h>

#include "eps.h"
#include "face.h"

using namespace std;

/** empty constructor */
Face::Face()
{
	clear();
} // Face

/** Face */
Face::Face(Edge* ab, Edge* bc, Edge* ca)
{
	clear();
	set(ab, bc, ca);
} // Face

#if 0
/** copy constructor */
Face::Face(const Face& face)
{
	clear();
	_normal = face._normal;
	for (int i=0; i<3; i++) {
		_vertex[i]   = face._vertex[i];
		_edge[i]     = face._edge[i];
		_neighbor[i] = face._neighbor[i];
//		_show[i]     = face._show[i];
	}
	calc();
} // Face
#endif

//	// clear links with faces
//	for (int i=0; i<3; i++)
//		if (_edge[i]) _edge[i]->remove(this);

/** clear */
void Face::clear()
{
	_processed = 0;
	memset(_vertex,   0, sizeof(_vertex));
	memset(_edge,     0, sizeof(_edge));
//	memset(_neighbor, 0, sizeof(_neighbor));
} // clear

/** set */
void Face::set(Edge* ab, Edge* bc, Edge* ca)
{
	_processed = 0;
//	_stat = Status_Unknown;

	// assign new edges
	_edge[0] = ab;
	_edge[1] = bc;
	_edge[2] = ca;

	for (int i=0; i<3; i++)
		_edge[i]->add(this);

	assert((ab==NULL&&bc==NULL&&ca==NULL) ||
	       (ab!=NULL&&bc!=NULL&&ca!=NULL));

	calc();
} // set

/** calc */
void Face::calc()
{
	/* find the correct vertices */
	for (int i=0; i<3; i++) {
		const Edge* e1 = _edge[i];
		const Edge* e2 = _edge[Next(i,3)];
		// Vertex[i] should be the not-common vertex of e1:e2 from e1
		if (e1->b==e2->a || e1->b==e2->b)
			_vertex[i] = e1->a;
		else {
			assert(e1->a==e2->a || e1->a==e2->b);
			_vertex[i] = e1->b;
		}
	}
	//for (int i=0; i<3; i++)
	//	_vertex[i] = _edge[i]->a;

	_normal = AB() ^ AC();
	double d = _normal.normalize();

	// It is better to expand the AB() x AC() but leads to 
	// a big expression
	// ---> Very crude!!!!

	_errorNormal = Max(A().max(), B().max(), C().max());
	if (d>SMALL) _errorNormal /= d;
} // calc

/** Set the normal vector of the face. */
void Face::normal(const Vector& n)
{
	_normal = n;
	_normal.normalize();
} // normal

/** find edge sharing two vectors
 * @param a,b		vectors of edge
 * @param strict	if true, accept edge only in that order (a,b) else both ways (a,b) and (b,a)
 * @return	edge index or -1 on failure
 */
int Face::findEdge(const Vertex* a, const Vertex* b, bool strict) const
{
	/* strict searching */
	if (_vertex[0]==a && _vertex[1]==b) return 0;
	if (_vertex[1]==a && _vertex[2]==b) return 1;
	if (_vertex[2]==a && _vertex[0]==b) return 2;

	/* reverse searching */
	if (!strict) {
		if (_vertex[0]==b && _vertex[1]==a) return 0;
		if (_vertex[1]==b && _vertex[2]==a) return 1;
		if (_vertex[2]==b && _vertex[0]==a) return 2;
	}

	return -1;
} // findEdge

/** flip */
void Face::flip()
{
	/* swap B and C, keep A as starting vertex */
	Swap(_vertex[1], _vertex[2]);

	/* swap AB with CA */
	Swap(_edge[0], _edge[2]);
//	Swap(_neighbor[0], _neighbor[2]);

	/* flip normal */
	_normal = -_normal;
} // flip

#if 0
/* Möller Trumbore ray-triangle intersection algorithm
 * @param O	ray origin
 * @param D	ray direction
 * @param out	intersection distance
 */
bool Face::interectRay(const Vertex& O, const Vector& D, double* out, double eps )
{
	// Find vectors for two edges sharing V1
	Vector e1 = AB();
	Vector e2 = AC();

	// Begin calculating determinant - also used to calculate u parameter
	Vector P = D ^ e2;

	// if determinant is near zero, ray lies in plane of triangle
	double det = e1 * P;

	// NOT CULLING
	if (det>-eps && det<eps) return false;
	double inv_det = 1.0 / det;

	// calculate distance from V1 to ray origin
	Vector T = O - A();

	//Calculate u parameter and test bound
	double u = (T*P) * inv_det;
	//The intersection lies outside of the triangle
	if (u<0.0 || u>1.0) return false;

	// Prepare to test v parameter
	Vector Q = T ^ e1;

	// Calculate V parameter and test bound
	double v = (D*Q) * inv_det;

	// The intersection lies outside of the triangle
	if (v<0.0 || u+v>1.0) return false;

	double t = (e2*Q) * inv_det;

	if (t > eps) { //ray intersection
		*out = t;
		return true;
	}

	// No hit, no win
	return false;
} // interectRay

#endif

/** break-down of the function needed to calculate
 * the volume of a mesh
 */
double Face::volume() const
{
	return  A().x*(B().y*C().z - B().z*C().y) +
	        B().x*(C().y*A().z - C().z*A().y) +
	        C().x*(A().y*B().z - A().z*B().y);
} // volume

/** print face */
std::ostream& operator << (std::ostream& s, const Face& face)
{
	s << "Face: normal=" << face.normal();
	s << "\tVertex: " << face.A() << " " << face.B() << " " << face.C();
	return s;
} /* operator << */
