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
#include "tetra.h"

using namespace std;

/* ================================== Plane ================================= */
/** set */
void Plane::set(const double aa, const double bb, const double cc, const double dd)
{
	a = aa;
	b = bb;
	c = cc;
	d = dd;
	normalize();
} // set

/** create a plane from 3 vertices */
void Plane::make(const Vertex& A, const Vertex& B, const Vertex& C)
{
	Vector AB = B - A;
	Vector AC = C - A;
	Vector ABxAC = AB ^ AC;
	ABxAC.normalize();
	a = ABxAC.x;
	b = ABxAC.y;
	c = ABxAC.z;
	d = - ABxAC * A;
} // make

/** normalize */
void Plane::normalize()
{
	double n = 1.0 / sqrt(a*a + b*b + c*c);
	a *= n;
	b *= n;
	c *= n;
	d *= n;
} // normalize

/* ================================== Tetra ================================= */
const int Tetra::faceOrder[] = {12, 13, 23, 123};
const int Tetra::faceVertex[4][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}};

/** Tetra */
Tetra::Tetra(Vertex* a, Vertex* b, Vertex* c, Vertex* d)
{
	clear();
	_vertex[0] = a;
	_vertex[1] = b;
	_vertex[2] = c;
	_vertex[3] = d;
} // Tetra

/** clear */
void Tetra::clear()
{
	memset(_vertex,   0, sizeof(_vertex));
	memset(_plane,    0, sizeof(_plane));
	memset(_side,     0, sizeof(_side));
	memset(_neighbor, 0, sizeof(_neighbor));
} // clear

/** set */
void Tetra::set(Plane* ABC, Plane* ABD, Plane* ACD, Plane* BCD)
{
	_plane[0] = ABC;
	_plane[1] = ABD;
	_plane[2] = ACD;
	_plane[3] = BCD;

	Point bary = barycenter();
	for (int i=0; i<4; i++)
		_side[i] = _plane[i]->onPositive(bary);

#if _DEBUG>0
	const double eps = 1e-10;
	assert(ABC->on(A(),eps));
	assert(ABC->on(B(),eps));
	assert(ABC->on(C(),eps));

	assert(ABD->on(A(),eps));
	assert(ABD->on(B(),eps));
	assert(ABD->on(D(),eps));

	assert(ACD->on(A(),eps));
	assert(ACD->on(C(),eps));
	assert(ACD->on(D(),eps));

	assert(BCD->on(B(),eps));
	assert(BCD->on(C(),eps));
	assert(BCD->on(D(),eps));
#endif
} // set

/** @return surface from all faces or only from borders
 * @param border	set to true to get only from borders
 */
double Tetra::surface(bool border) const
{
	double S = 0.0;
	for (int i=0; i<4; i++)
		if (!border || _neighbor[i]==NULL) {
			Vertex& A = *_vertex[faceVertex[i][0]];
			Vertex& B = *_vertex[faceVertex[i][1]];
			Vertex& C = *_vertex[faceVertex[i][2]];
			S += Abs((B-A).cross(C-A).length())/2.0;
		}
	return S;
} // surface

/** @return tetrahedra volume */
double Tetra::volume() const
{
	double v = 0.0;
	double s = -1.0;
	for (int i=0; i<4; i++) {
		Vertex& A = *_vertex[faceVertex[i][0]];
		Vertex& B = *_vertex[faceVertex[i][1]];
		Vertex& C = *_vertex[faceVertex[i][2]];
		v += s*(A.x*(B.y*C.z - B.z*C.y) +
		        B.x*(C.y*A.z - C.z*A.y) +
		        C.x*(A.y*B.z - A.z*B.y));
		s = -s;
	}
	return v;
} // volume

/** check if current tetra is neighbor of t and assign it
 * @param t	neighbor tetra to check
 */
bool Tetra::checkNeighbor(Tetra* t)
{
	int ta, tb, tc, td;
	if ((ta=t->has(_vertex[0]))>=0) {		  // ABC ABD ACD
		if ((tb=t->has(_vertex[1]))>=0) {	  // ABC, ABD
			if (_neighbor[0]==NULL &&
			   (tc=t->has(_vertex[2]))>=0) { // ABC [0]
				_neighbor[0] = t;
				t->_neighbor[t->faceNumber(ta,tb,tc)] = this;
				return true;
			} else
			if (_neighbor[1]==NULL &&
			   (td=t->has(_vertex[3]))>=0) { // ABD [1]
				//assert();
				_neighbor[1] = t;
				t->_neighbor[t->faceNumber(ta,tb,td)] = this;
				return true;
			}
		} else
		if (_neighbor[2]==NULL &&
		   (tc=t->has(_vertex[2]))>=0 &&
		   (td=t->has(_vertex[3]))>=0) { // ACD [2]
			//assert(_neighbor[2]==NULL);
			_neighbor[2] = t;
			t->_neighbor[t->faceNumber(ta,tc,td)] = this;
			return true;
		}
	} else
	if (_neighbor[3]==NULL &&
	   (tb=t->has(_vertex[1]))>=0 &&
	   (tc=t->has(_vertex[2]))>=0 &&
	   (td=t->has(_vertex[3]))>=0) { // BCD [3]
		//assert(_neighbor[3]==NULL);
		_neighbor[3] = t;
		t->_neighbor[t->faceNumber(tb,tc,td)] = this;
		return true;
	}
	return false;
} // checkNeighbor

/** return face number based on the tree index */
int Tetra::faceNumber(int a, int b, int c)
{
	// bubble sort numbers
	if (b>c) Swap(b,c);
	if (a>b) Swap(a,b);
	if (b>c) Swap(b,c);
	int abc = a*100 + b*10 + c;
	for (int i=0; i<4; i++)
		if (faceOrder[i]==abc) return i;
	assert(false);
	return -1;		// it should never reach here
} // faceNumber

/** print tetra */
std::ostream& operator << (std::ostream& s, const Tetra&)
{
	s << "Tetra:";
//	for (int i=0; i<4; i++) {
//		s << "\t" << tetra.face(i);
//		if (i!=3) s << endl;
//	}
	return s;
} /* operator << */
