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

#ifndef __TETRA_H
#define __TETRA_H

#include <assert.h>
#include "face.h"
#include "array.h"

/* ============================= Plane ================================ */
class Plane {
protected:
	double	a, b, c, d;		/** plane coefficients	*/

private:
	int	_processed;		/** processed flag			*/

public:
	Plane() : a(0.), b(0.), c(0.), d(0.) {}
	Plane(const double aa, const double bb, const double cc, const double dd) :
		a(aa), b(bb), c(cc), d(dd) { normalize(); }
	Plane(const Vertex& A, const Vertex& B, const Vertex& C)	{ make(A, B, C); }

	void	set(const double aa, const double bb, const double cc, const double dd);
	void	make(const Vertex& A, const Vertex& B, const Vertex& C);
	void	normalize();

	/** @return distance of a point from plane	*/
	double	distance(const Point& p)	const
			{ return a*p.x + b*p.y + c*p.z + d; }

	/** @return true if point p lies on plane */
	bool	on(const Point& p, const double eps)	const
			{ return Eq0(distance(p), eps); }

	/** @return true if point p lies on positive side */
	bool	onPositive(const Point& p) const { return distance(p) >= 0.0; }

	/** set/get processed flag */
	int	processed()		 const	{ return _processed; }
	void	processed(int p)		{ _processed = p; }
}; // class Plane

/* ============================= Tetra ================================ */
/**
 * Class implementing the basic tetrahedro
 *
 *                      v
 *                    .
 *                  ,/
 *                 /
 *              2
 *            ,/|`\
 *          ,/  |  `\
 *        ,/    '.   `\
 *      ,/       |     `\
 *    ,/         |       `\
 *   0-----------'.--------1 --> u
 *    `\.         |      ,/
 *       `\.      |    ,/
 *          `\.   '. ,/
 *             `\. |/
 *                `3
 *                   `\.
 *                      ` w
 */
class Tetra {
protected:
	// FIXME what do we really need?
	Vertex*	_vertex[4];		/** pointers to vertices		*/
	Plane*	_plane[4];		/** pointers to planes			*/
	bool	_side[4];		/** positive/negative side of plane	*/
	Tetra*	_neighbor[4];		/** neighbor tetrahedra			*/
	int	_processed;		/** processed flag			*/

private:
const	static int faceOrder[4];
const	static int faceVertex[4][3];

public:
	Tetra()		{ clear(); }
	Tetra(Vertex* a, Vertex* b, Vertex* c, Vertex* d);
//	Tetra(const Tetra&);
	~Tetra() {};

	void	set(Plane* ABC, Plane* ABD, Plane* ACD, Plane* BCD);
	Point	barycenter()	const { return 0.25*(*_vertex[0] + *_vertex[1] + *_vertex[2] + *_vertex[3]); }

const	Vertex&	A()		const	{ return *_vertex[0]; }
const	Vertex&	B()		const	{ return *_vertex[1]; }
const	Vertex&	C()		const	{ return *_vertex[2]; }
const	Vertex&	D()		const	{ return *_vertex[3]; }

	bool	checkNeighbor(Tetra* t);

	/** @return surface only of borders */
	double	surface(bool border)	const;
	double	volume() const;

#if 0
	/** @return normal */
//const	Vector&	normal()		const	{ return _normal; }

	/** @return point face-plane distance */
	// the plane equation is   n*(r-A) = n*r - n*A = 0
	// distance is given as n*(p-A)/|n| = (n*p-n*A)/|n|
	// normal n is already normalized |n| = 1
	double	distance(const Point& p) const { return normal()*(p - A()); }

//	bool	show(const int i)	const	{ return _show[i]; }
	bool	interectRay(const Vertex& O, const Vector& D, double* out, double eps=1e-7);
#endif

	int	nneighbors()		const {
			int n=0;
			for (int i=0; i<4; i++)
				if (_neighbor[i]) n++;
			return n;
		}

	/** set/get processed flag */
	int	processed()		 const	{ return _processed; }
	void	processed(int p)		{ _processed = p; }

protected:
	void	clear();
	int	has(const Vertex* v) {
			for (int i=0; i<4; i++)
				if (v==_vertex[i]) return i;
			return -1;
		}
	int	faceNumber(int a, int b, int c);

	/** Print face. */
friend	std::ostream& operator << (std::ostream&, const Tetra& tetra);
//friend class TetraMesh;
}; // Tetra

typedef Array<Plane*>	PlaneArray;
typedef Array<Tetra*>	TetraArray;
//typedef int	(*TetraFunc) (Tetra*, void*);

#endif
