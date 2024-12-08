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
 * Date:        20-Jun-2012
 */

#ifndef __EDGE_H
#define __EDGE_H

#include "array.h"
#include "vertex.h"

class Face;
class Mesh;

enum EdgeLocation {
	EDGE_OUT = -1,
	EDGE_ON  =  0,
	EDGE_A   =  1,
	EDGE_B   =  2
};

/* ============================== Edge ================================ */
/** Edge class, each edge contains two pointers to the vertices
 * but do not own the data
 * It owns the list of pointers to faces that share the edge
 *
 * WARNING for the moment we only treat CLOSED objects
 *         so every edge belongs only to two faces
 */
class Edge {
public:
	Vertex*	a;		/** pointer to 1st vertex		*/
	Vertex*	b;		/** pointer to 2nd vertex		*/
	Face*	fa;		/** first face of edge			*/
	Face*	fb;		/** second face of edge			*/
	bool	 show;		/** edge is visible			*/
//	List	 faces;		/** list of faces sharing that edge	*/
//	void	*reserved;	/** reserved space for operations	*/

public:
	Edge() : a(NULL), b(NULL), fa(NULL), fb(NULL), show(false) {}
	Edge(Vertex* aa, Vertex* bb) { set(aa,bb);}

	/** set vertices */
	void	set(Vertex* aa, Vertex* bb, bool s=true) {
			// force aa to be always the smallest pointer
			if (aa<bb) {
				a = aa;
				b = bb;
			} else {
				b = aa;
				a = bb;
			}
			fa = fb = NULL;
			show = s;
		}

const	Vertex&	A()	const	{ return *a; }
const	Vertex&	B()	const	{ return *b; }

	/** add a face to the edge, only 2 faces are allowed
	 * @return false in case of failure
	 */
	bool	add(Face* f)	{
			if (!fa) {
				fa = f;
				return true;
			}
			assert(fa!=f);
			if (!fb) {
				fb = f;
				return true;
			}
			assert(fb!=f);
			return false;
		}

	/** remove a link with face f */
	void	remove(const Face* f) {
			if (fa==f) fa=NULL;
			if (fb==f) fb=NULL;
		}

	/** return owner faces */
	Face*	fA()	const	{ return fa; }
	Face*	fB()	const	{ return fb; }

	bool	isValid()		const { return a!=NULL && b!=NULL; }
	bool	isDegenerate()		const { return a==b; }

	/** check if edge is similar with the edge defined
	 * by the two vertices
	 * @param aa first vertex
	 * @param bb second vertex */
	bool	same(const Vertex* aa, const Vertex* bb) const {
//			return (a==aa && b==bb);
			if (a==aa && b==bb) return true;
			if (a==bb && b==aa) return true;
			return false;
		}

	/** compare the edges
	 * @param e edge to compare with
	 * @return true if they share the same vertices */
	bool	same(const Edge& e)	const { return same(e.a, e.b); }

	bool operator == (const Edge& e) const { return this->same(e); }
	bool operator != (const Edge& e) const { return !(this->same(e)); }

	/** compare useful for sorted arrays, order edge with pointer values */
static int	compare(Edge* const& a, Edge* const& b) {
			if (a->a > b->a)
				return  1;
			else
			if (a->a < b->a)
				return -1;
			else
			if (a->b > b->b)
				return  1;
			else
			if (a->b < b->b)
				return -1;
			else
			return 0;
		}

	friend class Face;
	friend class Mesh;
}; // class Edge;

typedef Array<Edge*>	EdgeArray;
typedef int	(*EdgeFunc)(Edge*, void*);

#endif
