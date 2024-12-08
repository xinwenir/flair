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
 * Date:	1-Jun-2016
 */

#ifndef __MESH_H
#define __MESH_H

#include "bbox.h"
#include "edge.h"
#include "face.h"
#include "list.h"
#include "array.h"
#include "vector.h"
#include "vertex.h"
#include "matrix4.h"
#include "geometry.h"

enum MeshFlagType {
		Mesh_AllFaces,
		Mesh_PositivePid,
		Mesh_NegativePid,
		Mesh_EqualPid,
};

typedef List<GBody*>	BodyList;

/* ============================== Mesh ================================ */
class Mesh {
private:
	VertexArray	_vertices;	/** array of vertices	*/
	EdgeArray	_edges;		/** array of edges	*/
	FaceArray	_faces;		/** array of faces	*/
	BBox		_bbox;		/** bounding box	*/
private:
	double		eps;		/* precision for operations */

public:
	/** constructor*/
	Mesh();

	/** copy constructor */
	Mesh(const Mesh& m);

	/** destructor */
	~Mesh()	{ free(); };

	/** clean up mesh */
	void	free();
	void	freeVertices();
	void	freeEdges();
	void	freeFaces();
	void	resize(int nv, int ne, int nf);
	void	allocateVertices(int nv);
	bool	isEmpty()		const	{ return nvertices() == 0; }

	// --- Vertex operations ---
	int	nvertices()		const	{ return _vertices.count(); }
	Vertex&	vertex(int i)			{ return *_vertices[i]; }
const	Vertex&	vertex(int i)		const	{ return *_vertices[i]; }
	Vertex*	add(const Vertex& v);
	Vertex*	add(Vertex* v);
	int	findVertex(const Vertex&) const;
//	int	vertexIndex(Vertex*)	const	{ return _vertices.find(v); }

	// --- Edge operations ---
	int	nedges()		const	{ return _edges.count(); }
	Edge*	edge(int i)			{ return _edges[i]; }
const	Edge*	edge(int i)		const	{ return _edges[i]; }
	Edge*	add(Vertex* a, Vertex* b, bool show=true);
	Edge*	add(int a, int b, bool show=true)
			{ return add(_vertices[a], _vertices[b], show); }

	/** Find an edge that contains the two vertices
	 * @param a first vertex of the edge
	 * @param b second vertex of the edge
	 * @return pointer to the edge, NULL otherwise */
	int	findEdge(Vertex *a, Vertex *b) const {
			Edge e(a,b);
			return _edges.find(&e);
		}

	// --- Face operations ---
	int	nfaces()		const	{ return _faces.count(); }
	Face*	face(int i)			{ return _faces[i]; }
const	Face*	face(int i)		const	{ return _faces[i]; }
//	Face*	add(Face* f)			{ _faces.add(f); return f; }
	Face*	add(Face* f)			{ _faces << f; return f; }
	Face*	add(Edge* AB, Edge* BC, Edge* CA){ return add(new Face(AB, BC, CA)); }
	Face*	add(Vertex* A, Vertex* B, Vertex* C, bool ab=true, bool bc=true, bool ca=true);
	Face*	add(int a, int b, int c, bool ab=true, bool bc=true, bool ca=true)
				{ return add(_vertices[a], _vertices[b], _vertices[c], ab, bc, ca); }

	void	forEachVertex(VertexFunc, void* arg=NULL);
	void	forEachEdge(EdgeFunc, void* arg=NULL);
	void	forEachFace(FaceFunc, void* arg=NULL);

const	BBox&	bbox()			const	{ return _bbox; }
	void	calcBbox();

	/* --- bodies fitting --- */
	void	fit(BodyList& list, const double ang);

	/* operations */
	bool	isClosed();
	bool	isOrientable();
	void	makeOrientable();
	int	problematicEdges() const;
	double	surface() const;
	double	volume() const;

	void	flip();
	void	process();
	void	transform(const Matrix4& matrix);

	size_t	memory() const;

protected:
	int	fitPlanes(BodyList& list, int minFaces, int& pid);
	int	fitCylinders(BodyList& list, int minFaces, int& pid);

	void	clearProcessedFlag(int b=0, enum MeshFlagType type=Mesh_AllFaces);
	bool	addNeighborsWithSimilarNormal(Face* f, FaceArray& fara, const double dot, int pid);
	bool	addNeighborsCylinder(Face* f, FaceArray& fara, bool& first, double* dot, Vector* axis, int pid, double acc);
	void	faces2Vertices(FaceArray& fara, VertexArray& pts) const;

public:
static	void	createParallelepiped(Mesh& mesh, const Vertex& P,
			const Vector& X, const Vector& Y, const Vector& Z);

	/** print mesh */
friend	std::ostream& operator << (std::ostream&, Mesh& mesh);
}; // Mesh

#endif
