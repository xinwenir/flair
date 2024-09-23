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
		Plane
 * SUBSTITUTE GOODS OR SERVICES, LOSS OF USE, DATA OR PROFITS,
 * OR BUSINESS INTERRUPTION, HOWEVER CAUSED AND ON ANY THEORY
 * OF CONTRACT, WARRANTY, TORT (INCLUDING NEGLIGENCE), PRODUCT
 * LIABILITY OR OTHERWISE, ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGES.
 *
 * Author:	Vasilis.Vlachoudis@cern.ch
 * Date:	11-Oct-2016
 */

#ifndef __TETRAMESH_H
#define __TETRAMESH_H

#include "bbox.h"
#include "list.h"
#include "tetra.h"
#include "matrix4.h"

typedef List<Plane*>	PlaneList;
typedef List<Tetra*>	TetraList;

/* =========================== TetraMesh ============================== */
class TetraMesh {
protected:
	VertexArray	_vertices;	/** array of vertices		*/
	PlaneArray	_planes;	/** array of planes		*/
	TetraArray	_tetrahedra;	/** array of tetrahedra		*/
	BBox		_bbox;		/** bounding box		*/
	double		eps;		/** accuracy of operations	*/

private:
	// temporary arrays used during build up
	Array<PlaneList> vertexPlanes;	/* planes for each vertex	*/
	Array<TetraList> vertexTetra;	/* tetra for each vertex	*/
	int		processed;

public:
	TetraMesh() : eps(FLOATSMALL), processed(0) {}
	~TetraMesh()	{ free(); };

	/** clean up mesh */
	void	free();
	void	freeVertices();
	void	freePlanes();
	void	freeTetrahedra();
//	void	resize(int nv, int ne, int nf);
	void	allocateVertices(int nv);
//	bool	isEmpty()		const	{ return nvertices() == 0; }

	// --- Vertex operations ---
	int	nvertices()		const	{ return _vertices.count(); }
	Vertex*	vertex(int i)			{ return _vertices[i]; }
const	Vertex*	vertex(int i)		const	{ return _vertices[i]; }
	Vertex*	add(const Vertex& v);

	// --- Planes operations ---
	int	nplanes()		const	{ return _planes.count(); }
	Plane*	plane(int i)			{ return _planes[i]; }
const	Plane*	plane(int i)		const	{ return _planes[i]; }
	Plane*	add(int a, int b, int c);

	// --- Tetra operations ---
	int	ntetrahedra()		const	{ return _tetrahedra.count(); }
	Tetra*	tetra(int i)			{ return _tetrahedra[i]; }
const	Tetra*	tetra(int i)		const	{ return _tetrahedra[i]; }
	Tetra*	add(Tetra* t)			{ _tetrahedra << t; return t; }
	Tetra*	add(int a, int b, int c, int d);
//	Tetra*	add(Edge* AB, Edge* BC, Edge* CA){ return add(new Tetra(AB, BC, CA)); }
//	Tetra*	add(Vertex* A, Vertex* B, Vertex* C, bool ab=true, bool bc=true, bool ca=true);

//	void	forEachVertex(VertexFunc, void* arg=NULL);
//	void	forEachEdge(EdgeFunc, void* arg=NULL);
//	void	forEachFace(FaceFunc, void* arg=NULL);
//	void	forEachTetra(TetraFunc, void* arg=NULL);

const	BBox&	bbox()			const	{ return _bbox; }
	void	calcBbox();

//	/* operations */
	double	surface() const;
	double	volume() const;

	void	process();
	void	transform(const Matrix4& matrix);

	size_t	memory() const;

	/** print tetramesh */
friend	std::ostream& operator << (std::ostream&, TetraMesh& tetramesh);
}; // TetraMesh

#endif
