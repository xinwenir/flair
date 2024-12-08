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
 * Date:	11-Oct-2016
 */

//#define DUMP


#include <string>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <iostream>

#include "eps.h"
#include "tetramesh.h"

using namespace std;

/** free */
void TetraMesh::free()
{
	freeTetrahedra();
	freePlanes();
	freeVertices();
} // free

/** freeVertices */
void TetraMesh::freeVertices()
{
	for (int i=0; i<nvertices(); i++) delete _vertices[i];
	_vertices.clear();
	_bbox.reset();
} // freeVertices

/** freePlanes */
void TetraMesh::freePlanes()
{
	for (int i=0; i<nplanes(); i++) delete _planes[i];
	_planes.clear();
} // freePlanes

/** freeTetrahedra */
void TetraMesh::freeTetrahedra()
{
	for (int i=0; i<ntetrahedra(); i++) delete _tetrahedra[i];
	_tetrahedra.clear();
} // freeTetrahedra

#if 0
/** allocate initial arrays for vertices, edges, and faces */
void TetraMesh::resize(int nv, int ne, int nf)
{
	if (nv>0) _vertices.resize(nv);
	if (ne>0) _edges.resize(ne);
	if (nf>0) _faces.resize(nf);
} // resize
#endif

/** allocate initial number of vertices
 * DANGEROUS
 */
void TetraMesh::allocateVertices(int nv)
{
	if (nvertices() == nv) return;
	if (nv < nvertices()) {
		assert(0);	// FIXME not tested yet
		// delete excess
		for (int i=nv; i<nvertices(); i++)
			delete _vertices[i];
		_vertices.resize(nv);
	} else {
		int old = nvertices();
		_vertices.resize(nv);
		for (int i=old; i<nv; i++)
			_vertices << new Vertex();
	}
} // allocateVertices

/** create a new vertex and add it
 * @param v to add
 * @return the index of vertex
 */
Vertex* TetraMesh::add(const Vertex& v)
{
	_bbox.add(v);
	Vertex *vv = new Vertex(v);
	_vertices << vv;
	return vv;
} // add

/** add a new tetra */
Tetra* TetraMesh::add(int a, int b, int c, int d)
{
	Vertex* A = vertex(a);
	Vertex* B = vertex(b);
	Vertex* C = vertex(c);
	Vertex* D = vertex(d);

	// Add planes
	Plane* ABC = add(a,b,c);
	Plane* ABD = add(a,b,d);
	Plane* ACD = add(a,c,d);
	Plane* BCD = add(b,c,d);

	if (!vertexTetra.count()) vertexTetra.allocate(nvertices());

	Tetra* t = new Tetra(A,B,C,D);
	t->set(ABC, ABD, ACD, BCD);

	// find neighbors
	int n = 0;	// neighbor count
	ListIterator<Tetra*>	iter(vertexTetra[a]);
	while (iter && n<4) {
		Tetra* vt = iter++;
		if (vt->checkNeighbor(t)) n++;
	}
	if (n<4) {
		iter.reset(vertexTetra[b]);
		while (iter && n<4) {
			Tetra* vt = iter++;
			if (vt->checkNeighbor(t)) n++;
		}
		if (n<4) {
			iter.reset(vertexTetra[c]);
			while (iter && n<4) {
				Tetra* vt = iter++;
				if (vt->checkNeighbor(t)) n++;
			}
			if (n<4) {
				iter.reset(vertexTetra[d]);
				while (iter && n<4) {
					Tetra* vt = iter++;
					if (vt->checkNeighbor(t)) n++;
				}
			}
		}
	}
	if (n<4)
		for (int i=0; i<ntetrahedra(); i++)
			tetra(i)->checkNeighbor(t);

	// create caching lists
	vertexTetra[a] << t;
	vertexTetra[b] << t;
	vertexTetra[c] << t;
	vertexTetra[d] << t;

	// add the tetra
	_tetrahedra << t;
	return t;
} // add

/**  add a new plane or return an existing one that passes from 
 * the three vertices
 * @param A,B,C	vertex from where plane lies
 * @return plane pointer inside _planes array
 */
Plane* TetraMesh::add(int a, int b, int c)
{
	Vertex& A = *vertex(a);
	Vertex& B = *vertex(b);
	Vertex& C = *vertex(c);

	if (!vertexPlanes.count()) vertexPlanes.allocate(nvertices());

	// FIXME we should keep a list of vertices on each plane to make a better fitting!
	//find if plane exists
	const double e=1e-10;

	processed++;
	// First check on existing planes known from any of the vertices
	ListIterator<Plane*>	iter(vertexPlanes[a]);
	while (iter) {
		Plane* p = iter++;
		if (p->processed() == processed) continue;
		p->processed(processed);
		if (p->on(A, e) &&
		    p->on(B, e) &&
		    p->on(C, e)) {
			if (vertexPlanes[b].find(p) == -1) vertexPlanes[b] << p;
			if (vertexPlanes[c].find(p) == -1) vertexPlanes[c] << p;
			return p;
		}
	}
	iter.reset(vertexPlanes[b]);
	while (iter) {
		Plane* p = iter++;
		if (p->processed() == processed) continue;
		p->processed(processed);
		if (p->on(A, e) &&
		    p->on(B, e) &&
		    p->on(C, e)) {
			if (vertexPlanes[a].find(p) == -1) vertexPlanes[a] << p;
			if (vertexPlanes[c].find(p) == -1) vertexPlanes[c] << p;
			return p;
		}
	}
	iter.reset(vertexPlanes[c]);
	while (iter) {
		Plane* p = iter++;
		if (p->processed() == processed) continue;
		p->processed(processed);
		if (p->on(A, e) &&
		    p->on(B, e) &&
		    p->on(C, e)) {
			if (vertexPlanes[a].find(p) == -1) vertexPlanes[a] << p;
			if (vertexPlanes[b].find(p) == -1) vertexPlanes[b] << p;
			return p;
		}
	}

	for (int i=0; i<nplanes(); i++) {
		if (plane(i)->processed() == processed) continue;
		if (plane(i)->on(A, e) &&
		    plane(i)->on(B, e) &&
		    plane(i)->on(C, e))
			return plane(i);
	}

	Plane* p = new Plane(A,B,C);

	// add cached lists
	vertexPlanes[a] << p;
	vertexPlanes[b] << p;
	vertexPlanes[c] << p;

	_planes << p;
	return p;
} // findPlane

#if 0
/** for each vertex call function func with argument arg
* @param func function to call
* @param arg argument to pass */
void TetraMesh::forEachVertex(VertexFunc func, void* arg)
{
	ArrayIterator<Vertex*> iter(_vertices);
	while (iter)
		if ((func)(iter++, arg)) break;
} // forEachVertex

/** for each edge call function func with argument arg
* @param func function to call
* @param arg argument to pass */
void TetraMesh::forEachEdge(EdgeFunc func, void* arg)
{
	ArrayIterator<Edge*> iter(_edges);
	while (iter)
		if ((func)(iter++, arg)) break;
} // forEachEdge

/** for each face call function func with argument arg
* @param func function to call
* @param arg argument to pass */
void TetraMesh::forEachFace(FaceFunc func, void* arg)
{
	ArrayIterator<Face*> iter(_faces);
	while (iter)
		if ((func)(iter++, arg)) break;
} // forEachFace
#endif

/** calcBbox */
void TetraMesh::calcBbox()
{
	_bbox.reset();
	for (int i=0; i<nvertices(); i++)
		_bbox.add(*vertex(i));
} // calcBbox

/** process mesh */
void TetraMesh::process()
{
} // process

/** transform */
void TetraMesh::transform(const Matrix4& matrix)
{
	for (int i=0; i<nvertices(); i++)
		*vertex(i) = matrix * *vertex(i);
} // transform

/* @return surface of the mesh */
double TetraMesh::surface() const
{
	double s = 0.0;
	for (int i=0; i<ntetrahedra(); i++)
		s += tetra(i)->surface(true);
	return s;
} // surface

/* @return volume of the mesh */
double TetraMesh::volume() const
{
	double vol = 0.0;
	for (int i=0; i<ntetrahedra(); i++)
		vol += tetra(i)->volume();
	return vol/6.0;
} // volume

/** @return memory used by TetraMesh */
size_t TetraMesh::memory() const
{
	return  sizeof(*this)
	      + _vertices.memory()
	      + _tetrahedra.memory()
	      + sizeof(Vertex) * nvertices()
	      + sizeof(Tetra) * ntetrahedra();
} // memory

/** Print mesh. */
ostream& operator << (ostream& s, TetraMesh& tetramesh)
{
	s << "TetraMesh" << endl;
	s << " vertices: " << tetramesh.nvertices() << endl;
//	for (int i=0; i<tetramesh.nvertices(); i++)
//		s << "\t" << i << ": " << tetramesh.vertex(i) << endl;
	s << "   planes: " << tetramesh.nplanes() << endl;
	s << "    tetra: " << tetramesh.ntetrahedra() << endl;
//	for (int i=0; i<tetramesh.ntetrahedra(); i++) {
//		Tetra* t = tetramesh.tetra(i);
//		cout << endl;
//	}
	int neighbors = 0;
	int borders   = 0;
	for (int i=0; i<tetramesh.ntetrahedra(); i++) {
		int n = tetramesh.tetra(i)->nneighbors();
		neighbors += n;
		borders   += 4 - n;
	}
	neighbors /= 2;
	s << "neighbors: " << neighbors << endl;
	s << "  borders: " << borders << endl;
	s << "  surface: " << tetramesh.surface() << endl;
	s << "   volume: " << tetramesh.volume() << endl;
	s << "     bbox: " << tetramesh.bbox() << endl;
	return s;
} /* operator << */
