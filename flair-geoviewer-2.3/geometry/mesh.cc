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

//#define DUMP

#include <string>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <iostream>

#include "eps.h"
#include "mesh.h"
#include "gbody.h"
#include "matrix.h"
#include "matrix2.h"
#include "matrix3.h"

using namespace std;

static bool	fitPlane(VertexArray& verts, double* a, double* b, double* c, double* d, double* err, const double eps);
static BodyType plane2Fluka(double a, double b, double c, double d, double* what);
static bool	fitCircle(VertexArray& pts, int ia, int ib, double* a, double* b, double* r, double* err);
//static bool	fitSphere(VertexArray& pts, double *x, double *y, double *z, double *r, double *err, const double eps);
//static bool	fitQuadratic(VertexArray& pts, double C[10], double *err, const double eps);
//static bool	fitEllipticalCylinder(VertexArray& pts, double *y, double *z, double *Ry, double *Rz, double *err);
#if _DEBUG>1
static bool	writeArray2Stl(const char*, FaceArray&);
#endif

/** constructor*/
Mesh::Mesh() : eps(FLOATSMALL)
{
	_vertices.compare(Vertex::compare);
	_edges.compare(Edge::compare);
} // Mesh

/** copy constructor */
Mesh::Mesh(const Mesh& m) :
	_vertices(m._vertices),
	_edges(m._edges),
	_faces(m._faces),
	eps(m.eps)
{
	_vertices.compare(Vertex::compare);
	_edges.compare(Edge::compare);
} // Mesh

/** free */
void Mesh::free()
{
	freeFaces();
	freeEdges();
	freeVertices();
} // free

/** freeVertices */
void Mesh::freeVertices()
{
	for (int i=0; i<nvertices(); i++) delete _vertices[i];
	_vertices.clear();
	_bbox.reset();
} // freeVertices

/** freeEdges */
void Mesh::freeEdges()
{
	for (int i=0; i<nedges(); i++) delete _edges[i];
	_edges.clear();
} // freeEdges

/** freeFaces */
void Mesh::freeFaces()
{
	for (int i=0; i<nfaces(); i++) delete _faces[i];
	_faces.clear();
} // freeFaces

/** allocate initial arrays for vertices, edges, and faces */
void Mesh::resize(int nv, int ne, int nf)
{
	if (nv>0) _vertices.resize(nv);
	if (ne>0) _edges.resize(ne);
	if (nf>0) _faces.resize(nf);
} // resize

/** allocate initial number of vertices
 * DANGEROUS
 */
void Mesh::allocateVertices(int nv)
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
Vertex* Mesh::add(const Vertex& v)
{
	int idx = findVertex(v);
	if (idx>=0) return _vertices[idx];

	_bbox.add(v);
	Vertex *vv = new Vertex(v);
	_vertices << vv;
	return vv;
} // add

/** add vertex pointer
 * @param v to add
 * @return the index of vertex
 */
Vertex* Mesh::add(Vertex* v)
{
	assert(false);
	// FIXME can create memory leaks!
	// Who is owning the v????
#if 0
	int idx = vertexIndex(v);
	if (idx>=0) return v;

	_bbox.add(*v);
	_vertices << v;
#endif
	return v;
} // add

/**
 * @param v to search for
 * @return the index of closest vertex, -1 on failure
 */
int Mesh::findVertex(const Vertex& v) const
{
	Vertex lower(v.x-eps, v.y-eps, v.z-eps);

	// find lower location in z
	int i= _vertices.search(&lower);
	if (i>=nvertices()) return -1;
	if (i<0) i++;

	// scan all points until vertex.z+eps for the first closest match
	Vertex upper(v.x+eps, v.y+eps, v.z+eps);
	for (; i<nvertices(); i++) {
		if (vertex(i) > upper)
			break;
		if (vertex(i).isClose(v, eps))
			return i;
	}
	return -1;
} // findVertex

/** add a new edge */
Edge* Mesh::add(Vertex* a, Vertex* b, bool show)
{
	assert(a!=NULL && b!=NULL);
	int i = findEdge(a, b);
	if (i>=0) {
		assert(edge(i)->show == show);
		return edge(i);
	}
	Edge* e = new Edge(a, b);
	_edges << e;
	e->show = show;
	return e;
} // add

/** add face by vertex pointers */
Face* Mesh::add(Vertex* A, Vertex* B, Vertex* C, bool ab, bool bc, bool ca)
{
	Edge *AB = add(A, B, ab);
	Edge *BC = add(B, C, bc);
	Edge *CA = add(C, A, ca);

	return add(new Face(AB, BC, CA));
} // add

/** for each vertex call function func with argument arg
* @param func function to call
* @param arg argument to pass */
void Mesh::forEachVertex(VertexFunc func, void* arg)
{
	ArrayIterator<Vertex*> iter(_vertices);
	while (iter)
		if ((func)(iter++, arg)) break;
} // forEachVertex

/** for each edge call function func with argument arg
* @param func function to call
* @param arg argument to pass */
void Mesh::forEachEdge(EdgeFunc func, void* arg)
{
	ArrayIterator<Edge*> iter(_edges);
	while (iter)
		if ((func)(iter++, arg)) break;
} // forEachEdge

/** for each face call function func with argument arg
* @param func function to call
* @param arg argument to pass */
void Mesh::forEachFace(FaceFunc func, void* arg)
{
	ArrayIterator<Face*> iter(_faces);
	while (iter)
		if ((func)(iter++, arg)) break;
} // forEachFace

/** calcBbox */
void Mesh::calcBbox()
{
	_bbox.reset();
	for (int i=0; i<nvertices(); i++)
		_bbox.add(vertex(i));
} // calcBbox

/** process mesh */
void Mesh::process()
{
} // process

/** flip faces of mesh */
void Mesh::flip()
{
	for (int i=0; i<nfaces(); i++)
		_faces[i]->flip();
} // flip

/** transform */
void Mesh::transform(const Matrix4& matrix)
{
	for (int i=0; i<nvertices(); i++)
		vertex(i) = matrix * vertex(i);
} // transform

/** @return true if every edge is shared by two faces */
bool Mesh::isClosed()
{
	for (int i=0; i<nedges(); i++) {
		const Edge* e = _edges[i];
		if (!e->fA() || !e->fB()) return false;
	}
	return true;
} // isClosed

/** @return number of problematic edges (with only one face) as owner */
int Mesh::problematicEdges() const
{
	int count = 0 ;
	for (int i=0; i<nedges(); i++) {
		const Edge* e = _edges[i];
		if (!e->fA() || !e->fB()) count++;
	}
	return count;
} // problematicEdges

/** @return true if it is a closed mesh and every edge is transverse
 *               one way from one face and the opposite from the other face
 */
bool Mesh::isOrientable()
{
	for (int i=0; i<_faces.size(); i++) {
		Face* f = face(i);
		for (int j=0; j<3; j++) {
			const Face* n = f->neighbor(j);
			if (n == NULL) return false;
			const Vertex* a = f->vertex(j);
			const Vertex* b = f->vertex(Next(j,3));
			if (n->findEdge(b,a,true)<0) return false;
		}
	}
	return true;
} // isOrientable

/** makeOrientable */
void Mesh::makeOrientable()
{
	for (int i=0; i<_faces.size(); i++) {
		Face* f = face(i);
		for (int j=0; j<3; j++) {
			Face* n = f->neighbor(j);
			assert(n != NULL);
			const Vertex* a = f->vertex(j);
			const Vertex* b = f->vertex(Next(j,3));
			if (n->findEdge(b,a,true)<0) n->flip();
		}
	}
} // makeOrientable

/* @return surface of the mesh */
double Mesh::surface() const
{
	double s = 0.0;
	for (int i=0; i<nfaces(); i++)
		s += face(i)->surface();
	return s;
} // surface

/* @return volume of the mesh */
double Mesh::volume() const
{
	double vol = 0.0;
	for (int i=0; i<nfaces(); i++)
		vol += face(i)->volume();
	return vol/6.0;
} // volume

/** clear processed flag from all faces */
void Mesh::clearProcessedFlag(int b, enum MeshFlagType type)
{
	for (int i=0; i<nfaces(); i++)
		switch (type) {
			case Mesh_PositivePid:
				if (face(i)->processed()>0)
					face(i)->processed(b);
				break;
			case Mesh_NegativePid:
				if (face(i)->processed()<0)
					face(i)->processed(b);
				break;
			case Mesh_EqualPid:
				if (face(i)->processed() == b)
					face(i)->processed(0);
				break;
			default:
				face(i)->processed(b);
		}
} // clearProcessedFlag

/** @return a list of all body cards necessary for the mesh */
void Mesh::fit(BodyList& list, const double ang)
{
	FaceArray	fara;	// list of connected faces with similar normal
	VertexArray	pts;	// vertices of fara
//	pts.compare(Cmp);	// sort everything with pointers

	double dot = cos(ang);

	// strict check for planes
	clearProcessedFlag(0);
	int pid = 0;	// processed id

	// Scan for planes
	fitPlanes(list, 4, pid);

	// Scan for cylinders
	clearProcessedFlag(0, Mesh_PositivePid);
	int cnt=0;
	for (int i=0; i<nfaces(); i++) {
		Face* f = face(i);
		if (f->processed()) continue;	// ignore processed face

		pid++;	// increase processed id
		fara.clear();	// clear face array
		fara << f;	// add first face
		double d = dot;
		Vector axis;
		bool first = true;
		const double accuracy = 0.001;
		addNeighborsCylinder(f, fara, first, &d, &axis, pid, accuracy);
		if (fara.size() >= 5) {
			faces2Vertices(fara, pts);

			double a, b, r, err=accuracy;
			int found = -1;
			for (int j=0; j<3; j++) {
				int jn = (j+1)%3;
				double aj, bj, rj, errj;
				if (fitCircle(pts, j, jn, &aj, &bj, &rj, &errj)) {
					if (errj < err) {
						a = aj;
						b = bj;
						r = rj;
						err = errj;
						found = j;
					}
				}
			}
			if (found>=0) {
				cnt++;
#if _DEBUG>1
				printf("\nCylinder %d\n",cnt);
				printf("Best Cylinder %d-%d\n",found,(found+1)%3);
				cout << i << " Faces:" << fara.count() << endl;
				cout << "axis=" << axis << endl;
				printf("a=%g b=%g r=%g err=%g\n",a,b,r,err);
#endif
				int ia = found;
				int ib = (ia+1)%3;

				// flag them as negative
				for (int j=0; j<fara.size(); j++)
					fara[j]->processed(-pid);

				// find all potentially not processed faces
				// that might belong to the same plane
				for (int j=0; j<nfaces(); j++) {
					Face* fj = face(j);
					if (fj->processed()<0) continue;	// ignore processed face

					// All 3 points should be part of the circle
					bool ok = true;
					for (int k=0; k<3; k++) {
						double x = (*fj->vertex(k))[ia];
						double y = (*fj->vertex(k))[ib];
						double e = Sqr(x-a) + Sqr(y-b) - r*r;
						if (e>eps*(r+Abs(a)+Abs(b))) {
							ok = false;
							break;
						}
					}
					if (!ok) continue;

					fj->processed(-pid);
					fara << fj;
				}
				cout << i << " Added Faces:" << fara.count() << endl;

#if _DEBUG>1
				char filename[100];
				sprintf(filename,"cylinder_%03d.stl",cnt);
				writeArray2Stl(filename, fara);
#endif
			} else {
				// Could become smarter by recursively looking the neighbors of the faces
				// the opposite of addNeighborsCylinder
				clearProcessedFlag(pid, Mesh_EqualPid);
			}
		} else
			clearProcessedFlag(pid, Mesh_EqualPid);
	}
#if 0
		// Sphere
		double x, y, z, r, errsph=INFINITE;
		if (fitSphere(pts, &x, &y, &z, &r, &errsph, eps)) {
			cout << "# SPH " << x
			       << "  " << y
			       << "  " << z
			       << "  " << r
			       << " err=" << errsph << endl;
		}

		// Quadratic
		double C[10], errqua=INFINITE;
		if (fitQuadratic(pts, C, &errqua, eps)) {
//			if (errqua < eps) {
				cout << "QUA ";
				for (int k=0; k<10; k++)
					cout << " " << C[k];
				cout << endl;
			cout << "Err=" << errqua << endl;
//			}
			flagProcessedFaces(fara, 2);		// mark as 2
		}

		if (errsph < eps) {
			cout << "SPH " << x
			       << "  " << y
			       << "  " << z
			       << "  " << r
			       << " err=" << errsph << endl;

			double what[4];
			what[0] = x;
			what[1] = y;
			what[2] = z;
			what[3] = r;
			GBody* body = GBody::newBody("sphere",SPHbody);
			body->setWhat(what, NULL);
//#ifdef _DUMP
			cout << GBody::typeStr(SPHbody);
			for (int kk=0; kk<4; kk++)
				cout << "  " << setprecision(16) << what[kk];
			cout << endl;
//#endif
			list.append(body);
			flagProcessedFaces(fara, 2);		// mark as 2
		}
	}

	// check for planes (strict check)
//	fitPlanes(list, 2);

	// 3rd pass: create planes for remaining faces

	// 4th pass: check if possible to merge bodies to a macro body
	string type="";
	Array< Array<Face> > areas = getNeighbourAreasWithSimilarNormal(_faces);
	//check what figure it is :plan/cone/sphere/cylinder/etc..
	for (int i=0; i<areas.count(); i++)
	{
		//for plan
		Face firstFace = areas.get(i).get(0);
		Vertex vertex = firstFace.getVertexIndices()[0];
		BodyDefinition bd_pla = getPlan(areas.get(i), vertex);
		type = bd_pla.getType();
		if (type.compare("") != 0){
			result.push_back(bd_pla);
		}
	}

	//for sphere
	BodyDefinition bd_sph = getSphere(_vertices);
	type = bd_sph.getType();
	if (type.compare("") != 0) {
		result.push_back(bd_sph);
	}

	/*
	//for cylinder
	BodyDefinition bd = getCylinder(_vertices);
	type = bd.getType();
	if (type.compare("") != 0)
		//result.push_back(bd);
	*/

	//for eliptical cylinder
	BodyDefinition bd_cyl = getEllipticalCylinder(_vertices);
	type = bd_cyl.getType();
	if (type.compare("") != 0)
	{
		result.push_back(bd_cyl);
	}

	//print result
	for (vector<BodyDefinition>::iterator it = result.begin(); it!=result.end(); ++it)
	{
		cout<<*it<<endl;
	}
	return result;
#endif
} // fit

/** fitPlanes, using a strict check
 * @param list		list to add bodies to
 * @param pid		processed id of planes
 * @param minFaces	minimum number of faces to be considered as plane
 */
int Mesh::fitPlanes(BodyList& list, int minFaces, int& pid)
{
	FaceArray	fara;	// list of connected faces with similar normal
	VertexArray	pts;	// vertices of fara
//	pts.compare(Cmp);	// sort everything with pointers

	int count = 0;		// count number of planes found
	const double dotplane = 1.0 - eps;		// close to 1.0 with the mesh accuracy
#if _DEBUG>1
	int cnt = 0;
#endif
	// First scan for planes with many faces
	for (int i=0; i<nfaces(); i++) {
		Face* f = face(i);
		if (f->processed()) continue;	// ignore processed face

		pid++;	// increase processed id
		fara.clear();	// clear face array
		fara << f;	// add first face
		addNeighborsWithSimilarNormal(f, fara, dotplane, pid);
		if (fara.size()<minFaces) continue;

		faces2Vertices(fara, pts);

		// Plane fitting
		double a, b, c, d, err;
		if (!fitPlane(pts, &a, &b, &c, &d, &err, eps)) continue;
#if _DEBUG>1
		cout << "PLA " << a
		       << "  " << b
		       << "  " << c
		       << "  " << d
		       << " err=" << err << endl;
#endif
		// Convert to FLUKA plane
		double what[6];
		BodyType bt = plane2Fluka(a,b,c,d,what);
#if _DEBUG>1
		const char* plas[] = {"YZP", "XZP", "XYP", "PLA" };
		cout << plas[bt];
		for (int ii=0; ii<(bt==PLAbody?6:1); ii++)
			cout << " " << what[ii];
		cout << endl;
#endif
		// flag them as negative
		for (int j=0; j<fara.size(); j++)
			fara[j]->processed(-pid);

		// find all potentially not processed faces
		// that might belong to the same plane
		for (int j=0; j<nfaces(); j++) {
			Face* fj = face(j);
			if (fj->processed()<0) continue;	// ignore processed face

			// Check normal first
			if (Abs(f->normal() * fj->normal()) < dotplane) continue;

			// check if belongs to plane
			Vertex& A = f->A();
			Vertex& B = fj->A();
			if (Abs(f->normal() * (A-B)) > eps * Max(1.0, A.max(), B.max())) continue;

			fj->processed(-pid);
			fara << fj;
		}

		char name[9];
		nameNumber(name, sizeof(name), "plane", list.count()+1, 2);
		GBody* body = GBody::newBody(name,bt);
		body->setWhat(what, NULL);
		body->id(pid);
		list << body;
		count++;
#if _DEBUG>1
		char filename[100];
		sprintf(filename,"plane_%03d.stl",++cnt);
		writeArray2Stl(filename, fara);
#endif
	}
	return count;
} // fitPlanes

/** addNeighborsWithSimilarNormal */
bool Mesh::addNeighborsWithSimilarNormal(Face* f, FaceArray& fara, const double dot, int pid)
{
	assert(f->processed()==0);
	f->processed(pid);
	bool added = false;
	for (int i=0; i<3; i++) {
		Face* n = f->neighbor(i);
		if (n->processed()==0) {
			double d = f->normal() * n->normal();
			if (d<dot) continue;
			fara << n;
			addNeighborsWithSimilarNormal(n, fara, dot, pid);
			added = true;
		}
	}
	return added;
} // addNeighborsWithSimilarNormal

/** addNeighborsCylinder */
bool Mesh::addNeighborsCylinder(Face* f, FaceArray& fara, bool& first, double* dot, Vector* axis, int pid, double acc)
{
	// Add two types of dot products
	// 1. planar
	// 2. with a fixed value (unknown in the beginning)
	assert(f->processed()==0);
	f->processed(pid);
	const double dotplane = 1.0 - eps;		// close to 1.0 with the mesh accuracy
	bool added = false;
	for (int i=0; i<3; i++) {
		Face* n = f->neighbor(i);
		if (n->processed()==0) {
			double d = f->normal() * n->normal();
			if (d<dotplane) {
				//cout << "d=" << d << " dot=" << *dot << " first=" << first << " a*n=" << *axis * n->normal() << endl;
				if (d<*dot-acc) continue;
				if (first) {
					*dot  = d;
					*axis = (f->normal() ^ n->normal()).unit();
					first = false;
				} else
				//if (Abs(d-*dot) > acc || Abs(*axis * n->normal()) > acc)
				//	continue;
				if (Abs(d-*dot) > acc)
					continue;
			}
			fara << n;
			addNeighborsCylinder(n, fara, first, dot, axis, pid, acc);
			added = true;
		}
	}
	return added;
} // addNeighborsCylinder

/** faces2Vertices given a list of faces
 * @param fara	input array with pointers to faces
 * @param pts	ouput array with pointers to vectors
 *		FIXME could be sorted array in Z for faster lookup
 */
void Mesh::faces2Vertices(FaceArray& fara, VertexArray& pts) const
{
	pts.clear();

	DUMP(cout << endl << "# Mesh::faces2Vertices: faces=" << fara.size() << endl);
	ArrayIterator<Face*> iter(fara);
	while (iter) {
		Face* f = iter++;
		for (int i=0; i<3; i++) {
			Vertex* v = f->vertex(i);
			if (pts.find(v)<0) {
				DUMP(cout << v->x << ' ' << v->y << ' ' << v->z << endl);
				pts << v;
			}
		}
	}
	DUMP(cout << "# Mesh::faces2Vertices: pts=" << pts.size() << endl);
} // faces2Vertices

#if _DEBUG>1
/** writeArray2Stl */
static bool writeArray2Stl(const char* filename, FaceArray& fara)
{
	char buf[80];
	std::ofstream stream(filename, ios::binary);

	strncpy(buf, filename, sizeof(buf));
	stream.write(buf, sizeof(buf));

	dword triangles = fara.size();
	stream.write(reinterpret_cast<char*>(&triangles), sizeof(triangles));

	float triangle[4*3];
	triangle[0] = 0.0;	// normal
	triangle[1] = 0.0;
	triangle[2] = 0.0;

	for (int i=0; i<fara.size(); i++) {
		Face* f = fara[i];
		for (int j=0; j<3; j++) {
			const Vertex* v = f->vertex(j);
			triangle[3+j*3 + 0] = (float)v->x;
			triangle[3+j*3 + 1] = (float)v->y;
			triangle[3+j*3 + 2] = (float)v->z;
		}
		stream.write(reinterpret_cast<char*>(triangle), sizeof(triangle));

		word attr = 0;
		stream.write(reinterpret_cast<char*>(&attr), sizeof(attr));
	}

	stream.close();
	return true;
} // writeArray2Stl
#endif

/** _planeError */
static double _planeError(VertexArray& pts, double a, double b, double c, double d)
{
	double err = 0.0;
	for (int i=0; i<pts.size(); i++)
		err += Sqr(a*pts[i]->x + b*pts[i]->y + c*pts[i]->z + d);

	return sqrt(err/(Sqr(a)+Sqr(b)+Sqr(c))) / (double)pts.size();
} // _planeError

/* fitPlane from a list of vectors
 * @param pts		a list of points that form the plane
 * @param a,b,c,d	return plane parameters
 * @param eps		accuracy of the fit
 * @return		true if fit was possible
 */
static bool fitPlane(VertexArray& pts, double* a, double* b, double* c, double* d, double* err, const double eps)
{
	double Sx, Sy, Sz;
	double Sx2, Sy2, Sz2;
	double Sxy, Sxz, Syz;
	double Vx, Vy, Vz;

	// initialize variables
	*a = *b = *c = 0.0; *d = 1.0;

	// we need minimum 3 vertices
	int n = pts.count();
	if (n < 3) return false;

	// First do statistics with points
	Sx  = Sy  = Sz  = 0.0;
	Sx2 = Sy2 = Sz2 = 0.0;
	Sxy = Sxz = Syz = 0.0;
	for (int i=0; i<n; i++) {
		Sx  += pts[i]->x;
		Sy  += pts[i]->y;
		Sz  += pts[i]->z;

		Sx2 += Sqr(pts[i]->x);
		Sy2 += Sqr(pts[i]->y);
		Sz2 += Sqr(pts[i]->z);

		Sxy += pts[i]->x * pts[i]->y;
		Sxz += pts[i]->x * pts[i]->z;
		Syz += pts[i]->y * pts[i]->z;
	}
	Sx /= (double)n;		// mean value
	Sy /= (double)n;
	Sz /= (double)n;
	Vx = Sx2/(double)n - Sqr(Sx);	// variance
	Vy = Sy2/(double)n - Sqr(Sy);
	Vz = Sz2/(double)n - Sqr(Sz);

	// Count zero variances
	int nv = (int)Eq0(Vx, eps);
	if (Eq0(Vy, eps)) nv++;
	if (Eq0(Vz, eps)) nv++;
	if (nv) {
		if (nv>1) return false;
		// Planes parallel to axes
		// Try the solution of x=Xo or y=Yo or z=Zo
		if (Eq0(Vx, eps)) {
			*a = 1.0;
			*d = -Sx;
			*err = _planeError(pts, *a, *b, *c, *d);
			return true;
		} else
		if (Eq0(Vy, eps)) {
			*b = 1.0;
			*d = -Sy;
			*err = _planeError(pts, *a, *b, *c, *d);
			return true;
		} else {
			*c = 1.0;
			*d = -Sz;
			*err = _planeError(pts, *a, *b, *c, *d);
			return true;
		}
	}

	/*
	 * Try a generic solution
	 *  z = ax + by + d    <=>  ax + by -z + d = 0
	 *  assuming c=-1
	 *  it can only fail on ax + by + d = 0
	 *
	 *  / Sx2    Sxy    Sx \       / Sxz \
	 *  | Sxy    Sy2    Sy | * X = | Syz |
	 *  \ Sx     Sy     n  /       \ Sz  /
	 */
	Matrix3 A;
	Vertex B;

	A(0,0) = Sx2; A(0,1) = Sxy; A(0,2) = Sx;          B.x = Sxz;
	A(1,0) = Sxy; A(1,1) = Sy2; A(1,2) = Sy;          B.y = Syz;
	A(2,0) = Sx;  A(2,1) = Sy;  A(2,2) = (double)n;   B.z = Sz;

	cout << "A=" << endl << A << endl << "det=" << A.det() << endl;
	if (A.inverse(eps)) {
		Vector X = A*B;
		cout << "X=" << X << endl;
		*a = X.x;
		*b = X.y;
		*c =-1.0;
		*d = X.z;
		*err = _planeError(pts, *a, *b, *c, *d);
		return true;
	}

	/* Try a solution where c=0
	 * y = ax + d  <=>   ax -y +d = 0
	 *
	 *  / Sx2    Sx \       / Sxy \
	 *  |           | * X = |     |
	 *  \ Sx     n  /       \ Sy  /
	 */
	Matrix2 A2;
	Vertex  B2;

	A(0,0) = Sx2; A(0,1) = Sx;          B.x = Sxy;
	A(1,0) = Sx;  A(1,1) = (double)n;   B.y = Sy;

	cout << "A2=" << endl << A2 << endl << "det=" << A2.det() << endl;
	if (A2.inverse(eps)) {
		Vector X = A2*B2;
		*a = X.x;
		*b =-1.0;
		*c = 0.0;
		*d = X.y;
		*err = _planeError(pts, *a, *b, *c, *d);
		return true;
	}
	return false;

} // fitPlane

/** Convert a generic plane to a more appropriate form
 * @param a,b,c,d	plane parameters
 * @return a ("PLA", "XYP", "XZP", "YZP", [x|y|z|(nx,ny,nz,px,py,pz)], err
 */
static BodyType plane2Fluka(double a, double b, double c, double d, double *what)
{
	if (Eq0(a,FLOATSMALL)) a = 0.0;
	if (Eq0(b,FLOATSMALL)) b = 0.0;
	if (Eq0(c,FLOATSMALL)) c = 0.0;
	if (Eq0(d,FLOATSMALL)) d = 0.0;

	// prefer the positive solution
	if (a<=FLOATSMALL && b<=FLOATSMALL && c<=FLOATSMALL) {
		a = -a;
		b = -b;
		c = -c;
		d = -d;
	}

	for (int i=0; i<6; i++) what[i] = 0.0;

	if (Eq0(a,FLOATSMALL) && Eq0(b,FLOATSMALL) && !Eq0(c,FLOATSMALL)) {
		what[0] = -d/c;
		if (Eq0(what[0],FLOATSMALL)) what[0] = 0.0;
		return XYPbody;
	} else
	if (Eq0(a,FLOATSMALL) && !Eq0(b,FLOATSMALL) && Eq0(c,FLOATSMALL)) {
		what[0] = -d/b;
		if (Eq0(what[0],FLOATSMALL)) what[0] = 0.0;
		return XZPbody;
	} else
	if (!Eq0(a,FLOATSMALL) && Eq0(b,FLOATSMALL) && Eq0(c,FLOATSMALL)) {
		what[0] = -d/a;
		if (Eq0(what[0],FLOATSMALL)) what[0] = 0.0;
		return YZPbody;
	}

	// find a nice position
	// search intersections of the plane with XY, YZ, ZX axis planes
	// or the closest distance
	// for the most simpler representation
	// FIXME we could have used from the list of vertices a vertex
	// along the plane with the minimum coordinates (smallest distance to center)
	char str[256];
	int minlen = 999999;
	double abc = a+b+c;
	if (!Eq0(abc, FLOATSMALL)) {
		what[3] = what[4] = what[5] = -d/abc;
		sprintf(str,"%.15g %.15g %.15g", what[3], what[4], what[5]);
		minlen = strlen(str);
	}

	if (!Eq0(a, FLOATSMALL)) {
		sprintf(str,"%.15g 0.0 0.0", -d/a);
		int len = strlen(str);
		if (len < minlen) {
			minlen = len;
			what[3] = -d/a;
			what[4] = 0.0;
			what[5] = 0.0;
		}
	}

	if (!Eq0(b, FLOATSMALL)) {
		sprintf(str,"0.0 %.15g 0.0", -d/b);
		int len = strlen(str);
		if (len < minlen) {
			minlen = len;
			what[3] = 0.0;
			what[4] = -d/b;
			what[5] = 0.0;
		}
	}

	if (!Eq0(c, FLOATSMALL)) {
		sprintf(str,"0.0 0.0 %.15g", -d/c);
		int len = strlen(str);
		if (len < minlen) {
			what[3] = 0.0;
			what[4] = 0.0;
			what[5] = -d/c;
		}
	} // FIXME check that MUST have a solution!!!!!!!

	// normalize normal
	double s = sqrt(a*a + b*b + c*c);
	what[0] = a / s;
	what[1] = b / s;
	what[2] = c / s;

	for (int i=0; i<6; i++)
		if (Eq0(what[i],FLOATSMALL))
			what[i] = 0.0;

	return PLAbody;
} // convert2FlukaPlane

/** fit circle using indices ia and ib of points
 * @return true on success
 */
static bool fitCircle(VertexArray& pts, int ia, int ib, double* a, double* b, double* r, double* err)
{
	/*
	 * (x-xc)^2 + (y-yc)^2 = R^2
	 *
	 * x^2 + y^2 - 2xc x - 2yc y - R^2 + xc^2 + yc^2 = 0
	 *
	 * 1*(xc^2+yc^2) + 2*x*xc
	 */
	int n=pts.size();
	if (n<3) return false;

	// solution matrices
	Matrix A(n,3);
	Matrix B(n,1);
	Matrix X;

	// populate matrices
	for (int i=0; i<n; i++) {
		double x = (*pts[i])[ia];
		double y = (*pts[i])[ib];
		A(i,0) = 2.0*x;
		A(i,1) = 2.0*y;
		A(i,2) = 1.0;
		B(i,0) = x*x + y*y;
	}

	if (solveOverDetermined(A,B,X)) {
		*a = X(0,0);
		*b = X(1,0);
		double c = X(2,0);
		*r = sqrt(Sqr(*a) + Sqr(*b) + c);
		if (*r<0.0) return false;

		// find error
		*err = 0.0;
		for (int i=0; i<n; i++) {
			double x = (*pts[i])[ia];
			double y = (*pts[i])[ib];
			*err += Sqr( Sqr(x - *a)
			           + Sqr(y - *b)
				   - Sqr(*r));
		}
		*err /= (double)n;
		*err = sqrt(*err);
		return true;
	} else
		return false;
} // fitCircle

#if 0
/** fitSphere */
static bool fitSphere(VertexArray& pts, double *x, double *y, double *z, double *r, double *err, const double /*eps*/)
{
	/*
	 * (x-xc)^2 + (y-yc)^2 + (z-zc)^2 - R^2 = 0
	 *
	 * x^2 + y^2 + z^2 - 2xc x - 2yc y - 2zc z - R^2 + xc^2 + yc^2 + zc^2 = 0
	 */
	int n=pts.size();
	if (n<5) return false;

	// solution matrices
	Matrix A(n,4);
	Matrix B(n,1);
	Matrix X;
	// populate matrices
	for (int i=0; i<n; i++) {
		A(i,0) = 2.0*pts[i]->x;
		A(i,1) = 2.0*pts[i]->y;
		A(i,2) = 2.0*pts[i]->z;
		A(i,3) =-1.0;
		B(i,0) = pts[i]->x*pts[i]->x + pts[i]->y*pts[i]->y + pts[i]->z*pts[i]->z;
	}

	if (solveOverDetermined(A,B,X)) {
		*x = X(0,0);
		*y = X(1,0);
		*z = X(2,0);
		double c = X(3,0);
		*r = sqrt(-c + Sqr(*x) + Sqr(*y) + Sqr(*z));
		if (*r<0.0) return false;

		// find error
		*err = 0.0;
		for (int i=0; i<n; i++) {
			*err += Sqr( Sqr(pts[i]->x - *x)
			           + Sqr(pts[i]->y - *y)
			           + Sqr(pts[i]->z - *z)
				   - Sqr(*r));
		}
		*err /= (double)n;
		*err = sqrt(*err);
		return true;
	} else
		return false;

} // fitSphere

/** coef
 * @param p	vector
 * @param i	return ith coefficient of vector
 * @return coefficient of vector p
 */
static double coef(const Vertex& p, const int i)
{
	switch (i) {
		case 0:	return p.x * p.x;
		case 1:	return p.y * p.y;
		case 2:	return p.z * p.z;

		case 3:	return p.x * p.y;
		case 4:	return p.x * p.z;
		case 5:	return p.y * p.z;

		case 6:	return p.x;
		case 7:	return p.y;
		case 8:	return p.z;

		default:
			return 1.0;
	}
} // coef

/** fitQuadratic */
static bool fitQuadratic(VertexArray& pts, double C[10], double *err, const double eps)
{
	int n=pts.size() ;
	if (n<9) return false;

	// Initialize variables
	for (int j=0; j<10; j++) C[j] = 0.0;

	// Make statistics to see which variables are not zero
	double S[9], S2[9], V[9];
	for (int i=0; i<9; i++) S[i] = S2[i] = V[i] = 0.0;

	for (int i=0; i<n; i++)
		for (int j=0; j<9; j++) {
			double c = coef(*pts[i],j);
			S[j] += c;
			S2[j] += Sqr(c);
		}

	// Calculate variance
	for (int j=0; j<9; j++) {
		S[j]  /= (double)n;
		S2[j] /= (double)n;
		V[j] = S2[j] - Sqr(S[j]);
		cout << j << " Mean=" << S[j] << "   Var=" << V[j] << endl;
	}

	// Count non-zero terms
	int nv = 0;
	for (int j=0; j<9; j++)
		if (Eq0(V[j],eps)) nv++;

	// Assume C != 0.0
	// Create a matrix with the non-zero terms
	Matrix A(n, 9-nv);
	Matrix B(n, 1);
	Matrix X;

	// populate matrices
	C[9] = 1.0;
	for (int i=0; i<n; i++) {
		int col = 0;
		for (int j=0; j<9; j++) {
			if (Eq0(V[j],eps)) continue;
			A(i,col) = coef(*pts[i], j);
			col++;
		}
		B(i,0) = -1;
	}
	if (solveOverDetermined(A,B,X)) {
		int col = 0;
		for (int j=0; j<9; j++) {
			if (Eq0(V[j],eps)) continue;
			C[j] = X(col,0);
			col++;
		}
		goto _FINDERR;
	}

	// In case of non-convergence
	// from C[9] to C[1] set to 0.0 progressively
	// and check for a solution
//	C[9] = 0.0;
//	for (int k=8; k>=1; k--) {
//		if (Eq0(V[k],eps)) continue;
//		C[k] = 1.0;
//		A.make(n, 9-nv-k); //???
//		C[k] = 0.0;
//	}

	return false;
_FINDERR:
	*err = 0.0;
	for (int i=0; i<n; i++) {
		double q=0.0;
		for (int j=0; j<10; j++)
			q += C[j]*coef(*pts[i],j);
cout << i << " q=" << q << endl;
		*err += Sqr(q);
	}
	*err /= (double)n;
	*err = sqrt(*err);
	return true;
} // fitQuadratic

/** fitEllipticalCylinder*/
static bool fitEllipticalCylinder(VertexArray& pts, double *centre1, double *centre2, double *ray1, double *ray2, double *err)
{
	int n=pts.count();
	if (n<4) return false;
	Matrix AX(n,4);
	Matrix AY(n,4);
	Matrix AZ(n,4);
	Matrix B(n,1);

	for (int i=0; i<n; i++)
	{
		AX(i,0) = pts[i]->y*pts[i]->y;
		AX(i,1) = pts[i]->z*pts[i]->z;
		AX(i,2) = -2.0*pts[i]->y;
		AX(i,3) = -2.0*pts[i]->z;

		AY(i,0) = pts[i]->x*pts[i]->x;
		AY(i,1) = pts[i]->z*pts[i]->z;
		AY(i,2) = -2.0*pts[i]->x;
		AY(i,3) = -2.0*pts[i]->z;

		AZ(i,0) = pts[i]->x*pts[i]->x;
		AZ(i,1) = pts[i]->y*pts[i]->y;
		AZ(i,2) = -2.0*pts[i]->x;
		AZ(i,3) = -2.0*pts[i]->y;

		B(i,0) = 1.0;
	}
/*
	cout << "AX" << AX << endl;
	cout << "AY" << AY << endl;
	cout << "AZ" << AZ << endl;
*/
	double a=0,b=0,x=0,y=0,z=0,/*c=0,*/Ry=0,Rz=0,Rx=0;
	Matrix result(4,1);
	//"XEC"

	if (solveOverDetermined(AX,B,result)){
		a = result(0,0);
		b = result(1,0);
		y = result(2,0);
		z = result(3,0);
		if( a!=0 && b!=0){
			y/=a; z/=b;
			Ry = sqrt(y*y + b/a * z*z + 1.0/a);
			Rz = sqrt(a/b)*Ry;
			*centre1=y;
			*centre2=z;
			*ray1=Ry;
			*ray2=Rz;
			// find error
			for (int i=0; i<n; i++)
				*err += (*ray2)*(*ray2)*(Sqr(pts[i]->y - (*centre1))) +
					(*ray1)*(*ray1)*(Sqr(pts[i]->z - (*centre2)));
			*err = *err-n*(*ray2)*(*ray2)*(*ray1)*(*ray1);
			if ((*err)>(*ray1)) goto _YEC;
			cout<<"XEC ";
			return true;
		}
	}
	_YEC:
	if (solveOverDetermined(AY,B,result)){
		a = result(0,0);
		b = result(1,0);
		x = result(2,0);
		z = result(3,0);
		if( a!=0 && b!=0){
			x/=a; z/=b;
			Rx = sqrt(x*x + b/a * z*z + 1.0/a);
			Rz = sqrt(a/b)*Rx;
			*centre1=x;
			*centre2=z;
			*ray1=Rx;
			*ray2=Rz;
			// find error
			*err=0;
			for (int i=0; i<n; i++)
				*err += (*ray2)*(*ray2)*(Sqr(pts[0]->x - (*centre1))) +
					(*ray1)*(*ray1)*(Sqr(pts[0]->z - (*centre2)));
			*err = *err-n*(*ray2)*(*ray2)*(*ray1)*(*ray1);
			if ((*err)>(*ray1)) goto _ZEC;
			cout<<"YEC ";
			return true;
		}
	}
	_ZEC:
	if (solveOverDetermined(AZ,B,result)){
		a = result(0,0);
		b = result(1,0);
		x = result(2,0);
		y = result(3,0);
		if( a!=0 && b!=0){
			x/=a; y/=b;
			Rx = sqrt(x*x + b/a * y*y + 1.0/a);
			Ry = sqrt(a/b)*Rx;
			*centre1=x;
			*centre2=y;
			*ray1=Rx;
			*ray2=Ry;
			cout<<"ZEC ";
			// find error
			*err=0;
			for (int i=0; i<n; i++)
				*err += (*ray2)*(*ray2)*(Sqr(pts[i]->x - (*centre1))) +
					(*ray1)*(*ray1)*(Sqr(pts[i]->y - (*centre2)));
			*err = *err-n*(*ray2)*(*ray2)*(*ray1)*(*ray1);
			return true;
		}
	}
	return false;
} // fitEllipticalCylinder
#endif

/** @return memory used by Mesh */
size_t Mesh::memory() const
{
	return  sizeof(*this)
	      + _vertices.memory()
	      + _edges.memory()
	      + _faces.memory()
	      + sizeof(Vertex) * nvertices()
	      + sizeof(Edge)  * nedges()
	      + sizeof(Face)  * nfaces();
} // memory

/** createParallelepiped */
void Mesh::createParallelepiped(Mesh& mesh, const Vertex& P, const Vector& X, const Vector& Y, const Vector& Z)
{
	/*         Z
	 *         |
	 *         4------------7
	 *        /|           /|
	 *       5------------6 |
	 *       | |          | |
	 *       | |          | |
	 *       | | O        | |
	 *       | 0----------|-3--> Y
	 *       |/           |/
	 *       1------------2
	 *      /
	 *     X
	 */
	mesh.resize(8, 2*5+8, 12);
	mesh.allocateVertices(8);

	// Add vertices (without checking for the same!
	mesh.vertex(0) = P;
	mesh.vertex(1) = P + X;
	mesh.vertex(2) = P + X + Y;
	mesh.vertex(3) = P     + Y;

	mesh.vertex(4) = P         + Z;
	mesh.vertex(5) = P + X     + Z;
	mesh.vertex(6) = P + X + Y + Z;
	mesh.vertex(7) = P + Y     + Z;

	// Add faces
	mesh.add(0, 2, 1, false, true, true );
	mesh.add(0, 3, 2, true,  true, false);

	mesh.add(0, 5, 4, false, true, true );
	mesh.add(0, 1, 5, true,  true, false);

	mesh.add(1, 6, 5, false, true, true );
	mesh.add(1, 2, 6, true,  true, false);

	mesh.add(2, 7, 6, false, true, true );
	mesh.add(2, 3, 7, true,  true, false);

	mesh.add(3, 4, 7, false, true, true );
	mesh.add(3, 0, 4, true,  true, false);

	mesh.add(4, 6, 7, false, true, true );
	mesh.add(4, 5, 6, true,  true, false);

	mesh.calcBbox();
	mesh.process();
	assert(mesh.isClosed());
	assert(mesh.isOrientable());
#if _DEBUG>2
	cout << "BOX Mesh ";
	cout << " isClosed=" << mesh.isClosed();
	cout << " isOrientable=" << mesh.isOrientable() << endl;
	cout << "BOX volume=" << mesh.volume() << endl;
#endif
} // createParallelepiped

/** Print mesh. */
ostream& operator << (ostream& s, Mesh& mesh)
{
	s << "Mesh" << endl;
	s << "    vertices: " << mesh.nvertices() << endl;
	for (int i=0; i<mesh.nvertices(); i++)
		s << "\t" << i << ": " << mesh.vertex(i) << endl;
	s << "    edges: " << mesh.nedges() << endl;
	for (int i=0; i<mesh.nedges(); i++) {
		Edge* e = mesh.edge(i);
		int a = mesh.findVertex(*e->a);
		int b = mesh.findVertex(*e->b);
		s << "\t" << i << ": [" << a << ", " << b << "]" << endl;
	}
	s << "    faces: " << mesh.nfaces() << endl;
	for (int i=0; i<mesh.nfaces(); i++) {
		Face* f = mesh.face(i);
		s << "\t" << i << ": " << *f << "\tbend=";
		for (int j=0; j<3; j++) {
			Face* fn = f->neighbor(j);
			if (fn)
				cout << " " << f->normal() * fn->normal();
			else
				cout << " ? ";
		}
		cout << endl;
	}
	return s;
} /* operator << */
