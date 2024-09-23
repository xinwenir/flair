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

#include <iostream>

#include "bmath.h"
#include "vbody.h"
#include "viewport.h"

using namespace std;

/* ============================== VBody ================================ */
/* Describe a convex body with up to 6 intersecting quads */

VBody::VBody(GBody *abody)
{
#ifdef THREAD
	pthread_mutex_init(&mutex, NULL);
#endif
	for (int i=0; i<BODYCONICS; i++) {
		V[i].delta(16);
		V[i].compare(Vertex2D::compare);
	}
	init(abody);
} // VBody

/** ~VBody */
VBody::~VBody()
{
#ifdef THREAD
	pthread_mutex_destroy(&mutex);
#endif
	delVertices();
	_body=NULL;
} /* ~VBody */

/** init */
void VBody::init(GBody *abody)
{
	_body       = abody;

	nC          = 0;
	visible     = false;
	notref      = false;
	_generation = -1;
	location    = LOCATION_OVERLAP;	// default

	for (int i=0; i<BODYQUADS; i++) acc[i] = 1.0;
	delVertices();
} // init

/** delVertices */
void VBody::delVertices()
{
	for (int i=0; i<BODYCONICS; i++)
		V[i].clear();
} // delVertices

/* makeConics
 * Create viewport conics of the body and initializes location information
 * @param view	viewport
 */
void VBody::makeConics(const ViewPort& view)
{
	// clear settings
	delVertices();
	nC     = 0;
	notref = false;
	location = LOCATION_OVERLAP;	// default

	if (type() == ERRbody) {
		location = LOCATION_OUTSIDE;
		return;
	}

	// Create transformed quadratics
	Quad QN[BODYQUADS];
	for (int i=0; i<nQ(); i++) {
		QN[i] = Q(i);
		QN[i].transform(view.matrix());
		QN[i].normalize();
#if defined(_DUMP) && _DEBUG>1
		if (show()) {
			cout << "+++ VBody::makeConics:transform " << name() << ":" << i << endl;
			cout << "*-*\t Q[" << i << "] = " << Q(i) << endl;
			cout << "*-*\t QN["<< i << "] = " << QN[i] << endl;
			cout << "*-*\t M=" << endl << view.matrix() << endl;
		}
#endif
	}

	int i;
	double x = view.Uofs();
	double y = view.Vofs();

	for (i=0; i<nQ(); i++) acc[i] = 1.0;

#if _DEBUG>1
//	if (!strcmp(name(),"parabola"))
//		cout << "break\n";
#endif

	// Create conics
	for (i=0; i<nQ(); i++) {
		const Quad& QNi = QN[i];
		Conic conic(QNi.Cxx, QNi.Cxy/2.0, QNi.Cyy, QNi.Cx/2.0, QNi.Cy/2.0, QNi.C);
//fprintf(stderr,"Conic.%s.%d: %.22lg %.22lg %.22lg %.22lg %.22lg %.22lg %s\n",
//	name(), i,
//	conic.a, conic.h, conic.b, conic.g, conic.f, conic.c, conic.typeStr());
		DUMP(if (show()) cout << "+++ " << QNi << endl);
		DUMP(if (show()) cout << "+++ " << conic << endl);

		if (conic.type()==CONIC_LINES) {
			Conic c1, c2;
			conic.splitLines(&c1, &c2);
			c1.parametric();
			c2.parametric();
			DUMP(if (show()) cout << "+1+ " << c1 << endl);
			DUMP(if (show()) cout << "+2+ " << c2 << endl);

			bool in1 = view.inside(c1);
			bool in2 = view.inside(c2);

			// if any of the two is inside then ok
			if (in1 || in2) {
				if (c1.equal(c2,SMALL)) {
					// Conics are the same, so we are tangential
					// to the quad(s) decide with the normal
					// if we are in or out
					location = LOCATION_OUTSIDE;
					nC = 0;
					return;
				}
				if (in1) {
					c2q[nC] = i;
					C[nC++] = c1;
				}
				if (in2) {
					c2q[nC] = i;
					C[nC++] = c2;
				}
			} else {
				double q = QNi(x,y);
				double gx, gy, gz, g;
				QNi.grad(x, y, &gx, &gy, &gz);
				g = sqrt(gx*gx + gy*gy + gz*gz);
				if (g > SMALL) q /= g;
				DUMP(if (show()) cout << "\tq=" << q << endl);
				DUMP(if (show()) cout << "\tacc=" << QNi.acc(x,y,Conic::CONICPREC) << endl);
				if (q > QNi.acc(x,y,Conic::CONICPREC)) {
					location = LOCATION_OUTSIDE;
					nC = 0;
					DUMP(if (show()) cout << "\tlocation="
						<< location << " q=" << QNi(x,y) << endl);
					return;
				}
			}
		} else
		if (conic.type()==CONIC_LINE      || conic.type()==CONIC_ELLIPSE ||
		    conic.type()==CONIC_HYPERBOLA || conic.type()==CONIC_PARABOLA) {
			conic.parametric();
			if (view.inside(conic)) {
				c2q[nC] = i;
				C[nC++] = conic;
			} else {
				double q = QNi(x,y);
				double gx, gy, gz, g;
				QNi.grad(x, y, &gx, &gy, &gz);
				g = sqrt(gx*gx + gy*gy + gz*gz);
				if (g > SMALL) q /= g;
				DUMP(if (show()) cout << "\tq=" << q << endl);
				DUMP(if (show()) cout << "\tacc=" << QNi.acc(x,y,Conic::CONICPREC) << endl);
				if (q > QNi.acc(x,y,Conic::CONICPREC)) {
					location = LOCATION_OUTSIDE;
					nC = 0;
					DUMP(if (show()) cout << "\tlocation="
						<< location << " q=" << QNi(x,y) << endl);
					return;
				}
			}
		} else {
			DUMP(if (show()) cout << "\tq=" << QNi(x,y) << endl);
			double q = QNi(x,y);
			double gx, gy, gz, g;
			QNi.grad(x, y, &gx, &gy, &gz);
			g = sqrt(gx*gx + gy*gy + gz*gz);
			if (g > SMALL) q /= g;
			if (q > QNi.acc(x,y,Conic::CONICPREC)) {
				location = LOCATION_OUTSIDE;
				nC = 0;
				DUMP(if (show()) cout << "\tlocation=" << location << " q=" << QNi(x,y) << endl);
				return;
			}
		}
	}
} // makeConics

/** Intersect all the conics of a body
 * @param view	viewport
 */
void VBody::intersectSelf(const ViewPort& view)
{
	for (int i=0; i<nC; i++) {
		const Conic& aconic = C[i];
		if (aconic.type()==CONIC_ELLIPSE || aconic.type()==CONIC_HYPERBOLA) {
			// For ellipse add starting and ending point
			double x, y;
			aconic.getXY(PI, &x, &y);
			// XXX FIXME possible bug if it is outside the body
			// while at half distance there is another intersecting body!!!
			double rx = view.uv2x(x,y);
			double ry = view.uv2y(x,y);
			double rz = view.uv2z(x,y);
			if (view.inside(x,y) && inside2D(rx, ry, rz,
					-view.matrix(0,2), -view.matrix(1,2), -view.matrix(2,2),
					c2q[i]))
			{
				V[i].add(Vertex2D(-PI, x, y, this));
				V[i].add(Vertex2D( PI, x, y, this));
			}
		}
		for (int j=i+1; j<nC; j++) {
			const Conic& bconic = C[j];
			Vector2D pts[4];
			int n = aconic.intersect(bconic, pts);
			for (int k=0; k<n; k++) {
				double x = pts[k].x;
				double y = pts[k].y;
				double rx = view.uv2x(x,y);
				double ry = view.uv2y(x,y);
				double rz = view.uv2z(x,y);
				if (view.inside(x,y) &&
					inside2D(rx,ry,rz,
						-view.matrix(0,2), -view.matrix(1,2), -view.matrix(2,2),
						c2q[i], c2q[j]))
				{
					double t = aconic.getT(x,y);
					V[i].add(Vertex2D(t,x,y, this));
					t = bconic.getT(x,y);
					V[j].add(Vertex2D(t,x,y, this));
				}
			}
		}
	}
} // intersectSelf

/* bodyIntersections
 * Intersect a body with itself + viewport limits
 * if there is no intersection it updates the location information
 * @param view	viewport
 */
void VBody::intersectViewport(const ViewPort& view)
{
	if (!nC) return;
	// Intersect with the window
	intersectConic(view.ulow,  view);
	intersectConic(view.uhigh, view);
	intersectConic(view.vlow,  view);
	intersectConic(view.vhigh, view);
} // intersectViewport

/* updateLocation
 * if a body has no vertices means it will be either inside or outside the
 * viewport, update this information
 */
void VBody::updateLocation(const ViewPort& view)
{
	// No conics and location already determined from makeConics
	if (!nC && location != LOCATION_OVERLAP)
		return;

	// Check if body lies outside the plotting window
	for (int i=0; i<nC; i++)
		if (!V[i].empty()) return;

	// No intersection found, remove conics
	nC = 0;

	// find location of body, check several points on screen
	double dx = view.relWidth() / 4.0;
	double dy = view.relHeight() / 4.0;

	// 3 non-linear points should be sufficient
	// 1 or 2 might end-up on a degenerated conic of line or point
	// example a cylinder tangential to the plot will be only one line
	double x,y,z;
	view.origin(&x,&y,&z);

	assert(location == LOCATION_OVERLAP);
	if (!inside2D(x, y, z,
			-view.matrix(0,2), -view.matrix(1,2), -view.matrix(2,2))) {
		location = LOCATION_OUTSIDE;
		return;
	}

	double rdx = view.uv2x(dx,dy);
	double rdy = view.uv2y(dx,dy);
	double rdz = view.uv2z(dx,dy);

	if (!inside2D(rdx, rdy, rdz,
			-view.matrix(0,2), -view.matrix(1,2), -view.matrix(2,2))) {
		location = LOCATION_OUTSIDE;
		return;
	}

	rdx = view.uv2x(dx,-dy);
	rdy = view.uv2y(dx,-dy);
	rdz = view.uv2z(dx,-dy);

	if (!inside2D(rdx, rdy, rdz,
			-view.matrix(0,2), -view.matrix(1,2), -view.matrix(2,2))) {
		location = LOCATION_OUTSIDE;
		return;
	}

	// set to inside
	location = LOCATION_INSIDE;

	DUMP(if (show()) cout << "\tlocation=" << location << endl);
} // updateBodyLocation

/** Intersect the conics of two bodies and keep vertices inside the view
 * @param abody	first body
 * @param view	viewport
 */
bool VBody::intersectBody(VBody* abody, const ViewPort& view)
{
	assert(abody);
	assert(this != abody);

	bool intersected = false;

#if _DUMP && _DEBUG>1
	if (show() && abody->show()) {
		cout << "VBody::intersectBody" << endl;
		cout << "\tname=" << name()
		<< " abody->name=" << abody->name() << endl;
	}
#endif
	// Loop over all conics of this
	for (int i=0; i<nC; i++) {
		// with all conics of abody
		for(int j=0; j<abody->nC; j++) {
			Vector2D pts[4];

			int n = C[i].intersect(abody->C[j], pts);

			for (int k=0; k<n; k++) {
				double x  = pts[k].x;
				double y  = pts[k].y;
				double rx = view.uv2x(x,y);
				double ry = view.uv2y(x,y);
				double rz = view.uv2z(x,y);
				if (view.inside(x,y) &&
				    inside2D(rx,ry,rz,
						-view.matrix(0,2), -view.matrix(1,2), -view.matrix(2,2),
						c2q[i]) &&
				    abody->inside2D(rx,ry,rz,
						-view.matrix(0,2), -view.matrix(1,2), -view.matrix(2,2),
						abody->c2q[j]))
				{
					addVertex(i, x,y, abody);
					abody->addVertex(j, x,y, this);
					intersected = true;
				}
			}
		}
	}
	return intersected;
} // intersectBody

/** intersect body with a conic (normally from the viewport) */
void VBody::intersectConic(const Conic& conic, const ViewPort& view)
{
	Vector2D pts[4];

	for (int i=0; i<nC; i++) {
		int n = C[i].intersect(conic, pts);
		for (int j=0; j<n; j++) {
			double x  = pts[j].x;
			double y  = pts[j].y;
			double rx = view.uv2x(x,y);
			double ry = view.uv2y(x,y);
			double rz = view.uv2z(x,y);
			if (view.inside(x,y) &&
			    inside2D(rx,ry,rz, -view.matrix(0,2), -view.matrix(1,2), -view.matrix(2,2), c2q[i])) {
				double t = C[i].getT(x,y);
				V[i].add(Vertex2D(t,x,y, NULL));
			}
		}
	}
} // intersectConic

/** removeInvalidVertices
 * remove vertices that belong to invalid bodies and invalidate
 * This is tricky because removing a vertex in a segment does change the
 * region checks performed by scanSegments
 */
void VBody::removeInvalidVertices()
{
	assert(valid());
	for (int i=0; i<nC; i++) {
		Array<Vertex2D>& VV = V[i];
		for (int j=VV.size()-1; j>=0; j--) {
			Vertex2D& vertex = VV[j];
			if (vertex.body && vertex.body->invalid()) {
				VV.erase(j);
				// Invalidate touching segments
				if (j<VV.size())
					VV[j].invalid = true;
				if (j>0)
					VV[j-1].invalid = true;
			}
		}
	}
} // removeInvalidVertices

/** markInvalidVertices
 * mark segments touching invalid vertices as invalid
 * This is tricky because adding a new vertex in a segment does change the
 * region checks performed by scanSegments
 */
void VBody::markInvalidVertices(const VBody *invalidBody)
{
	//assert(valid());
	for (int i=0; i<nC; i++) {
		Array<Vertex2D>& VV = V[i];
		for (int j=VV.size()-1; j>=0; j--) {
			const Vertex2D& vertex = VV[j];
			if (vertex.body == invalidBody) {
				// Invalidate touching segments
				if (j<VV.size())
					VV[j].invalid = true;
				if (j+1<VV.size())
					VV[j+1].invalid = true;
			}
		}
	}
} // markInvalidVertices

/** remove problematic vertices from the vertex list
 * Vertex2D should be in the order of intersection as
 * vertex-list:   [self|window ] ...(other bodies)... [self|window]
 * For hyperbola accept it twice
 *
 * XXX FIXME Clean up bodies or the reference to bodies when they are not visible
 *   e.g. nC>2 and all vertex are closer than a certain visible distance!!!
 *   or nC=2 like infinite cylinders or parallel bodies with the intersection to the windows
 *           closer than a certain distance
 * XXX Maybe this check can be done when the body is constructed, even to avoid checking any
 *     intersection with it since it wont be visible!
 */
void VBody::removeWrongVertices()
{
	// for all conics
	for (int i=0; i<nC; i++) {
		if (V[i].size()<2) continue;
		if (C[i].type()==CONIC_HYPERBOLA) continue;
		Array<Vertex2D>& VV = V[i];

		// Find first valid point, intersection with self or viewport
		int j, jp=0;
		double tp = -INFINITE;
		for (j=0; j<VV.size(); j++) {
			Vertex2D& v = VV[j];
			if (!Eq0(v.t-tp, SMALL*(Abs(v.t)+1.0))) {
				jp = j;
				tp = v.t;
			}
			if (v.body==NULL || v.body==this) break;
		}
		if (j==VV.size()) {
			// delete everything!!!!
			VV.clear();
			continue;
		}

		// delete everything from 0 to jp
		if (jp>0) VV.erase(0, jp);

		// Do the same in reverse order
		tp = INFINITE;
		for (j=jp=VV.size()-1; j>=0; j--) {
			Vertex2D& v = VV[j];
			if (!Eq0(v.t-tp, SMALL*(Abs(v.t)+1.0))) {
				jp = j;
				tp = v.t;
			}
			if (v.body==NULL || v.body==this) break;
		}

		// delete everything from jp to end
		VV.erase(jp+1, VV.size());
	}
} // removeWrongVertices

/**
 * Calculate the limits of the accuracy to be used for each quad
 */
void VBody::calculateAccuracy(const ViewPort& view)
{
#define MAXACC	100000.0
	// Initial values
	for (int i=0; i<nQ(); i++)
		acc[i] = SMALL;

	// For each conic
	for (int i=0; i<nC; i++) {
		Array<Vertex2D>& VV = V[i];
		if (VV.empty()) {
			acc[i] = 10.0;
			continue;
		}

		int q = c2q[i]; // Quad index
		bool isHyperbola = C[i].type()==CONIC_HYPERBOLA;

		// Find limits
		double t0 = VV[0].t;
		double tp = t0;

		// calculate in the middle of the segment
		for (int j=1; j<VV.size(); j++) {
			Vertex2D const &v = VV[j];
			if (Eq0(v.t - tp, SMALL*Abs(v.t+SMALL)))
				continue;
			if (!(isHyperbola &&	 // Skip the segment between the two hyperbolas
			      ((tp<-PIover2 && v.t>-PIover2) ||
			       (tp< PIover2 && v.t> PIover2)))) {
				double x, y;
				C[i].getXY((v.t + tp)/2.0, &x, &y);
				double rx = view.uv2x(x,y) - SMALL*view.matrix(0,2);
				double ry = view.uv2y(x,y) - SMALL*view.matrix(1,2);
				double rz = view.uv2z(x,y) - SMALL*view.matrix(2,2);
				acc[q] = Max(acc[q], accuracy(q, rx, ry, rz));
#if _DEBUG>1
				if (accuracy(q, rx, ry, rz) > MAXACC) {
					cout << endl << "ERROR:" << endl;
					cout << "   body: " << name() << "  c=" << i << "  Q=" << q <<  endl;
					cout << "         " << *_body << endl;
					cout << "   quad[" << q << "]=" << Q(q) << endl;
					cout << "   conic[" << i << "]=" << C[i] << endl;
					cout << "   matrix=" << endl << view.matrix() << endl;
					cout << "   pos   = " << rx << " " << ry << " " << rz << endl;
					cout << "   acc   = " << accuracy(q, rx, ry, rz) << endl;
					cout << "   Q()   = " << Q(q)(rx,ry,rz) << endl;
					cout << "   accQ()= " << Q(q).acc(rx,ry,rz, SMALL2) << endl;
				}
#endif
			}
		}
	}
	// Final adjustment
	for (int i=0; i<nQ(); i++) {
		// Allow some margin (x10.0) since all the test will
		// be done slightly outside the viewing plane
		if (acc[i] > MAXACC) {
			if (acc[i] > 10.0*MAXACC)
				cerr << "ERROR: ";
			else
				cerr << "Warning: ";
			cerr << "VBody " << name() << ":" << i << " accuracy=" << acc[i] << endl;
		}
		acc[i] = 64.0 * Range(1.0, acc[i], MAXACC);
#if _DEBUG>1 & defined(_DUMP)
		if (show()) cout << "Accuracy: " << name() << "." << i << " = " << acc[i] << endl;
#endif
	}
} // calculateAccuracy

#if 0
/**
 * Find if a point (x,y) is on the body surface
 * @param x,y	position
 * @param d	distance from the boundary
 * Returns: true  on the surface
 * XXX FIXME XXX problematic routine on RCC!
 * have to check if the point is inside the body on a non-displayed segment
 */
bool VBody::isOn(const double x, const double y, const double d) const
{
	bool on = false;
	for (int i=0; i<nQ(); i++) {
// XXX if Conic::splitLines is corrected I can use the C[i].adist
		double q = Q(i).adist(x,y);
		if (q > d) {
			return false;	// outside
		} else
		if (q >= -d)
			on = true;
	}
	return on;
} // isOn
#endif

/** 2D bounding box on viewer plane
 * @return bounding box of body
 * @WARNING z coordinates are meaningless
 */
BBox VBody::bbox2D() const
{
	BBox bb;
	for (int i=0; i<nC; i++) {
		if (V[i].size()<2) continue;
		bool isLine = C[i].type()==CONIC_LINE;
		const Array<Vertex2D>& VV = V[i];
		double tp = 0.0;
		for (int k=0; k<VV.size(); k++) {
			const Vertex2D& v = VV[k];
			bb.add(v.x, v.y, 0.0);
			if (!isLine && v.t-tp > 1.0) {
				// split it into several parts at least 6
				double tstep = (v.t - tp) / 6.0;
				for (double t=tp+tstep; t<v.t; t += tstep) {
					double x,y;
					C[i].getXY(t,&x,&y);
					bb.add(x, y, 0.0);
				}
			}
			tp = v.t;
		}
	}
	return bb;
} // bbox2D

/** @return memory used by body */
size_t VBody::memory() const
{
	size_t mem = sizeof(VBody);
	for (int i=0; i<nC; i++)
		mem += V[i].memory();
	return mem;
} // memory

/** operator << */
ostream& operator << (ostream& s, const VBody& body)
{
	s << "VBody: " << body.name() << "\t[" << body.typeStr() << "]" << endl;
	for (int i=0; i<body.nQ(); i++)
		s << "   Quad #" << (i+1) << endl << body.Q(i) << endl;
	return s;
} /* operator << */
