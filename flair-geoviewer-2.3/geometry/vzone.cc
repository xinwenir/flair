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

#include <assert.h>
#include <iostream>
#include <ostream>

#include "os.h"
#include "geo.h"
#include "bbox.h"
#include "obbox.h"

#include "vzone.h"
#include "kernel.h"
#include "engine.h"

using namespace std;

#if defined(_DUMP) && _DEBUG>1
	#define RNAME	"R1"
	#define XCHK	3.35680660661
	#define YCHK	-0.866272672673
	#define ACHK	0.1
//	#define CHECK(X) if(!strcmp(region->name(),RNAME)) X
//	#define CHECK(X) if(!strcmp(region->name(),RNAME) && Abs(x-XCHK)<ACHK && Abs(y-YCHK)<ACHK) X
//	#define CHECK(X) if(Abs(x-XCHK)<ACHK && Abs(y-YCHK)<ACHK) X
	#define CHECK(X)	X
#else
	#define CHECK(X)
#endif

/** init */
void VZone::init(GZone* z, VRegion *r)
{
	_zone       = z;
	_region     = r;
	_generation = -1;
	location    = false;
	in          = false;
} // init

/** vBody */
inline VBody* VZone::vbody(const GBody* body) const
{
	return _region->_kernel->getBody(body);
} // vbody

/** cbody */
inline CBody* VZone::cbody(GeometryEngine* engine, const GBody* body)
{
	return engine->getBody(body);
} // cbody

/** updateLocation
 *
 * find the location of the Zone on the Z=0 plane and within the viewport
 * coordinates
 *
 * FIXME: Improvement on STD expressions we can ignore directly the zones
 *        with a negative result, and remove the bodies
 *        Probably add a reference count to the bodies!
 *
 * We have 3 possible locations inside(1) / outside(0) / unknown(?:2)
 *		  a + b
 *	0 + 0 = 0	? + 0 = 0
 *	0 + 1 = 0	? + 1 = ?
 *	1 + 0 = 0	0 + ? = 0
 *	1 + 1 = 1	1 + ? = ?
 *
 *		  a - b
 *	0 - 0 = 0	? - 0 = ?
 *	0 - 1 = 0	? - 1 = 0
 *	1 - 0 = 1	0 - ? = 0
 *	1 - 1 = 0	1 - ? = ?
 *
 *		  a | b
 *	0 | 0 = 0	? | 0 = ?
 *	0 | 1 = 1	? | 1 = 1
 *	1 | 0 = 1	0 | ? = ?
 *	1 | 1 = 1	1 | ? = 1
 */
void VZone::updateLocation()
{
#if 0 // Find location using the obb
	int viewLoc = obbox()->intersectWithPlane(
			_region->viewer->view.matrix,
			_region->viewer->view.minX,
			_region->viewer->view.maxX,
			_region->viewer->view.minY,
			_region->viewer->view.maxY);

	location = false;
	if (viewLoc == 1) {
		location = true;
		in = true;
		return;
	} else if (viewLoc == 0) {
		location = true;
		in = false;
		return;
	}

#else
	location = false;
	if (size()==0) {
		location = true;
		in = false;
		return;
	}
	if (!rpn()) {
		int  prod = 1;	// result of product
		int  bloc;	// body location
		int  i=0;
		while (i<size()) {
			const GBody *body = gexpr(i++);
			if (body == &GBody::tnull) break;

			bloc = vbody(body)->location;
			if (prod<=1 && bloc<=1)
				prod = prod && bloc;
			else
			if (prod==0 || bloc==0)
				prod = 0;
			else
				prod = 2;
			if (prod==0) {
				location = true;
				in = false;
				return;
			}
		}
		while (i<size()) {
			const GBody *body = gexpr(i++);
			assert(body != &GBody::tnull);

			bloc = vbody(body)->location;
			if (prod<=1 && bloc<=1)
				prod = prod && !bloc;
			else
			if (prod==0 || bloc==1)
				prod = 0;
			else
				prod = 2;
			if (prod==0) {
				location = true;
				in = false;
				return;
			}
		}

		// New zone or end...
		if (prod==1) {
			location = true;
			in = true;
		}
	} else { // RPN
		int stack[100], b;
		int  sp = -1;	/* stack pointer */

		for (int i=0; i<size(); i++) {
			const GBody *body = gexpr(i);
			if (body == &GBody::tplus) {
				if (sp<1) {	// Error recovery
					location = true;
					in = false;
					return;
				}
				b = stack[sp--];
				if (stack[sp]<=1 && b<=1)
					stack[sp] = stack[sp] && b;
				else
				if (stack[sp]==0 || b==0)
					stack[sp] = 0;
				else
					stack[sp] = 2;
			} else
			if (body == &GBody::tminus) {
				if (sp<1) {	// Error recovery
					location = true;
					in = false;
					return;
				}
				b = stack[sp--];
				if (stack[sp]<=1 && b<=1)
					stack[sp] = stack[sp] && !b;
				if (stack[sp]==0 || b==1)
					stack[sp] = 0;
				else
					stack[sp] = 2;
			} else
			if (body == &GBody::tunion) {
				if (sp<1) {	// Error recovery
					location = true;
					in = false;
					return;
				}
				b = stack[sp--];
				if (stack[sp]<=1 && b<=1)
					stack[sp] = stack[sp] || b;
				if (stack[sp]==1 || b==1)
					stack[sp] = 1;
				else
					stack[sp] = 2;
			} else
			if (body == &GBody::tuniverse)
				stack[++sp] = 1;
			else
				stack[++sp] = vbody(body)->location;
		}
		if (sp != 0) {	// Error recovery
			location = true;
			in = false;
			return;
		}
		if (stack[0] <= 1) {
			location = true;
			in = stack[0];
		}
	}
#endif
#if _DEBUG>1
//	if (location)
//		cout << "Region: " << region->name()
//		     << "\tlocation=" << location << "\tin=" << in << endl;
#endif
} // updateLocation

/** intersectRay
 * @param x,y,z		location of ray
 * @param dx,dy,dz	normalized direction of ray
 * @param tmin		I/O minimum intersection point
 * @param tmax		maximum to scan
 * @return body that an intersection in found value with the minimum (*tmin)
 *
 * NOTE: Similar logic has the method Geometry::intersectRayUndefinedRegion
 *
 */
CBody* VZone::intersectRay(GeometryEngine* engine,
			   const double  x, const double  y, const double  z,
			   const double dx, const double dy, const double dz,
			   double *tmin, const double tmax) const
{
	while (true) {
		CBody *tbody = NULL;
		double t = tmax;

		for (int i=0; i<size(); i++) {
			const GBody *body = gexpr(i);
			if (body->isOperator()) continue;
			CBody *cb = cbody(engine, body);
			if (cb->intersectRay(x,y,z, dx,dy,dz)) {
				if (InRangeOpen(*tmin, cb->tmin, t)) {
					t = cb->tmin * (1.0+SMALL3D);
					tbody = cb;
				} else
				if (InRangeOpen(*tmin, cb->tmax, t)) {
					t = cb->tmax * (1.0+SMALL3D);
					tbody = cb;
				}
			}
		}
		if (!tbody) {
			*tmin = tmax;
			return NULL;
		} else {
			// move it a bit
			// FIXME
			// It could be improved by getting rid of th SMALL3D
			// and find the next intersection (like the above loop
			// and check if it is insideRay or not
			t *= (1.0+SMALL3D);
			*tmin = t;
			if (!insideRay(engine, x,y,z, dx,dy,dz, t))
				return tbody;
		}
	}
} // intersectRay

/** _inside2D */
bool VZone::_inside2D(GeometryEngine* engine,
		      const double  x, const double  y, const double  z,
		      const double dx, const double dy, const double dz) const
{
	if (size()==0) return false;
	CHECK( cout << endl << "VZone: " << *this
			<< " pos= " << x << " " << y << " " << z
			<< " dir= " << dx << " " << dy << " " << dz << endl);
	if (!rpn()) {
		int  i=0;
		// first plus terms (if any)
		while (i<size()) {
			const GBody *body = gexpr(i++);
			if (body == &GBody::tnull) break;
#ifdef EXPERIMENTAL
			CHECK(int checked=cbody(engine, body)->isChecked());
			CHECK(cout << "   + " << body->name()
				     << " in["<<checked<<"]="
				     << cbody(engine, body)->inside2D(engine,x,y,z,dx,dy,dz) << endl);
			if (!cbody(engine, body)->inside2D(engine,x,y,z,dx,dy,dz)) return false;
#else
			if (!cbody(engine, body)->inside2D(x,y,z,dx,dy,dz)) return false;
#endif
		}
		// then negative terms (if any)
		while (i<size()) {
			const GBody *body = gexpr(i++);
#ifdef EXPERIMENTAL
			CHECK(int checked=cbody(engine, body)->isChecked());
			CHECK(cout << "   - " << body->name()
				   << " in["<<checked<<"]="
				   <<  cbody(engine, body)->inside2D(engine,x,y,z,dx,dy,dz) << endl);
			if (cbody(engine, body)->inside2D(engine,x,y,z,dx,dy,dz)) return false;
#else
			if (cbody(engine, body)->inside2D(x,y,z,dx,dy,dz)) return false;
#endif
		}
		CHECK(cout << "   | prod=1" << endl);
		return true;
	} else { // RPN
		bool stack[100], b;
		int  sp = -1;	/* stack pointer */

		for (int i=0; i<size(); i++) {
			const GBody *body = gexpr(i);
			if (body == &GBody::tplus) {
				if (!sp) return false;
				assert(sp>0);
				b = stack[sp--];
				stack[sp] = stack[sp] && b;
				CHECK(cout << "  +  " << stack[sp] << endl);
			} else
			if (body == &GBody::tminus) {
				if (!sp) return false;
				assert(sp>0);
				b = stack[sp--];
				stack[sp] = stack[sp] && !b;
				CHECK(cout << "  -  " << stack[sp] << endl);
			} else
			if (body == &GBody::tunion) {
				if (!sp) return false;
				assert(sp>0);
				b = stack[sp--];
				stack[sp] = stack[sp] || b;
				CHECK(cout << "  |  " << stack[sp] << endl);
			} else
			if (body == &GBody::tuniverse) {
				stack[++sp] = true;
				CHECK(cout << "  @  " << stack[sp] << endl);
			} else {
				CHECK(int checked=cbody(engine, body)->isChecked());
#ifdef EXPERIMENTAL
				stack[++sp] = cbody(engine, body)->inside2D(engine,x,y,z,dx,dy,dz);
#else
				stack[++sp] = cbody(engine, body)->inside2D(x,y,z,dx,dy,dz);
#endif
				if (sp>=100) return false;
				assert(sp<100);
				CHECK(cout << " --> "<<body->name()
					   << " in["<<checked<<"]="
					   << "=" << stack[sp] << ":\t";
					for (int j=0; j<=sp; j++) cout << stack[j] << ' ';
					cout << endl);
			}
		}
		CHECK(cout << region()->name() << "=" << stack[0] << endl << endl);
		if (sp) return false;
		assert(sp==0);
		return stack[0];
	}
} // _inside2D

/** _inside */
bool VZone::inside(GeometryEngine* engine,
		      const double  x, const double  y, const double  z,
		      const double dx, const double dy, const double dz) const
{
	if (size()==0) return false;
	CHECK( cout << endl << "VZone: " << *this
			<< " pos= " << x << " " << y << " " << z
			<< " dir= " << dx << " " << dy << " " << dz << endl);
	if (!rpn()) {
		int  i=0;
		// first plus terms (if any)
		while (i<size()) {
			const GBody *body = gexpr(i++);
			if (body == &GBody::tnull) break;
			if (!cbody(engine, body)->inside(x,y,z,dx,dy,dz)) return false;
		}
		// then negative terms (if any)
		while (i<size()) {
			const GBody *body = gexpr(i++);
			if (cbody(engine, body)->inside(x,y,z,dx,dy,dz)) return false;
		}
		CHECK(cout << "   | prod=1" << endl);
		return true;
	} else { // RPN
		bool stack[100], b;
		int  sp = -1;	/* stack pointer */
//#if NEW
#if 1
		CBody* bstack[100];	// NULL if already pushed to stack, else needs evaluation of inside

		for (int i=0; i<size(); i++) {
			const GBody *body = gexpr(i);
			// a+b = a && b
			if (body == &GBody::tplus) {
				if (!sp) return false;
				assert(sp>0);
				// check if any of the two is known and false
				if ((bstack[sp]==NULL   and !stack[sp]) or		// b
				    (bstack[sp-1]==NULL and !stack[sp-1])) {		// a
					// (a or b)=false => a+b = false
					sp--;
					stack[sp]  = false;
					bstack[sp] = NULL;
				} else {
					// retrieve b
					if (bstack[sp])
						b = bstack[sp--]->inside(x,y,z,dx,dy,dz);
					else
						b = stack[sp--];
					if (b) {
						// b = true => a+b = a
						if (bstack[sp]) {
							stack[sp]  = bstack[sp]->inside(x,y,z,dx,dy,dz);
							bstack[sp] = NULL;
						} // else already there no need to do anything
					} else {
						// b = false => a+b = false
						stack[sp]  = false;
						bstack[sp] = NULL;
					}
				}
				CHECK(cout << "  +  " << stack[sp] << endl);
			} else
			// a-b = a && !b
			if (body == &GBody::tminus) {
				if (!sp) return false;
				assert(sp>0);
				// check if any of the two is known and false
				if ((bstack[sp]==NULL   and  stack[sp]) or		//!b
				    (bstack[sp-1]==NULL and !stack[sp-1])) {		// a
					// (a or !b)=false => a-b = false
					sp--;
					stack[sp]  = false;
					bstack[sp] = NULL;
				} else {
					// retrieve b
					if (bstack[sp])
						b = bstack[sp--]->inside(x,y,z,dx,dy,dz);
					else
						b = stack[sp--];
					if (!b) {
						// b=false => a-b = a
						if (bstack[sp]) {
							stack[sp]  = bstack[sp]->inside(x,y,z,dx,dy,dz);
							bstack[sp] = NULL;
						} // else already there no need to do anything
					} else {
						// b=true => a-b = false
						stack[sp]  = false;
						bstack[sp] = NULL;
					}
				}
				CHECK(cout << "  -  " << stack[sp] << endl);
			} else
			// a | b
			if (body == &GBody::tunion) {
				if (!sp) return false;
				assert(sp>0);
				// check if any of the two is known and true
				if ((bstack[sp]==NULL   and  stack[sp]) or		// b
				    (bstack[sp-1]==NULL and  stack[sp-1])) {		// a
					// (a or b)=true => a|b = true
					sp--;
					stack[sp]  = true;
					bstack[sp] = NULL;
				} else {
					// retrieve b
					if (bstack[sp])
						b = bstack[sp--]->inside(x,y,z,dx,dy,dz);
					else
						b = stack[sp--];
					if (!b) {
						// b=false => a|b = a
						if (bstack[sp]) {
							stack[sp]  = bstack[sp]->inside(x,y,z,dx,dy,dz);
							bstack[sp] = NULL;
						} // else already there no need to do anything
					} else {
						// b=true => a|b = true
						stack[sp]  = true;
						bstack[sp] = NULL;
					}
				}
				CHECK(cout << "  |  " << stack[sp] << endl);
			} else
			if (body == &GBody::tuniverse) {
				stack[++sp] = true;
				bstack[sp]  = NULL;
				CHECK(cout << "  @  " << stack[sp] << endl);
			} else {
				CHECK(int checked=cbody(engine, body)->isChecked());
				bstack[++sp] = cbody(engine, body);
				if (sp>=100) return false;
				assert(sp<100);
				CHECK(cout << " --> "<<body->name()
					   << " in["<<checked<<"]="
					   << "=" << stack[sp] << ":\t";
					for (int j=0; j<=sp; j++) cout << stack[j] << ' ';
					cout << endl);
			}
		}
#else
		for (int i=0; i<size(); i++) {
			const GBody *body = gexpr(i);
			if (body == &GBody::tplus) {
				if (!sp) return false;
				assert(sp>0);
				b = stack[sp--];
				stack[sp] = stack[sp] && b;
				CHECK(cout << "  +  " << stack[sp] << endl);
			} else
			if (body == &GBody::tminus) {
				if (!sp) return false;
				assert(sp>0);
				b = stack[sp--];
				stack[sp] = stack[sp] && !b;
				CHECK(cout << "  -  " << stack[sp] << endl);
			} else
			if (body == &GBody::tunion) {
				if (!sp) return false;
				assert(sp>0);
				b = stack[sp--];
				stack[sp] = stack[sp] || b;
				CHECK(cout << "  |  " << stack[sp] << endl);
			} else
			if (body == &GBody::tuniverse) {
				stack[++sp] = true;
				CHECK(cout << "  @  " << stack[sp] << endl);
			} else {
				CHECK(int checked=cbody(engine, body)->isChecked());
				stack[++sp] = cbody(engine, body)->inside(x,y,z,dx,dy,dz);
				if (sp>=100) return false;
				assert(sp<100);
				CHECK(cout << " --> "<<body->name()
					   << " in["<<checked<<"]="
					   << "=" << stack[sp] << ":\t";
					for (int j=0; j<=sp; j++) cout << stack[j] << ' ';
					cout << endl);
			}
		}
#endif
		CHECK(cout << region()->name() << "=" << stack[0] << endl << endl);
		if (sp) return false;
		assert(sp==0);
		return stack[0];
	}
} // inside

/** insideRay */
bool VZone::insideRay(GeometryEngine* engine,
		      const double  x, const double  y, const double  z,
		      const double dx, const double dy, const double dz,
		      const double  t) const
{
	if (!rpn()) {
		int  i=0;
		// first plus terms (if any)
		while (i<size()) {
			const GBody *body = gexpr(i++);
			if (body == &GBody::tnull) break;
			if (!cbody(engine, body)->insideRay(x,y,z,dx,dy,dz, t)) return false;
		}
		// then negative terms (if any)
		while (i<size()) {
			const GBody *body = gexpr(i++);
			if (cbody(engine, body)->insideRay(x,y,z,dx,dy,dz, t)) return false;
		}
		return true;
	} else { // RPN
		bool stack[100], b;
		int  sp = -1;   /* stack pointer */
//#if NEW
#if 1
		CBody* bstack[100];	// NULL if already pushed to stack, else needs evaluation of inside

		for (int i=0; i<size(); i++) {
			const GBody *body = gexpr(i);
			// a+b = a && b
			if (body == &GBody::tplus) {
				assert(sp>0);
				if (!sp) return false;
				// check if any of the two is known and false
				if ((bstack[sp]==NULL   and !stack[sp]) or		// b
				    (bstack[sp-1]==NULL and !stack[sp-1])) {		// a
					// (a or b)=false => a+b = false
					sp--;
					stack[sp]  = false;
					bstack[sp] = NULL;
				} else {
					// retrieve b
					if (bstack[sp])
						b = bstack[sp--]->insideRay(x,y,z,dx,dy,dz,t);
					else
						b = stack[sp--];
					if (b) {
						// b = true => a+b = a
						if (bstack[sp]) {
							stack[sp]  = bstack[sp]->insideRay(x,y,z,dx,dy,dz,t);
							bstack[sp] = NULL;
						} // else already there no need to do anything
					} else {
						// b = false => a+b = false
						stack[sp]  = false;
						bstack[sp] = NULL;
					}
				}
			} else
			// a-b = a && !b
			if (body == &GBody::tminus) {
				assert(sp>0);
				if (!sp) return false;
				// check if any of the two is known and false
				if ((bstack[sp]==NULL   and  stack[sp]) or		//!b
				    (bstack[sp-1]==NULL and !stack[sp-1])) {		// a
					// (a or !b)=false => a-b = false
					sp--;
					stack[sp]  = false;
					bstack[sp] = NULL;
				} else {
					// retrieve b
					if (bstack[sp])
						b = bstack[sp--]->insideRay(x,y,z,dx,dy,dz,t);
					else
						b = stack[sp--];
					if (!b) {
						// b=false => a-b = a
						if (bstack[sp]) {
							stack[sp]  = bstack[sp]->insideRay(x,y,z,dx,dy,dz,t);
							bstack[sp] = NULL;
						} // else already there no need to do anything
					} else {
						// b=true => a-b = false
						stack[sp]  = false;
						bstack[sp] = NULL;
					}
				}
			} else
			// a | b
			if (body == &GBody::tunion) {
				assert(sp>0);
				if (!sp) return false;
				// check if any of the two is known and true
				if ((bstack[sp]==NULL   and  stack[sp]) or		// b
				    (bstack[sp-1]==NULL and  stack[sp-1])) {		// a
					// (a or b)=true => a|b = true
					sp--;
					stack[sp]  = true;
					bstack[sp] = NULL;
				} else {
					// retrieve b
					if (bstack[sp])
						b = bstack[sp--]->insideRay(x,y,z,dx,dy,dz,t);
					else
						b = stack[sp--];
					if (!b) {
						// b=false => a|b = a
						if (bstack[sp]) {
							stack[sp]  = bstack[sp]->insideRay(x,y,z,dx,dy,dz,t);
							bstack[sp] = NULL;
						} // else already there no need to do anything
					} else {
						// b=true => a|b = true
						stack[sp]  = true;
						bstack[sp] = NULL;
					}
				}
			} else
			if (body == &GBody::tuniverse) {
				stack[++sp] = true;
				bstack[sp]  = NULL;
			} else {
				bstack[++sp] = cbody(engine, body);
				if (sp>=100) return false;
				assert(sp<100);
			}
		}
#else
		for (int i=0; i<size(); i++) {
			const GBody *body = gexpr(i);
			if (body == &GBody::tplus) {
				//assert(sp>=1);
				if (sp<1) return false;
				b = stack[sp--];
				stack[sp] = stack[sp] && b;
			} else if (body == &GBody::tminus) {
				//assert(sp>=1);
				if (sp<1) return false;
				b = stack[sp--];
				stack[sp] = stack[sp] && !b;
			} else if (body == &GBody::tunion) {
				//assert(sp>=1);
				if (sp<1) return false;
				b = stack[sp--];
				stack[sp] = stack[sp] || b;
			} else if (body == &GBody::tuniverse)
				stack[++sp] = true;
			else {
				stack[++sp] = cbody(engine, body)->insideRay(x,y,z,dx,dy,dz,t);
				if (sp>=100) return false;
				assert(sp<100);
			}
		}
#endif
		if (sp!=0) return false;
		return stack[0];
	}
} // insideRay

/** operator << */
ostream& operator << (ostream& s, const VZone& zone)
{
	s << zone.name() << ":" << zone.id();
	return s;
} /* operator << */
