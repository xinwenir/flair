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
 * Date:	15-Mar-2012
 */

#include <sstream>

#include <stdlib.h>
#include <string.h>

#include "geo.h"
#include "token.h"
#include "gbody.h"
#include "bmath.h"
#include "obbox.h"
#include "gzone.h"
#include "bstring.h"
#include "gregion.h"
#include "geometry.h"

class Quad;

using namespace std;

#if defined(_DUMP) && _DEBUG>1
#	define CHECK(X)	X
#else
#	define CHECK(X)
#endif

/** FLUKA zone definition */
GZone::GZone(GRegion *reg, const bool ref) :
	_id(-1),
	_type(EXPR_PRODUCT),
	_generation(0),
	_bodyref(ref),
	_hash(0),
	cached_obbox(NULL),
	region(reg)
{
} // GZone

/** @return regions name */
const char* GZone::name() const
{
	return region?region->name() : "zone";
} // name

/** clear */
void GZone::clear()
{
	removeAllReferences();
	// set to null all expression
	expr.clear();
	_type = EXPR_PRODUCT;
	_hash = 0;
	clearOBB();
} // clear

void GZone::clearOBB()
{
	if (cached_obbox) {
		delete cached_obbox;
		cached_obbox = NULL;
	}
} // clearOBB

/** addReference */
void GZone::addReference(GBody* body)
{
	if (!_bodyref) return;
	if (!body->isOperator() && !body->hasZone(this))
		body->addZone(this);
} // addReference

/** addAllReferences */
void GZone::addAllReferences()
{
	if (!_bodyref) return;
	for (int i=0; i<size(); i++) {
		GBody* body = expr[i];
		if (!body->isOperator() && !body->hasZone(this))
			body->addZone(this);
	}
} // addAllReferences

/** removeAllReferences */
void GZone::removeAllReferences()
{
	if (!_bodyref) return;
	// Remove the reference to this zone from all bodies
	for (int i=0; i<size(); i++) {
		GBody* body = expr[i];
		if (!body->isOperator())
			body->delZone(this);
	}
} // removeAllReferences

/** add
 * The expression can either be in Reverse Polish Notation (RPN) mode or
 * STanDard mode with unions of intersection/subtractions
 * @param token		token to be added
 * @param body		body pointer to add
 * @param addref	add reference to bodies zones
 * @return		false in case of error, true otherwise
 */
bool GZone::add(const char *token, GBody *body, bool addref)
{
	assert(strlen(token)<sizeof(char16)-1);

	_hash = 0;
	clearOBB();

	if (body) {
		if (addref) addReference(body);
		expr.add(body);
		//body->nextGeneration();
	} else
	if (strlen(token)>1)
		return false;
	else
	if (token[0] == '$') {
		assert(_type==EXPR_PRODUCT);
		expr.add(&GBody::tnull);
	} else
	if (token[0] == '+') {
		assert(_type!=EXPR_PRODUCT);
		expr.add(&GBody::tplus);
	} else
	if (token[0] == '-') {
		assert(_type!=EXPR_PRODUCT);
		expr.add(&GBody::tminus);
	} else
	if (token[0] == '|') {
		assert(_type!=EXPR_PRODUCT);
		expr.add(&GBody::tunion);
	} else
	if (token[0] == '@') {
		assert(_type==EXPR_RPN);
		expr.add(&GBody::tuniverse);
	} else
	if (token[0] == '(') {
		assert(_type==EXPR_NORMAL);
		expr.add(&GBody::tleft);
	} else
	if (token[0] == ')') {
		assert(_type==EXPR_NORMAL);
		expr.add(&GBody::tright);
	} else
		return false;

	return true;
} // add

/** addPlus */
bool GZone::addPlus(GBody *body)
{
	if (type() != EXPR_PRODUCT) return false;
	expr.insert(_nullIndex(), body);
	addReference(body);
	return true;
} // addPlus

/** addMinus */
bool GZone::addMinus(GBody *body)
{
	if (type() != EXPR_PRODUCT) return false;
	addReference(body);
	expr.add(body);
	return true;
} // addMinus

/** @return hash value of expression
 * WARNING it has to be the same as in Viewer_zone("has",...) python binding
 */
dword GZone::hash() const
{
	if (_hash>0) return _hash;
	_hash = 0;
	for (int i=0; i<size(); i++) {
		GBody const* body = expr[i];
		_hash += (_hash<<5) + i;
		_hash += (_hash<<5) + body->hash();
	}
	return _hash;
} // hash

/** inside */
bool GZone::inside(const double  x, const double  y, const double  z,
		  const double dx, const double dy, const double dz) const
{
	if (size()==0) return false;
	CHECK( cout << endl << "GZone: " << *this
			<< " pos= " << x << " " << y << " " << z
			<< " dir= " << dx << " " << dy << " " << dz << endl);
	if (!rpn()) {
		int  i=0;
		// first plus terms (if any)
		while (i<size()) {
			const GBody* body = expr[i++];
			if (body==&GBody::tnull) break;
			CHECK(cout << "   + " << body->name()
				     << " in=" << body->inside(x,y,z,dx,dy,dz) << endl);
			if (!body->inside(x,y,z,dx,dy,dz)) return false;
		}
		// then negative terms (if any)
		while (i<size()) {
			const GBody* body = expr[i++];
			CHECK(cout << "   - " << body->name()
				   << " in=" << body->inside(x,y,z,dx,dy,dz) << endl);
			if (body->inside(x,y,z,dx,dy,dz)) return false;
		}
		CHECK(cout << "   | prod=1" << endl);
		return true;
	} else { // RPN
		bool stack[100], b;
		int  sp = -1;	/* stack pointer */

		for (int i=0; i<size(); i++) {
			const GBody* body = expr[i];
			if (body == &GBody::tplus) {
				if (sp<1) return false;
				b = stack[sp--];
				stack[sp] = stack[sp] && b;
				CHECK(cout << "  +  " << stack[sp] << endl);
			} else
			if (body == &GBody::tminus) {
				if (sp<1) return false;
				b = stack[sp--];
				stack[sp] = stack[sp] && !b;
				CHECK(cout << "  -  " << stack[sp] << endl);
			} else
			if (body == &GBody::tunion) {
				if (sp<1) return false;
				b = stack[sp--];
				stack[sp] = stack[sp] || b;
				CHECK(cout << "  |  " << stack[sp] << endl);
			} else
			if (body == &GBody::tuniverse) {
				stack[++sp] = true;
				CHECK(cout << "  @  " << stack[sp] << endl);
			} else {
				stack[++sp] = body->inside(x,y,z,dx,dy,dz);
				if (sp>=100) return false;
				CHECK(cout << name() << "."<<body->name()<<"="<<stack[sp] << ":\t";
					for (int j=0; j<=sp; j++) cout << stack[j] << ' ';
					cout << endl);
			}
		}
		CHECK(cout << name() << "=" << stack[0] << endl << endl);
		if (sp != 0) return false;
		return stack[0];
	}
} // inside

/** inside_threshold
 * @param v	location to search
 * @param ignore_a, ignore_b, ignore_c quads to ignore in the computation
 * @return:
 *	true if the point is inside
 *	false if outside
 */
bool GZone::insideThreshold(const Vector &v, const Quad *ignore_a,
			const Quad *ignore_b, const Quad *ignore_c) const
{
	assert(!rpn());

	// Positive terms
	int i = 0;
	for (; i<size(); i++) {
		GBody const* gbody = expr[i];
		if (gbody==&GBody::tnull) break;
		if (!gbody->inside(v.x, v.y, v.z, ignore_a, ignore_b, ignore_c))
			return false;
	}

	i++; // Skip the NULL term

	// Negative terms
	for (; i<size(); i++) {
		GBody const* gbody = expr[i];
		if (gbody==&GBody::tnull) break;
		if (!gbody->outside(v.x, v.y, v.z, ignore_a, ignore_b, ignore_c))
			return false;
	}

	return true;
} // insideThreshold

/** _bboxFromPlanes */
void GZone::_bboxFromPlanes(BBox& bb) const
{
	// Produce every possible combination of planes
	for (int i=0; i<size(); i++) {
		GBody const* abody = expr[i];
		if (abody==&GBody::tnull) continue;
		for (int ii=0; ii<abody->nQ(); ii++) {
			const Quad& pa = abody->Q(ii);

			int j, jj;
			if (ii<abody->nQ()-1) {
				j  = i;
				jj = ii+1;
			} else {
				j  = i+1;
				jj = 0;
			}

			for (; j<size(); j++) {
				GBody const* bbody = expr[j];
				if (bbody==&GBody::tnull) continue;
				for (; jj<bbody->nQ(); jj++) {
					const Quad& pb = bbody->Q(jj);

					int k, kk;
					if (jj<bbody->nQ()-1) {
						k  = j;
						kk = jj+1;
					} else {
						k  = j+1;
						kk = 0;
					}

					for (; k<size(); k++) {
						GBody const* cbody = expr[k];
						if (cbody==&GBody::tnull) continue;
						for (; kk<cbody->nQ(); kk++) {
							const Quad& pc = cbody->Q(kk);
							// Solve system and find vertices, add them to the bb
							Matrix3 m;
							m(0,0) = pa.Cx; m(0,1) = pa.Cy; m(0,2) = pa.Cz;
							m(1,0) = pb.Cx; m(1,1) = pb.Cy; m(1,2) = pb.Cz;
							m(2,0) = pc.Cx; m(2,1) = pc.Cy; m(2,2) = pc.Cz;
							if (m.inverse()) {
								Vector v = m.mult(-pa.C, -pb.C, -pc.C);
								if (insideThreshold(v, &pa, &pb, &pc))
									bb.add(v.x, v.y, v.z);
							}

						}
						kk = 0;
					}
				}
				jj = 0;
			}
		}
	}
} // _bboxFromPlanes

/** guess of 3D bounding box
 * @return bounding box of zone
 */
BBox GZone::bbox() const
{
	BBox bb;
	if (!rpn()) {
#if 1
		// Check if expression contains only planes!
		// FIXME should accept RPP, BOX, WED, RAW, ARB also
		bool onlyPlanes = true;
		for (int i=0; i<size(); i++) {
			GBody const* gbody = expr[i];
			if (gbody==&GBody::tnull) continue;
			//if (gbody->type()>RAWbody && gbody->type()!=ARBbody) {
			if (gbody->type()>PLAbody && gbody->type()!=ARBbody) {
				onlyPlanes = false;
				break;
			}
		}
		if (onlyPlanes) {
			_bboxFromPlanes(bb);
			return bb;
		}
#endif

		// Mixed expression
		bb.infinite();

		// FIRST PASS
		// Compute a rough bounding box using the bounding boxes of the bodies
		bool negative = false;
		for (int i=0; i<size(); i++) {
			GBody const* gbody = expr[i];

			// Detect negative term
			if (gbody==&GBody::tnull) {
				negative = true;
				continue;
			}

			if (!negative)
				bb += gbody->bbox();
			else {
				// Only subtract bodies for which
				// inside bbox = outside bbox
				// FIXME Difference should ONLY work if there is no
				//       arbitrary rotation (different than 0,90,180,270deg)
				if (!gbody->hasMatrix() &&
						(gbody->type() <= XYPbody || gbody->type() == RPPbody))
					bb -= gbody->bbox();
			}
			if (!bb.isValid()) return bb;
		}

		// SECOND PASS
		// Find all planes, intersect them with the bounding box and
		// create a new bounding box containing these vertices +
		// bounding box vertices if they are inside the zone

		bool changed;
		do {
			// Find all planes in bodies
			negative = false;
			changed = false;
			double vol = bb.volume();
			for (int i=0; i<size(); i++) {
				GBody const* gbody = expr[i];
				// Detect negative term
				if (gbody==&GBody::tnull) {
					negative = true;
					continue;
				}
				if (gbody->type() <= PLAbody)
					if (bb.intersectPlane(gbody->Q(0), negative)) {
						if (!Eq(vol, bb.volume(), vol*SMALL3D))
							changed = true;
					}
			}
		} while (changed);
	} else {
		assert(0);
	}
	return bb;
} // bbox

/** guess of 3D oriented bounding box
 * updates _obbox
 */
OBBox* GZone::updateOBB() const
{
	OBBox* _obbox = new OBBox();
	if (!rpn()) {
		int pass;
		// Iterate body number times so each body has a change to
		// improve the BB
		for (pass=0; pass<size()-1; pass++) {
			bool negative = false;
			//double vol = _obbox.volume();
			for (int i=0; i<size(); i++) {
				GBody const* gbody = expr[i];

				// Detect negative term
				if (gbody==&GBody::tnull) {
					if (!_obbox->isValid())
						return _obbox;
					negative = true;
					continue;
				}

				if (!negative) {
					if (_obbox->isValid())
						_obbox->Intersect(*gbody->obbox(false));
					else
						*_obbox = *gbody->obbox(false);
				} else {
					// Subtract the inner bounding box
					_obbox->Difference(*gbody->obbox(true));
				}
				if (!_obbox->isValid())
					return _obbox;
			}
			// If the volume didn't change from last time, assume it
			// is stable
			//if (Eq0(_obbox.volume() - vol, SMALL2))
			//break;
		}
		// All bodies had a chance to improve the BB
		// if (pass == size()-1) cout << "Exhausted bb stabilization passes, must be optimal" << endl;
	} else { // RPN
#if 0
		// FIXME: review!
		OBBox aux;
		const int stack_size = 100;
		OBBox *stack[stack_size], *b;

		int  sp = -1;   /* stack pointer */

		for (int i=0; i<size(); i++) {
			const GBody *body = expr[i++];
			if (body == GZone::tplus) {
				if (sp<1) break;
				b = stack[sp];
				sp--;
				aux = OBBox(*stack[sp]);
				stack[sp] = &aux;
				stack[sp]->Intersect(*b);
			} else if (body == GZone::tminus) {
				if (sp<1) break;
				b = stack[sp];
				sp--;
				aux = OBBox(*stack[sp]);
				stack[sp] = &aux;
				stack[sp]->Difference(*b);
			} else if (body == GZone::tunion) {
				if (sp<1) break;
				b = stack[sp];
				sp--;
				aux = OBBox(*stack[sp]);
				stack[sp] = &aux;
				stack[sp]->Union(*b);
			} else if (body == GZone::tuniverse) {
				++sp;
				aux.infinite();
				stack[sp] = &aux;
			} else {
				++sp;
				stack[sp] = body->obbox(false);
				if (sp>=100) break;
				assert(sp<100);
			}
		}
		if (sp==0)
			*_obbox = *stack[0];
#endif
	}
	return _obbox;
} // updateOBB

//#define DUMP
/** optimize the zone's expression only in Product format
 * @return true if any optimization took place
 */
bool GZone::optimize()
{
	if (type() != EXPR_PRODUCT) return false;

	// Remove all references
	removeAllReferences();

	// 1st find negative terms starting index
	int negative = _nullIndex()+1;

#ifdef _DUMP
	cout << "Zone::optimize\nIn= ";
	for (int i=0; i<size(); i++)
		if (expr[i]) cout << (i<negative?" +":" -") << expr[i]->name() << "[" << expr[i]->typeStr() << "]";
	cout << endl;
#endif

	// 2nd remove all duplicate terms:
	//	+A +A -> +A
	//	+A -A ->  0
	//	-A -A -> -A
	for (int i=0; i<size(); i++) {
		GBody const* body = expr[i];
		if (body==&GBody::tnull) continue;
		for (int j=i+1; j<size(); j++) {
			if (expr[j]==NULL) continue;
			if (body != expr[j]) continue;
			if (i<negative && j<negative)
				expr[j] = NULL;			// +A +A
			else
			if (i>=negative && j>=negative)
				expr[j] = NULL;			// -A -A
			else {
				size(0);			// +A -A
				return true;
			}
		}
	}

	// TODO XXX Find bounding box for better location estimation
	//BBox bb = bbox();

	// Geometrical optimization
	for (int i=0; i<size(); i++) {
		GBody const* A = expr[i];
		if (A==NULL || A->isOperator()) continue;
		bool plusA = i<negative;
		for (int j=i+1; j<size(); j++) {
			GBody const* B = expr[j];
			if (B==NULL || B->isOperator()) continue;
			bool plusB = j<negative;

			// Find location of A vs B
			Location loc = A->bbLocationWrt(B);	// Check first bounding boxes
			if (loc != LOCATION_OUTSIDE)		// if overlap then check better
				loc = A->locationWrt(B);
#ifdef _DUMP
			cout <<"A="<<A->name()<<" B="<<B->name() << " loc=" << loc << endl;
			cout << "\tAbb="<<A->bbox() << endl;
			cout << "\tBbb="<<B->bbox() << endl;
#endif
			switch (loc) {
				case LOCATION_OUTSIDE:
					// +A +B = 0
					// +A -B = +A
					// -A +B = +B
					// -A -B = same
					if (plusA) {
						if (plusB) {	// +A +B -> 0
							size(0);
							DUMP(cout << "Out= -none-" << endl);
							return true;
						} else		// +A -B -> +A
							expr[j] = NULL;
					} else
						if (plusB)	// -A +B -> +B
							expr[i] = NULL;
					break;

				case LOCATION_AinB:
					// +A +B = +A
					// +A -B = 0
					// -A +B = same
					// -A -B = -B
					if (plusA) {
						if (plusB)	// +A +B -> +A
							expr[j] = NULL;
						else {		// +A -B -> 0
							size(0);
							DUMP(cout << "Out= -none-" << endl);
							return true;
						}
					} else
						if (!plusB)	// -A -B -> -B
							expr[i] = NULL;
					break;

				case LOCATION_BinA:
					// +A +B = +B
					// +A -B = same
					// -A +B = 0
					// -A -B = -A
					if (plusA) {
						if (plusB)	// +A +B -> +B
							expr[i] = NULL;
					} else {
						if (plusB) {	// -A +B -> 0
							size(0);
							DUMP(cout << "Out= -none-" << endl);
							return true;
						} else		// -A -B -> -A
							expr[j] = NULL;
					}
					break;

				default: /* do nothing */;
			}
			// If A==NULL go to next
			if (expr[i]==NULL) break;
		}
	}

#ifdef _DUMP
	cout << "Mid=";
	for (int i=0; i<size(); i++)
		if (expr[i]) cout << (i<negative?" +":" -") << expr[i]->name() << "[" << expr[i]->typeStr() << "]";
	cout << endl;
#endif

	// Finally: fix by removing nulls
	// WARNING reference to bodies
	bool opt = false;
	int pos = 0;
	for (int i=0; i<negative-1; i++) {
		GBody* body = expr[i];
		if (body)
			expr[pos++] = body;
		else
			opt = true;
	}
	expr[pos++] = &GBody::tnull;
	for (int i=negative; i<size(); i++) {
		GBody* body = expr[i];
		if (body)
			expr[pos++] = body;
		else
			opt = true;
	}

	size(pos);		// Fix size
	addAllReferences();	// restore references

#ifdef _DUMP
	cout << "Out=";
	for (int i=0; i<size(); i++)
		if (expr[i]) cout << (i<negative?" +":" -") << expr[i]->name() << "[" << expr[i]->typeStr() << "]";
	cout << endl;
#endif
	return opt;
} // optimize

/** _optimizeProducts */
bool GZone::_optimizeProducts()
{
	assert(type()==EXPR_RPN);

	// Product target should be in the form
	// ... t1 t2 [+-] [ti [+-]]* ...
	//     i               iend
	int i = 0;
	bool optimized = false;
	while (i<expr.count()-2) {
		if (expr[i]->isOperator()) {
			i++;
			continue;
		}
		// Find as many repetition of ti [+-]
		int iend = i + 1;
		while (iend<expr.count()-1) {
			if (!(expr[iend]->isOperator()) && _isPlus(expr[iend+1]))
				iend += 2;
			else
				break;
		}
		if (iend != i+1) {
			// Product found
			GZone product(NULL,false);
			product.type(EXPR_PRODUCT);
			// add first + term (Universe? @)
			product.add(expr[i]->name(), expr[i]);
			for (int j=i+1; j<iend-1; j++)
				if (expr[j+1]==&GBody::tplus)
					product.add(expr[j]->name(), expr[j]);
			product.add("$",NULL);
			for (int j=i+1; j<iend-1; j++)
				if (expr[j+1]==&GBody::tminus)
					product.add(expr[j]->name(), expr[j]);

			cout << "Product=" << product << endl;

			if (product.optimize()) {
				// Replace expression with optimized one!
				cout << "Optimized=" << product << endl;
				bool plus=true;
				bool first=true;
				int ii = i;
				for (int j=0; j<product.size(); j++) {
					GBody *body = product[j];
					if (body==&GBody::tnull)
						plus = false;
					else
					if (first) {
						expr[ii++] = body;
						first = false;
					} else {
						expr[ii++] = body;
						expr[ii++] = (plus? &GBody::tplus : &GBody::tminus);
					}
				}
				// delete terms from ii-iend
				optimized = true;
				expr.erase(ii,iend);
				cout << "Expr=" << *this << endl;
			}
		}
		i = iend;
	}
	return optimized;
} // _optimizeProducts

/** tokenize a string to expr */
void GZone::parse(Geometry& geometry, const char *expstr)
{
	stringstream stream(expstr);
	Token token(stream);
	token.next();
	TokenType tok;

	type(EXPR_NORMAL);
	do {
		tok = token();
		if (tok != eof_tok)
			add(token.string(), geometry.getBody(token.string()));
		token.next();
	} while (tok);
} // parse

/** return token input and output
 * @param ip input priority
 * @return output priority
 */
//               Operator In Out
//priorities = {"(" : [ 99,  0],
//		"|" : [  1,  2],
//		"+" : [  3,  4],
//		"-" : [  3,  4],
//		")" : [  0, 99],
//		" " : [101,100]}
int GZone::_priority(const GBody* body, int* ip)
{
	// Find priorities
	if (body==&GBody::tleft) {
		if (ip) *ip = 99;
		return 0;
	} else
	if (body==&GBody::tunion) {
		if (ip) *ip =  1;
		return 2;
	} else
	if (body==&GBody::tplus || body==&GBody::tminus) {
		if (ip) *ip =  3;
		return 4;
	} else
	if (body==&GBody::tright) {
		if (ip) *ip =  0;
		return 99;
	} else {
		if (ip) *ip = 101;
		return 100;
	}
} // _priority

/** @return depth of the rpn expression */
bool GZone::depth()
{
	int d = 0;
	for (int i=0; i<size(); i++) {
		if (expr[i]->isOperator())
			d--;
		else
			d++;
	}
	return d==1;
} // depth

/* ---------------------------------------------------------------------- *
 * exp2rpn
 *
 * Convert the FLUKA Boolean expression to Reverse Polish Notation (RPN)
 * Since the normalization routine does not accept the +/- as signs to
 * objects, the routine is converting the leading - to @- (Universe minus)
 * where the special symbol @ is treated as the universe.
 * Furthermore the leading + is ignored as well as the leading | which is
 * accepted by fluka.
 *
 * WARNING: No check is done for the correctness of the expression apart
 *          from the parenthesis nesting
 *
 * ie.
 *        A+B         -> A B +
 *       (A+B)|C      -> A B + C |
 *       (A|B)|C+D|E  -> A B | C D + | E |
 *       -A           -> @ A -
 *       -(-A)        -> @ @ A - -
 *
 * The routine is using the same array for returning the Reverse Polish
 * expression, since the format is more compact.
 * This is generally true apart one case -A -> @ A -
 *
 * Priorities are treated as:
 *     Operator  Priority         In     Out
 *     --------  --------        ---     ---
 *       |       lower             1       2
 *       +       high              3       4
 *       -       high              3       4
 *       (       higher           99       0
 *       )       higher            0      99
 *       object  highest         101     100
 *
 * Algorithm
 *  Consider the expression as a train moving on a railroad with a
 *  T-shape, where each token is one wagon
 *
 *                      <-  (A|B)|C+D
 *      ------------.   .------------
 *      RPN-End      \ /      Exp-End
 *                    |
 *                   S|
 *                   t|
 *                   a|
 *                   c|
 *                   k|
 *
 *  Each wagon to move from the Exp-End to the RPN-End it has to make
 *  a stop first in the Stack-End. Before entering in the stack, the
 *  priority-IN will be checked against the objects in the stack.
 *  All top-most objects currently present in the stack with
 *  higher priority-OUT will be transfered from the stack to the
 *  RPN-End. Apart from the opening parenthesis ( which is discarded.
 *
 *  Example:
 *  (1)                         A+B|C (1)
 *  (2) A                        +B|C (2)
 *  (3) A                         B|C (3)
 *  (4) A                          |C (4)
 *  (5) A B +                       C (5)
 *  (6) A B + C
 *  (7) A B + C |
 *      ------------.   .------------
 *      RPN-End      \ /      Exp-End
 *                    |
 *                   S|
 *                   t|
 *                   a|
 *                   c|     B   C
 *                   k| A + + | | |
 *                      1 2 3 4 5 6 7
 */
void GZone::exp2rpn()
{
	bool newproduct = true;
	int i = 0;
	int m = 0;
	Array<GBody*> stack;

	assert(type()==EXPR_NORMAL);

	while (i < size()) {
		GBody* body = expr[i];

		//# Check for special leading chars
		if (newproduct && (body==&GBody::tunion || body==&GBody::tplus)) {
			newproduct = (body==&GBody::tunion);
			i++;
			continue;
		}

		if (newproduct && body==&GBody::tminus) {
			expr.insert(i, &GBody::tnull);	// insert space in ith position
			body = &GBody::tuniverse;
		}

		newproduct = (body==&GBody::tleft || body==&GBody::tunion);

		int ip;
		_priority(body,&ip);

		// Remove from the stack everything with higher priority
		while (!stack.empty() && ip<_priority(stack.tail())) {
			if (stack.tail() != &GBody::tleft) {
				expr[m] = stack.tail();
				m++;
			}
			stack.pop();
		}

		//# Push it into the stack
		if (body != &GBody::tright)
			stack.append(body);
		else {
			// Should be an opening parenthesis
			if (stack.empty()) {
				throw 0;
				//raise CSGException("Unbalanced parenthesis")
			}
			stack.pop();
		}
		i++;
	}

	// Empty Stack
	while (!stack.empty()) {
		if (stack.tail() == &GBody::tleft) {
			throw 0;
			//CSGException("Unbalanced parenthesis")
		}

		if (m >= expr.size())
			expr.append(stack.tail());
		else
			expr[m] = stack.tail();
		stack.pop();
		m++;
	}

	// Delete unwanted items
	expr.erase(m, expr.size());
	//for i in xrange(len(expr)-1,m-1,-1):
	//	del expr[i]

	type(EXPR_RPN);
} // exp2rpn

/** ---------------------------------------------------------------------- *
 * rpn2exp
 *
 * Convert a NORMALIZED Reverse Polish notation to a standard expression
 * WARNING: The routine expects an expression where for each UNION
 * operator the right-sub-expression is a product while the left can be
 * UNION or a product
 * ---------------------------------------------------------------------- */
void GZone::rpn2exp()
{
	assert(type() == EXPR_RPN);

	Array<GBody*> newexpr(expr.size());
	Array<GBody*> plus;
	Array<GBody*> minus;
//	plus.compare(GBody::compare);
//	minus.compare(GBody::compare);

	int i = 0;
	int nstack = 0;
	bool endprod = false;
	while (i<expr.size()) {
		GBody* body = expr[i];
		bool lastPlus = true;

		if (body->isOperator())
			nstack--;
		else {
			// First term is always a plus
			// .. peek then next operator to check for sign
			if (plus.empty() || expr[i+1] == &GBody::tplus) {
				plus.append(expr[i]);
				lastPlus = true;
			} else
			if (expr[i+1] == &GBody::tminus) {
				minus.append(expr[i]);
				lastPlus = false;
			}
			nstack++;
		}

		if (nstack==0)
			endprod = true;
		else
		if (nstack==3) {
			if (lastPlus)
				plus.pop();
			else
				minus.pop();
			i -= 2;
			endprod = true;
		} else
		if (body==&GBody::tunion) {
			i -= 2;
			endprod = true;
		} else
		if (i==expr.top())
			endprod = true;

		if (endprod) {
			//optZone(plus,minus)
//			if (!plus.empty() || !minus.empty())
//				zones.append((plus,minus));
			bool something_added = false;
			for (int ii=0; ii<plus.count(); ii++)
				if (plus[ii]->name()[0] != ' ') {
					newexpr.append(&GBody::tplus);
					newexpr.append(plus[ii]);
					something_added = true;
				}
			for (int ii=0; ii<minus.count(); ii++)
				if (minus[ii]->name()[0] != ' ') {
					newexpr.append(&GBody::tminus);
					newexpr.append(minus[ii]);
					something_added = true;
				}
			if (something_added)
				newexpr.append(&GBody::tunion);
			plus.clear();
			minus.clear();
			nstack  = 0;
			endprod = false;
		}
		i++;
	}

	// Remove last |
	if (!newexpr.empty()) newexpr.pop();

	// Copy back expression
	expr.clear();
	expr.allocate(newexpr.count());
	memcpy(&expr[0], &newexpr[0], sizeof(GBody*)*newexpr.count());

	//# Remove duplicates of products
//	rmDoubles(zones);

#if 0
	// Reconstruct expression
	expr = []
	for plus,minus in zones:
		if len(expr)>0 and expr[-1]!="|": expr.append("|")
		//# Fill the new array
		for j in plus:
			expr.append("+")
			expr.append(j)
		for j in minus:
			expr.append("-")
			expr.append(j)
#endif
	type(EXPR_EXPANDED);
//	return expr
} // rpn2exp

/** --------------------------------------------------------------------- *
 * rpnorm
 *
 * Normalize a CG expression given in Reverse Polish Notation.
 * Normalized CG expression is an expression given as sum (Boolean OR) of
 * products (Boolean intersection or subtraction).
 * The normalization (expansion of parenthesis and operator priorities)
 * should be performed by recursively calling the RPNRULE subroutine.
 * Since Fortran-77 doesn't have recursion, call the RPNRULE for every
 * operator starting from the right-most one, until no rule is found.
 *
 * WARNING: fix references
 * ----------------------------------------------------------------------
 */
void GZone::rpnorm()
{
	// Loop until there is no any extra change needed
	// Scan to find the first operators
	bool changed;
	_optimizeProducts();
	do {
		changed = false;
		int i = expr.count()-1;
		while (i >= 4) {
			if (expr[i]->isOperator()) {
				int l = expr.count();
				int rule = _rpnrule(i);
				if (rule>0) {
					changed = true;
					i += expr.count()-l+1;
				}
			}
			i--;
		}

		// Find possible products and optimize them
		if (changed) _optimizeProducts();
	} while (changed);
} // rpnorm

/** --------------------------------------------------------------------- *
 * rpnrule
 *
 * Find a matching rule and apply it on the sub-expression starting from
 * the N position in the Reverse Polish Notation
 *
 * An expression is in normal form when all the parenthesis are expanded
 * and the expression is described as a sum (UNIONS) of products
 * (INTERSECTIONS and/or SUBTRACTIONS)
 *
 * An expression can be converted to normal form by repeatedly applying
 * the following set of production rules to the expression and then to its
 * sub-expressions:
 *
 *    Normal Form                        Reverse Polish Notation
 *  1. X-(Y|Z) -> (X-Y)-Z                X Y Z | -  ->  X Y - Z -
 *  2. X+(Y|Z) -> (X+Y)|(X+Z)            X Y Z | +  ->  X Y + X Z + |
 *  3. X-(Y+Z) -> (X-Y)|(X-Z)            X Y Z + -  ->  X Y - X Z - |
 *  4. X+(Y+Z) -> (X+Y)+Z                X Y Z + +  ->  X Y + Z +
 *  5. X-(Y-Z) -> (X-Y)|(X+Z)            X Y Z - -  ->  X Y - X Z + |
 *  6. X+(Y-Z) -> (X+Y)-Z                X Y Z - +  ->  X Y + Z -
 *  7. X|(Y|Z) -> (X|Y)|Z                X Y Z | |  ->  X Y | Z |
 *  8. (X-Y)+Z -> (X+Z)-Y                X Y - Z +  ->  X Z + Y -
 *  9. (X|Y)-Z -> (X-Z)|(Y-Z)            X Y | Z -  ->  X Z - Y Z - |
 * 10. (X|Y)+Z -> (X+Z)|(Y+Z)            X Y | Z +  ->  X Z + Y Z + |
 * X,Y, and Z here match both primitives or sub-expressions.
 * ---------------------------------------------------------------------- */
int GZone::_rpnrule(int n)
{
	// Reset rule
	int rule = 0;

	// Top-most operator
	GBody* op = expr[n];
	if (!op->isOperator()) return rule;

	// Right operator
	GBody* rop = expr[n-1];

	// Find left and right sub-trees
	int ll, lr;
	_subTerms(n, &ll, &lr);

	// Left operator
	GBody* lop = &GBody::tnull;
	if (lr>0) lop = expr[lr-1];

	// Find Rule
	if (op==&GBody::tminus && rop==&GBody::tunion)
		rule = 1;
	else
	if (op==&GBody::tplus  && rop==&GBody::tunion)
		rule = 2;
	else
	if (op==&GBody::tminus && rop==&GBody::tplus)
		rule = 3;
	else
	if (op==&GBody::tplus  && rop==&GBody::tplus)
		rule = 4;
	else
	if (op==&GBody::tminus && rop==&GBody::tminus)
		rule = 5;
	else
	if (op==&GBody::tplus  && rop==&GBody::tminus)
		rule = 6;
	else
	if (op==&GBody::tunion && rop==&GBody::tunion)
		rule = 7;
	else
	if (op==&GBody::tplus  && lop==&GBody::tminus)
		rule = 8;
	else
	if (op==&GBody::tminus && lop==&GBody::tunion)
		rule = 9;
	else
	if (op==&GBody::tplus  && lop==&GBody::tunion)
		rule = 10;
	else
		return rule;

	// limits of subexpressions X,Y,Z
	int Xu, Xl, /*Yu,*/ Yl, Zu, Zl;
	int L;		// aux variables

	// Find sub expressions X Y Z
	if (rule<=7) {    // X op (Y rop Z)
		Xu = lr-1;
		Xl = ll;
		_subTerms(n-1, &ll, &lr);
		//Yu = lr-1;
		Yl = ll;
		Zu = n-2;
		Zl = lr;
	} else {	// (X lop Y) op Z
		Zu = n-1;
		Zl = lr;
		L  = lr-1;
		_subTerms(L, &ll, &lr);
		Xu = lr-1;
		Xl = ll;
		//Yu = L-1;
		Yl = lr;
	}

//	printf("Rule: %d\n",rule);
//	printf(" X[%d,%d,] = \n",Xl,Xu); //,"]=",expr[Xl:Xu+1]
//	printf(" Y[%d,%d,] = \n",Yl,Yu); //,"]=",expr[Yl:Yu+1]
//	printf(" Z[%d,%d,] = \n",Zl,Zu); //,"]=",expr[Zl:Zu+1]

	// Expand the rule
	switch (rule) {
		//# 1. X-(Y|Z) -> (X-Y)-Z	 X Y Z | -  ->  X Y - Z -
		case 1:
			// Leave X Y
			// Insert a - operator after Y
			expr.insert(Zl, &GBody::tminus);
			// Chop length by 1
			expr.erase(Zu+2);
			// Change the last operator to -
			expr[Zu+2] = &GBody::tminus;
			break;

		// 2. X+(Y|Z) -> (X+Y)|(X+Z)     X Y Z + |  ->  X Y + X Z + |
		case 2:
			// Leave X Y
			// Insert a + operator after Y
			expr.insert(Zl, &GBody::tplus);
			// Copy X after the + operator
			_copy(Zl+1, Xl, Xu-Xl+1);
			Zu += Xu-Xl+2;
			// Change last 2 operators to + |
			expr[Zu+1] = &GBody::tplus;
			expr[Zu+2] = &GBody::tunion;
			break;

		// 3. X-(Y+Z) -> (X-Y)|(X-Z)     X Y Z + -  ->  X Y - X Z - |
		case 3:
			// Leave X Y
			// Insert a - operator after Y
			expr.insert(Zl, &GBody::tminus);
			// Copy X after the - operator
			_copy(Zl+1, Xl, Xu-Xl+1);
			Zu += Xu-Xl+2;
			// Change last 2 operators to - |
			expr[Zu+1] = &GBody::tminus;
			expr[Zu+2] = &GBody::tunion;
			break;

		// 4. X+(Y+Z) -> (X+Y)+Z	 X Y Z + +  ->  X Y + Z +
		case 4:
			// Leave X Y
			// Insert a + operator after Y
			expr.insert(Zl, &GBody::tplus);
			// Chop length by 1
			expr.erase(Zu+2);
			// Change the last operator to +
			expr[Zu+2] = &GBody::tplus;
			break;

		// 5. X-(Y-Z) -> (X-Y)|(X+Z)     X Y Z - -  ->  X Y - X Z + |
		case 5:
			// Leave X Y
			// Insert a - operator after Y
			expr.insert(Zl, &GBody::tminus);
			// Copy X after the - operator
			_copy(Zl+1, Xl, Xu-Xl+1);
			Zu += Xu-Xl+2;
			// Change last 2 operators to + |
			expr[Zu+1] = &GBody::tplus;
			expr[Zu+2] = &GBody::tunion;
			break;

		// 6. X+(Y-Z) -> (X+Y)-Z	 X Y Z - +  ->  X Y + Z -
		case 6:
			// Leave X Y
			// Insert a + operator after Y
			expr.insert(Zl, &GBody::tplus);
			// Chop length by 1
			expr.erase(Zu+2);
			// Change the last operator to -
			expr[Zu+2] = &GBody::tminus;
			break;

		// 7. X|(Y|Z) -> (X|Y)|Z	 X Y Z | |  ->  X Y | Z |
		case 7:
			// Leave X Y
			// Insert a | operator after Y
			expr.insert(Zl, &GBody::tunion);
			// Chop length by 1
			expr.erase(Zu+2);
			// Change the last operator to |
			expr[Zu+2] = &GBody::tunion;
			break;

		// 8. (X-Y)+Z -> (X+Z)-Y	 X Y - Z +  ->  X Z + Y -
		case 8:
			// Leave X
			// Copy "Z +" after X
			L = Zu-Zl+2;
			_copy(Yl, Zl, Zu-Zl+2);
			// Delete old "Z +"
			expr.erase(Zl+L,Zl+L+L);	// XXX to be verified!
			break;

		// 9. (X|Y)-Z -> (X-Z)|(Y-Z)     X Y | Z -  ->  X Z - Y Z - |
		//10. (X|Y)+Z -> (X+Z)|(Y+Z)     X Y | Z +  ->  X Z + Y Z + |
		case 9:
		case 10:
			// Leave X
			// Copy "Z -" or "Z +" after X
			L = Zu-Zl+2;
			// Correct Z position
			_copy(Yl, Zl, L);
			Zl += L;
			Zu += L;
			// Delete the | in front of Z
			expr.erase(Zl-1);
			// Add | at the end
			expr.insert(Zu+1, &GBody::tunion);
			break;
	}

	return rule;
} // _rpnrule

/** --------------------------------------------------------------------- *
 * subTerms
 *
 * This routine returns the pointers in the RPN terms array of the
 * starting point of the left sub-expression LOWLEFT and right
 * sub-expression LOWRIGHT given the pointer of the expr-operator NTX.
 * The searching is performed by scanning from right to left the number
 * of operators and objects pushed into the stack
 *
 * EXP          (expr-left) op (expr-right)
 *
 * RPN ...... | expr-left | expr-right | op | .......
 *           Lowleft      LowRight  op-1 ntx
 * ---------------------------------------------------------------------- */
void GZone::_subTerms(int n, int* lowLeft, int* lowRight)
{
	int nop = 0;
	*lowRight = 0;
	*lowLeft  = 0;

	while (n>=0) {
		if (expr[n]->isOperator())
			nop++;
		else
			nop--;

		n--;
		if (nop==0) {
			if (*lowRight==0) {
				*lowRight = n+1;
				nop++;
				continue;
			} else {
				*lowLeft = n+1;
				return;
			}
		}
	}
} // _subTerms

/** copy a subexpression
 *  @param dst		insert position
 *  @param src		starting position
 *  @param length	length of subexpression
 */
void GZone::_copy(int dst, int src, int length)
{
	// Create space by shifting everything by length
	expr.shiftby(dst, length);

	// Correct starting position if needed
	int nf;
	if (src > dst)
		nf = src + length;
	else
		nf = src;

	// copy terms
	memcpy(&expr[dst], &expr[nf], sizeof(GBody*)*length);
} // _copy

/** memory */
size_t GZone::memory() const
{
	return sizeof(GZone); // + _maxSize * sizeof(GBody*);
} // memory

/** operator << */
ostream& operator << (ostream& s, const GZone& zone)
{
	if (zone.region!=NULL)
		s << zone.name() << "-" << zone.id() << ":";

	bool plus = true;
	for (int i=0; i<zone.size(); i++) {
		GBody const* body = zone[i];
		switch (zone.type()) {
			case EXPR_PRODUCT:
				if (body==&GBody::tnull)
					plus = false;
				else {
					s << (plus?" +" :" -");
					s << body->name();
				}
				break;

			case EXPR_RPN:
				s << ' ' << body->name();
				break;

			case EXPR_EXPANDED:
			case EXPR_NORMAL:
				if (body==&GBody::tplus || body==&GBody::tminus)
					s << ' ' << body->name();
				else
				if (body==&GBody::tunion)
					s << ' ' << body->name();
				else
				if (body==&GBody::tleft)
					s << body->name() << ' ';
				else
				if (body==&GBody::tright)
					s << ' ' << body->name();
				else
					s << body->name();
				break;
		}
	}
	return s;
} /* operator << */
