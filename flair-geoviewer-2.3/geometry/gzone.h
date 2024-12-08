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

#ifndef __GZONE_H
#define __GZONE_H

#include <assert.h>
#include <iostream>

#include "os.h"
#include "bbox.h"
#include "gbody.h"
#include "vector.h"

class GRegion;
class OBBox;
//class BBox;
class Quad;
class Geometry;

enum ExprType {
	EXPR_PRODUCT    = 0,	// Product expects [Plus-Bodies]... <NullBody> [Minus-Bodies]...
	EXPR_RPN        = 1,	// Reverse polish notation with body, +, -, |, @(Universe)
	EXPR_EXPANDED   = 2,	// Expanded standard expression
	EXPR_NORMAL     = 3,	// Normal mathematical expression with parenthesis
};

/** Zone class. Zone is a region in space defined:
 * 1. STD (normal form) with only with intersection or subtraction
 * 2. RPN form with an Reverse Polish Notion expression
 *
 */
class GZone {
private:
	int	_id;		// id of zone
	ExprType _type;		// expression type
	int	_generation;	// Zone generation
	bool	_bodyref;	// bodies has reference to this
	Array<GBody*> expr;	// expression or rpn
mutable dword	_hash;		// hash cache
mutable OBBox*	cached_obbox;

public:
	GRegion	*region;	// parent region

public:
	GZone(GRegion *reg=NULL, const bool ref=true);
	~GZone()				{ clear(); }

	void	clear();

	// Id
	void	id(const int i)			{ _id = i; }
	int	id()			const	{ return _id; }
const	char*	name() const;

	/** set region mode RPN or STanDard */
	void	rpn(const bool m)		{ _type = (ExprType)(int)m; }
	bool	rpn() const			{ return _type==EXPR_RPN; }
	ExprType type() const			{ return _type; }
	void	type(ExprType t)		{ _type = t; }
	//void	rpn(const bool m)		{ _rpn = m; }
	//bool	rpn() const			{ return _rpn; }

	int	size()	const			{ return expr.count(); }
	bool	resize(const int n)		{ return expr.resize(n); }
//	int	size()	const			{ return _n; }
//	bool	resize(const int n);
	GBody*	operator [] (const int i)	{ return expr[i]; }
const	GBody*	operator [] (const int i) const	{ return expr[i]; }

	void	addReference(GBody* body);
	void	addAllReferences();
	void	removeAllReferences();

	bool	add(const char *token, GBody* body, bool addref=true);
	bool	addPlus(GBody *body);
	bool	addMinus(GBody *body);
	dword	hash() const;

	// Generation
	int	generation()		const	{ return _generation; }
	void	nextGeneration(int newGen) {
			_generation = newGen;
			clearOBB();
		}

	// bounding box
	BBox	bbox2D() const;
	BBox	bbox() const;
	OBBox*	obbox() const {
		if (!cached_obbox) cached_obbox = updateOBB();
		return cached_obbox;
	};

	// Operations
	bool	optimize();

	bool	inside( const double  x, const double  y, const double  z,
			const double dx, const double dy, const double dz) const;
	bool	inside( const Vector& p, const Vector& d) const
			{ return inside(p.x, p.y, p.z, d.x, d.y, d.z); }
	bool	insideThreshold(const Vector &v, const Quad *ignore_a,
				const Quad *ignore_b, const Quad *ignore_c) const;

	bool	intersectBbox(const double  x, const double  y, const double  z,
			  const double dx, const double dy, const double dz,
			  const double tmin, const double tmax ) const
			{ return bbox().intersect(x,y,z,dx,dy,dz,tmin,tmax); }

	// CSG operations
	void	parse(Geometry& geometry, const char *expstr);
	void	exp2rpn();
	void	rpn2exp();
	void	rpnorm();
	bool	depth();

	size_t	memory() const;

protected:
	void	size(const int n)		{ assert(n<=expr.capacity()); expr.forceCount(n); }
	void	clearOBB();
	OBBox*	updateOBB() const;

private:
	int	_nullIndex() const	{ assert(type()==EXPR_PRODUCT); return expr.find(&GBody::tnull); }
	int	_rpnrule(int n);
	void	_subTerms(int n, int* lowLeft, int* lowRight);
	void	_copy(int dst, int src, int length);
static	bool	_isPlus(const GBody* t)		{ return t==&GBody::tplus || t==&GBody::tminus;}
static	int	_priority(const GBody* body, int* ip=NULL);
	bool	_optimizeProducts();
	void	_bboxFromPlanes(BBox& bbox) const;

friend class VZone;
friend class CZone;
}; // GZone

std::ostream& operator << (std::ostream&, const GZone&);

#endif
