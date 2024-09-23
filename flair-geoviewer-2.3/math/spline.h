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
 * Date:	6-Feb-2013
 */

#ifndef __SPLINE_H
#define __SPLINE_H

#include "array.h"
#include "vector.h"

/** Abstract class of SplineNodes */
class BaseSplineNode {
private:
	double	_x;
public:
	BaseSplineNode(double a=0.0) : _x(a) {}
virtual ~BaseSplineNode() {}

	/* set/get */
	void	x(double a)		{ _x = a; }
	double	x()		const	{ return _x; }

virtual	int	size()		const	= 0;
virtual	double	y(int i)	const	= 0;
virtual	void	y(int i, double v)	= 0;
}; // BaseSplineNode

/** Spline Node */
class SplineNode : public BaseSplineNode {
private:
	double	_y;

public:
	SplineNode(double a=0.0, double v=0.0) : BaseSplineNode(a), _y(v) {}
virtual ~SplineNode() {}

virtual	int	size()		const	{ return 1; }
virtual	double	y()	const		{ return _y; }
virtual	double	y(int)	const		{ return _y; }
virtual	void	y(int, double v)	{ _y = v; }
}; // SplineNode

/** Vector Node */
class VectorSplineNode : public BaseSplineNode {
private:
	Vector	_y;

public:
	VectorSplineNode(double a, Vector v) : BaseSplineNode(a), _y(v) {}
virtual ~VectorSplineNode() {}

virtual	int	size()		const	{ return 3; }
virtual	double	y(int i)	const	{ return _y[i]; }
virtual	void	y(int i, double v)	{ _y[i] = v; }
}; // VectorSplineNode

/** Cardinal Spline class */
class CardinalSpline {
private:
	double	A;
	double	matrix[4][4];
	Array<BaseSplineNode*>	pts;
public:
	CardinalSpline(const double a=0.5);
	~CardinalSpline() { /* delete pts!!!! */ }

	void	add(BaseSplineNode* n)		{ pts.add(n); }
	BaseSplineNode* point(int i)		{ return pts[i]; }
	int	size()				{ return pts.size(); }

	void	calcMatrix(const double a);
	double	spline(int k, const double t, BaseSplineNode* result);
}; // CardinalSpline

#endif
