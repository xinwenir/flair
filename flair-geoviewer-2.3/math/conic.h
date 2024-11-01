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
 * Date:	04-Feb-2010
 */

#ifndef __CONIC_H
#define __CONIC_H

#include <ostream>

#include "bmath.h"
#include "vector.h"

class Matrix3;

// Conic types
enum ConicType {
	CONIC_DEGENERATE,
	CONIC_POINT,
	CONIC_LINE,
	CONIC_LINES,
	CONIC_ELLIPSE,
	CONIC_PARABOLA,
	CONIC_HYPERBOLA,
	CONIC_NONE
};

/**
 * Generalized conic line
 * S = ax^2 + 2hxy + by^2 + 2gx + 2fy + c = 0
 *
 *        /x\T  / a h g \   /x\
 *    S = |y| * | h b f | * |y| = 0
 *        \1/   \ g f c /   \1/
 *
 * Invariants
 *        | a h g |
 *    D = | h b f |
 *        | g f c |
 *
 *    I = a + b
 *    J = a*b - h^2
 *    K = ca - g^2 + bc - f^2
 *
 * relative invariants
 *    I' = r*I
 *    J' = r^2*J
 *    D' = r^3*D
 */
class Conic {
public:
	double a, h, b, g, f, c;
private:
	ConicType _type;
	double D, I, J, K;
	double c1, c2, c3;
	double c4, c5, c6;

static	const char* _typeStr[];

public:
static	const double CONICPREC;

public:
	Conic() { set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0); }
	Conic(const double aa, const double ah, const double ab,
	      const double ag, const double af, const double ac)
		{ set(aa, ah, ab, ag, af, ac); }

	void	set(const double aa, const double ah, const double ab,
		    const double ag, const double af, const double ac);
	void	copy(const Conic& src);

	ConicType type() const	{ return _type; }
	void	matrix(Matrix3* m) const;

	/** @return the rotation of axes to eliminate the h component */
	double	theta() const	{ return 0.5 * atan2(2.0*h, a - b); }
	void	rotate(const double theta);
	void	translate2Origin(double* dx, double* dy) const;
	void	translate(const double dx, const double dy);
	void	splitLines(Conic *conic1, Conic *conic2);

	int	intersect(const Conic& conic, Vector2D pts[4]) const;

	void	parametric();
	double	operator () (const double x, const double y) const
		{ return (a*x + 2.0*(h*y + g))*x + (b*y + 2.0*f)*y + c; }

	double	getT(double x, double y) const;
	void	getXY(const double t, double *x, double *y) const;
	int	getY(const double x, double *y1, double *y2) const;

	void	grad(const double x, const double y, double *dx, double *dy) const;
	double	adist(const double x, const double y) const;
	bool	equal(const Conic& conic, const double acc) const;

const	char*	typeStr()		const	{ return _typeStr[_type]; }
static	const char* typeStr(ConicType t)	{ return _typeStr[t]; }

static	std::ostream& _Sx(std::ostream& s, const double a, const char *suffix);

private:
#if _DEBUG>1
	int	_intersect(const Conic& conic, Vector2D pts[4]) const;
#endif
	int	_intersectLine(const Conic& line, Vector2D pts[4]) const;
	bool	addPoint(const Conic& conic, int *n, Vector2D pts[4], double x, double y) const;

friend	std::ostream& operator << (std::ostream&, const Conic&);
}; // Conic

std::ostream& operator << (std::ostream&, const Conic&);

#endif
