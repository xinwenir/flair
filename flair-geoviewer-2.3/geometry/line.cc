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
 *
 * Date:	16-Feb-2016
 */

#include "os.h"
#include "eps.h"
#include "line.h"

/* --- pointLineDistance --- */
bool pointLineDistance( const Point& point,
		const Point& start, const Point& end,
		double *distance )
{
	Vector dir = end-start;
	double u = ((point-start) * dir) / dir.length2();
	if (u<0.0 || u>1.0)
		return false;   // closest point does not fall within the dir segment

	if (distance!=NULL) {
		Point intersection = start + u*dir;
		*distance = (point-intersection).length();
	}
	return true;
} // PointLineDistance

/**
 * Case 1. parallel lines
 *	return PointLineDistance of p1 wrt to line2
 * Case 2. skew lines
 * The distance between two skew lines with equations
 *      y = p1 + z1*s
 *      y = p2 + z2*s
 *   is given as
 *          |(p2-p1) . (z1 x z2)|
 *      D = ---------------------
 *              |z1 x z2|
 * @return the minimum distance from two lines
 */
double lineLineDistance(const Point& p1, const Vector& z1,
			const Point& p2, const Vector& z2)
{
	Vector axb = z1 ^ z2;

	double len_axb = axb.length();
	if (ISSMALL(len_axb)) { // parallel lines
		// return point line distance
		double u = ((p1-p2) * z2) / z2.length2();
		Point intersection = p2 + u*z2;
		return (p1-intersection).length();
	} else	// skew lines
		return Abs((p2-p1) * axb) / Abs(len_axb);
} // lineLineDistance

/* --- LineLineIntersect --- */
bool lineLineIntersect(
		const Point& p1, const Point& p2,
		const Point& p3, const Point& p4,
		Point *pa, Point *pb,
		double *mua,  double *mub, const double eps)
{
	// based on http://astronomy.swin.edu.au/~pbourke/geometry/lineline3d/
	Vector p13 = p1 - p3;
	Vector p43 = p4 - p3;

	if (Abs(p43.x) < eps && Abs(p43.y) < eps && Abs(p43.z) < eps)
		return false;

	Vector p21 = p2 - p1;

	if (Abs(p21.x) < eps && Abs(p21.y) < eps && Abs(p21.z) < eps)
		return false;

	double d1343 = p13 * p43;	// dot product
	double d4321 = p43 * p21;
	double d1321 = p13 * p21;
	double d4343 = p43 * p43;
	double d2121 = p21 * p21;
	double denom = d2121*d4343 - d4321*d4321;

	if (Abs(denom) < eps)
		return false;

	*mua = (d1343*d4321 - d1321*d4343) / denom;
	*mub = (d1343 + d4321*(*mua)) / d4343;

	if (pa) *pa = p1 + *mua * p21;
	if (pb) *pb = p3 + *mub * p43;

	return true;
} // lineLineIntersect
