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

#include <stdlib.h>
#include <ostream>

#include "bmath.h"
#include "vector.h"
#include "matrix4.h"

using namespace std;

const Vector Vector::O (0.0, 0.0, 0.0);
const Vector Vector::Xo(1.0, 0.0, 0.0);
const Vector Vector::Yo(0.0, 1.0, 0.0);
const Vector Vector::Zo(0.0, 0.0, 1.0);

double Vector::_epsilon = 1E-10;

#if 1
/* --- operator *= --- */
Vector& Vector::operator *=(const Matrix4& m)
{
	return *this = m * (*this);
} // operator *=

/* --- transform --- */
Vector& Vector::transform(const Matrix4& m)
{
	return *this = m * (*this);
} // transform
#endif

/* --- rotateX --- */
void Vector::rotateX(double ang)
{
	double s, c;
	bsincos(ang, &s, &c);
	double yy = y;
	y = c*yy - s*z;
	z = s*yy + c*z;
} // rotateX

/* --- rotateY --- */
void Vector::rotateY(double ang)
{
	double s, c;
	bsincos(ang, &s, &c);
	double zz = z;
	z = c*zz - s*x;
	x = s*zz + c*x;
} // rotateY

/* --- rotateY --- */
void Vector::rotateZ(double ang)
{
	double s, c;
	bsincos(ang, &s, &c);
	double xx = x;
	x = c*xx - s*y;
	y = s*xx + c*y;
} // rotateZ

#if 0
/* --- rotate --- */
void Vector::rotate(double angle, const Vector& axis)
{
	Rotate trans(ang,axis);
	operator*=(trans);
} // rotate
#endif

/* --- rotateUz --- */
void Vector::rotateUz(Vector& newUzVector)
{
	// newUzVector must be normalized !

	double u1 = newUzVector.x;
	double u2 = newUzVector.y;
	double u3 = newUzVector.z;
	double up = u1*u1 + u2*u2;

	if (Eq0(up,1.0E-20)) {
		up = sqrt(up);
		double px = x,  py = y,  pz = z;
		x = (u1*u3*px - u2*py + u1*up*pz)/up;
		y = (u2*u3*px + u1*py + u2*up*pz)/up;
		z = (u3*u3*px -    px + u3*up*pz)/up;
	}
	else if (u3 < 0.0) {
		x = -x;		// phi=0  theta=pi
		z = -z;
	} else {
		// nop for the moment
	};
} // rotateUz

/** return direction of a unit vector
 * @return	0 - unknown
 *		[+/-]1 - X
 *		[+/-]2 - Y
 *		[+/-]3 - Z
 */
int Vector::direction(const double acc) const
{
	if (Eq(x, 1.0, acc))
		return  1;
	if (Eq(x,-1.0, acc))
		return -1;

	if (Eq(y, 1.0, acc))
		return  2;
	if (Eq(y,-1.0, acc))
		return -2;

	if (Eq(z, 1.0, acc))
		return  3;
	if (Eq(z,-1.0, acc))
		return -3;

	return 0;
} // direction

/* compare - order vectors in ascending z first then y then x order,
 *  useful for sorted arrays
 */
int Vector::compare(const Vector* a, const Vector* b)
{
	if (a->z > b->z)
		return  1;
	else
	if (a->z < b->z)
		return -1;
	else
	if (a->y > b->y)
		return  1;
	else
	if (a->y < b->y)
		return -1;
	else
	if (a->x > b->x)
		return  1;
	else
	if (a->x < b->x)
		return -1;
	return 0;
} // compare

/* --- operator << --- */
/**
 * output a vector as an string to an output stream
 */
ostream& operator << (ostream& s, const Vector& q)
{
	return s << "[" << q.x << ", " << q.y << ", " << q.z << "]";
} /* operator << */
