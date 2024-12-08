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
 * Date:	19-Feb-2013
 */

#include "random.h"

/**
 * return a random angle sin and cos from 0..2PI
 * @param s	random sin
 * @param c	random cosine
 */
void Random::sincos(double *s, double *c)
{
	double x,y,x2,y2,d;
	do {
		x  = 2.0*real() - 1.0;
		y  = 2.0*real() - 1.0;
		x2 = Sqr(x);
		y2 = Sqr(y);
		d  = x2 + y2;
	} while (d > 1.0);

	*s = 2.0*x*y / d;
	*c = (x2 - y2) / d;
} // sincos

/**
 * @return a normal distributed random number generator
 * based on particle data book
 */
double Random::normal()
{
	if (_normalSave>=0.0) {
		double s = _normalSave;
		_normalSave = -1.0;
		return s;
	}

	double s, c;
	::sincos(PI2*real(), &s, &c);

	double sq = sqrt(-2.0*log(real()));
	_normalSave = s * sq;

	return c * sq;
} // normal

/**
 * return 2 normal distributed random numbers
 * @param z1,z2
 */
void Random::normal(double *z1, double *z2)
{
#if 0
	// Slower variant
	double r1 = PI2*real();
	double r2 = sqrt(-2.0*log(real()));
	*z1 = sin(r1) * r2;
	*z2 = cos(r2) * r2;
#endif
	double u1, u2;
	double r2;

	do {	// calculate sin/cos like above method sincos()
		u1 = 2.0*real() - 1.0;
		u2 = 2.0*real() - 1.0;
		r2 = u1*u1 + u2*u2;
	} while (r2 > 1.0);

	r2 = sqrt(-2.0*log(r2)/r2);
	*z1 = u1 * r2;
	*z2 = u2 * r2;
} // normal

/** sample uniform x,y coordinates inside a disc or radius r
 * @param r	radius of disc
 * @param x,y	return position
 */
void Random::disc(const double r, double *x, double *y)
{
	double R = r * sqrt(real());
	double s,c;
	sincos(&s, &c);
	*x = R * c;
	*y = R * s;
} // disc

/** sample uniform x,y,z coordinates inside a sphere of radius r
 * @param r	radius of sphere
 * @param x,y,z	return position
 */
void Random::sphere(const double r, double *x, double *y, double *z)
{
	double R = r * pow(real(), 1.0/3.0);

	double sphi, cphi;
	sincos(&sphi, &cphi);

	double cth = real(-1.0, 1.0);
	double sth = sqrt(1.0 - cth*cth);

	*x = R * sth * cphi;
	*y = R * sth * sphi;
	*z = R * cth;
} // sphere

/**
 * create an isotropic random vector of unit length using a
 * rejection technique
 * @return *this
 */
Vector Random::vector()
{
	double x, y, z;
	double len;

	do {
		x = 2.0*real() - 1.0;
		y = 2.0*real() - 1.0;
		z = 2.0*real() - 1.0;
		len = x*x + y*y + z*z;
	} while (len > 1.0);

	len = sqrt(1.0/len);
	x *= len;
	y *= len;
	z *= len;
	return Vector(x,y,z);
} // vector
