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

#include <math.h>
#include <assert.h>
#include <string.h>

#include <iomanip>
#include <iostream>

#include "os.h"
#include "bmath.h"
#include "matrix3.h"

using namespace std;

double Matrix3::_identity[9] = {
	1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0,
};

/** transpose */
void Matrix3::transpose()
{
	for (int i=0; i<3; i++)
		for (int j=i+1; j<3; j++) {
			double t = (*this)(i,j);
			(*this)(i, j) = (*this)(j, i);
			(*this)(j, i) = t;
		}
} // transpose

/** operator * *
Matrix3& Matrix3::operator *= (const Matrix3& mat)
{
	// FIXME Optimize
	return (*this) = (*this)*mat;
} * operator *= */

/** negate matrix */
void Matrix3::negate()
{
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			(*this)(i,j) = -(*this)(i,j);
} // negate

/** return trace of matrix (sum of diagonal elements) */
double Matrix3::trace() const
{
	double t = 0.0;
	for (int i=0; i<3; i++)
		t += m(i,i);
	return t;
} // trace

/** invert the matrix in place and
 * @return true on success, false otherwise
 */
bool Matrix3::inverse(const double eps)
{
	Matrix3 inv;
	double D = det();
	if (Eq0(D, eps)) return false;
	D = 1.0/D;

	inv(0,0) =  ((*this)(1,1)*(*this)(2,2) - (*this)(2,1)*(*this)(1,2)) * D;
	inv(0,1) = -((*this)(0,1)*(*this)(2,2) - (*this)(0,2)*(*this)(2,1)) * D;
	inv(0,2) =  ((*this)(0,1)*(*this)(1,2) - (*this)(0,2)*(*this)(1,1)) * D;
	inv(1,0) = -((*this)(1,0)*(*this)(2,2) - (*this)(1,2)*(*this)(2,0)) * D;
	inv(1,1) =  ((*this)(0,0)*(*this)(2,2) - (*this)(0,2)*(*this)(2,0)) * D;
	inv(1,2) = -((*this)(0,0)*(*this)(1,2) - (*this)(1,0)*(*this)(0,2)) * D;
	inv(2,0) =  ((*this)(1,0)*(*this)(2,1) - (*this)(2,0)*(*this)(1,1)) * D;
	inv(2,1) = -((*this)(0,0)*(*this)(2,1) - (*this)(2,0)*(*this)(0,1)) * D;
	inv(2,2) =  ((*this)(0,0)*(*this)(1,1) - (*this)(1,0)*(*this)(0,1)) * D;
	copy(inv);
	return true;
} // inverse

/** rotate matrix around an axis
 * @param angle
 * @param axis 0=X, 1=Y, 2=Z
 * */
void Matrix3::rotate(const double angle, const int axis)
{
	double s, c;
	identity();

	bsincos(angle, &s, &c);

	int m1 = ((axis+1)%3)+1;
	int m2 = m1%3;
	m1 -= 1;

	(*this)(m1,m1) =  c;
	(*this)(m2,m2) =  c;
	(*this)(m1,m2) = -s;
	(*this)(m2,m1) =  s;
} // rotate

/** rotate matrix around a vector */
void Matrix3::rotate(const double angle, double x, double y, double z)
{
	double s, c;
	identity();

	bsincos(angle, &s, &c);

	double invl = 1.0/sqrt(x*x + y*y + z*z);
	x *= invl;
	y *= invl;
	z *= invl;

	double c1 = 1.0 - c;
	(*this)(0,0) = x*x + (1-x*x)*c;
	(*this)(0,1) = x*y*c1 - z*s;
	(*this)(0,2) = x*z*c1 + y*s;

	(*this)(1,0) = x*y*c1 + z*s;
	(*this)(1,1) = y*y + (1-y*y)*c;
	(*this)(1,2) = y*z*c1 - x*s;

	(*this)(2,0) = x*z*c1 - y*s;
	(*this)(2,1) = y*z*c1 + x*s;
	(*this)(2,2) = z*z + (1-z*z)*c;
} // rotate

/** matrix print out */
ostream& operator << (ostream& s, const Matrix3& matrix)
{
	for (int r=0; r<matrix.rows(); r++) {
		if (r==0)
			s << "/";
		else
		if (r+1==matrix.rows())
			s << "\\";
		else
			s << "|";
		for (int c=0; c<matrix.columns(); c++) {
			s << " " << setw(16) << setprecision(10) << matrix(r,c);
		}
		if (r==0)
			s << " \\" << endl;
		else
		if (r+1==matrix.rows())
			s << " /" << endl;
		else
			s << " |" << endl;
	}
	return s;
} /* operator << */
