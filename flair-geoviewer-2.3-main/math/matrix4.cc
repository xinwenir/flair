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

#include <stdio.h>

#include <iomanip>
#include <iostream>

#include "os.h"
#include "bmath.h"
#include "matrix4.h"

using namespace std;

double Matrix4::_identity[16] = {
	1.0, 0.0, 0.0, 0.0,
	0.0, 1.0, 0.0, 0.0,
	0.0, 0.0, 1.0, 0.0,
	0.0, 0.0, 0.0, 1.0
};

/** transpose */
void Matrix4::transpose()
{
	for (int i=0; i<4; i++)
		for (int j=i+1; j<4; j++) {
			double t = (*this)(i,j);
			(*this)(i, j) = (*this)(j, i);
			(*this)(j, i) = t;
		}
} // identity

/** operator * *
Matrix4& Matrix4::operator *= (const Matrix4& mat)
{
	// FIXME Optimize
	return (*this) = (*this)*mat;
} * operator *= */

/** negate matrix */
void Matrix4::negate()
{
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			(*this)(i,j) = -(*this)(i,j);
} // negate

/** det *
double Matrix4::det() const
{
	if (4 == 2)
		return m(0,0)*m(1,1) - m(1,0)*m(0,1);
	else
	if (4 == 3)
		return m(0,0)*(m(1,1)*m(2,2) - m(2,1)*m(1,2))
		     - m(0,1)*(m(1,0)*m(2,2) - m(2,0)*m(1,2))
		     + m(0,2)*(m(1,0)*m(2,1) - m(2,0)*m(1,1));
	else
		return 0.0;
} * det */

/** return trace of matrix (sum of diagonal elements) */
double Matrix4::trace() const
{
	double t = 0.0;
	for (int i=0; i<4; i++)
		t += m(i,i);
	return t;
} // trace

/** inverse the matrix. */
void Matrix4::inverse()
{
	Matrix4 inv;
	invertMatrix(_data, inv._data);
	copy(inv);
} // inverse

/** rotate matrix around an axis
 * @param angle
 * @param axis 0=X, 1=Y, 2=Z
 *
 * WARNING possibly a call to fix should be made after
 */
void Matrix4::rotate(const double angle, const int axis)
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
void Matrix4::rotate(const double angle, double x, double y, double z)
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

//
// From Mesa-2.2\src\glu\project.c
//
//
// Compute the inverse of a 4x4 matrix.  Contributed by scotter@lafn.org
//
#define MAT(m,r,c) (m)[(c)*4+(r)]
/* Here's some shorthand converting standard (row,column) to index. */
#define m11 MAT(m,0,0)
#define m12 MAT(m,0,1)
#define m13 MAT(m,0,2)
#define m14 MAT(m,0,3)
#define m21 MAT(m,1,0)
#define m22 MAT(m,1,1)
#define m23 MAT(m,1,2)
#define m24 MAT(m,1,3)
#define m31 MAT(m,2,0)
#define m32 MAT(m,2,1)
#define m33 MAT(m,2,2)
#define m34 MAT(m,2,3)
#define m41 MAT(m,3,0)
#define m42 MAT(m,3,1)
#define m43 MAT(m,3,2)
#define m44 MAT(m,3,3)

void Matrix4::invertMatrixGeneral(const double *m, double *out)
{
	double det;
	double d12, d13, d23, d24, d34, d41;
	double tmp[16]; /* Allow out == in. */

	/* Inverse = adjoint / det. (See linear algebra texts.)*/
	/* pre-compute 2x2 dets for last two rows when computing */
	/* cofactors of first two rows. */
	d12 = (m31*m42-m41*m32);
	d13 = (m31*m43-m41*m33);
	d23 = (m32*m43-m42*m33);
	d24 = (m32*m44-m42*m34);
	d34 = (m33*m44-m43*m34);
	d41 = (m34*m41-m44*m31);

	tmp[0] =  (m22 * d34 - m23 * d24 + m24 * d23);
	tmp[1] = -(m21 * d34 + m23 * d41 + m24 * d13);
	tmp[2] =  (m21 * d24 + m22 * d41 + m24 * d12);
	tmp[3] = -(m21 * d23 - m22 * d13 + m23 * d12);

	/* Compute determinant as early as possible using these cofactors. */
	det = m11 * tmp[0] + m12 * tmp[1] + m13 * tmp[2] + m14 * tmp[3];

	/* Run singularity test. */
	if (det == 0.0) {
		/* printf("invert_matrix: Warning: Singular matrix.\n"); */
		memcpy(out,_identity,16*sizeof(double));
	} else {
		double invDet = 1.0 / det;

		/* Compute rest of inverse. */
		tmp[0] *= invDet;
		tmp[1] *= invDet;
		tmp[2] *= invDet;
		tmp[3] *= invDet;

		tmp[4] = -(m12 * d34 - m13 * d24 + m14 * d23) * invDet;
		tmp[5] =  (m11 * d34 + m13 * d41 + m14 * d13) * invDet;
		tmp[6] = -(m11 * d24 + m12 * d41 + m14 * d12) * invDet;
		tmp[7] =  (m11 * d23 - m12 * d13 + m13 * d12) * invDet;

		/* Pre-compute 2x2 dets for first two rows when computing */
		/* cofactors of last two rows. */
		d12 = m11*m22-m21*m12;
		d13 = m11*m23-m21*m13;
		d23 = m12*m23-m22*m13;
		d24 = m12*m24-m22*m14;
		d34 = m13*m24-m23*m14;
		d41 = m14*m21-m24*m11;

		tmp[8]  =  (m42 * d34 - m43 * d24 + m44 * d23) * invDet;
		tmp[9]  = -(m41 * d34 + m43 * d41 + m44 * d13) * invDet;
		tmp[10] =  (m41 * d24 + m42 * d41 + m44 * d12) * invDet;
		tmp[11] = -(m41 * d23 - m42 * d13 + m43 * d12) * invDet;
		tmp[12] = -(m32 * d34 - m33 * d24 + m34 * d23) * invDet;
		tmp[13] =  (m31 * d34 + m33 * d41 + m34 * d13) * invDet;
		tmp[14] = -(m31 * d24 + m32 * d41 + m34 * d12) * invDet;
		tmp[15] =  (m31 * d23 - m32 * d13 + m33 * d12) * invDet;

		memcpy(out, tmp, 16*sizeof(double));
	}
} // invertMatrixGeneral

//
// Invert matrix m.  This algorithm contributed by Stephane Rehel
// <rehel@worldnet.fr>
//
void Matrix4::invertMatrix(const double *m, double *out)
{
	register double det;
	double tmp[16]; /* Allow out == in. */

	if( m41 != 0.0 || m42 != 0.0 || m43 != 0.0 || m44 != 1.0 ) {
		invertMatrixGeneral(m, out);
		return;
	}

	/* Inverse = adjoint / det. (See linear algebra texts.)*/
	tmp[0]= m22 * m33 - m23 * m32;
	tmp[1]= m23 * m31 - m21 * m33;
	tmp[2]= m21 * m32 - m22 * m31;

	/* Compute determinant as early as possible using these cofactors. */
	det= m11 * tmp[0] + m12 * tmp[1] + m13 * tmp[2];

	/* Run singularity test. */
	if (det == 0.0) {
		/* printf("invert_matrix: Warning: Singular matrix.\n"); */
		memcpy( out, _identity, 16*sizeof(double) );
	} else {
		double d12, d13, d23, d24, d34, d41;
		register double im11, im12, im13, im14;

		det= 1. / det;

		/* Compute rest of inverse. */
		tmp[0] *= det;
		tmp[1] *= det;
		tmp[2] *= det;
		tmp[3]  = 0.;

		im11= m11 * det;
		im12= m12 * det;
		im13= m13 * det;
		im14= m14 * det;
		tmp[4] = im13 * m32 - im12 * m33;
		tmp[5] = im11 * m33 - im13 * m31;
		tmp[6] = im12 * m31 - im11 * m32;
		tmp[7] = 0.;

		/* Pre-compute 2x2 dets for first two rows when computing */
		/* cofactors of last two rows. */
		d12 = im11*m22 - m21*im12;
		d13 = im11*m23 - m21*im13;
		d23 = im12*m23 - m22*im13;
		d24 = im12*m24 - m22*im14;
		d34 = im13*m24 - m23*im14;
		d41 = im14*m21 - m24*im11;

		tmp[8]  =  d23;
		tmp[9]  = -d13;
		tmp[10] = d12;
		tmp[11] = 0.0;

		tmp[12] = -(m32 * d34 - m33 * d24 + m34 * d23);
		tmp[13] =  (m31 * d34 + m33 * d41 + m34 * d13);
		tmp[14] = -(m31 * d24 + m32 * d41 + m34 * d12);
		tmp[15] =  1.0;

		memcpy(out, tmp, 16*sizeof(double));
	}
} // invertMatrix
#undef m11
#undef m12
#undef m13
#undef m14
#undef m21
#undef m22
#undef m23
#undef m24
#undef m31
#undef m32
#undef m33
#undef m34
#undef m41
#undef m42
#undef m43
#undef m44
#undef MAT

#define SMALL 1e-10

/** fix matrix to be a valid rotation matrix */
void Matrix4::fix01()
{
	// round to 0.0, 1.0 and -1.0 if needed with precision
	// of 1e-10 since 1e-11**2 = limit of precision
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
//			if (Eq0((*this)(i,j), SMALL))
//				(*this)(i,j) = 0.0;
//			else
//			if (Eq((*this)(i,j), 1.0, SMALL)) {

			// Something that should never be done for floating point
			// however up to phi=1.06e-8, cos(phi)==1.0 exactly
			// and sin(phi)>0.0, so cos(phi)**2 + sin(phi)**2 > 1.0
			if ((*this)(i,j) == 1.0) {
//				(*this)(i,j) = 1.0;
				// Set the other two to zero
				(*this)(i,(j+1)%3) = 0.0;
				(*this)(i,(j+2)%3) = 0.0;
				break;
			} else
//			if (Eq((*this)(i,j), -1.0, SMALL)) {
			if ((*this)(i,j) == -1.0) {
//				(*this)(i,j) = -1.0;
				// Set the other two to zero
				(*this)(i,(j+1)%3) = 0.0;
				(*this)(i,(j+2)%3) = 0.0;
				break;
			}
		}
	}
} // fix

/** fix matrix to be a valid rotation matrix */
void Matrix4::fix()
{
	fix01();

	// Try to correct (*this) for possible numerical precision problems
	// normalize vectors
	double len = 1.0 / sqrt(Sqr((*this)(0,0)) + Sqr((*this)(0,1)) + Sqr((*this)(0,2)));
	(*this)(0,0) *= len;
	(*this)(0,1) *= len;
	(*this)(0,2) *= len;

	len = 1.0 / sqrt(Sqr((*this)(1,0)) + Sqr((*this)(1,1)) + Sqr((*this)(1,2)));
	(*this)(1,0) *= len;
	(*this)(1,1) *= len;
	(*this)(1,2) *= len;

	// FIXME Check orthogonality of row0 and row1
	// force orthogonality row[2] = row[0] x row[1]
	(*this)(2,0) = (*this)(0,1)*(*this)(1,2) - (*this)(1,1)*(*this)(0,2),
	(*this)(2,1) = (*this)(0,2)*(*this)(1,0) - (*this)(1,2)*(*this)(0,0),
	(*this)(2,2) = (*this)(0,0)*(*this)(1,1) - (*this)(1,0)*(*this)(0,1);

	len = 1.0 / sqrt(Sqr((*this)(2,0)) + Sqr((*this)(2,1)) + Sqr((*this)(2,2)));
	(*this)(2,0) *= len;
	(*this)(2,1) *= len;
	(*this)(2,2) *= len;

	// force last row
	(*this)(3,0) = 0.0;
	(*this)(3,1) = 0.0;
	(*this)(3,2) = 0.0;
	(*this)(3,3) = 1.0;

	fix01();
} // fix

/** matrix print out */
ostream& operator << (ostream& s, const Matrix4& matrix)
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
			s << " " << setw(22) << setprecision(16) << matrix(r,c);
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
