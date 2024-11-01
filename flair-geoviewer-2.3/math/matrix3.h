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

#ifndef __MATRIX3_H
#define __MATRIX3_H

#include <ostream>
#include <string.h>

#include "os.h"
#include "vector.h"

////////////////////// Matrix3 ////////////////////////
/**
 * Matrix3 3x3 array
 */
class Matrix3 {
private:
	double	_data[9];			// data
static	double	_identity[9];

public:
	Matrix3()			{}
	Matrix3(int type)		{ if (type==0) zero(); else identity(); }
	// XXX Warning data are uninitialized XXX
	Matrix3(const Matrix3& mat)	{ copy(mat); }

	// Get/Set data
	int	rows()		const { return 3; }
	int	columns()	const { return 3; }

	/** @return data at [row, col] */
	double	operator () (const int row, const int col) const {
			assert(row<rows() && col<columns());
			return _data[row*columns()+col];
		}
	double&	operator () (const int row, const int col) {
			assert(row<columns() && col<columns());
			return _data[row*columns()+col];
		}

	/** copy matrix from src
	 * @param src matrix to copy from
	 */
	void	copy(const Matrix3& src)	{ memcpy(_data, src._data, sizeof(_data)); }
	bool	isIdentity()	const	{ return !memcmp(_data, _identity, sizeof(_data)); }
	Matrix3& operator =(const Matrix3& matrix) { copy(matrix); return *this; }

	void	zero()		{ bzero(_data, sizeof(_data)); }
	void	identity()	{ memcpy(_data, _identity, sizeof(_data)); }
	void	transpose();
	void	T() {transpose(); }

inline	void	multiply(const Matrix3& A, const Matrix3& B);
inline	void	multiplyT(const Matrix3& A, const Matrix3& B);
//	Matrix3& operator *= (const Matrix3&);
	void	negate();

	bool	inverse(const double eps=Vector::_epsilon);
//	bool	inverse(const Matrix3& src) { invertMatrix(src._data, _data);}

	double	trace() const;

	Vector	operator *(const Vector &v) const {
			return Vector(v.x*m(0,0) + v.y*m(0,1) + v.z*m(0,2),
				      v.x*m(1,0) + v.y*m(1,1) + v.z*m(1,2),
				      v.x*m(2,0) + v.y*m(2,1) + v.z*m(2,2));
		}
	Vector	mult(const double x, const double y, const double z) const {
			return Vector(x*m(0,0) + y*m(0,1) + z*m(0,2),
				      x*m(1,0) + y*m(1,1) + z*m(1,2),
				      x*m(2,0) + y*m(2,1) + z*m(2,2));
		}
	Vector	mult(const double x, const double y) const {
			return Vector(x*m(0,0) + y*m(0,1),
				      x*m(1,0) + y*m(1,1),
				      x*m(2,0) + y*m(2,1));
		}
	Vector	mult(const Vector& v) const { return mult(v.x, v.y, v.z); }

	void	scale(const double sx, const double sy=0.0, const double sz=0.0) {
			identity();
			(*this)(0,0) = sx;
			(*this)(1,1) = (sy==0.0?sx:sy);
			(*this)(2,2) = (sz==0.0?sx:sz);
		}

	void	rotate(const double angle, const int axis);
	void	rotate(const double angle, const double x, const double y,
				const double z);
	void	rotate(const double angle, const Vector& axis)
			{ rotate(angle, axis.x, axis.y, axis.z); }

	void	rotX(const double angle)	{rotate(angle,0);}
	void	rotY(const double angle)	{rotate(angle,1);}
	void	rotZ(const double angle)	{rotate(angle,2);}

	void	make(const Vector& vx, const Vector& vy, const Vector& vz) {
			identity();
			for (int i=0; i<3; i++) {
				(*this)(0,i) = vx[i];
				(*this)(1,i) = vy[i];
				(*this)(2,i) = vz[i];
			}
		}

	/** @return determinate */
	double	det() const {
			return	  m(0,0)*m(1,1)*m(2,2)
				+ m(0,1)*m(1,2)*m(2,0)
				+ m(0,2)*m(1,0)*m(2,1)
				- m(0,2)*m(1,1)*m(2,0)
				- m(0,1)*m(1,0)*m(2,2)
				- m(0,0)*m(1,2)*m(2,1);
		}

	void	fix();

private:
	double	m(const int row, const int col) const {return (*this)(row,col);}
}; // Matrix3

std::ostream& operator << (std::ostream&, const Matrix3&);

/** multiply */
inline void Matrix3::multiply(const Matrix3& A, const Matrix3& B)
{
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			(*this)(i,j) =	A(i,0)*B(0,j) +
					A(i,1)*B(1,j) +
					A(i,2)*B(2,j);
} /* multiply */

/** multiplyT
 * multiply the transpose of A with B: A.T() * B
 */
inline void Matrix3::multiplyT(const Matrix3& A, const Matrix3& B)
{
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			(*this)(i,j) =	A(0,i)*B(0,j) +
					A(1,i)*B(1,j) +
					A(2,i)*B(2,j);
} /* multiplyT */

#endif
