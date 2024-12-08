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
 * Date:	14-Jun-2013
 */

#ifndef __MATRIX2_H
#define __MATRIX2_H

#include <ostream>
#include <string.h>

#include "os.h"
#include "vector.h"

////////////////////// Matrix2 ////////////////////////
/**
 * Matrix2 2x2 array
 */
class Matrix2 {
private:
	double	_data[4];			// data
static	double	_identity[4];

public:
	Matrix2()			{}
	Matrix2(int type)		{ if (type==0) zero(); else identity(); }
	Matrix2(const Matrix2& mat)	{ copy(mat); }

	// Get/Set data
	int	rows()		const { return 2; }
	int	columns()	const { return 2; }

	/** @return data at [row, col] */
	double	operator () (const int row, const int col) const {
			assert(row<rows() && col<columns());
			return _data[row*rows()+col];
		}
	double&	operator () (const int row, const int col) {
			assert(row<rows() && col<columns());
			return _data[row*rows()+col];
		}

	/** copy matrix from src
	 * @param src matrix to copy from
	 */
	void	copy(const Matrix2& src)	{ memcpy(_data, src._data, sizeof(_data)); }
	bool	isIdentity()	const	{ return !memcmp(_data, _identity, sizeof(_data)); }
	Matrix2& operator =(const Matrix2& matrix) { copy(matrix); return *this; }

	void	zero()		{ bzero(_data, sizeof(_data)); }
	void	identity()	{ memcpy(_data, _identity, sizeof(_data)); }
	void	transpose();
	void	T() {transpose(); }

inline	void	multiply(const Matrix2& A, const Matrix2& B);
inline	void	multiplyT(const Matrix2& A, const Matrix2& B);
//	Matrix2& operator *= (const Matrix2&);
	void	negate();

	bool	inverse(const double eps=Vector::_epsilon);
//	bool	inverse(const Matrix2& src) { invertMatrix(src._data, _data);}

	double	trace() const	{ return (*this)(0,0) + (*this)(1,1); }

	Vector	operator *(const Vector &v) const {
			return Vector(v.x*m(0,0) + v.y*m(0,1),
				      v.x*m(1,0) + v.y*m(1,1),
				      v.z);
		}
	Vector	mult(const double x, const double y) const {
			return Vector(x*m(0,0) + y*m(0,1),
				      x*m(1,0) + y*m(1,1),
				      0.0);
		}
	Vector	mult(const Vector& v) const { return mult(v.x, v.y); }

	void	scale(const double sx, const double sy=0.0) {
			identity();
			(*this)(0,0) = sx;
			(*this)(1,1) = (sy==0.0?sx:sy);
		}

	void	rotate(const double angle);

	/** @return determinate */
	double	det()	const { return m(0,0)*m(1,1) - m(0,1)*m(1,0); }

	void	fix();

private:
	double	m(const int row, const int col) const {return (*this)(row,col);}
}; // Matrix2

std::ostream& operator << (std::ostream&, const Matrix2&);

/** multiply */
inline void Matrix2::multiply(const Matrix2& A, const Matrix2& B)
{
	(*this)(0,0) = A(0,0)*B(0,0) + A(0,1)*B(1,0);
	(*this)(0,1) = A(0,0)*B(0,1) + A(0,1)*B(1,1);

	(*this)(1,0) = A(1,0)*B(0,0) + A(1,1)*B(1,0);
	(*this)(1,1) = A(1,0)*B(0,1) + A(1,1)*B(1,1);
} /* multiply */

/** multiplyT
 * multiply the transpose of A with B: A.T() * B
 */
inline void Matrix2::multiplyT(const Matrix2& A, const Matrix2& B)
{
	(*this)(0,0) = A(0,0)*B(0,0) + A(1,0)*B(1,0);
	(*this)(0,1) = A(0,0)*B(0,1) + A(1,0)*B(1,1);

	(*this)(1,0) = A(0,1)*B(0,0) + A(1,1)*B(1,0);
	(*this)(1,1) = A(0,1)*B(0,1) + A(1,1)*B(1,1);
} /* multiplyT */
#endif
