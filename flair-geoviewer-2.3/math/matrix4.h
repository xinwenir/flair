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

#ifndef __MATRIX4_H
#define __MATRIX4_H

#include <ostream>
#include <assert.h>
#include <string.h>

#include "os.h"
#include "point.h"
#include "vector.h"

////////////////////// Matrix4 ////////////////////////
/**
 * Matrix4 4x4 array
 * WARNING data are uninitialized on constructor
 */
class Matrix4 {
private:
	double	_data[16];		// data
static	double	_identity[16];

public:
	Matrix4()			{}
	Matrix4(int type)		{ if (type) identity(); else zero(); }
	Matrix4(const Matrix4& mat)	{ copy(mat); }

	// Get/Set data
	int	rows()		const { return 4; }
	int	columns()	const { return 4; }

	/** @return data at [row, col] */
	double	operator () (const int row, const int col) const {
			assert(row<rows() && col<columns());
			return _data[row*columns()+col];
		}
	double&	operator () (const int row, const int col) {
			assert(row<rows() && col<columns());
			return _data[row*columns()+col];
		}

	/** copy matrix from src
	 * @param src matrix to copy from
	 */
	void	copy(const Matrix4& src)	{ memcpy(_data, src._data, sizeof(_data)); }
	bool	isIdentity()	const	{ return !memcmp(_data, _identity, sizeof(_data)); }
	Matrix4& operator =(const Matrix4& matrix) { copy(matrix); return *this; }

	void	zero()		{ bzero(_data, sizeof(_data)); }
	void	identity()	{ memcpy(_data, _identity, sizeof(_data)); }
	void	transpose();
	void	T() {transpose(); }

inline	void	multiply(const Matrix4& A, const Matrix4& B);
inline	void	multiplyT(const Matrix4& A, const Matrix4& B);
	Matrix4 operator *(const Matrix4& B)	{ Matrix4 r; r.multiply(*this, B); return r; }
	Matrix4& operator *=(const Matrix4& B)	{ Matrix4 r; r.multiply(*this, B); copy(r); return *this; }
	void	negate();

	//double det() const;
	void	inverse();
	void	inverse(const Matrix4& src) { invertMatrix(src._data, _data);}

	double	trace() const;

	// transform point in place coordinates (x,y,z) = Matrix * (x,y,z)
	void	transform(double *x, double *y, double *z) const {
			double ox, oy, oz;
			ox = *x*m(0,0) + *y*m(0,1) + *z*m(0,2) + m(0,3);
			oy = *x*m(1,0) + *y*m(1,1) + *z*m(1,2) + m(1,3);
			oz = *x*m(2,0) + *y*m(2,1) + *z*m(2,2) + m(2,3);
			*x = ox;
			*y = oy;
			*z = oz;
		}

	// transform in place vector V = Matrix * V
//	void	transform(Point& p)  const { transform(&p.x, &p.y, &p.z); }

	// FIXME remove the translation! Point 
	Vector	operator *(const Vector &v) const {
			return Vector(v.x*m(0,0) + v.y*m(0,1) + v.z*m(0,2),  // + m(0,3),
				      v.x*m(1,0) + v.y*m(1,1) + v.z*m(1,2),  // + m(1,3),
				      v.x*m(2,0) + v.y*m(2,1) + v.z*m(2,2)); // + m(2,3));
		}

	Point	operator *(const Point &v) const {
			return Point(v.x*m(0,0) + v.y*m(0,1) + v.z*m(0,2) + m(0,3),
				     v.x*m(1,0) + v.y*m(1,1) + v.z*m(1,2) + m(1,3),
				     v.x*m(2,0) + v.y*m(2,1) + v.z*m(2,2) + m(2,3));
		}

	Point	mult(const double x, const double y, const double z) const {
			return Vector(x*m(0,0) + y*m(0,1) + z*m(0,2) + m(0,3),
				      x*m(1,0) + y*m(1,1) + z*m(1,2) + m(1,3),
				      x*m(2,0) + y*m(2,1) + z*m(2,2) + m(2,3));
		}
	/** multiply with a vector - no translation is included - */
	Vector	multVector(const double x, const double y, const double z) const {
			return Vector(x*m(0,0) + y*m(0,1) + z*m(0,2),
				      x*m(1,0) + y*m(1,1) + z*m(1,2),
				      x*m(2,0) + y*m(2,1) + z*m(2,2));
		}

	Point	mult(const double x, const double y) const {
			return Vector(x*m(0,0) + y*m(0,1) + m(0,3),
				      x*m(1,0) + y*m(1,1) + m(1,3),
				      x*m(2,0) + y*m(2,1) + m(2,3));
		}
	Vector	multVector(const double x, const double y) const {
			return Vector(x*m(0,0) + y*m(0,1),
				      x*m(1,0) + y*m(1,1),
				      x*m(2,0) + y*m(2,1));
		}


	/** multiply with a vector - no translation is included - */
	Vector	multVector(const Vector &v) const {
			return Vector(v.x*m(0,0) + v.y*m(0,1) + v.z*m(0,2),
				      v.x*m(1,0) + v.y*m(1,1) + v.z*m(1,2),
				      v.x*m(2,0) + v.y*m(2,1) + v.z*m(2,2));
		}

//	Point	mult(const Point& p)  const { return mult(p.x, p.y, p.z); }
//	Vector	mult(const Vector& v) const { return mult(v.x, v.y, v.z); }

	// transform vector in place coordinates (x,y,z) = Matrix * (x,y,z) (NO TRANSLATION)
	void	transformVector(double *x, double *y, double *z) const {
			double ox, oy, oz;
			ox = *x*m(0,0) + *y*m(0,1) + *z*m(0,2);
			oy = *x*m(1,0) + *y*m(1,1) + *z*m(1,2);
			oz = *x*m(2,0) + *y*m(2,1) + *z*m(2,2);
			*x = ox;
			*y = oy;
			*z = oz;
		}

//	void	transform(Vector& v) const { transformVector(&v.x, &v.y, &v.z); }
	void	transform(Vector& v) const { transform(&v.x, &v.y, &v.z); }

	// Create translation matrix
	void	translate(const double x, const double y, const double z) {
			identity();
			(*this)(0,3) = x;
			(*this)(1,3) = y;
			(*this)(2,3) = z;
		}
	void	translate(const Vector& v) {
			identity();
			(*this)(0,3) = v.x;
			(*this)(1,3) = v.y;
			(*this)(2,3) = v.z;
		}
	void	scale(const double sx, const double sy=0.0, const double sz=0.0) {
			identity();
			(*this)(0,0) = sx;
			(*this)(1,1) = (sy==0.0?sx:sy);
			(*this)(2,2) = (sz==0.0?sx:sz);
		}

	// Create rotation matrix
	void	rotate(const double angle, const int axis);
	void	rotate(const double angle, const double x, const double y, const double z);
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

	/** @return determinate of the rotation 3x3 matrix only */
	double	det3() {
			return	  m(0,0)*m(1,1)*m(2,2)
				+ m(0,1)*m(1,2)*m(2,0)
				+ m(0,2)*m(1,0)*m(2,1)
				- m(0,2)*m(1,1)*m(2,0)
				- m(0,1)*m(1,0)*m(2,2)
				- m(0,0)*m(1,2)*m(2,1);
		}

	void	fix();

private:
	void	fix01();
	double	m(const int row, const int col) const {return (*this)(row,col);}
static	void	invertMatrixGeneral(const double *, double *);
static	void	invertMatrix(const double *, double *);
}; // Matrix4

std::ostream& operator << (std::ostream&, const Matrix4&);

/** multiply */
inline void Matrix4::multiply(const Matrix4& A, const Matrix4& B)
{
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			(*this)(i,j) =	A(i,0)*B(0,j) +
					A(i,1)*B(1,j) +
					A(i,2)*B(2,j) +
					A(i,3)*B(3,j);
} /* multiply */

/** multiplyT
 * multiply the transpose of A with B: A.T() * B
 */
inline void Matrix4::multiplyT(const Matrix4& A, const Matrix4& B)
{
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			(*this)(i,j) =	A(0,i)*B(0,j) +
					A(1,i)*B(1,j) +
					A(2,i)*B(2,j) +
					A(3,i)*B(3,j);
} /* multiplyT */

#endif
