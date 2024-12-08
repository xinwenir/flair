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
#include "matrix.h"

using namespace std;

/* maximum size for matrix inversion */
#define MAXSIZE 10

/** make */
void Matrix::make(int r, int c)
{
	assert(r>0 && c>=0);
	if (c==0) c=r;
	if (_rows!=r || _cols!=c) {
		if (_data) delete [] _data;
		_rows = r;
		_cols = c;
		_data = new double[_rows*_cols];
	}
} // make

/** copy matrix from src
 * @param src matrix to copy from
 */
void Matrix::copy(const Matrix& src)
{
	make(src._rows, src._cols);
	memcpy(_data, src._data, _rows*_cols*sizeof(double));
} // copy

/** zero */
void Matrix::zero()
{
	for (int i=0; i<_rows; i++)
		for (int j=0; j<_cols; j++)
			(*this)(i,j) = 0.0;
} // zero

/** identity */
void Matrix::identity()
{
	for (int i=0; i<_rows; i++)
		for (int j=0; j<_cols; j++)
			(*this)(i,j) = i==j? 1.0 : 0.0;
} // identity

/** transpose */
void Matrix::transpose()
{
	if (_rows==_cols)
		for (int i=0; i<_rows; i++)
			for (int j=i+1; j<_cols; j++) {
				double t = (*this)(i,j);
				(*this)(i, j) = (*this)(j, i);
				(*this)(j, i) = t;
		}
	else {
		// FIXME could be done in place
		double *newdata = new double[_rows*_cols];
		for (int i=0; i<_cols; i++)
			for (int j=0; j<_rows; j++)
				newdata[i*_rows+j] = (*this)(j,i);
		delete [] _data;
		_data = newdata;
		Swap(_rows, _cols);
	}
} // transpose

/** transpose */
void Matrix::transpose(const Matrix& A)
{
	make(A._cols, A._rows);
	for (int i=0; i<_rows; i++)
		for (int j=0; j<_cols; j++)
			(*this)(i,j) = A(j,i);
} // transpose

/** multiply */
void Matrix::multiply(const Matrix& A, const Matrix& B)
{
	assert(A._cols == B._rows);
	make(A._rows, B._cols);

	for (int i=0; i<A._rows; i++)
		for (int j=0; j<B._cols; j++) {
			double sum = 0.0;
			for (int k=0; k<A._cols; k++)
				sum += A(i,k) * B(k,j);
			(*this)(i,j) = sum;
		}
} // multiply

/** multiplyT
 * multiply the transpose of A with B: A.T() * B
 */
void Matrix::multiplyT(const Matrix& A, const Matrix& B)
{
	assert(A._rows == B._rows);
	make(A._cols, B._cols);

	for (int i=0; i<A._cols; i++)
		for (int j=0; j<B._cols; j++) {
			double sum = 0.0;
			for (int k=0; k<A._rows; k++)
				sum += A(k,i) * B(k,j);
			(*this)(i,j) = sum;
		}
} // multiplyT

/** multiply
 * multiply matrix in place with B
 */
void Matrix::multiply(const Matrix& B)
{
	Matrix C(_rows, B._cols);

	int i,j,k;
	for (i=0; i<_rows; i++) {
		for (j=0; j<B._cols; j++) {
			double sum = 0.0;
			for (k=0; k<B._rows; k++)
				sum += (*this)(i,k) * B(k,j);
			C(i,j) = sum;
		}
	}

	_cols = B._cols;
	delete [] _data;
	_data = C._data;
	C._data = NULL;
} // multiply

/** negate matrix */
void Matrix::negate()
{
	for (int i=0; i<_rows; i++)
		for (int j=0; j<_cols; j++)
			(*this)(i,j) = -(*this)(i,j);
} // negate

/** det */
double Matrix::det() const
{
	assert (_rows == _cols);

	if (_rows == 2)
		return m(0,0)*m(1,1) - m(1,0)*m(0,1);
	else
	if (_rows == 3)
		return m(0,0)*(m(1,1)*m(2,2) - m(2,1)*m(1,2))
		     - m(0,1)*(m(1,0)*m(2,2) - m(2,0)*m(1,2))
		     + m(0,2)*(m(1,0)*m(2,1) - m(2,0)*m(1,1));
	else
		return 0.0;
} // det

/** return trace of matrix (sum of diagonal elements) */
double Matrix::trace() const
{
	assert(_rows==_cols);
	double t = 0.0;
	for (int i=0; i<_rows; i++)
		t += m(i,i);
	return t;
} // trace

/** clean up almost zero elements of matrix */
void Matrix::cleanupZero(const double eps)
{
	int siz = _rows * _cols;
	for (int i=0; i<siz; i++)
		if (Eq0(_data[i], eps))
			_data[i] = 0.0;
} // cleanupZero

/** return the inverse of a matrix. */
void Matrix::inverse()
{
	assert(_rows==_cols);
	int	i, j, indx[MAXSIZE];
	double	d, col[MAXSIZE];

	// make a copy of matrix
	Matrix mat(*this);
	mat.ludcmp(indx, &d);		/* matrix lu decomposition */

	for (j=0; j<_rows; j++) {	/* matrix inversion */
		for (i=0; i<_rows; i++)
			col[i] = 0.0;
		col[j] = 1.0;
		mat.lubksb(indx, col);
		for (i=0; i<_rows; i++)
			(*this)(i,j) = col[i];
	}
} // inverse

/**
 * lubksb(Mat4 a, int *indx, Flt b[])
 * backward substitution
 *	int	*indx	row permutation record
 *	double	b[]	right hand vector (?)
 */
void Matrix::lubksb(int *indx, double b[])
{
	int	i, j, ii=-1;
	double	sum;

	for (i=0; i<_rows; i++) {
		int ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii>=0)
			for (j=ii; j<=i-1; j++)
				sum -= m(i,j) * b[j];
		else
		if (sum != 0.0)
			ii = i;
		b[i] = sum;
	}
	for (i=_rows-1; i>=0; i--) {
		sum = b[i];
		for (j=i+1; j<_rows; j++)
			sum -= m(i,j) * b[j];
		b[i] = sum/m(i,i);
	}
} // lubksb

/**
 * ludcmp(Mat4 a, int *indx, Flt *d)
 * LU decomposition.
 * Parameteres
 *	Mat4	a		input matrix. gets thrashed
 *	int	*indx		row permutation record
 *	Flt	*d		+/- 1.0 even or odd # of row interchanges
 */
void Matrix::ludcmp(int *indx, double *d)
{
	double	vv[MAXSIZE];               /* implicit scale for each row */
	double	big, dum, sum, tmp;

	*d = 1.0;
	for (int i=0; i<_rows; i++) {
		big = 0.0;
		for (int j=0; j<_rows; j++) {
			if ((tmp=Abs(m(i,j))) > big)
				big = tmp;
		}
		if (big == 0.0) {
			//cerr << "error ludcmp(): singular matrix found..." << endl;
			throw 0;
		}
		vv[i] = 1.0/big;
	}
	for (int j=0; j<_rows; j++) {
		for (int i=0; i<j; i++) {
			sum = m(i,j);
			for (int k=0; k<i; k++)
				sum -= m(i,k) * m(k,j);
			(*this)(i,j) = sum;
		}
		big  = 0.0;
		int imax = 0;
		for (int i=j; i<_rows; i++) {
			sum = m(i,j);
			for (int k=0; k<j; k++)
				sum -= m(i,k)*m(k,j);
			(*this)(i,j) = sum;
			if ((dum=vv[i]*Abs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (int k=0; k<_rows; k++) {
				dum = m(imax,k);
				(*this)(imax,k) = m(j,k);
				(*this)(j,k) = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (m(j,j) == 0.0)
			(*this)(j,j) = 1.0E-20;      /* can be 0.0 also... */
		if (j != _rows-1) {
			dum = 1.0/m(j,j);
			for (int i=j+1; i<_rows; i++)
				(*this)(i,j) *= dum;
		}
	}
} // ludcmp

/** matrix print out */
ostream& operator << (ostream& s, const Matrix& matrix)
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

/** solve over determined system A*X = B
 * @param A	input matrix
 * @param B	input matrix
 * @param X	return solution
 * @return true on success
 */
bool solveOverDetermined(Matrix& A, Matrix& B, Matrix& X)
{
	Matrix ATA, ATB;
	ATA.multiplyT(A, A);
	ATB.multiplyT(A, B);
	ATA.cleanupZero(Vector::_epsilon);
	ATB.cleanupZero(Vector::_epsilon);
	try {
		ATA.inverse();
	} catch (...) {
		return false;
	}
	X.multiply(ATA,ATB);
	return true;
} // solveOverDetermined
