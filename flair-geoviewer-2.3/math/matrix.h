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

#ifndef __MATRIX_H
#define __MATRIX_H

#include <ostream>

#include "os.h"
#include "vector.h"

////////////////////// Matrix ////////////////////////
/**
 * Matrix class
 */
class Matrix {
private:
	int	 _rows;		/** Matrix rows	*/
	int	 _cols;		/** Matrix columns	*/
	double	*_data;		/** allocated data	*/

public:
	Matrix() : _rows(0), _cols(0), _data(NULL) {};
	Matrix(int r, int c=0)	: _rows(0), _cols(0), _data(NULL) { make(r,c); };
	Matrix(const Matrix& mat) : _rows(0), _cols(0), _data(NULL) { copy(mat); }
	~Matrix()	{ if (_data) delete [] _data; };

	void	make(int r, int c=0);

	// Get/Set data
	int	rows()		const { return _rows; }
	int	columns()	const { return _cols; }

	double	operator () (int r, int c) const {
			assert(r<_rows && c<_cols);
			return _data[r*_cols+c];
		}
	double&	operator () (int r, int c) {
			assert(r<_rows && c<_cols);
			return _data[r*_cols+c];
		}
	double	m(int r, int c)	const {return (*this)(r,c);}

	void	copy(const Matrix& src);
	Matrix&	operator =(const Matrix& matrix) { copy(matrix); return *this;}
	void	identity();
	void	zero();
	void	transpose();
	void	transpose(const Matrix& A);
	void	T() {transpose(); }

	void	multiply(const Matrix& A, const Matrix& B);
	void	multiplyT(const Matrix& A, const Matrix& B);
	void	multiply(const Matrix& B);
	/** operator *= */
	Matrix&	operator *= (const Matrix& B)	{ multiply(B); return *this; }
	void	negate();

	double	det() const;
	void	inverse();
	void	inverse(const Matrix& src) {
			copy(src);
			inverse();
		}

	double	trace() const;

	void	cleanupZero(const double eps);
private:
	void	lubksb(int *indx, double b[]);
	void	ludcmp(int *indx, double *d);
}; // Matrix

bool solveOverDetermined(Matrix& A, Matrix& B, Matrix& X);

std::ostream& operator << (std::ostream&, const Matrix&);
#endif
