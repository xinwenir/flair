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

#include <math.h>
#include <assert.h>
#include <string.h>

#include <iomanip>
#include <iostream>

#include "os.h"
#include "bmath.h"
#include "matrix2.h"

using namespace std;

double Matrix2::_identity[4] = {
	1.0, 0.0,
	0.0, 1.0
};

/** transpose */
void Matrix2::transpose()
{
	double t = (*this)(0,1);
	(*this)(1,0) = (*this)(0,1);
	(*this)(0,1) = t;
} // transpose

/** operator * *
Matrix2& Matrix2::operator *= (const Matrix2& mat)
{
	// FIXME Optimize
	return (*this) = (*this)*mat;
} * operator *= */

/** negate matrix */
void Matrix2::negate()
{
	(*this)(0,0) = -(*this)(0,0);
	(*this)(0,1) = -(*this)(0,1);
	(*this)(1,0) = -(*this)(1,0);
	(*this)(1,1) = -(*this)(1,1);
} // negate

/** invert the matrix in place and
 * @return true on success, false otherwise
 */
bool Matrix2::inverse(const double eps)
{
	double D = det();
	if (Eq0(D, eps)) return false;
	D = 1.0/D;

	double a = (*this)(0,0);
	double b = (*this)(0,1);
	double c = (*this)(1,0);
	double d = (*this)(1,1);

	(*this)(0,0) =  d*D;
	(*this)(0,1) = -b*D;
	(*this)(1,0) = -c*D;
	(*this)(1,1) =  a*D;

	return true;
} // inverse

/** rotate matrix around an axis
 * @param angle
 * */
void Matrix2::rotate(const double angle)
{
	double s, c;
	bsincos(angle, &s, &c);

	(*this)(0,0) =  c;
	(*this)(0,1) = -s;
	(*this)(1,0) =  s;
	(*this)(1,1) =  c;
} // rotate

/** matrix print out */
ostream& operator << (ostream& s, const Matrix2& matrix)
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
