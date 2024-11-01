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
 * author:	Vasilis.Vlachoudis@cern.ch
 * Date:	06-Feb-2013
 */

#include <iostream>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "spline.h"

using namespace std;

/** CardinalSpline */
CardinalSpline::CardinalSpline(const double a)
{
	calcMatrix(a);
} // CardinalSpline

/** calcMatrix */
void CardinalSpline::calcMatrix(const double a)
{
	A = a;
	matrix[0][0] =    -A; matrix[0][1] = 2.0-A; matrix[0][2] = A-2.0;     matrix[0][3] =    A;
	matrix[1][0] = 2.0*A; matrix[1][1] = A-3.0; matrix[1][2] = 3.0-2.0*A; matrix[1][3] =   -A;
	matrix[2][0] =    -A; matrix[2][1] =   0.0; matrix[2][2] = A;         matrix[2][3] =  0.0;
	matrix[3][0] =   0.0; matrix[3][1] =   1.0; matrix[3][2] = 0.0;       matrix[3][3] =  0.0;
} // calcMatrix

/** spline */
double CardinalSpline::spline(int k, const double t, BaseSplineNode* result)
{
	double T[4], R[4];

	T[0] = t*t*t;    T[1] = t*t;    T[2] = t;   T[3] = 1.0;

	for (int i=0; i<4; i++) {
		R[i] = 0.0;
		for (int j=0; j<4; j++)
			R[i] += T[j] * matrix[j][i];
	}
	double x =  0.0;

	for (int i=0; i<4; i++) {
		int idx = Range(0, k+i-1, pts.size()-1);
		x += R[i] * pts[idx]->x();
	}
	result->x(x);

	for (int j=0; j<result->size(); j++) {
		double y = 0.0;
		for (int i=0; i<4; i++) {
			int idx = Range(0, k+i-1, pts.size()-1);
			y += R[i] * pts[idx]->y(j);
		}
		result->y(j, y);
	}
	return x;
} // spline
