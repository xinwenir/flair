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
 *
 * Date:	16-Feb-2016
 */
#ifndef __LINE_H
#define __LINE_H

#include "point.h"

/** Calculate the shortest distance between a point and a line
 * @param point to find distance
 * @param lineStart start point of line
 * @param lineEnd end point of line
 * @param distance if not NULL return the distance
 * @return TRUE is the shortest distance is in the segment between the
 *	lineStart lineEnd
 */
bool pointLineDistance( const Point& point,
		const Point& start, const Point& end,
		double *distance=NULL);

/** @return the minimum distance between two lines */
double lineLineDistance(const Point& p1, const Vector& z1,
			const Point& p2, const Vector& z2);

/** Calculate the line segment PaPb that is the shortest route between
 * two lines P1P2 and P3P4. Calculate also the values of mua and mub where
 *    Pa = P1 + mua (P2 - P1)
 *    Pb = P3 + mub (P4 - P3)
 * @param p1,p2 points defining the first line segment
 * @param p3,p4 points defining the second line segment
 * @param pa,pb if not NULL return the position of intersection on both lines
 * @param mua,mub return the distance from p1,p3 respectively
 * @param eps the required accuracy of the operation
 * @return false if no solution exists.
 */
bool lineLineIntersect(
		const Point& p1, const Point& p2,
		const Point& p3, const Point& p4,
		Point *pa, Point *pb,
		double *mua,  double *mub,
		const double eps=Vector::_epsilon);
#endif
