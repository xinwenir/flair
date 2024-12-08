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
 * Date:	05-Apr-2011
 *
 * TODO
 * GObject has to be split in two GObject and VObject
 * otherwise we will have problems when parallel drawings are done
 */

#ifndef __XDRAW_H
#define __XDRAW_H

#include <X11/X.h>
#include <X11/Xutil.h>

class ViewPort;
class Matrix4;
class Point;

/* General function prototypes */
bool clipSegment(int *x1, int *y1, int *x2, int *y2,
		const int left, const int top, const int right, const int bottom);
void XDrawLine3D(Display* display, Drawable drawable, GC gc,
		const ViewPort& view,
		const Point& a, const Point& b);
void XDrawAxes(Display* display, Drawable drawable, GC gc,
		const int x, const int y, const int size, const Matrix4& matrix, bool xyz=true);
void XDrawTrackball(Display* display, Drawable drawable, GC gc,
		const int x, const int y, const int radius, const Matrix4& matrix);
#endif
