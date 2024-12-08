/*
 *
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
 * Date:	13-Aug-2013
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
//#include <string.h>

#include "geo.h"
#include "xdraw.h"
#include "matrix4.h"
#include "painter.h"
#include "viewport.h"

using namespace std;

//static const char dottedPattern[2] = {1,1};
//static const char dashedPattern[2] = {3,3};

/* ============================== static functions =========================== */
bool clipSegment(int *x1, int *y1, int *x2, int *y2,
		const int left, const int top, const int right, const int bottom)
{
	union	_OutCodeUnion ocu1, ocu2;
	bool	in;
	bool	out;

	/* Initialize 4-bit codes */
	ocu1.outcodes = 0;
	ocu1.ocs.code0 = (*x1 < left);
	ocu1.ocs.code1 = (*y1 < top);
	ocu1.ocs.code2 = (*x1 > right);
	ocu1.ocs.code3 = (*y1 > bottom);

	ocu2.outcodes = 0;
	ocu2.ocs.code0 = (*x2 < left);
	ocu2.ocs.code1 = (*y2 < top);
	ocu2.ocs.code2 = (*x2 > right);
	ocu2.ocs.code3 = (*y2 > bottom);

	in  = ((ocu1.outcodes | ocu2.outcodes) == 0);
	out = ((ocu1.outcodes & ocu2.outcodes) != 0);

	while (!out && !in) {
		if (ocu1.outcodes==0) {		/* Swap endpoints if necessary so */
			Swap(*x1,*x2);
			Swap(*y1,*y2);
			Swap(ocu1.outcodes,ocu2.outcodes);
		}

		if (ocu1.ocs.code0) {		/* Clip left */
			*y1 += (int)((long)((*y2-*y1) * (left-*x1))/(*x2-*x1));
			*x1 = left;
		}
		else
		if (ocu1.ocs.code1) {		/* Clip above */
			*x1 += (int)((long)((*x2-*x1) * (top-*y1))/(*y2-*y1));
			*y1 = top;
		}
		else
		if (ocu1.ocs.code2) {		/* Clip right */
			*y1 += (int)((long)((*y2-*y1) * (right-*x1))/(*x2-*x1));
			*x1 = right;
		}
		else
		if (ocu1.ocs.code3) {		/* Clip below */
			*x1 += (int)((long)((*x2-*x1) * (bottom-*y1))/(*y2-*y1));
			*y1 = bottom;
		}

		/* update for (*x1,*y1) */
		ocu1.outcodes = 0;
		ocu1.ocs.code0 = (*x1 < left);
		ocu1.ocs.code1 = (*y1 < top);
		ocu1.ocs.code2 = (*x1 > right);
		ocu1.ocs.code3 = (*y1 > bottom);

		in  = ((ocu1.outcodes | ocu2.outcodes) == 0); /* update */
		out = ((ocu1.outcodes & ocu2.outcodes) != 0); /* 4-bit codes */
	}
	return in;
} // clipSegment

/** XDrawLine3D
 * draw a pseudo 3D line
 */
void XDrawLine3D(Display *display, Drawable drawable, GC gc,
		const ViewPort& view,
		const Point& a, const Point& b)
{
	XGCValues	gcValues;	// context values

	int xx1 = view.u2i(a.x);
	int yy1 = view.v2j(a.y);
	int xx2 = view.u2i(b.x);
	int yy2 = view.v2j(b.y);

	// It should take care of the perspective viewports
	gcValues.line_style = LineSolid;
	XChangeGC(display, gc, GCLineStyle, &gcValues);
	XDrawLine(display, drawable, gc, xx1, yy1, xx2, yy2);
#if 0
	if (Eq0(z1,SMALL) && Eq0(z2,SMALL)) {
		gcValues.line_style = LineSolid;
		gcValues.line_width = 1;
		XChangeGC(display, gc, GCLineWidth|GCLineStyle, &gcValues);
		XDrawLine(display, drawable, gc, xx1, yy1, xx2, yy2);
	} else
	if (z1>=-SMALL && z2>=-SMALL) {
		gcValues.line_style = LineSolid;
		gcValues.line_width = 0;
		XChangeGC(display, gc, GCLineWidth|GCLineStyle, &gcValues);
		XDrawLine(display, drawable, gc, xx1, yy1, xx2, yy2);
	} else
	if (z1<=SMALL && z2<=SMALL) {
		gcValues.line_style = LineOnOffDash;
		gcValues.line_width = 0;
		XChangeGC(display, gc, GCLineWidth|GCLineStyle, &gcValues);
		XSetDashes(display, gc, 0, dottedPattern, SIZE(dottedPattern));
		XDrawLine(display, drawable, gc, xx1, yy1, xx2, yy2);
	} else {
		double ww = z1 / (z2-z1);
		int xx0 = view.u2i(x1 - (x2-x1) * ww);
		int yy0 = view.v2j(y1 - (y2-y1) * ww);

		gcValues.line_style = LineSolid;
		gcValues.line_width = 0;
		XChangeGC(display, gc, GCLineWidth|GCLineStyle, &gcValues);

		if (z1>=-SMALL) {
			XDrawLine(display, drawable, gc, xx1, yy1, xx0, yy0);

			gcValues.line_style = LineOnOffDash;
			XChangeGC(display, gc, GCLineStyle, &gcValues);
			XSetDashes(display, gc, 0, dottedPattern, SIZE(dottedPattern));

			XDrawLine(display, drawable, gc, xx0, yy0, xx2, yy2);
		} else {
			XDrawLine(display, drawable, gc, xx2, yy2, xx0, yy0);

			gcValues.line_style = LineOnOffDash;
			XChangeGC(display, gc, GCLineStyle, &gcValues);
			XSetDashes(display, gc, 0, dottedPattern, SIZE(dottedPattern));

			XDrawLine(display, drawable, gc, xx0, yy0, xx1, yy1);
		}
	}
#endif
} // XDrawLine3D

/** draw a 3D axes system */
void XDrawAxes(Display *display, Drawable drawable, GC gc,
		const int x, const int y, const int size,
		const Matrix4& matrix, bool xyz)
{
	/* X - axis */
	double r = (double)size;
	int xe = x + Round(r * matrix(0,0));
	int ye = y - Round(r * matrix(0,1));
	int color;
	if (xyz)
		color = matrix(0,2)>=0.0? COLOR_X : COLOR_X_DARK;
	else
		color = matrix(0,2)>=0.0? COLOR_U : COLOR_U_DARK;

	XSetForeground(display, gc, color);
	XDrawLine(display, drawable, gc, x, y, xe, ye);

	/* Y - axis */
	xe = x + Round(r * matrix(1,0));
	ye = y - Round(r * matrix(1,1));
	if (xyz)
		color = matrix(1,2)>=0.0? COLOR_Y : COLOR_Y_DARK;
	else
		color = matrix(1,2)>=0.0? COLOR_V : COLOR_V_DARK;
	XSetForeground(display, gc, color);
	XDrawLine(display, drawable, gc, x, y, xe, ye);

	/* Z - axis */
	xe = x + Round(r * matrix(2,0));
	ye = y - Round(r * matrix(2,1));
	if (xyz)
		color = matrix(2,2)>=0.0? COLOR_Z : COLOR_Z_DARK;
	else
		color = matrix(2,2)>=0.0? COLOR_W : COLOR_W_DARK;
	XSetForeground(display, gc, color);
	XDrawLine(display, drawable, gc, x, y, xe, ye);
} // XDrawAxes

/** draw a rotated ellipse */
void XDrawEllipse(Display *display, Drawable drawable, GC gc,
		const int x, const int y,
		const int major, const int minor,
		double rotation , double from, double to)
{
#define SEGMENTS	32
	XPoint pts[SEGMENTS];
	double step = (from-to) / (double)(SEGMENTS-1);
	double crot, srot;
	sincos(rotation, &srot, &crot);

	double angle = from;
	for (int i=0; i<SEGMENTS; i++, angle += step) {
		double s, c;
		sincos(angle, &s, &c);
		double xx = (double)major * c;
		double yy = (double)minor * s;
		pts[i].x = x + (int)( xx*crot + yy*srot);
		pts[i].y = y + (int)(-xx*srot + yy*crot);
	}
	XDrawLines(display, drawable, gc, pts, SEGMENTS, CoordModeOrigin);
} // XDrawEllipse

/** draw a 3D track ball */
void XDrawTrackball(Display *display, Drawable drawable, GC gc,
		const int x, const int y, const int radius,
		const Matrix4& matrix)
{
	dword col;
//	col=0x7F0000;
//	for (int i=0; i<3; i++, col >>= 8) {
//		XSetForeground(display, gc, col);
//		double angle = atan2(matrix(i,1), matrix(i,0)) + PIover2;
//		XDrawEllipse(display, drawable, gc, x, y, radius, radius*matrix(i,2), angle, PI, PI2);
//	}
	col = 0xFF0000;
	for (int i=0; i<3; i++, col >>= 8) {
		XSetForeground(display, gc, col);
		double angle = atan2(matrix(i,1), matrix(i,0)) + PIover2;
		XDrawEllipse(display, drawable, gc, x, y, radius, Round(radius*matrix(i,2)), angle, 0.0, PI);
	}

	XSetForeground(display, gc, 0xAF00AF);

	XDrawArc(display, drawable, gc,
		x-radius, y-radius, 2*radius, 2*radius, 0, 360*64);

	XDrawAxes(display, drawable, gc, x, y, radius, matrix);
} // XDrawTrackball
