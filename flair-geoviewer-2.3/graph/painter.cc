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
 * Date:	16-Jul-2010
 */

#include <wchar.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <iostream>

#include "os.h"
#include "bmath.h"
#include "color.h"
#include "painter.h"

// if wchar_t is 4 bytes use the wmemset for faster background setting
#if WCHAR_MAX > 0xffffu
//	wchar_t is 4 bytes
#	define WCHAR4
#endif

#define MINSIZE		10
#define MAXSIZE		3000
#define MINSCALE	1e-7
#define MAXSCALE	1e14

using namespace std;

/* ============================= Painter ================================= */
/** Painter */
Painter::Painter(const int w, const int h) :
	_width(0),
	_height(0),
	_dataSize(0),
	_maxSize(0),
	_data(NULL)
{
	background(0x707070);
	init(w,h);
} // Painter

/** ~Painter */
Painter::~Painter()
{
	if (_data) delete [] _data;
} /* ~Painter */

/** init
 * @param w	width of image
 * @param h	height of image
 */
void Painter::init(const int w, const int h)
{
	_width  = Range(MINSIZE, w, MAXSIZE);
	_height = Range(MINSIZE, h, MAXSIZE);
	resetClip();

	_dataSize = _width*_height*(int)sizeof(*_data);
	if (_dataSize > _maxSize) {
		if (_width*_height < 1<<23) {
			// good for typical 1920x1280 screen
			// round to the next power of 2
			// to avoid too many resizings
			int ds = 1;
			while (_dataSize!=0) {
				_dataSize >>= 1;
				ds <<= 1;
			}
			_dataSize = ds;
		}
		_maxSize = _dataSize;
		if (_data) delete [] _data;
		_data = new dword[_maxSize];
	}
	clear();
} // init

/** reset clipping region */
void Painter::resetClip()
{
	_clip.left   = 0;
	_clip.top    = 0;
	_clip.right  = _width-1;
	_clip.bottom = _height-1;
} // resetClip

/** set clipping region. If values are negative then consider them from the other side
 * @param l	X-left coordinate
 * @param t	Y-top coordinate
 * @param r	X-right coordinate
 * @param b	Y-bottom coordinate
 */
void Painter::clip(int l, int t, int r, int b)
{
	if (l<0) l = _width  + l;
	if (t<0) t = _height + t;
	if (r<0) r = _width  + r;
	if (b<0) b = _height + b;
	_clip.left   = Range(0, l, _width-1);
	_clip.top    = Range(0, t, _height-1);
	_clip.right  = Range(0, r, _width-1);
	_clip.bottom = Range(0, b, _height-1);
} // clip

/** clip line coordinates to the appropriate clipping region
 * @param x1,y1	starting point of line
 * @param x2,y2	end point of line
 * @return true if a fraction of the line falls inside th clipping region
 */
bool Painter::clipLine(int *x1, int *y1, int *x2, int *y2) const
{
	union _OutCodeUnion ocu1, ocu2;
	bool   in;
	bool   out;

	/* Initialize 4-bit codes */
	ocu1.outcodes = 0;
	ocu1.ocs.code0 = (*x1 < _clip.left);
	ocu1.ocs.code1 = (*y1 < _clip.top);
	ocu1.ocs.code2 = (*x1 > _clip.right);
	ocu1.ocs.code3 = (*y1 > _clip.bottom);

	ocu2.outcodes = 0;
	ocu2.ocs.code0 = (*x2 < _clip.left);
	ocu2.ocs.code1 = (*y2 < _clip.top);
	ocu2.ocs.code2 = (*x2 > _clip.right);
	ocu2.ocs.code3 = (*y2 > _clip.bottom);

	in  = ((ocu1.outcodes | ocu2.outcodes) == 0);
	out = ((ocu1.outcodes & ocu2.outcodes) != 0);

	while (!out && !in) {
		if (ocu1.outcodes==0) {		/* Swap endpoints if necessary so */
			Swap(*x1,*x2);
			Swap(*y1,*y2);
			Swap(ocu1.outcodes,ocu2.outcodes);
		}

		if (ocu1.ocs.code0) {		/* Clip left */
			*y1 += (int)((long)((*y2-*y1) * (_clip.left-*x1))/(*x2-*x1));
			*x1 = _clip.left;
		}
		else
		if (ocu1.ocs.code1) {		/* Clip above */
			*x1 += (int)((long)((*x2-*x1) * (_clip.top-*y1))/(*y2-*y1));
			*y1 = _clip.top;
		}
		else
		if (ocu1.ocs.code2) {		/* Clip right */
			*y1 += (int)((long)((*y2-*y1) * (_clip.right-*x1))/(*x2-*x1));
			*x1 = _clip.right;
		}
		else
		if (ocu1.ocs.code3) {		/* Clip below */
			*x1 += (int)((long)((*x2-*x1) * (_clip.bottom-*y1))/(*y2-*y1));
			*y1 = _clip.bottom;
		}

		/* update for (*x1,*y1) */
		ocu1.outcodes = 0;
		ocu1.ocs.code0 = (*x1 < _clip.left);
		ocu1.ocs.code1 = (*y1 < _clip.top);
		ocu1.ocs.code2 = (*x1 > _clip.right);
		ocu1.ocs.code3 = (*y1 > _clip.bottom);

		in  = ((ocu1.outcodes | ocu2.outcodes) == 0); /* update */
		out = ((ocu1.outcodes & ocu2.outcodes) != 0); /* 4-bit codes */
	}
	return in;
}  /* clipLine */

/** clear image
 * @param byte	color of byte to fill 0x00 or 0xff
 */
void Painter::clear()
{
	if (_clip.left==0 && _clip.top==0 &&
	    _clip.right==_width-1 && _clip.bottom==_height-1) {
#ifdef WCHAR4
		wmemset((wchar_t*)_data, (wchar_t)_background, _dataSize>>2);
#else
		// dirty trick
		_data[0] = _background;
		memcpy(_data+1, _data, _dataSize-sizeof(dword));
#endif
	} else {
#ifdef WCHAR4
		int mx = (_clip.right - _clip.left);
#else
		int mx = (_clip.right - _clip.left) * (int)sizeof(dword);
#endif
		int my = _clip.bottom - _clip.top;
		dword *src = pixelPtr(_clip.left, _clip.top);
		for (int i=0; i<my; i++) {
#ifdef WCHAR4
			wmemset((wchar_t*)src, (wchar_t)_background, mx);
#else
			src[0] = _background;
			memcpy(src+1, src, mx-sizeof(dword));
#endif
			src += _width;
		}
	}
} // clear

/** move image by (dx,dy)
 * @param dx,dy	pixels to move image
 */
void Painter::move(const int dx, const int dy)
{
	dword *src, *dst;
	int mx = (_width - Abs(dx))*(int)sizeof(dword);
	int my = _height - Abs(dy);

	if (dy >= 0) {
		if (dx>=0) {
			src = pixelPtr( 0, _height-1-dy);
			dst = pixelPtr(dx, _height-1);
		} else {	/* dx<0 */
			src = pixelPtr(-dx, _height-1-dy);
			dst = pixelPtr(  0, _height-1);
		}
		for (int j=0; j<my; j++) {
			memmove(dst, src, mx);
			dst -= _width;
			src -= _width;
		}
	} else {
		if (dx>=0) {
			src = pixelPtr( 0,-dy);
			dst = pixelPtr(dx,  0);
		} else {	/* dx<0 */
			src = pixelPtr(-dx,-dy);
			dst = pixelPtr(  0,  0);
		}
		for (int j=0; j<my; j++) {
			memmove(dst, src, mx);
			dst += _width;
			src += _width;
		}
	}
} // move

/** unclippedLine
 * @param x1,y1	starting point of line
 * @param x2,y2	end point of line
 * @param color	color of line
 */
void Painter::unclippedLine(int x1, int y1, int x2, int y2, const dword color)
{
	int   d, dx, dy;
	int   Ainc, Binc;
	int   x, y;

	dx = Abs(x1-x2);
	dy = Abs(y1-y2);

	if (dx > dy) {
		if (x1 > x2) {		/* force x1<x2 */
			Swap(x1,x2);
			Swap(y1,y2);
		}
		int yinc = (y2>y1)? 1: -1;	/* determine increment for y */
		d = 2*dy - dx;
		Ainc = 2*(dy-dx);
		Binc = 2*dy;

		x = x1;
		y = y1;
		pixel(x,y,color);

		for (x=x1+1; x<=x2; x++) {
			if (d >= 0) {		/* setpixel A */
				y += yinc;
				d += Ainc;
			} else
				d += Binc;	/* setpixel B */
			pixel(x,y,color);
		}
	} else {	/* DO THE SAME BUT FOR Y */
		if (y1 > y2) {		/* force y1<y2 */
			Swap(y1,y2);
			Swap(x1,x2);
		}
		int xinc = (x2>x1)? 1: -1;		/* determine increment for x */
		d = 2*dx - dy;
		Ainc = 2*(dx-dy);
		Binc = 2*dx;

		x = x1;
		y = y1;
		pixel(x,y,color);

		for (y=y1+1; y<=y2; y++) {
			if (d >= 0) {		/* setpixel A */
				x += xinc;
				d += Ainc;
			} else
				d += Binc;	/* setpixel B */
			pixel(x,y,color);
		}
	}
} // unclippedLine

/** unclippedThickLine
 * @param x1,y1	starting point of line
 * @param x2,y2	end point of line
 * @param d	line thickness in pixels
 * @param color	color of line
 */
void Painter::unclippedThickLine(int x1, int y1, int x2, int y2, const int w, const dword color)
{
	int   d, dx, dy;
	int   Ainc, Binc;
	int   x, y;

	dx = Abs(x1-x2);
	dy = Abs(y1-y2);

	if (!dx && !dy) {
		// fill a square
		int w2 = w/2;
		for (y=y1-w2; y<=y1+w2; y++)
			for (x=x1-w2; x<=x1+w2; x++)
				pixel(x,y,color);
		return;
	}

	if (dx > dy) {
		if (x1 > x2) {		/* force x1<x2 */
			Swap(x1,x2);
			Swap(y1,y2);
		}
		int l = w * isqrt(Sqr(dx)+Sqr(dy)) / dx / 2;
		int yinc = (y2>y1)? 1: -1;	/* determine increment for y */
		d = 2*dy - dx;
		Ainc = 2*(dy-dx);
		Binc = 2*dy;

		x = x1;
		y = y1;
		for (int yy=y-l; yy<=y+l; yy++)
			pixel(x,yy,color);

		for (x=x1+1; x<=x2; x++) {
			if (d >= 0) {		/* setpixel A */
				y += yinc;
				d += Ainc;
			} else
				d += Binc;	/* setpixel B */
			for (int yy=y-l; yy<=y+l; yy++)
				pixel(x,yy,color);
		}
	} else {	/* DO THE SAME BUT FOR Y */
		if (y1 > y2) {		/* force y1<y2 */
			Swap(y1,y2);
			Swap(x1,x2);
		}
		int l = w * isqrt(Sqr(dx)+Sqr(dy)) / dy / 2;
		int xinc = (x2>x1)? 1: -1;		/* determine increment for x */
		d = 2*dx - dy;
		Ainc = 2*(dx-dy);
		Binc = 2*dx;

		x = x1;
		y = y1;
		for (int xx=x-l; xx<=x+l; xx++)
			pixel(xx,y,color);

		for (y=y1+1; y<=y2; y++) {
			if (d >= 0) {		/* setpixel A */
				x += xinc;
				d += Ainc;
			} else
				d += Binc;	/* setpixel B */
			for (int xx=x-l; xx<=x+l; xx++)
				pixel(xx,y,color);
		}
	}
} // unclippedThickLine

/** draw pixel at (x,y) with color using intensity from 0..1
 * @param x,y		position to draw
 * @param color		color to use
 * @param intensity	intensity of the color*intensity
 */
void Painter::pixel(const int x, const int y, dword color, double intensity)
{
	Color32 c; c.val = color;
	Color32 p; p.val = pixel(x,y);

	c.rgb.red   = (int)((double)c.rgb.red   * intensity + (double)p.rgb.red   * (1.0-intensity));
	c.rgb.green = (int)((double)c.rgb.green * intensity + (double)p.rgb.green * (1.0-intensity));
	c.rgb.blue  = (int)((double)c.rgb.blue  * intensity + (double)p.rgb.blue  * (1.0-intensity));
	pixel(x,y,c.val);
} // pixel

/** unClippedLineAntialias
 * @param x1,y1	starting point of line
 * @param x2,y2	end point of line
 * @param color	color of line
 *
 * FIXME DOESN'T WORK PROPERLY
 */
void Painter::unclippedLineAntialias(double x1, double y1, double x2, double y2, const dword color)
{
	double dx = x2 - x1;
	double dy = y2 - y1;
	if (dx==0.0 && dy==0.0) return;

	if (Abs(dx) < Abs(dy)) {
		Swap(x1, y1);
		Swap(x2, y2);
		Swap(dx, dy);
	}

	if (x2 < x1) {
		Swap(x1,x2);
		Swap(y1,y2);
	}

	double gradient = dy / dx;

	// handle first endpoint
	double xend = Round(x1);
	double yend = y1 + gradient * (xend - x1);
	double xgap = RFrac(x1 + 0.5);
	int xpxl1 = (int)xend;			// this will be used in the main loop
	int ypxl1 = Int(yend);
	pixel(xpxl1, ypxl1,   color, RFrac(yend)*xgap);
	pixel(xpxl1, ypxl1+1, color,  Frac(yend)*xgap);
	double intery = yend + gradient;	// first y-intersection for the main loop

	// handle second endpoint
	xend = Round(x2);
	yend = y2 + gradient * (xend - x2);
	xgap = Frac(x2 + 0.5);
	int xpxl2 = (int)xend;			// this will be used in the main loop
	int ypxl2 = Int(yend);
	pixel(xpxl2, ypxl2,   color, RFrac(yend)*xgap);
	pixel(xpxl2, ypxl2+1, color,  Frac(yend)*xgap);

	// main loop
	for (int x=xpxl1+1; x<=xpxl2-1; x++) {
		pixel(x, Int(intery),   color, RFrac(intery));
		pixel(x, Int(intery)+1, color,  Frac(intery));
		intery += gradient;
	}
} // unclippedLineAntialias

/** draw rectangle
 * @param x1,y1	upper-left corner
 * @param x2,y2 lower-right corner
 * @param color	color to draw
 */
bool Painter::rectangle(int x1, int y1, int x2, int y2, const dword color)
{
	bool b = line(x1,y1, x2,y1, color);
	b |= line(x2,y1, x2,y2, color);
	b |= line(x2,y2, x1,y2, color);
	b |= line(x1,y2, x1,y1, color);
	return b;
} // rectangle

/** fill rectangle
 * @param x1,y1	upper-left corner
 * @param x2,y2 lower-right corner
 * @param color	color to draw
 */
void Painter::fillRect(int x1, int y1, int x2, int y2, const dword color)
{
	clipPoint(&x1, &y1);
	clipPoint(&x2, &y2);

	dword *src = pixelPtr(x1,y1);
	for (int i=x1; i<=x2; i++)
		*src++ = color;

	int w = x2-x1+1;
	src -= w;
	dword *dst = src + _width;
	w *= (int)sizeof(*_data);
	for (int j=y1+1; j<=y2; j++) {
		memcpy(dst, src, w);
		dst += _width;
	}
} // fillRect

/** lighten or darken rectangle
 * @param x1,y1	upper-left corner
 * @param x2,y2 lower-right corner
 * @param level	to shift colors
 */
void Painter::levelShiftRect(int x1, int y1, int x2, int y2, int level)
{
	if (level==0) return;

	clipPoint(&x1, &y1);
	clipPoint(&x2, &y2);

	if (level > 0) {
		for (int j=y1; j<=y2; j++) {
			dword *ptr = pixelPtr(x1,j);
			for (int i=x1; i<=x2; i++, ptr++)
				*ptr = Lighten(*ptr, level);
		}
	} else {
		level = -level;
		for (int j=y1; j<=y2; j++) {
			dword *ptr = pixelPtr(x1,j);
			for (int i=x1; i<=x2; i++, ptr++)
				*ptr = Darken(*ptr, level);
		}
	}
} // levelShiftRect

/** draw bitmap
 * @param x,y	upper-left corner where to draw
 * @param w,h	width and height of the bitmap
 * @param bitmap	bitmap data of size w*h
 * @return	true if it draw something, else false
 */
bool Painter::drawBitmap(int x, int y, int w, int h, const dword *bitmap)
{
	int xs=x,   ys=y;
	int xe=x+w-1, ye=y+h-1;

	// Clip starting and ending point
	//
	clipPoint(&xs, &ys);
	clipPoint(&xe, &ye);

	if (xe-xs<1 || ye-ys<1) return false;

	const Color32 *src = (const Color32*)bitmap;
	dword *dst   = pixelPtr(xs,ys);

	// skip pixels outside
	src += (ys-y)*w + (xs-x);

	int dw = xe-xs+1;
	w -= dw;
	int Dw = _width-dw;

	//bool found = false;
	for (int j=ys; j<=ye; j++) {
		for (int i=xs; i<=xe; i++) {
			if (src->rgb.alpha != 0xFF)
				*dst++ = src->val;
			else
				dst++;
			src++;
		}
		src +=  w;
		dst += Dw;
	}
	return true;
} // drawBitmap

/** flood fill similar colors
 * @param x,y	starting point
 * @param color	color to fill
 */
void Painter::fill(int x, int y, const dword color, const dword color2, const FillType type)
{
#define STACKSIZE 500
	int stackX[STACKSIZE];
	int stackY[STACKSIZE];

	dword *pixelup, *pixeldown;

	if (!insideClip(x,y)) return;
	dword startColor = pixel(x,y);
	if (startColor==color || (type && startColor==color2)) return;

	int tos = 0;					// top of stack
	for (;;) {					// do forever
		while (pixel(x,y)!=startColor) {	// pick an empty point
			if (!tos) return;
			x = stackX[--tos];
			y = stackY[  tos];
		}

		int xrem = x+1;
		dword* pix = pixelPtr(x,y);			// get pointer

		/* scan line */
		if (_clip.top<y)
			pixelup = pix - _width;		// up one line
		else
			pixelup = NULL;

		if (y<_clip.bottom)
			pixeldown = pix + _width;	// down one line
		else
			pixeldown = NULL;

		int startup   = (pixelup && *pixelup==startColor);
		int startdown = (pixeldown && *pixeldown==startColor);

		int up   = 0;
		int down = 0;
		while (*pix == startColor) {
			switch (type) {
				case FILL_FLOOD:
					*pix = color;
					break;
				case FILL_DOTS:
					*pix = (((x+y)&3) || (y&1))? color : color2;
					break;
				case FILL_HASH:
					*pix = (x+y)&7? color : color2;
					break;
				case FILL_HASHR:
					*pix = (x-y)&7? color : color2;
					break;
				case FILL_X:
					*pix = (((x+y)&3) && (x-y)&3)? color : color2;
					break;
				case FILL_X2:
					*pix = (((x+y)&7) && (x-y)&7)? color : color2;
					break;
				default:
					assert(0);
			}

			int pup = (pixelup && *pixelup==startColor);
			if (!up && pup && tos<STACKSIZE) {
				stackX[tos  ] = x;
				stackY[tos++] = y-1;
			}
			up = pup;

			int pdown = (pixeldown && *pixeldown==startColor);
			if (!down && pdown && tos<STACKSIZE) {
				stackX[tos  ] = x;
				stackY[tos++] = y+1;
			}
			down = pdown;
			if (x==_clip.left) break;
			x--;
			pix--;
			if (pixelup) pixelup--;
			if (pixeldown) pixeldown--;
		}

		if (xrem>_clip.right) continue;
		x = xrem;

		pix = pixelPtr(x,y);

		/* scan line */
		if (_clip.top<y)
			pixelup = pix - _width;		// up one line
		else
			pixelup = NULL;

		if (y<_clip.bottom)
			pixeldown = pix + _width;	// down one line
		else
			pixeldown = NULL;

		up   = startup;
		down = startdown;
		while (*pix == startColor) {
			switch (type) {
				case FILL_FLOOD:
					*pix = color;
					break;
				case FILL_DOTS:
					*pix = (((x+y)&3) || (y&1))? color : color2;
					break;
				case FILL_HASH:
					*pix = (x+y)&7? color : color2;
					break;
				case FILL_HASHR:
					*pix = (x-y)&7? color : color2;
					break;
				case FILL_X:
					*pix = (((x+y)&3) && (x-y)&3)? color : color2;
					break;
				case FILL_X2:
					*pix = (((x+y)&7) && (x-y)&7)? color : color2;
					break;
				default:
					assert(0);
			}

			int pup = (pixelup && *pixelup==startColor);
			if (!up && pup && tos<STACKSIZE) {
				stackX[tos  ] = x;
				stackY[tos++] = y-1;
			}
			up = pup;

			int pdown = (pixeldown && *pixeldown==startColor);
			if (!down && pdown && tos<STACKSIZE) {
				stackX[tos  ] = x;
				stackY[tos++] = y+1;
			}
			down = pdown;
			if (x==_clip.right) break;
			x++;
			pix++;
			if (pixelup) pixelup++;
			if (pixeldown) pixeldown++;
		}
	}
} // fill

/** measure string */
int Painter::measuref(BFont& font, const char *fmt, ...)
{
	char	text[1024];
	va_list	ap;
	if (fmt==NULL) return 0;
	va_start(ap, fmt);
	vsprintf(text, fmt, ap);
	va_end(ap);

	return measure(font, text);
} // measuref

/** print formatted string
 * @return right most pixel of the string
 */
int Painter::printf(BFont& font, int x, int y, const dword color, const char *fmt, ...)
{
	char	text[1024];
	va_list	ap;
	if (fmt==NULL) return x;
	va_start(ap, fmt);
	vsprintf(text, fmt, ap);
	va_end(ap);

	for (uint8_t *ch=(uint8_t*)text; *ch; ch++)
		x += font.draw(*this, x, y, color, *ch);

	return x;
} // printf

/** vertical print formatted string
 * @return right most pixel of the string
 */
int Painter::vprintf(BFont& font, int x, int y, const dword color, const char *fmt, ...)
{
	char	text[1024];
	va_list	ap;
	if (fmt==NULL) return x;
	va_start(ap, fmt);
	vsprintf(text, fmt, ap);
	va_end(ap);

	for (uint8_t *ch=(uint8_t*)text; *ch; ch++)
		y -= font.vdraw(*this, x, y, color, *ch);

	return y;
} // printf

/** print formatted string with shadow
 * @return right most pixel of the string
 */
int Painter::printfShadow(BFont& font, int x, int y, int ofs,
			const dword color, const dword shadow, const char *fmt, ...)
{
	char	text[1024];
	va_list	ap;
	if (fmt==NULL) return x;
	va_start(ap, fmt);
	vsprintf(text, fmt, ap);
	va_end(ap);

	for (uint8_t *ch=(uint8_t*)text; *ch; ch++) {
		     font.draw(*this, x+ofs, y+ofs, shadow, *ch);
		x += font.draw(*this, x, y, color, *ch);
	}

	return x;
} // printfShadow

/** print formatted string with shadow
 * @return right most pixel of the string
 */
int Painter::printfOutline(BFont& font, int x, int y,
			const dword color, const dword outline, const char *fmt, ...)
{
	char	text[1024];
	va_list	ap;
	if (fmt==NULL) return x;
	va_start(ap, fmt);
	vsprintf(text, fmt, ap);
	va_end(ap);

	for (uint8_t *ch=(uint8_t*)text; *ch; ch++)
		x += font.drawOutline(*this, x, y, color, outline, *ch);

	return x;
} // printfShadow

/** operator << *
ostream& operator << (ostream& s, const Painter& painter)
{
	s << "Painter" << endl;
	s << "\tWindow  [" << painter.width << ", " << painter.height << "]" << endl;
	s << "\tMemory  dataSize=" << painter.dataSize << " Max=" << painter.maxSize << endl;
	s << "\tClip    (" << painter.clip().left << ", " << painter.clip().top << ") - ("
			   << painter.clip().right << ", " << painter.clip().bottom << ")" << endl;
	s << "\tBackground = " << painter.background << endl;
	return s;
} * operator << */

/* -------------------- BFont ----------------------- */
BFont::~BFont()
{
	clean();
} /* ~BFont */

/** clean */
void BFont::clean()
{
	_name.clear();
	if (imageData) {
		delete [] imageData;
		imageData = NULL;
	}
} // clean

/** load a font image
 * Loads A TGA File Into Memory
 * Expects to be compression:none, TopLeft format (byte 17=0x20)
 * @return true on success
 */
bool BFont::load(const char *filename)
{
#if 0
	// Uncompressed TGA Header
	typedef struct {
	0	uint8_t  identSize;	// size of ID field that follows 18 byte header (0 usually)
	1	uint8_t  colorMap;		// type of color map 0=none, 1=has palette
	2	uint8_t  imageType;	// type of image
					// 0=none,1=indexed,2=rgb,3=grey,+8=rle packed

	3	short colorMapStart;	// first color map entry in palette
	5	short colorMapLength;	// number of colors in palette
	7	uint8_t  colorMapBits;	// number of bits per palette entry 15,16,24,32

	// Alignemnt issues on 64bit for the xstart is in odd position!!!
	8	short xstart;		// image x origin
	10	short ystart;		// image y origin
	12	short width;		// image width in pixels
	14	short height;		// image height in pixels
	16	uint8_t  bits;		// image bits per pixel 8,16,24,32
	17	uint8_t  descriptor;	// image descriptor bits (vh flip bits)
	} TGA_Header;
	TGA_Header tga;
#endif
	uint8_t header[18];

	clean();
	assert(sizeof(short)==2);
	FILE *file = fopen(filename, "rb");
	if (file==NULL) return false;

	if (fread(&header, 1, sizeof(header), file) != sizeof(header)) {
		fclose(file);
		return false;
	}

	// It has to be grey, 8 bit
	if (header[1]  != 0 ||	// colourMap
	    header[2]  != 3 ||	// imageType
	    header[3]  != 0 ||	// colorMap start
	    header[4]  != 0 ||	//  -//-
	    header[5]  != 0 ||	// map length
	    header[6]  != 0 ||	//  -//-
	    header[7]  != 0 ||	// map bits
	    header[16] != 8) {	// xstart (should be 0)
		fclose(file);
		return false;
	}

	if (header[0]) fseek(file, header[0], SEEK_CUR);	// skip title

	imageWidth  = (int)header[12] + ((int)header[13]<<8);
	imageHeight = (int)header[14] + ((int)header[15]<<8);

	_width  = imageWidth  / 16;
	_height = imageHeight / 16;

	// Calculate The Memory Required For The TGA Data
	size_t imageSize = imageWidth * imageHeight * header[16]/8;

	// Reserve Memory To Hold The TGA Data
	imageData= new uint8_t[imageSize];

	if (imageData==NULL || fread(imageData, 1, imageSize, file)!=imageSize) {
		if(imageData!=NULL)
			delete [] imageData;
		imageData = NULL;
		fclose(file);
		return false;
	}
	fclose(file);

	// Copy name
	_name = filename;

	build();

	return true;
} // load

/** set */
void BFont::set(const char* n, int w, int h, uint8_t* data)
{
	clean();
	_name   = n;
	_width  = w;
	_height = h;
	imageData = new uint8_t[_width*_height];
	memcpy(imageData, data, _width*_height);
	build();
} // set

/** set */
void BFont::set(const char* n, int w, int h, dword* data)
{
	clean();
	_name   = n;
	_width  = w;
	_height = h;

	imageData = new uint8_t[_width*_height];

	uint8_t *bptr = imageData;
	dword *dptr = data;

	FILE *f = fopen("font.gray","wb");
	for (int i=0; i<_width*_height; i++, bptr++, dptr++) {
		*bptr = *dptr? 0xFF : 0x00;
		fputc(*bptr,f);
	}
	fclose(f);
	build();
} // set

/** build font */
void BFont::build()
{
	if (imageData==NULL) return;

	// Convert to black/white
	uint8_t *ptr = imageData;
	for (size_t i=0; i<(size_t)(_width*_height); i++)
		if (*ptr < 0x7F)
			*ptr = 0;
		else
			*ptr = 0xFF;

	// scan fonts to find the length of each character
	for (int ch=0; ch<256; ch++) {
		ptr = charPtr(ch);
		int maxlen = 0;
		for (int j=0; j<_height; j++) {
			for (int i=_width-2; i>=2; i--)
				if (ptr[i]) {
					maxlen = Max(maxlen, i);
					break;
			}
			ptr += imageWidth;
		}
		charLength[ch] = (uint8_t)(maxlen+2);
		// space
		if (ch==' ') charLength[ch] = (uint8_t)(_width/2);
	}
} // build

/** @return string length in pixels */
int BFont::measure(const char *str)
{
	int size = 0;
	for (const uint8_t *ch=(const uint8_t*)str; *ch; ch++)
		size += charLength[*ch];

	return size;
} // measure

/** draw a character into a viewport
 * @param painter	painter to use
 * @param x,y	location where to draw character
 * @param color	color to use
 * @param ch	character to draw
 */
int BFont::draw(Painter& painter, int x, int y, const dword color, const uint8_t ch)
{
	uint8_t  *src;
	dword *dst;

	if (imageData==NULL) return 0;

	int w = charLength[ch];

	src = charPtr(ch);
	dst = painter.pixelPtr(x,y);
	for (int j=0; j<_height; j++) {
		int xp = x;
		for (int i=0; i<w; i++) {
			if (painter.insideClip(xp,y) && *src != 0) *dst = color;
			xp++;
			src++;
			dst++;
		}
		y++;
		src += imageWidth - w;
		dst += painter.width() - w;
	}
	return w;
} // draw

/** draw a vertical character into a viewport
 * @param painter	painter to use
 * @param x,y	location where to draw character
 * @param color	color to use
 * @param ch	character to draw
 */
int BFont::vdraw(Painter& painter, int x, int y, const dword color, const uint8_t ch)
{
	uint8_t  *src;
	dword *dst;

	if (imageData==NULL) return 0;

	int w = charLength[ch];

	src = charPtr(ch);
	dst = painter.pixelPtr(x,y);
	for (int j=0; j<_height; j++) {
		int yp = y;
		for (int i=0; i<w; i++) {
			if (painter.insideClip(x,yp) && *src != 0) *dst = color;
			yp--;
			src++;
			dst -= painter.width();
		}
		x++;
		src += imageWidth - w;
		dst = painter.pixelPtr(x,y);
	}
	return w;
} // vdraw

/** draw a character outline
 * @param painter	painter to use
 * @param x,y	location where to draw character
 * @param color	color to use
 * @param outline	color to use
 * @param ch	character to draw
 */
int BFont::drawOutline(Painter& painter, int x, int y,
			const dword color, const dword outline, const uint8_t ch)
{
	uint8_t  *src;
	dword *dst;

	if (imageData==NULL) return 0;

	int w = charLength[ch];

	src = charPtr(ch);
	dst = painter.pixelPtr(x,y);
	int yp = y;
	for (int j=0; j<_height; j++) {
		int xp = x;
		for (int i=0; i<w; i++) {
			if (painter.insideClip(xp,yp)) {
				if (*src)
					*dst = color;
				else
				if (src[1] || src[-1])
					*dst = outline;
			}
			xp++;
			src++;
			dst++;
		}
		yp++;
		src += imageWidth - w;
		dst += painter.width() - w;
	}
	return w;
} // drawOutline
