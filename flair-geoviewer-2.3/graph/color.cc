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

// if wchar_t is 4 bytes use the wmemset for faster background setting
#if WCHAR_MAX > 0xffffu
//	wchar_t is 4 bytes
#	define WCHAR4
#endif

#define MINSIZE		10
#define MAXSIZE		3000
#define MINSCALE	1e-7
#define MAXSCALE	1e14

// --- Declare Basic Colors ---
const Color Color::Black(0.0,0.0,0.0);
const Color Color::Red(1.0,0.0,0.0);
const Color Color::Green(0.0,1.0,0.0);
const Color Color::Blue(0.0,0.0,1.0);
const Color Color::Magenta(1.0,0.0,1.0);
const Color Color::Turquoise(0.0,1.0,1.0);
const Color Color::Yellow(1.0,1.0,0.0);
const Color Color::White(1.0,1.0,1.0);


using namespace std;

/** adjust color level inside the range lw..hg */
dword ColorLevel(const dword col, const Color32& lw, const Color32& hg)
{
	Color32 color;
	color.val = col;

	color.rgb.red   = (uint8_t)((((int)color.rgb.red  *((int)hg.rgb.red   - (int)lw.rgb.red)  )>>8)
			+ (int)lw.rgb.red);
	color.rgb.green = (uint8_t)((((int)color.rgb.green*((int)hg.rgb.green - (int)lw.rgb.green))>>8)
			+ (int)lw.rgb.green);
	color.rgb.blue  = (uint8_t)((((int)color.rgb.blue *((int)hg.rgb.blue  - (int)lw.rgb.blue) )>>8)
			+ (int)lw.rgb.blue);

	return color.val;
} // ColorLevel

/** alpha color blending */
dword alphaBlend(const Color32 &colA, const Color32 &colB)
{
	Color32 result;
	uint8_t alphaA      = colA.rgb.alpha;
	uint8_t alphaBnegA  = byteMul(colB.rgb.alpha, 255 - alphaA);

	result.rgb.red   = byteMul(colA.rgb.red,   alphaA) + byteMul(colB.rgb.red,   alphaBnegA);
	result.rgb.green = byteMul(colA.rgb.green, alphaA) + byteMul(colB.rgb.green, alphaBnegA);
	result.rgb.blue  = byteMul(colA.rgb.blue,  alphaA) + byteMul(colB.rgb.blue,  alphaBnegA);

	result.rgb.alpha = 0;

	return result.val;
} // alphaBlend

/* ========================= HSV ================================== *

			      ^ Value
		       _______|______
		       \      |      /   Side face of hexcone HSV
			 \    |    /
			   \  |  /
			     \|/
		    Green          Yellow               /
			-------------                 /
		      / \           / \             /
		    /     \       /     \         /  ^  Hue
		  /         \   /         \           \
	  Cyan  /_____________x_____________\   Red  --------->
		\       White  \            /                Saturation
		  \         /    \        /
		    \     /        \    /
		      \ /            \/          Top face of hexcone HSV
		  Blue  ~~~~~~~~~~~~~~ Magenta    color model

 * ================================================================ */

/** rgb2hsv - convert an RGB color to HSV
 *  Output:  h in [0,360],  s and v in [0,1], except if s = 0
 *           then h is undefined which is defined as a constant -1
 *           outside the interval [0,360]
 * @param	col	Color32 to convert to HSV
 * @param	h	return hue [0..360]
 * @param	s	return saturation [0..1]
 * @param	v	return value [0..1]
 */
void rgb2hsv(const Color32& col, double *h, double *s, double *v)
{
	int maximum = Max(col.rgb.red, col.rgb.green, col.rgb.blue);
	int minimum = Min(col.rgb.red, col.rgb.green, col.rgb.blue);

	/* value */
	*v = (double)maximum/255.0;

	if (maximum != 0)
		/* saturation */
		*s = (double)(maximum - minimum) / (double)maximum;
	else
		*s = 0.0;

	if (*s == 0.0)
		*h = -1.0;		/* hue, undefined */
	else {
		/* rc measures "distance" of color from red */
		double rc = (double)(maximum - col.rgb.red)   / (double)(maximum - minimum);
		double gc = (double)(maximum - col.rgb.green) / (double)(maximum - minimum);
		double bc = (double)(maximum - col.rgb.blue)  / (double)(maximum - minimum);

		if (col.rgb.red == maximum)
			*h = bc - gc;		/* resulting color between yellow and magenta */
		else
		if (col.rgb.green == maximum)
			*h = 2.0 + rc - bc;	/* resulting color between cyan and yellow */
		else
			*h = 4.0 + gc - rc;	/* resulting color between magenta and cyan */

		*h *= 60.0;			/* convert to degrees */

		if (*h < 0.0) *h += 360.0;	/* make non negative */
	}  /* chromatic case */
} // rgb2hsv

/** hsv2rgb - convert an HSV color to RGB
 * @param	h in [0..360]
 * @param	s in [0..1]
 * @param	v in [0..1], except if s = 0.0
 * @return	rgb value of color in dword
 */
dword hsv2rgb(double h, double s, double v)
{
	double r, g, b;

	if (s==0.0) {			/* achromatic color there is no hue */
		if (h==-1.0)		/* undefined */
			r = g = b = v;	/* this is the achromatic case */
		 else
			r = g = b = 0;	/* error if s = 0 and h has a value */
	} else {			/* chromatic color: there is a hue */
		if (h==360.0) h = 0.0;
		h /= 60.0;		/* h is now in [0,6] */

		int i = (int)h;		/* largest integer <= h */
		double f = h - (double)i;/* fractional part of h */
		double p = v * (1.0 - s);
		double q = v * (1.0 - s*f);
		double t = v * (1.0 - s*(1.0 - f));

		switch (i) {
			case  0:
				r = v;
				g = t;
				b = p;
				break;

			case  1:
				r = q;
				g = v;
				b = p;
				break;

			case  2:
				r = p;
				g = v;
				b = t;
				break;

			case  3:
				r = p;
				g = q;
				b = v;
				break;

			case  4:
				r = t;
				g = p;
				b = v;
				break;

			case  5:
				r = v;
				g = p;
				b = q;
				break;

			default:
				r = g = b = 0;
		}  /* end of switch */
	}

	return RGB(int(r*255.0), int(g*255.0), int(b*255.0));
} // hsv2rgb

/* =============================== Color ================================== */
/* ----------------------------------------------------------------- *
 *  setHSV(h,s,v);                                                 *
 *  Input:   h in [0,360] or undefined (-1), s and v in [0,1]        *
 * ----------------------------------------------------------------- */
void Color::setHSV(float h, float s, float v)
{
	if (Eq0(s, _CEPSILON)) {		// achromatic color there is no hue
		if (Eq(h, -1.0, _CEPSILON))	// undefined
			// this is the achromatic case
			_red = _green = _blue = v;
//		else
//			cerr << WARNINGSTR
//				<< "Color::setHSV invalid values" << endl;
	} else {				// chromatic color: there is a hue
		if (Eq(h, 360.0, _CEPSILON)) h = 0.0;
		h /= 60.0;			// h is now in [0,6]

		int i = (int)h;			// largest integer <= h
		float f = h - (float)i;		// fractional ray of h
		float p = v * (1-s);
		float q = v * (1 - s*f);
		float t = v * (1 - s*(1-f));

		switch (i) {
			case  0:
				_red   = v;
				_green = t;
				_blue  = p;
				break;

			case  1:
				_red   = q;
				_green = v;
				_blue  = p;
				break;

			case  2:
				_red   = p;
				_green = v;
				_blue  = t;
				break;

			case  3:
				_red   = p;
				_green = q;
				_blue  = v;
				break;

			case  4:
				_red   = t;
				_green = p;
				_blue  = v;
				break;

			case  5:
				_red   = v;
				_green = p;
				_blue  = q;
				break;

			default: ;
//				cerr << WARNINGSTR <<
//				"Color::setHSV invalid values" << endl;
		}  /* end of switch */
	}
} // setHSV

/* ========================= HSV ================================== *

			      ^ Value
		       _______|______
		       \      |      /   Side face of hexcone HSV
			 \    |    /
			   \  |  /
			     \|/
		    Green          Yellow               /
			-------------                 /
		      / \           / \             /
		    /     \       /     \         /  ^  Hue
		  /         \   /         \           \
	  Cyan  /_____________x_____________\   Red  --------->
		\       White  \            /                Saturation
		  \         /    \        /
		    \     /        \    /
		      \ /            \/          Top face of hexcone HSV
		  Blue  ~~~~~~~~~~~~~~ Magenta    color model

 * ================================================================ */


/* ================================================================= *
 *  getHSV(float& h, float& s, float& v);                                *
 *                                                                   *
 *  Input:   r,g,b, each in [0,1]                                    *
 *  Output:  h in [0,360],  s and v in [0,1], except if s = 0        *
 *           then h is undefined which is defined as a constant -1   *
 *           outside the interval [0,360]                            *
 * ================================================================= */
void Color::getHSV(float& h, float& s, float& v)
{
	float maximum = Max(_red,_green,_blue);
	float minimum = Min(_red,_green,_blue);

	v = maximum;	/* value */

	if (!Eq0(maximum, _CEPSILON))
		s = (maximum - minimum) / maximum;	// saturation
	else
		s = 0;

	if (Eq0(s, _CEPSILON))
		h = -1;		// hue, undefined
	else {
		// rc measures "distance" of color from red
		float rc = (maximum - _red) / (maximum - minimum);
		float gc = (maximum - _green) / (maximum - minimum);
		float bc = (maximum - _blue)  / (maximum - minimum);

		if (Eq(_red, maximum, _CEPSILON))
			// resulting color between yellow and magenta
			h = bc - gc;
		else
		if (Eq(_green, maximum, _CEPSILON))
			// resulting color between cyan and yellow
			h = 2 + rc - bc;
		else
			// resulting color between magenta and cyan
			h = 4 + gc - rc;

		h *= 60;		// convert to degrees

		if (h < 0) h += 360;	// make non negative
	}  /* chromatic case */
} // getHSV

/* ============================== Color3D ================================= */
void Color3D::operator () (dword c)
{
	Color32 col;
	double h, s, v;

	col.val = color = c & 0xFFFFFF;
	rgb2hsv(col, &h, &s, &v);
	dark   = hsv2rgb(h, s, (1.0 + v)/2.0);
	bright = hsv2rgb(h, s, 1.0);
} /* operator () */
