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
 * Date:	14-Jan-2014
 */

#ifndef __COLOR_H
#define __COLOR_H

#include <string>
#include <ostream>
#include <cstdint>

#include "os.h"
#define _CEPSILON		1E-7

#define COLOR_BLACK		0x00000000
#define COLOR_WHITE		0x00FFFFFF
#define COLOR_TRANSPARENT	0xFFFFFFFF
#define COLOR_ALPHA		0xFF000000

/* int breakup into color uint8_ts for little endian*/
typedef union {
	struct {
//#if CPU_ENDIAN == LITTLE_ENDIAN
		uint8_t blue;
		uint8_t green;
		uint8_t red;
		uint8_t alpha;
//#else
//		byte alpha;
//		byte red;
//		byte green;
//		byte blue;
//#endif
	} rgb;
	dword val;
} Color32;

inline int Red(const int col) {
//#if CPU_ENDIAN == LITTLE_ENDIAN
	return (col>>16)&0xff;
//#else
//	return col&0xff;
//#endif
}
inline int Green(const int col) {
	return (col>>8)&0xff;
}
inline int Blue(const int col) {
//#if CPU_ENDIAN == LITTLE_ENDIAN
	return col&0xff;
//#else
//	return (col>>16)&0xff;
//#endif
}

inline int Alpha(const dword col) {
	Color32 color;
	color.val = col;
	return color.rgb.alpha;
}

/** darken by color in the range 0..level */
inline dword Darken(const dword col, const int level) {
		Color32 color;
		color.val = col;
		color.rgb.red   = (uint8_t)(((int)color.rgb.red   * level)>>8);
		color.rgb.green = (uint8_t)(((int)color.rgb.green * level)>>8);
		color.rgb.blue  = (uint8_t)(((int)color.rgb.blue  * level)>>8);
		return color.val;
	}

/** lighten by color in the range level..255*/
inline dword Lighten(const dword col, const int level) {
		Color32 color;
		color.val = col;
		color.rgb.red   = (uint8_t)((((int)color.rgb.red   * (256-level)) >> 8) + level);
		color.rgb.green = (uint8_t)((((int)color.rgb.green * (256-level)) >> 8) + level);
		color.rgb.blue  = (uint8_t)((((int)color.rgb.blue  * (256-level)) >> 8) + level);
		return color.val;
	}

/** intensity */
inline int Intensity(const dword col){
	return ((77*Red(col) + 151*Green(col) + 28*Blue(col))>>8);
	//return (int)(0.3*(double)Red(col) + 0.59*(double)Green(col) + 0.11*(double)Blue(col));
}

/** make RGB color */
inline dword RGB(const int red, const int green, const int blue)
		{ return red<<16 | green<<8 | blue; }

/** make RGBA color */
inline dword RGBA(const int red, const int green, const int blue, const int alpha)
		{ return alpha<<24 | red<<16 | green<<8 | blue; }

/* adjust color level */
dword ColorLevel(const dword col, const Color32& lw, const Color32& hg);

/* alphaBlend - transparent addition of colors */
dword alphaBlend(const Color32& colA, const Color32& colB);

/* rgb2hsv */
void rgb2hsv(const Color32& col, double *h, double *s, double *v);

/* hsv2rgb */
dword hsv2rgb(double h, double s, double v);

/* ============================== Color32 ============================= */
//class Color32 {
//public:
//public:
//}; // Color32

/* =============================== Color ============================== */
// Colors are using SINGLE precision numbers "float"
// No need for the extra digits
// WARNING: Colors can have values from 0.0 to +inf (more than 1.0)
class Color {
private:
	// The components.
	float	_red;	/** red	ranging from 0..1	*/
	float	_green;	/** gree			*/
	float	_blue;	/** blue			*/

public:	// Some basic colors
static const Color	Black;
static const Color	Red;
static const Color	Green;
static const Color	Blue;
static const Color	Magenta;
static const Color	Turquoise;
static const Color	Yellow;
static const Color	White;

public:
	Color(const float R=0.0, const float G=0.0, const float B=0.0) { set(R,G,B); }
	Color(const double R, const double G, const double B) { set(R,G,B); }
	Color(const int R, const int G, const int B) { set(R,G,B); }
	Color(const dword color) { set(color); }
	Color(const Color32& color) { set(color); }

	// The components in cartesian coordinate system.
	float red()	const { return _red; }
	float green()	const { return _green; }
	float blue()	const { return _blue; }

	void red(const float R)		{ _red   = Max(0.0f, R); }
	void green(const float G)	{ _green = Max(0.0f, G); }
	void blue(const float B)	{ _blue  = Max(0.0f, B); }
	void red(const int R)		{ _red   = Max(0.0f, (float)R/255.0f); }
	void green(const int G)		{ _green = Max(0.0f, (float)G/255.0f); }
	void blue(const int B)		{ _blue  = Max(0.0f, (float)B/255.0f); }

	// Get components by index
	float operator () (int) const;

	// return color as dword
	dword	color32()	const	{
			Color32 color;
			color.rgb.red   = (uint8_t)Min(255, (int)(255.0*red()));
			color.rgb.green = (uint8_t)Min(255, (int)(255.0*green()));
			color.rgb.blue  = (uint8_t)Min(255, (int)(255.0*blue()));
			return color.val;
		}

	// Set the components in cartesian coordinate system.
	void set(const Color& color) {
			_red   = color._red;
			_green = color._green;
			_blue  = color._blue;
		}

	void set(const float R, const float G, const float B) {
			_red   = Max(0.0f, R);
			_green = Max(0.0f, G);
			_blue  = Max(0.0f, B);
		}

	void set(const double R, const double G, const double B) {
			_red   = Max(0.0, R);
			_green = Max(0.0, G);
			_blue  = Max(0.0, B);
		}

	void set(const int R, const int G, const int B) {
			_red   = Max(0.0, (float)R/255.0);
			_green = Max(0.0, (float)G/255.0);
			_blue  = Max(0.0, (float)B/255.0);
		}

	void set(const Color32& color) {
			_red   = (float)color.rgb.red  /255.0;
			_green = (float)color.rgb.green/255.0;
			_blue  = (float)color.rgb.blue /255.0;
		}

	void set(const dword color) {
			Color32 col;
			col.val = color;
			set(col);
		}

	Color& operator = (const Color& p) {
			_red   = p._red;
			_green = p._green;
			_blue  = p._blue;
			return *this;
		}

	// Comparisons
	bool operator == (const Color& c) const {
			return (Eq(_red,   c._red,   _CEPSILON) &&
				Eq(_green, c._green, _CEPSILON) &&
				Eq(_blue,  c._blue,  _CEPSILON)) ? true : false;
		}

	bool operator != (const Color& c) const {
			return (Eq(_red,   c._red,   _CEPSILON) &&
				Eq(_green, c._green, _CEPSILON) &&
				Eq(_blue,  c._blue,  _CEPSILON)) ? false : true;
		}

	// Addition.
	// No range limitation since colors can have any positive value!
	Color& operator += (const Color& p) {
			_red   += p._red;
			_green += p._green;
			_blue  += p._blue;
			return *this;
		}

	// Subtraction.
	Color& operator -= (const Color& p) {
			_red   = Max(0.0f, _red   - p._red);
			_green = Max(0.0f, _green - p._green);
			_blue  = Max(0.0f, _blue  - p._blue);
			return *this;
		}

	// Scaling with real numbers.
	Color& operator *= (float a) {
			_red   *= a;
			_green *= a;
			_blue  *= a;
			return *this;
		}

	// Scaling with real numbers.
	Color& operator *= (double a) {
			_red   *= a;
			_green *= a;
			_blue  *= a;
			return *this;
		}

	// Scaling with real numbers.
	Color& operator /= (float a) {
			float inva = 1.0/a;
			_red   *= inva;
			_green *= inva;
			_blue  *= inva;
			return *this;
		}

	// Scaling with another color.
	Color& operator *= (const Color& color) {
			_red   *= color._red;
			_green *= color._green;
			_blue  *= color._blue;
			return *this;
		}

	void	blend(Color& color, float alpha) {
			_red   = red()  *alpha + color.red()  *(1.0-alpha);
			_green = green()*alpha + color.green()*(1.0-alpha);
			_blue  = blue() *alpha + color.blue() *(1.0-alpha);
		}

	// Color transformations
	void	setHSV(float, float, float);
	void	getHSV(float&, float&, float&);
	float	gray();
	float	grey()	const { return 0.3*_red  + 0.59*_green + 0.11*_blue; };
}; // Color

// Output to a stream.
inline std::ostream& operator << (std::ostream& os, const Color& col) {
		os << "[" << col.red() << ", " << col.green() << ", " << col.blue() << "]";
		return os;
	}

// Addition of 3-colors
inline Color operator + (const Color& a, const Color& b) {
		return Color(a.red()+b.red(), a.green()+b.green(), a.blue()+b.blue());
	}

// Subtraction of 3-colors.
inline Color operator - (const Color& a, const Color& b) {
		return Color(a.red()-b.red(), a.green()-b.green(), a.blue()-b.blue());
	}

// Scaling of a color with a real number
inline Color operator * (const Color& p, float a) {
		return Color(a*p.red(), a*p.green(), a*p.blue());
	}

inline Color operator * (float a, const Color& p) {
		return Color(a*p.red(), a*p.green(), a*p.blue());
	}

inline Color operator * (double a, const Color& p) {
		return Color(a*p.red(), a*p.green(), a*p.blue());
	}

// Scaling of a color with a real number
inline Color operator / (const Color& p, float a) {
		return Color(p.red()/a, p.green()/a, p.blue()/a);
	}

inline Color operator / (float a, const Color& p) {
		return Color(a/p.red(), a/p.green(), a/p.blue());
	}

// Scaling with another color
inline Color operator * (const Color& a, const Color& b) {
		return Color(a.red()*b.red(), a.green()*b.green(), a.blue()*b.blue());
	}

/** Color3D that holds 4 variations of the same color for
 * pseudo 3D drawing
 */
class Color3D {
public:
	dword	color;			/** default color		*/
	dword	dark;			/** darker color		*/
	dword	bright;			/** brighter color		*/
public:
	Color3D() : color(0),  dark(0), bright(0) {}
	Color3D(dword c)		{ (*this)(c); }
	void	operator ()(dword c);
	dword	operator ()()	const	{ return color; }
}; // Color3D

#endif
