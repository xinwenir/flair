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

#ifndef __PAINTER_H
#define __PAINTER_H

#include <string>
#include <stdint.h>

#include "os.h"

enum FillType {
	FILL_FLOOD,
	FILL_DOTS,
	FILL_HASH,
	FILL_HASHR,
	FILL_X,
	FILL_X2
};

/* rectangular clipping region */
struct ClipRegion {
	int  left;	/* X-upper left corner	*/
	int  top;	/* Y-upper left corner	*/
	int  right;	/* X-lower right corner	*/
	int  bottom;	/* Y-lower right corner	*/
};

/* out codes for clipping lines */
union _OutCodeUnion {		/* outcodes are represented as bit fields */
	struct {
		unsigned code0 : 1;	/* x < left	*/
		unsigned code1 : 1;	/* y < top	*/
		unsigned code2 : 1;	/* x > right	*/
		unsigned code3 : 1;	/* y > bottom	*/
	}  ocs;
	int outcodes;
};

/** integer point class */
class IPoint {
public:
	word	x, y;
public:
	IPoint() : x(0), y(0) {}
	IPoint(int i,int j) : x(i), y(j) {}
};

class Painter;
class ViewPort;

/* BFont class */
class BFont {
private:
	std::string	_name;
	int	imageWidth;		/** image width			*/
	int	imageHeight;		/** image height		*/
	int	_width;			/** fixed width of each char	*/
	int	_height;		/** height of each char		*/
	uint8_t*	imageData;		/** Font image data		*/
	uint8_t	charLength[256];	/** var-width of each character	*/

public:
	BFont() : _width(-1), _height(-1), imageData(NULL) {}
	BFont(const char *filename) :
			_width(-1),
			_height(-1),
			imageData(NULL)
			{ load(filename); }
	~BFont();

const	std::string& name()	const	{ return _name; }
	int	width()		const	{ return _width; }
	int	height()	const	{ return _height; }

	bool	load(const char *filename);
	void	set(const char* n, int w, int h, uint8_t* data);
	void	set(const char* n, int w, int h, dword* data);

	int	draw(Painter& painter, int x, int y, const dword color, const uint8_t ch);
	int	vdraw(Painter& painter, int x, int y, const dword color, const uint8_t ch);
	int	drawOutline(Painter& painter, int x, int y,
			const dword color, const dword outline, const uint8_t ch);

	int	measure(const char *str);

protected:
	uint8_t*	charPtr(const int ch) const
			{ return imageData + ((ch>>4)*_height*imageWidth) + (ch&0x0F)*_width; }
	void	build();
	void	clean();
}; // BFont

/** Painter base class
 * basic painter class, provide means of clipping lines
 * Holds the data of an Image ready to be used as XImage
 * Default color depth is 24, therefore bits per pixel 32
 */
class Painter {
protected:
	int	_width;		/** width in pixels		*/
	int	_height;		/** height in pixels		*/
	int	_dataSize;	/** width*height*sizeof(*data)	*/
	int	_maxSize;	/** allocated memory size bytes	*/
	dword*	_data;		/** pointer to allocated memory	*/
	dword	_background;	/** background (gray) color	*/
	ClipRegion _clip;	/** clip region			*/

public:
	Painter(const int w, const int h);
	~Painter();

	void	init(const int w, const int h);

	// Data
	void	dataNull()		{ _data = NULL; }	// WARNING needed to deallocate from X11

	dword	*data()			{ return _data; }
	int	width()		const	{ return _width; }
	int	height()	const	{ return _height; }

	/** @param byte	background byte (only gray scale are allowed) */
	void	background(dword b)	{ _background = b; }
	dword	background()	const	{ return _background; }

	// Clipping
	void	resetClip();
	void	clip(int l, int t, int r, int b);
const	ClipRegion& clip()	const	{ return _clip; }

	bool	insideClip(const int x, const int y) const {
			return _clip.left<=x && x<=_clip.right &&
			       _clip.top <=y && y<=_clip.bottom;
		}
	void	clipPoint(int *x, int *y) const {
			*x = Range(_clip.left, *x, _clip.right);
			*y = Range(_clip.top,  *y, _clip.bottom);
		}
	bool	clipLine(int *x1, int *y1, int *x2, int *y2) const;

	// Drawing
	void	clear();
	bool	hasPixel()	const { return true; }
	void	pixel(const int x, const int y, dword color)
			{ _data[y*_width+x] = color; }
	void	pixel(const int x, const int y, dword color, double intensity);
	dword	pixel(const int x, const int y) const
			{ return _data[y*_width+x]; }
	dword	operator()(const int x, const int y)	const	{ return pixel(x,y); }
	dword	*pixelPtr(const int x, const int y)
			{ return _data + y*_width+x; }

	void	move(const int dx, const int dy);

	/** draw a line
	 * @param x1,y1	starting point of line
	 * @param x2,y2	end point of line
	 * @param color	color of line
	 * @return true if line is inside clip area
	 */
	bool	line(int x1, int y1, int x2, int y2, const dword color) {
			if (clipLine(&x1, &y1, &x2, &y2)) {
				unclippedLine(x1,y1,x2,y2,color);
				return true;
			}
			return false;
		}
	bool	lineThick(int x1, int y1, int x2, int y2, const int w, const dword color) {
			if (clipLine(&x1, &y1, &x2, &y2)) {
				unclippedThickLine(x1,y1, x2,y2, w, color);
				return true;
			}
			return false;
		}
	bool	lineAntialias(int x1, int y1, int x2, int y2, const dword color) {
			if (clipLine(&x1, &y1, &x2, &y2)) {
				unclippedLineAntialias(x1,y1, x2,y2, color);
				return true;
			}
			return false;
		}
	bool	drawBitmap(int, int, int, int, const dword *);

	bool	rectangle(int, int, int, int, const dword);
	void	fillRect(int, int, int, int, const dword);
	void	levelShiftRect(int, int, int, int, int);
	void	fill(int, int, const dword, const dword = 0, const FillType = FILL_FLOOD);


	int	measure(BFont &font, const char *str) { return font.measure(str); }
	int	measuref(BFont &font, const char *, ...);
	int	printf(BFont&, int, int, const dword, const char *, ...);
	int	vprintf(BFont&, int, int, const dword, const char *, ...);
	int	printfShadow(BFont&, int, int, int, const dword, const dword, const char *, ...);
	int	printfOutline(BFont&, int, int, const dword, const dword, const char *, ...);
	void	drawchar(BFont& font, int x, int y, const dword color, const char ch)
			{ font.draw(*this, x, y, color, ch); }
	void	vdrawchar(BFont& font, int x, int y, const dword color, const char ch)
			{ font.vdraw(*this, x, y, color, ch); }
	size_t	memory() const { return sizeof(Painter) + _maxSize; }

protected:
	void	unclippedLine(int x1, int y1, int x2, int y2, const dword color);
	void	unclippedThickLine(int x1, int y1, int x2, int y2, const int w, const dword color);
	void	unclippedLineAntialias(double x1, double y1, double x2, double y2, const dword color);

friend class ViewPort;
}; // Painter

#endif
