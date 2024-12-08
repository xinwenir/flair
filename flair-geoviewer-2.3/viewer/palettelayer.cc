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
 */

#include <math.h>

#include "viewer.h"
#include "painter.h"
#include "palette.h"
#include "geometry.h"
#include "viewport.h"
#include "palettelayer.h"

#define MARGIN	4
#define FORMAT "%g"

static const double LOG10[10] = {
			0.0,
			0.0,
			0.301029995663981,
			0.477121254719662,
			0.602059991327962,
			0.698970004336019,
			0.778151250383644,
			0.845098040014257,
			0.903089986991944,
			0.954242509439325 };

/** PaletteLayer */
PaletteLayer::PaletteLayer(const Geometry& g, GeometryKernel& k, GeometryViewer& v) : Layer(g, k, v)
{
	for (int i=0; i<MAXPALETTE; i++)
		_display[i] = false;
	font.load(DEFAULT_FONT);
} // PaletteLayer

/** draw */
void PaletteLayer::draw(Painter& painter)
{
	int w = painter.width();
	for (int i=MAXPALETTE-1; i>=0; i--)
		if (display(i)) {
			int dw = drawPalette(painter, i, w);
			w -= dw + MARGIN;
		}
} // draw

/** draw palette on painter
 * @param painter to draw to
 * @param id index of palette
 * @param width palette width in pixels
 * @return width of palette drawn
 */
int PaletteLayer::drawPalette(Painter& painter, int id, int width)
{
	assert(id>=0 && id<MAXPALETTE);
	Palette& palette = _palette[id];

	palette.checkLimits();

	// Find vertical limits
	int fy  = font.height();
	int fy2 = fy/2;

	int h     = painter.height() - 4*fy - 2*MARGIN;
	int ylow  = painter.height() - 3*fy + MARGIN;
	int yhigh = ylow - h;

	double y     = (double)ylow;
	double ystep = (double)h / (double)palette.size();

	// Find maximum text width
	double cbrange = palette.range();
	int textwidth = 0;
	if (palette.log()) {
		// Log
		double min = (double)((int)palette.min() - 1);
		double max = (double)((int)palette.max() + 1);
		double vv = pow(10.0,min);
		for (double v=min; v<=max; v++, vv*=10.0) {
			if (InRange(palette.min(), v, palette.max()))
				textwidth = Max(textwidth,
					painter.measuref(font, FORMAT, vv));
		}
	} else {
		// Linear
		double r = 10.0;
		while (true) {
			int steps = Round(cbrange/r);
			if (steps < 2)
				r /= 10.0;
			else
			if (steps >= 20)
				r *= 10.0;
			else break;
		}
		double min = (double)(int)(palette.min() / r) * r;
		for (double v=min; v<=palette.max(); v+=r) {
			if (v<palette.min() || v>palette.max()) continue;
			textwidth = Max(textwidth,
				painter.measuref(font, FORMAT, v));
		}
	}

	if (!_label[id].empty())
		textwidth += 1 + font.height();

	// Find horizontal limits
	int cbwidth = 15;
	int x1 = width - cbwidth - 2*MARGIN - textwidth;
	int x2 = x1 + cbwidth;
	int x3 = width - MARGIN;

	// Draw background if any
	if (viewer.textBackgroundLevel)
		painter.levelShiftRect(
				x1-MARGIN, yhigh-fy2-MARGIN,
				x3,        ylow +fy2+MARGIN,
				viewer.textBackgroundLevel);

	// Draw color box
	if (palette.invert())
		if (palette.interpolate()) {
			double dp = palette.max()-palette.min();
			for (int yy=yhigh; yy<=ylow; yy++) {
				double value = (ylow-yy)*dp/(double)h + palette.min();
				painter.line(x1, yy, x2, yy, palette.color(value));
			}
		} else
			for (int i=palette.size()-1; i>=0; i--, y-=ystep)
				painter.fillRect(x1, Round(y-ystep), x2, Round(y), palette[i]);
	else {
		if (palette.interpolate()) {
			double dp = palette.max()-palette.min();
			for (int yy=ylow; yy>=yhigh; yy--) {
				double value = (ylow-yy)*dp/(double)h + palette.min();
				painter.line(x1, yy, x2, yy, palette.color(value));
			}
		} else
			for (int i=0; i<palette.size(); i++, y-=ystep)
				painter.fillRect(x1, Round(y-ystep), x2, Round(y), palette[i]);
	}

	painter.rectangle(x1, Round(y), x2, ylow, 0);

	// Draw tics
	cbrange = palette.range();
	if (palette.log()) {
		// Log
		double min = (double)((int)palette.min() - 1);
		double max = (double)((int)palette.max() + 1);
		double vv = pow(10.0,min);
		for (double v=min; v<=max; v++, vv*=10.0) {
			if (InRange(palette.min(), v, palette.max())) {
				int ypos = ylow - (int)((double)h*(v-palette.min())/cbrange);
				painter.line(x2-cbwidth/2, ypos, x2, ypos, 0);
				painter.printf(font, x2+MARGIN, ypos-fy2-1,
						geometry.gridTextColor, FORMAT, vv);

				for (int i=2; i<10; i++) {
					double vl = v + LOG10[i];
					if (vl<palette.min() || vl>palette.max()) continue;
					ypos = ylow - (int)((double)h*(vl-palette.min())/cbrange);
					painter.line(x2-cbwidth/4, ypos, x2, ypos, 0);
				}
			}
		}
	} else {
		// Linear
		double r = 10.0;
		while (true) {
			int steps = Round(cbrange/r);
			if (steps < 2)
				r /= 10.0;
			else
			if (steps >= 20)
				r *= 10.0;
			else break;
		}
		double min = (double)(int)(palette.min() / r) * r;
		for (double v=min; v<=palette.max(); v+=r) {
			if (v<palette.min() || v>palette.max()) continue;
			int ypos = ylow - (int)((double)h*(v-palette.min())/cbrange);
			painter.line(x2-cbwidth/2, ypos, x2, ypos, 0);
			painter.printf(font, x2+MARGIN, ypos-fy2-1,
					geometry.gridTextColor, FORMAT, v);
		}
	}

	// Draw label
	if (!_label[id].empty()) {
		int xpos = x3 - font.height();
		int ypos = (ylow + yhigh + font.measure(_label[id].c_str()))/2;
		painter.vprintf(font, xpos, ypos, geometry.gridTextColor, _label[id].c_str());
	}

	return x3 - x1 + MARGIN;
} // drawPalette
