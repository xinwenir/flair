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

//#include "bmath.h"
#include "viewer.h"
#include "matrix4.h"
#include "painter.h"
#include "geometry.h"
#include "viewport.h"
#include "decorationlayer.h"

#define MARGIN	4
#if 0
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
#endif

//				  U        V        W,       X        Y        Z
static const int axisColor[] = {COLOR_U, COLOR_V, COLOR_W, COLOR_X, COLOR_Y, COLOR_Z };

/** draw */
void DecorationLayer::draw(Painter& painter, const dword mask)
{
	if (mask&DRAW_GRID && showGrid)
		drawGrid(painter);

	if (mask&DRAW_AXES && showAxes)
		drawAxes(painter);

	if (mask&DRAW_TITLE && showTitle)
		drawTitle(painter);
} // draw

/** drawGrid
 * FIXME: The need to access very often to the fields of the viewport makes me
 * think that big parts of this method should be moved to the ViewPort class
 */
void DecorationLayer::drawGrid(Painter& painter)
{
	dword *ptr;

	grid_du = 100.0;
	double w = (double)painter.width() / view().Sx();

	/* calculation horizontal grid (vertical lines) */
	if (w<grid_du)
		while (w<grid_du) grid_du /= 10.0;
	else {
		while (w>grid_du) grid_du *= 10.0;
		if (grid_du>w) grid_du /= 10.0;
	}

	double stepi = grid_du*view().Sx();
	if (painter.width() / stepi <= 2)
		grid_du /= 10.0;
	stepi = grid_du*view().Sx();

	/* calculation vertical grid (horizontal lines) */
	w = (double)painter.height() / view().Sy();
	grid_dv = grid_du;
	if (w<grid_dv) {
		while (w<grid_dv) grid_dv /= 10.0;
	} else {
		while (w>grid_dv) grid_dv *= 10.0;
		if (grid_dv>w) grid_dv /= 10.0;
	}

	double stepj = grid_dv*view().Sy();
	if (painter.height() / stepj <= 2) {
		grid_dv /= 10.0;
		stepj  = grid_dv*view().Sy();
	}

	// Make grids comparable
	if (painter.width() / stepi > 4 * painter.height() / stepj) {
		grid_dv /= 10.0;
		stepj  = grid_dv*view().Sy();
	} else
	if (4 * painter.width() / stepi < painter.height() / stepj) {
		grid_du /= 10.0;
		stepi  = grid_du*view().Sx();
	}

	// Find starting position
	double u = view().i2u(painter.clip().left);
	double v = view().j2v(painter.clip().top);
	double x = view().uv2x(u,v);
	double y = view().uv2y(u,v);
	double z = view().uv2z(u,v);
	double tmp;
	double starti, startj;

	// find starting u position
	switch (gridU) {
		case 'X':
			grid_dx = grid_du;
			grid_x = (double)((long long)(x/grid_du))*grid_du;
			view().xyz2uv(grid_x, y, z, &grid_u, &tmp);
			break;
		case 'Y':
			grid_dy = grid_du;
			grid_y = (double)((long long)(y/grid_du))*grid_du;
			view().xyz2uv(x, grid_y, z, &grid_u, &tmp);
			break;
		case 'Z':
			grid_dz = grid_du;
			grid_z = (double)((long long)(z/grid_du))*grid_du;
			view().xyz2uv(x, y, grid_z, &grid_u, &tmp);
			break;
		default:
			grid_u = (double)((long long)((u-view().Uofs())/grid_du))*grid_du + view().Uofs();
	}

	// convert to pixel
	starti = view().u2id(grid_u);
	int brk=0;
	while (starti<painter.clip().left) {
		starti += stepi;
		if (++brk > 100) {
			starti = painter.clip().left;
			break;
		}
	}

	// find starting v position
	switch (gridV) {
		case 'X':
			grid_dx = grid_dv;
			grid_x = (double)((long long)(x/grid_dv))*grid_dv;
			view().xyz2uv(grid_x, y, z, &tmp, &grid_v);
			break;
		case 'Y':
			grid_dy = grid_dv;
			grid_y = (double)((long long)(y/grid_dv))*grid_dv;
			view().xyz2uv(x, grid_y, z, &tmp, &grid_v);
			break;
		case 'Z':
			grid_dz = grid_dv;
			grid_z = (double)((long long)(z/grid_dv))*grid_dv;
			view().xyz2uv(x, y, grid_z, &tmp, &grid_v);
			break;
		default:
			grid_v = (double)((long long)((v-view().Vofs())/grid_dv))*grid_dv + view().Vofs();
	}
	// convert to pixel
	startj = view().v2jd(grid_v);
	brk = 0;
	while (startj<painter.clip().top) {
		startj += stepj;
		if (++brk > 100) {
			startj = painter.clip().top;
			break;
		}
	}

	// Draw grid lines
	int label=-1;
	int labely = painter.height() - gridFont.height() - 1;

	// For horizontal axis
	for (double i=starti; i<=painter.clip().right; i+=stepi) {
		int iu = Round(i);
		ptr = painter.pixelPtr(iu, painter.clip().top);
		for (int j=painter.clip().top; j<=painter.clip().bottom; j+=geometry.gridStep) {
			// ignore the j%(ystep)==0
			*ptr = Darken(*ptr, gridLevel);
			ptr += geometry.gridStep*painter.width();
		}
		if (label+3 < i) {
			double num;
			u = view().i2u(iu);
			switch (gridU) {
				case 'X':
					num = view().uv2x(u,v);
					break;
				case 'Y':
					num = view().uv2y(u,v);
					break;
				case 'Z':
					num = view().uv2z(u,v);
					break;
				default:
					num = u - view().Uofs();
			}
			num = (double)LongRound(num / grid_du)*grid_du;
			if (viewer.textBackgroundLevel) {
				painter.levelShiftRect(
					iu,
					labely-1,
					iu + 2 + painter.measuref(gridFont, "%.12g", num),
					labely + gridFont.height()+1,
					viewer.textBackgroundLevel);
			}
			label = painter.printf(gridFont, iu+1, labely,
					geometry.gridTextColor, "%.12g", num);
		}
	}
	painter.drawchar(gridFont,
			painter.clip().right-gridFont.width()-1, labely+1,
			0, gridU);
	painter.drawchar(gridFont,
			painter.clip().right-gridFont.width()-2, labely,
			axisColor[gridU-'U'], gridU);

	// For vertical axis
	u = view().i2u(painter.clip().left);
	for (double j=startj; j<=painter.clip().bottom; j+=stepj) {
		int jv = Round(j);
		ptr = painter.pixelPtr(painter.clip().left, jv);
		for (int i=painter.clip().left; i<=painter.clip().right; i+=geometry.gridStep) {
			*ptr = Darken(*ptr, gridLevel);
			ptr += geometry.gridStep;
		}
		double num;
		v = view().j2v(jv);
		switch (gridV) {
			case 'X':
				num = view().uv2x(u,v);
				break;
			case 'Y':
				num = view().uv2y(u,v);
				break;
			case 'Z':
				num = view().uv2z(u,v);
				break;
			default:
				num = v - view().Vofs();
		}
		num = (double)LongRound(num / grid_dv)*grid_dv;
		if (viewer.textBackgroundLevel)
			painter.levelShiftRect(
					0,
					jv-gridFont.height()+2,
					painter.measuref(gridFont, "%.12g", num),
					jv,
					viewer.textBackgroundLevel);
		painter.printf(gridFont, 1, jv-gridFont.height()+1,
				geometry.gridTextColor, "%.12g", num);
	}
	painter.drawchar(gridFont, 2, 2, 0, gridV);
	painter.drawchar(gridFont, 1, 1, axisColor[gridV-'U'], gridV);
} // drawGrid

/** drawAxes */
void DecorationLayer::drawAxes(Painter& painter)
{
	int axisX0 = geometry.axisLen+2;
	int axisY0 = painter.height() - axisX0;

	/* make an transparent circle */
	int ymin = Max(axisY0 - geometry.axisLen, 0) - axisY0;
	int ymax = Min(axisY0 + geometry.axisLen, painter.height()) - axisY0;

	for (int y=ymin; y<ymax; y+= 1) {
//		assert(InRange(0,y+axisY0,painter.height()-1));
		int dx = (int)isqrt((geometry.axisLen-y)*(geometry.axisLen+y));
		int xmin = Max(axisX0-dx, 0) - axisX0;
		int xmax = Min(axisX0+dx, painter.width()) - axisX0;
		dword *ptr = painter.pixelPtr(axisX0+xmin, axisY0+y);
		for (int x=xmin; x<xmax; x++, ptr++)
			*ptr = Darken(*ptr, gridLevel);
	}

	/* X - axis */
	double r = (double)geometry.axisLen;
	int xe = axisX0 + Round(r * view().matrix(0,0));
	int ye = axisY0 - Round(r * view().matrix(0,1));
	int color = view().matrix(0,2)>=0.0? axisColor[3] : axisColor[3] & 0x7F7F7F;
	painter.line(axisX0, axisY0, xe, ye, color);
	painter.drawchar(gridFont, xe+1, ye, color, 'x');

	/* Y - axis */
	xe = axisX0 + Round(r * view().matrix(1,0));
	ye = axisY0 - Round(r * view().matrix(1,1));
	color = view().matrix(1,2)>=0.0? axisColor[4] : axisColor[4] & 0x7F7F7F;
	painter.line(axisX0, axisY0, xe, ye, color);
	painter.drawchar(gridFont, xe+1, ye, color, 'y');

	/* Z - axis */
	xe = axisX0 + Round(r * view().matrix(2,0));
	ye = axisY0 - Round(r * view().matrix(2,1));
	color = view().matrix(2,2)>=0.0? axisColor[5] : axisColor[5] & 0x7F7F7F;
	painter.line(axisX0, axisY0, xe, ye, color);
	painter.drawchar(gridFont, xe+1, ye, color, 'z');
} // drawAxes

/** drawTitle */
void DecorationLayer::drawTitle(Painter& painter)
{
	int x = painter.width() -MARGIN-1 -painter.measure(viewer.font, viewer.title());
	painter.printf(viewer.font, x+1, MARGIN, 0, viewer.title());	// shadow
	painter.printf(viewer.font, x, MARGIN-1, geometry.titleColor, viewer.title());
} // drawTitle

/** drawMessage */
void DecorationLayer::drawMessage(Painter& painter, const char *msg, dword color, int line)
{
	int x = painter.width()/2 - painter.measure(viewer.font, msg)/2;
	int y = MARGIN + line * viewer.font.height();
	painter.printfShadow(viewer.font, x, y, 1, color, 0, msg);
} // drawMessage

/** drawProgress */
void DecorationLayer::drawProgress(Painter& painter, int percent, const char *msg)
{
	// Show progress bar3D
	int green = (percent*200)/100;
	int red   = 200 - green;
	dword color = RGB(red, green, 0);
	// Draw percentage bar + status
	painter.fillRect(0, MARGIN-1,
			painter.width()-1, MARGIN+viewer.font.height(),
			geometry.gridTextColor);
	painter.fillRect(0, MARGIN-1,
			((painter.width()-1)*percent)/100, MARGIN+viewer.font.height(),
			color);
	painter.printf(viewer.font, 36, MARGIN,geometry.titleColor,
			"%d%%  %s", percent, msg);
} // drawProgress
