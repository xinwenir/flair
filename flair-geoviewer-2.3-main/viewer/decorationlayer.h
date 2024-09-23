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

#ifndef __DECORATION_LAYER_H
#define __DECORATION_LAYER_H

#include "os.h"
#include "layer.h"

class Painter;
class ViewPort;
class Geometry;
class GeometryViewer;

/**
 * DecorationLayer
 *
 * This class draws all the decorations of a viewport and provides an interface
 * to show messages and the progress bar.
 *
 *  The draw method draws the following:
 *  - Grid
 *  - Axes
 *  - Title
 */
class DecorationLayer : public Layer {
public:
	/* flags */
	bool	showAxes;		/** draw axes			*/
	bool	showTitle;		/** flag ON/OFF for title	*/
	bool	showGrid;		/** flag ON/OFF for grid	*/

	/* grid */
	double	grid_x,  grid_y,  grid_z;/** real positions of grid	*/
	double	grid_dx, grid_dy, grid_dz;/** real positions of grid	*/
	double	grid_u,  grid_v;	/** starting location of grid	*/
	double	grid_du, grid_dv;	/** step of grid		*/
	char	gridU;			/** axis to plot		*/
	char	gridV;			/** axis to plot		*/
	int	gridLevel;		/** grid level			*/
	BFont	gridFont;		/** drawing font for grid/axes	*/

public:
	DecorationLayer(const Geometry& g, GeometryKernel& k, GeometryViewer& v) : Layer(g, k, v) {
			showTitle = true;
			showAxes  = true;
			showGrid  = true;

			grid_du = grid_dv = 100.0;
			gridU = 'U';
			gridV = 'V';
			gridLevel = 200;
			gridFont.load(DEFAULT_FONT);
		}

	void draw(Painter& painter, const dword mask);
	void drawMessage(Painter& painter, const char *msg, dword color, int line);
	void drawProgress(Painter& painter, int percent, const char *msg);

private:
	void drawAxes(Painter& painter);
	void drawGrid(Painter& painter);
	void drawTitle(Painter& painter);
}; // DecorationLayer

#endif
