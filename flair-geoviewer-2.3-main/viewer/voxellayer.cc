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

#include "voxel.h"
#include "viewer.h"
#include "geometry.h"
#include "voxellayer.h"

/** draw */
void VoxelLayer::draw(Painter& painter)
{
	Point	r(0.0, 0.0, 0.0);
	dword *ptr = painter.data();
	int H = painter.height();	// save variables that might change
	int W = painter.width();	// during the scanning...

	for (int j=0; j<H; j++) {
		if (stop()) break;

		r.y = view().j2v(j);

		for (int i=0; i<W; i++, ptr++) {
			if (*ptr!=geometry.voxelColor && *ptr!=hashColor)
				continue;

			r.x = view().i2u(i);
			Point p = view().matrix() * r;

			bool ok;
			dword color = kernel.voxel.color(p, &ok);
			if (!ok) continue;
			// FIXME transparent color in voxels
			//if (color==COLOR_TRANSPARENT && show3D)
			//	color = geometry.transparentColor;

			if (*ptr == hashColor)
				*ptr = Darken(color, level) | FLAG_VOXEL;
			else
				*ptr = color | FLAG_VOXEL;
		}
	}
#ifdef EXPERIMENTAL
	drawStructure(painter);
#endif
} // draw

#ifdef EXPERIMENTAL
/** drawStructure */
void VoxelLayer::drawStructure(Painter& )
{
	for (int i=0; i<geometry.voxel.nroi(); i++) {
		if (kernel.voxel.roiShow(i)) {
		}
	}
} // drawStructure
#endif
