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

#include "viewer.h"
#include "usrbin.h"
#include "geometry.h"
#include "usrbinlayer.h"

UsrbinLayer::UsrbinLayer(const Geometry& g, GeometryKernel& k, GeometryViewer& v) : Layer(g, k, v)
{
	for (int i=0; i<MAXUSRBIN; i++) {
		_alpha[i]   = 0;
		_palette[i] = 0;
	}
} // UsrbinLayer

/** draw
 * FIXME 1. Can be optimized by introducing a bounding box around the image and calculate
 *          first the intersection of the 6 planes
 *       2. Modified to check for non-lines, in order to draw also the regions
 */
void UsrbinLayer::draw(Painter& painter)
{
	Point	r(0.0, 0.0, 0.0);
	int H = painter.height();	// save variables that might change
	int W = painter.width();	// during the scanning...
	dword *ptr = painter.data();	// data pointer

//	Matrix4 rot[MAXUSRBIN];
	int maxusrbin = -1;
	for (int usrbinidx=0; usrbinidx < MAXUSRBIN; usrbinidx++) {
		Usrbin &usrbin = usrbins[usrbinidx];
		if (!usrbin.hasData() && !usrbin.checker()) continue;
//	FIXME might be caching the rotation*viewerport could be helpful
//		if (usrbin.rotdefi)
//			rot[usrbinidx].multiply(*geometry.rotdefi(usrbin.rotdefi),
//						view().matrix);
		maxusrbin = usrbinidx;
	}
	if (maxusrbin<0) return;

	for (int j=0; j<H; j++) {
		if (stop()) break;
		r.y = view().jc2v(j);

		for (int i=0; i<W; i++, ptr++) {
			dword flag = *ptr & FLAG_INFOMASK;
			if (flag==0) continue;
			if (viewer.d3.show && flag!=FLAG_REGION) continue;
			r.x = view().ic2u(i);

			for (int usrbinidx=0; usrbinidx <= maxusrbin; usrbinidx++) {
				Usrbin &usrbin = usrbins[usrbinidx];
				if (!usrbin.hasData() && !usrbin.checker()) continue;
				Point p = view().matrix() * r;

				bool ok;
				double value = usrbin.getData(p.x, p.y, p.z, &ok);
				if (!ok) continue;

				Color32 pix, usr;
				pix.val = *ptr;
				pix.rgb.alpha = alpha(usrbinidx);

				int pal = _palette[usrbinidx];
				if (usrbin.checker()) {
					if (value<0.5)
						usr.val = viewer.palette[pal].first();
					else
						usr.val = viewer.palette[pal].last();
				} else {
					dword color = viewer.palette[pal].color(value);
					if (color == COLOR_TRANSPARENT)
						break;
					else
						usr.val = color;
				}
				usr.rgb.alpha = 0xFF;
				*ptr = alphaBlend(pix, usr) | FLAG_REGION;
			}
		}
	}
} // draw

/** shade
 * @param hit	input 3D hit position
 * @param color output color at hit position
 * @return 0xFF if no usrbin is found at hit position otherwise 0..255 alpha level
 */
uint8_t UsrbinLayer::shade(const Point& hit, Color& color) const
{
	// move a bit the hit position to avoid boundary effects
	for (int i=0; i<MAXUSRBIN; i++) {
		const Usrbin& bin = usrbins[i];
		if (!bin.hasData() && !bin.checker()) continue;

		bool ok;
		double value = bin.getData(hit.x, hit.y, hit.z, &ok);
		if (ok) {
			Color usrcolor;
//			color->rgb.alpha = alpha(i);
			int pal = _palette[i];
			dword c;
			if (bin.checker()) {
				if (value<0.5)
					c = viewer.palette[pal].first();
				else
					c = viewer.palette[pal].last();
			} else {
				c = viewer.palette[pal].color(value);
				if (c == COLOR_TRANSPARENT) break;
			}
			usrcolor.set(c);
//			usr.rgb.alpha = 0xFF;
			color.blend(usrcolor, (float)alpha(i)/255.0);
//			color->val = alphaBlend(*color, usr);
			return alpha(i);
		}
	}
	return 0xFF;
} // shade

/** regionColorFromUsrbin
 * Set color of all regions using the USRBIN values
 * if not defined set to 0
 */
void UsrbinLayer::regionColorFromUsrbin(int idx)
{
	if (!usrbins[idx].hasData()) return;
	int pal = _palette[idx];
	for (int i = 0; i < kernel.regions.size(); i++) {
		VRegion* region = kernel.getRegion(i);
		region->color(viewer.palette[pal].first());
		bool ok;
		double value = usrbins[idx].getData(region->id(),&ok);
		if (ok) region->color(viewer.palette[pal].color(value));
	}
} // regionColorFromUsrbin
