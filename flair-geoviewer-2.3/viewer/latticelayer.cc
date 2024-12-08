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
#include "geometry.h"
#include "latticelayer.h"

/** draw */
void LatticeLayer::draw(Painter& painter)
{
	int	H = painter.height();	// save variables that might change
	int	W = painter.width();	// during the scanning...
	Ray	ray;
	VZone*	repZone=NULL;	// replica zone

	double dx = view().uv2dx(1.0,0.0);
	double dy = view().uv2dy(1.0,0.0);
	double dz = view().uv2dz(1.0,0.0);

	GeometryEngine *eng = kernel.engine();

	dword *ptr = painter.data();
	for (int j=0; j<H; j++) {
		if (stop()) break;
		double v = view().jc2v(j);
		for (int i=0; i<W; i++, ptr++) {
			if (*ptr!=geometry.latticeColor && *ptr!=viewer.lattice.hashColor)
				continue;
			double u = view().ic2u(i);

			/* scan to find end of lattice row */
			dword *ps = ptr;
			int is = i;
			for (; i<W; i++, ptr++)
				if (*ptr!=geometry.latticeColor && *ptr!=viewer.lattice.hashColor)
					break;
			if (i==W) ptr--;	// avoid double counting at the end of the line

			double tmax = view().ic2u(i) - u;
			//double u2 = u + tmax/2.0;

			// Find zone and transformation
			double x = view().uv2x(u,v);
			double y = view().uv2y(u,v);
			double z = view().uv2z(u,v);

			geometry.lockRead();;
			eng->incBodyCheckId();
			repZone = eng->whereRay(x,y,z, dx,dy,dz, tmax/2.0, repZone);
			geometry.unlockRead();;

			if (repZone == NULL) continue;

			// prepare the ray
			ray.init();
			RaySegment segment;
			segment.pos.set(x, y, z);
			segment.dir.set(dx, dy, dz);
			segment.zone = repZone;
			segment.tmin = SMALL3D;
			segment.tmax = tmax;
			ray.push(segment);

			dword color;
			bool end;
			int k;
			do {
				if (ray.hitZone() == NULL) {
					// try to recover zone
					geometry.lockRead();;
					ray.segment().zone = eng->whereRay(
								ray.segment().pos,
								ray.segment().dir,
								ray.segment().tmin);
					geometry.unlockRead();;
					color = geometry.errorColor;
				} else {
					if (ray.hitRegion()->transparent() && viewer.d3.show)
						color = geometry.transparentColor();
					else
						color = ray.hitRegion()->color();
				}
				dword color2 = Darken(color, level);

				if (ray.hitZone() != NULL) {
					end = eng->intersectRay(&ray, true);
					k   = view().u2ic(u+ray.T())-is;
				} else {
					end = true;
					k   = i-is;
				}
				while (k>0 && ps<ptr) {
					*ps = (*ps==viewer.lattice.hashColor)?color2:color;
					ps++;
					is++;
					k--;
				}
			} while (!end);
			// fill remaining pixels
			while (ps<ptr) *ps++ = color;
		}
	}
} // draw
