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

#include <stdlib.h>

#include "viewer.h"
#include "geometry.h"
#include "imagelayer.h"

/** ImageLayer */
ImageLayer::ImageLayer(const Geometry& g, GeometryKernel& k, GeometryViewer& v) : Layer(g, k, v)
{
	// background image
	_data    = NULL;
	_size    = 0;
	_width   = 0;
	_height  = 0;
	_alpha   = 0x7F;
	_minW    = -INFINITE;
	_maxW    = INFINITE;
	_visible = false;
	R.identity();
	M.identity();
} // ImageLayer

/** ~ImageLayer */
ImageLayer::~ImageLayer()
{
	if (_size>0) delete [] _data;
} /* ~ImageLayer */

/** data */
bool ImageLayer::data(dword *dat, size_t siz)
{
	if (siz==0) {
		// set the pointer DO NOT COPY
		if (_size>0) delete [] _data;
		_data = dat;
		_size = 0;
		return true;
	} else
	if (siz != _size) {
		if (_data) delete [] _data;
		try {
			_data = new dword[siz];
			_size = siz;
		} catch ( ... ) {
			_data = NULL;
			_size = 0;
			return false;
		}
	}
	memcpy(_data, dat, siz);
	return true;
} // data

/** matrix
 * set image Matrix and prepare the inverse matrix and check for visibility
 */
void ImageLayer::matrix(const Matrix3 *r, const Matrix3 *m)
{
	if (r != NULL) R.copy(*r);
	if (m != NULL) M.copy(*m);

	// Calculate inverse and image matrices
	//_viewMatrix.multiply(_matrix, view().matrix());

	/* Check matrix and view matrix are on the same projection */
	Vector imgW(view().invMatrix(2,0), view().invMatrix(2,1), view().invMatrix(2,2));
	Vector w(M(2,0), M(2,1), M(2,2));
	_visible = Abs(imgW*w)/w.length() > 0.9999;
} // matrix

/** colorRange */
void ImageLayer::colorRange(dword low, dword high)
{
	Color32 lw, hg;
	lw.val = low;
	hg.val = high;

	dword *ptr = _data;
	for (size_t i=0; i<size(); i++, ptr++)
		*ptr = ColorLevel(*ptr, lw, hg);
} // colorRange

/** findBackground */
void ImageLayer::findBackground()
{
	int sum[3*256];
	dword color[3*256];
	memset(sum, 0, sizeof(sum));
	memset(color, 0, sizeof(color));
	dword* d = _data;
	for (size_t i=0; i<size(); i++) {
		int s = (*d & 0xff) + (*d>>8 & 0xff) + (*d>>16 & 0xff);
		sum[s]++;
		color[s] = *d;
		d++;
	}
	int mi = 0;
	for (int i=0; i<3*256; i++)
		if (sum[i]>sum[mi])
			mi = i;

	_background = color[mi];
} // findBackground

/** draw
 * FIXME 1. Can be optimized by introducing a bounding box around the image and calculate
 *          first the intersection of the 6 planes
 *       2. Modified to check for non-lines, in order to draw also the regions
 *       3. Rotation of the USRBIN
 */
void ImageLayer::draw(Painter& painter)
{
	if (!_visible) return;

	dword *ptr = painter.data();
	int H = painter.height();	// save variables that might change
	int W = painter.width();	// during the scanning...

	for (int j=0; j<H; j++) {
		if (stop()) break;
		double v = view().jc2v(j);
		for (int i=0; i<W; i++, ptr++) {
			if ((*ptr & FLAG_INFOMASK) == 0) continue;
			double u = view().ic2u(i);

			Point p = view().uv2xyz(u,v);
			double ui = M(0,0)*p.x + M(0,1)*p.y + M(0,2)*p.z;
			double vi = M(1,0)*p.x + M(1,1)*p.y + M(1,2)*p.z;

			double xi = R(0,0)*ui + R(0,1)*vi + R(0,2);
			double yi = R(1,0)*ui + R(1,1)*vi + R(1,2);

			int pi = (int)xi;
			int pj = (int)yi;

			if (pi<0 || pi>=_width || pj<0 || pj>=_height) continue;
			int loc = pj*_width + pi;

			Color32	pixel, pixelR, pixelD, pixelDR;
			pixel.val = _data[loc];
			if (pi < _width-1)
				pixelR.val = _data[loc+1];
			else
				pixelR.val = pixel.val;

			if (pj < _height-1) {
				pixelD.val = _data[loc+_width];
				if (pi < _width-1)
					pixelDR.val = _data[loc+_width+1];
				else
					pixelDR.val = pixelR.val;
			} else {
				pixelD.val = pixel.val;
				pixelDR.val = pixelR.val;
			}

			double a = xi - (double)pi;
			double b = yi - (double)pj;
			double a1 = 1.0 - a;
			double b1 = 1.0 - b;

			// Bilinear smoothing
			int red   = (int)(a1*b1*(double)pixel.rgb.red +
			                  a *b1*(double)pixelR.rgb.red +
			                  a1*b *(double)pixelD.rgb.red +
			                  a *b *(double)pixelDR.rgb.red);
			int green = (int)(a1*b1*(double)pixel.rgb.green +
			                  a *b1*(double)pixelR.rgb.green +
			                  a1*b *(double)pixelD.rgb.green +
			                  a *b *(double)pixelDR.rgb.green);
			int blue  = (int)(a1*b1*(double)pixel.rgb.blue +
			                  a *b1*(double)pixelR.rgb.blue +
			                  a1*b *(double)pixelD.rgb.blue +
			                  a *b *(double)pixelDR.rgb.blue);

			Color32 pix, img;
			pix.val = *ptr;
			pix.rgb.alpha = _alpha;
			img.val = RGBA(red,green,blue,0xFF);
			*ptr = alphaBlend(pix, img);
		}
	}
} // draw
