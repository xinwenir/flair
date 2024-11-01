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

#ifndef __IMAGE_LAYER_H
#define __IMAGE_LAYER_H

#include "layer.h"
#include "matrix3.h"

// Drawing type
#define	SHOW_OFF	0
#define	SHOW_LATE	1
#define SHOW_PROMPT	2

/* ================================ ImageLayer ============================== */
class ImageLayer : public Layer {
protected:
	dword	*_data;			/** image data			*/
	dword	 _size;			/** allocated size		*/
	int	 _width;		/** dimensions			*/
	int	 _height;
	int	 _alpha;		/** image transparency		*/
	double	 _minW;			/** minimum and maximum w	*/
	double	 _maxW;			/** minimum and maximum w	*/
	dword	 _background;		/** background color (guessed)	*/
	Matrix3	 R;			/** matrix (plane)->image	*/
	Matrix3	 M;			/** matrix (x,y,z)->plane	*/
private:
	bool	 _visible;		/** if image is visible		*/
//	Matrix4	 _viewMatrix;		/** viewing matrix of image	*/

public:
	ImageLayer(const Geometry& g, GeometryKernel& k, GeometryViewer& v);
	~ImageLayer();

	/* set/get */
	bool	data(dword *dat, size_t size);
	void	size(const int w, const int h) {
			_width  = w;
			_height = h;
		}
	size_t	size()		const		{ return _width*_height; }
	void	alpha(const int a)		{ _alpha = a; }
	void	matrix(const Matrix3 *r=NULL, const Matrix3 *m=NULL);
	void	colorRange(dword low, dword high);
	void	findBackground();

	void	draw(Painter& painter);
}; // ImageLayer

#endif
