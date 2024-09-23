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

#ifndef __USRBIN_LAYER_H
#define __USRBIN_LAYER_H

#include "usrbin.h"
#include "layer.h"

/* =============================== UsrbinLayer ============================== */
class UsrbinLayer : public Layer {
protected:
	Usrbin	usrbins[MAXUSRBIN];	/** usrbin objects		*/
	int	_alpha[MAXUSRBIN];	/** usrbin transparency		*/
	int	_palette[MAXUSRBIN];	/** palette index		*/

public:
	UsrbinLayer(const Geometry& g, GeometryKernel& k, GeometryViewer& v);

	/* set/get */
	int	alpha(int i)		const	{ return _alpha[i]; }
	void	alpha(int i, int a)		{ _alpha[i] = a; }

	int	palette(int i)		const	{ return _palette[i]; }
	void	palette(int i, int p)		{ _palette[i] = p; }

	Usrbin&	operator [](const int i)	{ return usrbins[i]; }

	void	draw(Painter& painter);
	uint8_t	shade(const Point& hit, Color& color) const;
	void	regionColorFromUsrbin(int idx);
}; // UsrbinLayer

#endif
