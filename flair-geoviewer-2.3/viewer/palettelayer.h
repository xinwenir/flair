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

#ifndef __PALETTE_LAYER_H
#define __PALETTE_LAYER_H

#include "os.h"
#include "layer.h"
#include "palette.h"

class Painter;
class ViewPort;

/**
 * PaletteLayer
 *
 * This class draws the palette
 */
class PaletteLayer : public Layer {
private:
	/* color scale */
	Palette	_palette[MAXPALETTE];		/** color palette		*/
	bool	_display[MAXPALETTE];		/** display palette		*/
	std::string	_label[MAXPALETTE];	/** palette label		*/
	int	_default;			/** default palette		*/
//	bool	cbVertical;			/** orientation			*/
public:
	BFont	font;				/** color palette font		*/

public:
	PaletteLayer(const Geometry& g, GeometryKernel& k, GeometryViewer& v);

	/* set/get */
	bool	display(int i)		const	{ return _display[i]; }
	void	display(int i, bool a)		{ _display[i] = a; }

	int	def()			const	{ return _default; }
	void	def(int a)			{ _default = a; }

	void	label(int idx, std::string l)	{ _label[idx] = l; }
	std::string label(int idx)	const	{ return _label[idx]; }

	void	draw(Painter& painter);

	Palette& operator [](const int i)	{ return _palette[i]; }

private:
	int	drawPalette(Painter& painter, int id, int width);
}; // PaletteLayer

#endif
