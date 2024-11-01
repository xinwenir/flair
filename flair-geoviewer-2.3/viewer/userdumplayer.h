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

#ifndef __MGDRAW_LAYER_H
#define __MGDRAW_LAYER_H

#include <string>
#include <sys/stat.h>

#include "layer.h"
#include "userdump.h"

#define NEG_PARTICLES	6			// negative particles
#define NPARTICLES	64+NEG_PARTICLES

/* ============================= UserDumpLayer ============================== */
class UserDumpLayer : public Layer {
protected:
	UserDump	userdump;	// userdump class
	time_t		_mtime;
	std::string	_filename;	// filename

public:
	int	start;			// starting source event
	int	n;			// how many events to show

	dword	_color[NPARTICLES];	// particle color transparent = no show
	double	_emin[NPARTICLES];	// energy limits to show
	double	_emax[NPARTICLES];

public:
	UserDumpLayer(const Geometry& g, GeometryKernel& k, GeometryViewer& v);
	~UserDumpLayer();

	void	cleanup();

	/* set/get */
const	std::string&	filename()	const	{ return _filename; }
	bool	open(const char *name=NULL);
	void	reset();

	void	display(int particle)		{ _color[particle+NEG_PARTICLES] &= ~COLOR_ALPHA; }
	void	hide(int particle)		{ _color[particle+NEG_PARTICLES] |=  COLOR_ALPHA; }
	void	color(int particle, dword c)	{ _color[particle+NEG_PARTICLES] = c; }
	dword	color(int particle)	const	{
			if (particle>=-NEG_PARTICLES && particle<NPARTICLES)
				return _color[particle+NEG_PARTICLES];
			else
				return COLOR_ALPHA;
		}

	void	draw(Painter& painter);
}; // UserDumpLayer

#endif
