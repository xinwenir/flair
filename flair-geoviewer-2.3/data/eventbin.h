/*
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
 * Author:      Wioletta.Kozlowska@cern.ch
 * based on the T.Boehlen version from 2011
 * Last change: 16-Dec-2016
 */
#ifndef EVENTBIN_H
#define EVENTBIN_H

#include <string>
#include <vector>

#include "usrbin.h"
#include "fortran.h"

#define FTNVGET(V,T,B)	{V=*(T*)B; B+=sizeof(T);}
//TTB: return V=Variable, T=Type, B=buffer,  return S=String, L=Length

class Eventbin {
private:
	int	_detector;		// detector loaded
	UsrbinType _type;		// type of eventbin
	int	_score;			// quantity scored
public:
	int	nx, ny, nz;
	double	xlow, ylow, zlow;
	double	xhigh, yhigh, zhigh;
	double	dx, dy, dz;
	double	xofs,  yofs,  zofs;	// offset of usrbin
	double	x0,  y0;		// for R-Z-Phi binning
	float	weight;
	int	nhits;

private:
	FortranFile	file;
	std::string	filename;

public:
	Eventbin() {};
	Eventbin(std::string fn, int det=0) { load(fn,det); };
	~Eventbin() {cleanup();};

	void cleanup();
	int load(std::string fn, int det=0);

	/* get/set */
	int	detector()	const	{ return _detector; }
	UsrbinType type()	const	{ return _type; }
	int	score()		const	{ return _score; }

	int eventbini(const double x) const {return int((x-xlow)/dx);}
	int eventbinj(const double y) const {return int((y-ylow)/dy);}
	int eventbink(const double z) const {return int((z-zlow)/dz);}

	int readHeader();
	char* readEvent(char *buffer, unsigned sizeBuffer);
}; // class Eventbin

#endif
