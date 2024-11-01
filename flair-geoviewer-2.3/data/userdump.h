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
 * Date:	19-Oct-2011
 */

#ifndef __MGDRAW_H
#define __MGDRAW_H

#include <stdio.h>

#include "os.h"
#include "array.h"
#include "fortran.h"
#include "matrix4.h"

struct UserDumpTrackPos {
	float	x, y, z;
};
struct UserDumpSourceParticle {
	int	particle;
	float	energy;
	float	weight;
	float	x, y, z;
	float	tx, ty, tz;
};

#define MGDRAW_TRACK	0x01
#define MGDRAW_ENERGY	0x02
#define MGDRAW_SOURCE	0x04

class UserDump {
private:
	FortranFile	file;

public:
	int	type;				// type of event
	union {
		struct {
			int	ndum;		// tracking data>0, energy=0, source particle<0
			int	mdum;		// code
			int	jdum;		// jtrack (particle type)
			float	edum;		// energy
			float	wdum;		// weight
		} generic;			// generic read record

		struct {
			int	ntrack;
			int	mtrack;
			int	particle;
			float	energy;
			float	weight;
			float	ctrack;		// curved path
			//float	cx, cy, cz;
		} track;			// track record

		struct {
			int	_dummy;
			int	icode;
			int	particle;
			float	energy;
			float	weight;

			float	x, y, z;
			float	edep;
		} energy;			// energy deposition record

		struct {
			int	ncase;
			int	npfluka;
			int	nstmax;
			float	energy_sum;
			float	weight_sum;
		} source;			// source particles record
	} event;

	Array<UserDumpTrackPos>		track;
	Array<UserDumpSourceParticle>	stack;

public:
	UserDump() {}
//	~UserDump() {}

	bool	open(const char *filename);
	operator int()		{ return file; }
	int	seek(int pos)	{ return file.seek(pos, SEEK_SET); }

	bool	readEvent(int filter);

private:
	bool	readTracking();
	bool	readEnergy();
	bool	readSource();
//	void	cleanup();
}; // UserDump
#endif
