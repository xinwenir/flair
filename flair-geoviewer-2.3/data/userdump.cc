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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>

#include "userdump.h"

using namespace std;

///////////////////////////////// UserDump ////////////////////////////////////
/** open */
bool UserDump::open(const char *filename)
{
	return file.open(filename,"rb");
} // open

/** Read or skip next event from mgread structure */
bool UserDump::readEvent(int filter)
{
	char	buffer[30];
	int	length;

	// Read header
	length = file.read(buffer,sizeof(buffer));
	if (length!=20) return false;

	FortranParser parser(buffer, length);
	event.generic.ndum = parser.readInt();
	event.generic.mdum = parser.readInt();
	event.generic.jdum = parser.readInt();
	event.generic.edum = parser.readFloat();
	event.generic.wdum = parser.readFloat();

	/* check event */
	if (event.generic.ndum > 0) {		// Tracking event
		type = MGDRAW_TRACK;
		if (filter & MGDRAW_TRACK)
			return readTracking();
		else
			file.skip();
	} else
	if (event.generic.ndum == 0) {	// Energy deposition event
		type = MGDRAW_ENERGY;
		if (filter & MGDRAW_TRACK)
			return readEnergy();
		else
			file.skip();
	} else {
		type = MGDRAW_SOURCE;
		event.source.ncase = -event.generic.ndum;	// Make positive
		if (filter & MGDRAW_SOURCE)
			return readSource();
		else
			file.skip();
	}
	return true;
} // readEvent

/** readTracking */
bool UserDump::readTracking()
{
	/* read ntrack points */
	int size = sizeof(UserDumpTrackPos) * (event.track.ntrack+1) +
		   sizeof(float) * event.track.mtrack +
		   sizeof(float);

	if (!file.mustBe(size)) return false;
	track.clear();
	for (int i=0; i<=event.track.ntrack; i++) {
		UserDumpTrackPos pos;
		file.readraw(&pos, sizeof(pos), 1);
		track.append(pos);
	}
	for (int i=0; i<event.track.mtrack; i++) {
		float dh;
		file.readraw(&dh, sizeof(dh), 1);
	}
	file.readraw(&event.track.ctrack, sizeof(float), 1);
	if (!file.mustBe(size)) return false;

	return true;
} // readTracking

/** readEnergy */
bool UserDump::readEnergy()
{
	char buffer[16];

	if (file.read(buffer,sizeof(buffer)) != 16) return false;

	FortranParser parser(buffer);
	event.energy.x    = parser.readFloat();
	event.energy.y    = parser.readFloat();
	event.energy.z    = parser.readFloat();
	event.energy.edep = parser.readFloat();
	return true;
} // readEnergy

/** readSource */
bool UserDump::readSource()
{
	/* read npfluka particles */
	int size = sizeof(UserDumpSourceParticle) * event.source.npfluka;
	if (!file.mustBe(size)) return false;

	stack.clear();

	for (int i=0; i<event.source.npfluka; i++) {
		UserDumpSourceParticle part;
		size_t sz = fread(&part, sizeof(part), 1, file.handle);
		if (sz != 1) return false;
		stack.append(part);
	}
	return file.mustBe(size);
} // readSource
