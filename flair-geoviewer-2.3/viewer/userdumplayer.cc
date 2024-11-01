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
 * Author"	Vasilis.Vlachoudis@cern.ch
 */

#include "random.h"
#include "viewer.h"
#include "geometry.h"
#include "userdumplayer.h"

UserDumpLayer::UserDumpLayer(const Geometry& g, GeometryKernel& k, GeometryViewer& v)
	: Layer(g, k, v),
	  _mtime(0)
{
	reset();
} // UserDumpLayer

UserDumpLayer::~UserDumpLayer()
{
	cleanup();
} // ~UserDumpLayer

/** cleanup */
void UserDumpLayer::cleanup()
{
	_filename.clear();
	_mtime = 0;
} // cleanup

/** reset */
void UserDumpLayer::reset()
{
	Random random(3141592654L);
	for (int i=0; i<NPARTICLES; i++)
		_color[i] = random.integer();

	// fix some colors
	_color[ 1+NEG_PARTICLES] = 0xFF0000;	// proton
	_color[ 3+NEG_PARTICLES] = 0x0000FF;	// electron
	_color[ 4+NEG_PARTICLES] = 0x00FFFF;	// positron
	_color[ 7+NEG_PARTICLES] = 0xFFFF00;	// photon
	_color[ 8+NEG_PARTICLES] = 0x00FF00;	// neutron
} // reset

/** open filename for userdump plotting */
bool UserDumpLayer::open(const char *name)
{
	struct stat st;

	if (_filename.compare(name))
		cleanup();

	// Check if information is valid
	if (stat(name, &st)) {
		cleanup();
		return true;
	}
	if ( st.st_mtime==_mtime) return false;

	cleanup();
	_filename = name;
	if (userdump.open(name)) return false;

	// Remember modification time
	_mtime = st.st_mtime;

	return true;
} // filename

/** draw */
void UserDumpLayer::draw(Painter& painter)
{
	if (!userdump) return;
	userdump.seek(0);

	int i = 0;

	while (i<start && userdump.readEvent(0)) {
		if (stop()) return;
		if (userdump.type == MGDRAW_SOURCE) i++;
	}

	i = n+1;
	while (i>0 && userdump.readEvent(-1)) {
		if (stop()) return;
		if (userdump.type == MGDRAW_TRACK) {
			dword col = color(userdump.event.track.particle);
			if ((col&COLOR_ALPHA) == COLOR_ALPHA) continue;

			UserDumpTrackPos& prev = userdump.track[0];
			double up = view().xyz2u(prev.x, prev.y, prev.z);
			double vp = view().xyz2v(prev.x, prev.y, prev.z);
			for (int j=1; j<userdump.track.size(); j++) {
				UserDumpTrackPos& pos = userdump.track[j];
				double u = view().xyz2u(pos.x, pos.y, pos.z);
				double v = view().xyz2v(pos.x, pos.y, pos.z);
				double cup = up;	// clip coordinates
				double cvp = vp;
				double cu  = u;
				double cv  = v;
				view().clipLine(&cup, &cvp, &cu, &cv);
				painter.line(view().u2i(cup),view().v2j(cvp),
				             view().u2i(cu), view().v2j(cv),
				             col);
				up = u;
				vp = v;
			}
		} else
		if (userdump.type == MGDRAW_SOURCE) i--;
	}
} // draw
