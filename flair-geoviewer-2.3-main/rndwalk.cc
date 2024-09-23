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
 * Date:	04-Feb-2010
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>


#include <fstream>
#include <iostream>

#include "os.h"
#include "nox.h"
#include "timer.h"
#include "random.h"
#include "particle.h"

using namespace std;

/** main */
int main(int ac, char *av[])
{
	if (ac!=3) {
		cerr << "syntax: " << av[0] << " <file.nox> <events>" << endl;
		return 1;
	}
	fpetrap();
	NoxGeometry world(av[1]);

	GeometryEngine* engine = world.viewer.kernel().engine();

	Random random(1234567890);
	int nevents = atoi(av[2]);

	Timer timer;
	timer.start();

#if 0
	for (int i=0; i<nevents; i++) {
		double x, y, z;
		//random.normal(&x,&y);
		//random.sincos(&x,&y);
		//random.disc(2.0,&x,&y);
		random.sphere(2.0,&x,&y,&z);
		cout << x << " " << y << " " << z << endl;
	}
#endif

	for (int i=0; i<nevents; i++) {
		Caloron particle;
		Vector pos = random.sphere(2.0);
		cout << pos.x << " " << pos.y << " " << pos.z << endl;
		VZone* zone = NULL;

		while (true) {
			Vector dir = random.vector();		// random direction
			double step = random.normal();		// random step

			engine->incBodyCheckId();
			zone = engine->whereRay(pos, dir, SMALL3D, zone);

			particle.init();
			RaySegment segment(pos, dir, zone);
			segment.tmax = step;
			particle.push(segment);

			// find next intersection
			if (engine->intersectRay(&particle, true)) {
				pos = particle.hit();
				cout << pos.x << " " << pos.y << " " << pos.z << endl;
			} else {
				break;
			}
		}
		cout << endl << endl;
	}

	timer.stop();
	cout << "# Timer1:" << timer << endl;
} // main
