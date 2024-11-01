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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "os.h"
#include "nox.h"

#ifdef MEM
#include "memory.h"
#endif

#include "mesh.h"

//#define _DRAW_2D
#define _DRAW_3D
//#define _DRAW_BBOX

using namespace std;

/** main */
int main(int ac, char *av[])
{
#ifdef MEM
	Memory::init();
#endif
	if (ac<3) {
		printf("%s <file.nox> <file-out.rgb> [<width:640> <height:480>]\n", av[0]);
		return 1;
	}
	fpetrap();

	NoxGeometry nox;
//	nox.kernel.initThreads(0);

	nox.load(av[1]);
	nox.dump = true;
	int width=640, height=480;
	if (ac>3) width  = atoi(av[3]);
	if (ac>4) height = atoi(av[4]);
	nox.viewer.resize(width, height);
	nox.viewer.view().calcWindow(4.0);
	//nox.viewer.xray = 100;
	//nox.viewer.supersampling = 4;

	nox.viewer.font.load("fonts/fixed8x13.tga");
	nox.viewer.decoration.gridFont.load("fonts/fixed8x13.tga");
	nox.viewer.palette.font.load("fonts/fixed8x13.tga");

	cout << endl << "Processing:" << endl;
	//nox.geometry.updateOBBs();
	//nox.viewer.d2.project();

#if _DEBUG>0
	nox.geometry.printMemory();
	nox.viewer.printMemory();
#endif

	dword mask = DRAW_CLEAR;

#ifdef _DRAW_2D
	nox.viewer.d2.fillRegions = true;
	nox.viewer.lattice.show   = true;
	mask |= DRAW_SEGMENTS;
	mask |= DRAW_REGIONS;
	mask |= DRAW_LATTICES;
	//mask |= DRAW_VOXEL;
#endif

#ifdef _DRAW_3D
	nox.viewer.d3.show        = true;
	nox.viewer.lattice.show   = true;
	nox.viewer.voxel.show     = true;
	//nox.viewer.d3.drawEdges   = true;
	mask |= DRAW_3D;
#endif

#ifdef _DRAW_BBOX
	mask |= DRAW_BBOX;
#endif
	mask |= DRAW_GRID;
	mask |= DRAW_AXES;

	cout << endl << "Processing:" << endl;
	nox.viewer.d2.project();
	nox.viewer.draw(mask, true);

#if _DEBUG>0
	nox.geometry.printMemory();
	nox.viewer.printMemory();
#endif

	printf("Writing file %s\n", av[2]);
	FILE *fout = fopen(av[2],"wb");
	dword *pixel = nox.viewer.painter.data();
	for (int i=0; i<width*height; i++, pixel++) {
		Color32 p;
		p.val = *pixel;
		fputc(p.rgb.red,   fout);
		fputc(p.rgb.green, fout);
		fputc(p.rgb.blue,  fout);
	}
	fclose(fout);

#ifdef MEM
	Memory::fini();
#endif
} // main
