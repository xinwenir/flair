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

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "os.h"
#include "nox.h"

#ifdef MEM
#include "memory.h"
#endif

#define ROTSTEP 180.0/12.0
//#define ROTSTEP	15.0
//#define ROTSTEP	1.0

//#define _DRAW_2D
#define _DRAW_3D
//#define _DRAW_BBOX

using namespace std;

static Display *display;
static int	screen;
static Visual  *visual;
static int	depth;
static Window	window;
static GC	gc;
//static Pixmap	pixmap;
//static int	pixmap_width=0, pixmap_height=0;
static unsigned black, white;
static XImage	*ximage;

NoxGeometry nox;

/** initX */
void initX()
{
	/* use the information from the environment variable DISPLAY
	   to create the X connection:
	*/
	display = XOpenDisplay(NULL);
	screen  = DefaultScreen(display);
	visual  = DefaultVisual(display, screen);
	depth   = DefaultDepth(display, screen);
	black   = BlackPixel(display, screen),	/* get color black */
	white   = WhitePixel(display, screen);  /* get color white */
	assert(depth>16);

	/* once the display is initialized, create the window.
	   This window will be have be 640 pixels across and 480 down.
	   It will have the foreground white and background black
	*/
	window = XCreateSimpleWindow(display, DefaultRootWindow(display),
			0, 0, nox.viewer.width(), nox.viewer.height(),
			5, black, white);

	/* here is where some properties of the window can be set.
	   The third and fourth items indicate the name which appears
	   at the top of the window and the name of the minimized window
	   respectively.
	*/
	XSetStandardProperties(display, window, "XGeoview","xgeoview",None,NULL,0,NULL);

	/* this routine determines which types of input are allowed in
	   the input.  see the appropriate section for details...
	*/
	XSelectInput(display, window, ExposureMask|ButtonPressMask|KeyPressMask|StructureNotifyMask);

	/* clear the window and bring it on top of the other windows */
	XClearWindow(display, window);
	XMapRaised(display, window);

	XGCValues values;
	values.function = GXcopy;
	values.plane_mask = AllPlanes;
	values.foreground = black;
	values.background = white;
	gc = XCreateGC(display, window,
		GCFunction | GCPlaneMask | GCForeground | GCBackground,
		&values);

	ximage = XCreateImage(display, visual, depth, ZPixmap, 0,
			(char*)nox.viewer.painter.data(),
			nox.viewer.width(), nox.viewer.height(), 32, 0);

} // initX

/** main */
int main(int ac, char *av[])
{
	Matrix4 mat;
	mat.identity();
	nox.kernel.initThreads(0);

	if (ac<2) {
		printf("%s <file.nox> [<width:640> <height:480>]\n", av[0]);
		return 1;
	}
	fpetrap();
	nox.load(av[1]);
	nox.dump = true;
	int width=640, height=480;
	if (ac>2) width  = atoi(av[2]);
	if (ac>3) height = atoi(av[3]);
	nox.viewer.resize(width, height);
	nox.viewer.view().calcWindow(4.0);
	//nox.viewer.xray = 100;
	//nox.viewer.supersampling = 4;

	nox.viewer.font.load("fonts/fixed8x13.tga");
	nox.viewer.decoration.gridFont.load("fonts/fixed8x13.tga");
	nox.viewer.palette.font.load("fonts/fixed8x13.tga");

	cout << endl << "Processing:" << endl;
	//nox.geometry.updateOBBs();
	nox.viewer.d2.project();

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

	initX();
	int configCount = 0;


	while (1) {
		XEvent event;		/* the XEvent declaration !!! */
		KeySym key;		/* a dealie-bob to handle KeyPress Events */
		char text[255];		/* a char buffer for KeyPress Events */

		/* get the next event and stuff it into our event variable.
		   Note:  only events we set the mask for are detected!
		*/
		XNextEvent(display, &event);

		if (event.type==Expose && event.xexpose.count==0) {
			/* the window was exposed redraw it! */
//			if (configCount>2)
//				configCount--;
//			else
			if (width!=nox.viewer.width() || height!=nox.viewer.height()) {
				nox.viewer.resize(width, height);
				ximage->width  = nox.viewer.width();
				ximage->height = nox.viewer.height();
				ximage->data   = (char*)nox.viewer.painter.data();
				ximage->bytes_per_line = 0;
				XInitImage(ximage);
				nox.viewer.draw(mask, true);
//				XClearWindow(display, window);
			}
			XPutImage(display, window, gc, ximage, 0, 0, 0, 0,
				nox.viewer.width(), nox.viewer.height());
//			XFlush(display);
		} else
		if (event.type==KeyPress) {
			int len = XLookupString(&event.xkey,text,255,&key,0);
			/* use the XLookupString routine to convert the invent
			   KeyPress data into regular text.  Weird but necessary...
			*/
			double u = nox.viewer.view().Uofs();	// old x
			double v = nox.viewer.view().Vofs();	// old y
			switch (text[0]) {
				case '=':
					nox.viewer.view().zoom(nox.viewer.view().zoom()*1.2);
					cout << "Zoom=" << nox.viewer.view().zoom() << endl;
					if (nox.viewer.view().zoom() > 1000.0) {
						nox.viewer.view().calcWindow(4.0);
						nox.viewer.d2.project();
					}
					break;

				case '-':
					nox.viewer.view().zoom(nox.viewer.view().zoom()/1.2);
					cout << "Zoom=" << nox.viewer.view().zoom() << endl;
					break;

				case 'h':
					u -= (nox.viewer.width()/10) / nox.viewer.view().Sx();
					//u -= nox.viewer.view().extX / nox.viewer.view().zoom()/5.0;
					nox.viewer.view().offset(u,v);
					break;

				case 'j':
					v -= (nox.viewer.height()/10) / nox.viewer.view().Sy();
					nox.viewer.view().offset(u,v);
					break;

				case 'k':
					v += (nox.viewer.height()/10) / nox.viewer.view().Sy();
					nox.viewer.view().offset(u,v);
					break;

				case 'l':
					u += (nox.viewer.width()/10) / nox.viewer.view().Sx();
					nox.viewer.view().offset(u,v);
					break;

				case 'H':
					u -= (nox.viewer.width()/100) / nox.viewer.view().Sx();
					nox.viewer.view().offset(u,v);
					break;

				case 'J':
					v -= (nox.viewer.height()/100) / nox.viewer.view().Sy();
					nox.viewer.view().offset(u,v);
					break;

				case 'K':
					v += (nox.viewer.height()/100) / nox.viewer.view().Sy();
					nox.viewer.view().offset(u,v);
					break;

				case 'L':
					u += (nox.viewer.width()/100) / nox.viewer.view().Sx();
					nox.viewer.view().offset(u,v);
					break;

				case 'n':
					mat.translate(0.0, 0.0, (nox.viewer.width()/50) / nox.viewer.view().Sx());
					nox.viewer.view().transform(mat);
					nox.viewer.d2.project();
					break;

				case 'm':
					mat.translate(0.0, 0.0, -(nox.viewer.width()/50) / nox.viewer.view().Sx());
					nox.viewer.view().transform(mat);
					nox.viewer.d2.project();
					break;

				case 's': // rotate left
					mat.rotY(RAD(-ROTSTEP), 1e-5);
					nox.viewer.view().transform(mat);
					nox.viewer.d2.project();
					break;

				case 'd': // rotate right
					mat.rotY(RAD(ROTSTEP), 1e-5);
					nox.viewer.view().transform(mat);
					nox.viewer.d2.project();
					break;

				case 'e': // rotate up
					mat.rotX(RAD(-ROTSTEP), 1e-5);
					nox.viewer.view().transform(mat);
					nox.viewer.d2.project();
					break;

				case 'x': // rotate down
					mat.rotX(RAD(ROTSTEP), 1e-5);
					nox.viewer.view().transform(mat);
					nox.viewer.d2.project();
					break;

				case 'w': // rotate aclk
					mat.rotZ(RAD(-ROTSTEP), 1e-5);
					nox.viewer.view().transform(mat);
					nox.viewer.d2.project();
					break;

				case 'r': // rotate clk
					mat.rotZ(RAD(ROTSTEP), 1e-5);
					nox.viewer.view().transform(mat);
					nox.viewer.d2.project();
					break;

				// check image shift
//				case 's':
//					nox.viewer.painter->move(-50, 0);
//					nox.viewer.painter->setClip(-51, 0, -1, -1);
//					x += (int)(50.0 / nox.viewer.view().Sx());
//					nox.viewer.view().setOrigin(x,y);
//					break;
//
//				case 'd':
//					nox.viewer.painter->move(50, 0);
//					nox.viewer.painter->setClip(50, 0, -1, -1);
//					break;

//				case 'e':
//					nox.viewer.painter->move(0, -50);
//					nox.viewer.painter->setClip(0, -50, -1, -1);
//					break;

//				case 'x':
//					nox.viewer.painter->move(0, 50);
//					nox.viewer.painter->setClip(0, 50, -1, -1);
//					break;

				case 'c':
					nox.viewer.view().offset(0.,0.);
					break;

				case 'z': // zoom test
					cout << "zoom test\n" ;
					cout << "view=" << nox.viewer.view() << endl;
					nox.viewer.view().calcWindow(2.0);
					cout << "\nafter\n";
					cout << "view=" << nox.viewer.view() << endl;
					nox.viewer.d2.project();
					break;

				case 'q':
					nox.viewer.stopThread();
					goto END;

				case 'b': nox.viewer.stopThread();
					break;

				default:
					printf("You pressed: ");
					for (int i=0; i<len; i++)
						printf("%02X ",text[i]);
					printf("\n");
			}
			if (nox.viewer.view().invalidWindow()) {
				nox.viewer.view().calcWindow(4.0);
				cout << "Invalid window" << endl;
				nox.viewer.d2.project();
				cout << "Transform" << endl;
			}
			nox.viewer.draw(mask, true);
			cout << nox.viewer.view() << endl;
			cout << "BasisX: " << nox.viewer.view().matrix(0,0)
					   << " " << nox.viewer.view().matrix(1,0)
					   << " " << nox.viewer.view().matrix(2,0) << endl;
			cout << "BasisY: " << nox.viewer.view().matrix(0,1)
					   << " " << nox.viewer.view().matrix(1,1)
					   << " " << nox.viewer.view().matrix(2,1) << endl;
			nox.viewer.painter.resetClip();
			XPutImage(display, window, gc, ximage, 0, 0, 0, 0,
				nox.viewer.width(), nox.viewer.height());
		} else
		if (event.type==ButtonPress) {
			/* tell where the mouse Button was Pressed */
			double x, y, z;
			nox.viewer.view().ij2xyz(event.xbutton.x, event.xbutton.y, &x, &y, &z);
			printf("Location [%d,%d] = [%g, %g, %g]\n",
				event.xbutton.x, event.xbutton.y, x, y, z);
		} else
		if (event.type==ConfigureNotify) {
			width  = event.xconfigure.width;
			height = event.xconfigure.height;
			configCount++;
		}
	}
END:
	nox.viewer.painter.dataNull();
	XDestroyImage(ximage);
	//XFreePixmap(display, pixmap);
	XFreeGC(display, gc);
	XDestroyWindow(display,window);
	XCloseDisplay(display);
	return 0;
} // main
