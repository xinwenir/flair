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
#include <unistd.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/extensions/XShm.h>

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glext.h>

#include "os.h"
#include "nox.h"
#include "glmesh.h"

#ifdef MEM
#include "memory.h"
#endif

#define ROTSTEP 180.0/12.0
//#define ROTSTEP	15.0
//#define ROTSTEP	1.0

//#define _DRAW_2D
#define _DRAW_3D
//#define _DRAW_BBOX

#define XK_STATE_SHIFT	0x01
#define XK_STATE_CTRL	0x04
#define XK_STATE_ALT	0x08

using namespace std;

static Display* display;
static int	depth;
static Window	window;
static GC	gc;
static unsigned black, white;
static XImage*	ximage;
static XVisualInfo* visinfo;
static Pixmap   pixmap;
static GLXContext context;
static GLXPixmap glxpixmap;

NoxGeometry nox;

static int sngBuf[] = {
		GLX_RGBA,
		GLX_RED_SIZE,    8,
		GLX_GREEN_SIZE,  8,
		GLX_BLUE_SIZE,   8,
		GLX_DEPTH_SIZE, 12,
//		GLX_DOUBLEBUFFER,
		None };
//static int dblBuf[] = {
//		GLX_RGBA,
//		GLX_RED_SIZE,    1,
//		GLX_GREEN_SIZE,  1,
//		GLX_BLUE_SIZE,   1,
////		GLX_DEPTH_SIZE, 12,
//		GLX_DOUBLEBUFFER,
//		None };

typedef struct {
	GLfloat x, y, z;
	GLfloat r, g, b, a;
} GLVertex;

/** drawPoint */
void drawPoint(GLVertex v1, GLfloat size)
{
	glPointSize(size);
	glBegin(GL_POINTS);
		glColor4f(v1.r, v1.g, v1.b, v1.a);
		glVertex3f(v1.x, v1.y, v1.z);
	glEnd();
} // drawPoint

/** drawPointsDemo */
void drawPointsDemo(int /*width*/, int /*height*/)
{
	GLfloat size = 5.0f;
	for (GLfloat x = 0.0f; x <= 1.0f; x += 0.2f, size += 5) {
		GLVertex v1 = { x, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f };
		drawPoint(v1, size);
	}
} // drawPointsDemo

/** drawLineSegment */
void drawLineSegment(GLVertex v1, GLVertex v2, GLfloat width)
{
	glLineWidth(width);
	glBegin(GL_LINES);
		glColor4f(v1.r, v1.g, v1.b, v1.a);
		glVertex3f(v1.x, v1.y, v1.z);
		glColor4f(v2.r, v2.g, v2.b, v2.a);
		glVertex3f(v2.x, v2.y, v2.z);
	glEnd();
} // drawLineSegment

/** drawGrid */
void drawGrid(GLfloat width, GLfloat height, GLfloat grid_width)
{
	//horizontal lines
	for(float i=-height; i<height; i+=grid_width){
		GLVertex v1 = {-width, i, 0.0f, 0.8f, 0.8f, 0.8f, 1.0f};
		GLVertex v2 = { width, i, 0.0f, 0.8f, 0.8f, 0.8f, 1.0f};
		drawLineSegment(v1, v2, 0.0f);
	}

	//vertical lines
	for(float i=-width; i<width; i+=grid_width){
		GLVertex v1 = {i, -height, 0.0f, 0.8f, 0.8f, 0.8f, 1.0f};
		GLVertex v2 = {i,  height, 0.0f, 0.8f, 0.8f, 0.8f, 1.0f};
		drawLineSegment(v1, v2, 0.0f);
	}
} // drawGrid

/** drawLineDemo */
void drawLineDemo()
{
	//draw a simple grid
	drawGrid(5.0f, 1.0f, 0.1f);
	//define the vertices and colors of the line segments
	GLVertex v1 = {-5.0f,  0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.7f};
	GLVertex v2 = { 5.0f,  0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.7f};
	GLVertex v3 = { 0.0f,  1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.7f};
	GLVertex v4 = { 0.0f, -1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.7f};
	//draw the line segments
	drawLineSegment(v1, v2, 0.1f);
	drawLineSegment(v3, v4, 0.1f);
} // drawLineDemo

/** drawTriangle */
void drawTriangle(GLVertex v1, GLVertex v2, GLVertex v3)
{
	glBegin(GL_TRIANGLES);
	glColor4f( v1.r, v1.g, v1.b, v1.a);
	glVertex3f(v1.x, v1.y, v1.z);
	glColor4f( v2.r, v2.g, v2.b, v2.a);
	glVertex3f(v2.x, v2.y, v2.z);
	glColor4f( v3.r, v3.g, v3.b, v3.a);
	glVertex3f(v3.x, v3.y, v3.z);
	glEnd();
} // drawTriangle

/** drawTriangleDemo */
void drawTriangleDemo()
{
	//Triangle Demo
	GLVertex v1 = { 0.0f,  0.8f, 0.0f, 1.0f, 0.0f, 0.0f, 0.6f};
	GLVertex v2 = {-1.0f, -0.8f, 0.0f, 0.0f, 1.0f, 0.0f, 0.6f};
	GLVertex v3 = { 1.0f, -0.8f, 0.0f, 0.0f, 0.0f, 1.0f, 0.6f};
	drawTriangle(v1, v2, v3);
} // drawTriangleDemo

/** initX */
void initX()
{
	/* use the information from the environment variable DISPLAY
	   to create the X connection:
	*/
	display = XOpenDisplay(NULL);
	int screen     = DefaultScreen(display);
	Visual* visual = DefaultVisual(display, screen);
	depth   = DefaultDepth(display, screen);
	black   = BlackPixel(display, screen),	/* get color black */
	white   = WhitePixel(display, screen);  /* get color white */
	assert(depth>16);
//	visinfo = glXChooseVisual(display, screen, dblBuf);
	visinfo = glXChooseVisual(display, screen, sngBuf);

	// Last argument is Direct rendering through False=through X11
	// Newer X-servers disable the indirect rendering!
	context = glXCreateContext(display, visinfo, None, True);

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
	values.function   = GXcopy;
	values.plane_mask = AllPlanes;
	values.foreground = black;
	values.background = white;
	gc = XCreateGC(display, window,
		GCFunction | GCPlaneMask | GCForeground | GCBackground,
		&values);

	ximage = XCreateImage(display, visual, depth, ZPixmap, 0,
			(char*)nox.viewer.painter.data(),
			nox.viewer.width(), nox.viewer.height(), 32, 0);
	pixmap = XCreatePixmap(display, window, ximage->width, ximage->height, depth);
	glxpixmap = glXCreateGLXPixmap(display, visinfo, pixmap);
	glXMakeCurrent(display, glxpixmap, context);
} // initX

/** draw */
void draw()
{
	ViewPort& view = nox.viewer.view();

	XPutImage(display, pixmap, gc, ximage, 0, 0, 0, 0, ximage->width, ximage->height);
	//XShmPutImage(display, pixmap, gc, ximage, 0, 0, 0, 0, width, height, False);

	glViewport(0, 0, view.width(), view.height());
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);

	/* remove back faces */
//	glDisable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

	glShadeModel(GL_SMOOTH);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);

	/* light */
	static GLfloat light_pos[4] = {5.0, 20.0, 10.0, 1.0 };
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
	glEnable(GL_LIGHT0);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	/* Setup projection */
	glMatrixMode(GL_PROJECTION);

//	glLoadIdentity();
//	float aspect = (float)ximage->width / (float)ximage->height;
//	glOrtho(-aspect, aspect, -1.f, 1.f, 1.f, -1.f);
//	drawLineDemo();
//	//glRotatef((float)glfwGetTime() * 50.f, 0.f, 0.f, 1.f);
//	drawPointsDemo(view.width(), view.height());
//	drawTriangleDemo();

	glLoadIdentity();
	float dw = view.imageWidth();
	float dh = view.imageHeight();
	glOrtho(-dw/2., dw/2., -dh/2, dh/2, 0.f, 1000.f);

	double matrix[4][4];
	memset(matrix, 0, sizeof(matrix));
	for (int j=0; j<3; j++)
		for (int i=0; i<3; i++)
			matrix[j][i] = view.matrix(j,i);
	matrix[3][3] = 1.0;
	glMultMatrixd((GLdouble*)&matrix);

	double x,y,z;
	view.origin(&x, &y, &z);

	glTranslated(-x,-y,-z);

	/* Setup geometry */
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Scan geometry and plot all meshes
	for (int i=2; i<nox.geometry.bodies.size(); i++) {
		GBody* body = nox.geometry.bodies[i];
		cout << *body << endl;
		GLMesh glmesh(&body->mesh);
		glmesh.draw(GL_RENDER);
	}

//	usleep(1000);	// FIXME stupid delay to wait to finish drawing
			// For some strange reason the Flush/Finish, nothing works

	glFlush();	// flush contents. Doesn't ensure that the graph is completed
//	glFinish();	// wait to finish and then flush
//	glXWaitGL();	// wait X11 to finish
//	glXSwapBuffers(display, glxpixmap);

	// ---------- Display the pixmap ----------
	XCopyArea(display, pixmap, window, gc, 0, 0, ximage->width, ximage->height, 0, 0);
//
} // draw

/** main */
int main(int ac, char *av[])
{
	Matrix4 mat;
	mat.identity();
//	nox.viewer.kernel().initThreads(0);

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

	ViewPort& view = nox.viewer.view();

	view.calcWindow(4.0);
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
		KeySym keysym;		/* a dealie-bob to handle KeyPress Events */
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
			if (width!=view.width() || height!=view.height()) {
				nox.viewer.resize(width, height);
				ximage->width  = view.width();
				ximage->height = view.height();
				ximage->data   = (char*)nox.viewer.painter.data();
				ximage->bytes_per_line = 0;
				XInitImage(ximage);
				nox.viewer.draw(mask, true);
				//XClearWindow(display, window);
				glXDestroyGLXPixmap(display, glxpixmap);
				XFreePixmap(display, pixmap);
				pixmap = XCreatePixmap(display, window, ximage->width, ximage->height, depth);
				glxpixmap = glXCreateGLXPixmap(display, visinfo, pixmap);
				glXMakeCurrent(display, glxpixmap, context);
			}
			draw();
		} else
		if (event.type==KeyPress) {
			// ---------- Display the pixmap ----------
			int len = XLookupString(&event.xkey,text,255,&keysym,0);
			/* use the XLookupString routine to convert the
			   KeyPress data into regular text.  Weird but necessary...
			*/
			//printf("key=%d sym=%d state=%d text=%s\n", event.xkey.keycode, keysym, event.xkey.state, text);

			double u = view.Uofs();	// old x
			double v = view.Vofs();	// old y
			switch (keysym) {
				case XK_equal:
					view.zoom(view.zoom()*1.2);
					cout << "Zoom=" << view.zoom() << endl;
					if (view.zoom() > 1000.0) {
						view.calcWindow(4.0);
						nox.viewer.d2.project();
					}
					break;

				case XK_minus:
					view.zoom(view.zoom()/1.2);
					cout << "Zoom=" << view.zoom() << endl;
					break;

				case XK_h:
				case XK_Left:
					if (event.xkey.state & XK_STATE_CTRL) {
						mat.rotY(RAD(-ROTSTEP));
						view.transform(mat);
						nox.viewer.d2.project();
					} else {
						if (event.xkey.state & XK_STATE_SHIFT)
							u -= (view.width()/100) / view.Sx();
						else
							u -= (view.width()/10) / view.Sx();
						view.offset(u,v);
					}
					break;

				case XK_Right:
					if (event.xkey.state & XK_STATE_CTRL) {
						mat.rotY(RAD(ROTSTEP));
						view.transform(mat);
						nox.viewer.d2.project();
					} else {
						if (event.xkey.state & XK_STATE_SHIFT)
							u += (view.width()/100) / view.Sx();
						else
							u += (view.width()/10) / view.Sx();
						view.offset(u,v);
					}
					break;

				case XK_Up:
					if (event.xkey.state & XK_STATE_CTRL) {
						mat.rotX(RAD(-ROTSTEP));
						view.transform(mat);
						nox.viewer.d2.project();
					} else {
						if (event.xkey.state & XK_STATE_SHIFT)
							v += (view.height()/100) / view.Sy();
						else
							v += (view.height()/10) / view.Sy();
						view.offset(u,v);
					}
					break;

				case XK_Down:
					if (event.xkey.state & XK_STATE_CTRL) {
						mat.rotX(RAD(ROTSTEP));
						view.transform(mat);
						nox.viewer.d2.project();
					} else {
						if (event.xkey.state & XK_STATE_SHIFT)
							v -= (view.height()/100) / view.Sy();
						else
							v -= (view.height()/10) / view.Sy();
						view.offset(u,v);
					}
					break;

				case XK_Page_Down:
					if (event.xkey.state & XK_STATE_CTRL) {
						mat.rotZ(RAD(-ROTSTEP));
					} else {
						mat.translate(0.0, 0.0, (view.width()/50) / view.Sx());
					}
					view.transform(mat);
					nox.viewer.d2.project();
					break;

				case XK_Page_Up:
					if (event.xkey.state & XK_STATE_CTRL) {
						mat.rotZ(RAD(ROTSTEP));
					} else {
						mat.translate(0.0, 0.0, -(view.width()/50) / view.Sx());
					}
					view.transform(mat);
					nox.viewer.d2.project();
					break;


				case XK_l:
				case XK_r:
					if (event.xkey.state & XK_STATE_CTRL) // redraw
						nox.viewer.draw(mask, true);
					break;

				case XK_H:
					u -= (view.width()/100) / view.Sx();
					view.offset(u,v);
					break;

				case XK_J:
					v -= (view.height()/100) / view.Sy();
					view.offset(u,v);
					break;

				case XK_K:
					v += (view.height()/100) / view.Sy();
					view.offset(u,v);
					break;

				case XK_L:
					u += (view.width()/100) / view.Sx();
					view.offset(u,v);
					break;

				// check image shift
//				case XK_s:
//					nox.viewer.painter->move(-50, 0);
//					nox.viewer.painter->setClip(-51, 0, -1, -1);
//					x += (int)(50.0 / view.Sx());
//					view.setOrigin(x,y);
//					break;
//
//				case XK_d:
//					nox.viewer.painter->move(50, 0);
//					nox.viewer.painter->setClip(50, 0, -1, -1);
//					break;

//				case XK_e:
//					nox.viewer.painter->move(0, -50);
//					nox.viewer.painter->setClip(0, -50, -1, -1);
//					break;

//				case XK_x:
//					nox.viewer.painter->move(0, 50);
//					nox.viewer.painter->setClip(0, 50, -1, -1);
//					break;

				case XK_c:
					view.offset(0.,0.);
					break;

				case XK_z: // zoom test
					cout << "zoom test\n" ;
					cout << "view=" << view << endl;
					view.calcWindow(2.0);
					cout << "\nafter\n";
					cout << "view=" << view << endl;
					nox.viewer.d2.project();
					break;

				case XK_q:
					nox.viewer.stopThread();
					goto END;

				case XK_b:
					nox.viewer.stopThread();
					break;

				default:
					printf("You pressed: ");
					for (int i=0; i<len; i++)
						printf("%02X ",text[i]);
					printf("\n");
			}
			if (view.invalidWindow()) {
				view.calcWindow(4.0);
				cout << "Invalid window" << endl;
				nox.viewer.d2.project();
				cout << "Transform" << endl;
			}
//			nox.viewer.draw(mask, true);
			cout << view << endl;
			cout << "BasisX: " << view.matrix(0,0)
					   << " " << view.matrix(1,0)
					   << " " << view.matrix(2,0) << endl;
			cout << "BasisY: " << view.matrix(0,1)
					   << " " << view.matrix(1,1)
					   << " " << view.matrix(2,1) << endl;

			//view.moveOriginTo0();

			nox.viewer.painter.resetClip();
			// ---------- Display the pixmap ----------
			draw();
		} else
		if (event.type==ButtonPress) {
			/* tell where the mouse Button was Pressed */
			double x, y, z;
			view.ij2xyz(event.xbutton.x, event.xbutton.y, &x, &y, &z);
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
	glXDestroyGLXPixmap(display, glxpixmap);
	XFreePixmap(display, pixmap);
	XFreeGC(display, gc);
	XDestroyWindow(display,window);
	XCloseDisplay(display);
	return 0;
} // main
