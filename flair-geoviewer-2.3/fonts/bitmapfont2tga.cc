/**
 * $Id$
 *
 * Author:	Vasilis.Vlachoudis@cern.ch
 * Date:	12-Dec-2012
 */

#include <os.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>

int main(int ac, char *av[])
{
	if (ac!=3) {
		printf("syntax: %s <font-name> <filename.tga>\n",av[0]);
		return 1;
	}

	Display* display = XOpenDisplay(NULL);
	int screen = DefaultScreen(display);

	XFontStruct* fontinfo = XLoadQueryFont(display, av[1]);
	if (fontinfo == NULL) {
		printf("Invalid font \'%s\' specified", av[1]);
		return 2;
	}
	int width  = fontinfo->max_bounds.width;
	int height = fontinfo->max_bounds.ascent + fontinfo->max_bounds.descent;

	printf("Font= %s\n", av[1]);
	printf("File= %s\n", av[2]);
	printf("Size= %dx%d\n",16*width, 16*height);

	Pixmap pixmap = XCreatePixmap(display,
				XRootWindow(display,screen),
				16*width, 16*height,
				DefaultDepth(display,screen));

	XGCValues values;
	values.font       = fontinfo->fid;
	values.function   = GXcopy;
	values.plane_mask = AllPlanes;
	values.foreground = BlackPixel(display,screen);
	values.background = BlackPixel(display,screen);

	GC gc = XCreateGC(display,
		pixmap,
		GCFont |
		GCFunction | GCPlaneMask | GCForeground | GCBackground,
		&values);

	XFillRectangle(display,pixmap,gc,0,0,16*width,16*height);

	XSetForeground(display, gc, WhitePixel(display,screen));
	for (int row=0; row<16; row++) {
		int y = row*height + fontinfo->max_bounds.ascent;
		for (int col=0; col<16; col++) {
			int x = col*width;
			char str[2];
			str[0] = row<<4 | col;
			str[1] = 0;
			XDrawImageString(display,pixmap,gc,x,y,str,1);
		}
	}
	XFreeFont(display, fontinfo);

	XImage* image = XGetImage(display, pixmap,
		0, 0, 16*width, 16*height, AllPlanes,ZPixmap);

//	font.setMap(fontpath, 16*width, 16*height, (dword*)image->data);
	int *ptr = (int*)image->data;

	FILE *f = fopen(av[2],"wb");

	byte header[18];
	memset(header,0,sizeof(header));
	header[2] = 3;
	*(short*)(header+12) = 16*width;
	*(short*)(header+14) = 16*height;
	header[16] = 8;
	header[17] = 0x20;
	fwrite(header,1,sizeof(header),f);

	for (int i=0; i<16*16*width*height; i++, ptr++)
		fputc(*ptr?0xFF:0, f);

	fclose(f);
	XDestroyImage(image);
	XFreeGC(display, gc);
	XFreePixmap(display, pixmap);
	return 0;
} // main
