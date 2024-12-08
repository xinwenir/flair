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
 * Date:	27-Mar-2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "histogram2d.h"

/** Histogram2D */
Histogram2D::Histogram2D(const char *atitle,
				int nx, double lowx, double highx,
				int ny, double lowy, double highy)
	: xbins(0), xlow(0.0), xhigh(1.0),
	  ybins(0), ylow(0.0), yhigh(1.0),
	  h(NULL), eh(NULL)
{
	if (atitle) title(atitle);
	set( nx, lowx, highx, ny, lowy, highy);
} // Histogram2D

/** Histogram2D */
Histogram2D::Histogram2D(const Histogram2D& hist)
	: xbins(0), xlow(0.0), xhigh(1.0),
	  ybins(0), ylow(0.0), yhigh(1.0),
	  h(NULL), eh(NULL)
{
	copy(hist);
} // Histogram2D

/** copy */
void Histogram2D::copy(const Histogram2D& hist)
{
	free();
	if (hist.title()) title(hist.title());
	set(hist.xbins, hist.xlow, hist.xhigh,
	    hist.ybins, hist.ylow, hist.yhigh);
	memcpy(h, hist.h, sizeof(double)*xbins*ybins);
	memcpy(eh, hist.eh, sizeof(double)*xbins*ybins);

	_entries = hist._entries;
	_total   = hist._total;
	_sum     = hist._sum;
} // copy

/** free */
void Histogram2D::free()
{
	if (h)  delete [] h;
	if (eh) delete [] eh;
	h  = NULL;
	eh = NULL;
	xbins = ybins = 0;
} // free

/** set */
void Histogram2D::set(	int nx, double lowx, double highx,
			int ny, double lowy, double highy)
{
	if (xbins != nx || ybins != ny) {
		free();
		xbins = nx;
		ybins = ny;
	}
	setXLimits(lowx, highx);
	setYLimits(lowy, highy);

	if (xbins && ybins && h==NULL) {
		h  = new double[xbins*ybins];
		eh = new double[xbins*ybins];
	}
	clear();
} // set

/** clear */
void Histogram2D::clear()
{
	_entries = 0;
	_total   = 0.0;
	_sum     = 0.0;
	memset(h,  0, sizeof(h[0]) *xbins*ybins);
	memset(eh, 0, sizeof(eh[0])*xbins*ybins);
} // clear

/** setXLimits */
void Histogram2D::setXLimits(double low, double high)
{
	xlow  = low;
	xhigh = high;
	if (xbins) xstep = (xhigh-xlow) / (double)xbins;
} // setXLimits

/** setYLimits */
void Histogram2D::setYLimits(double low, double high)
{
	ylow  = low;
	yhigh = high;
	if (ybins) ystep = (yhigh-ylow) / (double)ybins;
} // setYLimits

/** save */
bool Histogram2D::save(FILE* fout)
{
	fprintf(fout, "#xbins %d\n",  xbins);
	fprintf(fout, "#xmin  %lg\n", xlow);
	fprintf(fout, "#xmax  %lg\n", xhigh);
	fprintf(fout, "#xstep %lg\n", xstep);

	fprintf(fout, "#ybins %d\n",  ybins);
	fprintf(fout, "#ymin  %lg\n", ylow);
	fprintf(fout, "#ymax  %lg\n", yhigh);
	fprintf(fout, "#ystep %lg\n", ystep);

	double x = xlow;
	for (int i=0, ptr=0; i<xbins; i++, x += xstep) {
		double y = ylow;
		for (int j=0; j<ybins; j++, ptr++, y += ystep)
			fprintf(fout, "%lg %lg %lg %lg\n", x, y, h[ptr], eh[ptr]);
		fprintf(fout, "\n");
	}
	return true;
} // save

/** load */
bool Histogram2D::load(FILE* fin)
{
	int    nx=0, ny=0;
	double lowx=0.0, lowy=1.0;
	double highx=0.0, highy=1.0;
	char	line[256];
	free();

	int ptr = 0;
	while (!feof(fin)) {
		if (fgets(line, sizeof(line), fin)==NULL) break;
		int l = strlen(line);
		// chop new line at end
		if (l && line[l-1]=='\n') line[--l] = 0;
		if (line[0]==0)
			continue;
		else
		if (line[0]=='#') {
			char var[16], value[128];
			sscanf(line+1,"%s %s",var,value);
			for (char *p=var; *p; ++p) *p = tolower(*p);
			if (!strcmp(var, "bins") || !strcmp(var, "xbins"))
				nx = atoi(value);
			else
			if (!strcmp(var, "xmin"))
				lowx = atof(value);
			else
			if (!strcmp(var, "xmax"))
				highx = atof(value);
			else
			if (!strcmp(var, "ybins"))
				ny = atoi(value);
			else
			if (!strcmp(var, "ymin"))
				lowy = atof(value);
			else
			if (!strcmp(var, "ymax"))
				highy = atof(value);
		} else {
			if (nx*ny==0)
				fprintf(stderr, "Cannot load histogram xbins=%d ybins=%d",nx,ny);
			else
			if (ptr==0)
				set(nx, lowx, highx, ny, lowy, highy);

			double x,y,v,e=0.0;
			sscanf(line,"%lg %lg %lg %lg",&x,&y,&v,&e);
			h[ptr] = v;
			eh[ptr] = e;
			ptr++;
		}
	}
	return ptr==xbins*ybins;
} // load
