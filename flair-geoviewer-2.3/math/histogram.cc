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
 * Date:	02-Feb-2015
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "os.h"
#include "bmath.h"
#include "histogram.h"

/** Histogram */
Histogram::Histogram(const char *atitle, int bins, double low, double high)
	: xbins(0), xlow(0.0), xhigh(1.0), h(NULL), eh(NULL)
{
	if (atitle) title(atitle);
	set(bins, low, high);
} // Histogram

/** Histogram */
Histogram::Histogram(const Histogram& hist)
	: xbins(0), h(NULL), eh(NULL)
{
	if (hist.title()) title(hist.title());
	set(hist.xbins, hist.xlow, hist.xhigh);
	memcpy(h, hist.h, sizeof(double)*xbins);
	memcpy(eh, hist.eh, sizeof(double)*xbins);

	_entries = hist._entries;
	_under   = hist._under;
	_over    = hist._over;
	_total   = hist._total;
	_sum     = hist._sum;
} // Histogram

/** ~Histogram */
Histogram::~Histogram()
{
	free();
} // ~Histogram

/** free */
void Histogram::free()
{
	if (h)  delete [] h;
	if (eh) delete [] eh;
	h  = NULL;
	eh = NULL;
	xbins = 0;
} // free

/** set */
void Histogram::set(int n, double low, double high)
{
	if (n != xbins) {
		free();
		xbins = n;
	}
	setLimits(low,high);

	if (xbins && h==NULL) {
		h  = new double[xbins];
		eh = new double[xbins];
	}
	clear();
} // set

/** setLimits */
void Histogram::setLimits(double low, double high)
{
	xlow  = low;
	xhigh = high;
	if (xbins) xstep = (xhigh-xlow) / (double)xbins;
} // setLimits

/** clear */
void Histogram::clear()
{
	_entries = 0;
	_under   = 0.0;
	_over    = 0.0;
	_total   = 0.0;
	_sum     = 0.0;
	memset(h,  0, sizeof(h[0]) *xbins);
	memset(eh, 0, sizeof(eh[0])*xbins);
} // clear

/** set */
void Histogram::set(int bin, double content)
{
	_entries++;
	if (bin<0) {
		_total -= _under;
		_under  = content;
		_total += content;
	} else
	if (bin>=xbins) {
		_total -= _over;
		_over   = content;
		_total += content;
	} else {
		_total -= h[bin];
		h[bin]  = content;
		_total += content;
	}
} // set

/** fill histogram */
int Histogram::fill(const double x, const double w)
{
	_entries++;
	_total += w;
	_sum   += x*w;
	if (x<xlow) {
		_under += w;
		return -1;
	} else
	if (x>xhigh) {
		_over  += w;
		return -1;
	} else {
		h[geti(x)] += w;
		return geti(x);
	}
} // fill

/** normalize histogram with f value
 * @param f normalization factor
 */
void Histogram::norm(const double f)
{
	assert(xbins);
	_under *= f;
	_over  *= f;
	_total *= f;
	_sum   *= f;
	for (int i=0; i<xbins; i++) {
		h[i]  *= f;
		eh[i] *= f;
	}
} // norm

/** convert the histogram to cumulative one loop: h[i+1] += h[i]
 *  ignoring the _over and _under!
 * @param reverse from max to 0.0
 */
void Histogram::cumulative(bool reverse)
{
	assert(xbins);
	if (!reverse) {
		h[0] += _under;
		for (int i=1; i<xbins; i++)
			h[i] += h[i-1];
			// errors?
		_over += h[xbins-1];
	} else {
		h[xbins-1] += _over;
		for (int i=xbins-2; i>=0; i--)
			h[i] += h[i+1];
			// errors?
		_under += h[0];
	}
} // cumulative

/** convert a log10(histogram) to isolethargic */
void Histogram::isolethargic()
{
	norm(pow(10.0,xstep/2.0) / pow(10.0,xstep-1.0));
} // isolethargic

/** save */
void Histogram::save(FILE *fout, int columns, bool lastdup)
{
	fprintf(fout,"# title   %s\n",title());
	fprintf(fout,"# xmin    %.10g\n",xlow);
	fprintf(fout,"# xmax    %.10g\n",xhigh);
	fprintf(fout,"# xbins   %d\n"   ,xbins);
	fprintf(fout,"# xstep   %.10g\n",xstep);
	fprintf(fout,"# n       %d\n"   ,_entries);
	fprintf(fout,"# weight  %.10g\n",_total);
	fprintf(fout,"# under   %.10g\n",_under);
	fprintf(fout,"# over    %.10g\n",_over);
	fprintf(fout,"# total   %.10g\n",_total);
	fprintf(fout,"# sum     %.10g\n",_sum);

	double x = xlow;
	for (int i=0; i<xbins; i++) {
		double xh = x + xstep;
		switch (columns) {
			case 3:
				fprintf(fout,"%.10g %.10g %.10g\n",x,h[i],eh[i]);
				break;

			case 4:
				fprintf(fout,"%.10g %.10g %.10g %.10g\n",x,xh,h[i],eh[i]);
				break;

			default:	/* 2 */
				fprintf(fout,"%.10g %.10g\n",x,h[i]);
		}
		x = xh;
	}

	/* correction for gnuplot. Duplicate the last line */
	if (lastdup) {
		int i = xbins-1;
		switch (columns) {
			case 3:
				fprintf(fout,"%.10g %.10g %.10g\n",x,h[i],eh[i]);
				break;

			case 4:
				fprintf(fout,"%.10g %.10g %.10g %.10g\n",x,x,h[i],eh[i]);
				break;

			default:	/* 2 */
				fprintf(fout,"%.10g %.10g\n",x,h[i]);
		}
	}
} // save

/** load */
bool Histogram::load(FILE *)
{
	assert(0);
	return true;
} // load
