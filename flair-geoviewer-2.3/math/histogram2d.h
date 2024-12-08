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

#ifndef __HISTOGRAM2D_H
#define __HISTOGRAM2D_H

#include <string>
#include <assert.h>

class Histogram2D {
public:
	std::string _title;	/** histogram title		*/

	int	xbins;		/** number of x bins		*/
	double	xlow;		/** lower x limit		*/
	double	xhigh;		/** higher x limit		*/
	double	xstep;		/** x step size			*/

	int	ybins;		/** number of y bins		*/
	double	ylow;		/** lower y limit		*/
	double	yhigh;		/** higher y limit		*/
	double	ystep;		/** y step size			*/

	int	_entries;	/** number of entries		*/
	double	_total;		/** total sum of weights	*/
	double	_sum;		/** integral			*/

	double	*h;		/** histogram values		*/
	double	*eh;		/** errors			*/

public:
	Histogram2D(const char *atitle=NULL,
				int nx=0, double lowx=0.0, double highx=1.0,
				int ny=0, double lowy=0.0, double highy=1.0);
	Histogram2D(const Histogram2D& hist);
	~Histogram2D()	{ free(); }

	/* initialize */
	void	free();
	void	clear();
	void	copy(const Histogram2D& src);

	/* get/set */
const	char*	title()		const	{ return _title.c_str(); }
//	std::string title()	const	{ return _title; }
	void	title(const char *t)	{ _title = t; }
	void	title(const std::string& t)	{ _title = t; }

	void	set(int nx=0, double lowx=0.0, double highx=1.0,
		    int ny=0, double lowy=0.0, double highy=1.0);
	void	setXLimits(double low, double high);
	void	setYLimits(double low, double high);

	int	entries()	const	{ return _entries; }
	int	total()		const	{ return _total; }

//	double	width()		const	{ return xstep; }
//	double	mean()		const	{ return _sum / _total; }

//	/* fill histogram */
//	int	geti(const double x)	const	{ return (int)((x - xlow) / xstep); }
//	double	getx(const int i)	const	{ return (double)i*xstep + xlow; }
//	double	center(const int i)	const	{ return xlow + i*xstep + 0.5*xstep; }
//	int	fill(const double x, const double w=1.0);

//	void	set(int bin, double content);
	double	get(int i) const	{ assert(i>=0 && i<xbins*ybins); return h[i]; }
	double	operator[](int i) const	{ assert(i>=0 && i<xbins*ybins); return h[i]; }
	double&	operator[](int i)	{ assert(i>=0 && i<xbins*ybins); return h[i]; }

//	void	norm(const double f);
//	void	scale(const double f)	{ norm(f); }
//	void	isolethargic();
//	void	cumulative(bool reverse=false);

//	/* I/O */
	bool	save(FILE *fout);
	bool	save(const char* filename) {
			FILE *f = fopen(filename,"w");
			if (f) {
				bool rc = save(f);
				fclose(f);
				return rc;
			} else
				return false;
		}
	bool	load(FILE *fin);
	bool	load(const char* filename) {
			FILE *f = fopen(filename,"r");
			if (f) {
				bool rc = load(f);
				fclose(f);
				return rc;
			} else
				return false;
		}
}; // Histogram2D

class H2D : public Histogram2D {
public:
	H2D(const char *atitle=NULL,
		int nx=0, double lowx=0.0, double highx=1.0,
		int ny=0, double lowy=0.0, double highy=1.0)
		: Histogram2D(atitle, nx,lowx,highx, ny,lowy,highy) {}
	H2D(const H2D& hist) : Histogram2D(hist) {}
}; // class H2D

#endif
