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

#ifndef __HISTOGRAM_H
#define __HISTOGRAM_H

#include <string>

class Histogram {
private:
	std::string _title;	/** histogram title		*/

	int	xbins;		/** number of bins		*/
	double	xlow;		/** lower x limit		*/
	double	xhigh;		/** higher x limit		*/
	double	xstep;		/** x step size			*/

	int	_entries;	/** number of entries		*/
	double	_under;		/** weights below low x limit	*/
	double	_over;		/** weights above high x limit	*/
	double	_total;		/** total sum of weights	*/
	double	_sum;		/** integral			*/

	double	*h;		/** histogram values		*/
	double	*eh;		/** errors			*/

public:
	Histogram(const char *atitle=NULL, int bins=0, double low=0.0, double high=1.0);
	Histogram(const Histogram& hist);
	~Histogram();

	/* initialize */
	void	free();
	void	clear();

	/* get/set */
const	char*	title()		const	{ return _title.c_str(); }
//	std::string title()	const	{ return _title; }
	void	title(const char *t)	{ _title = t; }
	void	title(const std::string& t)	{ _title = t; }

	void	set(int n, double low, double high);
	void	setLimits(double low, double high);

	int	nbins()		const	{ return xbins; }
//	int	nbinsx()	const	{ return xbins; }
	int	entries()	const	{ return _entries; }
	int	under()		const	{ return _under; }
	int	over()		const	{ return _over; }
	int	total()		const	{ return _total; }

	double	width()		const	{ return xstep; }
	double	mean()		const	{ return _sum / _total; }

	/* fill histogram */
	int	geti(const double x)	const	{ return (int)((x - xlow) / xstep); }
	double	getx(const int i)	const	{ return (double)i*xstep + xlow; }
	double	center(const int i)	const	{ return xlow + i*xstep + 0.5*xstep; }
	int	fill(const double x, const double w=1.0);

	void	set(int bin, double content);
	double	get(int i) const	{ assert(i>=0 && i<xbins); return h[i]; }
	double	operator[](int i) const	{ assert(i>=0 && i<xbins); return h[i]; }
	double&	operator[](int i)	{ assert(i>=0 && i<xbins); return h[i]; }

	void	norm(const double f);
	void	scale(const double f)	{ norm(f); }
	void	isolethargic();
	void	cumulative(bool reverse=false);

	/* I/O */
	void	save(FILE *fout, int columns=4, bool lastdup=false);
	bool	load(FILE *fin);
}; // Histogram

class H1D : public Histogram {
public:
	H1D(const char *atitle=NULL, int bins=0, double low=0.0, double high=1.0)
		: Histogram(atitle, bins, low, high) {}
	H1D(const H1D& hist) : Histogram(hist) {}
}; // class H1D

#endif
