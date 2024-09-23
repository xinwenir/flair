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
 * Date:	26-Jul-2005
 */
#ifndef __USRBIN_H
#define __USRBIN_H

#include <sys/stat.h>
#include <string>

#include "os.h"
#include "matrix4.h"

#define USRBIN_RUNTIT	80
#define	USRBIN_RUNTIM	32
#define USRBIN_TITUSB	10

enum UsrbinType {
	Usrbin_XYZ_point      =  0,
	Usrbin_XYZ            = 10,

	Usrbin_RPhiZ_point    =  1,
	Usrbin_RPhiZ          = 11,

	Usrbin_region_point   =  2,
	Usrbin_region         = 12,

	Usrbin_Xsym_point     =  3,
	Usrbin_Xsym           = 13,

	Usrbin_Ysym_point     =  4,
	Usrbin_Ysym           = 14,

	Usrbin_Zsym_point     =  5,
	Usrbin_Zsym           = 15,

	Usrbin_XYZsym_point   =  6,
	Usrbin_XYZsym         = 16,

	Usrbin_RPhiZsym_point =  7,
	Usrbin_RPhiZsym       = 17,

	Usrbin_special_point  =  8,
	Usrbin_special        = 18
}; // UsrbinType

/*
 * If the standard system is used
 *	X = horizontal
 *	Y = vertical
 *	Z = horizontal-beam direction
 * Then the transformation that has to be performed are
 *	X -> -X		[X]  Y   Z	Size: -1.0
 *	Y ->  Z		 X   Y  [Z]	Size:  1.0
 *	Z ->  Y		 X  [Y]  Z	Size:  1.0
 */
class Usrbin {
private:
	int	_detector;		// detector loaded
	UsrbinType _type;		// type of usrbin
	int	_score;			// quantity scored
public:
	int	nx,    ny,    nz;
	double	xlow,  ylow,  zlow;
	double	xhigh, yhigh, zhigh;
	double	dx,    dy,    dz;
//	Point	offset;
	double	xofs,  yofs,  zofs;	// offset of usrbin
	double	x0,  y0;		// for R-Z-Phi binning
	double	min, max;
	char	runtitle[USRBIN_RUNTIT+1],
		runtime[USRBIN_RUNTIM+1];
	char	titusb[USRBIN_TITUSB+1];

	int	lntzer;			// not zero flag
	float	bkusbn;
	float	b2usbn;
	float	tcusbn;

	float	weipri;
	int	ncase;
	int	mcase;
	int	nbatch;

private:
	double	_norm;
	bool	_logscale;		// is logarithmic histogram
	bool	_checker;
	bool	_hasMatrix;		// if rotation matrix is set
	Matrix4	_matrix;

	int	nynx;
	time_t	_mtime;
	std::string	filename;

	float	*data;
//	float	*error;
//	double	 _total;

static	char*	_typeStr[];

public:
	Usrbin()			: data(NULL) { cleanup(); }
	Usrbin(const char* fn, int det)	: data(NULL) { load(fn,det); }
	~Usrbin()	{ cleanup(); }

	void	cleanup();
	bool	load(const char* fn, int det);
	bool	save(const char* fn=NULL);

	/* create usrbin */
	void	create(int det, UsrbinType atype, int ascore,
			const double axlow, const double axhigh, const int anx,
			const double aylow, const double ayhigh, const int any,
			const double azlow, const double azhigh, const int anz);
	void	checker(const UsrbinType atype,
			double axlow, double axhigh, const int anx,
			double aylow, double ayhigh, const int any,
			double azlow, double azhigh, const int anz);
	bool	checker()	const	{ return _checker; }

	/* get/set */
	int	detector()	const	{ return _detector; }
	UsrbinType type()	const	{ return _type; }
	int	score()		const	{ return _score; }

	bool	hasData()	const	{ return data != NULL; }

	double	xmin()		const	{
			if (_type==Usrbin_XYZsym_point || _type==Usrbin_XYZsym ||
			    _type==Usrbin_Xsym_point   || _type==Usrbin_Xsym)
				return -xhigh;
			else
				return xlow;
		}
	double	ymin()		const	{
			if (_type==Usrbin_XYZsym_point || _type==Usrbin_XYZsym ||
			    _type==Usrbin_Ysym_point   || _type==Usrbin_Ysym)
				return -yhigh;
			else
				return ylow;
		}
	double	zmin()		const	{
			if (_type==Usrbin_XYZsym_point   || _type==Usrbin_XYZsym ||
			    _type==Usrbin_Zsym_point     || _type==Usrbin_Zsym   ||
			    _type==Usrbin_RPhiZsym_point || _type==Usrbin_RPhiZsym)
				return -zhigh;
			else
				return zlow;
		}
	double	xmax()			const	{ return xhigh; }
	double	ymax()			const	{ return yhigh; }
	double	zmax()			const	{ return zhigh; }

	void	norm(const double n);
	double	norm()			const	{ return _norm; }
	void	normalize(const double n);

	bool	convert(bool tolog);
	bool	logscale()		const	{ return _logscale; }

	// WARNING no transformation is applied!!!
	int	usrbini(const double x)	const	{ return Int((x - xlow)/dx); }
	int	usrbinj(const double y)	const	{ return Int((y - ylow)/dy); }
	int	usrbink(const double z)	const	{ return Int((z - zlow)/dz); }

	double	usrbinx(const int i)	const	{ return xlow + i*dx; }
	double	usrbiny(const int j)	const	{ return ylow + j*dy; }
	double	usrbinz(const int k)	const	{ return zlow + k*dz; }

	double	usrbincx(const int i)	const	{ return xlow + i*dx + dx/2.0; }
	double	usrbincy(const int j)	const	{ return ylow + j*dy + dy/2.0; }
	double	usrbincz(const int k)	const	{ return zlow + k*dz + dz/2.0; }

	// with transformation
	bool	xyz2ijk(double x, double y, double z, int* i, int* j, int* k) const;

	double	getData(double x, double y, double z, bool *ok) const;
	double	getData(const int region, bool *ok) const;
	double	getDataRPhiZ(double r, double phi, double z, bool *ok) const;
	double	get(const int ptr) const	{ return data[ptr]; }
	double	get(int i, int j, int k) const	{ return (double)data[k*nynx + j*nx + i]; }

	/** @return real value of usrbin at location */
	double	get(double x, double y, double z, bool *ok) const {
			double v = getData(x,y,z,ok);
			if (*ok && logscale())
				v = pow(10.0,v);
			return v;
		}
	/** @return real value of usrbin for the region */
	double	get(const int region, bool *ok) const {
			double v = getData(region, ok);
			if (*ok && logscale())
				v = pow(10.0,v);
			return v;
		}

	/** set a value to a cell */
	void	set(int i, int j, int k, double value) {
			assert(!_logscale && !_checker && data);
			data[k*nynx + j*nx + i] = value;
		}

	/** increase cell value by amount v */
	void	add(int i, int j, int k, double value) {
			assert(!_logscale && !_checker && data);
			data[k*nynx + j*nx + i] += value;
			/*_total += value;*/
		}
	void	add(double x, double y, double z, double value) {
			assert(!_logscale && !_checker && data);
			int i,j,k;
			if (xyz2ijk(x,y,z, &i, &j, &k)) add(i,j,k,value);
		}
	double	add(Point start, Point stop, double value=1.0, bool tracklength=true);

	void	offset(const Point& ofs);
	void	scanMinMax();

	void	clearMatrix()			{ _hasMatrix = false; }
	bool	hasMatrix()		const	{ return _hasMatrix; }
	void	matrix(const Matrix4& M)	{ _matrix.copy(M); _hasMatrix=true; }
const	Matrix4& matrix()		const	{ return _matrix; }

private:
	// returned normalized value
	double	_normValue(double value) const {
			if (_logscale)
				return value + _norm;
			else
				return value * _norm;
		}
}; // Usrbin

#endif
