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

#ifndef __VOXEL_H
#define __VOXEL_H

#include <string.h>
#include <assert.h>
#include <sys/stat.h>

#include "os.h"
#include "eps.h"
#include "array.h"
#include "color.h"
#include "matrix4.h"

class VVoxel;

/* ============================== ROICombination ============================ */
// Internal class holding the ROI combinations
class ROICombination {
public:
	int	 length;	// length of possible rt-structures
	short	*roi;		// structure list

public:
	ROICombination() : length(0), roi(NULL) {}
	~ROICombination()	{ if (roi) delete [] roi; }

	short	operator[](int i) const	{ assert(i<length); return roi[i]; }

	void	parse(void *buffer, int buflen);
	void	print()	{
			for (int i=0; i<length; i++)
				std::cout << ':' << roi[i];
			std::cout << std::endl;
		}
	size_t	memory() const	{ return sizeof(ROICombination) + length * sizeof(short); }
}; // class ROICombination

/* =============================== ROIPlanarSlice =========================== */
class ROIPlanarSlice {
public:
	int	 n;		// number of X,Y pairs
	Point2D *pts;		// points
public:
	ROIPlanarSlice(int an) : n(0), pts(NULL) { alloc(an); }
	~ROIPlanarSlice() { if (pts) delete [] pts; }

	void	alloc(int an) {
			if (an==n) return;
			n = an;
			if (pts) delete [] pts;
			pts = new Point2D[n];
		}
}; // ROIPlanarSlice

/* ================================ ROIPlanar ============================ */
class ROIPlanar {
public:
	int		 nz;
	ROIPlanarSlice	**slice;
public:
	ROIPlanar() : nz(0), slice(NULL) {};
	~ROIPlanar();

	void	alloc(int anz);
	void	parse(void *buffer, int buflen);
		operator int()	const	{ return slice==NULL; }

	ROIPlanarSlice&	addZslice(int z, int n);
}; // ROIPlanar

/* ================================== GVoxel ================================ */
/** GVoxel, voxel file structure */
class GVoxel {
public:
	char	title[84];
	int	no;			// number of unique regions
	int	mo;			// number of unique voxel/hu id (+1 for outer shell)
	int	nx,    ny,    nz;
	double	dx,    dy,    dz;
	double	xlow,  ylow,  zlow;	// after setting *low call the calcLimits()
	double	xhigh, yhigh, zhigh;	// ... to set the *high

	Array<std::string> cards;

private:
	word		*_data;		// Voxel id
	word		*_kreg;		// voxel to region conversion _kreg[_data[]]

	int		nynx;

	time_t		_mtime;
	std::string	filename;

	bool		_hasMatrix;	// if rotation matrix is set
	Matrix4		_matrix;	// rotation matrix

	// Structure information in the voxel
	int		 _nroi;		// number of roistructures [in the file 1-based, in memory 0-based]
	word		*_roi;		// structure memory equal to data memory
	std::string	*_roiName;	// array of roi structure names, size=nroi
	dword		*_roiColor;	// array of roi colors
	int		 _roiCombN;	// number of combinations
	int		 _roiMaxComb;	// maximum length of combination
	ROICombination	*_roiComb;	// double array of roiCombN x nroi
	int		 _roiPlanarN;	// number of closed planar structures
	ROIPlanar	*_roiPlanar;	// ROI structures in polygonal form

public:
	GVoxel();
	~GVoxel()		{ cleanup(); }

	bool	load(const char *fn, bool loadCards=false);
	bool	save(const char *fn=NULL);
	void	cleanup();
	void	calcLimits();

	int	get(const int ptr) const
			{ return (int)_kreg[_data[ptr]]; }
	int	get(const int i, const int j, const int k) const
			{ return (int)_kreg[_data[k*nynx + j*nx + i]]; }
	int	get(const double x, const double y, const double z) const;
	int	get(const Point& r) const { return get(r.x, r.y, r.z); }

	int	operator () (const int i, const int j, const int k) const
			{ return get(i,j,k); }
	int	operator () (const double x, const double y, const double z) const
			{ return get(x,y,z); }
	int	operator () (const Point& r) const
			{ return get(r.x, r.y, r.z); }

	word	kreg(const int i) const	{ return _kreg[i]; }
	word	addKreg(word region);
	word	maxRegion() const;

	word	data(const int ptr) const
			{ return _data[ptr]; }
	word	data(const int i, const int j, const int k) const
			{ return _data[k*nynx + j*nx + i]; }
	word	data(const double x, const double y, const double z) const {
				int i,j,k;
				if (voxelijk(x,y,z,&i,&j,&k))
					return data(i,j,k);
				else
					return 0xFFFF;
			}
	void	data(const int i, const int j, const int k, word value)
			{ _data[k*nynx + j*nx + i] = value; }

	double	voxelx(int i)	const	{ return xlow + i*dx; }
	double	voxely(int j)	const	{ return ylow + j*dy; }
	double	voxelz(int k)	const	{ return zlow + k*dz; }

	double	voxelcx(int i)	const	{ return xlow + i*dx + dx/2.0; }
	double	voxelcy(int j)	const	{ return ylow + j*dy + dy/2.0; }
	double	voxelcz(int k)	const	{ return zlow + k*dz + dz/2.0; }

	// WARNING no transformation is applied!!!
	int	voxeli(const double x) const	{ return Int((x - xlow)/dx); }
	int	voxelj(const double y) const	{ return Int((y - ylow)/dy); }
	int	voxelk(const double z) const	{ return Int((z - zlow)/dz); }

	bool	voxelijk(double x, double y, double z, int *i, int *j, int *k) const;

	bool	intersectRay(const Point& p, const Vector &d, double *tmin) const;
	Vector	normal(const Point& r, const Vector& d) const;

	void	clearMatrix()			{ _hasMatrix = false; }
	bool	hasMatrix()		const	{ return _hasMatrix; }
	void	matrix(const Matrix4& M)	{ _matrix.copy(M); _hasMatrix=true;}
const	Matrix4& matrix()		const	{ return _matrix; }

	/* roi structures */
	int		nroi()		const	{ return _nroi; }
	word		roi(const int ptr)	{ return _roi[ptr]; }
	word		roi(int i, int j, int k){ return _roi[k*nynx + j*nx + i]; }
	std::string	roiName(int i)	const	{ return _roiName[i]; }
	dword		roiColor(int i)	const	{ return _roiColor[i]; }
	ROICombination&	roiComb(const int ptr)	{ return _roiComb[_roi[ptr]]; }
	ROICombination&	roiComb(int i, int j, int k) { return _roiComb[_roi[k*nynx + j*nx + i]]; }

	size_t	size()		const	{ return nynx*nz; }
	size_t	memory()	const	{
			size_t mem = 0;
			for (int i=0; i<_nroi; i++)
				mem += _roiName[i].length() + sizeof(std::string);
			for (int i=0; i<_roiCombN; i++)
				mem += _roiComb[i].memory();
			return sizeof(GVoxel)
				+ size()*sizeof(word)
				+ (mo+1)*sizeof(word)
				+ _nroi*sizeof(dword)	// roiColor
				+ 
				+ mem;
				// FIXME roi structure size
		}

	void	prefetch()	const	{
			__builtin_prefetch(_kreg);
			__builtin_prefetch(_data);
		}

friend class VVoxel;
}; // GVoxel

/* ================================== VVoxel ================================ */
/** Viewer voxel structure */
class VVoxel {
private:
	GVoxel&	 _voxel;	// pointer to voxel file
	dword	*_color;	// Color assigned to material color[_kreg[data[]]]
	int	 _no;		// number of colors allocated
	int	 _roiAlpha;	// alpha blend color with color of roi
	bool	*_roiShow;	// roi structures which are visible

private:
	int	 _lastRoi;	// last roi hit

public:
	VVoxel(GVoxel& v) : _voxel(v), _color(NULL), _no(0), _roiAlpha(128), _roiShow(NULL), _lastRoi(0) {}
	~VVoxel()	{ cleanup(); }
	void	cleanup();
	void	allocate();

	GVoxel&	voxel()	const	{ return _voxel; }

	int	no() const	{ return _voxel.no; }
	int	mo() const	{ return _voxel.mo; }

	int	nx() const	{ return _voxel.nx; }
	int	ny() const	{ return _voxel.ny; }
	int	nz() const	{ return _voxel.nz; }

	double	dx() const	{ return _voxel.dx; }
	double	dy() const	{ return _voxel.dy; }
	double	dz() const	{ return _voxel.dz; }

	double	xlow() const	{ return _voxel.xlow; }
	double	ylow() const	{ return _voxel.ylow; }
	double	zlow() const	{ return _voxel.zlow; }

	double	xhigh() const	{ return _voxel.xhigh; }
	double	yhigh() const	{ return _voxel.yhigh; }
	double	zhigh() const	{ return _voxel.zhigh; }

	double	voxelx(int i)	const	{ return _voxel.voxelx(i); }
	double	voxely(int j)	const	{ return _voxel.voxely(j); }
	double	voxelz(int k)	const	{ return _voxel.voxelz(k); }

	int	voxeli(double x) const	{ return _voxel.voxeli(x); }
	int	voxelj(double y) const	{ return _voxel.voxelj(y); }
	int	voxelk(double z) const	{ return _voxel.voxelk(z); }

	int	get(const int i, const int j, const int k) const
			{ return _voxel.get(i,j,k); }
	int	get(double x, double y, double z) const
			{ return _voxel.get(x,y,z); }
	int	get(const Point& r) const { return get(r.x, r.y, r.z); }

	int	operator () (const int i, const int j, const int k) const { return get(i,j,k); }
	int	operator () (const double x, const double y, const double z) const { return get(x,y,z); }
	int	operator () (const Point& r) const { return get(r.x, r.y, r.z); }

	void	color(int reg, dword c);
	dword	color(int reg)	const { assert(reg>=0 && reg<=_no); return _color[reg]; }
	dword	color(const int i, const int j, const int k) const {
			if (_color == NULL) return COLOR_TRANSPARENT;
			if (i<0 || j<0 || k<0 ||
			    i>=_voxel.nx || j>=_voxel.ny || k>=_voxel.nz)
				return COLOR_TRANSPARENT;
			else
				return _color[_voxel.get(i,j,k)];
		}
	dword	color(const double x, const double y, const double z, bool *ok) const;
	dword	color(const Point& r, bool *ok) const
			{ return color(r.x, r.y, r.z, ok); }
	dword	roiColor(const int reg);

	bool	roiShow(int i)	const	{ return _roiShow!=NULL?_roiShow[i] : false; }
	void	roiShow(int i, bool b);
	void	roiShowClear();

	int	roiAlpha()	const	{ return _roiAlpha; }
	void	roiAlpha(int a)		{ _roiAlpha = a; }

	dword	material(int ptr);

	bool	intersectRay(const Point& p, const Vector &d, double *tmin);

	bool	hasMatrix()	const	{ return _voxel._hasMatrix; }
const	Matrix4& matrix()	const	{ return _voxel._matrix; }

	size_t	memory()	const	{ return sizeof(VVoxel) + (_voxel.mo+1)*sizeof(dword); }
}; // VVoxel
#endif
