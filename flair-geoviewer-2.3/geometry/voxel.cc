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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>

#include <iostream>
#include <fstream>

#include "geo.h"
#include "voxel.h"
#include "color.h"
#include "fortran.h"

using namespace std;

/* =========================== ROICombination =============================== */
/** parse buffer for combinations */
void ROICombination::parse(void *buffer, int buflen)
{
	FortranParser parser(buffer, buflen);
	length = parser.readInt();
	roi    = new short[length];
	for (int i=0; i<length; i++)
		roi[i] = parser.readWord();
} // parse

/* =========================== ROICombination =============================== */
/** ROIPlanar */
ROIPlanar::~ROIPlanar()
{
	for (int i=0; i<nz; i++)
		if (slice[i]) delete slice[i];
	delete [] slice;
} // ROIPlanar

/** allocate structure slices */
void ROIPlanar::alloc(int anz)
{
	assert(nz==0);
	nz = anz;
	slice = new ROIPlanarSlice*[nz];
	bzero(slice, sizeof(ROIPlanarSlice*)*nz);
} // alloc

/** addZslice */
ROIPlanarSlice& ROIPlanar::addZslice(int z, int n)
{
	assert(slice[z]==NULL);
	slice[z] = new ROIPlanarSlice(n);
	return *slice[z];
} // addZslice

/* =============================== GVoxel =================================== */
/** GVoxel */
GVoxel::GVoxel()
	: _data(NULL), _kreg(NULL),
	  _roi(NULL), _roiName(NULL), _roiColor(NULL), _roiComb(NULL), _roiPlanar(NULL)
{
	nx = ny = nz = 0;
	no = mo = 0;
	_mtime  = 0;

	_nroi = _roiCombN = 0;
} // GVoxel

/** cleanup */
void GVoxel::cleanup()
{
	if (_data)	delete [] _data;
	if (_kreg)	delete [] _kreg;
	if (_roiName)	delete [] _roiName;
	if (_roiColor)	delete [] _roiColor;
	if (_roiComb)	delete [] _roiComb;
	if (_roiPlanar)	delete [] _roiPlanar;
	if (_roi)	delete [] _roi;
	filename.clear();
	_data      = NULL;
	_kreg      = NULL;
	_roiName   = NULL;
	_roiColor  = NULL;
	_roiComb   = NULL;
	_roiPlanar = NULL;
	_roi       = NULL;
	nx = ny = nz = 0;
	no = mo = 0;
	_mtime  = 0;
	_hasMatrix = false;
	_nroi = _roiCombN = 0;
} // cleanup

/** @return maximum region in kreg */
word GVoxel::maxRegion() const
{
	word maxreg = 0;
	for (int i=0; i<=mo; i++)
		maxreg = Max(maxreg, _kreg[i]);
	return maxreg;
} // maxRegion

/** addKreg */
word GVoxel::addKreg(word region)
{
	assert(_kreg);
	word maxreg = maxRegion();
	if (region > maxreg)
		region = maxreg + 1;

	// Increase the kreg array
	mo++;
	no++;
	word* kregnew = new word[mo+1];
	memcpy(kregnew, _kreg, sizeof(word)*mo);
	kregnew[mo] = region;

	delete [] _kreg;
	_kreg = kregnew;

	return region;
} // addKreg

/** load
 * @return	true on success
 *		false on failure
 */
bool GVoxel::load(const char *fn, bool loadCards)
{
	char	buffer[5120];
	int	length, readin, buflen;
	struct stat st;
	FortranParser parser;

	// Check if information is valid
	if (stat(fn, &st)) {
		cleanup();
		return false;
	}
	if (!filename.compare(fn) && st.st_mtime==_mtime) return true;

	cleanup();

	FortranFile f(fn, "rb");
	if (!f) return false;

	// Remember modification name and time
	filename = fn;
	_mtime = st.st_mtime;

	/* read title */
	memset(title, 0, sizeof(title));
	if (f.read(title,sizeof(title))<0) {
		cerr << "ERROR reading voxel title" << endl;
		goto ERROR;
	}
	FortranParser::strip(title);

	/* read arrays header  */
//	if ((buflen=f.read(buffer,sizeof(buffer)))<0 || (buflen!=20 && buflen!=32)) {
	if ((buflen=f.read(buffer,sizeof(buffer)))<0 || (buflen!=20 && buflen!=36)) {
		cerr << "ERROR reading voxel header" << endl;
		goto ERROR;
	}

	parser(buffer, buflen);
	nx = parser.readInt();
	ny = parser.readInt();
	nz = parser.readInt();
	no = parser.readInt();
	mo = parser.readInt();
	if (parser.more()) {
		_nroi       = parser.readInt();
		_roiCombN   = parser.readInt();
		_roiMaxComb = parser.readInt();
		_roiPlanarN = parser.readInt();
	}

	/* read voxel dimensions */
	if ((buflen=f.read(buffer,sizeof(buffer))) != 24) {
		cerr << "ERROR reading voxel dimensions" << endl;
		goto ERROR;
	}

	parser(buffer);
	dx = parser.readDouble();
	dy = parser.readDouble();
	dz = parser.readDouble();
	/* allocate arrays */
#if _DEBUG>1
	printf("VOXEL\tTitle: %s\n",title);
	printf("\tnx=%d ny=%d nz=%d\n",nx,ny,nz);
	printf("\tdx=%g dy=%g dz=%g\n",dx,dy,dz);
	printf("\tDx=%g Dy=%g Dz=%g\n",(double)nx*dx,(double)ny*dy,(double)nz*dz);
	printf("\tno=%d mo=%d\n",no,mo);
#endif
	nynx = ny*nx;

	/* Read data */
	_data  = new word[size()];
	length = size()*sizeof(word);
	readin = f.read(_data, length);
	if (readin != length) {
		fprintf(stderr,"ERROR wrong data length %d expected %d\n",readin,length);
		goto ERROR;
	}

	/* Read region definitions */
	_kreg = new word[mo+1];
	_kreg[0] = 0;
	readin = f.read(_kreg+1, mo*sizeof(word));
	if (readin != mo*(int)sizeof(word)) {
		fprintf(stderr,"ERROR wrong kreg length %d expected %d\n",readin,mo*(int)sizeof(word));
		goto ERROR;
	}

	/* ROI structures */
	if (_nroi) {
		/* Read the names */
		_roiName  = new std::string[_nroi+1];
		_roiColor = new dword[_nroi+1];
		for (int i=1; i<=_nroi; i++) {
			char name[64];
			//memset(buffer, 0, sizeof(buffer));
			if (f.read(buffer,sizeof(buffer)) != 68) {
				cerr << "ERROR reading voxel ROI names" << endl;
				goto ERROR;
			}
			parser(buffer);
			_roiColor[i] = parser.readInt();
			parser.read(name,sizeof(name));
			_roiName[i] = name;
			DUMP(if (buffer[0]) cout << "ROI: " << i+1 << ": " << _roiName[i] << ' ' << _roiColor[i] << endl;)
		}

		/* Read the combination dictionary */
		_roiComb = new ROICombination[_roiCombN+1];
		for (int i=1; i<=_roiCombN; i++) {
			if ((buflen=f.read(buffer,sizeof(buffer)))<0) {
				cerr << "ERROR reading voxel ROI combination dictionary" << endl;
				goto ERROR;
			}
			// roiComb[0] is the basic combination with nothing ink:W
			_roiComb[i].parse(buffer, buflen);
			DUMP(printf("ROICMB %d ",i); _roiComb[i].print();)
		}

		/* Read the roi struct data */
		_roi = new word[size()];
		length = size()*sizeof(word);
		readin = f.read(_roi, length);
		if (readin != length) {
			fprintf(stderr,"ERROR wrong data length %d expected %d\n",readin,length);
			goto ERROR;
		}

//		memset(stat,0,sizeof(stat));
//		for (int i=0; i<length; i++)
//			stat[_roi[i]]++;
//		for (int i=0; i<_roiCombN; i++)
//			cout << i << " : " << stat[i] << endl;

		/* Read the planar information */
		if (_roiPlanarN)
			_roiPlanar = new ROIPlanar[_nroi+1];
		for (int i=0; i<_roiPlanarN; i++) {
			if ((buflen=f.read(buffer,sizeof(buffer))) != 8) {
				cerr << "ERROR reading planar header" << endl;
				goto ERROR;
			}
			parser(buffer, buflen);
			/*int r =*/ parser.readInt();
			int n = parser.readInt();
#if _DEBUG>2
			cout << "ROI:" << r;
			cout << " r=" << r << endl;
			cout << " n=" << n << endl;
#endif

//			_roiPlanar[r].alloc(nz);

			for (int j=0; j<n; j++) {
//				buflen=f.read(buffer,sizeof(buffer));
				int buflen2;
				//				byte* buffer2 = f.readBuffer(&buflen2);
				uint8_t* buffer2 = f.readBuffer(&buflen2);				if (buffer2==NULL) {
					cerr << "ERROR reading planar information" << endl;
					goto ERROR;
				}
				parser(buffer2, buflen2);
				/*int zi =*/ parser.readInt();
				/*int m  =*/ parser.readInt();

				//assert(buflen == m*2*sizeof(float) + 2*sizeof(int));
#if _DEBUG>2
				cout << '[' << i << "] z=" << zi;
				cout << " m=" << m;
				for (int k=0; k<m; k++) {
					float x = parser.readFloat();
					float y = parser.readFloat();
					cout << " (" << x << ',' << y << ')';
				}
				cout << endl;
#endif
				delete [] buffer2;
			}
		}
	}

	/* XXX maybe read internal cards */
	if (loadCards) {
		char card[81];
		while (f) {
			f.read(card, sizeof(card));
			card[80] = 0;
			cards.add(card);
			DUMP(cout << "Voxel: " << card << endl;)
		}
	}
	return true;

ERROR:
	fprintf(stderr,"ERROR reading voxel file \"%s\"\n",filename.c_str());
	if (_data) delete [] _data;
	_data = NULL;
	return false;
} // load

/** save
 * @return	true on success
 *		false on failure
 */
bool GVoxel::save(const char* fn)
{
	char	buffer[512];
	if (fn) filename = fn;

	FortranFile f(filename.c_str(), "wb");
	if (!f) return false;

	FortranParser parser(buffer, sizeof(buffer));

	// Title
	parser.write(title, 80);
	f.write(parser);

	// Dimensions
	parser.reset();
	parser.write(nx);
	parser.write(ny);
	parser.write(nz);
	parser.write(no);
	parser.write(mo);
//	if (_nroi) {
//		parser.write(_nroi);
//		parser.write(_roiCombN);
//		parser.write(_roiMaxComb);
//	}
	f.write(parser);

	// Size
	parser.reset();
	parser.write(dx);
	parser.write(dy);
	parser.write(dz);
	f.write(parser);

	// write data
	f.write(_data, size()*sizeof(word));

	// write _kreg (skip first)
	f.write(_kreg+1, mo*sizeof(word));

	for (int i=0; i<cards.count(); i++) {
		parser.reset();
		parser.write(cards[i].c_str(), 80);
		f.write(parser);
	}

	f.close();
	return true;
} // save

/** calculate voxel limits */
void GVoxel::calcLimits()
{
	xhigh = voxelx(nx);
	yhigh = voxely(ny);
	zhigh = voxelz(nz);
} // calcLimits

/** voxelijk */
bool GVoxel::voxelijk(double x, double y, double z, int *i, int *j, int *k) const
{
	if (hasMatrix()) matrix().transform(&x,&y,&z);
	// be more flexible on the limits
	double chk = dx*SMALL;

	if (Abs(x-xlow) < chk)
		*i = 0;
	else
	if (Abs(x-xhigh) < chk)
		*i = nx-1;
	else {
		*i = voxeli(x);
		if (*i<0 || *i>=nx)
			return false;
	}

	chk = dy*SMALL;
	if (Abs(y-ylow) < chk)
		*j = 0;
	else
	if (Abs(y-yhigh) < chk)
		*j = ny-1;
	else {
		*j = voxelj(y);
		if (*j<0 || *j>=ny)
			return false;
	}

	chk = dz*SMALL;
	if (Abs(z-zlow) < chk)
		*k = 0;
	else
	if (Abs(z-zhigh) < chk)
		*k = nz-1;
	else {
		*k = voxelk(z);
		if (*k<0 || *k>=nz)
			return false;
	}
	return true;
} // voxelijk

/** get _kreg of data of location
 * @param x,y,z	location to return value
 * @return _kreg[data[x,y,z]] at location or -1 when error
 */
int GVoxel::get(double x, double y, double z) const
{
	int i, j, k;

	if (hasMatrix()) matrix().transform(&x,&y,&z);

	// be more flexible on the limits
	double chk = dx*SMALL;

	if (Abs(x-xlow) < chk)
		i = 0;
	else
	if (Abs(x-xhigh) < chk)
		i = nx-1;
	else {
		i = voxeli(x);
		if (i<0 || i>=nx)
			return -1;
	}

	chk = dy*SMALL;
	if (Abs(y-ylow) < chk)
		j = 0;
	else
	if (Abs(y-yhigh) < chk)
		j = ny-1;
	else {
		j = voxelj(y);
		if (j<0 || j>=ny)
			return -1;
	}

	chk = dz*SMALL;
	if (Abs(z-zlow) < chk)
		k = 0;
	else
	if (Abs(z-zhigh) < chk)
		k = nz-1;
	else {
		k = voxelk(z);
		if (k<0 || k>=nz)
			return -1;
	}

	return get(i,j,k);
} // get

/** intersectRay find the distance to the next voxel with different
 * material
 * @param p	initial position
 * @param d	direction
 * @param tmin	minimum distance of intersection
 * @return true if an intersection is found
 */
bool GVoxel::intersectRay(const Point& /*p*/, const Vector& /*d*/, double* /*tmin*/) const
{
	// FIXME can be cached for faster searching
	// position, direction and step
	assert(0);	// should use the VVoxel one for the moment
	return false;
} // intersectRay

/** return normal vector
 * @param r	location of normal
 * @return	returned normal vector
 */
Vector GVoxel::normal(const Point& p, const Vector& d) const
{
	/* find the closest boundary */
	double remX, remY, remZ;

	Point pos;
	Vector dir;
	if (hasMatrix()) {
		// Convert to voxel coordinates
		pos = matrix() * p;
		dir = matrix() * d;
//		dir = matrix().multVector(d);
	} else {
		pos = p;
		dir = d;
	}

	/* First on X */
	double di = (pos.x-xlow)/dx;
	int i = Round(di);
	if (i>=0 && i<=nx)
		remX = Abs((double)i - di);
	else
		remX = INFINITE;

	/* Then on Y */
	di = (pos.y-ylow)/dy;
	i = Round(di);
	if (i>=0 && i<=ny)
		remY = Abs((double)i - di);
	else
		remY = INFINITE;

	/* Finally on Z */
	di = (pos.z-zlow)/dz;
	i = Round(di);
	if (i>=0 && i<=nz)
		remZ = Abs((double)i - di);
	else
		remZ = INFINITE;

	if (remX<remY && remX<remZ) {
		if (dir.x>0.0)
			return -Vector::Xo;
		else
			return Vector::Xo;
	}

	if (remY<remX && remY<remZ) {
		if (dir.y>0.0)
			return -Vector::Yo;
		else
			return Vector::Yo;
	}

	if (dir.z>0.0)
		return -Vector::Zo;
	else
		return Vector::Zo;
} // normal

///////////////////////////////// VVoxel ////////////////////////////////////
/** cleanup */
void VVoxel::cleanup()
{
	roiShowClear();
	if (_color) delete [] _color;
	_color    = NULL;
	_no       = 0;
} // cleanup

/** allocate color array if needed */
void VVoxel::allocate()
{
	if (_color==NULL || _no != no()) {
		cleanup();
		_no = no();
		_color = new dword[_no+1];
		memset(_color, 0xFF, (_no+1)*sizeof(dword));
	}
} // allocate

/** get color of location
 * @param x,y,z	location to return value
 * @return value at (x,y,z) location
 */
dword VVoxel::color(const double x, const double y, const double z, bool *ok) const
{
	if (_color == NULL) return COLOR_TRANSPARENT;
	int d = _voxel.get(x,y,z);
	if (d<0) {
		*ok = false;
		return COLOR_TRANSPARENT;
	} else {
		*ok = true;
		return _color[d];
	}
} // color

/** set color ith (Base 0) */
void VVoxel::color(int reg, dword c)
{
	if (!_color) allocate();
	if (reg<0 || reg>_no) return;
	_color[reg] = c;
} // color

/** roiColor */
dword VVoxel::roiColor(const int reg)
{
	if (_lastRoi==0 || _roiShow==NULL || _roiAlpha==0)
		return color(reg);

	Color32 cmat, croi;
	cmat.val = color(reg);
	cmat.rgb.alpha = 0xff;
	croi.val = _voxel._roiColor[_lastRoi];
	croi.rgb.alpha = _roiAlpha;
	return alphaBlend(croi, cmat);
} // roiColor

/** roiShowClear */
void VVoxel::roiShowClear()
{
	if (_roiShow)	delete [] _roiShow;
	_roiShow = NULL;
	_lastRoi = 0;
} // roiShowClear

/** roiShow */
void VVoxel::roiShow(int idx, bool b)
{
	if (idx<0 || idx>_voxel._nroi) return;
	if (_roiShow == NULL) {
		_roiShow = new bool[_voxel._nroi+1];
		memset(_roiShow, 1, sizeof(bool)*(_voxel._nroi+1));
	}
	_roiShow[idx] = b;
} // roiShow

/** return material at position ptr.
 * For the moment return the color at position ptr
 */
inline dword VVoxel::material(int ptr)
{
	dword mat = color(_voxel.get(ptr));
	if (_roiShow==NULL || mat == COLOR_TRANSPARENT) return mat;

	int comb = _voxel._roi[ptr];
	if (comb==0) return COLOR_TRANSPARENT; //mat;
	//_lastRoi = 0;
	// FIXME maybe I should keep the original structure....
	for (int i=0; i<_voxel._roiComb[comb].length; i++) {
		int roi = Abs(_voxel._roiComb[comb][i]);
		//int roi = -_voxel._roiComb[comb][i];
		if (roi>0 && !_roiShow[roi]) {
			_lastRoi = roi;
			return mat;
		}
	}
	return COLOR_TRANSPARENT;
} // material

/** intersectRay find the distance to the next voxel with different
 * material
 * @param p	initial position
 * @param d	direction
 * @param tmin	minimum distance of intersection
 * @return true if an intersection is found
 */
bool VVoxel::intersectRay(const Point& p, const Vector &d, double *tmin)
{
	if (_color == NULL) return false;
	//_voxel.prefetch();
	//__builtin_prefetch(_voxel._kreg);
	//__builtin_prefetch(_voxel._data);

	// FIXME can be cached for faster searching
	// position, direction and step
	Point pos = p + d*(*tmin);
	Vector dir;
	if (hasMatrix()) {
		// Convert to voxel coordinates
		pos = matrix() * pos;
		dir = matrix() * d;
//		dir = matrix().multVector(d);
	} else
		dir = d;

	// Find location
	int i = voxeli(pos.x);
	if (i<0) {
		if (dir.x<=0.0) return false;
		else i = 0;	// correct for possible rounding error
	} else
	if (i>=nx()) {
		if (dir.x>=0.0) return false;
		else i = nx()-1;
	}

	int j = voxelj(pos.y);
	if (j<0) {
		if (dir.y<=0.0) return false;
		else j = 0;	// correct for possible rounding error
	} else
	if (j>=ny()) {
		if (dir.y>=0.0) return false;
		else j = ny()-1;
	}

	int k = voxelk(pos.z);
	if (k<0) {
		if (dir.z<=0.0) return false;
		else k = 0;	// correct for possible rounding error
	} else
	if (k>=nz()) {
		if (dir.z>=0.0) return false;
		else k = nz()-1;
	}

	int ptr = k*_voxel.nynx + j*_voxel.nx + i;
	dword mat = material(ptr);
	if (mat != COLOR_TRANSPARENT) return true;

	double tx,    ty,    tz;		// distance to next intersection with x/y/z plane
	double tdx,   tdy,   tdz;		// distance increase for every x/y/z plane
	int    stepx, stepy, stepz;		// step direction -1/0/+1
	#define ptrx stepx
	int           ptry,  ptrz;		// ptr change ptrx=stepx
	int    outx,  outy,  outz;		// ending condition

	if (dir.x > SMALL) {
		tx    = (voxelx(i+1) - pos.x) / dir.x;
		tdx   = dx() / dir.x;
		stepx = 1;
		outx  = nx();
	} else
	if (dir.x < -SMALL) {
		tx    = (voxelx(i) - pos.x) / dir.x;
		tdx   = -dx()/dir.x;
		stepx = -1;
		outx  = -1;
	} else {
		tx    = INFINITE;
		tdx   = 0.0;
		stepx = 0;
		outx  = 0;
	}

	if (dir.y > SMALL) {
		ty    = (voxely(j+1) - pos.y) / dir.y;
		tdy   = dy() / dir.y;
		stepy = 1;
		ptry  = nx();
		outy  = ny();
	} else
	if (dir.y < -SMALL) {
		ty    = (voxely(j) - pos.y) / dir.y;
		tdy   = -dy()/dir.y;
		stepy = -1;
		ptry  = -nx();
		outy  = -1;
	} else {
		ty    = INFINITE;
		tdy   = 0.0;
		stepy = 0;
		ptry  = 0;
		outy  = 0;
	}

	if (dir.z > SMALL) {
		tz    = (voxelz(k+1) - pos.z) / dir.z;
		tdz   = dz() / dir.z;
		stepz = 1;
		ptrz  = _voxel.nynx;
		outz  = nz();
	} else
	if (dir.z < -SMALL) {
		tz    = (voxelz(k) - pos.z) / dir.z;
		tdz   = -dz()/dir.z;
		stepz = -1;
		ptrz  = -_voxel.nynx;
		outz  = -1;
	} else {
		tz    = INFINITE;
		tdz   = 0.0;
		stepz = 0;
		ptrz  = 0;
		outz  = 0;
	}

	// Loop forever until we find a different material
	while (true) {
		if (tx<ty && tx<tz) {
			//__builtin_prefetch(_voxel._data+ptr+ptrx);
			i += stepx;
			if (i==outx) return false;
			ptr += ptrx;
			if (material(ptr) != mat) {
				*tmin += tx;
				return true;
			}
			tx += tdx;
		} else
		if (ty<tz) {
			//__builtin_prefetch(_voxel._data+ptr+ptry);
			j += stepy;
			if (j==outy) return false;
			ptr += ptry;
			if (material(ptr) != mat) {
				*tmin += ty;
				return true;
			}
			ty += tdy;
		} else {
			//__builtin_prefetch(_voxel._data+ptr+ptrz);
			k += stepz;
			if (k==outz) return false;
			ptr += ptrz;
			if (material(ptr) != mat) {
				*tmin += tz;
				return true;
			}
			tz += tdz;
		}
	}
} // intersectRay
