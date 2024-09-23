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

#include "os.h"
#include "eps.h"
#include "bmath.h"
#include "usrbin.h"
#include "fortran.h"

using namespace std;

#define TOLERANCE 1e-7

char* Usrbin::_typeStr[] = {
		"XYZ_point",
		"RPhiZ_point",
		"region_point",
		"|X|_point",
		"|Y|_point",
		"|Z|_point",
		"|XYZ|_point",
		"RPhi|Z|_point",
		"special_point",
		"",
		"XYZ",
		"RPhiZ",
		"region",
		"|X|",
		"|Y|",
		"|Z|",
		"|XYZ|",
		"RPhi|Z|",
		"special",
		""
	};

///////////////////////////////// Usrbin ////////////////////////////////////
/** cleanup */
void Usrbin::cleanup()
{
	filename.clear();
	if (data) delete [] data;
	data      = NULL;
	_detector = 0;
	nx = ny = nz = 0;
	x0 = y0  = 0.0;
	xofs = yofs = zofs = 0.0;
	runtitle[0] = 0;
	runtime[0]  = 0;
	titusb[0]   = 0;
	_norm       = 1.0;
	_logscale   = false;
	_checker    = false;
	_mtime      = 0;
	_hasMatrix  = false;

	lntzer = 0;
	bkusbn = 0.0;
	b2usbn = 0.0;
	tcusbn = 0.0;

	weipri = 0.0;
	ncase  = 0;
	mcase  = 0;
	nbatch = 0;
//	_total = 0.0;
} // cleanup

/** load
 * @return	true on success
 *		false on failure
 */
bool Usrbin::load(const char *fn, int det)
{
	char buffer[512];
	int  i, length;
	int  totSize = 0;
	struct stat st;

	// Check if information is valid
	if (stat(fn, &st)) {
		cleanup();
		return false;
	}
	if (!filename.compare(fn) &&
	     st.st_mtime==_mtime && det==_detector) return true;

	cleanup();

	FortranFile f(fn,"rb");
	if (!f) return false;

	filename = fn;

	// Remember modification time
	_mtime = st.st_mtime;

	length = f.read(buffer,sizeof(buffer)); //pbuffer = buffer;
	FortranParser parser(buffer, length);
	parser.read(runtitle, USRBIN_RUNTIT);
	parser.read(runtime,  USRBIN_RUNTIM);
	weipri = parser.readFloat();
	ncase  = parser.readInt();
	if (parser.more()) mcase  = parser.readInt();
	if (parser.more()) nbatch = parser.readInt();

	for (i=0; i<det; i++) {
		/* Read detector information
			WRITE (1) NB,TITUSB(NB),ITUSBN(NB),IDUSBN(NB),
				XLOW(NB),XHIGH(NB),NXBIN(NB),DXUSBN(NB),
				YLOW(NB),YHIGH(NB),NYBIN(NB),DYUSBN(NB),
				ZLOW(NB),ZHIGH(NB),NZBIN(NB),DZUSBN(NB),
				LNTZER(NB),BKUSBN(NB),B2USBN(NB),TCUSBN(NB)
		*/
		length = f.read(buffer, sizeof(buffer));
		if (length!=86) break;
		if (length<0) goto ERROR;
		parser(buffer);

		parser.skipInt();
		parser.read(titusb,USRBIN_TITUSB);
		_type  = (UsrbinType)parser.readInt();
		_score = parser.readInt();

		xlow  = (double)parser.readFloat();
		xhigh = (double)parser.readFloat();
		nx    =         parser.readInt();
		dx    = (double)parser.readFloat();

		ylow  = (double)parser.readFloat();
		yhigh = (double)parser.readFloat();
		ny    =         parser.readInt();
		dy    = (double)parser.readFloat();

		zlow  = (double)parser.readFloat();
		zhigh = (double)parser.readFloat();
		nz    =         parser.readInt();
		dz    = (double)parser.readFloat();

		if (_type!=Usrbin_region_point  && _type!=Usrbin_region &&
		    _type!=Usrbin_special_point && _type!=Usrbin_special) {
			// Re-normalize dx,dy,dz to get more precision from float
			dx  = (xhigh-xlow) / (double)nx;
			dy  = (yhigh-ylow) / (double)ny;
			dz  = (zhigh-zlow) / (double)nz;
		}

		lntzer = parser.readInt();
		bkusbn = parser.readFloat();
		b2usbn = parser.readFloat();
		tcusbn = parser.readFloat();
		/*
		printf("LNTZER %d\n", lntzer);
		printf("BKUSBN %g\n", bkusbn);
		printf("B2USBN %g\n", b2usbn);
		printf("TCUSBN %g\n", tcusbn);
		*/

		nynx = ny*nx;
		totSize = nynx*nz;

		/* skip detector */
		if (i<det-1) f.seek(8+totSize*sizeof(float), SEEK_CUR);
	}

	if (det==i && totSize>0) {
		data   = new float[totSize];
		length = totSize*sizeof(float);
		if (f.read(data, length) != length) goto ERROR;

		scanMinMax();

		printf("Usrbin: %d %s\n",det, titusb);
		printf("\tType:%d [%s]\tScore:%d\n",_type,_typeStr[_type],_score);
		printf("\tX: [%g .. %g] x %d (%g)\n",xlow, xhigh, nx, dx);
		printf("\tY: [%g .. %g] x %d (%g)\n",ylow, yhigh, ny, dy);
		printf("\tZ: [%g .. %g] x %d (%g)\n",zlow, zhigh, nz, dz);
		printf("\tMin=%-10.5g\tLogMin=%-10.5g\n",min,log10(min));
		printf("\tMax=%-10.5g\tLogMax=%-10.5g\n",max,log10(max));
	}
	_detector = det;
	return true;

ERROR:
	fprintf(stderr,"ERROR reading usrbin file \"%s\"\n",filename.c_str());
	cleanup();
	return true;
} // load

/** save
 * @return	true on success
 *		false on failure
 */
bool Usrbin::save(const char* fn)
{
	char	buffer[512];
	if (fn) filename = fn;

	FortranFile f(filename.c_str(), "wb");
	if (!f) return false;

	FortranParser parser(buffer, sizeof(buffer));

	// Write title
	parser.write(runtitle, USRBIN_RUNTIT);
	parser.write(runtime,  USRBIN_RUNTIM);
	parser.write(weipri);
	parser.write(ncase);
	parser.write(mcase);
	parser.write(nbatch);
	f.write(parser);

	// Write detector header
	parser.reset();
	parser.write(_detector);		// detector number
	parser.write(titusb, USRBIN_TITUSB);	// usrbin detector title
	parser.write((int)_type);		// type
	parser.write(_score);			// score

	parser.write((float)xlow);
	parser.write((float)xhigh);
	parser.write(nx);
	parser.write((float)dx);

	parser.write((float)ylow);
	parser.write((float)yhigh);
	parser.write(ny);
	parser.write((float)dy);

	parser.write((float)zlow);
	parser.write((float)zhigh);
	parser.write(nz);
	parser.write((float)dz);

	parser.write(lntzer);
	parser.write(bkusbn);
	parser.write(b2usbn);
	parser.write(tcusbn);
	assert(parser.length() == 86);
	f.write(parser);

	// Write data
	int len = nynx*nz*sizeof(float);
	if (f.write(data, len) != len)
		fprintf(stderr,"ERROR writing usrbin file \"%s\"\n",filename.c_str());

	return true;
} // save

/** create a new usrbin and initialize matrix to size ax*ay*az
 */
void Usrbin::create(int det, UsrbinType atype, int ascore,
			const double axlow, const double axhigh, const int anx,
			const double aylow, const double ayhigh, const int any,
			const double azlow, const double azhigh, const int anz)
{
	cleanup();

	_detector = det;
	_type  = atype;
	_score = ascore;

	nx    = Max(1,anx);
	ny    = Max(1,any);
	nz    = Max(1,anz);
	nynx  = nx*ny;
	data  = new float[nynx*nz];

	xlow  = axlow;
	ylow  = aylow;
	zlow  = azlow;

	xhigh = axhigh;
	yhigh = ayhigh;
	zhigh = azhigh;

	dx  = (xhigh-xlow) / (double)nx;
	dy  = (yhigh-ylow) / (double)ny;
	dz  = (zhigh-zlow) / (double)nz;

	assert(xlow < xhigh);
	assert(ylow < yhigh);
	assert(zlow < zhigh);
} // create

/** checker */
void Usrbin::checker(const UsrbinType atype,
			double axlow, double axhigh, const int anx,
			double aylow, double ayhigh, const int any,
			double azlow, double azhigh, const int anz)
{
	cleanup();

	_checker = true;

	if (axlow > axhigh) Swap(axlow, axhigh);

	_type = atype;
	xlow  = axlow;
	xhigh = axhigh;
	nx    = Max(1,anx);
	ny    = Max(1,any);

	dx = (xhigh - xlow)/(double)nx;

	if (_type==Usrbin_RPhiZ_point    || _type==Usrbin_RPhiZ ||
	    _type==Usrbin_RPhiZsym_point || _type==Usrbin_RPhiZsym) {
		// For R-Z-Phi use -PI to PI as limits
		ylow  = -PI;
		yhigh =  PI;
		x0    = aylow;
		y0    = ayhigh;
	} else {
		if (aylow > ayhigh) Swap(aylow, ayhigh);
		ylow  = aylow;
		yhigh = ayhigh;
		x0    = 0.0;
		y0    = 0.0;
	}
	dy = (yhigh - ylow)/(double)ny;

	if (azlow > azhigh) Swap(azlow, azhigh);
	zlow  = azlow;
	zhigh = azhigh;
	nz    = Max(1,anz);
	dz    = (zhigh - zlow)/(double)nz;

#if _DEBUG>1
	printf("Usrbin::setChecker it=%d\n",_type);
	printf("\tx:[%g..%g]x%d\n",xlow,xhigh,nx);
	printf("\ty:[%g..%g]x%d\n",ylow,yhigh,ny);
	printf("\tz:[%g..%g]x%d\n",zlow,zhigh,nz);
#endif

	if (data) {
		delete [] data;
		data = NULL;
	}
} // setChecker

/** offset */
void Usrbin::offset(const Point& ofs)
{
	xofs = ofs.x;
	yofs = ofs.y;
	zofs = ofs.z;
} // offset

/** scanMinMax */
void Usrbin::scanMinMax()
{
	/* scan for min and maximum */
	min =  1.0e30;
	max = -1.0e30;
	for (int j=0; j<nynx*nz; j++) {
		if (data[j]<min) min = data[j];
		if (data[j]>max) max = data[j];
	}
} // scanMinMax

/** norm */
void Usrbin::norm(const double n)
{
	if (_logscale)
		_norm = log10(n);
	else
		_norm = n;
} // norm

/** normalize usrbin by multiplying by n */
void Usrbin::normalize(double n)
{
	assert(!_checker && data);

	int totSize = nynx*nz;
	if (_logscale) {
		n = log10(n);
		for (int i=0; i<totSize; i++)
			data[i] += n;
	} else
		for (int i=0; i<totSize; i++)
			data[i] *= n;
} // normalize

/** convert */
bool Usrbin::convert(bool tolog)
{
	if (tolog == _logscale) return false;
	if (data==NULL) return true;
	bool err = false;
	int totSize = nynx*nz;
	for (int i=0; i<totSize; i++)
		if (tolog) {
			if (data[i]>0.0)
				data[i] = log10(data[i]);
			else {
				data[i] = -1E+30;
				err = true;
			}
		} else
			data[i] = pow(10.0,data[i]);
	if (tolog)
		_norm = log10(_norm);
	else
		_norm = pow(10.0,_norm);
	_logscale = tolog;
	return err;
} // convert

/** convert x,y,z coordinates to i,j,k in usrbin
 * @param x,y,z	location to return value
 * @return true if inside USRBIN, false otherwise
 */
inline bool Usrbin::xyz2ijk(double x, double y, double z, int* i, int* j, int* k) const
{
	double	r, phi;

	if (!_checker && !data) return false;

	if (hasMatrix()) matrix().transform(&x,&y,&z);
	x -= xofs;
	y -= yofs;
	z -= zofs;

	switch (_type) {
		case Usrbin_XYZ_point:		// Cartesian no symmetry
		case Usrbin_XYZ:
			break;

		case Usrbin_RPhiZ_point:	// R-Phi-Z
		case Usrbin_RPhiZ:
			x  -= x0;
			y  -= y0;
			r   = hypot(x,y);
			phi = atan2(y,x);
			x   = r;
			y   = phi;
			break;

		case Usrbin_Xsym_point:		// X-symmetry
		case Usrbin_Xsym:
			if (x<0.0) x = -x;
			break;

		case Usrbin_Ysym_point:		// Y-symmetry
		case Usrbin_Ysym:
			if (y<0.0) y = -y;
			break;

		case Usrbin_Zsym_point:		// Z-symmetry
		case Usrbin_Zsym:
			if (z<0.0) z = -z;
			break;

		case Usrbin_XYZsym_point:	// X,Y,Z-symmetry
			if (x<0.0) x = -x;
			if (y<0.0) y = -y;
			if (z<0.0) z = -z;
			break;

		case Usrbin_RPhiZsym_point:	// R-Phi-Z, Z-symmetry
		case Usrbin_RPhiZsym:
			x  -= x0;
			y  -= y0;
			r   = hypot(x,y);
			phi = atan2(y,x);
			if (z<0.0) z = -z;
			x   = r;
			y   = phi;
			break;

		default:	// special (ignored)
			return 0.0;
	}

	// be more flexible on the limits
	if (Abs(x-xlow)<Abs(x)*TOLERANCE)
		*i = 0;
	else
	if (Abs(x-xhigh)<Abs(x)*TOLERANCE)
		*i = nx-1;
	else {
		*i = usrbini(x);
		if (*i<0 || *i>=nx)
			return false;
	}

	if (Abs(y-ylow)<Abs(y)*TOLERANCE)
		*j = 0;
	else
	if (Abs(y-yhigh)<Abs(y)*TOLERANCE)
		*j = ny-1;
	else {
		*j = usrbinj(y);
		if (*j<0 || *j>=ny)
			return false;
	}

	if (Abs(z-zlow)<Abs(z)*TOLERANCE)
		*k = 0;
	else
	if (Abs(z-zhigh)<Abs(z)*TOLERANCE)
		*k = nz-1;
	else {
		*k = usrbink(z);
		if (*k<0 || *k>=nz)
			return false;
	}
	return true;
} // xyz2ijk

/** get stored value (log10 in case of log)
 * @param x,y,z	location to return value
 * @return value at (x,y,z) location
 */
double Usrbin::getData(const double x, const double y, const double z, bool *ok) const
{
	int i,j,k;

	if (xyz2ijk(x,y,z, &i, &j, &k)) {
		*ok = true;
		if (_checker)
			return ((i+j+k)&1? 0.0 : 1.0);
		else
			return _normValue((double)data[k*nynx + j*nx + i]);
	} else {
		*ok = false;
		return 0.0;
	}
} // getData

/** get stored value (log10 in case of log) by region (ONLY for binning per region)
 * @param region
 * @return value for region
 */
double Usrbin::getData(const int region, bool *ok) const
{
	*ok = false;

	if (!data) return 0.0;

	if (_type!=Usrbin_region_point && _type!=Usrbin_region) return 0.0;

	// First set...
	int rstart = (int)xlow  - 1;
	int rend   = (int)xhigh - 1;
	int rstep  = (int)dx;
	//int rlen   = (int)nx;

	if (rstart < region && region > rend) return 0.0;

	if ((region - rstart) % rstep) return 0.0;

	*ok = true;
	return _normValue((double)data[(region - rstart) / rstep]);
	// FIXME should check second and third set
} // getData

/** get stored value (log10 in case of log)
 * @param r,phi,z	location to return value
 * @return value at (r,phi,z) location
 */
double Usrbin::getDataRPhiZ(double r, double phi, double z, bool *ok) const
{
	*ok = false;

	assert( _type==Usrbin_RPhiZ_point    || _type==Usrbin_RPhiZ ||
		_type==Usrbin_RPhiZsym_point || _type==Usrbin_RPhiZsym);
	assert(!_checker && data);

	if (_type==Usrbin_RPhiZsym_point || _type==Usrbin_RPhiZsym)
		if (z<0.0) z = -z;

	int i, j, k;
	// be more flexible on the limits
	if (Abs(r-xlow)<Abs(r)*TOLERANCE)
		i = 0;
	else
	if (Abs(r-xhigh)<Abs(r)*TOLERANCE)
		i = nx-1;
	else {
		i = Int((r-xlow)/dx);
		if (i<0 || i>=nx)
			return 0.0;
	}

	if (Abs(phi-ylow)<Abs(phi)*TOLERANCE)
		j = 0;
	else
	if (Abs(phi-yhigh)<Abs(phi)*TOLERANCE)
		j = ny-1;
	else {
		j = Int((phi-ylow)/dy);
		if (j<0 || j>=ny)
			return 0.0;
	}

	if (Abs(z-zlow)<Abs(z)*TOLERANCE)
		k = 0;
	else
	if (Abs(z-zhigh)<Abs(z)*TOLERANCE)
		k = nz-1;
	else {
		k = Int((z-zlow)/dz);
		if (k<0 || k>=nz)
			return 0.0;
	}

	*ok = true;
	return _normValue(data[k*nynx + j*nx + i]);
} // getDataRPhiZ

/** distribute amount 'value' to all the usrbin cells that the line start-end will cross.
 * The amount scored per cell is proportional to the distance travelled in the cell.
 * @param start, stop   starting and ending points of the line to distribute the amount value.
 *                      warning points can lie also outside the usrbin
 * @param value		amount to score (times the weight of the particle)
 * @param tracklength	if true score the amount as tracklength (do not divide by the total length)
 * @return the amount that was scored
 */
double Usrbin::add(Point start, Point stop, double value, bool tracklength)
{
	assert(!_logscale && !_checker && data);

	if (hasMatrix()) {
		// Convert to usrbin coordinates
		matrix().transform(start);
		matrix().transform(stop);
	}
	start.x -= xofs;	// FIXME convert offset to Point
	start.y -= yofs;
	start.z -= zofs;

	stop.x -= xofs;
	stop.y -= yofs;
	stop.z -= zofs;

	Vector dir = stop - start;
	double totallength = dir.normalize();
	if (!tracklength) value /= totallength;		// already divide by length

	//cout << "total length=" << totallength << endl;

	// for the moment it works only with XYZ
	assert(_type==Usrbin_XYZ_point || _type==Usrbin_XYZ);

	double t = 0.0;					// distance travelled to far

	// move to boundary
	if (start.x < xmin()) {
		if (dir.x<=0.0) return 0.0;		// going out ignore
		double f = (xmin() - start.x) / dir.x;	// distance to travel
		t += f;					// increase distance traveled
		start.x  = xmin();			// new position
		start.y += f*dir.y;
		start.z += f*dir.z;
	}
	if (start.y < ymin()) {
		if (dir.y<=0.0) return 0.0;
		double f = (ymin() - start.y) / dir.y;
		t += f;
		start.x += f*dir.x;
		start.y  = ymin();
		start.z += f*dir.z;
	}
	if (start.z < zmin()) {
		if (dir.z<=0.0) return 0.0;
		double f = (zmin() - start.z) / dir.z;
		t += f;
		start.x += f*dir.x;
		start.y += f*dir.y;
		start.z  = zmin();
	}

	if (start.x > xmax()) {
		if (dir.x>=0.0) return 0.0;
		double f = (xmax() - start.x) / dir.x;
		t += f;
		start.x  = xmax();
		start.y += f*dir.y;
		start.z += f*dir.z;
	}
	if (start.y > ymax()) {
		if (dir.y>=0.0) return 0.0;
		double f = (ymax() - start.y) / dir.y;
		t += f;
		start.x += f*dir.x;
		start.y  = ymax();
		start.z += f*dir.z;
	}
	if (start.z > zmax()) {
		if (dir.z>=0.0) return 0.0;
		double f = (zmax() - start.z) / dir.z;
		t += f;
		start.x += f*dir.x;
		start.y += f*dir.y;
		start.z  = zmax();
	}
	//cout << "starting length=" << t << endl;

	// Find starting location
	int i = usrbini(start.x);
	if (i<0) i = 0;	// correct for starting rounding error
	else
	if (i>=nx) i = nx-1;

	int j = usrbinj(start.y);
	if (j<0) j = 0;
	else
	if (j>=ny) j = ny-1;

	int k = usrbink(start.z);
	if (k<0) k = 0;
	else
	if (k>=nz) k = nz-1;

	// Moving variables
	int ptr = k*nynx + j*nx + i;
	double tx,    ty,    tz;		// distance to next intersection with x/y/z plane
	double tdx,   tdy,   tdz;		// distance increase for every x/y/z plane
	int    stepx, stepy, stepz;		// step direction -1/0/+1
	#define ptrx stepx
	int           ptry,  ptrz;		// ptr change ptrx=stepx
	int    outx,  outy,  outz;		// ending condition

	if (dir.x > SMALL) {
		tx    = (usrbinx(i+1) - start.x) / dir.x;
		tdx   = dx / dir.x;
		stepx = 1;
		outx  = nx;
	} else
	if (dir.x < -SMALL) {
		tx    = (usrbinx(i) - start.x) / dir.x;
		tdx   = -dx / dir.x;
		stepx = -1;
		outx  = -1;
	} else {
		tx    = INFINITE;
		tdx   = 0.0;
		stepx = 0;
		outx  = 0;
	}

	if (dir.y > SMALL) {
		ty    = (usrbiny(j+1) - start.y) / dir.y;
		tdy   = dy / dir.y;
		stepy = 1;
		ptry  = nx;
		outy  = ny;
	} else
	if (dir.y < -SMALL) {
		ty    = (usrbiny(j) - start.y) / dir.y;
		tdy   = -dy / dir.y;
		stepy = -1;
		ptry  = -nx;
		outy  = -1;
	} else {
		ty    = INFINITE;
		tdy   = 0.0;
		stepy = 0;
		ptry  = 0;
		outy  = 0;
	}

	if (dir.z > SMALL) {
		tz    = (usrbinz(k+1) - start.z) / dir.z;
		tdz   = dz / dir.z;
		stepz = 1;
		ptrz  = nynx;
		outz  = nz;
	} else
	if (dir.z < -SMALL) {
		tz    = (usrbinz(k) - start.z) / dir.z;
		tdz   = -dz / dir.z;
		stepz = -1;
		ptrz  = -nynx;
		outz  = -1;
	} else {
		tz    = INFINITE;
		tdz   = 0.0;
		stepz = 0;
		ptrz  = 0;
		outz  = 0;
	}

	// increase by initial starting length
	tx += t;
	ty += t;
	tz += t;
	double amount = 0.0;	// amount to score
	// Loop forever until we run out of length
	while (true) {
		//cout << "pos=" << i << ", " << j << ", " << k << " t=" << t << " total=" << totallength << endl;
		if (tx<ty && tx<tz) {
			i += stepx;
			if (i==outx) return amount;
			ptr += ptrx;
			if (tx < totallength) {
				double dt = tx - t;		// length to score
				data[ptr] += dt*value;
				amount    += dt*value;
				t = tx;
			} else {
				double dt = totallength - t;
				data[ptr] += dt*value;
				amount    += dt*value;
				return amount;
			}
			tx += tdx;
		} else
		if (ty<tz) {
			j += stepy;
			if (j==outy) return amount;
			ptr += ptry;
			if (ty < totallength) {
				double dt = ty - t;		// length to score
				data[ptr] += dt*value;
				amount    += dt*value;
				t = ty;
			} else {
				double dt = totallength - t;
				data[ptr] += dt*value;
				amount    += dt*value;
				return amount;
			}
			ty += tdy;
		} else {
			k += stepz;
			if (k==outz) return amount;
			ptr += ptrz;
			if (tz < totallength) {
				double dt = tz - t;		// length to score
				data[ptr] += dt*value;
				amount    += dt*value;
				t = tz;
			} else {
				double dt = totallength - t;
				data[ptr] += dt*value;
				amount    += dt*value;
				return amount;
			}
			tz += tdz;
		}
	}
} // add
