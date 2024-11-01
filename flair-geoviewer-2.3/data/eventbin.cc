// ********************************************
// Author:  wioletta.kozlowska@cern.ch
// based on the T.Boehlen version from 2011
// Modified wioletta.kozlowska@cern.ch
// Version: 3.1
// Last change: 19/09/2016
// *********************************************
//
#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "fortran.h"
#include "eventbin.h"

using namespace std;

int Eventbin::load(string fn, int /*det*/)
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Evetbin::load"<<endl;
	cout<<"---------------------------------------------------------"<<endl;
	filename = fn;
	file.open(filename.c_str(),"rb");
	if (!file) cerr << "ERROR: PBMatrix cannot open file " << filename << endl;
	if (readHeader())	cerr<<"ERROR: PBMatrix reading Eventbin"<<endl;
	return 0;
} // load

void Eventbin::cleanup()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Evenbin::clean"<<endl;
	cout<<"---------------------------------------------------------"<<endl;
	filename.clear();
	nx        = 0;
	ny        = 0;
	nz        = 0;
	xlow      = 0;
	ylow      = 0;
	zlow      = 0;
	nhits     = 0;
	weight    = 0;
	_detector = 0;
	file.close();
} // cleanup

int Eventbin::readHeader()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Eventbin::readHeader"<<endl;
	cout<<"---------------------------------------------------------"<<endl;
	char	runtitle[USRBIN_RUNTIT+1];
	char	runtime[USRBIN_RUNTIM+1];
	char	titusb[USRBIN_TITUSB+1];
	char	buffer[512];
	int	length;

	__attribute__((unused)) int	itusbn;		// binning type
	__attribute__((unused)) int	lntzer;
	__attribute__((unused)) float	bkusbn;
	__attribute__((unused)) float	b2usbn;
	__attribute__((unused)) float	tcusbn;

	// First buffer WRITE (1) RUNTIT*80,RUNTIM*32,SNGL(WCTOT),NCTOT
	length = file.read(buffer, sizeof(buffer));
	FortranParser parser(buffer, length);
	parser.read(runtitle, USRBIN_RUNTIT);
	parser.read(runtime,  USRBIN_RUNTIM);

	//length = fortranRead(file,buffer,sizeof(buffer)); //TTB for each detector read 1 fortran block
	length = file.read(buffer, sizeof(buffer));
	if (length<0) {
		fprintf(stderr,"Error reading input file !\n");
		exit(-1);	// FIXME
	}

	// parse buffer
	parser(buffer);

	_detector = parser.readInt();
	parser.read(titusb, USRBIN_TITUSB);
	itusbn   = parser.readInt();
	_score   = parser.readInt();
	printf("detector %i titusb %s itusb %i score %i \n", _detector, titusb, itusbn, _score);

	xlow  = parser.readFloat();
	xhigh = parser.readFloat();
	nx    = parser.readInt();
	dx    = parser.readFloat();
	printf(" xlow %f xhigh %f nx %i dx%f \n", xlow, xhigh, nx, dx);

	ylow  = parser.readFloat();
	yhigh = parser.readFloat();
	ny    = parser.readInt();
	dy    = parser.readFloat();
	printf(" ylow %f yhigh %f ny %i dy%f \n", ylow, yhigh, ny, dy);

	zlow  = parser.readFloat();
	zhigh = parser.readFloat();
	nz    = parser.readInt();
	dz    = parser.readFloat();
	printf(" zlow %f zhigh %f nz %i dz %f \n", zlow, zhigh, nz, dz);

	lntzer = parser.readInt();
	bkusbn = parser.readFloat();
	b2usbn = parser.readFloat();
	tcusbn = parser.readFloat();

	return 0;
} // readHeader

/*--------------------------------------------------------*/
char* Eventbin::readEvent(char *buffer, unsigned sizeBuffer)
{
	//READ EVENT HEADER
	//length = fortranRead(file,Buffer,SizeBuffer); //TTB for each detector read 1 fortran block
	int length = file.read(buffer,sizeBuffer);		//TTB for each detector read 1 fortran block

	if (length<0) {
		fprintf(stderr,"1:Error reading input file!\n");
		exit(-1);	// FIXME
	}

	// parse buffer
	FortranParser parser(buffer, length);

	__attribute__((unused))	int nbx,ncase;
	nbx   = parser.readInt();
	ncase = parser.readInt();		// ncases nb of rund cases
	float acase = parser.readFloat();	// acase -required in WHAT(4) eventbin
	int   mcase = parser.readInt();		// mcase??

	weight=acase+1.E+9*(float)mcase;
//	cout<<"DEBUG: nbx,"<<nbx<<" ncase, "<<ncase<<" acase, "<<acase<<" mcase "<<mcase<<endl;

	length = file.read(buffer,sizeBuffer); //TTB for each detector read 1 fortran block
	if (length<0) {
		cerr << "Error reading input file!" << endl;
		exit(-1);	// FIXME
	}

	parser(buffer, length);
	nhits = parser.readInt();

	return parser.getPtr();
} // readEvent
