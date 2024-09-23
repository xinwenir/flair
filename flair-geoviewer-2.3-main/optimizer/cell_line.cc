//********************************************
// Author:  wioletta.kozlowska@cern.ch
// based on the T.Boehlen version from 2011
// Modified wioletta.kozlowska@cern.ch
// Version: 3.1
// Last change: 23/08/2016
//********************************************

#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "cell_line.h"

using namespace std;

Cell_Line::Cell_Line(unsigned int n, string fileName)
{
	switch (n) {
		//case 1: see below
		case 0 :
			//read from file
			readCellLine(fileName);
		case 2 :
			//HSG - FURUSAWA nAlphaX=0.313D0 Gy-1  nBetaX=0.0615D0 Gy-2   DT=10.5D0 Gy (adjusted to fit the exp data, as usual)
			BioXRayNominal.cellName = "HSG Furusawa";
			BioXRayNominal.nAlpha   = 0.313;
			BioXRayNominal.nBeta    = 0.0615;
			BioXRayNominal.nCut     = 10.5;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "HSG Furusawa";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			//nSmax=nAlphaX+2.D0*nBetaX*DT
			break;

		case 3 :
			//CHORDOMA (LEM-IV -- HIT version) nAlphaX= 0.0081d0 nBetaX = 0.0033d0 DT    = 30.0d0
			BioXRayNominal.cellName = "Chordoma LEM-IV";
			BioXRayNominal.nAlpha   = 0.0081;
			BioXRayNominal.nBeta    = 0.0033;
			BioXRayNominal.nCut     = 30.;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "Chordoma LEM-IV";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 4 :
			//LIVER Dt30 nAlphaX= 0.01d0 nBetaX = 0.00067d0 DT    = 30.0d0
			//nAlpha/nBeta = 15.0 +/- 2.0 Gy, nAlpha = 0.010 +/- 0.001 Gy (-1), T(d) = 128 +/- 12 day.
			//http://www.ncbi.nlm.nih.gov/pubmed/18262101
			BioXRayNominal.cellName = "Liver Dt30Gy";
			BioXRayNominal.nAlpha   = 0.01;
			BioXRayNominal.nBeta    = 0.00067;
			BioXRayNominal.nCut     = 30.;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "Liver Dt30Gy";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 5 :
			//PROSTATE Dt20 nAlphaX= 0.15d0 nBetaX = 0.0682d0 DT    = 20.0d0
			//nAlpha/nBeta = 2.2  Gy
			//info from S. Brons
			BioXRayNominal.cellName = "Prostate Dt20Gy";
			BioXRayNominal.nAlpha   = 0.15;
			BioXRayNominal.nBeta    = 0.0682;
			BioXRayNominal.nCut     = 20.;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "Prostate Dt20Gy";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 6 :
			//V79 nAlphax=0.129d0  nBetax  =0.049d0 dt   =4.0d0
			//nAlpha/nBeta = 2.63265  Gy
			BioXRayNominal.cellName = "V79 Dt4Gy";
			BioXRayNominal.nAlpha   = 0.129;
			BioXRayNominal.nBeta    = 0.049;
			BioXRayNominal.nCut     = 4.;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "V79 Dt4Gy";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 7 :
			//PROSTATE Dt20 nAlphaX= 0.15d0 nBetaX = 0.0682d0 DT    = 30.0d0
			//nAlpha/nBeta = 2.2  Gy
			//info from S. Brons
			BioXRayNominal.cellName = "Prostate Dt30Gy";
			BioXRayNominal.nAlpha   = 0.15;
			BioXRayNominal.nBeta    = 0.0682;
			BioXRayNominal.nCut     = 30.;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "Prostate Dt30Gy";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 8 :
			//LIVER Dt40
			//nAlphaX= 0.01d0 nBetaX = 0.00067d0 DT    = 40.0d0
			//nAlpha/nBeta = 15.0 +/- 2.0 Gy, nAlpha = 0.010 +/- 0.001 Gy (-1), T(d) = 128 +/- 12 day.
			//http://www.ncbi.nlm.nih.gov/pubmed/18262101
			BioXRayNominal.cellName = "Liver Dt40Gy";
			BioXRayNominal.nAlpha   = 0.01;
			BioXRayNominal.nBeta    = 0.00067;
			BioXRayNominal.nCut     = 40.;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "Liver Dt40Gy";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 9 :
			//Normal Tissue Dt20
			BioXRayNominal.cellName = "Normal Tissue Dt20Gy";
			BioXRayNominal.nAlpha   = 0.045;
			BioXRayNominal.nBeta    = 0.030;
			BioXRayNominal.nCut     = 20.;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "Normal Tissue Dt20Gy";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 10 :
			//Normal Tissue Dt20
			BioXRayNominal.cellName = "a_x=0.173,b_x=0.032,dt=15";
			BioXRayNominal.nAlpha   = 0.173;
			BioXRayNominal.nBeta    = 0.032;
			BioXRayNominal.nCut     = 15.;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "a_x=0.173,b_x=0.032,dt=15";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

			////////////////////////////////////////////////////
		case 20075 :
			BioXRayNominal.cellName = "HSG Furusawa";
			BioXRayNominal.nAlpha   = 0.313;
			BioXRayNominal.nBeta    = 0.0615;
			BioXRayNominal.nCut     = 10.5;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "HSG Furusawa nAlpha-25%";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha*0.75;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 20125 :
			BioXRayNominal.cellName = "HSG Furusawa";
			BioXRayNominal.nAlpha   = 0.313;
			BioXRayNominal.nBeta    = 0.0615;
			BioXRayNominal.nCut     = 10.5;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "HSG Furusawa nAlpha+25%";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha*1.25;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 21075 :
			BioXRayNominal.cellName = "HSG Furusawa";
			BioXRayNominal.nAlpha   = 0.313;
			BioXRayNominal.nBeta    = 0.0615;
			BioXRayNominal.nCut     = 10.5;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "HSG Furusawa nBeta-25%";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta*0.75;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 21125 :
			BioXRayNominal.cellName = "HSG Furusawa";
			BioXRayNominal.nAlpha   = 0.313;
			BioXRayNominal.nBeta    = 0.0615;
			BioXRayNominal.nCut     = 10.5;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "HSG Furusawa nBeta+25%";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta*1.25;
			BioXRay.nCut     = calcSmax();
			BioXRay.nSmax    = calcSmax();
			break;

		case 22002 : //HSG with nAlpha/nBeta=2
			//HSG - FURUSAWA nAlphaX=0.123D0 Gy-1  nBetaX=0.0615D0 Gy-2   DT=10.5D0 Gy (adjusted to fit the exp data, as usual)
			BioXRayNominal.cellName = "HSG with a/b=2";
			BioXRayNominal.nAlpha   = 0.123;
			BioXRayNominal.nBeta    = 0.0615;
			BioXRayNominal.nCut     = 10.5;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "HSG with a/b=2";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		case 22010 : //HSG with nAlpha/nBeta=10
			//HSG - FURUSAWA nAlphaX=0.615D0 Gy-1  nBetaX=0.0615D0 Gy-2   DT=10.5D0 Gy (adjusted to fit the exp data, as usual)
			BioXRayNominal.cellName = "HSG with a/b=10";
			BioXRayNominal.nAlpha   = 0.615;
			BioXRayNominal.nBeta    = 0.0615;
			BioXRayNominal.nCut     = 10.5;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "HSG with a/b=10";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
			break;

		default: //this is case 1;
			//CHORDOMA LEM-I nAlphaX= 0.1 Gy^-1, nBetaX=0.05Gy^-2  and Dt = 30 Gy !LEM I
			BioXRayNominal.cellName = "Chordoma LEM-I";
			BioXRayNominal.nAlpha   = 0.1;
			BioXRayNominal.nBeta    = 0.05;
			BioXRayNominal.nCut     = 30.;
			BioXRayNominal.nSmax    = calcNominalSmax();

			BioXRay.cellName = "Chordoma LEM-I";
			BioXRay.nAlpha   = BioXRayNominal.nAlpha;
			BioXRay.nBeta    = BioXRayNominal.nBeta;
			BioXRay.nCut     = BioXRayNominal.nCut;
			BioXRay.nSmax    = calcSmax();
	}
	//nAlphaDivBeta=BioXRay.nAlpha/BioXRay.nBeta;
} // Cell_Line

int Cell_Line::readCellLine(string fileName)
{
	string Line;
	string Dummy;

	cout<<"Reading Cell Line info from File "<<fileName<<endl;
	cout.flush();
	ifstream FileIn(fileName.c_str());

	if (!FileIn.good()) {
		cerr<<"\nERROR: Optimizer::readCellLine(): I/O error!\n"<<endl;
		exit(-1);	// FIXME throw an exception but not exit!!!
	}

	while(FileIn.good()) {
		//CHORDOMA LEM-I nAlphaX= 0.1 Gy^-1, nBetaX=0.05Gy^-2  and Dt = 30 Gy !LEM I
		FileIn>>Dummy>>	BioXRayNominal.cellName;
		FileIn>>Dummy>>	BioXRayNominal.nAlpha;
		FileIn>>Dummy>> BioXRayNominal.nBeta;
		FileIn>>Dummy>> BioXRayNominal.nCut;
		BioXRayNominal.nSmax=calcNominalSmax();

		BioXRay.cellName=BioXRayNominal.cellName;
		BioXRay.nAlpha=BioXRayNominal.nAlpha;
		BioXRay.nBeta=BioXRayNominal.nBeta;
		BioXRay.nCut=BioXRayNominal.nCut;
		BioXRay.nSmax=calcSmax();
	}

	//nAlphaDivBeta=BioXRay.nAlpha/BioXRay.nBeta;
	cout<< *this;
	return true;
} // readCellLine

ostream& operator<<(ostream& out, const Cell_Line& data)
{
	out<< "Bio Cell Line:\t"<<data.BioXRayNominal.cellName<<"\nAlpha value:\t"<<data.BioXRayNominal.nAlpha
	   <<"\nBeta value:\t"<<data.BioXRayNominal.nBeta<<"\nDose Cut:\t"<<data.BioXRayNominal.nCut
	   <<"\nMaxSurvival:\t"<<data.BioXRayNominal.nSmax<<endl;
	return out;
} // operator<<
