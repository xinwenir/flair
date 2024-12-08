//********************************************
// Author:  wioletta.kozlowska@cern.ch
// based on the T.Boehlen version from 2011
// Modified wioletta.kozlowska@cern.ch
// Version: 3.1
// Last change: 23/08/2016
//********************************************

#ifndef CELL_LINE_H
#define CELL_LINE_H

#include <string>

/*
 Note: General struture for saving
		LEM RBE Cell info
 Parameters: alpha, sqrt(beta),
			CutDose, SurivivalMax
			Name
 */
typedef struct {
	double nAlpha;
	double nBeta;
	double nCut;
	double nSmax;//maximum slope
	std::string  cellName;
} BIO_PARAMETERS;

/*
 Cell_Line class
 reeading the external files, saving biased paraemetrs
 */
class Cell_Line {
private:
	//store here the possibly biased parameters
	BIO_PARAMETERS BioXRay;

	//store here the nominal parameters
	BIO_PARAMETERS BioXRayNominal; //TODO: Is it necessary?

public:

	Cell_Line(unsigned int n=1,std::string fileName="");
	virtual ~Cell_Line() {};

	int readCellLine(std::string fileName);
	friend std::ostream& operator<<(std::ostream& out, const Cell_Line& data);

	double nAlpha() {return BioXRay.nAlpha;};
	double nBeta() {return BioXRay.nBeta;};
	double nCut()  {return BioXRay.nCut;};
	double nSmax() {return BioXRay.nSmax;};
	std::string cellName() {return BioXRay.cellName;};
	double nAlphaDivBeta() {return (BioXRay.nAlpha/BioXRay.nBeta);};		//helper variable

	double nAlphaNominal() {return BioXRayNominal.nAlpha;};
	double nBetaNominal() {return BioXRayNominal.nBeta;};
	double nCutNominal()  {return BioXRayNominal.nCut;};
	double nSmaxNominal() {return BioXRayNominal.nSmax;};
	std::string cellNameNominal() {return BioXRayNominal.cellName;};

private:
	inline double calcSmax() { //nSmax=nAlphaX+2.D0*nBetaX*DT
		return (BioXRay.nAlpha+2.*BioXRay.nBeta*BioXRay.nCut);
	};

	inline double calcNominalSmax() { //nSmax=nAlphaX+2.D0*nBetaX*DT
		return (BioXRayNominal.nAlpha+2.*BioXRayNominal.nBeta*BioXRayNominal.nCut);
	};
};
 #endif
