// ********************************************
// Author:  wioletta.kozlowska@cern.ch
// based on the T.Boehlen version from 2011
// Modified wioletta.kozlowska@cern.ch
// Version: 3.1
// Last change: 24/08/2016
// *********************************************
#ifndef PBMATRIX_H
#define PBMATRIX_H

#include <string>
#include <vector>

class Roi;
class Eventbin;
class PencilBeam;

typedef struct {
	std::vector<int>   RoiId;            //Voxel index in ROI
	std::vector<float> Dose;             //Biological Dose [GyE] or Absorbed Dose [Gy]
	std::vector<float> Rbe;              //only used if isBiologicalOptimization
	std::vector<float> AlphaMeanDose;    //AlphaMean Dose
	std::vector<float> SqrtBetaMeanDose; //sqrtBeta dose from fluka
	std::vector<float> Let;              //Dose weighted LET
} DOSE_GRID;

typedef struct {
	int*   VoxelIndex;                   //Voxel index in DOSE_GRID
	float* DosePerPrimary;               //WARNING: in case of an error matrix
	                                     //the following entry is replaced by the respective standard error
	float* AlphaMean;                    //Alpha Mean Dose
	float* SqrtBetaMean;                 //Sqrt Beta Mean Dose
	float* LetMean;                      //Let

	//     number of pencil beam, tot number of hit voxels, tot number of events for this pencil beam,
	int PbIndex;				//ID of the PB
	int VoxelHit;				//nb. of voxels with energy deposition
	int VoxelHitAboveThreshold;	//nb. of voxels with energy deposition above threshold
	int Events;					//nb. of primaries used for PB simulation

	int TempVoxelHitBeta;
} PB_MATRIX;

class PBMatrix {
public:
	std::vector<std::vector<PB_MATRIX> >  box;
	std::vector< std::vector<PB_MATRIX> > BoxError;

public:
	PBMatrix() {}
//	PBMatrix(PencilBeam*, Roi* , std::string ,
//		 std::string, std::string, std::string,
//		 double, double, bool,
//		 bool, std::vector<DOSE_GRID>);

	~PBMatrix() { clean(); }

	void clean() {
			clean(box);
			clean(BoxError);
		}

	//TODO we do not need so many arguments
	int readBinaryMatrix(PencilBeam *PB, Roi *roi,
			     std::string vFileDose, std::string vFileAlpha,
			     std::string vFileBeta, std::string vFileLet,
			     double nMatrixMinDoseThreshold, double nRbeFixed,
			     bool isLetOptimization, bool isBiologicalOptimization,
			     std::vector<DOSE_GRID> fDoseGrid);//TODO: in dev mode!

	int readErrorMatrix(std::string FileName, int RoiBoxId,
			    PencilBeam *PB, Roi *roi,
			    double nMatrixMinDoseThreshold, double nRbeFixed);

private:
	int clean(std::vector< std::vector<PB_MATRIX> >& Struct);
}; // class PBMatrix

#endif
