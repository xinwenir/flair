//********************************************
// Author:  wioletta.kozlowska@cern.ch
// based on the T.Boeahlen version from 2011
// Modified wioletta.kozlowska@cern.ch
// March 2016 - code refactoring
// Version: 3.1
// Last change: 27/04/2016
//********************************************

#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

#include <math.h>
#include <stdlib.h>

#include "vector.h"
#include "histogram.h"
#include "optimizer.h"

#define VERSION "3.1"

using namespace std;
/*--------------------------------------------------------*/
void Optimizer::init(string in, string out, bool IsBioOpt, int ncells, string vCellFile)
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::constructor long"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	cout<<"###################################################"<<endl;
	cout<<"### Optimizer called  #############################"<<endl;
	cout<<"###################################################"<<endl;
	cout<<"  Version: "<<VERSION<<endl;
	cout<<"  Compilation date: "<<__DATE__<<" time: "<<__TIME__<<endl;
	cout<<"  Execution date & time: "<<CurrentDateTime()<<endl;
	cout<<"###################################################"<<endl;

	isInitialized=false;

	//*** Initialize all variables ***//

	fDataIn=in;
	fDataOut=out;
	nIterationMax=10000;

	PB.box.clear();//Set of pencil beams from different vPbFiles objects
	ROI.box.clear(); //set of roi
	Matrix.clean();

	isBiologicalOptimization=IsBioOpt;
	nOptimizationAlgorithm=3;
	CellLine=Cell_Line(ncells, vCellFile);

	nConvergenceCriterium=0.01;
	nScaleOptimizationStep=1.0;
	nMatrixMinDoseThreshold=0.0;
	nMinParticlesPerPb=0;

	isInitialEvaluation=false;
	isCheckPb=false;
	isLetOptimization=false;
	//Vodoo programming. It is not necessary since the container clear the content by itsels
	//TODO: Chek if in the containers are smart pointers to be cleaned (new)
	vFileDose.clear();
	vFileAlpha.clear();
	vFileBeta.clear();
	vFileLet.clear();
	vFileErrorMatrix.clear();
	vFileErrorMatrixRoiBoxId.clear();

	FileNewIntensities="";
	isReadNewIntensities=false;
	isReadNewIntensitiesBefore=false;
	nScaleIntensities=1.;

	nRbeFactor=1.;
	isRbeFixed=false;
	nRbeFixed=1.;
	nTargetRbeMean=0.;
	isRbeInPtvFixed=false;
	isPhotonEquivSurvival=false;
	isSelectiveEvaluation=false;

	nIterationActionInterval=20;//int(nIterationMax/50);

	//clear stuff
	CleanCost();

	cout<<"### Clean-up ######################################"<<endl;
	//output directory //dirty way
	//ensure that directory is existing
	string comm("mkdir -pv");
	  comm+=string(fDataOut);
	system(comm.c_str());
	//clean it
	comm="rm -v "+string(fDataOut)+"/*.res ";
	system(comm.c_str());
	cout<<"Remaining files which were not deleted:"<<endl;
	comm=string("ls ");
	  comm+=string(fDataOut);
	  comm+=string("/");
	system(comm.c_str());
	cout<<"###################################################"<<endl;
} // init

/*--------------------------------------------------------*/
int Optimizer::Clean(vector<DOSE_GRID> *Struct)
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::Clean"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	for (unsigned int i0=Struct->size();i0--;) {
		(*Struct)[i0].RoiId.clear();
		(*Struct)[i0].Dose.clear();
		(*Struct)[i0].Rbe.clear();
		(*Struct)[i0].AlphaMeanDose.clear();
		(*Struct)[i0].SqrtBetaMeanDose.clear();
	}
	Struct->clear();

	return 0;
} // Clean

/*--------------------------------------------------------*/
int Optimizer::CleanCost()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::CleanCost"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	fCostFunction.clear();
	fCostFunctionTime.clear();

	return 0;
} // CleanCost

/*
 * This is the main routine
 * Launches the initialization and optimization
 */
/*--------------------------------------------------------*/
int Optimizer::RunOptimizer()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::RunOptimizer"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	//create and clean dose grid
	InitializeDoseGrid();

	//get the matrix with the dose and bio-parameters for a given PB for each voxel in the ROI
	if (!vFileDose.empty()) {
		Matrix.readBinaryMatrix(&PB, &ROI, vFileDose, vFileAlpha, vFileBeta, vFileLet,
				nMatrixMinDoseThreshold, nRbeFixed, isLetOptimization,
				isBiologicalOptimization, fDoseGrid);
	} else {
		cerr<<"ERROR: Optimizer::RunOptimizer(): Dose File is empty. ";
	}

	//check sizes
	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		if (Matrix.box[iRoiId].size()!=PB.box.size()) {
			cerr<<"ERROR: Optimizer::RunOptimizer(): Matrix.box["<<iRoiId<<"].size()!=fPencilBeam.size(): "<<Matrix.box[iRoiId].size()<<" "<<PB.box.size()<<endl;// shall calculate actually the number of pencil beams? not the number of files
			exit(-1);	// FIXME
		}
	}

	// for (int i=0; i<PB.box.size();i++) {
	cout<<"BinaryMAtrixMaxdoseperprimary"<<PB.box[0].MaxDosePerPrimary<<endl;//}

	//now we have everything initiated from the files
	isInitialized=true;

	if (isReadNewIntensities&&isReadNewIntensitiesBefore) PB.readNewIntensities(FileNewIntensities); //read new intensities from a file BEFORE multiplying PBs

	if (isCheckPb) CheckPencilBeams(); //find out if peaked in a PTV or not and switch them off

	if ((nOptimizationAlgorithm==3)||(nOptimizationAlgorithm==4)||(nOptimizationAlgorithm==6)||(nOptimizationAlgorithm==8)) ComputeInverseMapping();//do index cross linking only if really needed (to safe memory)

	if (isReadNewIntensities&&!isReadNewIntensitiesBefore) PB.readNewIntensities(FileNewIntensities); //read new intensities from a file AFTER multiplying PBs

	if (fabs(nScaleIntensities-1.)>1.E-10) {//scale all particles numbers (to account for example for lateral scaling/occupation factor for 1D -> 3D)
		//apply when from 1D opt to 3D results: (2/1)^2=4 (i.e. dose per primary per area unit)
		PB.scaleIntensities(nScaleIntensities);
	}

	PrintInitial();

	if (isBiologicalOptimization&&isRbeFixed) {
		cerr<<"ERROR: Optimizer::RunOptimizer(): isBiologicalOptimization&&fIsRbeFixed: This is probably not what you want! Run Stopped!"<<endl;
		exit(-1);	// FIXME
	}

	//TODO: PreOptimize();//improve initial PB intensities

	Optimize();
	PrintFinal();

	//compute and plot statistical errors
	if (vFileErrorMatrix.size()>0) {
		cout<<"### Evaluation of Statistical Errors ##############"<<endl;

		//get the error matrix with the dose and bio-parameters for a given PB for each voxel in the ROI
		//  if (vFileMatrix.size()==0) {//read the standard file
		//	Matrix.ReadErrorMatrix(string(fDataIn+"/forwardpCT001_TPS_MATRIXERROR.TXT"),0,PB, ROI, nMatrixMinDoseThreshold, nRbeFixed);
		// } else {//read a file list
			for (unsigned int i0=0;i0<vFileErrorMatrix.size();i0++) {
				Matrix.readErrorMatrix(vFileErrorMatrix[i0],vFileErrorMatrixRoiBoxId[i0], &PB, &ROI,nMatrixMinDoseThreshold, nRbeFixed);
			}
			//}

		ComputeDoseError();
		Evaluate();

		cout<<"###################################################"<<endl;
	}

	//compute and plot equivalent photon survival
	if (isBiologicalOptimization&&isPhotonEquivSurvival) {
		cout<<"### Equivalent Photon Survival Evaluation #########"<<endl;
		ComputePhotonEquivSurvival();
		Evaluate();
		cout<<"###################################################"<<endl;
	}

	if (isSelectiveEvaluation) {
		//do some selective dose plotting with optimized results
		cout<<"### Selective Evaluation ##########################"<<endl;
		if (isBiologicalOptimization) {
			ComputeSelectiveBioDose();
		} else {
			ComputeSelectiveAbsDose();
		}
		Evaluate();
		cout<<"###################################################"<<endl;
	}
	return 0;
} // RunOptimizer

/*--------------------------------------------------------*/
int Optimizer::WriteDvh(string FileName, vector<string> *Tag, vector<int> *Type, vector<H1D*> *DvhD, vector<H1D*> *DvhI)
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::WriteDvh"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	cout<<"WRITE <"<<FileName<<"> ...";
	cout.flush();
	ofstream FileOut( FileName.c_str() ); //ios::trunc

	if (!FileOut.is_open()) {
		cerr<<"\nERROR: Optimizer::WriteDvh(): I/O error!\n"<<endl;
		exit(-1);	// FIXME
	}

	if (isBiologicalOptimization||isRbeFixed) {
		FileOut<<"# DVHs are in <RBE-weighted dose [Gy (RBE)]>. ";
	} else {
		FileOut<<"# DVHs are in <Dose [Gy]>. ";
	}
	FileOut<<"Next lines give <number of PTVs and OARs> and their <name tags>"<<endl;

	//in PENCIL.RES also the "total intensity" is given in the first line
	FileOut<<setw(8)<<DvhD->size()<<endl;// number of DVHs (i.e. number of PTVs and OARs) TODO
	for (unsigned int i0=0;i0<DvhD->size();i0++) {
		FileOut<<setw(8)<<Type->at(i0)<<" "<<Tag->at(i0)<<endl;
	}

	FileOut<<setw(10)<<"# Bin nb."<<setw(16)<<"Bin centre[Gy/Gy (RBE)]";

	for (unsigned int iRoi=0;iRoi<DvhD->size();iRoi++) {
		//DVH Volume(%) 2N columns with DoseDvhDif (%) and DoseDvhInt(%)
		FileOut<<setw(16)<<"diff. DVH[%]"<<setw(16)<<"int. DVH[%]";
	}
	FileOut<<endl;

	for (int iBin=1;iBin<=DvhD->at(0)->nbins();iBin++) {
		FileOut<<setw(10)<<iBin<<" "<<setw(16)<<DvhD->at(0)->center(iBin-1);
		for (unsigned int iRoi=0;iRoi<DvhD->size();iRoi++) {
			//DVH Volume(%) 2N columns with DoseDvhDif (%) and DoseDvhInt(%)
//		cout<<"\n  DvhBincenter "<<DvhD->at(0)->GetBinCenter(iBin)<<" "<< iBin;		cout<<" Dvhbincontent "<<DvhD->at(0)->GetBinContent(iBin);
		FileOut<<setw(16)<<DvhD->at(iRoi)->get(iBin-1)<<setw(16)<<DvhI->at(iRoi)->get(iBin-1);

		}
		FileOut<<endl;
	}

	FileOut.close();
	cout<<" DONE"<<endl;
	cout<<" Total number of PTVs and OARs for plan: "<<DvhD->size()<<endl;

	return 0;
} // WriteDvh

/*
 * Create dose grid
 */
/*--------------------------------------------------------*/
int Optimizer::InitializeDoseGrid()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::InitializeDoseGrid"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	//	cout<<"Initialize Dose Grid ...";
	cout.flush();

	for (unsigned int iRoiId=0;iRoiId<ROI.box.size();iRoiId++) {

		DOSE_GRID CurDoseGrid;
		CurDoseGrid.RoiId.clear();
		CurDoseGrid.Dose.clear();

		CurDoseGrid.Rbe.clear();
		CurDoseGrid.AlphaMeanDose.clear();
		CurDoseGrid.SqrtBetaMeanDose.clear();
		CurDoseGrid.Let.clear();

		cout<<" for Roi Id "<<ROI.box[iRoiId].RoiId<<" ... ";

		//do indexing
		for (int i0=0;i0<ROI.box[iRoiId].VoxelIndexMax+1;i0++) {
			CurDoseGrid.RoiId.push_back(-1); //initialize to no PTV or OAR assigned
			CurDoseGrid.Dose.push_back(0.);

			if (isBiologicalOptimization) {
				CurDoseGrid.Rbe.push_back(0.);
				CurDoseGrid.AlphaMeanDose.push_back(0.);
				CurDoseGrid.SqrtBetaMeanDose.push_back(0.);
			}
			if (isLetOptimization) {
				CurDoseGrid.Let.push_back(0.);
			}
		}

		//fill the RoiId to associate dose grid voxels with respective ROI
		for (unsigned int i0=0;i0<ROI.box[iRoiId].VoxelIndex.size();i0++) {
			int idx = ROI.box[iRoiId].VoxelIndex[i0];
			if ((idx<(signed int)CurDoseGrid.RoiId.size())&&(idx>0)) {
				// remove it as it is not suited for overlapping rois
				if (CurDoseGrid.RoiId[idx]==-1) {
					CurDoseGrid.RoiId[idx]=i0;
				} else {
					cerr<<"ERROR: Optimizer::InitializeDoseGrid(): Corrupt ROI-file! Double assignment of ROI! Voxel index: "<<idx<<endl;
					exit(-1);	// FIXME
				}
			} else {//index out-of-bounds
				cerr<<"ERROR: Optimizer::InitializeDoseGrid(): ROI entry: "<<i0<<": Index out-of-bounds: "<<idx<<" Maximum index: "<<CurDoseGrid.RoiId.size()<<endl;
				exit(-1);	// FIXME
			}
		}
		fDoseGrid.push_back(CurDoseGrid);
	}

	cout<<" DONE"<<endl;
	return 0;
} // InitializeDoseGrid

/*
 * reset dose grid to 0
 */
/*--------------------------------------------------------*/
int Optimizer::ResetDoseGrid()
{	//TODO MR: loop over ROI number

	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ResetDose Grid"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=ROI.box[iRoiId].VoxelIndexMax+1;i0--;) {
			fDoseGrid[iRoiId].Dose[i0]=0.;
			if (isBiologicalOptimization) {
				fDoseGrid[iRoiId].Rbe[i0]=0.;
				fDoseGrid[iRoiId].AlphaMeanDose[i0]=0.;
				fDoseGrid[iRoiId].SqrtBetaMeanDose[i0]=0.;
			}
		}
	}
	return 0;
} // ResetDoseGrid

/*
 * This has to be extra to compute the dose-weighted LET using the abs dose form the grid and resetting it
 */
/*--------------------------------------------------------*/
int Optimizer::ResetLetGrid()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ResetLetGrid"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (isLetOptimization) {
		for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
			for (unsigned int i0=ROI.box[iRoiId].VoxelIndexMax+1;i0--;) {
				fDoseGrid[iRoiId].Let[i0]=0.;
			}
		}
	} else {
		cerr<<"ERROR: Optimizer::ResetLetGrid(): Method called without having LET optimization activated!"<<endl;
		exit(-1);	// FIXME
	}

	return 0;
} // ResetLetGrid

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
int Optimizer::ComputeInverseMapping()
{	//TODO MR: loop over ROI number for both dosegrid and roi.matrixindex
	//if really necessary do the index cross linking
	//who cares, its just another huge array
	//cout<<"In ComputeInverseMapping function"<<vPbFiles[0].PB.BoxSet[0].Intensity.back()<<endl;
	//	cout<<"Compute inverse mapping ...";

	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeInverseMapping"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	cout.flush();

	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		cout<<" "<<iRoiId<<" ... "; cout.flush();

		for (unsigned int iVox=0;iVox<ROI.box[iRoiId].VoxelIndex.size();iVox++) {

			ROI.box[iRoiId].MatrixIndex.push_back(new int[PB.box.size()]) ;//fill Matrix voxel index for each PB
			//initialize
			for (unsigned int iPb=PB.box.size();iPb--;) {
				ROI.box[iRoiId].MatrixIndex[iVox][iPb]=-1;
			}
		}
			cout<<ROI.box[iRoiId].VoxelIndex.size()<<endl;
			for (unsigned int iPb=0;iPb<Matrix.box[iRoiId].size();iPb++) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					if (fDoseGrid[iRoiId].RoiId[Matrix.box[iRoiId][iPb].VoxelIndex[iVox]]>=0) {//voxel is a PTV or OAR voxel
						ROI.box[iRoiId].MatrixIndex[ fDoseGrid[iRoiId].RoiId[Matrix.box[iRoiId][iPb].VoxelIndex[iVox]] ][Matrix.box[iRoiId][iPb].PbIndex]=iVox;//fill Matrix voxel index
				}
			}
		}
	}

	cout<<" DONE"<<endl;
	return 0;
} // ComputeInverseMapping

/*--------------------------------------------------------*/
int Optimizer::CheckPencilBeams()
{	//TODO: Check with the ROI
	cout<<"### Check pencil beams ############################"<<endl;
	//check that max dose of a PB is inside a PTV
	int CountDeActivatedPb=0;
	for (unsigned int i0=PB.box.size();i0--;) {
		bool IsInPtv=false;
		//look in ROI with max dose if it is in a PTV
		for (unsigned int i1=ROI.box[PB.box[i0].MaxDoseRoiBoxId].VoxelIndex.size();i1--;) {
			if ((ROI.box[PB.box[i0].MaxDoseRoiBoxId].Id>0)&&(PB.box[i0].MaxDoseVoxelIndex==ROI.box[PB.box[i0].MaxDoseRoiBoxId].VoxelIndex[i1])) {
				IsInPtv=true;
				break;
			}
		}

		if (!IsInPtv) {//deactivate beams which peak not in a PTV
			cout<<"  PB #"<<i0<<" is not peaked in a PTV. Peak position is: ";
			if (PB.box[i0].MaxDoseVoxelIndex!=-1) {PrintVoxelIndex(PB.box[i0].MaxDoseVoxelIndex,PB.box[i0].MaxDoseRoiBoxId);}
			else{ cout<<"not known!"<<endl; }
			PB.box[i0].Active=false;
			PB.box[i0].Intensity[0]=0.;
			CountDeActivatedPb++;
		}
	}

	cout<<"  De-activated "<<CountDeActivatedPb<<" pencil beams which are not peaked in a PTV"<<endl;
	cout<<"###################################################"<<endl;

	return 0;
} // CheckPencilBeams

/*--------------------------------------------------------*/
int Optimizer::PreOptimize()
{
	cout<<"### Pre-Optimize ##################################"<<endl;
	cout<<" NOT DONE YET!"<<endl;

	PreOptimizeDRbe();

	cout<<"###################################################"<<endl;
	exit(-1);	// FIXME

	return 0;
} // PreOptimize

/*--------------------------------------------------------*/
int Optimizer::PreOptimizeDRbe()
{
	cout<<"### Pre-Optimize ##################################"<<endl;
	cout<<"  Back-to-front optimization of dose and RBE"<<endl;
	cout<<"###################################################"<<endl;
	exit(-1);	// FIXME

	return 0;
} // PreOptimizeDRbe

/*--------------------------------------------------------*/
int Optimizer::Optimize()
{
	//cout<<"In Optimize  function"<<vPbFiles[0].PB.BoxSet[0].Intensity.back()<<endl;

	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::Optimize"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::Optimize()"<<endl;
		exit(-1);	// FIXME
	}

	//initial values
	fCurrentIteration=0;//flag initial
	ComputeDose();
	ComputeCost();
	if (isInitialEvaluation) Evaluate();

	double CostChange=1.;
	unsigned int iLoop=0;
	cout<<"### Starting Optimization #########################"<<endl;
	cout <<setw(4)<<0<<setw(9)<<fCostFunctionTime[0]<<"s"<<" Cost: "<<setw(6)<<fCostFunction[0]<<endl;

	while((CostChange>nConvergenceCriterium)&&(iLoop<(unsigned int)nIterationMax)) {
		iLoop++;fCurrentIteration=iLoop;
		cout.flush();
		ComputeStep();//compute our next step (and check that minimum particle numbers per beam are fullfilled)
		ComputeDose();//compute dose, RBE, etc.
		ComputeCost();//compute cost function for given intensities

		CostChange=fabs((fCostFunction[iLoop]-fCostFunction[iLoop-1])/max(1.E-20,fCostFunction[iLoop-1]));
		if (fCostFunction[iLoop-1]-fCostFunction[iLoop]<0.) cout<<"### WARNING: Divergent cost function! ###"<<endl;
		cout <<setw(4)<<iLoop<<setw(9)<<fCostFunctionTime[iLoop]<<"s"<<" Cost: "<<setw(6)<<fCostFunction[iLoop]<<" Change: "<<CostChange*100<<"%"<<" Time: "<<fCostFunctionTime[iLoop]-fCostFunctionTime[iLoop-1]<<"s"<<endl;

		if (iLoop%nIterationActionInterval==0) {
			cout<<"### Current state ##################################"<<endl;
			Evaluate();//all
			cout<<"####################################################"<<endl;
		}

	}
	cout<<"### Optimization finished #########################"<<endl;

	fCurrentIteration=-1;//flag final

	return 0;
} // Optimize

/*
 * Computing the dose for given PB intensities
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeDose()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeDose"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeDose()"<<endl;
		exit(-1);	// FIXME
	}

	if (isBiologicalOptimization) {
		if (isLetOptimization) {//compute mean LET
			ComputeAbsDose();//this is need for the weighting of the dose-weighted LET
			ComputeLet();
		}
		ComputeBioDose();
	} else {
		ComputeAbsDose();
		if (isLetOptimization) {//compute mean LET
			ComputeLet();
		}
	}

	return 0;
} // ComputeDose

/*--------------------------------------------------------*/
int Optimizer::ComputeAbsDose()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeAbsDose"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeAbsDose()!"<<endl;
		exit(-1);	// FIXME
	}

	//initialize
	ResetDoseGrid();

	//  cout<<"DEBUG: Matrix.box.size(): "<<Matrix.box.size()<<" PB.box.size():"<<PB.box.size()<<endl;
	for (unsigned int iPb=PB.box.size();iPb--;) {
		if (PB.box[iPb].Active) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					fDoseGrid[iRoiId].Dose[ Matrix.box[iRoiId][iPb].VoxelIndex[iVox] ]+=
						Matrix.box[iRoiId][iPb].DosePerPrimary[iVox]*PB.box[iPb].Intensity.back();
				}
			}
		}

		if (PB.box[iPb].Intensity.back()<0.) {//just checking
			cerr<<"ERROR: Optimizer::ComputeAbsDose(): Particle number below zero: "<< PB.box[iPb].Intensity.back()<<" PB: "<<iPb<<endl;
			exit(-1);	// FIXME
		}
	}

	return 0;
} // ComputeAbsDose

/*--------------------------------------------------------*/
int Optimizer::ComputeBioDose()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeBioDose"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeBioDose()"<<endl;
		exit(-1);	// FIXME
	}

	//initialize
	ResetDoseGrid();

	//c           add here bio calc for each voxel:
	//c      avg_alpha(I)=alpha(JP,I)*DD*WJ(JP)+avg_alpha(I)
	//c in the end divide by sum DD*WJ(JP)
	//c      avg_beta(I)=sqrtbeta(JP,I)*DD*WJ(JP)+avg_beta(I)
	//c in the end divide by sum DD*WJ(JP) and then everything pow2

	//compute dose in voxel and mean alpha and sqrt(beta)
	for (unsigned int iPb=PB.box.size();iPb--;) {
		if (PB.box[iPb].Active) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					int idx= Matrix.box[iRoiId][iPb].VoxelIndex[iVox];
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];

					//calc abs. dose
					fDoseGrid[iRoiId].Dose[idx]+=DosePerPrimary *PB.box[iPb].Intensity.back();
					//calc dose-weighted average (see Krämer and Scholz 2006, Kanai et al. 1997, Zaider and Rossi 1980)
					//calc alpha mean in voxel
					fDoseGrid[iRoiId].AlphaMeanDose[idx]+=
							Matrix.box[iRoiId][iPb].AlphaMean[iVox]*DosePerPrimary*PB.box[iPb].Intensity.back();
					//calc sqrt beta mean in voxel
					fDoseGrid[iRoiId].SqrtBetaMeanDose[idx]+=
							Matrix.box[iRoiId][iPb].SqrtBetaMean[iVox]*DosePerPrimary*PB.box[iPb].Intensity.back();
				}
			}
		}

		if (PB.box[iPb].Intensity.back()<0.) {//just checking
			cerr<<"ERROR: Optimizer::ComputeBioDose(): Particle number below zero: "<< PB.box[iPb].Intensity.back()<<" PB: "<<iPb<<endl;
			exit(-1);	// FIXME
		}
	}

	const double nBeta=CellLine.nBeta();
	const double nAlphaDivBeta=CellLine.nAlphaDivBeta();

	//now compute RBE and RBE-weighted dose
	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=fDoseGrid[iRoiId].Dose.size();i0--;) {
			if (fDoseGrid[iRoiId].Dose[i0]>0.) {
				//TODO: introduce the Dt cut parameter (is this really necessary? then also need "slope max")
				//calculate survival
				double NegLnS_ion= fDoseGrid[iRoiId].AlphaMeanDose[i0]+pow(fDoseGrid[iRoiId].SqrtBetaMeanDose[i0],2);
				fDoseGrid[iRoiId].Rbe[i0]=(sqrt( NegLnS_ion/nBeta+pow(0.5*nAlphaDivBeta,2) )
						- 0.5*nAlphaDivBeta)/max(fDoseGrid[iRoiId].Dose[i0],(float)1.E-20);
				//assign RBE-weighted dose
				fDoseGrid[iRoiId].Rbe[i0]=pow(fDoseGrid[iRoiId].Rbe[i0],nRbeFactor);//do exponential RBE scaling
				fDoseGrid[iRoiId].Dose[i0]=fDoseGrid[iRoiId].Rbe[i0]*fDoseGrid[iRoiId].Dose[i0];
			} else if (fDoseGrid[iRoiId].Dose[i0]<0.) {
				cerr<<"ERROR: Optimizer::ComputeBioDose(): dose below zero: "<<fDoseGrid[iRoiId].Dose[i0]<<" RBE: "<<fDoseGrid[iRoiId].Rbe[i0]<<" Voxel: "<<i0<<endl;
				exit(-1);	// FIXME
			}
		}
	}
	return 0;
} // ComputeBioDose

//GO ON HERE: TEST ComputeLet and add Let plotting and then optimization

/*--------------------------------------------------------*/
// Computes dose-weighted average LET in each voxel
// WARNING: Assumes that the updated AbsDose is stored in the dose grid (and NO BioDose)!!!
//(i.e. if doing a bio-opt call first ComputeAbsDose then ComputeLet and then ComputeBioDose)
int Optimizer::ComputeLet()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeLet"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeLet()!"<<endl;
		exit(-1);	// FIXME
	}

	if (!isLetOptimization) {
		cerr<<"ERROR: Optimizer::ComputeLet(): called, but no LET optimization requested!"<<endl;
		exit(-1);	// FIXME
	}

	//initialize
	ResetLetGrid();

	//c           add here bio calc for each voxel:
	//c      avg_alpha(I)=alpha(JP,I)*DD*WJ(JP)+avg_alpha(I)
	//c in the end divide by sum DD*WJ(JP)
	//c      avg_beta(I)=sqrtbeta(JP,I)*DD*WJ(JP)+avg_beta(I)
	//c in the end divide by sum DD*WJ(JP) and then everything pow2

	//  cout<<"DEBUG: Matrix.box.size(): "<<Matrix.box.size()<<" fPencilBeam.size():"<<fPencilBeam.size()<<endl;
	for (unsigned int iPb=PB.box.size();iPb--;) {
		if (PB.box[iPb].Active) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					//dose-weighted LET average in voxel: SUM_AllPBs( LETFromPB*DoseFromPB ) /TotalDoseInVoxel
					int idx= Matrix.box[iRoiId][iPb].VoxelIndex[iVox];
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];

					if (fDoseGrid[iRoiId].Dose[idx]>0.) {
						fDoseGrid[iRoiId].Let[idx]+=
						Matrix.box[iRoiId][iPb].LetMean[iVox]
							* DosePerPrimary*PB.box[iPb].Intensity.back()
							/fDoseGrid[iRoiId].Dose[idx];
					}
				}
			}
		}

		if (PB.box[iPb].Intensity.back()<0.) {//just checking
			cerr<<"ERROR: Optimizer::ComputeLet(): Particle number below zero: "<< PB.box[iPb].Intensity.back()<<" PB: "<<iPb<<endl;
			exit(-1);	// FIXME
		}
	}
	return 0;
} // ComputeLet

/*
 * Computing the dose for given PB intensities
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeDoseError()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeDoseError"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeDoseError()"<<endl;
		exit(-1);	// FIXME
	}

	if (isBiologicalOptimization) {
		ComputeBioDoseError();
	} else {
		ComputeAbsDoseError();
	}

	return 0;
} // ComputeDoseError

/*--------------------------------------------------------*/
int Optimizer::ComputeAbsDoseError()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeAbsDoseError"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeAbsDoseError()!"<<endl;
		exit(-1);	// FIXME
	}

	cout<<"Calculate Absorbed Dose Error"<<endl;
//	ResetDoseGrid();
	fCurrentIteration=-4; //flag dose grid and plots!

	cout<<"Initialize Dose Error Grid ...";
	cout.flush();

	vector< vector<double> > DoseError;
	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		vector<double> CurDoseError;
		CurDoseError.clear();
		cout<<" "<<iRoiId<<" ... ";

		//do indexing
		for (unsigned int i0=0;i0<fDoseGrid[iRoiId].Dose.size();i0++) {
			CurDoseError.push_back(0.);
		}

		DoseError.push_back(CurDoseError);
	}
	cout<<" DONE"<<endl;

	//initialize
	ResetDoseGrid();

	//	double CurDose;
	//	double CurDoseError2;
	//	double DosePerPrimary;

	for (unsigned int iPb=PB.box.size();iPb--;) {//loop PB
		if (PB.box[iPb].Active) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {//loop ROI box
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {//loop on voxels hit by PB
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];
					double CurDose      = DosePerPrimary*PB.box[iPb].Intensity.back();
				//compute error^2 :
					double CurDoseError2 = (Matrix.BoxError[iRoiId][iPb].DosePerPrimary[iVox] - pow(DosePerPrimary,2) )/Matrix.BoxError[iRoiId][iPb].Events;
					fDoseGrid[iRoiId].Dose[ Matrix.box[iRoiId][iPb].VoxelIndex[iVox] ]+=CurDose;
					DoseError[iRoiId][ Matrix.box[iRoiId][iPb].VoxelIndex[iVox] ]+=CurDoseError2*pow(PB.box[iPb].Intensity.back(),2);
				}
			}
		}

		if (PB.box[iPb].Intensity.back()<0.) {//just checking
			cerr<<"ERROR: Optimizer::ComputeAbsDoseError(): Particle number below zero: "<< PB.box[iPb].Intensity.back()<<" PB: "<<iPb<<endl;
			exit(-1);	// FIXME
		}
	}

	//take sqrt of results and express as relative error (fraction) in the DOSEGRID!!!
	int CountPtvVoxels=0;
	double DoseMeanErrorInPTV(0.),DoseRmsErrorInPTV(0.),DoseMaxErrorInPTV(0.),DoseMinErrorInPTV(1000.);
	for (unsigned int iRoiId=(unsigned int)ROI.box.size();iRoiId--;) {
		for (unsigned int i0=fDoseGrid[iRoiId].Dose.size();i0--;) {
			fDoseGrid[iRoiId].Dose[ i0 ]=
				sqrt(DoseError[iRoiId][ i0 ])/max((double)fDoseGrid[iRoiId].Dose[ i0 ] ,1.E-30) ;
			//if voxel in a PTV compute average/min/max
			if (ROI.box[iRoiId].Id>0) {//a PTV
					DoseMeanErrorInPTV+=fDoseGrid[iRoiId].Dose[ i0 ];
					DoseRmsErrorInPTV+=pow(fDoseGrid[iRoiId].Dose[ i0 ],2);
					CountPtvVoxels++;
					//find min and max
					if (DoseMaxErrorInPTV<fDoseGrid[iRoiId].Dose[ i0 ])
						DoseMaxErrorInPTV=fDoseGrid[iRoiId].Dose[ i0 ];
					if (DoseMinErrorInPTV>fDoseGrid[iRoiId].Dose[ i0 ])
						DoseMinErrorInPTV=fDoseGrid[iRoiId].Dose[ i0 ];
			}
		}
	}
	DoseMeanErrorInPTV/=(double)CountPtvVoxels;
	DoseRmsErrorInPTV=sqrt(DoseRmsErrorInPTV)/(double)CountPtvVoxels;

	cout<<"### Statistical Uncertainty Abs. Dose in PTV ######"<<endl;
	cout<<"  Percentage relative to the actual abs. dose in the PTV: deltaDose/Dose*100."<<endl;
	cout<<"  Mean uncert. dose : "<<DoseMeanErrorInPTV*100.<<" %"<<endl;
	cout<<"  RMS of mean uncert. dose : "<<DoseRmsErrorInPTV*100.<<" %"<<endl;
	cout<<"  Min uncert. dose  : "<<DoseMinErrorInPTV*100.<<" %"<<endl;
	cout<<"  Max uncert. dose  : "<<DoseMaxErrorInPTV*100.<<" %"<<endl;
	cout<<"###################################################"<<endl;

	return 0;
} // ComputeAbsDoseError

/*--------------------------------------------------------*/
int Optimizer::ComputeBioDoseError()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeBioDoseError"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeBioDoseError()"<<endl;
		exit(-1);	// FIXME
	}

	cout<<"Calculate Biological Dose Error"<<endl;
//	ResetDoseGrid();
	fCurrentIteration=-4; //flag dose grid and plots!

	cout<<"Initialize Dose Error Grid ...";
	cout.flush();

	vector< vector<double> > DoseErrorSq,AlphaErrorSq,SqrtBetaErrorSq;
	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		vector<double> CurDoseError;
		CurDoseError.clear();

		cout<<" "<<iRoiId<<" ... ";

		//do indexing
		for (unsigned int i0=0;i0<fDoseGrid[iRoiId].Dose.size();i0++) {
			CurDoseError.push_back(0.);
		}

		DoseErrorSq.push_back(CurDoseError);//fill with dummies
		AlphaErrorSq.push_back(CurDoseError);
		SqrtBetaErrorSq.push_back(CurDoseError);
	}
	cout<<" DONE"<<endl;

	//initialize
	ResetDoseGrid();

	//c           add here bio calc for each voxel:
	//c      avg_alpha(I)=alpha(JP,I)*DD*WJ(JP)+avg_alpha(I)
	//c in the end divide by sum DD*WJ(JP)
	//c      avg_beta(I)=sqrtbeta(JP,I)*DD*WJ(JP)+avg_beta(I)
	//c in the end divide by sum DD*WJ(JP) and then everything pow2

	//compute dose in voxel and mean alpha and sqrt(beta)
	for (unsigned int iPb=PB.box.size();iPb--;) {
		if (PB.box[iPb].Active) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					int idx=Matrix.box[iRoiId][iPb].VoxelIndex[iVox];
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];

					//calc abs. dose
					fDoseGrid[iRoiId].Dose[idx]+=DosePerPrimary*PB.box[iPb].Intensity.back();
					//calc dose-weighted average (see Krämer and Scholz 2006, Kanai et al. 1997, Zaider and Rossi 1980)
					//calc alpha mean in voxel
					fDoseGrid[iRoiId].AlphaMeanDose[idx]+=
							Matrix.box[iRoiId][iPb].AlphaMean[iVox]*DosePerPrimary*PB.box[iPb].Intensity.back();
					//calc sqrt beta mean in voxel
					fDoseGrid[iRoiId].SqrtBetaMeanDose[idx]+=
							Matrix.box[iRoiId][iPb].SqrtBetaMean[iVox]*DosePerPrimary*PB.box[iPb].Intensity.back();

					//compute the respective sum of the error^2 quantities
					DoseErrorSq[iRoiId][idx] +=pow(PB.box[iPb].Intensity.back(),2)*
							( Matrix.BoxError[iRoiId][iPb].DosePerPrimary[iVox] - pow(DosePerPrimary,2) );

					AlphaErrorSq[iRoiId][idx]+=pow(PB.box[iPb].Intensity.back(),2)*
							( Matrix.BoxError[iRoiId][iPb].AlphaMean[iVox] - pow(Matrix.box[iRoiId][iPb].AlphaMean[iVox],2) );

					SqrtBetaErrorSq[iRoiId][idx]+=pow(PB.box[iPb].Intensity.back(),2)*
							( Matrix.BoxError[iRoiId][iPb].SqrtBetaMean[iVox] - pow(Matrix.box[iRoiId][iPb].SqrtBetaMean[iVox],2) );

				}
			}
		}

		if (PB.box[iPb].Intensity.back()<0.) {//just checking
			cerr<<"ERROR: Optimizer::ComputeBioDose(): Particle number below zero: "<< PB.box[iPb].Intensity.back()<<" PB: "<<iPb<<endl;
			exit(-1);	// FIXME
		}
	}

	int CountPtvVoxels=0;
	double DoseMeanErrorInPTV(0.),DoseRmsErrorInPTV(0.),DoseMaxErrorInPTV(0.),DoseMinErrorInPTV(1000.);
	const double nBeta=CellLine.nBeta();
	const double nAlphaDivBeta=CellLine.nAlphaDivBeta();

	//now compute RBE and RBE-weighted dose together with relative error
	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=fDoseGrid[iRoiId].Dose.size();i0--;) {
			if (fDoseGrid[iRoiId].Dose[i0]>0.) {
				//TODO: introduce the Dt cut parameter (is this really necessary? then also need "slope max")
				//calculate survival
				double NegLnS_ion= fDoseGrid[iRoiId].AlphaMeanDose[i0]+pow(fDoseGrid[iRoiId].SqrtBetaMeanDose[i0],2);
				fDoseGrid[iRoiId].Rbe[i0]=(sqrt( NegLnS_ion/nBeta+pow(0.5*nAlphaDivBeta,2) )
						- 0.5*nAlphaDivBeta)/max(fDoseGrid[iRoiId].Dose[i0],(float)1.E-20);
				//assign RBE-weighted dose
				fDoseGrid[iRoiId].Rbe[i0]=pow(fDoseGrid[iRoiId].Rbe[i0],nRbeFactor);//do exponential RBE scaling
				fDoseGrid[iRoiId].Dose[i0]=fDoseGrid[iRoiId].Rbe[i0]*fDoseGrid[iRoiId].Dose[i0];

				//now compute bio dose error
				double CurAbsDose=fDoseGrid[iRoiId].Dose[i0]/fDoseGrid[iRoiId].Rbe[i0];
				double Denominator=sqrt( NegLnS_ion/nBeta+pow(0.5*nAlphaDivBeta,2) );
				double BioDoseError=0.25*pow( ( fDoseGrid[iRoiId].AlphaMeanDose[i0]/CurAbsDose +
					2.*pow(fDoseGrid[iRoiId].SqrtBetaMeanDose[i0],2)/CurAbsDose )/
					Denominator,2) // (dBD/dAD)^2
					* DoseErrorSq[iRoiId][i0] // deltaAD^2
					+0.25*pow( CurAbsDose / Denominator ,2) // (dBD/dAlpha)^2
					*AlphaErrorSq[iRoiId][i0] // deltaAlpha^2
					+pow( fDoseGrid[iRoiId].SqrtBetaMeanDose[i0]*CurAbsDose/Denominator ,2) //(dBD/dSqrtBeta)^2
					*SqrtBetaErrorSq[iRoiId][i0]; // deltaSqrtBeta^2
				//store RELATIVE error of bio dose in bio dose array
				fDoseGrid[iRoiId].Dose[i0]=sqrt(BioDoseError)/max((double)fDoseGrid[iRoiId].Dose[ i0 ] ,1.E-30);

				//if voxel in a PTV compute average/min/max
				if (ROI.box[iRoiId].Id>0) {//a PTV
					DoseMeanErrorInPTV+=fDoseGrid[iRoiId].Dose[ i0 ];
					DoseRmsErrorInPTV+=pow(fDoseGrid[iRoiId].Dose[ i0 ],2);
					CountPtvVoxels++;
					//find min and max
					if (DoseMaxErrorInPTV<fDoseGrid[iRoiId].Dose[ i0 ])
						DoseMaxErrorInPTV=fDoseGrid[iRoiId].Dose[ i0 ];
					if (DoseMinErrorInPTV>fDoseGrid[iRoiId].Dose[ i0 ])
						DoseMinErrorInPTV=fDoseGrid[iRoiId].Dose[ i0 ];
				}
			} else if (fDoseGrid[iRoiId].Dose[i0]<0.) {
				cerr<<"ERROR: Optimizer::ComputeBioDoseError(): dose below zero: "<<fDoseGrid[iRoiId].Dose[i0]<<" RBE: "<<fDoseGrid[iRoiId].Rbe[i0]<<" Voxel: "<<i0<<endl;
				exit(-1); // FIXME
			}
			//    cout<<"DEBUG: DoseGrid:"<<i0<<" BioDose:"<<fDoseGrid[iRoiId].Dose[i0]<<"Gy (RBE) RBE:"<<fDoseGrid[iRoiId].Rbe[i0]<<endl;
		}
	}

	DoseMeanErrorInPTV/=(double)CountPtvVoxels;
	DoseRmsErrorInPTV=sqrt(DoseRmsErrorInPTV)/(double)CountPtvVoxels;

	cout<<"### Statistical Uncertainty Bio. Dose in PTV ######"<<endl;
	cout<<"  Percentage relative to the prescribed bio. dose in the PTV"<<endl;
	cout<<"  Mean uncert. dose : "<<DoseMeanErrorInPTV*100.<<" %"<<endl;
	cout<<"  RMS of mean uncert. dose : "<<DoseRmsErrorInPTV*100.<<" %"<<endl;
	cout<<"  Min uncert. dose  : "<<DoseMinErrorInPTV*100.<<" %"<<endl;
	cout<<"  Max uncert. dose  : "<<DoseMaxErrorInPTV*100.<<" %"<<endl;
	cout<<"###################################################"<<endl;

	return 0;
} // ComputeBioDoseError

/*
 * Calculate photon equivalent survival to assess misestimations in case of photons
 * to be used with care
 */
/*--------------------------------------------------------*/
int Optimizer::ComputePhotonEquivSurvival()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputePhotonEquivSurvival"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputePhotonEquivSurvival()"<<endl;
		exit(-1); // FIXME
	}

	cout<<"Calculate Photon Equivalent Survival"<<endl;
	cout<<"WARNING: Optimizer::ComputePhotonEquivSurvival(): in case you really want to calculate the `equivalent photon survival` you must use the NOMINAL ion alpha-beta values (in the FLUKA run) - NOT the varied ones - but give the varied photon parameters to the optimizer!"<<endl;

	fCurrentIteration=-3; //flag photon equivalent survival plots!

	//notice that by recalculating Dose and RBE here the selective plotting won't work for photon equivalent survival!
	//maybe it could work if one uses the doses calculated before, but we don't mess with this for the moment

	//initialize
	ResetDoseGrid();

	for (unsigned int iPb=PB.box.size();iPb--;) {
		if (PB.box[iPb].Active) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					int idx=Matrix.box[iRoiId][iPb].VoxelIndex[iVox];
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];					//calc abs. dose
					fDoseGrid[iRoiId].Dose[idx]+=DosePerPrimary
							*PB.box[iPb].Intensity.back();
					//calc dose-weighted average (see Krämer and Scholz 2006, Kanai et al. 1997, Zaider and Rossi 1980)
					//calc alpha mean in voxel
					fDoseGrid[iRoiId].AlphaMeanDose[idx]+=
							Matrix.box[iRoiId][iPb].AlphaMean[iVox]*DosePerPrimary*PB.box[iPb].Intensity.back();
					//calc sqrt beta mean in voxel
					fDoseGrid[iRoiId].SqrtBetaMeanDose[idx]+=
							Matrix.box[iRoiId][iPb].SqrtBetaMean[iVox]*DosePerPrimary*PB.box[iPb].Intensity.back();
				}
			}
		}
	}

	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=fDoseGrid[iRoiId].Dose.size();i0--;) {
			if (fDoseGrid[iRoiId].Dose[i0]>0.) {
				//assign RBE-weighted dose
				//TODO: introduce the Dt cut parameter (is this really necessary? then also need "slope max")
				double NegLnS_ion= fDoseGrid[iRoiId].AlphaMeanDose[i0]+pow(fDoseGrid[iRoiId].SqrtBetaMeanDose[i0],2); //compute -ln(S_ion)

				//compute the photon equivalent survival
				//now we want to verify what would happen if we have the same bio parameter variation for photons (i.e. alpha,beta,Dt) with the same dose distribution
				//for this we need to: calculate RBE with normal ion and photon values (so we have the same dose coverage in the PTV as for an equivalent photon treatment) - then with modified photon parameters we calculate survival
				//in other words: this calculates the biological dose (RBE-weighted in Gy (RBE)) as usual, then supposes it to be absorbed dose (in Gy) and calculates survival with the possibly modified PHOTON parameters! (alpha_photon, beta_photon)
				//use here the nominal alpha/beta NOT the possibly modified ones

				fDoseGrid[iRoiId].Rbe[i0]=(sqrt( NegLnS_ion/CellLine.nBetaNominal()+pow(0.5*CellLine.nAlphaNominal()/CellLine.nBetaNominal(),2) ) - 0.5*CellLine.nAlphaNominal()/CellLine.nBetaNominal())/max(fDoseGrid[iRoiId].Dose[i0],(float)1.E-20);

				fDoseGrid[iRoiId].Rbe[i0]=pow(fDoseGrid[iRoiId].Rbe[i0],nRbeFactor);//introduce here a possible RBE variation
				fDoseGrid[iRoiId].Dose[i0]=fDoseGrid[iRoiId].Rbe[i0]*fDoseGrid[iRoiId].Dose[i0];//store bio dose in the dose grid

				//now compute survival with the possibly modified photon alpha/beta
				//reverse the RBE calculation line to include also possible systematic misestimations from RBE, see above: fRbeFactor
				//use the BIOLOGIGAL DOSE stored in Dose[i0]
				//      Survival[i0]= exp(-CellLine.BioXRay.Alpha*Dose[i0] - CellLine.BioXRay.Beta*pow(Dose[i0],2));
				//TODO: CHECK do a check that NegLnS_ion == the new stuff if fRbeFactor=1
				//        Survival[i0]= exp(-NegLnS_ion);//this is the normal survival

			} else if (fDoseGrid[iRoiId].Dose[i0]<0.) {
				cerr<<"ERROR: Optimizer::ComputePhotonEquivSurvival(): dose below zero: "<<fDoseGrid[iRoiId].Dose[i0]<<" Voxel: "<<i0<<endl;
				exit(-1);	// FIXME
			}
		}
	}
	return 0;
} // ComputePhotonEquivSurvival

/*
 * The following methods specifies the selection of PBs for selective dose calculations
 */
/*--------------------------------------------------------*/
bool Optimizer::IsPbSelected(PENCIL_BEAM* Pb)
{
	//cout<<"In IsPbSelected function"<<vPbFiles[0].PB.BoxSet[0].Intensity.back()<<endl;
	if (Pb->Z==6) {
		return true;
	} else {
		return false;
	}
} // IsPbSelected

/*
 * The following versions can be used to calculate the dose for a restricted amount of PB.
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeSelectiveAbsDose()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeSelectiveAbsDose"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeSelectiveAbsDose()"<<endl;
		exit(-1);	// FIXME
	}

	cout<<"Calculate Selective Absorbed Dose"<<endl;
	ResetDoseGrid();
	fCurrentIteration=-2; //flag dose grid and plots!

	for (unsigned int iPb=PB.box.size();iPb--;) {
		if ((PB.box[iPb].Active)and (IsPbSelected(&PB.box[iPb]))) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					fDoseGrid[iRoiId].Dose[ Matrix.box[iRoiId][iPb].VoxelIndex[iVox] ]+=
						Matrix.box[iRoiId][iPb].DosePerPrimary[iVox]*PB.box[iPb].Intensity.back();
									}
			}
		}
	}
	return 0;
} // ComputeSelectiveAbsDose

/*
 * The following versions can be used to calculate the bio dose for a restricted amount of PB.
 * WARNING: for the moment we calculate the bio-dose as if only these PBs were present (and no additional dose contributions from other PB)
 * this is the way as one would use to verify the for example a bio experiment with separate beams
 * but it does not give the fractions of biodose directly from each field
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeSelectiveBioDose()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeSelectiveBioDose"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeSelectiveBioDose()!"<<endl;
		exit(-1);	// FIXME
	}

	cout<<"Calculate Selective Biological Dose"<<endl;
	fCurrentIteration=-2; //flag dose grid and plots!

	//initialize
	ResetDoseGrid();


	for (unsigned int iPb=PB.box.size();iPb--;) {
		if ((PB.box[iPb].Active) and (IsPbSelected(&PB.box[iPb]))) {//our condition
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					int idx=Matrix.box[iRoiId][iPb].VoxelIndex[iVox];
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];
					//calc abs. dose
					fDoseGrid[iRoiId].Dose[idx]+=DosePerPrimary*PB.box[iPb].Intensity.back();
					//calc dose-weighted average (see Krämer and Scholz 2006, Kanai et al. 1997, Zaider and Rossi 1980)
					//calc alpha mean in voxel
					fDoseGrid[iRoiId].AlphaMeanDose[idx]+=
						Matrix.box[iRoiId][iPb].AlphaMean[iVox]*DosePerPrimary*PB.box[iPb].Intensity.back();
					//calc sqrt beta mean in voxel
					fDoseGrid[iRoiId].SqrtBetaMeanDose[idx]+=
						Matrix.box[iRoiId][iPb].SqrtBetaMean[iVox]*DosePerPrimary*PB.box[iPb].Intensity.back();
				}
			}
		}
	}

	const double nBeta=CellLine.nBeta();
	const double nAlphaDivBeta=CellLine.nAlphaDivBeta();

	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=fDoseGrid[iRoiId].Dose.size();i0--;) {
			if (fDoseGrid[iRoiId].Dose[i0]>0.) {
				//assign RBE-weighted dose
				double NegLnS_ion= fDoseGrid[iRoiId].AlphaMeanDose[i0]+pow(fDoseGrid[iRoiId].SqrtBetaMeanDose[i0],2);
				fDoseGrid[iRoiId].Rbe[i0]=(sqrt( NegLnS_ion/nBeta+pow(0.5*nAlphaDivBeta,2) ) - 0.5*nAlphaDivBeta)/max(fDoseGrid[iRoiId].Dose[i0],(float)1.E-10);
				fDoseGrid[iRoiId].Rbe[i0]=pow(fDoseGrid[iRoiId].Rbe[i0],nRbeFactor);
				fDoseGrid[iRoiId].Dose[i0]=fDoseGrid[iRoiId].Rbe[i0]*fDoseGrid[iRoiId].Dose[i0];
			} else if (fDoseGrid[iRoiId].Dose[i0]<0.) {
				cerr<<"ERROR: Optimizer::ComputeSelectiveBioDose(): dose below zero: "<<fDoseGrid[iRoiId].Dose[i0]<<" Voxel: "<<i0<<endl;
				exit(-1);	// FIXME
			}
		}
	}
	return 0;
} // ComputeSelectiveBioDose

/*
 * Compute mean RBE from all PTV voxels
 */
/*--------------------------------------------------------*/
double Optimizer::ComputeRbeMean()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeRbeMean"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeRbeMean()!"<<endl;
		exit(-1);	// FIXME
	}

	double RbeMean=0;
	int CountVoxel=0;

	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=ROI.box[iRoiId].VoxelIndex.size();i0--;) {
			if (ROI.box[iRoiId].Id>0) {//a PTV
				RbeMean+=fDoseGrid[iRoiId].Rbe[ROI.box[iRoiId].VoxelIndex[i0]];
				CountVoxel++;
			}
		}
	}

	if (CountVoxel!=0) {
		RbeMean/=(double)CountVoxel;
	} else {
		cerr<<"ERROR: Optimizer::ComputeRbeMean(): CountVoxel==0! "<<RbeMean<<endl;
		exit(-1);	// FIXME
	}

	if (Verbose) cout<<"RBE_mean="<<RbeMean<<endl;

	return RbeMean;
} // ComputeRbeMean

/*
 * Compute <mean LET> from all PTV voxels
 * WARNING: it is <mean LET> over all voxels NOT <dose-weighted mean LET>!
 * This is supposed to be used only for converging to a constant <dose-weighted mean LET> in all PTV voxels, so the above fact shouldn't matter.
 */
/*--------------------------------------------------------*/
double Optimizer::ComputeLetMean()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeLetMean"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeLetMean()!"<<endl;
		exit(-1);	// FIXME
	}

	if (!isLetOptimization) {
		cerr<<"ERROR: Optimizer::ComputeLetMean(): called, but no LET optimization requested!"<<endl;
		exit(-1);	// FIXME
	}

	double LetMean=0;
	double CountVoxel=0;

//	cout<<"DEBUG: LET in voxels:";
	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=ROI.box[iRoiId].VoxelIndex.size();i0--;) {
			if (ROI.box[iRoiId].Id>0) {//a PTV
				LetMean+=fDoseGrid[iRoiId].Let[ROI.box[iRoiId].VoxelIndex[i0]];
				CountVoxel++;
			}
		}
	}

	if (CountVoxel!=0) {
		LetMean/=(double)CountVoxel;
	} else {
		cerr<<"ERROR: Optimizer::ComputeLetMean(): CountVoxel==0! "<<LetMean<<endl;
		exit(-1);	// FIXME
	}

	if (Verbose) cout<<"LET_voxel_mean="<<LetMean<<" (note: Is <mean LET> over all PTV voxels NOT <dose-weighted mean LET>!)"<<endl;

	return LetMean;
} // ComputeLetMean

/*--------------------------------------------------------*/
int Optimizer::ComputeCost()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeCost"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	double Chi2=0.;

	switch (nOptimizationAlgorithm) {
		case 2: Chi2=ComputeCostDoseRbe();
			break;
		case 3: Chi2=ComputeCostDose();
			break;
		case 4: Chi2=ComputeCostDoseRbe();
			break;
		case 5: Chi2=ComputeCostDoseRbe();
			break;
		case 6: Chi2=ComputeCostDoseRbe();
			break;
		case 7: Chi2=ComputeCostRbe();
			break;
		case 8: Chi2=ComputeCostDoseRbe();
			break;
		case 9: Chi2=ComputeCostDoseLet();
			break;
		default:
			Chi2=ComputeCostDose();//this is case 1;
	}

	fCostFunction.push_back(Chi2);
	fCostFunctionTime.push_back(((double) clock())/CLOCKS_PER_SEC);

	return 0;
} // ComputeCost

/*
 * Compute cost function based on (bio-)dose differences
 */
/*--------------------------------------------------------*/
double Optimizer::ComputeCostDose()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeCostDose"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeCostDose()!"<<endl;
		exit(-1);	// FIXME
	}

	double Chi2=0.;
	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=ROI.box[iRoiId].VoxelIndex.size();i0--;) {
			int idx = ROI.box[iRoiId].VoxelIndex[i0];
			double Weight = ROI.box[iRoiId].Weight[i0];

			if (ROI.box[iRoiId].Id>0) {//a PTV
				Chi2+= Weight*pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2)
						/pow(max(ROI.box[iRoiId].TargetDose,1.E-10),2);

			} else if (ROI.box[iRoiId].Id<0) {//OAR
				if (fDoseGrid[iRoiId].Dose[idx]>ROI.box[iRoiId].TargetDose) {//heaviside function
					Chi2+= Weight*pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2)
						/pow(max(ROI.box[iRoiId].TargetDose,1.E-10),2);
				}
			} else {
				cerr<<"ERROR: Optimizer::ComputeCostDose(): Entry: "<<i0<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
				exit(-1);	// FIXME
			}
		}
	}

	fCurrentChi2Dose=Chi2;
	fCurrentChi2Rbe=-1.;

	return Chi2;
} // ComputeCostDose

/*
 * Compute cost function based on (bio-)dose differences and RBE-difference (in PTVs only)
 * for dual goal optimization
 */
/*--------------------------------------------------------*/
double Optimizer::ComputeCostDoseRbe()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeCostDoseRbe"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeCostDoseRbe()!"<<endl;
		exit(-1);	// FIXME
	}
	if (!isBiologicalOptimization) {
		cerr<<"ERROR: Optimizer::ComputeCostDoseRbe(): To be called only for bio optimizations!"<<endl;
		exit(-1);	// FIXME
	}

	if (!isRbeInPtvFixed) nTargetRbeMean=ComputeRbeMean();

	double Chi2=0.;
	double Chi2_1=0.;//optimization goal 1
	double Chi2_2=0.;//optimization goal 2



	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=ROI.box[iRoiId].VoxelIndex.size();i0--;) {
			int idx = ROI.box[iRoiId].VoxelIndex[i0];
			double Weight = ROI.box[iRoiId].Weight[i0];
			//WARNING: QAD: set voxel-specific target RBE!
			//		nTargetRbeMean=ROI.box[iRoiId].TargetRbe; if (fCurrentIteration==1) cout<<"WARNING: QAD: set voxel-specific target RBE!"<<endl;
			//TODO: discriminate OAR and PTV
			if (ROI.box[iRoiId].Id>0) {//a PTV
				Chi2+= Weight*pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2) //chi2 for dose
					/pow(max(ROI.box[iRoiId].TargetDose,1.E-10),2);
				Chi2+= Weight*pow(fDoseGrid[iRoiId].Rbe[idx]-nTargetRbeMean,2)//chi2 for RBE
					/pow(max(nTargetRbeMean,1.E-10),2);

				Chi2_1+= Weight*pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2) //chi2 for dose
					/pow(max(ROI.box[iRoiId].TargetDose,1.E-10),2);
				Chi2_2+= Weight*pow(fDoseGrid[iRoiId].Rbe[idx]-nTargetRbeMean,2)//chi2 for RBE
					/pow(max(nTargetRbeMean,1.E-10),2);
			} else if (ROI.box[iRoiId].Id<0) {//OAR
				if (fDoseGrid[iRoiId].Dose[idx]>ROI.box[iRoiId].TargetDose) {//heaviside function
					Chi2+= Weight*pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2)
						/pow(max(ROI.box[iRoiId].TargetDose,1.E-10),2);
					Chi2_1+= Weight*pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2)
						/pow(max(ROI.box[iRoiId].TargetDose,1.E-10),2);
				}

			} else {
				cerr<<"ERROR: Optimizer::ComputeCostDoseRbe(): Entry: "<<i0<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
				exit(-1);	// FIXME
			}
					}
	}

	//  cout<<"Chi2_DoseRbe="<<Chi2<<endl;

	if (Verbose) cout<<"ComputeCostDoseRbe(): Chi2_DoseRbe="<<Chi2<<" Chi2_Dose="<<Chi2_1<<" Chi2_Rbe="<<Chi2_2<<" Current target mean RBE:"<<nTargetRbeMean<<" Actual mean RBE in PTV:"<<ComputeRbeMean()<<endl;

	fCurrentChi2Dose=Chi2_1;
	fCurrentChi2Rbe =Chi2_2;

	return Chi2;
} // ComputeCostDoseRbe

/*
 * Compute cost function based on the RBE-difference only (in PTVs only)
 */
/*--------------------------------------------------------*/
double Optimizer::ComputeCostRbe()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeCostRbe"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeCostRbe()!"<<endl;
		exit(-1);	// FIXME
	}

	if (!isBiologicalOptimization) {
		cerr<<"ERROR: Optimizer::ComputeCostRbe(): To be called only for bio optimizations!"<<endl;
		exit(-1);	// FIXME
	}

	if (!isRbeInPtvFixed) nTargetRbeMean=ComputeRbeMean();

	double Chi2=0.;

	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=ROI.box[iRoiId].VoxelIndex.size();i0--;) {
			int idx = ROI.box[iRoiId].VoxelIndex[i0];
			double Weight =ROI.box[iRoiId].Weight[i0];
			//TODO: discriminate OAR and PTV
			if (ROI.box[iRoiId].Id>0) {//a PTV
				Chi2+= Weight*pow(fDoseGrid[iRoiId].Rbe[idx]-nTargetRbeMean,2)//chi2 for RBE
					/pow(max(nTargetRbeMean,1.E-10),2);
			} else if (ROI.box[iRoiId].Id<0) {//OAR
				if (fDoseGrid[iRoiId].Dose[idx]>ROI.box[iRoiId].TargetDose) {//heaviside function
					//do nothing
				}
			} else {
				cerr<<"ERROR: Optimizer::ComputeCostRbe(): Entry: "<<i0<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
				exit(-1);	// FIXME
			}
		}
	}

	cout<<"ComputeCostRbe(): Chi2_Rbe="<<Chi2<<" Current target mean RBE:"<<nTargetRbeMean<<" Actual mean RBE in PTV:"<<ComputeRbeMean()<<endl;

	fCurrentChi2Dose=-1.;
	fCurrentChi2Rbe=Chi2;

	return Chi2;
} // ComputeCostRbe

/*
 * Compute cost function based on (bio-)dose differences and RBE-difference (in PTVs only)
 * for dual goal optimization
 */
/*--------------------------------------------------------*/
double Optimizer::ComputeCostDoseLet()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeCostDoseLet"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeCostDoseLet()!"<<endl;
		exit(-1);	// FIXME
	}

	if (!isLetOptimization) {
		cerr<<"ERROR: Optimizer::ComputeCostDoseLet(): To be called only for LET optimizations!"<<endl;
		exit(-1);	// FIXME
	}

	if (!isLetInPtvFixed) nTargetLetMean=ComputeLetMean();
	//  fLetMean=1.7;cout<<"WARNING: Fixed target mean LET of:"<<fTargetRbeMean<<endl;

	double Chi2=0.;
	double Chi2_1=0.;//optimization goal 1
	double Chi2_2=0.;//optimization goal 2


	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=ROI.box[iRoiId].VoxelIndex.size();i0--;) {
			int idx = ROI.box[iRoiId].VoxelIndex[i0];
			double Weight =ROI.box[iRoiId].Weight[i0];
			if (ROI.box[iRoiId].Id>0) {
				Chi2+= Weight*pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2) //chi2 for dose
						/pow(max(ROI.box[iRoiId].TargetDose,1.E-10),2);
				Chi2+= Weight*pow(fDoseGrid[iRoiId].Let[idx]-nTargetLetMean,2)//chi2 for LET
						/pow(max(nTargetLetMean,1.E-10),2);

				Chi2_1+= Weight*pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2) //chi2 for dose
						/pow(max(ROI.box[iRoiId].TargetDose,1.E-10),2);
				Chi2_2+= Weight*pow(fDoseGrid[iRoiId].Let[idx]-nTargetRbeMean,2)//chi2 for RBE
						/pow(max(nTargetLetMean,1.E-10),2);

			} else if (ROI.box[iRoiId].Id<0) {//OAR
				if (fDoseGrid[iRoiId].Dose[idx]>ROI.box[iRoiId].TargetDose) {//heaviside function
					Chi2+= Weight*pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2)
						/pow(max(ROI.box[iRoiId].TargetDose,1.E-10),2);

					Chi2_1+= Weight*pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2)
						/pow(max(ROI.box[iRoiId].TargetDose,1.E-10),2);

				//no prescribed LET goal for OAR possible right now!
				}

			} else {
				cerr<<"ERROR: Optimizer::ComputeCostDoseLet(): Entry: "<<i0<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
				exit(-1);	// FIXME
			}
		}
	}

	if (Verbose) cout<<"ComputeCostDoseLet(): Chi2_DoseLet="<<Chi2<<" Chi2_Dose="<<Chi2_1<<" Chi2_Let="<<Chi2_2<<" Current target mean LET:"<<nTargetLetMean<<" Actual mean LET in PTV:"<<ComputeLetMean()<<endl;

	fCurrentChi2Dose=Chi2_1;
	fCurrentChi2Let =Chi2_2;

	return Chi2;
		return 0;
} // ComputeCostDoseLet

/*--------------------------------------------------------*/
//cost like in FORTRAN programm FMCTPSDEV.f
double Optimizer::ComputeCostOld()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeCostOld"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeCostOld()!"<<endl;
		exit(-1);	// FIXME
	}

	double Chi2=0.;

	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=ROI.box[iRoiId].VoxelIndex.size();i0--;) {
			int idx = ROI.box[iRoiId].VoxelIndex[i0];
			Chi2+= pow(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose,2);
		}
	}

	cout<<"Chi2_Old="<<Chi2<<endl;

	fCurrentChi2Dose=Chi2;
	fCurrentChi2Rbe=-1.;

	return Chi2;
} // ComputeCostOld

/*--------------------------------------------------------*/
int Optimizer::ComputeStep()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeStep"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::ComputeStep()!"<<endl;
		exit(-1);	// FIXME
	}

	switch (nOptimizationAlgorithm) {
		case 2: ComputeStepDRbeD();	//dose-difference and RBE-difference TESTING!
			cout<<"WARNING: Optimization algorithm is still being tested and NOT FOR PRODUCTION!"<<endl;
			break;
		case 3: ComputeStepPg();	//plain-gradient: dose OK!
			break;
		case 4: ComputeStepPgRbe();	//plain-gradient: dose and RBE TESTING/OK!
			break;
		case 5: ComputeStepDRbeD2();	//dose-difference and RBE-difference TESTING/OK!
			//cout<<"WARNING: Optimization algorithm is still being tested and NOT FOR PRODUCTION!"<<endl;
			break;
		case 6 : ComputeStepPgRbe2();	//plain-gradient: dose and RBE TESTING!
			cout<<"WARNING: Optimization algorithm is still being tested and NOT FOR PRODUCTION!"<<endl;
			break;
		case 7: ComputeStepPgRbeOnly();	//plain-gradient: RBE in PTVs only OK!
			break;
		case 8: //hybrid PB //plain-gradient hybrid: dose and RBE TESTING!
			//compute step with the larger Chi2!
			//           if (fCurrentChi2Rbe>fCurrentChi2Dose) ComputeStepPgRbeOnly();
			//           else ComputeStepPG();
			//or define step:
			cout<<"WARNING: Optimization algorithm is still being tested and NOT FOR PRODUCTION!"<<endl;
			ComputeStepPgRbeOnly();
			ComputeBioDose();
			ComputeStepPg();
			break;
		case 9: ComputeStepDLetD();//dose-difference and LET-difference TESTING!
			cout<<"WARNING: Optimization algorithm is still being tested and NOT FOR PRODUCTION!"<<endl;
			break;
		default:
			ComputeStepD();//this is case 1; //dose-difference OK!
	}

	//discard PB with low particle numbers (if threshold is set)
	DiscardLowIntensityPb();

	return 0;
} // ComputeStep

/*
 * check if PB particle number falls below a given threshold, if so: set to zero
 * to respect lower accelerator limit
 */
/*--------------------------------------------------------*/
int Optimizer::DiscardLowIntensityPb()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::DiscardLowIntensityPb"<<endl;
	cout<<"---------------------------------------------------------"<<endl;


	for (unsigned int iPb=PB.box.size();iPb--;) {
		if (PB.box[iPb].Active) {
			//TODO: is it better for optimization to deactivate the zero-PB also?
			if (PB.box[iPb].Intensity.back()<nMinParticlesPerPb) {
				PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=0.;
			}
		}
	}

	//For taking into account the min threshold for particle: the value used @CNAO (for moment) for the proton beams is 5E5.
	//At HIT for carbon ion it is: 10000=1E4
	return 0;
} // DiscardLowIntensityPb

/*
 * Compute step with dose-difference (D) scaling
 * à la Lomax et al.
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeStepD()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeStepD"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	//compute the current gradient
	for (unsigned int iPb=PB.box.size();iPb--;) { //for each PB
		if (PB.box[iPb].Active) {
			double Up(0.),Down(0.);
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) { //for each voxel receiving dose by current PB
					int idx = Matrix.box[iRoiId][iPb].VoxelIndex[iVox];
					double Weight =ROI.box[iRoiId].Weight[iVox];
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];
					if ((fDoseGrid[iRoiId].RoiId[idx]>=0)//voxel is a PTV or OAR voxel
					   &&(fDoseGrid[iRoiId].Dose[idx]>0.)//current voxel dose>0
					   &&(ROI.box[iRoiId].TargetDose>0.)) {//target voxel dose>0
						if (ROI.box[iRoiId].Id>0) {//a PTV
							//sum_voxel(weight*DosePerParticle^2*wanted dose/current dose) / sum_voxel(weight*DosePerParticle^2)
							Up+=  Weight*pow(DosePerPrimary,2)
								*ROI.box[iRoiId].TargetDose/fDoseGrid[iRoiId].Dose[idx];
							Down+=Weight*pow(DosePerPrimary,2);

						} else if (ROI.box[iRoiId].Id<0) {//an OAR
							//OLD WITH ERROR ?!          } else if (ROI.box.Id[Matrix.box[iPb].VoxelIndex[iVox]]<0) {//an OAR
							if (fDoseGrid[iRoiId].Dose[idx]>ROI.box[iRoiId].TargetDose) {//heaviside function
								Up+=  Weight*pow(DosePerPrimary,2)
									*ROI.box[iRoiId].TargetDose/fDoseGrid[iRoiId].Dose[idx];
								Down+=Weight*pow(DosePerPrimary,2);
							}
						} else {
							cerr<<"ERROR: Optimizer::Optimizer::ComputeStepDD(): Entry: "<<iVox<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
							exit(-1);	// FIXME
						}
					}

				}//end of voxel loop
				//	cout<<"Pushing back intensity "<<PB.box[iPb].Intensity.back()<<" nScaleOptimizationStep* "<<(nScaleOptimizationStep*(Up/Down-1.)+1.)<<endl;
			}//end of RoiBoxId

			//assign new intensity
			if (Down>0.) {

				//rescale PB intensity:
				PB.box[iPb].Intensity.push_back( (nScaleOptimizationStep*(Up/Down-1.)+1.)*PB.box[iPb].Intensity.back() );

			} else {
				//			cout<<"WARNING: Optimizer::ComputeStepD(): PB: "<<iPb<<" Down: "<<Down<<", weights==0 OR DosePerPrimary==0, de-activate PB!"<<endl;
				PB.box[iPb].Active=false;
				PB.box[iPb].Intensity.push_back(0.);
			}
		} else {//PB  not active
			PB.box[iPb].Intensity.push_back(0.);
		}

	}// EO PB loop

	return 0;
} // ComputeStepD

/*
 * Compute step with the Plain Gradient (PG) method
 * for description see for instance A. Gemmel et al. 2008
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeStepPg()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeStepPg"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	double GradChi2[PB.box.size()];

	for (unsigned int iPb=PB.box.size();iPb--;) {

		GradChi2[iPb]=0.;

		if (PB.box[iPb].Active) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					int idx = Matrix.box[iRoiId][iPb].VoxelIndex[iVox];
					double Weight = ROI.box[iRoiId].Weight[iVox];
					float DosePerPrimary = Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];
					float TargetDose=ROI.box[iRoiId].TargetDose;//target voxel dose>0
					double Dose=fDoseGrid[iRoiId].Dose[idx]; //current voxel dose>0

					if ((fDoseGrid[iRoiId].RoiId[idx]>=0)//voxel is a PTV or OAR voxel
					   &&(Dose>0.) &&(TargetDose>0.)) {
						//now compute gradient of Chi2
						if (ROI.box[iRoiId].Id>0) {
							GradChi2[iPb]+=-2.*(Weight * (TargetDose - Dose)
									    / pow(TargetDose,2)
									    * DosePerPrimary
									    * ((isBiologicalOptimization)?fDoseGrid[iRoiId].Rbe[idx] :1.) //if bio opt: multiply by current RBE
									    );
						} else if ((ROI.box[iRoiId].Id<0) &&(Dose>TargetDose)) {//heavy func
								GradChi2[iPb]+=-2.*(Weight * (TargetDose - Dose)
										    / pow(TargetDose,2)
										    * DosePerPrimary
										    * ((isBiologicalOptimization)?fDoseGrid[iRoiId].Rbe[idx] :1.) //if bio opt: multiply by current RBE
										    );
						} else {
							cerr<<"ERROR: Optimizer::Optimizer::ComputeStepPG(): Entry: "<<iVox<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
							exit(-1);	// FIXME
						}
					}
				}//end of voxel loop
			}//end of RoiBoId loop
		}//EO PB  active
	}// EO PB loop

	//gradient is calculated - now compute the good stepsize
	double StepSizeUp=0.;
	double StepSizeDown=0.;

	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int iVoxX=ROI.box[iRoiId].VoxelIndex.size();iVoxX--;) {//for each ROI voxel
			double Weight=ROI.box[iRoiId].Weight[iVoxX];

			//According to Gemmel et al. loop ONLY over all PTV voxels. Instead we loop also OAR? -> No, can lead to divergences!
			// Why can't we add OARs? (besides the fact that they are not beeing defined for the point of the heaviside function and the point changes itself -> but can be approximated by old point of heaviside kick-in?) -> could lead to instabilities ... but not sure

			if (((ROI.box[iRoiId].Id>0))&&(ROI.box[iRoiId].TargetDose>0.) ) {//target voxel dose>0
				double RRR(0);
				for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
					if ((PB.box[iPb].Active)&&
						//altenatively one could change the loop structure and create large StepSizeUp + StepSizeDown arrays -> is probably better but needs to change sum squared!!!!
					   (ROI.box[iRoiId].MatrixIndex[iVoxX][iPb]>=0)) {//if dose deposited in ROI voxel by current PB
						RRR+= -GradChi2[iPb]*Matrix.box[iRoiId][iPb].DosePerPrimary[ROI.box[iRoiId].MatrixIndex[iVoxX][iPb]] * ((isBiologicalOptimization)? fDoseGrid[iRoiId].Rbe[ROI.box[iRoiId].VoxelIndex[iVoxX]]:1. ); //if bio opt: multiply by current RBE
					}//EO PB active
				}//EO PB loop
				//calculated without OAR:
				//CORRECT STEP-SIZE

				StepSizeUp+=  Weight//weight
					* (ROI.box[iRoiId].TargetDose - fDoseGrid[iRoiId].Dose[ROI.box[iRoiId].VoxelIndex[iVoxX]]) //(target dose - current dose)
					/ pow(ROI.box[iRoiId].TargetDose,2) // div by (target dose)^2
					* RRR;

				StepSizeDown+= Weight * pow( RRR / ROI.box[iRoiId].TargetDose ,2);  // weight*(R/ target dose )^2
			}//EO PTV only if
		}//EO ROI voxel loop
	}//EO RoiBoxId loop

	double StepSize(0.);
	if (StepSizeDown>0.) StepSize=StepSizeUp/StepSizeDown;

	//assign new intensities
	for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
		if (PB.box[iPb].Active) {
			//nScaleOptimizationStep: according to A.Gemmel for ~3Gy (RBE) -> 0.5 is good for <1Gy (RBE) use 0.25
			PB.box[iPb].Intensity.push_back(PB.box[iPb].Intensity.back()-GradChi2[iPb]*nScaleOptimizationStep*1.0*StepSize);
			//don't allow negative particle numbers
			if (PB.box[iPb].Intensity.back()<0.) {
				PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=0.;
			}
		} else {//PB  not active
			PB.box[iPb].Intensity.push_back(0.);
		}
	}// EO PB loop

	return 0;
} // ComputeStepPG

/*
 * Compute step with dose-difference (D) + RBE-difference scaling
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeStepDRbeD()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeStepDRbePg"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	double RbeMean=nTargetRbeMean;

	for (unsigned int iPb=PB.box.size();iPb--;) {
		if (PB.box[iPb].Active) {
			double Up(0.),Down(0.);
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					int idx = Matrix.box[iRoiId][iPb].VoxelIndex[iVox];
					double Weight=ROI.box[iRoiId].Weight[iVox];
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];
					if ((fDoseGrid[iRoiId].RoiId[idx]>=0)//voxel is a PTV or OAR voxel
					   &&(fDoseGrid[iRoiId].Dose[idx]>0.)//current voxel dose>0
					   &&(ROI.box[iRoiId].TargetDose>0.)//target voxel dose>0
					   &&(fDoseGrid[iRoiId].Rbe[idx]>0.)) {//RBE in voxel >0
						if (ROI.box[iRoiId].Id>0) {//a PTV
							//sum_voxel(weight*DosePerParticle^2*wanted dose/current dose) / sum_voxel(weight*DosePerParticle^2)
							Up+=  Weight*pow(DosePerPrimary,2)
								*ROI.box[iRoiId].TargetDose/fDoseGrid[iRoiId].Dose[idx]//dose-difference part
								//TODO: maybe a stronger RBE weighting is necessary e.g. squared
								// OR: decouple up and down in dose and RBE goal part: up_dose/down_dose*up_rbe/down_rbe
								//this would allow to monitor also the interference of the two goals
								//TODO: do a new function starting from this one to test some stuff -> why not use the normal chi2 part (see PG method) and define only the gradient differently?
								*RbeMean/fDoseGrid[iRoiId].Rbe[idx];//RBE-difference part
							Down+=Weight*pow(DosePerPrimary,2);
						} else if ((ROI.box[iRoiId].Id<0)&&//an OAR
							 (fDoseGrid[iRoiId].Dose[idx]>ROI.box[iRoiId].TargetDose)) {//heavy func
							Up+=  Weight*pow(DosePerPrimary,2)
								*ROI.box[iRoiId].TargetDose/fDoseGrid[iRoiId].Dose[idx];
							Down+=Weight*pow(DosePerPrimary,2);
						} else {
							cerr<<"ERROR: Optimizer::Optimizer::ComputeStepDDRbeD(): Entry: "<<iVox<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
							exit(-1);	// FIXME
						}
					}
				}//end of voxel loop
			}//end of RoiBoxId

			//assign new intensity
			if (Down>0.) {
				//rescale PB intensity:
				PB.box[iPb].Intensity.push_back( (nScaleOptimizationStep*(Up/Down-1.)+1.)*PB.box[iPb].Intensity.back() );
			} else {
				cout<<"WARNING: Optimizer::ComputeStepDDRbeD(): PB: "<<iPb<<" Down: "<<Down<<", keep old PB intensity!"<<endl;
				PB.box[iPb].Intensity.push_back(PB.box[iPb].Intensity.back());
			}
		} else {//PB  not active
			PB.box[iPb].Intensity.push_back(0.);
		}
	}
	return 0;
} // ComputeStepDRbePg

/*
 * Compute step with the Plain Gradient method
 * for dose and RBE
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeStepPgRbe()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeStepPgRbe"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	double RbeMean=nTargetRbeMean;
	double GradChi2[PB.box.size()];

	const double nBeta=CellLine.nBeta();
	const double nAlpha=CellLine.nAlpha();

	//compute the current gradient
	for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
		GradChi2[iPb]=0.;
		if (PB.box[iPb].Active) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {//for each voxel receiving dose by current PB
					double Weight=ROI.box[iRoiId].Weight[iVox];
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];
					int idx=Matrix.box[iRoiId][iPb].VoxelIndex[iVox];

					if ((fDoseGrid[iRoiId].RoiId[idx]>=0)//voxel is a PTV or OAR voxel
					   &&(fDoseGrid[iRoiId].Dose[idx]>0.)//current voxel dose>0
					   &&(ROI.box[iRoiId].TargetDose>0.)) {//target voxel dose>0
						//now compute gradient of Chi2
						if (ROI.box[iRoiId].Id>0) {//a PTV
							//compute grad(biodose) and grad(RBE) for a  VOXEL and PB
							double NegLnS_ion= fDoseGrid[iRoiId].AlphaMeanDose[idx]
								+pow(fDoseGrid[iRoiId].SqrtBetaMeanDose[idx],2);
							double GradBioDose= ( Matrix.box[iRoiId][iPb].AlphaMean[iVox]*DosePerPrimary
									      + 2.*fDoseGrid[iRoiId].SqrtBetaMeanDose[idx]
									      * Matrix.box[iRoiId][iPb].SqrtBetaMean[iVox]*DosePerPrimary
									      ) / ( 2. *nBeta*sqrt( NegLnS_ion/nBeta + pow(nAlpha/2./nBeta,2) ) );

							double GradRbe;
							if (fDoseGrid[iRoiId].Dose[idx]>0.) {
								GradRbe = ( GradBioDose - fDoseGrid[iRoiId].Rbe[idx] * DosePerPrimary )
									/ fDoseGrid[iRoiId].Dose[idx] * fDoseGrid[iRoiId].Rbe[idx];//the whole line is: abs dose
							} else {
								GradRbe = 0.; //assume zero for the zero-dose limiting case the initial RBE is given by RBE_alpha=alpha_ion/alpha_photon
							}
							//dose part
							GradChi2[iPb]+=-2.*( Weight //weight
									     * (ROI.box[iRoiId].TargetDose
										- fDoseGrid[iRoiId].Dose[idx])
									     / pow(ROI.box[iRoiId].TargetDose,2)
									     *( DosePerPrimary
										* ((isBiologicalOptimization)?fDoseGrid[iRoiId].Rbe[idx] :1.) //if bio opt: multiply by current RBE
							 //Here one could add the Grad(RBE) term: Grad(BioDose)=RBE*Grad(Dose) + Dose*Grad(RBE)
							 //but don't since according to A. Gemmel et al. 2008 not necessary
							 + fDoseGrid[iRoiId].Dose[idx] * GradBioDose ) //in case
							);

							//RBE part
							GradChi2[iPb]+=-2.*(Weight
							     * (RbeMean - fDoseGrid[iRoiId].Rbe[idx])
							     / pow(RbeMean,2) // div by (target RBE)^2
							     * GradRbe //grad(RBE)
							);
						} else if ((ROI.box[iRoiId].Id<0)&&//an OAR
							 (fDoseGrid[iRoiId].Dose[idx]>ROI.box[iRoiId].TargetDose)) {//heavy func
								GradChi2[iPb]+=-2.*(Weight
										    * (ROI.box[iRoiId].TargetDose
										       - fDoseGrid[iRoiId].Dose[idx])
										    / pow(ROI.box[iRoiId].TargetDose,2)
										    * DosePerPrimary
										    * ((isBiologicalOptimization)?fDoseGrid[iRoiId].Rbe[idx] :1.) //if bio opt: multiply by current RBE
								  //Here one could add the Grad(RBE) term: Grad(BioDose)=RBE*Grad(Dose) + Dose*Grad(RBE)
								  //but don't since according to A. Gemmel et al. 2008 not necessary
								  //                               * GradBioDose//in case
								);
						} else {
							cerr<<"ERROR: Optimizer::Optimizer::ComputeStepPGRBE(): Entry: "<<iVox<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
							exit(-1);	// FIXME
						}
					}
				}//end of voxel loop
			}//end of RoiBoxId
		}//EO PB  active
	}// EO PB loop

	//TODO: Our step is not correctly calculated yet!!!
	//We should calculate taking into account also the RBE gradient??? Or what?

	//gradient is calculated - now compute the good stepsize
	double StepSizeUp=0.;
	double StepSizeDown=0.;

	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int iVoxX=ROI.box[iRoiId].VoxelIndex.size();iVoxX--;) {//for each ROI voxel
			if ((ROI.box[iRoiId].Id>0)//a PTV, loop ONLY over all PTV voxels
					&&(ROI.box[iRoiId].TargetDose>0.) ) {//target voxel dose>0
				double RRR(0);
				double Weight = ROI.box[iRoiId].Weight[iVoxX];
				for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[ROI.box[iRoiId].MatrixIndex[iVoxX][iPb]];
					if ((PB.box[iPb].Active)&&
					   (ROI.box[iRoiId].MatrixIndex[iVoxX][iPb]>=0)) {//if dose deposited in ROI voxel by current PB
						RRR+= -GradChi2[iPb]*DosePerPrimary
							* ((isBiologicalOptimization)? fDoseGrid[iRoiId].Rbe[ROI.box[iRoiId].VoxelIndex[iVoxX]]:1. ); //if bio opt: multiply by current RBE
					}//EO PB active
				}//EO PB loop
				//calculated without OAR:
				//CORRECT STEP-SIZE
				StepSizeUp+= Weight
					* (ROI.box[iRoiId].TargetDose - fDoseGrid[iRoiId].Dose[ROI.box[iRoiId].VoxelIndex[iVoxX]]) //(target dose - current dose)
					/ pow(ROI.box[iRoiId].TargetDose,2) // div by (target dose)^2
					* RRR;
				StepSizeDown+= Weight * pow( RRR / ROI.box[iRoiId].TargetDose ,2);  // weight*(R/ target dose )^2
			}//EO PTV only if
		}//EO ROI voxel loop
	}//EO RoiBoxId loop

	double StepSize(0.);
	if (StepSizeDown>0.) StepSize=StepSizeUp/StepSizeDown;

	//assign new intensities
	for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
		if (PB.box[iPb].Active) {
			//nScaleOptimizationStep: according to A.Gemmel for ~3Gy (RBE) -> 0.5 is good for <1Gy (RBE) use 0.25
			PB.box[iPb].Intensity.push_back(PB.box[iPb].Intensity.back()-GradChi2[iPb]*nScaleOptimizationStep*1.0*StepSize);
			//don't allow negative particle numbers
			if (PB.box[iPb].Intensity.back()<0.) {
				PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=0.;
			}

		} else {//PB  not active
			PB.box[iPb].Intensity.push_back(0.);
		}
	}// EO PB loop
	return 0;
} // ComputeStepPgRbe

/*
 * TEST Compute step with dose-difference (DD) + RBE-difference scaling
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeStepDRbeD2()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeStepDRbeD2"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	double RbeMean=nTargetRbeMean;
	double RbeRelImportance=1.;

	const double nBeta=CellLine.nBeta();
	const double nAlphaDivBeta=CellLine.nAlphaDivBeta();

	for (unsigned int iPb=PB.box.size();iPb--;) {
		if (PB.box[iPb].Active) {

			double OtherUp(0.),OtherDown(0.);
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					double Weight=ROI.box[iRoiId].Weight[iVox];
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];
					int idx=Matrix.box[iRoiId][iPb].VoxelIndex[iVox];
					if ((fDoseGrid[iRoiId].RoiId[idx]>=0)//voxel is a PTV or OAR voxel
							&&(fDoseGrid[iRoiId].Dose[idx]>0.)//current voxel dose>0
							&&(ROI.box[iRoiId].TargetDose>0.)//target voxel dose>0
							&&(fDoseGrid[iRoiId].Rbe[idx]>0.)) {//RBE in voxel >0
						//WARNING: QAD: set voxel-specific target RBE!
						RbeMean=ROI.box[iRoiId].TargetRBE; if (fCurrentIteration==1) cout<<"WARNING: QAD: use voxel-specific target RBE for stepping!"<<endl;
						//compute absorbed dose in current voxel
						double CurrentDose=fDoseGrid[iRoiId].Dose[idx]/max(fDoseGrid[iRoiId].Rbe[idx],(float)1.E-20);
						//compute RBE for current PB assuming that all current dose is only by this beam
						double NegLnS_ion= Matrix.box[iRoiId][iPb].AlphaMean[iVox]*CurrentDose
							+ pow(Matrix.box[iRoiId][iPb].SqrtBetaMean[iVox]*CurrentDose,2);
						double RbePb=(sqrt( NegLnS_ion/nBeta+pow(0.5*nAlphaDivBeta,2) ) - 0.5*nAlphaDivBeta)/CurrentDose;

						if (ROI.box[iRoiId].Id>0) {
							OtherUp+=Weight*(
									 (RbeMean-fDoseGrid[iRoiId].Rbe[idx])
									 * (RbePb-RbeMean)*DosePerPrimary
									 *(RbeRelImportance*fCurrentChi2Rbe/(fCurrentChi2Dose+RbeRelImportance*fCurrentChi2Rbe))
									 +(fCurrentChi2Dose/(fCurrentChi2Dose+RbeRelImportance*fCurrentChi2Rbe))*
									 DosePerPrimary*RbePb*
									 (ROI.box[iRoiId].TargetDose//Bio dose part
									  - fDoseGrid[iRoiId].Dose[idx])

									 )*DosePerPrimary*RbePb*PB.box[iPb].Intensity.back()//damping factor
								/fDoseGrid[iRoiId].Dose[idx]
								;
							OtherDown+=Weight *(
									    pow(DosePerPrimary*RbePb,2) //dose part
									    *(fCurrentChi2Dose/(fCurrentChi2Dose+RbeRelImportance*fCurrentChi2Rbe))
									    +pow((RbePb-RbeMean)*DosePerPrimary,2) //RBE part
									    *(RbeRelImportance*fCurrentChi2Rbe/(fCurrentChi2Dose+RbeRelImportance*fCurrentChi2Rbe))
									    );

						} else if ((ROI.box[iRoiId].Id<0)&&
							 (fDoseGrid[iRoiId].Dose[idx]>ROI.box[iRoiId].TargetDose)) {
							OtherUp+=Weight*(
									 (fCurrentChi2Dose/(fCurrentChi2Dose+RbeRelImportance*fCurrentChi2Rbe))*
									 DosePerPrimary*RbePb* //Bio dose part
									 (ROI.box[iRoiId].TargetDose
									  - fDoseGrid[iRoiId].Dose[idx])
									 )*DosePerPrimary*RbePb*PB.box[iPb].Intensity.back()//damping factor
								/fDoseGrid[iRoiId].Dose[idx]
								;
							OtherDown+=Weight *(
									    (fCurrentChi2Dose/(fCurrentChi2Dose+RbeRelImportance*fCurrentChi2Rbe))*
									    pow(DosePerPrimary*RbePb,2) //bio dose part
									    );
						} else {
							cerr<<"ERROR: Optimizer::Optimizer::ComputeStepDDRbeD2(): Entry: "<<iVox<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
							exit(-1);	// FIXME
						}
					}

				}//end of voxel loop
			}//end of RoiBoxId

			if (OtherDown>0.) {
				PB.box[iPb].Intensity.push_back( PB.box[iPb].Intensity.back()+nScaleOptimizationStep*OtherUp/OtherDown );
			} else {
				cout<<"WARNING: Optimizer::ComputeStepDDRbeD2(): PB: "<<iPb<<" OtherUp: "<<OtherUp<<" OtherDown: "<<OtherDown<<", keep old PB intensity!"<<endl;
				PB.box[iPb].Intensity.push_back(PB.box[iPb].Intensity.back());
			}
		} else {//PB  not active
			PB.box[iPb].Intensity.push_back(0.);
		}
	}
	return 0;
} // ComputeStepDDRbeD2

/*
 * TEST Compute step with the Plain Gradient method
 * for dose and RBE
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeStepPgRbe2()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeStepPgRbe2"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (Verbose) cout<<"Enter Optimizer::ComputeStepPGRBE2()"<<endl;

	double RbeMean=nTargetRbeMean;

	double GradChi2[PB.box.size()];

	const double nBeta = CellLine.nBeta();
	const double nAlpha= CellLine.nAlpha();

	//compute the current gradient
	for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
		GradChi2[iPb]=0.;

		if (PB.box[iPb].Active) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {//for each voxel receiving dose by current PB
					double Weight=ROI.box[iRoiId].Weight[iVox];
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];
					int idx=Matrix.box[iRoiId][iPb].VoxelIndex[iVox];

					if ((fDoseGrid[iRoiId].RoiId[idx]>=0)//voxel is a PTV or OAR voxel
					   &&(fDoseGrid[iRoiId].Dose[idx]>0.)//current voxel dose>0
					   &&(ROI.box[iRoiId].TargetDose>0.)) {//target voxel dose>0
						//now compute gradient of Chi2
						if (ROI.box[iRoiId].Id>0) {//a PTV
							//compute grad(biodose) and grad(RBE) for a given VOXEL and PB
							double NegLnS_ion= fDoseGrid[iRoiId].AlphaMeanDose[idx]
								+pow(fDoseGrid[iRoiId].SqrtBetaMeanDose[idx],2);
							double GradBioDose= ( Matrix.box[iRoiId][iPb].AlphaMean[iVox]*DosePerPrimary
									      + 2.*fDoseGrid[iRoiId].SqrtBetaMeanDose[idx]
									      * Matrix.box[iRoiId][iPb].SqrtBetaMean[iVox]*DosePerPrimary
									      ) / ( 2. *nBeta*sqrt( NegLnS_ion/nBeta + pow(nAlpha/2./nBeta,2) ) );
							double GradRbe;
							if (fDoseGrid[iRoiId].Dose[idx]>0.) {
								GradRbe = ( GradBioDose - fDoseGrid[iRoiId].Rbe[idx] * DosePerPrimary )
									/ fDoseGrid[iRoiId].Dose[idx] * fDoseGrid[iRoiId].Rbe[idx];//the whole line is: abs dose
							} else {
								GradRbe = 0.; //assume zero for the zero-dose limiting case the initial RBE is given by RBE_alpha=alpha_ion/alpha_photon
							}
							//dose part
							GradChi2[iPb]+=-2.*(Weight //weight
									    * (ROI.box[iRoiId].TargetDose
									       - fDoseGrid[iRoiId].Dose[idx]) //(target dose - current dose)
									    / pow(ROI.box[iRoiId].TargetDose,2) // div by (target dose)^2
									    *( DosePerPrimary//absorbed dose per primary (the grad(dose)*RBE)
									       * ((isBiologicalOptimization)?fDoseGrid[iRoiId].Rbe[idx] :1.) //if bio opt: multiply by current RBE
										//Here one could add the Grad(RBE) term: Grad(BioDose)=RBE*Grad(Dose) + Dose*Grad(RBE)
							 //but don't since according to A. Gemmel et al. 2008 not necessary
										+ fDoseGrid[iRoiId].Dose[idx] * GradBioDose ) //in case
									     );
							//RBE part
							GradChi2[iPb]+=-2.*( Weight //weight
									     * (RbeMean - fDoseGrid[iRoiId].Rbe[idx]) //(target RBE - current RBE)
									     / pow(RbeMean,2) // div by (target RBE)^2
									     * GradRbe //grad(RBE)
									     );

						} else if ((ROI.box[iRoiId].Id<0)&&
							 (fDoseGrid[iRoiId].Dose[idx]>ROI.box[iRoiId].TargetDose)) {//heavy function
							GradChi2[iPb]+=-2.*(Weight //weight
									    * (ROI.box[iRoiId].TargetDose
									       - fDoseGrid[iRoiId].Dose[idx]) //(target dose - current dose)
									    / pow(ROI.box[iRoiId].TargetDose,2) // div by (target dose)^2
									    * DosePerPrimary //absorbed dose per primary
									    * ((isBiologicalOptimization)?fDoseGrid[iRoiId].Rbe[idx] :1.) //if bio opt: multiply by current RBE
									    //Here one could add the Grad(RBE) term: Grad(BioDose)=RBE*Grad(Dose) + Dose*Grad(RBE)
									    //but don't since according to A. Gemmel et al. 2008 not necessary
									    //                               * GradBioDose//in case
									    );
						} else {
							cerr<<"ERROR: Optimizer::Optimizer::ComputeStepPGRBE(): Entry: "<<iVox<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
							exit(-1);	// FIXME
						}
					}

				}//end of voxel loop
			}//end of RoiBoxId

		}//EO PB  active
	}// EO PB loop

	//TODO: Our step is not correctly calculated yet!!!
	//gradient is calculated - now compute the good stepsize
	//zero search algorithm

	// intialize the step value
	double StepSize(1.);

	//CALCULATE HERE HESSE MATRIX in diagonal approximation AND DO Newton's root finding algorithm

	//assign new intensities
	for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
		if (PB.box[iPb].Active) {
			//nScaleOptimizationStep: according to A.Gemmel for ~3Gy (RBE) -> 0.5 is good for <1Gy (RBE) use 0.25 *nScaleOptimizationStep
			PB.box[iPb].Intensity.push_back(PB.box[iPb].Intensity.back()-GradChi2[iPb]*StepSize);
			//don't allow negative particle numbers
			if (PB.box[iPb].Intensity.back()<0.) {
				PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=0.;
			}
		} else {//PB  not active
			PB.box[iPb].Intensity.push_back(0.);
		}
	}// EO PB loop

	//compute Chi2:
	double InitialChi2=fCostFunction.back();
	double CurrentChi2=ComputeCostDoseRbe();
	double PreviousChi2=InitialChi2;
	double MinChi2=InitialChi2;

	double MinChi2Step=0.;

	int Counter=0;
	double CurrentChange=.9;
	if (Verbose) cout<<"DEBUG: START STEP SEARCH: "<<Counter<<" Curr Chi2:"<<CurrentChi2<<" Prev Chi2:"<<PreviousChi2<<" Initial Chi2:"<<InitialChi2<<" Stepsize:"<<StepSize<<endl;

	bool FlagNegParticles=false;

	while((fabs(CurrentChi2-PreviousChi2)/PreviousChi2>nConvergenceCriterium*1.)||(Counter<20)) {

		if (CurrentChi2>PreviousChi2) CurrentChange=-0.5*CurrentChange;

		StepSize*=1.+CurrentChange;
		cout<<"DEBUG: Counter: "<<Counter<<" Step size: "<<StepSize<<" Chi2:"<<CurrentChi2<<endl;

		FlagNegParticles=false;
		for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
			if (PB.box[iPb].Active) {
				PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-2]-GradChi2[iPb]*(StepSize);
				//don't allow negative particle numbers
				if (PB.box[iPb].Intensity.back()<0.) {
					PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=0.;
				}
			}
		}// EO PB loop

		//TODO: add a small change for step random scaling: (1-[0-1))/(1-[0-1)) test if better chi2, if not take the one from before

		PreviousChi2=CurrentChi2;
		ComputeBioDose();
		CurrentChi2=ComputeCostDoseRbe();

		if (MinChi2>CurrentChi2) {//new best step found
			if (!FlagNegParticles) {
				MinChi2=CurrentChi2;
				MinChi2Step=StepSize;
			} else {
				cout<<"WARNING: New Chi2 minimum: "<<MinChi2<<" could not be set as it contained negative particle numbers!"<<endl;
			}
		}

		Counter++;
		if (Verbose) cout<<"  DEBUG: "<<Counter<<" Curr Chi2:"<<CurrentChi2<<" Prev Chi2:"<<PreviousChi2<<" Initial Chi2:"<<InitialChi2<<" Stepsize:"<<StepSize<<" Cur Change:"<<CurrentChange<<endl;
	}
	if (Verbose) cout<<"DEBUG: STOP STEP SEARCH"<<endl;

	if (MinChi2Step==0.) {
		cout<<"WARNING: MinChi2Step==0: endless loop!"<<endl;
	}

	//now set the min step-size we found
	for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
		if (PB.box[iPb].Active) {
			PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-2]-GradChi2[iPb]*(MinChi2Step);
			//don't allow negative particle numbers
			if (PB.box[iPb].Intensity.back()<0.) {
				PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=0.;
			}
		}
	}// EO PB loop
	return 0;
} // ComputeStepPgRbe2

/*
 * Compute with Plain Gradient method a flat RBE
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeStepPgRbeOnly()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeStepPgRbeOnly"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	double RbeMean=nTargetRbeMean;
	double GradChi2[PB.box.size()];

	const double nBeta=CellLine.nBeta();
	const double nAlpha=CellLine.nAlpha();

	//compute the current gradient
	for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
		GradChi2[iPb]=0.;

		if (PB.box[iPb].Active) {
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {//for each voxel receiving dose by current PB
					double Weight = ROI.box[iRoiId].Weight[iVox] ;
					float DosePerPrimary=Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];
					int idx =Matrix.box[iRoiId][iPb].VoxelIndex[iVox];
					if ((fDoseGrid[iRoiId].RoiId[idx]>=0)//voxel is a PTV or OAR voxel
							&&(fDoseGrid[iRoiId].Dose[idx]>0.)//current voxel dose>0
							&&(ROI.box[iRoiId].TargetDose>0.)) {//target voxel dose>0
						//now compute gradient of Chi2
						if (ROI.box[iRoiId].Id>0) {//a PTV
							//compute grad(biodose) and grad(RBE) for a given VOXEL and PB
							double NegLnS_ion= fDoseGrid[iRoiId].AlphaMeanDose[idx]
								+pow(fDoseGrid[iRoiId].SqrtBetaMeanDose[idx],2);
							double GradBioDose= ( Matrix.box[iRoiId][iPb].AlphaMean[iVox]*DosePerPrimary
									      + 2.*fDoseGrid[iRoiId].SqrtBetaMeanDose[idx]
									      * Matrix.box[iRoiId][iPb].SqrtBetaMean[iVox]*DosePerPrimary
									      ) / ( 2. *nBeta*sqrt( NegLnS_ion/nBeta + pow(nAlpha/2./nBeta,2) ) );
							double GradRbe;
							if (fDoseGrid[iRoiId].Dose[idx]>0.) {
								GradRbe = ( GradBioDose - fDoseGrid[iRoiId].Rbe[idx] * DosePerPrimary )
									/ fDoseGrid[iRoiId].Dose[idx] * fDoseGrid[iRoiId].Rbe[idx];//the whole line is: abs dose
							} else {
								GradRbe = 0.; //assume zero for the zero-dose limiting case the initial RBE is given by RBE_alpha=alpha_ion/alpha_photon
							}

							//RBE part
							GradChi2[iPb]+=-2.*( Weight//weight
									     * (RbeMean - fDoseGrid[iRoiId].Rbe[idx]) //(target RBE - current RBE)
									     / pow(RbeMean,2) // div by (target RBE)^2
									     * GradRbe //grad(RBE)
									     );

						} else if ((ROI.box[iRoiId].Id<0)&&
							 (fDoseGrid[iRoiId].Dose[idx]>ROI.box[iRoiId].TargetDose)) {//heaviside function
						} else {
							//Do nothing
							cerr<<"ERROR: Optimizer::Optimizer::ComputeStepPGRBE(): Entry: "<<iVox<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
							exit(-1);	// FIXME
						}
					}
				}//end of voxel loop
			}//end of RoiBoxId
		}//EO PB  active
	}// EO PB loop

	//TODO: Our step is not correctly calculated yet!!!
	//gradient is calculated - now compute the good stepsize
	//zero search algorithm

	// intialize the step value
	double StepSize(100.);

	//assign new intensities
	for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
		if (PB.box[iPb].Active) {
			//nScaleOptimizationStep: according to A.Gemmel for ~3Gy (RBE) -> 0.5 is good for <1Gy (RBE) use 0.25 *nScaleOptimizationStep
			PB.box[iPb].Intensity.push_back(PB.box[iPb].Intensity.back()-GradChi2[iPb]*StepSize);

			//don't allow negative particle numbers
			if (PB.box[iPb].Intensity.back()<0.) {
				PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=0.;
			}
		} else {//PB  not active
			PB.box[iPb].Intensity.push_back(0.);
		}
	}// EO PB loop

	//compute Chi2:
	double InitialChi2=fCostFunction.back();
	double CurrentChi2=ComputeCostRbe();
	double PreviousChi2=InitialChi2;
	double MinChi2=InitialChi2;

	double MinChi2Step=0.;

	int Counter=0;
	double CurrentChange=.9;
	if (Verbose) cout<<"DEBUG: START STEP SEARCH: "<<Counter<<" Curr Chi2:"<<CurrentChi2<<" Prev Chi2:"<<PreviousChi2<<" Initial Chi2:"<<InitialChi2<<" Stepsize:"<<StepSize<<endl;

	bool FlagNegParticles=false;

	while((fabs(CurrentChi2-PreviousChi2)/PreviousChi2>nConvergenceCriterium*1.)||(Counter<9)) {
		if (CurrentChi2>PreviousChi2) CurrentChange=-0.5*CurrentChange;
		//    PreviousStepSize=StepSize;
		StepSize*=1.+CurrentChange;

		FlagNegParticles=false;
		for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
			if (PB.box[iPb].Active) {
				PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-2]-GradChi2[iPb]*(StepSize);
				//don't allow negative particle numbers
				if (PB.box[iPb].Intensity.back()<0.) {
					PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=0.;
				}
			}
		}// EO PB loop

		//TODO: add a small change for step random scaling: (1-[0-1))/(1-[0-1)) test if better chi2, if not take the one from before

		PreviousChi2=CurrentChi2;
		ComputeBioDose();
		CurrentChi2=ComputeCostRbe();

		if (MinChi2>CurrentChi2) {//new best step found
			if (!FlagNegParticles) {
				MinChi2=CurrentChi2;
				MinChi2Step=StepSize;
			} else {
				cout<<"WARNING: New Chi2 minimum: "<<MinChi2<<" could not be set as it contained negative particle numbers!"<<endl;
			}
		}

		Counter++;
		if (Verbose) cout<<"  DEBUG: "<<Counter<<" Curr Chi2:"<<CurrentChi2<<" Prev Chi2:"<<PreviousChi2<<" Initial Chi2:"<<InitialChi2<<" Stepsize:"<<StepSize<<" Cur Change:"<<CurrentChange<<endl;
	}
	if (Verbose) cout<<"DEBUG: STOP STEP SEARCH"<<endl;

	if (MinChi2Step==0.) {
		cout<<"WARNING: MinChi2Step==0: endless loop!"<<endl;
	}

	//now set the min step-size we found
	for (unsigned int iPb=PB.box.size();iPb--;) {//for each PB
		if (PB.box[iPb].Active) {
			PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-2]-GradChi2[iPb]*(MinChi2Step);
			//don't allow negative particle numbers
			if (PB.box[iPb].Intensity.back()<0.) {
				PB.box[iPb].Intensity[PB.box[iPb].Intensity.size()-1]=0.;
			}
		}
	}// EO PB loop
	return 0;
} // ComputeStepPgRbeOnly

/*
 * TEST Compute step with dose-difference (DD) + LET-difference scaling
 * Dose-difference + LET-difference is based on ComputeStepDRbeD2() method!!!
 */
/*--------------------------------------------------------*/
int Optimizer::ComputeStepDLetD()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::ComputeStepDLetD"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	double LetMean=nTargetLetMean;
	//TODO implement floating constant LET
	//LetMean=1.8;cout<<"WARNING: Fixed target mean LET of:"<<LetMean<<endl;

	double LetRelImportance=1.;

	const double nBeta=CellLine.nBeta();
	const double nAlphaDivBeta=CellLine.nAlphaDivBeta();

	for (unsigned int iPb=PB.box.size();iPb--;) {
		if (PB.box[iPb].Active) {
			double OtherUp(0.),OtherDown(0.);
			for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
				for (unsigned int iVox=Matrix.box[iRoiId][iPb].VoxelHitAboveThreshold;iVox--;) {
					double Weight = ROI.box[iRoiId].Weight[iVox];//weight
					float DosePerPrimary = Matrix.box[iRoiId][iPb].DosePerPrimary[iVox];//absorbed dose per primary
					float TargetDose=ROI.box[iRoiId].TargetDose;
					double Dose=fDoseGrid[iRoiId].Dose[Matrix.box[iRoiId][iPb].VoxelIndex[iVox]]; //(target dose - current dose)
					float Let=fDoseGrid[iRoiId].Let[Matrix.box[iRoiId][iPb].VoxelIndex[iVox]];

					if ((fDoseGrid[iRoiId].RoiId[Matrix.box[iRoiId][iPb].VoxelIndex[iVox]]>=0)//voxel is a PTV or OAR voxel
							&&(Dose>0.)//current voxel dose>0
							&&(TargetDose>0.)//target voxel dose>0
							&&(Let>0.)) {//LET in voxel >0

						//compute absorbed dose in current voxel
						double CurrentDose=Dose/max(Let,(float)1.E-20);
						//compute LET for current PB assuming that all current dose is only by this beam
						double NegLnS_ion= Matrix.box[iRoiId][iPb].AlphaMean[iVox]*CurrentDose
								+ pow(Matrix.box[iRoiId][iPb].SqrtBetaMean[iVox]*CurrentDose,2);
						double LetPb=(sqrt( NegLnS_ion/nBeta+pow(0.5*nAlphaDivBeta,2) ) - 0.5*nAlphaDivBeta)/CurrentDose;

						if (ROI.box[iRoiId].Id>0) {//a PTV
							OtherUp+=Weight*(
									 (LetMean-Let) //Let part
									 * (LetPb-LetMean)*DosePerPrimary
									 *(LetRelImportance*fCurrentChi2Let/(fCurrentChi2Dose+LetRelImportance*fCurrentChi2Let))
									 +
									 (fCurrentChi2Dose/(fCurrentChi2Dose+LetRelImportance*fCurrentChi2Let))*
									 DosePerPrimary*LetPb*
									 (TargetDose//Bio dose part
									  - Dose)
									 )*DosePerPrimary*LetPb*PB.box[iPb].Intensity.back()//damping factor
								/Dose;
							OtherDown+=Weight *(
									    pow(DosePerPrimary*LetPb,2) //dose part
									    *(fCurrentChi2Dose/(fCurrentChi2Dose+LetRelImportance*fCurrentChi2Let))
									    +pow((LetPb-LetMean)*DosePerPrimary,2) //LET part
									    *(LetRelImportance*fCurrentChi2Let/(fCurrentChi2Dose+LetRelImportance*fCurrentChi2Let))
									    );
						} else if ((ROI.box[iRoiId].Id<0)&&(Dose>TargetDose)) {//heavy func
							OtherUp+=Weight*(
									 (fCurrentChi2Dose/(fCurrentChi2Dose+LetRelImportance*fCurrentChi2Let))*
									 DosePerPrimary*LetPb* //Bio dose part
									 (TargetDose - Dose)
									 )*DosePerPrimary*LetPb*PB.box[iPb].Intensity.back()//damping factor
								/Dose;
							OtherDown+=Weight*(
									   (fCurrentChi2Dose/(fCurrentChi2Dose+LetRelImportance*fCurrentChi2Let))*
									   pow(DosePerPrimary*LetPb,2) //bio dose part
									   );
						} else {
							cerr<<"ERROR: Optimizer::Optimizer::ComputeStepDLetD(): Entry: "<<iVox<<" PTV/OAR value not allowed:"<<ROI.box[iRoiId].Id<<"!"<<endl;
							exit(-1);	// FIXME
						}
					}
				}//end of voxel loop
			}//end of RoiBoxId

			if (OtherDown>0.) {
				PB.box[iPb].Intensity.push_back( PB.box[iPb].Intensity.back()+nScaleOptimizationStep*OtherUp/OtherDown );
			} else {
				cout<<"WARNING: Optimizer::ComputeStepDLetD(): PB: "<<iPb<<" OtherUp: "<<OtherUp<<" OtherDown: "<<OtherDown<<", keep old PB intensity!"<<endl;
				PB.box[iPb].Intensity.push_back(PB.box[iPb].Intensity.back());
			}
		} else {//PB  not active
			PB.box[iPb].Intensity.push_back(0.);
		}
	}
	return 0;
} // ComputeStepDLetD

/*
 * Evaluate the dose distribution
 */
/*--------------------------------------------------------*/
int Optimizer::Evaluate()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::Evaluate"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::Evaluate()!"<<endl;
		exit(-1);	// FIXME
	}

	//store stuff
	SavePb();
	//TODO:  SaveDose();

	//print stuff
	PrintEvaluation();
	PlotDvh();
	return 0;
} // Evaluate

/*--------------------------------------------------------*/
void Optimizer::SavePb()
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Optimizer::SavePb"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	ostringstream comm;
	if (fCurrentIteration==-1) {
		PB.writePB(string(fDataOut+"/pb_final.res"),&PB.box); //TODO this is done in very naive way
	} else if (fCurrentIteration==-2) {//selective dose
		//NOT WORKING, just does the same WritePencilBeams(char*::Format("%s/pb_select.res",fDataOut.Data()));
	} else {
		comm.str("");
		comm.clear();
		comm<<fDataOut<<"/pb_"<<fCurrentIteration<<".res";
		PB.writePB((comm.str()).c_str(),&PB.box);
		comm.str("");
		comm.clear();
		comm<<fDataOut<<"/pb_"<<fCurrentIteration<<"_id";
	}
} // SavePb

/*--------------------------------------------------------*/
void Optimizer::PrintEvaluation()
{
	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::PrintEvaluation()"<<endl;
		exit(-1);	// FIXME
	}

	//compute min,max,average deviation in each PTV and OAR in each ROI
	vector<int>    Type,Entries;
	vector<double> DevAvg,DevMin,DevMax;

	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=ROI.box[iRoiId].VoxelIndex.size();i0--;) {
			bool IsRegistered=false;
			int idx = ROI.box[iRoiId].VoxelIndex[i0];
			for (unsigned int i1=Type.size();i1--;) {
				if (Type[i1]==ROI.box[iRoiId].Id) {
					IsRegistered=true;
					Entries[i1]++;
					DevAvg[i1]+=fabs(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose)/ROI.box[iRoiId].TargetDose;
					if (DevMax[i1]<fDoseGrid[iRoiId].Dose[idx]/ROI.box[iRoiId].TargetDose) DevMax[i1]=fDoseGrid[iRoiId].Dose[idx]/ROI.box[iRoiId].TargetDose;
					if (DevMin[i1]>fDoseGrid[iRoiId].Dose[idx]/ROI.box[iRoiId].TargetDose) DevMin[i1]=fDoseGrid[iRoiId].Dose[idx]/ROI.box[iRoiId].TargetDose;
					break;
				}
			}
			if (!IsRegistered) {
				Type.push_back(ROI.box[iRoiId].Id);
				Entries.push_back(1);
				DevAvg.push_back( fabs(fDoseGrid[iRoiId].Dose[idx]-ROI.box[iRoiId].TargetDose)/ROI.box[iRoiId].TargetDose );
				DevMax.push_back( fDoseGrid[iRoiId].Dose[idx]/ROI.box[iRoiId].TargetDose );
				DevMin.push_back( fDoseGrid[iRoiId].Dose[idx]/ROI.box[iRoiId].TargetDose );
			}
		}
	}

	for (unsigned int i1=Type.size();i1--;) {
		DevAvg[i1]/=max((double)Entries[i1],1.E-10);
	}

	//print
	cout<<"### Evaluation ####################################"<<endl;
	cout<<"  Evaluation for current dose distribution:"<<endl;
	for (unsigned int i1=0;i1<Type.size();i1++) {
		if (Type[i1]>0) {
			cout<<"  PTV: "<<setw(3)<<Type[i1]<<" Deviation Avg: "<<setw(6)<<DevAvg[i1]*100.
					<<"% Min: "<<setw(6)<<(DevMin[i1]-1.)*100.
					<<"% Max: "<<setw(6)<<(DevMax[i1]-1.)*100.<<"%"<<endl;
		}
	}
	for (unsigned int i1=0;i1<Type.size();i1++) {
		if (Type[i1]<0) {
			cout<<"  OAR: "<<setw(3)<<Type[i1]<<" Deviation Avg: "<<setw(6)<<DevAvg[i1]*100.
					<<"% Min: "<<setw(6)<<(DevMin[i1]-1.)*100.
					<<"% Max: "<<setw(6)<<(DevMax[i1]-1.)*100.<<"%"<<endl;
		}
	}

	cout<<"###################################################"<<endl;
} // PrintEvaluation

void Optimizer::PlotDvh()
{
	ostringstream comm;

	double MaxDosePercent=200.;
	int DoseBinNumber=(int)(MaxDosePercent*10.);

	//do differential and cumulative DVHs for multiple PTVs and OARs

	vector<int>    Type,Voxel;
	vector<string>  Tag;
	vector<double> TargetDose;
	vector<H1D*>    DvhD,DvhI;
	for (unsigned int iRoiId=ROI.box.size();iRoiId--;) {
		for (unsigned int i0=ROI.box[iRoiId].VoxelIndex.size();i0--;) {
			bool IsRegistered=false;
			int idx = ROI.box[iRoiId].VoxelIndex[i0];
			for (unsigned int i1=Type.size();i1--;) {
				if (Type[i1]==ROI.box[iRoiId].Id) {
					IsRegistered=true;
					Voxel[i1]++;
					DvhD[i1]->fill(fDoseGrid[iRoiId].Dose[idx]);//absolute abs./bio. dose
				//	cout <<DvhD[i1]->Fill(fDoseGrid[iRoiId].Dose[idx])<<"     ";
					break;
				}
			}
			if (!IsRegistered) {
				Type.push_back(ROI.box[iRoiId].Id);
				Tag.push_back(ROI.box[iRoiId].Tag);
				Voxel.push_back(1);
				comm.str("");
				comm.clear();
				comm<<"DvhDifferential"<<Type.back();
				DvhD.push_back(new H1D((comm.str()).c_str(),DoseBinNumber,0.,MaxDosePercent/100.*ROI.box[0].TargetDose));//assume here that first ROI given is PTV (i.e.: iRoiId=0 and i0=0)

				comm.str("");
				comm.clear();
				comm<<"DvhCumulative"<<Type.back();
				DvhI.push_back(new H1D((comm.str()).c_str(),DoseBinNumber,0.,MaxDosePercent/100.*ROI.box[0].TargetDose));//assume here that first ROI given is PTV (i.e.: iRoiId=0 and i0=0)

				DvhD.back()->fill(fDoseGrid[iRoiId].Dose[idx]);
//				DvhI.push_back(new H1D(DvhD[i1]));
				TargetDose.push_back(ROI.box[iRoiId].TargetDose);

				//Add here EUD computation
				//Compute here mean survival in PTV/OAR
				//				MeanSurvival+=XXX;
				//MeanSurvivalWeight+=1.;
			}
		}
	}

	for (unsigned int i1=Type.size();i1--;) {//loop over PTVs/OARs
		DvhD[i1]->scale(100./DvhD[i1]->entries());//normalize
		//compute cumulative DVH
		for (int i2=DoseBinNumber-1;i2>=0;i2--) {
			for (int  i3=0;i3<=i2;i3++) {
				DvhI[i1]->set(i3,DvhI[i1]->get(i3)+DvhD[i1]->get(i2));
			}
		}

		double Dmean,D1(-1.),D99(-1.),D5(-1.),D95(-1.);
		double D1rel(-1.),D99rel(-1.),D5rel(-1.),D95rel(-1.);

		cout<<"GetMean "<<DvhD[i1]->mean();
		Dmean=DvhD[i1]->mean()/max(TargetDose[i1],1.E-10)*100.;//from differential DVH

		for (int  i2=0;i2<DoseBinNumber;i2++) {
			//compute: near max/min and  inhomogeneity coefficients (IC1 and IC5) IC1=(D1-D99)/Dmean

			if (DvhI[i1]->get(i2)>99.) D99=(DvhI[i1]->center(i2)+DvhI[i1]->width()/2.) /max(TargetDose[i1],1.E-10)*100.;//store right edge of bin
			if (DvhI[i1]->get(i2)>95.) D95=(DvhI[i1]->center(i2)+DvhI[i1]->width()/2.) /max(TargetDose[i1],1.E-10)*100.;//store right edge
			if (DvhI[i1]->get(i2)>05.) D5 =(DvhI[i1]->center(i2)+DvhI[i1]->width()/2.) /max(TargetDose[i1],1.E-10)*100.;//store right edge
			if (DvhI[i1]->get(i2)>01.) D1 =(DvhI[i1]->center(i2)+DvhI[i1]->width()/2.) /max(TargetDose[i1],1.E-10)*100.;//store right edge
		}

		//compute: RELATIVE (to mean dose) near max/min and RELATIVE inhomogeneity coefficients (IC1 and IC5) IC1=(D1-D99)/Dmean
		D99rel=D99/Dmean*100.;
		D95rel=D95/Dmean*100.;
		D5rel=D5/Dmean*100.;
		D1rel=D1/Dmean*100.;

		cout<<endl;cout<<endl;
		cout<<"Quantities for PTV/OAR number: "<<i1<<endl;
		cout<<"  Quantities relative to TARGET mean PTV/OAR dose ("<<TargetDose[i1]<<"Gy (RBE))"<<endl;
		cout<<"  PTV/OAR:"<<setw(4)<<Type[i1]<<" Dose_mean = "<<Dmean<<"% Dose_meandev = "<<Dmean-100.<<"%"<<endl;
		cout<<"  Near maximum doses: D1  = "<<D1-100. <<"% D5  = "<<D5-100.<<"%"<<endl;
		cout<<"  Near minimum doses: D99 = "<<D99-100.<<"% D95 = "<<D95-100.<<"%"<<endl;
		cout<<"  Inhomogeneity coefficient IC1=(D1-D99)/Dmean = "<<(D1-D99)/Dmean*100.<<"%"<<endl;
		cout<<"  Inhomogeneity coefficient IC5=(D5-D95)/Dmean = "<<(D5-D95)/Dmean*100.<<"%"<<endl;
		cout<<endl;
		cout<<"  Quantities relative to ACTUAL mean PTV/OAR dose ("<<Dmean/100.*TargetDose[i1]<<"Gy (RBE))"<<endl;
		cout<<"  Rel. near maximum doses: D1  = "<<D1rel-100. <<"% D5  = "<<D5rel-100.<<"%"<<endl;
		cout<<"  Rel. near minimum doses: D99 = "<<D99rel-100.<<"% D95 = "<<D95rel-100.<<"%"<<endl;
		cout<<"  Rel. inhomogeneity coefficient IC1=(D1-D99)/Dmean = "<<(D1-D99)/Dmean*100.<<"%"<<endl;
		cout<<"  Rel. inhomogeneity coefficient IC5=(D5-D95)/Dmean = "<<(D5-D95)/Dmean*100.<<"%"<<endl;
	}

	//save DVHs to file
	if (fCurrentIteration==-1) {
		comm.str("");
		comm.clear();
		comm<<fDataOut<<"/Dvh_final.res";
		WriteDvh((comm.str()),&Tag,&Type,&DvhD,&DvhI);
	} else if (fCurrentIteration==-2) {
		comm.str("");
		comm.clear();
		comm<<fDataOut<<"/Dvh_select.res";
		WriteDvh((comm.str()),&Tag,&Type,&DvhD,&DvhI);
	} else {
		comm.str("");
		comm.clear();
		comm<<fDataOut<<"/Dvh_"<<fCurrentIteration<<".res";
		WriteDvh((comm.str()),&Tag,&Type,&DvhD,&DvhI);
	}

	// clean up
	for (unsigned i=0; i<DvhD.size(); i++) delete DvhD[i];
	for (unsigned i=0; i<DvhI.size(); i++) delete DvhI[i];
} // PlotDvh

/*--------------------------------------------------------*/
void Optimizer::PrintInitial()
{
	// cout<<"In PrintInitial function"<<vPbFiles[0].PB.BoxSet[0].Intensity.back()<<endl;
	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::PrintInitial()!"<<endl;
		exit(-1);	// FIXME
	}

	cout<<"### Initialized ###################################"<<endl;
	cout<<"  Optimizer was initialized with configuration:"<<endl;
	cout<<"  ROI boxes read:            "<<ROI.box.size()<<endl;
	cout<<"  ROI tag(s):  ";
	for (unsigned int iRoiId=0;iRoiId<ROI.box.size();iRoiId++) {
		cout<<"<"<<ROI.box[iRoiId].Tag<<"> ";
		cout<<"<"<<ROI.box[iRoiId].Id<<"> ";
	}
	cout<<endl;
	cout<<"  LET optimization:  "<<(isLetOptimization==0?"False":"True")<<endl;
	cout<<"  Biological optimization:   "<<(isBiologicalOptimization==0?"False":"True")<<endl;
	if (isBiologicalOptimization) {
	  cout<<"  Assuming for X-ray parameters cell line <"<<CellLine.cellName()<<">"<<endl;
	}
	cout<<"  RBE is fixed:              "<<(isRbeFixed==0?"False":"True");
	if (isRbeFixed) { cout <<" (RBE = "<<nRbeFixed<<") "<<endl;}
	else{cout<<endl;}
	cout<<"  RBE exponential scaling factor: "<<nRbeFactor<<endl;
	cout<<"  Optimization algorithm:    ";
	switch(nOptimizationAlgorithm) {
	case 2 : cout<<"<Dose-difference + RBE-difference>"<<endl;
	if (!isBiologicalOptimization) {
		cerr<<"ERROR: Optimizer::PrintInitial(): for this algorithm bio optimization has to be switched on!"<<endl;
		exit(-1);	// FIXME
	}
	break;
	case 3 : cout<<"<Plain Gradient>"<<endl;
	break;
	case 4 : cout<<"<Plain Gradient for dose and RBE>"<<endl;
	if (!isBiologicalOptimization) {
		cerr<<"ERROR: Optimizer::PrintInitial(): for this algorithm bio optimization has to be switched on!"<<endl;
		exit(-1);	// FIXME
	}
	break;
	case 5 : cout<<"<Dose-difference + RBE-difference 2>"<<endl;
	if (!isBiologicalOptimization) {
		cerr<<"ERROR: Optimizer::PrintInitial(): for this algorithm bio optimization has to be switched on!"<<endl;
		exit(-1);	// FIXME
	}
	break;
	case 6 : cout<<"<Plain Gradient for dose and RBE 2>"<<endl;
	if (!isBiologicalOptimization) {
		cerr<<"ERROR: Optimizer::PrintInitial(): for this algorithm bio optimization has to be switched on!"<<endl;
		exit(-1);	// FIXME
	}
	break;

	case 7 : cout<<"<Plain Gradient for RBE only>"<<endl;
	if (!isBiologicalOptimization) {
		cerr<<"ERROR: Optimizer::PrintInitial(): for this algorithm bio optimization has to be switched on!"<<endl;
		exit(-1);	// FIXME
	}
	break;

	case 8 : cout<<"<Plain Gradient hybrid switching Dose/RBE only>"<<endl;
	if (!isBiologicalOptimization) {
		cerr<<"ERROR: Optimizer::PrintInitial(): for this algorithm bio optimization has to be switched on!"<<endl;
		exit(-1);	// FIXME
	}
	break;
	case 9 : cout<<"TODO: Add LET OPTIMIZATION!"<<endl;//TODO
		cout<<"<Dose-difference + LET-difference 2>"<<endl;
		exit(-1);	// FIXME
	if (!isLetOptimization) {
		cerr<<"ERROR: Optimizer::PrintInitial(): for this algorithm LET optimization has to be switched on!"<<endl;
		exit(-1);	// FIXME
	}
	default: cout<<"<Dose-difference>"<<endl;//this is case 1;
	}
	cout<<"  Convergence criterium:     "<<nConvergenceCriterium*100<<"%"<<endl;
	cout<<"  Threshold min nb. of particles per PB: "<<nMinParticlesPerPb<<endl;
	cout<<"  Max iterations:            "<<nIterationMax<<endl;
	cout<<"  Iteration action interval: "<<nIterationActionInterval<<" step(s)"<<endl;
	cout<<"  Scaling of optim. step:    "<<nScaleOptimizationStep<<endl;
	cout<<"  Data input:                "<<fDataIn<<endl;
	cout<<"  Data output:               "<<fDataOut<<endl;
	cout<<"  Time for initialization:   "<<((double) clock())/CLOCKS_PER_SEC/60.<<" min"<<endl;
	cout<<"###################################################"<<endl;
} // PrintInitial

/*--------------------------------------------------------*/
void Optimizer::PrintFinal()
{
	// cout<<"In PrintFinal  function"<<vPbFiles[0].PB.BoxSet[0].Intensity.back()<<endl;
	if (!isInitialized) {
		cerr<<"ERROR: Initialize before calling Optimizer::PrintFinal()!"<<endl;
		exit(-1);	// FIXME
	}
	Evaluate();
} // PrintFinal

/*--------------------------------------------------------*/
void Optimizer::PrintVoxelIndex(int Index, int RoiBoxId)
{
	// cout<<"In PrintVoxelIndex  function"<<vPbFiles[0].PB.BoxSet[0].Intensity.back()<<endl;
	Vector Pos=	ROI.index2Position(Index);
	VECTOR_INT Bin=	ROI.index2Bin(Index);
	cout<<" RoiBoxId: "<<RoiBoxId<<" "<<" VoxelIndex: "<<Index<<" "<<ROI.position2Index(Pos,RoiBoxId)
	    <<" Position: "<<Pos.x<<" "<<Pos.y<<" "<<Pos.z<<" Bins: "<<Bin.V[0]<<" "<<Bin.V[1]<<" "<<Bin.V[2]<<endl;
} // PrintVoxelIndex

/*--------------------------------------------------------*/
void Optimizer::PrintLoadedPb()
{
	// cout<<"In PrintLoadedPb  function"<<vPbFiles[0].PB.BoxSet[0].Intensity.back()<<endl;
	cout<<"### Loaded PBs ####################################"<<endl;
	cout<<"Total: "<<PB.box.size()<<endl;
	for (unsigned int i0=0;i0<PB.box.size();i0++) {
		//save the pencil beams with the optimized intensities
		cout<<setw(5)<<i0
		    <<" Pos:"<<setw(8)<<PB.box[i0].Position.x<<" "<<setw(8)<<PB.box[i0].Position.y
		    <<" "<<setw(8)<<PB.box[i0].Position.z
		    <<" Mom:"<<setw(8)<<PB.box[i0].Energy
		    <<" Z:"<<setw(3)<<PB.box[i0].Z<<" A:"<<setw(3)<<PB.box[i0].A
		    <<" Spot:"<<setw(5)<<PB.box[i0].SpotSigma[0]<<setw(5)<<PB.box[i0].SpotSigma[1]
		    <<" Intens:"<<setw(8)<<PB.box[i0].Intensity.back()
		<<" Act:"<<PB.box[i0].Active<<endl;
	}
	cout<<"###################################################"<<endl;
} // PrintLoadedPb

/*--------------------------------------------------------*/
void Optimizer::PrintLoadedPbMatrix()
{
	//cout<<"In PrintLoadedPbMatrix  function"<<vPbFiles[0].PB.BoxSet[0].Intensity.back()<<endl;
	cout<<"### Loaded PB matrices ############################"<<endl;

	for (unsigned int iRoiId=0;iRoiId<ROI.box.size();iRoiId++) {
		cout<<"ROI BOX ID: "<<iRoiId <<" Matrix size: "<<Matrix.box[iRoiId].size()<<endl;
		for (unsigned int i0=0;i0<Matrix.box[iRoiId].size();i0++) {

			cout<<setw(5)<<i0<<"  PbIndex:"<<setw(5)<<Matrix.box[iRoiId][i0].PbIndex
					<<" VoxelHit:"<<setw(8)<<Matrix.box[iRoiId][i0].VoxelHit<<" VoxelHit>Thres#:"<<setw(8)<<Matrix.box[iRoiId][i0].VoxelHitAboveThreshold
			    <<" Prim:"<<setw(8)<<Matrix.box[iRoiId][i0].Events
					<<" ME(0): VoxelID:"<<setw(8)<<Matrix.box[iRoiId][i0].VoxelIndex[0]<<" D:"<<setw(8)<<Matrix.box[iRoiId][i0].DosePerPrimary[0]
				    <<" A:"<<setw(8)<<Matrix.box[iRoiId][i0].AlphaMean[0]<<" SqrtB:"<<setw(8)<<Matrix.box[iRoiId][i0].SqrtBetaMean[0]<<endl;

		}
	}
	cout<<"###################################################"<<endl;
} // PrintLoadedPbMatrix

/*--------------------------------------------------------*/
//given \sum_i^N x_i, \sum_i^N x^2 and N, this function returns in X the Average and in X2 the error
int Optimizer::CalcAvgErr(double *X, double *X2, int N)
{
	if (N<=0) {
		*X=0.;
		*X2=-1.;
		return -1;
	} else if (N==1) {
		*X2=-1.;
		return 1;
	} else {
		*X=(*X)/float(N);
		*X2=sqrt( fabs( (*X2)/float(N)- (*X)*(*X) )/float(N-1) );
		return 0;
	}
} // CalcAvgErr

/*--------------------------------------------------------*/
//given \sum_i^N x_i, \sum_i^N x^2 and N, this function returns in X the Average and in X2 the standard deviation
int Optimizer::CalcAvgStd(double *X, double *X2, int N)
{
	if (N<=0) {
		*X=0.;
		*X2=-1.;
		return -1;
	} else {
		*X=(*X)/float(N);
		*X2=sqrt( fabs( (*X2)/float(N)- (*X)*(*X) ) );
		return 0;
	}
} // CalcAvgStd

const std::string Optimizer::CurrentDateTime()
{
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	// Visit http://www.cplusplus.com/reference/clibrary/ctime/strftime/
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
} // CurrentDateTime
