//********************************************
// Author: T.Boehlen@gmail.com
// Modified: wioletta.kozlowska@cern.ch
// Created: some time 2011
// Modified: January 2016 - root extraction
//********************************************
#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <string>

#include "roi.h"
#include "pbmatrix.h"
#include "cell_line.h"
#include "pencilbeam.h"

class H1D;

class Optimizer {
public:
	bool	Verbose;

private:
	bool	isDev;                    //Turn on dev version

	std::string  fDataIn;             //Input Data localisation
	std::string  fDataOut;            //Output Data localisation
	bool   isInitialized;             //Main variable checking initialization of the variables

	int    nIterationMax;             //Max number of iterations

	//------------------------------------ Box collectors-----------------------------------------------------------//
	PencilBeam PB;
	Roi ROI;
	PBMatrix Matrix;
//	PBMatrix MatrixError;

	int    fCurrentIteration;
	double fCurrentChi2Dose;
	double fCurrentChi2Rbe;
	double fCurrentChi2Let;

	//storing and optimizing

	std::vector<double> fCostFunction;
	std::vector<double> fCostFunctionTime;//store full CPU time (also including initialization)
	std::vector<DOSE_GRID> fDoseGrid;

	Cell_Line CellLine;

	//----------------------------------------------------------------------------------------------------------------//
	//**************************** Running options. Shall be defined before Run ****************************************
	//---------------------------------- General Variables --------------------------------------------//
	bool   isBiologicalOptimization;             // Biological Optimization
	int    nOptimizationAlgorithm;               // Optimization Algorithm
	int    nCellLine;                            // Cell Line (for x-ray parameters)

	double nConvergenceCriterium;                // Convergence Criterion TODO: What does it mean
	double nScaleOptimizationStep;               // Scale the optimization step
	double nMatrixMinDoseThreshold;              // Minimum dose threshold for the dose matrix elements
	int nMinParticlesPerPb;                   // Minimum number of particles per pencil beam

	bool   isInitialEvaluation;                  // Do an initial evaluation TODO:What does it mean
	bool   isCheckPb;                            // Do an initial check of the PBs
	bool   isLetOptimization;

	//---------------------------------------- Files vectors -----------------------------------------//
	// Fluka binary dose file format  (Warning: you also need to supply the corresponding pencil beam file(s),
	// and possibly alpha, beta and LET files in the same order!)  TODO:Whatever it means
	std::string vFileDose;

	// New FLUKA binary alpha file (Warning: you also need to supply the corresponding pencil beam,
	// dose and beta file(s) in the same order!)
	std::string vFileAlpha; //TODO: Alpha File class

	// New FLUKA binary beta file  (Warning: you also need to supply the corresponding pencil beam,
	// dose and alpha file(s) in the same order!)
	std::string vFileBeta; //TODO: Beta File Class

	// New FLUKA binary LET file (Warning: you also need to supply the corresponding pencil beam,
	// dose beta file(s) and possibly alpha and beta files in the same order!)
	std::string vFileLet; //TODO: LET File Class

	// New error matrix file and its corresponding ROI box ID.
	// TODO: WATCH OUT: works only with "old format" for the moment!
	// Don't use with -mindosethres, -mirrorz and -3dfrom1d !
	// it needs the quantities dose^2  alpha^2 and beta per primary
	std::vector<std::string> vFileErrorMatrix; //TODO: Error Matrix File
	// ERROR matrix file corresponding RoiBoxId (which is starting from 0)
	std::vector<int>         vFileErrorMatrixRoiBoxId; //TODO: Error Matrix File

	std::string  FileNewIntensities;          // Read new intensities from a file TODO:structure of the file?
	bool   isReadNewIntensities;              // Read new intensities from a file yes no? TODO: what is it for?
	bool   isReadNewIntensitiesBefore;        // Read new intensities from a file before multiplying TODO: what is it for?
	double nScaleIntensities;                 // Scale intensities

	//------------------------------------ Factors and scaling values -----------------------------------------//
	double nRbeFactor;                        // This is for exponential scaling of RBE
	                                          //Set a systematic bias (scaling) of RBE (RBE^(factor))
	bool   isRbeFixed;                        // Set a fixed rbe TODO:What is it for?
	double nRbeFixed;                         // This is forcing a fixed RBE (i.e. for protons set to 1.1)
	double nTargetRbeMean;                    // Set a target mean rbe in the PTV in PTV
	double nTargetLetMean;
	bool isRbeInPtvFixed;
	bool isLetInPtvFixed;
	bool   isPhotonEquivSurvival;             // Compute photon equivalent survival
	bool   isSelectiveEvaluation;             // Do selective evaluation

	bool   isOptimize1d;                      // Do optimization only in 1D along (x=y=0)
	int    nIterationActionInterval;          // Set Action Interval TODO:What does it mean

	//----------------------------------------------Plot---------------------------------------------------//
	//*** Plotting ***
	bool   isPlot;                            // Activate plotting
	int    nPlotSlices;                       // Plot all slices
	int    nPlotRoiTest;                      // Plot all ROI slices
	bool   isPlotColour;                      // Plot in colour
	std::string  PlottingFormat;                   // Set plotting format

public :
	Optimizer() {}
	~Optimizer() {}

	//Launch the initialization and optimization
	int RunOptimizer();

	//initialization functions
	void init(std::string in, std::string out, bool IsBioOpt, int ncells, std::string vfileCell);

private:
	//Cleaning member structures
	int Clean(std::vector<DOSE_GRID> *Struct);//clear completely (to be re-initialized afterwards)
	int CleanCost();

	int ManipulateRoi(int RoiBoxId);
	int InitializeDoseGrid();
	int ResetDoseGrid();//needed after each step
	int ResetLetGrid();
	int ComputeInverseMapping();

	int CheckPencilBeams();

	//Launch the optimization
	int PreOptimize();
	int PreOptimizeDRbe();
	int Optimize1d();
	int Optimize();

	int ComputeDose();
	int ComputeBioDose();
	int ComputeAbsDose();
	int ComputeLet();

	int ComputeDoseError();
	int ComputeBioDoseError();
	int ComputeAbsDoseError();

	int ComputePhotonEquivSurvival();

	//For evaluating dose distributions produced by a sub-set of pencil-beams (e.g. for single field plotting)
	bool IsPbSelected(PENCIL_BEAM* Pb);
	int ComputeSelectiveAbsDose();
	int ComputeSelectiveBioDose();

	double ComputeRbeMean();
	double ComputeLetMean();

	//compute cost functions
	int ComputeCost();
	double ComputeCostDose();
	double ComputeCostDoseRbe();
	double ComputeCostRbe();
	double ComputeCostDoseLet();
	double ComputeCostOld();//à la FORTRAN

	//compute the next optimization step
	int ComputeStep();//generic calling of the computation of the next step
	int DiscardLowIntensityPb();// check if PB particle number falls below a given threshold, if so "zero"-it
	//single objective (absorbed OR biological dose)
	int ComputeStepD();//à la Lomax et al.
	int ComputeStepPg();//Plain Gradient algorithm for description see for instance A. Gemmel et al. 2008
	//dual objective (dose AND mean RBE in PTV)
	int ComputeStepDRbeD();//Dose-difference + RBE-difference
	int ComputeStepPgRbe();//Plain Gradient method for dose and RBE
	int ComputeStepDRbeD2();//Dose-difference + RBE-difference 2 (slightly modified)
	int ComputeStepPgRbe2();//Plain Gradient method for dose and RBE (slightly modified)
	int ComputeStepPgRbeOnly();//Compute with Plain Gradient method a flat RBE
	int ComputeStepDLetD();//Dose-difference + LET-difference based on ComputeStepDRbeD2()!!!

	/*--------------------------------------------------------*/

	//Evaluation of quantities from current iteration step
	int Evaluate();
	void SavePb();
	void PrintEvaluation();
	//ON	void PlotCost();
	void PlotDvh();

	//Printing stuff
	void PrintInitial();
	void PrintFinal();
	int PrintResults();
	void PrintVoxelIndex(int Index, int RoiBoxId);
	void PrintLoadedPb();
	void PrintLoadedPbMatrix();

	//store results
	int WriteDvh(	std::string FileName,
			std::vector<std::string> *Tag,
			std::vector<int> *Type,
			std::vector<H1D*> *DvhD,
			std::vector<H1D*> *DvhI);
	//TODO: int WriteDose(char* FileName);

	//Old stuff
	//given \sum_i^N x_i, \sum_i^N x^2 and N, this function returns in X the Average and in X2 the error
	int CalcAvgErr(double *X, double *X2, int N);

	//given \sum_i^N x_i, \sum_i^N x^2 and N, this function returns in X the Average and in X2 the standard deviation
	int CalcAvgStd(double *X, double *X2, int N);

	void PressEnter();

	const std::string CurrentDateTime();

public:
	//*---------------------------------------------------------------------------------------------------------------------*
	//***************************** Setters for initialization variables**************************************************//

	void setIsBiologicalOptimization(bool IsBiologicalOptimization) {
			isBiologicalOptimization=IsBiologicalOptimization;
		}

	void setIsLetOptimization(bool IsLetOptimization) {
			isLetOptimization=IsLetOptimization;
		}

	void setOptimizationAlgorithm(int OptimizationAlgorithm) {
			nOptimizationAlgorithm=OptimizationAlgorithm;
			switch(nOptimizationAlgorithm) {//check if LET computation is needed or not
				case 9 :
					isLetOptimization=true;
					break;
				default:
					isLetOptimization=false;
			}
		}

	void setCellLine(int n) {
			nCellLine = n;
		}

	void setConvergenceCriterium(double Criterium) {
			if (Criterium>0.)
				nConvergenceCriterium=Criterium;
			else
				std::cerr<<"ERROR: SetConvergenceCriterium(): Wrong value:"<<Criterium<<std::endl;
		}

	void setScaleOptimizationStep(double ScaleOptimizationStep) {
			if (ScaleOptimizationStep>0.)
				nScaleOptimizationStep=ScaleOptimizationStep;
			else
				std::cerr<<"ERROR: SetScaleOptimizationStep(): Wrong value:"<<ScaleOptimizationStep<<std::endl;
		}

	void setMatrixMinDoseThreshold(double MatrixMinDoseThreshold) {
			nMatrixMinDoseThreshold=MatrixMinDoseThreshold;
		}

	void setMinParticlesPerPb(int MinParticlesPerPb) {
			nMinParticlesPerPb=MinParticlesPerPb;
		}

	void setIsInitialEvaluation(bool IsInitialEvaluation) {
			isInitialEvaluation=IsInitialEvaluation;
		}

	void setIsCheckPb(bool IsCheckPb) {
			isCheckPb=IsCheckPb;
		}

	void setRoiFiles(std::string File) {
			if (File.empty())
				ROI = Roi("ROI.txt");
			else
				ROI = Roi(File);
		}

	void setPencilBeamFiles(std::string File) {
			if (File.empty())
				PB= PencilBeam("raster.txt");
			else
				PB= PencilBeam(File);
		}

	void setVxlFiles(std::string File, float nXmin, float nYmin, float nZmin,
			 std::string evtFile, std::vector<int> numRoi,
			 std::vector<int> isPtv, std::vector<float> TargetDose,
			 std::vector<float> TargetRBE, std::vector<float> TargetLET) {
			if (File.empty())
				ROI= Roi("ROI.txt");
			else
				ROI= Roi(File, nXmin, nYmin, nZmin, evtFile, numRoi, isPtv,
					TargetDose, TargetRBE, TargetLET);
		}

	void setDoseFiles(std::string File) {
			vFileDose=File;
		}

	void setAlphaFiles(std::string File) {
			vFileAlpha=File;
		}

	void setBetaFiles(std::string File) {
			vFileBeta=File;
		}

	void setLetFiles(std::string File) {
				vFileLet=File;
			}

	void setErrorMatrixFiles(std::vector<std::string> Files, std::vector<int> RoiBoxId) { //TODO: Error Matrix Files
			vFileErrorMatrix=Files;
			vFileErrorMatrixRoiBoxId=RoiBoxId;
		}

	void doReadNewIntensities(std::string FNewIntensities,bool IsReadNewIntensitiesBefore) {
			isReadNewIntensities=true;
			isReadNewIntensitiesBefore=IsReadNewIntensitiesBefore;
			FileNewIntensities=FNewIntensities;
			std::cout<<FNewIntensities<<std::endl;
		}

	void setScaleIntensities(double ScaleIntensities) {
			nScaleIntensities=ScaleIntensities;
		}

	void setRbeFactor(double RbeFactor) {
			nRbeFactor=RbeFactor;
		}

	void setRbeFixed(double RbeFixed) {
			nRbeFixed=RbeFixed;
			if (fabs(nRbeFixed-1.)>0.0001) {
				isRbeFixed=true;
			}
		}

	void setRbeInPtv(double RbeInPtv) {
			nTargetRbeMean=RbeInPtv;
			if (fabs(RbeInPtv)>0.0001) {
				isRbeInPtvFixed=true;
			}
		}

	void setLetInPtv(double LetInPtv) {//target LET
			nTargetLetMean=LetInPtv;
			if (fabs(LetInPtv)>0.0001) {
				isLetInPtvFixed=true;
			}
		}

	void setPhotonEquivSurvival(bool PhotonEquivSurvival) {
		       isPhotonEquivSurvival=PhotonEquivSurvival;
		}

	void setSelectiveEvaluation(bool SelectiveEvaluation) {
			isSelectiveEvaluation=SelectiveEvaluation;
		}

	void setIsOptimize1d(bool IsOptimize1d) {
			isOptimize1d=IsOptimize1d;
		}

	void setIterationActionInterval(int IterationActionInterval) {
			if (IterationActionInterval>0) {
				nIterationActionInterval=IterationActionInterval;
			} else {
				std::cerr<<"ERROR: SetIteractionActionInterval(): Wrong input <=0!"<<std::endl;
				exit(-1); // FIXME
			}
		}

	void setDoPlot(bool DoPlot) {
			isPlot=DoPlot;
		}

	void setPlotSlices(int PlotSlices) {
			if ((PlotSlices>=0)&&(PlotSlices<=3)) {
				nPlotSlices=PlotSlices;
			} else {
				std::cerr<<"ERROR: SetPlotSlices(): Wrong input value!"<<std::endl;
				exit(-1); // FIXME
			}
		}

	void setPlotRoiTest(int PlotRoiTest) {
			if ((PlotRoiTest>=0)&&(PlotRoiTest<=3)) {
				nPlotRoiTest=PlotRoiTest;
			} else {
				std::cerr<<"ERROR: SetPlotRoiTest(): Wrong input value!"<<std::endl;
				exit(-1); // FIXME
			}
		}

	void setPlotColour(bool DoPlotColour) {
			isPlotColour=DoPlotColour;
		}

	/*--------------------------------------------------------------------------------------------------------------*/
	//********************************************** Getters Functions **********************************************
	bool	getIsBiologicalOptimization() const {return isBiologicalOptimization;}
	int	getOptimizationAlgorithm() const {return nOptimizationAlgorithm;}
	int	getCellLine() const {return nCellLine;}

	double	getConvergenceCriterium() const	{ return nConvergenceCriterium; }
	double	getScaleOptimizationStep() const	{ return nScaleOptimizationStep; }
	double	getMatrixMinDoseThreshold() const{ return nMatrixMinDoseThreshold; }
	int	getMinParticlesPerPb()	const { return nMinParticlesPerPb; }

	bool	getIsInitialEvaluation() const{ return isInitialEvaluation; }
	bool	getIsCheckPb()		const { return isCheckPb; }

	std::string getDoseFiles()	const { return vFileDose; }
	std::string getAlphaFiles()	const { return vFileAlpha; }
	std::string getBetaFiles()	const { return vFileBeta; }
	std::string getLetFiles()	const { return vFileLet; }

	std::vector<std::string> getErrorMatrixFiles()	const { return vFileErrorMatrix;}; //TODO: Change Error MAtrix Files into class
	std::vector<int> getErrorMatrixFilesRoiBoxId()	const { return vFileErrorMatrixRoiBoxId;}
	bool	getIsReadNewIntensities()		const { return isReadNewIntensities;}
	bool	getIsReadNewIntensitiesBefore()		const { return isReadNewIntensitiesBefore;}
	std::string getReadNewIntensities()		const { return FileNewIntensities;}
	double	getScaleIntensities()			const { return nScaleIntensities;}

//	double	getPrescribedPtvDose()	const {return nPrescribedPtvDose;}
	double	getRbeFactor()		const {return nRbeFactor;}
	double	getRbeFixed()		const {return nRbeFixed;}
	double	getTargetRbeMean()	const {return nTargetRbeMean;}

	bool	getPhotonEquivSurvival() const {return isPhotonEquivSurvival;}
	bool	getSelectiveEvaluation() const {return isSelectiveEvaluation;}

	bool	getIsOptimize1d()	const {return isOptimize1d;}
	int	getIterationActionInterval() const {return nIterationActionInterval;}
}; // class Optimizer

#endif
