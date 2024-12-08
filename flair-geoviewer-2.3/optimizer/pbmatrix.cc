// ********************************************
// Author:  wioletta.kozlowska@cern.ch
// based on the T.Boehlen version from 2011
// Modified wioletta.kozlowska@cern.ch
// Version: 3.1
// Last change: 24/08/2016
// *********************************************

#include <iomanip>
#include <iostream>

#include "roi.h"
#include "pbmatrix.h"
#include "pencilbeam.h"

using namespace std;

/*--------------------------------------------------------*/
int PBMatrix::clean(vector< vector<PB_MATRIX> >& Struct)
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"PBMAtrix::clean"<<endl;
	cout<<"---------------------------------------------------------"<<endl;
	for (unsigned i0=0; i0<Struct.size(); i0++) {
		for (unsigned i1=0; i1<Struct[i0].size(); i1++) {
			if (Struct[i0][i1].VoxelIndex)     delete [] Struct[i0][i1].VoxelIndex;
			if (Struct[i0][i1].DosePerPrimary) delete [] Struct[i0][i1].DosePerPrimary;
			if (Struct[i0][i1].AlphaMean)      delete [] Struct[i0][i1].AlphaMean;
			if (Struct[i0][i1].SqrtBetaMean)   delete [] Struct[i0][i1].SqrtBetaMean;
			if (Struct[i0][i1].LetMean)        delete [] Struct[i0][i1].LetMean;
		}
		Struct[i0].clear();
	}
	Struct.clear();
	return 0;
} // clean

/*--------------------------------------------------------*/
int PBMatrix::readBinaryMatrix(PencilBeam* PB, Roi* roi, string vFileDose, string vFileAlpha,
		string vFileBeta, string vFileLet, double nMatrixMinDoseThreshold,
		double nRbeFixed, bool isLetOptimization, bool isBiologicalOptimization,
		vector<DOSE_GRID> DoseGrid)
{
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"PBMatrix::readBinaryMatrix"<<endl;
	cout<<"---------------------------------------------------------"<<endl;

	//check which files are existing
	bool IsDose  = (!vFileDose.empty());
	bool IsAlpha = (!vFileAlpha.empty());
	bool IsBeta  = (!vFileBeta.empty());
	bool IsLet   = (!vFileLet.empty());

	if (!IsDose) {
		cerr<<"\nERROR: PBMatrix::ReadBinaryMatrix(): No dose to read!"
		    <<" vFileDose has no entries for  (<"<<vFileDose<<")!\n"<<endl;
		exit(-1);	// FIXME
	}

	if (isBiologicalOptimization&&(!IsAlpha || !IsBeta)) {
		cerr<<"\nERROR: PBMatrix::ReadBinaryMatrix(): "
		    <<"Biological optimization requested! Please provide Alpha and Beta files\n"<<endl;
		exit(-1);	// FIXME
	}

	if (isLetOptimization&&(!IsLet)) {
		cerr<<"\nERROR:PBMatrix::ReadBinaryMatrix(): "
		    <<"LET optimization requested! Please provide LET files\n"<<endl;
		exit(-1);	// FIXME
	}

	if ((IsAlpha && !IsDose) ||	// FIXME no reason for IsDose is protected above
	    (IsBeta  && !IsDose) ||
	    (IsLet   && !IsDose)) {
		cerr<<"\nERROR: PBMatrix::ReadBinaryMatrix(): "
		    <<"You are missing file for alpha, beta and LET!\n"<<endl;
		exit(-1);	// FIXME
	}

	if (PB->getFilePB().empty()) { //check it
		cerr<<"\nERROR: PBMatrix::ReadBinaryMatrix(): Call first ReadPencilBeam()!\n"<<endl;
		exit(-1);	// FIXME
	}

	//check if RoiBoxId already exists, if not, create extra entry(ies)
	if (box.size()<=(unsigned)roi->box.size()) {
		vector<PB_MATRIX> DummyPbMatrix;
		for (unsigned i0=box.size(); i0<(unsigned)roi->box.size(); i0++) {
			box.push_back(DummyPbMatrix);
		}
	}

	//figure out the number of PB to read from the current file
	for (unsigned RoiBoxId=0; RoiBoxId<roi->box.size(); RoiBoxId++) {
		int PbToBeRead=0;
		if (PbToBeRead==(int)box[RoiBoxId].size()) {
			PbToBeRead=PB->NumPB();
		}

		if (PbToBeRead==0) {
			cerr<<"\nERROR: PBMatrix::ReadBinaryMatrix(): No match of PbToBeRead! Something went wrong\n"<<endl;
			cerr<<setw(3)<<" PB per file:"<<PB->NumPB()<<endl;
			exit(-1);	// FIXME
		}

		//TODO: check that PbToBeRead is correct!!
		cout<<"READ FOR ROI BOX "<<RoiBoxId<<" AND "<<PbToBeRead<<" PBs MATRIX OF:"<<endl;

		cout<<"  <dose> FROM <"<<vFileDose<<endl;
		Eventbin Dose(vFileDose);
		Eventbin Alpha, Beta, Let;
		if (IsAlpha) {
			Alpha.load(vFileAlpha);
			cout<<endl<<"  <alpha> FROM <"<<vFileAlpha<<">";
		}
		if (IsBeta) {
			Beta.load(vFileBeta);
			cout<<endl<<", <beta>  FROM <"<<vFileBeta<<">";
		}
		if (IsLet) {
			Let.load(vFileLet);
			cout<<endl<<", <LET>   FROM <"<<vFileLet<<">";
		}
		cout<<" ...";
		cout.flush();

		//////////////
		char	*pbufferedDose=NULL, *pbufferedAlpha=NULL, *pbufferedBeta=NULL, *pbufferedLet=NULL;
		char	BufferedDose[200000];
		char    BufferedAlpha[200000], BufferedBeta[200000], BufferedLet[200000];

		int CountPbRead=0;
		int CountElements=0;
		int CountReducedElements=0;
		int CountDeactivatedPb=0;
		double OverallMaxDoseInVoxel=0.;
		int CountSkipBetaPerFile=0;

		for (int i0=0;i0<PbToBeRead;i0++) {
			PB_MATRIX CurPbMatrix;

			// set to null everything
			memset(&CurPbMatrix, 0, sizeof(CurPbMatrix));

			//store in matrix
			CurPbMatrix.PbIndex=box[RoiBoxId].size();

			//Read event and return pointer
			pbufferedDose=Dose.readEvent(BufferedDose, sizeof(BufferedDose));
			CurPbMatrix.VoxelHit=Dose.nhits;
			CurPbMatrix.Events=(int)Dose.weight;
			if (IsAlpha) {
				pbufferedAlpha=Alpha.readEvent(BufferedAlpha, sizeof(BufferedAlpha));
				if (CurPbMatrix.VoxelHit!=Alpha.nhits ||
					CurPbMatrix.Events!=(int)Alpha.weight) {
					cerr<<"\nERROR: PBMatrix::ReadBinaryMatrix(): "
					<<"Prior matrix entries (DOSE) are not the same! VoxelHit: "
					<<CurPbMatrix.VoxelHit<<"!="<<Alpha.nhits<<" or "
					<<CurPbMatrix.Events<<"!="<<(int)Alpha.weight<<" !\n"<<endl;
					exit(-1);	// FIXME
				}
			}
			if (IsBeta) {
				pbufferedBeta=Beta.readEvent(BufferedBeta,sizeof(BufferedBeta));
				if (CurPbMatrix.Events!=(int)Beta.weight) {
						cerr<<"\nERROR: PBMatrix::ReadBinaryMatrix(): "
						<<"Prior matrix entries (DOSE) do not have the same weight: "
						<<CurPbMatrix.Events<<"!="<<(int)Beta.weight<<" ! VoxelHit (DOSE): "
						<<CurPbMatrix.VoxelHit<<" and now: "<<Beta.nhits<<"\n"<<endl;
						exit(-1);	// FIXME
				}
				//if ((Verbose>0)&&((*PbMatrix).VoxelHit!=nhits)) {
				//	cerr<<"\nWARNING: Optimizer::ReadEventbinEvent(): VoxelHit (DOSE): "
				//	<<(*PbMatrix).VoxelHit<<" and now (BETA): "
				//	<<nhits<<" This is probably due to beta=0.0 for high LET in the model which were skipped!\n"<<endl;
				//}
				CurPbMatrix.TempVoxelHitBeta=Beta.nhits;//this is in case of less beta entries compared to dose!!!
			}
			if (IsLet) {
				pbufferedLet=Let.readEvent(BufferedLet,sizeof(BufferedLet));
				if (CurPbMatrix.VoxelHit!=Let.nhits ||
						CurPbMatrix.Events!=(int)Let.weight) {
						cerr<<"\nERROR: PBMatrix::ReadBinaryMatrix(): "
						<<"Prior matrix entries (DOSE) are not the same! VoxelHit: "
						<<CurPbMatrix.VoxelHit<<"!="<<Let.nhits<<" or "
						<<CurPbMatrix.Events<<"!="<<(int)Let.weight<<" !\n"<<endl;
						exit(-1);	// FIXME
				}
			}

			if ( (unsigned)CurPbMatrix.PbIndex > PB->NumPB() || CurPbMatrix.PbIndex<0 ) {
				cerr<<"\nERROR: PBMatrix::ReadBinaryMatrix(): CurPbMatrix.PbIndex>fPencilBeam.size(): "
				    <<CurPbMatrix.PbIndex<<">"<<PB->NumPB()<<"!\n"
				    <<CurPbMatrix.PbIndex<<endl;
				exit(-1);	// FIXME
			}

			CountPbRead++;

			//Prepare arrays
			CurPbMatrix.VoxelIndex= new int[CurPbMatrix.VoxelHit];
			CurPbMatrix.DosePerPrimary= new float[CurPbMatrix.VoxelHit];
			//To save time for non-bio calcs: only if fIsBiologicalOptimization or fIsLetOptimization
			if (IsAlpha) CurPbMatrix.AlphaMean= new float[CurPbMatrix.VoxelHit];
			if (IsBeta) CurPbMatrix.SqrtBetaMean= new float[CurPbMatrix.VoxelHit];
			if (IsLet) CurPbMatrix.LetMean= new float[CurPbMatrix.VoxelHit];

			int weight = CurPbMatrix.Events;//get weight from matrix again!

			int iFilled=0;
			double MaxDoseInVoxel=PB->MaxDosePerPrimary(CurPbMatrix.PbIndex);//set to the highest former one
			unsigned int MaxDoseVoxelIndex=PB->MaxDoseVoxelIndex(CurPbMatrix.PbIndex);
			unsigned int MaxDoseRoiBoxId=PB->MaxDoseRoiBoxId(CurPbMatrix.PbIndex);

			int dosevoxid;//=Dose.dosevoxid;
			float dose;//=Dose.dose;

			int alphavoxid=0, betavoxid=0,letvoxid=0;
			float alpha=0.0,beta=0.0,let=0.0;

			bool readbeta=true;//in case that a beta was skipped since it is equal to zero.
			int  countskipbeta=0;
			//cout<<"CurPbMatrix.VoxelHit "<<CurPbMatrix.VoxelHit<<endl;
			for (int iLine=0;iLine<CurPbMatrix.VoxelHit;iLine++) {
				//WATCHOUT: we assume that FLUKA provides with data in the form:
				//(FLUKA Voxel index)
				// Dose (in bin=Voxel volume) for all events(=weight) [Gy*cm^3]
				// | Dose*Alpha_Mean [cm^3] | Dose*Sqrt(Beta)_Mean [cm^3/Gy] | Dose*LET [Gy*GeV/cm]
				//then we need to convert for the optimizer to:
				// Dose/Primary [Gy] | Alpha_Mean [Gy] | Sqrt(Beta)_Mean [1/Gy^2] | LET_Mean [keV/um]

				FTNVGET(dosevoxid,int,pbufferedDose);
				FTNVGET(dose,float,pbufferedDose);

				if (IsAlpha&&IsBeta) {
					FTNVGET(alphavoxid,int,pbufferedAlpha);
					FTNVGET(alpha,float,pbufferedAlpha);
					alpha/=dose;

					if (readbeta) {//the beta from before was matched so read a new beta
						FTNVGET(betavoxid,int,pbufferedBeta);
						FTNVGET(beta,float,pbufferedBeta);
						//WATCH OUT: We can't do beta/=dose if values don't belong to the same voxid!!!
						//so do it later and only do:
						beta/=weight*roi->VoxelVolume;
					}
				}
				if (IsLet) {
					FTNVGET(letvoxid,int,pbufferedLet);
					FTNVGET(let,float,pbufferedLet);
					let*=1E+2/dose; // GeV/cm -> keV/um  (1E+6/1E+4)
				}

				dose/=weight*roi->VoxelVolume; //divide by weight and voxel volume to obtain: (Dose/Primary)

				//CHECK if dose and dosevoxid have valid values:
				if ( (dose<0.0)||(dosevoxid<=0)||(dosevoxid> (int)DoseGrid[RoiBoxId].RoiId.size()) ) {
					cout<<"\nERROR: PBMatrix::ReadBinaryMatrix(): PB#:"<<i0<<" Line:"<<iLine
					    <<" PbMatId:"<<iFilled <<" DG VoxelIndex out-of-bounds: 0< "<<dosevoxid
					<<" <"<<DoseGrid[RoiBoxId].RoiId.size()<<" or negative Dose/prim: "<<dose<<" or weight "<<weight<<" "<<endl;
					exit(-1);	// FIXME
				}

				if (MaxDoseInVoxel<dose) {//new max dose/primary ?
					MaxDoseInVoxel=dose;
					MaxDoseVoxelIndex=dosevoxid;
					MaxDoseRoiBoxId=RoiBoxId;
				}
				if (OverallMaxDoseInVoxel<dose) OverallMaxDoseInVoxel=dose;

				if (dose>nMatrixMinDoseThreshold) {//drop elements, if below threshold
					CurPbMatrix.VoxelIndex[iFilled]=dosevoxid;//do indexing
					CurPbMatrix.DosePerPrimary[iFilled]=dose*nRbeFixed;//multiply with a possible constant RBE scaling factor

					if (IsAlpha&&IsBeta) {
						if (dosevoxid!=alphavoxid) {
							cerr<<"ERROR: PBMatrix::ReadBinaryMatrix(): Alpha voxid not the same!"
							    <<endl;
							exit(-1);	// FIXME
						}
						CurPbMatrix.AlphaMean[iFilled]=alpha;

						if (dosevoxid!=betavoxid) {
							CurPbMatrix.SqrtBetaMean[iFilled]=0.0;
							readbeta=false;//try to match the current betavoxid with the next entry
							countskipbeta++;
						} else {
							CurPbMatrix.SqrtBetaMean[iFilled]=beta/dose;
							//do division by dose here, see above
							readbeta=true;
						}
					}
					if (IsLet) {
						if (dosevoxid!=letvoxid) {
							cerr<<"ERROR: PBMatrix::ReadBinaryMatrix(): LET voxid not the same!"
							    <<endl;
							exit(-1);	// FIXME
						}
						CurPbMatrix.LetMean[iFilled]=let;
					}

					iFilled++;
					CountElements++;
				} else {
					CountReducedElements++;
				}
			}//voxel loop

			if (IsBeta) {
				CountSkipBetaPerFile+=countskipbeta;
				if (CurPbMatrix.VoxelHit-CurPbMatrix.TempVoxelHitBeta!=countskipbeta) {
					cerr<<"\nERROR: Optimizer::ReadBinaryMatrix(): "
					    <<"CurPbMatrix.VoxelHit-CurPbMatrix.TempVoxelHitBeta!=countskipbeta, i.e.: "
					    << CurPbMatrix.VoxelHit <<"-"<<CurPbMatrix.TempVoxelHitBeta<<"="
					    <<CurPbMatrix.VoxelHit-CurPbMatrix.TempVoxelHitBeta<<"!="<<countskipbeta<<" !"
					    <<endl;
				exit(-1);	// FIXME
				}
			}

			CurPbMatrix.VoxelHitAboveThreshold=iFilled;

			box[RoiBoxId].push_back(CurPbMatrix);

			//store max dose voxel
			PB->setMaxDosePerPrimary(CurPbMatrix.PbIndex,MaxDoseInVoxel);
			PB->setMaxDoseVoxelIndex(CurPbMatrix.PbIndex,MaxDoseVoxelIndex);
			PB->setMaxDoseRoiBoxId(CurPbMatrix.PbIndex,MaxDoseRoiBoxId);

			//BUGFIX: Take out the following if{} (can give a problem if many ROIs are used!)
			//if no dosel was hit: de-activate PB
			//if (CurPbMatrix.VoxelHitAboveThreshold==0) {
			//	fPencilBeam[CurPbMatrix.PbIndex].Active=false;
			//	fPencilBeam[CurPbMatrix.PbIndex].Intensity[0]=0.;
			//	CountDeactivatedPb++;
			//	//if (fVerbose>0)
			//	cout<<"\n WARNING: No ROI voxel hit, de-activate PB#"<<CurPbMatrix.PbIndex<<endl;
			//}
		}//PB loop

		if (IsBeta&&CountSkipBetaPerFile!=0) {
			cout<<"\nWARNING: Optimizer::ReadBinaryMatrix(): Assumed "
			    <<CountSkipBetaPerFile<<" beta-values to be zero in this matrix!"<<endl;
		}

		cout<<" DONE"<<endl;
		cout<<"### Matrix Summary ################################"<<endl;
		cout<<"  ROI box ID: "<<roi->box[RoiBoxId].Tag<<endl;
		cout<<"  Read matrices of "<<CountPbRead<<" pencil beams"<<endl;
		cout<<"  Total matrices read: "<<box[RoiBoxId].size()<<endl;
		cout<<"  Matrix elements read:   "<<CountElements<<endl;
		cout<<"  Matrix reduced by: "<<CountReducedElements<<" elements (="
		    <<(double)CountReducedElements/(double)(CountReducedElements+CountElements)*100.
		    <<"% smaller)"<<endl;
		cout<<"  Matrix maximum dose per primary: "<<OverallMaxDoseInVoxel<<" Gy/primary"<<endl;
		cout<<"  Matrix minimum dose threshold:   "<<nMatrixMinDoseThreshold<<" Gy/primary"<<endl;
		cout<<"  De-activated pencil beams: "<<CountDeactivatedPb<<"/"<<PbToBeRead<<endl;
		if (IsLet)   cout<<"  LET read!"<<endl;
		if (IsAlpha) cout<<"  Alpha read!"<<endl;
		if (IsBeta)  cout<<"  Beta read!"<<endl;
		cout<<"###################################################"<<endl;
	}
	return 0;
} // readBinaryMatrix

/*--------------------------------------------------------*/
int PBMatrix::readErrorMatrix(string FileName, int RoiBoxId, PencilBeam *PB, Roi *roi,
			double nMatrixMinDoseThreshold, double nRbeFixed)
{
	//cout<<"In ReadErrorMatrix  function"<<vPbFiles[0].PB.BoxSet[0].Intensity.back()<<endl;
	PB_MATRIX CurPbMatrixError;
	double Help[4];

	if (PB->getFilePB().empty()) {
		cerr<<"\nERROR: PBMatrix::ReadErrorMatrix(): Call first ReadPencilBeam()!\n"<<endl;
		exit(-1);	// FIXME
	}

	//check if RoiBoxId already exists, if not, create extra entry(ies)
	if (BoxError.size()<=(unsigned int)RoiBoxId) {
		vector<PB_MATRIX> DummyPbMatrix;
		for (unsigned int i0=BoxError.size();i0<=(unsigned int)RoiBoxId;i0++) {
			BoxError.push_back(DummyPbMatrix);
		}
	}

	int PbNumberBefore=BoxError[RoiBoxId].size();

	//red out the number of PB to read from the current file
	int PbToBeRead=0;
	if (PbToBeRead==(int)BoxError[RoiBoxId].size()) {
		PbToBeRead=PB->NumPB();
	}
	if (PbToBeRead==0) {
		cerr<<"\nERROR: PBMatrix::ReadMatrixError(): No match of PbToBeRead! Something went wrong\n"<<endl;
		cerr<<setw(3)<<" PB per file:"<<PB->NumPB()<<endl;
		exit(-1);	// FIXME
	}

	cout<<"READ "<<PbToBeRead<<" PB FROM <"<<FileName<<"> ...";
	cout.flush();
	ifstream FileIn(FileName.c_str() );

	if (!FileIn.good()) {
		cerr<<"\nERROR: PBMatrix::ReadMatrixError(): I/O error!\n"<<endl;
		exit(-1);	// FIXME
	}

	int CountPbRead=0;
	int CountElements=0;
	int CountReducedElements=0;
	int CountDeactivatedPb=0;
	double OverallMaxDoseInVoxel=0.;
	while((FileIn.good())&&(CountPbRead<PbToBeRead)) {
		memset(&CurPbMatrixError, 0, sizeof(CurPbMatrixError));
		FileIn>>CurPbMatrixError.PbIndex>>CurPbMatrixError.VoxelHit>>CurPbMatrixError.Events;
		CurPbMatrixError.PbIndex=BoxError[RoiBoxId].size();//this assumes that the corresponding pencil beams were read in the same order!
		if ( (CurPbMatrixError.PbIndex>(int)PB->box.size())||(CurPbMatrixError.PbIndex<0) ) {//usually ">=" but since last entry is often messed up: ">"
			cerr<<"\nERROR: PBMatrix::ReadMatrix(): CurPbMatrixError.PbIndex>PB.box.size(): "<<CurPbMatrixError.PbIndex<<">"<<PB->box.size()<<"!\n"<<endl;
			exit(-1);	// FIXME
		}

		if (PbNumberBefore+CountPbRead!=CurPbMatrixError.PbIndex) {//just to check
			cerr<<"ERROR: PBMatrix::ReadMatrixError(): PbNumberBefore+CountPbRead!=CurPbMatrixError.PbIndex "<<PbNumberBefore<<"+"<<CountPbRead<<" "<<CurPbMatrixError.PbIndex<<endl;
			exit(-1);	// FIXME
		}
		CountPbRead++;

		//Prepare arrays
		CurPbMatrixError.VoxelIndex= new int[CurPbMatrixError.VoxelHit];
		CurPbMatrixError.DosePerPrimary= new float[CurPbMatrixError.VoxelHit];
		//TODO: to save time for non-bio calcs:  if (isBiologicalOptimization) {
		CurPbMatrixError.AlphaMean= new float[CurPbMatrixError.VoxelHit];
		CurPbMatrixError.SqrtBetaMean= new float[CurPbMatrixError.VoxelHit];

		int iFilled=0;
		double MaxDoseInVoxel=PB->box[CurPbMatrixError.PbIndex].MaxDosePerPrimary;//set to the highest former one
		int MaxDoseVoxelIndex=PB->box[CurPbMatrixError.PbIndex].MaxDoseVoxelIndex;
		int MaxDoseRoiBoxId=PB->box[CurPbMatrixError.PbIndex].MaxDoseRoiBoxId;
		for (int iLine=0;iLine<CurPbMatrixError.VoxelHit;iLine++) {
			//WATCHOUT: we need from the matrix:
			//Voxel index in ROI | E((Dose*(Voxel volume))^2) /Primary  [Gy*cm^3] | E(Alpha_Mean^2) [1/Gy] | E(Sqrt(Beta)_Mean^2) [1/Gy]
			FileIn>>Help[0]>>Help[1]>>Help[2]>>Help[3];

			Help[1]=Help[1]/(roi->VoxelVolume*roi->VoxelVolume);//divide per voxel volume^2 to obtain: (E(Dose^2)/Primary)

			//THIS IS NOT WORKING FOR THE MOMENT
			if (MaxDoseInVoxel<Help[1]) {//new max dose/primary ?
				MaxDoseInVoxel=Help[1];
				MaxDoseVoxelIndex=(int)Help[0];
				MaxDoseRoiBoxId=RoiBoxId;
			}
			if (OverallMaxDoseInVoxel<Help[1]) OverallMaxDoseInVoxel=Help[1];

//FOR THE MOMENT DON'T			if (Help[1]>nMatrixMinDoseThreshold) {//drop elements, if below threshold
			CurPbMatrixError.VoxelIndex[iFilled]=(int)Help[0];//do indexing
			CurPbMatrixError.DosePerPrimary[iFilled]=Help[1]*pow(nRbeFixed,2);//multiply with a possible fixed scaling factor
			CurPbMatrixError.AlphaMean[iFilled]=Help[2];
			CurPbMatrixError.SqrtBetaMean[iFilled]=Help[3];
			iFilled++;
			CountElements++;

		}//voxel loop
		CurPbMatrixError.VoxelHitAboveThreshold=iFilled;

		BoxError[RoiBoxId].push_back(CurPbMatrixError);

		//store max dose voxel
		PB->box[CurPbMatrixError.PbIndex].MaxDosePerPrimary=MaxDoseInVoxel;

		PB->box[CurPbMatrixError.PbIndex].MaxDoseVoxelIndex=MaxDoseVoxelIndex;
		PB->box[CurPbMatrixError.PbIndex].MaxDoseRoiBoxId=MaxDoseRoiBoxId;
	}

	FileIn.close();
	cout<<" DONE"<<endl;

	cout<<"### Error Matrix Summary ##########################"<<endl;
	cout<<"  ROI box ID: "<<RoiBoxId<<endl;
	cout<<"  Read matrices of "<<CountPbRead<<" pencil beams"<<endl;
	cout<<"  Total matrices read: "<<BoxError[RoiBoxId].size()<<endl;
	cout<<"  Matrix elements read:   "<<CountElements<<endl;
	cout<<"  Matrix reduced by: "<<CountReducedElements<<" elements (="<<(double)CountReducedElements/(double)(CountReducedElements+CountElements)*100.<<"% smaller)"<<endl;
	cout<<"  Matrix maximum dose per primary: "<<OverallMaxDoseInVoxel<<" Gy/primary"<<endl;
	cout<<"  Matrix minimum dose threshold:   "<<nMatrixMinDoseThreshold<<" Gy/primary"<<endl;
	cout<<"  De-activated pencil beams: "<<CountDeactivatedPb<<"/"<<PbToBeRead<<endl;
	cout<<"###################################################"<<endl;

	return 0;
} // readErrorMatrix
