 // ********************************************
// Author:  wioletta.kozlowska@cern.ch
// based on the T.Boehlen version from 2011
// Modified wioletta.kozlowska@cern.ch
// Version: 3.1
// Last change: 19/09/2016
// *********************************************

#include <stdlib.h>

#include "roi.h"

using namespace std;

/** Roi */
Roi::Roi(string filename)
{
	fileROI=filename;
	RoiNb=0;
	readRoi(filename); //RoiBoxId
} // Roi

/** Roi */
Roi::Roi(string filename, float nXmin, float nYmin, float nZmin,
	 string evtBin, vector<int> numRoi, vector<int> isPTV,
	 vector<float> targetDose, vector<float> targetRBE,vector<float> targetLET)
{
	vxlFile = filename;
	evtFile = evtBin;
	for (unsigned int idx=0; idx<numRoi.size(); idx++) {
		ROI CurrentRoi;
		CurrentRoi.RoiId=numRoi[idx];
		CurrentRoi.Id=isPTV[idx];
		CurrentRoi.TargetDose=targetDose[idx];
		CurrentRoi.TargetRBE=targetRBE[idx];
		CurrentRoi.TargetLET=targetLET[idx];
		CurrentRoi.VoxelIndexMax=0;
		box.push_back(CurrentRoi);
	}
	RoiNb=numRoi.size();
	readVxl(vxlFile, nXmin, nYmin, nZmin, evtFile); //RoiBoxId,
} // Roi


/** Roi */
Roi::~Roi()
{
	for (unsigned j=0; j<box.size(); j++)
		for (unsigned i=0; i<box[j].MatrixIndex.size(); i++)
			delete [] box[j].MatrixIndex[i];
} // Roi

/** clean */
int Roi::clean(vector<ROI> *Struct)
{
	for (unsigned int i0=0; i0<Struct->size(); i0++) {
		(*Struct)[i0].VoxelIndexMax=0;
		(*Struct)[i0].Id=0;
		(*Struct)[i0].VoxelIndex.clear();
		(*Struct)[i0].TargetDose=0;
		(*Struct)[i0].TargetRBE=0;
		(*Struct)[i0].TargetLET=0;
		(*Struct)[i0].Weight.clear();
		//TODO: add free and clear matrix index
	}
	Struct->clear();
	return 0;
} // clean

/** readVxl */
int Roi::readVxl(string filename, float nXmin, float nYmin, float nZmin, string evtBin)
{
	cout<<"Reading Roi from Voxel File "<<filename<<"..."<<endl;
	for (unsigned idx=0;idx<box.size() ;idx++ ) {
		cout<<"For Roi Id no "<<box[idx].RoiId<<"  "<<endl;
	}

	ifstream FileIn(filename.c_str() );

	cout<<"READ <"<<filename<<"> ...";
	cout.flush();

	if (!FileIn.good()) {
		cerr<<"\nERROR: Roi::ReadRoi(): I/O error!\n"<<endl;
		exit(-1);	// FIXME
	}

	Eventbin* Binning = new Eventbin(evtBin);
	GVoxel* Voxel = new GVoxel();

	BinNb.set(Binning->nx, Binning->ny, Binning->nz);
	Min.set(Binning->xlow, Binning->ylow, Binning->zlow);
	Max.set(Binning->xhigh, Binning->yhigh, Binning->zhigh);
	Delta.set((Max.x-Min.x)/((double)BinNb.x),
		  (Max.y-Min.y)/((double)BinNb.y),
		  (Max.z-Min.z)/((double)BinNb.z));

	if (!Voxel->load(filename.c_str())) {
		 cerr<<"ERROR: loading voxel file"<<filename<<endl;
	}

	Voxel->calcLimits();

	cout<<"voxel nx "<<Voxel->nx<<" voxel ny "<<Voxel->ny
		<<" voxel nz "<<Voxel->nz<<endl;

	Voxel->xlow=nXmin;
	Voxel->ylow=nYmin;
	Voxel->zlow=nZmin;
	cout<<"voxel xlow "<<Voxel->xlow<<" voxel ylow "<<Voxel->ylow
		<<" voxel zlow "<<Voxel->zlow<<endl;
	cout.flush();


	cout<<"binning nx "<<Binning->nx<<" binning ny "<<Binning->ny
	    <<"binning nz "<<Binning->nz<<endl;
	for (int idx=0; idx<RoiNb; idx++) {
		set<int> roinb;
		box[idx].VoxelIndexMax=Binning->nx*Binning->ny*Binning->nz;
		 for (int k=0; k<Voxel->nz; k++) {
			 double z  = Voxel->voxelcz(k);
			 double zz = Binning->eventbink(z);
			 for (int j=0; j<Voxel->ny; j++) {
				 double y  = Voxel->voxelcy(j);
				 double yy = Binning->eventbinj(y);

				 for (int i=0; i<Voxel->nx; i++) {
					 double x  = Voxel->voxelcx(i);
					 double xx = Binning->eventbini(x);
					 const ROICombination& comb = Voxel->roiComb(i,j,k);
					 for (int m=0; m<comb.length; m++) {
						int dummy = Abs(comb[m]);
						 if (dummy==box[idx].RoiId) {
							 int binnb=xx+(yy)*Binning->nx+(zz)*Binning->nx*Binning->ny;
							roinb.insert(binnb);
						 }
					 }
				 }
			 }
		 }
		 copy(roinb.begin(), roinb.end(), back_inserter(box[idx].VoxelIndex));

		 for (unsigned int id=0;id<box[idx].VoxelIndex.size();id++) {
			box[idx].Weight.push_back(1.);
		 }
	}
	VoxelVolume=Delta.x*Delta.y*Delta.z;
	cout<<VoxelVolume<<" ";

	delete(Voxel);
	delete(Binning);
	cout<<" DONE"<<endl;

	cout<<"### ROI Summary ###################################"<<endl;
	cout<<"  Voxel dimension:("<<Delta.x<<"cm,"<<Delta.y<<"cm,"<<Delta.z<<"cm)"<<endl;
	cout<<"  Usrbin minimum: ("<<Min.x<<"cm,"<<Min.y<<"cm,"<<Min.z<<"cm)"<<endl;
	cout<<"  Usrbin maximum: ("<<Max.x<<"cm,"<<Max.y<<"cm,"<<Max.z<<"cm)"<<endl;
	cout<<"  Usrbin bins:    ("<<BinNb.x<<","<<BinNb.y<<","<<BinNb.z<<")"<<endl;
	cout<<"  Voxel volume:    "<<VoxelVolume<<"cm³"<<endl;
	cout<<"  Roi to be read:  "<<RoiNb<<endl;
	for (int idx=0;idx<RoiNb;idx++) {
		cout<<"###################################################"<<endl;
		cout<<"  ROI ID:          "<<box[idx].RoiId<<endl;
		cout<<"  ROI tag:         "<<box[idx].Tag<<endl;
		cout<<"  isPTV[1]:        "<<box[idx].Id<<endl;
		cout<<"  TargetDose:      "<<box[idx].TargetDose<<" Gy"<<endl;
		cout<<"  TargetRbe:       "<<box[idx].TargetRBE<<" Gy"<<endl;
		cout<<"  TargetLET:       "<<box[idx].TargetLET<<endl;
		cout<<"  ROI voxels:      "<<box[idx].VoxelIndex.size()<<endl;
	}
	cout<<"###################################################"<<endl;
	return 0;
} // readVxl

/*--------------------------------------------------------*/
/** readRoi */
int Roi::readRoi(string FileName)
{
	cout<<"Reading Roi from File "<<FileName<<"..."<<endl;
	cout.flush();

	ifstream FileIn(FileName.c_str() );

	cout<<"READ <"<<FileName<<"> ...";
	cout.flush();

	if (!FileIn.good()) {
		cerr<<"\nERROR: Roi::ReadRoi(): I/O error!\n"<<endl;
		exit(-1);	// FIXME
	}

	string Line;
	double Help[5];
	int VoxelToBeRead;

	if (FileIn.is_open()) {
			while (FileIn.good() && getline(FileIn,Line))
				if (Line[0]!='*'&&!(Line.empty()))
					break;
			RoiNb=0;
			RoiNb=atoi(Line.c_str());
			FileIn>>Help[0]>>Help[1]>>Help[2];
			BinNb.set(Help[0],Help[1],Help[2]);
			FileIn>>Help[0]>>Help[1]>>Help[2];
			Min.set(Help[0],Help[1],Help[2]);
			FileIn>>Help[0]>>Help[1]>>Help[2];
			Max.set(Help[0],Help[1],Help[2]);

			int CountRoi=0;
			while(FileIn.good()&& (CountRoi<RoiNb)) {
				ROI CurrentRoi;

				//stringstream buffer(Line);
				FileIn>>CurrentRoi.Tag>>CurrentRoi.RoiId>>CurrentRoi.Id
				>>CurrentRoi.TargetDose>>CurrentRoi.TargetRBE>>CurrentRoi.TargetLET;
				FileIn>>VoxelToBeRead;

				//read file body
				int CountVoxel=0;
				//int MaxVoxel=0;
				while(FileIn.good()&&(CountVoxel<VoxelToBeRead)) {
				// <voxelindex> <voxelweight>
					FileIn>>Help[0]>>Help[1];
					CurrentRoi.VoxelIndex.push_back((int)Help[0]);  //do indexing (dose grid index)
					CountVoxel++;
					CurrentRoi.Weight.push_back(int(Help[1]));
					//MaxVoxel=max(int(Help[0]),MaxVoxel);
				}
				CurrentRoi.VoxelIndexMax=BinNb.x*BinNb.y*BinNb.z;//MaxVoxel;
				CurrentRoi.VoxelNb=CountVoxel;
				CountRoi++;
				cout<< CountVoxel<<endl;;
				box.push_back(CurrentRoi);
			}
	}
	FileIn.close();

	Delta.set((Max.x-Min.x)/((double)BinNb.x),
		  (Max.y-Min.y)/((double)BinNb.y),
		  (Max.z-Min.z)/((double)BinNb.z));

	//this is from old roi - check which one is correct

	VoxelVolume=Delta.x*Delta.y*Delta.z;
	cout<<" DONE"<<endl;

	cout<<"### ROI Summary ###################################"<<endl;
	cout<<"  Voxel dimension:("<<Delta.x<<"cm,"<<Delta.y<<"cm,"<<Delta.z<<"cm)"<<endl;
	cout<<"  Usrbin minimum: ("<<Min.x<<"cm,"<<Min.y<<"cm,"<<Min.z<<"cm)"<<endl;
	cout<<"  Usrbin maximum: ("<<Max.x<<"cm,"<<Max.y<<"cm,"<<Max.z<<"cm)"<<endl;
	cout<<"  Usrbin bins:    ("<<BinNb.x<<","<<BinNb.y<<","<<BinNb.z<<")"<<endl;
	cout<<"  Voxel volume:    "<<VoxelVolume<<"cm³"<<endl;
	cout<<"  Roi to be read:  "<<RoiNb<<endl;
	for (int idx=0;idx<RoiNb;idx++) {
		cout<<"###################################################"<<endl;
		cout<<"  ROI ID:          "<<box[idx].RoiId<<endl;
		cout<<"  ROI tag:         "<<box[idx].Tag<<endl;
		cout<<"  isPTV[1]:        "<<box[idx].Id<<endl;
		cout<<"  TargetDose:      "<<box[idx].TargetDose<<" Gy"<<endl;
		cout<<"  TargetRbe:       "<<box[idx].TargetRBE<<" Gy"<<endl;
		cout<<"  TargetLET:       "<<box[idx].TargetLET<<endl;
		cout<<"  ROI voxels:      "<<box[idx].VoxelIndex.size()<<endl;
	}
	cout<<"###################################################"<<endl;
	return 0;
} // readRoi

/*--------------------------------------------------------*/
/** index2Position */
Vector Roi::index2Position(int VoxelIndex)
{	//was inline
	VECTOR_INT Bin;

	//TODO: add check here

	Bin=index2Bin(VoxelIndex);

	return Vector(
			Min.x+((double)Bin.V[0]-0.5)*Delta.x,
			Min.y+((double)Bin.V[1]-0.5)*Delta.y,
			Min.z+((double)Bin.V[2]-0.5)*Delta.z
	);
} // index2Position

/*--------------------------------------------------------*/
int Roi::position2Index(Vector Position, int Id)
{	//TEST IT
	int CurrentBin[3];

	//We do rounding here
	CurrentBin[0]=  (int)( (Position.x-Min.x)/Delta.x );
	CurrentBin[1]=  (int)( (Position.y-Min.y)/Delta.y );
	CurrentBin[2]=  (int)( (Position.z-Min.z)/Delta.z );

	if ( isInsideRoi(CurrentBin) ) {
		return bin2Index(CurrentBin, Id);
	}
	else{
	  cerr<<"ERROR: Optimizer::Position2Index(): Index out-of-bounds in ROI ID "<<Id<<": "<<CurrentBin[0]<<" "<<CurrentBin[1]<<" "<<CurrentBin[2]<<" Max values: "<<BinNb.x<<" "<<BinNb.y<<" "<<BinNb.z<<endl;
		cerr<<"  Position: "<<Position.x<<" "<<Position.y<<" "<<Position.z<<endl;
		exit(-1);	// FIXME
	}
} // position2Index

/*--------------------------------------------------------*/
int Roi::bin2Index(int CurrentBin[3], int Id)
{	//TEST IT
	int VoxelIndex;

	if ( isInsideRoi(CurrentBin) ) {
	  VoxelIndex= CurrentBin[2]*BinNb.x*BinNb.y+CurrentBin[1]*BinNb.x+CurrentBin[0]+1;
		return VoxelIndex;
	} else {
		cerr<<"Optimizer::Bin2Index(): Index out-of-bounds in ROI ID "<<Id<<": "<<CurrentBin[0]<<" "<<CurrentBin[1]<<" "<<CurrentBin[2]<<" Max values: "<<BinNb.x<<" "<<BinNb.y<<" "<<BinNb.z<<endl;
			exit(-1);	// FIXME
	}
} // bin2Index

/*--------------------------------------------------------*/
/** index2Bin */
VECTOR_INT Roi::index2Bin(int VoxelIndex)
{	 // TEST IT
	VECTOR_INT Bin;

	//TODO: add check here
	Bin.V[2]=(int)( (VoxelIndex-1)/(BinNb.x*BinNb.y) )+1;
	Bin.V[1]=(int)( (VoxelIndex-1- (Bin.V[2]-1)*(BinNb.x*BinNb.y))/BinNb.x )+1;
	Bin.V[0]=VoxelIndex- ( (Bin.V[2]-1)*(BinNb.x*BinNb.y) + (Bin.V[1]-1)*BinNb.x );

	return Bin;
} // index2Bin

/*--------------------------------------------------------*/
bool Roi::isInsideRoi(int Bin[3])
{	// TEST IT needs to be redone definitely
	// WARNING: the binning start with 1 (not 0!)???
	if ( (Bin[0]>=0)&&(Bin[0]<BinNb.x) &&
	    (Bin[1]>=0)&&(Bin[1]<BinNb.y) &&
	    (Bin[2]>=0)&&(Bin[2]<BinNb.z))
		return true;
	else
		return false;
} // isInsideRoi

//not used TODO: check it
/*
 * Manipulates voxel weights and target (bio-)doses of ROI
 */
/*--------------------------------------------------------
int Roi::manipulateRoi(int RoiBoxId)
{
	//cout<<"In ManipulateRoi  function"<<PBFiles[0].PBBoxSet[0].Intensity.back()<<endl;
	cout<<"Manipulate ROI ...";
	cout.flush();

	double NewWeight=1.;
	double NewDose=3.;//Gy (RBE)

	//WARNING: works only for boxes at the moment
	Vector Min=setRoi[RoiBoxId].Max;
	Vector Max=setRoi[RoiBoxId].Min;

	//find min and max bin in Z
	for (unsigned int i0=0;i0<setRoi[RoiBoxId].VoxelIndex.size();i0++) {
		Vector Pos=index2Position(setRoi[RoiBoxId].VoxelIndex[i0],RoiBoxId);

		if (setRoi[RoiBoxId].Id[i0]>0) {//if PTV
			if (Min.z>Pos.z) Min.setZ(Pos.z);
			if (Max.z<Pos.z) Max.setZ(Pos.z);
		}

		if (setRoi[RoiBoxId].Id[i0]>0) setRoi[RoiBoxId].TargetDose[i0]=NewDose;//set new target dose
	}

	//set new weight
	for (unsigned int i0=0;i0<setRoi[RoiBoxId].VoxelIndex.size();i0++) {
		Vector Pos=index2Position(setRoi[RoiBoxId].VoxelIndex[i0],RoiBoxId);

		if (setRoi[RoiBoxId].Id[i0]>0) {//if PTV
			if (fabs(Min.z-Pos.z)<0.01) setRoi[RoiBoxId].Weight[i0]=NewWeight;
			if (fabs(Max.z-Pos.z)<0.01) setRoi[RoiBoxId].Weight[i0]=NewWeight;
		}
	}

	cout<<" DONE"<<endl;

	cout<<"  Set a target dose of "<<NewDose<<"Gy (RBE) in all PTVs"<<endl;
	cout<<"  Set a weight of "<<NewWeight<<" in all min and max PTV voxels in z (only working for boxes)."<<endl;
	return 0;
} // manipulateRoi
*/
