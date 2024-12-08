// ********************************************
// Author:  wioletta.kozlowska@cern.ch
// based on the T.Boehlen version from 2011
// Modified wioletta.kozlowska@cern.ch
// Version: 3.1
// Last change: 24/08/2016
// *********************************************

#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

#include <stdlib.h>

#include "pencilbeam.h"

using namespace std;

//#TODO: Set up a values fixed on the input?
int PR=15;

PencilBeam::PencilBeam(string File)
{
	FilePB=File;
	TotalParticleNb=0;
	readPB(FilePB);
} // PencilBeam

PencilBeam::~PencilBeam(void)
{
	box.clear(); //it clears vector<vector<int>> only problem with *int
} // ~PencilBeam

/*--------------------------------------------------------------------*/
bool PencilBeam::readPB(string FileName)
{
	cout<<"Reading Pencil Beam from File "<<FileName<<"..."<<endl;
	cout.flush();

	ifstream FileIn;
	FileIn.open(FileName.c_str() );

	cout<<"READ <"<<FileName<<"> ...";
	cout.flush();

	if (!FileIn.good()) {
		cerr<<"\nERROR: PencilBeam::readPB(): I/O error!\n"<<endl;
		exit(-1);	// FIXME
	}

	setFilePB(FileName);

	string Line;
	string Dummy;
	double Help[9];
	PENCIL_BEAM CurrentPb;

	double nparticles=0;
	unsigned int PbNbBefore=NumPB();

	vector<unsigned int> Z,A;
	double EnMin(1.E30),EnMax(0.);

	if (FileIn.is_open()) {
		unsigned int card=0;

		while(FileIn.good()&&(getline(FileIn,Line))) {

			if (Line[0]=='*'||Line.empty())  continue;
			istringstream buffer(Line); //destructor necessary?

			//* - Spot Energy GeV/u, Momentum spread, Angular Spread,
			//* SpotSigmaX size, SpotSigma Y size, SpotWeight,
			//* xPos, yPos ,zPos,cosx,cosy,cosz, Xx, Xy, Xz,
			//* Atomic Number, Mass Number, BeamID
			// in 3 lines!

				/*buffer>>CurrentPb.Energy>>CurrentPb.MomSpread
			      >>CurrentPb.AngSpread>>CurrentPb.SpotSigma[0]
			      >>CurrentPb.SpotSigma[1]>>CurrentPb.InitialIntensity
			      >>Help[0]>>Help[1]>>Help[2]>>Help[3]>>Help[4]>>Help[5]
			      >>Help[6]>>Help[7]>>Help[8]>>CurrentPb.Z>>CurrentPb.A
			      >>CurrentPb.SpotId;*/

			getline(buffer,Dummy, ',');

			if (card%3==0) {

				/*FileIn>>Dummy>>a>>CurrentPb.Energy>>CurrentPb.MomSpread
			      >>CurrentPb.AngSpread>>CurrentPb.SpotSigma[0]
				      >>CurrentPb.SpotSigma[1]>>CurrentPb.InitialIntensity;*/

				getline(buffer,Dummy, ',');
				CurrentPb.Energy=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				CurrentPb.MomSpread=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				CurrentPb.AngSpread=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				CurrentPb.SpotSigma[0]=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				CurrentPb.SpotSigma[1]=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				CurrentPb.InitialIntensity=atof(Dummy.c_str());
				card++;

			}
			else if (card%3==1) {
				//		getline(buffer,Dummy, ',');
				getline(buffer,Dummy, ',');
				Help[0]=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				Help[1]=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				Help[2]=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				Help[3]=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				Help[4]=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				Help[5]=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				Help[6]=atof(Dummy.c_str());
				//buffer>>Help[0]>>Help[1]>>Help[2]>>Help[3]>>Help[4]>>Help[5]>>Help[6];
				card++;
			}
			else{
				//	getline(buffer,Dummy, ',');
				getline(buffer,Dummy, ',');
				Help[7]=atof(Dummy.c_str());
				getline(buffer,Dummy, ',');
				CurrentPb.Z=atoi(Dummy.c_str());
				getline(buffer,Dummy, ',');
				CurrentPb.A=atoi(Dummy.c_str());
				getline(buffer,Dummy, ',');
				CurrentPb.SpotId=atoi(Dummy.c_str());
				buffer>>Help[7]>>Help[8]>>CurrentPb.Z>>CurrentPb.A>>CurrentPb.SpotId;
				card++;


				nparticles+=CurrentPb.InitialIntensity;

				CurrentPb.Position.set(Help[0], Help[1], Help[2]);
				CurrentPb.DirCosines.set(Help[3], Help[4], Help[5]);
				CurrentPb.CSCosines.set(Help[6], Help[7], Help[8]);
				CurrentPb.Energy*=-1; // FLUKA requirements
				CurrentPb.Active=true;

				CurrentPb.MaxDoseVoxelIndex=-1;//fill later
				CurrentPb.MaxDoseRoiBoxId=-1;
				CurrentPb.MaxDosePerPrimary=-1.;

				CurrentPb.Intensity.push_back(CurrentPb.InitialIntensity);

				if (CurrentPb.Energy>EnMax) EnMax=CurrentPb.Energy;
				if (CurrentPb.Energy<EnMin) EnMin=CurrentPb.Energy;

				bool IsNewIon=true; //TODO: For what?

				for (unsigned int i0=Z.size();i0--;) {
					if ((Z[i0]==CurrentPb.Z)&&(A[i0]==CurrentPb.A))
						IsNewIon=false;
				}
				if (IsNewIon) {
					Z.push_back(CurrentPb.Z); A.push_back(CurrentPb.A);
				}

				addPB(CurrentPb);
			}
		};

		FileIn.close();
	};

	setTotPartNb(nparticles);

	cout<<"DONE"<<endl;

	cout<<"### Pencil Beam Summary ###########################"<<endl;
	cout<<"  Read "<<NumPB()-PbNbBefore<<" pencil beams"<<endl;
	cout<<"  Total number of pencil beams: "<<NumPB()<<endl;
	cout<<"  Initial total particle number: "<<getTotalParticleNb()<<endl;
	cout<<"  Ions of type: "<<endl;
	cout<<"  Z:";
	for (unsigned int i0=0;i0<Z.size();i0++) {
		cout<<setw(4)<<Z[i0]<<" ";
	}
	cout<<endl;
	cout<<"  A:";
	for (unsigned int i0=0;i0<A.size();i0++) {
		cout<<setw(4)<<A[i0]<<" ";
	}
	cout<<endl;
	cout<<"  Energy Min: "<<EnMin<<" GeV, Max: "<<EnMax<<" GeV"<<endl;
	cout<<"  WARNING: line above is now kinetic energy for p: MeV and for ions: MeV/u!"<<endl;
	cout<<"###################################################"<<endl;

	return 1;
} // readPB

/*--------------------------------------------------------*/
bool PencilBeam::writePB(string FileName, vector<PENCIL_BEAM>* ModifPB)
{
	cout<<"Write Pencil Beam to File <"<<FileName<<"> ...";
	cout.flush();

	ofstream FileOut( FileName.c_str() );
	double nparticles=0;

	//TODO: PB needs cleaning in case or erasing some PB spots?
	//calculate total number of particles in the plan
	for (unsigned int i0=(*ModifPB).size();i0--;) {// we are not sure if it is the same nb of pencil beams..
		setPB(&(*ModifPB)[i0],i0);

		if (Intensity(i0)>=0.) {
			nparticles+=Intensity(i0);
		} else {
			cerr<<"\nERROR: Optimizer::WritePencilBeams(): negative PB particle number found: PB: "<<i0<<" Particle Nb: "<<Intensity(i0)<<" !\n"<<endl;
			exit(-1);	// FIXME
		}
	}

	if (!FileOut.is_open()) {
		cerr<<"\nERROR: Optimizer::WritePencilBeams(): I/O error!\n"<<endl;
		exit(-1);	// FIXME
	}

	FileOut<<"*********** Optimized RTPLAN Data for Beam Source Routine ********************\n"
	       <<"* Spot Energy [GeV/u], Momentum spread [GeV/c], Angular divergence [rad],"
	       <<"* Spot width in X dir [cm], Spot width in Y dis [cm], Spot Weight,"
	       <<"* Spot Pos in X dir [cm], Spot Pos in Y dir [cm] ,Spot Pos in Z dir [cm], "
	       <<"* Direction cosine X, Direction cosine Y, Direction cosine Z,"
	       <<"* Xx, Xy, Xz, Atomic Number, Mass Number, Beam ID " <<endl;

	for (unsigned int i0=0;i0<NumPB();i0++) {
		if ( (!isActive(i0))&&(Intensity(i0)!=0.) ) {
			cout<<"\nWARNING: Optimizer::WritePencilBeams(): Something went wrong: PB#"
			    <<i0<<" is inactive but has an intensity of:"<<Intensity(i0)<<"\n"<<endl;
		}
		if ( (isActive(i0))&&(Intensity(i0)!=0.))continue;
		//saving the pencil beams with the optimized intensities

			//* - Spot Energy GeV/u, Momentum spread, Angular Spread,
			//* SpotSigmaX size, SpotSigma Y size, SpotWeight,
			//* xPos, yPos ,zPos,cosx,cosy,cosz, Xx, Xy, Xz,
			//* Atomic Number, Mass Number, BeamID

		FileOut<<"USRICALL  , "<<setw(PR)<<box[i0].Energy*(-1)<<" " //FLUKA requirements
		       <<setw(PR)<<MomSpread(i0)<<" "
		       <<setw(PR)<<AngSpread(i0)<<" "
		       <<setw(PR)<<SigmaX(i0)<<" "
		       <<setw(PR)<<SigmaY(i0)<<" "
		       <<setw(PR)<<Intensity(i0)<<" \n"
		       <<"USRICALL  , "<<setw(PR)<<Position(i0).x<<" "
		       <<setw(PR)<<Position(i0).y<<" "
		       <<setw(PR)<<Position(i0).z<<" "
		       <<setw(PR)<<DirCosines(i0).x<<" "
		       <<setw(PR)<<DirCosines(i0).y<<" "
		       <<setw(PR)<<DirCosines(i0).z<<" &\n"
		       <<"USRICALL  , "<<setw(PR)<<CSCosines(i0).x<<" "
		       <<setw(PR)<<CSCosines(i0).y<<" "
		       <<setw(PR)<<CSCosines(i0).z<<" "
		       <<setw(2)<<Z(i0)<<" "
		       <<setw(2)<<A(i0)<<" "
		       <<setw(1)<<BeamId(i0)<<" &&"<<endl;

	}

	FileOut.close();
	cout<<" DONE"<<endl;
	cout<<" Total number of particles for plan: "<<getTotalParticleNb()<<endl;

	return 0;
} // writePB

int PencilBeam::clean(vector<PENCIL_BEAM> *Struct)
{
	for (unsigned int i0=0; i0<Struct->size(); i0++) {
		(*Struct)[i0].Intensity.clear();
	}

	Struct->clear();

	return 1;
} // clean

/*--------------------------------------------------------*/
int PencilBeam::readNewIntensities(string FileName)
{
	cout<<"Reading New Intensities from File <"<<FileName<<"> ..."<<endl;
	cout.flush();

	PencilBeam tempPB=PencilBeam(FileName);

	int PbNb=tempPB.NumPB();

	if (PbNb!=(int)box.size()) {
		cerr<<"\nERROR: Optimizer::readNewIntensities(): PbNb!=fPencilBeam.size()"<<PbNb<<"!="<<(int)box.size()<<"\n"<<endl;
		exit(-1);	// FIXME
	}

	//FIXME: Search throught Id, not position in the matrix
//	int CountPb=0;
	while(PbNb) {
		PbNb--;
		addIntensity(PbNb,tempPB.Intensity(PbNb));
//		CountPb++;
	}

	cout<<" DONE"<<endl;

	cout<<"### Read new intensities ##########################"<<endl;
	cout<<"  Total number of pencil beams: "<<box.size()<<endl;
	cout<<"  Set "<<PbNb<<" new pencil beam instensities"<<endl;
	cout<<"###################################################"<<endl;

	return 0;
} // readNewIntensities

int PencilBeam::scaleIntensities(double ScaleIntensities)
{
	if (ScaleIntensities<0.) {
		cerr<<"ERROR: Scaling by negative number not permitted: "<<ScaleIntensities<<"!"<<endl;
		exit(-1);	// FIXME
	}

	for (unsigned int i0=box.size();i0--;) {
		addIntensity(i0,box[i0].Intensity.back()*ScaleIntensities);
	}

	cout<<"### Scale intensities #############################"<<endl;
	cout<<"  Scaling by factor "<<ScaleIntensities<<" "<<endl;
	cout<<"###################################################"<<endl;

	return 0;
} // scaleIntensities
