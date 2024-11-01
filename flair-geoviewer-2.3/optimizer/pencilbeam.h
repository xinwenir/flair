// ********************************************
// Author:  wioletta.kozlowska@cern.ch
// based on the T.Boehlen version from 2011
// Modified wioletta.kozlowska@cern.ch
// Version: 3.2
// Last change: 13/11/2016
// *********************************************
#ifndef PENCILBEAM_H
#define PENCILBEAM_H

#include <string>
#include <vector>

#include "vector.h"

typedef struct {
	//Information from RTPlan
	word	BeamId;			// Tag for different beams
	unsigned SpotId;		// Tag for spot id
	word	Z,A;			// Atomic and Mass Number
	Vector	Position;		// x,y,z in the Isocenter from RTPLAN
	Vector	DirCosines;		// cosx, cosy, cosz directional cosines
	Vector	CSCosines;		// Xx, Xy, Xz
	double	InitialIntensity;	// Spot Weight
	double	Energy;			// Spot Energy
	double	MomSpread;		// Momentum Spread
	double	AngSpread;		// Angular Spread
	double	SpotSigma[2];		// FWHM for beam at extraction point

	std::vector<double> Intensity;	// Intensity [number of primaries]
	bool   Active;			// Is the beam Active? -> int>0
	int    MaxDoseVoxelIndex;	// Voxel Index of Max dose
	int    MaxDoseRoiBoxId;		// In which ROI?
	double MaxDosePerPrimary;	// in Gy/primary
} PENCIL_BEAM;

class PencilBeam {
private:
	std::string	FilePB;
	double		TotalParticleNb;

public:
	std::vector<PENCIL_BEAM> box;
	PencilBeam(std::string File);
	PencilBeam() {};
	virtual ~PencilBeam();

	bool readPB          (std::string fileName);
	bool writePB         (std::string fileName, std::vector<PENCIL_BEAM>* vectorPB);

	int clean           (std::vector<PENCIL_BEAM>* vectorPB);
	int readNewIntensities(std::string FileName);
	int scaleIntensities(double ScaleIntensities);

	/*
	 * The following methods specifies the selection of PBs for selective dose calculations
	 */
	bool isPbSelected(unsigned idPB) {
			if (box[idPB].Z==6) {
				return true;
			} else {
				return false;
			}
		};

	//set
	void setPB (PENCIL_BEAM *PB, unsigned idPB) { box[idPB]=*PB; }

	void setFilePB (std::string File) { FilePB = File; }

	void setTotPartNb (double ParticleNb) { TotalParticleNb = ParticleNb; }

	void setIntensity (unsigned idPB, unsigned idInt, double Intensity) {
			box[idPB].Intensity[idInt] = Intensity;	 //try and catch?
		}

	void setMaxDoseVoxelIndex (unsigned idPB, int VxlID) {
			box[idPB].MaxDoseVoxelIndex = VxlID;
		}

	void setMaxDoseRoiBoxId (unsigned idPB, unsigned RoiID) {
			box[idPB].MaxDoseRoiBoxId = RoiID;
		}

	void setMaxDosePerPrimary (unsigned idPB, double DosePP) {
			box[idPB].MaxDosePerPrimary = DosePP;
		}

	void activate (unsigned idPB) { box[idPB].Active = true; }

	inline void deactivate (unsigned idPB) { box[idPB].Active = false; }

	//add
	void addPB (PENCIL_BEAM SinglePB) { box.push_back(SinglePB); }

	void addIntensity (unsigned idPB, double Intensity) {
			box[idPB].Intensity.push_back(Intensity); //try and catch?
		}

	//get
	std::vector<PENCIL_BEAM> *getPBBox() { return &box; }

	PENCIL_BEAM *getSinglePB (unsigned idPB) {
		return &box[idPB]; }

	std::string getFilePB () {
		return FilePB; }

	double getTotalParticleNb () {
		return TotalParticleNb;}

	double Intensity (unsigned idPB) {
			return box[idPB].Intensity.back(); //try and catch
		}

	double Intensity (unsigned idPB, unsigned idInt) {
			return box[idPB].Intensity[idInt]; //try and catch
		}

	double InitIntensity (unsigned idPB) {
			return box[idPB].InitialIntensity;
		}

	bool isActive (unsigned idPB) {
			return box[idPB].Active;
		}

	unsigned NumPB () {
			return box.size();
		}

	unsigned BeamId (unsigned idPB) {
			return box[idPB].BeamId;
		}

	unsigned SpotId (unsigned idPB) {
			return box[idPB].SpotId;
		}

	unsigned A (unsigned idPB) {
			return box[idPB].A;
		}

	unsigned Z (unsigned idPB) {
			return box[idPB].Z;
		}

	Vector Position (unsigned idPB) {
			return box[idPB].Position;
		}

	Vector DirCosines (unsigned idPB) {
			return box[idPB].DirCosines;
		}

	Vector CSCosines (unsigned idPB) {
			return box[idPB].CSCosines;
		}

	double Energy (unsigned idPB) {
			return box[idPB].Energy;
		}

	double MomSpread (unsigned idPB) {
			return box[idPB].MomSpread;
		}

	double AngSpread (unsigned idPB) {
			return box[idPB].AngSpread;
		}

	double SigmaX (unsigned idPB) {
			return box[idPB].SpotSigma[0];
		}

	double SigmaY (unsigned idPB) {
			return box[idPB].SpotSigma[1];
		}

	int MaxDoseVoxelIndex (unsigned idPB) {
			return box[idPB].MaxDoseVoxelIndex;
		}

	int MaxDoseRoiBoxId (unsigned idPB) {
			return box[idPB].MaxDoseRoiBoxId;
		}

	double MaxDosePerPrimary (unsigned idPB) {
			return box[idPB].MaxDosePerPrimary;
		}

	//some more getters are missing
	//clear function in missing
}; // class PencilBeam

#endif
