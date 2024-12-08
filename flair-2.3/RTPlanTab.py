# 

#
# Copyright and User License
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright 2006-2019 CERN and INFN
# 
#
# Please consult the LICENSE file for the license 
#
# DISCLAIMER
# ~~~~~~~~~~
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY, OF
# SATISFACTORY QUALITY, AND FITNESS FOR A PARTICULAR PURPOSE
# OR USE ARE DISCLAIMED. THE COPYRIGHT HOLDERS AND THE
# AUTHORS MAKE NO REPRESENTATION THAT THE SOFTWARE AND
# MODIFICATIONS THEREOF, WILL NOT INFRINGE ANY PATENT,
# COPYRIGHT, TRADE SECRET OR OTHER PROPRIETARY RIGHT.
#
# LIMITATION OF LIABILITY
# ~~~~~~~~~~~~~~~~~~~~~~~
# THE COPYRIGHT HOLDERS AND THE AUTHORS SHALL HAVE NO
# LIABILITY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
# CONSEQUENTIAL, EXEMPLARY, OR PUNITIVE DAMAGES OF ANY
# CHARACTER INCLUDING, WITHOUT LIMITATION, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES, LOSS OF USE, DATA OR PROFITS,
# OR BUSINESS INTERRUPTION, HOWEVER CAUSED AND ON ANY THEORY
# OF CONTRACT, WARRANTY, TORT (INCLUDING NEGLIGENCE), PRODUCT
# LIABILITY OR OTHERWISE, ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGES.
#
# Author: Wioletta.Kozlowska@cern.ch
# Date:	16-Aug-2015



import os
import math
import time
try:
	import pydicom as dicom
except ImportError:
	import dicom
from numpy import mean

import undo
import Dicom
import Input
import Ribbon
import tkFlair
import tkExtra
import FlairRibbon
import bFileDialog

from log import say

try:
	from tkinter import *
except ImportError:
	from tkinter import *

import math

#maximum energy for tabulation in the BEAM card
MAX_EN = 0.3

#===============================================================================
# RT plan processing
#===============================================================================
class RTPlan:
	#-----------------------------------------------------------------------
	def __init__(self, dataset):
		#Patient Info
		self.patientName = self.patientAge = self.patientSex = "No information"

		#RTPlan info
		self.planName = self.planDate = self.doseType = self.planIntent = \
			self.appStatus = self.reviewDate = self.reviewerName = \
			self.studyDate = "No information"
		self.patientSetup = "HFS"

		#Beam Info
		self.bidx	 = {}	# index to beam name#
		self.bName	 = {}	# Beam Name to idx
		self.bWeight	 = {}	# Beam Weight
		self.bType	 = {}	# Planned Beam Type
		self.radType	 = {}	# Radiation Type. Description of ions not supported
		self.sMode	 = {}	# Scanning Mode
		self.numCP	 = {}	# Number of Control Points
		self.vSAxis	 = {}	# Virtual Source-Axis Distance
		self.totWeight	 = {}	# Final Cumulative Meterset Weight
		self.unit	 = {}	# Unit Number of Particles
		self.isocenter	 = {}

		self.gantryAngle    = {} # Gantry Angle		 - No Gantry Rotation Supported
		self.patientAngle   = {} # Patient Angle	 - No Patient Rotation Supported
		self.gantryPAngle   = {} # Gantry Pitch Angle	 - No Gantry Rotation Supported
		self.tableTPAngle   = {} # Table Top Pitch Angle - No Table Rotation Supported
		self.tableTRAngle   = {} # Table Top Roll Angle  - No Table Rotation Supported

		# no beam limiting devices
		if isinstance(dataset,str):
			try:
				dataset = dicom.read_file(dataset, force=True)
				self.dataset = dataset
			except:
				return
		self.fillRTPlanData(dataset)

	#-----------------------------------------------------------------------
	def fillRTPlanData(self,dataset):
		""" Filling information from RTPlan"""
		try:
			self.patientName = str(dataset.PatientsName.family_name)+", " +str(dataset.PatientsName.given_name)
		except AttributeError:
			self.patientName = "No information"

		try:
			self.patientAge = str(dataset.PatientsAge)
		except AttributeError:
			self.patientAge = "No information"

		try:
			self.patientSex = str(dataset.PatientsSex)
		except AttributeError:
			self.patientSex = "No information"

		try:
			self.planName = str(dataset.RTPlanName)
		except AttributeError:
			self.planName = "No information"

		try:
			RTdate = dataset.RTPlanDate
			RTtime = dataset.RTPlanTime
			self.planDate = "%s.%s.%s %s:%s" % \
				(RTdate[: 4], RTdate[ 4: 6], RTdate[ 6:],
				 RTtime[:-4], RTtime[-4:-2])
		except AttributeError:
			self.planDate = "No information"

		try:
			self.doseType = str(dataset.DoseType)
		except AttributeError:
			self.doseType = "No information"

		try:
			self.planIntent = str(dataset.PlanIntent)
		except AttributeError:
			self.planIntent = "No information"

		try:
			self.appStatus = str(dataset.ApprovalStatus)
		except AttributeError:
			self.appStatus = "No information"

		try:
			Revdate = dataset.ReviewDate
			Revtime = dataset.ReviewTime
			self.reviewDate = "%s.%s.%s %s:%s" % \
				(Revdate[: 4], Revdate[ 4: 6], Revdate[ 6:],
				 Revtime[:-4], Revtime[-4:-2])
		except AttributeError:
			self.reviewDate = "No information"

		try:
			self.reviewerName = str(dataset.ReviewerName.family_name)+", " +str(dataset.ReviewerName.given_name)
		except AttributeError:
			self.reviewerName = "No information"

		try:
			Studydate = dataset.StudyDate
			Studytime = dataset.StudyTime
			self.studyDate = "%s.%s.%s %s:%s" % \
				(Studydate[: 4], Studydate[ 4: 6], Studydate[ 6:],
				 Studytime[:-4], Studytime[-4:-2])
		except AttributeError:
			self.studyDate = "No information"

		try:
			self.beamNum = len(dataset.IonBeamSequence)
		except AttributeError:
			self.beamNum = 0

		try:
			self.patientSetup = str(dataset.PatientSetupSequence.PatientPosition)
		except AttributeError:
			self.patientSetup = "HFS"

		for beamId in range(self.beamNum):
			try:
				beam = dataset.IonBeamSequence[beamId]
			except AttributeError:
				return

			try:
				self.bidx[beamId] = int(beam.BeamNumber)
			except AttributeError:
				self.bidx[beamId] = "No information"

			try:
				self.numCP[beamId] = int(beam.NumberofControlPoints)
			except AttributeError:
				self.numCP[beamId] = 0

			try:
				self.bName[beamId] = str(beam.BeamName)
			except AttributeError:
				self.bName[beamId] = "No information"

			try:
				self.bWeight[beamId] = float(beam.FinalCumulativeMetersetWeight)
			except AttributeError:
				self.bWeight[beamId] = 0

			try:
				self.bType[beamId] = str(beam.BeamType)
			except AttributeError:
				self.bType[beamId] = "No information"

			try:
				self.radType[beamId] = str(beam.RadiationType)
			except AttributeError:
				self.radType[beamId] = "No information"

			try:
				self.sMode[beamId] = str(beam.ScanMode)
			except AttributeError:
				self.sMode[beamId] = "No information"

			try:
				icps = beam.IonControlPointSequence[0]
			except:
				return

			try:
				self.gantryAngle[beamId] =  self.correctAng(icps.GantryAngle)
			except AttributeError:
				self.gantryAngle[beamId] = 0

			try:
				self.gantryPAngle[beamId] = self.correctAng(icps.GantryPitchAngle)
			except AttributeError:
				self.gantryPAngle[beamId] = 0

			try:
				self.patientAngle[beamId] = self.correctAng(icps.PatientSupportAngle)
			except AttributeError:
				self.patientAngle[beamId] = 0

			try:
				self.tableTPAngle[beamId] = self.correctAng(icps.TableTopPitchAngle)
			except AttributeError:
				self.tableTPAngle[beamId] = 0

			try:
				self.tableTRAngle[beamId] = self.correctAng(icps.TableTopRollAngle)
			except AttributeError:
				self.tableTRAngle[beamId] = 0

			try:
				self.isocenter[beamId] = icps.IsocenterPosition
			except AttributeError:
				self.isocenter[beamId] = [0,0,0]

	#		try:
	#			self.refPatientSetup[int(beam.BeamNumber)] = beam.ReferencedPatientSetupNumber
	#		except AttributeError:
	#			self.refPatientSetup[int(beam.BeamNumber)] = 0

			try:
				self.vSAxis[beamId] = beam.VirtualSourceAxisDistances
			except AttributeError:
				self.vSAxis[beamId] = "No information"

			try:
				self.totWeight[beamId] = beam.FinalCumulativeMetersetWeight
			except AttributeError:
				self.totWeight[beamId] = "No information"
			try:
				self.unit[beamId] = str(beam.PrimaryDosimeterUnit)
			except AttributeError:
				self.unit[beamId] = "None"

		#no beam limiting devices

	#-----------------------------------------------------------------------
	def correctAng(self, ang):
		""" Correcting angle values, if angle do not exist return 0"""
		if ang == [] or ang == 0.0:
			return 0
		return ang

	#-----------------------------------------------------------------------
	def beamData(self, dataset, beamID, dMomentum=0, dAngular=0, printed=0, dFWHMx=0, dFWHMy=0):
		""" Exporting Beam parameters to the consisent table"""

		try:
			IonSeq=dataset.IonBeamSequence[beamID]
			IonCP=IonSeq.IonControlPointSequence
		except:
			say("Beam File Error",
				"No Beam Data File was created")
			return None

		unit = self.unit[beamID]

		# FIXME: If there is no NP adapt the table with MU
		if unit!= 'NP':
			say("Beam File Error",
				"Primary Dosimeter Unit is not the Number of Particles!!!")
			return None

		#Radiation type
		RadType  = self.radType[beamID]
		if RadType == "PROTON":
			Particle="PROTON"
		elif RadType == "ION":
			try:
				AtNb = int(IonSeq.RadiationAtomicNumber)
				MsNb = int(IonSeq.RadiationMassNumber)
				Particle == -(MsNb*100000+AtNb*100)
			except AttributeError:
				say("Beam File Error",
				"No information on Ion Type")
				return None
		else:
			say("Beam File Error",
			    "No information on Ion Type")
			AtNb, MsNb = 1, 1
			Particle="PROTON"

		NumCP	  = self.numCP[beamID]
		VSAxis	  = self.vSAxis[beamID]
		TotWeight = self.totWeight[beamID]

		GAngle	  = self.gantryAngle[beamID]
		GPAngle   = self.gantryPAngle[beamID]
		PAngle	  = self.patientAngle[beamID]
		TTPAngle  = self.tableTPAngle[beamID]
		TTRAngle  = self.tableTRAngle[beamID]

		PartWeight = 0
		BeamMet   = 1
		CPTab	  = []

		try:
			BeamMet = dataset.FractionGroupSequence[0].ReferencedBeamSequence[beamID].BeamMeterset
		except AttributeError:
			BeamMet = 1

		PartWeight = BeamMet/TotWeight

		#fill the table with control point info
		SumWeight = 0 # Control sum weight

		for i in range(NumCP):
			xPos,yPos	 = 0, 0
			SSWeight, NumSSP = 0, 0		# Number of Scan Spot Positions
			try:
				CP = IonSeq.IonControlPointSequence[i]

				#not really the thruth.. because this is the spot size at isocenter
				if (not dFWHMx):
					xSSSize = - CP.ScanningSpotSize[0]/10
				else:
					xSSSize = dFWHMx

				if (not dFWHMy):
					ySSSize = - CP.ScanningSpotSize[0]/10
				else:
					ySSSize = dFWHMy

				BEnergy = float(CP.NominalBeamEnergy)/1000 # Nominal Beam Energy in GeV/n
				xPos=[x/10 for x in CP.ScanSpotPositionMap[0::2]] # Scan Pos x in mm
				yPos=[x/10 for x in CP.ScanSpotPositionMap[1::2]] # Scan Pos y in mm
				zPos=mean(VSAxis)/10
				NumSSP = CP.NumberOfScanSpotPositions
				if NumSSP ==1:
					SSweight = [CP.ScanSpotMetersetWeights]
				else:
					SSweight = CP.ScanSpotMetersetWeights
			except AttributeError:
				say("Beam File Error", "No Control Points Exist")
				return None

			# check is this scanning spot is correct
			if abs(sum(SSweight))<0.001:
				continue

			SumWeight += abs(sum(SSweight)) #Control sum weight

			#Solve all angles
			f1 = lambda a,b:a*math.cos(math.radians(ang))-math.sin(math.radians(ang))*b
			f2 = lambda a,b:a*math.cos(math.radians(ang))+math.sin(math.radians(ang))*b
			for spot in range(NumSSP):

				posx, posy = 0,0
				posz = abs(zPos)

				# Gantry Angle (X,Z)
				ang  = GAngle
				cosx = f2(xPos[spot],-abs(zPos))
				cosy = yPos[spot]
				cosz = f1(-abs(zPos),xPos[spot])
				posx,posz = f2(posx,posz),f1(posz,posx)

				#Gantry Pitch Angle (Y,Z)
				ang  = GPAngle
				cosy,cosz = f1(cosy,cosz),f2(cosz,cosy)
				posy,posz = f1(posy,posz),f2(posz,posy)

				#Patient Angle (X,Y)
				ang  = PAngle
				cosx,cosy = f2(cosx,cosy),f1(cosy,cosx)
				posx,posy = f2(posx,posy),f1(posy,posx)

				#Table Top Pitch Angle (Y,Z)

				ang  = TTPAngle
				cosy,cosz = f2(cosy,cosz),f1(cosz,cosy)
				posy,posz = f2(posy,posz),f1(posz,posy)

				#Table Top Roll Angle (X,Z)
				ang  = TTRAngle
				cosx,cosz = f1(cosx,cosz),f2(cosz,cosx)
				posx,posz = f1(posx,posz),f2(posz,posx)

				r = math.sqrt(math.pow(cosx,2)+math.pow(cosy,2)+math.pow(cosz,2))

				cosx = cosx/r
				cosy = cosy/r
				cosz = cosz/r

				# TODO: give the user chance to change it?
				Xx, Xy, Xz = 0, 0, 1

			# fill in the Control Point Table
				CPTab.append([Particle,-BEnergy, SSweight[spot]*PartWeight, dMomentum, dAngular,printed,
					      posx, posy, posz, xSSSize, ySSSize,
					      cosx, cosy, cosz, Xx, Xy, Xz,
					      self.bidx[beamID]])

		if ((SumWeight-TotWeight)/SumWeight>0.0001):
			say("Beam File Error", "Beam Meterset Weight not correct")
		#	return

		return CPTab

	#-----------------------------------------------------------------------
	def exportData(self, CPTab, filename, bidx ):
                """Exporting data to the external file"""

#		out_fname =""
		# Output file
                try:
                        out_fname = "RTPLANbeam_%s_%d.inp" % (filename,bidx)
                        out = open(out_fname,"w+")
                except:
                        say("Beam File Error",
				"No Beam Data File was created")
                        return

		# Print to output file
                try:
                        print("********************* RTPLAN Data for Beam Spot Cards ************************",file=out)
                        print("*Card Spotbeam",file=out)
                        print("* 1)particle id 2)beam laboratory kinetic energy/n [GeV/n], 3)beam spot weight,",file=out)
                        print("* 4)momentum FWHM spread (GeV/c), 5)beam angular divergence (mrad), 6)print beam spot info,",file=out)
                        print("*Card Spotpos",file=out)
                        print("* 1)x-position of the beam spot, 2)y-position of the beam spot, 3)z-position of the beam spot, ",file=out)
                        print("* 4)gaussian beam FWHM with in x-direction (cm), 5)gaussian beam FWHM width in y-direction (cm),",file=out)
                        print("*Card Spotdir",file=out)
                        print("* 1)x-cosine of the beam spot, 2)y-cosine of the beam spot, 3)z-cosine of the beam spot,",file=out)
                        print("* 4)x-cosine of ht beam x-axis, 5)y-cosine of the beam x-axis, 6)z-cosine of the beam x-axis ",file=out)

                        beamidx=1

                        for spot in CPTab:
                                print("SPOTBEAM",", ".join(map(str,spot[0:6]  )),   int(beamidx), sep=", ", file=out)
                                print("SPOTPOS", ", ".join(map(str,spot[6:11] )),"",int(beamidx), sep=", ", file=out)
                                print("SPOTDIR", ", ".join(map(str,spot[11:17])),   int(beamidx), sep=", ", file=out)
                                beamidx+=1
                        out.close()
                except:                        
                        say("Beam File Error", "No data written to file")
                        return

                return out_fname

	#-----------------------------------------------------------------------
	def beams2Files(self, dataset, filename, dMom=0, dAng=0, dPrinted=0, dFWHMx=0, dFWHMy=0):
                """Arranging the beam to correct files"""

                data  =[""]
                files =[]
                spots =[]

		#create beam TP files
                for idx in range(len(self.bidx)):
                        data.append(self.beamData(dataset,idx,dMom,dAng,dPrinted,dFWHMx,dFWHMy))

                data[0]  = [item for sublist in data for item in sublist]

                for idx in range(len(data)):
                        files.append(self.exportData(data[idx], filename, idx))
                        spots.append(len(data[idx]))

                return files,spots

#===============================================================================
# RTPlan Tab class
#===============================================================================
class RTPlanTab(FlairRibbon.FlairTab):
	#-----------------------------------------------------------------------
	def __init__(self, master, page, **kw):
		FlairRibbon.FlairTab.__init__(self, master, page, **kw)

		self.filesel = None
		self.beamsel = None

		# ------------------ Plan files ---------------
		pane = PanedWindow(self, orient=HORIZONTAL, opaqueresize=1)
		pane.pack(fill=BOTH, expand=YES)

		filesframe =tkExtra.LabelFrame(pane, text = "Choose RT file")
		pane.add(filesframe, stretch="always")

		self.rtList = tkExtra.MultiListbox(filesframe,
						(("RTPlan", 20, None),
						 ("Date",   10, None)))
		self.rtList.pack(anchor=NW, side=TOP, fill=BOTH, expand=YES)
		self.rtList.bind("<Return>",	      self.execute)
		self.rtList.bind("<KP_Enter>",	      self.execute)
		self.rtList.bind("<F2>",	      self.execute)
		self.rtList.bind("<<ListboxSelect>>", self.fileSelect)
		self.rtList.sortAssist = None

		beamframe =tkExtra.LabelFrame(pane, text = "Choose Beam Number")
		pane.add(beamframe, stretch="always")

		self.beamList = tkExtra.MultiListbox(beamframe,
						(("Beam",	10, None),
						 ("# points",	10, None),
						 ("Gantry Ang", 10, None),
						 ("Patient Ang",10, None)))
		self.beamList.pack(anchor=NW, side=TOP, fill=BOTH, expand =Y)
		self.beamList.bind("<<ListboxSelect>>", self.beamSelect)
		self.beamList.sortAssist = None


		#-------------------------Information frame ---------------------
		inputFrame =tkExtra.LabelFrame(pane, text = "Information")
		pane.add(inputFrame, stretch="always")

		rowL, colL = 0,0

		infoframe=tkExtra.LabelFrame(inputFrame, text="RTPLan Information")
		infoframe.grid(row=rowL, column=colL,sticky=EW)
		# ---
		row = 0
		Label(infoframe,text="Patient Info:", foreground="DarkRed").grid(row=row, column=0, sticky=W)

		row += 1
		Label(infoframe, text="Patient's Name:").grid(row=row, column=0, sticky=E)
		self.patname = Label(infoframe, justify=LEFT, anchor=W,
				width=32,
			        foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.patname.grid(row=row, column=1, sticky=EW)

		row += 1
		Label(infoframe, text="Patient's Age:").grid(row=row, column=0, sticky=E)
		self.patage = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.patage.grid(row=row, column=1, sticky=EW)

		row += 1
		Label(infoframe, text="Patient's Sex:").grid(row=row, column=0, sticky=E)
		self.patsex = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.patsex.grid(row=row, column=1, sticky=EW)

		# ---
		row += 5
		Label(infoframe,text="RT Plan Info:", foreground="DarkRed").grid(row=row, column=0, sticky=W)

		row +=1
		Label(infoframe, text="RT Plan Name:").grid(row=row, column=0, sticky=E)
		self.planname = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.planname.grid(row=row, column=1, sticky=EW)

		row += 1
		Label(infoframe, text="Dose Type:").grid(row=row, column=0, sticky=E)
		self.dosetype = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.dosetype.grid(row=row, column=1, sticky=EW)

		row += 1
		Label(infoframe, text="Plan Intent:").grid(row=row, column=0, sticky=E)
		self.planintent = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.planintent.grid(row=row, column=1, sticky=EW)

		row += 1
		Label(infoframe, text="Approval Status:").grid(row=row, column=0, sticky=E)
		self.appstatus = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.appstatus.grid(row=row, column=1, sticky=EW)

		row += 1
		Label(infoframe, text="Review Date").grid(row=row, column=0, sticky=E)
		self.reviewdate = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.reviewdate.grid(row=row, column=1, sticky=EW)

		row += 1
		Label(infoframe, text="Reviewer Name").grid(row=row, column=0, sticky=E)
		self.reviewername = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.reviewername.grid(row=row, column=1, sticky=EW)

		# ---
		row += 5
		Label(infoframe,text="General Beam Info:", foreground="DarkRed").grid(row=row, column=0, sticky=W)

		row +=1
		Label(infoframe, text="Beam Name:").grid(row=row, column=0, sticky=E)
		self.beamname = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.beamname.grid(row=row, column=1, sticky=EW)

		row +=1
		Label(infoframe, text="Beam Weight:").grid(row=row, column=0, sticky=E)
		self.beamdescr = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.beamdescr.grid(row=row, column=1, sticky=EW)

		row +=1
		Label(infoframe, text="Beam Type:").grid(row=row, column=0, sticky=E)
		self.beamtype = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.beamtype.grid(row=row, column=1, sticky=EW)

		row +=1
		Label(infoframe, text="Radiation Type:").grid(row=row, column=0, sticky=E)
		self.radtype = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.radtype.grid(row=row, column=1, sticky=EW)

		row +=1
		Label(infoframe, text="Scan Mode:").grid(row=row, column=0, sticky=E)
		self.scanmode = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.scanmode.grid(row=row, column=1, sticky=EW)

		row +=1
		Label(infoframe, text="# Control Points:").grid(row=row, column=0, sticky=E)
		self.numcp = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.numcp.grid(row=row, column=1, sticky=EW)

		# ---
		row += 5
		Label(infoframe,text="Rotations Info:", foreground="DarkRed").grid(row=row, column=0, sticky=W)

		row +=1
		Label(infoframe, text="Gantry Angle:").grid(row=row, column=0, sticky=E)
		self.gangle = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.gangle.grid(row=row, column=1, sticky=EW)

		row +=1
		Label(infoframe, text="Gantry Pitch Angle:").grid(row=row, column=0, sticky=E)
		self.gpangle = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.gpangle.grid(row=row, column=1, sticky=EW)

		row +=1
		Label(infoframe, text="Patient Support Angle:").grid(row=row, column=0, sticky=E)
		self.pangle = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.pangle.grid(row=row, column=1, sticky=EW)

		row +=1
		Label(infoframe, text="Table Top Pitch Angle:").grid(row=row, column=0, sticky=E)
		self.ttpangle = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.ttpangle.grid(row=row, column=1, sticky=EW)

		row +=1
		Label(infoframe, text="Table Roll Angle:").grid(row=row, column=0, sticky=E)
		self.ttrangle = Label(infoframe, justify=LEFT, anchor=W,
				foreground=tkFlair._ILABEL_FOREGROUND_COLOR,
				background=tkFlair._ILABEL_BACKGROUND_COLOR)
		self.ttrangle.grid(row=row, column=1, sticky=EW)

		infoframe.grid_columnconfigure(1, weight=1)

		rowL+=1

		#--------------------------Additional Info to fill ---------------
		fillframe=tkExtra.LabelFrame(inputFrame, text="Additional Information")
		fillframe.grid(row=rowL, column=colL,sticky=EW)

		# ------------------Additional data for SPOTBEAM -----------------

		#TODO: Save the latest values!
		#FIXME: This is not very useful for all beam energies, how to apply as a % or what?
		col, row = 0,0
		Label(fillframe, text = "Momentum spread [GeV/c]:",).grid(
				row=row, column=col, sticky=E)

		col += 1
		self.momentum = tkExtra.FloatEntry(fillframe, width=10, background="White")
		self.momentum.insert(END,0)
		self.momentum.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.momentum, "Additional momentum spread GeV/c")

		col += 1
		self.momentumGauss=IntVar()
		c1=Checkbutton(fillframe, text="Gaussian spread", variable=self.momentumGauss
				  , onvalue=TRUE )
		c1.grid(row=row, column=col, sticky=W)
		c1.select()

		col = 0
		row +=1
		Label(fillframe, text = "Angular divergence [mrad]:").grid(
				row=row, column=col, sticky=E)
		col +=1
		self.angular = tkExtra.FloatEntry(fillframe, width=10, background="White")
		self.angular.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.angular, "Beam angular divergence [mrad]")

		col += 1
		self.angularGauss=IntVar()
		c2=Checkbutton(fillframe, text="Gaussian spread", variable=self.angularGauss
				  , onvalue=TRUE )
		c2.grid(row=row, column=col, sticky=W)
		c2.select()

		# ------------------Additional data for SPOTPOS -----------------

		col =0
		row += 1
		Label(fillframe, text = "Beam FWHM in x-direction[cm]:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.fwhmx = tkExtra.FloatEntry(fillframe, width=10, background="White")
		self.fwhmx.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.fwhmx, "Beam FWHM(gaussian beam)/width(rectangular beam) in x direction from accelerator in [cm]")

		col += 1
		self.xGauss=IntVar()
		c1=Checkbutton(fillframe, text="Gaussian spread", variable=self.xGauss
				  , onvalue=TRUE)
		c1.grid(row=row, column=col, sticky=W)
		c1.select()

		col = 0
		row += 1
		Label(fillframe, text = "Beam FWHM in y-direction[cm]:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.fwhmy = tkExtra.FloatEntry(fillframe, width=10, background="White")
		self.fwhmy.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.fwhmy, "Beam FWHM(gaussian beam)/width(rectangular beam) in y direction from accelerator in [cm]")

		col +=1
		self.yGauss=IntVar()
		c2=Checkbutton(fillframe, text="Gaussian spread", variable=self.yGauss
				  , onvalue=TRUE)
		c2.grid(row=row, column=col, sticky=W)
		c2.select()

	#-----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairTab.createRibbon(self)

		# --------------Choose correct RTDOSE-----------------
		groupDose = Ribbon.LabelGroup(self.ribbon, "Choose RTDose")
		groupDose.label["background"] = Ribbon._BACKGROUND_GROUP2
		groupDose.pack(side=LEFT, fill=Y, padx=0, pady=0)

		col, row = 0,1
		Label(groupDose.frame, text="RTDOSE:",
			 background=Ribbon._BACKGROUND).grid(
				 row=row, column=col, sticky=EW)
		row += 1
		self.rtdose = tkExtra.Combobox(groupDose.frame, width=30, foreground="DarkBlue")
		self.rtdose.grid(row=row, column=col, sticky=EW)
		tkExtra.Balloon.set(self.rtdose, "Choose RTDOSE input file")

		row +=1
		self.printed=IntVar()
		c1=Checkbutton(groupDose.frame, text="Print info", variable=self.printed
				  , onvalue=TRUE)
		c1.grid(row=row, column=col, sticky=W)
		tkExtra.Balloon.set(c1, "Print all beam spot info")

		# -------------------Process RTPlan -------------------
		group = Ribbon.LabelGroup(self.ribbon, "Process RTPlan")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		col, row = 0,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["gear32"],
				text="Create",
				compound=TOP,
				command=self.execute,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, sticky=E)
		tkExtra.Balloon.set(b, "Create RTPLAN input file")

		# ---
#				col+=1
#		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["clean32"],
#				text="Clean",
#				compound=TOP,
#				command=self.clear,
#				background=Ribbon._BACKGROUND)
#		b.grid(row=row, column=col, sticky=E)
#		tkExtra.Balloon.set(b, "Remove all generated files")

		self._updated = False

	# ----------------------------------------------------------------------
	def dicom(self):	return self.page.dicom

	# ----------------------------------------------------------------------
	def refresh(self):
		di = self.page.dicominfo
		if di is None:
			self.clear()
			return

		self.dataset = None

		self.fwhmx.set("text")
		self.fwhmy.set("text")
		self.momentum.set("text")
		self.angular.set("text")
		self.filesel = 0
		self.beamsel = 0
		self.rtdose.set(di.get("rtdose","").strip())

		# Fill a list of RTDOSE
		self.rtdose.delete(0,END)
		self.rtdose.insert(END,"")
		for di in self.project.dicoms:
				if di.modality=="RTDOSE":
					self.rtdose.insert(END,di.name)

	# ----------------------------------------------------------------------
	def activate(self):
		di = self.page.selectedDicom()

		if di is None:
			self.clear()
			return
		try:
			fn = os.path.join(self.page.dicominfo["directory"], di.files[0][0])
			self.dataset = dicom.read_file(fn,
						defer_size=256,
						force=True)
		except dicom.filereader.InvalidDicomError:
			self.dataset = None
			return

		if self.dataset:
			self.dicom().init(self.dataset)

		if self.dicom().modality == "RTPLAN":
			self.filesel = 0
			self.beamsel = 0
			self.rtList.labels(["RTPlan","Date"])
			self.updateRTPlan(di.files)

			self.rtdose.delete(0,END)
			self.rtdose.insert(END,"")
			for di in self.project.dicoms:
				if di.modality=="RTDOSE":
					self.rtdose.insert(END,di.name)
			return
		else:
			self.flair.notify("No DICOM file selected",
				"Please select RTPlan file first")
			self.clear()
			return

	# ----------------------------------------------------------------------
	def clear(self):
		""" Clear table info"""
		self.dataset=None
		self.rtList.delete(0,END)
		self.beamList.delete(0,END)
		self.rtdose.delete(0,END)
		self.patname["text"]	 = ""
		self.patage["text"]	 = ""
		self.patsex["text"]	 = ""
		self.planname["text"]	 = ""
		self.dosetype["text"]	 = ""
		self.planintent["text"]  = ""
		self.appstatus["text"]	 = ""
		self.reviewdate["text"]  = ""
		self.reviewername["text"]= ""
		self.beamname["text"]	 = ""
		self.beamdescr["text"]	 = ""
		self.beamtype["text"]	 = ""
		self.radtype["text"]	 = ""
		self.numcp["text"]	 = ""
		self.gangle["text"]	 = ""
		self.gpangle["text"]	 = ""
		self.pangle["text"]	 = ""
		self.ttpangle["text"]	 = ""
		self.ttrangle["text"]	 = ""
		self.scanmode["text"]	 = ""

	#-----------------------------------------------------------------------
	def fileSelect(self, event=None):
		"""Update file selection"""

		self.filesel = self.rtList.curselection()[0]

		fn= os.path.join(self.page.dicominfo["directory"],self.rtList.get(self.filesel)[0])
		self.dataset = dicom.read_file(fn,
						defer_size=256,
						force=True)
		self.rtplan = RTPlan(self.dataset)

		self.updateBeam()

	#-----------------------------------------------------------------------
	def beamSelect(self, event=None):
		"""Update beam selection"""

		self.beamsel = self.beamList.curselection()[0]
		self.updateInfo()

	#-----------------------------------------------------------------------
	def updateRTPlan(self,files):
		"""Update the RTPlan Info file"""

		self.clear()

		for fnz in files:
			try:
				fn = os.path.join(self.page.dicominfo["directory"], fnz[0])
				self.dataset = dicom.read_file(fn,
						defer_size=256,
						force=True)
			except:
				self.dataset = None
				return

			if self.dataset:
				self.rtplan = RTPlan(self.dataset)
				self.rtList.insert(END,(fnz[0],self.rtplan.planDate))

		self.rtList.selection_set(self.filesel)
		fn = os.path.join(self.page.dicominfo["directory"],self.rtList.get(self.filesel)[0])

		self.dataset = dicom.read_file(fn,
						defer_size=256,
						force=True)
		self.rtplan = RTPlan(self.dataset)

		self.updateBeam()

	#-----------------------------------------------------------------------
	def updateBeam(self):
		"""Update RTPlan info in the table"""

		self.beamList.delete(0,END)
		self.patname["text"]	 = self.rtplan.patientName
		self.patage["text"]	 = self.rtplan.patientAge
		self.patsex["text"]	 = self.rtplan.patientSex
		self.planname["text"]	 = self.rtplan.planName
		self.dosetype["text"]	 = self.rtplan.doseType
		self.planintent["text"]  = self.rtplan.planIntent
		self.appstatus["text"]	 = self.rtplan.appStatus
		self.reviewdate["text"]  = self.rtplan.reviewDate
		self.reviewername["text"]= self.rtplan.reviewerName

		for beam in list(self.rtplan.bidx):
			self.beamList.insert(END,(self.rtplan.bidx[beam],
						  self.rtplan.numCP[beam],
						  self.rtplan.gantryAngle[beam],
						  self.rtplan.patientAngle[beam]))
		if not self.beamsel:
			self.beamsel = 0
		self.beamList.selection_set(self.beamsel)
		self.updateInfo()

	#-----------------------------------------------------------------------
	def updateInfo(self):
		"""Update Beam info in the table"""

		self.beamname["text"]  = self.rtplan.bName[int(self.beamsel)]
		self.beamdescr["text"] = self.rtplan.bWeight[int(self.beamsel)]
		self.beamtype["text"]  = self.rtplan.bType[int(self.beamsel)]
		self.radtype["text"]   = self.rtplan.radType[int(self.beamsel)]
		self.numcp["text"]	 = self.rtplan.numCP[int(self.beamsel)]
		self.gangle["text"]	= self.rtplan.gantryAngle[int(self.beamsel)]
		self.gpangle["text"]   = self.rtplan.gantryPAngle[int(self.beamsel)]
		self.pangle["text"]	= self.rtplan.patientAngle[int(self.beamsel)]
		self.ttpangle["text"]  = self.rtplan.tableTPAngle[int(self.beamsel)]
		self.ttrangle["text"]  = self.rtplan.tableTRAngle[int(self.beamsel)]
		self.scanmode["text"]  = self.rtplan.sMode[int(self.beamsel)]

	#-----------------------------------------------------------------------
	def findRTDose(self):
		"""Find and open RTSTRUCT dataset"""

		rtuid = self.rtdose.get()
		if rtuid:
			for di in self.project.dicoms:
				if di["uid"] == rtuid:
					try:
						fn = os.path.join(di["directory"],
								  di.files[0][0])
						self.dose=dicom.read_file(fn,defer_size=256,
									  force=True)
					except:
						self.flair.notify("RTDose File",
						  "Problem with opening RTDose File")
						self.dose = None
					break
		else:
			try:
				self.scan4ReferencedDose()
				fn = os.path.join(self.doseInfo["directory"],
						  self.doseInfo.files[0][0])
				self.dose = dicom.read_file(fn, defer_size=256,
								force=True)
			except :
				self.flair.notify("RTDose File",
						  "Corresponding RTDose File do not exist")
				self.dose = None


	#-----------------------------------------------------------------------
	# Processing of input, creating files
	#-----------------------------------------------------------------------

	def execute(self):
		"""Add all required information into the input file"""

		inp=self.flair.project.input
		undoinfo = self.flair.undolistRefresh()
		self.page.get()

		self.findRTDose()

		self._inccards = []
		self.makePatRD(self.rtplan, undoinfo, inp)
		self.makeUSRBIN(self.rtplan, undoinfo, inp)
		self.makeBeamFile(self.rtplan, undoinfo, inp)
		self.makeDef(self.rtplan, undoinfo, inp)
		self.expandIncludeFiles(self._inccards)

		self.flair.addUndolistRefresh(undoinfo,"Update input file from RTPLAN")


		self.flair.notify("RTPLAN",
				"Input File was modified:\n"\
				"VOXEL card modified\n" \
				"USRBIN scoring is added\n"\
				"SOURCE added with RTPlan information\n" \
				"BEAM card modified")
		self.flair.refresh("card")
		return True


	#-----------------------------------------------------------------------
	def makeDef(self,rtplan,undoinfo, inp):
		"""Update the beam numbers define the preprocessor cards"""

		pos = inp.bestPosition("GLOBAL")
		undoinfo.append(self.flair.addCardUndo(
				Input.Card("#define", ["full_treatment"]),pos,False))
		for i in range(rtplan.beamNum):
			pos+=1
			card=	Input.Card("#define", ["only_beam"+str(i+1)])
			undoinfo.append(self.flair.addCardUndo(card,pos,False))
			card.disable()

	#-----------------------------------------------------------------------
	def makePatRD(self,rtplan, undoinfo,inp):
		"""Update the info of the patient position"""

		pos = inp.bestPosition("ROT-DEFI")

		#defining the isocenter
		isocenter=[]

		for i in range(rtplan.beamNum):
			if isocenter != [rtplan.isocenter[i]]:
				isocenter.append(rtplan.isocenter[i])

		for i in range(len(isocenter)): # assuming that is only one
			pos=self._isocenter(i,isocenter[i],undoinfo,pos)
			if ((i>0) and i==len(isocenter)-1):
				undoinfo.append(self.flair.addCardUndo(
					Input.Card("#endif",[]),pos,False))
				pos+=1

		#defining patient setup - only one available!
		setup = [rtplan.patientSetup]
		pos  = self._patposition(0,setup[0],undoinfo,pos)

		#For FLUKA release- - 2011.1c.5---------------------------------
		#undoinfo.append(self.flair.addCardUndo(Input.Card("ROT-DEFI", ["nDICOM", 100, 0, 0,
		#		(isocenter[0][0]/10), (isocenter[0][1]/10), (isocenter[0][2]/10)],
		#		"Isocenter translation"),pos, False))
		#pos+=1
		#---------------------------------------------------------------

		#modifying Voxel
		card = inp["VOXELS"]

		if card:
			card = card[0]
			undoinfo.append(self.flair.setWhatUndo(card,4,"DICOM"))
			undoinfo.append(self.flair.setCommentUndo(card,"Voxel with DICOM->RT rotation"))
		else:
			self.flair.notify("No Voxel Card",
			"Add Voxel Card from Voxel Tab first")

	#-----------------------------------------------------------------------
	def makeUSRBIN(self,rtplan, undoinfo, inp):
		"""Create the USRBIN files data"""

		card = inp["VOXELS"]
		if card:
			card = card[0]
			voxel= Input.Voxel(str(card.what(0))+".vxl")

			xmin = float(card.what(1))
			ymin = float(card.what(2))
			zmin = float(card.what(3))

			if self.dose is not None:

				try:
					imSpac = [float(self.dose.PixelSpacing[0]/10),
						  float(self.dose.PixelSpacing[1]/10),
						  float(self.dose.GridFrameOffsetVector[1]/10)]
					imNum  = [int(self.dose.Columns),
						  int(self.dose.Rows),
						  int(len(self.dose.GridFrameOffsetVector))]
					imPos  = [float(self.dose.ImagePositionPatient[0]/10),
						  float(self.dose.ImagePositionPatient[1]/10),
						  float(self.dose.ImagePositionPatient[2]/10)]
				except AttributeError:
					imSpac = [xnum,ynum,znum]
					say("Dicom File Error",
						"RTDOSE Internal Error")

				[xmin, ymin, zmin] = list(map(lambda x,y: x-y/2, imPos, imSpac))
				[xnum, ynum, znum] = imNum
				xmax = xmin+xnum*imSpac[0]
				ymax = ymin+ynum*imSpac[1]
				zmax = zmin+znum*imSpac[2]

			else:
				xnum = int(voxel.nx)
				ynum = int(voxel.ny)
				znum = int(voxel.nz)
				xmax = xmin+xnum*voxel.dx
				ymax = ymin+ynum*voxel.dy
				zmax = zmin+znum*voxel.dz
		else:
			self.flair.notify("No Voxel Card",
					  "Add Voxel Card from Voxel Tab first")
			return

		pos = inp.bestPosition("USRBIN")

		X = [xmin,xmax,xnum]
		Y = [ymin,ymax,ynum]
		Z = [zmin,zmax,znum]

		pos = self._usrbinCards(X,Y,Z, undoinfo, pos)

	#-----------------------------------------------------------------------
	def makeBeamFile(self,rtplan, undoinfo, inp):
		"""Create the beam data files"""

		dMomentum = self.momentum.getfloat()
		if(self.momentumGauss.get()):
			dMomentum *= -1

		dAngular = self.angular.getfloat()
		if(self.angularGauss.get()):
			dAngular *= -1

		dFWHMx = self.fwhmx.getfloat()
		if(self.xGauss.get()):
			dFWHMx *= -1

		dFWHMy = self.fwhmy.getfloat()
		if(self.yGauss.get()):
			dFWHMy *= -1

		if(self.printed.get()):
			dPrinted = 1
		else:
			dPrinted = 0

		fn = self.rtList.get(self.filesel)[0]

		filename,nSpots = rtplan.beams2Files(self.dataset,(fn).split('.')[0],
							 dMomentum, dAngular, dPrinted, dFWHMx, dFWHMy)

		if not filename:
			self.flair.warning("No File",
			"Beam Profile Files were not created, check your DICOM")
			return None

		card = inp["DEFAULTS"][0]
		card.setSdum("HADROTHE")

		# set the maximum energy for tabulation
		#TODO: Find maximum energy in the beam spots
		card = inp["BEAM"][0]
		card.setWhat(1,(float(-MAX_EN)))
		pos  = inp.bestPosition("SOURCE")

		undoinfo.append(self.flair.addCardUndo(
			(Input.Card("FREE")),pos, False))
		pos+=1
		for i in range(len(filename)):
			pos = self._includeCards(i,filename[i], nSpots[i], undoinfo, pos)
			if ((i==len(filename)-1) and (i>0)):
				undoinfo.append(self.flair.addCardUndo(
					Input.Card("#endif",[]),pos,False))
				pos+=1

		undoinfo.append(self.flair.addCardUndo(
			Input.Card("FIXED"),pos, False))

	#-----------------------------------------------------------------------
	def _patposition(self, iTrasf, setup, undoinfo, pos):
		""" Adding ROT-DEFI card based on VOXEL to GANTRY coordinate system"""

		if setup == "HFS":
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 100, 0, 90],
				   "Voxel Rotation"),pos, False))
			pos+=1

		#For FLUKA release - 2011.1c.5----------------------------------
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 100, 0, -90],
		#		   "Voxel Rotation"),pos, False))
		#	pos+=1
		#---------------------------------------------------------------

		elif setup == "FFS":
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 100, 0, -90],
				   "Voxel Rotation"),pos, False))
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 200, 0, 180],
				   "Voxel Rotation"),pos+1, False))
			pos+=2

		#For FLUKA release - 2011.1c.5 --------------------------------
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 100, 0, 90],
		#		   "Voxel Rotation"),pos, False))
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 200, 0, -180],
		#		   "Voxel Rotation"),pos+1, False))
		#	pos+=2
		#---------------------------------------------------------------

		elif setup == "HFP":
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 100, 0, 90],
				   "Voxel Rotation"),pos, False))
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 200, 0, 180],
				   "Voxel Rotation"),pos+1, False))
			pos+=2

		#For FLUKA release - 2011.1c.5----------------------------------
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 100, 0, -90],
		#		   "Voxel Rotation"),pos, False))
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 200, 0, -180],
		#		   "Voxel Rotation"),pos+1, False))
		#	pos+=2
		#---------------------------------------------------------------

		elif setup == "FFP":
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 100, 0, -90],
				   "Voxel Rotation"),pos, False))
			pos+=1

		#For FLUKA release - 2011.1c.5---------------------------------
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 100, 0, 90],
		#		   "Voxel Rotation"),pos, False))
		#	pos+=1
		#---------------------------------------------------------------

		elif setup == "HFDL":
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 100, 0, 90],
				   "Voxel Rotation"),pos, False))
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 200, 0, 90],
				   "Voxel Rotation"),pos+1, False))
			pos+=2

		#For FLUKA release - 2011.1c.5----------------------------------
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 100, 0, -90],
		#		   "Voxel Rotation"),pos, False))
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 200, 0, -90],
		#		   "Voxel Rotation"),pos+1, False))
		#	pos+=2
		#---------------------------------------------------------------

		elif setup == "HFDR":
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 100, 0, 90],
				   "Voxel Rotation"),pos, False))
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 200, 0, -90],
				   "Voxel Rotation"),pos+1, False))
			pos+=2

		#For FLUKA releaseFor - 2011.1c.5-------------------------------
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 100, 0, -90],
		#		   "Voxel Rotation"),pos, False))
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 200, 0, 90],
		#		   "Voxel Rotation"),pos+1, False))
		#	pos+=2
		#---------------------------------------------------------------

		elif setup == "FFDL":
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 100, 0, 90],
				   "Voxel Rotation"),pos, False))
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 200, 0, -90],
				   "Voxel Rotation"),pos+1, False))
			pos+=2
		#For FLUKA release - 2011.1c.5----------------------------------
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 100, 0, -90],
		#		   "Voxel Rotation"),pos, False))
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 200, 0, 90],
		#		   "Voxel Rotation"),pos+1, False))
		#	pos+=2
		#---------------------------------------------------------------

		elif setup == "FFDR":
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 100, 0, -90],
				   "Voxel Rotation"),pos, False))
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("ROT-DEFI", ["DICOM", 200, 0, 90],
				   "Voxel Rotation"),pos+1, False))
			pos+=2
		#For FLUKA release - 2011.1c.5----------------------------------
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 100, 0, 90],
		#		   "Voxel Rotation"),pos, False))
		#	undoinfo.append(self.flair.addCardUndo(
		#		Input.Card("ROT-DEFI", ["nDICOM", 200, 0, -90],
		#		   "Voxel Rotation"),pos+1, False))
		#	pos+=2
		#---------------------------------------------------------------

		return pos

	#-----------------------------------------------------------------------
	def _isocenter(self, iTrasf, isocenter, undoinfo, pos):
		"""Updating Voxel files - beam isocenter"""

		if iTrasf ==1:
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("#if", ["only_beam"+str(iTrasf)],
					   "Defines isocenter setup for beam"+str(iTrasf)),pos, False))
			pos+=1
		elif iTrasf >1:
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("#elif", ["only_beam"+str(iTrasf),],
					   "Defines isocenter setup for beam"+str(iTrasf)),pos, False))
			pos+=1

		undoinfo.append(self.flair.addCardUndo(Input.Card("ROT-DEFI", ["DICOM", 100, 0, 0,
				-(isocenter[0]/10), -(isocenter[1]/10), -(isocenter[2]/10)],
				"Isocenter translation"),pos, False))
		#---------------------------------------------------------------

		return pos+1

	#-----------------------------------------------------------------------
	def _usrbinCards(self, X, Y, Z, undoinfo, pos):
		"""Adding USRBIN Cards"""

		undoinfo.append(self.flair.addCardUndo(
			Input.Card("USRBIN", ["Dose-H2O", 10,"DOSE-H2O",
				-22,X[1],Y[1],Z[1],X[0],Y[0],Z[0],X[2],Y[2],Z[2]],
				"Dose-H2O"),pos, False))

		card = Input.Card("USRBIN", ["All-Part", 10,"ALL-PART",
				-23,X[1],Y[1],Z[1],X[0],Y[0],Z[0],X[2],Y[2],Z[2]],
				"All particles")
		undoinfo.append(self.flair.addCardUndo(card,pos+1, False))

		#For FLUKA release - 2011.1c.5----------------------------------
		# card.disable()
		#---------------------------------------------------------------

		card=Input.Card("USRBIN", ["Dose", 10,"DOSE",
				-24,X[1],Y[1],Z[1],X[0],Y[0],Z[0],X[2],Y[2],Z[2]],
				"Dose")
		undoinfo.append(self.flair.addCardUndo(card,pos+2, False))
		#For FLUKA release - 2011.1c.5-----------------------------------
		#card.disable()
		#----------------------------------------------------------------

		undoinfo.append(self.flair.addCardUndo(
			Input.Card("ROTPRBIN", ["","","-DICOM","","DOSE-H2O","DOSE"],
				   "ROTBIN"),pos+3, False))

		#For FLUKA release - 2011.1c.5-----------------------------------
		#undoinfo.append(self.flair.addCardUndo(
		#	Input.Card("ROTPRBIN", ["","","nDICOM","","Dbeam","Dbeam"],
		#		   "ROTBIN"),pos+3, False))
		#----------------------------------------------------------------

		return pos+4

	#-------------------------------------------------------------------------
	def _includeCards(self, iTrasf, filename, nSpots, undoinfo, pos):
		"""Adding source file for beam description """

		if iTrasf == 0:
			card = Input.Card("#if", ["full_treatment",],
					   "Defines beam source file for full treatment")
			undoinfo.append(self.flair.addCardUndo(card, pos, False))
			card = Input.Card("SPECSOUR", ["BEAMSPOT", nSpots],
				   "Card to intialize special beam - SPOTBEAM")
			undoinfo.append(self.flair.addCardUndo(card, pos+1, False))
			card = Input.Card("#include", [str(filename)])
			self._inccards.append(card)
			undoinfo.append(self.flair.addCardUndo(card,pos+2, False))
			undoinfo.append(self.flair.addCardUndo(
				Input.Card("#endinclude",[]),pos+3, False))

		elif iTrasf > 0:
			card = Input.Card("#elif", ["only_beam"+str(iTrasf),],
					   "Defines beam source file for beam"+str(iTrasf))
			undoinfo.append(self.flair.addCardUndo(card, pos, False))
			card = Input.Card("SPECSOUR", ["BEAMSPOT", nSpots],
				   "Card to intialize special beam - SPOTBEAM")
			undoinfo.append(self.flair.addCardUndo(card, pos+1, False))
			card = Input.Card("#include", [str(filename)])
			self._inccards.append(card)
			card.disable()
			undoinfo.append(self.flair.addCardUndo(card,pos+2, False))
			card = Input.Card("#endinclude",[])
			card.disable()
			undoinfo.append(self.flair.addCardUndo(card,pos+3, False))

		return pos+4

	#-----------------------------------------------------------------------
	# FIXME for the moment it doesn't go through the undo/redo
	#-----------------------------------------------------------------------
	def expandIncludeFiles(self, cards):
		self.flair.project.input.renumber()
		for card in cards:
			if card.enable:
				self.flair.project.input.include(card)

	#-----------------------------------------------------------------------
	# FIXME: what to do with non correct dicom files?
	#-----------------------------------------------------------------------
	def scan4ReferencedDose(self):
		self.doseInfo = None
		try:
			self.doseInfo = self.page.findDicomFromSOPInstanceUID(
					self.dataset.ReferencedDoseSequence[0].ReferencedSOPInstanceUID)
		except AttributeError:
			self.doseInfo = None
			say("Dicom File Error",
				"RTDOSE with correct UID does not exist")
		return
