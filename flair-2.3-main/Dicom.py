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
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	14-Aug-2012

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import sys
import math
import time
import struct
import string
import getopt
import decimal
import fortran
import fnmatch
from log import say

import bmath
import Input
import tkDialogs

from operator import itemgetter

try:
	from PIL import Image
	from PIL import ImageDraw
	has_pil = True
except ImportError:
	has_pil = False

try:
	import numpy
except ImportError:
	numpy = None

try:
	import pydicom as dicom
except ImportError:
	try:
		import dicom
	except ImportError:
		dicom = None

HOUNSFIELD_MIN = -2000
HOUNSFIELD_MAX =  4096

MODALITY = {
	"AS"       : "Angioscopy",
	"AU"       : "Audio",
	"BI"       : "Biomagnetic imaging",
	"CD"       : "Color flow Doppler",
	"CF"       : "Cinefluorography (retired)",
	"CP"       : "Culposcopy",
	"CR"       : "Computed Radiography",
	"CS"       : "Cystoscopy",
	"CT"       : "Computed Tomography",
	"DD"       : "Duplex Doppler",
	"DF"       : "Digital fluoroscopy (retired)",
	"DG"       : "Diaphanography",
	"DM"       : "Digital microscopy",
	"DS"       : "Digital Subtraction Angiography (retired)",
	"DX"       : "Digital Radiography",
	"EC"       : "Echocardiography",
	"ECG"      : "Electrocardiography",
	"EPS"      : "Cardiac Electrophysiology",
	"ES"       : "Endoscopy",
	"FA"       : "Fluorescein angiography",
	"FS"       : "Fundoscopy",
	"GM"       : "General Microscopy",
	"HC"       : "Hard Copy",
	"HD"       : "Hemodynamic Waveform",
	"IO"       : "Intra-oral Radiography",
	"IVUS"     : "Intravascular Ultrasound",
	"KO"       : "Key Object Selection",
	"LP"       : "Laparoscopy",
	"LS"       : "Laser surface scan",
	"MA"       : "Magnetic resonance angiography",
	"MG"       : "Mammography",
	"MR"       : "Magnetic Resonance",
	"MS"       : "Magnetic resonance spectroscopy",
	"NM"       : "Nuclear Medicine",
	"OCT"      : "Optical Coherence Tomography",
	"OP"       : "Ophthalmic Photography",
	"OPM"      : "Ophthalmic Mapping",
	"OPR"      : "Ophthalmic Refraction",
	"OPV"      : "Ophthalmic Visual Field",
	"OT"       : "Other",
	"PR"       : "Presentation State",
	"PT"       : "Positron emission tomography (PET)",
	"PX"       : "Panoramic X-Ray",
	"REG"      : "Registration",
	"RF"       : "Radio Fluoroscopy",
	"RG"       : "Radiographic imaging (conventional film/screen)",
	"RTDOSE"   : "Radiotherapy Dose",
	"RTIMAGE"  : "Radiotherapy Image",
	"RTPLAN"   : "Radiotherapy Plan",
	"RTRECORD" : "RT Treatment Record",
	"RTSTRUCT" : "Radiotherapy Structure Set",
	"SEG"      : "Segmentation",
	"SM"       : "Slide Microscopy",
	"SMR"      : "Stereometric Relationship",
	"SR"       : "SR Document",
	"ST"       : "Single-photon emission computed tomography (SPECT)",
	"TG"       : "Thermography",
	"US"       : "Ultrasound",
	"VF"       : "Videofluorography (retired)",
	"XA"       : "X-Ray Angiography",
	"XC"       : "External-camera Photography",
}

# ----------------------------------------------------------------------
# Draw line for every pixel that it is intercepted, accepting
# floating points as pixel position
# ----------------------------------------------------------------------
def drawLine(image, x1, y1, x2, y2, c):
	dx = x2-x1
	dy = y2-y1
	length = math.sqrt(dx**2 + dy**2)
	if length < 1e-10: return
	dirx = dx / length
	diry = dy / length

	i  = int(x1)
	j  = int(y1)

	if dirx > 1e-10:
		tx    = (float(i+1) - x1)/ dirx	# distance to next cell
		tdx   = 1.0 / dirx
		stepx = 1
	elif dirx < -1e-10:
		tx    = (float(i) - x1)/ dirx	# distance to next cell
		tdx   = -1.0 / dirx
		stepx = -1
	else:
		tx    = 1e10
		tdx   = 0.0
		stepx = 0

	if diry > 1e-10:
		ty    = (float(j+1) - y1)/ diry	# distance to next cell
		tdy   = 1.0 / diry
		stepy = 1
	elif diry < -1e-10:
		ty    = (float(j) - y1)/ diry	# distance to next cell
		tdy   = -1.0 / diry
		stepy = -1
	else:
		ty    = 1e10
		tdy   = 0.0
		stepy = 0

	try:
		image.putpixel((i,j),c)
	except IndexError:
		pass
	while tx<=length or ty<=length:
		if tx<ty:
			i += stepx
			tx += tdx
		else:
			j += stepy
			ty += tdy
		if i<0 or i>=image.size[0] or j<0 or j>=image.size[1]: break
		try:
			image.putpixel((i,j),c)
		except IndexError:
			pass

	i2 = int(x2)
	j2 = int(y2)
	try:
		image.putpixel((i2,j2),c)
	except IndexError:
		pass

# ----------------------------------------------------------------------
# Return the header of a dicom for filling the list
# ----------------------------------------------------------------------
def dicomHeader(uid,filename):
	try:
		dataset = dicom.read_file(filename,
			force=True,
			defer_size=256,
			stop_before_pixels=True)
	except:
		return (uid,"?","?","?","?","?","?","?")

	try:
		modality = dataset.Modality
	except:
		modality = 'dicom'

	try:
		acqdate = dataset.StudyDate
		acqtime = dataset.StudyTime
		datetime = "%s.%s.%s %s:%s"% \
			(acqdate[: 4], acqdate[ 4: 6], acqdate[ 6:],
			 acqtime[:-4], acqtime[-4:-2])
	except AttributeError:
		datetime = "?"

	if "PatientIdentityRemoved" in dataset:
		name = "-patient identity removed-"
	elif "PatientsName" in dataset:
		pat_name = dataset.PatientsName
		name = pat_name.family_name + ", " + pat_name.given_name
	else:
		name = "none"

	try:
		institution_name = dataset.InstitutionName
	except AttributeError:
		institution_name = "none"

	try:
		study_description = dataset.StudyDescription
	except AttributeError:
		study_description = "none"

	try:
		age = dataset.PatientsAge
	except AttributeError:
		try:
			birth = dataset.PatientBirthDate[0:4]
			age = str(time.gmtime()[0] - int(birth))
		except:
			age = "?"

	try:
		sex = dataset.PatientsSex
	except AttributeError:
		sex = "?"

	return (uid,
		  modality,
		  datetime,
		  name,
		  age,
		  sex,
		  institution_name,
		  study_description)

#===============================================================================
# Scan a directory for dicom sets and cache the results
#===============================================================================
class Scan:
	directory = ""
	datasets  = {}

	# ----------------------------------------------------------------------
	# Scan directory for unique dicom files
	# ----------------------------------------------------------------------
	@staticmethod
	def scan(master, path):
		apath = os.path.abspath(path)
		if apath == Scan.directory: return

		dlg = tkDialogs.ProgressDialog(master, "Scanning directory: %s"%(path))
		Scan.directory = apath
		Scan.datasets.clear()

		try:
			listdir = os.listdir(apath)
			dlg.setLimits(0.0, len(listdir))
			for i,fn in enumerate(listdir):
				filename = os.path.join(path, fn)
				try:
					dataset = dicom.read_file(filename,
						force=True,
						defer_size=256,
						stop_before_pixels=True)
				except:
					continue
				try:
					uid = dataset.SeriesInstanceUID
				except AttributeError:
					continue

				#if "SliceLocation" in dataset:
				#	z = float(dataset.SliceLocation)
				if "ImagePositionPatient" in dataset:
					z = float(dataset.ImagePositionPatient[2])
				else:
					try:
						z = float(dataset.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0].ImagePositionPatient[2])
					except:
						z = 0.0

				try:
					sop = dataset.SOPInstanceUID
				except:
					sop = None

				try:
					modality = dataset.Modality
				except:
					modality = 'dicom'

				if uid in Scan.datasets:
					# append file
					Scan.datasets[uid][1].append((fn,z,sop))
				else:
					try:
						frames = int(dataset.NumberOfFrames)
					except AttributeError:
						frames = 0

					Scan.datasets[uid] = [frames, [(fn,z,sop)], modality]

				# Show progress
				if dlg.show(i,"Scanning: "+fn):
					break
				del dataset
		except OSError:
			pass

		# Sort lists
		for (frames,files,modality) in list(Scan.datasets.values()):
			files.sort(key=lambda x: Scan.datasets[uid][1])

		dlg.stop()

	# ----------------------------------------------------------------------
	@staticmethod
	def _cmpZ(a, b):
		return cmp(a[1], b[1])

	# ----------------------------------------------------------------------
	@staticmethod
	def rescan(master, path):
		Scan.directory = ""
		Scan.scan(master, path)

	# ----------------------------------------------------------------------
	@staticmethod
	def __getitem__(uid):
		return Scan.datasets[uid]

	# ----------------------------------------------------------------------
	@staticmethod
	def headerUID(uid):
		return dicomHeader(uid,os.path.join(Scan.directory, Scan.files(uid)[0][0]))

	# ----------------------------------------------------------------------
	@staticmethod
	def frames(uid):
		return Scan.datasets[uid][0]

	# ----------------------------------------------------------------------
	@staticmethod
	def files(uid):
		return Scan.datasets[uid][1]

	# ----------------------------------------------------------------------
	@staticmethod
	def modality(uid):
		return Scan.datasets[uid][2]

#===============================================================================
# Hounsfield units class
#===============================================================================
class Hounsfield:
	prev   = -9999

	#-----------------------------------------------------------------------
	def __init__(self, humax, mat, crho_min, crho_max, cdedx_min, cdedx_max):
		self.humin  = Hounsfield.prev
		self.humax  = humax
		Hounsfield.prev = humax
		self.mat    = mat
		self.reg    = 0
		self.crho_min  = crho_min
		self.crho_max  = crho_max
		self.cdedx_min = cdedx_min
		self.cdedx_max = cdedx_max

		if self.humax > self.humin+1:
			self._rho_scale  = (self.crho_max  - self.crho_min)  / float(self.humax - self.humin - 1)
			self._dedx_scale = (self.cdedx_max - self.cdedx_min) / float(self.humax - self.humin - 1)
		else:
			self._rho_scale  = 0.0
			self._dedx_scale = 0.0

#		if abs(self.crho(self.humin) - self.crho_min) > 0.001:
#			print "error1", self
#		if abs(self.crho(self.humax-1) - self.crho_max) > 0.001:
#			print
#			print "error2", self
#			print "interpolated:",self.crho(self.humax-1)
#			print "max:",self.crho_max

	#-----------------------------------------------------------------------
	def __str__(self):
		return "[%5d .. %5d)	%s  rho:[%g .. %g]  dedx:[%g .. %g]" % \
			(self.humin, self.humax, self.mat,
			 self.crho_min, self.crho_max,
			 self.cdedx_min, self.cdedx_max)

	#-----------------------------------------------------------------------
	def str(self):
		return "%d %s  %.10g %.10g  %.10g  %.10g" % \
			(self.humax, self.mat,
			 self.crho_min, self.crho_max,
			 self.cdedx_min, self.cdedx_max)

	#-----------------------------------------------------------------------
	def crho(self, h):
		return (h - self.humin) * self._rho_scale  + self.crho_min

	#-----------------------------------------------------------------------
	def cdedx(self, h):
		return (h - self.humin) * self._dedx_scale + self.cdedx_min

	#-----------------------------------------------------------------------
	@staticmethod
	def reset():
		Hounsfield.prev = -9999

#===============================================================================
# Dicom helper class
#===============================================================================
class Dicom:
	#-----------------------------------------------------------------------
	def __init__(self, filename=None, log=None):
		self.setLog(log)
		self.clear()
		self.clearHounsfield()
		self.xfrom = None
		self.yfrom = None
		self.zfrom = None
		self.xto   = None
		self.yto   = None
		self.zto   = None

		self.sliceFrom = 0
		self.sliceTo   = sys.maxsize
		self.rowFrom   = 0
		self.rowTo     = sys.maxsize
		self.colFrom   = 0
		self.colTo     = sys.maxsize
		self.field     = 0

		if filename is not None:
			self.init(filename)

	#-----------------------------------------------------------------------
	def clear(self):
		self.uid      = None
		self.title    = ""
		self.slice    = False
		self._slices  = []
		self._scanned = False
		self.modality = "?"

		self.x  = self.y  = self.z  = 0.0
		self.dx = self.dy = self.dz = 0.0
		self.nx = self.ny = self.nz = 0

	#-----------------------------------------------------------------------
	def clearHounsfield(self):
		Hounsfield.reset()
		self.hounsfield = []
		self.corrrho    = []
		self.corrdedx   = []
		self.correction = {}
		self.material   = Input.Input()

	# ----------------------------------------------------------------------
	def setLog(self, log):
		self._log = log

	# ----------------------------------------------------------------------
	def log(self, s):
		if self._log is None:
			sys.stdout.write("%s\n"%(s))
		else:
			self._log(s)

	# ----------------------------------------------------------------------
	def scanned(self):	return self._scanned

	# ----------------------------------------------------------------------
	# @return true if they have changed
	# ----------------------------------------------------------------------
	def setLimits(self, xfrom=None, xto=None, yfrom=None, yto=None, zfrom=None, zto=None):
		self.xfrom = xfrom
		self.yfrom = yfrom
		self.zfrom = zfrom
		self.xto = xto
		self.yto = yto
		self.zto = zto
		return self.calcLimits()

	# ----------------------------------------------------------------------
	def calcLimits(self):
		sf = self.sliceFrom
		st = self.sliceTo

		rf = self.rowFrom
		rt = self.rowTo

		cf = self.colFrom
		ct = self.colTo

		#print "calcLimits: x",self.xfrom, self.xto, self.x, self.dx
		#print "calcLimits: y",self.yfrom, self.yto, self.y, self.dy
		#print "calcLimits: z",self.zfrom, self.zto, self.z, self.dz
		try:    self.sliceFrom = min(max(0, int(math.floor((float(self.zfrom)-self.z)/self.dz))),self.nz-1)
		except: self.sliceFrom = 0

		try:    self.sliceTo   = min(max(0, int(math.floor((float(self.zto)-self.z)/self.dz))),self.nz-1)
		except: self.sliceTo   = self.nz-1

		try:    self.rowFrom   = min(max(0, int(math.floor((float(self.yfrom)-self.y)/self.dy))),self.ny-1)
		except: self.rowFrom   = 0

		try:    self.rowTo     = min(max(0, int(math.floor((float(self.yto)-self.y)/self.dy))),self.ny-1)
		except: self.rowTo     = self.ny-1

		try:    self.colFrom   = min(max(0, int(math.floor((float(self.xfrom)-self.x)/self.dx))),self.nx-1)
		except: self.colFrom   = 0

		try:    self.colTo     = min(max(0, int(math.floor((float(self.xto)-self.x)/self.dx))),self.nx-1)
		except: self.colTo     = self.nx-1

		#print "calcLimits:", self.colFrom,self.colTo, self.rowFrom, self.rowTo, self.sliceFrom, self.sliceTo
		return sf!=self.sliceFrom or st!=self.sliceTo or \
		       rf!=self.rowFrom   or rt!=self.rowTo   or \
		       cf!=self.colFrom   or ct!=self.colTo

	#-----------------------------------------------------------------------
	def columns(self): return min(self.nx, self.colTo  -self.colFrom+1)
	def rows(self):    return min(self.ny, self.rowTo  -self.rowFrom+1)
	def slices(self):  return min(self.nz, self.sliceTo-self.sliceFrom+1)

	def xmin(self):	return self.x + (self.colFrom-0.5)  *self.dx
	def ymin(self): return self.y + (self.rowFrom-0.5)  *self.dy
	def zmin(self): return self.z + (self.sliceFrom-0.5)*self.dz

	def xmax(self):	return self.x + (self.columns()+ self.colFrom-0.5)  *self.dx
	def ymax(self):	return self.y + (self.rows()   + self.rowFrom-0.5)  *self.dy
	def zmax(self):	return self.z + (self.slices() + self.sliceFrom-0.5)*self.dz

	#-----------------------------------------------------------------------
	# Read Hounsfield material and ranges from a file
	#-----------------------------------------------------------------------
	def readHounsfield(self, filename):
		Hounsfield.reset()
		del self.hounsfield[:]

		try: f = open(filename,"r")
		except: return False

		lines = []
		for line in f:
			line = line.strip()
			if len(line)==0 or line[0] in ("#","*"): continue
			elem     = line.split()
			humax    = int(elem[0])
			lines.append((humax, elem))

		lines.sort(key=itemgetter(0))

		for humax,elem in lines:
			material = elem[1]
			try:	crho_min = float(elem[2])
			except: crho_min = 1.0
			try:	crho_max = float(elem[3])
			except: crho_max = 1.0
			try:	cdedx_min = float(elem[4])
			except: cdedx_min = 1.0
			try:	cdedx_max = float(elem[5])
			except: cdedx_max = 1.0
			hu = Hounsfield(humax, material, crho_min, crho_max, cdedx_min, cdedx_max)
			self.hounsfield.append(hu)

		f.close()
		self.interpolateCorrections()

	#-----------------------------------------------------------------------
	def interpolateCorrections(self):
		self.correction.clear()
		for hu in self.hounsfield:
			# Interpolate correction
			if hu.humax > hu.humin:
				for h in range(hu.humin, hu.humax):
					self.correction[h] = (-hu.cdedx(h), hu.crho(h))
			else:
					self.correction[hu.humin] = (-hu.cdedx_min, hu.crho_min)

	#-----------------------------------------------------------------------
	def readMaterials(self, filename):
		self.material.read(filename)

	#-----------------------------------------------------------------------
	def readCorrection(self, filename):
		self.correction = {}
		try: f = open(filename,"r")
		except: return False
		for line in f:
			line = line.strip()
			if len(line)==0 or line[0]=="*": continue
			elem = list(map(float,line.split()))
			self.correction[int(elem[0])] = elem[1:]
		f.close()

	#-----------------------------------------------------------------------
	# Initialize dicom information by reading the first slice
	#-----------------------------------------------------------------------
	def init(self, filename, nz=0, z=0.0):
		self.title    = ""
		self.modality = ""
		self._slices  = []
		self.values   = {}
		self.m = 1.0
		self.b = 0.0
		self.minimum =  99999999999
		self.maximum = -99999999999
		self.uid = None
		self.dx  = 0.0
		self.dy  = 0.0
		self.dz  = 0.0
		self.nx  = 0
		self.ny  = 0
		self.nz  = nz
		self.x   = 0.0
		self.y   = 0.0
		self.z   = z
		self.slice    = False
		self._scanned = False
		self.hasPixel = False

		# Load one slice to get dimensions
		if isinstance(filename,str):
			try:
				dataset = dicom.read_file(filename, force=True)
			except:
				self.log("Error cannot open file %s"%(filename))
				return False
		else:
			dataset = filename

		# Very slow for RTSTRUCT
#		dataset.decode()

		self.uid = dataset.SeriesInstanceUID

		if "PatientIdentityRemoved" in dataset:
			name = "-patient identity removed-"
		elif "PatientsName" in dataset:
			pat_name = dataset.PatientsName
			name = pat_name.family_name + ", " + pat_name.given_name
		else:
			name = "none"
		self.title = "%-80s"%(name)

		try:
			self.modality = dataset.Modality
		except:
			self.modality = 'dicom'

		if "PixelData" not in dataset:
			self.hasPixel = False
			#self.log("Error cannot handle DICOM without PixelData")
			return False
		else:
			self.hasPixel = True

		if "PixelSpacing" not in dataset:
			# Most probably new 3D format
			#import pdb; pdb.set_trace()

			self.dx = self.dy = None
			self.dz = 0.0
			found = False

			# Try first in SharedFunctionalGroupsSequence
			if "SharedFunctionalGroupsSequence" in dataset:
				found = True
				frame = dataset.SharedFunctionalGroupsSequence[0]
				try:
					self.dx, self.dy = list(map(float,frame.PixelMeasuresSequence[0].PixelSpacing))
					self.dz = float(frame.PixelMeasuresSequence[0].SliceThickness)
				except:
					pass

			# Then in PerFrameFunctionalGroupsSequence
			if "PerFrameFunctionalGroupsSequence" in dataset:
				found = True
				frame = dataset.PerFrameFunctionalGroupsSequence[0]
				if "PixelMeasuresSequence" in frame:
					try:
						self.dx, self.dy = list(map(float,frame.PixelMeasuresSequence[0].PixelSpacing))
					except:
						pass
					if "SliceThickness" in frame.PixelMeasuresSequence[0]:
						self.dz = float(frame.PixelMeasuresSequence[0].SliceThickness)
					elif "SliceThickness" in frame[0x2005140f][0]:	# Private tag?
						self.dz = float(frame[0x2005140f][0].SliceThickness)

			if not found:
				self.log("Unknown format: No PixelSpacing or" \
					" PerFrameFunctionalGroupsSequence or" \
					" SharedFunctionalGroupsSequence found in DICOM")
				return False

			if self.dx is None or self.dy is None:
				self.log("Unknown dicom format, PixelSpacing is not found.")
				return False

		else:
			# in mm
			self.dx, self.dy = list(map(float,dataset.PixelSpacing))
			try:
				self.dz = float(dataset.SliceThickness)
			except:
				self.dz = 0.0

			frame = None

		#if dataset.RescaleType != "HU":
		#	self.log("Error cannot handle DICOM with rescale type different than HU")
		#	return False

		self.m = 1.0
		self.b = 0.0

		if frame is not None:
			try:
				self.m = float(frame.PixelValueTransformationSequence[0].RescaleSlope)
				self.b = float(frame.PixelValueTransformationSequence[0].RescaleIntercept)
			except AttributeError:
				pass
		elif self.modality == "RTDOSE":
			try :
				self.m = float(dataset.DoseGridScaling)
			except AttributeError:
				pass
		else:
			try:
				self.m = float(dataset.RescaleSlope)
				self.b = float(dataset.RescaleIntercept)
			except AttributeError:
				pass

		if frame is not None:
			try:
				self.x, self.y, z = list(map(float,frame.PlanePositionSequence[0].ImagePositionPatient))
			except:
				self.x = self.y = z = 0.0
		else:
			try:
				self.x, self.y, z = list(map(float,dataset.ImagePositionPatient))
			except:
				self.x = self.y = z = 0.0

		# get limits
		self.nx = int(dataset.Columns)
		self.ny = int(dataset.Rows)
		if 'NumberOfFrames' in dataset:
			# 3D data
			self.nz = int(dataset.NumberOfFrames)
			self.slice = False
			self.z  = z
			self.sliceZ = 0.0
		else:
			# 2D slice data
			self.slice = True
			self.sliceZ = z
			if 'SliceLocation' not in dataset and 'NumberOfSlices' not in dataset:
				self.log("Error cannot handle DICOM without both SliceLocation and NumberOfFrames %s" \
					%(name))
				return False

		# Calculate the limits for the first time
		self.calcLimits()
		if isinstance(filename,str): del dataset
		return True

	#-----------------------------------------------------------------------
	# Scan slice for unique Hounsfield values
	#-----------------------------------------------------------------------
	def _scanSlice(self, filename, s=None, update=None):
		try: dataset = dicom.read_file(filename, force=True)
		except: return False

		if self.uid != dataset.SeriesInstanceUID:
			self.log("Different Series Instance UID found")
			return False

		# Scan for unique Hounsfield regions
		if self.slice:
			if s is None or self.sliceFrom <= s <= self.sliceTo:
				self.__scanSliceData(dataset.pixel_array)
			try:
				self._slices.append((float(dataset.ImagePositionPatient[2]), filename))
			except:
				self._slices.append((float(dataset.SliceLocation), filename))

		else:
			try:
				# FIXME Only for DOSE maps
				zoff = dataset.GridFrameOffsetVector
			except:
				zoff = [x*self.dz for x in range(len(dataset.pixel_array))]

			for s,data_slice in enumerate(dataset.pixel_array):
				self._slices.append((float(zoff[s]), filename))
				if s is None or self.sliceFrom <= s <= self.sliceTo:
					self.__scanSliceData(data_slice)
				if update and update(s): break

		del dataset
		return True

	#-----------------------------------------------------------------------
	# Scan data to create a histogram of unique Hounsfield values
	#-----------------------------------------------------------------------
	def __scanSliceData(self, data):
		if self.rowFrom>0 or self.rowTo<self.ny-1 or \
		   self.colFrom>0 or self.colTo<self.nx-1:
			array = data[self.rowFrom:self.rowTo+1, self.colFrom:self.colTo+1]
		else:
			array = data
		mi = array.min()
		ma = array.max()
		hist,val = numpy.histogram(array, ma-mi+1, (mi-0.5,ma+0.5))
		for i in range(hist.size):
			if hist[i] > 0:
				v = int(math.ceil(val[i]))	# closest int is the upper one
				if v not in self.values:
					self.values[v] = hist[i]
				else:
					self.values[v] += hist[i]

	#-----------------------------------------------------------------------
	def _scanSliceEnd(self):
		self._scanned = True

		keys = sorted(self.values.keys())
		self.minimum = int(self.m*float(min(keys)) + self.b)
		self.maximum = int(self.m*float(max(keys)) + self.b)

		# Sorted list of existing Hounsfield units in all slices
		# Each Hounsfield is one region
		# XXX WARNING the conversion should be the same with the _writeVoxelSlice ###
		self.reg2hu = [int(self.m*float(val)+self.b) for val in keys]
		self.reg2hu.sort()

		# DEBUG
#		if Input._developer:
#			say("Limits:",self.minimum, self.maximum)
#			self.printHistogram()
#			say("\nreg2hu")
#			for i,r in enumerate(self.reg2hu):
#				say(i,r)

		# inverse table
		self.kreg = [0]*(self.maximum - self.minimum)

		# Keep first Unit (normally -1024) as the master volume organ=0
		for i,h in enumerate(self.reg2hu[1:]):
			self.kreg[h-self.minimum-1] = i+1

#		if Input._developer:
#			say("\nkreg")
#			for i,r in enumerate(self.kreg):
#				say(i,r)

		# Sort slices
		self._slices = sorted(self._slices, key=lambda x:x[0])
		if self.nz==0: self.nz = len(self._slices)

		# Correct slice thickness, use the average from the last SliceLocation to the first
		self.dz = (self._slices[-1][0] - self._slices[0][0]) / (len(self._slices)-1)

		prev = self._slices[0][0] - self.dz
		err = ""
		dzs = []
		for l,fn in self._slices:
			dzs.append(l-prev)
			if abs((l-prev) - self.dz) > 1E-5:
				self.log("\nSlice: %s"%(fn))
				self.log("Difference from previous: %g"%(l-prev))
				self.log("Report slice thickness:   %g"%(self.dz))
				if not err:
					err = "Warning: Slice is missing or wrong thickness."
				err +=  "\nLocation z= %s\n" \
					"Difference from previous: %s\n" \
					"Slice thickness: %s" \
					% (str(l),str(l-prev), str(self.dz))
			prev = l
		return err

	#-----------------------------------------------------------------------
	# Scan directory for all files with a specific pattern and sort them
	# in slice order
	#-----------------------------------------------------------------------
	def scanSlices(self, path, slices):
		self.init(os.path.join(path,slices[0]))
		for fn in slices:
			if not self._scanSlice(os.path.join(path,fn)):
				return False
		return self._scanSliceEnd()

	#-----------------------------------------------------------------------
	def assignMaterials2Regions(self, field=0):
		self.field = field
		self.assignmat = []
		for i,h in enumerate(self.reg2hu):
			for hounsfield in self.hounsfield:
				if hounsfield.humin <= h < hounsfield.humax:
					self.assignmat.append(hounsfield.mat)
					break
			else:
				self.log("WARNING: Out of range HU= %d assigned last material"%(h))
				self.assignmat.append(self.hounsfield[-1].mat)
		return True

	#-----------------------------------------------------------------------
	def _write(self, data):
		fortran.write(self._f, data)

	#-----------------------------------------------------------------------
	def _writeStr(self, data):
		self._write(data.encode("utf-8"))

	#-----------------------------------------------------------------------
	def _writeStruct(self, fmt, *kw):
		fortran.write(self._f, struct.pack(fmt, *kw))

	#-----------------------------------------------------------------------
	# Initialize the voxel file
	# Write the header
	#-----------------------------------------------------------------------
	def _writeVoxelInit(self, filename):
		try: self._f = open(filename,"wb")
#		try: self._f = open(filename,"w")
		except: return False

		# voxel data
		self._data = []

		# initialize the rtstuct arrays
		self._roiname   = {}	# index=ROINumber, value=name,volume
		self._roicomb   = {}	# index=:RTSTRUCT[:RTSTRUCT] (Unique combination) value=sequential index
		self._roidata   = []	# matrix of seq-index
		self._roicolor  = {}	# color of each structure
		self._roiplanar = {}	# planar information for CLOSED_PLANAR

		self._tdraw = 0.0
		self._tproc = 0.0
		return True

	#-----------------------------------------------------------------------
	# Populate the self._data array with the conversion
	# to voxel from a single slice or a complete 3D
	#-----------------------------------------------------------------------
	def _writeVoxelSlice(self, filename, s=None, update=None):
		if s is not None and not self.sliceFrom <= s <= self.sliceTo: return
		dataset = dicom.read_file(filename, force=True)
		if self.slice:
			self.__writeVoxelSliceData(dataset.pixel_array)
		else:
			for s,data_slice in enumerate(dataset.pixel_array):
				if s is None or self.sliceFrom <= s <= self.sliceTo:
					self.__writeVoxelSliceData(data_slice)
				if update and update(s): break
		del dataset

	#-----------------------------------------------------------------------
	def __writeVoxelSliceData(self, data):
		if self.rowFrom>0 or self.rowTo<self.ny-1 or \
		   self.colFrom>0 or self.colTo<self.nx-1:
			array = data[self.rowFrom:self.rowTo+1, self.colFrom:self.colTo+1]
		else:
			array = data

		array2 = numpy.uint16(self.m*array + self.b - self.minimum)
		self._data.append(bytes(array2.data))

	#-----------------------------------------------------------------------
	# Find the ROI (RTSTRUCT) and create the 3D matrix as well the
	# correspondence table
	#-----------------------------------------------------------------------
	def _roiStructSlice(self, dataset, sop):
		if not has_pil: return
		rois   = []
		images = []

		# Prepare table with index, names, volume
		for roi in dataset.StructureSetROISequence:
			try:
				self._roiname[int(roi.ROINumber)] = str(roi.ROIName)
			except:
				pass

		# Loop over all ROIs
		t = time.time()
		for sequence in dataset.ROIContourSequence:
			# Loop over all contours to find SOP
			draw = None
			roi = int(sequence.ReferencedROINumber)
			planar = []

			if roi not in self._roicolor:
				try:
					color = sequence.ROIDisplayColor
					# Swap color order to be used in geoviewer
					self._roicolor[roi] = int(color[0])<<16 | int(color[1])<<8 | int(color[2])
				except AttributeError:
					self._roicolor[roi] = 128

			for contour in sequence.ContourSequence:
				if contour.ContourImageSequence[0].ReferencedSOPInstanceUID != sop:
					continue

				# Ignore non-planar (for the moment)
				if contour.ContourGeometricType != "CLOSED_PLANAR" or \
				   contour.NumberOfContourPoints == 1:
					continue

				# We might have more than one contour per slice
				if draw is None:
					# First time create the image/draw
					rois.append(int(sequence.ReferencedROINumber))
					image = Image.new("L", (self.columns(),self.rows()))
					images.append(image)
					draw = ImageDraw.Draw(image)

				# FIXME cropped images!!!
				coords = list(map(float, contour.ContourData))
				X  = [x-self.x+self.dx/2.0 for x in coords[ ::3]]	# WARNING center of voxel
				Y  = [y-self.y+self.dy/2.0 for y in coords[1::3]]
				zi = int(coords[2] - self.z / self.dz)	# z-index of structure

				# FIXME me correct the origin!
				planar.append((zi,list(zip(X,Y))))

				iX = [int(x/self.dx) for x in X]
				iY = [int(y/self.dy) for y in Y]
				draw.polygon(list(zip(iX,iY)), fill=0x80)

				# draw outline with our routine to correctly mark all
				# possible crossing voxels for the outline, based on the
				# true floating point position of the structure
				x1 = X[-1]/self.dx
				y1 = Y[-1]/self.dy
				for x,y in zip(X,Y):
					x2 = x/self.dx
					y2 = y/self.dy
					drawLine(image, x1, y1, x2, y2, 0xff)
					x1,y1 = x2,y2

				#image.save("ROI%02d_z%03d.png"%(roi, zi))

			if planar: self._roiplanar[roi] = planar

			del draw
		self._tdraw += time.time() - t

		roimg = [(roi,img.getdata())
				for roi,img in zip(rois,images)]

		# Process images
		t = time.time()
		ptr  = 0
		data = ""
		# XXX XXX XXX XXX XXX XXX XXX
		# FIXME cropped images!!!
		# XXX XXX XXX XXX XXX XXX XXX
		# Warning ZERO based
		self._roicomb[""] = 0	# first element
		for j in range(self.rows()):
			col = []
			for i in range(self.columns()):
				# rois appearing in pixel as unique string :roi[:roi]
				pixel = ""
				for roi,img in roimg:
					if   img[ptr]==0x80: pixel += ":%d"%(roi)
					elif img[ptr]==0xff: pixel += ":-%d"%(roi)
				ptr += 1
				# try:except is faster than if None or setdefault
				try:
					idx = self._roicomb[pixel]
				except KeyError:
					# add new id
					self._roicomb[pixel] = idx = len(self._roicomb)
				col.append(idx)
			data += struct.pack("=%dH"%(len(col)), *col)
		self._roidata.append(data)
		del data
		del images
		self._tproc += time.time() - t

	#-----------------------------------------------------------------------
	# Finalize the file by writing the data, materials, etc..
	#-----------------------------------------------------------------------
	def _writeVoxelEnd(self):
		# write title
		self._writeStr(self.title)

		# memory size
		header = struct.pack("=5i", self.columns(), self.rows(), self.slices(),
						len(self.reg2hu)-1, len(self.kreg))
		if len(self._roicomb)>1:
			# Number of structures (consecutive)
			header += struct.pack("=i", max(self._roiname.keys()))

			# Number of combinations
			header += struct.pack("=i", len(self._roicomb)-1)	# Minus the 0 index

			# Maximum length of combinations
			maxcomb = max([x.count(":") for x in list(self._roicomb.keys())])
			header += struct.pack("=i", maxcomb)	# Minus the 0 index

			# Number of CLOSED_PLANAR structures
			header += struct.pack("=i", len(self._roiplanar))

		self._write(header)

		# voxel dimension
		self._writeStruct("=3d", 0.1*self.dx, 0.1*self.dy, 0.1*self.dz)

		# data
		self._write(b"".join(self._data))
		del self._data

		# hu -> region matrix
		self._writeStruct("=%dH"%(len(self.kreg)), *self.kreg)

		# rt struct
		self.__writeVoxelRTstruct()

		# embed material definition
		self.writeMaterials()

		# assignmat matrix
		self.writeAssignmats()

		# correction factors
		self.writeCorrFactors()

		self._f.close()
		return True

	#-----------------------------------------------------------------------
	# Write the RTSTRUCT if they are there
	#-----------------------------------------------------------------------
	def __writeVoxelRTstruct(self):
		if len(self._roicomb)<=1: return

		# Write structure information
		for i in range(1,max(self._roiname.keys())+1):
			name  = self._roiname.get(i,"")
			color = self._roicolor.get(i,0)
			self._writeStruct("=i64s", color, name+" "*(64-len(name)))
			#self._write(name+" "*(64-len(name)))

		# Write indexes (WARNING skip 0 which is empty)
		maxcomb = max([x.count(":") for x in list(self._roicomb.keys())]) # maximum combination
		for i in range(1,len(self._roicomb)):
			# Find the key that has this index
			for key,value in list(self._roicomb.items()):
				if value == i:
					# Found break the key into the ROI sequence
					rois = sorted(map(int,key[1:].split(":")))
					lenrois = len(rois)
					rois.extend([0]*(maxcomb-lenrois))
					say("*-*", i,lenrois,len(rois),rois)
					self._writeStruct("=i%dh"%(len(rois)),lenrois,*rois)
					#self._writeStruct("=i%dH"%(len(rois)),lenrois,*rois)
					break

		# Write the data
		self._write(" ".join(self._roidata,""))

		#say("RTSTRUCT len=",len(string.join(self._roidata,"")))
		#say("Dictionary len=",len(self._roicomb),"max=",\)
		#	max([len(x.split(":")) for x in self._roicomb.keys()])

		###################
		# Write planar information
		if self._roiplanar:
			for roi in sorted(self._roiplanar.keys()):
				planar = self._roiplanar[roi]
				self._writeStruct("=ii", roi, len(planar))
				for zi,xy in planar:
					rec = struct.pack("=ii",zi,len(xy))
					for x,y in xy:
						rec += struct.pack("=ff",x,y)
					self._write(rec)
		###################

		# Cleanup memory
		del self._roiplanar
		del self._roicomb
		del self._roidata

	#-----------------------------------------------------------------------
	# Write USRBIN header
	#-----------------------------------------------------------------------
	def _writeUsrbinInit(self, filename):
		try: self._f = open(filename,"wb")
		except: return False

		# write title
		self._writeStruct("=80s32sfiii", ((self.title+80*" ").encode()),((time.asctime()+32*" ").encode()),0.0,0,0,0)

		# write dicom as first usrbin file
		self._writeStruct("=i10siiffifffifffififff",
				1, ((self.modality+" "*10)[:10]).encode(), 0, 208,
				self.xmin()/10.0, self.xmax()/10.0, self.columns(), self.dx/10.0,
				self.ymin()/10.0, self.ymax()/10.0, self.rows(),    self.dy/10.0,
				self.zmin()/10.0, self.zmax()/10.0, self.slices(),  self.dz/10.0,
				0, 0.0, 0.0, 0.0)

		self._data = []
		self._min =  1e10
		self._max = -1e10
		return True

	#-----------------------------------------------------------------------
	# Process a single slice of USRBIN
	#-----------------------------------------------------------------------
	def _writeUsrbinSlice(self, filename, s=None, update=None):
		if s is not None and not self.sliceFrom <= s <= self.sliceTo: return
		dataset = dicom.read_file(filename, force=True)

		m = self.m
		b = self.b

#		if frame is not None:
#			try:
#				m = float(frame.PixelValueTransformationSequence[0].RescaleSlope)
#				b = float(frame.PixelValueTransformationSequence[0].RescaleIntercept)
#			except AttributeError:
#				pass
		if self.modality == "RTDOSE":
			try :
				m = float(dataset.DoseGridScaling)
			except AttributeError:
				pass
		else:
			try:
				m = float(dataset.RescaleSlope)
				b = float(dataset.RescaleIntercept)
			except AttributeError:
				pass

		if self.slice:
			self.__writeUsrbinSliceData(dataset.pixel_array, m, b)
		else:
			for s,data_slice in enumerate(dataset.pixel_array):
				if s is None or self.sliceFrom <= s <= self.sliceTo:
					self.__writeUsrbinSliceData(data_slice, m, b)
				if update and update(s): break
		del dataset

	#-----------------------------------------------------------------------
	# append data
	#-----------------------------------------------------------------------
	def __writeUsrbinSliceData(self, data, m, b):
		if self.rowFrom>0 or self.rowTo<self.ny-1 or \
		   self.colFrom>0 or self.colTo<self.nx-1:
			array = data[self.rowFrom:self.rowTo+1, self.colFrom:self.colTo+1]
		else:
			array = data
		array2 = numpy.float32(m*array + b)
		self._min = min(self._min, array2.min())
		self._max = max(self._max, array2.max())
		self._data.append(bytes(array2.data))

	#-----------------------------------------------------------------------
	# Finalize USRBIN by writing the data
	#-----------------------------------------------------------------------
	def _writeUsrbinEnd(self):
		# data
		self._write(b"".join(self._data))
		del self._data
		self._f.close()
		say("Limits:", self._min, self._max)
		return True

	#-----------------------------------------------------------------------
	# Convert slices to voxels
	#-----------------------------------------------------------------------
	def writeVoxel(self, filename):
		if not self._writeVoxelInit(filename): return False
		for l,fn in self._slices:
			self._writeVoxelSlice(fn)

		return self._writeVoxelEnd()

	#-----------------------------------------------------------------------
	# Embed materials inside voxel file
	#-----------------------------------------------------------------------
	def writeMaterials(self):
		for tag in ["MATERIAL", "COMPOUND", "MAT-PROP"]:
			for card in self.material[tag]:
				if card.ignore(): continue
				scard = card.toStr(Input.FORMAT_SINGLE)
				for line in scard.splitlines():
#					self._writeStr("%-80s"%(line[:80]))
#					self._writeStruct("=80s",(line[:80]))
#					self._writeStruct("=80s",(line[:80]+b'\0'*80))
					self._writeStruct("=80s",(line+b' '*80))

		return True

	#-----------------------------------------------------------------------
	# Write assignmats
	#-----------------------------------------------------------------------
	def writeAssignmats(self):
		if self.field:
			extra="%20s%8d.0"%("",self.field)
		else:
			extra=""
		self._writeStr("ASSIGNMA  %10s%10s%-50s" % ("VACUUM","VOXEL",extra))
		for i,asmat in enumerate(self.assignmat):
			self._writeStr("ASSIGNMA  %10s%10s%-50s" % \
				(asmat[:8],Input.Voxel.regionName(i+1),extra))

	#-----------------------------------------------------------------------
	# Write correction factors
	#-----------------------------------------------------------------------
	def writeCorrFactors(self):
		if not self.correction: return
		card = Input.Card("CORRFACT",[""])
		for i,hu in enumerate(self.reg2hu):
			try:
				for j,c in enumerate(self.correction[hu]):
					card.setWhat(j+1, c)
				card.setWhat(4, Input.Voxel.regionName(i+1))
				scard = card.toStr()
#				self._writeStr("%-80s"%(scard))
				self._writeStruct("=80s",(scard+b' '*80))

			except:
				pass

	#-----------------------------------------------------------------------
	def printHistogram(self):
		say("Dicom histogram")
		for v in sorted(self.values.keys()):
			say(v, int(self.m*float(v) + self.b), self.values[v])

#===============================================================================
# RT structure processing
#===============================================================================
class RTStruct:
	#-----------------------------------------------------------------------
	def __init__(self, dataset, ctinfo=None):
		self.idx     = {}	# index to name
		self.names   = {}	# name to index
		self.volumes = bmath.DefaultDict((0.,0.,0.))	# volumes per index

		if isinstance(dataset,str):
			try:
				dataset = dicom.read_file(dataset, force=True)
			except:
				return

		if ctinfo is not None:
			ds = Dicom(ctinfo)
			self.x  = ds.x
			self.y  = ds.y
			self.dx = ds.dx
			self.dy = ds.dy
			self.dz = ds.dz
			self.nx = ds.nx
			self.ny = ds.ny
		else:
			self.x  = self.y  = None
			self.dx = self.dy = self.dz = None
			self.nx = self.ny = None

		self.findNames(dataset)
		self.calcVolumes(dataset)

	#-----------------------------------------------------------------------
	def findNames(self, dataset):
		for roi in dataset.StructureSetROISequence:
			self.names[int(roi.ROINumber)] = str(roi.ROIName)
			self.idx[str(roi.ROIName)]     = int(roi.ROINumber)
			try:
				volume = float(roi.ROIVolume)
			except:
				volume = 0.0
			self.volumes[int(roi.ROINumber)] = (volume,0.,0.)

	#-----------------------------------------------------------------------
	# Calculate real volumes based on coordinates
	#-----------------------------------------------------------------------
	def calcVolumes(self, dataset):
		# Loop over all sequences
		dz = 1e99
		for sequence in dataset.ROIContourSequence:
			# Loop over all contours
			zmin =  1e99
			zmax = -1e99
			zp   = None
			Sv   = 0
			Sp   = 0.0
			if self.nx is not None:
				image = Image.new("L",(self.nx, self.ny))
				draw = ImageDraw.Draw(image)
			else:
				image = None

			for contour in sequence.ContourSequence:
				# Ignore non-planar (for the moment)
				if contour.ContourGeometricType != "CLOSED_PLANAR" or \
				   contour.NumberOfContourPoints == 1:
					continue

				# Coordinates
				coords = list(map(float, contour.ContourData))
				X = coords[ ::3]
				Y = coords[1::3]
				z = coords[2]
				zmin = min(zmin,z)
				zmax = max(zmax,z)

				# Calculate volume from polygon
				S = 0.0
				for i in range(len(X)):
					x = X[i]
					y = Y[i]
					j = i+1
					if j>=len(X): j=0
					xn = X[j]
					yn = Y[j]
					S += (x+xn) * (y-yn)

				# Voxel-based volume
				if image:
					X  = [x-self.x+self.dx/2.0 for x in X]	# WARNING center of voxel
					Y  = [y-self.y+self.dy/2.0 for y in Y]
					iX = [int(x/self.dx) for x in X]
					iY = [int(y/self.dy) for y in Y]

					# clear image
					draw.rectangle(((0,0),(self.nx,self.ny)), fill=0, outline=0)
					draw.polygon(list(zip(iX,iY)), fill=1)

					# draw outline with our routine to correctly mark all
					# possible crossing voxels for the outline, based on the
					# true floating point position of the structure
					x1 = X[-1]/self.dx
					y1 = Y[-1]/self.dy
					for x,y in zip(X,Y):
						x2 = x/self.dx
						y2 = y/self.dy
						drawLine(image, x1, y1, x2, y2, 1)
						x1,y1 = x2,y2

					# Some all entries with weight=1
					Sv += image.histogram()[1]
				Sp += abs(S)

				# Find next z
				if zp is not None and z != zp:
					dz = min(dz, z - zp)
				zp = z

			if dz is not None:
				if image:
					#print "***",sequence.ReferencedROINumber,Sv
					Sv *= self.dx*self.dy*self.dz
				Sp = Sp * dz / 2.0
				Rv = self.volumes[int(sequence.ReferencedROINumber)][0]
				self.volumes[int(sequence.ReferencedROINumber)] = (Rv,Sv,Sp)

#===============================================================================
def usage(rc):
	say("DICOM to voxel")
	say("syntax: Dicom.py [options]")
	say("\t-?\t\thelp")
	say("\t-h\t--hounsfield=<file>\tload file with Hounsfield units to material conversion [no-default]")
	say("\t-m\t--material=<file>\tload fluka input file for material definitions [default:none]")
	say("\t-c\t--correction=<file>\tload file with density correction factors [default:none]")
	say("\t-d\t--dir=<path>\tpath where dicom files are [default: .]")
	say("\t-p\t--pattern=<filepattern>\tfiles to load e.g CT0*.dcm [default: *.dcm]")
	say("\t-o\t--output=<voxelfile>\toutput voxel file [default: output.vxl]")
	sys.exit(rc)

if __name__ == "__main__":
	if dicom is None:
		say("pydicom module cannot be loaded")

	if len(sys.argv)==1:
		usage(0)

	try:
		optlist, args = getopt.getopt(sys.argv[1:],
			'?h:m:c:d:p:o:',
			['help', 'hounsfield=','material=', 'correction=', 'dir=', 'pattern=',"output="])
	except getopt.GetoptError:
		usage(1)

	Input.init()
	dd = Dicom()

	path    = "."
	pattern = "*.dcm"
	output  = "dicom.vxl"

	for opt, val in optlist:
		if opt in ("-?", "--help"):
			usage(0)
		elif opt in ("-h", "--hounsfield="):
			say("Hounsfield:", val)
			dd.readHounsfield(val)
		elif opt in ("-m", "--material="):
			say("Materials:", val)
			dd.readMaterials(val)
		elif opt in ("-c", "--correction"):
			say("Correction:",val)
			dd.readCorrection(val)
		elif opt in ("-d","--dir="):
			say("Path:", val)
			path = val
		elif opt in ("-p","--pattern="):
			say("Pattern:", val)
			pattern = val
		elif opt in ("-o","--output="):
			say("Output:", val)
			output = val
		else:
			say("Unknown option",opt,val)
			usage(2)

	first = True
	for fn in os.listdir(path):
		if not fnmatch.fnmatch(fn, pattern): continue
		filename = os.path.join(path,fn)
		if first:
			dd.init(filename)
			first = False
		dd._scanSlice(filename)
	dd._scanSliceEnd()

#	dd.scanSlices(path,pattern) # "/clueet/home/mchin/DICOM/lung_petct","T33001_CT*.dcm"
	dd.assignMaterials2Regions()
	dd.writeVoxel(output)
