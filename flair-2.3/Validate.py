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
# Date:	17-May-2004

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import sys
import pdb
from log import say

import bmath
import Input
import Project

#===============================================================================
# Validate input
#===============================================================================
class Validate:
	def __init__(self, project, log=None):
		self.project = project
		self.input   = project.input
		self._log    = log
		self.skip    = []	# errors/warnings to skip
		self.init()

	# ----------------------------------------------------------------------
	def init(self):
		# List of cards with errors/warnings
		self.warnings = []
		self.errors   = []
#		self.fatal    = []

	# ----------------------------------------------------------------------
	def log(self, s):
		if self._log is None:
			sys.stdout.write("%s\n"%(s))
		else:
			self._log(s)

	# ----------------------------------------------------------------------
	def write(self, s):
		self.log(s)

	# ----------------------------------------------------------------------
	def addWarning(self, card, msg):
		self.warnings.append((card,msg))
		self.write("\nWarning: %s\n%s"%(msg,str(card)))

	# ----------------------------------------------------------------------
	def addError(self, card, msg):
		self.errors.append((card,msg))
		self.write("\nERROR: %s\n%s"%(msg,str(card)))

	#-----------------------------------------------------------------------
	# Check specific run for logical errors
	#-----------------------------------------------------------------------
	def check(self, run=None):
		if run is not None:
			defines = run.defines
		else:
			defines = None

		# We have to check the problematic cards!!!
		for card,err in self.input.preprocess(defines):
			self.addError(card,err)

		# Check cards out of order
		self._checkOutOfOrder()

		# Check for start..stop transformations
		# self._checkTransformations()

		# Assign the region properties
		self.input.regionProperties()

		# Assignment
		# - Bodies not used
		# - Regions not assigned properties
		#   . Material
		#   . EMF threshold (if EMF cut is defined)
		#   . BIAS if bias is used

		# Physics check
		# Thresholds
		# - limits
		# - photon/electron in PART-THRES
		# - EMF-CUT (PROD without TRANSPORT or vice versa)
		#           EMF-CUT with no material/region definition
		# - maybe all region ranges if they are missing

		# Scoring
		# - Units mixing
		# - RESNUCLEi with w/wo evolution on same unit
		# - USRBDX of region that do not share any common body

	#-----------------------------------------------------------------------
	# Check for cards out of logical order
	#-----------------------------------------------------------------------
	def _checkOutOfOrder(self):
		header    = ("TITLE", "GLOBAL", "DEFAULTS")
		present   = {}	# Cards present in input

		bodies    = {}	# lists for checking for unused cards
		regions   = {}
		materials = {}
		assignmat = {}
		rotdefi   = {}

		def exist(tag):		return present.get(tag, 0)>0
		def count(tag):		return present.get(tag, 0)
		def mark(tag):
			if tag in present:
				present[tag] += 1
			else:
				present[tag] = 1
		def beforeGeometry():	return not exist("GEOBEGIN")
		def inGeometry():	return exist("GEOBEGIN") and not exist("GEOEND")
		def afterGeometry():	return exist("GEOEND")
		# Assuming we re in geometry
		def inBodies():		return count("END")==0 and exist("body")
		def inRegions():	return count("END")==1 and exist("REGION")
		def inLattice():	return count("END")==2
#		def inVolume():
		def addGeometry(card):
			if card.type()==Input.Card.BODY:
				pass
			else:
				pass

		say()
		say("-"*70)

		# Check cards
		for card in self.input.cardlist:
			if card.ignore(): continue
			tag = card.tag
			if tag[0] == "#": continue	# Preprocessor skip

			# FIXME Special $transform cards??????

			# ---- Check for header cards ----
			say(">>>",tag)
			if tag in header:
				if present.get(tag,False):
					self.addWarning(card, "Duplicate: %s"%(tag))

				elif exist("normal"):
					self.addWarning(card,
						"Out of order: GLOBAL,DEFAULTS. Must be the first cards")

			# ---- Before geometry ----
			elif beforeGeometry() or afterGeometry():
				if card.isGeo():
					self.addError(card,
						"Out of order: %s outside the GEOBEGIN...GEOEND"%(tag))
					addGeometry(card)

				elif tag == "GEOBEGIN" and afterGeometry():
					self.addError(card, "Duplicate: %s"%(tag))

				elif tag == "GEOEND":
					self.addError(card, "Out of order: %s"%(tag))

			# ---- Check for geometry ----
			elif inGeometry():
				# VOXELS must be the first in the geometry
				if tag=="VOXELS":
					if exist(tag):
						self.addError(card, "Duplicate: %s"%(tag))

					elif inBodies() or inRegions() or inLattice():
						self.addError(card,"Out of order: VOXELS must be immediately after GEOBEGIN")

				# From lattice we can continue on lattice or with GEOEND (or with VOLUMES?)
				if inLattice() and tag not in ("LATTICE","GEOEND"):
					self.addError(card, "Out of order: %s, LATTICE or GEOEND expected here"%(tag))

				elif inRegions() and tag not in ("REGION", "END"):
					self.addError(card, "Out of order: %s, REGION or END expected here"%(tag))

				elif card.type() == Input.Card.BODY:
					mark("body")

# accept $cards...
				elif inBodies() and (card.type != Input.Card.BODY and tag!="END"):
					self.addError(card, "Out of order: %s, body card or END expected here"%(tag))

				elif count("END") > 2:
					self.addError(card, "Too many END")

			# Mark cards
			if tag not in header: mark("normal")
			mark(tag)

"""
	order = [("#&define"  , _P_DEFINE_layout4 ),
		( "#&undef"   , _P_UNDEF_layout4  ),
		( "#inc&lude" , _P_INCLUDE_layout4),
		( "#&if"      , _P_IF_layout4     ),
		( "#ifdef"    , _P_IFDEF_layout4  ),
		( "#ifndef"   , _P_IFNDEF_layout4 ),
		( "#eli&f"    , _P_ELIF_layout4   ),
		( "#&else"    , _P_ELSE_layout4   ),
		( "#e&ndif"   , _P_ENDIF_layout4  ),
		( "#endinclude", _empty_layout4   ),
		( "TITLE"     , _TITLE_layout4    ),
		( "FIXED"     , _empty_layout4    ),
		( "FREE"      , _empty_layout4    ),
		( "GLOBAL"    , _GLOBAL_layout4   ),
		( "DEFAULTS"  , _DEFAULTS_layout4 ),
		( "&BEAM"     , _BEAM_layout4     ),
		( "BEAM&POS"  , _BEAMPOS_layout4  ),
		( "BEAM&AXES" , _BEAMAXES_layout4 ),
		( "HI-PROPE"  , _HI_PROPE_layout4 ),
		( "SO&URCE"   , _SOURCE_layout4   ),
		( "SPECSOUR"  , _SPECSOUR_layout4 ),
		( "GCR-SPE"   , _GCR_SPE_layout4  ),
		( "POLARIZA"  , _POLARIZA_layout4 ),
		( "DISCARD"   , _DISCARD_layout4  ),
		( "BME"       , _BME_layout4      ),
		( "DELTARAY"  , _DELTARAY_layout4 ),
		( "DPMJET"    , _DPMJET_layout4   ),
		( "EMFFLUO"   , _EMFFLUO_layout4  ),
		( "EMFRAY"    , _EMFRAY_layout4   ),
		( "EVENTYPE"  , _EVENTYPE_layout4 ),
		( "IONFLUCT"  , _IONFLUCT_layout4 ),
		( "IONTRANS"  , _IONTRANS_layout4 ),
		( "MULSOPT"   , _MULSOPT_layout4  ),
		( "MUPHOTON"  , _MUPHOTON_layout4 ),
		( "MYRQMD"    , _MYRQMD_layout4   ),
		( "OPT-PROD"  , _OPT_PROD_layout4 ),
		( "PAIRBREM"  , _PAIRBREM_layout4 ),
		( "PHOTONUC"  , _PHOTONUC_layout4 ),
		( "PHYSICS"   , _PHYSICS_layout4  ),
		( "RQMD"      , _RQMD_layout4     ),

		( "!point"    , _POINT_layout4    ),
		( "!arrow"    , _ARROW_layout4    ),
		( "!ruler"    , _RULER_layout4    ),
		( "!light"    , _LIGHT_layout4    ),
		( "!mesh"     , _MESH_layout4     ),

		( "GEOBEGIN"  , _GEOBEGIN_layout4 ),
		( "VOXELS"    , _VOXELS_layout4   ),
		( "&RPP"      , _RPP_layout4      ),
		( "&TET"      , _TET_layout4      ), 
		( "&BOX"      , _BOX_layout4      ),
		( "&SPH"      , _SPH_layout4      ),
		( "R&CC"      , _RCC_layout4      ),
		( "R&EC"      , _REC_layout4      ),	
		( "&TRC"      , _TRC_layout4      ),
		( "&TRX"      , _TRXYZ_layout4    ),
		( "&TRY"      , _TRXYZ_layout4    ),
		( "&TRZ"      , _TRXYZ_layout4    ),
		( "E&LL"      , _ELL_layout4      ),
		( "&WED"      , _WED_layout4      ),
		( "RAW"       , _RAW_layout4      ),
		( "&ARB"      , _ARB_layout4      ),
		( "&PLA"      , _PLA_layout4      ),
		( "&QUA"      , _QUA_layout4      ),
		( "&YZP"      , _YZP_layout4      ),
		( "X&ZP"      , _XZP_layout4      ),
		( "&XYP"      , _XYP_layout4      ),
		( "XCC"       , _XCC_layout4      ),
		( "YCC"       , _YCC_layout4      ),
		( "ZCC"       , _ZCC_layout4      ),
		( "XEC"       , _XEC_layout4      ),
		( "YEC"       , _YEC_layout4      ),
		( "ZEC"       , _ZEC_layout4      ),
		( "END"       , _empty_layout4	  ),
		( "$start_expansion", _START_EXPANSION_layout4 ),
		( "$start_translat",  _START_TRANSLAT_layout4  ),
		( "$start_transform", _START_TRANSFORM_layout4 ),
		( "$end_expansion",   _END_EXPANSION_layout4   ),
		( "$end_translat",    _END_TRANSLAT_layout4    ),
		( "$end_transform",   _END_TRANSFORM_layout4   ),
		( "REGION"    , _REGION_layout4   ),
		( "LATTICE"   , _LATTICE_layout4  ),
		( "GEO&END"   , _GEOEND_layout4   ),
		( "&MATERIAL" , _MATERIAL_layout4 ),
		( "COMPOUND"  , _COMPOUND_layout4 ),
		( "LOW-MAT"   , _LOW_MAT_layout4  ),
		( "&OPT-PROP" , _OPT_PROP_layout4 ),
		( "MAT-&PROP" , _MAT_PROP_layout4 ),
		( "CORRFACT"  , _CORRFACT_layout4 ),
		( "STERNHEI"  , _STERNHEI_layout4 ),
		( "ASSIGNMA"  , _ASSIGNMA_layout4 ),
		( "EMF"       , _EMF_layout4      ),
		( "EMF&CUT"   , _EMFCUT_layout4   ),
		( "EMFFIX"    , _EMFFIX_layout4   ),
		( "ELCFIELD"  , _ELCFIELD_layout4 ),
		( "EXPTRANS"  , _EXPTRANS_layout4 ),
		( "&FLUKAFIX" , _FLUKAFIX_layout4 ),
		( "&LOW-NEUT" , _LOW_NEUT_layout4 ),
		( "MCSTHRES"  , _MCSTHRES_layout4 ),
		( "MGNFIELD"  , _MGNFIELD_layout4 ),
		( "&PART-THR" , _PART_THR_layout4 ),
		( "STEPSIZE"  , _STEPSIZE_layout4 ),
		( "THRESHOL"  , _THRESHOL_layout4 ),
		( "TIME-&CUT" , _TIME_CUT_layout4 ),
		( "BIASING"   , _BIASING_layout4  ),
		( "EMF-BIAS"  , _EMF_BIAS_layout4 ),
		( "LAM-BIAS"  , _LAM_BIAS_layout4 ),
		( "LOW-BIAS"  , _LOW_BIAS_layout4 ),
		( "LOW-DOWN"  , _LOW_DOWN_layout4 ),
		( "WW-FACTO"  , _WW_FACTO_layout4 ),
		( "WW-PROFI"  , _WW_PROFI_layout4 ),
		( "WW-THRES"  , _WW_THRES_layout4 ),
		( "BAMJET"    , _BAMJET_layout4   ),
		( "DPM-PARA"  , _DPM_PARA_layout4 ),
		( "EVXTEST"   , _EVXTEST_layout4  ),
		( "DCYSCORE"  , _DCYSCORE_layout4 ),
		( "DCYTIMES"  , _DCYTIMES_layout4 ),
		( "IRRPROFI"  , _IRRPROFI_layout4 ),
		( "RADDECAY"  , _RADDECAY_layout4 ),
		( "ROT-&DEFI" , _ROT_DEFI_layout4 ),
		( "ROTPRBIN"  , _ROTPRBIN_layout4 ),
		( "AUXSCORE"  , _AUXSCORE_layout4 ),
		( "TCQUENCH"  , _TCQUENCH_layout4 ),
		( "DETECT"    , _DETECT_layout4   ),
		( "USR&BIN"   , _USRBIN_layout4   ),
		( "USRBD&X"   , _USRBDX_layout4   ),
		( "USR&COLL"  , _USRCOLL_layout4  ),
		( "USR&TRACK" , _USRTRACK_layout4 ),
		( "USR&YIELD" , _USRYIELD_layout4 ),
		( "&EVENTBIN" , _EVENTBIN_layout4 ),
		( "EVENTDAT"  , _EVENTDAT_layout4 ),
		( "&USERDUMP" , _USERDUMP_layout4 ),
		( "USER&WEIG" , _USERWEIG_layout4 ),
		( "&RESNUCLE" , _RESNUCLE_layout4 ),
		( "&SCORE"    , _SCORE_layout4    ),
		( "OPEN"      , _OPEN_layout4     ),
		( "USRGCALL"  , _USRGCALL_layout4 ),
		( "USRICALL"  , _USRICALL_layout4 ),
		( "USROCALL"  , _USROCALL_layout4 ),
		( "PLOTGEOM"  , _PLOTGEOM_layout4 ),
		( "&RANDOMIZ" , _RANDOMIZ_layout4 ),
		( "&START"    , _START_layout4    ),
		( "STOP"      , _STOP_layout4     )]
"""

#-------------------------------------------------------------------------------
if __name__ == "__main__":
	p = Project.Project()
	p.load(sys.argv[1])
	v = Validate(p)
	v.check()
