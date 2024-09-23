# $Id: makefile 3727 2016-02-18 08:33:29Z bnv $
#
# Copyright and User License
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright Vasilis.Vlachoudis@cern.ch for the
# European Organization for Nuclear Research (CERN)
#
# All rights not expressly granted under this license are reserved.
#
# Installation, use, reproduction, display of the
# software ("flair"), in source and binary forms, are
# permitted free of charge on a non-exclusive basis for
# internal scientific, non-commercial and non-weapon-related
# use by non-profit organizations only.
#
# For commercial use of the software, please contact the main
# author Vasilis.Vlachoudis@cern.ch for further information.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.
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
# Date:	07-Apr-2016
#

LIBRARY = lib$(NAME).a
SOLIB   = $(NAME).$(SO)
OBJECTS = $(SOURCES:.cc=.$(OBJ))
DEPS    = $(SOURCES:.cc=.d)

# **************************
#  CXX Flags
# **************************
# _DEBUG can be 0, 1 or 2
# NDEBUG is used by assert()
CXXFLAGS = -D$(SYSTEM)

# **************************
#  Debug + Optimization
# **************************
ifeq ($(DEBUG),yes)
	# With debug info.
	CXXFLAGS += -g
	CXXFLAGS += -D_DEBUG=2
else
	# Optimised version
#ifeq ($(CLANG),yes)
ifeq ($(COMP),clang)
#	CXXFLAGS += -O3
	CXXFLAGS += -O
	CXXFLAGS += -ftrapv
else
	CXXFLAGS += -O4
endif
	CXXFLAGS += -funroll-loops
#	Floating point, prefer sse and if not then 387
#	CXXFLAGS += -ffast-math
#	CXXFLAGS += -ftree-vectorizer-verbose=5
#	CXXFLAGS += -ftree-vectorize
#	CXXFLAGS += -mfpmath=sse,387
#	CXXFLAGS += -mfpmath=387
	# XXX for the beta release XXX
	CXXFLAGS += -g
	CXXFLAGS += -DNDEBUG
	CXXFLAGS += -D_DEBUG=0
	#CXXFLAGS += -mtune=native
	#CXXFLAGS += -march=native
endif

# **************************
#  Defines
# **************************
ifeq ($(THREAD),yes)
	CXXFLAGS += -DTHREAD
endif

ifeq ($(MEM),yes)
	CXXFLAGS += -DMEM
endif
# Use the DUMP=yes to enable logging
ifeq ($(DUMP),yes)
	CXXFLAGS += -D_DUMP
endif

ifeq ($(EXP),yes)
	CXXFLAGS += -DEXPERIMENTAL
endif

ifeq ($(OPENGL),yes)
	CXXFLAGS += -DOPENGL
endif

ifeq ($(STAT),yes)
	CXXFLAGS += -D_STAT
endif

ifeq ($(PROFILE),yes)
	CXXFLAGS += -pg
endif
ifeq ($(PROF),yes)
	CXXFLAGS += -pg
endif

ifeq ($(OPENMP),yes)
	CXXFLAGS += -fopenmp
endif

# **************************
#  Warnings
# **************************
#CXXFLAGS += -W
CXXFLAGS += -Wall
CXXFLAGS += -Wcast-align
CXXFLAGS += -Wcast-qual
#CXXFLAGS += -Wconversion
CXXFLAGS += -Wextra
CXXFLAGS += -Wformat
CXXFLAGS += -Wpointer-arith
CXXFLAGS += -Wredundant-decls
CXXFLAGS += -Wshadow
CXXFLAGS += -Wno-write-strings
#ifeq ($(CXX),c++)
#	CXXFLAGS += -Wconversion-null
#endif
#CXXFLAGS += -Wfloat-equal
#CXXFLAGS += -Winline
#CXXFLAGS += -Wno-deprecated
#CXXFLAGS += -Werror


# **************************
#  Warnings
# **************************
# XXX important for floating point speed up. Look fedora-14 release notes
#CXXFLAGS += -fexcess-precision=fast	# for gcc 4.5
# Signed arithmetic overflow of +,-,* wraps using 2-complement
CXXFLAGS  += -fwrapv
CXXFLAGS  += -fstack-protector
ALIASING   = -fstrict-aliasing
NOALIASING = -fno-strict-aliasing

ifeq ($(SYSTEM),UNIX)
	CXXFLAGS += -fpic
else ifeq ($(SYSTEM),ANDROID)
	CXXFLAGS += -fpic
	CXXFLAGS += -fexceptions
#	CXXFLAGS += -frtti
endif

ifeq ($(SYSTEM),Darwin)
	LDFLAGS  = -dynamiclib
else
	SOFLAGS  = -shared
	LDFLAGS += -Wl,-O4
endif

# **************************
# Include directories
# **************************
CXXFLAGS += $(INCDIR)

# **************************
# Link flags
# **************************
LDFLAGS += $(LIBDIR) $(LIBS)

# **************************
# Generic compilation of any program
# **************************
%: %.cc $(LIBRARY)
	$(CXX) -o $@ $< $(ARCH) $(CXXFLAGS) -I. $(CPPFLAGS) $(LDFLAGS)

# **************************
# Compiling of cc programs
# **************************
%.$(OBJ): %.cc
	$(CXX) -c $< $(ARCH) $(CXXFLAGS) $(ALIASING)

# **************************
# Dependencies rule
# **************************
%.d: %.cc
	$(CXX) -MM $(INCDIR) $< > $@

#.PHONY: doc
#doc:
#	doxygen doxygen.cfg 2>&1 | tee doxygen.out
