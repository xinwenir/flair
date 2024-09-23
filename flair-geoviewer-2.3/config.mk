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

# ******************************
# Default variables
# ******************************
ifndef SYSTEM
	SYSTEM = UNIX
endif

ifndef THREAD
	THREAD = yes
endif

ifndef MEM
	MEM = no
endif

# ******************************
# File extensions
# ******************************
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
	SYSTEM = $(UNAME)
	SO  = so
	OBJ = o
	EXE =
else ifneq (,$(findstring CYGWIN,$(UNAME)))
	SYSTEM = CYGWIN
	SO  = dll
	OBJ = o
	EXE = .exe
else
	SO  = so
	OBJ = o
	EXE =
endif

# ******************************
# Compilers
# ******************************
ifeq ($(COMP),clang)
#ifeq ($(CLANG),yes)
	CXX := clang++
else
	CXX := c++
endif

# ******************************
# Architecture
# ******************************
ifndef TARGET
	TARGET := $(shell uname -m)
endif

ifeq ($(TARGET),i386)
	ARCH      = -m32

else ifeq ($(TARGET),i586)
	ARCH      = -m32

else ifeq ($(TARGET),i686)
	ARCH      = -m32

else ifeq ($(TARGET),x86_64)
	ARCH      = -m64

else ifeq ($(TARGET),amd64)
	ARCH      = -m64

else ifeq ($(TARGET),arm64)
	ARCH      = -m64

else ifeq ($(TARGET),android)
	ANDNDK    = /opt/android-ndk
	ANDOS     = 18
	ANDVER    = 4.8
	ANDARCH   = arm
	#ANDBIT    = linux-x86
	#ANDBIT    = linux-x86_64
	#ANDBIN   := $(ANDNDK)/toolchains/$(ANDARCH)-linux-androideabi-$(ANDVER)/prebuilt/$(ANDBIT)/bin
	ANDBIN    = $(shell ls -d $(ANDNDK)/toolchains/$(ANDARCH)-linux-androideabi-$(ANDVER)/prebuilt/*)/bin
	ANDROOT  := $(ANDNDK)/platforms/android-$(ANDOS)/arch-$(ANDARCH)
	ANDARM   := $(ANDBIN)/$(ANDARCH)-linux-androideabi
	ANDSTL   := $(ANDNDK)/sources/cxx-stl/stlport
	ARCH      = --sysroot=$(ANDROOT)
	CXX      := $(ANDARM)-c++
	PYX      := no
	SYSTEM   := ANDROID
	THREAD   := no
else
	#$(shell echo "Error TARGET architecture not found")
endif
