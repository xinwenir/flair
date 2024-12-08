/*
 * $Id$
 *
 * Copyright and User License
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright 2006-2019 CERN and INFN
 * 
 *
 * Please consult the LICENSE file for the license 
 *
 * DISCLAIMER
 * ~~~~~~~~~~
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
 * NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY, OF
 * SATISFACTORY QUALITY, AND FITNESS FOR A PARTICULAR PURPOSE
 * OR USE ARE DISCLAIMED. THE COPYRIGHT HOLDERS AND THE
 * AUTHORS MAKE NO REPRESENTATION THAT THE SOFTWARE AND
 * MODIFICATIONS THEREOF, WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADE SECRET OR OTHER PROPRIETARY RIGHT.
 *
 * LIMITATION OF LIABILITY
 * ~~~~~~~~~~~~~~~~~~~~~~~
 * THE COPYRIGHT HOLDERS AND THE AUTHORS SHALL HAVE NO
 * LIABILITY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
 * CONSEQUENTIAL, EXEMPLARY, OR PUNITIVE DAMAGES OF ANY
 * CHARACTER INCLUDING, WITHOUT LIMITATION, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES, LOSS OF USE, DATA OR PROFITS,
 * OR BUSINESS INTERRUPTION, HOWEVER CAUSED AND ON ANY THEORY
 * OF CONTRACT, WARRANTY, TORT (INCLUDING NEGLIGENCE), PRODUCT
 * LIABILITY OR OTHERWISE, ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGES.
 *
 * Author:	Vasilis.Vlachoudis@cern.ch
 * Date:	19-Feb-2013
 */

#ifndef __UNITS_H
#define __UNITS_H

// Definition of units for easier reading

// Mathematical units
const double fwhm   = 2.354820045030949327014; // 2.*sqrt(2.*log(2.))

// Prefixes
const double yotta  = 1e24;
const double zetta  = 1e21;
const double exa    = 1e18;
const double peta   = 1e15;
const double tera   = 1e12;
const double giga   = 1e9;
const double mega   = 1e6;
const double kilo   = 1e3;
const double hecto  = 1e2;
const double deca   = 1e1;

const double deci   = 1e-1;
const double centi  = 1e-2;
const double milli  = 1e-3;
const double micro  = 1e-6;
const double nano   = 1e-9;
const double pico   = 1e-12;
const double femto  = 1e-15;
const double atto   = 1e-18;
const double zepto  = 1e-21;
const double yocto  = 1e-24;

// Distance
const double nm     = 0.1E-6;
const double um     = 0.1E-3;
const double mm     = 0.1;
const double cm     = 1.0;
const double dm     = 10.0;
const double m      = 100.0;
const double km     = 100.0e3;

// Imperial system distances
const double inch   = 2.54;
const double feet   = 30.48;
const double ft     = 30.48;
const double mile   = 160934.4;
const double mi     = 160934.4;

// Time
const double ns     = 1.0e-9;
const double us     = 1.0e-6;
const double ms     = 0.001;
const double s      = 1.0;
const double min    = 60.0;
const double hour   = 3600.0;
const double day    = 86400.0;
const double week   = 7.0*86400.0;
const double month  = 365.25/12.0*86400.0;
const double year   = 365.25*86400.0;

// Energy;
const double eV     = 1e-9;
const double keV    = 1e-6;
const double MeV    = 1e-3;
const double GeV    = 1.0;
const double TeV    = 1e3;
const double PeV    = 1e6;
const double J      = 1.0/1.60217646e-10;

// Angle
const double deg    = PI/180.0;
const double rad    = 1.0;
const double mrad   = 0.001;

// Physics
const double a      = 1.0 / 137.035989561;	// fine constant
const double h      = 4.135667516e-15 *eV*s;	// eV s
const double c      = 299792458e2*cm/s;		// cm/s
const double Na     = 6.0221367e23;		// Avogadro number
const double qe     = 1.60217646e-19;		// electron charge
const double re     = 2.8179409183694872e-13;	// electron radius

// Masses
const double amu    = 0.93149432*GeV;
const double amuC12 = 0.93123890*GeV;
const double amugr  = 1.6605402E-24;
const double Mp     = 0.93827231*GeV;		// Proton mass
const double Mn     = 0.93956563*GeV;		// Neutron mass
const double Me     = 0.510998910e-3*GeV;	// Electron mass

#endif
