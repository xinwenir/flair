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
 */

#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "os.h"
#include "geo.h"
#include "ray.h"
#include "units.h"

/*
enum ParticleType {
	PNeutron,
	PProton,
	PGamma,
	PElectron,
	PPositron,
	PDeuteron,
	PTriton,
	PAlpha,
	PIon,
	PPionPlus,
	PPionMinus,
};
*/

/* =============================== Particle ================================= */
class Particle : public Ray {
public:
//	ParticleType	type;
//	PType	type;			// particle type
//	double	mass;
//	double	energy;

public:
	Particle()	{}
virtual const char*	name()	const = 0;
virtual double		mass()	const = 0;
}; // Particle

/* =============================== Caloron ================================== *
 * Fictitious particle transfering heat
 */
class Caloron : public Particle {
public:
virtual const char*	name()	const { return "Caloron"; }
virtual double		mass()	const { return 0.0; }
} ;

/* =============================== Neutron ================================== */
class Neutron : public Particle {
public:
virtual const char*	name()	const { return "Neutron"; }
virtual double		mass()	const { return Mn; }
}; // Neutron

#endif
