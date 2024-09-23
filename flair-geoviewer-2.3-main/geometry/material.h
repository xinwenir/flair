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
 * Date:        26 Feb 2013
 */

#ifndef __MATERIAL_H
#define __MATERIAL_H

#include <iostream>

#include "geo.h"
#include "painter.h"

#define MAT_CACHE         (0x0001)
#define MAT_NO_ANTIALIAS  (0x0002)
#define MAT_TRANSFORM     (0x0004)

#define MAT_TM_DIFF       (0x0008)
#define MAT_TM_SPEC       (0x0010)
#define MAT_TM_TRANS      (0x0020)
#define MAT_TM_AMB        (0x0040)
#define MAT_TM_MAPPING    (0x0078)        /* all bits for mapping */

/* ============================== Material ============================ */
class Material {
protected:
	int		_id;		// material id
	char16		_name;		// Material name

	// Physics parameters
	double		_density;

	int		_Z;		// atomic number
	int		_A;		// atomic mass
	int		_weight;	// atomic weight
//	Array		_elements;


	// Optical properties
	dword		_diff;		// diffuse color
	double		_spec;		// specular (reflected) color
	double		_shine;		// specular spot exponent
	double		_ior;		// index of refraction
	double		_fuzz;		// Surface fuzz
//	Textmap		tm_diff;
//	Textmap		tm_spec;
//	dword		_amb;		// ambient color
//	Textmap		tm_amb;
//	dword		_cshine;	// spec spot color
//	dword		_trans;		// Transparent component
//	Textmap		tm_trans;

	dword		_flag;		// is this material valid for shadow caching

//	struct t_texture        *tex;   // ptr for color texture
//	struct t_bump           *bump;  // ptr for surface normal texture
//	Matrix  matrix;                 // transformation matrix

private:
	bool		_used;		// if this material is being used

public:
	Material();
	Material(const char *);
	~Material() {}

	// Access Data
	void	name(const char *aname);
const	char*	name()			const { return _name; }

	void	id(const int i)		      { _id = i; }
	int	id()			const { return _id; }

	double	density()		const { return _density; }
	void	density(double d)	      { _density = d; }

	int	Z()			const { return _Z; }
	void	Z(int a)		      { _Z = a; }

	int	A()			const { return _A; }
	void	A(int a)		      { _A = a; }

	int	N()			const { return _A-_Z; }

	double	weight()		const { return _weight; }
	void	weight(double a)	      { _weight = a; }

//	dword	ambient()		const { return _amb; }
//	void	ambient(dword a)	      { _amb = a; }

	dword	diffuse()		const { return _diff; }
	dword	color()			const { return _diff; }
	void	diffuse(dword a)	      { _diff = a; }

	double	specular()		const { return _spec; }
	void	specular(double a)	      { _spec = Range(0.0,a,1.0); }

	double	shine()			const { return _shine; }
	void	shine(double a)		      { _shine = Max(0.0,a); }
//	dword	shineColor()		const { return _cshine; }
//	void	shineColor(dword b)	      { _cshine = b; }
//	void	shine(double a, dword b)      { _shine = a; _cshine = b; }

//	dword	alpha()			const { return _trans; }
//	dword	transparent()		const { return _trans; }
//	void	transparent(dword a)	      { _trans = a; }

	double	ior()			const { return _ior; }
	void	ior(double a)		      { _ior  = a; }

	double	fuzz()			const { return _fuzz; }
	void	fuzz(double a)		      { _fuzz  = Range(0.0,a,3.0); }

	bool	used()			const { return _used; }
	void	used()			      { _used = true; }

	// Flags
	int	flags()			const { return _flag; }
	void	noAntialias()		      { _flag |= MAT_NO_ANTIALIAS; }
	void	antialias()		      { _flag &= ~MAT_NO_ANTIALIAS; }
}; // Material

std::ostream& operator << (std::ostream&, const Material&);

#endif
