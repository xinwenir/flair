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

#include <stdlib.h>
#include <string.h>

#include "material.h"

using namespace std;

/* -------------- Material ---------------- */
Material::Material()
{
	_name[0] = 0;
	_density = 0.0;
	_Z       = 0;
	_A       = 0;
	_weight  = 0.0;

	_diff   = 0x000000;
//	_amb    = 0x000000;
//	_cshine = 0x000000;
//	_trans  = 0x000000;
	_spec   = 0.0;
	_shine  = 0.0;
	_ior    = 1.0;
	_fuzz   = 0.0;
	_flag   = MAT_CACHE;
	_used   = false;
} // Material

/* -------------- Material ---------------- */
Material::Material(const char *n)
{
	name(n);
	_density = 0.0;
	_Z       = 0;
	_A       = 0;
	_weight  = 0.0;

	_diff   = 0x000000;
//	_amb    = 0x000000;
//	_cshine = 0x000000;
//	_trans  = 0x000000;
	_spec   = 0.0;
	_shine  = 0.0;
	_ior    = 1.0;
	_fuzz   = 0.0;
	_flag   = MAT_CACHE;
	_used   = false;
} // Material

/** name - set name of body */
void Material::name(const char *aname)
{
	strncpy(_name, aname, sizeof(_name));
	_name[sizeof(_name)-1] = 0;
} // name

/* -------------- operator << ---------------- */
ostream& operator << (ostream& os, const Material& mat)
{
	os << "material \"" << mat.name() << "\" {" << endl;
	os << "\tdiffuse=" << mat.diffuse() << endl;
	os << "\tspecular=" << mat.specular() << endl;
//	os << "\tambient=" << mat.ambient() << endl;
//	os << "\ttransparent=" << mat.transparent() << endl;
//	if (RED(mat.shineColor())>0)
//		os << "\tshine=" << mat.shineColor() << endl;
//	else
//		os << "\tshine=" << mat.shineExp() << endl;
	os << "\tior=" << mat.ior() << endl;
	os << "\tfuzz=" << mat.fuzz() << endl;
	if (mat.flags() & MAT_CACHE)
		os << "\tCACHE" << endl;
	if (mat.flags() & MAT_NO_ANTIALIAS)
		os << "\tNO_ANTIALIAS" << endl;
	os << '}';
	return os;
} /* operator << */
