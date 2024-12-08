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
 *
 */

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "os.h"
#include "eps.h"
#include "bmath.h"
#include "vector.h"
#include "random.h"

/** sample random position and direction on a unit sphere to generate
 * a uniform and isotropic fluence
 * @param pos	position to return
 * @param dir	direction to return
 */
void sflood(Random& random, Vector *pos, Vector *dir)
{
	// Random position direction
	double costh0 = 2.0 * random() - 1.0;
	double sinth0 = sqrt((1.0-costh0) * (1.0+costh0));
	double sinph0, cosph0;
	random.sincos(&sinph0, &cosph0);

	// These are the position vector components
	pos->set(cosph0*sinth0, sinph0*sinth0, costh0);

	// Now sample the direction vector
	double rndspe = random();
	double snth   = sqrt(rndspe);
	double csth   = sqrt(1.0 - rndspe);
	double snph, csph;
	random.sincos(&snph, &csph);

	double tr = -csth;
	double tp = -snth * snph;
	double tt = -snth * csph;

	if ( sinth0 > TOOSMALL ) {
		dir->x = pos->x * tr + ( pos->x * pos->z * tt - pos->y * tp ) / sinth0;
		dir->y = pos->y * tr + ( pos->y * pos->z * tt + pos->x * tp ) / sinth0;
		dir->z = pos->z * tr - sinth0 * tt;
	} else
		dir->set(tt, tp, pos->z * tr);
} // sflood
