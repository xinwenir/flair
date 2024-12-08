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

#ifndef __RANDOM_H
#define __RANDOM_H

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "os.h"
#include "vector.h"

#if defined(__MACH__) || defined(CYGWIN) || defined(ANDROID)
	#define HAS_RAND48
#endif

/* =============================== Vector ============================= */
/** Random class, generate random distributions using the drand48
 * C library
 */
class Random {
private:
	long int	_seed;			// random number seed
	double		_normalSave;		// second normal number saved
#ifndef HAS_RAND48
	struct drand48_data	rnd;
#endif

public:
	Random(long int s=0) : _normalSave(-1.0) { seed(s!=0?s:time(NULL)); }

	/** set seed value */
	void	seed(long int s) {
			_seed = s;
#ifdef HAS_RAND48
			srand48(_seed);
#else
			srand48_r(_seed, &rnd);
#endif
		}

	/** @return seed value */
	long int seed()		const	{ return _seed; }

	/** @return an integer random number [0, 2^31) */
	long int integer() {
#ifdef HAS_RAND48
			return lrand48();
#else
			long int x;
			lrand48_r(&rnd, &x);
			return x;
#endif
		}

	/** @return an integer random number [0, imax) */
	long int integer(long int imax) {
#ifdef HAS_RAND48
			return lrand48()%imax;
#else
			long int x;
			lrand48_r(&rnd, &x);
			return x%imax;
#endif
		}

	/** @return a real random number [0, 1) */
	double	real() {
#ifdef HAS_RAND48
			return drand48();
#else
			double x;
			drand48_r(&rnd, &x);
			return x;
#endif
		}

	/** @return a real random number [0, max) */
	double	real(const double max)				{ return max*real(); }

	/** @return a real random number [min, max) */
	double	real(const double min, const double max)	{ return (max-min)*real() + min; }
	double	uniform(const double min, const double max)	{ return (max-min)*real() + min; }

	/** @return a real random number [0, 1) */
	double operator ()()					{ return real(); }

	/** @return a real random number [0, max) */
	double operator ()(const double max)			{ return max*real(); }

	/** @return a real random number [min, max) */
	double operator ()(const double min, const double max)	{ return (max-min)*real() + min; }

	/** @return exponential spectrum exp(-t/a) to sample the
	 *  collision distance or random decay
	 */
	double	exp(const double a)				{ return -log(real())/a; }

	void	sincos(double *s, double *c);
	double	normal();
	void	normal(double *z1, double *z2);

	void	disc(const double r, double *x, double *y);
	void	sphere(const double r, double *x, double *y, double *z);
	Vector	sphere(const double r) {
			double x, y, z;
			sphere(r, &x, &y, &z);
			return Vector(x,y,z);
		}

	Vector	vector();
}; // Random

#endif
