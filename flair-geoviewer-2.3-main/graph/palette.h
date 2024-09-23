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

#ifndef __PALETTE_H
#define __PALETTE_H

#define _MAXCOLORS	256

/* ============================= Palette ============================== */
class Palette {
private:
	/* color scale */
	double	_min;			/** minimum value		*/
	double	_max;			/** maximum value		*/
	double	_step;			/** stepping log or linear	*/
	int	_n;			/** number of colors in palette */
	bool	_log;			/** logarithmic color box	*/
	bool	_interpolate;		/** interpolate colors		*/
	bool	_invert;		/** invert palette		*/
	bool	_alphamin;		/** transparent below minimum   */
	bool	_alphamax;		/** transparent above maximum   */
	dword	_palette[_MAXCOLORS];	/** color palette		*/

public:
	Palette()			{ reset(); }
	void	reset();

	/* set/get */
	double	min()		const	{ return _min; }
	void	min(double a)		{ _min = a; }

	double	max()		const	{ return _max; }
	void	max(double a)		{ _max = a; }

	double	range()		const	{ return _max - _min; }

	bool	log()		const	{ return _log; }
	void	log(bool a)		{ _log = a; }

	bool	interpolate()	const	{ return _interpolate; }
	void	interpolate(bool a)	{ _interpolate = a; }

	bool	invert()	const	{ return _invert; }
	void	invert(bool a)		{ _invert = a; }

	int	size()		const	{ return _n; }
	void	size(int a)		{ _n = Min(_MAXCOLORS,a); }

	bool	alphamin()	const	{ return _alphamin; }
	void	alphamin(bool a)	{ _alphamin = a; }

	bool	alphamax()	const	{ return _alphamax; }
	void	alphamax(bool a)	{ _alphamax = a; }

	dword	operator [](const int i) const	{ return _palette[i]; }
	dword&	operator [](const int i)	{ return _palette[i]; }

	void	init();
	void	checkLimits();
	dword	first()	const		{ return _palette[0]; }
	dword	last()	const		{ return _palette[_n-1]; }
	dword	color(const double value) const;
}; // Palette

#endif
