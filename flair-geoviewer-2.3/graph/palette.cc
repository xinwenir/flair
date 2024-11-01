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

#include <math.h>

#include <iostream>

#include "os.h"
#include "eps.h"
#include "color.h"
#include "palette.h"

/** reset */
void Palette::reset()
{
	// color scale
	_log         = true;
	_interpolate = false;
	_min         = -4;
	_max         =  4;
	_invert      = false;
	_alphamin    = false;
	_alphamax    = false;
	_n           = 2;
	_palette[0]  = 0;
	_palette[1]  = 0XFFFFFF;
	init();
} // reset

/** init */
void Palette::init()
{
	_n    = Range(2, _n, _MAXCOLORS);
	_step = (_max - _min) / (double)(_interpolate? _n-1 : _n);
} // init

/** checkLimits */
void Palette::checkLimits()
{
	if (_min > _max) {
		Swap(_min, _max);
		_invert = !_invert;
		if (_step < 0.0) _step = -_step;
	}
	if (_step <= SMALL) {
		_step = 1.0;
		_max = _min + 4.0;
	}
} // checkLimits

/** color */
dword Palette::color(const double value) const
{
	if (value<=_min) {
		if (_alphamin) return COLOR_TRANSPARENT;
		return _invert? _palette[_n-1] : _palette[0];
	}
	if (value>=_max) {
		if (_alphamax) return COLOR_TRANSPARENT;
		return _invert? _palette[0] : _palette[_n-1];
	}

//	std::cout << value << std::endl;
	double div = (value-_min) / _step;
	int ic = (int)div;
	if (_invert) {
		div = (double)_n - div;
		ic  = _n - ic - 1;
		if (_interpolate) {
			div--;
			ic--;
		}
	}
	if (_interpolate) {
		int x  = (int)(255.0*(div - (double)ic));
		int x1 = 256 - x;

		Color32 lo, hi, col;
		if (_invert) {
			Swap(x,x1);
			lo.val = _palette[ic+1];
			hi.val = _palette[ic];
		} else {
			lo.val = _palette[ic];
			hi.val = _palette[ic+1];
		}

		col.val = 0;
		col.rgb.red   = (uint8_t)(((int)lo.rgb.red  *x1 + (int)hi.rgb.red  *x)>>8);
		col.rgb.green = (uint8_t)(((int)lo.rgb.green*x1 + (int)hi.rgb.green*x)>>8);
		col.rgb.blue  = (uint8_t)(((int)lo.rgb.blue *x1 + (int)hi.rgb.blue *x)>>8);

		return col.val;
	} else
		return _palette[ic];
} // color
