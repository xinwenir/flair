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
 * Date:	04-Feb-2002
 */

#ifndef __TIMER_H
#define __TIMER_H

#include <iostream>

#include "os.h"

//////////////// Timer /////////////////
/** A cronometer class for both system and real time */
class Timer {
protected:
	bool	active;
	int	_laps;
	double	startRealTime,
		stopRealTime;
	double	startSystemTime,
		stopSystemTime;
public:
	Timer(bool on=false) : active(false), _laps(0),
			startRealTime(0.0),
			stopRealTime(0.0),
			startSystemTime(0.0),
			stopSystemTime(0.0)
		{ if (on) start(); }
	~Timer() {}

	/** reset timer */
	void reset()	{ active = false; }

	/** start measuring */
	void start();

	/** restart timer */
	void restart()	{ start(); }

	/** mark a lap */
	void lap();

	/** stop timer */
	void stop();

	/** resume timer */
	void resume();

	/** @return true if timer is active */
	bool working()	const	{ return active; }

	/** @return number of laps */
	int laps()	const	{ return _laps; }

	/** @return real time in seconds */
	double elapsed()	{ lap(); return realTime(); }

	/** @return real time in seconds */
	double realTime() const
		{ return (stopRealTime - startRealTime); }

	/** @return system time in seconds */
	double systemTime() const
		{ return (stopSystemTime - startSystemTime); }

	/** wait until real time t seconds is passed
	 * @param t time to wait
	 */
	void delayRealTime(double t);

	/** wait until system time t seconds is passed
	 * @param t time to wait
	 */
	void delaySystemTime(double t);
}; // Timer

/** Output to a stream the timer */
inline std::ostream& operator << (std::ostream& os, Timer& tm)
{
	if (tm.working()) tm.elapsed();
	os << "Real: " << tm.realTime();
	os << "  System: " << tm.systemTime() << " (sec)";
	return os;
}

#endif
