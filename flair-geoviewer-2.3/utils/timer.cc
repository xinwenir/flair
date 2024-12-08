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

#include <sys/time.h>
#include <unistd.h>

#include "os.h"
#include "timer.h"

/* ------------ start -------------- */
void Timer::start()
{
	struct timeval tv;
	struct timezone tz;

	gettimeofday(&tv,&tz);
	startRealTime   = (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
	startSystemTime = (double)clock() / (double)CLOCKS_PER_SEC;

	active = true;
	_laps  = 0;
} // start

/* ------------ lap -------------- */
void Timer::lap()
{
	struct timeval tv;
	struct timezone tz;

	gettimeofday(&tv,&tz);
	stopRealTime   = (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
	stopSystemTime = (double)clock() / (double)CLOCKS_PER_SEC;
	_laps++;
} // lap

/* ------------ stop -------------- */
void Timer::stop()
{
	struct timeval tv;
	struct timezone tz;

	gettimeofday(&tv,&tz);
	stopRealTime = (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
	stopSystemTime = (double)clock() / (double)CLOCKS_PER_SEC;

	active = false;
} // stop

/* ------------ resume -------------- */
void Timer::resume()
{
	struct timeval tv;
	struct timezone tz;

	double rt = realTime();		// remember real and system time
	double st = systemTime();

	gettimeofday(&tv,&tz);
	startRealTime   = (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
	startSystemTime = (double)clock() / (double)CLOCKS_PER_SEC;

	startRealTime   -= rt;		// remove it from the start time
	startSystemTime -= st;

	active = true;
} // resume

/* ------------ delayRealTime -------------- */
void Timer::delayRealTime(double t)
{
	struct timeval tv;
	struct timezone tz;
	double	ela;

	start();
	do {
		gettimeofday(&tv,&tz);
		ela = (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
	} while (ela<t);
	active = false;
} // delayRealTime

/* ------------ delaySystemTime -------------- */
void Timer::delaySystemTime(double t)
{
	struct timeval tv;
	struct timezone tz;
	double	ela;

	start();
	do {
		gettimeofday(&tv,&tz);
		ela = (double)clock() / (double)CLOCKS_PER_SEC;
	} while (ela<t);
	active = false;
} // delaySystemTime
