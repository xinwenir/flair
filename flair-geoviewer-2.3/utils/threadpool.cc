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
 * Based on original code from David Sinuela Pastor siu_coders@gmail.com
 * Modified:	Vasilis.Vlachoudis@cern.ch
 * Date:	15-Jan-2013
 *
 */

#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include "threadpool.h"

/* -------------------- ThreadPool ----------------------- */
/** ThreadPool */
ThreadPool::ThreadPool() :
	_nthreads(0),
	activeJobs(0),
	requestThreadExit(false),
	stopped(false),
	threads(NULL),
	feeder(NULL)
{
	pthread_mutex_init(&mutex,     NULL);
	pthread_cond_init(&condFeeder, NULL);
	pthread_cond_init(&condEnd,    NULL);
} // ThreadPool

/** ~ThreadPool */
ThreadPool::~ThreadPool()
{
	end();
	pthread_mutex_destroy(&mutex);
	pthread_cond_destroy(&condFeeder);
	pthread_cond_destroy(&condEnd);
} /* ~ThreadPool */

/** init
 * @return true on success
 *	   false if it didn't succeed in creating all threads
 */
bool ThreadPool::init(int n)
{
	if (n == _nthreads) return true;
	end();
	_nthreads = n;
	if (_nthreads==0) return true;
	threads = new _ThreadData[_nthreads];
	memset(threads, 0, _nthreads*sizeof(pthread_t));

	for (int i=0; i<_nthreads; i++) {
		threads[i].pool = this;
		threads[i].id   = i;
		int ret = pthread_create(&(threads[i].thread),
				NULL,
				ThreadPool::threadExecute,
				(void*)&(threads[i]));
		if (ret) {
			fprintf(stderr,"ERROR: cannot create thread\n");
			_nthreads = i-1;
			return false;
		}
	}
	return true;
} // init

/** end all running jobs before closing the threadpool */
void ThreadPool::end()
{
	if (!threads) return;
	requestThreadExit = true;
	stop();
	bool allFinished;
	do {
		// Send condFeeder to wakeup one thread
		lock();
		pthread_cond_broadcast(&condFeeder);
		unlock();
		// Check if all threads are finished
		allFinished = true;
		for (int i=0; i<_nthreads; i++) {
			if (threads[i].id>=0)
				allFinished = false;
		}
		// Make a micro-sleep
		usleep(1);
	} while (!allFinished);

	delete [] threads;
	threads   = NULL;
	_nthreads = 0;
	requestThreadExit = false;
} // end

/** execute a new job with the feeder f
 * @param f	feeder to feed the job
 * @return true if it was stopped by the user, false otherwise
 */
bool ThreadPool::execute(ThreadPoolFeeder* f)
{
	assert(_nthreads>0);

	// cannot execute while there is a stop signal pending
	if (stopped) return false;

	lock();
	assert(feeder==NULL && activeJobs==0);
	feeder = f;
	pthread_cond_broadcast(&condFeeder);
	while (feeder || activeJobs>0) {
		// wake up every few s to avoid deadlocks :(
		struct timespec ts;
		getTimeout(&ts);
		pthread_cond_timedwait(&condEnd, &mutex, &ts);
	}
	bool st = stopped;	// remember stop condition to return
	stopped = false;	// reset stopped to permit next execution
	unlock();
	return st;
} // execute

/** stop/kill the running job if any */
void ThreadPool::stop()
{
	lock();
	// Enable stop request even if we are not running
	stopped = true;
	if (feeder==NULL && activeJobs==0) {
		unlock();
		return;
	}
	feeder = NULL;
	pthread_cond_broadcast(&condFeeder);
	while (feeder || activeJobs>0) {
		// wake up every sec to avoid deadlocks :(
		struct timespec ts;
		getTimeout(&ts);
		pthread_cond_timedwait(&condEnd, &mutex, &ts);
	}
	unlock();
} // stop

/** decActiveJobs */
void ThreadPool::decActiveJobs()
{
	lock();
	if (activeJobs>0) activeJobs--;
	unlock();
} // decActiveJobs

/** the main loop providing work on the threads.
 */
ThreadPoolWorker* ThreadPool::fetchWork(int thread_id)
{
	/* all threads but one should be stopped at this point
	 * waiting for:
	 * - either condFeeder event to provide new workers
	 * - or next worker available
	 */
	lock();
	while (true) {
		while (feeder==NULL) {
			/* broadcast if there are no jobs running */
			if (activeJobs==0) pthread_cond_broadcast(&condEnd);
			pthread_cond_wait(&condFeeder, &mutex);
			if (requestThreadExit) {
				unlock();
				return NULL;
			}
		}

		ThreadPoolWorker* worker = feeder? feeder->feed(thread_id) : NULL;
		if (worker) {
			activeJobs++;
			unlock();
			return worker;
		}

		// set feeder to null
		feeder = NULL;
	}
} // fetchWork

/** the work function, waiting for the pool to provide a worker to run */
void* ThreadPool::threadExecute(void *param)
{
	_ThreadData *threaddata = (_ThreadData*)param;
	ThreadPoolWorker* worker;

	while ((worker=threaddata->pool->fetchWork(threaddata->id))) {
		(*worker)();
		threaddata->pool->decActiveJobs();
	}
	threaddata->id = -1;	// mark that is dead
	return NULL;		// never comes here but keeps compiler happy
} // threadExecute

/** getTimeout */
void ThreadPool::getTimeout(struct timespec* ts)
{
	// Solution coming from https://gist.github.com/jbenet/1087739
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
	clock_serv_t cclock;
	mach_timespec_t mts;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &mts);
	mach_port_deallocate(mach_task_self(), cclock);
	ts->tv_sec  = mts.tv_sec;
	ts->tv_nsec = mts.tv_nsec;
#else
	clock_gettime(CLOCK_REALTIME, ts);
#endif
	ts->tv_nsec += 100000000L;	// increase by 100ms
	if (ts->tv_nsec>1000000000L) {	// overflow >1s = 1e9ns
		ts->tv_nsec -= 1000000000L;
		ts->tv_sec++;
	}
} // getTimeout
