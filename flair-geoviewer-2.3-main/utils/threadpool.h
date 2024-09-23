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
#ifndef __THREADPOOL_H
#define __THREADPOOL_H

#include <pthread.h>
#include <semaphore.h>

class ThreadPool;
class ThreadPoolFeeder;

/* ========================= ThreadPoolWorker ========================= */
/** ThreadPoolWorker
 * A work unit for the TreadPool
 */
class ThreadPoolWorker {
protected:
	ThreadPoolFeeder*	_feeder;	// to recover easier the information...
public:
	ThreadPoolWorker(ThreadPoolFeeder* f=NULL) : _feeder(f) {}
	~ThreadPoolWorker() {}		// Destroy

virtual	void	operator()() = 0;		// Perform the work
virtual	void	feeder(ThreadPoolFeeder* f)	{ _feeder = f; }
virtual	ThreadPoolFeeder* feeder()		{ return _feeder; }
}; // ThreadPoolWorker

/* ========================= ThreadPoolFeeder ========================= */
/** ThreadPoolFeeder
 * Worker provider unit for the threadpool.
 * It has to be subclassed to send work to the threadpool
 */
class ThreadPoolFeeder {
protected:
	ThreadPool&	pool;

public:
	ThreadPoolFeeder(ThreadPool& p) : pool(p) {}
	virtual ~ThreadPoolFeeder() {}
	virtual	ThreadPoolWorker* feed(int threadId) = 0;
}; // ThreadPoolFeeder

/* ========================= _ThreadPoolData ========================== */
/** Internal structure for linking the running threads with the threadpool */
struct _ThreadData {
	int		id;
	pthread_t	thread;
	ThreadPool*	pool;
}; // _ThreadData

/* ============================ ThreadPool ============================ */
/** ThreadPool */
class ThreadPool {
private:
	int	_nthreads;			// Number of threads
	int	activeJobs;			// Number of active jobs
	bool	requestThreadExit;		// request threads to end
	bool	stopped;			// flag exited normally or stopped

	_ThreadData*		threads;	// Thread data
	ThreadPoolFeeder*	feeder;		// Feeder providing workers
	pthread_mutex_t		mutex;
	pthread_cond_t		condFeeder;
	pthread_cond_t		condEnd;

public:
	ThreadPool();
	~ThreadPool();

	int	nthreads()	const	{ return _nthreads; }

	bool	init(int n);
	void	end();

	bool	execute(ThreadPoolFeeder* f);
	void	stop();

	ThreadPoolWorker* fetchWork(int thread_id);

	bool	isRunning()	{ return feeder!=NULL || activeJobs>0; }
	void	reset()		{ stopped = false; }

static	int	num_cores()	{ return sysconf( _SC_NPROCESSORS_ONLN ); }

private:
	int	lock()		{ return pthread_mutex_lock(&mutex); }
	int	unlock()	{ return pthread_mutex_unlock(&mutex); }
	void	decActiveJobs();

static	void*	threadExecute(void *param);
static	void	getTimeout(struct timespec* ts);
}; // ThreadPool

#endif
