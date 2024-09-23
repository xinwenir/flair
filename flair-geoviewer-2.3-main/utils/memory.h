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
 * Date:        2 Mar 2014
 */

#ifndef __MEMORY_H
#define __MEMORY_H

#include <pthread.h>
#include <new>

#define _MAGIC	0xDECAFFEE
#define _MAGIC3	0xDE
#define _MAGIC2	0xCA
#define _MAGIC1	0xFF
#define _MAGIC0	0xEE

/*
 * This file provides some debugging functions for memory allocation
 */
struct Memory {
private:
	const char*	name;
	Memory*		next;
	Memory*		prev;
	size_t		size;
	size_t		line;
	unsigned	magic;
	unsigned char	data[sizeof(unsigned)];

public:
static	void	init();
static	bool	fini();

static	void*	malloc(size_t size, const char* filename, const size_t lineno);
static	void*	realloc(void *ptr, size_t size);
static	void	free(void *ptr);

static	void	dump();
static	void	coredump();

static	void	check();
static	size_t	allocated()		{ return total; }
static	size_t	items()			{ return n; }
static	double	resident();

private:
	void	setTail() {
			data[size]   = _MAGIC0;
			data[size+1] = _MAGIC1;
			data[size+2] = _MAGIC2;
			data[size+3] = _MAGIC3;
		}
	bool	checkTail() {
			return	data[size]   == _MAGIC0 &&
				data[size+1] == _MAGIC1 &&
				data[size+2] == _MAGIC2 &&
				data[size+3] == _MAGIC3;
		}
	void	print(int count);

private:
static	void	lock()		{ pthread_mutex_lock(&mutex);   }
static	void	unlock()	{ pthread_mutex_unlock(&mutex); }
static	pthread_mutex_t	mutex;	// locking mutex when running through threadpool

static	bool	_init;
static	Memory*	head;
static	size_t	total;
static	size_t	n;
}; // Memory

inline void* operator new(size_t size, const char* filename, const size_t lineno) {
		return Memory::malloc(size, filename, lineno);
	}
inline void* operator new[](size_t size, const char* filename, const size_t lineno) {
		return Memory::malloc(size, filename, lineno);
	}

inline void  operator delete(void* ptr)   noexcept { Memory::free(ptr); }
inline void  operator delete[](void* ptr) noexcept { Memory::free(ptr); }

/* C++ -Wsized-deallocation requests to implement this */
inline void  operator delete(void* ptr, std::size_t)   noexcept { Memory::free(ptr); }
inline void  operator delete[](void* ptr, std::size_t) noexcept { Memory::free(ptr); }

#if 1
inline void* operator new(size_t size)		{ return Memory::malloc(size, "", 0); }
inline void* operator new[](size_t size)	{ return Memory::malloc(size, "[]", 0); }
#else
#undef new
//#define new	new(__FILE__,__LINE__)
#define NEW	new(__FILE__,__LINE__)
#define new	NEW
//#define DELETE	delete
#endif

#endif
