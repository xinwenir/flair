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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>

#include <string>
#include <ostream>
#include <fstream>

#include "memory.h"

typedef unsigned char byte;
typedef unsigned dword;

Memory*	Memory::head  = NULL;
size_t	Memory::total = 0L;
bool	Memory::_init = false;
size_t	Memory::n     = 0;
pthread_mutex_t	Memory::mutex;

/** init */
void Memory::init()
{
	pthread_mutex_init(&mutex, NULL);
	_init = true;
} // init

/** fini */
bool Memory::fini()
{
	pthread_mutex_destroy(&mutex);
	if (allocated()) {
		dump();
		coredump();
		return true;
	}
	return false;
} // fini

/** malloc */
void* Memory::malloc(size_t size, const char* filename, const size_t lineno)
{
	if (!_init) init();
	Memory* mem = (Memory*)::malloc(sizeof(Memory)+size);
	if (mem) {
		/* Create the memory header */
		mem->magic = _MAGIC;
		mem->name  = filename;
		mem->line  = lineno;
		mem->size  = size;

		lock();
		mem->next  = NULL;
		mem->prev  = head;
		if (head)
			head->next = mem;
		head = mem;
		total += size;
		n++;
		unlock();

		/* Mark also the END of data */
//		*(dword*)(mem->data+mem->size) = _MAGIC;
		mem->setTail();

#ifdef _DEBUG
//	if (size==144) {
//		printf("Debug 144\n");
////		coredump();
//	}
#endif

		return (void *)(mem->data);
	} else {
		fprintf(stderr,"Not enough memory to allocate object %s:%zu size=%zu\n",
				filename,lineno,size);
		return NULL;
	}
} // malloc

/** realloc */
void* Memory::realloc(void *ptr, size_t size)
{
	/* find our header */
	Memory* mem = (Memory *)((char *)ptr - (sizeof(Memory)-sizeof(dword)));

	/* check if the memory is valid */
	if (mem->magic != _MAGIC) {
		fprintf(stderr,"Memory::realloc: PREFIX Magic number doesn't match of object %p!\n",ptr);
		dump();
		coredump();
	}

	if (!mem->checkTail()) {
		fprintf(stderr,"Memory::realloc: SUFFIX Magic number doesn't match of object %p!\n",ptr);
		dump();
		coredump();
	}

	lock();
	total -= mem->size;
	n--;
	bool ishead = (mem==head);
	mem = (Memory *)::realloc(mem,size+sizeof(Memory));
	if (mem==NULL) {
		unlock();
		fprintf(stderr,"Not enough memory to allocate object %s:%zu size=%zu\n",
				mem->name,mem->line,size);
		return NULL;
	}

	if (ishead) head = mem;
	mem->size = size;
	total += size;
	n++;

	Memory* other = mem->prev;
	if (other)	other->next = mem;
	other = mem->next;
	if (other)	other->prev = mem;
	unlock();

	/* Mark also the new END of data */
	//*(dword *)(mem->data+mem->size) = _MAGIC;
	mem->setTail();

	return (void *)(mem->data);
} // realloc

/** free */
void Memory::free(void *ptr)
{
	/* find our header */
	Memory* mem = (Memory *)((char *)ptr - (sizeof(Memory)-sizeof(dword)));

	if (mem->magic != _MAGIC) {
		fprintf(stderr,"Memory::free: PREFIX Magic number doesn't match of object %p!\n",ptr);
		dump();
		coredump();
	}
	if (!mem->checkTail()) {
		fprintf(stderr,"Memory::free: SUFFIX Magic number doesn't match!\n");
		dump();
		coredump();
	}

	/* Remove the _MAGIC number, just to catch invalid entries */
	mem->magic = 0L;

	lock();
	Memory* mem_prev = mem->prev;
	Memory* mem_next = mem->next;
	total -= mem->size;
	n--;
	bool ishead = (mem==head);

	::free(mem);

	if (mem_next) mem_next->prev = mem_prev;
	if (mem_prev) mem_prev->next = mem_next;
	if (ishead) head = mem_prev;
	unlock();
} // free

/** print */
void Memory::print(int count)
{
	fputs((magic==_MAGIC)?"  ":"??",stderr);

	fprintf(stderr,"%3d %5zu %p %s:%zu\t\"",
		count, size, data, name, line);
	for (int i=0; i<10; i++)
		fprintf(stderr,"%c",
			isprint(data[i])? data[i]: '.');
	fprintf(stderr,"\" ");
	for (int i=0; i<10; i++)
		fprintf(stderr,"%02X ",data[i]);
	fprintf(stderr,"\n");
} // print

/** dump */
void Memory::dump()
{
	fprintf(stderr,"\nMemory dump: items: %zu  mem: %zu b\n",items(), allocated());

#if 0
	int y = 0;
#endif
	int count = 0;

	lock();
	Memory* mem = head;
	while (mem) {
		mem->print(count++);
		mem = mem->prev;
#if 0
		if (++y==15) {
			if (getchar()=='q') exit(0);
			y = 0;
		}
#endif
	}
	unlock();
	fprintf(stderr,"\n");
} // dump

/** coredump */
void Memory::coredump()
{
	raise(SIGSEGV);
} // coredump

/** check */
void Memory::check()
{
	Memory*	mem;
	int	i=0;

	lock();
	for (mem=head; mem; mem = mem->prev,i++) {
		if (mem->magic != _MAGIC) {
			fprintf(stderr,"PREFIX Magic number doesn't match! ID=%d\n",i);
			mem->print(i);
			dump();
			coredump();
		}
		if (!mem->checkTail()) {
			fprintf(stderr,"SUFFIX Magic number doesn't match! ID=%d\n",i);
			mem->print(i);
			dump();
			coredump();
		}
	}
	unlock();
} // check

/** return process memory as reported from the system */
double Memory::resident()
{
	//double vm_usage     = 0.0;
	double resident_set = 0.0;

	// the two fields we want
	size_t vsize;
	long rss;
	std::string ignore;
	std::ifstream ifs("/proc/self/stat", std::ios_base::in);
	ifs >> ignore >> ignore >> ignore >> ignore >> ignore
	    >> ignore >> ignore >> ignore >> ignore >> ignore
	    >> ignore >> ignore >> ignore >> ignore >> ignore
	    >> ignore >> ignore >> ignore >> ignore >> ignore
		>> ignore >> ignore >> vsize >> rss;

	// in case x86-64 is configured to use 2MB pages
	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
	//vm_usage = vsize / 1024.0;
	resident_set = rss * page_size_kb;

	return resident_set;
} // resident
