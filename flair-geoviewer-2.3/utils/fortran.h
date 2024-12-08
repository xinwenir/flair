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
 * Date:	26-Jul-2005
 */

#ifndef __FORTRAN_H
#define __FORTRAN_H

#include <stdio.h>
#include <string.h>
#include <cstdint>

#include "os.h"

/* ============================== FortranParser ============================ */
class FortranParser {
private:
	char*	buffer;		//** buffer pointer
	char*	ptr;		//** current location pointer
	char*	ptrEnd;		//** end of buffer pointer

public:
	FortranParser()	: buffer(NULL), ptr(NULL), ptrEnd(NULL) {}
	FortranParser(void *buf, int len=0)		{ set(buf, len); }
	FortranParser(const FortranParser& p) : buffer(p.buffer), ptr(p.ptr), ptrEnd(p.ptrEnd) {}
	void	operator()(void *buf, int len=0)	{ set(buf, len); }
	void*	operator()()				{ return (void*)buffer; }
	void	set(void *buf, int len=0) {
			ptr = buffer = (char*)buf;
			ptrEnd = ptr+len;
		}
	void	reset()		{ ptr = buffer; }
	int	length() const	{ return ptr - buffer; }
	bool	more()		{ return ptr<ptrEnd; }
	char*	getPtr()	{ return ptr; }

	void	skipWord()	{ ptr += sizeof(word); }
	void	skipInt()	{ ptr += sizeof(int); }
	void	skipFloat()	{ ptr += sizeof(float); }
	void	skipDouble()	{ ptr += sizeof(double); }

	int	read(char *s, int len) {
			memcpy(s,ptr,len);
			s[len]=0;
			ptr += len;
			return strip(s);
		}
	void	write(const char *s, int size) {
			memset(ptr, ' ', size);		// initialize to spaces
			memcpy(ptr, s, strlen(s));
			ptr += size;
		}
#ifdef ALIASING
	int	readWord() {
			char num[2];
			num[0] = *ptr++; num[1] = *ptr++;
			return *reinterpret_cast<word*>(&num);
		}
	int	readInt() {
			char num[4];
			memcpy(&num,ptr,4);
			return *reinterpret_cast<int*>(&num);
		}
	float	readFloat() {
			char num[4];
			memcpy(&num,ptr,4);
			return *reinterpret_cast<float*>(&num);
		}
	double	readDouble() {
			char num[8];
			memcpy(&num,ptr,8);
			return *reinterpret_cast<double*>(&num);
		}

	void	write(const word a) {
			memcpy(ptr, &a, sizeof(a));
			ptr += sizeof(a);
		}
	void	write(const int a) {
			memcpy(ptr, &a, sizeof(a));
			ptr += sizeof(a);
		}
	void	write(const float a) {
			memcpy(ptr, &a, sizeof(a));
			ptr += sizeof(a);
		}
	void	write(const double a) {
			memcpy(ptr, &a, sizeof(a));
			ptr += sizeof(a);
		}
#else
	int	readWord()	{ word   a=*reinterpret_cast<word*>(ptr);   ptr+=sizeof(word);   return a; }
	int	readInt()	{ int    a=*reinterpret_cast<int*>(ptr);    ptr+=sizeof(int);    return a; }
	float	readFloat()	{ float  a=*reinterpret_cast<float*>(ptr);  ptr+=sizeof(float);  return a; }
	double	readDouble()	{ double a=*reinterpret_cast<double*>(ptr); ptr+=sizeof(double); return a; }

	void	write(const word a) {
			*reinterpret_cast<word*>(ptr) = a;
			ptr += sizeof(word);
		}
	void	write(const int a) {
			*reinterpret_cast<int*>(ptr) = a;
			ptr += sizeof(int);
		}
	void	write(const float a) {
			*reinterpret_cast<float*>(ptr) = a;
			ptr += sizeof(float);
		}
	void	write(const double a) {
			*reinterpret_cast<double*>(ptr) = a;
			ptr += sizeof(double);
		}
#endif

static	int	strip(char *line);
}; // FortranParser

/* ================================ FortranFile ============================= */
class FortranFile {
public:
	FILE	*handle;	// fortran handle

public:
	FortranFile()		: handle(NULL)	{}
	FortranFile(const char *filename, const char *mode)	{ handle = fopen(filename, mode); }
	FortranFile(FILE *af)	: handle(af)	{}

	~FortranFile()	{ if (handle) fclose(handle); }

	operator int()	{ return handle!=NULL && !feof(handle); }

	bool	open(const char *filename, const char *mode) {
			handle = fopen(filename, mode);
			return handle!=NULL;
		}
	void	close() {
			if (handle) fclose(handle);
			handle = NULL;
		}

	int	seek(long offset, int whence=SEEK_SET)	{ return fseek(handle, offset, whence); }
	size_t	tell()					{ return ftell(handle); }
	int	getpos(fpos_t* position)		{ return fgetpos(handle, position); }
	int	setpos(fpos_t* position)		{ return fsetpos(handle, position); }

	int	readraw(void *buffer, int s, int n)	{ return fread(buffer, s, n, handle); }
	int	read(void *buffer, const int maxSize);
	uint8_t*	readBuffer(int* length);

	int	write(void *buffer, const int size);
	int	write(FortranParser& parser)		{ return write(parser(), parser.length()); }

	int	backspace();
	int	skip();
	int	blockSize();
	bool	mustBe(int length);
}; // class FortranFile

#endif
