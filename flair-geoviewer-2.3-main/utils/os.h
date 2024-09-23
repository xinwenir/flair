/*
 * $Id$
 *
 * Copyright 1996-2000+ Vasilis Vlachoudis. All Rights Reserved.
 *
 * This software is the proprietary information of Vasilis Vlachoudis
 * Use is subject to license terms.
 */

#ifndef __OS_H
#define __OS_H

/* ----------- Definitions depending on the OS system ---------- */

/*
 * The following definitions should be supplied from the
 *
 * Operating system
 *	WCE	- For a PalmTop version Windows CE
 *	MSDOS	- Compiled for a MSDOS operating system
 *	MSWIN	- Compiled for a WIN32 system
 *	XWIN	- For a unix with X11
 *	CURSES	- For a unix with Curses
 *
 * Compiler
 *	__BORLANDC__	For the Borland C++ compiler
 *	GCC	- GNU C compiler
 *	MSC_VER	- Microsoft C
 *
 * General
 *	OS	- Operating system name
 *	SHELL	- Command Shell
 *	PATHSEP	- Path separator
 *	FILESEP - Directory separator
 *	CPU_ALIGN	- Align the code in 4-bytes (used for RISC machines)
 *	CPU_ENDIAN	- LITTLE_ENDIAN or BIG_ENDIAN cpu type
 *
 *	USE_UNICODE	- Use the unicode functions instead
 *	UNICODE
 *	HAS_STDIO
 *	HAS_STRING
 *	HAS_CTYPE
 *	HAS_XTOY
 *	HAS_SIGNAL
 */

#ifndef BIG_ENDIAN
	#define LITTLE_ENDIAN	1234
	#define BIG_ENDIAN	4321
#endif

/* ======== Operating system specifics ========= */
#if defined(WCE)
	#include <windows.h>

	#define	OS		"WINCE"
	#define	SHELL		"SHELL"
	#define	FILESEP		'\\'
	#define	PATHSEP		';'
	#define	CPU_ALIGN	1

	#pragma warning(disable : 4018)	/* signed unsigned comparison warning */

	#ifndef CPU_ENDIAN
		#ifdef _WIN32_WCE_EMULATION
			#define	CPU_ENDIAN LITTLE_ENDIAN
		#else
			#define	CPU_ENDIAN BIG_ENDIAN
		#endif
	#endif

#elif defined(MSDOS)
	#define	OS		"MSDOS"
	#define	SHELL		"COMSPEC"
	#define	FILESEP		'\\'
	#define	PATHSEP		';'

	#define	HAS_STDIO
	#define	HAS_STRING	/* For BorlandC use special Strings */
	#define	HAS_CTYPE
	#define	HAS_XTOY
	#define	HAS_SIGNAL

	#ifndef CPU_ENDIAN
		#define	CPU_ENDIAN	LITTLE_ENDIAN
	#endif

#elif defined(MSWIN)
	#include <tchar.h>
	#include <windows.h>

	#define	OS		"WIN"
	#define	SHELL		"COMSPEC"
	#define	FILESEP		'\\'
	#define	PATHSEP		';'

	#define	HAS_STDIO
	#define	HAS_STRING
	#define	HAS_CTYPE
	#define	HAS_XTOY

	#define HAS_ABS

	#pragma warning(disable : 4018)	/* signed unsigned comparison warning */

	#ifndef CPU_ENDIAN
		#define	CPU_ENDIAN	LITTLE_ENDIAN
	#endif

#elif	defined(XWIN) || defined(CURSES) || \
	defined(UNIX) || defined(Darwin) || defined(CYGWIN) || \
	defined(unix) || defined(__unix__) || defined(__unix) || \
	defined(ANDROID) || defined(__APPLE__)

	#if defined(UNICODE)
		#include <wchar.h>
	#endif

	#include <unistd.h>

	#define	OS		"UNIX"
	#define	SHELL		"bash"
	#define	FILESEP		'/'
	#define	PATHSEP		':'

	#define HAS_STRING
	#define HAS_XTOY
	#define HAS_CTYPE
	#define HAS_STDIO

	#define HAS_ABS

	#define	_INCLUDE_POSIX_SOURCE
	/* LITTLE or BIG ENDIAN depends on the CPU */
	#ifndef CPU_ENDIAN
		#if  defined(__i386__) || defined(__x86_64__) || defined(__INTEL__)
			#define	CPU_ENDIAN	LITTLE_ENDIAN
//#if _DEBUG>1
//			#include <sys/ptrace.h>
//			#define BREAKON(x)	if ((x) && ptrace(PTRACE_TRACEME,0, NULL) == -1) __asm__("int $3")
//#else
//#endif
//			#define BREAKON(x)
		#else
			#define	CPU_ENDIAN	BIG_ENDIAN
		#endif
	#endif

#else
	#error "Please define the operating system"
#endif

#ifdef MEM
#	include "memory.h"
//#	include "memtrack.h"
#endif


/* ----------- Other definitions --------------- */
#ifdef WCE
	#define	__CDECL	__cdecl
#else
	#define	__CDECL
#endif

#ifndef __BORLANDC__
	#define huge
#endif

/* --------------- NON UNICODE ----------------- */
#if !defined(TCHAR)
	#if defined(UNICODE)
		#define	TCHAR		wchar_t
		#define	LPTSTR		wchar_t*
		#define	LPCTSTR		const wchar_t*
	#else
		#define	TCHAR		char
		#define	LPTSTR		char*
		#define	LPCTSTR		const char*
	#endif
	#ifndef TEXT
		#define	TEXT(x)		(x)
	#endif
#endif

#ifndef INVALID_HANDLE_VALUE
	#define	INVALID_HANDLE_VALUE	NULL
#endif

/* ------------- type definitions -------------- */
//typedef unsigned char	byte;
typedef unsigned short	word;
typedef unsigned int	dword;
typedef unsigned long	ddword;

typedef signed char	int8;
typedef signed short	int16;
typedef signed int	int32;
typedef signed long	int64;
#ifndef bool
//	typedef int		Bool;
#ifndef __cplusplus
	typedef int             bool;
#endif
#else
//	typedef int             boolean;
#endif

/* -------- comonly used definitions ----------- */
#ifndef FALSE
	#define FALSE 0
	#define TRUE  1
#endif
#ifndef OFF
	#define OFF   0
	#define ON    1
#endif
#ifndef NO
	#define NO    0
	#define YES   1
#endif

#define MAX_REAL	1e308

#ifndef PI
	#define PI	3.141592653589793238462643383279
#endif
#define PI2		6.283185307179586476925286766559
#define PIover2		1.570796326794896619231321691640
#define PIover180	0.0174532925199433
#define PIunder180	57.2957795130823
#define SqrPI
#define SqrtPI		1.772453850905516027298167483341
#define SqrtPI2		2.506628274631000502415765284811

#define Sqrt2		1.414213562373095048801688724210
#define Sqrt3		1.732050807568877293527446341506
#define Sqrt5		2.236067977499789696409173668731
#define Sqrt6		2.449489742783178098197284074706
#define Sqrt7		2.645751311064590590501615753639
#define Sqrt12		3.464101615137754587054892683012

#define Ln2		0.693147180559945309417232121458
#define Ln10		2.302585092994045684017991454684

#define KB		1024
#define MB		KB*KB
#define GB		KB*MB

/* --- BIT definitions --- */
#define MASK1      0x0001
#define MASK2      0x0003
#define MASK3      0x0007
#define MASK4      0x000f
#define MASK5      0x001f
#define MASK6      0x003f
#define MASK7      0x007f
#define MASK8      0x00ff
#define MASK9      0x01ff
#define MASK10     0x03ff
#define MASK11     0x07ff
#define MASK12     0x0fff
#define MASK13     0x1fff
#define MASK14     0x3fff
#define MASK15     0x7fff
#define MASK16     0xffff
#define MASK17     0x1ffffL
#define MASK18     0x3ffffL
#define MASK19     0x7ffffL
#define MASK20     0xfffffL
#define MASK21     0x1fffffL
#define MASK22     0x3fffffL
#define MASK23     0x7fffffL
#define MASK24     0xffffffL
#define MASK25     0x1ffffffL
#define MASK26     0x3ffffffL
#define MASK27     0x7ffffffL
#define MASK28     0xfffffffL
#define MASK29     0x1fffffffL
#define MASK30     0x3fffffffL
#define MASK31     0x7fffffffL
#define MASK32     0xffffffffL

#define BIT0       0x0001
#define BIT1       0x0002
#define BIT2       0x0004
#define BIT3       0x0008
#define BIT4       0x0010
#define BIT5       0x0020
#define BIT6       0x0040
#define BIT7       0x0080
#define BIT8       0x0100
#define BIT9       0x0200
#define BIT10      0x0400
#define BIT11      0x0800
#define BIT12      0x1000
#define BIT13      0x2000
#define BIT14      0x4000
#define BIT15      0x8000
#define BIT16      0x10000L
#define BIT17      0x20000L
#define BIT18      0x40000L
#define BIT19      0x80000L
#define BIT20      0x100000L
#define BIT21      0x200000L
#define BIT22      0x400000L
#define BIT23      0x800000L
#define BIT24      0x1000000L
#define BIT25      0x2000000L
#define BIT26      0x4000000L
#define BIT27      0x8000000L
#define BIT28      0x10000000L
#define BIT29      0x20000000L
#define BIT30      0x40000000L
#define BIT31      0x80000000L

/* --- bytes, words --- */
#ifndef LOBYTE
	#define LOBYTE(x)	((byte)(x))
	#define HIBYTE(x)	((byte)(((word)(x)>>8) & MASK8))
	#define LOWORD(x)	((word)(x))
	#define HIWORD(x)	((word)(((dword)(x)>>16) & MASK16))
	#define MAKEWORD(l,h)	((word)(((byte)(l)) | ((word)((byte)(h))) << 8))
	#define MAKEDWORD(l,h)	((dword)(((word)(l)) | ((dword)((word)(h))) << 16))
#endif

/* ------------------ comonly used macros -------------------- */
#define BOOL(a)		((a)?1:0)
#define BIT(n)		(1L << (n))
#define GETBIT(a,n)	BOOL((a) & BIT(n))
#define SETBIT(n,i)	((n) |= BIT(i))
#define CLRBIT(n,i)	((n) &= ~BIT(i))
#define TESTBIT(n,i)	(((n) & BIT(i)) != 0)

#define SHIFT(a,x)	((x>=0) ? ((a) << (x)) : ((a) >> (x)))
#define SUBBIT(a,p,x)	((dword)(((dword) a) << (32-p-x)) >> (32-x))

#define BISDIGIT(c)	((c)>='0' && (c)<='9')
#define	BISPRINT(c)	((c)>=' ' && (c)<='~' && (c)!='\'')
#define BISSPACE(c)	((c==0x09) || (c==0x0D) || (c==0x20))
#define HEXVAL(x)	(((x)>='A')?((((x)>='a')?((x)&(0xdf)):(x))-'A'+10):((x)-'0'))
#define CTL(a)		((a) & 0x1f)
#ifndef ESC
#define	ESC		((char)27)
#endif
#ifndef TAB
#define	TAB		((char)9)
#endif

#define	MOD(a,b)	((a)-((int)((a)/(b)))*(b))

#define SWAP(a,b)	a ^= b ^= a ^= b;

#define SIZE(p)		(sizeof(p) / sizeof(p[0]))
#define ABS(a)		(((a)<0)?-(a):(a))

#ifndef MAX
	#define MAX(a,b)	(((a)>(b))?(a):(b))
	#define MIN(a,b)	(((a)<(b))?(a):(b))
#endif

#define MAX3(a,b,c)	MAX(MAX(a,b),c)
#define MAX4(a,b,c,d)	MAX(MAX3(a,b,c),d)

#define MIN3(a,b,c)	MIN(MIN(a,b),c)
#define MIN4(a,b,c,d)	MIN(MIN3(a,b,c),d)

#define RANGE(a,x,b)	((x)<(a)?(a): ((x)>(b)?(b):(x)))
#define IN_RANGE(a,x,b)	(((a) <= (x)) && ((x) <= (b)))

#define ODD(n)		((n)&1)
#define EVEN(n)		(!((n)&1))

#define TRUNC(a)	((a)>=0? (int)(a) : (int)((a)-1))
#define ROUND(a)	((a)>=0? (int)((a)+0.5): -(int)(0.5-(a)))
#define CEILING(a)	((a)==(int)(a)? (int)(a): \
			(a)>0? (int)(1+(int)(a)): -(int)(1+(int)(-(a))))
#define FRAC(x)		((x)-(float)((int)(x)))
#define SIGN(x)		((x)>0? 1:(((x)==0)?0:-1))
#define SQR(x)		((x)*(x))

#define RAD(a)		((a)*PIover180)
#define DEG(a)		((a)*PIunder180)

#define COSD(a)		cos(RAD(a))
#define SIND(a)		sin(RAD(a))
#define TAND(a)		tan(RAD(a))

#define LOG(x)		(log(x)*ONEoverLn10)
#define EXP10(x)	exp((x)*Ln10)

#define RANDOM()	((double)rand()/(double)RAND_MAX)

/* -------------- Enable code ------------------ */
#ifdef _DUMP
	#define DUMP(X)		X
#else
	#define DUMP(X)
#endif

#ifdef EXPERIMENTAL
	#define EXP(X)		X
#else
	#define EXP(X)
#endif

#ifdef _STAT
	#define STATS(X)	X
#else
	#define STATS(X)
#endif

#if _DEBUG>1
	#define DEBUG(X)	X
#else
	#define DEBUG(X)
#endif

/* -------------------- I/O -------------------- */
#if defined(HAS_STDIO)
	/* --- Use the standard I/O --- */
	#define	STDIN		stdin
	#define	STDOUT		stdout
	#define	STDERR		stderr

	#define	FILEP		FILE*
#if !defined(FOPEN)
	#define	FOPEN		fopen
#endif
	#define	FEOF		feof
	#define	FTELL		ftell
	#define	FSEEK		fseek
	#define	FFLUSH		fflush
	#define	FCLOSE		fclose
	#define	FPUTC		fputc
	#define	FPUTS		fputs
	#define	FGETC		fgetc
#if !defined(FWRITE)
	#define	FWRITE		fwrite
#endif
	#define	PRINTF		printf
	#define	PUTS		puts
	#define	PUTCHAR		putchar

	#ifdef MSWIN
		#define	GETCWD		_getcwd
		#define	CHDIR		_chdir
	#else
		#define	GETCWD		getcwd
		#define	CHDIR		chdir
	#endif

#elif defined(USE_UNICODE)
	#define	FILEP		HANDLE
	#define	FOPEN		CreateFile
	#define	FCLOSE		CloseHandle
	#define	FWRITE(p,n,m,f)	{int tmp; WriteFile(f,p,(n)*(m),&tmp,NULL);}

	#define	FGETC		fgetc
	#define	FPUTC		fputc
	#define	FTELL(f)	SetFilePointer(f,0,0,FILE_BEGIN)
	#define	FSEEK(f,p,m)	SetFilePointer(f,p,0,FILE_BEGIN)
#else
	/* --- Use the home made I/O --- */
	#define	STDIN		NULL
	#define	STDOUT		NULL
	#define	STDERR		NULL

	#define	FILEP		BFILE*
	#define	FOPEN		Bfopen
	#define	FCLOSE		Bfclose
	#define	FEOF		Bfeof
	#define	FTELL		Bftell
	#define	FSEEK		Bfseek
	#define	FFLUSH		Bfflush
	#define	FPUTC		Bfputc
	#define	FPUTS		Bfputs
	#define	FGETC		Bfgetc

	#define	PUTS		Bputs
	#define	PUTINT		Bputint
	#define	PUTCHAR		Bputch
	#define	PRINTF		Bprintf
	#define	GETCHAR		WReadKey

	#define	GETCWD		Bgetcwd
	#define	CHDIR		Bchdir
#endif

/* ---------------- Memory Ops ------------------- */
#if defined(__BORLANDC__)	/* For the HUGE/LARGE memory model */
	#define	MEMMOVE		memmove
	#define	MEMCPY		_fmemcpy
	#define	MEMCMP		_fmemcmp
	#define	MEMCHR		_fmemchr
	#define	MEMSET		_fmemset
#else
	#define	MEMMOVE		memmove
	#define	MEMCPY		memcpy
	#define	MEMCMP		memcmp
	#define	MEMCHR		memchr
	#define	MEMSET		memset
#endif

/* ----------------- Strings --------------------- */
#if defined(UNIX) ||defined(unix) || defined(__unix__) || defined(__unix)
	#if defined(UNICODE)
		#define _tcscmp		wcscmp
		#define _tcsncmp		wcsncmp
		#define _tcscpy		wcscpy
		#define _tcslen		wcslen
	#else
		#define _tcscmp		strcmp
		#define _tcsncmp		strncmp
		#define _tcscpy		strcpy
		#define _tcslen		strlen
		#define _ttoi		atoi
	#endif
#endif

#if defined(__BORLANDC__)	/* For the HUGE/LARGE memory model */
	#define	STRCPY		_fstrcpy
	#define	STRCMP		_fstrcmp
	#define	STRNCMP		_fstrncmp
	#define	STRCAT		_fstrcat
	#define	STRCHR		_fstrchr
	#define	STRLEN		_fstrlen
	#define	STRSTR		_fstrstr
	#define	STRNCPY		_fstrncpy
	#define	SPRINTF		sprintf
	#define	MKTEMP		mktemp
#elif defined(USE_UNICODE)
	#define	STRCPY		wcscpy
	#define	STRCMP		wcscmp
	#define	STRNCMP		wcsncmp
	#define	STRNCPY		wcsncpy
	#define	STRCAT		wcscat
	#define	STRCHR		wcschr
	#define	STRLEN		wcslen
	#define	STRDUP		wcsdup
	#define	SPRINTF		swprintf
#elif defined(HAS_STRING)
	#define	STRCPY		strcpy
	#define	STRCMP		strcmp
	#define	STRNCMP		strncmp
	#define	STRNCPY		strncpy
	#define	STRCAT		strcat
	#define	STRCHR		strchr
	#define	STRLEN		strlen
	#define	STRSTR		strstr
	#define	STRUPR		strupr
	#define	STRDUP		strdup
	#define	SPRINTF		sprintf
	#define	MKTEMP		mktemp
#else
	#define	STRCPY		Bstrcpy
	#define	STRCMP		Bstrcmp
	#define	STRNCMP		Bstrncmp
	#define	STRCAT		Bstrcat
	#define	STRCHR		Bstrchr
	#define	STRLEN		Bstrlen
	#define	STRSTR		Bstrstr
	#define	STRUPR		Bstrupr
	#define	SPRINTF		Bsprintf
#endif

/* ----------------- Ctype ------------------------- */
#if defined(USE_UNICODE)
	#define	TOUPPER		towupper
	#define	TOLOWER		towlower
	#define	ISDIGIT		iswdigit
	#define	ISPRINT		iswprint
#elif defined(HAS_CTYPE)
	#define	TOUPPER		toupper
	#define	TOLOWER		tolower
	#define	ISSPACE		isspace
	#define	ISDIGIT		isdigit
	#define	ISPRINT		isprint
	#define	ISXDIGIT	isxdigit
	#define	ISALPHA		isalpha
#else
	#define	ISSPACE		Bisspace
	#define	ISDIGIT		Bisdigit
	#define	ISXDIGIT	Bisxdigit
	#define	ISALPHA		Bisalpha
#endif

/* ----------------- Conversions ------------------- */
#ifdef USE_UNICODE
	#define	STRTOL		wcstol
#elif defined(HAS_XTOY)
	#define	LTOA		ltoa
	#define	GCVT		gcvt
	#define	FCVT		fcvt
	#define	ECVT		ecvt
	#define	STRTOL		strtol
#elif defined(WCE)
	#define	LTOA		_ltoa
	#define	GCVT		_gcvt
	#define	FCVT		_fcvt
	#define	ECVT		_ecvt
#else
	#error	"ERROR HAS_XTOY: No home made conversions!"
#endif

#define CLEAR(x)	memset(&(x), 0, sizeof(x))

/* ------------------- Signal ---------------------- */
#ifndef SIGNAL
#ifdef HAS_SIGNAL
	#define	SIGNAL		signal
#else
	#define	SIGNAL		WSignal
#endif
#endif

#ifdef __cplusplus
/* --- Min --- */
template <class T>
inline const T& Min(const T& a, const T& b) {
	const T& retval = (a < b ? a : b);
	return retval;
}

/* --- Min --- */
template <class T>
inline const T& Min(const T& a, const T& b, const T& c) {
	const T& retval = Min(Min(a,b),c);
	return retval;
}

/* --- Min --- */
template <class T>
inline const T& Min(const T& a, const T& b, const T& c, const T& d) {
	const T& retval = Min(Min(a,b,c),d);
	return retval;
}

/* --- Max --- */
template <class T>
inline const T& Max(const T& a, const T& b) {
	const T& retval = (a > b ? a : b);
	return retval;
}

/* --- Max --- */
template <class T>
inline const T& Max(const T& a, const T& b, const T& c) {
	const T& retval = Max(Max(a,b),c);
	return retval;
}

/* --- Max --- */
template <class T>
inline const T& Max(const T& a, const T& b, const T& c, const T& d) {
	const T& retval = Max(Max(a,b,c),d);
	return retval;
}

/** interpolate a value [0,x,1] in the limits [a,b]
 * same as Map(x,0,1,a,b)
 */
template <class T>
inline const T& Interpolate(const T& a, const T& x, const T& b) {
	const T& retval = (b-a)*x + a;
	return retval;
}

/** remap a value from [in_min, in_max] to [out_min, out_max] */
template <class T>
inline const T& Map(const T& x, const T& in_min, const T& in_max, const T& out_min, const T& out_max) {
	const T& retval = (x-in_min) * (out_max-out_min) / (in_max-in_min) + out_min;
	return retval;
}

/* --- Range --- */
template <class T>
inline const T& Range(const T& a, const T& x, const T& b) {
	const T& retval = (x<a ? a: (x>b? b : x ));
	return retval;
}

/* --- InRange --- */
template <class T>
inline bool InRange(const T& a, const T& x, const T& b) {
	return ((a <= x) && (x <= b));
}

/* --- InRangeOpen --- */
template <class T>
inline bool InRangeOpen(const T& a, const T& x, const T& b) {
	return ((a < x) && (x < b));
}

/* --- Swap --- */
template <class T>
inline void Swap(T& a, T& b) {
	T tmp = a;
	a = b;
	b = tmp;
}

/* --- Sign --- */
template <typename T>
inline int Sign(const T val) {
	return (T(0) < val) - (val < T(0));
}

/* --- Int --- */
template <class T>
inline int Int(const T& x) {
	return x>=0? (int)x : (int)x - 1;
}

/* --- Int --- */
template <class T>
inline long long LongInt(const T& x) {
	return x>=0? (long long)x : (long long)x - 1;
}

/* --- Round --- */
inline int Round(const double& a)
{
	return Int(a+0.5);
} /* Round */

/* --- LongRound --- */
inline long long LongRound(const double& a)
{
	return LongInt(a+0.5);
} /* LongRound */

/* --- Fractional --- */
template <class T>
inline T Frac(const T& x) {
	return x - (T)Int(x);
}

/* --- RFractional --- */
template <class T>
inline T RFrac(const T& x) {
	return 1.0 - (x - (T)Int(x));
}

/* --- Sqr --- */
template <class T>
inline T Sqr(const T& x) {
	return x*x;
}

/* --- Sqrt --- *
 * safe sqrt for negative numbers
 */
template <class T>
inline T Sqrt(const T& x) {
	return x>0.0? sqrt(x): 0.0;
}

/* --- Cube --- */
template <class T>
inline T Cube(const T& x) {
	return x*x*x;
}

/* --- Abs --- */
template <class T>
inline T Abs(const T& a) {
	return (a<0) ? -a : a;
}
#ifdef HAS_ABS
	#include <math.h>
	// Partial specializations
	template <>
	inline double Abs(const double& a) {
		return fabs(a);
	}
	template <>
	inline float Abs(const float& a) {
		return fabsf(a);
	}
#endif

/* --- Eq0 --- */
inline bool Eq0(const double& n, const double &eps)
{
	return -eps<=n && n<=eps;
} /* Eq0 */

/* --- Cmp0 --- */
inline int Cmp0(const double& n, const double &eps)
{
	if (-eps<=n && n<=eps)
		return 0;

	if (n<0.0)
		return -1;

	return 1;
} /* Cmp0 */

/* --- Eq --- */
inline bool Eq(const double &a, const double &b, const double &eps)
{
	return Abs(a-b) <= eps;
} /* Eq */

/* --- EqRelative --- */
inline bool EqRelative(const double &a, const double &b, const double &eps)
{
	double v;
	double aa = Abs(a);
	double bb = Abs(b);
	double ab = a-b;
	if (ab<0.0) ab = -ab;

	if (aa>ab)
		v = aa;
	else
		v = bb;

	if (v<1.0) {
		if (ab <= eps)
			return TRUE;
	} else
	if (ab <= v*eps)
		return TRUE;

	return FALSE;
} /* Eq */

/* --- Cmp --- */
inline int Cmp(const double &a, const double &b)
{
	if (a < b)
		return -1;
	else
	if (a > b)
		return  1;
	else
		return  0;
} /* Cmp */

/* --- Cmp --- */
inline int Cmp(const double &a, const double &b, const double &eps)
{
	double v;
	double aa = Abs(a);
	double bb = Abs(b);
	double ab = a-b;
	if (ab<0.0) ab = -ab;

	if (aa>bb)
		v = aa;
	else
		v = bb;

	if (v<1.0) {
		if (ab <= eps)
			return 0;
	} else
	if (ab <= v*eps)
		return 0;

	if (a<b)
		return -1;
	return 1;
} /* Cmp */

/* --- Next --- *
 * return next index up to a maximum number n then go back to 0
 */
inline int Next(int i, const int n)
{
	if (++i == n)
		return 0;
	return i;
} /* Next */

/* --- Prev --- *
 * return prev index up to a minimum number 0 then go back to n-1
 */
inline int Prev(int i, const int n)
{
	if (--i == -1)
		return n-1;
	return i;
} /* Prev */
#endif

#endif
