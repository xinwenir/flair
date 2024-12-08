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
 * Date:	31-Mar-2010
 */

#include <string>
#include <ostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "bstring.h"

using namespace std;

/** string hash algorithm from http://www.cse.yorku.ca/~oz/hash.html
 djb2
 this algorithm (k=33) was first reported by dan bernstein many years ago in
 comp.lang.c. another version of this algorithm (now favored by bernstein) uses
 xor: hash(i) = hash(i - 1) * 33 ^ str[i]; the magic of number 33 (why it works
 better than many other constants, prime or not) has never been adequately
 explained.
*/
dword hash_djb2(const char *str)
{
	dword hash = 5381;
	int c;

	while ((c=*(const uint8_t*)str++))
		hash += (hash << 5) + c; /* hash * 33 + c */

	return hash;
} // hash_djb2

/** create a name000number format fitting inside length
 * with a name and a number with minimum digits
 * @param output	output string of size
 * @param size		size of output string
 * @param name		name as a prefix
 * @param number	number as a suffix
 * @param digits	minimum number of digits to append
 * @return output
 */
char* nameNumber(char* output, int size, const char* name, int number, int digits)
{
	strncpy(output, name, size);
	output[size-1] = '\0';
	char num[16];
	sprintf(num,"%0*d",digits,number);

	int slen = strlen(output);
	int nlen = strlen(num);

	if (slen + nlen < size-1)
		strcpy(output+slen, num);
	else {
		int pos = size-1-nlen;
		if (pos<=0)
			strcpy(output+1, num+1-pos);
		else
			strcpy(output+pos, num);
	}
	return output;
} // nameNumber

/**
bool AlmostEqual2sComplement(float A, float B, int maxUlps)
{
	// Make sure maxUlps is non-negative and small enough that the
	// default NAN won't compare as equal to anything.
	assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
	int aInt = *(int*)&A;
	// Make aInt lexicographically ordered as a twos-complement int
	if (aInt < 0) aInt = 0x80000000 - aInt;
	// Make bInt lexicographically ordered as a twos-complement int
	int bInt = *(int*)&B;
	if (bInt < 0) bInt = 0x80000000 - bInt;
	int intDiff = abs(aInt - bInt);
	if (intDiff <= maxUlps) return true;
	return false;
}
*/
