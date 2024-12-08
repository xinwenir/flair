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

#include <stdio.h>
#include <string.h>

#include "fortran.h"

/** strip - strip a fortran character block from trailing spaces
 * @param line	fortran character string to be stripped
 */
int FortranParser::strip(char *line)
{
	int length = (int)strlen(line)-1;
	while (length>=0 && line[length]==' ')
		length--;
	line[length+1] = 0;
	return length;
} // strip

/** read a fortran block
 * @param buffer	buffer where data should be stored
 * @param maxSize	maximum length of buffer
 * @return length of data read
 */
int FortranFile::read(void *buffer, const int maxSize)
{
	int length, length2;
	if (fread(&length,sizeof(length),1,handle)<=0)
		return -1;
	if (length>maxSize) {
		fseek(handle,sizeof(length),SEEK_CUR);
		return 0;
	}
	if (fread(buffer,length,1,handle)<=0)
		return -1;
	if (fread(&length2,sizeof(length2),1,handle)<=0)
		return -1;
	if (length!=length2)
		return -1;
	return length;
} // read

/** allocate memory (with new) and read a fortran record
 * WARNING: The calling procedure is responsible for freeing the memory!!
 * @param size	return size of the buffer allocated
 */
uint8_t* FortranFile::readBuffer(int* length)
{
	if (fread(length,sizeof(*length),1,handle)<=0)
		return NULL;

	uint8_t* buffer = new uint8_t[*length];
	if (buffer==NULL) return NULL;

	if (fread(buffer,*length,1,handle)<=0) {
		*length = -1;
		delete [] buffer;
		return NULL;
	}

	int length2;
	if (fread(&length2,sizeof(length2),1,handle)<=0) {
		*length = -1;
		delete [] buffer;
		return NULL;
	}
	if (*length!=length2) {
		*length = -1;
		delete [] buffer;
		return NULL;
	}
	return buffer;
} // readBuffer

/** write
 * @param buffer	write a buffer into file
 * @param size		length of data to write
 * @return length of data written (normally size)
 * */
int FortranFile::write(void *buffer, const int size)
{
	if (fwrite(&size, sizeof(size), 1, handle) != 1) return -1;
	if (fwrite(buffer, size, 1, handle) != 1) return -1;
	if (fwrite(&size, sizeof(size), 1, handle) != 1) return -1;
	return size;
} // write

/** skip a fortran block
 * @return length of data skipped
 */
int FortranFile::skip()
{
	int length, length2;
	if (fread(&length,sizeof(length),1,handle)<=0)
		return -1;
	fseek(handle, length, SEEK_CUR);
	if (fread(&length2,sizeof(length2),1,handle)<=0)
		return -1;
	if (length!=length2)
		return -1;
	return length;
} // skip

/** backspace a fortran block
 * @return length of data skipped
 */
int FortranFile::backspace()
{
	int length;
	if (ftell(handle)<=8) return 0;
	fseek(handle, -4, SEEK_CUR);
	if (fread(&length,sizeof(length),1,handle)<=0)
		return -1;
	fseek(handle, -(length+8), SEEK_CUR);
	return length;
} // backspace

/** read fortran block size
 * @return size of block data
 */
int FortranFile::blockSize()
{
	int length;
	if (fread(&length,sizeof(length),1,handle)<=0)
		return -1;
	fseek(handle, -4, SEEK_CUR);
	return length;
} // blockSize

/** fortran size mustbe
 * @param length size of block
 * @return size of block data
 */
bool FortranFile::mustBe(int length)
{
	int length2;
	if (fread(&length2,sizeof(length2),1,handle)<=0)
		return false;
	return length==length2;
} // mustBe
