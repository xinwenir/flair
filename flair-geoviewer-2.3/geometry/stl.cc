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
 * Date:	08-Aug-2012
 */

#include <fstream>
#include <sys/stat.h>

#include "stl.h"
#include "token.h"

using namespace std;

/** read a mesh from an STL file, ASCII or binary
 * @param name to open
 * @param mesh to populate
 * @return true if everything ok
 */
bool STL::read(const char *name, Mesh& mesh)
{
	if (!open(name, StreamRead, true)) return false;
	stream.seekg(80, ifstream::beg); //skip header

	int nf;			// read number of faces
	stream.read(reinterpret_cast<char*>(&nf), sizeof(int));
	mesh.free();

	struct stat filestatus;
	stat(name, &filestatus);

	// only binary format contains information on the number of faces,
	// so the file size can be calculated
	if ((nf*50 + 84) == filestatus.st_size )
		return _readBinary(mesh);	// it is a binary file
	else {
		close();
		open(name, StreamRead, false);
		return _readAscii(mesh);	// it is an ASCII file
	}
} // read

/** read mesh from acsii stl.file
 * @return true on success, false otherwise
 *
 * solid name
 * _Face_Begin____________________________________
 * facet normal n1 n2 n3
 *	outer loop
 *		vertex v1 v1 v1
 *		vertex v2 v2 v2
 *		vertex v3 v3 v3
 *	endloop
 * endfacet
 * _Face_End______________________________________
 * endsolid name
 */
bool STL::_readAscii(Mesh& mesh)
{
	char buf[256];

	Token token(stream, true, true);
	token.next();

	// Could be that we have many solids.... (not implemented at the moment)
	token.mustbe("SOLID",1);
	stream.getline(buf, sizeof(buf));	// get remaining line as name
//	mesh.name(buf);
	token.next();

	while (token()==ident_tok && !strcmp(token.string(),"FACET")) {
		token.next();
		if (!token.mustbe("NORMAL",1))    goto _ERROR;
		if (!token.mustbe(number_tok, 1)) goto _ERROR;	/* skip normal */
		if (!token.mustbe(number_tok, 1)) goto _ERROR;
		if (!token.mustbe(number_tok, 1)) goto _ERROR;

		if (!token.mustbe("OUTER",1)) goto _ERROR;

		if (!token.mustbe("LOOP",1)) goto _ERROR;

		Vertex* V[3];
		for (int i=0; i<3; i++) {
			if (!token.mustbe("VERTEX",1))  goto _ERROR;
			Vertex v;
			v.x = token.getNumber();
			v.y = token.getNumber();
			v.z = token.getNumber();
			V[i] = mesh.add(v);
		}
		mesh.add(V[0],V[1],V[2]);

		if (!token.mustbe("ENDLOOP",1))  goto _ERROR;
		if (!token.mustbe("ENDFACET",1))  goto _ERROR;
	}
	if (!token.mustbe("ENDSOLID",1))  goto _ERROR;

//	mesh.process();

	close();
	return true;

_ERROR:
	close();
	return false;
} // _readAscii

/** Read mesh form binary stl.file */
/** read mesh from binary stl.file
 * @return true on success, false otherwise
 *
 * format:
 * Byte	Data_type	Description
 * 80	ASCII		Header
 * 4	unsigned int	Number of faces
 * _Face_Begin____________________________________
 * 4	float		i for normal
 * 4	float		j
 * 4	float		k
 *
 * 4	float		x for vertex 1
 * 4	float		y
 * 4	float		z
 *
 * 4	float		x for vertex 2
 * 4	float		y
 * 4	float		z
 *
 * 4	float		x for vertex 3
 * 4	float		y
 * 4	float		z
 *
 * 2	unsigned integer	Atribute byte count(not used)
 * _Face_End______________________________________
 */
bool STL::_readBinary(Mesh& mesh)
{
	stream.seekg(80, ifstream::beg);
	//read name
	int size = stream.tellg();
	char *memblock = new char [size];
	stream.seekg(0, ios::beg);
	stream.read(memblock, size);
//	mesh.name(memblock);
	delete[] memblock;
	//done reading name

	int numFaces;
	stream.read(reinterpret_cast<char*>(&numFaces), 4);
	mesh.resize(0,0,numFaces);

	for (int i=0; i<numFaces; i++) {
		float n[3];
		stream.read(reinterpret_cast<char*>(n), sizeof(n));

		float v[3*3];
		stream.read(reinterpret_cast<char*>(v), sizeof(v));

		Vertex* A = mesh.add(Vertex(v[0],v[1],v[2]));
		Vertex* B = mesh.add(Vertex(v[3],v[4],v[5]));
		Vertex* C = mesh.add(Vertex(v[6],v[7],v[8]));
		mesh.add(A,B,C);

		// skip byte count
		unsigned short byte_cnt;
		stream.read(reinterpret_cast<char*>(&byte_cnt), 2);
	}
	close();

	mesh.process();
	return true;
} // _readBinary
