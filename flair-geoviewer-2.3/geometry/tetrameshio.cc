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
 * Date:	11-Oct-2016
 */

#include <fstream>
#include <assert.h>
#include <sys/stat.h>

#include "token.h"
#include "tetrameshio.h"

using namespace std;

/** read a mesh from an MESH file, ASCII or binary
 * @param mesh to populate
 * @return true if everything ok
 */
bool TetraMeshIO::read(TetraMesh& mesh)
{
	char dummy[12];
	if (!open(StreamRead, false)) return false;
	stream.read(dummy, 11); dummy[11] = 0;
	if (!strcmp(dummy, "$MeshFormat"))
		return _readAscii(mesh);
	else {
		close();
		open(StreamRead, true);
		return _readBinary(mesh);
	}
} // read

/** read mesh from acsii .msh file
 * @return true on success, false otherwise
 */
bool TetraMeshIO::_readAscii(TetraMesh& mesh)
{
//	char buf[256];

	stream.seekg(0, ifstream::beg);		 // go to start

	Token token(stream, true);
	token.setAlpha('$');
	token.next();

	// Could be that we have many solids.... (not implemented at the moment)
	int datasize, filetype;
	int n;
	double version;

	if (!token.mustbe("$MeshFormat",1)) goto _ERROR;

	version = token.getNumber();
	filetype = token.getInteger();
	datasize = token.getInteger();

	cout << "Version=" << version << " Type=" << filetype << " size=" << datasize << endl;
	if (!token.mustbe("$EndMeshFormat",1)) goto _ERROR;

	if (!token.mustbe("$Nodes", 1)) goto _ERROR;
	n = token.getInteger();
	for (int i=0; i<n; i++) {
		int j = token.getInteger();
		assert(i+1 == j);
		if (i+1 != j) goto _ERROR;
		Vertex v;
		v.x = token.getNumber();
		v.y = token.getNumber();
		v.z = token.getNumber();
		mesh.add(v);
		// FIXME how to ensure that we don't have duplicates!
	}
	if (!token.mustbe("$EndNodes", 1)) goto _ERROR;
	if (!token.mustbe("$Elements", 1)) goto _ERROR;
	n = token.getInteger();
	for (int i=0; i<n; i++) {
		int j = token.getInteger();
		assert(i+1 == j);
		if (i+1 != j) goto _ERROR;
		int etype  = token.getInteger();
		int ntags  = token.getInteger();
		/*int entity =*/ token.getInteger();
		// skip remaining tags
		for (int k=1; k<ntags; k++) token.getInteger();

		if (etype==4) {
			// tetrahedron
			int a = token.getInteger()-1;
			int b = token.getInteger()-1;
			int c = token.getInteger()-1;
			int d = token.getInteger()-1;
			mesh.add(a,b,c,d);
		} else {
			token.skipLine();
		}
	}
	if (!token.mustbe("$EndElements", 1)) goto _ERROR;

#if 0
//	mesh.process();
#endif

	close();
	return true;

_ERROR:
	close();
	return false;
} // _readAscii

/** Read mesh form binary gmsh file */
/** read mesh from binary stl.file
 * @return true on success, false otherwise
 */
bool TetraMeshIO::_readBinary(TetraMesh&)
{
	close();
	return false;
} // _readBinary
