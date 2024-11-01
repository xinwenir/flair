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
 */
#ifndef __GLMESH_H
#define __GLMESH_H

#include <GL/gl.h>

#include "os.h"
#include "mesh.h"
#include "color.h"
//#include "geometry.h"
//#include "globject.h"

/* --- GLMesh --- */
class GLMesh { //: GLObject {
private:
	Mesh*	mesh;		/** Mesh to draw		*/
	bool	fill;		/** Fill Faces			*/

public:
static	bool	drawNormals;	/** Draw Normals		*/
static	bool	drawVertices;	/** Draw vertices		*/
static	bool	drawEdges;	/** Draw Edges			*/
static	bool	drawFaces;	/** Draw Faces			*/

static	double	normalLength;	/** Normal vector length	*/
static	double	foldCos;	/** Cosine that face is fold	*/

static	Color	normalColor;	/** Normal vector color		*/
static	Color	edgeColor;	/** Edge color			*/
static	Color	vertexColor;	/** Vertex color		*/

public:
	GLMesh(Mesh *amesh=NULL) : //: GLObject(),
		fill(TRUE)
		{ mesh = amesh; } //setColor(0.8, 0.8, 0.8); }

	virtual ~GLMesh() {}

	void setMesh(Mesh *amesh)
		{ mesh = amesh; }

	virtual void draw(GLenum mode);
}; // class GLMesh

#endif
