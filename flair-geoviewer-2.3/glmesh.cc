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
#include <iostream>

#include "os.h"
#include "edge.h"
#include "face.h"
#include "glmesh.h"

#define FACE_NORMAL	0
#define FACE_FOLD	1
#define FACE_ORPHAN	2

double	GLMesh::foldCos       = 0.9;
double	GLMesh::normalLength  = 1.0;

bool	GLMesh::drawNormals  = true;
bool	GLMesh::drawVertices = false;
bool	GLMesh::drawEdges    = true;
bool	GLMesh::drawFaces    = true;

Color	GLMesh::normalColor(34, 221, 221);
Color	GLMesh::vertexColor(255, 187, 255);
Color	GLMesh::edgeColor(.7, .3, 0.);

using namespace std;

typedef struct {
	double	maxcos;
	Color	color;
} TFaceData;

/* --- faceType --- */
static int faceType(Face* /*face*/, double /*maxcos*/)
{
#if 0
	for (int i=0; i<3; i++) {
		Edge *e = face->edge(i);
		if (e->faceCount()==1)
			return FACE_ORPHAN;
		else if (e->faceCount()==2) {
			Face *nf = e->face(0);
			if (nf==face) nf = e->face(1);
			if (face->normal().dot(nf->normal()) < maxcos)
				//cout << face->normal().dot(nf->normal()) << endl;
				return FACE_FOLD;
		}
	}
#endif
	return FACE_NORMAL;
} // faceType

/* --- drawFace --- */
static int drawFace(Face *face, void *data)
{
	switch (faceType(face, ((TFaceData*)data)->maxcos)) {
		case FACE_NORMAL: {
				glColor3f(0.8, 0.8, 0.8);
				//Color color  = ((TFaceData*)data)->color;
				//GLCOLOR(color);
			}
			break;

		case FACE_FOLD:
			glColor3f(0., 0., 1.);
			break;

		case FACE_ORPHAN:
			glColor3f(1., 0., 0.);
			break;
	}

	glNormal3d( face->normal().x, face->normal().y, face->normal().z );
//	glVertex3d( face->A().x, face->A().y, face->A().z );

	glVertex3d( face->A().x, face->A().y, face->A().z );
	glVertex3d( face->B().x, face->B().y, face->B().z );
	glVertex3d( face->C().x, face->C().y, face->C().z );
	return false;
} // drawFace

/* --- drawNormal --- */
static int drawNormal(Face* face, void* data)
{
	double nl = *(double*)data;
//	glNormal3dv( face->normal().get() );
//	glNormal3d( face->normal().x, face->normal().y, face->normal().z );
	Vector center = face->center();
	glVertex3d( center.x, center.y, center.z );
	Vector end = center + face->normal() * nl;
	glVertex3d( end.x, end.y, end.z );
	return false;
} // drawNormal

/* --- drawVertex --- */
#if 0
static void drawVertex(Vertex *vertex, void *)
{
//	if (bSelected)
//		glColor3f(1.,.5,.5);
//	else
		glColor3f(0.8, 1.,0.8);
	glNormal3dv( face->normal()->get() );
	glVertex3dv( vertex->get() );
} // drawVertex
#endif

/* --- drawEdge --- */
static int drawEdge(Edge *edge, void *)
{
	//glNormal3dv( face->normal().get() );

	if (edge->show)
		glColor3f(1.0, 0.0, 0.0);
	else
		glColor3f(0.0, 0.5, 0.0);

	glVertex3d( edge->a->x, edge->a->y, edge->a->z );
	glVertex3d( edge->b->x, edge->b->y, edge->b->z );

	return false;
} // drawEdge

#if 0
/* --- drawVertexName --- */
static void drawVertexName(Vertex *v, void *data)
{
/*
	 int *name=(int*)data[0];
	 v->name=name[0];
	 glLoadName(name[0]);
	 name[0]++;
	 glBegin(GL_POINTS);
	 glVertex3d(p->x,p->y,p->z);
	 glEnd();
*/
} // drawVertexName
#endif

/* --- draw --- */
void GLMesh::draw(GLenum mode)
{
	GLfloat	polyfactor = 1.0;
	GLfloat polyunits  = 1.0;

	if (mode==GL_RENDER) {
		glEnable(GL_LIGHTING);
	} else
	if (mode==GL_SELECT) {
	}

//	if (fill)
//		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//	else
//		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);


	if (drawFaces) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(polyfactor,polyunits);
		glBegin(GL_TRIANGLES);
		TFaceData fd;
		fd.maxcos = foldCos;
//		fd.color  = color();
		mesh->forEachFace(drawFace, &fd);
		glEnd();
		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	glDisable(GL_LIGHTING);

	if (drawNormals) {
//		GLCOLOR(normalColor);
		glBegin(GL_LINES);
		mesh->forEachFace(drawNormal, &normalLength);
		glEnd();
	}

//	if (drawVertices) {
//		glPointSize(6.);
//		glColor3f(1., 0., 1.);
//		glBegin(GL_POINTS);
//		mesh->forEachVertex(drawVertex);
//		glEnd();
//		glPointSize(1.);
//	}

	if (drawEdges) {
		glBegin(GL_LINES);
		mesh->forEachEdge(drawEdge);
		glEnd();
	}
} // draw
