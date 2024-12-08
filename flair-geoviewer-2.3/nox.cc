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
 * Date:        19 Feb 2013
 */

#include <fstream>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "os.h"
#include "nox.h"

#define _DRAW_2D true
#define _DRAW_3D false
#define _DRAW_BBOX false

using namespace std;

/** NoxGeometry */
NoxGeometry::NoxGeometry(const char *filename)
	: kernel(geometry),
	  viewer(geometry,kernel),
	  dump(false)
{
	geometry.backgroundColor(0x707070);
	if (filename) load(filename);
} // NoxGeometry

/** parse matrix */
void NoxGeometry::parseMatrix(Matrix4& matrix)
{
	matrix.identity();
	for (int row=0; row<4; row++)
		for (int col=0; col<4; col++) {
			const char *token = strtok(NULL," ");
			if (!token) break;
			matrix(row,col) = atof(token);
		}
} // parseMatrix

/** parseBody */
void NoxGeometry::parseBody(const Matrix4* matrix)
{
	GBody* body;
	char *type = strtok(NULL," ");
	char *name = strtok(NULL," ");

	double what[30];
	memset(what, 0, sizeof(what));

	if (!strcmp(type,"VOXELS")) {
		if (!geometry.voxel.load(name)) {
			cerr << "ERROR Unable to load voxel " << name << endl;
		}
		geometry.voxel.xlow = atof(strtok(NULL," "));
		geometry.voxel.ylow = atof(strtok(NULL," "));
		geometry.voxel.zlow = atof(strtok(NULL," "));
		geometry.voxel.calcLimits();

		// create a voxel body
		body = geometry.addBody("VOXEL", "RPP");
		what[0] = geometry.voxel.xlow;
		what[1] = geometry.voxel.xhigh;
		what[2] = geometry.voxel.ylow;
		what[3] = geometry.voxel.yhigh;
		what[4] = geometry.voxel.zlow;
		what[5] = geometry.voxel.zhigh;

	} else {
		body = geometry.addBody(name, type);

		int i=0;
		const char *token;
		while ((token = strtok(NULL," ")))
			what[i++] = atof(token);
	}

	char err[128];
	err[0] = 0;
	body->setWhat(what,err);
	if (err[0]) cerr << err << endl;
	body->create();
	if (matrix) {
		body->matrix(*matrix);
		body->transform();
	}
	if (dump) cout << "Body: " << name << " " << type << endl;
} // parseBody

/** parseRegion */
void NoxGeometry::parseRegion()
{
	char *name = strtok(NULL," ");
	char err[128];
	GRegion *gregion = geometry.addRegion(name);
	if (dump) cout << "Region: " << name << endl;
	gregion->type((RegionType)atoi(strtok(NULL," ")));
	//	gregion->rotdefi = atoi(strtok(NULL," "));

	region_material_t regmat;
	regmat.region_id = gregion->id();
	regmat.color = atoi(strtok(NULL," "));
	regmat.alpha = atoi(strtok(NULL," "));
	region_materials.add(regmat);
#if _DRAW_BBOX
	gregion->show = BIT_BBOX;
#endif
	const char *token;

	while ((token=strtok(NULL," "))) {
		geometry.addZone(gregion,(bool)!strcmp(token,"RPN"));
		while ((token=strtok(NULL," "))) {
			if (token[0]==':')
				break;
			else
			if (geometry.add2exp(gregion, token, err))
				cerr << err << endl;
		}
		if (token==NULL) break;
	}
} // parseRegion

/** parseVoxel */
void NoxGeometry::parseVoxel()
{
	char *cmd = strtok(NULL," ");

	if (!strcmp(cmd,"color")) {
		int   i     = atoi(strtok(NULL," "));
		dword color = (dword)atoi(strtok(NULL," "));

		viewer.voxel().color(i-1, color);
	}
} // parseVoxel

/** parse */
void NoxGeometry::parse(FILE *input)
{
	char str[32000];
	bool hasMatrix = false;
	Matrix4 matrix;

	while (fgets(str,sizeof(str),input)) {
		str[strlen(str)-1] = 0;
		if (str[0] == 0 || str[0] == '#') continue;

		const char *token = strtok(str," ");

		//		if (!strcmp(token,"rotdefi")) {
		//			int idx = atoi(strtok(NULL," "));
		//			parseMatrix(matrix);
		//			geometry.rotdefi(idx,matrix);
		//		} else
		if (!strcmp(token,"matrix")) {
			parseMatrix(matrix);
			hasMatrix = true;
			if (dump) cout << "\tMatrix" << endl;
		} else
		if (!strcmp(token,"body")) {
			parseBody(hasMatrix? &matrix : NULL);
			hasMatrix = false;
		} else
		if (!strcmp(token,"region"))
			parseRegion();
		else
		if (!strcmp(token,"voxel"))
			parseVoxel();
		else
		if (!strcmp(token,"viewmatrix")) {
			parseMatrix(matrix);
			if (dump) cout << "View Matrix:" << endl << matrix << endl;
			viewer.matrix(matrix);
		} else
		if (!strcmp(token,"extends")) {
			double x1 = atof(strtok(NULL," "));
			double y1 = atof(strtok(NULL," "));
			double x2 = atof(strtok(NULL," "));
			double y2 = atof(strtok(NULL," "));
			if (dump) cout	<< "Extends [" << x1 << ", " << y1
					<< "] - [" << x2 << ", " << y2 << "]" << endl;
			viewer.window(x1,y1,x2,y2);
		} else
		if (!strcmp(token,"origin")) {
			double x = atof(strtok(NULL," "));
			double y = atof(strtok(NULL," "));
			double z = atof(strtok(NULL," "));
			if (dump) cout	<< "Origin [" << x << ", " << y
					<< ", " << z << "]" << endl;
			viewer.origin(x,y,z);
		} else
		if (!strcmp(token,"zoom")) {
			double z = atof(strtok(NULL," "));
			if (dump) cout	<< "Zoom " << z << endl;
			viewer.view().zoom(z);
		}
	}
	kernel.derive();
	assignRegionMaterials();	// must be called after derive
} // parse

/** assignRegionMaterials */
void NoxGeometry::assignRegionMaterials()
{
	//assign region materials
	for(int i=0; i<region_materials.size(); i++) {
		const region_material_t& regmat = region_materials[i];
		VRegion* region = kernel.getRegion(regmat.region_id);
		if (region) {
			region->color(regmat.color);
			region->alpha(regmat.alpha);
		}
	}
} // assignRegionMaterials

/** load */
bool NoxGeometry::load(const char *filename)
{
	FILE *input = fopen(filename, "r");
	if (input == NULL) {
		cerr << "Error opening file " << filename << endl;
		return false;
	}
	parse(input);
	fclose(input);
	return true;
} // load

/** cleanup */
void NoxGeometry::cleanup()
{
	geometry.cleanup();
} // cleanup
