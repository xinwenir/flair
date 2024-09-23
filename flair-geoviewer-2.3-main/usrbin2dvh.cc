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
 * Date:	02-Feb-2015
 */

#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

#include <iostream>

#include "os.h"
#include "voxel.h"
#include "usrbin.h"
#include "histogram.h"

#ifdef MEM
#include "memory.h"
#endif

using namespace std;

/** global variables */
char	*prgname;
Usrbin	usrbin;
GVoxel	voxel;
bool	same;		// if usrbin has same dimensions as the voxel
double	relative = 0.0;
double	norm     = 1.0;

class DVH {
public:
	int	id;
	int	volume;
	double	maxdose;
	H1D	histogram;
public:
	DVH()	: id(0), volume(0.), maxdose(0.) {}

	/* compare maxdose */
	void	cmpDose(double d)	{ maxdose = max(maxdose, d); }
}; // DVH histogram information

DVH	*dvh = NULL;


//int	*volumes=NULL;	// array with ROI structure volume information (should move to GVoxel)
//double	*maxdose=NULL;	// array with maximum dose information

/** check if usrbin and voxel have the same binning */
bool checkUsrbinNVoxel()
{
	const double eps = 1e-7;

	if (usrbin.nx!=voxel.nx) return false;
	if (usrbin.ny!=voxel.ny) return false;
	if (usrbin.nz!=voxel.nz) return false;

	if (Eq(usrbin.xlow, voxel.xlow, eps)) return false;
	if (Eq(usrbin.xlow, voxel.xlow, eps)) return false;
	if (Eq(usrbin.xlow, voxel.xlow, eps)) return false;

	if (Eq(usrbin.xhigh, voxel.xhigh, eps)) return false;
	if (Eq(usrbin.xhigh, voxel.xhigh, eps)) return false;
	if (Eq(usrbin.xhigh, voxel.xhigh, eps)) return false;

	return true;
} // checkUsrbinNVoxel

/** scan4Structures: scan voxel for structures return volumes and maximum dose per structure */
void scan4Structures()
{
	/** crap */
	dvh = new DVH[voxel.nroi()+1];

	for (int k=0; k<voxel.nz; k++) {
		double z = voxel.voxelcz(k);
		for (int j=0; j<voxel.ny; j++) {
			double y = voxel.voxelcy(j);
			for (int i=0; i<voxel.nx; i++) {
				double x = voxel.voxelcx(i);
				const ROICombination& comb = voxel.roiComb(i,j,k);
				double dose;
				bool ok;
				if (same)
					dose = norm * usrbin.get(i,j,k);
				else
					dose = norm * usrbin.getData(x,y,z,&ok);

				for (int m=0; m<comb.length; m++) {
					int roi = Abs(comb[m]);
					dvh[roi].volume++;
					dvh[roi].cmpDose(dose);
				}
			}
		}
	}

	/* prepare histograms */
	dvh[0].histogram.title("void");
	dvh[0].histogram.set(100, 0.0, dvh[0].maxdose);
	for (int i=1; i<=voxel.nroi(); i++) {
		dvh[i].histogram.title(voxel.roiName(i));
		dvh[i].histogram.set(100, 0.0, dvh[i].maxdose);
	}

	/* print maximum doses and values */
	printf("ROI\tVolume\tMax Dose\n");
	double dxdydz = voxel.dx * voxel.dy * voxel.dz;
	for (int i=0; i<=voxel.nroi(); i++) {
		string name;
		if (i==0)
			name = "void";
		else
			name = voxel.roiName(i);
		if (name.length()==0) continue;
		printf("%2d\t%-24s\t%.5g\t%.5g\n", i, name.c_str(),
		       dxdydz*(double)dvh[i].volume, dvh[i].maxdose*norm); 
	}
} // scan4Structures

/** fillDVH */
void fillDVH()
{
	// prepare histograms
	dvh[0].histogram.title("void");
	dvh[0].histogram.set(100, 0.0, dvh[0].maxdose);
	for (int i=1; i<=voxel.nroi(); i++) {
		dvh[i].histogram.title(voxel.roiName(i));
		if (dvh[i].maxdose > 0.0)
			dvh[i].histogram.set(100, 0.0, dvh[i].maxdose);
	}

	// Fill the histograms
	for (int k=0; k<voxel.nz; k++) {
		double z = voxel.voxelcz(k);
		for (int j=0; j<voxel.ny; j++) {
			double y = voxel.voxelcy(j);
			for (int i=0; i<voxel.nx; i++) {
				double x = voxel.voxelcx(i);
				const ROICombination& comb = voxel.roiComb(i,j,k);
				double dose;
				bool ok;
				if (same)
					dose = norm * usrbin.get(i,j,k);
				else
					dose = norm * usrbin.getData(x,y,z,&ok);

				if (dose>0.0)
					for (int m=0; m<comb.length; m++)
						dvh[Abs(comb[m])].histogram.fill(dose);
			}
		}
	}

	// convert to cumulative
	for (int i=0; i<=voxel.nroi(); i++)
		if (dvh[i].maxdose > 0.0) {
			if (relative!=0.0)
				//dvh[i].histogram.setLimits(0., 100.0);
			  dvh[i].histogram.setLimits(0., dvh[i].maxdose/relative*100.0);
			dvh[i].histogram.cumulative(true);
			dvh[i].histogram.norm(100.0/dvh[i].histogram.total());
		}
} // fillDVH

/** writeDVH */
void writeDVH(FILE *fout)
{
	/* First write volumes and maxdose */
	/* print maximum doses and values */
	double dxdydz = voxel.dx * voxel.dy * voxel.dz;

	int det = 1;
	fprintf(fout,"# Detector: %d volumes, maxdose\n",det++);
	for (int i=0; i<voxel.nroi()+1; i++)
		fprintf(fout,"%2d \"%s\" %.5g %.5g\n", i, dvh[i].histogram.title(), dxdydz*(double)dvh[i].volume, dvh[i].maxdose);
	fprintf(fout,"\n\n");

	for (int i=0; i<voxel.nroi()+1; i++) {
		if (dvh[i].histogram.title()[0]==0 || dvh[i].maxdose==0.0) continue;
		fprintf(fout,"# Detector: %d %s\n",det++, dvh[i].histogram.title());
		dvh[i].histogram.save(fout,4);
		if (i<voxel.nroi()) fprintf(fout,"\n\n");
	}
} // writeDVH

/* --- usage --- */
void usage()
{
	fprintf(stderr,"syntax: %s [options] <infile>\n",prgname);
	fprintf(stderr,"options:\n");
	fprintf(stderr,"\t-h?, --help\t\tThis help page\n");
	fprintf(stderr,"\t-d, --detector #\tDetector index number starting from 1\n");
	fprintf(stderr,"\t-n, --norm #\tNormalization factor\n");
	fprintf(stderr,"\t-o, --output file\tA tab_lis file with the dvh histograms\n");
	fprintf(stderr,"\t-r, --relative #\tRelative to target dose (Default absolute)\n");
	fprintf(stderr,"\t-u, --usrbin #\t\tUsrbin file to load\n");
	fprintf(stderr,"\t-v, --voxel #\t\tVoxel file containing the ROI STRUCTures\n");
	fprintf(stderr,"\t-x, --xvoxel #\tX voxel origin (cm)\n");
	fprintf(stderr,"\t-y, --yvoxel #\tY voxel origin (cm)\n");
	fprintf(stderr,"\t-z, --zvoxel #\tZ voxel origin (cm)\n");
	fprintf(stderr,"author: Vasilis.Vlachoudis@cern.ch (2015)\n");
	exit(0);
} // usage

/** main */
int main(int ac, char *av[])
{
	prgname = av[0];
	string	usrbin_file, voxel_file;
	int	detector=0;
	char*	fileout=NULL;
	FILE*	fout=stdout;

	while (true) {
		int option_index = 0;
		static struct option long_options[] = {
			{"help",    0, 0, 'h'},
			{"detector",1, 0, 'd'},
			{"norm",    1, 0, 'n'},
			{"output",  1, 0, 'o'},
			{"relative",1, 0, 'r'},
			{"usrbin",  1, 0, 'u'},
			{"voxel",   1, 0, 'v'},
			{"xvoxel",  1, 0, 'x'},
			{"yvoxel",  1, 0, 'y'},
			{"zvoxel",  1, 0, 'z'},
			{0, 0, 0, 0}
		};

		int c = getopt_long  (ac, av, "h?d:n:o:r:u:v:x:y:z:",
			long_options, &option_index);
		if (c == -1)
			break;

		switch (c) {
			case 0:
				printf ("option %s", long_options[option_index].name);
				if (optarg)
					printf (" with arg %s", optarg);
				printf ("\n");
				usage();
				break;

			case '?':
			case 'h':
				usage();
				break;

			case 'd':
				detector = atoi(optarg);
				break;

			case 'n':
				norm = atof(optarg);
				break;

			case 'o':
				fileout = strdup(optarg);
				break;

			case 'r':
				relative = atof(optarg);
				break;

			case 'u':
				usrbin_file = optarg;
				break;

			case 'v':
				voxel_file = optarg;
				break;

			case 'x':
				voxel.xlow = atof(optarg);
				break;

			case 'y':
				voxel.ylow = atof(optarg);
				break;

			case 'z':
				voxel.zlow = atof(optarg);
				break;

			default:
				cerr << "?? getopt returned character code " << c << endl;
				usage();
		}
	}
	voxel.calcLimits();
//	if (optind < ac) {
//		if ((fin=fopen(av[optind],"rb"))==NULL) {
//			fprintf(stderr,"Error opening input file %s\n",av[2]);
//			return 2;
//		}
//	}
	if (usrbin_file.length()==0) {
		cerr << "ERROR: usrbin file was not provided\n" << endl;
		return 1;
	}
	if (detector <= 0) {
		cerr << "ERROR: usrbin detector index was not provided\n" << endl;
		return 2;
	}
	if (voxel_file.length()==0) {
		cerr << "ERROR: voxel file was not provided\n" << endl;
		return 3;
	}

	if (!usrbin.load(usrbin_file.c_str(), detector)) {
		cerr << "ERROR: loading usrbin detector " << usrbin_file << endl;
		return 4;
	}
	if (!voxel.load(voxel_file.c_str())) {
		cerr << "ERROR: loading voxel file " << voxel_file << endl;
		return 5;
	}
	if (voxel.nroi()==0) {
		cerr << "ERROR: voxel do not contain ROI STRUCT information" << endl;
		return 6;
	}

	same = checkUsrbinNVoxel();
	scan4Structures();
	fillDVH();

	if (fileout) {
		fout = fopen(fileout,"w");
		if (fout==NULL) {
			cerr << "ERROR: opening output file " << optarg << endl;
			return 4;
		}
	}
	writeDVH(fout);

	if (dvh) delete [] dvh;

	if (fout != stdout) fclose(fout);
	if (fileout) free(fileout);
} // main
