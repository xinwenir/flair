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
	int		id;
	int		volume;
	double		maxdose;
	Histogram	histogram;
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

void background()
{
	//OK it can be done better
	int p =0;
	int c =0;
	for (int k=0; k<usrbin.nz; k++) {
		double kk = voxel.voxelk(usrbin.usrbincz(k));
		for (int j=0; j<usrbin.ny; j++) {
			double jj = voxel.voxelj(usrbin.usrbincy(j));
			for (int i=0; i<usrbin.nx; i++) {
				double ii = voxel.voxeli(usrbin.usrbincx(i));
				const ROICombination& comb = voxel.roiComb(ii,jj,kk);
				c++;
				if (comb.length==0) {
					usrbin.set(i,j,k,0);
					p++;
				}
		        }
		}
	}
	cout << "I have changed "<<p<<" bins out of "<<c <<"\n";
}

/* --- usage --- */
void usage()
{
	fprintf(stderr,"syntax: %s [options] <infile>\n",prgname);
	fprintf(stderr,"options:\n");
	fprintf(stderr,"\t-h?, --help\t\tThis help page\n");
	fprintf(stderr,"\t-d, --detector #\tDetector index number starting from 1\n");
	fprintf(stderr,"\t-n, --norm #\tNormalization factor\n");
	fprintf(stderr,"\t-o, --output file\tA tab_lis file with the dvh histograms\n");
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
			{"usrbin",  1, 0, 'u'},
			{"voxel",   1, 0, 'v'},
			{"xvoxel",  1, 0, 'x'},
			{"yvoxel",  1, 0, 'y'},
			{"zvoxel",  1, 0, 'z'},
			{0, 0, 0, 0}
		};

		int c = getopt_long  (ac, av, "h?d:n:o:u:v:x:y:z:",
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

	background();
        if (fileout) {
		fout = fopen(fileout,"w");
		if (fout==NULL) {
			cerr << "ERROR: opening output file " << optarg << endl;
			return 4;
		}
	}
	if (usrbin.save(fileout)) cout << "File saved " << fileout << endl;

	if (fout != stdout) fclose(fout);
	if (fileout) free(fileout);


} // main
