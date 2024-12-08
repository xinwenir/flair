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
GVoxel	voxel;

int	freq[0xFFFF];

/** replace */
void replace(word rtstruct, word hu)
{
	memset(freq, 0, sizeof(freq));
	// build a frequency table
	for (size_t ptr=0; ptr<voxel.size(); ptr++) {
		word w = voxel.data(ptr);
		freq[w]++;
	}

	word region = voxel.addKreg(hu);
	hu = voxel.mo;

	cout << endl << "Frequencies (Before)" << endl;
	for (word l=0; l<=voxel.mo; l++)
		cout << l << " " << freq[l] << endl;

	for (int k=0; k<voxel.nz; k++)
		for (int j=0; j<voxel.ny; j++)
			for (int i=0; i<voxel.nx; i++) {
				word roi = voxel.roi(i,j,k);
				if (roi == rtstruct) {
					word w = voxel.data(i,j,k);
					freq[w]--;
					voxel.data(i,j,k,hu);
				}
			}

	cout << endl << "Frequencies (After)" << endl;
	for (word l=0; l<=voxel.mo; l++)
		cout << l << " " << freq[l] << endl;

	int ii=0, jj=0, kk=0;
	int count=0;
	for (word l=0; l<=voxel.mo; l++)
		if (freq[l]==0) {
			voxel.data(ii,jj,kk, l);
			count++;
			if (++ii == voxel.nx) {
				ii = 0;
				if (++jj == voxel.ny) {
					jj = 0;
					kk++;
				}
			}
		}

	memset(freq, 0, sizeof(freq));
	// build a frequency table
	for (int k=0; k<voxel.nz; k++)
		for (int j=0; j<voxel.ny; j++)
			for (int i=0; i<voxel.nx; i++) {
				word w = voxel.data(i,j,k);
				freq[w]++;
			}

	cout << endl << "Frequencies (AfterCorrection)" << endl;
	for (word l=0; l<=voxel.mo; l++)
		cout << l << " " << freq[l] << endl;

	cout << "Corrected: " << count << endl;
	cout << "Added hu:" << voxel.mo << " region:" << region << endl;
} // replace

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
	string	input_file, output_file;
	word	rtstruct = 0xFFFF;
	word	hu       = 0xFFFF;

	while (true) {
		int option_index = 0;
		static struct option long_options[] = {
			{"help",    0, 0, '?'},
			{"input",   1, 0, 'i'},
			{"output",  1, 0, 'o'},
			{"struct",  1, 0, 's'},
			{"hu",      1, 0, 'h'},
			{0, 0, 0, 0}
		};

		int c = getopt_long  (ac, av, "?i:o:s:h:",
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
				usage();
				break;

			case 'i':
				input_file = optarg;
				break;

			case 'o':
				output_file = optarg;
				break;

			case 's':
				rtstruct = atoi(optarg);
				break;

			case 'h':
				hu = atoi(optarg);
				break;

			default:
				cerr << "?? getopt returned character code " << c << endl;
				usage();
		}
	}
	if (rtstruct==0xFFFF || hu==0xFFFF) {
		cerr << "ERROR: didn't supply structure and hu to replace" << endl;
		return 1;
	}

	if (!input_file.length()) {
		cerr << "ERROR: input file not provided" << endl;
		return 1;
	}

	if (!output_file.length()) {
		cerr << "ERROR: output file not provided" << endl;
		return 1;
	}

	if (!voxel.load(input_file.c_str(), true)) {
		cerr << "ERROR: loading voxel file " << input_file << endl;
		return 2;
	}

	if (voxel.nroi()==0) {
		cerr << "ERROR: voxel do not contain ROI STRUCT information" << endl;
		return 3;
	}

	replace(rtstruct, hu);

	voxel.save(output_file.c_str());
} // main
