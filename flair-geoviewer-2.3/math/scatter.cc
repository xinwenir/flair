/**
 * $Id$
 *
 * Author:	Vasilis.Vlachoudis@cern.ch
 * Date:	07-Mar-2016
 */

#include <stdio.h>

#include "os.h"
#include "scatter.h"

/** Scatter */
Scatter::~Scatter()
{
//	for (int i=0; i<data.size(); i++)
//		delete data[i];
} // Scatter

/** load */
bool Scatter::load(const char* filename)
{
	FILE* fin = fopen(filename,"r");
	if (fin==NULL) return false;

	while (true) {
		double x,y;
		if (!fscanf(fin,"%lg %lg",&x, &y)) break;
		if (feof(fin)) break;
		data.add(XYPair(x,y));
	}
	fclose(fin);
//	for (int i=0; i<data.count(); i++)
//		printf("%d %g %g\n",i,data[i].x,data[i].y);
	return true;
} // load

/** @return linear interpolated value of y */
double Scatter::interpolate(const double x) const
{
	int sp = data.search(x);
	if (sp<0 || sp>=data.count()) return 0.0;
	int si = sp+1;

	/* interpolate value */
	double xlow = data[sp].x;
	double ylow = data[sp].y;
	return ylow + (x-xlow) * (data[si].y-ylow) / (data[si].x-xlow);
} // interpolate

/** @return log interpolated value of y */
double Scatter::interpolateLog(const double x) const
{
	int sp = data.search(x);
	if (sp<0 || sp>=data.count()) return 0.0;
	int si = sp+1;

	/* interpolate value */
	double xlow  = log(data[sp].x);
	double ylow  = log(data[sp].y);
	double xhigh = log(data[si].x);
	double yhigh = log(data[si].y);

	return exp(ylow + (log(x)-xlow) * (yhigh-ylow) / (xhigh-xlow));
} // interpolate
