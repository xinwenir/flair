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
 */

#include <string.h>

#include "viewer.h"
#include "geometry.h"
#include "exportlayer.h"

#include "dxfexport.h"
#include "svgexport.h"
#include "exportbase.h"

using namespace std;

/** draw */
void ExportLayer::draw(const dword)
{
} // draw

/** exportSegment
 * @param body	body of which to draw the segment
 * @param ic	index of conic to draw
 * @param type	segment type
 * @param ts,xs,ys	start coordinates of segment
 * @param te,xe,ye	end coordinates of segment
 * @param step	I/O step used in drawing (<0.0 to reset)
 */
void ExportLayer::exportSegment(ExportBase& exporter, VBody *body,
			const int ic, const dword color,
			const double ts, const double xs, const double ys,
			const double te, const double xe, const double ye,
			double *step)
{
	double x1,y1, x2,y2;
	const Conic& conic = body->C[ic];

	if (ISZERO(Sqr(xe-xs) + Sqr(ye-ys)) && ISEQ(ts,te))
		return;

	// draw lines
	x1 = xs;
	y1 = ys;
	if (conic.type() == CONIC_LINE) {
		x2 = xe;
		y2 = ye;
		exporter.line(x1,y1,x2,y2,color,body->name());
	} else {
		Path path;
		path.add(Vector2D(x1,y1));
		// initial guess
		if (*step<=0.0) *step = MAXSTEP;
		if (*step > te-ts) *step = te-ts;
		double tp = ts;
		do {
			double t = tp + *step;
			if (t>te) {
				t = te;
				*step = te - ts;
			}
			double x, y;
			conic.getXY(t, &x, &y);
			x2 = x;
			y2 = y;
			if (Abs(x2-x1)<=2 && Abs(y2-y1)<=2) {
				path.add(Vector2D(x2,y2));
				x1 = x2;
				y1 = y2;
				tp = t;
				*step *= 2.0;
				continue;
			}

			// find at half way what is the distance from a line
			while (1) {
				double xh, yh;
				double th = tp + *step/2.0;
				assert(th<=te);
				conic.getXY(th, &xh, &yh);
				double x3 = xh;
				double y3 = yh;
				if (Abs(x3-x1)<=2 && Abs(y3-y1)<=2) {
					path.add(Vector2D(x3,y3));
					x1 = x3;
					y1 = y3;
					tp = th;
					break;
				}
				// distance of line 1,2
				double l2 = Sqr((double)(x2-x1)) + Sqr((double)(y2-y1));
				// perpendicular distance of midpoint 3 wrt to line 1,2
				double d2 = Sqr((double)(x3-x1)*(double)(y2-y1)
					      - (double)(y3-y1)*(double)(x2-x1))
					      / (double)l2;
				if (d2<=0.01 && l2<50.0) {
					path.add(Vector2D(x2,y2));
					x1 = x2;
					y1 = y2;
					tp = t;
					*step *= 2.0;
					if (*step > MAXSTEP) *step = MAXSTEP;
					break;
				} else
				if (d2<=0.1) {
					path.add(Vector2D(x2,y2));
					x1 = x2;
					y1 = y2;
					tp = t;
					break;
				} else {
					if (*step <= MINSTEP+SMALL) {
						path.add(Vector2D(x2,y2));
						x1 = x2;
						y1 = y2;
						tp = t;
						break;
					}
					x2 = x3;	// half
					y2 = y3;
					t  = th;
					*step /= 2.0;
					if (*step < MINSTEP) {
						*step = MINSTEP;
						break;
					}
				}
			}
		} while (tp < te - *step);
		x2 = xe;
		y2 = ye;
		path.add(Vector2D(x2,y2));
		exporter.polyline(path, color, body->name());
	}
} // exportSegment

/** export geometry to file */
void ExportLayer::exportSegments(ExportBase& exporter, VBody *body)
{
	dword	type;
	double	tp, xp, yp;
	double	t, x, y;
	double	step;
	const dword mask = SEGMENT_REGION | SEGMENT_ERROR;

	for (int i=0; i < body->nC; i++) {
		if (body->V[i].size()<2) continue;

		ArrayIterator<Vertex2D> viter(body->V[i]);

		// Initial vertex
		Vertex2D* v = &(viter++);
		tp = v->t;
		xp = v->x;
		yp = v->y;

		// Next vertex
		v = &(viter++);
		type = v->type&mask;
		t = v->t;
		x = v->x;
		y = v->y;

		step = -1.0;

		while (viter) {
			// Check next vertex
			v = &(viter++);
			if ((v->type&mask) == type) {
				t = v->t;
				x = v->x;
				y = v->y;
				continue;
			}

			if (type&mask)
				exportSegment( exporter, body, i,
					(type==SEGMENT_REGION? geometry.regionColor : geometry.errorColor),
					tp, xp, yp, t, x, y, &step);

			tp = t;
			xp = x;
			yp = y;
			type = v->type;
			t = v->t;
			x = v->x;
			y = v->y;
		}
		if (type&mask)
			exportSegment( exporter, body, i,
				(type==SEGMENT_REGION? geometry.regionColor : geometry.errorColor),
				tp, xp, yp, t, x, y, &step);
	}
} // exportSegments

/** exportDXF
 * @param mask	of options to draw
 * @return	projection state that was drawn
 */
void ExportLayer::exportDXF(const char *filename)
{
	// create exporter
	DXFExport exporter(filename);
	if( !exporter ) {
		kernel.error("Geometry could not be exported to DXF.\n");
#if _DEBUG>1
		cerr << "Geometry could not be exported to DXF." << endl;
		cerr << "Output stream " << filename << " could not be openend!" << endl;
#endif
		return;
	}

	if (viewer.state()==PROJECTION_FINISHED)
		for (int i=0; i<geometry.bodies.count(); i++)
			exportSegments(exporter, kernel.getBody(i));
} // exportDXF

/** exportSVG
 * @param mask	of options to draw
 * @return	projection state that was drawn
 */
void ExportLayer::exportSVG(const char *filename)
{
	// create exporter
	SVGExport exporter(filename);
	if( !exporter ) {
		kernel.error("Geometry could not be exported to SVG.\n");
#if _DEBUG>1
		cerr << "Geometry could not be exported to SVG." << endl;
		cerr << "Output stream " << filename << " could not be opened!" << endl;
#endif
		return;
	}

	if (viewer.state()==PROJECTION_FINISHED)
		for (int i=0; i<geometry.bodies.count(); i++)
			exportSegments(exporter, kernel.getBody(i));
} // exportSVG

/** export image to filename */
void ExportLayer::operator()(const char *filename)
{
	// a bit stupid search!
	if (strstr(filename,".dxf") || strstr(filename,".DXF"))
		exportDXF(filename);
	else
	if (strstr(filename,".svg") || strstr(filename,".SVG"))
		exportSVG(filename);
} /* operator() */
