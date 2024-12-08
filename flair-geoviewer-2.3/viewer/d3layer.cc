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
 *
 * TODO
 * - draw with hash lines internal of bodies (cut by clipbody)
 */

#include <iomanip>

#include "bbox.h"
#include "array.h"
#include "obbox.h"
#include "viewer.h"
#include "random.h"
#include "d3layer.h"
#include "geometry.h"

#ifdef THREAD

/* =========================== Body3DWorker =========================== */
/** operator() */
void Body3DWorker::operator()()
{
	int sample = 1;
	for (int sx=0; sx<d3->_samples; sx++, sample++) {
		double u = d3->view().i2u((double)ii + ((double)sx+0.5)/(double)d3->_samples);
		for (int sy=0; sy<d3->_samples; sy++) {
			double v = d3->view().j2v((double)jj + ((double)sy+0.5)/(double)d3->_samples);
			zone = d3->drawPixel(engine, ray, ptr, W,H, ii,jj,step,sample, u,v, zone);
		}
	}
} // operator

/* =========================== Body3DFeeder =========================== */
/** allocate */
void Body3DFeeder::allocate()
{
	if (nworkers != pool.nthreads()) {
		if (workers) delete [] workers;
		workers  = new Body3DWorker[pool.nthreads()];
		nworkers = pool.nthreads();
	}
} // allocate

/** reset */
void Body3DFeeder::reset(Painter* p, Ray& ray)
{
	pool.reset();
	allocate();

	// painter width&height
	painter = p;
	W = painter->width();
	H = painter->height();

	// initialize workers
	for (int k=0; k<nworkers; k++) {
		workers[k].d3     = d3;
		workers[k].engine = d3->kernel.engine(k);
		workers[k].ray    = ray;
		workers[k].zone   = NULL;
		workers[k].W      = W;
		workers[k].H      = H;
	}
} // reset

/** init */
void Body3DFeeder::init(int astep)
{
	// nested loops initialization
	step = astep;
	ii   = 0;
	jj   = 0;
	ptr = pline = painter->data();
} // init

/** loop - for next pixel
 * @return true if we have finished, false otherwise
 */
bool Body3DFeeder::loop()
{
	// (ii) loop
	ii  += step;
	ptr += step;
	if (ii>=W) {
		// (jj) loop
		jj     += step;
		if (jj>=H) return true;

		pline += step*W;

		ii = 0;
		ptr = pline;
	}
	return false;
} // loop

/** feed */
ThreadPoolWorker* Body3DFeeder::feed(int threadId)
{
	// Check for ending condition
	if (d3->viewer.stop()) return NULL;

	// update worker
	workers[threadId].ii   = ii;
	workers[threadId].jj   = jj;
	workers[threadId].ptr  = ptr;
	workers[threadId].step = step;

	if (loop()) return NULL;
	return (ThreadPoolWorker*)&workers[threadId];
} // feed
#endif

/* ============================== D3Layer ============================= */
D3Layer::D3Layer(const Geometry& g, GeometryKernel& k, GeometryViewer& v) :
		Layer(g, k, v)
#ifdef THREAD
		, feeder(k.threadpool, this)
#endif
{
	deflights     = false;
	ambient       = 0x30;		// ambient intensity
	xray          = 0;
	antialias     = 1;
	drawEdges     = false;
	shadows       = true;
	skip_1stblack = false;
	_step         = 32;		// initial guess
	_maxDepth     = 0;
	_usrbinastexture = true;
	drawTime(15);
} // D3Layer

/** nextIntersection
 * @return false if the ray is ended or on error, true otherwise
 */
bool D3Layer::nextIntersection(GeometryEngine* eng, Ray *ray)
{
	Point hit;
INT:	if (eng->intersectRay(ray, false))
		// End of ray, maybe we should return sky color?
		return false;

	RaySegment &segment = ray->segment();
	VZone *zone = segment.zone;

	if (zone == NULL) {
		ray->error = true;
		return false;
	}

	// Set default normal the opposite direction of the beam
	ray->normal = -ray->segment().dir;

	int n = ray->n;

	if (ray->project_hit)
		/* do nothing */;
	else
	if (zone->region()->type() == REGION_BLACKHOLE) // Blackhole
		return false;
	else
	if (zone->region()->type() == REGION_VOXEL) {
		// Inside voxel region
		hit = segment.hit(SMALL3D3);
		ray->voxelreg = viewer.voxel().get(hit);
		if (ray->voxelreg<0 || viewer.voxel().color(ray->voxelreg)==COLOR_TRANSPARENT) {
			// End of the voxel or precision error
			ray->moveby(SMALL3D3);
			goto INT;		// treat as transparent
		} else {
			if (segment.body && ray->clip_hit) {
				ray->normal = segment.body->body()->normal(hit);
				if (ray->normal*segment.dir > 0.0)
					ray->normal.negate();
			} else
			if (!ray->start())
				ray->normal = geometry.voxel.normal(hit, segment.dir);
		}
	} else {
		// Normal region
		// Check if we are at the entry of a lattice
		// in order to avoid unnecessary transformations
		while (n>0 && ray->segment(n).body == NULL) n--;
		if (ray->segment(n).body != NULL) {
			hit = ray->segment(n).hit(-SMALL3D3);

			// Find body normal
			ray->normal = ray->segment(n).body->body()->normal(hit);

			// check that the normal is opposite to the direction
			if (ray->normal*ray->segment(n).dir > 0.0)
				ray->normal.negate();
		}
	}
	// Transform normal due to lattice if any
	ray->normal.normalize();

	if (n>0) {
		STATS(stats.shade_correctRotdefiNormal_count++);
		for (; n>=0; n--)
		  //			if (ray->segment(n).rotdefi)
//				ray->normal = geometry.invRotdefi(
//						ray->segment(n).rotdefi).multVector(
//							ray->normal);
			if (ray->segment(n).region)
				ray->normal = ray->segment(n).invMatrix() * ray->normal;
		//				ray->normal = geometry.invRotdefi(ray->segment(n).rotdefi) * ray->normal;
	}

	return true;
} // nextIntersection

/** shade
 * returns the final color of the ray
 */
Color D3Layer::shade(GeometryEngine* eng, Ray *ray)
{
	Point hit;
	Color color;
	const Material* material = NULL;

	STATS(stats.shade_call_count++);

INT:	if (!nextIntersection(eng, ray))
		return color;	// on error or blackbody

	// Get region/voxel color
	if (ray->voxelreg > 0)
		color.set(viewer.voxel().roiColor(ray->voxelreg));	// voxel color
	else {
		color.set(ray->hitRegion()->color());			// region color
		material = ray->hitRegion()->region()->material();
	}

	// Check for USRBIN
	if (viewer.usrbin.show) {
		if (usrbinAsTexture() || ray->project_hit) {
//			Color32 col;
//			col.val = 0;
			ray->project_alpha = viewer.usrbin.shade(ray->hit(SMALL3D3), color);
//			color.set(col);
			if (ray->project_hit && ray->project_alpha == 0xFF) {
				ray->project_hit = false;
				ray->moveby(SMALL3D3);
				goto INT;
			}
		}
	}

	// If we are in the start of the plane do not apply any shadow
	if (ray->start()) return color;

	if (material && material->fuzz()>0.0) {
		ray->normal += material->fuzz() *
			Vector(eng->random(-1.0,1.0), eng->random(-1.0,1.0), eng->random(-1.0,1.0));
		ray->normal.normalize();
	}

	// All light rays should start from outside the object
	// Move a bit backwards
	// Should be bigger than what is inside intersectRay
	hit = ray->hit(-SMALL3D3);

	// Calculate illumination
	double intensity = 0.0;
	for (int i=0; i<Min(ray->lights, geometry.lights); i++) {
		const Light& L = geometry.light[i];	// Geometry Light properties
		VLight& VL = light[i];			// Viewer light coordinate & direction
		Vector dir;
		double dist, power;
		switch (L.type) {
			case LIGHT_SUN:
				dir   = -VL.dir;
				dist  = INFINITE;
				power = L.power;
				if (L.shadow && ray->shadow) {
					double dot = view().matrix(0,2)*dir.x +
						     view().matrix(1,2)*dir.y +
						     view().matrix(2,2)*dir.z;
					if (dot>SMALL3D3 && view().projection==Projection_Orthographic) {
						dist = (view().matrix(0,2)*(view().matrix(0,3)-hit.x) +
							view().matrix(1,2)*(view().matrix(1,3)-hit.y) +
							view().matrix(2,2)*(view().matrix(2,3)-hit.z)) /
							dot;
						if (dist < SMALL3D3) continue;
					}
				}
				break;

			case LIGHT_OMNI:
			case LIGHT_SPOT:
				dir  = VL.pos - hit;
				dist = dir.normalize();
				if (dist<SMALL3D3) continue;
				if (dist>L.distance) continue;
				if (L.falloff) {
					if (L.falloff==1)
						power = L.power * 100.0 / dist;
					else
						power = L.power * 10000.0 / Sqr(dist);
					if (power<0.001) continue;
				} else
					power = L.power;
				if (L.shadow && ray->shadow) {
					double dot = view().matrix(0,2)*dir.x +
						     view().matrix(1,2)*dir.y +
						     view().matrix(2,2)*dir.z;
					if (dot>SMALL3D3) {
						double pdist = (view().matrix(0,2)*(view().matrix(0,3)-hit.x) +
							        view().matrix(1,2)*(view().matrix(1,3)-hit.y) +
							        view().matrix(2,2)*(view().matrix(2,3)-hit.z)) /
							        dot;
						dist = Min(dist, pdist);
						if (dist < SMALL3D3) continue;
					}
				}
				break;

			default:
				continue;
		}
		double dot = dir * ray->normal;
		if (dot>0.00001) {
			if (L.shadow && ray->shadow) {
				// Find intersection with light
				Ray lightRay;
				lightRay.use_clip = ray->use_clip;
				RaySegment new_segment(hit, dir, ray->prevZone());
				new_segment.tmin = SMALL3D3*ray->segment().tmin;
				new_segment.tmax = dist;
				lightRay.push(new_segment);

				int oldBodyCheckId = eng->incBodyCheckId();
				if (!eng->intersectRay(&lightRay, false)) {
					// Check if we hit a black hole or an error
					VZone *z = lightRay.hitZone();
					if (z!=NULL && (z->region()->type()==REGION_NORMAL ||
					    z->region()->type()==REGION_VOXEL)) {
						eng->bodyCheckId(oldBodyCheckId);
						continue;
					}
				}
				eng->bodyCheckId(oldBodyCheckId);
			}
			intensity += dot * power;
			if (L.spec && material && material->shine()>0.0) { //if (L.spec>0.001) {
				dot = ray->dir() * ray->normal;
				Vector reflect = (-2.0*dot)*ray->normal + ray->dir();
				double spec = dir * reflect;
				if (spec > 0.0)
					intensity += power*pow(spec, material->shine());
			}
		}
	}

#if 0
	int log2depth;
	frexp(hit.z, &log2depth);
#endif
	color *= intensity + (float)ambient/255.0;
	if (ray->depth < ray->max_depth &&
	    material!=NULL &&
	    material->specular() > 0.0 ) {
		double dot = ray->normal * ray->dir();
		if (dot<0.0) {
			Vector reflect = (-2.0*dot)*ray->normal + ray->dir();
			// Find intersection with light
			Ray reflRay(*ray);
			reflRay.push(RaySegment(hit, reflect, ray->prevZone()));
			int oldBodyCheckId = eng->incBodyCheckId();
			color += shade(eng, &reflRay) * material->specular();
			eng->bodyCheckId(oldBodyCheckId);
		}
	}
	return color;
} // shade

/** shadeXray
 * Computes the color of a pixel
 * for transparent objects it invokes iteratively shade while the color contribution
 * is above the threshold
 */
dword D3Layer::shadeXray(GeometryEngine* eng, Ray* ray,
			const double u, const double v,
			int alpha, VRegion* last_region)
{
	// alpha=255 first level

	VZone *zone = ray->segment(0).zone;
	while (true) {
		Color32 thisColor, nextColor, finalColor;
		int nextAlpha;
		int transparency;
		bool shadowsChanged;

		// Remove shadows for deeper levels
		if (alpha!=255 && ray->shadow) {
			ray->shadow = false;
			shadowsChanged = true;
		} else
			shadowsChanged = false;

		thisColor.val = shade(eng, ray).color32();

		// restore shadows
		if (shadowsChanged) ray->shadow = true;

		if (alpha==255 && _edgeDetect && !ray->project_hit) {
			// remember checkid to restore in case of no edge
			int oldBodyCheckId = eng->bodyCheckId();
			if (edgePixel(eng, *ray, u,v, zone)) return 0;
			eng->bodyCheckId(oldBodyCheckId);
		}
		thisColor.rgb.alpha = 0;

		if (ray->ended() ||
		   !ray->hitZone() ||
		    ray->hitRegion()->type() == REGION_BLACKHOLE)
			return geometry._backgroundColor;
		else
		if (ray->hitRegion() == last_region && \
		    ray->hitRegion()->type() != REGION_VOXEL) {
			ray->moveby(SMALL3D3);
			ray->skip_current = true;
			continue;
		}

		if (ray->project_hit)
			transparency = ray->project_alpha;
		else
			transparency = (xray>0)? xray: ray->hitRegion()->alpha();

		nextAlpha = byteMul(alpha, transparency);
		if (nextAlpha < 1)
			return thisColor.val;

		ray->skip_transparent = !ray->skip_transparent;
		nextColor.val = shadeXray(eng, ray, u, v, nextAlpha, ray->hitRegion());
		nextColor.rgb.alpha = 255;
		thisColor.rgb.alpha = 255 - transparency;

		finalColor.val = alphaBlend(thisColor, nextColor);
		// "Fix" color.rgb.alpha
		finalColor.rgb.alpha = 0;

		return finalColor.val;
	}
} // shadeXray

/** edgePixel
 * @return true if pixel is on edge
 */
bool D3Layer::edgePixel(GeometryEngine* eng, Ray& ray, const double u, const double v, VZone* zone)
{
	// Sample few rays to detect change in depth, normal, region
	Point  hit[5];		// hit position
	Point  projhit[5];	// hit position projected to the hit-plane of center pixel
	Vector normal[5];	// normal at hit position
	CBody* body[5];

	// Do not return edge for zone=NULL or BLACKBODY
	if (!ray.hitZone() || ray.hitRegion()->type() == REGION_BLACKHOLE) return false;

	double depth = ray.T();

	hit[0]    = projhit[0] = ray.hit();
	normal[0] = ray.normal;
	body[0]   = ray.hitBody();

	Ray edgeray;
	edgeray.use_clip      = ray.use_clip;
	edgeray.skip_1stblack = ray.skip_1stblack;

	int cid = 1;

	double du = view().Dx();
	double dv = view().Dy();

	double x=0.0,y=0.0,z=0.0;

	if (view().projection == Projection_Perspective)
		view().cameraPosition(&x, &y, &z);

	double uu = u - du/2.0;
	for (int i=0; i<2; i++) {
		double vv = v - dv/2.0;
		for (int j=0; j<2; j++) {
			double dx, dy, dz;
			if (view().projection!=Projection_Perspective) {
				x = view().uv2x(uu,vv);
				y = view().uv2y(uu,vv);
				z = view().uv2z(uu,vv);
			}
			view().rayDirection(uu,vv, &dx,&dy,&dz);

			eng->incBodyCheckId();
			zone = eng->whereRay(x,y,z, dx,dy,dz, SMALL3D, zone);
			edgeray.init();
			edgeray.push(RaySegment(x,y,z, dx,dy,dz, zone));

			nextIntersection(eng, &edgeray);

			if (edgeray.hitZone()) {
				if (edgeray.hitRegion()->type() == REGION_VOXEL) {
					if (edgeray.voxelreg != ray.voxelreg) return true;
				} else
				if (edgeray.hitRegion() != ray.hitRegion()) return true;
			}

			hit[cid]     = edgeray.hit();
			projhit[cid] = edgeray.dir() * depth + edgeray.pos();
			normal[cid]  = edgeray.normal;
			body[cid]    = edgeray.hitBody();

			for (int k=0; k<cid; k++) {
				double dot = normal[cid] * normal[k];
				if (body[k]==body[cid] && dot<0.5) return true;
				else
				if (body[k]!=body[cid] && dot<0.999999999) return true;
				else
				if (dot>0.99999) {
					double cosf = normal[cid]*ray.dir();
					double be = (hit[cid]-projhit[k]).length2();
					double bd = (hit[cid]-hit[k]).length2() * Sqr(cosf);
					if (bd > 2.0 * be) return true;
				}
			}
			cid++;
			vv += dv;
		}
		uu += du;
	}

	return false;
} // edgePixel

/** drawPixel */
VZone* D3Layer::drawPixel(GeometryEngine* eng, Ray& ray,
			dword *ptr, int W, int H,
			int i, int j, int step, int sample,
			double u, double v,
			VZone* zone)
{
	double x,y,z;
	double dx,dy,dz;

	view().rayPosition(u,v, &x, &y, &z);
	view().rayDirection(u,v, &dx, &dy, &dz);

	eng->incBodyCheckId();
	zone = eng->whereRay(x,y,z, dx,dy,dz, SMALL3D, zone);

	// prepare the ray
	ray.init();
	ray.push(RaySegment(x,y,z, dx,dy,dz, zone));

	dword color = shadeXray(eng, &ray, u, v);

	int selected = 0;
	dword color2 = 0;
	if (ray.error && _showErrors && ((i+j)&15)==0) {
		color = geometry.errorColor;
		kernel.setError();
	} else
	if (ray.hitZone() && ray.hitRegion()->show() & BIT_SELECT) {
		Color32 col;
		col.val = ray.hitRegion()->color();
		col.rgb.red   = col.rgb.red  >0x7F? 0:0xFF;
		col.rgb.green = col.rgb.green>0x7F? 0:0xFF;
		col.rgb.blue  = col.rgb.blue >0x7F? 0:0xFF;
		color2 = col.val | FLAG_3D;
		selected = 1;
	}

	// FIXME not the best way but works
	// Maybe I should add a selected flag inside the zone
	if (geometry.editRegion().nzones()) {
		Point hit = ray.hit();
		geometry.lockEdit();
		if (geometry.editRegion().inside(hit.x,hit.y,hit.z, dx,dy,dz)) {
			color2   = 0xFF00FF | FLAG_3D;
			selected = 2;
		}
		geometry.unlockEdit();
	}

	color |= FLAG_3D;
	if (color & FLAG_TRANSPARENT)
		if ((i^j)&8)
			color = Darken(color, 200);

	dword *p = ptr;
	if (step) {
		for (int jj=0; jj<step; jj++, p+=W) {
			if (j+jj>=H) break;
			int step2 = Min(step, W-i);
			for (int ii=0; ii<step2; ii++) {
				dword flag = p[ii] & FLAG_INFOMASK;
				if (flag && flag!=FLAG_REGION && flag!=FLAG_ERROR) {
					Color32 new_color;
					Color32 mean;

					if ((selected==1 && (!(((i+ii+j+jj)&3) || ((j+jj)&1)))) ||
					    (selected==2 && ((i+ii+j+jj)&1)))
						new_color.val = color2;
					else
						new_color.val = color;
					mean.val = p[ii];

					Color32 delta;
					delta.rgb.red   = (new_color.rgb.red   - mean.rgb.red)  / sample;
					delta.rgb.green = (new_color.rgb.green - mean.rgb.green)/ sample;
					delta.rgb.blue  = (new_color.rgb.blue  - mean.rgb.blue) / sample;

					new_color.rgb.red   = mean.rgb.red   + delta.rgb.red;
					new_color.rgb.green = mean.rgb.green + delta.rgb.green;
					new_color.rgb.blue  = mean.rgb.blue  + delta.rgb.blue;

					p[ii] = new_color.val;
				}
			}
		}
	}

	return zone;
} // drawPixel

/** draw */
void D3Layer::draw(Painter& painter, bool checkTime)
{
	// FIXME ray... Should move to workers
	Ray ray;

	STATS(stats.reset());

	// Check for the default lights
//XXX	if (geometry.lights==0 || deflights) geometry.defaultLights();

	clock_t starttime, endtime;
	if (checkTime && _drawTime) {
		starttime   = clock();
		endtime     = starttime + _drawTime;	// find ending time
		ray.shadow  = false;
		ray.lights  = 2;

		_edgeDetect = false;
		_showErrors = false;
		_samples    = 1;
	} else {
		starttime   = 0;
		endtime     = 0;
		ray.shadow  = shadows;
		ray.lights  = geometry.lights;

		_edgeDetect = drawEdges;
		_showErrors = viewer.showErrors;
		_samples    = antialias;
	}
	ray.use_clip      = kernel.clipBodyCount()    > 0;
	ray.use_project   = kernel.projectBodyCount() > 0;
	ray.max_depth     = _maxDepth;
	ray.skip_1stblack = skip_1stblack;

	// Convert lights from relative to absolute position and direction
	for (int l=0; l<geometry.lights; l++) {
		if (geometry.light[l].relative) {
			light[l].pos = view().matrix() * geometry.light[l].pos;	// Point multiplication
			light[l].dir = view().matrix() * geometry.light[l].dir;	// Vector multiplication
//			light[l].dir = view().matrix().multVector(geometry.light[l].dir);
		} else {
			light[l].pos = geometry.light[l].pos;
			light[l].dir = geometry.light[l].dir;
		}
	}

	geometry.lockRead();
#ifdef THREAD
	if (kernel.nthreads()) {
		if (kernel.nGeometryBodies()) {
			assert(!kernel.isRunning());
			feeder.reset(&painter, ray);
			for (int step=_step; step>0; step>>=1) {
				feeder.init(step);
				if (kernel.threadpool.execute(&feeder)) {
//std::cout << "Stop " << (double)(clock() - starttime)/(double)CLOCKS_PER_SEC << "s  check=" << checkTime << std::endl;
					break;
				}
				// Stop if we have eaten up our time
				if (endtime && clock()>endtime) break;
			}
		}
	} else
#endif
	{
		int H = painter.height();	// save variables that might change
		int W = painter.width();	// during the scanning...
		VZone *zone = NULL;

		// scan in steps
		for (int step=_step; step>0; step>>=1) {
			dword *pline = painter.data();
			for (int j=0; j<H; j+=step) {
				dword *ptr = pline;
				for (int i=0; i<W; i+=step, ptr+=step) {
					dword flag = *ptr & FLAG_INFOMASK;
					if (flag==0 || flag==FLAG_REGION) continue;
					int sample = 1;
					for (int sx=0; sx<_samples; sx++, sample++) {
						double u = view().i2u((double)i + (sx+0.5)/(double)_samples);
						for (int sy=0; sy<_samples; sy++) {
							if (stop()) {
								geometry.unlockRead();
								return;
							}
							double v = view().j2v((double)j + (sy+0.5)/(double)_samples);
							zone = drawPixel(engine(), ray,
									ptr, W, H,
									i, j, step, sample,
									u, v,
									zone);
						}
					}
				}
				pline += step*W;
			}
			// Stop if we have eaten up our time
			if (endtime && clock()>endtime) break;
		}
	}
	geometry.unlockRead();

	// adaptive step adjust _step if needed
	if (endtime) {
		clock_t dt = clock() - starttime;
		if (dt > _drawTime*4) {		// _step goes with the square (surface)
			_step <<= 1;
			if (_step>64) _step = 64;
		} else
		if (dt < _drawTime/4) {
			_step >>= 1;
			if (_step==0) _step = 1;
		}
	}
//std::cout << "End  " << (double)(clock() - starttime)/(double)CLOCKS_PER_SEC << "s  check=" << checkTime << std::endl;

	STATS(std::cout << "Stats\n" << stats << std::endl);
} // draw

/** draw the projection of a 3D line */
bool D3Layer::draw3Dline(Painter& painter,
			const Point& a, const Point& b,
			const Color3D& color)
{
	Point ra, rb;

	if (view().clipLine3D(a, b, &ra, &rb)) {
		dword col;
//		if (selected)
//			color = geometry->wireframeColor.select;
//		else
		if (ra.z>SMALL3D && rb.z>SMALL3D) {
			// FIXME normally should be broken into two parts >0 and <0
			col = color.bright;
		} else
		if (ra.z<-SMALL3D && rb.z<-SMALL3D)
			col = color.color;
		else
			col = color.dark;

		int x1 = view().u2i(ra.x);
		int y1 = view().v2j(ra.y);
		int x2 = view().u2i(rb.x);
		int y2 = view().v2j(rb.y);

		if (painter.line(x1,y1, x2,y2, col))
			return true;
	}
	return false;
} // draw3Dline

/** drawWireframe */
void D3Layer::drawWireframe(Painter& painter)
{
	geometry.lockRead();
	for (int i=0; i<kernel.nGeometryBodies(); i++) {
		VBody* body = kernel.getBody(i);
		if (body->show() & BIT_WIREFRAME)
			drawWireframe(painter, body);
	}
	geometry.unlockRead();
} // drawWireframe

/** drawWireframe */
void D3Layer::drawWireframe(Painter& painter, VBody *body)
{
	const Color3D& color = (body->show() & BIT_SELECT)?
			geometry.select3DColor : geometry.wireframeColor;

	/* let's use brute force */
	GBody* gbody = body->body();
	Mesh& mesh = gbody->mesh;
	for (int i=0; i<mesh.nedges(); i++){
		const Edge* edge = mesh.edge(i);
		if (edge->show && draw3Dline(painter, edge->A(), edge->B(), color))
			body->visible = true;
	}

	for (int n=0; n<gbody->nodes(); n++) {
		Point r;
		view().xyz2uvw3D(gbody->node(n), &r);
		if (show || Eq0(r.z,SMALL3D)) {
			int x = view().u2i(r.x);
			int y = view().v2j(r.y);
			if (painter.rectangle(x-1,y-1,x+1,y+1, color.color))
				body->visible = true;
		}
	}

#ifdef OPENGL
#if 0
	// FIXME VERY PRIMITIVE
	gbody->mesh.clearProcessedFlag();
	for (int i=0; i<gbody->mesh.faceCount(); i++) {
		Face* f = gbody->mesh.face(i);
		for (int j=0; j<3; j++) {
			if (!f->show(j)) continue;
			Face* fn = f->neighbor(j);
			if (fn && fn->processed()) continue;
			int jn = (j+1)%3;
			Point* A = f->vertex(j);
			Point* B = f->vertex(jn);
			draw3Dline(painter, *A, *B, color);
		}
		f->processed(true);
	}
#endif
#endif
} // drawWireframe

/** drawBBox */
void D3Layer::drawBBox(Painter& painter, const BBox& bbox, const Color3D& color)
{
	if (!bbox.isValid()) return;

	for(int i=0; i<bbox.edges(); i++) {
		int a, b;
		bbox.edge(i, &a, &b);
		draw3Dline(painter, bbox.vertex(a), bbox.vertex(b), color);
	}
} // drawBBox

/** drawOBBox */
void D3Layer::drawOBBox(Painter& painter, const OBBox* bbox, const Color3D& color)
{
	if (!bbox->isValid()) return;

	for(int i=0; i<bbox->edges(); i++) {
		int a, b;
		bbox->edge(i, &a, &b);
		draw3Dline(painter, bbox->vertex(a), bbox->vertex(b), color);
	}
} // drawOBBox

/** drawBodiesBBox */
void D3Layer::drawBodiesBBox(Painter& painter)
{
	geometry.lockRead();;
	for (int i=0; i<kernel.nGeometryBodies();i++) {
		VBody* body = kernel.getBody(i);
		if (body->show() & BIT_BBOX) {
			drawBBox(painter, body->body()->bbox(),
					geometry.bodyBBoxOutColor);
			drawOBBox(painter, body->body()->obbox(false),
					geometry.bodyBBoxOutColor);
			drawOBBox(painter, body->body()->obbox(true),
					geometry.bodyBBoxInColor);
		}
	}
	geometry.unlockRead();;
} // drawBodiesBBox

/** drawRegionsBBox */
void D3Layer::drawRegionsBBox(Painter& painter)
{
	geometry.lockRead();;
	// FIXME scan regions before if BBOX is selected
	ArrayIterator<GRegion*> iter_regions(geometry.regions);
	while (iter_regions) {
		GRegion *region = iter_regions++;
		if (region->show & BIT_BBOX) {
			drawBBox(painter,  region->bbox(),  geometry.regionBBoxColor);
			drawOBBox(painter, region->obbox(), geometry.zoneBBoxColor);
		}
	}
#if 0
	geometry.lockEdit();
	if (geometry.editRegion().show&BIT_BBOX && geometry.editRegion().nzones()) {
		ArrayIterator<GZone*> ziter(geometry.editRegion().zones());
		while (ziter) {
			GZone* z = ziter++;
			drawBBox(painter,  z->bbox(),  geometry.zoneBBoxColor);
			drawOBBox(painter, z->obbox(), geometry.zoneBBoxColor);
		}
	}
	geometry.unlockEdit();
#endif
	geometry.unlockRead();;
} // drawRegionsBBox

/* ========================== D3LayerStats ========================== */
#ifdef _STAT
/** reset drawing statistics */
void D3LayerStats::reset()
{
	shade_call_count                 = 0;
	shade_correctRotdefiNormal_count = 0;
} // reset

/** operator << */
std::ostream & operator<<(std::ostream &os, const D3LayerStats &gv)
{
	os << " shade call count:        " << gv.shade_call_count                 << std::endl;
	os << " shade corrected normals: " << gv.shade_correctRotdefiNormal_count << std::endl;
	return os;
} /* operator << */
#endif
