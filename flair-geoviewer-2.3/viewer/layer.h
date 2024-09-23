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

#ifndef __LAYER_H
#define __LAYER_H

class Geometry;
class GeometryViewer;
class GeometryKernel;
class GeometryEngine;
class Painter;
class ViewPort;

class Layer {
public:
const	Geometry&	geometry;
	GeometryKernel&	kernel;
	GeometryViewer&	viewer;

public:
	int	show;		/** show current layer */

public:
	Layer(const Geometry& g, GeometryKernel& k, GeometryViewer& v);
	ViewPort&	view()	const;
	bool		stop()	const;

	GeometryEngine*	engine() const;
#ifdef THREAD
	GeometryEngine*	engine(int i) const;
#endif
}; // Layer
#endif
