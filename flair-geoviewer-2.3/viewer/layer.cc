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

#include "layer.h"

#include "engine.h"
#include "kernel.h"
#include "viewer.h"

// XXX XXX XXX
// FIXME why I cannot inline them?
// maybe make them as MACROS
// XXX XXX XXX

/** Layer */
Layer::Layer(const Geometry& g, GeometryKernel& k, GeometryViewer& v) :
	geometry(g),
	kernel(k),
	viewer(v),
	show(0)
{
} // Layer

/** view */
ViewPort& Layer::view() const
{
	return kernel.view;
} // view

/** engine */
GeometryEngine* Layer::engine() const
{
	return kernel.engine();
} // engine

#ifdef THREAD
/** engine */
GeometryEngine* Layer::engine(int i) const
{
	return kernel.engine(i);
} // engine
#endif

/** stop */
bool Layer::stop() const
{
	return viewer.stop();
} // stop
