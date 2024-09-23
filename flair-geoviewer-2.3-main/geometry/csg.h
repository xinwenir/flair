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

#ifndef __CSG_H
#define __CSG_H

#include <iostream>

#include "gbody.h"
#include "array.h"

class GBody;
class Geometry;

/* =================================== CSG ================================== */
/** Constructive solid geometry class */
class CSG {
public:
	Array<GBody*>	expr;
private:
	Geometry&	geometry;

public:
	CSG(Geometry& g) : geometry(g) {}

	~CSG() {}

	void	clear()		{ expr.clear(); }
	void	exp2rpn();
	void	rpn2exp();
	void	tokenize(const char *expstr);
	void	rpnorm();

	bool	depth();
private:
	int	_rpnrule(int n);
	void	_subTerms(int n, int* lowLeft, int* lowRight);
	void	_copy(int dst, int src, int length);
static	bool	_isPlus(const GBody* t)		{ return t==&GBody::tplus || t==&GBody::tminus;}
static	int	_priority(const GBody* body, int* ip=NULL);
}; // CSG

std::ostream& operator << (std::ostream&, const CSG&);

#endif
