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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sstream>

#include "csg.h"
#include "geometry.h"
#include "token.h"

using namespace std;

/* =================================== CSG ================================== */
/** return token input and output
 * @param ip input priority
 * @return output priority
 */
//               Operator In Out
//priorities = {"(" : [ 99,  0],
//		"|" : [  1,  2],
//		"+" : [  3,  4],
//		"-" : [  3,  4],
//		")" : [  0, 99],
//		" " : [101,100]}
int CSG::_priority(const GBody* body, int* ip)
{
	// Find priorities
	if (body==&GBody::tleft) {
		if (ip) *ip = 99;
		return 0;
	} else
	if (body==&GBody::tunion) {
		if (ip) *ip =  1;
		return 2;
	} else
	if (body==&GBody::tplus || body==&GBody::tminus) {
		if (ip) *ip =  3;
		return 4;
	} else
	if (body==&GBody::tright) {
		if (ip) *ip =  0;
		return 99;
	} else {
		if (ip) *ip = 101;
		return 100;
	}
} // _priority

/** @return depth of the rpn expression */
bool CSG::depth()
{
	int d = 0;
	for (int i=0; i<expr.count(); i++) {
		if (expr[i]->isOperator())
			d--;
		else
			d++;
	}
	return d==1;
} // depth

/* ---------------------------------------------------------------------- *
 * exp2rpn
 *
 * Convert the FLUKA Boolean expression to Reverse Polish Notation (RPN)
 * Since the normalization routine does not accept the +/- as signs to
 * objects, the routine is converting the leading - to @- (Universe minus)
 * where the special symbol @ is treated as the universe.
 * Furthermore the leading + is ignored as well as the leading | which is
 * accepted by fluka.
 *
 * WARNING: No check is done for the correctness of the expression apart
 *          from the parenthesis nesting
 *
 * ie.
 *        A+B         -> A B +
 *       (A+B)|C      -> A B + C |
 *       (A|B)|C+D|E  -> A B | C D + | E |
 *       -A           -> @ A -
 *       -(-A)        -> @ @ A - -
 *
 * The routine is using the same array for returning the Reverse Polish
 * expression, since the format is more compact.
 * This is generally true apart one case -A -> @ A -
 *
 * Priorities are treated as:
 *     Operator  Priority         In     Out
 *     --------  --------        ---     ---
 *       |       lower             1       2
 *       +       high              3       4
 *       -       high              3       4
 *       (       higher           99       0
 *       )       higher            0      99
 *       object  highest         101     100
 *
 * Algorithm
 *  Consider the expression as a train moving on a railroad with a
 *  T-shape, where each token is one wagon
 *
 *                      <-  (A|B)|C+D
 *      ------------.   .------------
 *      RPN-End      \ /      Exp-End
 *                    |
 *                   S|
 *                   t|
 *                   a|
 *                   c|
 *                   k|
 *
 *  Each wagon to move from the Exp-End to the RPN-End it has to make
 *  a stop first in the Stack-End. Before entering in the stack, the
 *  priority-IN will be checked against the objects in the stack.
 *  All top-most objects currently present in the stack with
 *  higher priority-OUT will be transfered from the stack to the
 *  RPN-End. Apart from the opening parenthesis ( which is discarded.
 *
 *  Example:
 *  (1)                         A+B|C (1)
 *  (2) A                        +B|C (2)
 *  (3) A                         B|C (3)
 *  (4) A                          |C (4)
 *  (5) A B +                       C (5)
 *  (6) A B + C
 *  (7) A B + C |
 *      ------------.   .------------
 *      RPN-End      \ /      Exp-End
 *                    |
 *                   S|
 *                   t|
 *                   a|
 *                   c|     B   C
 *                   k| A + + | | |
 *                      1 2 3 4 5 6 7
 */
void CSG::exp2rpn()
{
	bool newproduct = true;
	int i = 0;
	int m = 0;
	Array<GBody*> stack;

	while (i < expr.count()) {
		GBody* body = expr[i];

		//# Check for special leading chars
		if (newproduct && (body==&GBody::tunion || body==&GBody::tplus)) {
			newproduct = (body==&GBody::tunion);
			i++;
			continue;
		}

		if (newproduct && body==&GBody::tminus) {
			expr.insert(i, &GBody::tnull);	// insert space in ith position
			body = &GBody::tuniverse;
		}

		newproduct = (body==&GBody::tleft || body==&GBody::tunion);

		int ip;
		_priority(body,&ip);

		// Remove from the stack everything with higher priority
		while (!stack.empty() && ip<_priority(stack.tail())) {
			if (stack.tail() != &GBody::tleft) {
				expr[m] = stack.tail();
				m++;
			}
			stack.pop();
		}

		//# Push it into the stack
		if (body != &GBody::tright)
			stack.append(body);
		else {
			// Should be an opening parenthesis
			if (stack.empty()) {
				throw 0;
				//raise CSGException("Unbalanced parenthesis")
			}
			stack.pop();
		}
		i++;
	}

	// Empty Stack
	while (!stack.empty()) {
		if (stack.tail() == &GBody::tleft) {
			throw 0;
			//CSGException("Unbalanced parenthesis")
		}

		if (m >= expr.count())
			expr.append(stack.tail());
		else
			expr[m] = stack.tail();
		stack.pop();
		m++;
	}

	// Delete unwanted items
	expr.erase(m, expr.count());
	//for i in xrange(len(expr)-1,m-1,-1):
	//	del expr[i]
} // exp2rpn

/** ---------------------------------------------------------------------- *
 * rpn2exp
 *
 * Convert a NORMALIZED Reverse Polish notation to a standard expression
 * WARNING: The routine expects an expression where for each UNION
 * operator the right-sub-expression is a product while the left can be
 * UNION or a product
 * ---------------------------------------------------------------------- */
void CSG::rpn2exp()
{
	Array<GBody*> newexpr(expr.count());
	Array<GBody*> plus;
	Array<GBody*> minus;
//	plus.compare(GBody::compare);
//	minus.compare(GBody::compare);

	int i = 0;
	int nstack = 0;
	bool endprod = false;
	while (i<expr.count()) {
		GBody* body = expr[i];
		bool lastPlus = true;

		if (body->isOperator())
			nstack--;
		else {
			// First term is always a plus
			// .. peek then next operator to check for sign
			if (plus.empty() || expr[i+1] == &GBody::tplus) {
				plus.append(expr[i]);
				lastPlus = true;
			} else
			if (expr[i+1] == &GBody::tminus) {
				minus.append(expr[i]);
				lastPlus = false;
			}
			nstack++;
		}

		if (nstack==0)
			endprod = true;
		else
		if (nstack==3) {
			if (lastPlus)
				plus.pop();
			else
				minus.pop();
			i -= 2;
			endprod = true;
		} else
		if (body==&GBody::tunion) {
			i -= 2;
			endprod = true;
		} else
		if (i==expr.top())
			endprod = true;

		if (endprod) {
			//optZone(plus,minus)
//			if (!plus.empty() || !minus.empty())
//				zones.append((plus,minus));
			bool something_added = false;
			for (int ii=0; ii<plus.count(); ii++)
				if (plus[ii]->name()[0] != ' ') {
					newexpr.append(&GBody::tplus);
					newexpr.append(plus[ii]);
					something_added = true;
				}
			for (int ii=0; ii<minus.count(); ii++)
				if (minus[ii]->name()[0] != ' ') {
					newexpr.append(&GBody::tminus);
					newexpr.append(minus[ii]);
					something_added = true;
				}
			if (something_added)
				newexpr.append(&GBody::tunion);
			plus.clear();
			minus.clear();
			nstack  = 0;
			endprod = false;
		}
		i++;
	}

	// Remove last |
	if (!newexpr.empty()) newexpr.pop();

	// Copy back expression
	expr.clear();
	expr.allocate(newexpr.count());
	memcpy(&expr[0], &newexpr[0], sizeof(GBody*)*newexpr.count());

	//# Remove duplicates of products
//	rmDoubles(zones);

#if 0
	// Reconstruct expression
	expr = []
	for plus,minus in zones:
		if len(expr)>0 and expr[-1]!="|": expr.append("|")
		//# Fill the new array
		for j in plus:
			expr.append("+")
			expr.append(j)
		for j in minus:
			expr.append("-")
			expr.append(j)
#endif
//	return expr
} // rpn2exp

/** --------------------------------------------------------------------- *
 * rpnorm
 *
 * Normalize a CG expression given in Reverse Polish Notation.
 * Normalized CG expression is an expression given as sum (Boolean OR) of
 * products (Boolean intersection or subtraction).
 * The normalization (expansion of parenthesis and operator priorities)
 * should be performed by recursively calling the RPNRULE subroutine.
 * Since Fortran-77 doesn't have recursion, call the RPNRULE for every
 * operator starting from the right-most one, until no rule is found.
 * ----------------------------------------------------------------------
 */
void CSG::rpnorm()
{
	// Loop until there is no any extra change needed
	// Scan to find the first operators
	bool changed;
	do {
		changed = false;
		int i = expr.count()-1;
		while (i >= 4) {
			if (expr[i]->isOperator()) {
				int l = expr.count();
				int rule = _rpnrule(i);
				if (rule>0) {
					changed = true;
					i += expr.count()-l+1;
				}
			}
			i--;
		}

		// Find possible products and optimize them
		if (changed) {
			// Product target should be in the form
			// ... t1 t2 [+-] [ti [+-]]* ...
			//     i               iend
			i = 0;
			while (i<expr.count()-2) {
				if (expr[i]->isOperator()) {
					i++;
					continue;
				}
				// Find as many repetition of ti [+-]
				int iend = i + 1;
				while (iend<expr.count()-1) {
					if (!(expr[iend]->isOperator()) && _isPlus(expr[iend+1]))
						iend += 2;
					else
						break;
				}
				if (iend != i+1) {
					// Product found
					cout << "\tProduct= ";
					for (int ii=i; ii<iend; ii++)
						cout << expr[ii]->name() << ' ';
					cout << endl;
				}
				i = iend;
			}
		}
	} while (changed);
} // rpnorm

/** --------------------------------------------------------------------- *
 * rpnrule
 *
 * Find a matching rule and apply it on the sub-expression starting from
 * the N position in the Reverse Polish Notation
 *
 * An expression is in normal form when all the parenthesis are expanded
 * and the expression is described as a sum (UNIONS) of products
 * (INTERSECTIONS and/or SUBTRACTIONS)
 *
 * An expression can be converted to normal form by repeatedly applying
 * the following set of production rules to the expression and then to its
 * sub-expressions:
 *
 *    Normal Form                        Reverse Polish Notation
 *  1. X-(Y|Z) -> (X-Y)-Z                X Y Z | -  ->  X Y - Z -
 *  2. X+(Y|Z) -> (X+Y)|(X+Z)            X Y Z | +  ->  X Y + X Z + |
 *  3. X-(Y+Z) -> (X-Y)|(X-Z)            X Y Z + -  ->  X Y - X Z - |
 *  4. X+(Y+Z) -> (X+Y)+Z                X Y Z + +  ->  X Y + Z +
 *  5. X-(Y-Z) -> (X-Y)|(X+Z)            X Y Z - -  ->  X Y - X Z + |
 *  6. X+(Y-Z) -> (X+Y)-Z                X Y Z - +  ->  X Y + Z -
 *  7. X|(Y|Z) -> (X|Y)|Z                X Y Z | |  ->  X Y | Z |
 *  8. (X-Y)+Z -> (X+Z)-Y                X Y - Z +  ->  X Z + Y -
 *  9. (X|Y)-Z -> (X-Z)|(Y-Z)            X Y | Z -  ->  X Z - Y Z - |
 * 10. (X|Y)+Z -> (X+Z)|(Y+Z)            X Y | Z +  ->  X Z + Y Z + |
 * X,Y, and Z here match both primitives or sub-expressions.
 * ---------------------------------------------------------------------- */
int CSG::_rpnrule(int n)
{
	// Reset rule
	int rule = 0;

	// Top-most operator
	GBody* op = expr[n];
	if (!op->isOperator()) return rule;

	// Right operator
	GBody* rop = expr[n-1];

	// Find left and right sub-trees
	int ll, lr;
	_subTerms(n, &ll, &lr);

	// Left operator
	GBody* lop = &GBody::tnull;
	if (lr>0) lop = expr[lr-1];

	// Find Rule
	if (op==&GBody::tminus && rop==&GBody::tunion)
		rule = 1;
	else
	if (op==&GBody::tplus  && rop==&GBody::tunion)
		rule = 2;
	else
	if (op==&GBody::tminus && rop==&GBody::tplus)
		rule = 3;
	else
	if (op==&GBody::tplus  && rop==&GBody::tplus)
		rule = 4;
	else
	if (op==&GBody::tminus && rop==&GBody::tminus)
		rule = 5;
	else
	if (op==&GBody::tplus  && rop==&GBody::tminus)
		rule = 6;
	else
	if (op==&GBody::tunion && rop==&GBody::tunion)
		rule = 7;
	else
	if (op==&GBody::tplus  && lop==&GBody::tminus)
		rule = 8;
	else
	if (op==&GBody::tminus && lop==&GBody::tunion)
		rule = 9;
	else
	if (op==&GBody::tplus  && lop==&GBody::tunion)
		rule = 10;
	else
		return rule;

	// limits of subexpressions X,Y,Z
	int Xu, Xl, /*Yu,*/ Yl, Zu, Zl;
	int L;		// aux variables

	// Find sub expressions X Y Z
	if (rule<=7) {    // X op (Y rop Z)
		Xu = lr-1;
		Xl = ll;
		_subTerms(n-1, &ll, &lr);
		//Yu = lr-1;
		Yl = ll;
		Zu = n-2;
		Zl = lr;
	} else {	// (X lop Y) op Z
		Zu = n-1;
		Zl = lr;
		L  = lr-1;
		_subTerms(L, &ll, &lr);
		Xu = lr-1;
		Xl = ll;
		//Yu = L-1;
		Yl = lr;
	}

//	printf("Rule: %d\n",rule);
//	printf(" X[%d,%d,] = \n",Xl,Xu); //,"]=",expr[Xl:Xu+1]
//	printf(" Y[%d,%d,] = \n",Yl,Yu); //,"]=",expr[Yl:Yu+1]
//	printf(" Z[%d,%d,] = \n",Zl,Zu); //,"]=",expr[Zl:Zu+1]

	// Expand the rule
	switch (rule) {
		//# 1. X-(Y|Z) -> (X-Y)-Z	 X Y Z | -  ->  X Y - Z -
		case 1:
			// Leave X Y
			// Insert a - operator after Y
			expr.insert(Zl, &GBody::tminus);
			// Chop length by 1
			expr.erase(Zu+2);
			// Change the last operator to -
			expr[Zu+2] = &GBody::tminus;
			break;

		// 2. X+(Y|Z) -> (X+Y)|(X+Z)     X Y Z + |  ->  X Y + X Z + |
		case 2:
			// Leave X Y
			// Insert a + operator after Y
			expr.insert(Zl, &GBody::tplus);
			// Copy X after the + operator
			_copy(Zl+1, Xl, Xu-Xl+1);
			Zu += Xu-Xl+2;
			// Change last 2 operators to + |
			expr[Zu+1] = &GBody::tplus;
			expr[Zu+2] = &GBody::tunion;
			break;

		// 3. X-(Y+Z) -> (X-Y)|(X-Z)     X Y Z + -  ->  X Y - X Z - |
		case 3:
			// Leave X Y
			// Insert a - operator after Y
			expr.insert(Zl, &GBody::tminus);
			// Copy X after the - operator
			_copy(Zl+1, Xl, Xu-Xl+1);
			Zu += Xu-Xl+2;
			// Change last 2 operators to - |
			expr[Zu+1] = &GBody::tminus;
			expr[Zu+2] = &GBody::tunion;
			break;

		// 4. X+(Y+Z) -> (X+Y)+Z	 X Y Z + +  ->  X Y + Z +
		case 4:
			// Leave X Y
			// Insert a + operator after Y
			expr.insert(Zl, &GBody::tplus);
			// Chop length by 1
			expr.erase(Zu+2);
			// Change the last operator to +
			expr[Zu+2] = &GBody::tplus;
			break;

		// 5. X-(Y-Z) -> (X-Y)|(X+Z)     X Y Z - -  ->  X Y - X Z + |
		case 5:
			// Leave X Y
			// Insert a - operator after Y
			expr.insert(Zl, &GBody::tminus);
			// Copy X after the - operator
			_copy(Zl+1, Xl, Xu-Xl+1);
			Zu += Xu-Xl+2;
			// Change last 2 operators to + |
			expr[Zu+1] = &GBody::tplus;
			expr[Zu+2] = &GBody::tunion;
			break;

		// 6. X+(Y-Z) -> (X+Y)-Z	 X Y Z - +  ->  X Y + Z -
		case 6:
			// Leave X Y
			// Insert a + operator after Y
			expr.insert(Zl, &GBody::tplus);
			// Chop length by 1
			expr.erase(Zu+2);
			// Change the last operator to -
			expr[Zu+2] = &GBody::tminus;
			break;

		// 7. X|(Y|Z) -> (X|Y)|Z	 X Y Z | |  ->  X Y | Z |
		case 7:
			// Leave X Y
			// Insert a | operator after Y
			expr.insert(Zl, &GBody::tunion);
			// Chop length by 1
			expr.erase(Zu+2);
			// Change the last operator to |
			expr[Zu+2] = &GBody::tunion;
			break;

		// 8. (X-Y)+Z -> (X+Z)-Y	 X Y - Z +  ->  X Z + Y -
		case 8:
			// Leave X
			// Copy "Z +" after X
			L = Zu-Zl+2;
			_copy(Yl, Zl, Zu-Zl+2);
			// Delete old "Z +"
			expr.erase(Zl+L,Zl+L+L);	// XXX to be verified!
			break;

		// 9. (X|Y)-Z -> (X-Z)|(Y-Z)     X Y | Z -  ->  X Z - Y Z - |
		//10. (X|Y)+Z -> (X+Z)|(Y+Z)     X Y | Z +  ->  X Z + Y Z + |
		case 9:
		case 10:
			// Leave X
			// Copy "Z -" or "Z +" after X
			L = Zu-Zl+2;
			// Correct Z position
			_copy(Yl, Zl, L);
			Zl += L;
			Zu += L;
			// Delete the | in front of Z
			expr.erase(Zl-1);
			// Add | at the end
			expr.insert(Zu+1, &GBody::tunion);
			break;
	}

	return rule;
} // _rpnrule

/** --------------------------------------------------------------------- *
 * subTerms
 *
 * This routine returns the pointers in the RPN terms array of the
 * starting point of the left sub-expression LOWLEFT and right
 * sub-expression LOWRIGHT given the pointer of the expr-operator NTX.
 * The searching is performed by scanning from right to left the number
 * of operators and objects pushed into the stack
 *
 * EXP          (expr-left) op (expr-right)
 *
 * RPN ...... | expr-left | expr-right | op | .......
 *           Lowleft      LowRight  op-1 ntx
 * ---------------------------------------------------------------------- */
void CSG::_subTerms(int n, int* lowLeft, int* lowRight)
{
	int nop = 0;
	*lowRight = 0;
	*lowLeft  = 0;

	while (n>=0) {
		if (expr[n]->isOperator())
			nop++;
		else
			nop--;

		n--;
		if (nop==0) {
			if (*lowRight==0) {
				*lowRight = n+1;
				nop++;
				continue;
			} else {
				*lowLeft = n+1;
				return;
			}
		}
	}
} // _subTerms

/** copy a subexpression
 *  @param dst		insert position
 *  @param src		starting position
 *  @param length	length of subexpression
 */
void CSG::_copy(int dst, int src, int length)
{
	// Create space by shifting everything by length
	expr.shiftby(dst, length);

	// Correct starting position if needed
	int nf;
	if (src > dst)
		nf = src + length;
	else
		nf = src;

	// copy terms
	memcpy(&expr[dst], &expr[nf], sizeof(GBody*)*length);
} // _copy

#if 0
//# ----------------------------------------------------------------------
def splitZones(expr):
	"""split a tokenized expression into zones"""
	zones = []
	zone  = []
	depth = 0
	for token in expr:
		if token=="|" and depth==0:
			if zone:
				zones.append(zone)
				zone = []
			continue

		elif token == "(":
			depth += 1

		elif token == ")":
			depth -= 1
			if depth<0:
				raise CSGException("Too many closing parenthesis")

		zone.append(token)

	if zone: zones.append(zone)
	return zones

//# ----------------------------------------------------------------------
def toString(expr):
	"""return a string version of a tokenize expression"""
	s    = ""
	prev = "("
	for token in expr:
		if prev == "(":
			s += token

		elif token == "|":
			s += " | "

		elif token in ("+", "-"):
			s += " %s"%(token)

		else:
			s += token

		prev = token

	return s

//# ----------------------------------------------------------------------
//# Subroutine:	optZone
//# Author:	Vasilis.Vlachoudis@cern.ch
//# Date:		20/4/2004
//#
//# Optimize a product, by removing the universe @, duplicated terms like
//# A+A, A-A, -A-A, and finally calling the geometrical optimizations and
//# sorting by name the primitives in both the plus and minus terms.
//#
//# The product is described by 2 arrays the PLUS,NPLUS and MINUS,NMINUS
//# with all the plus and minus terms of the product
//#
//# WARNING: It doesn't delete the term from the arrays but changes the
//# name to space.
//# ----------------------------------------------------------------------
def optZone(plus,minus):
	//# Remove Universe @ from PLUS
	//# and
	//# remove duplicate terms A+A=A, A-A={}
	for i in xrange(len(plus)):
		if plus[i]=='@': plus[i]=None

		if plus[i]!=None:
			for j in xrange(i+1,len(plus)):
				if plus[i]==plus[j]: plus[j]=None
			for j in xrange(len(minus)):
				if plus[i]==minus[j]:
					//# Remove everything
					del plus[:]
					del minus[:]
					return

	//# Discard product if universe @ exists in MINUS
	//# Check for duplicates in MINUS like -A-A=-A
	for i in xrange(len(minus)):
		if minus[i]=='@':
			//# Remove everything
			del plus[:]
			del minus[:]
			return
		if minus[i]!=None:
			for j in xrange(i+1,len(minus)):
				if minus[i]==minus[j]: minus[j]=None

	//# Perform the Geometrical optimization in the product
	//#call OptGeo(nplus,plus,nminus,minus)

	//# Peform a bubble sort on product terms
	plus.sort()
	minus.sort()

//# ----------------------------------------------------------------------
//# Subroutine:	rmDoubles
//# Author:	Vasilis.Vlachoudis@cern.ch
//# Date:		20/4/2004
//#
//# Remove duplicates of products A+B|A+B|C+D  ->  A+B|C+D
//# ----------------------------------------------------------------------
def rmDoubles(zones):
	i = -1
	while i<len(zones)-1:
		i += 1
		plus1, minus1 = zones[i]
		j = i
		while j<len(zones)-1:
			j += 1
			plus2, minus2 = zones[j]

			//# Check if they are the same
			if len(plus1)!=len(plus2) or len(minus1)!=len(minus2): continue
			diffplus = -1
			for k in range(len(plus1)):
				if plus1[k] != plus2[k]:
					diffplus = k
					break
			//#if not same: continue

			diffminus = -1
			for k in range(len(minus1)):
				if minus1[k] != minus2[k]:
					diffminus = k
					break
			//#if not same: continue

			if diffplus==-1 and diffminus==-1:
				del zones[j]
				j -= 1
				continue

//# ----------------------------------------------------------------------
//# Subroutine:	split
//# Author:	Vasilis.Vlachoudis@cern.ch
//# Date:		18/5/2004
//#
//# Break products into lists
//# ----------------------------------------------------------------------
def split(expr):
	brk=[]
	expr=expr[:]	//# Make a copy of the expression
	while 1:
		try:
			p = expr.index("|")
			brk.append(expr[0:p])
			del expr[0:p+1]
		except:
			brk.append(expr)
			return brk

//# ----------
if __name__=="__main__":
	import sys, string, pprint
	expr = tokenize(string.join(sys.argv[1:]))

	print "Zones="
	pprint.pprint(splitZones(expr))

	print "ZoneString="
	pprint.pprint(map(toString, splitZones(expr)))

	print "ORG=",expr

	exp2rpn(expr)
	print "RPN=",expr

	rpnorm(expr)
	print "RPNORM=",expr

	expnorm = rpn2exp(expr)
	print "NORM=",expnorm

//#	list = breakprod(expnorm)
	expr = string.join(expnorm)
	expr = expr.replace("+ ","+")
	expr = expr.replace("- ","-")
	print "EXP=",expr
#endif

/** operator << */
ostream& operator << (ostream& s, const CSG& csg)
{
	for (int i=0; i<csg.expr.count(); i++)
		s << csg.expr[i]->name() << ' ';
	return s;
} /* operator << */
