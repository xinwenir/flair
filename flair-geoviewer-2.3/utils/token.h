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

#ifndef __TOKEN_H
#define __TOKEN_H

#include <ostream>
#include <iostream>
#include <cstdint>
#include "os.h"

enum TokenType {
	eof_tok,		//  0
	semicolon_tok,		//  1 ;
	ident_tok,		//  2 identifier
	string_tok,		//  3 "string"
	file_string_tok,	//  4 <string> file string, to be searched
	function_tok,		//  5 ident(
	number_tok,		//  6 0.12345789E-5
	comma_tok,		//  7 ,
	colon_tok,		//  8 :
	equal_tok,		//  9 =
	plus_tok,		// 10 +
	minus_tok,		// 11 -
	mult_tok,		// 12 *
	div_tok,		// 13 /
	not_tok,		// 14 ^
	or_tok,			// 15 |
	and_tok,		// 16 &
	left_brace_tok,		// 17 {
	right_brace_tok,	// 18 }
	left_bracket_tok,	// 19 [
	right_bracket_tok,	// 20 ]
	left_paren_tok,		// 21 (
	right_paren_tok,	// 22 )
	dollar_tok		// 23 $
};

#define ERR_OK	0
#define ERR_FILE_NOT_FOUND	1
#define ERR_INVALID_CHAR	2
#define ERR_STRING_EXPECTED	3
#define ERR_NUMBER_EXPECTED	4
#define ERR_INTEGER_EXPECTED	5
#define ERR_EXTRA_DATA		6
#define ERR_INVALID_SYNTAX	7
#define ERR_OBJECT_EXISTS	10
#define ERR_INVALID_BOOLEAN	12
#define	ERR_OBJECT_NOTFOUND	16
#define ERR_RIGHT_PAREN_EXPECTED 18

/*-------------------------------------------- */
/*                Token                       */
/*-------------------------------------------- */
class Token {
protected:
	std::istream&	stream;
	bool		newline;
	int		linenum;
	int		colpos;
	TokenType	cur_token;
	bool		token_pushed;

	double		num;
	long		integer;
	char		str[128];
	bool		extended_alpha[256];

public:
	bool		ignore_case;
	bool		ignore_newline;

public:
	Token(std::istream& s, bool inl=false, bool ic=false)
		: stream(s) { init(inl,ic); }
	~Token() {}

	void	printErrorLine(std::ostream& os);
virtual	void	error(const int);
	void	setAlpha(char c)	{ extended_alpha[(uint8_t)c] = true; }

	void	push();
	TokenType next();
	TokenType operator ()() const	{ return cur_token; }
	TokenType current() const	{ return cur_token; }

	double	number() const		{ return num; }
const	char*	string() const		{ return str; }
const	char*	ident()	 const		{ return str; }
	char	character() const	{ return str[0]; }
	void	upper();		// Translate str to uppercase
	bool	cmp(const char *);	// abbreviated comparison

	void	skipSemicolon();
	void	skipLine();
	void	nextSkipEqual();

	bool	mustbe(const TokenType, const int);
	bool	mustbe(const char*, const int);
	bool	mustbe(const char, const int);
	void	getString(char *, const int);
	int	getInteger();
	bool	getBoolean();
	double	getNumber();
	bool	getUpperIdent(const int errnum);

private:
	void	init(bool inl=false, bool ic=false);
}; // Token

std::ostream& operator << (std::ostream& s, const Token& token);

#define MUSTBE_SEMICOLON	token.mustbe(semicolon_tok,ERR_EXTRA_DATA)
#endif
