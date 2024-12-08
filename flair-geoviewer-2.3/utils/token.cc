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

#include <fstream>

#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "os.h"
#include "token.h"

using namespace std;

const static char* tokenType[] = {
			"EOF",
			";",
			"ident",
			"string",
			"file",
			"function",
			"number",
			",",
			":",
			"=",
			"+",
			"-",
			"*",
			"/",
			"^",
			"|",
			"&",
			"{",
			"}",
			"[",
			"]",
			"(",
			")",
			"$"
		};

/** init */
void Token::init(bool inl, bool ic)
{
	newline      = true;
	linenum      = 0;
	colpos       = 0;
	cur_token    = string_tok;
	num          = 0.0;
	str[0]       = 0;
	ignore_newline = inl;
	ignore_case  = ic;
	token_pushed = false;
	memset(extended_alpha, 0, sizeof(extended_alpha));
} // init

/** push token if not needed */
void Token::push()
{
	token_pushed = true;
} // push

/** next
 * Gets the next token from the input stream and
 * returns it's value.  If the previous token has been
 * pushed back, then it is returned again.
 * @return next token
 */
TokenType Token::next()
{
	int	i,c;

	if (token_pushed) {
		token_pushed = false;
		return cur_token;
	}

	i = 0;
	str[0] = 0;

	while (newline) {
		newline = false;
		c = stream.get();	colpos++;
		if (c=='#') {
			stream >> linenum;
			//stream >> currentfile;
			while ((c=stream.get())!='\n') colpos++;
			newline = true;
			colpos=0;
		} else {
			stream.putback(c);
			colpos--;
		}
	}

	/* ignore white space */
	while (isspace(c=stream.get())) {
		colpos++;
		if (c == '\n') {
			linenum++;
			newline = true;
			colpos = 0;
			if (!ignore_newline) {
				cur_token = semicolon_tok;
				return cur_token;
			}
		}
	}
	colpos++;

	/* check for end of file */
	str[0] = c;	// store in string
	str[1] = 0;
	if (isalpha(c) || extended_alpha[(uint8_t)c]) {         /* must be a keyword */
		do {
			str[i++] = c;
			c = stream.get();
			colpos++;
		} while (isalnum(c) || c=='_');
		str[i] = 0;
		if (c=='(')
			cur_token = function_tok;
		else {
			if (ignore_case) upper();
			cur_token = ident_tok;
			// push back the character that doesn't belong
			stream.putback(c);
			colpos--;
		}
		return cur_token;
	} else
	if (isdigit(c) || c=='+' || c=='.' || c=='-') {
		do {
			str[i++] = c;
			c = stream.get();
			colpos++;
		} while (isdigit(c) || c=='+' || c=='.'
				|| c=='-' || c=='e' || c=='E');
		stream.putback(c);	// push back the character that doesn't belong
		colpos--;

		str[i] = 0;
		num = atof(str);
		cur_token = number_tok;
		return cur_token;
	} else
	switch (c) {
		case EOF:
			cur_token = eof_tok;
			return cur_token;

		case '{':
			cur_token = left_brace_tok;
			return cur_token;

		case '}':
			cur_token = right_brace_tok;
			return cur_token;

		case '[':
			cur_token = left_bracket_tok;
			return cur_token;

		case ']':
			cur_token = right_bracket_tok;
			return cur_token;

		case '(':
			cur_token = left_paren_tok;
			return cur_token;

		case ')':
			cur_token = right_paren_tok;
			return cur_token;

		case ';':
			cur_token = semicolon_tok;
			return cur_token;

		case ':':
			cur_token = colon_tok;
			return cur_token;

		case '+':
			cur_token = plus_tok;
			return cur_token;

		case '-':
			if (c=='-')
				cur_token = minus_tok;
			else
				cur_token = plus_tok;
			// peek next to see if it is a number
			c = stream.get();
			stream.putback(c);
			if (!isdigit(c))
				return cur_token;
			if (cur_token == minus_tok)
				c = '-';
			else
				c = '+';
			// try a number, go below
			break;

		case '*':
			cur_token = mult_tok;
			return cur_token;

		case '/':
			cur_token = div_tok;
			return cur_token;

		case '^':
			cur_token = not_tok;
			return cur_token;

		case '|':
			cur_token = or_tok;
			return cur_token;

		case '=':
			cur_token = equal_tok;
			return cur_token;

		case ',':
			cur_token = comma_tok;
			return cur_token;

		case '$':
			cur_token = comma_tok;
			return dollar_tok;

		case '\"':
			while((c = stream.get()) != '\"') {
				str[i++] = c;
				colpos++;
			}
			str[i] = 0;
			cur_token = string_tok;
			return cur_token;

		case '<':
			while((c = stream.get()) != '>') {
				str[i++] = c;
				colpos++;
			}
			str[i] = 0;
			cur_token = file_string_tok;
			return cur_token;
	}

	// if we get down here something is really wrong
	cerr << "\nERROR: Found the character: " << (char)c << endl;
//		<< " Hex 0x" << hex << c << dec << endl;
	Token::error(ERR_INVALID_CHAR);
	return cur_token;	// never gets here but keeps compiler happy
} // next

/** printErrorLine
 * if an error is found print current line with a pointer to the wrong token
 */
void Token::printErrorLine(std::ostream& os)
{
	// Move file pointer "colpos" back and display current line
	stream.seekg(-(long)colpos-1,ios::cur);
	os.width(5);
	os << linenum << " *-* ";
	os.width();
	int col = 0;
	int pos = 0;
	int colerr = 0;
	while (1) {
		int c = stream.get(); pos++;
		if (pos==colpos) colerr = col;
		if (c=='\n' || c==EOF) break;
		if (c=='\t') {
			c = 8-col%8;
			for (int i=0; i<c; i++) os << ' ';
			col += c;
		} else {
			os << (char)c;
			col++;
		}
	}
	colerr += 10;
	os << endl;

	for (int i=0; i<colerr; i++) os << ' ';
	os << '^' << endl;
} // printErrorLine

/** error
 * @param errnum error number to be printed
 */
void Token::error(const int errnum)
{
	printErrorLine(cerr);
	cerr << "Error " << errnum << " parsing line " << linenum << endl;
//	assert(0);
//	exit(errnum);
} // error

/** skipSemicolon */
void Token::skipSemicolon()
{
	while (cur_token==semicolon_tok) next();
} // skipSemicolon

/** skipLine */
void Token::skipLine()
{
	int c;
	do {
		c = stream.get();
	} while (c!='\n' && c!=-1);
	next();
} // skipLine

/** nextSkipEqual */
void Token::nextSkipEqual()
{
	next();
	if (cur_token == colon_tok || cur_token==equal_tok)
		next();
} // nextSkipEqual

/** mustbe
 * @param tok next token must be tok otherwise display error
 * @param errnum error number to display
 */
bool Token::mustbe(const TokenType tok, const int errnum)
{
	if (cur_token != tok)
		Token::error(errnum);
	else {
		next();
		return true;
	}
	return false;
} // mustbe

/** mustbe
 * @param tok next token must be tok otherwise display error
 * @param errnum error number to display
 */
bool Token::mustbe(const char* tokstr, const int errnum)
{
	if (cur_token != ident_tok)
		Token::error(errnum);
	else
	if (strcmp(ident(), tokstr))
		Token::error(errnum);
	else {
		next();
		return true;
	}
	return false;
} // mustbe

/** mustbe
 * @param tok next token must be tok otherwise display error
 * @param errnum error number to display
 */
bool Token::mustbe(const char tokchar, const int errnum)
{
	if (character() != tokchar)
		Token::error(errnum);
	else {
		next();
		return true;
	}
	return false;
} // mustbe

/**  upper - convert to uppercase
 */
void Token::upper()
{
	for(char *c=str; *c; c++)
		*c = toupper(*c);
} // upper

/** cmp
 * Perform an ABBREVIATED comparison
 * a star '*' marks the start of abbreviation
 * @param teststr	abbreviation to check token
 * @return true/false if it matches the abbreviation
 */
bool Token::cmp(const char *teststr)
{
	char	*c = str;

	// until the star must be exactly equal
	while (*c && *teststr && *c == *teststr) {
		c++;
		teststr++;
	}
	// skip the star
	if (*teststr=='*')
		teststr++;
	else
	if (*teststr || *c)
		return false;
	else
		return true;

	// compare abbreviated
	while (*c && *teststr && *c == *teststr) {
		c++;
		teststr++;
	}
	if (*c == 0)
		return true;
	else
		return false;
} // cmp

/** getString, get string from file
 * @param getstr	string to return
 * @param maxsize	maximum size of getstr
 */
void Token::getString(char *getstr, const int maxsize)
{
	if (cur_token == string_tok || cur_token == ident_tok) {
		strncpy(getstr, str, maxsize);
		next();
	} else
		Token::error(ERR_STRING_EXPECTED);
} // getString

/** getInteger
 * @return next integer token
 */
int Token::getInteger()
{
	if (cur_token == number_tok) {
		int n = (int)num;
		if (!Eq((double)n, num, 1.0e-15))
			Token::error(ERR_INTEGER_EXPECTED);
		next();
		return n;
	} else
		Token::error(ERR_INTEGER_EXPECTED);
	return 0;	// Never gets here but keeps compiler happy
} // getInteger

/** getBoolean
 * @return next boolean token
 */
bool Token::getBoolean()
{
	if (cur_token == number_tok) {
		int	n = (int)num;
		if ((!Eq((double)n, num, 1.0e-15)) || (n!=0 && n!=1))
			Token::error(ERR_INVALID_BOOLEAN);
		next();
		return n;
	} else
	if (cur_token == string_tok || cur_token == ident_tok) {
		upper();
		if (cmp("T*RUE") || cmp("ON") || cmp("YES")) {
			next();
			return true;
		}
		if (cmp("F*ALSE") || cmp("OFF") || cmp("NO")) {
			next();
			return false;
		}
	}
	Token::error(ERR_INVALID_BOOLEAN);
	return 0;	// Never gets here but keeps compiler happy
} // getBoolean

/** getNumber
 * @return next floating point number
 */
double Token::getNumber()
{
	if (cur_token == number_tok) {
		double	n = num;
		next();
		return n;
	} else
		Token::error(ERR_NUMBER_EXPECTED);
	return 0.0;	// Never gets here but keeps compiler happy
} // getNumber

/** getUpperIdent
 * @param errnum	error number to display in case of problem
 * @return next identifier as uppercase
 */
bool Token::getUpperIdent(const int errnum)
{
	skipSemicolon();
	if (current() == right_brace_tok) {
		next();
		return false;
	}

	if (current() != ident_tok)
		Token::error(errnum);

	upper();
	return true;
} // getUpperIdent

/** print token */
std::ostream& operator << (std::ostream& s, const Token& token)
{
	s << "Token ";
	s << " Type=" << tokenType[token()];
	s << " Value=\"" << token.string() << '"';
	return s;
} /* operator << */
