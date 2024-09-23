# $Id$
# $Log: check.py,v $
# Revision 1.3  2010/02/17 08:48:09  bnv
# PreRelease 0.8.0
#
# Revision 1.2  2007/06/01 07:25:59  bnv
# pipe: changed to popen to avoid blocking the stdout
#
# Revision 1.1  2006/09/29 15:16:35  bnv
# *** empty log message ***
#
#
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	17-May-2004

import pdb
import bmath

def __checkNumber(n,l=10,e=False):
	m = bmath.format(n,l,e)
	print((n, m, len(m), l, e))

def __checkFormating():
	__checkNumber("9.99999999999999999999999998",10)

	__checkNumber("0",10)
	__checkNumber("000000000",10)
	__checkNumber("0.0",10)
	__checkNumber("0.00000000",10)
	__checkNumber("0.00000000e0000",10)
	__checkNumber("0.00000000e1000",10)
	__checkNumber("0.00000000e-100",10)
	__checkNumber("-000000000",10)
	__checkNumber("-0.00000000",10)
	__checkNumber("-0.00000000e0000",10)
	__checkNumber("-0.00000000e1000",10)
	__checkNumber("-0.00000000e-100",10)
	print()

	__checkNumber("1",10)
	__checkNumber("10",10)
	__checkNumber("100",10)
	__checkNumber("1000",10)
	__checkNumber("1000000",10)
	__checkNumber("1000000000",10)
	__checkNumber("1000000000000",10)
	print()

	__checkNumber("0.1",10)
	__checkNumber("0.01",10)
	__checkNumber("0.001",10)
	__checkNumber("0.0001",10)
	__checkNumber("0.000001",10)
	__checkNumber("0.000000001",10)
	__checkNumber("0.000000000001",10)
	print()

	__checkNumber("0.100000",10)
	__checkNumber("0.0100000",10)
	__checkNumber("0.00100000",10)
	__checkNumber("0.000100000",10)
	__checkNumber("0.00000100000",10)
	__checkNumber("0.00000000100000",10)
	__checkNumber("0.00000000000100000",10)
	print()

	__checkNumber("0.100006",10)
	__checkNumber("0.01000006",10)
	print()

	__checkNumber("-1",10)
	__checkNumber("-10",10)
	__checkNumber("-100",10)
	__checkNumber("-1000",10)
	__checkNumber("-1000000",10)
	__checkNumber("-1000000000",10)
	__checkNumber("-1000000000000",10)
	print()

	__checkNumber("-0.1",10)
	__checkNumber("-0.01",10)
	__checkNumber("-0.001",10)
	__checkNumber("-0.0001",10)
	__checkNumber("-0.000001",10)
	__checkNumber("-0.000000001",10)
	__checkNumber("-0.000000000001",10)
	print()

	__checkNumber("123.456789012345678e13",10)
	__checkNumber("-123.456789012345678e-13",10)
	__checkNumber("+123.456789012345678e-13",10)
	__checkNumber("123456789012345678",10)
	__checkNumber("-123456789012345678",10)
	__checkNumber("+123456789012345678",10)
	__checkNumber("12345678901234567800000000000000",10)
	__checkNumber("-0.00000000000123456789012345678",10)
	__checkNumber("+0.00000000000123456789012345678",10)
	print()

	__checkNumber("000000000",10,True)
	__checkNumber("0.00000000",10,True)
	__checkNumber("0.00000000e0000",10,True)
	__checkNumber("0.00000000e1000",10,True)
	__checkNumber("0.00000000e-100",10,True)
	__checkNumber("-000000000",10,True)
	__checkNumber("-0.00000000",10,True)
	__checkNumber("-0.00000000e0000",10,True)
	__checkNumber("-0.00000000e1000",10,True)
	__checkNumber("-0.00000000e-100",10,True)
	print()

	__checkNumber("+1",10)
	__checkNumber("+10",10)
	__checkNumber("+100",10)
	__checkNumber("+1000",10)
	__checkNumber("+1000000",10)
	__checkNumber("+1000000000",10)
	__checkNumber("+1000000000000",10)
	print()

	__checkNumber("+0.1",10)
	__checkNumber("+0.01",10)
	__checkNumber("+0.001",10)
	__checkNumber("+0.0001",10)
	__checkNumber("+0.000001",10)
	__checkNumber("+0.000000001",10)
	__checkNumber("+0.000000000001",10)
	print()

	__checkNumber("000000000",10,True)
	__checkNumber("0.00000000",10,True)
	__checkNumber("0.00000000e0000",10,True)
	__checkNumber("0.00000000e1000",10,True)
	__checkNumber("0.00000000e-100",10,True)
	__checkNumber("-000000000",10,True)
	__checkNumber("-0.00000000",10,True)
	__checkNumber("-0.00000000e0000",10,True)
	__checkNumber("-0.00000000e1000",10,True)
	__checkNumber("-0.00000000e-100",10,True)
	print()

	__checkNumber("1",10,True)
	__checkNumber("10",10,True)
	__checkNumber("100",10,True)
	__checkNumber("1000",10,True)
	__checkNumber("1000000",10,True)
	__checkNumber("1000000000",10,True)
	__checkNumber("1000000000000",10,True)
	print()

	__checkNumber("0.1",10,True)
	__checkNumber("0.01",10,True)
	__checkNumber("0.001",10,True)
	__checkNumber("0.0001",10,True)
	__checkNumber("0.000001",10,True)
	__checkNumber("0.000000001",10,True)
	__checkNumber("0.000000000001",10,True)
	print()

	__checkNumber("0.100000",10,True)
	__checkNumber("0.0100000",10,True)
	__checkNumber("0.00100000",10,True)
	__checkNumber("0.000100000",10,True)
	__checkNumber("0.00000100000",10,True)
	__checkNumber("0.00000000100000",10,True)
	__checkNumber("0.00000000000100000",10,True)
	print()

	__checkNumber("0.100006",10,True)
	__checkNumber("0.01000006",10,True)
	print()

	__checkNumber("-1",10,True)
	__checkNumber("-10",10,True)
	__checkNumber("-100",10,True)
	__checkNumber("-1000",10,True)
	__checkNumber("-1000000",10,True)
	__checkNumber("-1000000000",10,True)
	__checkNumber("-1000000000000",10,True)
	print()

	__checkNumber("-0.1",10,True)
	__checkNumber("-0.01",10,True)
	__checkNumber("-0.001",10,True)
	__checkNumber("-0.0001",10,True)
	__checkNumber("-0.000001",10,True)
	__checkNumber("-0.000000001",10,True)
	__checkNumber("-0.000000000001",10,True)
	print()

	__checkNumber("123.4567890123456781238971239e13",10,True)
	__checkNumber("-123.4567890123456712389712938e-13",10,True)
	__checkNumber("12345678901234567812989127389123",10,True)
	__checkNumber("-12345678901234567812389712389",10,True)
	__checkNumber("123450000008901234567812989127389123",10,True)
	__checkNumber("-123450000008901234567812389712389",10,True)
	__checkNumber("-0.000000000000008901234567812389712389",10,True)
	__checkNumber("12345678901234567800000000000000",10,True)
	__checkNumber("-0.00000000000123456789012345678",10,True)
	print()

if __name__ == "__main__":
	__checkFormating()
	#inp = Input()
	#inp.verbose = True
	#inp.read(sys.argv[1])
	#inp.write("test.inp")
