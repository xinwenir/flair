@ECHO OFF
REM $Id: flair.bat,v 1.5 2010/02/23 10:34:20 bnv Exp $
REM shell script to launch the program
REM Author:	Vasilis.Vlachoudis@cern.ch
REM Date:	13-Jun-2006
set PYTHONPATH=c:\Documents And Settings\bnv\prg\physics\fluka\flair\lib
c:\soft\python24\python flair.py %1
