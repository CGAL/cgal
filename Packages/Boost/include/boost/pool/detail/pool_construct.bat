@echo off
rem
rem Copyright (C) 2000, 2001 Stephen Cleary
rem
rem This file can be redistributed and/or modified under the terms found
rem  in "copyright.html"
rem This software and its documentation is provided "as is" without express or
rem  implied warranty, and with no claim as to its suitability for any purpose.
rem
rem See http://www.boost.org for updates, documentation, and revision history.
rem

rem Check for Windows NT
if %OS%==Windows_NT goto NT

rem Not NT - run m4 as normal, then exit
m4 -P -E -DNumberOfArguments=%1 pool_construct.m4 > pool_construct.inc
goto end

rem DJGPP programs (including m4) running on Windows/NT do NOT support long
rem  file names (see the DJGPP v2 FAQ, question 8.1)
rem Note that the output doesn't have to be a short name because it's an
rem  argument to the command shell, not m4.
:NT
m4 -P -E -DNumberOfArguments=%1 < pool_construct.m4 > pool_construct.inc

:end