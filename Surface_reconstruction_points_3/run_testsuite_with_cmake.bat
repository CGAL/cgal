@echo off

rem Double-clickable version of run_testsuite_with_cmake on Windows/Cygwin

rem Path to Cygwin is hard-coded
PATH=c:\Applis\cygwin\bin;C:\Program Files\Microsoft Visual Studio 8\Common7\IDE;Z:\src\cgal\SVNROOT\trunk\Testsuite\test;%PATH%

sh -c "run_testsuite_with_cmake ; cat error.txt"
pause

