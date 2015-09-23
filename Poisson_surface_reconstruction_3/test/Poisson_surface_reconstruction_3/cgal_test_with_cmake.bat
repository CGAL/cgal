@echo off

rem Double-clickable version of cgal_test_with_cmake on Windows/Cygwin/VisualC++

rem Path to Cygwin is hard-coded
PATH=c:\Applis\cygwin\bin;%PATH%

rem Add path to VisualC++
if "%CMAKE_GENERATOR%"=="-GVisual Studio 8 2005"       PATH=%VS80COMNTOOLS%\..\IDE;%PATH%
if "%CMAKE_GENERATOR%"=="-GVisual Studio 8 2005 Win64" PATH=%VS80COMNTOOLS%\..\IDE;%PATH%
if "%CMAKE_GENERATOR%"=="-GVisual Studio 9 2008"       PATH=%VS90COMNTOOLS%\..\IDE;%PATH%
if "%CMAKE_GENERATOR%"=="-GVisual Studio 9 2008 Win64" PATH=%VS90COMNTOOLS%\..\IDE;%PATH%

rem Call shell script
sh ./cgal_test_with_cmake

pause

rem Compare logs with Beyond Compare
BC2.exe ./cgal_test_with_cmake.log.bak ./cgal_test_with_cmake.log
