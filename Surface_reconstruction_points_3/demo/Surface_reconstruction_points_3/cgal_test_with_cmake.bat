@echo off

rem Double-clickable version of cgal_test_with_cmake on Windows/Cygwin/VisualC++

rem Path to Cygwin and VisualC++ is hard-coded
PATH=c:\Applis\cygwin\bin;C:\Program Files\Microsoft Visual Studio 8\Common7\IDE;%PATH%

sh ./cgal_test_with_cmake
pause

rem Compare logs with Beyond Compare
"C:\Program Files\Beyond Compare 2\BC2.exe" ./cgal_test_with_cmake.log.bak ./cgal_test_with_cmake.log
