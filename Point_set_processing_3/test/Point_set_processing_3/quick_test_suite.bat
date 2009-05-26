@echo off

rem Double-clickable version of quick_test_suite.sh on Windows/Cygwin

rem Path to Cygwin is hard-coded
PATH=%PATH%;c:\Applis\cygwin\bin

sh ./quick_test_suite.sh
pause

rem Compare logs with Beyond Compare
"C:\Program Files\Beyond Compare 2\BC2.exe" ./quick_test_suite.log.bak ./quick_test_suite.log
