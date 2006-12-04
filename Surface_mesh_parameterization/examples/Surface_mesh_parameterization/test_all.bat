@echo off

rem Double-clickable version of test_all.sh on Windows/Cygwin
rem sh.exe must be in %PATH%

sh ./test_all.sh 2>&1 | tee test/test_all.log
pause
