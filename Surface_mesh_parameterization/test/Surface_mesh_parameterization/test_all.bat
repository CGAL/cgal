@echo off

rem Double-clickable version of test_all.sh on Windows/Cygwin

rem Path to Cygwin is hard-coded
PATH=%PATH%;c:\Applis\cygwin\bin

sh ./test_all.sh 2>&1 | tee test_all.log
pause
