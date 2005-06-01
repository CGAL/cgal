@echo off

rem Double-clickable version of test_all.sh on Windows/Cygwin

sh .\test_all.sh | tee test_all.log
pause
