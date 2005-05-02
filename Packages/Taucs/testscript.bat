REM This is an automatically-generated script
del testscript.log
echo TAUCS TEST LOG >  testscript.log
echo ============== >> testscript.log
echo Win32          >> testscript.log
echo ============== >> testscript.log
echo =============== >> testscript.log
echo = test_cilk_snmf = >> testscript.log
call configure in=progs\test_cilk_snmf.c %*
nmake /F build\%TAUCS_LASTCONF%\makefile 
bin\%TAUCS_LASTCONF%\test_cilk_snmf >> testscript.log
if errorlevel 1 goto :error_test_cilk_snmf
echo = TEST PASSED (test_cilk_snmf) >> testscript.log
goto :next_test_cilk_snmf
:error_test_cilk_snmf
echo = TEST FAILED (test_cilk_snmf) >> testscript.log
:next_test_cilk_snmf
echo =============== >> testscript.log
echo =============== >> testscript.log
echo = test_complex = >> testscript.log
call configure in=progs\test_complex.c %*
nmake /F build\%TAUCS_LASTCONF%\makefile 
bin\%TAUCS_LASTCONF%\test_complex >> testscript.log
if errorlevel 1 goto :error_test_complex
echo = TEST PASSED (test_complex) >> testscript.log
goto :next_test_complex
:error_test_complex
echo = TEST FAILED (test_complex) >> testscript.log
:next_test_complex
echo =============== >> testscript.log
echo =============== >> testscript.log
echo = test_linsolve = >> testscript.log
call configure in=progs\test_linsolve.c %*
nmake /F build\%TAUCS_LASTCONF%\makefile 
bin\%TAUCS_LASTCONF%\test_linsolve >> testscript.log
if errorlevel 1 goto :error_test_linsolve
echo = TEST PASSED (test_linsolve) >> testscript.log
goto :next_test_linsolve
:error_test_linsolve
echo = TEST FAILED (test_linsolve) >> testscript.log
:next_test_linsolve
echo =============== >> testscript.log
echo =============== >> testscript.log
echo = test_memory = >> testscript.log
call configure in=progs\test_memory.c %*
nmake /F build\%TAUCS_LASTCONF%\makefile 
bin\%TAUCS_LASTCONF%\test_memory >> testscript.log
if errorlevel 1 goto :error_test_memory
echo = TEST PASSED (test_memory) >> testscript.log
goto :next_test_memory
:error_test_memory
echo = TEST FAILED (test_memory) >> testscript.log
:next_test_memory
echo =============== >> testscript.log
echo =============== >> testscript.log
echo = test_notimer = >> testscript.log
call configure in=progs\test_notimer.c %*
nmake /F build\%TAUCS_LASTCONF%\makefile 
bin\%TAUCS_LASTCONF%\test_notimer >> testscript.log
if errorlevel 1 goto :error_test_notimer
echo = TEST PASSED (test_notimer) >> testscript.log
goto :next_test_notimer
:error_test_notimer
echo = TEST FAILED (test_notimer) >> testscript.log
:next_test_notimer
echo =============== >> testscript.log
echo =============== >> testscript.log
echo = test_stack = >> testscript.log
call configure in=progs\test_stack.c %*
nmake /F build\%TAUCS_LASTCONF%\makefile 
bin\%TAUCS_LASTCONF%\test_stack >> testscript.log
if errorlevel 1 goto :error_test_stack
echo = TEST PASSED (test_stack) >> testscript.log
goto :next_test_stack
:error_test_stack
echo = TEST FAILED (test_stack) >> testscript.log
:next_test_stack
echo =============== >> testscript.log
echo =============== >> testscript.log
echo = test_timer = >> testscript.log
call configure in=progs\test_timer.c %*
nmake /F build\%TAUCS_LASTCONF%\makefile 
bin\%TAUCS_LASTCONF%\test_timer >> testscript.log
if errorlevel 1 goto :error_test_timer
echo = TEST PASSED (test_timer) >> testscript.log
goto :next_test_timer
:error_test_timer
echo = TEST FAILED (test_timer) >> testscript.log
:next_test_timer
echo =============== >> testscript.log
