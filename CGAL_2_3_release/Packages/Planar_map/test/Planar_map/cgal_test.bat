@echo off
set cfiles=test01 test02 test03 test04
set adv_cfiles=test05
set dfiles=1 2 3 4 5 6
set cext=.C
set dpre=DATA\
set din=.ix.
set rpre=RESULT\
set rin=.ix.
set testers=1 2
set strategies=1 2 3
set bboxes=1
set adv_bboxes=2
set make=nmake -nologo -f makefile.ms
set cmp=fc
set OUTPUTFILE=error.txt

rem help page
if "%1"=="-h" goto help
if "%1"=="-help" goto help
if "%1"=="?" goto help
if "%1"=="/?" goto help

rem execute and log on report
if not "%5"=="" goto execute

rem fix bounding box strategy
if not "%4"=="" goto bbox

rem fix stategy
if not "%3"=="" goto strategy

rem fix tester
if not "%2"=="" goto tester

rem fix file
if not "%1"=="" goto file

:test
echo. > %OUTPUTFILE%
rem set EXTRA_FLAGS_ORG="%EXTRA_FLAGS%"

set advanced=false
for %%f in (%cfiles%) do if exist %%f%cext% call cgal_test %%f
set advanced=true
for %%f in (%adv_cfiles%) do if exist %%f%cext% call cgal_test %%f
goto end

:file
for %%t in (%testers%) do if exist %dpre%%1%din%* call cgal_test %1 %%t
goto end

:tester
for %%s in (%strategies%) do call cgal_test %1 %2 %%s
goto end

:strategy
if "%advanced%"=="false" for %%b in (%bboxes%) do cgal_test %1 %2 %3 %%b
if "%advanced%"=="true" for %%b in (%adv_bboxes%) do cgal_test %1 %2 %3 %%b
goto end

:bbox
rem set EXTRA_FLAGS=%EXTRA_FLAGS% -DTESTR#%2 -DSTRATEGY#%3 -DBBOX#%4

:make
%make% clean %1.exe EXTRA_FLAGS="%EXTRA_FLAGS% -DTESTR#%2 -DSTRATEGY#%3 -DBBOX#%4"
if errorlevel==1 goto makefail

:makesucc
echo "compilation of %1 with TESTR=%2 STRATEGY=%3 BBOX=%4 succeeded" >> %OUTPUTFILE%

:ready
for %%d in (%dfiles%) do if exist %dpre%%1%din%%%d call cgal_test %1 %2 %3 %4 %%d
goto end

:makefail
echo "ERROR: compilation of %1 with TESTR=%2 STRATEGY=%3 BBOX=%4 and %1 failed" >> %OUTPUTFILE%
goto end

:execute
set goodresult="%rpre%%1%rin%%5.%2.%4"
%1 < %dpre%%1%din%%5 > aresult
if not exist %goodresult% goto nottested
%cmp% aresult %goodresult% >> %OUTPUTFILE%
goto succeed

:nottested
echo "ERROR: %1 with TESTR=%2 STRATEGY=%3 BBOX=%4 was not tested on %goodresult%" >> %OUTPUTFILE%
copy aresult %goodresult%
goto end

:succeed
echo  "execution of %1 TESTR=%2 STRATEGY=%3 BBOX=%4 with %goodresult% succesful">> %OUTPUTFILE%
goto done

:failed
echo "ERROR: %1 < %dpre%%1%din%%2%3 did not run well" >> %OUTPUTFILE%
goto done

:done
rem if exist %1.exe del %1.exe
rem if exist %1.obj del %1.obj
goto end

:help
echo Syntax:
echo    cgal_test [testfile]
echo.
echo Abstract:
echo    The utility tests your package again each *%cext% file.
echo    make is called for each of the different TESTR,STRATEGY and
echo    BBOX values and then to executed with the data files in
echo    %dpre%, comparing the results to %rpre%.
echo.
echo    Passing one argument to cgal_test, e.g. cgal_test tst20 will 
echo    selectively test only tst20%cext%.
goto end

:end
rem set EXTRA_FLAGS="%EXTRA_FLAGS_ORG%"
rem set EXTRA_FLAGS_ORG=
set make=
set cmp=
set cpre=
set cfiles=
set cext=
set dpre=
set din=
set rpre=
set rin=
set strategies=
set bboxes=
set testers=
set goodresult=
set advaced=

