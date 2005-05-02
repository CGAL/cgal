@echo off

set OSTYPE=win32

nmake /Fconfigurator\makefile.win32

configurator\configurator win32 %* > .lastconf

for /F "usebackq" %%i IN (`type .lastconf`) DO set TAUCS_LASTCONF=%%i

rem LS 04/2005: was "del .lastconf /a:h"
del .lastconf

goto :EOF

