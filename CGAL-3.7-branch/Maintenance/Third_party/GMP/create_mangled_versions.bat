@echo OFF

set ROOT=%~dp0

if "%1"=="" goto :ERROR

set BASEDIR=%1%

cd %BASEDIR%

if EXIST gmp.h (
  set GMP_OR_MPFR=gmp
) else (
  if EXIST mpfr.h set GMP_OR_MPFR=mpfr
)

if NOT DEFINED GMP_OR_MPFR goto :ERROR

if NOT EXIST mangled_binaries mkdir mangled_binaries

if EXIST build.vc8 call :PROCESS_BUILD_FOLDER vc8
if EXIST build.vc9 call :PROCESS_BUILD_FOLDER vc9

cd ..

goto :EOF

:PROCESS_BUILD_FOLDER
: %1 <- vc(8|9)

cd build.%1

for /D %%i in (lib_%GMP_OR_MPFR%*) do (

  call :PROCESS_LIB_FOLDER %%i %1
)
   
cd ..

goto :EOF

:PROCESS_LIB_FOLDER
: %1 <- lib_gmp|mpfr_* folder name
: %2 <- vc(8|9)

cd %1
if EXIST Win32 call :PROCESS_PALTFORM_FOLDER Win32 %2
if EXIST x64   call :PROCESS_PALTFORM_FOLDER x64   %2
cd ..

goto :EOF

:PROCESS_PALTFORM_FOLDER
: %1 <- (Win32|x64)
: %2 <- vc(8|9)

cd %1
if EXIST Debug   call :PROCESS_VARIANT_FOLDER Debug   %2 -gd
if EXIST Release call :PROCESS_VARIANT_FOLDER Release %2
cd ..

goto :EOF

:PROCESS_VARIANT_FOLDER
: %1 <- (Debug|Release)
: %2 <- vc(8|9)
: %3 <- (-gd|)

cd %1
if EXIST %GMP_OR_MPFR%.lib call :COPY_MANGLED %2 lib %3
if EXIST %GMP_OR_MPFR%.pdb call :COPY_MANGLED %2 pdb %3
cd ..

goto :EOF


:COPY_MANGLED
: %1 <- vc(8|9)
: %2 <- lib|pdb
: %3 <- (-gd|)
  
set TOOLSET=%1%
set EXT=%2%
set VARIANT=%3%

set MANGLED=%GMP_OR_MPFR%-%TOOLSET%-mt%VARIANT%.%EXT%
    
echo Copying %GMP_OR_MPFR%.%EXT% as %MANGLED% in %cd% and %ROOT%\mangled_binaries
    
copy /Y %GMP_OR_MPFR%.%EXT% %MANGLED%
    
copy /-Y %MANGLED% %ROOT%\mangled_binaries

goto :EOF

:ERROR

echo Usage:
echo.
echo   create_mangled_versions folder_with_gmp_or_mpfr
echo.
echo This batch script scans the given directory searching for gmp or mpfr libs
echo copying each such lib with a special mangled name based on
echo the subdirectory path where the .lib was found.
echo.
echo For example, if it finds the file
echo.
echo   gmp_4.2.4/build.vc9/lib_gmp_gc/Win32/Debug/gmp.lib
echo.
echo it creates a copy of it named 
echo.
echo   gmp_vc90_mt_gd.lib
echo.
