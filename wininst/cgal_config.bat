@echo off

rem %1: compiler
rem %2: compiler options (optional)
rem %3: LEDA libs
rem %4: LEDA includes

rem -----------------------------------------------------
rem 		Only Windows!
rem -----------------------------------------------------

if "%windir%"  == "" goto :unknos

rem -----------------------------------------------------
rem		compiler settings
rem -----------------------------------------------------

set cc=unknown

if "%1" == "msc"	set cc=msvc
if "%1" == "bcc"	set cc=bcc

if %cc% == unknown goto :usage

set ccopt=
set ledaroot=%2
set ledaincl=%3

if %cc% == bcc	goto :bccconf
if %cc% == msvc	goto :msvcconf

:msvcconf

set cclibo=unknown
if "%2" == "ml"		set cclibo=ML
if "%2" == "mt"		set cclibo=MT
if "%2" == "md"		set cclibo=MD
if "%2" == "mld"	set cclibo=MLd
if "%2" == "mtd"	set cclibo=MTd
if "%2" == "mdd"	set cclibo=MDd
if %cclibo% == unknown set cclibo=ML

set ccopt=unknown
if "%2" == "ml"		set ccopt=-ML
if "%2" == "mt"		set ccopt=-MT
if "%2" == "md"		set ccopt=-MD
if "%2" == "mld"	set ccopt=-MLd -Z7
if "%2" == "mtd"	set ccopt=-MTd -Z7
if "%2" == "mdd"	set ccopt=-MDd -Z7
if %cclibo% == ML	set ccopt=-ML

set cxx=cl -nologo -GX 
set make=nmake -nologo /S
set extralibs=kernel32.lib
set cgallibpref=-LIBPATH:
goto :ledaconf

:bccconf
set cclibo=nd
if "%2" == "d"		set ccopt=-v
if "%2" == "d"		set cclibo=d

set cxx=bcc32
rem set cxx=bcc32 -nologo -- BCC 5.5.1 does not have -nologo option!
set make=make -s -N
set extralibs=
set cgallibpref=-L

rem --------------------------------------------------------
rem	5.4 or 5.5 ?
rem --------------------------------------------------------
if exist tmpt0.cpp del tmpt0.cpp > NUL
echo #if ( __BORLANDC__ == 0x540 ) > tmpt0.cpp
echo #error 1 >> tmpt0.cpp
echo #endif >> tmpt0.cpp
%cxx% -c tmpt0.cpp >NUL
if ERRORLEVEL 1 goto :bcc54
if exist include\CGAL\config\bcc\cctype del include\CGAL\config\bcc\cctype >NUL
goto :ledaconf

:bcc54
copy include\CGAL\config\bcc\cctype.txt include\CGAL\config\bcc\cctype >NUL
goto :ledaconf


rem --------------------------------------------------------
rem	LEDA, or no LEDA ? And if yes, then where?
rem --------------------------------------------------------

:ledaconf
rem --------------------------------------------------------
rem	checking if the compiler works...
rem --------------------------------------------------------
if exist tmpt0.cpp del tmpt0.cpp > NUL
echo int main() {return 0;} > tmpt0.cpp
%cxx% -c tmpt0.cpp >NUL
if ERRORLEVEL 1 goto :nocxx


if "%ccopt%" == "" goto :ledacnfa
if "%ledaroot%" == "" goto :noleda

set ledaroot=%3
set ledaincl=%4

rem --------------------------------------------------------
rem	LEDA dirs
rem --------------------------------------------------------

:ledacnfa
if "%ledaroot%" == "" goto :noleda
if "%ledaincl%" == "" goto :dfltinc
goto :conf1

:dfltinc
set ledaincl=%ledaroot%\incl

:conf1
echo.
echo Configuring with LEDA libs in %ledaroot% and headers in %ledaincl%
echo.

if not exist %ledaroot%\libL.lib goto :noledar
if not exist %ledaincl%\LEDA\basic.h goto :noledai

rem --------------------------------------------------------
rem	checking LEDA_STD_HEADERS
rem --------------------------------------------------------
if exist tmpt0.cpp del tmpt0.cpp >NUL
copy winutils\src\ledatst.cpp tmpt0.cpp >NUL
%cxx% -I%ledaincl% -DLEDA_PREFIX -c tmpt0.cpp >NUL
if ERRORLEVEL 1 goto :noledast

set ledain=leda
set ledaflag=-DCGAL_USE_LEDA -DLEDA_PREFIX -I$(LEDA_INCL_DIR)
set ledasupport="SUPPORTED"
set ledalink=%cgallibpref%%ledaroot%
set ledalibs=libP.lib libG.lib libL.lib libD3.lib libGeoW.lib libW.lib
set cgalwindirpath=
goto :done

rem --------------------------------------------------------
rem	no LEDA
rem --------------------------------------------------------

:noleda
echo.
echo Configuring without LEDA
echo.
set ledain=noleda
set ledaroot=
set ledaincl=
set ledaflag=
set ledalink=
set ledasupport="NOT SUPPORTED"
set ledalibs=CGALWin.lib
set cgalwindirpath=%cgallibpref%$(CGAL_ROOT)\lib\win32\%cc%\%cclibo%
goto :done


rem --------------------------------------------------------
rem	error messages...
rem --------------------------------------------------------

:usage
echo.
echo Usage: cgal_config cc [comp. opts] [LEDA libs dir] [LEDA includes dir]
echo.
echo Possible compilers (cc):
echo msc    : Microsoft Visual C++ 6.0	(opts: ml mld md mt mdd mtd)
echo bcc    : Borland C++ Builder 4, C++ 5.4, 5.5	(opts: d)
echo.
echo   examples: "> cgal_config msc mld g:\LEDA"
echo             "> cgal_config msc g:\LEDA"
echo             "> cgal_config msc mt"
goto :quit

:noledar
echo.
echo LEDA libs are not present in %ledaroot% !
echo.
goto :quit

:noledai
echo.
echo LEDA headers are not present in %ledaincl% !
echo.
goto :quit


:noledast
echo. 
echo LEDA appear to be not configured with LEDA_STD_HEADERS !
echo In this case it does not work with CGAL. Details in
echo CGAL Installation guide and in
echo  http://www.mpi-sb.mpg.de/LEDA/download/windows
echo.
goto :quit

:unknos
echo.
echo This does not appear to be a Windows 9* or Windows NT machine!
echo Quitting...
goto :quit

:nocxx
echo.
echo compiler call "%cxx% -c tmpt0.cpp" failed. Here is tmpt0.cpp:
echo        int main() {return 0;}
echo Probably the compiler cannot be found in your PATH.
echo.
goto :quit

rem --------------------------------------------------------
rem	Dumping header of the makefile
rem --------------------------------------------------------
:done

rem --------------------------------------------------------
rem 	Building the necessary tools
rem --------------------------------------------------------

if exist winutils\bin\pwd.exe goto :donetls
echo Building tools...
if exist winutils\bin goto :nombindi
mkdir winutils\bin >NUL

:nombindi

cd winutils\src
%make% all > NUL
cd ..\..
echo Done.
echo.
:donetls

echo # This file contains CGAL makefile.mak settings for the following platform: > makefile.mak
echo # OS: %PROCESSOR_ARCHITECTURE% %OS% >> makefile.mak
echo # COMPILER:        %cc% >> makefile.mak
echo # GMP:             supported >> makefile.mak
echo # LEDA:            %ledasupport% >> makefile.mak
echo #---------------------------------------------------------------------# >> makefile.mak
echo #                    installation directory >> makefile.mak
echo #---------------------------------------------------------------------#>> makefile.mak 
winutils\bin\pwd CGAL_ROOT=>> makefile.mak
echo #---------------------------------------------------------------------# >> makefile.mak
echo #                    os/compiler description >> makefile.mak
echo #---------------------------------------------------------------------# >> makefile.mak
echo CGAL_OS_COMPILER = %cc%>> makefile.mak
echo.  >> makefile.mak
echo # LEDA include directory *** >> makefile.mak
echo LEDA_INCL_DIR = %ledaincl% >> makefile.mak
echo.  >> makefile.mak
echo # LEDA libs directory *** >> makefile.mak
echo LEDA_LIB_DIR = %ledaroot% >> makefile.mak
echo.  >> makefile.mak
echo # LEDA-specific compilation flags *** >> makefile.mak
echo LE_CXXFLAGS = %ledaflag% >> makefile.mak
echo.  >> makefile.mak
echo # LEDA-specific linking flags *** >> makefile.mak
echo LE_LIB_DIR = %ledalink% >> makefile.mak
echo LE_LIBS_LIST = %ledalibs% >> makefile.mak
echo.  >> makefile.mak
echo # CGAL_WINDOW library path *** >> makefile.mak
echo CGALWIN_LIB_DIR_PATH = %cgalwindirpath% >> makefile.mak
echo. >> makefile.mak
rem --------------------------------------------------------
rem	setting compiler/linker options
rem --------------------------------------------------------
echo # *** Extra compiler flags  *** >> makefile.mak
echo CUSTOM1_CXXFLAGS = %ccopt% >> makefile.mak
echo # *** Extra linker flags  *** >> makefile.mak
echo CUSTOM1_LDFLAGS = %ccopt% >> makefile.mak

rem --------------------------------------------------------
rem	attaching header and the right makefile.mak
rem --------------------------------------------------------
copy makefile.mak + winutils\make\%cc%\makefile makefile.mak > NUL

rem --------------------------------------------------------
rem	creating make_lib batchfile
rem --------------------------------------------------------
echo rem batchfile to create CGAL libs > make_lib.bat
echo. >> make_lib.bat
echo rem this file is autmatically generated, do not edit. >> make_lib.bat
echo. >> make_lib.bat
rem --------------------------------------------------------
rem	setting CGAL_MAKEFILE there
rem --------------------------------------------------------
winutils\bin\pwd set CGAL_MAKEFILE= \makefile.mak >> make_lib.bat
echo. >> make_lib.bat
rem --------------------------------------------------------
rem	cleaning up
rem --------------------------------------------------------
echo if exist lib\%cc%\ goto :direxi  >> make_lib.bat
echo mkdir lib\%cc% >> make_lib.bat
echo :direxi >> make_lib.bat
echo if exist lib\%cc%\CGAL.lib del lib\%cc%\CGAL.lib  >> make_lib.bat
rem --------------------------------------------------------
rem	building lib
rem --------------------------------------------------------
echo cd src >> make_lib.bat
echo %make% -f makefile_lib.mak >> make_lib.bat
echo cd .. >> make_lib.bat
echo @echo off >> make_lib.bat
echo echo Done. Now you can build examples and demo's >> make_lib.bat
echo echo.  >> make_lib.bat
echo echo Use make_demo and make_examples batch files  >> make_lib.bat
echo echo.  >> make_lib.bat
echo echo Or change to examples\"desired example"  >> make_lib.bat
echo echo and type %make% -f makefile.mak all >> make_lib.bat
echo echo or, change to demo\"desired demo"  >> make_lib.bat
echo echo     and type %make% -f makefile.mak all >> make_lib.bat
echo echo.  >> make_lib.bat 

rem ------------------------------------------------------------------
rem 		creating make_demo batchfile
rem ------------------------------------------------------------------
echo rem batchfile to compile and link CGAL demos > make_demo.bat
echo.>> make_demo.bat
echo rem this file is autmatically generated, do not edit. >> make_demo.bat
echo.>> make_demo.bat
echo set MMAKE=%make% -f makefile.mak >> make_demo.bat
echo.>> make_demo.bat
copy make_demo.bat + winutils\make\%cc%\%ledain%\demo.bat  make_demo.bat >NUL
echo @echo off >> make_demo.bat
echo set MMAKE=>> make_demo.bat
echo echo. >> make_demo.bat
echo echo Now you can run the demos is the usual way. >> make_demo.bat
echo. >> make_demo.bat


rem ------------------------------------------------------------------
rem 		creating make_examples batchfile
rem ------------------------------------------------------------------
echo rem batchfile to compile and link CGAL examples > make_examples.bat
echo. >> make_examples.bat
echo rem this file is autmatically generated, do not edit. >> make_examples.bat
echo. >> make_examples.bat
echo set MMAKE=%make% -f makefile.mak >> make_examples.bat
echo.>> make_examples.bat
copy make_examples.bat + winutils\make\%cc%\%ledain%\ex.bat make_examples.bat >NUL
echo @echo off >> make_examples.bat
echo set MMAKE=>> make_examples.bat
echo echo. >> make_examples.bat
echo echo Now you can run the examples is the usual way. >> make_examples.bat
echo. >> make_examples.bat


echo.
echo Finished configuring. Now run make_lib to build libs
echo.
goto :quit

:quit
if exist tmpt0.cpp del tmpt0.cpp > NUL
set cc=
set ccopt=
set ledaroot=
set ledaincl=
set ledaflag=
set ledasupport=
set ledalink=
set libpref=
set ledalibs=
set cxx=
set make=
set extralibs=
set ledain=
set cgalwindirpath=
set cclibo=
