# Microsoft Developer Studio Project File - Name="cgalwin" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=cgalwin - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "cgalwin.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "cgalwin.mak" CFG="cgalwin - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cgalwin - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "cgalwin - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "cgalwin - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /TP /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /TP /c
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"lib\msvc\visual\Release\cgalwin.lib"

!ELSEIF  "$(CFG)" == "cgalwin - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /TP /GZ /c
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /I "$(CGAL)\stlport" /I "$(CGAL)\auxilary\wingmp\gmp-2.0.2" /I "$(CGAL)\include\cgal\config\msvc" /I "$(CGAL)\include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "CGAL_USE_GMP" /D "CGAL_USE_CGAL_HEADERS" /YX /FD /TP /GZ /c
# ADD BASE RSC /l 0x40c /d "_DEBUG"
# ADD RSC /l 0x40c /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"lib\msvc\visual\Debug\cgalwin.lib"

!ENDIF 

# Begin Target

# Name "cgalwin - Win32 Release"
# Name "cgalwin - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE="$(CGAL)\src\CGALWin\_base_panel.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\CGALWin\_base_window.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\CGALWin\_basic.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\CGALWin\_color.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\CGALWin\_file.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\CGALWin\_file_panel.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\CGALWin\_string_manip.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\CGALWin\_window.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\CGALWin\mswin\_x_basic.C"
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# End Target
# End Project
