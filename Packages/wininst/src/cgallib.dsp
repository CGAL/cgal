# Microsoft Developer Studio Project File - Name="cgallib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=cgallib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "cgallib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "cgallib.mak" CFG="cgallib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cgallib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "cgallib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "cgallib - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "NDEBUG" /D "_MBCS" /D "_LIB" /D "WIN32" /YX /FD /TP /c
# ADD CPP /nologo /W3 /GX /GR /O2 /I "$(CGALROOT)\stlport" /I "$(CGALROOT)\auxilary\wingmp\gmp-4.0.1" /I "$(CGALROOT)\include\cgal\config\msvc6" /I "$(CGALROOT)\include" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "CGAL_USE_GMP" /YX /FD /TP /c
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"../lib/msvc6/cgal.lib"

!ELSEIF  "$(CFG)" == "cgallib - Win32 Debug"

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
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /TP /D /GZ "WIN32" /c
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /I "$(CGALROOT)\stlport" /I "$(CGALROOT)\auxilary\wingmp\gmp-4.0.1" /I "$(CGALROOT)\include\cgal\config\msvc6" /I "$(CGALROOT)\include" /D "_DEBUG" /D "CGAL_USE_GMP" /D "WIN32" /D "_MBCS" /D "_LIB" /FR /YX /FD /TP /GZ /c
# ADD BASE RSC /l 0x40c /d "_DEBUG"
# ADD RSC /l 0x40c /i "$(CGALROOT)\stlport" /i "$(CGALROOT)\include\cgal\config\msvc6" /i "$(CGALROOT)\auxilary\wingmp\gmp-4.0.1" /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"../lib/msvc6/cgal.lib"

!ENDIF 

# Begin Target

# Name "cgallib - Win32 Release"
# Name "cgallib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "C;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\aff_transformation_tags.C
# End Source File
# Begin Source File

SOURCE=.\assertions.C
# End Source File
# Begin Source File

SOURCE=.\Bbox_2_intersections.C
# End Source File
# Begin Source File

SOURCE=.\Bbox_3_intersections.C
# End Source File
# Begin Source File

SOURCE=.\cgal_logo.C
# End Source File
# Begin Source File

SOURCE=.\Color.C
# End Source File
# Begin Source File

SOURCE=.\File_header_extended_OFF.C
# End Source File
# Begin Source File

SOURCE=.\File_header_OFF.C
# End Source File
# Begin Source File

SOURCE=.\File_scanner_OFF.C
# End Source File
# Begin Source File

SOURCE=.\File_writer_inventor.C
# End Source File
# Begin Source File

SOURCE=.\File_writer_OFF.C
# End Source File
# Begin Source File

SOURCE=.\File_writer_VRML_2.C
# End Source File
# Begin Source File

SOURCE=.\File_writer_wavefront.C
# End Source File
# Begin Source File

SOURCE=.\Geomview_stream.C
# End Source File
# Begin Source File

SOURCE=.\Interval_arithmetic.C
# End Source File
# Begin Source File

SOURCE=.\io.C
# End Source File
# Begin Source File

SOURCE=.\MP_Float.C
# End Source File
# Begin Source File

SOURCE=.\optimisation_basic.C
# End Source File
# Begin Source File

SOURCE=.\Origin.C
# End Source File
# Begin Source File

SOURCE=.\Random.C
# End Source File
# Begin Source File

SOURCE=.\Real_timer.C
# End Source File
# Begin Source File

SOURCE=.\Timer.C
# End Source File
# Begin Source File

SOURCE=.\Triangulation_3.C
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# End Target
# End Project
