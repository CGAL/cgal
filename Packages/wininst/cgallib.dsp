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
# ADD CPP /nologo /W3 /GX /O2 /I "$(CGAL)\stlport" /I "$(CGAL)\auxilary\wingmp\gmp-2.0.2" /I "$(CGAL)\include\cgal\config\msvc" /I "$(CGAL)\include" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "CGAL_USE_GMP" /YX /FD /TP /c
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"lib\msvc\visual\Release\cgal.lib"

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
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /I "$(CGAL)\stlport" /I "$(CGAL)\auxilary\wingmp\gmp-2.0.2" /I "$(CGAL)\include\cgal\config\msvc" /I "$(CGAL)\include" /D "_DEBUG" /D "CGAL_USE_GMP" /D "WIN32" /D "_MBCS" /D "_LIB" /FR /YX /FD /TP /GZ /c
# SUBTRACT CPP /X
# ADD BASE RSC /l 0x40c /d "_DEBUG"
# ADD RSC /l 0x40c /i "$(CGAL)\include" /i "$(CGAL)\include\cgal" /i "$(CGAL)\stlport" /i "$(CGAL)\include\cgal\config\msvc" /i "$(CGAL)\auxilary\wingmp\gmp-2.0.2" /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"lib\msvc\visual\Debug\cgal.lib"

!ENDIF 

# Begin Target

# Name "cgallib - Win32 Release"
# Name "cgallib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "C;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE="$(CGAL)\src\aff_transformation_tags.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\assertions.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\Bbox_2.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\Bbox_2_intersections.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\Bbox_3_intersections.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\Color.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\File_header_extended_OFF.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\File_header_OFF.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\File_scanner_OFF.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\File_writer_inventor.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\File_writer_OFF.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\File_writer_VRML_2.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\File_writer_wavefront.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\Geomview_stream.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\Interval_arithmetic.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\io.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\MP_Float.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\optimisation_basic.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\Origin.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\Random.C"
# End Source File
# Begin Source File

SOURCE="$(CGAL)\src\Triangulation_3.C"
# End Source File
# Begin Source File

SOURCE=.\src\Interval_arithmetic\workaround_4_ms.c

!IF  "$(CFG)" == "cgallib - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputDir=.\src\Interval_arithmetic
OutDir=.\Release
InputPath=.\src\Interval_arithmetic\workaround_4_ms.c
InputName=workaround_4_ms

"$(OutDir)\$(InputName).obj" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	copy $(InputDir)\$(InputName).obj $(OutDir)\$(InputName).obj

# End Custom Build

!ELSEIF  "$(CFG)" == "cgallib - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputDir=.\src\Interval_arithmetic
OutDir=.\Debug
InputPath=.\src\Interval_arithmetic\workaround_4_ms.c
InputName=workaround_4_ms

"$(OutDir)\$(InputName).obj" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	copy $(InputDir)\$(InputName).obj $(OutDir)\$(InputName).obj

# End Custom Build

!ENDIF 

# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# End Target
# End Project
