# Microsoft Developer Studio Project File - Name="cgalQt" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=cgalQt - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "CGALQt.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "CGALQt.mak" CFG="cgalQt - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cgalQt - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "cgalQt - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "cgalQt - Win32 Release"

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
# ADD CPP /nologo /W3 /GX /O2 /I "$(CGALROOT)\stlport" /I "$(CGALROOT)\auxilary\wingmp\gmp-2.0.2" /I "$(CGALROOT)\include\cgal\config\msvc" /I "$(CGALROOT)\include" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "CGAL_USE_GMP" /YX /FD /TP /c
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"lib\msvc\visual\Release\cgal.lib"

!ELSEIF  "$(CFG)" == "cgalQt - Win32 Debug"

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
# ADD RSC /l 0x40c /i "$(CGALROOT)\include" /i "$(CGALROOT)\include\cgal" /i "$(CGALROOT)\stlport" /i "$(CGALROOT)\include\cgal\config\msvc" /i "$(CGALROOT)\auxilary\wingmp\gmp-2.0.2" /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"../../lib/msvc6/CGALQt.lib"

!ENDIF 

# Begin Target

# Name "cgalQt - Win32 Release"
# Name "cgalQt - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "C;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\Qt_widget.C
# End Source File
# Begin Source File

SOURCE=.\Qt_widget_layer.C
# End Source File
# Begin Source File

SOURCE=.\Qt_widget_standard_toolbar.C
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\include\CGAL\IO\Qt_widget.h

!IF  "$(CFG)" == "cgalQt - Win32 Release"

!ELSEIF  "$(CFG)" == "cgalQt - Win32 Debug"

# Begin Custom Build
InputPath=..\..\include\CGAL\IO\Qt_widget.h

"Qt_widget.moc" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	$(QTDIR)/bin/moc.exe $(CGALROOT)/include/CGAL/IO/Qt_widget.h -o Qt_widget.moc

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\include\CGAL\IO\Qt_widget_layer.h

!IF  "$(CFG)" == "cgalQt - Win32 Release"

!ELSEIF  "$(CFG)" == "cgalQt - Win32 Debug"

# Begin Custom Build
InputPath=..\..\include\CGAL\IO\Qt_widget_layer.h

"Qt_widget_layer.moc" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	$(QTDIR)/bin/moc.exe $(CGALROOT)/include/CGAL/IO/Qt_widget_layer.h -o Qt_widget_layer.moc

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\include\CGAL\IO\Qt_widget_standard_toolbar.h

!IF  "$(CFG)" == "cgalQt - Win32 Release"

!ELSEIF  "$(CFG)" == "cgalQt - Win32 Debug"

# Begin Custom Build
InputPath=..\..\include\CGAL\IO\Qt_widget_standard_toolbar.h

"Qt_widget_standard_toolbar.moc" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	$(QTDIR)/bin/moc.exe $(CGALROOT)/include/CGAL/IO/Qt_widget_standard_toolbar.h -o Qt_widget_standard_toolbar.moc

# End Custom Build

!ENDIF 

# End Source File
# End Group
# End Target
# End Project
