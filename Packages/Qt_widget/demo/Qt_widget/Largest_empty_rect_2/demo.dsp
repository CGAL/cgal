# Microsoft Developer Studio Project File - Name="demo" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

CFG=demo - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "demo.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "demo.mak" CFG="demo - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "demo - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "demo - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "demo - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386

!ELSEIF  "$(CFG)" == "demo - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /I "$(QTDIR)\include" /I "$(CGALROOT)\stlport" /I "$(CGALROOT)\include\cgal\config\msvc6" /I "$(CGALROOT)\auxiliary\wingmp\gmp-4.0.1" /I "$(CGALROOT)\include" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "CGAL_USE_QT" /D "QT_DLL" /D "UNICODE" /D "QT_THREAD_SUPPORT" /FR /YX /FD /GZ /c /Tp
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x40c /d "_DEBUG"
# ADD RSC /l 0x40c /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib imm32.lib winmm.lib wsock32.lib imm32.lib wsock32.lib winmm.lib $(QTDIR)\lib\qt-mt230nc.lib $(QTDIR)\lib\qtmain.lib $(CGALROOT)\lib\msvc6\CGAL.lib $(CGALROOT)\lib\msvc6\CGALQt.lib /nologo /subsystem:windows /debug /machine:I386 /nodefaultlib:"libc" /nodefaultlib:"msvcrt.lib" /nodefaultlib:"libcmtd" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "demo - Win32 Release"
# Name "demo - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\largest_empty_rect_2.C
# End Source File
# Begin Source File

SOURCE=.\Qt_widget_move_list_point.C
# End Source File
# Begin Source File

SOURCE=.\Qt_widget_toolbar.C
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\Qt_widget_move_list_point.h

!IF  "$(CFG)" == "demo - Win32 Release"

!ELSEIF  "$(CFG)" == "demo - Win32 Debug"

# Begin Custom Build
InputPath=.\Qt_widget_move_list_point.h

"Qt_widget_move_list_point.moc" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	$(QTDIR)\bin\moc.exe -o "Qt_widget_move_list_point.moc" "Qt_widget_move_list_point.h"

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Qt_widget_toolbar.h

!IF  "$(CFG)" == "demo - Win32 Release"

!ELSEIF  "$(CFG)" == "demo - Win32 Debug"

# Begin Custom Build
InputPath=.\Qt_widget_toolbar.h

BuildCmds= \
	$(QTDIR)\bin\moc.exe -o "Qt_widget_toolbar.moc" "Qt_widget_toolbar.h" \
	$(QTDIR)\bin\moc.exe -o largest_empty_rect_2.moc largest_empty_rect_2.C \
	

"Qt_widget_toolbar.moc" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
   $(BuildCmds)

"largest_empty_rect_2.moc" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ENDIF 

# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
