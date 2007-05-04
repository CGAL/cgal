;============================
; Copyright 2007 GeometryFactory (France)
; Author: Andreas Fabri (andreas.fabri@geometryfactrory.com), Fernando Cacciola (fernando.cacciola@geometryfactrory.com)
;============================
; Some portions of this file have been extracted/derived from "boost,.nsi", the Boost Windows Installer, contributed by www.boost-consulting.org.
;
; Copyright 2006 Daniel Wallin
; Copyright 2006 Eric Niebler
; Distributed under the Boost Software License, Version 1.0. (See
; accompanying file LICENSE_1_0.txt or copy at
; http://www.boost.org/LICENSE_1_0.txt)
;============================

;Include Modern UI

!include "MUI.nsh"
!include "WriteEnvStr.nsh"
!include "StrFunc.nsh"
!include "Sections.nsh"
!include "LogicLib.nsh"

!include "TextLog.nsh"
!include "script_cgal_3_3.nsh"


;-------------------------------------------------------------------------------------------------------
;
;                                          -= SOURCES =-
;
;
; The following defintions specify the source folders for the files to install. 
; ALL the variants for the precompiled libraries for ALL dependencies must exist in the source folders.
;
; These are:
;   LIB-vc71-mt-gd.lib 
;   LIB-vc71-mt-s.lib
;   LIB-vc71-mt-sgd.lib
;   LIB-vc71-mt.lib
;   LIB-vc71-s.lib
;   LIB-vc71-sgd.lib
;   LIB-vc80-mt-gd.lib
;   LIB-vc80-mt-s.lib
;   LIB-vc80-mt-sgd.lib
;   LIB-vc80-mt.lib
;
; where LIB is: cgal,CGALcore++,gmp,mpfr,taucs and zlib
;
; Since the files for third-party libraries MPFR-GMP,TAUCS and ZLIB are accessed during NSIS compilation
; only to add them to the installer package, this script searches for files related to a given third party library
; in the single defined here (so all of .h, .lib, .txt etc must be in the same folder)
;--------------------------------

!define CGAL_SRC  "CGAL-3.3"
!define GMP_SRC   "gmp-4.1.2-gladman"
!define TAUCS_SRC "taucs_precompiled"
!define ZLIB_SRC  "zlib123"


;--------------------------------
; General
;--------------------------------

  ;Name and file
  Name "CGAL-3.3"
  OutFile "CGAL-3.3-Setup.exe"

  ;Default installation folder
  InstallDir "$PROGRAMFILES\CGAL-3.3"

  ;Get installation folder from registry if available
  InstallDirRegKey HKCU "Software\CGAL-3.3" ""

;--------------------------------
; Variables
;--------------------------------

  Var MUI_TEMP
  Var STARTMENU_FOLDER
  
;--------------------------------
; Interface Settings
;--------------------------------

  !define MUI_ICON  "cgal-16.ico"
  !define MUI_UNICON  "cgal-16.ico"
  !define MUI_HEADERIMAGE
  !define MUI_HEADERIMAGE_BITMAP_NOSTRETCH

  !define MUI_HEADERIMAGE_BITMAP "cgal_very_small_FFFFFF.bmp" ; optional
 
  !define MUI_ABORTWARNING

  !define MUI_WELCOMEFINISHPAGE_BITMAP zirkel.bmp
  !define MUI_WELCOMEFINISHPAGE_BITMAP_NOSTRETCH

;  !define MUI_PAGE_HEADER_TEXT "page header text"

  !define MUI_WELCOMEPAGE_TEXT "This installs CGAL-3.3 on your machine, precompiled with .Net 2003 and 2005 (VC 7.1 and 8.0)."


  !define MUI_FINISHPAGE_TITLE "Installation Not Finished Yet!!!"

  !define MUI_FINISHPAGE_TEXT "CGAL needs part of the Boost Library which you must download yourself (from www.boost.org or www.boost-consulting.com/download.html).\r\n\r\nThe vcproj files in the example and demo subdirectories do not specify where the header files and libraries for CGAL, boost, GMP, zlib, Qt, and Taucs are located.  \r\n\r\nFor .net 2003 (VC++ 7.1) the installer already modified the Developer Studio settings for the CGAL libraries and the third party libraries GMP, Taucs, and Zlib.\r\n\r\nFor Boost and Qt you can add these paths for Developer Studio in Tools->Options->Projects->VC++Directories in the tab 'Show Directories' for 'Include Files' and 'Library Files'.\r\n\r\nVC 8.0 users must go to the CGAL-3.3/src directory and recompile the libraries."

  !define MUI_FINISHPAGE_LINK "More information about CGAL and Visual C++"
  !define MUI_FINISHPAGE_LINK_LOCATION http://www.cgal.org/platforms_frame.html
;--------------------------------
; Pages
;--------------------------------

  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_LICENSE "${CGAL_SRC}\LICENSE"
  
  ; A page were the user can specify a default variant configuration (taken from the boost installer)
  Page custom defaultVariantsPage 
  
  !insertmacro MUI_PAGE_COMPONENTS
  
  ; A page were the user can check/uncheck the enviroment variables
  ; used to specify paths in vproj files to be added.
  Page custom envarsPage 
  
  !insertmacro MUI_PAGE_DIRECTORY
  
  ;Start Menu Folder Page Configuration
  !define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU" 
  !define MUI_STARTMENUPAGE_REGISTRY_KEY "CGAL-3.3" 
  !define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Start Menu Folder"
  
  !insertmacro MUI_PAGE_STARTMENU Application $STARTMENU_FOLDER
  
  !insertmacro MUI_PAGE_INSTFILES
  
  
  !insertmacro MUI_PAGE_FINISH

  !insertmacro MUI_UNPAGE_WELCOME
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES
  !insertmacro MUI_UNPAGE_FINISH

;--------------------------------
; Languages

  !insertmacro MUI_LANGUAGE "English"

;--------------------------------
; Sections
;--------------------------------


;--------------------------------
Section "!Main CGAL" MAIN_Idx

  SectionIn RO 
  SetOutPath "$INSTDIR\auxiliary"
  FILE /r "${CGAL_SRC}\auxiliary\*.*"
  SetOutPath "$INSTDIR\config"
  FILE /r "${CGAL_SRC}\config\*.*"
  SetOutPath "$INSTDIR\include"
  FILE /r "${CGAL_SRC}\include\*.*"
  SetOutPath "$INSTDIR\make"
  FILE /r "${CGAL_SRC}\make\*.*"
  SetOutPath "$INSTDIR\src"
  FILE /r "${CGAL_SRC}\src\*.*"
  SetOutPath "$INSTDIR\examples"
  FILE /r "${CGAL_SRC}\examples\*.*"
  SetOutPath "$INSTDIR\scripts"
  FILE /r "${CGAL_SRC}\scripts\*.*"

  SetOutPath "$INSTDIR"
  FILE "${CGAL_SRC}\Changes"
  FILE "${CGAL_SRC}\Install"
  FILE "${CGAL_SRC}\INSTALL.win32"
  FILE "${CGAL_SRC}\install_cgal"
  FILE "${CGAL_SRC}\License"
  FILE "${CGAL_SRC}\LICENSE.LGPL"
  FILE "${CGAL_SRC}\LICENSE.QPL"
  FILE "${CGAL_SRC}\LICENSE.Free_Use"
  FILE "${CGAL_SRC}\Readme"
  FILE ".\cgal.ico"

  !insertmacro MUI_STARTMENU_WRITE_BEGIN  Application  
    ;Create shortcuts
    CreateDirectory "$SMPROGRAMS\$STARTMENU_FOLDER"

    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  
  !insertmacro MUI_STARTMENU_WRITE_END
    ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"
SectionEnd
;--------------------------------

;--------------------------------
Section "CGAL Demos" DEMOS_Idx

  SetOutPath "$INSTDIR\demo"
  FILE /r "${CGAL_SRC}\demo\*.*"
SectionEnd

;--------------------------------

;--------------------------------
; Multi Variant Sections
; Each of the sections below is a group enclosing all the variants for a given set of precomp libraries
; NOTE: The variant selection code uses the trailing "libs" in the group name to identify components.
;       DO NOT change the trailing "libs" in the section name.
;
${MultiVariantSection} "CGAL precomp libs"                      Install_CGAL_libs CGAL_Idx
${MultiVariantSection} "GMP and MPFR headers and precomp libs"  Install_GMP_MPFR  GMP_Idx
${MultiVariantSection} "TAUCS headers and precomp libs"         Install_TAUCS     TAUCS_Idx
${MultiVariantSection} "Zlib headers and precomp libs"          Install_ZLIB      ZLIB_Idx
;--------------------------------

;--------------------------------
;Uninstaller Section

Section "Uninstall"

  ;ADD YOUR OWN FILES HERE...

  Delete "$INSTDIR\Uninstall.exe"

  RMDir /r "$INSTDIR"
  
  !insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP
    
  Delete "$SMPROGRAMS\$MUI_TEMP\Uninstall.lnk"
  
  ;Delete empty start menu parent diretories
  StrCpy $MUI_TEMP "$SMPROGRAMS\$MUI_TEMP"
 
  startMenuDeleteLoop:
	ClearErrors
    RMDir $MUI_TEMP
    GetFullPathName $MUI_TEMP "$MUI_TEMP\.."
    
    IfErrors startMenuDeleteLoopDone
  
    StrCmp $MUI_TEMP $SMPROGRAMS startMenuDeleteLoopDone startMenuDeleteLoop
  startMenuDeleteLoopDone:

  DeleteRegKey /ifempty HKCU "Software\CGAL-3.3"

SectionEnd




;--------------------------------
;Descriptions

  ;Language strings
  LangString DESC_MAIN  ${LANG_ENGLISH} "The main components of the CGAL Library."
  LangString DESC_DEMOS ${LANG_ENGLISH} "The CGAL demos, for which you will need Qt 3 in order to build them."
  LangString DESC_CGAL  ${LANG_ENGLISH} "The precompiled CGAL libraries."
  LangString DESC_GMP   ${LANG_ENGLISH} "The precompiled GMP and MPFR libraries, which provide exact number types."
  LangString DESC_TAUCS ${LANG_ENGLISH} "The precompiled TAUCS library which provides a solver for sparse matrices that can be used together with the surface parametrization package."
  LangString DESC_ZLIB  ${LANG_ENGLISH} "The precompiled ZLIB library which provides compression algorithms that are used in examples of the surface mesher package."

  ;Assign language strings to sections
  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
    !insertmacro MUI_DESCRIPTION_TEXT ${MAIN_Idx}  $(DESC_MAIN)
    !insertmacro MUI_DESCRIPTION_TEXT ${DEMOS_Idx} $(DESC_DEMOS)
    !insertmacro MUI_DESCRIPTION_TEXT ${CGAL_Idx}  $(DESC_CGAL)
    !insertmacro MUI_DESCRIPTION_TEXT ${GMP_Idx}   $(DESC_GMP)
    !insertmacro MUI_DESCRIPTION_TEXT ${TAUCS_Idx} $(DESC_TAUCS)
    !insertmacro MUI_DESCRIPTION_TEXT ${ZLIB_Idx}  $(DESC_ZLIB)
  !insertmacro MUI_FUNCTION_DESCRIPTION_END





;--------------------------------
;Uninstaller Section
;--------------------------------

Section "Uninstall"

  ;ADD YOUR OWN FILES HERE...

  Delete "$INSTDIR\Uninstall.exe"

  RMDir "$INSTDIR"

  DeleteRegKey /ifempty HKCU "Software\CGAL-3.3"

SectionEnd


;--------------------------------
; Functions
;--------------------------------

Function .onInit

	${LogSetFileName} "install_log.txt"
	${LogSetOn}
  ${LogText} "CGAL Installer started"

  # the plugins dir is automatically deleted when the installer exits
  InitPluginsDir
  File /oname=$PLUGINSDIR\splash.bmp ".\CGAL.bmp"
  advsplash::show 1000 600 400 -1 $PLUGINSDIR\splash
  
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "default_variants.ini"
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "enviroment_variables.ini"

  #Call initSelectionFlags

FunctionEnd

Function .onSelChange
    ${LogText} "onSelChange. Old selected_libs=$selected_libs"
    ClearErrors
    StrCpy $0 0 ; Section index
    StrCpy $1 0 ; Lib index
  next:
    SectionGetText $0 $2
    ${LogText} "Section $0 is '$2'"
    IfErrors bail
    StrCpy $3 $2 "" -4
    StrCmp $3 "libs" 0 not_lib
    ${LogText} "  Section $0 is a library"
    StrCpy $3 $selected_libs 1 $1 ; $3 == old flag
    ${LogText} "  Old selection state $3"
    SectionGetFlags $0 $4 ; $4 == flag
    IntOp $5 $4 & 65
    ${LogText} "  New selection state $5"
    StrCmp $5 0 not_true
    StrCpy $5 1
  not_true:
    StrCmp $3 $5 not_toggled 0
    ${LogText} "  Selection state CHANGED"
    StrCpy $6 $selected_libs $1 ; Before
    IntOp $7 $1 + 1
    StrCpy $7 $selected_libs "" $7 ; After
    StrCpy $selected_libs "$6$5$7"
    StrCmp $5 1 0 not_selected
    ; -- New library was selected, select default variants
    Push $0
    call SelectDefaultVariants
  not_selected:
  not_toggled:
    IntOp $1 $1 + 1
  not_lib:
    IntOp $0 $0 + 1
    goto next
  bail:
FunctionEnd

Function defaultVariantsPage

    !insertmacro MUI_HEADER_TEXT "Select Default Variants" "Choose which binary variants that will be default selected for libraries."

    !insertmacro MUI_INSTALLOPTIONS_INITDIALOG "default_variants.ini"
    Pop $0
    !insertmacro MUI_INSTALLOPTIONS_SHOW
    Pop $0

    !insertmacro MUI_INSTALLOPTIONS_READ $0 "default_variants.ini" "Field 4" "State"
    ${If} $0 == 0
        !insertmacro MUI_INSTALLOPTIONS_READ $0 "default_variants.ini" "Field 5" "State"
    ${EndIf}

    IntOp $0 $0 !
    StrCpy $no_default_compilers $0
    StrCpy $0 6
    StrCpy $1 0
    ClearErrors
  next:
    !insertmacro MUI_INSTALLOPTIONS_READ $2 "default_variants.ini" "Field $0" "State"
    IfErrors bail
    IntOp $1 $1 || $2
    IntOp $0 $0 + 1
    goto next
  bail:
    IntOp $1 $1 !
    StrCpy $no_default_variants $1
 
    ${LogText} "no_default_compilers=$no_default_compilers" 
    ${LogText} "no_default_variants=$no_default_variants" 
    
    call initSelectionFlags
FunctionEnd

# Disables the env var checkbox # FN
!macro DisableEnvStrCB FN
  !insertmacro MUI_INSTALLOPTIONS_WRITE "enviroment_variables.ini" "Field ${FN}" "State" "0"
  !insertmacro MUI_INSTALLOPTIONS_WRITE "enviroment_variables.ini" "Field ${FN}" "Flags" "DISABLED"
!macroend

!macro ProcessEnvStrCB ALLUSERS FN VALUE
  !insertmacro MUI_INSTALLOPTIONS_READ $8 "enviroment_variables.ini" "Field ${FN}" "State"
  !insertmacro MUI_INSTALLOPTIONS_READ $9 "enviroment_variables.ini" "Field ${FN}" "Text"
  ${If} $8 = 1
    ${If} ${ALLUSERS} = 1
      ${LogText} "Setting enviroment variable $9='${VALUE}' for All Users."
    ${Else}
      ${LogText} "Setting enviroment variable $9='${VALUE}' for Current User Only."
    ${Endif}
  ${EndIF}
!macroend

Function envarsPage

  !insertmacro MUI_HEADER_TEXT "Setting Enviroment Variables" "Choose which enviroment variables to set and for which users."
  
  # INITIALIZATION - Disables checkboxes corresponding to env vars that cannot be set
  # either because they are already defined in the system or because the corresponding library
  # was not installed.
  ${For} $1 5 15
  
    # The text of the checkbox matches exactly the enviroment variable to be set
    !insertmacro MUI_INSTALLOPTIONS_READ $2 "enviroment_variables.ini" "Field $1" "Text"
    
    # Disable the checkbox if the envar is alreadyc defined
    ReadEnvStr $3 $2
    ${If} $3 != ""
      !insertmacro DisableEnvStrCB $1
    ${Endif}
    
    # Disable checkboxes corresponding to optional components that might have not be installed
    StrCpy $4 $2 1
        
    ${Select} $4
      ${Case} "C" # checkboxes for CGAL libs
        ${Unless}   ${SectionIsSelected} ${CGAL_Idx}
        ${AndUnless} ${SectionIsPartiallySelected} ${CGAL_Idx}
          !insertmacro DisableEnvStrCB $1
        ${EndUnless}
      ${Case2} "M" "G" # checkboxes for MPFR and GMP
        ${Unless}    ${SectionIsSelected} ${GMP_Idx}
        ${AndUnless} ${SectionIsPartiallySelected} ${GMP_Idx}
          !insertmacro DisableEnvStrCB $1
        ${EndUnless}
      ${Case} "T" # checkboxes for TAUCS
        ${Unless}    ${SectionIsSelected} ${TAUCS_Idx}
        ${AndUnless} ${SectionIsPartiallySelected} ${TAUCS_Idx}
          !insertmacro DisableEnvStrCB $1
        ${EndUnless}
      ${Case} "Z" # checkboxes for ZLIB
        ${Unless}    ${SectionIsSelected} ${ZLIB_Idx}
        ${AndUnless} ${SectionIsPartiallySelected} ${ZLIB_Idx}
          !insertmacro DisableEnvStrCB $1
        ${EndUnless}
    ${EndSelect}
  ${Next}
  
  !insertmacro MUI_INSTALLOPTIONS_INITDIALOG "enviroment_variables.ini"
  Pop $0

  !insertmacro MUI_INSTALLOPTIONS_SHOW
  Pop $0

  # PROCESSING - Installs selected enviroment variables
  
  !insertmacro MUI_INSTALLOPTIONS_READ $1 "enviroment_variables.ini" "Field 2" "State" # $1=Is ALL USERS selected
  
  !insertmacro ProcessEnvStrCB $1 5  "$INSTDIR\include"    
  !insertmacro ProcessEnvStrCB $1 6  "$INSTDIR\lib"
  !insertmacro ProcessEnvStrCB $1 7  "boost-root"    
  !insertmacro ProcessEnvStrCB $1 8  "$INSTDIR\auxiliary\gmp\include"
  !insertmacro ProcessEnvStrCB $1 9  "$INSTDIR\auxiliary\gmp\lib"    
  !insertmacro ProcessEnvStrCB $1 10 "$INSTDIR\auxiliary\gmp\include"    
  !insertmacro ProcessEnvStrCB $1 11 "$INSTDIR\auxiliary\gmp\lib"    
  !insertmacro ProcessEnvStrCB $1 12 "$INSTDIR\auxiliary\taucs\include"    
  !insertmacro ProcessEnvStrCB $1 13 "$INSTDIR\auxiliary\taucs\lib"    
  !insertmacro ProcessEnvStrCB $1 14 "$INSTDIR\auxiliary\zlib\include"    
  !insertmacro ProcessEnvStrCB $1 15 "$INSTDIR\auxiliary\zlib\lib"    

FunctionEnd
