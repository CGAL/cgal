;============================
; Copyright 2007 GeometryFactory (France)
; Author: Andreas Fabri (andreas.fabri@geometryfactrory.com), Fernando Cacciola (fernando.cacciola@geometryfactrory.com)
;============================
; Some portions of this file have been derived from "boost.nsi", the Boost Windows Installer, contributed by www.boost-consulting.org.
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
; ALL the variants for the precompiled libraries for ALL dependencies must exist.
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
; where LIB is: cgal,CGALcore++,gmp,mpfr
;
; For MPFR/GMP, ALL files (.h, .lib, .txt etc) must be in a single folder.
;--------------------------------

!define CGAL_SRC  "CGAL-3.3"
!define GMP_SRC   "gmp-4.1.2-gladman"


;--------------------------------
; General
;--------------------------------

  ;Name and file
  Name "CGAL-3.3"
  OutFile "CGAL-3.3-Setup.exe"

  ;Default installation folder
  InstallDir "$PROGRAMFILES\CGAL-3.3a"

  ;Get installation folder from registry if available
  InstallDirRegKey HKCU "Software\CGAL-3.3" ""
  
  BrandingText "The CGAL Project and GeometryFactory"

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
 
  !define MUI_FINISHPAGE_NOAUTOCLOSE
  
  !define MUI_ABORTWARNING

  !define MUI_WELCOMEFINISHPAGE_BITMAP Zirkel.bmp 
  !define MUI_WELCOMEFINISHPAGE_BITMAP_NOSTRETCH

  !define MUI_COMPONENTSPAGE_SMALLDESC
  
;  !define MUI_PAGE_HEADER_TEXT "page header text"

  !define MUI_WELCOMEPAGE_TEXT "This installs CGAL-3.3 on your machine, precompiled with .Net 2003 and 2005 (VC 7.1 and 8.0).\r\nThe project files for building the library itself are provided for both VC7.1 and VC8.0 separatedly, but for examples and demos, only the VC7.1 project files are provided since all its settings are compatible with VC8.0 (thus you can just open them with Visual Studio 2005 and follow the conversion wizard)."


  !define MUI_FINISHPAGE_TITLE "Installation Not Finished Yet!!!"

  !define MUI_FINISHPAGE_TEXT "CGAL needs part of the Boost Library which you must download yourself (from www.boost.org or www.boost-consulting.com/download.html).\r\nMost demos need the Qt 3 library which you must also download yourself (from www.trolltech.com).\r\nThe vcproj files in the example, demo amd src subdirectories all use the enviroment variables CGALROOT and BOOSTROOT to locate CGAL and Boost resp., but for other third-party libraries like Qt3, TAUCS and ZLib, you must add the corresponding include/lib paths manually."

  !define MUI_FINISHPAGE_SHOWREADME "$INSTDIR\INSTALL.win32.txt"
  !define MUI_FINISHPAGE_SHOWREADME_TEXT "Read the full installation notes for further instructions"
  
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

!ifndef TestingOnly
  SectionIn RO 
  SetOutPath "$INSTDIR\auxiliary"
  File /r "${CGAL_SRC}\auxiliary\*.*"
  SetOutPath "$INSTDIR\config"
  File /r "${CGAL_SRC}\config\*.*"
  SetOutPath "$INSTDIR\include"
  File /r "${CGAL_SRC}\include\*.*"
  SetOutPath "$INSTDIR\make"
  File /r "${CGAL_SRC}\make\*.*"
  SetOutPath "$INSTDIR\src"
  File /r "${CGAL_SRC}\src\*.*"
  SetOutPath "$INSTDIR\examples"
  File /r "${CGAL_SRC}\examples\*.*"
  SetOutPath "$INSTDIR\scripts"
  File /r "${CGAL_SRC}\scripts\*.*"

  SetOutPath "$INSTDIR"
  File "${CGAL_SRC}\Changes"
  File "${CGAL_SRC}\Install"
  File "${CGAL_SRC}\INSTALL.win32.txt"
  File "${CGAL_SRC}\install_cgal"
  File "${CGAL_SRC}\License"
  File "${CGAL_SRC}\LICENSE.LGPL"
  File "${CGAL_SRC}\LICENSE.QPL"
  File "${CGAL_SRC}\LICENSE.Free_Use"
  File "${CGAL_SRC}\Readme"
  File ".\cgal.ico"
!endif

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

!ifndef TestingOnly
  SetOutPath "$INSTDIR\demo"
  File /r "${CGAL_SRC}\demo\*.*"
!endif
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
  LangString DESC_MAIN    ${LANG_ENGLISH} "The main components of the CGAL Library."
  LangString DESC_DEMOS   ${LANG_ENGLISH} "The CGAL demos, for which you will need Qt 3 in order to build them."
  LangString DESC_CGAL    ${LANG_ENGLISH} "The precompiled CGAL libraries."
  LangString DESC_GMP     ${LANG_ENGLISH} "The precompiled GMP and MPFR libraries, which provide exact number types."
  LangString DESC_EVARSET ${LANG_ENGLISH} "(this enviroment variable already exists in your system)"

  ;Assign language strings to sections
  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
    !insertmacro MUI_DESCRIPTION_TEXT ${MAIN_Idx}  $(DESC_MAIN)
    !insertmacro MUI_DESCRIPTION_TEXT ${DEMOS_Idx} $(DESC_DEMOS)
    !insertmacro MUI_DESCRIPTION_TEXT ${CGAL_Idx}  $(DESC_CGAL)
    !insertmacro MUI_DESCRIPTION_TEXT ${GMP_Idx}   $(DESC_GMP)
  !insertmacro MUI_FUNCTION_DESCRIPTION_END


;--------------------------------
; Functions
;--------------------------------

Function .onInit

	${LogSetFileName} "$INSTDIR\windows_install.log"
	${LogSetOn}

  # the plugins dir is automatically deleted when the installer exits
  InitPluginsDir
  File /oname=$PLUGINSDIR\splash.bmp ".\CGAL.bmp"
  advsplash::show 1000 600 400 -1 $PLUGINSDIR\splash
  
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "default_variants.ini"
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "enviroment_variables.ini"

FunctionEnd

Function .onSelChange
    ClearErrors
    StrCpy $0 0 ; Section index
    StrCpy $1 0 ; Lib index
  next:
    SectionGetText $0 $2
    IfErrors bail
    StrCpy $3 $2 "" -4
    StrCmp $3 "libs" 0 not_lib
    StrCpy $3 $selected_libs 1 $1 ; $3 == old flag
    SectionGetFlags $0 $4 ; $4 == flag
    IntOp $5 $4 & 65
    StrCmp $5 0 not_true
    StrCpy $5 1
  not_true:
    StrCmp $3 $5 not_toggled 0
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

    !insertmacro MUI_HEADER_TEXT "Select Default Variants" "Choose the default variants to autoselect in the next page."
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

    call initSelectionFlags
FunctionEnd

# Disables the env var checkbox # FN and textbox # FN+1
!macro UncheckEnvStrCheckbox FN
  !insertmacro MUI_INSTALLOPTIONS_WRITE "enviroment_variables.ini" "Field ${FN}" "State" "0"
!macroend

# Disables the env var checkbox # FN
!macro SetEnvStrValueSlot FN VAL
  !insertmacro MUI_INSTALLOPTIONS_WRITE "enviroment_variables.ini" "Field ${FN}" "State" "${VAL}"
!macroend

!macro SetEnvStrLabel FN VAL
  !insertmacro MUI_INSTALLOPTIONS_WRITE "enviroment_variables.ini" "Field ${FN}" "Text" "${VAL}"
!macroend

!macro ProcessEnvStrCB ALLUSERS VAR FN

  # ${ALLUSERS} is 0 or 1
  # ${VAR} is the env var to set
  # ${FN} is the filed number of the checkbox corresponding to the env var
  
  # $6  is the filed number of the textbox corresponding to the env var value 
  # $7  is the state of checkbox
  # $8  is the value of the env var
  
  IntOp $6 ${FN} + 1
  
  !insertmacro MUI_INSTALLOPTIONS_READ $7 "enviroment_variables.ini" "Field ${FN}" "State" 
  !insertmacro MUI_INSTALLOPTIONS_READ $8 "enviroment_variables.ini" "Field $6"    "State"
  
  ${If} $7 = 1 # checkbox selected - set the env var
    ${If} ${ALLUSERS} = 1 
      ${LogText} "Setting enviroment variable ${VAR}='$8' for All Users."
      !define ALL_USERS
      !ifndef TestingOnly
        ${WriteEnvStr} ${VAR} $8
      !endif
    ${Else}
      ${LogText} "Setting enviroment variable ${VAR}='$8' for Current User Only."
      !undef ALL_USERS
      !ifndef TestingOnly
        ${WriteEnvStr} ${VAR} $8
      !endif
    ${Endif}
  ${EndIF}
!macroend

Function envarsPage

  !insertmacro MUI_HEADER_TEXT "Setting Enviroment Variables" "Choose whether to set or not the following enviroment variables"
  
  ReadEnvStr $1 "CGALROOT"   # $1 = existing value for CGALROOT
  
  ${If} $1 != ""
      !insertmacro UncheckEnvStrCheckbox 6
      !insertmacro SetEnvStrValueSlot    7 $1
      !insertmacro SetEnvStrLabel        8 $(DESC_EVARSET)
 ${Else}
      !insertmacro SetEnvStrValueSlot 7 $INSTDIR
  ${Endif}
  
  ReadEnvStr $2 "BOOSTROOT"  # $2 = existing value for BOOSTROOT
  ${If} $2 != ""
      !insertmacro UncheckEnvStrCheckbox 10
      !insertmacro SetEnvStrValueSlot    11 $2
      !insertmacro SetEnvStrLabel        12 $(DESC_EVARSET)
  ${Else}
      !insertmacro SetEnvStrValueSlot 11 "guess of boost"
  ${Endif}
  
  !insertmacro MUI_INSTALLOPTIONS_INITDIALOG "enviroment_variables.ini"
  Pop $0

  !insertmacro MUI_INSTALLOPTIONS_SHOW
  Pop $0

  # PROCESSING - Installs selected enviroment variables
  
  !insertmacro MUI_INSTALLOPTIONS_READ $3 "enviroment_variables.ini" "Field 2" "State" # $3=Is ALL USERS selected
  
  !insertmacro ProcessEnvStrCB $3 "CGALROOT"  6
  !insertmacro ProcessEnvStrCB $3 "BOOSTROOT" 10
FunctionEnd
