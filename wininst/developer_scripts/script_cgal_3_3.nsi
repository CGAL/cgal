;===========================
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


!include "MUI.nsh"
!include "WriteEnvStr.nsh"
!include "StrFunc.nsh"
!include "Sections.nsh"
!include "LogicLib.nsh"
!include "Locate.nsh"
!include "StrRep.nsh"
!include "ReplaceInFile.nsh"


!include "script_cgal_3_3.nsh"

!ifdef DebugLog
!include "TextLog.nsh"
!EndIf

!define CGAL_SRC  "CGAL-3.3"
;!define FTP_SRC   "http://www.geometryfactory.com/precompiled_libs/"
!define FTP_SRC   "ftp://ftp.mpi-sb.mpg.de/pub/outgoing/CGAL/precompiled_libs/"

;--------------------------------
; General
;--------------------------------

  ;Name and file
  Name "GAL-3.3"
  
  !ifdef FetchLocal
  OutFile "CGAL-3.3-Full-Setup.exe"
  !else
  OutFile "CGAL-3.3-Setup.exe"
  !endif

  ;Default installation folder
  InstallDir "$PROGRAMFILES\CGAL-3.3"

  ;Get installation folder from registry if available
  InstallDirRegKey HKCU "Software\CGAL-3.3" ""
  
  BrandingText "The CGAL Project and GeometryFactory - Installer created with NSIS."

  VIProductVersion "3.3.0.0"
  VIAddVersionKey "ProductName"     "CGAL Windows Installer"
  VIAddVersionKey "CompanyName"     "The CGAL Project and GeometryFactory"
  VIAddVersionKey "LegalCopyright"  "© The CGAL Project and GeometryFactory"
  VIAddVersionKey "FileDescription" "Windows Installer for CGAL"
  VIAddVersionKey "FileVersion"     "3.3"
  
;--------------------------------
; Variables
;--------------------------------

  Var MUI_TEMP
  Var STARTMENU_FOLDER
  Var DoFixupProjectFiles
  Var SetCGALROOT
  Var SetBOOSTROOT
  Var SetEnvAllUsers
  Var IsGmpInstalled
  
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

  !define MUI_WELCOMEPAGE_TEXT "This installs CGAL-3.3 on your machine, precompiled with .Net 2003 and 2005 (VC 7.1 and 8.0).\r\nThe project files for building the library itself are provided for both VC7.1 and VC8.0 separatedly, but for examples and demos, only the VC7.1 project files are provided since all its settings are compatible with VC8.0 (thus you can just open them with Visual Studio 2005 and follow the conversion wizard)."


  !define MUI_FINISHPAGE_TITLE "Installation Not Finished Yet!!!"

  !define MUI_FINISHPAGE_TEXT "CGAL needs part of the Boost Library which you must download yourself (from www.boost.org or www.boost-consulting.com/download.html).\r\nMost demos need the Qt 3 library which you must also download yourself (from www.trolltech.com).\r\nThe vcproj files in the example, demo amd src subdirectories all use the environment variables CGALROOT and BOOSTROOT to locate CGAL and Boost resp., but for other third-party libraries like Qt3 and ZLib, you must add the corresponding include/lib paths manually."

  !define MUI_FINISHPAGE_SHOWREADME "$INSTDIR\INSTALL.win32.txt"
  !define MUI_FINISHPAGE_SHOWREADME_TEXT "Read the full installation notes for further instructions"
  
  !define MUI_FINISHPAGE_LINK "More information about CGAL and Visual C++"
  !define MUI_FINISHPAGE_LINK_LOCATION http://www.cgal.org/platforms_frame.html
  
;--------------------------------
; Pages
;--------------------------------

  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_LICENSE "${CGAL_SRC}\LICENSE"
  
  ; A page where the user can specify a default variant configuration (taken from the boost installer)
  Page custom defaultVariantsPage 
  
  !insertmacro MUI_PAGE_COMPONENTS
 
  !insertmacro MUI_PAGE_DIRECTORY
  
  ; A page where the user can check/uncheck the environment variables
  ; used to specify paths in vproj files to be added.
  Page custom envarsPage 
  
  ; A page where the user can decide not to remove CGAL_USE_GMP from project files.
  ; This page is only shown if the MPFR/GMP component has not been installed.
  Page custom fixupPage  

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

!ifndef SkipFiles
  SectionIn RO 
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
  SetOutPath "$INSTDIR\lib"
  File /nonfatal "${CGAL_SRC}\lib\README"

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

!ifndef SkipFiles
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
${MultiVariantSection} "CGAL precompiled libs"          Install_CGAL_bin     CGAL_LIB_Idx
${MultiVariantSection} "GMP and MPFR precompiled libs"  Install_GMP_MPFR_bin GMP_LIB_Idx
;--------------------------------

;--------------------------------
Section /o "LAPACK and TAUCS precompiled libs" TAUCS_LIB_Idx

!ifndef SkipFiles
  SetOutPath "$INSTDIR\auxiliary\taucs\include"
  File /r "${CGAL_SRC}\auxiliary\taucs\include\*.*"
!endif

!ifndef FetchLocal
    SetOutPath "$INSTDIR\auxiliary\taucs\lib"
    File "${CGAL_SRC}\auxiliary\taucs\lib\README"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "libatlas.lib.zip"           "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "libcblas.lib.zip"           "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "libf77blas.lib.zip"         "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "liblapack.lib.zip"          "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "libmetis-vc71-mt-s.lib.zip" "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "libmetis-vc71-mt.lib.zip"   "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "libmetis-vc71-s.lib.zip"    "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "libtaucs-vc71-mt-s.lib.zip" "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "libtaucs-vc71-mt.lib.zip"   "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "libtaucs-vc71-s.lib.zip"    "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "libtstatlas.lib.zip"        "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "vcf2c-vc71-mt-s.lib.zip"    "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "vcf2c-vc71-mt.lib.zip"      "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "vcf2c-vc71-s.lib.zip"       "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "vcf2c-vc71-mt-sgd.lib.zip"  "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "vcf2c-vc71-mt-gd.lib.zip"   "$INSTDIR\auxiliary\taucs\lib"
	!insertmacro DownloadFile "auxiliary/TAUCS-CGAL-3.3/" "vcf2c-vc71-sgd.lib.zip"     "$INSTDIR\auxiliary\taucs\lib"
!else
  !ifndef SkipFiles
    SetOutPath "$INSTDIR\auxiliary\taucs\lib"
    File "${CGAL_SRC}\auxiliary\taucs\lib\README"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\libatlas.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\libcblas.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\libf77blas.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\liblapack.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\libmetis-vc71-mt-s.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\libmetis-vc71-mt.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\libmetis-vc71-s.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\libtaucs-vc71-mt-s.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\libtaucs-vc71-mt.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\libtstatlas.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\vcf2c-vc71-mt-s.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\vcf2c-vc71-mt.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\vcf2c-vc71-s.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\vcf2c-vc71-mt-sgd.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\vcf2c-vc71-mt-gd.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\taucs\lib\vcf2c-vc71-sgd.lib"
  !endif  
!endif  
SectionEnd
;--------------------------------

Section "-unzip" UNZIP_Idx
  
  ${locate::Open} "$INSTDIR" "/D=0 /X=zip" $0
  ${If} $0 != 0
    ${Do}
  	  ${locate::Find} $0 $1 $2 $3 $4 $5 $6
      ${If} "$1" != ""
        ZipDLL::extractall $1 $2
		Pop $7
		${If} "$7" == "success"
          Delete $1
		${EndIf}
      ${EndIf}
    ${LoopUntil} "$1" == ""
  ${EndIf}  
  ${locate::Close} $0
  ${locate::Unload}
  
SectionEnd

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
  LangString DESC_MAIN      ${LANG_ENGLISH} "The main components of the CGAL Library."
  LangString DESC_DEMOS     ${LANG_ENGLISH} "The CGAL demos, for which you need Qt 3 in order to build them."
  LangString DESC_CGAL_LIB  ${LANG_ENGLISH} "The precompiled CGAL libraries."
  LangString DESC_GMP_LIB   ${LANG_ENGLISH} "The precompiled GMP and MPFR libraries (needed for exact constructions)."
  LangString DESC_TAUCS_LIB ${LANG_ENGLISH} "The precompiled LAPACK and TAUCS libraries."
  LangString DESC_ENVSET    ${LANG_ENGLISH} "already set"

  ;Assign language strings to sections
  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
    !insertmacro MUI_DESCRIPTION_TEXT ${MAIN_Idx}      $(DESC_MAIN)
    !insertmacro MUI_DESCRIPTION_TEXT ${DEMOS_Idx}     $(DESC_DEMOS)
    !insertmacro MUI_DESCRIPTION_TEXT ${CGAL_LIB_Idx}  $(DESC_CGAL_LIB)
    !insertmacro MUI_DESCRIPTION_TEXT ${GMP_LIB_Idx}   $(DESC_GMP_LIB)
    !insertmacro MUI_DESCRIPTION_TEXT ${TAUCS_LIB_Idx} $(DESC_TAUCS_LIB)
  !insertmacro MUI_FUNCTION_DESCRIPTION_END


;--------------------------------
; Functions
;--------------------------------

Function .onInit

  !ifdef DebugLog
  ${LogSetFileName} "$INSTDIR\windows_install_log.txt"
  ${LogSetOn}
  !endif	

  # the plugins dir is automatically deleted when the installer exits
  InitPluginsDir
  File /oname=$PLUGINSDIR\splash.bmp ".\CGAL.bmp"
  advsplash::show 1000 600 400 -1 $PLUGINSDIR\splash
  
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "default_variants.ini"
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "environment_variables.ini"
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "fixup_projects.ini"
  
  !insertmacro SelectSection ${UNZIP_Idx}
  
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
  
  !insertmacro SelectSection ${UNZIP_Idx}
  
FunctionEnd

Function .onInstSuccess
  !insertmacro SelectSection ${UNZIP_Idx}
  
  ${If} "$SetBOOSTROOT" != ""
    !insertmacro SetEnvStr $SetEnvAllUsers "BOOSTROOT"  $SetBOOSTROOT
  ${EndIf}
  
  ${If} "$SetCGALROOT" != ""
    !insertmacro SetEnvStr $SetEnvAllUsers "CGALROOT"  $SetCGALROOT
  ${EndIf}
  
  ${If} "$DoFixupProjectFiles" == "y"
    Call FixupProjectFiles
  ${EndIf}
FunctionEnd

Function defaultVariantsPage

    !insertmacro MUI_HEADER_TEXT "Select Default Variants" "Choose the default variants to autoselect in the next page."
    !insertmacro MUI_INSTALLOPTIONS_INITDIALOG "default_variants.ini"
    !insertmacro MUI_INSTALLOPTIONS_SHOW

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
  !insertmacro MUI_INSTALLOPTIONS_WRITE "environment_variables.ini" "Field ${FN}" "State" "0"
!macroend

# Disables the env var checkbox # FN
!macro SetEnvStrValueSlot FN VAL
  !insertmacro MUI_INSTALLOPTIONS_WRITE "environment_variables.ini" "Field ${FN}" "State" "${VAL}"
!macroend

!macro SetEnvStrLabel FN VAL
  !insertmacro MUI_INSTALLOPTIONS_WRITE "environment_variables.ini" "Field ${FN}" "Text" "${VAL}"
!macroend


Function envarsPage

  Push $1
  Push $2
  Push $3
  Push $4
  
  !insertmacro MUI_HEADER_TEXT "Setting Environment Variables" "Choose whether to set or not the following environment variables"
  
  ReadEnvStr $1 "CGALROOT"   # $1 = existing value for CGALROOT
  
  !insertmacro SetEnvStrValueSlot 7 $INSTDIR
  !insertmacro SetEnvStrValueSlot 7 $INSTDIR
 
  ${If} $1 != ""
      StrCpy $3 "($(DESC_ENVSET):  $1 )"    
      !insertmacro UncheckEnvStrCheckbox 6
      !insertmacro SetEnvStrLabel        8 $3
  ${Endif}
  
  ReadEnvStr $2 "BOOSTROOT"  # $2 = existing value for BOOSTROOT
  ${If} $2 != ""
      StrCpy $4 "($(DESC_ENVSET):  $2 )"    
      !insertmacro UncheckEnvStrCheckbox 10
      !insertmacro SetEnvStrValueSlot    11 $2
      !insertmacro SetEnvStrLabel        12 $4
  ${Else}
      Call FindBoostFolder
	  ${If} "$FoundBoostFolder" == ""
	    StrCpy $FoundBoostFolder "(insert here the root folder of the boost libraries)"
	  ${EndIf}
      !insertmacro SetEnvStrValueSlot 11 $FoundBoostFolder
  ${Endif}
  
  !insertmacro MUI_INSTALLOPTIONS_INITDIALOG "environment_variables.ini"

  !insertmacro MUI_INSTALLOPTIONS_SHOW_RETURN
  Pop $0
  ${If} "$0" = "success"
    # PROCESSING - Installs selected environment variables
    
    !insertmacro MUI_INSTALLOPTIONS_READ $SetEnvAllUsers "environment_variables.ini" "Field 2" "State" # $3=Is ALL USERS selected
    
    !insertmacro MUI_INSTALLOPTIONS_READ $3 "environment_variables.ini" "Field 6" "State" # CGALROOT checkbox
    !insertmacro MUI_INSTALLOPTIONS_READ $4 "environment_variables.ini" "Field 7" "State" # CGALROOT value
    ${If} $3 = 1 
      StrCpy $SetCGALROOT $4
    ${EndIF}

    !insertmacro MUI_INSTALLOPTIONS_READ $3 "environment_variables.ini" "Field 10" "State" # BOOSTROOT checkbox
    !insertmacro MUI_INSTALLOPTIONS_READ $4 "environment_variables.ini" "Field 11" "State" # BOOSTROOT value 
    ${If} $3 = 1 
      StrCpy $SetBOOSTROOT $4
    ${EndIF}
  ${EndIf}
  
  Pop $4
  Pop $3
  Pop $2
  Pop $1
  
FunctionEnd

Function fixupPage

  !insertmacro SelectSection ${UNZIP_Idx}

  ${Unless}    ${SectionIsSelected}          ${GMP_LIB_Idx}
  ${AndUnless} ${SectionIsPartiallySelected} ${GMP_LIB_Idx}
    !insertmacro MUI_HEADER_TEXT "Customizing Project Files" "Customize the installed project files"
    
    !insertmacro MUI_INSTALLOPTIONS_WRITE "fixup_projects.ini" "Field 2" "State" "1"
    !insertmacro MUI_INSTALLOPTIONS_WRITE "fixup_projects.ini" "Field 3" "Text"  "(because MPFR/GMP has not been installed)"

    !insertmacro MUI_INSTALLOPTIONS_INITDIALOG "fixup_projects.ini"

    !insertmacro MUI_INSTALLOPTIONS_SHOW_RETURN
    Pop $0
    ${If} "$0" = "success"
    
      !insertmacro MUI_INSTALLOPTIONS_READ $1 "fixup_projects.ini" "Field 2" "State" # $1=Remove CGAL_USE_GMP flags
      
      ${If} $1 == 1
        StrCpy $DoFixupProjectFiles "y"
      ${EndIf}
      
    ${EndIf}
  ${EndUnless}
FunctionEnd
