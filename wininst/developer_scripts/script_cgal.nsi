;===========================
; Copyright 2007, 2008 GeometryFactory (France)
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
!include "Sections.nsh"
!include "LogicLib.nsh"
!include "Locate.nsh"
!include "WriteEnvStr.nsh"
!include "EnvVarUpdate.nsh"
!include "x64.nsh"

;!define DebugLog

!ifdef DebugLog
  !include "TextLog.nsh"
!endif

!include "script_cgal.nsh"

!define CGAL_SRC  "CGAL-4.8"
!define FTP_SRC   "https://cgal.geometryfactory.com/CGAL/precompiled_libs/"

;--------------------------------
; General
;--------------------------------

  ;Name and file
  Name "${CGAL_SRC}"
  
  !ifdef FetchLocal
  OutFile "${CGAL_SRC}-Full-Setup.exe"
  !else
  OutFile "${CGAL_SRC}-Setup.exe"
  !endif

  ;Default installation folder: C:\dev\CGAL-4.8
  ; See also .onInit
  Installdir ""


  ;Get installation folder from registry if available
  InstallDirRegKey HKCU "Software\${CGAL_SRC}" ""
  
  BrandingText "The CGAL Project and GeometryFactory - Installer created with NSIS."

  VIProductVersion "4.8.0.0"
  VIAddVersionKey "ProductName"     "CGAL Windows Installer"
  VIAddVersionKey "CompanyName"     "The CGAL Project and GeometryFactory"
  VIAddVersionKey "LegalCopyright"  "© The CGAL Project and GeometryFactory"
  VIAddVersionKey "FileDescription" "Windows Installer for CGAL"
  VIAddVersionKey "FileVersion"     "4.8"
  
;--------------------------------
; Variables
;--------------------------------

  Var SetCGAL_DIR
  Var RegLoc
  Var Add_GMP_LIB_DIR_to_PATH  
  
;--------------------------------
; Interface Settings
;--------------------------------

  !define MUI_ICON  "cgal.ico"
  !define MUI_UNICON  "cgal.ico"
  !define MUI_HEADERIMAGE
  !define MUI_HEADERIMAGE_BITMAP_NOSTRETCH
  !define MUI_HEADERIMAGE_UNBITMAP_NOSTRETCH

  !define MUI_HEADERIMAGE_BITMAP "cgal_very_small_FFFFFF.bmp" ; optional
  !define MUI_HEADERIMAGE_UNBITMAP "cgal_very_small_FFFFFF.bmp" ; optional
 
  !define MUI_FINISHPAGE_NOAUTOCLOSE
  
  !define MUI_ABORTWARNING

  !define MUI_WELCOMEFINISHPAGE_BITMAP Zirkel.bmp 
  !define MUI_WELCOMEFINISHPAGE_BITMAP_NOSTRETCH
  !define MUI_UNWELCOMEFINISHPAGE_BITMAP Zirkel.bmp 
  !define MUI_UNWELCOMEFINISHPAGE_BITMAP_NOSTRETCH

  !define MUI_COMPONENTSPAGE_SMALLDESC

  !define MUI_WELCOMEPAGE_TEXT "This downloads ${CGAL_SRC} to your machine."

  !define MUI_FINISHPAGE_TITLE "Downloading finished"

  !define MUI_FINISHPAGE_TEXT "You have downloaded CGAL successfully. Please continue the installation, reading the installation instructions."

  !define MUI_FINISHPAGE_LINK "Installation instructions"
  
  !define MUI_FINISHPAGE_LINK_LOCATION "http://www.cgal.org/download/windows.html"
  
;--------------------------------
; Pages
;--------------------------------

  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_LICENSE "${CGAL_SRC}\LICENSE"
     
  !insertmacro MUI_PAGE_COMPONENTS
 
  ; A page where the user can specify a default variant configuration (taken from the boost installer)
  Page custom VariantsPage

  !insertmacro MUI_PAGE_DIRECTORY
  
  ; A page where the user can check/uncheck the environment variables
  ; used to specify paths in vcproj files to be added.
  Page custom envarsPage 
  
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
  SetOutPath "$INSTDIR\auxiliary"
  File /nonfatal /r "${CGAL_SRC}\auxiliary\*.*"
  SetOutPath "$INSTDIR\cmake"
  File /r "${CGAL_SRC}\cmake\*.*"
  SetOutPath "$INSTDIR\config"
  File /r "${CGAL_SRC}\config\*.*"
  SetOutPath "$INSTDIR\doc_html"
  File /r "${CGAL_SRC}\doc_html\*.*"
  SetOutPath "$INSTDIR\include"
  File /r "${CGAL_SRC}\include\*.*"
  SetOutPath "$INSTDIR\scripts"
  File /r "${CGAL_SRC}\scripts\*.*"
  SetOutPath "$INSTDIR\src"
  File /r "${CGAL_SRC}\src\*.*"
  SetOutPath "$INSTDIR\demo\icons"
  File /r "${CGAL_SRC}\demo\icons\*.*"
  SetOutPath "$INSTDIR\demo\resources"
  File /r "${CGAL_SRC}\demo\resources\*.*"

  SetOutPath "$INSTDIR"
  File "${CGAL_SRC}\AUTHORS"
  File "${CGAL_SRC}\CHANGES"
  File "${CGAL_SRC}\CMakeLists.txt"
  File "${CGAL_SRC}\INSTALL.md"
  File "${CGAL_SRC}\LICENSE"
  File "${CGAL_SRC}\LICENSE.FREE_USE"
  File "${CGAL_SRC}\LICENSE.LGPL"
  File "${CGAL_SRC}\LICENSE.GPL"
  File "${CGAL_SRC}\VERSION"
  File ".\cgal.ico"
!endif

  ; Write uninstall informations
  ;   http://nsis.sourceforge.net/Add_uninstall_information_to_Add/Remove_Programs
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "DisplayName" "${CGAL_SRC} -- Computational Geometry Algorithms Library, version 4.8"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "UninstallString" "$\"$INSTDIR\Uninstall.exe$\""
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "QuietUninstallString" "$\"$INSTDIR\Uninstall.exe$\" /S"

  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "NoRepair" 1
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "InstallLocation" "$\"$INSTDIR$\""
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "DisplayIcon" "$\"$INSTDIR\cgal.ico$\""
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "Publisher" "The CGAL Project and GeometryFactory"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "URLInfoAbout" "http://www.cgal.org/"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "DisplayedVersion" "4.8.0"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
                   "CGALUninstallRegLoc" "$RegLoc"

  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"
SectionEnd
;--------------------------------

;--------------------------------
Section "CGAL Examples and Demos" SAMPLES_Idx

!ifndef SkipFiles
  SetOutPath "$INSTDIR\examples"
  File /r "${CGAL_SRC}\examples\*.*"
  SetOutPath "$INSTDIR\demo"
  File /r "${CGAL_SRC}\demo\*.*"
!endif
SectionEnd

; Download and install GMP and MPFR binaries.
; Depend only on the platform (one variant per platform)
Section "GMP and MPFR precompiled libs"  GMP_LIB_Idx
  !ifndef FetchLocal
    !insertmacro Install_GMP_MPFR_bin "$Platform"
  !endif
SectionEnd


;--------------------------------


Section /o "HTML Manuals" DOC_Idx
  !ifndef FetchLocal
    !insertmacro DownloadFileFrom "https://cgal.geometryfactory.com/" "CGAL/4.8/Manual/" "cgal_manual.zip"  "$INSTDIR\doc_html"
  !endif  
SectionEnd

Section "-Unzip"

  ${locate::Open} "$INSTDIR" "/D=0 /X=zip" $0
  ${If} $0 <> 0
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

  ReadRegStr $RegLoc HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}" \
    "CGALUninstallRegLoc"

  DeleteRegKey /ifempty HKCU "Software\${CGAL_SRC}"
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${CGAL_SRC}"

  ${un.EnvVarUpdate} $0 "PATH" "R" $RegLoc "$INSTDIR\auxiliary\gmp\lib"
  ${If} $RegLoc == HKLM
    ${un.DeleteEnvStr} "CGAL_DIR" 1
  ${Else}
    ${un.DeleteEnvStr} "CGAL_DIR" 0
  ${EndIf}

SectionEnd




;--------------------------------
;Descriptions

  ;Language strings
  LangString DESC_MAIN      ${LANG_ENGLISH} "The main components of the CGAL Library."
  LangString DESC_SAMPLES   ${LANG_ENGLISH} "The CGAL demos and examples, for which you need Qt 5 in order to build them (and Qt 3 for some)."
  LangString DESC_GMP_LIB   ${LANG_ENGLISH} "The precompiled GMP and MPFR libraries (needed for exact constructions)."
  LangString DESC_DOC       ${LANG_ENGLISH} "The HTML manuals."
  LangString DESC_ENVSET    ${LANG_ENGLISH} "already set"

  ;Assign language strings to sections
  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
    !insertmacro MUI_DESCRIPTION_TEXT ${MAIN_Idx}      $(DESC_MAIN)
    !insertmacro MUI_DESCRIPTION_TEXT ${SAMPLES_Idx}   $(DESC_SAMPLES)
    !insertmacro MUI_DESCRIPTION_TEXT ${GMP_LIB_Idx}   $(DESC_GMP_LIB)
    !insertmacro MUI_DESCRIPTION_TEXT ${DOC_Idx}       $(DESC_DOC)
  !insertmacro MUI_FUNCTION_DESCRIPTION_END


;--------------------------------
; Functions
;--------------------------------

Function .onInit

  # Setup the default installation dir
  ${If} $InstDir == "" ; /D= was not used on the command line
      StrCpy $InstDir "C:\dev\${CGAL_SRC}"
  ${EndIf}

  !ifdef DebugLog
  ${LogSetFileName} "${CGAL_SRC}_install_log.txt"
  ${LogSetOn}
  !endif	
  
  StrCpy $Platform "win32"

  # the plugins dir is automatically deleted when the installer exits
  InitPluginsDir
  File /oname=$PLUGINSDIR\splash.bmp ".\CGAL.bmp" 
  advsplash::show 1000 600 400 -1 $PLUGINSDIR\splash
  
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "variants.ini"
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "environment_variables.ini"
 
FunctionEnd


Function .onInstSuccess

  ${If} $SetCGAL_DIR != ""
    ; RegLoc can be either HKLM (all users) or HKCU (current user).
    ${If} $RegLoc == HKLM
      ${WriteEnvStr} "CGAL_DIR"  $SetCGAL_DIR 1
    ${Else}
      ${WriteEnvStr} "CGAL_DIR"  $SetCGAL_DIR 0
    ${Endif}
  ${EndIf}
  
  ${If} $Add_GMP_LIB_DIR_to_PATH = 1
    ; Append "$INSTDIR\auxiliary\gmp\lib" to the PATH.
    ; RegLoc can be either HKLM (all users) or HKCU (current user).
    ; The return value goes to $0
    ${EnvVarUpdate} $0 "PATH" "A" $RegLoc "$INSTDIR\auxiliary\gmp\lib"
  ${EndIf}
  
FunctionEnd

Function VariantsPage

    !insertmacro MUI_HEADER_TEXT "Select platform" "Choose the platform for precompiled libraries."
    !insertmacro MUI_INSTALLOPTIONS_INITDIALOG "variants.ini"
    !insertmacro MUI_INSTALLOPTIONS_SHOW

  !insertmacro MUI_INSTALLOPTIONS_READ $0 "variants.ini" "Field 1" "State"
  ${If} $0 = 1
    StrCpy $Platform "win32"
  ${Else}
    StrCpy $Platform "x64"
  ${Endif}

FunctionEnd

# Disables the env var checkbox # FN and textbox # FN+1
!macro UncheckEnvStrCheckbox FN
  !insertmacro MUI_INSTALLOPTIONS_WRITE "environment_variables.ini" "Field ${FN}" "State" "0"
!macroend

!macro CheckEnvStrCheckbox FN
  !insertmacro MUI_INSTALLOPTIONS_WRITE "environment_variables.ini" "Field ${FN}" "State" "1"
!macroend

!macro DisableEnvStrCheckbox FN
  !insertmacro MUI_INSTALLOPTIONS_WRITE "environment_variables.ini" "Field ${FN}" "Flags" "DISABLED"
!macroend

!macro EnableEnvStrCheckbox FN
  !insertmacro MUI_INSTALLOPTIONS_WRITE "environment_variables.ini" "Field ${FN}" "Flags" ""
!macroend

# Disables the env var checkbox # FN
!macro SetEnvStrValueSlot FN VAL
  !insertmacro MUI_INSTALLOPTIONS_WRITE "environment_variables.ini" "Field ${FN}" "State" "${VAL}"
!macroend

!macro SetEnvStrLabel FN VAL
  !insertmacro MUI_INSTALLOPTIONS_WRITE "environment_variables.ini" "Field ${FN}" "Text" "${VAL}"
!macroend


Function envarsPage

  Push $0
  Push $1
  Push $2
  Push $3
  Push $4
  Push $5

  
  !insertmacro MUI_HEADER_TEXT "Setting Environment Variables" "Choose whether to set or not the following environment variables"
  
  ReadEnvStr $1 "CGAL_DIR"   # $1 = existing value for CGAL_DIR
  
  !insertmacro SetEnvStrValueSlot 7 $INSTDIR
  !insertmacro SetEnvStrValueSlot 7 $INSTDIR
 
  ${If} $1 != ""
      StrCpy $3 "($(DESC_ENVSET):  $1 )"    
      !insertmacro UncheckEnvStrCheckbox 6
      !insertmacro SetEnvStrLabel        8 $3
  ${Endif}

  SectionGetText ${GMP_LIB_Idx} $2 
  
  SectionGetFlags ${GMP_LIB_Idx} $1 
  IntOp $2 $1 & ${SF_SELECTED}
  
  ${If} $2 == 0
    !insertmacro UncheckEnvStrCheckbox 9
    !insertmacro DisableEnvStrCheckbox 9
  ${Else}
    !insertmacro CheckEnvStrCheckbox  9
    !insertmacro EnableEnvStrCheckbox 9
  ${Endif}
  
  !insertmacro MUI_INSTALLOPTIONS_INITDIALOG "environment_variables.ini"

  !insertmacro MUI_INSTALLOPTIONS_SHOW_RETURN
  Pop $0
  ${If} "$0" = "success"
    # PROCESSING - Installs selected environment variables
    
    !insertmacro MUI_INSTALLOPTIONS_READ $3 "environment_variables.ini" "Field 2" "State" # $3=Is ALL USERS selected

    ${If} $3 == 1
      StrCpy $RegLoc "HKLM"
    ${Else}  
      StrCpy $RegLoc "HKCU"
    ${EndIf}    
    
    !insertmacro MUI_INSTALLOPTIONS_READ $3 "environment_variables.ini" "Field 6" "State" # CGAL_DIR checkbox
    !insertmacro MUI_INSTALLOPTIONS_READ $4 "environment_variables.ini" "Field 7" "State" # CGAL_DIR value
    ${If} $3 == 1 
      StrCpy $SetCGAL_DIR $4
    ${EndIF}

    !insertmacro MUI_INSTALLOPTIONS_READ $5 "environment_variables.ini" "Field 9" "State" # Add to PATH checkbox
    ${If} $5 == 1 
      StrCpy $Add_GMP_LIB_DIR_to_PATH 1
    ${EndIF}
    
  ${EndIf}
  
  Pop $5
  Pop $4
  Pop $3
  Pop $2
  Pop $1
  Pop $0
  
FunctionEnd
