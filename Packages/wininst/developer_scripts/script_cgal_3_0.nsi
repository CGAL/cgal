; The name of the installer
Name "CgalInstaller"

; The file to write
OutFile "CGAL3.0Setup.exe"

; The default installation directory
InstallDir $PROGRAMFILES\CGAL-3.0

!define MUI_PRODUCT "CGAL" ;Define your own software name here
!define MUI_VERSION "3.0" ;Define your own software version here
!define MUI_ICON "Icons\cgal.ico"
!define MUI_UNICON "Icons\cgal-uninstall2.ico"


!include "MUI.nsh"

;Remember the Start Menu Folder
!define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU" 
!define MUI_STARTMENUPAGE_REGISTRY_KEY "Software\${MUI_PRODUCT}{MUI_VERSION}" 
!define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Start Menu Folder"

!define TEMP $R0
XPStyle on

;--------------------------------
;Modern UI Configuration
  !define MUI_WELCOMEPAGE
  !define MUI_CUSTOMPAGECOMMANDS
  !define MUI_LICENSEPAGE
  !define MUI_COMPONENTSPAGE
  !define MUI_DIRECTORYPAGE
  !define MUI_STARTMENUPAGE
  !define MUI_ABORTWARNING
  !define MUI_UNINSTALLER
  !define MUI_UNCONFIRMPAGE
  ;Modern UI System
  ;!insertmacro MUI_SYSTEM


;--------------------------------
;Languages

  ;!insertmacro MUI_LANGUAGE "Bulgarian"
  !insertmacro MUI_LANGUAGE "Czech"
  !insertmacro MUI_LANGUAGE "SimpChinese"
  !insertmacro MUI_LANGUAGE "TradChinese"
  !insertmacro MUI_LANGUAGE "Danish"
  !insertmacro MUI_LANGUAGE "Dutch"
  !insertmacro MUI_LANGUAGE "English"
  !insertmacro MUI_LANGUAGE "French"
  !insertmacro MUI_LANGUAGE "German"
  !insertmacro MUI_LANGUAGE "Greek"
  !insertmacro MUI_LANGUAGE "Hungarian"
  !insertmacro MUI_LANGUAGE "Italian"
  !insertmacro MUI_LANGUAGE "Japanese"
  !insertmacro MUI_LANGUAGE "Macedonian"
  !insertmacro MUI_LANGUAGE "Polish"
  !insertmacro MUI_LANGUAGE "PortugueseBR"
  !insertmacro MUI_LANGUAGE "Romanian"
  !insertmacro MUI_LANGUAGE "Russian"
  !insertmacro MUI_LANGUAGE "Slovak"
  !insertmacro MUI_LANGUAGE "Spanish"
  !insertmacro MUI_LANGUAGE "Thai"
  !insertmacro MUI_LANGUAGE "Ukrainian"
  
;--------------------------------


;--------------------------------
;Data
  
  LicenseData /LANG=${LANG_ENGLISH} "License.txt"
  LicenseData /LANG=${LANG_FRENCH} "License.txt"
  LicenseData /LANG=${LANG_GERMAN} "License.txt"
  LicenseData /LANG=${LANG_SPANISH} "License.txt"
  LicenseData /LANG=${LANG_SIMPCHINESE} "License.txt"
  LicenseData /LANG=${LANG_TRADCHINESE} "License.txt"
  LicenseData /LANG=${LANG_JAPANESE} "License.txt"
  LicenseData /LANG=${LANG_ITALIAN} "License.txt"
  LicenseData /LANG=${LANG_DUTCH} "License.txt"
  LicenseData /LANG=${LANG_DANISH} "License.txt"
  LicenseData /LANG=${LANG_GREEK} "License.txt"
  LicenseData /LANG=${LANG_RUSSIAN} "License.txt"
  LicenseData /LANG=${LANG_PORTUGUESEBR} "License.txt"
  LicenseData /LANG=${LANG_POLISH} "License.txt"
  LicenseData /LANG=${LANG_UKRAINIAN} "License.txt"
  LicenseData /LANG=${LANG_CZECH} "License.txt"
  LicenseData /LANG=${LANG_SLOVAK} "License.txt"
  ;LicenseData /LANG=${LANG_BULGARIAN} "License.txt"
  LicenseData /LANG=${LANG_HUNGARIAN} "License.txt"
  LicenseData /LANG=${LANG_THAI} "License.txt"
  LicenseData /LANG=${LANG_ROMANIAN} "License.txt"
  LicenseData /LANG=${LANG_MACEDONIAN} "License.txt"




;--------------------------------
;Pages
  !insertmacro MUI_PAGECOMMAND_WELCOME
  !insertmacro MUI_PAGECOMMAND_LICENSE
  !insertmacro MUI_PAGECOMMAND_COMPONENTS
  !insertmacro MUI_PAGECOMMAND_DIRECTORY
  !insertmacro MUI_PAGECOMMAND_STARTMENU
  !insertmacro MUI_PAGECOMMAND_INSTFILES

;--------------------------------
;Installer Sections


Section "CGAL library" SecCopyUI
  SectionIn RO 
  SetOutPath "$INSTDIR\auxiliary"
  FILE /r "CGAL3.0\auxiliary\*.*"
  SetOutPath "$INSTDIR\config"
  FILE /r "CGAL3.0\config\*.*"
  SetOutPath "$INSTDIR\include"
  FILE /r "CGAL3.0\include\*.*"
  SetOutPath "$INSTDIR\lib"
  FILE /r "CGAL3.0\lib\*.*"
  SetOutPath "$INSTDIR\make"
  FILE /r "CGAL3.0\make\*.*"
  SetOutPath "$INSTDIR\src"
  FILE /r "CGAL3.0\src\*.*"
;  SetOutPath "$INSTDIR\winutils"
;  FILE /r "CGAL3.0\winutils\*.*"
  SetOutPath "$INSTDIR\examples"
  FILE /r "CGAL3.0\examples\*.*"
  SetOutPath "$INSTDIR\demo"
  FILE /r "CGAL3.0\demo\*.*"
  SetOutPath "$INSTDIR\scripts"
  FILE /r "CGAL3.0\scripts\*.*"
  SetOutPath "$INSTDIR\doc_html"
  FILE /r "CGAL3.0\doc_html\*.*"
  SetOutPath "$INSTDIR\doc_pdf"
  FILE /r "CGAL3.0\doc_pdf\*.*"
  SetOutPath "$INSTDIR\doc_ps"
  FILE /r "CGAL3.0\doc_ps\*.*"



  ;the next will not be part of the public release
;  SetOutPath "$INSTDIR\developer_scripts"
;  FILE /r "CGAL3.0\developer_scripts\*.*"
;  SetOutPath "$INSTDIR\doc_tex"
;  FILE /r "CGAL3.0\doc_tex\*.*"
;  SetOutPath "$INSTDIR\test"
;  FILE /r "CGAL3.0\test\*.*"

  SetOutPath "$INSTDIR"
  FILE "CGAL3.0\Changes"
  FILE "CGAL3.0\Install"
  FILE "CGAL3.0\INSTALL.win32"
  FILE "CGAL3.0\Install_win32.chm"
  FILE "CGAL3.0\install_cgal"
  FILE "CGAL3.0\License"
  FILE "CGAL3.0\LICENSE.LGPL"
  FILE "CGAL3.0\LICENSE.QPL"
  FILE "CGAL3.0\Readme"
  FILE "Icons\cgal.ico"
  ;Read a value from an InstallOptions INI File
  ;!insertmacro MUI_INSTALLOPTIONS_READ ${TEMP} "ioCGAL.ini" "Field 2" "State"
  ;StrCmp ${TEMP} "0" "" +2
    ;Checked
    ;MessageBox MB_OK "A MessageBox..."
 

  !insertmacro MUI_STARTMENU_WRITE_BEGIN    
    ;Create shortcuts
    CreateDirectory "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}"
    CreateDirectory "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\DOC"
    CreateShortCut "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\DOC\installation_pdf.lnk" "$INSTDIR\doc_pdf\installation.pdf"
    CreateShortCut "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\DOC\installation_ps.lnk" "$INSTDIR\doc_ps\installation.ps"
    CreateShortCut "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\DOC\contents.lnk" "$INSTDIR\doc_html\installation\contents.html"

    CreateShortCut "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\Install_win32.lnk" "$INSTDIR\Install_win32.chm" 
    CreateShortCut "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\README.lnk" "$INSTDIR\Readme"
    CreateShortCut "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  
    ;Write shortcut location to the registry (for Uninstaller)
    WriteRegStr HKCU "Software\${MUI_PRODUCT}${MUI_VERSION}" "Start Menu Folder" "${MUI_STARTMENUPAGE_VARIABLE}"
  !insertmacro MUI_STARTMENU_WRITE_END
  
;  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall" "${MUI_PRODUCT}${MUI_VERSION}"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${MUI_PRODUCT}${MUI_VERSION}" "DisplayName" "${MUI_PRODUCT}${MUI_VERSION}"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${MUI_PRODUCT}${MUI_VERSION}" "UninstallString" "$INSTDIR\Uninstall.exe"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${MUI_PRODUCT}${MUI_VERSION}" "DisplayIcon" "$INSTDIR\cgal.ico"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${MUI_PRODUCT}${MUI_VERSION}" "DisplayVersion" "${MUI_VERSION}"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${MUI_PRODUCT}${MUI_VERSION}" "Publisher" "Geometry Factory"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${MUI_PRODUCT}${MUI_VERSION}" "HelpLink" "mailto:contact@cgal.org"
;  WriteRegStr HKCU "Software\${MUI_PRODUCT}${MUI_VERSION}" "Installer Language" 

  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"   

SectionEnd




Section "Set CGALROOT" ROOTSec5
  Call GetWindowsVersion
  Pop $R0
  Push $R1
  StrCpy $R1 "You are running on Windows $R0.$\r$\nI configure CGALROOT to be $INSTDIR" 100 0
  MessageBox MB_OK $R1
  Pop $R1
  StrCmp $R0 "ME" lbl_me
  StrCmp $R0 "98" lbl_ninetyeight
  StrCmp $R0 "98" lbl_ninetyfive
  Goto lbl_ntlike

  lbl_me:
    Goto lbl_autoexec
  Goto lbl_done
  lbl_ninetyeight:
    Goto lbl_autoexec
  Goto lbl_done
  lbl_ninetyfive:
    Goto lbl_autoexec
  Goto lbl_done
  lbl_ntlike:
    WriteRegStr HKCU "Environment" "CGALROOT" $INSTDIR
    MessageBox MB_YESNO "Do you want to set it for every user on this computer?$\r$\nYou need administrator rights to do it." IDNO no
    WriteRegStr HKEY_LOCAL_MACHINE "SYSTEM\CurrentControlSet\Control\Session Manager\Environment" "CGALROOT" "$INSTDIR"
    MessageBox MB_OK "You must logout and login to have CGALROOT set in your environment.$\r$\nAlso you should copy $INSTDIR\auxiliary\wingmp\gmp-4.1.2\gmp.dll on your PATH. (eg. C:\Windows\System32)"
    no:
    MessageBox MB_OK "You must logout and login to have CGALROOT set in your environment.$\r$\nAlso you should copy $INSTDIR\auxiliary\wingmp\gmp-4.1.2\gmp.dll on your PATH. (eg. C:\Windows\System32)"

  Goto lbl_done
  lbl_autoexec:
    FileOpen $1 "$1\autoexec.bat" a
    FileSeek $1 0 END
    FileWrite $1 "$\r$\nSET CGALROOT=$INSTDIR$\r$\n"
    FileClose $1
    MessageBox MB_OK "You will have to reboot before the CGALROOT is set in your environment.$\r$\nAlso you should copy $INSTDIR\auxiliary\wingmp\gmp-4.1.2\gmp.dll on your PATH. (eg. C:\Windows\System32)"
    SetRebootFlag true
  Goto lbl_done


  lbl_done:
  	
SectionEnd



;Display the Finish header
;Insert this macro after the sections if you are not using a finish page
!insertmacro MUI_SECTIONS_FINISHHEADER




;--------------------------------
;Installer Functions

Function .onInit

  # the plugins dir is automatically deleted when the installer exits
  InitPluginsDir
  File /oname=$PLUGINSDIR\splash.bmp "LOGO\CGAL.bmp"
  advsplash::show 1000 600 400 -1 $PLUGINSDIR\splash

  ;Language selection

  ;Font
  Push Tahoma
  Push 8

  ;!insertmacro MUI_LANGDLL_PUSH "Bulgarian"

  Push ""
  Push ${LANG_CZECH}
  Push Czech
  Push ${LANG_SIMPCHINESE}
  Push SimpChinese
  Push ${LANG_TRADCHINESE}
  Push TradChinese
  Push ${LANG_DANISH}
  Push Danish
  Push ${LANG_GERMAN}
  Push German
  Push ${LANG_ENGLISH}
  Push English
  Push ${LANG_SPANISH}
  Push Spanish
  Push ${LANG_FRENCH}
  Push French
  Push ${LANG_GREEK}
  Push Greek
  Push ${LANG_ITALIAN}
  Push Italian  
  Push ${LANG_JAPANESE}
  Push Japanese
  Push ${LANG_MACEDONIAN}
  Push Macedonian
  Push ${LANG_HUNGARIAN}
  Push Hungarian
  Push ${LANG_DUTCH}
  Push Dutch
  Push ${LANG_POLISH}
  Push Polish
  Push ${LANG_PORTUGUESE}
  Push Portuguese
  Push ${LANG_ROMANIAN}
  Push Romanian
  Push ${LANG_RUSSIAN}
  Push Russian
  Push ${LANG_SLOVAK}
  Push Slovak
  Push ${LANG_THAI}
  Push Thai
  Push ${LANG_UKRAINIAN}
  Push Ukrainian

  Push A ; A means auto count languages
  ; for the auto count to work the first empty push (Push "") must remain

  LangDLL::LangDialog "Installer Language" "Please select a language."

  Pop $LANGUAGE
  StrCmp $LANGUAGE "cancel" 0 +2
    Abort


FunctionEnd

 
 ; GetWindowsVersion
 ;
 ; Based on Yazno's function, http://yazno.tripod.com/powerpimpit/
 ; Returns on top of stack
 ;
 ; Windows Version (95, 98, ME, NT x.x, 2000, XP, .NET Server)
 ; or
 ; '' (Unknown Windows Version)
 ;
 ; Usage:
 ;   Call GetWindowsVersion
 ;   Pop $R0
 ;   ; at this point $R0 is "NT 4.0" or whatnot

 Function GetWindowsVersion
   Push $R0
   Push $R1
   ReadRegStr $R0 HKLM "SOFTWARE\Microsoft\Windows NT\CurrentVersion" CurrentVersion
   StrCmp $R0 "" 0 lbl_winnt
   ; we are not NT.
   ;ReadRegStr $R0 HKLM SOFTWARE\Microsoft\Windows\CurrentVersion VersionNumber

   StrCpy $R1 $R0 1
   StrCmp $R1 '4' 0 lbl_error

   StrCpy $R1 $R0 3

   StrCmp $R1 '4.0' lbl_win32_95
   StrCmp $R1 '4.9' lbl_win32_ME lbl_win32_98

   lbl_win32_95:
     StrCpy $R0 '95'
   Goto lbl_done

   lbl_win32_98:
     StrCpy $R0 '98'
   Goto lbl_done

   lbl_win32_ME:
     StrCpy $R0 'ME'
   Goto lbl_done

   lbl_winnt:

     StrCpy $R1 $R0 1

     StrCmp $R1 '3' lbl_winnt_x
     StrCmp $R1 '4' lbl_winnt_x

     StrCpy $R1 $R0 3

     StrCmp $R1 '5.0' lbl_winnt_2000
     StrCmp $R1 '5.1' lbl_winnt_XP
     StrCmp $R1 '5.2' lbl_winnt_dotNET lbl_error

     lbl_winnt_x:
       StrCpy $R0 "NT $R0" 6
     Goto lbl_done

     lbl_winnt_2000:
       Strcpy $R0 '2000'
     Goto lbl_done

     lbl_winnt_XP:
       Strcpy $R0 'XP'
     Goto lbl_done

     lbl_winnt_dotNET:
       Strcpy $R0 '.NET Server'
     Goto lbl_done

   lbl_error:
     Strcpy $R0 ''
   lbl_done:
   Pop $R1
   Exch $R0
 FunctionEnd




; special uninstall section.
Section "Uninstall"
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${MUI_PRODUCT}${MUI_VERSION}"
  DeleteRegKey HKLM "SOFTWARE\${MUI_PRODUCT}${MUI_VERSION}"
  DeleteRegKey HKCU "SOFTWARE\${MUI_PRODUCT}${MUI_VERSION}"

  RMDir /r "$INSTDIR"

  Delete "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\DOC\Installation_ps.lnk"
  Delete "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\DOC\Installation_pdf.lnk"
  Delete "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\DOC\contents.lnk"
  Delete "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\Install_win32.lnk"  
  Delete "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\README.lnk"
  Delete "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\Uninstall.lnk"
  RMDir "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}\DOC" ;Only if empty, so it won't delete other 
  RMDir "$SMPROGRAMS\${MUI_PRODUCT}${MUI_VERSION}" ;Only if empty, so it won't delete other shortcuts

  IfFileExists "$INSTDIR" 0 NoErrorMsg
   MessageBox MB_OK "Note: $INSTDIR could not be removed!" IDOK 0 ; skipped if file doesn't exist
NoErrorMsg:
  ;Display the Finish header
  !insertmacro MUI_UNFINISHHEADER

SectionEnd

;--------------------------------
;Uninstaller Functions

Function un.onInit
  ;Get language from registry
  ReadRegStr $LANGUAGE HKCU "Software\${MUI_PRODUCT}" "Installer Language"  
FunctionEnd

