#
# Taken from http://nsis.sourceforge.net/Setting_Environment_Variables
# User handling modified by Fernando Cacciola
#
!ifndef _WriteEnvStr_nsh
!define _WriteEnvStr_nsh
 
#
#  Macro definition added by Fernando Cacciola
#
!define WriteEnvStr "!insertmacro WriteEnvStr" 
!macro WriteEnvStr name value all_users
  Push ${name}
  Push ${value}
  Push ${all_users}
  Call WriteEnvStr
!macroend

!include WinMessages.nsh
 
!define WriteEnvStr_RegKey_AllUsers 'HKLM "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"'
	   
!define WriteEnvStr_RegKey_CurrentUser 'HKCU "Environment"'
 
#
# WriteEnvStr - Writes an environment variable
# Note: Win9x systems requires reboot
#
# Example:
#  Push "HOMEDIR"           # name
#  Push "C:\New Home Dir\"  # value
#  Push 1 (all users) or 0 (current user only) 
#  Call WriteEnvStr
#
Function WriteEnvStr
  Exch $2 ; $2 all users?
  Exch
  Exch $1 ; $1 has environment variable value
  Exch 2
  Exch $0 ; $0 has environment variable name
  Push $3
 
  Call IsNT
  Pop $3
  StrCmp $3 1 WriteEnvStr_NT
    ; Not on NT
    StrCpy $3 $WINDIR 2 ; Copy drive of windows (c:)
    FileOpen $3 "$3\autoexec.bat" a
    FileSeek $3 0 END
    FileWrite $3 "$\r$\nSET $0=$1$\r$\n"
    FileClose $3
    SetRebootFlag true
    Goto WriteEnvStr_done
 
  WriteEnvStr_NT:
    StrCmp $2 1 AllUsers
      WriteRegExpandStr ${WriteEnvStr_RegKey_CurrentUser} $0 $1
	  Goto Written
	AllUsers:
      WriteRegExpandStr ${WriteEnvStr_RegKey_AllUsers} $0 $1
	Written:  
      SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000
 
  WriteEnvStr_done:
    Pop $3
    Pop $2
    Pop $0
    Pop $1
FunctionEnd
 
!ifndef IsNT_KiCHiK
!define IsNT_KiCHiK
 
#
# [un.]IsNT - Pushes 1 if running on NT, 0 if not
#
# Example:
#   Call IsNT
#   Pop $0
#   StrCmp $0 1 +3
#     MessageBox MB_OK "Not running on NT!"
#     Goto +2
#     MessageBox MB_OK "Running on NT!"
#
!macro IsNT UN
Function ${UN}IsNT
  Push $0
  ReadRegStr $0 HKLM \
    "SOFTWARE\Microsoft\Windows NT\CurrentVersion" CurrentVersion
  StrCmp $0 "" 0 IsNT_yes
  ; we are not NT.
  Pop $0
  Push 0
  Return
 
  IsNT_yes:
    ; NT!!!
    Pop $0
    Push 1
FunctionEnd
!macroend
!insertmacro IsNT ""
 
!endif ; IsNT_KiCHiK
 
!endif ; _WriteEnvStr_nsh