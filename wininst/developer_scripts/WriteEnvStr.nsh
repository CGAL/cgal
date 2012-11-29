#
# Taken from http://nsis.sourceforge.net/Setting_Environment_Variables
# User handling modified by Fernando Cacciola
# Laurent Rineau added un.DeleteEnvStr and adapted it to handle user.
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

!define un.DeleteEnvStr "!insertmacro un.DeleteEnvStr" 
!macro un.DeleteEnvStr name all_users
  Push ${name}
  Push ${all_users}
  Call un.DeleteEnvStr
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
 
#
# un.DeleteEnvStr - Removes an environment variable
# Note: Win9x systems requires reboot
#
# Example:
#  Push "HOMEDIR"           # name
#  Push 1 (all users) or 0 (current user only) 
#  Call un.DeleteEnvStr
#
Function un.DeleteEnvStr
  Exch $1 ; $1 all users?
  Exch
  Exch $0 ; $0 now has the name of the variable
  Push $2
  Push $3
  Push $4
  Push $5
  Push $6
 
  Call un.IsNT
  Pop $2
  StrCmp $2 1 DeleteEnvStr_NT
    ; Not on NT
    StrCpy $2 $WINDIR 2
    FileOpen $2 "$2\autoexec.bat" r
    GetTempFileName $5
    FileOpen $3 $5 w
    StrCpy $0 "SET $0="
    SetRebootFlag true
 
    DeleteEnvStr_dosLoop:
      FileRead $2 $4
      StrLen $6 $0
      StrCpy $6 $4 $6
      StrCmp $6 $0 DeleteEnvStr_dosLoop
      StrCmp $6 "" DeleteEnvStr_dosLoopEnd
      FileWrite $3 $4
      Goto DeleteEnvStr_dosLoop
 
    DeleteEnvStr_dosLoopEnd:
      FileClose $3
      FileClose $2
      StrCpy $2 $WINDIR 2
      Delete "$2\autoexec.bat"
      CopyFiles /SILENT $5 "$2\autoexec.bat"
      Delete $5
      Goto DeleteEnvStr_done
 
  DeleteEnvStr_NT:
    StrCmp $1 1 DelAllUsers
    DeleteRegValue ${WriteEnvStr_RegKey_CurrentUser} $0
    Goto DelWritten
  DelAllUsers:
    DeleteRegValue ${WriteEnvStr_RegKey_AllUsers} $0
  DelWritten:
    SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} \
      0 "STR:Environment" /TIMEOUT=5000
 
  DeleteEnvStr_done:
    Pop $6
    Pop $5
    Pop $4
    Pop $3
    Pop $2
    Pop $1
    Pop $0
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
!insertmacro IsNT "un."
 
!endif ; IsNT_KiCHiK
 
!endif ; _WriteEnvStr_nsh
