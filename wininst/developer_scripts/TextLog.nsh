# TextLog.nsh v1.1 - 2005-12-26
# Written by Mike Schinkel [http://www.mikeschinkel.com/blog/]
 
Var /GLOBAL __TextLog_FileHandle
Var /GLOBAL __TextLog_FileName
Var /GLOBAL __TextLog_State
 
!define LogMsg '!insertmacro LogMsgCall'
!macro LogMsgCall _text
        Call LogSetOn
	Push "${_text}"
	Call LogText
	Call LogSetOff
!macroend
 
 
!define LogText '!insertmacro LogTextCall'
!macro LogTextCall _text
	Push "${_text}"
	Call LogText
!macroend
 
Function LogText
	Exch $0   	; pABC -> 0ABC
	FileWrite $__TextLog_FileHandle "$0$\r$\n"
	Pop $0		; 0ABC -> ABC
FunctionEnd
 
!define LogSetFileName '!insertmacro LogSetFileNameCall'
!macro LogSetFileNameCall _filename
	Push "${_filename}"
	Call LogSetFileName
!macroend
 
Function LogSetFileName
	Exch $0   	; pABC -> 0ABC
	StrCpy $__TextLog_FileName "$0"
	StrCmp $__TextLog_State "open" +1 +3
	Call LogSetOff
	Call LogSetOn
	Pop $0		; 0ABC -> ABC
FunctionEnd
 
!define LogSetOn '!insertmacro LogSetOnCall'
!macro LogSetOnCall
	Call LogSetOn
!macroend
 
Function LogSetOn
	StrCmp $__TextLog_FileName "" +1 AlreadySet
	StrCpy $__TextLog_FileName "$INSTDIR\install.log"
AlreadySet:
	StrCmp $__TextLog_State "open" +2
	FileOpen $__TextLog_FileHandle  "$__TextLog_FileName"  a
        FileSeek $__TextLog_FileHandle 0 END
	StrCpy $__TextLog_State "open"
FunctionEnd
 
!define LogSetOff '!insertmacro LogSetOffCall'
!macro LogSetOffCall
 	Call LogSetOff
!macroend
 
Function LogSetOff
	StrCmp $__TextLog_State "open" +1 +2
	FileClose $__TextLog_FileHandle
	StrCpy $__TextLog_State ""
FunctionEnd