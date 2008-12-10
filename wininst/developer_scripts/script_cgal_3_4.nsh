;============================
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

${StrStr}

;!define SkipFiles
;!define SkipSetEnvVar
;!define SkipDownload
;!define FetchLocal
;!define DebugLog
!define ViaFTP

Var no_default_compilers
Var no_default_variants
Var selected_libs
Var FoundBoostFolder
Var Platform

;--------------------------------
; Macros
;--------------------------------


!define MultiVariantSection "!insertmacro MultiVariantSection"

; Expands to a Section Group named "SecName" which contains all library variants.
; For each variant, the macro "Handler" is expanded with the variant name as argument
!macro MultiVariantSection SecName Handler Platform Idx
  SectionGroup "${SecName}" ${Idx}
    SectionGroup "VC8.0"
      Section /o "Multithread Debug"
        !insertmacro "${Handler}" "${Platform}" "vc80-mt-gd"
      SectionEnd
      Section /o "Multithread"
        !insertmacro "${Handler}" "${Platform}" "vc80-mt"
      SectionEnd
      Section /o "Multithread, static runtime"
        !insertmacro "${Handler}" "${Platform}" "vc80-mt-s"
      SectionEnd
      Section /o "Multithread Debug, static runtime"
        !insertmacro "${Handler}" "${Platform}" "vc80-mt-sgd"
      SectionEnd
    SectionGroupEnd
    SectionGroup "VC9.0"
      Section /o "Multithread Debug"
        !insertmacro "${Handler}" "${Platform}" "vc90-mt-gd"
      SectionEnd
      Section /o "Multithread"
        !insertmacro "${Handler}" "${Platform}" "vc90-mt"
      SectionEnd
      Section /o "Multithread, static runtime"
        !insertmacro "${Handler}" "${Platform}" "vc90-mt-s"
      SectionEnd
      Section /o "Multithread Debug, static runtime"
        !insertmacro "${Handler}" "${Platform}" "vc90-mt-sgd"
      SectionEnd
    SectionGroupEnd
  SectionGroupEnd
!macroend

!ifdef ViaFTP
  !define DownloadOK      "OK"
  !define DownloadAborted "cancel"
!else
  !define DownloadOK      "success"
  !define DownloadAborted "cancel"
!endif

!macro DownloadFile SRC_FOLDER FILE TGT
!ifndef SkipDownload
  !ifdef ViaFTP
    inetc::get ${FTP_SRC}${SRC_FOLDER}${FILE} ${TGT}\${FILE}
  !else
    NSISdl::download ${FTP_SRC}${SRC_FOLDER}${FILE} ${TGT}\${FILE}
  !endif  
  Pop $0
  ${If}     "$0" == "OK"
    DetailPrint "${FILE} downloaded successfully."
  ${ElseIf} "$0" == "URL Parts Error"
    DetailPrint "${FILE} downloaded successfully."
  ${ElseIf} "$0" == "Terminated"
    DetailPrint "${FILE} download CANCELLED."
  ${ElseIf} "$0" == "Cancelled"
    DetailPrint "${FILE} download CANCELLED."
  ${Else}  
    MessageBox MB_OK "Unable to download ${FTP_SRC}${SRC_FOLDER}${FILE}. Error: $0"
    DetailPrint "ERROR $0: Unable to download ${FTP_SRC}${SRC_FOLDER}${FILE}."
  ${Endif}
!endif	
!macroend

!macro Install_PDB_if_debug_variant HANDLER PLATFORM VARIANT
  ${StrStr} $R0 ${VARIANT} "gd"
  ${If} "$R0" != ""
    !insertmacro "${HANDLER}" "${PLATFORM}" "${VARIANT}"
  ${EndIf}  
!macroend

!macro Install_DLL_if_dynamic_variant HANDLER PLATFORM VARIANT
  ${StrStr} $R0 ${VARIANT} "s"
  ${If} "$R0" == ""
    !insertmacro "${HANDLER}" "${PLATFORM}" "${VARIANT}"
  ${EndIf}  
!macroend

!macro Install_GMP_MPFR_libs PLATFORM VARIANT
!ifndef FetchLocal
    !insertmacro DownloadFile "auxiliary/${PLATFORM}/GMP/4.2.4/"  "gmp-${VARIANT}.lib.zip"  "$INSTDIR\auxiliary\gmp\lib"
    !insertmacro DownloadFile "auxiliary/${PLATFORM}/MPFR/2.3.2/" "mpfr-${VARIANT}.lib.zip" "$INSTDIR\auxiliary\gmp\lib"
!else
  !ifndef SkipFiles
    SetOutPath "$INSTDIR\auxiliary\gmp\lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\gmp-${VARIANT}.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\mpfr-${VARIANT}.lib"
  !endif  
!endif
!macroend

!macro Install_GMP_MPFR_dlls PLATFORM VARIANT
!ifndef FetchLocal
    !insertmacro DownloadFile "auxiliary/${PLATFORM}/GMP/4.2.4/"  "gmp-${VARIANT}.dll.zip"  "$INSTDIR\auxiliary\gmp\lib"
    !insertmacro DownloadFile "auxiliary/${PLATFORM}/MPFR/2.3.2/" "mpfr-${VARIANT}.dll.zip" "$INSTDIR\auxiliary\gmp\lib"
!else
  !ifndef SkipFiles
    SetOutPath "$INSTDIR\auxiliary\gmp\lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\gmp-${VARIANT}.dll"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\mpfr-${VARIANT}.dll"
  !endif  
!endif
!macroend

!macro Install_GMP_MPFR_pdbs PLATFORM VARIANT
!ifndef FetchLocal
    !insertmacro DownloadFile "auxiliary/${PLATFORM}/GMP/4.2.4/"  "gmp-${VARIANT}.pdb.zip"  "$INSTDIR\auxiliary\gmp\lib"
    !insertmacro DownloadFile "auxiliary/${PLATFORM}/MPFR/2.3.2/" "mpfr-${VARIANT}.pdb.zip" "$INSTDIR\auxiliary\gmp\lib"
!else
  !ifndef SkipFiles
    SetOutPath "$INSTDIR\auxiliary\gmp\lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\gmp-${VARIANT}.pdb"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\mpfr-${VARIANT}.pdb"
  !endif  
!endif
!macroend

!macro Install_GMP_MPFR_bin PLATFORM VARIANT

  ; Headers are not VARIANT dependent so we include this only once, but here since
  ; we want to download headers only if at least one lib variant was selected.
  ${If} $IsGmpInstalled = 0
    StrCpy $IsGmpInstalled 1 
    SetOutPath "$INSTDIR\auxiliary\gmp\lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\README"
    SetOutPath "$INSTDIR\auxiliary\gmp\include"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\include\README"
    !ifndef FetchLocal
        !insertmacro DownloadFile "auxiliary/${PLATFORM}/GMP/4.2.4/"  "gmp.h.zip"  "$INSTDIR\auxiliary\gmp\include"
        !insertmacro DownloadFile "auxiliary/${PLATFORM}/MPFR/2.3.2/" "mpfr.h.zip" "$INSTDIR\auxiliary\gmp\include"
    !else
      !ifndef SkipFiles
        SetOutPath "$INSTDIR\auxiliary\gmp\include"
        File /r "${CGAL_SRC}\auxiliary\gmp\include\*.*"
        SetOutPath "$INSTDIR\auxiliary\gmp\lib"
        File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\README"
      !endif
    !endif
  ${Endif}
  
  !insertmacro Install_GMP_MPFR_libs                                "${PLATFORM}" "${VARIANT}"
  !insertmacro Install_PDB_if_debug_variant Install_GMP_MPFR_pdbs   "${PLATFORM}" "${VARIANT}"
  !insertmacro Install_DLL_if_dynamic_variant Install_GMP_MPFR_dlls "${PLATFORM}" "${VARIANT}"
!macroend


!macro SetEnvStr ALLUSERS VAR VALUE
  # ${ALLUSERS} is 0 or 1
  # ${VAR}      is the env var to set
  # ${VALUE}    is the env var value
  ${WriteEnvStr} ${VAR} ${VALUE} ${ALLUSERS}
!macroend

;--------------------------------
; Functions
;--------------------------------

Function initSelectionFlags

    StrCpy $selected_libs ""
    ClearErrors
    StrCpy $0 0
    
    
  next:
    SectionGetText $0 $1
    IfErrors bail
    StrCpy $2 $1 "" -4
    StrCmp $2 "libs" 0 not_lib
    Push $0
    call SelectDefaultVariants
    StrCpy $selected_libs "$selected_libs1"
  not_lib:
    IntOp $0 $0 + 1
    goto next
  bail:
FunctionEnd

; Stack 0: compiler name
; Stack 1: variant name
; Stack 2: section

Function MaybeSelectVariant
    Exch $2
    ; c, v, r2
    Exch
    ; c, r2, v
    Exch $1
    ; c, r2, r1
    Exch
    ; c, r1, r2
    Exch 2
    ; r2, r1, c
    Exch $0
    ; r2, r1, r0
    Exch 2
    ; r0, r1, r2

    Push $3
    Push $4

        ${If} $0 == "VC8.0"
            !insertmacro MUI_INSTALLOPTIONS_READ $3 "variants.ini" "Field 5" "State"
        ${Else}
            !insertmacro MUI_INSTALLOPTIONS_READ $3 "variants.ini" "Field 6" "State"
        ${EndIf}

    StrCpy $5 0

    ${If} $3 != 0
    
        StrCpy $3 7
      next:
        !insertmacro MUI_INSTALLOPTIONS_READ $4 "variants.ini" "Field $3" "Text"
        IfErrors bail
        ${If} $4 == $1
            !insertmacro MUI_INSTALLOPTIONS_READ $5 "variants.ini" "Field $3" "State"
            goto bail
        ${EndIf}
        IntOp $3 $3 + 1
        goto next
      bail:
    ${EndIf}

    ${If} $5 == 0
      !insertmacro UnselectSection $2
    ${Else}
      !insertmacro SelectSection $2
    ${EndIf}

    Pop $4
    Pop $3
    Pop $2
    Pop $1
    Pop $0
FunctionEnd

; Stack 0: top level section index
Function SelectDefaultVariants
    Exch $0
    Push $1
    Push $2
    Push $3
    Push $4

    IntOp $0 $0 + 1

    StrCpy $1 0 ; Last section was group end
    StrCpy $4 "" ; Current compiler
 next:
    SectionGetFlags $0 $2
    IfErrors bail
    IntOp $3 $2 & ${SF_SECGRPEND}
    StrCmp $3 0 not_end
    StrCmp $1 0 0 bail ; two groups in a row means we are backing out
 not_end:
    StrCpy $1 $3
    IntOp $3 $2 & 6
    StrCmp $3 0 0 not_variant
    SectionGetText $0 $2
    Push $4
    Push $2
    Push $0
    call MaybeSelectVariant
    goto variant
 not_variant:
    SectionGetText $0 $4
 variant:
    IntOp $0 $0 + 1
    goto next
 bail:
    Pop $4
    Pop $3
    Pop $2
    Pop $1
    Pop $0
FunctionEnd


Function FindBoostFolder
  Push $0
  Push $1
  Push $2
  Push $3
  Push $4
  Push $5
  Push $6
  
  ${locate::Open} "$PROGRAMFILES\boost" "/F=0 /D=1 /M=boost_*" $0
	${If} $0 != 0
    ${DoUntil} "$FoundBoostFolder" != ""
  	  ${locate::Find} $0 $1 $2 $3 $4 $5 $6
      ${If} "$2" != ""
        StrCpy $FoundBoostFolder $1
      ${EndIf}
    ${Loop}
  ${EndIf}  
  ${locate::Close} $0
  ${locate::Unload}
  
  Pop $6
  Pop $5
  Pop $4
  Pop $3
  Pop $2
  Pop $1
  Pop $0
FunctionEnd


