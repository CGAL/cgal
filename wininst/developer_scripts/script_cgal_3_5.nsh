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

;!define SkipFiles
;!define SkipSetEnvVar
;!define SkipDownload
!define ViaFTP

Var Platform
Var IsGmpInstalled

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

!macro DownloadFileFrom SERVER SRC_FOLDER FILE TGT
!ifndef SkipDownload
  !ifdef ViaFTP
    inetc::get ${SERVER}${SRC_FOLDER}${FILE} ${TGT}\${FILE}
  !else
    NSISdl::download ${SERVER}${SRC_FOLDER}${FILE} ${TGT}\${FILE}
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
    MessageBox MB_OK "Unable to download ${SERVER}${SRC_FOLDER}${FILE}. Error: $0"
    DetailPrint "ERROR $0: Unable to download ${SERVER}${SRC_FOLDER}${FILE}."
  ${Endif}
!endif	
!macroend

!macro DownloadFile SRC_FOLDER FILE TGT
  !insertmacro DownloadFileFrom ${FTP_SRC} ${SRC_FOLDER} ${FILE} ${TGT}
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
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/GMP/4.2.4/"  "gmp-${VARIANT}.lib.zip"  "$INSTDIR\auxiliary\gmp\lib"
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/MPFR/2.3.2/" "mpfr-${VARIANT}.lib.zip" "$INSTDIR\auxiliary\gmp\lib"
!macroend

!macro Install_GMP_MPFR_dlls PLATFORM VARIANT
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/GMP/4.2.4/"  "gmp-${VARIANT}.dll.zip"  "$INSTDIR\auxiliary\gmp\lib"
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/MPFR/2.3.2/" "mpfr-${VARIANT}.dll.zip" "$INSTDIR\auxiliary\gmp\lib"
!macroend

!macro Install_GMP_MPFR_pdbs PLATFORM VARIANT
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/GMP/4.2.4/"  "gmp-${VARIANT}.pdb.zip"  "$INSTDIR\auxiliary\gmp\lib"
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/MPFR/2.3.2/" "mpfr-${VARIANT}.pdb.zip" "$INSTDIR\auxiliary\gmp\lib"
!macroend

!macro Install_GMP_MPFR_bin PLATFORM VARIANT

  ; Headers are not VARIANT dependent so we include this only once, but here since
  ; we want to download headers only if at least one lib variant was selected.
  ${If} $IsGmpInstalled = 0
    StrCpy $IsGmpInstalled 1 
    !insertmacro DownloadFile "auxiliary/${PLATFORM}/GMP/4.2.4/"  "gmp.h.zip"  "$INSTDIR\auxiliary\gmp\include"
    !insertmacro DownloadFile "auxiliary/${PLATFORM}/MPFR/2.3.2/" "mpfr.h.zip" "$INSTDIR\auxiliary\gmp\include"
  ${Endif}
  
  !insertmacro Install_GMP_MPFR_libs                                "${PLATFORM}" "${VARIANT}"
  !insertmacro Install_PDB_if_debug_variant Install_GMP_MPFR_pdbs   "${PLATFORM}" "${VARIANT}"
  !insertmacro Install_DLL_if_dynamic_variant Install_GMP_MPFR_dlls "${PLATFORM}" "${VARIANT}"
!macroend


;--------------------------------
; Functions
;--------------------------------
 
; Given a section ($2) implicitely corresponding
; to a certain compiler ($0) and variant choice ($2)
; select it or unselect it based on the user choices in the variants page
Function __MaybeSelectVariant
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

    ; If the corresponding compiler+variant is not found in the variant page
    ; the section is unselected
    StrCpy $5 0

    ${If} $3 <> 0 ; Is the compiler selected?
    
      ; variants are the fields 7 to 10
      ${For} $3 7 10
      
        !insertmacro MUI_INSTALLOPTIONS_READ $4 "variants.ini" "Field $3" "Text"
        
        ${If} $4 == $1 ; Is this variant field the one we are looking for?
        
          ; Found the variant field. Read the state and exit the loop
          !insertmacro MUI_INSTALLOPTIONS_READ $5 "variants.ini" "Field $3" "State"
          
          goto break ; 
          
        ${EndIf}
        
      ${Next}
      
      break:
      
    ${EndIf}

    ${If} $5 = 0
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

!macro _MaybeSelectVariant Compiler Variant Sec

  Push "${Compiler}"
  Push "${Variant}"
  Push "${Sec}"
  
  call __MaybeSelectVariant
  
!macroend

!define MaybeSelectVariant "!insertmacro _MaybeSelectVariant"
