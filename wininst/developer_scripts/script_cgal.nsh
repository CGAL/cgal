;============================
; Copyright 2007, 2008, 2009 GeometryFactory (France)
; Authors: Andreas Fabri (andreas.fabri@geometryfactrory.com),
;          Fernando Cacciola (fernando.cacciola@geometryfactrory.com),
;          Laurent Rineau (laurent.rineau@geometryfactory.com)
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
Var IsTAUCSInstalled

;--------------------------------
; Macros
;--------------------------------


!define MultiVariantSection "!insertmacro MultiVariantSection"

; Expands to a Section Group named "SecName" which contains all library variants.
; For each variant, the macro "Handler" is expanded with the variant name as argument
!macro MultiVariantSection SecName Handler Platform Idx
  SectionGroup "${SecName}" ${Idx}
    SectionGroup "VC10.0"
      Section /o "Multithread Debug"
        !insertmacro "${Handler}" "${Platform}" "vc100-mt-gd"
      SectionEnd
      Section /o "Multithread"
        !insertmacro "${Handler}" "${Platform}" "vc100-mt"
      SectionEnd
      Section /o "Multithread, static runtime"
        !insertmacro "${Handler}" "${Platform}" "vc100-mt-s"
      SectionEnd
      Section /o "Multithread Debug, static runtime"
        !insertmacro "${Handler}" "${Platform}" "vc100-mt-sgd"
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
  !ifdef DebugLog
    ${LogMsg} "Downloadimg ${SERVER}${SRC_FOLDER}${FILE} into ${TGT}\${FILE}"
  !endif	
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

!macro Install_LAPACK_TAUCS_libs PLATFORM VARIANT
  ; Headers are not VARIANT dependent so we include this only once, but here since
  ; we want to download headers only if at least one lib variant was selected.
  ${If} $IsTAUCSInstalled = 0
    StrCpy $IsTAUCSInstalled 1 
    
    !insertmacro DownloadFile "auxiliary/$Platform/TAUCS-CGAL-3.6/"  "taucs.h.zip"               "$INSTDIR\auxiliary\taucs\include"
    !insertmacro DownloadFile "auxiliary/$Platform/TAUCS-CGAL-3.6/"  "taucs_private.h.zip"       "$INSTDIR\auxiliary\taucs\include"
    !insertmacro DownloadFile "auxiliary/$Platform/TAUCS-CGAL-3.6/"  "taucs_config_tests.h.zip"  "$INSTDIR\auxiliary\taucs\include"
    !insertmacro DownloadFile "auxiliary/$Platform/TAUCS-CGAL-3.6/"  "taucs_config_build.h.zip"  "$INSTDIR\auxiliary\taucs\include"
    !insertmacro DownloadFile "auxiliary/$Platform/TAUCS-CGAL-3.6/"  "blaswrap.h.zip"            "$INSTDIR\auxiliary\taucs\include"

    ${If} "$Platform" == "win32"
      !insertmacro DownloadFile "auxiliary/win32/TAUCS-CGAL-3.6/" "common.zip" "$INSTDIR\auxiliary\taucs\lib"
    ${Endif}
    
  ${Endif}
  
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/TAUCS-CGAL-3.6/"  "libs-${VARIANT}.zip"  "$INSTDIR\auxiliary\taucs\lib"
!macroend

!macro Install_GMP_MPFR_bin PLATFORM
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/GMP/5.0.1/"  "gmp-all.zip"  "$INSTDIR\auxiliary\gmp"
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/MPFR/3.0.0/" "mpfr-all.zip" "$INSTDIR\auxiliary\gmp"
!macroend

!macro _MaybeSelectVariant Compiler Variant Sec1 Sec2

  Push $0
  Push $1
  
  ${If} "${Compiler}" == "VC10.0"
      !insertmacro MUI_INSTALLOPTIONS_READ $0 "variants.ini" "Field 5" "State"
  ${Else}
      !insertmacro MUI_INSTALLOPTIONS_READ $0 "variants.ini" "Field 6" "State"
  ${EndIf}

  ; If the corresponding compiler+variant is not found in the variant page
  ; the section is unselected
  StrCpy $1 0

  ${If} $0 <> 0 ; Is the compiler selected?
  
    ; variants are the fields 7 to 10
    ${For} $0 7 10
    
      !insertmacro MUI_INSTALLOPTIONS_READ $1 "variants.ini" "Field $0" "Text"
      
      ${If} "$1" == "${Variant}" ; Is this variant field the one we are looking for?
      
        ; Found the variant field. Read the state and exit the loop
        !insertmacro MUI_INSTALLOPTIONS_READ $1 "variants.ini" "Field $0" "State"
        
        goto break_${Sec1} ; 
        
      ${EndIf}
      
    ${Next}
    
    break_${Sec1}:
    
  ${EndIf}

  ${If} $1 = 0
    !insertmacro UnselectSection ${Sec1}
    !insertmacro UnselectSection ${Sec2}
  ${Else}
    !insertmacro SelectSection ${Sec1}
    !insertmacro SelectSection ${Sec2}
  ${EndIf}

  Pop $1
  Pop $0
    
!macroend

!define MaybeSelectVariant "!insertmacro _MaybeSelectVariant"
