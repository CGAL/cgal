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

;--------------------------------
; Macros
;--------------------------------


!define MultiVariantSection "!insertmacro MultiVariantSection"

; Expands to a Section Group named "SecName" which contains all library variants.
; For each variant, the macro "Handler" is expanded with the variant name as argument
!macro MultiVariantSection SecName Handler Idx
  SectionGroup "${SecName}" ${Idx}
    SectionGroup "VC7.1"
      Section /o "Multithread Debug"
        !insertmacro "${Handler}" "vc71-mt-gd"
      SectionEnd
      Section /o "Multithread"
        !insertmacro "${Handler}" "vc71-mt"
      SectionEnd
      Section /o "Multithread, static runtime"
        !insertmacro "${Handler}" "vc71-mt-s"
      SectionEnd
      Section /o "Multithread Debug, static runtime"
        !insertmacro "${Handler}" "vc71-mt-sgd"
      SectionEnd
      Section /o "Single thread, static runtime"
        !insertmacro "${Handler}" "vc71-s"
      SectionEnd
      Section /o "Single thread Debug, static runtime"
        !insertmacro "${Handler}" "vc71-sgd"
      SectionEnd
    SectionGroupEnd
    SectionGroup "VC8.0"
      Section /o "Multithread Debug"
        !insertmacro "${Handler}" "vc80-mt-gd"
      SectionEnd
      Section /o "Multithread"
        !insertmacro "${Handler}" "vc80-mt"
      SectionEnd
      Section /o "Multithread, static runtime"
        !insertmacro "${Handler}" "vc80-mt-s"
      SectionEnd
      Section /o "Multithread Debug, static runtime"
        !insertmacro "${Handler}" "vc80-mt-sgd"
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

!macro Install_CGAL_libs VARIANT
!ifndef FetchLocal
  !insertmacro DownloadFile "CGAL/3.3/" "cgal-${VARIANT}.lib.zip"         "$INSTDIR\lib"
  !insertmacro DownloadFile "CGAL/3.3/" "CGALcore++-${VARIANT}.lib.zip"   "$INSTDIR\lib"
  !insertmacro DownloadFile "CGAL/3.3/" "CGALimageIO-${VARIANT}.lib.zip"  "$INSTDIR\lib"
  !insertmacro DownloadFile "CGAL/3.3/" "CGALPDB-${VARIANT}.lib.zip"      "$INSTDIR\lib"
  !insertmacro DownloadFile "CGAL/3.3/" "cgalqt-${VARIANT}.lib.zip"       "$INSTDIR\lib"
!else
  !ifndef SkipFiles
    SetOutPath "$INSTDIR\lib"
    File /nonfatal "${CGAL_SRC}\lib\cgal-${VARIANT}.lib"
    File /nonfatal "${CGAL_SRC}\lib\CGALcore++-${VARIANT}.lib"
    File /nonfatal "${CGAL_SRC}\lib\CGALimageIO-${VARIANT}.lib"
    File /nonfatal "${CGAL_SRC}\lib\CGALPDB-${VARIANT}.lib"
  !endif  
!endif  
!macroend

!macro Install_CGAL_pdbs VARIANT
!ifndef FetchLocal
  !insertmacro DownloadFile "CGAL/3.3/" "cgal-${VARIANT}.pdb.zip"         "$INSTDIR\lib"
  !insertmacro DownloadFile "CGAL/3.3/" "cgalcore++-${VARIANT}.pdb.zip"   "$INSTDIR\lib"
  !insertmacro DownloadFile "CGAL/3.3/" "cgalimageio-${VARIANT}.pdb.zip"  "$INSTDIR\lib"
  !insertmacro DownloadFile "CGAL/3.3/" "cgalpdb-${VARIANT}.pdb.zip"      "$INSTDIR\lib"
  !insertmacro DownloadFile "CGAL/3.3/" "cgalqt-${VARIANT}.pdb.zip"       "$INSTDIR\lib"
!else
  !ifndef SkipFiles
    SetOutPath "$INSTDIR\lib"
    File /nonfatal "${CGAL_SRC}\lib\cgal-${VARIANT}.pdb"
    File /nonfatal "${CGAL_SRC}\lib\cgalcore++-${VARIANT}.pdb"
    File /nonfatal "${CGAL_SRC}\lib\cgalimageio-${VARIANT}.pdb"
    File /nonfatal "${CGAL_SRC}\lib\cgalpdb-${VARIANT}.pdb"
  !endif  
!endif  
!macroend

!macro Install_PDB_if_debug_variant HANDLER VARIANT
  ${StrStr} $R0 ${VARIANT} "gd"
  ${If} "$R0" != ""
    !insertmacro "${HANDLER}" ${VARIANT}
  ${EndIf}  
!macroend

!macro Install_CGAL_bin VARIANT
  !insertmacro Install_CGAL_libs ${VARIANT}
  !insertmacro Install_PDB_if_debug_variant Install_CGAL_pdbs ${VARIANT}
!macroend

!macro Install_GMP_MPFR_libs VARIANT
!ifndef FetchLocal
    !insertmacro DownloadFile "auxiliary/GMP/4.2.1/"  "gmp-${VARIANT}.lib.zip"  "$INSTDIR\auxiliary\gmp\lib"
    !insertmacro DownloadFile "auxiliary/MPFR/2.2.1/" "mpfr-${VARIANT}.lib.zip" "$INSTDIR\auxiliary\gmp\lib"
!else
  !ifndef SkipFiles
    SetOutPath "$INSTDIR\auxiliary\gmp\lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\gmp-${VARIANT}.lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\mpfr-${VARIANT}.lib"
  !endif  
!endif
!macroend

!macro Install_GMP_MPFR_pdbs VARIANT
!ifndef FetchLocal
    !insertmacro DownloadFile "auxiliary/GMP/4.2.1/"  "gmp-${VARIANT}.pdb.zip"  "$INSTDIR\auxiliary\gmp\lib"
    !insertmacro DownloadFile "auxiliary/MPFR/2.2.1/" "mpfr-${VARIANT}.pdb.zip" "$INSTDIR\auxiliary\gmp\lib"
!else
  !ifndef SkipFiles
    SetOutPath "$INSTDIR\auxiliary\gmp\lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\gmp-${VARIANT}.pdb"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\mpfr-${VARIANT}.pdb"
  !endif  
!endif
!macroend

!macro Install_GMP_MPFR_bin VARIANT
  ${If} $IsGmpInstalled = 0
    StrCpy $IsGmpInstalled 1 
  !ifndef SkipFiles
    SetOutPath "$INSTDIR\auxiliary\gmp\include"
    File /r "${CGAL_SRC}\auxiliary\gmp\include\*.*"
    SetOutPath "$INSTDIR\auxiliary\gmp\lib"
    File /nonfatal "${CGAL_SRC}\auxiliary\gmp\lib\README"
  !endif
  ${Endif}

  !insertmacro Install_GMP_MPFR_libs ${VARIANT}
  !insertmacro Install_PDB_if_debug_variant Install_GMP_MPFR_pdbs ${VARIANT}
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

        ${If} $0 == "VC7.1"
            !insertmacro MUI_INSTALLOPTIONS_READ $3 "default_variants.ini" "Field 4" "State"
        ${Else}
            !insertmacro MUI_INSTALLOPTIONS_READ $3 "default_variants.ini" "Field 5" "State"
        ${EndIf}

    StrCpy $5 0

    ${If} $3 != 0
    
        StrCpy $3 6
      next:
        !insertmacro MUI_INSTALLOPTIONS_READ $4 "default_variants.ini" "Field $3" "Text"
        IfErrors bail
        ${If} $4 == $1
            !insertmacro MUI_INSTALLOPTIONS_READ $5 "default_variants.ini" "Field $3" "State"
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

Function FixupProjectFile
  Exch $0
  DetailPrint "Removing CGAL_USE_GMP from $0"
  !insertmacro ReplaceInFile $0 "CGAL_USE_GMP" ""
  Pop $0
FunctionEnd

Function FixupProjectFiles
  Push $0
  Push $1
  Push $2
  Push $3
  Push $4
  Push $5
  Push $6
  
  ${locate::Open} "$INSTDIR" "/D=0 /X=vcproj" $0
	${If} $0 != 0
    ${Do}
  	  ${locate::Find} $0 $1 $2 $3 $4 $5 $6
      ${If} "$1" != ""
        Push $1
        Call FixupProjectFile
      ${EndIf}
    ${LoopUntil} "$1" == ""
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


