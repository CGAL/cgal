;============================
; Copyright 2007 GeometryFactory (France)
; Author: Andreas Fabri (andreas.fabri@geometryfactrory.com), Fernando Cacciola (fernando.cacciola@geometryfactrory.com)
;============================
; Some portions of this file have been extracted/derived from "boost,.nsi", the Boost Windows Installer, contributed by www.boost-consulting.org.
;
; Copyright 2006 Daniel Wallin
; Copyright 2006 Eric Niebler
; Distributed under the Boost Software License, Version 1.0. (See
; accompanying file LICENSE_1_0.txt or copy at
; http://www.boost.org/LICENSE_1_0.txt)
;============================

Var no_default_compilers
Var no_default_variants
Var selected_libs

;--------------------------------
; Macros
;--------------------------------

!macro AddConfigFlag CFLAG ADDED
  StrCmp $ADDED "y" Set
    return
  Set:  
    Push "#define CGAL_USE_${CFLAG} 1$\r$\n"
    Call AppendLineToCompilerConfig
    StrCpy $ADDED "y"
!macroend

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

!macro InstallThirdPartyLib SRC_MISC SRC_INC SRC_LIB TGT
  SetOutPath "$INSTDIR\auxiliary\${TGT}"
  FILE /r "${SRC_MISC}"
  SetOutPath "$INSTDIR\auxiliary\${TGT}\include"
  FILE /r "${SRC_INC}"
  SetOutPath "$INSTDIR\auxiliary\${TGT}\lib"
  FILE /r "${SRC_LIB}"
!macroend

!macro Install_CGAL_libs VARIANT
  
  SetOutPath "$INSTDIR\lib"
  FILE /r "${CGAL_SRC}\lib\cgal-${VARIANT}.lib"
  FILE /r "${CGAL_SRC}\lib\CGALcore++-${VARIANT}.lib"
!macroend

!macro Install_GMP_MPFR VARIANT
  !insertmacro InstallThirdPartyLib "${GMP_SRC}\README" "${GMP_SRC}\*.h" "${GMP_SRC}\*-${VARIANT}.lib" "gmp"
!macroend

!macro Install_TAUCS VARIANT
  !insertmacro InstallThirdPartyLib "${TAUCS_SRC}\LICENSE" "${TAUCS_SRC}\*.h" "${TAUCS_SRC}\*-${VARIANT}.lib" "taucs"
!macroend

!macro Install_ZLIB VARIANT
  !insertmacro InstallThirdPartyLib "${ZLIB_SRC}\*.txt" "${TAUCS_SRC}\*.h" "${TAUCS_SRC}\*-${VARIANT}.lib" "zlib"
!macroend

;--------------------------------
; Functions
;--------------------------------


# Appends "Line" to "compiler_config.h"
# Usage:
#   Push Line
#   Call  AppendLineToCompilerConfig
Function AppendLineToCompilerConfig
  Push "$INSTDIR\include\CGAL\config\msvc\CGAL\compiler_config.h"
  Call AppendLineToFile
FunctionEnd

# Appends "Line" to "Fle"
# Usage:
#   Push Line
#   Push File
#   Call  AppendLineToFile
Function AppendLineToFile
  Exch $0 # file
  Exch
  Exch $1 # line
  Push $2 # handle 
  FileOpen $2 $0 a
  IfErrors error
  FileSeek $2 0 END
  FileWrite $2 $1
  FileClose $2
  Goto done
  error:
    DetailPrint "ERROR: Unable to add line:$\r$\n  $1 $\r$\nto file:$\r$\n  $0 $\r$\n"
  done:
  Pop $2
  Pop $1
  Pop $0  
FunctionEnd

Function initSelectionFlags

    ${LogText} "initSelectionFlags..."
  
    StrCpy $selected_libs ""
    ClearErrors
    StrCpy $0 0
    
    
  next:
    SectionGetText $0 $1
    IfErrors bail
    ${LogText} "Section $0 is '$1'"
    StrCpy $2 $1 "" -4
    StrCmp $2 "libs" 0 not_lib
    ${LogText} "Section $0 is a libray"
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

    ${LogText} "MaybeSelectVariant comp=$0 variant=$1 section=$2"
    
    Push $3
    Push $4

        ${If} $0 == "VC7.1"
            !insertmacro MUI_INSTALLOPTIONS_READ $3 "default_variants.ini" "Field 4" "State"
        ${Else}
            !insertmacro MUI_INSTALLOPTIONS_READ $3 "default_variants.ini" "Field 5" "State"
        ${EndIf}

    StrCpy $5 0

    ${If} $3 != 0
      ${LogText} "$0 is checked. Looking for variant in selection."
    
        StrCpy $3 6
      next:
        !insertmacro MUI_INSTALLOPTIONS_READ $4 "default_variants.ini" "Field $3" "Text"
        IfErrors bail
        ${If} $4 == $1
            !insertmacro MUI_INSTALLOPTIONS_READ $5 "default_variants.ini" "Field $3" "State"
            ${LogText} "variant found. selected=$5"
            goto bail
        ${EndIf}
        IntOp $3 $3 + 1
        goto next
      bail:
        ${If} $3 == 12
          ${LogText} "ERROR!! variant NOT found in selection"
        ${EndIf} 
    ${EndIf}

    ${If} $5 == 0
      !insertmacro SelectSection $2
    ${Else}
      !insertmacro UnselectSection $2
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

    ${LogText} "SelectDefaultVariants for section $0"
    
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


