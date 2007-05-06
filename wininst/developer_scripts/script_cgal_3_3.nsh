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

;!define TestingOnly

Var no_default_compilers
Var no_default_variants
Var selected_libs

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

!macro InstallThirdPartyLib SRC_MISC SRC_INC SRC_LIB TGT
!ifndef TestingOnly
  SetOutPath "$INSTDIR\auxiliary\${TGT}"
  File /r "${SRC_MISC}"
  SetOutPath "$INSTDIR\auxiliary\${TGT}\include"
  File /r "${SRC_INC}"
  SetOutPath "$INSTDIR\auxiliary\${TGT}\lib"
  File /r "${SRC_LIB}"
!endif
!macroend

!macro Install_CGAL_libs VARIANT
!ifndef TestingOnly
  SetOutPath "$INSTDIR\lib"
  File /r "${CGAL_SRC}\lib\cgal-${VARIANT}.lib"
  File /r "${CGAL_SRC}\lib\CGALcore++-${VARIANT}.lib"
!endif  
!macroend

!macro Install_GMP_MPFR VARIANT
  !insertmacro InstallThirdPartyLib "${GMP_SRC}\README" "${GMP_SRC}\*.h" "${GMP_SRC}\*-${VARIANT}.lib" "gmp"
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


