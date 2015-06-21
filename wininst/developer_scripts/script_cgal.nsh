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

!macro Install_GMP_MPFR_bin PLATFORM
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/GMP/5.0.1/"  "gmp-all-CGAL-3.9.zip"  "$INSTDIR\auxiliary\gmp"
  !insertmacro DownloadFile "auxiliary/${PLATFORM}/MPFR/3.0.0/" "mpfr-all-CGAL-3.9.zip" "$INSTDIR\auxiliary\gmp"
!macroend

