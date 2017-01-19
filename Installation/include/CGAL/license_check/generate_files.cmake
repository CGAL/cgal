file(STRINGS package_list.txt package_list)

foreach(LINE ${package_list})
  string(REGEX MATCH "^[^ ]*" PACKAGE_CODE_NAME ${LINE})
  string(REPLACE "${PACKAGE_CODE_NAME}" "" PACKAGE_PRETTY_NAME ${LINE})
  string(STRIP "${PACKAGE_PRETTY_NAME}" PACKAGE_PRETTY_NAME)
  configure_file(gpl.h.in ${PACKAGE_CODE_NAME}.h @ONLY NEWLINE_STYLE LF)
endforeach()
