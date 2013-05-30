project(Documentation NONE)

macro(subdirlist result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
        set(dirlist ${dirlist} ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
ENDMACRO()

function(configure_doxygen_package CGAL_PACKAGE_NAME)
  set(CGAL_PACKAGE_DIR ${CMAKE_SOURCE_DIR}/${CGAL_PACKAGE_NAME})
  subdirlist(SUBPACKAGES ${CGAL_PACKAGE_DIR}/doc)
  foreach(package ${SUBPACKAGES})
    set(CGAL_PACKAGE_NAME ${package}) # variable used for configuration
    if(EXISTS ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/doxygen.in)
      message("Configuring doxygen for: ${CGAL_PACKAGE_NAME}")
      set(CGAL_DOC_PACKAGE_DEFAULTS ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}_defaults.dxy)
      file(REMOVE ${CGAL_DOC_PACKAGE_DEFAULTS})
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "@INCLUDE = ${CGAL_DOC_DOXY_DEFAULT}\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "OUTPUT_DIRECTORY = ${CGAL_DOC_OUTPUT_DIR}/${CGAL_PACKAGE_NAME}\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "IMAGE_PATH = ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/fig\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "EXAMPLE_PATH = ${CGAL_PACKAGE_DIR}/examples\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_TAGFILE = ${CGAL_DOC_TAG_DIR}/${CGAL_PACKAGE_NAME}.tag\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_PATH = ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_PATH += ${CGAL_PACKAGE_DIR}/include/\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_INC_PATH = ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_INC_PATH += ${CGAL_PACKAGE_DIR}/include/\n")
      
      # figure out the dependencies
      file(STRINGS ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/dependeds.txt DEPENDENCIES)
      message("${CGAL_PACKAGE_NAME} Depends on ${DEPENDENCIES}")
      foreach(depend ${DEPENDENCIES})
        file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "TAGFILES += ${CGAL_DOC_TAG_DIR}/${depend}.tag=${CGAL_DOC_OUTPUT_DIR}/${depend}/html\n")
      endforeach()
      
      configure_file(${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/doxygen.in ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}.dxy)

      
      add_custom_target(${CGAL_PACKAGE_NAME}_doxy
        ${DOXYGEN_EXECUTABLE} ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}.dxy
        DEPENDS CGAL_doc) # list actual dependencies here... 
      
      add_custom_target(${CGAL_PACKAGE_NAME}_doc
        # to be implemented
        # ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/html_output_post_processing.py --package path_to_package
        # depend because we need the styefiles and collective resources
        DEPENDS ${CGAL_PACKAGE_NAME}_doxy)
    endif()
  endforeach()
endfunction()

find_package(Doxygen)
find_package(PythonInterp 2.7)

if(DOXYGEN_FOUND)
  # set up the directories and variables
  set(CGAL_DOC_MATHJAX_LOCATION "http://cdn.mathjax.org/mathjax/latest")

  set(CGAL_DOC_OUTPUT_DIR "${CMAKE_BINARY_DIR}/doc_output")
  file(MAKE_DIRECTORY "${CGAL_DOC_OUTPUT_DIR}")
  set(CGAL_DOC_LOG_DIR "${CMAKE_BINARY_DIR}/doc_log")
  file(MAKE_DIRECTORY "${CGAL_DOC_LOG_DIR}")
  set(CGAL_DOC_TAG_DIR "${CMAKE_BINARY_DIR}/doc_tags")
  file(MAKE_DIRECTORY "${CGAL_DOC_TAG_DIR}")
  set(CGAL_DOC_DXY_DIR "${CMAKE_BINARY_DIR}/doc_dxy")
  file(MAKE_DIRECTORY "${CGAL_DOC_DXY_DIR}")

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CGAL_DOC_DXY_DIR}/BaseDoxyfile)
  set(CGAL_DOC_HEADER ${CGAL_DOC_DXY_DIR}/header.html)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/header.html ${CGAL_DOC_HEADER} @ONLY)
  set(CGAL_DOC_HEADER_PACKAGE ${CGAL_DOC_DXY_DIR}/header_package.html)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/header_package.html ${CGAL_DOC_HEADER_PACKAGE} @ONLY)
  
  set(CGAL_DOC_DOXY_DEFAULT "${CGAL_DOC_DXY_DIR}/BaseDoxyfile")

  # pkglist_filter gets the path to the pkglist_filter of this source
  # directory.
  if(WIN32)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pkglist_filter.bat ${CMAKE_BINARY_DIR}/pkglist_filter)
  else()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pkglist_filter ${CMAKE_BINARY_DIR}/pkglist_filter)
  endif()

  # the main package
  message("Configuring doxygen for: CGAL main package")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CGAL.dxy.in ${CGAL_DOC_DXY_DIR}/CGAL.dxy)
  add_custom_target(CGAL_doc
    ${DOXYGEN_EXECUTABLE} ${CGAL_DOC_DXY_DIR}/CGAL.dxy
    )
  
  foreach(package ${CGAL_CONFIGURED_PACKAGES_NAMES})
    configure_doxygen_package(${package})
  endforeach()

  add_custom_target(html_postprocess
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/html_output_post_processing.py --output ${CGAL_DOC_OUTPUT_DIR}
    )
  
  # add_custom_target(doc) dummy to depend on all docs, require to run
  # twice or something
endif()
