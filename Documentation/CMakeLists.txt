project(Documentation NONE)

# Visual Studio users might appreciate this
# set_property(GLOBAL PROPERTY USE_FOLDERS ON)

macro(subdirlist result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
      set(dirlist ${dirlist} ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

function(configure_doxygen_package CGAL_PACKAGE_NAME)
  set(CGAL_PACKAGE_DIR ${CMAKE_SOURCE_DIR}/${CGAL_PACKAGE_NAME})
  subdirlist(SUBPACKAGES ${CGAL_PACKAGE_DIR}/doc)
  foreach(package ${SUBPACKAGES})
    set(CGAL_PACKAGE_NAME ${package}) # variable used for configuration
    if(EXISTS ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/Doxyfile.in)
      set(CGAL_DOC_PACKAGE_DEFAULTS ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}_defaults.dxy)
      file(REMOVE ${CGAL_DOC_PACKAGE_DEFAULTS})
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "@INCLUDE = ${CGAL_DOC_DOXY_DEFAULT}\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "OUTPUT_DIRECTORY = ${CGAL_DOC_OUTPUT_DIR}/${CGAL_PACKAGE_NAME}\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "HTML_OUTPUT = .\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "SHOW_FILES = NO\n")

      if(CGAL_DOC_CREATE_LOGS)
        file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "WARN_LOGFILE = ${CGAL_DOC_LOG_DIR}/${CGAL_PACKAGE_NAME}.log\n")
      endif()

      if (CGAL_DOC_RELEASE)
        file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_TODOLIST = NO\n")
      endif()

      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "IMAGE_PATH = ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/fig\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "EXAMPLE_PATH = ${CGAL_PACKAGE_DIR}/examples\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_TAGFILE = ${CGAL_DOC_TAG_DIR}/${CGAL_PACKAGE_NAME}.tag\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_PATH = ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_PATH += ${CGAL_PACKAGE_DIR}/include/\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_INC_PATH = ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_INC_PATH += ${CGAL_PACKAGE_DIR}/include/\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "ALIASES += \"cgalPkgDescriptionBegin{2}=\\details \"\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "ALIASES += \"cgalPkgManuals{2}=<BR>\"\n")
      
      # figure out the dependencies
      if(EXISTS ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/dependencies)
        file(STRINGS ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/dependencies DEPENDENCIES)
      else()
        set(DEPENDENCIES "")
      endif()
      
      foreach(depend ${DEPENDENCIES})
        file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} 
          "TAGFILES += ${CGAL_DOC_TAG_DIR}/${depend}.tag=../${depend}\n")
      endforeach()
      
      configure_file(${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/Doxyfile.in 
        ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}.dxy)
      
      # TODO we also want to run html_post_process per package as a custom_command with POST_BUILD
      add_custom_target(${CGAL_PACKAGE_NAME}_doc
        ${DOXYGEN_EXECUTABLE} ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}.dxy)

      set_target_properties(${CGAL_PACKAGE_NAME}_doc PROPERTIES FOLDER Documentation/Packages)

      # don't do this for now
      # foreach(depend ${DEPENDENCIES})
      #   add_dependencies(${CGAL_PACKAGE_NAME}_doc
      #     ${depend}_doc)
      # endforeach()
      
    endif()
  endforeach()
endfunction()

option(BUILD_DOC "Build the doxygen documentation" OFF)

if (BUILD_DOC)

find_package(Doxygen)
find_package(PythonInterp 2.6.7)

if(DOXYGEN_FOUND)
  # set up the directories and variables
  set(CGAL_DOC_MATHJAX_LOCATION "http://cdn.mathjax.org/mathjax/latest"
      CACHE STRING "The location of MathJax to be used for the documentation.")

  set(CGAL_DOC_OUTPUT_DIR "${CMAKE_BINARY_DIR}/doc_output")
  file(MAKE_DIRECTORY "${CGAL_DOC_OUTPUT_DIR}")
  
  option(CGAL_DOC_CREATE_LOGS "Write log files or print warnings and errors to cerr when building the doc." OFF)

  if(CGAL_DOC_CREATE_LOGS)
    set(CGAL_DOC_LOG_DIR "${CMAKE_BINARY_DIR}/doc_log")
    file(MAKE_DIRECTORY "${CGAL_DOC_LOG_DIR}")
  endif()

  set(CGAL_DOC_TAG_DIR "${CMAKE_BINARY_DIR}/doc_tags")
  file(MAKE_DIRECTORY "${CGAL_DOC_TAG_DIR}")
  set(CGAL_DOC_DXY_DIR "${CMAKE_BINARY_DIR}/doc_dxy")
  file(MAKE_DIRECTORY "${CGAL_DOC_DXY_DIR}")

  set(CGAL_DOC_HEADER ${CGAL_DOC_DXY_DIR}/header.html)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/header.html ${CGAL_DOC_HEADER} @ONLY)
  set(CGAL_DOC_HEADER_PACKAGE ${CGAL_DOC_DXY_DIR}/header_package.html)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/header_package.html ${CGAL_DOC_HEADER_PACKAGE} @ONLY)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/BaseDoxyfile.in ${CGAL_DOC_DXY_DIR}/BaseDoxyfile)
  
  set(CGAL_DOC_DOXY_DEFAULT "${CGAL_DOC_DXY_DIR}/BaseDoxyfile")

  # pkglist_filter gets the path to the pkglist_filter of this source
  # directory.
  if(WIN32)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pkglist_filter.bat ${CMAKE_BINARY_DIR}/pkglist_filter)
  else()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pkglist_filter ${CMAKE_BINARY_DIR}/pkglist_filter)
  endif()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pkglist_filter.py ${CMAKE_BINARY_DIR}/pkglist_filter.py)

  foreach(package ${CGAL_CONFIGURED_PACKAGES_NAMES})
    configure_doxygen_package(${package})
  endforeach()

  
  # hard-coded. this is not recognized as a package by the
  # build-system, because it doesn't have a package_description subdir.
  configure_doxygen_package("Miscellany")

  # do the main package manually, it isn't part of ${CGAL_CONFIGURED_PACKAGES_NAMES}
  configure_doxygen_package("Documentation")
  set_target_properties(Documentation_doc PROPERTIES FOLDER Documentation)

  # hard-code that Documentation_doc depends on all packages
  file(STRINGS ${CMAKE_SOURCE_DIR}/Documentation/doc/Documentation/dependencies DEPENDENCIES)
  foreach(depend ${DEPENDENCIES})
    add_dependencies(Documentation_doc ${depend}_doc)
  endforeach()

  if(PYTHONINTERP_FOUND)
    add_custom_target(Documentation_postprocess
      ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/html_output_post_processing.py --output ${CGAL_DOC_OUTPUT_DIR} --resources ${CMAKE_CURRENT_SOURCE_DIR}/resources
      DEPENDS Documentation_doc)

    if(CGAL_DOC_CREATE_LOGS)
      add_custom_target(Documentation_test
        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testsuite.py --output-dir ${CGAL_DOC_OUTPUT_DIR} --doc-log-dir ${CGAL_DOC_LOG_DIR}
        DEPENDS Documentation_doc)

      add_custom_target(Documentation_test_publish
        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testsuite.py --output-dir ${CGAL_DOC_OUTPUT_DIR} --doc-log-dir ${CGAL_DOC_LOG_DIR} --publish ${CGAL_DOC_PUBLISH_DIR} --do-copy-results
        DEPENDS Documentation_doc)

      set(CGAL_DOC_VERSION_H "/tmp/version.h" CACHE PATH "Path to CGAL/version.h")

      add_custom_target(Documentation_testsuite_publish
        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testsuite.py --output-dir ${CGAL_DOC_OUTPUT_DIR} --doc-log-dir ${CGAL_DOC_LOG_DIR} --publish ${CGAL_DOC_PUBLISH_DIR} --do-copy-results --cgal-version ${CGAL_DOC_VERSION_H} --version-to-keep 10
        DEPENDS Documentation_postprocess)
    endif()
  endif()
  
  
  # add_custom_target(doc) dummy to depend on all docs, require to run
  # twice or something
endif()

endif(BUILD_DOC)
