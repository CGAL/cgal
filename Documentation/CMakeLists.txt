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

      if(CGAL_DOC_CREATE_LOGS)
        file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "WARN_LOGFILE = ${CGAL_DOC_LOG_DIR}/${CGAL_PACKAGE_NAME}.log\n")
      endif()

      if (CGAL_DOC_RELEASE)
        file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_TODOLIST = NO\n")
        file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_BUGLIST = NO\n")
      endif()

      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "IMAGE_PATH = ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME}/fig\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "EXAMPLE_PATH = ${CGAL_PACKAGE_DIR}/examples\n")
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_TAGFILE = ${CGAL_DOC_TAG_GEN_DIR}/${CGAL_PACKAGE_NAME}.tag\n")
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
      # target that build the doc and put the tag file in the generation directory CGAL_DOC_TAG_GEN_DIR
      add_custom_target(${CGAL_PACKAGE_NAME}_internal_doxygen_run
        ${DOXYGEN_EXECUTABLE} ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}.dxy)

      set_target_properties(${CGAL_PACKAGE_NAME}_internal_doxygen_run PROPERTIES FOLDER Documentation/Packages)

      #target moving the tag file into the read directory CGAL_DOC_TAG_DIR
      if ( "${CGAL_PACKAGE_NAME}" STREQUAL "Documentation" )
        set(CGAL_PACKAGE_TAGFILE "Manual.tag")
      else()
        set(CGAL_PACKAGE_TAGFILE "${CGAL_PACKAGE_NAME}.tag")
      endif()
      add_custom_target(${CGAL_PACKAGE_NAME}_copy_doc_tags
                        ${CMAKE_COMMAND} -E copy
                          "${CGAL_DOC_TAG_GEN_DIR}/${CGAL_PACKAGE_TAGFILE}"
                          "${CGAL_DOC_TAG_DIR}/${CGAL_PACKAGE_TAGFILE}"
      )

      #add the doc target doing both the doc generation and then the tag file copy
      add_custom_target(${CGAL_PACKAGE_NAME}_doc
        ${DOXYGEN_EXECUTABLE} ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}.dxy
        COMMAND
        ${CMAKE_COMMAND} -E copy
                          "${CGAL_DOC_TAG_GEN_DIR}/${CGAL_PACKAGE_TAGFILE}"
                          "${CGAL_DOC_TAG_DIR}/${CGAL_PACKAGE_TAGFILE}"
      )

      #     ${depend}_doc)
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

  #we use two directories for the generation/reading of tag files to prevent issues
  #if the targets are built in parallel
  set(CGAL_DOC_TAG_GEN_DIR "${CMAKE_BINARY_DIR}/doc_gen_tags")
  file(MAKE_DIRECTORY "${CGAL_DOC_TAG_GEN_DIR}")
  set(CGAL_DOC_TAG_DIR "${CMAKE_BINARY_DIR}/doc_tags")
  file(MAKE_DIRECTORY "${CGAL_DOC_TAG_DIR}")
  set(CGAL_DOC_DXY_DIR "${CMAKE_BINARY_DIR}/doc_dxy")
  file(MAKE_DIRECTORY "${CGAL_DOC_DXY_DIR}")

  set(CGAL_DOC_HEADER ${CGAL_DOC_DXY_DIR}/header.html)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/header.html ${CGAL_DOC_HEADER} @ONLY)
  set(CGAL_DOC_HEADER_PACKAGE ${CGAL_DOC_DXY_DIR}/header_package.html)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/header_package.html ${CGAL_DOC_HEADER_PACKAGE} @ONLY)

  #Generate the bibtex file
  #this is needed also for BaseDoxyfile.in
  if ( ${CMAKE_VERSION} VERSION_GREATER "2.8.10")
    string(TIMESTAMP CGAL_BUILD_YEAR4 "%Y")
    string(TIMESTAMP CGAL_BUILD_YEAR2 "%y")
    string(TIMESTAMP CGAL_BUILD_MONTH "%m")
  else()
    MESSAGE(WARNING "You're using an old version of CMake, date in bibtex files will be incorrect")
    set(CGAL_BUILD_YEAR4 "2000")
    set(CGAL_BUILD_YEAR2 "00")
    set(CGAL_BUILD_MONTH "1")
  endif()

  if (${CGAL_BUILD_MONTH} GREATER 6)
    SET(CGAL_RELEASE_YEAR_ID "${CGAL_BUILD_YEAR2}b")
  else()
    SET(CGAL_RELEASE_YEAR_ID "${CGAL_BUILD_YEAR2}a")
  endif()

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/resources/how_to_cite_cgal.bib.in
  ${CMAKE_BINARY_DIR}/how_to_cite_cgal.bib)

  find_program(BIBTEX2HTML_EXECUTABLE bibtex2html DOC "Path to bibtex2html")
  find_program(SH_EXECUTABLE sh DOC "Path to sh")
  if (BIBTEX2HTML_EXECUTABLE AND SH_EXECUTABLE)
    set(ENV{openout_any} "a") # may be needed with TexLive
    execute_process(COMMAND ${SH_EXECUTABLE}
                            ${CMAKE_CURRENT_SOURCE_DIR}/scripts/makebiblio
                            ${CMAKE_BINARY_DIR}/how_to_cite_cgal.bib
                            ${BIBTEX2HTML_EXECUTABLE}
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    OUTPUT_VARIABLE makebiblio_OUTPUT
                    ERROR_VARIABLE makebiblio_ERROR
                    RESULT_VARIABLE makebiblio_RESULT) #bibtex2html also prints info in stderr
    if(makebiblio_RESULT GREATER 0)
      message(WARNING "
= makebiblio returned the error number ${makebiblio_RESULT}

== Its standard output was:
${makebiblio_OUTPUT}
== Its error output was:
${makebiblio_ERROR}
"
        )
    endif()
  else()
    MESSAGE(WARNING "bibtex2html or sh missing, using default how_to_cite files")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/resources/how_to_cite.html.default ${CMAKE_BINARY_DIR}/how_to_cite.html COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/resources/how_to_cite_cgal.txt.default ${CMAKE_BINARY_DIR}/how_to_cite_cgal.txt COPYONLY)
  endif()
  ## end bibtex configuration

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/BaseDoxyfile.in ${CGAL_DOC_DXY_DIR}/BaseDoxyfile)
  
  set(CGAL_DOC_DOXY_DEFAULT "${CGAL_DOC_DXY_DIR}/BaseDoxyfile")

  # pkglist_filter gets the path to the pkglist_filter of this source
  # directory.
  if(WIN32)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pkglist_filter.bat ${CMAKE_BINARY_DIR}/pkglist_filter.bat)
  endif()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pkglist_filter ${CMAKE_BINARY_DIR}/pkglist_filter)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pkglist_filter.py ${CMAKE_BINARY_DIR}/pkglist_filter.py)

  foreach(package ${CGAL_CONFIGURED_PACKAGES_NAMES})
    configure_doxygen_package(${package})
  endforeach()

  # hard-coded. this is not recognized as a package by the
  # build-system, because it doesn't have a package_description subdir.
  configure_doxygen_package("Miscellany")

  # Add a custom target "doc"
  add_custom_target(doc)
  add_custom_target(doc_pre)
  add_custom_target(doc_post)

  # do the main package manually, it isn't part of ${CGAL_CONFIGURED_PACKAGES_NAMES}
  configure_doxygen_package("Documentation")
  set_target_properties(Documentation_doc PROPERTIES FOLDER Documentation)

  # hard-code that doc depends on all packages
  file(STRINGS ${CMAKE_SOURCE_DIR}/Documentation/doc/Documentation/dependencies DEPENDENCIES)
  foreach(depend ${DEPENDENCIES})
    add_dependencies(doc_pre ${depend}_internal_doxygen_run)
    add_dependencies(doc_post ${depend}_copy_doc_tags)
    add_dependencies(${depend}_copy_doc_tags doc_pre)
  endforeach()
  add_dependencies(doc_pre Documentation_internal_doxygen_run)
  add_dependencies(doc_post Documentation_copy_doc_tags)
  add_dependencies(Documentation_copy_doc_tags doc_pre)

  #total level doc dependencies
  add_dependencies(doc doc_post)

  if(PYTHONINTERP_FOUND)
    add_custom_target(doc_with_postprocessing
      ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/html_output_post_processing.py --output ${CGAL_DOC_OUTPUT_DIR} --resources ${CMAKE_CURRENT_SOURCE_DIR}/resources
    )
    add_dependencies(doc_with_postprocessing doc)

    if(CGAL_DOC_CREATE_LOGS)
      add_custom_target(Documentation_test
        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testsuite.py --output-dir ${CGAL_DOC_OUTPUT_DIR} --doc-log-dir ${CGAL_DOC_LOG_DIR}
      )
      add_dependencies(Documentation_test doc)

      if (BIBTEX2HTML_EXECUTABLE AND SH_EXECUTABLE)
        add_custom_target(Documentation_test_publish
          ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testsuite.py --output-dir ${CGAL_DOC_OUTPUT_DIR} --doc-log-dir ${CGAL_DOC_LOG_DIR} --publish ${CGAL_DOC_PUBLISH_DIR} --do-copy-results
        )
        add_dependencies(Documentation_test_publish doc)

        set(CGAL_DOC_VERSION_H "/tmp/version.h" CACHE PATH "Path to CGAL/version.h")

        add_custom_target(doc_and_publish_testsuite
          ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testsuite.py --output-dir ${CGAL_DOC_OUTPUT_DIR} --doc-log-dir ${CGAL_DOC_LOG_DIR} --publish ${CGAL_DOC_PUBLISH_DIR} --do-copy-results --cgal-version ${CGAL_DOC_VERSION_H} --version-to-keep 10
        )
        add_dependencies(doc_and_publish_testsuite doc_with_postprocessing)
      endif()
    endif()
  endif()
  
  
  # add_custom_target(doc) dummy to depend on all docs, require to run
  # twice or something
endif()

endif(BUILD_DOC)
