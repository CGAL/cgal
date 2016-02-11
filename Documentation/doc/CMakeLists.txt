project(Documentation NONE)

option(BUILD_DOC "Build the doxygen documentation" OFF)
if(NOT BUILD_DOC)
  return()
endif()

find_package(Doxygen)
find_package(PythonInterp 2.6.7)

if(NOT DOXYGEN_FOUND)
  message(WARNING "Cannot build the documentation without Doxygen!")
  return()
endif()

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
  if(CGAL_BRANCH_BUILD)
    set(CGAL_PACKAGE_DIR ${CMAKE_SOURCE_DIR}/${CGAL_PACKAGE_NAME})
  else()
    set(CGAL_PACKAGE_DIR ${CMAKE_SOURCE_DIR})
  endif()

  set(CGAL_PACKAGE_DOC_DIR ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME})
  set(CGAL_PACKAGE_INCLUDE_DIR ${CGAL_PACKAGE_DIR}/include/)

  if(NOT EXISTS ${CGAL_PACKAGE_DOC_DIR}/Doxyfile.in)
    return()
  endif()

  set(CGAL_DOC_PACKAGE_DEFAULTS ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}_defaults.dxy)
  file(REMOVE ${CGAL_DOC_PACKAGE_DEFAULTS})
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "@INCLUDE = ${CGAL_DOC_DOXY_DEFAULT}\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "OUTPUT_DIRECTORY = ${CGAL_DOC_OUTPUT_DIR}/${CGAL_PACKAGE_NAME}\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "HTML_OUTPUT = .\n")

  if(CGAL_DOC_CREATE_LOGS)
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "WARN_LOGFILE = ${CGAL_DOC_LOG_DIR}/${CGAL_PACKAGE_NAME}.log\n")
  endif()

  if(CGAL_DOC_RELEASE)
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_TODOLIST = NO\n")
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_BUGLIST = NO\n")
  endif()

  if(EXISTS "${CGAL_PACKAGE_DIR}/examples")
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "EXAMPLE_PATH = ${CGAL_PACKAGE_DIR}/examples\n")
  endif()
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_TAGFILE = ${CGAL_DOC_TAG_GEN_DIR}/${CGAL_PACKAGE_NAME}.tag\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_PATH = ${CGAL_PACKAGE_DOC_DIR}/\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_PATH += ${CGAL_PACKAGE_DIR}/include/\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_INC_PATH = ${CGAL_PACKAGE_DOC_DIR}/\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "STRIP_FROM_INC_PATH += ${CGAL_PACKAGE_DIR}/include/\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "ALIASES += \"cgalPkgDescriptionBegin{2}=\\details \"\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "ALIASES += \"cgalPkgManuals{2}=<BR>\"\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "INPUT = ${CGAL_PACKAGE_DOC_DIR}\n")
  if(NOT EXISTS "${CGAL_PACKAGE_DOC_DIR}/CGAL")
    # This package has in-source documentation.
    if(CGAL_BRANCH_BUILD) # use the entire include subdir
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "INPUT += ${CGAL_PACKAGE_DIR}/include/\n")
    else() # use the filelist
      if(EXISTS "${CGAL_PACKAGE_DOC_DIR}/filelist.txt")
        file(STRINGS "${CGAL_PACKAGE_DOC_DIR}/filelist.txt" CGAL_PKG_FILES)
        foreach(pkg_file ${CGAL_PKG_FILES})
          file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "INPUT += ${CGAL_PACKAGE_DIR}/include/${pkg_file}\n")
        endforeach()
      endif()
    endif()
  endif()

  # IMAGE_PATH is set by default. For Documentation, we generate the extra path using packages.txt
  set(IMAGE_PATHS "${CGAL_PACKAGE_DOC_DIR}/fig")
  set(DEPENDENCIES "")
  # figure out the dependencies, using the file `dependencies` or `packages.txt` for Documentation
  if (${CGAL_PACKAGE_NAME} STREQUAL "Documentation")
    set(PackageFile ${CGAL_PACKAGE_DOC_DIR}/packages.txt)
    file(STRINGS ${PackageFile} EntriesAsList REGEX "^\\\\package_listing{.+}$")
    foreach(Line ${EntriesAsList})
      string(REGEX REPLACE "^\\\\package_listing{(.+)}$" "\\1" PKG "${Line}")
      list(APPEND DEPENDENCIES ${PKG})

      if(CGAL_BRANCH_BUILD)
        SET(IMG_DIR "${CMAKE_SOURCE_DIR}/${PKG}/doc/${PKG}/fig")
      else()
        SET(IMG_DIR "${CMAKE_SOURCE_DIR}/doc/${PKG}/fig")
      endif()

      if(EXISTS ${IMG_DIR})
        list(APPEND IMAGE_PATHS ${IMG_DIR})
      endif()
    endforeach()
  else()
    if(EXISTS ${CGAL_PACKAGE_DOC_DIR}/dependencies)
      file(STRINGS ${CGAL_PACKAGE_DOC_DIR}/dependencies DEPENDENCIES)
    endif()
  endif()

  foreach(depend ${DEPENDENCIES})
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} 
      "TAGFILES += ${CGAL_DOC_TAG_DIR}/${depend}.tag=../${depend}\n")
  endforeach()

  foreach(image_path ${IMAGE_PATHS})
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
      "IMAGE_PATH += ${image_path}\n")
  endforeach()

  configure_file(${CGAL_PACKAGE_DOC_DIR}/Doxyfile.in ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}.dxy)

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
endfunction()

# set up the directories and variables
set(CGAL_DOC_MATHJAX_LOCATION "https://cdn.mathjax.org/mathjax/latest"
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

set(CGAL_DOC_RESOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/resources")
set(CGAL_DOC_BIBLIO_DIR "${CMAKE_CURRENT_LIST_DIR}/biblio")
set(CGAL_DOC_SCRIPT_DIR "${CMAKE_CURRENT_LIST_DIR}/scripts")

set(CGAL_DOC_HEADER ${CGAL_DOC_DXY_DIR}/header.html)
configure_file(${CGAL_DOC_RESOURCE_DIR}/header.html ${CGAL_DOC_HEADER} @ONLY)
set(CGAL_DOC_HEADER_PACKAGE ${CGAL_DOC_DXY_DIR}/header_package.html)
configure_file(${CGAL_DOC_RESOURCE_DIR}/header_package.html ${CGAL_DOC_HEADER_PACKAGE} @ONLY)

#Generate the bibtex file
#this is needed also for BaseDoxyfile.in
string(TIMESTAMP CGAL_BUILD_YEAR4 "%Y")
string(TIMESTAMP CGAL_BUILD_YEAR2 "%y")
string(TIMESTAMP CGAL_BUILD_MONTH "%m")

if(${CGAL_BUILD_MONTH} GREATER 6)
  set(CGAL_RELEASE_YEAR_ID "${CGAL_BUILD_YEAR2}b")
else()
  set(CGAL_RELEASE_YEAR_ID "${CGAL_BUILD_YEAR2}a")
endif()

# Introduce our own version variable. This way we avoid tying this too
# closely to the convoluted versioning code and can adapt without a
# huge diff.
set(CGAL_DOC_VERSION ${CGAL_CREATED_VERSION_NUM})

## generate how_to_cite files
if(PYTHONINTERP_FOUND)
  execute_process(COMMAND ${PYTHON_EXECUTABLE}
                          ${CGAL_DOC_SCRIPT_DIR}/generate_how_to_cite.py
                          ${CMAKE_SOURCE_DIR}
                          ${CMAKE_BINARY_DIR}
                          "${CGAL_BRANCH_BUILD}"
                  RESULT_VARIABLE GENERATE_HOW_TO_CITE_RESULT)
  if(NOT GENERATE_HOW_TO_CITE_RESULT EQUAL "0")
    message(FATAL_ERROR "generate_how_to_cite did not exit correctly")
  endif()

  configure_file(${CMAKE_BINARY_DIR}/how_to_cite_cgal.bib.in ${CMAKE_BINARY_DIR}/how_to_cite_cgal.bib)
  configure_file(${CMAKE_BINARY_DIR}/how_to_cite_cgal.txt.in ${CMAKE_BINARY_DIR}/how_to_cite_cgal.txt)
  configure_file(${CMAKE_BINARY_DIR}/how_to_cite.html.in ${CMAKE_BINARY_DIR}/how_to_cite.html)
endif()

configure_file(${CGAL_DOC_RESOURCE_DIR}/BaseDoxyfile.in ${CGAL_DOC_DXY_DIR}/BaseDoxyfile)

set(CGAL_DOC_DOXY_DEFAULT "${CGAL_DOC_DXY_DIR}/BaseDoxyfile")

# pkglist_filter gets the path to the pkglist_filter of this source
# directory.
if(WIN32)
  configure_file(${CGAL_DOC_SCRIPT_DIR}/pkglist_filter.bat ${CMAKE_BINARY_DIR}/pkglist_filter.bat)
endif()
configure_file(${CGAL_DOC_SCRIPT_DIR}/pkglist_filter ${CMAKE_BINARY_DIR}/pkglist_filter)
configure_file(${CGAL_DOC_SCRIPT_DIR}/pkglist_filter.py ${CMAKE_BINARY_DIR}/pkglist_filter.py)

set(CGAL_DOC_PACKAGES "")
if(CGAL_BRANCH_BUILD)
  foreach(pkg ${CGAL_CONFIGURED_PACKAGES_NAMES} "Miscellany")
    if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/${pkg}/doc/${pkg})
      list(APPEND CGAL_DOC_PACKAGES ${pkg})
    endif()
    if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/${pkg}/examples)
      set(CGAL_${pkg}_EXAMPLE_DIR "${CMAKE_SOURCE_DIR}/${pkg}/examples")
    endif()
    if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/${pkg}/include)
      set(CGAL_${pkg}_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/${pkg}/include")
    endif()
    if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/${pkg}/include)
      set(CGAL_${pkg}_DEMO_DIR "${CMAKE_SOURCE_DIR}/${pkg}/demo")
    endif()
  endforeach()
else()
  subdirlist(CGAL_DOC_PACKAGES_TMP ${CMAKE_SOURCE_DIR}/doc)
  set(CGAL_DOC_PACKAGES)
  foreach(package ${CGAL_DOC_PACKAGES_TMP})
    if(EXISTS "${CMAKE_SOURCE_DIR}/doc/${package}/PackageDescription.txt")
      list(APPEND CGAL_DOC_PACKAGES "${package}")
    endif()
  endforeach()

  subdirlist(CGAL_EXAMPLE_PACKAGES ${CMAKE_SOURCE_DIR}/examples)
  subdirlist(CGAL_DEMO_PACKAGES ${CMAKE_SOURCE_DIR}/demo)
  list(APPEND CGAL_EXAMPLE_PACKAGES "BGL") # manually add BGL to the list of packages
  foreach(pkg ${CGAL_EXAMPLE_PACKAGES})
    set(CGAL_${pkg}_EXAMPLE_DIR "${CMAKE_SOURCE_DIR}/examples")
    # This might not be entirely correct, but there is no reliable,
    # easy way to get a list of all package names when not in a
    # branch build.
    set(CGAL_${pkg}_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
  endforeach()
  foreach(pkg ${CGAL_DEMO_PACKAGES})
    set(CGAL_${pkg}_DEMO_DIR "${CMAKE_SOURCE_DIR}/demo")
  endforeach()
endif()

foreach(package ${CGAL_DOC_PACKAGES})
  configure_doxygen_package(${package})
endforeach()

# Add a custom target "doc"
add_custom_target(doc)
add_custom_target(doc_pre)
add_custom_target(doc_post)

# do the main package manually, it isn't part of ${CGAL_CONFIGURED_PACKAGES_NAMES}
configure_doxygen_package("Documentation")
set_target_properties(Documentation_doc PROPERTIES FOLDER Documentation)

# hard-code that doc depends on all packages
foreach(depend ${CGAL_DOC_PACKAGES})
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
  set(CGAL_DOC_TESTSUITE_SCRIPT "${CGAL_DOC_SCRIPT_DIR}/testsuite.py")

  add_custom_target(doc_with_postprocessing
    ${PYTHON_EXECUTABLE} ${CGAL_DOC_SCRIPT_DIR}/html_output_post_processing.py --output ${CGAL_DOC_OUTPUT_DIR} --resources ${CGAL_DOC_RESOURCE_DIR}
  )
  add_dependencies(doc_with_postprocessing doc)

  if(CGAL_DOC_CREATE_LOGS)
    add_custom_target(Documentation_test
      ${PYTHON_EXECUTABLE} ${CGAL_DOC_TESTSUITE_SCRIPT} --output-dir ${CGAL_DOC_OUTPUT_DIR} --doc-log-dir ${CGAL_DOC_LOG_DIR}
    )
    add_dependencies(Documentation_test doc)

    add_custom_target(Documentation_test_publish
      ${PYTHON_EXECUTABLE} ${CGAL_DOC_TESTSUITE_SCRIPT} --output-dir ${CGAL_DOC_OUTPUT_DIR} --doc-log-dir ${CGAL_DOC_LOG_DIR} --publish ${CGAL_DOC_PUBLISH_DIR} --do-copy-results
    )
    add_dependencies(Documentation_test_publish doc)

    add_custom_target(doc_and_publish_testsuite
      ${PYTHON_EXECUTABLE} ${CGAL_DOC_TESTSUITE_SCRIPT} --output-dir ${CGAL_DOC_OUTPUT_DIR} --doc-log-dir ${CGAL_DOC_LOG_DIR} --publish ${CGAL_DOC_PUBLISH_DIR} --do-copy-results --cgal-version "${CGAL_FULL_VERSION}" --version-to-keep 10
    )
    add_dependencies(doc_and_publish_testsuite doc_with_postprocessing)

  endif()
endif()
