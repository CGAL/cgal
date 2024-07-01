cmake_minimum_required(VERSION 3.12...3.29)
project(Documentation NONE)

# Check whether this cmake script is the top level one
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  # decide if this is a branch build
  if(NOT EXISTS "${CMAKE_SOURCE_DIR}/../VERSION")
    set(CGAL_BRANCH_BUILD
        ON
        CACHE INTERNAL "CGAL Branch build")
  endif()

  if(CGAL_BRANCH_BUILD)
    set(CGAL_ROOT "${CMAKE_SOURCE_DIR}/../..")
  else()
    set(CGAL_ROOT "${CMAKE_SOURCE_DIR}/..")
  endif()
  get_filename_component(CGAL_ROOT ${CGAL_ROOT} REALPATH)
else()
  option(BUILD_DOC "Build the doxygen documentation" OFF)
  if(NOT BUILD_DOC)
    return()
  endif()
  set(CGAL_ROOT "${CMAKE_SOURCE_DIR}")
endif()

find_package(Doxygen REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)

if (NOT Python3_EXECUTABLE)
  message(FATAL_ERROR "Cannot build the documentation without Python3!")
  return()
endif()

message(STATUS  ${Python3_EXECUTABLE})

if(NOT DOXYGEN_FOUND)
  message(WARNING "Cannot build the documentation without Doxygen!")
  return()
endif()

#starting from cmake 3.9 the usage of DOXYGEN_EXECUTABLE is deprecated
if(TARGET Doxygen::doxygen)
  get_property(
    DOXYGEN_EXECUTABLE
    TARGET Doxygen::doxygen
    PROPERTY IMPORTED_LOCATION)
endif()

# Visual Studio users might appreciate this
# set_property(GLOBAL PROPERTY USE_FOLDERS ON)

macro(subdirlist result curdir)
  file(
    GLOB children
    RELATIVE ${curdir}
    ${curdir}/*)
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
    set(CGAL_PACKAGE_DIR ${CGAL_ROOT}/${CGAL_PACKAGE_NAME})
  else()
    set(CGAL_PACKAGE_DIR ${CGAL_ROOT})
  endif()

  set(CGAL_PACKAGE_DOC_DIR ${CGAL_PACKAGE_DIR}/doc/${CGAL_PACKAGE_NAME})
  set(CGAL_PACKAGE_INCLUDE_DIR ${CGAL_PACKAGE_DIR}/include/)

  if(NOT EXISTS ${CGAL_PACKAGE_DOC_DIR}/Doxyfile.in)
    return()
  endif()
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${CGAL_PACKAGE_DOC_DIR}/Doxyfile.in)

  set(CGAL_DOC_PACKAGE_DEFAULTS
      ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}_defaults.dxy)
  file(REMOVE ${CGAL_DOC_PACKAGE_DEFAULTS})
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
       "@INCLUDE = ${CGAL_DOC_DOXY_DEFAULT}\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
       "OUTPUT_DIRECTORY = ${CGAL_DOC_OUTPUT_DIR}/${CGAL_PACKAGE_NAME}\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "HTML_OUTPUT = .\n")

  if(CGAL_DOC_CREATE_LOGS)
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
         "WARN_LOGFILE = ${CGAL_DOC_LOG_DIR}/${CGAL_PACKAGE_NAME}.log\n")
  endif()

  if(CGAL_DOC_RELEASE)
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_TODOLIST = NO\n")
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_BUGLIST = NO\n")
  endif()

  if(EXISTS "${CGAL_PACKAGE_DIR}/examples")
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
         "EXAMPLE_PATH = ${CGAL_PACKAGE_DIR}/examples\n")
  endif()
  if(CGAL_GENERATE_XML)
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "GENERATE_XML = YES\n")
  endif()

  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
       "GENERATE_TAGFILE = ${CGAL_DOC_TAG_GEN_DIR}/${CGAL_PACKAGE_NAME}.tag\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
       "STRIP_FROM_PATH = ${CGAL_PACKAGE_DOC_DIR}/\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
       "STRIP_FROM_PATH += ${CGAL_PACKAGE_DIR}/include/\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
       "STRIP_FROM_INC_PATH = ${CGAL_PACKAGE_DOC_DIR}/\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
       "STRIP_FROM_INC_PATH += ${CGAL_PACKAGE_DIR}/include/\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
       "ALIASES += \"cgalPkgDescriptionBegin{2}=\\details \"\n")
  file(
    APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
    "ALIASES += \"cgalPkgManuals{2}=<div class=\\\"PkgManuals\\\"></div>\"\n")
  file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "INPUT = ${CGAL_PACKAGE_DOC_DIR}\n")
  if(NOT EXISTS "${CGAL_PACKAGE_DOC_DIR}/CGAL")
    # This package has in-source documentation.
    if(CGAL_BRANCH_BUILD) # use the entire include subdir
      file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
           "INPUT += ${CGAL_PACKAGE_DIR}/include/\n")
    else() # use the filelist
      if(EXISTS "${CGAL_PACKAGE_DOC_DIR}/filelist.txt")
        file(STRINGS "${CGAL_PACKAGE_DOC_DIR}/filelist.txt" CGAL_PKG_FILES)
        foreach(pkg_file ${CGAL_PKG_FILES})
          file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
               "INPUT += ${CGAL_PACKAGE_DIR}/include/${pkg_file}\n")
        endforeach()
      endif()
    endif()
  endif()
  if(EXISTS "${CGAL_PACKAGE_DIR}/include/CGAL/${CGAL_PACKAGE_NAME}/internal")
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
        "EXCLUDE += ${CGAL_PACKAGE_DIR}/include/CGAL/${CGAL_PACKAGE_NAME}/internal\n")
  endif()

  # IMAGE_PATH is set by default. For Documentation, we generate the extra path using packages.txt
  set(IMAGE_PATHS "${CGAL_PACKAGE_DOC_DIR}/fig")
  set(DEPENDENCIES "")
  # figure out the dependencies, using the file `dependencies` or `packages.txt` for Documentation
  if(${CGAL_PACKAGE_NAME} STREQUAL "Documentation")
    set(PackageFile ${CGAL_PACKAGE_DOC_DIR}/packages.txt)
    file(STRINGS ${PackageFile} EntriesAsList REGEX "^\\\\package_listing{.+}$")
    foreach(Line ${EntriesAsList})
      string(REGEX REPLACE "^\\\\package_listing{(.+)}$" "\\1" PKG "${Line}")
      list(APPEND DEPENDENCIES ${PKG})

      if(CGAL_BRANCH_BUILD)
        set(IMG_DIR "${CGAL_ROOT}/${PKG}/doc/${PKG}/fig")
      else()
        set(IMG_DIR "${CGAL_ROOT}/doc/${PKG}/fig")
      endif()

      if(EXISTS ${IMG_DIR})
        list(APPEND IMAGE_PATHS ${IMG_DIR})
      endif()
    endforeach()
  else()
    if(EXISTS ${CGAL_PACKAGE_DOC_DIR}/dependencies)
      set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${CGAL_PACKAGE_DOC_DIR}/dependencies)
      file(STRINGS ${CGAL_PACKAGE_DOC_DIR}/dependencies DEPENDENCIES)
    endif()
  endif()

  foreach(depend ${DEPENDENCIES})
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS}
         "TAGFILES += ${CGAL_DOC_TAG_DIR}/${depend}.tag=../${depend}\n")
  endforeach()

  foreach(image_path ${IMAGE_PATHS})
    file(APPEND ${CGAL_DOC_PACKAGE_DEFAULTS} "IMAGE_PATH += ${image_path}\n")
  endforeach()

  configure_file(${CGAL_PACKAGE_DOC_DIR}/Doxyfile.in
                 ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}.dxy)

  # TODO we also want to run html_post_process per package as a custom_command with POST_BUILD
  # target that build the doc and put the tag file in the generation directory CGAL_DOC_TAG_GEN_DIR
  add_custom_target(
    ${CGAL_PACKAGE_NAME}_internal_doxygen_run
    ${DOXYGEN_EXECUTABLE} ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}.dxy)

  set_target_properties(${CGAL_PACKAGE_NAME}_internal_doxygen_run
                        PROPERTIES FOLDER Documentation/Packages)

  #target moving the tag file into the read directory CGAL_DOC_TAG_DIR
  if("${CGAL_PACKAGE_NAME}" STREQUAL "Documentation")
    set(CGAL_PACKAGE_TAGFILE "Manual.tag")
  else()
    set(CGAL_PACKAGE_TAGFILE "${CGAL_PACKAGE_NAME}.tag")
  endif()
  add_custom_target(
    ${CGAL_PACKAGE_NAME}_copy_doc_tags
    ${CMAKE_COMMAND} -E copy "${CGAL_DOC_TAG_GEN_DIR}/${CGAL_PACKAGE_TAGFILE}"
    "${CGAL_DOC_TAG_DIR}/${CGAL_PACKAGE_TAGFILE}")

  #add the doc target doing both the doc generation and then the tag file copy
  add_custom_target(
    ${CGAL_PACKAGE_NAME}_doc
    ${DOXYGEN_EXECUTABLE} ${CGAL_DOC_DXY_DIR}/${CGAL_PACKAGE_NAME}.dxy
    COMMAND
      ${CMAKE_COMMAND} -E copy
      "${CGAL_DOC_TAG_GEN_DIR}/${CGAL_PACKAGE_TAGFILE}"
      "${CGAL_DOC_TAG_DIR}/${CGAL_PACKAGE_TAGFILE}")

  #     ${depend}_doc)
  # don't do this for now
  # foreach(depend ${DEPENDENCIES})
  #   add_dependencies(${CGAL_PACKAGE_NAME}_doc
  #     ${depend}_doc)
  # endforeach()
endfunction()

# set up the directories and variables
if(CGAL_DOC_MATHJAX_LOCATION)
  set(CGAL_DOC_MATHJAX_LOCATION_FULL_OPTION_LINE
      "MATHJAX_RELPATH = ${CGAL_DOC_MATHJAX_LOCATION}")
else()
  set(CGAL_DOC_MATHJAX_LOCATION
      "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2")
endif()

set(CGAL_DOC_OUTPUT_DIR "${CMAKE_BINARY_DIR}/doc_output")
file(MAKE_DIRECTORY "${CGAL_DOC_OUTPUT_DIR}")

option(
  CGAL_DOC_CREATE_LOGS
  "Write log files or print warnings and errors to cerr when building the doc."
  OFF)

if(CGAL_DOC_CREATE_LOGS)
  set(CGAL_DOC_LOG_DIR "${CMAKE_BINARY_DIR}/doc_log")
  file(MAKE_DIRECTORY "${CGAL_DOC_LOG_DIR}")
endif()

option(CGAL_EXTRACT_ALL_NO_DETAILED_IF_EMPTY
       "Use CGAL special doxygen setting EXTRACT_ALL_NO_DETAILED_IF_EMPTY." ON)
if(CGAL_EXTRACT_ALL_NO_DETAILED_IF_EMPTY)
  set(CGAL_OPT_EXTRACT_ALL_NO_DETAILED_IF_EMPTY
      "EXTRACT_ALL_NO_DETAILED_IF_EMPTY = YES")
else()
  # The default is NO, so we could leave it out, but it is better to have a commented out placeholder
  # this will work for versions with and without the setting.
  set(CGAL_OPT_EXTRACT_ALL_NO_DETAILED_IF_EMPTY
      "#EXTRACT_ALL_NO_DETAILED_IF_EMPTY = NO")
endif()

option(CGAL_NO_ADDITIONAL_DETAILS
       "Use CGAL special doxygen setting NO_ADDITIONAL_DETAILS." ON)
if(CGAL_NO_ADDITIONAL_DETAILS)
  set(CGAL_OPT_NO_ADDITIONAL_DETAILS
      "NO_ADDITIONAL_DETAILS = YES")
else()
  # The default is NO, so we could leave it out, but it is better to have a commented out placeholder
  # this will work for versions with and without the setting.
  set(CGAL_OPT_NO_ADDITIONAL_DETAILS
      "#NO_ADDITIONAL_DETAILS = NO")
endif()

#we use two directories for the generation/reading of tag files to prevent issues
#if the targets are built in parallel
set(CGAL_DOC_TAG_GEN_DIR "${CMAKE_BINARY_DIR}/doc_gen_tags")
file(MAKE_DIRECTORY "${CGAL_DOC_TAG_GEN_DIR}")
set(CGAL_DOC_TAG_DIR "${CMAKE_BINARY_DIR}/doc_tags")
file(MAKE_DIRECTORY "${CGAL_DOC_TAG_DIR}")
set(CGAL_DOC_DXY_DIR "${CMAKE_BINARY_DIR}/doc_dxy")
file(MAKE_DIRECTORY "${CGAL_DOC_DXY_DIR}")

#Setting the resource directory depending on the version of doxygen
set(CGAL_DOC_RESOURCE_DIR_DEFAULT "${CMAKE_CURRENT_LIST_DIR}/resources/1.10.0")

#  first look if resources for the specific doxygen version is available, fallback
#  on the default otherwise

#select only the version number (not the commit hash)
string(REPLACE " " ";" DOXYGEN_VERSION ${DOXYGEN_VERSION})
list(GET DOXYGEN_VERSION 0 DOXYGEN_VERSION)

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/resources/${DOXYGEN_VERSION}")
  set(CGAL_DOC_RESOURCE_DIR
      "${CMAKE_CURRENT_LIST_DIR}/resources/${DOXYGEN_VERSION}")
else()
  set(CGAL_DOC_RESOURCE_DIR "${CGAL_DOC_RESOURCE_DIR_DEFAULT}")
endif()

set(CGAL_DOC_BIBLIO_DIR "${CMAKE_CURRENT_LIST_DIR}/biblio")
set(CGAL_DOC_SCRIPT_DIR "${CMAKE_CURRENT_LIST_DIR}/scripts")

set(CGAL_DOC_HEADER ${CGAL_DOC_DXY_DIR}/header.html)
configure_file(${CGAL_DOC_RESOURCE_DIR}/header.html ${CGAL_DOC_HEADER} @ONLY)
set(CGAL_DOC_HEADER_PACKAGE ${CGAL_DOC_DXY_DIR}/header_package.html)
configure_file(${CGAL_DOC_RESOURCE_DIR}/header_package.html
               ${CGAL_DOC_HEADER_PACKAGE} @ONLY)

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

if(NOT CGAL_CREATED_VERSION_NUM)
  if(CGAL_BRANCH_BUILD)
    include(${CGAL_ROOT}/CGALConfigVersion.cmake)

    if(CGAL_BUGFIX_VERSION AND CGAL_BUGFIX_VERSION GREATER 0)
      set(CGAL_CREATED_VERSION_NUM
          "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}.${CGAL_BUGFIX_VERSION}")
    else()
      set(CGAL_CREATED_VERSION_NUM
          "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")
    endif()
  else()
    if(EXISTS "${CGAL_ROOT}/doc/public_release_name")
      file(STRINGS "${CGAL_ROOT}/doc/public_release_name" CGAL_VERSION_LINE)
      string(REGEX REPLACE "CGAL-" "" CGAL_CREATED_VERSION_NUM
                           "${CGAL_VERSION_LINE}")
    else()
      #read version.h and get the line with CGAL_VERSION
      file(STRINGS "${CGAL_ROOT}/include/CGAL/version.h" CGAL_VERSION_LINE
           REGEX "CGAL_VERSION ")
      #extract release id
      string(REGEX MATCH "[0-9]+\\.[0-9]+\\.?[0-9]*" CGAL_CREATED_VERSION_NUM
                   "${CGAL_VERSION_LINE}")
    endif()
  endif()
endif()

# Introduce our own version variable. This way we avoid tying this too
# closely to the convoluted versioning code and can adapt without a
# huge diff.
set(CGAL_DOC_VERSION ${CGAL_CREATED_VERSION_NUM})

## generate how_to_cite files
execute_process(
  COMMAND ${Python3_EXECUTABLE} ${CGAL_DOC_SCRIPT_DIR}/generate_how_to_cite.py
          ${CGAL_ROOT} ${CMAKE_BINARY_DIR} "${CGAL_BRANCH_BUILD}"
  RESULT_VARIABLE GENERATE_HOW_TO_CITE_RESULT)
if(NOT GENERATE_HOW_TO_CITE_RESULT EQUAL "0")
  message(FATAL_ERROR "generate_how_to_cite did not exit correctly")
endif()

configure_file(${CMAKE_BINARY_DIR}/how_to_cite_cgal.bib.in
               ${CMAKE_BINARY_DIR}/how_to_cite_cgal.bib)
configure_file(${CMAKE_BINARY_DIR}/how_to_cite_cgal.txt.in
               ${CMAKE_BINARY_DIR}/how_to_cite_cgal.txt)
configure_file(${CMAKE_BINARY_DIR}/how_to_cite.html.in
               ${CMAKE_BINARY_DIR}/how_to_cite.html)

#
configure_file(${CGAL_DOC_RESOURCE_DIR}/BaseDoxyfile.in
               ${CGAL_DOC_DXY_DIR}/BaseDoxyfile)

set(CGAL_DOC_DOXY_DEFAULT "${CGAL_DOC_DXY_DIR}/BaseDoxyfile")

if (BE_QUIET)
  file(APPEND ${CGAL_DOC_DOXY_DEFAULT} "WARN_LOGFILE=doxygen.log\n")
endif()

# pkglist_filter gets the path to the pkglist_filter of this source
# directory.
if(WIN32)
  configure_file(${CGAL_DOC_SCRIPT_DIR}/pkglist_filter.bat
                 ${CMAKE_BINARY_DIR}/pkglist_filter.bat)
endif()
configure_file(${CGAL_DOC_SCRIPT_DIR}/pkglist_filter
               ${CMAKE_BINARY_DIR}/pkglist_filter)
configure_file(${CGAL_DOC_SCRIPT_DIR}/pkglist_filter.py
               ${CMAKE_BINARY_DIR}/pkglist_filter.py)

set(CGAL_DOC_PACKAGES "")
if(CGAL_BRANCH_BUILD)
  file(
    GLOB CGAL_ALL_PACKAGES
    RELATIVE ${CGAL_ROOT}
    "${CGAL_ROOT}/*")
  foreach(pkg ${CGAL_ALL_PACKAGES})
    #detect packages with documentation to process
    if(IS_DIRECTORY ${CGAL_ROOT}/${pkg})
      if(EXISTS "${CGAL_ROOT}/${pkg}/doc/${pkg}/PackageDescription.txt")
        list(APPEND CGAL_DOC_PACKAGES ${pkg})
      endif()
    endif()
    #set up directories of packages
    if(IS_DIRECTORY ${CGAL_ROOT}/${pkg}/examples)
      set(CGAL_${pkg}_EXAMPLE_DIR "${CGAL_ROOT}/${pkg}/examples")
    endif()
    if(IS_DIRECTORY ${CGAL_ROOT}/${pkg}/include)
      set(CGAL_${pkg}_INCLUDE_DIR "${CGAL_ROOT}/${pkg}/include")
    endif()
    if(IS_DIRECTORY ${CGAL_ROOT}/${pkg}/demo)
      set(CGAL_${pkg}_DEMO_DIR "${CGAL_ROOT}/${pkg}/demo")
    endif()
  endforeach()
else()
  subdirlist(CGAL_DOC_PACKAGES_TMP ${CGAL_ROOT}/doc)
  set(CGAL_DOC_PACKAGES)
  foreach(package ${CGAL_DOC_PACKAGES_TMP})
    if(EXISTS "${CGAL_ROOT}/doc/${package}/PackageDescription.txt")
      list(APPEND CGAL_DOC_PACKAGES "${package}")
    endif()
  endforeach()

  subdirlist(CGAL_EXAMPLE_PACKAGES ${CGAL_ROOT}/examples)
  subdirlist(CGAL_DEMO_PACKAGES ${CGAL_ROOT}/demo)
  list(APPEND CGAL_EXAMPLE_PACKAGES "BGL"
  )# manually add BGL to the list of packages
  foreach(pkg ${CGAL_EXAMPLE_PACKAGES})
    set(CGAL_${pkg}_EXAMPLE_DIR "${CGAL_ROOT}/examples")
    # This might not be entirely correct, but there is no reliable,
    # easy way to get a list of all package names when not in a
    # branch build.
    set(CGAL_${pkg}_INCLUDE_DIR "${CGAL_ROOT}/include")
  endforeach()
  foreach(pkg ${CGAL_DEMO_PACKAGES})
    set(CGAL_${pkg}_DEMO_DIR "${CGAL_ROOT}/demo")
  endforeach()
endif()

#special cases
foreach(pkg "Mesher_level")
  if(CGAL_BRANCH_BUILD)
    set(CGAL_${pkg}_INCLUDE_DIR "${CGAL_ROOT}/${pkg}/include")
  else()
    set(CGAL_${pkg}_INCLUDE_DIR "${CGAL_ROOT}/include")
  endif()
endforeach()

option(CGAL_BUILD_THREE_DOC "Build the documentation of the Three package" OFF)

if(NOT CGAL_BUILD_THREE_DOC)
  list(REMOVE_ITEM CGAL_DOC_PACKAGES "Three")
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


set(CGAL_DOC_TESTSUITE_SCRIPT "${CGAL_DOC_SCRIPT_DIR}/testsuite.py")

add_custom_target(
  doc_with_postprocessing
  ${Python3_EXECUTABLE} ${CGAL_DOC_SCRIPT_DIR}/html_output_post_processing.py
  --output ${CGAL_DOC_OUTPUT_DIR} --resources ${CGAL_DOC_RESOURCE_DIR})
add_dependencies(doc_with_postprocessing doc)

if(CGAL_DOC_CREATE_LOGS)
  add_custom_target(
    Documentation_test
    ${Python3_EXECUTABLE} ${CGAL_DOC_TESTSUITE_SCRIPT} --output-dir
    ${CGAL_DOC_OUTPUT_DIR} --doc-log-dir ${CGAL_DOC_LOG_DIR})
  add_dependencies(Documentation_test doc)

  add_custom_target(
    Documentation_test_publish
    ${Python3_EXECUTABLE}
    ${CGAL_DOC_TESTSUITE_SCRIPT}
    --output-dir
    ${CGAL_DOC_OUTPUT_DIR}
    --doc-log-dir
    ${CGAL_DOC_LOG_DIR}
    --publish
    ${CGAL_DOC_PUBLISH_DIR}
    --do-copy-results)
  add_dependencies(Documentation_test_publish doc)

  add_custom_target(
    doc_and_publish_testsuite
    ${Python3_EXECUTABLE}
    ${CGAL_DOC_TESTSUITE_SCRIPT}
    --output-dir
    ${CGAL_DOC_OUTPUT_DIR}
    --doc-log-dir
    ${CGAL_DOC_LOG_DIR}
    --publish
    ${CGAL_DOC_PUBLISH_DIR}
    --do-copy-results
    --cgal-version
    "${CGAL_FULL_VERSION}"
    --version-to-keep
    10)
  add_dependencies(doc_and_publish_testsuite doc_with_postprocessing)

endif()
