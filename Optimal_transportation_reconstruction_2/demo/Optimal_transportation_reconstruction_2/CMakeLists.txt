# This is the CMake script for compiling the Optimal_transportation_reconstruction_2 demo.

cmake_minimum_required(VERSION 3.1...3.15)
project(Optimal_transportation_reconstruction_2_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Include this package's headers first
include_directories(BEFORE ./ ./include)

# Find CGAL and CGAL Qt5
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

# Find Qt5 itself
find_package(Qt5 5.4 QUIET COMPONENTS OpenGL)
if(Qt5_FOUND)
  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif(Qt5_FOUND)

# Find CImg
find_path(
  CIMG_INCLUDE_DIR
  NAMES CImg.h
  HINTS ENV CIMG_INC_DIR
  DOC "Path to the header of the CImg library")

if(CIMG_INCLUDE_DIR)
  message(
    STATUS "CImg library found, the demo can load point set from image files.")
else()
  message(
    STATUS
      "CImg library was not found, the demo will not be able to load point set from image files. "
      "Try setting the environment variable CIMG_INC_DIR to point to the path of the directory containing CImg.h."
  )
endif()

if(CGAL_Qt5_FOUND AND Qt5_FOUND)

  set(SRCS glviewer.cpp scene.cpp Otr2_demo.cpp window.cpp render.cpp)

  set(CGAL_Qt5_MOC_FILES moc_dialog_options.cxx moc_glviewer.cxx moc_window.cxx)

  set(UIS pwsrec.ui options.ui)

  qt5_wrap_ui(UI_FILES ${UIS})

  include(AddFileDependencies)

  qt5_generate_moc("window.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_window.cxx")
  add_file_dependencies(moc_window.cxx "${CMAKE_CURRENT_SOURCE_DIR}/window.h")

  qt5_generate_moc("glviewer.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_glviewer.cxx")
  add_file_dependencies(moc_glviewer.cxx
                        "${CMAKE_CURRENT_SOURCE_DIR}/glviewer.h")

  qt5_generate_moc("dialog_options.h"
                   "${CMAKE_CURRENT_BINARY_DIR}/moc_dialog_options.cxx")
  add_file_dependencies(moc_dialog_options.cxx
                        "${CMAKE_CURRENT_SOURCE_DIR}/dialog_options.h")

  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES pwsrec.qrc)

  add_executable(Otr2_demo ${SRCS} ${CGAL_Qt5_MOC_FILES} ${UI_FILES}
                           ${CGAL_Qt5_RESOURCE_FILES})

  target_link_libraries(Otr2_demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Gui
                                          Qt5::OpenGL)

  # Link with pthread if necessary
  if(CIMG_INCLUDE_DIR)
    target_compile_definitions(Otr2_demo PRIVATE -DCGAL_USE_CIMG)
    target_include_directories(Otr2_demo PRIVATE ${CIMG_INCLUDE_DIR})
    # Is pthread around? If yes, we need to link against it
    set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
    find_package(Threads)
    if(CMAKE_USE_PTHREADS_INIT)
      target_link_libraries(Otr2_demo PRIVATE ${CMAKE_THREAD_LIBS_INIT})
    endif()
  endif()

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Otr2_demo)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Otr2_demo)
else(
  CGAL_Qt5_FOUND
  AND Qt5_FOUND)

  set(OTR2_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(OTR2_MISSING_DEPS "the CGAL Qt5 library, ${OTR2_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(OTR2_MISSING_DEPS "Qt5.4, ${OTR2_MISSING_DEPS}")
  endif()

  message(
    STATUS
      "NOTICE: This demo requires ${OTR2_MISSING_DEPS}and will not be compiled."
  )

endif(
  CGAL_Qt5_FOUND
  AND Qt5_FOUND)
