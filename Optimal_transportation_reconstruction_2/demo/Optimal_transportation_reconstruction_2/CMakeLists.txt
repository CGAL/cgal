# This is the CMake script for compiling the Optimal_transportation_reconstruction_2 demo.

cmake_minimum_required(VERSION 3.12...3.29)
project(Optimal_transportation_reconstruction_2_Demo)

# Find CGAL and CGAL Qt6
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

# Find Qt6 itself
find_package(Qt6 QUIET COMPONENTS Widgets)

# Find CImg
find_path(
  CIMG_INCLUDE_DIR
  NAMES CImg.h
  HINTS ENV CIMG_INC_DIR
  DOC "Path to the header of the CImg library")

if(CIMG_INCLUDE_DIR)
  message(STATUS "NOTICE: CImg library found, the demo can load point set from image files.")
else()
  message(STATUS "CImg library was not found, the demo will not be able to load point set from image files. "
                 "Try setting the environment variable CIMG_INC_DIR to point to the path of the directory containing CImg.h.")
endif()

if(CGAL_Qt6_FOUND AND Qt6_FOUND)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)

  add_executable(Otr2_demo glviewer.cpp scene.cpp Otr2_demo.cpp window.cpp render.cpp dialog_options.cpp
                           pwsrec.ui options.ui pwsrec.qrc)

  target_link_libraries(Otr2_demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt6)

  # Link with pthread if necessary
  if(CIMG_INCLUDE_DIR)
    target_compile_definitions(Otr2_demo PRIVATE -DCGAL_USE_CIMG)
    target_include_directories(Otr2_demo PRIVATE ${CIMG_INCLUDE_DIR})
    # Is pthread around? If yes, we need to link against it
    set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
    find_package(Threads QUIET)
    if(CMAKE_USE_PTHREADS_INIT)
      target_link_libraries(Otr2_demo PRIVATE ${CMAKE_THREAD_LIBS_INIT})
    endif()
  endif()

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Otr2_demo)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Otr2_demo)
else()

  set(OTR2_MISSING_DEPS "")

  if(NOT CGAL_Qt6_FOUND)
    set(OTR2_MISSING_DEPS "the CGAL Qt6 library, ${OTR2_MISSING_DEPS}")
  endif()

  if(NOT Qt6_FOUND)
    set(OTR2_MISSING_DEPS "Qt6, ${OTR2_MISSING_DEPS}")
  endif()

  message("NOTICE: This demo requires ${OTR2_MISSING_DEPS} and will not be compiled.")

endif()
