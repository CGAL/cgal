cmake_minimum_required(VERSION 3.1...3.15)
project(CGAL_ipelets_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

macro(remove_leading_zero var)
  string(SUBSTRING "${${var}}" 0 1 ONECHAR)
  string(COMPARE EQUAL "${ONECHAR}" "0" ISZERO)
  if(${ISZERO})
    string(SUBSTRING "${${var}}" 1 1 ONECHAR)
    set(${var} ${ONECHAR})
  endif()
endmacro()

#path where to build libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen_support)
if(NOT TARGET CGAL::Eigen_support)
  message(
    STATUS
      "NOTICE: This project requires the Eigen library, and will not be compiled."
  )
  return()
endif()

find_package(IPE 6)

if(IPE_FOUND)
  include_directories(BEFORE ${IPE_INCLUDE_DIR})

  if(${IPE_VERSION} EQUAL "7")
    set(WITH_IPE_7 ON)
  elseif(${IPE_VERSION} EQUAL "6")
    set(WITH_IPE_7 OFF)
  else()
    message(
      "-- Error: ${IPE_VERSION} is not a supported version of IPE (only 6 and 7 are)."
    )
    set(IPE_FOUND FALSE)
  endif()
  # starting ipe 7.2.1, a compiler with c++11 must be used to compile ipelets
  if(${IPE_VERSION} EQUAL "7"
     AND ${IPE_MINOR_VERSION_1} GREATER "1"
     AND ${IPE_MINOR_VERSION_2} GREATER "0")
    message(
      STATUS
        "Starting from Ipe 7.2.1 a compiler with c++11 support must be used")

    list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_generalized_initializers has_cpp11)
    if(has_cpp11 LESS 0)
      message(
        STATUS
          "NOTICE: This demo requires a C++11 compiler and will not be compiled."
      )
      return()
    endif()

    # Use C++11 for this directory and its sub-directories.
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

  endif()
endif()

if(IPE_FOUND AND IPE_VERSION)
  if(WITH_IPE_7)
    add_definitions(-DCGAL_USE_IPE_7)
  endif()

  message("-- Using IPE version ${IPE_VERSION} compatibility.")

  #setting installation directory
  get_filename_component(IPE_LIBRARY_DIR ${IPE_LIBRARIES} PATH)
  if(IPE_FOUND AND NOT IPELET_INSTALL_DIR)
    if(WITH_IPE_7)
      remove_leading_zero(IPE_MINOR_VERSION_1)
      remove_leading_zero(IPE_MINOR_VERSION_2)
      set(INSTALL_PATHS
          "${IPE_LIBRARY_DIR}/ipe/7.${IPE_MINOR_VERSION_1}.${IPE_MINOR_VERSION_2}/ipelets/"
      )
      find_path(
        IPELET_INSTALL_DIR
        NAMES libgoodies.lua goodies.lua
        PATHS ${INSTALL_PATHS}
        DOC "The folder where ipelets will be installed" ENV IPELETPATH)
    else()
      foreach(VER RANGE 28 40)
        string(REPLACE XX ${VER} PATHC
                       "${IPE_LIBRARY_DIR}/ipe/6.0preXX/ipelets/")
        set(INSTALL_PATHS ${INSTALL_PATHS} ${PATHC})
      endforeach()
      set(INSTALL_PATHS ${INSTALL_PATHS} ${PATHC})
      set(INSTALL_PATHS ${INSTALL_PATHS} /usr/lib64/ipe/6.0/ipelets)
      set(INSTALL_PATHS ${INSTALL_PATHS} /usr/lib/ipe/6.0/ipelets)

      find_library(
        IPELET_INSTALL_DIR_FILES
        NAMES align
        PATHS ${INSTALL_PATHS} ENV IPELETPATH)
      if(IPELET_INSTALL_DIR_FILES)
        get_filename_component(IPELET_INSTALL_DIR ${IPELET_INSTALL_DIR_FILES}
                               PATH)
      endif()
    endif()
  endif()

  set(CGAL_IPELETS ${CGAL_IPELETS})
  set(CGAL_IPELETS ${CGAL_IPELETS} alpha_shapes)
  set(CGAL_IPELETS ${CGAL_IPELETS} arrangement)
  set(CGAL_IPELETS ${CGAL_IPELETS} bbox_restriction)
  set(CGAL_IPELETS ${CGAL_IPELETS} diagrams)
  set(CGAL_IPELETS ${CGAL_IPELETS} hilbert_sort)
  set(CGAL_IPELETS ${CGAL_IPELETS} hull)
  set(CGAL_IPELETS ${CGAL_IPELETS} generator)
  set(CGAL_IPELETS ${CGAL_IPELETS} mesh_2)
  set(CGAL_IPELETS ${CGAL_IPELETS} minkowski)
  set(CGAL_IPELETS ${CGAL_IPELETS} mst)
  set(CGAL_IPELETS ${CGAL_IPELETS} multi_delaunay)
  set(CGAL_IPELETS ${CGAL_IPELETS} multi_regular)
  set(CGAL_IPELETS ${CGAL_IPELETS} partition)
  set(CGAL_IPELETS ${CGAL_IPELETS} pca)
  set(CGAL_IPELETS ${CGAL_IPELETS} skeleton)
  set(CGAL_IPELETS ${CGAL_IPELETS} svdlinf)
  set(CGAL_IPELETS ${CGAL_IPELETS} triangulation)
  set(CGAL_IPELETS ${CGAL_IPELETS} circle_pencils)
  set(CGAL_IPELETS ${CGAL_IPELETS} hyperbolic)
  set(CGAL_IPELETS ${CGAL_IPELETS} distance)
  if(CGAL_Core_FOUND)
    set(CGAL_IPELETS ${CGAL_IPELETS} cone_spanners)
  endif()

  if(IPELET_INSTALL_DIR)
    message(STATUS "Set Ipelets install dir: ${IPELET_INSTALL_DIR}")
  endif()

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)

  foreach(IPELET ${CGAL_IPELETS})
    add_library(CGAL_${IPELET} MODULE ${IPELET}.cpp)
    add_to_cached_list(CGAL_EXECUTABLE_TARGETS CGAL_${IPELET})
    target_link_libraries(CGAL_${IPELET} PRIVATE CGAL::CGAL CGAL::Eigen_support
                                                 ${IPE_LIBRARIES})
    if(IPELET_INSTALL_DIR)
      install(TARGETS CGAL_${IPELET} DESTINATION ${IPELET_INSTALL_DIR})
      if(WITH_IPE_7)
        install(FILES ./lua/libCGAL_${IPELET}.lua
                DESTINATION ${IPELET_INSTALL_DIR}) #only for ipe 7
      endif()
    endif()
    cgal_add_compilation_test(CGAL_${IPELET})
  endforeach(IPELET)
  if(CGAL_Core_FOUND)
    target_link_libraries(CGAL_cone_spanners PRIVATE CGAL::CGAL_Core
                                                     CGAL::Eigen_support)
  endif()
  #example in doc not installed
  add_library(simple_triangulation MODULE simple_triangulation.cpp)
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS simple_triangulation)
  target_link_libraries(simple_triangulation CGAL::Eigen_support
                        ${IPE_LIBRARIES})
  cgal_add_compilation_test(simple_triangulation)

else()
  message(
    STATUS
      "NOTICE: This program requires the Ipe include files and library, and will not be compiled."
  )
endif()
