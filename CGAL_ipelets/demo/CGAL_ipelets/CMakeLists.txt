cmake_minimum_required(VERSION 3.1...3.20)
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
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message(
    STATUS
      "NOTICE: This project requires the Eigen library, and will not be compiled."
  )
  return()
endif()

find_package(IPE 7)

if(IPE_FOUND)
  if ( NOT ${IPE_VERSION} EQUAL "7")
    message("-- Error: ${IPE_VERSION} is not a supported version of IPE (only 7 is).")
    set(IPE_FOUND FALSE)
  endif()
endif()


if(IPE_FOUND AND IPE_VERSION)
  message("-- Using IPE version ${IPE_VERSION} compatibility.")

  #setting installation directory
  get_filename_component(IPE_LIBRARY_DIR ${IPE_LIBRARIES} PATH)
  if (IPE_FOUND AND NOT IPELET_INSTALL_DIR)
      remove_leading_zero(IPE_MINOR_VERSION_1)
      remove_leading_zero(IPE_MINOR_VERSION_2)
      set(INSTALL_PATHS "${IPE_LIBRARY_DIR}/ipe/7.${IPE_MINOR_VERSION_1}.${IPE_MINOR_VERSION_2}/ipelets/")
      find_path(IPELET_INSTALL_DIR
                    NAMES libgoodies.lua goodies.lua
                    PATHS ${INSTALL_PATHS}
                    DOC "The folder where ipelets will be installed"
                    ENV IPELETPATH
                   )
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
    target_include_directories(CGAL_${IPELET} BEFORE PRIVATE ${IPE_INCLUDE_DIR})

    add_to_cached_list(CGAL_EXECUTABLE_TARGETS CGAL_${IPELET})
    target_link_libraries(CGAL_${IPELET} PRIVATE CGAL::CGAL CGAL::Eigen3_support
                                                 ${IPE_LIBRARIES})
    if(IPELET_INSTALL_DIR)
      install(TARGETS CGAL_${IPELET} DESTINATION ${IPELET_INSTALL_DIR})
      install(FILES ./lua/libCGAL_${IPELET}.lua DESTINATION ${IPELET_INSTALL_DIR}) #only for ipe 7
    endif ()
    cgal_add_compilation_test(CGAL_${IPELET})
  endforeach(IPELET)
  if(CGAL_Core_FOUND)
    target_link_libraries(CGAL_cone_spanners PRIVATE CGAL::CGAL_Core
                                                     CGAL::Eigen3_support)
  endif()
  #example in doc not installed
  add_library(simple_triangulation MODULE simple_triangulation.cpp)
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS simple_triangulation)
  target_link_libraries(simple_triangulation CGAL::Eigen3_support
                        ${IPE_LIBRARIES})
  target_include_directories(simple_triangulation BEFORE PRIVATE ${IPE_INCLUDE_DIR})
  if (WITH_IPE_7)
      target_compile_definitions(simple_triangulation PRIVATE CGAL_USE_IPE_7)
    endif()
  cgal_add_compilation_test(simple_triangulation)

else()
  message(STATUS "NOTICE: This program requires the Ipe include files and library, and will not be compiled.")
endif()
