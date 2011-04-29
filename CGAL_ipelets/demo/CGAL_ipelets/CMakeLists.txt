cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()


#path where to build libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

find_package(CGAL QUIET COMPONENTS Core )
set( IPE_VERSION "AUTODETECT" CACHE STRING "The version of Ipe. Only 6 and 7 are supported.")

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  find_package(IPE)
  if ( IPE_FOUND )
    include_directories(BEFORE ../../include)
    include_directories(include/CGAL_ipelets)
    include_directories(${IPE_INCLUDE_DIR})
    link_directories( ${IPE_LIBRARY_DIR} )
    
    if (WITH_IPE_7)
      add_definitions(-DCGAL_USE_IPE_7)
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
    set(CGAL_IPELETS ${CGAL_IPELETS} multi_delaunay)
    set(CGAL_IPELETS ${CGAL_IPELETS} multi_regular)
    set(CGAL_IPELETS ${CGAL_IPELETS} partition)
    set(CGAL_IPELETS ${CGAL_IPELETS} pca)    
    set(CGAL_IPELETS ${CGAL_IPELETS} skeleton)
    set(CGAL_IPELETS ${CGAL_IPELETS} triangulation)


    foreach(IPELET ${CGAL_IPELETS})
      add_library(CGAL_${IPELET} MODULE ${IPELET}.cpp)
      add_to_cached_list(CGAL_EXECUTABLE_TARGETS CGAL_${IPELET})
      target_link_libraries( CGAL_${IPELET} ipe)
      if ( IPELET_INSTALL_DIR )
        install(TARGETS CGAL_${IPELET} DESTINATION ${IPELET_INSTALL_DIR})
        if (WITH_IPE_7)
          install(FILES ./lua/libCGAL_${IPELET}.lua DESTINATION ${IPELET_INSTALL_DIR}) #only for ipe 7
        endif()
      endif ()
    endforeach(IPELET)
    #example in doc not installed
    add_library(simple_triangulation MODULE simple_triangulation.cpp)
    target_link_libraries( simple_triangulation ipe )
    add_to_cached_list(CGAL_EXECUTABLE_TARGETS simple_triangulation)

  else()
    message(STATUS "NOTICE: This program requires the Ipe include files and library, and will not be compiled.")
  endif()
else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
    
endif()
