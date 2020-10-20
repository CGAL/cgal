cmake_minimum_required(VERSION 3.1...3.15)
project( FastEnvelope_Examples )

set(FAST_ENVELOPE_BUILD_DIR "" CACHE PATH "Path to fast-evelope build directory")

find_package(CGAL REQUIRED)

find_package(Eigen3 REQUIRED 3.2.0) #(requires 3.2.0 or greater)
include(CGAL_Eigen_support)


create_single_source_cgal_program( "fast.cpp" )

if (FAST_ENVELOPE_BUILD_DIR)
  message(STATUS "Using ${FAST_ENVELOPE_BUILD_DIR} as build directory of fast-evelope")

  include_directories("${FAST_ENVELOPE_BUILD_DIR}/include")
  link_directories ( "${FAST_ENVELOPE_BUILD_DIR}/lib"  "${FAST_ENVELOPE_BUILD_DIR}"  "${FAST_ENVELOPE_BUILD_DIR}/tbb"   )

  find_package(OpenMP)
  if (OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()

  create_single_source_cgal_program( "fastE.cpp" )
  target_link_libraries( fastE PUBLIC CGAL::Eigen_support)
  target_link_libraries( fastE PUBLIC FastEnvelope IndirectPredicates geogram)

else()
  message(STATUS "Cmake variable FAST_ENVELOPE_BUILD_DIR is not defined fastE will not be built")
endif()
