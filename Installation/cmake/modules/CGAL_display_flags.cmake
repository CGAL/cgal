if(CMAKE_CONFIGURATION_TYPES)
  # Multi-configuration CMake generator,
  message( STATUS "Multi-configuration CMake generator: cannot display flags" )
  return()
endif()

message( STATUS "Build type: ${CMAKE_BUILD_TYPE}" )

string( TOUPPER "${CMAKE_BUILD_TYPE}" CGAL_BUILD_TYPE_UPPER )

message( STATUS "USING CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CGAL_BUILD_TYPE_UPPER}}'" )

message( STATUS "USING EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${CGAL_BUILD_TYPE_UPPER}}'" )
