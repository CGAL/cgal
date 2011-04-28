# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Envelope_3_demo )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core Qt3 )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  
  find_package(Qt3-patched QUIET )
  # FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
  # that it can be used together with FindQt4: all its variables are prefixed
  # by "QT3_" instead of "QT_".
  
  if(CGAL_Qt3_FOUND AND QT3_FOUND)
  
    include( Qt3Macros-patched )
    qt3_automoc(  envelope_3.cpp )

    # Make sure the compiler can find generated .moc files
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  
    include_directories( ${QT3_INCLUDE_DIR} )

  include_directories (BEFORE ../../include)

    add_executable  (envelope_3  envelope_3.cpp)
  
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS envelope_3 )
  
  
    # Link the executable to CGAL and third-party libraries
    target_link_libraries(envelope_3 ${QT3_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  else()
  
    message(STATUS "NOTICE: This demo requires Qt3 and the CGAL Qt3 library, and will not be compiled.")
  
  endif()
  
else()
  
    message(STATUS "NOTICE: This demo requires the CGAL library, and will not be compiled.")
  
endif()

