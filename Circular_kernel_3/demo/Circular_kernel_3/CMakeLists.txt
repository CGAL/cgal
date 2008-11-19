# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Circular_kernel_3_demo )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Core Qt3 )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  
  find_package(Qt3-patched )
  # FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
  # that it can be used together with FindQt4: all its variables are prefixed
  # by "QT3_" instead of "QT_".
  
  if(CGAL_Qt3_FOUND AND QT3_FOUND)
  
    include( Qt3Macros-patched )
    qt3_automoc(  Circular_kernel_3.cpp )

    # Make sure the compiler can find generated .moc files
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  
    include_directories( ${QT3_INCLUDE_DIR} )

    add_executable  (Circular_kernel_3  Circular_kernel_3.cpp)
  
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Circular_kernel_3 )
  
  
    # Link the executable to CGAL and third-party libraries
    target_link_libraries(Circular_kernel_3 ${QT3_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  else()
  
    message(STATUS "NOTICE: This demo requires Qt3 and the CGAL Qt3 library, and will not be compiled.")
  
  endif()
  
else()
  
    message(STATUS "NOTICE: This demo requires the CGAL library, and will not be compiled.")
  
endif()

