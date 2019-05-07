# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Stream_support_Tests )

cmake_minimum_required(VERSION 3.1)

find_package(CGAL QUIET)
if ( CGAL_FOUND )
  set(3MF_INCLUDE_DIR "" CACHE PATH "Path to lib3MF headers")
  set(3MF_LIBRARY_DIR "" CACHE PATH "Path to lib3MF library files")
  if(IS_DIRECTORY "${3MF_INCLUDE_DIR}/Common")
    set(3MF_FOUND true)
  endif()
  
  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    if ( "${cppfile}" STREQUAL "test_3mf_to_sm.cpp" )
      
      if(3MF_FOUND)
        if (WIN32)
          set(LSUFFIX "dll")
          set(LSUFFIXOUT ".dll")
        elseif(UNIX AND NOT APPLE)
          set(LSUFFIX "so")
          set(LSUFFIXOUT ".so.1")
        else()
          set(LSUFFIX "dylib")
          set(LSUFFIXOUT ".1.dylib")
        endif()
        include_directories(${3MF_INCLUDE_DIR})
        link_directories(${3MF_LIBRARY_DIR})
        create_single_source_cgal_program( "${cppfile}" )
        
        if (WIN32)
          target_link_libraries(test_3mf_to_sm PRIVATE lib3MF)
        else()
          # Unix prefixes the name of the library with "lib" anyway
          target_link_libraries(test_3mf_to_sm PRIVATE 3MF)
          if(APPLE)
            target_link_directories(test_3mf_to_sm PRIVATE ${3MF_LIBRARY_DIR}/..)
          endif()
        endif()

      else()
        message(STATUS "NOTICE: This program requires the lib3MF library, and will not be compiled.")
      endif()
    else()
      create_single_source_cgal_program( "${cppfile}" )
    endif()
  endforeach()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()


