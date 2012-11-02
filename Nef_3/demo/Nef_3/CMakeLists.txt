# This is the CMake script for compiling a CGAL application.

project (Nef_3_Demo)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL COMPONENTS Qt3)
include(${CGAL_USE_FILE})

find_package(Qt3-patched )
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".

macro(create_single_source_cgal_qt3_program first )

  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${first})
    set( all ${CMAKE_CURRENT_SOURCE_DIR}/${first} )
    
    foreach( i ${ARGN} )
      set( all ${all} ${CMAKE_CURRENT_SOURCE_DIR}/${i} ) 
    endforeach()
    
    get_filename_component(exe_name ${first} NAME_WE)
    
    QT3_AUTOMOC( ${all} )
    
    # Make sure the compiler can find generated .moc files
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
    include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
    
    include_directories( ${QT3_INCLUDE_DIR} )
    
    add_executable  (${exe_name} ${all})
  
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${exe_name} )
    
    # Link the executable to CGAL and third-party libraries
    if ( CGAL_AUTO_LINK_ENABLED )    
      target_link_libraries(${exe_name} ${CGAL_3RD_PARTY_LIBRARIES} ${QT3_LIBRARIES} )
    else()
      target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT3_LIBRARIES})
    endif()
    
  endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${first})
  
endmacro()

if ( CGAL_FOUND AND QT3_FOUND AND CGAL_Qt3_FOUND )

  include(Qt3Macros-patched)
  include(CGAL_CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_qt3_program( nef_3_extended.cpp    )
  create_single_source_cgal_qt3_program( nef_3_filtered.cpp    )
  create_single_source_cgal_qt3_program( nef_3_homogeneous.cpp )
  create_single_source_cgal_qt3_program( visual_hull.cpp       )
  create_single_source_cgal_qt3_program( visualization_SM.cpp  )
  create_single_source_cgal_qt3_program( visualization_SNC.cpp )

else( CGAL_FOUND AND QT3_FOUND AND CGAL_Qt3_FOUND )

  message(STATUS "NOTICE: This demo requires CGAL and Qt3, and will not be compiled.")

endif( CGAL_FOUND AND QT3_FOUND AND CGAL_Qt3_FOUND )
