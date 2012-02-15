project (CGAL_Qt4)

# Minimal version of CMake:
cmake_minimum_required(VERSION 2.6.2)

# Tested version:
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

if ( NOT CGAL_FOUND )   
  find_package(CGAL REQUIRED)
endif()

include(${CGAL_USE_FILE})

set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTOPENGL TRUE )
find_package(Qt4 QUIET)

find_package(OpenGL QUIET )

if( QT4_FOUND )

  if( OPENGL_FOUND )
    
    include_directories (BEFORE ../../include)
    
    include(${QT_USE_FILE})

    get_dependency_version(OPENGL OpenGL)
    
    get_dependency_version(QT4 QT)

    if(COMMAND add_config_flag)
      set( CGAL_HAS_QT4 TRUE ) 
      add_config_flag( CGAL_HAS_QT4 )
    endif()
    
    message( STATUS "OpenGL include:     ${OPENGL_INCLUDE_DIR}" )
    message( STATUS "OpenGL libraries:   ${OPENGL_LIBRARIES}" )
    message( STATUS "OpenGL definitions: ${OPENGL_DEFINITIONS}" )
    message( STATUS "Qt4 include:        ${QT_INCLUDE_DIR}" )
    message( STATUS "Qt4 libraries:      ${QT_LIBRARIES}" )
    message( STATUS "Qt4 definitions:    ${QT_DEFINITIONS}" )
    message( STATUS "moc executable:     ${QT_MOC_EXECUTABLE}" )
    message( STATUS "uic executable:     ${QT_UIC_EXECUTABLE}" )

    cache_set(CGAL_Qt4_3RD_PARTY_INCLUDE_DIRS ${QT_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR} )
    cache_set(CGAL_Qt4_3RD_PARTY_LIBRARIES    ${QT_LIBRARIES}   ${OPENGL_LIBRARIES}   )
    cache_set(CGAL_Qt4_3RD_PARTY_DEFINITIONS  ${QT_DEFINITIONS} ${OPENGL_DEFINITIONS} )

    set(mocfiles "")
    set(RESOURCE_FILES "")

    foreach (package ${CGAL_CONFIGURED_PACKAGES} )
      file(GLOB PACKAGE_QTMOC_FILES "${package}/src/CGALQt4/*.qtmoc.cmake")
      foreach(package_qtmoc_file ${PACKAGE_QTMOC_FILES})
         # includes 'moccing' for sources/headers in package + collects lists of moc-files for dependency (to properly build the lib)
        include(${package_qtmoc_file}) 
#        message(STATUS QTMOC------------FILE: ${package_qtmoc_file})      
      endforeach()
    endforeach()
  
    foreach(mocfile ${mocfiles})
      list(APPEND additional_files ${mocfile})
    endforeach()
    foreach(resfile ${RESOURCE_FILES})
      list(APPEND additional_files ${resfile})
    endforeach()

#    message(STATUS "Additional input files: ${additional_files}")

    build_cgal_library( CGAL_Qt4 CGALQt4 "${additional_files}")

    if ( NOT CGAL_IS_CURRENT_SCRIPT_TOP_LEVEL )
      add_dependencies( CGAL_Qt4 CGAL )
    endif()    

    target_link_libraries( CGAL_Qt4 ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_Qt4_3RD_PARTY_LIBRARIES} )

    add_definitions ( ${CGAL_3RD_PARTY_DEFINITIONS} )

    if($ENV{CGAL_FAKE_PUBLIC_RELEASE})
      add_definitions( -DCGAL_FAKE_PUBLIC_RELEASE )
    endif()
    
  else()
    message( STATUS "CGAL_Qt4 needs OpenGL, cannot be configured.")
  endif()  
  
else()
  message( STATUS "CGAL_Qt4 needs Qt4, cannot be configured.")
endif()
