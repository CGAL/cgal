message("Configuring libCGAL_Qt5")

if(WIN32)
message("Qt5 on Windows needs Windows SDK.")
set(CMAKE_LIBRARY_PATH "C:\\Program Files (x86)\\Windows Kits\\8.1\\Lib\\winv6.3\\um\\x64")
endif()

set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTOPENGL TRUE )

find_package(OpenGL QUIET )

cmake_policy(SET CMP0020 NEW)

FIND_PACKAGE(Qt5Core QUIET)
FIND_PACKAGE(Qt5Gui QUIET)
FIND_PACKAGE(Qt5Widgets QUIET)
FIND_PACKAGE(Qt5OpenGL QUIET)

list(APPEND QT_INCLUDE_DIR ${Qt5Core_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5OpenGL_INCLUDE_DIRS}) 
list(APPEND QT_LIBRARIES ${Qt5Core_LIBRARIES} ${Qt5Gui_LIBRARIES} ${Qt5Widgets_LIBRARIES} ${Qt5OpenGL_LIBRARIES})
list(APPEND QT_DEFINITIONS ${Qt5Core_DEFINITIONS} ${Qt5Gui_DEFINITIONS} ${Qt5Widgets_DEFINITIONS} ${Qt5OpenGL_DEFINITIONS})



if( Qt5Widgets_FOUND AND Qt5Gui_FOUND AND Qt5Widgets_FOUND AND Qt5OpenGL_FOUND )

  if( OPENGL_FOUND )
  
    include_directories (BEFORE ../../include)
	
	#Marche pas car vide
    #include(${QT_USE_FILE})

    get_dependency_version(OPENGL)
    
    get_dependency_version(QT5 QT)

    if(COMMAND add_config_flag)
      set( CGAL_HAS_QT5 TRUE ) 
      add_config_flag( CGAL_HAS_QT5 )
    endif()
    
    message( STATUS "OpenGL include:     ${OPENGL_INCLUDE_DIR}" )
    message( STATUS "OpenGL libraries:   ${OPENGL_LIBRARIES}" )
    message( STATUS "OpenGL definitions: ${OPENGL_DEFINITIONS}" )
    message( STATUS "Qt5 include:        ${QT_INCLUDE_DIR}" )
    message( STATUS "Qt5 libraries:      ${QT_LIBRARIES}")
    message( STATUS "Qt5 definitions:    ${QT_DEFINITIONS}" )
    message( STATUS "moc executable:     ${QT_MOC_EXECUTABLE}" )
    message( STATUS "uic executable:     ${QT_UIC_EXECUTABLE}" )

    cache_set(CGAL_Qt5_3RD_PARTY_INCLUDE_DIRS ${QT_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR} )
    cache_set(CGAL_Qt5_3RD_PARTY_LIBRARIES    ${QT_LIBRARIES}   ${OPENGL_LIBRARIES}   )
    cache_set(CGAL_Qt5_3RD_PARTY_DEFINITIONS  ${QT_DEFINITIONS} ${OPENGL_DEFINITIONS} -DNOMINMAX)

    set(mocfiles "")
    set(RESOURCE_FILES "")

    foreach (package ${CGAL_CONFIGURED_PACKAGES} )
      file(GLOB PACKAGE_QTMOC_FILES "${package}/src/CGAL_Qt5/*.qtmoc.cmake")  
      foreach(package_qtmoc_file ${PACKAGE_QTMOC_FILES})
         # includes 'moccing' for sources/headers in package + collects lists of moc-files for dependency (to properly build the lib)
        include(${package_qtmoc_file}) 
        #message(STATUS QTMOC------------FILE: ${package_qtmoc_file})      
      endforeach()
    endforeach()
	
    foreach(mocfile ${mocfiles})
      list(APPEND additional_files ${mocfile})
    endforeach()
    foreach(resfile ${RESOURCE_FILES})
      list(APPEND additional_files ${resfile})
    endforeach()

#    message(STATUS "Additional input files: ${additional_files}")

    use_essential_libs()

    include_directories( SYSTEM ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_Qt5_3RD_PARTY_INCLUDE_DIRS} )

    link_directories    ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )

    collect_cgal_library( CGAL_Qt5 "${additional_files}")

    add_dependencies( CGAL_Qt5 CGAL )
	
    target_link_libraries( CGAL_Qt5 ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_Qt5_3RD_PARTY_LIBRARIES} )

    add_definitions ( ${CGAL_3RD_PARTY_DEFINITIONS} ${CGAL_Qt5_3RD_PARTY_DEFINITIONS} )
	
    if($ENV{CGAL_FAKE_PUBLIC_RELEASE})
      add_definitions( -DCGAL_FAKE_PUBLIC_RELEASE )
    endif()
    
    message("libCGAL_Qt5 is configured")

  else()
    message( STATUS "libCGAL_Qt5 needs OpenGL, cannot be configured.")
  endif()  
  
else()
  message( STATUS "libCGAL_Qt5 needs Qt5, cannot be configured.")
endif()
