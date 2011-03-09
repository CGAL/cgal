project (CGAL_Qt4)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5 FATAL_ERROR)
   
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
  
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

if ( NOT CGAL_FOUND )   
  find_package(CGAL REQUIRED)
  include(${CGAL_USE_FILE})
endif()

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

      # moc files that are compiled directly as cpp files
    qt4_wrap_cpp(mocfiles ../../include/CGAL/Qt/GraphicsViewNavigation.h 
                          ../../include/CGAL/Qt/DemosMainWindow.h
                          ../../include/CGAL/Qt/GraphicsItem.h
                          ../../include/CGAL/Qt/GraphicsViewInput.h)

    # qrc files (resources files, that contain icons, at least)
    qt4_add_resources ( RESOURCE_FILES ../../demo/resources/CGAL.qrc ../../demo/icons/Input.qrc ../../demo/icons/File.qrc ../../demo/icons/Triangulation_2.qrc)

    set( CGAL_Qt4_BASENAME CGAL_Qt4 )

    add_library( CGAL_Qt4 debug.cpp 
			  GraphicsViewNavigation.cpp
			  GraphicsViewPolylineInput.cpp
			  DemosMainWindow.cpp
			  utility.cpp
                          ${mocfiles} 
			  ${RESOURCE_FILES} )
    
    if ( NOT CGAL_IS_CURRENT_SCRIPT_TOP_LEVEL )
      add_dependencies( CGAL_Qt4 CGAL )
    endif()    

    target_link_libraries( CGAL_Qt4 ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_Qt4_3RD_PARTY_LIBRARIES} )

    add_definitions ( ${CGAL_3RD_PARTY_DEFINITIONS} )

    if($ENV{CGAL_FAKE_PUBLIC_RELEASE})
      add_definitions( -DCGAL_FAKE_PUBLIC_RELEASE )
    endif()
    
    if(CGAL_SOVERSION AND CGAL_SONAME_VERSION)
      set_target_properties(CGAL_Qt4 PROPERTIES VERSION "${CGAL_SOVERSION}" SOVERSION "${CGAL_SONAME_VERSION}")
    endif()

    if (NOT CGAL_AUTO_LINK_ENABLED )
    
      set_target_properties( CGAL_Qt4 PROPERTIES OUTPUT_NAME "${CGAL_Qt4_BASENAME}" )
      get_target_property( CGAL_Qt4_LIBRARY_FULL_PATH CGAL_Qt4 LOCATION  )

      set(CGAL_Qt4_LIBRARY "${CGAL_Qt4_LIBRARY_FULL_PATH}" CACHE STRING "The CGAL_Qt4 library" FORCE)
      mark_as_advanced(CGAL_Qt4_LIBRARY)

    else()
    
      set(CGAL_Qt4_LIBRARY)
      
      set(CGAL_Qt4_RELEASE_MANGLED_NAME        "${CGAL_Qt4_BASENAME}-${CGAL_TOOLSET}-mt" ) 
      set(CGAL_Qt4_DEBUG_MANGLED_NAME          "${CGAL_Qt4_BASENAME}-${CGAL_TOOLSET}-mt-gd" ) 
      set(CGAL_Qt4_MINSIZEREL_MANGLED_NAME     "${CGAL_Qt4_BASENAME}-${CGAL_TOOLSET}-mt-o" ) 
      set(CGAL_Qt4_RELWITHDEBINFO_MANGLED_NAME "${CGAL_Qt4_BASENAME}-${CGAL_TOOLSET}-mt-g" ) 
      
      set_target_properties( CGAL_Qt4 PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_Qt4_DEBUG_MANGLED_NAME}" 
                                                 RELEASE_OUTPUT_NAME        "${CGAL_Qt4_RELEASE_MANGLED_NAME}" 
                                                 MINSIZEREL_OUTPUT_NAME     "${CGAL_Qt4_MINSIZEREL_MANGLED_NAME}" 
                                                 RELWITHDEBINFO_OUTPUT_NAME "${CGAL_Qt4_RELWITHDEBINFO_MANGLED_NAME}" 
                           )
      
      if ( HAS_CFG_INTDIR )                      
        add_custom_command(TARGET CGAL_Qt4 POST_BUILD COMMAND copy /Y \"$(TargetPath)\" \"$(SolutionDir)lib\" )
        add_custom_command(TARGET CGAL_Qt4 POST_BUILD COMMAND if exist \"$(TargetDir)$(TargetName).pdb\" copy /Y \"$(TargetDir)$(TargetName).pdb\" \"$(SolutionDir)lib\" )
      endif()

    endif()

    if ( CGAL_INSTALL_LIB_DIR )
      install(TARGETS CGAL_Qt4 DESTINATION ${CGAL_INSTALL_LIB_DIR} )
    endif()

  else()
    message( STATUS "CGAL_Qt4 needs OpenGL, cannot be configured.")
  endif()  
  
else()
  message( STATUS "CGAL_Qt4 needs Qt4, cannot be configured.")
endif()
