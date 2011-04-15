project (CGAL_Qt3)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5 FATAL_ERROR)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
  
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

if ( NOT CGAL_FOUND )
  find_package(CGAL REQUIRED)
  include(${CGAL_USE_FILE})
endif()

find_package(OpenGL QUIET )

find_package(Qt3-patched QUIET )
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".

if( QT3_FOUND )

  if( OPENGL_FOUND )

    get_dependency_version(OPENGL)
    get_dependency_version(QT QT3)
    
    if(COMMAND add_config_flag)
      set( CGAL_HAS_QT3 TRUE ) 
      add_config_flag( CGAL_HAS_QT3 )
    endif()
    
    message( STATUS "OpenGL include:     ${OPENGL_INCLUDE_DIR}" )
    message( STATUS "OpenGL include:     ${OPENGL_INCLUDE_DIR}" )
    message( STATUS "OpenGL libraries:   ${OPENGL_LIBRARIES}" )
    message( STATUS "OpenGL definitions: ${OPENGL_DEFINITIONS}" )
    message( STATUS "Qt3 include:        ${QT3_INCLUDE_DIR}" )
    message( STATUS "Qt3 libraries:      ${QT3_LIBRARIES}" )
    message( STATUS "Qt3 definitions:    ${QT3_DEFINITIONS}" )
    message( STATUS "moc executable:     ${QT3_MOC_EXECUTABLE}" )
    message( STATUS "uic executable:     ${QT3_UIC_EXECUTABLE}" )
    
    cache_set(CGAL_Qt3_3RD_PARTY_INCLUDE_DIRS ${QT3_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR} )
    cache_set(CGAL_Qt3_3RD_PARTY_LIBRARIES    ${QT3_LIBRARIES}   ${OPENGL_LIBRARIES}   )
    cache_set(CGAL_Qt3_3RD_PARTY_DEFINITIONS  ${QT3_DEFINITIONS} ${OPENGL_DEFINITIONS} )
    
    include(Qt3Macros-patched)
       
    include_directories( ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_Qt3_3RD_PARTY_INCLUDE_DIRS} )
       
    # Make sure the compiler can find generated .moc files
    include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})

    # Make sure QT3_AUTOMOC compiler can find CGALQt include files
    include_directories (BEFORE ${CMAKE_SOURCE_DIR}/include/CGAL/IO)

    QT3_AUTOMOC( Qt_widget.cpp
                 Qt_widget_layer.cpp
                 Qt_widget_standard_toolbar.cpp
                 Qt_help_window.cpp
                 Qt_widget_history.cpp
                 Qt_widget_xpm_icons.cpp
                 Qt_widget_OpenGL.cpp
               )  

    QT3_GENERATE_MOC(${CMAKE_SOURCE_DIR}/include/CGAL/Kinetic/IO/internal/Qt_core.h
                     ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_core.moc
                     )
                          
    QT3_GENERATE_MOC(${CMAKE_SOURCE_DIR}/include/CGAL/Kinetic/IO/internal/Qt_timer.h
                     ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_timer.moc
                    )
                        
    QT3_GENERATE_MOC(${CMAKE_SOURCE_DIR}/include/CGAL/Kinetic/IO/internal/Qt_widget_2_core.h
                     ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_widget_2_core.moc
                    )
                        
    QT3_GENERATE_MOC(${CMAKE_SOURCE_DIR}/include/CGAL/Kinetic/IO/internal/Qt_window_2.h
                     ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_window_2.moc
                     )
                    
    add_file_dependencies(all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_core.moc )
    add_file_dependencies(all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_timer.moc )
    add_file_dependencies(all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_widget_2_core.moc )
    add_file_dependencies(all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_window_2.moc )
    add_file_dependencies(all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/Qt_widget.moc ) 
    add_file_dependencies(all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/Qt_widget_layer.moc )
    add_file_dependencies(all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/Qt_widget_standard_toolbar.moc )
    add_file_dependencies(all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/Qt_help_window.moc )
    add_file_dependencies(all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/Qt_widget_history.moc )
    add_file_dependencies(all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/Qt_widget_OpenGL.moc )

    add_library (CGAL_Qt3 "all_files.cpp")
    
    add_dependencies( CGAL_Qt3 CGAL )
    
    if(CGAL_SOVERSION AND CGAL_SONAME_VERSION)
      set_target_properties(CGAL_Qt3 PROPERTIES VERSION "${CGAL_SOVERSION}" SOVERSION "${CGAL_SONAME_VERSION}")
    endif()
    
    target_link_libraries( CGAL_Qt3 ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_Qt3_3RD_PARTY_LIBRARIES} )
        
    add_definitions ( ${CGAL_3RD_PARTY_DEFINITIONS} ${CGAL_Qt3_3RD_PARTY_DEFINITIONS} )
    
    set( CGAL_Qt3_BASENAME CGAL_Qt3 )
    
    if (NOT CGAL_AUTO_LINK_ENABLED )

      set_target_properties( CGAL_Qt3 PROPERTIES OUTPUT_NAME "${CGAL_Qt3_BASENAME}" )
      
      get_target_property( CGAL_Qt3_LIBRARY_FULL_PATH CGAL_Qt3 LOCATION  )
      
      set(CGAL_Qt3_LIBRARY "${CGAL_Qt3_LIBRARY_FULL_PATH}" CACHE STRING "The CGAL_Qt3 library" FORCE)
      
      mark_as_advanced(CGAL_Qt3_LIBRARY)
      
    else()

      set(CGAL_Qt3_LIBRARY)
      
      set(CGAL_Qt3_RELEASE_MANGLED_NAME        "${CGAL_Qt3_BASENAME}-${CGAL_TOOLSET}-mt" ) 
      set(CGAL_Qt3_DEBUG_MANGLED_NAME          "${CGAL_Qt3_BASENAME}-${CGAL_TOOLSET}-mt-gd" ) 
      set(CGAL_Qt3_MINSIZEREL_MANGLED_NAME     "${CGAL_Qt3_BASENAME}-${CGAL_TOOLSET}-mt-o" ) 
      set(CGAL_Qt3_RELWITHDEBINFO_MANGLED_NAME "${CGAL_Qt3_BASENAME}-${CGAL_TOOLSET}-mt-g" ) 
      
      set_target_properties( CGAL_Qt3 PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_Qt3_DEBUG_MANGLED_NAME}" 
                                                 RELEASE_OUTPUT_NAME        "${CGAL_Qt3_RELEASE_MANGLED_NAME}" 
                                                 MINSIZEREL_OUTPUT_NAME     "${CGAL_Qt3_MINSIZEREL_MANGLED_NAME}" 
                                                 RELWITHDEBINFO_OUTPUT_NAME "${CGAL_Qt3_RELWITHDEBINFO_MANGLED_NAME}" 
                           )

      if ( HAS_CFG_INTDIR )                      
        add_custom_command(TARGET CGAL_Qt3 POST_BUILD COMMAND copy /Y \"$(TargetPath)\" \"$(SolutionDir)lib\" )
        add_custom_command(TARGET CGAL_Qt3 POST_BUILD COMMAND if exist \"$(TargetDir)$(TargetName).pdb\" copy /Y \"$(TargetDir)$(TargetName).pdb\" \"$(SolutionDir)lib\" )
      endif()
      
    endif()

    if ( CGAL_INSTALL_LIB_DIR )
      install(TARGETS CGAL_Qt3 DESTINATION ${CGAL_INSTALL_LIB_DIR} )
    endif()

  else()
    message( STATUS "CGAL_Qt3 needs OpenGL, cannot be configured.")
  endif()  
else()
  message( STATUS "CGAL_Qt3 needs Qt3, cannot be configured.")
endif()

   

