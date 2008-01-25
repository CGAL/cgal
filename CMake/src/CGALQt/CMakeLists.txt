# src/CGALQt/CMakeLists.txt
#
# The following variables must be set from the parent script:
# The first two cannot be empty
#
# CGAL_BINARY_DIR
# CGAL_LIB_INSTALL_DIR
#
# CGAL_INCLUDE_DIRS
# CGAL_DEFINITIONS
# CGAL_3RD_PARTY_INCLUDE_DIRS
# CGAL_3RD_PARTY_LIBRARIES
# CGAL_3RD_PARTY_DEFINITIONS
# CGAL_LIB_SHARED_LINKER_FLAGS
# CGAL_VERSION
# CGAL_MAJOR_VERSION
# CGAL_MINOR_VERSION

find_package(Qt3)
if(QT_FOUND)

    include(Qt3Macros)

    add_definitions( "-DCGAL_USE_QT" )    
    
    assert(CGAL_BINARY_DIR)
    assert(CGAL_LIB_INSTALL_DIR)
       
    # Make sure the compiler can find generated .moc files
    include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})

    # Make sure QT3_AUTOMOC compiler can find CGALQt include files
    include_directories (BEFORE ${CGAL_SOURCE_DIR}/include/CGAL/IO)

    QT3_AUTOMOC( Qt_widget.cpp
                 Qt_widget_layer.cpp
                 Qt_widget_standard_toolbar.cpp
                 Qt_help_window.cpp
                 Qt_widget_history.cpp
                 Qt_widget_xpm_icons.cpp
                 Qt_widget_OpenGL.cpp
               )  

    QT3_GENERATE_MOC(${CGAL_SOURCE_DIR}/include/CGAL/Kinetic/IO/internal/Qt_core.h
                     ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_core.moc
                     )
                          
    QT3_GENERATE_MOC(${CGAL_SOURCE_DIR}/include/CGAL/Kinetic/IO/internal/Qt_timer.h
                     ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_timer.moc
                    )
                        
    QT3_GENERATE_MOC(${CGAL_SOURCE_DIR}/include/CGAL/Kinetic/IO/internal/Qt_widget_2_core.h
                     ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_widget_2_core.moc
                    )
                        
    QT3_GENERATE_MOC(${CGAL_SOURCE_DIR}/include/CGAL/Kinetic/IO/internal/Qt_window_2.h
                     ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_window_2.moc
                     )
                    
    add_file_dependencies(Kinetic_Qt_core.cpp          ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_core.moc )
    add_file_dependencies(Kinetic_Qt_timer.cpp         ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_timer.moc )
    add_file_dependencies(Kinetic_Qt_widget_2_core.cpp ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_widget_2_core.moc )
    add_file_dependencies(Kinetic_Qt_window_2.cpp      ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_window_2.moc )

    add_library (CGAL_QT "all_files.cpp")

    if (NOT AUTO_LINK_ENABLED )

      set_target_properties( CGAL_QT PROPERTIES OUTPUT_NAME "CGALQt" )
      
    else()

      set(CGAL_QT_RELEASE_MANGLED_NAME        "cgalqt${TOOLSET}-mt" ) 
      set(CGAL_QT_DEBUG_MANGLED_NAME          "cgalqt${TOOLSET}-mt-gd" ) 
      set(CGAL_QT_MINSIZEREL_MANGLED_NAME     "cgalqt${TOOLSET}-mt-o" ) 
      set(CGAL_QT_RELWITHDEBINFO_MANGLED_NAME "cgalqt${TOOLSET}-mt-g" ) 
      
      set_target_properties( CGAL PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_QT_DEBUG_MANGLED_NAME}" 
                                             RELEASE_OUTPUT_NAME        "${CGAL_QT_RELEASE_MANGLED_NAME}" 
                                             MINSIZEREL_OUTPUT_NAME     "${CGAL_QT_MINSIZEREL_MANGLED_NAME}" 
                                             RELWITHDEBINFO_OUTPUT_NAME "${CGAL_QT_RELWITHDEBINFO_MANGLED_NAME}" 
                           )
      
      add_custom_command(TARGET CGAL POST_BUILD COMMAND move /Y $(TargetDir)$(TargetName)$(TargetExt) ../../lib)
      add_custom_command(TARGET CGAL POST_BUILD COMMAND if exist "$(TargetDir)$(TargetName).pdb." move /Y $(TargetDir)$(TargetName).pdb ../../lib)
      
    endif()

    install(TARGETS CGAL_QT DESTINATION ${CGAL_LIB_INSTALL_DIR} )
endif()

   

