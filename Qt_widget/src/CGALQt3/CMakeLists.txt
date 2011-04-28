project (CGAL_Qt3)

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

    # Make sure QT3_AUTOMOC compiler can find CGALQt3 include files
    include_directories( BEFORE ${CGAL_INCLUDE_DIRS} )
       
    # Make sure the compiler can find generated .moc files
    include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})

    # needed for dependency (collected in package-specific CMakeLists.${packagename}.qtmoc)
    set(CGAL_Qt3_MOC_FILES "")  

    foreach (package ${CGAL_CONFIGURED_PACKAGES} )
      file(GLOB PACKAGE_QTMOC_FILES "${package}/src/CGALQt3/*.qtmoc.cmake")
      foreach(package_qtmoc_file ${PACKAGE_QTMOC_FILES})
         # includes 'moccing' for sources/headers in package + collects lists of moc-files for dependency (to properly build the lib)
        include(${package_qtmoc_file}) 
#        message(STATUS QTMOC------------FILE: ${package_qtmoc_file})      
      endforeach()
    endforeach()
  
#    message(STATUS QT3-mocfiles: ${CGAL_Qt3_MOC_FILES})

    build_cgal_library( CGAL_Qt3 CGALQt3 "${CGAL_Qt3_MOC_FILES}")

    foreach(mocfile ${CGAL_Qt3_MOC_FILES})
      add_file_dependencies(${CMAKE_CURRENT_BINARY_DIR}/all_files.cpp ${CMAKE_CURRENT_BINARY_DIR}/${mocfile})
    endforeach()

    add_dependencies( CGAL_Qt3 CGAL )

    if (CGAL_SOVERSION AND CGAL_SONAME_VERSION)
      set_target_properties(CGAL_Qt3 PROPERTIES VERSION "${CGAL_SOVERSION}" SOVERSION "${CGAL_SONAME_VERSION}")
    endif()
    
    target_link_libraries( CGAL_Qt3 ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_Qt3_3RD_PARTY_LIBRARIES} )
        
    add_definitions( ${CGAL_3RD_PARTY_DEFINITIONS} ${CGAL_Qt3_3RD_PARTY_DEFINITIONS} )
    
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

