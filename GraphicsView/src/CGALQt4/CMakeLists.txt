# src/CGALQt4/CMakeLists.txt
#
# The following variables must be set from the parent script:
# The first two cannot be empty
#
# CGAL_BINARY_DIR
# CGAL_LIB_INSTALL_DIR
#
# CGAL_INCLUDE_DIRS
# GMP_INCLUDE_DIRS
# CGAL_DEFINITIONS
# CGAL_3RD_PARTY_INCLUDE_DIRS
# CGAL_3RD_PARTY_LIBRARIES
# CGAL_3RD_PARTY_DEFINITIONS
# CGAL_LIB_SHARED_LINKER_FLAGS
# CGAL_VERSION
# CGAL_MAJOR_VERSION
# CGAL_MINOR_VERSION
# CGAL_SOVERSION

project (CGALQt4)
   
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTOPENGL TRUE )
find_package(Qt4 QUIET)

include_directories (BEFORE ../../include)

if(QT4_FOUND) 
  include(${QT_USE_FILE})

  get_dependency_version(QT4 Qt4)

  # moc files that are compiled directly as cpp files
  qt4_generate_moc( ../../include/CGAL/Qt/GraphicsViewNavigation.h GraphicsViewNavigation_moc.cpp )
  qt4_generate_moc( ../../include/CGAL/Qt/DemosMainWindow.h DemosMainWindow_moc.cpp )
  qt4_generate_moc( ../../include/CGAL/Qt/GraphicsItem.h GraphicsItem_moc.cpp )
  qt4_generate_moc( ../../include/CGAL/Qt/GraphicsViewInput.h GraphicsViewInput_moc.cpp )

  # qrc files (resources files, that contain icons, at least)
  qt4_add_resources ( RESOURCE_FILES ../../demo/resources/CGAL.qrc ../../demo/icons/Input.qrc ../../demo/icons/File.qrc ../../demo/icons/Triangulation_2.qrc)

  set( CGAL_QT4_BASENAME CGAL-Qt4 )

  add_library( CGAL_QT4 debug.cpp GraphicsViewInput_moc.cpp GraphicsViewNavigation.cpp GraphicsViewNavigation_moc.cpp GraphicsItem_moc.cpp GraphicsViewPolylineInput.cpp DemosMainWindow.cpp DemosMainWindow_moc.cpp ${RESOURCE_FILES} )
  target_link_libraries( CGAL_QT4 ${QT_LIBRARIES} )
  if ( NOT AUTO_LINK_ENABLED )    
    # If auto-link is not enabled (on Linux for example), link with libCGAL explicitly.
    target_link_libraries( CGAL_QT4 CGAL )
  endif()

  set_target_properties(CGAL_QT4 PROPERTIES VERSION "${CGAL_VERSION}" SOVERSION "${CGAL_SOVERSION}")

  if (NOT AUTO_LINK_ENABLED )
    set_target_properties( CGAL_QT4 PROPERTIES OUTPUT_NAME "${CGAL_QT4_BASENAME}" )
    get_target_property( CGAL_QT4_LIBRARY_FULL_PATH CGAL_QT4 LOCATION  )

    set(CGAL_QT4_LIBRARY "${CGAL_QT4_LIBRARY_FULL_PATH}" CACHE STRING "The CGAL-ImageIO library")
    mark_as_advanced(CGAL_QT4_LIBRARY)

  else(NOT AUTO_LINK_ENABLED)
  
    set(CGAL_QT4_LIBRARY)
    
    set(CGAL_QT4_RELEASE_MANGLED_NAME        "${CGAL_QT4_BASENAME}-${TOOLSET}-mt" ) 
    set(CGAL_QT4_DEBUG_MANGLED_NAME          "${CGAL_QT4_BASENAME}-${TOOLSET}-mt-gd" ) 
    set(CGAL_QT4_MINSIZEREL_MANGLED_NAME     "${CGAL_QT4_BASENAME}-${TOOLSET}-mt-o" ) 
    set(CGAL_QT4_RELWITHDEBINFO_MANGLED_NAME "${CGAL_QT4_BASENAME}-${TOOLSET}-mt-g" ) 
    
    set_target_properties( CGAL_QT4 PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_QT4_DEBUG_MANGLED_NAME}" 
                                               RELEASE_OUTPUT_NAME        "${CGAL_QT4_RELEASE_MANGLED_NAME}" 
                                               MINSIZEREL_OUTPUT_NAME     "${CGAL_QT4_MINSIZEREL_MANGLED_NAME}" 
                                               RELWITHDEBINFO_OUTPUT_NAME "${CGAL_QT4_RELWITHDEBINFO_MANGLED_NAME}" 
                         )
    
    if ( HAS_CFG_INTDIR )                      
      add_custom_command(TARGET CGAL_QT4 POST_BUILD COMMAND copy /Y "$(TargetPath)" "$(SolutionDir)lib" )
      add_custom_command(TARGET CGAL_QT4 POST_BUILD COMMAND if exist "$(TargetDir)$(TargetName).pdb" copy /Y "$(TargetDir)$(TargetName).pdb" "$(SolutionDir)lib" )
    endif()

  endif(NOT AUTO_LINK_ENABLED)

endif(QT4_FOUND)
