# src/Core/CMakeLists.txt
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

project (CGALQt4)
   
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

set( QT_USE_QTMAIN   TRUE )
find_package(Qt4 QUIET)

include_directories (BEFORE ../../include)

if(QT4_FOUND) 
  include(${QT_USE_FILE})

  # a moc file that is included in a .cpp file
  qt4_generate_moc( ../../include/CGAL/IO/QtNavigation.h QtNavigation.moc )
  macro_add_file_dependencies( QtNavigation.cpp QtNavigation.moc )

  # two moc files that are compiled directly as cpp files
  qt4_generate_moc( ../../include/CGAL/IO/QtGraphicsItem.h QtGraphicsItem_moc.cpp )
  qt4_generate_moc( ../../include/CGAL/IO/QtInput.h QtInput_moc.cpp )

  # qrc files (resources files, that contain icons, at least)
  qt4_add_resources ( RESOURCE_FILES ../../demo/icons/Input.qrc ../../demo/icons/File.qrc ../../demo/icons/Logos.qrc ../../demo/icons/Triangulation_2.qrc)

  set( CGAL_QT4_BASENAME CGAL-Qt4 )

  add_library( CGAL_QT4 QtInput_moc.cpp QtNavigation.cpp QtGraphicsItem_moc.cpp QtPolylineInput.cpp ${RESOURCE_FILES} )
  target_link_libraries( CGAL_QT4 ${QT_LIBRARIES} )
  if ( NOT AUTO_LINK_ENABLED )    
    # If auto-link is not enabled (on Linux for example), link with libCGAL explicitly.
    target_link_libraries( CGAL_QT4 CGAL )
  endif()

  set_target_properties( CGAL_QT4 PROPERTIES OUTPUT_NAME "${CGAL_QT4_BASENAME}" )

endif(QT4_FOUND)
