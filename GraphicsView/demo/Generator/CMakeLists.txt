# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (GeneratorDemo)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL COMPONENTS Qt4)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(Qt4)

include_directories (BEFORE ../../include)

if ( CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND )

  include(${QT_USE_FILE})
  #--------------------------------
  # Demo: Generator_2
  #--------------------------------
  # UI files (Qt Designer files)
  qt4_wrap_ui( DT_UI_FILES Generator_2.ui )

  # qrc files (resources files, that contain icons, at least)
  qt4_add_resources ( DT_RESOURCE_FILES ./Generator_2.qrc )

  # use the Qt MOC preprocessor on classes that derives from QObject
  qt4_generate_moc( "Generator_2.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Generator_2.moc" )

  # The executable itself.
  add_executable  ( Generator_2 Generator_2.cpp Generator_2.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Generator_2 )

  # Link with Qt libraries
  target_link_libraries( Generator_2 ${QT_LIBRARIES} )
  # And with CGAL libraries
  target_link_libraries( Generator_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

else()

  message(STATUS "NOTICE: This demo requires CGAL, and Qt4, and will not be compiled.")

endif()
