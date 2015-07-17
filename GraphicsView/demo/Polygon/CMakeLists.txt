# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (GraphicsView)

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL COMPONENTS Qt5 Core)

include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)

include_directories (BEFORE ../../include)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )

  add_definitions(-DQT_NO_KEYWORDS)

  if( CGAL_Core_FOUND)
    add_definitions(-DCGAL_USE_CORE)
  endif()

  #--------------------------------
  # Demo: Polygon_2
  #--------------------------------
  # UI files (Qt Designer files)
  qt5_wrap_ui( DT_UI_FILES Polygon_2.ui )

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( DT_RESOURCE_FILES ./Polygon_2.qrc )

  # use the Qt MOC preprocessor on classes that derives from QObject
  qt5_generate_moc( Polygon_2.cpp "${CMAKE_CURRENT_BINARY_DIR}/Polygon_2.moc" )

#  add_library( CGAL SHARED IMPORTED ) 
#  SET_PROPERTY(TARGET CGAL PROPERTY IMPORTED_LOCATION ${CGAL_LIBRARY} )

  # The executable itself.
  add_executable  ( Polygon_2 Polygon_2.cpp Polygon_2.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

  qt5_use_modules(Polygon_2 Xml Script OpenGL Svg)

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polygon_2 )


  # Link with Qt libraries
  target_link_libraries( Polygon_2 ${QT_LIBRARIES} )
  # And with CGAL libraries
  target_link_libraries( Polygon_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

else()

  message(STATUS "NOTICE: This demo requires CGAL, CGAL_Core, and Qt5, and will not be compiled.")

endif()
