# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Largest_empty_rect_2_Demo)

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

find_package(CGAL COMPONENTS Qt5)

include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)

include_directories (BEFORE ../../include)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )

  add_definitions(-DQT_NO_KEYWORDS)
  
  #--------------------------------
  # Demo: Largest_empty_rectangle_2
  #--------------------------------
  # UI files (Qt Designer files)
  qt5_wrap_ui( DT_UI_FILES Largest_empty_rectangle_2.ui )

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ./Largest_empty_rectangle_2.qrc )

  # use the Qt MOC preprocessor on classes that derives from QObject
  qt5_generate_moc( Largest_empty_rectangle_2.cpp Largest_empty_rectangle_2.moc )

  # The executable itself.
  add_executable  ( Largest_empty_rectangle_2 Largest_empty_rectangle_2.cpp Largest_empty_rectangle_2.moc ${DT_UI_FILES}  ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES} )

  qt5_use_modules(Largest_empty_rectangle_2 Xml Script OpenGL Svg)

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Largest_empty_rectangle_2 )

  # Link with Qt libraries
  target_link_libraries( Largest_empty_rectangle_2 ${QT_LIBRARIES} )
  # And with CGAL libraries
  target_link_libraries( Largest_empty_rectangle_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

else()

  message(STATUS "NOTICE: This demo requires CGAL, and Qt5, and will not be compiled.")

endif()
