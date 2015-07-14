# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Polyline_simplification_2)

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL COMPONENTS Qt5)

include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Widgets Svg)

include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../../include )
include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )

add_definitions(-DQT_NO_KEYWORDS)


# UI files (Qt Designer files)
qt5_wrap_ui( CDT_UI_FILES Polyline_simplification_2.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( CDT_RESOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/Polyline_simplification_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt5_generate_moc( ${CMAKE_CURRENT_SOURCE_DIR}/Polyline_simplification_2.cpp Polyline_simplification_2.moc )

# The executable itself.
add_executable  ( Polyline_simplification_2 ${CMAKE_CURRENT_SOURCE_DIR}/Polyline_simplification_2.cpp Polyline_simplification_2.moc ${CDT_UI_FILES} ${CDT_RESOURCE_FILES} )

qt5_use_modules(Polyline_simplification_2 Xml Script OpenGL Widgets Svg)
# Link with Qt libraries
target_link_libraries( Polyline_simplification_2 ${QT_LIBRARIES} )

# Link with CGAL
target_link_libraries( Polyline_simplification_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polyline_simplification_2 )

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
