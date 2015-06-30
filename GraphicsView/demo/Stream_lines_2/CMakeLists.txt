# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Stream_lines_2_demo)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL REQUIRED COMPONENTS Qt5)

include(${CGAL_USE_FILE})


find_package(Qt5 COMPONENTS Xml Script OpenGL)
include_directories (BEFORE ../../include)


if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )

add_definitions(-DQT_NO_KEYWORDS)


#--------------------------------
# UI files (Qt Designer files)
qt5_wrap_ui( DT_UI_FILES Stream_lines_2.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( DT_RESOURCE_FILES ./Stream_lines_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt5_generate_moc( Stream_lines_2.cpp Stream_lines_2.moc )

# The executable itself.
add_executable  ( Stream_lines_2 Stream_lines_2.cpp Stream_lines_2.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

qt5_use_modules(Stream_lines_2 Xml Script OpenGL)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Stream_lines_2 )

# Link with Qt libraries
target_link_libraries( Stream_lines_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Stream_lines_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
