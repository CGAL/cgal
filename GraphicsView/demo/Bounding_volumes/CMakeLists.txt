# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Bounding_volumes_Demo)

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

find_package(CGAL COMPONENTS Qt5)

include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)
include_directories (BEFORE ../../include)
include_directories (BEFORE ./include)
include_directories (BEFORE ../../../Matrix_search/include)


if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )

add_definitions(-DQT_NO_KEYWORDS)

#----------------------------------------------
# The "Bounding volumes" demo: Bounding_volumes
#----------------------------------------------
# UI files (Qt Designer files)
qt5_wrap_ui( DT_UI_FILES Bounding_volumes.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ./Bounding_volumes.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt5_generate_moc( "Bounding_volumes.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Bounding_volumes.moc" )

# The executable itself.
add_executable  ( Bounding_volumes Bounding_volumes.cpp Bounding_volumes.moc ${DT_UI_FILES} ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})

qt5_use_modules(Bounding_volumes Xml Script OpenGL Svg)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Bounding_volumes )

# Link with Qt libraries
target_link_libraries( Bounding_volumes ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Bounding_volumes ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
