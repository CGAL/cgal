# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Triangulation_2)

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL COMPONENTS Qt5)

include(${CGAL_USE_FILE})


find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)

include_directories (BEFORE ../../include)
include_directories (BEFORE ./include)
include_directories (BEFORE ../../../Number_types/include)

include_directories (BEFORE ../../../Triangulation_2/include)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )

add_definitions(-DQT_NO_KEYWORDS)



#--------------------------------
# The "constrained Delaunay" demo: Constrained_Delaunay_triangulation_2
#--------------------------------
# UI files (Qt Designer files)
qt5_wrap_ui( CDT_UI_FILES Constrained_Delaunay_triangulation_2.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( CDT_RESOURCE_FILES ./Constrained_Delaunay_triangulation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt5_generate_moc( "Constrained_Delaunay_triangulation_2.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Constrained_Delaunay_triangulation_2.moc" )

# The executable itself.
add_executable  ( Constrained_Delaunay_triangulation_2 Constrained_Delaunay_triangulation_2.cpp Constrained_Delaunay_triangulation_2.moc ${CDT_UI_FILES} ${CDT_RESOURCE_FILES} )

qt5_use_modules(Constrained_Delaunay_triangulation_2 Xml Script OpenGL Svg)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Constrained_Delaunay_triangulation_2 )

# Link with Qt libraries
target_link_libraries( Constrained_Delaunay_triangulation_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Constrained_Delaunay_triangulation_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

#--------------------------------
# The "Delaunay" demo: Delaunay_triangulation_2
#--------------------------------
# UI files (Qt Designer files)
qt5_wrap_ui( DT_UI_FILES Delaunay_triangulation_2.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( DT_RESOURCE_FILES ./Delaunay_triangulation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt5_generate_moc( "Delaunay_triangulation_2.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Delaunay_triangulation_2.moc" )

# The executable itself.
add_executable  ( Delaunay_triangulation_2 Delaunay_triangulation_2.cpp Delaunay_triangulation_2.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

qt5_use_modules(Delaunay_triangulation_2 Xml Script OpenGL Svg)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Delaunay_triangulation_2 )

# Link with Qt libraries
target_link_libraries( Delaunay_triangulation_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Delaunay_triangulation_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

#--------------------------------
# The "Regular" demo: Regular_triangulation_2
#--------------------------------
# UI files (Qt Designer files)
qt5_wrap_ui( DT_UI_FILES Regular_triangulation_2.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( DT_RESOURCE_FILES ./Regular_triangulation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt5_generate_moc( "Regular_triangulation_2.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Regular_triangulation_2.moc" )

# The executable itself.
add_executable  ( Regular_triangulation_2 Regular_triangulation_2.cpp Regular_triangulation_2.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

qt5_use_modules(Regular_triangulation_2 Xml Script OpenGL Svg)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Regular_triangulation_2 )

# Link with Qt libraries
target_link_libraries( Regular_triangulation_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Regular_triangulation_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
