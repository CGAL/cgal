# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Hyperbolic_Triangulation_2_demo)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL COMPONENTS Qt4)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

#find_package(Qt4)

include_directories (BEFORE ../../../GraphicsView/demo/Triangulation_2)
include_directories (BEFORE ../../include)
# demos contain their own headers
include_directories (BEFORE include)

if ( CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND )

include(${QT_USE_FILE})

#--------------------------------
# The "Delaunay" demo: Delaunay_triangulation_2
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( DT_UI_FILES ../../../GraphicsView/demo/Triangulation_2/Delaunay_triangulation_2.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( DT_RESOURCE_FILES ../../../GraphicsView/demo/Triangulation_2/Delaunay_triangulation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( "Hyperbolic_Delaunay_triangulation_2_demo.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Hyperbolic_Delaunay_triangulation_2_demo.moc" )

# The executable itself.
add_executable  ( Hyperbolic_Delaunay_triangulation_2_demo Hyperbolic_Delaunay_triangulation_2_demo.cpp Hyperbolic_Delaunay_triangulation_2_demo.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Hyperbolic_Delaunay_triangulation_2_demo )

# Link with Qt libraries
target_link_libraries( Hyperbolic_Delaunay_triangulation_2_demo ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Hyperbolic_Delaunay_triangulation_2_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

#--------------------------------
# The "hyperbolic" demo: Hyperbolic_translations_2
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( DT_UI_FILES ui/Hyperbolic_translations_2.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( DT_RESOURCE_FILES resources/Hyperbolic_translations_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( "Hyperbolic_translations_2_demo.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Hyperbolic_translations_2_demo.moc" )

# The executable itself.
add_executable  ( Hyperbolic_translations_2_demo Hyperbolic_translations_2_demo.cpp Hyperbolic_translations_2_demo.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Hyperbolic_translations_2_demo )

# Link with Qt libraries
target_link_libraries( Hyperbolic_translations_2_demo ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Hyperbolic_translations_2_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

#--------------------------------
# The "Periodic Delaunay hyperbolic " demo:
#--------------------------------

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( "Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.moc" )

# The executable itself.
add_executable  ( Periodic_2_Delaunay_hyperbolic_triangulation_2_demo Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.cpp  Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Periodic_2_Delaunay_hyperbolic_triangulation_2_demo )

# Link with Qt libraries
target_link_libraries( Periodic_2_Delaunay_hyperbolic_triangulation_2_demo ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Periodic_2_Delaunay_hyperbolic_triangulation_2_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt4, and will not be compiled.")

endif()
