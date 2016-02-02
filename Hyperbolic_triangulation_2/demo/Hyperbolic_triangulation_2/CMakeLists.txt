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

project( Hyperbolic_triangulation_2 )


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Qt5Core Qt5Giu Qt5Widgets )

# Find the QtWidgets library
find_package(Qt5Core)
find_package(Qt5Gui)
find_package(Qt5Widgets)

# Add the include directories for the Qt 5 Widgets module to
# the compile lines.
include_directories(${Qt5Core_INCLUDE_DIRS})
include_directories(${Qt5Gui_INCLUDE_DIRS})
include_directories(${Qt5Widgets_INCLUDE_DIRS})

# Use the compile definitions defined in the Qt 5 Widgets module
add_definitions(${Qt5Core_DEFINITIONS})
add_definitions(${Qt5Gui_DEFINITIONS})
add_definitions(${Qt5Widgets_DEFINITIONS})

# Add compiler flags for building executables (-fPIE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Core_EXECUTABLE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Gui_EXECUTABLE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")

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

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

#--------------------------------
# The "Periodic Delaunay hyperbolic " demo:
#--------------------------------

qt5_generate_moc( Hyperbolic_Delaunay_triangulation_2_demo.cpp Hyperbolic_Delaunay_triangulation_2_demo.moc )
create_single_source_cgal_program( "Hyperbolic_Delaunay_triangulation_2_demo.cpp" )
target_link_libraries(Hyperbolic_Delaunay_triangulation_2_demo CGAL CGAL_Qt5 Qt5::Core Qt5::Gui Qt5::Widgets)

################


qt5_generate_moc( Hyperbolic_Dirichlet_region_2_demo.cpp Hyperbolic_Dirichlet_region_2_demo.moc )
create_single_source_cgal_program( "Hyperbolic_Dirichlet_region_2_demo.cpp" )
target_link_libraries(Hyperbolic_Dirichlet_region_2_demo CGAL CGAL_Qt5 Qt5::Core Qt5::Gui Qt5::Widgets)

################


qt5_generate_moc( Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.cpp Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.moc )
create_single_source_cgal_program( "Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.cpp" )
target_link_libraries(Periodic_2_Delaunay_hyperbolic_triangulation_2_demo CGAL CGAL_Qt5 Qt5::Core Qt5::Gui Qt5::Widgets)

################

qt5_generate_moc( Hyperbolic_translations_2_demo.cpp Hyperbolic_translations_2_demo.moc )
create_single_source_cgal_program( "Hyperbolic_translations_2_demo.cpp" )
target_link_libraries(Hyperbolic_translations_2_demo CGAL CGAL_Qt5 Qt5::Core Qt5::Gui Qt5::Widgets)




# Link with Qt libraries
target_link_libraries( Hyperbolic_Dirichlet_region_2_demo ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Hyperbolic_Dirichlet_region_2_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt4, and will not be compiled.")

endif()
