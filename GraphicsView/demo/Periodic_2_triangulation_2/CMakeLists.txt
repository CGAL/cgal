project (Periodic_2_triangulation_2_Demo)

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

find_package(CGAL COMPONENTS Qt5)

include(${CGAL_USE_FILE})



find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)
include_directories (BEFORE ../../../include)
include_directories (BEFORE ./include)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )

add_definitions(-DQT_NO_KEYWORDS)


#--------------------------------
# The "2D Periodic triangulation" demo: Periodic_2_triangulation_2
#--------------------------------
# UI files (Qt Designer files)
qt5_wrap_ui( DT_UI_FILES Periodic_2_triangulation_2.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ./Periodic_2_triangulation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt5_generate_moc( Periodic_2_Delaunay_triangulation_2.cpp Periodic_2_triangulation_2.moc )

# find header files for projects that can show them

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()
file(GLOB headers "*.h")
file(GLOB QT_headers "include/CGAL/Qt/*.h")
file(GLOB P2T2_headers "../../../include/CGAL/*.h")
SOURCE_GROUP("P2T2" FILES ${P2T2_headers})
SOURCE_GROUP("QT" FILES ${QT_headers})

# The executable itself.
add_executable  ( Periodic_2_Delaunay_triangulation_2 
  Periodic_2_Delaunay_triangulation_2.cpp
  Periodic_2_triangulation_2.moc 
  ${DT_UI_FILES}  ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES}
  ${headers} ${QT_headers} ${P2T2_headers})

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Periodic_2_Delaunay_triangulation_2 )

target_link_libraries( Periodic_2_Delaunay_triangulation_2 PRIVATE
  CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Gui )

include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
cgal_add_compilation_test(Periodic_2_Delaunay_triangulation_2)
else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
