# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (GraphicsView)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )


find_package(Qt4       REQUIRED)
include(${QT_USE_FILE})
include_directories (BEFORE ../../include)
link_directories(../../src/CGALQt4)

#--------------------------------
# The "constrained Delaunay" demo: constrained_Delaunay_triangulation
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( UI_FILES_CDT Constrained_Delaunay_triangulation_MainWindow.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( RESOURCE_FILES_CDT ./Constrained_Delaunay_triangulation_MainWindow.qrc ../icons/Input.qrc ../icons/File.qrc ../icons/Triangulation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_automoc( Constrained_Delaunay_triangulation_MainWindow.cpp )
qt4_automoc( DemosMainWindow.cpp )

# The executable itself.
add_executable  ( constrained_Delaunay_triangulation constrained_Delaunay_triangulation.cpp Constrained_Delaunay_triangulation_MainWindow.cpp DemosMainWindow.cpp ${UI_FILES_CDT} ${RESOURCE_FILES_CDT} )

# Link with Qt libraries
target_link_libraries( constrained_Delaunay_triangulation ${QT_LIBRARIES} )
if ( NOT AUTO_LINK_ENABLED )    
  # If auto-link is not enabled (on Linux for example), link with libCGAL explicitly.
  target_link_libraries( constrained_Delaunay_triangulation CGAL )
endif()
target_link_libraries( constrained_Delaunay_triangulation CGAL-Qt4 )

#--------------------------------
# The "Delaunay" demo: Delaunay_triangulation
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( UI_FILES_DT Delaunay_triangulation_MainWindow.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( RESOURCE_FILES_DT ./Delaunay_triangulation_MainWindow.qrc ../icons/Input.qrc ../icons/File.qrc ../icons/Triangulation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_automoc( Delaunay_triangulation_MainWindow.cpp )

# The executable itself.
add_executable  ( Delaunay_triangulation Delaunay_triangulation.cpp Delaunay_triangulation_MainWindow.cpp DemosMainWindow.cpp ${UI_FILES_DT} ${RESOURCE_FILES_DT} )

# Link with Qt libraries
target_link_libraries( Delaunay_triangulation ${QT_LIBRARIES} )
if ( NOT AUTO_LINK_ENABLED )    
  # If auto-link is not enabled (on Linux for example), link with libCGAL explicitly.
  target_link_libraries( Delaunay_triangulation CGAL )
endif()
target_link_libraries( Delaunay_triangulation CGAL-Qt4 )


#----------------------------------------
# minimal graphics view executable: min 
#----------------------------------------
add_executable  ( min EXCLUDE_FROM_ALL min.cpp )
qt4_generate_moc( min.cpp min.moc )
target_link_libraries( min ${QT_LIBRARIES} )

