# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (GraphicsView)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED)

include(${CGAL_USE_FILE})
include(CreateSingleSourceCGALQt4Program)

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )


find_package(Qt4       REQUIRED)
include(${QT_USE_FILE})
include_directories (BEFORE ../../include)

#--------------------------------
# The "delaunay" demo: delaunay 
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( UI_FILES MainWindow.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( RESOURCE_FILES ./MainWindow.qrc ../icons/Input.qrc ../icons/File.qrc ../icons/Logos.qrc ../icons/Triangulation_2.qrc)

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( ../../include/CGAL/IO/QtNavigation.h QtNavigation.moc )
qt4_generate_moc( ../../include/CGAL/IO/QtGraphicsItem.h QtGraphicsItem.moc )
qt4_generate_moc( ../../include/CGAL/IO/QtInput.h QtInput.moc )
macro_add_file_dependencies( QNavigation.cpp QtNavigation.moc )
macro_add_file_dependencies( QInput.cpp QtInput.moc )
macro_add_file_dependencies( QGraphicsItem_2.cpp QtGraphicsItem.moc )
qt4_automoc( MainWindow.cpp )

# The executable itself.
add_executable  ( delaunay main.cpp MainWindow.cpp QNavigation.cpp QGraphicsItem_2.cpp QInput.cpp QPolylineInput_2_non_templated_base.cpp ${UI_FILES} ${RESOURCE_FILES} )

# Link with Qt libraries
target_link_libraries( delaunay ${QT_LIBRARIES} )
if ( NOT AUTO_LINK_ENABLED )    
  # If auto-link is not enabled (on Linux for example), link with libCGAL explicitly.
  target_link_libraries( delaunay CGAL )
endif()

#----------------------------------------
# minimal graphics view executable: min 
#----------------------------------------
add_executable  ( min min.cpp )
qt4_generate_moc( min.cpp min.moc )
target_link_libraries( min ${QT_LIBRARIES} )

