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
# The demo: Circular_kernel_2
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( DT_UI_FILES Circular_kernel_2.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( DT_RESOURCE_FILES ./Circular_kernel_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( Circular_kernel_2.cpp Circular_kernel_2.moc )

# The executable itself.
add_executable  ( Circular_kernel_2 Circular_kernel_2.cpp ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

# Link with Qt libraries
target_link_libraries( Circular_kernel_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Circular_kernel_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})