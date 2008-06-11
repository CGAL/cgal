# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (GraphicsView)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALQt4Program)

  set( QT_USE_QTXML    TRUE )
  set( QT_USE_QTMAIN   TRUE )
  set( QT_USE_QTSCRIPT  TRUE )


  find_package(QT4       REQUIRED)
  include(${QT_USE_FILE})
  include_directories (BEFORE ../../include)

  qt4_wrap_ui( UI_FILES MainWindow.ui )
  qt4_add_resources ( RESOURCE_FILES ./MainWindow.qrc ../icons/Input.qrc ../icons/File.qrc ../icons/Logos.qrc ../icons/Triangulation_2.qrc)
  qt4_automoc( MainWindow.cpp QNavigation.cpp QNavigation2.cpp QGraphicsItem_2.cpp QInput.cpp QDatastructure.cpp )

  add_executable  ( delaunay main.cpp MainWindow.cpp QNavigation.cpp QNavigation2.cpp QDatastructure.cpp QGraphicsItem_2.cpp QInput.cpp ${UI_FILES} ${RESOURCE_FILES} )
  target_link_libraries( delaunay ${QT_LIBRARIES} )
endif()
