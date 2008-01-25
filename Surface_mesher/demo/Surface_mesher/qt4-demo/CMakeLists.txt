set ( prj Surface_mesher_Qt4_Demo )

project ( ${prj} )

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_BINARY_DIR}/cmake/modules" )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules"  )

# QGLViwer needs Qt4 configured with QtOpenGL and QtXml support

set( QT_USE_QTOPENGL TRUE )
set( QT_USE_QTXML    TRUE )
set(QT_USE_QTMAIN    TRUE )

find_package(CGAL      REQUIRED )
find_package(Qt4       REQUIRED )
find_package(QGLViewer REQUIRED )

if ( CGAL_FOUND AND QT_FOUND AND QGLVIEWER_FOUND )

  include( ${QT_USE_FILE})
  include(UseCGAL)
  
  # Make sure the compiler can find local include files
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include")
    include_directories (BEFORE include)
  endif()
  
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
  
  include_directories( ${QGLVIEWER_INCLUDE_DIR} )

  set( sources surface_mesher.cpp viewer.cpp mainwindow.cpp polyhedral_surface.cpp )
  
  qt4_automoc( ${sources} )
  
  qt4_wrap_ui( uis mainwindow.ui optionsdialog.ui )
    
  qt4_add_resources( surface_mesher.qrc ) 

  add_definitions ( -DCGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE )
  
  add_executable  ( ${prj}  ${sources} ${uis} )
  
  add_dependencies( ${prj} CGAL)
    
  # Link the executable to CGAL and third-party libraries
  if ( NOT AUTO_LINK_ENABLED )    
    target_link_libraries( ${prj} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  endif()
  
  target_link_libraries( ${prj} ${QT_LIBRARIES} ${QGLVIEWER_LIBRARY} )
  
endif()
