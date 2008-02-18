set ( prj Surface_mesher_Qt4_Demo )

project ( ${prj} )

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )

# QGLViwer needs Qt4 configured with QtOpenGL and QtXml support

set( QT_USE_QTOPENGL TRUE )
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )

find_package(CGAL      REQUIRED )
find_package(Qt4       REQUIRED )
find_package(QGLViewer REQUIRED )
find_package(OpenGL    REQUIRED )

if ( CGAL_FOUND AND QT_FOUND AND QGLVIEWER_FOUND )

  include( ${QT_USE_FILE})
  include(UseCGAL)
  
  # Make sure the compiler can find local include files
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include")
    include_directories (BEFORE include)
  endif()

  foreach(INC_DIR ../../../include ../../../../Data_structure_for_queries_3/include ../../../../Marching_cube/include ../../../../CGALimageIO/include)
    if (EXISTS ${INC_DIR})
      include_directories (BEFORE ${INC_DIR})
    endif()
  endforeach()
  
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
  
  include_directories( ${QGLVIEWER_INCLUDE_DIR} )

  set( sources colorlisteditor.cpp isovalues_list.cpp mainwindow.cpp surface_mesher.cpp viewer.cpp volume.cpp ../../../../Marching_cube/src/mc/ply.c)
  
  qt4_automoc( ${sources} )

  qt4_generate_moc(surface.h moc_surface.cpp)

  set( sources ${sources} moc_surface.cpp)
  
  qt4_wrap_ui( uis ui/isovalues_list.ui ui/mainwindow.ui ui/meshing_bar.ui ui/optionsdialog.ui )
    
  qt4_add_resources( qrc_sources isovalues_list.qrc surface_mesher.qrc  ) 

  add_executable  ( ${prj}  ${sources} ${uis} ${qrc_sources} )
  
#  add_dependencies( ${prj} CGAL CGAL_IMAGEIO)
    
  # Link the executable to CGAL and third-party libraries
  if ( NOT AUTO_LINK_ENABLED )    
    target_link_libraries( ${prj} CGAL CGAL-ImageIO ${CGAL_3RD_PARTY_LIBRARIES})
  endif()
  
  target_link_libraries( ${prj} ${QT_LIBRARIES} ${QGLVIEWER_LIBRARY} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )
  
endif()
