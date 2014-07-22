#set ( prj Surface_mesher )

project ( Surface_mesher )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

set(PACKAGE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)

# Add several CGAL packages to the include and link paths,
# if they lie in ${PACKAGE_ROOT}/.
foreach(INC_DIR ${PACKAGE_ROOT}/include ${PACKAGE_ROOT}/../Mesh_2/include ${PACKAGE_ROOT}/../Data_structure_for_queries_3/include ${PACKAGE_ROOT}/../Marching_cube/include ${PACKAGE_ROOT}/../CGAL_ImageIO/include ${PACKAGE_ROOT}/../GraphicsView/include)
  if (EXISTS ${INC_DIR})
    include_directories (BEFORE ${INC_DIR})
  endif()
endforeach()
foreach(LIB_DIR ${PACKAGE_ROOT}/../CGAL_ImageIO/src/CGAL_ImageIO ${PACKAGE_ROOT}/../GraphicsView/src/CGALQt4)
  if (EXISTS ${LIB_DIR})
    link_directories (${LIB_DIR})
  endif()
endforeach()

include_directories( ./ )

#New for Qt5 version !

# QGLViwer needs Qt4 configured with QtOpenGL and QtXml support

set( QT_USE_QTOPENGL TRUE )
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )

option (USE_QT5 "Use Qt5 version instead of Qt4" ON)

if(USE_QT5)

	find_package(CGAL COMPONENTS ImageIO Qt5)
	include(${CGAL_USE_FILE})
	find_package(Qt5)

else (USE_QT5)
	find_package(CGAL COMPONENTS ImageIO Qt4)

	include(${CGAL_USE_FILE})

	find_package(Qt4)
	include(${QT_USE_FILE})

	if(CGAL_Qt4_FOUND AND QT4_FOUND)
		message("Qt4 found")
		set(QT4 TRUE)
	endif()
endif (USE_QT5)

include(QtChoice)

if ( CGAL_FOUND AND (QT4 OR QT5) AND CGAL_ImageIO_FOUND)
  include( ${CGAL_USE_FILE} )

  find_package(QGLViewer )
  find_package(OpenGL    )
  
  if ( QGLVIEWER_FOUND )
        
    include_directories( ${QGLVIEWER_INCLUDE_DIR} )
    add_definitions(${QGLVIEWER_DEFINITIONS})

    set( sources Raw_image_dialog.cpp colorlisteditor.cpp values_list.cpp mainwindow.cpp Surface_mesher.cpp viewer.cpp volume.cpp )

    if(EXISTS ${PACKAGE_ROOT}/../Marching_cube)
      option(SURFACE_MESH_DEMO_USE_MARCHING_CUBE "Embed a marching cube implementation in the Surface Mesh demo." ON)
      mark_as_advanced(SURFACE_MESH_DEMO_USE_MARCHING_CUBE)
    endif()
    
    if(EXISTS ${PACKAGE_ROOT}/include/CGAL/Polyhedral_surface_3.h)
      option(SURFACE_MESH_DEMO_WITH_POLYHEDRAL_SURFACE "Compile the support for polyhedral surfaces." OFF)
      mark_as_advanced(SURFACE_MESH_DEMO_WITH_POLYHEDRAL_SURFACE)
    endif()

    option(SURFACE_MESH_DEMO_VERBOSE "Set this option if you want the Surface Mesh demo to display messages on standard output." OFF)
    mark_as_advanced(SURFACE_MESH_DEMO_VERBOSE)

    if(SURFACE_MESH_DEMO_VERBOSE)
      add_definitions(-DCGAL_SURFACE_MESHER_VERBOSE)
    endif()

    if(SURFACE_MESH_DEMO_WITH_POLYHEDRAL_SURFACE)
      set(sources ${sources} polyhedral_surface.cpp)
    else(SURFACE_MESH_DEMO_WITH_POLYHEDRAL_SURFACE)
      add_definitions(-DCGAL_DO_NOT_USE_POLYHEDRAL_SURFACE)
    endif(SURFACE_MESH_DEMO_WITH_POLYHEDRAL_SURFACE)

    if(SURFACE_MESH_DEMO_USE_MARCHING_CUBE)
      set(sources ${sources} ${PACKAGE_ROOT}/../Marching_cube/src/mc/ply.c)
      add_definitions(-DCGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE)
    endif()
	
    qt_automoc( ${sources} )

    qt_generate_moc( "surface.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_surface.cpp" )

    set( sources ${sources} moc_surface.cpp)
    
    qt_wrap_ui( uis ui/values_list.ui ui/mainwindow.ui ui/optionsdialog.ui ui/raw_image.ui )
    
    qt_add_resources( qrc_sources values_list.qrc surface_mesher.qrc  ) 

    add_executable  ( Surface_mesher  
		      ${sources} 
                      ${uis} 
                      ${qrc_sources} )
    
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Surface_mesher )

    # Link the executable to CGAL and third-party libraries
    target_link_libraries( Surface_mesher 
			   ${CGAL_LIBRARIES} 				   
			   ${CGAL_3RD_PARTY_LIBRARIES})
    
    target_link_libraries( Surface_mesher 
		           ${QT_LIBRARIES} 
			   ${QGLVIEWER_LIBRARIES})

    target_link_libraries( Surface_mesher 
			   ${OPENGL_gl_LIBRARY} 
		   	   ${OPENGL_glu_LIBRARY} )
    
  else( QGLVIEWER_FOUND)
    message(STATUS "NOTICE: This demo needs libQGLViewer, and will not be compiled.")
  endif( QGLVIEWER_FOUND)
else(CGAL_FOUND AND (QT4 OR QT5) AND CGAL_ImageIO_FOUND)
  if(RUNNING_CGAL_AUTO_TEST)
    # Just to avoid a warning from CMake if that variable is set on the command line...
  endif()
  message(STATUS "NOTICE: This demo needs Qt (version 4 or 5), and will not be compiled.")
endif(CGAL_FOUND AND (QT4 OR QT5) AND CGAL_ImageIO_FOUND)
