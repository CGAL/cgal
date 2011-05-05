set ( prj Surface_mesher )

project ( ${prj} )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set(PACKAGE_ROOT ../..)

# Add several CGAL packages to the include and link paths,
# if they lie in ${PACKAGE_ROOT}/.
foreach(INC_DIR ${PACKAGE_ROOT}/include ${PACKAGE_ROOT}/../Mesh_2/include ${PACKAGE_ROOT}/../Data_structure_for_queries_3/include ${PACKAGE_ROOT}/../Marching_cube/include ${PACKAGE_ROOT}/../CGALimageIO/include ${PACKAGE_ROOT}/../GraphicsView/include)
  if (EXISTS ${INC_DIR})
    include_directories (BEFORE ${INC_DIR})
  endif()
endforeach()
foreach(LIB_DIR ${PACKAGE_ROOT}/../CGALimageIO/src/CGALimageIO ${PACKAGE_ROOT}/../GraphicsView/src/CGALQt4)
  if (EXISTS ${LIB_DIR})
    link_directories (${LIB_DIR})
  endif()
endforeach()

include_directories( ./ )

# QGLViwer needs Qt4 configured with QtOpenGL and QtXml support

set( QT_USE_QTOPENGL TRUE )
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )

find_package(CGAL COMPONENTS ImageIO Qt4)
if ( CGAL_FOUND AND CGAL_Qt4_FOUND)
  include( ${CGAL_USE_FILE} )

  find_package(Qt4 )
  find_package(QGLViewer )
  find_package(OpenGL    )
  
  if ( QGLVIEWER_FOUND )

    include( ${QT_USE_FILE})
        
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
    
    qt4_automoc( ${sources} )

    qt4_generate_moc( "surface.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_surface.cpp" )

    set( sources ${sources} moc_surface.cpp)
    
    qt4_wrap_ui( uis ui/values_list.ui ui/mainwindow.ui ui/optionsdialog.ui ui/raw_image.ui )
    
    qt4_add_resources( qrc_sources values_list.qrc surface_mesher.qrc  ) 

    add_executable  ( ${prj}  ${sources} ${uis} ${qrc_sources} )
    
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${prj} )

    # Link the executable to CGAL and third-party libraries
    target_link_libraries( ${prj} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
    
    target_link_libraries( ${prj} ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )
    
  else( QGLVIEWER_FOUND)
    message(STATUS "NOTICE: This demo needs libQGLViewer, and will not be compiled.")
  endif( QGLVIEWER_FOUND)
else(CGAL_FOUND AND CGAL_Qt4_FOUND)
  message(STATUS "NOTICE: This demo needs Qt4, and will not be compiled.")
endif(CGAL_FOUND AND CGAL_Qt4_FOUND)
