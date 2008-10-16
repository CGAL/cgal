cmake_minimum_required(VERSION 2.4.5)

set ( prj Surface_mesher )

project ( ${prj} )

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(PACKAGE_ROOT ../../..)


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

option(WITH_POLYHEDRAL_SURFACE "Compile the support for polyhedral surfaces." ON)

# QGLViwer needs Qt4 configured with QtOpenGL and QtXml support

set( QT_USE_QTOPENGL TRUE )
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )

find_package(CGAL      REQUIRED COMPONENTS ImageIO Qt4)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  find_package(Qt4       REQUIRED )
  find_package(QGLViewer REQUIRED )
  find_package(OpenGL    REQUIRED )
  find_package(ZLIB      REQUIRED )
  find_package(VTK       QUIET)

  find_path(bimap_DIR boost/bimap/bimap.hpp ${Boost_INCLUDE_DIRS} DOC "The directory containing the bimap subdirectory.")

  if(EXISTS ${bimap_DIR}/boost/bimap/bimap.hpp)
    if("${bimap_DIR}" STREQUAL "${Boost_INCLUDE_DIRS}")
    else()
      include_directories(${bimap_DIR})
      add_definitions(-DCGAL_USE_BOOST_BIMAP)
    endif()
  endif(EXISTS ${bimap_DIR}/boost/bimap/bimap.hpp)

  if( VTK_FOUND AND EXISTS ${VTK_USE_FILE})
#    if( ${VTK_MAJOR_VERSION} STRGREATER "5")
      message(STATUS "VTK-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}.${VTK_BUILD_VERSION} found in ${VTK_DIR}.")
      set(VTK_BUILD_SETTINGS_FILE "")
      include( ${VTK_USE_FILE} )
      add_definitions(-DCGAL_USE_VTK)
#    endif( ${VTK_MAJOR_VERSION} GREATER 5)
    SET(PRJ_VTK_LIBS vtkImaging vtkIO)
  else( VTK_FOUND AND EXISTS ${VTK_USE_FILE})
    message(STATUS "VTK not found.")
  endif( VTK_FOUND  AND EXISTS ${VTK_USE_FILE})

  add_definitions(-DCGAL_USE_ZLIB)

  if ( QT_FOUND AND QGLVIEWER_FOUND )

    include( ${QT_USE_FILE})
    
    # Make sure the compiler can find local include files
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include")
      include_directories (BEFORE include)
    endif()

    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
    include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
    
    include_directories( ${QGLVIEWER_INCLUDE_DIR} )

    set( sources colorlisteditor.cpp values_list.cpp mainwindow.cpp Surface_mesher.cpp viewer.cpp volume.cpp )

    if(WITH_POLYHEDRAL_SURFACE)
      set(sources ${sources} polyhedral_surface.cpp ${PACKAGE_ROOT}/../Marching_cube/src/mc/ply.c)
    else(WITH_POLYHEDRAL_SURFACE)
      add_definitions(-DCGAL_DO_NOT_USE_POLYHEDRAL_SURFACE)
    endif(WITH_POLYHEDRAL_SURFACE)
    
    qt4_automoc( ${sources} )

    qt4_generate_moc(surface.h moc_surface.cpp)

    set( sources ${sources} moc_surface.cpp)
    
    qt4_wrap_ui( uis ui/values_list.ui ui/mainwindow.ui ui/optionsdialog.ui ui/raw_image.ui )
    
    qt4_add_resources( qrc_sources values_list.qrc surface_mesher.qrc  ) 

    add_executable  ( ${prj}  ${sources} ${uis} ${qrc_sources} )
    
    # Link the executable to CGAL and third-party libraries
    target_link_libraries( ${prj} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
    
    target_link_libraries( ${prj} ${QT_LIBRARIES} ${QGLVIEWER_LIBRARY} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} ${PRJ_VTK_LIBS} )
    
  endif()
endif(CGAL_FOUND)
