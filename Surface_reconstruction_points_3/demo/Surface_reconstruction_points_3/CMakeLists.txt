# This is the CMake script for compiling the CGAL Point Set demo.

project( Point_set_demo )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

#option(POINT_SET_DEMO_ENABLE_FORWARD_DECL "In the Point Set demo, enable " OFF)
#mark_as_advanced(POINT_SET_DEMO_ENABLE_FORWARD_DECL)

# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

# Find CGAL and CGAL Qt4
find_package(CGAL COMPONENTS Qt4)
include( ${CGAL_USE_FILE} )

# Find Qt4 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(QT4_FOUND)
  include(${QT_USE_FILE})
  find_package(QGLViewer )
endif(QT4_FOUND)

# Find BLAS, LAPACK and TAUCS (optional), for Poisson
find_package(TAUCS)

# Find BLAS and LAPACK only (optional), for Jet Fitting
find_package(LAPACK)

# Find Eigen3 (requires 3.1.0 or greater)
find_package(Eigen3 3.1.0)

# Find Glew (optional), for splatting
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/GlSplat/cmake)
find_package(GLEW)

if(CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  # Link with BLAS, LAPACK and TAUCS (optional), for Poisson
  if(TAUCS_FOUND)
    include( ${TAUCS_USE_FILE} )
  else(TAUCS_FOUND)
    message(STATUS "NOTICE: TAUCS library is not found.")
  endif(TAUCS_FOUND)
  
  if(EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
  else()
    message(STATUS "NOTICE: Eigen library is not found.")
  endif()

  # Link with BLAS and LAPACK only (optional), for Jet Fitting
  if(LAPACK_FOUND)
    include( ${LAPACK_USE_FILE} )
  else(LAPACK_FOUND)
    if (NOT EIGEN3_FOUND)
      message(STATUS "NOTICE: Nor Eigen 3.1 (or greater) nor LAPACK library were found. Normal estimation and smoothing will not be available.")
    endif()
  endif(LAPACK_FOUND)

  if(GLEW_FOUND)
    include_directories ( ${GLEW_INCLUDE_DIR} )
    add_definitions(-DCGAL_GLEW_ENABLED)
  else(GLEW_FOUND)
    message(STATUS "NOTICE: GLEW library is not found. Splat rendering will not be available.")
  endif(GLEW_FOUND)

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    # Use /EHa option to catch stack overflows.
    # Note: TAUCS needs a stack >= 2MB. CGAL default is 10MB.
    string(REGEX REPLACE "/EH[asc]*" "/EHa" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

    # Use /FR to turn on IntelliSense
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /FR")

    # Allow Windows applications to use up to 3GB of RAM
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

    # Turn off stupid VC++ warnings
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4311 /wd4800 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")

    # Prints new compilation options
    message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
    message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
    message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
    message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
  endif()

  # Temporary debugging stuff
  ADD_DEFINITIONS( "-DDEBUG_TRACE" ) # turn on traces

  qt4_wrap_ui( UI_FILES MainWindow.ui)

  include(AddFileDependencies)

  qt4_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp" )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt4_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt4_generate_moc( "Scene.h" "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp" )
  add_file_dependencies( Scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" )

  qt4_add_resources ( RESOURCE_FILES Point_set_demo.qrc )

  qt4_automoc(Scene_item.cpp
              Scene_plane_item.cpp
              Point_set_scene_item.cpp
              Scene_polyhedron_item.cpp)

  # AUXILIARY LIBRARIES

  # put plugins (which are shared libraries) at the same location as
  # executable files

  set(LIBRARY_OUTPUT_PATH ${RUNTIME_OUTPUT_PATH})

  add_library(PS_demo_scene_item SHARED
    Scene_item.cpp Scene_item.moc
    Scene_item_with_display_list.cpp
    Polyhedron_demo_plugin_helper.cpp)
  target_link_libraries(PS_demo_scene_item ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${GLEW_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})

  add_library(PS_demo_scene_basic_objects SHARED
    Scene_plane_item.cpp Scene_plane_item.moc)
  target_link_libraries(PS_demo_scene_basic_objects PS_demo_scene_item ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})

  add_library(PS_demo_scene_polyhedron_item SHARED
    Scene_polyhedron_item.cpp Scene_polyhedron_item.moc)
  target_link_libraries(PS_demo_scene_polyhedron_item PS_demo_scene_item)

  if(GLEW_FOUND)
    qt4_add_resources(gl_splat_rc GlSplat/glsplat.qrc)
    add_library(gl_splat SHARED
      GlSplat/GlSplat.cpp GlSplat/Shader.cpp ${gl_splat_rc})
    target_link_libraries(gl_splat ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${GLEW_LIBRARIES})
  endif(GLEW_FOUND)

  add_library(point_set SHARED
    Point_set_scene_item.cpp Point_set_scene_item.moc)
  target_link_libraries(point_set PS_demo_scene_item)

  foreach(lib PS_demo_scene_item PS_demo_scene_basic_objects PS_demo_scene_polyhedron_item point_set)
    add_to_cached_list(CGAL_EXECUTABLE_TARGETS ${lib})
  endforeach()

  if(GLEW_FOUND)
    target_link_libraries( point_set gl_splat ${GLEW_LIBRARIES} )
  endif(GLEW_FOUND)

  add_definitions(-DQT_STATICPLUGIN)

#  if(POINT_SET_DEMO_ENABLE_FORWARD_DECL)
    add_definitions(-DUSE_FORWARD_DECL)
    add_executable  ( Point_set_demo
      ${UI_FILES}
      MainWindow.cpp
      Point_set_demo.cpp
      Viewer.cpp
      Scene.cpp
      MainWindow_moc.cpp
      Scene_moc.cpp
      Viewer_moc.cpp
      ${RESOURCE_FILES} )
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Point_set_demo )
#  else(POINT_SET_DEMO_ENABLE_FORWARD_DECL)
#    add_file_dependencies( Point_set_demo.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
#                                            "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp"
#					    "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
#    add_executable  ( Point_set_demo Point_set_demo.cpp ${UI_FILES} ${RESOURCE_FILES} )
#  endif(POINT_SET_DEMO_ENABLE_FORWARD_DECL)


  # Link with Qt libraries
  target_link_libraries( Point_set_demo ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( Point_set_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Point_set_demo ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

  # Link with the scene_item libraries
  target_link_libraries( Point_set_demo PS_demo_scene_item PS_demo_scene_polyhedron_item point_set )

  if(GLEW_FOUND)
    target_link_libraries( Point_set_demo gl_splat ${GLEW_LIBRARIES} )
  endif(GLEW_FOUND)

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Point_set_demo )


  ###########
  # PLUGINS #
  ###########
  remove_definitions(-DQT_STATICPLUGIN)

  macro(point_set_demo_plugin plugin_name plugin_implementation_base_name)
    list_split(option ARGN_TAIL ${ARGN} )
    if(NOT ${option} STREQUAL "EXCLUDE_FROM_ALL")
      set(other_sources ${ARGN})
      set(option "")
    else()
      set(other_sources ${ARGN_TAIL})
    endif()
    qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" ${plugin_implementation_base_name}.moc )
    add_file_dependencies( ${plugin_implementation_base_name}.moc "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" )

    add_library(${plugin_name} MODULE ${option} ${plugin_implementation_base_name}.moc ${plugin_implementation_base_name}.cpp ${other_sources})
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} )
    # Link with scene_item
    target_link_libraries( ${plugin_name} PS_demo_scene_item)
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  endmacro(point_set_demo_plugin)

  if(TAUCS_FOUND OR EIGEN3_FOUND)
    qt4_wrap_ui( POISSON_UI_FILES PS_demo_poisson_plugin.ui)
    point_set_demo_plugin(PS_demo_poisson_plugin
      PS_demo_poisson_plugin
      PS_demo_poisson_plugin_cgal_code.cpp
      ${POISSON_UI_FILES})
    target_link_libraries(PS_demo_poisson_plugin PS_demo_scene_polyhedron_item point_set)
  else()
    message(STATUS "NOTICE: Neither TAUCS nor Eigen 3.1 (or greater) libraries have been found. Poisson reconstruction will not be available.")
  endif()

  point_set_demo_plugin(PS_demo_inside_out_plugin PS_demo_inside_out_plugin)
  target_link_libraries(PS_demo_inside_out_plugin PS_demo_scene_polyhedron_item point_set)

  point_set_demo_plugin(PS_demo_off_plugin PS_demo_off_plugin)
  target_link_libraries(PS_demo_off_plugin PS_demo_scene_polyhedron_item point_set)

  point_set_demo_plugin(PS_demo_xyz_plugin PS_demo_xyz_plugin)
  target_link_libraries(PS_demo_xyz_plugin PS_demo_scene_polyhedron_item point_set)

  qt4_wrap_ui(SIMPLIFICATION_UI_FILES PS_demo_simplification_plugin.ui)
  point_set_demo_plugin(PS_demo_simplification_plugin
    PS_demo_simplification_plugin
    ${SIMPLIFICATION_UI_FILES})
  target_link_libraries(PS_demo_simplification_plugin point_set)

  point_set_demo_plugin(PS_demo_local_spacing_plugin PS_demo_local_spacing_plugin)
  target_link_libraries(PS_demo_local_spacing_plugin point_set)

  point_set_demo_plugin(PS_demo_average_spacing_plugin PS_demo_average_spacing_plugin)
  target_link_libraries(PS_demo_average_spacing_plugin point_set)

  if(LAPACK_FOUND OR EIGEN3_FOUND)
    point_set_demo_plugin(PS_demo_smoothing_plugin PS_demo_smoothing_plugin)
    target_link_libraries(PS_demo_smoothing_plugin point_set)
  endif()

  if(LAPACK_FOUND)
    qt4_wrap_ui( NORMAL_UI_FILES PS_demo_normal_estimation_plugin.ui)
    point_set_demo_plugin(PS_demo_normal_estimation_plugin
      PS_demo_normal_estimation_plugin
      ${NORMAL_UI_FILES})
    target_link_libraries(PS_demo_normal_estimation_plugin point_set)
  endif()

  qt4_wrap_ui( CLEANING_UI_FILES PS_demo_cleaning_plugin.ui)
  point_set_demo_plugin(PS_demo_cleaning_plugin
    PS_demo_cleaning_plugin
    ${CLEANING_UI_FILES})
  target_link_libraries(PS_demo_cleaning_plugin point_set)

else (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(POINT_SET_DEMO_MISSING_DEPS "")

  if(NOT CGAL_Qt4_FOUND)
    set(POINT_SET_DEMO_MISSING_DEPS "the CGAL Qt4 library, ${POINT_SET_DEMO_MISSING_DEPS}")
  endif()

  if(NOT QT4_FOUND)
    set(POINT_SET_DEMO_MISSING_DEPS "Qt4, ${POINT_SET_DEMO_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(POINT_SET_DEMO_MISSING_DEPS "OpenGL, ${POINT_SET_DEMO_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(POINT_SET_DEMO_MISSING_DEPS "QGLViewer, ${POINT_SET_DEMO_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${POINT_SET_DEMO_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
