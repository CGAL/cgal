cmake_minimum_required(VERSION 3.16...3.31)
project(Lab_Demo)
include(FeatureSummary)


add_compile_definitions(CGAL_NO_DEPRECATED_CODE)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

#Defines flags to emulate windows behavior for linking error generation
if(CMAKE_CXX_COMPILER_ID STREQUAL Clang
   OR CMAKE_COMPILER_IS_GNUCC
   OR CMAKE_COMPILER_IS_GNUCXX)
  if(UNIX OR APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
  endif()
  if(UNIX AND NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -z defs")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -z defs")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -z defs")
  endif()
endif()

# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Include this package's headers first
include_directories(BEFORE ./ ./include ./CGAL_demo)
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}")

# Find CGAL
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6 ImageIO)
set_package_properties(CGAL PROPERTIES TYPE REQUIRED)

# Find Qt6 itself
find_package(Qt6 QUIET
             COMPONENTS OpenGLWidgets Widgets Qml
             OPTIONAL_COMPONENTS WebSockets Network)

set_package_properties(
  Qt6 PROPERTIES
  TYPE REQUIRED
  PURPOSE "Enables the 3D Features, for GUI and visualization."
  DESCRIPTION "To find this package, it should be sufficient to fill the Qt6_DIR variable with: <Qt_dir>/<Qt_version>/<Compilator>/lib/cmake/Qt6")

if(Qt6_FOUND)
  add_compile_definitions(QT_NO_KEYWORDS)
endif(Qt6_FOUND)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

set_package_properties(
  Eigen3 PROPERTIES
  DESCRIPTION "A library for linear algebra."
  PURPOSE "Required for most plugins (Meshing, Mesh and Point Set Processing, etc.).")

if(NOT TARGET CGAL::Eigen3_support)
  message(STATUS "NOTICE: Eigen was not found.")
endif()

find_package(METIS QUIET)
include(CGAL_METIS_support)
set_package_properties(
  METIS PROPERTIES
  DESCRIPTION "A library for graph partitioning."
  PURPOSE "Required for the partition plugin.")

# Activate concurrency?
option(LAB_DEMO_ACTIVATE_CONCURRENCY "Enable concurrency" ON)
if(LAB_DEMO_ACTIVATE_CONCURRENCY)
  find_package(TBB QUIET)
  include(CGAL_TBB_support)
  if(NOT TARGET CGAL::TBB_support)
    message(STATUS "NOTICE: Intel TBB was not found. Bilateral smoothing and WLOP plugins are faster if TBB is linked.")
  endif()
endif()

#find libssh for scene sharing
find_package(LibSSH QUIET)
set_package_properties(
  LibSSH PROPERTIES
  DESCRIPTION "A library implementing the SSH protocol on client and server side. "
  PURPOSE "Required for loading (saving) a scene to (from) a distant server.")
if(NOT LIBSSH_FOUND)
  message(STATUS "NOTICE: The SSH features will be disabled.")
endif()

# Activate concurrency ? (turned OFF by default)
option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)

# And add -DCGAL_CONCURRENT_MESH_3 if that option is ON
if(CGAL_ACTIVATE_CONCURRENT_MESH_3 OR "$ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3}")
  add_compile_definitions(CGAL_CONCURRENT_MESH_3)
  if(NOT TBB_FOUND)
    find_package(TBB REQUIRED)
    include(CGAL_TBB_support)
    if(NOT TARGET CGAL::TBB_support)
      message(STATUS "NOTICE: Intel TBB was not found. Mesh_3 is faster if TBB is linked.")
    endif()
  endif()
else(CGAL_ACTIVATE_CONCURRENT_MESH_3 OR "$ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3}")
  option(LINK_WITH_TBB "Link with TBB anyway so we can use TBB timers for profiling" ON)
  if(LINK_WITH_TBB)
    find_package(TBB)
    include(CGAL_TBB_support)
  endif(LINK_WITH_TBB)
endif()

set_package_properties(
  TBB PROPERTIES
  DESCRIPTION "A library for parallel programming."
  PURPOSE "Plugins such as Mesh_3, Bilateral smoothing, and WLOP are faster if TBB is linked.")

if(CGAL_Qt6_FOUND AND Qt6_FOUND)
  qt6_wrap_ui(MainWindowUI_files MainWindow.ui)
  qt6_wrap_ui(SubViewerUI_files SubViewer.ui)
  qt6_wrap_ui(statisticsUI_FILES Statistics_on_item_dialog.ui)
  qt6_wrap_ui(FileLoaderDialogUI_files FileLoaderDialog.ui)
  qt6_wrap_ui(Show_point_dialogUI_FILES Show_point_dialog.ui)
  qt6_wrap_ui(PreferencesUI_FILES Preferences.ui Details.ui SSH_dialog.ui)
  qt6_wrap_ui(Show_point_dialogUI_FILES Show_point_dialog.ui)
  qt6_wrap_ui(ViewerUI_FILES LightingDialog.ui)
  qt6_generate_moc("File_loader_dialog.h"
                   "${CMAKE_CURRENT_BINARY_DIR}/File_loader_dialog_moc.cpp")

  include(${CMAKE_CURRENT_SOURCE_DIR}/CGALlab_macros.cmake)

  qt6_add_resources(CGAL_Qt6_RESOURCE_FILES CGALlab.qrc)
  find_path(
    CGAL_THREE_HEADERS_PATH
    NAMES CGAL/Three/Scene_item.h
    HINTS ${CGAL_INCLUDE_DIRS}
    NO_DEFAULT_PATH NO_CMAKE_FIND_ROOT_PATH
    DOC "Path to CGAL/Three/Scene_item.h")

  if(CGAL_THREE_HEADERS_PATH)
    qt6_generate_moc("${CGAL_THREE_HEADERS_PATH}/CGAL/Three/Viewer_interface.h"
                     "${CMAKE_CURRENT_BINARY_DIR}/Viewer_interface_moc.cpp")
    qt6_generate_moc("${CGAL_THREE_HEADERS_PATH}/CGAL/Three/Scene_item.h"
                     "${CMAKE_CURRENT_BINARY_DIR}/Scene_item_moc.cpp")
    qt6_generate_moc("${CGAL_THREE_HEADERS_PATH}/CGAL/Three/Scene_group_item.h"
                     "${CMAKE_CURRENT_BINARY_DIR}/Scene_group_item_moc.cpp")
    qt6_generate_moc(
      "${CGAL_THREE_HEADERS_PATH}/CGAL/Three/Scene_item_rendering_helper.h"
      "${CMAKE_CURRENT_BINARY_DIR}/Scene_item_rendering_helper_moc.cpp")
    qt6_generate_moc("${CGAL_THREE_HEADERS_PATH}/CGAL/Three/TextRenderer.h"
                     "${CMAKE_CURRENT_BINARY_DIR}/TextRenderer_moc.cpp")
  else()
    message(FATAL_ERROR "Cannot find <CGAL/Three/Scene_item.h>")
  endif()

  unset(CGAL_THREE_HEADERS_PATH CACHE)

  # AUXILIARY LIBRARIES

  # put plugins (which are shared libraries) at the same location as
  # executable files
  set(CGAL_LAB_DEMO_PLUGINS_DIR "${RUNTIME_OUTPUT_PATH}")
  set(LIBRARY_OUTPUT_PATH "${CGAL_LAB_DEMO_PLUGINS_DIR}")

  add_library(
    demo_framework SHARED
    Scene.cpp
    Viewer.cpp
    Three.cpp
    ${ViewerUI_FILES}
    ${CGAL_Qt6_RESOURCE_FILES}
    ${CGAL_Qt6_MOC_FILES}
    Viewer_interface_moc.cpp
    Scene_item_moc.cpp
    Scene_item.cpp
    Scene_item_with_properties.cpp
    Triangle_container.cpp
    Edge_container.cpp
    Point_container.cpp
    Scene_group_item.cpp
    Scene_group_item_moc.cpp
    TextRenderer.cpp
    TextRenderer_moc.cpp
    Scene_item_rendering_helper.cpp
    Scene_item_rendering_helper_moc.cpp
    Primitive_container.cpp
    CGALlab_plugin_helper.cpp
    CGAL_double_edit.cpp)
  target_link_libraries(demo_framework PRIVATE Qt6::OpenGLWidgets Qt6::Widgets Qt6::Gui
                                       PUBLIC CGAL::CGAL_Qt6)
  if(TARGET Qt6::WebSockets)
    target_link_libraries(demo_framework PRIVATE Qt6::WebSockets)
    message(STATUS "Qt6WebSockets was found. Using WebSockets is therefore possible.")
  endif()

  # This custom target is useless. It is used only as a flag to
  # detect that the test has already been created.
  add_custom_target("compilation_of__demo_framework")
  cgal_add_compilation_test(demo_framework)
  if(TEST "compilation of  demo_framework")
    set_property(TEST "compilation of  demo_framework"
      APPEND PROPERTY LABELS "CGAL_build_system" "Installation" "${PROJECT_NAME}")
    set_property(TEST "compilation of  demo_framework"
      APPEND PROPERTY FIXTURES_SETUP "demo_framework_SetupFixture")
    set_property(TEST "compilation of  demo_framework"
      APPEND PROPERTY DEPENDS "check build system" "compilation of  CGAL_Qt6_moc_and_resources")
  endif()

  # Let's define `three_EXPORT` during the compilation of `demo_framework`,
  # in addition of `demo_framework_EXPORT` (defined automatically by
  # CMake). That is to deal with the visibility of symbols of
  # `Three.h`/`Three.cpp`.
  target_compile_definitions(demo_framework PRIVATE three_EXPORTS=1)
  target_compile_definitions(demo_framework PRIVATE -DCGAL_USE_Qt6_RESOURCES)

  add_library(scene_basic_objects SHARED Scene_plane_item.cpp
                                         Scene_spheres_item.cpp)
  target_link_libraries(
    scene_basic_objects PRIVATE demo_framework CGAL::CGAL_Qt6 Qt6::OpenGLWidgets
                                Qt6::Gui Qt6::Widgets)
  add_library(scene_color_ramp SHARED Color_ramp.cpp)
  target_link_libraries(scene_color_ramp PRIVATE Qt6::Core)

  add_library(scene_callback_signaler SHARED Callback_signaler.cpp)
  target_link_libraries(scene_callback_signaler PRIVATE Qt6::Core)

  add_library(point_dialog SHARED Show_point_dialog.cpp Show_point_dialog.ui
                                  ${Show_point_dialogUI_FILES})
  target_link_libraries(point_dialog PRIVATE Qt6::OpenGLWidgets Qt6::Gui
                                            Qt6::Widgets)
  cgal_add_compilation_test(point_dialog)

  macro(add_item item_name)
    add_library(${item_name} SHARED ${ARGN})
    target_link_libraries(
      ${item_name} PRIVATE demo_framework CGAL::CGAL_Qt6 Qt6::OpenGLWidgets Qt6::Gui
                           Qt6::Widgets)
    cgal_add_compilation_test(${item_name})
    add_to_cached_list(CGAL_EXECUTABLE_TARGETS ${item_name})
    CGAL_install_hooks()
  endmacro(add_item)

  add_item(scene_triangulation_3_item Scene_triangulation_3_item.cpp)
  target_link_libraries(scene_triangulation_3_item PRIVATE scene_basic_objects scene_edit_box_item)

  add_item(scene_c3t3_item Scene_c3t3_item.cpp)
  target_link_libraries(
    scene_c3t3_item
    PUBLIC
    scene_triangulation_3_item
    PRIVATE
    scene_surface_mesh_item scene_polygon_soup_item
    scene_basic_objects ${TBB_LIBRARIES})

  add_item(scene_tetrahedra_item Scene_tetrahedra_item.cpp)
  target_link_libraries(scene_tetrahedra_item PRIVATE scene_c3t3_item)


  # Support for precompiled headers, for Mesh_3 (since CMake 3.16)
  add_library(c3t3_type_pch OBJECT c3t3_type_pch_empty_source.cpp C3t3_type.h)
  target_link_libraries(c3t3_type_pch PRIVATE CGAL::CGAL_Qt6 Qt6::OpenGLWidgets Qt6::Gui)
  set_property(TARGET c3t3_type_pch PROPERTY POSITION_INDEPENDENT_CODE TRUE)
  if(TARGET CGAL::TBB_support)
    target_link_libraries(c3t3_type_pch PRIVATE CGAL::TBB_support)
  endif()
  # target_precompile_headers(c3t3_type_pch PRIVATE [["C3t3_type.h"]])

  function(CGAL_Lab_target_use_c3t3_type target)
    target_link_libraries(${target} PRIVATE c3t3_type_pch)
    # target_precompile_headers(${target} REUSE_FROM c3t3_type_pch)
    if(TARGET CGAL::TBB_support)
      target_link_libraries(${target} PRIVATE CGAL::TBB_support)
    endif()
  endfunction(CGAL_Lab_target_use_c3t3_type)

  CGAL_Lab_target_use_c3t3_type(scene_c3t3_item)
  CGAL_Lab_target_use_c3t3_type(scene_triangulation_3_item)
  CGAL_Lab_target_use_c3t3_type(scene_tetrahedra_item)

  add_item(scene_aff_transformed_item Plugins/PCA/Scene_aff_transformed_item.cpp)

  add_item(scene_aff_transformed_point_set_item Plugins/PCA/Scene_aff_transformed_point_set_item.cpp)
  target_link_libraries(scene_aff_transformed_point_set_item PRIVATE scene_points_with_normal_item
                                                                     scene_aff_transformed_item)
  add_item(scene_aff_transformed_polygon_soup_item Plugins/PCA/Scene_aff_transformed_polygon_soup_item.cpp)
  target_link_libraries(scene_aff_transformed_polygon_soup_item PRIVATE scene_polygon_soup_item
                                                                        scene_aff_transformed_item)
  add_item(scene_aff_transformed_surface_mesh_item Plugins/PCA/Scene_aff_transformed_surface_mesh_item.cpp)
  target_link_libraries(scene_aff_transformed_surface_mesh_item PRIVATE scene_surface_mesh_item
                                                                        scene_aff_transformed_item)

  add_item(scene_edit_box_item Plugins/PCA/Scene_edit_box_item.cpp)

  add_item(scene_image_item Scene_image_item.cpp)
  add_item(scene_surface_mesh_item Scene_surface_mesh_item.cpp)

  if(TBB_FOUND)
    target_link_libraries(scene_surface_mesh_item PRIVATE CGAL::TBB_support)
  endif()
  # special

  add_item(scene_item_decorator Scene_polyhedron_item_decorator.cpp)
  target_link_libraries(scene_item_decorator PRIVATE scene_surface_mesh_item)

  add_item(scene_k_ring_selection
           Plugins/PMP/Scene_facegraph_item_k_ring_selection.cpp)
  target_link_libraries(scene_k_ring_selection PRIVATE scene_surface_mesh_item)

  add_item(scene_selection_item Scene_polyhedron_selection_item.cpp)
  target_link_libraries(scene_selection_item
    PUBLIC
    scene_surface_mesh_item
    scene_item_decorator
    scene_k_ring_selection)
  if(TBB_FOUND)
    target_link_libraries(scene_selection_item PRIVATE CGAL::TBB_support)
  endif()

  add_item(scene_shortest_path_item
           Plugins/Surface_mesh/Scene_polyhedron_shortest_path_item.cpp)
  target_link_libraries(
    scene_shortest_path_item
    PRIVATE scene_item_decorator scene_surface_mesh_item scene_polylines_item)

  add_item(scene_movable_sm_item Plugins/AABB_tree/Scene_movable_sm_item.cpp)

  if(TARGET CGAL::Eigen3_support)
    add_item(scene_textured_item Scene_textured_surface_mesh_item.cpp
             texture.cpp)
    target_link_libraries(scene_textured_item PRIVATE CGAL::Eigen3_support)
    add_item(scene_mcf_item Plugins/PMP/Scene_mcf_item.cpp)
    target_link_libraries(scene_mcf_item PRIVATE CGAL::Eigen3_support)
  endif()

  add_item(scene_implicit_function_item Scene_implicit_function_item.cpp)
  target_link_libraries(scene_implicit_function_item PRIVATE scene_color_ramp)

  add_item(scene_polygon_soup_item Scene_polygon_soup_item.cpp)
  target_link_libraries(scene_polygon_soup_item PRIVATE scene_surface_mesh_item)

  add_item(scene_nef_polyhedron_item Scene_nef_polyhedron_item.cpp)
  target_link_libraries(scene_nef_polyhedron_item
                        PRIVATE scene_surface_mesh_item)

  add_item(scene_points_with_normal_item Scene_points_with_normal_item.cpp)
  if(TARGET CGAL::Eigen3_support)
    target_link_libraries(scene_points_with_normal_item
                          PRIVATE CGAL::Eigen3_support)
  endif()

  if(TBB_FOUND)
    target_link_libraries(scene_points_with_normal_item
                          PRIVATE CGAL::TBB_support)
  endif()

  add_item(scene_polylines_item Scene_polylines_item.cpp)
  target_link_libraries(
    scene_polylines_item PRIVATE scene_basic_objects
                                 scene_points_with_normal_item)

  add_item(scene_lcc_item Scene_lcc_item.cpp)

  foreach(lib demo_framework scene_basic_objects scene_color_ramp
              scene_polygon_soup_item scene_nef_polyhedron_item)
    add_to_cached_list(CGAL_EXECUTABLE_TARGETS ${lib})
  endforeach()

  add_compile_definitions(USE_FORWARD_DECL)
  add_library(
    cgal_lab SHARED
    MainWindow.cpp
    CGAL_Lab.cpp
    File_loader_dialog_moc.cpp
    Use_ssh.cpp
    ${CGAL_Qt6_RESOURCE_FILES}
    ${CGAL_Qt6_MOC_FILES}
    ${FileLoaderDialogUI_files}
    ${MainWindowUI_files}
    ${PreferencesUI_FILES}
    ${statisticsUI_FILES}
    ${SubViewerUI_files})
  target_link_libraries(
    cgal_lab PRIVATE demo_framework point_dialog Qt6::Widgets
             PUBLIC CGAL::CGAL Qt6::Gui Qt6::OpenGLWidgets Qt6::Qml
                           )
  if(LIBSSH_FOUND)
    target_compile_definitions(cgal_lab PRIVATE -DCGAL_USE_SSH)
    target_link_libraries(cgal_lab PRIVATE ${LIBSSH_LIBRARIES})
    target_include_directories(cgal_lab SYSTEM PRIVATE ${LIBSSH_INCLUDE_DIR})
  endif() #libssh
  if(TARGET Qt6::WebSockets)
    target_compile_definitions(cgal_lab PRIVATE -DCGAL_USE_WEBSOCKETS)
    target_compile_definitions(demo_framework PRIVATE -DCGAL_USE_WEBSOCKETS)
    target_link_libraries(cgal_lab PRIVATE Qt6::WebSockets)
  endif()
  cgal_add_compilation_test(cgal_lab)

  #
  # the `CGALlab` executable
  #   required as a dependency by all plugins
  #
  add_executable(CGALlab CGALlab.cpp)
  target_link_libraries(CGALlab PRIVATE cgal_lab)
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS CGALlab)
  cgal_add_compilation_test(CGALlab)

  #
  # PLUGINS
  #

  file(
    GLOB DEMO_PLUGINS
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/Plugins/"
    "${CMAKE_CURRENT_SOURCE_DIR}/Plugins/*")
  foreach(SUB_DIR ${DEMO_PLUGINS})
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Plugins/${SUB_DIR}")
  endforeach()

  add_subdirectory(implicit_functions)

  add_dependencies(CDT_3  IO_surface_meshes detect_sharp_edges_plugin polylines_io_plugin selection_plugin)
  add_dependencies(Mesh_3 IO_surface_meshes)

  #
  # EXECUTABLES
  #
  add_executable(CGAL_Mesh_3 Mesh_3.cpp)
  add_dependencies(CGAL_Mesh_3 Mesh_3 IO_surface_meshes)
  target_link_libraries(CGAL_Mesh_3 PRIVATE cgal_lab)
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS CGAL_Mesh_3)

  add_executable(CGAL_CDT_3 CGAL_Lab_CDT_3.cpp)
  add_dependencies(CGAL_CDT_3 CDT_3 IO_surface_meshes)
  target_link_libraries(CGAL_CDT_3 PRIVATE cgal_lab)
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS CGAL_CDT_3)

  add_executable(CGAL_Classification Classification.cpp)
  add_dependencies(CGAL_Classification Classification)
  target_link_libraries(CGAL_Classification PRIVATE cgal_lab)
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS CGAL_Classification)

  add_executable(CGAL_PMP PMP.cpp)
  add_dependencies(CGAL_PMP PMP)
  target_link_libraries(CGAL_PMP PRIVATE cgal_lab)
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS CGAL_PMP)

  #WS Server
  if(TARGET Qt6::WebSockets AND TARGET Qt6::Network)
    add_executable(WS_server Server_ws.cpp)
    target_link_libraries(WS_server PRIVATE Qt6::WebSockets Qt6::Widgets Qt6::Network)
    message(STATUS "Qt6WebSockets was found. Using WebSockets is therefore possible.")
  endif()
  #
  # Exporting
  #
  if(TARGET CGAL_Qt6)
    export(
      TARGETS CGAL CGAL_Qt6 CGAL_ImageIO
      FILE cgal_lab_targets.cmake
      NAMESPACE Lab_)
  endif()
  if(TARGET CGAL_Qt6_moc_and_resources)
    export(
      TARGETS CGAL_Qt6_moc_and_resources
      NAMESPACE Lab_ APPEND
      FILE cgal_lab_targets.cmake)
  endif()
  export(
    TARGETS demo_framework
            scene_surface_mesh_item
            scene_points_with_normal_item
            scene_color_ramp
            scene_implicit_function_item
            scene_polylines_item
            scene_basic_objects
            scene_selection_item
            scene_item_decorator
            scene_k_ring_selection
    NAMESPACE Lab_ APPEND
    FILE cgal_lab_targets.cmake)

  configure_file(CGALlabConfig.cmake.in
                 CGALlabConfig.cmake)
  #TO DO script the activation of all the plugins.

else(CGAL_Qt6_FOUND AND Qt6_FOUND)

  set(Lab_MISSING_DEPS "")

  if(NOT CGAL_Qt6_FOUND)
    set(Lab_MISSING_DEPS "the CGAL Qt6 library, ${Lab_MISSING_DEPS}")
  endif()

  if(NOT Qt6_FOUND)
    set(Lab_MISSING_DEPS "Qt6, ${Lab_MISSING_DEPS}")
  endif()

  message("NOTICE: XX This demo requires ${Lab_MISSING_DEPS} and will not be compiled.")

endif(CGAL_Qt6_FOUND AND Qt6_FOUND)

feature_summary(
  WHAT REQUIRED_PACKAGES_NOT_FOUND
  INCLUDE_QUIET_PACKAGES
  DESCRIPTION "NOTICE: Missing required packages that prevent the demo from being compiled:"
  QUIET_ON_EMPTY
  VAR NotFound_REQ_PACKAGES)
if(NOT ${NotFound_REQ_PACKAGES} STREQUAL "")
  message("${NotFound_REQ_PACKAGES}")
endif()

feature_summary(
  WHAT OPTIONAL_PACKAGES_NOT_FOUND INCLUDE_QUIET_PACKAGES
  DESCRIPTION "Missing optional packages:"
  QUIET_ON_EMPTY
  VAR NotFound_PACKAGES)
if(NOT ${NotFound_PACKAGES} STREQUAL "")
  message(STATUS "${NotFound_PACKAGES}")
endif()

feature_summary(
  WHAT DISABLED_FEATURES
  DESCRIPTION "Missing Packages Components:"
  QUIET_ON_EMPTY
  VAR Disabled_Features)
if(NOT ${Disabled_Features} STREQUAL "")
  message(STATUS "${Disabled_Features}")
endif()
