# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.
# cmake ../ -DCMAKE_BUILD_TYPE=Debug

cmake_minimum_required(VERSION 3.1...3.15)
project(Linear_cell_complex_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()
if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

## To add expensive tests
# add_definitions("-DCGAL_CHECK_EXPENSIVE")
# add_definitions("-Wall -Wextra")

## For profilling with gprof
# add_definitions("-pg")
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

## To add an is_valid test after each operation (only in debug mode)
# add_definitions("-DCGAL_CMAP_TEST_VALID_REMOVALS")
# add_definitions("-DCGAL_CMAP_TEST_VALID_CONTRACTIONS")
# add_definitions("-DCGAL_CMAP_TEST_VALID_INSERTIONS")

# Option allowing to profile each operation of the demo (cout times on stdout).
add_definitions(-DCGAL_PROFILE_LCC_DEMO)

##################
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)

if(NOT
   (CGAL_Qt5_FOUND
    AND Qt5_FOUND))

  message(STATUS "NOTICE: This demo requires CGAL, "
                 "and Qt5, and will not be compiled.")

else()

  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)

  # ui file, created wih Qt Designer
  qt5_wrap_ui(uis MainWindow.ui CreateMesh.ui CreateMenger.ui
              CreateSierpinskiCarpet.ui CreateSierpinskiTriangle.ui)

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES ./Linear_cell_complex_3.qrc)

  add_executable(
    Linear_cell_complex_3_demo
    Linear_cell_complex_3_demo.cpp
    MainWindow.cpp
    Viewer.cpp
    Linear_cell_complex_3_subdivision.cpp
    Linear_cell_complex_pqq_subdivision.cpp
    typedefs.h
    import_moka.h
    ${uis}
    ${CGAL_Qt5_RESOURCE_FILES}
    ${CGAL_Qt5_MOC_FILES})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Linear_cell_complex_3_demo)

  target_link_libraries(Linear_cell_complex_3_demo
                        PUBLIC CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Gui Qt5::OpenGL)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Linear_cell_complex_3_demo)
endif()
