# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.
# cmake ../ -DCMAKE_BUILD_TYPE=Debug

cmake_minimum_required(VERSION 3.12...3.31)
project(Linear_cell_complex_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

## To add expensive tests
# add_compile_definitions(CGAL_CHECK_EXPENSIVE)
# add_definitions("-Wall -Wextra")
# add_definitions(-fsanitize=address)
# add_link_options(-fsanitize=address)
# add_compile_definitions(_GLIBCXX_DEBUG)

## For profilling with gprof
# add_compile_definitions("-pg")
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

## To add an is_valid test after each operation (only in debug mode)
# add_compile_definitions(CGAL_CMAP_TEST_VALID_REMOVALS)
# add_compile_definitions(CGAL_CMAP_TEST_VALID_CONTRACTIONS)
# add_compile_definitions(CGAL_CMAP_TEST_VALID_INSERTIONS)

# Option allowing to profile each operation of the demo (cout times on stdout).
add_compile_definitions(CGAL_PROFILE_LCC_DEMO)

add_compile_definitions(CMAP_WITH_INDEX) # to use cc with index (handle otherwise)

##################
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 COMPONENTS Widgets OpenGL)

if(NOT CGAL_Qt6_FOUND OR NOT Qt6_FOUND)

  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")

else()

  add_compile_definitions(CGAL_USE_BASIC_VIEWER QT_NO_KEYWORDS)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  qt_add_executable(
    Linear_cell_complex_3_demo
    Linear_cell_complex_3_demo.cpp
    Linear_cell_complex_3.qrc
    MainWindow.cpp
    Viewer.cpp
    Linear_cell_complex_3_subdivision.cpp
    Linear_cell_complex_pqq_subdivision.cpp
    typedefs.h
    import_moka.h
    MainWindow.h
    Viewer.h)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Linear_cell_complex_3_demo)

  target_link_libraries(Linear_cell_complex_3_demo
                        PRIVATE CGAL::CGAL CGAL::CGAL_Basic_viewer)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Linear_cell_complex_3_demo)
endif()
