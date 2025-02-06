cmake_minimum_required(VERSION 3.12...3.29)


if(NOT CGAL_MODULES_DIR)
  find_package(CGAL REQUIRED)
endif()
include(${CGAL_MODULES_DIR}/CGALHelpers.cmake)
CGAL_handle_subdirectories(test tests)
