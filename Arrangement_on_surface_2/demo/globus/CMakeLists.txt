# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Globus_demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# General
set(GLOBUS_MODULES_REL_DIR cmake/modules)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/${GLOBUS_MODULES_REL_DIR})

find_package(CGAL QUIET COMPONENTS Qt5 OPTIONAL_COMPONENTS Core)
find_package(Qt5 QUIET COMPONENTS Gui Widgets)
# find_package(FileGDBAPI REQUIRED)

################################################################################
include(FetchContent)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
FetchContent_MakeAvailable(json)
# find_package(nlohmann_json 3.2.0 REQUIRED)
################################################################################

if (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND)
  include(${CGAL_USE_FILE})
  add_compile_definitions(QT_NO_KEYWORDS)
  include_directories( BEFORE ./ )
  # include_directories(${FileGDBAPI_INCLUDE_DIR})

  # Arrangement package includes
  add_definitions(-DQT_NO_KEYWORDS)
  option(COMPILE_UTILS_INCREMENTALLY
    "Compile files in Utils directory incrementally, or compile them all as a unit. \
    Incremental compilation will be better for development and consume less \
    memory while compiling but will take longer to compile."
    OFF)

  set(UTILS_SOURCE_FILES "")

  if (COMPILE_UTILS_INCREMENTALLY)
    set(UTILS_COMPILE_FILES ${UTILS_SOURCE_FILES})
  else()
    set(UTILS_CPP_FILES_INCLUDES "")
    foreach(utils_src IN LISTS UTILS_SOURCE_FILES)
      string(APPEND UTILS_CPP_FILES_INCLUDES "#include \"${utils_src}\"\n")
    endforeach()
    file(WRITE "${CMAKE_BINARY_DIR}/CombinedUtils.cpp" ${UTILS_CPP_FILES_INCLUDES})
    set(UTILS_COMPILE_FILES "${CMAKE_BINARY_DIR}/CombinedUtils.cpp")
  endif()


  qt5_wrap_ui(globus_uis
              Globus_window.ui)

  qt5_wrap_cpp(CGAL_Qt5_MOC_FILES
               Globus_window.h)

  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES Globus.qrc)

  add_executable(globus
                 globus.cpp
                 Globus_window.cpp
                 ${UTILS_COMPILE_FILES}
                 ${globus_uis}
                 ${CGAL_Qt5_RESOURCE_FILES}
                 ${CGAL_Qt5_MOC_FILES})

  # target_link_libraries(globus ${FileGDBAPI_LIBRARY})
  target_link_libraries(globus Qt5::Core Qt5::Gui Qt5::Widgets)
  target_link_libraries(globus CGAL::CGAL CGAL::CGAL_Qt5)
  target_link_libraries(globus nlohmann_json::nlohmann_json)
  if(CGAL_Core_FOUND)
    target_link_libraries(globus CGAL::CGAL_Core)
  endif()

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS globus)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(globus)

  # set_property(TARGET globus PROPERTY INSTALL_RPATH "${FileGDBAPI_LIBRARY_DIR}")
  set_property(TARGET globus PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

  set_target_properties(globus PROPERTIES LINK_FLAGS "-Wl,--disable-new-dtags")

else()
  set(MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(MISSING_DEPS "CGAL, ${MISSING_DEPS}")
  endif()
  if(NOT CGAL_Qt5_FOUND)
    set(MISSING_DEPS "the CGAL Qt5 library, ${MISSING_DEPS}")
  endif()
  if(NOT Qt5_FOUND)
    set(MISSING_DEPS "Qt5, ${MISSING_DEPS}")
  endif()
  message(STATUS
          "NOTICE: This demo requires ${MISSING_DEPS} and will not be compiled.")
endif()
