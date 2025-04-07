cmake_minimum_required(VERSION 3.12...3.31)
project(Periodic_2_triangulation_2_Demo)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets)
include_directories(BEFORE ./include)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  add_compile_definitions(QT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)

  # UI files (Qt Designer files)
  qt6_wrap_ui(DT_UI_FILES Periodic_2_triangulation_2.ui)

  # qrc files (resources files, that contain icons, at least)
  qt6_add_resources(CGAL_Qt6_RESOURCE_FILES ./Periodic_2_triangulation_2.qrc)

  # find header files for projects that can show them
  file(GLOB headers "*.h")
  file(GLOB QT_headers "include/CGAL/Qt/*.h")
  file(GLOB P2T2_headers "../../../include/CGAL/*.h")
  source_group("P2T2" FILES ${P2T2_headers})
  source_group("QT" FILES ${QT_headers})

  # The executable itself.
  qt_add_executable(
    Periodic_2_Delaunay_triangulation_2
    Periodic_2_Delaunay_triangulation_2.cpp
    ${DT_UI_FILES}
    ${CGAL_Qt6_RESOURCE_FILES}
    ${CGAL_Qt6_MOC_FILES}
    ${headers}
    ${QT_headers}
    ${P2T2_headers})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS
                     Periodic_2_Delaunay_triangulation_2)

  target_link_libraries(Periodic_2_Delaunay_triangulation_2
                        PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Periodic_2_Delaunay_triangulation_2)
else()

  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")

endif()
