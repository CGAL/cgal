cmake_minimum_required(VERSION 3.1...3.23)
project(Property_map_Tests)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("test_property_map.cpp")
create_single_source_cgal_program("dynamic_property_map.cpp")
create_single_source_cgal_program("dynamic_properties_test.cpp")
create_single_source_cgal_program("kernel_converter_properties_test.cpp")
create_single_source_cgal_program("test_Property_container.cpp")

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  message(STATUS "Found OpenMesh")
  include(UseOpenMesh)

  target_link_libraries(dynamic_properties_test PRIVATE ${OPENMESH_LIBRARIES})
  target_compile_definitions(dynamic_properties_test PRIVATE -DCGAL_USE_OPENMESH)
else()
  message(STATUS "NOTICE: Tests will not use OpenMesh.")
endif()
