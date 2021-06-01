cmake_minimum_required(VERSION 3.1...3.20)
project(Surface_mesh_shortest_path_Tests)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)

include(${CGAL_USE_FILE})

find_package(LEDA QUIET)

include_directories(BEFORE "include")

create_single_source_cgal_program("Surface_mesh_shortest_path_test_1.cpp")
create_single_source_cgal_program("Surface_mesh_shortest_path_test_2.cpp")
create_single_source_cgal_program("Surface_mesh_shortest_path_test_3.cpp")
create_single_source_cgal_program("Surface_mesh_shortest_path_test_4.cpp")
create_single_source_cgal_program("Surface_mesh_shortest_path_test_5.cpp")
create_single_source_cgal_program("Surface_mesh_shortest_path_test_6.cpp")
create_single_source_cgal_program("Surface_mesh_shortest_path_traits_test.cpp")

# Link with Boost.ProgramOptions (optional)
find_package(Boost QUIET COMPONENTS program_options)
if(Boost_PROGRAM_OPTIONS_FOUND)
  if(CGAL_Core_FOUND OR LEDA_FOUND)
    create_single_source_cgal_program("TestMesh.cpp")
    if(TARGET Boost::filesystem)
      target_link_libraries(TestMesh PRIVATE Boost::program_options)
    else()
      target_link_libraries(TestMesh PRIVATE ${Boost_PROGRAM_OPTIONS_LIBRARY})
    endif()
  else()
    message(
      STATUS
        "NOTICE: Example TestMesh.cpp requires the CGAL_Core library (or LEDA) and will not be compiled."
    )
  endif()
else()
  message(
    STATUS
      "NOTICE: Example TestMesh.cpp requires boost program_option and will not be compiled."
  )
endif()
