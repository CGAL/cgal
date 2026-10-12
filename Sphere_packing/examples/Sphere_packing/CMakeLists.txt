cmake_minimum_required(VERSION 3.12...3.31)
project(Sphere_packing_Examples)

find_package(CGAL REQUIRED)
find_package(ProtoSphere)
include(CGAL_ProtoSphere_support)

if (TARGET CGAL::ProtoSphere)
  create_single_source_cgal_program("sphere_packing.cpp")
  target_link_libraries(sphere_packing PRIVATE CGAL::ProtoSphere_support )

  # Windows: copy the ProtoSphere DLL (Debug or Release, matching the build) next to the executable
  if(WIN32)
    add_custom_command(TARGET sphere_packing POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              $<TARGET_FILE:ProtoSphere::ProtoSphereExtended>
              $<TARGET_FILE_DIR:sphere_packing>)
  endif()
else()
  message(STATUS "NOTICE: The Sphere_packing package requires the ProtoSphere library will not be compiled.")
endif()
