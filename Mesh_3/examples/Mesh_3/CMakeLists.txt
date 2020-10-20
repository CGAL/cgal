cmake_minimum_required(VERSION 3.1...3.15)
project(Mesh_3_Examples)

add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if(CGAL_MESH_3_VERBOSE)
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

find_package(CGAL REQUIRED COMPONENTS ImageIO)
find_package(Boost)

option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)

if(CGAL_ACTIVATE_CONCURRENT_MESH_3 OR ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3})
  add_definitions(-DCGAL_CONCURRENT_MESH_3)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
endif()

# Use Eigen
find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include(CGAL_Eigen_support)

find_package(VTK QUIET COMPONENTS vtkImagingGeneral vtkIOImage NO_MODULE)
if(VTK_FOUND)
  if(VTK_USE_FILE)
    include(${VTK_USE_FILE})
  endif()
  if("${VTK_VERSION_MAJOR}" GREATER "5" OR VTK_VERSION VERSION_GREATER 5)
    message(STATUS "VTK found")
    if(TARGET VTK::IOImage)
      set(VTK_LIBRARIES VTK::ImagingGeneral VTK::IOImage)
    endif()
  else()
    message(STATUS "VTK version 6.0 or greater is required")
  endif()
else()
  message(STATUS "VTK was not found")
endif()

create_single_source_cgal_program("mesh_hybrid_mesh_domain.cpp")
target_link_libraries(mesh_hybrid_mesh_domain PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_implicit_sphere.cpp")
target_link_libraries(mesh_implicit_sphere PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_implicit_sphere_variable_size.cpp")
target_link_libraries(mesh_implicit_sphere_variable_size
                      PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_two_implicit_spheres_with_balls.cpp")
target_link_libraries(mesh_two_implicit_spheres_with_balls
                      PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_implicit_domains_2.cpp"
                                  "implicit_functions.cpp")
target_link_libraries(mesh_implicit_domains_2 PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_cubes_intersection.cpp")
target_link_libraries(mesh_cubes_intersection PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_cubes_intersection_with_features.cpp")
target_link_libraries(mesh_cubes_intersection_with_features
                      PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_implicit_domains.cpp"
                                  "implicit_functions.cpp")
target_link_libraries(mesh_implicit_domains PUBLIC CGAL::Eigen_support)

find_package( ITT QUIET )
if(TARGET ITT::ITT)
  target_link_libraries(mesh_implicit_domains PRIVATE ITT::ITT)
  target_compile_definitions(mesh_implicit_domains PRIVATE CGAL_MESH_3_USE_INTEL_ITT)
endif()

create_single_source_cgal_program( "mesh_polyhedral_domain.cpp" )
target_link_libraries(mesh_polyhedral_domain PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_polyhedral_domain_sm.cpp")
target_link_libraries(mesh_polyhedral_domain_sm PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program(
  "mesh_polyhedral_domain_with_surface_inside.cpp")
target_link_libraries(mesh_polyhedral_domain_with_surface_inside
                      PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("remesh_polyhedral_surface.cpp")
target_link_libraries(remesh_polyhedral_surface PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("remesh_polyhedral_surface_sm.cpp")
target_link_libraries(remesh_polyhedral_surface_sm PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_polyhedral_domain_with_features.cpp")
target_link_libraries(mesh_polyhedral_domain_with_features
                      PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_polyhedral_domain_with_features_sm.cpp")
target_link_libraries(mesh_polyhedral_domain_with_features_sm
                      PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program(
  "mesh_polyhedral_domain_with_lipschitz_sizing.cpp")
target_link_libraries(mesh_polyhedral_domain_with_lipschitz_sizing
                      PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_polyhedral_complex.cpp")
target_link_libraries(mesh_polyhedral_complex PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("mesh_polyhedral_complex_sm.cpp")
target_link_libraries(mesh_polyhedral_complex_sm PUBLIC CGAL::Eigen_support)

if(TARGET CGAL::CGAL_ImageIO)
  if(VTK_FOUND AND ("${VTK_VERSION_MAJOR}" GREATER "5" OR VTK_VERSION
                                                          VERSION_GREATER 5))
    add_executable(mesh_3D_gray_vtk_image mesh_3D_gray_vtk_image.cpp)
    target_link_libraries(
      mesh_3D_gray_vtk_image
      PUBLIC CGAL::Eigen_support CGAL::CGAL CGAL::CGAL_ImageIO
             ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${VTK_LIBRARIES})
    cgal_add_test(mesh_3D_gray_vtk_image)
    add_to_cached_list(CGAL_EXECUTABLE_TARGETS mesh_3D_gray_vtk_image)
  endif()

  create_single_source_cgal_program("mesh_3D_gray_image.cpp")
  target_link_libraries(mesh_3D_gray_image PUBLIC CGAL::Eigen_support)

  create_single_source_cgal_program("mesh_3D_gray_image_multiple_values.cpp")
  target_link_libraries(mesh_3D_gray_image_multiple_values
                        PUBLIC CGAL::Eigen_support)

  create_single_source_cgal_program("mesh_3D_image_with_features.cpp")
  target_link_libraries(mesh_3D_image_with_features PUBLIC CGAL::Eigen_support)

  if(CGAL_ImageIO_USE_ZLIB)
    create_single_source_cgal_program("mesh_optimization_example.cpp")
    target_link_libraries(mesh_optimization_example PUBLIC CGAL::Eigen_support)

    create_single_source_cgal_program("mesh_optimization_lloyd_example.cpp")
    target_link_libraries(mesh_optimization_lloyd_example
                          PUBLIC CGAL::Eigen_support)

    create_single_source_cgal_program("mesh_3D_image.cpp")
    target_link_libraries(mesh_3D_image PUBLIC CGAL::Eigen_support)

    create_single_source_cgal_program(
      "mesh_3D_image_with_custom_initialization.cpp")
    target_link_libraries(mesh_3D_image_with_custom_initialization
                          PUBLIC CGAL::Eigen_support)

    create_single_source_cgal_program("mesh_3D_image_variable_size.cpp")
    target_link_libraries(mesh_3D_image_variable_size
                          PUBLIC CGAL::Eigen_support)
  else()
    message(
      STATUS
        "NOTICE: The examples mesh_3D_image.cpp, mesh_3D_image_variable_size.cpp, mesh_optimization_example.cpp and mesh_optimization_lloyd_example.cpp need CGAL_ImageIO to be configured with ZLIB support, and will not be compiled."
    )
  endif()

else()
  message(
    STATUS
      "NOTICE: Some examples need the CGAL_ImageIO library, and will not be compiled."
  )
endif()

if(CGAL_ACTIVATE_CONCURRENT_MESH_3 AND TARGET CGAL::TBB_support)
  foreach(
    target
    mesh_3D_image
    mesh_3D_image_variable_size
    mesh_3D_image_with_custom_initialization
    mesh_3D_image_with_features
    mesh_implicit_sphere
    mesh_implicit_sphere_variable_size
    mesh_optimization_example
    mesh_optimization_lloyd_example
    mesh_polyhedral_complex
    mesh_polyhedral_complex_sm
    mesh_polyhedral_domain
    mesh_polyhedral_domain_sm
    mesh_polyhedral_domain_with_features
    mesh_polyhedral_domain_with_features_sm
    mesh_polyhedral_domain_with_lipschitz_sizing
    mesh_two_implicit_spheres_with_balls)
    if(TARGET ${target})
      target_link_libraries(${target} PUBLIC CGAL::TBB_support)
    endif()
  endforeach()
endif()
