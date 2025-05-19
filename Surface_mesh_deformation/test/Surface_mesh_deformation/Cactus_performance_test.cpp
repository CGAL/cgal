#include "Surface_mesh_deformation_test_commons.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh_deformation.h>

typedef CGAL::Simple_cartesian<double> Kernel;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>  Polyhedron;

typedef CGAL::Surface_mesh_deformation<Polyhedron, CGAL::Default, CGAL::Default,
  CGAL::ORIGINAL_ARAP> Deform_mesh_arap;

typedef CGAL::Surface_mesh_deformation<Polyhedron, CGAL::Default, CGAL::Default,
  CGAL::SPOKES_AND_RIMS> Deform_mesh_spoke;

int main()
{
  Polyhedron mesh_1;
  read_to_polyhedron(CGAL::data_file_path("meshes/cactus.off"), mesh_1);
  Polyhedron mesh_2 = mesh_1;

  init_indices(mesh_1);
  init_indices(mesh_2);

  Deform_mesh_arap deform_mesh_arap(mesh_1);
  Deform_mesh_spoke deform_mesh_spoke(mesh_2);

  const int deformation_iteration = 500;
  const double x = -0.55; const double y = -0.50; const double z = -0.0;

  std::cerr << "ORIGINAL_ARAP performance: " << std::endl;
  preprocess_and_deform(deform_mesh_arap,
    "data/cactus_roi.txt",
    "data/cactus_handle.txt",
    CGAL::Simple_cartesian<double>::Vector_3(x, y, z),
    deformation_iteration);

  std::cerr << "SPOKES_AND_RIMS performance: " << std::endl;
  preprocess_and_deform(deform_mesh_spoke,
    "data/cactus_roi.txt",
    "data/cactus_handle.txt",
    CGAL::Simple_cartesian<double>::Vector_3(x, y, z),
    deformation_iteration);

  std::cerr << "Save deformed models" << std::endl;
  std::ofstream output("data/cactus_deformed_arap.off");
  output << mesh_1;
  output.close();
  output.open("data/cactus_deformed_spokes.off");
  output << mesh_2;
  std::cerr << "All done!" << std::endl;
  return EXIT_SUCCESS;
}

