// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex/IO/VTK.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/Mesh_smoothing_3/mesh_representations.h>
#include <CGAL/Mesh_smoothing_3/Mesh_smoothing_3.h>
#include <CGAL/Mesh_smoothing_3/Projectors.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using LCC = CGAL::Linear_cell_complex_for_combinatorial_map<3, 3>;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_3<Surface_mesh, K>;

namespace {

[[maybe_unused]] void test_lcc_mixed_element_smooth()
{
  std::cout << std::string(MS3_EXAMPLE_MESH_DIRECTORY) + "volume.vtk" << std::endl;
  std::string filename = std::string(MS3_EXAMPLE_MESH_DIRECTORY) + "volume.vtk";
  std::string boundary_filename = std::string(MS3_EXAMPLE_MESH_DIRECTORY) + "target.obj";
  LCC lcc;

  [[maybe_unused]] bool res = CGAL::IO::read_VTK(filename.c_str(), lcc);
  assert(res);

  assert(lcc.is_valid());

  std::cout << "Loaded LCC: ";
  lcc.display_characteristics(std::cout);
  std::cout << '\n';

  CGAL::IO::write_VTK("lcc_initial.vtk", lcc);


  Surface_mesh poly;
  CGAL::IO::read_polygon_mesh(boundary_filename, poly);
  CGAL::Polygon_mesh_processing::triangulate_faces(poly);
  Mesh_domain mesh_domain(poly);
  CGAL::Mesh_smoothing_3::Polyhedral_mesh_domain_projector<Mesh_domain> projector(mesh_domain);

  using ConcurrencyTag = CGAL::Mesh_smoothing_3::Parallel_if_available_tag;
  unsigned max_nb_metric_evaluations = 10000;
  double time_limit = 0.;
  bool verbose = true;


  using LCC_mixed_mesh = CGAL::Mesh_smoothing_3::cgal_types::LCC_mixed_mesh<LCC>;
  using Mixed_mesh = CGAL::Mesh_smoothing_3::helper_structures::Mixed_mesh_wrapper<LCC_mixed_mesh>;

  LCC_mixed_mesh lcc_mixed_mesh(lcc);
  Mixed_mesh mixed_mesh(lcc_mixed_mesh);

  using Lcc_surface = CGAL::Mesh_smoothing_3::cgal_types::LCC_surface_mesh<LCC>;
  Lcc_surface surface(lcc);

  using Curve_network = CGAL::Mesh_smoothing_3::default_structures::Empty_edge_network<Mixed_mesh::Vertex_descriptor>;

  CGAL::Mesh_smoothing_3::Mesh_smoother<Mixed_mesh, Lcc_surface, Curve_network, ConcurrencyTag> 
  smoother(mixed_mesh, surface);

  smoother.set_predicates_mode(CGAL::Mesh_smoothing_3::Parameters::STRONG_ENFORCEMENT);
  smoother.set_verbose(verbose);
  smoother.set_maximum_running_time(time_limit);
  smoother.set_maximum_number_of_metric_evaluations(max_nb_metric_evaluations);

  smoother.set_boundary_query([&](std::vector<Point_3> const& pts, Lcc_surface::Surface_patch_index patch_face) {
      auto proj = projector.patch_face_projection_plane(patch_face, pts);
      return std::make_tuple(proj.origin(), proj.vector(), 1.);
  });

  auto result = smoother.run();

  std::cout << "Number of inverted elements: " << result.nb_invalid_elements << std::endl;
  std::cout << "Number of vertex updates: " << result.nb_vertex_updates << std::endl;
  std::cout << "Number of metric evaluations: " << result.nb_metric_evaluations << std::endl;
  std::cout << "Pre-processing time: " << result.pre_processing_time << std::endl;
  std::cout << "Smoothing time: " << result.optimization_time << std::endl;

  CGAL::IO::write_VTK("lcc_smoothed.vtk", lcc);

  std::cout << "Done" << std::endl;
}

} // namespace

int main()
{
  // not run because the input mesh is too large
  // This test is used as a compilation checks for code that should be 
  // developed in later releases
  // test_lcc_mixed_element_smooth();
  return EXIT_SUCCESS;
}
