#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel     EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel       EPECK;

void test_stitch_boundary_cycles(const char* fname)
{
  typedef CGAL::Surface_mesh<EPICK::Point_3>                    Mesh;

  std::ifstream input(fname);
  Mesh mesh;
  if (!input || !(input >> mesh)){
    std::cerr << "Error: can not read file.";
    return;
  }

  std::size_t res = PMP::internal::stitch_boundary_cycles(mesh);
  assert(res > 0);
  assert(is_valid(mesh));
}

template <typename K>
void test_polyhedron(const char* fname)
{
  std::cout << "Testing Polyhedron_3 " << fname << "..." << std::flush;

  typedef CGAL::Polyhedron_3<K>                                 Polyhedron;

  std::ifstream input(fname);
  Polyhedron poly;
  if (!input || !(input >> poly)){
    std::cerr << "Error: can not read file.";
    return;
  }

  assert(poly.size_of_vertices() > 0);

  PMP::stitch_borders(poly);
  poly.normalize_border();

  assert(poly.is_valid(false, 5));
  std::cout << "OK\n";
}

template <typename K>
void test_polyhedron_cc(const char* fname)
{
  std::cout << "Testing Polyhedron_3 " << fname << "..." << std::flush;

  typedef CGAL::Polyhedron_3<K>                                 Polyhedron;

  std::ifstream input(fname);
  Polyhedron poly;
  if (!input || !(input >> poly)){
    std::cerr << "Error: can not read file.";
    return;
  }

  assert(poly.size_of_vertices() > 0);

  PMP::stitch_borders(poly, params::apply_per_connected_component(true));
  poly.normalize_border();

  assert(poly.is_valid(false, 5));
  std::cout << "OK\n";
}

void test_surface_mesh(const char* fname)
{
  std::cout << "Testing Surface_mesh " << fname << "..." << std::flush;

  typedef CGAL::Surface_mesh<EPICK::Point_3>                    Mesh;

  std::ifstream input(fname);
  Mesh mesh;
  if (!input || !(input >> mesh)){
    std::cerr << "Error: can not read file.";
    return;
  }

  PMP::stitch_borders(mesh);
  assert(is_valid_polygon_mesh(mesh));

  std::cout << "OK\n";
}

void test_surface_mesh_cc(const char* fname)
{
  std::cout << "Testing Surface_mesh " << fname << "..." << std::flush;

  typedef CGAL::Surface_mesh<EPICK::Point_3>                    Mesh;

  std::ifstream input(fname);
  Mesh mesh;
  if (!input || !(input >> mesh)){
    std::cerr << "Error: can not read file.";
    return;
  }
  
  PMP::stitch_borders(mesh, params::apply_per_connected_component(true));
  assert(is_valid(mesh));

  std::cout << "OK\n";
}

void bug_test()
{
  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_3 Point_3;
  CGAL::Surface_mesh<Point_3> tm;

  CGAL::make_triangle(Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0), tm);
  CGAL::make_triangle(Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0), tm);
  CGAL::make_triangle(Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0), tm);
  CGAL::make_triangle(Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0), tm);

  CGAL::Polygon_mesh_processing::stitch_borders(tm);
}

int main()
{
  test_stitch_boundary_cycles("data_stitching/boundary_cycle.off");
  test_stitch_boundary_cycles("data_stitching/boundary_cycle_2.off");

  test_polyhedron<EPECK>("data_stitching/full_border.off");
  test_polyhedron<EPICK>("data_stitching/full_border.off");
  test_polyhedron<EPICK>("data_stitching/full_border_quads.off");
  test_polyhedron<EPICK>("data_stitching/half_border.off");
  test_polyhedron<EPICK>("data_stitching/mid_border.off");
  test_polyhedron<EPICK>("data_stitching/multiple_incidence.off");
  test_polyhedron<EPICK>("data_stitching/incidence_3.off");
  test_polyhedron<EPICK>("data_stitching/incoherent_patch_orientation.off");
  test_polyhedron<EPICK>("data_stitching/non_stitchable.off");
  test_polyhedron<EPICK>("data_stitching/deg_border.off");
  test_polyhedron<EPICK>("data_stitching/two_patches.off");
  test_polyhedron<EPICK>("data_stitching/non_manifold.off");
  test_polyhedron<EPICK>("data_stitching/non_manifold2.off");
  test_polyhedron_cc<EPICK>("data_stitching/nm_cubes.off");

  test_surface_mesh("data_stitching/full_border.off");
  test_surface_mesh("data_stitching/full_border_quads.off");
  test_surface_mesh("data_stitching/half_border.off");
  test_surface_mesh("data_stitching/mid_border.off");
  test_surface_mesh("data_stitching/multiple_incidence.off");
  test_surface_mesh("data_stitching/incidence_3.off");
  test_surface_mesh("data_stitching/incoherent_patch_orientation.off");
  test_surface_mesh("data_stitching/non_stitchable.off");
  test_surface_mesh("data_stitching/deg_border.off");
  test_surface_mesh("data_stitching/non_manifold.off");
  test_surface_mesh_cc("data_stitching/nm_cubes.off");

  bug_test();

  return EXIT_SUCCESS;
}
