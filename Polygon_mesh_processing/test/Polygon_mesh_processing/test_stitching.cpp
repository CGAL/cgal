#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <iostream>
#include <fstream>
#include <set>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

void test_polyhedron(const char* fname)
{
  typedef CGAL::Polyhedron_3<K> Polyhedron;

  std::cout << "Testing Polyhedron_3 " << fname << "..." << std::flush;
  std::ifstream input(fname);
  Polyhedron poly;
  if (!input || !(input >> poly)){
    std::cerr << "Error: can not read file.";
    return;
  }

  CGAL_assertion(poly.size_of_vertices() > 0);
  
  CGAL::Polygon_mesh_processing::stitch_borders(poly);

  CGAL_assertion(poly.is_valid(false, 5));
  std::cout << "OK\n";
}

void test_surface_mesh(const char* fname)
{
  typedef K::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> Mesh;

  std::cout << "Testing Surface_mesh " << fname << "..." << std::flush;
  std::ifstream input(fname);
  Mesh m;
  if (!input || !(input >> m)){
    std::cerr << "Error: can not read file.";
    return;
  }
  
  CGAL::Polygon_mesh_processing::stitch_borders(m);
  //todo : add a validity test

  std::cout << "OK\n";
}

int main()
{
  test_polyhedron("data_stitching/full_border.off");
  test_polyhedron("data_stitching/full_border_quads.off");
  test_polyhedron("data_stitching/half_border.off");
  test_polyhedron("data_stitching/mid_border.off");
  test_polyhedron("data_stitching/multiple_incidence.off");
  test_polyhedron("data_stitching/incidence_3.off");

  test_surface_mesh("data_stitching/full_border.off");
  test_surface_mesh("data_stitching/full_border_quads.off");
  test_surface_mesh("data_stitching/half_border.off");
  test_surface_mesh("data_stitching/mid_border.off");
  test_surface_mesh("data_stitching/multiple_incidence.off");
  test_surface_mesh("data_stitching/incidence_3.off");
}

