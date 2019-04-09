#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Point_3                                              Point_3;

typedef std::vector<std::size_t>                                Polygon;
typedef CGAL::Surface_mesh<Point_3>                             Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int, char**)
{
  // First, construct a polygon soup with some problems
  std::vector<Point_3> points;
  std::vector<Polygon> polygons;

  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(1,0,0));
  points.push_back(Point_3(0,1,0));
  points.push_back(Point_3(-1,0,0));
  points.push_back(Point_3(0,-1,0));
  points.push_back(Point_3(0,1,0)); // duplicate point
  points.push_back(Point_3(0,-2,0)); // unused point

  Polygon p;
  p.push_back(0); p.push_back(1); p.push_back(2);
  polygons.push_back(p);

  // degenerate face
  p.clear();
  p.push_back(0); p.push_back(0); p.push_back(0);
  polygons.push_back(p);

  p.clear();
  p.push_back(0); p.push_back(1); p.push_back(4);
  polygons.push_back(p);

  // duplicate face with different orientation
  p.clear();
  p.push_back(0); p.push_back(4); p.push_back(1);
  polygons.push_back(p);

  p.clear();
  p.push_back(0); p.push_back(3); p.push_back(5);
  polygons.push_back(p);

  // degenerate face
  p.clear();
  p.push_back(0); p.push_back(3); p.push_back(0);
  polygons.push_back(p);

  p.clear();
  p.push_back(0); p.push_back(3); p.push_back(4);
  polygons.push_back(p);

  // pinched and degenerate face
  p.clear();
  p.push_back(0); p.push_back(1); p.push_back(2); p.push_back(3);
  p.push_back(4); p.push_back(3); p.push_back(2); p.push_back(1);
  polygons.push_back(p);

  PMP::repair_polygon_soup(points, polygons);
  PMP::orient_polygon_soup(points, polygons);

  Mesh mesh;
  PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);

  std::cout << "Mesh has " << num_vertices(mesh) << " vertices and " << num_faces(mesh) << " faces" << std::endl;

  assert(num_vertices(mesh) == 5);
  assert(num_faces(mesh) == 4);

  return 0;
}
