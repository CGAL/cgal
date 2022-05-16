#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/generators.h>

#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef K::Point_3 Point_3;

int main()
{
  Surface_mesh tet, tri;

  CGAL::make_tetrahedron(Point_3(0,0,0), Point_3(10,0,0), Point_3(0, 10,0), Point_3(0,0,10), tet);

  CGAL::make_triangle(Point_3(-1, -1, 1), Point_3(11, -1, 1), Point_3(-1, 11,1), tri);

  PMP::clip(tet, tri);
  std::cout << tet << std::endl;

  return EXIT_SUCCESS;
}
