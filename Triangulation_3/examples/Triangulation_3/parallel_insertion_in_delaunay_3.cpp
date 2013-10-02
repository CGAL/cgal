#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// Delaunay T3
typedef CGAL::Triangulation_data_structure_3< 
  CGAL::Triangulation_vertex_base_3<K>, 
  CGAL::Triangulation_cell_base_3<K>, 
  CGAL::Parallel_tag>                          Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Triangulation;

typedef Triangulation::Point          Point;

int main()
{
  CGAL::Random_points_in_cube_3<Point> rnd(1.);

  // Construction from a vector of 1000000 points
  std::vector<Point> V;
  V.reserve(1000000);
  for (int i = 0; i != 1000000; ++i)
    V.push_back(*rnd++);
  
  // Construct the locking data-structure, using the bounding-box of the points
  Triangulation::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1., -1., -1., 1., 1., 1.), 50);
  // Contruct the triangulation in parallel
  Triangulation T(V.begin(), V.end(), &locking_ds);

  return 0;
}
