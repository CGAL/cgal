#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>

#include <iostream>
#include <fstream>
#include <vector>

int main()
{
#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

  // Delaunay T3
  typedef CGAL::Triangulation_data_structure_3<
            CGAL::Triangulation_vertex_base_3<K>,
            CGAL::Delaunay_triangulation_cell_base_3<K>,
            CGAL::Parallel_tag>                               Tds;
  typedef CGAL::Delaunay_triangulation_3<K, Tds>              Triangulation;

  typedef Triangulation::Point                                Point;

  const int NUM_INSERTED_POINTS = 5000;

  CGAL::Random_points_in_cube_3<Point> rnd(1.);

  // Construction from a vector of 1,000,000 points
  std::vector<Point> V;
  V.reserve(NUM_INSERTED_POINTS);
  for (int i = 0; i != NUM_INSERTED_POINTS; ++i)
    V.push_back(*rnd++);

  // Construct the locking data-structure, using the bounding-box of the points
  Triangulation::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1., -1., -1., 1., 1., 1.), 50);

  // Construct the triangulation in parallel
  Triangulation T(V.begin(), V.end(), &locking_ds);
  assert(T.is_valid());

#endif //CGAL_LINKED_WITH_TBB

  return 0;
}
