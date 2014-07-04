#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <fstream>
#include <vector>

int main()
{
#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
  typedef Traits::Bare_point                                  Point;

  // Regular T3
  typedef CGAL::Triangulation_data_structure_3< 
    CGAL::Triangulation_vertex_base_3<Traits>, 
    CGAL::Regular_triangulation_cell_base_3<Traits>, 
    CGAL::Parallel_tag>                                       Tds;

  typedef CGAL::Regular_triangulation_3<Traits, Tds>          Rt;
  typedef Rt::Vertex_handle                                   Vertex_handle;

  const int NUM_INSERTED_POINTS = 5000;

  CGAL::Random_points_in_cube_3<Point> rnd(1.);

  // Construction from a vector of 1,000,000 points
  std::vector<Point> V;
  V.reserve(NUM_INSERTED_POINTS);
  for (int i = 0; i != NUM_INSERTED_POINTS; ++i)
    V.push_back(*rnd++);
  
  // Construct the locking data-structure, using the bounding-box of the points
  Rt::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1., -1., -1., 1., 1., 1.), 50);
  // Contruct the triangulation in parallel
  std::cerr << "Construction and insertion" << std::endl;
  Rt rtr(V.begin(), V.end(), &locking_ds);

  assert(rtr.is_valid());

  std::cerr << "Remove" << std::endl;
  // Remove the first 1/10 vertices
  std::vector<Vertex_handle> vertices_to_remove;
  Rt::Finite_vertices_iterator vit = rtr.finite_vertices_begin();
  for (int i = 0 ; i < NUM_INSERTED_POINTS/10 ; ++i)
    vertices_to_remove.push_back(vit++);
  // Parallel remove
  rtr.remove(vertices_to_remove.begin(), vertices_to_remove.end());
  
  assert(rtr.is_valid());

#endif //CGAL_LINKED_WITH_TBB

  return 0;
}
