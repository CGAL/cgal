// file: examples/Triangulation_3/example_regular.C

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <cassert>
#include <vector>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;

typedef Traits::RT                                          Weight;
typedef Traits::Bare_point                                  Point;
typedef Traits::Weighted_point                              Weighted_point;

typedef CGAL::Regular_triangulation_3<Traits>               Rt;

typedef Rt::Vertex_iterator                                 Vertex_iterator;
typedef Rt::Vertex_handle                                   Vertex_handle;

int main()
{
  Rt T;

  // insertion of points on a 3D grid
  std::vector<Vertex_handle> V;

  for (int z=0 ; z<5 ; z++)
    for (int y=0 ; y<5 ; y++)
      for (int x=0 ; x<5 ; x++) {
	  Point p(x, y, z);
          Weight w = (x+y-z*y*x)*2.0; // let's say this is the weight.
	  Weighted_point wp(p, w);
	  V.push_back(T.insert(wp));
      }

  assert( T.is_valid() );
  assert( T.dimension() == 3 );

  std::cout << "Number of vertices : " << T.number_of_vertices() << std::endl;

  return 0;
}
