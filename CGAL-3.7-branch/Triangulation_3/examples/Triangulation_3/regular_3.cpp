#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K>  Traits;

typedef Traits::RT                                          Weight;
typedef Traits::Bare_point                                  Point;
typedef Traits::Weighted_point                              Weighted_point;

typedef CGAL::Regular_triangulation_3<Traits>               Rt;

typedef Rt::Vertex_iterator                                 Vertex_iterator;
typedef Rt::Vertex_handle                                   Vertex_handle;

int main()
{
  // generate points on a 3D grid
  std::vector<Weighted_point> P;

  int number_of_points = 0;

  for (int z=0 ; z<5 ; z++)
    for (int y=0 ; y<5 ; y++)
      for (int x=0 ; x<5 ; x++) {
	  Point p(x, y, z);
          Weight w = (x+y-z*y*x)*2.0; // let's say this is the weight.
	  P.push_back(Weighted_point(p, w));
          ++number_of_points;
      }

  Rt T;

  // insert all points in a row (this is faster than one insert() at a time).
  T.insert (P.begin(), P.end());

  assert( T.is_valid() );
  assert( T.dimension() == 3 );

  std::cout << "Number of vertices : " << T.number_of_vertices() << std::endl;

  // removal of all vertices
  int count = 0;
  while (T.number_of_vertices() > 0) {
      T.remove (T.finite_vertices_begin());
      ++count;
  }

  assert( count == number_of_points );

  return 0;
}
