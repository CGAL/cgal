#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>

#include <vector>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay;
typedef Delaunay::Point Point;

int main()
{
  // generating points on a grid.
  std::vector<Point> P;

  for (int z=0 ; z<20 ; z++)
    for (int y=0 ; y<20 ; y++)
      for (int x=0 ; x<20 ; x++)
	  P.push_back(Point(x,y,z));

  // building their Delaunay triangulation.
  Delaunay T(P.begin(), P.end());

  assert( T.number_of_vertices() == 8000 );

  // performing nearest vertex queries to a series of random points,
  // which is a case where the Fast_location policy is beneficial.
  for (int i=0; i<10000; ++i)
    T.nearest_vertex(Point(CGAL::default_random.get_double(0, 20),
			   CGAL::default_random.get_double(0, 20),
			   CGAL::default_random.get_double(0, 20)));

  return 0;
}
