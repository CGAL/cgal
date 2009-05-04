#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_filtered_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>

#include <vector>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_3_triangulation_filtered_traits_3<K> GT;

typedef CGAL::Periodic_3_Delaunay_triangulation_3<GT> Delaunay;
typedef Delaunay::Point                               Point;
typedef Delaunay::Cell_handle                         Cell_handle;
typedef Delaunay::Facet                               Facet;

int main()
{
  Delaunay T;
  CGAL::Random_points_in_cube_3<Point> rnd(0.5);
  GT::Vector_3 v(0.5,0.5,0.5);

  // First, make sure the triangulation is 3D.
  T.insert(Point(0,0,0));
  T.insert(Point(.1,0,0));
  T.insert(Point(0,.1,0));
  T.insert(Point(0,0,.1));

  // Gets the conflict region of 100 random points
  // in the Delaunay tetrahedralization
  for (int i = 0; i != 100; ++i) {
    Point p = (*rnd++)+v;

    // Locate the point
    Delaunay::Locate_type lt;
    int li, lj;
    Cell_handle c = T.locate(p, lt, li, lj);
    if (lt == Delaunay::VERTEX)
      continue; // Point already exists

    // Get the cells that conflict with p in a vector V,
    // and a facet on the boundary of this hole in f.
    std::vector<Cell_handle> V;
    Facet f;

    T.find_conflicts(p, c,
		     CGAL::Oneset_iterator<Facet>(f), // Get one boundary facet
		     std::back_inserter(V));          // Conflict cells in V
    }

  std::cout << "Final triangulation has " << T.number_of_vertices()
            << " vertices." << std::endl;

  return 0;
}
