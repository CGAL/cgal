#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_filtered_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>

#include <vector>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_triangulation_filtered_traits_2<K> GT;

typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT> Delaunay;
typedef Delaunay::Point                               Point;
typedef Delaunay::Face_handle                         Face_handle;

int main()
{
  Delaunay T;
  CGAL::Random_points_in_iso_rectangle_2<Point> rnd(Point(0, 0), Point(1, 1));

  // First, make sure the triangulation is 2D.
  T.insert(Point(0, 0));
  T.insert(Point(.1, 0));
  T.insert(Point(0, .1));

  // Gets the conflict region of 100 random points
  // in the Delaunay tetrahedralization
  for (int i = 0; i != 100; ++i)
    {
      Point p = (*rnd++);

      // Locate the point
      Delaunay::Locate_type lt;
      int li;
      Face_handle f = T.locate(p, lt, li);
      if (lt == Delaunay::VERTEX)
        continue; // Point already exists

      // Get the cells that conflict with p in a vector V,
      // and a facet on the boundary of this hole in f.
      std::vector<Face_handle> V;

      T.get_conflicts(p,
                      std::back_inserter(V), // Conflict cells in V
                      f);
    }

  std::cout << "Final triangulation has " << T.number_of_vertices()
            << " vertices." << std::endl;

  return 0;
}
