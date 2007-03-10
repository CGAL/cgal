#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>

#include <vector>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_3<K>        Delaunay;
typedef Delaunay::Point                          Point;
typedef Delaunay::Cell_handle                    Cell_handle;
typedef Delaunay::Facet                          Facet;

int main()
{
  Delaunay T;
  CGAL::Random_points_in_sphere_3<Point> rnd;

  // First, make sure the triangulation is 3D.
  T.insert(Point(0,0,0));
  T.insert(Point(1,0,0));
  T.insert(Point(0,1,0));
  T.insert(Point(0,0,1));

  assert(T.dimension() == 3);

  // Inserts 100 random points if and only if their insertion
  // in the Delaunay tetrahedralization conflicts with
  // an even number of cells.
  for (int i = 0; i != 100; ++i) {
    Point p = *rnd++;

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

    if ((V.size() & 1) == 0)  // Even number of conflict cells ?
      T.insert_in_hole(p, V.begin(), V.end(), f.first, f.second);
  }

  std::cout << "Final triangulation has " << T.number_of_vertices()
            << " vertices." << std::endl;

  return 0;
}
