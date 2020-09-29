#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#define CGAL_P2T2_USE_COMPACT_OFFSET

#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/point_generators_2.h>

#include <cassert>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;

typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT>       Delaunay;
typedef Delaunay::Point                                     Point;
typedef Delaunay::Face_handle                               Face_handle;

typedef K::Iso_rectangle_2                                  Iso_rectangle_2;

int main(int, char**)
{
  Delaunay T(Iso_rectangle_2(-0.5, -0.5, 0.5, 0.5));
  CGAL::Random_points_in_iso_rectangle_2<Point> rnd(Point(0, 0), Point(1, 1));

  std::cout << "1" << std::endl;
  T.insert(Point(0, 0));

  std::cout << "2" << std::endl;
  T.insert(Point(.1, .2));

  std::cout << "3" << std::endl;
  T.insert(Point(.3, .1));

  std::cout << "number of vertices: " << T.tds().number_of_vertices() << std::endl;
  std::cout << "number of faces: " << T.tds().number_of_faces() << std::endl;

  T.draw_to_OFF("P2T2.off");

  CGAL_assertion(T.is_valid());
  std::cin.get();

  // Gets the conflict region of 100 random points
  // in the Delaunay tetrahedralization
  for(int i = 0; i != 100; ++i)
  {
    Point p = (*rnd++);

    // Locate the point
    Delaunay::Locate_type lt;
    int li;
    Face_handle f = T.locate(p, lt, li);
    if(lt == Delaunay::VERTEX)
      continue; // Point already exists

    // Get the faces that conflict with p in a vector V,
    // and a facet on the boundary of this hole in f.
    std::vector<Face_handle> faces_in_conflict;
    T.find_conflicts(p, f, std::back_inserter(faces_in_conflict));
    std::cout << faces_in_conflict.size() << " faces in conflict" << std::endl;
  }

  std::cout << "Final triangulation has " << T.number_of_vertices() << " vertices." << std::endl;

  return 0;
}
