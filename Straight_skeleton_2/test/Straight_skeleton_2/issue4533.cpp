#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>
#include <CGAL/Polygon_2.h>

#include <memory>

#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;

typedef K::Point_2                                                  Point;
typedef CGAL::Polygon_2<K>                                          Polygon_2;
typedef CGAL::Straight_skeleton_2<K>                                Ss;

typedef std::shared_ptr<Ss>                                       SsPtr;

int main()
{
  Polygon_2 poly;

  poly.push_back(Point(     0,     0));
  poly.push_back(Point(  2000,  8000));
  poly.push_back(Point( 10000, 10000));
  poly.push_back(Point(  2000, 12000));
  poly.push_back(Point(     0, 20000));
  poly.push_back(Point( -2000, 12000));
  poly.push_back(Point(-10000, 10000));
  poly.push_back(Point( -2000,  8000));

  assert(poly.is_simple());
  assert(poly.is_counterclockwise_oriented());

  SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

  CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
  draw(*iss);

  assert(iss->size_of_vertices() == 9);
  assert(iss->size_of_halfedges() == 32);
  assert(iss->size_of_faces() == 8);
  assert(iss->is_valid());

  return EXIT_SUCCESS;
}
