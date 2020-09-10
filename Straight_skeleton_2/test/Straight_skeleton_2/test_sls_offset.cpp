#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Polygon_offset_builder_2.h>
#include "print.h"

#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_polygon_2.h>

#include <boost/shared_ptr.hpp>

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel            EPECK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  EPECK_w_sqrt;

template <typename K>
void test_offset_polygon_with_hole()
{
  std::cout << "Test Polygon with Hole, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::Point_2                                        Point;
  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;
  typedef boost::shared_ptr<Polygon_with_holes_2>                    Polygon_with_holes_2_ptr;
  typedef std::vector<Polygon_with_holes_2_ptr>                      Polygon_with_holes_2_ptr_container;

  Polygon_2 outer, hole;
  outer.push_back(Point(-45, -45));
  outer.push_back(Point(  5, -45));
  outer.push_back(Point(  5,   5));
  outer.push_back(Point(-45,   5));

  hole.push_back(Point( -5, -30));
  hole.push_back(Point(  0, -30));
  hole.push_back(Point(  0, -35));
  hole.push_back(Point( -5, -35));

  Polygon_with_holes_2 poly(outer);
  poly.add_hole(hole);

  // Generic case, offset of the outer boundary and of the hole are distinct
  std::cout << "Interior skeleton and offset, value: 2" << std::endl;
  Polygon_with_holes_2_ptr_container offset_poly_with_holes
    = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(2, poly);

  for (const auto& offp : offset_poly_with_holes)
    print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 4);
  assert(offset_poly_with_holes[0]->number_of_holes() == 1);
  assert(offset_poly_with_holes[0]->holes_begin()->size() == 4);

  // The offset value is such that
  // - skeleton bisectors and offset edge overlap
  // - switch from 2 CCs to 1 CC
  std::cout << "Interior skeleton and offset, value: 2.5" << std::endl;
  offset_poly_with_holes = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(2.5, poly);

  for (const auto& offp : offset_poly_with_holes)
    print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 4);
  assert(offset_poly_with_holes[0]->number_of_holes() == 0);

  // Same, but using a different skeleton kernel
  std::cout << "Interior skeleton and offset, value: 2.5 (EPICK)" << std::endl;
  offset_poly_with_holes = create_interior_skeleton_and_offset_polygons_with_holes_2(2.5, poly,
                                                                                     K(), EPICK());

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 4);
  assert(offset_poly_with_holes[0]->number_of_holes() == 0);

  // Single CC in the offset
  std::cout << "Interior skeleton and offset, value: 3.5" << std::endl;
  offset_poly_with_holes = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(3.5, poly);

  for (const auto& offp : offset_poly_with_holes)
    print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 8);
  assert(offset_poly_with_holes[0]->number_of_holes() == 0);

  // Outer skeleton and offset
  std::cout << "Outer skeleton and offset, value: 0.1" << std::endl;
  offset_poly_with_holes = create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, outer);

  // Same, but different skeleton kernel
  std::cout << "Outer skeleton and offset, value: 0.1" << std::endl;
  offset_poly_with_holes = create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, outer,
                                                                                     K(), EPICK());
}

template <typename K>
void test_offset_multiple_CCs()
{
  typedef typename K::FT                                                       FT;
  typedef typename K::Point_2                                                  Point_2;
  typedef CGAL::Polygon_2<K>                                                   Contour;
  typedef boost::shared_ptr<Contour>                                           ContourPtr;
  typedef std::vector<ContourPtr>                                              Contour_sequence;

  typedef CGAL::Straight_skeleton_2<K>                                         Ss;

  typedef CGAL::Straight_skeleton_builder_traits_2<K>                          Skeleton_builder_traits;
  typedef CGAL::Straight_skeleton_builder_2<Skeleton_builder_traits, Ss>       Skeleton_builder;

  typedef CGAL::Polygon_offset_builder_traits_2<K>                             Offset_builder_traits;
  typedef CGAL::Polygon_offset_builder_2<Ss, Offset_builder_traits, Contour>   Offset_builder;

  Point_2 pts[] = { Point_2(  0,   0),
                    Point_2(700,   0),
                    Point_2(700, 600),
                    Point_2(  0, 600),
                    Point_2(  0, 300),
                    Point_2(300, 300),
                    Point_2(300, 400),
                    Point_2(600, 400),
                    Point_2(600, 100),
                    Point_2(300, 100),
                    Point_2(300, 200),
                    Point_2(  0, 200)
                  };
  std::vector<Point_2> input(pts, pts+12);

  FT offset = 50;
  boost::optional<FT> margin = CGAL::compute_outer_frame_margin(input.begin(), input.end(), offset);
  assert(margin);

  // Get the bbox of the polygon
  CGAL::Bbox_2 bbox = CGAL::bbox_2(input.begin(), input.end());

  // Compute the boundaries of the frame
  double m = CGAL::to_double(*margin);
  double fxmin = bbox.xmin() - m;
  double fxmax = bbox.xmax() + m;
  double fymin = bbox.ymin() - m;
  double fymax = bbox.ymax() + m;

  // Create the rectangular frame
  Point_2 frame[4] = { Point_2(fxmin, fymin),
                       Point_2(fxmax, fymin),
                       Point_2(fxmax, fymax),
                       Point_2(fxmin, fymax)
                     };

  Skeleton_builder ssb;

  ssb.enter_contour(frame, frame+4);
  ssb.enter_contour(input.rbegin(), input.rend());

  boost::shared_ptr<Ss> ss = ssb.construct_skeleton();
  assert(ss);

  Contour_sequence offset_contours;
  Offset_builder ob(*ss);
  ob.construct_offset_contours(offset, std::back_inserter(offset_contours));

  for(const ContourPtr p : offset_contours)
    print_polygon(*p);

  assert(offset_contours.size() == 3);
}

int main(int, char**)
{
  test_offset_polygon_with_hole<EPICK>();
  test_offset_polygon_with_hole<EPECK>();
  test_offset_polygon_with_hole<EPECK_w_sqrt>();

  test_offset_multiple_CCs<EPICK>();
  test_offset_multiple_CCs<EPECK>();
  test_offset_multiple_CCs<EPECK_w_sqrt>();
}
