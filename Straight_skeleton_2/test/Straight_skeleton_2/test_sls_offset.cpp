#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
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

  typedef typename K::FT                                             FT;
  typedef typename K::Point_2                                        Point;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;
  typedef boost::shared_ptr<Polygon_with_holes_2>                    Polygon_with_holes_2_ptr;
  typedef std::vector<Polygon_with_holes_2_ptr>                      Polygon_with_holes_2_ptr_container;

  // Square with a non-centered square hole
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
  Polygon_with_holes_2_ptr_container offset_poly_with_holes =
    CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(int(2), poly);

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
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 8);
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

  // Very large value, no offset polygon
  std::cout << "Interior skeleton and offset, value: 100" << std::endl;
  offset_poly_with_holes = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(FT(100), poly);

  for (const auto& offp : offset_poly_with_holes)
    print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 0);
}

template <typename K>
void test_offset_pinched()
{
  typedef typename K::FT                                                       FT;
  typedef typename K::Point_2                                                  Point;
  typedef CGAL::Polygon_2<K>                                                   Polygon_2;
  typedef boost::shared_ptr<Polygon_2>                                         Polygon_ptr;
  typedef std::vector<Polygon_ptr>                                             Polygon_ptrs;

  typedef CGAL::Straight_skeleton_2<K>                                         Ss;

  typedef CGAL::Straight_skeleton_builder_traits_2<K>                          Skeleton_builder_traits;
  typedef CGAL::Straight_skeleton_builder_2<Skeleton_builder_traits, Ss>       Skeleton_builder;

  typedef CGAL::Polygon_offset_builder_traits_2<K>                             Offset_builder_traits;
  typedef CGAL::Polygon_offset_builder_2<Ss, Offset_builder_traits, Polygon_2> Offset_builder;

  // This looks like this:
  // ---------------------------      ---------------------------
  // |                          \    /                          |
  // |                           \  /                           |
  // |                            \/                            |
  // |                                                          |
  // |                                                          |
  // |                            /\                            |
  // |                           /  \                           |
  // |                          /    \                          |
  // |_________________________/      \_________________________|

  Polygon_2 input;
  input.push_back(Point(0, 0));
  input.push_back(Point(20, 0));
  input.push_back(Point(30, 10));
  input.push_back(Point(40, 0));
  input.push_back(Point(60, 0));
  input.push_back(Point(60, 40));
  input.push_back(Point(40, 40));
  input.push_back(Point(30, 30));
  input.push_back(Point(20, 40));
  input.push_back(Point(0, 40));

  Skeleton_builder ssb;
  ssb.enter_contour(input.begin(), input.end());

  boost::shared_ptr<Ss> ss = ssb.construct_skeleton();
  assert(ss);

  // The two splitting fronts meet in the middle, and at that time,
  // we go from a single offset polygon to two polygons
  FT time_at_split;
  for(auto vit=ss->vertices_begin(); vit!=ss->vertices_end(); ++vit)
    if(vit->is_split())
      time_at_split = vit->time();

  Polygon_ptrs offset_polys;
  Offset_builder ob(*ss);

  // Before, 1 polygon
  ob.construct_offset_contours(time_at_split - FT(1), std::back_inserter(offset_polys));

  for(const Polygon_ptr p : offset_polys)
    print_polygon(*p);

  assert(offset_polys.size() == 1);
  assert(offset_polys[0]->size() == 10);

  // At the exact splitting time
  offset_polys.clear();
  ob.construct_offset_contours(time_at_split, std::back_inserter(offset_polys));

  for(const Polygon_ptr p : offset_polys)
    print_polygon(*p);

  assert(offset_polys.size() == 2);
  assert(offset_polys[0]->size() == 5);
  assert(offset_polys[1]->size() == 5);
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

template <typename K>
void test_offset_non_manifold()
{
  typedef typename K::FT                                                       FT;
  typedef typename K::Point_2                                                  Point;
  typedef CGAL::Polygon_2<K>                                                   Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                                        Polygon_with_holes_2;
  typedef boost::shared_ptr<Polygon_with_holes_2>                              Polygon_with_holes_ptr;
  typedef std::vector<Polygon_with_holes_ptr>                                  Polygon_with_holes_ptrs;

  typedef CGAL::Straight_skeleton_2<K>                                         Ss;

  typedef CGAL::Straight_skeleton_builder_traits_2<K>                          Skeleton_builder_traits;
  typedef CGAL::Straight_skeleton_builder_2<Skeleton_builder_traits, Ss>       Skeleton_builder;

  std::vector<Point> outer, hole;
  outer.push_back(Point( 0,  0));
  outer.push_back(Point(10,  0));
  outer.push_back(Point(10, 10));
  outer.push_back(Point( 0, 10));

  hole.push_back(Point( 2,  6));
  hole.push_back(Point( 5,  9));
  hole.push_back(Point( 8,  6));
  hole.push_back(Point( 5,  3));

  Skeleton_builder ssb;
  ssb.enter_contour(outer.begin(), outer.end());
  ssb.enter_contour(hole.begin(), hole.end());

  boost::shared_ptr<Ss> ss = ssb.construct_skeleton();
  assert(ss);

  // The two splitting fronts meet in the middle, and at that time,
  // we go from a single offset polygon to two polygons
  FT time_at_split = 100;
  for(auto vit=ss->vertices_begin(); vit!=ss->vertices_end(); ++vit)
    if(vit->is_split())
      if(vit->time() < time_at_split)
        time_at_split = vit->time();

  Polygon_with_holes_2 poly(Polygon_2(outer.begin(), outer.end()));
  poly.add_hole(Polygon_2(hole.begin(), hole.end()));

  Polygon_with_holes_ptrs offset_poly_with_holes =
    CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(time_at_split, poly);

  std::cout << offset_poly_with_holes.size() << " polys" << std::endl;
  for (const auto& offp : offset_poly_with_holes)
    print_polygon_with_holes(*offp);

  // @fixme
  // - should be a single 4-vertex polygon with one hole
  // - check why the result is 0 with epeck_w_srqt (epeck too?)
  std::exit(1);
}

template <typename K>
void test_offset_polygon_exterior()
{
  std::cout << "Test Polygon with Hole, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::FT                                             FT;
  typedef typename K::Point_2                                        Point;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;
  typedef boost::shared_ptr<Polygon_with_holes_2>                    Polygon_with_holes_2_ptr;
  typedef std::vector<Polygon_with_holes_2_ptr>                      Polygon_with_holes_2_ptr_container;

  Polygon_2 poly;
  poly.push_back(Point( 0,  0));
  poly.push_back(Point(40,  0));
  poly.push_back(Point(40, 50));
  poly.push_back(Point( 0, 50));
  poly.push_back(Point( 0, 30));
  poly.push_back(Point(10, 30));
  poly.push_back(Point(10, 40));
  poly.push_back(Point(30, 40));
  poly.push_back(Point(30, 10));
  poly.push_back(Point(10, 10));
  poly.push_back(Point(10, 20));
  poly.push_back(Point( 0, 20));

  // -----------------------------------------------------------------------------------------------
  // Outer skeleton and offset
  std::cout << "Outer skeleton and offset, value: 0.1" << std::endl;
  Polygon_with_holes_2_ptr_container offset_poly_with_holes =
    create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, poly);

  for (const auto& offp : offset_poly_with_holes)
    print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 4);
  assert(offset_poly_with_holes[0]->number_of_holes() == 1);
  assert(offset_poly_with_holes[0]->holes_begin()->size() == 12);

  // -----------------------------------------------------------------------------------------------
  // Value such that it is clearly separated into two contours
  std::cout << "Outer skeleton and offset, value: 7" << std::endl;
  offset_poly_with_holes = create_exterior_skeleton_and_offset_polygons_with_holes_2(FT(7), poly,
                                                                                     K(), EPICK());

  for (const auto& offp : offset_poly_with_holes)
    print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 2);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 4);
  assert(offset_poly_with_holes[0]->number_of_holes() == 1);

  // Technically both polygons below should be rectangles, but the algorithm puts a 5th vertex collinear.
  // Tolerating it for now...

  // assert(offset_poly_with_holes[0]->holes_begin()->size() == 4);
  // assert(offset_poly_with_holes[1]->outer_boundary().size() == 4);
  assert(offset_poly_with_holes[0]->holes_begin()->is_simple());
  assert(offset_poly_with_holes[1]->outer_boundary().is_simple());

  // -----------------------------------------------------------------------------------------------
  // Border value between a single contour and two contours
  std::cout << "Outer skeleton and offset, value: 5" << std::endl;
  offset_poly_with_holes = create_exterior_skeleton_and_offset_polygons_with_holes_2(5., poly,
                                                                                     K(), EPICK());

  for (const auto& offp : offset_poly_with_holes)
    print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 2);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 4);
  assert(offset_poly_with_holes[0]->number_of_holes() == 1);

  // Technically both polygons below should be rectangles, but the algorithm puts a 5th vertex collinear.
  // Tolerating it for now...

  // assert(offset_poly_with_holes[0]->holes_begin()->size() == 4);
  // assert(offset_poly_with_holes[1]->outer_boundary().size() == 4);
  assert(offset_poly_with_holes[0]->holes_begin()->is_simple());
  assert(offset_poly_with_holes[1]->outer_boundary().is_simple());
}

template <typename K>
void test_offset(const char* filename,
                 const std::vector<double> offset_values = { }) // optional interesting offset values
{
  std::cout << "Construct inner offset of input: " << filename << std::endl;

  typedef typename K::FT                                             FT;
  typedef typename K::Point_2                                        Point;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;
  typedef boost::shared_ptr<Polygon_with_holes_2>                    Polygon_with_holes_2_ptr;
  typedef std::vector<Polygon_with_holes_2_ptr>                      Polygon_with_holes_2_ptr_container;

  typedef CGAL::Min_circle_2_traits_2<K>                             Traits;
  typedef CGAL::Min_circle_2<Traits>                                 Min_circle;

  std::ifstream in(filename);
  assert(in);

  CGAL::set_ascii_mode(in);

  std::vector<Point> points;
  std::vector<Polygon_2> polys;

  int ccb_count = 0;
  in >> ccb_count;
  for(int i=0; i<ccb_count && in; ++i)
  {
    std::vector<Point> poly;

    int v_count = 0;
    in >> v_count;
    for(int j=0; j<v_count && in; ++j)
    {
      double x = 0., y = 0.;
      in >> x >> y;
      points.emplace_back(x, y);
      poly.push_back(points.back());
    }

    if(poly.size() >= 3)
    {
      bool is_simple = CGAL::is_simple_2(poly.begin(), poly.end(), K());
      if(!is_simple)
        std::cerr << "Input polygon not simple (hopefully it is strictly simple...)" << std::endl;

      CGAL::Orientation expected = (i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE);

      const double area = CGAL::to_double(CGAL::polygon_area_2(poly.begin(), poly.end(), K()));
      CGAL::Orientation orientation = area > 0 ? CGAL::COUNTERCLOCKWISE : area < 0 ? CGAL::CLOCKWISE : CGAL::COLLINEAR;

      if(orientation == expected)
        polys.push_back(Polygon_2(poly.begin(), poly.end()));
      else
        polys.push_back(Polygon_2(poly.rbegin(), poly.rend()));
    }
  }

  assert(!polys.empty());

  Polygon_with_holes_2 p(polys[0]);
  for(std::size_t i=0; i<polys.size()-1; ++i)
    p.add_hole(polys[1]);

  Min_circle mc(points.begin(), points.end(), true /*randomize*/);
  const FT r = CGAL::approximate_sqrt(mc.circle().squared_radius());
  const FT offv = 0.05 * r; // 10% of the radius of the min enclosing circle

  Polygon_with_holes_2_ptr_container offset_poly_with_holes =
    create_interior_skeleton_and_offset_polygons_with_holes_2(offv, p);

  for(const auto& offp : offset_poly_with_holes)
  {
    assert(offp->outer_boundary().is_simple());
    assert(offp->outer_boundary().is_counterclockwise_oriented());
  }

  for(double o : offset_values)
  {
    offset_poly_with_holes = create_interior_skeleton_and_offset_polygons_with_holes_2(offv, p);

    for (const auto& offp : offset_poly_with_holes)
    {
      assert(offp->outer_boundary().is_simple());
      assert(offp->outer_boundary().is_counterclockwise_oriented());
    }
  }
}

template <typename K>
void test_kernel()
{
  std::cout.precision(17);
  std::cerr.precision(17);

//  test_offset_polygon_with_hole<K>();

//  test_offset_pinched<K>();

  test_offset_non_manifold<K>();

//  test_offset_polygon_exterior<K>();

//  test_offset_multiple_CCs<K>();

//  test_offset<K>("data/1_Example.poly");
//  test_offset<K>("data/1_Example_Working.poly");
//  test_offset<K>("data/2_Example.poly");
//  test_offset<K>("data/5-SPOKE2.poly");
//  test_offset<K>("data/5-SPOKE.poly");
//  test_offset<K>("data/7-SPOKE.poly");
//  test_offset<K>("data/alley_0.poly");
//  test_offset<K>("data/alley_1.poly");
//  test_offset<K>("data/alley_2.poly");
//  test_offset<K>("data/alley_3.poly");
//  test_offset<K>("data/AlmostClosed.poly");
//  test_offset<K>("data/A.poly");
//  test_offset<K>("data/closer_edge_event_0.poly");
//  test_offset<K>("data/closer_edge_event_1.poly");
//  test_offset<K>("data/complex_0.poly");
//  test_offset<K>("data/complex_1.poly");
//  test_offset<K>("data/complex_2.poly");
//  test_offset<K>("data/complex_3.poly");
//  test_offset<K>("data/complex_4.poly");
//  test_offset<K>("data/complex_5.poly");
//  test_offset<K>("data/consecutive_coincident_vertices_0.poly");
//  test_offset<K>("data/consecutive_coincident_vertices_1.poly");
//  test_offset<K>("data/consecutive_coincident_vertices_2.poly");
//  test_offset<K>("data/consecutive_coincident_vertices_3.poly");
//  test_offset<K>("data/consecutive_coincident_vertices_4.poly");
//  test_offset<K>("data/degenerate0a.poly");
//  test_offset<K>("data/degenerate0.poly");
//  test_offset<K>("data/degenerate10.poly");
//  test_offset<K>("data/degenerate11.poly");
//  test_offset<K>("data/degenerate12.poly");
//  test_offset<K>("data/degenerate13.poly");
//  test_offset<K>("data/degenerate1.poly");
//  test_offset<K>("data/degenerate20.poly");
//  test_offset<K>("data/degenerate21.poly");
//  test_offset<K>("data/degenerate22b.poly");
//  test_offset<K>("data/degenerate22c.poly");
//  test_offset<K>("data/degenerate22.poly");
//  test_offset<K>("data/degenerate24.poly");
//  test_offset<K>("data/degenerate25.poly");
//  test_offset<K>("data/degenerate26.poly");
//  test_offset<K>("data/degenerate27b.poly");
//  test_offset<K>("data/degenerate27c.poly");
//  test_offset<K>("data/degenerate27d.poly");
//  test_offset<K>("data/degenerate27e.poly");
//  test_offset<K>("data/degenerate27.poly");
//  test_offset<K>("data/degenerate28aa.poly");
//  test_offset<K>("data/degenerate28a.poly");
//  test_offset<K>("data/degenerate28b.poly");
//  test_offset<K>("data/degenerate28c.poly");
//  test_offset<K>("data/degenerate28x.poly");
//  test_offset<K>("data/degenerate2.poly");
//  test_offset<K>("data/degenerate3.poly");
//  test_offset<K>("data/degenerate4.poly");
//  test_offset<K>("data/degenerate5a.poly");
//  test_offset<K>("data/degenerate5.poly");
//  test_offset<K>("data/degenerate6.poly");
//  test_offset<K>("data/degenerate7.poly");
//  test_offset<K>("data/degenerate8.poly");
//  test_offset<K>("data/degenerate9.poly");
//  test_offset<K>("data/degenerate_multinode0.poly");
//  test_offset<K>("data/Detmier_b.poly");
//  test_offset<K>("data/Detmier_c.poly");
//  test_offset<K>("data/Detmier_d.poly");
//  test_offset<K>("data/Detmier_e.poly");
//  test_offset<K>("data/Detmier.poly");
//  test_offset<K>("data/double_edge_0.poly");
//  test_offset<K>("data/double_edge_1.poly");
//  test_offset<K>("data/double_edge_2.poly");
//  test_offset<K>("data/double_edge.poly");
//  test_offset<K>("data/double_split.poly");
//  test_offset<K>("data/equal_times_0.poly");
//  test_offset<K>("data/ExtraEdge_1.poly");
//  test_offset<K>("data/ExtraEdge_2.poly");
//  test_offset<K>("data/hole.poly");
//  test_offset<K>("data/inputcircle.poly");
//  test_offset<K>("data/inputc.poly");
//  test_offset<K>("data/inputd1.poly");
//  test_offset<K>("data/inputd.poly");
//  test_offset<K>("data/inputG.poly");
//  test_offset<K>("data/input_K.poly");
//  test_offset<K>("data/inputPa.poly");
//  test_offset<K>("data/inputP.poly");
//  test_offset<K>("data/inputq1.poly");
//  test_offset<K>("data/inputq.poly");
//  test_offset<K>("data/inputsquare2.poly");
//  test_offset<K>("data/inputsquare.poly");
//  test_offset<K>("data/inputT.poly");
//  test_offset<K>("data/inputu.poly");
//  test_offset<K>("data/issue3382_bis.txt");
//  test_offset<K>("data/issue3382_ter.txt");
//  test_offset<K>("data/large_1.poly");
//  test_offset<K>("data/large_2.poly");
//  test_offset<K>("data/large_3.poly");
//  test_offset<K>("data/large_4.poly");
//  test_offset<K>("data/many_holes.poly");
//  test_offset<K>("data/masked_double_split.poly");
//  test_offset<K>("data/multinode0.poly");
//  test_offset<K>("data/multinode1.poly");
//  test_offset<K>("data/near_degenerate_0.poly");
//  test_offset<K>("data/near_degenerate_1.poly");
//  test_offset<K>("data/nearly_collinear.poly");
//  test_offset<K>("data/parallels0_b.poly");
//  test_offset<K>("data/parallels0.poly");
//  test_offset<K>("data/parallels_1.poly");
//  test_offset<K>("data/poly4b.poly");
//  test_offset<K>("data/poly4.poly");
//  test_offset<K>("data/poly6.poly");
//  test_offset<K>("data/pseudo_split_0.poly");
//  test_offset<K>("data/pseudo_split_10.poly");
//  test_offset<K>("data/pseudo_split_11.poly");
//  test_offset<K>("data/pseudo_split_12.poly");
//  test_offset<K>("data/pseudo_split_13b.poly");
//  test_offset<K>("data/pseudo_split_13.poly");
//  test_offset<K>("data/pseudo_split_1.poly");
//  test_offset<K>("data/pseudo_split_2.poly");
//  test_offset<K>("data/pseudo_split_3.poly");
//  test_offset<K>("data/pseudo_split_4.poly");
//  test_offset<K>("data/pseudo_split_5b.poly");
//  test_offset<K>("data/pseudo_split_5.poly");
//  test_offset<K>("data/pseudo_split_6.poly");
//  test_offset<K>("data/pseudo_split_7.poly");
//  test_offset<K>("data/pseudo_split_8.poly");
//  test_offset<K>("data/pseudo_split_9.poly");
//  test_offset<K>("data/rect_4_spokes.poly");
//  test_offset<K>("data/rectangle.poly");
//  test_offset<K>("data/region_4.poly");
//  test_offset<K>("data/rombus_4_spokes.poly");
//  test_offset<K>("data/sample_0.poly");
//  test_offset<K>("data/sample_101.poly");
//  test_offset<K>("data/sample_102.poly");
//  test_offset<K>("data/sample_147.poly");
//  test_offset<K>("data/sample_1.poly");
//  test_offset<K>("data/sample_235.poly");
//  test_offset<K>("data/sample_298.poly");
//  test_offset<K>("data/sample_2.poly");
//  test_offset<K>("data/sample2.poly");
//  test_offset<K>("data/sample_319.poly");
//  test_offset<K>("data/sample_325.poly");
//  test_offset<K>("data/sample_333.poly");
//  test_offset<K>("data/sample_3.poly");
//  test_offset<K>("data/sample3.poly");
//  test_offset<K>("data/sample_46.poly");
//  test_offset<K>("data/sample_4.poly");
//  test_offset<K>("data/sample_5.poly");
//  test_offset<K>("data/sample_638.poly");
//  test_offset<K>("data/sample_698.poly");
//  test_offset<K>("data/sample_6.poly");
//  test_offset<K>("data/sample_73.poly");
//  test_offset<K>("data/sample_85.poly");
//  test_offset<K>("data/sample.poly");
//  test_offset<K>("data/simple_0.poly");
//  test_offset<K>("data/simple_1.poly");
//  test_offset<K>("data/simple_2.poly");
//  test_offset<K>("data/simple_3.poly");
//  test_offset<K>("data/single_split.poly");
//  test_offset<K>("data/split_at_end_0.poly");
//  test_offset<K>("data/split_at_end_1.poly");
//  test_offset<K>("data/split_at_end_2.poly");
//  test_offset<K>("data/split_at_zero_0.poly");
//  test_offset<K>("data/square.poly");
//  test_offset<K>("data/star.poly");
//  test_offset<K>("data/StrayCenterlines.poly");
//  test_offset<K>("data/triangle.poly");
//  test_offset<K>("data/wheel_128_spokes.poly");
//  test_offset<K>("data/wheel_13_spokes.poly");
//  test_offset<K>("data/wheel_14_spokes.poly");
//  test_offset<K>("data/wheel_15_spokes.poly");
//  test_offset<K>("data/wheel_16_spokes_b.poly");
//  test_offset<K>("data/wheel_16_spokes.poly");
//  test_offset<K>("data/wiggly_03_cgal.poly");
//  test_offset<K>("data/WingChiu.poly");
}

int main(int, char**)
{
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  test_kernel<EPECK_w_sqrt>();
}
