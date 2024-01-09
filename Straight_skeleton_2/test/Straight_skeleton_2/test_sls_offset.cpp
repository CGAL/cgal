#define CGAL_SLS_TEST_SPEED_THINGS_UP_FOR_THE_TESTSUITE
#define CGAL_ENABLE_DISABLE_ASSERTIONS_AT_RUNTIME

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/draw_straight_skeleton_2.h>

#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Polygon_offset_builder_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>

#include <memory>

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel            EPECK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  EPECK_w_sqrt;

template <typename K>
void test_API()
{
  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;

  typedef CGAL::Polygon_2<EPICK>                                     Polygon_2_EPICK;
  typedef CGAL::Polygon_with_holes_2<EPICK>                          Polygon_with_holes_2_EPICK;

  Polygon_2 p;
  Polygon_with_holes_2 pwh;

  std::vector< std::shared_ptr<Polygon_2> > res;
  std::vector< std::shared_ptr<Polygon_2_EPICK> > res_EPICK;
  std::vector< std::shared_ptr<Polygon_with_holes_2> > res_w;
  std::vector< std::shared_ptr<Polygon_with_holes_2_EPICK> > res_w_EPICK;

  // First kernel is the offset construction (and thus output kernel), second kernel is the skeleton construction

  // simple interior, no holes
  res_EPICK = create_interior_skeleton_and_offset_polygons_2(0.1, p) ;
  res_EPICK = create_interior_skeleton_and_offset_polygons_2(0.1, p, EPICK()) ;
  res_EPICK = create_interior_skeleton_and_offset_polygons_2(0, p, EPICK(), K()) ;
  res = create_interior_skeleton_and_offset_polygons_2(0.1, p, K()) ;
  res = create_interior_skeleton_and_offset_polygons_2(0, p, K(), EPICK()) ;
  res = create_interior_skeleton_and_offset_polygons_2(0.1, p, K(), K()) ;

  // simple interior, holes
  res_EPICK = create_interior_skeleton_and_offset_polygons_2(0.1, pwh) ;
  res_EPICK = create_interior_skeleton_and_offset_polygons_2(0.1, pwh, EPICK()) ;
  res_EPICK = create_interior_skeleton_and_offset_polygons_2(0, pwh, EPICK(), K()) ;
  res = create_interior_skeleton_and_offset_polygons_2(0.1, pwh, K()) ;
  res = create_interior_skeleton_and_offset_polygons_2(0, pwh, K(), EPICK()) ;
  res = create_interior_skeleton_and_offset_polygons_2(0.1, pwh, K(), K()) ;

  // simple exterior, no holes
  res_EPICK = create_exterior_skeleton_and_offset_polygons_2(0.1, p) ;
  res_EPICK = create_exterior_skeleton_and_offset_polygons_2(0.1, p, EPICK()) ;
  res_EPICK = create_exterior_skeleton_and_offset_polygons_2(0, p, EPICK(), K()) ;
  res = create_exterior_skeleton_and_offset_polygons_2(0.1, p, K()) ;
  res = create_exterior_skeleton_and_offset_polygons_2(0, p, K(), EPICK()) ;
  res = create_exterior_skeleton_and_offset_polygons_2(0.1, p, K(), K()) ;

  // simple exterior, holes
  res_EPICK = create_exterior_skeleton_and_offset_polygons_2(0.1, pwh) ;
  res_EPICK = create_exterior_skeleton_and_offset_polygons_2(0.1, pwh, EPICK()) ;
  res_EPICK = create_exterior_skeleton_and_offset_polygons_2(0, pwh, EPICK(), K()) ;
  res = create_exterior_skeleton_and_offset_polygons_2(0.1, pwh, K()) ;
  res = create_exterior_skeleton_and_offset_polygons_2(0, pwh, K(), EPICK()) ;
  res = create_exterior_skeleton_and_offset_polygons_2(0.1, pwh, K(), K()) ;

  // Same, but the result has holes --------------------

  // arranged interior, no holes
  res_w_EPICK = create_interior_skeleton_and_offset_polygons_with_holes_2(0.1, p) ;
  res_w_EPICK = create_interior_skeleton_and_offset_polygons_with_holes_2(0.1, p, EPICK()) ;
  res_w_EPICK = create_interior_skeleton_and_offset_polygons_with_holes_2(0, p, EPICK(), K()) ;
  res_w = create_interior_skeleton_and_offset_polygons_with_holes_2(0.1, p, K()) ;
  res_w = create_interior_skeleton_and_offset_polygons_with_holes_2(0, p, K(), EPICK()) ;
  res_w = create_interior_skeleton_and_offset_polygons_with_holes_2(0.1, p, K(), K()) ;

  // arranged interior, holes
  res_w_EPICK = create_interior_skeleton_and_offset_polygons_with_holes_2(0.1, pwh) ;
  res_w_EPICK = create_interior_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, EPICK()) ;
  res_w_EPICK = create_interior_skeleton_and_offset_polygons_with_holes_2(0, pwh, EPICK(), K()) ;
  res_w = create_interior_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, K()) ;
  res_w = create_interior_skeleton_and_offset_polygons_with_holes_2(0, pwh, K(), EPICK()) ;
  res_w = create_interior_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, K(), K()) ;

  // arranged exterior, no holes
  res_w_EPICK = create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, p) ;
  res_w_EPICK = create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, p, EPICK()) ;
  res_w_EPICK = create_exterior_skeleton_and_offset_polygons_with_holes_2(0, p, EPICK(), K()) ;
  res_w = create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, p, K()) ;
  res_w = create_exterior_skeleton_and_offset_polygons_with_holes_2(0, p, K(), EPICK()) ;
  res_w = create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, p, K(), K()) ;

  // arranged exterior, holes
  res_w_EPICK = create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, pwh) ;
  res_w_EPICK = create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, EPICK()) ;
  res_w_EPICK = create_exterior_skeleton_and_offset_polygons_with_holes_2(0, pwh, EPICK(), K()) ;
  res_w = create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, K()) ;
  res_w = create_exterior_skeleton_and_offset_polygons_with_holes_2(0, pwh, K(), EPICK()) ;
  res_w = create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, K(), K()) ;
}

template <typename K, typename StraightSkeleton>
bool is_valid(const std::shared_ptr<StraightSkeleton>& ss)
{
  typedef typename StraightSkeleton::Traits::Point_2 Point;

  if(!ss)
    return false;

  if(!ss->is_valid())
    return false;

  assert(static_cast<std::size_t>(std::distance(ss->vertices_begin(), ss->vertices_end())) ==
           static_cast<std::size_t>(ss->size_of_vertices()));

  std::set<Point> unique_vertices;
  for(auto vit=ss->vertices_begin(); vit!=ss->vertices_end(); ++vit)
    unique_vertices.insert(vit->point());

  std::cout << unique_vertices.size() << " unique vertices (" << ss->size_of_vertices() << ")" << std::endl;

  // Can't guarantee that the embedding is correct with EPICK
  if(!std::is_same<K, EPICK>::value)
  {
    if(unique_vertices.size() != ss->size_of_vertices())
      return false;

    for(auto hit=ss->halfedges_begin(); hit!=ss->halfedges_end(); ++hit)
      if(hit->vertex()->point() == hit->opposite()->vertex()->point())
        return false;
  }

  return true;
}

// Test what happens when the offset is a single point
template <typename K>
void test_offset_square()
{
  std::cout << " --- Test Square, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::Point_2                                                  Point;
  typedef CGAL::Polygon_2<K>                                                   Polygon_2;
  typedef std::shared_ptr<Polygon_2>                                         Polygon_ptr;
  typedef std::vector<Polygon_ptr>                                             Polygon_ptr_container;

  typedef CGAL::Straight_skeleton_2<K>                                         Ss;

  typedef CGAL::Straight_skeleton_builder_traits_2<K>                          Skeleton_builder_traits;
  typedef CGAL::Straight_skeleton_builder_2<Skeleton_builder_traits, Ss>       Skeleton_builder;

  std::vector<Point> square;
  square.push_back(Point(0, 0));
  square.push_back(Point(1, 0));
  square.push_back(Point(1, 1));
  square.push_back(Point(0, 1));

  Skeleton_builder ssb;
  ssb.enter_contour(square.begin(), square.end());

  std::shared_ptr<Ss> ss = ssb.construct_skeleton();
  assert(is_valid<K>(ss));

  Polygon_ptr_container offset_polys =
    CGAL::create_interior_skeleton_and_offset_polygons_2(0.5, Polygon_2(square.begin(), square.end()), K());

  std::cout << offset_polys.size() << " polygons" << std::endl;
  for(const auto& offp : offset_polys)
    CGAL::Straight_skeletons_2::IO::print_polygon(*offp);

  assert(offset_polys.empty());
}

// Test what happens when the offset is a single point
template <typename K>
void test_offset_four_square_holes()
{
  std::cout << " --- Test square with four holes, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::Point_2                                        Point;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;
  typedef std::shared_ptr<Polygon_with_holes_2>                    Polygon_with_holes_2_ptr;
  typedef std::vector<Polygon_with_holes_2_ptr>                      Polygon_with_holes_2_ptr_container;

  Polygon_2 outer, hole1, hole2, hole3, hole4;
  outer.push_back(Point( 0,  0));
  outer.push_back(Point(10,  0));
  outer.push_back(Point(10, 10));
  outer.push_back(Point(0,  10));

  hole1.push_back(Point(1, 1));
  hole1.push_back(Point(1, 4.5));
  hole1.push_back(Point(4.5, 4.5));
  hole1.push_back(Point(4.5, 1));

  hole2.push_back(Point(5.5, 1));
  hole2.push_back(Point(5.5, 4.5));
  hole2.push_back(Point(9, 4.5));
  hole2.push_back(Point(9, 1));

  hole3.push_back(Point(1, 5.5));
  hole3.push_back(Point(1, 9));
  hole3.push_back(Point(4.5, 9));
  hole3.push_back(Point(4.5, 5.5));

  hole4.push_back(Point(5.5, 5.5));
  hole4.push_back(Point(5.5, 9));
  hole4.push_back(Point(9, 9));
  hole4.push_back(Point(9, 5.5));

  Polygon_with_holes_2 poly(outer);
  poly.add_hole(hole1);
  poly.add_hole(hole2);
  poly.add_hole(hole3);
  poly.add_hole(hole4);

  // Generic case, offset of the outer boundary and of the hole are distinct
  std::cout << "Interior skeleton and offset, value: 2" << std::endl;
  Polygon_with_holes_2_ptr_container offset_poly_with_holes =
    CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(0.5, poly, K());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  // @todo
//  assert(offset_poly_with_holes.size() == 1);
//  assert(offset_poly_with_holes[0]->outer_boundary().size() == 4);
//  assert(offset_poly_with_holes[0]->number_of_holes() == 1);
//  assert(offset_poly_with_holes[0]->holes_begin()->size() == 4);
}

// Test what happens when part an offset polygon collapses
template <typename K>
void test_offset_L()
{
  std::cout << " --- Test L, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::Point_2                                                  Point;
  typedef CGAL::Polygon_2<K>                                                   Polygon_2;
  typedef std::shared_ptr<Polygon_2>                                         Polygon_ptr;
  typedef std::vector<Polygon_ptr>                                             Polygon_ptr_container;

  typedef CGAL::Straight_skeleton_2<K>                                         Ss;

  typedef CGAL::Straight_skeleton_builder_traits_2<K>                          Skeleton_builder_traits;
  typedef CGAL::Straight_skeleton_builder_2<Skeleton_builder_traits, Ss>       Skeleton_builder;

  std::vector<Point> L;
  L.push_back(Point( 0  ,  0));
  L.push_back(Point( 0.2,  0));
  L.push_back(Point(10  ,  0));
  L.push_back(Point(10  , 10));
  L.push_back(Point( 4  , 10));
  L.push_back(Point( 4  ,  2));
  L.push_back(Point( 0.2,  2));
  L.push_back(Point( 0  ,  2));

  Skeleton_builder ssb;
  ssb.enter_contour(L.begin(), L.end());

  std::shared_ptr<Ss> ss = ssb.construct_skeleton();
  assert(is_valid<K>(ss));

  Polygon_ptr_container offset_polys =
    CGAL::create_interior_skeleton_and_offset_polygons_2(int(1), Polygon_2(L.begin(), L.end()), K());

  std::cout << offset_polys.size() << " polygons" << std::endl;
//  for(const auto& offp : offset_polys)
//    CGAL::Straight_skeletons_2::IO::print_polygon(*offp);

  assert(offset_polys.size() == 1);
  assert(offset_polys[0]->size() == 4);
}

template <typename K>
void test_offset_polygon_with_hole()
{
  std::cout << " --- Test Polygon with Hole, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::FT                                             FT;
  typedef typename K::Point_2                                        Point;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;
  typedef std::shared_ptr<Polygon_with_holes_2>                    Polygon_with_holes_2_ptr;
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
    CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(int(2), poly, K());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 4);
  assert(offset_poly_with_holes[0]->number_of_holes() == 1);
  assert(offset_poly_with_holes[0]->holes_begin()->size() == 4);

  // The offset value is such that
  // - skeleton bisectors and offset edge overlap
  // - switch from 2 CCs to 1 CC
  std::cout << "Interior skeleton and offset, value: 2.5" << std::endl;
  offset_poly_with_holes = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(2.5, poly, K());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 8);
  assert(offset_poly_with_holes[0]->number_of_holes() == 0);

  // Same, but using a different skeleton kernel
  std::cout << "Interior skeleton and offset, value: 2.5 (EPICK)" << std::endl;
  offset_poly_with_holes = create_interior_skeleton_and_offset_polygons_with_holes_2(2.5, poly, K(), EPICK());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 8);
  assert(offset_poly_with_holes[0]->number_of_holes() == 0);

  // Single CC in the offset
  std::cout << "Interior skeleton and offset, value: 3.5" << std::endl;
  offset_poly_with_holes = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(3.5, poly, K());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 8);
  assert(offset_poly_with_holes[0]->number_of_holes() == 0);

  // Very large value, no offset polygon
  std::cout << "Interior skeleton and offset, value: 100" << std::endl;
  offset_poly_with_holes = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(FT(100), poly, K());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 0);
}

template <typename K>
void test_offset_pinched()
{
  std::cout << " --- Test Pinched, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::FT                                                       FT;
  typedef typename K::Point_2                                                  Point;
  typedef CGAL::Polygon_2<K>                                                   Polygon_2;
  typedef std::shared_ptr<Polygon_2>                                         Polygon_ptr;
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

  std::shared_ptr<Ss> ss = ssb.construct_skeleton();
  assert(is_valid<K>(ss));

  // The two splitting fronts meet in the middle, and at that time,
  // we go from a single offset polygon to two polygons
  FT time_at_split = 1;
  for(auto vit=ss->vertices_begin(); vit!=ss->vertices_end(); ++vit)
    if(vit->is_split())
      time_at_split = vit->time();

  Polygon_ptrs offset_polys;
  Offset_builder ob(*ss);

  // Before, 1 polygon
  ob.construct_offset_contours(time_at_split - FT(1), std::back_inserter(offset_polys));

//  for(const Polygon_ptr p : offset_polys)
//    print_polygon(*p);

  assert(offset_polys.size() == 1);
  assert(offset_polys[0]->size() == 10);

  // At the exact splitting time
  offset_polys.clear();
  ob.construct_offset_contours(time_at_split, std::back_inserter(offset_polys));

//  for(const Polygon_ptr p : offset_polys)
//    print_polygon(*p);

  // This splits into two polygons, but might have just created a pinched polygon too...
//  assert(offset_polys.size() == 2);
//  assert(offset_polys[0]->size() == 5);
//  assert(offset_polys[1]->size() == 5);
}

template <typename K>
void test_offset_multiple_CCs()
{
  std::cout << " --- Test Multi CC, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::FT                                                       FT;
  typedef typename K::Point_2                                                  Point_2;
  typedef CGAL::Polygon_2<K>                                                   Contour;
  typedef std::shared_ptr<Contour>                                           ContourPtr;
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

  const FT offset = 50;
  std::optional<FT> margin = CGAL::compute_outer_frame_margin(input.begin(), input.end(), offset);
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

  std::shared_ptr<Ss> ss = ssb.construct_skeleton();
  assert(is_valid<K>(ss));

  Contour_sequence offset_contours;
  Offset_builder ob(*ss);
  ob.construct_offset_contours(offset, std::back_inserter(offset_contours));

//  for(const ContourPtr p : offset_contours)
//    print_polygon(*p);

  assert(offset_contours.size() == 3);
}

template <typename K>
void test_offset_non_manifold()
{
  std::cout << " --- Test Non Manifold #1, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::FT                                                       FT;
  typedef typename K::Point_2                                                  Point;
  typedef CGAL::Polygon_2<K>                                                   Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                                        Polygon_with_holes_2;
  typedef std::shared_ptr<Polygon_with_holes_2>                              Polygon_with_holes_ptr;
  typedef std::vector<Polygon_with_holes_ptr>                                  Polygon_with_holes_ptrs;

  typedef CGAL::Straight_skeleton_2<K>                                         Ss;

  typedef CGAL::Straight_skeleton_builder_traits_2<K>                          Skeleton_builder_traits;
  typedef CGAL::Straight_skeleton_builder_2<Skeleton_builder_traits, Ss>       Skeleton_builder;

  // Square with a diamond
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

  std::shared_ptr<Ss> ss = ssb.construct_skeleton();
  assert(is_valid<K>(ss));

  // The two splitting fronts meet in the middle, and at that time,
  // we go from a single offset polygon to two polygons
  std::set<FT> times;

  FT time_at_split = 100;
  for(auto vit=ss->vertices_begin(); vit!=ss->vertices_end(); ++vit)
  {
    times.insert(vit->time());
    if(vit->is_split())
      if(vit->time() < time_at_split)
        time_at_split = vit->time();
  }

  Polygon_with_holes_2 poly(Polygon_2(outer.begin(), outer.end()));
  poly.add_hole(Polygon_2(hole.begin(), hole.end()));

  // Basic, one polygon with a hole
  Polygon_with_holes_ptrs offset_poly_with_holes =
    CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(0.1, poly, K());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  // The way the algorithm currently works, this creates a non-simple polygon
  // and not a square offset with a diamond-hole tangent to its border
  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 4);
  assert(offset_poly_with_holes[0]->number_of_holes() == 1);
  assert(offset_poly_with_holes[0]->holes_begin()->size() == 4);

  // Outer and hole polygons merge into a single polygon
  std::cout << "time_at_split: " << time_at_split << std::endl;
  offset_poly_with_holes =
    CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(time_at_split, poly, K());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  // The way the algorithm currently works, this sometimes creates a non-simple (but strictly simple) polygon
  // and not a square offset with a diamond-hole tangent to its border
  assert(offset_poly_with_holes.size() == 1);
//  assert(offset_poly_with_holes[0]->outer_boundary().size() == 9);
//  assert(offset_poly_with_holes[0]->number_of_holes() == 0);
}

template <typename K>
void test_offset_non_manifold_2()
{
  std::cout << " --- Test Non Manifold #2, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::FT                                                       FT;
  typedef typename K::Point_2                                                  Point;
  typedef CGAL::Polygon_2<K>                                                   Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                                        Polygon_with_holes_2;
  typedef std::shared_ptr<Polygon_with_holes_2>                              Polygon_with_holes_ptr;
  typedef std::vector<Polygon_with_holes_ptr>                                  Polygon_with_holes_ptrs;

  typedef CGAL::Straight_skeleton_2<K>                                         Ss;

  typedef CGAL::Straight_skeleton_builder_traits_2<K>                          Skeleton_builder_traits;
  typedef CGAL::Straight_skeleton_builder_2<Skeleton_builder_traits, Ss>       Skeleton_builder;

  std::vector<Point> outer, hole;
  outer.push_back(Point(15, 10));
  outer.push_back(Point(20, 20));
  outer.push_back(Point(25, 10));
  outer.push_back(Point(25,  0));
  outer.push_back(Point(40,  0));
  outer.push_back(Point(40, 45));
  outer.push_back(Point( 0, 45));
  outer.push_back(Point( 0,  0));
  outer.push_back(Point(15,  0));

  hole.push_back(Point( 5,  5));
  hole.push_back(Point( 5, 40));
  hole.push_back(Point(35, 40));
  hole.push_back(Point(35,  5));
  hole.push_back(Point(30,  5));
  hole.push_back(Point(30, 25));
  hole.push_back(Point(10, 25));
  hole.push_back(Point(10,  5));

  Skeleton_builder ssb;
  ssb.enter_contour(outer.begin(), outer.end());
  ssb.enter_contour(hole.begin(), hole.end());

  std::shared_ptr<Ss> ss = ssb.construct_skeleton();
  assert(is_valid<K>(ss));

  // Similar to the previous function, a split event happens and at that particular time,
  // the offset of the hole and of the outer border merge into a single polygon
  FT time_at_split = 100;
  for(auto vit=ss->vertices_begin(); vit!=ss->vertices_end(); ++vit)
  {
    if(vit->is_split())
      if(vit->time() < time_at_split)
        time_at_split = vit->time();
  }

  Polygon_with_holes_2 poly(Polygon_2(outer.begin(), outer.end()));
  poly.add_hole(Polygon_2(hole.begin(), hole.end()));

  // Before, one polygon with a hole
  Polygon_with_holes_ptrs offset_poly_with_holes =
    CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(0.1, poly, K());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 9);
  assert(offset_poly_with_holes[0]->number_of_holes() == 1);
  assert(offset_poly_with_holes[0]->holes_begin()->size() == 8);

  // At the split
  std::cout << "time_at_split: " << time_at_split << std::endl;
  offset_poly_with_holes =
    CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(time_at_split, poly, K());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  // Current implementation will give a pinched polygon, but it could be a polygon with one incident hole
  assert(offset_poly_with_holes.size() == 1);
//  assert(offset_poly_with_holes[0]->outer_boundary().size() == 18);
//  assert(offset_poly_with_holes[0]->number_of_holes() == 0);
}

template <typename K>
void test_offset_polygon_exterior()
{
  std::cout << " --- Test Polygon exterior, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::FT                                             FT;
  typedef typename K::Point_2                                        Point;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;
  typedef std::shared_ptr<Polygon_with_holes_2>                    Polygon_with_holes_2_ptr;
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
    create_exterior_skeleton_and_offset_polygons_with_holes_2(0.1, poly, K());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->outer_boundary().size() == 12);
  assert(offset_poly_with_holes[0]->number_of_holes() == 0);

  // -----------------------------------------------------------------------------------------------
  // Value such that it is clearly separated into two contours
  std::cout << "Outer skeleton and offset, value: 7" << std::endl;
  offset_poly_with_holes =
    create_exterior_skeleton_and_offset_polygons_with_holes_2(FT(7), poly, K(), EPICK());

  // for(const auto& offp : offset_poly_with_holes)
  //   CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() == 1);
  assert(offset_poly_with_holes[0]->number_of_holes() == 1);

  // Technically both polygons below should be rectangles, but the algorithm puts a 5th vertex collinear.
  // Tolerating it for now...

  assert(offset_poly_with_holes[0]->holes_begin()->size() >= 4);
  assert(offset_poly_with_holes[0]->outer_boundary().size() >= 4);
  assert(offset_poly_with_holes[0]->holes_begin()->is_simple());
  assert(offset_poly_with_holes[0]->outer_boundary().is_simple());

  // -----------------------------------------------------------------------------------------------
  // Border value between a single contour and two contours
  std::cout << "Outer skeleton and offset, value: 5" << std::endl;
  offset_poly_with_holes =
    create_exterior_skeleton_and_offset_polygons_with_holes_2(5., poly, K(), EPICK());

  std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
//  for(const auto& offp : offset_poly_with_holes)
//    CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

  assert(offset_poly_with_holes.size() >= 1);
  assert(offset_poly_with_holes[0]->number_of_holes() == 1);

  // Technically both polygons below should be rectangles, but the algorithm puts a 5th vertex collinear.
  // Tolerating it for now...

  assert(offset_poly_with_holes[0]->holes_begin()->size() >= 4);
  assert(offset_poly_with_holes[0]->outer_boundary().size() >= 4);
  assert(offset_poly_with_holes[0]->holes_begin()->is_simple());
  assert(offset_poly_with_holes[0]->outer_boundary().is_simple());
}

template <typename K>
void test_offset_polygon_with_holes_exterior()
{
  std::cout << " --- Test Polygon exterior, kernel: " << typeid(K).name() << std::endl;

  typedef typename K::Point_2                                        Point;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;
  typedef std::shared_ptr<Polygon_with_holes_2>                    Polygon_with_holes_2_ptr;
  typedef std::vector<Polygon_with_holes_2_ptr>                      Polygon_with_holes_2_ptr_container;

  Polygon_2 outer ;
    outer.push_back( Point( 10.0, 10.0) ) ;
    outer.push_back( Point(-10.0, 10.0) ) ;
    outer.push_back( Point(-10.0, -10.0) ) ;
    outer.push_back( Point(10.0, -10.0) ) ;

  Polygon_2 hole ;
    hole.push_back( Point(5.0,5.0) ) ;
    hole.push_back( Point(5.0,-5.0) ) ;
    hole.push_back( Point(-5.0,-5.0) ) ;
    hole.push_back( Point(-5.0,5.0) ) ;

  Polygon_with_holes_2 pwh(outer) ;
  pwh.add_hole( hole ) ;

  Polygon_with_holes_2_ptr_container offset_poly_with_holes_1 =
    CGAL::create_exterior_skeleton_and_offset_polygons_with_holes_2(1., pwh, K(), EPICK());
  assert(offset_poly_with_holes_1.size()==1);
  assert(offset_poly_with_holes_1[0]->number_of_holes()==1);

  Polygon_with_holes_2_ptr_container offset_poly_with_holes_2 =
    CGAL::create_exterior_skeleton_and_offset_polygons_with_holes_2(5., pwh, K(), EPICK());
  assert(offset_poly_with_holes_2.size()==1);
  assert(offset_poly_with_holes_2[0]->number_of_holes()==0);
}

template <typename K>
void test_offset(const char* filename,
                 const bool skeleton_only = false)
{
  std::cout << "Construct inner offset of input: " << filename
            << ", Kernel = " << typeid(K).name() << std::endl;

  typedef typename K::FT                                             FT;
  typedef typename K::Point_2                                        Point;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;
  typedef std::shared_ptr<Polygon_with_holes_2>                    Polygon_with_holes_2_ptr;
  typedef std::vector<Polygon_with_holes_2_ptr>                      Polygon_with_holes_2_ptr_container;

  typedef CGAL::Straight_skeleton_2<K>                               Ss;

  typedef CGAL::Straight_skeleton_builder_traits_2<K>                Ss_builder_traits;
  typedef CGAL::Straight_skeleton_builder_2<Ss_builder_traits, Ss>   Skeleton_builder;

  typedef CGAL::Min_circle_2_traits_2<K>                             Traits;
  typedef CGAL::Min_circle_2<Traits>                                 Min_circle;

  std::ifstream in(filename);
  assert(in);

  CGAL::IO::set_ascii_mode(in);

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
  Skeleton_builder ssb;
  ssb.enter_contour(polys[0].begin(), polys[0].end());

  for(std::size_t i=0; i<polys.size()-1; ++i)
  {
    p.add_hole(polys[i+1]);
    ssb.enter_contour(polys[i+1].begin(), polys[i+1].end());
  }

  std::shared_ptr<Ss> ss = ssb.construct_skeleton();
  assert(is_valid<K>(ss));

  if(skeleton_only)
  {
    CGAL::draw(*ss);
    return;
  }

  std::set<FT> offset_times;
  for(auto vit=ss->vertices_begin(); vit!=ss->vertices_end(); ++vit)
  {
    if(vit->time() != 0)
    {
      // if the kernel is not exact, avoid trouble
      if(std::is_same<K, EPECK_w_sqrt>::value)
        offset_times.insert(vit->time());
      else
        offset_times.insert(0.99 * vit->time());
    }
  }

  Min_circle mc(points.begin(), points.end(), true /*randomize*/);
  const FT r = CGAL::approximate_sqrt(mc.circle().squared_radius());
  const FT offt = 0.05 * r; // 10% of the radius of the min enclosing circle
  offset_times.insert(offt);

  // Offset construction
  int i = 0;
  for(const FT& ot : offset_times)
  {
    std::cout << "Offset #" << i++ << " = " << ot << std::endl;
    Polygon_with_holes_2_ptr_container offset_poly_with_holes =
      CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(ot, p, K());

    std::cout << offset_poly_with_holes.size() << " polygons with holes" << std::endl;
    for(const auto& offp : offset_poly_with_holes)
      CGAL::Straight_skeletons_2::IO::print_polygon_with_holes(*offp);

    CGAL::set_use_assertions(false);
    for(const auto& offp : offset_poly_with_holes){
      (void)offp;
      assert(offp->outer_boundary().is_counterclockwise_oriented());
    }
    CGAL::set_use_assertions(true);

#ifdef CGAL_SLS_TEST_SPEED_THINGS_UP_FOR_THE_TESTSUITE
    if(i > 2)
      break;
#endif
  }
}

template <typename K>
void test_kernel()
{
  std::cout.precision(17);
  std::cerr.precision(17);

#ifndef CGAL_SLS_TEST_SPEED_THINGS_UP_FOR_THE_TESTSUITE
  // test_API<K>();
#endif

  // Artificial data
  test_offset_square<K>();
  test_offset_four_square_holes<K>();
  test_offset_L<K>();
  test_offset_polygon_with_hole<K>();
  test_offset_pinched<K>();
  test_offset_non_manifold<K>();
  test_offset_non_manifold_2<K>();
  test_offset_polygon_exterior<K>();
  test_offset_polygon_with_holes_exterior<K>();
  test_offset_multiple_CCs<K>();

  // Real data
  test_offset<K>("data/1_Example.poly");
  test_offset<K>("data/1_Example_Working.poly");
  test_offset<K>("data/2_Example.poly");
  test_offset<K>("data/5-SPOKE2.poly");
  test_offset<K>("data/5-SPOKE.poly");
  test_offset<K>("data/7-SPOKE.poly");
  test_offset<K>("data/AlmostClosed.poly");
  test_offset<K>("data/A.poly");
  test_offset<K>("data/closer_edge_event_0.poly");
  test_offset<K>("data/closer_edge_event_1.poly");
  test_offset<K>("data/consecutive_coincident_vertices_0.poly");
  test_offset<K>("data/consecutive_coincident_vertices_1.poly");
  test_offset<K>("data/consecutive_coincident_vertices_2.poly");
  test_offset<K>("data/consecutive_coincident_vertices_3.poly");
  test_offset<K>("data/consecutive_coincident_vertices_4.poly");
  test_offset<K>("data/degenerate0a.poly");
  test_offset<K>("data/degenerate0.poly");
  test_offset<K>("data/degenerate2.poly");
  test_offset<K>("data/degenerate3.poly");
  test_offset<K>("data/degenerate4.poly");
  test_offset<K>("data/degenerate5a.poly");
  test_offset<K>("data/degenerate5.poly");
  test_offset<K>("data/degenerate6.poly");
  test_offset<K>("data/degenerate7.poly");
  test_offset<K>("data/degenerate8.poly");
  test_offset<K>("data/degenerate9.poly");
  test_offset<K>("data/degenerate10.poly");
  test_offset<K>("data/degenerate11.poly");
  test_offset<K>("data/degenerate12.poly");
  test_offset<K>("data/degenerate13.poly");
  test_offset<K>("data/degenerate1.poly");
  test_offset<K>("data/degenerate20.poly");
  test_offset<K>("data/degenerate21.poly");
  test_offset<K>("data/degenerate22.poly");
  test_offset<K>("data/degenerate22b.poly");
  test_offset<K>("data/degenerate22c.poly");
  test_offset<K>("data/degenerate24.poly");
  test_offset<K>("data/degenerate25.poly");
  test_offset<K>("data/degenerate26.poly");
  test_offset<K>("data/degenerate27.poly");
  test_offset<K>("data/degenerate27b.poly");
  test_offset<K>("data/degenerate27c.poly");
  test_offset<K>("data/degenerate27d.poly");
  test_offset<K>("data/degenerate27e.poly");
  test_offset<K>("data/degenerate28a.poly");
  test_offset<K>("data/degenerate28aa.poly");
  test_offset<K>("data/degenerate28b.poly");
  test_offset<K>("data/degenerate28c.poly");
  test_offset<K>("data/degenerate28x.poly");
  test_offset<K>("data/degenerate_multinode0.poly");
  test_offset<K>("data/Detmier_b.poly");
  test_offset<K>("data/Detmier_c.poly");
  test_offset<K>("data/Detmier_d.poly");
  test_offset<K>("data/Detmier_e.poly");
  test_offset<K>("data/Detmier.poly");
  test_offset<K>("data/double_edge_0.poly");
  test_offset<K>("data/double_edge_1.poly");
  test_offset<K>("data/double_edge_2.poly");
  test_offset<K>("data/double_edge.poly");
  test_offset<K>("data/double_split.poly");
  test_offset<K>("data/equal_times_0.poly");
  test_offset<K>("data/ExtraEdge_1.poly");
  test_offset<K>("data/ExtraEdge_2.poly");
  test_offset<K>("data/hole.poly");
  test_offset<K>("data/inputcircle.poly");
  test_offset<K>("data/inputsquare2.poly");
  test_offset<K>("data/inputsquare.poly");
  test_offset<K>("data/many_holes.poly");
  test_offset<K>("data/masked_double_split.poly");
  test_offset<K>("data/multinode0.poly");
  test_offset<K>("data/multinode1.poly");
  test_offset<K>("data/near_degenerate_0.poly");
  test_offset<K>("data/near_degenerate_1.poly");
  test_offset<K>("data/nearly_collinear.poly");
  test_offset<K>("data/parallels0.poly");
  test_offset<K>("data/parallels0_b.poly");
  test_offset<K>("data/parallels_1.poly");
  test_offset<K>("data/poly4b.poly");
  test_offset<K>("data/poly4.poly");
  test_offset<K>("data/poly6.poly");
  test_offset<K>("data/pseudo_split_0.poly");
  test_offset<K>("data/pseudo_split_10.poly");
  test_offset<K>("data/pseudo_split_11.poly");
  test_offset<K>("data/pseudo_split_12.poly");
  test_offset<K>("data/pseudo_split_13b.poly");
  test_offset<K>("data/pseudo_split_13.poly");
  test_offset<K>("data/pseudo_split_1.poly");
  test_offset<K>("data/pseudo_split_2.poly");
  test_offset<K>("data/pseudo_split_3.poly");
  test_offset<K>("data/pseudo_split_4.poly");
  test_offset<K>("data/pseudo_split_5b.poly");
  test_offset<K>("data/pseudo_split_5.poly");
  test_offset<K>("data/pseudo_split_6.poly");
  test_offset<K>("data/pseudo_split_7.poly");
  test_offset<K>("data/pseudo_split_8.poly");
  test_offset<K>("data/pseudo_split_9.poly");
  test_offset<K>("data/rect_4_spokes.poly");
  test_offset<K>("data/rectangle.poly");
  test_offset<K>("data/region_4.poly");
  test_offset<K>("data/rombus_4_spokes.poly");

#ifndef CGAL_SLS_TEST_SPEED_THINGS_UP_FOR_THE_TESTSUITE
  test_offset<K>("data/sample.poly");
  test_offset<K>("data/sample_0.poly");
  test_offset<K>("data/sample_1.poly");
  test_offset<K>("data/sample_2.poly");
  test_offset<K>("data/sample2.poly");
  test_offset<K>("data/sample_3.poly");
  test_offset<K>("data/sample3.poly");
  test_offset<K>("data/sample_4.poly");
  test_offset<K>("data/sample_5.poly");
  test_offset<K>("data/sample_6.poly");
  test_offset<K>("data/sample_46.poly");
  test_offset<K>("data/sample_73.poly");
  test_offset<K>("data/sample_85.poly");
  test_offset<K>("data/sample_101.poly");
  test_offset<K>("data/sample_102.poly");
  test_offset<K>("data/sample_147.poly");
  test_offset<K>("data/sample_235.poly");
  test_offset<K>("data/sample_298.poly");
  test_offset<K>("data/sample_319.poly");
  test_offset<K>("data/sample_325.poly");
  test_offset<K>("data/sample_333.poly");
  test_offset<K>("data/sample_638.poly");
  test_offset<K>("data/sample_698.poly");
  test_offset<K>("data/simple_0.poly");
  test_offset<K>("data/simple_1.poly");
  test_offset<K>("data/simple_2.poly");
  test_offset<K>("data/simple_3.poly");
  test_offset<K>("data/single_split.poly");
  test_offset<K>("data/split_at_end_0.poly");
  test_offset<K>("data/split_at_end_1.poly");
  test_offset<K>("data/split_at_end_2.poly");
  test_offset<K>("data/split_at_zero_0.poly");
  test_offset<K>("data/square.poly");
  test_offset<K>("data/star.poly");
  test_offset<K>("data/StrayCenterlines.poly");
  test_offset<K>("data/triangle.poly");
  test_offset<K>("data/wheel_13_spokes.poly");
  test_offset<K>("data/wheel_14_spokes.poly");
  test_offset<K>("data/wheel_15_spokes.poly");
  test_offset<K>("data/wheel_16_spokes_b.poly");
  test_offset<K>("data/wheel_16_spokes.poly");
  test_offset<K>("data/wiggly_03_cgal.poly");
  test_offset<K>("data/WingChiu.poly");

  // Below is particularly long
  test_offset<K>("data/alley_0.poly");
  test_offset<K>("data/alley_1.poly");
  test_offset<K>("data/alley_2.poly");
  test_offset<K>("data/alley_3.poly");
  test_offset<K>("data/inputc.poly");
  test_offset<K>("data/inputd1.poly");
  test_offset<K>("data/inputd.poly");
  test_offset<K>("data/inputG.poly");
  test_offset<K>("data/input_K.poly");
  test_offset<K>("data/inputPa.poly");
  test_offset<K>("data/inputP.poly");
  test_offset<K>("data/inputq1.poly");
  test_offset<K>("data/inputq.poly");
  test_offset<K>("data/inputT.poly");
  test_offset<K>("data/inputu.poly");
  test_offset<K>("data/large_1.poly");
  test_offset<K>("data/large_2.poly");
  test_offset<K>("data/large_3.poly");
  test_offset<K>("data/large_4.poly");
  test_offset<K>("data/complex_0.poly");
  test_offset<K>("data/complex_1.poly");
  test_offset<K>("data/complex_2.poly");
  test_offset<K>("data/complex_3.poly");
  test_offset<K>("data/complex_4.poly");
  test_offset<K>("data/complex_5.poly");
  test_offset<K>("data/wheel_128_spokes.poly");
#endif
}

int main(int, char**)
{
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  test_kernel<EPECK_w_sqrt>();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
