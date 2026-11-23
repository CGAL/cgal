#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/Straight_skeleton_2/IO/print.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/random_polygon_2.h>

#include <CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/create_weighted_offset_polygons_from_polygon_with_holes_2.h>

#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

namespace SS = CGAL::CGAL_SS_i;

typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel            EPECK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  EPECK_w_sqrt;

// Artificial vertices (see `Straight_skeleton_builder_2::EnforceSimpleConnectedness()`) are on the bisectors
// that they split: their times must be (strictly) between the times of their neighbors along the bisector,
// or equal if the bisector is horizontal, and they must be exactly collinear with these neighbors
// if the kernel has exact constructions.
template <typename K, typename StraightSkeleton>
void check_artificial_vertices(const StraightSkeleton& ss)
{
  using Vertex_const_handle = typename StraightSkeleton::Vertex_const_handle;
  using Halfedge_const_handle = typename StraightSkeleton::Halfedge_const_handle;

  CGAL::Straight_skeleton_builder_traits_2<K> traits;

  // contour vertices have time 0, and other vertices have trisegments
  auto compare_times = [&](Vertex_const_handle u, Vertex_const_handle v) -> CGAL::Comparison_result
  {
    if(u->is_contour() || v->is_contour())
      return (u->is_contour() == v->is_contour()) ? CGAL::EQUAL : (u->is_contour() ? CGAL::SMALLER : CGAL::LARGER);
    return static_cast<CGAL::Comparison_result>(traits.compare_ss_event_times_2_object()(u->trisegment(), v->trisegment()));
  };

  int artificial_vertices = 0;
  for(Vertex_const_handle v = ss.vertices_begin(); v != ss.vertices_end(); ++v)
  {
    if(v->is_contour() || v->has_infinite_time() || !v->trisegment())
      continue;
    if(v->trisegment()->e0().id() != v->trisegment()->e1().id()) // marker of artificial vertices
      continue;

    ++artificial_vertices;

    // the two pieces of the split bisector are the incident halfedges that are not artificial bisectors
    std::vector<Vertex_const_handle> neighbors;
    Halfedge_const_handle h = v->halfedge(), done = h;
    do
    {
      assert(h->vertex() == v);
      if(h->face() != h->opposite()->face())
        neighbors.push_back(h->opposite()->vertex());
      h = h->next()->opposite();
    }
    while(h != done);

    assert(neighbors.size() == 2);
    if(neighbors.size() != 2)
      continue;

    const CGAL::Comparison_result c0 = compare_times(neighbors[0], v);
    const CGAL::Comparison_result c1 = compare_times(v, neighbors[1]);
    if(c0 != c1)
      std::cerr << "Error: artificial vertex V" << v->id() << " is not between V" << neighbors[0]->id()
                << " and V" << neighbors[1]->id() << " (" << c0 << " " << c1 << ")" << std::endl;
    assert(c0 == c1);

    if constexpr(!std::is_same<K, EPICK>::value)
    {
      assert(CGAL::are_strictly_ordered_along_line(neighbors[0]->point(), v->point(), neighbors[1]->point()));
    }
  }

  std::cout << artificial_vertices << " artificial vertices" << std::endl;
}

template <typename K>
void test_kernel(const int hole_n, const int hole_nv, CGAL::Random& rnd)
{
  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
  using Straight_skeleton_2_ptr = std::shared_ptr<Straight_skeleton_2>;

  auto generate_random_polygon = [&](CGAL::Random& rnd) -> Polygon_2
  {
    typedef CGAL::Random_points_in_square_2<Point_2> Point_generator;

    Polygon_2 poly;
    CGAL::random_polygon_2(hole_nv, std::back_inserter(poly), Point_generator(0.25, rnd));
    return poly;
  };

  // each hole is in a square of size 1
  std::vector<Point_2> ob = { Point_2(-hole_n-1, -0.5),
                              Point_2( hole_n+1, -0.5),
                              Point_2( hole_n+1,  0.5),
                              Point_2(-hole_n-1,  0.5) };
  Polygon_2 obp(ob.begin(), ob.end());
  Polygon_with_holes_2 pwh(obp);

  std::cout << "pwh.outer_boundary() = " << pwh.outer_boundary() << std::endl;

  std::vector<std::vector<FT> > weights(1);

  // tiny weight (far-reaching) for vertical sides
  weights[0].push_back(rnd.get_double(0.05, 0.5));
  weights[0].push_back(rnd.get_double(1, 10));

  // half the time, use the same weight on the left and right span as to get a nicer bisector
  // in the middle, which increases the likelihood of multiple splits being required
  if(rnd.get_int(0, 2) == 1)
    weights[0].push_back(rnd.get_double(0.05, 0.5));
  else
    weights[0].push_back(weights[0][0]);
  weights[0].push_back(rnd.get_double(1, 10));

  for(int hole_i=0; hole_i<hole_n; ++hole_i)
  {
    Polygon_2 poly = generate_random_polygon(rnd);
    std::vector<Point_2> hole_pts = poly.container();
    double y_nudge = rnd.get_double(-0.25, 0.25);
    for(Point_2& pt : hole_pts)
      pt = Point_2(pt.x() + hole_i + 1, pt.y() + y_nudge);
    Polygon_2 hole(hole_pts.begin(), hole_pts.end());
    if(hole.is_counterclockwise_oriented())
      hole.reverse_orientation();
    pwh.add_hole(hole);
    weights.push_back(std::vector<FT>(hole.size(), 10));

    // same but for negative x
    poly = generate_random_polygon(rnd);
    hole_pts = poly.container();
    y_nudge = rnd.get_double(-0.25, 0.25);
    for(Point_2& pt : hole_pts)
      pt = Point_2(pt.x() - hole_i - 1, pt.y() + y_nudge);
    hole = Polygon_2(hole_pts.begin(), hole_pts.end());
    if(hole.is_counterclockwise_oriented())
      hole.reverse_orientation();
    pwh.add_hole(hole);
    weights.push_back(std::vector<FT>(hole.size(), 10));
  }

  std::cout << "PWH:\n" << pwh << std::endl;
  std::cout << "Weights:" << std::endl;
  for(const std::vector<FT>& ws : weights)
  {
    for(FT w : ws)
      std::cout << w << " ";
    std::cout << std::endl;
  }

//  CGAL::draw(pwh);

  Straight_skeleton_2_ptr ss_ptr = CGAL::create_interior_weighted_straight_skeleton_2(pwh, weights, K());
  assert(ss_ptr);
  if(!ss_ptr)
  {
    std::cerr << "Error: failed to create straight skeleton" << std::endl;
    return;
  }

//  CGAL::draw(*ss_ptr);

  check_artificial_vertices<K>(*ss_ptr);
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  int hole_n = (argc > 1) ? std::atoi(argv[1]) : 2;
  int hole_nv = (argc > 2) ? std::atoi(argv[2]) : 10;
  int seed = (argc > 3) ? std::atoi(argv[3]) : std::time(nullptr);

  CGAL::Random rnd(seed);

  std::cout << "Seed is " << rnd.get_seed() << std::endl;
  std::cout << 2*hole_n << " holes of size " << hole_nv << std::endl;

  test_kernel<EPICK>(hole_n, hole_nv, rnd);
  test_kernel<EPECK>(hole_n, hole_nv, rnd);
  test_kernel<EPECK_w_sqrt>(hole_n, hole_nv, rnd);

  return EXIT_SUCCESS;
}
