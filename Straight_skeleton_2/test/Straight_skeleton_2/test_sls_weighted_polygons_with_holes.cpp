#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Straight_skeleton_2/IO/print.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/draw_surface_mesh.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/random_polygon_2.h>

#include <CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/create_weighted_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/extrude_skeleton.h>

#include <boost/shared_ptr.hpp>

#include <iostream>

namespace SS = CGAL::CGAL_SS_i;

typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel            EPECK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  EPECK_w_sqrt;

template <typename K>
void test_kernel(const int hole_n, const int hole_nv, CGAL::Random& rnd)
{
  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;
  using Vector_2 = typename K::Vector_2;
  using Point_3 = typename K::Point_3;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
  using Straight_skeleton_2_ptr = boost::shared_ptr<Straight_skeleton_2>;

  using Mesh = CGAL::Surface_mesh<Point_3>;

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

  Mesh sm;
  bool success = extrude_skeleton(pwh, sm, CGAL::parameters::weights(weights));
  assert(success);
  if(!success)
  {
    std::cerr << "Error: failed to extrude skeleton" << std::endl;
    return;
  }

  std::cout << num_vertices(sm) << " vertices and " << num_faces(sm) << " faces" << std::endl;

//  CGAL::draw(sm);
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
