#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Straight_skeleton_2/IO/print.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_surface_mesh.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/random_polygon_2.h>

#include <CGAL/create_weighted_straight_skeleton_2.h>
#include <CGAL/create_weighted_offset_polygons_2.h>
#include <CGAL/extrude_skeleton.h>

#include <boost/shared_ptr.hpp>

#include <iostream>

namespace SS = CGAL::CGAL_SS_i;

typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel            EPECK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  EPECK_w_sqrt;

template <typename K>
void test_API()
{
  typedef typename K::FT                                             FT;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;

  typedef CGAL::Straight_skeleton_2<EPICK>                           Straight_skeleton_EPICK;
  typedef boost::shared_ptr<Straight_skeleton_EPICK>                 Straight_skeleton_Ptr_EPICK;

  typedef CGAL::Straight_skeleton_2<K>                               Straight_skeleton;
  typedef boost::shared_ptr<Straight_skeleton>                       Straight_skeleton_Ptr;

  std::vector<std::vector<FT> > weights;

  Polygon_2 p;
  Straight_skeleton_Ptr_EPICK ss_epick = CGAL::create_interior_weighted_straight_skeleton_2(p, weights);
  Straight_skeleton_Ptr ss = CGAL::create_interior_weighted_straight_skeleton_2(p, weights, K());
  ss_epick = CGAL::create_exterior_weighted_straight_skeleton_2(double(1.01), p, weights);
  ss = CGAL::create_exterior_weighted_straight_skeleton_2(int(2), p, weights, K());

  Polygon_with_holes_2 pwh;
  ss_epick = CGAL::create_interior_weighted_straight_skeleton_2(pwh, weights);
  ss = CGAL::create_interior_weighted_straight_skeleton_2(pwh, weights, K());
  ss_epick = CGAL::create_exterior_weighted_straight_skeleton_2(double(1.01), p, weights);
  ss = CGAL::create_exterior_weighted_straight_skeleton_2(int(2), p, weights, K());
}

template <typename K>
void test_kernel(const int polygon_nv, CGAL::Random& rnd)
{
  using FT = typename K::FT;
  using Point_2 = typename  K::Point_2;
  using Vector_2 = typename  K::Vector_2;
  using Point_3 = typename K::Point_3;

  using Polygon_2 = CGAL::Polygon_2<K>;

  using Straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
  using Straight_skeleton_2_ptr = boost::shared_ptr<Straight_skeleton_2>;

  using Mesh = CGAL::Surface_mesh<Point_3>;

  void (*dummy_ptr)() = &test_API<K>;

  typedef CGAL::Random_points_in_square_2<Point_2> Point_generator;
  Polygon_2 pol;
  CGAL::random_polygon_2(polygon_nv, std::back_inserter(pol), Point_generator(0.25, rnd));

  std::vector<std::vector<FT> > weights(1);
  for(int i=0; i<polygon_nv; ++i)
    weights[0].push_back(rnd.get_double(1, 10));

  std::cout << "Weights:" << std::endl;
  for(const std::vector<FT>& ws : weights)
  {
    for(FT w : ws)
      std::cout << w << " ";
    std::cout << std::endl;
  }

  CGAL::draw(pol);

  auto ss_ptr = CGAL::create_interior_weighted_straight_skeleton_2(pol, weights);
  assert(ss_ptr);
  if(!ss_ptr)
  {
    std::cerr << "Error: failed to create straight skeleton" << std::endl;
    return;
  }

  CGAL::draw(*ss_ptr);

  ss_ptr = CGAL::create_exterior_weighted_straight_skeleton_2(0.1, pol, weights);
  assert(ss_ptr);
  if(!ss_ptr)
  {
    std::cerr << "Error: failed to create straight skeleton" << std::endl;
    return;
  }

  CGAL::draw(*ss_ptr);

  Mesh sm;
  bool success = extrude_skeleton(pol, sm, CGAL::parameters::weights(weights));
  assert(success);
  if(!success)
  {
    std::cerr << "Error: failed to extrude skeleton" << std::endl;
    return;
  }

  std::cout << num_vertices(sm) << " vertices and " << num_faces(sm) << " faces" << std::endl;

  CGAL::draw(sm);
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const int polygon_nv = (argc > 1) ? std::atoi(argv[1]) : 10;
  const int seed = (argc > 2) ? std::atoi(argv[2]) : std::time(nullptr);

  CGAL::Random rnd(seed);
  std::cout << "Seed is " << rnd.get_seed() << std::endl;

  test_kernel<EPICK>(polygon_nv, rnd);
  test_kernel<EPECK>(polygon_nv, rnd);
  test_kernel<EPECK_w_sqrt>(polygon_nv, rnd);

  return EXIT_SUCCESS;
}
