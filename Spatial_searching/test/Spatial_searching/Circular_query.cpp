// test whether circular queries are computed correctly for random data
//
// 1) generate list of query points using report_all
// 2) remove and check reported points from these list
// 3) check if no remaining points should have been reported

#include "Point_with_info.h"

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_sphere.h>

#include <CGAL/algorithm.h>
#include <CGAL/iterator.h>
#include <CGAL/Origin.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>
#include <CGAL/use.h>

#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <cassert>
#include <iostream>
#include <list>
#include <vector>

typedef CGAL::Simple_cartesian<double>                                 K;
typedef K::FT                                                          FT;
typedef K::Point_2                                                     Point;
typedef K::Vector_2                                                    Vector;

typedef CGAL::Random_points_in_square_2<Point>                         Random_points_iterator;
typedef CGAL::Random_points_on_circle_2<Point>                         Random_points_on_circle_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator>                N_Random_points_iterator;
typedef CGAL::Search_traits_2<K>                                       Traits;

//for Point_with_info
typedef Point_with_info_helper<Point>::type                            Point_with_info;
typedef Point_property_map<Point>                                      Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,Traits>      Traits_with_info;

template <class SearchTraits>
void run_with_fuzziness(std::list<Point> all_points, // intentional copy
                        const Point& center,
                        const FT radius,
                        const FT fuzziness,
                        CGAL::Random& rnd)
{
  typedef CGAL::Fuzzy_sphere<SearchTraits> Fuzzy_circle;

  std::cout << "test with center: " << center << " radius: " << radius << " eps: " << fuzziness << "... ";

  const FT inner_radius = radius - fuzziness;
  const FT sq_inner_radius = (inner_radius < FT(0)) ? FT(-1) : inner_radius * inner_radius;
  const FT outer_radius = radius + fuzziness;
  const FT sq_outer_radius = outer_radius * outer_radius;

  // add a bunch of points on the boundary of the inner approximation
  // (in general all the constructed points won't be exactly on the boundary,
  // but as long as some are, that's okay)
  std::size_t N = all_points.size();
  Random_points_on_circle_iterator rpocit(inner_radius, rnd);
  for(std::size_t i=0, max=N/10; i<max; ++i)
    all_points.push_back(*rpocit++ + Vector(CGAL::ORIGIN, center));

  // same for the outer approximation boundary
  Random_points_on_circle_iterator rpocit2(outer_radius, rnd);
  for(std::size_t i=0, max = N/10; i<max; ++i)
    all_points.push_back(*rpocit2++ + Vector(CGAL::ORIGIN, center));

  typedef CGAL::Kd_tree<SearchTraits> Tree;

  // Insert also the N points in the tree
  Tree tree(
    boost::make_transform_iterator(all_points.begin(), Create_point_with_info<typename SearchTraits::Point_d>()),
    boost::make_transform_iterator(all_points.end(), Create_point_with_info<typename SearchTraits::Point_d>())
  );

  Fuzzy_circle default_range(typename SearchTraits::Point_d(center), radius);
  std::list<typename SearchTraits::Point_d> result;
  tree.search(std::back_inserter(result), default_range);

  typedef std::vector<typename SearchTraits::Point_d> V;
  V vec;
  vec.resize(result.size());
  typename V::iterator it = tree.search(vec.begin(), default_range);
  CGAL_USE(it);
  assert(it == vec.end());
  result.clear();

  tree.search(CGAL::Emptyset_iterator(), Fuzzy_circle(center, radius) ); //test compilation when Point != Traits::Point_d

  // range searching
  Fuzzy_circle approximate_range(typename SearchTraits::Point_d(center), radius, fuzziness);
  tree.search(std::back_inserter(result), approximate_range);
  std::cout << result.size() << " hits... Verifying correctness...";

  // test the results
  for(typename std::list<typename SearchTraits::Point_d>::iterator pt=result.begin(); (pt != result.end()); ++pt)
  {
    // a point with distance d to the center can only be reported if d <= r + eps
    bool is_correct = (CGAL::squared_distance(center, get_point(*pt)) <= sq_outer_radius);
    if(!is_correct)
      std::cout << get_point(*pt) << " at distance = " << CGAL::squared_distance(center, get_point(*pt))
                << " should not have been reported (max is: " << sq_outer_radius << ")" << std::endl;
    assert(is_correct);
    all_points.remove(get_point(*pt));
  }

  // note: if eps > r, the squared radius is set to '-1'. We cannot have missed
  // any point because there is no inner approximation.
  if(sq_inner_radius >= FT(0))
  {
    for(std::list<Point>::const_iterator pt=all_points.begin(); (pt != all_points.end()); ++pt)
    {
      // all points with a distance d <= r - eps must have been reported
      bool is_correct = (CGAL::squared_distance(center,*pt) > sq_inner_radius);
      if(!is_correct)
        std::cout << "missed " << *pt << " with distance = " << CGAL::squared_distance(center,*pt) << std::endl;

      assert(is_correct);
    }
  }

  std::cout << "done" << std::endl;
}

template <class SearchTraits>
void run(std::list<Point> all_points, // intentional copy
         CGAL::Random& rnd)
{
  Point center(0.25, 0.25);

  // add some interesting points
  all_points.push_back(center);

  run_with_fuzziness<SearchTraits>(all_points, center, 0. /*radius*/, 0. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, center, 0.25 /*radius*/, 0. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, center, 0.25 /*radius*/, 0.125 /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, center, 0.25 /*radius*/, 0.25 /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, center, 0.25 /*radius*/, 1. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, center, 1. /*radius*/, 0. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, center, 1. /*radius*/, 0.25 /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, center, 10. /*radius*/, 0. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, center, 10. /*radius*/, 10. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, center, 10. /*radius*/, 100. /*fuzziness*/, rnd);
}

int main()
{
  CGAL::Set_ieee_double_precision pfr;

  CGAL::Random rnd = CGAL::get_default_random();
  std::cout << "seed: " << rnd.get_seed() << std::endl;

  // generator for random data points in the square ( (-1,-1), (1,1) )
  Random_points_iterator rpit(1.0, rnd);

  // construct list containing N random points
  const int N=1000;
  std::list<Point> all_points(N_Random_points_iterator(rpit,0),
                              N_Random_points_iterator(N));

  run<Traits>(all_points, rnd);
  run<Traits_with_info>(all_points, rnd);

  return 0;
}
