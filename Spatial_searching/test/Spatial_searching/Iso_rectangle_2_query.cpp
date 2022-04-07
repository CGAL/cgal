// test whether box queries are computed correctly for random data
//
// 1) generate list of query points using report_all
// 2) remove and check reported points from these list
// 3) check if no remaining points should have been reported

#include <CGAL/Simple_cartesian.h>

#include "Point_with_info.h"

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_iso_box.h>

#include <CGAL/iterator.h>
#include <CGAL/point_generators_2.h>

#include <cassert>
#include <vector>
#include <iostream>

typedef CGAL::Simple_cartesian<double>                                K;
typedef K::FT                                                         FT;
typedef K::Point_2                                                    Point;
typedef K::Vector_2                                                   Vector;
typedef K::Iso_rectangle_2                                            Iso_rectangle;

typedef CGAL::Random_points_in_square_2<Point>                        Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator>               N_Random_points_iterator;
typedef CGAL::Search_traits_2<K>                                      Traits;
typedef Point_with_info_helper<Point>::type                           Point_with_info;
typedef Point_property_map<Point>                                     Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,Traits>     Traits_with_info;

template <class SearchTraits>
void run_with_fuzziness(std::list<Point> all_points,
                        const Point& p, const Point& q,
                        const FT fuzziness,
                        CGAL::Random& rnd)
{
  typedef CGAL::Fuzzy_iso_box<SearchTraits> Fuzzy_box;

  std::cout << "test with box: [" << p << " || " << q << "] and eps: " << fuzziness << "... ";

  // test the results of the approximate query
  Iso_rectangle inner_ic(p + fuzziness*Vector(1,1), q - fuzziness*Vector(1,1));
  Iso_rectangle outer_ic(p - fuzziness*Vector(1,1), q + fuzziness*Vector(1,1));

  // If the fuziness is greater than half of the largest dimension of the box,
  // then the inner box does not exist
  const FT max_box_edge_length = (std::max)(q[1] - p[1], q[0] - p[0]);
  const bool is_inner_c_empty = (fuzziness > 0.5 * max_box_edge_length);
  if(is_inner_c_empty)
    std::cout << " (empty inner box)... ";

  // Insert a bunch of points on the boundary of the inner approximation
  std::size_t N = all_points.size();

  if(!is_inner_c_empty)
  {
    for(std::size_t i=0, max=N/10; i<max; ++i)
    {
      double x = rnd.get_double(-1., 1.);
      all_points.push_back(Point(x, p.y() + fuzziness));
      all_points.push_back(Point(p.x() + fuzziness, x));
    }
  }

  // same for the outer boundary
  for(std::size_t i=0, max=N/10; i<max; ++i)
  {
    double x = rnd.get_double(-1., 1.);
    all_points.push_back(Point(x, q.y() + fuzziness));
    all_points.push_back(Point(q.x() + fuzziness, x));
  }

  // Insert also the N points in the tree
  CGAL::Kd_tree<SearchTraits> tree(
    boost::make_transform_iterator(all_points.begin(), Create_point_with_info<typename SearchTraits::Point_d>()),
    boost::make_transform_iterator(all_points.end(), Create_point_with_info<typename SearchTraits::Point_d>())
  );

  tree.search(CGAL::Emptyset_iterator(), Fuzzy_box(p,q)); //test compilation when Point != Traits::Point_d

  typename SearchTraits::Point_d pp(p);
  typename SearchTraits::Point_d qq(q);

  // approximate range searching
  std::list<typename SearchTraits::Point_d> result;
  Fuzzy_box approximate_range(pp, qq, fuzziness);
  tree.search(std::back_inserter(result), approximate_range);
  std::cout << result.size() << " hits... Verifying correctness...";

  for (typename std::list<typename SearchTraits::Point_d>::iterator pt=result.begin(); (pt != result.end()); ++pt)
  {
    // a point can only be reported if it is in the outer box
    bool is_correct = outer_ic.has_on_bounded_side(get_point(*pt)) || outer_ic.has_on_boundary(get_point(*pt));
    if(!is_correct)
      std::cout << get_point(*pt) << " should have not been reported" << std::endl;
    assert(is_correct);
    all_points.remove(get_point(*pt));
  }

  // nothing to test if the inner box is empty because everything is on the unbounded side
  if(!is_inner_c_empty)
  {
    for (std::list<Point>::iterator pt=all_points.begin(); (pt != all_points.end()); ++pt)
    {
      // all points that have not been reported must be outside the inner box
      bool is_correct = inner_ic.has_on_unbounded_side(*pt);
      if(!is_correct)
        std::cout << *pt << " should have been reported" << std::endl;

      assert(is_correct);
    }
  }

  std::cout << "done" << std::endl;
}

template <class SearchTraits>
void run(std::list<Point> all_points, // intentional copy
         CGAL::Random& rnd)
{
  // bigger box
  Point p0(-10., -10.);
  Point q0( 10.,  10.);

  // a subset
  Point p1(-CGAL_PI/10., -CGAL_PI/10.);
  Point q1( CGAL_PI/10.,  CGAL_PI/10.);

  // another subset
  Point p2(0.1, 0.2);
  Point q2(0.3, 0.4);

  // degenerate
  Point p3(0., 0.);
  Point q3(0., 0.);

  run_with_fuzziness<SearchTraits>(all_points, p0, q0, 0. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, p0, q0, 0.1 /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, p0, q0, 1. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, p0, q0, 10. /*fuzziness*/, rnd);

  run_with_fuzziness<SearchTraits>(all_points, p1, q1, 0. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, p1, q1, 0.1 /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, p1, q1, 1. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, p1, q1, 10. /*fuzziness*/, rnd);

  run_with_fuzziness<SearchTraits>(all_points, p2, q2, 0. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, p2, q2, 0.1 /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, p2, q2, 0.4 /*fuzziness*/, rnd);

  run_with_fuzziness<SearchTraits>(all_points, p3, q3, 0. /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, p3, q3, 0.33 /*fuzziness*/, rnd);
  run_with_fuzziness<SearchTraits>(all_points, p3, q3, 1. /*fuzziness*/, rnd);
}

int main()
{
  const int N=10000;

  CGAL::Random rnd = CGAL::get_default_random();
  std::cout << "seed: " << rnd.get_seed() << std::endl;

  // generator for random data points in the square ( (-1,-1), (1,1) )
  Random_points_iterator rpit(1.0, rnd);

  // construct list containing N random points
  std::list<Point> all_points(N_Random_points_iterator(rpit,0),
                              N_Random_points_iterator(N));

  // add some interesting points
  all_points.push_back(Point(0., 0.));
  all_points.push_back(Point(-CGAL_PI/10.+0.1, -CGAL_PI/10.+0.1));
  all_points.push_back(Point(1., 1.));
  all_points.push_back(Point(0., 1.));
  all_points.push_back(Point(0.3, 0.4));
  all_points.push_back(Point(0.2, 0.3));
  all_points.push_back(Point(0., 0.1));

  run<Traits>(all_points, rnd);
  run<Traits_with_info>(all_points, rnd);

  return 0;
}


