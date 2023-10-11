#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Fuzzy_iso_box.h>
#include "Point_with_info.h"

#include <cassert>
#include <vector>
#include <iostream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Iso_cuboid_3 Iso_cuboid;
typedef CGAL::Random_points_in_sphere_3<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_3<K>                        Traits;
//for Point_with_info
typedef Point_with_info_helper<Point>::type                                          Point_with_info;
typedef Point_property_map<Point>                                                    Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,Traits>                    Traits_with_info;

const int N=10000;

template <class SearchTraits>
void run(std::list<Point> all_points)
{
  typedef CGAL::Fuzzy_iso_box<SearchTraits> Fuzzy_box;

  // Insert also the N points in the tree
  CGAL::Kd_tree<SearchTraits> tree(
    boost::make_transform_iterator(all_points.begin(),Create_point_with_info<typename SearchTraits::Point_d>()),
    boost::make_transform_iterator(all_points.end(),Create_point_with_info<typename SearchTraits::Point_d>())
  );

  Point p(0.1, 0.2, 0.3);
  Point q(0.3, 0.5, 0.4);

  typename SearchTraits::Point_d pp(p);
  typename SearchTraits::Point_d qq(q);

  Iso_cuboid exact_ic(p,q);

  Fuzzy_box exact_range(pp,qq);

  std::list<typename SearchTraits::Point_d> result;
  // Searching the box r exactly
  tree.search( std::back_inserter( result ), exact_range);

  // test the results of the exact query
  std::list<Point> copy_all_points(all_points);

  for (typename std::list<typename SearchTraits::Point_d>::iterator pt=result.begin(); (pt != result.end()); ++pt) {
    assert(! exact_ic.has_on_unbounded_side(get_point(*pt)));
    copy_all_points.remove(get_point(*pt));
  }

  for (std::list<Point>::iterator pt=copy_all_points.begin(); (pt != copy_all_points.end()); ++pt) {
    assert(exact_ic.has_on_unbounded_side(*pt));
  }


  result.clear();
  // approximate range searching using value 0.05 for fuzziness parameter
  Fuzzy_box approximate_range(pp,qq,0.05);
  Iso_cuboid inner_ic(p+ 0.05*Vector(1,1,1),q-0.05*Vector(1,1,1));
  Iso_cuboid outer_ic(p- 0.05*Vector(1,1,1),q+0.05*Vector(1,1,1));

  tree.search(std::back_inserter( result ), approximate_range);
  // test the results of the approximate query
  for (typename std::list<typename SearchTraits::Point_d>::iterator pt=result.begin(); (pt != result.end()); ++pt) {
    // a point we found may be slightly outside the isocuboid
    assert(! outer_ic.has_on_unbounded_side(get_point(*pt)));
    all_points.remove(get_point(*pt));
  }

  for (std::list<Point>::iterator pt=all_points.begin(); (pt != all_points.end()); ++pt) {
    assert(inner_ic.has_on_unbounded_side(*pt));
  }
  std::cout << "done" << std::endl;
}

int main() {
  // generator for random data points in the square ( (-1,-1), (1,1) )
  Random_points_iterator rpit( 1.0);

  // construct list containing N random points
  std::list<Point> all_points(N_Random_points_iterator(rpit,0),
                              N_Random_points_iterator(N));

  run<Traits>(all_points);
  run<Traits_with_info>(all_points);

  return 0;
}

