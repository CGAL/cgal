// file: Iso_rectangle_2_query.C

#include <CGAL/Simple_cartesian.h>
#include <cassert>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/iterator.h>
#include "Point_with_info.h"

#include <vector>
#include <iostream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point;
typedef K::Vector_2 Vector;
typedef K::Iso_rectangle_2 Iso_rectangle;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K> Traits;
typedef Point_with_info_helper<Point>::type                                   Point_with_info;
typedef Point_property_map<Point>                                             Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,Traits>             Traits_with_info;

template <class SearchTraits>
void run(std::list<Point> all_points)
{
  typedef CGAL::Fuzzy_iso_box<SearchTraits> Fuzzy_box;
  
  // Insert also the N points in the tree
  CGAL::Kd_tree<SearchTraits> tree(
    boost::make_transform_iterator(all_points.begin(),Create_point_with_info<typename SearchTraits::Point_d>()),
    boost::make_transform_iterator(all_points.end(),Create_point_with_info<typename SearchTraits::Point_d>())
  );
  
  Point p(0.1, 0.2);
  Point q(0.3, 0.5);
  typename SearchTraits::Point_d pp(p);
  typename SearchTraits::Point_d qq(q);
  
  
  Iso_rectangle exact_ic(p,q);

  Fuzzy_box exact_range(pp,qq);

  std::list<typename SearchTraits::Point_d> result;
  // Searching the box r exactly
  tree.search( std::back_inserter( result ), exact_range);
  
  tree.search(CGAL::Emptyset_iterator(), Fuzzy_box(p,q) ); //test compilation when Point != Traits::Point_d
  
  // test the results of the exact query
  std::list<Point> copy_all_points(all_points);
  for (typename std::list<typename SearchTraits::Point_d>::iterator pt=result.begin(); (pt != result.end()); ++pt) {
    assert(! exact_ic.has_on_unbounded_side(get_point(*pt)) || exact_ic.has_on_boundary(get_point(*pt)));
    copy_all_points.remove(get_point(*pt));
  }
  
  for (std::list<Point>::iterator pt=copy_all_points.begin(); (pt != copy_all_points.end()); ++pt) {
    assert(exact_ic.has_on_unbounded_side(*pt) || exact_ic.has_on_boundary(*pt));
  }


  result.clear();
  // approximate range searching using value 0.1 for fuzziness parameter
  Fuzzy_box approximate_range(pp,qq,0.05);
  Iso_rectangle inner_ic(p+ 0.05*Vector(1,1),q-0.05*Vector(1,1));
  Iso_rectangle outer_ic(p- 0.05*Vector(1,1),q+0.05*Vector(1,1));

  tree.search(std::back_inserter( result ), approximate_range);
  // test the results of the approximate query
  for (typename std::list<typename SearchTraits::Point_d>::iterator pt=result.begin(); (pt != result.end()); ++pt) {
    // a point we found may be slighlty outside the isorectangle
    assert(! outer_ic.has_on_unbounded_side(get_point(*pt)) || outer_ic.has_on_boundary(get_point(*pt)));
    all_points.remove(get_point(*pt));
  }
  
  for (std::list<Point>::iterator pt=all_points.begin(); (pt != all_points.end()); ++pt) {
    assert(inner_ic.has_on_unbounded_side(*pt) || inner_ic.has_on_boundary(*pt));
  }
  std::cout << "done" << std::endl;  
}

int main() {

  const int N=10000;
  
  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit( 1.0);
  
  // construct list containing N random points
  std::list<Point> all_points(N_Random_points_iterator(rpit,0),
			      N_Random_points_iterator(N));
  
  run<Traits>(all_points);
  run<Traits_with_info>(all_points);
  
  return 0;
}


