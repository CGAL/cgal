// file          : test/Spatial_searching/Circular_query.C
// test whether circular queries are computed correctly for random data
// 
// 1) generate list of query points using report_all  
// 2) remove and check reported points from these list
// 3) check if no remaining points should have been reported

#include <CGAL/Cartesian.h>
#include <cassert>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include "Point_with_info.h"

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K>                Traits;
//for Point_with_info
typedef Point_with_info_helper<Point>::type                                          Point_with_info;
typedef CGAL::Search_traits_adapter<Point_with_info,Point_property_map,Traits>         Traits_with_info;


template <class Traits>
void run(std::list<Point> all_points){
  typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_circle;
  typedef CGAL::Kd_tree<Traits> Tree;  
  
  // Insert also the N points in the tree
  Tree tree(
    boost::make_transform_iterator(all_points.begin(),Create_point_with_info<typename Traits::Point_d>()),
    boost::make_transform_iterator(all_points.end(),Create_point_with_info<typename Traits::Point_d>())
  );

  // define exact circular range query  (fuzziness=0)
  Point center(0.2, 0.2);
  Fuzzy_circle exact_range(typename Traits::Point_d(center), 0.2);
    
  std::list<typename Traits::Point_d> result;
  tree.search(std::back_inserter( result ), exact_range);
  
  // test the results of the exact query
  std::list<Point> copy_all_points(all_points);
  for (typename std::list<typename Traits::Point_d>::iterator pt=result.begin(); (pt != result.end()); ++pt) {
    assert(CGAL::squared_distance(center,get_point(*pt))<=0.04);
    copy_all_points.remove(get_point(*pt));
  }
  
  for (std::list<Point>::iterator pt=copy_all_points.begin(); (pt != copy_all_points.end()); ++pt) {
    if(CGAL::squared_distance(center,*pt)<=0.04){
      std::cout << "we missed " << *pt << " with distance = " << CGAL::squared_distance(center,*pt) << std::endl;
    }
    assert(CGAL::squared_distance(center,*pt)>0.04);
  }


  result.clear();
  // approximate range searching using value 0.1 for fuzziness parameter
  Fuzzy_circle approximate_range(typename Traits::Point_d(center), 0.2, 0.1);

  tree.search(std::back_inserter( result ), approximate_range);
  // test the results of the approximate query
  for (typename std::list<typename Traits::Point_d>::iterator pt=result.begin(); (pt != result.end()); ++pt) {
    // a point we found may be slighlty outside the circle
    assert(CGAL::squared_distance(center,get_point(*pt))<=0.09);
    all_points.remove(get_point(*pt));
  }
  
  for (std::list<Point>::iterator pt=all_points.begin(); (pt != all_points.end()); ++pt) {
    if(CGAL::squared_distance(center,*pt)<=0.01){
      std::cout << "we missed " << *pt << " with distance = " << CGAL::squared_distance(center,*pt) << std::endl;
    }
    assert(CGAL::squared_distance(center,*pt)> 0.01);
  }
  std::cout << "done" << std::endl;  
}

int main() {

  const int N=1000;

  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit( 1.0);

  // construct list containing N random points
  std::list<Point> all_points(N_Random_points_iterator(rpit,0),
			      N_Random_points_iterator(N));
  
  run<Traits>(all_points);
  run<Traits_with_info>(all_points);

  return 0;
}

