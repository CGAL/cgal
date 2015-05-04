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
#include <CGAL/iterator.h>
#include "Point_with_info.h"

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K>                Traits;
//for Point_with_info
typedef Point_with_info_helper<Point>::type                                          Point_with_info;
typedef Point_property_map<Point>                                                    Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,Traits>                    Traits_with_info;


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
  Point center(0.25, 0.25);
  Fuzzy_circle exact_range(typename Traits::Point_d(center), 0.25);
    
  std::list<typename Traits::Point_d> result;
  tree.search(std::back_inserter( result ), exact_range);

  typedef std::vector<typename Traits::Point_d> V;
  V vec;
  vec.resize(result.size());
  typename V::iterator it = tree.search(vec.begin(), exact_range);
  assert(it == vec.end());
 
  tree.search(CGAL::Emptyset_iterator(), Fuzzy_circle(center, 0.25) ); //test compilation when Point != Traits::Point_d
  
  // test the results of the exact query
  std::list<Point> copy_all_points(all_points);
  for (typename std::list<typename Traits::Point_d>::iterator pt=result.begin(); (pt != result.end()); ++pt) {
    assert(CGAL::squared_distance(center,get_point(*pt))<=0.0625);
    copy_all_points.remove(get_point(*pt));
  }
  
  for (std::list<Point>::iterator pt=copy_all_points.begin(); (pt != copy_all_points.end()); ++pt) {
    if(CGAL::squared_distance(center,*pt)<=0.0625){
      std::cout << "we missed " << *pt << " with distance = " << CGAL::squared_distance(center,*pt) << std::endl;
    }
    assert(CGAL::squared_distance(center,*pt)>0.0625);
  }


  result.clear();
  // approximate range searching using value 0.125 for fuzziness parameter
  Fuzzy_circle approximate_range(typename Traits::Point_d(center), 0.25, 0.125);

  tree.search(std::back_inserter( result ), approximate_range);
  // test the results of the approximate query
  for (typename std::list<typename Traits::Point_d>::iterator pt=result.begin(); (pt != result.end()); ++pt) {
    // a point with distance d to the center may be reported if d <= r + eps
    assert(CGAL::squared_distance(center,get_point(*pt))<=0.140625); // (0.25 + 0.125)²
    all_points.remove(get_point(*pt));
  }
  
  for (std::list<Point>::iterator pt=all_points.begin(); (pt != all_points.end()); ++pt) {
    if(CGAL::squared_distance(center,*pt)<=0.015625){ // 0.125²
      std::cout << "we missed " << *pt << " with distance = " << CGAL::squared_distance(center,*pt) << std::endl;
    }
    //all points with a distance d <= r-eps must be reported
    assert(CGAL::squared_distance(center,*pt)> 0.015625);
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

