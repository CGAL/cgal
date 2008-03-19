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
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K> Traits;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_circle;
typedef CGAL::Kd_tree<Traits> Tree;
  
int main() {

  const int N=1000;

  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit( 1.0);

  // construct list containing N random points
  std::list<Point> all_points(N_Random_points_iterator(rpit,0),
			      N_Random_points_iterator(N));
  
  // Insert also the N points in the tree
  Tree tree(all_points.begin(),all_points.end());

  // define exact circular range query  (fuzziness=0)
  Point center(0.2, 0.2);
  Fuzzy_circle exact_range(center, 0.2);
    
  std::list<Point> result;
  tree.search(std::back_inserter( result ), exact_range);
  
  // test the results of the exact query
  std::list<Point> copy_all_points(all_points);
  std::list<Point>::iterator pt;
  for (pt=result.begin(); (pt != result.end()); ++pt) {
    assert(CGAL::squared_distance(center,*pt)<=0.04);
    copy_all_points.remove(*pt);
  }
  
  for (pt=copy_all_points.begin(); (pt != copy_all_points.end()); ++pt) {
    if(CGAL::squared_distance(center,*pt)<=0.04){
      std::cout << "we missed " << *pt << " with distance = " << CGAL::squared_distance(center,*pt) << std::endl;
    }
    assert(CGAL::squared_distance(center,*pt)>0.04);
  }


  result.clear();
  // approximate range searching using value 0.1 for fuzziness parameter
  Fuzzy_circle approximate_range(center, 0.2, 0.1);

  tree.search(std::back_inserter( result ), approximate_range);
  // test the results of the approximate query
  for (pt=result.begin(); (pt != result.end()); ++pt) {
    // a point we found may be slighlty outside the circle
    assert(CGAL::squared_distance(center,*pt)<=0.09);
    all_points.remove(*pt);
  }
  
  for (pt=all_points.begin(); (pt != all_points.end()); ++pt) {
    if(CGAL::squared_distance(center,*pt)<=0.01){
      std::cout << "we missed " << *pt << " with distance = " << CGAL::squared_distance(center,*pt) << std::endl;
    }
    assert(CGAL::squared_distance(center,*pt)> 0.01);
  }
  std::cout << "done" << std::endl;
  return 0;
}

