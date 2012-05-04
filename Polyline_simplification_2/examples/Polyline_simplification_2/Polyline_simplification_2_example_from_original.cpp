// Recommended kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
 
// The simplification data structure
#include <CGAL/Polyline_simplification_2.h>

// Stop predicate
#include <CGAL/Polyline_simplification_2/Stop_below_count_ratio_threshold.h>

// Cost functions
#include <CGAL/Polyline_simplification_2/Scaled_squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2 Point ;

typedef CGAL::Simplify_polylines_2<K> PS;

int main( int argc, char** argv ) 
{
  PS ps ;
  
  Point points[] = { Point(0,1)
                   , Point(1,2)
                   , Point(2,1)
                   , Point(3,3)
                   , Point(4,1)
                   , Point(5,4)
                   , Point(6,1)
                   , Point(3,0)
                   }  ;
       
  // Insert polygon into simplification class
  PS::Closed_polyline_handle poly = ps.insert_polygon(points,points+8);
  
  std::cout << "Before simplification" << std::endl ;
  std::cout << std::endl ;

  // Iterate over the vertices of the inserted polygon before simplification
  PS::Closed_polyline::Circulator head = poly->circulator();
  PS::Closed_polyline::Circulator cit = head;
  
  do
  {
     std::cout << cit->point() << std::endl ;
  }
  while ( ++ cit != head ) ;
  
  // Define the stop predicate to finish when the number of vertices drops 
  // below half the initial amount
  CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold stop(0.5);
  
  CGAL::Polyline_simplification_2::Scaled_squared_distance_cost<double> cost1;
  CGAL::Polyline_simplification_2::Squared_distance_cost<double>        cost2;
  
  // Proceed with the simplification   
  ps.simplify(stop, cost1) ;

  std::cout << "Using the scaled squared distance" << std::endl ;
  std::cout << std::endl ;

  // Iterate over the now simplified polyline showing the remaining vertices
  head = poly->circulator();
  cit = head;
  
  do
  {
     std::cout << cit->point() << std::endl ;
  }
  while ( ++ cit != head ) ;

  // Start again but using the other cost function
  ps.simplify_original(stop, cost2) ;

  std::cout << "Using the absolute squared distance" << std::endl ;
  std::cout << std::endl ;

  // Iterate over the now simplified polyline showing the remaning vertices
  head = poly->circulator();
  cit = head;
  
  do
  {
     std::cout << cit->point() << std::endl ;
  }
  while ( ++ cit != head ) ;

  return 0 ;      
}


