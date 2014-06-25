
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/random_convex_hull_in_disc_2.h>
#include <iostream>


  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_2                                          Point_2;
  typedef CGAL::Polygon_2<K>                                  Polygon_2;

int
main( )
{


  Polygon_2 p;
  int n( 1000);
  boost::random::mt19937 gen;

  // build random hull from n random points in a disc:
  random_convex_hull_in_disc_2(n,1.0,gen,std::back_inserter(p),K());

  // check convexity:
  if ( ! p.is_convex()) {
    std::cerr << "ERROR: polygon is not convex." << std::endl;
    return 1;
  }
  return 0;
} // int main( )

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

