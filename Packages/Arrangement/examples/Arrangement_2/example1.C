// examples/Arrangement_2/example1.C
// ---------------------------------

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;

typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;

typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits>                Arr_2;

int main() 
{
  Arr_2 arr;

  //insertion of the curves
  arr.insert(Curve_2(Point_2(0, 0), Point_2(1, 1)));
  arr.insert(Curve_2(Point_2(0, 1), Point_2(1, 0))); 
  
  //traversal of the curves
  Arr_2::Curve_iterator cit;
  Arr_2::Edge_iterator eit;
  for (cit = arr.curve_node_begin(); cit != arr.curve_node_end(); ++cit) 
  {
    std::cout << std::endl << "Curve level:" << std::endl << cit->curve()
              << std::endl ;
    std::cout << "Edge level:" << std::endl;

    //traversal of the edges of the current curve
    for (eit = cit->edges_begin(); eit != cit->edges_end(); ++eit) 
    {
      std::cout << eit->x_curve() << std::endl ;
    }
  }

  return 0;
}
