// file: examples/Arrangement_2/example4.C

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits>                Arr_2;

// A simple function that splits a segment into 2
void my_split_f(const Curve_2 & cv, std::list<Curve_2> & l) 
{
  Point_2 s = cv.source(); // Uses the knowledge of the curve functions
  Point_2 t = cv.target();
  Point_2 m1 = s + (t - s) / 2.0;
  l.push_back(Curve_2(s, m1));
  l.push_back(Curve_2(m1, t));
}

typedef void (*SPLIT_FUNC)(const Curve_2 & cv, std::list<Curve_2> & l);

int main() 
{
  std::vector<SPLIT_FUNC> func_vec;
  func_vec.push_back(&my_split_f);
  Arr_2 arr;
   
  // Insertion with user-defined function
  Arr_2::Curve_iterator cit = arr.insert(Curve_2(Point_2(0, 0),Point_2(6, 6)),
                                         func_vec.begin(), func_vec.end());
   
  // Regular insertion                                  
  cit = arr.insert(Curve_2(Point_2(0, 4), Point_2(6, 4))); 
   
  // Traversal of the curves
  Arr_2::Edge_iterator eit;
  for (cit = arr.curve_node_begin(); cit != arr.curve_node_end(); ++cit) 
  {
    std::cout << std::endl << "Curve level:" << std::endl << cit->curve()
              << std::endl ;
    std::cout << "Edge level:" << std::endl;
    for (eit = cit->edges_begin(); eit != cit->edges_end(); ++eit) 
    {
      std::cout << eit->x_curve() << std::endl ;
    }
  }

  return 0;
}
