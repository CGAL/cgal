//examples/Arrangement_2/example1.C
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Arrangement_2.h>

// We use here double instead of leda_real to enable compilation without LEDA.
// This is not recommended generally.
// Read more in the README file or in the manual.
typedef double                                        NT;
typedef CGAL::Cartesian<NT>                           R;
typedef CGAL::Arr_segment_exact_traits<R>             Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Arr_base_node<Curve>                    Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>              Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >   Arr_2;

int main() {
  
  Arr_2 arr;

  //insertion of the curves
  Arr_2::Curve_iterator cit=arr.insert(Curve(Point(0,0),Point(1,1)));
  cit=arr.insert(Curve(Point(0,1),Point(1,0))); 
  
  //traversal of the curves
  Arr_2::Edge_iterator eit;
  for (cit=arr.curve_node_begin(); cit!=arr.curve_node_end(); ++cit) {
    std::cout << std::endl << "Curve level:" << std::endl << cit->curve()
	      << std::endl ;
    std::cout << "Edge level:" << std::endl;
    for (eit=cit->edges_begin(); eit!=cit->edges_end(); ++eit) {
      std::cout << eit->curve() << std::endl ;
    }
  }

  return 0;
}
