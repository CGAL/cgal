// examples/Arrangement_2/example7.C
// ---------------------------------

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
typedef Traits::Point_2                                 Point;
typedef Traits::Curve_2                                 Curve;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;

// A global variable to keep track of when we inserted the edges.
int global_int = 0;

struct Info 
{ 
  Info() : d(global_int) {} 
  // Initializes the private instance with the global variable
  
  int num() {return d;}
  int d;
};

struct My_base_node : public CGAL::Arr_base_node<Curve, X_monotone_curve_2> 
{
  typedef CGAL::Arr_base_node<Curve, X_monotone_curve_2> Base;
  My_base_node() : Base(), info() {}

  Info info; // Will be initialized with the default ctr.
};


struct My_halfedge : public CGAL::Arr_2_halfedge_base<My_base_node> {
  typedef CGAL::Arr_2_halfedge_base<My_base_node> Base;
  My_halfedge() : Base(), info() {}

  Info info; // Will be initialized with the default ctr.
};

typedef CGAL::Pm_dcel<CGAL::Arr_2_vertex_base<Point>,  
                      My_halfedge, 
                      CGAL::Arr_2_face_base>            Dcel;

typedef CGAL::Arrangement_2<Dcel,Traits,My_base_node>   Arr_2;

int main() 
{
  Arr_2 arr;
  
  arr.set_update(false); // The insertions won't update the planar map.

  // Insertion of the curves
  Arr_2::Curve_iterator cit=arr.insert(Curve(Point(0,0),Point(1,1)));
  cit=arr.insert(Curve(Point(0,1),Point(1,0))); 
  CGAL_assertion(arr.number_of_halfedges()==0);
 
  // Change the global variable so we will see the changes
  // all Info structs will have 1 in them from now on.
  global_int=1;

  arr.set_update(true);
  CGAL_assertion(arr.number_of_halfedges() == 8);

  // Traversal of the curves
  Arr_2::Edge_iterator eit;
  for (cit = arr.curve_node_begin(); cit!=arr.curve_node_end(); ++cit) 
  {
    std::cout << std::endl << "Curve level info:" << std::endl 
              << cit->info.num() << std::endl ;
    std::cout << "Edge level info:" << std::endl;
    for (eit = cit->edges_begin(); eit!=cit->edges_end(); ++eit) 
    {
      std::cout << eit->info.num() << std::endl ;
    }
  }
  return 0;
}
