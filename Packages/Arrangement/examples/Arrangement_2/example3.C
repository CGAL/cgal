//example3.C   
#include <CGAL/basic.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_circles_real_traits.h> 
#include <CGAL/Arrangement_2.h>

//#include <CGAL/leda_real.h>
//typedef leda_real                                     NT;

// We use here double instead of leda_real to enable compilation without LEDA.
// This is not recommended generally.
// Read more in the README file or in the manual.
typedef double                                        NT;

typedef CGAL::Arr_circles_real_traits<NT>             Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Arr_base_node<Curve>                     Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>               Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >    Arr_2;

using namespace std;

int main() {
   Arr_2 arr;  

   //2 ccw circles with radius 5 and center (0,0) and (6,0) resp.
   Arr_2::Curve_iterator cit=arr.insert(Curve(0,0,25));
   cit=arr.insert(Curve(6,0,25)); 

   //upward vertical ray shooting
   Arr_2::Locate_type lt;
#ifndef CGAL_NO_ASSERTIONS // in order to avoid warnings
   Arr_2::Halfedge_handle e=
#endif
     arr.vertical_ray_shoot(Point(-1,0),lt,true);

   CGAL_assertion(e->source()->point()==Point(3,4)); 
   CGAL_assertion(e->target()->point()==Point(-5,0));

   return 0;
}
