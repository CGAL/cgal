//examples/Arrangement_2/example9.C

// Define shorter names to please linker (g++/egcs)
#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Seg_traits;
typedef CGAL::Arr_polyline_traits_2<Seg_traits>         Traits;

typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;

typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits>                Arr_2;

int main() 
{
  Arr_2                arr;
  std::vector<Point_2> pts;

  // We insert two intersecting squares
  pts.push_back(Point_2(  0, 50));
  pts.push_back(Point_2( 50, 50));
  pts.push_back(Point_2( 50,  0));
  pts.push_back(Point_2(  0,  0));
  pts.push_back(Point_2(  0, 50));
  arr.insert (Curve_2(pts.begin(), pts.end())); 
  
  pts.clear();
  pts.push_back(Point_2( 25, 75));
  pts.push_back(Point_2( 75, 75));
  pts.push_back(Point_2( 75, 25));
  pts.push_back(Point_2( 25, 25));
  pts.push_back(Point_2( 25, 75));
  arr.insert (Curve_2(pts.begin(), pts.end())); 
  
  // Upward vertical ray shooting
  // the edge <25,50>-<50,50> is supposed to be have been created
  Arr_2::Locate_type lt;
  Arr_2::Halfedge_handle e=arr.vertical_ray_shoot(Point_2(30,30),lt,true);
  
  CGAL_assertion(e->source()->point()==Point_2(50,50)); 
  CGAL_assertion(e->target()->point()==Point_2(25,50));
   
  // We expect <50,25> to be an intersection (of the polylines)
  e = arr.locate(Point_2(50,25), lt);
  CGAL_assertion(lt == Arr_2::VERTEX);
  
  return 0;
}
