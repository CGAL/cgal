//examples/Arrangement_2/example10.C

// Define shorter names to please linker (g++)
#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>

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

  // Curve #1, not x monotone.
  pts.push_back(Point_2(  0,  0));
  pts.push_back(Point_2( 10, 10));
  pts.push_back(Point_2(  0, 20));
  arr.insert (Curve_2(pts.begin(), pts.end())); 
  
  // Curve #2, x monotone.
  pts.clear();
  pts.push_back(Point_2(100,  0));
  pts.push_back(Point_2(150, 50));
  pts.push_back(Point_2(200,  0));
  arr.insert (Curve_2(pts.begin(), pts.end()));

  // Curve #1 is broken into two edges. Point_2 (10,10) turns into a vertex.
  Arr_2::Locate_type lt;
  arr.locate(Point_2(10, 10), lt);
  CGAL_assertion(lt == Arr_2::VERTEX);

  return 0;
}
