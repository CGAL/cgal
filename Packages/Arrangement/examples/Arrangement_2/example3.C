// file: examples/Arrangement_2/example3.C   

#include "short_names.h"

#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
// To enable compilation without core:
int main ()
{
  return (0);
}

#else

#include <CGAL/Cartesian.h>
#include <CORE/BigInt.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_conic_traits_2.h> 
#include <CGAL/Arrangement_2.h>

typedef CORE::BigInt                                CfNT;
typedef CGAL::Cartesian<CfNT>                       Int_kernel;
typedef Int_kernel::Point_2                         Int_point_2;
typedef Int_kernel::Circle_2                        Int_circle_2;

typedef CORE::Expr                                      CoNT;
typedef CGAL::Cartesian<CoNT>                           Alg_kernel;

typedef CGAL::Arr_conic_traits_2<Int_kernel,Alg_kernel> Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arr_2_default_dcel<Traits_2>              Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits_2>              Arr_2;

int main()
{
  Arr_2 arr;  
              
  // 2 ccw circles with radius 5 and center (0,0) and (6,0) resp.
  Int_circle_2  c1 (Int_point_2(0,0), 5*5);
  Int_circle_2  c2 (Int_point_2(6,0), 5*5);

  Arr_2::Curve_iterator cit = arr.insert(Curve_2 (c1));
  cit = arr.insert(Curve_2 (c2)); 

  // upward vertical ray shooting
  Arr_2::Locate_type lt;
  Arr_2::Halfedge_handle e = arr.vertical_ray_shoot(Point_2(-1, 0), lt, true);

  CGAL_assertion(e->source()->point() == Point_2(3, 4)); 
  CGAL_assertion(e->target()->point() == Point_2(-5, 0));

  return 0;
}

#endif
