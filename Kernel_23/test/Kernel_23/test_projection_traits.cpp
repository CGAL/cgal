#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Projection_traits_xy_3<Epic> K;

typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;


int main()
{
  Point_2 p3(1,0,0), q3(1,1,0), r3(1,2,0);
  Vector_2 v3(1, 0, 0), w3(0,1,0);

  Epic::Point_2 p2(1,0), q2(1,1), r2(1,2);
  Epic::Vector_2 v2(1, 0), w2(0,1);

  K k;

  assert( k.compute_scalar_product_2_object()(v3, w3) ==
          v2 * w2 );

  assert( k.collinear_2_object()(p3,q3,r3) ==
          CGAL::collinear(p2,q2,r2) );

  assert( k.collinear_are_ordered_along_line_2_object()(p3,q3,r3) ==
          CGAL::collinear_are_ordered_along_line(p2,q2,r2) );

  assert( k.compute_squared_length_2_object()(v3) ==
          v2.squared_length() );

 return 0;
}
