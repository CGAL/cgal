#include <CGAL/Cartesian.h>
#include <CGAL/Circular_kernel.h>
#include <CGAL/Algebraic_kernel_2_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Filtered_interval_circular_kernel.h>
#include <CGAL/intersections.h>
#include <iostream>

typedef CGAL::MP_Float RT;
typedef CGAL::Quotient<RT> NT1;
typedef CGAL::Cartesian<NT1> Linear_k1;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1> Algebraic_k1;
typedef CGAL::Circular_kernel_2<Linear_k1, Algebraic_k1> CircularKernel;
typedef CGAL::Filtered_interval_circular_kernel<CircularKernel> CK;

  CK ck;

#include <CGAL/_test_circles_predicates.h>
#include <CGAL/_test_circles_constructions.h>
#include <CGAL/_test_circles_extention.h>
  
int main() {

  _test_circle_predicat(ck);
  _test_circle_construct(ck);
  _test_circle_bbox(ck);
  _test_circular_arc_bbox(ck);
  _test_has_on(ck);

  return 0;
}
