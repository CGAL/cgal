#include <CGAL/basic.h>
#include <CGAL/Exact_circular_kernel.h>
#include <CGAL/intersections.h>
#include <iostream>

typedef CGAL::Exact_circular_kernel_2 CK;

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
