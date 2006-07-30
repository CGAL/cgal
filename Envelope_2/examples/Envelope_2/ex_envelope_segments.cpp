#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Env_default_diagram_1.h>
#include <CGAL/envelope_2.h>

#include <list>
#include <iostream>

typedef CGAL::Gmpq                                      NT;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Segment_2;
typedef CGAL::Env_default_diagram_1<Traits_2>           Diagram_1;

int main ()
{
  Traits_2               traits;
  std::list<Segment_2>   segments;
  Diagram_1              diag;

  segments.push_back (Segment_2 (Point_2 (1, 1), Point_2 (6, 3)));
  segments.push_back (Segment_2 (Point_2 (5, 3), Point_2 (9, 1)));
  segments.push_back (Segment_2 (Point_2 (2, 2), Point_2 (5, 2)));
  segments.push_back (Segment_2 (Point_2 (6, 1), Point_2 (8, 4)));
  
  lower_envelope_2 (traits,
                    segments.begin(), segments.end().
                    diag);

  Diagram_1::Edge_const_handle    e = diag.leftmost();
  Diagram_1::Vertex_const_handle  v;

  while (e != diag.rightmost())
  {
    if (! e->is_empty())
      std::cout << e->curve() << "  ";
    else
      std::cout << "[empty]" << "  ";

    v = e->right();
    std::cout << "(" << v->point() << ")  ";

    e = v->right();
  }
  std::cout << "[empty]" << std::endl;

  return (0);
}
