//! \file examples/Example_2/ex_envelope_segments.cpp
// Constructing the lower envelope of a set of segments.

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Envelope_diagram_1.h>
#include <CGAL/envelope_2.h>

#include <list>
#include <iostream>

typedef CGAL::Gmpq                                      Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef Segment_traits_2::X_monotone_curve_2            Segment_2;
typedef CGAL::Arr_curve_data_traits_2<Segment_traits_2,
                                      char>             Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::X_monotone_curve_2                    Labeled_segment_2;
typedef CGAL::Envelope_diagram_1<Traits_2>              Diagram_1;

int main ()
{
  // Consrtuct the input segments and label them 'A' ... 'H'.
  std::list<Labeled_segment_2>   segments;

  segments.push_back (Labeled_segment_2 (Segment_2 (Point_2 (0, 1),
                                                    Point_2 (2, 3)), 'A'));
  segments.push_back (Labeled_segment_2 (Segment_2 (Point_2 (1, 2),
                                                    Point_2 (4, 5)), 'B'));
  segments.push_back (Labeled_segment_2 (Segment_2 (Point_2 (1, 5),
                                                    Point_2 (7, 2)), 'C'));
  segments.push_back (Labeled_segment_2 (Segment_2 (Point_2 (4, 2),
                                                    Point_2 (6, 4)), 'D'));
  segments.push_back (Labeled_segment_2 (Segment_2 (Point_2 (8, 3),
                                                    Point_2 (8, 6)), 'E'));
  segments.push_back (Labeled_segment_2 (Segment_2 (Point_2 (9, 2),
                                                    Point_2 (12, 4)), 'F'));
  segments.push_back (Labeled_segment_2 (Segment_2 (Point_2 (10, 2),
                                                    Point_2 (12, 1)), 'G'));
  segments.push_back (Labeled_segment_2 (Segment_2 (Point_2 (11, 0),
                                                    Point_2 (11, 5)), 'H'));

  // Compute the minimization diagram that represents their lower envelope.
  Diagram_1              min_diag;

  lower_envelope_x_monotone_2 (segments.begin(), segments.end(),
                               min_diag);

  // Print the minimization diagram.
  Diagram_1::Edge_const_handle     e = min_diag.leftmost();
  Diagram_1::Vertex_const_handle   v;
  Diagram_1::Curve_const_iterator  cit;

  while (e != min_diag.rightmost())
  {
    std::cout << "Edge:";
    if (! e->is_empty())
    {
      for (cit = e->curves_begin(); cit != e->curves_end(); ++cit)
        std::cout << ' ' << cit->data();
    }
    else
      std::cout << " [empty]";
    std::cout << std::endl;

    v = e->right();
    std::cout << "Vertex (" << v->point() << "):";
    for (cit = v->curves_begin(); cit != v->curves_end(); ++cit)
      std::cout << ' ' << cit->data();
    std::cout << std::endl;

    e = v->right();
  }
  CGAL_assertion (e->is_empty());
  std::cout << "Edge: [empty]" << std::endl;

  return (0);
}
