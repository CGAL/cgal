//! \file examples/Example_2/ex_envelope_circles.cpp
// Constructing the envelopes of a set of circles using the circle-segment
// traits.

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Envelope_diagram_1.h>
#include <CGAL/envelope_2.h>

typedef CGAL::Gmpq                                    Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef Kernel::Point_2                               Kernel_point_2;
typedef Kernel::Circle_2                              Circle_2;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>     Traits_2;
typedef Traits_2::Curve_2                             Curve_2;
typedef CGAL::Envelope_diagram_1<Traits_2>            Diagram_1;

/*! Print the given envelope diagram. */
void print_diagram (const Diagram_1& diag)
{
  Diagram_1::Edge_const_handle     e = diag.leftmost();
  Diagram_1::Vertex_const_handle   v;

  while (e != diag.rightmost())
  {
    std::cout << "Edge: ";
    if (! e->is_empty())
    {
      Circle_2      circ = e->curve().supporting_circle();
      std::cout << " (x - " << CGAL::to_double(circ.center().x()) << ")^2 +"
                << " (y - " << CGAL::to_double(circ.center().y()) << ")^2 = "
                << CGAL::to_double(circ.squared_radius()) << std::endl;
    }
    else
      std::cout << " [empty]" << std::endl;

    v = e->right();
    std::cout << "Vertex (" << CGAL::to_double(v->point().x()) << ' '
              << CGAL::to_double(v->point().y()) << ')' << std::endl;

    e = v->right();
  }
  CGAL_assertion (e->is_empty());
  std::cout << "Edge: [empty]" << std::endl;

  return;
}

/*! The main program. */
int main ()
{
  // Create four input circles.
  Curve_2          circles[4];

  circles[0] = Circle_2 (Kernel_point_2 (1, 3), CGAL::square(2));
  circles[1] = Circle_2 (Kernel_point_2 (4, 5), CGAL::square(4));
  circles[2] = Circle_2 (Kernel_point_2 (5, 1), CGAL::square(1));
  circles[3] = Circle_2 (Kernel_point_2 (6, 7), CGAL::square(2));

  // Compute the minimization diagram that represents their lower envelope.
  Diagram_1              min_diag;

  lower_envelope_2 (&(circles[0]), &(circles[4]), min_diag);
  print_diagram (min_diag);

  // Compute the maximization diagram that represents the upper envelope.
  Diagram_1              max_diag;

  upper_envelope_2 (&(circles[0]), &(circles[4]), max_diag);
  print_diagram (max_diag);

  return (0);
}
