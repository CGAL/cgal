//! \file examples/Envelope_3/ex_envelope_triangles.cpp
// Constructing the lower and the upper envelope of a set of triangles.

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Env_triangle_traits_3.h>
#include <CGAL/Env_surface_data_traits_3.h>
#include <CGAL/envelope_3.h>
#include <iostream>
#include <list>

typedef CGAL::Exact_rational                             Number_type;
typedef CGAL::Cartesian<Number_type>                     Kernel;
typedef CGAL::Env_triangle_traits_3<Kernel>              Traits_3;
typedef Kernel::Point_3                                  Point_3;
typedef Traits_3::Surface_3                              Triangle_3;
typedef CGAL::Env_surface_data_traits_3<Traits_3, char>  Data_traits_3;
typedef Data_traits_3::Surface_3                         Data_triangle_3;
typedef CGAL::Envelope_diagram_2<Data_traits_3>          Envelope_diagram_2;

/* Auxiliary function - print the features of the given envelope diagram. */
void print_diagram (const Envelope_diagram_2& diag)
{
  // Go over all arrangement faces.
  Envelope_diagram_2::Face_const_iterator            fit;
  Envelope_diagram_2::Ccb_halfedge_const_circulator  ccb;
  Envelope_diagram_2::Surface_const_iterator         sit;

  for (fit = diag.faces_begin(); fit != diag.faces_end(); ++fit)
  {
    // Print the face boundary.
    if (fit->is_unbounded())
    {
      std::cout << "[Unbounded face]";
    }
    else
    {
      // Print the vertices along the outer boundary of the face.
      ccb = fit->outer_ccb();
      std::cout << "[Face]  ";
      do
      {
        std::cout << '(' << ccb->target()->point() << ")  ";
        ++ccb;
      } while (ccb != fit->outer_ccb());
    }

    // Print the labels of the triangles that induce the envelope on this face.
    std::cout << "-->  " << fit->number_of_surfaces()
              << " triangles:";

    for (sit = fit->surfaces_begin(); sit != fit->surfaces_end(); ++sit)
      std::cout << ' ' << sit->data();
    std::cout << std::endl;
  }

  // Go over all arrangement edges.
  Envelope_diagram_2::Edge_const_iterator            eit;

  for (eit = diag.edges_begin(); eit != diag.edges_end(); ++eit)
  {
    // Print the labels of the triangles that induce the envelope on this edge.
    std::cout << "[Edge]  (" << eit->source()->point()
              << ")  (" << eit->target()->point()
              << ")  -->  " << eit->number_of_surfaces()
              << " triangles:";

    for (sit = eit->surfaces_begin(); sit != eit->surfaces_end(); ++sit)
      std::cout << ' ' << sit->data();
    std::cout << std::endl;
  }

  return;
}

/* The main program: */
int main ()
{
  // Construct the input triangles, makred A and B.
  std::list<Data_triangle_3>   triangles;

  triangles.push_back (Data_triangle_3 (Triangle_3 (Point_3 (0, 0, 0),
                                                    Point_3 (0, 6, 0),
                                                    Point_3 (5, 3, 4)),
                                        'A'));
  triangles.push_back (Data_triangle_3 (Triangle_3 (Point_3 (6, 0, 0),
                                                    Point_3 (6, 6, 0),
                                                    Point_3 (1, 3, 4)),
                                        'B'));
  
  // Compute and print the minimization diagram.
  Envelope_diagram_2      min_diag;

  CGAL::lower_envelope_3 (triangles.begin(), triangles.end(),
                          min_diag);

  std::cout << std::endl << "The minimization diagram:" << std::endl;
  print_diagram (min_diag);

  // Compute and print the maximization diagram.
  Envelope_diagram_2      max_diag;

  CGAL::upper_envelope_3 (triangles.begin(), triangles.end(),
                          max_diag);

  std::cout << std::endl << "The maximization diagram:" << std::endl;
  print_diagram (max_diag);

  return (0);
}
