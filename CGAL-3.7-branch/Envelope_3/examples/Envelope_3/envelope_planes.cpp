//! \file examples/Envelope_3/ex_envelope_planes.cpp
// Constructing the lower and the upper envelope of a set of planes.

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Env_plane_traits_3.h>
#include <CGAL/envelope_3.h>
#include <iostream>
#include <list>

typedef CGAL::Gmpq                                       Number_type;
typedef CGAL::Cartesian<Number_type>                     Kernel;
typedef Kernel::Plane_3                                  Plane_3;
typedef CGAL::Env_plane_traits_3<Kernel>                 Traits_3;
typedef Traits_3::Surface_3                              Surface_3;
typedef CGAL::Envelope_diagram_2<Traits_3>               Envelope_diagram_2;


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

    // Print the vertices along the outer boundary of the face.
    ccb = fit->outer_ccb();
    std::cout << "[Face]  ";
    do
    {
      if(!ccb->is_fictitious())
        std::cout << '(' << ccb->curve() << ") ";
      ++ccb;
    } while (ccb != fit->outer_ccb());

    // Print the planes that induce the envelope on this face.
    std::cout << "-->  " << fit->number_of_surfaces()
              << " planes:";

    for (sit = fit->surfaces_begin(); sit != fit->surfaces_end(); ++sit)
      std::cout << ' ' << sit->plane();
    std::cout << std::endl;
  }

  return;
}

/* The main program: */
int main ()
{
  // Construct the input planes.
  std::list<Surface_3>   planes;

  planes.push_back (Surface_3(Plane_3(0, -1, 1, 0)));
  planes.push_back (Surface_3(Plane_3(-1, 0, 1, 0)));
  planes.push_back (Surface_3(Plane_3(0, 1 , 1, 0)));
  planes.push_back (Surface_3(Plane_3(1, 0, 1,  0)));

  // Compute and print the minimization diagram.
  Envelope_diagram_2      min_diag;

  CGAL::lower_envelope_3 (planes.begin(), planes.end(), min_diag);

  std::cout << std::endl << "The minimization diagram:" << std::endl;
  print_diagram (min_diag);

  // Compute and print the maximization diagram.
  Envelope_diagram_2      max_diag;

  CGAL::upper_envelope_3 (planes.begin(), planes.end(), max_diag);

  std::cout << std::endl << "The maximization diagram:" << std::endl;
  print_diagram (max_diag);

  return (0);
}
