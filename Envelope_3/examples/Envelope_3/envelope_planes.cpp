//! \file examples/Envelope_3/ex_envelope_planes.cpp
// Constructing the lower and the upper envelope of a set of planes.

#include <iostream>
#include <list>

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Env_plane_traits_3.h>
#include <CGAL/envelope_3.h>

using Number_type = CGAL::Exact_rational;
using Kernel = CGAL::Cartesian<Number_type>;
using Plane_3 = Kernel::Plane_3;
using Traits_3 = CGAL::Env_plane_traits_3<Kernel>;
using Surface_3 = Traits_3::Surface_3;
using Envelope_diagram_2 = CGAL::Envelope_diagram_2<Traits_3>;

/* Auxiliary function - print the features of the given envelope diagram. */
void print_diagram(const Envelope_diagram_2& diag) {
  // Go over all arrangement faces.
  for (auto fit = diag.faces_begin(); fit != diag.faces_end(); ++fit) {
    // Print the face boundary.

    // Print the vertices along the outer boundary of the face.
    auto ccb = fit->outer_ccb();
    std::cout << "[Face] ";
    do if (!ccb->is_fictitious()) std::cout << '(' << ccb->curve() << ") ";
    while (++ccb != fit->outer_ccb());

    // Print the planes that induce the envelope on this face.
    std::cout << "--> " << fit->number_of_surfaces() << " planes:";

    for (auto sit = fit->surfaces_begin(); sit != fit->surfaces_end(); ++sit)
      std::cout << ' ' << sit->plane();
    std::cout << std::endl;
  }
}

/* The main program: */
int main() {
  // Construct the input planes.
  std::list<Surface_3> planes;

  planes.push_back(Surface_3(Plane_3(0, -1, 1, 0)));
  planes.push_back(Surface_3(Plane_3(-1, 0, 1, 0)));
  planes.push_back(Surface_3(Plane_3(0, 1 , 1, 0)));
  planes.push_back(Surface_3(Plane_3(1, 0, 1,  0)));

  // Compute and print the minimization diagram.
  Envelope_diagram_2 min_diag;
  CGAL::lower_envelope_3(planes.begin(), planes.end(), min_diag);
  std::cout << std::endl << "The minimization diagram:" << std::endl;
  print_diagram(min_diag);

  // Compute and print the maximization diagram.
  Envelope_diagram_2 max_diag;
  CGAL::upper_envelope_3(planes.begin(), planes.end(), max_diag);
  std::cout << std::endl << "The maximization diagram:" << std::endl;
  print_diagram (max_diag);

  return 0;
}
