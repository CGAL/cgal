#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>

#include <iostream>
#include <vector>

int main(int, char**)
{
  typedef CORE::Expr                                                          NT;
  typedef CGAL::Cartesian<NT>                                                 Kernel;
  typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel,
                                        CGAL::Hyperbolic_octagon_translation> Traits;
  typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>        Triangulation;

  Triangulation tr;
  assert(tr.is_valid());
  if(!tr.is_valid())
    return EXIT_FAILURE;

  std::cout << "triangulation works!" << std::endl;
  std::cout << "nb of vertices: " << tr.number_of_vertices() << std::endl;
  std::cout << "nb of faces: " << tr.number_of_faces() << std::endl;

  return EXIT_SUCCESS;
}
