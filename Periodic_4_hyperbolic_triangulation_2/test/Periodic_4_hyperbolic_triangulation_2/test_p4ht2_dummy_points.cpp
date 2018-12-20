#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/internal/Periodic_4_hyperbolic_triangulation_dummy_14.h>
#include <CGAL/Hyperbolic_octagon_translation.h>

#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/determinant.h>

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

typedef CORE::Expr                                                              NT;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT>                              AK;   // Algebraic kernel
typedef CGAL::Cartesian<NT>                                                     BK;   // Basic kernel
typedef CGAL::Circular_kernel_2<BK, AK>                                         Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<
                Kernel, CGAL::Hyperbolic_octagon_translation>                   Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef Triangulation::Face_iterator                                            Face_iterator;

int main(int, char**)
{
  Triangulation tr;

  std::cout << "Triangulation successfully initialized with dummy points!" << std::endl << "---------------------------------------------" << std::endl;
  std::cout << "Number of vertices:                  " << tr.number_of_vertices() << std::endl;
  std::cout << "Number of faces:                     " << tr.number_of_faces() << std::endl;
  std::cout << "Number of edges:                     " << tr.number_of_edges() << std::endl;
  std::cout << "Expected edges (by Euler relation):  " << tr.number_of_vertices() + tr.number_of_faces() + 2 << std::endl;

  assert(tr.is_valid(true));
  if(!tr.is_valid(true))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
