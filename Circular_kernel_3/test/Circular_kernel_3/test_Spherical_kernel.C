#include <CGAL/Cartesian.h>
#include <CGAL/Spherical_kernel_3.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/_test_sphere_predicates.h>
#include <CGAL/_test_sphere_constructions.h>
#include <CGAL/Polynomials_1_3.h>
#include <CGAL/Polynomials_2_3.h>
#include <CGAL/Polynomials_for_line_3.h>

int main()
{
  typedef CGAL::Quotient< CGAL::MP_Float >                    FT;
  typedef CGAL::Cartesian<FT>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_spheres_2_3<FT>          Algebraic_k1;
  typedef CGAL::Spherical_kernel_3<Linear_k1,Algebraic_k1>    SK1;
  SK1 sk1;
  _test_spherical_kernel_predicates(sk1);
  _test_spherical_kernel_construct(sk1); 
  return 0;
}
