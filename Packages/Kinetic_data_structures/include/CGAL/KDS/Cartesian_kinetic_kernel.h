#ifndef KINETIC_KERNEL_H
#define KINETIC_KERNEL_H
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/internal/Kernel/Cartesian_kinetic_kernel_base.h>

//#include <CGAL/KDS_internals/kernel_defines.h>
/*#define FOD(UC, lc, d) typedef internal::Cartesian_##lc##_##d<This> UC##_##d;\
UC##_##d lc##_##d##_object() const {\
return UC##_##d();\
}*/
/*#define FODW(UC, lc, d) typedef internal::Cartesian_##lc##_##d<Weighted_point_##d> UC##_##d;\
UC##_##d lc##_##d##_object() const {\
return UC##_##d();\
}*/

CGAL_KDS_BEGIN_NAMESPACE

//! A kinetic kernel using cartesian coordinates
/*!  It takes a PolynomialKernel as a template parameter. The
  PolynomialKernel is used to define the Motion_function and the
  Certificate_function.
*/
template <class Polynomial_k>
class Cartesian_kinetic_kernel: 
  public internal::Cartesian_kinetic_kernel_base<Polynomial_k, 
						 Cartesian_kinetic_kernel<Polynomial_k> > {
  typedef internal::Cartesian_kinetic_kernel_base<Polynomial_k, 
						  Cartesian_kinetic_kernel<Polynomial_k> > P;
public:
  Cartesian_kinetic_kernel(Polynomial_k pk): P(pk){}
  Cartesian_kinetic_kernel(){}
};

CGAL_KDS_END_NAMESPACE

//#include <CGAL/KDS_internals/kernel_undefs.h>
#endif
