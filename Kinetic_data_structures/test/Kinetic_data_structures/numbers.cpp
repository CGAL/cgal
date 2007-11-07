#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Kinetic/basic.h>
#include <CGAL/Testsuite/assert.h>
#include <limits>
#include <CGAL/Kinetic/Exact_simulation_traits.h>

template <class NT>
void check_nt(NT v) {
  NT zero(0);
  CGAL_test_assert(v <= zero || v > zero);                  // break symmetry due to compiler warning
  if (std::numeric_limits<NT>::has_infinity) {
    NT inf= std::numeric_limits<NT>::infinity();
    //NT inf2=inf;
    NT ninf= -inf;
    CGAL_test_assert(inf > v);
    CGAL_test_assert(ninf < v);
    CGAL_test_assert(inf > zero);
    CGAL_test_assert(ninf < zero);
    CGAL_test_assert(inf > ninf);
    CGAL_test_assert(ninf < inf);
    //CGAL_test_assert(inf2==inf);
  }
}

int main(int, char *[])
{
  check_nt(0.0);
  check_nt(CGAL::Kinetic::Default_field_nt(1));
  
  typedef CGAL::Kinetic::Exact_simulation_traits Tr;
  Tr tr(0,1000000);
  Tr::Simulator::Function_kernel::Function fn= tr.kinetic_kernel_object().function_kernel_object().construct_function_object()(1,0,-2);
  check_nt(Tr::Simulator::Time(-1));
 
  if (CGAL::Kinetic::internal::audit_failures__ != 0) return EXIT_FAILURE;
  else return EXIT_SUCCESS;
}
