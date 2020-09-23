#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>

template<class AK>
void test_get_arithmetic_kernel(){
  CGAL_assertion_code(typedef typename AK::Integer Integer);
  {
    CGAL_assertion_code(typedef CGAL::Polynomial<Integer> POLY;)
    CGAL_assertion_code(typedef typename CGAL::Get_arithmetic_kernel<POLY>::Arithmetic_kernel AK_;)
    CGAL_static_assertion((boost::is_same<AK,AK_>::value));
  }{
    CGAL_assertion_code(typedef CGAL::Polynomial<CGAL::Polynomial<Integer> > POLY;)
    CGAL_assertion_code(typedef typename CGAL::Get_arithmetic_kernel<POLY>::Arithmetic_kernel AK_;)
    CGAL_static_assertion((boost::is_same<AK,AK_>::value));
  }
}

int main(){
#ifdef CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
   test_get_arithmetic_kernel<CGAL::Arithmetic_kernel>();
#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
    return 0;
}


