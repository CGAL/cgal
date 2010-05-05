#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>

template<class AT>
void test_get_arithmetic_kernel(){
  typedef CGAL::Arithmetic_kernel AK;
  typedef AK::Integer Integer;
  {
    typedef CGAL::Polynomial<Integer> POLY;
    typedef CGAL::Get_arithmetic_kernel<POLY>::Arithmetic_kernel AK_;
    BOOST_STATIC_ASSERT((boost::is_same<AK,AK_>::value));
  }{
    typedef CGAL::Polynomial<CGAL::Polynomial<Integer> > POLY;
    typedef CGAL::Get_arithmetic_kernel<POLY>::Arithmetic_kernel AK_;
    BOOST_STATIC_ASSERT((boost::is_same<AK,AK_>::value));
  }
}

int main(){ 
#ifdef CGAL_USE_LEDA
   test_get_arithmetic_kernel<CGAL::LEDA_arithmetic_kernel>();
#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_CORE
    test_get_arithmetic_kernel<CGAL::CORE_arithmetic_kernel>();
#endif // CGAL_USE_CORE
    return 0;
}


