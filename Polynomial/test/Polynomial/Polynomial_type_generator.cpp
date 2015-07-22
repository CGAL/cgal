#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/use.h>

int main(){
  CGAL_assertion_code(typedef CGAL::Polynomial<int>        Poly_int_1;)
  CGAL_assertion_code(typedef CGAL::Polynomial<Poly_int_1> Poly_int_2;)
  CGAL_assertion_code(typedef CGAL::Polynomial<Poly_int_2> Poly_int_3;)
  
  {
    typedef CGAL::Polynomial_type_generator<int,1>::Type Polynomial;
    CGAL_USE_TYPE(Polynomial);
    CGAL_static_assertion((::boost::is_same<Polynomial, Poly_int_1>::value)); 
  } 
  {
    typedef CGAL::Polynomial_type_generator<int,2>::Type Polynomial;
    CGAL_USE_TYPE(Polynomial);
    CGAL_static_assertion((::boost::is_same<Polynomial, Poly_int_2>::value)); 
  } 
  {
    typedef CGAL::Polynomial_type_generator<int,3>::Type Polynomial;
    CGAL_USE_TYPE(Polynomial);
    CGAL_static_assertion((::boost::is_same<Polynomial, Poly_int_3>::value)); 
  }
}
