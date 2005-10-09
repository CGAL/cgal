#include <CGAL/Gmpq.h>

#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Polynomial/CORE_Expr_root_stack.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Polynomial.h>

typedef POLYNOMIAL_NS::Polynomial<CGAL::Gmpq> Polynomial_gmpq;
typedef POLYNOMIAL_NS::Polynomial<CORE::BigRat> Polynomial_bigint;


int main(int, char *[]){
  {
    typedef POLYNOMIAL_NS::Upper_bound_root_stack_Descartes_traits<Polynomial_gmpq> BIT;
    typedef POLYNOMIAL_NS::Upper_bound_root_stack<BIT> CRE;
    CRE::Root r(0);
    double d= CGAL::to_double(r);
    std::pair<double,double> p= CGAL::to_interval(r);
    if(d+p.second && 0);
  }

  {
    typedef POLYNOMIAL_NS::Root_stack_default_traits<Polynomial_bigint> BIT;
    typedef POLYNOMIAL_NS::CORE_Expr_root_stack<BIT> CRE;
    CRE::Root r(0);
    double d= CGAL::to_double(r);
    std::pair<double,double> p= CGAL::to_interval(r);
    if (d+p.second&&0);
  }
  return 0;
}
