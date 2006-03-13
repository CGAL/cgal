#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Polynomial.h>

#ifdef CGAL_USE_CORE
#include <CGAL/Polynomial/CORE_Expr_root_stack.h>
#endif

typedef CGAL_POLYNOMIAL_NS::Polynomial<CGAL::POLYNOMIAL::Default_field_nt> Polynomial_ft;


#ifdef CGAL_USE_CORE
typedef CGAL_POLYNOMIAL_NS::Polynomial<CORE::BigRat> Polynomial_bigint;
#endif

int main(int, char *[])
{
    {
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_Descartes_traits<Polynomial_ft> BIT;
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<BIT> CRE;
        CRE::Root r(0);
        double d= CGAL::to_double(r);
        std::pair<double,double> p= CGAL::to_interval(r);
        if((d+p.second>0) && 0);
    }

#ifdef CGAL_USE_CORE
    {
      typedef CGAL_POLYNOMIAL_NS::CORE_Expr_root_stack CRE;
      CRE::Root r(0);
      double d= CGAL::to_double(r);
      std::pair<double,double> p= CGAL::to_interval(r);
      if ((d+p.second>0)&&0);
    }
#endif
    return 0;
}
