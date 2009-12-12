#include <CGAL/basic.h>


#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/Numeric_root_stack.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Sturm_root_stack.h>
#include <CGAL/Polynomial/Sturm_root_stack_traits.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>
#include <CGAL/Polynomial/polynomial_converters.h>

#ifdef CGAL_POLYNOMIAL_USE_GSL
#include <CGAL/Polynomial/internal/GSL_numeric_solver.h>
#endif

typedef CGAL_POLYNOMIAL_NS::Polynomial<double> Pd;
typedef CGAL_POLYNOMIAL_NS::Root_stack_default_traits<Pd> Ddt;
typedef CGAL_POLYNOMIAL_NS::Polynomial<double> Polynomial_double;
typedef CGAL_POLYNOMIAL_NS::Polynomial<CGAL::POLYNOMIAL::Default_field_nt> Polynomial_ft;


template <class K, class P>
void solve(K k, P p, double lb, double ub, const char *name){
  typename CGAL::POLYNOMIAL::Polynomial_converter<P, typename K::Function> pc;
  typename K::Function f= pc(p);
  typename K::Root rlb(lb), rub(ub);
  typename K::Root_container rc= k.root_container_object(f, rlb, rub);
  std::cout << name <<"(" << f <<  "): ";
  for (typename K::Root_container::iterator it= rc.begin(); it != rc.end(); ++it){
    std::cout << *it << " ";
  }
  std::cout << std::endl << std::endl;
}


int main(int , char **)
{
  std::cout << "Enter a polynomial like 1+2*x+3*x^2\n";
  while (true) {
    typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Ddt,
      CGAL_POLYNOMIAL_NS::internal::Turkowski_numeric_solver> NRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    K::Function input;
    std::cout << "polynomial >> " << std::flush;
    std::cin >> input;
    if (!std::cin) break;
    std::cout << "lower bound on roots >> " << std::flush;
    double lb;
    std::cin >> lb;
    if (!std::cin) break;
    std::cout << "upper bound on roots >> " << std::flush;
    double ub;
    std::cin >> ub;
    if (!std::cin) break;


    solve(k, input, lb, ub, "Turkowski");


#ifdef CGAL_POLYNOMIAL_USE_GSL
    {
      typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Ddt, CGAL_POLYNOMIAL_NS::internal::GSL_numeric_solver> NRE;
      typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
      K k;
      solve(k, input, lb, ub, "GSL");
    }
#endif
#ifdef CGAL_USE_TNT
    {
      typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Ddt,
	CGAL_POLYNOMIAL_NS::internal::JAMA_numeric_solver> NRE;
      typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
      K k;
      solve(k, input, lb, ub, "JAMA");
    }
#endif
#ifdef CGAL_POLYNOMIAL_USE_GSL
    {
        typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Ddt,
            CGAL_POLYNOMIAL_NS::internal::GSL_cleaned_numeric_solver> NRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
        K k;
	solve(k, intput, "CleanGSL");
    }
#endif

#ifdef CGAL_USE_TNT
    {
        typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Ddt,
            CGAL_POLYNOMIAL_NS::internal::JAMA_cleaned_numeric_solver> NRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
        K k;
	solve(k, input, lb, ub, "CleanTNT");
    }
#endif
    {
      typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Ddt,
	CGAL_POLYNOMIAL_NS::internal::Turkowski_cleaned_numeric_solver> NRE;
      typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
      K k;
      solve(k, input, lb, ub, "CleanTurk");
    }

    {
      typedef CGAL_POLYNOMIAL_NS::Sturm_root_stack_traits<Polynomial_ft> RET;
      typedef CGAL_POLYNOMIAL_NS::Sturm_root_stack<RET> RE;
      typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_ft, RE> K;
      K k;
      solve(k, input, lb, ub, "Sturm");
    }

  };

    return EXIT_SUCCESS;
}
