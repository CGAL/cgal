#include <CGAL/basic.h>


#include <CGAL/Polynomial/Default_filtering_traits.h>
#include <CGAL/Polynomial/Filtered_kernel.h>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/Numeric_root_stack.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Sturm_root_stack.h>
#include <CGAL/Polynomial/Sturm_root_stack_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack_filtered_Descartes_traits.h>
#include <CGAL/Polynomial/internal/Filtered_function.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>
#include <CGAL/Polynomial/polynomial_converters.h>

#ifdef CGAL_POLYNOMIAL_USE_GSL
#include <CGAL/Polynomial/internal/GSL_numeric_solver.h>
#endif

typedef CGAL_POLYNOMIAL_NS::Polynomial<double> Pd;
typedef CGAL_POLYNOMIAL_NS::Root_stack_default_traits<Pd> Ddt;
typedef CGAL_POLYNOMIAL_NS::Polynomial<double> Polynomial_double;
typedef CGAL_POLYNOMIAL_NS::Polynomial<CGAL::Gmpq> Polynomial_gmpq;


template <class K, class P> 
void solve(K k, P p, const char *name){
  typename CGAL::POLYNOMIAL::Polynomial_converter<P, typename K::Function> pc;
  typename K::Function f= pc(p);
  typename K::Root_container rc= k.root_container_object(f);
  std::cout << name <<"(" << f <<  "): ";
  for (typename K::Root_container::iterator it= rc.begin(); it != rc.end(); ++it){
    std::cout << *it << " ";
  }
  std::cout << std::endl << std::endl;
}


int main(int argc, char *argv[])
{
  while (true) {
    typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Ddt,
      CGAL_POLYNOMIAL_NS::internal::Turkowski_numeric_solver> NRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    K::Function input;
    std::cout << ">> " << std::flush;
    std::cin >> input;
    if (!std::cin) break;
    
    
    solve(k, input, "Turkowski");


#ifdef CGAL_POLYNOMIAL_USE_GSL
    {
      typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Ddt, CGAL_POLYNOMIAL_NS::internal::GSL_numeric_solver> NRE;
      typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
      K k;
      solve(k, input, "GSL");
    }
#endif
#ifdef CGAL_USE_TNT
    {
      typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Ddt,
	CGAL_POLYNOMIAL_NS::internal::JAMA_numeric_solver> NRE;
      typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
      K k;
      solve(k, input, "JAMA");
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
	solve(k, input, "CleanTNT");
    }
#endif
    {
      typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Ddt,
	CGAL_POLYNOMIAL_NS::internal::Turkowski_cleaned_numeric_solver> NRE;
      typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
      K k;
      solve(k, input, "CleanTurk");
    }
    
    {
      typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_Descartes_traits<Polynomial_gmpq> BIT;
      typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<BIT> CRE;
      typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_gmpq, CRE> K;
      
      K k;
      solve(k, input, "Descartes");
    }
    {
      typedef CGAL_POLYNOMIAL_NS::Default_filtering_traits<CGAL::Gmpq> FT;
      typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_filtered_Descartes_traits<FT> DT;
      typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<DT> RE;
      typedef CGAL_POLYNOMIAL_NS::Filtered_kernel<FT, RE> K;
      K k;
      solve(k, input, "DescartesFiltered");
    }
    {
      typedef CGAL_POLYNOMIAL_NS::Default_filtering_traits<CGAL::Gmpq> FT;
      typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_filtered_Descartes_traits<FT> DT;
      typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<DT> RE;
      typedef CGAL_POLYNOMIAL_NS::Filtered_kernel<FT, RE> K;
      K k;
      solve(k, input, "SturmFiltered");
    }

    {
      typedef CGAL_POLYNOMIAL_NS::Sturm_root_stack_traits<Polynomial_gmpq> RET;
      typedef CGAL_POLYNOMIAL_NS::Sturm_root_stack<RET> RE;
      typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_gmpq, RE> K;
      K k;
      solve(k, input, "Sturm");
    }

  };

    return EXIT_SUCCESS;
};
