#include <iostream>
#include <cstdlib>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/Numeric_root_stack.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>

#ifdef CGAL_POLYNOMIAL_USE_GSL
#include <CGAL/Polynomial/internal/GSL_numeric_solver.h>
#endif

#include "Check_solver.h"


bool verbose=true;
typedef CGAL_POLYNOMIAL_NS::Polynomial<double> Pd;
typedef CGAL_POLYNOMIAL_NS::Root_stack_default_traits<Pd> Dt;



int main(int argc, char* argv[])
{
  if ( argc > 1 ) {
    int is_verbose = atoi(argv[1]);
    if ( is_verbose == 0 ) {
      verbose = false;
    } else verbose = true;
  }


#ifdef CGAL_POLYNOMIAL_USE_GSL
  {
    if (verbose) std::cout <<"GSL________________________________________\n";
    else std::cout << "GSL &\t";
    typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt, CGAL_POLYNOMIAL_NS::internal::GSL_numeric_solver> NRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.all();
    std::cout << std::endl;
  }
#endif 

  {
    if (verbose) std::cout <<"JAMA______________________________________\n";
    else std::cout << "JAMA &\t";
    typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt,
      CGAL_POLYNOMIAL_NS::internal::JAMA_numeric_solver> NRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.all();
    std::cout << std::endl;
  }

  {
    if (verbose) std::cout <<"Turk______________________________________\n";
    else std::cout << "Turk &\t";
    typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt,
      CGAL_POLYNOMIAL_NS::internal::Turkowski_numeric_solver> NRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.all();
    std::cout << std::endl;
  }


#ifdef CGAL_POLYNOMIAL_USE_GSL
  {
    if (verbose) std::cout <<"CleanGSL__________________________________\n";
    else std::cout << "CleanGSL &\t";
    typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt, 
      CGAL_POLYNOMIAL_NS::internal::GSL_cleaned_numeric_solver> NRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.cleaned();
    std::cout << std::endl;
  }
#endif 

  {
    if (verbose) std::cout <<"CleanJAMA________________________________\n";
    else std::cout << "CleanJAMA &\t";
    typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt,
      CGAL_POLYNOMIAL_NS::internal::JAMA_cleaned_numeric_solver> NRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.cleaned();
    std::cout << std::endl;
  }

  {
    if (verbose) std::cout <<"CleanTurk________________________________\n";
    else std::cout << "CleanTurk &\t";
    typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt,
      CGAL_POLYNOMIAL_NS::internal::Turkowski_cleaned_numeric_solver> NRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.cleaned();
    std::cout << std::endl;
  }



  /*{
    if (verbose) std::cout <<"Descartes__________________________________\n";
    else std::cout << "Descartes &\t";
    typedef CGAL_POLYNOMIAL_NS::Upper_bound_enumerator_Descartes_traits<Pd> Dt;
    typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_enumerator<Dt> NRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.all();
    std::cout << std::endl;
    }*/


  



  

  return 0;
}
