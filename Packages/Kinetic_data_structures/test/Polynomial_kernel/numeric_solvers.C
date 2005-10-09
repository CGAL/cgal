#include <iostream>
#include <cstdlib>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/Numeric_root_stack.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>

#include "Check_solver.h"


bool verbose=true;
typedef POLYNOMIAL_NS::Polynomial<double> Pd;
typedef POLYNOMIAL_NS::Root_stack_default_traits<Pd> Dt;


template <class Solver_traits>
class Nongsl_root_stack {
  typedef Nongsl_root_stack<Solver_traits> This;

  template <class Fn>
  void initialize(const Fn &f, double lb, double ub) {
    std::vector<double> c(f.degree()+1);
    for (int i=0; i<= f.degree(); ++i){
      c[i]= CGAL::to_double(f[i]);
    }
    POLYNOMIAL_NS::internal::jama_polynomial_compute_roots(&*c.begin(),
							   &*c.begin()+ f.degree()+1, 
							   lb, ub, roots_);
  }
  /*void initialize(const POLYNOMIAL_NS::Polynomial<double> &f, double lb, double ub) {
    polynomial_compute_roots(&*f.begin(), &*f.begin()+ f.degree()+1, lb, ub, CLEAN, roots_);
    }*/
public:
  typedef double Root;
  typedef typename Solver_traits::Function Function;
  typedef Solver_traits Traits;
  Nongsl_root_stack(const typename Solver_traits::Function &f, Root lb, Root ub, const Solver_traits&){
    initialize(f, lb, ub);
  }
  
  Nongsl_root_stack(){};

  void pop() {
    Polynomial_precondition(!roots_.empty());
    roots_.pop_back();
  }

  const Root& top() const {
    Polynomial_precondition(!roots_.empty());
    return roots_.back();
  }

  bool empty() const {
    return roots_.empty();
  }

protected:
  std::vector<double> roots_;
};



int main(int argc, char* argv[])
{
  if ( argc > 1 ) {
    int is_verbose = atoi(argv[1]);
    if ( is_verbose == 0 ) {
      verbose = false;
    } else verbose = true;
  }



  {
    if (verbose) std::cout <<"GSL________________________________________\n";
    else std::cout << "GSL &\t";
    typedef POLYNOMIAL_NS::Numeric_root_stack<Dt> NRE;
    typedef POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.all();
    std::cout << std::endl;
  }

  {
    if (verbose) std::cout <<"NonGSL_____________________________________\n";
    else std::cout << "NonGSL &\t";
    typedef Nongsl_root_stack<Dt> NRE;
    typedef POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.all();
    std::cout << std::endl;
  }



  /*{
    if (verbose) std::cout <<"Descartes__________________________________\n";
    else std::cout << "Descartes &\t";
    typedef POLYNOMIAL_NS::Upper_bound_enumerator_Descartes_traits<Pd> Dt;
    typedef POLYNOMIAL_NS::Upper_bound_root_enumerator<Dt> NRE;
    typedef POLYNOMIAL_NS::Kernel<Pd, NRE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.all();
    std::cout << std::endl;
    }*/


  



  

  return 0;
}
