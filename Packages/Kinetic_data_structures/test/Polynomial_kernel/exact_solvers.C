#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE

#include <iostream>
#include <cstdlib>

#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/Filtered_kernel.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/internal/Filtered_function.h>
#include <CGAL/Polynomial/CORE_Expr_root_stack.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack_filtered_Descartes_traits.h>
#include <CGAL/Polynomial/Sturm_root_stack.h>
#include <CGAL/Polynomial/Default_filtering_traits.h>
#include <CGAL/Polynomial/Sturm_root_stack_traits.h>

#include <CGAL/Gmpq.h>
#include <CORE/BigInt.h>

#include "Check_solver.h"

typedef CGAL_POLYNOMIAL_NS::Polynomial<double> Polynomial_double;
typedef CGAL_POLYNOMIAL_NS::Polynomial<CGAL::Gmpq> Polynomial_gmpq;
typedef CGAL_POLYNOMIAL_NS::Polynomial<CORE::BigRat> Polynomial_bigint;

bool verbose=true;


int main(int argc, char* argv[])
{
  if ( argc > 1 ) {
    int is_verbose = atoi(argv[1]);
    if ( is_verbose == 0 ) {
      verbose = false;
    } else verbose = true;
  }


  {
    if (verbose) std::cout <<"Descartes_filtered__________________________\n";
    else std::cout << "Descartes filtered&\t";
    typedef CGAL_POLYNOMIAL_NS::Default_filtering_traits<CGAL::Gmpq> FT;
    typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_filtered_Descartes_traits<FT> DT;
    typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<DT> RE;
    typedef CGAL_POLYNOMIAL_NS::Filtered_kernel<FT, RE> K;
    K k;
    Check_solver<K > cg(k, verbose);
    cg.all();
    cg.exact();
    std::cout << std::endl;
  }
 {
    if (verbose) std::cout <<"Sturm_filtered__________________________\n";
    else std::cout << "Sturm filtered&\t";
    typedef CGAL_POLYNOMIAL_NS::Default_filtering_traits<CGAL::Gmpq> FT;
    typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_filtered_Descartes_traits<FT> DT;
    typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<DT> RE;
    typedef CGAL_POLYNOMIAL_NS::Filtered_kernel<FT, RE> K;
    K k;
    Check_solver<K > cg(k, verbose);
    cg.all();
    cg.exact();
    std::cout << std::endl;
 }

  {
    if (verbose) std::cout <<"Descartes_exact_____________________________\n";
    else std::cout << "Descartes&\t";
    typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_Descartes_traits<Polynomial_gmpq> BIT;
    typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<BIT> CRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_gmpq, CRE> K;

    K k;
    Check_solver<K> cc(k,verbose);
    cc.all();
    cc.exact();
    if (!verbose) std::cout << " -- &";
    std::cout << std::endl;
  }

  {
    if (verbose) std::cout <<"Sturm_exact_________________________________\n";
    else std::cout << "Sturm&\t";
    typedef CGAL_POLYNOMIAL_NS::Sturm_root_stack_traits<Polynomial_gmpq> RET;
    typedef CGAL_POLYNOMIAL_NS::Sturm_root_stack<RET> RE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_gmpq, RE> K;
    K k;
    Check_solver<K > cg(k,verbose);
    cg.all();
    cg.exact();
    std::cout << std::endl;
  }
  

  /*{
    if (verbose) std::cout << "Filtered_Sturm( Field, Gmpq )______________________\n";
    else std::cout << "F S Gmpq & ";
    typedef CGAL_POLYNOMIAL_NS::Sturm_root_stack_traits<Polynomial_gmpq> RET;
    typedef CGAL_POLYNOMIAL_NS::Sturm_lazy_solver<RET> RE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_gmpq, RE> K;
    Check_solver<Filtered_Sturm_solver_f_gmpq > cs;
    cs.all();
    //cs.square_free();
    //cs.wilkinson();
    //cs.mignotte();
    //cs.small_intervals();
    //cs.non_simple();
    std::cout << std::endl;
    }*/
 



 /*{
    if (verbose) std::cout <<"Filtered Descartes__________________________\n";
    else std::cout << "Descartes filtered&\t";
    Check_solver<Filtered_Descartes_filtered_kernel > cg;
    cg.all();
    std::cout << std::endl;
    }*/

  
 {
    if (verbose) std::cout << "CORE_______________________________________\n";
    else std::cout << "CORE & ";
    typedef CGAL_POLYNOMIAL_NS::Root_stack_default_traits<Polynomial_bigint> BIT;
    typedef CGAL_POLYNOMIAL_NS::CORE_Expr_root_stack<BIT> CRE;
    typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_bigint, CRE> K;
    K k;
    Check_solver<K > cc(k,verbose);
    cc.square_free();
    //cc.wilkinson();
    //cc.mignotte();
    cc.small_intervals();
    //cc.exact();
    if (!verbose) std::cout << " -- &";
    std::cout << std::endl;
  }

#if 0
  {
   if (verbose) std::cout <<"Bezier_exact_________________________________\n";
    else std::cout << "Bezier&\t";
    Check_solver<Bezier_exact_kernel > cg;
    cg.all();
    std::cout << std::endl;
  }
#endif

  
  /*{
    if (verbose) std::cout <<"Descartes_filtered___________________________\n";
    else std::cout << "Descartes filtered&\t";
    Check_solver<Descartes_filtered > cg;
    cg.all();
    std::cout << std::endl;
    }*/

  /*{
    if (verbose) std::cout <<"Descartes_exact_fi___________________________\n";
    else std::cout << "Descates exact filtered_interval&\t";
    Check_solver<Descartes_exact_filtered_interval > cg;
    cg.all();
    std::cout << std::endl;
    }*/

  
  /*{
    if (verbose) std::cout <<"Filtered_descartes_filtered___________________\n";
    else std::cout << "Filtered Descartes filtered&\t";
    Check_solver<Filtered_Descartes_filtered > cg;
    cg.all();
    std::cout << std::endl;
    }*/
  
  /*{
    if (verbose) std::cout <<"Filtered_Descartes_exact_____________________\n";
    else std::cout << "Filtered Descartes exact&\t";
    Check_solver<Filtered_Descartes_exact > cg;
    cg.all();
    std::cout << std::endl;
    }*/

  /*{
    if (verbose) std::cout <<"Bezier_exact_________________________________\n";
    else std::cout << "Bezier exact&\t";
    Check_solver<Bezier_exact > cg;
    cg.all();
    std::cout << std::endl;
    }*/

  // too slow
  /*std::cout << "Sturm( Ring )_____________________________\n";
    Check_solver<Sturm_solver_ring_i > csr;
    csr.all();*/

  
  return EXIT_SUCCESS;
}
