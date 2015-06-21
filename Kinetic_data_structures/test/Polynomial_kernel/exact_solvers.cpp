#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE

#include <iostream>
#include <cstdlib>

#include <CGAL/Polynomial/Kernel.h>
//#include <CGAL/Polynomial/Filtered_kernel.h>
#include <CGAL/Polynomial/Polynomial.h>
//#include <CGAL/Polynomial/internal/Filtered_function.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
//#include <CGAL/Polynomial/Upper_bound_root_stack.h>
//#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
//#include <CGAL/Polynomial/Upper_bound_root_stack_filtered_Descartes_traits.h>
#include <CGAL/Polynomial/Sturm_root_stack.h>
//#include <CGAL/Polynomial/Default_filtering_traits.h>
#include <CGAL/Polynomial/Sturm_root_stack_traits.h>

#if 0
#include <CGAL/Polynomial/Lazy_upper_bound_root_stack.h>
#endif

#ifdef CGAL_USE_CORE
#include <CGAL/Polynomial/CORE_kernel.h>
#include <CGAL/Polynomial/CORE_Expr_root_stack.h>


#endif

#include "Check_solver.h"

typedef CGAL_POLYNOMIAL_NS::Polynomial<double> Polynomial_double;
typedef CGAL_POLYNOMIAL_NS::Polynomial<CGAL::POLYNOMIAL::Default_field_nt> Polynomial_ft;

bool verbose=true;

int main(int argc, char* argv[])
{
    if ( argc > 1 ) {
        int is_verbose = atoi(argv[1]);
        if ( is_verbose == 0 ) {
            verbose = false;
        } else verbose = true;
    }

    std::cout << "Begin." << std::endl;

#if 0
    {
      CORE::BigFloat cs[9];
      cs[0]=CORE::BigFloat(-2295485086.0);
      cs[1]=CORE::BigFloat(2072822157.0);
      cs[2]=CORE::BigFloat(116461914.2);
      cs[3]=CORE::BigFloat(-116175036.500);
      cs[4]=CORE::BigFloat(-10063149.8700);
      cs[5]=CORE::BigFloat(-196007.034400);
      cs[6]=CORE::BigFloat(3460.88600000);
      cs[7]=CORE::BigFloat(136.910039600);
      cs[8]=CORE::BigFloat(1.0);

      CORE::Polynomial<CORE::BigFloat> p(8, cs);
      //std::cout << p << std::endl;
      Polynomial<CORE::BigFloat> temp(p);
      //std::cout << temp << std::endl;
      Polynomial<CORE::BigFloat> pp=p;
      pp.differentiate(); 
      //std::cout << pp << std::endl;
      Polynomial<CORE::BigFloat> pg = gcd(p, temp.differentiate());

      CORE::BigFloat c;
      Polynomial<CORE::BigFloat> prem=p;
      Polynomial<CORE::BigFloat> pquo= prem.pseudoRemainder(pp, c);
      std::cout << "quo: " << pquo << std::endl;

      //std::cout << R << std::endl;
      CGAL::POLYNOMIAL::internal::CORE_polynomial cp(p);
      CGAL::POLYNOMIAL::internal::CORE_polynomial cpp(pp);
      CGAL::POLYNOMIAL::internal::CORE_polynomial cpg(pg);
      CGAL::POLYNOMIAL::internal::CORE_polynomial cprem(prem);
      CGAL::POLYNOMIAL::internal::CORE_polynomial cpquo(pquo);
      std::cout << "P: " << cp << std::endl;
      std::cout << "P': " << cpp<< std::endl;
      std::cout << "gcd: " << cpg<< std::endl;
      std::cout << "Rem: " << cprem << std::endl;
      std::cout << "Quo: " << cpquo << std::endl;
      std::cout << "C: " << c << std::endl;
    }
#endif

#if 0
    {
        if (verbose) std::cout <<"Descartes_lazy______________________________\n";
        else std::cout << "DL&\t";
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_Descartes_traits<Polynomial_ft> BIT;
        typedef CGAL_POLYNOMIAL_NS::Lazy_upper_bound_root_stack<BIT> CRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_ft, CRE> K;

        K k;
        Check_solver<K> cc(k,verbose);
        cc.all();
        cc.exact();
        if (!verbose) std::cout << " -- &";
        std::cout << std::endl;
    }
#endif

#if 0
    
    {
        if (verbose) std::cout <<"Descartes_exact_____________________________\n";
        else std::cout << "Descartes&\t";
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_Descartes_traits<Polynomial_ft> BIT;
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<BIT> CRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_ft, CRE> K;

        K k;
        Check_solver<K, True> cc(k,verbose);
        cc.all();
        cc.exact();
        //if (!verbose) std::cout << " -- &";
        std::cout << std::endl;
    }
    {
        if (verbose) std::cout <<"Descartes_filtered__________________________\n";
        else std::cout << "Descartes filtered&\t";
        typedef CGAL_POLYNOMIAL_NS::Default_filtering_traits<CGAL::POLYNOMIAL::Default_field_nt> FT;
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_filtered_Descartes_traits<FT> DT;
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<DT> RE;
        typedef CGAL_POLYNOMIAL_NS::Filtered_kernel<FT, RE> K;
        K k;
        Check_solver<K > cg(k, verbose);
        cg.all();
        cg.exact();
        std::cout << std::endl;
    }
    /*{
        if (verbose) std::cout <<"Sturm_filtered__________________________\n";
        else std::cout << "Sturm filtered&\t";
        typedef CGAL_POLYNOMIAL_NS::Default_filtering_traits<CGAL::Ft> FT;
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_filtered_Descartes_traits<FT> DT;
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<DT> RE;
        typedef CGAL_POLYNOMIAL_NS::Filtered_kernel<FT, RE> K;
        K k;
        Check_solver<K > cg(k, verbose);
        cg.all();
        cg.exact();
        std::cout << std::endl;
	}*/
#endif
#ifdef CGAL_USE_GMP
    {
        if (verbose) std::cout <<"Sturm_exact_________________________________\n";
        else std::cout << "Sturm&\t";
        typedef CGAL_POLYNOMIAL_NS::Sturm_root_stack_traits<Polynomial_ft> RET;
        typedef CGAL_POLYNOMIAL_NS::Sturm_root_stack<RET> RE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_ft, RE> K;
        K k;
        Check_solver<K > cg(k,verbose);
        cg.all();
        cg.exact();
        std::cout << std::endl;
	//K::Root rt;
	//std::cout << CGAL::to_interval(rt).first << std::endl;
    }
#endif

#ifdef CGAL_USE_CORE
    {
      
      if (verbose) std::cout << "CORE_______________________________________\n";
      else std::cout << "CORE & ";
      typedef CGAL_POLYNOMIAL_NS::CORE_kernel K;
      K k;
      Check_solver<K > cc(k,verbose);
      cc.all();
      cc.square_free();
      
      cc.mignotte();
      cc.small_intervals();
      cc.non_simple();
      
      //cc.wilkinson();
      
      //cc.wilkinson();
      //cc.mignotte();
      //cc.small_
      cc.exact();
      if (!verbose) std::cout << " -- &";
      std::cout << std::endl;
    }
#endif



#if 0
    {
        if (verbose) std::cout <<"Bezier_exact_________________________________\n";
        else std::cout << "Bezier&\t";
        Check_solver<Bezier_exact_kernel > cg;
        cg.all();
        std::cout << std::endl;
    }
#endif

    return EXIT_SUCCESS;
}
