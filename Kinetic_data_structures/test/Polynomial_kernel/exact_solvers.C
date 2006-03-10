//#define CGAL_CHECK_EXACTNESS
//#define CGAL_CHECK_EXPENSIVE

#include <iostream>
#include <cstdlib>

#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/Filtered_kernel.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/internal/Filtered_function.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack_filtered_Descartes_traits.h>
#include <CGAL/Polynomial/Sturm_root_stack.h>
#include <CGAL/Polynomial/Default_filtering_traits.h>
#include <CGAL/Polynomial/Sturm_root_stack_traits.h>

#if 0
#include <CGAL/Polynomial/Lazy_upper_bound_root_stack.h>
#endif

#include <CGAL/Gmpq.h>

#ifdef CGAL_USE_CORE
#include <CGAL/Polynomial/CORE_kernel.h>
#include <CGAL/Polynomial/CORE_Expr_root_stack.h>


#endif

#include "Check_solver.h"

typedef CGAL_POLYNOMIAL_NS::Polynomial<double> Polynomial_double;
typedef CGAL_POLYNOMIAL_NS::Polynomial<CGAL::Gmpq> Polynomial_gmpq;

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

    {
      CORE::BigFloat c[3];
      c[0]= CORE::BigFloat(1.0);
      c[1]= CORE::BigFloat(3.0);
      c[2]= CORE::BigFloat(-4.0);
      CORE::Polynomial<CORE::BigFloat> bf(2, c);
      CORE::BFInterval bfi(CORE::BigFloat(-.5), CORE::BigFloat(-.2));
      std::cout << bf << std::endl;
      CORE::Expr e(bf, bfi);
      std::cout << e << std::endl;
    }

    if (1) {
      /*CORE::Polynomial<CORE::BigFloat> p("-2+20t-50t^2+1t^50", 't');
      CORE::Sturm<CORE::BigFloat> s(p, true);
      std::cout << p << std::endl;
      std::cout << s.numberOfRoots() << std::endl;*/
      CORE::BigRat c[17];
      c[0]= CORE::BigRat("86966370606641/4503599627370496");
      c[1]= CORE::BigRat("9813373613357677/18014398509481984");
      c[2]= CORE::BigRat("3924498795192991/18014398509481984");
      c[3]= CORE::BigRat("-156857429476936177/576460752303423488");
      c[4]= CORE::BigRat("28036759770576529/72057594037927936");
      c[5]= CORE::BigRat("-3478688861042611/9007199254740992");
      c[6]= CORE::BigRat("4395898718207/35184372088832");
      c[7]= CORE::BigRat("-4078546949307093/9007199254740992");
      c[8]= CORE::BigRat("4063208312463087/9007199254740992");
      c[9]= CORE::BigRat("-1076999133570829/4503599627370496");
      c[10]= CORE::BigRat("-21830467801783653/72057594037927936");
      c[11]= CORE::BigRat("813684203309077/4503599627370496");
      c[12]= CORE::BigRat("19160047808947179/36028797018963968");
      c[13]= CORE::BigRat("-22915804375471/140737488355328");
      c[14]= CORE::BigRat("9001478429603223/36028797018963968");
      c[15]= CORE::BigRat("77422038291455461/144115188075855872");
      c[16]= CORE::BigRat("-3556900689495543/72057594037927936");
      CORE::Polynomial<CORE::BigRat> p( 16, c);
      std::cout << p << std::endl;
      CORE::Sturm<CORE::BigRat> s(p);
      std::cout << "done"  << std::endl;
      
    }
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
        if (verbose) std::cout <<"Descartes_lazy______________________________\n";
        else std::cout << "DL&\t";
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_Descartes_traits<Polynomial_gmpq> BIT;
        typedef CGAL_POLYNOMIAL_NS::Lazy_upper_bound_root_stack<BIT> CRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_gmpq, CRE> K;

        K k;
        Check_solver<K> cc(k,verbose);
        cc.all();
        cc.exact();
        if (!verbose) std::cout << " -- &";
        std::cout << std::endl;
    }
#endif


    
    {
        if (verbose) std::cout <<"Descartes_exact_____________________________\n";
        else std::cout << "Descartes&\t";
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack_Descartes_traits<Polynomial_gmpq> BIT;
        typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_stack<BIT> CRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Polynomial_gmpq, CRE> K;

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
    /*{
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
	}*/

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
