// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/basic.h>

#include <vector>
#include <sstream>
#include <cstddef>
#include <cassert>

#include <CGAL/Polynomial/sturm_habicht_sequence.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>


template<typename ArithmeticKernel> void test_routine() {

  // All results are computed with the sres_stha.mpl worksheet
  typedef typename ArithmeticKernel::Integer NT;
  typedef CGAL::Polynomial<NT> PNT_1;
  typedef CGAL::Polynomial_traits_d<PNT_1> PNT_1_traits;
  typedef CGAL::Polynomial<PNT_1> PNT_2;
  typedef CGAL::Polynomial_traits_d<PNT_2> PNT_2_traits;
  
  { // test for the principal Sturm-Habicht sequence
    // f=x^3*(x+1)*(x-1)*(x-2)^2
    PNT_1 f(NT(0),NT(0),NT(0),NT(-4),NT(4),NT(3),NT(-4),NT(1));
    std::vector<NT> stha;
    CGAL::principal_sturm_habicht_sequence
        (f,std::back_inserter(stha));
        
    assert(stha.size()==8);
    assert(stha[0] == NT(0));
    assert(stha[1] == NT(0));
    assert(stha[2] == NT(0));
    assert(stha[3] == NT(864));
    assert(stha[4] == NT(324));
    assert(stha[5] == NT(54));
    assert(stha[6] == NT(7));
    assert(stha[7] == NT(1));
    assert(CGAL::number_of_real_roots
             (stha.begin(),stha.end()) == 4);
    assert(CGAL::number_of_real_roots(f)==4);
    
    stha.clear();
    
    CGAL::internal::prs_principal_sturm_habicht_sequence<PNT_1_traits>
        (f,std::back_inserter(stha));
        
    assert(stha.size()==8);
    assert(stha[0] == NT(0));
    assert(stha[1] == NT(0));
    assert(stha[2] == NT(0));
    assert(stha[3] == NT(864));
    assert(stha[4] == NT(324));
    assert(stha[5] == NT(54));
    assert(stha[6] == NT(7));
    assert(stha[7] == NT(1));
    assert(CGAL::number_of_real_roots
             (stha.begin(),stha.end()) == 4);
    
    stha.clear();
    
    CGAL::internal::bezout_principal_sturm_habicht_sequence<PNT_1_traits>
        (f,std::back_inserter(stha));
        
    assert(stha.size()==8);
    assert(stha[0] == NT(0));
    assert(stha[1] == NT(0));
    assert(stha[2] == NT(0));
    assert(stha[3] == NT(864));
    assert(stha[4] == NT(324));
    assert(stha[5] == NT(54));
    assert(stha[6] == NT(7));
    assert(stha[7] == NT(1));
    assert(CGAL::number_of_real_roots
             (stha.begin(),stha.end()) == 4);

    stha.clear();

    typename 
        PNT_1_traits::Principal_sturm_habicht_sequence()
        (f,std::back_inserter(stha));
        
    assert(stha.size()==8);
    assert(stha[0] == NT(0));
    assert(stha[1] == NT(0));
    assert(stha[2] == NT(0));
    assert(stha[3] == NT(864));
    assert(stha[4] == NT(324));
    assert(stha[5] == NT(54));
    assert(stha[6] == NT(7));
    assert(stha[7] == NT(1));
    assert(CGAL::number_of_real_roots
             (stha.begin(),stha.end()) == 4);
  }
  
  { // Example for the defective case
    // f:=x^8+x^2
    PNT_1 f(NT(0),NT(0),NT(1),NT(0),NT(0),NT(0),NT(0),NT(0),NT(1));

    std::vector<NT> stha;
    CGAL::principal_sturm_habicht_sequence
        (f,std::back_inserter(stha));
        
    assert(stha.size()==9);
    assert(stha[0] == NT(0));
    assert(stha[1] == NT(-93312));
    assert(stha[2] == NT(-62208));
    assert(stha[3] == NT(0));
    assert(stha[4] == NT(0));
    assert(stha[5] == NT(0));
    assert(stha[6] == NT(0));
    assert(stha[7] == NT(8));
    assert(stha[8] == NT(1));
    assert(CGAL::number_of_real_roots
               (stha.begin(),stha.end()) == 1);

    stha.clear();
    
    CGAL::internal::prs_principal_sturm_habicht_sequence<PNT_1_traits> 
        (f,std::back_inserter(stha));
        
    assert(stha.size()==9);
    assert(stha[0] == NT(0));
    assert(stha[1] == NT(-93312));
    assert(stha[2] == NT(-62208));
    assert(stha[3] == NT(0));
    assert(stha[4] == NT(0));
    assert(stha[5] == NT(0));
    assert(stha[6] == NT(0));
    assert(stha[7] == NT(8));
    assert(stha[8] == NT(1));
    assert(CGAL::number_of_real_roots
               (stha.begin(),stha.end()) == 1);

    stha.clear();

    CGAL::internal::bezout_principal_sturm_habicht_sequence<PNT_1_traits> 
        (f,std::back_inserter(stha));
        
    assert(stha.size()==9);
    assert(stha[0] == NT(0));
    assert(stha[1] == NT(-93312));
    assert(stha[2] == NT(-62208));
    assert(stha[3] == NT(0));
    assert(stha[4] == NT(0));
    assert(stha[5] == NT(0));
    assert(stha[6] == NT(0));
    assert(stha[7] == NT(8));
    assert(stha[8] == NT(1));
    assert(CGAL::number_of_real_roots
               (stha.begin(),stha.end()) == 1);
    
  }
 

  { // test for the complete Sturm-Habicht sequence
    // f:=(2*x+1)^3*(x+1)*(x-1)
    PNT_1 f(NT(-1),NT(-6),NT(-11),NT(-2),NT(12),NT(8));
    PNT_1 fx(f);
    fx.diff();

    std::vector<PNT_1> stha;

    CGAL::sturm_habicht_sequence
        (f,std::back_inserter(stha));

    assert(stha.size()==6);

    assert(stha[5]==f);
    assert(stha[4]==fx);
    assert(stha[3]==PNT_1(NT(1024),NT(5568),NT(9984),NT(5888)));
    assert(stha[2]==PNT_1(NT(55296),NT(221184),NT(221184)));
    assert(stha[1].is_zero());
    assert(stha[0].is_zero());

    stha.clear();

    typename PNT_1_traits::Sturm_habicht_sequence()
        (f,std::back_inserter(stha));

    assert(stha.size()==6);

    assert(stha[5]==f);
    assert(stha[4]==fx);
    assert(stha[3]==PNT_1(NT(1024),NT(5568),NT(9984),NT(5888)));
    assert(stha[2]==PNT_1(NT(55296),NT(221184),NT(221184)));
    assert(stha[1].is_zero());
    assert(stha[0].is_zero());

  }

  { // test for cofactors
    // f:=(2*x+1)^3*(x+1)*(x-1)
    PNT_1 f(NT(-1),NT(-6),NT(-11),NT(-2),NT(12),NT(8));
    PNT_1 fx = CGAL::differentiate(f);

    std::vector<PNT_1> stha_check,stha,co_f,co_fx;

    CGAL::sturm_habicht_sequence
        (f,std::back_inserter(stha_check));
    CGAL::sturm_habicht_sequence_with_cofactors
        (f,
         std::back_inserter(stha),
         std::back_inserter(co_f),
         std::back_inserter(co_fx)
        );
    
    int n = static_cast<int>(stha.size());
    assert(n == static_cast<int>(stha_check.size()));
    for(int i = 0; i < n; i++) {
      assert(stha[i]==stha_check[i]);
      assert(stha[i]==co_f[i]*f + co_fx[i]*fx);
    }

    stha.clear();
    co_f.clear();
    co_fx.clear();
    typename PNT_1_traits::Sturm_habicht_sequence_with_cofactors()
        (f,
         std::back_inserter(stha),
         std::back_inserter(co_f),
         std::back_inserter(co_fx)
        );
    
    n = static_cast<int>(stha.size());
    assert(n == static_cast<int>(stha_check.size()));
    for(int i = 0; i < n; i++) {
      assert(stha[i]==stha_check[i]);
      assert(stha[i]==co_f[i]*f + co_fx[i]*fx);
    }

  }

  { // Cofactors, defective case
    // f:=x^8+x^2
    PNT_1 f(NT(0),NT(0),NT(1),NT(0),NT(0),NT(0),NT(0),NT(0),NT(1));
    PNT_1 fx = CGAL::differentiate(f);

    std::vector<PNT_1> stha_check,stha,co_f,co_fx;


    CGAL::sturm_habicht_sequence
        (f,std::back_inserter(stha_check));
    CGAL::sturm_habicht_sequence_with_cofactors
        (f,
         std::back_inserter(stha),
         std::back_inserter(co_f),
         std::back_inserter(co_fx)
        );
    
    int n = static_cast<int>(stha.size());
    assert(n == static_cast<int>(stha_check.size()));
    for(int i = 0; i < n; i++) {
      assert(stha[i]==stha_check[i]);
      assert(stha[i]==co_f[i]*f + co_fx[i]*fx);
    }


  }

  { // Test for bivariate polynomials
    // f:=y^4 + 6*y^3 + 9*y^2 + (9*x^3 + 3*x^2 + 2*x + 5)*y + (2*x^2 + 5*x + 7)
    std::stringstream ss("P[4(0,P[2(0,7)(1,5)(2,2)])(1,P[3(0,5)(1,2)(2,3)(3,9)])(2,P[0(0,9)])(3,P[0(0,6)])(4,P[0(0,1)])]");
    PNT_2 f;
    ss >> f;

    {
      std::vector<PNT_1> stha,costha;
      CGAL::internal::first_two_sturm_habicht_coefficients<PNT_2_traits>
          (f,std::back_inserter(stha),std::back_inserter(costha));
      
      assert(stha.size()==5);
      assert(costha.size()==4);
      
      NT c1[13] = { NT(-81383),NT(-113448),NT(-123840),NT(-223936),NT(-272142),
                    NT(-14700),NT(419954),NT(156816),NT(-247131),NT(-498636),
                    NT(-275562),NT(-236196),NT(-177147) };
      PNT_1 p1(c1,c1+13);
      assert(stha[0] == p1);
      assert(stha[1] == PNT_1(NT(-828),NT(-1008),NT(-864),NT(-1728),NT(-1620),
                                NT(-1944),NT(-2916)) );
      assert(stha[2] == PNT_1(NT(36)) );
      assert(stha[3] == PNT_1(NT(4)) );
      assert(stha[4] == PNT_1(NT(1)) );
      
      
      
      assert(costha[0] == PNT_1(NT(-2742),NT(-2592),NT(-1788),NT(-2100),
                                  NT(-1638),NT(108),NT(1458)) );
      assert(costha[1] == PNT_1(NT(48),NT(-24),NT(-36),NT(-108)) );
      assert(costha[2] == PNT_1(NT(18)) );
      assert(costha[3] == PNT_1(NT(6)) );
    }
    {
      std::vector<PNT_1> stha,costha;
      CGAL::internal::prs_first_two_sturm_habicht_coefficients<PNT_2_traits>
          (f,std::back_inserter(stha),std::back_inserter(costha));
      
      assert(stha.size()==5);
      assert(costha.size()==4);
      
      NT c1[13] = { NT(-81383),NT(-113448),NT(-123840),NT(-223936),NT(-272142),
                    NT(-14700),NT(419954),NT(156816),NT(-247131),NT(-498636),
                    NT(-275562),NT(-236196),NT(-177147) };
      PNT_1 p1(c1,c1+13);
      assert(stha[0] == p1);
      assert(stha[1] == PNT_1(NT(-828),NT(-1008),NT(-864),NT(-1728),NT(-1620),
                                NT(-1944),NT(-2916)) );
      assert(stha[2] == PNT_1(NT(36)) );
      assert(stha[3] == PNT_1(NT(4)) );
      assert(stha[4] == PNT_1(NT(1)) );
      
      
      
      assert(costha[0] == PNT_1(NT(-2742),NT(-2592),NT(-1788),NT(-2100),
                                  NT(-1638),NT(108),NT(1458)) );
      assert(costha[1] == PNT_1(NT(48),NT(-24),NT(-36),NT(-108)) );
      assert(costha[2] == PNT_1(NT(18)) );
      assert(costha[3] == PNT_1(NT(6)) );
    }
    {
      std::vector<PNT_1> stha,costha;
      CGAL::internal::bezout_first_two_sturm_habicht_coefficients<PNT_2_traits>
        (f,std::back_inserter(stha),
         std::back_inserter(costha));
      
      assert(stha.size()==5);
      assert(costha.size()==4);
      
      NT c1[13] = { NT(-81383),NT(-113448),NT(-123840),NT(-223936),NT(-272142),
                    NT(-14700),NT(419954),NT(156816),NT(-247131),NT(-498636),
                    NT(-275562),NT(-236196),NT(-177147) };
      PNT_1 p1(c1,c1+13);
      assert(stha[0] == p1);
      assert(stha[1] == PNT_1(NT(-828),NT(-1008),NT(-864),NT(-1728),NT(-1620),
                                NT(-1944),NT(-2916)) );
      assert(stha[2] == PNT_1(NT(36)) );
      assert(stha[3] == PNT_1(NT(4)) );
      assert(stha[4] == PNT_1(NT(1)) );
      
      
      
      assert(costha[0] == PNT_1(NT(-2742),NT(-2592),NT(-1788),NT(-2100),
                                  NT(-1638),NT(108),NT(1458)) );
      assert(costha[1] == PNT_1(NT(48),NT(-24),NT(-36),NT(-108)) );
      assert(costha[2] == PNT_1(NT(18)) );
      assert(costha[3] == PNT_1(NT(6)) );
    }
    
    
  }
  


}

int main() {

#ifdef CGAL_USE_LEDA
    test_routine<CGAL::LEDA_arithmetic_kernel>();
#else
    std::cout << "LEDA tests skipped" << std::endl;
#endif

#ifdef CGAL_USE_CORE
    test_routine<CGAL::CORE_arithmetic_kernel>();
#else
    std::cout << "CORE tests skipped" << std::endl;
#endif
    return 0;
}
