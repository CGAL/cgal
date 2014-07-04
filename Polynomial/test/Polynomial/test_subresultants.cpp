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
#include <CGAL/Polynomial/subresultants.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>

template<typename Poly_> Poly_ from_string(const char* s) {
  std::stringstream ss(s);
  Poly_ f;
  ss >> f;
  return f;
}

template<typename ArithmeticKernel>
void test_routine() {
  typedef ArithmeticKernel Arithmetic_kernel;

  typedef typename Arithmetic_kernel::Rational Rational;
  typedef typename Arithmetic_kernel::Integer Integer;
    

  typedef CGAL::Polynomial<Integer> Poly_int1;
  typedef CGAL::Polynomial_traits_d<Poly_int1> Poly_int1_traits;

  typedef CGAL::Polynomial<Poly_int1> Poly_int2;
  typedef CGAL::Polynomial_traits_d<Poly_int2> Poly_int2_traits;

  typedef CGAL::Polynomial<Rational> Poly_rat1;
  typedef CGAL::Polynomial_traits_d<Poly_rat1> Poly_rat1_traits;

  {
    //Example for the regular case:
    Poly_int1 f(-5,-2,3,-6,-7,3,-2,4);
    Poly_int1 g(0,-1,5,-7,5);
    std::vector<Poly_int1> sres;
    CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
        (f,g,std::back_inserter(sres));
    assert(sres.size()==5);
    assert(sres[4]==Integer(25)*g);
    assert(sres[3]==Poly_int1(-3125,-1768,4970,-9451));
    assert(sres[2]==Poly_int1(206535,-262340,252423));
    assert(sres[1]==Poly_int1(Integer(602925),Integer(657683)));
    assert(sres[0]==Poly_int1(Integer(4474810)));
    sres.clear();
    CGAL::internal::bezout_polynomial_subresultants<Poly_int1_traits>
        (f,g,std::back_inserter(sres));
    assert(sres.size()==5);
    assert(sres[4]==Integer(25)*g);
    assert(sres[3]==Poly_int1(-3125,-1768,4970,-9451));
    assert(sres[2]==Poly_int1(206535,-262340,252423));
    assert(sres[1]==Poly_int1(Integer(602925),Integer(657683)));
    assert(sres[0]==Poly_int1(Integer(4474810)));
    sres.clear();
    CGAL::polynomial_subresultants(f,g,std::back_inserter(sres));
    assert(sres.size()==5);
    assert(sres[4]==Integer(25)*g);
    assert(sres[3]==Poly_int1(-3125,-1768,4970,-9451));
    assert(sres[2]==Poly_int1(206535,-262340,252423));
    assert(sres[1]==Poly_int1(Integer(602925),Integer(657683)));
    assert(sres[0]==Poly_int1(Integer(4474810)));
    sres.clear();
    typename Poly_int1_traits::Polynomial_subresultants()
        (f,g,std::back_inserter(sres));
    assert(sres.size()==5);
    assert(sres[4]==Integer(25)*g);
    assert(sres[3]==Poly_int1(-3125,-1768,4970,-9451));
    assert(sres[2]==Poly_int1(206535,-262340,252423));
    assert(sres[1]==Poly_int1(Integer(602925),Integer(657683)));
    assert(sres[0]==Poly_int1(Integer(4474810)));
    std::vector<Integer> psres;
    CGAL::internal::prs_principal_subresultants<Poly_int1_traits>
        (f,g,std::back_inserter(psres));
    assert(psres.size()==5);
    assert(psres[4]==125);
    assert(psres[3]==-9451);
    assert(psres[2]==252423);
    assert(psres[1]==657683);
    assert(psres[0]==4474810);
    psres.clear();
    CGAL::internal::bezout_principal_subresultants<Poly_int1_traits>
        (f,g,std::back_inserter(psres));
    assert(psres.size()==5);
    assert(psres[4]==125);
    assert(psres[3]==-9451);
    assert(psres[2]==252423);
    assert(psres[1]==657683);
    assert(psres[0]==4474810);
    psres.clear();
    CGAL::principal_subresultants(f,g,std::back_inserter(psres));
    assert(psres.size()==5);
    assert(psres[4]==125);
    assert(psres[3]==-9451);
    assert(psres[2]==252423);
    assert(psres[1]==657683);
    assert(psres[0]==4474810);
    psres.clear();
    typename Poly_int1_traits::Principal_subresultants()
        (f,g,std::back_inserter(psres));
    assert(psres.size()==5);
    assert(psres[4]==125);
    assert(psres[3]==-9451);
    assert(psres[2]==252423);
    assert(psres[1]==657683);
    assert(psres[0]==4474810);
  }

  {    
    // Defective exampke
    Poly_int1 f(1,-2,0,0,0,0,1);
    std::vector<Poly_int1> sres;
    CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
      (f,CGAL::differentiate(f),std::back_inserter(sres));
    assert(sres.size()==6);
    assert(sres[5]==CGAL::differentiate(f));
    assert(sres[4]==Poly_int1(Integer(36),Integer(-60)));
    assert(sres[3]==Poly_int1(0));
    assert(sres[2]==Poly_int1(0));
    assert(sres[1]==Poly_int1(Integer(-36000),Integer(60000)));
    assert(sres[0]==Poly_int1(Integer(-153344)));
    sres.clear();
    CGAL::internal::bezout_polynomial_subresultants<Poly_int1_traits>
      (f,CGAL::differentiate(f),std::back_inserter(sres));
    assert(sres.size()==6);
    assert(sres[5]==CGAL::differentiate(f));
    assert(sres[4]==Poly_int1(Integer(36),Integer(-60)));
    assert(sres[3]==Poly_int1(0));
    assert(sres[2]==Poly_int1(0));
    assert(sres[1]==Poly_int1(Integer(-36000),Integer(60000)));
    assert(sres[0]==Poly_int1(Integer(-153344)));
    std::vector<Integer> psres;
    CGAL::internal::prs_principal_subresultants<Poly_int1_traits>
        (f,CGAL::differentiate(f),std::back_inserter(psres));
    assert(psres.size()==6);
    assert(psres[5]==6);
    assert(psres[4]==0);
    assert(psres[3]==0);
    assert(psres[2]==0);
    assert(psres[1]==60000);
    assert(psres[0]==-153344);
    psres.clear();
    CGAL::internal::bezout_principal_subresultants<Poly_int1_traits>
        (f,CGAL::differentiate(f),std::back_inserter(psres));
    assert(psres.size()==6);
    assert(psres[5]==6);
    assert(psres[4]==0);
    assert(psres[3]==0);
    assert(psres[2]==0);
    assert(psres[1]==60000);
    assert(psres[0]==-153344);
    

  }
    {
    //Example for rational values:
    Poly_rat1 f(-5,-2,3,-6,-7,3,-2,4);
    Poly_rat1 g(0,-1,5,-7,5);
    std::vector<Poly_rat1> sres;
    CGAL::internal::prs_polynomial_subresultants<Poly_rat1_traits>
        (f,g,std::back_inserter(sres));
    assert(sres.size()==5);
    assert(sres[4]==Rational(25)*g);
    assert(sres[3]==Poly_rat1(-3125,-1768,4970,-9451));
    assert(sres[2]==Poly_rat1(206535,-262340,252423));
    assert(sres[1]==Poly_rat1(Rational(602925),Rational(657683)));
    assert(sres[0]==Poly_rat1(Rational(4474810)));
    sres.clear();
    CGAL::internal::bezout_polynomial_subresultants<Poly_rat1_traits>
        (f,g,std::back_inserter(sres));
    assert(sres.size()==5);
    assert(sres[4]==Rational(25)*g);
    assert(sres[3]==Poly_rat1(-3125,-1768,4970,-9451));
    assert(sres[2]==Poly_rat1(206535,-262340,252423));
    assert(sres[1]==Poly_rat1(Rational(602925),Rational(657683)));
    assert(sres[0]==Poly_rat1(Rational(4474810)));
    sres.clear();
    CGAL::polynomial_subresultants(f,g,std::back_inserter(sres));
    assert(sres.size()==5);
    assert(sres[4]==Rational(25)*g);
    assert(sres[3]==Poly_rat1(-3125,-1768,4970,-9451));
    assert(sres[2]==Poly_rat1(206535,-262340,252423));
    assert(sres[1]==Poly_rat1(Rational(602925),Rational(657683)));
    assert(sres[0]==Poly_rat1(Rational(4474810)));
    std::vector<Rational> psres;
    CGAL::internal::prs_principal_subresultants<Poly_rat1_traits>
        (f,g,std::back_inserter(psres));
    assert(psres.size()==5);
    assert(psres[4]==125);
    assert(psres[3]==-9451);
    assert(psres[2]==252423);
    assert(psres[1]==657683);
    assert(psres[0]==4474810);
    psres.clear();
    CGAL::internal::bezout_principal_subresultants<Poly_rat1_traits>
        (f,g,std::back_inserter(psres));
    assert(psres.size()==5);
    assert(psres[4]==125);
    assert(psres[3]==-9451);
    assert(psres[2]==252423);
    assert(psres[1]==657683);
    assert(psres[0]==4474810);
    psres.clear();
    CGAL::principal_subresultants(f,g,std::back_inserter(psres));
    assert(psres.size()==5);
    assert(psres[4]==125);
    assert(psres[3]==-9451);
    assert(psres[2]==252423);
    assert(psres[1]==657683);
    assert(psres[0]==4474810);
    psres.clear();
  }

  {
    // Another test included on 10 Jan 2007 - situation comes from arrangement
    // computation of curves:
    
    // f:=y^4 + (x)*y^2 + (2*x^4 + (-1)*x^3);
    Poly_int2 f=from_string<Poly_int2>("P[4(0,P[4(3,-1)(4,2)])(2,P[1(1,1)])(4,P[0(0,1)])]");
    // g:=g:=eval(10000*f,[y=2*y,x=x-1/10]);
    Poly_int2 g=from_string<Poly_int2>("P[4(0,P[4(0,12)(1,-380)(2,4200)(3,-18000)(4,20000)])(2,P[1(0,-4000)(1,40000)])(4,P[0(0,160000)])]");
    std::vector<Poly_int2> sres;
    CGAL::internal::prs_polynomial_subresultants<Poly_int2_traits>
        (f,g,std::back_inserter(sres));
    assert(sres.size()==5);
    // Computed with MAPLE
    assert(sres[4]==g);
    assert(sres[3]==from_string<Poly_int2>("P[2(0,P[4(0,12)(1,-380)(2,4200)(3,142000)(4,-300000)])(2,P[1(0,-4000)(1,-120000)])]"));
    assert(sres[2]==from_string<Poly_int2>("P[2(0,P[5(0,-48000)(1,80000)(2,28800000)(3,-1072000000)(4,-15840000000)(5,36000000000)])(2,P[2(0,16000000)(1,960000000)(2,14400000000)])]"));
    assert(sres[1]==from_string<Poly_int2>("P[0(0,P[9(0,576000)(1,172800000)(2,5326400000)(3,-158512000000)(4,-5164000000000)(5,24705600000000)(6,615472000000000)(7,912480000000000)(8,-9864000000000000)(9,10800000000000000)])]"));
    assert(sres[0]==from_string<Poly_int2>("P[0(0,P[16(0,20736)(1,11197440)(2,1559232000)(3,5760000)(4,-3426163040000)(5,-9736288000000)(6,2377866144000000)(7,-1780931200000000)(8,-427278798400000000)(9,-507616640000000000)(10,31454608000000000000)(11,83909222400000000000)(12,-697197584000000000000)(13,-919113600000000000000)(14,9138960000000000000000)(15,-15336000000000000000000)(16,8100000000000000000000)])]"));
    sres.clear();
    CGAL::internal::bezout_polynomial_subresultants<Poly_int2_traits>
        (f,g,std::back_inserter(sres));
    assert(sres.size()==5);
    // Computed with MAPLE
    assert(sres[4]==g);
    assert(sres[3]==from_string<Poly_int2>("P[2(0,P[4(0,12)(1,-380)(2,4200)(3,142000)(4,-300000)])(2,P[1(0,-4000)(1,-120000)])]"));
    assert(sres[2]==from_string<Poly_int2>("P[2(0,P[5(0,-48000)(1,80000)(2,28800000)(3,-1072000000)(4,-15840000000)(5,36000000000)])(2,P[2(0,16000000)(1,960000000)(2,14400000000)])]"));
    assert(sres[1]==from_string<Poly_int2>("P[0(0,P[9(0,576000)(1,172800000)(2,5326400000)(3,-158512000000)(4,-5164000000000)(5,24705600000000)(6,615472000000000)(7,912480000000000)(8,-9864000000000000)(9,10800000000000000)])]"));
    assert(sres[0]==from_string<Poly_int2>("P[0(0,P[16(0,20736)(1,11197440)(2,1559232000)(3,5760000)(4,-3426163040000)(5,-9736288000000)(6,2377866144000000)(7,-1780931200000000)(8,-427278798400000000)(9,-507616640000000000)(10,31454608000000000000)(11,83909222400000000000)(12,-697197584000000000000)(13,-919113600000000000000)(14,9138960000000000000000)(15,-15336000000000000000000)(16,8100000000000000000000)])]"));
    sres.clear();
    CGAL::polynomial_subresultants(f,g,std::back_inserter(sres));
    assert(sres.size()==5);
    // Computed with MAPLE
    assert(sres[4]==g);
    assert(sres[3]==from_string<Poly_int2>("P[2(0,P[4(0,12)(1,-380)(2,4200)(3,142000)(4,-300000)])(2,P[1(0,-4000)(1,-120000)])]"));
    assert(sres[2]==from_string<Poly_int2>("P[2(0,P[5(0,-48000)(1,80000)(2,28800000)(3,-1072000000)(4,-15840000000)(5,36000000000)])(2,P[2(0,16000000)(1,960000000)(2,14400000000)])]"));
    assert(sres[1]==from_string<Poly_int2>("P[0(0,P[9(0,576000)(1,172800000)(2,5326400000)(3,-158512000000)(4,-5164000000000)(5,24705600000000)(6,615472000000000)(7,912480000000000)(8,-9864000000000000)(9,10800000000000000)])]"));
    assert(sres[0]==from_string<Poly_int2>("P[0(0,P[16(0,20736)(1,11197440)(2,1559232000)(3,5760000)(4,-3426163040000)(5,-9736288000000)(6,2377866144000000)(7,-1780931200000000)(8,-427278798400000000)(9,-507616640000000000)(10,31454608000000000000)(11,83909222400000000000)(12,-697197584000000000000)(13,-919113600000000000000)(14,9138960000000000000000)(15,-15336000000000000000000)(16,8100000000000000000000)])]"));
  }
  
  { // Test for f.degree() < g.degree()
    Poly_int1 f = from_string<Poly_int1>("P[10(2,-4)(3,-73)(6,97)(7,-62)(10,-56)]");
    Poly_int1 g = from_string<Poly_int1>("P[15(5,-75)(6,-17)(8,71)(11,-44)(14,80)(15,-82)]");
    
    std::vector<Poly_int1> sres;
    CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
      (f,g,std::back_inserter(sres));
    assert(sres.size()==11);
    assert(sres[10]==from_string<Poly_int1>("P[10(2,-39337984)(3,-717918208)(6,953946112)(7,-609738752)(10,-550731776)]"));
    assert(sres[9]==from_string<Poly_int1>("P[9(2,-305262755840)(3,-4966105776128)(4,10840151891968)(5,-5962969628672)(6,6702091010048)(7,-22436989575168)(8,19712814514176)(9,-3099911815168)]"));
    assert(sres[8]==from_string<Poly_int1>("P[8(2,9680208765108224)(3,156728928112275456)(4,-360058057701515264)(5,152421565222248448)(6,-176106510157656064)(7,746064242175737856)(8,-579306262853885952)]"));
    assert(sres[7]==from_string<Poly_int1>("P[7(2,-1486764415342623952896)(3,-26708453688140199103488)(4,8054371431842323176960)(5,3491478546579489027072)(6,38677404054075120663552)(7,-43611415672712033527296)]"));
    assert(sres[6]==from_string<Poly_int1>("P[6(2,9979660953345863413381568)(3,193906324691706110661441648)(4,213219021894852075965695648)(5,362883385829743140658019008)(6,-93331232763027671680558576)]"));
    assert(sres[5]==from_string<Poly_int1>("P[5(2,-70907443741718133356344565824)(3,-1389117178389007568663227054224)(4,-1747565962047585993090805107776)(5,-2771071359404167222865326731904)]"));
    assert(sres[4]==from_string<Poly_int1>("P[4(2,-715876861001937085090535662991104)(3,-13330409138559036348802401660667328)(4,-5258395198219675575175450753910144)]"));
    assert(sres[3]==from_string<Poly_int1>("P[3(2,-2842398393513781389004679512270277696)(3,-51817649494648561412362002061296745104)]"));
    assert(sres[2]==from_string<Poly_int1>("P[2(2,-46088453959546967380998470482217286656)]"));
    assert(sres[1]==Poly_int1(Integer(0)));
    assert(sres[0]==Poly_int1(Integer(0)));
    sres.clear();

    CGAL::internal::bezout_polynomial_subresultants<Poly_int1_traits>
      (f,g,std::back_inserter(sres));
    assert(sres.size()==11);
    assert(sres[10]==from_string<Poly_int1>("P[10(2,-39337984)(3,-717918208)(6,953946112)(7,-609738752)(10,-550731776)]"));
    assert(sres[9]==from_string<Poly_int1>("P[9(2,-305262755840)(3,-4966105776128)(4,10840151891968)(5,-5962969628672)(6,6702091010048)(7,-22436989575168)(8,19712814514176)(9,-3099911815168)]"));
    assert(sres[8]==from_string<Poly_int1>("P[8(2,9680208765108224)(3,156728928112275456)(4,-360058057701515264)(5,152421565222248448)(6,-176106510157656064)(7,746064242175737856)(8,-579306262853885952)]"));
    assert(sres[7]==from_string<Poly_int1>("P[7(2,-1486764415342623952896)(3,-26708453688140199103488)(4,8054371431842323176960)(5,3491478546579489027072)(6,38677404054075120663552)(7,-43611415672712033527296)]"));
    assert(sres[6]==from_string<Poly_int1>("P[6(2,9979660953345863413381568)(3,193906324691706110661441648)(4,213219021894852075965695648)(5,362883385829743140658019008)(6,-93331232763027671680558576)]"));
    assert(sres[5]==from_string<Poly_int1>("P[5(2,-70907443741718133356344565824)(3,-1389117178389007568663227054224)(4,-1747565962047585993090805107776)(5,-2771071359404167222865326731904)]"));
    assert(sres[4]==from_string<Poly_int1>("P[4(2,-715876861001937085090535662991104)(3,-13330409138559036348802401660667328)(4,-5258395198219675575175450753910144)]"));
    assert(sres[3]==from_string<Poly_int1>("P[3(2,-2842398393513781389004679512270277696)(3,-51817649494648561412362002061296745104)]"));
    assert(sres[2]==from_string<Poly_int1>("P[2(2,-46088453959546967380998470482217286656)]"));
    assert(sres[1]==Poly_int1(Integer(0)));
    assert(sres[0]==Poly_int1(Integer(0)));
    sres.clear();
  }
  { // Another example where sign switches are necessary
    Poly_int1 f = from_string<Poly_int1>("P[3(0,42)(1,-40)(2,-7)(3,-10)]");
    Poly_int1 g = from_string<Poly_int1>("P[5(0,74)(1,6)(2,-92)(3,75)(4,23)(5,-50)]");
    std::vector<Poly_int1> sres;
    CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
      (f,g,std::back_inserter(sres));
    assert(sres.size()==4);
    assert(sres[3]==from_string<Poly_int1>("P[3(0,-420)(1,400)(2,70)(3,100)]"));
    assert(sres[2]==from_string<Poly_int1>("P[2(0,-1058480)(1,688000)(2,698080)]"));
    assert(sres[1]==from_string<Poly_int1>("P[1(0,-22577275200)(1,28253151360)]"));
    assert(sres[0]==from_string<Poly_int1>("P[0(0,-103066942158720)]"));
    sres.clear();
    CGAL::internal::bezout_polynomial_subresultants<Poly_int1_traits>
      (f,g,std::back_inserter(sres));
    assert(sres.size()==4);
    assert(sres[3]==from_string<Poly_int1>("P[3(0,-420)(1,400)(2,70)(3,100)]"));
    assert(sres[2]==from_string<Poly_int1>("P[2(0,-1058480)(1,688000)(2,698080)]"));
    assert(sres[1]==from_string<Poly_int1>("P[1(0,-22577275200)(1,28253151360)]"));
    assert(sres[0]==from_string<Poly_int1>("P[0(0,-103066942158720)]"));
  }  
  { // Test cases with constant polynomials
    {
      Poly_int1 f(Integer(8));
      Poly_int1 g=from_string<Poly_int1>("P[3(0,43)(1,40)(2,-19)(3,2)]");
      std::vector<Poly_int1> sres;

      CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
        (f,g,std::back_inserter(sres));
      assert(sres.size()==1);
      //assert(sres[0].degree()==0);
      assert(sres[0]==Integer(8*8*8));
      sres.clear();

      CGAL::internal::bezout_polynomial_subresultants<Poly_int1_traits>
        (f,g,std::back_inserter(sres));
      assert(sres.size()==1);
      //assert(sres[0].degree()==0);
      assert(sres[0]==Integer(8*8*8));
      std::vector<Integer> psres;

      CGAL::internal::prs_principal_subresultants<Poly_int1_traits>
        (g,f,std::back_inserter(psres));
      assert(psres.size()==1);
      assert(psres[0]==Integer(8*8*8));
      psres.clear();

      CGAL::internal::bezout_principal_subresultants<Poly_int1_traits>
        (g,f,std::back_inserter(psres));
      assert(psres.size()==1);
      assert(psres[0]==Integer(8*8*8));
    }
    {
       Poly_int1 f(Integer(0));
       Poly_int1 g=from_string<Poly_int1>("P[1(0,1)(1,1)]");
       std::vector<Poly_int1> sres;
       CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
        (f,g,std::back_inserter(sres));
       assert(sres.size()==1);
       assert(CGAL::is_zero(sres[0]));
       sres.clear();

      CGAL::internal::bezout_polynomial_subresultants<Poly_int1_traits>
        (f,g,std::back_inserter(sres));
      assert(sres.size()==1);
      assert(CGAL::is_zero(sres[0]));
    }
    {
      Poly_int1 f(Integer(7));
      Poly_int1 g(Integer(-12));
      
      std::vector<Integer> sres;
      CGAL::internal::prs_principal_subresultants<Poly_int1_traits>
        (f,g,std::back_inserter(sres));
      assert(sres.size()==1);
      assert(sres[0]==Integer(1));
      sres.clear();
      
      CGAL::internal::bezout_principal_subresultants<Poly_int1_traits>
        (f,g,std::back_inserter(sres));
      assert(sres.size()==1);
      assert(sres[0]==Integer(1));
    }
    {
      Poly_int2 f=from_string<Poly_int2>("P[0(0,P[0(0,122)])]");
      Poly_int2 g=from_string<Poly_int2>("P[0(0,P[0(0,0)])]");
      std::vector<Poly_int2> sres;
      CGAL::internal::prs_polynomial_subresultants<Poly_int2_traits>
        (f,g,std::back_inserter(sres));
      assert(sres.size()==1);
      assert(sres[0].is_zero());
      sres.clear();
      
      CGAL::internal::bezout_polynomial_subresultants<Poly_int2_traits>
        (f,g,std::back_inserter(sres));
      assert(sres.size()==1);
      assert(sres[0].is_zero());
    }
    
  }
  { // Test for two linear polynomials:
    Poly_int1 f = from_string<Poly_int1>("P[1(0,-5)(1,6)]");
    Poly_int1 g = from_string<Poly_int1>("P[1(0,-3)(1,4)]");
    std::vector<Poly_int1> sres;
    CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
      (f,g,std::back_inserter(sres));
    assert(sres.size()==2);
    assert(sres[1]==g);
    assert(sres[0]==from_string<Poly_int1>("P[0(0,2)]"));
    sres.clear();
    CGAL::internal::bezout_polynomial_subresultants<Poly_int1_traits>
      (f,g,std::back_inserter(sres));
    assert(sres.size()==2);
    assert(sres[1]==g);
    assert(sres[0]==from_string<Poly_int1>("P[0(0,2)]"));
  }
  { // Test for cofactors, univariate
    Poly_int1 f = from_string<Poly_int1>("P[6(0,246)(1,100)(2,197)(3,14)(4,136)(5,191)(6,207)]");
    Poly_int1 g = from_string<Poly_int1>("P[5(0,100)(1,394)(2,42)(3,544)(4,955)(5,1242)]");
    std::vector<Poly_int1> sres_check,sres,coP,coQ;
    CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
      (f,g,std::back_inserter(sres_check));
    CGAL::polynomial_subresultants_with_cofactors
      (f,g,
       std::back_inserter(sres), 
       std::back_inserter(coP),
       std::back_inserter(coQ));
    assert(sres.size()==sres_check.size());
    assert(sres.size()==coP.size());
    assert(sres.size()==coQ.size());
    for(int i=0;i < static_cast<int>(sres.size()); i++) {
      assert(sres[i]==sres_check[i]);
      assert(sres[i]==coP[i]*f + coQ[i]*g);
    }
  }

    { // Test for cofactors, bivariate
    Poly_int2 f = from_string<Poly_int2>("P[7(0,P[7(0,72)(1,113)(2,238)(3,75)(4,77)(5,149)(6,34)(7,75)])(1,P[6(0,113)(1,69)(2,94)(3,171)(4,148)(5,103)(6,42)])(2,P[5(0,103)(1,233)(2,16)(3,131)(4,2)(5,156)])(3,P[4(0,249)(1,194)(2,156)(3,81)(4,176)])(4,P[3(0,39)(1,140)(2,134)(3,190)])(5,P[2(0,117)(1,249)(2,158)])(6,P[1(0,109)(1,46)])(7,P[0(0,236)])]");
    Poly_int2 g = from_string<Poly_int2>("P[6(0,P[6(0,113)(1,69)(2,94)(3,171)(4,148)(5,103)(6,42)])(1,P[5(0,206)(1,466)(2,32)(3,262)(4,4)(5,312)])(2,P[4(0,747)(1,582)(2,468)(3,243)(4,528)])(3,P[3(0,156)(1,560)(2,536)(3,760)])(4,P[2(0,585)(1,1245)(2,790)])(5,P[1(0,654)(1,276)])(6,P[0(0,1652)])]");
    std::vector<Poly_int2> sres_check,sres,coP,coQ;
    CGAL::internal::prs_polynomial_subresultants<Poly_int2_traits>
      (f,g,std::back_inserter(sres_check));
    CGAL::polynomial_subresultants_with_cofactors
      (f,g,
       std::back_inserter(sres), 
       std::back_inserter(coP),
       std::back_inserter(coQ));
    assert(sres.size()==sres_check.size());
    assert(sres.size()==coP.size());
    assert(sres.size()==coQ.size());
    for(int i=0;i < static_cast<int>(sres.size()); i++) {
      assert(sres[i]==sres_check[i]);
      assert(sres[i]==coP[i]*f + coQ[i]*g);
    }
    }

  { // Test for trivariate
      typedef CGAL::Polynomial<Poly_int2> Poly_int3;
      typedef CGAL::Polynomial_traits_d<Poly_int3> Poly_int3_traits;
      Poly_int3 f = from_string<Poly_int3>("P[6(0,P[6(0,P[6(2,3)(4,-3)(6,1)])(1,P[1(1,-2)])(2,P[4(0,3)(2,-5)(4,3)])(4,P[2(0,-3)(2,3)])(6,P[0(0,1)])])(3,P[1(0,P[0(0,2)])(1,P[1(1,-2)])])(6,P[0(0,P[0(0,1)])])]");
      Poly_int3 g = CGAL::differentiate(f);
      std::vector<Poly_int3> sres_check,sres,coP,coQ;
      CGAL::internal::prs_polynomial_subresultants<Poly_int3_traits>
          (f,g,std::back_inserter(sres_check));
      CGAL::polynomial_subresultants_with_cofactors
          (f,g,
           std::back_inserter(sres), 
           std::back_inserter(coP),
           std::back_inserter(coQ));
      assert(sres.size()==sres_check.size());
      assert(sres.size()==coP.size());
      assert(sres.size()==coQ.size());
      for(int i=static_cast<int>(sres.size())-1;i>=0 ; i--) {
	assert(sres[i]==sres_check[i]);
	assert(sres[i]==coP[i]*f + coQ[i]*g);
      }
  }
    
    Poly_int2 x=from_string<Poly_int2>("P[0(0,P[1(1,1)])]");
    Poly_int2 y=from_string<Poly_int2>("P[1(1,P[0(0,1)])]");
    { // Bug reported by Eric Berberich, 15.03.2010 -> (half) fixed
      Poly_int2 f = y*y + x*x-2;
      Poly_int2 g = -y*y + x*x;
      std::vector<Poly_int2> sres,coP,coQ,sres_check;
      CGAL::internal::prs_polynomial_subresultants<Poly_int2_traits>
	(f,g,std::back_inserter(sres_check));      
      CGAL::polynomial_subresultants_with_cofactors
	(f,g,
	 std::back_inserter(sres),
	 std::back_inserter(coP),
	 std::back_inserter(coQ));
      assert(sres.size()==3);
      assert(sres[2]==coP[2]*f + coQ[2]*g);
      assert(sres[1]==coP[1]*f + coQ[1]*g);
      assert(sres[0]==coP[0]*f + coQ[0]*g);
      CGAL::set_pretty_mode(std::cout);
      assert(sres[2]==sres_check[2]);
      assert(sres[1]==sres_check[1]);
      assert(sres[0]==sres_check[0]);
    }
    Poly_int1 t=from_string<Poly_int1>("P[1(1,1)]");   
    {
      Poly_int1 f = 3*CGAL::ipower(t,10) + 7*CGAL::ipower(t,4)-13;
      Poly_int1 g = (-5)*CGAL::ipower(t,10) + 11*CGAL::ipower(t,6)-17*t;
      std::vector<Poly_int1> sres,coP,coQ,sres_check;
      CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
	(f,g,std::back_inserter(sres_check));      
      CGAL::polynomial_subresultants_with_cofactors
	(f,g,
	 std::back_inserter(sres),
	 std::back_inserter(coP),
	 std::back_inserter(coQ));
      CGAL::set_pretty_mode(std::cout);
      for(int i=static_cast<int>(sres.size()-1);i>=0;i--) {
	assert(sres[i]==coP[i]*f + coQ[i]*g);
	assert(sres_check[i]==sres[i]);
      }
    }

    {
      Poly_int1 f = 11*t;
      for(int i =0; i < 25; i++) {
	f=t*f+((613*i+1225)%47);
      }
      Poly_int1 g = -23*t;
      for(int i =0; i < 25; i++) {
	g=t*g+((397*i+2423)%59);
      }
      std::vector<Poly_int1> sres,coP,coQ,sres_check;
      CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
	(f,g,std::back_inserter(sres_check));      
      CGAL::polynomial_subresultants_with_cofactors
	(f,g,
	 std::back_inserter(sres),
	 std::back_inserter(coP),
	 std::back_inserter(coQ));
      CGAL::set_pretty_mode(std::cout);
      for(int i=static_cast<int>(sres.size()-1);i>=0;i--) {
	assert(sres[i]==coP[i]*f + coQ[i]*g);
	assert(sres_check[i]==sres[i]);
      }
    }

    {
      Poly_int1 f = 771*CGAL::ipower(t,7);
      for(int i =0; i < 31; i+=5) {
	f=t*f+((881*i-14154)%31);
      }
      Poly_int1 g = 919*CGAL::ipower(t,7);
      for(int i =0; i < 31; i+=4) {
	g=t*g+((129*i+1151)%59);
      }
      std::vector<Poly_int1> sres,coP,coQ,sres_check;
      CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
	(f,g,std::back_inserter(sres_check));      
      CGAL::polynomial_subresultants_with_cofactors
	(f,g,
	 std::back_inserter(sres),
	 std::back_inserter(coP),
	 std::back_inserter(coQ));
      CGAL::set_pretty_mode(std::cout);
      for(int i=static_cast<int>(sres.size()-1);i>=0;i--) {
	assert(sres[i]==coP[i]*f + coQ[i]*g);
	assert(sres_check[i]==sres[i]);
      }
    }
   
    

    {
      Poly_int1 f=3*t*t*t*t+t*t*t-t+1;
      Poly_int1 g=5*t*t+7;
     
      std::vector<Poly_int1> sres,coP,coQ,sres_check;
      
      CGAL::internal::prs_polynomial_subresultants<Poly_int1_traits>
	(f,g,std::back_inserter(sres_check));      

      CGAL::polynomial_subresultants_with_cofactors
	(f,g,
	 std::back_inserter(sres),
	 std::back_inserter(coP),
	 std::back_inserter(coQ));
      for(int i = 0; i < static_cast<int>(sres.size());i++) {
	assert(sres_check[i]==sres[i]);
	assert(sres[i]==coP[i]*f + coQ[i]*g);
      }

    }
    { // Bug reported by Eric Berberich, 26.03.2010 -> Fixed

      Poly_int2 f=(CGAL::ipower(x,16) + 8*CGAL::ipower(x,14) + 28*CGAL::ipower(x,12) + 56*CGAL::ipower(x,10) + 70*CGAL::ipower(x,8) + 56*CGAL::ipower(x,6) + 28*CGAL::ipower(x,4) + 8*CGAL::ipower(x,2) + 1)*CGAL::ipower(y,16) + (8*CGAL::ipower(x,16) + 40*CGAL::ipower(x,14) + 72*CGAL::ipower(x,12) + 40*CGAL::ipower(x,10) + (-40)*CGAL::ipower(x,8) + (-72)*CGAL::ipower(x,6) + (-40)*CGAL::ipower(x,4) + (-8)*CGAL::ipower(x,2))*CGAL::ipower(y,14) + (28*CGAL::ipower(x,16) + 72*CGAL::ipower(x,14) + 40*CGAL::ipower(x,12) + (-24)*CGAL::ipower(x,10) + (-32)*CGAL::ipower(x,8) + (-40)*CGAL::ipower(x,6) + (-40)*CGAL::ipower(x,4) + (-8)*CGAL::ipower(x,2) + 4)*CGAL::ipower(y,12) + (56*CGAL::ipower(x,16) + 40*CGAL::ipower(x,14) + (-24)*CGAL::ipower(x,12) + 376*CGAL::ipower(x,10) + 328*CGAL::ipower(x,8) + (-424)*CGAL::ipower(x,6) + (-360)*CGAL::ipower(x,4) + 8*CGAL::ipower(x,2))*CGAL::ipower(y,10) + (70*CGAL::ipower(x,16) + (-40)*CGAL::ipower(x,14) + (-32)*CGAL::ipower(x,12) + 328*CGAL::ipower(x,10) + 1458*CGAL::ipower(x,8) + 320*CGAL::ipower(x,6) + 84*CGAL::ipower(x,4) + (-48)*CGAL::ipower(x,2) + 4)*CGAL::ipower(y,8) + (56*CGAL::ipower(x,16) + (-72)*CGAL::ipower(x,14) + (-40)*CGAL::ipower(x,12) + (-424)*CGAL::ipower(x,10) + 320*CGAL::ipower(x,8) + (-128)*CGAL::ipower(x,6) + 48*CGAL::ipower(x,4) + (-16)*CGAL::ipower(x,2))*CGAL::ipower(y,6) + (28*CGAL::ipower(x,16) + (-40)*CGAL::ipower(x,14) + (-40)*CGAL::ipower(x,12) + (-360)*CGAL::ipower(x,10) + 84*CGAL::ipower(x,8) + 48*CGAL::ipower(x,6) + 24*CGAL::ipower(x,4))*CGAL::ipower(y,4) + (8*CGAL::ipower(x,16) + (-8)*CGAL::ipower(x,14) + (-8)*CGAL::ipower(x,12) + 8*CGAL::ipower(x,10) + (-48)*CGAL::ipower(x,8) + (-16)*CGAL::ipower(x,6))*CGAL::ipower(y,2) + (CGAL::ipower(x,16) + 4*CGAL::ipower(x,12) + 4*CGAL::ipower(x,8));
     
      Poly_int2 g=(CGAL::ipower(x,2))*CGAL::ipower(y,9) + ((-3)*CGAL::ipower(x,2) + (-1)*x + (-1))*CGAL::ipower(y,8) + (2*CGAL::ipower(x,4) + (-1)*CGAL::ipower(x,3) + 2*CGAL::ipower(x,2) + 3*x + 3)*CGAL::ipower(y,7) + (CGAL::ipower(x,5) + (-6)*CGAL::ipower(x,4) + (-1)*CGAL::ipower(x,3) + CGAL::ipower(x,2) + (-1)*x + (-2))*CGAL::ipower(y,6) + (CGAL::ipower(x,6) + (-2)*CGAL::ipower(x,5) + CGAL::ipower(x,4) + 6*CGAL::ipower(x,3) + 2*CGAL::ipower(x,2) + (-3)*x + (-2))*CGAL::ipower(y,5) + (2*CGAL::ipower(x,7) + (-3)*CGAL::ipower(x,6) + 11*CGAL::ipower(x,4) + (-1)*CGAL::ipower(x,3) + (-2)*CGAL::ipower(x,2) + 2*x + 3)*CGAL::ipower(y,4) + ((-1)*CGAL::ipower(x,7) + CGAL::ipower(x,6) + 5*CGAL::ipower(x,5) + (-10)*CGAL::ipower(x,4) + (-11)*CGAL::ipower(x,3) + (-9)*CGAL::ipower(x,2) + (-1))*CGAL::ipower(y,3) + (CGAL::ipower(x,9) + (-3)*CGAL::ipower(x,7) + (-6)*CGAL::ipower(x,5) + (-3)*CGAL::ipower(x,4) + 8*CGAL::ipower(x,3) + 12*CGAL::ipower(x,2))*CGAL::ipower(y,2) + ((-1)*CGAL::ipower(x,8) + (-1)*CGAL::ipower(x,7) + 3*CGAL::ipower(x,6) + 4*CGAL::ipower(x,5) + 5*CGAL::ipower(x,4) + (-4)*CGAL::ipower(x,2))*y;
      std::vector<Poly_int2> sres,coP,coQ,sres_check;
      
      CGAL::internal::prs_polynomial_subresultants<Poly_int2_traits>
	(f,g,std::back_inserter(sres_check));      

      CGAL::polynomial_subresultants_with_cofactors
	(f,g,
	 std::back_inserter(sres),
	 std::back_inserter(coP),
	 std::back_inserter(coQ));
      for(int i = 0; i < static_cast<int>(sres.size());i++) {
	assert(sres[i]==sres_check[i]);
	assert(sres[i]==coP[i]*f + coQ[i]*g);
      }
    }

      



    { // bug reported by Eric Berberich
      Poly_int2 x = CGAL::shift(Poly_int2(1),1,0);
      Poly_int2 y = CGAL::shift(Poly_int2(1),1,1);
      Poly_int2 P = y*y + (x*x + (-2));
      Poly_int2 Q = (-1)*y*y + (x*x);
      std::vector<Poly_int2> res,coP,coQ;
      CGAL::polynomial_subresultants_with_cofactors(
          P,Q,
          std::back_inserter(res),
          std::back_inserter(coP),
          std::back_inserter(coQ)); 
    }

  return;
}


int main(){

#ifdef CGAL_USE_LEDA
    test_routine<CGAL::LEDA_arithmetic_kernel>();
#else
    std::cerr << "LEDA tests skipped!" << std::endl;
#endif
#ifdef CGAL_USE_CORE
    test_routine<CGAL::CORE_arithmetic_kernel>();
#else
    std::cerr << "CORE tests skipped!" << std::endl;
#endif
    return 0;
} 
