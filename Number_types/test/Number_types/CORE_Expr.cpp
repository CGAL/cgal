#include <iostream>
#include <CGAL/config.h>
#ifdef CGAL_USE_CORE
#include <sstream>
#include <list>
#include <cassert>
#include <cstdlib>

#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>

void precision_bug()
{
  typedef CORE::Expr FT;
  //computation with forced evaluation of sqrtD1D2
  {
    FT zero(0);
    FT ipx(0.3080532378), iqx(0.3080629044), irx(0.3080725711);
    FT a1(0.1282376364 - 0.1282279698);
    FT a2(0.1282473031 - 0.1282376364);
    FT b1 = ipx - iqx;
    FT b2 = iqx - irx;
    FT n1 = a1 * a1 + b1 * b1;
    FT n2 = a2 * a2 + b2 * b2;
    FT D1D2 = n1 * n2;
    FT sqrtD1D2 = CGAL::sqrt(D1D2);
    FT a1a2b1b2 = a1 * a2 + b1 * b2;
    FT uz =  sqrtD1D2 -  a1a2b1b2 ;

    sqrtD1D2.approx(53,1075); //force evaluation of sqrtD1D2
    uz.approx(53,1075);

    assert(!uz.isZero());
  }
  //computation without forced evaluation of sqrtD1D2
  {
    FT zero(0);
    FT ipx(0.3080532378), iqx(0.3080629044), irx(0.3080725711);
    FT a1(0.1282376364 - 0.1282279698);
    FT a2(0.1282473031 - 0.1282376364);
    FT b1 = ipx - iqx;
    FT b2 = iqx - irx;
    FT n1 = a1 * a1 + b1 * b1;
    FT n2 = a2 * a2 + b2 * b2;
    FT D1D2 = n1 * n2;
    FT sqrtD1D2 = CGAL::sqrt(D1D2);
    FT a1a2b1b2 = a1 * a2 + b1 * b2;
    FT uz =  sqrtD1D2 -  a1a2b1b2 ;

    assert(!uz.isZero());
  }
  std::cout << "precision bug OK\n";
}

void test_istream()
{
  std::list<std::string> strings;
  strings.push_back(std::string("0.231262"));
  strings.push_back(std::string("-15.123534563"));
  strings.push_back(std::string("12345"));
  std::list< std::pair<double,double> > intervals;
  intervals.push_back(std::make_pair(0.231261,0.231263));
  intervals.push_back(std::make_pair(-15.123534564,-15.123534562));
  intervals.push_back(std::make_pair(12345,12345));

  std::list< std::pair<double,double> >::iterator it_inter=intervals.begin();
  for (std::list<std::string>::iterator it=strings.begin();it!=strings.end();++it,++it_inter)
  {
    std::stringstream ss;
    ss << *it;

    CORE::Expr e;
    ss >> e;

    if( e < it_inter->first || e > it_inter->second ){
      std::cerr << "ERROR" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  std::cout << "test istream OK\n";
}

template <class FT>
void test_MSB_bug()
{
  FT px1(4.7320508075688767), py1(4.0000000000000000), pz1(1.9999999999999969);
  FT qx1(4.7320508075688767), qy1(4.0000000000000009), qz1(1.9999999999999880);
  FT rx1(5.1624261825907327), ry1(4.0000000000000000), rz1(2.0000000000000009);

  FT px2(4.7320508075688767), py2(4.0000000000000000), pz2(2.0000000000000044);
  FT qx2(4.2679491924311224), qy2(4.2679491924311233), qz2(2.0000000000000013);
  FT rx2(4.7320508075688767), ry2(4.0000000000000009), rz2(1.9999999999999880);

  FT rqy1 = qy1-ry1; // this operation will trigger a -infinity lMSB()
  FT rpx1 = px1-rx1;
  FT rpy1 = py1-ry1;
  FT rpz1 = pz1-rz1;
  FT rqx1 = qx1-rx1;
  FT rqz1 = qz1-rz1;

  FT b1 = rpz1*rqx1 - rqz1*rpx1;
  FT c1 = rpx1*rqy1 - rqx1*rpy1;

  FT rpx2 = px2-rx2;
  FT rpy2 = py2-ry2; // this operation will trigger a -infinity lMSB()
  FT rpz2 = pz2-rz2;
  FT rqx2 = qx2-rx2;
  FT rqy2 = qy2-ry2;
  FT rqz2 = qz2-rz2;
  FT b2 = rpz2*rqx2 - rqz2*rpx2;
  FT c2 = rpx2*rqy2 - rqx2*rpy2;
  FT res = b1*c2-c1*b2;

  assert(res != 0);
}

int main() {
    precision_bug();
    test_istream();
    test_MSB_bug<CORE::BigRat>();
    test_MSB_bug<CORE::Expr>();

    typedef CORE::Expr NT;
    typedef CGAL::Field_with_root_of_tag Tag;
    typedef CGAL::Tag_true Is_exact;
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(-15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6), NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(-15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(-15));

    CGAL::test_real_embeddable<NT>();

  return 0;
}

#if 0

#include <cstdlib>
#include <sstream>

template <class NT>
void test_io(const NT& x){
    NT tmp;
    std::ostringstream os;
    os << x;
    std::istringstream is(os.str());
    is >> tmp;
    assert( x == tmp );
}

/*
COEFF Coefficient type of Polynomial
REAL  a FieldWithSqrt
RATIONAL a numbertype representing the rational numbers
Z a numbertype representing Z (needed for Descartes)
*/
template <class COEFF,class REAL,class RATIONAL,class Z>
void algebraic_real_test()
{
    typedef COEFF Coeff_NT;
    typedef REAL real_NT;
    typedef RATIONAL rat_NT;
    typedef Z Integer;
    typedef typename CGAL::Coercion_traits< Coeff_NT, rat_NT>::Type Type;

    typedef NiX::Algebraic_real<Coeff_NT,real_NT,rat_NT> ALGNUM;
    typedef NiX::Polynomial<Coeff_NT> Poly;

    ::NiX::Residue::set_current_prime(29);
    typename ::NiX::NT_traits<real_NT>::Sqrt real_sqrt;

    // general test of comparable functionality
    NiX::test_real_comparable<ALGNUM>();
    // test of constructors
    Poly P_00(Coeff_NT(0));                   // zero polynomial
    Poly P_01(Coeff_NT(1));                   // constant polynomial
    Poly P_1(Coeff_NT(-1),Coeff_NT(1));       //(x-1)
    Poly P_2(Coeff_NT(-2),Coeff_NT(1));       //(x-2)
    Poly P_3(Coeff_NT(-3),Coeff_NT(1));       //(x-3)
    Poly P_4(Coeff_NT(-4),Coeff_NT(1));       //(x-4)
    Poly P_12=P_1*P_2;    //(x-1)(x-2)
    Poly P_123=P_1*P_2*P_3;    //(x-1)(x-2)(x-3)
    Poly P_s2(Coeff_NT(-2),Coeff_NT(0),Coeff_NT(1)); //(x^2-2)
    Poly P_s3(Coeff_NT(-3),Coeff_NT(0),Coeff_NT(1)); //(x^2-3)
    Poly P_s5(-Coeff_NT(5),Coeff_NT(0),Coeff_NT(1));
    Poly P_s10(-Coeff_NT(10),Coeff_NT(0),Coeff_NT(1));
    Poly P_s30(-Coeff_NT(30),Coeff_NT(0),Coeff_NT(1));
    Poly P_s25= P_s2*P_s5;
    Poly P_s2510= P_s2*P_s5*P_s10;
    Poly P_s530= P_s5 * P_s30;

    ALGNUM tmp;
    ALGNUM tmp1,tmp2;

    real_NT real, real1, real_2;

    rat_NT m;
    real_NT mm;
    // general constructors;
    // default
    // tmp = IS_RATIONAL = 0
    tmp = ALGNUM();
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()==0);
    // from int
    tmp = ALGNUM(1);
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()==1);

    tmp = ALGNUM(5);
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()==5);

    // from Field
    // from int
    tmp = ALGNUM(rat_NT(0));
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()==0);

    tmp = ALGNUM(rat_NT(1));
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()==1);

    tmp = ALGNUM(rat_NT(5)/ rat_NT(2));
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()== rat_NT(5)/ rat_NT(2));

    // general constructor
    // tmp = 1
    tmp = ALGNUM(P_1,-2,+2);
    if ((LiS::Compare_types< rat_NT, Type  >::same_type)) {
        NiX_test(tmp.is_rational());
        NiX_test(tmp.type()==NiX::IS_RATIONAL);
        NiX_test(tmp==rat_NT(1));
        NiX_test(tmp.rational()==1);
    } else {
        NiX_test(!tmp.is_rational());
        NiX_test(tmp.type()==NiX::IS_REAL);
        NiX_test(tmp.real()==1);
        NiX_test(tmp==rat_NT(1));
        NiX_test(tmp.type()==NiX::IS_RATIONAL);
        NiX_test(tmp.is_rational());
        NiX_test(tmp.rational()==1);
    }
    tmp = ALGNUM(P_1,1,1);
    if ((LiS::Compare_types< rat_NT, Type  >::same_type)) {
        NiX_test(tmp.is_rational());
        NiX_test(tmp.type()==NiX::IS_RATIONAL);
        NiX_test(tmp==rat_NT(1));
        NiX_test(tmp.rational()==1);
    } else {
        NiX_test(tmp.is_rational());
        NiX_test(tmp.type()==NiX::IS_RATIONAL);
        NiX_test(tmp.real()==1);
        NiX_test(tmp==rat_NT(1));
        NiX_test(tmp.is_rational());
        NiX_test(tmp.rational()==1);
    }

    // tmp IS_REAL == sqrt(2);
    tmp = ALGNUM(P_s2,1,2);
    NiX_test(tmp.type()==NiX::IS_REAL);
    NiX_test(tmp.real()==real_sqrt(real_NT(2)));

    // special constructors
    // from int
    tmp = ALGNUM(2);
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()==rat_NT(2));
    //from rat_NT
    tmp = ALGNUM(rat_NT(2));
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()==rat_NT(2));
    //from polynomial but rational
    tmp = ALGNUM(P_123,rat_NT(2),rat_NT(2));
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()==rat_NT(2));

    // member functions
    // tmp IS_GENERAL == 2;
    tmp = ALGNUM(P_123,rat_NT(3)/2,rat_NT(5)/2);
    NiX_test(tmp.type()==NiX::IS_GENERAL);
    NiX_test(tmp.polynomial()==P_123);
    NiX_test(tmp.low()==rat_NT(3)/2);
    NiX_test(tmp.high()==rat_NT(5)/2);
    NiX_test(tmp.sign_at_low()==P_123.sign_at(rat_NT(3)/2));

    //unary operator -
    // tmp IS_GENERAL == 2;
    tmp = - ALGNUM(P_123,rat_NT(3)/2,rat_NT(5)/2);
    NiX_test( tmp == ALGNUM(-2) );

    // refine
    tmp = ALGNUM(P_123,rat_NT(3)/2,rat_NT(5)/2);
    tmp.refine();
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()==rat_NT(2));
    // tmp IS_GENERAL = sqrt 2
    tmp = ALGNUM(P_s2*P_3,rat_NT(1),rat_NT(2));
    tmp.refine();
    NiX_test(tmp.low()==rat_NT(1));
    NiX_test(tmp.high()==rat_NT(3)/2);

    // strong_refine
    // tmp IS_GENERAL = sqrt 2
    tmp = ALGNUM(P_s2*P_3,rat_NT(1),rat_NT(2));
    m = rat_NT(3)/2;
    tmp.strong_refine(m);
    NiX_test(m < tmp.low() || tmp.high() < m);

    // tmp IS_GENERAL = sqrt 2
    tmp = ALGNUM(P_s2*P_3,rat_NT(1),rat_NT(2));
    mm = real_NT(3)/2;
    tmp.strong_refine(mm);
    NiX_test(tmp.low()!=mm);
    NiX_test(tmp.high()!=mm);
    mm = real_sqrt(real_NT(2));


    // refine_to(a,b)
    // tmp IS_GENERAL = sqrt 2
    tmp = ALGNUM(P_s2*P_4,rat_NT(0),rat_NT(3));
    NiX_test(tmp.type() == NiX::IS_GENERAL);
    tmp.refine_to(rat_NT(1), rat_NT(2));
    NiX_test(tmp.low()  >= rat_NT(1));
    NiX_test(tmp.high() <= rat_NT(2));

    // tmp IS_REAL = sqrt 2
    tmp = ALGNUM(P_s2,rat_NT(0),rat_NT(3));
    NiX_test(tmp.type() == NiX::IS_REAL);
    tmp.refine_to(rat_NT(1), rat_NT(2));
    NiX_test(tmp.low()  >= rat_NT(1));
    NiX_test(tmp.high() <= rat_NT(2));

    // compare(rat)
    // tmp IS_GENERAL = sqrt 2
    tmp = ALGNUM(P_s2*P_3,rat_NT(1),rat_NT(2));
    m = rat_NT(1);
    NiX_test(tmp.compare(m)==1);
    m = rat_NT(2);
    NiX_test(tmp.compare(m)==-1);
    mm = real_NT(1);
    NiX_test(tmp.compare(mm)==1);
    mm = real_NT(2);
    NiX_test(tmp.compare(mm)==-1);

    tmp1 = ALGNUM(P_3, rat_NT(1), rat_NT(4));
    tmp2 = ALGNUM(rat_NT(1));
    NiX_test(tmp1.compare(tmp2) == 1);
    NiX_test(tmp1.low() != rat_NT(1));

    tmp1 = ALGNUM(P_3, rat_NT(0), rat_NT(4));
    tmp2 = ALGNUM(rat_NT(1));
    NiX_test(tmp1.compare(tmp2) == 1);
    NiX_test(tmp1.low() != rat_NT(0));

    // tmp IS_GENERAL = 3
    tmp = ALGNUM(P_s2*P_3,rat_NT(2),rat_NT(4));
    m = rat_NT(3);
    NiX_test(tmp.compare(m)==0);
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.rational()==rat_NT(3));
    NiX_test(tmp.polynomial().degree() == 1);
    NiX_test(tmp.polynomial().evaluate(Coeff_NT(3)) == Coeff_NT(0));

    // tmp1/2 IS_GENERAL = 3
    tmp1 = ALGNUM(P_s2*P_3,rat_NT(2),rat_NT(4));
    tmp2 = ALGNUM(P_s2*P_3,rat_NT(2),rat_NT(4));
    NiX_test(tmp1.compare(tmp2)==0);

    tmp1 = ALGNUM(P_123,rat_NT(1)/2,rat_NT(3)/2);
    tmp2 = ALGNUM(P_s2*P_3,rat_NT(1),rat_NT(2));
    NiX_test(tmp1.compare(tmp2)==-1);

    // tmp1 GENERAL = sqrt(5)
    tmp1 = ALGNUM(P_s530,rat_NT(2),rat_NT(3));
    // tmp2 GENERAL = sqrt(5)
    tmp2 = ALGNUM(P_s2510,rat_NT(2),rat_NT(3));
    NiX_test(tmp1.compare(tmp2)==0);
    NiX_test(tmp1.type()==NiX::IS_REAL);
    NiX_test(tmp2.type()==NiX::IS_REAL);

    // compare (real)
    tmp1 = ALGNUM(P_s530,rat_NT(2),rat_NT(3));
    real1 = NiX::sqrt(real_NT(6));
    NiX_test(tmp1.compare(real1) == CGAL::SMALLER );
    real1 = NiX::sqrt(real_NT(5));
    NiX_test(tmp1.compare(real1) == CGAL::EQUAL );
    real1 = NiX::sqrt(real_NT(4));
    NiX_test(tmp1.compare(real1) == CGAL::LARGER );

    // compare_distinct()
    tmp1 = ALGNUM(P_s530, rat_NT(2), rat_NT(3)); // sqrt(5)  = 2.236...
    tmp2 = ALGNUM(P_s530, rat_NT(5), rat_NT(6)); // sqrt(30) = 5.477...
    NiX_test(tmp1.compare_distinct(tmp2) == CGAL::SMALLER);
    NiX_test(tmp2.compare_distinct(tmp1) == CGAL::LARGER);

    //member functions
    // is_root_of
    tmp1 = ALGNUM(P_s2510,rat_NT(1)/2,rat_NT(3)/2);
    NiX_test(tmp1.is_root_of(P_s530*P_s2));
    tmp1 = ALGNUM(P_s2510,rat_NT(1)/2,rat_NT(3)/2);
    NiX_test(!tmp1.is_root_of(P_s530));

    //rational_between
    {
        rat_NT r;
        tmp1 = ALGNUM(P_s2,rat_NT(1),rat_NT(2)); //sqrt2
        tmp2 = ALGNUM(P_s3,rat_NT(1),rat_NT(3)); //sqrt3
        r = tmp1.rational_between(tmp2);
        NiX_test(tmp1.compare(r)==CGAL::SMALLER);
        NiX_test(tmp2.compare(r)==CGAL::LARGER);

        r = tmp2.rational_between(tmp1);
        NiX_test(tmp1.compare(r)==CGAL::SMALLER);
        NiX_test(tmp2.compare(r)==CGAL::LARGER);
    }

    // to_double()
    tmp = ALGNUM(P_1*P_3*P_4, rat_NT(0), rat_NT(2));
    NiX_test(fabs(tmp.to_double() - 1.0) < 1e-10);
    tmp = ALGNUM(P_1*P_3, rat_NT(0), rat_NT(2));
    NiX_test(fabs(tmp.to_double() - 1.0) < 1e-10);
    tmp = ALGNUM(P_1, rat_NT(0), rat_NT(2));
    NiX_test(fabs(tmp.to_double() - 1.0) < 1e-10);

    // input/output
    //test rational input
    {
        tmp = ALGNUM(2);
        std::ostringstream os;
        os << tmp ;
        std::istringstream is(os.str());
        is >> tmp2;
        NiX_test_msg( ALGNUM(2) == tmp2, "IO_TEST failed");
    }
    //test general input
    {
        ALGNUM a1,b1,c1,a2,b2,c2;
        a1 = ALGNUM(P_s2*P_3,rat_NT(1),rat_NT(2));
        b1 = ALGNUM(rat_NT(2));
        c1 = ALGNUM(P_s2,-2,-1);

        test_io(a1);
        test_io(b1);
        test_io(c1);

        std::ostringstream os;
        os << a1 << b1 << c1;
        std::istringstream is(os.str());
        is>> a2 >> b2 >> c2;
        NiX_test_msg( a1 == a2, "IO_TEST failed");
        NiX_test_msg( b1 == b2, "IO_TEST failed");
        NiX_test_msg( c1 == c2, "IO_TEST failed");
    }
    // output benchmark (code coverage)
    {
        ALGNUM a1,b1,c1,a2,b2,c2;
        a1 = ALGNUM(P_s2*P_3,rat_NT(1),rat_NT(2));
        b1 = ALGNUM(rat_NT(2));
        c1 = ALGNUM(P_s2,-2,-1);

        test_io(a1);
        test_io(b1);
        test_io(c1);

        std::ostringstream os;
        LiS::set_benchmark_mode(os);
        os << a1 << b1 << c1;
    }


    // test for Handle with union
    {
        typedef
            NiX::Algebraic_real
            <Coeff_NT,real_NT,rat_NT,::LiS::Handle_with_union> Int;
        Int i(5);
        Int j(5);
        Int k(6);
        NiX_test( ! i.identical( j));
        NiX_test( ! i.identical( k));
        NiX_test( ! j.identical( k));
        NiX_test( i == j);
        NiX_test( ! (i == k));
        NiX_test( i.identical( j));
        NiX_test( ! i.identical( k));
        NiX_test( ! j.identical( k));
    }
    // test for Handle without union
    {
        typedef
            NiX::Algebraic_real
            <Coeff_NT,real_NT,rat_NT,::LiS::Handle_without_union> Int;
        Int i(5);
        Int j(5);
        Int k(6);
        NiX_test( ! i.identical( j));
        NiX_test( ! i.identical( k));
        NiX_test( ! j.identical( k));
        NiX_test( i == j);
        NiX_test( ! (i == k));
        NiX_test( ! i.identical( j));
        NiX_test( ! i.identical( k));
        NiX_test( ! j.identical( k));
    }


    // static member function conjugate
    {
        std::vector<ALGNUM> roots;
        // for the empty case
        ALGNUM::conjugate(roots.begin(),roots.end());

        roots.push_back(ALGNUM(P_s2510,1,2));
        roots.push_back(ALGNUM(P_s2510,-2,-1));
        roots.push_back(ALGNUM(P_s2510,2,3));
        roots.push_back(ALGNUM(P_s2510,-3,-2));
        roots.push_back(ALGNUM(P_s2510,3,4));
        roots.push_back(ALGNUM(P_s2510,-4,3));
        ALGNUM::conjugate(roots.begin(),roots.end());
        NiX_test(roots[0]==ALGNUM(P_s25,1,2));
        NiX_test(roots[0].polynomial()==P_s25);
        NiX_test(roots[1].polynomial()==P_s25);
        NiX_test(roots[2].polynomial()==P_s25);
        NiX_test(roots[3].polynomial()==P_s25);
        NiX_test(roots[4].polynomial()==P_s10);
        NiX_test(roots[5].polynomial()==P_s10);

    }
    // to_Interval
    {
        ALGNUM TMP;
        typename NiX::NT_traits<ALGNUM>::To_Interval to_Interval;
        NiX_test(NiX::in(25.0,to_Interval(ALGNUM(25))));
        NiX_test(NiX::in(sqrt(2),to_Interval(ALGNUM(P_s2,1,2))));
        NiX_test(NiX::in(sqrt(2),to_Interval(ALGNUM(P_s2510,1,2))));
        NiX_test(NiX::in(-sqrt(2),to_Interval(ALGNUM(P_s2510,-2,-1))));
        NiX_test(NiX::in(sqrt(5),to_Interval(ALGNUM(P_s2510,2,3))));
        NiX_test(NiX::in(-sqrt(5),to_Interval(ALGNUM(P_s2510,-3,-2))));
        NiX_test(NiX::in(sqrt(10),to_Interval(ALGNUM(P_s2510,3,4))));
        NiX_test(NiX::in(-sqrt(10),to_Interval(ALGNUM(P_s2510,-4,-3))));
    }
    //simplify
    {
        // just a syntax check
        ALGNUM(P_s2510,1,2).simplify();
    }
}



/*
COEFF Coefficient type of Polynomial
REAL  a FieldWithSqrt
RATIONAL a numbertype representing the rational numbers
Z a numbertype representing Z (needed for Descartes)
*/
template <class COEFF,class REAL,class RATIONAL,class Z>
void algebraic_real_test_for_set_rational(bool set_rational)
{
    typedef COEFF Coeff_NT;
    typedef REAL real_NT;
    typedef RATIONAL rat_NT;
    typedef Z Integer;
    typedef typename NiX::Coercion_traits< Coeff_NT, rat_NT>::Type Type;

    typedef NiX::Algebraic_real<Coeff_NT,real_NT,rat_NT> ALGNUM;
    typedef NiX::Polynomial<Coeff_NT> Poly;

    ::NiX::Residue::set_current_prime(29);

    // general test of comparable functionality
    NiX::test_real_comparable<ALGNUM>();
    // test of constructors
    Poly P_1(Coeff_NT(-1),Coeff_NT(1));       //(x-1)
    Poly P_3(Coeff_NT(-3),Coeff_NT(1));       //(x-3)
    Poly P_13=P_1*P_3;    //(x-1)(x-3)

    ALGNUM tmp;

    // general constructor for linear polynomial
    // tmp = 1
    tmp = ALGNUM(P_1,-2,+2);
    if (set_rational) {
        NiX_test(tmp.is_rational());
        NiX_test(tmp.type()==NiX::IS_RATIONAL);
        NiX_test(tmp==rat_NT(1));
        NiX_test(tmp.rational()==1);
    } else {
        NiX_test(!tmp.is_rational());
        NiX_test(tmp.type()==NiX::IS_REAL);
        NiX_test(tmp.real()==1);
        NiX_test(tmp==rat_NT(1));
        NiX_test(tmp.type()==NiX::IS_RATIONAL);
        NiX_test(tmp.is_rational());
        NiX_test(tmp.rational()==1);
    }

    tmp = ALGNUM(P_13,0,2);
    NiX_test(!tmp.is_rational());
    NiX_test(tmp.type()==NiX::IS_REAL);
    NiX_test(tmp.real()==1);
    NiX_test(tmp==rat_NT(1));
    NiX_test(tmp.type()==NiX::IS_RATIONAL);
    NiX_test(tmp.is_rational());
    NiX_test(tmp.rational()==1);
}

template <class AT>
void algebraic_real_test_at(){
/*
COEFF Coefficient type of Polynomial
REAL  a FieldWithSqrt
RATIONAL a numbertype representing the rational numbers
Z a numbertype representing Z (needed for Descartes)
template< class COEFF, class REAL, class RATIONAL, class Z>
algebraic_real_test()
*/


    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;
    typedef typename AT::Field_with_sqrt Real;
    {
        algebraic_real_test <Integer, Real, Rational,Integer>();
        algebraic_real_test_for_set_rational
            <Integer, Real, Rational, Integer> (true);
    }{
        algebraic_real_test <Rational, Real, Rational,Integer>();
        algebraic_real_test_for_set_rational
            <Rational, Real, Rational, Integer>(true);
    }{
        typedef Rational ROOT_NT;
        typedef NiX::Sqrt_extension<Rational,ROOT_NT> EXT;
        algebraic_real_test<EXT, Real, Rational, Integer>();
        algebraic_real_test_for_set_rational
            <EXT, Real , Rational, Integer>(false);
    }{
        typedef Integer ROOT_NT;
        typedef NiX::Sqrt_extension<Integer,ROOT_NT > EXT;
        algebraic_real_test<EXT, Real, Rational, Integer>();
        algebraic_real_test_for_set_rational
            <EXT, Real , Rational, Integer>(false);
    }
}


int main(){
#ifdef CGAL_USE_LEDA
     algebraic_real_test_at<NiX::LEDA_arithmetic_traits>();
#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_CORE
     algebraic_real_test_at<NiX::CORE_arithmetic_traits>();
#endif // CGAL_USE_CORE
     return 0;
}
#endif

#else // CGAL_USE_CORE
int main() { return 0; }
#endif // CGAL_USE_CORE
//EOF
