
#include <CGAL/use.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/convert_to_bfi.h>
#include <CGAL/Bigfloat_interval_traits.h>

#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Test/_test_coercion_traits.h>

#include <cstdlib>
#include <sstream>
#include <iostream>
#include <cassert>

#include <CGAL/disable_warnings.h>

// TODO: Included from EXACUS/NumeriX/include/NiX/number_type_utils.h
namespace CGAL {
template <class NT , class RT>
inline
void convert_to(const NT& x, RT& r){
    typedef CGAL::Coercion_traits<NT,RT> CT;
    typedef typename CT::Type Type;
    CGAL_USE_TYPE(Type);
    static_assert(::std::is_same<Type,RT>::value);
    r = typename CT::Cast()(x);
}
} //namespace CGAL
// end included from number_type_utils.h

template <class NT>
void test_io(const NT& x){
    NT tmp;
    std::ostringstream os;
    os << x;
    std::istringstream is(os.str());
    is >> tmp;
    assert( x == tmp );
}

template <class EXT>
void constructor_test(){
    typedef typename EXT::NT NT;
    typedef typename EXT::ROOT ROOT;
    assert(EXT()==EXT((NT)0));
    assert(EXT()==EXT((NT)0,(NT)0,(ROOT)7));
    assert(EXT().is_extended()==false);
    // from int
    assert(EXT( 1)==EXT((NT) 1,(NT)0,(ROOT)7));
    assert(EXT(-1)==EXT((NT)-1,(NT)0,(ROOT)7));
    assert(EXT( 7)==EXT((NT) 7,(NT)0,(ROOT)7));
    assert(EXT(-7)==EXT((NT)-7,(NT)0,(ROOT)7));
    assert(EXT( 1)+EXT(-1)==EXT());
    assert(EXT(1).is_extended()==false);
    // from NT
    assert(EXT((NT) 1)==EXT((NT) 1,(NT)0,(ROOT)7));
    assert(EXT((NT)-1)==EXT((NT)-1,(NT)0,(ROOT)7));
    assert(EXT((NT) 7)==EXT((NT) 7,(NT)0,(ROOT)7));
    assert(EXT((NT)-7)==EXT((NT)-7,(NT)0,(ROOT)7));
    assert(EXT((NT) 1)+EXT((NT)-1)==EXT());
    assert(EXT((NT)1).is_extended()==false);
    // general constructor
    const EXT x(NT(3),NT(4),ROOT(17));
    assert(x.a0()==NT(3));
    assert(x.a1()==NT(4));
    assert(x.root()==ROOT(17));
    assert(x.is_extended()==true);
}

template <class EXT>
void other_functions_test(){
    typedef typename EXT::NT NT;
    typedef typename EXT::ROOT ROOT;

    // is_zero
    assert( EXT(NT( 0),NT( 0),ROOT(17)).is_zero());
    assert(!EXT(NT( 3),NT( 0),ROOT(17)).is_zero());
    assert(!EXT(NT(-3),NT( 0),ROOT(17)).is_zero());
    assert(!EXT(NT( 0),NT( 3),ROOT(17)).is_zero());
    assert(!EXT(NT( 0),NT(-3),ROOT(17)).is_zero());
    assert(!EXT(NT(-3),NT(-3),ROOT(17)).is_zero());
    assert( EXT(NT( 0)).is_zero());
    assert(!EXT(NT(-5)).is_zero());
    assert(!EXT(NT( 5)).is_zero());
    assert( EXT(NT( 0),NT( 0),ROOT(4)).is_zero());
    assert(!EXT(NT( 6),NT( 3),ROOT(4)).is_zero());
    assert( EXT(NT( 6),NT(-3),ROOT(4)).is_zero());

    // abs
    {
        EXT x  = EXT(NT(-3),NT(-3),ROOT(17));
        EXT ax = EXT(NT( 3),NT( 3),ROOT(17));
        assert(x.abs() == ax);
    }{
        EXT x  = EXT(NT( 3),NT(-3),ROOT(17));
        EXT ax = EXT(NT(-3),NT( 3),ROOT(17));
        assert(x.abs() == ax);
    }{
        EXT x  = EXT(NT(-3),NT( 3),ROOT(17));
        EXT ax = EXT(NT(-3),NT( 3),ROOT(17));
        assert(x.abs() == ax);
    }{
        EXT x  = EXT(NT( 3),NT( 3),ROOT(17));
        EXT ax = EXT(NT( 3),NT( 3),ROOT(17));
        assert(x.abs() == ax);
    }
    // division / integral division in case root is a square
    {
        EXT x(4,6,4);
        EXT y(2,1,4);
        assert(x/y ==EXT(4));
    }{
        EXT x(4);
        EXT y(2,1,4);
        assert(x/y ==EXT(1));
    }
}

template <class EXT>
void io_test(){
    test_io(EXT(0,0,5));
    test_io(EXT(0,2,5));
    test_io(EXT(4,0,5));
    test_io(EXT(-4,2,5));
    test_io(EXT(4,-2,5));
    test_io(EXT(-4,-2,5));
    test_io(EXT(5));

    std::ostringstream os;
    EXT(0).output_maple(os); // no real test just to instantiate code
}

template<class NT, class ROOT,class REAL,class ACDE_TAG>
void convert_to_real(){
  typedef CGAL::Sqrt_extension<NT,ROOT,ACDE_TAG> EXT1;

    {
        EXT1 ext((NT)3 ,(NT)5  ,(ROOT)17);
        REAL real= REAL(3)+REAL(5)*CGAL_NTS sqrt(REAL(17));
        REAL tmp;
        CGAL::convert_to(ext,tmp);
        assert(tmp==real);
    }{
        EXT1 ext((NT)3);
        REAL real= REAL(3);
        REAL tmp;
        CGAL::convert_to(ext,tmp);
        assert(tmp==real);
    }
}

template<class NT, class ROOT,class BFI,class ACDE_TAG>
void convert_to_bfi(){
  typedef CGAL::Sqrt_extension<NT,ROOT,ACDE_TAG> EXT1;
    typename CGAL::Bigfloat_interval_traits<BFI>::Get_precision get_precision;
    typename CGAL::Bigfloat_interval_traits<BFI>::Set_precision set_precision;
    long old_precision = get_precision();

    for(int i = 0; i < 3; i++){
      long precision = old_precision;
      for(int p = 0; p < 3; p++){
        set_precision(precision);
        {
          EXT1 ext((NT)3 ,(NT)5  ,(ROOT)17);
          BFI interval= BFI(3)+BFI(5)*CGAL_NTS sqrt(BFI(17));
          BFI tmp = CGAL::convert_to_bfi(ext);
          assert(CGAL::overlap(tmp,interval));
        }{
          EXT1 ext((NT)3);
          BFI interval= BFI(3);
          BFI tmp = CGAL::convert_to_bfi(ext);
          assert(CGAL::overlap(tmp,interval));
        }
        precision*=2;
      }
    }
    set_precision(old_precision);
}


template < class AT, class ACDE_TAG>
void sqrt_ext_pretty_output_test(){
    typedef typename AT::Integer Integer;
    typedef CGAL::Sqrt_extension<Integer,Integer>  EXT1;
    typedef CGAL::Sqrt_extension<EXT1,Integer> EXT2;
    CGAL_USE_TYPE(EXT2);
    // TEST without Parens_as_product_tag
    {
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << EXT1(-1,0,2);
        assert( ss.str()=="-1");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << EXT1(0,0,2);
        assert( ss.str() == "0");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << EXT1(1,0,2);
        assert( ss.str() == "1");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << EXT1(1,2,2);
        assert( ss.str() == "1+2*sqrt(2)");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << EXT1(1,-2,2);
        assert( ss.str() == "1+(-2)*sqrt(2)");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << EXT1(-1,2,2);
        assert( ss.str() == "-1+2*sqrt(2)");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << EXT1(-1,-2,2);
        assert( ss.str() == "-1+(-2)*sqrt(2)");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << EXT1(0,2,2);
        assert( ss.str()=="2*sqrt(2)");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << EXT1(0,-2,2);
        assert( ss.str()=="(-2)*sqrt(2)");
    }
// TEST with Parens_as_product_tag
    {
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(EXT1(-1,0,2),CGAL::Parens_as_product_tag());
        assert( ss.str()=="(-1)");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(EXT1(0,0,2),CGAL::Parens_as_product_tag());
        assert( ss.str() == "0");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(EXT1(1,0,2),CGAL::Parens_as_product_tag());
        assert( ss.str() == "1");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(EXT1(1,2,2),CGAL::Parens_as_product_tag());
        assert( ss.str() == "(1+2*sqrt(2))");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(EXT1(1,-2,2),CGAL::Parens_as_product_tag());
        assert( ss.str() == "(1+(-2)*sqrt(2))");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(EXT1(-1,2,2),CGAL::Parens_as_product_tag());
        assert( ss.str() == "(-1+2*sqrt(2))");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(EXT1(-1,-2,2),CGAL::Parens_as_product_tag());
        assert( ss.str() == "(-1+(-2)*sqrt(2))");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(EXT1(0,2,2),CGAL::Parens_as_product_tag());
        assert( ss.str()=="2*sqrt(2)");
    }{
        std::stringstream ss;
        CGAL::IO::set_pretty_mode(ss);
        ss << CGAL::IO::oformat(EXT1(0,-2,2),CGAL::Parens_as_product_tag());
        assert( ss.str()=="(-2)*sqrt(2)");
    }
}


template <class EXT>
void to_double_test(){
    typedef typename EXT::NT NT;
    typedef typename EXT::ROOT ROOT;
    typedef CGAL::Real_embeddable_traits<EXT> NTT;

    typename NTT::To_double to_double;
    typename NTT::To_interval to_interval;

    {
        NT a0(1);
        NT a1(2);
        ROOT root(4);
        EXT ext(a0,a1,root);

        assert(to_double(ext)==5.0);
        assert(to_interval(ext)==CGAL_NTS to_interval(5.0) );
    }{
        NT a0(1);
        NT a1(-2);
        ROOT root(4);
        EXT ext(a0,a1,root);

        assert(to_double(ext)==-3.0);
        assert(to_interval(ext)==CGAL_NTS to_interval(-3.0) );
    }
}

//This test is dedicated to the comparison of numbers from different extensions
template <class EXT>
void test_compare(){
  typedef typename EXT::NT NT;
  typedef typename EXT::ROOT ROOT;
  typedef typename EXT::ACDE_TAG ACDE_TAG;

  EXT a(NT(1),NT(-5),ROOT(7));
  EXT b(NT(2),NT(8),ROOT(5));
  EXT c(NT(3),NT(0),ROOT(3));

  if(ACDE_TAG::value){
    assert(a.compare(b)==CGAL::SMALLER);
    assert(a.compare(c)==CGAL::SMALLER);
    assert(b.compare(c)==CGAL::LARGER );
    assert(a.compare(a)==CGAL::EQUAL);
    assert(b.compare(b)==CGAL::EQUAL);
    assert(c.compare(c)==CGAL::EQUAL);
  }
  assert(a.compare(a)==CGAL::EQUAL);
  assert(b.compare(b)==CGAL::EQUAL);
  assert(c.compare(c)==CGAL::EQUAL);


  assert(a.compare(b,false)==CGAL::SMALLER);
  assert(a.compare(c,false)==CGAL::SMALLER);
  assert(b.compare(c,false)==CGAL::LARGER );
  assert(a.compare(a,false)==CGAL::EQUAL);
  assert(b.compare(b,false)==CGAL::EQUAL);
  assert(c.compare(c,false)==CGAL::EQUAL);
  assert(a.compare(a,true)==CGAL::EQUAL);
  assert(b.compare(b,true)==CGAL::EQUAL);
  assert(c.compare(c,true)==CGAL::EQUAL);

  assert(CGAL::compare(a,b,false)==CGAL::SMALLER);
  assert(CGAL::compare(a,c,false)==CGAL::SMALLER);
  assert(CGAL::compare(b,c,false)==CGAL::LARGER );
  assert(CGAL::compare(a,a,false)==CGAL::EQUAL);
  assert(CGAL::compare(b,b,false)==CGAL::EQUAL);
  assert(CGAL::compare(c,c,false)==CGAL::EQUAL);
  assert(CGAL::compare(a,a,true)==CGAL::EQUAL);
  assert(CGAL::compare(b,b,true)==CGAL::EQUAL);
  assert(CGAL::compare(c,c,true)==CGAL::EQUAL);
}

template <class NT, class ROOT, class Algebraic_type,class ACDE_TAG>
void general_test(){
    typedef typename CGAL::Algebraic_structure_traits<NT>::Is_exact Is_exact;

    typedef CGAL::Sqrt_extension<NT,ROOT,ACDE_TAG> EXT1;
    EXT1 a1((NT)3 ,(NT)5  ,(ROOT)17);
    EXT1 b1((NT)-2,(NT)3  ,(ROOT)17);
    EXT1 c1((NT)7 ,(NT)-11,(ROOT)17);
    CGAL::test_algebraic_structure<EXT1,Algebraic_type,Is_exact>(a1,b1,c1);
    CGAL::test_real_embeddable<EXT1>();
    CGAL::test_implicit_interoperable<int,EXT1,EXT1>();
    CGAL::test_implicit_interoperable<NT ,EXT1,EXT1>();

    typedef CGAL::Sqrt_extension<EXT1,ROOT> EXT2;
    EXT2 a2(b1 ,a1  ,(ROOT)19);
    EXT2 b2(c1,(EXT1)0 ,(ROOT)19);
    EXT2 c2((EXT1)0 ,(EXT1)-11,(ROOT)19);
    CGAL::test_algebraic_structure<EXT2,Algebraic_type,Is_exact>(a2,b2,c2);
    CGAL::test_real_embeddable<EXT2>();

    typedef CGAL::Sqrt_extension<EXT2,ROOT> EXT3;
    EXT3 a3(b2 ,a2  ,(ROOT)23);
    EXT3 b3(c2,(EXT2)0 ,(ROOT)23);
    EXT3 c3((EXT2)0 ,(EXT2)-11,(ROOT)23);
    CGAL::test_algebraic_structure<EXT3,Algebraic_type,Is_exact>(a3,b3,c3);
    CGAL::test_real_embeddable<EXT3>();

    constructor_test<EXT1>();
    constructor_test<EXT2>();
    constructor_test<EXT3>();

    other_functions_test<EXT1>();
    other_functions_test<EXT2>();
    other_functions_test<EXT3>();

    io_test<EXT1>();
    io_test<EXT2>();
    io_test<EXT3>();

    to_double_test<EXT1>();

    test_compare<EXT1>();
    test_compare<EXT2>();
    test_compare<EXT3>();
}



template <class AT, class ACDE_TAG>
void fraction_traits_test(){
    //TEST Type traits ROOT of type INT
    typedef typename AT::Integer  INT;
    typedef typename AT::Rational RAT;

    typedef CGAL::Sqrt_extension<RAT     ,INT> RAT1_EXT;
    typedef CGAL::Sqrt_extension<INT     ,INT> INT1_EXT;
    typedef CGAL::Sqrt_extension<RAT1_EXT,INT> RAT2_EXT;
    typedef CGAL::Sqrt_extension<INT1_EXT,INT> INT2_EXT;
    typedef CGAL::Sqrt_extension<RAT2_EXT,INT> RAT3_EXT;
    typedef CGAL::Sqrt_extension<INT2_EXT,INT> INT3_EXT;

    typedef CGAL::Fraction_traits<RAT1_EXT> RAT_FT1;
    typedef CGAL::Fraction_traits<RAT2_EXT> RAT_FT2;
    typedef CGAL::Fraction_traits<RAT3_EXT> RAT_FT3;

    typedef CGAL::Fraction_traits<INT1_EXT> INT_FT1;
    typedef CGAL::Fraction_traits<INT2_EXT> INT_FT2;
    typedef CGAL::Fraction_traits<INT3_EXT> INT_FT3;
    CGAL_USE_TYPE(INT_FT3);

    // RAT_FTs decomposable
    assert((std::is_same< typename RAT_FT1::Is_fraction,
                    CGAL::Tag_true>::value));
    assert((std::is_same< typename RAT_FT2::Is_fraction,
                    CGAL::Tag_true>::value));
    assert((std::is_same< typename RAT_FT3::Is_fraction,
                    CGAL::Tag_true>::value));
    // RAT_FTi Numerator_type == INTi_EXT
    assert((std::is_same< typename RAT_FT1::Numerator_type,
                    INT1_EXT>::value));
    assert((std::is_same< typename RAT_FT2::Numerator_type,
                    INT2_EXT>::value));
    assert((std::is_same< typename RAT_FT3::Numerator_type,
                    INT3_EXT>::value));
    // RAT_FTi Denomiantor == INT
    assert((std::is_same< typename RAT_FT1::Denominator_type,
                    INT>::value));
    assert((std::is_same< typename RAT_FT2::Denominator_type,
                    INT>::value));
    assert((std::is_same< typename RAT_FT3::Denominator_type,
                    INT>::value));

    // INT_FTs not decomposable
    assert((std::is_same< typename INT_FT1::Is_fraction,
                    CGAL::Tag_false>::value));
    assert((std::is_same< typename INT_FT2::Is_fraction,
                    CGAL::Tag_false>::value));
    assert((std::is_same< typename INT_FT2::Is_fraction,
                    CGAL::Tag_false>::value));
    {
        typedef CGAL::Sqrt_extension<RAT,RAT>   RAT_RAT_EXT;
        typedef CGAL::Fraction_traits<RAT_RAT_EXT> RAT_RAT_FT;
        assert((std::is_same< typename RAT_RAT_FT::Is_fraction,
                        CGAL::Tag_false>::value));
    }{
        typedef CGAL::Sqrt_extension<INT1_EXT,INT1_EXT> INT_nEXT;
        typedef CGAL::Fraction_traits<INT_nEXT> INT_nEXT_FT;
        assert((std::is_same< typename INT_nEXT_FT::Is_fraction,
                        CGAL::Tag_false>::value));
    }

    {
        // Semantic for RAT1_EXT
        RAT1_EXT a(RAT(3)/RAT(10),RAT(4)/RAT(15),INT(5));
        typename RAT_FT1::Decompose decompose;
        typename RAT_FT1::Numerator_type  num;
        typename RAT_FT1::Denominator_type den;
        decompose(a,num,den);
        assert(num==typename RAT_FT1::Numerator_type(INT(9),INT(8),INT(5)));
        assert(den==typename RAT_FT1::Denominator_type(INT(30)));
        typename RAT_FT1::Compose compose;
        assert(a==compose(num,den));
        {
            RAT1_EXT a(RAT(3)/RAT(10));
            decompose(a,num,den);
            assert(num==typename RAT_FT1::Numerator_type(INT(3)));
            assert(den==typename RAT_FT1::Denominator_type(INT(10)));
            assert(!num.is_extended());
        }
    }
    {
        // Semantic RAT2_EXT
        RAT1_EXT a0(RAT(3)/RAT(10),RAT(4)/RAT(15),INT(5));
        RAT1_EXT a1(RAT(7)/RAT(3),RAT(5)/RAT(7),INT(5));
        RAT2_EXT a(a0,a1,INT(7));
        typename RAT_FT2::Decompose decompose;
        typename RAT_FT2::Numerator_type   num;
        typename RAT_FT2::Denominator_type den;
        decompose(a,num,den);
        assert(num==INT2_EXT(INT1_EXT(INT(63),INT(56),INT(5)),
                        INT1_EXT(INT(490),INT(150),INT(5)),
                        INT(7)));
        assert(den==typename RAT_FT2::Denominator_type(INT(210)));
        typename RAT_FT2::Compose compose;
        assert(a==compose(num,den));
    }


}

// this includes Scalar_factor_traits
// this includes Algebraic_extension_traits

template <class AT, class ACDE_TAG>
void scalar_factor_traits_test(){
    {
        typedef typename AT::Integer Integer;
        typedef CGAL::Sqrt_extension<Integer,Integer> EXT1;
        typedef CGAL::Sqrt_extension<EXT1,Integer>    EXT2;

        {
            typedef CGAL::Scalar_factor_traits<EXT1> SFT;
            typedef typename SFT::Scalar Scalar;
            typename SFT::Scalar_factor sfac;
            typename SFT::Scalar_div    sdiv;

            assert( (std::is_same<Scalar, Integer>::value) );

            assert((sfac(EXT1(0)))==Integer(0));
            assert((sfac(EXT1(3)))==Integer(3));
            assert((sfac(EXT1(0,3,2)))==Integer(3));
            assert((sfac(EXT1(6,9,2)))==Integer(3));
            EXT1 tmp;
            tmp = EXT1(3); sdiv(tmp,sfac(tmp));  assert(tmp==EXT1(1));
            tmp = EXT1(0,3,2); sdiv(tmp,sfac(tmp));  assert(tmp==EXT1(0,1,2));
            tmp = EXT1(6,9,2); sdiv(tmp,sfac(tmp));  assert(tmp==EXT1(2,3,2));
        }
        {
            typedef CGAL::Scalar_factor_traits<EXT2> SFT;
            typedef typename SFT::Scalar Scalar;
            typename SFT::Scalar_factor sfac;
            typename SFT::Scalar_div    sdiv;

            assert( (std::is_same<Scalar, Integer>::value) );

            assert((sfac(EXT2(0)))==Integer(0));
            assert((sfac(EXT2(3)))==Integer(3));
            assert((sfac(EXT2(EXT1(0),EXT1(3),2)))==Integer(3));
            assert((sfac(EXT2(EXT1(6),EXT1(9),2)))==Integer(3));
            EXT2 tmp;
            tmp = EXT2(3); sdiv(tmp,sfac(tmp));  assert(tmp==EXT2(1));
            tmp = EXT2(0,3,2); sdiv(tmp,sfac(tmp));  assert(tmp==EXT2(0,1,2));
            tmp = EXT2(6,9,2); sdiv(tmp,sfac(tmp));  assert(tmp==EXT2(2,3,2));
        }
        {
            typedef CGAL::Scalar_factor_traits<EXT1> SFT;
            typedef typename SFT::Scalar Scalar;
            typename SFT::Scalar_factor sfac;

            assert( (std::is_same<Scalar, Integer>::value) );

            assert((sfac(EXT1(0)))                 ==Integer(0));
            assert((sfac(EXT1(9),Integer(15)))     ==Integer(3));
            assert((sfac(EXT1(0,9,2),Integer(15))) ==Integer(3));
            assert((sfac(EXT1(18,9,2),Integer(15)))==Integer(3));
        }
    }
}

template <class AT,class ACDE_TAG>
void test_algebraic_extension_traits(){
    typedef typename AT::Integer  INT;
    typedef typename AT::Rational RAT;

    typedef CGAL::Sqrt_extension<RAT     ,INT,ACDE_TAG> RAT1_EXT;
    typedef CGAL::Sqrt_extension<INT     ,INT,ACDE_TAG> INT1_EXT;
    typedef CGAL::Sqrt_extension<RAT1_EXT,INT,ACDE_TAG> RAT2_EXT;
    typedef CGAL::Sqrt_extension<INT1_EXT,INT,ACDE_TAG> INT2_EXT;

    // normalisation factor
    typedef CGAL::Algebraic_extension_traits<RAT1_EXT> RAT1_EXT_ANT;
    typedef CGAL::Algebraic_extension_traits<INT1_EXT> INT1_EXT_ANT;
    typedef CGAL::Algebraic_extension_traits<RAT2_EXT> RAT2_EXT_ANT;
    typedef CGAL::Algebraic_extension_traits<INT2_EXT> INT2_EXT_ANT;

    {
        INT1_EXT a(INT(3),INT(5),INT(7));
        typename INT1_EXT_ANT::Normalization_factor normalization_factor;
        INT1_EXT r=normalization_factor(a);
        assert(r.a0()    ==INT( 3));
        assert(r.a1()    ==INT(-5));
        assert(r.root()  ==INT( 7));
        assert((r*a).a1()==INT( 0));
    }{
        RAT1_EXT a(RAT(3)/RAT(5),RAT(5)/RAT(11),INT(7));
        typename RAT1_EXT_ANT::Normalization_factor normalization_factor;
        RAT1_EXT r=normalization_factor(a);
        assert((r*a).a1()==RAT( 0));
    }{
        INT1_EXT a0(4,3,7);
        INT1_EXT a1(8,9,7);
        INT      root(5);
        INT2_EXT a(a0,a1,root);
        typename INT2_EXT_ANT::Normalization_factor normalization_factor;
        INT2_EXT r=normalization_factor(a);

        assert((r*a).a1()==INT1_EXT(0));
        assert((r*a).a0().a1()==INT(0));
    }{
        RAT1_EXT a0(4,3,7);
        RAT1_EXT a1(8,9,7);
        INT      root(5);
        RAT2_EXT a(a0,a1,root);
        typename RAT2_EXT_ANT::Normalization_factor normalization_factor;
        RAT2_EXT r=normalization_factor(a);
        assert((r*a).a1()==   RAT1_EXT(0));
        assert((r*a).a0().a1()==RAT(0));
    }{
        INT1_EXT a(INT(3));
        typename INT1_EXT_ANT::Normalization_factor normalization_factor;
        INT1_EXT r=normalization_factor(a);
        assert(r.a0()    ==INT(1));
        assert(r.a1()    ==INT(0));
        assert(r.root()  ==INT(0));
        assert(!(r*a).is_extended());
    }

    // denomiantor for algebraic integers
    {
        typedef typename AT::Integer  Integer;
        typedef CGAL::Sqrt_extension<Integer,Integer,ACDE_TAG> Extn_1;
        typedef CGAL::Algebraic_extension_traits<Extn_1> ANT;
        typename ANT::Denominator_for_algebraic_integers dfai;
        Extn_1 ext(1,2,5);
        assert(dfai(ext)==Extn_1(20));
    }{
        typedef typename AT::Integer  Integer;
        typedef CGAL::Sqrt_extension<Integer,Integer,ACDE_TAG> Extn_1;
        typedef CGAL::Sqrt_extension<Extn_1, Integer,ACDE_TAG> Extn_2;
        typedef CGAL::Algebraic_extension_traits<Extn_2> ANT;
        typename ANT::Denominator_for_algebraic_integers dfai;
        {
            Extn_1 a0(1);
            Extn_1 a1(2);
            Integer root(5);
            Extn_2 ext(a0,a1,root);
            assert(dfai(ext)==Extn_2(20));

        }{
            Extn_1 a0(1,2,5);
            Extn_1 a1(2,3,5);
            Integer root(7);
            Extn_2 ext(a0,a1,root);
            assert(dfai(ext)==Extn_2(20 *28));
        }
    }{
        typedef typename AT::Integer  Integer;
        typedef CGAL::Sqrt_extension<Integer,Integer,ACDE_TAG>  Extn_1;
        typedef CGAL::Sqrt_extension<Extn_1, Extn_1,ACDE_TAG>   Extn_2;
        typedef CGAL::Algebraic_extension_traits<Extn_2> ANT;
        typename ANT::Denominator_for_algebraic_integers dfai;
        {
            Extn_1 a0(1);
            Extn_1 a1(2);
            Extn_1 root(5);
            Extn_2 ext(a0,a1,root);
            assert(dfai(ext)==Extn_2(20));

        }{
            Extn_1 a0(1,2,5);
            Extn_1 a1(2,3,5);
            Extn_1 root(7);
            Extn_2 ext(a0,a1,root);
            assert(dfai(ext)==Extn_2(20 *28));
        }{
            Extn_1 a0(1,2,5);
            Extn_1 a1(2,3,5);
            Extn_1 root(4,5,5);
            Extn_2 ext(a0,a1,root);
            assert(dfai(ext)==Extn_2(20)*Extn_2(4)*Extn_2(Extn_1(4,5,5)));
        }
        {

            std::vector<Extn_2> vec;
            Extn_1 root(Extn_1(2,3,5));
            vec.push_back(Extn_2(0));
            assert(dfai(vec.begin(),vec.end())== Extn_2(1));
            vec.push_back(Extn_2(root));
            assert(dfai(vec.begin(),vec.end())== Extn_2(20));
            vec.push_back(Extn_2(Extn_1(1),Extn_1(3),root));
            assert(dfai(vec.begin(),vec.end())== Extn_2(20)*Extn_2(root)*Extn_2(4));
        }

    }
}

template<class AT, class ACDE_TAG>
void test_get_arithmetic_kernel(){
  typedef typename AT::Integer Integer;
  typedef typename AT::Rational Rational;
  {
    typedef CGAL::Sqrt_extension<Integer,Integer,ACDE_TAG> EXT;
    typedef typename CGAL::Get_arithmetic_kernel<EXT>::Arithmetic_kernel AT_;
    CGAL_USE_TYPE(AT_);
    static_assert(std::is_same<AT,AT_>::value);
  } {
    typedef CGAL::Sqrt_extension<Rational,Integer,ACDE_TAG> EXT;
    typedef typename CGAL::Get_arithmetic_kernel<EXT>::Arithmetic_kernel AT_;
    CGAL_USE_TYPE(AT_);
    static_assert(std::is_same<AT,AT_>::value);
  } {
    typedef CGAL::Sqrt_extension<Rational,Rational,ACDE_TAG> EXT;
    typedef typename CGAL::Get_arithmetic_kernel<EXT>::Arithmetic_kernel AT_;
    CGAL_USE_TYPE(AT_);
    static_assert(std::is_same<AT,AT_>::value);
  }
}

template <class AT, class ACDE_TAG>
void sqrt_extension_test(){
    CGAL_SNAP_ARITHMETIC_KERNEL_TYPEDEFS(AT);

    general_test<Integer,Integer,CGAL::Integral_domain_tag,ACDE_TAG>();
    general_test<Rational,Integer,CGAL::Field_tag,ACDE_TAG>();

    convert_to_real<Integer,Integer,Field_with_sqrt,ACDE_TAG>();

    convert_to_bfi<Integer,Integer,typename AT::Bigfloat_interval,ACDE_TAG>();


    sqrt_ext_pretty_output_test<AT,ACDE_TAG>();
    fraction_traits_test<AT,ACDE_TAG>();

    scalar_factor_traits_test<AT,ACDE_TAG>();
    test_algebraic_extension_traits<AT,ACDE_TAG>();

    test_get_arithmetic_kernel<AT,ACDE_TAG>();
}

#include <CGAL/Number_types/internal/Exact_type_selector.h>
void test_nt_converter()
{
  typedef CGAL::internal::Exact_field_selector<int>::Type NT;
  typedef CGAL::Sqrt_extension<double,double> Source;
  typedef CGAL::Sqrt_extension<NT,NT> Target;

  CGAL::NT_converter<Source,Target> converter;

  Source s;
  Target t=converter(s);
}

