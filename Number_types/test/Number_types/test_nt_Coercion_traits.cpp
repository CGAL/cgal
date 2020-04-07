// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================


#include <CGAL/use.h>
#include <cassert>
#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Coercion_traits.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Interval_nt.h>

#include <CGAL/Test/_test_coercion_traits.h>

typedef CGAL::Interval_nt<false> Interval;

void test_coercion_from_to(CGAL::Null_tag, CGAL::Null_tag){}
template<class A> void test_coercion_from_to(A, CGAL::Null_tag){}
template<class B> void test_coercion_from_to(CGAL::Null_tag, B){}
template<class A, class B> void test_coercion_from_to(A, B){
    CGAL::test_explicit_interoperable_from_to<A,B>();
}

template <class AT>
void check_type_coercion_at(){
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;
    typedef typename AT::Field_with_sqrt Field_with_sqrt;
    typedef typename AT::Field_with_kth_root Field_with_kth_root;
    typedef typename AT::Field_with_root_of  Field_with_root_of;
    typedef typename AT::Bigfloat Bigfloat;
    typedef typename AT::Bigfloat_interval Bigfloat_interval;


    test_coercion_from_to(int(),Integer());
    test_coercion_from_to(short(),Integer());
    test_coercion_from_to(Integer(),Integer());

    test_coercion_from_to(int(),Rational());
    test_coercion_from_to(short(),Rational());
    test_coercion_from_to(float(),Rational());
    test_coercion_from_to(double(),Rational());
    test_coercion_from_to(Integer(),Rational());
    test_coercion_from_to(Bigfloat(),Rational());
    test_coercion_from_to(Rational(),Rational());

    test_coercion_from_to(int(),Field_with_sqrt());
    test_coercion_from_to(short(),Field_with_sqrt());
    test_coercion_from_to(float(),Field_with_sqrt());
    test_coercion_from_to(double(),Field_with_sqrt());
    test_coercion_from_to(Integer(),Field_with_sqrt());
    test_coercion_from_to(Bigfloat(),Field_with_sqrt());
    test_coercion_from_to(Rational(),Field_with_sqrt());
    test_coercion_from_to(Field_with_sqrt(),Field_with_sqrt());

    test_coercion_from_to(int(),Field_with_kth_root());
    test_coercion_from_to(short(),Field_with_kth_root());
    test_coercion_from_to(float(),Field_with_kth_root());
    test_coercion_from_to(double(),Field_with_kth_root());
    test_coercion_from_to(Integer(),Field_with_kth_root());
    test_coercion_from_to(Bigfloat(),Field_with_kth_root());
    test_coercion_from_to(Rational(),Field_with_kth_root());
    test_coercion_from_to(Field_with_sqrt(),Field_with_kth_root());
    test_coercion_from_to(Field_with_kth_root(),Field_with_kth_root());


    test_coercion_from_to(int(),Field_with_root_of());
    test_coercion_from_to(short(),Field_with_root_of());
    test_coercion_from_to(float(),Field_with_root_of());
    test_coercion_from_to(double(),Field_with_root_of());
    test_coercion_from_to(Integer(),Field_with_root_of());
    test_coercion_from_to(Bigfloat(),Field_with_root_of());
    test_coercion_from_to(Rational(),Field_with_root_of());
    test_coercion_from_to(Field_with_sqrt(),Field_with_root_of());
    test_coercion_from_to(Field_with_kth_root(),Field_with_root_of());
    test_coercion_from_to(Field_with_root_of(),Field_with_root_of());


    test_coercion_from_to(int(),Bigfloat());
    test_coercion_from_to(short(),Bigfloat());
    test_coercion_from_to(float(),Bigfloat());
    test_coercion_from_to(double(),Bigfloat());
    test_coercion_from_to(Integer(),Bigfloat());
    test_coercion_from_to(Bigfloat(),Bigfloat());


    test_coercion_from_to(int(),Bigfloat_interval());
    test_coercion_from_to(short(),Bigfloat_interval());
    test_coercion_from_to(float(),Bigfloat_interval());
    test_coercion_from_to(double(),Bigfloat_interval());
    test_coercion_from_to(Integer(),Bigfloat_interval());
    test_coercion_from_to(Bigfloat(),Bigfloat_interval());
    test_coercion_from_to(Rational(),Bigfloat_interval());
    test_coercion_from_to(Field_with_sqrt(),Bigfloat_interval());
    test_coercion_from_to(Field_with_kth_root(),Bigfloat_interval());
    test_coercion_from_to(Field_with_root_of(),Bigfloat_interval());
    test_coercion_from_to(Bigfloat_interval(),Bigfloat_interval());

}


template <class AT>
void AT_coercion_test_for_cgal_types_int(){
    typedef typename AT::Integer Integer;
    CGAL::test_explicit_interoperable_from_to<int    ,Integer>();
    CGAL::test_explicit_interoperable_from_to<Integer,Integer>();
    CGAL::test_explicit_interoperable_from_to<Integer,Interval>();

    // Quotient
    typedef CGAL::Quotient<Integer> Quotient;

    CGAL::test_explicit_interoperable_from_to<int     ,Quotient>();
    CGAL::test_explicit_interoperable_from_to<double  ,Quotient>();
    CGAL::test_explicit_interoperable_from_to<Integer ,Quotient>();
    CGAL::test_explicit_interoperable_from_to<Quotient,Quotient>();
    CGAL::test_explicit_interoperable_from_to<Quotient,Interval>();

    typedef CGAL::Sqrt_extension<Integer, Integer> Extn_1;
    CGAL::test_explicit_interoperable_from_to<int    ,Extn_1>();
    CGAL::test_explicit_interoperable_from_to<Integer,Extn_1>();
    CGAL::test_explicit_interoperable_from_to<Extn_1 ,Extn_1>();
    CGAL::test_explicit_interoperable_from_to<Extn_1,Interval>();

    typedef CGAL::Sqrt_extension<Extn_1,  Extn_1 > Extn_2n;
    CGAL::test_explicit_interoperable_from_to<int    ,Extn_2n>();
    CGAL::test_explicit_interoperable_from_to<Integer,Extn_2n>();
    CGAL::test_explicit_interoperable_from_to<Extn_1 ,Extn_2n>();
    CGAL::test_explicit_interoperable_from_to<Extn_2n,Extn_2n>();
    CGAL::test_explicit_interoperable_from_to<Extn_2n,Interval>();

    typedef CGAL::Sqrt_extension<Extn_1,  Integer> Extn_2d;
    CGAL::test_explicit_interoperable_from_to<int    ,Extn_2d>();
    CGAL::test_explicit_interoperable_from_to<Integer,Extn_2d>();
    CGAL::test_explicit_interoperable_from_to<Extn_1 ,Extn_2d>();
    CGAL::test_explicit_interoperable_from_to<Extn_2d,Extn_2d>();
    CGAL::test_explicit_interoperable_from_to<Extn_2d,Interval>();
}

template <class AT>
void AT_coercion_test_for_cgal_types_rat(){

    AT_coercion_test_for_cgal_types_int<AT>();

    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;
    typedef typename AT::Bigfloat_interval Bigfloat_interval;

    CGAL_USE_TYPE(Bigfloat_interval);
    CGAL_static_assertion(!(::boost::is_same<Integer, CGAL::Null_tag>::value));
    CGAL_static_assertion(!(::boost::is_same<Rational, CGAL::Null_tag>::value));
    CGAL_static_assertion(!(::boost::is_same<Bigfloat_interval, CGAL::Null_tag>::value));

    CGAL::test_explicit_interoperable_from_to<int     ,Rational>();
    CGAL::test_explicit_interoperable_from_to<double  ,Rational>();
    CGAL::test_implicit_interoperable_from_to<Integer ,Rational>();
    CGAL::test_implicit_interoperable_from_to<Rational,Rational>();
    CGAL::test_explicit_interoperable_from_to<Rational,Interval>();


    typedef CGAL::Sqrt_extension<Integer  , Integer> Extn_1;
    typedef CGAL::Sqrt_extension<Extn_1   , Extn_1 > Extn_2n;
    CGAL_USE_TYPE(Extn_2n);
    typedef CGAL::Sqrt_extension<Extn_1   , Integer> Extn_2d;
    CGAL_USE_TYPE(Extn_2d);

    typedef CGAL::Sqrt_extension<Rational  , Integer> Extn_rat_int;
    CGAL::test_explicit_interoperable_from_to<int         ,Extn_rat_int>();
//  CGAL::test_explicit_interoperable_from_to<double      ,Extn_rat_int>();
    CGAL::test_explicit_interoperable_from_to<Integer     ,Extn_rat_int>();
    CGAL::test_explicit_interoperable_from_to<Rational    ,Extn_rat_int>();
    CGAL::test_explicit_interoperable_from_to<Extn_1      ,Extn_rat_int>();
    CGAL::test_explicit_interoperable_from_to<Extn_rat_int,Extn_rat_int>();

    typedef CGAL::Sqrt_extension<Rational  , Rational> Extn_rat_1;
    CGAL::test_explicit_interoperable_from_to<int         ,Extn_rat_1>();
//  CGAL::test_explicit_interoperable_from_to<double      ,Extn_rat_1>();
    CGAL::test_explicit_interoperable_from_to<Integer     ,Extn_rat_1>();
    CGAL::test_explicit_interoperable_from_to<Rational    ,Extn_rat_1>();
//  CGAL::test_explicit_interoperable_from_to<Extn_1      ,Extn_rat_1>();
    CGAL::test_explicit_interoperable_from_to<Extn_rat_1,Extn_rat_1>();

    typedef CGAL::Sqrt_extension<Extn_rat_1, Extn_rat_1 > Extn_rat_2n;
    CGAL::test_explicit_interoperable_from_to<int         ,Extn_rat_2n>();
//  CGAL::test_explicit_interoperable_from_to<double      ,Extn_rat_2n>();
    CGAL::test_explicit_interoperable_from_to<Integer     ,Extn_rat_2n>();
    CGAL::test_explicit_interoperable_from_to<Rational    ,Extn_rat_2n>();
    CGAL::test_explicit_interoperable_from_to<Extn_rat_1  ,Extn_rat_2n>();
    CGAL::test_explicit_interoperable_from_to<Extn_rat_2n ,Extn_rat_2n>();

    typedef CGAL::Sqrt_extension<Extn_rat_1, Rational> Extn_rat_2d;
    CGAL::test_explicit_interoperable_from_to<int         ,Extn_rat_2d>();
//  CGAL::test_explicit_interoperable_from_to<double      ,Extn_rat_2d>();
    CGAL::test_explicit_interoperable_from_to<Integer     ,Extn_rat_2d>();
    CGAL::test_explicit_interoperable_from_to<Rational    ,Extn_rat_2d>();
    CGAL::test_explicit_interoperable_from_to<Extn_rat_1  ,Extn_rat_2d>();
    CGAL::test_explicit_interoperable_from_to<Extn_rat_2d ,Extn_rat_2d>();
}


template <class AT>
void AT_coercion_test_for_cgal_types_fws(){
    AT_coercion_test_for_cgal_types_rat<AT>();

    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;
    typedef typename AT::Bigfloat_interval Bigfloat_interval;
    typedef typename AT::Field_with_sqrt Real;

    CGAL_USE_TYPE(Bigfloat_interval);
    CGAL_static_assertion(!(::boost::is_same<Integer, CGAL::Null_tag>::value));
    CGAL_static_assertion(!(::boost::is_same<Rational, CGAL::Null_tag>::value));
    CGAL_static_assertion(!(::boost::is_same<Bigfloat_interval, CGAL::Null_tag>::value));
    CGAL_static_assertion(!(::boost::is_same<Real, CGAL::Null_tag>::value));


    typedef CGAL::Sqrt_extension<Integer  , Integer> Extn_1;
    typedef CGAL::Sqrt_extension<Extn_1   , Extn_1 > Extn_2n;
    typedef CGAL::Sqrt_extension<Extn_1   , Integer> Extn_2d;

    typedef CGAL::Sqrt_extension<Rational  , Integer>     Extn_rat_int;
    typedef CGAL::Sqrt_extension<Rational  , Rational>    Extn_rat_1;
    typedef CGAL::Sqrt_extension<Extn_rat_1, Extn_rat_1 > Extn_rat_2n;
    typedef CGAL::Sqrt_extension<Extn_rat_1, Rational>    Extn_rat_2d;

    CGAL::test_explicit_interoperable_from_to<int         ,Real>();
    CGAL::test_explicit_interoperable_from_to<double      ,Real>();
    CGAL::test_explicit_interoperable_from_to<Integer     ,Real>();
    CGAL::test_explicit_interoperable_from_to<Rational    ,Real>();
    CGAL::test_explicit_interoperable_from_to<Extn_1      ,Real>();
    CGAL::test_explicit_interoperable_from_to<Extn_2n     ,Real>();
    CGAL::test_explicit_interoperable_from_to<Extn_2d     ,Real>();
    CGAL::test_explicit_interoperable_from_to<Extn_rat_int,Real>();
    CGAL::test_explicit_interoperable_from_to<Extn_rat_1  ,Real>();
    CGAL::test_explicit_interoperable_from_to<Extn_rat_2n ,Real>();
    CGAL::test_explicit_interoperable_from_to<Extn_rat_2d ,Real>();

  // direct casts
  // Integer;
    CGAL::test_explicit_interoperable_from_to<Integer,Interval>();

  // Rational
    CGAL::test_explicit_interoperable_from_to<Rational,Interval>();

  // Interval
    CGAL::test_explicit_interoperable_from_to<Interval,Interval>();
}

// TODO: We have no NT compactified, Matrix, Vector, etc. yet
/*void test_compactified(){
    typedef NiX::Compactified<int>  Compactified_int;
    typedef NiX::Compactified<long> Compactified_long;

    Compactified_int zero(0);
    Compactified_int one(1);
    Compactified_int p_infty(NiX::PLUS_INFTY);
    Compactified_int m_infty(NiX::MINUS_INFTY);

    {
        typedef NiX::Coercion_traits<Compactified_int,Compactified_int> CT;
        CT::Cast cast;
        assert(cast(one) == Compactified_int(1));
        assert(cast(p_infty) == Compactified_int(NiX::PLUS_INFTY));
        assert(cast(m_infty) == Compactified_int(NiX::MINUS_INFTY));
        assert(cast(zero) == Compactified_int(0));
    }{
        typedef NiX::Coercion_traits<Compactified_int,Compactified_long> CT;
        CT::Cast cast;
        assert(cast(Compactified_long(1)) == Compactified_long(1));
        assert(cast(one) == Compactified_long(1));
        assert(cast(p_infty) == Compactified_long(NiX::PLUS_INFTY));
        assert(cast(m_infty) == Compactified_long(NiX::MINUS_INFTY));
        assert(cast(zero) == Compactified_long(0));
    }{
        typedef NiX::Coercion_traits<Compactified_long,Compactified_int> CT;
        CT::Cast cast;
        assert(cast(Compactified_long(1)) == Compactified_long(1));
        assert(cast(one) == Compactified_long(1));
        assert(cast(p_infty) == Compactified_long(NiX::PLUS_INFTY));
        assert(cast(m_infty) == Compactified_long(NiX::MINUS_INFTY));
        assert(cast(zero) == Compactified_long(0));
    }
}

void test_tendency(){
    typedef NiX::Tendency<int>  Tendency_int;
    typedef NiX::Tendency<long> Tendency_long;

    Tendency_int zero(0);

    Tendency_int one_me (1,Tendency_int::MINUS_EPSILON);
    Tendency_int one_mee(1,Tendency_int::MINUS_SQUARED_EPSILON);
    Tendency_int one    (1,Tendency_int::NONE);
    Tendency_int one_pee(1,Tendency_int::PLUS_SQUARED_EPSILON);
    Tendency_int one_pe (1,Tendency_int::PLUS_EPSILON);

    {
        typedef NiX::Coercion_traits<Tendency_int,Tendency_int> CT;
        CT::Cast cast;
        assert(cast(zero) == Tendency_int(0));
        assert(cast(one) == Tendency_int(1));
        assert(cast(one_me)  == Tendency_int(1,Tendency_int::MINUS_EPSILON));
        assert(cast(one_mee) ==Tendency_int(1,Tendency_int::MINUS_SQUARED_EPSILON));
        assert(cast(one)     ==Tendency_int(1,Tendency_int::NONE));
        assert(cast(one_pee) ==Tendency_int(1,Tendency_int::PLUS_SQUARED_EPSILON));
        assert(cast(one_pe)  ==Tendency_int(1,Tendency_int::PLUS_EPSILON));

    }{
        typedef NiX::Coercion_traits<Tendency_int,Tendency_long> CT;
        CT::Cast cast;
        assert(cast(Tendency_long(0))    == Tendency_long(0));
        assert(cast(zero) == Tendency_long(0));
        assert(cast(one) == Tendency_long(1));
        assert(cast(one_me)  == Tendency_long(1,Tendency_long::MINUS_EPSILON));
        assert(cast(one_mee) ==Tendency_long(1,Tendency_long::MINUS_SQUARED_EPSILON));
        assert(cast(one)     ==Tendency_long(1,Tendency_long::NONE));
        assert(cast(one_pee) ==Tendency_long(1,Tendency_long::PLUS_SQUARED_EPSILON));
        assert(cast(one_pe)  ==Tendency_long(1,Tendency_long::PLUS_EPSILON));

    }{
        typedef NiX::Coercion_traits<Tendency_long,Tendency_int> CT;
        CT::Cast cast;
        assert(cast(Tendency_long(0))    == Tendency_long(0));
        assert(cast(zero)    == Tendency_long(0));
        assert(cast(one)     == Tendency_long(1));
        assert(cast(one_me)  ==Tendency_long(1,Tendency_long::MINUS_EPSILON));
        assert(cast(one_mee) ==Tendency_long(1,Tendency_long::MINUS_SQUARED_EPSILON));
        assert(cast(one)     ==Tendency_long(1,Tendency_long::NONE));
        assert(cast(one_pee) ==Tendency_long(1,Tendency_long::PLUS_SQUARED_EPSILON));
        assert(cast(one_pe)  ==Tendency_long(1,Tendency_long::PLUS_EPSILON));
    }
}

void test_matrix_d(){
    typedef NiX::Matrix_d<int > M_int;
    typedef NiX::Matrix_d<long> M_long;

    std::vector<int> input;
    for(int i= 0; i< 4 ;i++) input.push_back(i);
    M_int::Element_range erange_int ;
    M_long::Element_range erange_long ;

    M_int  m_int (erange_int ,2,2,input.begin(),input.end());
    M_long m_long(erange_long,2,2,input.begin(),input.end());

    {
        typedef NiX::Coercion_traits<M_int,M_int> CT;
        CT::Cast cast;
        assert(cast(m_int) == m_int);
    }{
        typedef NiX::Coercion_traits<M_long,M_int> CT;
        CT::Cast cast;
        assert(cast(m_int) == m_long);
        assert(cast(m_long) == m_long);
    }{
        typedef NiX::Coercion_traits<M_int,M_long> CT;
        CT::Cast cast;
        assert(cast(m_int) == m_long);
        assert(cast(m_long) == m_long);
    }
}

void test_vector_d(){
    typedef NiX::Vector_d<int > V_int;
    typedef NiX::Vector_d<long> V_long;


    V_int  v_int (2); v_int [0] = 0;v_int [1] = 3;
    V_long v_long(2); v_long[0] = 0;v_long[1] = 3;


    {
        typedef NiX::Coercion_traits<V_int,V_int> CT;
        CT::Cast cast;
        assert(cast(v_int) == v_int);
    }{
        typedef NiX::Coercion_traits<V_long,V_int> CT;
        CT::Cast cast;
        assert(cast(v_int) == v_long);
        assert(cast(v_long) == v_long);
    }{
        typedef NiX::Coercion_traits<V_int,V_long> CT;
        CT::Cast cast;
        assert(cast(v_int) == v_long);
        assert(cast(v_long) == v_long);
    }
}*/

int main(){
    CGAL::test_implicit_interoperable_from_to<short,int>();
    CGAL::test_implicit_interoperable_from_to<short,long>();
    CGAL::test_implicit_interoperable_from_to<short,long long>();
    CGAL::test_implicit_interoperable_from_to<short,float>();
    CGAL::test_implicit_interoperable_from_to<short,double>();
    CGAL::test_implicit_interoperable_from_to<short,long double>();

    CGAL::test_implicit_interoperable_from_to<int,long>();
    CGAL::test_implicit_interoperable_from_to<int,long long>();
    CGAL::test_implicit_interoperable_from_to<int,float>();
    CGAL::test_implicit_interoperable_from_to<int,double>();
    CGAL::test_implicit_interoperable_from_to<int,long double>();


    CGAL::test_implicit_interoperable_from_to<long,long long>();
    CGAL::test_implicit_interoperable_from_to<long,float>();
    CGAL::test_implicit_interoperable_from_to<long,double>();
    CGAL::test_implicit_interoperable_from_to<long,long double>();


    CGAL::test_implicit_interoperable_from_to<long long,float>();
    CGAL::test_implicit_interoperable_from_to<long long,double>();
    CGAL::test_implicit_interoperable_from_to<long long,long double>();


    CGAL::test_implicit_interoperable_from_to<float,double>();
    CGAL::test_implicit_interoperable_from_to<float,long double>();


    CGAL::test_implicit_interoperable_from_to<double,long double>();

    CGAL::test_implicit_interoperable_from_to<short,short>();
    CGAL::test_implicit_interoperable_from_to<int,int>();
    CGAL::test_implicit_interoperable_from_to<long,long>();
    CGAL::test_implicit_interoperable_from_to<long long,long long>();
    CGAL::test_implicit_interoperable_from_to<float,float>();
    CGAL::test_implicit_interoperable_from_to<double,double>();
    CGAL::test_implicit_interoperable_from_to<long double,long double>();

//    test_compactified();
//    test_tendency();
//    test_matrix_d();
//    test_vector_d();


#ifdef CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
    // since it might be GMP kernel this can only do the rat test.
    AT_coercion_test_for_cgal_types_rat<CGAL::Arithmetic_kernel>();
#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  check_type_coercion_at<CGAL::LEDA_arithmetic_kernel>();
  AT_coercion_test_for_cgal_types_fws<CGAL::LEDA_arithmetic_kernel>();
#endif // CGAL_USE_LEDA

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  check_type_coercion_at<CGAL::CORE_arithmetic_kernel>();
  AT_coercion_test_for_cgal_types_fws<CGAL::CORE_arithmetic_kernel>();
#endif // CGAL_USE_CORE

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  check_type_coercion_at<CGAL::GMP_arithmetic_kernel>();
  AT_coercion_test_for_cgal_types_rat<CGAL::GMP_arithmetic_kernel>();
#endif // CGAL_USE_GMP

  return EXIT_SUCCESS;
}
