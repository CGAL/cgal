// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)

// Test program for the Root_of_2 class (not other models of the concept yet).

#include <iostream>
#include <cassert>

#include <CGAL/Random.h>

#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Root_of_traits.h>
#include <iomanip>

#include <CGAL/Test/_test_real_embeddable.h>

#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#  include <CGAL/Gmpq.h>
#endif

#ifdef CGAL_USE_GMPXX
#  include <CGAL/gmpxx.h>
#endif

#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
#endif

#ifdef CGAL_USE_CORE
#  include <CGAL/CORE_Expr.h>
#endif

#include <CGAL/disable_warnings.h>

// We should put these in a nested namespace and use it.
using CGAL::compare; // for double
using CGAL::sign;    // for double

CGAL::Random rnd;

int test_loops;

// Generate random integral values of 7 bits.
template < typename NT >
NT my_rand()
{
  return rnd.get_int(-64, 63);
}


//----------------------

template<class T,class ACDE_TAG,class FP_TAG>
bool is_RO2_class(const CGAL::Sqrt_extension<T,T,ACDE_TAG,FP_TAG>& ){ return true;}
template<class T>
bool is_RO2_class(const T& ){ return false;}

template<class T,class ACDE_TAG,class FP_TAG>
CGAL::Sqrt_extension<T,T,ACDE_TAG,FP_TAG> conjugate(const CGAL::Sqrt_extension<T,T,ACDE_TAG,FP_TAG>& R){
  return R.conjugate();
}

template<class T>
T conjugate(const T&){ return T(); }


template<class T>
T inverse_helper(const T& R){ return (T) 1/R; }

template<class T>
bool is_smaller_helper(const T& R){ return R.exact().a1() <= 0; }

template<class T,class Tag1,class Tag2>
bool is_smaller_helper(const CGAL::Sqrt_extension<T,T,Tag1,Tag2>& R){ return R.a1()<= 0;}

template < typename RT >
struct bracket {
  template < typename T >
  RT operator()(const T& R, int i) const { return bracket<typename RT::ET>().operator()(R.exact(),i); }
  template <class FT,class Tag1,class Tag2>
  RT operator()(const CGAL::Sqrt_extension<FT,FT,Tag1,Tag2>& R, int i) const {
    typedef CGAL::Rational_traits< FT > Rational;
    assert((i>=0) & (i<3));
    Rational r;
    const RT r1 = r.numerator(R.a0());
    const RT d1 = r.denominator(R.a0());
    const RT r2 = r.numerator(R.a1());
    const RT d2 = r.denominator(R.a1());
    const RT r3 = r.numerator(R.root());
    const RT d3 = r.denominator(R.root());
    if(i == 0) {
      return (CGAL_NTS square(d2)) * d3;
    }
    if(i == 1) {
      return -2 * (CGAL_NTS square(d2)) * d3 * r1;
    }
    // i == 2
    return ((CGAL_NTS square(d2)) * d3 * (CGAL_NTS square(r1))) -
           ((CGAL_NTS square(d1)) * r3 * (CGAL_NTS square(r2)));
  }
};


//----------------------------
template < class Root, class RT >
Root create_root_helper(RT a, RT b, Root *)
{
   return Root(a)/b;
}

template <class RT, class Tag1, class Tag2 >
CGAL::Sqrt_extension<double,double,Tag1,Tag2> create_root_helper(RT a, RT b,CGAL::Sqrt_extension<double,double,Tag1,Tag2> *)
{
   return CGAL::Sqrt_extension<double,double,Tag1,Tag2>(a)/b;
}

template < class T, class RT,class Tag1, class Tag2 >
CGAL::Sqrt_extension<T,T,Tag1,Tag2> create_root_helper(RT a, RT b, CGAL::Sqrt_extension<T,T,Tag1,Tag2> *)
{
  typename CGAL::Fraction_traits<T>::Compose comp;
  return CGAL::Sqrt_extension<T,T,Tag1,Tag2>( comp(a, b) );
}

template < class Root, class RT >
Root create_root(RT a, RT b)
{
   return create_root_helper(a, b, (Root*) nullptr);
}

// Generate random Root_of_2 of degree 1
template < typename Root,typename RT >
Root my_rand_root_1()
{
  Root r;
  do {
    RT a = my_rand<RT>();
    if (a == 0)
      continue;
    r = create_root<Root>(my_rand<RT>(), a);
  } while (! is_valid(r));
  return r;
}

// Generate random Root_of_2 of degree 2
template < typename Root,typename RT >
Root my_rand_root()
{
  do {
    RT a = my_rand<RT>();
    RT b = my_rand<RT>();
    RT c = my_rand<RT>();
    if (a == 0 && b == 0)
      return Root(c);
    if (a == 0)
      return create_root<Root>(c, b);
    if (b*b-4*a*c >= 0)
      return CGAL::make_root_of_2(a, b, c, rnd.get_bool());
  } while (true);
}

// Compares the comparison results of two roots with
// the comparison results of their double approximations.
template < typename Root >
bool
test_compare(const Root &r1, const Root &r2)
{
  if (CGAL_NTS compare(r1, r2) == compare(to_double(r1), to_double(r2)))
    return true;

  std::cout << "ERROR in comparing :" << std::endl << " r1 = ";
  print(std::cout, r1);
  std::cout << " [approx = " << to_double(r1) << " ]\n" << " r2 = ";
  print(std::cout, r2);
  std::cout << " [approx = " << to_double(r2) << " ]\n"
            << " comparison gives : " << compare(r1, r2) << " for exact, but "
            << compare(to_double(r1), to_double(r2))
            << " for double approximation" << std::endl;
  return false;
}

template < typename RT, typename Root >
bool
test_to_interval(const Root &r1)
{
  std::pair<double, double> the_interval = CGAL_NTS to_interval(r1);
  if(!((to_double(r1) >= the_interval.first) && (to_double(r1) <= the_interval.second))) {
    std::cout << bracket<RT>()(r1,2) << " " << bracket<RT>()(r1,1) << " " << bracket<RT>()(r1,0) << " " << is_smaller_helper(r1) << std::endl;
    std::cout << std::setprecision (18) << to_double(r1) << std::endl;
    std::cout << "[" << std::setprecision (18) << the_interval.first << "," << std::setprecision (18) << the_interval.second << "]";
  }
  return (to_double(r1) >= the_interval.first) && (to_double(r1) <= the_interval.second);
}




template < typename Root,typename RT,typename FT >
bool
test_root_of_g()
{
  CGAL::test_real_embeddable<Root>();

  std::cout << "  Testing zeros" << std::endl;
  Root zero1=CGAL::make_root_of_2((RT) 1,(RT)0,(RT)0, true);
  Root zero2=CGAL::make_root_of_2((RT)1,(RT)0,(RT)0, false);
  Root zero3(0);

  if (is_RO2_class(zero1)) assert(conjugate(zero1) == zero2);
  if (is_RO2_class(zero2)) assert(conjugate(zero2) == zero1);
  if (is_RO2_class(zero3)) assert(conjugate(zero3) == zero3);
  assert(is_valid(zero1));
  assert(is_valid(zero2));
  assert(is_valid(zero3));
  assert(CGAL_NTS to_double(zero1) == 0.0);
  assert(CGAL_NTS to_double(zero2) == 0.0);
  assert(CGAL_NTS to_double(zero3) == 0.0);
  assert(test_to_interval<RT>(zero1));
  assert(test_to_interval<RT>(zero2));
  assert(test_to_interval<RT>(zero3));
  assert(CGAL_NTS sign(zero1) == 0);
  assert(CGAL_NTS sign(zero2) == 0);
  assert(CGAL_NTS sign(zero3) == 0);
  assert(CGAL_NTS compare(zero1, zero1) == 0);
  assert(CGAL_NTS compare(zero2, zero2) == 0);
  assert(CGAL_NTS compare(zero1, zero2) == 0);
  assert(CGAL_NTS compare(zero2, zero1) == 0);
  assert(zero2 == zero1);
  assert(zero2 == zero3);
  assert(zero1 == zero3);
  assert(zero2 <= zero1);
  assert(zero2 <= zero3);
  assert(zero1 <= zero3);
  assert(zero2 >= zero1);
  assert(zero2 >= zero3);
  assert(zero1 >= zero3);
  assert(!(zero2 != zero1));
  assert(!(zero2 != zero3));
  assert(!(zero1 != zero3));
  assert(!(zero2 < zero1));
  assert(!(zero2 < zero3));
  assert(!(zero1 < zero3));
  assert(!(zero2 > zero1));
  assert(!(zero2 > zero3));
  assert(!(zero1 > zero3));

  assert(CGAL_NTS compare(zero2, -zero1) == 0);
  assert(CGAL_NTS compare(zero1, -zero2) == 0);

  assert(CGAL_NTS compare(zero1, zero3) == 0);
  assert(CGAL_NTS compare(zero2, zero3) == 0);

  std::cout << "  Testing ones" << std::endl;
  Root one1=CGAL::make_root_of_2 ((RT)-1,(RT) 0,(RT) 1, false);
  Root mone1=CGAL::make_root_of_2((RT)-1,(RT) 0,(RT) 1, true);
  Root one2  = Root(1);
  Root mone2 = Root(-1);
  Root one3  = create_root<Root>((RT) 1,(RT) 1);
  Root mone3 = create_root<Root>(-1, 1);

  if (is_RO2_class(one1)) assert(conjugate(one1) == mone1);

  //It is not true that those must hold
  //assert(one2.conjugate() == mone2);
  //assert(one3.conjugate() == mone3);
  assert(is_valid(one1));
  assert(is_valid(mone1));
  assert(is_valid(one2));
  assert(is_valid(mone2));
  assert(CGAL_NTS to_double(one1)  ==  1.0);
  assert(CGAL_NTS to_double(mone1) == -1.0);
  assert(CGAL_NTS to_double(one2)  ==  1.0);
  assert(CGAL_NTS to_double(mone2) == -1.0);
  assert(test_to_interval<RT>(one1));
  assert(test_to_interval<RT>(mone1));
  assert(test_to_interval<RT>(one2));
  assert(test_to_interval<RT>(mone2));
  assert(test_to_interval<RT>(one3));
  assert(test_to_interval<RT>(mone3));

  assert(CGAL_NTS sign( one1) > 0);
  assert(CGAL_NTS sign( one2) > 0);
  assert(CGAL_NTS sign(mone1) < 0);
  assert(CGAL_NTS sign(mone2) < 0);
  assert(CGAL_NTS compare( one1,  one1) == 0);
  assert(CGAL_NTS compare(mone1, mone1) == 0);
  assert(CGAL_NTS compare(mone1,  one1) < 0);
  assert(CGAL_NTS compare( one1, mone1) > 0);

  assert(CGAL_NTS compare( one2,  one2) == 0);
  assert(CGAL_NTS compare(mone2, mone2) == 0);
  assert(CGAL_NTS compare(mone2,  one2) < 0);
  assert(CGAL_NTS compare( one2, mone2) > 0);

  assert(CGAL_NTS compare( one1, -mone1) == 0);
  assert(CGAL_NTS compare(-one1,  mone1) == 0);

  assert(CGAL_NTS compare( one1,  one2) == 0);
  assert(CGAL_NTS compare(mone1, mone2) == 0);
  assert(CGAL_NTS compare( one1,  one3) == 0);
  assert(CGAL_NTS compare(mone1, mone3) == 0);
  assert(CGAL_NTS compare( one2,  one3) == 0);
  assert(CGAL_NTS compare(mone2, mone3) == 0);

  assert(CGAL_NTS compare( one1, zero1) > 0);
  assert(CGAL_NTS compare(mone1, zero2) < 0);

  std::cout << "  Testing output" << std::endl;
  std::cout << "  zero = ";
  print(std::cout, zero1);
  std::cout << " approx =  " << zero1 << std::endl;

  std::cout << "  zero = ";
  print(std::cout, Root(0));
  std::cout << " approx =  " << Root(0) << std::endl;

  std::cout << "  one  = ";
  print(std::cout, one1);
  std::cout << " approx =  " << one1 << std::endl;

  std::cout << "  mone = ";
  print(std::cout, mone1);
  std::cout << " approx =  " << mone1 << std::endl;

  std::cout << "  Testing degree 1" << std::endl;
  assert(CGAL_NTS compare( Root(0),  Root(0)) == 0);
  assert(CGAL_NTS compare( Root(1),  Root(1)) == 0);
  assert(CGAL_NTS compare(Root(-1), Root(-1)) == 0);
  assert(CGAL_NTS compare( Root(0),  Root(1)) < 0);
  assert(CGAL_NTS compare( Root(1), Root(-1)) > 0);

  assert(CGAL_NTS compare( create_root<Root>(-4, 2), Root(-2)) == 0);
  assert(CGAL_NTS compare( create_root<Root>(-4, 2), Root(-1)) < 0);
  assert(CGAL_NTS compare( create_root<Root>(-4, 2), Root(-3)) > 0);

  std::cout << "  Testing degree 1 and 2" << std::endl;
  Root rone=CGAL::make_root_of_2((RT)1, (RT)0, (RT)-1, false);
  Root rmone=CGAL::make_root_of_2((RT)1, (RT)0,(RT) -1, true);
  if (is_RO2_class(rone)) assert(conjugate(rone) == rmone);
  assert(test_to_interval<RT>(rone));
  assert(test_to_interval<RT>(rmone));
  // Compare the two roots of the above polynomial (-1 and 1),
  // succesively with -2, -1, 0, 1, 2.
  assert(CGAL_NTS compare(Root(-2), rmone) < 0);
  assert(CGAL_NTS compare(Root(-2), rone)  < 0);
  assert(CGAL_NTS compare(rmone, Root(-2)) > 0);
  assert(CGAL_NTS compare(rone, Root(-2))  > 0);
  assert(CGAL_NTS compare(Root(-1), rmone) == 0);
  assert(CGAL_NTS compare(Root(-1), rone)  < 0);
  assert(CGAL_NTS compare(rmone, Root(-1)) == 0);
  assert(CGAL_NTS compare(rone, Root(-1))  > 0);
  assert(CGAL_NTS compare(Root(0), rmone)  > 0);
  assert(CGAL_NTS compare(Root(0), rone)   < 0);
  assert(CGAL_NTS compare(rmone, Root(0))  < 0);
  assert(CGAL_NTS compare(rone, Root(0))   > 0);
  assert(CGAL_NTS compare(Root(1), rmone)  > 0);
  assert(CGAL_NTS compare(Root(1), rone)   == 0);
  assert(CGAL_NTS compare(rmone, Root(1))  < 0);
  assert(CGAL_NTS compare(rone, Root(1))   == 0);
  assert(CGAL_NTS compare(Root(2), rmone)  > 0);
  assert(CGAL_NTS compare(Root(2), rone)   > 0);
  assert(CGAL_NTS compare(rmone, Root(2))  < 0);
  assert(CGAL_NTS compare(rone, Root(2))   < 0);

// All compare tests are wrong since the double
// cannot handle all of them

  std::cout << "  Testing random roots constructed by constants" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = Root(my_rand<RT>());
    Root r2 = Root(my_rand<RT>());
    assert(r1 == r1);
    assert(test_to_interval<RT>(r1));
    assert(r2 == r2);
    assert(test_to_interval<RT>(r2));
//    assert(test_compare(r1, r2));
  }

  std::cout << "  Testing random roots of degree 1" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = my_rand_root_1<Root,RT>();
    Root r2 = my_rand_root_1<Root,RT>();
    assert(r1 == r1);
    assert(test_to_interval<RT>(r1));
    assert(r2 == r2);
    assert(test_to_interval<RT>(r2));
//    assert(test_compare(r1, r2));
  }

  std::cout << "  Testing random roots of degree 2" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = my_rand_root<Root,RT>();
    Root r2 = my_rand_root<Root,RT>();
    assert(r1 == r1);
    assert(test_to_interval<RT>(r1));
    assert(r2 == r2);
    assert(test_to_interval<RT>(r2));
//    assert(test_compare(r1, r2));
  }

  std::cout << "  Testing random roots of degree 1 and 2" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = my_rand_root_1<Root,RT>();
    Root r2 = my_rand_root<Root,RT>();
    assert(r1 == r1);
    assert(test_to_interval<RT>(r1));
    assert(r2 == r2);
    assert(test_to_interval<RT>(r2));
//    assert(test_compare(r1, r2));
  }

  std::cout << "  Testing squares of random roots of degree 2" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = my_rand_root<Root,RT>();
    Root r2 = my_rand_root<Root,RT>();
    Root r1_sqr = CGAL_NTS square(r1);
    Root r2_sqr = CGAL_NTS square(r2);
    assert(test_to_interval<RT>(r1));
    assert(test_to_interval<RT>(r2));
    assert(test_to_interval<RT>(r1_sqr));
    assert(test_to_interval<RT>(r2_sqr));
//   double dr1 = to_double(r1);
//   double dr2 = to_double(r2);
//    assert(CGAL_NTS compare(r1_sqr, r2_sqr) == compare(dr1*dr1, dr2*dr2));
  }

  std::cout << "  Testing addition/subtraction of Sqrt_extension<NT,NT> with NT"
            << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    RT r    = my_rand<RT>();
    Root r1 = my_rand_root<Root,RT>();
    Root r2 = r1 - r;
    Root r3 = r1 + r;
    Root r4 = r - r1;
    Root r5 = r + r1;
    assert(r1 == r1);
    assert(r2 == r2);
    assert(r3 == r3);
    assert(r4 == r4);
    assert(r5 == r5);
    assert(test_to_interval<RT>(r1));
    assert(test_to_interval<RT>(r2));
    assert(test_to_interval<RT>(r3));
    assert(test_to_interval<RT>(r4));
    assert(test_to_interval<RT>(r5));
  //  assert(test_compare(r1, r2));
  //  assert(test_compare(r1, r3));
  //  assert(test_compare(r1, r4));
  //  assert(test_compare(r1, r5));
  //  assert(CGAL_NTS compare(r1, r2) ==   (int) CGAL_NTS sign(r));
  //  assert(CGAL_NTS compare(r1, r3) == - (int) CGAL_NTS sign(r));
  //  assert(CGAL_NTS compare(r1, r5) == - (int) CGAL_NTS sign(r));
  //  assert(CGAL_NTS compare(r, r4)  == (int) CGAL_NTS sign(r1));
  }

  std::cout << "  Testing multiplication of Sqrt_extension<NT,NT> with NT" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    RT r    = my_rand<RT>();
    Root r1 = my_rand_root<Root,RT>();
    Root r2 = r1 * r;
    Root r3 = r * r1;
    int n = rnd.get_int(0, 63);
    FT rn(n);
    if (r != 0){
      Root r4 = r2 / r;
      assert(r4 == r4);
      assert(test_to_interval<RT>(r4));
      assert(r4 == r1);
    }
    if (n != 0){
      Root r5 = my_rand_root<Root,RT>();
      Root r6 = r5 * rn;
      Root r4 = r6 / rn;
      assert(r4 == r4);
      assert(test_to_interval<RT>(r4));
      assert(r4 == r5);
    }
    assert(r1 == r1);
    assert(r2 == r2);
    assert(r3 == r3);
    assert(test_to_interval<RT>(r1));
    assert(test_to_interval<RT>(r2));
    assert(test_to_interval<RT>(r3));
//    assert(test_compare(r1, r2));
//    assert(test_compare(r1, r3));
//    assert(test_compare(r2, r3));
    assert(r2 == r3);
    if (r > 0) {
//      assert(CGAL_NTS compare(r1, r2) == compare(1, r) * sign(r1));
//      assert(CGAL_NTS compare(r1, r3) == compare(1, r) * sign(r1));
    } else if (r < 0) {
//      assert(CGAL_NTS compare(r1, -r2) == compare(1, -r) * sign(r1));
//      assert(CGAL_NTS compare(r1, -r3) == compare(1, -r) * sign(r1));
    } else {
      assert(r2 == Root(0));
      assert(r3 == Root(0));
    }
  }

  std::cout << "  Testing the inverse of random roots of degree 2" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = my_rand_root<Root,RT>();
    Root r2 = my_rand_root<Root,RT>();
    while(r1 == 0) r1 = my_rand_root<Root,RT>();
    while(r2 == 0) r2 = my_rand_root<Root,RT>();
    Root r1_inv = inverse_helper(r1);
    Root r2_inv = inverse_helper(r2);
    assert(test_to_interval<RT>(r1));
    assert(test_to_interval<RT>(r2));
    assert(test_to_interval<RT>(r1_inv));
    assert(test_to_interval<RT>(r2_inv));
//    double dr1 = to_double(r1);
//    double dr2 = to_double(r2);
//    assert(CGAL_NTS compare(r1_inv, r2_inv) == compare(1.0/dr1, 1.0/dr2));
  }

  std::cout << "  Testing make_sqrt(RT)" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    RT r1 = my_rand<RT>();
    RT r2 = my_rand<RT>();
    while(r1 < 0) r1 = my_rand<RT>();
    while(r2 < 0) r2 = my_rand<RT>();
    Root sqr_r1 = CGAL::make_sqrt(r1);
    Root sqr_r2 = CGAL::make_sqrt(r2);
    assert(test_to_interval<RT>(sqr_r1));
    assert(test_to_interval<RT>(sqr_r2));
//    double dr1 = to_double(r1);
//    double dr2 = to_double(r2);
//    assert(CGAL_NTS compare(sqr_r1, sqr_r2) == compare(std::sqrt(dr1), std::sqrt(dr2)));
  }

  std::cout << "  Testing Sqrt_extension<FT,FT>" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    int n = rnd.get_int(0, 63);
    typedef CGAL::Rational_traits<FT> Rat_traits;
    FT r(n);
    Root r1(r);
    Root r2(n);
    Root r3=create_root<Root>(Rat_traits().numerator(r), Rat_traits().denominator(r));
    assert(r1 == r1);
    assert(r2 == r2);
    assert(r3 == r3);
    assert(r1 == r2);
    assert(r2 == r3);
    assert(r1 == r3);
  }
  return true;
}

template<typename Root >
bool
test_root_of(){
  typedef typename Root::NT FT;
  typedef typename CGAL::Fraction_traits<FT>::Numerator_type RT;
  return test_root_of_g<Root,RT,FT>();
}

// ALL THOSE TESTS ARE BASED ON COMPARING THE DOUBLE RESULT
// WITH THE EXACT RESULT
// THAT IS NOT NECESSARILY CORRECT

int main(int argc, char **argv) {

  test_loops = argc == 2 ? atoi(argv[1]) : 1000;

  bool result = true;

  // Sqrt_extension requires a FT as template parameter
  //std::cout << "Testing Sqrt_extension<double,double>" << std::endl;
  //result = result && test_root_of_g<CGAL::Sqrt_extension<double,double,CGAL::Tag_true,CGAL::Tag_true>,double,double>();

  std::cout << "Testing Sqrt_extension with Quotient<MP_Float>" << std::endl;
  result = result && test_root_of<CGAL::Sqrt_extension<CGAL::Quotient<CGAL::MP_Float>,CGAL::Quotient<CGAL::MP_Float>,CGAL::Tag_true,CGAL::Tag_true> >();

  std::cout << "Testing Lazy_exact_nt<MP_Float>'s RootOf_2 " << std::endl;
  result = result &&
      test_root_of_g<CGAL::Root_of_traits<CGAL::Lazy_exact_nt<CGAL::MP_Float> >
        ::RootOf_2,CGAL::Lazy_exact_nt<CGAL::MP_Float>,CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > >();

#ifdef CGAL_USE_GMP
  std::cout << "Testing Sqrt_extension with Gmpq" << std::endl;
  result = result && test_root_of<CGAL::Sqrt_extension<CGAL::Gmpq,CGAL::Gmpq,CGAL::Tag_true,CGAL::Tag_true> >();

  std::cout << "Testing Lazy_exact_nt<Gmpz>'s RootOf_2 " << std::endl;
  result = result &&
      test_root_of_g<CGAL::Root_of_traits<CGAL::Lazy_exact_nt<CGAL::Gmpz> >
        ::RootOf_2,CGAL::Lazy_exact_nt<CGAL::Gmpz>,CGAL::Lazy_exact_nt<CGAL::Gmpq> >();
#endif

#ifdef CGAL_USE_GMPXX
  std::cout << "Testing Sqrt_extension with mpq_class" << std::endl;
  result = result && test_root_of<CGAL::Sqrt_extension<mpq_class,mpq_class,CGAL::Tag_true,CGAL::Tag_true > >();
#endif

  if (result) {
    std::cout << "OK" << std::endl;
    return 0;
  } else {
    std::cout << "ERROR" << std::endl;
    return -1;
  }
}
