// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     :   Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#include <cassert>
#include <algorithm>

//#include <CGAL/_test_basic.h>

#include <CGAL/_test_algebraic_kernel_1.h>
#include <CGAL/array.h>

#ifndef CGAL_TEST_ALGEBRAIC_KERNEL_2_H
#define CGAL_TEST_ALGEBRAIC_KERNEL_2_H

// For convenience
#define pow(a,b) CGAL::ipower(a,b)

namespace CGAL {

template< class AlgebraicKernel_2  >
void test_algebraic_kernel_2(const AlgebraicKernel_2& ak_2) {

  typedef AlgebraicKernel_2 AK_2;

  // AK_2 is also AK_1, thus test it:
  CGAL::test_algebraic_kernel_1<AK_2>(ak_2);

  typedef typename AK_2::Coefficient Coefficient;
  typedef typename AK_2::Bound Bound;
  typedef std::pair<Bound,Bound> BInterval;
  typedef std::array<Bound, 4> BArray;
  typedef typename AK_2::Polynomial_1 Polynomial_1;
  typedef typename AK_2::Polynomial_2 Polynomial_2;
  typedef typename AK_2::Algebraic_real_1 Algebraic_real_1;
  typedef typename AK_2::Algebraic_real_2 Algebraic_real_2;
  typedef typename AK_2::size_type size_type;
  typedef typename AK_2::Multiplicity_type Multiplicity_type;

  typedef CGAL::Polynomial_traits_d<Polynomial_1> Polynomial_traits_1;
  typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;


  #define CGAL_GET_FTOR(Name,name)                        \
  typedef typename AK_2::Name Name;      \
  const Name name = ak_2.name##_object();

  CGAL_GET_FTOR(Construct_algebraic_real_1,construct_algebraic_real_1);
  CGAL_GET_FTOR(Compare_1,compare_1);
  CGAL_GET_FTOR(Is_zero_at_1,is_zero_at_1);
  CGAL_GET_FTOR(Solve_1,solve_1);

  CGAL_GET_FTOR(Construct_algebraic_real_2,construct_algebraic_real_2);
  CGAL_GET_FTOR(Is_square_free_2,is_square_free_2);
  CGAL_GET_FTOR(Make_square_free_2,make_square_free_2);
  CGAL_GET_FTOR(Square_free_factorize_2,square_free_factorize_2);
  CGAL_GET_FTOR(Is_coprime_2,is_coprime_2);
  CGAL_GET_FTOR(Make_coprime_2,make_coprime_2);
  CGAL_GET_FTOR(Solve_2,solve_2);
  CGAL_GET_FTOR(Number_of_solutions_2,number_of_solutions_2);
  CGAL_GET_FTOR(Sign_at_2,sign_at_2);
  CGAL_GET_FTOR(Is_zero_at_2,is_zero_at_2);
  CGAL_GET_FTOR(Compute_x_2,compute_x_2);
  CGAL_GET_FTOR(Compute_y_2,compute_y_2);
  CGAL_GET_FTOR(Compute_polynomial_x_2,compute_polynomial_x_2);
  CGAL_GET_FTOR(Compute_polynomial_y_2,compute_polynomial_y_2);
  CGAL_GET_FTOR(Isolate_x_2,isolate_x_2);
  CGAL_GET_FTOR(Isolate_y_2,isolate_y_2);
  CGAL_GET_FTOR(Isolate_2,isolate_2);
  CGAL_GET_FTOR(Compare_x_2,compare_x_2);
  CGAL_GET_FTOR(Compare_y_2,compare_y_2);
  CGAL_GET_FTOR(Compare_xy_2,compare_xy_2);
  CGAL_GET_FTOR(Bound_between_x_2,bound_between_x_2);
  CGAL_GET_FTOR(Bound_between_y_2,bound_between_y_2);
  CGAL_GET_FTOR(Approximate_absolute_x_2,approximate_absolute_x_2);
  CGAL_GET_FTOR(Approximate_absolute_y_2,approximate_absolute_y_2);
  CGAL_GET_FTOR(Approximate_relative_x_2,approximate_relative_x_2);
  CGAL_GET_FTOR(Approximate_relative_y_2,approximate_relative_y_2);
#undef CGAL_GET_FTOR

  #define CGAL_CHECK_UFUNCTION(Name,AT,RT)                        \
  {                                                             \
    typedef typename Name::argument_type AT_;                   \
    typedef typename Name::result_type   RT_;                   \
    CGAL_USE_TYPE(AT_);                                           \
    CGAL_USE_TYPE(RT_);                                           \
    {static_assert(::std::is_same<AT,AT_>::value);}  \
    {static_assert(::std::is_same<RT,RT_>::value);}  \
  }
  #define CGAL_CHECK_BFUNCTION(Name,AT1,AT2,RT)                           \
  {                                                                     \
    typedef typename Name::first_argument_type AT1_;                    \
    typedef typename Name::second_argument_type AT2_;                   \
    typedef typename Name::result_type   RT_;                           \
    CGAL_USE_TYPE(AT1_);                                                 \
    CGAL_USE_TYPE(AT2_);                                                \
    CGAL_USE_TYPE(RT_);                                                 \
    {static_assert(::std::is_same<AT1,AT1_>::value);}        \
    {static_assert(::std::is_same<AT2,AT2_>::value);}        \
    {static_assert(::std::is_same<RT,RT_>::value);}          \
  }


  static_assert(::std::is_same
                         <Algebraic_real_2,
                         typename Construct_algebraic_real_2::result_type>
                           ::value);
  CGAL_CHECK_UFUNCTION(Is_square_free_2,Polynomial_2,bool);
  CGAL_CHECK_UFUNCTION(Make_square_free_2,Polynomial_2,Polynomial_2);
  // TODO: missing check for Square_free_factorize_2
  CGAL_CHECK_BFUNCTION(Is_coprime_2,Polynomial_2,Polynomial_2,bool);
  static_assert(::std::is_same
                         <bool,typename Make_coprime_2::result_type>::value);
  CGAL_CHECK_BFUNCTION(Number_of_solutions_2,Polynomial_2,Polynomial_2,
                       size_type);
  CGAL_CHECK_UFUNCTION(Compute_x_2,Algebraic_real_2,Algebraic_real_1);
  CGAL_CHECK_UFUNCTION(Compute_y_2,Algebraic_real_2,Algebraic_real_1);
  CGAL_CHECK_UFUNCTION(Compute_polynomial_x_2,Algebraic_real_2,Polynomial_1);
  CGAL_CHECK_UFUNCTION(Compute_polynomial_y_2,Algebraic_real_2,Polynomial_1);
  CGAL_CHECK_BFUNCTION(Isolate_x_2,Algebraic_real_2,Polynomial_1,BInterval);
  CGAL_CHECK_BFUNCTION(Isolate_y_2,Algebraic_real_2,Polynomial_1,BInterval);
  static_assert(::std::is_same
                         < BArray,typename Isolate_2::result_type>::value);
  CGAL_CHECK_BFUNCTION(Sign_at_2,Polynomial_2,Algebraic_real_2,Sign);
  CGAL_CHECK_BFUNCTION(Is_zero_at_2,Polynomial_2,Algebraic_real_2,bool);
  CGAL_CHECK_BFUNCTION(Compare_x_2,Algebraic_real_2,Algebraic_real_2,Sign);
  CGAL_CHECK_BFUNCTION(Compare_y_2,Algebraic_real_2,Algebraic_real_2,Sign);
  CGAL_CHECK_BFUNCTION(Compare_xy_2,Algebraic_real_2,Algebraic_real_2,Sign);
  CGAL_CHECK_BFUNCTION(Bound_between_x_2,Algebraic_real_2,
                       Algebraic_real_2,Bound);
  CGAL_CHECK_BFUNCTION(Bound_between_y_2,Algebraic_real_2,
                       Algebraic_real_2,Bound);
  CGAL_CHECK_BFUNCTION(Approximate_absolute_x_2,Algebraic_real_2,
                       int,BInterval);
  CGAL_CHECK_BFUNCTION(Approximate_absolute_y_2,Algebraic_real_2,
                       int,BInterval);
  CGAL_CHECK_BFUNCTION(Approximate_relative_x_2,Algebraic_real_2,
                       int,BInterval);
  CGAL_CHECK_BFUNCTION(Approximate_relative_y_2,Algebraic_real_2,
                       int,BInterval);
#undef CGAL_CHECK_BFUNCTION
#undef CGAL_CHECK_UFUNCTION

  typename Polynomial_traits_1::Construct_polynomial construct_polynomial_1;
  std::pair<CGAL::Exponent_vector,Coefficient> coeffs_t[1]
    = {std::make_pair(CGAL::Exponent_vector(1),Coefficient(1))};
  Polynomial_1 t=construct_polynomial_1(coeffs_t,coeffs_t+1);

  typename Polynomial_traits_2::Construct_polynomial construct_polynomial_2;
  std::pair<CGAL::Exponent_vector,Coefficient> coeffs_x[1]
    = {std::make_pair(CGAL::Exponent_vector(1,0),Coefficient(1))};
  Polynomial_2 x=construct_polynomial_2(coeffs_x,coeffs_x+1);
  std::pair<CGAL::Exponent_vector,Coefficient> coeffs_y[1]
    = {std::make_pair(CGAL::Exponent_vector(0,1),Coefficient(1))};
  Polynomial_2 y=construct_polynomial_2(coeffs_y,coeffs_y+1);
  Polynomial_2 one=construct_polynomial_2(Coefficient(1));

  {
    // Construct_algebraic_real_2
    Algebraic_real_2 ar1 = construct_algebraic_real_2(1,1);
    Algebraic_real_2 ar2 = construct_algebraic_real_2(Bound(1),Bound(1));
    Algebraic_real_2 ar3
      = construct_algebraic_real_2(Coefficient(2),Coefficient(-2));
    Algebraic_real_2 ar4
      = construct_algebraic_real_2(construct_algebraic_real_1(Coefficient(2)),
                                   construct_algebraic_real_1(Bound(-5)));
    // 3y+2x+2
    Polynomial_2 pol1_2 = y*3+x*2+2;
    // x^2+xy
    Polynomial_2 pol2_2 = pow(x,2)+x*y;
    Algebraic_real_2 ar5
      = construct_algebraic_real_2(pol1_2,pol2_2,0); // (0,-2/3)

    Algebraic_real_2 ar6
      = construct_algebraic_real_2(pol1_2,pol2_2,1); // (2,-2)

    Algebraic_real_2 ar7
      = construct_algebraic_real_2(pol1_2,pol2_2,
                                   Bound(-3),Bound(3),Bound(-1),Bound(100));

    // Comparisons
    assert(compare_xy_2(ar1,ar2)==CGAL::EQUAL);
    assert(compare_xy_2(ar1,ar3)==CGAL::SMALLER);
    assert(compare_xy_2(ar2,ar3)==CGAL::SMALLER);
    assert(compare_y_2(ar1,ar3)==CGAL::LARGER);
    assert(compare_x_2(ar3,ar4)==CGAL::EQUAL);
    assert(compare_y_2(ar3,ar4)==CGAL::LARGER);
    assert(compare_xy_2(ar3,ar4)==CGAL::LARGER);
    assert(compare_xy_2(ar3,ar6)==CGAL::EQUAL);
    assert(compare_x_2(ar1,ar5)==CGAL::LARGER);
    assert(compare_x_2(ar5,0)==CGAL::EQUAL);
    assert(compare_x_2(ar5,Bound(0))==CGAL::EQUAL);
    assert(compare_x_2(ar5,Coefficient(0))==CGAL::EQUAL);
    assert(compare_x_2(ar5,construct_algebraic_real_1(0))==CGAL::EQUAL);
    assert(compare_x_2(0,ar5)==CGAL::EQUAL);
    assert(compare_x_2(Bound(0),ar5)==CGAL::EQUAL);
    assert(compare_x_2(Coefficient(0),ar5)==CGAL::EQUAL);
    assert(compare_x_2(construct_algebraic_real_1(0),ar5)==CGAL::EQUAL);
    assert(compare_y_2(ar5,0)==CGAL::SMALLER);
    assert(compare_y_2(ar5,Bound(0))==CGAL::SMALLER);
    assert(compare_y_2(ar5,Coefficient(0))==CGAL::SMALLER);
    assert(compare_y_2(ar5,construct_algebraic_real_1(0))
                   ==CGAL::SMALLER);
    assert(compare_y_2(0,ar5)==CGAL::LARGER);
    assert(compare_y_2(Bound(0),ar5)==CGAL::LARGER);
    assert(compare_y_2(Coefficient(0),ar5)==CGAL::LARGER);
    assert(compare_y_2(construct_algebraic_real_1(0),ar5)==CGAL::LARGER);
    assert(compare_xy_2(ar5,0,0)==CGAL::SMALLER);
    assert(compare_xy_2(ar5,Bound(0),Bound(0))==CGAL::SMALLER);
    assert(compare_xy_2(ar5,Coefficient(0),Coefficient(0))
                   ==CGAL::SMALLER);
    assert(compare_xy_2(ar5,construct_algebraic_real_1(0),
                                construct_algebraic_real_1(0))==CGAL::SMALLER);
    assert(compare_xy_2(0,0,ar5)==CGAL::LARGER);
    assert(compare_xy_2(Bound(0),Bound(0),ar5)==CGAL::LARGER);
    assert(compare_xy_2(Coefficient(0),Coefficient(0),ar5)
                   ==CGAL::LARGER);
    assert(compare_xy_2(construct_algebraic_real_1(0),
                                construct_algebraic_real_1(0),
                                ar5)==CGAL::LARGER);
    assert(compare_y_2(0,ar5)==CGAL::LARGER);
    assert(compare_y_2(-1,ar5)==CGAL::SMALLER);
    assert(compare_xy_2(ar5,-1,1000)==CGAL::LARGER);
    assert(compare_xy_2(ar5,ar7)==CGAL::EQUAL);

    // Compute_x/y_2
    Algebraic_real_1 ar1_x = compute_x_2(ar1);
    assert(compare_1(1,ar1_x)==CGAL::EQUAL);
    assert(compare_1(ar1_x,1)==CGAL::EQUAL);
    Algebraic_real_1 ar6_y = compute_y_2(ar6);
    assert(compare_1(Bound(-2),ar6_y)==CGAL::EQUAL);
    assert(compare_1(ar6_y,Bound(-2))==CGAL::EQUAL);

    // SignAt_2, IsZeroAt_2
    assert(sign_at_2(pol1_2,ar4)==CGAL::NEGATIVE);
    assert(! is_zero_at_2(pol1_2,ar4));
    assert(sign_at_2(pol2_2,ar7)==CGAL::ZERO);
    assert(is_zero_at_2(pol2_2,ar7));
    assert(sign_at_2(pol2_2,ar1)==CGAL::POSITIVE);
    assert(! is_zero_at_2(pol1_2,ar4));

    // ComputePolynomial
    Polynomial_1 pol1_x = compute_polynomial_x_2(ar1);
    assert(is_zero_at_1(pol1_x,compute_x_2(ar1)));
    Polynomial_1 pol1_y = compute_polynomial_y_2(ar1);
    assert(is_zero_at_1(pol1_y,compute_y_2(ar1)));

    // Is_square_free, Make_square_free
    assert(is_square_free_2(pol1_2));
    assert(! is_square_free_2(pol1_2*pol1_2));
    Polynomial_2 pol1_2_prime = make_square_free_2(pol1_2*pol1_2);
    assert(CGAL::canonicalize(pol1_2)
                   ==CGAL::canonicalize(pol1_2_prime));
    // Is_coprime, Make_coprime
    assert(is_coprime_2(pol1_2,pol2_2));
    assert(is_coprime_2(pol1_2*Coefficient(3),pol2_2*Coefficient(3)));
    assert(! is_coprime_2(pol1_2,pol1_2*pol2_2));
    Polynomial_2 g,q1,q2;
    bool check = make_coprime_2(pol1_2,pol1_2*pol2_2,g,q1,q2);
    assert(! check);
    assert(CGAL::canonicalize(q1)
                   ==CGAL::canonicalize(one));
    assert(CGAL::canonicalize(q2)==CGAL::canonicalize(pol2_2));
    assert(CGAL::canonicalize(g)==CGAL::canonicalize(pol1_2));
    check = make_coprime_2(pol2_2,pol1_2,g,q1,q2);
    assert(check);
    assert(CGAL::canonicalize(q1)==CGAL::canonicalize(pol2_2));
    assert(CGAL::canonicalize(q2)==CGAL::canonicalize(pol1_2));
    assert(CGAL::total_degree(g)==0);
    // Test coprime also for factors in x only
    Polynomial_2 pol3_2 = pow(x,2)-3;
    assert(! is_coprime_2(pol3_2*pol1_2, pol3_2*pol3_2*pol2_2));
    check = make_coprime_2(pol3_2*pol1_2, pol3_2*pol3_2*pol2_2,g,q1,q2);
    assert(! check);
    assert(CGAL::canonicalize(g)==CGAL::canonicalize(pol3_2));
    assert(CGAL::canonicalize(q1)==CGAL::canonicalize(pol1_2));
    assert(CGAL::canonicalize(q2)==CGAL::canonicalize(pol3_2*pol2_2));
    // Square_free_factorize_2
    std::vector<std::pair<Polynomial_2,Multiplicity_type> > sqfr_factors;
    Polynomial_2 sqfr_fac_test_pol
      = pol1_2*CGAL::ipower(pol2_2,3)*CGAL::ipower(pol3_2,7);
    square_free_factorize_2(sqfr_fac_test_pol,
                            std::back_inserter(sqfr_factors));
    Polynomial_2 sqfr_fac_checksum=one;
    for(unsigned int i=0;i<sqfr_factors.size();i++) {
      assert(is_square_free_2(sqfr_factors[i].first));
      sqfr_fac_checksum
        *=CGAL::ipower(sqfr_factors[i].first, sqfr_factors[i].second);
    }
    assert(CGAL::canonicalize(sqfr_fac_checksum) ==
                   CGAL::canonicalize(sqfr_fac_test_pol));

    // BoundBetween
    Bound b = bound_between_x_2(ar1,ar3);
    assert(compare_1(compute_x_2(ar1),b)==CGAL::SMALLER);
    assert(compare_1(b,compute_x_2(ar3))==CGAL::SMALLER);
    b = bound_between_x_2(ar3,ar1);
    assert(compare_1(compute_x_2(ar1),b)==CGAL::SMALLER);
    assert(compare_1(b,compute_x_2(ar3))==CGAL::SMALLER);
    b = bound_between_y_2(ar5,ar6);
    assert(compare_1(compute_y_2(ar5),b)==CGAL::LARGER);
    assert(compare_1(b,compute_y_2(ar6))==CGAL::LARGER);
    b = bound_between_y_2(ar6,ar5);
    assert(compare_1(compute_y_2(ar5),b)==CGAL::LARGER);
    assert(compare_1(b,compute_y_2(ar6))==CGAL::LARGER);

  }

  {
    // Solve_2, Number_of_Solutions_2
    assert(number_of_solutions_2(one,x)==0);
    assert(number_of_solutions_2(one,y)==0);
    assert(number_of_solutions_2(y,x)==1);

    Polynomial_2 pol1 = pow(x,3)*pow(y,3)+4*x*y+pow(y,8)-x;
    Polynomial_2 pol2 = -pow(x,10)+6*y*pow(x,5)-10*x+5*y-3;
    Polynomial_2 pol3 = y*(x+2);

    assert(is_square_free_2(pol1));
    assert(is_square_free_2(pol2));
    assert(is_coprime_2(pol1,pol2));
    assert(number_of_solutions_2(pol1,pol2)==4);
    std::vector<std::pair<Algebraic_real_2,Multiplicity_type> > roots,roots2;
    solve_2(pol1,pol2,std::back_inserter(roots));
    assert(roots.size()==4);
    for(size_type i = 0; i < 4; i++) {
      assert(roots[i].second==1);
    }
    solve_2(pol1,pol2,Bound(-7),Bound(-1),Bound(-123),Bound(1),
            std::back_inserter(roots2));
    assert(roots2.size()==1);
    solve_2(pol1,pol2,Bound(-7),Bound(17),Bound(1),Bound(2),
            std::back_inserter(roots2));
    assert(roots2.size()==2); // sum of the roots so far
    solve_2(pol1,pol2,Bound(-1),Bound(8),Bound(-9),Bound(20),
            std::back_inserter(roots2));
    assert(roots2.size()==4); // sum of the roots so far
    std::sort(roots.begin(),roots.end());
    std::sort(roots2.begin(),roots2.end());
    assert(std::equal(roots.begin(),roots.end(),roots2.begin()));

    assert(number_of_solutions_2(pol1,pol2*pol3)==7);
    roots2.clear();
    solve_2(pol1,pol2*pol3,Bound(-2),Bound(0),Bound(0),Bound(1),
            std::back_inserter(roots2));
    assert(roots2.size()==5);

    // Isolate_2
    Polynomial_2 circle = pow(x,2)+pow(y,2)-1;
    Algebraic_real_2 ar1 = construct_algebraic_real_2(10*x+1,10*y-1,0);
    BArray box = isolate_2(ar1,circle);
    assert(box[0]>-1);
    assert(box[1]<1);
    assert(box[2]>-1);
    assert(box[3]<1);
    assert(sign_at_2(circle,construct_algebraic_real_2(box[0],box[2]))
                   == CGAL::NEGATIVE);
    assert(sign_at_2(circle,construct_algebraic_real_2(box[1],box[2]))
                   == CGAL::NEGATIVE);
    assert(sign_at_2(circle,construct_algebraic_real_2(box[0],box[3]))
                   == CGAL::NEGATIVE);
    assert(sign_at_2(circle,construct_algebraic_real_2(box[1],box[3]))
                   == CGAL::NEGATIVE);
    isolate_2(ar1,one);
    box = isolate_2(ar1,pol1,pol2*pol3);
    roots2.clear();
    solve_2(pol1,pol2*pol3,box[0],box[1],box[2],box[3],
            std::back_inserter(roots2));
    assert(roots2.size()==0);
    roots.clear();
    roots2.clear();
    solve_2(pol1,pol2*pol3,std::back_inserter(roots));
    assert(roots.size()==7);
    for(size_type i = 0; i < 7; i++) {
      box=isolate_2(roots[i].first,pol1,pol2*pol3);
      solve_2(pol1,pol2*pol3,box[0],box[1],box[2],box[3],
              std::back_inserter(roots2));
      assert(roots2.size()==static_cast<size_t>(i+1));
    }
    std::sort(roots.begin(),roots.end());
    std::sort(roots2.begin(),roots2.end());
    assert(std::equal(roots.begin(),roots.end(),roots2.begin()));

    // Isolate_x_2, and Isolate_y_2
    Polynomial_1 q = pow(t,8) - 7*pow(t,6)+3*pow(t,3)-2*pow(t,2)+1;
    assert(!is_zero_at_1(q,compute_x_2(ar1)));
    BInterval bounds = isolate_x_2(ar1,q);
    assert(bounds.first<=bounds.second);
    std::vector<Algebraic_real_1> roots_1,roots2_1,roots3_1;
    solve_1(q,false,bounds.first,bounds.second,std::back_inserter(roots_1));
    std::vector<std::pair<Algebraic_real_1,Multiplicity_type> > root_pairs_1;
    solve_1(q,bounds.first,bounds.second,std::back_inserter(root_pairs_1));
    assert(roots_1.size()==0);
    solve_1(q,false,std::back_inserter(roots_1));
    assert(roots_1.size()==4);
    for(size_type i=0;i<4;i++) {
      Algebraic_real_2 c_ar
        = construct_algebraic_real_2(roots_1[i],roots_1[3-i]);
      bounds=isolate_x_2(c_ar,q);
      solve_1(q,false,bounds.first,bounds.second,std::back_inserter(roots2_1));
      bounds=isolate_y_2(c_ar,q);
      solve_1(q,false,bounds.first,bounds.second,std::back_inserter(roots3_1));
    }
    assert(roots2_1.size()==roots3_1.size());
    assert(roots_1.size()==roots3_1.size());
    assert(std::equal(roots_1.begin(),roots_1.end(),roots2_1.begin()));
    assert(std::equal(roots_1.begin(),roots_1.end(),
                              roots3_1.rbegin()));

    // Approximate_absolute_2, Approximate_relative_2
    Polynomial_2 pol4 = pow(x,5)+pow(y,2)-y+y*x-1;
    Polynomial_2 pol5 = pow(y,2)-pow(x,3);
    // absolute
    Algebraic_real_2 ar2
      = construct_algebraic_real_2(pol4,pol5,
                                   Bound(0),Bound(2),Bound(0),Bound(3));
    for(int i=1;i<1024;i*=2) {
      Bound dist = CGAL::ipower(Bound(2),i);
      bounds = approximate_absolute_x_2(ar2,i);
      assert(bounds.first<=bounds.second);
      assert(compare_x_2(bounds.first,ar2)!=CGAL::LARGER);
      assert(compare_x_2(bounds.second,ar2)!=CGAL::SMALLER);
      // NOTE: The following check is not sufficient for ensuring
      //       the postcondition, but if it fails, the postcondition
      //       cannot be satisfied
      assert(dist*(bounds.second-bounds.first)<=Bound(2));
      // For the default implementation, we can test this - this
      // is sufficient for ensuring the postcondition, but other
      // kernels might not satisfy this property
      assert(dist*(bounds.second-bounds.first)<=Bound(1));
      bounds = approximate_absolute_y_2(ar2,i);
      assert(bounds.first<=bounds.second);
      assert(compare_y_2(bounds.first,ar2)!=CGAL::LARGER);
      assert(compare_y_2(bounds.second,ar2)!=CGAL::SMALLER);
      // NOTE: The following check is not sufficient for ensuring
      //       the postcondition, but if it fails, the postcondition
      //       cannot be satisfied
      assert(dist*(bounds.second-bounds.first)<=Bound(2));
      // For the default implementation, we can test this - this
      // is sufficient for ensuring the postcondition, but other
      // kernels might not satisfy this property
      assert(dist*(bounds.second-bounds.first)<=Bound(1));
    }

    // relative
    ar2 = construct_algebraic_real_2(pol4,pol5,
                                     Bound(0),Bound(2),Bound(-3),Bound(0));
    for(int i=1;i<1024;i*=2) {
      Bound dist = CGAL::ipower(Bound(2),i);
      bounds = approximate_relative_x_2(ar2,i);
      Bound min=(CGAL::min)(CGAL::abs(bounds.first),CGAL::abs(bounds.second));
      Bound max=(CGAL::max)(CGAL::abs(bounds.first),CGAL::abs(bounds.second));
      assert(bounds.first<=bounds.second);
      assert(CGAL::sign(bounds.first)==CGAL::sign(bounds.second));
      assert(compare_x_2(bounds.first,ar2)!=CGAL::LARGER);
      assert(compare_x_2(bounds.second,ar2)!=CGAL::SMALLER);
      // NOTE: The following check is not sufficient for ensuring
      //       the postcondition, but if it fails, the postcondition
      //       cannot be satisfied
      assert(dist*(bounds.second-bounds.first)<=Bound(2)*max);
      // For the default implementation, we can test this - this
      // is sufficient for ensuring the postcondition, but other
      // kernels might not satisfy this property
      assert(dist*(bounds.second-bounds.first)<=Bound(1)*min);
      bounds = approximate_relative_y_2(ar2,i);
      assert(bounds.first<=bounds.second);
      assert(compare_y_2(bounds.first,ar2)!=CGAL::LARGER);
      assert(compare_y_2(bounds.second,ar2)!=CGAL::SMALLER);
      // NOTE: The following check is not sufficient for ensuring
      //       the postcondition, but if it fails, the postcondition
      //       cannot be satisfied
      assert(dist*(bounds.second-bounds.first)<=Bound(2)*max);
      // For the default implementation, we can test this - this
      // is sufficient for ensuring the postcondition, but other
      // kernels might not satisfy this property
      assert(dist*(bounds.second-bounds.first)<=Bound(1)*min);
    }

  }




}

} //namespace CGAL

#endif
