// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

#include <cassert>
#include <CGAL/use.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_coercion_traits.h>
#include "_test_ak_polynomial_traits_d.h"


// Test for the Algebraic_kernel syntax
#ifndef CGAL_TEST_ALGEBRAIC_KERNEL_1_H
#define CGAL_TEST_ALGEBRAIC_KERNEL_1_H

namespace CGAL{

template <class AlgebraicKernel_d_1>
void test_algebraic_kernel_1(const AlgebraicKernel_d_1& ak_1){
  typedef AlgebraicKernel_d_1 Algebraic_kernel_d_1;
  CGAL_USE_TYPE(Algebraic_kernel_d_1);

  typedef typename AlgebraicKernel_d_1::Coefficient Coefficient;
  typedef typename AlgebraicKernel_d_1::Polynomial_1 Polynomial_1;
  typedef typename AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1;
  typedef typename AlgebraicKernel_d_1::Bound Bound;
  typedef std::pair<Bound,Bound> BInterval;
  typedef CGAL::Polynomial_traits_d<Polynomial_1> PT;

  {
    // check Coefficient
    typedef Algebraic_structure_traits<Coefficient> AST;
    typedef typename AST::Is_exact Is_exact;
    typedef typename AST::Algebraic_category Algebraic_category;
    test_algebraic_structure< Coefficient,Algebraic_category,Is_exact>();
    test_real_embeddable<Coefficient>();
  }
  {
    // check Polynomial_1
    typedef Algebraic_structure_traits<Coefficient> CAST;
    typedef typename CAST::Is_exact Is_exact;
    typedef typename CAST::Algebraic_category Algebraic_category;
    CGAL_USE_TYPE(Algebraic_category);
    test_ak_polynomial_traits_d<Is_exact,PT>()(PT());

    // test not possible due to bug in test_algebraic_structure
    // div(3,2)=3/2 != 0 in case of Polynomial<Rational>
    //  typedef Algebraic_structure_traits<Polynomial_1> AST;
    //  typedef typename AST::Algebraic_category Algebraic_category;
    //  test_algebraic_structure< Polynomial_1,Algebraic_category,Tag_true>();
  }
  {
    // check Algebraic_real_1
    test_real_embeddable<Algebraic_real_1>();
  }

  {
    typedef Algebraic_structure_traits<Bound> AST;
    typedef typename AST::Is_exact Is_exact;
    typedef typename AST::Algebraic_category Algebraic_category;
    test_algebraic_structure< Bound,Algebraic_category,Is_exact>();
    test_real_embeddable<Bound>();
  }

  test_explicit_interoperable_from_to<int, Coefficient>();
  test_explicit_interoperable_from_to<int, Bound>();

  // interoperability was removed from algebraic real concept
  //test_explicit_interoperable_from_to<int        , Algebraic_real_1>();
  //test_explicit_interoperable_from_to<Bound      , Algebraic_real_1>();
  //test_explicit_interoperable_from_to<Coefficient, Algebraic_real_1>();

#define CGAL_GET_FTOR(Name,name)                        \
  typedef typename AlgebraicKernel_d_1::Name Name;      \
  const Name name = ak_1.name##_object();


  CGAL_GET_FTOR(Construct_algebraic_real_1,construct_algebraic_real_1);
  CGAL_GET_FTOR(Compute_polynomial_1,compute_polynomial_1);
  CGAL_GET_FTOR(Is_square_free_1,is_square_free_1);
  CGAL_GET_FTOR(Make_square_free_1,make_square_free_1);
  CGAL_GET_FTOR(Square_free_factorize_1,square_free_factorize_1);
  CGAL_GET_FTOR(Is_coprime_1,is_coprime_1);
  CGAL_GET_FTOR(Make_coprime_1,make_coprime_1);
  CGAL_GET_FTOR(Solve_1,solve_1);
  CGAL_GET_FTOR(Number_of_solutions_1,number_of_solutions_1);
  CGAL_GET_FTOR(Sign_at_1,sign_at_1);
  CGAL_GET_FTOR(Is_zero_at_1,is_zero_at_1);
  CGAL_GET_FTOR(Compare_1,compare_1);
  CGAL_GET_FTOR(Bound_between_1,bound_between_1);
  CGAL_GET_FTOR(Approximate_absolute_1,approximate_absolute_1);
  CGAL_GET_FTOR(Approximate_relative_1,approximate_relative_1);
#undef CGAL_GET_FTOR

#define CGAL_CHECK_UFUNCTION(Name,AT,RT)                          \
  {                                                               \
    typedef typename Name::argument_type AT_;                     \
    typedef typename Name::result_type   RT_;                     \
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
    CGAL_USE_TYPE(AT1_);                                                \
    CGAL_USE_TYPE(AT2_);                                                \
    CGAL_USE_TYPE(RT_);                                                 \
    {static_assert(::std::is_same<AT1,AT1_>::value);}      \
    {static_assert(::std::is_same<AT2,AT2_>::value);}      \
    {static_assert(::std::is_same<RT,RT_>::value);}        \
  }

  // TODO: missing check for Construct_algebraic_real_1
  CGAL_CHECK_UFUNCTION(Compute_polynomial_1,Algebraic_real_1,Polynomial_1);
  CGAL_CHECK_UFUNCTION(Is_square_free_1,Polynomial_1,bool);
  CGAL_CHECK_UFUNCTION(Make_square_free_1,Polynomial_1,Polynomial_1);
  // TODO: missing check for Square_free_factorize_1
  CGAL_CHECK_BFUNCTION(Is_coprime_1,Polynomial_1,Polynomial_1,bool);
  // TODO: missing check for Make_coprime_1
  // TODO: missing check for Solve_1
  CGAL_CHECK_UFUNCTION(Number_of_solutions_1,Polynomial_1,int);
  CGAL_CHECK_BFUNCTION(Sign_at_1,Polynomial_1,Algebraic_real_1,Sign);
  CGAL_CHECK_BFUNCTION(Is_zero_at_1,Polynomial_1,Algebraic_real_1,bool);
  CGAL_CHECK_BFUNCTION(Compare_1,Algebraic_real_1,Algebraic_real_1,Sign);
  CGAL_CHECK_BFUNCTION(Bound_between_1,Algebraic_real_1,Algebraic_real_1,Bound);
  CGAL_CHECK_BFUNCTION(Approximate_absolute_1,Algebraic_real_1,int,BInterval);
  CGAL_CHECK_BFUNCTION(Approximate_relative_1,Algebraic_real_1,int,BInterval);
#undef CGAL_CHECK_BFUNCTION
#undef CGAL_CHECK_UFUNCTION

  Polynomial_1 x = typename PT::Shift()(Polynomial_1(1),1);
  {
    assert( is_square_free_1(ipower((x-1),1)));
    assert(!is_square_free_1(ipower((x-1),2)));
  }
  {
    assert( make_square_free_1(ipower(5*(x-1),2))==ipower((x-1),1));
  }
  {
    std::list< std::pair<Polynomial_1,int> > factors;
    square_free_factorize_1((x-1)*(x-2)*(x-2),std::back_inserter(factors));
    assert(factors.size()==2);
    assert(factors.front() != factors.back());
    assert(
        factors.front() == std::make_pair((x-1),1) ||
        factors.front() == std::make_pair((x-2),2) );
    assert(
        factors.back()  == std::make_pair((x-1),1) ||
        factors.back()  == std::make_pair((x-2),2) );
  }

  assert( is_coprime_1((x-1),(x-2)));
  assert(!is_coprime_1((x-1)*(x-2),(x-1)*(x-3)));

  {
    Polynomial_1 a,b,c,d,e;
    a = (x-1)*(x-2);
    b = (x-1)*(x-3);
    assert(!make_coprime_1(a,b,c,d,e));
    assert( c == (x-1) ); // gcd
    assert( d == (x-2) );
    assert( e == (x-3) );

    a = (x-1);
    b = (x-2);
    assert( make_coprime_1(a,b,c,d,e) );
    assert( c == (1) ); // gcd
    assert( d == (x-1) );
    assert( e == (x-2) );
  }

  {
    // solve_1 for OI::value_type == std::pair<Algebraic_real_1,int>
    typedef  std::list<std::pair<Algebraic_real_1,unsigned int> > ROOTS;
    ROOTS roots;
    std::back_insert_iterator<ROOTS> biit =
      solve_1((x-1)*(x-2)*(x-2),std::back_inserter(roots));

    assert(roots.size()==2);
    assert(roots.front() != roots.back());
    assert(
        roots.front() == std::make_pair(Algebraic_real_1(1),(unsigned int)1) ||
        roots.front() == std::make_pair(Algebraic_real_1(2),(unsigned int)2) );
    assert(
        roots.back()  == std::make_pair(Algebraic_real_1(1),(unsigned int)1) ||
        roots.back()  == std::make_pair(Algebraic_real_1(2),(unsigned int)2) );

    solve_1((x-1)*(x-2)*(x-2),biit); // use iterator again
    assert(roots.size()==4);
  }
  {
    // Compute_polynomial
    typedef  std::vector<std::pair<Algebraic_real_1,unsigned int> > ROOTS;
    ROOTS roots;
    Polynomial_1 p1 = (x-1)*(x-2)*(x-2);
    // std::back_insert_iterator<ROOTS> biit =
      solve_1(p1,std::back_inserter(roots));
    Algebraic_real_1 ar = roots[1].first;
    Polynomial_1 p2 = compute_polynomial_1(ar);
    assert(!is_coprime_1(p1,p2));
    assert(is_zero_at_1(p2,ar));

  }
  {
    // solve_1 for OI::value_type == std::pair<Algebraic_real_1>
    typedef  std::list<Algebraic_real_1 > ROOTS;
    ROOTS roots;
    std::back_insert_iterator<ROOTS> biit =
      solve_1((x-1)*(x-2)*(x-2),false,std::back_inserter(roots));

    assert(roots.size()==2);
    assert(roots.front() != roots.back());
    assert(
        roots.front() == Algebraic_real_1(1) ||
        roots.front() == Algebraic_real_1(2) );
    assert(
        roots.back()  == Algebraic_real_1(1) ||
        roots.back()  == Algebraic_real_1(2) );

    solve_1((x-1)*(x-2),true,biit); // use iterator again
    assert(roots.size()==4);
  }
  {
    // number_of_solutions
    assert(3 == number_of_solutions_1((x-1)*(x-2)*(x-3)));
  }
  {
    assert(sign_at_1(x*0,Algebraic_real_1(0)) == ZERO);

    assert(sign_at_1(x-1,Algebraic_real_1(0)) == NEGATIVE);
    assert(sign_at_1(x-1,Algebraic_real_1(1)) == ZERO);
    assert(sign_at_1(x-1,Algebraic_real_1(2)) == POSITIVE);

    std::vector<Algebraic_real_1> roots;
    solve_1(x*x-2,true,std::back_inserter(roots));
    assert(sign_at_1((x+1),roots[0]) == CGAL::NEGATIVE);

    Polynomial_1 f = (x*x-2)*(x*x-2);
    assert(sign_at_1(f,roots[0]) == CGAL::ZERO);
  }
  {
    std::list<Algebraic_real_1> roots;
    solve_1((x*x-3),true,std::back_inserter(roots));
    Algebraic_real_1 root = (CGAL::min)(roots.front(),roots.back());
    assert(sign_at_1(x*x-2,root) == POSITIVE);
    assert(sign_at_1(x*x-3,root) == ZERO);
    assert(sign_at_1((x*x-3)*(x-4),root) == ZERO);
    assert(sign_at_1(x*x-4,root) == NEGATIVE);
  }
  {
    assert(is_zero_at_1(x*0,Algebraic_real_1(0)) == true);

    assert(is_zero_at_1(x-1,Algebraic_real_1(0)) == false);
    assert(is_zero_at_1(x-1,Algebraic_real_1(1)) == true);
    assert(is_zero_at_1(x-1,Algebraic_real_1(2)) == false);

    std::list<Algebraic_real_1> roots;
    solve_1((x*x-3),true,std::back_inserter(roots));
    Algebraic_real_1 root = (CGAL::min)(roots.front(),roots.back());
    assert(is_zero_at_1(x*x-2,root) == false);
    assert(is_zero_at_1(x*x-3,root) == true);
    assert(is_zero_at_1((x*x-3)*(x-4),root) == true);
    assert(is_zero_at_1(x*x-4,root) == false);
  }
  {
    assert(compare_1(Algebraic_real_1( 1),Algebraic_real_1( 2)) == SMALLER);
    assert(compare_1(Algebraic_real_1( 2),Algebraic_real_1( 2)) == EQUAL);
    assert(compare_1(Algebraic_real_1( 3),Algebraic_real_1( 2)) == LARGER);
    assert(compare_1(Algebraic_real_1(-1),Algebraic_real_1(-2)) == LARGER);
    assert(compare_1(Algebraic_real_1(-2),Algebraic_real_1(-2)) == EQUAL);
    assert(compare_1(Algebraic_real_1(-3),Algebraic_real_1(-2)) == SMALLER);
  }
  {
    Bound bound = bound_between_1(Algebraic_real_1(1),Algebraic_real_1(2));
    assert(compare_1(bound,Algebraic_real_1(1)) == LARGER  );
    assert(compare_1(bound,Algebraic_real_1(2)) == SMALLER );
  }

  CGAL::IO::set_pretty_mode(std::cerr);

  // Approximations
  bool all_right = true;
  std::cout << "start test_approximation (takes a while) " << std::flush;
  Coefficient c = CGAL::ipower(Coefficient(2), 2500) + 1;
  // we choose coefficients: small, large, (close to) power of two
  std::vector< int > coeffs;
  coeffs.push_back(1);
  coeffs.push_back(13);
  //coeffs.push_back(255);
  //coeffs.push_back(499);
  //coeffs.push_back(512);
  //coeffs.push_back(10000);
  //coeffs.push_back(3);
  //coeffs.push_back(7);
  //coeffs.push_back(64);
  //coeffs.push_back(100);
  coeffs.push_back(1023);
  //coeffs.push_back(4096);
  std::vector< int > precs;
  //precs.push_back(0);
  //precs.push_back(1);
  //precs.push_back(2);
  //precs.push_back(4);
  precs.push_back(8);
  //precs.push_back(13);
  //  precs.push_back(1023);
  //  precs.push_back(2048);
  //precs.push_back(53);
  //precs.push_back(3);
  //precs.push_back(64);
  //precs.push_back(106);
  //precs.push_back(424);
  for (typename std::vector< int >::const_iterator c0i = coeffs.begin();
       c0i != coeffs.end(); c0i++) {
    for (typename std::vector< int >::const_iterator c2i = coeffs.begin();
         c2i != coeffs.end(); c2i++) {
      // we basically test a quadratic polynomial (with chosen small and large
      // quadratic and constant coefficient, which is disturbed by a root close to zero).
      //Polynomial_1 poly((*c2i*x*x - *c0i) * (c*x-1));
      Polynomial_1 poly((*c2i*x*x - *c0i) * (c*x-1));
      std::list<Algebraic_real_1> roots;
      solve_1(poly,true,std::back_inserter(roots));
      for (typename std::vector< int >::const_iterator pi = precs.begin();
           pi != precs.end(); pi++) {
        // all three roots are approximated with various precisions
        long p = *pi;
        { // Approximate_absolute_1 with positive p
          for (typename std::list< Algebraic_real_1 >::const_iterator rit = roots.begin();
               rit != roots.end(); rit++) {
            BInterval bi = approximate_absolute_1(*rit,p);
            assert(compare_1(bi.first ,*rit) != LARGER );
            assert(compare_1(bi.second,*rit) != SMALLER);
            assert(CGAL::sign(bi.second - bi.first) != NEGATIVE);
            if (!((bi.second - bi.first) * (p == 0 ? Bound(1) : ipower(Bound(2),p-1))
                  <= (p == 0 ? Bound(2) : Bound(1)) )) {
              all_right = false;
              std::cerr << "ERROR: Approximate_absolute_1 fails for prec = " << p
                        << " of this root: " << *rit << std::endl;
            }
          }
        }
        { // Approximate_absolute_1 with negative p
          for (typename std::list< Algebraic_real_1 >::const_iterator rit = roots.begin();
               rit != roots.end(); rit++) {
            BInterval bi = approximate_absolute_1(*rit,-p);
            assert(compare_1(bi.first ,*rit) != LARGER );
            assert(compare_1(bi.second,*rit) != SMALLER);
            assert(CGAL::sign(bi.second - bi.first) != NEGATIVE);
            if (!((bi.second - bi.first) <= ipower(Bound(2),1-(-p)) )) {
              all_right = false;
              std::cerr << "ERROR: Approximate_absolute_1 fails for prec = " << -p
                        << " of this root: " << *rit << std::endl;
            }
          }
        }
        { // Approximate_relative_1 with positive p
          for (typename std::list< Algebraic_real_1 >::const_iterator rit = roots.begin();
               rit != roots.end(); rit++) {
            BInterval bi = approximate_relative_1(*rit,p);
            assert(compare_1(bi.first ,*rit) != LARGER );
            assert(compare_1(bi.second,*rit) != SMALLER);
            assert(CGAL::sign(bi.second - bi.first) != NEGATIVE);
            if (!((bi.second - bi.first) * (p == 0 ? Bound(1) : ipower(Bound(2),p-1))
                  <= (p == 0 ? Bound(2) : Bound(1)) * (CGAL::max)(abs(bi.first),abs(bi.second)))) {
              all_right = false;
              std::cerr << "ERROR: Approximate_relative_1 fails for prec = " << p
                        << " of this root: " << *rit << std::endl;

            }
          }
        }
        { // Approximate_relative_1 with negative p
          for (typename std::list< Algebraic_real_1 >::const_iterator rit = roots.begin();
               rit != roots.end(); rit++) {
            BInterval bi = approximate_relative_1(*rit,-p);
            assert(compare_1(bi.first ,*rit) != LARGER );
            assert(compare_1(bi.second,*rit) != SMALLER);
            assert(CGAL::sign(bi.second - bi.first) != NEGATIVE);
            if (!((bi.second - bi.first) <=
                  ipower(Bound(2),1-(-p)) * (CGAL::max)(abs(bi.first),abs(bi.second)))) {
              all_right = false;
              std::cerr << "ERROR: Approximate_relative_1 fails for prec = " << -p
                        << " of this root: " << *rit << std::endl;
            }
          }
        }
      }
    }
  }
  assert(all_right); // some approximation was not good enough
  std::cout << " ok" << std::endl;
  {

#define CGAL_TEST_ALGEBRAIC_REAL_IO(_f)         \
    alg1=_f;                                    \
    ss<<CGAL::IO::oformat(alg1);                \
    assert(ss.good());                          \
    ss>>CGAL::IO::iformat(alg2);                \
    assert(!ss.fail());                         \
    ss.clear();                                 \
    assert(alg1==alg2)
    // Note: after the reading ss>>CGAL::IO::iformat(alg2) the state of ss can
    // have the eofbit. The C++ norm says if one tries to write to a stream
    // with eofbit, then the failbit will be set. That is why one must
    // clear the iostate with ss.clear().


    Algebraic_real_1 alg1,alg2;
    std::stringstream ss;
    CGAL::IO::set_ascii_mode(ss);

    // test construction from int, Coefficient and Bound
    CGAL_TEST_ALGEBRAIC_REAL_IO(construct_algebraic_real_1(int(2)));
    CGAL_TEST_ALGEBRAIC_REAL_IO(construct_algebraic_real_1(Coefficient(2)));
    CGAL_TEST_ALGEBRAIC_REAL_IO(construct_algebraic_real_1(Bound(2)));

  // construction by index
    Polynomial_1 x = CGAL::shift(Polynomial_1(1),1); // the monom x
    CGAL_TEST_ALGEBRAIC_REAL_IO(construct_algebraic_real_1(x*x-2,1));

    // construction by isolating interval
    CGAL_TEST_ALGEBRAIC_REAL_IO(construct_algebraic_real_1(x*x-2,Bound(0),Bound(2)));
#undef CGAL_TEST_ALGEBRAIC_REAL_IO
  }
}

} // namespace CGAL

#endif //CGAL_TEST_ALGEBRAIC_KERNEL_1_H
