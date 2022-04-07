// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-2.1-only
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#include <CGAL/config.h>

#if defined(CGAL_USE_GMP) && defined(CGAL_USE_MPFI) && defined(CGAL_USE_RS)

#include "include/CGAL/_test_algebraic_kernel_1.h"
#ifdef CGAL_RS_TEST_LOG_TIME
#include <ctime>
#endif

// default RS_AK_1
#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>

// different isolators
#include <CGAL/RS/rs2_isolator_1.h>
#ifdef CGAL_USE_RS3
#include <CGAL/RS/rs23_k_isolator_1.h>
#endif

// different refiners
#include <CGAL/RS/bisection_refiner_1.h>
#ifdef CGAL_USE_RS3
#include <CGAL/RS/rs3_refiner_1.h>
#include <CGAL/RS/rs3_k_refiner_1.h>
#endif

template <class AK_>
int test_ak(){
  typedef AK_                                           AK;
  typedef typename AK::Polynomial_1                     Polynomial_1;
  //typedef typename AK::Coefficient                      Coefficient;
  typedef typename AK::Bound                            Bound;
  typedef typename AK::Algebraic_real_1                 Algebraic_real_1;
  typedef typename AK::Multiplicity_type                Multiplicity_type;

#ifdef CGAL_RS_TEST_LOG_TIME
  clock_t start=clock();
#endif

  AK ak; // an algebraic kernel object
  CGAL::test_algebraic_kernel_1<AK>(ak); // we run standard tests first

  typename AK::Solve_1 solve_1 = ak.solve_1_object();
  Polynomial_1 x = CGAL::shift(Polynomial_1(1),1);
  int returnvalue=0;

  // variant using a bool indicating a square free polynomial
  // multiplicities are not computed
  std::vector<Algebraic_real_1> roots;
  solve_1(x*x-2,true, std::back_inserter(roots));
  if(roots.size()!=2){
    returnvalue-=1;
    std::cerr<<"error 1: the number of roots of x^2-2 must be 2"<<
      std::endl;
  }
  if(-1.42>=roots[0] || -1.41<=roots[0] ||
      1.41>=roots[1] || 1.42<=roots[1]){
    returnvalue-=2;
    std::cerr<<"error 2: the roots of x^2-2 are wrong"<<std::endl;
  }
  roots.clear();

  // variant for roots in a given range of a square free polynomial
  solve_1((x*x-2)*(x*x-3),true, Bound(0),Bound(10),
      std::back_inserter(roots));
  if(roots.size()!=2){
    returnvalue-=4;
    std::cerr<<"error 3: the number of roots of (x^2-2)*(x^2-3)"<<
      " between 0 and 10 must be 2"<<std::endl;
  }
  if(1.41>=roots[0] || 1.42<=roots[0] ||
      1.73>=roots[1] || 1.74<=roots[1]){
    returnvalue-=8;
    std::cerr<<"error 4: the roots of (x^2-2)*(x^2-3)"<<
      " between 0 and 10 are wrong"<<std::endl;
  }
  roots.clear();

  // variant computing all roots with multiplicities
  std::vector<std::pair<Algebraic_real_1,Multiplicity_type> > mroots;
  solve_1((x*x-2), std::back_inserter(mroots));
  if(mroots.size()!=2){
    returnvalue-=16;
    std::cerr<<"error 5: the number of roots of x^2-2 must be 2"<<
      std::endl;
  }
  if(-1.42>=mroots[0].first || -1.41<=mroots[0].first ||
      1.41>=mroots[1].first || 1.42<=mroots[1].first){
    returnvalue-=32;
    std::cerr<<"error 6: the roots of x^2-2 are wrong"<<std::endl;
  }
  if(mroots[0].second!=1 && mroots[1].second!=1){
    returnvalue-=64;
    std::cerr<<"error 7: the multiplicities of the"<<
      " roots of x^2-2 are wrong"<<std::endl;
  }
  mroots.clear();

  // variant computing roots with multiplicities for a range
  solve_1((x*x-2)*(x*x-3),Bound(0),Bound(10),std::back_inserter(mroots));
  if(mroots.size()!=2){
    returnvalue-=128;
    std::cerr<<"error 8: the number of roots of (x^2-2)*(x^2-3)"<<
      " between 0 and 10 must be 2"<<std::endl;
  }
  if(1.41>=mroots[0].first || 1.42<=mroots[0].first ||
      1.73>=mroots[1].first || 1.74<=mroots[1].first){
    returnvalue-=256;
    std::cerr<<"error 9: the roots of (x^2-2)*(x^2-3) are wrong"<<
      std::endl;
  }
  if(mroots[0].second!=1 && mroots[1].second!=1){
    returnvalue-=512;
    std::cerr<<"error 10: the multiplicities of the roots of"<<
      " (x^2-2)*(x^2-3) are wrong"<<std::endl;
  }

  typename AK::Number_of_solutions_1 nos_1 = ak.number_of_solutions_1_object();
  if(nos_1(x*x*x-2)!=1){
    returnvalue-=1024;
    std::cerr<<"error 11: x^3-2 must have only one root"<<std::endl;
  }

#ifdef CGAL_RS_TEST_LOG_TIME
  std::cerr<<"*** test time: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<
          " seconds"<<std::endl;
#endif

  return returnvalue;
}

int main(){
        // We'll test three different RS-based univariate AK's:
        // - the default RS one,
        // - one with RS2 functions only, and
        // - one with both RS2 and RS3 functions.

        // the default RS kernel
        typedef CGAL::Algebraic_kernel_rs_gmpz_d_1              AK_default;

        typedef CGAL::Polynomial<CGAL::Gmpz>                    Polynomial;
        typedef CGAL::Gmpfr                                     Bound;

        // the RS2-only kernel
        typedef CGAL::RS2::RS2_isolator_1<Polynomial,Bound>     RS2_isolator;
        typedef CGAL::Bisection_refiner_1<Polynomial,Bound>     B_refiner;
        typedef CGAL::RS_AK1::Algebraic_kernel_1<Polynomial,
                                                 Bound,
                                                 RS2_isolator,
                                                 B_refiner>     AK_RS2;

#ifdef CGAL_USE_RS3
        // the RS2/RS3 kernel
        typedef CGAL::RS23_k_isolator_1<Polynomial,Bound>       K_isolator;
        typedef CGAL::RS3::RS3_k_refiner_1<Polynomial,Bound>    RS3_k_refiner;
        typedef CGAL::RS_AK1::Algebraic_kernel_1<Polynomial,
                                                 Bound,
                                                 K_isolator,
                                                 RS3_k_refiner> AK_RS2_RS3;
#endif // CGAL_USE_RS3

        // test all and return the result
        long result=0;
        std::cerr<<"*** testing default RS AK_1:";
        result+=test_ak<AK_default>();
        std::cerr<<"*** testing RS2 AK_1:";
        result+=(2048*test_ak<AK_RS2>());
#ifdef CGAL_USE_RS3
        std::cerr<<"*** testing RS2/RS3 k_AK_1:";
        result+=(4096*test_ak<AK_RS2_RS3>());
#endif // CGAL_USE_RS3
        std::cerr<<"*** result of the tests (should be 0): "<<result<<std::endl;
        return result;
}

#else
int main(){
        return 0;
}
#endif
