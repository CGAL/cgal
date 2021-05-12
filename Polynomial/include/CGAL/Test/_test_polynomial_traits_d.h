// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial

#ifndef CGAL_TEST_POLYNOMIAL
#define CGAL_TEST_POLYNOMIAL



#include <iostream>
#include <cassert>

#include <CGAL/algorithm.h>
#include <CGAL/use.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/ipower.h>
#include <CGAL/Random.h>
#include <CGAL/Polynomial_traits_d.h>
#include <cmath>


namespace CGAL {


namespace Test_Pol {


static CGAL::Random my_rnd(346); // some seed

#define CGAL_SNAP_CGALi_TRAITS_D(T)                             \
  typedef T PT;                                                 \
  typedef typename PT::Polynomial_d          Polynomial_d;      \
  typedef typename PT::Coefficient_type           Coeff;             \
  typedef typename PT::Innermost_coefficient_type ICoeff;            \
  typedef CGAL::Polynomial_traits_d<Coeff> PTC;                 \
  typedef CGAL::Exponent_vector Exponent_vector;                \
  typedef std::pair< CGAL::Exponent_vector , ICoeff > Monom;\
 CGAL_USE_TYPE(Polynomial_d);\
 CGAL_USE_TYPE(PTC);\
  CGAL_USE_TYPE(Exponent_vector);\
  CGAL_USE_TYPE(Monom)


#define ASSERT_IS_NULL_FUNCTOR(T)                                       \
  CGAL_static_assertion((boost::is_same<T,CGAL::Null_functor >::value))



template <class Polynomial_d_>
Polynomial_d_
generate_sparse_random_polynomial_(int max_degree = 10){
  typedef CGAL::Polynomial_traits_d<Polynomial_d_> PT;
  CGAL_SNAP_CGALi_TRAITS_D(PT);
  typedef typename PT::Construct_polynomial Constructor;

  typedef CGAL::Exponent_vector Exponent_vector;
  typedef std::pair< CGAL::Exponent_vector , ICoeff > Monom;
  typedef std::vector< Monom > Monom_rep;

  int range = 20;
  int number_of_variables = PT::d;
  double md = max_degree+1;
  int number_of_coeffs =
    (CGAL::min)(number_of_variables * (int)std::ceil(std::log(md))+1,4);

  Polynomial_d result(0);
  for(int i = 0; i < number_of_coeffs; i++){
    CGAL::Exponent_vector exps;
    for(int j = 0; j < PT::d; j++){
      exps.push_back(my_rnd.get_int(0,max_degree));
    }
    ICoeff c = ICoeff(my_rnd.get_int(-range,range));
    Monom_rep monom_rep;
    monom_rep.push_back(Monom(exps,c));
    result += Constructor()(monom_rep.begin(), monom_rep.end());
  }

  return result;
}

// generates a random polynomial that contains at least one coefficient
// that effects the outer most variable.
// Several test in this file rely on this fact.
// The other case is usually handled as a special case.
template <class Polynomial_d_>
Polynomial_d_
generate_sparse_random_polynomial(int max_degree = 10){
  Polynomial_d_ p;
  while(p.degree()==0){
    p = generate_sparse_random_polynomial_<Polynomial_d_>(max_degree);
  }
  return p;
}

template <class Polynomial_traits_d> class Construct_test_polynomial {

  typedef Polynomial_traits_d PT;
  typedef typename PT::Polynomial_d          Polynomial_d;
  typedef typename PT::Coefficient_type           Coeff;
  typedef typename PT::Innermost_coefficient_type ICoeff;
  typedef CGAL::Polynomial_traits_d<Coeff> PTC;
  typedef CGAL::Exponent_vector Exponent_vector;
  typedef std::pair< CGAL::Exponent_vector , ICoeff > Monom;

  typedef typename PT::Construct_polynomial Constructor;

public:
  Polynomial_d operator ()(const Coeff a0, const Coeff a1)const {
    std::list<Coeff> coeffs;
    coeffs.push_back(Coeff(a0));
    coeffs.push_back(Coeff(a1));
    return Constructor()(coeffs.begin(),coeffs.end());
  }
  Polynomial_d operator ()(const Coeff a0, const Coeff a1, const Coeff a2)const {
    std::list<Coeff> coeffs;
    coeffs.push_back(Coeff(a0));
    coeffs.push_back(Coeff(a1));
    coeffs.push_back(Coeff(a2));
    return Constructor()(coeffs.begin(),coeffs.end());
  }
  Polynomial_d operator ()(const Coeff a0, const Coeff a1, const Coeff a2,
      const Coeff a3)const {
    std::list<Coeff> coeffs;
    coeffs.push_back(Coeff(a0));
    coeffs.push_back(Coeff(a1));
    coeffs.push_back(Coeff(a2));
    coeffs.push_back(Coeff(a3));
    return Constructor()(coeffs.begin(),coeffs.end());
  }
  Polynomial_d operator ()(const Coeff a0, const Coeff a1, const Coeff a2,
      const Coeff a3, const Coeff a4, const Coeff a5)const {
    std::list<Coeff> coeffs;
    coeffs.push_back(Coeff(a0));
    coeffs.push_back(Coeff(a1));
    coeffs.push_back(Coeff(a2));
    coeffs.push_back(Coeff(a3));
    coeffs.push_back(Coeff(a4));
    coeffs.push_back(Coeff(a5));
    return Constructor()(coeffs.begin(),coeffs.end());
  }
};

template <class Polynomial_traits_d>
void test_construct_polynomial(const Polynomial_traits_d&){
  std::cerr << "start test_construct_polynomial ";
  std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);

  const int d = PT::d;
  {
    typedef typename PT::Construct_polynomial Constructor;
    int array[] = {1,2,3,4};
    Polynomial_d p = Constructor()(array,array+4);
    assert(p==Construct_test_polynomial<Polynomial_traits_d>()(Coeff(1),Coeff(2),Coeff(3),Coeff(4)));
    assert(Constructor()(array,array) == Polynomial_d(0));
  }
  {
    typedef typename PT::Construct_polynomial Constructor;
    Coeff array[] = {Coeff(1),Coeff(2),Coeff(3),Coeff(4)};
    Polynomial_d p = Constructor()(array,array+4);
    assert(p==Construct_test_polynomial<Polynomial_traits_d>()(Coeff(1),Coeff(2),Coeff(3),Coeff(4)));
    assert(Constructor()(array,array) == Polynomial_d(0));
  }
  {
    typedef typename PT::Construct_polynomial Constructor;
    ICoeff array[] = {ICoeff(1),ICoeff(2),ICoeff(3),ICoeff(4)};
    Polynomial_d p = Constructor()(array,array+4);
    assert(p==Construct_test_polynomial<Polynomial_traits_d>()(Coeff(1),Coeff(2),Coeff(3),Coeff(4)));
    assert(Constructor()(array,array) == Polynomial_d(0));
  }
  { // Construct_polynomial
    typedef typename PT::Construct_polynomial Constructor;
    CGAL_static_assertion(
        !(boost::is_same< Constructor , CGAL::Null_functor >::value));
    typedef typename Constructor::result_type result_type;
    CGAL_static_assertion(
        (boost::is_same< result_type , Polynomial_d >::value));
    CGAL_USE_TYPE(result_type);
    typedef typename PT::Shift Shift;
    typedef typename PT::Evaluate Evaluate;
    Shift shift = Shift();
        Evaluate evaluate = Evaluate();

    Polynomial_d empty = Constructor()();
    assert(empty == Polynomial_d(0));
    assert(empty != Polynomial_d(1));
    assert(Constructor()(3) == Polynomial_d(3));
    assert(Constructor()(4) != Polynomial_d(3));
    assert(Constructor()(Coeff(2)) != Polynomial_d(1));
    assert(Constructor()(Coeff(2)) == Polynomial_d(2));
    assert(Constructor()(ICoeff(3)) != Polynomial_d(1));
    assert(Constructor()(ICoeff(4)) == Polynomial_d(4));

    // construct via iterator range
    Polynomial_d x = shift(Polynomial_d(1), 1,(d-1));
    Polynomial_d result = 0+1*x+2*x*x+3*x*x*x;

    std::list<Coeff> coeffs;
    assert(Constructor()(coeffs.begin(),coeffs.end()) == Constructor()(0));
    for(int i = 0; i<4;i++){coeffs.push_back(Coeff(i));}
    Polynomial_d compare = Constructor()(coeffs.begin(),coeffs.end());
        assert(compare == result);
        assert(evaluate(compare,Coeff(0))== Coeff(0));
        assert(evaluate(compare,Coeff(1))== Coeff(6));
        assert(evaluate(compare,Coeff(2))== Coeff(34));
        assert(evaluate(compare,Coeff(-1))== Coeff(-2));

    typedef std::list< Monom > Monom_list;
    typedef std::vector< Monom > Monom_vec;

    // construct from InputIterator
    Monom_list monom_list;
    assert(Constructor()(monom_list.begin(),monom_list.end()) ==
        Constructor()(0));
    CGAL::Random rnd(7);
    for(int j = 0; j < 2; j++){
      CGAL::Exponent_vector exps;
      for(int i = 0; i < d; i++){
        exps.push_back(j+i*5);
      }
      monom_list.push_back(Monom(exps,ICoeff(j+1)));
    };

    Monom_vec monom_vec(monom_list.begin(),monom_list.end());
    CGAL::cpp98::random_shuffle(monom_vec. begin(),monom_vec. end());

    Polynomial_d p1 = Constructor()(monom_vec. begin(),
        monom_vec. begin()+((monom_vec. end()- monom_vec. begin())/2));
    Polynomial_d p2 = Constructor()(monom_vec. begin()+
        ((monom_vec. end()- monom_vec. begin())/2), monom_vec. end());
    Polynomial_d p  = Constructor()(monom_vec. begin(), monom_vec. end());
    assert(p == p1+p2);

    assert(Constructor()(monom_vec. begin(),monom_vec. end())
        == Constructor()(monom_vec.rbegin(),monom_vec.rend()));
    // test with boolean flag is_sorted
    assert(Constructor()(monom_vec. begin(),monom_vec. end(),false)
        == Constructor()(monom_vec.rbegin(),monom_vec.rend(),false));
  }
  std::cerr << " ok "<< std::endl;
}

template< class Polynomial_traits_d >
void test_get_coefficient(const Polynomial_traits_d&) {

  std::cerr << "start test_get_coefficient ";
  std::cerr.flush();

  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Get_coefficient Get_coeff;
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;

  const Get_coeff get_coeff = Get_coeff();
  (void) get_coeff;

  Polynomial_d p = Construct_test_polynomial()(Coeff(1), Coeff(2), Coeff(3));
  assert(get_coeff(p, 0) == Coeff(1));
  assert(get_coeff(p, 1) == Coeff(2));
  assert(get_coeff(p, 2) == Coeff(3));
  assert(get_coeff(p, 3) == Coeff(0));

  std::cerr << " ok" << std::endl;
}

template< class Polynomial_traits_d >
void test_get_innermost_coefficient(const Polynomial_traits_d&) {
  std::cerr << "start test_get_innermost_coefficient ";
  std::cerr.flush();

  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Get_innermost_coefficient Get_innermost_coeff;
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;

  const Get_innermost_coeff get_innermost_coeff = Get_innermost_coeff();
  (void) get_innermost_coeff;

  Polynomial_d p = Construct_test_polynomial()(Coeff(1), Coeff(2), Coeff(3));

  Exponent_vector ev;

  for(int i = 0; i < PT::d-1; ++i) {
    ev.push_back(0);
  }

  ev.push_back(0);
  assert(get_innermost_coeff(p, ev) == ICoeff(1));

  ev.pop_back();
  ev.push_back(1);
  assert(get_innermost_coeff(p, ev) == ICoeff(2));

  ev.pop_back();
  ev.push_back(2);
  assert(get_innermost_coeff(p, ev) == ICoeff(3));

  ev.pop_back();
  ev.push_back(3);
  assert(get_innermost_coeff(p, ev) == ICoeff(0));

  std::cerr << " ok" << std::endl;
}

template <class Polynomial_traits_d>
void test_get_monom_representation(const Polynomial_traits_d&){
  std::cerr << "start test_get_monom_representation ";
  std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef std::vector< Monom > Monom_rep;
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Monomial_representation Gmr;
  const Gmr gmr = Gmr();

  {
    Polynomial_d zero= Constructor()(0);
    Monom_rep monom_rep;
    gmr(zero,std::back_inserter(monom_rep));
    assert(monom_rep.size()==1);
    assert(Constructor()(monom_rep.begin(),monom_rep.end()) == zero);
  }
  {
    Polynomial_d x = CGAL::shift(Constructor()(1),1,PT::d-1);
    Polynomial_d p = x*x-1;
    Monom_rep monom_rep;
    gmr(p,std::back_inserter(monom_rep));
    assert(monom_rep.size()==2);
    assert(Constructor()(monom_rep.begin(),monom_rep.end()) == p);
  }



  for (int i = 0; i < 5 ; i++){
    Polynomial_d p,q;
    Monom_rep monom_rep;
    p = generate_sparse_random_polynomial<Polynomial_d>();
    gmr(p,std::back_inserter(monom_rep));
    q = Constructor()(monom_rep.begin(), monom_rep.end());
    assert(q == p);
    CGAL::cpp98::random_shuffle(monom_rep.begin(), monom_rep.end());
    q = Constructor()(monom_rep.begin(), monom_rep.end());
    assert(q == p);
  }
  std::cerr << " ok "<< std::endl;
}




template <class Polynomial_traits_d>
void test_swap(const Polynomial_traits_d&){
  std::cerr << "start test_swap "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  int d = PT::d;
  typedef typename Polynomial_traits_d::Swap Swap;
  const Swap swap = Swap();

  Polynomial_d zero = Constructor()(0);
  assert(swap(zero,0,d-1) == zero);

  for(int i = 0; i < 5; i++){
    int i1 = my_rnd.get_int(0,d);
    int i2 = my_rnd.get_int(0,d);
    int i3 = my_rnd.get_int(0,d);

    Polynomial_d p,q;
    p = generate_sparse_random_polynomial<Polynomial_d>();
    q = swap(p,i1,i2);
    assert(p!=q || i1 == i2);
    q = swap(q,i1,i2);
    assert(p == q);

    if(i1 != i2 && i1 != i3 && i2 != i3){
      q = swap(q,i1,i2);
      assert(p!=q || i1 == i2);
      q = swap(q,i1,i3);
      q = swap(q,i1,i2);
      q = swap(q,i2,i3);
      assert(p == q);
    }
  }
  for(int i = 0; i < 5; i++){
    int n = my_rnd.get_int(0,d);
    int m = my_rnd.get_int(0,d);
    Polynomial_d p,q;
    p = generate_sparse_random_polynomial<Polynomial_d>();
    q = generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d pq_1, pq_2;
    pq_1= p*q;
    p = swap(p,n,m);
    q = swap(q,n,m);
    pq_2 = swap(p*q,n,m);
    assert(pq_1 == pq_2);
  }
  std::cerr << " ok "<< std::endl;

}



template <class Polynomial_traits_d>
void test_permute(const Polynomial_traits_d&){
  std::cerr << "start test_permutate "; std::cerr.flush();
  const int d=4;
  typedef typename Polynomial_traits_d::Innermost_coefficient_type ICoeff;
  typedef typename Polynomial_traits_d:: template Rebind<ICoeff,d>::Other PT4;
  typedef typename PT4::Construct_polynomial Constructor;
  typedef typename PT4::Permute Permutate;
  typedef typename PT4::Shift Shift;
  typedef typename PT4::Polynomial_d Polynomial_d;
  typedef typename std::vector<int>::iterator Iterator;

  const Permutate permutate = Permutate();

  Polynomial_d x  = Shift()(Constructor()(1),1,0);
  Polynomial_d y  = Shift()(Constructor()(1),1,1);
  Polynomial_d z  = Shift()(Constructor()(1),1,2);
  Polynomial_d w1 = Shift()(Constructor()(1),1,3);

  std::vector<int> change(d);
  Iterator changeb = change.begin();
  Iterator changee = change.end();
  int lauf = 0;
  for (Iterator itlauf = changeb; itlauf != changee; ++itlauf){
    change[lauf]= lauf;
    lauf++ ;
  }
  Polynomial_d zero = Constructor()(0);
  Polynomial_d one = Constructor()(1);
  assert(permutate (zero,changeb,changee)==zero);
  assert(permutate (one,changeb,changee)!=zero);

  //1. test //
  Polynomial_d Orig = y*z*z*w1*w1*w1;
  Polynomial_d Orig3 = x*x*y*y + z*z*z*w1*w1*w1*w1*x + x*w1*w1;
  change [0] = 1;change [1] = 3;change [2] = 0;change [3] = 2;
  Polynomial_d Res =  permutate  (Orig,changeb, changee);
  Polynomial_d Res3 = permutate  (Orig3,changeb, changee);
  assert(Res != Orig);
  assert(Res3!= Orig3);
  change [0] = 2; change [1] = 0; change [2] = 3; change [3] = 1;
  assert(permutate (Res,changeb,changee)==Orig);
  assert(permutate (Res3,changeb,changee)==Orig3);
  //2. test//
  change [0] = 1; change [1] = 0; change [2] = 3; change [3] = 2;
  Res =  permutate  (Orig,changeb, changee);
  Res3 =  permutate  (Orig3,changeb, changee);
  change [0] = 0; change [1] = 3; change [2] = 2; change [3] = 1;
  Res =  permutate  (Res,changeb, changee);
  Res3 =  permutate  (Res3,changeb, changee);
  change [0] = 3; change [1] = 0; change [2] = 1; change [3] = 2;
  assert(permutate (Orig,changeb,changee)==Res);
  assert(permutate (Orig3,changeb,changee)==Res3);
  //3. test//
  Polynomial_d Should = x*y*y*y*z*z;
  Polynomial_d Should3 = w1*w1*x*x+z*z*z*y*y*y*y*w1+w1*y*y;
  change [0] = 3; change [1] = 0; change [2] = 2; change [3] = 1;
  assert(permutate (Orig,changeb,changee)==Should);
  assert(permutate (Orig3,changeb,changee)==Should3);

  std::cerr << " ok "<< std::endl;
}



template <class Polynomial_traits_d>
void test_move(const Polynomial_traits_d&){
  std::cerr << "start test_move "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename Polynomial_traits_d::Move Move;
  const Move move = Move();(void) move;
  typename Polynomial_traits_d::Swap swap; (void) swap;

  //std::cout << "start_test ----------- "<< d << std::endl;
  for(int i = 0; i < 5; i++){
    int n = (0 == PT::d-1) ? 0 : my_rnd.get_int(0,PT::d-1);  // as [0, 0) is not well defined
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
    if(n <= PT::d-2){
      assert(move(p,n,n+1) == swap(p,n,n+1));
    }
    if(n <= PT::d-3){
      assert(move(p,n,n+2) == swap(swap(p,n,n+1),n+1,n+2));
    }
    if(n >= 1){
      assert(move(p,n-1,n) == swap(p,n-1,n));
    }
    if(n >= 2){
      assert(move(p,n-2,n) == swap(swap(p,n-1,n),n-1,n-2));
    }
  }

  std::cerr << " ok "<< std::endl;

}

template <class Polynomial_traits_d>
void test_degree(const Polynomial_traits_d&){
  std::cerr << "start test_degree "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);

  typedef typename PT::Construct_polynomial Constructor;
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Degree Degree;
  const Degree degree = Degree();
  (void) degree;

  Polynomial_d p;
  p= Constructor()(Coeff(0));
  assert(degree(p) == 0);
  assert(degree(p,0) == 0);
  p= Constructor()(Coeff(1));
  assert(degree(p) == 0);
  assert(degree(p,0) == 0);
  p= Construct_test_polynomial()(Coeff(1),Coeff(2));
  assert(degree(p) == 1);
  p= Constructor()(Coeff(0));
  assert(degree(p,(PT::d-1)) == 0);
  p= Constructor()(Coeff(1));
  assert(degree(p,(PT::d-1)) == 0);
  p= Construct_test_polynomial()(Coeff(1),Coeff(2));
  assert(degree(p,(PT::d-1)) == 1);
  std::cerr << " ok "<< std::endl;
}


//       Total_degree;
template <class Polynomial_traits_d>
void test_total_degree(const Polynomial_traits_d&){

  std::cerr << "start test_total_degree "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Total_degree Total_degree;
  const Total_degree total_degree = Total_degree();
  for(int i = 0; i < 5; i++){
    Polynomial_d p =generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q =generate_sparse_random_polynomial<Polynomial_d>();
    int tdp = total_degree(p); (void) tdp;
    int tdq = total_degree(q); (void) tdq;
    assert(total_degree(p*q) == tdp+tdq);
  }
  std::cerr << " ok "<< std::endl;
}
// //       Leading_coefficient;
template <class Polynomial_traits_d>
void test_leading_coefficient(const Polynomial_traits_d&){

  std::cerr << "start test_leading_coefficient "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Leading_coefficient Lcoeff;
  const Lcoeff lcoeff =Lcoeff();
  for(int i = 0; i < 5; i++){
    Polynomial_d p =generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q =generate_sparse_random_polynomial<Polynomial_d>();
    Coeff lcoeffp = lcoeff(p);
    Coeff lcoeffq = lcoeff(q);
    assert(lcoeff(p*q) == lcoeffp*lcoeffq);
  }
  std::cerr << " ok "<< std::endl;
}

template< class Polynomial_traits_d >
void test_innermost_leading_coefficient(const Polynomial_traits_d&) {

  std::cerr << "start test_innermost_leading_coefficient "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);

  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Innermost_leading_coefficient Ilcoeff;
  const Ilcoeff ilcoeff = Ilcoeff(); (void) ilcoeff;

  Polynomial_d p= Construct_test_polynomial()(Coeff(1), Coeff(2), Coeff(3));
  assert(ilcoeff(p) == ICoeff(3));

  p = generate_sparse_random_polynomial<Polynomial_d>();
  typename PT::Degree_vector degree_vector;  (void) degree_vector;
  typename PT::Get_innermost_coefficient icoeff;  (void) icoeff;
  assert(ilcoeff(p) == icoeff(p,degree_vector(p)));

  std::cerr << " ok" << std::endl;
}

// //       Univariate_content;
template <class Polynomial_traits_d>
void test_univariate_content(const Polynomial_traits_d&){

  std::cerr << "start test_univariate_content "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Univariate_content Univariate_content;
  const Univariate_content univariate_content =Univariate_content();
  assert(univariate_content(Constructor()(0)) == Coeff(0));
  assert(univariate_content(Constructor()(1)) == Coeff(1));
  assert(
      univariate_content(Constructor()(2))
      ==
      CGAL::integral_division(Coeff(2), CGAL::unit_part(Coeff(2))));
  for(int i = 0; i < 5; i++){
    Polynomial_d p =generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q =generate_sparse_random_polynomial<Polynomial_d>();
    Coeff ucontentp = univariate_content(p);
    Coeff ucontentq = univariate_content(q);
    assert(univariate_content(p*q) == ucontentp*ucontentq);
  }

  std::cerr << " ok "<< std::endl;
}
// //       Multivariate_content;
template <class Polynomial_traits_d>
void test_multivariate_content(const Polynomial_traits_d&){
  std::cerr << "start test_multivariate_content ";
  std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Multivariate_content Mcontent;
  const Mcontent mcontent = Mcontent();

  assert(mcontent(Constructor()(0)) ==
      CGAL::integral_division(ICoeff(0),CGAL::unit_part(ICoeff(0))));
  assert(mcontent(Constructor()(1)) ==
      CGAL::integral_division(ICoeff(1), CGAL::unit_part(ICoeff(1))));
  assert(mcontent(Constructor()(2)) ==
      CGAL::integral_division(ICoeff(2), CGAL::unit_part(ICoeff(2))));
  assert(mcontent(Constructor()(-2)) ==
      CGAL::integral_division(ICoeff(2), CGAL::unit_part(ICoeff(2))));

  for(int i = 0; i < 5; i++){
    Polynomial_d p =generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q =generate_sparse_random_polynomial<Polynomial_d>();

    ICoeff content_p = mcontent(p);
    ICoeff content_q = mcontent(q);

    assert(mcontent(p*q) == content_p*content_q);
    assert(mcontent(p*Polynomial_d(2))
        == content_p*
        CGAL::integral_division(ICoeff(2), CGAL::unit_part(ICoeff(2))));
    p = CGAL::integral_division(p,content_p);
    assert(mcontent(p) ==
        CGAL::integral_division(ICoeff(1), CGAL::unit_part(ICoeff(1))));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Shift;
template <class Polynomial_traits_d>
void test_shift(const Polynomial_traits_d&){

  std::cerr << "start test_shift "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Shift Shift;
  const Shift shift= Shift(); (void) shift;
  typename PT::Swap  swap;  (void) swap;
  for(int i = 0; i < 5; i++){
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q = p*CGAL::ipower(Construct_test_polynomial()(Coeff(0),
            Coeff(1)),5);
    assert(shift(p,5) == q);
  }
  int d = Polynomial_traits_d::d;
  if( d > 1){
    assert(shift(Constructor()(1),1,0) !=  shift(Constructor()(1),1,1));
    assert(shift(Constructor()(1),1,0) !=  shift(Constructor()(1),1,1));
    assert(shift(Constructor()(1),1,0) !=  shift(Constructor()(1),1,d-1));
    assert(shift(Constructor()(1),1,0) !=  shift(Constructor()(1),1,d-1));

    assert(shift(Constructor()(1),1,0) == swap(shift(Constructor()(1),1,1),0,1));
    assert(shift(Constructor()(1),1,0) == swap(shift(Constructor()(1),1,1),0,1));

    assert(shift(Constructor()(1),1,0) ==
        swap(shift(Constructor()(1),1,d-1),0,d-1));
    assert(shift(Constructor()(1),1,0) ==
        swap(shift(Constructor()(1),1,d-1),0,d-1));

  }
  std::cerr << " ok "<< std::endl;
}
// //       Negate;
template <class Polynomial_traits_d>
void test_negate(const Polynomial_traits_d&){
  std::cerr << "start test_negate "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Negate Negate;
  typename PT::Swap swap;

  const Negate negate = Negate();
  assert(negate(Constructor()(0)) == Constructor()(0));
  assert(negate(Constructor()(2)) == Constructor()(2));

  Polynomial_d p = Construct_test_polynomial()(Coeff(1),Coeff(2),Coeff(3));
  assert(negate(p) == Construct_test_polynomial()(Coeff(1),Coeff(-2),Coeff(3)));;

  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q = p;
    int n = (0 == PT::d-1) ? 0 : my_rnd.get_int(0,PT::d-1);  // as [0, 0) is not well defined
    int m = (0 == PT::d-1) ? 0 : my_rnd.get_int(0,PT::d-1);  // as [0, 0) is not well defined
    (void) n;  (void) m; (void) negate; (void) swap;
    assert(negate(swap(p,n,m),n) == swap(negate(p,m),n,m));
  }
  std::cerr << " ok "<< std::endl;
}
//      Invert;
template <class Polynomial_traits_d>
void test_invert(const Polynomial_traits_d&){
  std::cerr << "start test_invert "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Invert Invert;
  typename PT::Swap swap;
  const Invert invert = Invert(); (void) invert;
  (void) swap;
  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
    std::vector<Coeff> coeffs (p.begin(),p.end());
    p = invert(p);
    std::vector<Coeff> rcoeffs (p.begin(),p.end());
    assert(coeffs.size() >= rcoeffs.size());
    for (unsigned int i = 0; i < rcoeffs.size(); i++){
      assert(rcoeffs[i] == coeffs[coeffs.size()-i-1]);
    }
    int n = (0 == PT::d-1) ? 0 : my_rnd.get_int(0,PT::d-1);  // as [0, 0) is not well defined

    assert(invert(p,n) == swap(invert(swap(p,n,PT::d-1)),n,PT::d-1));
  }
  std::cerr << " ok "<< std::endl;
}
// //       Translate;
template <class Polynomial_traits_d>
void test_translate(const Polynomial_traits_d&){
  std::cerr << "start test_translate "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Translate Translate;
  typename PT::Evaluate evaluate;
  typename PT::Move move;
  const Translate translate = Translate(); (void) translate;
  (void) evaluate;
  (void) move;

  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
    assert(evaluate(translate(p,ICoeff(5)),ICoeff(3))
        == evaluate(p,ICoeff(8)));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Translate_homogeneous;
template <class Polynomial_traits_d>
void test_translate_homongenous(const Polynomial_traits_d&){
  std::cerr << "start test_translate_homongenous "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Translate_homogeneous Transh;
  typename PT::Canonicalize canonicalize;
  typename PT::Evaluate_homogeneous evh;
  const Transh transh = Transh(); (void) transh;
  (void) canonicalize;
  (void) evh;
  //typename PT::Move move;
  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p,q1,q2;
    p = generate_sparse_random_polynomial<Polynomial_d>();
    q1 = transh(transh(p,ICoeff(5),ICoeff(3)),ICoeff(3),ICoeff(2)) ;
    q2 = transh(p,ICoeff(19),ICoeff(6)) ;
    assert(canonicalize(q1) != canonicalize(p)) ;
    assert(canonicalize(q2) != canonicalize(p)) ;
    assert(canonicalize(q1) == canonicalize(q2));

    assert(
        evh(p,ICoeff(19),ICoeff(6)) == evh(q2,ICoeff(0),ICoeff(1)));

    q1 = transh(transh(p,ICoeff(5),ICoeff(3),0),ICoeff(3),ICoeff(2),0) ;
    q2 = transh(p,ICoeff(19),ICoeff(6),0) ;

    assert(canonicalize(q1) != canonicalize(p)) ;
    assert(canonicalize(q2) != canonicalize(p)) ;
    assert(canonicalize(q1) == canonicalize(q2));
  }
  std::cerr << " ok "<< std::endl;
}

template< class Polynomial_traits_d>
void test_scale(const Polynomial_traits_d&) {
  (std::cerr << "start test_scale ").flush();

  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;

  typedef typename PT::Scale Scale;
  const Scale scale = Scale();(void) scale;
  Polynomial_d p = Construct_test_polynomial()(Coeff(1), Coeff(2), Coeff(3));

  assert(
      scale(p, ICoeff(2)) == Construct_test_polynomial()(Coeff(1), Coeff(4),
          Coeff(12)));

  std::cerr << " ok" << std::endl;
}

// //       Scale_homogeneous;
template <class Polynomial_traits_d>
void test_scale_homogeneous(const Polynomial_traits_d&){
  std::cerr << "start test_scale_homogeneous "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Scale_homogeneous Scaleh;
  typename PT::Canonicalize canonicalize;
  const Scaleh scaleh= Scaleh(); (void) scaleh;
  (void) canonicalize;
  //typename PT::Move move;
  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p,q1,q2;
    p = generate_sparse_random_polynomial<Polynomial_d>();
    q1 = scaleh(scaleh(p,ICoeff(5),ICoeff(3)),ICoeff(3),ICoeff(2)) ;
    q2 = scaleh(p,ICoeff(15),ICoeff(6)) ;
    assert(canonicalize(q1) != canonicalize(p)) ;
    assert(canonicalize(q2) != canonicalize(p)) ;
    assert(canonicalize(q1) == canonicalize(q2));
      q1 = scaleh(scaleh(p,ICoeff(5),ICoeff(3),0),ICoeff(3),ICoeff(2),0) ;
      q2 = scaleh(p,ICoeff(15),ICoeff(6),0) ;

      assert(canonicalize(q1) != canonicalize(p)) ;
      assert(canonicalize(q2) != canonicalize(p)) ;
      assert(canonicalize(q1) == canonicalize(q2));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Differentiate;
template <class Polynomial_traits_d>
void test_differentiate (const Polynomial_traits_d&){
  std::cerr << "start test_differentiate "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Differentiate Diff;
  typename PT::Swap swap;
  const Diff diff= Diff(); (void) diff;
  (void) swap;

  assert(diff(Constructor()(0)) == Constructor()(0));
  assert(diff(Constructor()(1)) == Constructor()(0));
  assert(diff(Construct_test_polynomial()(Coeff(1), Coeff(2))) ==
      Polynomial_d(2));

  for(int i = 0 ; i < 5 ; i++){
    int n = (0 == PT::d-1) ? 0 : my_rnd.get_int(0,PT::d-1);  // as [0, 0) is not well defined
    Polynomial_d p,pd;
    p = generate_sparse_random_polynomial<Polynomial_d>();
    pd = diff(p,n);
    assert(pd == swap(diff(swap(p,n,PT::d-1)),n,PT::d-1));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Is_square_free;
template <class Polynomial_traits_d>
void test_is_square_free(const Polynomial_traits_d&){
  std::cerr << "start test_is_square_free "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Is_square_free Is_square_free;
  typedef typename PT::Make_square_free Make_square_free;


  const Is_square_free is_square_free = Is_square_free();
  const Make_square_free make_square_free = Make_square_free();

  assert(true == is_square_free(Constructor()(0)));
  assert(true == is_square_free(Constructor()(1)));
  assert(true == is_square_free(Constructor()(2)));

  //typename PT::Canonicalize canonicalize;
  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p;
    p = generate_sparse_random_polynomial<Polynomial_d>(3);
    assert(is_square_free(make_square_free(p)));
    assert(!is_square_free(p*p));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Make_square_free;
template <class Polynomial_traits_d>
void test_make_square_free(const Polynomial_traits_d&){
  std::cerr << "start test_make_square_free "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Make_square_free Make_square_free;
  typename PT::Leading_coefficient lcoeff;
  typename PT::Univariate_content_up_to_constant_factor ucontent_utcf;
  typename PT::Integral_division_up_to_constant_factor idiv_utcf;


  const Make_square_free make_square_free = Make_square_free();
  assert(Constructor()(0) == make_square_free(Constructor()(0)));
  assert(Constructor()(1) == make_square_free(Constructor()(1)));
  assert(Constructor()(1) == make_square_free(Constructor()(2)));

  //typename PT::Canonicalize canonicalize;
  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p;
    p = generate_sparse_random_polynomial<Polynomial_d>(3);

    p = idiv_utcf(p, Constructor()(ucontent_utcf(p)));
    p = make_square_free(p);
    Coeff f_intern = lcoeff(p);

    f_intern = typename PTC::Make_square_free()(f_intern);

    //std::cout <<"f: " << f_intern << std::endl;

    assert(p * f_intern == make_square_free(p*p*f_intern));
    assert(p * f_intern == make_square_free(p*f_intern*f_intern));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Square_free_factorize;
template <class Polynomial_traits_d>
void test_square_free_factorize(const Polynomial_traits_d&){
  std::cerr << "start test_square_free_factorize "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef CGAL::Algebraic_structure_traits<Polynomial_d> AST;
  typename AST::Integral_division idiv;
  typedef typename PT::Square_free_factorize Sqff;
  typename PT::Canonicalize canonicalize;
  typename PT::Total_degree total_degree;
  typename PT::Leading_coefficient lcoeff;
  (void) idiv;
  const Sqff sqff = Sqff(); (void) sqff;
  (void) canonicalize;

  for(int i = 0; i < 5; i++){
    Polynomial_d f1 = generate_sparse_random_polynomial<Polynomial_d>(3);
    Polynomial_d f2 = generate_sparse_random_polynomial<Polynomial_d>(3);
    f2 = Constructor()(lcoeff(f2));
    Polynomial_d p = f1*f1*f2*f2*f2;
    if( !CGAL::is_zero(p)){
      std::vector<std::pair<Polynomial_d,int> > fac_mul_pairs;
      sqff(p, std::back_inserter(fac_mul_pairs));
      std::size_t n = fac_mul_pairs.size();
      for (std::size_t j = 0; j < n; j++){
        Polynomial_d factor = fac_mul_pairs[j].first;
        assert( total_degree(factor) > 0);
        int multi = fac_mul_pairs[j].second;
        for (int k = 0; k < multi; k++){
          p = idiv(p,factor);
        }
      }
      assert(CGAL::is_one(canonicalize(p)));
    }
  }

  typename PT::Innermost_leading_coefficient ileading_coeff;
  typename PT::Multivariate_content multivariate_content;
  (void) ileading_coeff;
  (void) multivariate_content;
  Polynomial_d p = generate_sparse_random_polynomial< Polynomial_d >(2);
  p *= p*generate_sparse_random_polynomial< Polynomial_d >(2);
  std::vector< std::pair< Polynomial_d,int> > fac_mul_pairs;

  ICoeff alpha;


  sqff(p, std::back_inserter(fac_mul_pairs), alpha);

  assert(alpha
      == CGAL::unit_part(ileading_coeff(p)) * multivariate_content(p));

  //std::cerr << std::endl;
  Polynomial_d rec_p = Constructor()(alpha);
  for( unsigned int i = 0; i < fac_mul_pairs.size() ; i ++){
    Polynomial_d factor = fac_mul_pairs[i].first;
    int mult = fac_mul_pairs[i].second;
    assert(factor == canonicalize(factor));
    rec_p *= CGAL::ipower(factor,mult);
  }
  assert(p == rec_p);

  std::cerr << "ok"<< std::endl;

}
// //       Pseudo_division;
template <class Polynomial_traits_d>
void test_pseudo_division(const Polynomial_traits_d&){
  std::cerr << "start test_pseudo_division "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Pseudo_division Pdiv;
  const Pdiv pdiv = Pdiv();
  for(int i = 0; i < 10; i++){
    Polynomial_d f = generate_sparse_random_polynomial<Polynomial_d>(3);
    Polynomial_d g = generate_sparse_random_polynomial<Polynomial_d>(2);
    Coeff D;
    Polynomial_d q,r;
    if (!CGAL::is_zero(q)){
      pdiv(f,g,q,r,D);
      assert(f*Constructor()(D) == g*q+r);
    }
  }
  std::cerr << " ok "<< std::endl;
}

// //       Pseudo_division_remainder;
template <class Polynomial_traits_d>
void test_pseudo_division_remainder(const Polynomial_traits_d&){
  std::cerr << "start test_pseudo_division_remainder "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Pseudo_division pdiv;
  typedef typename PT::Pseudo_division_remainder Pdiv_r;
  (void) pdiv;
  const Pdiv_r pdiv_r = Pdiv_r(); (void) pdiv_r;

  for(int i = 0; i < 10; i++){
    Polynomial_d f = generate_sparse_random_polynomial<Polynomial_d>(3);
    Polynomial_d g = generate_sparse_random_polynomial<Polynomial_d>(2);
    Coeff D;
    Polynomial_d q,r;
    if (!CGAL::is_zero(q)){
      pdiv(f,g,q,r,D);
      assert(r == pdiv_r(f,g));
    }
  }
  std::cerr << " ok "<< std::endl;
}

// //       Pseudo_division_quotient;
template <class Polynomial_traits_d>
void test_pseudo_division_quotient(const Polynomial_traits_d&){
  std::cerr << "start test_pseudo_division_quotient "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Pseudo_division pdiv;
  typedef typename PT::Pseudo_division_quotient Pdiv_q;
  (void) pdiv;
  const Pdiv_q pdiv_q = Pdiv_q(); (void) pdiv_q;
  for(int i = 0; i < 10; i++){
    Polynomial_d f = generate_sparse_random_polynomial<Polynomial_d>(3);
    Polynomial_d g = generate_sparse_random_polynomial<Polynomial_d>(2);
    Coeff D;
    Polynomial_d q,r;
    if (!CGAL::is_zero(q)){
      pdiv(f,g,q,r,D);
      assert(q == pdiv_q(f,g));
    }
  }
  std::cerr << " ok "<< std::endl;
}

// //       Gcd_up_to_constant_factor;
template <class Polynomial_traits_d>
void test_gcd_up_to_constant_factor(const Polynomial_traits_d&){
  std::cerr << "start test_gcd_up_to_constant_factor "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Gcd_up_to_constant_factor Gcd_utcf;
  const Gcd_utcf gcd_utcf = Gcd_utcf(); (void) gcd_utcf;

  assert(
      Constructor()(0) == gcd_utcf(Constructor()(0),Constructor()(0)));
  assert(
      Constructor()(1) == gcd_utcf(Constructor()(1),Constructor()(0)));
  assert(
      Constructor()(1) == gcd_utcf(Constructor()(0),Constructor()(1)));
  assert(
      Constructor()(1) == gcd_utcf(Constructor()(1),Constructor()(1)));
  assert(
      Constructor()(1) == gcd_utcf(Constructor()(-1),Constructor()(-1)));
  assert(
      Constructor()(1) == gcd_utcf(Constructor()(2),Constructor()(2)));

  std::cerr << " ok "<< std::endl;
}

// //       Integral_division_up_to_constant_factor;
template <class Polynomial_traits_d>
void test_integral_division_up_to_constant_factor(const Polynomial_traits_d&){
  std::cerr << "start test_integral_division_up_to_constant_factor ";
  std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Integral_division_up_to_constant_factor Idiv_utcf;
  typename PT::Canonicalize canonicalize;
  const Idiv_utcf idiv_utcf = Idiv_utcf(); (void) idiv_utcf;
  (void) canonicalize;
  assert(
      Constructor()(0) == idiv_utcf(Constructor()(0),Constructor()(1)));
  assert(
      Constructor()(1) == idiv_utcf(Constructor()(1),Constructor()(1)));
  assert(
      Constructor()(1) == idiv_utcf(Constructor()(2),Constructor()(1)));

  for(int i = 0; i < 5; i++){
    Polynomial_d p,q;
    p = generate_sparse_random_polynomial<Polynomial_d>(3);
    q = generate_sparse_random_polynomial<Polynomial_d>(2);
    assert(canonicalize(p) == idiv_utcf(p*q,q));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Content_up_to_constant_factor;
template <class Polynomial_traits_d>
void test_univariate_content_up_to_constant_factor(const Polynomial_traits_d&){
  std::cerr << "start test_univariate_content_up_to_constant_factor ";
  std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Univariate_content_up_to_constant_factor Ucontent_utcf;
  typename PT::Integral_division_up_to_constant_factor idiv_utcf;
  typename PT::Leading_coefficient lcoeff;
  typename PT::Canonicalize canonicalize;

  const Ucontent_utcf ucontent_utcf = Ucontent_utcf(); (void) ucontent_utcf;
  (void) idiv_utcf;
  (void) lcoeff;
  (void) canonicalize;

  assert(Coeff(0) == ucontent_utcf(Constructor()(0)));
  assert(Coeff(1) == ucontent_utcf(Constructor()(1)));
  assert(Coeff(1) == ucontent_utcf(Constructor()(2)));
  assert(Coeff(1) == ucontent_utcf(Constructor()(-2)));

  for(int i = 0; i < 5; i++){
    Polynomial_d p,q;
    p = generate_sparse_random_polynomial<Polynomial_d>(3);
    Coeff content = ucontent_utcf(p);
    p = idiv_utcf(p,Constructor()(content));
    assert(Coeff(1) == ucontent_utcf(p));
    Coeff lc = lcoeff(p);
    p = p*lc*lc;
    assert(canonicalize(Constructor()(lc*lc)) == ucontent_utcf(p));

  }

  std::cerr << " ok "<< std::endl;
}


// //       Square_free_factorize_up_to_constant_factor;
template <class Polynomial_traits_d>
void test_square_free_factorize_up_to_constant_factor(const Polynomial_traits_d&)
{
  std::cerr << "start test_square_free_factorize_up_to_constant_factor ";
  std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Integral_division_up_to_constant_factor idiv_utcf;
  typedef typename PT::Square_free_factorize_up_to_constant_factor Sqff_utcf;
  typename PT::Canonicalize canonicalize;

  (void) idiv_utcf;
  const Sqff_utcf sqff_utcf = Sqff_utcf(); (void) sqff_utcf;
  (void) canonicalize;

  for(int i = 0; i < 5; i++){
    Polynomial_d f1 = generate_sparse_random_polynomial<Polynomial_d>(2);
    Polynomial_d f2 = generate_sparse_random_polynomial<Polynomial_d>(2);
    Polynomial_d p = f1*f1*f2;
    std::vector<std::pair<Polynomial_d,int> > fac_mul_pairs;
    sqff_utcf(p, std::back_inserter(fac_mul_pairs)) ;

    std::size_t n = fac_mul_pairs.size();

    for (std::size_t j = 0; j < n; j++){
      Polynomial_d factor = fac_mul_pairs[j].first;
      assert(factor == canonicalize(factor));
      int multi = fac_mul_pairs[j].second;
      for (int k = 0; k < multi; k++){
        p = idiv_utcf(p,factor);
      }
    }
    assert(CGAL::is_one(canonicalize(p)));
  }

  std::cerr << "ok"<< std::endl;
}





// //       Evaluate;
template <class Polynomial_traits_d>
void test_evaluate(const Polynomial_traits_d&){
  std::cerr << "start test_evaluate "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Evaluate Evaluate;
  typename PT::Move move;
  const Evaluate evaluate = Evaluate(); (void) evaluate;
  (void) move;
  assert(evaluate(Constructor()(0),Coeff(0)) == Coeff(0));
  assert(evaluate(Constructor()(1),Coeff(0)) == Coeff(1));
  assert(evaluate(Constructor()(2),Coeff(5)) == Coeff(2));

  assert( evaluate(Construct_test_polynomial()(Coeff(3),Coeff(2)),Coeff(0)) ==
      Coeff(3));
  assert( evaluate(Construct_test_polynomial()(Coeff(3),Coeff(2)),Coeff(1)) ==
      Coeff(5));
  assert( evaluate(Construct_test_polynomial()(Coeff(3),Coeff(2)),Coeff(2)) ==
      Coeff(7));

  std::cerr << " ok "<< std::endl;
}

// //       Evaluate_homogeneous;
template <class Polynomial_traits_d>
void test_evaluate_homogeneous(const Polynomial_traits_d&){
  std::cerr << "start test_evaluate_homogeneous "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Evaluate_homogeneous Evh;
  const Evh evh = Evh(); (void) evh;

  assert(evh(Constructor()(0),Coeff(0),Coeff(1)) == Coeff(0));
  assert(evh(Constructor()(1),Coeff(0),Coeff(2)) == Coeff(1));
  assert(evh(Constructor()(2),Coeff(5),Coeff(3)) == Coeff(2));

  assert(evh(Construct_test_polynomial()(Coeff(3),Coeff(2)) , Coeff(0),Coeff(1))
      == Coeff(3));
  assert(evh(Construct_test_polynomial()(Coeff(3),Coeff(2)) , Coeff(1),Coeff(1))
      == Coeff(5));
  assert(evh(Construct_test_polynomial()(Coeff(3),Coeff(2)) , Coeff(2),Coeff(3))
      == Coeff(9+4));

  std::cerr << " ok "<< std::endl;
}

template< class Polynomial_traits_d >
void test_is_zero_at(const Polynomial_traits_d&) {
  std::cerr << "start test_is_zero_at ";
  std::cerr.flush();

  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Is_zero_at Is_zero_at;
  const Is_zero_at is_zero_at = Is_zero_at();(void) is_zero_at;

  Polynomial_d p = Construct_test_polynomial()(Coeff(-1), Coeff(0), Coeff(1));

  std::list< ICoeff > cv;
  for(int i = 0; i < PT::d-1; ++i)
    cv.push_back(ICoeff(0));

  cv.push_back(ICoeff(0));
  assert(!is_zero_at(p, cv.begin(), cv.end()));

  cv.pop_back();
  cv.push_back(ICoeff(1));
  assert(is_zero_at(p, cv.begin(), cv.end()));

  cv.pop_back();
  cv.push_back(ICoeff(-1));
  assert(is_zero_at(p, cv.begin(), cv.end()));

  std::cerr << " ok" << std::endl;
}

template< class Polynomial_traits_d >
void test_is_zero_at_homogeneous(const Polynomial_traits_d&) {
  std::cerr << "start test_is_zero_at_homogeneous ";
  std::cerr.flush();

  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Is_zero_at_homogeneous Is_zero_at_homogeneous;
  const Is_zero_at_homogeneous is_zero_at_homogeneous = Is_zero_at_homogeneous();
  (void) is_zero_at_homogeneous;

  Polynomial_d p = Construct_test_polynomial()(Coeff(-1), Coeff(0), Coeff(4));

  std::list< ICoeff > cv;
  for(int i = 0; i < PT::d-1; ++i)
    cv.push_back(ICoeff(0));

  for(int v = 2; v < 5; ++v) {
    cv.push_back(ICoeff(v));
    cv.push_back(ICoeff(2*v));
    assert(is_zero_at_homogeneous(p, cv.begin(), cv.end()));
    cv.pop_back();
    cv.pop_back();

    cv.push_back(ICoeff(v));
    cv.push_back(ICoeff(2*-v));
    assert(is_zero_at_homogeneous(p, cv.begin(), cv.end()));
    cv.pop_back();
    cv.pop_back();
  }

  std::cerr << " ok" << std::endl;
}

template< class Polynomial_traits_d >
void test_sign_at(const Polynomial_traits_d&) {
  std::cerr << "start test_sign_at ";
  std::cerr.flush();

  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Sign_at Sign_at;
  const Sign_at sign_at = Sign_at();(void) sign_at;

  Polynomial_d p= Construct_test_polynomial()(Coeff(-1), Coeff(0), Coeff(1));

  std::list< ICoeff > cv;
  for(int i = 0; i < PT::d-1; ++i)
    cv.push_back(ICoeff(0));

  cv.push_back(ICoeff(0));
  assert(sign_at(p, cv.begin(), cv.end()) == CGAL::NEGATIVE);

  cv.pop_back();
  cv.push_back(ICoeff(1));
  assert(sign_at(p, cv.begin(), cv.end()) == CGAL::ZERO);

  cv.pop_back();
  cv.push_back(ICoeff(-1));
  assert(sign_at(p, cv.begin(), cv.end()) == CGAL::ZERO);

  cv.pop_back();
  cv.push_back(ICoeff(2));
  assert(sign_at(p, cv.begin(), cv.end()) == CGAL::POSITIVE);

  std::cerr << " ok" << std::endl;
}

template< class Polynomial_traits_d >
void test_sign_at_homogeneous(const Polynomial_traits_d&) {
  std::cerr << "start test_sign_at_homogeneous ";
  std::cerr.flush();

  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename PT::Sign_at_homogeneous Sign_at_homogeneous;
  const Sign_at_homogeneous sign_at_homogeneous = Sign_at_homogeneous();
  (void) sign_at_homogeneous;

  Polynomial_d p= Construct_test_polynomial()(Coeff(-1), Coeff(0), Coeff(1));

  std::list< ICoeff > cv;
  for(int i = 0; i < PT::d-1; ++i)
    cv.push_back(ICoeff(0));

  for(int v = 1; v < 5; ++v) {
    cv.push_back(ICoeff(0));
    cv.push_back(ICoeff(v));
    assert(
        sign_at_homogeneous(p, cv.begin(), cv.end()) == CGAL::NEGATIVE);

    cv.pop_back();
    cv.pop_back();
    cv.push_back(ICoeff(v));
    cv.push_back(ICoeff(v));
    assert(
        sign_at_homogeneous(p, cv.begin(), cv.end()) == CGAL::ZERO);

    cv.pop_back();
    cv.pop_back();
    cv.push_back(ICoeff(-v));
    cv.push_back(ICoeff(v));
    assert(
        sign_at_homogeneous(p, cv.begin(), cv.end()) == CGAL::ZERO);

    cv.pop_back();
    cv.pop_back();
    cv.push_back(ICoeff(v+1));
    cv.push_back(ICoeff(v));
    assert(
        sign_at_homogeneous(p, cv.begin(), cv.end()) == CGAL::POSITIVE);

    cv.pop_back();
    cv.pop_back();
  }

  std::cerr << " ok" << std::endl;
}

template< class Polynomial_traits_d>
void test_compare(const Polynomial_traits_d&) {
  (std::cerr << "start test compare ").flush();

  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);

  typedef typename PT::Compare Compare;
  typedef typename PT::Construct_polynomial Constructor;
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  const Compare compare = Compare(); (void) compare;

  Polynomial_d p0 = Constructor()(Coeff(0));
  Polynomial_d pp2 = Constructor()(Coeff(2));
  Polynomial_d pm2 = Constructor()(Coeff(-2));
  Polynomial_d pp1p2 = Construct_test_polynomial()(Coeff(1), Coeff(2));
  Polynomial_d pm1m2 = Construct_test_polynomial()(Coeff(-1), Coeff(-2));

  assert(compare(p0, p0) == CGAL::EQUAL);
  assert(compare(pp2, pp2) == CGAL::EQUAL);
  assert(compare(pm2, pm2) == CGAL::EQUAL);
  assert(compare(pp1p2, pp1p2) == CGAL::EQUAL);
  assert(compare(pm1m2, pm1m2) == CGAL::EQUAL);

  assert(compare(p0, pp2) == CGAL::SMALLER);
  assert(compare(p0, pm2) == CGAL::LARGER);
  assert(compare(pp2, pm2) == CGAL::LARGER);
  assert(compare(pm1m2, pp1p2) == CGAL::SMALLER);

  assert(compare(pp1p2, pp2) == CGAL::LARGER);
  assert(compare(pm1m2, pp2) == CGAL::SMALLER);
  assert(compare(pp1p2, pm2) == CGAL::LARGER);
  assert(compare(pm1m2, pm2) == CGAL::SMALLER);

  std::cerr << " ok" << std::endl;
}


// //       Resultant;
template <class Polynomial_traits_d>
void test_resultant(const Polynomial_traits_d&){
  std::cerr << "start test_resultant "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef Construct_test_polynomial<PT> Construct_test_polynomial;
  typedef typename  PT::Resultant Resultant;
  typename  PT::Move move;
  const Resultant resultant = Resultant();(void) resultant;
  (void) move;
  {
    Polynomial_d A = Constructor()(0);
    Polynomial_d B = Constructor()(0);
    assert(resultant(A,B) == Coeff(0));
  }{
    Polynomial_d A = Constructor()(4);
    Polynomial_d B = Constructor()(8);
    assert(resultant(A,B) == Coeff(1));
  }{
    Polynomial_d f = Construct_test_polynomial()(Coeff(2),Coeff(7),Coeff(1),
        Coeff(8),Coeff(1),Coeff(8));
    Polynomial_d g = Construct_test_polynomial()(Coeff(3),Coeff(1),Coeff(4),
        Coeff(1),Coeff(5),Coeff(9));

    assert(resultant(f,g) == Coeff(230664271L)); // Maple

    Polynomial_d h = Construct_test_polynomial()(Coeff(3),Coeff(4),Coeff(7),
        Coeff(7));
    Polynomial_d fh = f*h;
    Polynomial_d gh = g*h;
    assert(resultant(fh,gh) == Coeff(0));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Canonicalize;
template <class Polynomial_traits_d>
void test_canonicalize(const Polynomial_traits_d&){
  std::cerr << "start test_canonicalize "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Canonicalize Canonicalize;
  const Canonicalize canonicalize = Canonicalize();(void) canonicalize;

  assert(Constructor()(0) == canonicalize(Constructor()(0)));
  assert(Constructor()(1) == canonicalize(Constructor()(1)));
  assert(Constructor()(1) == canonicalize(Constructor()(2)));
  assert(Constructor()(1) == canonicalize(Constructor()(-2)));

  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>(3);
    Polynomial_d q = generate_sparse_random_polynomial<Polynomial_d>(3);
    assert(canonicalize(p)*canonicalize(q) == canonicalize(p*q));
  }
  std::cerr << " ok "<< std::endl;
}
// //       Substitute;
template <class Polynomial_traits_d>
void test_substitute(const Polynomial_traits_d&){
  std::cerr << "start test_substitute "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Substitute Substitute;
  typedef typename PT::Innermost_coefficient_type Innermost_coefficient_type;

  const Substitute substitute = Substitute(); (void) substitute;
  std::list<Innermost_coefficient_type> list;
  for(int i = 0; i < PT::d; i++){
    list.push_back(Innermost_coefficient_type(i));
  }
  assert(Innermost_coefficient_type(0)
      == substitute(Constructor()(0),list.begin(),list.end()));
  assert(Innermost_coefficient_type(1)
      == substitute(Constructor()(1),list.begin(),list.end()));
  assert(Innermost_coefficient_type(2)
      == substitute(Constructor()(2),list.begin(),list.end()));
  assert(Innermost_coefficient_type(-2)
      == substitute(Constructor()(-2),list.begin(),list.end()));

  for(int i = 0; i< 5; i++){
    typedef typename PT
      :: template Rebind<Innermost_coefficient_type,2>::Other PT_2;
    typedef typename PT_2::Polynomial_d Polynomial_2;
    std::vector<Polynomial_2> vector1,vector2;
    for(int j = 0; j < PT::d; j++){
      vector1.push_back(
          generate_sparse_random_polynomial<Polynomial_2>(3));
    }
    vector2=vector1;
    std::swap(vector2[0],vector2[PT::d-1]);
    Polynomial_d p
      = generate_sparse_random_polynomial<Polynomial_d>(3);
    assert( substitute(p,vector1.begin(),vector1.end()) ==
        substitute(typename PT::Swap()(p,0,PT::d-1),
            vector2.begin(),vector2.end()));
  }

  std::cerr << " ok "<< std::endl;
}

// //       Substitute;
template <class Polynomial_traits_d>
void test_substitute_homogeneous(const Polynomial_traits_d&){
  std::cerr << "start test_substitute_homogeneous "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef typename PT::Construct_polynomial Constructor;
  typedef typename PT::Substitute_homogeneous Substitute_homogeneous;
  const Substitute_homogeneous substitute_homogeneous = Substitute_homogeneous();
  (void) substitute_homogeneous;
  typedef typename PT::Innermost_coefficient_type Innermost_coefficient_type;


  std::list<Innermost_coefficient_type> list;
  for(int i = 0; i < PT::d; i++){
    list.push_back(Innermost_coefficient_type(i));
  }
  list.push_back(Innermost_coefficient_type(3));
  assert(Innermost_coefficient_type(0)
      == substitute_homogeneous(Constructor()(0),list.begin(),list.end()));
  assert(Innermost_coefficient_type(1)
      == substitute_homogeneous(Constructor()(1),list.begin(),list.end()));
  assert(Innermost_coefficient_type(2)
      == substitute_homogeneous(Constructor()(2),list.begin(),list.end()));
  assert(Innermost_coefficient_type(-2)
      == substitute_homogeneous(Constructor()(-2),list.begin(),list.end()));

  for(int i = 0; i< 2; i++){
    typedef typename PT
      :: template Rebind<Innermost_coefficient_type,2>::Other PT_2;
    typedef typename PT_2::Polynomial_d Polynomial_2;
    std::vector<Polynomial_2> vec1,vec2;
    for(int j = 0; j < PT::d+1; j++){
      vec1.push_back(
          generate_sparse_random_polynomial<Polynomial_2>(3));
    }
    vec2=vec1;
    std::swap(vec2[0],vec2[PT::d-1]);
    Polynomial_d p
      = generate_sparse_random_polynomial<Polynomial_d>(3);
    assert( substitute_homogeneous(p,vec1.begin(),vec1.end()) ==
        substitute_homogeneous(typename PT::Swap()(p,0,PT::d-1),
            vec2.begin(),vec2.end()));
  }

  std::cerr << " ok "<< std::endl;
}

template <class PT >
void test_construct_coefficient_const_iterator_range(const PT&) {
  std::cerr << "start test_construct_coefficient_const_iterator_range ";
  std::cerr.flush();

  typedef typename PT::Polynomial_d                Polynomial_d;
  typedef typename PT::Coefficient_const_iterator  CCIterator;

  typedef typename PT::Construct_coefficient_const_iterator_range Coeff_range;
  typename PT::Degree                  degree;
  typename PT::Get_coefficient         coeff;
  const Coeff_range coeff_range = Coeff_range();
  Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();


  CCIterator it = coeff_range(p).first;
  for(int i = 0; i <= degree(p); i++){
    assert(*it == coeff(p,i));
    it++;
  }
  assert(coeff_range(p).second == it);

  std::cerr << " ok "<< std::endl;
}

template <class PT>
void test_construct_innermost_coefficient_const_iterator_range(const PT&) {
  std::cerr << "start test_construct_innermost_coefficient_const_iterator_range ";
  std::cerr.flush();
  typedef typename PT::Innermost_coefficient_type NT;

  typedef typename PT:: template Rebind< NT, 1 >::Other PT_1;
  typedef typename PT:: template Rebind< NT, 2 >::Other PT_2;
  typedef typename PT:: template Rebind< NT, 3 >::Other PT_3;

  typedef typename PT_1::Polynomial_d Polynomial_1;
  typedef typename PT_2::Polynomial_d Polynomial_2;
  typedef typename PT_3::Polynomial_d Polynomial_3;

  Polynomial_1
    p1(NT( 1), NT( 2), NT( 3)),
    p2(NT( 4), NT( 5), NT( 6)),
    p3(NT( 7), NT( 8), NT( 9)),
    p4(NT(10), NT(11), NT(12)),
    p5(NT(13), NT(14), NT(15)),
    p6(NT(16), NT(17), NT(18)),
    p7(NT(19), NT(20), NT(21)),
    p8(NT(22), NT(23), NT(24)),
    p9(NT(25), NT(26), NT(27));

  Polynomial_2
    q1(p1, p2, p3),
    q2(p4, p5, p6),
    q3(p7, p8, p9);

  Polynomial_3 r(q1, q2, q3);

  int i;


  typename PT_1::Innermost_coefficient_const_iterator it1;(void) it1;
  typedef typename PT_1::
    Construct_innermost_coefficient_const_iterator_range Range1;
  const Range1 range1 = Range1();
  typename PT_2::Innermost_coefficient_const_iterator it2;(void) it2;
  typedef typename PT_2::
    Construct_innermost_coefficient_const_iterator_range Range2;
  const Range2 range2 = Range2();
  typename PT_3::Innermost_coefficient_const_iterator it3;(void) it3;
  typedef typename PT_3::
    Construct_innermost_coefficient_const_iterator_range Range3;
  const Range3 range3 = Range3();

  (void) range1;
  (void) range2;
  (void) range3;

  for (i = 1, it1 = (range1(p1).first); i <= 3; ++i, ++it1)
    assert(*it1 == i);
  assert(it1 == range1(p1).second);
  for (i = 1, it2 = range2(q1).first; i <= 9; ++i, ++it2)
    assert(*it2 == i);
  assert(it2 == range2(q1).second);
  for (i = 1, it3 = range3(r).first; i <= 27; ++i, ++it3)
    assert(*it3 == i);
  assert(it3 == range3(r).second);

  std::cerr << " ok "<< std::endl;
}




template< class PT >
void test_fundamental_functors(const PT& traits){
  std::cout << "\n start test for dimension: "
            << PT::d << std::endl;

  // Construction
  test_construct_polynomial(traits);

  // Gets
  test_get_coefficient(traits);
  test_get_innermost_coefficient(traits);
  test_get_monom_representation(traits);
  test_leading_coefficient(traits);
  test_innermost_leading_coefficient(traits);

  test_degree(traits);
  test_total_degree(traits);

  // modifier
  test_swap(traits);
  test_move(traits);

  test_substitute(traits);
  test_substitute_homogeneous(traits);

  test_shift(traits);
  test_negate(traits);
  test_invert(traits);
  test_translate(traits);
  test_translate_homongenous(traits);
  test_scale(traits);
  test_scale_homogeneous(traits);

  test_differentiate(traits);
  test_make_square_free(traits);
  test_is_square_free(traits);
  test_canonicalize(traits);

  // evaluates (sign depends on real embeddable)
  test_evaluate(traits);
  test_evaluate_homogeneous(traits);
  test_is_zero_at(traits);
  test_is_zero_at_homogeneous(traits);

  // pseudo division
  test_pseudo_division(traits);
  test_pseudo_division_remainder(traits);
  test_pseudo_division_quotient(traits);

  // utcf functions
  test_gcd_up_to_constant_factor(traits);
  test_integral_division_up_to_constant_factor(traits);
  test_univariate_content_up_to_constant_factor(traits);
  test_square_free_factorize_up_to_constant_factor(traits);

  // resultant
  test_resultant(traits);

}

template< class PT >
void test_real_embeddable_functors(const PT& traits, CGAL::Tag_true){
  test_sign_at(traits);
  test_sign_at_homogeneous(traits);
  test_compare(traits);
}

template< class PT >
void test_real_embeddable_functors(const PT&, CGAL::Tag_false){
  // Since Innermost_coefficient_type is not RealEmbeddable the following functors
  // should be CGAL::Null_functor.
  ASSERT_IS_NULL_FUNCTOR(typename PT::Sign_at);
  ASSERT_IS_NULL_FUNCTOR(typename PT::Sign_at_homogeneous);
  ASSERT_IS_NULL_FUNCTOR(typename PT::Compare);

}

// test functors depending on the Algebraic_category of ICoeff
template< class PT >
void
test_ac_icoeff_functors(const PT&, CGAL::Integral_domain_without_division_tag){
  ASSERT_IS_NULL_FUNCTOR(typename PT::Multivariate_content);
  //ASSERT_IS_NULL_FUNCTOR(typename PT::Interpolate);
}

template< class PT >
void test_ac_icoeff_functors(
    const PT& traits, CGAL::Unique_factorization_domain_tag){
  test_multivariate_content(traits);
  // ASSERT_IS_NULL_FUNCTOR(typename PT::Interpolate);
}
template< class PT >
void test_ac_icoeff_functors(const PT& traits, CGAL::Field_tag){
  test_multivariate_content(traits);
//  test_interpolate(traits);
}

// test functors depending on the Algebraic_category of Coefficient_type
template< class PT >
void test_ac_poly_functors(const PT&, CGAL::Integral_domain_without_division_tag){
  ASSERT_IS_NULL_FUNCTOR(typename PT::Univariate_content);
  ASSERT_IS_NULL_FUNCTOR(typename PT::Square_free_factorize);
}

template< class PT >
void test_ac_poly_functors(const PT& traits, CGAL::Unique_factorization_domain_tag){
  test_univariate_content(traits);
  test_square_free_factorize(traits);
}

template< class PT >
void test_polynomial_traits_d(const PT& traits){
  typedef typename PT::Polynomial_d           Polynomial_d;
  typedef typename PT::Innermost_coefficient_type  ICoeff;
  typedef typename PT::Coefficient_type            Coeff;
  CGAL_USE_TYPE(Coeff);
  test_fundamental_functors(traits);

  typedef typename CGAL::Algebraic_structure_traits<ICoeff> AST_IC;
  test_ac_icoeff_functors(traits, typename AST_IC::Algebraic_category());

  typedef typename CGAL::Algebraic_structure_traits<Polynomial_d> AST_Poly;
  test_ac_poly_functors(traits, typename AST_Poly::Algebraic_category());

  typedef typename CGAL::Real_embeddable_traits<ICoeff> RET_IC;
  typedef typename RET_IC::Is_real_embeddable Is_real_embeddable;
  test_real_embeddable_functors(traits, Is_real_embeddable());
  test_construct_coefficient_const_iterator_range(traits);
}

template <class  PT>
void test_rebind(const PT& /*traits*/){

  typedef typename PT::Innermost_coefficient_type IC;

{
  CGAL_assertion_code( const int dimension = 1; );
  typedef typename PT:: template Rebind<IC,1>::Other PT_IC_1;
  CGAL_USE_TYPE(PT_IC_1);
  CGAL_static_assertion((boost::is_same< typename PT_IC_1::Innermost_coefficient_type,
            IC>::value));
  CGAL_static_assertion((PT_IC_1::d==dimension));
}
{
  CGAL_assertion_code( const int dimension = 2; );
  typedef typename PT:: template Rebind<IC,2>::Other PT_IC_2;
  CGAL_USE_TYPE(PT_IC_2);
  CGAL_static_assertion((boost::is_same< typename PT_IC_2::Innermost_coefficient_type,
            IC>::value));
  CGAL_static_assertion((PT_IC_2::d==dimension));
}
{
  CGAL_assertion_code( const int dimension = 3; )
  typedef typename PT:: template Rebind<IC,3>::Other PT_IC_3;
  CGAL_USE_TYPE(PT_IC_3);
  CGAL_static_assertion((boost::is_same< typename PT_IC_3::Innermost_coefficient_type,
          IC>::value));
  CGAL_static_assertion((PT_IC_3::d==dimension));
}
{
  typedef typename PT:: template Rebind<IC,1>::Other PT_IC_1;
  typedef typename PT:: template Rebind<IC,2>::Other PT_IC_2;
  typedef typename PT:: template Rebind<IC,3>::Other PT_IC_3;

  CGAL_USE_TYPE(PT_IC_1);
  CGAL_USE_TYPE(PT_IC_2);
  CGAL_USE_TYPE(PT_IC_3);
  CGAL_assertion_code(typedef typename  PT_IC_1::Polynomial_d Poly1;)
  CGAL_assertion_code(typedef typename  PT_IC_2::Polynomial_d Poly2;)

  CGAL_static_assertion((boost::is_same< typename PT_IC_1::Coefficient_type,
          IC>::value));
  CGAL_static_assertion((boost::is_same< typename PT_IC_2::Coefficient_type,
          Poly1>::value));
  CGAL_static_assertion((boost::is_same< typename PT_IC_3::Coefficient_type,
          Poly2>::value));

}

#ifdef CGAL_USE_LEDA
{
  typedef CGAL::LEDA_arithmetic_kernel AT;
  typedef typename AT::Integer Integer;
  typedef typename AT::Rational Rational;
  const int dimension = 4; CGAL_USE(dimension);
  typedef typename PT:: template Rebind<Integer,4>::Other PT_Integer_4;
  CGAL_USE_TYPE(PT_Integer_4);
  typedef typename PT:: template Rebind<Rational,4>::Other PT_Rational_4;
  CGAL_USE_TYPE(PT_Rational_4);
  CGAL_static_assertion((boost::is_same< typename PT_Integer_4::Innermost_coefficient_type,
          Integer>::value));
  CGAL_static_assertion((boost::is_same< typename PT_Rational_4::Innermost_coefficient_type,
          Rational>::value));
  CGAL_static_assertion((PT_Integer_4::d==dimension));
  CGAL_static_assertion((PT_Rational_4::d==dimension));
}
#endif
#ifdef CGAL_USE_CORE
{
  typedef CGAL::CORE_arithmetic_kernel AT;
  typedef typename AT::Integer Integer;
  typedef typename AT::Rational Rational;
  CGAL_assertion_code( const int dimension = 4; )
  typedef typename PT:: template Rebind<Integer,4>::Other PT_Integer_4;
  typedef typename PT:: template Rebind<Rational,4>::Other PT_Rational_4;
  CGAL_USE_TYPE(PT_Integer_4);
  CGAL_USE_TYPE(PT_Rational_4);
  CGAL_static_assertion((boost::is_same< typename PT_Integer_4::Innermost_coefficient_type,
          Integer>::value));
  CGAL_static_assertion((boost::is_same< typename PT_Rational_4::Innermost_coefficient_type,
          Rational>::value));
  CGAL_static_assertion((PT_Integer_4::d==dimension));
  CGAL_static_assertion((PT_Rational_4::d==dimension));
}
#endif
{
  CGAL_assertion_code( const int dimension = 4; )
  typedef typename PT:: template Rebind<int,4>::Other PT_Integer_4;
  typedef typename PT:: template Rebind<double,4>::Other PT_Rational_4;
  CGAL_USE_TYPE(PT_Integer_4);
  CGAL_USE_TYPE(PT_Rational_4);
  CGAL_static_assertion((boost::is_same< typename PT_Integer_4::Innermost_coefficient_type,
          int>::value));
  CGAL_static_assertion((boost::is_same< typename PT_Rational_4::Innermost_coefficient_type,
          double>::value));
  CGAL_static_assertion((PT_Integer_4::d==dimension));
  CGAL_static_assertion((PT_Rational_4::d==dimension));
}
}



template< class PT >
void test_multiple_dimensions(const PT& traits) {
    test_rebind(traits);

    typedef typename PT::Innermost_coefficient_type IC;
    typedef typename PT:: template Rebind<IC,1>::Other PT_IC_1;
    typedef typename PT:: template Rebind<IC,2>::Other PT_IC_2;
    typedef typename PT:: template Rebind<IC,3>::Other PT_IC_3;

    test_permute(traits);
    test_construct_innermost_coefficient_const_iterator_range(traits);

    Test_Pol::test_polynomial_traits_d(PT_IC_1());
    Test_Pol::test_polynomial_traits_d(PT_IC_2());
    Test_Pol::test_polynomial_traits_d(PT_IC_3());
}
}//Namespace Test_Pol

template< class PT >
void test_polynomial_traits_d(const PT& traits){
  Test_Pol::test_polynomial_traits_d(traits);
}

} // namespace CGAL

#endif
