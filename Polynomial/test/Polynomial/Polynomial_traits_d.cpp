// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hemmer <hemmer@informatik.uni-mainz.de> 
//                 Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//
// ============================================================================

#include <iostream>
#include <CGAL/basic.h>

#include <CGAL/Sqrt_extension.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/ipower.h>

#include <cassert>

#include <CGAL/Test/_test_basic.h>

#include <CGAL/Random.h>

static CGAL::Random my_rnd(346); // some seed 

#define CGAL_SNAP_CGALi_TRAITS_D(T)                             \
  typedef T PT;                                                 \
  typedef typename PT::Polynomial_d          Polynomial_d;      \
  typedef typename PT::Coefficient           Coeff;             \
  typedef typename PT::Innermost_coefficient ICoeff;            \
  typedef CGAL::Polynomial_traits_d<Coeff> PTC;                 \
  typedef CGAL::Exponent_vector Exponent_vector;                \
  typedef std::pair< CGAL::Exponent_vector , ICoeff > Monom;    \
  typedef std::vector< Monom > Monom_rep;    

#define ASSERT_IS_NULL_FUNCTOR(T)                                       \
  BOOST_STATIC_ASSERT((boost::is_same<T,CGAL::Null_functor >::value))   

template <class Polynomial_d_>
Polynomial_d_
generate_sparse_random_polynomial(int max_degree = 10){
  typedef CGAL::Polynomial_traits_d<Polynomial_d_> PT;
  CGAL_SNAP_CGALi_TRAITS_D(PT);
  typename PT::Construct_polynomial construct; 


  typedef CGAL::Exponent_vector Exponent_vector;
  typedef std::pair< CGAL::Exponent_vector , ICoeff > Monom;
  typedef std::vector< Monom > Monom_rep;

  int range = 20;
  int number_of_variables = PT::d;
  int number_of_coeffs = 
    CGAL::min(number_of_variables * (int)ceil(log(max_degree+1))+1,4);
    
  Polynomial_d result; 
  for(int i = 0; i < number_of_coeffs; i++){
    CGAL::Exponent_vector exps(PT::d);
    for(int j = 0; j < PT::d; j++){
      exps[j]=my_rnd.get_int(0,max_degree);
    }
    ICoeff c = ICoeff(my_rnd.get_int(-range,range));
    Monom_rep monom_rep;
    monom_rep.push_back(Monom(exps,c));
    result += construct(monom_rep.begin(), monom_rep.end());
  }
    
  return result;
}


template <class Polynomial_traits_d>
void test_construct_polynomial(const Polynomial_traits_d&){
  std::cerr << "start test_construct_polynomial "; 
  std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
    
  const int d = PT::d;
  { // Construct_polynomial
    typedef typename PT::Construct_polynomial Constructor;
    BOOST_STATIC_ASSERT(
        !(boost::is_same< Constructor , CGAL::Null_functor >::value));
    typedef typename Constructor::result_type result_type;
    BOOST_STATIC_ASSERT(
        (boost::is_same< result_type , Polynomial_d >::value));
    Constructor construct;

    assert(Polynomial_d() == construct()); 
    assert(Polynomial_d(3) == construct(3)); 
    assert(Polynomial_d(Coeff(3)) == construct(Coeff(3))); 
    assert(Polynomial_d(ICoeff(3)) == construct(ICoeff(3))); 

    assert(construct(Coeff(2)) != Polynomial_d(0)); 
    assert(construct(Coeff(2)) == Polynomial_d(2)); 
    assert(construct(Coeff(0),Coeff(1)) != Polynomial_d(1)); 
    assert(construct() 
        == construct(Coeff(0)));
    assert(construct(Coeff(1)) 
        == construct(Coeff(1),Coeff(0)));
    assert(construct(Coeff(2),Coeff(1)) 
        == construct(Coeff(2),Coeff(1),Coeff(0)));
    assert(construct(Coeff(3),Coeff(2),Coeff(1)) 
        == construct(Coeff(3),Coeff(2),Coeff(1),Coeff(0)));
    assert(construct(Coeff(3),Coeff(2),Coeff(1)) 
        != construct(Coeff(3),Coeff(2),Coeff(1),Coeff(1)));
    // construct via iterator range
    std::vector<Coeff> coeffs;
    assert(construct(coeffs.begin(),coeffs.end()) == construct(0));
    for(int i = 0; i<4;i++){coeffs.push_back(Coeff(i));}
    assert(construct(coeffs.begin(),coeffs.end())
        == construct(Coeff(0),Coeff(1),Coeff(2),Coeff(3)));
        
    Monom_rep monom_rep;
    assert(
        construct(monom_rep.begin(),monom_rep.end()) == construct(0));
    CGAL::Random rnd(7);
    for(int j = 0; j < 2; j++){
      CGAL::Exponent_vector exps(d);
      for(int i = 0; i < d; i++){
        exps[i]=j+i*5;
      }
      monom_rep.push_back(Monom(exps,ICoeff(j+1)));
    }
    //std::cout<<"\n"<<std::endl;
    //std::cout<<"monom_rep_1: "; print_monom_rep(monom_rep); 
       
    std::random_shuffle(monom_rep. begin(),monom_rep. end());
    Polynomial_d p1 = construct(monom_rep. begin(),
        monom_rep. begin()+((monom_rep. end()- monom_rep. begin())/2));
    Polynomial_d p2 = construct(monom_rep. begin()+
        ((monom_rep. end()- monom_rep. begin())/2), monom_rep. end());
    Polynomial_d p  = construct(monom_rep. begin(), monom_rep. end());
    assert(p == p1+p2);
    
    assert(construct(monom_rep. begin(),monom_rep. end()) 
        == construct(monom_rep.rbegin(),monom_rep.rend()));  
  }
  std::cerr << " ok "<< std::endl; 
}

template< class Polynomial_traits_d >
void test_get_coefficient(const Polynomial_traits_d&) {
  std::cerr << "start test_get_coefficient ";
  std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
    
  typename PT::Construct_polynomial construct;
  typename PT::Get_coefficient get_coeff;
    
  Polynomial_d p = construct(Coeff(1), Coeff(2), Coeff(3));
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
    
  typename PT::Construct_polynomial construct;
  typename PT::Get_innermost_coefficient get_innermost_coeff;
    
  Polynomial_d p = construct(Coeff(1), Coeff(2), Coeff(3));
    
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
    
  typename PT::Construct_polynomial construct;
  typename PT::Get_monom_representation gmr;
        
  for (int i = 0; i < 5 ; i++){
    Polynomial_d p,q;
    Monom_rep monom_rep;
    p = generate_sparse_random_polynomial<Polynomial_d>();
    gmr(p,std::back_inserter(monom_rep));
    q = construct(monom_rep.begin(), monom_rep.end());
    assert(q == p);
    std::random_shuffle(monom_rep.begin(), monom_rep.end());
    q = construct(monom_rep.begin(), monom_rep.end());
    assert(q == p);
  }
  std::cerr << " ok "<< std::endl; 
}




template <class Polynomial_traits_d>
void test_swap(const Polynomial_traits_d&){
  std::cerr << "start test_swap "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  int d = PT::d;
  typename Polynomial_traits_d::Swap swap;
    
  //std::cout << "start_test ----------- "<< d << std::endl; 
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
void test_move(const Polynomial_traits_d&){
  std::cerr << "start test_move "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename Polynomial_traits_d::Move move;
  typename Polynomial_traits_d::Swap swap;
    
  //std::cout << "start_test ----------- "<< d << std::endl; 
  for(int i = 0; i < 5; i++){
    int n = my_rnd.get_int(0,PT::d-1);
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

  typename PT::Construct_polynomial construct;
  typename PT::Degree degree;

  CGAL::Test_functor_arity< typename PT::Degree >()(1);
    
  Polynomial_d p; 
  p= construct(Coeff(0));
  assert(degree(p) == 0);
  p= construct(Coeff(1));
  assert(degree(p) == 0);
  p= construct(Coeff(1),Coeff(2));
  assert(degree(p) == 1);

    
  p= construct(Coeff(0));
  assert(degree(p,(PT::d-1)) == 0);
  p= construct(Coeff(1));
  assert(degree(p,(PT::d-1)) == 0);
  p= construct(Coeff(1),Coeff(2));
  assert(degree(p,(PT::d-1)) == 1);
  std::cerr << " ok "<< std::endl; 
}


//       Total_degree;
template <class Polynomial_traits_d>
void test_total_degree(const Polynomial_traits_d&){

  std::cerr << "start test_total_degree "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Total_degree total_degree;
  for(int i = 0; i < 5; i++){
    Polynomial_d p =generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q =generate_sparse_random_polynomial<Polynomial_d>();
    int tdp = total_degree(p);
    int tdq = total_degree(q);
    assert(total_degree(p*q) == tdp+tdq);    
  }
  std::cerr << " ok "<< std::endl; 
}
// //       Leading_coefficient;
template <class Polynomial_traits_d>
void test_leading_coefficient(const Polynomial_traits_d&){
    
  std::cerr << "start test_leading_coefficient "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Leading_coefficient lcoeff;
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
    
  typename PT::Innermost_leading_coefficient ilcoeff;
  Polynomial_d p(Coeff(1), Coeff(2), Coeff(3));
    
  assert(ilcoeff(p) == ICoeff(3));
    
  std::cerr << " ok" << std::endl;
}

// //       Univariate_content;
template <class Polynomial_traits_d>
void test_univariate_content(const Polynomial_traits_d&){
    
  std::cerr << "start test_univariate_content "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Univariate_content univariate_content;
    
  assert(univariate_content(Polynomial_d(0)) == Coeff(0));
  assert(univariate_content(Polynomial_d(1)) == Coeff(1));
  assert(
      univariate_content(Polynomial_d(2)) 
      == 
      CGAL::integral_division(Coeff(2), CGAL::unit_part(Coeff(2))));
    
    
  for(int i = 0; i < 5; i++){
    Polynomial_d p =generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q =generate_sparse_random_polynomial<Polynomial_d>();
    Coeff ucontentp = univariate_content(p);
    Coeff ucontentq = univariate_content(q);
    assert(univariate_content(p*q) == ucontentp*ucontentq);     
  }
  for(int i = 0; i < 5; i++){
    Polynomial_d p =generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q =generate_sparse_random_polynomial<Polynomial_d>();
    Coeff ucontentp = univariate_content(p,0);
    Coeff ucontentq = univariate_content(q,0);
    assert(univariate_content(p*q,0) == ucontentp*ucontentq);    
  }
  std::cerr << " ok "<< std::endl; 
}
// //       Multivariate_content;
template <class Polynomial_traits_d>
void test_multivariate_content(const Polynomial_traits_d&){
  std::cerr << "start test_multivariate_content ";
  std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Multivariate_content mcontent;
    
  assert(mcontent(Polynomial_d(0)) == 
      CGAL::integral_division(ICoeff(0),CGAL::unit_part(ICoeff(0))));
  assert(mcontent(Polynomial_d(1)) == 
      CGAL::integral_division(ICoeff(1), CGAL::unit_part(ICoeff(1))));
  assert(mcontent(Polynomial_d(2)) == 
      CGAL::integral_division(ICoeff(2), CGAL::unit_part(ICoeff(2))));
  assert(mcontent(Polynomial_d(-2)) == 
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

// //       Multivariate_content;
template <class Polynomial_traits_d>
void test_interpolate(const Polynomial_traits_d&){
  std::cerr << "start test_interpolate "; std::cerr.flush();
  typedef Polynomial_traits_d PT_d; 
  typedef typename PT_d::Innermost_coefficient ICoeff;
  typedef typename PT_d::Polynomial_d Polynomial_d;
  typedef typename PT_d:: template Rebind<ICoeff,1>::Other PT_1;
  typedef typename PT_1::Polynomial_d Polynomial_1;
    
  typename PT_d::Interpolate interpolate;
  typename PT_d::Evaluate eval; 
    
  for(int i = 0; i < 5; i++){
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>(i);
       
    Polynomial_1 m(1); 
    Polynomial_d u(0);
    for (int j = 0; j <= i; j++){ 
      Polynomial_1 m1 = m; 
      Polynomial_d u1 = u; 
      Polynomial_1 m2 = Polynomial_1(ICoeff(-j),ICoeff(2));
      Polynomial_d u2 = eval(p,ICoeff(j)/ICoeff(2));
      interpolate(m1,u1,m2,u2,m,u);
    }
    assert(u == p);
  }
    
  for(int i = 0; i < 5; i++){
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>(i);
       
    Polynomial_1 m(1); 
    Polynomial_d u(0);
    for (int j = 0; j <= i; j++){ 
      Polynomial_1 m2 = m; 
      Polynomial_d u2 = u; 
      Polynomial_1 m1 = Polynomial_1(ICoeff(-j),ICoeff(2));
      Polynomial_d u1 = eval(p,ICoeff(j)/ICoeff(2));
      interpolate(m1,u1,m2,u2,m,u);
    }
    assert(u == p);
  }
  std::cerr << " ok "<< std::endl; 
}

// //       Shift;
template <class Polynomial_traits_d>
void test_shift(const Polynomial_traits_d&){
    
  std::cerr << "start test_shift "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Shift shift; 
  typename PT::Construct_polynomial construct;
  for(int i = 0; i < 5; i++){ 
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q = p*CGAL::ipower(construct(Coeff(0),Coeff(1)),5);
    assert(shift(p,5) == q);   
  }
  std::cerr << " ok "<< std::endl; 
}
// //       Negate;
template <class Polynomial_traits_d>
void test_negate(const Polynomial_traits_d&){
  std::cerr << "start test_negate "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Construct_polynomial construct;
  typename PT::Negate negate;
  typename PT::Swap swap;
    
  assert(negate(Polynomial_d(0)) == Polynomial_d(0));
  assert(negate(Polynomial_d(2)) == Polynomial_d(2));

  Polynomial_d p = construct(Coeff(1),Coeff(2),Coeff(3));
  assert(negate(p) == construct(Coeff(1),Coeff(-2),Coeff(3)));;
    
  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q = p;
    int n = my_rnd.get_int(0,PT::d-1);
    int m = my_rnd.get_int(0,PT::d-1);
    assert(negate(swap(p,n,m),n) == swap(negate(p,m),n,m));
  }
  std::cerr << " ok "<< std::endl;
}
//      Invert;
template <class Polynomial_traits_d>
void test_invert(const Polynomial_traits_d&){
  std::cerr << "start test_invert "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Invert invert; 
  typename PT::Swap swap; 
  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
    std::vector<Coeff> coeffs (p.begin(),p.end());
    p = invert(p);
    std::vector<Coeff> rcoeffs (p.begin(),p.end());
    assert(coeffs.size() >= rcoeffs.size());
    for (unsigned int i = 0; i < rcoeffs.size(); i++){
      assert(rcoeffs[i] == coeffs[coeffs.size()-i-1]);
    }   
    int n = my_rnd.get_int(0,PT::d-1);
        
    assert(
        invert(p,n) == swap(invert(swap(p,n,PT::d-1)),n,PT::d-1));
  }
  std::cerr << " ok "<< std::endl;
}
// //       Translate;
template <class Polynomial_traits_d>
void test_translate(const Polynomial_traits_d&){
  std::cerr << "start test_translate "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Translate translate;
  typename PT::Evaluate evaluate;
  typename PT::Move move;
  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
    assert(evaluate(translate(p,ICoeff(5)),ICoeff(3)) 
        == evaluate(p,ICoeff(8)));
    assert(evaluate(translate(p,ICoeff(5),0),ICoeff(3),0) 
        == evaluate(move(p,0,PT::d-1),ICoeff(8)));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Translate_homogeneous;
template <class Polynomial_traits_d>
void test_translate_homongenous(const Polynomial_traits_d&){
  std::cerr << "start test_translate_homongenous "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Translate_homogeneous transh;
  typename PT::Canonicalize canonicalize;
  typename PT::Evaluate_homogeneous evh;
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
    
  typename PT::Scale scale;
  Polynomial_d p(Coeff(1), Coeff(2), Coeff(3));
    
  assert(
      scale(p, ICoeff(2)) == Polynomial_d(Coeff(1), Coeff(4), Coeff(12)));
    
  std::cerr << " ok" << std::endl;
}

// //       Scale_homogeneous;
template <class Polynomial_traits_d>
void test_scale_homogeneous(const Polynomial_traits_d&){
  std::cerr << "start test_scale_homogeneous "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Scale_homogeneous scaleh;
  typename PT::Canonicalize canonicalize;
  //typename PT::Move move;
  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p,q1,q2;
    p = generate_sparse_random_polynomial<Polynomial_d>(); 
    q1 = scaleh(scaleh(p,ICoeff(5),ICoeff(3)),ICoeff(3),ICoeff(2)) ;
    q2 = scaleh(p,ICoeff(15),ICoeff(6)) ;
    assert(canonicalize(q1) != canonicalize(p)) ;
    assert(canonicalize(q2) != canonicalize(p)) ;
    assert(canonicalize(q1) == canonicalize(q2));
        
    q1 = scaleh(scaleh(p,ICoeff(5),ICoeff(3),1),ICoeff(3),ICoeff(2),0) ;  
    q2 = scaleh(p,ICoeff(15),ICoeff(6),0) ;
        
    assert(canonicalize(q1) != canonicalize(p)) ;
    assert(canonicalize(q2) != canonicalize(p)) ;
    assert(canonicalize(q1) == canonicalize(q2)); 
  }
  std::cerr << " ok "<< std::endl;
}

// //       Derivative;
template <class Polynomial_traits_d>
void test_derivative(const Polynomial_traits_d&){
  std::cerr << "start test_derivative "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
    
  typename PT::Derivative diff;
  typename PT::Swap swap;
    
  assert(diff(Polynomial_d(0)) == Polynomial_d(0));
  assert(diff(Polynomial_d(1)) == Polynomial_d(0));
  assert(diff(Polynomial_d(Coeff(1), Coeff(2))) == Polynomial_d(2));
    
  for(int i = 0 ; i < 5 ; i++){
    int n = my_rnd.get_int(0,PT::d-1);
    Polynomial_d p,pd;
    p = generate_sparse_random_polynomial<Polynomial_d>(); 
    pd = diff(p,n);
    assert(pd == swap(diff(swap(p,n,PT::d-1)),n,PT::d-1));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Make_square_free;
template <class Polynomial_traits_d>
void test_make_square_free(const Polynomial_traits_d&){
  std::cerr << "start test_make_square_free "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Make_square_free make_square_free;
  typename PT::Leading_coefficient lcoeff;
  typename PT::Univariate_content_up_to_constant_factor ucontent_utcf;
  typename PT::Integral_division_up_to_constant_factor idiv_utcf;
    
  assert(Polynomial_d(0) == make_square_free(Polynomial_d(0)));
  assert(Polynomial_d(1) == make_square_free(Polynomial_d(1)));
  assert(Polynomial_d(1) == make_square_free(Polynomial_d(2)));
    
  //typename PT::Canonicalize canonicalize; 
  for(int i = 0 ; i < 5 ; i++){
    Polynomial_d p;
    p = generate_sparse_random_polynomial<Polynomial_d>(3); 
        
    p = idiv_utcf(p, ucontent_utcf(p));
    p = make_square_free(p);
    Coeff f_intern = lcoeff(p);
        
    f_intern = typename PTC::Make_square_free()(f_intern);
       
    //std::cout <<"f: " << f_intern << std::endl;
        
    assert(p * f_intern == make_square_free(p*p*f_intern));
    assert(p * f_intern == make_square_free(p*f_intern*f_intern));
  }
  std::cerr << " ok "<< std::endl;
}


// //       Square_free_factorization;
template <class Polynomial_traits_d>
void test_square_free_factorization(const Polynomial_traits_d&){
  std::cerr << "start test_square_free_factorization "; std::cerr.flush(); 
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef CGAL::Algebraic_structure_traits<Polynomial_d> AST;
  typename AST::Integral_division idiv;
  typename PT::Square_free_factorization sqff;
  typename PT::Canonicalize canonicalize;

  for(int i = 0; i < 5; i++){
    Polynomial_d f1 = generate_sparse_random_polynomial<Polynomial_d>(2);
    Polynomial_d f2 = generate_sparse_random_polynomial<Polynomial_d>(2);    
    Polynomial_d p = f1*f1*f2;
    std::vector<std::pair<Polynomial_d,int> > fac_mul_pairs;
    sqff(p, std::back_inserter(fac_mul_pairs));
    int n = fac_mul_pairs.size();
    for (int j = 0; j < n; j++){
      Polynomial_d factor = fac_mul_pairs[j].first;
      int multi = fac_mul_pairs[j].second;
      for (int k = 0; k < multi; k++){
        p = idiv(p,factor);
      }
    }
    assert(CGAL::is_one(canonicalize(p)));
  }
    
  
  typename PT::Innermost_leading_coefficient ileading_coeff;
  typename PT::Multivariate_content multivariate_content;
  Polynomial_d p = generate_sparse_random_polynomial< Polynomial_d >(2);
  p *= p*generate_sparse_random_polynomial< Polynomial_d >(2);    
  std::vector< std::pair< Polynomial_d,int> > fac_mul_pairs;

  ICoeff alpha;
 
  
  sqff(p, std::back_inserter(fac_mul_pairs), alpha);

  assert(alpha 
      == CGAL::unit_part(ileading_coeff(p)) * multivariate_content(p));
    
  //std::cerr << std::endl;
  Polynomial_d rec_p(alpha);
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
  typename PT::Pseudo_division pdiv;
  for(int i = 0; i < 10; i++){
    Polynomial_d f = generate_sparse_random_polynomial<Polynomial_d>(3);
    Polynomial_d g = generate_sparse_random_polynomial<Polynomial_d>(2);    
    Coeff D;
    Polynomial_d q,r;
    pdiv(f,g,q,r,D);
    assert(f*Polynomial_d(D) == g*q+r);
  }
  std::cerr << " ok "<< std::endl;
}

// //       Pseudo_division_remainder;
template <class Polynomial_traits_d>
void test_pseudo_division_remainder(const Polynomial_traits_d&){
  std::cerr << "start test_pseudo_division_remainder "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Pseudo_division pdiv;
  typename PT::Pseudo_division_remainder pdiv_r;
  for(int i = 0; i < 10; i++){
    Polynomial_d f = generate_sparse_random_polynomial<Polynomial_d>(3);
    Polynomial_d g = generate_sparse_random_polynomial<Polynomial_d>(2);    
    Coeff D;
    Polynomial_d q,r;
    pdiv(f,g,q,r,D);
    assert(r == pdiv_r(f,g));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Pseudo_division_quotient;
template <class Polynomial_traits_d>
void test_pseudo_division_quotient(const Polynomial_traits_d&){
  std::cerr << "start test_pseudo_division_quotient "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Pseudo_division pdiv;
  typename PT::Pseudo_division_quotient pdiv_q;
  for(int i = 0; i < 10; i++){
    Polynomial_d f = generate_sparse_random_polynomial<Polynomial_d>(3);
    Polynomial_d g = generate_sparse_random_polynomial<Polynomial_d>(2);    
    Coeff D;
    Polynomial_d q,r;
    pdiv(f,g,q,r,D);
    assert(q == pdiv_q(f,g));
  }
  std::cerr << " ok "<< std::endl;
}

// //       Gcd_up_to_constant_factor;
template <class Polynomial_traits_d>
void test_gcd_up_to_constant_factor(const Polynomial_traits_d&){
  std::cerr << "start test_gcd_up_to_constant_factor "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Gcd_up_to_constant_factor gcd_utcf;
    
  assert(
      Polynomial_d(0) == gcd_utcf(Polynomial_d(0),Polynomial_d(0)));
  assert(
      Polynomial_d(1) == gcd_utcf(Polynomial_d(1),Polynomial_d(0)));
  assert(
      Polynomial_d(1) == gcd_utcf(Polynomial_d(0),Polynomial_d(1)));
  assert(
      Polynomial_d(1) == gcd_utcf(Polynomial_d(1),Polynomial_d(1)));
  assert(
      Polynomial_d(1) == gcd_utcf(Polynomial_d(-1),Polynomial_d(-1)));
  assert(
      Polynomial_d(1) == gcd_utcf(Polynomial_d(2),Polynomial_d(2)));
    
  std::cerr << " ok "<< std::endl;
}

// //       Integral_division_up_to_constant_factor;
template <class Polynomial_traits_d>
void test_integral_division_up_to_constant_factor(const Polynomial_traits_d&){
  std::cerr << "start test_integral_division_up_to_constant_factor "; 
  std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Integral_division_up_to_constant_factor idiv_utcf;
  typename PT::Canonicalize canonicalize;
  assert(
      Polynomial_d(0) == idiv_utcf(Polynomial_d(0),Polynomial_d(1)));
  assert(
      Polynomial_d(1) == idiv_utcf(Polynomial_d(1),Polynomial_d(1)));
  assert(
      Polynomial_d(1) == idiv_utcf(Polynomial_d(2),Polynomial_d(1)));
    
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
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Univariate_content_up_to_constant_factor ucontent_utcf;
  typename PT::Integral_division_up_to_constant_factor idiv_utcf;
  typename PT::Leading_coefficient lcoeff;
    
  typename PT::Canonicalize canonicalize;
  assert(Coeff(0) == ucontent_utcf(Polynomial_d(0)));
  assert(Coeff(1) == ucontent_utcf(Polynomial_d(1)));
  assert(Coeff(1) == ucontent_utcf(Polynomial_d(2)));
  assert(Coeff(1) == ucontent_utcf(Polynomial_d(-2)));
    
  for(int i = 0; i < 5; i++){
    Polynomial_d p,q;
    p = generate_sparse_random_polynomial<Polynomial_d>(3);
    Coeff content = ucontent_utcf(p);
    p = idiv_utcf(p,Polynomial_d(content));
    assert(Coeff(1) == ucontent_utcf(p));
    Coeff lc = lcoeff(p);
    p = p*lc*lc;
    assert(canonicalize(Polynomial_d(lc*lc)) == ucontent_utcf(p));
        
  }
    
  std::cerr << " ok "<< std::endl;
}


// //       Square_free_factorization_up_to_constant_factor;
template <class Polynomial_traits_d>
void test_square_free_factorization_up_to_constant_factor(const Polynomial_traits_d&){
  std::cerr << "start test_square_free_factorization_up_to_constant_factor "; 
  std::cerr.flush(); 
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Integral_division_up_to_constant_factor idiv_utcf;
  typename PT::Square_free_factorization_up_to_constant_factor sqff_utcf;
  typename PT::Canonicalize canonicalize;


    for(int i = 0; i < 5; i++){
    Polynomial_d f1 = generate_sparse_random_polynomial<Polynomial_d>(2);
    Polynomial_d f2 = generate_sparse_random_polynomial<Polynomial_d>(2);    
    Polynomial_d p = f1*f1*f2;
    std::vector<std::pair<Polynomial_d,int> > fac_mul_pairs;
    sqff_utcf(p, std::back_inserter(fac_mul_pairs)) ;
     
    int n = fac_mul_pairs.size();
   
    for (int j = 0; j < n; j++){
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
  typename PT::Evaluate evaluate;
  typename PT::Move move;
  assert(evaluate(Polynomial_d(0),ICoeff(0)) == Coeff(0));
  assert(evaluate(Polynomial_d(1),ICoeff(0)) == Coeff(1));
  assert(evaluate(Polynomial_d(2),ICoeff(5)) == Coeff(2));

  assert( evaluate(Polynomial_d(Coeff(3),Coeff(2)),ICoeff(0)) == Coeff(3));
  assert( evaluate(Polynomial_d(Coeff(3),Coeff(2)),ICoeff(1)) == Coeff(5));
  assert( evaluate(Polynomial_d(Coeff(3),Coeff(2)),ICoeff(2)) == Coeff(7));
    
  for(int i = 0; i < 5; i++){
    int n = my_rnd.get_int(0,PT::d-1);
    Polynomial_d p,q;
    p = generate_sparse_random_polynomial<Polynomial_d>();
    assert(evaluate(p,ICoeff(3),n) 
        == evaluate(move(p,n,PT::d-1),ICoeff(3)));
  }
    
  std::cerr << " ok "<< std::endl;
}

// //       Evaluate_homogeneous;
template <class Polynomial_traits_d>
void test_evaluate_homogeneous(const Polynomial_traits_d&){
  std::cerr << "start test_evaluate_homogeneous "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Evaluate_homogeneous evh;
  
  assert(evh(Polynomial_d(0),ICoeff(0),ICoeff(1)) == Coeff(0));
  assert(evh(Polynomial_d(1),ICoeff(0),ICoeff(2)) == Coeff(1));
  assert(evh(Polynomial_d(2),ICoeff(5),ICoeff(3)) == Coeff(2));

  assert(evh( Polynomial_d(Coeff(3),Coeff(2)) , ICoeff(0),ICoeff(1)) 
      == Coeff(3));
  assert(evh( Polynomial_d(Coeff(3),Coeff(2)) , ICoeff(1),ICoeff(1)) 
      == Coeff(5));
  assert(evh( Polynomial_d(Coeff(3),Coeff(2)) , ICoeff(2),ICoeff(3)) 
      == Coeff(9+4));
    
  std::cerr << " ok "<< std::endl;
}

template< class Polynomial_traits_d >
void test_is_zero_at(const Polynomial_traits_d&) {
  std::cerr << "start test_is_zero_at "; 
  std::cerr.flush();
    
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
    
  typename PT::Is_zero_at is_zero_at;
    
  Polynomial_d p(Coeff(-1), Coeff(0), Coeff(1));
    
  std::vector< ICoeff > cv;
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
    
  typename PT::Is_zero_at_homogeneous is_zero_at_homogeneous;
    
  Polynomial_d p(Coeff(-1), Coeff(0), Coeff(4));
    
  std::vector< ICoeff > cv;
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
    
  typename PT::Sign_at sign_at;
    
  Polynomial_d p(Coeff(-1), Coeff(0), Coeff(1));
    
  std::vector< ICoeff > cv;
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
    
  typename PT::Sign_at_homogeneous sign_at_homogeneous;
    
  Polynomial_d p(Coeff(-1), Coeff(0), Coeff(1));
    
  std::vector< ICoeff > cv;
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
    
  typename PT::Compare compare;
    
  Polynomial_d p0(Coeff(0));
  Polynomial_d pp2(Coeff(2));
  Polynomial_d pm2(Coeff(-2));
  Polynomial_d pp1p2(Coeff(1), Coeff(2));
  Polynomial_d pm1m2(Coeff(-1), Coeff(-2));
    
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
  typename  PT::Resultant resultant;
  typename  PT::Move move;
  {
    Polynomial_d A(0);
    Polynomial_d B(0);
    assert(resultant(A,B) == Coeff(0));
  }{
    Polynomial_d A(4);
    Polynomial_d B(8);
    assert(resultant(A,B) == Coeff(1));
  }{
    Polynomial_d f(Coeff(2),Coeff(7),Coeff(1),Coeff(8),Coeff(1),Coeff(8));
    Polynomial_d g(Coeff(3),Coeff(1),Coeff(4),Coeff(1),Coeff(5),Coeff(9));
        
    assert(resultant(f,g) == Coeff(230664271L)); // Maple
            
    Polynomial_d h(Coeff(3),Coeff(4),Coeff(7),Coeff(7));
    Polynomial_d fh(f*h);
    Polynomial_d gh(g*h);
    assert(resultant(fh,gh) == Coeff(0));
  } 
  for(int i = 0 ; i < 5 ; i++){
    int n = my_rnd.get_int(0,PT::d-1);
    Polynomial_d p,q;
    p = generate_sparse_random_polynomial<Polynomial_d>(3);
    q = generate_sparse_random_polynomial<Polynomial_d>(3);
    assert(resultant(p,q,n) == resultant(move(p,n),move(q,n)));
  }
    
  std::cerr << " ok "<< std::endl;    
}

// //       Canonicalize;
template <class Polynomial_traits_d>
void test_canonicalize(const Polynomial_traits_d&){
  std::cerr << "start test_canonicalize "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Canonicalize canonicalize;

  assert(Polynomial_d(0) == canonicalize(Polynomial_d(0)));
  assert(Polynomial_d(1) == canonicalize(Polynomial_d(1)));
  assert(Polynomial_d(1) == canonicalize(Polynomial_d(2)));
  assert(Polynomial_d(1) == canonicalize(Polynomial_d(-2)));

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
  typename PT::Substitute substitute;
  typedef typename PT::Innermost_coefficient Innermost_coefficient;
    
    
  std::vector<Innermost_coefficient> vec;
  for(int i = 0; i < PT::d; i++){
    vec.push_back(Innermost_coefficient(i));
  }
  assert(Innermost_coefficient(0) 
      == substitute(Polynomial_d(0),vec.begin(),vec.end()));
  assert(Innermost_coefficient(1) 
      == substitute(Polynomial_d(1),vec.begin(),vec.end()));
  assert(Innermost_coefficient(2) 
      == substitute(Polynomial_d(2),vec.begin(),vec.end()));
  assert(Innermost_coefficient(-2) 
      == substitute(Polynomial_d(-2),vec.begin(),vec.end())); 
    
  for(int i = 0; i< 5; i++){
    typedef typename PT
      :: template Rebind<Innermost_coefficient,2>::Other PT_2;
    typedef typename PT_2::Polynomial_d Polynomial_2;
    std::vector<Polynomial_2> vec1,vec2;
    for(int j = 0; j < PT::d; j++){
      vec1.push_back(
          generate_sparse_random_polynomial<Polynomial_2>(3));
    }
    vec2=vec1;
    std::swap(vec2[0],vec2[PT::d-1]);
    Polynomial_d p 
      = generate_sparse_random_polynomial<Polynomial_d>(3); 
    assert( substitute(p,vec1.begin(),vec1.end()) == 
        substitute(typename PT::Swap()(p,0,PT::d-1),
            vec2.begin(),vec2.end()));
  }   
    
  std::cerr << " ok "<< std::endl;
}

// //       Substitute;
template <class Polynomial_traits_d>
void test_substitute_homogeneous(const Polynomial_traits_d&){
  std::cerr << "start test_substitute_homogeneous "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Substitute_homogeneous substitute_homogeneous;
  typedef typename PT::Innermost_coefficient Innermost_coefficient;
    
    
  std::vector<Innermost_coefficient> vec;
  for(int i = 0; i < PT::d; i++){
    vec.push_back(Innermost_coefficient(i));
  }
  vec.push_back(Innermost_coefficient(3));
  assert(Innermost_coefficient(0) 
      == substitute_homogeneous(Polynomial_d(0),vec.begin(),vec.end()));
  assert(Innermost_coefficient(1) 
      == substitute_homogeneous(Polynomial_d(1),vec.begin(),vec.end()));
  assert(Innermost_coefficient(2) 
      == substitute_homogeneous(Polynomial_d(2),vec.begin(),vec.end()));
  assert(Innermost_coefficient(-2) 
      == substitute_homogeneous(Polynomial_d(-2),vec.begin(),vec.end())); 
    
  for(int i = 0; i< 2; i++){
    typedef typename PT
      :: template Rebind<Innermost_coefficient,2>::Other PT_2;
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


// #############

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

  test_derivative(traits);
  test_make_square_free(traits);
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
  test_square_free_factorization_up_to_constant_factor(traits);
  
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
void test_real_embeddable_functors(const PT& traits, CGAL::Tag_false){
  // Since Innermost_coefficient is not RealEmbeddable the following functors
  // should be CGAL::Null_functor.   
  ASSERT_IS_NULL_FUNCTOR(typename PT::Sign_at);
  ASSERT_IS_NULL_FUNCTOR(typename PT::Sign_at_homogeneous);
  ASSERT_IS_NULL_FUNCTOR(typename PT::Compare);
  
}

// test functors depending on the Algebraic_category of ICoeff
template< class PT >
void test_ac_icoeff_functors(
    const PT& traits, CGAL::Integral_domain_without_division_tag){
  ASSERT_IS_NULL_FUNCTOR(typename PT::Multivariate_content);
  ASSERT_IS_NULL_FUNCTOR(typename PT::Interpolate);
}

template< class PT >
void test_ac_icoeff_functors(
    const PT& traits, CGAL::Unique_factorization_domain_tag){
  test_multivariate_content(traits);
  ASSERT_IS_NULL_FUNCTOR(typename PT::Interpolate);
}
template< class PT >
void test_ac_icoeff_functors(const PT& traits, CGAL::Field_tag){
  test_multivariate_content(traits);
  test_interpolate(traits);
}

// test functors depending on the Algebraic_category of Coefficient 
template< class PT >
void test_ac_poly_functors(const PT& traits, CGAL::Integral_domain_without_division_tag){
  ASSERT_IS_NULL_FUNCTOR(typename PT::Univariate_content); 
  ASSERT_IS_NULL_FUNCTOR(typename PT::Square_free_factorization); 
}

template< class PT >
void test_ac_poly_functors(const PT& traits, CGAL::Unique_factorization_domain_tag){
  test_univariate_content(traits);  
  test_square_free_factorization(traits);
}




template< class PT >
void test_polynomial_traits_d(const PT& traits){
  typedef typename PT::Polynomial_d           Polynomial_d;
  typedef typename PT::Innermost_coefficient  ICoeff;
  typedef typename PT::Coefficient            Coeff;

  test_fundamental_functors(traits);

  typedef typename CGAL::Algebraic_structure_traits<ICoeff> AST_IC;
  test_ac_icoeff_functors(traits, typename AST_IC::Algebraic_category());

  typedef typename CGAL::Algebraic_structure_traits<Polynomial_d> AST_Poly;
  test_ac_poly_functors(traits, typename AST_Poly::Algebraic_category());

  typedef typename CGAL::Real_embeddable_traits<ICoeff> RET_IC;
  typedef typename RET_IC::Is_real_embeddable Is_real_embeddable; 
  test_real_embeddable_functors(traits, Is_real_embeddable());
  
}

template< class InnermostCoefficient >
void test_multiple_dimensions() {
  {
    typedef CGAL::Polynomial< InnermostCoefficient > Polynomial_1;
       
    const int dimension                    = 1;
    typedef Polynomial_1                   Polynomial_d;
    typedef InnermostCoefficient           Coefficient;
    typedef InnermostCoefficient           Innermost_coefficient;
        
    typedef CGAL::Polynomial_traits_d<Polynomial_d> PT;
    typedef CGAL::Algebraic_structure_traits<Innermost_coefficient> AST_IC;
    typedef typename AST_IC::Algebraic_category Algebraic_category;
                
    BOOST_STATIC_ASSERT(
        (boost::is_same<typename PT::Polynomial_d,Polynomial_d>::value));
    BOOST_STATIC_ASSERT(
        (boost::is_same< typename PT::Coefficient, Coefficient>::value));
    BOOST_STATIC_ASSERT((boost::is_same< typename PT::Innermost_coefficient, 
            Innermost_coefficient>::value));
    BOOST_STATIC_ASSERT((PT::d == dimension));
    test_polynomial_traits_d(PT());
  }
  {
    typedef CGAL::Polynomial< InnermostCoefficient > Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;

    const int dimension                  = 2;
    typedef Polynomial_2                   Polynomial_d;
    typedef Polynomial_1                   Coefficient;
    typedef InnermostCoefficient           Innermost_coefficient;
        
    typedef CGAL::Polynomial_traits_d<Polynomial_d> PT;
    typedef CGAL::Algebraic_structure_traits<Innermost_coefficient> AST_IC;
    typedef typename AST_IC::Algebraic_category Algebraic_category;
        
    BOOST_STATIC_ASSERT(
        (boost::is_same<typename PT::Polynomial_d,Polynomial_d>::value));
    BOOST_STATIC_ASSERT(
        (boost::is_same< typename PT::Coefficient, Coefficient>::value));
    BOOST_STATIC_ASSERT((boost::is_same< typename PT::Innermost_coefficient, 
            Innermost_coefficient>::value));
    BOOST_STATIC_ASSERT((PT::d == dimension));
    test_polynomial_traits_d(PT());
  }
  {
    typedef CGAL::Polynomial< InnermostCoefficient > Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;

    const int dimension                  = 3;
    typedef Polynomial_3                   Polynomial_d;
    typedef Polynomial_2                   Coefficient;
    typedef InnermostCoefficient           Innermost_coefficient;
  
    typedef CGAL::Algebraic_structure_traits<Innermost_coefficient> AST_IC;
    typedef typename AST_IC::Algebraic_category Algebraic_category;
    typedef CGAL::Polynomial_traits_d<Polynomial_d> PT;
        
        
    BOOST_STATIC_ASSERT(
        (boost::is_same<typename PT::Polynomial_d,Polynomial_d>::value));
    BOOST_STATIC_ASSERT(
        (boost::is_same< typename PT::Coefficient, Coefficient>::value));
    BOOST_STATIC_ASSERT((boost::is_same< typename PT::Innermost_coefficient, 
            Innermost_coefficient>::value));
    BOOST_STATIC_ASSERT((PT::d == dimension));
    test_polynomial_traits_d(PT());
  }   
} 

template < typename AK>
void test_rebind(){
  typedef typename AK::Integer Integer; 
  typedef typename AK::Rational Rational;
  typedef CGAL::Polynomial<Integer> Poly_int_1;                
  typedef CGAL::Polynomial<Poly_int_1> Poly_int_2;              
  typedef CGAL::Polynomial<Poly_int_2> Poly_int_3;              
  typedef CGAL::Polynomial<Rational> Poly_rat_1;               
  typedef CGAL::Polynomial<Poly_rat_1> Poly_rat_2;              
  typedef CGAL::Polynomial<Poly_rat_2> Poly_rat_3;         
    
  typedef CGAL::Polynomial_traits_d<Poly_int_1> PT_int_1;
  typedef CGAL::Polynomial_traits_d<Poly_int_2> PT_int_2;
  typedef CGAL::Polynomial_traits_d<Poly_int_3> PT_int_3;
  typedef CGAL::Polynomial_traits_d<Poly_rat_1> PT_rat_1;
  typedef CGAL::Polynomial_traits_d<Poly_rat_2> PT_rat_2;
  typedef CGAL::Polynomial_traits_d<Poly_rat_3> PT_rat_3;

  typedef typename PT_int_1:: template Rebind<Integer,1>::Other PT_int_1_;
  typedef typename PT_int_3:: template Rebind<Integer,2>::Other PT_int_2_;
  typedef typename PT_rat_3:: template Rebind<Integer,3>::Other PT_int_3_;
  typedef typename PT_int_1:: template Rebind<Rational,1>::Other PT_rat_1_;
  typedef typename PT_rat_2:: template Rebind<Rational,2>::Other PT_rat_2_;
  typedef typename PT_int_2:: template Rebind<Rational,3>::Other PT_rat_3_;
    
  BOOST_STATIC_ASSERT((boost::is_same<PT_int_1_,PT_int_1>::value));
  BOOST_STATIC_ASSERT((boost::is_same<PT_int_2_,PT_int_2>::value));
  BOOST_STATIC_ASSERT((boost::is_same<PT_int_3_,PT_int_3>::value));
  BOOST_STATIC_ASSERT((boost::is_same<PT_rat_1_,PT_rat_1>::value));
  BOOST_STATIC_ASSERT((boost::is_same<PT_rat_2_,PT_rat_2>::value));
  BOOST_STATIC_ASSERT((boost::is_same<PT_rat_3_,PT_rat_3>::value));

  BOOST_STATIC_ASSERT((!boost::is_same<PT_rat_3_,PT_rat_2>::value));
}

template < typename AT> 
void test_AT(){
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);


  std::cerr << std::endl;
  std::cerr << 
    "Test for coefficient type CGAL::Modular" 
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  test_multiple_dimensions< CGAL::Modular >();   


  typedef typename AT::Integer Integer;
  typedef typename AT::Rational Rational; 

  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Integer" << std::endl;
  std::cerr << "--------------------------------------" << std::endl;
  test_multiple_dimensions<Integer>();

  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Rational" << std::endl;
  std::cerr << "---------------------------------------" << std::endl;
  test_multiple_dimensions<Rational>();
    
  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Sqrt_extension< Integer, Integer >" 
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  test_multiple_dimensions< CGAL::Sqrt_extension< Integer, Integer > >();    

  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Sqrt_extension< Rational, Integer >"
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  test_multiple_dimensions< CGAL::Sqrt_extension< Rational, Integer > >();    

  std::cerr << std::endl;
  std::cerr << 
    "Test for coefficient type Sqrt_extension< Rational, Rational >" 
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  test_multiple_dimensions< CGAL::Sqrt_extension< Rational, Rational > >();   



    
  test_rebind<AT>();
}    
  

int main(){

#if 1    
  test_AT<CGAL::Arithmetic_kernel>();
#else
#ifdef CGAL_USE_LEDA
  {        
    typedef CGAL::LEDA_arithmetic_kernel AT;
    test_AT<AT>();
  }
#endif
#ifdef CGAL_USE_CORE
  {    
    typedef CGAL::CORE_arithmetic_kernel AT;
    test_AT<AT>();
  }
#endif
#endif
  return 0;
}
