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
#include <cassert>

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/ipower.h>
#include <CGAL/Random.h>
#include <cmath>

static CGAL::Random my_rnd(346); // some seed 

#define CGAL_SNAP_CGALi_TRAITS_D(T)                             \
  typedef T PT;                                                 \
  typedef typename PT::Polynomial_d          Polynomial_d;      \
  typedef typename PT::Coefficient_type           Coeff;             \
  typedef typename PT::Innermost_coefficient_type ICoeff;            \
  typedef CGAL::Polynomial_traits_d<Coeff> PTC;                 \
  typedef CGAL::Exponent_vector Exponent_vector;                \
  typedef std::pair< CGAL::Exponent_vector , ICoeff > Monom;    


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
  double md = max_degree+1;
  int number_of_coeffs = 
    (CGAL::min)(number_of_variables * (int)std::ceil(std::log(md))+1,4);
    
  Polynomial_d result; 
  for(int i = 0; i < number_of_coeffs; i++){
    CGAL::Exponent_vector exps;
    for(int j = 0; j < PT::d; j++){
      exps.push_back(my_rnd.get_int(0,max_degree));
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
    std::list<Coeff> coeffs;
    assert(construct(coeffs.begin(),coeffs.end()) == construct(0));
    for(int i = 0; i<4;i++){coeffs.push_back(Coeff(i));}
    assert(construct(coeffs.begin(),coeffs.end())
        == construct(Coeff(0),Coeff(1),Coeff(2),Coeff(3)));
        

    typedef std::list< Monom > Monom_list;    
    typedef std::vector< Monom > Monom_vec;    

    // construct from InputIterator
    Monom_list monom_list;
    assert(construct(monom_list.begin(),monom_list.end()) == construct(0));
    CGAL::Random rnd(7);
    for(int j = 0; j < 2; j++){
      CGAL::Exponent_vector exps;
      for(int i = 0; i < d; i++){
        exps.push_back(j+i*5);
      }
      monom_list.push_back(Monom(exps,ICoeff(j+1)));
    };
    
    Monom_vec monom_vec(monom_list.begin(),monom_list.end());
    std::random_shuffle(monom_vec. begin(),monom_vec. end());
      
    Polynomial_d p1 = construct(monom_vec. begin(),
        monom_vec. begin()+((monom_vec. end()- monom_vec. begin())/2));
    Polynomial_d p2 = construct(monom_vec. begin()+
        ((monom_vec. end()- monom_vec. begin())/2), monom_vec. end());
    Polynomial_d p  = construct(monom_vec. begin(), monom_vec. end());
    assert(p == p1+p2);
    
    assert(construct(monom_vec. begin(),monom_vec. end()) 
        == construct(monom_vec.rbegin(),monom_vec.rend()));
    // test with boolean flag is_sorted 
    assert(construct(monom_vec. begin(),monom_vec. end(),false) 
        == construct(monom_vec.rbegin(),monom_vec.rend(),false));  
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
  (void) construct;
  (void) get_coeff;
    
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
  (void) get_innermost_coeff;
  
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
  typedef std::vector< Monom > Monom_rep;    

  typename PT::Construct_polynomial construct;
  typename PT::Get_monom_representation gmr;
  
  {
    Polynomial_d zero(0);
    Monom_rep monom_rep;
    gmr(zero,std::back_inserter(monom_rep));
    assert(monom_rep.size()==1);
    assert(construct(monom_rep.begin(),monom_rep.end()) == zero);
  }
  {
    Polynomial_d x = CGAL::shift(Polynomial_d(1),1,PT::d-1);
    Polynomial_d p = x*x-1;
    Monom_rep monom_rep;
    gmr(p,std::back_inserter(monom_rep));
    assert(monom_rep.size()==2);
    assert(construct(monom_rep.begin(),monom_rep.end()) == p);
  }
    
 
        
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
    
  Polynomial_d zero; 
  assert(swap(zero,0,d-1) == zero);
  
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
void test_permute(const Polynomial_traits_d&){
  std::cerr << "start test_permutate "; std::cerr.flush();

  const int d=4;
  typedef typename Polynomial_traits_d::Innermost_coefficient_type ICoeff;
  typedef typename Polynomial_traits_d:: template Rebind<ICoeff,d>::Other PT4;
  typename PT4::Permute permutate;
  typedef typename PT4::Shift Shift;
  typedef typename PT4::Polynomial_d Polynomial_d;
  typedef typename std::vector<int>::iterator Iterator;
  
  Polynomial_d x  = Shift()(Polynomial_d(1),1,0); 
  Polynomial_d y  = Shift()(Polynomial_d(1),1,1); 
  Polynomial_d z  = Shift()(Polynomial_d(1),1,2); 
  Polynomial_d w1 = Shift()(Polynomial_d(1),1,3); 
     
  std::vector<int> change(d);
  Iterator changeb = change.begin();
  Iterator changee = change.end();
  int lauf = 0;
  for (Iterator itlauf = changeb; itlauf != changee; ++itlauf){
    change[lauf]= lauf;
    lauf++ ;
  }
  Polynomial_d zero (0);
  Polynomial_d one (1); 
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
  typename Polynomial_traits_d::Move move; (void) move;
  typename Polynomial_traits_d::Swap swap; (void) swap;
    
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
  (void) construct;
  (void) degree;
    
  Polynomial_d p; 
  p= construct(Coeff(0));
  assert(degree(p) == 0);
  assert(degree(p,0) == 0);
  p= construct(Coeff(1));
  assert(degree(p) == 0);
  assert(degree(p,0) == 0);
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
    
  typename PT::Innermost_leading_coefficient ilcoeff; (void) ilcoeff;
  Polynomial_d p(Coeff(1), Coeff(2), Coeff(3));
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

// //       Shift;
template <class Polynomial_traits_d>
void test_shift(const Polynomial_traits_d&){
    
  std::cerr << "start test_shift "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Shift shift;  (void) shift;
  typename PT::Swap  swap;  (void) swap;
  typename PT::Construct_polynomial construct;
  for(int i = 0; i < 5; i++){ 
    Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
    Polynomial_d q = p*CGAL::ipower(construct(Coeff(0),Coeff(1)),5);
    assert(shift(p,5) == q);   
  }
  int d = Polynomial_traits_d::d;
  if( d > 1){ 
    assert(shift(Polynomial_d(1),1,0) !=  shift(Polynomial_d(1),1,1));
    assert(shift(Polynomial_d(1),1,0) !=  shift(Polynomial_d(1),1,1));
    assert(shift(Polynomial_d(1),1,0) !=  shift(Polynomial_d(1),1,d-1));
    assert(shift(Polynomial_d(1),1,0) !=  shift(Polynomial_d(1),1,d-1));
    
    assert(shift(Polynomial_d(1),1,0) ==  swap(shift(Polynomial_d(1),1,1),0,1));
    assert(shift(Polynomial_d(1),1,0) ==  swap(shift(Polynomial_d(1),1,1),0,1));
    assert(shift(Polynomial_d(1),1,0) ==  swap(shift(Polynomial_d(1),1,d-1),0,d-1));
    assert(shift(Polynomial_d(1),1,0) ==  swap(shift(Polynomial_d(1),1,d-1),0,d-1));
    
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
  typename PT::Invert invert; 
  typename PT::Swap swap; 
  (void) invert; 
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
    int n; n = my_rnd.get_int(0,PT::d-1);
        
    assert(invert(p,n) == swap(invert(swap(p,n,PT::d-1)),n,PT::d-1));
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
  (void) translate;
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
  typename PT::Translate_homogeneous transh;
  typename PT::Canonicalize canonicalize;
  typename PT::Evaluate_homogeneous evh;
  (void) transh;
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
    
  typename PT::Scale scale;
  (void) scale;
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
  (void) scaleh;
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
void test_derivative(const Polynomial_traits_d&){
  std::cerr << "start test_derivative "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
    
  typename PT::Differentiate diff;
  typename PT::Swap swap;
  (void) diff;
  (void) swap;
    
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


// //       Square_free_factorize;
template <class Polynomial_traits_d>
void test_square_free_factorize(const Polynomial_traits_d&){
  std::cerr << "start test_square_free_factorize "; std::cerr.flush(); 
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typedef CGAL::Algebraic_structure_traits<Polynomial_d> AST;
  typename AST::Integral_division idiv;
  typename PT::Square_free_factorize sqff;
  typename PT::Canonicalize canonicalize;
  typename PT::Total_degree total_degree;
  typename PT::Leading_coefficient lcoeff;
  (void) idiv;
  (void) sqff;
  (void) canonicalize;

  for(int i = 0; i < 5; i++){
    Polynomial_d f1 = generate_sparse_random_polynomial<Polynomial_d>(2);
    Polynomial_d f2 = generate_sparse_random_polynomial<Polynomial_d>(2);  
    Polynomial_d f3 = generate_sparse_random_polynomial<Polynomial_d>(2);    
    f3 = Polynomial_d(lcoeff(f3));
    Polynomial_d p = f1*f1*f2*f3*f3*f3;
    std::vector<std::pair<Polynomial_d,int> > fac_mul_pairs;
    sqff(p, std::back_inserter(fac_mul_pairs));
    int n = fac_mul_pairs.size();
    assert(n >= 3 
        || total_degree(f1) == 0 || total_degree(f2) == 0 || total_degree(f3) == 0);
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
    if (!CGAL::is_zero(q)){
      pdiv(f,g,q,r,D);
      assert(f*Polynomial_d(D) == g*q+r);
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
  typename PT::Pseudo_division_remainder pdiv_r;
  (void) pdiv;
  (void) pdiv_r;

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
  typename PT::Pseudo_division_quotient pdiv_q;
  (void) pdiv_q;
  (void) pdiv;
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
  typename PT::Gcd_up_to_constant_factor gcd_utcf;
  (void) gcd_utcf;
    
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
  (void) idiv_utcf;
  (void) canonicalize;
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
  
  typename PT::Univariate_content_up_to_constant_factor ucontent_utcf;
  typename PT::Integral_division_up_to_constant_factor idiv_utcf;
  typename PT::Leading_coefficient lcoeff;
  typename PT::Canonicalize canonicalize; 
  
  (void) ucontent_utcf;
  (void) idiv_utcf;
  (void) lcoeff;
  (void) canonicalize;

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


// //       Square_free_factorize_up_to_constant_factor;
template <class Polynomial_traits_d>
void test_square_free_factorize_up_to_constant_factor(const Polynomial_traits_d&){
  std::cerr << "start test_square_free_factorize_up_to_constant_factor "; 
  std::cerr.flush(); 
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Integral_division_up_to_constant_factor idiv_utcf;
  typename PT::Square_free_factorize_up_to_constant_factor sqff_utcf;
  typename PT::Canonicalize canonicalize;
  
  (void) idiv_utcf;
  (void) sqff_utcf;
  (void) canonicalize;

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
  (void) evaluate;
  (void) move;
  assert(evaluate(Polynomial_d(0),Coeff(0)) == Coeff(0));
  assert(evaluate(Polynomial_d(1),Coeff(0)) == Coeff(1));
  assert(evaluate(Polynomial_d(2),Coeff(5)) == Coeff(2));

  assert( evaluate(Polynomial_d(Coeff(3),Coeff(2)),Coeff(0)) == Coeff(3));
  assert( evaluate(Polynomial_d(Coeff(3),Coeff(2)),Coeff(1)) == Coeff(5));
  assert( evaluate(Polynomial_d(Coeff(3),Coeff(2)),Coeff(2)) == Coeff(7));
      
  std::cerr << " ok "<< std::endl;
}

// //       Evaluate_homogeneous;
template <class Polynomial_traits_d>
void test_evaluate_homogeneous(const Polynomial_traits_d&){
  std::cerr << "start test_evaluate_homogeneous "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Evaluate_homogeneous evh;
  (void) evh;
  
  assert(evh(Polynomial_d(0),Coeff(0),Coeff(1)) == Coeff(0));
  assert(evh(Polynomial_d(1),Coeff(0),Coeff(2)) == Coeff(1));
  assert(evh(Polynomial_d(2),Coeff(5),Coeff(3)) == Coeff(2));

  assert(evh( Polynomial_d(Coeff(3),Coeff(2)) , Coeff(0),Coeff(1)) 
      == Coeff(3));
  assert(evh( Polynomial_d(Coeff(3),Coeff(2)) , Coeff(1),Coeff(1)) 
      == Coeff(5));
  assert(evh( Polynomial_d(Coeff(3),Coeff(2)) , Coeff(2),Coeff(3)) 
      == Coeff(9+4));
    
  std::cerr << " ok "<< std::endl;
}

template< class Polynomial_traits_d >
void test_is_zero_at(const Polynomial_traits_d&) {
  std::cerr << "start test_is_zero_at "; 
  std::cerr.flush();
    
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
    
  typename PT::Is_zero_at is_zero_at;
  (void) is_zero_at;
    
  Polynomial_d p(Coeff(-1), Coeff(0), Coeff(1));
    
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
    
  typename PT::Is_zero_at_homogeneous is_zero_at_homogeneous;
  (void) is_zero_at_homogeneous;
    
  Polynomial_d p(Coeff(-1), Coeff(0), Coeff(4));
    
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
    
  typename PT::Sign_at sign_at;
  (void) sign_at;
    
  Polynomial_d p(Coeff(-1), Coeff(0), Coeff(1));
    
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
    
  typename PT::Sign_at_homogeneous sign_at_homogeneous;
  (void) sign_at_homogeneous;
    
  Polynomial_d p(Coeff(-1), Coeff(0), Coeff(1));
    
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
    
  typename PT::Compare compare; (void) compare;
    
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
  (void) resultant;
  (void) move;
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
  std::cerr << " ok "<< std::endl;    
}

// //       Canonicalize;
template <class Polynomial_traits_d>
void test_canonicalize(const Polynomial_traits_d&){
  std::cerr << "start test_canonicalize "; std::cerr.flush();
  CGAL_SNAP_CGALi_TRAITS_D(Polynomial_traits_d);
  typename PT::Canonicalize canonicalize;
  (void) canonicalize;

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
  typename PT::Substitute substitute;  (void) substitute;
  typedef typename PT::Innermost_coefficient_type Innermost_coefficient_type;
    
    
  std::list<Innermost_coefficient_type> list;
  for(int i = 0; i < PT::d; i++){
    list.push_back(Innermost_coefficient_type(i));
  }
  assert(Innermost_coefficient_type(0) 
      == substitute(Polynomial_d(0),list.begin(),list.end()));
  assert(Innermost_coefficient_type(1) 
      == substitute(Polynomial_d(1),list.begin(),list.end()));
  assert(Innermost_coefficient_type(2) 
      == substitute(Polynomial_d(2),list.begin(),list.end()));
  assert(Innermost_coefficient_type(-2) 
      == substitute(Polynomial_d(-2),list.begin(),list.end())); 
    
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
  typename PT::Substitute_homogeneous substitute_homogeneous;
  (void) substitute_homogeneous;
  typedef typename PT::Innermost_coefficient_type Innermost_coefficient_type;
    
    
  std::list<Innermost_coefficient_type> list;
  for(int i = 0; i < PT::d; i++){
    list.push_back(Innermost_coefficient_type(i));
  }
  list.push_back(Innermost_coefficient_type(3));
  assert(Innermost_coefficient_type(0) 
      == substitute_homogeneous(Polynomial_d(0),list.begin(),list.end()));
  assert(Innermost_coefficient_type(1) 
      == substitute_homogeneous(Polynomial_d(1),list.begin(),list.end()));
  assert(Innermost_coefficient_type(2) 
      == substitute_homogeneous(Polynomial_d(2),list.begin(),list.end()));
  assert(Innermost_coefficient_type(-2) 
      == substitute_homogeneous(Polynomial_d(-2),list.begin(),list.end())); 
    
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

template <class PT >
void test_coefficient_const_iterator(const PT&) {

  typedef typename PT::Polynomial_d                Polynomial_d;
  typedef typename PT::Coefficient_type            Coefficient; 
  typedef typename PT::Coefficient_const_iterator  CCIterator;
  
  typename PT::Construct_coefficient_const_iterator_range coeff_range;
  typename PT::Degree                  degree;
  typename PT::Get_coefficient         coeff;
  
  Polynomial_d p = generate_sparse_random_polynomial<Polynomial_d>();
  
  
  CCIterator it = coeff_range(p).first;
  for(int i = 0; i <= degree(p); i++){
    assert(*it == coeff(p,i));
    it++;
  }
  assert(coeff_range(p).second == it);
}



template <class PT>
void test_innermost_coefficient_const_iterator(const PT&) {
  
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


  typename PT_1::Innermost_coefficient_const_iterator it1; (void) it1;
  typename PT_1::Construct_innermost_coefficient_const_iterator_range range1; 
  typename PT_2::Innermost_coefficient_const_iterator it2; (void) it2;
  typename PT_2::Construct_innermost_coefficient_const_iterator_range range2; 
  typename PT_3::Innermost_coefficient_const_iterator it3; (void) it3;
  typename PT_3::Construct_innermost_coefficient_const_iterator_range range3;
  
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
}


template< class PT >
void test_polynomial_traits_d(const PT& traits){
  typedef typename PT::Polynomial_d           Polynomial_d;
  typedef typename PT::Innermost_coefficient_type  ICoeff;
  typedef typename PT::Coefficient_type            Coeff;

  test_fundamental_functors(traits);

  typedef typename CGAL::Algebraic_structure_traits<ICoeff> AST_IC;
  test_ac_icoeff_functors(traits, typename AST_IC::Algebraic_category());

  typedef typename CGAL::Algebraic_structure_traits<Polynomial_d> AST_Poly;
  test_ac_poly_functors(traits, typename AST_Poly::Algebraic_category());

  typedef typename CGAL::Real_embeddable_traits<ICoeff> RET_IC;
  typedef typename RET_IC::Is_real_embeddable Is_real_embeddable; 
  test_real_embeddable_functors(traits, Is_real_embeddable()); 
  test_coefficient_const_iterator(traits);
  test_innermost_coefficient_const_iterator(traits);
}

template< class InnermostCoefficient_type >
void test_multiple_dimensions() {
  {
    typedef CGAL::Polynomial< InnermostCoefficient_type > Polynomial_1;
       
    const int dimension                    = 1;
    typedef Polynomial_1                   Polynomial_d;
    typedef InnermostCoefficient_type           Coefficient_type;
    typedef InnermostCoefficient_type           Innermost_coefficient_type;
        
    typedef CGAL::Polynomial_traits_d<Polynomial_d> PT;
    typedef CGAL::Algebraic_structure_traits<Innermost_coefficient_type> AST_IC;
    typedef typename AST_IC::Algebraic_category Algebraic_category;
                
    BOOST_STATIC_ASSERT(
        (boost::is_same<typename PT::Polynomial_d,Polynomial_d>::value));
    BOOST_STATIC_ASSERT(
        (boost::is_same< typename PT::Coefficient_type, Coefficient_type>::value));
    BOOST_STATIC_ASSERT((boost::is_same< typename PT::Innermost_coefficient_type, 
            Innermost_coefficient_type>::value));
    BOOST_STATIC_ASSERT((PT::d == dimension));
    test_polynomial_traits_d(PT());
  }
  {
    typedef CGAL::Polynomial< InnermostCoefficient_type > Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;

    const int dimension                  = 2;
    typedef Polynomial_2                   Polynomial_d;
    typedef Polynomial_1                   Coefficient_type;
    typedef InnermostCoefficient_type           Innermost_coefficient_type;
        
    typedef CGAL::Polynomial_traits_d<Polynomial_d> PT;
    typedef CGAL::Algebraic_structure_traits<Innermost_coefficient_type> AST_IC;
    typedef typename AST_IC::Algebraic_category Algebraic_category;
        
    BOOST_STATIC_ASSERT(
        (boost::is_same<typename PT::Polynomial_d,Polynomial_d>::value));
    BOOST_STATIC_ASSERT(
        (boost::is_same< typename PT::Coefficient_type, Coefficient_type>::value));
    BOOST_STATIC_ASSERT((boost::is_same< typename PT::Innermost_coefficient_type, 
            Innermost_coefficient_type>::value));
    BOOST_STATIC_ASSERT((PT::d == dimension));
    test_polynomial_traits_d(PT());
  }
  {
    typedef CGAL::Polynomial< InnermostCoefficient_type > Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;

    const int dimension                  = 3;
    typedef Polynomial_3                   Polynomial_d;
    typedef Polynomial_2                   Coefficient_type;
    typedef InnermostCoefficient_type           Innermost_coefficient_type;
  
    typedef CGAL::Algebraic_structure_traits<Innermost_coefficient_type> AST_IC;
    typedef typename AST_IC::Algebraic_category Algebraic_category;
    typedef CGAL::Polynomial_traits_d<Polynomial_d> PT;
        
        
    BOOST_STATIC_ASSERT(
        (boost::is_same<typename PT::Polynomial_d,Polynomial_d>::value));
    BOOST_STATIC_ASSERT(
        (boost::is_same< typename PT::Coefficient_type, Coefficient_type>::value));
    BOOST_STATIC_ASSERT((boost::is_same< typename PT::Innermost_coefficient_type, 
            Innermost_coefficient_type>::value));
    BOOST_STATIC_ASSERT((PT::d == dimension));
    test_polynomial_traits_d(PT());
  }   
  {
  	typedef CGAL::Polynomial< InnermostCoefficient_type > Polynomial_1;
  	typedef CGAL::Polynomial_traits_d<Polynomial_1> PT; 
  	test_permute(PT());
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

  std::cerr << std::endl;
  std::cerr << 
    "Test for coefficient type CGAL::Residue" 
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  // Enforce IEEE double precision before using modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);
  test_multiple_dimensions< CGAL::Residue >();   
}    
  

int main(){

    // Set wrong rounding mode to test modular arithmetic 
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);

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

  return 0;
}
