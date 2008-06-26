// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber  <mkerber@mpi-inf.mpg.de>
//                 Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef AcX_DEBUG_PRINT
#define AcX_DEBUG_PRINT 1
#endif

#if AcX_DEBUG_PRINT
#define AcX_DSTREAM(str) std::cout << str;
#else
#define AcX_DSTREAM(str) 
#endif


#ifndef BITSTREAM_USES_E08_TREE
#define BITSTREAM_USES_E08_TREE 0
#endif

#include <CGAL/basic.h>

#include <sstream>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_pure.h>

#include <CGAL/Algebraic_curve_kernel_2/Bitstream_descartes_at_x/Bitstream_descartes_bfs.h>

#include <CGAL/Algebraic_curve_kernel_2/Bitstream_descartes_at_x/Bitstream_descartes_traits_on_vert_line.h>

// Traits class for Integer coefficients
template <class ArithmeticKernel>
class Bitstream_descartes_rndl_tree_traits_from_Integer_coeff {

public:
  typedef ArithmeticKernel Arithmetic_kernel;

  // input coefficients of polynomial
  typedef typename Arithmetic_kernel::Integer Coefficient; 
  
  // type for internal computations
  typedef typename Arithmetic_kernel::Integer Integer; 

  // interval bdrys
  typedef typename Arithmetic_kernel::Rational Boundary; 


  class Boundary_creator {
      
  public:
      
    Boundary_creator() {}

    Boundary operator() (Integer x,long p) {
      Integer num=x, denom;
      if(p <0) {
          denom = CGAL::ipower(Integer(2),-p);
      }
      else {
          num*=CGAL::ipower(Integer(2),p);
          denom=1;
      }
      return Boundary(num,denom);
    } 
  };

  Bitstream_descartes_rndl_tree_traits_from_Integer_coeff(int state = 0)
  { }

  // Integer approximation to 2^p*x with abs error <= 1  (not(!) 1/2)
  class Approximator {
  public:
    Approximator() {
    }
    Integer operator() (Coefficient x, long p) {
      if      (p >  0)  return x << p;
      else if (p <  0)  return x >> -p;
      else   /*p == 0*/ return x;
    }
  };
  Approximator approximator_object() const { return Approximator(); }

  // We can just take this from NT_traits, but need to supply ..._object()
    typedef typename CGAL::CGALi::Real_embeddable_extension<Coefficient>
    ::Floor_log2_abs Lower_bound_log2_abs;
  Lower_bound_log2_abs lower_bound_log2_abs_object() const {
    return Lower_bound_log2_abs();
  }

  // For the upper bound approximator, we just the NT_traits again
  class Upper_bound_log2_abs_approximator {

  public:

    bool initial_upper_bound(Coefficient ai,
			     long& ub_log2,
			     bool& is_certainly_zero) {
      return improve_upper_bound(ai,ub_log2,is_certainly_zero);
    }
    
    bool improve_upper_bound(Coefficient ai,
			     long& ub_log2,
			     bool& is_certainly_zero) {
      if(ai==Coefficient(0)) {
	is_certainly_zero=true;
      }
      else {
	is_certainly_zero=false;
	ub_log2=CGAL::CGALi::ceil_log2_abs(ai);
      }
      return true;
    }

  };
  
  Upper_bound_log2_abs_approximator 
  upper_bound_log2_abs_approximator_object() {
    return Upper_bound_log2_abs_approximator();
  }

  // We can just take these from NT_traits
  typedef typename CGAL::Real_embeddable_traits<Integer>::Sign Sign;
    typedef typename CGAL::CGALi::Real_embeddable_extension<Integer>
    ::Ceil_log2_abs Ceil_log2_abs_Integer;
  typedef typename CGAL::CGALi::Real_embeddable_extension<Integer>
    ::Ceil_log2_abs Ceil_log2_abs_long;
};


// A simple model of the EventRefinement concept: 
// Uses a vector of Algebraic reals
template <typename AlgReal>
class Event_refinement {

  typedef AlgReal Alg_real;

  typedef typename AlgReal::Rational Rational;

private:
  std::vector<Alg_real> events;

public:
  Event_refinement() {}
  void add(Alg_real a) {events.push_back(a);}
  Rational lower_boundary(int i){return events[i].low();}
  Rational upper_boundary(int i){return events[i].high();}
  void refine(int i) {events[i].refine();}
};

template<typename ArithmeticKernel>
void test_routine() {
  
  typedef ArithmeticKernel Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Integer Integer;
  typedef typename Arithmetic_kernel::Rational Rational;
  
  typedef CGAL::Polynomial<Integer> Poly_int1;
  typedef CGAL::Polynomial<Poly_int1> Poly_int2;
  typedef CGAL::CGALi::Algebraic_real_pure<Integer,Rational>
      Algebraic_real;

  typedef 
    Bitstream_descartes_rndl_tree_traits_from_Integer_coeff<Arithmetic_kernel>
    Traits;
  
  typedef CGAL::CGALi::Bitstream_descartes_bfs<Traits> Bitstream_descartes;

  Traits traits;

  

  { // Test for the classical Descartes method
    std::stringstream ss("P[15(0,59738427711)(1,300038251861)(2,-471253844514)(3,538575935875)(4,22286946912)(5,548111721525)(6,-379185895352)(7,296681325489)(8,-464256140044)(9,410194800463)(10,232977578849)(11,-376080509486)(12,521721411895)(13,-100316773723)(14,-171187873598)(15,-189253202432)]");
    Poly_int1 f;
    ss >> f;

    // In Maple: f := 59738427711-189253202432*x^15-171187873598*x^14-100316773723*x^13+521721411895*x^12-376080509486*x^11+232977578849*x^10+410194800463*x^9-464256140044*x^8+296681325489*x^7-379185895352*x^6+548111721525*x^5+22286946912*x^4+538575935875*x^3-471253844514*x^2+300038251861*x

    // We expect 3 roots, at -1.176, -.154 and 1.168
    
    CGAL::CGALi::Square_free_descartes_tag t;

    Bitstream_descartes descartes(t, f);

    assert(descartes.degree_of_gcd() == 0);
    assert(descartes.square_free_part() == f);
    
    int n = descartes.number_of_real_roots();
    
    assert(n==3);

    for(int i = 0; i<n;i++) {
      assert(descartes.is_certainly_simple_root(i));
      assert(! descartes.is_certainly_multiple_root(i));
    }

    for(int i =0; i<n; i++) {
      while((descartes.right_boundary(i)-descartes.left_boundary(i))
	    > Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }
    assert(Rational(-118,100)<descartes.left_boundary(0));
    assert(Rational(-117,100)>descartes.right_boundary(0));
    assert(Rational(-155,1000)<descartes.left_boundary(1));
    assert(Rational(-154,1000)>descartes.right_boundary(1));
    assert(Rational(116,100)<descartes.left_boundary(2));
    assert(Rational(117,100)>descartes.right_boundary(2));
  }

  { // Test for the square free case, created with tree from outside

    typedef typename Bitstream_descartes::Bitstream_tree Bitstream_tree;

    std::stringstream ss("P[13(1,1)(4,1)(6,95)(8,77)(11,-51)(13,31)]");

    Poly_int1 f;

    ss >> f;
#if BITSTREAM_USES_E08_TREE
    Bitstream_tree tree(1,f.begin(),f.end(),
                        typename Bitstream_tree::Monomial_basis_tag());
#else
    Bitstream_tree tree(-2,2,0,f.begin(),f.end(),
                        typename Bitstream_tree::Monomial_basis_tag());
#endif

    // We expect 3 roots, at -1.597, -.388 and 0
    
    CGAL::CGALi::Square_free_descartes_tag t;

    Bitstream_descartes descartes(t, f, tree);
    
    assert(descartes.degree_of_gcd() == 0);
    assert(descartes.square_free_part() == f);

    int n = descartes.number_of_real_roots();

    assert(n==3);

    for(int i = 0; i<n;i++) {
      assert(descartes.is_certainly_simple_root(i));
      assert(! descartes.is_certainly_multiple_root(i));
    }

    for(int i =0; i<n; i++) {
      while((descartes.right_boundary(i)-descartes.left_boundary(i))
	    > Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }
    assert(Rational(-160,100)<descartes.left_boundary(0));
    assert(Rational(-159,100)>descartes.right_boundary(0));
    assert(Rational(-389,1000)<descartes.left_boundary(1));
    assert(Rational(-387,1000)>descartes.right_boundary(1));
    assert(Rational(-1,100)<descartes.left_boundary(2));
    assert(Rational(1,100)>descartes.right_boundary(2));
  }
    


  { // Test for the m-k-Descartes method

    // Polnomial with one multiple root:
    
    std::stringstream ss("P[11(0,-51960)(1,158454)(2,3015726)(3,-22833405)(4,64277882)(5,-86502719)(6,58622397)(7,-260172)(8,-77833332)(9,85923423)(10,-37885401)(11,1874259)]");
    // This is a polynomial with a 4-fold root at 1/3, and 3 more roots at
    // approximate positions -1.036, -0.104, 17.764

    // In Maple:g := -51960-86502719*x^5+158454*x-22833405*x^3+3015726*x^2+64277882*x^4+1874259*x^11-77833332*x^8+85923423*x^9-37885401*x^10+58622397*x^6-260172*x^7 

    Poly_int1 f;
    ss >> f;
    
    CGAL::CGALi::M_k_descartes_tag t;
    
    // We expect 4 real roots (m), and the gcd of g and g' is 3 (k)
    Bitstream_descartes descartes(t, f, 4, 3);

    assert(descartes.degree_of_gcd() == 3);
    try {
      Poly_int1 sq_f = descartes.square_free_part();
      assert(false);
    } catch (CGAL::CGALi::Virtual_method_exception ex) {
      // Expected
    }

    int n = descartes.number_of_real_roots();
    
    assert(n==4);

    for(int i = 0; i<n;i++) {
      if(i!=2) {
	assert(descartes.is_certainly_simple_root(i));
	assert(! descartes.is_certainly_multiple_root(i));
      }
      else{ 
	assert(! descartes.is_certainly_simple_root(i));
      }
    }
    for(int i =0; i<n; i++) {
      while(descartes.right_boundary(i)-descartes.left_boundary(i)
	    >Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }
    assert(Rational(-104,100)<descartes.left_boundary(0));
    assert(Rational(-103,100)>descartes.right_boundary(0));
    assert(Rational(-106,1000)<descartes.left_boundary(1));
    assert(Rational(-104,1000)>descartes.right_boundary(1));
    assert(Rational(33,100)<descartes.left_boundary(2));
    assert(Rational(34,100)>descartes.right_boundary(2));
    assert(Rational(1776,100)<descartes.left_boundary(3));
    assert(Rational(1777,100)>descartes.right_boundary(3));
    
  }

  { // Another test for the m-k-method

    // Test for the m-k-Descartes method

    // Polnomial with one multiple root:
    
    // In Maple:f := y^3 - y^2  -2*y 
    Poly_int1 f(0,-2,-1,1);
    
    CGAL::CGALi::M_k_descartes_tag t;
    
    // We expect 3 real roots (m), and the degree of gcd of g and g' is 0 (k)
    Bitstream_descartes descartes(t, f, 3, 0);

    assert(descartes.degree_of_gcd() == 0);
    try {
      Poly_int1 sq_f = descartes.square_free_part();
      assert(false);
    } catch (CGAL::CGALi::Virtual_method_exception ex) {
      // Expected
    }

    int n = descartes.number_of_real_roots();
    
    assert(n==3);

    for(int i = 0; i<n;i++) {
      assert(! descartes.is_certainly_multiple_root(i));
    }
    for(int i =0; i<n; i++) {
      while(descartes.right_boundary(i)-descartes.left_boundary(i)
	    >Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }
    assert(Rational(-101,100)<descartes.left_boundary(0));
    assert(Rational(-99,100)>descartes.right_boundary(0));
    assert(Rational(-1,1000)<descartes.left_boundary(1));
    assert(Rational(1,1000)>descartes.right_boundary(1));
    assert(Rational(199,100)<descartes.left_boundary(2));
    assert(Rational(201,100)>descartes.right_boundary(2));
    
  }
  
  /*
  { // Test for the m-ki-Descartes method

    // Polnomial with one multiple root:
    
    std::stringstream ss("P[11(0,-51960)(1,158454)(2,3015726)(3,-22833405)(4,64277882)(5,-86502719)(6,58622397)(7,-260172)(8,-77833332)(9,85923423)(10,-37885401)(11,1874259)]");
    // This is a polynomial with a 4-fold root at 1/3, and 3 more roots at
    // approximate positions -1.036, -0.104, 17.764

    // In Maple:g := -51960-86502719*x^5+158454*x-22833405*x^3+3015726*x^2+64277882*x^4+1874259*x^11-77833332*x^8+85923423*x^9-37885401*x^10+58622397*x^6-260172*x^7 

    Poly_int1 f;
    ss >> f;
    
    CGAL::CGALi::M_ki_descartes_tag t;
    
    // We expect 4 real roots (m), and the gcd of g and g' is 3 (k)
    Bitstream_descartes descartes(
            t, f, CGAL::CGALi::div_utcf(f, CGAL::CGALi::gcd_utcf(f, CGAL::diff(f)))
    );
    
    int n = descartes.number_of_real_roots();
    
    //std::cout << "num: " << n << std::endl;

    assert(n == 4);

    for(int i = 0; i<n;i++) {
      if(i!=2) {
	assert(descartes.is_certainly_simple_root(i));
	assert(! descartes.is_certainly_multiple_root(i));
        //std::cout << i << "s: " << descartes.multiplicity_of_root(i) << std::endl;
        assert(descartes.multiplicity_of_root(i) == 1);
      }
      else{ 
	assert(! descartes.is_certainly_simple_root(i));
        //std::cout << i << "m: " << descartes.multiplicity_of_root(i) << std::endl;
        assert(descartes.multiplicity_of_root(i) == 4);
      }
    }
    for(int i =0; i<n; i++) {
      while(descartes.right_boundary(i)-descartes.left_boundary(i)
	    >Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }
    assert(Rational(-104,100)<descartes.left_boundary(0));
    assert(Rational(-103,100)>descartes.right_boundary(0));
    assert(Rational(-106,1000)<descartes.left_boundary(1));
    assert(Rational(-104,1000)>descartes.right_boundary(1));
    assert(Rational(33,100)<descartes.left_boundary(2));
    assert(Rational(34,100)>descartes.right_boundary(2));
    assert(Rational(1776,100)<descartes.left_boundary(3));
    assert(Rational(1777,100)>descartes.right_boundary(3));
    
  }
  
  //std::cout << std::endl;

  { // Another test for the m-ki-method

    // In Maple:f := y^3 - y^2  -2*y 
    Poly_int1 f(0,-2,-1,1);
    
    CGAL::CGALi::M_ki_descartes_tag t;
    
    // We expect 3 real roots (m), and the degree of gcd of g and g' is 0 (k)
    Bitstream_descartes descartes(
            t, f, CGAL::CGALi::div_utcf(f, CGAL::CGALi::gcd_utcf(f, CGAL::diff(f)))
    );

    int n = descartes.number_of_real_roots();
    
    //std::cout << "num: " << n << std::endl;

    assert(n==3);

    for(int i = 0; i<n;i++) {
        assert(! descartes.is_certainly_multiple_root(i));
        //std::cout << i << "s: " << descartes.multiplicity_of_root(i) << std::endl;
        assert(descartes.multiplicity_of_root(i) == 1);
    }
    for(int i =0; i<n; i++) {
      while(descartes.right_boundary(i)-descartes.left_boundary(i)
	    >Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }
    assert(Rational(-101,100)<descartes.left_boundary(0));
    assert(Rational(-99,100)>descartes.right_boundary(0));
    assert(Rational(-1,1000)<descartes.left_boundary(1));
    assert(Rational(1,1000)>descartes.right_boundary(1));
    assert(Rational(199,100)<descartes.left_boundary(2));
    assert(Rational(201,100)>descartes.right_boundary(2));
    
  }
  
  //std::cout << std::endl;
  
  { // Test for the m-ki-Descartes method
      
    // Polnomial with three multiple (1(3), 2(2), -4(4)) and two 
    // (10 and -8) single roots:
    
      Poly_int1 fa(Integer(1), Integer(-1));
      Poly_int1 fb(Integer(2), Integer(-1));
      Poly_int1 fc(Integer(-4), Integer(-1));
      Poly_int1 fd(Integer(10), Integer(-1));
      Poly_int1 fe(Integer(-8), Integer(-1));
      
      Poly_int1 f = fa*fa*fa * fb*fb * fc*fc*fc*fc * fd * fe;
    
      CGAL::CGALi::M_ki_descartes_tag t;
      
      Bitstream_descartes descartes(
              t, f, CGAL::CGALi::div_utcf(f, CGAL::CGALi::gcd_utcf(f, CGAL::diff(f)))
      );
      
      int n = descartes.number_of_real_roots();
      
      //std::cout << "num: " << n << std::endl;
      
      assert(n == 5);

      for(int i = 0; i<n;i++) {
          if(i == 0 || i == 4) {
              assert(descartes.is_certainly_simple_root(i));
              //std::cout << i << "s: " << descartes.multiplicity_of_root(i) << std::endl;
              assert(! descartes.is_certainly_multiple_root(i));
              assert(descartes.multiplicity_of_root(i) == 1);
          }
          else{ 
              assert(! descartes.is_certainly_simple_root(i));
              //std::cout << i << "m: " << descartes.multiplicity_of_root(i) << std::endl;
              if (i == 1) {
                  assert(descartes.multiplicity_of_root(i) == 4);
              }
              if (i == 2) {
                  assert(descartes.multiplicity_of_root(i) == 3);
              }
              if (i == 3) {
                  assert(descartes.multiplicity_of_root(i) == 2);
              }
          }
      }

      for(int i =0; i<n; i++) {
          while(descartes.right_boundary(i)-descartes.left_boundary(i)
                >Rational(1,10000)) 
          {
              descartes.refine_interval(i);
          }
      }

      assert(Rational(-85,10)<descartes.left_boundary(0));
      assert(Rational(-75,10)>descartes.right_boundary(0));
      assert(Rational(-41,10)<descartes.left_boundary(1));
      assert(Rational(-39,10)>descartes.right_boundary(1));
      assert(Rational(1,10)<descartes.left_boundary(2));
      assert(Rational(11,10)>descartes.right_boundary(2));
      assert(Rational(19,100)<descartes.left_boundary(3));
      assert(Rational(21,10)>descartes.right_boundary(3));
      assert(Rational(9990,1000)<descartes.left_boundary(4));
      assert(Rational(10001,1000)>descartes.right_boundary(4));
  }

  //std::cout << std::endl;

  { // Test for the m-ki-Descartes method
      
    // Polnomial with three multiple (1(3), 2(2), -4(4)) and two 
    // (10 and -8) single roots:

      Poly_int1 f0(Integer(-20), Integer(-1));
      Poly_int1 f1(Integer(-8), Integer(-1));
      Poly_int1 f2(Integer(-4), Integer(-1));
      Poly_int1 f3(Integer(1), Integer(-1));
      Poly_int1 f4(Integer(2), Integer(-1));
      Poly_int1 f5(Integer(10), Integer(-1));
      Poly_int1 f6(Integer(103), Integer(-1));
      Poly_int1 f7(Integer(140), Integer(-1));
      
      Poly_int1 f = f0*f0*f0 * f1*f1*f1 * f2*f2 * f3*f3*f3*f3 * f4 * f5 * 
          f6*f6*f6*f6*f6*f6 * f7;
      
      CGAL::CGALi::M_ki_descartes_tag t;
      
      Bitstream_descartes descartes(
              t, f, CGAL::CGALi::div_utcf(f, CGAL::CGALi::gcd_utcf(f, CGAL::diff(f)))
      );
      
      int n = descartes.number_of_real_roots();
      
      //std::cout << "num: " << n << std::endl;
      
      assert(n == 8);

      for(int i = 0; i<n;i++) {
          switch (i) {
          case 4: case 5: case 7:
              assert(descartes.is_certainly_simple_root(i));
              //std::cout << i << "s: " << descartes.multiplicity_of_root(i) << std::endl;
              assert(! descartes.is_certainly_multiple_root(i));
              assert(descartes.multiplicity_of_root(i) == 1);
              break;
          case 0:
              assert(! descartes.is_certainly_simple_root(i));
              assert( descartes.is_certainly_multiple_root(i));
              assert(descartes.multiplicity_of_root(i) == 3);
              //std::cout << i << "s: " << descartes.multiplicity_of_root(i) << std::endl;
              break;
          case 1:
              assert(! descartes.is_certainly_simple_root(i));
              assert( descartes.is_certainly_multiple_root(i));
              assert(descartes.multiplicity_of_root(i) == 3);
              //std::cout << i << "s: " << descartes.multiplicity_of_root(i) << std::endl;
              break;
          case 2:
              assert(! descartes.is_certainly_simple_root(i));
              assert( descartes.is_certainly_multiple_root(i));
              assert(descartes.multiplicity_of_root(i) == 2);
              //std::cout << i << "s: " << descartes.multiplicity_of_root(i) << std::endl;
              break;
          case 3:
              assert(! descartes.is_certainly_simple_root(i));
              assert( descartes.is_certainly_multiple_root(i));
              assert(descartes.multiplicity_of_root(i) == 4);
              //std::cout << i << "s: " << descartes.multiplicity_of_root(i) << std::endl;
              break;
          case 6:
              assert(! descartes.is_certainly_simple_root(i));
              assert( descartes.is_certainly_multiple_root(i));
              assert(descartes.multiplicity_of_root(i) == 6);
              //std::cout << i << "s: " << descartes.multiplicity_of_root(i) << std::endl;
              break;
          default:
              assert(false);
          }
      }
      
      for(int i =0; i<n; i++) {
          while(descartes.right_boundary(i)-descartes.left_boundary(i)
                >Rational(1,10000)) 
          {
              descartes.refine_interval(i);
          }
      }

      assert(Rational(-2002,100)<descartes.left_boundary(0));
      assert(Rational(-1997,100)>descartes.right_boundary(0));
      assert(Rational(-85,10)<descartes.left_boundary(1));
      assert(Rational(-75,10)>descartes.right_boundary(1));
      assert(Rational(-41,10)<descartes.left_boundary(2));
      assert(Rational(-39,10)>descartes.right_boundary(2));
      assert(Rational(950,1000)<descartes.left_boundary(3));
      assert(Rational(1001,1000)>descartes.right_boundary(3));
      assert(Rational(1999,1000)<descartes.left_boundary(4));
      assert(Rational(2001,1000)>descartes.right_boundary(4));
      assert(Rational(9999,1000)<descartes.left_boundary(5));
      assert(Rational(10001,1000)>descartes.right_boundary(5));
      assert(Rational(10200,100)<descartes.left_boundary(6));
      assert(Rational(10400,100)>descartes.right_boundary(6));
      assert(Rational(13000,100)<descartes.left_boundary(7));
      assert(Rational(14100,100)>descartes.right_boundary(7));
  }
  
  */

  { // Test for the Exchange-Descartes method

    // Polnomial with one multiple root:
    
    std::stringstream ss("P[11(0,-51960)(1,158454)(2,3015726)(3,-22833405)(4,64277882)(5,-86502719)(6,58622397)(7,-260172)(8,-77833332)(9,85923423)(10,-37885401)(11,1874259)]");
    // This is a polynomial with a 4-fold root at 1/3, and 3 more roots at
    // approximate positions -1.036, -0.104, 17.764

    // In Maple:g := -51960-86502719*x^5+158454*x-22833405*x^3+3015726*x^2+64277882*x^4+1874259*x^11-77833332*x^8+85923423*x^9-37885401*x^10+58622397*x^6-260172*x^7 

    Poly_int1 f;
    ss >> f;
    
    CGAL::CGALi::Exchange_descartes_tag t;
    
    Poly_int1 sq_free_f = CGAL::CGALi::div_utcf(f, CGAL::CGALi::gcd_utcf(f, CGAL::diff(f)));

    // We expect 4 real roots (m), and the gcd of g and g' is 3 (k)
    Bitstream_descartes descartes(t, f, sq_free_f);

    assert(descartes.degree_of_gcd() == 3);
    assert(descartes.square_free_part() == sq_free_f);
    
    int n = descartes.number_of_real_roots();
    
    //std::cout << "num: " << n << std::endl;

    assert(n == 4);

    for(int i = 0; i<n;i++) {
      if(i!=2) {
	assert(! descartes.is_certainly_multiple_root(i));
      }
      else{ 
	assert(! descartes.is_certainly_simple_root(i));
	assert(descartes.is_certainly_multiple_root(i));
      }
    }
    for(int i =0; i<n; i++) {
      while(descartes.right_boundary(i)-descartes.left_boundary(i)
	    >Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }
    assert(Rational(-104,100)<descartes.left_boundary(0));
    assert(Rational(-103,100)>descartes.right_boundary(0));
    assert(Rational(-106,1000)<descartes.left_boundary(1));
    assert(Rational(-104,1000)>descartes.right_boundary(1));
    assert(Rational(33,100)<descartes.left_boundary(2));
    assert(Rational(34,100)>descartes.right_boundary(2));
    assert(Rational(1776,100)<descartes.left_boundary(3));
    assert(Rational(1777,100)>descartes.right_boundary(3));
    
  }
  
  { // Another test for the Exchange-method

    // In Maple:f := y^3 - y^2  -2*y 
    Poly_int1 f(0,-2,-1,1);
    
    CGAL::CGALi::Exchange_descartes_tag t;
    
    Poly_int1 sq_free_f = CGAL::CGALi::div_utcf(f,CGAL::CGALi::gcd_utcf(f, CGAL::diff(f)));

    // We expect 3 real roots (m), and the degree of gcd of g and g' is 0 (k)
    Bitstream_descartes descartes(t, f, sq_free_f);

    assert(descartes.degree_of_gcd() == 0);
    assert(descartes.square_free_part() == sq_free_f);

    int n = descartes.number_of_real_roots();
    
    //std::cout << "num: " << n << std::endl;

    assert(n==3);

    for(int i =0; i<n; i++) {
      while(descartes.right_boundary(i)-descartes.left_boundary(i)
	    >Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }
    assert(Rational(-101,100)<descartes.left_boundary(0));
    assert(Rational(-99,100)>descartes.right_boundary(0));
    assert(Rational(-1,1000)<descartes.left_boundary(1));
    assert(Rational(1,1000)>descartes.right_boundary(1));
    assert(Rational(199,100)<descartes.left_boundary(2));
    assert(Rational(201,100)>descartes.right_boundary(2));
    
  }
  
  //std::cout << std::endl;
  
  { // Test for the Exchange-Descartes method
      
    // Polnomial with three multiple (1(3), 2(2), -4(4)) and two 
    // (10 and -8) single roots:
    
      Poly_int1 fa(Integer(1), Integer(-1));
      Poly_int1 fb(Integer(2), Integer(-1));
      Poly_int1 fc(Integer(-4), Integer(-1));
      Poly_int1 fd(Integer(10), Integer(-1));
      Poly_int1 fe(Integer(-8), Integer(-1));
      
      Poly_int1 f = fa*fa*fa * fb*fb * fc*fc*fc*fc * fd * fe;
    
      Poly_int1 sq_free_f = CGAL::CGALi::div_utcf(f, CGAL::CGALi::gcd_utcf(f, CGAL::diff(f)));

      CGAL::CGALi::Exchange_descartes_tag t;
      
      Bitstream_descartes descartes(t, f, sq_free_f);
      
      assert(descartes.degree_of_gcd() == 6);
      assert(descartes.square_free_part() == sq_free_f);

      int n = descartes.number_of_real_roots();
      
      //std::cout << "num: " << n << std::endl;
      
      assert(n == 5);

      for(int i = 0; i<n;i++) {
          if(i == 0 || i == 4) {
              assert(! descartes.is_certainly_multiple_root(i));
              //std::cout << i << "s: " << descartes.multiplicity_of_root(i) << std::endl;
          }
          else{ 
              assert(! descartes.is_certainly_simple_root(i));
              //std::cout << i << "m: " << descartes.multiplicity_of_root(i) << std::endl;
          }
      }

      for(int i =0; i<n; i++) {
          while(descartes.right_boundary(i)-descartes.left_boundary(i)
                >Rational(1,10000)) 
          {
              descartes.refine_interval(i);
          }
      }

      assert(Rational(-85,10)<descartes.left_boundary(0));
      assert(Rational(-75,10)>descartes.right_boundary(0));
      assert(Rational(-41,10)<descartes.left_boundary(1));
      assert(Rational(-39,10)>descartes.right_boundary(1));
      assert(Rational(1,10)<descartes.left_boundary(2));
      assert(Rational(11,10)>descartes.right_boundary(2));
      assert(Rational(19,100)<descartes.left_boundary(3));
      assert(Rational(21,10)>descartes.right_boundary(3));
      assert(Rational(9990,1000)<descartes.left_boundary(4));
      assert(Rational(10001,1000)>descartes.right_boundary(4));
  }

  //std::cout << std::endl;


  //std::cout << std::endl;

  { // Test for the backshear method
    typedef Event_refinement<Algebraic_real> Ev_refinement;
    Ev_refinement event_refinement;
    // A polynomial with 3 real roots
    Poly_int1 f(19,21,-28,-22);
    event_refinement.add(Algebraic_real(f,Rational(-2),Rational(-1)));
    event_refinement.add(Algebraic_real(f,Rational(-8,10),Rational(-5,10)));
    event_refinement.add(Algebraic_real(f,Rational(0),Rational(1)));
    
    std::stringstream ss("P[14(0,75449)(1,359917)(2,299972)(3,-1003302)(4,-1857360)(5,42392)(6,2091751)(7,1825115)(8,840268)(9,-840578)(10,-2327140)(11,-1130116)(12,705936)(13,746328)(14,170368)]");
    Poly_int1 g;
    ss >> g;
    // g = f^3*h, where h is a polynomial with 3 simple roots 
    // g := 75449+42392*x^5+359917*x-1003302*x^3+299972*x^2-1857360*x^4-1130116*x^11+840268*x^8-840578*x^9-2327140*x^10+705936*x^12+2091751*x^6+1825115*x^7+746328*x^13+170368*x^14
    
    CGAL::CGALi::Backshear_descartes_tag t;

    // We expect 3 event roots, and 3 non-event roots
    Bitstream_descartes descartes(t, g, 6, 3, event_refinement);

    try {
      int k = descartes.degree_of_gcd();
      (void)k;
      assert(false);
    } catch (CGAL::CGALi::Virtual_method_exception ex) {
      // Expected
    }
    try {
      Poly_int1 sq_f = descartes.square_free_part();
      assert(false);
    } catch (CGAL::CGALi::Virtual_method_exception ex) {
      // Expected
    }

    int n = descartes.number_of_real_roots();

    assert(n==6);

    for(int i =0; i<n; i++) {
      while(descartes.right_boundary(i)-descartes.left_boundary(i)
	    >Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }

    // The non-event roots are expected at -0.953, -0.685, 1.009. The events at
    // -1.527, -0.635, 0.890 (approximate)

    assert(! descartes.is_certainly_simple_root(0));
    assert(descartes.is_certainly_multiple_root(0));
    assert(Rational(-1528,1000)<descartes.left_boundary(0));
    assert(Rational(-1527,1000)>descartes.right_boundary(0));

    assert(! descartes.is_certainly_multiple_root(1));
    assert(Rational(-954,1000)<descartes.left_boundary(1));
    assert(Rational(-953,1000)>descartes.right_boundary(1));

    assert(! descartes.is_certainly_multiple_root(2));
    assert(Rational(-686,1000)<descartes.left_boundary(2));
    assert(Rational(-685,1000)>descartes.right_boundary(2));
    
    assert(! descartes.is_certainly_simple_root(3));
    assert(descartes.is_certainly_multiple_root(3));
    assert(Rational(-636,1000)<descartes.left_boundary(3));
    assert(Rational(-635,1000)>descartes.right_boundary(3));

    assert(! descartes.is_certainly_simple_root(4));
    assert(descartes.is_certainly_multiple_root(4));
    assert(Rational(889,1000)<descartes.left_boundary(4));
    assert(Rational(891,1000)>descartes.right_boundary(4));

    assert(! descartes.is_certainly_multiple_root(5));
    assert(Rational(1009,1000)<descartes.left_boundary(5));
    assert(Rational(1010,1000)>descartes.right_boundary(5));

  }
  { // Another test for the square free method with tree from outside

    // More complicated polynomials needed - change the traits class
    typedef 
        CGAL::CGALi::
        Bitstream_descartes_traits_on_vert_line<Poly_int1,
                                                Algebraic_real,
                                                Integer>
      Traits_2;
  
    typedef CGAL::CGALi::Bitstream_descartes_bfs<Traits_2> Bitstream_descartes;

    typedef typename Bitstream_descartes::Bitstream_tree Bitstream_tree;

    std::stringstream ss("P[4(0,P[4(0,-87)(1,47)(2,43)(3,-88)(4,5)])(1,P[3(0,-90)(1,92)(2,-48)(3,13)])(2,P[2(0,-91)(1,53)(2,-10)])(3,P[1(0,-28)(1,-82)])(4,P[0(0,71)])]");

    // In MAPLE: f:=-87+47*x-90*y+43*x^2+92*x*y-91*y^2-88*x^3-48*x^2*y+53*x*y^2-28*y^3+5*x^4+13*x^3*y-10*x^2*y^2-82*x*y^3+71*y^4;

    Poly_int2 f;
    
    ss >> f;

    Poly_int1 r = CGAL::CGALi::resultant(f,CGAL::diff(f));

    Algebraic_real alpha(r,-2,-1);

    Traits_2 traits(alpha);

    CGAL::CGALi::M_k_descartes_tag mk;

    Bitstream_descartes m_k_descartes(mk,f,1,1,traits);

    assert(m_k_descartes.number_of_real_roots()==1);

    // The root is an "opening" x-extreme point

    // Refine the y-coordinate (just to make it harder...)
    while(m_k_descartes.right_boundary(0) - m_k_descartes.left_boundary(0) 
           > Rational(1,100000)) {
      m_k_descartes.refine_interval(0);
    }    

    // Get copy of the tree
    Bitstream_tree tree = m_k_descartes.get_tree().make_unique();


    // Now, find a value inside the isolating interval of alpha, on the right
    // of alpha
    Rational right_bound = alpha.high();
    
    while(alpha.high()==right_bound) {
      alpha.refine();
    }

    Rational on_right = (alpha.high() + right_bound) / 2;

    Algebraic_real beta(on_right);

    // New traits class for this point
    Traits_2 new_traits(beta);

    CGAL::CGALi::Square_free_descartes_tag sq_free;

    Bitstream_descartes sq_free_descartes(sq_free, f, tree, new_traits);

    assert(sq_free_descartes.number_of_real_roots() == 2);

    for( int i = 0 ; i<sq_free_descartes.number_of_real_roots(); i++ ) {
        // Both roots must be inside the isolating interval of the critical point
        assert(sq_free_descartes.left_boundary(i) > m_k_descartes.left_boundary(0));
        assert(sq_free_descartes.right_boundary(i) < m_k_descartes.right_boundary(0));

    }
    
  }

  { // Test for the Inverse transform Descartes
   
    // More complicated polynomials needed - change the traits class
      typedef 
          CGAL::CGALi::
          Bitstream_descartes_traits_on_vert_line<Poly_int1,
                                                  Algebraic_real,
                                                  Integer>
          Traits_2;

      typedef CGAL::CGALi::Bitstream_descartes_bfs<Traits_2> Bitstream_descartes;

      std::stringstream ss("P[5(0,P[6(0,-80)(1,-34)(2,-45)(3,58)(4,72)(5,63)(6,77)])(1,P[5(0,-99)(1,-45)(2,-2)(3,43)(4,-43)(5,-46)])(2,P[4(0,45)(1,-77)(2,-12)(3,8)(4,27)])(3,P[3(0,-38)(1,-62)(2,78)(3,-64)])(4,P[2(0,95)(1,58)(2,-84)])(5,P[1(0,67)(1,-85)])]");
      
      // In MAPLE: f:=-80+8*x^3*y^2-45*x*y-34*x-77*x*y^2+43*x^3*y-62*x*y^3+58*x*y^4-46*x^5*y+27*x^4*y^2-84*x^2*y^4-85*y^5*x+77*x^6-64*y^3*x^3+67*y^5-12*x^2*y^2-99*y-2*x^2*y-43*x^4*y-45*x^2+95*y^4+63*x^5+72*x^4+45*y^2+58*x^3-38*y^3+78*x^2*y^3

      typename CGAL::Fraction_traits<Rational>::Compose compose;
      
      Algebraic_real x(compose(Integer(67),Integer(85)));

      Traits_2 traits(x);

      Poly_int2 f;
      
      ss >> f;
  
      Poly_int2 g(f);
      
      g.reversal();

      CGAL::CGALi::Square_free_descartes_tag sq_free;

      Bitstream_descartes inv_descartes(sq_free, g, traits);

      assert(inv_descartes.number_of_real_roots() == 3);

      
      CGAL::CGALi::Inverse_transform_descartes_tag inv;

      Bitstream_descartes descartes(inv,
                                    f,
                                    inv_descartes);
      
      int n = descartes.number_of_real_roots();
      
      assert(n == 2);

    for(int i = 0; i<n;i++) {
	assert(! descartes.is_certainly_multiple_root(i));
	assert(descartes.is_certainly_simple_root(i));
    }
    for(int i =0; i<n; i++) {
      while(descartes.right_boundary(i)-descartes.left_boundary(i)
	    >Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }
    assert(Rational(-274,1000)<descartes.left_boundary(0));
    assert(Rational(-272,1000)>descartes.right_boundary(0));

    assert(Rational(160,100)<descartes.left_boundary(1));
    assert(Rational(161,100)>descartes.right_boundary(1));
  }

  { // Another test for the Inverse transform Descartes
   
    // More complicated polynomials needed - change the traits class
      typedef 
          CGAL::CGALi::
          Bitstream_descartes_traits_on_vert_line<Poly_int1, 
                                                  Algebraic_real,
                                                  Integer>
          Traits_2;

      typedef CGAL::CGALi::Bitstream_descartes_bfs<Traits_2> Bitstream_descartes;

      std::stringstream ss("P[6(0,P[7(0,-34)(1,73)(2,-51)(3,46)(4,-46)(5,5)(6,-34)(7,82)])(1,P[6(0,-20)(1,-78)(2,-71)(3,65)(4,86)(5,25)(6,96)])(2,P[5(0,73)(1,13)(2,10)(3,-51)(4,-31)(5,-27)])(3,P[4(0,-94)(1,93)(3,60)(4,84)])(4,P[3(0,99)(1,64)(2,-62)(3,-56)])(5,P[2(0,18)(1,-70)(2,-26)])(6,P[1(0,-17)(1,-71)])]");
      
      // In MAPLE: f:=96*x^6*y-34-51*x^3*y^2+84*y^3*x^4-78*x*y+73*x-56*x^3*y^4-27*x^5*y^2+82*x^7+13*x*y^2+65*x^3*y+93*x*y^3+64*x*y^4+25*x^5*y-71*x*y^6-31*x^4*y^2-62*x^2*y^4-70*y^5*x-26*y^5*x^2-34*x^6+60*y^3*x^3+18*y^5+10*x^2*y^2-20*y-71*x^2*y-17*y^6+86*x^4*y-51*x^2+99*y^4+5*x^5-46*x^4+73*y^2+46*x^3-94*y^3

      typename CGAL::Fraction_traits<Rational>::Compose compose;
      
      Algebraic_real x(compose(Integer(-17),Integer(71)));

      Traits_2 traits(x);

      Poly_int2 f;
      
      ss >> f;
  
      Poly_int2 g(f);

      g.translate(Integer(2));

      g.reversal();

      CGAL::CGALi::Square_free_descartes_tag sq_free;

      Bitstream_descartes inv_descartes(sq_free, g, traits);

      assert(inv_descartes.number_of_real_roots() == 4);
      
      for(int i =0; i<inv_descartes.number_of_real_roots(); i++) {
      while(inv_descartes.right_boundary(i)-inv_descartes.left_boundary(i)
	    >Rational(1,10000)) 
	{
	  inv_descartes.refine_interval(i);
	}
    }
 
      CGAL::CGALi::Inverse_transform_descartes_tag inv;

      Bitstream_descartes descartes(inv,
                                    f,
                                    inv_descartes,
                                    2);
      
      int n = descartes.number_of_real_roots();
      
      assert(descartes.inverse_transform_isolator().
                       number_of_real_roots() == 4);

      assert(n == 3);

    for(int i = 0; i<n;i++) {
	assert(! descartes.is_certainly_multiple_root(i));
	assert(descartes.is_certainly_simple_root(i));
    }
    for(int i =0; i<n; i++) {
      while(descartes.right_boundary(i)-descartes.left_boundary(i)
	    >Rational(1,10000)) 
	{
	  descartes.refine_interval(i);
	}
    }

    assert(Rational(-358,100)<descartes.left_boundary(0));
    assert(Rational(-357,100)>descartes.right_boundary(0));

    assert(Rational(-573,1000)<descartes.left_boundary(1));
    assert(Rational(-572,1000)>descartes.right_boundary(1));

    assert(Rational(972,1000)<descartes.left_boundary(2));
    assert(Rational(973,1000)>descartes.right_boundary(2));
  }



}



int main(int argc,char** argv) {
  
#ifndef CGAL_USE_LEDA
#ifndef LiS_HAVE_CORE
  std::cerr << "This tests requires LEDA and/or CORE" << std::endl;
  return 1;
#endif
#endif
#ifdef CGAL_USE_LEDA
  // LEDA TEST
  test_routine<CGAL::LEDA_arithmetic_kernel> ();
#else
  std::cerr << "LEDA tests skipped" << std::endl;
#endif
#ifdef CGAL_USE_CORE
  // CORE TEST
  test_routine<CGAL::CORE_arithmetic_kernel> ();
#else
  std::cerr << "CORE tests skipped" << std::endl;
#endif
  return 0;
  
}
