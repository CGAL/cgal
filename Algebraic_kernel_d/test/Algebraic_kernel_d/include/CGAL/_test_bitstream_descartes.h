// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Polynomial/include/CGAL/Polynomial.h $
// $Id: Polynomial.h 47254 2008-12-06 21:18:27Z afabri $
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Kerber  <mkerber@mpi-inf.mpg.de>
//                 Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

// #include <CGAL/Algebraic_curve_kernel_2/flags.h>

#include <sstream>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_d_1.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel_at_alpha.h>

#include <CGAL/Algebraic_kernel_d_1.h>

namespace CGAL {

namespace internal {

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
  Rational lower_bound(int i){return events[i].low();}
  Rational upper_bound(int i){return events[i].high();}
  void refine(int i) {events[i].refine();}
};

template<typename ArithmeticKernel>
void test_bitstream_descartes() {

  typedef ArithmeticKernel Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Integer Integer;
  typedef typename Arithmetic_kernel::Rational Rational;

  typedef typename CGAL::Polynomial_type_generator<Integer,1>::Type Poly_int1;
  typedef typename CGAL::Polynomial_type_generator<Integer,2>::Type Poly_int2;
  typedef CGAL::internal::Algebraic_real_d_1<Integer,Rational>
      Algebraic_real;

  typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits
      <CGAL::internal::Bitstream_coefficient_kernel<Integer> > Traits;

  typedef CGAL::internal::Bitstream_descartes<Traits> Bitstream_descartes;

  Traits traits;



  { // Test for the classical Descartes method
    std::stringstream ss("P[15(0,59738427711)(1,300038251861)(2,-471253844514)(3,538575935875)(4,22286946912)(5,548111721525)(6,-379185895352)(7,296681325489)(8,-464256140044)(9,410194800463)(10,232977578849)(11,-376080509486)(12,521721411895)(13,-100316773723)(14,-171187873598)(15,-189253202432)]");
    Poly_int1 f;
    ss >> f;

    // In Maple: f := 59738427711-189253202432*x^15-171187873598*x^14-100316773723*x^13+521721411895*x^12-376080509486*x^11+232977578849*x^10+410194800463*x^9-464256140044*x^8+296681325489*x^7-379185895352*x^6+548111721525*x^5+22286946912*x^4+538575935875*x^3-471253844514*x^2+300038251861*x

    // We expect 3 roots, at -1.176, -.154 and 1.168

    CGAL::internal::Square_free_descartes_tag t;

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
      while((descartes.right_bound(i)-descartes.left_bound(i))
            > Rational(1,10000))
        {
          descartes.refine_interval(i);
        }
    }
    assert(Rational(-118,100)<descartes.left_bound(0));
    assert(Rational(-117,100)>descartes.right_bound(0));
    assert(Rational(-155,1000)<descartes.left_bound(1));
    assert(Rational(-154,1000)>descartes.right_bound(1));
    assert(Rational(116,100)<descartes.left_bound(2));
    assert(Rational(117,100)>descartes.right_bound(2));
  }

  { // Test for the square free case, created with tree from outside

    typedef typename Bitstream_descartes::Bitstream_tree Bitstream_tree;

    std::stringstream ss("P[13(1,1)(4,1)(6,95)(8,77)(11,-51)(13,31)]");

    Poly_int1 f;

    ss >> f;
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
    Bitstream_tree tree(1,f.begin(),f.end(),
                        typename Bitstream_tree::Monomial_basis_tag());
#else
    Bitstream_tree tree(-2,2,0,f.begin(),f.end(),
                        typename Bitstream_tree::Monomial_basis_tag());
#endif

    // We expect 3 roots, at -1.597, -.388 and 0

    CGAL::internal::Square_free_descartes_tag t;

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
      while((descartes.right_bound(i)-descartes.left_bound(i))
            > Rational(1,10000))
        {
          descartes.refine_interval(i);
        }
    }
    assert(Rational(-160,100)<descartes.left_bound(0));
    assert(Rational(-159,100)>descartes.right_bound(0));
    assert(Rational(-389,1000)<descartes.left_bound(1));
    assert(Rational(-387,1000)>descartes.right_bound(1));
    assert(Rational(-1,100)<descartes.left_bound(2));
    assert(Rational(1,100)>descartes.right_bound(2));
  }



  { // Test for the m-k-Descartes method

    // Polnomial with one multiple root:

    std::stringstream ss("P[11(0,-51960)(1,158454)(2,3015726)(3,-22833405)(4,64277882)(5,-86502719)(6,58622397)(7,-260172)(8,-77833332)(9,85923423)(10,-37885401)(11,1874259)]");
    // This is a polynomial with a 4-fold root at 1/3, and 3 more roots at
    // approximate positions -1.036, -0.104, 17.764

    // In Maple:g := -51960-86502719*x^5+158454*x-22833405*x^3+3015726*x^2+64277882*x^4+1874259*x^11-77833332*x^8+85923423*x^9-37885401*x^10+58622397*x^6-260172*x^7

    Poly_int1 f;
    ss >> f;

    CGAL::internal::M_k_descartes_tag t;

    // We expect 4 real roots (m), and the gcd of g and g' is 3 (k)
    Bitstream_descartes descartes(t, f, 4, 3);

    assert(descartes.degree_of_gcd() == 3);
    try {
      Poly_int1 sq_f = descartes.square_free_part();
      assert(false);
    } catch (const CGAL::internal::Virtual_method_exception&) {
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
      while(descartes.right_bound(i)-descartes.left_bound(i)
            >Rational(1,10000))
        {
          descartes.refine_interval(i);
        }
    }
    assert(Rational(-104,100)<descartes.left_bound(0));
    assert(Rational(-103,100)>descartes.right_bound(0));
    assert(Rational(-106,1000)<descartes.left_bound(1));
    assert(Rational(-104,1000)>descartes.right_bound(1));
    assert(Rational(33,100)<descartes.left_bound(2));
    assert(Rational(34,100)>descartes.right_bound(2));
    assert(Rational(1776,100)<descartes.left_bound(3));
    assert(Rational(1777,100)>descartes.right_bound(3));

  }

  { // Another test for the m-k-method

    // Test for the m-k-Descartes method

    // Polnomial with one multiple root:

    // In Maple:f := y^3 - y^2  -2*y
    Poly_int1 f(0,-2,-1,1);

    CGAL::internal::M_k_descartes_tag t;

    // We expect 3 real roots (m), and the degree of gcd of g and g' is 0 (k)
    Bitstream_descartes descartes(t, f, 3, 0);

    assert(descartes.degree_of_gcd() == 0);
    try {
      Poly_int1 sq_f = descartes.square_free_part();
      assert(false);
    } catch (const CGAL::internal::Virtual_method_exception&) {
      // Expected
    }

    int n = descartes.number_of_real_roots();

    assert(n==3);

    for(int i = 0; i<n;i++) {
      assert(! descartes.is_certainly_multiple_root(i));
    }
    for(int i =0; i<n; i++) {
      while(descartes.right_bound(i)-descartes.left_bound(i)
            >Rational(1,10000))
        {
          descartes.refine_interval(i);
        }
    }
    assert(Rational(-101,100)<descartes.left_bound(0));
    assert(Rational(-99,100)>descartes.right_bound(0));
    assert(Rational(-1,1000)<descartes.left_bound(1));
    assert(Rational(1,1000)>descartes.right_bound(1));
    assert(Rational(199,100)<descartes.left_bound(2));
    assert(Rational(201,100)>descartes.right_bound(2));

  }

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

    CGAL::internal::Backshear_descartes_tag t;

    // We expect 3 event roots, and 3 non-event roots
    Bitstream_descartes descartes(t, g, 6, 3, event_refinement);

    try {
      int k = descartes.degree_of_gcd();
      (void)k;
      assert(false);
    } catch (const CGAL::internal::Virtual_method_exception&) {
      // Expected
    }
    try {
      Poly_int1 sq_f = descartes.square_free_part();
      assert(false);
    } catch (const CGAL::internal::Virtual_method_exception&) {
      // Expected
    }

    int n = descartes.number_of_real_roots();

    assert(n==6);

    for(int i =0; i<n; i++) {
      while(descartes.right_bound(i)-descartes.left_bound(i)
            >Rational(1,10000))
        {
          descartes.refine_interval(i);
        }
    }

    // The non-event roots are expected at -0.953, -0.685, 1.009. The events at
    // -1.527, -0.635, 0.890 (approximate)

    assert(! descartes.is_certainly_simple_root(0));
    assert(descartes.is_certainly_multiple_root(0));
    assert(Rational(-1528,1000)<descartes.left_bound(0));
    assert(Rational(-1527,1000)>descartes.right_bound(0));

    assert(! descartes.is_certainly_multiple_root(1));
    assert(Rational(-954,1000)<descartes.left_bound(1));
    assert(Rational(-953,1000)>descartes.right_bound(1));

    assert(! descartes.is_certainly_multiple_root(2));
    assert(Rational(-686,1000)<descartes.left_bound(2));
    assert(Rational(-685,1000)>descartes.right_bound(2));

    assert(! descartes.is_certainly_simple_root(3));
    assert(descartes.is_certainly_multiple_root(3));
    assert(Rational(-636,1000)<descartes.left_bound(3));
    assert(Rational(-635,1000)>descartes.right_bound(3));

    assert(! descartes.is_certainly_simple_root(4));
    assert(descartes.is_certainly_multiple_root(4));
    assert(Rational(889,1000)<descartes.left_bound(4));
    assert(Rational(891,1000)>descartes.right_bound(4));

    assert(! descartes.is_certainly_multiple_root(5));
    assert(Rational(1009,1000)<descartes.left_bound(5));
    assert(Rational(1010,1000)>descartes.right_bound(5));

  }
  { // Another test for the square free method with tree from outside

    // More complicated polynomials needed - change the traits class

    typedef CGAL::Algebraic_kernel_d_1<Integer> AK_1;

    typedef CGAL::internal::Bitstream_coefficient_kernel_at_alpha
        <AK_1> Bitstream_coefficient_kernel;

    typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits
        <Bitstream_coefficient_kernel > Traits_2;

    typedef CGAL::internal::Bitstream_descartes<Traits_2> Bitstream_descartes;

    typedef typename Bitstream_descartes::Bitstream_tree Bitstream_tree;

    std::stringstream ss("P[4(0,P[4(0,-87)(1,47)(2,43)(3,-88)(4,5)])(1,P[3(0,-90)(1,92)(2,-48)(3,13)])(2,P[2(0,-91)(1,53)(2,-10)])(3,P[1(0,-28)(1,-82)])(4,P[0(0,71)])]");

    // In MAPLE: f:=-87+47*x-90*y+43*x^2+92*x*y-91*y^2-88*x^3-48*x^2*y+53*x*y^2-28*y^3+5*x^4+13*x^3*y-10*x^2*y^2-82*x*y^3+71*y^4;

    AK_1 ak_1;

    Poly_int2 f;

    ss >> f;

    Poly_int1 r = CGAL::internal::resultant(f,CGAL::differentiate(f));

    Algebraic_real alpha(r,-2,-1);

    Bitstream_coefficient_kernel bitstream_coefficient_kernel(&ak_1,alpha);

    Traits_2 traits(bitstream_coefficient_kernel);

    CGAL::internal::M_k_descartes_tag mk;

    Bitstream_descartes m_k_descartes(mk,f,1,1,traits);

    assert(m_k_descartes.number_of_real_roots()==1);

    // The root is an "opening" x-extreme point

    // Refine the y-coordinate (just to make it harder...)
    while(m_k_descartes.right_bound(0) - m_k_descartes.left_bound(0)
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

    Bitstream_coefficient_kernel new_bitstream_coefficient_kernel(&ak_1,beta);

    Traits_2 new_traits(new_bitstream_coefficient_kernel);

    CGAL::internal::Square_free_descartes_tag sq_free;

    Bitstream_descartes sq_free_descartes(sq_free, f, tree, new_traits);

    assert(sq_free_descartes.number_of_real_roots() == 2);

    for( int i = 0 ; i<sq_free_descartes.number_of_real_roots(); i++ ) {
        // Both roots must be inside the isolating interval of the critical point
        assert(sq_free_descartes.left_bound(i) > m_k_descartes.left_bound(0));
        assert(sq_free_descartes.right_bound(i) < m_k_descartes.right_bound(0));

    }

  }
}

} // namespace internal

} //namespace CGAL
