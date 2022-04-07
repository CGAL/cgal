// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
//
// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : test/polynomial_gcd.C
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Dominik Huelse <dominik.huelse@gmx.de>
//                 Tobias Reithmann <treith@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file CGAL/polynomial_gcd.C
  test for computing the gcd for polynomials
*/


#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Random.h>
#include <CGAL/gen_polynomials.h>
#include <cassert>
#include <cstdlib>


// #define WITH_OUTPUT 1
CGAL::Random my_random(4711);

template<class POLY>
void gcd_test(const POLY& f, const POLY& g, const POLY& d) {
    POLY tmp = CGAL::gcd(f, g);
#ifdef WITH_OUTPUT
    ::CGAL::IO::set_pretty_mode(std::cout);
    std::cout << "\ngcd: ";
    std::cout << "\nf(x) = " << f;
    std::cout << "\ng(x) = " << g;
    std::cout << "\ngcd(f,g) = " << tmp << "\n";
    std::cout << "\nd        = " << d << "\n";
    std::cout << "\nd/unit_part(d)= "<<d/d.unit_part()<<"\n";

#endif

    assert( tmp == d/d.unit_part() );
}


template<class POLY>
void gcd_utcf_test(const POLY& f, const POLY& g, const POLY& d) {
    typedef CGAL::Polynomial_traits_d<POLY> PT;
    POLY tmp = typename PT::Gcd_up_to_constant_factor()(f, g);

#ifdef WITH_OUTPUT
    std::cout << "\ngcd_utcf: ";
    std::cout << "\nf(x) = " << f;
    std::cout << "\ng(x) = " << g;
    std::cout << "\ngcd_utcf(f,g) = " << tmp << "\n";
    std::cout << "\nd        = " << typename PT::Canonicalize()(d) << "\n";
#endif
    assert( tmp == typename PT::Canonicalize()(d) );
}

template<class AT>
void univariate_polynomial_test() {
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;

    typedef CGAL::Sqrt_extension<Integer  ,Integer> int_EXT_1;
    typedef CGAL::Sqrt_extension<int_EXT_1,Integer> int_EXT_2;
    typedef CGAL::Sqrt_extension<Rational ,Integer> rat_EXT_1;
    typedef CGAL::Sqrt_extension<rat_EXT_1,Integer> rat_EXT_2;

    {
        // testing univariate polynomials with integer coefficients
        typedef CGAL::Polynomial<Integer> int_POLY;
        typedef Integer NT;
        // special cases:
        gcd_test(int_POLY(0),int_POLY(0),int_POLY(0));
        gcd_test(int_POLY(0),int_POLY(NT(4),NT(2)),int_POLY(NT(4),NT(2)));
        gcd_test(int_POLY(NT(4),NT(2)),int_POLY(0),int_POLY(NT(4),NT(2)));

        gcd_utcf_test(int_POLY(0),int_POLY(0),int_POLY(0));
        gcd_utcf_test(int_POLY(0),int_POLY(NT(4),NT(2)),int_POLY(NT(2),NT(1)));
        gcd_utcf_test(int_POLY(NT(4),NT(2)),int_POLY(0),int_POLY(NT(2),NT(1)));

        int_POLY a, b, c, d, e, ac, bc;
        a = int_POLY(Integer(1), Integer(-5), Integer(1));
        b = int_POLY(Integer(-7), Integer(3));
        c = int_POLY(Integer(5), Integer(4));
        d = int_POLY(Integer(5), Integer(4));
        e = int_POLY(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_test(ac, bc, c);
        // there's no gcd
        gcd_test(ac, b, e);
        // g is already a divisor of f
        gcd_test(ac, c, c);
        // testing gcd_utcf
        gcd_utcf_test(ac, bc, d);

        int_POLY  f, g, h, result;
        int l;

        // random polynomials with integer coefficients
        for(l=0;l<2;l++){
            do{
                f = CGAL::internal::rand_Poly_int<Integer>
                    (my_random.get_int(10,1000));
                g = CGAL::internal::rand_Poly_int<Integer>
                    (my_random.get_int(10,1000));
                result = CGAL::gcd(f, g);
            }
            while(result != int_POLY(1));
            d = CGAL::internal::rand_Poly_int<Integer>(my_random.get_int(10,1000));
            int_POLY p1 = f*d;
            int_POLY p2 = g*d;
            gcd_utcf_test(p1, p2, d);
            gcd_test(p1, p2, d);
        }


    }{
        // testing univariate polynomials with rational coefficients
        typedef CGAL::Polynomial<Rational> rat_POLY;
        rat_POLY a, b, c, d, e, ac, bc;
        a = rat_POLY(Rational(1), Rational(-5), Rational(1));
        b = rat_POLY(Rational(-7), Rational(3));
        c = rat_POLY(Rational(6), Rational(4));
        d = rat_POLY(Rational(3,2), Rational(1));
        e = rat_POLY(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_test(ac, bc, d);
        // there's no gcd
        gcd_test(a, b, e);
        // g is already a divisor of f
        gcd_test(ac, c, d);
        // testing gcd_utcf
        gcd_utcf_test(ac, bc, d);

        rat_POLY  f, g, h, result;
        int l;

        // random polynomials with rational coefficients
        for(l=0;l<2;l++){
            do{
                f = rat_POLY(Rational(my_random.get_int(10,1000),
                                my_random.get_int(10,1000)),
                        Rational(my_random.get_int(10,1000),
                                my_random.get_int(10,1000)),
                        Rational(my_random.get_int(10,1000),
                                my_random.get_int(10,1000)),
                        Rational(my_random.get_int(10,1000),
                                my_random.get_int(10,1000)));
                g = rat_POLY(Rational(my_random.get_int(10,1000),
                                my_random.get_int(10,1000)),
                        Rational(my_random.get_int(10,1000),
                                my_random.get_int(10,1000)),
                        Rational(my_random.get_int(10,1000),
                                my_random.get_int(10,1000)));
                result = CGAL::gcd(f, g);
            }
            while(result != rat_POLY(1));

            d = rat_POLY(Rational(my_random.get_int(10,1000),
                            my_random.get_int(10,1000)),
                    Rational(my_random.get_int(10,1000),
                            my_random.get_int(10,1000)),
                    Rational(my_random.get_int(10,1000),
                            my_random.get_int(10,1000)));
            rat_POLY p1 = f*d;
            rat_POLY p2 = g*d;
            gcd_utcf_test(p1, p2, d);
            gcd_test(p1, p2, d);
        }

    }
    {
        // testing univariate polynomials with sqrt-extensions of
        // integers as coefficients

        typedef CGAL::Polynomial<int_EXT_1> int_EXT_1_POLY;
        int_EXT_1_POLY a, b, c, d, e, ac, bc;
        a = int_EXT_1_POLY(int_EXT_1(1), int_EXT_1(-5,-1,2), int_EXT_1(1));
        b = int_EXT_1_POLY(int_EXT_1(-7,-3,2), int_EXT_1(3,-1,2));
        c = int_EXT_1_POLY(int_EXT_1(3), int_EXT_1(2,1,2));
        d = int_EXT_1_POLY(int_EXT_1(6,-3,2), int_EXT_1(2));
        e = int_EXT_1_POLY(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_utcf_test(ac, bc, d);
        // there's no gcd
        gcd_utcf_test(a, b, e);
        // g is already a divisor of f
        gcd_utcf_test(ac, c, d);

        int_EXT_1_POLY  f, g, h, result;
        int l;
        Integer root;

        // random polynomials with sqrt-extensions of integers as coefficients
        for(l=0;l<2;l++){
            root = CGAL::abs(CGAL::internal::rand_int<Integer>
                    (my_random.get_int(10,1000)));
            do{
                f = CGAL::internal::rand_Poly_sqrt<int_EXT_1, Integer>
                    (my_random.get_int(10,1000), root);
                g = CGAL::internal::rand_Poly_sqrt<int_EXT_1, Integer>
                    (my_random.get_int(10,1000), root);
                result = CGAL::internal::gcd_utcf(f, g);
            }
            while(result != int_EXT_1_POLY(1));

            d = CGAL::internal::rand_Poly_sqrt<int_EXT_1, Integer>
                (my_random.get_int(10,1000), root);
            int_EXT_1_POLY p1 = f*d;
            int_EXT_1_POLY p2 = g*d;
            gcd_utcf_test(p1, p2, d);
        }

    }{
        // testing univariate polynomials with sqrt-extensions of
        // rationals as coefficients
        typedef CGAL::Polynomial<rat_EXT_1> rat_EXT_1_POLY;
        rat_EXT_1_POLY a, b, c, d, e, ac, bc;
        a = rat_EXT_1_POLY(rat_EXT_1(1), rat_EXT_1(-5,-1,2), rat_EXT_1(1));
        b = rat_EXT_1_POLY(rat_EXT_1(-7,-3,2), rat_EXT_1(3,-1,2));
        c = rat_EXT_1_POLY(rat_EXT_1(3), rat_EXT_1(2,1,2));
        d = rat_EXT_1_POLY(rat_EXT_1(Rational(3),Rational(-3,2),
                        Integer(2)), rat_EXT_1(1));
        e = rat_EXT_1_POLY(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_test(ac, bc, d);
        // there's no gcd
        gcd_test(a, b, e);
        // g is already a divisor of f
        gcd_test(ac, c, d);
        // testing gcd_utcf
        gcd_utcf_test(ac, bc, d);
    }
    {
        // testing univariate polynomials with nested sqrt-extensions
        // of integers as coefficients
        typedef CGAL::Polynomial<int_EXT_2> int_EXT_2_POLY;
        int_EXT_2_POLY a, b, c, d, e, ac, bc;
        a = int_EXT_2_POLY(
                int_EXT_2(int_EXT_1(1,1,2), int_EXT_1(3,-1,2), 3),
                int_EXT_2(int_EXT_1(5,1,2)),
                int_EXT_2(int_EXT_1(3,1,2), int_EXT_1(1,1,2), 3));
        b = int_EXT_2_POLY(
                int_EXT_2(int_EXT_1(4), int_EXT_1(2,1,2), 3),
                int_EXT_2(int_EXT_1(7,1,2), int_EXT_1(2,-3,2), 3));
        c = int_EXT_2_POLY(
                int_EXT_2(int_EXT_1(3), int_EXT_1(-1,-3,2), 3),
                int_EXT_2(int_EXT_1(2,1,2)));
        d = int_EXT_2_POLY(
                int_EXT_2(int_EXT_1(6,-3,2), int_EXT_1(4,-5,2), 3),
                int_EXT_2(int_EXT_1(2)));
        e = int_EXT_2_POLY(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_utcf_test(ac, bc, d);
        // there's no gcd
        gcd_utcf_test(a, b, e);
        // g is already a divisor of f
        gcd_utcf_test(ac, c, d);


        int_EXT_2_POLY  f, g, h, result;
        int l;
        Integer root;

// random polynomials with nested sqrt-extensions of integers as coefficients
        for(l=0;l<2;l++){
            root = CGAL::abs(CGAL::internal::rand_int<Integer>
                    (my_random.get_int(10,1000)));
            do{
                f = int_EXT_2_POLY(
                        int_EXT_2(CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                                (my_random.get_int(10,1000),root),
                                CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                                (my_random.get_int(10,1000),root), root),
                        int_EXT_2(CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                                (my_random.get_int(10,1000),root),
                                CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                                (my_random.get_int(10,1000),root), root),
                        int_EXT_2(CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                                (my_random.get_int(10,1000),root),
                                CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                                (my_random.get_int(10,1000),root), root));
                g = int_EXT_2_POLY(
                        int_EXT_2(CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                                (my_random.get_int(10,1000),root),
                                CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                                (my_random.get_int(10,1000),root), root),
                        int_EXT_2(CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                                (my_random.get_int(10,1000),root),
                                CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                                (my_random.get_int(10,1000),root), root));
                result = CGAL::internal::gcd_utcf(f, g);
            }
            while(result.degree() != 0);

            d = int_EXT_2_POLY(
                    int_EXT_2(CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                            (my_random.get_int(1,10),root),
                            CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                            (my_random.get_int(1,10),root), root),
                    int_EXT_2(CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                            (my_random.get_int(1,10),root),
                            CGAL::internal::rand_sqrt<int_EXT_1, Integer>
                            (my_random.get_int(1,10),root), root));
            int_EXT_2_POLY p1 = f*d;
            int_EXT_2_POLY p2 = g*d;
            gcd_utcf_test(p1, p2, d);
        }


    }{
        // testing univariate polynomials with nested sqrt-extensions of
        // rationals as coefficients
        typedef CGAL::Polynomial<rat_EXT_2> rat_EXT_2_POLY;
        rat_EXT_2_POLY a, b, c, d, e, ac, bc;
        a = rat_EXT_2_POLY(
                rat_EXT_2(rat_EXT_1(1,1,2), rat_EXT_1(3,-1,2), 3),
                rat_EXT_2(rat_EXT_1(5,1,2)),
                rat_EXT_2(rat_EXT_1(3,1,2), rat_EXT_1(1,1,2), 3));
        b = rat_EXT_2_POLY(
                rat_EXT_2(rat_EXT_1(4), rat_EXT_1(2,1,2), 3),
                rat_EXT_2(rat_EXT_1(7,1,2), rat_EXT_1(2,-3,2), 3));
        c = rat_EXT_2_POLY(
                rat_EXT_2(rat_EXT_1(3), rat_EXT_1(-1,-3,2), 3),
                rat_EXT_2(rat_EXT_1(2,1,2)));
        d = rat_EXT_2_POLY(
                rat_EXT_2(rat_EXT_1(Rational(3),Rational(-3,2),Integer(2)),
                        rat_EXT_1(Rational(2),Rational(-5,2),Integer(2)), 3),
                rat_EXT_2(rat_EXT_1(1)));
        e = rat_EXT_2_POLY(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_test(ac, bc, d);
        // there's no gcd
        gcd_test(a, b, e);
        // g is already a divisor of f
        gcd_test(ac, c, d);
        // testing gcd_utcf
        gcd_utcf_test(ac, bc, d);
    }
}

template<class AT>
void bivariate_polynomial_test() {
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;

    typedef CGAL::Sqrt_extension<Integer  ,Integer> int_EXT_1;
    typedef CGAL::Sqrt_extension<int_EXT_1,Integer> int_EXT_2;
    typedef CGAL::Sqrt_extension<Rational ,Integer> rat_EXT_1;
    typedef CGAL::Sqrt_extension<rat_EXT_1,Integer> rat_EXT_2;
    {
        // testing bivariate polynomials with integer coefficients
        typedef CGAL::Polynomial<Integer> int_POLY_1;
        typedef CGAL::Polynomial<int_POLY_1> int_POLY_2;

        int_POLY_2 a, b, c, d, e, ac, bc;
        a = int_POLY_2(
                int_POLY_1(Integer(-3), Integer(-5)),
                int_POLY_1(Integer(1), Integer(0), Integer(1)));
        b = int_POLY_2(
                int_POLY_1(Integer(1), Integer(-7)),
                int_POLY_1(Integer(0), Integer(3)));
        c = int_POLY_2(
                int_POLY_1(Integer(0), Integer(4)),
                int_POLY_1(Integer(6)));
        d = int_POLY_2(
                int_POLY_1(Integer(0), Integer(2)),
                int_POLY_1(Integer(3)));
        e = int_POLY_2(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_test(ac, bc, c);
        // there's no gcd
        gcd_test(ac, b, e);
        // g is already a divisor of f
        gcd_test(ac, c, c);
        // testing gcd_utcf
        gcd_utcf_test(ac, bc, d);

        int_POLY_2  f, g, h, result;
        int l;

        // random bivariate polynomials with integers as coefficients
        for(l=0;l<2;l++){
            do{
                f = int_POLY_2(
                        CGAL::internal::rand_Poly_int<Integer>
                        (my_random.get_int(10,100)),
                        CGAL::internal::rand_Poly_int<Integer>
                        (my_random.get_int(10,100)),
                        CGAL::internal::rand_Poly_int<Integer>
                        (my_random.get_int(10,100)));
                g = int_POLY_2(
                        CGAL::internal::rand_Poly_int<Integer>
                        (my_random.get_int(10,100)),
                        CGAL::internal::rand_Poly_int<Integer>
                        (my_random.get_int(10,100)),
                        CGAL::internal::rand_Poly_int<Integer>
                        (my_random.get_int(10,100)),
                        CGAL::internal::rand_Poly_int<Integer>
                        (my_random.get_int(10,100)));
                result = CGAL::gcd(f, g);
            }
            while(result != int_POLY_2(1));

            d = int_POLY_2(
                    CGAL::internal::rand_Poly_int<Integer>
                    (my_random.get_int(10,100)),
                    CGAL::internal::rand_Poly_int<Integer>
                    (my_random.get_int(10,100)),
                    CGAL::internal::rand_Poly_int<Integer>
                    (my_random.get_int(10,100)),
                    CGAL::internal::rand_Poly_int<Integer>
                    (my_random.get_int(10,100)));
            int_POLY_2 p1 = f*d;
            int_POLY_2 p2 = g*d;
            gcd_utcf_test(p1, p2, d);
            gcd_test(p1, p2, d);
        }


    }{
        // testing bivariate polynomials with rationals as coefficients
        typedef CGAL::Polynomial<Rational> rat_POLY_1;
        typedef CGAL::Polynomial<rat_POLY_1> rat_POLY_2;
        rat_POLY_2 a, b, c, d, e, ac, bc;
        a = rat_POLY_2(
                rat_POLY_1(Rational(-3),Rational(-5)),
                rat_POLY_1(Rational(1),Rational(0),Rational(1)));
        b = rat_POLY_2(
                rat_POLY_1(Rational(1),Rational(-7)),
                rat_POLY_1(Rational(0),Rational(3)));
        c = rat_POLY_2(
                rat_POLY_1(Rational(0),Rational(4)),
                rat_POLY_1(Rational(6)));
        d = rat_POLY_2(
                rat_POLY_1(Rational(0),Rational(2,3)),
                rat_POLY_1(Rational(1)));
        e = rat_POLY_2(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_test(ac, bc, d);
        // there's no gcd
        gcd_test(a, b, e);
        // g is already a divisor of f
        gcd_test(ac, c, d);
        // testing gcd_utcf
        gcd_utcf_test(ac, bc, d);
    }{
        // testing bivariate polynomials with sqrt-extensions of
        // integers as coefficients
        typedef CGAL::Polynomial<int_EXT_1> int_EXT_1_POLY_1;
        typedef CGAL::Polynomial<int_EXT_1_POLY_1> int_EXT_1_POLY_2;
        int_EXT_1_POLY_2 a, b, c, d, e, ac, bc;
        a = int_EXT_1_POLY_2(
                int_EXT_1_POLY_1(int_EXT_1(-3), int_EXT_1(-5,-1,2)),
                int_EXT_1_POLY_1(int_EXT_1(1),int_EXT_1(0), int_EXT_1(1)));
        b = int_EXT_1_POLY_2(
                int_EXT_1_POLY_1(int_EXT_1(1), int_EXT_1(-7,-3,2)),
                int_EXT_1_POLY_1(int_EXT_1(0), int_EXT_1(3,-1,2)));
        c = int_EXT_1_POLY_2(
                int_EXT_1_POLY_1(int_EXT_1(0), int_EXT_1(4,1,2)),
                int_EXT_1_POLY_1(int_EXT_1(6,-1,2)));
        d = int_EXT_1_POLY_2(
                int_EXT_1_POLY_1(int_EXT_1(0), int_EXT_1(13,5,2)),
                int_EXT_1_POLY_1(int_EXT_1(17)));
        e = int_EXT_1_POLY_2(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_utcf_test(ac, bc, d);
        // there's no gcd
        gcd_utcf_test(a, b, e);
        // g is already a divisor of f
        gcd_utcf_test(ac, c, d);
    }{
        // testing bivariate polynomials with sqrt-extensions of
        // rationals as coefficients
        typedef CGAL::Polynomial<rat_EXT_1> rat_EXT_1_POLY_1;
        typedef CGAL::Polynomial<rat_EXT_1_POLY_1> rat_EXT_1_POLY_2;
        rat_EXT_1_POLY_2 a, b, c, d, e, ac, bc;
        a = rat_EXT_1_POLY_2(
                rat_EXT_1_POLY_1(rat_EXT_1(-3), rat_EXT_1(-5,-1,2)),
                rat_EXT_1_POLY_1(rat_EXT_1(1),rat_EXT_1(0),rat_EXT_1(1)));
        b = rat_EXT_1_POLY_2(
                rat_EXT_1_POLY_1(rat_EXT_1(1), rat_EXT_1(-7,-3,2)),
                rat_EXT_1_POLY_1(rat_EXT_1(0), rat_EXT_1(3,-1,2)));
        c = rat_EXT_1_POLY_2(
                rat_EXT_1_POLY_1(rat_EXT_1(0), rat_EXT_1(4,1,2)),
                rat_EXT_1_POLY_1(rat_EXT_1(6,-1,2)));
        d = rat_EXT_1_POLY_2(
                rat_EXT_1_POLY_1(rat_EXT_1(0), rat_EXT_1(Rational(13,17),
                                Rational(5,17),2)),
                rat_EXT_1_POLY_1(rat_EXT_1(1)));
        e = rat_EXT_1_POLY_2(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_test(ac, bc, d);
        // there's no gcd
        gcd_test(a, b, e);
        // g is already a divisor of f
        gcd_test(ac, c, d);
        // testing gcd_utcf
        gcd_utcf_test(ac, bc, d);
    }{
        // testing bivariate polynomials with nested sqrt-extensions of
        // integers as coefficients
        typedef CGAL::Polynomial<int_EXT_2> int_EXT_2_POLY_1;
        typedef CGAL::Polynomial<int_EXT_2_POLY_1> int_EXT_2_POLY_2;
        int_EXT_2_POLY_2 a, b, c, d, e, ac, bc;
        a = int_EXT_2_POLY_2(
                int_EXT_2_POLY_1(
                        int_EXT_2(int_EXT_1(2), int_EXT_1(-7), 3),
                        int_EXT_2(int_EXT_1(-5,-1,2))),
                int_EXT_2_POLY_1(
                        int_EXT_2(int_EXT_1(1,1,2), int_EXT_1(3,-1,2), 3),
                        int_EXT_2(int_EXT_1(0)),
                        int_EXT_2(int_EXT_1(3,1,2), int_EXT_1(1,1,2), 3)));
        b = int_EXT_2_POLY_2(
                int_EXT_2_POLY_1(
                        int_EXT_2(int_EXT_1(1,1,2), int_EXT_1(-3), 3),
                        int_EXT_2(int_EXT_1(4), int_EXT_1(2,1,2), 3)),
                int_EXT_2_POLY_1(
                        int_EXT_2(int_EXT_1(0)),
                        int_EXT_2(int_EXT_1(7,1,2), int_EXT_1(2,-3,2), 3)));
        c = int_EXT_2_POLY_2(
                int_EXT_2_POLY_1(
                        int_EXT_2(int_EXT_1(0)),
                        int_EXT_2(int_EXT_1(6), int_EXT_1(-2,-6,2), 3)),
                int_EXT_2_POLY_1(
                        int_EXT_2(int_EXT_1(4,2,2))));
        d = int_EXT_2_POLY_2(
                int_EXT_2_POLY_1(
                        int_EXT_2(int_EXT_1(0)),
                        int_EXT_2(int_EXT_1(6,-3,2), int_EXT_1(4,-5,2), 3)),
                int_EXT_2_POLY_1(
                        int_EXT_2(int_EXT_1(2))));
        e = int_EXT_2_POLY_2(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_utcf_test(ac, bc, d);
        // there's no gcd
        gcd_utcf_test(a, b, e);
        // g is already a divisor of f
        gcd_utcf_test(ac, c, d);
    }{
        // testing bivariate polynomials with nested sqrt-extensions of
        // rationals as coefficients
        typedef CGAL::Polynomial<rat_EXT_2> rat_EXT_2_POLY_1;
        typedef CGAL::Polynomial<rat_EXT_2_POLY_1> rat_EXT_2_POLY_2;
        rat_EXT_2_POLY_2 a, b, c, d, e, ac, bc;
        a = rat_EXT_2_POLY_2(
                rat_EXT_2_POLY_1(
                        rat_EXT_2(rat_EXT_1(2), rat_EXT_1(-7), 3),
                        rat_EXT_2(rat_EXT_1(-5,-1,2))),
                rat_EXT_2_POLY_1(
                        rat_EXT_2(rat_EXT_1(1,1,2), rat_EXT_1(3,-1,2), 3),
                        rat_EXT_2(rat_EXT_1(0)),
                        rat_EXT_2(rat_EXT_1(3,1,2), rat_EXT_1(1,1,2), 3)));
        b = rat_EXT_2_POLY_2(
                rat_EXT_2_POLY_1(
                        rat_EXT_2(rat_EXT_1(1,1,2), rat_EXT_1(-3), 3),
                        rat_EXT_2(rat_EXT_1(4), rat_EXT_1(2,1,2), 3)),
                rat_EXT_2_POLY_1(
                        rat_EXT_2(rat_EXT_1(0)),
                        rat_EXT_2(rat_EXT_1(7,1,2), rat_EXT_1(2,-3,2), 3)));
        c = rat_EXT_2_POLY_2(
                rat_EXT_2_POLY_1(
                        rat_EXT_2(rat_EXT_1(0)),
                        rat_EXT_2(rat_EXT_1(6), rat_EXT_1(-2,-6,2), 3)),
                rat_EXT_2_POLY_1(
                        rat_EXT_2(rat_EXT_1(4,2,2))));
        d = rat_EXT_2_POLY_2(
                rat_EXT_2_POLY_1(
                        rat_EXT_2(rat_EXT_1(0)),
                        rat_EXT_2(rat_EXT_1(Rational(3),Rational(-3,2),
                                        Integer(2)),
                                rat_EXT_1(Rational(2),Rational(-5,2),2),
                                Integer(3))),
                rat_EXT_2_POLY_1(
                        rat_EXT_2(rat_EXT_1(1))));
        e = rat_EXT_2_POLY_2(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_test(ac, bc, d);
        // there's no gcd
        gcd_test(a, b, e);
        // g is already a divisor of f
        gcd_test(ac, c, d);
        // testing gcd_utcf
        gcd_utcf_test(ac, bc, d);
    }
}

template<class AT>
void trivariate_polynomial_test() {
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;

//    typedef CGAL::Sqrt_extension<Integer  ,Integer> int_EXT_1;
//    typedef CGAL::Sqrt_extension<int_EXT_1,Integer> int_EXT_2;
//    typedef CGAL::Sqrt_extension<Rational ,Integer> rat_EXT_1;
//    typedef CGAL::Sqrt_extension<rat_EXT_1,Integer> rat_EXT_2;
    {
        // testing trivariate polynomials with integer coefficients
        typedef CGAL::Polynomial<Integer> int_POLY_1;
        typedef CGAL::Polynomial<int_POLY_1> int_POLY_2;
        typedef CGAL::Polynomial<int_POLY_2> int_POLY_3;
        int_POLY_3 a, b, c, d, e, ac, bc;
        a = int_POLY_3(
                int_POLY_2(
                        int_POLY_1(Integer(4), Integer(2)),
                        int_POLY_1(Integer(0), Integer(0), Integer(5))),
                int_POLY_2(
                        int_POLY_1(Integer(-1), Integer(7)),
                        int_POLY_1(Integer(-4), Integer(3))));
        b = int_POLY_3(
                int_POLY_2(
                        int_POLY_1(Integer(1), Integer(0), Integer(7)),
                        int_POLY_1(Integer(-3), Integer(2))),
                int_POLY_2(
                        int_POLY_1(Integer(0)),
                        int_POLY_1(Integer(-3), Integer(1)),
                        int_POLY_1(Integer(2), Integer(11))));
        c = int_POLY_3(
                int_POLY_2(
                        int_POLY_1(Integer(-15), Integer(3)),
                        int_POLY_1(Integer(0), Integer(6))),
                int_POLY_2(
                        int_POLY_1(Integer(12)),
                        int_POLY_1(Integer(-9), Integer(6))));
        d = int_POLY_3(
                int_POLY_2(
                        int_POLY_1(Integer(-5), Integer(1)),
                        int_POLY_1(Integer(0), Integer(2))),
                int_POLY_2(
                        int_POLY_1(Integer(4)),
                        int_POLY_1(Integer(-3), Integer(2))));
        e = int_POLY_3(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_test(ac, bc, c);
        // there's no gcd
        gcd_test(ac, b, e);
        // g is already a divisor of f
        gcd_test(ac, c, c);
        // testing gcd_utcf
        gcd_utcf_test(ac, bc, c);

    }{
        // testing trivariate polynomials with rational coefficients
        typedef CGAL::Polynomial<Rational> rat_POLY_1;
        typedef CGAL::Polynomial<rat_POLY_1> rat_POLY_2;
        typedef CGAL::Polynomial<rat_POLY_2> rat_POLY_3;
        rat_POLY_3 a, b, c, d, e, ac, bc;
        a = rat_POLY_3(
                rat_POLY_2(
                        rat_POLY_1(Rational(4), Rational(2)),
                        rat_POLY_1(Rational(0), Rational(0), Rational(5))),
                rat_POLY_2(
                        rat_POLY_1(Rational(-1), Rational(7)),
                        rat_POLY_1(Rational(-4), Rational(3))));
        b = rat_POLY_3(
                rat_POLY_2(
                        rat_POLY_1(Rational(1), Rational(0), Rational(7)),
                        rat_POLY_1(Rational(-3), Rational(2))),
                rat_POLY_2(
                        rat_POLY_1(Rational(0)),
                        rat_POLY_1(Rational(-3), Rational(1)),
                        rat_POLY_1(Rational(2), Rational(11))));
        c = rat_POLY_3(
                rat_POLY_2(
                        rat_POLY_1(Rational(-15), Rational(3)),
                        rat_POLY_1(Rational(0), Rational(6))),
                rat_POLY_2(
                        rat_POLY_1(Rational(12)),
                        rat_POLY_1(Rational(-9), Rational(6))));
        d = rat_POLY_3(
                rat_POLY_2(
                        rat_POLY_1(Rational(-5,2), Rational(1,2)),
                        rat_POLY_1(Rational(0), Rational(1))),
                rat_POLY_2(
                        rat_POLY_1(Rational(2)),
                        rat_POLY_1(Rational(-3,2), Rational(1))));
        e = rat_POLY_3(1);
        ac = a*c; bc = b*c;

        // a real gcd exists
        gcd_test(ac, bc, d);
        // there's no gcd
        gcd_test(a, b, e);
        // g is already a divisor of f
        gcd_test(ac, c, c);
        // testing gcd_utcf
        gcd_utcf_test(ac, bc, c);
    }
}



template <class AT>
void polynomial_gcd_test() {
    ::CGAL::IO::set_pretty_mode(std::cout);
    std::cout<<" univariate "<<std::endl;
    univariate_polynomial_test<AT>();
    std::cout<<" bivariate "<<std::endl;
    bivariate_polynomial_test<AT>();
    std::cout<<" trivariate "<<std::endl;
    trivariate_polynomial_test<AT>();
}

int main(){
  // This is the wrong rounding mode for modular arithmetic by intention
  CGAL::Protect_FPU_rounding<> pfr(CGAL_FE_UPWARD);
#ifdef CGAL_USE_LEDA
    polynomial_gcd_test<CGAL::LEDA_arithmetic_kernel>();
#endif // CGAL_USE_LEDA
#ifdef CGAL_USE_CORE
    polynomial_gcd_test<CGAL::CORE_arithmetic_kernel>();
#endif // CGAL_USE_CORE

}


// EOF
