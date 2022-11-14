// Copyright (c) 2007  INRIA Sophia-Antipolis (France), Max-Planck-Institute
// Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sylvain Pion, Michael Hemmer, Alexander Kobel

#ifndef CGAL_RESIDUE_TYPE_H
#define CGAL_RESIDUE_TYPE_H

#include <CGAL/tss.h>
#include <CGAL/FPU.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/double.h>
#include <iostream>

#include <cfloat>
#include <cmath>

#include <boost/operators.hpp>


namespace CGAL {

class Residue;

Residue operator + (const Residue&);
Residue operator - (const Residue&);

std::ostream& operator << (std::ostream& os, const Residue& p);
std::istream& operator >> (std::istream& is, Residue& p);

/*! \ingroup CGAL_Modular_traits
 * \brief This class represents the Field Z mod p.
 *
 * This class uses the type double for representation.
 * Therefore the value of p is restricted to primes less than 2^26.
 * By default p is set to 67108859.
 *
 * It provides the standard operators +,-,*,/ as well as in&output.
 *
 * \see Modular_traits
 */
class Residue:
    boost::ordered_field_operators1< Residue,
    boost::ordered_field_operators2< Residue, int > >{

public:
  typedef Residue Self;
  typedef Residue NT;

private:

  static const double& get_static_CST_CUT()
  {
    static const double CST_CUT = std::ldexp( 3., 51 );
    return CST_CUT;
  }


  static int& prime_int_internal()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(int, prime_int, 67111067);
    return prime_int;
  }

  static inline int get_prime_int(){
    return prime_int_internal();
  }


  static double& prime_internal()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(double, prime, 67111067.0);
    return prime;
  }

  static inline double get_prime(){
    return prime_internal();
  }

  static double& prime_inv_internal()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(double, prime_inv, 0.000000014900672045640400859667452463541);
    return prime_inv;
  }

  static inline double get_prime_inv(){
    return prime_inv_internal();
  }




    /* Quick integer rounding, valid if a<2^51. for double */
    static inline
    double RES_round (double a){
      // call CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST)
      // before using modular arithmetic
      CGAL_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
      return ( (a + get_static_CST_CUT())  - get_static_CST_CUT());
    }

    /* Big modular reduction (e.g. after multiplication) */
    static inline
    double RES_reduce (double a){
      double result = a - get_prime() * RES_round(a * get_prime_inv());
      CGAL_postcondition(2*result <  get_prime());
      CGAL_postcondition(2*result > -get_prime());
      return result;
    }

    /* Little modular reduction (e.g. after a simple addition). */
    static inline
    double RES_soft_reduce (double a){
      double p = get_prime();
        double b = 2*a;
        return (b>p) ? a-p :
            ((b<-p) ? a+p : a);
    }


    /* -a */
    static inline
    double RES_negate(double a){
        return RES_soft_reduce(-a);
    }


    /* a*b */
    static inline
    double RES_mul (double a, double b){
        double c = a*b;
        return RES_reduce(c);
    }


    /* a+b */
    static inline
    double RES_add (double a, double b){
        double c = a+b;
        return RES_soft_reduce(c);
    }


    /* a^-1, using Bezout (extended Euclidian algorithm). */
    static inline
    double RES_inv (double ri1){
        CGAL_precondition (ri1 != 0.0);

        double bi = 0.0;
        double bi1 = 1.0;
        double ri = get_prime();
        double p, tmp, tmp2;

        Real_embeddable_traits<double>::Abs double_abs;
        while (double_abs(ri1) != 1.0)
        {
            p = RES_round(ri/ri1);
            tmp = bi - p * bi1;
            tmp2 = ri - p * ri1;
            bi = bi1;
            ri = ri1;
            bi1 = tmp;
            ri1 = tmp2;
        };

        return ri1 * RES_soft_reduce(bi1);        /* Quicker !!!! */
    }

    /* a/b */
    static inline
    double RES_div (double a, double b){
        return RES_mul(a, RES_inv(b));
    }

public:
    /*! \brief sets the current prime.
     *
     *  Note that you are going to change a static member!
     *  \pre p is prime, but we abstained from such a test.
     *  \pre 0 < p < 2^26
     *
     */
    static int
    set_current_prime(int p){
      int old_prime = get_prime_int();
      prime_int_internal() = p;
      prime_internal() = double(p);
      prime_inv_internal() =  1.0 / double(p);

      return old_prime;
    }

  /*! \brief returns the current prime.  */
    static int get_current_prime(){
      return get_prime_int();
    }

  int  get_value() const{
    CGAL_precondition(2*x_ <  get_prime());
    CGAL_precondition(2*x_ > -get_prime());
    return int(x_);
  }

private:
    double x_;

public:

    //! constructor of Residue, from int
    Residue(int n = 0){
        x_= RES_reduce(n);
    }

    //! constructor of Residue, from long

    Residue (long n) {
        x_= RES_soft_reduce (static_cast< double > (n % get_prime_int()));
    }

    //! constructor of Residue, from long long
    Residue (long long n) {
        x_= RES_soft_reduce (static_cast< double > (n % get_prime_int()));
    }

    //! Access operator for x, \c const
    const double& x() const { return x_; }
    //! Access operator for x
    double&       x()       { return x_; }

    Self& operator += (const Self& p2) {
        x() = RES_add(x(),p2.x());
        return (*this);
    }
    Self& operator -= (const Self& p2){
        x() = RES_add(x(),RES_negate(p2.x()));
        return (*this);
    }
    Self& operator *= (const Self& p2){
        x() = RES_mul(x(),p2.x());
        return (*this);
    }
    Self& operator /= (const Self& p2) {
        x() = RES_div(x(),p2.x());
        return (*this);
    }
    //
    Self& operator += (int p2) {
        x() = RES_add(x(),Residue(p2).x());
        return (*this);
    }
    Self& operator -= (int p2){
        x() = RES_add(x(),Residue(-p2).x());
        return (*this);
    }

    Self& operator *= (int p2){
        x() = RES_mul(x(),Residue(p2).x());
        return (*this);
    }

    Self& operator /= (int p2) {
        x() = RES_div(x(),Residue(p2).x());
        return (*this);
    }

    friend Self operator + (const Self&);
    friend Self operator - (const Self&);
};

inline Residue operator + (const Residue& p1)
{ return p1; }

inline Residue operator - (const Residue& p1){
    typedef Residue RES;
    Residue r;
    r.x() = RES::RES_negate(p1.x());
    return r;
}

inline bool operator == (const Residue& p1, const Residue& p2)
{ return ( p1.x()==p2.x() ); }
inline bool operator == (const Residue& p1, int p2)
{ return ( p1 == Residue(p2) ); }


inline bool operator < (const Residue& p1, const Residue& p2)
{ return ( p1.x() < p2.x() ); }
inline bool operator < (const Residue& p1, int p2)
{ return ( p1.x() < Residue(p2).x() ); }


// I/O
inline std::ostream& operator << (std::ostream& os, const Residue& p) {
    typedef Residue RES;
    os <<"("<< int(p.x())<<"%"<<RES::get_current_prime()<<")";
    return os;
}

inline std::istream& operator >> (std::istream& is, Residue& p) {
    char ch;
    int prime;

    is >> p.x();
    is >> ch;    // read the %
    is >> prime; // read the prime
    CGAL_precondition(prime==Residue::get_current_prime());
    return is;
}

} //namespace CGAL

#endif // CGAL_RESIDUE_TYPE_H
