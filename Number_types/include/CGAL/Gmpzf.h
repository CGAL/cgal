// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://gaertner@scm.gforge.inria.fr/svn/cgal/trunk/Number_types/include/CGAL/Gmpzf.h $
// $Id: Gmpzf.h 33782 2006-08-25 14:06:31Z gaertner $
// 
//
// Author(s)     : Bernd Gaertner <gaertner@inf.ethz.ch>

#ifndef CGAL_GMPZF_H
#define CGAL_GMPZF_H

// includes
#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <gmp.h>
#include <mpfr.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpz.h>

#include <boost/operators.hpp>

CGAL_BEGIN_NAMESPACE

struct Gmpzf_rep // as in Gmpz.h
{
// FIXME : bug if ~() is called before an mpz_init*() is called.
// not a problem in practice, but not nice.
// maybe the mpz_init_set* functions should move back to Gmpz_rep.
// But then we should use the Storage_traits::construct/get...

  mpz_t mpZ;

  Gmpzf_rep() {}
  ~Gmpzf_rep() { mpz_clear(mpZ); }

private:
  // Make sure it does not get accidentally copied.
  Gmpzf_rep(const Gmpzf_rep &);
  Gmpzf_rep & operator= (const Gmpzf_rep &);
};

// forward declarations
class Gmpz;

// class declaration
// =================
// This is an exact floating point type; it can represent numbers
// of the form m*2^e, where m is of type mpz_t and e (currently)
// of type long

class Gmpzf :
  Handle_for<Gmpzf_rep>,
    boost::ordered_euclidian_ring_operators1< Gmpzf
  , boost::ordered_euclidian_ring_operators2< Gmpzf, int
    > >
{
  typedef Handle_for<Gmpzf_rep> Base; 

public:
  // tags
  // ----
  typedef Tag_true   Has_gcd;
  typedef Tag_true   Has_division;
  typedef Tag_true   Has_sqrt;
  typedef Tag_true   Has_exact_ring_operations;
  typedef Tag_true   Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;

  // exponent type
  // --------------------
  typedef long Exponent; // may overflow, but if it does, the mantissa is
                         // potentially too large to be useful, anyway;
                         // still, repeated squaring of a power of two
                         // quickly brings this type to its limits...

private:
  // data members (mantissa is from Gmpzf_rep)
  // ----------------------------------------
  // Invariant: 
  // - number is in canonical form, i.e.(m,e) == 0 or m is odd
  Exponent  e; 

  // some static helpers
  static const int double_precision;  // std::numeric_limits<double>::digits
  static Gmpzf s;                     // for intermediate computations 

public:
  // access
  // ------
  const mpz_t& Gmpzf::man() const 
  {
    return Ptr()->mpZ;
  }

  mpz_t& Gmpzf::man() // actually, this shouldn't be public
  {
    return ptr()->mpZ;
  }

  const Gmpzf::Exponent& Gmpzf::exp() const
  {
    return e;
  }
 
  // construction
  // ------------

  Gmpzf::Gmpzf( )
    : e(0)
  {
    mpz_init(man());
    CGAL_postcondition(is_canonical());
  }


  Gmpzf::Gmpzf(const mpz_t z)
    : e(0)
  { 
    mpz_init_set(man(), z); 
    canonicalize();
  }  


  Gmpzf::Gmpzf(const Gmpz& n )
    : e(0)
  { 
    mpz_init_set(man(), n.mpz()); 
    canonicalize();
  }


  Gmpzf::Gmpzf( int i)
    : e(0)
  {   
    mpz_init_set_si( man(), i);  
    canonicalize();
  }


  Gmpzf::Gmpzf( long l)
    : e(0)
  {   
    mpz_init_set_si( man(), l);
    canonicalize();
  }


  Gmpzf::Gmpzf( double d)    
  {
    Protect_FPU_rounding<> P(CGAL_FE_TONEAREST);
    if (d == 0) {
      mpz_init (man());
      e = 0;
      return;
    }
    CGAL_assertion(is_finite(d) && is_valid(d));
    int exp;
    double x = std::frexp(d, &exp); // x in [1/2, 1], x*2^exp = d
    mpz_init_set_d (man(), std::ldexp( x, double_precision)); // an integer   
    e = exp - double_precision;
    canonicalize();
  }

  // arithmetics 
  // -----------
  Gmpzf operator-() const; 
  Gmpzf& operator+=( const Gmpzf& b);
  Gmpzf& operator+=( int i);
  Gmpzf& operator-=( const Gmpzf& b);
  Gmpzf& operator-=( int i);
  Gmpzf& operator*=( const Gmpzf& b);
  Gmpzf& operator*=( int i);
  Gmpzf& operator/= (const Gmpzf& b);
  Gmpzf& operator%= (const Gmpzf& b);
  Gmpzf& operator/= (int i);
  Gmpzf& operator%= (int i); 
  bool is_zero() const;
  Sign sign() const;
  Gmpzf exact_division(const Gmpzf& b) const;
  Gmpzf gcd (const Gmpzf& b) const;
  Gmpzf sqrt() const;
  Comparison_result compare (const Gmpzf &b) const;
 
private:
  void canonicalize();
  bool is_canonical() const;
  static void align ( const mpz_t*& a_aligned, const mpz_t*& b_aligned, 
		     Exponent& rexp, const Gmpzf& a, const Gmpzf& b);  
};

double to_double( const Gmpzf& a );
double to_double( const Quotient<Gmpzf > &q);
std::ostream& operator<< (std::ostream& os, const Gmpzf& a);
std::ostream& print (std::ostream& os, const Gmpzf& a);
std::istream&  operator>> ( std::istream& is, Gmpzf& a);
Comparison_result compare (Gmpzf &a, const Gmpzf &b);
bool operator<(const Gmpzf &a, const Gmpzf &b);
bool operator==(const Gmpzf &a, const Gmpzf &b);
bool operator<(const Gmpzf &a, int b);
bool operator==(const Gmpzf &a, int b);
bool operator>(const Gmpzf &a, int b);
Sign sign (const Gmpzf &a);
bool is_finite(const Gmpzf &);
bool is_valid(const Gmpzf &);
io_Operator io_tag(const Gmpzf &);
Gmpzf exact_division( const Gmpzf& a, const Gmpzf& b);
Gmpzf gcd ( const Gmpzf& a, const Gmpzf& b);
Gmpzf gcd ( const Gmpzf& a, int i);
Gmpzf div ( const Gmpzf& a, const Gmpzf& b);
Gmpzf sqrt (  const Gmpzf& b);

#if ! defined( CGAL_DONT_USE_LINK_PRAGMA) && defined( _MSC_VER )
    #pragma comment(lib, "gmp.lib")
    #pragma comment(lib, "mpfr.lib")
#endif 

CGAL_END_NAMESPACE

#endif // CGAL_GMPZF_H

// ===== EOF ==================================================================
