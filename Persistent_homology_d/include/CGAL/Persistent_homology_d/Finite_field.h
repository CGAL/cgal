// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Afra Zomorodian

// This package implements finite fields using integers.
// It is to be used for only small fields, where the order
// of the field is small.  For 32 bit machines, this means
// primes up to 2^16.  This is so that the largest number
// in the field is x < (2^16-1) and the largest product
// is x^2 < 2^32-1 and will not cause overflow.
// 
// WARNING:  NO overflow checks, or primality check done.
// only check is on div by 0.

#ifndef CGAL_PERSISTENT_HOMOLOGY_D_FINITE_FIELD_H
#define CGAL_PERSISTENT_HOMOLOGY_D_FINITE_FIELD_H

#include <iostream> // std::ostream

namespace CGAL {
namespace Persistent_homology_d {


template < int prime>
class Finite_field {
private:
  typedef Finite_field< prime> Self;
  unsigned int element;
public:
  Finite_field() {}
  // ANSI guarantees correctness for negative numbers for %, but
  // does not specify the behavior, so result could be negative, hence
  // the implementation here.
  Finite_field( int x) {
    element = ( unsigned int) ((( x % prime) + prime) % prime);    
  }
  // addition
  Finite_field operator+( const Finite_field& rhs) const
  {
    return Self(( element + rhs.element) % prime);
  }
  Finite_field& operator+=(const Finite_field& rhs) const
  {
    element = ( element + rhs.element) % prime;
    return *this;
  }

  // subtraction
  Finite_field operator-( const Finite_field& rhs) const
  {
    return Self(((( element - rhs.element) % prime) + prime) % prime);
  }
  Finite_field& operator-=( const Finite_field& rhs)
  {
    element = ((( element - rhs.element) % prime) + prime) % prime;
    return *this;
  }
  
  // multiplication
  Finite_field operator*( const Finite_field& rhs) const
  {
    return Self(( element * rhs.element) % prime);
  }
  Finite_field& operator*=( const Finite_field& rhs) 
  {
    element = ( element * rhs.element) % prime;
    return *this;
  }
  Finite_field operator/( const Finite_field& rhs) const
  {
    if( rhs.element == 0) {
      fprintf( stderr, "Finite_field:  division by 0.\n");
      exit(1);
    }
    return *this*rhs.inverse();
  }
  Finite_field& operator/=( const Finite_field& rhs) 
  {
    if( rhs.element == 0) {
      fprintf( stderr, "Finite_field:  division by 0.\n");
      exit(1);
    }
    *this = *this*rhs.inverse();
    return *this;
  }

  // integer assignment
  Finite_field& operator=( const unsigned int& rhs) 
  {
    element = (( rhs % prime) + prime) % prime;
    return *this;
  }
  
  // unary negation
  Finite_field operator-() const
  {
    return Self(( -element + prime) % prime);
  }

  // inverse
  Finite_field inverse() const
  {
    // To do this, calculate extended euclid(rhs, prime)
    // 2.142 Algorithm Computing multiplicative
    // inverses in Zn
    // See "Handbook of Applied Cryptography" by
    // Alfred J. Menezes et al page 71.
    // or Algorithm E in Knuth, Vol 1.
    // Follwoing Pate Williams' implementation online
    int a = element;
    int b = prime;
    int x, y, q, r, x1, x2, y1, y2;
    x2 = 1, x1 = 0, y2 = 0, y1 = 1;
    while (b > 0) {
      q = a / b, r = a - q * b;
      x = x2 - q * x1, y = y2 - q * y1;
      a = b, b = r;
      x2 = x1, x1 = x, y2 = y1, y1 = y;
    }
    // d = a, x = x2, y = y2;
    // x2 is the inverse
    return Self(x2);
  }
  
  // casting operator
  operator unsigned int () const { return element; }
};

template < unsigned int prime>
inline std::ostream &operator<< ( std::ostream &out,
				  const Finite_field< prime>& e)
{
  out << static_cast<unsigned int>( e);
}

} // namespace Persistent_homology_d
} // namespace CGAL


#endif // CGAL_PERSISTENT_HOMOLOGY_D_FINITE_FIELD_H
