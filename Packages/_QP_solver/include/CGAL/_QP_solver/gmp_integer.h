// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef GMP_INTEGER_H
#define GMP_INTEGER_H

// includes
#include <gmp.h>
#include <iostream>

namespace GMP {

// Class declaration
// =================
class Integer;

// Function declaration
// ====================
std::ostream&  operator << ( std::ostream&, const Integer&);


// Class interface
// ===============
class Integer {
  public:
    // construction
    Integer( );
    Integer( const Integer&);
    Integer( int);
    Integer( signed long);
    Integer( unsigned long);
    Integer( double);

    // destruction
    ~Integer( );

    // assignment
    Integer&  operator = ( const Integer&);

    // comparisons
    bool  operator == ( const Integer&) const;
    bool  operator != ( const Integer&) const;
    bool  operator <  ( const Integer&) const;
    bool  operator >  ( const Integer&) const;
    bool  operator <= ( const Integer&) const;
    bool  operator >= ( const Integer&) const;

    // arithmetic operations
    Integer  operator - ( ) const;

    Integer  operator + ( const Integer&) const;
    Integer  operator - ( const Integer&) const;
    Integer  operator * ( const Integer&) const;
    Integer  operator / ( const Integer&) const;

    // arithmetic assignment operations
    Integer&  operator += ( const Integer&);
    Integer&  operator -= ( const Integer&);
    Integer&  operator *= ( const Integer&);
    Integer&  operator /= ( const Integer&);

    // shift operations
    Integer  operator << ( unsigned long) const;
    Integer  operator >> ( unsigned long) const;
    
    // shift assignment operations
    Integer&  operator <<= ( unsigned long);
    Integer&  operator >>= ( unsigned long);
    
    // sign function
    int  sign( ) const;

    // conversion functions
    long    to_long  ( ) const;
    double  to_double( ) const;

    // access functions to the internal representation
    size_t         length( ) const;
    unsigned long  zeros ( ) const;

  private:
    // data members
    mpz_t  value;

    // friends
    friend  std::ostream&  operator << ( std::ostream&, const Integer&);
};

// ============================================================================

// Class implementation
// ====================

// construction
// ------------

// default constructor
inline
Integer::
Integer( )
{
    mpz_init( value);
}

// copy constructor
inline
Integer::
Integer( const Integer& integer)
{
    mpz_init_set( value, integer.value);
}

// constructor (int)
inline
Integer::
Integer( int i)
{
    mpz_init_set_si( value, i);
}

// constructor (signed long)
inline
Integer::
Integer( signed long sl)
{
    mpz_init_set_si( value, sl);
}

// constructor (unsigned long)
inline
Integer::
Integer( unsigned long ul)
{
    mpz_init_set_ui( value, ul);
}

// constructor (double)
inline
Integer::
Integer( double d)
{
    mpz_init_set_d ( value, d);
}


// destruction
// -----------

// destructor
inline
Integer::
~Integer( )
{
    mpz_clear( value);
}


// assignment
// ----------
inline
Integer&
Integer::
operator = ( const Integer& integer)
{
    mpz_set( value, integer.value);
    return *this;
}

// comparisons
// -----------

// equal
inline
bool
Integer::
operator == ( const Integer& integer) const
{
    return ( mpz_cmp( value, integer.value) == 0);
}

// not equal
inline
bool
Integer::
operator != ( const Integer& integer) const
{
    return ( mpz_cmp( value, integer.value) != 0);
}

// less
inline
bool
Integer::
operator < ( const Integer& integer) const
{
    return ( mpz_cmp( value, integer.value) < 0);
}

// greater
inline
bool
Integer::
operator > ( const Integer& integer) const
{
    return ( mpz_cmp( value, integer.value) > 0);
}

// less equal
inline
bool
Integer::
operator <= ( const Integer& integer) const
{
    return ( mpz_cmp( value, integer.value) <= 0);
}

// greater equal
inline
bool
Integer::
operator >= ( const Integer& integer) const
{
    return ( mpz_cmp( value, integer.value) >= 0);
}


// arithmetic operations
// ---------------------

// unary minus
inline
Integer
Integer::
operator - ( ) const
{
    Integer  result;
    mpz_neg( result.value, value);
    return result;
}

// addition
inline
Integer
Integer::
operator + ( const Integer& integer) const
{
    Integer  result;
    mpz_add( result.value, value, integer.value);
    return result;
}

// subtraction
inline
Integer
Integer::
operator - ( const Integer& integer) const
{
    Integer  result;
    mpz_sub( result.value, value, integer.value);
    return result;
}

// multiplication
inline
Integer
Integer::
operator * ( const Integer& integer) const
{
    Integer  result;
    mpz_mul( result.value, value, integer.value);
    return result;
}

// division (truncated)
inline
Integer
Integer::
operator / ( const Integer& integer) const
{
    Integer  result;
    mpz_tdiv_q( result.value, value, integer.value);
    return result;
}


// arithmetic assignment operations
// --------------------------------

// addition assignment
inline
Integer&
Integer::
operator += ( const Integer& integer)
{
    mpz_add( value, value, integer.value);
    return *this;
}

// subtraction assignment
inline
Integer&
Integer::
operator -= ( const Integer& integer)
{
    mpz_sub( value, value, integer.value);
    return *this;
}

// multiplication assignment
inline
Integer&
Integer::
operator *= ( const Integer& integer)
{
    mpz_mul( value, value, integer.value);
    return *this;
}

// division assignment
inline
Integer&
Integer::
operator /= ( const Integer& integer)
{
    mpz_tdiv_q( value, value, integer.value);
    return *this;
}


// shift operations
// ----------------

// left shift
inline
Integer
Integer::
operator << ( unsigned long i) const
{
    Integer  result;
    mpz_mul_2exp( result.value, value, i);
    return result;
}

// right shift
inline
Integer
Integer::
operator >> ( unsigned long i) const
{
    Integer  result;
    mpz_tdiv_q_2exp( result.value, value, i);
    return result;
}

    
// shift assignment operations
// ---------------------------

// left shift assignment
inline
Integer&
Integer::
operator <<= ( unsigned long i)
{
    mpz_mul_2exp( value, value, i);
    return *this;
}

// right shift assignment
inline
Integer&
Integer::
operator >>= ( unsigned long i)
{
    mpz_tdiv_q_2exp( value, value, i);
    return *this;
}


// sign function
// -------------
inline
int
Integer::
sign( ) const
{
    return mpz_sgn( value);
}


// conversion functions
// --------------------

// conversion to long
inline
long
Integer::
to_long( ) const
{
    return mpz_get_si( value);
}

// conversion to double
inline
double
Integer::
to_double( ) const
{
    return mpz_get_d( value);
}


// access functions to the internal representation
// -----------------------------------------------

// number of bits
inline
size_t
Integer::
length( ) const
{
    return mpz_sizeinbase( value, 2);
}

// number of trailing zeros
inline
unsigned long
Integer::
zeros ( ) const
{
    if ( mpz_sgn( value) == 0)
	return 0;
    else
	return mpz_scan1( value, 0);
}


// I/O
// ---

// output operator
inline
std::ostream&
operator << ( std::ostream& out, const Integer& integer)
{
    out << mpz_get_str( (char*)0, 10, integer.value);
    return out;
}


}; // namespace GMP

#endif // GMP_INTEGER_H

// ===== EOF ==================================================================
