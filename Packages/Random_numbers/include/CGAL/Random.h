// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Random.h
// package       : $CGAL_Package: Random_numbers $
// chapter       : Random Numbers Generator
//
// source        : web/Random.aw
// revision      : 2.5
// revision_date : 2001/03/21
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : INRIA Sophia-Antipolis
//
// implementation: Random Numbers Generator
// ============================================================================

#ifndef CGAL_RANDOM_H
#define CGAL_RANDOM_H

// includes
// --------
#ifndef CGAL_BASIC_H
#  include <CGAL/basic.h>
#endif

CGAL_BEGIN_NAMESPACE

// Class declaration
// =================
class Random;

// Class interface
// ===============
class Random {
  public:
    // types
    typedef  unsigned short  State[3];                  // 48 Bits
    
    // creation
    Random( );
    Random( long seed);
    Random( State state);
    
    // operations
    bool    get_bool  ( );
    int     get_int   ( int lower, int upper);
    double  get_double( double lower = 0.0, double upper = 1.0);
    
    int     operator () ( int upper);
    
    // state functions
    void       save_state(       State& state) const;
    void    restore_state( const State& state);
    
    // equality test
    bool  operator == ( const Random& rnd) const;

  private:
    // data members
    unsigned short  _state[3];                          // 48 Bits
};

// Global variables
// ================
extern  Random  default_random;

CGAL_END_NAMESPACE

// ============================================================================

// Class implementation (inline functions)
// =======================================
// includes
#ifndef CGAL_PROTECT_CSTDLIB
#  include <cstdlib>
#  define CGAL_PROTECT_CSTDLIB
#endif

CGAL_BEGIN_NAMESPACE

// operations
inline
bool
Random::
get_bool( )
{
    return static_cast<bool>( erand48( _state) < 0.5);
}

inline
int
Random::
get_int( int lower, int upper)
{
    return( lower + static_cast<int>(
        static_cast<double>( upper-lower) * erand48( _state)));
}

inline
double
Random::
get_double( double lower, double upper)
{
    return( lower + ( upper-lower) * erand48( _state));
}

inline
int
Random::
operator () ( int upper)
{
    return( get_int( 0, upper));
}

inline
bool
Random::
operator == ( const Random& rnd) const
{
    return( static_cast<bool>(
                ( _state[ 0] == rnd._state[ 0]) &&
                ( _state[ 1] == rnd._state[ 1]) &&
                ( _state[ 2] == rnd._state[ 2]) ) );
}

CGAL_END_NAMESPACE

#endif // CGAL_RANDOM_H

// ===== EOF ==================================================================
