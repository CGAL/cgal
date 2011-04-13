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
// revision      : $Revision$
// revision_date : $Date$
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
    // creation
    Random( );
    Random( unsigned int  seed);
    
    // operations
    bool    get_bool  ( );
    int     get_int   ( int lower, int upper);
    double  get_double( double lower = 0.0, double upper = 1.0);
    
    int     operator () ( int upper);
  private:
    // data members
    const double  rand_max_plus_1;
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
    return( static_cast< bool>( rand() & 1));
}

inline
int
Random::
get_int( int lower, int upper)
{
    return( lower + static_cast< int>(
      ( static_cast< double>( upper) - lower) * rand() / rand_max_plus_1));
}

inline
double
Random::
get_double( double lower, double upper)
{
    return( lower + ( ( upper-lower) * rand() / rand_max_plus_1));
}

inline
int
Random::
operator () ( int upper)
{
    return( get_int( 0, upper));
}

CGAL_END_NAMESPACE

#endif // CGAL_RANDOM_H

// ===== EOF ==================================================================
