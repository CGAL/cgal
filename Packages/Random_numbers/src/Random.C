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
// file          : src/Random.C
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

#include <CGAL/Random.h>

// additional includes
#ifndef CGAL_PROTECT_CTIME
#  include <ctime>
#  define CGAL_PROTECT_CTIME
#endif

CGAL_BEGIN_NAMESPACE

// Class implementation (continued)
// ================================

// constructors
Random::
Random( )
    : rand_max_plus_1( RAND_MAX+1.0)
{
    // get system's time
    time_t s;
    time( &s);
    unsigned int  seed = s;

    // initialize random numbers generator
    srand( seed);
}

Random::
Random( unsigned int  seed)
    : rand_max_plus_1( RAND_MAX+1.0)
{
    // initialize random numbers generator
    srand( seed);
}

// Global variables
// ================
Random  default_random;

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
