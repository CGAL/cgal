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
// release       : $CGAL_Revision: CGAL-wip $
// release_date  : $CGAL_Date$
//
// file          : src/Random.C
// package       : Random_numbers 2.2.3 (20 Mar 2001)
// chapter       : Random Numbers Generator
//
// source        : web/Random.aw
// revision      : 2.4
// revision_date : 2000/11/01 13:48:56
//
// author(s)     : Sven Schönherr
// maintainer    : Sven Schönherr <sven@inf.ethz.ch>
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
#ifndef CGAL_PROTECT_SYS_TIME_H
#  include <sys/time.h>
#  define CGAL_PROTECT_SYS_TIME_H
#endif

CGAL_BEGIN_NAMESPACE

// Class implementation (continued)
// ================================

// constructors
Random::
Random( )
{
    // get system's microseconds
    timeval tv;
    gettimeofday( &tv, NULL);
    unsigned long  ms = tv.tv_sec*1000000+tv.tv_usec;

    // initialize random numbers generator
    _state[ 0] = _state[ 2] = static_cast<unsigned short>( ms >> 16);
    _state[ 1] =              static_cast<unsigned short>( ms & 65535);
}

Random::
Random( long seed)
{
    // initialize random numbers generator
    _state[ 0] = _state[ 2] = static_cast<unsigned short>( seed >> 16);
    _state[ 1] =              static_cast<unsigned short>( seed & 65535);
}

Random::
Random( State state)
{
    // initialize random numbers generator
    _state[ 0] = state[ 0];
    _state[ 1] = state[ 1];
    _state[ 2] = state[ 2];
}

// state functions
void
Random::
save_state( State& state) const
{
    state[ 0] = _state[ 0];
    state[ 1] = _state[ 1];
    state[ 2] = _state[ 2];
}

void
Random::
restore_state( const State& state)
{
    _state[ 0] = state[ 0];
    _state[ 1] = state[ 1];
    _state[ 2] = state[ 2];
}

// Global variables
// ================
Random  default_random;

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
