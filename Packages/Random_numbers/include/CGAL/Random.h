// Copyright (c) 1997-2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>, Sylvain Pion

#ifndef CGAL_RANDOM_H
#define CGAL_RANDOM_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

class Random {
  public:
    // creation
    Random( );
    Random( unsigned int  seed);
    
    // operations
    bool    get_bool  ( );
    int     get_int   ( int lower, int upper);
    double  get_double( double lower = 0.0, double upper = 1.0);

    // Computes a random int value smaller than 2^b.
    // It's supposed to be fast, useful for randomized algorithms.
    // The distribution is not perfectly flat, but this is a sacrifice against
    // efficiency.
    template <int b>
    int get_bits()
    {
	CGAL_assertion(0<b && b<16);
        if (val == 0) {
            random_value = (421U * random_value + 2073U) % 32749U;
            val = random_value;
        }
        int ret = val & ((1<<b)-1);
        val >>= 1; // Shifting by b would be slightly better, but is slower.
        return ret;
    }

    int     operator () ( int upper);
  private:
    // data members
    const double  rand_max_plus_1;
    unsigned int random_value; // Current 15 bits random value.
    unsigned int val; // random_value shifted by used bits.
};

// Global variables
// ================
extern  Random  default_random;

CGAL_END_NAMESPACE

// ============================================================================

// Class implementation (inline functions)
// =======================================
// includes
#  include <cstdlib>

CGAL_BEGIN_NAMESPACE

// operations
inline
bool
Random::
get_bool( )
{
    return( static_cast< bool>( CGAL_CLIB_STD::rand() & 1));
}

inline
int
Random::
get_int( int lower, int upper)
{
    return( lower + static_cast< int>(
      ( static_cast< double>( upper) - lower) * 
      CGAL_CLIB_STD::rand() / rand_max_plus_1));
}

inline
double
Random::
get_double( double lower, double upper)
{
    return( lower + ( ( upper-lower) * 
		      CGAL_CLIB_STD::rand() / rand_max_plus_1));
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
