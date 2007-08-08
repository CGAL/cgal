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
// $URL$
// $Id$
// 
//
// Author(s)     : Sven Schönherr <sven@inf.ethz.ch>

#include <CGAL/Random.h>
#include <ctime>

CGAL_BEGIN_NAMESPACE

// Class implementation (continued)
// ================================

// constructors
Random::
Random( )
    : rand_max_plus_1( RAND_MAX+1.0), val(0)
{
    // get system's time
    std::time_t s;
    std::time( &s);
    unsigned int  seed = s;

    // initialize random numbers generator
    std::srand( seed);
    random_value = get_int(0, 1<<15);
}

Random::
Random( unsigned int  seed)
    : rand_max_plus_1( RAND_MAX+1.0), val(0)
{
    // initialize random numbers generator
    std::srand( seed);
    random_value = get_int(0, 1<<15);
}

// Global variables
// ================
Random  default_random;

CGAL_END_NAMESPACE
