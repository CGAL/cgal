// Copyright (c) 1997-2001
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <ctime>
#include <iostream>
#include <sstream>

namespace CGAL {

// Class implementation (continued)
// ================================

// constructors
CGAL_INLINE_FUNCTION
Random::
Random()
    :  val(0)
{
    // get system's time
    std::time_t s;
    std::time( &s);
    seed = (unsigned int)s;
    // initialize random numbers generator
    rng.seed(static_cast<boost::int32_t>(seed));
    random_value = get_int(0, 1<<15);
}

CGAL_INLINE_FUNCTION
Random::
Random(internal::Random_print_seed)
    :  val(0)
{
    // get system's time
    std::time_t s;
    std::time( &s);
    seed = (unsigned int)s;
    std::cerr << "CGAL::Random()::get_seed() = " << seed << std::endl;
    // initialize random numbers generator
    rng.seed(static_cast<boost::int32_t>(seed));
    random_value = get_int(0, 1<<15);
}

CGAL_INLINE_FUNCTION
Random::
Random( unsigned int  seed)
    : val(0), seed(seed)
{
    // initialize random numbers generator
    rng.seed(static_cast<boost::int32_t>(seed));
    random_value = get_int(0, 1<<15);
}


// seed
CGAL_INLINE_FUNCTION
unsigned int
Random::get_seed () const
{
  return seed;
}

// state
CGAL_INLINE_FUNCTION
void
Random::save_state( Random::State& state) const
{
  std::ostringstream os;
  os << rng;
  state = Random::State(os.str(),random_value, val, seed);
}

CGAL_INLINE_FUNCTION
void
Random::restore_state( const Random::State& state)
{
  std::istringstream is(state.rng);
  is >> rng;
  random_value = state.random_value;
  val = state.val;
  seed = state.seed;
}

} //namespace CGAL
