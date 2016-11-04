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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>, Sylvain Pion, Andreas Fabri

#ifndef CGAL_RANDOM_H
#define CGAL_RANDOM_H

#include <string>
#include <utility>
#include <CGAL/basic.h>
#include <CGAL/tss.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4244)
#endif
#include <boost/random/uniform_smallint.hpp>
#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

namespace CGAL {

  namespace internal { 
    struct Random_print_seed{};
  }

class Random {
public:
  // types
  
  struct State {
    std::string rng;
    unsigned int random_value, val, seed;
    
    State()
    {}
    
    State(std::string rng, 
          unsigned int random_value, 
          unsigned int val, 
          unsigned int seed)
      : rng(rng), random_value(random_value), val(val), seed(seed)
    {}
  };
  // creation
  CGAL_EXPORT Random( );
  CGAL_EXPORT Random( internal::Random_print_seed );
  CGAL_EXPORT Random( unsigned int  seed );
  
  // seed
  CGAL_EXPORT unsigned int get_seed ( ) const;
    
  // operations
  bool get_bool( )
  {
    return( static_cast< bool>( rng() & 1));
  }


  template <typename IntType>
  IntType
  uniform_smallint(IntType lower, IntType upper)
  {
    // uniform_smallint has a closed interval, CGAL a halfopen
    boost::uniform_smallint<IntType> dist(lower,upper-1);
    boost::variate_generator<boost::rand48&, boost::uniform_smallint<IntType> > generator(rng,dist);
    
    return generator();
  }

  template <typename IntType>
  IntType
  uniform_smallint(IntType lower)
  {
    return uniform_smallint<IntType>(lower,9);
  }

  template <typename IntType>
  IntType
  uniform_smallint()
  {
    return uniform_smallint<IntType>(0,9);
  }

  template <typename IntType>
  IntType
  uniform_int(IntType lower, IntType upper)
  {
    // uniform_int has a closed interval, CGAL a halfopen
    boost::uniform_int<IntType> dist(lower,upper);
    boost::variate_generator<boost::rand48&, boost::uniform_int<IntType> > generator(rng,dist);
    
    return generator();
  }


  template <typename IntType>
  IntType
  uniform_int(IntType lower)
  {
    return uniform_int<IntType>(lower,9);
  }

  template <typename IntType>
  IntType
  uniform_int()
  {
    return uniform_int<IntType>(0,9);
  }
 

  
  template <typename IntType>
  IntType
  operator () (IntType upper)
  {
    return uniform_int<IntType>(0, upper-1);
  }
 
  int
  get_int(int lower, int upper)
  {
    return uniform_int<int>(lower,upper-1);
  }


  template <typename RealType>
  RealType
  uniform_real( RealType lower, RealType upper)
  {
    // uniform_real as well as CGAL have a halfopen interval
    boost::uniform_real<RealType> dist(lower,upper);
    boost::variate_generator<boost::rand48&, boost::uniform_real<RealType> > generator(rng,dist);
    
    return generator();
  }


  template <typename RealType>
  RealType
  uniform_real( RealType lower)
  {
    return uniform_real<RealType>(lower, 1.0);
  }


  template <typename RealType>
  RealType
  uniform_real()
  {
    return uniform_real<RealType>(0.0, 1.0);
  }


  template <typename RealType>
  RealType
  uniform_01()
  {
    // uniform_01 as well as CGAL have a halfopen interval
    boost::uniform_01<RealType> dist;
    boost::variate_generator<boost::rand48&, boost::uniform_01<RealType> > generator(rng,dist);
    
    return generator();
  }


  double
  get_double( double lower = 0.0, double upper = 1.0)
  {
    return uniform_real<double>(lower, upper);
  }

    // state 
    CGAL_EXPORT void save_state( State& state) const;
    CGAL_EXPORT void restore_state( const State& state);

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

    
  bool    operator==(Random rd) const
  {
    return (rng == rd.rng) 
      && (random_value == rd.random_value)
      && (val == rd.val)
      && (seed == rd.seed);
  }

  private:
    // data members
    unsigned int random_value; // Current 15 bits random value.
    unsigned int val; // random_value shifted by used bits.
    unsigned int seed; 
    boost::rand48 rng;
};

inline Random& get_default_random()
{
#if (defined( CGAL_TEST_SUITE ) || defined( CGAL_PRINT_SEED )) && !defined(CGAL_HEADER_ONLY)
  internal::Random_print_seed rps;
  CGAL_STATIC_THREAD_LOCAL_VARIABLE(Random, default_random, rps);
#else
  CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Random, default_random);
#endif
  return default_random;
}

#ifndef CGAL_NO_DEPRECATED_CODE
  namespace { CGAL_DEPRECATED_UNUSED CGAL::Random& default_random = get_default_random(); }
#endif // CGAL_NO_DEPRECATED_CODE



} //namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Random_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_RANDOM_H
