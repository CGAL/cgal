// Copyright (c) 1998  
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
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_RANDOM_CONVEX_SET_TRAITS_2_H
#define CGAL_RANDOM_CONVEX_SET_TRAITS_2_H 1

#include <CGAL/Point_2.h>

namespace CGAL {

template < class Kernel >
struct Random_convex_set_traits_2 : public Kernel {

  typedef typename Kernel::Point_2      Point_2;
  typedef typename Kernel::Direction_2  Direction_2;
  typedef typename Kernel::FT           FT;

  Random_convex_set_traits_2() : _origin( ORIGIN)
  {}

  const Point_2 &
  origin() const
  { return _origin; }

  struct Max_coordinate
  : public std::unary_function< Point_2, FT >
  {
    FT
    operator()( const Point_2& p) const
    { 
      return (std::max)( CGAL_NTS abs( p.x()), CGAL_NTS abs( p.y())); 
    }
  };

  struct Sum
  : public std::binary_function< Point_2, Point_2, Point_2 >
  {
    Point_2
    operator()( const Point_2& p, const Point_2& q) const
    { return p + (q - ORIGIN); }
  };

  struct Scale
  : public std::binary_function< Point_2, FT, Point_2 >
  {
    Point_2
    operator()( const Point_2& p, const FT& k) const
    { return ORIGIN + (p - ORIGIN) * k; }
  };

  struct Angle_less
  : public std::binary_function< Point_2, Point_2, bool >
  {
    bool
    operator()( const Point_2& p, const Point_2& q) const
    {
      return Direction_2( p - ORIGIN) < Direction_2( q - ORIGIN);
    }
  };

private:
  Point_2 _origin;
};

template <class Kernel>
struct Random_convex_set_traits : public Random_convex_set_traits_2<Kernel>
{};

template < class OutputIterator, class Point_generator >
inline
OutputIterator
random_convex_set_2( std::size_t n,
                     OutputIterator o,
                     const Point_generator& pg)
{
  typedef typename Point_generator::value_type Point_2;
  return CGAL_random_convex_set_2(n, o, pg, reinterpret_cast<Point_2*>(0));
}

template < class OutputIterator, class Point_generator, class R >
inline
OutputIterator
CGAL_random_convex_set_2( std::size_t n,
                          OutputIterator o,
                          const Point_generator& pg,
                          Point_2< R >*)
{
  return random_convex_set_2(
    n, o, pg, Random_convex_set_traits_2< R >());
}


} //namespace CGAL

#endif // ! (CGAL_RANDOM_CONVEX_SET_TRAITS_2_H)
