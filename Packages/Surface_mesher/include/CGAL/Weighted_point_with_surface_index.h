// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POINT_WITH_SURFACE_INDEX_H
#define CGAL_POINT_WITH_SURFACE_INDEX_H

#include <CGAL/Point_traits.h>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

namespace CGAL {

template <class Weighted_point>
class Weighted_point_with_surface_index : public Weighted_point
{
  typedef Point_traits<Weighted_point> Point_traits;
  typedef typename Point_traits::Bare_point Bare_point;

  BOOST_STATIC_ASSERT((Is_weighted<Weighted_point>::value));
  BOOST_STATIC_ASSERT((::boost::is_same<typename Point_traits::Is_weighted,
                                        Tag_true>::value));
  BOOST_STATIC_ASSERT((::boost::is_same<Bare_point,
                                    typename Weighted_point::Point>::value));

public:
  Weighted_point_with_surface_index() : Weighted_point(), index(0) {}

  Weighted_point_with_surface_index(const Weighted_point& p)
    : Weighted_point(p), index(0) {}

  Weighted_point_with_surface_index(const Bare_point& bp)
    : Weighted_point(bp), index(0) {}

  Weighted_point_with_surface_index(const Bare_point& bp, 
                                    typename Weighted_point::Weight weight,
                                    int i) 
    : Weighted_point(bp,weight), index(i) {}

  Weighted_point_with_surface_index(const Weighted_point_with_surface_index& 
                                                                            pi)
    : Weighted_point(pi), index(pi.surface_index()) {}

//   template <typename RT>
//   Weighted_point_with_surface_index(const RT& x, const RT& y, const RT& z)
//     : Weighted_point(Weighted_point_traits().point(Bare_point(x, y, z))), index(0) {}

  int surface_index() const
  {
    return index;
  }

  void set_surface_index(const int i)
  {
    index = i;
  }
private:
  int index;
}; // end class Weighted_point_with_surface_index

template <class P>
struct Is_weighted< Weighted_point_with_surface_index<P> > 
  : public Is_weighted<P> {};

} // end namespace CGAL

#endif
