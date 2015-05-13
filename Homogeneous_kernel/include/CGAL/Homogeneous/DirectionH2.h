// Copyright (c) 1999  
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
// Author(s)     : Stefan Schirra
 
#ifndef CGAL_HOMOGENEOUS_DIRECTION_2_H
#define CGAL_HOMOGENEOUS_DIRECTION_2_H

#include <CGAL/array.h>
#include <CGAL/Handle_for.h>
#include <CGAL/kernel_config.h>
#include <CGAL/number_type_basic.h>

namespace CGAL {

template < class R_ >
class DirectionH2
{
  typedef DirectionH2<R_>                   Self;
  typedef typename R_::FT                   FT;
  typedef typename R_::RT                   RT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Segment_2            Segment_2;

  typedef cpp11::array<RT, 3>               Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:

  typedef R_                                    R;

  typedef const RT& Homogeneous_coordinate_type;

  const Self&
  rep() const
  {
    return static_cast<const Self& >(*this);
  }
  
   DirectionH2() {}

   DirectionH2(const RT& x, const RT& y)
      : base(CGAL::make_array(x, y, RT(1))) {}

   // TODO Not documented : should not exist, not used.
   // we should also change array<RT, 3> -> array<RT, 2>
   DirectionH2(const RT& x, const RT& y, const RT& w )
     : base( w > RT(0) ? CGAL::make_array(x, y, w)
                       : CGAL::make_array<RT>(-x, -y, -w) ) {}

    bool    operator==( const DirectionH2<R>& d) const;
    bool    operator!=( const DirectionH2<R>& d) const;


    Vector_2       to_vector() const;

    const RT & x() const { return CGAL::get_pointee_or_identity(base)[0]; }
    const RT & y() const { return CGAL::get_pointee_or_identity(base)[1]; }

    const RT & delta(int i) const;
    const RT & dx() const { return CGAL::get_pointee_or_identity(base)[0]; }
    const RT & dy() const { return CGAL::get_pointee_or_identity(base)[1]; }

};

template <class R >
CGAL_KERNEL_INLINE
bool
DirectionH2<R>::operator==( const DirectionH2<R>& d) const
{
  return (  ( x() * d.y() == y() * d.x() )
          &&( CGAL_NTS sign( x() ) == CGAL_NTS sign( d.x() ) )
          &&( CGAL_NTS sign( y() ) == CGAL_NTS sign( d.y() ) ) );
}

template <class R >
inline
bool
DirectionH2<R>::operator!=( const DirectionH2<R>& d) const
{ return !(*this == d); }

} //namespace CGAL

#include <CGAL/Homogeneous/predicates_on_directionsH2.h>

namespace CGAL {

template <class R >
CGAL_KERNEL_INLINE
typename DirectionH2<R>::Vector_2
DirectionH2<R>::to_vector() const
{ return Vector_2(dx(), dy()); }

} //namespace CGAL

#endif // CGAL_HOMOGENEOUS_DIRECTION_2_H
