// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_LINEH2_H
#define CGAL_LINEH2_H

#include <CGAL/kernel_config.h>
#include <CGAL/array.h>

namespace CGAL {

template < class R_ >
class LineH2
{
    typedef typename R_::FT                   FT;
    typedef typename R_::RT                   RT;
    typedef typename R_::Point_2              Point_2;
    typedef typename R_::Vector_2             Vector_2;
    typedef typename R_::Direction_2          Direction_2;
    typedef typename R_::Segment_2            Segment_2;
    typedef typename R_::Ray_2                Ray_2;
    typedef typename R_::Line_2               Line_2;

    typedef std::array<RT, 3>               Rep;
    typedef typename R_::template Handle<Rep>::type  Base;

    Base base;

public:

    typedef R_                                    R;

    LineH2() {}
    LineH2(const RT& a, const RT& b, const RT& c)
      : base(CGAL::make_array(a, b, c)) {}

    bool           operator==(const LineH2<R>& l) const;
    bool           operator!=(const LineH2<R>& l) const;

    const RT &     a() const { return get_pointee_or_identity(base)[0]; }
    const RT &     b() const { return get_pointee_or_identity(base)[1]; }
    const RT &     c() const { return get_pointee_or_identity(base)[2]; }

};

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
LineH2<R>::operator==(const LineH2<R>& l) const
{
  if (  (a() * l.c() != l.a() * c() )
      ||(b() * l.c() != l.b() * c() ) )
  {
      return false;
  }
  int sc  = static_cast<int>(CGAL_NTS sign(c()));
  int slc = static_cast<int>(CGAL_NTS sign(l.c()));
  if ( sc == slc )
  {
      if (sc == 0)
          return (  (a()*l.b() == b()*l.a() )
                  &&(CGAL_NTS sign(a() )== CGAL_NTS sign( l.a() ))
                  &&(CGAL_NTS sign(b() )== CGAL_NTS sign( l.b() )) );
      else
          return true;
  }
  else
      return false;
}

template < class R >
inline
bool
LineH2<R>::operator!=(const LineH2<R>& l) const
{ return !(*this == l); }

} //namespace CGAL

#endif // CGAL_LINEH2_H
