// Copyright (c) 2000  
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_TRIANGLE_3_H
#define CGAL_CARTESIAN_TRIANGLE_3_H

#include <CGAL/Handle_for.h>
#include <CGAL/array.h>

namespace CGAL {

template <class R_>
class TriangleC3
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Triangle_3           Triangle_3;

  typedef cpp11::array<Point_3, 3>          Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  TriangleC3() {}

  TriangleC3(const Point_3 &p, const Point_3 &q, const Point_3 &r)
    : base(CGAL::make_array(p, q, r)) {}

  bool       operator==(const TriangleC3 &t) const;
  bool       operator!=(const TriangleC3 &t) const;

  Plane_3    supporting_plane() const;

  bool       has_on(const Point_3 &p) const;
  bool       is_degenerate() const;

  const Point_3 & vertex(int i) const;
  const Point_3 & operator[](int i) const;

  FT         squared_area() const;
};

template < class R >
bool
TriangleC3<R>::operator==(const TriangleC3<R> &t) const
{
  if (CGAL::identical(base, t.base))
      return true;

  int i;
  for(i=0; i<3; i++)
    if ( vertex(0) == t.vertex(i) )
       break;

  return (i<3) && vertex(1) == t.vertex(i+1) && vertex(2) == t.vertex(i+2);
}

template < class R >
inline
bool
TriangleC3<R>::operator!=(const TriangleC3<R> &t) const
{
  return !(*this == t);
}

template < class R >
const typename TriangleC3<R>::Point_3 &
TriangleC3<R>::vertex(int i) const
{
  if (i<0) i=(i%3)+3;
  else if (i>2) i=i%3;
  return (i==0) ? get(base)[0] :
         (i==1) ? get(base)[1] :
                  get(base)[2];
}

template < class R >
inline
const typename TriangleC3<R>::Point_3 &
TriangleC3<R>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename TriangleC3<R>::FT
TriangleC3<R>::squared_area() const
{
  return internal::squared_area(vertex(0), vertex(1), vertex(2), R());
}

template < class R >
inline
typename TriangleC3<R>::Plane_3
TriangleC3<R>::supporting_plane() const
{
  return Plane_3(vertex(0), vertex(1), vertex(2));
}

template < class R >
inline
bool
TriangleC3<R>::
has_on(const typename TriangleC3<R>::Point_3 &p) const
{
  return R().has_on_3_object()
               (static_cast<const typename R::Triangle_3&>(*this), p);
}

template < class R >
bool
TriangleC3<R>::is_degenerate() const
{
  return collinear(vertex(0),vertex(1),vertex(2));
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_TRIANGLE_3_H
