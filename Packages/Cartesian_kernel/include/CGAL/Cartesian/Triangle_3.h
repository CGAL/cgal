// Copyright (c) 2000  Utrecht University (The Netherlands),
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_TRIANGLE_3_H
#define CGAL_CARTESIAN_TRIANGLE_3_H

#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class TriangleC3
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Triangle_3           Triangle_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef Threetuple<Point_3>                      Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  TriangleC3() {}

  TriangleC3(const Point_3 &p, const Point_3 &q, const Point_3 &r)
    : base(p, q, r) {}

  bool       operator==(const TriangleC3 &t) const;
  bool       operator!=(const TriangleC3 &t) const;

  Plane_3    supporting_plane() const;

  Triangle_3 transform(const Aff_transformation_3 &t) const
  {
    return TriangleC3<R>(t.transform(vertex(0)),
                t.transform(vertex(1)),
                t.transform(vertex(2)));
  }

  bool       has_on(const Point_3 &p) const;
  bool       is_degenerate() const;

  const Point_3 & vertex(int i) const;
  const Point_3 & operator[](int i) const;

  Bbox_3     bbox() const;
  
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
  return (i==0) ? get(base).e0 :
         (i==1) ? get(base).e1 :
                  get(base).e2;
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
  return CGALi::squared_area(vertex(0), vertex(1), vertex(2), R());
}

template < class R >
inline
typename TriangleC3<R>::Plane_3
TriangleC3<R>::supporting_plane() const
{
  return Plane_3(vertex(0), vertex(1), vertex(2));
}

template < class R >
Bbox_3
TriangleC3<R>::bbox() const
{
  typename R::Construct_bbox_3 construct_bbox_3;
  return construct_bbox_3(vertex(0)) 
    + construct_bbox_3(vertex(1)) 
    + construct_bbox_3(vertex(2));
}

template < class R >
inline
bool
TriangleC3<R>::
has_on(const typename TriangleC3<R>::Point_3 &p) const
{
  return R().has_on_3_object()
               (static_cast<const typename R::Triangle_3>(*this), p);
}

template < class R >
bool
TriangleC3<R>::is_degenerate() const
{
  return collinear(vertex(0),vertex(1),vertex(2));
}

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLEC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const TriangleC3<R> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0]  << t[1]  << t[2];
    default:
        os << "TriangleC3(" << t[0] <<  ", " << t[1] <<   ", " << t[2] <<")";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLEC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLEC3
template < class R >
std::istream &
operator>>(std::istream &is, TriangleC3<R> &t)
{
    typename R::Point_3 p, q, r;

    is >> p >> q >> r;

    if (is)
	t = TriangleC3<R>(p, q, r);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLEC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TRIANGLE_3_H
