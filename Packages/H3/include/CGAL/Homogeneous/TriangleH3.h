// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// Author(s)     : Stefan Schirra

#ifndef CGAL_TRIANGLEH3_H
#define CGAL_TRIANGLEH3_H

CGAL_BEGIN_NAMESPACE

template < class R_ >
class TriangleH3
  : public R_::template Handle<Threetuple<typename R_::Point_3> >::type
{
  typedef typename R_::RT                   RT;
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Tetrahedron_3        Tetrahedron_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef Threetuple<Point_3>                      rep;
  typedef typename R_::template Handle<rep>::type  base;

  const base& Base() const { return *this; }
  base& Base() { return *this; }

public:
  typedef R_                R;

  TriangleH3() {}

  TriangleH3(const Point_3 &p,
             const Point_3 &q,
             const Point_3 &r)
    : base(p, q, r) {}

  bool          operator==(const TriangleH3<R> &t) const;
  bool          operator!=(const TriangleH3<R> &t) const;

  Plane_3       supporting_plane() const;

  TriangleH3<R> transform(const Aff_transformation_3 &t) const;
  bool          has_on(const Point_3 &p) const;
  bool          nondegenerate_has_on(const Point_3 &p) const;
  bool          is_degenerate() const;

  const Point_3 & vertex(int i) const;
  const Point_3 & operator[](int i) const;

  FT       squared_area() const;

  Bbox_3   bbox() const;
};

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleH3<R>::operator==(const TriangleH3<R> &t) const
{
  int i;
  for(i = 0; (i< 3) && (vertex(0) != t.vertex(i) ); i++) {}
  if (i==3)
  {
      return false;
  }
  return ( vertex(1) == t.vertex(i+1) && vertex(2) == t.vertex(i+2) );
}

template < class R >
inline
bool
TriangleH3<R>::operator!=(const TriangleH3<R> &t) const
{ return !(*this == t); }

template < class R >
CGAL_KERNEL_INLINE
const typename TriangleH3<R>::Point_3 &
TriangleH3<R>::vertex(int i) const
{
  if (i<0) i=(i%3)+3;
  else if (i>2) i=i%3;
  return (i==0) ? get(Base()).e0 :
         (i==1) ? get(Base()).e1 :
                  get(Base()).e2;
}

template < class R >
inline
const typename TriangleH3<R>::Point_3 &
TriangleH3<R>::operator[](int i) const
{ return vertex(i); }

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename TriangleH3<R>::FT
TriangleH3<R>::squared_area() const
{ 
  return CGAL::squared_area(vertex(0), vertex(1), vertex(2), R());
}

template < class R >
CGAL_KERNEL_INLINE
typename TriangleH3<R>::Plane_3
TriangleH3<R>::supporting_plane() const
{ return Plane_3(vertex(0), vertex(1), vertex(2)); }

template < class R >
inline
Bbox_3
TriangleH3<R>::bbox() const
{ return vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox(); }

template < class R >
CGAL_KERNEL_INLINE
TriangleH3<R>
TriangleH3<R>::
transform(const typename TriangleH3<R>::Aff_transformation_3 &t) const
{
  return TriangleH3<R>(t.transform(vertex(0)),
                                t.transform(vertex(1)),
                                t.transform(vertex(2)));
}


#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLEH3
template < class R >
std::ostream &operator<<(std::ostream &os, const TriangleH3<R> &t)
{
  switch(os.iword(IO::mode))
  {
      case IO::ASCII :
          return os << t[0] << ' ' << t[1] << ' ' << t[2];
      case IO::BINARY :
          return os << t[0]  << t[1]  << t[2];
      default:
          os << "TriangleH3(" << t[0] <<  ", " << t[1] <<   ", " << t[2] <<")";
          return os;
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLEH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLEH3
template < class R >
std::istream &operator>>(std::istream &is, TriangleH3<R> &t)
{
  typename R::Point_3 p, q, r;
  is >> p >> q >> r;
  t = TriangleH3<R>(p, q, r);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLEH3

template < class R >
CGAL_KERNEL_INLINE
bool
TriangleH3<R>::
nondegenerate_has_on(const typename TriangleH3<R>::Point_3 &p) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  Plane_3 sup_pl = supporting_plane();
  if ( !sup_pl.has_on(p) )
  {
      return false;
  }
  Tetrahedron_3 tetrapak( vertex(0),
                          vertex(1),
                          vertex(2),
                          vertex(0) + sup_pl.orthogonal_vector());
  return tetrapak.has_on_boundary(p);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleH3<R>::has_on(const typename TriangleH3<R>::Point_3 &p) const
{
  if (!is_degenerate() )
  {
      return nondegenerate_has_on(p);
  }
  Point_3 minp( vertex(0) );
  Point_3 maxp( vertex(1) );
  if (lexicographically_xyz_smaller(vertex(1),vertex(0)) )
  {
      minp = vertex(1);
      maxp = vertex(0);
  }
  if (lexicographically_xyz_smaller(vertex(2),minp ) )
  {
      minp = vertex(2);
  }
  if (lexicographically_xyz_smaller(maxp, vertex(2)) )
  {
      maxp = vertex(2);
  }
  if (minp == maxp)
  {
      return (p == maxp);
  }
  Segment_3 s(minp,maxp);
  return s.has_on(p);
}

template < class R >
CGAL_KERNEL_INLINE
bool
TriangleH3<R>::is_degenerate() const
{ return collinear(vertex(0),vertex(1),vertex(2)); }

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGLEH3_H
