// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : include/CGAL/Homogeneous/TriangleH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_TRIANGLEH3_H
#define CGAL_TRIANGLEH3_H

CGAL_BEGIN_NAMESPACE

template < class R_ >
class TriangleH3
  : public R_::template Handle<Threetuple<typename R_::Point_3> >::type
{
CGAL_VC7_BUG_PROTECTED
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

public:
  typedef R_                R;

  TriangleH3()
    : base(rep()) {}

  TriangleH3(const Point_3 &p,
             const Point_3 &q,
             const Point_3 &r)
    : base(rep(p,q,r)) {}

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
  return (i==0) ? Ptr()->e0 :
         (i==1) ? Ptr()->e1 :
                  Ptr()->e2;
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
   Vector_3 v1 = vertex(1) - vertex(0);
   Vector_3 v2 = vertex(2) - vertex(0);
   Vector_3 v3 = cross_product(v1, v2);
   return (v3 * v3)/FT(4); 
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
