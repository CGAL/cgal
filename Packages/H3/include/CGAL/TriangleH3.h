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
// file          : TriangleH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_TRIANGLEH3_H
#define CGAL_TRIANGLEH3_H

#include <CGAL/predicates_on_pointsH3.h>
#include <CGAL/PlaneH3.h>
#include <CGAL/TetrahedronH3.h>

CGAL_BEGIN_NAMESPACE

template < class FT, class RT >
class TriangleH3 : public Handle_for< Threetuple< PointH3<FT,RT> > >
{
public:
  TriangleH3();
  TriangleH3(const PointH3<FT,RT> &p,
             const PointH3<FT,RT> &q,
             const PointH3<FT,RT> &r);

  bool          operator==(const TriangleH3<FT,RT> &t) const;
  bool          operator!=(const TriangleH3<FT,RT> &t) const;

  PlaneH3<FT,RT>
                supporting_plane() const;

  TriangleH3<FT,RT>
                transform(const Aff_transformationH3<FT,RT> &t) const;
  bool          has_on(const PointH3<FT,RT> &p) const;
  bool          nondegenerate_has_on(const PointH3<FT,RT> &p) const;
  bool          is_degenerate() const;


  PointH3<FT,RT> vertex(int i) const;
  PointH3<FT,RT> operator[](int i) const;

  Bbox_3   bbox() const;

};



template < class FT, class RT >
inline
TriangleH3<FT,RT>::TriangleH3()
 : Handle_for< Threetuple< PointH3<FT,RT> > >( Threetuple< PointH3<FT,RT> >())
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
TriangleH3<FT,RT>::TriangleH3(const PointH3<FT,RT> &p,
                              const PointH3<FT,RT> &q,
                              const PointH3<FT,RT> &r)
 : Handle_for<Threetuple<PointH3<FT,RT> > >(Threetuple<PointH3<FT,RT> >(p,q,r))
{}

template < class FT, class RT >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleH3<FT,RT>::operator==(const TriangleH3<FT,RT> &t) const
{
  int i;
  for(i = 0; (i< 3) && (vertex(0) != t.vertex(i) ); i++) {}
  if (i==3)
  {
      return false;
  }
  return ( vertex(1) == t.vertex(i+1) && vertex(2) == t.vertex(i+2) );
}

template < class FT, class RT >
inline
bool
TriangleH3<FT,RT>::operator!=(const TriangleH3<FT,RT> &t) const
{ return !(*this == t); }
template < class FT, class RT >
CGAL_KERNEL_INLINE
PointH3<FT,RT>
TriangleH3<FT,RT>::vertex(int i) const
{
  switch (i)
  {
      case 0:  return ptr->e0;
      case 1:  return ptr->e1;
      case 2:  return ptr->e2;
      default: return vertex(i%3);
  }
  // return PointH3<FT,RT>();
}
template < class FT, class RT >
inline
PointH3<FT,RT>
TriangleH3<FT,RT>::operator[](int i) const
{ return vertex(i); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
PlaneH3<FT,RT>
TriangleH3<FT,RT>::supporting_plane() const
{ return PlaneH3<FT,RT>(vertex(0), vertex(1), vertex(2)); }

template < class FT, class RT >
inline
Bbox_3
TriangleH3<FT,RT>::bbox() const
{ return vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox(); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
TriangleH3<FT,RT>
TriangleH3<FT,RT>::
transform(const Aff_transformationH3<FT,RT> &t) const
{
  return TriangleH3<FT,RT>(t.transform(vertex(0)),
                                t.transform(vertex(1)),
                                t.transform(vertex(2)));
}


#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLEH3
template < class FT, class RT >
std::ostream &operator<<(std::ostream &os, const TriangleH3<FT,RT> &t)
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
template < class FT, class RT >
std::istream &operator>>(std::istream &is, TriangleH3<FT,RT> &t)
{
  PointH3<FT,RT> p, q, r;
  is >> p >> q >> r;
  t = TriangleH3<FT,RT>(p, q, r);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLEH3

template < class FT, class RT >
CGAL_KERNEL_INLINE
bool
TriangleH3<FT,RT>::
nondegenerate_has_on(const PointH3<FT,RT> &p) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  PlaneH3<FT,RT> sup_pl = supporting_plane();
  if ( !sup_pl.has_on(p) )
  {
      return false;
  }
  TetrahedronH3<FT,RT> tetrapak( vertex(0),
                                      vertex(1),
                                      vertex(2),
                                      vertex(0) + sup_pl.orthogonal_vector());
  return tetrapak.has_on_boundary(p);
}

template < class FT, class RT >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleH3<FT,RT>::has_on(const PointH3<FT,RT> &p) const
{
  if (!is_degenerate() )
  {
      return nondegenerate_has_on(p);
  }
  PointH3<FT,RT> minp( vertex(0) );
  PointH3<FT,RT> maxp( vertex(1) );
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
  SegmentH3<FT,RT> s(minp,maxp);
  return s.has_on(p);
}

template < class FT, class RT >
CGAL_KERNEL_INLINE
bool
TriangleH3<FT,RT>::is_degenerate() const
{ return collinear(vertex(0),vertex(1),vertex(2)); }

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGLEH3_H
