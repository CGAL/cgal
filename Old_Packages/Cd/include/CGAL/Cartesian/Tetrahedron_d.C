// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/Tetrahedron_d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

#ifndef CGAL_CARTESIAN_TETRAHEDRON_D_C
#define CGAL_CARTESIAN_TETRAHEDRON_D_C

#include <CGAL/Cartesian/predicates_on_points_d.h>
#include <CGAL/Cartesian/solve_d.h>
#include <vector>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class R >
_Fourtuple< typename TetrahedronCd<R CGAL_CTAG>::Point_d >*  
TetrahedronCd<R CGAL_CTAG>::ptr() const
{
  return (_Fourtuple< Point_d >*)PTR;
}

template < class R >
TetrahedronCd<R CGAL_CTAG>::
TetrahedronCd()
{
  PTR = new _Fourtuple< Point_d >;
}


template < class R >
TetrahedronCd<R CGAL_CTAG>::
TetrahedronCd(const TetrahedronCd<R CGAL_CTAG> &t)
  : Handle(t)
{}


template < class R >
TetrahedronCd<R CGAL_CTAG>::
TetrahedronCd(const typename TetrahedronCd<R CGAL_CTAG>::Point_d &p,
              const typename TetrahedronCd<R CGAL_CTAG>::Point_d &q,
              const typename TetrahedronCd<R CGAL_CTAG>::Point_d &r,
              const typename TetrahedronCd<R CGAL_CTAG>::Point_d &s)
{
  CGAL_kernel_precondition( p.dimension() == q.dimension() );
  CGAL_kernel_precondition( p.dimension() == r.dimension() );
  CGAL_kernel_precondition( p.dimension() == s.dimension() );
  PTR = new _Fourtuple< Point_d >(p, q, r, s);
}

template < class R >
inline
TetrahedronCd<R CGAL_CTAG>::~TetrahedronCd()
{}

template < class R >
TetrahedronCd<R CGAL_CTAG> &
TetrahedronCd<R CGAL_CTAG>::
operator=(const TetrahedronCd<R CGAL_CTAG> &t)
{
  Handle::operator=(t);
  return *this;
}

template < class Point_d >
struct LessCd {
  // cannot reuse it from predicate_classes, because of
  // problems with file inclusions...
  bool operator() (Point_d const &p, Point_d const &q) {
      typename Point_d::const_iterator pi,qi;
      for (pi=p.begin(),qi=q.begin(); pi!=p.end(); ++pi,++qi) {
        if (*pi<*qi) return true;
        if (*qi<*pi) return false;
      }
      return false;
    }
};

template < class R >
bool
TetrahedronCd<R CGAL_CTAG>::
operator==(const TetrahedronCd<R CGAL_CTAG> &t) const
{
  if ( id() == t.id() ) return true;
  if ( orientation() != t.orientation() ) return false;

  std::vector< Point_d > V1;
  std::vector< Point_d > V2;
  std::vector< Point_d >::iterator uniq_end1;
  std::vector< Point_d >::iterator uniq_end2;
  int k;
  for ( k=0; k < 4; k++) V1.push_back( vertex(k));
  for ( k=0; k < 4; k++) V2.push_back( t.vertex(k));
  std::sort(V1.begin(), V1.end(), LessCd<Point_d>());
  std::sort(V2.begin(), V2.end(), LessCd<Point_d>());
  uniq_end1 = std::unique( V1.begin(), V1.end());
  uniq_end2 = std::unique( V2.begin(), V2.end());
  V1.erase( uniq_end1, V1.end());
  V2.erase( uniq_end2, V2.end());
  return V1 == V2;
}

template < class R >
inline
bool
TetrahedronCd<R CGAL_CTAG>::
operator!=(const TetrahedronCd<R CGAL_CTAG> &t) const
{
  return !(*this == t);
}

template < class R >
inline
long TetrahedronCd<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
inline
int
TetrahedronCd<R CGAL_CTAG>::dimension() const
{
  return vertex(0).dimension();
}

template < class R >
typename TetrahedronCd<R CGAL_CTAG>::Point_d
TetrahedronCd<R CGAL_CTAG>::
vertex(int i) const
{
  if (i<0) i=(i%4)+4;
  else if (i>3) i=i%4;
  switch (i)
    {
    case 0: return ptr()->e0;
    case 1: return ptr()->e1;
    case 2: return ptr()->e2;
    default: return ptr()->e3;
    }
}

template < class R >
inline
typename TetrahedronCd<R CGAL_CTAG>::Point_d
TetrahedronCd<R CGAL_CTAG>::
operator[](int i) const
{
  return vertex(i);
}

template < class R >
Orientation
TetrahedronCd<R CGAL_CTAG>::
orientation() const
{
  CGAL_kernel_precondition( dimension()==3 );
  Point_d v[4] = { vertex(0), vertex(1), vertex(2), vertex(3) };
  return CGAL::orientation(v+0, v+4, Cartesian_tag());
}

template < class R >
Oriented_side
TetrahedronCd<R CGAL_CTAG>::
oriented_side(const typename TetrahedronCd<R CGAL_CTAG>::Point_d &p) const
{
  Orientation o = orientation();
  if (o != ZERO)
    return Oriented_side(o * bounded_side(p));

  CGAL_kernel_assertion (!is_degenerate());
  return ON_ORIENTED_BOUNDARY;
}

template < class R >
Bounded_side
TetrahedronCd<R CGAL_CTAG>::
bounded_side(const typename TetrahedronCd<R CGAL_CTAG>::Point_d &p) const
{
  FT alpha, beta, gamma;
  CGAL_kernel_precondition( dimension()==3 );

  solve(vertex(1)-vertex(0), vertex(2)-vertex(0), vertex(3)-vertex(0),
             p - vertex(0), alpha, beta, gamma);
  if (   (alpha < FT(0)) || (beta < FT(0)) || (gamma < FT(0))
      || (alpha + beta + gamma > FT(1)) )
      return ON_UNBOUNDED_SIDE;

  if (   (alpha == FT(0)) || (beta == FT(0)) || (gamma == FT(0))
      || (alpha+beta+gamma == FT(1)) )
    return ON_BOUNDARY;

  return ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
TetrahedronCd<R CGAL_CTAG>::has_on_boundary
  (const typename TetrahedronCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_ORIENTED_BOUNDARY;
}

template < class R >
inline
bool
TetrahedronCd<R CGAL_CTAG>::has_on_positive_side
  (const typename TetrahedronCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
bool
TetrahedronCd<R CGAL_CTAG>::has_on_negative_side
  (const typename TetrahedronCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
TetrahedronCd<R CGAL_CTAG>::has_on_bounded_side
  (const typename TetrahedronCd<R CGAL_CTAG>::Point_d &p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
TetrahedronCd<R CGAL_CTAG>::has_on_unbounded_side
  (const typename TetrahedronCd<R CGAL_CTAG>::Point_d &p) const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
bool
TetrahedronCd<R CGAL_CTAG>::is_degenerate() const
{
  return (orientation() == ZERO);
}

/*
template < class R >
inline
Bbox_d
TetrahedronCd<R CGAL_CTAG>::bbox() const
{
  return vertex(0).bbox() + vertex(1).bbox()
       + vertex(2).bbox() + vertex(3).bbox();
}
*/

template < class R >
inline
TetrahedronCd<R CGAL_CTAG>
TetrahedronCd<R CGAL_CTAG>::transform
  (const typename TetrahedronCd<R CGAL_CTAG>::Aff_transformation_d &t) const
{
  return TetrahedronCd<R CGAL_CTAG>(t.transform(vertex(0)),
                           t.transform(vertex(1)),
                           t.transform(vertex(2)),
                           t.transform(vertex(3)));
}

#ifndef CGAL_NO_OSTREAM_INSERT_TETRAHEDRONCD
template < class R >
std::ostream &
operator<<(std::ostream &os, const TetrahedronCd<R CGAL_CTAG> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3];
    case IO::BINARY :
        return os << t[0]  << t[1]  << t[2] << t[3];
    default:
        os << "TetrahedronCd(" << t[0] <<  ", " << t[1] <<   ", " << t[2];
        os <<  ", " << t[3] << ")";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TETRAHEDRONCD

#ifndef CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRONCD
template < class R >
std::istream &
operator>>(std::istream &is, TetrahedronCd<R CGAL_CTAG> &t)
{
    typename TetrahedronCd<R CGAL_CTAG>::Point_d p, q, r, s;

    is >> p >> q >> r >> s;

    if (is)
        t = TetrahedronCd<R CGAL_CTAG>(p, q, r, s);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRONCD
 
CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_TETRAHEDRON_D_C
