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
// file          : include/CGAL/Cartesian/Tetrahedron_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_TETRAHEDRON_3_H
#define CGAL_CARTESIAN_TETRAHEDRON_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Fourtuple.h>
#include <CGAL/Cartesian/solve_3.h>
#include <vector>
#include <functional>

CGAL_BEGIN_NAMESPACE

template <class R_>
class TetrahedronC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle_for< Fourtuple< typename R_::Point_3 > >
{
public:
  typedef R_                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef TetrahedronC3<R CGAL_CTAG>       Self;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Plane_3              Plane_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef TetrahedronC3<R>                      Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Plane_3_base              Plane_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  TetrahedronC3();
  TetrahedronC3(const Self &t);
  TetrahedronC3(const Point_3 &p,
                const Point_3 &q,
                const Point_3 &r,
                const Point_3 &s);
  ~TetrahedronC3() {}

  Point_3    vertex(int i) const;
  Point_3    operator[](int i) const;

  bool       operator==(const Self &t) const;
  bool       operator!=(const Self &t) const;

  Bbox_3     bbox() const;

  Self       transform(const Aff_transformation_3 &t) const;

  Orientation    orientation() const;
  Oriented_side  oriented_side(const Point_3 &p) const;
  Bounded_side   bounded_side(const Point_3 &p) const;

  bool       has_on_boundary(const Point_3 &p) const;
  bool       has_on_positive_side(const Point_3 &p) const;
  bool       has_on_negative_side(const Point_3 &p) const;
  bool       has_on_bounded_side(const Point_3 &p) const;
  bool       has_on_unbounded_side(const Point_3 &p) const;

  bool       is_degenerate() const;

};

template < class R >
CGAL_KERNEL_CTOR_INLINE
TetrahedronC3<R CGAL_CTAG>::TetrahedronC3()
  : Handle_for<Fourtuple< typename R::Point_3> >( Fourtuple< typename R::Point_3>()) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
TetrahedronC3<R CGAL_CTAG>::TetrahedronC3(const TetrahedronC3<R CGAL_CTAG> &t)
  : Handle_for<Fourtuple< typename R::Point_3> >(t) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
TetrahedronC3<R CGAL_CTAG>::
TetrahedronC3(const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &p,
              const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &q,
              const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &r,
              const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &s)
  : Handle_for<Fourtuple< typename R::Point_3> >( Fourtuple< typename R::Point_3>(p, q, r, s)) {}

template < class Point_3 >
struct Less_xyzC3 {
  // cannot reuse it from predicate_classes, because of
  // problems with file inclusions...
  bool operator() (Point_3 const &p, Point_3 const &q) {
      return lexicographically_xyz_smaller_or_equal(p,q);
    }
};

template < class R >
bool
TetrahedronC3<R CGAL_CTAG>::
operator==(const TetrahedronC3<R CGAL_CTAG> &t) const
{
  if ( identical(t) ) return true;
  if ( orientation() != t.orientation() ) return false;

  std::vector< Point_3 > V1;
  std::vector< Point_3 > V2;
  typename std::vector< Point_3 >::iterator uniq_end1;
  typename std::vector< Point_3 >::iterator uniq_end2;
  int k;
  for ( k=0; k < 4; k++) V1.push_back( vertex(k));
  for ( k=0; k < 4; k++) V2.push_back( t.vertex(k));
  std::sort(V1.begin(), V1.end(), Less_xyzC3<Point_3>());
  std::sort(V2.begin(), V2.end(), Less_xyzC3<Point_3>());
  uniq_end1 = std::unique( V1.begin(), V1.end());
  uniq_end2 = std::unique( V2.begin(), V2.end());
  V1.erase( uniq_end1, V1.end());
  V2.erase( uniq_end2, V2.end());
  return V1 == V2;
}

template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::
operator!=(const TetrahedronC3<R CGAL_CTAG> &t) const
{
  return !(*this == t);
}

template < class R >
typename TetrahedronC3<R CGAL_CTAG>::Point_3
TetrahedronC3<R CGAL_CTAG>::
vertex(int i) const
{
  if (i<0) i=(i%4)+4;
  else if (i>3) i=i%4;
  switch (i)
    {
    case 0: return ptr->e0;
    case 1: return ptr->e1;
    case 2: return ptr->e2;
    default: return ptr->e3;
    }
}

template < class R >
inline
typename TetrahedronC3<R CGAL_CTAG>::Point_3
TetrahedronC3<R CGAL_CTAG>::
operator[](int i) const
{
  return vertex(i);
}

template < class R >
Orientation
TetrahedronC3<R CGAL_CTAG>::
orientation() const
{
  return CGAL::orientation(vertex(0), vertex(1), vertex(2), vertex(3));
}

template < class R >
Oriented_side
TetrahedronC3<R CGAL_CTAG>::
oriented_side(const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  Orientation o = orientation();
  if (o != ZERO)
    return Oriented_side(o * bounded_side(p));

  CGAL_assertion (!is_degenerate());
  return ON_ORIENTED_BOUNDARY;
}

template < class R >
Bounded_side
TetrahedronC3<R CGAL_CTAG>::
bounded_side(const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  FT alpha, beta, gamma;

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
TetrahedronC3<R CGAL_CTAG>::has_on_boundary
  (const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  return oriented_side(p) == ON_ORIENTED_BOUNDARY;
}

template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::has_on_positive_side
  (const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::has_on_negative_side
  (const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::has_on_bounded_side
  (const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::has_on_unbounded_side
  (const typename TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
bool
TetrahedronC3<R CGAL_CTAG>::is_degenerate() const
{
  Plane_3 plane(vertex(0), vertex(1), vertex(2));
  return (plane.is_degenerate()) ? true
                                 : plane.has_on_boundary(vertex(3));
}

template < class R >
inline
Bbox_3
TetrahedronC3<R CGAL_CTAG>::bbox() const
{
  return vertex(0).bbox() + vertex(1).bbox()
       + vertex(2).bbox() + vertex(3).bbox();
}

template < class R >
inline
TetrahedronC3<R CGAL_CTAG>
TetrahedronC3<R CGAL_CTAG>::transform
  (const typename TetrahedronC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  return TetrahedronC3<R CGAL_CTAG>(t.transform(vertex(0)),
                           t.transform(vertex(1)),
                           t.transform(vertex(2)),
                           t.transform(vertex(3)));
}

#ifndef CGAL_NO_OSTREAM_INSERT_TETRAHEDRONC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const TetrahedronC3<R CGAL_CTAG> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3];
    case IO::BINARY :
        return os << t[0]  << t[1]  << t[2] << t[3];
    default:
        os << "TetrahedronC3(" << t[0] <<  ", " << t[1] <<   ", " << t[2];
        os <<  ", " << t[3] << ")";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TETRAHEDRONC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRONC3
template < class R >
std::istream &
operator>>(std::istream &is, TetrahedronC3<R CGAL_CTAG> &t)
{
    typename TetrahedronC3<R CGAL_CTAG>::Point_3 p, q, r, s;

    is >> p >> q >> r >> s;

    t = TetrahedronC3<R CGAL_CTAG>(p, q, r, s);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRONC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TETRAHEDRON_3_H
