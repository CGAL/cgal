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
// file          : include/CGAL/Cartesian/Ray_d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

#include <CGAL/Cartesian/distance_computations_d.h>

#ifndef CGAL_CARTESIAN_RAY_D_C
#define CGAL_CARTESIAN_RAY_D_C

CGAL_BEGIN_NAMESPACE

template < class R >
_Twotuple< typename RayCd<R CGAL_CTAG>::Point_d > *
RayCd<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple< Point_d >*)PTR;
}

template < class R >
RayCd<R CGAL_CTAG>::RayCd()
{
  PTR = new _Twotuple< Point_d >;
}

template < class R >
RayCd<R CGAL_CTAG>::
RayCd(const RayCd<R CGAL_CTAG>  &r)
  : Handle((Handle&)r)
{}

template < class R >
RayCd<R CGAL_CTAG>::
RayCd(const typename RayCd<R CGAL_CTAG>::Point_d &sp,
      const typename RayCd<R CGAL_CTAG>::Point_d &secondp)
{
  CGAL_kernel_precondition( sp.dimension() == secondp.dimension() );
  PTR = new _Twotuple< Point_d >(sp,secondp);
}

template < class R >
RayCd<R CGAL_CTAG>::
RayCd(const typename RayCd<R CGAL_CTAG>::Point_d &sp,
      const typename RayCd<R CGAL_CTAG>::Direction_d &d)
{
  CGAL_kernel_precondition( sp.dimension() == d.dimension() );
  PTR = new _Twotuple< Point_d >(sp, sp + d.to_vector());
}


template < class R >
inline RayCd<R CGAL_CTAG>::~RayCd()
{}

template < class R >
RayCd<R CGAL_CTAG> &
RayCd<R CGAL_CTAG>::operator=(const RayCd<R CGAL_CTAG> &r)
{
  Handle::operator=(r);
  return *this;
}

template < class R >
inline bool RayCd<R CGAL_CTAG>::operator==(const RayCd<R CGAL_CTAG> &r) const
{
  if (ptr() == r.ptr()) return true; // identical
  return (source() == r.source()) && (direction() == r.direction());
}

template < class R >
inline bool RayCd<R CGAL_CTAG>::operator!=(const RayCd<R CGAL_CTAG> &r) const
{
  return !(*this == r);
}

template < class R >
inline
long RayCd<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
inline
typename RayCd<R CGAL_CTAG>::Point_d
RayCd<R CGAL_CTAG>::start() const
{
  return ptr()->e0;
}

template < class R >
inline
typename RayCd<R CGAL_CTAG>::Point_d
RayCd<R CGAL_CTAG>::source() const
{
  return ptr()->e0;
}

template < class R >
inline
typename RayCd<R CGAL_CTAG>::Point_d
RayCd<R CGAL_CTAG>::second_point() const
{
  return ptr()->e1;
}


template < class R >
CGAL_KERNEL_INLINE
typename RayCd<R CGAL_CTAG>::Point_d
RayCd<R CGAL_CTAG>::point(int i) const
{
  CGAL_kernel_precondition( i >= 0 );
  if (i == 0) return ptr()->e0;
  if (i == 1) return ptr()->e1;
  return source() + FT(i) * (second_point() - source());
}

template < class R >
inline
typename RayCd<R CGAL_CTAG>::Direction_d
RayCd<R CGAL_CTAG>::direction() const
{
  return Direction_d( second_point() - source() );
}

template < class R >
inline
typename RayCd<R CGAL_CTAG>::Line_d
RayCd<R CGAL_CTAG>::supporting_line() const
{
  return Line_d(*this);
}

template < class R >
inline
RayCd<R CGAL_CTAG>
RayCd<R CGAL_CTAG>::opposite() const
{
  return RayCd<R CGAL_CTAG>( source(), - direction() );
}

template < class R >
inline
RayCd<R CGAL_CTAG>
RayCd<R CGAL_CTAG>::
transform(const typename RayCd<R CGAL_CTAG>::Aff_transformation_d &t) const
{
  return RayCd<R CGAL_CTAG>(t.transform(source()),
                            t.transform(second_point()));
}

template < class R >
bool
RayCd<R CGAL_CTAG>::
has_on(const typename RayCd<R CGAL_CTAG>::Point_d &p) const
{
  return (p == source()) ||
         ( collinear(source(), p, second_point())
           && ( Direction_d(p - source()) == direction() ));
}

template < class R >
inline
bool
RayCd<R CGAL_CTAG>::is_degenerate() const
{
  return source() == second_point();
}

template < class R >
inline
bool
RayCd<R CGAL_CTAG>::
collinear_has_on(const typename RayCd<R CGAL_CTAG>::Point_d &p) const
{
  CGAL_kernel_exactness_precondition( collinear(source(), p, second_point()) );

  Comparison_result cx = compare_x(source(), second_point());
  if (cx != EQUAL)
    return cx != compare_x(p, source());

  Comparison_result cy = compare_y(source(), second_point());
  if (cy != EQUAL)
    return cy != compare_y(p, source());

  Comparison_result cz = compare_z(source(), second_point());
  if (cz != EQUAL)
    return cz != compare_z(p, source());

  return true; // p == source()
}


#ifndef CGAL_NO_OSTREAM_INSERT_RAYCd
template < class R >
std::ostream &
operator<<(std::ostream &os, const RayCd<R CGAL_CTAG> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r.start() << ' ' << r.direction();
    case IO::BINARY :
        return os<< r.start() << r.direction();
    default:
        return os << "RayCd(" << r.start() <<  ", " << r.direction() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_RAYCd

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAYCd
template < class R >
std::istream &
operator>>(std::istream &is, RayCd<R CGAL_CTAG> &r)
{
    typename RayCd<R CGAL_CTAG>::Point_d p;
    typename RayCd<R CGAL_CTAG>::Direction_d d;

    is >> p >> d;

    if (is)
        r = RayCd<R CGAL_CTAG>(p, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAYCd

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_RAY_d_C
