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
// file          : include/CGAL/Cartesian/Segment_d.C
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

#ifndef CGAL_CARTESIAN_SEGMENT_D_C
#define CGAL_CARTESIAN_SEGMENT_D_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Twotuple< typename SegmentCd<R CGAL_CTAG>::Point_d > *
SegmentCd<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple< Point_d >*)PTR;
}

template < class R >
SegmentCd<R CGAL_CTAG>::
SegmentCd()
{
  PTR = new _Twotuple< Point_d >;
}

template < class R >
SegmentCd<R CGAL_CTAG>::
SegmentCd(const SegmentCd<R CGAL_CTAG>  &s) :
  Handle((Handle&)s)
{}

template < class R >
SegmentCd<R CGAL_CTAG>::
SegmentCd(const typename SegmentCd<R CGAL_CTAG>::Point_d &sp,
          const typename SegmentCd<R CGAL_CTAG>::Point_d &ep)
{
  CGAL_kernel_precondition( sp.dimension() == ep.dimension() );
  PTR = new _Twotuple< Point_d >(sp, ep);
}

template < class R >
inline SegmentCd<R CGAL_CTAG>::~SegmentCd()
{}

template < class R >
SegmentCd<R CGAL_CTAG> &
SegmentCd<R CGAL_CTAG>::operator=(const SegmentCd<R CGAL_CTAG> &s)
{
  Handle::operator=(s);
  return *this;
}

template < class R >
inline
bool
SegmentCd<R CGAL_CTAG>::operator==(const SegmentCd<R CGAL_CTAG> &s) const
{
  if (ptr() == s.ptr()) return true; // identical
  return (source() == s.source())  && (target() == s.target());
}

template < class R >
inline
bool
SegmentCd<R CGAL_CTAG>::operator!=(const SegmentCd<R CGAL_CTAG> &s) const
{
  return !(*this == s);
}

template < class R >
long  SegmentCd<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
typename SegmentCd<R CGAL_CTAG>::Point_d
SegmentCd<R CGAL_CTAG>::start() const
{
  return ptr()->e0;
}

template < class R >
typename SegmentCd<R CGAL_CTAG>::Point_d
SegmentCd<R CGAL_CTAG>::end() const
{
  return ptr()->e1;
}

template < class R >
typename SegmentCd<R CGAL_CTAG>::Point_d
SegmentCd<R CGAL_CTAG>::source() const
{
  return ptr()->e0;
}

template < class R >
typename SegmentCd<R CGAL_CTAG>::Point_d
SegmentCd<R CGAL_CTAG>::target() const
{
  return ptr()->e1;
}

template < class R >
inline
typename SegmentCd<R CGAL_CTAG>::Point_d
SegmentCd<R CGAL_CTAG>::min() const
{
  return (lexicographically_d_smaller(source(),target())) ? source()
                                                          : target();
}

template < class R >
inline
typename SegmentCd<R CGAL_CTAG>::Point_d
SegmentCd<R CGAL_CTAG>::max() const
{
  return (lexicographically_d_smaller(source(),target())) ? target()
                                                          : source();
}

template < class R >
inline
typename SegmentCd<R CGAL_CTAG>::Point_d
SegmentCd<R CGAL_CTAG>::vertex(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
inline
typename SegmentCd<R CGAL_CTAG>::Point_d
SegmentCd<R CGAL_CTAG>::point(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
inline
typename SegmentCd<R CGAL_CTAG>::Point_d
SegmentCd<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
inline
typename SegmentCd<R CGAL_CTAG>::FT
SegmentCd<R CGAL_CTAG>::squared_length() const
{
  return squared_distance(target(), source());
}

template < class R >
inline
typename SegmentCd<R CGAL_CTAG>::Direction_d
SegmentCd<R CGAL_CTAG>::direction() const
{
  return Direction_d( target() - source() );
}

template < class R >
inline
typename SegmentCd<R CGAL_CTAG>::Line_d
SegmentCd<R CGAL_CTAG>::supporting_line() const
{
  return Line_d(*this);
}

template < class R >
inline
SegmentCd<R CGAL_CTAG>
SegmentCd<R CGAL_CTAG>::opposite() const
{
  return SegmentCd<R CGAL_CTAG>(target(), source());
}

template < class R >
inline
SegmentCd<R CGAL_CTAG>
SegmentCd<R CGAL_CTAG>::
transform(const typename SegmentCd<R CGAL_CTAG>::Aff_transformation_d &t) const
{
  return SegmentCd<R CGAL_CTAG>(t.transform(source()), t.transform(target()));
}

template < class R >
inline
bool
SegmentCd<R CGAL_CTAG>::is_degenerate() const
{
  return source() == target();
}

/*
template < class R >
inline
Bbox_d
SegmentCd<R CGAL_CTAG>::bbox() const
{
  return source().bbox() + target().bbox();
}
*/

template < class R >
bool SegmentCd<R CGAL_CTAG>::
has_on(const typename SegmentCd<R CGAL_CTAG>::Point_d &p) const
{
  return are_ordered_along_line(source(), p, target());
}

template < class R >
inline
bool
SegmentCd<R CGAL_CTAG>::
collinear_has_on(const typename SegmentCd<R CGAL_CTAG>::Point_d &p) const
{
  return collinear_are_ordered_along_line(source(), p, target());
}

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENTCd
template < class R >
std::ostream &
operator<<(std::ostream &os, const SegmentCd<R CGAL_CTAG> &s)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << s.source() << ' ' << s.target();
    case IO::BINARY :
        return os << s.source() << s.target();
    default:
        return os << "SegmentCd(" << s.source() <<  ", " << s.target() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENTCd

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENTCd
template < class R >
std::istream &
operator>>(std::istream &is, SegmentCd<R CGAL_CTAG> &s)
{
    typename SegmentCd<R CGAL_CTAG>::Point_d p, q;

    is >> p >> q;

    if (is)
        s = SegmentCd<R CGAL_CTAG>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENTCd

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_SEGMENT_D_C
