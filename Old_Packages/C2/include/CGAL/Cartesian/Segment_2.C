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
// file          : include/CGAL/Cartesian/Segment_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_SEGMENT_2_C
#define CGAL_CARTESIAN_SEGMENT_2_C

#include <CGAL/Cartesian/predicates_on_points_2.h>

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_CTOR_INLINE
SegmentC2<R CGAL_CTAG>::SegmentC2()
{
  new ( static_cast< void*>(ptr)) Twotuple<Point_2>();
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
SegmentC2<R CGAL_CTAG>::SegmentC2(const SegmentC2<R CGAL_CTAG>  &s)
  : Handle_for<Twotuple<typename R::Point_2> >(s)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
SegmentC2<R CGAL_CTAG>::
SegmentC2(const typename SegmentC2<R CGAL_CTAG>::Point_2 &sp,
          const typename SegmentC2<R CGAL_CTAG>::Point_2 &ep)
{
   new ( static_cast< void*>(ptr)) Twotuple<Point_2>(sp, ep);
}

template < class R >
inline
SegmentC2<R CGAL_CTAG>::~SegmentC2()
{}



template < class R >
inline
bool
SegmentC2<R CGAL_CTAG>::operator==(const SegmentC2<R CGAL_CTAG> &s) const
{
  return source() == s.source()  && target() == s.target();
}

template < class R >
inline
bool
SegmentC2<R CGAL_CTAG>::operator!=(const SegmentC2<R CGAL_CTAG> &s) const
{
  return !(*this == s);
}


template < class R >
inline
typename SegmentC2<R CGAL_CTAG>::Point_2
SegmentC2<R CGAL_CTAG>::start() const
{
  return ptr->e0;
}

template < class R >
inline
typename SegmentC2<R CGAL_CTAG>::Point_2
SegmentC2<R CGAL_CTAG>::end() const
{
  return ptr->e1;
}

template < class R >
inline
typename SegmentC2<R CGAL_CTAG>::Point_2
SegmentC2<R CGAL_CTAG>::source() const
{
  return ptr->e0;
}

template < class R >
inline
typename SegmentC2<R CGAL_CTAG>::Point_2
SegmentC2<R CGAL_CTAG>::target() const
{
  return ptr->e1;
}

template < class R >
CGAL_KERNEL_INLINE
typename SegmentC2<R CGAL_CTAG>::Point_2
SegmentC2<R CGAL_CTAG>::min() const
{
  return lexicographically_xy_smaller(source(),target()) ? source() : target();
}

template < class R >
CGAL_KERNEL_INLINE
typename SegmentC2<R CGAL_CTAG>::Point_2
SegmentC2<R CGAL_CTAG>::max() const
{
  return lexicographically_xy_smaller(source(),target()) ? target() : source();
}

template < class R >
CGAL_KERNEL_INLINE
typename SegmentC2<R CGAL_CTAG>::Point_2
SegmentC2<R CGAL_CTAG>::vertex(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
CGAL_KERNEL_INLINE
typename SegmentC2<R CGAL_CTAG>::Point_2
SegmentC2<R CGAL_CTAG>::point(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
inline
typename SegmentC2<R CGAL_CTAG>::Point_2
SegmentC2<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
CGAL_KERNEL_INLINE
typename SegmentC2<R CGAL_CTAG>::FT
SegmentC2<R CGAL_CTAG>::squared_length() const
{
  return squared_distance(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
typename SegmentC2<R CGAL_CTAG>::Direction_2
SegmentC2<R CGAL_CTAG>::direction() const
{
  return Direction_2( target() - source() );
}

template < class R >
inline
typename SegmentC2<R CGAL_CTAG>::Line_2
SegmentC2<R CGAL_CTAG>::supporting_line() const
{
  return Line_2(*this);
}

template < class R >
inline
SegmentC2<R CGAL_CTAG>
SegmentC2<R CGAL_CTAG>::opposite() const
{
  return SegmentC2<R CGAL_CTAG>(target(), source());
}

template < class R >
inline
SegmentC2<R CGAL_CTAG>
SegmentC2<R CGAL_CTAG>::
transform(const typename SegmentC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{
  return SegmentC2<R CGAL_CTAG>(t.transform(source()), t.transform(target()));
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_2
SegmentC2<R CGAL_CTAG>::bbox() const
{
  return source().bbox() + target().bbox();
}

template < class R >
inline
bool
SegmentC2<R CGAL_CTAG>::is_degenerate() const
{
  return (source() == target());
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentC2<R CGAL_CTAG>::is_horizontal() const
{
  return source().y() == target().y();
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentC2<R CGAL_CTAG>::is_vertical() const
{
  return source().x() == target().x();
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentC2<R CGAL_CTAG>::
has_on(const typename SegmentC2<R CGAL_CTAG>::Point_2 &p) const
{
  return are_ordered_along_line(source(), p, target());
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
SegmentC2<R CGAL_CTAG>::
collinear_has_on(const typename SegmentC2<R CGAL_CTAG>::Point_2 &p) const
{
    return collinear_are_ordered_along_line(source(), p, target());
}

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENTC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const SegmentC2<R CGAL_CTAG> &s)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << s.source() << ' ' << s.target();
    case IO::BINARY :
        return os << s.source() << s.target();
    default:
        return os << "SegmentC2(" << s.source() <<  ", " << s.target() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENTC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENTC2
template < class R >
std::istream &
operator>>(std::istream &is, SegmentC2<R CGAL_CTAG> &s)
{
    typename SegmentC2<R CGAL_CTAG>::Point_2 p, q;

    is >> p >> q;

    s = SegmentC2<R CGAL_CTAG>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENTC2

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_SEGMENT_2_C
