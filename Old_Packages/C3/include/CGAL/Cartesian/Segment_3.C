// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_SEGMENT_3_C
#define CGAL_CARTESIAN_SEGMENT_3_C

#include <CGAL/Cartesian/distance_computations_3.h>

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Twotuple< typename SegmentC3<R CGAL_CTAG>::Point_3 > *
SegmentC3<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple< Point_3 >*)PTR;
}

template < class R >
SegmentC3<R CGAL_CTAG>::
SegmentC3()
{
  PTR = new _Twotuple< Point_3 >;
}

template < class R >
SegmentC3<R CGAL_CTAG>::
SegmentC3(const SegmentC3<R CGAL_CTAG>  &s) :
  Handle((Handle&)s)
{}

template < class R >
SegmentC3<R CGAL_CTAG>::
SegmentC3(const typename SegmentC3<R CGAL_CTAG>::Point_3 &sp,
          const typename SegmentC3<R CGAL_CTAG>::Point_3 &ep)
{
  PTR = new _Twotuple< Point_3 >(sp, ep);
}

template < class R >
inline SegmentC3<R CGAL_CTAG>::~SegmentC3()
{}

template < class R >
SegmentC3<R CGAL_CTAG> &
SegmentC3<R CGAL_CTAG>::operator=(const SegmentC3<R CGAL_CTAG> &s)
{
  Handle::operator=(s);
  return *this;
}

template < class R >
inline
bool
SegmentC3<R CGAL_CTAG>::operator==(const SegmentC3<R CGAL_CTAG> &s) const
{
  return (source() == s.source())  && (target() == s.target());
}

template < class R >
inline
bool
SegmentC3<R CGAL_CTAG>::operator!=(const SegmentC3<R CGAL_CTAG> &s) const
{
  return !(*this == s);
}

template < class R >
long  SegmentC3<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::start() const
{
  return ptr()->e0;
}

template < class R >
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::end() const
{
  return ptr()->e1;
}

template < class R >
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::source() const
{
  return ptr()->e0;
}

template < class R >
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::target() const
{
  return ptr()->e1;
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::min() const
{
  return (lexicographically_xyz_smaller(source(),target())) ? source()
                                                            : target();
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::max() const
{
  return (lexicographically_xyz_smaller(source(),target())) ? target()
                                                            : source();
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::vertex(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::point(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::FT
SegmentC3<R CGAL_CTAG>::squared_length() const
{
  return squared_distance(target(), source());
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Direction_3
SegmentC3<R CGAL_CTAG>::direction() const
{
  return Direction_3( target() - source() );
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Line_3
SegmentC3<R CGAL_CTAG>::supporting_line() const
{
  return Line_3(*this);
}

template < class R >
inline
SegmentC3<R CGAL_CTAG>
SegmentC3<R CGAL_CTAG>::opposite() const
{
  return SegmentC3<R CGAL_CTAG>(target(), source());
}

template < class R >
inline
SegmentC3<R CGAL_CTAG>
SegmentC3<R CGAL_CTAG>::
transform(const typename SegmentC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  return SegmentC3<R CGAL_CTAG>(t.transform(source()), t.transform(target()));
}

template < class R >
inline
bool
SegmentC3<R CGAL_CTAG>::is_degenerate() const
{
  return source() == target();
}

template < class R >
inline
Bbox_3
SegmentC3<R CGAL_CTAG>::bbox() const
{
  return source().bbox() + target().bbox();
}

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENTC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const SegmentC3<R CGAL_CTAG> &s)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << s.source() << ' ' << s.target();
    case IO::BINARY :
        return os << s.source() << s.target();
    default:
        return os << "SegmentC3(" << s.source() <<  ", " << s.target() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENTC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENTC3
template < class R >
std::istream &
operator>>(std::istream &is, SegmentC3<R CGAL_CTAG> &s)
{
    typename SegmentC3<R CGAL_CTAG>::Point_3 p, q;

    is >> p >> q;

    s = SegmentC3<R CGAL_CTAG>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENTC3

template < class R >
bool SegmentC3<R CGAL_CTAG>::
has_on(const typename SegmentC3<R CGAL_CTAG>::Point_3 &p) const
{
  return are_ordered_along_line(source(), p, target());
}

template < class R >
inline
bool
SegmentC3<R CGAL_CTAG>::
collinear_has_on(const typename SegmentC3<R CGAL_CTAG>::Point_3 &p) const
{
  return collinear_are_ordered_along_line(source(), p, target());
}

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_SEGMENT_3_C
