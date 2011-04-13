// revision      : 2.8
// revision_date : 28 Oct 1999 
// author(s)     : Hervé Brönnimann

#ifndef CGAL_SEGMENT_D_H
#define CGAL_SEGMENT_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/SegmentH3.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Cartesian/Segment_d.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Line_d.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class Segment_d : public _R::Segment_d_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Segment_d_base  RSegment_d;

  Segment_d() : RSegment_d()
  {}
  Segment_d(const CGAL::Segment_d<R>& s) : RSegment_d(s)
  {}
  Segment_d(const CGAL::Point_d<R>& sp, const CGAL::Point_d<R>& ep)
    : RSegment_d(sp,ep)
  {}
  Segment_d(const RSegment_d&  s) : RSegment_d(s)
  {}

  CGAL::Segment_d<R>&   operator=(const CGAL::Segment_d<R>& s)
  {
    RSegment_d::operator=(s);
    return *this;
  }
  bool                 has_on(const CGAL::Point_d<R>& p) const
  { return RSegment_d::has_on(p); }
  bool                 operator==(const CGAL::Segment_d<R>& s) const
  { return RSegment_d::operator==(s); }
  bool                 operator!=(const CGAL::Segment_d<R>& s) const
  { return !(*this == s); }
  int                  id() const   /* XXX */
  { return (int) PTR ; }
  CGAL::Point_d<R>     start() const
  { return RSegment_d::start(); }
  CGAL::Point_d<R>     end() const
  { return RSegment_d::end(); }
  CGAL::Point_d<R>     source() const
  { return RSegment_d::source(); }
  CGAL::Point_d<R>     target() const
  { return RSegment_d::target(); }
  CGAL::Point_d<R>     min() const
  { return RSegment_d::min(); }
  CGAL::Point_d<R>     max() const
  { return RSegment_d::max(); }
  CGAL::Point_d<R>     vertex(int i) const
  { return RSegment_d::vertex(i); }
  CGAL::Point_d<R>     operator[](int i) const
  { return vertex(i); }
  FT                   squared_length() const
  { return RSegment_d::squared_length(); }
  CGAL::Direction_d<R> direction() const
  { return RSegment_d::direction(); }
  CGAL::Segment_d<R>  opposite() const
  { return CGAL::Segment_d<R>(target(),source()); }
  // CGAL::Segment_d<R>  transform(const CGAL::Aff_transformation_d<R>& t) const
  // { return RSegment_d::transform(t); }
  CGAL::Line_d<R>     supporting_line() const
  { return RSegment_d::supporting_line(); }
  bool                is_degenerate() const
  { return RSegment_d::is_degenerate(); }
  // Bbox_d         bbox() const
  // { return source().bbox() + target().bbox(); }
};


#ifndef NO_OSTREAM_INSERT_SEGMENT_D
template < class R>
std::ostream&
operator<<(std::ostream& os, const Segment_d<R>& s)
{
  typedef typename  R::Segment_d_base  RSegment_d;
  return os << (const RSegment_d& )s;
}
#endif // NO_OSTREAM_INSERT_SEGMENT_D

#ifndef NO_ISTREAM_EXTRACT_SEGMENT_D
template < class R>
std::istream&
operator>>(std::istream& is, Segment_d<R>& s)
{
  typedef typename  R::Segment_d_base  RSegment_d;
  return is >> (RSegment_d& )s;
}
#endif // NO_ISTREAM_EXTRACT_SEGMENT_D


CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_D_H
