// revision      : 2.8
// revision_date : 28 Oct 1999 
// author(s)     : Hervé Brönnimann

#ifndef CGAL_LINE_D_H
#define CGAL_LINE_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/LineHd.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Cartesian/Line_d.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Segment_d.h>
#include <CGAL/Point_d.h>
#include <CGAL/Ray_d.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class Line_d : public _R::Line_d_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Line_d_base           RLine_d;

  Line_d() : RLine_d()
  {}
  Line_d(const CGAL::Line_d<R>  & l) : RLine_d( ( const RLine_d&  )l)
  {}
  Line_d(const CGAL::Point_d<R> & p,
              const CGAL::Point_d<R> & q) : RLine_d(p,q)
  {}
  // conversion impl -> interface class
  Line_d(const RLine_d&  l) : RLine_d(l)
  {}
  Line_d(const CGAL::Segment_d<R> & s) : RLine_d( s )
  {}
  Line_d(const CGAL::Ray_d<R> & r) : RLine_d( r )
  {}
  Line_d(const CGAL::Point_d<R> & p,
              const CGAL::Direction_d<R> & d) : RLine_d( p, d )
  {}

  CGAL::Line_d<R>&     operator=(const CGAL::Line_d<R> & l)
  {
    RLine_d::operator=(l);
    return *this;
  }

  bool                operator==(const CGAL::Line_d<R> & l) const
  { return RLine_d::operator==(l); }

  bool                operator!=(const CGAL::Line_d<R> & l) const
  { return !(*this == l); }

  int                 id() const    /* XXX */
  { return (int) PTR; }

  CGAL::Plane_d<R>     perpendicular_plane(const CGAL::Point_d<R> & p) const
  { return RLine_d::perpendicular_plane(p); }

  CGAL::Line_d<R>      opposite() const
  { return RLine_d::opposite(); }

  CGAL::Point_d<R>     point() const
  { return RLine_d::point(); }

  CGAL::Point_d<R>     point(int i) const
  { return RLine_d::point(i); }

  CGAL::Point_d<R>     projection(const CGAL::Point_d<R>& p) const
  { return RLine_d::projection(p); }

  CGAL::Direction_d<R> direction() const
  { return RLine_d::direction(); }

  bool                has_on(const CGAL::Point_d<R>& p) const
  { return RLine_d::has_on(p); }

  bool                is_degenerate() const
  { return RLine_d::is_degenerate(); }

  // CGAL::Line_d<R> transform(const CGAL::Aff_transformation_d<R> & t) const
  // { return RLine_d::transform(t); }
};

#ifndef NO_OSTREAM_INSERT_LINE_D
template < class R >
std::ostream&
operator<<(std::ostream& os, const Line_d<R>& l)
{
  typedef typename  R::Line_d_base  RLine_d;
  return os << (const RLine_d& )l;
}
#endif // NO_OSTREAM_INSERT_LINE_D

#ifndef NO_ISTREAM_EXTRACT_LINE_D
template < class R >
std::istream&
operator>>(std::istream & is, Line_d<R> & p)
{
  typedef typename  R::Line_d_base  RLine_d;
  is >> ( RLine_d&  )p;
  return is;
}
#endif // NO_ISTREAM_EXTRACT_LINE_D

CGAL_END_NAMESPACE


#ifndef CGAL_PLANE_D_H
#include <CGAL/Plane_3.h>
#endif // CGAL_PLANE_3_H

#endif // CGAL_LINE_3_H
