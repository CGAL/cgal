// revision      : 2.8
// revision_date : 28 Oct 1999 
// author(s)     : Hervé Brönnimann

#ifndef CGAL_RAY_D_H
#define CGAL_RAY_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/RayH3.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Cartesian/Ray_d.h>
#endif // CGAL_CARTESIAN_H

CGAL_BEGIN_NAMESPACE

template <class _R>
class Ray_d : public _R::Ray_d_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Ray_d_base  RRay_d;

  Ray_d() : RRay_d()
  {}
  Ray_d(const CGAL::Ray_d<R>& r) : RRay_d(r)
  {}
  Ray_d(const RRay_d&  r) : RRay_d(r)
  {}
  Ray_d(const CGAL::Point_d<R>& sp,
            const CGAL::Point_d<R>& secondp)
    : RRay_d(sp, secondp)
  {}
  Ray_d(const CGAL::Point_d<R>& sp,
            const CGAL::Direction_d<R>& d)
    : RRay_d(sp, d)
  {}

  CGAL::Ray_d<R>&      operator=(const CGAL::Ray_d<R>& r)
  {
      RRay_d::operator=(r);
      return *this;
  }
  bool                operator==(const CGAL::Ray_d<R>& r) const
  { return RRay_d::operator==(r); }
  bool                operator!=(const CGAL::Ray_d<R>& r) const
  { return !(*this == r); }

  int                 id() const  /* XXX */
  { return (int)  PTR ; }

  CGAL::Point_d<R>     start() const
  { return RRay_d::start(); }
  CGAL::Point_d<R>     source() const
  { return RRay_d::source(); }
  CGAL::Point_d<R>     second_point() const
  { return RRay_d::second_point(); }
  CGAL::Point_d<R>     point(int i) const
  { return RRay_d::point(i); }
  CGAL::Direction_d<R> direction() const
  { return RRay_d::direction(); }
  CGAL::Line_d<R>      supporting_line() const
  { return RRay_d::supporting_line(); }
  CGAL::Ray_d<R>       opposite() const
  { return RRay_d::opposite(); }
  // CGAL::Ray_d<R> transform(const CGAL::Aff_transformation_d<R>& t) const
  // { return RRay_d::transform(t); }
  bool                is_degenerate() const
  { return RRay_d::is_degenerate(); }
  bool                has_on(const CGAL::Point_d<R>& p) const
  { return RRay_d::has_on(p); }
};

#ifndef NO_OSTREAM_INSERT_RAY_D
template < class R >
std::ostream&
operator<<(std::ostream& os, const Ray_d<R>& r)
{
  typedef typename  R::Ray_d_base  RRay_d;
  return os << (const RRay_d& )r;
}
#endif // NO_OSTREAM_INSERT_RAY_D

#ifndef NO_ISTREAM_EXTRACT_RAY_D
template < class R >
std::istream&
operator>>(std::istream& is, Ray_d<R>& r)
{
  typedef typename  R::Ray_d_base  RRay_d;
  return is >> (RRay_d& )r;
}
#endif // NO_ISTREAM_EXTRACT_RAY_D


CGAL_END_NAMESPACE

#include <CGAL/Line_d.h>

#endif // CGAL_RAY_D_H
