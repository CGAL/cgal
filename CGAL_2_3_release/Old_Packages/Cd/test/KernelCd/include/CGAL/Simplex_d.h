// revision      : 2.8
// revision_date : 28 Oct 1999 
// author(s)     : Hervé Brönnimann

#ifndef CGAL_SIMPLEX_D_H
#define CGAL_SIMPLEX_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/SimplexHd.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Cartesian/Simplex_d.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Plane_d.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class Simplex_d : public _R::Simplex_d_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Simplex_d_base        RSimplex_d;
  typedef typename CGAL::Point_d<R>         Point_d;

  Simplex_d() : RSimplex_d()
  {}
  Simplex_d(const CGAL::Simplex_d<R>& t) : RSimplex_d(t)
  {}
  Simplex_d(const RSimplex_d&  t) : RSimplex_d(t)
  {}
  Simplex_d(const Point_d& p, const Point_d& q, const Point_d& r)
  {
    Point_d v[3] = { p, q, r };
    *this = RSimplex_d(v+0,v+3);
  }
  Simplex_d(const Point_d& p, const Point_d& q, const Point_d& r,
            const Point_d& s)
  {
    Point_d v[4] = { p, q, r, s };
    *this = RSimplex_d(v+0,v+4);
  }
  template < class PointIterator >
  Simplex_d(const PointIterator &first, const PointIterator &last)
    : RSimplex_d(first,last)
  {}

  CGAL::Simplex_d<R>&
                     operator=(const CGAL::Simplex_d<R>& t)
                     {
                       RSimplex_d::operator=(t);
                       return *this;
                     }
  Point_d   vertex(int i) const
                     { return RSimplex_d::vertex(i); }
  Point_d   operator[](int i) const
                     { return vertex(i); }
  bool               operator==(const CGAL::Simplex_d<R>& t) const
                     { return RSimplex_d::operator==(t); }
  bool               operator!=(const CGAL::Simplex_d<R>& t) const
                     { return !(*this == t); }
  int                id() const    /* XXX */
                     { return (int)PTR ; }
  /*
  Bbox_d        bbox() const
                     {
                       return vertex(0).bbox() + vertex(1).bbox()
                            + vertex(2).bbox() + vertex(3).bbox();
                     }
  CGAL::Simplex_d<R>
                     transform(const CGAL::Aff_transformation_d<R>& t) const
                     {
                       return
                       CGAL::Simplex_d<R>(RSimplex_d::transform(t));
                     }
  */
  Orientation   orientation() const
                     { return RSimplex_d::orientation(); }
  Oriented_side oriented_side(const Point_d& p) const
                     { return RSimplex_d::oriented_side(p); }
  bool               has_on_positive_side(const Point_d& p) const
                     { return oriented_side(p) == ON_POSITIVE_SIDE; }
  bool               has_on_negative_side(const Point_d& p) const
                     { return oriented_side(p) == ON_NEGATIVE_SIDE; }
  Bounded_side  bounded_side(const Point_d& p) const
                     { return RSimplex_d::bounded_side(p); }
  bool               has_on_boundary(const Point_d& p) const
                     { return bounded_side(p) == ON_BOUNDARY; }
  bool               has_on_bounded_side(const Point_d& p) const
                     { return bounded_side(p) == ON_BOUNDED_SIDE; }
  bool               has_on_unbounded_side(const Point_d& p) const
                     { return bounded_side(p) == ON_UNBOUNDED_SIDE; }
  bool               is_degenerate() const
                     { return RSimplex_d::is_degenerate(); }
};

#ifndef NO_OSTREAM_INSERT_SIMPLEX_D
template < class R >
std::ostream&
operator<<(std::ostream& os, const Simplex_d<R>& t)
{
  typedef typename  R::Simplex_d_base  RSimplex_d;
  return os << (const RSimplex_d& )t;
}
#endif // NO_OSTREAM_INSERT_SIMPLEX_D

#ifndef NO_ISTREAM_EXTRACT_SIMPLEX_D
template < class R >
std::istream&
operator>>(std::istream& is, Simplex_d<R>& t)
{
  typedef typename  R::Simplex_d_base  RSimplex_d;
  return is >> (RSimplex_d& )t;
}
#endif // NO_ISTREAM_EXTRACT_SIMPLEX_D

CGAL_END_NAMESPACE

#endif // CGAL_SIMPLEX_D_H
