// revision      : 2.8
// revision_date : 28 Oct 1999 
// author(s)     : Hervé Brönnimann

#ifndef CGAL_TETRAHEDRON_D_H
#define CGAL_TETRAHEDRON_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/TetrahedronHd.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Cartesian/Tetrahedron_d.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Plane_d.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class Tetrahedron_d : public _R::Tetrahedron_d_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Tetrahedron_d_base  RTetrahedron_d;

  Tetrahedron_d() : RTetrahedron_d()
  {}
  Tetrahedron_d(const CGAL::Tetrahedron_d<R>& t) : RTetrahedron_d(t)
  {}
  Tetrahedron_d(const RTetrahedron_d&  t) : RTetrahedron_d(t)
  {}
  Tetrahedron_d(const CGAL::Point_d<R>& p,
                     const CGAL::Point_d<R>& q,
                     const CGAL::Point_d<R>& r,
                     const CGAL::Point_d<R>& s)
    : RTetrahedron_d(p,q,r,s)
  {}

  CGAL::Tetrahedron_d<R>&
                     operator=(const CGAL::Tetrahedron_d<R>& t)
                     {
                       RTetrahedron_d::operator=(t);
                       return *this;
                     }
  CGAL::Point_d<R>    vertex(int i) const
                     { return RTetrahedron_d::vertex(i); }
  CGAL::Point_d<R>    operator[](int i) const
                     { return vertex(i); }
  bool               operator==(const CGAL::Tetrahedron_d<R>& t) const
                     { return RTetrahedron_d::operator==(t); }
  bool               operator!=(const CGAL::Tetrahedron_d<R>& t) const
                     { return !(*this == t); }
  int                id() const    /* XXX */
                     { return (int)PTR ; }
  /*
  Bbox_d        bbox() const
                     {
                       return vertex(0).bbox() + vertex(1).bbox()
                            + vertex(2).bbox() + vertex(3).bbox();
                     }
  CGAL::Tetrahedron_d<R>
                     transform(const CGAL::Aff_transformation_d<R>& t) const
                     {
                       return
                       CGAL::Tetrahedron_d<R>(RTetrahedron_d::transform(t));
                     }
  */
  Orientation   orientation() const
                     { return RTetrahedron_d::orientation(); }
  Oriented_side oriented_side(const CGAL::Point_d<R>& p) const
                     { return RTetrahedron_d::oriented_side(p); }
  bool               has_on_positive_side(const CGAL::Point_d<R>& p) const
                     { return oriented_side(p) == ON_POSITIVE_SIDE; }
  bool               has_on_negative_side(const CGAL::Point_d<R>& p) const
                     { return oriented_side(p) == ON_NEGATIVE_SIDE; }
  Bounded_side  bounded_side(const CGAL::Point_d<R>& p) const
                     { return RTetrahedron_d::bounded_side(p); }
  bool               has_on_boundary(const CGAL::Point_d<R>& p) const
                     { return bounded_side(p) == ON_BOUNDARY; }
  bool               has_on_bounded_side(const CGAL::Point_d<R>& p) const
                     { return bounded_side(p) == ON_BOUNDED_SIDE; }
  bool               has_on_unbounded_side(const CGAL::Point_d<R>& p) const
                     { return bounded_side(p) == ON_UNBOUNDED_SIDE; }
  bool               is_degenerate() const
                     { return RTetrahedron_d::is_degenerate(); }
};

#ifndef NO_OSTREAM_INSERT_TETRAHEDRON_D
template < class R >
std::ostream&
operator<<(std::ostream& os, const Tetrahedron_d<R>& t)
{
  typedef typename  R::Tetrahedron_d_base  RTetrahedron_d;
  return os << (const RTetrahedron_d& )t;
}
#endif // NO_OSTREAM_INSERT_TETRAHEDRON_D

#ifndef NO_ISTREAM_EXTRACT_TETRAHEDRON_D
template < class R >
std::istream&
operator>>(std::istream& is, Tetrahedron_d<R>& t)
{
  typedef typename  R::Tetrahedron_d_base  RTetrahedron_d;
  return is >> (RTetrahedron_d& )t;
}
#endif // NO_ISTREAM_EXTRACT_TETRAHEDRON_D

CGAL_END_NAMESPACE

#endif  // CGAL_TETRAHEDRON_D_H
