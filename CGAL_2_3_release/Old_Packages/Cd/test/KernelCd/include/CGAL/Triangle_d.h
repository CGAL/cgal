// revision      : 
// revision_date : 
// author(s)     : Hervé Brönnimann

#ifndef CGAL_TRIANGLE_D_H
#define CGAL_TRIANGLE_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/TriangleH3.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Cartesian/Triangle_d.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Plane_d.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class Triangle_d : public _R::Triangle_d_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Triangle_d_base  RTriangle_d;

  Triangle_d() : RTriangle_d()
  {}
  Triangle_d(const CGAL::Triangle_d<R>& t) : RTriangle_d(t)
  {}
  Triangle_d(const RTriangle_d&  t) : RTriangle_d(t)
  {}
  Triangle_d(const CGAL::Point_d<R>& p,
                  const CGAL::Point_d<R>& q,
                  const CGAL::Point_d<R>& r)
    : RTriangle_d(p,q,r)
  {}

  CGAL::Triangle_d<R>& operator=(const CGAL::Triangle_d<R>& t)
  {
    RTriangle_d::operator=(t);
    return *this;
  }
  bool                operator==(const CGAL::Triangle_d<R>& t) const
                      { return RTriangle_d::operator==(t); }
  bool                operator!=(const CGAL::Triangle_d<R>& t) const
                      { return !(*this == t); }
  int                 id() const   /* XXX */
                      { return (int) PTR ; }
  CGAL::Plane_d<R>     supporting_plane() const
                      {
                        return
                        CGAL::Plane_d<R>(
                            RTriangle_d::supporting_plane());
                      }
  /*
  CGAL::Triangle_d<R>  transform(
                      const CGAL::Aff_transformation_d<R>& t) const
                      {
                        return
                        CGAL::Triangle_d<R>(RTriangle_d::transform( t ));
                      }
  */
  bool                has_on(const CGAL::Point_d<R>& p) const
                      { return RTriangle_d::has_on(p); }
  bool                is_degenerate() const
                      { return RTriangle_d::is_degenerate(); }
  CGAL::Point_d<R>     vertex(int i) const
                      { return RTriangle_d::vertex(i); }
  CGAL::Point_d<R>     operator[](int i) const
                      { return vertex(i); }
  /*
  Bbox_d         bbox() const
                      {
                        return vertex(0).bbox()
                             + vertex(1).bbox()
                             + vertex(2).bbox();
                      }
  */
};

#ifndef NO_OSTREAM_INSERT_TRIANGLE_D
template < class R >
std::ostream&
operator<<(std::ostream& os, const Triangle_d<R>& t)
{
  typedef typename  R::Triangle_d_base  RTriangle_d;
  return os << (const RTriangle_d& )t;
}
#endif // NO_OSTREAM_INSERT_TRIANGLE_D

#ifndef NO_ISTREAM_EXTRACT_TRIANGLE_D
template < class R >
std::istream&
operator>>(std::istream& is, Triangle_d<R>& t)
{
  typedef typename  R::Triangle_d_base  RTriangle_d;
  return is >> (RTriangle_d& )t;
}
#endif // NO_ISTREAM_EXTRACT_TRIANGLE_D

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGLE_D_H
