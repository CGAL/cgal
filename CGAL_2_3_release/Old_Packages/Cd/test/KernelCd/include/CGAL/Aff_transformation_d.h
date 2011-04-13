// revision      : 2.8
// revision_date : 28 Oct 1999 
// author(s)     : Hervé Brönnimann

#ifndef CGAL_AFF_TRANSFORMATION_D_H
#define CGAL_AFF_TRANSFORMATION_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/Aff_transformationHd.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Cartesian/Aff_transformation_d.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Point_d.h>
#include <CGAL/Vector_d.h>
#include <CGAL/Direction_d.h>
#include <CGAL/Plane_d.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class Aff_transformation_d : public _R::Aff_transformation_d_base
{
public:
  typedef _R                                R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Plane_d_base  RPlane_d;
  typedef typename R::Aff_transformation_d_base  RAff_transformation_d;

  // default constructor
  Aff_transformation_d() : RAff_transformation_d()
  {}
  // copy constructor
  Aff_transformation_d(const CGAL::Aff_transformation_d<R>& t)
    : RAff_transformation_d(t)
  {}
  // up cast constructor
  Aff_transformation_d(const RAff_transformation_d&  t)
    : RAff_transformation_d(t)
  {}
  // translation:
  Aff_transformation_d(const Translation tag,
                       const CGAL::Vector_d<R>& v)
    : RAff_transformation_d(tag, v)
  {}
  // scaling:
  Aff_transformation_d(const Scaling tag, int d,
                       const RT& s,
                       const RT& w= RT(1) )
    : RAff_transformation_d(tag, d, s, w)
  {}
  // the general case:
  Aff_transformation_d(
      const RT& m11, const RT& m12, const RT& m13, const RT& m14,
      const RT& m21, const RT& m22, const RT& m23, const RT& m24,
      const RT& m31, const RT& m32, const RT& m33, const RT& m34,
                                                   const RT& w= RT(1) )
  {
    RT m[9] = { m11, m12, m13, m21, m22, m23, m31, m32, m33 };
    RT v[3] = { m14, m24, m34 };
    *this = RAff_transformation_d(3, m+0, m+9, v+0, v+3, w);
  }
  Aff_transformation_d(
      const RT& m11, const RT& m12, const RT& m13,
      const RT& m21, const RT& m22, const RT& m23,
      const RT& m31, const RT& m32, const RT& m33,
                                                   const RT& w = RT(1) )
  {
    RT m[9] = { m11, m12, m13, m21, m22, m23, m31, m32, m33 };
    *this = RAff_transformation_d(3, m+0, m+9, w);
  }
  // dtor
  ~Aff_transformation_d()
  {}
  // transformations
  CGAL::Point_d<R>       transform(const CGAL::Point_d<R>& p) const
                        { return RAff_transformation_d::transform(p); }
  CGAL::Point_d<R>       operator()(const CGAL::Point_d<R>& p) const
                        { return RAff_transformation_d::transform(p); }
  CGAL::Vector_d<R>      transform(const CGAL::Vector_d<R>& v) const
                        { return RAff_transformation_d::transform(v); }
  CGAL::Vector_d<R>      operator()(const CGAL::Vector_d<R>& v) const
                        { return RAff_transformation_d::transform(v); }
  CGAL::Direction_d<R>   transform(const CGAL::Direction_d<R>& d) const
                        { return RAff_transformation_d::transform(d); }
  CGAL::Direction_d<R>   operator()(const CGAL::Direction_d<R>& d) const
                        { return RAff_transformation_d::transform(d); }
  CGAL::Plane_d<R>       transform(const CGAL::Plane_d<R>& pl) const
                        { return RAff_transformation_d::transform(pl); }
  CGAL::Plane_d<R>       operator()(const CGAL::Plane_d<R>& pl) const
                        { return transform(pl); }
  // further members
  CGAL::Aff_transformation_d<R>
                        inverse() const
                        { return RAff_transformation_d::inverse(); }
  int                   dimension() const
                        { return RAff_transformation_d::dimension(); }
  bool                  is_even() const
                        { return RAff_transformation_d::is_even(); }
  bool                  is_odd() const
                        { return !is_even(); }
  // access
  FT                    cartesian(int i, int j) const
                        { return RAff_transformation_d::cartesian(i,j); }
  RT                    homogeneous(int i, int j) const
                        { return RAff_transformation_d::homogeneous(i,j); }
  FT                    m(int i, int j) const
                        { return RAff_transformation_d::m(i,j); }
  RT                    hm(int i, int j) const
                        { return RAff_transformation_d::hm(i,j); }
  // composition
  CGAL::Aff_transformation_d<R>
                        operator*(const CGAL::Aff_transformation_d<R>& t) const
                        {
                          return
                          static_cast<const RAff_transformation_d&>(*this) *
                          static_cast<const RAff_transformation_d&>(t) ;
                        }
};

// I/O operators
#ifndef NO_OSTREAM_INSERT_AFF_TRANSFORMATION_D
template < class R >
std::ostream&
operator<<(std::ostream& os, const CGAL::Aff_transformation_d<R>& t)
{
  typedef typename   R::Aff_transformation_d_base  RAff_transformation_d;
  return os << static_cast<const RAff_transformation_d&>(t);
}
#endif // NO_OSTREAM_INSERT_AFF_TRANSFORMATION_D

#ifndef NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_D
template < class R >
std::istream&
operator>>(std::istream& is, CGAL::Aff_transformation_d<R>& t)
{
  typedef typename   R::Aff_transformation_d_base  RAff_transformation_d;
  return is >> static_cast<const RAff_transformation_d&>(t);
}
#endif // NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_D

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Aff_transformation_d.C>
#endif 

#endif // CGAL_AFF_TRANSFORMATION_D_H
