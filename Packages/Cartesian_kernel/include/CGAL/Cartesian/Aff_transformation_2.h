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
// file          : include/CGAL/Cartesian/Aff_transformation_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Lutz Kettner
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H

#include <cmath>
#include <CGAL/Cartesian/redefine_names_2.h>

CGAL_BEGIN_NAMESPACE

class Identity_transformation;
template < class R > class Aff_transformation_rep_baseC2;
template < class R > class Aff_transformation_repC2;
template < class R > class Translation_repC2;
template < class R > class Rotation_repC2;
template < class R > class Scaling_repC2;

CGAL_END_NAMESPACE

#include <CGAL/Cartesian/Aff_transformation_rep_2.h>
#include <CGAL/Cartesian/Translation_rep_2.h>
#include <CGAL/Cartesian/Rotation_rep_2.h>
#include <CGAL/Cartesian/Scaling_rep_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class Aff_transformationC2 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Aff_transformation_handle_2
{
public:
  typedef R_                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
  typedef Aff_transformation_rep_baseC2<R> Aff_t_base;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef Aff_transformationC2<R,Cartesian_tag> Self;
  typedef typename R::Point_2              Point_2;
  typedef typename R::Vector_2             Vector_2;
  typedef typename R::Direction_2          Direction_2;
  typedef typename R::Line_2               Line_2;
#else
  typedef Aff_transformationC2<R>          Self;
  typedef typename R::Point_2_base         Point_2;
  typedef typename R::Vector_2_base        Vector_2;
  typedef typename R::Direction_2_base     Direction_2;
  typedef typename R::Line_2_base          Line_2;
#endif // CGAL_CFG_NO_ADVANCED_KERNEL
   
  Aff_transformationC2()
  {
    PTR = new Aff_transformation_repC2<R>(FT(1), FT(0), FT(0), FT(1));
  }

  Aff_transformationC2(const Identity_transformation)
  {
    PTR = new Aff_transformation_repC2<R>(FT(1), FT(0), FT(0), FT(1));
  }

  Aff_transformationC2(const Translation, const Vector_2 &v)
  {
    PTR = new Translation_repC2<R>(v);
  }

  // Rational Rotation:
  Aff_transformationC2(const Rotation,
                       const Direction_2 &d,
                       const FT &num,
                       const FT &den = FT(1))
  {
    PTR = new Rotation_repC2<R>(d, num, den);
  }

  Aff_transformationC2(const Rotation,
                       const FT &sine,
                       const FT &cosine,
                       const FT &w = FT(1))
  {
    if (w != FT(1))
      PTR = new Rotation_repC2<R>(sine/w, cosine/w);
    else
      PTR = new Rotation_repC2<R>(sine, cosine);
  }

  Aff_transformationC2(const Scaling, const FT &s, const FT &w = FT(1))
  {
    if (w != FT(1))
      PTR = new Scaling_repC2<R>(s/w);
    else
      PTR = new Scaling_repC2<R>(s);
  }

  // The general case:
  // a 3x2 matrix for the operations combining rotation, scaling, translation
  Aff_transformationC2(const FT & m11, const FT & m12, const FT & m13,
                       const FT & m21, const FT & m22, const FT & m23,
                       const FT &w = FT(1))
  {
    if (w != FT(1))
      PTR = new Aff_transformation_repC2<R>(m11/w, m12/w, m13/w,
                                            m21/w, m22/w, m23/w);
    else
      PTR = new Aff_transformation_repC2<R>(m11, m12, m13,
                                            m21, m22, m23);
  }

  Aff_transformationC2(const FT & m11, const FT & m12,
                       const FT & m21, const FT & m22,
                       const FT &w = FT(1))
  {
    PTR = new Aff_transformation_repC2<R>(m11/w, m12/w, m21/w, m22/w);
  }

  Point_2
  transform(const Point_2 &p) const 
  { return ptr()->transform(p); } 

  Point_2
  operator()(const Point_2 &p) const
  { return transform(p); }

  Vector_2
  transform(const Vector_2 &v) const 
  { return ptr()->transform(v); }

  Vector_2
  operator()(const Vector_2 &v) const
  { return transform(v); } // FIXME : not compiled by the test-suite.

  Direction_2
  transform(const Direction_2 &d) const
  { return ptr()->transform(d); }

  Direction_2
  operator()(const Direction_2 &d) const
  { return transform(d); }

  Line_2
  transform(const Line_2 &l) const
  { return l.transform(*this); }

  Line_2
  operator()(const Line_2 &l) const
  { return transform(l); }

  Self inverse() const { return ptr()->inverse(); }

  bool is_even() const { return ptr()->is_even(); }
  bool is_odd() const { return ! (ptr()->is_even()); }

  FT cartesian(int i, int j) const { return ptr()->cartesian(i,j); }
  FT homogeneous(int i, int j) const { return cartesian(i,j); }
  FT m(int i, int j) const { return cartesian(i,j); }
  FT hm(int i, int j) const { return cartesian(i,j); }

  Self operator*(const Self &t) const
  {
    return (*ptr()) * (*t.ptr());
  }

  std::ostream &
  print(std::ostream &os) const;

private:
  Aff_t_base* ptr() const { return static_cast<Aff_t_base*>(PTR); }
  // FIXME : ptr() should be in Handle.
};

template < class R >
std::ostream&
Aff_transformationC2<R CGAL_CTAG>::print(std::ostream &os) const
{
  ptr()->print(os);
  return os;
}

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC2
template < class R >
std::ostream&
operator<<(std::ostream& os, const Aff_transformationC2<R CGAL_CTAG>& t)
{
  t.print(os);
  return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONC2
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H
