// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Aff_transformation_2.h
// source        : include/CGAL/Cartesian/Aff_transformation_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr, kettner@inf.fu-berlin.de
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#include <CGAL/Cartesian/redefine_names_2.h>
#endif

#ifndef CGAL_PROTECT_CMATH
#define CGAL_PROTECT_CMATH
#include <cmath>
#endif

#ifndef CGAL_CONFIG_H
#include <CGAL/config.h>
#endif // CGAL_CONFIG_H
#ifndef CGAL_HANDLE_H
#include <CGAL/Handle.h>
#endif // CGAL_HANDLE_H

CGAL_BEGIN_NAMESPACE

class Identity_transformation;
template < class R > class Aff_transformation_rep_baseC2;
template < class R > class Aff_transformation_repC2;
template < class R > class Translation_repC2;
template < class R > class Rotation_repC2;
template < class R > class Scaling_repC2;

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_2_H
#include <CGAL/Cartesian/Aff_transformation_rep_2.h>
#endif

#ifndef CGAL_CARTESIAN_TRANSLATION_REP_2_H
#include <CGAL/Cartesian/Translation_rep_2.h>
#endif

#ifndef CGAL_CARTESIAN_ROTATION_REP_2_H
#include <CGAL/Cartesian/Rotation_rep_2.h>
#endif

#ifndef CGAL_CARTESIAN_SCALING_REP_2_H
#include <CGAL/Cartesian/Scaling_rep_2.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class _R >
class Aff_transformationC2
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
public:
  typedef _R                               R;
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
   
  Aff_transformationC2();
  Aff_transformationC2(const Self &t);

  // Identity
  Aff_transformationC2(const Identity_transformation);

  // Translation:
  Aff_transformationC2(const Translation,
                       const Vector_2 &v);

  // Rational Rotation:
  Aff_transformationC2(const Rotation,
                       const Direction_2 &d,
                       const FT &num,
                       const FT &den = FT(1));

  Aff_transformationC2(const Rotation,
                       const FT &sine_rho,
                       const FT &cosine_rho,
                       const FT &hw = FT(1));

  // Scaling:
  Aff_transformationC2(const Scaling,
                       const FT &s,
                       const FT &w = FT(1));

  // The general case:
  Aff_transformationC2(const FT & m11,
                       const FT & m12,
                       const FT & m13,
                       const FT & m21,
                       const FT & m22,
                       const FT & m23,
                       const FT &w = FT(1));

  Aff_transformationC2(const FT & m11, const FT & m12,
                       const FT & m21, const FT & m22,
                       const FT &w = FT(1));

  ~Aff_transformationC2();

  Self &operator=(const Self &t);

  Point_2     transform(const Point_2 &p) const { return ptr()->transform(p); } 
  Point_2     operator()(const Point_2 &p) const { return transform(p); }

  Vector_2    transform(const Vector_2 &v) const { return ptr()->transform(v); }
  Vector_2    operator()(const Vector_2 &v) const { return transform(p); }

  Direction_2 transform(const Direction_2 &d) const { return ptr()->transform(d); }
  Direction_2 operator()(const Direction_2 &d) const { return transform(d); }

  Line_2      transform(const Line_2 &l) const { return l.transform(*this); }
  Line_2      operator()(const Line_2 &l) const {return transform(l); }

  Self        inverse() const { return ptr()->inverse(); }


  bool        is_even() const { return ptr()->is_even(); }
  bool        is_odd() const { return ! (ptr()->is_even()); }

  FT          cartesian(int i, int j) const { return ptr()->cartesian(i,j); }
  FT          homogeneous(int i, int j) const { return cartesian(i,j); }
  FT          m(int i, int j) const { return cartesian(i,j); }
  FT          hm(int i, int j) const { return cartesian(i,j); }

  Self operator*(const Self &t) const
  {
    return (*ptr()) * (*t.ptr());
  }

  std::ostream &print(std::ostream &os) const;

private:
  Aff_t_base* ptr() const { return  (Aff_t_base*)PTR; }
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_2_C
#include <CGAL/Cartesian/Aff_transformation_2.C>
#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_2_C
#endif

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H
