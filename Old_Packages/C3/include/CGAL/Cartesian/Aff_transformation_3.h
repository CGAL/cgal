// ==========================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Aff_transformation_3.h
// source        : include/CGAL/Cartesian/Aff_transformation_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_3_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

#include <cmath>

CGAL_BEGIN_NAMESPACE

class Identity_transformation;
template <class R> class Aff_transformation_rep_baseC3;
template <class R> class Aff_transformation_repC3;
template <class R> class Translation_repC3;
template <class R> class Scaling_repC3;

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_H
#include <CGAL/Cartesian/Aff_transformation_rep_3.h>
#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_H

#ifndef CGAL_CARTESIAN_CARTESIAN_TRANSLATION_REP_3_H
#include <CGAL/Cartesian/Translation_rep_3.h>
#endif

#ifndef CGAL_CARTESIAN_CARTESIAN_SCALING_REP_3_H
#include <CGAL/Cartesian/Scaling_rep_3.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class _R >
class Aff_transformationC3
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
  friend class PlaneC3<_R CGAL_CTAG>;

public:
  typedef _R                               R;
  typedef typename R::FT                   FT;
  typedef typename R::FT                   RT;
  typedef Aff_transformation_rep_baseC3<R> Aff_t_base;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef Aff_transformationC3<R,Cartesian_tag> Self;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Vector_3             Vector_3;
  typedef typename R::Direction_3          Direction_3;
  typedef typename R::Plane_3              Plane_3;
#else
  typedef Aff_transformationC3<R>          Self;
  typedef typename R::Point_3_base         Point_3;
  typedef typename R::Vector_3_base        Vector_3;
  typedef typename R::Direction_3_base     Direction_3;
  typedef typename R::Plane_3_base         Plane_3;
#endif

  Aff_transformationC3();
  // Aff_transformationC3(const Self &t); // Provided by default

  // Identity constructor:
  Aff_transformationC3(const Identity_transformation &);

  // Translation:
  Aff_transformationC3(const Translation,
                       const Vector_3 &v);

  // Scaling:
  Aff_transformationC3(const Scaling,
                       const FT &s,
                       const FT &w = FT(1));

  // General form: without translation
  Aff_transformationC3(const FT& m11, const FT& m12, const FT& m13,
                       const FT& m21, const FT& m22, const FT& m23,
                       const FT& m31, const FT& m32, const FT& m33,
                       const FT& w= FT(1));

  // General form: with translation
  Aff_transformationC3(const FT& m11, const FT& m12, const FT& m13,
                                                           const FT& m14,
                       const FT& m21, const FT& m22, const FT& m23,
                                                           const FT& m24,
                       const FT& m31, const FT& m32, const FT& m33,
                                                           const FT& m34,
                       const FT& w = FT(1));

  ~Aff_transformationC3();

  Point_3     transform(const Point_3 &p) const { return ptr()->transform(p); }
  Point_3     operator()(const Point_3 &p) const { return transform(p); }

  Vector_3    transform(const Vector_3 &v) const { return ptr()->transform(v); }
  Vector_3    operator()(const Vector_3 &v) const { return transform(v); }

  Direction_3 transform(const Direction_3 &d) const
                                              { return ptr()->transform(d); }
  Direction_3 operator()(const Direction_3 &d) const { return transform(d); }

  Plane_3     transform(const Plane_3& p) const { return p.transform(*this); }
  Plane_3     operator()(const Plane_3& p) const { return transform(l); }

  Self        inverse() const { return ptr()->inverse(); }
  
  bool        is_even() const { return ptr()->is_even(); }
  bool        is_odd() const { return  ! (ptr()->is_even()); }
  
  FT          cartesian(int i, int j) const { return ptr()->cartesian(i,j); }
  FT          homogeneous(int i, int j) const { return cartesian(i,j); }
  FT          m(int i, int j) const { return cartesian(i,j); }
  FT          hm(int i, int j) const { return cartesian(i,j); }

  Self operator*(const Self &t) const { return (*ptr()) * (*t.ptr()); }

protected:
  Self        transpose() const { return ptr()->transpose(); }

private:
  Aff_t_base*       ptr() const { return  (Aff_t_base*)PTR; }
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_3_C
#include <CGAL/Cartesian/Aff_transformation_3.C>
#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_3_C
#endif 

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_3_H
