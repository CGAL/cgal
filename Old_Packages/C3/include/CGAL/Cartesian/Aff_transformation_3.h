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

CGAL_BEGIN_NAMESPACE

class Identity;
template <class R> class Aff_transformation_rep_baseC3;
template <class R> class Aff_transformation_repC3;
template <class R> class Translation_repC3;
template <class R> class Scaling_repC3;

template <class R>
Aff_transformationC3<R CGAL_CTAG>
_general_transformation_composition (const Aff_transformation_rep_baseC3<R> &l,
                                     const Aff_transformation_rep_baseC3<R> &r );

template <class R>
Aff_transformationC3<R CGAL_CTAG>
operator*(const Aff_transformationC3<R CGAL_CTAG> &a,
          const Aff_transformationC3<R CGAL_CTAG> &b );


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

template <class R>
std::ostream &operator<<(std::ostream &os,
                         const Aff_transformationC3<R CGAL_CTAG> &t);

template < class _R >
class Aff_transformationC3
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
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

  friend std::ostream &operator<< CGAL_NULL_TMPL_ARGS(std::ostream &os,
                                 const Self &t);

  friend Self
  _general_transformation_composition CGAL_NULL_TMPL_ARGS(
                                 const Aff_t_base &l,
                                 const Aff_t_base &r );

  friend Self operator* CGAL_NULL_TMPL_ARGS(
                                 const Self &a,
                                 const Self &b );

  // default constructor:
  Aff_transformationC3();

  // Identity constructor:
  Aff_transformationC3(const Identity &);

  // Translation:
  Aff_transformationC3(const Translation,
                       const Vector_3 &v);

  // Scaling:
  Aff_transformationC3(const Scaling,
                       const FT &s,
                       const FT &w = FT(1));

  // General form:
  Aff_transformationC3(const FT& m11, const FT& m12,
                       const FT& m13,
                       const FT& m21,
                       const FT& m22,
                       const FT& m23,
                       const FT& m31,
                       const FT& m32,
                       const FT& m33,
                       const FT& w= FT(1));

  Aff_transformationC3(const FT& m11, const FT& m12,
                       const FT& m13, const FT& m14,
                       const FT& m21, const FT& m22,
                       const FT& m23, const FT& m24,
                       const FT& m31, const FT& m32,
                       const FT& m33, const FT& m34,
                       const FT& w = FT(1));

  ~Aff_transformationC3();

  Point_3     transform(const Point_3 &p) const;
  Point_3     operator()(const Point_3 &p) const;

  Vector_3    transform(const Vector_3 &v) const;
  Vector_3    operator()(const Vector_3 &v) const;

  Direction_3 transform(const Direction_3 &d) const;
  Direction_3 operator()(const Direction_3 &d) const;

  Plane_3     transform(const Plane_3& p) const;
  Plane_3     operator()(const Plane_3& p) const;

  Self        inverse() const;
  Self        transpose() const;
  Self        general_form() const;
  bool        is_even() const;
  bool        is_odd() const;
  FT          cartesian(int i, int j) const { return ptr()->cartesian(i,j); }
  FT          homogeneous(int i, int j) const { return cartesian(i,j); }
  FT          m(int i, int j) const { return cartesian(i,j); }
  FT          hm(int i, int j) const { return cartesian(i,j); }

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
