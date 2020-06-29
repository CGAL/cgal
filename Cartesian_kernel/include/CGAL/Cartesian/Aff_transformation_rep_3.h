// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_H

#include <ostream>
#include <CGAL/determinant.h>

namespace CGAL {

template < class R >
class Aff_transformation_rep_baseC3
  : public Ref_counted_virtual
{
public:
  typedef typename R::FT                   FT;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Vector_3             Vector_3;
  typedef typename R::Direction_3          Direction_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;

  virtual ~Aff_transformation_rep_baseC3(){}

  virtual Point_3     transform(const Point_3 &p) const = 0;
  virtual Vector_3    transform(const Vector_3 &v) const = 0;
  virtual Direction_3 transform(const Direction_3 &d) const = 0;

  virtual Aff_transformation_3 operator*(
                       const Aff_transformation_rep_baseC3 &t) const = 0;

  virtual Aff_transformation_3 compose(
                       const Translation_repC3<R> &t) const  = 0;

  virtual Aff_transformation_3 compose(
                       const Scaling_repC3<R> &t) const  = 0;

  virtual Aff_transformation_3 compose(
                       const Aff_transformation_repC3<R> &t) const  = 0;

  virtual Aff_transformation_3 inverse() const = 0;
  virtual Aff_transformation_3 transpose() const = 0;
  virtual bool                 is_even() const = 0;
  virtual FT                   cartesian(int i, int j) const = 0;
  virtual std::ostream         &print(std::ostream &os) const = 0;
};

template < class R >
class Aff_transformation_repC3
  : public Aff_transformation_rep_baseC3<R>
{
  friend class Translation_repC3<R>;
  friend class Scaling_repC3<R>;

public:
  typedef typename R::FT                                FT;
  typedef Aff_transformation_repC3<R>                   Self;
  typedef Aff_transformation_rep_baseC3<R>              Transformation_base_3;
  typedef Aff_transformation_repC3<R>                   Transformation_3;
  typedef Translation_repC3<R>                          Translation_3;
  typedef Scaling_repC3<R>                              Scaling_3;
  typedef typename Transformation_base_3::Point_3       Point_3;
  typedef typename Transformation_base_3::Vector_3      Vector_3;
  typedef typename Transformation_base_3::Direction_3   Direction_3;
  typedef typename Transformation_base_3::
                                   Aff_transformation_3 Aff_transformation_3;

  Aff_transformation_repC3()
  {}

  Aff_transformation_repC3(const FT& m11, const FT& m12, const FT& m13,
                           const FT& m21, const FT& m22, const FT& m23,
                           const FT& m31, const FT& m32, const FT& m33)
    : t11(m11), t12(m12), t13(m13), t14(FT(0)),
      t21(m21), t22(m22), t23(m23), t24(FT(0)),
      t31(m31), t32(m32), t33(m33), t34(FT(0))
  {}

  Aff_transformation_repC3(
     const FT& m11, const FT& m12, const FT& m13, const FT& m14,
     const FT& m21, const FT& m22, const FT& m23, const FT& m24,
     const FT& m31, const FT& m32, const FT& m33, const FT& m34
  )
    : t11(m11), t12(m12), t13(m13), t14(m14),
      t21(m21), t22(m22), t23(m23), t24(m24),
      t31(m31), t32(m32), t33(m33), t34(m34)
  {}

  virtual ~Aff_transformation_repC3()
  {}

  virtual Point_3 transform(const Point_3& p) const // FIXME : construction
  {
    typename R::Construct_point_3 construct_point_3;
    return construct_point_3(t11 * p.x() + t12 * p.y() + t13 * p.z() + t14,
                             t21 * p.x() + t22 * p.y() + t23 * p.z() + t24,
                             t31 * p.x() + t32 * p.y() + t33 * p.z() + t34);
  }

  // note that a vector is not translated
  virtual Vector_3 transform(const Vector_3& v) const // FIXME : construction
  {
    return Vector_3(t11 * v.x() + t12 * v.y() + t13 * v.z(),
                    t21 * v.x() + t22 * v.y() + t23 * v.z(),
                    t31 * v.x() + t32 * v.y() + t33 * v.z());
  }

  // note that a direction is not translated
  virtual Direction_3 transform(const Direction_3& dir) const
  { // FIXME : construction
    Vector_3 v = dir.to_vector();
    return Direction_3(t11 * v.x() + t12 * v.y() + t13 * v.z(),
                       t21 * v.x() + t22 * v.y() + t23 * v.z(),
                       t31 * v.x() + t32 * v.y() + t33 * v.z());
  }

  // Note that Aff_transformation is not defined yet,
  // so the following 6 functions have to be implemented
  // outside class body
  virtual Aff_transformation_3 inverse() const;
  virtual Aff_transformation_3 transpose() const;
  virtual Aff_transformation_3 operator*(const Transformation_base_3 &t) const;
  virtual Aff_transformation_3 compose(const Transformation_3 &t) const;
  virtual Aff_transformation_3 compose(const Translation_3 &t) const;
  virtual Aff_transformation_3 compose(const Scaling_3 &t) const;

  virtual bool is_even() const
  {
    return sign_of_determinant(t11, t12, t13,
                                  t21, t22, t23,
                                  t31, t32, t33) == POSITIVE;
  }

  virtual FT cartesian(int i, int j) const
  {
    switch (i)
    {
    case 0: switch (j)
            {
              case 0: return t11;
              case 1: return t12;
              case 2: return t13;
              default: return t14;
            }
    case 1: switch (j)
            {
              case 0: return t21;
              case 1: return t22;
              case 2: return t23;
              default: return t24;
            }
    case 2: switch (j)
            {
              case 0: return t31;
              case 1: return t32;
              case 2: return t33;
              default: return t34;
            }
    case 3: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(0);
              case 2: return FT(0);
              default: return FT(1);
            }
    }
  return FT(0);
  }

  virtual std::ostream &print(std::ostream &os) const
  {
    os <<"Aff_transformationC3("<<t11<<' '<<t12<<' '<<t13<<' '<<t14<<std::endl;
    os <<"                     "<<t21<<' '<<t22<<' '<<t23<<' '<<t24<<std::endl;
    os <<"                     "<<t31<<' '<<t32<<' '<<t33<<' '<<t34<<")";
    return os;
  }

private:
  FT   t11, t12, t13, t14; // FIXME : Wouldn't this be better with an array ?
  FT   t21, t22, t23, t24;
  FT   t31, t32, t33, t34;
};

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::inverse() const // FIXME : construction
{
  return Aff_transformation_3(
      determinant( t22, t23, t32, t33),         // i 11
     -determinant( t12, t13, t32, t33),         // i 12
      determinant( t12, t13, t22, t23),         // i 13
     -determinant( t12, t13, t14, t22, t23, t24, t32, t33, t34 ),

     -determinant( t21, t23, t31, t33),         // i 21
      determinant( t11, t13, t31, t33),         // i 22
     -determinant( t11, t13, t21, t23),         // i 23
      determinant( t11, t13, t14, t21, t23, t24, t31, t33, t34 ),

      determinant( t21, t22, t31, t32),         // i 31
     -determinant( t11, t12, t31, t32),         // i 32
      determinant( t11, t12, t21, t22),         // i 33
     -determinant( t11, t12, t14, t21, t22, t24, t31, t32, t34 ),

      determinant( t11, t12, t13, t21, t22, t23, t31, t32, t33 ));
}

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::
operator*(const Aff_transformation_rep_baseC3<R> &t) const
{
  return t.compose(*this);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::
compose(const Aff_transformation_repC3<R> &t) const // FIXME : construction
{
  return Aff_transformation_3(t.t11*t11 + t.t12*t21 + t.t13*t31,
                              t.t11*t12 + t.t12*t22 + t.t13*t32,
                              t.t11*t13 + t.t12*t23 + t.t13*t33,
                              t.t11*t14 + t.t12*t24 + t.t13*t34 + t.t14,

                              t.t21*t11 + t.t22*t21 + t.t23*t31,
                              t.t21*t12 + t.t22*t22 + t.t23*t32,
                              t.t21*t13 + t.t22*t23 + t.t23*t33,
                              t.t21*t14 + t.t22*t24 + t.t23*t34 + t.t24,

                              t.t31*t11 + t.t32*t21 + t.t33*t31,
                              t.t31*t12 + t.t32*t22 + t.t33*t32,
                              t.t31*t13 + t.t32*t23 + t.t33*t33,
                              t.t31*t14 + t.t32*t24 + t.t33*t34 + t.t34);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::
compose(const Translation_repC3<R> &t) const // FIXME : construction
{
  return Aff_transformation_3(t11,
                              t12,
                              t13,
                              t14 + t.translationvector_.x(),

                              t21,
                              t22,
                              t23,
                              t24 + t.translationvector_.y(),

                              t31,
                              t32,
                              t33,
                              t34 + t.translationvector_.z());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::
compose(const Scaling_repC3<R> &t) const // FIXME : construction
{
  return Aff_transformation_3(t.scalefactor_ * t11,
                              t.scalefactor_ * t12,
                              t.scalefactor_ * t13,
                              t.scalefactor_ * t14,

                              t.scalefactor_ * t21,
                              t.scalefactor_ * t22,
                              t.scalefactor_ * t23,
                              t.scalefactor_ * t24,

                              t.scalefactor_ * t31,
                              t.scalefactor_ * t32,
                              t.scalefactor_ * t33,
                              t.scalefactor_ * t34);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::transpose() const
{
  return Aff_transformation_3( t11, t21, t31, t14,
                               t12, t22, t32, t24,
                               t13, t23, t33, t34);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_H
