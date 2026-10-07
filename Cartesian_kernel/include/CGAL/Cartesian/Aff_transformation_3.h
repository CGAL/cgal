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
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_3_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_3_H

#include <cmath>
#include <CGAL/aff_transformation_tags.h>

namespace CGAL {

template < class R_ >
class Aff_transformationC3
{
  friend class PlaneC3<R_>; // FIXME: why ?

  typedef typename R_::FT                   FT;

  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  std::vector<FT> entries;
  int variant; // 0: identity, 1: scaling, 2: translation, 3: general
public:
  typedef R_                               R;

  Aff_transformationC3()
  : variant(0)
  {}

  Aff_transformationC3(const Identity_transformation)
  : variant(0)
  {}

  Aff_transformationC3(const Scaling, const FT &s, const FT &w = FT(1))
  : entries(1), variant(1)
  {
    if (w != FT(1))
      entries[0] = s/w;
    else
      entries[0] = s;
  }

  Aff_transformationC3(const Translation, const Vector_3 &v)
  : entries(3), variant(2)
  {
    entries[0] = v.x();
    entries[1] = v.y();
    entries[2] = v.z();
  }

  // General form: without translation
  Aff_transformationC3(const FT& m11, const FT& m12, const FT& m13,
                       const FT& m21, const FT& m22, const FT& m23,
                       const FT& m31, const FT& m32, const FT& m33)
  : entries(12), variant(3)
  {
     FT zero(0);
     entries = {m11, m12, m13, zero,
                m21, m22, m23, zero,
                m31, m32, m33, zero };
  }

  // General form: without translation
  Aff_transformationC3(const FT& m11, const FT& m12, const FT& m13,
                       const FT& m21, const FT& m22, const FT& m23,
                       const FT& m31, const FT& m32, const FT& m33,
                       const FT& w)
  : entries(12), variant(3)
  {
    FT zero(0);
    entries = {m11/w, m12/w, m13/w, zero,
               m21/w, m22/w, m23/w, zero,
               m31/w, m32/w, m33/w, zero};
  }

  // General form: with translation
  Aff_transformationC3(
              const FT& m11, const FT& m12, const FT& m13, const FT& m14,
              const FT& m21, const FT& m22, const FT& m23, const FT& m24,
              const FT& m31, const FT& m32, const FT& m33, const FT& m34)
  : entries(12), variant(3)
  {
    entries = {m11, m12, m13, m14,
               m21, m22, m23, m24,
               m31, m32, m33, m34 };
  }

  // General form: with translation
  Aff_transformationC3(
              const FT& m11, const FT& m12, const FT& m13, const FT& m14,
              const FT& m21, const FT& m22, const FT& m23, const FT& m24,
              const FT& m31, const FT& m32, const FT& m33, const FT& m34,
              const FT& w)
  : entries(12), variant(3)
  {
    entries = {m11/w, m12/w, m13/w, m14/w,
               m21/w, m22/w, m23/w, m24/w,
               m31/w, m32/w, m33/w, m34/w };
  }


  template <class T>
  T transform(const T &t) const
  { return t.transform(*this); }

template <class T>
  T operator()(const T &t) const
  { return t.transform(*this); }

  Point_3
  transform(const Point_3 &p) const
  { if(variant == 0) return p;
    if(variant == 1) return Point_3(entries[0]*p.x(), entries[0]*p.y(), entries[0]*p.z());
    if(variant == 2) return Point_3(p.x()+entries[0], p.y()+entries[1], p.z()+entries[2]);
    typename R::Construct_point_3 construct_point_3;
    return construct_point_3(entries[0] * p.x() + entries[1] * p.y() + entries[2]  * p.z() + entries[3] ,
                             entries[4] * p.x() + entries[5] * p.y() + entries[6]  * p.z() + entries[7] ,
                             entries[8] * p.x() + entries[9] * p.y() + entries[10] * p.z() + entries[11] );
  }

  Point_3
  operator()(const Point_3 &p) const
  { return transform(p); }

  Vector_3
  transform(const Vector_3 &v) const
  { if(variant == 0) return v;
    if(variant == 1) return Vector_3(entries[0]*v.x(), entries[0]*v.y(), entries[0]*v.z());
    if(variant == 2) return v;
    typename R::Construct_vector_3 construct_vector_3;
    return construct_vector_3(entries[0] * v.x() + entries[1] * v.y() + entries[2]  * v.z(),
                              entries[4] * v.x() + entries[5] * v.y() + entries[6]  * v.z(),
                              entries[8] * v.x() + entries[9] * v.y() + entries[10] * v.z());
   }

   Vector_3
   operator()(const Vector_3 &v) const
  { return transform(v); }

  Direction_3
  transform(const Direction_3 &d) const
  { if(variant == 0) return d;
    if(variant == 1) return d;
    if(variant == 2) return d;

    Vector_3 v = d.to_vector();
    return transform(v).direction();
   }

  Direction_3
  operator()(const Direction_3 &d) const
  { return transform(d); }


  Plane_3
  transform(const Plane_3& p) const
  {
    if(variant == 0) return p;
    if(variant == 1) return Plane_3(p.a(),p.b(),p.c(), p.d()*entries[0]);
    if(variant == 2) return Plane_3(p.a(),
                                    p.b(),
                                    p.c(),
                                    p.d()  - ( p.a()*entries[0] + p.b()*entries[1]  + p.c()*entries[2] ));
    if (is_even())
      return Plane_3(transform(p.point()),
                     transpose().inverse().transform(p.orthogonal_direction()));
    else
      return Plane_3(transform(p.point()),
                     - transpose().inverse().transform(p.orthogonal_direction()));
  }

  Plane_3
  operator()(const Plane_3& p) const
  { return transform(p); } // FIXME : not compiled by the test-suite !

  Aff_transformation_3 inverse() const {
    if(variant == 0) return *this;
    if(variant == 1) return Aff_transformation_3(SCALING, FT(1)/entries[0]);
    if(variant == 2) return Aff_transformation_3(TRANSLATION, -Vector_3(entries[0],entries[1], entries[2]));
    const FT & t11 = entries[0], &t12 = entries[1], &t13 = entries[2], &t14 = entries[3];
    const FT & t21 = entries[4], &t22 = entries[5], &t23 = entries[6], &t24 = entries[7];
    const FT & t31 = entries[8], &t32 = entries[9], &t33 = entries[10], &t34 = entries[11];
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

  bool is_even() const {
    if(variant == 0) return true;
    if(variant == 1) return true;
    if(variant == 2) return true;
    return sign_of_determinant(entries[0], entries[1], entries[2],
                               entries[4], entries[5], entries[6],
                               entries[8], entries[9], entries[10]) == POSITIVE;
   }

  bool is_odd() const { return  ! is_even(); }
  bool is_translation() const { return variant == 2; }
  bool is_scaling() const { return variant == 1; }

  bool has_rotation() const {
    if(variant == 3)
      return  !(is_zero(entries[0]) && is_zero(entries[1]) && is_zero(entries[2]) && is_zero(entries[4]) && is_zero(entries[6]) && is_zero(entries[8]) && is_zero(entries[9]));
  }


  FT cartesian(int i, int j) const {
    if(variant == 0){
      if (i!=j) return FT(0);
      if (i==3) return FT(1);
      return FT(1);
    }
    if(variant == 1){
      if (i!=j) return FT(0);
      if (i==3) return FT(1);
      return entries[0];
    }
    if(variant == 2){
      if (j==i) return FT(1);
      if (j==3) return entries[i];
      return FT(0);
    }
    if(i < 3)
      return entries[i*4+j];
    if(i == 3 && j == 3)
      return FT(1);
    return FT(0);
  }

  FT homogeneous(int i, int j) const { return cartesian(i,j); }
  FT m(int i, int j) const { return cartesian(i,j); }
  FT hm(int i, int j) const { return cartesian(i,j); }

  Aff_transformation_3 compose(const Aff_transformation_3 &t) const
  {
    if(variant == 0){
      return t;
    }
    if(variant == 1){
      if(t.variant == 0) return *this;
      if(t.variant == 1) return Aff_transformation_3(SCALING, entries[0]*t.entries[0]);
      if(t.variant == 2) return Aff_transformation_3(TRANSLATION, Vector_3(entries[0]*t.entries[0], entries[0]*t.entries[1], entries[0]*t.entries[2]));
      return Aff_transformation_3(entries[0]*t.entries[0], entries[0]*t.entries[1], entries[0]*t.entries[2], entries[0]*t.entries[3],
                                  entries[0]*t.entries[4], entries[0]*t.entries[5], entries[0]*t.entries[6], entries[0]*t.entries[7],
                                  entries[0]*t.entries[8], entries[0]*t.entries[9], entries[0]*t.entries[10], entries[0]*t.entries[11]);
    }
    if(variant == 2){
      if(t.variant == 0) return *this;
      if(t.variant == 1) return Aff_transformation_3(TRANSLATION, Vector_3(entries[0]*t.entries[0], entries[1]*t.entries[0], entries[2]*t.entries[0]));
      if(t.variant == 2) return Aff_transformation_3(TRANSLATION, Vector_3(entries[0]+t.entries[0], entries[1]+t.entries[1], entries[2]+t.entries[2]));
      return Aff_transformation_3(t.entries[0]+entries[0], t.entries[1]+entries[1], t.entries[2]+entries[2], t.entries[3],
                                  t.entries[4]+entries[0], t.entries[5]+entries[1], t.entries[6]+entries[2], t.entries[7],
                                  t.entries[8]+entries[0], t.entries[9]+entries[1], t.entries[10]+entries[2], t.entries[11]);
    }
      if(t.variant == 0) return *this;
      if(t.variant == 1)
        return Aff_transformation_3(entries[0]*t.entries[0], entries[1]*t.entries[0], entries[2]*t.entries[0], entries[3]*t.entries[0],
                                    entries[4]*t.entries[0], entries[5]*t.entries[0], entries[6]*t.entries[0], entries[7]*t.entries[0],
                                    entries[8]*t.entries[0], entries[9]*t.entries[0], entries[10]*t.entries[0], entries[11]*t.entries[0]);
      if(t.variant == 2)
        return Aff_transformation_3(entries[0]+entries[3]*t.entries[0], entries[1]+entries[3]*t.entries[1], entries[2]+entries[3]*t.entries[2], entries[3],
                                    entries[4]+entries[7]*t.entries[0], entries[5]+entries[7]*t.entries[1], entries[6]+entries[7]*t.entries[2], entries[7],
                                    entries[8]+entries[11]*t.entries[0], entries[9]+entries[11]*t.entries[1], entries[10]+entries[11]*t.entries[2], entries[11]);

      return Aff_transformation_3(t.entries[0]*entries[0] + t.entries[1]*entries[4] + t.entries[2]*entries[8],
                                  t.entries[0]*entries[1] + t.entries[1]*entries[5] + t.entries[2]*entries[9],
                                  t.entries[0]*entries[2] + t.entries[1]*entries[6] + t.entries[2]*entries[10],
                                  t.entries[0]*entries[3] + t.entries[1]*entries[7] + t.entries[2]*entries[11] + t.entries[3],

                                  t.entries[4]*entries[0] + t.entries[5]*entries[4] + t.entries[6]*entries[8],
                                  t.entries[4]*entries[1] + t.entries[5]*entries[5] + t.entries[6]*entries[9],
                                  t.entries[4]*entries[2] + t.entries[5]*entries[6] + t.entries[6]*entries[10],
                                  t.entries[4]*entries[3] + t.entries[5]*entries[7] + t.entries[6]*entries[11] + t.entries[7],

                                  t.entries[8]*entries[0] + t.entries[9]*entries[4] + t.entries[10]*entries[8],
                                  t.entries[8]*entries[1] + t.entries[9]*entries[5] + t.entries[10]*entries[9],
                                  t.entries[8]*entries[2] + t.entries[9]*entries[6] + t.entries[10]*entries[10],
                                  t.entries[8]*entries[3] + t.entries[9]*entries[7] + t.entries[10]*entries[11] + t.entries[11]);

  }

  Aff_transformation_3 operator*(const Aff_transformation_3 &t) const
  { return t.compose(*this); }

  std::ostream &
  print(std::ostream &os) const;

  bool operator==(const Aff_transformation_3 &t)const
  {
    for(int i=0; i<3; ++i)
      for(int j = 0; j< 4; ++j)
        if(cartesian(i,j)!=t.cartesian(i,j))
          return false;
    return true;
  }

  bool operator!=(const Aff_transformation_3 &t)const
  {
    return !(*this == t);
  }

protected:
  Aff_transformation_3  transpose() const {
    if(variant == 0) return *this;
    if(variant == 1) return *this;
    if(variant == 2) return *this;
    return Aff_transformation_3( entries[0], entries[4], entries[8], entries[3],
                                 entries[1], entries[5], entries[9], entries[7],
                                 entries[2], entries[6], entries[10], entries[11]);
  }
};


template < class R >
std::ostream&
Aff_transformationC3<R>::print(std::ostream &os) const
{
  if(variant == 0) {
    return os << "Identity";
  }
  if(variant == 1) {
    os << "Scaling: " << entries[0];
  }
  if(variant == 2) {
    return os << "Translation: " << entries[0] << " " << entries[1] << " " << entries[2];
  }
  return os << entries[0] << " " << entries[1] << " " << entries[2] << " " << entries[3] << "\n"
            << entries[4] << " " << entries[5] << " " << entries[6] << " " << entries[7] << "\n"
            << entries[8] << " " << entries[9] << " " << entries[10] << " " << entries[11];
  return os;
}

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC3
template < class R >
std::ostream&
operator<<(std::ostream &os, const Aff_transformationC3<R> &t)
{
  t.print(os);
  return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC3

} //namespace CGAL

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_3_H
