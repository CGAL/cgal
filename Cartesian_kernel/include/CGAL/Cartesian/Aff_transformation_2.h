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
// Author(s)     : Andreas Fabri, Lutz Kettner

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H

#include <CGAL/aff_transformation_tags.h>
#include <CGAL/determinant.h>
#include <CGAL/rational_rotation.h>

#include <ostream>
#include <vector>
#include <cmath>

namespace CGAL {

class Identity_transformation;


template < class R_ >
class Aff_transformationC2
{
  typedef typename R_::FT                   FT;

  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  enum Variant {identity = 0, scaling, translation, rotation, reflection, general};
  std::vector<FT> entries;
  Variant variant;

public:
  typedef R_                                R;

  Aff_transformationC2()
  : variant(identity)
  {}

  Aff_transformationC2(const Identity_transformation)
  : variant(identity)
  {}

  Aff_transformationC2(const Scaling, const FT &s, const FT &w = FT(1))
  : entries(1), variant(scaling)
  {
    if (w != FT(1))
      entries[0] = s/w;
    else
      entries[0] = s;
  }


  Aff_transformationC2(const Translation, const Vector_2 &v)
  : entries(2), variant(translation)
  {
    entries[0] = v.x();
    entries[1] = v.y();
  }


  // Rational Rotation:
  Aff_transformationC2(const Rotation,
                       const Direction_2 &d,
                       const FT &num,
                       const FT &den = FT(1))
  : entries(2), variant(rotation)
  {
    FT sine;
    FT cosine;
    FT w;

    rational_rotation_approximation(d.dx(),
                                    d.dy(),
                                    sine,
                                    cosine,
                                    w,
                                    num,
                                    den);
          entries[0] = sine/w;
          entries[1] = cosine/w;
  }

  Aff_transformationC2(const Rotation,
                       const FT &sine,
                       const FT &cosine,
                       const FT &w = FT(1))
  : entries(2), variant(rotation)
  {
          entries[0] = sine/w;
          entries[1] = cosine/w;
  }

  const FT& sinus() const { return entries[0];}
  const FT& cosinus() const { return entries[1];}

  Aff_transformationC2(const Reflection, const Line_2& l)
  : entries(4), variant(reflection)
  {
    Vector_2 t;
    if(l.a() == 0)
      t = -Vector_2(0, l.c()/l.b());
    else
      t = -Vector_2(l.c()/l.a(),0);
    entries[2]= t.x();
    entries[3] =  t.y();
    Vector_2 l_to_v = l.to_vector();
    FT scal = l_to_v.x(); //Projection of l_to_v on Ox. = |L|*cos(a)
    FT det = l_to_v.y();// = |L|*sin(a)
    entries[0] = 2*det*scal/l_to_v.squared_length(); //sin(2a) = 2*sin(a)*cos(a)
    FT sq_cos = scal*scal/l_to_v.squared_length(); //cos(a)*cos(a)
    entries[1] = 2*sq_cos-1;
  }
  const FT& tx()const { return entries[2];}
  const FT& ty()const { return entries[3];}

 // General form: without translation
  Aff_transformationC2(const FT& m11, const FT& m12,
                       const FT& m21, const FT& m22)
  : entries(6), variant(general)
  {
     FT zero(0);
     entries = {m11, m12, zero,
                m21, m22, zero};
  }

  // General form: without translation
  Aff_transformationC2(const FT& m11, const FT& m12,
                       const FT& m21, const FT& m22,
                       const FT& w)
  : entries(6), variant(general)
  {
    FT zero(0);
    entries = {m11/w, m12/w, zero,
               m21/w, m22/w, zero};
  }

  // General form: with translation
  Aff_transformationC2(
              const FT& m11, const FT& m12, const FT& m13,
              const FT& m21, const FT& m22, const FT& m23)
  : entries(6), variant(general)
  {
    entries = {m11, m12, m13,
               m21, m22, m23 };
  }

  // General form: with translation
  Aff_transformationC2(
              const FT& m11, const FT& m12, const FT& m13,
              const FT& m21, const FT& m22, const FT& m23,
              const FT& w)
  : entries(6), variant(general)
  {
    entries = {m11/w, m12/w, m13/w,
               m21/w, m22/w, m23/w};
  }


  Point_2
  transform(const Point_2&p) const
  { if(is_identity()) return p;
    typename R::Construct_point_2 construct_point_2;
    if(is_scaling()) return construct_point_2(entries[0]*p.x(), entries[0]*p.y());
    if(is_translation()) return construct_point_2(p.x()+entries[0], p.y()+entries[1]);
    if(is_rotation())
      return construct_point_2(cosinus() * p.x() - sinus() * p.y(),
                               sinus() * p.x() + cosinus() * p.y());
    if(is_reflection())
      return construct_point_2(
          cosinus()*p.x()+sinus()*p.y()-cosinus()*tx()-sinus()*ty()+tx(),
          sinus()*p.x()-cosinus()*p.y()-sinus()*tx()+cosinus()*ty()+ty());

    return construct_point_2(entries[0] * p.x() + entries[1] * p.y() + entries[2] ,
                             entries[3] * p.x() + entries[4] * p.y() + entries[5] );
  }

  Point_2
  operator()(const Point_2 &p) const
  { return transform(p); }

  Vector_2
  transform(const Vector_2 &v) const
  { if(is_identity()) return v;
    typename R::Construct_vector_2 construct_vector_2;

    if(is_scaling()) return construct_vector_2(entries[0]*v.x(), entries[0]*v.y());
    if(is_translation()) return v;
    if(is_rotation())
      return construct_vector_2(cosinus() * v.x() - sinus() * v.y(),
                                sinus() * v.x() + cosinus() * v.y());
    if(is_reflection())
      return construct_vector_2(cosinus()*v.x()+sinus()*v.y()-cosinus()*tx()-sinus()*ty()+tx(),
                                sinus()*v.x()-cosinus()*v.y()-sinus()*tx()+cosinus()*ty()+ty());
    return construct_vector_2(entries[0] * v.x() + entries[1] * v.y() ,
                              entries[3] * v.x() + entries[4] * v.y());
   }

  Vector_2
  operator()(const Vector_2 &v) const
  { return transform(v); } // FIXME : not compiled by the test-suite.


  Direction_2
  transform(const Direction_2 &d) const
  { if(is_identity()) return d;
    if(is_scaling()) return d;
    if(is_translation()) return d;
    if(is_rotation())
      return Direction_2(cosinus() * d.dx() - sinus() * d.dy(),
                         sinus() * d.dx() + cosinus() * d.dy());
    if(is_reflection())
      return transform(d.vector()).direction();

    Vector_2 v = d.to_vector();
    return transform(v).direction();
   }


  Direction_2
  operator()(const Direction_2 &d) const
  { return transform(d); } // FIXME : not compiled by the test-suite.


  Line_2
  transform(const Line_2 &l) const
  { assert(false);return l; }

  Line_2
  operator()(const Line_2 &l) const
  { return transform(l); }

  Aff_transformation_2 inverse() const
  {
   if(is_identity()) return *this;

    if(is_scaling()) return Aff_transformation_2(SCALING,FT(1)/entries[0]);
    if(is_translation())
      return Aff_transformation_2(TRANSLATION, Vector_2(-entries[0], -entries[1]));
    if(is_reflection())
      return Aff_transformation_2(cartesian(0,0), cartesian(0,1), cartesian(0,2),
                                  cartesian(1,0), cartesian(1,1), cartesian(1,2));;
    if(is_rotation())
      return Aff_transformation_2(ROTATION, - sinus(), cosinus(), FT(1));

   FT det = FT(1) / (entries[0] * entries[4] - entries[1] * entries[3]);
  return Aff_transformation_2(
    det * entries[4],    det * (-entries[1]), det * (entries[1]*entries[5]-entries[2]*entries[4]),
    det * (-entries[3]), det * entries[0] ,   det * (entries[2]*entries[3]-entries[0]*entries[5]));
   }

  bool is_even() const
  {
    if(is_identity()) return true;
    if(is_scaling()) return true;
    if(is_translation()) return true;
    if(is_reflection()) return true;
    if(is_rotation()) { return true;}
    return sign_of_determinant(entries[0], entries[1],
                               entries[3], entries[4]) == POSITIVE;
  }

  bool is_odd() const { return ! is_even(); }

  bool is_identity() const { return variant == identity; }
  bool is_translation() const { return variant == translation; }
  bool is_scaling() const { return variant == scaling; }
  bool is_reflection() const { return variant == reflection; }
  bool is_rotation() const { return variant == rotation; }
  bool is_general() const { return variant == general; }

  bool has_rotation() const {
    if(is_rotation() || is_reflection()) return true;
    return false;
  }

  FT cartesian(int i, int j) const
  {
    if(is_identity()){
      if (i!=j) return FT(0);
      if (i==2) return FT(1);
      return FT(1);
    }
    if(is_scaling()){
      if (i!=j) return FT(0);
      if (i==2) return FT(1);
      return entries[0];
    }
    if(is_translation()){
      if (j==i) return FT(1);
      if (j==2) return entries[i];
      return FT(0);
    }
    if(is_rotation()){
      switch (i)
      {
      case 0: switch (j)
            {
              case 0: return cosinus();
              case 1: return -sinus();
              default: return FT(0);
            }
    case 1: switch (j)
            {
              case 0: return sinus();
              case 1: return cosinus();
              default: return FT(0);
            }
    case 2: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(0);
              default: return FT(1);
            }
    }
    return FT(0);
    }


    if(is_reflection()){
      switch (i)
      {
      case 0: switch (j)
            {
              case 0: return cosinus();
              case 1: return sinus();
              default: return -cosinus()*entries[2]-sinus()*entries[3]+entries[2];
            }
    case 1: switch (j)
            {
              case 0: return sinus();
              case 1: return -cosinus();
              default: return -sinus()*entries[2]+cosinus()*entries[3]+entries[3];
            }
    case 2: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(0);
              default: return FT(1);
            }
    }
    return FT(0);
    }


    if(i < 2)
      return entries[i*3+j];
    if(i == 2 && j == 2)
      return FT(1);
    return FT(0);
  }

  FT homogeneous(int i, int j) const { return cartesian(i,j); }
  FT m(int i, int j) const { return cartesian(i,j); }
  FT hm(int i, int j) const { return cartesian(i,j); }

  Aff_transformation_2 operator*(const Aff_transformationC2 &t) const
  {
    if(is_identity()) return t;
    if(t.is_identity()) return *this;
    if(is_scaling() && t.is_scaling())
      return Aff_transformation_2(SCALING, entries[0]* t.entries[0]);
    if(is_translation() && t.is_translation())
      return Aff_transformation_2(TRANSLATION, Vector_2(entries[0] + t.entries[0],
                                                        entries[1] + t.entries[1]));
    return Aff_transformation_2(m(0,0)*t.m(0,0) + m(0,1)*t.m(1,0),
                              m(0,0)*t.m(0,1) + m(0,1)*t.m(1,1),
                              m(0,0)*t.m(0,2) + m(0,1)*t.m(1,2) + m(0,2),
                              m(1,0)*t.m(0,0) + m(1,1)*t.m(1,0),
                              m(1,0)*t.m(0,1) + m(1,1)*t.m(1,1),
                              m(1,0)*t.m(0,2) + m(1,1)*t.m(1,2) + m(1,2) );
  }

  std::ostream &
  print(std::ostream &os) const;

  bool operator==(const Aff_transformationC2 &t)const
  {
    for(int i=0; i<2; ++i)
      for(int j = 0; j< 3; ++j)
        if(cartesian(i,j)!=t.cartesian(i,j))
          return false;
    return true;
  }

  bool operator!=(const Aff_transformationC2 &t)const
  {
    return !(*this == t);
  }

};

template < class R >
std::ostream&
Aff_transformationC2<R>::print(std::ostream &os) const
{
  if(is_identity()) {
    return os << "Identity";
  }
  if(is_scaling()) {
    return os << "Scaling: " << entries[0];
  }
  if(is_translation()) {
    return os << "Translation: " << entries[0] << " " << entries[1] ;
  }
  if(is_rotation()) {
    return os << "Rotation: " << entries[0] << " " << entries[1] ;
  }

  if(is_reflection()) {
    return os << "Reflection: " << entries[0] << " " << entries[1] << " " << entries[2] << " " << entries[3] ;
  }
  return os << entries[0] << " " << entries[1] << " " << entries[2]  << "\n"
            << entries[3] << " " << entries[4] << " " << entries[5] ;
}


template < class R >
std::ostream&
operator<<(std::ostream& os, const Aff_transformationC2<R>& t)
{
  t.print(os);
  return os;
}


} //namespace CGAL

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H
