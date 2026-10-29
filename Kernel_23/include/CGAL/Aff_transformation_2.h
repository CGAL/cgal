// Copyright (c) 1999
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
// Author(s)     : Andreas Fabri, Stefan Schirra

#ifndef CGAL_AFF_TRANSFORMATION_2_H
#define CGAL_AFF_TRANSFORMATION_2_H

#include <CGAL/config.h>
#include <CGAL/Dimension.h>
#include <CGAL/aff_transformation_tags.h>

#include <iostream>

namespace CGAL {

template <class R_>
class Aff_transformation_2 : public R_::Kernel_base::Aff_transformation_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::FT                    FT;
  typedef typename R_::Line_2                Line_2;
  typedef typename R_::Direction_2           Direction_2;
  typedef typename R_::Vector_2              Vector_2;


  typedef Aff_transformation_2       Self;

public:

  typedef CGAL::Dimension_tag<2>            Ambient_dimension;


  typedef typename R_::Kernel_base::Aff_transformation_2 Rep;


  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }


  typedef  R_                               R;

  Aff_transformation_2() {}

  Aff_transformation_2(const Rep&  t)
    : Rep(t) {}

  Aff_transformation_2(Rep&& t)
    : Rep(std::move(t)) {}



  Aff_transformation_2(const Identity_transformation tag)
    : Rep(tag)
  {}

  Aff_transformation_2(const Translation tag, const Vector_2 &v)
    : Rep(tag, v)
  {}

  // Rational Rotation:
  Aff_transformation_2(const Rotation tag,
                       const Direction_2 &d,
                       const RT &num,
                       const RT &den = RT(1))
    : Rep(tag, d, num, den)
  {}

  Aff_transformation_2(const Rotation tag,
                       const RT &sin,
                       const RT &cos,
                       const RT &den = RT(1))
    : Rep(tag, sin, cos, den)
  {}

  Aff_transformation_2(const Reflection tag, const Line_2& l )
    : Rep(tag, l)
  {}

  Aff_transformation_2(const Scaling tag,
                       const RT &s,
                       const RT &w= RT(1))
    : Rep(tag, s, w)
  {}

  // The general case:
  Aff_transformation_2(const RT & m11,
                       const RT & m12,
                       const RT & m13,

                       const RT & m21,
                       const RT & m22,
                       const RT & m23,

                       const RT &w= RT(1))
    : Rep(m11, m12, m13,
          m21, m22, m23,
                         w)
  {}

  Aff_transformation_2(const RT & m11, const RT & m12,
                       const RT & m21, const RT & m22,
                       const RT &w = RT(1))
    : Rep(m11, m12,
          m21, m22,
                   w)
  {}
};


template < class R >
std::ostream &
operator<<(std::ostream &os, const CGAL::Aff_transformation_2<R> &t)
{
  typedef typename R::Kernel_base::Aff_transformation_2  RAff_transformation_2;
  return os << static_cast<const RAff_transformation_2&>(t);
}

template < class R >
std::istream &
operator>>(std::istream &is, CGAL::Aff_transformation_2<R> &t)
{
  typedef typename R::Kernel_base::Aff_transformation_2  RAff_transformation_2;
  return is >> static_cast<RAff_transformation_2&>(t);
}


} //namespace CGAL

#endif // CGAL_AFF_TRANSFORMATION_2_H
