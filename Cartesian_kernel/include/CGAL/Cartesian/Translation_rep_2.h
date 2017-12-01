// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_TRANSLATION_REP_2_H
#define CGAL_CARTESIAN_TRANSLATION_REP_2_H

#include <CGAL/Cartesian/Aff_transformation_rep_2.h>

namespace CGAL {

template < class R >
class Translation_repC2 : public Aff_transformation_rep_baseC2<R>
{
friend class Aff_transformation_repC2<R>;
friend class Rotation_repC2<R>;
friend class Scaling_repC2<R>;

public:
  typedef typename R::FT                         FT;
  typedef Aff_transformation_rep_baseC2<R>       Aff_t_base;
  typedef Aff_transformation_repC2<R>            Transformation;
  typedef Translation_repC2<R>                   Translation;
  typedef Rotation_repC2<R>                      Rotation;
  typedef Scaling_repC2<R>                       Scaling;
  typedef typename Aff_t_base::Point_2           Point_2;
  typedef typename Aff_t_base::Vector_2          Vector_2;
  typedef typename Aff_t_base::Direction_2       Direction_2;
  typedef typename Aff_t_base::Aff_transformation_2 Aff_transformation_2;


  Translation_repC2() {}

  Translation_repC2(const Vector_2 &tv)
    : translationvector_(tv)
  {}

  Point_2     transform(const Point_2 &p) const
  { 
    typename R::Construct_translated_point_2 translated_point;
    return translated_point(p, translationvector_); 
  }

  Vector_2    transform(const Vector_2 &v) const { return v; }
  Direction_2 transform(const Direction_2 &d) const { return d; }

  Aff_transformation_2 operator*(const Aff_t_base &t) const
  {
    return t.compose(*this);
  }

  Aff_transformation_2 compose(const Translation &t) const
  {
    return Aff_transformation_2(TRANSLATION,
                                translationvector_ + t.translationvector_);
  }

  Aff_transformation_2 compose(const Rotation &t) const
  {
    return Aff_transformation_2(t.cosinus_,
                                -t.sinus_,
                                t.cosinus_*translationvector_.x() -
                                t.sinus_*translationvector_.y(),

                                t.sinus_,
                                t.cosinus_,
                                t.sinus_*translationvector_.x() +
                                t.cosinus_*translationvector_.y());
  }

  Aff_transformation_2 compose(const Scaling &t) const
  {
    return Aff_transformation_2(t.scalefactor_,
                                FT(0),
                                t.scalefactor_*translationvector_.x(),

                                FT(0),
                                t.scalefactor_,
                                t.scalefactor_*translationvector_.y());
  }

  Aff_transformation_2 compose(const Transformation &t) const
  {
    return Aff_transformation_2(t.t11,
                                t.t12,
                                t.t11 * translationvector_.x()
                                + t.t12 * translationvector_.y()
                                + t.t13,

                                t.t21,
                                t.t22,
                                t.t21 * translationvector_.x()
                                + t.t22*translationvector_.y()
                                + t.t23);
  }

  Aff_transformation_2 inverse() const
  {
    return Aff_transformation_2(TRANSLATION, - translationvector_);
  }

  bool         is_even() const
  {
    return true;
  }

  FT cartesian(int i, int j) const
  {
    if (j==i) return FT(1);
    if (j==2) return translationvector_[i];
    return FT(0);
  }

  std::ostream &print(std::ostream &os) const
  {
    os << "Aff_transformationC2(VectorC2(" << translationvector_.x() << ", "
       << translationvector_.y()  <<  "))";
    return os;
  }

private:
  Vector_2   translationvector_;
};

} //namespace CGAL

#endif // CGAL_CARTESIAN_TRANSLATION_REP_2_H
