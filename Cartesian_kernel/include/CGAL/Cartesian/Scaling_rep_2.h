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

#ifndef CGAL_CARTESIAN_SCALING_REP_2_H
#define CGAL_CARTESIAN_SCALING_REP_2_H

namespace CGAL {

template < class R >
class Scaling_repC2: public Aff_transformation_rep_baseC2<R>
{
friend class Aff_transformation_repC2<R>;
friend class Translation_repC2<R>;
friend class Rotation_repC2<R>;
friend class Reflection_repC2<R>;

public:
  typedef Aff_transformation_rep_baseC2<R> Aff_t_base;
  typedef typename Aff_t_base::FT                FT;
  typedef typename Aff_t_base::Point_2           Point_2;
  typedef typename Aff_t_base::Vector_2          Vector_2;
  typedef typename Aff_t_base::Direction_2       Direction_2;
  typedef typename Aff_t_base::Aff_transformation_2   Aff_transformation_2;
  typedef Aff_transformation_repC2<R>      Transformation;
  typedef Translation_repC2<R>             Translation;
  typedef Rotation_repC2<R>                Rotation;
  typedef Reflection_repC2<R>              Reflection;
  typedef Scaling_repC2<R>                 Scaling;

  Scaling_repC2()
  {}

  Scaling_repC2(const FT &scalefactor) :
    scalefactor_(scalefactor)
  {}

  ~Scaling_repC2()
  {}

  Point_2      transform(const Point_2 &p) const
  {
    return Point_2(scalefactor_ * p.x(), scalefactor_ * p.y());
  }

  Vector_2      transform(const Vector_2 &p) const
  {
    return Vector_2(scalefactor_ * p.x(), scalefactor_ * p.y());
  }

  Direction_2  transform(const Direction_2 &d) const
  {
    return d;
  }

  Aff_transformation_2 operator*(const Aff_t_base &t) const
  {
   return t.compose(*this);
  }

  Aff_transformation_2 compose(const Translation &t) const
  {
    FT ft0(0);
    return Aff_transformation_2(scalefactor_,
                                ft0,
                                t.translationvector_.x(),
                                ft0,
                                scalefactor_,
                                t.translationvector_.y());
  }

  Aff_transformation_2 compose(const Rotation &t) const
  {
    return Aff_transformation_2(scalefactor_ * t.cosinus_,
                                scalefactor_ * -t.sinus_,
                                scalefactor_ * t.sinus_,
                                scalefactor_ * t.cosinus_);
  }

  Aff_transformation_2 compose(const Reflection &r) const
  {
    return Aff_transformation_2(scalefactor_*r.cosinus_, scalefactor_*r.sinus_,-r.cosinus_*r.t.x()-r.sinus_*r.t.y()+r.t.x(),
                                scalefactor_*r.sinus_, -scalefactor_*r.cosinus_, -r.sinus_*r.t.x()+r.cosinus_*r.t.y()-r.t.y());
  }

  Aff_transformation_2 compose(const Scaling &t) const
  {
    return Aff_transformation_2(SCALING, scalefactor_*t.scalefactor_);
  }

  Aff_transformation_2 compose(const Transformation &t) const
  {
    return Aff_transformation_2(scalefactor_ * t.t11,
                                scalefactor_ * t.t12,
                                t.t13,
                                scalefactor_ * t.t21,
                                scalefactor_ * t.t22,
                                t.t23);
  }

  Aff_transformation_2  inverse() const
  {
    return Aff_transformation_2(SCALING, FT(1)/scalefactor_);
  }

  bool is_even() const
  {
    return true;
  }

  FT cartesian(int i, int j) const
  {
    if (i!=j) return FT(0);
    return (i==2) ? FT(1) : scalefactor_;
  }

  std::ostream &print(std::ostream &os) const
  {
    os << "Aff_transformationC2(" << scalefactor_ <<  ")";
    return os;
  }

private:
  FT scalefactor_;
};

} //namespace CGAL

#endif // CGAL_CARTESIAN_SCALING_REP_2_H
