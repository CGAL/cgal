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

#ifndef CGAL_CARTESIAN_SCALING_REP_3_H
#define CGAL_CARTESIAN_SCALING_REP_3_H

namespace CGAL {

template < class R >
class Scaling_repC3 : public Aff_transformation_rep_baseC3<R>
{
  friend class Aff_transformation_repC3<R>;
  friend class Translation_repC3<R>;

public:
  typedef typename R::FT                                FT;
  typedef Aff_transformation_rep_baseC3<R>              Transformation_base_3;
  typedef Aff_transformation_repC3<R>                   Transformation_3;
  typedef Translation_repC3<R>                          Translation_3;
  typedef Scaling_repC3<R>                              Scaling_3;
  typedef typename Transformation_base_3::Point_3       Point_3;
  typedef typename Transformation_base_3::Vector_3      Vector_3;
  typedef typename Transformation_base_3::Direction_3   Direction_3;
  typedef typename Transformation_base_3::Plane_3       Plane_3;
  typedef typename Transformation_base_3::Aff_transformation_3
                                                        Aff_transformation_3;

  Scaling_repC3() {}
  Scaling_repC3(const FT &s) : scalefactor_(s) {}
  ~Scaling_repC3() override {}

  Point_3      transform(const Point_3 &p) const override
  {
    return Point_3(scalefactor_ * p.x(),
                   scalefactor_ * p.y(),
                   scalefactor_ * p.z());
  }

  Vector_3     transform(const Vector_3 &v) const override
  {
    return Vector_3(scalefactor_ * v.x(), scalefactor_ * v.y(),
                    scalefactor_ * v.z());
  }

  Direction_3  transform(const Direction_3 &d) const override
  {
    return d;
  }

  Plane_3  transform(const Plane_3 &p) const override
  {
    // direction ( which is (p.a(), p.b(), p.c())) does not change
    return Plane_3(p.a(),p.b(),p.c(), p.d()*scalefactor_);
  }

  Aff_transformation_3 operator*(const Transformation_base_3 &t) const override
  {
    return t.compose(*this);
  }

  Aff_transformation_3 compose(const Transformation_3 &t) const override
  {
    return Aff_transformation_3(scalefactor_ * t.t11,
                                scalefactor_ * t.t12,
                                scalefactor_ * t.t13,
                                t.t14,

                                scalefactor_ * t.t21,
                                scalefactor_ * t.t22,
                                scalefactor_ * t.t23,
                                t.t24,

                                scalefactor_ * t.t31,
                                scalefactor_ * t.t32,
                                scalefactor_ * t.t33,
                                t.t34);
  }

  Aff_transformation_3 compose(const Translation_3 &t) const override
  {
    FT ft0(0);
    return Aff_transformation_3(scalefactor_,
                                ft0,
                                ft0,
                                t.translationvector_.x(),

                                ft0,
                                scalefactor_,
                                ft0,
                                t.translationvector_.y(),

                                ft0,
                                ft0,
                                scalefactor_,
                                t.translationvector_.z());
  }

  Aff_transformation_3 compose(const Scaling_3 &t) const override
  {
    return Aff_transformation_3(SCALING, scalefactor_*t.scalefactor_);
  }

  Aff_transformation_3 inverse() const override
  {
    return Aff_transformation_3(SCALING, FT(1)/scalefactor_);
  }

  Aff_transformation_3 transpose() const override
  {
    return Aff_transformation_3(SCALING, scalefactor_);
  }

  bool is_even() const override
  {
    return true;
  }

  bool is_scaling() const override
  {
    return true;
  }

  FT cartesian(int i, int j) const override
  {
    if (i!=j) return FT(0);
    if (i==3) return FT(1);
    return scalefactor_;
  }

  std::ostream &print(std::ostream &os) const override
  {
    os << "Aff_transformationC3(" << scalefactor_ << ")";
    return os;
  }

private:
  FT   scalefactor_;
};

} //namespace CGAL

#endif // CGAL_CARTESIAN_SCALING_REP_3_H
