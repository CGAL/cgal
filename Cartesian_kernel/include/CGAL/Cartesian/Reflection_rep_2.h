// Copyright (c) 2018
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
// Author(s)     : Maxime Gimeno

#ifndef CGAL_CARTESIAN_REFLECTION_REP_2_H
#define CGAL_CARTESIAN_REFLECTION_REP_2_H
#include <cmath>
namespace CGAL {

template < class R >
class Reflection_repC2: public Aff_transformation_rep_baseC2<R>
{
  friend class Translation_repC2<R>;
  friend class Rotation_repC2<R>;
  friend class Scaling_repC2<R>;
  friend class Aff_transformation_repC2<R>;

public:
  typedef Aff_transformation_rep_baseC2<R> Aff_t_base;
  typedef typename Aff_t_base::FT                FT;
  typedef typename Aff_t_base::Point_2           Point_2;
  typedef typename Aff_t_base::Vector_2          Vector_2;
  typedef typename Aff_t_base::Direction_2       Direction_2;
typedef typename CGAL::Line_2<R>                 Line_2;
  typedef typename Aff_t_base::Aff_transformation_2   Aff_transformation_2;
  typedef Aff_transformation_repC2<R>      Transformation;
  typedef Reflection_repC2<R>              Reflection;
  typedef Scaling_repC2<R>                 Scaling;
  typedef Rotation_repC2<R>                Rotation;
  typedef Translation_repC2<R>             Translation;

  Reflection_repC2(const Line_2 &l)
  {
    if(l.a() == 0)
      t = -Vector_2(0, l.c()/l.b());
    else
      t = -Vector_2(l.c()/l.a(),0);

    Vector_2 l_to_v = l.to_vector();
    FT scal = l_to_v.x(); //Projection of l_to_v on Ox. = |L|*cos(a)
    FT det = l_to_v.y();// = |L|*sin(a)
    sinus_ = 2*det*scal/l_to_v.squared_length(); //sin(2a) = 2*sin(a)*cos(a)
    FT sq_cos = scal*scal/l_to_v.squared_length(); //cos(a)*cos(a)
    cosinus_ = 2*sq_cos-1;
  }

  ~Reflection_repC2()
  {}

  Point_2      transform(const Point_2 &p) const
  {
    return Point_2(
          cosinus_*p.x()+sinus_*p.y()-cosinus_*t.x()-sinus_*t.y()+t.x(),
          sinus_*p.x()-cosinus_*p.y()-sinus_*t.x()+cosinus_*t.y()+t.y());
  }

  Vector_2      transform(const Vector_2 &p) const
  {
    return Vector_2(
          cosinus_*p.x()+sinus_*p.y()-cosinus_*t.x()-sinus_*t.y()+t.x(),
          sinus_*p.x()-cosinus_*p.y()-sinus_*t.x()+cosinus_*t.y()+t.y());
  }

  Direction_2  transform(const Direction_2 &d) const
  {

    return transform(d.vector()).direction();
  }

  Aff_transformation_2 operator*(const Aff_t_base &t) const
  {
    return t.compose(*this);
  }

  Aff_transformation_2 compose(const Translation &tr) const
  {
    return Aff_transformation_2(cosinus_, sinus_, t13()+tr.translationvector_.x(),
                                sinus_, -cosinus_, t23()+tr.translationvector_.y());
  }

  Aff_transformation_2 compose(const Scaling &s) const
  {
    return Aff_transformation_2(s.scalefactor_ * cosinus_,
                                s.scalefactor_ * sinus_,
                                s.scalefactor_ * t13(),
                                s.scalefactor_ * sinus_,
                                -s.scalefactor_ * cosinus_,
                                s.scalefactor_ * t23());
  }

  Aff_transformation_2 compose(const Transformation &tr) const
  {
    return Aff_transformation_2(
          tr.t11*cosinus_+tr.t12*sinus_,
          tr.t11*sinus_-tr.t12*cosinus_,
          tr.t11*t13()+tr.t12*t23()+tr.t13,
          tr.t21*cosinus_+tr.t22*sinus_,
          tr.t21*sinus_-tr.t22*cosinus_,
          tr.t21*t13()+tr.t22*t23()+tr.t23);
  }

  Aff_transformation_2 compose(const Rotation &r) const
  {
    return Aff_transformation_2(
          r.cosinus_*cosinus_-r.sinus_*sinus_,
          r.cosinus_*sinus_+r.sinus_*cosinus_,
          r.cosinus_*t13()
          -r.sinus_*t23(),
          r.sinus_*cosinus_+r.cosinus_*sinus_,
          r.sinus_*sinus_-r.cosinus_*cosinus_,
          r.sinus_*t13()
          +r.cosinus_*t23());
  }

  Aff_transformation_2 compose(const Reflection &r) const
  {
    return Aff_transformation_2(
          cosinus_*r.cosinus_+sinus_*r.sinus_,
          r.cosinus_*sinus_-r.sinus_*cosinus_,
          r.cosinus_*(t13()-r.t.x()) + r.sinus_*(t23()-r.t.y())+r.t.x(),

          r.sinus_*cosinus_ - r.cosinus_*sinus_,
          r.sinus_*sinus_+r.cosinus_*cosinus_,
          r.sinus_*(t13()-r.t.x()) -r.cosinus_*(t23()-r.t.y())+r.t.y());
  }

  Aff_transformation_2  inverse() const
  {
    return Aff_transformation_2(cartesian(0,0), cartesian(0,1), cartesian(0,2),
                                cartesian(1,0), cartesian(1,1), cartesian(1,2));
  }

  bool is_even() const
  {
    return true;
  }

  FT cartesian(int i, int j) const
  {
    switch (i)
    {
    case 0: switch (j)
            {
              case 0: return cosinus_;
              case 1: return sinus_;
              default: return FT(0);
            }
    case 1: switch (j)
            {
              case 0: return sinus_;
              case 1: return -cosinus_;
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

  std::ostream &print(std::ostream &os) const
  {
    os << "Aff_transformationC2(" << sinus_ << ", " << cosinus_ <<  "; "<< t <<")";
    return os;
  }

  //convevience functions for composition
  FT t13()const
  {
    return FT(-cosinus_*t.x()-sinus_*t.y()+t.x());
  }
  FT t23()const
  {
    return FT(-sinus_*t.x()+cosinus_*t.y()+t.y());
  }

private:
  Vector_2 t;
  FT sinus_, cosinus_;
};

} //namespace CGAL

#endif // CGAL_CARTESIAN_REFLECTION_REP_2_H
