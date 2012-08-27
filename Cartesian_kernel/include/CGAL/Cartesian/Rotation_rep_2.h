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
// 
//
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_ROTATION_REP_2_H
#define CGAL_CARTESIAN_ROTATION_REP_2_H

#include <CGAL/rational_rotation.h>

namespace CGAL {

template < class R >
class Rotation_repC2: public Aff_transformation_rep_baseC2<R>
{
friend class Aff_transformation_repC2<R>;
friend class Translation_repC2<R>;
friend class Scaling_repC2<R>;

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
  typedef Scaling_repC2<R>                 Scaling;

  Rotation_repC2() {}

  Rotation_repC2(const FT &sinus, const FT &cosinus)
    : sinus_(sinus), cosinus_(cosinus) {}

  Rotation_repC2(const Direction_2 &d,
                 const FT &eps_num,
                 const FT &eps_den = FT(1))
  {
    FT sin_num;
    FT cos_num;
    FT denom;

    rational_rotation_approximation(d.dx(),
                                    d.dy(),
                                    sin_num,
                                    cos_num,
                                    denom,
                                    eps_num,
                                    eps_den);
    sinus_   = sin_num/denom;
    cosinus_ = cos_num/denom;
  }

  Point_2      transform(const Point_2 &p) const
  {
    return Point_2(cosinus_ * p.x() - sinus_ * p.y(),
                   sinus_ * p.x() + cosinus_ * p.y());
  }

  Vector_2     transform(const Vector_2 &v) const
  {
    return Vector_2(cosinus_ * v.x() - sinus_ * v.y(),
                    sinus_ * v.x() + cosinus_ * v.y());
  }

  Direction_2  transform(const Direction_2 &d) const
  {
    return Direction_2(cosinus_ * d.dx() - sinus_ * d.dy(),
                       sinus_ * d.dx() + cosinus_ * d.dy());
  }

  Aff_transformation_2 inverse() const
  {
    return Aff_transformation_2(ROTATION, - sinus_, cosinus_, FT(1));
  }

  Aff_transformation_2 operator*(const Aff_t_base &t) const
  {
    return t.compose(*this);
  }

  Aff_transformation_2 compose(const Translation &t) const
  {
    return Aff_transformation_2(cosinus_,
                                -sinus_,
                                t.translationvector_.x(),
                                sinus_,
                                cosinus_,
                                t.translationvector_.y());
  }

  Aff_transformation_2 compose(const Rotation &t) const
  {
    return Aff_transformation_2(ROTATION,
                                t.sinus_*cosinus_ + t.cosinus_*sinus_,
                                t.cosinus_*cosinus_-t.sinus_*sinus_ );
  }

  Aff_transformation_2 compose(const Scaling &t) const
  {
    return Aff_transformation_2(t.scalefactor_*cosinus_,
                                t.scalefactor_*-sinus_,
                                t.scalefactor_*sinus_,
                                t.scalefactor_*cosinus_);
  }

  Aff_transformation_2 compose(const Transformation &t) const
  {
    return Aff_transformation_2(cosinus_*t.t11 + sinus_*t.t12,
                                -sinus_*t.t11 + cosinus_*t.t12,
                                t.t13,
                                cosinus_*t.t21 + sinus_*t.t22,
                                -sinus_*t.t21 + cosinus_*t.t22,
                                t.t23);
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
              case 1: return -sinus_;
              case 2: return FT(0);
            }
    case 1: switch (j)
            {
              case 0: return sinus_;
              case 1: return cosinus_;
              case 2: return FT(0);
            }
    case 2: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(0);
              case 2: return FT(1);
            }
    }
    return FT(0);
  }

  std::ostream &print(std::ostream &os) const
  {
    os << "Aff_transformationC2(" << sinus_ << ", " << cosinus_ <<  ")";
    return os;
  }

private:
  FT sinus_, cosinus_;
};

} //namespace CGAL

#endif // CGAL_CARTESIAN_ROTATION_REP_2_H
