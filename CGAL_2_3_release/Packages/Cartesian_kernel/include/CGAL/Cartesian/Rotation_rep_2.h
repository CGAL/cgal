// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Rotation_rep_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_ROTATION_REP_2_H
#define CGAL_CARTESIAN_ROTATION_REP_2_H

#include <CGAL/rational_rotation.h>

CGAL_BEGIN_NAMESPACE

template < class R >
class Rotation_repC2: public Aff_transformation_rep_baseC2<R>
{
friend class Aff_transformation_repC2<R>;
friend class Translation_repC2<R>;
friend class Scaling_repC2<R>;

public:
  typedef Aff_transformation_rep_baseC2<R> Aff_t_base;
  typedef typename Aff_t_base::FT                FT;
  typedef typename Aff_t_base::RT                RT;
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
  { // FIXME : construction
    FT sin_num;
    FT cos_num;
    FT denom;

    rational_rotation_approximation(d.to_vector().x(),
                                    d.to_vector().y(),
                                    sin_num,
                                    cos_num,
                                    denom,
                                    eps_num,
                                    eps_den);
    sinus_   = sin_num/denom;
    cosinus_ = cos_num/denom;
  }

  Point_2      transform(const Point_2 &p) const
  { // FIXME : construction
    return Point_2(cosinus_ * p.x() - sinus_ * p.y(),
                   sinus_ * p.x() + cosinus_ * p.y());
  }

  Vector_2     transform(const Vector_2 &v) const
  { // FIXME : construction
    return Vector_2(cosinus_ * v.x() - sinus_ * v.y(),
                    sinus_ * v.x() + cosinus_ * v.y());
  }

  Direction_2  transform(const Direction_2 &d) const
  { // FIXME : construction
    Vector_2  v = d.to_vector();
    return Direction_2(cosinus_ * v.x() - sinus_ * v.y(),
                       sinus_ * v.x() + cosinus_ * v.y());
  }

  Aff_transformation_2 inverse() const
  { // FIXME : construction
    return Aff_transformation_2(ROTATION, - sinus_, cosinus_, FT(1));
  }

  Aff_transformation_2 operator*(const Aff_t_base &t) const
  {
    return t.compose(*this);
  }

  Aff_transformation_2 compose(const Translation &t) const
  { // FIXME : construction
    return Aff_transformation_2(cosinus_,
                                -sinus_,
                                t.translationvector_.x(),
                                sinus_,
                                cosinus_,
                                t.translationvector_.y());
  }

  Aff_transformation_2 compose(const Rotation &t) const
  { // FIXME : construction
    return Aff_transformation_2(ROTATION,
                                t.sinus_*cosinus_ + t.cosinus_*sinus_,
                                t.cosinus_*cosinus_-t.sinus_*sinus_ );
  }

  Aff_transformation_2 compose(const Scaling &t) const
  { // FIXME : construction
    return Aff_transformation_2(t.scalefactor_*cosinus_,
                                t.scalefactor_*-sinus_,
                                t.scalefactor_*sinus_,
                                t.scalefactor_*cosinus_);
  }

  Aff_transformation_2 compose(const Transformation &t) const
  { // FIXME : construction
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
  { // FIXME : construction
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

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_ROTATION_REP_2_H
