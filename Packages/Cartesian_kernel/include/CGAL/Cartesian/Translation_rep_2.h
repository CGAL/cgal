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
// file          : include/CGAL/Cartesian/Translation_rep_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_TRANSLATION_REP_2_H
#define CGAL_CARTESIAN_TRANSLATION_REP_2_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Translation_repC2 : public Aff_transformation_rep_baseC2<R>
{
friend class Aff_transformation_repC2<R>;
friend class Rotation_repC2<R>;
friend class Scaling_repC2<R>;

public:
  typedef typename R::FT                         FT;
  typedef typename R::RT                         RT;
  typedef Aff_transformation_rep_baseC2<R>       Aff_t_base;
  typedef Aff_transformation_repC2<R>            Transformation;
  typedef Translation_repC2<R>                   Translation;
  typedef Rotation_repC2<R>                      Rotation;
  typedef Scaling_repC2<R>                       Scaling;
  typedef typename Aff_t_base::Point_2           Point_2;
  typedef typename Aff_t_base::Vector_2          Vector_2;
  typedef typename Aff_t_base::Direction_2       Direction_2;
  typedef typename Aff_t_base::Aff_transformation_2 Aff_transformation_2;

  Translation_repC2()
  {}

  Translation_repC2(const Vector_2 &tv)
    : translationvector_(tv)
  {}

  Point_2     transform(const Point_2 &p) const
                            { return p + translationvector_; }
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

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TRANSLATION_REP_2_H
