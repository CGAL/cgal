// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Translation_rep_2.h
// source        : include/CGAL/Cartesian/Translation_rep_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr, Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


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
    : _translationvector(tv)
  {}

  Point_2     transform(const Point_2 &p) const
                            { return p + _translationvector; }
  Vector_2    transform(const Vector_2 &v) const { return v; }
  Direction_2 transform(const Direction_2 &d) const { return d; }

  Aff_transformation_2 operator*(const Aff_t_base &t) const
  {
    return t.compose(*this);
  }

  Aff_transformation_2 compose(const Translation &t) const
  {
    return Aff_transformation_2(TRANSLATION,
                                _translationvector + t._translationvector);
  }

  Aff_transformation_2 compose(const Rotation &t) const
  {
    return Aff_transformation_2(t._cosinus,
                                -t._sinus,
                                t._cosinus*_translationvector.x() -
                                t._sinus*_translationvector.y(),

                                t._sinus,
                                t._cosinus,
                                t._sinus*_translationvector.x() +
                                t._cosinus*_translationvector.y());
  }

  Aff_transformation_2 compose(const Scaling &t) const
  {
    return Aff_transformation_2(t._scalefactor,
                                FT(0),
                                t._scalefactor*_translationvector.x(),

                                FT(0),
                                t._scalefactor,
                                t._scalefactor*_translationvector.y());
  }

  Aff_transformation_2 compose(const Transformation &t) const
  {
    return Aff_transformation_2(t.t11,
                                t.t12,
                                t.t11 * _translationvector.x()
                                + t.t12 * _translationvector.y()
                                + t.t13,

                                t.t21,
                                t.t22,
                                t.t21 * _translationvector.x()
                                + t.t22*_translationvector.y()
                                + t.t23);
  }

  Aff_transformation_2 inverse() const
  {
    return Aff_transformation_2(TRANSLATION, - _translationvector);
  }

  bool         is_even() const
  {
    return true;
  }

  FT cartesian(int i, int j) const
  {
    if (j==i) return FT(1);
    if (j==2) return _translationvector[i];
    return FT(0);
  }

  std::ostream &print(std::ostream &os) const
  {
    os << "Aff_transformationC2(VectorC2(" << _translationvector.x() << ", "
       << _translationvector.y()  <<  "))";
    return os;
  }

private:
  Vector_2   _translationvector;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TRANSLATION_REP_2_H
