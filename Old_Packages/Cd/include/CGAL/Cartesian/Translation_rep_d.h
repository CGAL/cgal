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
// file          : include/CGAL/Cartesian/Translation_rep_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_TRANSLATION_REP_D_H
#define CGAL_CARTESIAN_TRANSLATION_REP_D_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Translation_repCd : public Aff_transformation_rep_baseCd<R>
{
  friend class Aff_transformation_repCd<R>;
  friend class Scaling_repCd<R>;

public:
  typedef typename R::FT                                FT;
  typedef typename R::RT                                RT;
  typedef Aff_transformation_rep_baseCd<R>              Transformation_base_d;
  typedef Aff_transformation_repCd<R>                   Transformation_d;
  typedef Translation_repCd<R>                          Translation_d;
  typedef Scaling_repCd<R>                              Scaling_d;
  typedef Homothecy_repCd<R>                            Homothecy_d;
  typedef Reflexion_repCd<R>                            Reflexion_d;
  typedef typename Transformation_base_d::Point_d       Point_d;
  typedef typename Transformation_base_d::Vector_d      Vector_d;
  typedef typename Transformation_base_d::Direction_d   Direction_d;
  typedef typename Transformation_base_d::Aff_transformation_d
	                                                Aff_transformation_d;

  Translation_repCd() {}
  Translation_repCd(const Vector_d &tv) : _translationvector(tv) {}
  virtual ~Translation_repCd() {}

  virtual Point_d     transform(const Point_d &p) const
  {
    return p + _translationvector;
  }

  virtual Vector_d    transform(const Vector_d &v) const
  {
    return v;
  }

  virtual Direction_d transform(const Direction_d &d) const
  {
    return d;
  }

  virtual Aff_transformation_d operator*(const Transformation_base_d &t) const
  {
    return t.compose(*this);
  }

  virtual Aff_transformation_d compose(const Transformation_d &t) const
  {
  }

  virtual Aff_transformation_d compose(const Translation_d &t) const
  {
    return Aff_transformation_d(TRANSLATION,
                                _translationvector + t._translationvector);
  }

  virtual Aff_transformation_d compose(const Scaling_d &t) const
  {
  }

  virtual Aff_transformation_d inverse() const
  {
  }

  virtual Aff_transformation_d transpose() const
  {
  }
  
  virtual bool is_even() const
  {
    return true;
  }

  virtual FT cartesian(int i, int j) const
  {
    if (j==i) return FT(1);
    if (j==dimension()) return _translationvector[i];
    return FT(0);
  }

private:
  Vector_d   _translationvector;
};
CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TRANSLATION_REP_D_H
