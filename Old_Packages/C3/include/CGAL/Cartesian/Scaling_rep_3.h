// ==========================================================================
//
// Copyright (c) 1998, 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Scaling_rep_3.h
// source        : include/CGAL/Cartesian/Scaling_rep_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef _MSC_VER
#define typename
#endif

#ifndef CGAL_CARTESIAN_SCALING_REP_3_H
#define CGAL_CARTESIAN_SCALING_REP_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
class Scaling_repC3 : public Aff_transformation_rep_baseC3<R>
{
  friend class Aff_transformation_repC3<R>;
  friend class Translation_repC3<R>;

public:
  typedef typename R::FT                                FT;
  typedef typename R::RT                                RT;
  typedef Aff_transformation_rep_baseC3<R>              Transformation_base_3;
  typedef Aff_transformation_repC3<R>                   Transformation_3;
  typedef Translation_repC3<R>                          Translation_3;
  typedef Scaling_repC3<R>                              Scaling_3;
  typedef typename Transformation_base_3::Point_3       Point_3;
  typedef typename Transformation_base_3::Vector_3      Vector_3;
  typedef typename Transformation_base_3::Direction_3   Direction_3;
  typedef typename Transformation_base_3::Aff_transformation_3
	                                                Aff_transformation_3;

  Scaling_repC3() {}
  Scaling_repC3(const FT &s) : _scalefactor(s) {}
  virtual ~Scaling_repC3() {}

  virtual Point_3      transform(const Point_3 &p) const
  {
    return Point_3(_scalefactor * p.x(),
                   _scalefactor * p.y(),
                   _scalefactor * p.z());
  }

  virtual Vector_3     transform(const Vector_3 &v) const
  {
    return Vector_3(_scalefactor * v.x(), _scalefactor * v.y(),
                    _scalefactor * v.z());
  }

  virtual Direction_3  transform(const Direction_3 &d) const
  {
    return d;
  }

  virtual Aff_transformation_3 operator*(const Transformation_base_3 &t) const
  {
    return t.compose(*this);
  }

  virtual Aff_transformation_3 compose(const Transformation_3 &t) const
  {
    return Aff_transformation_3(_scalefactor * t.t11,
                                _scalefactor * t.t12,
                                _scalefactor * t.t13,
                                t.t14,
				
                                _scalefactor * t.t21,
                                _scalefactor * t.t22,
                                _scalefactor * t.t23,
                                t.t24,
				
                                _scalefactor * t.t31,
                                _scalefactor * t.t32,
                                _scalefactor * t.t33,
                                t.t34);
  }

  virtual Aff_transformation_3 compose(const Translation_3 &t) const
  {
    FT ft0(0);
    return Aff_transformation_3(_scalefactor,
                                ft0,
                                ft0,
				t._translationvector.x(),
				
                                ft0,
				_scalefactor,
                                ft0,
				t._translationvector.y(),
				
                                ft0,
                                ft0,
				_scalefactor,
				t._translationvector.z());
  }

  virtual Aff_transformation_3 compose(const Scaling_3 &t) const
  {
    return Aff_transformation_3(SCALING, _scalefactor*t._scalefactor);
  }

  virtual Aff_transformation_3 inverse() const
  {
    return Aff_transformation_3(SCALING, FT(1)/_scalefactor);
  }

  virtual Aff_transformation_3 transpose() const
  {
    return Aff_transformation_3(SCALING, _scalefactor);
  }

  virtual bool            is_even() const
  {
    return true;
  }

  virtual FT cartesian(int i, int j) const
  {
    if (i!=j) return FT(0);
    if (i==3) return FT(1);
    return _scalefactor;
  }

  virtual std::ostream &print(std::ostream &os) const
  {
    FT ft0(0);
    os << "Aff_transformationC3(" << _scalefactor << ")";
    return os;
  }

private:
  FT   _scalefactor;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SCALING_REP_3_H
