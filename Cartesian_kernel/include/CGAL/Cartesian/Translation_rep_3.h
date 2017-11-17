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
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_TRANSLATION_REP_3_H
#define CGAL_CARTESIAN_TRANSLATION_REP_3_H

namespace CGAL {

template < class R >
class Translation_repC3 : public Aff_transformation_rep_baseC3<R>
{
  friend class Aff_transformation_repC3<R>;
  friend class Scaling_repC3<R>;

public:
  typedef typename R::FT                                FT;
  typedef Aff_transformation_rep_baseC3<R>              Transformation_base_3;
  typedef Aff_transformation_repC3<R>                   Transformation_3;
  typedef Translation_repC3<R>                          Translation_3;
  typedef Scaling_repC3<R>                              Scaling_3;
  typedef typename Transformation_base_3::Point_3       Point_3;
  typedef typename Transformation_base_3::Vector_3      Vector_3;
  typedef typename Transformation_base_3::Direction_3   Direction_3;
  typedef typename Transformation_base_3::Aff_transformation_3
	                                                Aff_transformation_3;

  Translation_repC3() {}
  Translation_repC3(const Vector_3 &tv) : translationvector_(tv) {}
  virtual ~Translation_repC3() {}

  virtual Point_3     transform(const Point_3 &p) const
  {
    return p + translationvector_;
  }

  virtual Vector_3    transform(const Vector_3 &v) const
  {
    return v;
  }

  virtual Direction_3 transform(const Direction_3 &d) const
  {
    return d;
  }

  virtual Aff_transformation_3 operator*(const Transformation_base_3 &t) const
  {
    return t.compose(*this);
  }

  virtual Aff_transformation_3 compose(const Transformation_3 &t) const
  {
    return Aff_transformation_3(t.t11,
                                t.t12,
				t.t13,
				t.t11 * translationvector_.x()
				+ t.t12 * translationvector_.y()
				+ t.t13 * translationvector_.z() + t.t14,
				
				t.t21,
                                t.t22,
				t.t23,
				t.t21 * translationvector_.x()
				+ t.t22 * translationvector_.y()
				+ t.t23 * translationvector_.z() + t.t24,
				
				t.t31,
                                t.t32,
				t.t33,
				t.t31 * translationvector_.x()
				+ t.t32 * translationvector_.y()
				+ t.t33 * translationvector_.z() + t.t34);
  }

  virtual Aff_transformation_3 compose(const Translation_3 &t) const
  {
    return Aff_transformation_3(TRANSLATION,
                                translationvector_ + t.translationvector_);
  }

  virtual Aff_transformation_3 compose(const Scaling_3 &t) const
  {
    FT ft0(0);
    return Aff_transformation_3(t.scalefactor_,
                                ft0,
				ft0,
				t.scalefactor_ * translationvector_.x(),
				
				ft0,
                                t.scalefactor_,
				ft0,
				t.scalefactor_ * translationvector_.y(),
				
				ft0,
                                ft0,
				t.scalefactor_,
				t.scalefactor_ * translationvector_.z());
  }

  virtual Aff_transformation_3 inverse() const
  {
    return Aff_transformation_3(TRANSLATION, - translationvector_);
  }

  virtual Aff_transformation_3 transpose() const
  {
    return Aff_transformation_3(TRANSLATION, translationvector_);
  }
  
  virtual bool is_even() const
  {
    return true;
  }

  virtual FT cartesian(int i, int j) const
  {
    if (j==i) return FT(1);
    if (j==3) return translationvector_[i];
    return FT(0);
  }

  virtual std::ostream &print(std::ostream &os) const
  {
    os << "Aff_transformationC3(VectorC3("<< translationvector_.x() << ","
       << translationvector_.y() << ","
       << translationvector_.z() << "))\n";
    return os;
  }

private:
  Vector_3   translationvector_;
};

} //namespace CGAL

#endif // CGAL_CARTESIAN_TRANSLATION_REP_3_H
