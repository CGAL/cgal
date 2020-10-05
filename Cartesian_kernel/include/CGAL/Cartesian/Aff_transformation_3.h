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
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_3_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_3_H

#include <cmath>
#include <CGAL/Handle_for_virtual.h>

namespace CGAL {

class Identity_transformation;
template <class R> class Aff_transformation_rep_baseC3;
template <class R> class Aff_transformation_repC3;
template <class R> class Translation_repC3;
template <class R> class Scaling_repC3;

} //namespace CGAL

#include <CGAL/Cartesian/Aff_transformation_rep_3.h>
#include <CGAL/Cartesian/Translation_rep_3.h>
#include <CGAL/Cartesian/Scaling_rep_3.h>

namespace CGAL {

template < class R_ >
class Aff_transformationC3
  : public Handle_for_virtual<Aff_transformation_rep_baseC3<R_> >
{
  friend class PlaneC3<R_>; // FIXME: why ?

  typedef typename R_::FT                   FT;
  typedef Aff_transformation_rep_baseC3<R_> Aff_t_base;

  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  using Handle_for_virtual<Aff_t_base>::initialize_with;
public:
  typedef R_                               R;

  Aff_transformationC3()
  {
    FT ft1(1), ft0(0);
    initialize_with(Aff_transformation_repC3<R>(ft1, ft0, ft0,
                                          ft0, ft1, ft0,
                                          ft0, ft0, ft1));
  }

  Aff_transformationC3(const Identity_transformation)
  {
    FT ft1(1), ft0(0);
    initialize_with(Aff_transformation_repC3<R>(ft1, ft0, ft0,
                                          ft0, ft1, ft0,
                                          ft0, ft0, ft1));
  }

  Aff_transformationC3(const Translation, const Vector_3 &v)
  {
    initialize_with(Translation_repC3<R>(v));
  }

  Aff_transformationC3(const Scaling, const FT &s, const FT &w = FT(1))
  {
    if (w != FT(1))
      initialize_with(Scaling_repC3<R>(s/w));
    else
      initialize_with(Scaling_repC3<R>(s));
  }

  // General form: without translation
  Aff_transformationC3(const FT& m11, const FT& m12, const FT& m13,
                       const FT& m21, const FT& m22, const FT& m23,
                       const FT& m31, const FT& m32, const FT& m33,
                       const FT& w = FT(1))
  {
    if (w != FT(1))
      initialize_with(Aff_transformation_repC3<R>(m11/w, m12/w, m13/w,
                                            m21/w, m22/w, m23/w,
                                            m31/w, m32/w, m33/w));
    else
      initialize_with(Aff_transformation_repC3<R>(m11, m12, m13,
                                            m21, m22, m23,
                                            m31, m32, m33));
  }

  // General form: with translation
  Aff_transformationC3(
              const FT& m11, const FT& m12, const FT& m13, const FT& m14,
              const FT& m21, const FT& m22, const FT& m23, const FT& m24,
              const FT& m31, const FT& m32, const FT& m33, const FT& m34,
              const FT& w = FT(1))
  {
    if (w != FT(1))
      initialize_with(Aff_transformation_repC3<R>(m11/w, m12/w, m13/w, m14/w,
                                            m21/w, m22/w, m23/w, m24/w,
                                            m31/w, m32/w, m33/w, m34/w));
    else
      initialize_with(Aff_transformation_repC3<R>(m11, m12, m13, m14,
                                            m21, m22, m23, m24,
                                            m31, m32, m33, m34));
  }

  Point_3
  transform(const Point_3 &p) const
  { return this->Ptr()->transform(p); }

  Point_3
  operator()(const Point_3 &p) const
  { return transform(p); }

  Vector_3
  transform(const Vector_3 &v) const
  { return this->Ptr()->transform(v); }

  Vector_3
  operator()(const Vector_3 &v) const
  { return transform(v); }

  Direction_3
  transform(const Direction_3 &d) const
  { return this->Ptr()->transform(d); }

  Direction_3
  operator()(const Direction_3 &d) const
  { return transform(d); }

  Plane_3
  transform(const Plane_3& p) const
  {
    if (is_even())
      return Plane_3(transform(p.point()),
                 transpose().inverse().transform(p.orthogonal_direction()));
    else
      return Plane_3(transform(p.point()),
               - transpose().inverse().transform(p.orthogonal_direction()));
  }

  Plane_3
  operator()(const Plane_3& p) const
  { return transform(p); } // FIXME : not compiled by the test-suite !

  Aff_transformation_3 inverse() const { return this->Ptr()->inverse(); }

  bool is_even() const { return this->Ptr()->is_even(); }
  bool is_odd() const { return  ! (this->Ptr()->is_even()); }

  FT cartesian(int i, int j) const { return this->Ptr()->cartesian(i,j); }
  FT homogeneous(int i, int j) const { return cartesian(i,j); }
  FT m(int i, int j) const { return cartesian(i,j); }
  FT hm(int i, int j) const { return cartesian(i,j); }

  Aff_transformation_3 operator*(const Aff_transformationC3 &t) const
  { return (*this->Ptr()) * (*t.Ptr()); }

  std::ostream &
  print(std::ostream &os) const;

  bool operator==(const Aff_transformationC3 &t)const
  {
    for(int i=0; i<3; ++i)
      for(int j = 0; j< 4; ++j)
        if(cartesian(i,j)!=t.cartesian(i,j))
          return false;
    return true;
  }

  bool operator!=(const Aff_transformationC3 &t)const
  {
    return !(*this == t);
  }

protected:
  Aff_transformation_3  transpose() const { return this->Ptr()->transpose(); }
};


template < class R >
std::ostream&
Aff_transformationC3<R>::print(std::ostream &os) const
{
  this->Ptr()->print(os);
  return os;
}

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC3
template < class R >
std::ostream&
operator<<(std::ostream &os, const Aff_transformationC3<R> &t)
{
  t.print(os);
  return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC3

} //namespace CGAL

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_3_H
