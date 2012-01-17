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
// Author(s)     : Andreas Fabri, Lutz Kettner

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H

#include <cmath>
#include <CGAL/Handle_for_virtual.h>

namespace CGAL {

class Identity_transformation;
template < class R > class Aff_transformation_rep_baseC2;
template < class R > class Aff_transformation_repC2;
template < class R > class Translation_repC2;
template < class R > class Rotation_repC2;
template < class R > class Scaling_repC2;

} //namespace CGAL

#include <CGAL/Cartesian/Aff_transformation_rep_2.h>
#include <CGAL/Cartesian/Translation_rep_2.h>
#include <CGAL/Cartesian/Rotation_rep_2.h>
#include <CGAL/Cartesian/Scaling_rep_2.h>

namespace CGAL {

template < class R_ >
class Aff_transformationC2
  : public Handle_for_virtual< Aff_transformation_rep_baseC2<R_> >
{
  typedef typename R_::FT                   FT;
  typedef Aff_transformation_rep_baseC2<R_> Aff_t_base;

  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  using Handle_for_virtual<Aff_t_base>::initialize_with;

public:
  typedef R_                                R;
   
  Aff_transformationC2()
  {
    initialize_with(Aff_transformation_repC2<R>(FT(1), FT(0), FT(0), FT(1)));
  }

  Aff_transformationC2(const Identity_transformation)
  {
    initialize_with(Aff_transformation_repC2<R>(FT(1), FT(0), FT(0), FT(1)));
  }

  Aff_transformationC2(const Translation, const Vector_2 &v)
  {
    initialize_with(Translation_repC2<R>(v));
  }

  // Rational Rotation:
  Aff_transformationC2(const Rotation,
                       const Direction_2 &d,
                       const FT &num,
                       const FT &den = FT(1))
  {
    initialize_with(Rotation_repC2<R>(d, num, den));
  }

  Aff_transformationC2(const Rotation,
                       const FT &sine,
                       const FT &cosine,
                       const FT &w = FT(1))
  {
    if (w != FT(1))
      initialize_with(Rotation_repC2<R>(sine/w, cosine/w));
    else
      initialize_with(Rotation_repC2<R>(sine, cosine));
  }

  Aff_transformationC2(const Scaling, const FT &s, const FT &w = FT(1))
  {
    if (w != FT(1))
      initialize_with(Scaling_repC2<R>(s/w));
    else
      initialize_with(Scaling_repC2<R>(s));
  }

  // The general case:
  // a 3x2 matrix for the operations combining rotation, scaling, translation
  Aff_transformationC2(const FT & m11, const FT & m12, const FT & m13,
                       const FT & m21, const FT & m22, const FT & m23,
                       const FT &w = FT(1))
  {
    if (w != FT(1))
      initialize_with(Aff_transformation_repC2<R>(m11/w, m12/w, m13/w,
                                                  m21/w, m22/w, m23/w));
    else
      initialize_with(Aff_transformation_repC2<R>(m11, m12, m13,
                                                  m21, m22, m23));
  }

  Aff_transformationC2(const FT & m11, const FT & m12,
                       const FT & m21, const FT & m22,
                       const FT &w = FT(1))
  {
    initialize_with(Aff_transformation_repC2<R>(m11/w, m12/w, m21/w, m22/w));
  }

  Point_2
  transform(const Point_2 &p) const 
  { return this->Ptr()->transform(p); } 

  Point_2
  operator()(const Point_2 &p) const
  { return transform(p); }

  Vector_2
  transform(const Vector_2 &v) const 
  { return this->Ptr()->transform(v); }

  Vector_2
  operator()(const Vector_2 &v) const
  { return transform(v); } // FIXME : not compiled by the test-suite.

  Direction_2
  transform(const Direction_2 &d) const
  { return this->Ptr()->transform(d); }

  Direction_2
  operator()(const Direction_2 &d) const
  { return transform(d); }

  Line_2
  transform(const Line_2 &l) const
  { return l.transform(*this); }

  Line_2
  operator()(const Line_2 &l) const
  { return transform(l); }

  Aff_transformation_2 inverse() const { return this->Ptr()->inverse(); }

  bool is_even() const { return this->Ptr()->is_even(); }
  bool is_odd() const { return ! (this->Ptr()->is_even()); }

  FT cartesian(int i, int j) const { return this->Ptr()->cartesian(i,j); }
  FT homogeneous(int i, int j) const { return cartesian(i,j); }
  FT m(int i, int j) const { return cartesian(i,j); }
  FT hm(int i, int j) const { return cartesian(i,j); }

  Aff_transformation_2 operator*(const Aff_transformationC2 &t) const
  {
    return (*(this->Ptr())) * (*t.Ptr());
  }

  std::ostream &
  print(std::ostream &os) const;
};

template < class R >
std::ostream&
Aff_transformationC2<R>::print(std::ostream &os) const
{
  this->Ptr()->print(os);
  return os;
}

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC2
template < class R >
std::ostream&
operator<<(std::ostream& os, const Aff_transformationC2<R>& t)
{
  t.print(os);
  return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC2

} //namespace CGAL

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_2_H
