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

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_2_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_2_H

#include <CGAL/determinant.h>
#include <CGAL/Handle_for_virtual.h>
#include <CGAL/Cartesian/Aff_transformation_2.h>

namespace CGAL {

template < class R >
class Aff_transformation_rep_baseC2
  : public Ref_counted_virtual
{
public:
  typedef typename R::FT                   FT;
  typedef typename R::Point_2              Point_2;
  typedef typename R::Vector_2             Vector_2;
  typedef typename R::Direction_2          Direction_2;
  typedef typename R::Aff_transformation_2 Aff_transformation_2;

  virtual ~Aff_transformation_rep_baseC2() {}

  virtual Point_2     transform(const Point_2 &p) const  = 0;
  virtual Vector_2    transform(const Vector_2 &v) const = 0;
  virtual Direction_2 transform(const Direction_2 &d) const=0;

  virtual Aff_transformation_2 operator*(
                       const Aff_transformation_rep_baseC2 &t) const = 0;

  virtual Aff_transformation_2 compose(
                       const Aff_transformation_repC2<R> &t) const  = 0;

  virtual Aff_transformation_2 compose(
                       const Translation_repC2<R> &t) const  = 0;

  virtual Aff_transformation_2 compose(
                       const Rotation_repC2<R> &t) const  = 0;

  virtual Aff_transformation_2 compose(
                       const Scaling_repC2<R> &t) const  = 0;

  virtual Aff_transformation_2 inverse() const  = 0;
  virtual bool                 is_even() const  = 0;
  virtual FT                   cartesian(int i, int j) const = 0;
  virtual std::ostream         &print(std::ostream &os) const = 0;
};

template < class R >
class Aff_transformation_repC2
  : public Aff_transformation_rep_baseC2<R>
{
public:
  typedef typename R::FT                              FT;
  typedef Aff_transformation_repC2<R>                 Self;
  typedef Aff_transformation_rep_baseC2<R>            Aff_t_base;
  typedef typename Aff_t_base::Point_2                Point_2;
  typedef typename Aff_t_base::Vector_2               Vector_2;
  typedef typename Aff_t_base::Direction_2            Direction_2;
  typedef typename Aff_t_base::Aff_transformation_2   Aff_transformation_2;

friend class Translation_repC2<R>;
friend class Rotation_repC2<R>;
friend class Scaling_repC2<R>;

  Aff_transformation_repC2()
  {}

  Aff_transformation_repC2( const FT& m11, const FT& m12,
                            const FT& m21, const FT& m22)
    : t11(m11), t12(m12), t13(0),
      t21(m21), t22(m22), t23(0)
  {}

  Aff_transformation_repC2( const FT& m11, const FT& m12, const FT& m13,
                            const FT& m21, const FT& m22, const FT& m23)
    : t11(m11), t12(m12), t13(m13),
      t21(m21), t22(m22), t23(m23)
  {}

  Point_2 transform(const Point_2& p) const
  {
    typename R::Construct_point_2 construct_point_2;
    return construct_point_2(t11 * p.x() + t12 * p.y() + t13,
			     t21 * p.x() + t22 * p.y() + t23);
  }

  // note that a vector is not translated
  Vector_2 transform(const Vector_2& v) const
  {
    return Vector_2(t11 * v.x() + t12 * v.y(),
                    t21 * v.x() + t22 * v.y());
  }

  // note that a direction is not translated
  Direction_2 transform(const Direction_2& dir) const
  {
    return Direction_2(t11 * dir.dx() + t12 * dir.dy(),
                       t21 * dir.dx() + t22 * dir.dy());
  }

  // Note that Aff_transformation is not defined yet,
  // so the following 6 functions have to be implemented later...
  Aff_transformation_2 inverse() const;
  Aff_transformation_2 operator*(const Aff_t_base &t) const;
  Aff_transformation_2 compose(const Self &t) const;
  Aff_transformation_2 compose(const Translation_repC2<R> &t) const;
  Aff_transformation_2 compose(const Rotation_repC2<R> &t) const;
  Aff_transformation_2 compose(const Scaling_repC2<R> &t) const;

  bool is_even() const
  {
    return sign_of_determinant(t11, t12, t21, t22) == POSITIVE;
  }

  FT cartesian(int i, int j) const
  {
    switch (i)
    {
    case 0: switch (j)
            {
              case 0: return t11;
              case 1: return t12;
              case 2: return t13;
            }
    case 1: switch (j)
            {
              case 0: return t21;
              case 1: return t22;
              case 2: return t23;
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
    os <<"Aff_transformationC2(" <<t11<<" "<<t12<<" "<<t13<<std::endl;
    os <<"                     " <<t21<<" "<<t22<<" "<<t23<<")";
    return os;
  }

private:
    FT   t11, t12, t13;
    FT   t21, t22, t23;
};

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
inverse() const
{
  FT det = FT(1) / (t11 * t22 - t12 * t21);
  return Aff_transformation_2(
    det * t22,    det * (-t12), det * (t12*t23-t13*t22),
    det * (-t21), det * t11 ,   det * (t13*t21-t11*t23));
}

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
operator*(const Aff_transformation_rep_baseC2<R> &t) const
{
  return t.compose(*this);
}
 
template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
compose(const Aff_transformation_repC2<R> &t) const
{
  return Aff_transformation_2(t.t11*t11 + t.t12*t21,
                              t.t11*t12 + t.t12*t22,
                              t.t11*t13 + t.t12*t23 + t.t13,
                              t.t21*t11 + t.t22*t21,
                              t.t21*t12 + t.t22*t22,
                              t.t21*t13 + t.t22*t23 + t.t23 );
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
compose(const Translation_repC2<R> &t) const
{
  return Aff_transformation_2(t11,
                              t12,
                              t13 + t.translationvector_.x(),
                              t21,
                              t22,
                              t23 + t.translationvector_.y());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
compose(const Rotation_repC2<R> &t) const
{
  return Aff_transformation_2(t.cosinus_*t11 - t.sinus_*t21,
                              t.cosinus_*t12 - t.sinus_*t22,
                              t.cosinus_*t13 - t.sinus_*t23,
                              t.sinus_*t11 + t.cosinus_*t21,
                              t.sinus_*t12 + t.cosinus_*t22,
                              t.sinus_*t13 + t.cosinus_*t23);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
compose(const Scaling_repC2<R> &t) const
{
   return Aff_transformation_2(t.scalefactor_ * t11,
                               t.scalefactor_ * t12,
                               t.scalefactor_ * t13,
                               t.scalefactor_ * t21,
                               t.scalefactor_ * t22,
                               t.scalefactor_ * t23);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_2_H
