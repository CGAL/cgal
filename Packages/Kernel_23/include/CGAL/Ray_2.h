// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_RAY_2_H
#define CGAL_RAY_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Ray_2 : public R_::Kernel_base::Ray_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Direction_2           Direction_2;
  typedef typename R_::Vector_2              Vector_2;
  typedef typename R_::Line_2                Line_2;
  typedef typename R_::Aff_transformation_2  Aff_transformation_2;

  typedef typename R_::Kernel_base::Ray_2    RRay_2;
public:

  typedef RRay_2 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_   R;

  Ray_2() 
    : RRay_2(typename R::Construct_ray_2()().rep())
   {}

  Ray_2(const RRay_2& r)
    : RRay_2(r) {}

  Ray_2(const Point_2 &sp, const Point_2 &secondp)
    : RRay_2(typename R::Construct_ray_2()(sp, secondp).rep()) {}

  Ray_2(const Point_2 &sp, const Direction_2 &d)
    : RRay_2(typename R::Construct_ray_2()(sp, d).rep()) {}

  Ray_2(const Point_2 &sp, const Vector_2 &v)
    : RRay_2(typename R::Construct_ray_2()(sp, v).rep()) {}

  Ray_2(const Point_2 &sp, const Line_2 &l)
    : RRay_2(typename R::Construct_ray_2()(sp, l).rep()) {}


  const Point_2 &     
  source() const
  {
    return R().construct_source_2_object()(*this);
  }


  const Point_2 &
  second_point() const
  {
    return R().construct_second_point_2_object()(*this);
  }


  Point_2
  point(int i) const
  {
    CGAL_kernel_precondition( i >= 0 );
    
    typename R::Construct_vector_2 construct_vector;
    typename R::Construct_scaled_vector_2 construct_scaled_vector;
    typename R::Construct_translated_point_2 construct_translated_point;
    if (i == 0) return source();
    if (i == 1) return second_point();
    return construct_translated_point(source(),
				      construct_scaled_vector(construct_vector(source(), 
									       second_point()),
							      FT(i)));
  }


  const Point_2 &
  start() const
  {
    return source();
  }



  bool is_horizontal() const
  {
    return R().equal_y_2_object()(source(), second_point());
  }

  bool is_vertical() const
  {
    return R().equal_x_2_object()(source(), second_point());
  }

  bool is_degenerate() const
  {
    return source() == second_point();
  }

  Direction_2
  direction() const
  {
    typename R::Construct_vector_2 construct_vector;
    typename R::Construct_direction_2 construct_direction;
    return construct_direction( construct_vector(source(), second_point()) );
  }


  Vector_2
  to_vector() const
  {
    typename R::Construct_vector_2 construct_vector;
    return construct_vector(source(), second_point());
  }

  bool
  has_on(const Point_2 &p) const
  {
    typename R::Construct_vector_2  construct_vector;
    return p == source()
      || R().collinear_2_object()(source(), p, second_point())
      && Direction_2(construct_vector( source(), p)) == direction();
  }



  bool
  collinear_has_on(const Point_2 &p) const
  {
    return R().collinear_has_on_2_object()(*this, p);
  }

  Ray_2
  opposite() const
  {
    return Ray_2( source(), - direction() );
  }


  Line_2
  supporting_line() const
  {
    return R().construct_line_2_object()(source(), second_point());
  }
  
  

  bool
  operator==(const Ray_2& r) const
  {
    return R().equal_2_object()(*this, r);
  }


  bool
  operator!=(const Ray_2& r) const
  {
    return !(*this == r);
  }

  Ray_2 
  transform(const Aff_transformation_2 &t) const
  {
    return rep().transform(t);
  }


};

#ifndef CGAL_NO_OSTREAM_INSERT_RAY_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Ray_2<R> &r)
{
  return os << r.rep();
}
#endif // CGAL_NO_OSTREAM_INSERT_RAY_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAY_2
template < class R >
std::istream &
operator>>(std::istream &is, Ray_2<R> &r)
{
  return is >> r.rep();
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAY_2



CGAL_END_NAMESPACE

#endif  // CGAL_RAY_2_H
