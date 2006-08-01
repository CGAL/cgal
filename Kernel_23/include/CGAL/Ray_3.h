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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri, Stefan Schirra
 
#ifndef CGAL_RAY_3_H
#define CGAL_RAY_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Ray_3 : public R_::Kernel_base::Ray_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Line_3                Line_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;
  typedef typename R_::Kernel_base::Ray_3    RRay_3;
public:

  typedef RRay_3 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  Ray_3() {}

  Ray_3(const RRay_3& r)
      : RRay_3(r) {}

  Ray_3(const Point_3& sp, const Point_3& secondp)
    : RRay_3(sp, secondp) {}

  Ray_3(const Point_3& sp, const Vector_3& v)
    : RRay_3(sp, v) {}

  Ray_3(const Point_3& sp, const Direction_3& d)
    : RRay_3(sp, d) {}

  Ray_3(const Point_3& sp, const Line_3& l)
    : RRay_3(sp, l) {}

  Ray_3 transform(const Aff_transformation_3 &t) const
  {
    return Ray_3(t.transform(this->source()),
                 t.transform(this->second_point()));
    // NB : Homogeneous used direction() instead of second_point().
  }

/*
  const Point_3 &   start() const;
  const Point_3 &   source() const
  {
      return get(base).e0;
  }

  Direction_3 direction() const;
  Vector_3    to_vector() const;
  Line_3      supporting_line() const;
  Ray_3       opposite() const;

  bool        is_degenerate() const;
  bool        collinear_has_on(const Point_3 &p) const;
*/

  Point_3 point(int i) const // TODO : use Qrt
  {
    return R().construct_point_on_3_object()(*this, i);
  }

  // FIXME : Use Qrt
  //typename Qualified_result_of<typename R_::Construct_source_3, Ray_3 >::type
  Point_3
  source() const
  {
    return R().construct_source_3_object()(*this);
  }

  Point_3 second_point() const // TODO : use Qrt
  {
    return R().construct_second_point_3_object()(*this);
  }

  bool has_on(const Point_3 &p) const
  {
    return R().has_on_3_object()(*this, p);
  }

  Direction_3
  direction() const
  {
    typename R::Construct_vector_3 construct_vector;
    typename R::Construct_direction_3 construct_direction;
    return construct_direction( construct_vector(source(), second_point()) );
  }

  Ray_3
  opposite() const
  {
    return Ray_3( source(), - direction() );
  }

  Vector_3
  to_vector() const
  {
    typename R::Construct_vector_3 construct_vector;
    return construct_vector(source(), second_point());
  }

  Line_3
  supporting_line() const
  {
    return R().construct_line_3_object()(source(), second_point());
  }

  bool is_degenerate() const
  {
    return R().is_degenerate_3_object()(*this);
  }

};

#ifndef CGAL_NO_OSTREAM_INSERT_RAY_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Ray_3<R>& r)
{
  typedef typename  R::Kernel_base::Ray_3  RRay_3;
  return os << static_cast<const RRay_3&>(r);
}
#endif // CGAL_NO_OSTREAM_INSERT_RAY_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAY_3
template < class R >
std::istream&
operator>>(std::istream& is, Ray_3<R>& r)
{
  typedef typename  R::Kernel_base::Ray_3  RRay_3;
  return is >> static_cast<RRay_3&>(r);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAY_3

CGAL_END_NAMESPACE

#endif // CGAL_RAY_3_H
