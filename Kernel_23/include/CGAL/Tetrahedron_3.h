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

#ifndef CGAL_TETRAHEDRON_3_H
#define CGAL_TETRAHEDRON_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Tetrahedron_3 : public R_::Kernel_base::Tetrahedron_3
{
  typedef typename R_::Point_3             Point_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;
  typedef typename R_::Kernel_base::Tetrahedron_3  RTetrahedron_3;
public:

  typedef RTetrahedron_3 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  Tetrahedron_3() {}

  Tetrahedron_3(const RTetrahedron_3& t)
      : RTetrahedron_3(t) {}

  Tetrahedron_3(const Point_3& p,
                const Point_3& q,
                const Point_3& r,
                const Point_3& s)
    : RTetrahedron_3(p,q,r,s) {}

  Tetrahedron_3 transform(const Aff_transformation_3 &t) const
  {
    return Tetrahedron_3(t.transform(this->vertex(0)),
                         t.transform(this->vertex(1)),
                         t.transform(this->vertex(2)),
                         t.transform(this->vertex(3)));
  }

  // FIXME TODO : Why doesn't Qrt work here ???
  //typename Qualified_result_of<typename R::Construct_vertex_3, Tetrahedron_3, int>::type
  Point_3
  vertex(int i) const
  {
    return R().construct_vertex_3_object()(*this,i);
  }

  //typename Qualified_result_of<typename R::Construct_vertex_3, Tetrahedron_3, int>::type
  Point_3
  operator[](int i) const
  {
    return vertex(i);
  }

  bool
  is_degenerate() const
  {
    return R().is_degenerate_3_object()(*this);
  }

  Orientation orientation() const
  {
    return R().orientation_3_object()(*this);
  }

  Bounded_side bounded_side(const Point_3 &p) const
  {
    return R().bounded_side_3_object()(*this, p);
  }

  Oriented_side oriented_side(const Point_3 &p) const
  {
    return R().oriented_side_3_object()(*this, p);
  }

  bool has_on_positive_side(const Point_3 &p) const
  {
    return R().has_on_positive_side_3_object()(*this, p);
  }

  bool has_on_negative_side(const Point_3 &p) const
  {
    return R().has_on_negative_side_3_object()(*this, p);
  }

  bool has_on_boundary(const Point_3 &p) const
  {
    return R().has_on_boundary_3_object()(*this, p);
  }

  bool has_on_bounded_side(const Point_3 &p) const
  {
    return R().has_on_bounded_side_3_object()(*this, p);
  }

  bool has_on_unbounded_side(const Point_3 &p) const
  {
    return R().has_on_unbounded_side_3_object()(*this, p);
  }

  typename Qualified_result_of<typename R::Compute_volume_3, Tetrahedron_3>::type
  volume() const
  {
    return R().compute_volume_3_object()(*this);
  }

};

#ifndef CGAL_NO_OSTREAM_INSERT_TETRAHEDRON_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Tetrahedron_3<R>& t)
{
  typedef typename  R::Kernel_base::Tetrahedron_3  RTetrahedron_3;
  return os << static_cast<const RTetrahedron_3&>(t);
}
#endif // CGAL_NO_OSTREAM_INSERT_TETRAHEDRON_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRON_3
template < class R >
std::istream&
operator>>(std::istream& is, Tetrahedron_3<R>& t)
{
  typedef typename  R::Kernel_base::Tetrahedron_3  RTetrahedron_3;
  return is >> static_cast<RTetrahedron_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRON_3

CGAL_END_NAMESPACE

#endif  // CGAL_TETRAHEDRON_3_H
