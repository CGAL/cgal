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

#ifndef CGAL_ISO_RECTANGLE_2_H
#define CGAL_ISO_RECTANGLE_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Iso_rectangle_2 : public R_::Kernel_base::Iso_rectangle_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::FT                    FT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Aff_transformation_2  Aff_transformation_2;
  typedef typename R_::Kernel_base::Iso_rectangle_2  RIso_rectangle_2;
public:
  typedef RIso_rectangle_2 Rep;

  typedef typename R_::Cartesian_coordinate_type Cartesian_coordinate_type;
  typedef typename R_::Homogeneous_coordinate_type Homogeneous_coordinate_type;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_   R;

  Iso_rectangle_2() 
    : RIso_rectangle_2(typename R::Construct_iso_rectangle_2()().rep())
  {}

  Iso_rectangle_2(const RIso_rectangle_2& r)
    : RIso_rectangle_2(r) {}

  Iso_rectangle_2(const Point_2 &p, const Point_2 &q)
    : RIso_rectangle_2(typename R::Construct_iso_rectangle_2()(p,q).rep()) {}

  Iso_rectangle_2(const Point_2 &left, const Point_2 &right,
                  const Point_2 &bottom, const Point_2 &top)
    : RIso_rectangle_2(typename R::Construct_iso_rectangle_2()(left, right, bottom, top).rep()) {}

  Iso_rectangle_2(const RT& min_hx, const RT& min_hy, 
                  const RT& max_hx, const RT& max_hy)
    : RIso_rectangle_2(typename R::Construct_iso_rectangle_2()(min_hx, min_hy, max_hx, max_hy).rep()) {}

  Iso_rectangle_2(const RT& min_hx, const RT& min_hy, 
                  const RT& max_hx, const RT& max_hy, const RT& hw)
    : RIso_rectangle_2(typename R::Construct_iso_rectangle_2()(min_hx, min_hy, max_hx, max_hy, hw).rep()) {}


  const Point_2 & min() const
  {
    return R().construct_min_vertex_2_object()(*this);
  }

  const Point_2 & max() const
  {
    return R().construct_max_vertex_2_object()(*this);
  }

  bool
  operator==(const Iso_rectangle_2 &i) const
  {
    return R().equal_2_object()(*this, i);
  }

  bool
  operator!=(const Iso_rectangle_2 &i) const
  {
    return ! (*this == i);
  }

  Point_2 
  vertex(int i) const
  {
    return R().construct_vertex_2_object()(*this,i);
  }

  Point_2 
  operator[](int i) const
  {
    return R().construct_vertex_2_object()(*this,i);
  }

  Cartesian_coordinate_type xmin() const
  {
    return R().compute_xmin_2_object()(*this);
  }

  Cartesian_coordinate_type xmax() const
  {
    return R().compute_xmax_2_object()(*this);
  }

  Cartesian_coordinate_type ymin() const
  {
    return R().compute_ymin_2_object()(*this);
  }

  Cartesian_coordinate_type ymax() const
  {
    return R().compute_ymax_2_object()(*this);
  }

  Cartesian_coordinate_type
  min_coord(int i) const
  {
    CGAL_kernel_precondition( i == 0 || i == 1 );
    if (i == 0)
      return xmin();
    else
      return ymin();
  }

  Cartesian_coordinate_type
  max_coord(int i) const
  {
    CGAL_kernel_precondition( i == 0 || i == 1 );
    if (i == 0)
      return xmax();
    else
      return ymax();
  }

  FT
  area() const
  {
    return R().compute_area_2_object()(*this);
  }


  bool
  has_on_boundary(const Point_2 &p) const
  {
    return R().has_on_boundary_2_object()(*this,p);
  }


  bool
  has_on_bounded_side(const Point_2 &p) const
  {
    return R().has_on_bounded_side_2_object()(*this,p);
  }


  bool
  has_on_unbounded_side(const Point_2 &p) const
  {
    return R().has_on_unbounded_side_2_object()(*this,p);
  }

  Bounded_side
  bounded_side(const Point_2 &p) const
  {
    return R().bounded_side_2_object()(*this,p);
  }


  bool
  is_degenerate() const
  {
    return R().is_degenerate_2_object()(*this);
  }

  Bbox_2
  bbox() const
  { 
    return R().construct_bbox_2_object()(*this);
  }

  Iso_rectangle_2 transform(const Aff_transformation_2 &t) const
  {
    //return Iso_rectangle_2(t.transform(min()), t.transform(max()));
    return rep().transform(t);
  }
};

#ifndef CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Iso_rectangle_2<R> &r)
{
  return  os << r.rep();
}
#endif // CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLE_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLE_2
template < class R >
std::istream &
operator>>(std::istream &is, Iso_rectangle_2<R> &r)
{
  return is >> r.rep();
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLE_2

CGAL_END_NAMESPACE

#endif // CGAL_ISO_RECTANGLE_2_H
