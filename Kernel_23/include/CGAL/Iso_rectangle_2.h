// Copyright (c) 1999  
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_ISO_RECTANGLE_2_H
#define CGAL_ISO_RECTANGLE_2_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Dimension.h>
#include <CGAL/result_of.h>

namespace CGAL {

template <class R_>
class Iso_rectangle_2 : public R_::Kernel_base::Iso_rectangle_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::FT                    FT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Aff_transformation_2  Aff_transformation_2;

  typedef Iso_rectangle_2                    Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Iso_rectangle_2>::value));

public:

  typedef Dimension_tag<2>  Ambient_dimension;
  typedef Dimension_tag<2>  Feature_dimension;

  typedef typename R_::Kernel_base::Iso_rectangle_2  Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_   R;

  Iso_rectangle_2() {}

  Iso_rectangle_2(const Rep& r)
    : Rep(r) {}

  Iso_rectangle_2(const Point_2 &p, const Point_2 &q, int)
    : Rep(typename R::Construct_iso_rectangle_2()(Return_base_tag(), p, q, 0)) {}

  Iso_rectangle_2(const Point_2 &p, const Point_2 &q)
    : Rep(typename R::Construct_iso_rectangle_2()(Return_base_tag(), p, q)) {}

  Iso_rectangle_2(const Point_2 &left, const Point_2 &right,
                  const Point_2 &bottom, const Point_2 &top)
    : Rep(typename R::Construct_iso_rectangle_2()(Return_base_tag(), left, right, bottom, top)) {}

  Iso_rectangle_2(const RT& min_hx, const RT& min_hy,
                  const RT& max_hx, const RT& max_hy)
    : Rep(typename R::Construct_iso_rectangle_2()(Return_base_tag(), min_hx, min_hy, max_hx, max_hy)) {}

  Iso_rectangle_2(const RT& min_hx, const RT& min_hy,
                  const RT& max_hx, const RT& max_hy, const RT& hw)
    : Rep(typename R::Construct_iso_rectangle_2()(Return_base_tag(), min_hx, min_hy, max_hx, max_hy, hw)) {}

  Iso_rectangle_2(const Bbox_2& bbox)
    : Rep(typename R::Construct_iso_rectangle_2()(Return_base_tag(), bbox.xmin(), bbox.ymin(), bbox.xmax(), bbox.ymax())) {}

  typename cpp11::result_of<typename R::Construct_min_vertex_2( Iso_rectangle_2 )>::type
  min BOOST_PREVENT_MACRO_SUBSTITUTION () const
  {
    return R().construct_min_vertex_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Construct_max_vertex_2( Iso_rectangle_2 )>::type
  max BOOST_PREVENT_MACRO_SUBSTITUTION () const
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


  typename cpp11::result_of<typename R::Construct_vertex_2( Iso_rectangle_2, int )>::type
  vertex(int i) const
  {
    return R().construct_vertex_2_object()(*this,i);
  }

  typename cpp11::result_of<typename R::Construct_vertex_2( Iso_rectangle_2, int )>::type
  operator[](int i) const
  {
    return R().construct_vertex_2_object()(*this,i);
  }

  typename cpp11::result_of<typename R::Compute_xmin_2( Iso_rectangle_2 )>::type
  xmin() const
  {
    return R().compute_xmin_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_xmax_2( Iso_rectangle_2 )>::type
  xmax() const
  {
    return R().compute_xmax_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_ymin_2( Iso_rectangle_2 )>::type
  ymin() const
  {
    return R().compute_ymin_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_ymax_2( Iso_rectangle_2 )>::type
  ymax() const
  {
    return R().compute_ymax_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_xmin_2( Iso_rectangle_2 )>::type
  min_coord(int i) const
  {
    CGAL_kernel_precondition( i == 0 || i == 1 );
    if (i == 0)
      return xmin();
    else
      return ymin();
  }

  typename cpp11::result_of<typename R::Compute_xmin_2( Iso_rectangle_2 )>::type
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
    // FIXME : We need a precondition like this!!!
    // CGAL_kernel_precondition(t.is_axis_preserving());
    return Iso_rectangle_2(t.transform(min  BOOST_PREVENT_MACRO_SUBSTITUTION ()), 
			   t.transform(max  BOOST_PREVENT_MACRO_SUBSTITUTION ()));
  }
};


template < class R >
std::ostream &
operator<<(std::ostream &os, const Iso_rectangle_2<R> &r)
{
  switch(get_mode(os)) {
  case IO::ASCII :
    return os << (r.min)() << ' ' << (r.max)();
  case IO::BINARY :
    return os << (r.min)() << (r.max)();
  default:
    return os << "Iso_rectangle_2(" << (r.min)() << ", " << (r.max)() << ")";
  }
}

template < class R >
std::istream &
operator>>(std::istream &is, Iso_rectangle_2<R> &r)
{
  typename R::Point_2 p, q;

  is >> p >> q;

  if (is)
    r = Iso_rectangle_2<R>(p, q);
  return is;
}

} //namespace CGAL

#endif // CGAL_ISO_RECTANGLE_2_H
