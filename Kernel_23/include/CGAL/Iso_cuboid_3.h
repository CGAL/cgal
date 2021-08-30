// Copyright (c) 1999
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
// Author(s)     : Stefan Schirra

#ifndef CGAL_ISO_CUBOID_3_H
#define CGAL_ISO_CUBOID_3_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R_>
class Iso_cuboid_3 : public R_::Kernel_base::Iso_cuboid_3
{
  typedef typename R_::RT                 RT;
  typedef typename R_::Point_3            Point_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;

  typedef Iso_cuboid_3                    Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Iso_cuboid_3>::value));

public:

  typedef Dimension_tag<3>  Ambient_dimension;
  typedef Dimension_tag<3>  Feature_dimension;

  typedef typename R_::Kernel_base::Iso_cuboid_3  Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                    R;

  Iso_cuboid_3() {}

  Iso_cuboid_3(const Rep&  r)
      : Rep(r) {}

  Iso_cuboid_3(const Point_3& p, const Point_3& q)
   : Rep(typename R::Construct_iso_cuboid_3()(Return_base_tag(), p,q)) {}

  Iso_cuboid_3(const Point_3& p, const Point_3& q, int)
   : Rep(typename R::Construct_iso_cuboid_3()(Return_base_tag(), p, q, 0)) {}

  Iso_cuboid_3(const Point_3 &left,   const Point_3 &right,
               const Point_3 &bottom, const Point_3 &top,
               const Point_3 &far_,   const Point_3 &close)
   : Rep(typename R::Construct_iso_cuboid_3()(Return_base_tag(), left, right, bottom,
                                                    top, far_, close)) {}

  Iso_cuboid_3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz,
               const RT& hw)
   : Rep(typename R::Construct_iso_cuboid_3()(Return_base_tag(), min_hx, min_hy, min_hz,
                                     max_hx, max_hy, max_hz, hw)) {}

  Iso_cuboid_3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz)
   : Rep(typename R::Construct_iso_cuboid_3()(Return_base_tag(), min_hx, min_hy, min_hz,
                                             max_hx, max_hy, max_hz)) {}

  Iso_cuboid_3(const Bbox_3& bbox)
   : Rep(typename R::Construct_iso_cuboid_3()(Return_base_tag(), bbox.xmin(), bbox.ymin(), bbox.zmin(),
                                                                 bbox.xmax(), bbox.ymax(), bbox.zmax())) {}

  decltype(auto)
  min BOOST_PREVENT_MACRO_SUBSTITUTION () const
  {
    return R().construct_min_vertex_3_object()(*this);
  }

  decltype(auto)
  max BOOST_PREVENT_MACRO_SUBSTITUTION () const
  {
    return R().construct_max_vertex_3_object()(*this);
  }

  decltype(auto)
  vertex(int i) const
  {
    return R().construct_vertex_3_object()(*this,i);
  }

  decltype(auto)
  operator[](int i) const
  {
    return R().construct_vertex_3_object()(*this,i);
  }

  decltype(auto)
  xmin() const
  {
    return R().compute_xmin_3_object()(*this);
  }

  decltype(auto)
  xmax() const
  {
    return R().compute_xmax_3_object()(*this);
  }

  decltype(auto)
  ymin() const
  {
    return R().compute_ymin_3_object()(*this);
  }

  decltype(auto)
  ymax() const
  {
    return R().compute_ymax_3_object()(*this);
  }

  decltype(auto)
  zmin() const
  {
    return R().compute_zmin_3_object()(*this);
  }

  decltype(auto)
  zmax() const
  {
    return R().compute_zmax_3_object()(*this);
  }

  decltype(auto)
  min_coord(int i) const
  {
    CGAL_kernel_precondition( i == 0 || i == 1 || i == 2 );
    if (i == 0)
       return xmin();
    else if (i == 1)
       return ymin();
    else
       return zmin();
  }

  decltype(auto)
  max_coord(int i) const
  {
    CGAL_kernel_precondition( i == 0 || i == 1 || i == 2 );
    if (i == 0)
       return xmax();
    else if (i == 1)
       return ymax();
    else
       return zmax();
  }

  bool
  has_on_bounded_side(const Point_3 &p) const
  {
    return R().has_on_bounded_side_3_object()(*this,p);
  }

  bool
  has_on_unbounded_side(const Point_3 &p) const
  {
    return R().has_on_unbounded_side_3_object()(*this,p);
  }

  bool
  has_on_boundary(const Point_3 &p) const
  {
    return R().has_on_boundary_3_object()(*this,p);
  }

  bool
  has_on(const Point_3 &p) const
  {
    return has_on_boundary(p);
  }

  Bounded_side
  bounded_side(const Point_3 &p) const
  {
    return R().bounded_side_3_object()(*this,p);
  }

  bool
  is_degenerate() const
  {
    return R().is_degenerate_3_object()(*this);
  }

  decltype(auto)
  volume() const
  {
    return R().compute_volume_3_object()(*this);
  }

  Bbox_3
  bbox() const
  {
    return R().construct_bbox_3_object()(*this);
  }

  Iso_cuboid_3
  transform(const Aff_transformation_3 &t) const
  {
    return Iso_cuboid_3(t.transform((this->min)()), t.transform((this->max)()));
  }

};


template < class R >
std::ostream &
operator<<(std::ostream& os, const Iso_cuboid_3<R>& r)
{
  switch(IO::get_mode(os)) {
  case IO::ASCII :
    return os << (r.min)() << ' ' << (r.max)();
  case IO::BINARY :
    return os << (r.min)() << (r.max)();
  default:
    return os << "Iso_cuboid_3(" << (r.min)() << ", " << (r.max)() << ")";
  }
}

template < class R >
std::istream &
operator>>(std::istream& is, Iso_cuboid_3<R>& r)
{
  typename R::Point_3 p, q;
  is >> p >> q;
  if (is)
      r = Iso_cuboid_3<R>(p, q);
  return is;
}

} //namespace CGAL

#endif // CGAL_ISO_CUBOID_3_H
