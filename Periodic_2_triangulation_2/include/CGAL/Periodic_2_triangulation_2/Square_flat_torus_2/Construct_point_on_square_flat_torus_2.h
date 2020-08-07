// Copyright (c) 1997-2013, 2017 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>,
//                 Mael Rouxel-Labbé

#ifndef CGAL_P2T2_CONSTRUCT_POINT_ON_SQUARE_FLAT_TORUS_2_H
#define CGAL_P2T2_CONSTRUCT_POINT_ON_SQUARE_FLAT_TORUS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_offset_2.h>

namespace CGAL {
namespace P2T2 {
namespace internal {

// the domain is necessarily an iso rectangle
template <typename K_,
          typename Offset_ = CGAL::Periodic_2_offset_2,
          typename Construct_point_2_base_ = typename K_::Construct_point_2>
class Construct_point_on_square_flat_torus_2
  : public Construct_point_2_base_
{
  typedef Construct_point_2_base_            Base;
  typedef K_                                 Kernel;

  typedef typename Kernel::Point_2           Point;
  typedef Offset_                            Offset;
  typedef typename Kernel::Iso_rectangle_2   Iso_rectangle_2;

public:
  Construct_point_on_square_flat_torus_2(const Iso_rectangle_2* dom, const Base& cp)
    : Base(cp), _dom(dom)
  { }

  using Base::operator();

  Point operator()(const Point& p, const Offset& o) const
  {
    return operator()(p.x() + (_dom->xmax() - _dom->xmin()) * o.x(),
                      p.y() + (_dom->ymax() - _dom->ymin()) * o.y());
  }

private:
  const Iso_rectangle_2* _dom;
};

} // namespace internal
} // namespace P2T2
} // namespace CGAL

#endif // CGAL_P2T2_CONSTRUCT_POINT_ON_SQUARE_FLAT_TORUS_2_H
