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

#ifndef CGAL_PERIODIC_2_CONSTRUCT_POINT_2_H
#define CGAL_PERIODIC_2_CONSTRUCT_POINT_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

namespace CGAL {

template < typename K_, typename Construct_point_2_base_>
class Periodic_2_construct_point_2
  : public Construct_point_2_base_
{
  typedef Construct_point_2_base_            Base;
  typedef K_                                 Kernel;

  typedef typename Kernel::Point_2           Point;
  typedef typename Kernel::Offset            Offset;
  typedef typename Kernel::Iso_rectangle_2   Iso_rectangle_2;

public:
  Periodic_2_construct_point_2(const Iso_rectangle_2* dom,
                               const Base& cp)
    : Base(cp), _dom(dom)
  { }

  using Base::operator();

  Point operator() ( const Point& p, const Offset& o ) const
  {
    return operator()(p.x() + (_dom->xmax() - _dom->xmin()) * o.x(),
                      p.y() + (_dom->ymax() - _dom->ymin()) * o.y());
  }

private:
  const Iso_rectangle_2* _dom;
};

} // namespace CGAL

#endif // CGAL_PERIODIC_2_CONSTRUCT_POINT_2_H
