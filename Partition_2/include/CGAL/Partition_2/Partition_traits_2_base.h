// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_PARTITION_TRAITS_2_BASE_H
#define CGAL_PARTITION_TRAITS_2_BASE_H

#include <CGAL/license/Partition_2.h>


#include <CGAL/Polygon_2.h>
#include <list>

namespace CGAL {

template <class Kernel_>
class Partition_traits_2_base
{
  private:
    typedef Kernel_                                 Kernel;
  public:
    typedef typename Kernel::Point_2                Point_2;
    typedef ::std::list<Point_2>                    Container;
    typedef CGAL::Polygon_2<Kernel, Container>      Polygon_2;
    typedef typename Kernel::Equal_2                Equal_2;
    typedef typename Kernel::Less_yx_2              Less_yx_2;
    typedef typename Kernel::Less_xy_2              Less_xy_2;
    typedef typename Kernel::Left_turn_2             Left_turn_2;
    typedef typename Kernel::Orientation_2          Orientation_2;
    typedef typename Kernel::Compare_y_2            Compare_y_2;
    typedef typename Kernel::Compare_x_2            Compare_x_2;

    Equal_2
    equal_2_object() const
    { return Equal_2(); }

    Less_yx_2
    less_yx_2_object() const
    { return Less_yx_2(); }

    Less_xy_2
    less_xy_2_object() const
    { return Less_xy_2(); }

    Left_turn_2
    left_turn_2_object() const
    { return Left_turn_2(); }

    Orientation_2
    orientation_2_object() const
    { return Orientation_2(); }

    Compare_y_2
    compare_y_2_object() const
    {  return Compare_y_2(); }

    Compare_x_2
    compare_x_2_object() const
    {  return Compare_x_2(); }

};

}

#endif // CGAL_PARTITION_TRAITS_2_BASE_H
