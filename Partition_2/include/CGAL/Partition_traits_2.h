// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_PARTITION_TRAITS_2_H
#define CGAL_PARTITION_TRAITS_2_H

#include <CGAL/license/Partition_2.h>


#include <CGAL/Polygon_2.h>
#include <CGAL/Partition_2/Partition_traits_2_base.h>
#include <CGAL/polygon_function_objects.h>
#include <list>

namespace CGAL {

template <class Kernel_>
class Partition_traits_2  : public Partition_traits_2_base<Kernel_>
{
  private:
    typedef Kernel_                                     Kernel;
    typedef Partition_traits_2<Kernel_>                 Self;
  
  public:
    typedef typename Kernel::Point_2                    Point_2;
    typedef ::std::list<Point_2>                        Container;
    typedef CGAL::Polygon_2<Kernel, Container>          Polygon_2;
    typedef typename Kernel::Less_yx_2                  Less_yx_2;
    typedef typename Kernel::Less_xy_2                  Less_xy_2;
    typedef typename Kernel::Left_turn_2                Left_turn_2;
    typedef typename Kernel::Orientation_2              Orientation_2;
    typedef typename Kernel::Compare_y_2                Compare_y_2;
    typedef typename Kernel::Compare_x_2                Compare_x_2;
    typedef CGAL::Is_convex_2<Self>                     Is_convex_2;
    typedef CGAL::Is_y_monotone_2<Self>                 Is_y_monotone_2;


    // needed by visibility graph and thus by optimal convex
    typedef typename Kernel::Collinear_are_ordered_along_line_2
                                            Collinear_are_ordered_along_line_2;
    typedef typename Kernel::Are_strictly_ordered_along_line_2
                                            Are_strictly_ordered_along_line_2;

    Collinear_are_ordered_along_line_2
    collinear_are_ordered_along_line_2_object() const
    { return Collinear_are_ordered_along_line_2(); }

    Are_strictly_ordered_along_line_2
    are_strictly_ordered_along_line_2_object() const
    { return Are_strictly_ordered_along_line_2(); }

    Is_convex_2
    is_convex_2_object(const Self& traits) const
    {  return Is_convex_2(traits); }

    Is_y_monotone_2
    is_y_monotone_2_object(const Self& traits) const
    {  return Is_y_monotone_2(traits); }

};

}

#endif // CGAL_PARTITION_TRAITS_2_H
