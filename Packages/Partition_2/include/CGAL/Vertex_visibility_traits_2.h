// Copyright (c) 2000  Max-Planck-Institute Saarbrucken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef VERTEX_VISIBILITY_TRAITS_2_H
#define VERTEX_VISIBILITY_TRAITS_2_H

namespace CGAL {

template <class Kernel_>
class Vertex_visibility_traits_2  
{
  public:
    typedef Kernel_                               Kernel;
    typedef typename Kernel::Point_2              Point_2; 
    typedef typename Kernel::Segment_2            Segment_2; 
    typedef typename Kernel::Ray_2                Ray_2; 
    typedef typename Kernel::Construct_segment_2  Construct_segment_2;
    typedef typename Kernel::Construct_ray_2      Construct_ray_2;
    typedef typename Kernel::Less_yx_2            Less_yx_2;
    typedef typename Kernel::Less_xy_2            Less_xy_2;
    typedef typename Kernel::Compare_x_2          Compare_x_2;
    typedef typename Kernel::Compare_y_2          Compare_y_2;
    typedef typename Kernel::Left_turn_2           Left_turn_2;
    typedef typename Kernel::Orientation_2        Orientation_2;
    typedef typename Kernel::Collinear_are_ordered_along_line_2
                                          Collinear_are_ordered_along_line_2;
    typedef typename Kernel::Are_strictly_ordered_along_line_2
                                          Are_strictly_ordered_along_line_2;

    Compare_x_2
    compare_x_2_object() const
    {  return Compare_x_2(); }

    Compare_y_2
    compare_y_2_object() const
    {  return Compare_y_2(); }

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

    Collinear_are_ordered_along_line_2
    collinear_are_ordered_along_line_2_object() const
    { return Collinear_are_ordered_along_line_2(); }

    Are_strictly_ordered_along_line_2
    are_strictly_ordered_along_line_2_object() const
    { return Are_strictly_ordered_along_line_2(); }

    Construct_segment_2
    construct_segment_2_object() const
    { return Construct_segment_2(); }

    Construct_ray_2
    construct_ray_2_object() const
    { return Construct_ray_2(); }

};

}

#endif // VERTEX_VISIBILITY_TRAITS_2_H
