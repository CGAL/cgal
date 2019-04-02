// Copyright (c) 2019  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_PAIR_PARTITION_TRAITS_2_H
#define CGAL_PAIR_PARTITION_TRAITS_2_H

#include <CGAL/license/Partition_2.h>

#include <CGAL/polygon_function_objects.h>
#include <CGAL/Polygon_2.h>
#include <list>

namespace CGAL {


  template <typename Pair, typename Functor>
  struct Pair_functor {
    Functor fct;

    Pair_functor()
      : fct()
    {}
    
    Pair_functor(const Functor fct)
      : fct(fct)
    {}

    typename Functor::result_type
    operator()(const Pair& p, const Pair& q)const
    {
      return fct(p.first,q.first);
    }
    
    typename Functor::result_type
    operator()(const Pair& p, const Pair& q, const Pair& r)const
    {
      return fct(p.first,q.first,r.first);
    }
  };

  template <typename Pair, typename K>
  struct Pair_compare_x_at_y_2 {
    
    typename K::Compare_x_at_y_2 fct;
    
    Pair_compare_x_at_y_2(typename K::Compare_x_at_y_2 fct)
      : fct(fct)
    {}
    
    typename K::Compare_x_at_y_2::result_type
    operator()(const Pair& p, const typename K::Line_2& line) const
    {
      return fct(p.first,line);
    }
  };

  template <typename Pair, typename K>
  struct Pair_collinear_are_ordered_along_line_2 {

    Pair_collinear_are_ordered_along_line_2()
    {}
    
    typename K::Collinear_are_ordered_along_line_2 fct;
    
    Pair_collinear_are_ordered_along_line_2(typename K::Collinear_are_ordered_along_line_2 fct)
      : fct(fct)
    {}
    
    typename K::Collinear_are_ordered_along_line_2::result_type
    operator()(const Pair& p, const Pair& q, const typename K::Point_2& r) const
    {
      return fct(p.first, q.first, r);
    }
    
    typename K::Collinear_are_ordered_along_line_2::result_type
    operator()(const Pair& p, const Pair& q, const Pair& r) const
    {
      return fct(p.first, q.first, r.first);
    }
  };

  
template <class Kernel_>
class Pair_partition_traits_2
{
private:
  typedef Kernel_                                     Kernel;
  typedef Pair_partition_traits_2<Kernel_>    Self;
  typedef Surface_mesh<typename Kernel::Point_2>      Surface_mesh;
  
public:

  Pair_partition_traits_2()
  {}

  
  typedef typename boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;

  typedef typename Kernel::FT FT;
  typedef typename boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
  typedef std::pair<typename Kernel::Point_2, vertex_descriptor> Point_2;
  typedef ::std::list<Point_2>                        Container;
  typedef CGAL::Polygon_2<Self, Container>          Polygon_2;

  
  typedef Pair_functor<Point_2, typename Kernel::Equal_2>                Equal_2;

  typedef Pair_functor<Point_2, typename Kernel::Less_yx_2>                  Less_yx_2;
  typedef Pair_functor<Point_2, typename Kernel::Less_xy_2>                  Less_xy_2;
  typedef Pair_functor<Point_2, typename Kernel::Left_turn_2>                 Left_turn_2;
  typedef Pair_functor<Point_2, typename Kernel::Orientation_2>              Orientation_2;
  typedef Pair_functor<Point_2, typename Kernel::Compare_y_2>                Compare_y_2;
  typedef Pair_functor<Point_2, typename Kernel::Compare_x_2>                Compare_x_2;
  typedef CGAL::Is_convex_2<Self>                     Is_convex_2;
  typedef CGAL::Is_y_monotone_2<Self>                 Is_y_monotone_2;

  // needed by Indirect_edge_compare, used in y_monotone and greene_approx
  typedef typename Kernel::Line_2                     Line_2;
  typedef Pair_functor<Point_2, typename Kernel::Construct_line_2>           Construct_line_2;
  typedef Pair_compare_x_at_y_2<Point_2,Kernel>            Compare_x_at_y_2;
  typedef typename Kernel::Is_horizontal_2            Is_horizontal_2;
  
#if 1
  // needed by visibility graph and thus by optimal convex
  typedef typename Kernel::Ray_2                      Ray_2; 
  typedef Pair_collinear_are_ordered_along_line_2<Point_2,Kernel>  Collinear_are_ordered_along_line_2;
  typedef Pair_functor<Point_2,typename Kernel::Are_strictly_ordered_along_line_2>
  Are_strictly_ordered_along_line_2;
  typedef typename Kernel::Intersect_2                Intersect_2;
  typedef typename Kernel::Assign_2                   Assign_2;
  typedef typename Kernel::Object_2                   Object_2;

  // needed by approx_convex (for constrained triangulation)
  // and optimal convex (for vis. graph)
#endif  
  typedef typename Kernel::Segment_2                  Segment_2;
#if 1 
  // needed by optimal convex (for vis. graph)
  typedef Pair_functor<Point_2,typename Kernel::Construct_segment_2>        Construct_segment_2;
  typedef Pair_functor<Point_2,typename Kernel::Construct_ray_2>            Construct_ray_2;
#endif
  Equal_2
  equal_2_object() const
  { return Equal_2(Kernel().equal_2_object()); }

  Orientation_2
  orientation_2_object() const
  { return Orientation_2(Kernel().orientation_2_object()); }



    Less_yx_2
    less_yx_2_object() const
    { return Less_yx_2(Kernel().less_yx_2_object()); }

    Less_xy_2
    less_xy_2_object() const
    { return Less_xy_2(Kernel().less_xy_2_object()); }

    Left_turn_2
    left_turn_2_object() const
    { return Left_turn_2(Kernel().left_turn_2_object()); }

    Compare_y_2
    compare_y_2_object() const
    {  return Compare_y_2(Kernel().compare_y_2_object()); }

    Compare_x_2
    compare_x_2_object() const
    {  return Compare_x_2(Kernel().compare_x_2_object()); }

  
  
  Construct_line_2
  construct_line_2_object() const
  {  return Construct_line_2(Kernel().construct_line_2_object()); }

  Compare_x_at_y_2
  compare_x_at_y_2_object() const
  { return Compare_x_at_y_2(Kernel().compare_x_at_y_2_object()); }


#if 1

  Construct_segment_2
  construct_segment_2_object() const
  { return Construct_segment_2(); }

  Construct_ray_2
  construct_ray_2_object() const
  { return Construct_ray_2(); }

  Collinear_are_ordered_along_line_2
  collinear_are_ordered_along_line_2_object() const
  { return Collinear_are_ordered_along_line_2(); }

  Are_strictly_ordered_along_line_2
  are_strictly_ordered_along_line_2_object() const
  { return Are_strictly_ordered_along_line_2(); }

#endif
  
  Is_horizontal_2
  is_horizontal_2_object() const
  {  return Is_horizontal_2(); }

  Is_convex_2
  is_convex_2_object(const Self& traits) const
  {  return Is_convex_2(traits); }

  Is_y_monotone_2
  is_y_monotone_2_object(const Self& traits) const
  {  return Is_y_monotone_2(traits); }

#if 1
  Intersect_2
  intersect_2_object() const
  {  return Intersect_2(); }

  Assign_2
  assign_2_object() const
  {  return Assign_2(); }
#endif  
};

}

#endif // CGAL_PAIR_PARTITION_TRAITS_2_H
