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

#ifndef CGAL_PARTITION_TRAITS_ADAPTER_2_H
#define CGAL_PARTITION_TRAITS_ADAPTER_2_H

#include <CGAL/license/Partition_2.h>

#include <boost/call_traits.hpp>

#include <CGAL/property_map.h>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/Polygon_2.h>
#include <list>

namespace CGAL {

  
template <class Base_traits,class PointPropertyMap>
class Partition_traits_adapter_2 : public Base_traits
{
private:
  typedef typename Base_traits::Kernel                                   Kernel;
  typedef Partition_traits_adapter_2<Base_traits,PointPropertyMap>       Self;
  
  PointPropertyMap ppmap;
public:

  Partition_traits_adapter_2(Base_traits base=Base_traits())
    : Base_traits(base)
  {}
  
  Partition_traits_adapter_2(const PointPropertyMap& ppmap,Base_traits base=Base_traits())
  : Base_traits(base),ppmap(ppmap)
  {}
  
  typedef typename Kernel::FT FT;
  typedef typename boost::property_traits<PointPropertyMap>::key_type Point_2;
  typedef typename boost::call_traits<Point_2>::param_type Arg_type;

  typedef ::std::list<Point_2>                      Container;
  typedef CGAL::Polygon_2<Self, Container>          Polygon_2;

  
  template <typename BaseFct>
  struct Pmap_fct : public BaseFct {
    Pmap_fct(const PointPropertyMap& ppmap, const BaseFct& base)
      : BaseFct(base),ppmap(ppmap)
    {}
    
    const PointPropertyMap& ppmap;

    typename BaseFct::result_type operator()(Arg_type p, Arg_type q) const {
      return static_cast<const typename BaseFct*>(this)->operator()(get(ppmap,p),get(ppmap,q));
    }
    
    typename BaseFct::result_type operator()(Arg_type p, Arg_type q, Arg_type r) const {
      return static_cast<const typename BaseFct*>(this)->operator()(get(ppmap,p),get(ppmap,q),get(ppmap,r));
    }
  };

    template <typename PointPropertyMap, typename K>
  struct Pmap_compare_x_at_y_2 {
    
    PointPropertyMap ppmap;
    typename K::Compare_x_at_y_2 fct;
    
    Pmap_compare_x_at_y_2(const PointPropertyMap& ppmap,typename K::Compare_x_at_y_2 fct)
      : ppmap(ppmap), fct(fct)
    {}
    
    typename K::Compare_x_at_y_2::result_type
    operator()(Arg_type p, const typename K::Line_2& line) const
    {
      return fct(get(ppmap,p),line);
    }
  };


  
  template <typename PointPropertyMap, typename K>
  struct Pmap_collinear_are_ordered_along_line_2 {

    PointPropertyMap ppmap;
    typename K::Collinear_are_ordered_along_line_2 fct;
    
    Pmap_collinear_are_ordered_along_line_2(const PointPropertyMap& ppmap,typename K::Collinear_are_ordered_along_line_2 fct)
      : ppmap(ppmap), fct(fct) 
    {}
    
    
    
    typename K::Collinear_are_ordered_along_line_2::result_type
    operator()(Arg_type p, Arg_type q, const typename K::Point_2& r) const
    {
      return fct(p.first, q.first, r);
    }
    
    typename K::Collinear_are_ordered_along_line_2::result_type
    operator()(Arg_type p, Arg_type q, Arg_type r) const
    {
      return fct(get(ppmap,p), get(ppmap,q), get(ppmap,r));
    }
  };

  
  typedef Pmap_fct<typename Kernel::Equal_2>                    Equal_2;
  typedef Pmap_fct<typename Kernel::Less_yx_2>                  Less_yx_2;
  typedef Pmap_fct<typename Kernel::Less_xy_2>                  Less_xy_2;
  typedef Pmap_fct<typename Kernel::Left_turn_2>                Left_turn_2;
  typedef Pmap_fct<typename Kernel::Orientation_2>              Orientation_2;
  typedef Pmap_fct<typename Kernel::Compare_y_2>                Compare_y_2;
  typedef Pmap_fct<typename Kernel::Compare_x_2>                Compare_x_2;
  typedef CGAL::Is_convex_2<Self>                                            Is_convex_2;
  typedef CGAL::Is_y_monotone_2<Self>                                        Is_y_monotone_2;

  // needed by Indirect_edge_compare, used in y_monotone and greene_approx
  typedef typename Kernel::Line_2                                            Line_2;
  typedef Pmap_fct<typename Kernel::Construct_line_2>           Construct_line_2;
  typedef Pmap_compare_x_at_y_2<PointPropertyMap,Kernel>                              Compare_x_at_y_2;
  typedef typename Kernel::Is_horizontal_2                                   Is_horizontal_2;
  
  // needed by visibility graph and thus by optimal convex
  typedef Pmap_collinear_are_ordered_along_line_2<PointPropertyMap,Kernel>   Collinear_are_ordered_along_line_2;
  typedef Pmap_fct<typename Kernel::Are_strictly_ordered_along_line_2>
                                                                             Are_strictly_ordered_along_line_2;

  // needed by approx_convex (for constrained triangulation)
  // and optimal convex (for vis. graph)

  typedef typename Kernel::Segment_2                                        Segment_2;

  // needed by optimal convex (for vis. graph)
  typedef Pmap_fct<typename Kernel::Construct_segment_2>        Construct_segment_2;

  Equal_2
  equal_2_object() const
  { return Equal_2(ppmap,static_cast<const Base_traits*>(this)->equal_2_object()); }

  Orientation_2
  orientation_2_object() const
  { return Orientation_2(ppmap,static_cast<const Base_traits*>(this)->orientation_2_object()); }

  Less_yx_2
  less_yx_2_object() const
  { return Less_yx_2(ppmap,static_cast<const Base_traits*>(this)->less_yx_2_object()); }
  
  Less_xy_2
  less_xy_2_object() const
  { return Less_xy_2(ppmap,static_cast<const Base_traits*>(this)->less_xy_2_object()); }
  
  Left_turn_2
  left_turn_2_object() const
  { return Left_turn_2(ppmap,static_cast<const Base_traits*>(this)->left_turn_2_object()); }
  
  Compare_y_2
  compare_y_2_object() const
  {  return Compare_y_2(ppmap,static_cast<const Base_traits*>(this)->compare_y_2_object()); }
  
  Compare_x_2
  compare_x_2_object() const
  {  return Compare_x_2(ppmap,static_cast<const Base_traits*>(this)->compare_x_2_object()); }
  
  
  Construct_line_2
  construct_line_2_object() const
  {  return Construct_line_2(ppmap,static_cast<const Base_traits*>(this)->construct_line_2_object()); }
  
  Compare_x_at_y_2
  compare_x_at_y_2_object() const
  { return Compare_x_at_y_2(ppmap,static_cast<const Base_traits*>(this)->compare_x_at_y_2_object()); }


  Construct_segment_2
  construct_segment_2_object() const
  { return Construct_segment_2(); }

  Collinear_are_ordered_along_line_2
  collinear_are_ordered_along_line_2_object() const
  { return Collinear_are_ordered_along_line_2(ppmap,static_cast<const Base_traits*>(this)->collinear_are_ordered_along_line_2_object()); }

  Are_strictly_ordered_along_line_2
  are_strictly_ordered_along_line_2_object() const
  { return Are_strictly_ordered_along_line_2(ppmap,static_cast<const Base_traits*>(this)->are_strictly_ordered_along_line_2_object()); }

  
  Is_horizontal_2
  is_horizontal_2_object() const
  {  return Is_horizontal_2(); }

  Is_convex_2
  is_convex_2_object(const Self& traits) const
  {  return Is_convex_2(traits); }

  Is_y_monotone_2
  is_y_monotone_2_object(const Self& traits) const
  {  return Is_y_monotone_2(traits); }



};

}

#endif // CGAL_PARTITION_TRAITS_ADAPTER_2_H
