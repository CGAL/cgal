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

#ifndef CGAL_PARTITION_TRAITS_2_H
#define CGAL_PARTITION_TRAITS_2_H

#include <CGAL/license/Partition_2.h>


#include <boost/call_traits.hpp>

#include <CGAL/property_map.h>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/Polygon_2.h>
#include <list>

namespace CGAL {

  template <class Base_traits, class PointPropertyMap = Identity_property_map<typename Base_traits::Point_2> >
class Partition_traits_2;

  template <typename BT, typename PM>
  struct Polygon_traits_getter{
    typedef Partition_traits_2<BT,PM> type;
  };

  template <typename BT>
  struct Polygon_traits_getter<BT,Identity_property_map<typename BT::Point_2> > {
    typedef BT type;
  };

template <class Base_traits, class PointPropertyMap>
class Partition_traits_2 : public Base_traits
{
private:
  typedef Base_traits                                              Kernel;
  typedef Partition_traits_2<Base_traits,PointPropertyMap> Self;

  PointPropertyMap ppmap;

public:

  Partition_traits_2(const Base_traits& base=Base_traits())
    : Base_traits(base)
  {}

  Partition_traits_2(const PointPropertyMap& ppmap,
                             const Base_traits& base=Base_traits())
  : Base_traits(base),ppmap(ppmap)
  {}

  typedef typename Kernel::FT FT;
  typedef typename boost::property_traits<PointPropertyMap>::key_type Point_2;
  typedef typename boost::call_traits<Point_2>::param_type Arg_type;

  typedef ::std::list<Point_2>                      Container;
  typedef typename Polygon_traits_getter<Base_traits,PointPropertyMap>::type PolygonTraits;
  typedef CGAL::Polygon_2<PolygonTraits, Container>          Polygon_2;

  template <typename BaseFct>
  struct Pmap_fct : public BaseFct {
    Pmap_fct(const PointPropertyMap& ppmap, const BaseFct& base)
      : BaseFct(base),ppmap(ppmap)
    {}

    const PointPropertyMap& ppmap;

    typename BaseFct::result_type operator()(const Arg_type& p, const Arg_type& q) const {
      return static_cast<const BaseFct*>(this)->operator()(get(ppmap,p),get(ppmap,q));
    }

    typename BaseFct::result_type operator()(const Arg_type& p, const Arg_type& q, const Arg_type& r) const {
      return static_cast<const BaseFct*>(this)->operator()(get(ppmap,p),get(ppmap,q),get(ppmap,r));
    }
  };


  typedef Pmap_fct<typename Kernel::Equal_2>                    Equal_2;
  typedef Pmap_fct<typename Kernel::Less_yx_2>                  Less_yx_2;
  typedef Pmap_fct<typename Kernel::Less_xy_2>                  Less_xy_2;
  typedef Pmap_fct<typename Kernel::Left_turn_2>                Left_turn_2;
  typedef Pmap_fct<typename Kernel::Orientation_2>              Orientation_2;
  typedef Pmap_fct<typename Kernel::Compare_y_2>                Compare_y_2;
  typedef Pmap_fct<typename Kernel::Compare_x_2>                Compare_x_2;
  typedef CGAL::Is_convex_2<Self>                               Is_convex_2;
  typedef CGAL::Is_y_monotone_2<Self>                           Is_y_monotone_2;


  // needed by visibility graph and thus by optimal convex
  typedef Pmap_fct<typename Kernel::Collinear_are_ordered_along_line_2>   Collinear_are_ordered_along_line_2;
  typedef Pmap_fct<typename Kernel::Are_strictly_ordered_along_line_2>
                                                                Are_strictly_ordered_along_line_2;


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


  Collinear_are_ordered_along_line_2
  collinear_are_ordered_along_line_2_object() const
  { return Collinear_are_ordered_along_line_2(ppmap,static_cast<const Base_traits*>(this)->collinear_are_ordered_along_line_2_object()); }

  Are_strictly_ordered_along_line_2
  are_strictly_ordered_along_line_2_object() const
  { return Are_strictly_ordered_along_line_2(ppmap,static_cast<const Base_traits*>(this)->are_strictly_ordered_along_line_2_object()); }


  Is_convex_2
  is_convex_2_object(const Self& traits) const
  {  return Is_convex_2(traits); }

  Is_y_monotone_2
  is_y_monotone_2_object(const Self& traits) const
  {  return Is_y_monotone_2(traits); }



};

}


#endif // CGAL_PARTITION_TRAITS_2_H
