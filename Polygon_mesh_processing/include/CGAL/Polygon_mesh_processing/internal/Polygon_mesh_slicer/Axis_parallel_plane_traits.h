// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Sebastien Loriot

#include <CGAL/enum.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/array.h>

#include <boost/optional.hpp>
#include <boost/variant.hpp>


#ifndef CGAL_INTERNAL_POLYGON_MESH_SLICER_AXIS_PARALLEL_PLANE_TRAITS_H
#define CGAL_INTERNAL_POLYGON_MESH_SLICER_AXIS_PARALLEL_PLANE_TRAITS_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>


namespace CGAL{
namespace Polygon_mesh_slicer_{

template <class Traits>
class Axis_parallel_plane_traits
{
  typedef typename Traits::FT FT;

  const Traits& m_traits;
  const int m_cst_coord; // 0, 1 or 2 indicates which coordinates is constant
  const FT m_value; // indicates the value of the constant coordinate

public:
  typedef typename Traits::Plane_3 Plane_3;

  Axis_parallel_plane_traits(int cst_coord, FT value, const Traits& traits)
    : m_traits(traits)
    , m_cst_coord(cst_coord)
    , m_value(value)
  {}

  struct Oriented_side_3
  {
    const int m_cst_coord;
    const FT m_value;
    const typename Traits::Construct_cartesian_const_iterator_3 m_coord_iterator;

    typedef Oriented_side result_type;

    Oriented_side_3(const Axis_parallel_plane_traits<Traits>& traits)
      : m_cst_coord(traits.m_cst_coord)
      , m_value(traits.m_value)
      , m_coord_iterator(traits.m_traits.construct_cartesian_const_iterator_3_object())
    {}

    result_type
    operator()(const typename Traits::Plane_3&, const typename Traits::Point_3& pt) const
    {
      if ( *( m_coord_iterator(pt)+m_cst_coord) > m_value ) return ON_POSITIVE_SIDE;
      if ( *( m_coord_iterator(pt)+m_cst_coord) < m_value ) return ON_NEGATIVE_SIDE;
      return ON_ORIENTED_BOUNDARY;
    }
  };

  struct Do_intersect_3{
    const int m_cst_coord;
    const FT m_value;

    typedef bool result_type;

    Do_intersect_3(int cst_coord, FT value)
      : m_cst_coord(cst_coord)
      , m_value(value)
    {}

    result_type
    operator()(const typename Traits::Plane_3&, const Bbox_3& bbox) const
    {
      return (bbox.min)(m_cst_coord) <= m_value &&
             (bbox.max)(m_cst_coord) >= m_value;
    }
  };

  struct Intersect_3{
    const int m_cst_coord;
    const FT m_value;
    const typename Traits::Construct_cartesian_const_iterator_3 m_coord_iterator;
    const typename Traits::Construct_point_3 m_point_3;
    const typename Traits::Construct_source_3 m_source_3;
    const typename Traits::Construct_target_3 m_target_3;

    typedef boost::variant<typename Traits::Point_3, typename Traits::Segment_3> Variant_type;
    typedef boost::optional< Variant_type > result_type;

    Intersect_3(const Axis_parallel_plane_traits<Traits>& traits)
      : m_cst_coord(traits.m_cst_coord)
      , m_value(traits.m_value)
      , m_coord_iterator(traits.m_traits.construct_cartesian_const_iterator_3_object())
      , m_point_3(traits.m_traits.construct_point_3_object())
      , m_source_3(traits.m_traits.construct_source_3_object())
      , m_target_3(traits.m_traits.construct_target_3_object())
    {}

    result_type
    operator()( const typename Traits::Plane_3&,
                const typename Traits::Segment_3& s) const
    {
      const typename Traits::Point_3& src = m_source_3(s);
      const typename Traits::Point_3& tgt = m_target_3(s);

      std::array<FT,3> src_coords = {{ *m_coord_iterator(src),
                                         *(m_coord_iterator(src)+1),
                                         *(m_coord_iterator(src)+2) }};
      std::array<FT,3> tgt_coords = {{ *m_coord_iterator(tgt),
                                         *(m_coord_iterator(tgt)+1),
                                         *(m_coord_iterator(tgt)+2) }};

      FT alpha = ( m_value - src_coords[m_cst_coord] ) / ( tgt_coords[m_cst_coord] - src_coords[m_cst_coord] );
      src_coords[m_cst_coord]=m_value;
      for (int i=1;i<3;++i)
      {
        int index = (m_cst_coord+i)%3;
        src_coords[index]+=(tgt_coords[index]-src_coords[index])*alpha;
      }

      return Variant_type( m_point_3(src_coords[0], src_coords[1], src_coords[2]) );
    }
  };

  Oriented_side_3 oriented_side_3_object() const
  {
    return Oriented_side_3(*this);
  }

  Do_intersect_3 do_intersect_3_object() const
  {
    return Do_intersect_3(m_cst_coord, m_value);
  }

  Intersect_3 intersect_3_object() const
  {
    return Intersect_3(*this);
  }
};

} } // end of namespace CGAL::Polygon_mesh_slicer_


#endif // CGAL_INTERNAL_POLYGON_MESH_SLICER_AXIS_PARALLEL_PLANE_TRAITS_H
