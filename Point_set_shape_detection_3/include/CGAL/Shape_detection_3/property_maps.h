// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot
//

#include <CGAL/license/Point_set_shape_detection_3.h>
#include <CGAL/Shape_detection_3/Shape_base.h>
#include <CGAL/Shape_detection_3.h>

#ifndef CGAL_SHAPE_DETECTION_3_PROPERTY_MAPS_H
#define CGAL_SHAPE_DETECTION_3_PROPERTY_MAPS_H

namespace CGAL {

namespace Shape_detection_3 {

/*!
   \ingroup PkgPointSetShapeDetection3

   Property map that associate a point index to its assigned shape
   found by a shape detection algorithm (such as
   `CGAL::Shape_detection_3::Efficient_RANSAC` or
   `CGAL::Shape_detection_3::Region_growing`).
 */
  template <typename Traits>
  class Point_to_shape_index_map
  {
    typedef CGAL::Shape_detection_3::Shape_base<Traits> Shape;
    boost::shared_ptr<std::vector<int> > m_indices;
    
  public:
    typedef std::size_t key_type; ///< Index of the point in the random access point range.
    typedef int value_type; ///< Index of the shape (-1 if the point is not assigned to any shape).
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    /// \cond SKIP_IN_MANUAL
    Point_to_shape_index_map () { }
    /// \endcond

    /*!
      Constructs a property map to map points to their associated shape.

      \note `shapes` must be a range of shapes detected using `points`.

      \tparam ShapeRange an `Iterator_range` with a bidirectional
      constant iterator type with value type
      `boost::shared_ptr<CGAL::Shape_detection_3::Shape_base<Traits> >`.
     */
    template <typename PointRange, typename ShapeRange>
    Point_to_shape_index_map (const PointRange& points,
                              const ShapeRange& shapes)
      : m_indices (new std::vector<int>(points.size(), -1))
    {
      int idx = 0;
      for (typename ShapeRange::const_iterator it = shapes.begin();
             it != shapes.end(); ++ it)
      {
        for (std::size_t j = 0; j < (*it)->indices_of_assigned_points ().size (); ++ j)
          (*m_indices)[(*it)->indices_of_assigned_points()[j]] = idx;
        ++ idx;
      }
    }

    inline friend value_type get (const Point_to_shape_index_map& pm, const key_type& k)
    {
      return (*(pm.m_indices))[k];
    }

  };

/*!
   \ingroup PkgPointSetShapeDetection3

   Property map that associates a detected plane object
   (`CGAL::Shape_detection_3::Plane`) to a `CGAL::Plane_3` object.
 */
  template <typename Traits>
  class Plane_map
  {
  public:
    typedef CGAL::Shape_detection_3::Plane<Traits> Plane_shape;
    typedef boost::shared_ptr<Plane_shape> key_type;
    typedef typename Traits::Plane_3 value_type;
    typedef value_type reference;
    typedef boost::read_write_property_map_tag category;

    inline friend reference get (const Plane_map&, const key_type& k)
    {
      return value_type(*k);
    }

    inline friend void put (const Plane_map&, const key_type& k, const value_type& v)
    {
      k->update(v);
    }

  };


} // namespace Shape_detection_3

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_3_PROPERTY_MAPS_H
