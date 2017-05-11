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
//
//
// Author(s)     : Simon Giraudot
//

/**
 * \ingroup PkgPointSetShapeDetection3
 *
 */


#ifndef CGAL_SHAPE_DETECTION_3_PROPERTY_MAPS_H
#define CGAL_SHAPE_DETECTION_3_PROPERTY_MAPS_H

namespace CGAL {

namespace Shape_detection_3 {

  template <typename ShapeDetectionTraits>
  class Point_to_shape_index_map
  {
    typedef CGAL::Shape_detection_3::Shape_base<ShapeDetectionTraits> Shape;
    boost::shared_ptr<std::vector<int> > m_indices;
    
  public:
    typedef std::size_t key_type;
    typedef int value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    template <typename ShapeRange>
    Point_to_shape_index_map (const typename ShapeDetectionTraits::Input_range& points,
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

  template <typename ShapeDetectionTraits>
  class Plane_map
  {
  public:
    typedef CGAL::Shape_detection_3::Plane<ShapeDetectionTraits> Plane_shape;
    typedef boost::shared_ptr<Plane_shape> key_type;
    typedef typename ShapeDetectionTraits::Plane_3 value_type;
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
