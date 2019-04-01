// Copyright (c) 2019 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
#define CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR_3/Data_structure.h>

namespace CGAL
{

template <typename GeomTraits>
class Kinetic_shape_reconstruction_3
{
public:

  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;

  typedef KSR_3::Data_structure<Kernel> Data;

private:

  Data m_data;

public:

  Kinetic_shape_reconstruction_3()
  {

  }


  template <typename PolygonRange, typename PolygonMap>
  void partition (const PolygonRange& polygons, PolygonMap polygon_map)
  {
    CGAL::Bbox_3 bbox;
    for (const auto& vt : polygons)
    {
      const std::vector<Point_3>& poly = get (polygon_map, vt);
      bbox += CGAL::bbox_3 (poly.begin(), poly.end());
    }

    m_data.add_bbox_as_polygons (bbox);
    
    m_data.add_polygons (polygons, polygon_map);
    m_data.compute_support_lines();
    m_data.make_polygons_intersection_free();
    m_data.initialize_queue();
    
    m_data.run();
  }

  
  template <typename PointRange, typename PointMap, typename VectorMap>
  void reconstruct (const PointRange& points, PointMap point_map, VectorMap normal_map)
  {

  }

  template <typename VertexOutputIterator, typename FacetOutputIterator>
  void output_partition_facets_to_polygon_soup (VertexOutputIterator vertices, FacetOutputIterator facets) const
  {
    for (std::size_t i = 0; i < m_data.number_of_vertices(); ++ i)
      *(vertices ++) = m_data.point(i);
    for (std::size_t i = 0; i < m_data.number_of_polygons(); ++ i)
      *(facets ++) = m_data.polygon(i);
  }

  template <typename OutputIterator>
  OutputIterator output_partition_cells_to_surface_meshes (OutputIterator output) const
  {

  }

};



} // namespace CGAL


#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
