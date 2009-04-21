// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s) : Laurent Saboret

#ifndef CGAL_RADIAL_ORIENT_NORMALS_H
#define CGAL_RADIAL_ORIENT_NORMALS_H

#include <CGAL/Orientable_normal_3.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/point_set_processing_assertions.h>

#include <math.h>
#ifndef M_PI
  #define M_PI       3.14159265358979323846
#endif

CGAL_BEGIN_NAMESPACE


/// Radial orientation of the normals of a point set.
/// Normals are oriented towards exterior of the point set.
/// This very fast method is intended to convex objects.
///
/// @commentheading Preconditions:
/// - VertexIterator is a model of ForwardIterator.
/// - VertexIndexMap is a model of boost::readable_property_map.
/// - VertexPointMap is a model of boost::readable_property_map.
/// - VertexNormalMap is a model of boost::lvalue_property_map.
/// - Normals must be unit vectors.

template<class VertexIterator, class VertexPointMap, class VertexIndexMap, class VertexNormalMap>
void
radial_orient_normals(VertexIterator first, ///< first input vertex.
                      VertexIterator beyond, ///< past-the-end input vertex.
                      VertexIndexMap vertex_index_map, ///< property map VertexIterator -> index.
                      VertexPointMap vertex_point_map, ///< property map VertexIterator -> Point_3.
                      VertexNormalMap vertex_normal_map) ///< property map VertexIterator -> Normal (in and out).
{
    CGAL_TRACE("Call radial_orient_normals()\n");

    // Input mesh's types
    typedef typename boost::property_traits<VertexPointMap>::value_type Point;
    typedef typename boost::property_traits<VertexNormalMap>::value_type Normal;
    typedef typename Normal::Vector Vector;
    typedef typename Normal::Geom_traits Geom_traits; // Kernel
    typedef typename Geom_traits::FT FT;

    // Precondition: at least one element in the container.
    CGAL_point_set_processing_precondition(first != beyond);

    // Find points barycenter.
    // Note: We should use CGAL::centroid() from PCA component.
    //       Unfortunately, it is not compatible with property maps.
    Vector sum = CGAL::NULL_VECTOR;
    int nb_points = 0;
    for (VertexIterator it = first; it != beyond; it++)
    {
      Point point = get(vertex_point_map, it);
      sum = sum + (point - CGAL::ORIGIN);
      nb_points++;
    }
    Point barycenter = CGAL::ORIGIN + sum / (FT)nb_points;

    // Radial orientation of the normals of a point set.
    // Normals are oriented towards exterior of the point set.
    for (VertexIterator it = first; it != beyond; it++)
    {
      Point point = get(vertex_point_map, it);

      // Radial vector towards exterior of the point set
      Vector vec1 = point - barycenter;

      // Point's normal
      Normal& normal2 = vertex_normal_map[it];
      Vector vec2 = normal2;

      double dot = vec1 * vec2;

      //         ->               ->
      // Orient vec2 parallel to vec1
      if (dot < 0)
        vec2 = -vec2;

      // Is orientation robust?
      //normal2 = Normal(vec2, true /* oriented */);
      bool oriented = (std::abs(dot) > std::cos(80.*M_PI/180.)); // oriented iff angle < 80 degrees
      normal2 = Normal(vec2, oriented);
    }

    long memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("End of radial_orient_normals()\n");
}


CGAL_END_NAMESPACE

#endif // CGAL_RADIAL_ORIENT_NORMALS_H

