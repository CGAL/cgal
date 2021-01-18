// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_TWO_VERTICES_PARAMETERIZER_3_H_INCLUDED
#define CGAL_SURFACE_MESH_PARAMETERIZATION_TWO_VERTICES_PARAMETERIZER_3_H_INCLUDED

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <boost/function_output_iterator.hpp>

#include <cfloat>
#include <climits>

/// \file Two_vertices_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

//
// Declaration
//

/// \ingroup PkgSurfaceMeshParameterizationBorderParameterizationMethods
///
/// The class `Two_vertices_parameterizer_3` parameterizes two extreme vertices
/// of a 3D surface.
/// This kind of border parameterization is used by free border parameterizations.
///
/// Implementation note:
/// To simplify the implementation, the border parameterizer knows only the
/// `TriangleMesh` class and does not know the parameterization algorithm
/// requirements or the kind of sparse linear system used.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
template< typename TriangleMesh_ >
class Two_vertices_parameterizer_3
{
// Public types
public:
  /// Triangle mesh type
  typedef TriangleMesh_                                            Triangle_mesh;

  typedef TriangleMesh_                                            TriangleMesh;

  /// Mesh vertex type
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;

  /// Mesh halfedge type
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor halfedge_descriptor;

// Private types
private:
  // Traits subtypes:
  typedef typename internal::Kernel_traits<Triangle_mesh>::PPM      PPM;
  typedef typename internal::Kernel_traits<Triangle_mesh>::Kernel   Kernel;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Point_3                                  Point_3;
  typedef typename Kernel::Vector_3                                 Vector_3;

  vertex_descriptor vxmin, vxmax;
  bool vertices_given;

// Public operations
public:
  // Default constructor, copy constructor and operator =() are fine.

  /// Constructor.
  Two_vertices_parameterizer_3()
    : vertices_given(false)
  { }

  /// Constructor where fixed vertices are provided.
  Two_vertices_parameterizer_3(vertex_descriptor v1, vertex_descriptor v2)
    : vxmin(v1), vxmax(v2), vertices_given(true)
  { }

  template <typename VertexContainer,
            typename VertexUVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize(const Triangle_mesh& mesh,
                          const VertexContainer& vertices,
                          VertexUVmap uvmap,
                          VertexIndexMap /* vimap */,
                          VertexParameterizedMap vpmap)
  {
    if(vertices_given) {
      bool found_min = false, found_max = false;
      for(vertex_descriptor vd : vertices) {
        if(vd == vxmin) {
          found_min = true;
          if(found_max) break;
        }

        if(vd == vxmax) {
          found_max = true;
          if(found_min) break;
        }
      }

      if(!found_min || !found_max) {
        std::cerr << "Error: Fixed vertices must be in the same connected component" << std::endl;
        return ERROR_NON_CONVEX_BORDER;
      }

      put(uvmap, vxmin, Point_2(0, 0.5));
      put(uvmap, vxmax, Point_2(1, 0.5));
      put(vpmap, vxmin, true);
      put(vpmap, vxmax, true);
      return OK;
    }

    const PPM ppmap = get(vertex_point, mesh);

    // Get mesh's bounding box
    double xmin = std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double zmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();

    for(vertex_descriptor vd : vertices) {
      const Point_3& position = get(ppmap,vd);

      xmin = (std::min)(position.x(), xmin);
      ymin = (std::min)(position.y(), ymin);
      zmin = (std::min)(position.z(), zmin);

      xmax = (std::max)(position.x(), xmax);
      ymax = (std::max)(position.y(), ymax);
      zmax = (std::max)(position.z(), zmax);
    }

    // Find longest bounding box axes
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    double dz = zmax - zmin;
    enum { X_AXIS, Y_AXIS, Z_AXIS } longest_axis, second_longest_axis;
    if(dx < dy && dx < dz) {
      if(dy > dz) {
        longest_axis        = Y_AXIS;
        second_longest_axis = Z_AXIS;
      } else {
        longest_axis        = Z_AXIS;
        second_longest_axis = Y_AXIS;
      }
    } else if(dy < dx && dy < dz) {
      if(dx > dz) {
        longest_axis        = X_AXIS;
        second_longest_axis = Z_AXIS;
      } else {
        longest_axis        = Z_AXIS;
        second_longest_axis = X_AXIS;
      }
    } else { // (dz < dx && dz < dy)
      if(dx > dy) {
        longest_axis        = X_AXIS;
        second_longest_axis = Y_AXIS;
      } else {
        longest_axis        = Y_AXIS;
        second_longest_axis = X_AXIS;
      }
    }

    Vector_3 V1,                // bounding box' longest axis
             V2;                // bounding box' 2nd longest axis
    double V1_min=0, V1_max=0;  // bounding box' dimensions along V1
    double V2_min=0, V2_max=0;  // bounding box' dimensions along V2
    switch (longest_axis)
    {
    case X_AXIS:
      V1 = Vector_3(1, 0, 0);
      V1_min = xmin;
      V1_max = xmax;
      break;
    case Y_AXIS:
      V1 = Vector_3(0, 1, 0);
      V1_min = ymin;
      V1_max = ymax;
      break;
    case Z_AXIS:
      V1 = Vector_3(0, 0, 1);
      V1_min = zmin;
      V1_max = zmax;
      break;
    default:
      CGAL_assertion(false);
    }
    switch (second_longest_axis)
    {
    case X_AXIS:
      V2 = Vector_3(1, 0, 0);
      V2_min = xmin;
      V2_max = xmax;
      break;
    case Y_AXIS:
      V2 = Vector_3(0, 1, 0);
      V2_min = ymin;
      V2_max = ymax;
      break;
    case Z_AXIS:
      V2 = Vector_3(0, 0, 1);
      V2_min = zmin;
      V2_max = zmax;
      break;
    default:
      CGAL_assertion(false);
    }

    // Project onto longest bounding box axes,
    // Set extrema vertices' (u,v) in unit square and mark them as "parameterized"
    double umin = std::numeric_limits<double>::infinity();
    double umax = -std::numeric_limits<double>::infinity();
    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();

    for(vertex_descriptor vd : vertices) {
      const Point_3& position = get(ppmap, vd);
      Vector_3 position_as_vector = position - Point_3(0, 0, 0);

      // coordinate along the bounding box' main axes
      double u = position_as_vector * V1;
      double v = position_as_vector * V2;

      // convert to unit square coordinates
      CGAL_assertion(V1_max > V1_min);
      CGAL_assertion(V2_max > V2_min);
      u = (u - V1_min) / (V1_max - V1_min);
      v = (v - V2_min) / (V2_max - V2_min);

      if(u < umin || (u==umin && v < vmin)) {
        vxmin = vd;
        umin = u;
        vmin = v;
      }
      if(u > umax || (u==umax && v > vmax)) {
        vxmax = vd;
        umax = u;
        vmax = v;
      }
    }

    put(uvmap, vxmin, Point_2(umin, vmin)); // useful only for vxmin and vxmax
    put(uvmap, vxmax, Point_2(umax, vmax)); // useful only for vxmin and vxmax
    put(vpmap, vxmin, true);
    put(vpmap, vxmax, true);

#ifdef DEBUG_TRACE
    std::cerr << "  map two vertices..." << std::endl;
#endif

    return OK;
  }

  /// maps two extreme vertices of the 3D mesh and mark them as <i>parameterized</i>.
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `Triangle_mesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  /// \tparam VertexParameterizedMap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         a Boolean as value type.
  ///
  /// \param mesh a triangulated surface.
  /// \param bhd a halfedge descriptor on the boundary of `mesh`.
  /// \param uvmap an instanciation of the class `VertexUVmap`.
  /// \param vimap an instanciation of the class `VertexIndexMap`.
  /// \param vpmap an instanciation of the class `VertexParameterizedMap`.
  ///
  /// \pre `mesh` must be a triangular mesh.
  /// \pre The vertices must be indexed (vimap must be initialized).
  ///
  template <typename VertexUVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize(const Triangle_mesh& mesh,
                          halfedge_descriptor bhd,
                          VertexUVmap uvmap,
                          VertexIndexMap vimap,
                          VertexParameterizedMap vpmap)
  {
    // Fill containers
    boost::unordered_set<vertex_descriptor> vertices;
    internal::Containers_filler<Triangle_mesh> fc(mesh, vertices);
    Polygon_mesh_processing::connected_component(
                                      face(opposite(bhd, mesh), mesh),
                                      mesh,
                                      boost::make_function_output_iterator(fc));

    return parameterize(mesh, vertices, uvmap, vimap, vpmap);
  }

  /// indicates if the border's shape is convex.
  /// Meaningless for free border parameterization algorithms.
  bool is_border_convex() const { return false; }
};

} // Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_TWO_VERTICES_PARAMETERIZER_3_H_INCLUDED
