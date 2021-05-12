// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_PARAMETERIZE_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_PARAMETERIZE_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>

#include <boost/property_map/property_map.hpp>

/// \file parameterize.h

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceMeshParameterizationMainFunction
///
/// computes a one-to-one mapping from a 3D triangle surface `mesh` to a
/// simple 2D domain.
/// The mapping is piecewise linear on the triangle mesh.
/// The result is a pair `(u,v)` of parameter coordinates for each vertex of the input mesh.
///
/// A one-to-one mapping may be guaranteed or not, depending on
/// the chosen Parameterizer algorithm.
///
/// \tparam TriangleMesh must be a model of `FaceGraph`.
/// \tparam Parameterizer must be a model of `Parameterizer_3`.
/// \tparam HD must be the halfedge_descriptor type corresponding to the graph
///         traits of TriangleMesh.
/// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
///         %Point_2 (type deduced from `TriangleMesh` using `Kernel_traits`)
///         as value type.
///
/// \param mesh a triangulated surface.
/// \param parameterizer a parameterizer.
/// \param bhd a halfedge descriptor on the boundary of `mesh`.
/// \param uvmap an instanciation of the class `VertexUVmap`.
///
/// \pre `mesh` must be a triangular mesh.
/// \pre The mesh border must be mapped onto a convex polygon
///   (for fixed border parameterizations).
/// \pre The vertices must be indexed (vimap must be initialized).
///
template <class TriangleMesh, class Parameterizer, class HD, class VertexUVmap>
Error_code parameterize(TriangleMesh& mesh,
                        Parameterizer parameterizer,
                        HD bhd,
                        VertexUVmap uvmap)
{
  CGAL_precondition(is_valid_polygon_mesh(mesh));
  CGAL_precondition(bhd != boost::graph_traits<TriangleMesh>::null_halfedge() && is_border(bhd, mesh));

  typedef CGAL::dynamic_vertex_property_t<int>                                 Vertex_int_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_int_tag>::type     Vertex_int_map;
  Vertex_int_map vimap = get(Vertex_int_tag(), mesh);
  internal::fill_index_map_of_cc(bhd, mesh, vimap);

  typedef CGAL::dynamic_vertex_property_t<bool>                                Vertex_bool_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_bool_tag>::type    Vertex_bool_map;
  Vertex_bool_map vpmap = get(Vertex_bool_tag(), mesh);

  return parameterizer.parameterize(mesh, bhd, uvmap, vimap, vpmap);
}

/// \ingroup  PkgSurfaceMeshParameterizationMainFunction
///
/// computes a one-to-one mapping from a 3D triangle surface `mesh` to a
/// 2D circle, using Floater Mean Value Coordinates algorithm.
/// A one-to-one mapping is guaranteed.
///
/// The mapping is piecewise linear on the input mesh triangles.
/// The result is a `(u,v)` pair of parameter coordinates for each vertex of the input mesh.
///
/// \tparam TriangleMesh must be a model of `FaceGraph`.
/// \tparam HD must be the halfedge_descriptor type corresponding to the graph
///         traits of TriangleMesh.
/// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
///         %Point_2 (type deduced from `TriangleMesh` using `Kernel_traits`)
///         as value type.
///
/// \param mesh a triangulated surface.
/// \param bhd a halfedge descriptor on the boundary of `mesh`.
/// \param uvmap an instanciation of the class `VertexUVmap`.
///
/// \pre `mesh` must be a triangular mesh.
/// \pre The vertices must be indexed (vimap must be initialized).
///
template <class TriangleMesh, class HD, class VertexUVmap>
Error_code parameterize(TriangleMesh& mesh,
                        HD bhd,
                        VertexUVmap uvmap)
{
  Mean_value_coordinates_parameterizer_3<TriangleMesh> parameterizer;
  return parameterize(mesh, parameterizer, bhd, uvmap);
}

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_PARAMETERIZE_H
