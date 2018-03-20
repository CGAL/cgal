// Copyright (c) 2005  INRIA (France).
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
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_PARAMETERIZE_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_PARAMETERIZE_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/Bool_property_map.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <boost/function_output_iterator.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

/// \file parameterize.h

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceParameterizationMainFunction
///
/// Compute a one-to-one mapping from a 3D triangle surface `mesh` to a
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
/// \param uvm an instanciation of the class `VertexUVmap`.
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
                        VertexUVmap uvm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  typedef boost::unordered_map<vertex_descriptor, int> Indices;
  Indices indices;
  CGAL::Polygon_mesh_processing::connected_component(
         face(opposite(bhd, mesh), mesh),
         mesh,
         boost::make_function_output_iterator(
         internal::Index_map_filler<TriangleMesh, Indices>(mesh, indices)));
  boost::associative_property_map<Indices> vipm(indices);

  boost::unordered_set<vertex_descriptor> vs;
  internal::Bool_property_map<boost::unordered_set<vertex_descriptor> > vpm(vs);

  return parameterizer.parameterize(mesh, bhd, uvm, vipm, vpm);
}

/// \ingroup  PkgSurfaceParameterizationMainFunction
///
/// Compute a one-to-one mapping from a 3D triangle surface `mesh` to a
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
/// \param uvm an instanciation of the class `VertexUVmap`.
///
/// \pre `mesh` must be a triangular mesh.
/// \pre The vertices must be indexed (vimap must be initialized).
///
template <class TriangleMesh, class HD, class VertexUVmap>
Error_code parameterize(TriangleMesh& mesh,
                        HD bhd,
                        VertexUVmap uvm)
{
  Mean_value_coordinates_parameterizer_3<TriangleMesh> parameterizer;
  return parameterize(mesh, parameterizer, bhd, uvm);
}

template <class TM, class SEM, class SVM>
class Seam_mesh;

template <class TM, class SEM, class SVM, class Parameterizer, class HD, class VertexUVmap>
Error_code parameterize(Seam_mesh<TM, SEM, SVM>& mesh,
                        Parameterizer parameterizer,
                        HD bhd,
                        VertexUVmap uvm)
{
  typedef typename boost::graph_traits<Seam_mesh<TM, SEM, SVM> >::vertex_descriptor vertex_descriptor;
  boost::unordered_set<vertex_descriptor> vs;
  internal::Bool_property_map<boost::unordered_set<vertex_descriptor> > vpm(vs);

  typedef boost::unordered_map<vertex_descriptor, int> Indices;
  Indices indices;
  CGAL::Polygon_mesh_processing::connected_component(
         face(opposite(bhd, mesh), mesh),
         mesh,
         boost::make_function_output_iterator(
         internal::Index_map_filler<Seam_mesh<TM, SEM, SVM>,
                                    Indices>(mesh, indices)));
  boost::associative_property_map<Indices> vipm(indices);

  return parameterizer.parameterize(mesh, bhd, uvm, vipm, vpm);
}

template <class TM, class SEM, class SVM, class HD, class VertexUVmap>
Error_code parameterize(Seam_mesh<TM, SEM, SVM>& mesh,
                        HD bhd,
                        VertexUVmap uvm)
{
  Mean_value_coordinates_parameterizer_3<Seam_mesh<TM, SEM, SVM> > parameterizer;
  return parameterize(mesh, parameterizer, bhd, uvm);
}

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_PARAMETERIZE_H
