// Copyright (c) 2026 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ange Clement

#ifndef CGAL_FEATURE_GRAPH_INTERNAL_DEBUGGING_METHODS_H
#define CGAL_FEATURE_GRAPH_INTERNAL_DEBUGGING_METHODS_H

#include <cstddef>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

#include <CGAL/license/Feature_graph.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/range/iterator_range.hpp>

#include <CGAL/Feature_graph/internal/feature_graph_operation.h>

namespace CGAL
{

namespace Feature_graph
{

namespace internal
{

namespace debug
{

template <typename K,
  typename PolygonMesh, typename MeshPointPropertyMap,
  typename SharpnessEstimator, typename ElementAndDomainTranslator,
  typename NormalEstimator>
void output_sharpness_on_surface_elements(
  const std::string& filename,
  const PolygonMesh& pmesh, const MeshPointPropertyMap& points,
  const SharpnessEstimator& sharpness_estimator, const ElementAndDomainTranslator& element_and_domain_translator,
  const NormalEstimator& normal_estimator)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  auto sharpness_estimator_domain = element_and_domain_translator.translate_domain(pmesh);
  std::vector<Point_3> out_points;
  std::vector<double> out_sharpness;
  double min_sharpness = 0;
  double max_sharpness = 0;
  std::set<vertex_descriptor> inserted_vertices;
  std::set<edge_descriptor> inserted_edges;
  for (const face_descriptor& face : faces(pmesh))
  {
    halfedge_descriptor face_edge = halfedge(face, pmesh);
    auto face_element_sharpness = sharpness_estimator.template operator()<Dimension_tag<2>>(element_and_domain_translator.translate_face_element(face, pmesh, points), sharpness_estimator_domain);
    if (face_element_sharpness < min_sharpness)
      min_sharpness = face_element_sharpness;
    if (face_element_sharpness > max_sharpness)
      max_sharpness = face_element_sharpness;
    Vector_3 barycenter(0,0,0);
    int number_of_edges = 0;
    halfedge_descriptor current_face_edge = face_edge;
    do
    {
      vertex_descriptor source_vertex = source(current_face_edge, pmesh);
      bool is_source_vertex_inserted = inserted_vertices.insert(source_vertex).second;
      if (is_source_vertex_inserted)
      {
        auto source_vertex_sharpness = sharpness_estimator.template operator()<Dimension_tag<0>>(element_and_domain_translator.translate_vertex_element(source_vertex, pmesh, points), sharpness_estimator_domain);
        out_points.push_back(get(points, source_vertex));
        out_sharpness.push_back(source_vertex_sharpness);
        if (source_vertex_sharpness < min_sharpness)
          min_sharpness = source_vertex_sharpness;
        if (source_vertex_sharpness > max_sharpness)
          max_sharpness = source_vertex_sharpness;
      }
      edge_descriptor current_as_edge = edge(current_face_edge, pmesh);
      bool is_edge_inserted = inserted_edges.insert(current_as_edge).second;
      if (is_edge_inserted)
      {
        Point_3 edge_point = CGAL::midpoint(get(points, source_vertex), get(points, target(current_face_edge, pmesh)));
        auto edge_sharpness = sharpness_estimator.template operator()<Dimension_tag<1>>(element_and_domain_translator.translate_halfedge_element(current_face_edge, pmesh, points), sharpness_estimator_domain);
        out_points.push_back(edge_point);
        out_sharpness.push_back(edge_sharpness);
        if (edge_sharpness < min_sharpness)
          min_sharpness = edge_sharpness;
        if (edge_sharpness > max_sharpness)
          max_sharpness = edge_sharpness;
      }
      barycenter += Vector_3(CGAL::ORIGIN, get(points, source_vertex));
      number_of_edges++;
      current_face_edge = next(current_face_edge, pmesh);
    } while (current_face_edge != face_edge);
    if (number_of_edges != 0)
      barycenter /= number_of_edges;
    out_points.push_back(CGAL::ORIGIN + barycenter);
    out_sharpness.push_back(face_element_sharpness);
  }
  std::ofstream points_output(filename);
  points_output << "ply\nformat ascii 1.0\n"
                    << "element vertex " << out_points.size() << "\n"
                    << "property double x\nproperty double y\nproperty double z\n"
                    << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
                    << "end_header\n";
  for (std::size_t i = 0; i < out_points.size(); i++)
  {
    const double& cur_sharpness = out_sharpness[i];
    if (cur_sharpness >= 0)
    {
      double scaled_sharpness = cur_sharpness;
      if (max_sharpness != 0)
        scaled_sharpness /= max_sharpness;
      int sharpness = static_cast<int>(scaled_sharpness * 255.0);
      std::clamp(sharpness, 0, 255);
      points_output << out_points[i] << " " << sharpness << " 0 0" << std::endl;
    }
    else
    {
      double scaled_sharpness = cur_sharpness;
      if (min_sharpness != 0)
        scaled_sharpness /= min_sharpness;
      int sharpness = static_cast<int>(scaled_sharpness * 255.0);
      std::clamp(sharpness, 0, 255);
      points_output << out_points[i] << " " << "0 0 " << sharpness << std::endl;
    }
  }
  points_output.close();
}

template <typename FeatureGraph,
  typename VertexPointMap, typename VertexDimensionMap,
  typename VertexSharpnessMap, typename VertexNextToHoleMap>
void output_feature_graph_properties(
  const std::string& filename,
  const FeatureGraph& graph,
  const VertexPointMap& vertex_point_map,
  const VertexDimensionMap& vertex_dimension_map,
  const VertexSharpnessMap& vertex_sharpness_map,
  const VertexNextToHoleMap& vertex_is_next_to_hole_map
)
{
  // [TODO] create single ply with points and edges (when CGAL-Lab can display them)
  typedef typename boost::graph_traits<FeatureGraph>::vertex_descriptor vertex_descriptor;
  std::ofstream points_output(filename+"_points.ply");
  std::ofstream edges_output(filename+"_edges.polylines.txt");
  points_output << "ply\nformat ascii 1.0\n"
                   << "element vertex " << num_vertices(graph) << "\n"
                   << "property double x\nproperty double y\nproperty double z\n"
                   << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
                   << "end_header\n";
  for (auto vertex : boost::make_iterator_range(vertices(graph)))
  {
    int dim = static_cast<float>(get(vertex_dimension_map, vertex)) * 255.0 / 3.0;
    int sharpness = static_cast<float>(get(vertex_sharpness_map, vertex)) * 255.0;
    int is_next_to_hole = static_cast<float>(get(vertex_is_next_to_hole_map, vertex)) * 255.0;
    points_output << get(vertex_point_map, vertex) << " " << dim << " " << sharpness << " " << is_next_to_hole << std::endl;
  }
  for (auto edge : boost::make_iterator_range(edges(graph)))
    edges_output << "2 " << get(vertex_point_map, source(edge, graph)) << " " << get(vertex_point_map, target(edge, graph)) << std::endl;
  points_output.close();
  edges_output.close();
}

} /* namespace debug */

} /* namespace internal */

} /* namespace Feature_graph */

} /* namespace CGAL */

#endif // CGAL_FEATURE_GRAPH_INTERNAL_DEBUGGING_METHODS_H