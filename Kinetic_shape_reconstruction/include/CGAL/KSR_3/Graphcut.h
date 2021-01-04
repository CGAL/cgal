// Copyright (c) 2019 GeometryFactory SARL (France).
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

#ifndef CGAL_KSR_3_GRAPHCUT_H
#define CGAL_KSR_3_GRAPHCUT_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/assertions.h>
#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/boost/graph/alpha_expansion_graphcut.h>

// Internal includes.
#include <CGAL/KSR/enum.h>
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR/property_map.h>
#include <CGAL/KSR_3/Data_structure.h>

namespace CGAL {
namespace KSR_3 {

  template<typename GeomTraits>
  class Graphcut {

  public:
    using Kernel = GeomTraits;

    using FT 		  = typename Kernel::FT;
    using Point_3 = typename Kernel::Point_3;
    using Indices = std::vector<std::size_t>;

    using Data_structure = KSR_3::Data_structure<Kernel>;
    using Volume_cell    = typename Data_structure::Volume_cell;

    using Visibility_label = KSR::Visibility_label;
		// using Alpha_expansion = CGAL::Alpha_expansion_boost_compressed_sparse_row_impl;

    Graphcut(
      const Data_structure& data,
      const FT graphcut_beta) :
    m_data(data),
    m_beta(graphcut_beta)
    { }

    void compute(std::vector<Volume_cell>& volumes) {

      // if (partition.empty()) return;
      // auto& pfaces = partition.faces;
      // auto& pedges = partition.edges;

			// compute_weights(pfaces);
			// compute_weights(pedges);

      // std::vector<Size_pair> edges;
      // std::vector<double> edge_weights;
			// set_graph_edges(pedges, edges, edge_weights);

			// std::vector< std::vector<double> > cost_matrix;
			// set_cost_matrix(pfaces, cost_matrix);

			// std::vector<std::size_t> labels;
			// set_initial_labels(pfaces, labels);

			// compute_graphcut(edges, edge_weights, cost_matrix, labels);
			// apply_new_labels(labels, pfaces);

      CGAL_assertion_msg(false, "TODO: FINISH GRAPHCUT!");
    }

  private:
    const Data_structure& m_data;
    const FT m_beta;

    // template<typename Object>
		// void compute_weights(
		// 	std::vector<Object>& objects) const {

		// 	FT sum = FT(0); FT max_value = FT(-1);
		// 	for (auto& object : objects) {
		// 		object.compute_weight();
		// 		const FT weight = CGAL::abs(object.weight);
		// 		sum += weight;
		// 		max_value = CGAL::max(weight, max_value);
		// 	}
		// 	CGAL_assertion(sum > FT(0));

		// 	if (m_use_max) {

		// 		for (auto& object : objects)
		// 			object.weight /= max_value;

		// 	} else {

		// 		for (auto& object : objects)
		// 			object.weight /= sum;
		// 	}
		// }

    // void set_graph_edges(
    //   const std::vector<Edge>& pedges,
    //   std::vector<Size_pair>& edges,
    //   std::vector<double>& edge_weights) const {

		// 	edges.clear();
		// 	edge_weights.clear();
		// 	for (const auto& pedge : pedges) {

		// 		const FT edge_weight = pedge.weight;
		// 		const auto& neighbors = pedge.neighbors;

		// 		const int idx1 = neighbors.first;
		// 		const int idx2 = neighbors.second;

		// 		// Boundary edges.
		// 		if (idx1 < 0 && idx2 >= 0)
		// 			continue;
		// 		if (idx2 < 0 && idx1 >= 0)
		// 			continue;

		// 		// Internal edges.
		// 		CGAL_assertion(idx1 >= 0);
		// 		const std::size_t id1 = static_cast<std::size_t>(idx1);
		// 		CGAL_assertion(idx2 >= 0);
		// 		const std::size_t id2 = static_cast<std::size_t>(idx2);

		// 		CGAL_assertion(edge_weight >= 0.0);
		// 		edges.push_back(std::make_pair(id1, id2));
		// 		edge_weights.push_back(get_graph_edge_cost(edge_weight));
		// 	}
		// }

		// double get_graph_edge_cost(const FT edge_weight) const {
		// 	return CGAL::to_double(m_beta * edge_weight);
		// }

    // void set_cost_matrix(
    //   const std::vector<Face>& pfaces,
    //   std::vector< std::vector<double> >& cost_matrix) const {

		// 	cost_matrix.clear();
		// 	cost_matrix.resize(2);
		// 	cost_matrix[0].resize(pfaces.size());
		// 	cost_matrix[1].resize(pfaces.size());

		// 	for (std::size_t i = 0; i < pfaces.size(); ++i) {
		// 		const auto& pface = pfaces[i];

		// 		const FT in = pface.inside;
		// 		const FT out = pface.outside;

		// 		CGAL_assertion(in >= FT(0) && in <= FT(1));
		// 		CGAL_assertion(out >= FT(0) && out <= FT(1));
		// 		CGAL_assertion((in + out) == FT(1));

		// 		const FT face_weight = pface.weight;
		// 		CGAL_precondition(face_weight >= FT(0));

		// 		const double cost_in = get_graph_face_cost(in, face_weight);
		// 		const double cost_out = get_graph_face_cost(out, face_weight);

		// 		cost_matrix[0][i] = cost_in;
		// 		cost_matrix[1][i] = cost_out;
		// 	}
		// }

		// double get_graph_face_cost(
    //   const FT face_prob, const FT face_weight) const {

		// 	const double weight = CGAL::to_double(face_weight);
		// 	const double value = (1.0 - CGAL::to_double(face_prob));
    //   return weight * value;
		// }

    // void set_initial_labels(
    //   const std::vector<Face>& pfaces,
    //   std::vector<std::size_t>& labels) const {

		// 	labels.clear();
		// 	labels.resize(pfaces.size());

		// 	for (std::size_t i = 0; i < pfaces.size(); ++i) {
		// 		if (pfaces[i].visibility == Visibility_label::INSIDE) labels[i] = 0;
		// 		else labels[i] = 1;
		// 	}
		// }

		// void compute_graphcut(
    //   const std::vector<Size_pair>& edges,
    //   const std::vector<double>& edge_weights,
    //   const std::vector< std::vector<double> >& cost_matrix,
    //   std::vector<std::size_t>& labels) const {

    //   Alpha_expansion graphcut;
    //   graphcut(edges, edge_weights, cost_matrix, labels);
    // }

		// void apply_new_labels(
    //   const std::vector<std::size_t>& labels,
		// 	std::vector<Face>& pfaces) const {

		// 	for (std::size_t i = 0; i < labels.size(); ++i) {
		// 		if (labels[i] == 0) pfaces[i].visibility = Visibility_label::INSIDE;
		// 		else pfaces[i].visibility = Visibility_label::OUTSIDE;
		// 	}
		// }
  };

} // KSR_3
} // CGAL

#endif // CGAL_KSR_3_GRAPHCUT_H
