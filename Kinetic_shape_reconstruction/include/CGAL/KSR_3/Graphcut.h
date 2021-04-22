// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_GRAPHCUT_H
#define CGAL_KSR_3_GRAPHCUT_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/assertions.h>
#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/boost/graph/alpha_expansion_graphcut.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Cartesian_converter.h>

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

    using FT 		  	 = typename Kernel::FT;
    using Point_3 	 = typename Kernel::Point_3;
		using Triangle_2 = typename Kernel::Triangle_2;
    using Indices 	 = std::vector<std::size_t>;

    using Data_structure = KSR_3::Data_structure<Kernel>;
    using Volume_cell    = typename Data_structure::Volume_cell;
		using PFace 				 = typename Data_structure::PFace;

    using Visibility_label = KSR::Visibility_label;

		using IK         = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Delaunay_2 = CGAL::Delaunay_triangulation_2<Kernel>;
		using Delaunay_3 = CGAL::Delaunay_triangulation_3<IK>;
		using Converter  = CGAL::Cartesian_converter<Kernel, IK>;

		struct Wrapper {
			PFace pface;
			FT weight = FT(0);
			std::pair<int, int> neighbors;
			bool is_boundary = false;
		};

    Graphcut(
      const Data_structure& data,
      const FT graphcut_beta) :
    m_data(data),
    m_beta(graphcut_beta)
    { }

    void compute(std::vector<Volume_cell>& volumes) {

      if (volumes.size() == 0) return;

			std::vector<Wrapper> wrappers;
			create_pface_wrappers(wrappers);

			compute_weights(wrappers);
			compute_weights(volumes);

      std::vector< std::pair<std::size_t, std::size_t> > edges;
      std::vector<double> edge_costs;
			set_graph_edges(wrappers, edges, edge_costs);

			std::vector< std::vector<double> > cost_matrix;
			set_cost_matrix(volumes, cost_matrix);

			std::vector<std::size_t> labels;
			set_initial_labels(volumes, labels);

			compute_graphcut(edges, edge_costs, cost_matrix, labels);
			apply_new_labels(labels, volumes);
    }

  private:
    const Data_structure& m_data;
    const FT m_beta;

		void create_pface_wrappers(
			std::vector<Wrapper>& wrappers) const {

			Wrapper wrapper;
			const auto& pface_neighbors = m_data.pface_neighbors();

			wrappers.clear();
			for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
				const auto pfaces = m_data.pfaces(i);
				for (const auto pface : pfaces) {
					wrapper.pface = pface;
					wrapper.is_boundary = (i < 6) ? true : false;
					CGAL_assertion(pface_neighbors.find(pface) != pface_neighbors.end());
					const auto& pair = pface_neighbors.at(pface);
					wrapper.neighbors = pair;
					wrappers.push_back(wrapper);
				}
			}
			CGAL_assertion(wrappers.size() > 6);
		}

		void compute_weights(
			std::vector<Wrapper>& wrappers) const {

			FT sum = FT(0);
			for (auto& wrapper : wrappers) {
				auto& weight = wrapper.weight;
				const auto& pface = wrapper.pface;

				Delaunay_2 tri;
				const auto pvertices = m_data.pvertices_of_pface(pface);
        for (const auto pvertex : pvertices) {
          CGAL_assertion(m_data.has_ivertex(pvertex));
          const auto ivertex = m_data.ivertex(pvertex);
          const auto& point = m_data.point_2(pface.first, ivertex);
          tri.insert(point);
        }

				weight = FT(0);
				for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
					const Triangle_2 triangle(
						fit->vertex(0)->point(),
						fit->vertex(1)->point(),
						fit->vertex(2)->point());
					weight += triangle.area();
				}
				sum += weight;
			}

			CGAL_assertion(sum > FT(0));
			for (auto& wrapper : wrappers) {
				wrapper.weight /= sum;
			}
		}

    void compute_weights(
			std::vector<Volume_cell>& volumes) const {

			FT sum = FT(0);
			const Converter converter;

			for (auto& volume : volumes) {
				auto& weight = volume.weight;
				const auto& pvertices = volume.pvertices;

				Delaunay_3 tri;
				for (const auto& pvertex : pvertices) {
					CGAL_assertion(m_data.has_ivertex(pvertex));
					const auto ivertex = m_data.ivertex(pvertex);
					tri.insert(converter(m_data.point_3(ivertex)));
				}

				weight = FT(0);
				for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit) {
        	const auto& tet = tri.tetrahedron(cit);
					weight += tet.volume();
				}
				sum += weight;
			}

			CGAL_assertion(sum > FT(0));
			for (auto& volume : volumes) {
				volume.weight /= sum;
			}
		}

    void set_graph_edges(
      const std::vector<Wrapper>& wrappers,
      std::vector< std::pair<std::size_t, std::size_t> >& edges,
      std::vector<double>& edge_costs) const {

			edges.clear();
			edge_costs.clear();
			for (const auto& wrapper : wrappers) {

				const FT edge_weight = wrapper.weight;
				const auto& neighbors = wrapper.neighbors;

				const int idx1 = neighbors.first;
				const int idx2 = neighbors.second;

				// Boundary edges.
				CGAL_assertion(idx1 >= 0 || idx2 >= 0);
				if (idx1 < 0 && idx2 >= 0) { continue; }
				if (idx2 < 0 && idx1 >= 0) { continue; }

				// Internal edges.
				CGAL_assertion(idx1 >= 0);
				const std::size_t id1 = static_cast<std::size_t>(idx1);
				CGAL_assertion(idx2 >= 0);
				const std::size_t id2 = static_cast<std::size_t>(idx2);

				edges.push_back(std::make_pair(id1, id2));
				edge_costs.push_back(compute_edge_cost(edge_weight));
			}
		}

		const double compute_edge_cost(const FT edge_weight) const {

			CGAL_assertion(m_beta 		 >= FT(0) && m_beta 		 <= FT(1));
			CGAL_assertion(edge_weight >= FT(0) && edge_weight <= FT(1));
			return CGAL::to_double(m_beta * edge_weight);
		}

    void set_cost_matrix(
      const std::vector<Volume_cell>& volumes,
      std::vector< std::vector<double> >& cost_matrix) const {

			cost_matrix.clear();
			cost_matrix.resize(2);
			cost_matrix[0].resize(volumes.size());
			cost_matrix[1].resize(volumes.size());

			for (std::size_t i = 0; i < volumes.size(); ++i) {
				const auto& volume = volumes[i];

				const FT in  = volume.inside;
				const FT out = volume.outside;

				CGAL_assertion(in  >= FT(0) && in  <= FT(1));
				CGAL_assertion(out >= FT(0) && out <= FT(1));
				CGAL_assertion((in + out) == FT(1));

				const FT face_weight = volume.weight;
				const double cost_in  = get_face_cost(in , face_weight);
				const double cost_out = get_face_cost(out, face_weight);

				cost_matrix[0][i] = cost_in;
				cost_matrix[1][i] = cost_out;
			}
		}

		const double get_face_cost(
      const FT face_prob, const FT face_weight) const {

			CGAL_assertion(face_prob   >= FT(0) && face_prob   <= FT(1));
			CGAL_assertion(face_weight >= FT(0) && face_weight <= FT(1));

			const double weight = CGAL::to_double(face_weight);
			const double value  = (1.0 - CGAL::to_double(face_prob));
      return weight * value;
		}

    void set_initial_labels(
      const std::vector<Volume_cell>& volumes,
      std::vector<std::size_t>& labels) const {

			labels.clear();
			labels.resize(volumes.size());

			for (std::size_t i = 0; i < volumes.size(); ++i) {
				const auto& volume = volumes[i];
				if (volume.visibility == Visibility_label::INSIDE) {
					labels[i] = 0;
				} else {
					labels[i] = 1;
				}
			}
		}

		void compute_graphcut(
      const std::vector< std::pair<std::size_t, std::size_t> >& edges,
      const std::vector<double>& edge_costs,
      const std::vector< std::vector<double> >& cost_matrix,
      std::vector<std::size_t>& labels) const {

      CGAL::alpha_expansion_graphcut(
				edges, edge_costs, cost_matrix, labels);
    }

		void apply_new_labels(
      const std::vector<std::size_t>& labels,
			std::vector<Volume_cell>& volumes) const {

			CGAL_assertion(volumes.size() == labels.size());
			for (std::size_t i = 0; i < labels.size(); ++i) {
				const std::size_t label = labels[i];
				auto& volume = volumes[i];

				if (label == 0) {
					volume.visibility = Visibility_label::INSIDE;
				} else {
					volume.visibility = Visibility_label::OUTSIDE;
				}
			}
		}
  };

} // KSR_3
} // CGAL

#endif // CGAL_KSR_3_GRAPHCUT_H
