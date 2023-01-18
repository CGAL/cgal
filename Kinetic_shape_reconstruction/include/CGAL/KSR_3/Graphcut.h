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
// Author(s)     : Simon Giraudot, Dmitry Anisimov

#ifndef CGAL_KSR_3_GRAPHCUT_H
#define CGAL_KSR_3_GRAPHCUT_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/assertions.h>
//#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/boost/graph/alpha_expansion_graphcut.h>
#include <CGAL/boost/graph/Alpha_expansion_MaxFlow_tag.h>
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

    using FT          = typename Kernel::FT;
    using Point_3    = typename Kernel::Point_3;
    using Triangle_2 = typename Kernel::Triangle_2;
    using Triangle_3 = typename Kernel::Triangle_3;
    using Indices    = std::vector<std::size_t>;

    using Data_structure = KSR_3::Data_structure<Kernel>;
    using Mesh           = typename Data_structure::Mesh;
    using Volume_cell    = typename Data_structure::Volume_cell;
    using PFace          = typename Data_structure::PFace;

    using Visibility_label = KSR::Visibility_label;

    using IK         = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Delaunay_2 = CGAL::Delaunay_triangulation_2<Kernel>;
    using Delaunay_3 = CGAL::Delaunay_triangulation_3<IK>;
    using Converter  = CGAL::Cartesian_converter<Kernel, IK>;

    struct Wrapper {
      PFace pface;
      FT weight = FT(0);
      std::pair<int, int> neighbors;
      std::size_t support_plane;
    };

    Graphcut(
      const Data_structure& data,
      const FT graphcut_beta) :
    m_data(data),
    m_beta(graphcut_beta)
    { }

    void compute(std::vector<Volume_cell>& volumes, std::size_t inliers) {

      if (volumes.size() == 0) return;

      std::vector<Wrapper> wrappers;
      create_pface_wrappers(wrappers);

      compute_weights(wrappers, inliers);
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
    FT m_total_area;

    void create_pface_wrappers(
      std::vector<Wrapper>& wrappers) const {

      Wrapper wrapper;
      const auto& pface_neighbors = m_data.pface_neighbors();

      wrappers.clear();
      for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
        const auto pfaces = m_data.pfaces(i);
        for (const auto pface : pfaces) {
          wrapper.pface = pface;
          wrapper.support_plane = i;
          CGAL_assertion(pface_neighbors.find(pface) != pface_neighbors.end());
          const auto& pair = pface_neighbors.at(pface);
          wrapper.neighbors = pair;
          wrappers.push_back(wrapper);
        }
      }
      CGAL_assertion(wrappers.size() > 6);
    }

    void compute_weights(
      std::vector<Wrapper>& wrappers, std::size_t inliers) {

      FT min = 10000000000;
      FT max = 0;

      std::vector<PFace> inside, boundary;

      m_total_area = 0;
      for (auto& wrapper : wrappers) {
        auto& weight = wrapper.weight;
        const auto& pface = wrapper.pface;

        if (wrapper.neighbors.first < 0 || wrapper.neighbors.second < 0) {
          boundary.push_back(pface);
        }
        else inside.push_back(pface);

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

        m_total_area += weight;
        min = (std::min<FT>)(min, weight);
        max = (std::max<FT>)(max, weight);
      }

      std::cout << "total area: " << m_total_area << " min: " << min << " max: " << max << " mean: " << (m_total_area / wrappers.size()) << std::endl;

      dump_pfaces(m_data, inside, "inside_pfaces.ply");
      dump_pfaces(m_data, boundary, "boundary_pfaces.ply");
      std::cout << inside.size() << " inside faces" << std::endl;
      std::cout << boundary.size() << " boundary faces" << std::endl;

      CGAL_assertion(m_total_area > FT(0));
      for (auto& wrapper : wrappers) {
        wrapper.weight = 2.0 * inliers * wrapper.weight / m_total_area;
      }
    }

    void compute_weights(
      std::vector<Volume_cell>& volumes) const {

      FT sum = FT(0);
      const Converter converter;

      std::size_t index = 0;

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
        index++;
      }

      CGAL_assertion(sum > FT(0));
      for (auto& volume : volumes) {
        volume.weight /= sum;
      }

      std::cout << volumes.size() << " volumes" << std::endl;
    }

    void set_graph_edges(
      const std::vector<Wrapper>& wrappers,
      std::vector< std::pair<std::size_t, std::size_t> >& edges,
      std::vector<double>& edge_costs) const {

      std::map<std::pair<std::size_t, std::size_t>, std::size_t> edges2index;

      edges.clear();
      edge_costs.clear();
      std::size_t internal = 0, external = 0;
      for (const auto& wrapper : wrappers) {

        const FT edge_weight = wrapper.weight;
        const auto& neighbors = wrapper.neighbors;

        const int idx1 = neighbors.first;
        const int idx2 = neighbors.second;

        std::size_t id1, id2;

        bool intern = false;

        // Boundary edges.
        CGAL_assertion(idx1 >= 0 || idx2 >= 0);
        if( (idx1 < 0 && idx2 >= 0) || (idx2 < 0 && idx1 >= 0)) {
          if (idx1 < 0) {
            id1 = wrapper.support_plane;
            id2 = static_cast<std::size_t>(idx2) + 6;
          }
          else {
            id1 = static_cast<std::size_t>(idx1) + 6;
            id2 = wrapper.support_plane;
          }
        }
        else {
          intern = true;
          // Internal edges.
          CGAL_assertion(idx1 >= 0);
          id1 = static_cast<std::size_t>(idx1) + 6;
          CGAL_assertion(idx2 >= 0);
          id2 = static_cast<std::size_t>(idx2) + 6;
        }

        if (id2 < id1) {
          std::size_t tmp = id1;
          id1 = id2;
          id2 = tmp;
        }

        std::pair<std::size_t, std::size_t> idx(id1, id2);

        auto it = edges2index.find(idx);
        if (it == edges2index.end()) {
          edges2index[idx] = edges.size();
          edges.push_back(std::make_pair(id1, id2));
          edge_costs.push_back(compute_edge_cost(edge_weight));
          if (intern)
            internal++;
          else
            external++;
        }
        else {
          edge_costs[edges2index[idx]] += compute_edge_cost(edge_weight);
        }

      }
      std::cout << internal << " internal " << external << " external edges" << std::endl;
    }

    double compute_edge_cost(const FT edge_weight) const {

      //CGAL_assertion(m_beta      >= FT(0) && m_beta      <= FT(1));
      //CGAL_assertion(edge_weight >= FT(0) && edge_weight <= FT(1));
      return CGAL::to_double(m_beta * edge_weight);
    }

    void set_cost_matrix(
      const std::vector<Volume_cell>& volumes,
      std::vector< std::vector<double> >& cost_matrix) const {

      cost_matrix.clear();
      cost_matrix.resize(2);
      cost_matrix[0].resize(volumes.size() + 6);
      cost_matrix[1].resize(volumes.size() + 6);

      FT min_in = 100000000, max_in = 0, mean_in = 0;
      FT min_out = 100000000, max_out = 0, mean_out = 0;
      int min_in_count = 100000000, max_in_count = -100000000;
      int min_out_count = 100000000, max_out_count = -100000000;

      // Searching minimum values
      for (std::size_t i = 0; i < volumes.size(); i++) {
        const auto& volume = volumes[i];
        min_in_count = (std::min<int>)(min_in_count, static_cast<int>(volume.inside_count));
        max_in_count = (std::max<int>)(max_in_count, static_cast<int>(volume.inside_count));
        min_out_count = (std::min<int>)(min_out_count, static_cast<int>(volume.outside_count));
        max_out_count = (std::max<int>)(max_out_count, static_cast<int>(volume.outside_count));
      }

      std::cout << "min in count: " << min_in_count << " max in count: " << max_in_count << std::endl;
      std::cout << "min out count: " << min_out_count << " max out count: " << max_out_count << std::endl;

      int minimum = (min_in_count < min_out_count) ? min_in_count : min_out_count;

      // Setting preferred outside label for bbox plane nodes
      // Order:
      // 0 zmin
      // 1 ymin
      // 2 xmax
      // 3 ymax
      // 4 xmin
      // 5 zmax

      for (std::size_t i = 0; i < 6; i++) {
        cost_matrix[0][i] = 100000000;
      }
      cost_matrix[0][0] = -100000000;

      for (std::size_t i = 0; i < 6; i++) {
        cost_matrix[1][i] = -cost_matrix[0][i];
      }

      // Using new data term
      if (true) {
        for (std::size_t i = 0; i < volumes.size(); i++) {
          cost_matrix[0][i + 6] = (volumes[i].outside_count - volumes[i].inside_count) * (1.0 - m_beta);
          cost_matrix[1][i + 6] = (volumes[i].inside_count - volumes[i].outside_count) * (1.0 - m_beta);
          //std::cout << i << ": " << cost_matrix[0][i + 6] << " " << cost_matrix[1][i + 6] << std::endl;
          mean_in += cost_matrix[0][i + 6];
          mean_out += cost_matrix[1][i + 6];
        }
      }
      else {
        for (std::size_t i = 0; i < volumes.size(); ++i) {
          const auto& volume = volumes[i];

          const FT in = volume.inside;
          const FT out = volume.outside;

          CGAL_assertion(in >= FT(0) && in <= FT(1));
          CGAL_assertion(out >= FT(0) && out <= FT(1));
          CGAL_assertion((in + out) == FT(1));

          const FT face_weight = volume.weight;
          const double cost_in = get_face_cost(in, face_weight);
          const double cost_out = get_face_cost(out, face_weight);

          cost_matrix[0][i + 6] = cost_in;
          cost_matrix[1][i + 6] = cost_out;

          min_in = (std::min<FT>)(min_in, cost_in);
          max_in = (std::max<FT>)(max_in, cost_in);
          mean_in += cost_in;

          min_out = (std::min<FT>)(min_out, cost_out);
          max_out = (std::max<FT>)(max_out, cost_out);
          mean_out += cost_out;

          min_in_count = (std::min<int>)(min_in_count, volume.inside_count);
          max_in_count = (std::max<int>)(max_in_count, volume.inside_count);
          min_out_count = (std::min<int>)(min_out_count, volume.outside_count);
          max_out_count = (std::max<int>)(max_out_count, volume.outside_count);

          //std::cout << volume.index << " in: " << cost_in << " before: " << in << " " << volume.inside_count << std::endl;
          //std::cout << " out: " << cost_out << " before " << out << " " << volume.outside_count << std::endl;
        }
      }

      mean_in /= FT(volumes.size());
      mean_out /= FT(volumes.size());

      //std::cout << "min in: " << min_in << " max in: " << max_in << "  mean: " << mean_in << std::endl;
      //std::cout << "min out: " << min_out << " max out: " << max_out << "  mean: " << mean_out << std::endl;
    }

    double get_face_cost(
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
      labels.resize(volumes.size() + 6);

      for (std::size_t i = 0; i < 6; i++) {
        labels[i] = 1;
      }

      for (std::size_t i = 0; i < volumes.size(); ++i) {
        const auto& volume = volumes[i];
        if (volume.visibility == Visibility_label::INSIDE) {
          labels[i + 6] = 0;
        } else {
          labels[i + 6] = 1;
        }
      }
    }

    void compute_graphcut(
      const std::vector< std::pair<std::size_t, std::size_t> >& edges,
      const std::vector<double>& edge_costs,
      const std::vector< std::vector<double> >& cost_matrix,
      std::vector<std::size_t>& labels) const {
      std::vector<std::size_t> tmp = labels;

      std::cout << "beta" << m_beta << std::endl;

      double min = 10000000000;
      double max = -min;
      for (std::size_t i = 0; i < edge_costs.size(); i++) {
        min = (std::min<double>)(edge_costs[i], min);
        max = (std::max<double>)(edge_costs[i], max);
      }

      std::cout << "edge costs" << std::endl;
      std::cout << "min: " << min << std::endl;
      std::cout << "max: " << max << std::endl;

      min = 1000000000;
      max = -min;

      for (std::size_t i = 6; i < cost_matrix[0].size(); i++) {
        min = (std::min<double>)(cost_matrix[0][i], min);
        max = (std::max<double>)(cost_matrix[0][i], max);
      }

      std::cout << "label costs" << std::endl;
      std::cout << "min: " << min << std::endl;
      std::cout << "max: " << max << std::endl;

/*
      CGAL::min_cut(
          edges, edge_costs, cost_matrix, labels, CGAL::parameters::implementation_tag(CGAL::Alpha_expansion_MaxFlow_tag()));

      bool difference = false;
      for (std::size_t i = 0; i < labels.size(); i++) {
        if (tmp[i] != labels[i]) {
          difference = true;
          break;
        }
      }
      std::cout << "Labels changed: " << difference << std::endl;
      */

      CGAL::alpha_expansion_graphcut(
            edges, edge_costs, cost_matrix, labels, CGAL::parameters::implementation_tag(CGAL::Alpha_expansion_MaxFlow_tag()));

      bool difference = false;
      for (std::size_t i = 0; i < labels.size(); i++) {
        if (tmp[i] != labels[i]) {
          difference = true;
          break;
        }
      }

      std::cout << "Labels changed: " << difference << std::endl;
    }

    void apply_new_labels(
      const std::vector<std::size_t>& labels,
      std::vector<Volume_cell>& volumes) const {

      CGAL_assertion((volumes.size() + 6) == labels.size());
      for (std::size_t i = 0; i < volumes.size(); ++i) {
        const std::size_t label = labels[i + 6];
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
