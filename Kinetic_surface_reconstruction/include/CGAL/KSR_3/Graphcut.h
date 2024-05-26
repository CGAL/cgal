// Copyright (c) 2023 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Florent Lafarge, Dmitry Anisimov, Simon Giraudot

#ifndef CGAL_KSR_3_GRAPHCUT_H
#define CGAL_KSR_3_GRAPHCUT_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

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
#include <CGAL/KSP/utils.h>
#include <CGAL/KSP/debug.h>
#include <CGAL/KSP_3/Data_structure.h>

namespace CGAL {
namespace KSR_3 {

#ifdef DOXYGEN_RUNNING
#else

  template<typename GeomTraits>
  class Graphcut {

  public:
    using FT          = typename GeomTraits::FT;

    Graphcut(
      const FT beta) : m_beta(beta) { }

    void solve(const std::vector<std::pair<std::size_t, std::size_t> >& edges, const std::vector<FT>& edge_weights, const std::vector<std::vector<double> > &cost_matrix, std::vector<std::size_t> &labels) {
      assert(edges.size() == edge_weights.size());
      assert(!cost_matrix.empty());
      labels.resize(cost_matrix[0].size());
      for (std::size_t i = 0; i < cost_matrix[0].size(); i++) {
        // Verify quadratic size
        assert(cost_matrix[0].size() == cost_matrix[1].size());
        labels[i] = (cost_matrix[0][i] > cost_matrix[1][0]) ? 1 : 0;
      }
      compute_graphcut(edges, edge_weights, cost_matrix, labels);
    }

  private:
    const FT m_beta;

    double get_face_cost(
      const FT face_prob, const FT face_weight) const {

      CGAL_assertion(face_prob   >= FT(0) && face_prob   <= FT(1));
      CGAL_assertion(face_weight >= FT(0) && face_weight <= FT(1));

      const double weight = CGAL::to_double(face_weight);
      const double value  = (1.0 - CGAL::to_double(face_prob));
      return weight * value;
    }

    void compute_graphcut(
      const std::vector< std::pair<std::size_t, std::size_t> >& edges,
      const std::vector<FT>& edge_costs,
      const std::vector< std::vector<double> >& cost_matrix,
      std::vector<std::size_t>& labels, bool verbose = false) const {
      std::vector<std::size_t> tmp = labels;

      if (verbose)
        std::cout << std::endl << "beta " << m_beta << std::endl;

      std::vector<FT> edge_costs_lambda(edge_costs.size());
      std::vector<std::vector<double> > cost_matrix_lambda(2);

      for (std::size_t i = 0; i < edge_costs.size(); i++)
        edge_costs_lambda[i] = edge_costs[i] * m_beta;

      for (std::size_t i = 0; i < cost_matrix.size(); i++) {
        cost_matrix_lambda[i].resize(cost_matrix[i].size());
        for (std::size_t j = 0; j < cost_matrix[i].size(); j++)
          cost_matrix_lambda[i][j] = cost_matrix[i][j] * (1.0 - m_beta);
      }

      if (verbose) {
        double min = (std::numeric_limits<double>::max)();
        double max = -min;
        double edge_sum = 0;
        for (std::size_t i = 0; i < edge_costs.size(); i++) {
          min = (std::min<double>)(edge_costs[i], min);
          max = (std::max<double>)(edge_costs[i], max);
          edge_sum += edge_costs[i] * m_beta;
        }

        std::cout << "edge costs" << std::endl;
        std::cout << "sum: " << edge_sum << std::endl;
        std::cout << "min: " << min << std::endl;
        std::cout << "max: " << max << std::endl;

        min = (std::numeric_limits<double>::max)();
        max = -min;

        std::size_t sum_inside = 0;
        std::size_t sum_outside = 0;

        for (std::size_t i = 6; i < cost_matrix[0].size(); i++) {
          sum_inside += static_cast<std::size_t>(cost_matrix[0][i]);
          sum_outside += static_cast<std::size_t>(cost_matrix[1][i]);
          min = (std::min<double>)(cost_matrix[0][i], min);
          min = (std::min<double>)(cost_matrix[1][i], min);
          max = (std::max<double>)(cost_matrix[0][i], max);
          max = (std::max<double>)(cost_matrix[1][i], max);
        }

        std::cout << "label costs" << std::endl;
        std::cout << "min: " << min << std::endl;
        std::cout << "max: " << max << std::endl;
        std::cout << "sum inside: " << sum_inside << std::endl;
        std::cout << "sum outside: " << sum_outside << std::endl;
      }

      CGAL::alpha_expansion_graphcut(edges, edge_costs_lambda, cost_matrix_lambda, labels);
    }
  };

#endif //DOXYGEN_RUNNING

} // KSR_3
} // CGAL

#endif // CGAL_KSR_3_GRAPHCUT_H
