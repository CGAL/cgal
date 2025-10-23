// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   algo/3d/ConvexVertexSplitter.h
 * author Gernot Walzl
 * date   2013-02-01
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONVEX_VERTEX_SPLITTER_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONVEX_VERTEX_SPLITTER_H

#include <CGAL/Straight_skeleton_3/internal/algorithm/vertex_splitters/Combinatorial_vertex_splitter.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_self_intersection.h>

#include <list>
#include <memory>
#include <string>
#include <vector>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class Convex_vertex_splitter
  : public Combi_vertex_splitter<Traits>
{
  using Base = Combi_vertex_splitter<Traits>;
  using Convex_vertex_splitter_sptr = std::shared_ptr<Convex_vertex_splitter>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;

private:
  using Transformation = algorithm::Polyhedron_transformation<Traits>;
  using Self_intersection = algorithm::Self_intersection<Traits>;

private:
  using vec2i = boost::shared_array<int>;
  using combi = std::vector<vec2i>;

public:
  Convex_vertex_splitter()
  {
    this->type_ = Base::CONVEX_VERTEX_SPLITTER;
    optimization_ = -1;
    ConfigurationSPtr config = Configuration::get_instance();
    if (config->is_loaded()) {
      std::string s_optimization = config->get_string("Algorithm", "convex_split_optimization");
      if (s_optimization.compare("max") == 0) {
        optimization_ = -1;
      } else if (s_optimization.compare("min") == 0) {
        optimization_ = 1;
      }
    }
  }

  virtual ~Convex_vertex_splitter() { /*intentionally does nothing*/ }

  static Convex_vertex_splitter_sptr create()
  {
    return std::make_shared<Convex_vertex_splitter>();
  }

  static int count_convex_edges(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    int result = 0;
    for (EdgeSPtr edge : polyhedron->edges()) {
      if (!edge->is_reflex()) {
        result++;
      }
    }
    return result;
  }

  virtual PolyhedronSPtr split_vertex(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_SS3_SPLITTER_TRACE_V(16, "\n> Splitting " << vertex->to_string());
    PolyhedronSPtr polyhedron = vertex->get_polyhedron();
    if (vertex->degree() <= 3) {
        return polyhedron;
    }
    vertex->sort();
    std::list<combi> combinations = Base::generate_all_combinations(vertex->degree());
    CGAL_SS3_SPLITTER_TRACE_V(16, combinations.size() << " combinations");

    combi combi_opt;
    PolyhedronSPtr poly_opt;
    PolyhedronSPtr poly_opt_offset;
    int num_convex_edges_opt = 0;

    std::list<combi>::iterator it_combi = combinations.begin();
    while (it_combi != combinations.end()) {
      combi combination = *it_combi++;
      CGAL_SS3_SPLITTER_TRACE_V(64, "-- Testing split-combination: " << Base::combi_to_string(combination));

      // don't take it out of the loop
      PolyhedronSPtr poly_c = Base::copy_vertex(vertex);
      CGAL_SS3_DEBUG_SPTR(poly_c);
      CGAL_assertion(poly_c->facets().size() == vertex->facets().size());
      poly_c->initialize_all_IDs();

      VertexSPtr vertex_c = poly_c->vertices().front();
      Base::split_vertex(vertex_c, combination);

      PolyhedronSPtr poly_c_offset = Base::copy_vertex(vertex);
      bool okShift = Transformation::shift_facets_deg1(poly_c_offset, -1.0);
      if (!okShift) {
        CGAL_SS3_SPLITTER_TRACE("Warning: failed to create offset of corner");
        continue;
      }

      auto update_optimal_combination = [&](const combi& combination,
                                            const PolyhedronSPtr& poly_c,
                                            const PolyhedronSPtr& poly_c_offset,
                                            int num_convex_edges)
      {
        if (!Self_intersection::has_self_intersecting_surface(poly_c_offset)) {
          CGAL_SS3_SPLITTER_TRACE("Valid split-combination found: " << Base::combi_to_string(combination));
          combi_opt = combination;
          poly_opt = poly_c;
          poly_opt_offset = poly_c_offset;
          num_convex_edges_opt = num_convex_edges;
        }
      };

      int num_convex_edges = count_convex_edges(poly_c_offset);

      if (!poly_opt) {
        update_optimal_combination(combination, poly_c, poly_c_offset, num_convex_edges);
        continue;
      }

      if (optimization_ < 0) {
        if (num_convex_edges > num_convex_edges_opt) {
          update_optimal_combination(combination, poly_c, poly_c_offset, num_convex_edges);
        }
      } else if (optimization_ > 0) {
        if (num_convex_edges < num_convex_edges_opt) {
          update_optimal_combination(combination, poly_c, poly_c_offset, num_convex_edges);
        }
      }
    }
    CGAL_assertion(combi_opt != combi());
    CGAL_SS3_SPLITTER_TRACE("Selected split-combination: " << Base::combi_to_string(combi_opt));
    Combi_vertex_splitter<Traits>::apply(poly_opt, vertex);
    return polyhedron;
  }

  virtual std::string to_string() const
  {
    std::string result("Convex_vertex_splitter(");
    if (optimization_ == -1) {
      result += "max";
    } else if (optimization_ == 1) {
      result += "min";
    }
    result += ")";
    return result;
  }

protected:
  int optimization_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONVEX_VERTEX_SPLITTER_H */

