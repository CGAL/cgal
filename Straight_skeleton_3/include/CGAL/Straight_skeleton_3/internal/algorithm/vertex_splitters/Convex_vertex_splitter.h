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
 * @file   algo/3d/ConvexVertexSplitter.h
 * @author Gernot Walzl
 * @date   2013-02-01
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONVEX_VERTEX_SPLITTER_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_CONVEX_VERTEX_SPLITTER_H

#include <CGAL/Straight_skeleton_3/internal/algorithm/vertex_splitters/Combinatorial_vertex_splitter.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_self_intersection.h>

#include <memory>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class ConvexVertexSplitter
  : public CombiVertexSplitter<Traits>
{
  using Base = CombiVertexSplitter<Traits>;
  using ConvexVertexSplitterSPtr = std::shared_ptr<ConvexVertexSplitter>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;

private:
  using PolyhedronTransformation = algorithm::PolyhedronTransformation<Traits>;
  using SelfIntersection = algorithm::SelfIntersection<Traits>;

private:
  using vec2i = boost::shared_array<int>;
  using combi = std::vector<vec2i>;

public:
  ConvexVertexSplitter()
  {
    this->type_ = Base::CONVEX_VERTEX_SPLITTER;
    optimization_ = -1;
    ConfigurationSPtr config = Configuration::getInstance();
    if (config->isLoaded()) {
      std::string s_optimization = config->getString("algo_3d_ConvexVertexSplitter", "optimization");
      if (s_optimization.compare("max") == 0) {
        optimization_ = -1;
      } else if (s_optimization.compare("min") == 0) {
        optimization_ = 1;
      }
    }
  }

  virtual ~ConvexVertexSplitter() { /*intentionally does nothing*/ }

  static ConvexVertexSplitterSPtr create()
  {
    return std::make_shared<ConvexVertexSplitter>();
  }

  static int countConvexEdges(PolyhedronSPtr polyhedron)
  {
    int result = 0;
    typename std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
      EdgeSPtr edge = *it_e++;
      if (!edge->isReflex()) {
        result++;
      }
    }
    return result;
  }

  virtual PolyhedronSPtr splitVertex(VertexSPtr vertex)
  {
    std::cout << "\n> Splitting " << vertex->toString() << std::endl;

    PolyhedronSPtr polyhedron = vertex->getPolyhedron();
    if (vertex->degree() <= 3) {
        return polyhedron;
    }
    vertex->sort();
    std::list<combi> combinations = Base::generateAllCombinations(vertex->degree());
    std::cout << combinations.size() << " combinations" << std::endl;

    combi combi_opt;
    PolyhedronSPtr poly_opt;
    PolyhedronSPtr poly_opt_offset;
    int num_convex_edges = 0;
    int num_convex_edges_opt = 0;

    std::list<combi>::iterator it_combi = combinations.begin();
    while (it_combi != combinations.end()) {
      combi combination = *it_combi++;
      // std::cout << "-- Testing split-combination: " << Base::combiToString(combination) << std::endl;

      // don't take it out of the loop
      PolyhedronSPtr poly_c = Base::copyVertex(vertex);
      CGAL_assertion(bool(poly_c));
      CGAL_assertion(poly_c->facets().size() == vertex->facets().size());
      poly_c->initializeAllIDs();

      VertexSPtr vertex_c = poly_c->vertices().front();
      Base::splitVertex(vertex_c, combination);

      PolyhedronSPtr poly_c_offset = PolyhedronTransformation::shiftFacets(poly_c, -1.0);
      if (!poly_c_offset) {
        std::cerr << "Warning: failed to create offset of corner" << std::endl;
        continue;
      }

      // std::cout << "= Base Polyhedron\n" << poly_c->toString() << std::endl;
      // std::cout << "= Shifted Polyhedron\n" << poly_c_offset->toString() << std::endl;

      // static int test_id = -1;
      // ++test_id;
      // IO::OBJFile::save("results/split_" + std::to_string(test_id) + ".obj", poly_c, false);
      // IO::OBJFile::save("results/split_" + std::to_string(test_id) + "_offset.obj", poly_c_offset, false);

      // poly_c->dumpEdges("results/last_convex_split_base");
      // poly_c_offset->dumpEdges("results/last_convex_split_offset");

      num_convex_edges = countConvexEdges(poly_c_offset);
      auto updateOptimalCombination = [&](const combi& combination,
                                          PolyhedronSPtr& poly_c,
                                          PolyhedronSPtr& poly_c_offset,
                                          int num_convex_edges)
      {
        if (!SelfIntersection::hasSelfIntersectingSurface(poly_c_offset)) {
          CGAL_SS3_SPLITTER_TRACE("Valid split-combination found: " << Base::combiToString(combination));
          combi_opt = combination;
          poly_opt = poly_c;
          poly_opt_offset = poly_c_offset;
          num_convex_edges_opt = num_convex_edges;
        }
      };

      if (!poly_opt) {
        updateOptimalCombination(combination, poly_c, poly_c_offset, num_convex_edges);
        continue;
      }

      if (optimization_ < 0) {
        if (num_convex_edges > num_convex_edges_opt) {
          updateOptimalCombination(combination, poly_c, poly_c_offset, num_convex_edges);
        }
      } else if (optimization_ > 0) {
        if (num_convex_edges < num_convex_edges_opt) {
          updateOptimalCombination(combination, poly_c, poly_c_offset, num_convex_edges);
        }
      }
    }
    CGAL_assertion(combi_opt != combi());
    CGAL_SS3_SPLITTER_TRACE("Selected split-combination: " << Base::combiToString(combi_opt));
    CombiVertexSplitter<Traits>::apply(poly_opt, vertex);
    return polyhedron;
  }

  virtual std::string toString() const
  {
    std::string result("ConvexVertexSplitter(");
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

