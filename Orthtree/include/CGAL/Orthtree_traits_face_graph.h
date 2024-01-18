// Copyright (c) 2023  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien LORIOT


#ifndef CGAL_ORTHREE_TRAITS_FACE_GRAPH_H
#define CGAL_ORTHREE_TRAITS_FACE_GRAPH_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Orthtree_traits_base_for_dimension.h>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/Dimension.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {

/*!
\ingroup PkgOrthtreeTraits

Traits class for the `Orthtree` class to be used to contruct a 3D octree around
a triangulated surface mesh. Each node of the octree will store all the faces of the
mesh intersected by its bounding box. The subdivision of the octree is controlled
by the nested class `Orthtree_traits_face_graph::Split_predicate_node_min_extent`
to which the minimal extend of a node should be provided.

\tparam TriangleMesh a model of `FaceListGraph` with all faces being triangles
\tparam VertexPointMap a property map associating points to the vertices of `TriangleMesh`

\todo check how to adapt to non regular splits (cubes vs rectangular cuboid)

\cgalModels{OrthtreeTraits}
\sa `CGAL::Orthtree_traits_base_for_dimension<GeomTraits, DimensionTag>`
*/
template <class TriangleMesh, class VertexPointMap>
struct Orthtree_traits_face_graph : public Orthtree_traits_base_for_dimension<
  typename Kernel_traits<typename boost::property_traits<VertexPointMap>::value_type>::type,
  Dimension_tag<3> > {

  Orthtree_traits_face_graph(const TriangleMesh& pm, VertexPointMap vpm)
    : m_pm(pm), m_vpm(vpm) {}

  /// \name Types
  /// @{

  using Self = Orthtree_traits_face_graph<TriangleMesh, VertexPointMap>;
  using Tree = Orthtree<Self>;

  using Point_d = typename Self::Point_d;
  using Dimension = typename Self::Dimension;
  using Bbox_d = typename Self::Bbox_d;
  using FT = typename Self::FT;
  using Cartesian_const_iterator_d = typename Self::Cartesian_const_iterator_d;

  using Node_data = std::vector<typename boost::graph_traits<TriangleMesh>::face_descriptor>;

  using Geom_traits = typename Kernel_traits<Point_d>::type;

  using Construct_root_node_bbox = std::function<Bbox_d()>;
  using Construct_root_node_contents = std::function<Node_data()>;
  using Distribute_node_contents = std::function<void(typename Tree::Node_index, Tree&, const Point_d&)>;

  /// @}

  /// \name Operations
  /// @{

  auto construct_root_node_bbox_object() const {
    return [&]() -> Bbox_d {

      std::array<FT, Dimension::value> min = {0.0, 0}, max = {0.0, 0};
      if (faces(m_pm).begin() != faces(m_pm).end()) {
        const Point_d& p = get(m_vpm, *vertices(m_pm).begin());
        min = {p.x(), p.y(), p.z()};
        max = min;
        for (auto v: vertices(m_pm)) {
          const Point_d& p_v = get(m_vpm, v);
          for (int i = 0; i < 3; ++i) {
            if (p_v[i] < min[i]) min[i] = p_v[i];
            if (p_v[i] > max[i]) max[i] = p_v[i];
          }
        }
      }

      return {std::apply(Self::construct_point_d_object(), min),
              std::apply(Self::construct_point_d_object(), max)};
    };
  }

  auto construct_root_node_contents_object() {
    return [&]() -> Node_data {
      return {faces(m_pm).begin(), faces(m_pm).end()};
    };
  }

  auto distribute_node_contents_object() {
    return [&](typename Tree::Node_index n, Tree& tree, const Point_d& center) -> void {
      Node_data& ndata = tree.data(n);
      auto traits = tree.traits();
      for (int i = 0; i < 8; ++i) {
        typename Tree::Node_index child = tree.child(n, i);
        Node_data& child_data = tree.data(child);
        Bbox_d bbox = tree.bbox(child);
        for (auto f : ndata) {
          typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
            h = halfedge(f, m_pm);
          typename Geom_traits::Triangle_3 t(get(m_vpm, source(h, m_pm)),
            get(m_vpm, target(h, m_pm)),
            get(m_vpm, target(next(h, m_pm), m_pm)));
          if (do_intersect(t, bbox))
            child_data.push_back(f);
        }
      }
      };
  }

  /// @}

  /// Recommended split predicate to pass to `Orthtree::refine()` function so
  /// that the octree is refined until a node is either empty or has an extent
  /// that would be smaller after split than the value provided to the constructor.
  class Split_predicate_node_min_extent {

    FT m_min_extent;

  public:

    /// constructor with `me` being the minimal value a node extent could be.
    Split_predicate_node_min_extent(FT me)
      : m_min_extent(me) {}

    /*!
      \brief returns `true` if `ni` should be split, `false` otherwise.
     */
    template <typename Node_index, typename Tree>
    bool operator()(Node_index ni, const Tree& tree) const {
      if (tree.data(ni).empty()) return false;

      Bbox_d bb = tree.bbox(ni);

      for (int i = 0; i < 3; ++i)
        if ((bb.max()[i] - bb.min()[i]) < 2 * m_min_extent)
          return false;
      return true;
    }
  };


private:

  const TriangleMesh& m_pm;
  VertexPointMap m_vpm;
};

} // end of CGAL namespace


#endif // CGAL_ORTHREE_TRAITS_FACE_GRAPH_H