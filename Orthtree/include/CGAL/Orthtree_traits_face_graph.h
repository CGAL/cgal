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

template <class PolygonMesh, class VPM>
struct Orthtree_traits_face_graph : public Orthtree_traits_base_for_dimension<
  typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type,
  Dimension_tag<3>
  // todo: it should be possible to determine the ambient dimension automatically, but this isn't working
//  Ambient_dimension<
//    typename boost::property_traits<VPM>::value_type,
//    typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type
//  >
> {

  Orthtree_traits_face_graph(const PolygonMesh& pm, VPM vpm)
    : m_pm(pm), m_vpm(vpm) {}

  using Self = Orthtree_traits_face_graph<PolygonMesh, VPM>;
  using Tree = Orthtree<Self>;

  using Point_d = typename Self::Point_d;
  using Dimension = typename Self::Dimension;
  using Bbox_d = typename Self::Bbox_d;
  using FT = typename Self::FT;
  using Cartesian_const_iterator_d = typename Self::Cartesian_const_iterator_d;

  // SL: these could be considered as built-in data and if the typedefs are not present, the tree have none
  using Node_data = std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor>;

  using Geom_traits = typename Kernel_traits<Point_d>::type;

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
      for (int i = 0; i < 8; ++i) {
        typename Tree::Node_index child = tree.child(n, i);
        Node_data& child_data = tree.data(child);
        Bbox_d bbox = tree.bbox(child);
        for (auto f: ndata) {
          typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
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

  class Split_predicate_node_min_extent {

    FT m_min_extent;

  public:

    Split_predicate_node_min_extent(FT me)
      : m_min_extent(me) {}

    /*!
      \brief returns `true` if `ni` should be split, `false` otherwise.
     */
    template <typename Node_index, typename Tree>
    bool operator()(Node_index ni, const Tree& tree) const {
      if (tree.data(ni).empty()) return false;

      Bbox_d bb = tree.bbox(ni);
      // TODO: we should get better version to get guarantees
      // TODO: as long as the bbox is cubic you can use depth and initial size to conclude.
      // todo (jackson): bbox is _not_ guaranteed to be cubic now, this may break in some very niche cases
      for (int i = 0; i < 3; ++i)
        if ((bb.max()[i] - bb.min()[i]) < 2 * m_min_extent)
          return false;
      return true;
    }
  };

  const PolygonMesh& m_pm;
  VPM m_vpm;
};

} // end of CGAL namespace


#endif // CGAL_ORTHREE_TRAITS_FACE_GRAPH_H
