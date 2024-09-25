// Copyright (c) 2023  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sven Oesau


#ifndef CGAL_ORTHREE_TRAITS_POLYGONS_H
#define CGAL_ORTHREE_TRAITS_POLYGONS_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Orthtree_traits_base.h>

namespace CGAL
{

template <class GeomTraits>
struct Orthtree_traits_polygons : public Orthtree_traits_base<GeomTraits, 3>
{
  Orthtree_traits_polygons(const std::vector<typename GeomTraits::Point_3>& points, const std::vector<std::vector<std::size_t> >& polygons, typename GeomTraits::FT bbox_dilation = 1.1)
    : m_points(points), bbox_dilation(bbox_dilation) {
    m_polygons.resize(polygons.size());
    for (std::size_t i = 0;i<polygons.size();i++) {
      m_polygons[i].first = i;
      m_polygons[i].second.resize(polygons[i].size());
      for (std::size_t j = 0; j < polygons[i].size(); j++)
        m_polygons[i].second[j] = points[polygons[i][j]];
    }
  }

  using Self = Orthtree_traits_polygons<GeomTraits>;
  using Tree = Orthtree<Self>;

  using Geom_traits = GeomTraits;
  using Point_d = typename GeomTraits::Point_3;
  using Dimension = Dimension_tag<3>;
  using Bbox_d = typename GeomTraits::Iso_cuboid_3;
  using FT = typename Geom_traits::FT;
  using Sphere_d = typename Geom_traits::Sphere_3;
  using Array = std::array<FT, 3>;
  using Cartesian_const_iterator_d = typename Geom_traits::Cartesian_const_iterator_3;

  using Node_data_element = std::vector<Point_d>;
  using Node_data = std::vector<std::pair<std::size_t, Node_data_element> >;

  struct Construct_bbox_d {
    Bbox_d operator()(const Array& min,
                      const Array& max) const {
      return Bbox_d(min[0], min[1], min[2], max[0], max[1], max[2]);
    }
  };

  struct Construct_point_d_from_array {
    Point_d operator()(const Array& array) const {
      return Point_d(array[0], array[1], array[2]);
    }
  };

  Construct_point_d_from_array construct_point_d_from_array_object() const { return Construct_point_d_from_array(); }
  Construct_bbox_d construct_bbox_d_object() const { return Construct_bbox_d(); }

  auto construct_root_node_bbox_object() const {
    return [&]() -> typename Self::Bbox_d {
    Array bbox_min = {(std::numeric_limits<double>::max)(), (std::numeric_limits<double>::max)(), (std::numeric_limits<double>::max)() }, bbox_max = { -(std::numeric_limits<double>::max)(), -(std::numeric_limits<double>::max)(), -(std::numeric_limits<double>::max)() };

    for (const std::pair<std::size_t, Node_data_element> &p : m_polygons) {
      const Node_data_element& poly = p.second;
      for (int i = 0; i < static_cast<int>(poly.size()); i++)
        for (int d = 0; d < Dimension::value; d++) {
          bbox_min[d] = (std::min)(bbox_min[d], poly[i][d]);
          bbox_max[d] = (std::max)(bbox_max[d], poly[i][d]);
        }
    }

    for (std::size_t d = 0; d < Dimension::value; d++) {
      const FT mid = (bbox_min[d] + bbox_max[d]) * 0.5;
      const FT side = (bbox_max[d] - mid) * bbox_dilation;
      bbox_min[d] = mid - side;
      bbox_max[d] = mid + side;
    }

    return { construct_point_d_from_array_object()(bbox_min),
            construct_point_d_from_array_object()(bbox_max) };
    };
  }

  Point_d interpolate(FT a, FT b, FT l, const Point_d pa, const Point_d pb) const {
    FT f = CGAL::abs((a - l) / (a - b));
    assert(f <= 1.0);
    return Point_d((1 - f) * pa.x() + f * pb.x(), (1 - f) * pa.y() + f * pb.y(), (1 - f) * pa.z() + f * pb.z());
  }

  void clip_polygons_plane(std::size_t dimension, FT mid, const Node_data& polys, Node_data& lower_polys, Node_data& upper_polys) const {
    if (polys.empty())
      return;

    for (const std::pair<std::size_t, Node_data_element> & p : polys) {
      Node_data_element lower, upper;
      const Node_data_element& poly = p.second;

      FT last = poly.back()[static_cast<int>(dimension)];
      bool last_lower = (last <= mid);
      bool last_upper = (mid <= last);
      std::size_t last_index = poly.size() - 1;

      for (std::size_t i = 0; i < poly.size(); i++) {
        FT d = poly[i][static_cast<int>(dimension)];
        bool clower = d <= mid;
        bool cupper = mid <= d;

        const Point_d& l = poly[last_index];
        const Point_d& c = poly[i];

        if (last_lower && clower)
          lower.push_back(poly[i]);

        if (last_upper && cupper)
          upper.push_back(poly[i]);

        // switched sides?
        if (last_upper && !cupper)
          upper.push_back(interpolate(d, last, mid, c, l));

        if (last_lower && !clower)
          lower.push_back(interpolate(d, last, mid, c, l));

        if (!last_upper && cupper) {
          upper.push_back(interpolate(d, last, mid, c, l));
          upper.push_back(poly[i]);
        }

        if (!last_lower && clower) {
          lower.push_back(interpolate(d, last, mid, c, l));
          lower.push_back(poly[i]);
        }

        last = d;
        last_index = i;
        last_upper = cupper;
        last_lower = clower;
      }
      if (!upper.empty())
        upper_polys.emplace_back(std::make_pair(p.first, upper));

      if (!lower.empty())
        lower_polys.emplace_back(std::make_pair(p.first, lower));
    }
  }

  template<class NodeIndex, class Tree>
  void distribute_node_contents(NodeIndex n, Tree &tree, const Point_d &center) const {
    Node_data& ndata = tree.data(n);

    Node_data Xsplits[2];
    Node_data XYsplits[4];

    clip_polygons_plane(0, center.x(), ndata, Xsplits[0], Xsplits[1]);

    clip_polygons_plane(1, center.y(), Xsplits[0], XYsplits[0], XYsplits[1]);
    clip_polygons_plane(1, center.y(), Xsplits[1], XYsplits[2], XYsplits[3]);

    clip_polygons_plane(2, center.z(), XYsplits[0], tree.data(tree.child(n, 0)), tree.data(tree.child(n, 4)));
    clip_polygons_plane(2, center.z(), XYsplits[1], tree.data(tree.child(n, 2)), tree.data(tree.child(n, 6)));
    clip_polygons_plane(2, center.z(), XYsplits[2], tree.data(tree.child(n, 1)), tree.data(tree.child(n, 5)));
    clip_polygons_plane(2, center.z(), XYsplits[3], tree.data(tree.child(n, 3)), tree.data(tree.child(n, 7)));
  }

  auto construct_root_node_contents_object() const {
    return [&]() -> typename Self::Node_data {
      return m_polygons;
    };
  }

  auto distribute_node_contents_object() {
    return [&](typename Tree::Node_index n, Tree& tree, const typename Self::Point_d& center) {
      CGAL_precondition(!tree.is_leaf(n));
      distribute_node_contents(n, tree, center);
    };
  }

  Node_data m_polygons;
  const std::vector<Point_d>& m_points;
  FT bbox_dilation;
};

} // end of CGAL namespace


#endif // CGAL_ORTHREE_TRAITS_FACE_GRAPH_H
