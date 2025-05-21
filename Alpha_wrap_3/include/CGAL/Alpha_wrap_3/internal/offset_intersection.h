// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_OFFSET_INTERSECTION_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_OFFSET_INTERSECTION_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/number_utils.h>

#include <boost/algorithm/clamp.hpp>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

template <typename AABBTree,
          typename AABBTraversalTraits>
struct AABB_tree_oracle_helper;

template <typename AABBTree, typename AABBTraversalTraits>
struct AABB_distance_oracle
{
  using FT = typename AABBTree::FT;
  using Point_3 = typename AABBTree::Point;

  using AABB_helper = AABB_tree_oracle_helper<AABBTree, AABBTraversalTraits>;

  AABB_distance_oracle(const AABBTree& tree) : tree(tree) { }

  FT operator()(const Point_3& p) const
  {
    return approximate_sqrt(AABB_helper::squared_distance(p, tree));
  }

private:
  const AABBTree& tree;
};

// @todo even with EPECK, the precision cannot be 0 (otherwise it will not converge),
// thus exactness is pointless. Might as well use a cheap kernel (e.g. SC<double>), as long
// as there exists a mechanism to catch when the cheap kernel fails to converge (iterations?
// see also Tr_3::locate() or Mesh_3::Robust_intersection_traits_3.h)
template <class Kernel, class DistanceOracle>
class Offset_intersection
{
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;

public:
  Offset_intersection(const DistanceOracle& oracle,
                      const FT& off,
                      const FT& prec,
                      const FT& lip)
    : dist_oracle(oracle), offset(off), precision(prec), lipschitz(lip)
  { }

  bool first_intersection(const Point_3& s,
                          const Point_3& t,
                          Point_3& output_pt)
  {
    source = s;
    target = t;
    seg_length = approximate_sqrt(squared_distance(s, t));
    seg_unit_v = (t - s) / seg_length;
    const Point_2 p0 { 0, dist_oracle(source) };
    const Point_2 p1 { seg_length, dist_oracle(target) };

    return recursive_dichotomic_search(p0, p1, output_pt);
  }

private:
  Point_3 source;
  Point_3 target;
  FT seg_length;
  Vector_3 seg_unit_v;
  DistanceOracle dist_oracle;
  FT offset;
  FT precision;
  FT lipschitz;

  template <class Point>
  bool recursive_dichotomic_search(const Point_2& s, const Point_2& t,
                                   Point& output_pt)
  {
    if(CGAL::abs(s.x() - t.x()) < precision)
    {
      if(CGAL::abs(s.y() - offset) < precision)
      {
        const FT x_clamp = boost::algorithm::clamp<FT>(s.x(), FT{0}, seg_length);
        output_pt = source + (seg_unit_v * x_clamp);
        return true;
      }

      return false;
    }

    const bool sign_s = (s.y() > offset);
    const bool sign_t = (t.y() > offset);
    const FT gs_a = (sign_s) ? -lipschitz : lipschitz;
    const FT gs_b = s.y() - (gs_a * s.x());
    const FT gt_a = (sign_t) ? lipschitz : -lipschitz;

    const FT gt_b = t.y() - (gt_a * t.x());
    FT ms = (offset - gs_b) / gs_a;
    FT mt = (offset - gt_b) / gt_a;

    // early exit if there is no intersection
    if(sign_s == sign_t)
    {
      FT ui = (gt_b - gs_b) / (gs_a - gt_a);
      const FT gs_ui = (gs_a * ui) + gs_b;
      if((sign_s && (gs_ui > offset)) || (!sign_s && (gs_ui < offset)))
      {
        if(CGAL::abs(s.y() - offset) < precision)
        {
          const FT x_clamp = boost::algorithm::clamp<FT>(s.x(), FT{0}, seg_length);
          output_pt = source + (seg_unit_v * x_clamp);
          return true;
        }
        else if(CGAL::abs(t.y() - offset) < precision)
        {
          const FT x_clamp = boost::algorithm::clamp<FT>(t.x(), FT{0}, seg_length);
          output_pt = source + (seg_unit_v * x_clamp);
          return true;
        }

        return false;
      }
      else
      {
        ms = boost::algorithm::clamp<FT>(ms, FT{0}, seg_length);
        ui = boost::algorithm::clamp<FT>(ui, FT{0}, seg_length);
        mt = boost::algorithm::clamp<FT>(mt, FT{0}, seg_length);
        const Point_2 ms_pt { ms, dist_oracle(source + (seg_unit_v * ms)) };
        const Point_2 ui_pt { ui, dist_oracle(source + (seg_unit_v * ui)) };
        const Point_2 mt_pt { mt, dist_oracle(source + (seg_unit_v * mt)) };

        if(CGAL::abs(ms_pt.y() - offset) < precision)
        {
          const FT x_clamp = boost::algorithm::clamp<FT>(ms_pt.x(), FT{0}, seg_length);
          output_pt = source + (seg_unit_v * x_clamp);
          return true;
        }
        else if(CGAL::abs(ui_pt.y() - offset) < precision)
        {
          const FT x_clamp = boost::algorithm::clamp<FT>(ui_pt.x(), FT{0}, seg_length);
          output_pt = source + (seg_unit_v * x_clamp);
          return true;
        }
        else if(CGAL::abs(mt_pt.y() - offset) < precision)
        {
          const FT x_clamp = boost::algorithm::clamp<FT>(mt_pt.x(), FT{0}, seg_length);
          output_pt = source + (seg_unit_v * x_clamp);
          return true;
        }

        return (recursive_dichotomic_search(ms_pt, ui_pt, output_pt) ||
                recursive_dichotomic_search(ui_pt, mt_pt, output_pt));
      }
    }
    else // there is an intersection
    {
      if(CGAL::abs(mt - ms) <= precision) // linear approximation
      {
        const FT fsft_a = (t.y() - s.y()) / (t.x() - s.x());
        const FT fsft_b = s.y() - fsft_a * s.x();
        FT m_fsft;
        if(fsft_a == FT{0})
        {
          if(CGAL::abs(s.y() - offset) < precision)
            m_fsft = s.x();
          else
            return false;
        }
        else
        {
          m_fsft = (offset - fsft_b) / fsft_a;
        }
        m_fsft = boost::algorithm::clamp<FT>(m_fsft, FT{0}, seg_length);
        output_pt = source + (seg_unit_v * m_fsft);
        return true;
      }
      else
      {
        FT m = (ms + mt) / FT{2};
        ms = boost::algorithm::clamp<FT>(ms, FT{0}, seg_length);
        m = boost::algorithm::clamp<FT>(m, FT{0}, seg_length);
        mt = boost::algorithm::clamp<FT>(mt, FT{0}, seg_length);

        const Point_2 ms_pt { ms, dist_oracle(source + (seg_unit_v * ms)) };
        const Point_2 m_pt { m, dist_oracle(source + (seg_unit_v * m)) };
        const Point_2 mt_pt { mt, dist_oracle(source + (seg_unit_v * mt)) };

        return (recursive_dichotomic_search(ms_pt, m_pt, output_pt) ||
                recursive_dichotomic_search(m_pt, mt_pt, output_pt));
      }
    }
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_OFFSET_INTERSECTION_H
