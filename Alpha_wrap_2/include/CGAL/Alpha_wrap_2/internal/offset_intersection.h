// Copyright (c) 2019-2022 Google LLC (USA).
// Copyright (c) 2025 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_ALPHA_WRAP_2_INTERNAL_OFFSET_INTERSECTION_H
#define CGAL_ALPHA_WRAP_2_INTERNAL_OFFSET_INTERSECTION_H

#include <CGAL/license/Alpha_wrap_2.h>

#include <CGAL/number_utils.h>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

template <typename AABBTree,
          typename AABBTraversalTraits>
struct AABB_tree_oracle_helper;

template <typename AABBTree, typename AABBTraversalTraits>
struct AABB_distance_oracle
{
  using FT = typename AABBTree::FT;
  using Point_2 = typename AABBTree::Point;

  using AABB_helper = AABB_tree_oracle_helper<AABBTree, AABBTraversalTraits>;

  AABB_distance_oracle(const AABBTree& tree) : tree(tree) { }

  FT operator()(const Point_2& p) const
  {
    return approximate_sqrt(AABB_helper::squared_distance(p, tree));
  }

public:
  const AABBTree& tree;
};

// @todo even with EPECK, the precision cannot be 0 (otherwise it will not converge),
// thus exactness is pointless. Might as well use a cheap kernel (e.g. SC<double>), as long
// as there exists a mechanism to catch when the cheap kernel fails to converge (iterations?)
template <class Kernel, class DistanceOracle>
class Offset_intersection
{
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Vector_2 = typename Kernel::Vector_2;

public:
  Offset_intersection(const DistanceOracle& oracle,
                      const FT& off,
                      const FT& prec,
                      const FT& lip)
    : dist_oracle(oracle), offset(off), precision(prec), lipschitz(lip),
      sq_offset_minus_precision(CGAL::square(offset - precision)),
      sq_offset_plus_precision(CGAL::square(offset + precision))
  {
    CGAL_assertion(offset > precision);
  }

  bool first_intersection(const Point_2& s,
                          const Point_2& t,
                          Point_2& output_pt)
  {
    return SQsphere_marching_search(s, t, output_pt);
  }

private:
  Point_2 source;
  Point_2 target;
  FT seg_length;
  Vector_2 seg_unit_v;
  DistanceOracle dist_oracle;
  FT offset;
  FT precision;
  FT lipschitz;

  FT sq_offset_minus_precision;
  FT sq_offset_plus_precision;

  bool sphere_marching_search(const Point_2& s,
                              const Point_2& t,
                              Point_2& output_pt)
  {
#ifdef CGAL_AW2_DEBUG_SPHERE_MARCHING
    std::cout << "Sphere march between " << s << " and " << t << std::endl;
#endif

    const FT sq_seg_length = squared_distance(s, t);
    const FT seg_length = approximate_sqrt(sq_seg_length);
    const Vector_2 seg_unit_v = (t - s) / seg_length;

    Point_2 current_pt = s;
    Point_2 closest_point = dist_oracle.tree.closest_point(current_pt);
    FT current_dist = approximate_sqrt(squared_distance(current_pt, closest_point)) - offset;

    for(;;)
    {
#ifdef CGAL_AW2_DEBUG_SPHERE_MARCHING
      std::cout << "current point " << current_pt << std::endl;
      std::cout << "current dist " << current_dist << std::endl;
#endif

      if(CGAL::abs(current_dist) < precision)
      {
        output_pt = current_pt;
        return true;
      }

      // use the previous closest point as a hint: it's an upper bound
      current_pt = current_pt + (current_dist * seg_unit_v);

      if(squared_distance(s, current_pt) > sq_seg_length)
        return false;

      closest_point = dist_oracle.tree.closest_point(current_pt, closest_point /*hint*/);
      current_dist = approximate_sqrt(squared_distance(current_pt, closest_point)) - offset;
    }

    return false;
  }

  bool SQsphere_marching_search(const Point_2& s,
                                const Point_2& t,
                                Point_2& output_pt)
  {
#ifdef CGAL_AW2_DEBUG_SPHERE_MARCHING
    std::cout << "Sphere march between " << s << " and " << t << std::endl;
#endif

    const FT sq_seg_length = squared_distance(s, t);
    const FT seg_length = approximate_sqrt(sq_seg_length);
    const Vector_2 seg_unit_v = (t - s) / seg_length;

    Point_2 current_pt = s;
    Point_2 closest_point = dist_oracle.tree.closest_point(current_pt);
    FT sq_current_dist = squared_distance(current_pt, closest_point);
    FT step = 0;

#ifdef CGAL_AW2_DEBUG_SPHERE_MARCHING
    std::cout << "bounds: " << sq_offset_minus_precision << " " << sq_offset_plus_precision << std::endl;
#endif

    for(;;)
    {
#ifdef CGAL_AW2_DEBUG_SPHERE_MARCHING
      std::cout << "current point " << current_pt << std::endl;
      std::cout << "current sq dist " << sq_current_dist << std::endl;
      std::cout << "closest point: " << closest_point << std::endl;
      std::cout << "sq dist to closest: " << sq_current_dist << std::endl;
#endif

      // abs(dist - offset) < epsilon
      if((sq_current_dist > sq_offset_minus_precision) &&
         (sq_current_dist < sq_offset_plus_precision))
      {
        output_pt = current_pt;
        return true;
      }

      step += (std::max)(approximate_sqrt(sq_current_dist) - offset, 2 * precision);
      CGAL_assertion(step > 0);
      current_pt = s + (step * seg_unit_v);

      if(squared_distance(s, current_pt) > sq_seg_length)
      {
#ifdef CGAL_AW2_DEBUG_SPHERE_MARCHING
        std::cout << "Next target farther than the segment's extremity: " << current_pt << std::endl;
#endif
        return false;
      }

      // the previous closest point gives an upper bound so it's a good hint
      // @todo
      // closest_point = dist_oracle.tree.closest_point(current_pt, closest_point /*hint*/);
      closest_point = dist_oracle.tree.closest_point(current_pt);
      sq_current_dist = squared_distance(current_pt, closest_point);
    }

    return false;
  }

  bool SQsphere_marching_search_pp(const Point_2& s,
                                   const Point_2& t,
                                   Point_2& output_pt)
  {
#ifdef CGAL_AW2_DEBUG_SPHERE_MARCHING
    std::cout << "Sphere march between " << s << " and " << t << std::endl;
#endif

    const FT seg_length = approximate_sqrt(squared_distance(s, t));
    const Vector_2 seg_unit_v = (t - s) / seg_length;

    Point_2 current_pt = s;
    Point_2 closest_point = dist_oracle.tree.closest_point(current_pt);
    FT sq_current_dist = squared_distance(current_pt, closest_point);
    FT step = 0;

    bool relaxing = true;
    FT w = 1.8; // over-extending factor

    for(;;)
    {
#ifdef CGAL_AW2_DEBUG_SPHERE_MARCHING
      std::cout << "current point " << current_pt << std::endl;
      std::cout << "current sq dist " << sq_current_dist << std::endl;
      std::cout << "bounds: " << sq_offset_minus_precision << " " << sq_offset_plus_precision << std::endl;
#endif

      // If abs(dist - offset) < precision, we're done
      if((sq_current_dist > sq_offset_minus_precision) &&
         (sq_current_dist < sq_offset_plus_precision))
      {
        output_pt = current_pt;
        return true;
      }

      const Point_2 previous_pt = current_pt;
      const Point_2 previous_hint = closest_point;
      const FT previous_radius = approximate_sqrt(sq_current_dist) - offset;
      const FT previous_step = step;

      const FT local_step = (std::max)(previous_radius, 2 * precision);

      if(relaxing)
      {
        step += w * local_step;
        w = 1.1; // take bigger and bigger steps
      }
      else
      {
        step += local_step;
      }

      CGAL_assertion(step > 0);

      // move to the next point on the segment
      current_pt = s + (step * seg_unit_v);

      // the previous closest point gives an upper bound so it's a good hint
      closest_point = dist_oracle.tree.closest_point(current_pt, closest_point /*hint*/);
      sq_current_dist = squared_distance(current_pt, closest_point);

      // check if we have over-relaxed (the sphere are disjoint)
      if(relaxing)
      {
        const FT centers_dist = approximate_sqrt(squared_distance(previous_pt, current_pt));
        const FT current_radius = approximate_sqrt(sq_current_dist) - offset;
        if(previous_radius + current_radius < centers_dist)
        {
#ifdef CGAL_AW2_DEBUG_SPHERE_MARCHING
          std::cout << "Over relaxed, reverting" << std::endl;
#endif

          // revert the last step, and no more relaxation
          relaxing = false;

          current_pt = previous_pt;
          closest_point = previous_hint;
          sq_current_dist = squared_distance(current_pt, closest_point);
          step = previous_step;

          continue;
        }
      }

      // check if we ran out of segment to test
      if(step > seg_length)
        return false;
    }

    return false;
  }
};

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_INTERNAL_OFFSET_INTERSECTION_H
