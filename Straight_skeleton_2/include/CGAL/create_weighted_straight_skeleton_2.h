// Copyright(c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL(www.cgal.org).
//

// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_CREATE_WEIGHTED_STRAIGHT_SKELETON_2_H
#define CGAL_CREATE_WEIGHTED_STRAIGHT_SKELETON_2_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Straight_skeleton_2/Polygon_iterators.h>
#include <CGAL/create_straight_skeleton_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel_traits.h>

#include <optional>
#include <memory>
#include <algorithm>
#include <iterator>
#include <vector>

namespace CGAL {

template <typename PointIterator, typename HoleIterator,
          typename WeightIterator, typename HoleWeightsIterator,
          typename K>
std::shared_ptr<Straight_skeleton_2<K> >
create_interior_weighted_straight_skeleton_2(PointIterator outer_contour_vertices_begin,
                                             PointIterator outer_contour_vertices_end,
                                             HoleIterator holes_begin,
                                             HoleIterator holes_end,
                                             WeightIterator outer_contour_weights_begin,
                                             WeightIterator outer_contour_weights_end,
                                             HoleWeightsIterator holes_weights_begin,
                                             HoleWeightsIterator holes_weights_end,
                                             const K&)
{
  using Skeleton = Straight_skeleton_2<K>;

  using Skeleton_build_traits = Straight_skeleton_builder_traits_2<K>;

  using Skeleton_builder = Straight_skeleton_builder_2<Skeleton_build_traits, Skeleton>;

  using Input_point = typename std::iterator_traits<PointIterator>::value_type;
  using Input_kernel = typename Kernel_traits<Input_point>::Kernel;

  CGAL_precondition(std::distance(outer_contour_vertices_begin, outer_contour_vertices_end) ==
                    std::distance(outer_contour_weights_begin, outer_contour_weights_end));
  CGAL_precondition(std::distance(holes_begin, holes_end) ==
                    std::distance(holes_weights_begin, holes_weights_end));

  Cartesian_converter<Input_kernel, K> point_converter;
  NT_converter<typename Input_kernel::FT, typename K::FT> weight_converter;

  Skeleton_builder ssb;
  ssb.enter_contour(outer_contour_vertices_begin, outer_contour_vertices_end, point_converter);
  ssb.enter_contour_weights(outer_contour_weights_begin, outer_contour_weights_end, weight_converter);

  for(HoleIterator hi = holes_begin; hi != holes_end && holes_weights_begin != holes_weights_end; ++hi, ++holes_weights_begin)
  {
    CGAL_precondition(std::distance(CGAL_SS_i::vertices_begin(*hi), CGAL_SS_i::vertices_end(*hi)) ==
                      std::distance(holes_weights_begin->begin(), holes_weights_begin->end()));
    ssb.enter_contour(CGAL_SS_i::vertices_begin(*hi), CGAL_SS_i::vertices_end(*hi), point_converter);
    ssb.enter_contour_weights(holes_weights_begin->begin(), holes_weights_begin->end(), weight_converter);
  }

  return ssb.construct_skeleton();
}

template <typename PointIterator, typename HoleIterator,
          typename WeightIterator, typename HoleWeightsIterator,
          typename Weights>
std::shared_ptr<Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
inline
create_interior_weighted_straight_skeleton_2(PointIterator outer_contour_vertices_begin,
                                             PointIterator outer_contour_vertices_end,
                                             HoleIterator holes_begin,
                                             HoleIterator holes_end,
                                             WeightIterator outer_contour_weights_begin,
                                             WeightIterator outer_contour_weights_end,
                                             HoleWeightsIterator holes_weights_begin,
                                             HoleWeightsIterator holes_weights_end)
{
  return create_interior_weighted_straight_skeleton_2(outer_contour_vertices_begin, outer_contour_vertices_end,
                                                      holes_begin, holes_end,
                                                      outer_contour_weights_begin, outer_contour_weights_end,
                                                      holes_weights_begin, holes_weights_end,
                                                      Exact_predicates_inexact_constructions_kernel());
}

template <typename PointIterator,
          typename WeightIterator,
          typename K>
std::shared_ptr<Straight_skeleton_2<K> >
inline
create_interior_weighted_straight_skeleton_2(PointIterator outer_contour_vertices_begin,
                                             PointIterator outer_contour_vertices_end,
                                             WeightIterator outer_contour_weights_begin,
                                             WeightIterator outer_contour_weights_end,
                                             const K& k)
{
  using InputPoint = typename std::iterator_traits<PointIterator>::value_type;
  using InputKernel = typename Kernel_traits<InputPoint>::Kernel;
  using InputFT = typename InputKernel::FT;

  std::vector<Polygon_2<InputKernel> > no_holes;
  std::vector<std::vector<InputFT> > no_hole_weights;

  return create_interior_weighted_straight_skeleton_2(outer_contour_vertices_begin,
                                                      outer_contour_vertices_end,
                                                      no_holes.begin(),
                                                      no_holes.end(),
                                                      outer_contour_weights_begin,
                                                      outer_contour_weights_end,
                                                      no_hole_weights.begin(),
                                                      no_hole_weights.end(),
                                                      k);
}

template <typename PointIterator,
          typename WeightIterator>
std::shared_ptr<Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
inline
create_interior_weighted_straight_skeleton_2(PointIterator outer_contour_vertices_begin,
                                             PointIterator outer_contour_vertices_end,
                                             WeightIterator outer_contour_weights_begin,
                                             WeightIterator outer_contour_weights_end)
{
  return create_interior_weighted_straight_skeleton_2(outer_contour_vertices_begin,
                                                      outer_contour_vertices_end,
                                                      outer_contour_weights_begin,
                                                      outer_contour_weights_end,
                                                      Exact_predicates_inexact_constructions_kernel());
}

template <typename Polygon,
          typename Weights,
          typename K>
std::shared_ptr<Straight_skeleton_2<K> >
inline
create_interior_weighted_straight_skeleton_2(const Polygon& out_contour,
                                             const Weights& weights,
                                             const K& k,
                                             std::enable_if_t<! CGAL_SS_i::has_Hole_const_iterator<Polygon>::value>* = nullptr)
{
  return create_interior_weighted_straight_skeleton_2(CGAL_SS_i::vertices_begin(out_contour),
                                                      CGAL_SS_i::vertices_end(out_contour),
                                                      weights.begin(),
                                                      weights.end(),
                                                      k);
}

template <typename Polygon,
          typename Weights>
std::shared_ptr<Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
inline
create_interior_weighted_straight_skeleton_2(const Polygon& out_contour,
                                             const Weights& weights)
{
  return create_interior_weighted_straight_skeleton_2(out_contour,
                                                      weights,
                                                      Exact_predicates_inexact_constructions_kernel());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// EXTERIOR

template <typename FT,
          typename PointIterator,
          typename WeightIterator,
          typename K>
std::shared_ptr<Straight_skeleton_2<K> >
create_exterior_weighted_straight_skeleton_2(const FT& max_offset,
                                             PointIterator vertices_begin,
                                             PointIterator vertices_end,
                                             WeightIterator weights_begin,
                                             WeightIterator weights_end,
                                             const K& k)
{
  CGAL_precondition(max_offset > 0);
  CGAL_precondition(std::distance(weights_begin, weights_end) == std::distance(vertices_begin, vertices_end));

  using Point_2 = typename std::iterator_traits<PointIterator>::value_type;
  using IK = typename Kernel_traits<Point_2>::Kernel;
  using IFT = typename IK::FT;

  std::shared_ptr<Straight_skeleton_2<K> > skeleton;

  // That's because we might not have FT == IK::FT (e.g. `double` and `Core`)
  // Note that we can also have IK != K (e.g. `Simple_cartesian<Core>` and `EPICK`)
  IFT offset = max_offset;

  // @todo This likely should be done in the kernel K rather than the input kernel(i.e. the same
  // converter stuff that is done in `create_partial_exterior_straight_skeleton_2`?).
  std::optional<IFT> margin = compute_outer_frame_margin(vertices_begin,
                                                           vertices_end,
                                                           weights_begin,
                                                           weights_end,
                                                           offset);

  if(margin)
  {
    const IFT lm = *margin;
    const Bbox_2 bbox = bbox_2(vertices_begin, vertices_end);

    const IFT fxmin = IFT(bbox.xmin()) - lm;
    const IFT fxmax = IFT(bbox.xmax()) + lm;
    const IFT fymin = IFT(bbox.ymin()) - lm;
    const IFT fymax = IFT(bbox.ymax()) + lm;

    Point_2 frame[4];
    frame[0] = Point_2(fxmin,fymin);
    frame[1] = Point_2(fxmax,fymin);
    frame[2] = Point_2(fxmax,fymax);
    frame[3] = Point_2(fxmin,fymax);

    CGAL_STSKEL_BUILDER_TRACE(2, "Frame:\n" << frame[0] << "\n" << frame[1] << "\n" << frame[2] << "\n" << frame[3]);

    typedef std::vector<Point_2> Hole;

    Hole poly(vertices_begin, vertices_end);
    std::reverse(poly.begin(), poly.end());

    std::vector<Hole> holes;
    holes.push_back(poly);

    // put a weight large enough such that frame edges are not relevant
    const FT frame_weight = FT(10) * *(std::max_element(weights_begin, weights_end));
    CGAL_STSKEL_BUILDER_TRACE(4, "Frame weight = " << frame_weight);

    std::vector<FT> lFrameWeights(4, frame_weight);
    std::vector<std::vector<FT> > lHoleWeights;
    lHoleWeights.emplace_back(weights_begin, weights_end);

    // If w[0] pointed to v_0, then when we reverse the polygon, the last polygon is pointing to v_{n-1}
    // but it is the edge v_0 v_{n-1}, which has the weight w_0.
    std::reverse(lHoleWeights[0].begin(), lHoleWeights[0].end());
    std::rotate(lHoleWeights[0].rbegin(), lHoleWeights[0].rbegin()+1, lHoleWeights[0].rend());

    skeleton = create_interior_weighted_straight_skeleton_2(frame, frame+4,
                                                            holes.begin(), holes.end(),
                                                            lFrameWeights.begin(), lFrameWeights.end(),
                                                            lHoleWeights.begin(), lHoleWeights.end(),
                                                            k);
  }

  return skeleton;
}

template <typename FT,
          typename PointIterator,
          typename WeightIterator>
std::shared_ptr<Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
inline
create_exterior_weighted_straight_skeleton_2(const FT& max_offset,
                                             PointIterator vertices_begin,
                                             PointIterator vertices_end,
                                             WeightIterator weights_begin,
                                             WeightIterator weights_end)
{
  return create_exterior_weighted_straight_skeleton_2(max_offset,
                                                      vertices_begin,
                                                      vertices_end,
                                                      weights_begin,
                                                      weights_end,
                                                      Exact_predicates_inexact_constructions_kernel());
}

template <typename FT,
          typename Polygon,
          typename Weights,
          typename K>
std::shared_ptr<Straight_skeleton_2<K> >
inline
create_exterior_weighted_straight_skeleton_2(const FT& max_offset,
                                             const Polygon& aPoly,
                                             Weights& weights,
                                             const K& k)
{
  return create_exterior_weighted_straight_skeleton_2(max_offset,
                                                      CGAL_SS_i::vertices_begin(aPoly),
                                                      CGAL_SS_i::vertices_end(aPoly),
                                                      weights.begin(),
                                                      weights.end(),
                                                      k);
}

template <typename FT,
          typename Weights,
          typename Polygon>
std::shared_ptr<Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
inline
create_exterior_weighted_straight_skeleton_2(const FT& max_offset,
                                             Weights& weights,
                                             const Polygon& aPoly)
{
  return create_exterior_weighted_straight_skeleton_2(max_offset,
                                                      aPoly,
                                                      weights,
                                                      Exact_predicates_inexact_constructions_kernel());
}

} // namespace CGAL

#endif // CGAL_CREATE_WEIGHTED_STRAIGHT_SKELETON_2_H
