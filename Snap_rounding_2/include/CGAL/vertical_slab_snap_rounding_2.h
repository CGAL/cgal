// Copyright (c) 2025 Geometry Factory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Léo Valque

#ifndef CGAL_VERTICAL_SLAB_SNAP_ROUNDING_2_H
#define CGAL_VERTICAL_SLAB_SNAP_ROUNDING_2_H

#include <CGAL/license/Snap_rounding_2.h>

#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
#include <iostream>
#endif

#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Double_grid_snap_rounding_traits_2.h>
#include <CGAL/intersection_2.h>

#include <CGAL/internal/Wrap_vertical_slab_snap_rounding_traits_2.h>

#include <boost/property_map/function_property_map.hpp>

#include <CGAL/Polygon_2.h>
#include <CGAL/internal/Snap_rounding_helpers.h>

#include <set>
#include <vector>


namespace CGAL {

namespace internal {

template <typename GeometryTraits_2, typename Points_, typename Polylines_, typename Allocator_ = CGAL_ALLOCATOR(int)>
class Snap_rounding_visitor :
    public Ss2::Default_visitor<Snap_rounding_visitor<GeometryTraits_2, Points_, Polylines_, Allocator_>,
                                GeometryTraits_2, Allocator_> {
public:
  using Gt2 = GeometryTraits_2;
  using Points = Points_;
  using Polylines = Polylines_;
  using Allocator = Allocator_;

  using P_idx = std::size_t;
  using L_idx = std::size_t;
  using Polyline = typename Polylines::value_type;

private:
  using Self = Snap_rounding_visitor<Gt2, Points, Polylines, Allocator>;
  using Base = Ss2::Default_visitor<Self, Gt2, Allocator>;

public:
  using Event = typename Base::Event;
  using Subcurve = typename Base::Subcurve;
  using Status_line_iterator = typename Subcurve::Status_line_iterator;
  using X_monotone_curve_2 = typename Gt2::X_monotone_curve_2;
  using Point_2 = typename Gt2::Point_2;
  using Surface_sweep_2 = typename Base::Surface_sweep_2;
  using FT = typename Gt2::FT;

  using Base_segment_2 = typename Gt2::Base_segment_2;

protected:
  Gt2& traits;
  Points& points;
  const Polylines& input_polylines;
  Polylines& output_polylines;
  std::vector< double >& round_bounds;

  // If pi is closed enough to sc, subdivide sc in the output and return true
  bool is_pi_closed_to_sc_and_subdivide(P_idx &pi, Subcurve* sc){
    auto point_at_x = traits.construct_point_at_x_on_segment_2_object();
    auto csq_dist_2 = traits.compare_squared_distance_2_object();
    auto round_bound = traits.compute_squared_round_bound_2_object();

    L_idx l_idx = sc->last_curve().polyline_indices[0]; // Use to get an input segment to compute the squared distance
    P_idx src = sc->last_curve().src;
    P_idx trg = sc->last_curve().trg;
    // Skip if this polyline already ends with the event point
    if(!output_polylines[l_idx].empty() && output_polylines[l_idx].back() == pi)
      return true;

    const auto &p = points[pi];
    const auto &pl = input_polylines[l_idx];
    Base_segment_2 seg = Base_segment_2(points[pl.front()], points[pl.back()]); // Get an input segment supporting the current one

    // (A+B)^2 <= 4*max(A^2,B^2)
    double bound = (std::max)({ round_bounds[pi],
                                round_bounds[src],
                                round_bounds[trg] });
    bound *= 4;

    if(possibly(csq_dist_2(p, seg, bound) != LARGER)){
      if constexpr (typename Gt2::Evaluation_tag()){
        auto evaluate = traits.evaluate_object();
        // We refine the pts to reduce the rounding shift and check again
        evaluate(points[pi]);
        evaluate(points[output_polylines[l_idx].back()]);
        // The following calls of exact evaluation affect 'seg' since there are its endpoints
        evaluate(points[src]);
        evaluate(points[trg]);
        // Update the bounds
        round_bounds[pi]  = round_bound(points[pi]);
        round_bounds[output_polylines[l_idx].back()] = round_bound(points[output_polylines[l_idx].back()]);
        round_bounds[src] = round_bound(points[src]);
        round_bounds[trg] = round_bound(points[trg]);
        bound = (std::max)({ round_bounds[pi],
                             round_bounds[output_polylines[l_idx].back()],
                             round_bounds[src],
                             round_bounds[trg] });
        bound *= 4;

        // Check if the point and the segment are still too close for a safe rounding
        if(csq_dist_2(p, seg, bound) == LARGER)
          return false;
      }

      // Check if the segment was not already subdivided by another point
      if(points[output_polylines[l_idx].back()].x() == points[pi].x())
        return false;

      // Create a point on seg at the same x coordinate than p
      points.push_back(point_at_x(seg, points[pi].x()));
      round_bounds.push_back(round_bound(points.back()));
      P_idx new_pi=points.size()-1;
      pi = new_pi; // Store the new_points in pi so that it can be tested against neighboring segments

      for(auto l_idx: sc->last_curve().polyline_indices){
        CGAL_assertion(points[output_polylines[l_idx].back()].x() != points[pi].x());
#ifdef CGAL_DOUBLE_2D_SNAP_FULL_VERBOSE
        std::cout << "Create point " << new_pi << " on " << l_idx << " due to proximity with " << pi << "_____________________________" << std::endl;
#endif
        // We insert it on the output
        output_polylines[l_idx].push_back(new_pi);
      }
      return true;
    }
    return false;
  }

public:
  Snap_rounding_visitor(Gt2 &traits_, Points& points_, const Polylines& in, Polylines& out, std::vector< double >& bounds) :
    traits(traits_),
    points(points_),
    input_polylines(in),
    output_polylines(out),
    round_bounds(bounds)
  {
    output_polylines.resize(input_polylines.size());
  }

  bool after_handle_event(Event* event, Status_line_iterator iter, bool /* flag */) {
    if (! event->is_closed()) return true;
    if (! event->has_left_curves() && ! event->has_right_curves()) return true;

    auto pi = event->point();

    // Insert the point in the output
    for (auto it = event->left_curves_begin(); it != event->left_curves_end(); ++it) {
      const Subcurve *sc = *it;
      for(auto l_idx: sc->last_curve().polyline_indices)
        output_polylines[l_idx].push_back(pi);
    }

    for (auto it = event->right_curves_begin(); it != event->right_curves_end(); ++it) {
      const Subcurve *sc = *it;
      for(auto l_idx: sc->last_curve().polyline_indices)
        if(output_polylines[l_idx].empty() || output_polylines[l_idx].back() != pi) // Maybe already by left curves
          output_polylines[l_idx].push_back(pi);
    }

    // Check segments above and create a subdivision point if they are too close for safe rounding
    auto pi_above = pi;
    auto above = iter;
    while(above != this->status_line_end() && is_pi_closed_to_sc_and_subdivide(pi_above, *above))
      ++above;

    // same with segments below
    auto below = iter;
    if(below != this->status_line_begin()){
      --below;
      while(is_pi_closed_to_sc_and_subdivide(pi, *below)){
        if(below == this->status_line_begin())
          break;
        --below;
      }
    }

    return true;
  }
};

/*
scans the vertices from left to right while maintaining the status line ordering of the segments.
Segments that are too close to a vertex are subdivided.
*/
template <class Traits, class PointsRange , class PolylinesRange>
void snap_rounding_scan(PointsRange &pts, PolylinesRange &polylines, const Traits &traits){
  auto round_bound = traits.compute_squared_round_bound_2_object();
  auto less_y_2 = traits.less_y_2_object();

  // Wrap the traits so they operate on point indices stored in `pts`.
  using Base_point_2 = typename Traits::Point_2;
  struct Get_point
  {
    const std::vector<Base_point_2>& pts;
    Base_point_2 operator()(const std::size_t& idx) const { return pts[idx]; }
  };
  using PMap = boost::function_property_map<Get_point, std::size_t, Base_point_2>;
  using Wrap_traits = Wrap_vertical_slab_snap_rounding_traits_2<Traits, PMap>;

  Get_point get{pts};
  PMap pm(get);
  Wrap_traits wrap_traits(traits, pm);

  using X_monotone_curve_2 = typename Wrap_traits::X_monotone_curve_2;
  using Visitor = Snap_rounding_visitor<Wrap_traits, PointsRange, PolylinesRange>;
  using Surface_sweep = Ss2::No_intersection_surface_sweep_2<Visitor>;

  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename PolylinesRange::iterator>::value_type>;

  std::vector< double > round_bounds;
  round_bounds.reserve(pts.size());
  for(std::size_t i=0; i!=pts.size(); ++i)
    round_bounds.push_back(round_bound(pts[i]));

  // Create the curves
  std::set<X_monotone_curve_2> curves;
  for(std::size_t j=0; j!=polylines.size(); ++j)
    for(std::size_t i=0; i!=polylines[j].size()-1; ++i){
      std::size_t src = polylines[j][i];
      std::size_t trg = polylines[j][i+1];
      auto it = curves.emplace(X_monotone_curve_2(src, trg)).first;
      it->add_polyline(j);
    }
  std::vector< Polyline > out;
  Visitor visitor(wrap_traits, pts, polylines, out, round_bounds);
  Surface_sweep surface_sweep(&wrap_traits, &visitor);
  surface_sweep.sweep(curves.begin(), curves.end());

  // We sort the point by y, evaluation by this operation ensure the y-order of the point is preserved after the rounding
  // Round value of y coordinates are modified in the case of a filter failure
  auto sort_pi=[&](std::size_t i, std::size_t j){
    return less_y_2(pts[i],pts[j]);
  };
  std::vector< std::size_t > indices(pts.size(),0);
  std::iota(indices.begin(),indices.end(),0);
  std::sort(indices.begin(),indices.end(),sort_pi);

  std::swap(polylines, out);
  return;
}

template <class Traits, class PointsRange , class PolylinesRange>
void merge_duplicate_points_in_polylines(PointsRange &pts, PolylinesRange &polylines, const Traits &traits)
{
  using Point_2 = typename Traits::Point_2;
  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename PolylinesRange::iterator>::value_type>;

  using Less_xy_2 = typename Traits::Less_xy_2;
  using Equal_2 = typename Traits::Equal_2;
  Equal_2 equal = traits.equal_2_object();

  auto Less_indices_xy_2=[&](std::size_t i, std::size_t j){
    return Less_xy_2()(pts[i], pts[j]);
  };

  std::vector< std::size_t > unique_points(pts.size());
  std::iota(unique_points.begin(), unique_points.end(), 0);
  std::sort(unique_points.begin(), unique_points.end(), Less_indices_xy_2);
  std::vector<Point_2> new_pts;
  std::vector<std::size_t> old_to_new_index(pts.size());
  for(std::size_t i=0; i!=pts.size(); ++i){
    if(i==0 || !equal(pts[unique_points[i]],pts[unique_points[i-1]]))
      new_pts.push_back(pts[unique_points[i]]);
    old_to_new_index[unique_points[i]]=new_pts.size()-1;
  }

  std::swap(pts, new_pts);
  for (Polyline& polyline : polylines) {
    std::vector<std::size_t> updated_polyline;
    updated_polyline.reserve(polylines.size()); // Potentially more than is necessary
    for (std::size_t i=0; i<polyline.size(); ++i) {
      std::size_t new_pi=old_to_new_index[polyline[i]];
      if(i==0 || (new_pi!=updated_polyline[updated_polyline.size()-1]))
        updated_polyline.push_back(new_pi);
    }
    std::swap(polyline, updated_polyline);
  }
}

// Some points may have collapsed on a vertical segment, we subdivide these vertical segments accordingly
template <class Traits, class PointsRange , class PolylinesRange>
void snap_post_process(PointsRange &pts, PolylinesRange &polylines, const Traits &traits)
{
  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename PolylinesRange::iterator>::value_type>;
  using Less_xy_2 = typename Traits::Less_xy_2;

  Less_xy_2 less_xy_2 = traits.less_xy_2_object();

  auto Less_indices_xy_2=[&](std::size_t i, std::size_t j){
    return less_xy_2(pts[i], pts[j]);
  };

  std::vector< std::size_t > p_sort_by_x(pts.size());
  std::iota(p_sort_by_x.begin(), p_sort_by_x.end(), 0);
  std::sort(p_sort_by_x.begin(), p_sort_by_x.end(), Less_indices_xy_2);

  for(Polyline &poly: polylines){
    std::vector<std::size_t> updated_polyline;
    updated_polyline.push_back(poly.front());
    for(std::size_t i=1; i!=poly.size(); ++i){
      if(pts[poly[i-1]].x()==pts[poly[i]].x()){
        std::vector< std::size_t >::iterator start, end;
        // Get all vertices between the two endpoints along x order
        if(Less_indices_xy_2(poly[i-1],poly[i])){
          start=std::upper_bound(p_sort_by_x.begin(), p_sort_by_x.end(), poly[i-1]);
          end=std::lower_bound(p_sort_by_x.begin(), p_sort_by_x.end(), poly[i]);
          for(auto it=start; it!=end; ++it){
            updated_polyline.push_back(*it);
          }
        } else {
          start=std::upper_bound(p_sort_by_x.begin(), p_sort_by_x.end(), poly[i]);
          end=std::lower_bound(p_sort_by_x.begin(), p_sort_by_x.end(), poly[i-1]);
          for(auto it=end; it!=start;){
            --it;
            updated_polyline.push_back(*it);
          }
        }
      }
      updated_polyline.push_back(poly[i]);
    }
    std::swap(poly, updated_polyline);
  }
}

template <class SegmentRange, class Traits, class PointsRange , class PolylinesRange>
void vertical_slab_snap_rounding_2_impl(const SegmentRange &segments, PointsRange &pts, PolylinesRange &polylines, const Traits &traits){
  using Point_2 = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;

  auto to_exact = traits.converter_to_exact_object();
  auto round    = traits.construct_rounded_point_2_object();

  std::vector< Segment_2 > convert_input;
  for(const auto& seg: segments)
    if(seg.source()!=seg.target())
      convert_input.push_back(to_exact(seg));
#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Solved intersections" << std::endl;
  std::cout << "do intersect? " << do_curves_intersect(convert_input.begin(), convert_input.end()) << std::endl;
#endif
  compute_intersection_polylines(convert_input.begin(), convert_input.end(), pts, polylines);

#ifdef CGAL_DOUBLE_2D_SNAP_FULL_VERBOSE
  std::cout << "Input points" << std::endl;
  std::size_t i=0;
  for(const auto &p: pts)
    std::cout << i++ << ": " << p << std::endl;
  std::cout << std::endl;
  std::cout << "Input polylines" << std::endl;
  i=0;
  for(const auto &pl: polylines){
    std::cout << i++ << ":";
    for(std::size_t pi: pl)
      std::cout << " " << pi;
    std::cout << std::endl;
  }
  std::cout << std::endl;
#endif

  snap_rounding_scan(pts, polylines, traits);
#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Round" << std::endl;
#endif
  for(Point_2 &p: pts)
    p=round(p);
#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Merging points that are collapsed together" << std::endl;
#endif
  merge_duplicate_points_in_polylines(pts, polylines, traits);
#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  // The algorithm prevents a vertex that goes through a segment but a vertex may lie on a horizontal/vertical segment after rounding
  std::cout << "Subdivide vertical segments with vertices on them" << std::endl;
#endif
  snap_post_process(pts, polylines, traits);

#ifdef CGAL_DOUBLE_2D_SNAP_FULL_VERBOSE
  {
    std::cout << "Output points" << std::endl;
    std::size_t i=0;
    for(const auto &p: pts)
      std::cout << i++ << ": " << p << std::endl;
    std::cout << std::endl;
    std::cout << "Output polylines" << std::endl;
    i=0;
    for(const auto &pl: polylines){
      std::cout << i++ << ":";
      for(std::size_t pi: pl)
        std::cout << " " << pi;
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
#endif
}

template <class PolygonRange, class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
void vertical_slab_snap_rounding_2_polygon(const PolygonRange  &polygons,
                                           OutputContainer &out,
                                           const NamedParameters &np = parameters::default_values())
{
  using Polygon_2 = typename std::iterator_traits<typename PolygonRange::iterator>::value_type;
  using InputKernel = typename Kernel_traits<typename Polygon_2::Point_2>::Kernel;
  using DefaultTraits = Double_grid_snap_rounding_traits_2<InputKernel>;
  using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                              NamedParameters,
                                                              DefaultTraits>::type;
  using Point_2 = typename Traits::Point_2;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const Traits &traits = choose_parameter(get_parameter(np, internal_np::geom_traits), DefaultTraits());

  auto from_exact= traits.converter_from_exact_object();

#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Change format to range of points and indices" << std::endl;
#endif
  std::vector< typename InputKernel::Segment_2 > input_segments;
  std::vector< Point_2 > pts;
  std::vector< std::vector< std::size_t> > polylines;

  // Store the indices of segment of a new polygon, segments between [ polygon_index[i] and polygon_index[i+1] [ belong to polygon i
  std::vector< std::size_t > polygon_indices;

  polygon_indices.reserve(std::distance(polygons.begin(), polygons.end()));
  for(const Polygon_2 &P: polygons){
    polygon_indices.push_back(input_segments.size());
    for(std::size_t i=0; i<P.size()-1; ++i)
      input_segments.emplace_back(P[i], P[i+1]);
    input_segments.emplace_back(P[P.size()-1], P[0]);
  }
  polygon_indices.push_back(input_segments.size());

  // Main algorithm
  internal::vertical_slab_snap_rounding_2_impl(input_segments, pts, polylines, traits);

#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Build output" << std::endl;
#endif

  // Reassemble the polygons
  for(std::size_t polygon_idx = 0; polygon_idx != polygon_indices.size()-1; ++polygon_idx){
    Polygon_2 P;
    std::size_t idx_start = polygon_indices[polygon_idx];
    std::size_t idx_end = polygon_indices[polygon_idx+1];
    std::size_t last_insert;
    for(std::size_t pl_idx = idx_start; pl_idx != idx_end; ++pl_idx){
      auto &pl = polylines[pl_idx];
      // The first element is not add, it is identical to the last

      bool go_forward;
      if(pl_idx == idx_start)
        go_forward = (pl.back() == polylines[pl_idx+1].front()) || (pl.back() == polylines[pl_idx+1].back());
      else
        go_forward = (last_insert == pl.front());

      // Add the element in forward direction
      if(go_forward){
        for(std::size_t i = 1; i != pl.size(); ++i)
          P.push_back(from_exact(pts[pl[i]]));
        last_insert = pl.back();

      // Add the element in backward direction
      } else {
        for(std::size_t i = pl.size()-1; i != 0; --i)
          P.push_back(from_exact(pts[pl[i-1]]));
        last_insert = pl.front();
      }
    }
    out.push_back(P);
  }
}

} // end of namespace internal

#if DOXYGEN_RUNNING
/**
* \ingroup Snap_rounding_vertical_slab_grp
*
* \brief subdivides and rounds a range of segments so that they are pairwise disjoint in their interiors.
*
* By default, each polyline of the output corresponds to an input segment. Consequently, duplicate segments may appear in the output, for instance when multiple input segments collapse.
* When the parameter `output_unique_segments` is set to `true`, the polylines are decomposed into individual segments (represented as polylines with two points), and duplicates are removed.
*
* @tparam SegmentRange model of a `ConstRange` whose iterator is model of `ForwardIterator` and whose value_type is `geom_traits::Segment_2`, where the type of geom_traits is detailed by `np::geom_traits`.
* @tparam OutputContainer model of the concept `BackInsertionSequence` whose value type is itself a model of the concepts `DefaultConstructible` and `BackInsertionSequence` whose value type is `geom_traits::Point_2`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param segments the input segment range
* \param out the output container
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `VerticalSlabSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `CGAL::Double_grid_snap_rounding_traits_2<Kernel>` where Kernel is deduced from the segment type, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
*   \cgalParamNBegin{output_unique_segments}
*     \cgalParamDescription{If set to true, the output polylines are unique pairs of distinct points represented a segment. As a result, the total number of output polylines may differ from the number of input segments.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{true}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template <class SegmentRange , class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
void vertical_slab_snap_rounding_2(const SegmentRange &segments,
                                   OutputContainer    &out,
                                   const NamedParameters &np = parameters::default_values());

/**
* \ingroup Snap_rounding_vertical_slab_grp
*
* \brief subdivides and rounds a range of polygons so that their boundary segments are pairwise disjoint in their interiors.
*
* If the input polygons are disjoint, the output polygons remain non-overlapping, although they may share vertices or edges.
* Each output polygon is free of self-intersections but may present pinched sections.
*
* @tparam PolygonRange model of a `ConstRange` whose iterator is model of `ForwardIterator` and whose value_type is model of `CGAL::Polygon_2`.
* @tparam OutputContainer model of the concept `BackInsertionSequence` whose value type is model of `CGAL::Polygon_2`.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param polygons the range of input polygons
* \param out the output container
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `VerticalSlabSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `CGAL::Double_grid_snap_rounding_traits_2<Kernel>` where Kernel is deduced from the point type of the polygons, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
* @warning a convex input polygon might no longer be convex after rounding.
*/
template <class PolygonRange, class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
void vertical_slab_snap_rounding_2(const PolygonRange  &polygons,
                                   OutputContainer     &out,
                                   const NamedParameters &np = parameters::default_values());

#else

template <class InputRange , class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
void vertical_slab_snap_rounding_2(const InputRange &inputs,
                                   OutputContainer  &out,
                                   const NamedParameters &np = parameters::default_values())
{
  using Input = std::remove_cv_t<typename std::iterator_traits<typename InputRange::iterator>::value_type>;

  if constexpr(internal::is_instance_of_Polygon_2< Input >){
    internal::vertical_slab_snap_rounding_2_polygon(inputs, out, np);
    return;
  } else {
    using Polyline = std::remove_cv_t<typename OutputContainer::value_type>;

    using InputKernel = typename Kernel_traits<std::remove_cv_t<typename std::iterator_traits<typename InputRange::iterator>::value_type>>::Kernel;
    using DefaultTraits = Double_grid_snap_rounding_traits_2<InputKernel>;
    using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                                NamedParameters,
                                                                DefaultTraits>::type;

    using Point_2 = typename Traits::Point_2;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    Traits traits = choose_parameter(get_parameter(np, internal_np::geom_traits), DefaultTraits());
    bool unique_segments = choose_parameter(get_parameter(np, internal_np::output_unique_segments), false);

    // auto to_exact=   traits.converter_to_exact_object();
    auto from_exact= traits.converter_from_exact_object();

    // Main algorithm
    std::vector<Point_2> pts;
    std::vector< std::vector< std::size_t> > polylines;
    internal::vertical_slab_snap_rounding_2_impl(inputs, pts, polylines, traits);

  #ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
    std::cout << "Build output" << std::endl;
  #endif

    if (unique_segments){
      // Output Segments while removing duplicate ones
      auto segment_2 = traits.construct_segment_2_object();
      std::set< std::pair<std::size_t,std::size_t> > set_out_segs;
      for(auto &poly: polylines)
        for(std::size_t i=1; i<poly.size(); ++i)
          set_out_segs.emplace((std::min)(poly[i-1],poly[i]),(std::max)(poly[i-1],poly[i]));
      for(auto &pair: set_out_segs){
        Polyline pl;
        pl.push_back(from_exact(pts[pair.first]));
        pl.push_back(from_exact(pts[pair.second]));
        out.push_back(std::move(pl));
        CGAL_assertion(pts[pair.first]!=pts[pair.second]);
      }
    } else {
      // Output polylines
      for(auto &poly: polylines){
        Polyline new_line;
        for(std::size_t pi: poly)
          new_line.push_back(from_exact(pts[pi]));
        out.push_back(std::move(new_line));
      }
    }
 }
}
#endif

} //namespace CGAL

#endif
