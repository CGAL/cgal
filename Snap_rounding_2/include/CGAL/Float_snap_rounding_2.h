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
// author(s)     : Léo Valque

#ifndef CGAL_FLOAT_SNAP_ROUNDING_2_H
#define CGAL_FLOAT_SNAP_ROUNDING_2_H

#include <CGAL/license/Snap_rounding_2.h>

#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
#include <iostream>
#endif

#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Float_snap_rounding_traits_2.h>
#include <CGAL/intersection_2.h>

#include <CGAL/internal/Wrap_float_snap_rounding_traits_2.h>

#include <boost/property_map/function_property_map.hpp>

#include <set>
#include <vector>
#include <type_traits>


namespace CGAL {

  namespace internal{

template <typename GeometryTraits_2, typename Points_, typename Polylines_, typename Allocator_ = CGAL_ALLOCATOR(int)>
class Snap_rounding_visitor :
    public Ss2::Default_visitor<Snap_rounding_visitor<GeometryTraits_2, Points_, Polylines_, Allocator_>,
                                GeometryTraits_2, Allocator_> {
public:
  using Geometry_traits_2 = GeometryTraits_2;
  using Points = Points_;
  using Polylines = Polylines_;
  using Allocator = Allocator_;

  using P_idx = std::size_t;
  using L_idx = std::size_t;
  using Polyline = typename Polylines::value_type;

private:
  using Gt2 = Geometry_traits_2;
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
  Gt2 traits;
  Points& points;
  const Polylines& input_polylines;
  Polylines& output_polylines;
  std::vector< double >& round_bounds;

  // If pi is closed enough to sc, subdivide sc in the output and return true
  bool is_pi_closed_to_sc_and_subdivide(P_idx pi, Subcurve* sc){
    auto point_at_x = traits.construct_point_at_x_on_segment_2_object();
    auto csq_dist_2 = traits.compare_squared_distance_2_object();
    auto round_bound = traits.compute_squared_round_bound_2_object();

    L_idx li = sc->last_curve().first.first;
    std::size_t i = sc->last_curve().first.second;

    P_idx src = input_polylines[li][i];
    P_idx trg = input_polylines[li][i+1];
    // Skip if this polyline already ends with the event point
    if(!output_polylines[li].empty() && output_polylines[li].back() == pi)
      return true;

    const auto &p = points[pi];
    const auto &pl = input_polylines[li];
    Base_segment_2 seg = Base_segment_2(points[pl.front()], points[pl.back()]); //get the input segment

    // (A+B)^2 <= 4*max(A^2,B^2)
    double bound = (std::max)({ round_bounds[pi],
                                round_bounds[src],
                                round_bounds[trg] });
    bound *= 4;

    if(possibly(csq_dist_2(p, seg, bound) != LARGER)){
      // if constexpr (std::is_same_v<CGAL::Exact_predicates_exact_constructions_kernel, typename G2t::Exact_type>){
      if(true){
        internal::Evaluate<FT> evaluate;
        // We refine the pts to reduce the rounding shift and check again
        evaluate(points[pi]);
        evaluate(points[output_polylines[li].back()]);
        // The two following call of exact act on seg variables since they appear in its DAG
        evaluate(points[src]);
        evaluate(points[trg]);
        // Update the bounds
        round_bounds[pi]  = round_bound(points[pi]);
        round_bounds[output_polylines[li].back()] = round_bound(points[output_polylines[li].back()]);
        round_bounds[src] = round_bound(points[src]);
        round_bounds[trg] = round_bound(points[trg]);
        bound = (std::max)({ round_bounds[pi],
                              round_bounds[output_polylines[li].back()],
                              round_bounds[src],
                              round_bounds[trg] });
        bound *= 4;

        // Check if the point and the segment are still too closed for a safe rounding
        if(csq_dist_2(p, seg, bound) == LARGER)
          return false;
      }

#ifdef CGAL_DOUBLE_2D_SNAP_FULL_VERBOSE
      std::cout << "Create point " << new_pi << " on " << li << " due to proximity with " << pi << "_____________________________" << std::endl;
      std::cout << new_pi <<": " << pts[new_pi] << std::endl;
      std::cout << li << ":";
      for(std::size_t i: pl)
        std::cout << " " << i;
      std::cout << std::endl;
#endif

      // Check if segment was not already subdivided by another point
      if(points[pl.back()].x() == points[pi].x())
        return false;

      // Create a point on seg at the same x coordinate than p
      points.push_back(point_at_x(seg, points[pi].x()));
      round_bounds.push_back(round_bound(points.back()));
      P_idx new_pi=points.size()-1;
      // We insert it on the output
      output_polylines[li].push_back(new_pi);
      pi = new_pi;
      return true;
    }
    return false;
  }

public:
  Snap_rounding_visitor(Gt2 &traits_, Points& points_, Polylines& in, Polylines& out, std::vector< double >& bounds) :
    traits(traits_),
    points(points_),
    input_polylines(in),
    output_polylines(out),
    round_bounds(bounds)
  {}

  /*!
   */
  template <typename CurveIterator>
  void sweep(CurveIterator begin, CurveIterator end) {
    output_polylines.resize(input_polylines.size());
    Surface_sweep_2* sl = this->surface_sweep();
    sl->sweep(begin, end);
  }

  /*!
   */
  bool after_handle_event(Event* event, Status_line_iterator iter, bool /* flag */) {
    if (! event->is_closed()) return true;
    if (! event->has_left_curves() && ! event->has_right_curves()) return true;

    auto pi = event->point().idx;
    Surface_sweep_2* sl = this->surface_sweep();

    // Insert the point in the output
    for (auto it = event->left_curves_begin(); it != event->left_curves_end(); ++it) {
      Subcurve *sc = *it;
      std::vector<Subcurve*> leaves;
      sc->all_leaves(std::back_inserter(leaves));
      for(auto l: leaves)
        output_polylines[l->last_curve().first.first].push_back(pi);
    }

    for (auto it = event->right_curves_begin(); it != event->right_curves_end(); ++it) {
      Subcurve *sc = *it;
      std::vector<Subcurve*> leaves;
      sc->all_leaves(std::back_inserter(leaves));
      for(auto l: leaves)
        if(output_polylines[l->last_curve().first.first].empty() || output_polylines[l->last_curve().first.first].back() != pi) // Maybe already by left curves
          output_polylines[l->last_curve().first.first].push_back(pi);
    }

    // Look above segments and creates a point if there too close for a safe rounding
    auto pi_above = pi;
    auto above = iter;
    while(above != sl->status_line_end() && is_pi_closed_to_sc_and_subdivide(pi_above, *above))
      ++above;

    // same with below segments
    auto below = iter;
    if(below != sl->status_line_begin()){
      --below;
      while(is_pi_closed_to_sc_and_subdivide(pi, *below)){
        if(below == sl->status_line_begin())
          break;
        --below;
      }
    }

    return true;
  }
};

/*
Scan the vertices from left to right while maintening the status line order of the segments.
Subdivide the segments if there are too close to a vertex
*/
template <class Concurrency_tag=Sequential_tag, class Traits, class PointsRange , class PolylinesRange>
void snap_rounding_scan(PointsRange &pts, PolylinesRange &polylines, const Traits &traits){

#ifdef CGAL_DOUBLE_2D_SNAP_FULL_VERBOSE
  std::cout << "Input points" << std::endl;
  std::size_t i=0;
  for(const Point_2 &p: pts)
    std::cout << i++ << ": " << p << std::endl;
  std::cout << std::endl;
  std::cout << "Input polylines" << std::endl;
  i=0;
  for(const Polyline &pl: polylines){
    std::cout << i++ << ":";
    for(std::size_t pi: pl)
      std::cout << " " << pi;
    std::cout << std::endl;
  }
  std::cout << std::endl;
#endif

// Wrap the traits with something supporting the data structure pts, polylines
  struct Key{
    Key(){}
    Key(std::size_t f, std::size_t s, std::size_t i):first(f), second(s), idx(i){};
    bool operator==(Key b) const{
      return idx==b.idx;
    }
    std::size_t first;
    std::size_t second;
    std::size_t idx;
  };
  auto get = [&](const Key &idx){ return pts[idx.idx]; };
  auto wrap_traits = make_wrap_float_snap_rounding_traits_2(traits, boost::make_function_property_map<Key>(get),[](){});

  auto round_bound = traits.compute_squared_round_bound_2_object();

  using Wrap_traits = decltype(wrap_traits);
  using X_monotone_curve_2 = typename Wrap_traits::X_monotone_curve_2;
  using Visitor = Snap_rounding_visitor<Wrap_traits, PointsRange, PolylinesRange>;
  using Surface_sweep = Ss2::No_intersection_surface_sweep_2<Visitor>;

  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename PolylinesRange::iterator>::value_type>;

  // auto less_xy_2 = traits.less_xy_2_object();
  auto less_y_2 = traits.less_y_2_object();

  std::vector< double > round_bound_pts;
  round_bound_pts.reserve(pts.size());
  for(std::size_t i=0; i!=pts.size(); ++i)
    round_bound_pts.push_back(round_bound(pts[i]));

  // Create the curves
  std::vector<X_monotone_curve_2> curves;
  for(std::size_t j=0; j!=polylines.size(); ++j)
    for(std::size_t i=0; i!=polylines[j].size()-1; ++i)
      curves.emplace_back(Key(j, i, polylines[j][i]), Key(j, i+1, polylines[j][i+1]));

  std::vector< Polyline > out;
  Visitor visitor(wrap_traits, pts, polylines, out, round_bound_pts);
  Surface_sweep surface_sweep(&wrap_traits, &visitor);
  visitor.sweep(curves.begin(), curves.end());

  // We sort the point by y to ensure the order of the point is preserved
  // Round value of y coordinates are indirectly modified in the case of a filter failure
  auto sort_pi=[&](std::size_t i, std::size_t j){
    return less_y_2(pts[i],pts[j]);
  };
  std::vector< std::size_t > indices(pts.size(),0);
  std::iota(indices.begin(),indices.end(),0);
  std::sort(indices.begin(),indices.end(),sort_pi);

  std::swap(polylines, out);
#ifdef CGAL_DOUBLE_2D_SNAP_FULL_VERBOSE
  std::cout << "Output points" << std::endl;
  i=0;
  for(const Point_2 &p: pts)
    std::cout << i++ << ": " << p << std::endl;
  std::cout << std::endl;
  std::cout << "Output polylines" << std::endl;
  i=0;
  for(const Polyline &pl: polylines){
    std::cout << i++ << ":";
    for(std::size_t pi: pl)
      std::cout << " " << pi;
    std::cout << std::endl;
  }
  std::cout << std::endl;
#endif
  return;
}

template <class Concurrency_tag=Sequential_tag, class Traits, class PointsRange , class PolylinesRange>
void merge_duplicate_points_in_polylines(PointsRange &pts, PolylinesRange &polylines, const Traits &traits)
{
  using Point_2 = typename Traits::Point_2;
  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename PolylinesRange::iterator>::value_type>;

  using Less_xy_2 = typename Traits::Less_xy_2;
  using Equal_2 = typename Traits::Equal_2;
  Equal_2 equal = traits.equal_2_object();

  auto Less_indexes_xy_2=[&](std::size_t i, std::size_t j){
    return Less_xy_2()(pts[i], pts[j]);
  };

  std::vector< std::size_t > unique_points(pts.size());
  std::iota(unique_points.begin(), unique_points.end(), 0);
  std::sort(unique_points.begin(), unique_points.end(), Less_indexes_xy_2);
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
    for (std::size_t i=0; i<polyline.size(); ++i) {
        std::size_t new_pi=old_to_new_index[polyline[i]];
        if(i==0 || (new_pi!=updated_polyline[updated_polyline.size()-1]))
          updated_polyline.push_back(new_pi);
        assert(new_pi<pts.size());
        assert(pts[new_pi]==new_pts[polyline[i]]);
    }
    std::swap(polyline, updated_polyline);
  }
}

// Some points may have collapsed on a vertical segment, we subdivide these vertical segments accordingly
template <class Concurrency_tag=Sequential_tag, class Traits, class PointsRange , class PolylinesRange>
void snap_post_process(PointsRange &pts, PolylinesRange &polylines, const Traits &traits)
{
  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename PolylinesRange::iterator>::value_type>;
  using Less_xy_2 = typename Traits::Less_xy_2;

  Less_xy_2 less_xy_2 = traits.less_xy_2_object();

  auto Less_indexes_xy_2=[&](std::size_t i, std::size_t j){
    return less_xy_2(pts[i], pts[j]);
  };

  std::vector< std::size_t > p_sort_by_x(pts.size());
  std::iota(p_sort_by_x.begin(), p_sort_by_x.end(), 0);
  std::sort(p_sort_by_x.begin(), p_sort_by_x.end(), Less_indexes_xy_2);

  for(Polyline &poly: polylines){
    std::vector<std::size_t> updated_polyline;
    updated_polyline.push_back(poly.front());
    for(std::size_t i=1; i!=poly.size(); ++i){
      if(pts[poly[i-1]].x()==pts[poly[i]].x()){
        std::vector< std::size_t >::iterator start, end;
        // Get all vertices between the two endpoints along x order
        if(Less_indexes_xy_2(poly[i-1],poly[i])){
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

template <class Concurrency_tag=Sequential_tag, class InputIterator, class Traits, class PointsRange , class PolylinesRange>
void double_snap_rounding_2_impl(InputIterator begin, InputIterator end, PointsRange &pts, PolylinesRange &polylines, const Traits &traits){
  using Point_2 = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;

  auto to_exact   = traits.converter_to_exact_object();
  // auto from_exact = traits.converter_from_exact_object();
  auto round      = traits.construct_rounded_point_2_object();

  std::vector< Segment_2 > convert_input;
  for(InputIterator it=begin; it!=end; ++it)
    if(it->source()!=it->target())
      convert_input.push_back(to_exact(*it));
  std::vector<Segment_2> segs;
#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Solved intersections" << std::endl;
  std::cout << "do intersect? " << do_curves_intersect(convert_input.begin(), convert_input.end()) << std::endl;
#endif
  compute_intersection_polylines(convert_input.begin(), convert_input.end(), pts, polylines);

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
  // The algorithm prevents the a vertex that goes through a segment but a vertex may lie on an horizontal/vertical segments after rounding
  std::cout << "Subdivide vertical segments with vertices on them" << std::endl;
#endif
  snap_post_process(pts, polylines, traits);
}

} // end of namespace internal

/**
* \ingroup PkgSnapRounding2Ref
*
* Subdivides and rounded a set of segments so that they are pairwise disjoint in their interiors.
* The output is a range of polyline with each polyline corresponding to an input segment.
*
* TODO Currently compute_subcurves have no visitor to obtain a polyline per input segment, thus
* currently the output contain one polyline per ubsegment of the arrangement
*
* @tparam InputIterator iterator over a range of `Segment_2`
* @tparam OutputContainer inserter over a range of `Polyline`. `Polyline` must be a type that provides a `push_back(Point_2)` function.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param begin,end the input segment range
* \param out the output inserter
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{concurrency_tag}
*     \cgalParamDescription{That template parameter enables to choose whether the algorithm is to be run in parallel, if CGAL::Parallel_tag
*                           is specified and CGAL has been linked with the Intel TBB library, or sequentially, otherwise.}
*     \cgalParamType{CGAL::Concurrency_tag}
*     \cgalParamDefault{CGAL::Sequential_tag}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `FloatSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `Float_snap_rounding_traits_2`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template <class InputIterator , class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
OutputContainer double_snap_rounding_2(InputIterator    begin,
                                       InputIterator    end,
                                       OutputContainer  out,
                                       const NamedParameters &np = parameters::default_values())
{
  using Concurrency_tag = typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t,
                                                              NamedParameters,
                                                              Sequential_tag>::type;

  using Polyline = std::remove_cv_t<typename OutputContainer::container_type::value_type>;

  using InputKernel = typename Kernel_traits<std::remove_cv_t<typename std::iterator_traits<InputIterator>::value_type>>::Kernel;
  using DefaultTraits = Float_snap_rounding_traits_2<InputKernel>;
  using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                              NamedParameters,
                                                              DefaultTraits>::type;

  using Point_2 = typename Traits::Point_2;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const Traits &traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Traits());

  // auto to_exact=   traits.converter_to_exact_object();
  auto from_exact= traits.converter_from_exact_object();

  // Main algorithm
  std::vector<Point_2> pts;
  std::vector< std::vector< std::size_t> > polylines;
  internal::double_snap_rounding_2_impl<Concurrency_tag>(begin, end, pts, polylines, traits);

#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Build output" << std::endl;
#endif

  // Output polylines
  for(auto &poly: polylines){
    Polyline new_line;
    for(std::size_t pi: poly)
      new_line.push_back(from_exact(pts[pi]));
    *out++ = new_line;
  }

  return out;
}

/**
* \ingroup PkgSnapRounding2Ref
*
* Given a range of segments, compute rounded subsegments that are pairwise disjoint in their interior, as induced by the input curves.
*
* @tparam Concurrency_tag That template parameter enables to choose whether the algorithm is to be run in
* parallel, if CGAL::Parallel_tag is specified and CGAL has been linked with the Intel TBB library, or sequentially, if CGAL::Sequential_tag - the default value - is specified.
* @tparam InputIterator iterator of a segment range
* @tparam OutputContainer inserter of a segment range
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param begin,end the input segment range
* \param out the output inserter
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{concurrency_tag}
*     \cgalParamDescription{That template parameter enables to choose whether the algorithm is to be run in parallel, if CGAL::Parallel_tag
*                           is specified and CGAL has been linked with the Intel TBB library, or sequentially, otherwise.}
*     \cgalParamType{CGAL::Concurrency_tag}
*     \cgalParamDefault{CGAL::Sequential_tag}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `FloatSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `Float_snap_rounding_traits_2`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template <class InputIterator , class OutputIterator, class NamedParameters = parameters::Default_named_parameters>
OutputIterator compute_snapped_subcurves_2(InputIterator     begin,
                                           InputIterator     end,
                                           OutputIterator    out,
                                           const NamedParameters &np = parameters::default_values())
{
  using Concurrency_tag = typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t,
                                                              NamedParameters,
                                                              Sequential_tag>::type;

  using InputKernel = typename Kernel_traits<std::remove_cv_t<typename std::iterator_traits<InputIterator>::value_type>>::Kernel;
  using DefaultTraits = Float_snap_rounding_traits_2<InputKernel>;
  using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                              NamedParameters,
                                                              DefaultTraits>::type;

  using Point_2 = typename Traits::Point_2;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const Traits &traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Traits());

  // auto to_exact=   traits.converter_to_exact_object();
  auto from_exact= traits.converter_from_exact_object();
  auto segment_2 = traits.construct_segment_2_object();

  // Main algorithm
  std::vector<Point_2> pts;
  std::vector< std::vector< std::size_t> > polylines;
  internal::double_snap_rounding_2_impl<Concurrency_tag>(begin, end, pts, polylines, traits);

#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Build output" << std::endl;
#endif

  // Output a range of segments while removing duplicate ones
  std::set< std::pair<std::size_t,std::size_t> > set_out_segs;
  for(auto &poly: polylines)
    for(std::size_t i=1; i<poly.size(); ++i)
      set_out_segs.emplace((std::min)(poly[i-1],poly[i]),(std::max)(poly[i-1],poly[i]));
  for(auto &pair: set_out_segs){
    *out++=from_exact(segment_2(pts[pair.first], pts[pair.second]));
    assert(pts[pair.first]!=pts[pair.second]);
  }

  return out;
}

/**
* \ingroup PkgSnapRounding2Ref
*
* Given a range of `Polygon_2`, compute rounded polygons such that their segments are either equal either disjoint in their interior, as induced by the input polygons.
* The polygons are intended to be non-intersecting, unless the named parameter `compute_intersections` is set to `true`.
* Any polygon is guaranteed to remain a Polygon in the output but may present pinched section or/and common vertices or segments with
* other polygons.
*
* @tparam InputIterator iterator of a CGAL::Polygon_2 range
* @tparam OutputContainer inserter of a CGAL::Polygon_2 range
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param begin,end the input polygon range
* \param out the output inserter
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{concurrency_tag}
*     \cgalParamDescription{That template parameter enables to choose whether the algorithm is to be run in parallel, if CGAL::Parallel_tag
*                           is specified and CGAL has been linked with the Intel TBB library, or sequentially, otherwise.}
*     \cgalParamType{CGAL::Concurrency_tag}
*     \cgalParamDefault{CGAL::Sequential_tag}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `FloatSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `Float_snap_rounding_traits_2`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
* @warning If an input polygon is convex, it might no longer be convex in the output of this function
*/
template <class InputIterator, class OutputIterator, class NamedParameters = parameters::Default_named_parameters>
void compute_snapped_polygons_2(InputIterator  begin,
                                InputIterator  end,
                                OutputIterator out,
                                const NamedParameters &np = parameters::default_values())
{
  using Concurrency_tag = typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t,
                                                              NamedParameters,
                                                              Sequential_tag>::type;

  using Polygon_2 = typename std::iterator_traits<InputIterator>::value_type;
  using InputKernel = typename Kernel_traits<typename Polygon_2::Point_2>::Kernel;
  using DefaultTraits = Float_snap_rounding_traits_2<InputKernel>;
  using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                              NamedParameters,
                                                              DefaultTraits>::type;
  using Point_2 = typename Traits::Point_2;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const Traits &traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Traits());

  auto from_exact= traits.converter_from_exact_object();

#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Change format to range of points and indexes" << std::endl;
#endif
  std::vector< typename InputKernel::Segment_2 > input_segments;
  std::vector< Point_2 > pts;
  std::vector< std::vector< std::size_t> > polylines;

  // Store the indexes of segment of a new polygon, segments between [ polygon_index[i] and polygon_index[i+1] [ belong to polygin i
  std::vector< std::size_t > polygon_indexes;

  polygon_indexes.reserve(std::distance(begin, end));
  for(InputIterator it=begin; it!=end; ++it){
    polygon_indexes.push_back(input_segments.size());
    const Polygon_2 &P = *it;
    for(std::size_t i=0; i<P.size()-1; ++i)
      input_segments.emplace_back(P[i], P[i+1]);
  }
  polygon_indexes.push_back(input_segments.size());

  // Main algorithm
  internal::double_snap_rounding_2_impl<Concurrency_tag>(input_segments.begin(), input_segments.end(), pts, polylines, traits);

#ifdef CGAL_DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Build output" << std::endl;
#endif

  // Reassemble the polygons
  for(std::size_t polygon_idx = 0; polygon_idx != polygon_indexes.size()-1; ++polygon_idx){
    Polygon_2 P;
    std::size_t idx_start = polygon_indexes[polygon_idx];
    std::size_t idx_end = polygon_indexes[polygon_idx+1];
    std::size_t last_insert;
    for(std::size_t pl_idx = idx_start; pl_idx != idx_end; ++pl_idx){
      auto &pl = polylines[pl_idx];
      // Add the first element
      if(pl_idx == idx_start)
        P.push_back(from_exact(pts[pl.front()]));

      // Add the element in forward direction
      if(pl_idx == idx_start || last_insert == pl.front())
        for(std::size_t i = 1; i != pl.size(); ++i)
          P.push_back(from_exact(pts[pl[i]]));

      // Add the element in backward direction
      else
        for(std::size_t i = pl.size(); i != 0; --i)
          P.push_back(from_exact(pts[pl[i-1]]));
    }
    *out++=P;
  }
}

/**
* \ingroup PkgSnapRounding2Ref
*
* Given a Polygon_2, compute rounded segments that are pairwise disjoint in their interior, as induced by the input polygon.
* The output is guarantee to be a Polygon but may present pinched section.
*
* @tparam Polygon_2 model of `CGAL::Polygon_2`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param P the input polygon
* \param out the output polygon
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{concurrency_tag}
*     \cgalParamDescription{That template parameter enables to choose whether the algorithm is to be run in parallel, if CGAL::Parallel_tag
*                           is specified and CGAL has been linked with the Intel TBB library, or sequentially, otherwise.}
*     \cgalParamType{CGAL::Concurrency_tag}
*     \cgalParamDefault{CGAL::Sequential_tag}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{a multiplier of the error value for boundary edges to preserve the boundaries}
*     \cgalParamType{double}
*     \cgalParamDefault{100}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
* \warning The convex property is not necessarly preserved
*/
template <class Polygon_2, class NamedParameters = parameters::Default_named_parameters>
void compute_snapped_polygon_2(const Polygon_2 &P, Polygon_2 &out, const NamedParameters &np = parameters::default_values())
{
  std::array<Polygon_2, 1> vec({P});
  std::vector<Polygon_2> out_vec;
  compute_snapped_polygons_2(vec.begin(), vec.end(), std::back_inserter(out_vec), np);
  out = out_vec[0];
}

} //namespace CGAL

#endif
