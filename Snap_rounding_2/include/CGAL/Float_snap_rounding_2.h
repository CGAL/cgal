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

#ifdef DOUBLE_2D_SNAP_VERBOSE
#include <iostream>
#endif

#include <CGAL/Surface_sweep_2_algorithms.h>

#include <CGAL/Float_snap_rounding_traits_2.h>

#include <CGAL/intersection_2.h>
#include <set>
#include <vector>

#include <CGAL/Named_function_parameters.h>

#include  <CGAL/mutex.h>

namespace CGAL {

namespace internal{
  /*
  Scan the vertices from left to right while maintening the y order of the segments.
  Subdivide the segments if there are too close to a vertex
  */
template <class Concurrency_tag=Sequential_tag, class Traits, class PointsRange , class PolylineRange>
void snap_rounding_scan(PointsRange &pts, PolylineRange &polylines, const Traits &traits){
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;

  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename PolylineRange::iterator>::value_type>;

  using Less_xy_2 = typename Traits::Less_xy_2;
  using Less_y_2 = typename Traits::Less_y_2;
  using Construct_segment_2 = typename Traits::Construct_segment_2;
  using Construct_point_at_x_on_segment_2 = typename Traits::Construct_point_at_x_on_segment_2;

  using Compare_squared_distance_2 = typename Traits::Compare_squared_distance_2;
  using Squared_round_bound_2      = typename Traits::Squared_round_bound_2;

  Less_xy_2 less_xy_2 = traits.less_xy_2_object();
  Less_y_2 less_y_2 = traits.less_y_2_object();
  Construct_segment_2 segment_2 = traits.construct_segment_2_object();
  Construct_point_at_x_on_segment_2 point_at_x = traits.construct_point_at_x_on_segment_2_object();

  Compare_squared_distance_2 csq_dist_2 = traits.compare_squared_distance_2_object();
  Squared_round_bound_2 round_bound = traits.squared_round_bound_2_object();

  // Precompute round bounds for points
  std::vector<double> round_bound_pts;
  round_bound_pts.reserve(pts.size());
  for(std::size_t i=0; i<pts.size(); ++i)
    round_bound_pts.emplace_back(round_bound(pts[i]));

  enum class EVENT_TYPE{INSERT, MIDDLE, REMOVE};
  struct Event{
    Event(size_t pi_, size_t li_, EVENT_TYPE t_):pi(pi_),li(li_),type(t_){}
    size_t pi;
    size_t li;
    EVENT_TYPE type;
  };

  auto event_comparator=[&](Event a, Event b){
    const Point_2 &pa=pts[a.pi];
    const Point_2 &pb=pts[b.pi];

    // We sort event lexicographically
    if(a.pi!=b.pi)
      return less_xy_2(pa,pb);

    // If two events use the same point, we place remove event first for a smaller y_order range
    if(a.type != b.type)
      return a.type > b.type;
    // Only to have a complete order
    return a.li < b.li;
  };

  auto pi_below_li=[&](size_t li, std::pair<size_t, size_t> pair_pi_li){
    std::size_t ind_pi_support = pair_pi_li.second;
    std::size_t pi = pair_pi_li.first;
    if(ind_pi_support==li) return false;

    const Polyline &pl=polylines[li];

    Orientation ori=COLLINEAR;
    if(pl.front()!=pi && pl.back()!=pi)
      ori=orientation(pts[pl.front()], pts[pl.back()], pts[pi]);
    // When collinear, pick an other vertex of the line of pi to decide
    if(ori==COLLINEAR){
      const Polyline &pi_support = polylines[ind_pi_support];
      std::size_t other_point = (pi_support.front()!=pi) ? pi_support.front() : pi_support.back();
      ori=orientation(pts[pl.front()], pts[pl.back()], pts[other_point]);
      assert(ori!=COLLINEAR);
    }
    return ori==RIGHT_TURN;
  };

#ifdef DOUBLE_2D_SNAP_FULL_VERBOSE
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

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Create the Event Queue" << std::endl;
#endif

  // Input polylines are supposed going from left to right
  std::vector<Event> event_queue;
  event_queue.reserve(2*polylines.size());
  for(size_t li=0; li<polylines.size(); ++li){
    Polyline &pl=polylines[li];
    event_queue.emplace_back(pl.front(), li, EVENT_TYPE::INSERT);
    for(size_t pi=1; pi<pl.size()-1; ++pi){
      event_queue.emplace_back(pl[pi], li, EVENT_TYPE::MIDDLE);
      assert(less_xy_2(pts[pl[pi-1]], pts[pl[pi]]));
    }
    event_queue.emplace_back(pl.back(), li, EVENT_TYPE::REMOVE);
    assert(less_xy_2(pts[pl.front()], pts[pl.back()]));
  }
#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Sort the event queue along x direction" << std::endl;
#endif
  std::sort(event_queue.begin(), event_queue.end(), event_comparator);
#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Goes through Event" << std::endl;
#endif
  // TODO Vector is suboptimal for arbitrary insertion (Actually negligeable in the running time)
  // Track order of the lines along y along the events
  std::vector< size_t > y_order;
  y_order.push_back(event_queue[0].li);
  for(size_t i=1; i<event_queue.size();++i)
  {
    Event event = event_queue[i];
#ifdef DOUBLE_2D_SNAP_FULL_VERBOSE
    std::cout << ((event.type==EVENT_TYPE::INSERT)?"Insert": "Remove") << " Event at point " << event.pi << " on line "<< event.li << std::endl;
#endif
    // Find the position of the event along the y direction
    std::vector<size_t>::iterator pos_it = y_order.begin();
    if(!y_order.empty())
      pos_it = std::lower_bound(y_order.begin(), y_order.end(), std::pair<size_t, size_t>(event.pi, event.li), pi_below_li);
    if(event.type==EVENT_TYPE::REMOVE){
      assert(*pos_it==event.li);
      y_order.erase(pos_it);
    }

    if(!y_order.empty()){
      auto close_event=[&](std::size_t pi, std::size_t li){
        if((pi==polylines[li].front()) || (pi==polylines[li].back()))
          return true;

        const Point_2 &p = pts[pi];
        Polyline &pl=polylines[li];
        Segment_2 seg=segment_2(pts[pl.front()], pts[pl.back()]);

        // (A+B)^2 <= 4*max(A^2,B^2)
        double bound=round_bound_pts[pi];
        bound = (std::max)(bound, round_bound_pts[pl.front()]);
        bound = (std::max)(bound, round_bound_pts[pl.back()]);
        bound*=4;

        if(possibly(csq_dist_2(p, seg, bound)!=CGAL::LARGER))
        {
          if constexpr(std::is_same_v<Exact_predicates_exact_constructions_kernel, typename Traits::Exact_type>){
            // We refine the pts to reduce the rounding shift and check again
            pts[pi].exact();
            pts[pl[pl.size()-2]].exact();
            // The two following call of exact act on seg variables since they appear in its DAG
            pts[pl.back()].exact();
            pts[pl.front()].exact();
            // Update the bounds
            round_bound_pts[pi]=round_bound(pts[pi]);
            round_bound_pts[pl[pl.size()-2]]=round_bound(pts[pl[pl.size()-2]]);
            round_bound_pts[pl.back()]=round_bound(pts[pl.back()]);
            round_bound_pts[pl.front()]=round_bound(pts[pl.front()]);
            bound=round_bound_pts[pi];
            bound = (std::max)(bound, round_bound_pts[pl[pl.size()-2]]);
            bound = (std::max)(bound, round_bound_pts[pl.back()]);
            bound*=4;

            // Check if the point and the segment are still too closed for a safe rounding
            if(csq_dist_2(p, seg, bound)==CGAL::LARGER)
              return false;
          }

          // Check if segment was not already subdivided by another point
          for(std::size_t i: pl)
            if(pts[i].x()==pts[pi].x())
              return false;

          // Create a point on seg at the same x coordinate than p
          pts.push_back(point_at_x(seg, pts[event.pi].x()));
          std::size_t new_pi=pts.size()-1;
          round_bound_pts.emplace_back(round_bound(pts[new_pi]));
          // We insert it on pl before the last vertex
          pl.insert(pl.end()-1, new_pi);
#ifdef DOUBLE_2D_SNAP_FULL_VERBOSE
          std::cout << "Create point " << new_pi << " on " << li << " due to proximity with " << pi << "_____________________________" << std::endl;
          std::cout << new_pi <<": " << pts[new_pi] << std::endl;
          std::cout << li << ":";
          for(std::size_t i: pl)
            std::cout << " " << i;
          std::cout << std::endl;
#endif
          return true;
        }
        return false;
      };

      if(event.pi != event_queue[i-1].pi)
      {
        // Look above segments and creates a point if there too close for a safe rounding
        auto above = pos_it;
        size_t pi=event.pi;
        while(above!=y_order.end() && close_event(pi,*above)){
          if((pi!=polylines[*above].front()) && (pi!=polylines[*above].back()))
            pi=pts.size()-1;
          ++above;
        }

        // same with below segments
        auto below = pos_it;
        pi=event.pi;
        if(below!=y_order.begin()){
          --below;
          while(close_event(pi,*(below))){
            if(below==y_order.begin())
              break;
            if((pi!=polylines[*below].front()) && (pi!=polylines[*below].back()))
              pi=pts.size()-1;
            --below;
          }
        }
      }
    }

    if(event.type==EVENT_TYPE::INSERT){
      y_order.insert(pos_it, event.li);
    }
  }

  // We sort the point by y to ensure the order of the point is preserved
  // Round value of y coordinates are indirectly modified in the case of a filter failure
  auto sort_pi=[&](std::size_t i, std::size_t j){
    return less_y_2(pts[i],pts[j]);
  };
  std::vector< std::size_t > indices(pts.size(),0);
  std::iota(indices.begin(),indices.end(),0);
  std::sort(indices.begin(),indices.end(),sort_pi);

#ifdef DOUBLE_2D_SNAP_FULL_VERBOSE
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
}

template <class Concurrency_tag=Sequential_tag, class Traits, class PointsRange , class PolylineRange>
void merge_duplicate_points_in_polylines(PointsRange &pts, PolylineRange &polylines, const Traits &traits)
{
  using Point_2 = typename Traits::Point_2;
  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename PolylineRange::iterator>::value_type>;

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
template <class Concurrency_tag=Sequential_tag, class Traits, class PointsRange , class PolylineRange>
void snap_post_process(PointsRange &pts, PolylineRange &polylines, const Traits &traits)
{
  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename PolylineRange::iterator>::value_type>;
  using Less_xy_2 = typename Traits::Less_xy_2;

  Less_xy_2 less_xy_2 = traits.less_xy_2_object();

  auto Less_indexes_xy_2=[&](std::size_t i, std::size_t j){
    return less_xy_2(pts[i], pts[j]);
  };

  std::set<std::size_t, decltype(Less_indexes_xy_2)> p_sort_by_x(Less_indexes_xy_2);
  for(std::size_t i=0; i!=pts.size(); ++i)
    p_sort_by_x.insert(i);

  using Iterator_set_x = typename std::set<std::size_t, decltype(Less_indexes_xy_2)>::iterator;

  for(Polyline &poly: polylines){
    std::vector<std::size_t> updated_polyline;
    updated_polyline.push_back(poly.front());
    for(std::size_t i=1; i!=poly.size(); ++i){
      if(pts[poly[i-1]].x()==pts[poly[i]].x()){
        Iterator_set_x start, end;
        // Get all vertices between the two endpoints along x order
        if(Less_indexes_xy_2(poly[i-1],poly[i])){
          start=p_sort_by_x.upper_bound(poly[i-1]);
          end=p_sort_by_x.lower_bound(poly[i]);
          for(auto it=start; it!=end; ++it){
            updated_polyline.push_back(*it);
          }
        } else {
          start=p_sort_by_x.upper_bound(poly[i]);
          end=p_sort_by_x.lower_bound(poly[i-1]);
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

template <class Concurrency_tag=Sequential_tag, class Traits, class PointsRange , class PolylineRange>
void double_snap_rounding_2_disjoint(PointsRange &pts, PolylineRange &polylines, const Traits &traits)
{
  using Point_2   = typename Traits::Point_2;
  using Construct_round_point_2    = typename Traits::Construct_round_point_2;

  Construct_round_point_2 round = traits.construct_round_point_2_object();

  snap_rounding_scan(pts, polylines, traits);
#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Round" << std::endl;
#endif
  for(Point_2 &p: pts)
    p=round(p);
#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Merging points that are collapsed together" << std::endl;
#endif
  merge_duplicate_points_in_polylines(pts, polylines, traits);
#ifdef DOUBLE_2D_SNAP_VERBOSE
  // The algorithm prevents the a vertex that goes through a segment but a vertex may lie on an horizontal/vertical segments after rounding
  std::cout << "Subdivide vertical segments with vertices on them" << std::endl;
#endif
  snap_post_process(pts, polylines, traits);
}

} // end of namespace internal

/**
* ingroup PkgSnapRounding2Ref
*
* Given a range of segments, compute rounded subsegments that are pairwise disjoint in their interiors, based on the input curves.
* The output is a range of polyline with each polyline corresponding to an input segment.
*
* TODO Currently compute_subcurves have no visitor to obtain a polyline per input segment, thus
* currently the output contain one polyline per ubsegment of the arrangement
*
* @tparam Concurrency_tag allows choosing whether the algorithm runs in parallel or sequentially.
* If `CGAL::Parallel_tag` is specified and CGAL is linked with the Intel TBB library, the algorithm will run in parallel.
* Otherwise, if `CGAL::Sequential_tag` is specified (the default), the algorithm will run sequentially.
* @tparam InputIterator iterator over a range of `Segment_2`
* @tparam OutputContainer inserter over a range of `Polyline`. `Polyline` must be a type that provides a `push_back(Point_2)` function.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
* \param begin,end the input segment range
* \param out the output container
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
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
*/
template <class InputIterator , class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
typename OutputContainer::iterator double_snap_rounding_2(InputIterator  	begin,
		                                                      InputIterator  	end,
		                                                      OutputContainer &out,
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
  using Segment_2 = typename Traits::Segment_2;
  using I2E = typename Traits::Converter_to_exact;
  using E2O = typename Traits::Converter_from_exact;
  using VectorIterator = typename std::vector<Segment_2>::iterator;

  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename OutputContainer::iterator>::value_type>;

  using Less_xy_2 = typename Traits::Less_xy_2;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const Traits &traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Traits());

  I2E to_exact=traits.converter_to_exact_object();
  E2O from_exact=traits.converter_from_exact_object();

  std::vector<Segment_2> convert_input;
  for(InputIterator it=begin; it!=end; ++it)
    convert_input.push_back(to_exact(*it));
  std::vector<Segment_2> segs;
#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Solved intersections" << std::endl;
  std::cout << "do intersect? " << do_curves_intersect(convert_input.begin(), convert_input.end()) << std::endl;
#endif
  compute_subcurves(convert_input.begin(), convert_input.end(), std::back_inserter(segs));

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Change format to range of points and indexes" << std::endl;
#endif
  std::set<Point_2> unique_point_set;
  std::map<Point_2, int> point_to_index;
  std::vector<Point_2> pts;
  std::vector< std::vector< std::size_t> > polylines;

  // Transform range of the segments in the range of points and polyline of indexes
  for(VectorIterator it=segs.begin(); it!=segs.end(); ++it)
  {
    const Point_2& p1 = it->source();
    const Point_2& p2 = it->target();

    if (unique_point_set.find(p1) == unique_point_set.end()) {
      unique_point_set.insert(p1);
      pts.push_back(p1);
      point_to_index[p1] = pts.size() - 1;
    }
    if (unique_point_set.find(p2) == unique_point_set.end()) {
      unique_point_set.insert(p2);
      pts.push_back(p2);
      point_to_index[p2] = pts.size() - 1;
    }
  }

  for(VectorIterator it=segs.begin(); it!=segs.end(); ++it)
  {
    std::size_t index1 = point_to_index[it->source()];
    std::size_t index2 = point_to_index[it->target()];
    if(Less_xy_2()(it->source(), it->target()))
      polylines.push_back({index1, index2});
    else
      polylines.push_back({index2, index1});
  }

  // Main algorithm
  internal::double_snap_rounding_2_disjoint<Concurrency_tag, Traits>(pts, polylines, traits);

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Build output" << std::endl;
#endif

  // Output polylines
  for(auto &poly: polylines){
    Polyline new_line;
    for(std::size_t pi: poly)
      new_line.push_back(from_exact(pts[pi]));
    out.push_back(new_line);
  }

  return out.end();
}

/**
* ingroup PkgSnapRounding2Ref
*
* Given a range of segments, compute rounded subsegments that are pairwise disjoint in their interior, as induced by the input curves.
*
* @tparam Concurrency_tag That template parameter enables to choose whether the algorithm is to be run in
* parallel, if CGAL::Parallel_tag is specified and CGAL has been linked with the Intel TBB library, or sequentially, if CGAL::Sequential_tag - the default value - is specified.
* @tparam InputIterator iterator of a segment range
* @tparam OutputContainer inserter of a segment range
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
* \param begin,end the input segment range
* \param out the output inserter
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
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
*/
template <class InputIterator , class OutputIterator, class NamedParameters = parameters::Default_named_parameters>
OutputIterator compute_snapped_subcurves_2(InputIterator  	 begin,
		                                       InputIterator  	 end,
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
  using Segment_2 = typename Traits::Segment_2;
  using I2E = typename Traits::Converter_to_exact;
  using E2O = typename Traits::Converter_from_exact;

  using SegmentRange = std::vector<Segment_2>;
  using SegmentRangeIterator = typename SegmentRange::iterator;

  using Less_xy_2 = typename Traits::Less_xy_2;
  using Construct_segment_2 = typename Traits::Construct_segment_2;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const Traits &traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Traits());

  I2E to_exact=traits.converter_to_exact_object();
  E2O from_exact=traits.converter_from_exact_object();

  Less_xy_2 less_xy_2 = traits.less_xy_2_object();
  Construct_segment_2 segment_2 = traits.construct_segment_2_object();

  std::vector<Segment_2> convert_input;
  for(InputIterator it=begin; it!=end; ++it)
    convert_input.push_back(to_exact(*it));
  std::vector<Segment_2> segs;
#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Solved intersections" << std::endl;
  std::cout << "do intersect? " << do_curves_intersect(convert_input.begin(), convert_input.end()) << std::endl;
#endif
  compute_subcurves(convert_input.begin(), convert_input.end(), std::back_inserter(segs));

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Change format to range of points and indexes" << std::endl;
#endif
  std::set<Point_2, Less_xy_2> unique_point_set(less_xy_2);
  std::map<Point_2, std::size_t> point_to_index;
  std::vector<Point_2> pts;
  std::vector< std::vector< std::size_t> > polylines;

  // Transform range of the segments in the range of points and polyline of indexes
  for(SegmentRangeIterator it=segs.begin(); it!=segs.end(); ++it)
  {
    const Point_2& p1 = it->source();
    const Point_2& p2 = it->target();

    if (unique_point_set.find(p1) == unique_point_set.end()) {
      unique_point_set.insert(p1);
      pts.push_back(p1);
      point_to_index[p1] = pts.size() - 1;
    }
    if (unique_point_set.find(p2) == unique_point_set.end()) {
      unique_point_set.insert(p2);
      pts.push_back(p2);
      point_to_index[p2] = pts.size() - 1;
    }
  }

  for(SegmentRangeIterator it=segs.begin(); it!=segs.end(); ++it)
  {
    std::size_t index1 = point_to_index[it->source()];
    std::size_t index2 = point_to_index[it->target()];
    if(less_xy_2(it->source(), it->target()))
      polylines.push_back({index1, index2});
    else
      polylines.push_back({index2, index1});
  }


  // Main algorithm
  internal::double_snap_rounding_2_disjoint<Concurrency_tag>(pts, polylines, traits);

#ifdef DOUBLE_2D_SNAP_VERBOSE
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
* ingroup PkgSnapRounding2Ref
*
* Given a range of `Polygon_2`, compute rounded polygons such that their segments are either equal either disjoint in their interior, as induced by the input polygons.
* The polygons are intended to be non-intersecting, unless the named parameter `compute_intersection` is set to `true`.
* Any polygon is guarantee to remain a Polygon in the output but may present pinched section or/and common vertices or segments with
* other polygons.
*
* @tparam InputIterator iterator of a CGAL::Polygon_2 range
* @tparam OutputContainer inserter of a CGAL::Polygon_2 range
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
* \param begin,end the input polygon range
* \param out the output inserter
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* \cgalNamedParamsBegin
*   \cgalParamNBegin{concurrency_tag}
*     \cgalParamDescription{That template parameter enables to choose whether the algorithm is to be run in parallel, if CGAL::Parallel_tag
*                           is specified and CGAL has been linked with the Intel TBB library, or sequentially, otherwise.}
*     \cgalParamType{CGAL::Concurrency_tag}
*     \cgalParamDefault{CGAL::Sequential_tag}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{compute_intersection}
*     \cgalParamDescription{Enable intersection computation between the polygons before performing snapping.}
*     \cgalParamType{boolean}
*     \cgalParamDefault{false}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `FloatSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `Float_snap_rounding_traits_2`}
*   \cgalParamNEnd
* @warning The convex property of the polygons is not necessarly preserved
*/
template <class InputIterator, class OutputIterator, class NamedParameters = parameters::Default_named_parameters>
void compute_snapped_polygons_2(InputIterator begin,
                                InputIterator end,
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
  using I2E = typename Traits::Converter_to_exact;
  using E2O = typename Traits::Converter_from_exact;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  Polygon_2 P=*begin;

  const Traits &traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Traits());
  const bool compute_intersections = choose_parameter(get_parameter(np, internal_np::compute_intersection), false);

  I2E to_exact=traits.converter_to_exact_object();
  E2O from_exact=traits.converter_from_exact_object();

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Change format to range of points and indexes" << std::endl;
#endif
  std::vector<Point_2> pts;
  std::vector< std::vector< std::size_t> > polylines;
  std::vector<size_t> polygon_index;

  if(compute_intersections){
    // TODO Need to compute_subcurves that output polylines to track the polygons
    assert(0);
  } else {
    // Index of the segments that introduced a new polygon
    polygon_index.reserve(std::distance(begin, end));
    for(InputIterator it=begin; it!=end; ++it){
      size_t index_start=polylines.size();
      polygon_index.push_back(polylines.size());
      for(const typename Polygon_2::Point_2 &p: it->vertices())
        pts.push_back(to_exact(p));
      for(size_t i=0; i<P.size()-1; ++i)
        polylines.push_back({i, i+1});
      polylines.push_back({pts.size()-1,index_start});
      assert(pts.size()==polylines.size());
    }
    polygon_index.push_back(polylines.size());
  }

  // Main algorithm
  internal::double_snap_rounding_2_disjoint<Concurrency_tag>(pts, polylines, traits);

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Build output" << std::endl;
#endif

  // Output a range of segments while removing duplicate ones
  for(std::size_t input_ind=0; input_ind<std::size_t(std::distance(begin,end)); ++input_ind){
    for(std::size_t pl_ind=polygon_index[input_ind]; pl_ind<polygon_index[input_ind+1]; ++pl_ind){
      std::vector<std::size_t> &poly = polylines[pl_ind];
      Polygon_2 P;
      for(std::size_t i=1; i<poly.size(); ++i)
        P.push_back(from_exact(pts[poly[i]]));
    }
    *out++=P;
  }
}


/**
* ingroup
*
* Given a Polygon_2, compute rounded segments that are pairwise disjoint in their interior, as induced by the input polygon.
* The output is guarantee to be a Polygon but may present pinched section.
*
* \tparam Polygon_2 model of CGAL::Polygon_2
* \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
* \param P the input polygon
* \param out the output polygon
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
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
* @warning The convex property is not necessarly preserved
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
