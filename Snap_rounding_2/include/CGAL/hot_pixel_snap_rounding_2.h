// Copyright (c) 2001  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// author(s)     : Eli Packer <elip@post.tau.ac.il>

#ifndef CGAL_HOT_PIXEL_SNAP_ROUNDING_2_H
#define CGAL_HOT_PIXEL_SNAP_ROUNDING_2_H

#include <CGAL/license/Snap_rounding_2.h>

#include <algorithm>
#include <cstddef>
#include <deque>
#include <iostream>
#include <iterator>
#include <list>
#include <memory>
#include <set>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <CGAL/enum.h>
#include <CGAL/Hot_pixel_snap_rounding_traits_2.h>
#include <CGAL/intersection_2.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/utility.h>
#include <CGAL/internal/Snap_rounding_helpers.h>
#include <CGAL/internal/Snap_rounding_kd_2.h>

namespace CGAL {

namespace SR_internal = Snap_Rounding_2::internal;

namespace Snap_Rounding_2::internal {

/*! the direction of a segment from its source to its target.
 */
enum SEG_Direction {
  SEG_UP_RIGHT, SEG_UP_LEFT, SEG_DOWN_RIGHT, SEG_DOWN_LEFT, SEG_UP, SEG_DOWN, SEG_LEFT, SEG_RIGHT, SEG_POINT_SEG
};

#ifdef CGAL_SR_DEBUG
inline int number_of_false_hp = 0;      // the number of rejected kd-tree candidates
#endif

/*! a segment represented by its two endpoints.
 */
template <typename Traits_>
class Segment_data {
private:
  using Traits = Traits_;
  using Segment_2 = typename Traits::Segment_2;
  using Point_2 = typename Traits::Point_2;

  Point_2 m_p;
  Point_2 m_q;

public:
  /*! constructs from two endpoints.
   */
  Segment_data(const Point_2& p, const Point_2& q) : m_p(p), m_q(q) {}

  /*! obtains the segment.
   */
  Segment_2 segment() const { return Segment_2(m_p, m_q); }

  /*! obtains the source.
   */
  const Point_2& source() const { return m_p; }

  /*! obtains the target.
   */
  const Point_2& target() const { return m_q; }

  /*! obtains the direction from the source to the target.
   * \param gt the traits.
   */
  SEG_Direction direction(const Traits& gt) const;
};

/*! a hot pixel, that is, a pixel that contains a vertex of the arrangement of the input segments.
 */
template <typename Traits_>
class Hot_pixel {
private:
  using Traits = Traits_;
  using NT = typename Traits::FT;
  using Segment_2 = typename Traits::Segment_2;
  using Point_2 = typename Traits::Point_2;
  using Segment_data = SR_internal::Segment_data<Traits>;

  Point_2 m_p;                  // the center of the hot pixel
  Point_2 m_p_left;
  Point_2 m_p_right;
  Point_2 m_p_down;
  Point_2 m_p_up;
  NT m_pixel_size = 0;
  Segment_2 m_right_seg;
  Segment_2 m_left_seg;
  Segment_2 m_top_seg;
  Segment_2 m_bot_seg;

public:
  /*! constructs the hot pixel that contains a point.
   * \param inp_point the point.
   * \param inp_pixel_size the pixel size.
   * \param gt the traits.
   */
  Hot_pixel(const Point_2& inp_point, NT inp_pixel_size, const Traits& gt);

  /*! obtains the center.
   */
  const Point_2& get_center() const { return m_p; }

  /*! obtains the center.
   * \param int_output indicates whether to express the center in the integer grid.
   * \param gt the traits.
   */
  Point_2 get_center(bool int_output, const Traits& gt) const;

  /*! determines whether a segment intersects the left edge of the pixel, following the boundary
   * rules of `intersect()`. `intersect_right()`, `intersect_bot()`, and `intersect_top()` handle the
   * other edges.
   * \param seg the segment.
   * \param seg_dir the direction of `seg`.
   * \param gt the traits.
   */
  bool intersect_left(const Segment_2& seg, SEG_Direction seg_dir, const Traits& gt) const;
  bool intersect_right(const Segment_2& seg, SEG_Direction seg_dir, const Traits& gt) const;
  bool intersect_bot(const Segment_2& seg, SEG_Direction seg_dir, const Traits& gt) const;
  bool intersect_top(const Segment_2& seg, SEG_Direction seg_dir, const Traits& gt) const;

  /*! determines whether a segment intersects the hot pixel. Depending on the direction of the
   * segment, some contacts with the boundary of the pixel do not count as intersections.
   * \param seg the segment.
   * \param seg_dir the direction of `seg`.
   * \param gt the traits.
   */
  bool intersect(const Segment_data& seg, SEG_Direction seg_dir, const Traits& gt) const;
};

/*! a comparator that orders the hot pixels a segment intersects along the direction of the segment.
 */
template <typename Traits_>
class Hot_pixel_dir_cmp {
private:
  using Hot_pixel = SR_internal::Hot_pixel<Traits_>;

  const Traits_* m_gt;
  SEG_Direction m_dir;

public:
  /*! constructs from the direction of the segment.
   * \param dir the direction of the segment.
   * \param gt the traits, which must outlive the comparator.
   */
  Hot_pixel_dir_cmp(SEG_Direction dir, const Traits_& gt) : m_gt(&gt), m_dir(dir) {}

  /*! determines whether `h1` precedes `h2` along the segment.
   */
  bool operator()(const Hot_pixel* h1, const Hot_pixel* h2) const;
};

/*! the hot-pixel snap-rounding algorithm.
 */
template <typename Traits_, typename OutputContainer>
class Hot_pixel_snap_rounding_2 {
private:
  using Traits = Traits_;
  using NT = typename Traits::FT;
  using Segment_2 = typename Traits::Segment_2;
  using Point_2 = typename Traits::Point_2;
  using Polyline_type = typename OutputContainer::value_type;
  using Hot_pixel = SR_internal::Hot_pixel<Traits>;
  using Segment_data = SR_internal::Segment_data<Traits>;
  using Hot_pixel_kd_tree = SR_internal::Hot_pixel_kd_tree<Traits, Hot_pixel*>;
  using Hot_pixel_dir_cmp = SR_internal::Hot_pixel_dir_cmp<Traits>;
  using Hot_pixel_set = std::set<Hot_pixel*, Hot_pixel_dir_cmp>;

public:
  using Segment_data_list = std::vector<Segment_data>;

  /*! constructs from the traits, which is copied.
   */
  explicit Hot_pixel_snap_rounding_2(const Traits& gt) : m_gt(gt) {}

  /*! computes the hot pixels, that is, the pixels that contain the endpoints of the segments and
   * their intersection points, and builds the kd-tree of the hot pixels.
   * \param pixel_size the pixel size.
   * \param seg_list the input segments.
   */
  void find_hot_pixels_and_create_kd_tree(NT pixel_size, const Segment_data_list& seg_list);

  /*! reroutes every segment through the centers of the hot pixels it intersects.
   * \param output_container the output; one polyline per input segment.
   * \param pixel_size the pixel size.
   * \param int_output indicates whether to express the output in the integer grid.
   * \param do_isr indicates whether to apply iterated snap rounding.
   * \param seg_list the input segments.
   * \pre `find_hot_pixels_and_create_kd_tree()` was called with `seg_list`.
   */
  void iterate(OutputContainer& output_container, NT pixel_size, bool int_output, bool do_isr,
               const Segment_data_list& seg_list);

private:
  Traits m_gt;
  std::deque<Hot_pixel> m_hot_pixels;
  std::unique_ptr<Hot_pixel_kd_tree> m_kd_tree;
  std::vector<Hot_pixel*> m_candidates;

  /*! obtains the hot pixels a segment intersects, ordered along the segment.
   */
  Hot_pixel_set find_intersected_hot_pixels(const Segment_data& seg);

  /*! appends to `seg_output` the centers of the hot pixels in `hot_pixels` but the first.
   */
  void reroute_sr(const Hot_pixel_set& hot_pixels, Polyline_type& seg_output, bool int_output);

  /*! appends to `seg_output` the centers of the hot pixels in `hot_pixels` but the first, where
   * the segment between every two consecutive centers is rerouted recursively.
   * \param first_time indicates whether `hot_pixels` are the hot pixels of an input segment.
   */
  void reroute_isr(const Hot_pixel_set& hot_pixels, Polyline_type& seg_output, bool first_time, bool int_output);
};

template <typename Traits_>
SEG_Direction Segment_data<Traits_>::direction(const Traits_& gt) const {
  auto compare_x = gt.compare_x_2_object();
  auto compare_y = gt.compare_y_2_object();

  Comparison_result cx = compare_x(m_p, m_q);
  Comparison_result cy = compare_y(m_p, m_q);

  if (cx == SMALLER) {
    if (cy == SMALLER) return SEG_UP_RIGHT;
    if (cy == EQUAL) return SEG_RIGHT;
    return SEG_DOWN_RIGHT;
  }
  if (cx == EQUAL) {
    if (cy == SMALLER) return SEG_UP;
    if (cy == EQUAL) return SEG_POINT_SEG;
    return SEG_DOWN;
  }
  if (cy == SMALLER) return SEG_UP_LEFT;
  if (cy == EQUAL) return SEG_LEFT;
  return SEG_DOWN_LEFT;
}

template <typename Traits_>
Hot_pixel<Traits_>::Hot_pixel(const Point_2& inp_point, NT inp_pixel_size, const Traits_& gt) :
  m_pixel_size(inp_pixel_size) {
  NT x, y;
  gt.snap_2_object()(inp_point, m_pixel_size, x, y);

  const NT half_pixel = m_pixel_size / NT(2);
  NT left_coord = x - half_pixel;
  NT right_coord = x + half_pixel;
  NT bottom_coord = y - half_pixel;
  NT top_coord = y + half_pixel;
  m_p = Point_2(x, y);
  m_p_left = Point_2(left_coord, y);
  m_p_right = Point_2(right_coord, y);
  m_p_down = Point_2(x, bottom_coord);
  m_p_up = Point_2(x, top_coord);

  auto construct_seg = gt.construct_segment_2_object();

  Point_2 lb(left_coord, bottom_coord);
  Point_2 rb(right_coord, bottom_coord);
  Point_2 lt(left_coord, top_coord);
  Point_2 rt(right_coord, top_coord);
  m_right_seg = construct_seg(rb, rt);
  m_left_seg = construct_seg(lb, lt);
  m_top_seg = construct_seg(lt, rt);
  m_bot_seg = construct_seg(lb, rb);
}

template <typename Traits_>
typename Hot_pixel<Traits_>::Point_2 Hot_pixel<Traits_>::get_center(bool int_output, const Traits_& gt) const {
  return (int_output) ? gt.integer_grid_point_2_object()(m_p, m_pixel_size) : m_p;
}

template <typename Traits_>
bool Hot_pixel<Traits_>::intersect_left(const Segment_2& seg, SEG_Direction seg_dir, const Traits_& gt) const {
  const auto result = intersection(seg, m_left_seg);
  if (! result) return false;
  const Point_2* p = std::get_if<Point_2>(&*result);
  if (p == nullptr) return true; // the segments overlap

  auto compare_y = gt.compare_y_2_object();
  auto construct_vertex = gt.construct_vertex_2_object();

  Comparison_result c_p = compare_y(*p, m_p_up);
  Comparison_result c_target = compare_y(construct_vertex(seg, 1), m_p_up);
  Comparison_result c_source = compare_y(construct_vertex(seg, 0), m_p_up);

  return (c_p != EQUAL || (seg_dir == SEG_UP_LEFT && c_source != EQUAL) ||
          (seg_dir == SEG_DOWN_RIGHT && c_target != EQUAL));
}

template <typename Traits_>
bool Hot_pixel<Traits_>::intersect_right(const Segment_2& seg, SEG_Direction seg_dir, const Traits_& gt) const {
  const auto result = intersection(seg, m_right_seg);
  if (! result) return false;
  const Point_2* p = std::get_if<Point_2>(&*result);
  if (p == nullptr) return false; // the segments overlap

  // bottom right point was checked in intersect_bot
  auto compare_y = gt.compare_y_2_object();
  auto construct_vertex = gt.construct_vertex_2_object();

  Comparison_result c1 = compare_y(*p, m_p_up);
  Comparison_result c2 = compare_y(*p, m_p_down);

  Point_2 src = construct_vertex(seg, 0);
  Point_2 trg = construct_vertex(seg, 1);
  Comparison_result c3 = compare_y(src, m_p_up);
  Comparison_result c4 = compare_y(trg, m_p_up);

  if (c1 == EQUAL)
    return ((seg_dir == SEG_UP_RIGHT && c3 != EQUAL) || (seg_dir == SEG_DOWN_LEFT && c4 != EQUAL));
  else if (c2 == EQUAL)
    return false; // was checked

  auto compare_x = gt.compare_x_2_object();
  Comparison_result c_target = compare_x(m_p_right, trg);
  Comparison_result c_source = compare_x(m_p_right, src);

  return (((seg_dir == SEG_LEFT || seg_dir == SEG_DOWN_LEFT || seg_dir == SEG_UP_LEFT) && c_target != EQUAL) ||
          ((seg_dir == SEG_RIGHT || seg_dir == SEG_DOWN_RIGHT || seg_dir == SEG_UP_RIGHT) && c_source != EQUAL));
}

template <typename Traits_>
bool Hot_pixel<Traits_>::intersect_bot(const Segment_2& seg, SEG_Direction seg_dir, const Traits_& gt) const {
  const auto result = intersection(seg, m_bot_seg);
  if (! result) return false;
  const Point_2* p = std::get_if<Point_2>(&*result);
  if (p == nullptr) return true; // the segments overlap

  auto compare_x = gt.compare_x_2_object();
  auto construct_vertex = gt.construct_vertex_2_object();

  Comparison_result c_p = compare_x(*p, m_p_right);
  Comparison_result c_target = compare_x(construct_vertex(seg, 1), m_p_right);
  Comparison_result c_source = compare_x(construct_vertex(seg, 0), m_p_right);

  return (c_p != EQUAL || (seg_dir == SEG_UP_LEFT && c_target != EQUAL) ||
          (seg_dir == SEG_DOWN_RIGHT && c_source != EQUAL));
}

template <typename Traits_>
bool Hot_pixel<Traits_>::intersect_top(const Segment_2& seg, SEG_Direction seg_dir, const Traits_& gt) const {
  const auto result = intersection(seg, m_top_seg);
  if (! result) return false;
  const Point_2* p = std::get_if<Point_2>(&*result);
  if (p == nullptr) return false; // the segments overlap

  auto compare_x = gt.compare_x_2_object();
  auto compare_y = gt.compare_y_2_object();
  auto construct_vertex = gt.construct_vertex_2_object();

  // corner points was checked in intersect_bot
  Comparison_result c1 = compare_x(*p, m_p_left);
  Comparison_result c2 = compare_x(*p, m_p_right);
  Comparison_result c3 = compare_y(construct_vertex(seg, 1), m_p_up);
  Comparison_result c4 = compare_y(construct_vertex(seg, 0), m_p_up);

  if (c1 == EQUAL || c2 == EQUAL)
    return false; // were checked
  else
    return (((seg_dir == SEG_DOWN || seg_dir == SEG_DOWN_LEFT || seg_dir == SEG_DOWN_RIGHT) && c3 != EQUAL) ||
            ((seg_dir == SEG_UP || seg_dir == SEG_UP_LEFT || seg_dir == SEG_UP_RIGHT) && c4 != EQUAL));
}

template <typename Traits_>
bool Hot_pixel<Traits_>::intersect(const Segment_data& seg, SEG_Direction seg_dir, const Traits_& gt) const {
  Segment_2 s = seg.segment();
  return (intersect_bot(s, seg_dir, gt) || intersect_left(s, seg_dir, gt) ||
          intersect_right(s, seg_dir, gt) || intersect_top(s, seg_dir, gt));
}

template <typename Traits_>
bool Hot_pixel_dir_cmp<Traits_>::operator()(const Hot_pixel* h1, const Hot_pixel* h2) const {
  auto compare_x = m_gt->compare_x_2_object();
  auto compare_y = m_gt->compare_y_2_object();

  Comparison_result cx = compare_x(h1->get_center(), h2->get_center());
  Comparison_result cy = compare_y(h1->get_center(), h2->get_center());

  // Point segment intersects only one pixel, thus ignored
  return ((m_dir == SEG_UP_RIGHT && (cx == SMALLER || (cx == EQUAL && cy == SMALLER))) ||
          (m_dir == SEG_UP_LEFT && (cx == LARGER || (cx == EQUAL && cy == SMALLER))) ||
          (m_dir == SEG_DOWN_RIGHT && (cx == SMALLER || (cx == EQUAL && cy == LARGER))) ||
          (m_dir == SEG_DOWN_LEFT && (cx == LARGER || (cx == EQUAL && cy == LARGER))) ||
          (m_dir == SEG_UP && cy == SMALLER) ||
          (m_dir == SEG_DOWN && cy == LARGER) ||
          (m_dir == SEG_LEFT && cx == LARGER) ||
          (m_dir == SEG_RIGHT && cx == SMALLER));
}

template <typename Traits_, typename OutputContainer>
void Hot_pixel_snap_rounding_2<Traits_, OutputContainer>::
find_hot_pixels_and_create_kd_tree(NT pixel_size, const Segment_data_list& seg_list) {
  if (seg_list.empty()) return;

  auto construct_seg = m_gt.construct_segment_2_object();
  std::vector<Segment_2> segments;
  segments.reserve(seg_list.size());
  for (const auto& seg : seg_list) segments.push_back(construct_seg(seg.source(), seg.target()));

  // get intersection points (with endpoints)
  std::list<Point_2> points;
  CGAL::compute_intersection_points(segments.begin(), segments.end(), std::back_inserter(points), true, m_gt);

  typename Hot_pixel_kd_tree::Point_saved_pair_list centers;
  centers.reserve(points.size());
  for (const auto& p : points) {
    Hot_pixel& hp = m_hot_pixels.emplace_back(p, pixel_size, m_gt);
    centers.emplace_back(hp.get_center(), &hp);
  }

  m_kd_tree = std::make_unique<Hot_pixel_kd_tree>(centers, pixel_size);
}

template <typename Traits_, typename OutputContainer>
typename Hot_pixel_snap_rounding_2<Traits_, OutputContainer>::Hot_pixel_set
Hot_pixel_snap_rounding_2<Traits_, OutputContainer>::find_intersected_hot_pixels(const Segment_data& seg) {
  const SEG_Direction seg_dir = seg.direction(m_gt);
  Hot_pixel_set hot_pixels(Hot_pixel_dir_cmp{seg_dir, m_gt});

  m_candidates.clear();
  m_kd_tree->get_intersecting_points(m_candidates, seg.segment(), m_gt);
  for (Hot_pixel* hp : m_candidates) {
    if (hp->intersect(seg, seg_dir, m_gt)) hot_pixels.insert(hp);
#ifdef CGAL_SR_DEBUG
    else ++number_of_false_hp;
#endif
  }
  return hot_pixels;
}

template <typename Traits_, typename OutputContainer>
void Hot_pixel_snap_rounding_2<Traits_, OutputContainer>::reroute_sr(const Hot_pixel_set& hot_pixels,
                                                                     Polyline_type& seg_output,
                                                                     bool int_output) {
  for (auto it = std::next(hot_pixels.begin()); it != hot_pixels.end(); ++it)
    seg_output.push_back((*it)->get_center(int_output, m_gt));
}

template <typename Traits_, typename OutputContainer>
void Hot_pixel_snap_rounding_2<Traits_, OutputContainer>::reroute_isr(const Hot_pixel_set& hot_pixels,
                                                                      Polyline_type& seg_output,
                                                                      bool first_time,
                                                                      bool int_output) {
  if (first_time || (hot_pixels.size() > 2)) {
    const auto last = std::prev(hot_pixels.end());
    for (auto it = hot_pixels.begin(); it != last; ++it) {
      Segment_data seg((*it)->get_center(), (*std::next(it))->get_center());
      reroute_isr(find_intersected_hot_pixels(seg), seg_output, false, int_output);
    }
  }
  else {
    // insert second hot pixel
    seg_output.push_back((*std::next(hot_pixels.begin()))->get_center(int_output, m_gt));
  }
}

template <typename Traits_, typename OutputContainer>
void Hot_pixel_snap_rounding_2<Traits_, OutputContainer>::iterate(OutputContainer& output_container,
                                                                  NT pixel_size,
                                                                  bool int_output, bool do_isr,
                                                                  const Segment_data_list& seg_list) {
  Polyline_type seg_output;
  for (const auto& seg : seg_list) {
    seg_output.clear();
    const Hot_pixel_set hot_pixels = find_intersected_hot_pixels(seg);
    if (hot_pixels.empty()) {
      // segment entirely inside a pixel
      seg_output.push_back(Hot_pixel(seg.source(), pixel_size, m_gt).get_center(int_output, m_gt));
    }
    else {
      seg_output.push_back((*hot_pixels.begin())->get_center(int_output, m_gt));
      // segments that intersect a single hot pixel are done
      if (hot_pixels.size() > 1) {
        if (do_isr) reroute_isr(hot_pixels, seg_output, true, int_output);
        else reroute_sr(hot_pixels, seg_output, int_output);
      }
    }
    output_container.push_back(std::move(seg_output));
  }
}

/*! snap rounds a range of segments; one polyline per input segment.
 * \param gt the traits; its type determines `Traits`.
 */
template <typename Traits, typename InputIterator, typename OutputContainer>
void hot_pixel_snap_rounding_2(InputIterator begin,
                               InputIterator end,
                               OutputContainer& output_container,
                               const Traits& gt,
                               typename Traits::NT pixel_size,
                               bool do_isr = true,
                               bool int_output = true) {
#ifdef CGAL_SR_DEBUG
  number_of_false_hp = 0;
#endif

  using Snap_rounding = Hot_pixel_snap_rounding_2<Traits, OutputContainer>;

  typename Snap_rounding::Segment_data_list seg_list;
  for (; begin != end; ++begin) seg_list.emplace_back(begin->source(), begin->target());

  output_container.clear();
  Snap_rounding sr(gt);
  sr.find_hot_pixels_and_create_kd_tree(pixel_size, seg_list);
  sr.iterate(output_container, pixel_size, int_output, do_isr, seg_list);

#ifdef CGAL_SR_DEBUG
  std::cout << "Overall number of false hot pixels in all the queries : " << number_of_false_hp << std::endl;
#endif
}

/*! snap rounds a range of polygons.
 */
template <typename PolygonRange, typename OutputContainer,
          typename NamedParameters = parameters::Default_named_parameters>
void hot_pixel_snap_rounding_2_polygon(const PolygonRange& polygons, OutputContainer& out,
                                       const NamedParameters& np = parameters::default_values()) {
  using Polygon_2 = typename std::iterator_traits<typename PolygonRange::iterator>::value_type;
  using Kernel = typename Kernel_traits<typename Polygon_2::Point_2>::Kernel;
  using DefaultTraits = Hot_pixel_snap_rounding_traits_2<Kernel>;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const auto gt = choose_parameter<DefaultTraits>(get_parameter(np, internal_np::geom_traits));
  using Traits = std::decay_t<decltype(gt)>;
  using Point_2 = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;

  std::vector<Segment_2> input_segments;
  std::vector<std::vector<Point_2>> polylines;

  // Store the indices of segment of a new polygon, segments between
  // [ polygon_index[i] and polygon_index[i+1] [ belong to polygon i
  std::vector< std::size_t > polygon_indices;

  polygon_indices.reserve(std::distance(polygons.begin(), polygons.end()));
  for (const Polygon_2& P: polygons){
    polygon_indices.push_back(input_segments.size());
    for (std::size_t i = 0; i < P.size()-1; ++i)
      input_segments.emplace_back(P[i], P[i+1]);
    input_segments.emplace_back(P[P.size()-1], P[0]);
  }
  polygon_indices.push_back(input_segments.size());

  auto pixel_size = choose_parameter(get_parameter(np, internal_np::pixel_size), 1.);
  bool do_isr = choose_parameter(get_parameter(np, internal_np::do_iterative_snap_rounding), true);
  bool int_output = choose_parameter(get_parameter(np, internal_np::use_grid_coordinates), true);

  // Main algorithm
  hot_pixel_snap_rounding_2(input_segments.begin(), input_segments.end(), polylines, gt,
                            pixel_size, do_isr, int_output);

  // Reassemble the polygons
  for (std::size_t polygon_idx = 0; polygon_idx != polygon_indices.size()-1; ++polygon_idx){
    Polygon_2 P;
    std::size_t idx_start = polygon_indices[polygon_idx];
    std::size_t idx_end = polygon_indices[polygon_idx+1];
    Point_2 last_insert;
    for (std::size_t pl_idx = idx_start; pl_idx != idx_end; ++pl_idx){
      const auto &pl = polylines[pl_idx];
      // The first element is not add, it is identical to the last

      bool go_forward;
      if (pl_idx == idx_start)
        go_forward = (pl.back() == polylines[pl_idx+1].front()) || (pl.back() == polylines[pl_idx+1].back());
      else
        go_forward = (last_insert == pl.front());

      // Add the element in forward direction
      if (go_forward){
        for (std::size_t i = 1; i != pl.size(); ++i)
          P.push_back(pl[i]);
        last_insert = pl.back();

      // Add the element in backward direction
      }
      else {
        for (std::size_t i = pl.size() - 1; i != 0; --i)
          P.push_back(pl[i-1]);
        last_insert = pl.front();
      }
    }
    out.push_back(P);
  }
}

} // namespace internal

#if DOXYGEN_RUNNING

/** \ingroup Snap_rounding_hot_pixel_grp
 *
 * \brief subdivides and rounds a range of segments so that they are pairwise disjoint in their interiors.
 *
 * By default, each polyline of the output corresponds to an input
 * segment. Consequently, duplicate segments may appear in the output, for
 * instance when multiple input segments collapse.  When the parameter
 * `output_unique_segments` is set to `true`, the polylines are decomposed into
 * individual segments (represented as polylines with two points), and duplicates
 * are removed.
 *
 * @tparam SegmentRange model of a `ConstRange` whose iterator is model of `ForwardIterator` and whose value_type
 *                      is `geom_traits::Segment_2`, where the type of geom_traits is detailed by `np::geom_traits`.
 * @tparam OutputContainer model of the concept `BackInsertionSequence` whose value type is itself a model of the
 *                         concepts `DefaultConstructible` and `BackInsertionSequence` whose value type is
 *                         `geom_traits::Point_2`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param segments the input segment range
 * \param out the output container
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{pixel_size}
 *     \cgalParamDescription{The size of the pixel. The plane will be tiled with square pixels of that width such that the origin is the center of a pixel.}
 *     \cgalParamType{`geom_traits::FT`}
 *     \cgalParamDefault{1}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{do_iterative_snap_rounding}
 *     \cgalParamDescription{determines whether to apply Iterative Snap Rounding, see the user manual for more details.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{true}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{use_grid_coordinates}
 *     \cgalParamDescription{If set to true, the output coordinates are expressed in the integer grid (pixel indices).
 *                           Otherwise, they are given in the input coordinate system.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{true}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{output_unique_segments}
 *     \cgalParamDescription{If set to true, the output polylines are unique pairs of distinct points represented a segment. As a result, the total number of output polylines may differ from the number of input segments.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{false}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{The traits class must be a model of `HotPixelSnapRoundingTraits_2`}
 *     \cgalParamDefault{an instance of `CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>` where Kernel is deduced from the segment type, using `CGAL::Kernel_traits`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 */
template <typename SegmentRange, typename OutputContainer,
          typename NamedParameters = parameters::Default_named_parameters>
void hot_pixel_snap_rounding_2(const SegmentRange& segments, OutputContainer& out,
                               const NamedParameters& np = parameters::default_values());

/** \ingroup Snap_rounding_hot_pixel_grp
 *
 * \brief subdivides and rounds a range of polygons so that their boundary
 * segments are pairwise disjoint in their interiors.  If the input polygons are
 * disjoint, the output polygons remain non-overlapping, although they may share
 * vertices or edges.  Each output polygon is free of self-intersections but may
 * present pinched sections.
 *
 * @tparam PolygonRange model of a `ConstRange` whose iterator is model of `ForwardIterator` and whose value_type
 *                      is model of `CGAL::Polygon_2`.
 * @tparam OutputContainer model of the concept `BackInsertionSequence` whose value type is model of `CGAL::Polygon_2`.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param polygons the range of input polygons
 * \param out the output container
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{pixel_size}
 *     \cgalParamDescription{The size of the pixel. The plane will be tiled with square pixels of that width such that the origin is the center of a pixel.}
 *     \cgalParamType{`geom_traits::FT`}
 *     \cgalParamDefault{1}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{do_iterative_snap_rounding}
 *     \cgalParamDescription{determines whether to apply Iterative Snap Rounding, see the user manual for more details.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{true}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{use_grid_coordinates}
 *     \cgalParamDescription{If set to true, the output coordinates are expressed in the integer grid (pixel indices).
 *                           Otherwise, they are given in the input coordinate system.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{true}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{The traits class must respect the concept of `HotPixelSnapRoundingTraits_2`}
 *     \cgalParamDefault{an instance of `CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>` where Kernel is deduced from the point type of the polygons, using `CGAL::Kernel_traits`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 * @warning a convex input polygon might no longer be convex after rounding.
 */
template <typename PolygonRange, typename OutputContainer,
          typename NamedParameters = parameters::Default_named_parameters>
void hot_pixel_snap_rounding_2(const PolygonRange& polygons,
                               OutputContainer& out,
                               const NamedParameters& np = parameters::default_values());

# else

template <typename InputRange, typename OutputContainer,
          typename NamedParameters = parameters::Default_named_parameters>
void hot_pixel_snap_rounding_2(const InputRange& inputs, OutputContainer& out,
                               const NamedParameters& np = parameters::default_values()) {
  using Input_iterator = typename InputRange::iterator;
  using Input_value = typename std::iterator_traits<Input_iterator>::value_type;
  using Input = std::remove_cv_t<Input_value>;

  if constexpr(internal::is_instance_of_Polygon_2<Input>) {
    SR_internal::hot_pixel_snap_rounding_2_polygon(inputs, out, np);
    return;
  }
  else {
    using SR_internal::hot_pixel_snap_rounding_2;
    using Polyline = std::remove_cv_t<typename OutputContainer::value_type>;

    using Kernel = typename Kernel_traits<Input>::Kernel;
    using DefaultTraits = Hot_pixel_snap_rounding_traits_2<Kernel>;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    const auto gt = choose_parameter<DefaultTraits>(get_parameter(np, internal_np::geom_traits));
    using Traits = std::decay_t<decltype(gt)>;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;

    auto pixel_size = choose_parameter(get_parameter(np, internal_np::pixel_size), 1.);
    bool do_isr = choose_parameter(get_parameter(np, internal_np::do_iterative_snap_rounding), true);
    bool int_output = choose_parameter(get_parameter(np, internal_np::use_grid_coordinates), true);
    bool unique_segments = choose_parameter(get_parameter(np, internal_np::output_unique_segments), false);

    if (unique_segments){
      // Output Segments while removing duplicate ones
      std::vector< std::vector<Point_2>> output_container;
      hot_pixel_snap_rounding_2(inputs.begin(), inputs.end(), output_container, gt, pixel_size, do_isr, int_output);

      // Build a set with lexicographic order of the segments
      auto comp = gt.compare_xy_2_object();
      auto less = gt.less_xy_2_object();
      auto pred = [&](const Segment_2& a, const Segment_2 &b){
        auto res = comp(a.source(), b.source());
        if (res != EQUAL)
          return res == SMALLER;
        return less(a.target(), b.target());
      };
      std::vector<Segment_2> out_segs;
      for (auto& poly: output_container)
        for (std::size_t i = 1; i < poly.size(); ++i)
          if (less(poly[i-1], poly[i]))
            out_segs.emplace_back(poly[i-1], poly[i]);
          else
            out_segs.emplace_back(poly[i], poly[i-1]);
      std::sort(out_segs.begin(), out_segs.end(), pred);
      auto last = std::unique(out_segs.begin(), out_segs.end());
      out_segs.erase(last, out_segs.end());
      for (Segment_2 &s: out_segs){
        Polyline pl;
        pl.push_back(std::move(s.source()));
        pl.push_back(std::move(s.target()));
        out.push_back(std::move(pl));
      }
    }
    else {
      // Output polylines
      hot_pixel_snap_rounding_2(inputs.begin(), inputs.end(), out, gt, pixel_size, do_isr, int_output);
    }
  }
}

#endif

} //namespace CGAL

#endif
