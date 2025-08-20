// Copyright (c) 2012
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel     <efifogel@gmail.com>
//            Shepard Liu	  <shepard0liu@gmail.com>

#ifndef CGAL_DRAW_ARRANGEMENT_2_H
#define CGAL_DRAW_ARRANGEMENT_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Random.h>
#include <CGAL/config.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/unordered_flat_map.h>
#include <CGAL/Draw_aos/type_utils.h>
#include <CGAL/Draw_aos/Arr_viewer.h>
#include "CGAL/Bbox_2.h"

namespace CGAL {

namespace draw_aos {

template <typename Arr, typename GSOptions>
class Draw_arr_tool
{
public:
  using Halfedge_const_handle = typename Arr::Halfedge_const_handle;
  using Vertex_const_handle = typename Arr::Vertex_const_handle;
  using Face_const_handle = typename Arr::Face_const_handle;
  using Ccb_halfedge_const_circulator = typename Arr::Ccb_halfedge_const_circulator;
  using Inner_ccb_const_iterator = typename Arr::Inner_ccb_const_iterator;
  using Outer_ccb_const_iterator = typename Arr::Outer_ccb_const_iterator;
  using Gt = typename Arr::Geometry_traits_2;
  using Point = typename Arr::Point_2;
  using X_monotone_curve = typename Arr::X_monotone_curve_2;

  /*! Construct
   */
  Draw_arr_tool(Arr& a_aos, CGAL::Graphics_scene& a_gs, const GSOptions& a_gso)
      : m_aos(a_aos)
      , m_gs(a_gs)
      , m_gso(a_gso) {}

  /// Add a face.
  void add_face(Face_const_handle face) {
    // std::cout << "add_face()\n";
    for(Inner_ccb_const_iterator it = face->inner_ccbs_begin(); it != face->inner_ccbs_end(); ++it) add_ccb(*it);

    if(!face->is_unbounded()) {
      for(Outer_ccb_const_iterator it = face->outer_ccbs_begin(); it != face->outer_ccbs_end(); ++it) {
        add_ccb(*it);
        draw_region(*it);
      }
    }
  }

  /// Add a Connected Component of the Boundary.
  void add_ccb(Ccb_halfedge_const_circulator circ) {
    // std::cout << "add_ccb()\n";
    auto curr = circ;
    do {
      auto new_face = curr->twin()->face();
      if(m_visited.find(new_face) != m_visited.end()) continue;
      m_visited[new_face] = true;
      add_face(new_face);
    } while(++curr != circ);
  }

  /// Draw a region.
  void draw_region(Ccb_halfedge_const_circulator circ) {
    // std::cout << "draw_region()\n";
    /* Check whether the traits has a member function called
     * approximate_2_object() and if so check whether the return type, namely
     * `Approximate_2` has an appropriate operator.
     *
     * C++20 supports concepts and `requires` expression; see, e.g.,
     * https://en.cppreference.com/w/cpp/language/constraints; thus, the first
     * condition above can be elegantly verified as follows:
     * constexpr bool has_approximate_2_object =
     *   requires(const Gt& traits) { traits.approximate_2_object(); };
     *
     * C++17 has experimental constructs called is_detected and
     * is_detected_v that can be used to achieve the same goal.
     *
     * For now we use C++14 features.
     */
    if(m_gso.colored_face(m_aos, circ->face()))
      m_gs.face_begin(m_gso.face_color(m_aos, circ->face()));
    else
      m_gs.face_begin();

    const auto* traits = this->m_aos.geometry_traits();
    auto ext = find_smallest(circ, *traits);
    auto curr = ext;

    do {
      // Skip halfedges that are "antenas":
      while(curr->face() == curr->twin()->face()) curr = curr->twin()->next();
      while(curr->face() == curr->twin()->face()) curr = curr->twin()->next();
      draw_region_impl1(*traits, curr);
      curr = curr->next();
    } while(curr != ext);

    m_gs.face_end();
  }

  /// Compile time dispatching

  ///
  template <typename T, typename A, std::enable_if_t<!has_approximate_point_v<T, A>, int> = 0>
  void draw_region_impl2(const T& /* traits */, const A& /* approximate */, Halfedge_const_handle curr) {
    draw_exact_region(curr);
  }

  ///
  template <typename T, typename A, std::enable_if_t<has_approximate_point_v<T, A>, int> = 0>
  auto draw_region_impl2(const T& /* traits */, const A& approx, Halfedge_const_handle curr) {
    draw_approximate_region(curr, approx);
  }

  /*! Draw a region, where the traits does not has approximate_2_object.
   */
  template <typename T, std::enable_if_t<!has_approximate_2_object_v<T> && !is_or_derived_from_agas_v<T>, int> = 0>
  void draw_region_impl1(const T& /* traits */, Halfedge_const_handle curr) {
    draw_exact_region(curr);
  }

  ///
  template <typename T, std::enable_if_t<has_approximate_2_object_v<T> && !is_or_derived_from_agas_v<T>, int> = 0>
  auto draw_region_impl1(const T& traits, Halfedge_const_handle curr) {
    using Approximate = typename Gt::Approximate_2;
    draw_region_impl2(traits, traits.approximate_2_object(), curr);
  }

  /*! Draw a geodesic region
   */
  template <typename T, std::enable_if_t<is_or_derived_from_agas_v<T>, int> = 0>
  void draw_region_impl1(const T& traits, Halfedge_const_handle curr) {
    //! \todo not implemented yet; for now, we just draw the boundaries using draw_curve_impl1()
    draw_curve_impl1(traits, curr->curve(), false, CGAL::IO::Color());
  }

  /*! Draw a region using approximate coordinates.
   * Call this member function only if the geometry traits is equipped with
   * the coordinate-approximation functionality of a curve.
   * This function must be inlined (e.g., a template) to enable the
   * compiled-time dispatching in the function `draw_region()`.
   */
  template <typename Approximate>
  void draw_approximate_region(Halfedge_const_handle curr, const Approximate& approx) {
    // std::cout << "draw_approximate_region()\n";
    std::vector<typename Gt::Approximate_point_2> polyline;
    double error(0.01); // TODO? (this->pixel_ratio());
    bool l2r = curr->direction() == ARR_LEFT_TO_RIGHT;
    approx(curr->curve(), error, std::back_inserter(polyline), l2r);
    if(polyline.empty()) return;
    auto it = polyline.begin();
    auto prev = it++;
    for(; it != polyline.end(); prev = it++) m_gs.add_point_in_face(*prev);
  }

  /*! Draw an exact curve.
   */
  template <typename XMonotoneCurve>
  void draw_exact_curve(const XMonotoneCurve& curve, bool colored, const CGAL::IO::Color& c) {
    const auto* traits = this->m_aos.geometry_traits();
    auto ctr_min = traits->construct_min_vertex_2_object();
    auto ctr_max = traits->construct_max_vertex_2_object();
    m_gs.add_segment(ctr_min(curve), ctr_max(curve));
    if(colored)
      m_gs.add_segment(ctr_min(curve), ctr_max(curve), c);
    else
      m_gs.add_segment(ctr_min(curve), ctr_max(curve));
  }

  /*! Draw a region in an exact manner.
   *  This fallback simply draws the curve in an exact manner (and even this is not guaranteed).
   */
  void draw_exact_region(Halfedge_const_handle curr) { draw_exact_curve(curr->curve(), false, CGAL::IO::Color()); }

  /// Add all faces.
  template <typename Traits>
  void add_faces(const Traits&) {
    for(auto it = m_aos.unbounded_faces_begin(); it != m_aos.unbounded_faces_end(); ++it) add_face(it);
  }

  /// Compile time dispatching

  /*! Draw a point using approximate coordinates.
   */
  template <typename Approximate>
  void draw_approximate_point(const Point& p, const Approximate& approx, bool colored, const CGAL::IO::Color& color) {
    if(colored)
      m_gs.add_point(approx(p), color);
    else
      m_gs.add_point(approx(p));
  }

  ///
  void draw_exact_point(const Point& p, bool colored, const CGAL::IO::Color& color) {
    if(colored)
      m_gs.add_point(p, color);
    else
      m_gs.add_point(p);
  }

  ///
  template <typename T, typename A, std::enable_if_t<!has_approximate_point_v<T, A>, int> = 0>
  void draw_point_impl2(
      const T& /* traits */, const A& /* approximate */, const Point& p, bool colored, const CGAL::IO::Color& c) {
    draw_exact_point(p, colored, c);
  }

  ///
  template <typename T, typename A, std::enable_if_t<has_approximate_point_v<T, A>, int> = 0>
  auto
  draw_point_impl2(const T& /* traits */, const A& approx, const Point& p, bool colored, const CGAL::IO::Color& c) {
    draw_approximate_point(p, approx, colored, c);
  }

  /*! Draw a point, where the traits does not has approximate_2_object.
   */
  template <typename T, std::enable_if_t<!has_approximate_2_object_v<T> && !is_or_derived_from_agas_v<T>, int> = 0>
  void draw_point_impl1(const T& /* traits */, const Point& p, bool colored, const CGAL::IO::Color& c) {
    draw_exact_point(p, colored, c);
  }

  /*! Draw a point, where the traits does have approximate_2_object.
   */
  template <typename T, std::enable_if_t<has_approximate_2_object_v<T> && !is_or_derived_from_agas_v<T>, int> = 0>
  auto draw_point_impl1(const T& traits, const Point& p, bool colored, const CGAL::IO::Color& c) {
    draw_point_impl2(traits, traits.approximate_2_object(), p, colored, c);
  }

  /*! Draw a geodesic point.
   */
  template <typename T, std::enable_if_t<is_or_derived_from_agas_v<T>, int> = 0>
  void draw_point_impl1(const T& traits, const Point& p, bool colored, const CGAL::IO::Color& color) {
    using Traits = T;
    using Ak = typename Traits::Approximate_kernel;
    using Approx_point_3 = typename Ak::Point_3;

    auto approx = traits.approximate_2_object();
    auto ap = approx(p);
    auto x = ap.dx();
    auto y = ap.dy();
    auto z = ap.dz();
    auto l = std::sqrt(x * x + y * y + z * z);
    Approx_point_3 p3(x / l, y / l, z / l);
    if(colored)
      m_gs.add_point(p3, color);
    else
      m_gs.add_point(p3);
  }

  /// Draw a point.
  void draw_point(const Point& p, bool colored, const CGAL::IO::Color& c) {
    const auto* traits = m_aos.geometry_traits();
    draw_point_impl1(*traits, p, colored, c);
  }

  ///
  template <typename Kernel, int AtanX, int AtanY>
  Halfedge_const_handle find_smallest(Ccb_halfedge_const_circulator circ,
                                      Arr_geodesic_arc_on_sphere_traits_2<Kernel, AtanX, AtanY> const&) {
    return circ;
  }

  /*! Find the halfedge incident to the lexicographically smallest vertex
   * along the CCB, such that there is no other halfedge underneath.
   */
  template <typename Traits>
  Halfedge_const_handle find_smallest(Ccb_halfedge_const_circulator circ, const Traits&) {
    // std::cout << "find_smallest()\n";
    const auto* traits = this->m_aos.geometry_traits();
    auto cmp_xy = traits->compare_xy_2_object();
    auto cmp_y = traits->compare_y_at_x_right_2_object();

    // Find the first halfedge directed from left to right
    auto curr = circ;
    do
      if(curr->direction() == CGAL::ARR_LEFT_TO_RIGHT) break;
    while(++curr != circ);
    Halfedge_const_handle ext = curr;

    // Find the halfedge incident to the lexicographically smallest vertex,
    //  such that there is no other halfedge underneath.
    do {
      // Discard edges not directed from left to right:
      if(curr->direction() != CGAL::ARR_LEFT_TO_RIGHT) continue;

      auto res = cmp_xy(curr->source()->point(), ext->source()->point());

      // Discard the edges inciden to a point strictly larger than the point
      // incident to the stored extreme halfedge:
      if(res == LARGER) continue;

      // Store the edge inciden to a point strictly smaller:
      if(res == SMALLER) {
        ext = curr;
        continue;
      }

      // The incident points are equal; compare the halfedges themselves:
      if(cmp_y(curr->curve(), ext->curve(), curr->source()->point()) == SMALLER) ext = curr;
    } while(++curr != circ);

    return ext;
  }

  /// Add all elements to be drawn.
  void add_elements() {
    // std::cout << "add_elements()\n";
    // std::cout << "ratio: " << this->pixel_ratio() << std::endl;
    m_visited.clear();

    if(m_aos.is_empty()) return;

    if(m_gso.are_faces_enabled()) add_faces(*(this->m_aos.geometry_traits()));

    // Add edges that do not separate faces.
    if(m_gso.are_edges_enabled()) {
      for(auto it = m_aos.edges_begin(); it != m_aos.edges_end(); ++it) {
        if(it->face() != it->twin()->face()) {
          if(m_gso.draw_edge(m_aos, it)) {
            if(m_gso.colored_edge(m_aos, it))
              draw_curve(it->curve(), true, m_gso.edge_color(m_aos, it));
            else
              draw_curve(it->curve(), false, CGAL::IO::Color());
          }
        }
      }
    }

    // Add all points
    if(m_gso.are_vertices_enabled()) {
      for(auto it = m_aos.vertices_begin(); it != m_aos.vertices_end(); ++it) {
        if(m_gso.colored_vertex(m_aos, it))
          draw_point(it->point(), true, m_gso.vertex_color(m_aos, it));
        else
          draw_point(it->point(), false, CGAL::IO::Color());
      }
    }

    m_visited.clear();
  }

  /*! Draw a curve using approximate coordinates.
   * Call this member function only of the geometry traits is equipped with
   * the coordinate-aproximation functionality of a curve.
   * This function must be inlined (e.g., a template) to enable the
   * compiled-time dispatching in the function `draw_curve()`.
   */
  template <typename XMonotoneCurve, typename Approximate>
  void draw_approximate_curve(const XMonotoneCurve& curve,
                              const Approximate& approx,
                              bool colored,
                              const CGAL::IO::Color& c) {
    // std::cout << "draw_approximate_curve\n";
    std::vector<typename Gt::Approximate_point_2> polyline;
    double error(0.01); // TODO? (this->pixel_ratio());
    approx(curve, error, std::back_inserter(polyline));
    if(polyline.empty()) return;
    auto it = polyline.begin();
    auto prev = it++;
    for(; it != polyline.end(); prev = it++) {
      if(colored)
        m_gs.add_segment(*prev, *it, c);
      else
        m_gs.add_segment(*prev, *it);
    }
  }

  ///
  template <typename T, typename A, std::enable_if_t<!has_approximate_point_v<T, A>, int> = 0>
  void draw_curve_impl2(const T& /* traits */,
                        const A& /* approximate */,
                        const X_monotone_curve& xcv,
                        bool colored,
                        const CGAL::IO::Color& c) {
    draw_exact_curve(xcv, colored, c);
  }

  ///
  template <typename T, typename A, std::enable_if_t<has_approximate_point_v<T, A>, int> = 0>
  auto draw_curve_impl2(
      const T& /* traits */, const A& approx, const X_monotone_curve& xcv, bool colored, const CGAL::IO::Color& c) {
    draw_approximate_curve(xcv, approx, colored, c);
  }

  /*! Draw a curve, where the traits does not has approximate_2_object.
   */
  template <typename T, std::enable_if_t<!has_approximate_2_object_v<T> && !is_or_derived_from_agas_v<T>, int> = 0>
  void draw_curve_impl1(const T& /* traits */, const X_monotone_curve& xcv, bool colored, const CGAL::IO::Color& c) {
    draw_exact_curve(xcv, colored, c);
  }

  /*! Draw a curve, where the traits does have approximate_2_object.
   */
  template <typename T, std::enable_if_t<has_approximate_2_object_v<T> && !is_or_derived_from_agas_v<T>, int> = 0>
  auto draw_curve_impl1(const T& traits, const X_monotone_curve& xcv, bool colored, const CGAL::IO::Color& c) {
    using Approximate = typename Gt::Approximate_2;
    draw_curve_impl2(traits, traits.approximate_2_object(), xcv, colored, c);
  }

  /*! Draw a geodesic curve
   */
  template <typename T, std::enable_if_t<is_or_derived_from_agas_v<T>, int> = 0>
  void draw_curve_impl1(const T& traits, const X_monotone_curve& xcv, bool colored, const CGAL::IO::Color& c) {
    // std::cout << "draw_curve (geodesic)\n";
    using Traits = T;
    using Kernel = typename Traits::Kernel;
    using Ak = typename Traits::Approximate_kernel;
    using Ap = typename Traits::Approximate_point_2;
    using Approx_point_3 = typename Ak::Point_3;

    auto approx = traits.approximate_2_object();
    std::vector<Ap> apoints;
    double error(0.01);
    approx(xcv, error, std::back_inserter(apoints));
    auto it = apoints.begin();
    auto x = it->dx();
    auto y = it->dy();
    auto z = it->dz();
    auto l = std::sqrt(x * x + y * y + z * z);
    Approx_point_3 prev(x / l, y / l, z / l);
    for(++it; it != apoints.end(); ++it) {
      auto x = it->dx();
      auto y = it->dy();
      auto z = it->dz();
      auto l = std::sqrt(x * x + y * y + z * z);
      Approx_point_3 next(x / l, y / l, z / l);
      if(colored)
        m_gs.add_segment(prev, next, c);
      else
        m_gs.add_segment(prev, next);
      prev = next;
    }
  }

  /// Draw a curve.
  template <typename XMonotoneCurve>
  void draw_curve(const XMonotoneCurve& curve, bool colored, const CGAL::IO::Color& c) {
    /* Check whether the traits has a member function called
     * approximate_2_object() and if so check whether the return type, namely
     * `Approximate_2` has an appropriate operator.
     *
     * C++20 supports concepts and `requires` expression; see, e.g.,
     * https://en.cppreference.com/w/cpp/language/constraints; thus, the first
     * condition above can be elegantly verified as follows:
     * constexpr bool has_approximate_2_object =
     *   requires(const Gt& traits) { traits.approximate_2_object(); };
     *
     * C++17 has experimental constructs called is_detected and
     * is_detected_v that can be used to achieve the same goal.
     *
     * For now we use C++14 features.
     */
#if 0
    if constexpr (std::experimental::is_detected_v<approximate_2_object_t, Gt>) {
      const auto* traits = this->m_aos.geometry_traits();
      auto approx = traits->approximate_2_object();
      draw_approximate_curve(curve, approx);
      return;
    }
    draw_exact_curve(curve);
#else
    const auto* traits = this->m_aos.geometry_traits();
    draw_curve_impl1(*traits, curve, colored, c);
#endif
  }

protected:
  Arr& m_aos;
  CGAL::Graphics_scene& m_gs;
  const GSOptions& m_gso;
  std::unordered_map<Face_const_handle, bool> m_visited;
};

template <typename Value1, typename Value2, typename Range1, typename Range2>
static auto map_from_pair_ranges(Range1 range1, Range2 range2) {
  CGAL_assertion_msg(range1.size() == range2.size(), "The two ranges must have the same size.");
  auto begin = boost::make_zip_iterator(boost::make_tuple(range1.begin(), range2.begin()));
  auto end = boost::make_zip_iterator(boost::make_tuple(range1.end(), range2.end()));
  auto tuple_to_pair = [](const auto& t) { return std::make_pair(boost::get<0>(t), boost::get<1>(t)); };
  return unordered_flat_map<Value1, Value2>(boost::make_transform_iterator(begin, tuple_to_pair),
                                            boost::make_transform_iterator(end, tuple_to_pair));
}

/*!
 * \brief tracking changes between an arrangement and its copy that will be later inserted to.
 *
 * \note tracks insertions only. If any other actions made(e.g. deletions, merging, etc), the state of the tracker
 * instance may become invalid.
 *
 * \tparam Arrangement
 */
template <typename Arrangement>
class Arr_insertion_tracker : Arr_observer<Arrangement>
{
  using Base = Arr_observer<Arrangement>;
  using Halfedge_handle = typename Arrangement::Halfedge_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_handle = typename Arrangement::Face_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Vertex_handle = typename Arrangement::Vertex_handle;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using X_monotone_curve_2 = typename Arrangement::X_monotone_curve_2;

protected:
  virtual void after_create_vertex(Vertex_handle v) override { m_vertex_map[v] = Vertex_const_handle(); }

  virtual void after_create_edge(Halfedge_handle e) override {
    m_halfedge_map[e] = Halfedge_const_handle();
    m_halfedge_map[e->twin()] = Halfedge_const_handle(); // twin is created as well
  }

  virtual void before_split_edge(Halfedge_handle e,
                                 Vertex_handle v,
                                 const X_monotone_curve_2& c1,
                                 const X_monotone_curve_2& c2) override {
    if(m_vertex_map.find(v) == m_vertex_map.end()) m_vertex_map[v] = Vertex_const_handle(); // v is newly created
  }

  virtual void after_split_edge(Halfedge_handle e1, Halfedge_handle e2) override {
    if(auto it = m_halfedge_map.find(e1); it == m_halfedge_map.end())
      m_halfedge_map[e2] = e1;
    else if(it->second == Halfedge_const_handle())
      m_halfedge_map[e2] = Halfedge_const_handle(); // e1 has no corresponding edge in the original arrangement
    else
      m_halfedge_map[e2] = it->second; // e1 is created by splitting an existing edge
  }

  virtual void after_split_face(Face_handle f1, Face_handle f2, bool) override {
    // Face cannot be created but by splitting an existing face.
    if(auto it = m_face_map.find(f1); it == m_face_map.end())
      m_face_map[f2] = f1;
    else
      m_face_map[f2] = it->second; // f1 is created by splitting an existing face
  }

public:
  Arr_insertion_tracker(Arrangement& arr)
      : Base(arr) {}

  /*!
   * \brief Query the original face of a given face.
   *
   * \param fh a valid face handle in the modified arrangement.
   * \return Face_const_handle
   */
  Face_const_handle original_face(Face_const_handle fh) const {
    auto it = m_face_map.find(fh);
    if(it == m_face_map.end()) return fh;
    return it->second; // new face from splitting an existing face
  }

  /*!
   * \brief Query the original halfedge of a given halfedge.
   *
   * \param heh a valid halfedge handle in the modified arrangement.
   * \return Halfedge_const_handle
   */
  Halfedge_const_handle original_halfedge(Halfedge_const_handle he) const {
    auto it = m_halfedge_map.find(he);
    if(it == m_halfedge_map.end()) return he;
    if(it->second == Halfedge_const_handle()) return Halfedge_const_handle(); // newly created halfedge
    return it->second;
  }

  /*!
   * \brief Query the original vertex of a given vertex.
   *
   * \param vh a valid vertex handle in the modified arrangement.
   * \return Vertex_const_handle
   */
  Vertex_const_handle original_vertex(Vertex_const_handle vh) const {
    auto it = m_vertex_map.find(vh);
    if(it == m_vertex_map.end()) return vh;
    if(it->second == Vertex_const_handle()) return Vertex_const_handle(); // newly created vertex
    return it->second;                                                    // it will never reach here.
  }

private:
  /*!
   * Maps tracking the changes between the original arrangement and modified arrangement.
   * The key is the current feature, and the value is the corresponding feature before modification.
   * If there is no entry about a feature, the corresponding feature is itself.
   * If the value is a invalid handle, it means that the feature is newly created and thus has no corresponding
   * feature in the original arrangement.
   */
  unordered_flat_map<Face_const_handle, Face_const_handle> m_face_map;
  unordered_flat_map<Halfedge_const_handle, Halfedge_const_handle> m_halfedge_map;
  unordered_flat_map<Vertex_const_handle, Vertex_const_handle> m_vertex_map;
};

void draw_unimplemented() {
  std::cerr << "Geometry traits type of arrangement is required to support approximation of Point_2 and "
               "X_monotone_curve_2. Traits on curved surfaces needs additional support for parameterization."
            << std::endl;
  exit(1);
}

template <typename Arrangement, typename GSOptions>
void draw_impl_planar(
    const Arrangement& arr, const GSOptions& gso, const char* title, Bbox_2 initial_bbox, QApplication& app) {
  Arr_viewer viewer(app.activeWindow(), arr, gso, title, initial_bbox);
  viewer.show();
  app.exec();
}

template <typename Arrangement, typename GSOptions>
void draw_impl_agas(
    const Arrangement& arr, const GSOptions& gso, const char* title, Bbox_2 initial_bbox, QApplication& app) {
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using X_monotone_curve_2 = typename Geom_traits::X_monotone_curve_2;
  using Direction_3 = typename Geom_traits::Direction_3;
  using Point_2 = typename Geom_traits::Point_2;
  using Agas_template_args = tmpl_args<Geom_traits>;

  Arrangement derived_arr(arr);
  auto vertex_map = map_from_pair_ranges<Vertex_const_handle, Vertex_const_handle>(derived_arr.vertex_handles(),
                                                                                   arr.vertex_handles());
  auto halfedge_map = map_from_pair_ranges<Halfedge_const_handle, Halfedge_const_handle>(derived_arr.halfedge_handles(),
                                                                                         arr.halfedge_handles());
  auto face_map =
      map_from_pair_ranges<Face_const_handle, Face_const_handle>(derived_arr.face_handles(), arr.face_handles());
  // setup tracker and insert the identification curve.
  Arr_insertion_tracker<Arrangement> tracker(derived_arr);
  X_monotone_curve_2 id_curve = arr.geometry_traits()->construct_x_monotone_curve_2_object()(
      Point_2(Direction_3(0, 0, -1), Point_2::MIN_BOUNDARY_LOC),
      Point_2(Direction_3(0, 0, 1), Point_2::MAX_BOUNDARY_LOC),
      Direction_3(Agas_template_args::atan_y, -Agas_template_args::atan_x, 0));
  insert(derived_arr, id_curve);

  // derived_gso proxies the call to the original gso
  GSOptions derived_gso(gso);
  derived_gso.draw_vertex = [&](const Arrangement&, const Vertex_const_handle& vh) {
    Vertex_const_handle original_vh = tracker.original_vertex(vh);
    if(original_vh == Vertex_const_handle() || vertex_map.find(original_vh) == vertex_map.end()) return false;
    return gso.draw_vertex(arr, vertex_map.at(original_vh));
  };
  derived_gso.colored_vertex = [&](const Arrangement&, const Vertex_const_handle& vh) {
    Vertex_const_handle original_vh = tracker.original_vertex(vh);
    if(original_vh == Vertex_const_handle() || vertex_map.find(original_vh) == vertex_map.end()) return false;
    return gso.colored_vertex(arr, vertex_map.at(original_vh));
  };
  derived_gso.vertex_color = [&](const Arrangement&, const Vertex_const_handle& vh) -> CGAL::IO::Color {
    Vertex_const_handle original_vh = tracker.original_vertex(vh);
    if(original_vh == Vertex_const_handle() || vertex_map.find(original_vh) == vertex_map.end())
      return CGAL::IO::Color();
    return gso.vertex_color(arr, vertex_map.at(original_vh));
  };
  derived_gso.draw_edge = [&](const Arrangement&, const Halfedge_const_handle& he) {
    Halfedge_const_handle original_he = tracker.original_halfedge(he);
    if(original_he == Halfedge_const_handle() || halfedge_map.find(original_he) == halfedge_map.end()) return false;
    return gso.draw_edge(arr, halfedge_map.at(original_he));
  };
  derived_gso.colored_edge = [&](const Arrangement&, const Halfedge_const_handle& he) {
    Halfedge_const_handle original_he = tracker.original_halfedge(he);
    if(original_he == Halfedge_const_handle() || halfedge_map.find(original_he) == halfedge_map.end()) return false;
    return gso.colored_edge(arr, halfedge_map.at(original_he));
  };
  derived_gso.edge_color = [&](const Arrangement&, const Halfedge_const_handle& he) -> CGAL::IO::Color {
    Halfedge_const_handle original_he = tracker.original_halfedge(he);
    if(original_he == Halfedge_const_handle() || halfedge_map.find(original_he) == halfedge_map.end())
      return CGAL::IO::Color();
    return gso.edge_color(arr, halfedge_map.at(original_he));
  };
  derived_gso.draw_face = [&](const Arrangement&, const Face_const_handle& fh) {
    Face_const_handle original_fh = tracker.original_face(fh);
    if(face_map.find(original_fh) == face_map.end()) return false;
    return gso.draw_face(arr, face_map.at(original_fh));
  };
  derived_gso.colored_face = [&](const Arrangement&, const Face_const_handle& fh) {
    Face_const_handle original_fh = tracker.original_face(fh);
    if(face_map.find(original_fh) == face_map.end()) return false;
    return gso.draw_face(arr, face_map.at(original_fh));
  };
  derived_gso.face_color = [&](const Arrangement&, const Face_const_handle& fh) -> CGAL::IO::Color {
    Face_const_handle original_fh = tracker.original_face(fh);
    if(face_map.find(original_fh) == face_map.end()) return CGAL::IO::Color();
    return gso.face_color(arr, face_map.at(original_fh));
  };

  Arr_viewer viewer(app.activeWindow(), derived_arr, derived_gso, title, initial_bbox);
  viewer.show();
  app.exec();
}

template <typename Arrangement, typename GSOptions, typename... Args>
void draw(const Arrangement& arr, const GSOptions& gso, Args&&... args) {
  using Geom_traits = typename Arrangement::Geometry_traits_2;

  if constexpr(!has_approximate_traits_v<Geom_traits>)
    return draw_unimplemented();
  else if constexpr(is_or_derived_from_agas_v<Geom_traits>)
    // Arrangements on curved surfaces require special handling. The identification curve must be present to make the
    // curved surface homeomorphic to a bounded plane.
    return draw_impl_agas(arr, gso, std::forward<Args>(args)...);
  else
    return draw_impl_planar(arr, gso, std::forward<Args>(args)...);
}

} // namespace draw_aos

/*!
 * \brief Draw an arrangement on surface.
 *
 * \tparam Arrangement
 * \tparam GSOptions
 * \param arr the arrangement to be drawn
 * \param gso graphics scene options
 * \param title title of the viewer window
 * \param initial_bbox parameter space bounding box to be shown intially. If empty, the approximate bounding box of the
 * arrangement is used. For arrangements induced by unbounded curves, the default initial bounding box is computed from
 * vertex coordinates.
 */
template <typename Arrangement, typename GSOptions>
void draw(const Arrangement& arr,
          const GSOptions& gso,
          const char* title = "2D Arrangement on Surface Viewer",
          Bbox_2 initial_bbox = Bbox_2(0, 0, 0, 0)) {
  Qt::init_ogl_context(4, 3);
  int argc;
  QApplication app(argc, nullptr);
  draw_aos::draw(arr, gso, title, initial_bbox, app);
}

/*!
 * \brief Draw an arrangement on surface with default graphics scene options. Faces are colored randomly.
 *
 * \tparam Arrangement
 * \param arr the arrangement to be drawn
 * \param title title of the viewer window
 * \param initial_bbox parameter space bounding box to be shown intially. If empty, the approximate bounding box of the
 * arrangement is used. For arrangements induced by unbounded curves, the default initial bounding box is computed from
 * vertex coordinates.
 */
template <typename Arrangement>
void draw(const Arrangement& arr,
          const char* title = "2D Arrangement on Surface Viewer",
          Bbox_2 initial_bbox = Bbox_2(0, 0, 0, 0)) {
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using GSOptions =
      CGAL::Graphics_scene_options<Arrangement, Vertex_const_handle, Halfedge_const_handle, Face_const_handle>;

  GSOptions gso;
  gso.enable_vertices();
  gso.draw_vertex = [](const Arrangement&, const Vertex_const_handle&) { return true; };
  gso.colored_vertex = [](const Arrangement&, const Vertex_const_handle&) { return true; };
  gso.vertex_color = [](const Arrangement&, const Vertex_const_handle& vh) -> CGAL::IO::Color {
    return CGAL::IO::Color(255, 0, 0);
  };
  gso.enable_edges();
  gso.draw_edge = [](const Arrangement&, const Halfedge_const_handle&) { return true; };
  gso.colored_edge = [](const Arrangement&, const Halfedge_const_handle&) { return true; };
  gso.edge_color = [](const Arrangement&, const Halfedge_const_handle& heh) -> CGAL::IO::Color {
    return CGAL::IO::Color(0, 0, 0);
  };
  gso.enable_faces();
  gso.draw_face = [](const Arrangement&, const Face_const_handle&) { return true; };
  gso.colored_face = [](const Arrangement&, const Face_const_handle&) { return true; };
  gso.face_color = [](const Arrangement&, const Face_const_handle& fh) -> CGAL::IO::Color {
    CGAL::Random random(std::size_t(fh.ptr()));
    return get_random_color(random);
  };

  draw(arr, gso, title, initial_bbox);
}

#define CGAL_ARR_TYPE CGAL::Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>

///
template <typename GeometryTraits_2, typename TopologyTraits, class GSOptions>
void add_to_graphics_scene(const CGAL_ARR_TYPE& aos, CGAL::Graphics_scene& graphics_scene, const GSOptions& gso) {
  draw_aos::Draw_arr_tool dar(aos, graphics_scene, gso);
  dar.add_elements();
}

///
template <typename GeometryTraits_2, typename TopologyTraits>
void add_to_graphics_scene(const CGAL_ARR_TYPE& aos, CGAL::Graphics_scene& graphics_scene) {
  CGAL::Graphics_scene_options<CGAL_ARR_TYPE, typename CGAL_ARR_TYPE::Vertex_const_handle,
                               typename CGAL_ARR_TYPE::Halfedge_const_handle, typename CGAL_ARR_TYPE::Face_const_handle>
      gso;
  // colored face?
  gso.colored_face = [](const CGAL_ARR_TYPE&, typename CGAL_ARR_TYPE::Face_const_handle) -> bool { return true; };

  // face color
  gso.face_color = [](const CGAL_ARR_TYPE&, typename CGAL_ARR_TYPE::Face_const_handle fh) -> CGAL::IO::Color {
    CGAL::Random random((unsigned int)(std::size_t)(&*fh));
    return get_random_color(random);
  };

  add_to_graphics_scene(aos, graphics_scene, gso);
}

/// Draw an arrangement on surface.
template <typename GeometryTraits_2, typename TopologyTraits, class GSOptions>
void draw_old(const CGAL_ARR_TYPE& aos,
              const GSOptions& gso,
              const char* title = "2D Arrangement on Surface Basic Viewer") {
  CGAL::Graphics_scene graphics_scene;
  add_to_graphics_scene(aos, graphics_scene, gso);
  draw_graphics_scene(graphics_scene, title);
}

/// Draw an arrangement on surface.
template <typename GeometryTraits_2, typename TopologyTraits>
void draw_old(const CGAL_ARR_TYPE& aos, const char* title = "2D Arrangement on Surface Basic Viewer") {
  CGAL::Graphics_scene graphics_scene;
  add_to_graphics_scene(aos, graphics_scene);
  draw_graphics_scene(graphics_scene, title);
}

} // namespace CGAL

#endif
