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
// Author(s): Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_DRAW_ARRANGEMENT_2_H
#define CGAL_DRAW_ARRANGEMENT_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/config.h>

#include <unordered_map>
#include <cstdlib>
#include <random>
#include <type_traits>

#include <CGAL/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Random.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

namespace CGAL {

namespace draw_function_for_arrangement_2 {

// ============================
// Detection idiom using void_t
// ============================

// Primary template: detection fails by default
template <typename, typename = std::void_t<>>
struct has_approximate_2_object : std::false_type {};

// Specialization: detection succeeds if decltype(T::approximate_2_object()) is valid
template <typename T>
struct has_approximate_2_object<T, std::void_t<decltype(std::declval<T>().approximate_2_object())>> :
    std::true_type
{};

// Convenience variable template (C++17)
template <typename T>
inline constexpr bool has_approximate_2_object_v = has_approximate_2_object<T>::value;

///
template<typename Arr, typename GSOptions>
class Draw_arr_tool {
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
  Draw_arr_tool(Arr& a_aos, CGAL::Graphics_scene& a_gs, const GSOptions& a_gso) :
    m_aos(a_aos), m_gs(a_gs),
    m_gso(a_gso)
  {}

  /// Add a face.
  void add_face(Face_const_handle face) {
    // std::cout << "add_face()\n";
    for (Inner_ccb_const_iterator it = face->inner_ccbs_begin();
         it != face->inner_ccbs_end(); ++it)
     add_ccb(*it);

    for (Outer_ccb_const_iterator it = face->outer_ccbs_begin();
         it != face->outer_ccbs_end(); ++it) {
      add_ccb(*it);
      draw_region(*it);
    }
  }

  /// Add a Connected Component of the Boundary.
  void add_ccb(Ccb_halfedge_const_circulator circ) {
    // std::cout << "add_ccb()\n";
    auto curr = circ;
    do {
      auto new_face = curr->twin()->face();
      if (m_visited.find(new_face) != m_visited.end()) continue;
      m_visited[new_face] = true;
      add_face(new_face);
    } while (++curr != circ);
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
    if (m_gso.colored_face(m_aos, circ->face()))
      m_gs.face_begin(m_gso.face_color(m_aos, circ->face()));
    else m_gs.face_begin();

    const auto* traits = this->m_aos.geometry_traits();
    auto ext = find_smallest(circ, *traits);
    auto curr = ext;

    do {
      // Skip halfedges that are "antenas":
      while (curr->face() == curr->twin()->face()) curr = curr->twin()->next();
      draw_region_impl1(*traits, curr);
      curr = curr->next();
    } while (curr != ext);

    m_gs.face_end();
  }

  /// Compile time dispatching
#if 0
  template <typename T, typename I = void>
  void draw_region_impl2(Halfedge_const_handle curr, T const&, long)
  { draw_exact_region(curr); }

  template <typename T, typename I>
  auto draw_region_impl2(Halfedge_const_handle curr, T const& approx, int) ->
    decltype(approx.template operator()<I>(X_monotone_curve{}, double{}, I{},
                                           bool{}), void())
  { draw_approximate_region(curr, approx); }
#endif

  /*! Draw a region, where the traits does not has approximate_2_object.
   */
  template <typename T, std::enable_if_t<! has_approximate_2_object_v<T>, int> = 0>
  void draw_region_impl1(T const& /* traits */, Halfedge_const_handle curr)
  { draw_exact_region(curr); }

  ///
  template <typename T, std::enable_if_t<has_approximate_2_object_v<T>, int> = 0>
  auto draw_region_impl1(T const& traits, Halfedge_const_handle curr) {
    using Approximate = typename Gt::Approximate_2;
    // draw_region_impl2<Approximate, int>(xcv, traits.approximate_2_object(), 0);
    draw_approximate_region(curr, traits.approximate_2_object());
  }

  /*! Draw a geodesic region
   */
  template <typename Kernel_, int AtanX, int AtanY>
  void draw_region_impl1
  (Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY> const& traits,
   Halfedge_const_handle curr)
  {
    if (! m_gso.draw_edge(m_aos, curr)) return;

    auto approx = traits.approximate_2_object();
    using Kernel = Kernel_;
    using Traits = Arr_geodesic_arc_on_sphere_traits_2<Kernel, AtanX, AtanY>;
    using Ak = typename Traits::Approximate_kernel;
    using Ap = typename Traits::Approximate_point_2;
    using Approx_point_3 = typename Ak::Point_3;

    std::vector<Ap> polyline;
    double error(0.01);
    bool l2r = curr->direction() == ARR_LEFT_TO_RIGHT;
    approx(curr->curve(), error, std::back_inserter(polyline), l2r);
    if (polyline.empty()) return;
    auto it = polyline.begin();
    auto x = it->dx();
    auto y = it->dy();
    auto z = it->dz();
    auto l = std::sqrt(x*x + y*y + z*z);
    Approx_point_3 prev(x/l, y/l, z/l);
    for (++it; it != polyline.end(); ++it) {
      auto x = it->dx();
      auto y = it->dy();
      auto z = it->dz();
      auto l = std::sqrt(x*x + y*y + z*z);
      Approx_point_3 next(x/l, y/l, z/l);

      if (m_gso.colored_edge(m_aos, curr))
        m_gs.add_segment(prev, next, m_gso.edge_color(m_aos, curr));
      else m_gs.add_segment(prev, next);

      prev = next;
      // m_gs.add_point_in_face(*prev);
    }
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
    if (polyline.empty()) return;
    auto it = polyline.begin();
    auto prev = it++;
    for (; it != polyline.end(); prev = it++) m_gs.add_point_in_face(*prev);
  }

  /*! Draw an exact curve.
   */
  template <typename XMonotoneCurve>
  void draw_exact_curve(const XMonotoneCurve& curve,
                        bool colored, const CGAL::IO::Color& c) {
    const auto* traits = this->m_aos.geometry_traits();
    auto ctr_min = traits->construct_min_vertex_2_object();
    auto ctr_max = traits->construct_max_vertex_2_object();
    m_gs.add_segment(ctr_min(curve), ctr_max(curve));
    if (colored) m_gs.add_segment(ctr_min(curve), ctr_max(curve), c);
    else m_gs.add_segment(ctr_min(curve), ctr_max(curve));
  }

  /*! Draw a region in an exact manner.
   *  This fallback simply draws the curve in an exact manner (and even this is not guaranteed).
   */
  void draw_exact_region(Halfedge_const_handle curr)
  { draw_exact_curve(curr->curve(), false, CGAL::IO::Color()); }

  /// Add all faces.
  template <typename Traits>
  void add_faces(const Traits&) {
    for (auto it = m_aos.unbounded_faces_begin(); it != m_aos.unbounded_faces_end(); ++it)
      add_face(it);
  }

  /// Add all faces.
  template <typename Kernel_, int AtanX, int AtanY>
  void add_faces(Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY> const&)
  { add_face(m_aos.faces_begin()); }

  /// Compile time dispatching
#if 0
  template <typename T>
  void draw_point_impl2(const Point& p, T const&, long) { m_gs.add_point(p); }

  template <typename T>
  auto draw_point_impl2(const Point& p, T const& approx, int) ->
    decltype(approx.operator()(p), void())
  { m_gs.add_point(approx(p)); }

  template <typename T>
  void draw_point_impl1(T const&, const Point& p, long) { m_gs.add_point(p); }
#endif

  /*! Draw a point using approximate coordinates.
   */
  template <typename Approximate>
  void draw_approximate_point(const Point& p, const Approximate& approx,
                              bool colored, const CGAL::IO::Color& color) {
    if (colored) m_gs.add_point(approx(p), color);
    else m_gs.add_point(approx(p));
  }

  ///
  void draw_exact_point(const Point& p, bool colored, const CGAL::IO::Color& color) {
    if (colored) m_gs.add_point(p, color);
    else m_gs.add_point(p);
  }

  /*! Draw a point, where the traits does not has approximate_2_object.
   */
  template <typename T, std::enable_if_t<! has_approximate_2_object_v<T>, int> = 0>
  void draw_point_impl1(T const& /* traits */, const Point& p,
                        bool colored, const CGAL::IO::Color& c)
  { draw_exact_point(p, colored, c); }

  /*! Draw a point, where the traits does have approximate_2_object.
   */
  template <typename T, std::enable_if_t<has_approximate_2_object_v<T>, int> = 0>
  auto draw_point_impl1(T const& traits, const Point& p,
                        bool colored, const CGAL::IO::Color& c) {
    using Approximate = typename Gt::Approximate_2;
    // draw_curve_impl2<Approximate, int>(p, traits.approximate_2_object(), 0);
    draw_approximate_point(p, traits.approximate_2_object(), colored, c);
  }

  /*! Draw a geodesic point.
   */
  template <typename Kernel_, int AtanX, int AtanY>
  void draw_point_impl1
  (Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY> const& traits,
   const Point& p,
   bool colored,
   const CGAL::IO::Color& color)
  {
    auto approx = traits.approximate_2_object();
    using Traits = Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY>;
    using Ak = typename Traits::Approximate_kernel;
    using Approx_point_3 = typename Ak::Point_3;
    auto ap = approx(p);
    auto x = ap.dx();
    auto y = ap.dy();
    auto z = ap.dz();
    auto l = std::sqrt(x*x + y*y + z*z);
    Approx_point_3 p3(x/l, y/l, z/l);
    if (colored) m_gs.add_point(p3, color);
    else m_gs.add_point(p3);
  }

  /// Draw a point.
  void draw_point(const Point& p, bool colored, const CGAL::IO::Color& c) {
    const auto* traits = m_aos.geometry_traits();
    draw_point_impl1(*traits, p, colored, c);
  }

  ///
  template <typename Kernel, int AtanX, int AtanY>
  Halfedge_const_handle
  find_smallest(Ccb_halfedge_const_circulator circ,
                Arr_geodesic_arc_on_sphere_traits_2<Kernel, AtanX, AtanY> const&)
  { return circ; }

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
    do if (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT) break;
    while (++curr != circ);
    Halfedge_const_handle ext = curr;

    // Find the halfedge incident to the lexicographically smallest vertex,
    //  such that there is no other halfedge underneath.
    do {
      // Discard edges not directed from left to right:
      if (curr->direction() != CGAL::ARR_LEFT_TO_RIGHT) continue;

      auto res = cmp_xy(curr->source()->point(), ext->source()->point());

      // Discard the edges inciden to a point strictly larger than the point
      // incident to the stored extreme halfedge:
      if (res == LARGER) continue;

      // Store the edge inciden to a point strictly smaller:
      if (res == SMALLER) {
        ext = curr;
        continue;
      }

      // The incident points are equal; compare the halfedges themselves:
      if (cmp_y(curr->curve(), ext->curve(), curr->source()->point()) == SMALLER)
        ext = curr;
    } while (++curr != circ);

    return ext;
  }

  /// Add all elements to be drawn.
  void add_elements() {
    // std::cout << "add_elements()\n";
    // std::cout << "ratio: " << this->pixel_ratio() << std::endl;
    m_visited.clear();

    if (m_aos.is_empty()) return;

    if (m_gso.are_faces_enabled()) add_faces(*(this->m_aos.geometry_traits()));

    // Add edges that do not separate faces.
    if (m_gso.are_edges_enabled()) {
      std::cout << "Edges not enabled\n";
      for (auto it = m_aos.edges_begin(); it != m_aos.edges_end(); ++it) {
        if (it->face() == it->twin()->face()) {
          if (m_gso.draw_edge(m_aos, it)) {
            if (m_gso.colored_edge(m_aos, it))
              draw_curve(it->curve(), true, m_gso.edge_color(m_aos, it));
            else draw_curve(it->curve(), false, CGAL::IO::Color());
          }
        }
      }
    }

    // Add all points
    if (m_gso.are_vertices_enabled()) {
      for (auto it = m_aos.vertices_begin(); it != m_aos.vertices_end(); ++it) {
        if (m_gso.colored_vertex(m_aos, it))
          draw_point(it->point(), true, m_gso.vertex_color(m_aos, it));
        else draw_point(it->point(), false, CGAL::IO::Color());
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
                              bool colored, const CGAL::IO::Color& c) {
    // std::cout << "draw_approximate_curve\n";
    std::vector<typename Gt::Approximate_point_2> polyline;
    double error(0.01); // TODO? (this->pixel_ratio());
    approx(curve, error, std::back_inserter(polyline));
    if (polyline.empty()) return;
    auto it = polyline.begin();
    auto prev = it++;
    for (; it != polyline.end(); prev = it++) {
      if (colored) m_gs.add_segment(*prev, *it, c);
      else m_gs.add_segment(*prev, *it);
    }
  }

  /*! Compile time dispatching
   */
  // template <typename T, typename I = void>
  // void draw_curve_impl2(const X_monotone_curve& xcv, T const&, long)
  // { draw_exact_curve(xcv); }

  // template <typename T, typename I>
  // auto draw_curve_impl2(const X_monotone_curve& xcv, T const& approx, int) ->
  //   decltype(approx.template operator()<I>(X_monotone_curve{}, double{}, I{},
  //                                            bool{}), void())
  // { draw_approximate_curve(xcv, approx); }

  /*! Draw a curve, where the traits does not has approximate_2_object.
   */
  template <typename T, std::enable_if_t<! has_approximate_2_object_v<T>, int> = 0>
  void draw_curve_impl1(T const& /* traits */, const X_monotone_curve& xcv,
                        bool colored, const CGAL::IO::Color& c)
  { draw_exact_curve(xcv, colored, c); }

  /*! Draw a curve, where the traits does have approximate_2_object.
   */
  template <typename T, std::enable_if_t<has_approximate_2_object_v<T>, int> = 0>
  auto draw_curve_impl1(T const& traits, const X_monotone_curve& xcv,
                        bool colored, const CGAL::IO::Color& c) {
    using Approximate = typename Gt::Approximate_2;
    // draw_curve_impl2<Approximate, int>(xcv, traits.approximate_2_object(), 0);
    draw_approximate_curve(xcv, traits.approximate_2_object(), colored, c);
  }

  /*! Draw a geodesic curve
   */
  template <typename Kernel_, int AtanX, int AtanY>
  void draw_curve_impl1
  (Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY> const& traits,
   const X_monotone_curve& xcv,
   bool colored, const CGAL::IO::Color& c)
  {
    // std::cout << "draw_curve (geodesic)\n";
    auto approx = traits.approximate_2_object();
    using Kernel = Kernel_;
    using Traits = Arr_geodesic_arc_on_sphere_traits_2<Kernel, AtanX, AtanY>;
    using Ak = typename Traits::Approximate_kernel;
    using Ap = typename Traits::Approximate_point_2;
    using Approx_point_3 = typename Ak::Point_3;
    std::vector<Ap> apoints;
    double error(0.01);
    approx(xcv, error, std::back_inserter(apoints));
    auto it = apoints.begin();
    auto x = it->dx();
    auto y = it->dy();
    auto z = it->dz();
    auto l = std::sqrt(x*x + y*y + z*z);
    Approx_point_3 prev(x/l, y/l, z/l);
    for (++it; it != apoints.end(); ++it) {
      auto x = it->dx();
      auto y = it->dy();
      auto z = it->dz();
      auto l = std::sqrt(x*x + y*y + z*z);
      Approx_point_3 next(x/l, y/l, z/l);
      if (colored) m_gs.add_segment(prev, next, c);
      else m_gs.add_segment(prev, next);
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

} // namespace draw_function_for_arrangement_2

#define CGAL_ARR_TYPE CGAL::Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>

///
template <typename GeometryTraits_2, typename TopologyTraits, class GSOptions>
void add_to_graphics_scene(const CGAL_ARR_TYPE& aos,
                           CGAL::Graphics_scene& graphics_scene,
                           const GSOptions& gso) {
  draw_function_for_arrangement_2::Draw_arr_tool dar(aos, graphics_scene, gso);
  dar.add_elements();
}

///
template <typename GeometryTraits_2, typename TopologyTraits>
void add_to_graphics_scene(const CGAL_ARR_TYPE& aos,
                           CGAL::Graphics_scene& graphics_scene) {
  CGAL::Graphics_scene_options<CGAL_ARR_TYPE,
                               typename CGAL_ARR_TYPE::Vertex_const_handle,
                               typename CGAL_ARR_TYPE::Halfedge_const_handle,
                               typename CGAL_ARR_TYPE::Face_const_handle>
    gso;
  gso.colored_face=[](const CGAL_ARR_TYPE&,
                      typename CGAL_ARR_TYPE::Face_const_handle) -> bool
  { return true; };

  gso.face_color=[](const CGAL_ARR_TYPE&,
                    typename CGAL_ARR_TYPE::Face_const_handle fh) -> CGAL::IO::Color
  {
    CGAL::Random random((unsigned int)(std::size_t)(&*fh));
    return get_random_color(random);
  };

  add_to_graphics_scene(aos, graphics_scene, gso);
}

/// Draw an arrangement on surface.
template <typename GeometryTraits_2, typename TopologyTraits, class GSOptions>
void draw(const CGAL_ARR_TYPE& aos, const GSOptions& gso,
          const char* title = "2D Arrangement on Surface Basic Viewer") {
  CGAL::Graphics_scene graphics_scene;
  add_to_graphics_scene(aos, graphics_scene, gso);
  draw_graphics_scene(graphics_scene, title);

}

/// Draw an arrangement on surface.
template <typename GeometryTraits_2, typename TopologyTraits>
void draw(const CGAL_ARR_TYPE& aos,
          const char* title = "2D Arrangement on Surface Basic Viewer") {
  CGAL::Graphics_scene graphics_scene;
  add_to_graphics_scene(aos, graphics_scene);
  draw_graphics_scene(graphics_scene, title);
}

#undef CGAL_ARR_TYPE

} // namespace CGAL

#endif
