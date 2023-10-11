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
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_DRAW_ARRANGEMENT_2_H
#define CGAL_DRAW_ARRANGEMENT_2_H

#include <CGAL/config.h>

#include <unordered_map>
#include <cstdlib>
#include <random>

#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <type_traits>
#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

namespace CGAL {

struct Default_color_generator {
  /*! Obtain color
   */
  template <typename HalfedgeHandle>
  CGAL::IO::Color operator()(HalfedgeHandle /* h */) {
    static std::random_device rd;
    static std::mt19937 rng(rd());
    static std::uniform_int_distribution<int> uni(0, 255);
    return CGAL::IO::Color(uni(rng), uni(rng), uni(rng));
  }
};

// Viewer class for`< Polygon_2
template <typename ArrangementOnSurface_2,
          typename ColorGenerator = Default_color_generator>
class Aos_2_basic_viewer_qt : public Basic_viewer_qt {
  using Aos = ArrangementOnSurface_2;
  using Color_generator = ColorGenerator;
  using Base = Basic_viewer_qt;
  using Gt = typename Aos::Geometry_traits_2;
  using Point = typename Aos::Point_2;
  using X_monotone_curve = typename Aos::X_monotone_curve_2;
  using Vertex_const_handle = typename Aos::Vertex_const_handle;
  using Halfedge_const_handle = typename Aos::Halfedge_const_handle;
  using Face_const_handle = typename Aos::Face_const_handle;
  using Ccb_halfedge_const_circulator =
    typename Aos::Ccb_halfedge_const_circulator;

public:
  /// Construct the viewer.
  /// @param arr the arrangement to view
  /// @param title the title of the window
  Aos_2_basic_viewer_qt(QWidget* parent, const Aos& aos,
                        Color_generator color_generator,
                        const char* title = "2D Arrangement Basic Viewer",
                        bool draw_vertices = false) :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, draw_vertices, true, true, false, false),
    m_aos(aos),
    m_color_generator(color_generator)
  {
    // mimic the computation of Camera::pixelGLRatio()
    auto bbox = bounding_box();
    CGAL::qglviewer::Vec minv(bbox.xmin(), bbox.ymin(), 0);
    CGAL::qglviewer::Vec maxv(bbox.xmax(), bbox.ymax(), 0);
    auto diameter = (maxv - minv).norm();
    m_pixel_ratio = diameter / m_height;
  }

  /*! Intercept the resizing of the window.
   */
  virtual void resizeGL(int width, int height) {
    CGAL::QGLViewer::resizeGL(width, height);
    m_width = width;
    m_height = height;
    CGAL::qglviewer::Vec p;
    auto ratio = camera()->pixelGLRatio(p);
    if (ratio != m_pixel_ratio) {
      m_pixel_ratio = ratio;
      add_elements();
    }
  }

  /*! Compute an approximation of the bounding box of a point.
   * \param[in] p the (exact) point.
   * Call this member function only if the geometry traits is equipped with
   * the coordinate-approximation functionality of a point coordinate.
   * This function must be inlined (e.g., a template) to enable the
   * compiled-time dispatching in the function `bounding_box()`.
   */
  template <typename Point, typename Approximate>
  CGAL::Bbox_2 approximate_bbox(const Point& p, const Approximate& approx) {
    auto x = approx(p, 0);
    auto y = approx(p, 1);
    return CGAL::Bbox_2(x, y, x, y);
  }

  /*! Obtain the bounding box of a point.
   * \param[in] p the point.
   * We assume that if the coordinate-approximation functionality is not
   * supported, the point supports the member function `bbox()`.
   * This function must be inlined (e.g., a template) to enable the
   * compiled-time dispatching in the function `bounding_box()`.
   */
  template <typename Point>
  CGAL::Bbox_2 exact_bbox(const Point& p) { return p.bbox(); }

  /*! Compile time dispatching
   */
#if 0
  template <typename T>
  void bounding_box_impl2(CGAL::Bbox_2& bbox, const Point& p, T const&, long)
  { bbox += exact_bbox(p); }

  template <typename T>
  auto bounding_box_impl2(CGAL::Bbox_2& bbox, const Point& p, T const& approx,
                          int) -> decltype(approx.operator()(p), void())
  { bbox += approximate_bbox(p, approx); }

  template <typename T>
  void bounding_box_impl1(CGAL::Bbox_2& bbox, const Point& p, T const&, long)
  { bbox += exact_bbox(p); }

  template <typename T>
  auto bounding_box_impl1(CGAL::Bbox_2& bbox, const Point& p, T const& traits,
                          int) ->
    decltype(traits.approximate_2_object(), void()) {
    using Approximate = typename Gt::Approximate_2;
    bounding_box_impl2<Approximate>(bbox, p, traits.approximate_2_object(), 0);
  }
#else
  template <typename T>
  void bounding_box_impl1(CGAL::Bbox_2& bbox, const Point& p, T const& traits,
                          int)
  { bbox += approximate_bbox(p, traits.approximate_2_object()); }
#endif

  /*! Compute the bounding box.
   */
  CGAL::Bbox_2 bounding_box() {
    CGAL::Bbox_2 bbox;
    const auto* traits = this->m_aos.geometry_traits();
    // At this point we assume that the arrangement is not open, and thus the
    // bounding box is defined by the vertices.
    for (auto it = m_aos.vertices_begin(); it != m_aos.vertices_end(); ++it)
      bounding_box_impl1(bbox, it->point(), *traits, 0);
    return bbox;
  }

  /*! Add all faces.
   */
  template <typename Traits>
  void add_faces(const Traits&) {
    for (auto it = m_aos.unbounded_faces_begin();
         it != m_aos.unbounded_faces_end(); ++it)
      add_face(it);
  }

  /*! Add all faces.
   */
  template <typename Kernel_, int AtanX, int AtanY>
  void
  add_faces(Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY> const&)
  { add_face(m_aos.faces_begin()); }

  /*! Add all elements to be drawn.
   */
  void add_elements() {
    // std::cout << "add_elements()\n";
    // std::cout << "ratio: " << this->pixel_ratio() << std::endl;
    clear();
    m_visited.clear();

    if (m_aos.is_empty()) return;
    add_faces(*(this->m_aos.geometry_traits()));

    // Add edges that do not separate faces.
    for (auto it = m_aos.edges_begin(); it != m_aos.edges_end(); ++it)
      if (it->face() == it->twin()->face()) draw_curve(it->curve());

    // Add all points
    for (auto it = m_aos.vertices_begin(); it != m_aos.vertices_end(); ++it)
      draw_point(it->point());

    m_visited.clear();
  }

  /*/ Obtain the pixel ratio
   */
  double pixel_ratio() const { return m_pixel_ratio; }

protected:
  template <typename Kernel, int AtanX, int AtanY>
  Halfedge_const_handle
  find_smallest(Ccb_halfedge_const_circulator circ,
                Arr_geodesic_arc_on_sphere_traits_2<Kernel, AtanX, AtanY> const&)
  { return circ; }

  /*! Find the halfedge incident to the lexicographically smallest vertex
   * along the CCB, such that there is no other halfedge underneath.
   */
  template <typename Traits>
  Halfedge_const_handle find_smallest(Ccb_halfedge_const_circulator circ,
                                      const Traits&) {
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
      if (cmp_y(curr->curve(), ext->curve(), curr->source()->point()) ==
          SMALLER)
        ext = curr;
    } while (++curr != circ);

    return ext;
  }

  /*! Draw a region using aproximate coordinates.
   * Call this member function only if the geometry traits is equipped with
   * the coordinate-approximation functionality of a curve.
   * This function must be inlined (e.g., a template) to enable the
   * compiled-time dispatching in the function `draw_region()`.
   */
  template <typename Approximate>
  void draw_approximate_region(Halfedge_const_handle curr,
                               const Approximate& approx) {
    // std::cout << "draw_approximate_region()\n";
    std::vector<typename Gt::Approximate_point_2> polyline;
    double error(this->pixel_ratio());
    bool l2r = curr->direction() == ARR_LEFT_TO_RIGHT;
    approx(curr->curve(), error, std::back_inserter(polyline), l2r);
    if (polyline.empty()) return;
    auto it = polyline.begin();
    auto prev = it++;
    for (; it != polyline.end(); prev = it++) {
      this->add_segment(*prev, *it);
      this->add_point_in_face(*prev);
    }
  }

  /*! Draw an exact curve.
   */
  template <typename XMonotoneCurve>
  void draw_exact_curve(const XMonotoneCurve& curve) {
    const auto* traits = this->m_aos.geometry_traits();
    auto ctr_min = traits->construct_min_vertex_2_object();
    auto ctr_max = traits->construct_max_vertex_2_object();
    this->add_segment(ctr_min(curve), ctr_max(curve));
  }

  /*! Draw an exact region.
   */
  void draw_exact_region(Halfedge_const_handle curr) {
    // this->add_point_in_face(curr->source()->point());
    draw_exact_curve(curr->curve());
  }

  /*! Compile time dispatching
   */
#if 0
  template <typename T, typename I = void>
  void draw_region_impl2(Halfedge_const_handle curr, T const&, long)
  { draw_exact_region(curr); }

  template <typename T, typename I>
  auto draw_region_impl2(Halfedge_const_handle curr, T const& approx, int) ->
    decltype(approx.template operator()<I>(X_monotone_curve{}, double{}, I{},
                                           bool{}), void())
  { draw_approximate_region(curr, approx); }

  template <typename T>
  void draw_region_impl1(Halfedge_const_handle curr, T const&, long)
  { draw_exact_region(curr); }

  template <typename T>
  auto draw_region_impl1(Halfedge_const_handle curr, T const& traits, int) ->
    decltype(traits.approximate_2_object(), void()) {
    using Approximate = typename Gt::Approximate_2;
    draw_region_impl2<Approximate, int>(curr, traits.approximate_2_object(), 0);
  }
#else
  template <typename T>
  void draw_region_impl1(Halfedge_const_handle curr, T const& traits, int)
  { draw_approximate_region(curr, traits.approximate_2_object()); }
#endif

  template <typename Kernel_, int AtanX, int AtanY>
  void draw_region_impl1
  (Halfedge_const_handle curr,
   Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY> const& traits,
   int) {
    // std::cout << "draw_region_impl1()\n";
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
      this->add_segment(prev, next);
      prev = next;
      // this->add_point_in_face(*prev);
    }
  }

  /*! Draw a region.
   */
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
    auto color = m_color_generator(circ->face());
    this->face_begin(color);

    const auto* traits = this->m_aos.geometry_traits();
    auto ext = find_smallest(circ, *traits);
    auto curr = ext;

    do {
      // Skip halfedges that are "antenas":
      while (curr->face() == curr->twin()->face()) curr = curr->twin()->next();
      draw_region_impl1(curr, *traits, 0);
      curr = curr->next();
    } while (curr != ext);

    this->face_end();
  }

  /*! Draw a curve using aproximate coordinates.
   * Call this member function only of the geometry traits is equipped with
   * the coordinate-aproximation functionality of a curve.
   * This function must be inlined (e.g., a template) to enable the
   * compiled-time dispatching in the function `draw_curve()`.
   */
  template <typename XMonotoneCurve, typename Approximate>
  void draw_approximate_curve(const XMonotoneCurve& curve,
                              const Approximate& approx) {
    std::vector<typename Gt::Approximate_point_2> polyline;
    double error(this->pixel_ratio());
    approx(curve, error, std::back_inserter(polyline));
    if (polyline.empty()) return;
    auto it = polyline.begin();
    auto prev = it++;
    for (; it != polyline.end(); prev = it++) this->add_segment(*prev, *it);
  }

  /*! Compile time dispatching
   */
#if 0
  template <typename T, typename I = void>
  void draw_curve_impl2(const X_monotone_curve& xcv, T const&, long)
  { draw_exact_curve(xcv); }

  template <typename T, typename I>
  auto draw_curve_impl2(const X_monotone_curve& xcv, T const& approx, int) ->
    decltype(approx.template operator()<I>(X_monotone_curve{}, double{}, I{},
                                           bool{}), void())
  { draw_approximate_curve(xcv, approx); }

  template <typename T>
  void draw_curve_impl1(const X_monotone_curve& xcv, T const&, long)
  { draw_exact_curve(xcv); }

  template <typename T>
  auto draw_curve_impl1(const X_monotone_curve& xcv, T const& traits, int) ->
    decltype(traits.approximate_2_object(), void()) {
    using Approximate = typename Gt::Approximate_2;
    draw_curve_impl2<Approximate, int>(xcv, traits.approximate_2_object(), 0);
  }
#else
  template <typename T>
  void draw_curve_impl1(const X_monotone_curve& xcv, T const& traits, int)
  { draw_approximate_curve(xcv, traits.approximate_2_object()); }
#endif

  template <typename Kernel_, int AtanX, int AtanY>
  void draw_curve_impl1
  (const X_monotone_curve& xcv,
   Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY> const& traits,
   int) {
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
      this->add_segment(prev, next);
      prev = next;
    }
  }

  /*! Draw a curve.
   */
  template <typename XMonotoneCurve>
  void draw_curve(const XMonotoneCurve& curve) {
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
    if constexpr (std::experimental::is_detected_v<approximate_2_object_t, Gt>)
    {
      const auto* traits = this->m_aos.geometry_traits();
      auto approx = traits->approximate_2_object();
      draw_approximate_curve(curve, approx);
      return;
    }
    draw_exact_curve(curve);
#else
    const auto* traits = this->m_aos.geometry_traits();
    draw_curve_impl1(curve, *traits, 0);
#endif
  }

  /*! Compile time dispatching
   */
#if 0
  template <typename T>
  void draw_point_impl2(const Point& p, T const&, long) { add_point(p); }

  template <typename T>
  auto draw_point_impl2(const Point& p, T const& approx, int) ->
    decltype(approx.operator()(p), void())
  { add_point(approx(p)); }

  template <typename T>
  void draw_point_impl1(const Point& p, T const&, long) { add_point(p); }

  template <typename T>
  auto draw_point_impl1(const Point& p, T const& traits, int) ->
    decltype(traits.approximate_2_object(), void()) {
    using Approximate = typename Gt::Approximate_2;
    draw_point_impl2<Approximate>(p, traits.approximate_2_object(), true);
  }
#else
  template <typename T>
  void draw_point_impl1(const Point& p, T const& traits, int)
  { add_point(traits.approximate_2_object()(p)); }
#endif

  template <typename Kernel_, int AtanX, int AtanY>
  void draw_point_impl1
  (const Point& p,
   Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY> const& traits,
   int) {
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
    add_point(p3);
  }

  /*! Draw a point.
   */
  void draw_point(const Point& p) {
    const auto* traits = m_aos.geometry_traits();
    draw_point_impl1(p, *traits, 0);
  }

  /*! Add a Connected Component of the Boundary.
   */
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

  /*! Add a face.
   */
  void add_face(Face_const_handle face) {
    // std::cout << "add_face()\n";
    using Inner_ccb_const_iterator = typename Aos::Inner_ccb_const_iterator;
    using Outer_ccb_const_iterator = typename Aos::Outer_ccb_const_iterator;

    for (Inner_ccb_const_iterator it = face->inner_ccbs_begin();
         it != face->inner_ccbs_end(); ++it)
      add_ccb(*it);

    for (Outer_ccb_const_iterator it = face->outer_ccbs_begin();
         it != face->outer_ccbs_end(); ++it) {
      add_ccb(*it);
      draw_region(*it);
    }
  }

  //!
  virtual void keyPressEvent(QKeyEvent* e) {
    // Test key pressed:
    //    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    //    if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton)) { ... }

    // Call: * add_elements() if the model changed, followed by
    //       * redraw() if some viewing parameters changed that implies some
    //                  modifications of the buffers
    //                  (eg. type of normal, color/mono)
    //       * update() just to update the drawing

    // Call the base method to process others/classicals key
    Base::keyPressEvent(e);
  }

protected:
  //! The window width in pixels.
  int m_width = CGAL_BASIC_VIEWER_INIT_SIZE_X;

  //! The window height in pixels.
  int m_height = CGAL_BASIC_VIEWER_INIT_SIZE_Y;

  //! The ratio between pixel and opengl units (in world coordinate system).
  double m_pixel_ratio = 1;

  //! The arrangement to draw.
  const Aos& m_aos;

  //! The color generator.
  Color_generator m_color_generator;

  std::unordered_map<Face_const_handle, bool> m_visited;
};

//! Basic viewer of a 2D arrangement.
template <typename ArrangementOnSurface_2,
          typename ColorGenerator = Default_color_generator>
class Aos_2_viewer_qt : public Aos_2_basic_viewer_qt<ArrangementOnSurface_2,
                                                     ColorGenerator> {
public:
  using Aos = ArrangementOnSurface_2;
  using Color_generator = ColorGenerator;
  using Base = Aos_2_basic_viewer_qt<Aos, Color_generator>;
  using Point = typename Aos::Point_2;
  using X_monotone_curve = typename Aos::X_monotone_curve_2;
  using Halfedge_const_handle = typename Aos::Halfedge_const_handle;
  using Face_const_handle = typename Aos::Face_const_handle;
  using Ccb_halfedge_const_circulator =
    typename Aos::Ccb_halfedge_const_circulator;

  /// Construct the viewer.
  /// @param arr the arrangement to view
  /// @param title the title of the window
  Aos_2_viewer_qt(QWidget* parent, const Aos& aos,
                  Color_generator color_generator,
                  const char* title = "2D Arrangement on Surface Basic Viewer",
                  bool draw_vertices = false) :
    Base(parent, aos, color_generator, title, draw_vertices)
  {}
};

/*! Draw an arrangement on surface.
 */
template <typename GeometryTraits_2, typename TopologyTraits>
void draw(const Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& aos,
          const char* title = "2D Arrangement on Surface Basic Viewer",
          bool draw_vertices = false) {
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (cgal_test_suite) return;
  using Gt = GeometryTraits_2;
  using Tt = TopologyTraits;
  using Aos = CGAL::Arrangement_on_surface_2<Gt, Tt>;
  using Viewer = Aos_2_viewer_qt<Aos, Default_color_generator>;

  CGAL::Qt::init_ogl_context(4,3);

  int argc = 1;
  const char* argv[2] = {"t2_viewer", nullptr};
  QApplication app(argc, const_cast<char**>(argv));
  Default_color_generator color_generator;
  Viewer mainwindow(app.activeWindow(), aos, color_generator, title,
                    draw_vertices);
  mainwindow.add_elements();
  mainwindow.show();
  app.exec();
}

/*! Draw an arrangement.
 */
template <typename GeometryTraits_2, typename Dcel>
void draw(const Arrangement_2<GeometryTraits_2, Dcel>& arr,
          const char* title = "2D Arrangement Basic Viewer",
          bool draw_vertices = false) {
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (cgal_test_suite) return;
  using Gt = GeometryTraits_2;
  using Arr = CGAL::Arrangement_2<Gt, Dcel>;
  using Viewer = Aos_2_viewer_qt<Arr, Default_color_generator>;

  CGAL::Qt::init_ogl_context(4,3);

  int argc = 1;
  const char* argv[2] = {"t2_viewer", nullptr};
  QApplication app(argc, const_cast<char**>(argv));
  Default_color_generator color_generator;
  Viewer mainwindow(app.activeWindow(), arr, color_generator, title,
                    draw_vertices);
  mainwindow.add_elements();
  mainwindow.show();
  app.exec();
}

/*! Draw an arrangement using a given color generator.
 */
template <typename GeometryTraits_2, typename Dcel,
          typename ColorGenerator>
void draw(const Arrangement_2<GeometryTraits_2, Dcel>& arr,
          ColorGenerator color_generator,
          const char* title = "2D Arrangement Basic Viewer",
          bool draw_vertices = false) {
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (cgal_test_suite) return;

  using Color_generator = ColorGenerator;
  using Gt = GeometryTraits_2;
  using Arr = CGAL::Arrangement_2<Gt, Dcel>;
  using Viewer = Aos_2_viewer_qt<Arr, Color_generator>;

  CGAL::Qt::init_ogl_context(4,3);

  int argc = 1;
  const char* argv[2] = {"t2_viewer", nullptr};
  QApplication app(argc, const_cast<char**>(argv));
  Viewer mainwindow(app.activeWindow(), arr, color_generator, title,
                    draw_vertices);
  mainwindow.add_elements();
  mainwindow.show();
  app.exec();
}

}

#endif
#endif
