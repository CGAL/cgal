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

#include <unordered_map>
#include <cstdlib>
#include <random>

#include <experimental/type_traits>
#include <boost/hana.hpp>

#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Arrangement_2.h>

namespace CGAL {

// Viewer class for`< Polygon_2
template <typename GeometryTraits_2, typename Dcel>
class Arr_2_basic_viewer_qt : public Basic_viewer_qt {
  typedef GeometryTraits_2                      Gt;
  typedef CGAL::Arrangement_2<Gt, Dcel>         Arr;
  typedef Basic_viewer_qt                       Base;
  typedef typename Arr::Point_2                 Point;
  typedef typename Arr::X_monotone_curve_2      X_monotone_curve;
  typedef typename Arr::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arr::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arr::Face_const_handle       Face_const_handle;
  typedef typename Arr::Ccb_halfedge_const_circulator
                                                Ccb_halfedge_const_circulator;

  template <typename T>
  using approximate_2_object_t =
    decltype(std::declval<T&>().approximate_2_object());

public:
  /// Construct the viewer.
  /// @param arr the arrangement to view
  /// @param title the title of the window
  Arr_2_basic_viewer_qt(QWidget* parent, const Arr& arr,
                        const char* title = "2D Arrangement Basic Viewer") :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false),
    m_arr(arr),
    m_uni(0, 255)
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

  //! Compute the bounding box
  CGAL::Bbox_2 bounding_box() {
    namespace bh = boost::hana;
    auto has_approximate_2_object =
      bh::is_valid([](auto&& x) -> decltype(x.approximate_2_object()){});

    // At this point we assume that the arrangement is not open, and thus the
    // bounding box is defined by the vertices.
    //! The bounding box
    CGAL::Bbox_2 bbox;
    for (auto it = m_arr.vertices_begin(); it != m_arr.vertices_end(); ++it) {
      bh::if_(has_approximate_2_object(Gt{}),
              [&](auto& x) {
                const auto* traits = this->m_arr.geometry_traits();
                auto approx = traits->approximate_2_object();
                auto has_operator =
                  bh::is_valid([](auto&& x) ->
                               decltype(x.operator()(Point{}, int{})){});
                bh::if_(has_operator(approx),
                        [&](auto& x)
                          { bbox += x.approximate_bbox(it->point(), approx); },
                        [&](auto& x) { bbox += x.exact_bbox(it->point()); }
                        )(*this);
              },
              [&](auto& x) { bbox += x.exact_bbox(it->point()); }
              )(*this);
    }
    return bbox;
  }

  //!
  void add_elements() {
    // std::cout << "ratio: " << this->pixel_ratio() << std::endl;
    clear();
    m_visited.clear();

    std::random_device rd;
    m_rng.seed(rd());

    if (m_arr.is_empty()) return;
    for (auto it = m_arr.unbounded_faces_begin();
         it != m_arr.unbounded_faces_end(); ++it)
      add_face(it);

    // Add edges that do not separe faces.
    for (auto it = m_arr.edges_begin(); it != m_arr.edges_end(); ++it)
      if (it->face() == it->twin()->face()) draw_curve(it->curve());

    // Add all points
    for (auto it = m_arr.vertices_begin(); it != m_arr.vertices_end(); ++it)
      draw_point(it->point());

    m_visited.clear();
  }

  /*/ Obtain the pixel ratio
   */
  double pixel_ratio() const { return m_pixel_ratio; }

protected:
  /*! Find the halfedge incident to the lexicographically smallest vertex
   * along the CCB, such that there is no other halfedge underneath.
   */
  Halfedge_const_handle find_smallest(Ccb_halfedge_const_circulator circ) {
    const auto* traits = this->m_arr.geometry_traits();
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

  //! Draw an exact curve.
  template <typename XMonotoneCurve>
  void draw_exact_curve(const XMonotoneCurve& curve) {
    const auto* traits = this->m_arr.geometry_traits();
    auto ctr_min = traits->construct_min_vertex_2_object();
    auto ctr_max = traits->construct_max_vertex_2_object();
    this->add_segment(ctr_min(curve), ctr_max(curve));
  }

  //! Draw an exact region.
  void draw_exact_region(Halfedge_const_handle curr) {
    this->add_point_in_face(curr->source()->point());
    draw_exact_curve(curr->curve());
  }

  //! Utility struct
  template <typename T>
  struct can_call_operator_curve {
    template <typename F>
    constexpr auto operator()(F&& f) ->
      decltype(f.template operator()<T>(X_monotone_curve{}, double{}, T{}, bool{})) {}
  };

  //! Draw a region.
  void draw_region(Ccb_halfedge_const_circulator circ) {
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
     * For now we use C++14 features and boost.hana instead.
     */
    namespace bh = boost::hana;
    auto has_approximate_2_object =
      bh::is_valid([](auto&& x) -> decltype(x.approximate_2_object()){});

    CGAL::IO::Color color(m_uni(m_rng), m_uni(m_rng), m_uni(m_rng));
    this->face_begin(color);

    const auto* traits = this->m_arr.geometry_traits();
    auto ext = find_smallest(circ);
    auto curr = ext;
    do {
      // Skip halfedges that are "antenas":
      while (curr->face() == curr->twin()->face()) curr = curr->twin()->next();

      bh::if_(has_approximate_2_object(Gt{}),
              [&](auto& x) {
                auto approx = traits->approximate_2_object();
                auto has_operator =
                  bh::is_valid(can_call_operator_curve<int>{});
                bh::if_(has_operator(approx),
                        [&](auto& x) { x.draw_approximate_region(curr, approx); },
                        [&](auto& x) { x.draw_exact_region(curr); }
                        )(*this);
              },
              [&](auto& x) { x.draw_exact_region(curr); }
              )(*this);
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

  //! Draw a curve.
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
     * For now we use C++14 features and boost.hana instead.
     */
#if 0
    if constexpr (std::experimental::is_detected_v<approximate_2_object_t, Gt>)
    {
      const auto* traits = this->m_arr.geometry_traits();
      auto approx = traits->approximate_2_object();
      draw_approximate_curve(curve, approx);
      return;
    }
    draw_exact_curve(curve);
#else
    namespace bh = boost::hana;
    auto has_approximate_2_object =
      bh::is_valid([](auto&& x) -> decltype(x.approximate_2_object()){});
    const auto* traits = this->m_arr.geometry_traits();
    bh::if_(has_approximate_2_object(Gt{}),
            [&](auto& x) {
              auto approx = traits->approximate_2_object();
              auto has_operator = bh::is_valid(can_call_operator_curve<int>{});
              bh::if_(has_operator(approx),
                      [&](auto& x) { x.draw_approximate_curve(curve, approx); },
                      [&](auto& x) { x.draw_exact_curve(curve); }
                      )(*this);
            },
            [&](auto& x) { x.draw_exact_curve(curve); }
            )(*this);
#endif
  }

  /*! Add an approximation of a point.
   * \param[in] p the (exact) point.
   * Call this member function only if the geometry traits is equipped with
   * the coordinate-approximation functionality of a point.
   * This function must be inlined (e.g., a template) to enable the
   * compiled-time dispatching in the function `draw_point()`.
   */
  template <typename Point, typename Approximate>
  void draw_approximate_point(const Point& p, const Approximate& approx)
  { this->add_point(approx(p)); }

  //!
  void draw_point(const Point& p) {
    namespace bh = boost::hana;
    auto has_approximate_2_object =
      bh::is_valid([](auto&& x) -> decltype(x.approximate_2_object()){});
    bh::if_(has_approximate_2_object(Gt{}),
            [&](auto& x) {
              const auto* traits = x.m_arr.geometry_traits();
              auto approx = traits->approximate_2_object();
              auto has_operator =
                bh::is_valid([](auto&& x) -> decltype(x.operator()(Point{})){});
              bh::if_(has_operator(approx),
                      [&](auto& x) { x.draw_approximate_point(p, approx); },
                      [&](auto& x) { x.add_point(p); }
                      )(*this);
            },
            [&](auto& x) { x.add_point(p); }
            )(*this);
  }

  //!
  void add_ccb(Ccb_halfedge_const_circulator circ) {
    auto curr = circ;
    do {
      auto new_face = curr->twin()->face();
      if (m_visited.find(new_face) != m_visited.end()) continue;
      m_visited[new_face] = true;
      add_face(new_face);
    } while (++curr != circ);
  }

  //!
  void add_face(Face_const_handle face) {
    for (auto it = face->inner_ccbs_begin(); it != face->inner_ccbs_end(); ++it)
      add_ccb(*it);

    for (auto it = face->outer_ccbs_begin(); it != face->outer_ccbs_end(); ++it)
    {
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
  const Arr& m_arr;

  std::unordered_map<Face_const_handle, bool> m_visited;

  std::mt19937 m_rng;
  std::uniform_int_distribution<size_t> m_uni;  // guaranteed unbiased
};

//! Basic viewer of a 2D arrangement.
template <typename GeometryTraits_2, typename Dcel>
class Arr_2_viewer_qt : public Arr_2_basic_viewer_qt<GeometryTraits_2, Dcel> {
public:
  typedef GeometryTraits_2                      Gt;
  typedef CGAL::Arrangement_2<Gt, Dcel>         Arr;
  typedef Arr_2_basic_viewer_qt<Gt, Dcel>       Base;
  typedef typename Arr::Point_2                 Point;
  typedef typename Arr::X_monotone_curve_2      X_monotone_curve;
  typedef typename Arr::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arr::Face_const_handle       Face_const_handle;
  typedef typename Arr::Ccb_halfedge_const_circulator
                                                Ccb_halfedge_const_circulator;

  /// Construct the viewer.
  /// @param arr the arrangement to view
  /// @param title the title of the window
  Arr_2_viewer_qt(QWidget* parent, const Arr& arr,
                  const char* title = "2D Arrangement Basic Viewer") :
    Base(parent, arr, title)
  {}
};

//!
template <typename GeometryTraits_2, typename Dcel>
void draw(const Arrangement_2<GeometryTraits_2, Dcel>& arr,
          const char* title = "2D Arrangement Basic Viewer") {
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite) {
    typedef GeometryTraits_2            Gt;

    CGAL::Qt::init_ogl_context(4,3);

    int argc = 1;
    const char* argv[2] = {"t2_viewer", nullptr};
    QApplication app(argc, const_cast<char**>(argv));
    Arr_2_viewer_qt<Gt, Dcel> mainwindow(app.activeWindow(), arr, title);
    mainwindow.add_elements();
    mainwindow.show();
    app.exec();
  }
}

}

#endif
#endif
