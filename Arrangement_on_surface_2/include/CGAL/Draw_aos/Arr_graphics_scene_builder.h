// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Efi Fogel    <efifogel@gmail.com>
//             Shepard Liu  <shepard0liu@gmail.com>

#ifndef CGAL_DRAW_AOS_ARR_GRAPHICS_SCENE_BUILDER_H
#define CGAL_DRAW_AOS_ARR_GRAPHICS_SCENE_BUILDER_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/iterator/function_output_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Random.h>
#include <CGAL/assertions.h>
#include <CGAL/number_type_config.h>
#include <CGAL/unordered_flat_map.h>
#include <CGAL/Draw_aos/Arr_bounded_renderer.h>
#include <CGAL/Draw_aos/Arr_coordinate_converter.h>
#include <CGAL/Draw_aos/Arr_face_point_generator.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/type_utils.h>

/* This file contains the machinery that converts an arrangement into the
 * primitives of a `Graphics_scene`. It does not depend on Qt; it is used both
 * by the interactive viewer (see `Arr_viewer.h`), which rebuilds the scene
 * every time the camera changes, and by `add_to_graphics_scene()`, which builds
 * the scene once.
 */

namespace CGAL {
namespace draw_aos {

/*! \brief The parameter space of an arrangement.
 *
 * Provides the bounding box of the arrangement in the parameter space, the
 * approximation error that corresponds to a given resolution, and the mapping
 * from the parameter space to the (3D) space of the graphics scene.
 *
 * \tparam Arrangement
 */
template <typename Arrangement, typename = void>
class Arr_parameter_space;

// Specialization for planar arrangements
template <typename Arrangement>
class Arr_parameter_space<Arrangement,
                          std::enable_if_t<! is_or_derived_from_curved_surf_traits_v
                                           <typename Arrangement::Geometry_traits_2>>> {
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Point = typename Approx_traits::Point;
  using Local_point = Graphics_scene::Local_point;

public:
  explicit Arr_parameter_space(const Geom_traits&) {}

  /*! computes a subpixel-level approximation error based on the bounding box
   * and the number of pixels that span its width.
   */
  double approximation_error(const Bbox_2& bbox, int resolution) const
  { return bbox.x_span() / resolution; }

  /*! computes a parameter space bounding box that contains the arrangement
   * with some margin.
   *
   * \note For an arrangement induced by unbounded curves the bounding box fits
   * the bounded parts of the curves only. If the arrangement has no bounded
   * parts at all (e.g., it is induced by a single line), the bounding box is
   * the square [-1,1]x[-1,1].
   */
  Bbox_2 bbox(const Arrangement& arr) const {
    const auto& traits = *arr.geometry_traits();
    auto approx = traits.approximate_2_object();

    // Computes a rough bounding box from the vertices.
    Bbox_2 bbox;
    for (const auto& vh : arr.vertex_handles()) bbox += approx(vh->point()).bbox();
    if (is_empty(bbox)) return Bbox_2(-1, -1, 1, 1);

    // Computes a more precise bounding box from the edges. Avoid a zero error,
    // which arises when all vertices lie on a vertical (or horizontal) line.
    if constexpr (has_approximate_xcv_v<Geom_traits, typename Geom_traits::Approximate_2>) {
      double span = (std::max)(bbox.x_span(), bbox.y_span());
      double approx_error = (span > 0) ? span / 100 : 0.01;
      auto add = boost::make_function_output_iterator([&bbox](const Approx_point& pt) { bbox += pt.bbox(); });
      for (const auto& he : arr.edge_handles()) approx(he->curve(), approx_error, add);
    }

    // Places margin around the bounding box.
    double dx = bbox.x_span() * 0.1;
    double dy = bbox.y_span() * 0.1;
    bbox = Bbox_2(bbox.xmin() - dx, bbox.ymin() - dy, bbox.xmax() + dx, bbox.ymax() + dy);
    // Makes sure the bounding box is not degenerate.
    if (bbox.x_span() == 0) bbox += Bbox_2(bbox.xmin() - 1, bbox.ymin(), bbox.xmax() + 1, bbox.ymax());
    if (bbox.y_span() == 0) bbox += Bbox_2(bbox.xmin(), bbox.ymin() - 1, bbox.xmax(), bbox.ymax() + 1);
    return bbox;
  }

  //! converts a parameter space point to a point of the graphics scene.
  Local_point to_local_point(const Point& pt) const { return Local_point(pt.x(), pt.y(), 0.0); }

private:
  static bool is_empty(const Bbox_2& bbox) { return (bbox.xmin() > bbox.xmax()) || (bbox.ymin() > bbox.ymax()); }
};

// Specialization for arrangements on the sphere
template <typename Arrangement>
class Arr_parameter_space<Arrangement,
                          std::enable_if_t<is_or_derived_from_agas_v<typename Arrangement::Geometry_traits_2>>> {
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Point = typename Approx_traits::Point;
  using Local_point = Graphics_scene::Local_point;

public:
  explicit Arr_parameter_space(const Geom_traits& traits) : m_coords(traits) {}

  /*! computes the approximation error based on the bounding box and the
   * number of pixels that span its width.
   */
  double approximation_error(const Bbox_2& bbox, int resolution) const {
    // If crossing hemisphere
    if (bbox.x_span() >= CGAL_PI) return 1.0 / resolution;
    // Otherwise we evaluate the error bound with respect to the longest longitude arc
    double theta =
      std::abs(bbox.ymin() - CGAL_PI / 2.0) < std::abs(bbox.ymax() - CGAL_PI / 2.0) ? bbox.ymin() : bbox.ymax();
    return bbox.x_span() * std::sin(theta) / resolution;
  }

  //! obtains the entire parameter space.
  Bbox_2 bbox(const Arrangement&) const { return Bbox_2(0, 0, 2 * CGAL_PI, CGAL_PI); }

  //! converts a parameter space point to a point on the unit sphere.
  Local_point to_local_point(const Point& pt) const {
    auto approx_pt = m_coords.to_cartesian(pt);
    return Local_point(approx_pt.dx(), approx_pt.dy(), approx_pt.dz());
  }

private:
  Arr_coordinate_converter<Geom_traits> m_coords;
};

/*! constructs a map from the elements of one range to the corresponding (in
 * iteration order) elements of another range of the same size.
 */
template <typename Value1, typename Value2, typename Range1, typename Range2>
auto map_from_pair_ranges(Range1 range1, Range2 range2) {
  CGAL_assertion_msg(range1.size() == range2.size(), "The two ranges must have the same size.");
  auto begin = boost::make_zip_iterator(boost::make_tuple(range1.begin(), range2.begin()));
  auto end = boost::make_zip_iterator(boost::make_tuple(range1.end(), range2.end()));
  auto tuple_to_pair = [](const auto& t) { return std::make_pair(boost::get<0>(t), boost::get<1>(t)); };
  return unordered_flat_map<Value1, Value2>(boost::make_transform_iterator(begin, tuple_to_pair),
                                            boost::make_transform_iterator(end, tuple_to_pair));
}

/*! tracks changes between an arrangement and its copy that will be later inserted to.
 *
 * \note tracks insertions only. If any other actions made (e.g. deletions,
 * merging, etc), the state of the tracker instance may become invalid.
 *
 * \tparam Arrangement
 */
template <typename Arrangement>
class Arr_insertion_tracker : Arr_observer<Arrangement> {
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

  virtual void before_split_edge(Halfedge_handle /* e */, Vertex_handle v,
                                 const X_monotone_curve_2& /* c1 */,
                                 const X_monotone_curve_2& /* c2 */) override
  { if (m_vertex_map.find(v) == m_vertex_map.end()) m_vertex_map[v] = Vertex_const_handle(); }

  virtual void after_split_edge(Halfedge_handle e1, Halfedge_handle e2) override {
    if (auto it = m_halfedge_map.find(e1); it == m_halfedge_map.end())
      m_halfedge_map[e2] = e1;
    else if (it->second == Halfedge_const_handle())
      m_halfedge_map[e2] = Halfedge_const_handle(); // e1 has no corresponding edge in the original arrangement
    else
      m_halfedge_map[e2] = it->second; // e1 is created by splitting an existing edge
  }

  virtual void after_split_face(Face_handle f1, Face_handle f2, bool) override {
    // Face cannot be created but by splitting an existing face.
    if (auto it = m_face_map.find(f1); it == m_face_map.end())
      m_face_map[f2] = f1;
    else
      m_face_map[f2] = it->second; // f1 is created by splitting an existing face
  }

public:
  Arr_insertion_tracker(Arrangement& arr) : Base(arr) {}

  /*! queries the original face of a given face.
   *
   * \param fh a valid face handle in the modified arrangement.
   */
  Face_const_handle original_face(Face_const_handle fh) const {
    auto it = m_face_map.find(fh);
    if (it == m_face_map.end()) return fh;
    return it->second; // new face from splitting an existing face
  }

  /*! queries the original halfedge of a given halfedge.
   *
   * \param heh a valid halfedge handle in the modified arrangement.
   * \return the original halfedge or an invalid handle if `he` is new.
   */
  Halfedge_const_handle original_halfedge(Halfedge_const_handle he) const {
    auto it = m_halfedge_map.find(he);
    if (it == m_halfedge_map.end()) return he;
    return it->second;
  }

  /*! queries the original vertex of a given vertex.
   *
   * \param vh a valid vertex handle in the modified arrangement.
   * \return the original vertex or an invalid handle if `vh` is new.
   */
  Vertex_const_handle original_vertex(Vertex_const_handle vh) const {
    auto it = m_vertex_map.find(vh);
    if (it == m_vertex_map.end()) return vh;
    return it->second;
  }

private:
  /*! maps tracking the changes between the original arrangement and modified arrangement.
   * The key is the current feature, and the value is the corresponding feature before modification.
   * If there is no entry about a feature, the corresponding feature is itself.
   * If the value is a invalid handle, it means that the feature is newly created and thus has no corresponding
   * feature in the original arrangement.
   */
  unordered_flat_map<Face_const_handle, Face_const_handle> m_face_map;
  unordered_flat_map<Halfedge_const_handle, Halfedge_const_handle> m_halfedge_map;
  unordered_flat_map<Vertex_const_handle, Vertex_const_handle> m_vertex_map;
};

/*! \brief A copy of an arrangement on the sphere, into which the
 * identification curve is inserted, together with a mapping from the features
 * of the copy to the features of the original arrangement.
 *
 * The identification curve makes the sphere homeomorphic to a bounded plane,
 * which is required by the renderer. Features that are introduced by the
 * insertion of the identification curve and do not originate from features of
 * the original arrangement are mapped to invalid handles.
 *
 * \tparam Arrangement
 */
template <typename Arrangement>
class Arr_identified_arrangement {
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;

public:
  explicit Arr_identified_arrangement(const Arrangement& arr) : m_arr(arr) {
    // Map the features of the copy to the features of the original arrangement.
    auto vertex_map =
      map_from_pair_ranges<Vertex_const_handle, Vertex_const_handle>(m_arr.vertex_handles(), arr.vertex_handles());
    auto halfedge_map = map_from_pair_ranges<Halfedge_const_handle, Halfedge_const_handle>(m_arr.halfedge_handles(),
                                                                                           arr.halfedge_handles());
    auto face_map =
      map_from_pair_ranges<Face_const_handle, Face_const_handle>(m_arr.face_handles(), arr.face_handles());

    // Insert the identification curve while tracking the changes.
    Arr_insertion_tracker<Arrangement> tracker(m_arr);
    insert(m_arr, construct_identification_curve(*arr.geometry_traits()));

    // Compose the two mappings, so that each query takes a single lookup.
    for (auto vh : m_arr.vertex_handles()) {
      auto it = vertex_map.find(tracker.original_vertex(vh));
      if (it != vertex_map.end()) m_vertices.emplace(vh, it->second);
    }
    for (auto he : m_arr.halfedge_handles()) {
      auto it = halfedge_map.find(tracker.original_halfedge(he));
      if (it != halfedge_map.end()) m_halfedges.emplace(he, it->second);
    }
    // The tracker records only one of the two halfedges that result from a
    // split; obtain the other one through the twin.
    for (auto he : m_arr.halfedge_handles()) {
      if (m_halfedges.find(he) != m_halfedges.end()) continue;
      auto it = m_halfedges.find(he->twin());
      if (it != m_halfedges.end()) m_halfedges.emplace(he, it->second->twin());
    }
    for (auto fh : m_arr.face_handles()) {
      auto it = face_map.find(tracker.original_face(fh));
      if (it != face_map.end()) m_faces.emplace(fh, it->second);
    }
  }

  Arr_identified_arrangement(const Arr_identified_arrangement&) = delete;
  Arr_identified_arrangement& operator=(const Arr_identified_arrangement&) = delete;

  //! obtains the arrangement with the identification curve.
  const Arrangement& arrangement() const { return m_arr; }

  //! obtains the original vertex, or an invalid handle if there is none.
  Vertex_const_handle original(Vertex_const_handle vh) const { return find(m_vertices, vh); }

  //! obtains the original halfedge, or an invalid handle if there is none.
  Halfedge_const_handle original(Halfedge_const_handle he) const { return find(m_halfedges, he); }

  //! obtains the original face, or an invalid handle if there is none.
  Face_const_handle original(Face_const_handle fh) const { return find(m_faces, fh); }

private:
  template <typename Map, typename Handle>
  static Handle find(const Map& map, Handle h) {
    auto it = map.find(h);
    return (it == map.end()) ? Handle() : it->second;
  }

  Arrangement m_arr;
  unordered_flat_map<Vertex_const_handle, Vertex_const_handle> m_vertices;
  unordered_flat_map<Halfedge_const_handle, Halfedge_const_handle> m_halfedges;
  unordered_flat_map<Face_const_handle, Face_const_handle> m_faces;
};

/*! \brief Graphics scene options of an identified arrangement (see
 * `Arr_identified_arrangement`) that forward every query to the graphics
 * scene options of the original arrangement.
 *
 * Features that do not originate from features of the original arrangement
 * are not drawn.
 *
 * \tparam Arrangement
 * \tparam GSOptions the graphics scene options of the original arrangement.
 */
template <typename Arrangement, typename GSOptions>
class Arr_identified_gs_options {
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Identified_arrangement = Arr_identified_arrangement<Arrangement>;

public:
  Arr_identified_gs_options(const Identified_arrangement& ia, const Arrangement& arr, const GSOptions& gso) :
    m_ia(ia),
    m_arr(arr),
    m_gso(gso)
  {}

  bool are_vertices_enabled() const { return m_gso.are_vertices_enabled(); }
  bool are_edges_enabled() const { return m_gso.are_edges_enabled(); }
  bool are_faces_enabled() const { return m_gso.are_faces_enabled(); }

  bool draw_vertex(const Arrangement&, Vertex_const_handle vh) const {
    auto ovh = m_ia.original(vh);
    return (ovh != Vertex_const_handle()) && m_gso.draw_vertex(m_arr, ovh);
  }

  bool colored_vertex(const Arrangement&, Vertex_const_handle vh) const {
    auto ovh = m_ia.original(vh);
    return (ovh != Vertex_const_handle()) && m_gso.colored_vertex(m_arr, ovh);
  }

  CGAL::IO::Color vertex_color(const Arrangement&, Vertex_const_handle vh) const {
    auto ovh = m_ia.original(vh);
    return (ovh != Vertex_const_handle()) ? m_gso.vertex_color(m_arr, ovh) : CGAL::IO::Color();
  }

  bool draw_edge(const Arrangement&, Halfedge_const_handle he) const {
    auto ohe = m_ia.original(he);
    return (ohe != Halfedge_const_handle()) && m_gso.draw_edge(m_arr, ohe);
  }

  bool colored_edge(const Arrangement&, Halfedge_const_handle he) const {
    auto ohe = m_ia.original(he);
    return (ohe != Halfedge_const_handle()) && m_gso.colored_edge(m_arr, ohe);
  }

  CGAL::IO::Color edge_color(const Arrangement&, Halfedge_const_handle he) const {
    auto ohe = m_ia.original(he);
    return (ohe != Halfedge_const_handle()) ? m_gso.edge_color(m_arr, ohe) : CGAL::IO::Color();
  }

  bool draw_face(const Arrangement&, Face_const_handle fh) const {
    auto ofh = m_ia.original(fh);
    return (ofh != Face_const_handle()) && m_gso.draw_face(m_arr, ofh);
  }

  bool colored_face(const Arrangement&, Face_const_handle fh) const {
    auto ofh = m_ia.original(fh);
    return (ofh != Face_const_handle()) && m_gso.colored_face(m_arr, ofh);
  }

  CGAL::IO::Color face_color(const Arrangement&, Face_const_handle fh) const {
    auto ofh = m_ia.original(fh);
    return (ofh != Face_const_handle()) ? m_gso.face_color(m_arr, ofh) : CGAL::IO::Color();
  }

private:
  const Identified_arrangement& m_ia;
  const Arrangement& m_arr;
  GSOptions m_gso;
};

/*! \brief The arrangement that is actually rendered and its graphics scene
 * options.
 *
 * In the planar case this is the given arrangement and its options. In the
 * spherical case this is a copy of the given arrangement into which the
 * identification curve is inserted (see `Arr_identified_arrangement`), and
 * options that forward to the given options.
 */
template <typename Arrangement, typename GSOptions, typename = void>
class Arr_renderable {
public:
  using Gs_options = GSOptions;

  Arr_renderable(const Arrangement& arr, const GSOptions& gso) : m_arr(arr), m_gso(gso) {}

  const Arrangement& arrangement() const { return m_arr; }
  const Gs_options& gs_options() const { return m_gso; }

private:
  const Arrangement& m_arr;
  Gs_options m_gso;
};

template <typename Arrangement, typename GSOptions>
class Arr_renderable<Arrangement, GSOptions,
                     std::enable_if_t<is_or_derived_from_agas_v<typename Arrangement::Geometry_traits_2>>> {
public:
  using Gs_options = Arr_identified_gs_options<Arrangement, GSOptions>;

  Arr_renderable(const Arrangement& arr, const GSOptions& gso) : m_ia(arr), m_gso(m_ia, arr, gso) {}

  // m_gso refers to m_ia; thus, an object of this type must not be copied or moved.
  Arr_renderable(const Arr_renderable&) = delete;
  Arr_renderable& operator=(const Arr_renderable&) = delete;

  const Arrangement& arrangement() const { return m_ia.arrangement(); }
  const Gs_options& gs_options() const { return m_gso; }

private:
  Arr_identified_arrangement<Arrangement> m_ia;
  Gs_options m_gso;
};

/*! \brief Builds the graphics scene of an arrangement.
 *
 * The arrangement is approximated within a bounding box in the parameter space
 * with a given approximation error, and the vertices, edges, and faces of the
 * approximation are added to a graphics scene. The faces are triangulated.
 *
 * The builder holds a reference to the arrangement; the arrangement must not
 * be modified (or destroyed) while the builder is in use.
 *
 * \tparam Arrangement
 * \tparam GSOptions a model of the `GraphicsSceneOptions` concept.
 */
template <typename Arrangement, typename GSOptions>
class Arr_graphics_scene_builder {
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Point = typename Approx_traits::Point;
  using Local_point = Graphics_scene::Local_point;
  using Point_generator = Arr_face_point_generator<Arrangement>;
  using Renderable = Arr_renderable<Arrangement, GSOptions>;

  static_assert(has_approximate_traits_v<Geom_traits>,
                "The geometry traits must support the approximation of points and x-monotone curves.");

public:
  //! The default resolution (in pixels across the bounding box) for a static scene.
  static constexpr int default_resolution = 2000;

  Arr_graphics_scene_builder(const Arrangement& arr, const GSOptions& gso) :
    m_renderable(arr, gso),
    m_param_space(*m_renderable.arrangement().geometry_traits())
  {}

  //! obtains a bounding box in the parameter space that contains the arrangement.
  Bbox_2 bbox() const { return m_param_space.bbox(m_renderable.arrangement()); }

  //! computes the approximation error for a bounding box and a resolution.
  double approximation_error(const Bbox_2& bbox, int resolution = default_resolution) const
  { return m_param_space.approximation_error(bbox, resolution); }

  /*! adds the part of the arrangement within `bbox` to `gs`, where the curves
   * are approximated with the error `approx_error`.
   */
  void operator()(Graphics_scene& gs, const Bbox_2& bbox, double approx_error) const {
    const Arrangement& arr = m_renderable.arrangement();
    const auto& gso = m_renderable.gs_options();

    auto face_points = Point_generator()(arr, approx_error);
    Arr_render_context<Arrangement> ctx(arr, approx_error, face_points);
    Arr_bounded_renderer<Arrangement> renderer(ctx, bbox);
    auto cache = renderer.render();

    if (gso.are_faces_enabled()) add_faces(gs, arr, gso, cache);
    if (gso.are_edges_enabled()) add_edges(gs, arr, gso, cache, bbox);
    if (gso.are_vertices_enabled()) add_vertices(gs, arr, gso, cache, bbox);
  }

  //! adds the arrangement to `gs`; the bounding box and error are computed.
  void operator()(Graphics_scene& gs) const {
    Bbox_2 bb = bbox();
    operator()(gs, bb, approximation_error(bb));
  }

private:
  static bool contains(const Bbox_2& bbox, const Point& pt)
  { return bbox.xmin() <= pt.x() && pt.x() <= bbox.xmax() && bbox.ymin() <= pt.y() && pt.y() <= bbox.ymax(); }

  template <typename Gso, typename Cache>
  void add_faces(Graphics_scene& gs, const Arrangement& arr, const Gso& gso, const Cache& cache) const {
    std::vector<Local_point> points;
    for (const auto& [fh, tf] : cache.faces()) {
      if (tf.triangles.empty() || ! gso.draw_face(arr, fh)) continue;
      bool colored = gso.colored_face(arr, fh);
      auto color = colored ? gso.face_color(arr, fh) : CGAL::IO::Color();
      // Convert every point once; a point is typically shared by several triangles.
      points.clear();
      points.reserve(tf.points.size());
      for (const auto& pt : tf.points) points.push_back(m_param_space.to_local_point(pt));
      for (const auto& tri : tf.triangles) {
        if (colored) gs.face_begin(color);
        else gs.face_begin();
        for (const auto i : tri) gs.add_point_in_face(points[i]);
        gs.face_end();
      }
    }
  }

  template <typename Gso, typename Cache>
  void add_edges(Graphics_scene& gs, const Arrangement& arr, const Gso& gso, const Cache& cache,
                 const Bbox_2& bbox) const {
    for (const auto& [he, polyline] : cache.halfedges()) {
      // Draw every edge once.
      if (he->direction() == ARR_RIGHT_TO_LEFT || polyline.size() < 2 || ! gso.draw_edge(arr, he)) continue;
      bool colored = gso.colored_edge(arr, he);
      auto color = colored ? gso.edge_color(arr, he) : CGAL::IO::Color();
      // skip first two if starts with a sep point.
      std::size_t start_idx = Approx_traits::is_null(polyline.front()) ? 2 : 0;
      // skip last two if ends with a sep point.
      std::size_t end_idx = Approx_traits::is_null(polyline.back()) ? polyline.size() - 2 : polyline.size();
      // The polyline may consist of separator pairs only (e.g., when the curve
      // has no part inside the bounding box), in which case end_idx <= 1 or
      // end_idx < start_idx. Avoid computing end_idx - 1, which may underflow.
      for (std::size_t i = start_idx; i + 1 < end_idx; ++i) {
        const auto& src = polyline[i];
        const auto& tgt = polyline[i + 1];
        if (Approx_traits::is_null(src) || Approx_traits::is_null(tgt)) continue;
        if (! contains(bbox, src) || ! contains(bbox, tgt)) continue;
        if (colored)
          gs.add_segment(m_param_space.to_local_point(src), m_param_space.to_local_point(tgt), color);
        else
          gs.add_segment(m_param_space.to_local_point(src), m_param_space.to_local_point(tgt));
      }
    }
  }

  template <typename Gso, typename Cache>
  void add_vertices(Graphics_scene& gs, const Arrangement& arr, const Gso& gso, const Cache& cache,
                    const Bbox_2& bbox) const {
    for (const auto& [vh, pt] : cache.vertices()) {
      if (! contains(bbox, pt) || ! gso.draw_vertex(arr, vh)) continue;
      if (gso.colored_vertex(arr, vh)) gs.add_point(m_param_space.to_local_point(pt), gso.vertex_color(arr, vh));
      else gs.add_point(m_param_space.to_local_point(pt));
    }
  }

  Renderable m_renderable;
  Arr_parameter_space<Arrangement> m_param_space;
};

/*! obtains the default graphics scene options of an arrangement: vertices are
 * red, edges are black, and faces are colored randomly (but consistently).
 */
template <typename Arrangement>
auto default_gs_options() {
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using GSOptions =
    CGAL::Graphics_scene_options<Arrangement, Vertex_const_handle, Halfedge_const_handle, Face_const_handle>;

  GSOptions gso;
  gso.colored_vertex = [](const Arrangement&, const Vertex_const_handle&) { return true; };
  gso.vertex_color = [](const Arrangement&, const Vertex_const_handle&) { return CGAL::IO::Color(255, 0, 0); };
  gso.colored_edge = [](const Arrangement&, const Halfedge_const_handle&) { return true; };
  gso.edge_color = [](const Arrangement&, const Halfedge_const_handle&) { return CGAL::IO::Color(0, 0, 0); };
  gso.colored_face = [](const Arrangement&, const Face_const_handle&) { return true; };
  gso.face_color = [](const Arrangement&, const Face_const_handle& fh) {
    CGAL::Random random(static_cast<unsigned int>(reinterpret_cast<std::size_t>(fh.ptr())));
    return get_random_color(random);
  };
  return gso;
}

} // namespace draw_aos
} // namespace CGAL

#endif
